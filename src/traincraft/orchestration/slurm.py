"""Slurm + Apptainer executor: the same stages as dependency-chained HPC jobs.

For each enabled stage it renders an ``sbatch`` script that ``apptainer exec``s
the right image and runs ``traincraft stage <name> CONFIG``. The scripts are
chained with ``--dependency=afterok`` so the pipeline runs as a DAG on the
cluster. The science is untouched — this only wires the stage functions to Slurm.

Image map (overridable per stage in ``[orchestration.slurm.stages.*]``):
    sample -> traincraft-mlip.sif  (GPU/Booster, --nv)
    others -> traincraft-core.sif  (CPU)

Two orthogonal axes make this portable across clusters (``[orchestration.slurm]``):

* ``runtime`` — ``apptainer`` (our images) or ``native`` (binaries already on the
  host: site modules / conda / EasyBuild). ``native`` drops the container wrapper,
  so the same DAG runs on machines where bind-mounting our images is awkward
  (e.g. Cray, where the site MPI/binaries are the path of least resistance).
* ``mpi`` — the Slurm MPI plugin for the multi-node DFT launch: ``pmix`` for
  InfiniBand+Slurm (Leonardo), ``cray_shasta`` for Cray/Slingshot (LUMI has no
  pmix), ``pmi2`` as a portable fallback. There is no universal MPI plugin —
  ``srun --mpi=list`` on the target cluster is the ground truth.

The ``label`` stage runs ``traincraft`` in the *core* image (or natively) but
injects the DFT engine command (``TRAINCRAFT_AIMS_COMMAND`` / ``TRAINCRAFT_PW_COMMAND``)
so the heavy MPI binary launches under ``srun --mpi=<plugin>`` — keeping the
calculator plugins container-agnostic (DESIGN §20.3). With ``runtime=apptainer``
the command wraps FHI-aims in ``traincraft-dft.sif`` and Quantum ESPRESSO in
``traincraft-qe.sif``; with ``runtime=native`` it calls the bare binary. Both are
fully overridable via ``[orchestration.slurm].aims_command``/``pw_command``.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from ..config import TrainCraftConfig
from ..core import Workspace
from .stages import enabled_stages, workspace_for

logger = logging.getLogger(__name__)

_DEFAULT_IMAGE = {
    "geometry": "traincraft-core.sif",
    "sample": "traincraft-mlip.sif",
    "select": "traincraft-core.sif",
    "label": "traincraft-core.sif",
    "dataset": "traincraft-core.sif",
}
_DEFAULT_GPUS = {"sample": 1}

_JOB_ID_RE = re.compile(r"Submitted batch job (\d+)")


@dataclass
class StageJob:
    stage: str
    script: Path
    depends_on: str | None = None  # name of the stage this one waits for
    job_id: str | None = field(default=None)


def _bind_flags(binds: list[str]) -> str:
    return "".join(f" --bind {b}" for b in binds)


def _launch_prefix(mpi: str) -> str:
    """``srun`` with the cluster's MPI plugin (omitted when ``mpi='none'``)."""
    return "srun" if mpi == "none" else f"srun --mpi={mpi}"


def _dft_command(binary: str, image: str, slurm, runtime: str, mpi: str) -> str:
    """Compose the injected DFT engine command for the requested runtime.

    apptainer: ``srun --mpi=<plugin> apptainer exec <binds> <sif> <binary>``
    native:    ``srun --mpi=<plugin> <binary>``  (binary from host modules/conda)
    """
    prefix = _launch_prefix(mpi)
    if runtime == "native":
        return f"{prefix} {binary}"
    return f"{prefix} apptainer exec{_bind_flags(slurm.binds)} {slurm.sif_dir}/{image} {binary}"


def _stage_exec(stage, image, runtime, gpus, slurm, config_path) -> str:
    """The line that runs ``traincraft stage`` — wrapped in our image or bare."""
    inner = f"traincraft stage {stage} {config_path}"
    if runtime == "native":
        return inner
    nv = " --nv" if gpus else ""
    return f"apptainer exec{nv}{_bind_flags(slurm.binds)} {slurm.sif_dir}/{image} {inner}"


def _sbatch_header(stage, st, slurm, run_name, ws: Workspace) -> list[str]:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={run_name}-{stage}",
        f"#SBATCH --output={ws.root}/slurm/{stage}-%j.out",
        f"#SBATCH --error={ws.root}/slurm/{stage}-%j.err",
        f"#SBATCH --nodes={st.nodes}",
        f"#SBATCH --time={st.time}",
    ]
    if slurm.account:
        lines.append(f"#SBATCH --account={slurm.account}")
    if st.partition:
        lines.append(f"#SBATCH --partition={st.partition}")
    if st.ntasks is not None:
        lines.append(f"#SBATCH --ntasks={st.ntasks}")
    if st.cpus_per_task is not None:
        lines.append(f"#SBATCH --cpus-per-task={st.cpus_per_task}")
    gpus = st.gpus if st.gpus is not None else _DEFAULT_GPUS.get(stage)
    if gpus:
        lines.append(f"#SBATCH --gpus={gpus}")
    if st.mem:
        lines.append(f"#SBATCH --mem={st.mem}")
    if st.qos:
        lines.append(f"#SBATCH --qos={st.qos}")
    lines.extend(f"#SBATCH {extra}" for extra in st.extra_sbatch)
    return lines


def render_sbatch(stage: str, config: TrainCraftConfig, ws: Workspace, config_path: str) -> str:
    from ..config.models import SlurmStage

    slurm = config.orchestration.slurm
    st = slurm.stages.get(stage, SlurmStage())
    image = st.image or _DEFAULT_IMAGE[stage]
    gpus = st.gpus if st.gpus is not None else _DEFAULT_GPUS.get(stage)
    runtime = st.runtime or slurm.runtime
    mpi = st.mpi or slurm.mpi
    binds = _bind_flags(slurm.binds)

    body = ["", "set -euo pipefail"]
    body += [f"module load {m}" for m in slurm.modules]
    body += list(slurm.pre_commands)
    body += list(st.pre_commands)
    body += [f'export {k}="{v}"' for k, v in slurm.env.items()]
    body += [f'export {k}="{v}"' for k, v in st.env.items()]

    if stage == "label":
        if slurm.aims_command:
            aims = slurm.aims_command.format(binds=binds, sif_dir=slurm.sif_dir, mpi=mpi)
        else:
            aims = _dft_command("aims.x", slurm.aims_image, slurm, runtime, mpi)
        if slurm.pw_command:
            pw = slurm.pw_command.format(binds=binds, sif_dir=slurm.sif_dir, mpi=mpi)
        else:
            pw = _dft_command("pw.x", slurm.qe_image, slurm, runtime, mpi)
        body.append(f'export TRAINCRAFT_AIMS_COMMAND="{aims}"')
        body.append(f'export TRAINCRAFT_PW_COMMAND="{pw}"')

    exec_line = _stage_exec(stage, image, runtime, gpus, slurm, config_path)
    body += ["", exec_line, ""]

    header = _sbatch_header(stage, st, slurm, config.run.name, ws)
    return "\n".join(header + body) + "\n"


def render(
    config: TrainCraftConfig, config_path: str, ws: Workspace | None = None
) -> list[StageJob]:
    """Write one sbatch script per enabled stage; return the dependency-chained plan."""
    if config.orchestration is None or config.orchestration.slurm is None:
        raise ValueError("slurm executor needs an [orchestration.slurm] section")
    ws = ws or workspace_for(config)
    slurm_dir = ws.subdir("slurm")
    cfg_abs = str(Path(config_path).resolve())

    jobs: list[StageJob] = []
    prev: str | None = None
    for stage in enabled_stages(config):
        text = render_sbatch(stage, config, ws, cfg_abs)
        path = slurm_dir / f"{stage}.sbatch"
        path.write_text(text)
        jobs.append(StageJob(stage=stage, script=path, depends_on=prev))
        prev = stage
    return jobs


def submit(
    config: TrainCraftConfig,
    config_path: str,
    ws: Workspace | None = None,
    *,
    dry_run: bool = False,
) -> list[StageJob]:
    """Render scripts and (unless ``dry_run``) submit them with afterok chaining."""
    ws = ws or workspace_for(config)
    jobs = render(config, config_path, ws)

    if dry_run or shutil.which("sbatch") is None:
        if not dry_run:
            logger.warning(
                "sbatch not found; wrote scripts to %s without submitting", ws.root / "slurm"
            )
        return jobs

    prev_id: str | None = None
    for job in jobs:
        cmd = ["sbatch"]
        if prev_id is not None:
            cmd.append(f"--dependency=afterok:{prev_id}")
        cmd.append(str(job.script))
        out = subprocess.run(cmd, capture_output=True, text=True, check=True)
        match = _JOB_ID_RE.search(out.stdout)
        if not match:
            raise RuntimeError(f"could not parse job id from sbatch output: {out.stdout!r}")
        job.job_id = prev_id = match.group(1)
        logger.info("submitted %s as job %s (after %s)", job.stage, job.job_id, job.depends_on)
    return jobs
