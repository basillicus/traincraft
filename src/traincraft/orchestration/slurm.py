"""Slurm + Apptainer executor: the same stages as dependency-chained HPC jobs.

For each enabled stage it renders an ``sbatch`` script that ``apptainer exec``s
the right image and runs ``traincraft stage <name> CONFIG``. The scripts are
chained with ``--dependency=afterok`` so the pipeline runs as a DAG on the
cluster. The science is untouched — this only wires the stage functions to Slurm.

Image map (overridable per stage in ``[orchestration.slurm.stages.*]``):
    sample -> traincraft-mlip.sif  (GPU/Booster, --nv)
    others -> traincraft-core.sif  (CPU)

The ``label`` stage runs ``traincraft`` in the *core* image but injects the DFT
engine command (``TRAINCRAFT_AIMS_COMMAND`` / ``TRAINCRAFT_PW_COMMAND``) so the
heavy MPI binary runs in the *dft* image under ``srun`` — keeping the calculator
plugins container-agnostic (DESIGN §20.3).
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
_DEFAULT_AIMS_CMD = "srun apptainer exec{binds} {sif_dir}/traincraft-dft.sif aims.x"
_DEFAULT_PW_CMD = "srun apptainer exec{binds} {sif_dir}/traincraft-dft.sif pw.x"

_JOB_ID_RE = re.compile(r"Submitted batch job (\d+)")


@dataclass
class StageJob:
    stage: str
    script: Path
    depends_on: str | None = None  # name of the stage this one waits for
    job_id: str | None = field(default=None)


def _bind_flags(binds: list[str]) -> str:
    return "".join(f" --bind {b}" for b in binds)


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
    binds = _bind_flags(slurm.binds)

    body = ["", "set -euo pipefail"]
    body += [f"module load {m}" for m in slurm.modules]
    body += [f'export {k}="{v}"' for k, v in slurm.env.items()]
    body += [f'export {k}="{v}"' for k, v in st.env.items()]

    if stage == "label":
        aims = (slurm.aims_command or _DEFAULT_AIMS_CMD).format(binds=binds, sif_dir=slurm.sif_dir)
        pw = (slurm.pw_command or _DEFAULT_PW_CMD).format(binds=binds, sif_dir=slurm.sif_dir)
        body.append(f'export TRAINCRAFT_AIMS_COMMAND="{aims}"')
        body.append(f'export TRAINCRAFT_PW_COMMAND="{pw}"')

    nv = " --nv" if gpus else ""
    exec_line = (
        f"apptainer exec{nv}{binds} {slurm.sif_dir}/{image} "
        f"traincraft stage {stage} {config_path}"
    )
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
