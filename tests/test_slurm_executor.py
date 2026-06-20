from __future__ import annotations

from traincraft.config.models import (
    CrystalBuilder,
    DatasetConfig,
    EmtCalc,
    FhiAimsCalc,
    GeometryConfig,
    LabelingConfig,
    MdSampling,
    OrchestrationConfig,
    RunConfig,
    SelectionConfig,
    SlurmConfig,
    SlurmStage,
    TrainCraftConfig,
)
from traincraft.orchestration import submit_slurm


def _hpc_cfg(tmp_path):
    return TrainCraftConfig(
        run=RunConfig(name="hpc", outdir=str(tmp_path)),
        geometry=GeometryConfig(builder=CrystalBuilder(name="Cu", crystalstructure="fcc", a=3.6)),
        calculator=EmtCalc(),
        sampling=MdSampling(steps=10),
        selection=SelectionConfig(),
        labeling=LabelingConfig(calculator=FhiAimsCalc(xc="pbe")),
        dataset=DatasetConfig(),
        orchestration=OrchestrationConfig(
            engine="slurm",
            slurm=SlurmConfig(
                account="ACC_1",
                sif_dir="/work/sif",
                binds=["/scratch", "/work"],
                modules=["apptainer"],
                stages={
                    "sample": SlurmStage(image="traincraft-mlip.sif", partition="boost", gpus=1),
                    "label": SlurmStage(partition="dcgp", nodes=2, ntasks=224, time="06:00:00"),
                },
            ),
        ),
    )


def test_submit_dry_run_writes_chained_scripts(tmp_path):
    cfg = _hpc_cfg(tmp_path)
    jobs = submit_slurm(cfg, "myconfig.toml", dry_run=True)

    assert [j.stage for j in jobs] == ["geometry", "sample", "select", "label", "dataset"]
    # dependency chain: each waits on the previous
    assert [j.depends_on for j in jobs] == [None, "geometry", "sample", "select", "label"]
    assert all(j.job_id is None for j in jobs)  # nothing actually submitted
    for j in jobs:
        assert j.script.exists()


def test_sample_stage_uses_gpu_mlip_image(tmp_path):
    submit_slurm(_hpc_cfg(tmp_path), "myconfig.toml", dry_run=True)
    text = (tmp_path / "hpc" / "slurm" / "sample.sbatch").read_text()
    assert "traincraft-mlip.sif" in text
    assert "--nv" in text
    assert "#SBATCH --gpus=1" in text
    assert "#SBATCH --partition=boost" in text
    assert "traincraft stage sample" in text
    assert "--bind /scratch" in text and "--bind /work" in text


def test_label_stage_injects_dft_command_and_core_image(tmp_path):
    submit_slurm(_hpc_cfg(tmp_path), "myconfig.toml", dry_run=True)
    text = (tmp_path / "hpc" / "slurm" / "label.sbatch").read_text()
    # traincraft runs in core; FHI-aims runs in the dft image via the injected command
    assert "traincraft-core.sif traincraft stage label" in text
    assert "export TRAINCRAFT_AIMS_COMMAND=" in text
    assert "traincraft-dft.sif aims.x" in text
    # the open-source QE engine is injected too, pointing at its own image
    assert "export TRAINCRAFT_PW_COMMAND=" in text
    assert "traincraft-qe.sif pw.x" in text
    # default MPI plugin is pmix (InfiniBand+Slurm), launched under srun
    assert "srun --mpi=pmix apptainer exec" in text
    assert "#SBATCH --nodes=2" in text
    assert "#SBATCH --ntasks=224" in text
    assert "#SBATCH --account=ACC_1" in text


def test_geometry_stage_defaults_to_core_image(tmp_path):
    submit_slurm(_hpc_cfg(tmp_path), "myconfig.toml", dry_run=True)
    text = (tmp_path / "hpc" / "slurm" / "geometry.sbatch").read_text()
    assert "traincraft-core.sif" in text
    assert "--nv" not in text  # no GPU on the geometry stage


def _cray_native_cfg(tmp_path):
    """A Cray-style cluster (LUMI): native runtime, cray_shasta, no pmix/containers."""
    return TrainCraftConfig(
        run=RunConfig(name="lumi", outdir=str(tmp_path)),
        geometry=GeometryConfig(builder=CrystalBuilder(name="Cu", crystalstructure="fcc", a=3.6)),
        calculator=EmtCalc(),
        sampling=MdSampling(steps=10),
        selection=SelectionConfig(),
        labeling=LabelingConfig(calculator=FhiAimsCalc(xc="pbe")),
        dataset=DatasetConfig(),
        orchestration=OrchestrationConfig(
            engine="slurm",
            slurm=SlurmConfig(
                account="proj_465",
                runtime="native",
                mpi="cray_shasta",
                modules=["LUMI/24.03", "cray-mpich"],
                pre_commands=["source $HOME/tc-venv/bin/activate"],
                stages={
                    "label": SlurmStage(
                        partition="standard",
                        nodes=2,
                        ntasks=256,
                        pre_commands=["module load fhi-aims/240507"],
                    ),
                },
            ),
        ),
    )


def test_native_runtime_drops_apptainer_wrapper(tmp_path):
    submit_slurm(_cray_native_cfg(tmp_path), "myconfig.toml", dry_run=True)
    geo = (tmp_path / "lumi" / "slurm" / "geometry.sbatch").read_text()
    # native: the binary/orchestrator runs directly, no `apptainer exec`, no .sif
    assert "apptainer exec" not in geo
    assert ".sif" not in geo
    assert "traincraft stage geometry" in geo
    assert "source $HOME/tc-venv/bin/activate" in geo


def test_native_label_uses_cray_mpi_and_bare_binary(tmp_path):
    submit_slurm(_cray_native_cfg(tmp_path), "myconfig.toml", dry_run=True)
    text = (tmp_path / "lumi" / "slurm" / "label.sbatch").read_text()
    # Cray MPI plugin, bare site binary (no container), site FHI-aims module loaded
    assert 'export TRAINCRAFT_AIMS_COMMAND="srun --mpi=cray_shasta aims.x"' in text
    assert 'export TRAINCRAFT_PW_COMMAND="srun --mpi=cray_shasta pw.x"' in text
    assert "apptainer exec" not in text
    assert "module load fhi-aims/240507" in text
