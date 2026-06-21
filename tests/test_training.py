"""Tests for the MACE training wrapper.

MACE/torch are not installed in CI, and training takes hours, so nothing is
actually trained: we verify the *rendered* ``mace_run_train`` command, the
train/valid split + ``REF_*`` re-keying, command injection, and the manifest —
mirroring ``test_dft`` (inspect, don't run) and ``test_slurm_executor`` (render,
don't submit). A fake subprocess stands in for the real binary in the run path.
"""

from __future__ import annotations

import json

import pytest

from traincraft.config.models import MaceTrainer
from traincraft.core import RegistryError, Structure, Workspace, get
from traincraft.training import run_training, write_training_xyz
from traincraft.training.mace import (
    DEFAULT_TRAIN_COMMAND,
    TRAIN_COMMAND_ENV,
    _model_and_loss,
    _split,
    render_command,
    train_mace,
)

pytest.importorskip("ase")
from ase.build import bulk, molecule  # noqa: E402
from ase.io import read as ase_read  # noqa: E402


def _flags(argv: list[str]) -> dict[str, str | bool]:
    """Parse ``--flag value`` / ``--flag`` argv into a dict (store_true -> True)."""
    out: dict[str, str | bool] = {}
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok.startswith("--"):
            if i + 1 < len(argv) and not argv[i + 1].startswith("--"):
                out[tok[2:]] = argv[i + 1]
                i += 2
            else:
                out[tok[2:]] = True
                i += 1
        else:
            i += 1
    return out


def _labelled_frames(n=6):
    frames = []
    for i in range(n):
        atoms = bulk("Cu", "fcc", a=3.6 + 0.001 * i, cubic=True)
        s = Structure(atoms=atoms, properties={
            "energy": -3.5 - 0.01 * i,
            "forces": [[0.0, 0.0, 0.01 * i]] * len(atoms),
            "stress": [[0.1, 0, 0], [0, 0.1, 0], [0, 0, 0.1]],
        })
        frames.append(s)
    return frames


# --------------------------------------------------------------------------- registry
def test_trainer_registered():
    assert get("trainer", "mace") is train_mace


def test_unknown_trainer_raises():
    with pytest.raises(RegistryError):
        get("trainer", "nequip")


# ----------------------------------------------------------------- render: defaults
def test_render_defaults_follow_paper(tmp_path):
    cfg = MaceTrainer()
    argv = render_command(cfg, tmp_path / "train.xyz", tmp_path / "valid.xyz", tmp_path)
    f = _flags(argv)
    assert argv[0] == DEFAULT_TRAIN_COMMAND
    # paper-informed fine-tuning defaults
    assert f["E0s"] == "foundation"
    assert f["weight_decay"] == "0.0"
    assert f["ema"] is True and f["ema_decay"] == "0.995"
    assert f["energy_weight"] == "10.0" and f["forces_weight"] == "10.0"
    assert f["foundation_model"] == "medium"
    assert f["multiheads_finetuning"] == "True"
    # reference keys match write_training_xyz
    assert f["energy_key"] == "REF_energy" and f["forces_key"] == "REF_forces"
    # plain energy/forces -> no explicit model/loss (foundation decides)
    assert "model" not in f and "loss" not in f
    assert f["train_file"].endswith("train.xyz")
    assert f["valid_file"].endswith("valid.xyz")


def test_render_omits_valid_when_none(tmp_path):
    argv = render_command(MaceTrainer(), tmp_path / "train.xyz", None, tmp_path)
    assert "--valid_file" not in argv


# --------------------------------------------------------------- render: strategies
def test_naive_disables_multihead(tmp_path):
    argv = render_command(MaceTrainer(strategy="naive"), tmp_path / "t.xyz", None, tmp_path)
    f = _flags(argv)
    assert f["multiheads_finetuning"] == "False"
    assert "pt_train_file" not in f


def test_multihead_adds_replay(tmp_path):
    cfg = MaceTrainer(strategy="multihead", pt_train_file="mp", num_samples_pt=5000)
    f = _flags(render_command(cfg, tmp_path / "t.xyz", None, tmp_path))
    assert f["pt_train_file"] == "mp"
    assert f["num_samples_pt"] == "5000"
    assert "weight_pt" in f and "weight_ft" in f


def test_scratch_has_no_foundation(tmp_path):
    cfg = MaceTrainer(strategy="scratch", foundation_model=None,
                      hidden_irreps="128x0e + 128x1o", r_max=4.5)
    f = _flags(render_command(cfg, tmp_path / "t.xyz", None, tmp_path))
    assert "foundation_model" not in f
    assert "multiheads_finetuning" not in f
    assert f["hidden_irreps"] == "128x0e + 128x1o"
    assert f["r_max"] == "4.5"


# --------------------------------------------------------------- multi-head mapping
def test_model_and_loss_mapping():
    assert _model_and_loss({"energy", "forces"}) == (None, None)
    assert _model_and_loss({"dipole"}) == ("AtomicDipolesMACE", "dipole")
    assert _model_and_loss({"energy", "forces", "dipole"}) == (
        "EnergyDipolesMACE", "energy_forces_dipole")
    assert _model_and_loss({"energy", "forces", "dipole", "polarizability"}) == (
        "AtomicDielectricMACE", "dipole_polar")


def test_render_dielectric_heads(tmp_path):
    cfg = MaceTrainer(heads=["energy", "forces", "dipole", "polarizability"])
    f = _flags(render_command(cfg, tmp_path / "t.xyz", None, tmp_path))
    assert f["model"] == "AtomicDielectricMACE"
    assert f["loss"] == "dipole_polar"
    assert f["dipole_key"] == "REF_dipole"
    assert f["polarizability_key"] == "REF_polarizability"


def test_render_stress_head(tmp_path):
    f = _flags(render_command(
        MaceTrainer(heads=["energy", "forces", "stress"]), tmp_path / "t.xyz", None, tmp_path))
    assert f["stress_key"] == "REF_stress"
    assert "stress_weight" in f


# ------------------------------------------------------------- command injection
def test_command_injection(tmp_path, monkeypatch):
    injected = "srun --nv apptainer exec --bind /scratch traincraft-mlip.sif mace_run_train"
    monkeypatch.setenv(TRAIN_COMMAND_ENV, injected)
    argv = render_command(MaceTrainer(), tmp_path / "t.xyz", None, tmp_path)
    assert argv[:5] == ["srun", "--nv", "apptainer", "exec", "--bind"]
    # the injected launcher precedes all the --flags we append
    assert argv[:8] == [
        "srun", "--nv", "apptainer", "exec", "--bind", "/scratch",
        "traincraft-mlip.sif", "mace_run_train",
    ]


# ------------------------------------------------------------------ extra passthrough
def test_extra_overrides_default(tmp_path):
    cfg = MaceTrainer(extra={"lr": 0.001234, "scheduler": "cosine", "amsgrad": True})
    f = _flags(render_command(cfg, tmp_path / "t.xyz", None, tmp_path))
    assert f["lr"] == "0.001234"  # overrode the default --lr
    assert f["scheduler"] == "cosine"
    assert f["amsgrad"] is True
    # no duplicate --lr
    assert render_command(cfg, tmp_path / "t.xyz", None, tmp_path).count("--lr") == 1


# ----------------------------------------------------------------------- split + io
def test_split_deterministic_and_nonempty():
    frames = _labelled_frames(10)
    train, valid = _split(frames, 0.2, seed=0)
    assert len(train) == 8 and len(valid) == 2
    # deterministic for the same seed
    t2, v2 = _split(frames, 0.2, seed=0)
    assert [s.hash for s in valid] == [s.hash for s in v2]


def test_write_training_xyz_rekeys(tmp_path):
    s = Structure(atoms=molecule("H2O"), properties={
        "energy": -2078.0,
        "forces": [[0, 0, 0.1]] * 3,
        "dipole": [0.0, 0.0, 1.85],
        "polarizability": [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]],
    })
    path = write_training_xyz(tmp_path / "out.xyz", [s])
    atoms = ase_read(str(path))
    assert atoms.info["REF_energy"] == pytest.approx(-2078.0)
    assert "REF_forces" in atoms.arrays and atoms.arrays["REF_forces"].shape == (3, 3)
    assert atoms.info["REF_dipole"][2] == pytest.approx(1.85)
    assert len(atoms.info["REF_polarizability"]) == 9


# --------------------------------------------------------------------- dry run + run
def test_train_dry_run_writes_splits_and_manifest(tmp_path):
    ws = Workspace(tmp_path)
    job = ws.job("model")
    result = run_training(_labelled_frames(6), MaceTrainer(seed=0), job, dry_run=True)
    assert result.model_path is None
    assert result.n_train + result.n_valid == 6
    assert (job.dir / "train.xyz").exists()
    assert (job.dir / "valid.xyz").exists()
    manifest = json.loads((job.dir / "manifest.json").read_text())
    assert manifest["dry_run"] is True
    assert manifest["reference"].startswith("Tompa")
    assert manifest["command"][0] == DEFAULT_TRAIN_COMMAND


def test_train_runs_subprocess_and_finds_model(tmp_path, monkeypatch):
    ws = Workspace(tmp_path)
    job = ws.job("model")

    def fake_run(cmd, check, stdout, stderr):
        f = _flags(cmd)
        from pathlib import Path

        model_dir = Path(f["model_dir"])
        model_dir.mkdir(parents=True, exist_ok=True)
        (model_dir / f"{f['name']}.model").write_text("fake")
        class _R:  # noqa: N801
            returncode = 0
        return _R()

    monkeypatch.setattr("traincraft.training.mace.subprocess.run", fake_run)
    result = run_training(_labelled_frames(4), MaceTrainer(name="cu", seed=0), job)
    assert result.model_path is not None and result.model_path.name == "cu.model"
    manifest = json.loads((job.dir / "manifest.json").read_text())
    assert manifest["model_path"].endswith("cu.model")
    assert "wall_seconds" in manifest


# ------------------------------------------------------------------ stage wiring
def test_train_stage_runs_after_dataset(tmp_path, monkeypatch):
    from traincraft.config.models import (
        CrystalBuilder,
        DatasetConfig,
        EmtCalc,
        GeometryConfig,
        LabelingConfig,
        MdSampling,
        RunConfig,
        SelectionConfig,
        TrainCraftConfig,
    )
    from traincraft.orchestration import enabled_stages, run_stage, workspace_for

    cfg = TrainCraftConfig(
        run=RunConfig(name="trainwire", outdir=str(tmp_path), seed=1),
        geometry=GeometryConfig(builder=CrystalBuilder(
            name="Cu", crystalstructure="fcc", a=3.6, cubic=True, supercell=(2, 2, 2))),
        calculator=EmtCalc(),
        sampling=MdSampling(temperature=300.0, steps=40, interval=10),
        selection=SelectionConfig(budget=4),
        labeling=LabelingConfig(calculator=EmtCalc()),
        dataset=DatasetConfig(),
        training=MaceTrainer(name="cu", strategy="naive", seed=0),
    )
    assert enabled_stages(cfg)[-1] == "train"

    def fake_run(cmd, check, stdout, stderr):
        from pathlib import Path

        f = _flags(cmd)
        d = Path(f["model_dir"])
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{f['name']}.model").write_text("fake")
        class _R:  # noqa: N801
            returncode = 0
        return _R()

    monkeypatch.setattr("traincraft.training.mace.subprocess.run", fake_run)
    ws = workspace_for(cfg)
    for stage in enabled_stages(cfg):
        run_stage(stage, cfg, ws)
    assert (ws.root / "model" / "manifest.json").exists()
    assert (ws.root / "model" / "cu.model").exists()


# --------------------------------------------------------------------- config guards
def test_scratch_rejects_foundation():
    from pydantic import ValidationError

    with pytest.raises(ValidationError):
        MaceTrainer(strategy="scratch")  # foundation_model defaults to "medium"


def test_empty_heads_rejected():
    from pydantic import ValidationError

    with pytest.raises(ValidationError):
        MaceTrainer(heads=[])
