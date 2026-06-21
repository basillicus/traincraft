"""MACE trainer: a thin, container-agnostic wrapper over ``mace_run_train``.

``core`` prepares the train/valid extended-XYZ (with explicit ``REF_*`` keys),
renders the ``mace_run_train`` command, and shells out to it. The MACE/torch
stack lives in ``traincraft-mlip.sif`` (or a local ``mace`` env); this module
never imports torch — it only builds argv and runs a subprocess, so the package
stays importable with no ML stack and an HPC executor can inject
``srun --nv apptainer exec … traincraft-mlip.sif mace_run_train`` via
``$TRAINCRAFT_MACE_TRAIN_COMMAND`` (DESIGN §20.3).

Fine-tuning defaults follow Tompa et al., arXiv:2606.12704 (see ``MaceTrainer``):
foundation-consistent E0s, multihead replay against forgetting, ``weight_decay=0``
and high EMA for fine-tuning, constant energy-prioritised loss weights.

Multi-head property targets map onto MACE model types + losses:

    heads include …            --model               --loss
    energy/forces (+stress)    (foundation default)  (foundation default)
    dipole only                AtomicDipolesMACE     dipole
    energy/forces + dipole     EnergyDipolesMACE     energy_forces_dipole
    polarizability             AtomicDielectricMACE  dipole_polar

The dipole/polarizability model types are the ones MACE ships for dielectric
properties (AtomicDipolesMACE / EnergyDipolesMACE; AtomicDielectricMACE, cf.
MACE-MDP / mace-field). They evolve faster than the energy/forces path, so the
exact ``--model``/``--loss`` are overridable through ``MaceTrainer.extra`` —
verify them against your installed MACE version before a production polarizability
run, the same way the QE ``ph.x`` polarizability path is gated in ``dft.py``.
"""

from __future__ import annotations

import json
import logging
import os
import shlex
import subprocess
import time
from pathlib import Path

import numpy as np

from ..core import make_rng, register
from .base import REF_KEYS, TrainResult, write_training_xyz

logger = logging.getLogger(__name__)

# Command injection (see DESIGN §20.3); default to the bare console script.
TRAIN_COMMAND_ENV = "TRAINCRAFT_MACE_TRAIN_COMMAND"
DEFAULT_TRAIN_COMMAND = "mace_run_train"


def _model_and_loss(heads: set[str]) -> tuple[str | None, str | None]:
    """MACE ``--model`` / ``--loss`` for a head set (None = foundation default)."""
    if "polarizability" in heads:
        return "AtomicDielectricMACE", "dipole_polar"
    if "dipole" in heads:
        if heads & {"energy", "forces"}:
            return "EnergyDipolesMACE", "energy_forces_dipole"
        return "AtomicDipolesMACE", "dipole"
    return None, None  # plain energy/forces (+stress): let the foundation decide


def _split(frames, valid_fraction, seed):
    """Deterministic train/valid split. At least one frame lands in each side."""
    gen = make_rng(seed)
    idx = np.arange(len(frames))
    gen.shuffle(idx)
    n_valid = max(1, int(round(len(frames) * valid_fraction))) if len(frames) > 1 else 0
    valid_idx = set(idx[:n_valid].tolist())
    train = [frames[i] for i in range(len(frames)) if i not in valid_idx]
    valid = [frames[i] for i in range(len(frames)) if i in valid_idx]
    return train, valid


def render_command(cfg, train_file: Path, valid_file: Path | None, out: Path) -> list[str]:
    """Render the full ``mace_run_train`` argv for ``cfg`` (pure; easy to test)."""
    base = os.environ.get(TRAIN_COMMAND_ENV) or DEFAULT_TRAIN_COMMAND
    argv = shlex.split(base)

    heads = set(cfg.heads)
    args: list[str] = [
        "--name", cfg.name,
        "--train_file", str(train_file),
        "--model_dir", str(out),  # the .model lands directly in the stage dir
        "--checkpoints_dir", str(out / "checkpoints"),
        "--results_dir", str(out / "results"),
        "--log_dir", str(out / "logs"),
        # explicit reference keys — match write_training_xyz so labels are read
        "--energy_key", REF_KEYS["energy"],
        "--forces_key", REF_KEYS["forces"],
        "--energy_weight", str(cfg.energy_weight),
        "--forces_weight", str(cfg.forces_weight),
        "--E0s", cfg.e0s,
        "--lr", str(cfg.lr),
        "--weight_decay", str(cfg.weight_decay),
        "--max_num_epochs", str(cfg.max_num_epochs),
        "--batch_size", str(cfg.batch_size),
        "--device", cfg.device,
        "--default_dtype", cfg.default_dtype,
    ]
    if valid_file is not None:
        args += ["--valid_file", str(valid_file)]
    if "stress" in heads:
        args += ["--stress_key", REF_KEYS["stress"], "--stress_weight", str(cfg.stress_weight)]
    if "dipole" in heads:
        args += ["--dipole_key", REF_KEYS["dipole"], "--dipole_weight", str(cfg.dipole_weight)]
    if "polarizability" in heads:
        args += [
            "--polarizability_key", REF_KEYS["polarizability"],
            "--polarizability_weight", str(cfg.polarizability_weight),
        ]

    model, loss = _model_and_loss(heads)
    if model is not None:
        args += ["--model", model]
    if loss is not None:
        args += ["--loss", loss]

    # foundation / fine-tuning strategy
    if cfg.strategy != "scratch":
        args += ["--foundation_model", cfg.foundation_model]
        multihead = cfg.strategy == "multihead"
        args += ["--multiheads_finetuning", "True" if multihead else "False"]
        if multihead:
            if cfg.pt_train_file:
                args += ["--pt_train_file", cfg.pt_train_file]
            args += [
                "--num_samples_pt", str(cfg.num_samples_pt),
                "--weight_pt", str(cfg.weight_pt),
                "--weight_ft", str(cfg.weight_ft),
            ]
    else:
        if cfg.hidden_irreps:
            args += ["--hidden_irreps", cfg.hidden_irreps]
        args += ["--r_max", str(cfg.r_max)]

    if cfg.ema:
        args += ["--ema", "--ema_decay", str(cfg.ema_decay)]
    if cfg.swa:
        args += ["--swa"]
    if cfg.seed is not None:
        args += ["--seed", str(cfg.seed)]

    # passthrough: arbitrary flags win over our defaults (de-dup the long name).
    for key, value in cfg.extra.items():
        flag = f"--{key}"
        if flag in args:
            pos = args.index(flag)
            # remove the old flag and its value (store_true flags have no value)
            if pos + 1 < len(args) and not args[pos + 1].startswith("--"):
                del args[pos : pos + 2]
            else:
                del args[pos]
        if isinstance(value, bool):
            args += [flag] if value else []
        else:
            args += [flag, str(value)]

    return argv + args


def _found_model(out: Path, name: str) -> Path | None:
    """Locate the produced ``.model`` (MACE writes ``<name>.model`` in model_dir)."""
    direct = out / f"{name}.model"
    if direct.exists():
        return direct
    candidates = sorted(out.glob("*.model"))
    return candidates[0] if candidates else None


@register("trainer", "mace")
def train_mace(frames, cfg, job, *, dry_run: bool = False) -> TrainResult:
    """Split, export, render, and run ``mace_run_train``; return a TrainResult."""
    if not frames:
        raise ValueError("training needs at least one labelled frame")
    out = Path(job.dir)
    out.mkdir(parents=True, exist_ok=True)

    train, valid = _split(frames, cfg.valid_fraction, cfg.seed)
    train_file = write_training_xyz(out / "train.xyz", train)
    valid_file = write_training_xyz(out / "valid.xyz", valid) if valid else None

    command = render_command(cfg, train_file, valid_file, out)
    manifest = {
        "backend": cfg.type,
        "name": cfg.name,
        "foundation_model": cfg.foundation_model,
        "strategy": cfg.strategy,
        "heads": list(cfg.heads),
        "e0s": cfg.e0s,
        "n_train": len(train),
        "n_valid": len(valid),
        "command": command,
        "reference": "Tompa et al., arXiv:2606.12704",
    }

    log_path = out / "train.log"
    if dry_run:
        manifest["dry_run"] = True
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
        logger.info("training (dry run): %s", " ".join(command))
        return TrainResult(
            model_path=None, n_train=len(train), n_valid=len(valid),
            heads=list(cfg.heads), command=command, manifest=manifest,
        )

    logger.info("training: %s", " ".join(command))
    t0 = time.perf_counter()
    with log_path.open("w") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    manifest["wall_seconds"] = round(time.perf_counter() - t0, 2)

    model_path = _found_model(out, cfg.name)
    manifest["model_path"] = str(model_path) if model_path else None
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    if model_path is None:
        logger.warning("training finished but no .model found under %s/model", out)

    return TrainResult(
        model_path=model_path, n_train=len(train), n_valid=len(valid),
        heads=list(cfg.heads), command=command, log_path=log_path, manifest=manifest,
    )
