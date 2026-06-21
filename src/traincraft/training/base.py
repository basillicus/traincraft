"""Training subsystem: turn a labelled dataset into a trained MLIP.

A ``trainer`` is a registry plugin (``@register("trainer", name)``) that takes a
list of labelled :class:`~traincraft.core.Structure` frames plus a config and a
job directory, and returns a :class:`TrainResult` (model path + metrics +
manifest). The model interface is pluggable on purpose — adding a backend
(MatterSim/Orb/SevenNet/…) is one new file and one registry entry, exactly like
calculators and samplers (DESIGN §12).

This module owns two things every backend shares:

* :func:`run_training` — resolve the trainer by ``cfg.type`` and run it.
* :func:`write_training_xyz` — export frames to an extended-XYZ the trainer can
  read, **re-keying** TrainCraft's ``tc_*`` properties onto explicit reference
  keys (``REF_energy``/``REF_forces``/…). The trainer then tells the engine which
  keys to read, so labels never get silently dropped on the floor.

Heavy backends (torch, mace) are imported lazily inside the trainer factories,
mirroring ``calculators/potentials.py`` — the package imports with no ML stack.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ase.io import write

from ..core import Structure, get

# Reference-property keys written into the training extxyz. These are what the
# trainer passes to the backend as --energy_key/--forces_key/etc., so the schema
# is owned here in one place rather than guessed per backend.
REF_KEYS = {
    "energy": "REF_energy",
    "forces": "REF_forces",
    "stress": "REF_stress",
    "dipole": "REF_dipole",
    "polarizability": "REF_polarizability",
}


@dataclass
class TrainResult:
    """Outcome of a training run."""

    model_path: Path | None  # the trained .model (None for a dry run)
    n_train: int
    n_valid: int
    heads: list[str]
    command: list[str]  # the rendered backend command (argv)
    log_path: Path | None = None
    manifest: dict = field(default_factory=dict)


def _as_float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_array(value) -> np.ndarray | None:
    """Coerce a stored property (ndarray, list, or JSON string) to an array."""
    if value is None:
        return None
    if isinstance(value, str):
        import json

        try:
            value = json.loads(value)
        except json.JSONDecodeError:
            return None
    return np.asarray(value, dtype=float)


def write_training_xyz(path: str | Path, frames: list[Structure]) -> Path:
    """Write ``frames`` to extended-XYZ with explicit ``REF_*`` property keys.

    Energy/dipole/polarizability go to ``atoms.info``; forces to ``atoms.arrays``;
    stress is flattened to its 9 components (3×3 or Voigt-6 both accepted). Frames
    missing a label for a given property simply omit that key — MACE skips them
    per head rather than training on zeros.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    images = []
    for s in frames:
        atoms = s.atoms.copy()
        props = s.properties
        energy = _as_float(props.get("energy"))
        if energy is not None:
            atoms.info[REF_KEYS["energy"]] = energy
        forces = _as_array(props.get("forces"))
        if forces is not None and forces.shape == (len(atoms), 3):
            atoms.arrays[REF_KEYS["forces"]] = forces
        stress = _as_array(props.get("stress"))
        if stress is not None:
            atoms.info[REF_KEYS["stress"]] = stress.reshape(-1).tolist()
        dipole = _as_array(props.get("dipole"))
        if dipole is not None:
            atoms.info[REF_KEYS["dipole"]] = dipole.reshape(-1).tolist()
        pol = _as_array(props.get("polarizability"))
        if pol is not None:
            atoms.info[REF_KEYS["polarizability"]] = pol.reshape(-1).tolist()
        images.append(atoms)
    write(str(path), images, format="extxyz")
    return path


def run_training(frames: list[Structure], cfg, job, *, dry_run: bool = False) -> TrainResult:
    """Resolve the trainer for ``cfg.type`` and train on ``frames``.

    ``job`` is a :class:`~traincraft.core.Workspace`-style object exposing a
    ``dir`` attribute (the run directory). ``dry_run`` renders the command and
    writes the train/valid splits without launching the backend — used by tests
    and ``submit --dry-run``.
    """
    trainer = get("trainer", cfg.type)
    return trainer(frames, cfg, job, dry_run=dry_run)
