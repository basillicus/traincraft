"""Extended-XYZ IO with provenance + properties.

Properties live in ``info``/``arrays`` with a ``tc_`` prefix; provenance is
JSON-encoded so the labeling origin and level-of-theory survive a round-trip.
(Reimplements the useful part of the legacy ``scripts/qe2extxyz.py``.)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.io import read, write

from ..core import Provenance, Structure

_META_KEYS = {"tc_provenance", "tc_hash"}


def _to_atoms(s: Structure):
    atoms = s.atoms.copy()
    atoms.info["tc_provenance"] = json.dumps(s.provenance.to_dict())
    atoms.info["tc_hash"] = s.hash
    for key, value in s.properties.items():
        if value is None:
            continue
        if key == "forces":
            atoms.arrays["tc_forces"] = np.asarray(value)
        elif isinstance(value, (int, float, str)):
            atoms.info[f"tc_{key}"] = value
        else:
            atoms.info[f"tc_{key}"] = json.dumps(np.asarray(value).tolist())
    return atoms


def write_frames(path: str | Path, structures: list[Structure], append: bool = False) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    images = [_to_atoms(s) for s in structures]
    write(str(path), images, format="extxyz", append=append)
    return path


def read_frames(path: str | Path) -> list[Structure]:
    images = read(str(path), index=":", format="extxyz")
    if not isinstance(images, list):
        images = [images]
    out: list[Structure] = []
    for atoms in images:
        props: dict = {}
        if "tc_forces" in atoms.arrays:
            props["forces"] = np.asarray(atoms.arrays["tc_forces"])
        for key in list(atoms.info):
            if key.startswith("tc_") and key not in _META_KEYS:
                props[key[3:]] = atoms.info[key]
        raw = atoms.info.get("tc_provenance")
        prov = Provenance.from_dict(json.loads(raw)) if raw else Provenance()
        out.append(Structure(atoms=atoms, properties=props, provenance=prov))
    return out
