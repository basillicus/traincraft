"""DFT labeling stage — the expensive step the selection funnel protects.

Runs *after* selection on the chosen frames: attaches a (DFT) calculator and
computes energy/forces/stress, plus dipole and polarizability when the calculator
declares them. Each labeled frame is tagged ``origin="dft_labeled"`` with its
level of theory, so the costly data stays cleanly separable from the cheap,
ML-generated points.

Each frame is computed in its own sub-directory (``frame_0000/`` …) so file-IO
calculators (FHI-aims, QE) don't clobber one another and reruns are resumable.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import numpy as np

from .calculators import make_calculator
from .core import Structure

logger = logging.getLogger(__name__)


def _level_of_theory(calc_cfg) -> dict:
    data = calc_cfg.model_dump()
    return {"calculator": data.pop("type"), **data}


def label_frames(frames: list[Structure], calc_cfg, *, out_dir=None) -> list[Structure]:
    """Label ``frames`` with ``calc_cfg``; return new ``dft_labeled`` Structures.

    ``out_dir`` (a directory) enables per-frame work dirs and writes a
    ``manifest.json`` recording the level of theory and the property set.
    """
    calc = make_calculator(calc_cfg)
    extra_props = list(getattr(calc_cfg, "properties", []) or [])
    level = _level_of_theory(calc_cfg)
    out_dir = Path(out_dir) if out_dir is not None else None

    labeled: list[Structure] = []
    t0 = time.perf_counter()
    for i, s in enumerate(frames):
        atoms = s.to_ase(with_properties=False)
        if out_dir is not None and hasattr(calc, "directory"):
            frame_dir = out_dir / f"frame_{i:04d}"
            frame_dir.mkdir(parents=True, exist_ok=True)
            calc.directory = str(frame_dir)
        atoms.calc = calc

        props: dict = {
            "energy": float(atoms.get_potential_energy()),
            "forces": np.asarray(atoms.get_forces()),
        }
        if atoms.cell.rank == 3 and atoms.get_pbc().any():
            try:
                props["stress"] = np.asarray(atoms.get_stress())
            except Exception:  # calc/system may not support stress
                logger.debug("calculator %s produced no stress", calc_cfg.type)
        if "dipole" in extra_props:
            props["dipole"] = np.asarray(atoms.get_dipole_moment())
        if "polarizability" in extra_props:
            pol = getattr(calc, "results", {}).get("polarizability")
            if pol is not None:
                props["polarizability"] = np.asarray(pol)

        new = s.copy()
        new.properties = props
        new.provenance.origin = "dft_labeled"
        new.provenance.calculator = calc_cfg.type
        new.provenance.level_of_theory = level
        labeled.append(new)
        logger.info("labeled %d/%d (E = %.4f eV)", i + 1, len(frames), props["energy"])

    if out_dir is not None:
        manifest = {
            "level_of_theory": level,
            "calculator": calc_cfg.type,
            "properties": ["energy", "forces", "stress", *extra_props],
            "n_frames": len(labeled),
            "wall_seconds": round(time.perf_counter() - t0, 2),
            "origin": "dft_labeled",
        }
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    return labeled
