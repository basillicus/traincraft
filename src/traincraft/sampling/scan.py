"""Sampler: deterministic grid scan of one fragment's position/orientation.

Unlike ``md`` / ``monte_carlo`` / ``rattle`` (which sample *stochastically*),
``scan`` enumerates a regular **Cartesian product** of translations and/or
rotations applied to a chosen fragment, producing an exact grid of candidates —
ideal for potential-energy-surface scans, binding/dissociation curves and
orientation grids that MD/MC sample inefficiently.

The fragment to move is selected by ``tc_fragment`` id (the adsorbate/guest is
fragment 0; the substrate/tube is framework ``-1`` and never moves). Energies are
*not* computed here — the labeling stage attaches them, exactly as for the other
samplers.

!!! tip
    A scan is meant to keep *every* grid point, so pair it with a generous
    ``[selection] budget`` (or drop the ``diversity`` step) — otherwise the
    funnel will thin your grid.
"""

from __future__ import annotations

import itertools
import logging

import numpy as np

from ..core import Job, Provenance, Structure, register
from ..core.fragments import get_fragments

logger = logging.getLogger(__name__)

_AXES = {
    "x": np.array([1.0, 0.0, 0.0]),
    "y": np.array([0.0, 1.0, 0.0]),
    "z": np.array([0.0, 0.0, 1.0]),
}


def _unit(axis) -> np.ndarray:
    """Resolve 'x'/'y'/'z' or an explicit vector to a unit vector."""
    v = _AXES[axis].copy() if isinstance(axis, str) else np.asarray(axis, dtype=float)
    n = float(np.linalg.norm(v))
    if n == 0.0:
        raise ValueError("scan axis vector must be non-zero")
    return v / n


def _grid(spec) -> list[float]:
    """`steps` values evenly spaced over [start, stop] (inclusive of both ends)."""
    return [float(x) for x in np.linspace(spec.start, spec.stop, spec.steps)]


def _fragment_indices(structure: Structure, fragment) -> np.ndarray:
    """Indices of the atoms the scan moves (a fragment id, or every atom)."""
    atoms = structure.atoms
    n = len(atoms)
    if fragment == "all":
        return np.arange(n)
    frag = get_fragments(atoms)
    if frag is None:
        logger.info("scan: no fragments set on the structure; moving the whole system")
        return np.arange(n)
    mask = frag == int(fragment)
    if not mask.any():
        available = sorted({int(i) for i in frag if i != -1})
        raise ValueError(
            f"scan: fragment {fragment!r} has no atoms; available mobile fragments: "
            f"{available}. Use a builder that tags fragments, or set fragment = \"all\"."
        )
    return np.where(mask)[0]


@register("sampler", "scan")
def sample_scan(structure: Structure, calc, job: Job, cfg) -> list[Structure]:
    from scipy.spatial.transform import Rotation

    idx = _fragment_indices(structure, cfg.fragment)

    t_vals = _grid(cfg.translate) if cfg.translate is not None else [None]
    t_axis = _unit(cfg.translate.axis) if cfg.translate is not None else None
    r_vals = _grid(cfg.rotate) if cfg.rotate is not None else [None]
    r_axis = _unit(cfg.rotate.axis) if cfg.rotate is not None else None

    frames: list[Structure] = []
    for t, r in itertools.product(t_vals, r_vals):
        atoms = structure.atoms.copy()  # preserves tc_fragment (ASE per-atom array)
        pos = atoms.get_positions()
        if r is not None:
            centroid = pos[idx].mean(axis=0)
            rot = Rotation.from_rotvec(np.deg2rad(r) * r_axis)
            pos[idx] = rot.apply(pos[idx] - centroid) + centroid
        if t is not None:
            pos[idx] = pos[idx] + t * t_axis
        atoms.set_positions(pos)

        tag = ",".join(
            part for part in (
                f"t={t:.3f}" if t is not None else None,
                f"r={r:.1f}" if r is not None else None,
            ) if part is not None
        )
        frames.append(
            Structure.from_ase(
                atoms,
                provenance=Provenance(
                    origin="ml_sampled",
                    source=f"sampler:scan:{tag}",
                    parents=[structure.hash],
                ),
            )
        )

    logger.info(
        "scan: generated %d grid candidate(s) (fragment=%s, translate=%s, rotate=%s)",
        len(frames), cfg.fragment, cfg.translate is not None, cfg.rotate is not None,
    )
    return frames
