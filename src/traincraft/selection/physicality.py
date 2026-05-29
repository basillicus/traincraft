"""Selector: drop unphysical frames (atoms too close)."""

from __future__ import annotations

import numpy as np

from ..core import Structure, register


@register("selector", "physicality")
def select_physicality(frames: list[Structure], cfg) -> list[Structure]:
    out: list[Structure] = []
    for s in frames:
        atoms = s.atoms
        n = len(atoms)
        if n < 2:
            out.append(s)
            continue
        d = atoms.get_all_distances(mic=bool(any(atoms.pbc)))
        iu = np.triu_indices(n, k=1)
        if d[iu].min() >= cfg.min_distance:
            out.append(s)
    return out
