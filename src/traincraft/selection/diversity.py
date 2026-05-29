"""Selector: farthest-point sampling for diversity.

Uses a fixed-length pairwise-distance histogram descriptor (works for any atom
count). SOAP/ACE descriptors via the ``descriptors`` extra come later.
"""

from __future__ import annotations

import numpy as np

from ..core import Structure, register


def _descriptor(atoms, nbins: int = 24, rcut: float = 8.0) -> np.ndarray:
    n = len(atoms)
    if n < 2:
        return np.zeros(nbins)
    d = atoms.get_all_distances(mic=bool(any(atoms.pbc)))
    iu = np.triu_indices(n, k=1)
    dd = d[iu]
    hist, _ = np.histogram(dd, bins=nbins, range=(0.0, rcut))
    return hist.astype(float) / max(len(dd), 1)


@register("selector", "diversity")
def select_diversity(frames: list[Structure], cfg) -> list[Structure]:
    k = cfg.budget
    if k is None or k >= len(frames) or len(frames) == 0:
        return frames
    x = np.array([_descriptor(s.atoms) for s in frames])
    chosen = [0]
    dist = np.linalg.norm(x - x[0], axis=1)
    for _ in range(1, k):
        idx = int(np.argmax(dist))
        chosen.append(idx)
        dist = np.minimum(dist, np.linalg.norm(x - x[idx], axis=1))
    return [frames[i] for i in chosen]
