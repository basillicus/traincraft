"""Shared slab/surface helpers used by the ``slab`` and ``surface_*`` builders."""

from __future__ import annotations

import ase.build as ab
from ase import Atoms

# Map config facet strings to ase.build functions (named low-index facets).
FACET_FN = {
    "fcc111": "fcc111",
    "fcc100": "fcc100",
    "fcc110": "fcc110",
    "bcc110": "bcc110",
    "bcc100": "bcc100",
    "bcc111": "bcc111",
    "hcp0001": "hcp0001",
}


def build_named_slab(
    element: str,
    facet: str,
    size: tuple[int, int, int],
    vacuum: float,
    *,
    orthogonal: bool = False,
) -> Atoms:
    """Build a slab from a named low-index facet (e.g. ``fcc111``).

    ``vacuum`` is the *total* vacuum thickness; ASE adds half on each side.
    """
    fn = getattr(ab, FACET_FN[facet])
    return fn(element, size=size, vacuum=vacuum / 2, orthogonal=orthogonal)


def build_miller_slab(
    element: str,
    crystalstructure: str | None,
    miller: tuple[int, int, int],
    layers: int,
    vacuum: float,
    *,
    a: float | None = None,
    c: float | None = None,
    cubic: bool = False,
    size: tuple[int, int] = (1, 1),
    periodic: bool = False,
) -> Atoms:
    """Build a slab for arbitrary Miller indices by cleaving a bulk crystal.

    ``vacuum`` is the *total* vacuum thickness; ASE's ``surface`` adds half on
    each side.  ``size`` repeats the resulting slab in-plane.
    """
    lattice = ab.bulk(element, crystalstructure, a=a, c=c, cubic=cubic)
    slab = ab.surface(lattice, miller, layers, vacuum=vacuum / 2, periodic=periodic)
    if size != (1, 1):
        slab = slab.repeat((size[0], size[1], 1))
    return slab
