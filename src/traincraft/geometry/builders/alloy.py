"""Shared *solid solution* machinery: turn an ordered lattice into a mixed solid.

A "mixed solid" (random substitutional alloy) is the lattice analogue of a
molecular mixture: instead of packing several species into a region, we replace
a fraction of the host lattice sites with other elements, chosen at random with
the builder's seed. The same ``ratio`` idea — and the same largest-remainder
apportionment — used for molecular mixtures applies here, so ``crystal``,
``slab`` and the ``surface_*`` builders can all grow an alloy with one helper.
"""

from __future__ import annotations

from ase import Atoms

from ...core import make_rng
from .mixture import apportion


def apply_solid_solution(atoms: Atoms, components: list, seed: int | None) -> tuple[Atoms, dict]:
    """Randomly substitute lattice sites per ``components`` (AlloyComponent list).

    Each component replaces ``round(ratio * n_sites)`` randomly chosen host sites
    with its element; ratios are apportioned together so they never overlap or
    exceed the lattice. Returns ``(new_atoms, composition)`` where ``composition``
    maps element -> substituted site count.
    """
    if not components:
        return atoms, {}
    total_ratio = sum(c.ratio for c in components)
    if total_ratio > 1.0 + 1e-9:
        raise ValueError(
            f"alloy composition ratios sum to {total_ratio:.3f} > 1 — leave room "
            "for the host element"
        )

    atoms = atoms.copy()
    n = len(atoms)
    # Apportion the substituted-site budget (ratio of all sites) across components.
    budget = int(round(total_ratio * n))
    weights = [c.ratio for c in components]
    n_sub = apportion(weights, budget)

    rng = make_rng(seed)
    sites = rng.permutation(n)
    cursor = 0
    composition: dict[str, int] = {}
    symbols = list(atoms.get_chemical_symbols())
    for c, k in zip(components, n_sub, strict=True):
        for idx in sites[cursor : cursor + k]:
            symbols[int(idx)] = c.element
        composition[c.element] = composition.get(c.element, 0) + int(k)
        cursor += k
    atoms.set_chemical_symbols(symbols)
    return atoms, composition
