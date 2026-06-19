"""Transform: fix atoms with an ASE ``FixAtoms`` constraint.

Selectors are OR-ed together into one set of atom indices, then frozen:

  * ``indices``    — explicit atom indices;
  * ``elements``   — every atom of these chemical symbols (e.g. ``["Cu"]``);
  * ``fragments``  — every atom with these ``tc_fragment`` ids (``-1`` = substrate);
  * ``below_z``    — every atom with a Cartesian *z* below this value (Å) — the
    usual way to freeze the bottom layers of a slab.

This deliberately selects on the *final* structure, so it is the correct place
to (re)apply constraints **after** a builder that reorders atoms (e.g. Packmol),
fixing the legacy bug where index-based constraints were applied before packing
and then pointed at the wrong atoms.
"""

from __future__ import annotations

import numpy as np
from ase.constraints import FixAtoms

from ...core import Structure, register


def _selected_indices(atoms, cfg) -> set[int]:
    n = len(atoms)
    chosen: set[int] = set()

    if cfg.indices:
        for i in cfg.indices:
            if not 0 <= i < n:
                raise IndexError(f"constraint index {i} out of range (0..{n - 1})")
            chosen.add(int(i))

    if cfg.elements:
        wanted = set(cfg.elements)
        chosen.update(i for i, s in enumerate(atoms.get_chemical_symbols()) if s in wanted)

    if cfg.fragments is not None:
        frag = atoms.arrays.get("tc_fragment")
        if frag is None:
            raise ValueError("constraints: 'fragments' given but structure has no tc_fragment")
        wanted = set(cfg.fragments)
        chosen.update(i for i, f in enumerate(frag) if int(f) in wanted)

    if cfg.below_z is not None:
        z = atoms.get_positions()[:, 2]
        chosen.update(int(i) for i in np.nonzero(z < cfg.below_z)[0])

    return chosen


@register("transform", "constraints")
def transform_constraints(structure: Structure, cfg) -> Structure:
    out = structure.copy()
    idx = _selected_indices(out.atoms, cfg)
    if not idx:
        raise ValueError(
            "constraints transform selected no atoms; set one of "
            "'indices', 'elements', 'fragments', or 'below_z'"
        )

    out.atoms.set_constraint(FixAtoms(indices=sorted(idx)))
    out.provenance.transforms.append(f"constraints:fixed={len(idx)}")
    return out
