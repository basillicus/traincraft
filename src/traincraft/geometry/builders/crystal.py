"""Builder: bulk crystal (optionally with point defects).

Builds a periodic bulk cell via ``ase.build.bulk``, optionally repeats it into a
supercell, then applies point defects:

  * **vacancy**       — remove the atom at ``index``;
  * **substitution**  — swap the atom at ``index`` for ``element``;
  * **interstitial**  — insert an atom of ``element`` at ``position``.

Defects are applied in a stable order (substitutions, then vacancies, then
interstitials) so that the ``index`` fields always refer to the supercell
*before* any atoms are added or removed — making a list of defects predictable.
"""

from __future__ import annotations

import numpy as np
from ase.build import bulk

from ...core import Provenance, Structure, register


def _apply_defects(atoms, defects):
    """Return a new Atoms with the requested point defects applied."""
    atoms = atoms.copy()
    n = len(atoms)

    sub = [d for d in defects if d.kind == "substitution"]
    vac = [d for d in defects if d.kind == "vacancy"]
    inter = [d for d in defects if d.kind == "interstitial"]

    for d in sub:
        idx = d.index if d.index is not None else 0
        if not 0 <= idx < n:
            raise IndexError(f"substitution index {idx} out of range (0..{n - 1})")
        atoms[idx].symbol = d.element

    # Remove vacancies together so indices stay consistent with the supercell.
    if vac:
        drop = set()
        for d in vac:
            idx = d.index if d.index is not None else 0
            if not 0 <= idx < n:
                raise IndexError(f"vacancy index {idx} out of range (0..{n - 1})")
            drop.add(idx)
        keep = [i for i in range(n) if i not in drop]
        atoms = atoms[keep]

    for d in inter:
        pos = np.asarray(d.position, dtype=float)
        if not d.cartesian:
            pos = pos @ np.asarray(atoms.get_cell())  # fractional -> cartesian
        atoms += type(atoms)(d.element, positions=[pos], cell=atoms.get_cell(),
                             pbc=atoms.get_pbc())
    return atoms


@register("builder", "crystal")
def build_crystal(cfg) -> Structure:
    atoms = bulk(
        cfg.name,
        cfg.crystalstructure,
        a=cfg.a,
        b=cfg.b,
        c=cfg.c,
        cubic=cfg.cubic,
        orthorhombic=cfg.orthorhombic,
    )
    if tuple(cfg.supercell) != (1, 1, 1):
        atoms = atoms.repeat(tuple(cfg.supercell))

    extra: dict = {}
    if cfg.defects:
        atoms = _apply_defects(atoms, cfg.defects)
        extra["defects"] = [
            {"kind": d.kind, "index": d.index, "element": d.element}
            for d in cfg.defects
        ]

    return Structure.from_ase(
        atoms,
        provenance=Provenance(
            origin="generated",
            source=f"builder:crystal:{cfg.name}",
            extra=extra,
        ),
    )
