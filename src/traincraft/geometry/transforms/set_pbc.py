"""Transform: set periodic boundary conditions.

``pbc`` is either a single bool (applied to all three axes) or a 3-tuple of
bools. This only flips the periodicity flags; it does not move atoms or change
the cell. Pair it with ``vacuum`` if you need padding after going non-periodic.
"""

from __future__ import annotations

from ...core import Structure, register


@register("transform", "set_pbc")
def transform_set_pbc(structure: Structure, cfg) -> Structure:
    pbc = cfg.pbc if isinstance(cfg.pbc, bool) else tuple(cfg.pbc)
    out = structure.copy()
    out.atoms.set_pbc(pbc)
    out.provenance.transforms.append(f"set_pbc:{pbc}")
    return out
