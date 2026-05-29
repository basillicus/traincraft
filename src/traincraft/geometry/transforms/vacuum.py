"""Transform: add vacuum padding around the structure."""

from __future__ import annotations

from ...core import Structure, register


@register("transform", "vacuum")
def transform_vacuum(structure: Structure, cfg) -> Structure:
    out = structure.copy()
    out.atoms.center(vacuum=cfg.amount / 2)
    out.provenance.transforms.append("vacuum")
    return out
