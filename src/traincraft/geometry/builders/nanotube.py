"""Builder: carbon nanotube (ported from legacy ``gengeom.create_nanotube``)."""

from __future__ import annotations

from ase.build import nanotube

from ...core import Provenance, Structure, register


@register("builder", "nanotube")
def build_nanotube(cfg) -> Structure:
    cnt = nanotube(cfg.n, cfg.m, length=cfg.length, bond=cfg.bond)
    cnt.center(vacuum=cfg.vacuum / 2, axis=(0, 1))
    return Structure.from_ase(
        cnt,
        provenance=Provenance(
            origin="generated", source=f"builder:nanotube:{cfg.n}-{cfg.m}-l{cfg.length}"
        ),
    )
