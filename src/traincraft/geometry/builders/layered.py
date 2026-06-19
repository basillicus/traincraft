"""Builder: layered / 2D-stacked materials.

Stacks ``n_layers`` copies of a 2D monolayer (graphene, hBN, or an MX2 such as
MoS2) with a controllable interlayer spacing and stacking order:

  * ``stacking="AA"`` — layers sit directly on top of each other;
  * ``stacking="AB"`` — Bernal stacking; alternate layers are shifted by one
    nearest-neighbour vector ``(a1 + a2) / 3``;
  * ``twist`` (degrees) — each layer ``i`` is rotated by ``i * twist`` about the
    in-plane centroid, producing a moiré stack.

Periodicity note: an untwisted stack tiles its in-plane cell and is returned
periodic in xy.  A twisted stack is generally **incommensurate** with any small
cell, so a nonzero ``twist`` yields a finite, fully non-periodic flake padded
with vacuum — honest geometry with no boundary artefacts.
"""

from __future__ import annotations

import ase.build as ab
import numpy as np
from ase import Atoms

from ...core import Provenance, Structure, register


def _monolayer(cfg) -> Atoms:
    if cfg.material == "graphene":
        return ab.graphene(a=cfg.a) if cfg.a else ab.graphene()
    if cfg.material == "hbn":
        return ab.graphene(formula="BN", a=cfg.a or 2.5)
    if cfg.material == "mx2":
        return ab.mx2(formula=cfg.formula or "MoS2", a=cfg.a or 3.18)
    raise ValueError(f"unknown layered material {cfg.material!r}")


@register("builder", "layered")
def build_layered(cfg) -> Structure:
    base = _monolayer(cfg)
    nx, ny = cfg.size
    if (nx, ny) != (1, 1):
        base = base.repeat((nx, ny, 1))

    cell = np.asarray(base.get_cell())
    a1, a2 = cell[0], cell[1]
    pos0 = base.get_positions()
    pos0[:, 2] -= pos0[:, 2].min()  # baseline the layer at z = 0
    symbols = base.get_chemical_symbols()

    centroid = pos0[:, :2].mean(axis=0)
    ab_offset = ((a1 + a2) / 3.0)[:2]

    all_pos, all_sym = [], []
    for i in range(cfg.n_layers):
        p = pos0.copy()
        if cfg.twist:
            ang = np.deg2rad(i * cfg.twist)
            c, s = np.cos(ang), np.sin(ang)
            rot = np.array([[c, -s], [s, c]])
            p[:, :2] = (p[:, :2] - centroid) @ rot.T + centroid
        if cfg.stacking == "AB" and i % 2 == 1:
            p[:, :2] += ab_offset
        p[:, 2] += i * cfg.interlayer_spacing
        all_pos.append(p)
        all_sym += symbols

    positions = np.vstack(all_pos)

    if cfg.twist:
        atoms = Atoms(symbols=all_sym, positions=positions)
        atoms.set_pbc(False)
        atoms.center(vacuum=cfg.vacuum / 2)
    else:
        z_extent = float(positions[:, 2].max() - positions[:, 2].min())
        new_cell = np.array([a1, a2, [0.0, 0.0, z_extent + cfg.vacuum]])
        atoms = Atoms(
            symbols=all_sym, positions=positions, cell=new_cell, pbc=(True, True, False)
        )
        atoms.center(axis=2)

    return Structure.from_ase(
        atoms,
        provenance=Provenance(
            origin="generated",
            source=f"builder:layered:{cfg.material}-{cfg.n_layers}L",
            extra={
                "material": cfg.material,
                "n_layers": cfg.n_layers,
                "stacking": cfg.stacking,
                "twist": cfg.twist,
            },
        ),
    )
