"""Builder: intercalate guest atoms into the galleries of a layered host.

Builds a planar layered host (graphene or hBN — see note) via the ``layered``
builder, then inserts ``n_per_gallery`` guest atoms on an in-plane grid at the
mid-plane of each selected gallery. ``stage`` controls staging: guests go into
galleries whose index is a multiple of ``stage`` (stage-1 fills every gallery,
stage-2 every other, ...), matching the electrochemical staging convention.

The galleries are detected from the host's distinct layer z-planes, so the host
must have planar layers; ``mx2`` (three sub-planes per layer) is rejected. The
gallery may be opened up by ``gallery_expansion`` (Å added to the host spacing
around occupied galleries) to leave room for the guest before relaxation.
"""

from __future__ import annotations

import math

import numpy as np

from ...core import Provenance, Structure, register, set_fragments
from .layered import build_layered


def _grid_fracs(n: int) -> list[tuple[float, float]]:
    """`n` in-plane fractional points spread over a near-square grid in (0,1)."""
    side = math.ceil(math.sqrt(n))
    fracs = []
    for k in range(n):
        i, j = divmod(k, side)
        fracs.append(((i + 0.5) / side, (j + 0.5) / side))
    return fracs


@register("builder", "intercalation")
def build_intercalation(cfg) -> Structure:
    if cfg.host.material == "mx2":
        raise ValueError("intercalation host must have planar layers; 'mx2' is not supported")
    if cfg.host.twist:
        raise ValueError("intercalation needs a periodic (untwisted) host")

    host = build_layered(cfg.host)
    atoms = host.to_ase()
    n_host = len(atoms)
    cell = np.asarray(atoms.get_cell())
    a1, a2 = cell[0], cell[1]

    # Distinct layer z-planes (host is planar), low -> high.
    z = atoms.get_positions()[:, 2]
    layer_z = np.array(sorted({round(float(v), 2) for v in z}))
    if len(layer_z) != cfg.host.n_layers:
        raise RuntimeError(
            f"expected {cfg.host.n_layers} layer planes, found {len(layer_z)}; "
            "host is not planar as required"
        )

    galleries = [(layer_z[i] + layer_z[i + 1]) / 2.0 for i in range(len(layer_z) - 1)]
    targets = [g for k, g in enumerate(galleries) if k % cfg.stage == 0]

    fracs = _grid_fracs(cfg.n_per_gallery)
    guest_pos = []
    for z_mid in targets:
        for f1, f2 in fracs:
            xy = f1 * a1[:2] + f2 * a2[:2]
            guest_pos.append([xy[0], xy[1], z_mid])

    if guest_pos:
        guest = type(atoms)(
            [cfg.guest] * len(guest_pos),
            positions=np.asarray(guest_pos),
            cell=cell,
            pbc=atoms.get_pbc(),
        )
        atoms = atoms + guest

    # Optionally open the cell along c to make room for the guests.
    if cfg.gallery_expansion:
        new_cell = cell.copy()
        new_cell[2, 2] += cfg.gallery_expansion * len(targets)
        atoms.set_cell(new_cell, scale_atoms=False)

    # Fragment tagging: host = -1, all guests = 0 (a single mobile species).
    frag = np.full(len(atoms), -1, dtype=int)
    frag[n_host:] = 0
    set_fragments(atoms, frag)

    return Structure.from_ase(
        atoms,
        provenance=Provenance(
            origin="generated",
            source=f"builder:intercalation:{cfg.guest}@{cfg.host.material}",
            extra={
                "guest": cfg.guest,
                "n_guests": len(guest_pos),
                "stage": cfg.stage,
                "galleries_filled": len(targets),
            },
        ),
    )
