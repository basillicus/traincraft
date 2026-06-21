"""Builder: a carbon nanotube randomly filled with molecules ("fillMyTubes").

The original driver for this project: take a CNT and stuff molecules inside it,
then sample/label to study (e.g.) confined-phase Raman. Packmol packs the guests
into a cylinder coaxial with the tube and just narrower than its wall; the tube
carbons are tagged framework (``tc_fragment == -1``) and every guest molecule
gets its own fragment id, so the MC sampler and the ``constraints`` transform can
address tube and guests separately.

The guest may be a single molecule (``molecule_name``/``smiles``/``file``) or a
*mixture* (a ``species`` list with counts or ratios) — the same species/packing
machinery used by the ``liquid`` and ``surface_packing`` builders (see
:mod:`.mixture`). Only the Packmol *region* differs (``inside cylinder``).
"""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.build import nanotube

from ...config.models import Species
from ...core import FRAMEWORK, Provenance, Structure, register, set_fragments
from .mixture import resolve_mixture, run_packmol, tag_mixture


def _tube_axis_and_radius(cnt: Atoms) -> tuple[float, float, float]:
    """Return (cx, cy, R): the tube's xy axis centre and mean carbon radius."""
    pos = cnt.get_positions()
    cx, cy = pos[:, 0].mean(), pos[:, 1].mean()
    radii = np.hypot(pos[:, 0] - cx, pos[:, 1] - cy)
    return float(cx), float(cy), float(radii.mean())


def _guest_species(cfg) -> list[Species]:
    """The single-molecule shortcut becomes a one-element mixture."""
    if cfg.species:
        return list(cfg.species)
    return [
        Species(
            molecule_name=cfg.molecule_name,
            smiles=cfg.smiles,
            file=cfg.file,
            count=cfg.n_molecules,
        )
    ]


@register("builder", "filled_nanotube")
def build_filled_nanotube(cfg) -> Structure:
    cnt = nanotube(cfg.n, cfg.m, length=cfg.length, bond=cfg.bond)
    cnt.center(vacuum=cfg.vacuum / 2, axis=(0, 1))

    cx, cy, radius = _tube_axis_and_radius(cnt)
    length_z = float(cnt.get_cell()[2][2])
    pack_radius = radius - cfg.radial_margin
    pack_length = length_z - 2 * cfg.axial_margin
    if pack_radius <= 0 or pack_length <= 0:
        raise ValueError(
            f"tube too small to fill: usable radius {pack_radius:.2f} Å, length "
            f"{pack_length:.2f} Å. Widen the tube (n/m) or lengthen it (length), "
            "or reduce radial_margin/axial_margin."
        )

    resolved = resolve_mixture(_guest_species(cfg), cfg.n_molecules)
    seed = cfg.seed if cfg.seed is not None else 12345
    # inside cylinder a1 a2 a3 d1 d2 d3 r l — base centre, axis direction, radius, length
    region = (
        f"inside cylinder {cx:.4f} {cy:.4f} {cfg.axial_margin:.4f} 0. 0. 1. "
        f"{pack_radius:.4f} {pack_length:.4f}"
    )
    items = [(atoms, n) for atoms, n, _, _ in resolved]
    guests = run_packmol(items, region, cfg.tolerance, seed, need="filled_nanotube")

    # Combine: tube first (framework), then the packed guests (one fragment each).
    system = cnt + guests
    system.set_cell(cnt.get_cell())
    system.set_pbc((False, False, cfg.pbc))

    frag = np.empty(len(system), dtype=int)
    frag[: len(cnt)] = FRAMEWORK
    n_frag, fragment_smiles, fragment_species, _ = tag_mixture(frag, resolved, len(cnt))
    set_fragments(system, frag)

    extra: dict = {
        "n_molecules": int(n_frag),
        "tube_radius": round(radius, 4),
        "pack_radius": round(pack_radius, 4),
    }
    if fragment_species:
        extra["fragment_species"] = fragment_species
    if fragment_smiles:
        extra["fragment_smiles"] = fragment_smiles

    return Structure.from_ase(
        system,
        provenance=Provenance(
            origin="generated",
            source=f"builder:filled_nanotube:{cfg.n}-{cfg.m}-l{cfg.length}",
            extra=extra,
        ),
    )
