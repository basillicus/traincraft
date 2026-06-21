"""Builders: molecule(s) on a crystalline surface.

Two builders:
  surface_adsorbate  — single adsorbate via ase.build.add_adsorbate (zero extra deps)
  surface_packing    — N-molecule (or mixed) coverage via the Packmol binary

Either may grow a *mixed-solid* (alloy) slab via ``composition`` (see
:mod:`.alloy`); ``surface_packing`` may pack a *mixture* of adsorbate species
(see :mod:`.mixture`).

Fragment convention:
  tc_fragment == -1   substrate atoms (immobile)
  tc_fragment >= 0    each adsorbate molecule is a separate mobile fragment
"""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.build import add_adsorbate

from ...config.models import Species
from ...core import Provenance, Structure, register, set_fragments

# re-exported for back-compat; canonical home is now _adsorbate.py
from ._adsorbate import _canonical_smiles, _resolve_adsorbate  # noqa: F401
from ._common import build_named_slab
from .alloy import apply_solid_solution
from .mixture import resolve_mixture, run_packmol, tag_mixture


def _build_slab(element: str, facet: str, size: tuple, vacuum: float) -> Atoms:
    return build_named_slab(element, facet, size, vacuum)


def _adsorbate_species(cfg) -> list[Species]:
    """The single-adsorbate shortcut becomes a one-element mixture."""
    if getattr(cfg, "species", None):
        return list(cfg.species)
    return [
        Species(
            molecule_name=cfg.molecule_name,
            smiles=cfg.smiles,
            file=cfg.file,
            count=getattr(cfg, "n_molecules", 1),
        )
    ]


@register("builder", "surface_adsorbate")
def build_surface_adsorbate(cfg) -> Structure:
    slab = _build_slab(cfg.element, cfg.facet, cfg.size, cfg.vacuum)
    slab, composition = apply_solid_solution(slab, cfg.composition, cfg.seed)
    n_slab = len(slab)

    ads = _resolve_adsorbate(cfg.molecule_name, cfg.smiles, cfg.file)
    n_ads = len(ads)

    kwargs = {"height": cfg.height, "position": cfg.site}
    if cfg.offset is not None:
        kwargs["offset"] = cfg.offset
    add_adsorbate(slab, ads, **kwargs)
    slab.info.pop("adsorbate_info", None)  # not extxyz-serialisable

    # Fragment tagging: substrate=-1, adsorbate=0
    frag = np.full(n_slab + n_ads, -1, dtype=int)
    frag[n_slab:] = 0
    set_fragments(slab, frag)

    extra: dict = {}
    if composition:
        extra["composition"] = composition
    if cfg.smiles is not None:
        canonical = _canonical_smiles(cfg.smiles)
        extra["smiles"] = canonical
        extra["fragment_smiles"] = {"0": canonical}

    return Structure.from_ase(
        slab,
        provenance=Provenance(
            origin="generated",
            source=f"builder:surface_adsorbate:{cfg.element}-{cfg.facet}",
            extra=extra,
        ),
    )


@register("builder", "surface_packing")
def build_surface_packing(cfg) -> Structure:
    slab = _build_slab(cfg.element, cfg.facet, cfg.size, cfg.vacuum)
    slab, composition = apply_solid_solution(slab, cfg.composition, cfg.seed)
    n_slab = len(slab)

    resolved = resolve_mixture(_adsorbate_species(cfg), cfg.n_molecules)

    # Determine packing region: above slab top, within the xy cell footprint.
    cell = slab.get_cell()
    slab_top = float(slab.get_positions()[:, 2].max())
    z_lo = slab_top + cfg.gap
    z_hi = z_lo + cfg.region_height
    x_hi = float(np.linalg.norm(cell[0]))
    y_hi = float(np.linalg.norm(cell[1]))
    region = f"inside box 0.0 0.0 {z_lo:.4f} {x_hi:.4f} {y_hi:.4f} {z_hi:.4f}"

    seed = cfg.seed if cfg.seed is not None else 12345
    items = [(atoms, n) for atoms, n, _, _ in resolved]
    packed = run_packmol(items, region, cfg.tolerance, seed, need="surface_packing")

    combined = slab + packed

    frag = np.full(len(combined), -1, dtype=int)
    n_frag, fragment_smiles, fragment_species, _ = tag_mixture(frag, resolved, n_slab)
    set_fragments(combined, frag)

    extra: dict = {"n_molecules": int(n_frag)}
    if composition:
        extra["composition"] = composition
    if fragment_species:
        extra["fragment_species"] = fragment_species
    if fragment_smiles:
        extra["fragment_smiles"] = fragment_smiles

    return Structure.from_ase(
        combined,
        provenance=Provenance(
            origin="generated",
            source=f"builder:surface_packing:{cfg.element}-{cfg.facet}",
            extra=extra,
        ),
    )
