"""Builder: liquids, mixtures, and confined bulk via Packmol.

Packs one or more molecular :class:`Species` into a periodic box. The box is
either given explicitly (``box = [a, b, c]`` in Å) or derived from a target mass
``density`` (g/cm³) for a cubic cell. Amounts may be absolute (``count``) or
relative (``ratio`` + a total ``n_molecules``). Each packed molecule becomes its
own ``tc_fragment`` (numbered globally across species), so the MC sampler and
constraints can address them individually.

The species/packing/tagging machinery is shared with the ``surface_packing`` and
``filled_nanotube`` builders (see :mod:`.mixture`).
"""

from __future__ import annotations

import numpy as np

from ...core import Provenance, Structure, register, set_fragments
from .mixture import resolve_mixture, run_packmol, tag_mixture

_AMU_G = 1.66053906660e-24  # grams per atomic mass unit


def _cubic_edge_from_density(total_mass_amu: float, density: float) -> float:
    """Cubic box edge (Å) holding ``total_mass_amu`` at ``density`` g/cm³."""
    volume_cm3 = (total_mass_amu * _AMU_G) / density
    edge_cm = volume_cm3 ** (1.0 / 3.0)
    return edge_cm * 1.0e8  # cm -> Å


@register("builder", "liquid")
def build_liquid(cfg) -> Structure:
    # Resolve the mixture: [(atoms, count, label, canonical_smiles), ...].
    resolved = resolve_mixture(cfg.species, cfg.n_molecules)

    # Box: explicit, or cubic from target density.
    if cfg.box is not None:
        box = tuple(float(x) for x in cfg.box)
    elif cfg.density is not None:
        total_mass = sum(float(a.get_masses().sum()) * n for a, n, _, _ in resolved)
        edge = _cubic_edge_from_density(total_mass, cfg.density)
        box = (edge, edge, edge)
    else:
        raise ValueError("liquid builder needs either 'box' or 'density'")

    seed = cfg.seed if cfg.seed is not None else 12345
    # Inset packing by `tolerance` from each wall so molecules don't clash across
    # the periodic boundary (Packmol itself is not PBC-aware).
    lo, m = cfg.tolerance, cfg.tolerance
    region = (
        f"inside box {lo:.4f} {lo:.4f} {lo:.4f} "
        f"{box[0] - m:.4f} {box[1] - m:.4f} {box[2] - m:.4f}"
    )
    items = [(atoms, n) for atoms, n, _, _ in resolved]
    packed = run_packmol(items, region, cfg.tolerance, seed, need="liquid builder")

    packed.set_cell(box)
    packed.set_pbc(cfg.pbc)

    frag = np.empty(len(packed), dtype=int)
    n_frag, fragment_smiles, fragment_species, _ = tag_mixture(frag, resolved, 0)
    set_fragments(packed, frag)

    extra: dict = {"n_molecules": int(n_frag), "box": [round(x, 4) for x in box]}
    if fragment_species:
        extra["fragment_species"] = fragment_species
    if fragment_smiles:
        extra["fragment_smiles"] = fragment_smiles

    return Structure.from_ase(
        packed,
        provenance=Provenance(origin="generated", source="builder:liquid", extra=extra),
    )
