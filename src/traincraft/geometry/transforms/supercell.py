"""Transform: build a supercell by repeating the cell."""

from __future__ import annotations

from ...core import Provenance, Structure, register


@register("transform", "supercell")
def transform_supercell(structure: Structure, cfg) -> Structure:
    atoms = structure.atoms.repeat(tuple(cfg.repeat))
    prov = Provenance.from_dict(structure.provenance.to_dict())
    prov.transforms.append(f"supercell:{tuple(cfg.repeat)}")
    return Structure(atoms=atoms, properties={}, provenance=prov)
