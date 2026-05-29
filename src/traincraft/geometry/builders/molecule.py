"""Builder: a single molecule from an ASE g2 name (SMILES path: geometry chunk)."""

from __future__ import annotations

from ...core import Provenance, Structure, register


@register("builder", "molecule")
def build_molecule(cfg) -> Structure:
    if cfg.name:
        from ase.build import molecule

        atoms = molecule(cfg.name)
        tag = cfg.name
    elif cfg.smiles:  # pragma: no cover - Phase 2
        raise NotImplementedError("SMILES builder (RDKit) lands in the geometry chunk")
    else:
        raise ValueError("molecule builder needs 'name' or 'smiles'")
    atoms.center(vacuum=cfg.vacuum / 2)
    return Structure.from_ase(
        atoms,
        provenance=Provenance(origin="generated", source=f"builder:molecule:{tag}"),
    )
