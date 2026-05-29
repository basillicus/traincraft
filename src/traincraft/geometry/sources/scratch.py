"""Source: build a structure from scratch via ASE primitives."""

from __future__ import annotations

from ...core import Provenance, Structure, register


@register("source", "scratch")
def source_scratch(cfg) -> Structure:
    if cfg.molecule:
        from ase.build import molecule

        atoms = molecule(cfg.molecule)
        tag = f"molecule:{cfg.molecule}"
    elif cfg.bulk:
        from ase.build import bulk

        atoms = bulk(cfg.bulk)
        tag = f"bulk:{cfg.bulk}"
    else:
        raise ValueError("scratch source needs 'molecule' or 'bulk'")
    return Structure.from_ase(
        atoms, provenance=Provenance(origin="generated", source=f"source:scratch:{tag}")
    )
