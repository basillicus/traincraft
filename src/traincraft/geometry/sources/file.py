"""Source: read a structure from any ASE-readable file."""

from __future__ import annotations

from ase.io import read

from ...core import Provenance, Structure, register


@register("source", "file")
def source_file(cfg) -> Structure:
    atoms = read(cfg.path)
    return Structure.from_ase(
        atoms, provenance=Provenance(origin="generated", source=f"source:file:{cfg.path}")
    )
