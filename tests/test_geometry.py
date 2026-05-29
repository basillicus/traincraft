from __future__ import annotations

from traincraft.config.models import (
    GeometryConfig,
    MoleculeBuilder,
    NanotubeBuilder,
    VacuumTransform,
)
from traincraft.geometry import build_geometry


def test_nanotube_builder():
    s = build_geometry(GeometryConfig(builder=NanotubeBuilder(n=5, m=0, length=1)))
    assert len(s.atoms) > 0
    assert s.provenance.source.startswith("builder:nanotube")


def test_molecule_builder_with_transform():
    g = GeometryConfig(
        builder=MoleculeBuilder(name="H2O"),
        transforms=[VacuumTransform(amount=8.0)],
    )
    s = build_geometry(g)
    assert len(s.atoms) == 3
    assert "vacuum" in s.provenance.transforms
