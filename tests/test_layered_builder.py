from __future__ import annotations

import numpy as np

from traincraft.config.models import GeometryConfig, LayeredBuilder
from traincraft.geometry import build_geometry


def test_layered_graphene_periodic_stack():
    s = build_geometry(
        GeometryConfig(builder=LayeredBuilder(material="graphene", n_layers=3))
    )
    assert len(s.atoms) == 2 * 3  # 2 atoms/primitive cell x 3 layers
    pbc = s.atoms.get_pbc()
    assert pbc[0] and pbc[1] and not pbc[2]


def test_layered_interlayer_spacing():
    spacing = 4.0
    s = build_geometry(
        GeometryConfig(builder=LayeredBuilder(
            material="graphene", n_layers=2, interlayer_spacing=spacing, stacking="AA"
        ))
    )
    z = np.sort(np.unique(np.round(s.atoms.get_positions()[:, 2], 3)))
    # AA: two flat layers; their z separation equals the spacing.
    assert np.isclose(z[-1] - z[0], spacing, atol=1e-6)


def test_layered_ab_offsets_layers():
    aa = build_geometry(GeometryConfig(builder=LayeredBuilder(
        material="graphene", n_layers=2, stacking="AA")))
    ab = build_geometry(GeometryConfig(builder=LayeredBuilder(
        material="graphene", n_layers=2, stacking="AB")))
    # AB shifts the second layer in-plane, so xy positions differ from AA.
    assert not np.allclose(
        np.sort(aa.atoms.get_positions()[:, 0]),
        np.sort(ab.atoms.get_positions()[:, 0]),
    )


def test_layered_twist_is_nonperiodic_flake():
    s = build_geometry(
        GeometryConfig(builder=LayeredBuilder(
            material="graphene", n_layers=2, twist=21.8, size=(3, 3)
        ))
    )
    assert not s.atoms.get_pbc().any()
    assert s.provenance.extra["twist"] == 21.8


def test_layered_hbn_and_mx2():
    hbn = build_geometry(GeometryConfig(builder=LayeredBuilder(material="hbn", n_layers=2)))
    assert sorted(set(hbn.atoms.get_chemical_symbols())) == ["B", "N"]

    mos2 = build_geometry(GeometryConfig(builder=LayeredBuilder(
        material="mx2", formula="MoS2", n_layers=2)))
    assert sorted(set(mos2.atoms.get_chemical_symbols())) == ["Mo", "S"]
