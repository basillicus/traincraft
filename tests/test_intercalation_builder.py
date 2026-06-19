from __future__ import annotations

import pytest

from traincraft.config.models import GeometryConfig, IntercalationBuilder, LayeredBuilder
from traincraft.geometry import build_geometry


def test_intercalation_fills_every_gallery():
    s = build_geometry(
        GeometryConfig(
            builder=IntercalationBuilder(
                host=LayeredBuilder(material="graphene", n_layers=3),
                guest="Li",
                n_per_gallery=2,
            )
        )
    )
    symbols = s.atoms.get_chemical_symbols()
    # 3 layers -> 2 galleries -> 2 * n_per_gallery guests.
    assert symbols.count("Li") == 4
    frag = s.atoms.arrays["tc_fragment"]
    assert int((frag == -1).sum()) == 2 * 3  # graphene host: 2 atoms/cell x 3 layers
    assert int((frag == 0).sum()) == 4  # all guests = one fragment
    assert s.provenance.extra["galleries_filled"] == 2


def test_intercalation_staging_skips_galleries():
    s = build_geometry(
        GeometryConfig(
            builder=IntercalationBuilder(
                host=LayeredBuilder(material="graphene", n_layers=3),
                guest="Li",
                n_per_gallery=1,
                stage=2,
            )
        )
    )
    # 2 galleries, stage 2 -> only gallery index 0 is filled.
    assert s.atoms.get_chemical_symbols().count("Li") == 1


def test_intercalation_rejects_mx2_host():
    with pytest.raises(ValueError, match="planar"):
        build_geometry(
            GeometryConfig(
                builder=IntercalationBuilder(
                    host=LayeredBuilder(material="mx2", formula="MoS2", n_layers=2)
                )
            )
        )


def test_intercalation_rejects_twisted_host():
    with pytest.raises(ValueError, match="untwisted"):
        build_geometry(
            GeometryConfig(
                builder=IntercalationBuilder(
                    host=LayeredBuilder(material="graphene", n_layers=2, twist=10.0)
                )
            )
        )
