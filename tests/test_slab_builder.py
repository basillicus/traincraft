from __future__ import annotations

import numpy as np
import pytest

from traincraft.config.models import GeometryConfig, SlabBuilder
from traincraft.core.fragments import FRAMEWORK, get_fragments
from traincraft.geometry import build_geometry


def test_slab_named_facet():
    s = build_geometry(
        GeometryConfig(builder=SlabBuilder(element="Cu", facet="fcc111", size=(2, 2, 4)))
    )
    assert len(s.atoms) == 2 * 2 * 4
    pbc = s.atoms.get_pbc()
    assert pbc[0] and pbc[1] and not pbc[2]
    assert s.provenance.source == "builder:slab:Cu-fcc111"


def test_slab_is_all_framework():
    s = build_geometry(
        GeometryConfig(builder=SlabBuilder(element="Cu", facet="fcc111", size=(2, 2, 3)))
    )
    frag = get_fragments(s.atoms)
    assert frag is not None
    assert np.all(frag == FRAMEWORK)
    assert s.n_fragments == 0


def test_slab_miller_indices():
    s = build_geometry(
        GeometryConfig(builder=SlabBuilder(
            element="Cu", miller=(1, 1, 0), crystalstructure="fcc", a=3.6, layers=6
        ))
    )
    assert len(s.atoms) > 0
    pbc = s.atoms.get_pbc()
    assert pbc[0] and pbc[1] and not pbc[2]
    assert "Cu-110" in s.provenance.source


def test_slab_miller_clears_facet_default():
    cfg = SlabBuilder(element="Cu", miller=(2, 1, 1), crystalstructure="fcc", a=3.6)
    assert cfg.facet is None


def test_slab_requires_facet_or_miller():
    with pytest.raises(ValueError, match="facet"):
        SlabBuilder(element="Cu", facet=None)
