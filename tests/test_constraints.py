from __future__ import annotations

import numpy as np
import pytest

from traincraft.config.models import (
    ConstraintsTransform,
    GeometryConfig,
    SlabBuilder,
    SurfaceAdsorbateBuilder,
)
from traincraft.geometry import build_geometry


def _slab():
    return SlabBuilder(element="Cu", facet="fcc111", size=(2, 2, 3))


def _fixed_indices(structure):
    (con,) = structure.atoms.constraints
    return set(int(i) for i in con.index)


def test_constraints_by_elements_fixes_all_cu():
    s = build_geometry(
        GeometryConfig(builder=_slab(), transforms=[ConstraintsTransform(elements=["Cu"])])
    )
    assert _fixed_indices(s) == set(range(len(s.atoms)))
    assert any(t.startswith("constraints") for t in s.provenance.transforms)


def test_constraints_by_indices():
    s = build_geometry(
        GeometryConfig(builder=_slab(), transforms=[ConstraintsTransform(indices=[0, 3, 5])])
    )
    assert _fixed_indices(s) == {0, 3, 5}


def test_constraints_below_z_fixes_bottom_layer():
    base = build_geometry(GeometryConfig(builder=_slab()))
    z = base.atoms.get_positions()[:, 2]
    threshold = z.min() + 0.5  # only the lowest layer sits below this
    expected = set(int(i) for i in np.nonzero(z < threshold)[0])

    s = build_geometry(
        GeometryConfig(builder=_slab(), transforms=[ConstraintsTransform(below_z=threshold)])
    )
    assert _fixed_indices(s) == expected
    assert 0 < len(expected) < len(base.atoms)


def test_constraints_by_fragment_fixes_substrate():
    # surface_adsorbate tags substrate as fragment -1, adsorbate as 0.
    builder = SurfaceAdsorbateBuilder(element="Cu", molecule_name="CO", size=(2, 2, 3))
    base = build_geometry(GeometryConfig(builder=builder))
    n_substrate = int((base.atoms.arrays["tc_fragment"] == -1).sum())

    s = build_geometry(
        GeometryConfig(builder=builder, transforms=[ConstraintsTransform(fragments=[-1])])
    )
    assert len(_fixed_indices(s)) == n_substrate


def test_constraints_index_out_of_range():
    with pytest.raises(IndexError):
        build_geometry(
            GeometryConfig(builder=_slab(), transforms=[ConstraintsTransform(indices=[10_000])])
        )


def test_constraints_selecting_nothing_raises():
    with pytest.raises(ValueError, match="selected no atoms"):
        build_geometry(
            GeometryConfig(builder=_slab(), transforms=[ConstraintsTransform(below_z=-1e6)])
        )
