from __future__ import annotations

import numpy as np
import pytest

from traincraft.config.models import (
    CrystalBuilder,
    GeometryConfig,
    RotateTransform,
    SetPbcTransform,
    SlabBuilder,
    StrainTransform,
)
from traincraft.geometry import build_geometry


def _cu_bulk():
    return CrystalBuilder(name="Cu", crystalstructure="fcc", a=3.6, cubic=True)


# ---------------------------------------------------------------------------
# strain

def test_strain_hydrostatic_scales_volume():
    base = build_geometry(GeometryConfig(builder=_cu_bulk()))
    v0 = base.atoms.get_volume()
    s = build_geometry(
        GeometryConfig(builder=_cu_bulk(), transforms=[StrainTransform(hydrostatic=0.02)])
    )
    assert np.isclose(s.atoms.get_volume(), v0 * 1.02**3, rtol=1e-6)
    assert any(t.startswith("strain") for t in s.provenance.transforms)


def test_strain_voigt_uniaxial():
    base = build_geometry(GeometryConfig(builder=_cu_bulk()))
    c0 = np.asarray(base.atoms.get_cell())
    s = build_geometry(
        GeometryConfig(builder=_cu_bulk(),
                       transforms=[StrainTransform(voigt=(0.05, 0.0, 0.0, 0, 0, 0))])
    )
    c1 = np.asarray(s.atoms.get_cell())
    # x lattice vector stretched by 5%, others unchanged.
    assert np.isclose(np.linalg.norm(c1[0]), np.linalg.norm(c0[0]) * 1.05, rtol=1e-6)


def test_strain_requires_one_form():
    with pytest.raises(ValueError, match="exactly one"):
        StrainTransform()
    with pytest.raises(ValueError, match="exactly one"):
        StrainTransform(hydrostatic=0.01, voigt=(0, 0, 0, 0, 0, 0))


def test_strain_on_zero_cell_raises():
    # scratch source builds a molecule with no cell (unlike the molecule builder,
    # which centers it in a vacuum box).
    from traincraft.config.models import ScratchSource
    with pytest.raises(ValueError, match="periodic cell"):
        build_geometry(
            GeometryConfig(source=ScratchSource(molecule="H2O"),
                           transforms=[StrainTransform(hydrostatic=0.01)])
        )


# ---------------------------------------------------------------------------
# rotate

def test_rotate_preserves_atom_count_and_geometry():
    from traincraft.config.models import MoleculeBuilder
    base = build_geometry(GeometryConfig(builder=MoleculeBuilder(name="H2O")))
    s = build_geometry(
        GeometryConfig(builder=MoleculeBuilder(name="H2O"),
                       transforms=[RotateTransform(angle=90, axis="z")])
    )
    assert len(s.atoms) == len(base.atoms)
    # Rotation is rigid: pairwise distances are invariant.
    assert np.allclose(
        np.sort(base.atoms.get_all_distances().ravel()),
        np.sort(s.atoms.get_all_distances().ravel()),
    )


# ---------------------------------------------------------------------------
# set_pbc

def test_set_pbc_bool():
    s = build_geometry(
        GeometryConfig(builder=SlabBuilder(element="Cu", facet="fcc111", size=(2, 2, 3)),
                       transforms=[SetPbcTransform(pbc=True)])
    )
    assert s.atoms.get_pbc().all()


def test_set_pbc_tuple():
    s = build_geometry(
        GeometryConfig(builder=SlabBuilder(element="Cu", facet="fcc111", size=(2, 2, 3)),
                       transforms=[SetPbcTransform(pbc=(False, False, False))])
    )
    assert not s.atoms.get_pbc().any()
