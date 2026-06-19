from __future__ import annotations

import pytest

from traincraft.config.models import CrystalBuilder, DefectSpec, GeometryConfig
from traincraft.geometry import build_geometry


def _cu(**kw):
    return CrystalBuilder(
        name="Cu", crystalstructure="fcc", a=3.6, cubic=True, **kw
    )


def test_crystal_bulk_periodic():
    s = build_geometry(GeometryConfig(builder=_cu()))
    assert len(s.atoms) == 4  # conventional fcc cell
    assert s.atoms.get_pbc().all()
    assert s.provenance.source == "builder:crystal:Cu"


def test_crystal_supercell():
    s = build_geometry(GeometryConfig(builder=_cu(supercell=(2, 2, 2))))
    assert len(s.atoms) == 4 * 8


def test_crystal_vacancy_removes_one_atom():
    base = build_geometry(GeometryConfig(builder=_cu(supercell=(2, 2, 2))))
    s = build_geometry(
        GeometryConfig(builder=_cu(supercell=(2, 2, 2),
                                   defects=[DefectSpec(kind="vacancy", index=0)]))
    )
    assert len(s.atoms) == len(base.atoms) - 1
    assert s.provenance.extra["defects"][0]["kind"] == "vacancy"


def test_crystal_substitution():
    s = build_geometry(
        GeometryConfig(builder=_cu(supercell=(2, 2, 2),
                                   defects=[DefectSpec(kind="substitution", index=0,
                                                       element="Ni")]))
    )
    assert "Ni" in s.atoms.get_chemical_symbols()
    assert s.atoms.get_chemical_symbols().count("Ni") == 1


def test_crystal_interstitial_adds_atom():
    base = build_geometry(GeometryConfig(builder=_cu(supercell=(2, 2, 2))))
    s = build_geometry(
        GeometryConfig(builder=_cu(supercell=(2, 2, 2),
                                   defects=[DefectSpec(kind="interstitial", element="H",
                                                       position=(0.5, 0.5, 0.5))]))
    )
    assert len(s.atoms) == len(base.atoms) + 1
    assert "H" in s.atoms.get_chemical_symbols()


def test_crystal_combined_defects_stable_indexing():
    # vacancy at 0 + substitution at 1 + interstitial: all index refer to supercell
    s = build_geometry(
        GeometryConfig(builder=_cu(supercell=(2, 2, 2), defects=[
            DefectSpec(kind="substitution", index=1, element="Ni"),
            DefectSpec(kind="vacancy", index=0),
            DefectSpec(kind="interstitial", element="H", position=(0.1, 0.1, 0.1)),
        ]))
    )
    syms = s.atoms.get_chemical_symbols()
    assert syms.count("Ni") == 1
    assert "H" in syms
    assert len(s.atoms) == 4 * 8 - 1 + 1


def test_crystal_vacancy_index_out_of_range():
    with pytest.raises(IndexError):
        build_geometry(
            GeometryConfig(builder=_cu(defects=[DefectSpec(kind="vacancy", index=999)]))
        )


def test_defect_validation_requires_element():
    with pytest.raises(ValueError, match="substitution"):
        DefectSpec(kind="substitution", index=0)
    with pytest.raises(ValueError, match="interstitial"):
        DefectSpec(kind="interstitial", element="H")
