"""Tests for the shared mixture (Species) + solid-solution (alloy) machinery."""

from __future__ import annotations

import shutil

import numpy as np
import pytest

from traincraft.config.models import (
    AlloyComponent,
    CrystalBuilder,
    FilledNanotubeBuilder,
    GeometryConfig,
    LiquidBuilder,
    SlabBuilder,
    Species,
    SurfacePackingBuilder,
)
from traincraft.core import FRAMEWORK
from traincraft.geometry import build_geometry
from traincraft.geometry.builders.alloy import apply_solid_solution
from traincraft.geometry.builders.mixture import apportion, species_counts

_HAS_PACKMOL = shutil.which("packmol") is not None


# --- apportionment (pure, no deps) -----------------------------------------

def test_apportion_sums_exactly():
    assert sum(apportion([1, 1, 1], 10)) == 10
    assert sum(apportion([0.5, 0.3, 0.2], 7)) == 7
    assert apportion([1, 1], 4) == [2, 2]


def test_apportion_largest_remainder():
    # 0.5/0.3/0.2 of 7 -> 3.5/2.1/1.4 -> floors 3/2/1 (=6), leftover to largest frac (0.5)
    assert apportion([0.5, 0.3, 0.2], 7) == [4, 2, 1]


def test_species_counts_count_mode():
    sp = [Species(molecule_name="H2O", count=8), Species(molecule_name="CH4", count=2)]
    assert species_counts(sp, None) == [8, 2]


def test_species_counts_bare_defaults_to_one():
    sp = [Species(molecule_name="H2O"), Species(molecule_name="CH4")]
    assert species_counts(sp, None) == [1, 1]


def test_species_counts_ratio_mode():
    sp = [Species(molecule_name="H2O", ratio=3), Species(molecule_name="CH4", ratio=1)]
    assert species_counts(sp, 20) == [15, 5]


def test_species_counts_ratio_needs_total():
    sp = [Species(molecule_name="H2O", ratio=1)]
    with pytest.raises(ValueError, match="needs a total"):
        species_counts(sp, None)


# --- config validation ------------------------------------------------------

def test_species_count_xor_ratio():
    with pytest.raises(ValueError, match="not both"):
        Species(molecule_name="H2O", count=2, ratio=1.0)


def test_no_mixing_count_and_ratio_in_builder():
    with pytest.raises(ValueError, match="don't mix"):
        FilledNanotubeBuilder(
            species=[
                Species(molecule_name="H2O", count=2),
                Species(molecule_name="CH4", ratio=1.0),
            ]
        )


def test_single_shortcut_xor_species_list():
    with pytest.raises(ValueError, match="not both"):
        SurfacePackingBuilder(
            element="Cu", molecule_name="CO", species=[Species(molecule_name="H2O")]
        )


def test_ratio_liquid_needs_n_molecules():
    with pytest.raises(ValueError, match="n_molecules"):
        LiquidBuilder(species=[Species(molecule_name="H2O", ratio=1.0)], box=(10, 10, 10))


# --- solid solution (alloy), no Packmol ------------------------------------

def test_apply_solid_solution_counts_and_seed():
    from ase.build import bulk

    atoms = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((3, 3, 3))  # 108 sites
    comps = [AlloyComponent(element="Au", ratio=0.25)]
    out, comp = apply_solid_solution(atoms, comps, seed=1)
    n_au = int((np.array(out.get_chemical_symbols()) == "Au").sum())
    assert n_au == 27  # round(0.25 * 108)
    assert comp == {"Au": 27}
    # deterministic with the same seed, different with another
    out2, _ = apply_solid_solution(atoms, comps, seed=1)
    out3, _ = apply_solid_solution(atoms, comps, seed=2)
    assert out.get_chemical_symbols() == out2.get_chemical_symbols()
    assert out.get_chemical_symbols() != out3.get_chemical_symbols()


def test_alloy_ratios_must_not_exceed_one():
    from ase.build import bulk

    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    comps = [AlloyComponent(element="Au", ratio=0.7), AlloyComponent(element="Ag", ratio=0.5)]
    with pytest.raises(ValueError, match="sum to"):
        apply_solid_solution(atoms, comps, seed=1)


def test_crystal_builder_alloy_provenance():
    s = build_geometry(
        GeometryConfig(
            builder=CrystalBuilder(
                name="Cu", cubic=True, supercell=(2, 2, 2),
                composition=[AlloyComponent(element="Au", ratio=0.5)], seed=3,
            )
        )
    )
    assert "composition" in s.provenance.extra
    syms = np.array(s.atoms.get_chemical_symbols())
    assert (syms == "Au").sum() == s.provenance.extra["composition"]["Au"]


def test_slab_builder_alloy():
    s = build_geometry(
        GeometryConfig(
            builder=SlabBuilder(
                element="Cu", facet="fcc111", size=(3, 3, 4),
                composition=[AlloyComponent(element="Ni", ratio=0.33)], seed=7,
            )
        )
    )
    assert "Ni" in np.array(s.atoms.get_chemical_symbols())
    assert s.provenance.extra["composition"]["Ni"] > 0


# --- end-to-end packed mixtures (need Packmol) -----------------------------

@pytest.mark.skipif(not _HAS_PACKMOL, reason="packmol binary not installed")
def test_liquid_ratio_mixture():
    s = build_geometry(
        GeometryConfig(
            builder=LiquidBuilder(
                species=[
                    Species(molecule_name="H2O", ratio=3),
                    Species(molecule_name="CH4", ratio=1),
                ],
                n_molecules=8,
                box=(16.0, 16.0, 16.0),
                seed=1,
            )
        )
    )
    # 6 waters (3 atoms) + 2 methane (5 atoms) = 18 + 10 = 28 atoms, 8 fragments
    assert s.provenance.extra["n_molecules"] == 8
    frag = s.atoms.arrays["tc_fragment"]
    assert sorted(set(int(f) for f in frag)) == list(range(8))
    species_map = s.provenance.extra["fragment_species"]
    assert sum(v == "H2O" for v in species_map.values()) == 6
    assert sum(v == "CH4" for v in species_map.values()) == 2


@pytest.mark.skipif(not _HAS_PACKMOL, reason="packmol binary not installed")
def test_filled_nanotube_mixture():
    s = build_geometry(
        GeometryConfig(
            builder=FilledNanotubeBuilder(
                n=12, m=12, length=6,
                species=[
                    Species(molecule_name="H2O", count=4),
                    Species(molecule_name="CH4", count=2),
                ],
                seed=2,
            )
        )
    )
    atoms = s.atoms
    frag = atoms.arrays["tc_fragment"]
    n_guest = 4 * 3 + 2 * 5  # 4 H2O + 2 CH4
    n_tube = len(atoms) - n_guest
    assert (frag[:n_tube] == FRAMEWORK).all()  # tube is framework, listed first
    assert (frag[n_tube:] >= 0).all()  # guests are mobile fragments
    assert s.provenance.extra["n_molecules"] == 6
    species_map = s.provenance.extra["fragment_species"]
    assert sum(v == "H2O" for v in species_map.values()) == 4
    assert sum(v == "CH4" for v in species_map.values()) == 2


@pytest.mark.skipif(not _HAS_PACKMOL, reason="packmol binary not installed")
def test_surface_packing_mixture_on_alloy():
    s = build_geometry(
        GeometryConfig(
            builder=SurfacePackingBuilder(
                element="Cu", facet="fcc111", size=(4, 4, 4),
                composition=[AlloyComponent(element="Au", ratio=0.5)],
                species=[
                    Species(molecule_name="CO", ratio=1),
                    Species(molecule_name="H2O", ratio=1),
                ],
                n_molecules=4,
                seed=5,
            )
        )
    )
    frag = s.atoms.arrays["tc_fragment"]
    assert int((frag == FRAMEWORK).sum()) == 4 * 4 * 4
    assert s.provenance.extra["n_molecules"] == 4
    assert "Au" in np.array(s.atoms.get_chemical_symbols())
    assert "composition" in s.provenance.extra
