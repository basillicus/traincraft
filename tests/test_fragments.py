from __future__ import annotations

import numpy as np
import pytest
from ase.build import molecule

from traincraft.core.fragments import (
    FRAMEWORK,
    fragment_ids,
    fragment_mask,
    get_fragments,
    infer_fragments,
    set_fragments,
)
from traincraft.core.provenance import Provenance
from traincraft.core.structure import Structure
from traincraft.datasets.io import read_frames, write_frames


def _h2o():
    return molecule("H2O")


def _co2():
    return molecule("CO2")


# --- set / get round-trip ---------------------------------------------------

def test_set_get_roundtrip():
    atoms = _h2o()
    arr = np.array([0, 0, 0])
    set_fragments(atoms, arr)
    result = get_fragments(atoms)
    assert result is not None
    np.testing.assert_array_equal(result, arr)


def test_get_returns_none_when_unset():
    atoms = _h2o()
    assert get_fragments(atoms) is None


def test_set_wrong_length_rejected():
    atoms = _h2o()  # 3 atoms
    with pytest.raises(ValueError, match="shape"):
        set_fragments(atoms, [0, 0])  # wrong length


def test_fragment_ids_excludes_framework():
    atoms = _h2o()
    set_fragments(atoms, [-1, 0, 0])
    assert fragment_ids(atoms) == [0]


def test_fragment_mask():
    atoms = _h2o()
    set_fragments(atoms, [0, 0, 1])
    mask = fragment_mask(atoms, 0)
    np.testing.assert_array_equal(mask, [True, True, False])


def test_fragment_mask_no_array_raises():
    atoms = _h2o()
    with pytest.raises(ValueError, match="no fragment array"):
        fragment_mask(atoms, 0)


# --- ase.Atoms.copy() preserves the array -----------------------------------

def test_copy_preserves_fragments():
    atoms = _h2o()
    set_fragments(atoms, [0, 0, 0])
    copy = atoms.copy()
    result = get_fragments(copy)
    assert result is not None
    np.testing.assert_array_equal(result, [0, 0, 0])


# --- infer_fragments ---------------------------------------------------------

def test_infer_two_isolated_molecules():
    """Two molecules placed far apart should yield 2 distinct component ids."""
    from ase import Atoms

    h2o = _h2o()
    co2 = _co2()

    # Place CO2 far from H2O (no cell, so use large position offset).
    co2_pos = co2.get_positions() + np.array([20.0, 0.0, 0.0])
    combined = h2o + Atoms(
        symbols=co2.get_chemical_symbols(),
        positions=co2_pos,
    )

    result = infer_fragments(combined)
    ids = sorted(int(i) for i in np.unique(result) if i != FRAMEWORK)
    assert len(ids) == 2, f"expected 2 fragments, got {ids}"


def test_infer_framework_mask():
    """Atoms in framework_mask are forced to FRAMEWORK (-1)."""
    atoms = _h2o()  # H2O: O bonded to each H, but H-H not bonded
    mask = np.array([True, False, False])  # mark O as framework
    result = infer_fragments(atoms, framework_mask=mask)
    assert result[0] == FRAMEWORK
    # The two H atoms are only bonded via O; removing O leaves 2 isolated components.
    mobile_ids = [int(i) for i in result[1:] if i != FRAMEWORK]
    assert len(mobile_ids) == 2
    assert len(set(mobile_ids)) == 2  # distinct ids


def test_infer_framework_mask_wrong_length():
    atoms = _h2o()
    with pytest.raises(ValueError, match="framework_mask"):
        infer_fragments(atoms, framework_mask=np.array([True, False]))


# --- extxyz round-trip -------------------------------------------------------

def test_extxyz_roundtrip_preserves_fragments(tmp_path):
    atoms = _h2o()
    arr = np.array([0, 0, 0], dtype=int)
    set_fragments(atoms, arr)

    s = Structure(
        atoms=atoms,
        provenance=Provenance(origin="generated"),
    )

    path = tmp_path / "out.extxyz"
    write_frames(path, [s])
    loaded = read_frames(path)

    assert len(loaded) == 1
    recovered = get_fragments(loaded[0].atoms)
    assert recovered is not None, "tc_fragment array lost in extxyz round-trip"
    np.testing.assert_array_equal(recovered, arr)


# --- Structure helpers -------------------------------------------------------

def test_structure_n_fragments():
    atoms = _h2o()
    set_fragments(atoms, [0, 0, 0])
    s = Structure(atoms=atoms, provenance=Provenance())
    assert s.n_fragments == 1


def test_structure_set_fragments_helper():
    atoms = _h2o()
    s = Structure(atoms=atoms, provenance=Provenance())
    s.set_fragments([0, 0, 0])
    assert s.fragments is not None
    np.testing.assert_array_equal(s.fragments, [0, 0, 0])
