from __future__ import annotations

import numpy as np
import pytest
from ase.build import bulk, molecule

from traincraft.core import Structure
from traincraft.core.converter import (
    ase_to_pymatgen,
    ase_to_rdkit,
    pymatgen_to_ase,
    rdkit_to_ase,
)

# ---------------------------------------------------------------------------
# pymatgen

def test_ase_to_pymatgen_structure_for_periodic():
    pytest.importorskip("pymatgen")
    cu = bulk("Cu", "fcc", a=3.6, cubic=True)
    pmg = ase_to_pymatgen(cu)
    assert type(pmg).__name__ == "Structure"
    assert len(pmg) == len(cu)


def test_ase_to_pymatgen_molecule_for_nonperiodic():
    pytest.importorskip("pymatgen")
    pmg = ase_to_pymatgen(molecule("H2O"))
    assert type(pmg).__name__ == "Molecule"


def test_pymatgen_roundtrip_preserves_cell():
    pytest.importorskip("pymatgen")
    cu = bulk("Cu", "fcc", a=3.6, cubic=True)
    back = pymatgen_to_ase(ase_to_pymatgen(cu))
    assert len(back) == len(cu)
    assert np.allclose(back.get_cell(), cu.get_cell())


# ---------------------------------------------------------------------------
# rdkit

def test_ase_to_rdkit_perceives_bonds():
    pytest.importorskip("rdkit")
    mol = ase_to_rdkit(molecule("H2O"))
    assert mol.GetNumAtoms() == 3
    assert mol.GetNumBonds() == 2


def test_rdkit_to_ase_roundtrip():
    pytest.importorskip("rdkit")
    h2o = molecule("H2O")
    back = rdkit_to_ase(ase_to_rdkit(h2o))
    assert sorted(back.get_chemical_symbols()) == sorted(h2o.get_chemical_symbols())


def test_ase_to_rdkit_rejects_periodic():
    pytest.importorskip("rdkit")
    with pytest.raises(ValueError, match="non-periodic"):
        ase_to_rdkit(bulk("Cu", "fcc", a=3.6, cubic=True))


# ---------------------------------------------------------------------------
# Structure wiring

def test_structure_to_from_pymatgen():
    pytest.importorskip("pymatgen")
    cu = bulk("Cu", "fcc", a=3.6, cubic=True)
    s = Structure.from_ase(cu)
    pmg = s.to_pymatgen()
    assert len(Structure.from_pymatgen(pmg).atoms) == len(cu)


def test_structure_to_from_rdkit():
    pytest.importorskip("rdkit")
    s = Structure.from_ase(molecule("CH4"))
    mol = s.to_rdkit()
    assert mol.GetNumAtoms() == 5
    assert len(Structure.from_rdkit(mol).atoms) == 5
