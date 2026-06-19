from __future__ import annotations

import numpy as np
import pytest

from traincraft.config.models import MaterialsProjectSource, OptimadeSource, PubchemSource
from traincraft.geometry.sources.providers import _atoms_from_optimade

# --- OPTIMADE JSON -> Atoms parser (no network) ----------------------------

def test_atoms_from_optimade_periodic():
    attrs = {
        "species_at_sites": ["Na", "Cl"],
        "species": [
            {"name": "Na", "chemical_symbols": ["Na"]},
            {"name": "Cl", "chemical_symbols": ["Cl"]},
        ],
        "cartesian_site_positions": [[0.0, 0.0, 0.0], [2.8, 2.8, 2.8]],
        "lattice_vectors": [[5.6, 0, 0], [0, 5.6, 0], [0, 0, 5.6]],
        "dimension_types": [1, 1, 1],
    }
    atoms = _atoms_from_optimade(attrs)
    assert atoms.get_chemical_symbols() == ["Na", "Cl"]
    assert atoms.get_pbc().all()
    assert np.allclose(np.diag(np.asarray(atoms.get_cell())), [5.6, 5.6, 5.6])
    assert np.allclose(atoms.get_positions()[1], [2.8, 2.8, 2.8])


def test_atoms_from_optimade_molecular_no_cell():
    attrs = {
        "species_at_sites": ["H", "H"],
        "species": [{"name": "H", "chemical_symbols": ["H"]}],
        "cartesian_site_positions": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        "dimension_types": [0, 0, 0],
    }
    atoms = _atoms_from_optimade(attrs)
    assert not atoms.get_pbc().any()
    assert len(atoms) == 2


# --- config validation (no network) ----------------------------------------

def test_pubchem_needs_exactly_one_identifier():
    with pytest.raises(ValueError, match="exactly one"):
        PubchemSource()
    with pytest.raises(ValueError, match="exactly one"):
        PubchemSource(name="water", cid=962)
    PubchemSource(name="water")  # ok


def test_optimade_and_mp_minimal_configs():
    assert OptimadeSource(base_url="https://example.org/optimade/v1").type == "optimade"
    assert MaterialsProjectSource(material_id="mp-149").type == "materials_project"
