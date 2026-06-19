from __future__ import annotations

import shutil

import numpy as np
import pytest

from traincraft.config.models import GeometryConfig, LiquidBuilder, LiquidSpecies
from traincraft.geometry import build_geometry
from traincraft.geometry.builders.liquid import _cubic_edge_from_density

_HAS_PACKMOL = shutil.which("packmol") is not None


# --- validation + math (no Packmol needed) ---------------------------------

def test_liquid_needs_box_xor_density():
    with pytest.raises(ValueError, match="exactly one"):
        LiquidBuilder(species=[LiquidSpecies(molecule_name="H2O")])
    with pytest.raises(ValueError, match="exactly one"):
        LiquidBuilder(species=[LiquidSpecies(molecule_name="H2O")], box=(10, 10, 10), density=1.0)


def test_liquid_species_needs_one_source():
    with pytest.raises(ValueError, match="exactly one"):
        LiquidSpecies(molecule_name="H2O", smiles="O")


def test_cubic_edge_from_density_matches_water():
    # 33 water molecules (~18 amu each) at 1 g/cm^3 -> ~10 Å cube.
    edge = _cubic_edge_from_density(33 * 18.015, 1.0)
    assert np.isclose(edge, 10.0, atol=0.3)


# --- end-to-end packing (needs the Packmol binary) -------------------------

@pytest.mark.skipif(not _HAS_PACKMOL, reason="packmol binary not installed")
def test_liquid_packs_water_box():
    s = build_geometry(
        GeometryConfig(
            builder=LiquidBuilder(
                species=[LiquidSpecies(molecule_name="H2O", count=8)],
                box=(12.0, 12.0, 12.0),
                seed=1,
            )
        )
    )
    assert len(s.atoms) == 8 * 3  # 8 waters
    assert s.atoms.get_pbc().all()
    assert np.allclose(np.diag(np.asarray(s.atoms.get_cell())), [12.0, 12.0, 12.0])
    frag = s.atoms.arrays["tc_fragment"]
    assert sorted(set(int(f) for f in frag)) == list(range(8))
    assert s.provenance.extra["n_molecules"] == 8
