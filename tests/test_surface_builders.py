from __future__ import annotations

import shutil

import numpy as np
import pytest

from traincraft.config.models import SurfaceAdsorbateBuilder, SurfacePackingBuilder
from traincraft.core.fragments import FRAMEWORK, get_fragments
from traincraft.geometry.builders.surface import (
    build_surface_adsorbate,
    build_surface_packing,
)

# ---------------------------------------------------------------------------
# surface_adsorbate

def test_surface_adsorbate_cu_co():
    cfg = SurfaceAdsorbateBuilder(
        element="Cu", facet="fcc111", size=(2, 2, 3), molecule_name="CO", site="ontop"
    )
    s = build_surface_adsorbate(cfg)

    frag = get_fragments(s.atoms)
    assert frag is not None

    n_framework = int(np.sum(frag == FRAMEWORK))
    n_mobile = int(np.sum(frag >= 0))

    # slab = 2*2*3 = 12 Cu atoms; CO = 2 atoms
    assert n_framework == 12
    assert n_mobile == 2
    assert s.n_fragments == 1


def test_surface_adsorbate_pbc():
    cfg = SurfaceAdsorbateBuilder(
        element="Cu", facet="fcc111", size=(2, 2, 3), molecule_name="CO"
    )
    s = build_surface_adsorbate(cfg)
    pbc = s.atoms.get_pbc()
    assert pbc[0] and pbc[1]  # periodic in xy
    assert not pbc[2]          # non-periodic in z


def test_surface_adsorbate_provenance():
    cfg = SurfaceAdsorbateBuilder(
        element="Fe", facet="bcc110", size=(2, 2, 3), molecule_name="H2"
    )
    s = build_surface_adsorbate(cfg)
    assert s.provenance.origin == "generated"
    assert "surface_adsorbate" in s.provenance.source
    assert "Fe" in s.provenance.source


def test_surface_adsorbate_smiles(tmp_path):
    pytest.importorskip("rdkit", reason="rdkit not installed")
    cfg = SurfaceAdsorbateBuilder(
        element="Cu", facet="fcc111", size=(2, 2, 3), smiles="O"
    )
    s = build_surface_adsorbate(cfg)
    assert s.n_fragments == 1
    assert "fragment_smiles" in s.provenance.extra


def test_surface_adsorbate_file(tmp_path):
    co_file = tmp_path / "co.xyz"
    co_file.write_text("2\nCO\nC 0.0 0.0 0.0\nO 0.0 0.0 1.2\n")
    cfg = SurfaceAdsorbateBuilder(
        element="Cu", facet="fcc111", size=(2, 2, 3), file=str(co_file)
    )
    s = build_surface_adsorbate(cfg)
    assert s.n_fragments == 1


# ---------------------------------------------------------------------------
# surface_packing

@pytest.mark.skipif(
    shutil.which("packmol") is None, reason="packmol binary not on PATH"
)
def test_surface_packing_n_fragments():
    cfg = SurfacePackingBuilder(
        element="Cu", facet="fcc111", size=(3, 3, 4),
        molecule_name="CO", n_molecules=3,
    )
    s = build_surface_packing(cfg)

    frag = get_fragments(s.atoms)
    assert frag is not None
    assert s.n_fragments == 3

    n_framework = int(np.sum(frag == FRAMEWORK))
    expected_slab = 3 * 3 * 4
    assert n_framework == expected_slab


@pytest.mark.skipif(
    shutil.which("packmol") is None, reason="packmol binary not on PATH"
)
def test_surface_packing_substrate_immobile():
    cfg = SurfacePackingBuilder(
        element="Cu", facet="fcc111", size=(3, 3, 4),
        molecule_name="CO", n_molecules=2,
    )
    s = build_surface_packing(cfg)
    frag = get_fragments(s.atoms)
    substrate_pos = s.atoms.get_positions()[frag == FRAMEWORK]
    # All substrate atoms should be below the packing region.
    assert substrate_pos.shape[0] == 3 * 3 * 4
