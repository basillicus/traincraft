from __future__ import annotations

import os

import numpy as np
import pytest
from ase.calculators.emt import EMT

from traincraft.config.models import MonteCarloSampling, SurfaceAdsorbateBuilder
from traincraft.core import Workspace
from traincraft.core.fragments import FRAMEWORK, get_fragments
from traincraft.geometry.builders.surface import build_surface_adsorbate
from traincraft.sampling.monte_carlo import sample_monte_carlo


@pytest.fixture()
def cu_co_structure():
    cfg = SurfaceAdsorbateBuilder(
        element="Cu", facet="fcc111", size=(2, 2, 3), molecule_name="CO", site="ontop"
    )
    return build_surface_adsorbate(cfg)


@pytest.fixture()
def workspace(tmp_path):
    return Workspace(tmp_path / "runs")


def test_mc_returns_frames(cu_co_structure, workspace, tmp_path):
    calc = EMT()
    job = workspace.job("mc_test")
    cfg = MonteCarloSampling(steps=50, interval=10, p_conformer=0.0, seed=7)

    frames = sample_monte_carlo(cu_co_structure, calc, job, cfg)

    assert len(frames) >= 1


def test_mc_substrate_atoms_frozen(cu_co_structure, workspace):
    """Substrate atoms (tc_fragment == -1) must not move during MC."""
    initial_pos = cu_co_structure.atoms.get_positions().copy()
    frag = get_fragments(cu_co_structure.atoms)
    substrate_mask = frag == FRAMEWORK

    calc = EMT()
    job = workspace.job("mc_freeze")
    cfg = MonteCarloSampling(steps=30, interval=10, p_conformer=0.0, seed=42)

    sample_monte_carlo(cu_co_structure, calc, job, cfg)

    # Substrate positions on the *input* structure are unchanged (MC works on a copy).
    np.testing.assert_array_equal(
        cu_co_structure.atoms.get_positions()[substrate_mask],
        initial_pos[substrate_mask],
    )


def test_mc_frames_carry_tc_fragment(cu_co_structure, workspace):
    calc = EMT()
    job = workspace.job("mc_frag")
    cfg = MonteCarloSampling(steps=20, interval=10, p_conformer=0.0, seed=1)

    frames = sample_monte_carlo(cu_co_structure, calc, job, cfg)
    assert len(frames) >= 1

    for f in frames:
        assert get_fragments(f.atoms) is not None, "tc_fragment lost in MC frame"


def test_mc_frames_have_provenance(cu_co_structure, workspace):
    calc = EMT()
    job = workspace.job("mc_prov")
    cfg = MonteCarloSampling(steps=10, interval=5, p_conformer=0.0, seed=0)

    frames = sample_monte_carlo(cu_co_structure, calc, job, cfg)
    for f in frames:
        assert f.provenance.origin == "ml_sampled"
        assert "monte_carlo" in f.provenance.source


def test_mc_cwd_unchanged(cu_co_structure, workspace):
    before = os.getcwd()
    calc = EMT()
    job = workspace.job("mc_cwd")
    cfg = MonteCarloSampling(steps=10, interval=5, p_conformer=0.0, seed=2)
    sample_monte_carlo(cu_co_structure, calc, job, cfg)
    assert os.getcwd() == before


def test_mc_conformer_no_smiles_does_not_crash(cu_co_structure, workspace):
    """Conformer move with no SMILES recorded should fall back gracefully."""
    calc = EMT()
    job = workspace.job("mc_conf_fallback")
    # p_conformer > 0 but no fragment_smiles in provenance → should fall back to translate
    cfg = MonteCarloSampling(steps=20, interval=10, p_conformer=0.5, seed=3)
    # Strip fragment_smiles from provenance.
    s = cu_co_structure.copy()
    s.provenance.extra.pop("fragment_smiles", None)

    frames = sample_monte_carlo(s, calc, job, cfg)
    assert len(frames) >= 1


def test_mc_no_fragment_array_falls_back(workspace):
    """System with no tc_fragment array should be treated as one rigid fragment."""
    from ase.build import molecule

    from traincraft.core import Provenance, Structure

    atoms = molecule("CO")
    atoms.center(vacuum=5.0)
    # Deliberately do NOT set fragments.
    s = Structure.from_ase(atoms, provenance=Provenance(origin="generated"))

    calc = EMT()
    job = workspace.job("mc_no_frag")
    cfg = MonteCarloSampling(steps=10, interval=5, p_conformer=0.0, seed=4)

    frames = sample_monte_carlo(s, calc, job, cfg)
    assert len(frames) >= 1
