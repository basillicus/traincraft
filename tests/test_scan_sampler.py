"""Tests for the deterministic grid `scan` sampler."""

from __future__ import annotations

import numpy as np
import pytest

from traincraft.config.models import ScanAxis, ScanSampling, SurfaceAdsorbateBuilder
from traincraft.core import Workspace
from traincraft.core.fragments import FRAMEWORK, get_fragments
from traincraft.geometry.builders.surface import build_surface_adsorbate
from traincraft.sampling.scan import sample_scan


@pytest.fixture()
def cu_co_structure():
    cfg = SurfaceAdsorbateBuilder(
        element="Cu", facet="fcc111", size=(2, 2, 3), molecule_name="CO", site="ontop"
    )
    return build_surface_adsorbate(cfg)


@pytest.fixture()
def job(tmp_path):
    return Workspace(tmp_path / "runs").job("scan")


# --- config validation ------------------------------------------------------

def test_scan_needs_a_grid():
    with pytest.raises(ValueError, match="at least one of"):
        ScanSampling()  # neither translate nor rotate


def test_scan_axis_steps_positive():
    with pytest.raises(ValueError, match="steps >= 1"):
        ScanAxis(start=0.0, stop=1.0, steps=0)


# --- enumeration ------------------------------------------------------------

def test_translate_grid_count_and_displacement(cu_co_structure, job):
    cfg = ScanSampling(fragment=0, translate=ScanAxis(axis="z", start=0.0, stop=2.0, steps=5))
    frames = sample_scan(cu_co_structure, None, job, cfg)

    assert len(frames) == 5  # 5 grid points
    frag = get_fragments(cu_co_structure.atoms)
    ads = frag == 0
    base_z = cu_co_structure.atoms.get_positions()[ads, 2]
    # last frame is the adsorbate lifted by +2.0 Å in z
    lifted_z = frames[-1].atoms.get_positions()[ads, 2]
    np.testing.assert_allclose(lifted_z, base_z + 2.0, atol=1e-9)


def test_cartesian_product_of_translate_and_rotate(cu_co_structure, job):
    cfg = ScanSampling(
        fragment=0,
        translate=ScanAxis(axis="z", start=0.0, stop=1.0, steps=3),
        rotate=ScanAxis(axis="x", start=0.0, stop=180.0, steps=4),
    )
    frames = sample_scan(cu_co_structure, None, job, cfg)
    assert len(frames) == 3 * 4  # full grid


def test_substrate_never_moves(cu_co_structure, job):
    frag = get_fragments(cu_co_structure.atoms)
    sub = frag == FRAMEWORK
    base = cu_co_structure.atoms.get_positions()[sub].copy()

    cfg = ScanSampling(fragment=0, translate=ScanAxis(start=0.0, stop=3.0, steps=4))
    for f in sample_scan(cu_co_structure, None, job, cfg):
        np.testing.assert_array_equal(f.atoms.get_positions()[sub], base)


def test_frames_keep_fragment_tags(cu_co_structure, job):
    cfg = ScanSampling(fragment=0, rotate=ScanAxis(axis="z", start=0.0, stop=90.0, steps=2))
    frames = sample_scan(cu_co_structure, None, job, cfg)
    for f in frames:
        assert get_fragments(f.atoms) is not None


def test_unknown_fragment_raises(cu_co_structure, job):
    cfg = ScanSampling(fragment=9, translate=ScanAxis(start=0.0, stop=1.0, steps=2))
    with pytest.raises(ValueError, match="no atoms; available mobile fragments"):
        sample_scan(cu_co_structure, None, job, cfg)


def test_registered():
    from traincraft.core import get

    assert get("sampler", "scan") is not None
