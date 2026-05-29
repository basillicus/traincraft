from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.build import molecule

from traincraft import Structure
from traincraft.config.models import SelectionConfig
from traincraft.selection import run_funnel

_rng = np.random.default_rng(0)


def _perturbed(scale: float) -> Structure:
    a = molecule("H2O")
    a.positions = a.positions + _rng.normal(0.0, scale, a.positions.shape)
    return Structure.from_ase(a)


def test_dedup_collapses_identical():
    a = molecule("H2O")
    frames = [Structure.from_ase(a) for _ in range(3)]
    out = run_funnel(frames, SelectionConfig(steps=["dedup"], budget=None))
    assert len(out) == 1


def test_physicality_filters_overlapping_atoms():
    good = Structure.from_ase(molecule("H2O"))
    bad = Structure.from_ase(Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.1]]))
    out = run_funnel(
        [good, bad], SelectionConfig(steps=["physicality"], budget=None, min_distance=0.7)
    )
    assert len(out) == 1


def test_diversity_reduces_to_budget():
    frames = [_perturbed(0.15) for _ in range(10)]
    out = run_funnel(frames, SelectionConfig(steps=["diversity"], budget=3))
    assert len(out) == 3
