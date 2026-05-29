"""Seeded RNG factory.

Replaces the legacy global ``set_seed``: a seed is threaded explicitly into the
jobs that need it, so parallel tasks have independent, reproducible streams.
"""

from __future__ import annotations

import numpy as np


def make_rng(seed: int | None = None) -> np.random.Generator:
    return np.random.default_rng(seed)
