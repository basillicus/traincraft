"""Training: turn a labelled dataset into a trained MLIP (MACE-centric).

Importing the package registers the built-in trainer backends. ``run_training``
resolves the backend from config and trains; ``write_training_xyz`` is the shared
extxyz exporter (TrainCraft ``tc_*`` properties → MACE ``REF_*`` keys).
"""

from __future__ import annotations

from . import mace  # noqa: F401  (import registers the "mace" trainer)
from .base import TrainResult, run_training, write_training_xyz

__all__ = ["TrainResult", "run_training", "write_training_xyz"]
