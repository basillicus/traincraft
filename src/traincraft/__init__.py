"""TrainCraft: modular MLIP dataset generation & MACE-centric active learning.

Importing the package registers all built-in plugins (sources, builders,
transforms, calculators, samplers, selectors) and exposes a curated API so the
pieces are usable standalone in your own scripts — the CLI is just a thin shell.
"""

from __future__ import annotations

from .calculators import make_calculator
from .config import TrainCraftConfig, load_config, loads_config
from .core import (
    Provenance,
    Result,
    Structure,
    Workspace,
    available,
    get,
    make_rng,
    register,
)
from .datasets import Dataset, read_frames, write_frames
from .geometry import build_geometry
from .orchestration import run_pipeline
from .sampling import run_sampling
from .selection import run_funnel

__version__ = "0.1.0"

__all__ = [
    "Dataset",
    "Provenance",
    "Result",
    "Structure",
    "TrainCraftConfig",
    "Workspace",
    "available",
    "build_geometry",
    "get",
    "load_config",
    "loads_config",
    "make_calculator",
    "make_rng",
    "read_frames",
    "register",
    "run_funnel",
    "run_pipeline",
    "run_sampling",
    "write_frames",
]
