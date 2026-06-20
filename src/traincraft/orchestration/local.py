"""Local engine: serial wiring of the pipeline.

The science lives in pure stage functions (see ``stages.py``); this only orders
them in one process. The Slurm executor (``slurm.py``) wires the same stages as
dependency-chained jobs — identical science, swappable engine.
"""

from __future__ import annotations

import logging

import numpy as np

from ..config import TrainCraftConfig
from ..datasets import Dataset
from .stages import (
    stage_dataset,
    stage_geometry,
    stage_label,
    stage_sample,
    stage_select,
    workspace_for,
)

logger = logging.getLogger(__name__)


def run_pipeline(config: TrainCraftConfig) -> dict:
    if config.geometry is None:
        raise ValueError("config has no [geometry] section: nothing to build")

    ws = workspace_for(config)
    if config.run.seed is not None:
        np.random.seed(config.run.seed)

    stage_geometry(config, ws)
    candidates = stage_sample(config, ws)
    selected = stage_select(config, ws)
    labeled = stage_label(config, ws)
    stage_dataset(config, ws)

    dataset_path = None
    if config.dataset is not None:
        dataset_path = str(Dataset(ws.root / config.dataset.path).path)

    return {
        "workspace": str(ws.root),
        "n_candidates": len(candidates),
        "n_selected": len(selected),
        "n_labeled": len(labeled) if config.labeling is not None else 0,
        "dataset": dataset_path,
    }
