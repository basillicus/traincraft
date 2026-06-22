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
    stage_train,
    workspace_for,
)

logger = logging.getLogger(__name__)


def run_pipeline(config: TrainCraftConfig, *, force: bool = False) -> dict:
    if config.geometry is None:
        raise ValueError("config has no [geometry] section: nothing to build")

    ws = workspace_for(config)
    if config.run.seed is not None:
        np.random.seed(config.run.seed)

    stage_geometry(config, ws, force=force)
    candidates = stage_sample(config, ws, force=force)
    selected = stage_select(config, ws, force=force)
    labeled = stage_label(config, ws, force=force)
    stage_dataset(config, ws, force=force)
    stage_train(config, ws, force=force)

    dataset_path = None
    if config.dataset is not None:
        dataset_path = str(Dataset(ws.root / config.dataset.path).path)
    model_dir = str(ws.root / "model") if config.training is not None else None

    return {
        "workspace": str(ws.root),
        "n_candidates": len(candidates),
        "n_selected": len(selected),
        "n_labeled": len(labeled) if config.labeling is not None else 0,
        "dataset": dataset_path,
        "model": model_dir,
    }
