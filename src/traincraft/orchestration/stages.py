"""Pipeline as discrete, resumable stages over workspace artifacts.

Each stage reads the latest upstream artifact, does its work, and writes its own
artifact — so a stage can run **standalone** (e.g. inside a container on an HPC
node) and the next stage picks up from disk. This is what the Slurm executor
dispatches as separate, dependency-chained jobs; the local engine just calls them
in order in one process.

Artifacts (under the run workspace):
    geometry -> structures/initial.extxyz
    sample   -> candidates/candidates.extxyz
    select   -> selected/selected.extxyz
    label    -> labeled_dft/labeled.extxyz  (+ manifest.json, frame_*/)
    dataset  -> <dataset.path>/ (a Dataset)
    train    -> model/ (mace_run_train: <name>.model + manifest.json)
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from ..calculators import make_calculator
from ..config import TrainCraftConfig
from ..core import Workspace
from ..datasets import Dataset, read_frames, write_frames
from ..geometry import build_geometry
from ..labeling import label_frames
from ..sampling import run_sampling
from ..selection import run_funnel
from ..training import run_training

logger = logging.getLogger(__name__)

# Stages in pipeline order. dataset is terminal for frames; train consumes the
# dataset (or labelled frames) and emits a model rather than a frames artifact.
STAGE_ORDER = ("geometry", "sample", "select", "label", "dataset", "train")
_PRODUCERS = ("geometry", "sample", "select", "label")
_ARTIFACT = {
    "geometry": ("structures", "initial.extxyz"),
    "sample": ("candidates", "candidates.extxyz"),
    "select": ("selected", "selected.extxyz"),
    "label": ("labeled_dft", "labeled.extxyz"),
}


def workspace_for(config: TrainCraftConfig) -> Workspace:
    return Workspace(Path(config.run.outdir) / config.run.name)


def enabled_stages(config: TrainCraftConfig) -> list[str]:
    """Stages actually requested by the config (geometry is always present)."""
    stages = ["geometry"]
    for name, section in (
        ("sample", config.sampling),
        ("select", config.selection),
        ("label", config.labeling),
        ("dataset", config.dataset),
        ("train", config.training),
    ):
        if section is not None:
            stages.append(name)
    return stages


def _artifact_path(ws: Workspace, stage: str) -> Path:
    sub, name = _ARTIFACT[stage]
    return ws.root / sub / name


def _latest_upstream(ws: Workspace, stage: str):
    """Frames from the most recent existing artifact before ``stage``."""
    idx = _PRODUCERS.index(stage) if stage in _PRODUCERS else len(_PRODUCERS)
    for upstream in reversed(_PRODUCERS[:idx]):
        path = _artifact_path(ws, upstream)
        if path.exists():
            return read_frames(path)
    raise FileNotFoundError(
        f"stage {stage!r} has no upstream artifact in {ws.root}; run earlier stages first"
    )


# --- individual stages ------------------------------------------------------

def stage_geometry(config, ws, force=False):
    out = _artifact_path(ws, "geometry")
    if out.exists() and not force:
        logger.info("geometry: using cached %s", out)
        return read_frames(out)
    if config.geometry is None:
        raise ValueError("config has no [geometry] section: nothing to build")
    structure = build_geometry(config.geometry)
    ws.subdir("structures")
    write_frames(out, [structure])
    return [structure]


def stage_sample(config, ws, force=False):
    if config.sampling is None:
        return _latest_upstream(ws, "sample")
    out = _artifact_path(ws, "sample")
    if out.exists() and not force:
        logger.info("sample: using cached %s", out)
        return read_frames(out)
    if config.calculator is None:
        raise ValueError("[sampling] needs a [calculator]")
    structure = _latest_upstream(ws, "sample")[0]
    calc = make_calculator(config.calculator)
    job = ws.job("candidates")
    candidates = run_sampling(structure, calc, job, config.sampling)
    write_frames(out, candidates)
    return candidates


def stage_select(config, ws, force=False):
    if config.selection is None:
        return _latest_upstream(ws, "select")
    out = _artifact_path(ws, "select")
    if out.exists() and not force:
        logger.info("select: using cached %s", out)
        return read_frames(out)
    selected = run_funnel(_latest_upstream(ws, "select"), config.selection)
    ws.subdir("selected")
    write_frames(out, selected)
    return selected


def stage_label(config, ws, force=False):
    if config.labeling is None:
        return _latest_upstream(ws, "label")
    out = _artifact_path(ws, "label")
    if out.exists() and not force:
        logger.info("label: using cached %s", out)
        return read_frames(out)
    job = ws.job("labeled_dft")
    frames = _latest_upstream(ws, "label")
    labeled = label_frames(frames, config.labeling.calculator, out_dir=job.dir)
    write_frames(out, labeled)
    return labeled


def stage_dataset(config, ws, force=False):
    frames = _latest_upstream(ws, "dataset")
    if config.dataset is None:
        return frames
    ds = Dataset(ws.root / config.dataset.path)
    ds.append(frames)
    ds.write()
    logger.info("dataset written: %s (%d frames)", ds.path, len(ds))
    return frames


def stage_train(config, ws, force=False):
    if config.training is None:
        return _latest_upstream(ws, "train")
    job = ws.job("model")
    manifest = job.dir / "manifest.json"
    if manifest.exists() and not force:
        logger.info("train: using cached %s", manifest)
        return _training_frames(config, ws)
    frames = _training_frames(config, ws)
    result = run_training(frames, config.training, job)
    logger.info(
        "trained '%s' on %d frames (model: %s)",
        config.training.name, result.n_train, result.model_path,
    )
    return frames


def _training_frames(config, ws) -> list:
    """Frames to train on: the dataset artifact if built, else the latest upstream."""
    if config.dataset is not None:
        ds_path = Dataset(ws.root / config.dataset.path).path
        if ds_path.exists():
            return read_frames(ds_path)
    return _latest_upstream(ws, "train")


_DISPATCH = {
    "geometry": stage_geometry,
    "sample": stage_sample,
    "select": stage_select,
    "label": stage_label,
    "dataset": stage_dataset,
    "train": stage_train,
}


def run_stage(name: str, config, ws=None, *, force=False):
    """Run a single named stage (seeding RNG so standalone runs are reproducible)."""
    if name not in _DISPATCH:
        raise ValueError(f"unknown stage {name!r}; choose from {STAGE_ORDER}")
    ws = ws or workspace_for(config)
    if config.run.seed is not None:
        np.random.seed(config.run.seed)
    return _DISPATCH[name](config, ws, force=force)
