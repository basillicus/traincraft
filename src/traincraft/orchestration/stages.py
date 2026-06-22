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

**Staleness (smart caching).** Every artifact gets a sidecar ``<artifact>.key``
holding a hash of the inputs that produced it: the stage's own config section
(plus the run seed where it matters) **and the upstream stage's key**. On a
rerun a stage is skipped only if its artifact exists *and* that key still
matches. Because each key folds in the upstream key, editing the geometry config
changes the geometry key, which changes the sample key, which changes the select
key, and so on — so exactly the affected stages (and everything downstream)
recompute, while unrelated/expensive stages stay cached. ``force=True`` ignores
the keys and recomputes regardless.
"""

from __future__ import annotations

import hashlib
import json
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


def _latest_upstream_path(ws: Workspace, stage: str) -> Path | None:
    """Path of the most recent existing producer artifact before ``stage``."""
    idx = _PRODUCERS.index(stage) if stage in _PRODUCERS else len(_PRODUCERS)
    for upstream in reversed(_PRODUCERS[:idx]):
        path = _artifact_path(ws, upstream)
        if path.exists():
            return path
    return None


def _latest_upstream(ws: Workspace, stage: str):
    """Frames from the most recent existing artifact before ``stage``."""
    path = _latest_upstream_path(ws, stage)
    if path is None:
        raise FileNotFoundError(
            f"stage {stage!r} has no upstream artifact in {ws.root}; run earlier stages first"
        )
    return read_frames(path)


# --- staleness keys ("receipts") -------------------------------------------

def _dump(model) -> dict | None:
    """A JSON-able, deterministic snapshot of a config sub-model (or None)."""
    return model.model_dump(mode="json") if model is not None else None


def _hash(payload) -> str:
    blob = json.dumps(payload, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def _stage_inputs(config, stage: str) -> dict:
    """The config that determines a stage's output (its own section + seed)."""
    seed = config.run.seed
    return {
        "geometry": {"geometry": _dump(config.geometry), "seed": seed},
        "sample": {
            "sampling": _dump(config.sampling),
            "calculator": _dump(config.calculator),
            "seed": seed,
        },
        "select": {"selection": _dump(config.selection), "seed": seed},
        "label": {"labeling": _dump(config.labeling)},
        "dataset": {"dataset": _dump(config.dataset)},
        "train": {"training": _dump(config.training)},
    }[stage]


def _key_file(artifact: Path) -> Path:
    return artifact.with_name(artifact.name + ".key")


def _read_key(artifact: Path) -> str | None:
    kf = _key_file(artifact)
    return kf.read_text().strip() if kf.exists() else None


def _write_key(artifact: Path, key: str) -> None:
    _key_file(artifact).write_text(key)


def _upstream_key(config, ws: Workspace, stage: str) -> str | None:
    """Key of the artifact this stage consumes (None if it has no producer yet)."""
    if stage == "train" and config.dataset is not None:
        return _read_key(Dataset(ws.root / config.dataset.path).path)
    path = _latest_upstream_path(ws, stage)
    return _read_key(path) if path is not None else None


def _stage_key(config, ws: Workspace, stage: str) -> str:
    """Hash of this stage's own inputs folded together with its upstream key."""
    return _hash({"self": _stage_inputs(config, stage), "up": _upstream_key(config, ws, stage)})


def _fresh(artifact: Path, key: str) -> bool:
    """True if the artifact exists and was produced with this exact key."""
    return artifact.exists() and _read_key(artifact) == key


# --- individual stages ------------------------------------------------------

def stage_geometry(config, ws, force=False):
    out = _artifact_path(ws, "geometry")
    key = _stage_key(config, ws, "geometry")
    if not force and _fresh(out, key):
        logger.info("geometry: up to date, skipping (%s)", out)
        return read_frames(out)
    if config.geometry is None:
        raise ValueError("config has no [geometry] section: nothing to build")
    structure = build_geometry(config.geometry)
    ws.subdir("structures")
    write_frames(out, [structure])
    _write_key(out, key)
    return [structure]


def stage_sample(config, ws, force=False):
    if config.sampling is None:
        return _latest_upstream(ws, "sample")
    out = _artifact_path(ws, "sample")
    key = _stage_key(config, ws, "sample")
    if not force and _fresh(out, key):
        logger.info("sample: up to date, skipping (%s)", out)
        return read_frames(out)
    if config.calculator is None:
        raise ValueError("[sampling] needs a [calculator]")
    structure = _latest_upstream(ws, "sample")[0]
    calc = make_calculator(config.calculator)
    job = ws.job("candidates")
    candidates = run_sampling(structure, calc, job, config.sampling)
    write_frames(out, candidates)
    _write_key(out, key)
    return candidates


def stage_select(config, ws, force=False):
    if config.selection is None:
        return _latest_upstream(ws, "select")
    out = _artifact_path(ws, "select")
    key = _stage_key(config, ws, "select")
    if not force and _fresh(out, key):
        logger.info("select: up to date, skipping (%s)", out)
        return read_frames(out)
    selected = run_funnel(_latest_upstream(ws, "select"), config.selection)
    ws.subdir("selected")
    write_frames(out, selected)
    _write_key(out, key)
    return selected


def stage_label(config, ws, force=False):
    if config.labeling is None:
        return _latest_upstream(ws, "label")
    out = _artifact_path(ws, "label")
    key = _stage_key(config, ws, "label")
    if not force and _fresh(out, key):
        logger.info("label: up to date, skipping (%s)", out)
        return read_frames(out)
    job = ws.job("labeled_dft")
    frames = _latest_upstream(ws, "label")
    labeled = label_frames(frames, config.labeling.calculator, out_dir=job.dir)
    write_frames(out, labeled)
    _write_key(out, key)
    return labeled


def stage_dataset(config, ws, force=False):
    frames = _latest_upstream(ws, "dataset")
    if config.dataset is None:
        return frames
    ds = Dataset(ws.root / config.dataset.path)
    key = _stage_key(config, ws, "dataset")
    if not force and _fresh(ds.path, key):
        logger.info("dataset: up to date, skipping (%s)", ds.path)
        return frames
    ds.append(frames)
    ds.write()
    _write_key(ds.path, key)
    logger.info("dataset written: %s (%d frames)", ds.path, len(ds))
    return frames


def stage_train(config, ws, force=False):
    if config.training is None:
        return _latest_upstream(ws, "train")
    job = ws.job("model")
    manifest = job.dir / "manifest.json"
    key = _stage_key(config, ws, "train")
    if not force and _fresh(manifest, key):
        logger.info("train: up to date, skipping (%s)", manifest)
        return _training_frames(config, ws)
    frames = _training_frames(config, ws)
    result = run_training(frames, config.training, job)
    _write_key(manifest, key)
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
