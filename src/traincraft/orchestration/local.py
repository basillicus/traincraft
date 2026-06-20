"""Local engine: serial wiring of the pipeline.

This is the seam where a QuACC/Covalent adapter later replaces the loop with a
parallel DAG. The science lives in pure functions; this only orders them.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from ..calculators import make_calculator
from ..config import TrainCraftConfig
from ..core import Workspace
from ..datasets import Dataset, write_frames
from ..geometry import build_geometry
from ..labeling import label_frames
from ..sampling import run_sampling
from ..selection import run_funnel

logger = logging.getLogger(__name__)


def run_pipeline(config: TrainCraftConfig) -> dict:
    if config.geometry is None:
        raise ValueError("config has no [geometry] section: nothing to build")

    ws = Workspace(Path(config.run.outdir) / config.run.name)
    if config.run.seed is not None:
        np.random.seed(config.run.seed)

    # --- geometry ---------------------------------------------------------
    structure = build_geometry(config.geometry)
    write_frames(ws.subdir("structures") / "initial.extxyz", [structure])
    candidates = [structure]

    # --- sampling ---------------------------------------------------------
    if config.sampling is not None:
        if config.calculator is None:
            raise ValueError("[sampling] needs a [calculator]")
        calc = make_calculator(config.calculator)
        job = ws.job("candidates")
        candidates = run_sampling(structure, calc, job, config.sampling)
        write_frames(job.path("candidates.extxyz"), candidates)

    # --- selection (the funnel, before any expensive labeling) -----------
    selected = candidates
    if config.selection is not None:
        selected = run_funnel(candidates, config.selection)
        write_frames(ws.subdir("selected") / "selected.extxyz", selected)

    # --- labeling (DFT on the selected frames) ----------------------------
    labeled = selected
    if config.labeling is not None:
        job = ws.job("labeled_dft")
        labeled = label_frames(selected, config.labeling.calculator, out_dir=job.dir)
        write_frames(job.path("labeled.extxyz"), labeled)
        job.mark_done()
        logger.info("labeled %d frames -> %s", len(labeled), job.dir)

    # --- dataset ----------------------------------------------------------
    dataset_path = None
    if config.dataset is not None:
        ds = Dataset(ws.root / config.dataset.path)
        ds.append(labeled)
        dataset_path = str(ds.write())
        logger.info("dataset written: %s (%d frames)", ds.path, len(ds))

    summary = {
        "workspace": str(ws.root),
        "n_candidates": len(candidates),
        "n_selected": len(selected),
        "n_labeled": len(labeled) if config.labeling is not None else 0,
        "dataset": dataset_path,
    }
    return summary
