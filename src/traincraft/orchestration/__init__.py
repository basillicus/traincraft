"""Orchestration engines. ``local`` + ``slurm`` now; QuACC/Covalent later."""

from __future__ import annotations

from .local import run_pipeline
from .slurm import render as render_slurm
from .slurm import submit as submit_slurm
from .stages import STAGE_ORDER, enabled_stages, run_stage, workspace_for

__all__ = [
    "STAGE_ORDER",
    "enabled_stages",
    "render_slurm",
    "run_pipeline",
    "run_stage",
    "submit_slurm",
    "workspace_for",
]
