"""Orchestration engines. ``local`` now; QuACC/Covalent adapters later."""

from __future__ import annotations

from .local import run_pipeline

__all__ = ["run_pipeline"]
