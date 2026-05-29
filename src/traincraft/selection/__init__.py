"""Selection funnel + built-in selectors."""

from __future__ import annotations

from . import dedup, diversity, physicality  # noqa: F401  (registers selectors)
from .base import Selector, run_funnel

__all__ = ["Selector", "run_funnel"]
