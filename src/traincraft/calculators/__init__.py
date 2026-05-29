"""Calculators: cheap potentials (force field / semiempirical / MLIP)."""

from __future__ import annotations

from . import potentials  # noqa: F401  (registers built-in calculators)
from .base import CalculatorFactory, make_calculator

__all__ = ["CalculatorFactory", "make_calculator"]
