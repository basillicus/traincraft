"""Calculator factory protocol + resolver.

A factory takes a calculator config model and returns an ASE calculator. Each is
registered under kind ``"calculator"`` with declared ``capabilities`` (the
properties it can produce), so a workflow can be validated up front.
"""

from __future__ import annotations

from typing import Protocol

from ..core import get


class CalculatorFactory(Protocol):
    def __call__(self, cfg) -> object:  # returns an ase.calculators Calculator
        ...


def make_calculator(cfg):
    """Resolve and build the ASE calculator named by ``cfg.type``."""
    factory = get("calculator", cfg.type)
    return factory(cfg)
