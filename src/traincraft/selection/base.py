"""Selection funnel: run ordered selectors, cheap -> expensive.

This is the layer that stops the pipeline from spending DFT on redundant or
unphysical frames. Selectors are plugins resolved by name from ``cfg.steps``.
"""

from __future__ import annotations

from typing import Protocol

from ..core import Structure, get


class Selector(Protocol):
    def __call__(self, frames: list[Structure], cfg) -> list[Structure]:
        ...


def run_funnel(frames: list[Structure], cfg) -> list[Structure]:
    out = list(frames)
    for step in cfg.steps:
        out = get("selector", step)(out, cfg)
    if cfg.budget is not None and len(out) > cfg.budget:
        out = out[: cfg.budget]
    return out
