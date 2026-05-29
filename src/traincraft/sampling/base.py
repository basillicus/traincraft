"""Sampler protocol + resolver.

A sampler turns one input structure into many candidate frames using a cheap
calculator. It writes artifacts into the job's absolute directory and returns
``list[Structure]`` — it does NOT decide what gets DFT-labeled (that is the
selection layer's job).
"""

from __future__ import annotations

from typing import Protocol

from ..core import Job, Structure, get


class Sampler(Protocol):
    def __call__(self, structure: Structure, calc, job: Job, cfg) -> list[Structure]:
        ...


def run_sampling(structure: Structure, calc, job: Job, cfg) -> list[Structure]:
    return get("sampler", cfg.type)(structure, calc, job, cfg)
