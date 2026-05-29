"""Sampler: Metropolis Monte Carlo (interface placeholder).

Designed signature is molecule-aware so the full implementation (rigid-body
translation/rotation of molecules + RDKit conformer generation, the primary tool
for molecules on surfaces) lands in the next chunk without changing the API.
"""

from __future__ import annotations

from ..core import Job, Structure, register


@register("sampler", "monte_carlo")
def sample_monte_carlo(structure: Structure, calc, job: Job, cfg) -> list[Structure]:
    raise NotImplementedError(
        "Monte Carlo sampler (rigid-body moves + RDKit conformers) lands in the next chunk"
    )
