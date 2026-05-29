"""Sampler: HiPhive rattle (ported from legacy ``gen_rattled_geometries``).

HiPhive is an optional dependency; install via the ``sampling`` extra/feature.
"""

from __future__ import annotations

import numpy as np

from ..core import Job, Provenance, Structure, register


@register("sampler", "rattle")
def sample_rattle(structure: Structure, calc, job: Job, cfg) -> list[Structure]:
    try:
        from hiphive.structure_generation import (
            generate_mc_rattled_structures,
            generate_rattled_structures,
        )
    except ImportError as exc:  # pragma: no cover - optional dep
        raise ImportError(
            "rattle sampling needs HiPhive: install the 'sampling' extra/feature"
        ) from exc

    atoms = structure.atoms
    frames: list[Structure] = []
    for _ in range(cfg.n_structures):
        seed = int(np.random.randint(2**31 - 1))
        if cfg.method == "mc":
            rattled = generate_mc_rattled_structures(
                atoms, 1, 0.25 * cfg.std, cfg.min_distance, seed=seed
            )
        else:
            rattled = generate_rattled_structures(atoms, 1, cfg.std, seed=seed)
        for a in rattled:
            frames.append(
                Structure.from_ase(
                    a,
                    provenance=Provenance(
                        origin="ml_sampled",
                        source=f"sampler:rattle:{cfg.method}",
                        parents=[structure.hash],
                    ),
                )
            )
    return frames
