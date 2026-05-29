"""Transform: Gaussian perturbation of atomic positions."""

from __future__ import annotations

import numpy as np

from ...core import Structure, register


@register("transform", "perturb")
def transform_perturb(
    structure: Structure, cfg, rng: np.random.Generator | None = None
) -> Structure:
    rng = rng if rng is not None else np.random.default_rng()
    out = structure.copy()
    out.atoms.positions = out.atoms.positions + rng.normal(
        0.0, cfg.stddev, out.atoms.positions.shape
    )
    out.properties = {}  # geometry changed: stale properties dropped
    out.provenance.transforms.append(f"perturb:{cfg.stddev}")
    return out
