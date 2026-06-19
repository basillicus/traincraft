# Write a Custom Sampler

The `Sampler` protocol is straightforward: accept a structure, a calculator,
a job directory, and a sampling config; return a list of structures.

---

## Step 1 — Write the sampler

Create `src/traincraft/sampling/my_sampler.py`:

```python
from __future__ import annotations

import logging

from ..core import Provenance, Structure, register

logger = logging.getLogger(__name__)


@register("sampler", "my_sampler")
def run_my_sampler(
    structure: Structure,
    calc,               # ASE calculator
    job,                # traincraft.core.Job
    cfg,                # config model instance
) -> list[Structure]:
    """Generate structures by simple random displacement."""
    import numpy as np
    rng = np.random.default_rng(cfg.seed)

    frames = []
    for i in range(cfg.n_structures):
        s = structure.copy()
        s.atoms.positions += rng.normal(0, cfg.std, s.atoms.positions.shape)
        # optionally: run a single-point calculation
        # s.atoms.calc = calc
        # s.properties["energy"] = s.atoms.get_potential_energy()
        prov = s.provenance
        prov.transforms.append(f"my_sampler:step={i}")
        frames.append(s)

    logger.info("my_sampler: generated %d frames", len(frames))
    return frames
```

---

## Step 2 — Config model

Add to `src/traincraft/config/models.py`:

```python
class MySampling(TCModel):
    type: Literal["my_sampler"] = "my_sampler"
    n_structures: int = 10
    std: float = 0.05
    seed: int | None = None
```

Add `MySampling` to the `SamplingConfig` union.

---

## Step 3 — Register on import

Add to `src/traincraft/sampling/__init__.py`:

```python
from . import md, monte_carlo, my_sampler, rattle
```

---

## Step 4 — Use in TOML

```toml
[sampling]
type         = "my_sampler"
n_structures = 20
std          = 0.08
seed         = 42
```

---

## Tips

!!! tip "Use the Job directory for large outputs"
    Write trajectory files into `job.path("trajectory.traj")` — the `Job`
    object provides an absolute path inside the run's workspace directory.

!!! tip "Fragment-aware sampling"
    Access the fragment array via `structure.fragments` to implement moves
    that respect the substrate/adsorbate division, matching the MC sampler pattern.
