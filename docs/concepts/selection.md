# Selection Funnel

The selection funnel is TrainCraft's quality-control layer. It sits between
sampling and DFT labeling, ensuring that only geometrically valid and
informationally diverse frames get expensive computation.

---

## The funnel in the pipeline

```
   Sampling output (hundreds of frames)
            │
            ▼
   ┌─────────────────┐
   │  physicality    │  ← filter: min interatomic distance
   └────────┬────────┘
            │  survivors
            ▼
   ┌─────────────────┐
   │     dedup       │  ← filter: exact hash deduplication
   └────────┬────────┘
            │  survivors
            ▼
   ┌─────────────────┐
   │   diversity     │  ← select: farthest-point sampling
   └────────┬────────┘
            │  budget-capped selection
            ▼
   DFT labeling (or dataset append)
```

---

## Stage: `physicality`

**What it does:** drops any frame where the minimum interatomic distance falls
below `min_distance` (in Å).

**Why:** MD runs occasionally produce bad steps — atoms fly through each other,
constraints fail, or the timestep was too large. These frames have unphysical
geometries that would poison the training set.

**Implementation:** `O(n²)` distance check using ASE's distance matrix.
Fast enough for any practical system size.

```python
from traincraft.selection.physicality import PhysicalitySelector

sel = PhysicalitySelector()
good = sel.select(frames, cfg)   # drops frames with d_min < cfg.min_distance
```

---

## Stage: `dedup`

**What it does:** removes structurally identical frames using a content hash.

**The hash:** `sha1(atomic_numbers + positions_rounded_4dp + cell_rounded_4dp + pbc)`

**Why:** when you run the pipeline multiple times (different temperatures,
different seeds, or restarting from a checkpoint), you accumulate duplicates.
Training on duplicate frames wastes compute and can bias the model.

Dedup also runs at the `Dataset` level — appending to a dataset never creates
duplicate entries.

---

## Stage: `diversity`

**What it does:** selects the `budget` most diverse frames using
**farthest-point sampling (FPS)** over a descriptor.

**The descriptor:** currently a composition + pair-distance histogram:
- Atom counts per element
- Binned histogram of all pairwise distances

This is cheap to compute (no SOAP/ACE expansion) and still captures the key
structural features that distinguish frames.

**FPS algorithm:**

1. Start with a random seed frame.
2. At each step, pick the frame that is *farthest* in descriptor space from
   all already-selected frames.
3. Repeat until `budget` frames are selected.

This guarantees that the selected set covers descriptor space as uniformly as
possible — you avoid spending DFT budget on frames that are all "the same".

```
Before diversity:  ● ● ● ●●●● ●●● ● ●●
                   (clumped in low-energy region)

After FPS:         ★      ★      ★     ★
                   (evenly spread across configuration space)
```

---

## Selector protocol

All selectors implement a simple protocol:

```python
class Selector(Protocol):
    def select(
        self,
        frames: list[Structure],
        cfg: SelectionConfig,
    ) -> list[Structure]: ...
```

Registered with:
```python
@register("selector", "diversity")
class DiversitySelector:
    def select(self, frames, cfg): ...
```

---

## Configuring the funnel

```toml
[selection]
steps        = ["physicality", "dedup", "diversity"]
budget       = 50
min_distance = 0.8
```

- **Order matters.** Run cheap filters first (`physicality`, `dedup`) to reduce
  the input to `diversity` — FPS is O(n²) in the surviving frame count.
- **`budget`** is applied at the very end. If fewer frames survive the filters,
  all are kept.
- You can **drop any stage** by removing it from `steps`. An empty `steps = []`
  passes all frames through unchanged.
- Set `budget = null` to remove the budget cap.

---

## Extending the funnel

Add a custom selector by registering it:

```python
from traincraft.core import register
from traincraft.core.structure import Structure
from traincraft.config.models import SelectionConfig

@register("selector", "energy_window")
class EnergyWindowSelector:
    def select(self, frames: list[Structure], cfg: SelectionConfig) -> list[Structure]:
        lo, hi = cfg.extra.get("e_min", -inf), cfg.extra.get("e_max", inf)
        return [f for f in frames if lo <= f.properties.get("energy", 0) <= hi]
```

Then add `"energy_window"` to `steps` in your TOML.
