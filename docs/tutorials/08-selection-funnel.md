# Tutorial 8 · The Selection Funnel

**What you'll learn:** what each stage of the selection funnel does, how to
configure it, and when to use which combination of stages.

**Prerequisites:** Tutorial 1. No extra dependencies.

---

## Why a funnel?

Generating thousands of candidate frames with MD or MC is cheap. Labeling
them with DFT is expensive (often hours per frame on an HPC cluster). The
selection funnel is the bridge: it removes bad and redundant frames *before*
any DFT is dispatched.

```
raw frames
  → 1. physicality  (< 1 ms/frame — drops unphysical geometries)
  → 2. dedup        (< 1 ms/frame — drops exact duplicates)
  → 3. diversity    (~1 ms/frame  — keeps the most diverse subset)
  → budget cap      → DFT labeling
```

Each stage is independent and reorderable.

---

## `physicality` — remove unphysical frames

```toml
[selection]
steps        = ["physicality"]
min_distance = 0.7    # (1)
```

1. **`min_distance`** — minimum allowed interatomic distance in Å. Any frame
   with *any* pair of atoms closer than this is dropped.

This catches:
- MD runs that went unstable (atoms flying through each other)
- MC moves that created overlapping atoms
- Geometry optimisations that didn't converge

**Typical values:**

| System | Suggested `min_distance` |
|---|---|
| Organic molecules | 0.7 Å |
| Metal surfaces | 1.0 Å |
| Dense crystals | 1.2–1.5 Å |

---

## `dedup` — remove exact duplicates

```toml
[selection]
steps = ["physicality", "dedup"]
```

Dedup compares frames by their **content hash** — a SHA-1 of atomic numbers,
positions (rounded to 4 decimal places), cell, and PBC flags. Two frames with
identical hashes are structurally identical up to floating-point precision.

This matters more than it sounds: if you restart a run, or combine outputs from
multiple sources, dedup ensures you never label the same geometry twice.

!!! tip "Dedup across runs"
    The `Dataset` class also deduplicates when you append to it. So even if you
    run the pipeline 10 times with different seeds, you won't accumulate duplicate
    frames in your final dataset.

---

## `diversity` — farthest-point sampling

```toml
[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 20    # (1)
```

1. **`budget`** — the maximum number of frames to keep. Diversity picks the
   `budget` most mutually-distant frames in descriptor space.

The descriptor is currently a **composition + geometry histogram**: atomic
species counts and pair-distance histograms. It's cheap to compute (no
expensive SOAP/ACE) and works well for diverse-structure selection.

```
       ┌──────────────────────────────────┐
       │         ALL FRAMES               │
       │  ●  ●  ●  ●  ●  ●  ●  ●  ●  ●  │
       │  ●  ●  ●  ●  ●  ●  ●  ●  ●  ●  │
       └──────────────────────────────────┘
                       ↓ diversity (budget=5)
       ┌──────────────────────────────────┐
       │  ★           ★            ★      │
       │        ★             ★           │
       └──────────────────────────────────┘
```

Farthest-point sampling (FPS) guarantees that no two selected frames are
closer in descriptor space than half the selection radius — maximising coverage.

---

## Configuring the full funnel

```toml
[selection]
steps        = ["physicality", "dedup", "diversity"]   # (1)
budget       = 50                                       # (2)
min_distance = 0.8                                      # (3)
```

1. **`steps`** — the ordered list. You can drop any stage you don't need.
   Stages run left-to-right; each stage's output is the next stage's input.

2. **`budget`** — applied at the very end, after all filter stages. If fewer
   frames survive the filters than the budget, all surviving frames are kept.

3. **`min_distance`** — used by the `physicality` stage.

---

## Example: aggressive diversity filtering

Start with many frames, keep only the most diverse 5%:

```toml
[sampling]
type        = "md"
temperature = 800.0   # high T → many diverse candidates
steps       = 2000
interval    = 10      # → 200 candidate frames

[selection]
steps        = ["physicality", "dedup", "diversity"]
budget       = 10     # keep only 10 of 200 = 5%
min_distance = 0.7
```

---

## Example: no diversity (just cleanup)

When you've already done diversity selection externally, just apply the cheap
filters:

```toml
[selection]
steps        = ["physicality", "dedup"]
budget       = null    # no budget cap — keep everything that passes
min_distance = 1.0
```

Setting `budget = null` (or simply omitting it) removes the budget cap.

---

## Using the funnel from Python

You can apply the funnel programmatically to any list of `Structure` objects:

```python
from traincraft import run_funnel, read_frames
from traincraft.config.models import SelectionConfig

frames = read_frames("all_candidates.extxyz")

cfg = SelectionConfig(
    steps=["physicality", "dedup", "diversity"],
    budget=20,
    min_distance=0.8,
)
selected = run_funnel(frames, cfg)
print(f"{len(selected)} frames selected from {len(frames)}")
```

---

## Summary

| Stage | What it removes | Cost |
|---|---|---|
| `physicality` | Unphysical / crashed geometries | Very fast (distance check) |
| `dedup` | Exact structural duplicates | Fast (hash lookup) |
| `diversity` | Redundant near-duplicates | Moderate (O(n²) FPS) |

**Next:** [Tutorial 9](09-interop.md) — converting between ASE, pymatgen, and RDKit.
