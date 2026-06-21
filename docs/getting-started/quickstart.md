# Quick Start

This page gets you from a blank directory to a complete, running TrainCraft
pipeline in about five minutes. No heavy dependencies required — just a working
Python install and TrainCraft.

---

## Step 1 — Generate a starter config

```bash
traincraft new my_first_run.toml
```

Open `my_first_run.toml`. It looks like this:

```toml title="my_first_run.toml (generated)"
[run]
name   = "traincraft_run"
outdir = "runs"
seed   = 42

[geometry.builder]
type   = "nanotube"
n      = 8
m      = 0
length = 1
vacuum = 6.0

[calculator]
type = "emt"

[sampling]
type        = "md"
temperature = 500.0
steps       = 1000
interval    = 20

[selection]
steps        = ["physicality", "dedup", "diversity"]
budget       = 50
min_distance = 0.7

[dataset]
path = "dataset"
```

Each section maps to one stage of the pipeline. Let's understand each one briefly:

| Section | What it does |
|---|---|
| `[run]` | Names the run and sets the output directory and random seed |
| `[geometry]` | Specifies what structure to build |
| `[calculator]` | Which energy/force calculator to use for sampling |
| `[sampling]` | How to explore configuration space |
| `[selection]` | A funnel that keeps only the most useful frames |
| `[dataset]` | Where to write the final `extxyz` dataset |

---

## Step 2 — Validate the config

Before running anything expensive, validate your config:

```bash
traincraft validate my_first_run.toml
```

```
OK: config is valid
  run:    traincraft_run
  stages: geometry, calculator, sampling, selection, dataset
```

This catches typos in field names, wrong types, and logical errors
(like a `sampling` section without a `calculator`).

---

## Step 3 — Run it

```bash
traincraft run my_first_run.toml
```

You'll see logging as each stage completes:

```
2026-05-30 10:00:01 [INFO] MD: T=500.0 steps=1000 interval=20
2026-05-30 10:00:05 [INFO] MD finished: 51 frames sampled
2026-05-30 10:00:05 [INFO] dataset written: runs/traincraft_run/dataset.extxyz (50 frames)
Done:
  workspace:    runs/traincraft_run
  n_candidates: 51
  n_selected:   50
  dataset:      runs/traincraft_run/dataset.extxyz
```

---

## Step 4 — Inspect the output

The workspace directory has a predictable layout:

```
runs/
  traincraft_run/
    structures/
      initial.extxyz      ← the structure as-built (before sampling)
    candidates/
      candidates.extxyz   ← all MD frames (before selection)
    selected/
      selected.extxyz     ← frames that passed the funnel
    dataset.extxyz        ← final dataset (hash-deduped, appended across runs)
```

Read the dataset with ASE or TrainCraft directly:

```python
from traincraft import read_frames

frames = read_frames("runs/traincraft_run/dataset.extxyz")
print(f"{len(frames)} frames")
print(frames[0].provenance.source)    # "builder:nanotube:8-0-l1"
print(frames[0].provenance.transforms) # ["supercell:(1,1,1)", ...]
```

Or with ASE directly:

```python
from ase.io import read

atoms_list = read("runs/traincraft_run/dataset.extxyz", index=":")
print(atoms_list[0].info["tc_provenance"])
```

---

## What's next?

You've just run a complete geometry → sampling → selection → dataset pipeline.
The result is a set of diverse, physically valid frames ready to be labeled with DFT.

- **[Tutorial 1](../tutorials/01-first-dataset.md)** — a deep dive into exactly this
  example, explaining every config option in detail.
- **[Tutorial 3](../tutorials/03-surfaces.md)** — molecules on surfaces with Monte
  Carlo sampling (more scientifically interesting for most use cases).
- **[Tutorial 10](../tutorials/10-training.md)** — once you have a labeled dataset,
  *train a MACE model* on it: add a `[training]` section and `traincraft run` does
  the rest.
- **[Config Schema](../reference/config.md)** — the complete reference for every
  TOML field.

!!! tip "Using TrainCraft as a library"
    You don't have to use the CLI — every stage is a pure Python function:

    ```python
    import traincraft as tc

    cfg = tc.load_config("my_first_run.toml")
    structure = tc.build_geometry(cfg.geometry)
    calc = tc.make_calculator(cfg.calculator)
    # ... etc.
    ```

    See the [Python API reference](../reference/api.md) for details.
