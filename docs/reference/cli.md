# CLI Reference

The `traincraft` CLI is installed as a console script when you install the package.

```bash
traincraft --help
```

---

## `traincraft run`

Run the complete workflow declared in a TOML config file.

```bash
traincraft run CONFIG [--force]
```

**Arguments / options:**

| Argument / option | Description |
|---|---|
| `CONFIG` | Path to a `.toml` file (must exist and be readable) |
| `--force`, `-f` | Recompute **every** stage, ignoring cached artifacts (a clean slate) |

**What it does:**

1. Loads and validates `CONFIG` via pydantic.
2. Creates the workspace directory (`runs/<name>/`).
3. Runs each active stage in order: geometry → sampling → selection → labeling →
   dataset → training (only the sections you declared run).
4. Prints a summary to stdout.

!!! note "Slurm routing"
    If `[orchestration].engine = "slurm"`, `run` instead renders and submits the
    stages as dependency-chained Slurm jobs (see [`submit`](#traincraft-submit)).

**Example:**

```bash
traincraft run examples/07_co_on_cu_mc.toml
```

```
2026-05-30 10:00:01 [INFO] MD: T=500.0 steps=100 interval=20
2026-05-30 10:00:02 [INFO] dataset written: runs/07_co_on_cu_mc/dataset.extxyz (3 frames)
Done:
  workspace:    runs/07_co_on_cu_mc
  n_candidates: 6
  n_selected:   3
  dataset:      runs/07_co_on_cu_mc/dataset.extxyz
```

!!! tip "Smart caching — rerunning only redoes what changed"
    Each stage writes a sidecar `<artifact>.key` fingerprinting the inputs that
    produced it: its own config section (plus the run seed) **and** the upstream
    stage's key. On a rerun a stage is skipped only if its artifact exists *and*
    that fingerprint still matches.

    So if you **edit the config and rerun**, the stages whose settings changed —
    and everything downstream of them — recompute automatically, while unrelated
    or expensive stages (e.g. DFT labeling) stay cached. Rerunning an *unchanged*
    config is a no-op that finishes in seconds.

    Pass `--force` to ignore the fingerprints and recompute everything from
    scratch. (Per stage, the same escape hatch is `traincraft stage NAME CONFIG
    --force`.)

---

## `traincraft stage`

Run a **single** pipeline stage standalone, reading the previous stage's artifact
from the workspace and writing its own. This is what the HPC executor dispatches
as separate jobs — and it's how you train on an existing dataset.

```bash
traincraft stage NAME CONFIG [--force]
```

**Arguments / options:**

| Argument | Description |
|---|---|
| `NAME` | One of: `geometry`, `sample`, `select`, `label`, `dataset`, `train` |
| `CONFIG` | Path to a `.toml` file |
| `--force` | Recompute even if the stage's artifact already exists |

**Example — train on a dataset you already built:**

```bash
pixi run -e mace traincraft stage train my_run.toml
```

```
stage 'train': 42 frames
```

The `train` stage writes `runs/<name>/model/<name>.model` and a `manifest.json`.
See [Tutorial 10 · Training](../tutorials/10-training.md).

---

## `traincraft submit`

Render (and, unless `--dry-run`, submit) the workflow as dependency-chained
Slurm jobs via Apptainer. Requires an `[orchestration.slurm]` section.

```bash
traincraft submit CONFIG [--dry-run]
```

| Option | Description |
|---|---|
| `--dry-run` | Write the sbatch scripts but don't submit |

The GPU stages (`sample`, `train`) run in `traincraft-mlip.sif` with `--nv`; the
DFT `label` stage injects the engine command. See [Run on HPC](../how-to/hpc.md).

---

## `traincraft validate`

Validate a config file without running anything.

```bash
traincraft validate CONFIG
```

**Arguments:**

| Argument | Description |
|---|---|
| `CONFIG` | Path to a `.toml` file |

**What it does:** loads and validates the config, prints which stages are active,
and exits. Use this to catch typos and logical errors before a long run.

**Example:**

```bash
traincraft validate examples/07_co_on_cu_mc.toml
```

```
OK: config is valid
  run:    07_co_on_cu_mc
  stages: geometry, calculator, sampling, selection, dataset
```

If the config is invalid:

```bash
traincraft validate broken.toml
```

```
Error: 1 validation error for TrainCraftConfig
geometry -> builder -> type
  Input should be 'nanotube', 'molecule', ... [type=literal_error, ...]
```

---

## `traincraft new`

Write a starter config file to a given path.

```bash
traincraft new PATH
```

**Arguments:**

| Argument | Description |
|---|---|
| `PATH` | Where to write the starter `.toml` file |

**What it does:** generates a commented starter config covering all major sections
and writes it to `PATH`. Refuses to overwrite an existing file.

**Example:**

```bash
traincraft new my_run.toml
cat my_run.toml
```

---

## `traincraft plugins`

List all registered plugins grouped by kind.

```bash
traincraft plugins
```

**Output:**

```
source:    file, materials_project, optimade, pubchem, scratch, smiles, url
builder:   crystal, intercalation, layered, liquid, molecule, nanotube, slab, surface_adsorbate, surface_packing
transform: constraints, perturb, rotate, set_pbc, strain, supercell, vacuum
calculator: emt, fhi_aims, mace, qe, tblite, xtb
sampler:   md, monte_carlo, rattle
selector:  dedup, diversity, physicality
trainer:   mace
```

!!! note "Plugin availability depends on your environment"
    `smiles` and `surface_packing` appear only when RDKit/Packmol are installed.
    `mace` and `xtb` appear only when their respective packages are installed.
