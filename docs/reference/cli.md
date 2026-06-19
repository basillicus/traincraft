# CLI Reference

The `traincraft` CLI is installed as a console script when you install the package.

```bash
traincraft --help
```

---

## `traincraft run`

Run the complete workflow declared in a TOML config file.

```bash
traincraft run CONFIG
```

**Arguments:**

| Argument | Description |
|---|---|
| `CONFIG` | Path to a `.toml` file (must exist and be readable) |

**What it does:**

1. Loads and validates `CONFIG` via pydantic.
2. Creates the workspace directory (`runs/<name>/`).
3. Runs each active stage in order: geometry → sampling → selection → dataset.
4. Prints a summary to stdout.

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

!!! tip "Rerunning is safe"
    Each run writes to a named workspace. Running the same config twice writes
    to the same workspace and **appends** to the dataset (deduplicating by hash).
    Nothing is ever overwritten.

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
source:    file, scratch, smiles, url
builder:   crystal, layered, molecule, nanotube, slab, surface_adsorbate, surface_packing
transform: perturb, rotate, set_pbc, strain, supercell, vacuum
calculator: emt, mace, tblite, xtb
sampler:   md, monte_carlo, rattle
selector:  dedup, diversity, physicality
```

!!! note "Plugin availability depends on your environment"
    `smiles` and `surface_packing` appear only when RDKit/Packmol are installed.
    `mace` and `xtb` appear only when their respective packages are installed.
