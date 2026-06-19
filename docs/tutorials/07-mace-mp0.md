# Tutorial 7 · MACE-MP0 Sampling

**What you'll learn:** how to use MACE foundation models — `mace-mp0` (universal,
for any periodic system) and `mace-off23` (organic molecules) — as the energy/force
engine driving MD or MC sampling.

**Prerequisites:** Tutorial 1. Requires `pixi install -e mace` (torch + mace-torch).

---

## Why use a foundation model?

EMT is fast but limited to metals. tblite/GFN-xTB is good for molecules but slow
for large periodic systems. MACE-MP0 covers **89 elements** with near-DFT accuracy
for a wide range of crystal and surface structures — making it ideal for generating
high-quality initial training data before any fine-tuning.

---

## Using MACE-MP0

```toml title="examples/06_mace_mp0_sampling.toml"
[run]
name = "cu_mace_mp0"
seed = 42

[geometry.source]
type = "scratch"
bulk = "Cu"           # ase.build.bulk("Cu")

[[geometry.transforms]]
type   = "supercell"
repeat = [2, 2, 2]

[calculator]
type          = "mace"
model         = "mace-mp0"    # (1)
device        = "cpu"         # (2)
default_dtype = "float32"     # (3)

[sampling]
type        = "md"
temperature = 500.0
steps       = 200
interval    = 20

[selection]
steps        = ["physicality", "dedup", "diversity"]
budget       = 5
min_distance = 1.5

[dataset]
path = "dataset"
```

```bash
pixi run -e mace example-06
```

1. **`model`** — which foundation model to use:

    | Value | Coverage | Best for |
    |---|---|---|
    | `mace-mp0` | 89 elements, periodic systems | Inorganic crystals, surfaces, alloys |
    | `mace-off23` | Organic molecules (CHNO...) | Drug-like molecules, soft matter |

2. **`device`** — `"cpu"` or `"cuda"`. MACE is 10–100× faster on a GPU.

3. **`default_dtype`** — `"float32"` (faster, default) or `"float64"` (more precise).

!!! note "First-run download"
    On first use, MACE downloads the model weights (~50–100 MB) and caches them
    locally. Subsequent runs use the cached version.

---

## Using a local fine-tuned model

After you've fine-tuned MACE on your DFT data, point TrainCraft to the checkpoint:

```toml
[calculator]
type       = "mace"
model      = "mace-mp0"     # (1)
model_path = "/path/to/my_finetuned.model"   # (2)
device     = "cuda"
```

1. **`model`** — used as a fallback label if `model_path` is set (for provenance).

2. **`model_path`** — absolute path to a `.model` checkpoint file produced by
   `mace_run_train`. When set, this overrides the foundation model download.

See [How-to: Use Your Own MACE Model](../how-to/own-mace-model.md) for more.

---

## MACE-MP0 with surfaces

MACE-MP0 works well for metallic slabs and adsorbate systems:

```toml
[run]
name = "co_cu111_mace"
seed = 0

[geometry.builder]
type          = "surface_adsorbate"
element       = "Cu"
facet         = "fcc111"
size          = [3, 3, 4]
vacuum        = 12.0
molecule_name = "CO"
site          = "ontop"
height        = 1.9

[calculator]
type          = "mace"
model         = "mace-mp0"
device        = "cpu"
default_dtype = "float32"

[sampling]
type          = "monte_carlo"
steps         = 200
temperature   = 500.0
interval      = 20
p_translate   = 0.6
p_rotate      = 0.4
p_conformer   = 0.0
max_translate = 0.3
max_rotate    = 20.0

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 8

[dataset]
path = "dataset"
```

---

## Performance tips

!!! tip "CPU parallelism"
    Set `OMP_NUM_THREADS` and `MKL_NUM_THREADS` before running to control
    how many CPU cores PyTorch uses. For a 4-core machine:
    ```bash
    OMP_NUM_THREADS=4 traincraft run my_run.toml
    ```

!!! tip "Batch MD vs single-frame MC"
    For pure exploration, MD is faster per-frame. For systems with many local
    minima (e.g., molecules on surfaces), MC with MACE gives better coverage.

!!! tip "dtype trade-off"
    `float32` is 2× faster than `float64` and accurate enough for dataset
    generation. Use `float64` only if you need sub-meV/atom energy accuracy.

---

## Summary

| Model | pip package | When to use |
|---|---|---|
| `mace-mp0` | `mace-torch` | Inorganic systems, surfaces, alloys |
| `mace-off23` | `mace-torch` | Organic molecules |
| Custom checkpoint | `mace-torch` | After fine-tuning on your DFT data |

**Next:** [Tutorial 8](08-selection-funnel.md) — understanding and customising the selection funnel.
