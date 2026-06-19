# Tutorial 1 · Your First Dataset

**What you'll learn:** how every section of a TrainCraft TOML config works,
and what happens at each stage of the pipeline.

**Prerequisites:** TrainCraft installed (core only — no heavy deps).

**Time:** ~10 minutes.

---

## The complete config

Create `examples/01_cnt_emt_md.toml` (or use the one already in the repo):

```toml title="examples/01_cnt_emt_md.toml"
[run]
name   = "01_cnt_emt_md"
outdir = "runs"
seed   = 42

[geometry.builder]
type   = "nanotube"
n      = 5
m      = 0
length = 1
vacuum = 6.0

[calculator]
type = "emt"

[sampling]
type        = "md"
temperature = 300.0
steps       = 50
interval    = 10
timestep    = 1.0

[selection]
steps        = ["physicality", "dedup", "diversity"]
budget       = 3
min_distance = 0.7

[dataset]
path = "dataset"
```

Run it:

```bash
traincraft run examples/01_cnt_emt_md.toml
```

Now let's understand every line.

---

## `[run]` — bookkeeping

```toml
[run]
name   = "01_cnt_emt_md"   # (1)
outdir = "runs"             # (2)
seed   = 42                 # (3)
```

1. **`name`** — the subdirectory created inside `outdir`. All outputs land in
   `runs/01_cnt_emt_md/`. Changing the name starts a fresh workspace; the old
   one is untouched, making reruns safe.

2. **`outdir`** — where workspaces live. Relative to the working directory when
   you call `traincraft run`.

3. **`seed`** — seeds NumPy's global RNG so every run is reproducible. Set to
   `null` to use a random seed.

---

## `[geometry]` — building the structure

```toml
[geometry.builder]
type   = "nanotube"  # (1)
n      = 5           # (2)
m      = 0           # (3)
length = 1           # (4)
vacuum = 6.0         # (5)
```

1. **`type`** — selects a registered builder. Available: `nanotube`, `molecule`,
   `crystal`, `slab`, `layered`, `surface_adsorbate`, `surface_packing`.

2. **`n`, `m`** — the chiral indices of the carbon nanotube. `(5,0)` is a
   zigzag CNT; `(5,5)` would be armchair.

3. **`length`** — number of unit cells along the tube axis.

4. **`vacuum`** — adds `vacuum/2` Å of vacuum on each side in the radial
   directions. This makes the cell non-periodic in xy and periodic in z (the
   tube axis).

Instead of a `builder`, you can use a `source` to load an existing file:

```toml
[geometry.source]
type = "file"
path = "my_structure.extxyz"
```

Or download from a URL:

```toml
[geometry.source]
type    = "url"
url     = "https://raw.githubusercontent.com/.../co.xyz"
format  = "xyz"
```

### Transforms

After the geometry is built, you can chain transforms. For example, to triple
the unit cell along the tube axis:

```toml
[[geometry.transforms]]
type   = "supercell"
repeat = [1, 1, 3]

[[geometry.transforms]]
type   = "perturb"
stddev = 0.05
```

Transforms are applied in order. Available transforms: `supercell`, `vacuum`,
`perturb`, `strain`, `rotate`, `set_pbc`.

---

## `[calculator]` — the energy/force engine

```toml
[calculator]
type = "emt"
```

EMT (Effective Medium Theory) is a simple force field built into ASE. It
requires no extra dependencies and works for metals — perfect for tests and
demonstrations.

For real science, switch to:

=== "tblite / GFN-xTB (semiempirical)"

    ```toml
    [calculator]
    type   = "tblite"
    method = "GFN2-xTB"  # or "GFN1-xTB"
    ```
    Great for organic molecules. Requires `pixi install -e science`.

=== "MACE-MP0 (foundation MLIP)"

    ```toml
    [calculator]
    type         = "mace"
    model        = "mace-mp0"
    device       = "cpu"
    default_dtype = "float32"
    ```
    Universal potential covering elements 1–89. Requires `pixi install -e mace`.
    See [Tutorial 7](07-mace-mp0.md).

---

## `[sampling]` — exploring configuration space

```toml
[sampling]
type        = "md"      # (1)
temperature = 300.0     # (2)
steps       = 50        # (3)
interval    = 10        # (4)
timestep    = 1.0       # (5)
```

1. **`type = "md"`** — Langevin NVT molecular dynamics via ASE. Other options:
   `rattle` (random displacements, good for periodic solids) and `monte_carlo`
   (rigid-body moves + conformer swaps, ideal for molecules on surfaces).

2. **`temperature`** — in Kelvin. Higher temperature explores further from
   equilibrium and generates more diverse frames.

3. **`steps`** — total MD steps. 50 steps × 1 fs timestep = 50 fs of dynamics.

4. **`interval`** — save a frame every `interval` steps.
   Here: 50 / 10 = 5 frames (plus the initial frame = 6 total).

5. **`timestep`** — in femtoseconds. 1–2 fs is typical.

!!! note "How many frames to generate?"
    Generate roughly 5–10× more candidates than your final budget. The selection
    funnel will pick the most diverse, physically valid subset.

---

## `[selection]` — the quality filter

This is where TrainCraft earns its name. Before any expensive DFT labeling,
the funnel removes bad and redundant frames:

```toml
[selection]
steps        = ["physicality", "dedup", "diversity"]  # (1)
budget       = 3                                       # (2)
min_distance = 0.7                                     # (3)
```

1. **`steps`** — the ordered list of filter stages:
   - `physicality`: drops frames where any two atoms are closer than
     `min_distance` Å (catches MD crashes, overlapping atoms).
   - `dedup`: removes exact duplicates by content hash.
   - `diversity`: farthest-point sampling over a histogram descriptor — keeps
     the most structurally diverse subset.

2. **`budget`** — the maximum number of frames to keep after all filters.
   Here we keep 3 out of 6 candidates.

3. **`min_distance`** — the physicality threshold in Å. 0.7 Å is conservative
   (almost any pair of atoms closer than this is unphysical).

!!! tip "Reorder the funnel for speed"
    Run `physicality` and `dedup` first — they're cheap. Only then run
    `diversity` (FPS), which is O(n²) in the number of surviving frames.

---

## `[dataset]` — persistent, hash-deduped storage

```toml
[dataset]
path = "dataset"
```

TrainCraft writes `runs/01_cnt_emt_md/dataset.extxyz`. If you run the pipeline
multiple times (e.g., after updating the geometry), frames are **appended** and
deduplicated by content hash — you never get exact duplicates across runs.

---

## Inspecting the results

```python
from traincraft import read_frames

frames = read_frames("runs/01_cnt_emt_md/dataset.extxyz")
print(f"Got {len(frames)} frames")

for f in frames:
    print(f.provenance.source, f.provenance.transforms)
    # → "builder:nanotube:5-0-l1"  []
```

The `extxyz` file is also readable by any ASE-compatible tool (OVITO, ASE, etc.):

```bash
ase gui runs/01_cnt_emt_md/dataset.extxyz
```

---

## Summary

| Stage | What happened |
|---|---|
| **geometry** | Built a (5,0) CNT with 6 Å vacuum in x/y |
| **sampling** | 50 steps of Langevin MD at 300 K → 6 candidate frames |
| **selection** | physicality → dedup → diversity kept 3 frames |
| **dataset** | Written to `dataset.extxyz` with provenance |

**Next:** [Tutorial 2](02-molecules.md) — building molecules from chemical names and SMILES strings.
