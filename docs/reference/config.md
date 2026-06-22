# Config Schema

Complete reference for every field in a TrainCraft TOML config file.
All fields are validated by pydantic v2. Unknown fields raise a `ValidationError`.

## Paths

Every filesystem path in the config (`outdir`, `dataset.path`, source/builder
`path`/`file`, `model_path`, `pt_train_file`, `slurm.sif_dir`, …) supports:

- **`~` expansion** — a leading `~` resolves to your home directory, so
  `outdir = "~/tests/runs"` lands in `$HOME/tests/runs`, not a literal `~`
  folder under the working directory.
- **Environment variables** — `$VAR` / `${VAR}` are expanded, e.g.
  `path = "$SCRATCH/in.xyz"`.

Relative *input* paths (e.g. `[geometry.source] path`) are resolved against the
config file's directory; relative *output* paths (`outdir`, `dataset.path`) are
left relative to the working directory. Expansion happens before this, so an
expanded `~`/`$VAR` path is absolute and is used as-is. Expansion uses the home
and environment of wherever the config is *parsed* — on a different machine,
prefer absolute paths or environment variables that exist there.

---

## `[run]`

Global settings for the pipeline run.

| Field | Type | Default | Description |
|---|---|---|---|
| `name` | `str` | `"traincraft_run"` | Name of the run; becomes the workspace subdirectory |
| `outdir` | `str` | `"runs"` | Root output directory (relative to working dir; `~`/`$VAR` expanded — see [Paths](#paths)) |
| `seed` | `int \| null` | `42` | RNG seed for reproducibility. `null` → random |

```toml
[run]
name   = "my_experiment"
outdir = "runs"
seed   = 42
```

---

## `[geometry]`

Specifies the initial structure. Exactly one of `source` or `builder` must be set.
Transforms are optional and applied in order.

```toml
[geometry]              # section header (optional if using sub-tables)

[geometry.source]       # OR [geometry.builder]
type = "..."

[[geometry.transforms]]
type = "..."
```

---

## Sources

### `type = "file"`

| Field | Type | Default | Description |
|---|---|---|---|
| `path` | `str` | required | Path to any ASE-readable file (extxyz, POSCAR, CIF, xyz, …) |

### `type = "scratch"`

| Field | Type | Default | Description |
|---|---|---|---|
| `molecule` | `str \| null` | `null` | ASE g2 molecule name (e.g. `"H2O"`, `"CO2"`) |
| `bulk` | `str \| null` | `null` | ASE bulk symbol (e.g. `"Cu"`, `"Si"`) |

Exactly one of `molecule` or `bulk` must be set.

### `type = "smiles"`

Requires RDKit (`pixi install -e science`).

| Field | Type | Default | Description |
|---|---|---|---|
| `smiles` | `str` | required | SMILES string (canonicalised automatically) |
| `n_conformers` | `int` | `1` | Number of ETKDG conformers to generate |
| `optimize` | `bool` | `true` | Run MMFF geometry optimisation |
| `seed` | `int \| null` | `null` | RNG seed for conformer embedding |
| `vacuum` | `float` | `6.0` | Vacuum padding in Å |

### `type = "url"`

| Field | Type | Default | Description |
|---|---|---|---|
| `url` | `str` | required | URL to download (http/https/file://) |
| `format` | `str \| null` | `null` | ASE format override (inferred from suffix if null) |
| `timeout` | `float` | `30.0` | Download timeout in seconds |

---

## Builders

### Shared blocks: `Species` and `AlloyComponent`

Two small blocks are reused across builders.

**`Species`** — one component of a *mixture*. Used by `liquid`,
`surface_packing` and `filled_nanotube` (as a `[[geometry.builder.species]]`
list). Provide exactly one *identity* and at most one *amount*:

| Field | Type | Default | Description |
|---|---|---|---|
| `molecule_name` | `str \| null` | `null` | ASE g2 name (identity) |
| `smiles` | `str \| null` | `null` | SMILES, built with RDKit (identity) |
| `file` | `str \| null` | `null` | Path to a structure file (identity) |
| `count` | `int \| null` | `null` | Absolute number of copies |
| `ratio` | `float \| null` | `null` | Relative share (apportioned from `n_molecules`) |

Within one mixture use a single amount style — all `count`, or all `ratio`
(then the builder's `n_molecules` is the total). A bare species is one copy
(count mode) or an equal share (ratio mode).

**`AlloyComponent`** — one substituent of a *random solid solution*. Used by
`crystal`, `slab`, `surface_adsorbate` and `surface_packing` (as a
`[[geometry.builder.composition]]` list):

| Field | Type | Default | Description |
|---|---|---|---|
| `element` | `str` | required | Element to substitute in |
| `ratio` | `float` | required | Fraction of host sites to replace, in (0, 1] |

Component ratios must sum to ≤ 1 (the remainder stays the host element);
substituted sites are picked with the builder's `seed`.

### `type = "nanotube"`

| Field | Type | Default | Description |
|---|---|---|---|
| `n` | `int` | `8` | Chiral index n |
| `m` | `int` | `0` | Chiral index m (`m=0` → zigzag; `n=m` → armchair) |
| `length` | `int` | `1` | Number of unit cells along the tube axis |
| `bond` | `float` | `1.42` | C–C bond length in Å |
| `vacuum` | `float` | `6.0` | Radial vacuum in Å |

### `type = "filled_nanotube"`

A carbon nanotube **randomly filled with molecules** (the original "fillMyTubes"
use case). Packmol packs `n_molecules` guests into a cylinder coaxial with the
tube and just narrower than its wall. The tube carbons are tagged framework
(`tc_fragment == -1`); each guest molecule is its own fragment. Requires the
Packmol binary (`pixi install -e science`).

| Field | Type | Default | Description |
|---|---|---|---|
| `n` | `int` | `10` | Chiral index n (default is roomy armchair) |
| `m` | `int` | `10` | Chiral index m |
| `length` | `int` | `6` | Unit cells along the tube axis |
| `bond` | `float` | `1.42` | C–C bond length in Å |
| `vacuum` | `float` | `8.0` | Radial vacuum in Å |
| `molecule_name` | `str \| null` | `null` | Single guest: ASE g2 name (e.g. `"H2O"`) |
| `smiles` | `str \| null` | `null` | Single guest: SMILES (needs RDKit) |
| `file` | `str \| null` | `null` | Single guest: path to a structure file |
| `species` | `list[Species]` | `[]` | A *mixture* of guests (alternative to the single-guest fields) |
| `n_molecules` | `int` | `4` | Single mode: copies; ratio mixture: total guests |
| `radial_margin` | `float` | `0.0` | Extra room taken off the packing cylinder (Å) |
| `axial_margin` | `float` | `1.5` | Inset at each tube end (Å; avoids PBC clashes) |
| `tolerance` | `float` | `2.0` | Packmol minimum separation (Å) |
| `pbc` | `bool` | `true` | Periodic along the tube axis |
| `seed` | `int \| null` | `null` | Packmol seed |

Use **either** a single guest (`molecule_name`/`smiles`/`file`) **or** a
`species` mixture — not both.

**No wall overlap.** The tube is handed to Packmol as a *fixed obstacle*, so
every guest atom is kept at least `tolerance` from the wall in true 3D — that is
what prevents overlap. The packing cylinder is only sized to keep guest atom
centres inside the carbon van der Waals shell, leaving the full cavity usable
(raise `tolerance` for more breathing room, or set `radial_margin` for extra).
If the tube is too narrow/short the builder raises a clear error — widen
(`n`/`m`), lengthen (`length`), or reduce the margins.

### `type = "molecule"`

| Field | Type | Default | Description |
|---|---|---|---|
| `name` | `str \| null` | `null` | ASE g2 molecule name |
| `smiles` | `str \| null` | `null` | SMILES string (Phase 2) |
| `vacuum` | `float` | `6.0` | Vacuum in Å |

### `type = "crystal"`

| Field | Type | Default | Description |
|---|---|---|---|
| `name` | `str` | required | Element symbol or compound (e.g. `"Cu"`, `"NaCl"`) |
| `crystalstructure` | `str \| null` | `null` | `fcc`, `bcc`, `hcp`, `diamond`, `rocksalt`, `wurtzite`, … |
| `a` | `float \| null` | `null` | Primary lattice constant in Å |
| `b` | `float \| null` | `null` | b lattice constant (orthorhombic) |
| `c` | `float \| null` | `null` | c lattice constant (HCP, tetragonal) |
| `cubic` | `bool` | `false` | Use the cubic conventional cell |
| `orthorhombic` | `bool` | `false` | Use the orthorhombic cell |
| `supercell` | `[int, int, int]` | `[1, 1, 1]` | Supercell repetitions |
| `defects` | `list[DefectSpec]` | `[]` | Point defects to introduce |
| `composition` | `list[AlloyComponent]` | `[]` | Random solid solution (mixed solid) |
| `seed` | `int \| null` | `null` | Seeds the random site substitution |

**`DefectSpec`:**

| Field | Type | Default | Description |
|---|---|---|---|
| `kind` | `"vacancy" \| "substitution" \| "interstitial"` | required | Defect type |
| `index` | `int \| null` | `null` (→ 0) | Atom index in the pre-defect supercell |
| `element` | `str \| null` | `null` | New element (substitution/interstitial) |
| `position` | `[float, float, float] \| null` | `null` | Insertion position (interstitial) |
| `cartesian` | `bool` | `false` | If true, `position` is Cartesian (Å); else fractional |

### `type = "slab"`

| Field | Type | Default | Description |
|---|---|---|---|
| `element` | `str` | required | Substrate element |
| `facet` | see below | `"fcc111"` | Named facet (ignored if `miller` is set) |
| `miller` | `[int, int, int] \| null` | `null` | Miller indices (overrides `facet`) |
| `crystalstructure` | `str \| null` | `null` | Needed for Miller mode |
| `a` | `float \| null` | `null` | Lattice constant for Miller mode |
| `c` | `float \| null` | `null` | c lattice constant for HCP |
| `cubic` | `bool` | `false` | Cubic cell for Miller mode |
| `layers` | `int` | `4` | Slab thickness in atomic layers (Miller mode) |
| `size` | `[int, int, int]` | `[3, 3, 4]` | `[nx, ny, nz]`; nz=layers in facet mode |
| `vacuum` | `float` | `12.0` | Total vacuum in Å |
| `orthogonal` | `bool` | `false` | Request an orthogonal cell (facet mode) |
| `periodic` | `bool` | `false` | Keep periodic along the surface normal (Miller mode) |
| `composition` | `list[AlloyComponent]` | `[]` | Random solid solution (mixed-solid slab) |
| `seed` | `int \| null` | `null` | Seeds the random site substitution |

**`facet` values:** `fcc111`, `fcc100`, `fcc110`, `bcc110`, `bcc100`, `bcc111`, `hcp0001`

### `type = "layered"`

| Field | Type | Default | Description |
|---|---|---|---|
| `material` | `"graphene" \| "hbn" \| "mx2"` | `"graphene"` | 2D crystal type |
| `formula` | `str \| null` | `null` | MX₂ formula (e.g. `"MoS2"`, `"WSe2"`) |
| `a` | `float \| null` | `null` | In-plane lattice constant override |
| `size` | `[int, int]` | `[1, 1]` | In-plane repetitions of the primitive cell |
| `n_layers` | `int` | `2` | Number of layers to stack |
| `interlayer_spacing` | `float` | `3.35` | Vertical layer spacing in Å |
| `stacking` | `"AA" \| "AB"` | `"AB"` | Stacking order |
| `twist` | `float` | `0.0` | Twist angle in degrees (non-zero → non-periodic flake) |
| `vacuum` | `float` | `15.0` | Total vacuum in Å |

### `type = "surface_adsorbate"`

| Field | Type | Default | Description |
|---|---|---|---|
| `element` | `str` | required | Substrate element |
| `facet` | `str` | `"fcc111"` | Named facet |
| `size` | `[int, int, int]` | `[3, 3, 4]` | Slab size |
| `vacuum` | `float` | `12.0` | Total vacuum in Å |
| `composition` | `list[AlloyComponent]` | `[]` | Random solid solution (mixed-solid slab) |
| `seed` | `int \| null` | `null` | Seeds the random site substitution |
| `molecule_name` | `str \| null` | `null` | ASE g2 adsorbate name |
| `smiles` | `str \| null` | `null` | SMILES adsorbate |
| `file` | `str \| null` | `null` | Path to adsorbate file |
| `site` | `str` | `"ontop"` | Initial site: `ontop`, `bridge`, `hollow`, `fcc`, `hcp` |
| `height` | `float` | `2.0` | Initial height above slab top (Å) |
| `offset` | `[float, float] \| null` | `null` | In-plane offset from the site |

Exactly one of `molecule_name`, `smiles`, or `file` must be set.

### `type = "surface_packing"`

All `surface_adsorbate` substrate parameters (including `composition` for a
mixed-solid slab), plus:

| Field | Type | Default | Description |
|---|---|---|---|
| `species` | `list[Species]` | `[]` | A *mixture* of adsorbates (alternative to a single molecule) |
| `n_molecules` | `int` | `4` | Single mode: copies; ratio mixture: total adsorbates |
| `tolerance` | `float` | `2.0` | Packmol intermolecular distance (Å) |
| `region_height` | `float` | `8.0` | Height of the packing box above the slab (Å) |
| `gap` | `float` | `2.0` | Clearance between slab top and packing box (Å) |
| `seed` | `int \| null` | `null` | Packmol seed / substitution seed |

Use **either** a single adsorbate (`molecule_name`/`smiles`/`file`) **or** a
`species` mixture — not both.

### `type = "liquid"`

Packs a `species` mixture into a periodic box (a liquid, a solvent blend, or
confined bulk). Requires the Packmol binary (`pixi install -e science`).

| Field | Type | Default | Description |
|---|---|---|---|
| `species` | `list[Species]` | required | The mixture to pack (≥ 1 species) |
| `n_molecules` | `int \| null` | `null` | Total for ratio-based mixtures |
| `box` | `[float, float, float] \| null` | `null` | Explicit cell edges in Å |
| `density` | `float \| null` | `null` | g/cm³; a cubic box is sized when `box` is unset |
| `tolerance` | `float` | `2.0` | Packmol min separation / wall inset (Å) |
| `pbc` | `bool` | `true` | Periodic cell |
| `seed` | `int \| null` | `null` | Packmol seed |

Set exactly one of `box` or `density`. Ratio-based mixtures need `n_molecules`.

---

## Transforms

### `type = "supercell"`

| Field | Type | Default | Description |
|---|---|---|---|
| `repeat` | `[int, int, int]` | `[1, 1, 1]` | Repetitions along each cell vector |

### `type = "vacuum"`

| Field | Type | Default | Description |
|---|---|---|---|
| `amount` | `float` | `6.0` | Total vacuum to add (Å); split equally on both sides |

### `type = "perturb"`

| Field | Type | Default | Description |
|---|---|---|---|
| `stddev` | `float` | `0.05` | Standard deviation of Gaussian displacements (Å) |

### `type = "strain"`

| Field | Type | Default | Description |
|---|---|---|---|
| `hydrostatic` | `float \| null` | `null` | Isotropic engineering strain (e.g. `0.02` = +2%) |
| `voigt` | `[float×6] \| null` | `null` | `(e_xx, e_yy, e_zz, e_yz, e_xz, e_xy)` |

Exactly one of `hydrostatic` or `voigt` must be set.

### `type = "rotate"`

| Field | Type | Default | Description |
|---|---|---|---|
| `angle` | `float` | required | Rotation angle in degrees |
| `axis` | `"x" \| "y" \| "z" \| [float×3]` | `"z"` | Rotation axis |
| `rotate_cell` | `bool` | `false` | Also rotate the periodic cell vectors |

### `type = "set_pbc"`

| Field | Type | Default | Description |
|---|---|---|---|
| `pbc` | `bool \| [bool×3]` | `true` | New PBC flags |

---

## `[calculator]`

### `type = "emt"`

No parameters. ASE's built-in force field — works for metals, zero deps.

### `type = "tblite"`

| Field | Type | Default | Description |
|---|---|---|---|
| `method` | `str` | `"GFN2-xTB"` | `"GFN2-xTB"` or `"GFN1-xTB"` |

### `type = "xtb"`

Same as `tblite`. Uses the `xtb-python` binding instead.

### `type = "mace"`

| Field | Type | Default | Description |
|---|---|---|---|
| `model` | `str` | `"mace-mp0"` | Foundation model: `"mace-mp0"` or `"mace-off23"` |
| `model_path` | `str \| null` | `null` | Path to a local `.model` checkpoint (overrides `model`) |
| `device` | `str` | `"cpu"` | `"cpu"` or `"cuda"` |
| `default_dtype` | `str` | `"float32"` | `"float32"` or `"float64"` |

---

## `[sampling]`

### `type = "md"`

Langevin NVT molecular dynamics.

| Field | Type | Default | Description |
|---|---|---|---|
| `temperature` | `float` | `500.0` | Temperature in K |
| `steps` | `int` | `1000` | Total MD steps |
| `interval` | `int` | `20` | Save a frame every N steps |
| `timestep` | `float` | `1.0` | Timestep in femtoseconds |
| `friction` | `float` | `0.02` | Langevin friction coefficient (1/fs) |

### `type = "rattle"`

HiPhive random displacement sampling. Requires `pixi install -e science`.

| Field | Type | Default | Description |
|---|---|---|---|
| `method` | `"mc" \| "standard"` | `"mc"` | Rattle method |
| `n_structures` | `int` | `10` | Number of rattled structures to generate |
| `std` | `float` | `0.12` | Displacement standard deviation (Å) |
| `min_distance` | `float` | `1.3` | Minimum interatomic distance (Å) |

### `type = "monte_carlo"`

Metropolis MC with rigid-body translate/rotate and optional conformer-swap moves.

| Field | Type | Default | Description |
|---|---|---|---|
| `steps` | `int` | `1000` | Total MC moves |
| `temperature` | `float` | `500.0` | Metropolis temperature (K) |
| `interval` | `int` | `20` | Save a snapshot every N accepted moves |
| `p_translate` | `float` | `0.5` | Probability of a translate move |
| `p_rotate` | `float` | `0.4` | Probability of a rotate move |
| `p_conformer` | `float` | `0.1` | Probability of a conformer-swap move |
| `max_translate` | `float` | `0.5` | Maximum displacement per translate move (Å) |
| `max_rotate` | `float` | `30.0` | Maximum rotation per rotate move (degrees) |
| `refresh_fragments` | `bool` | `false` | Re-infer fragments after each accepted move |
| `refresh_scale` | `float` | `1.2` | Bond-radius scale for `infer_fragments` |
| `seed` | `int \| null` | `null` | RNG seed |

### `type = "scan"`

**Deterministic** grid scan of one fragment's position/orientation — the exact,
combinatorial counterpart to the stochastic samplers above. It enumerates the
**Cartesian product** of a `translate` and/or `rotate` grid applied to the chosen
fragment, so it's the tool for potential-energy-surface scans, binding/
dissociation curves and orientation grids. Energies are attached by the
[`[labeling]`](#labeling) stage, like every sampler.

| Field | Type | Default | Description |
|---|---|---|---|
| `fragment` | `int \| "all"` | `0` | `tc_fragment` id to move (adsorbate/guest = `0`); `"all"` moves the whole system |
| `translate` | `ScanAxis \| null` | `null` | Grid of displacements (Å) along an axis, added to the fragment's position |
| `rotate` | `ScanAxis \| null` | `null` | Grid of rotations (degrees) about an axis through the fragment centroid |

At least one of `translate` / `rotate` is required. A **`ScanAxis`** is a small
table — `steps` values evenly spaced over `[start, stop]` (both ends inclusive):

| Field | Type | Default | Description |
|---|---|---|---|
| `axis` | `"x" \| "y" \| "z"` or `[x, y, z]` | `"z"` | Axis to translate along / rotate about |
| `start` | `float` | — | First grid value (Å for translate, degrees for rotate) |
| `stop` | `float` | — | Last grid value |
| `steps` | `int` | `5` | Number of grid points (so `translate.steps × rotate.steps` candidates) |

```toml
[sampling]
type = "scan"
fragment = 0
translate = { axis = "z", start = 0.0, stop = 3.0, steps = 13 }   # 13-point binding curve
rotate    = { axis = "x", start = 0.0, stop = 180.0, steps = 7 }  # × 7 orientations = 91 frames
```

!!! tip "Keep every grid point"
    A scan is meant to be exhaustive, so pair it with a generous
    [`[selection]`](#selection) `budget` (or drop the `diversity` step), or the
    funnel will thin your grid.

---

## `[selection]`

| Field | Type | Default | Description |
|---|---|---|---|
| `steps` | `list[str]` | `["physicality", "dedup", "diversity"]` | Ordered filter stages |
| `budget` | `int \| null` | `50` | Maximum frames to keep. `null` = no cap |
| `min_distance` | `float` | `0.7` | Physicality threshold (Å) |

---

## `[dataset]`

| Field | Type | Default | Description |
|---|---|---|---|
| `path` | `str` | `"dataset"` | Path (without extension) for the output extxyz dataset |
| `format` | `"extxyz"` | `"extxyz"` | Output format (only extxyz is currently supported) |

---

## `[labeling]`

Labels the **selected** frames with an (expensive) calculator — distinct from the
cheap `[calculator]` that drives sampling. Presence of this section enables the
`label` stage.

| Field | Type | Default | Description |
|---|---|---|---|
| `calculator` | calculator config | required | A `[labeling.calculator]` sub-table; usually `fhi_aims` or `qe` (any calculator works — `emt` stands in for DFT in demos) |

```toml
[labeling.calculator]
type             = "fhi_aims"
xc               = "pbe"
species_defaults = "tight"
properties       = ["polarizability"]   # E/F/stress always; + dipole/polarizability on request
```

See [Calculators & DFT Labeling](../concepts/calculators.md) for the DFT
calculator fields and the environment-injected run command.

---

## `[training]`

Trains a MACE model on the dataset (or the labelled frames). Presence of this
section enables the `train` stage, which runs after `dataset`. Defaults follow
[Tompa et al. (arXiv:2606.12704)](https://arxiv.org/abs/2606.12704); see the
[Training concept page](../concepts/training.md).

### `type = "mace"`

| Field | Type | Default | Description |
|---|---|---|---|
| `name` | `str` | `"mace_model"` | Output model name → `model/<name>.model` |
| `foundation_model` | `str \| null` | `"medium"` | Model to fine-tune from: `small`/`medium`/`large`, `mace-mp0`, `mace-off23`, or a path. `null` with `strategy="scratch"` trains from scratch |
| `strategy` | `"multihead" \| "naive" \| "scratch"` | `"multihead"` | Fine-tuning recipe. `multihead` = replay against forgetting |
| `heads` | `list[str]` | `["energy", "forces"]` | Properties to learn: any of `energy`, `forces`, `stress`, `dipole`, `polarizability`. Selects the MACE model type + loss |
| `e0s` | `str` | `"foundation"` | Isolated-atom energies: `"foundation"` (reuse), `"average"` (avoid), or a JSON `Z→energy` dict |
| `valid_fraction` | `float` | `0.1` | Held-out validation fraction (deterministic given `[run].seed`) |
| `pt_train_file` | `str \| null` | `null` | Replay data for multihead fine-tuning: `"mp"` or a path |
| `num_samples_pt` | `int` | `30000` | Number of replay samples |
| `weight_pt`, `weight_ft` | `float` | `1.0`, `1.0` | Replay-head / fine-tune-head weights |
| `energy_weight`, `forces_weight` | `float` | `10.0`, `10.0` | Loss weights (paper: energy-prioritised, constant) |
| `stress_weight`, `dipole_weight`, `polarizability_weight` | `float` | `1.0` | Per-head loss weights (used when that head is present) |
| `lr` | `float` | `1e-3` | Learning rate |
| `weight_decay` | `float` | `0.0` | Weight decay (paper: keep `0` for fine-tuning) |
| `ema` | `bool` | `true` | Use exponential moving average of weights |
| `ema_decay` | `float` | `0.995` | EMA decay (paper: `> 0.99` for fine-tuning) |
| `max_num_epochs` | `int` | `200` | Training epochs |
| `batch_size` | `int` | `10` | Batch size |
| `swa` | `bool` | `false` | Stochastic weight averaging (final-stage) |
| `hidden_irreps` | `str \| null` | `null` | Model size for `scratch` (e.g. `"128x0e + 128x1o"`) |
| `r_max` | `float` | `5.0` | Cutoff radius (Å) for `scratch` |
| `device` | `str` | `"cpu"` | `"cpu"` or `"cuda"` |
| `default_dtype` | `"float32" \| "float64"` | `"float64"` | Training precision |
| `seed` | `int \| null` | `null` | Trainer seed (also used for the split) |
| `extra` | `dict` | `{}` | Arbitrary `mace_run_train` flags (without `--`); **wins over** the defaults above |

```toml
[training]
type             = "mace"
name             = "cu_finetune"
foundation_model = "medium"
strategy         = "multihead"
heads            = ["energy", "forces", "stress"]
e0s              = "foundation"
pt_train_file    = "mp"
device           = "cuda"

# escape hatch: pass any raw mace_run_train flag (overrides a default)
[training.extra]
scheduler = "ReduceLROnPlateau"
```

!!! note "Multi-head model selection"
    `heads` chooses the MACE architecture automatically: `dipole` →
    `EnergyDipolesMACE`/`AtomicDipolesMACE`, `polarizability` →
    `AtomicDielectricMACE`. Override `--model`/`--loss`/keys via `[training.extra]`
    if your MACE version differs.
