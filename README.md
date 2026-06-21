# TrainCraft

Modular generation of training datasets for — and training of — machine-learned
interatomic potentials (MLIPs), with a MACE-centric **active-learning** loop that
only spends expensive DFT on structures that actually improve the model.

> **Status:** greenfield rewrite in progress (`rewrite/foundation-vertical-slice`,
> commit `d43d3b8`). Architecture: [`DESIGN.md`](DESIGN.md). Phased plan:
> [`ROADMAP.md`](ROADMAP.md). Legacy Raman-work code is preserved at git commit
> `e39647c` — this rewrite has **no backward compatibility** with it.

---

## 📚 Documentation

**📖 Hosted docs: <https://basillicus.github.io/traincraft/>** (built and deployed
from `main` by GitHub Actions).

The `docs/` folder is a **MkDocs (Material)** site — reading the raw `.md` files
directly misses the navigation, cross-links, and auto-generated API reference, so
read it as the rendered site above, or serve it locally:

```bash
pixi install -e docs      # one-time: MkDocs + Material + mkdocstrings
pixi run docs-serve       # live-reload server → open http://127.0.0.1:8000
```

To produce static HTML (e.g. to publish or browse offline):

```bash
pixi run docs-build       # strict build into ./site/  — open site/index.html
```

Good entry points once the site is up:

- **[Getting started → Full Workflow: From Scratch to HPC](docs/getting-started/end-to-end.md)**
  — the line-by-line guide from `git clone` to a DFT-labeled dataset on a cluster.
- **Concepts** — [the pipeline](docs/concepts/index.md), [calculators & DFT
  labeling](docs/concepts/calculators.md), [geometry](docs/concepts/geometry.md),
  [provenance](docs/concepts/provenance.md).
- **How-to → [Run on HPC (Slurm + Apptainer)](docs/how-to/hpc.md)** — containers,
  the `runtime`/`mpi` knobs, and worked Leonardo (PMIx) and LUMI (Cray) examples.

> The links above point at the source files for convenience; in the served site
> they appear in the left-hand navigation with full cross-linking and search.

---

## Quickstart

### 1 — Install

**With pixi (recommended — mixes conda-forge and PyPI in one lockfile):**

```bash
git clone https://github.com/basillicus/traincraft && cd traincraft

pixi install           # default env: core deps only (ASE, pydantic, typer)
pixi install -e dev    # + pytest, ruff, mypy
pixi install -e science  # + HiPhive rattle, tblite/GFN-xTB, geometry tools, dscribe
pixi install -e mace   # + torch + mace-torch (heavy; downloads model on first run)
```

**With pip/uv (pure-Python, no conda):**

```bash
uv pip install -e ".[dev]"               # core + dev tools
uv pip install -e ".[dev,sampling]"      # + HiPhive rattle
uv pip install -e ".[dev,semiempirical]" # + tblite/GFN-xTB
uv pip install -e ".[dev,mace]"          # + MACE + torch
```

### 2 — Run your first example

**Important: always run examples via pixi to ensure the right environment is used.**

```bash
# No heavy deps — uses ASE EMT force field (works in any env):
pixi run example-01                     # from the repo root

# Or equivalently:
pixi run -e dev traincraft run examples/01_cnt_emt_md.toml
```

> **Why `pixi run` and not bare `traincraft`?**  
> `pixi run` executes the command inside the pixi-managed environment. If you
> run `traincraft` directly from a shell that isn't inside a pixi env, the
> right packages (e.g. hiphive, tblite, mace-torch) may not be on the Python
> path even if `pixi install -e science` has been run.  
> Use `pixi shell -e science` to enter the environment interactively, then
> bare `traincraft` commands work as expected.

Examples that need optional dependencies:

```bash
# HiPhive rattle — install the science env first:
pixi install -e science
pixi run -e science example-02

# MACE-MP0 sampling — install the mace env first:
pixi install -e mace
pixi run -e mace example-06
```

This builds a (5,0) carbon nanotube, runs 50 steps of Langevin MD, filters
through the selection funnel (physicality → dedup → diversity), and writes 3
provenance-tagged frames to `runs/01_cnt_emt_md/dataset.extxyz`.

### 3 — Explore the output

```
runs/
└── 01_cnt_emt_md/
    ├── structures/initial.extxyz    ← generated geometry
    ├── candidates/candidates.extxyz ← all MD frames
    ├── candidates/md.traj           ← full ASE trajectory
    ├── selected/selected.extxyz     ← post-funnel, pre-label
    └── dataset.extxyz               ← the final dataset (with provenance)
```

Every frame in `dataset.extxyz` carries `tc_provenance` (origin, calculator,
parent hash, seed) and `tc_energy` / `tc_forces` where available.

### 4 — Validate or scaffold a config

```bash
traincraft validate examples/01_cnt_emt_md.toml  # check config, show stages
traincraft new my_run.toml                        # write a starter config
traincraft plugins                                # list all registered plugins
```

---

## Examples

| File | What it shows | Extra deps |
|------|---------------|------------|
| [`01_cnt_emt_md.toml`](examples/01_cnt_emt_md.toml) | Nanotube + EMT + MD | none |
| [`02_molecule_emt_rattle.toml`](examples/02_molecule_emt_rattle.toml) | Molecule + rattle (HiPhive) | `sampling` |
| [`03_molecule_from_file.toml`](examples/03_molecule_from_file.toml) | Load any ASE-readable file | none |
| [`04_nanotube_supercell_rattle.toml`](examples/04_nanotube_supercell_rattle.toml) | Supercell + rattle | `sampling` |
| [`05_selection_funnel_demo.toml`](examples/05_selection_funnel_demo.toml) | Selection funnel parameters | none |
| [`06_mace_mp0_sampling.toml`](examples/06_mace_mp0_sampling.toml) | MACE-MP0 as sampler | `mace` |
| [`07_co_on_cu_mc.toml`](examples/07_co_on_cu_mc.toml) | CO on Cu(111) + Metropolis MC | none |
| [`08_smiles_molecule.toml`](examples/08_smiles_molecule.toml) | Molecule from SMILES + MD | `science` |
| [`09_packing_on_surface.toml`](examples/09_packing_on_surface.toml) | Packmol surface coverage + MC | `science` |
| [`10_ethanol_conformers_on_cu.toml`](examples/10_ethanol_conformers_on_cu.toml) | Ethanol gauche/anti conformers on Cu(111) | `science` |
| [`11_butane_conformers_on_cu.toml`](examples/11_butane_conformers_on_cu.toml) | Butane conformer landscape, 3-molecule coverage | `science` |

---

## CLI reference

```bash
traincraft run      <config.toml>   # execute the full workflow
traincraft validate <config.toml>   # parse + validate only
traincraft new      <path.toml>     # write a starter config
traincraft plugins                  # list registered plugins by kind
```

---

## One TOML drives the whole workflow

A single validated TOML declares every stage. Presence of a section enables that
stage; its absence skips it. Adding a new geometry builder, calculator, sampler,
or selector is one self-registering file — no dispatcher to edit.

```toml
[run]
name   = "my_run"
outdir = "runs"
seed   = 42

[geometry.builder]
type   = "nanotube"   # or: molecule / crystal / surface / ...
n = 8; m = 0; length = 2

[calculator]
type   = "mace"       # or: emt / tblite / xtb
model  = "mace-mp0"

[sampling]
type        = "md"    # or: rattle / monte_carlo
temperature = 500.0
steps       = 1000

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 50

[dataset]
path = "dataset"
```

---

## Use as a library

Every stage is a pure function — use them directly in your own scripts:

```python
import traincraft as tc

# Run the full pipeline from a config file
cfg     = tc.load_config("examples/01_cnt_emt_md.toml")
summary = tc.run_pipeline(cfg)

# Or compose individual steps
structure = tc.build_geometry(cfg.geometry)
calc      = tc.make_calculator(cfg.calculator)
frames    = tc.run_sampling(structure, calc, job, cfg.sampling)
selected  = tc.run_funnel(frames, cfg.selection)
tc.write_frames("out.extxyz", selected)
```

---

## Package layout

```
src/traincraft/
  __init__.py            # curated public API
  cli.py                 # Typer shell (run / validate / new / plugins)
  config/                # pydantic v2 models + TOML loader
  core/                  # Structure, registry, Workspace, Result, provenance, rng
  geometry/
    sources/             # file, scratch  (SMILES/URL/MP — Phase 2)
    builders/            # nanotube, molecule  (crystal/surface/… — Phase 2)
    transforms/          # vacuum, supercell, perturb  (strain/rotate/… — Phase 2)
  calculators/
    potentials.py        # cheap calcs: emt (force field), tblite/xtb (semiempirical),
                         #              mace (MLIP)  — these are NOT all MLIPs
    dft.py               # QE + FHI-AIMS, E/F/stress/dipole/polarizability — Phase 3
  sampling/              # md, rattle, monte_carlo (rigid-body + RDKit conformers)
  selection/             # physicality, dedup, diversity (FPS), budget
                         #   uncertainty/committee — Phase 5
  datasets/              # extxyz IO, Dataset (dedup + filter by provenance)
                         #   health tooling (coverage, distributions) — Phase 4
  training/              # MACE fine-tune/train, multi-head — Phase 4
  validation/            # parity, learning curves, IR/Raman spectra — Phase 4
  active_learning/       # explore → select → label → train → converge — Phase 5
  orchestration/
    local.py             # serial engine (now)
                         # QuACC / Covalent / Parsl adapter — Phase 6
```

---

## Data organisation

Runs write a predictable directory tree. DFT-labeled frames are always in a
separate sub-directory so the expensive dataset is clean and shareable:

```
runs/<run-name>/
  structures/      # generated geometries (origin: generated)
  candidates/      # sampler output, all frames (origin: ml_sampled)
  selected/        # post-funnel frames, queued for labeling
  labeled_dft/     # DFT-labeled ← THE SHAREABLE DATASET  [Phase 3]
                   #   + manifest.json with level-of-theory and counts
  models/          # trained MACE checkpoints + metrics    [Phase 4]
  logs/
```

Every frame carries an `origin` tag in its provenance:
`generated` → `ml_sampled` → `ml_labeled` → `dft_labeled`.

---

## Environment (pixi)

Pixi mixes conda-forge (packmol, hiphive, tblite, QE) and PyPI (mace-torch,
dscribe, mdapackmol-fmt) in a single reproducible lockfile.

| Environment | Command | Contents |
|-------------|---------|----------|
| `default` | `pixi install` | core (ASE, pydantic, typer, tomlkit) |
| `dev` | `pixi install -e dev` | + pytest, ruff, mypy |
| `science` | `pixi install -e science` | + hiphive, tblite, geometry tools, dscribe |
| `mace` | `pixi install -e mace` | + torch, mace-torch |
| `docs` | `pixi install -e docs` | + mkdocs-material, mkdocstrings (see [Documentation](#-documentation)) |

---

## License

See [`LICENSE`](LICENSE). If you use the legacy TrainCraft code (commit
`e39647c`) in research, please cite Zenodo DOI `10.5281/zenodo.8174842` and
the packages it builds on.
