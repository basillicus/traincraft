---
hide:
  - toc
---

<div class="hero" markdown>

# TrainCraft

**Build training data for machine-learned interatomic potentials — systematically.**

TrainCraft is a modular Python toolkit for generating, selecting, and labeling
atomic structures to train MACE and other MLIPs. It covers the full pipeline:
geometry → sampling → selection → (DFT labeling) → training — with a clean
config file as the only glue.

[Get started in 5 minutes](getting-started/quickstart.md){ .md-button .md-button--primary }
[Browse the tutorials](tutorials/index.md){ .md-button }

</div>

---

## Why TrainCraft?

<div class="grid cards" markdown>

-   :material-molecule: **Every structure type covered**

    ---

    Bulk crystals, point defects, bare slabs, molecules on surfaces, 2D
    bilayers, moiré twists, SMILES-derived conformers — all from one TOML key.

-   :material-filter: **Principled dataset selection**

    ---

    A composable funnel — physicality → dedup → diversity (FPS) — removes
    unphysical frames and redundant near-duplicates before any expensive DFT.

-   :material-atom: **MACE-first, model-agnostic**

    ---

    Foundation models (`mace-mp0`, `mace-off23`) work out of the box. Swap in
    a local fine-tuned checkpoint with one config line.

-   :material-puzzle: **Plugin architecture**

    ---

    Every builder, calculator, sampler, and selector is a decorated function in
    the registry. Adding a new capability is one new file — no dispatcher to edit.

-   :material-history: **Provenance everywhere**

    ---

    Every frame records exactly how it was made: source, builder, transforms,
    calculator, seed. The `origin` tag keeps cheap and expensive data separable.

-   :material-cog: **Config is data**

    ---

    One TOML file drives the whole workflow. Validated by pydantic v2 —
    typos fail loudly. The same format a future workflow editor would emit.

</div>

---

## 30-second install

```bash
# Recommended: pixi (manages conda-forge + PyPI in one lockfile)
curl -fsSL https://pixi.sh/install.sh | sh
git clone https://github.com/basillicus/traincraft && cd traincraft
pixi install          # core dependencies
pixi install -e dev   # + pytest / ruff / mypy

# Alternative: pip / uv
pip install traincraft              # or: uv pip install traincraft
pip install "traincraft[geometry]"  # + rdkit, pymatgen, packmol
```

---

## 30-second example

Create `my_run.toml`:

```toml
[run]
name = "hello_traincraft"

[geometry.builder]
type   = "nanotube"
n      = 5
m      = 0
length = 1

[calculator]
type = "emt"

[sampling]
type        = "md"
temperature = 300.0
steps       = 200
interval    = 20

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 5

[dataset]
path = "dataset"
```

Run it:

```bash
traincraft run my_run.toml
```

```
Done:
  workspace: runs/hello_traincraft
  n_candidates: 11
  n_selected:   5
  dataset:      runs/hello_traincraft/dataset.extxyz
```

Five diverse, physically valid frames — ready to label with DFT. That's the whole loop.

---

## What's next?

=== "I want to learn"

    Start with [Tutorial 1: Your First Dataset](tutorials/01-first-dataset.md) — it walks
    through every section of the config in detail, explaining *why* each choice exists.

=== "I need a specific structure type"

    Jump to the tutorial for your system:

    - [Molecules & SMILES](tutorials/02-molecules.md) — organic molecules, RDKit conformers
    - [Molecules on Surfaces](tutorials/03-surfaces.md) — adsorbate coverage, MC sampling
    - [Crystals & Defects](tutorials/04-crystals-defects.md) — bulk, vacancies, substitutions
    - [Slabs & Strain](tutorials/05-slabs-strain.md) — surface models, mechanical deformation
    - [2D Materials](tutorials/06-2d-materials.md) — graphene, hBN, MoS₂, moiré stacks

=== "I want the full reference"

    See the [Config Schema](reference/config.md) for every TOML field, the
    [CLI reference](reference/cli.md) for all commands, or the
    [Python API](reference/api.md) for library usage.
