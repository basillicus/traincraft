# Installation

## Requirements

- Python **3.10 – 3.12** (3.13 support pending conda-forge coverage)
- An x86-64 Linux machine for the full science stack; macOS works for core deps

---

## Option A — pixi (recommended)

[pixi](https://pixi.sh) is a fast, cross-platform package manager that mixes
conda-forge packages (compiled tools like Packmol, tblite, RDKit) with PyPI
packages in one reproducible lockfile.

```bash
# Install pixi itself (one-liner)
curl -fsSL https://pixi.sh/install.sh | sh

# Clone the repo
git clone https://github.com/basillicus/traincraft && cd traincraft

# Core environment (ASE, numpy, pydantic, typer — runs examples 01–07)
pixi install

# Dev environment (+ pytest, ruff, mypy)
pixi install -e dev

# Science environment (+ RDKit, pymatgen, Packmol, hiphive, tblite, dscribe)
# Needed for examples 08–14, SMILES sources, slab builders, conformer MC
pixi install -e science

# MACE environment (+ torch + mace-torch — needed for MACE-MP0)
pixi install -e mace
```

!!! tip "Activating an environment"
    You don't need to activate the environment to run commands — just prefix
    with `pixi run -e <env>`. For an interactive shell, use `pixi shell -e science`.

---

## Option B — pip / uv

```bash
pip install traincraft                  # core only (zero optional deps)
pip install "traincraft[geometry]"      # + rdkit, pymatgen
pip install "traincraft[sampling]"      # + hiphive
pip install "traincraft[semiempirical]" # + tblite / GFN-xTB
pip install "traincraft[mace]"          # + mace-torch (needs torch)
pip install "traincraft[dev]"           # + pytest, ruff, mypy
```

With uv (faster):
```bash
uv pip install "traincraft[geometry,sampling,dev]"
```

!!! warning "Packmol with pip"
    Packmol is a Fortran binary available on conda-forge but not PyPI. If you
    need `surface_packing` or liquid-mixture builders, use the pixi `science`
    environment or install Packmol manually from
    [github.com/m3g/packmol](https://github.com/m3g/packmol).

---

## Verifying the installation

```bash
traincraft --help
```

You should see:

```
Usage: traincraft [OPTIONS] COMMAND [ARGS]...

  TrainCraft: modular MLIP dataset generation & active learning.

Commands:
  new       Write a starter config to PATH.
  plugins   List registered plugins by kind.
  run       Run the entire workflow declared in CONFIG (a single TOML).
  validate  Validate CONFIG and print the resolved stages.
```

Check what plugins are registered in your current environment:

```bash
traincraft plugins
```

```
source:    file, scratch, smiles, url
builder:   crystal, layered, molecule, nanotube, slab, surface_adsorbate, surface_packing
transform: perturb, rotate, set_pbc, strain, supercell, vacuum
calculator: emt, mace, tblite, xtb
sampler:   md, monte_carlo, rattle
selector:  dedup, diversity, physicality
```

!!! note "Plugin availability depends on installed extras"
    `smiles`, `surface_packing`, and conformer-based MC moves require
    the `science` environment. `mace` and `xtb` require their respective extras.

---

## Building the docs

The docs you're reading are built with [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/).
To build them locally:

```bash
pixi install -e docs          # one-time setup
pixi run -e docs docs-serve   # live-reload server at http://127.0.0.1:8000
pixi run -e docs docs-build   # build to site/
```
