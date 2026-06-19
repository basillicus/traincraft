# Tutorial 2 · Molecules & SMILES

**What you'll learn:** how to build molecules using ASE's built-in database,
SMILES strings, and local files — and how to sample their configuration space.

**Prerequisites:** Tutorial 1. For SMILES: `pixi install -e science` (RDKit).

---

## Option A — ASE g2 molecule database

The simplest way to get a molecule: look it up by its G2 database name.

```toml title="examples/02_molecule_emt_rattle.toml"
[run]
name = "water_rattle"
seed = 1

[geometry.builder]
type   = "molecule"
name   = "H2O"       # (1)
vacuum = 8.0          # (2)

[calculator]
type = "emt"

[sampling]
type         = "rattle"   # (3)
method       = "mc"
n_structures = 20
std          = 0.08
min_distance = 1.0

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 10

[dataset]
path = "dataset"
```

1. **`name`** — any molecule in [ASE's G2 database](https://wiki.fysik.dtu.dk/ase/ase/build/build.html#ase.build.molecule).
   Examples: `"H2O"`, `"CO2"`, `"CH4"`, `"NH3"`, `"C2H6"`.

2. **`vacuum`** — padding added around the molecule in all directions.
   8 Å is typical for isolated molecules.

3. **`rattle`** — HiPhive-based random displacement. Faster than MD for
   generating a diverse set of perturbed geometries (no time integration).
   Requires `pixi install -e science`.

---

## Option B — SMILES strings (RDKit)

For anything not in the G2 database, use a SMILES string. TrainCraft uses
RDKit's ETKDG algorithm to embed 3D coordinates and MMFF to optimize them.

```toml
[geometry.source]
type        = "smiles"
smiles      = "CCO"      # (1)
n_conformers = 3          # (2)
optimize    = true        # (3)
vacuum      = 8.0
```

1. **`smiles`** — standard SMILES notation. Ethanol: `"CCO"`. Toluene:
   `"Cc1ccccc1"`. Caffeine: `"Cn1cnc2c1c(=O)n(c(=O)n2C)C"`.

2. **`n_conformers`** — how many independent conformers to embed. Currently
   the first one is used as the initial structure; future versions will expose
   all conformers.

3. **`optimize`** — run MMFF geometry optimization after embedding (recommended).

!!! tip "Canonical SMILES"
    TrainCraft automatically canonicalises your SMILES via RDKit and stores
    the canonical form in the provenance. This means `"OCC"` and `"CCO"`
    produce the same provenance entry.

### Fragment tagging for SMILES

When a structure is built from a SMILES source, every atom is tagged as
**fragment 0** (a single mobile fragment). This means the Monte Carlo sampler
can rotate and translate the whole molecule as a rigid body.

```python
from traincraft.config.models import SmilesSource, GeometryConfig
from traincraft.geometry import build_geometry
from traincraft.core.fragments import get_fragments

s = build_geometry(GeometryConfig(source=SmilesSource(smiles="CCO", vacuum=8.0)))
frags = get_fragments(s.atoms)
print(frags)   # [0 0 0 0 0 0 0 0 0]  (all atoms are fragment 0)
```

---

## Option C — reading from a file

```toml
[geometry.source]
type = "file"
path = "my_molecule.xyz"   # any ASE-readable format
```

Or download directly from a URL:

```toml
[geometry.source]
type   = "url"
url    = "https://raw.githubusercontent.com/example/repo/main/ethanol.xyz"
format = "xyz"    # optional; inferred from the URL suffix if omitted
```

---

## Pairing molecules with MD sampling

While `rattle` is fast, MD gives more physically realistic trajectories for
flexible molecules:

```toml
[sampling]
type        = "md"
temperature = 500.0   # high T explores conformational space
steps       = 500
interval    = 25
timestep    = 0.5     # fs — shorter for light atoms (H)
```

!!! warning "EMT and organic molecules"
    EMT is parametrised for metals. For organic molecules, use `tblite`
    (GFN2-xTB) or `mace-off23` (organic foundation model):

    ```toml
    [calculator]
    type  = "tblite"
    method = "GFN2-xTB"
    ```

---

## Example: ethanol conformers

This generates a diverse set of ethanol geometries covering the OH and CC
rotational degrees of freedom:

```toml
[run]
name = "ethanol_conformers"
seed = 42

[geometry.source]
type     = "smiles"
smiles   = "CCO"
optimize = true
vacuum   = 8.0

[calculator]
type   = "tblite"
method = "GFN2-xTB"

[sampling]
type        = "md"
temperature = 800.0
steps       = 2000
interval    = 50
timestep    = 0.5

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 20

[dataset]
path = "dataset"
```

```bash
pixi run -e science traincraft run ethanol_conformers.toml
```

---

## Summary

| Method | When to use | Extra deps |
|---|---|---|
| `molecule` builder (g2 name) | Common molecules: H₂O, CO₂, CH₄ | None |
| `smiles` source | Any organic molecule | RDKit (`science` env) |
| `file` / `url` source | Pre-optimised geometries | None |

**Next:** [Tutorial 3](03-surfaces.md) — placing molecules on crystalline surfaces
and exploring them with Monte Carlo.
