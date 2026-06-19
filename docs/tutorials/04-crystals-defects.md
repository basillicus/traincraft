# Tutorial 4 · Crystals & Defects

**What you'll learn:** how to build bulk crystals with arbitrary structures and
lattice constants, create supercells, and introduce point defects (vacancies,
substitutions, interstitials).

**Prerequisites:** Tutorial 1. No extra dependencies.

---

## Building a bulk crystal

```toml title="examples/12_bulk_vacancy_md.toml"
[run]
name = "bulk_vacancy_md"
seed = 7

[geometry.builder]
type             = "crystal"
name             = "Cu"            # (1)
crystalstructure = "fcc"           # (2)
a                = 3.61            # (3)
cubic            = true            # (4)
supercell        = [2, 2, 2]       # (5)

[[geometry.builder.defects]]       # (6)
kind  = "vacancy"
index = 0

[calculator]
type = "emt"

[sampling]
type        = "md"
temperature = 600.0
steps       = 200
interval    = 20

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 8

[dataset]
path = "dataset"
```

```bash
pixi run example-12
```

### Crystal builder parameters

1. **`name`** — the element symbol (or compound formula for multi-element
   crystals, e.g. `"NaCl"`, `"GaAs"`, `"TiO2"`).

2. **`crystalstructure`** — the crystal structure type. Common values:
   `fcc`, `bcc`, `hcp`, `diamond`, `rocksalt`, `wurtzite`, `zincblende`.
   Leave `null` to let ASE infer it from the element.

3. **`a`** — the primary lattice constant in Å. Leave `null` to use ASE's
   built-in experimental value.

4. **`cubic`** — whether to use the cubic conventional cell (larger but
   orthogonal). Recommended for most MD runs.

5. **`supercell`** — repeat the primitive/conventional cell before adding
   defects. `[2,2,2]` → 8× the unit cell.

6. **`[[geometry.builder.defects]]`** — a list of point defects to introduce
   (see below).

---

## Point defects

### Vacancy

Remove one atom from the supercell:

```toml
[[geometry.builder.defects]]
kind  = "vacancy"
index = 0     # (1)
```

1. **`index`** — the atom index in the *pre-defect supercell* (0-based).
   Index 0 is the first atom.

### Substitution

Replace an atom with a different element:

```toml
[[geometry.builder.defects]]
kind    = "substitution"
index   = 0
element = "Ni"    # replace Cu with Ni
```

### Interstitial

Insert an extra atom at a specific position:

```toml
[[geometry.builder.defects]]
kind       = "interstitial"
element    = "H"
position   = [0.5, 0.5, 0.5]    # (1)
cartesian  = false               # (2)
```

1. **`position`** — where to insert the atom.

2. **`cartesian`** — `false` means the position is in fractional coordinates
   (relative to the supercell vectors). Set to `true` for Cartesian coordinates
   in Å.

---

## Combining multiple defects

Defects are applied in a **stable order**: substitutions first, then vacancies,
then interstitials. The `index` fields always refer to the supercell *before*
any atoms are added or removed.

```toml
[geometry.builder]
type             = "crystal"
name             = "Cu"
crystalstructure = "fcc"
a                = 3.61
cubic            = true
supercell        = [3, 3, 3]

[[geometry.builder.defects]]
kind    = "substitution"
index   = 0
element = "Ni"          # swap atom 0 for Ni

[[geometry.builder.defects]]
kind  = "vacancy"
index = 1               # remove atom 1 from the original supercell

[[geometry.builder.defects]]
kind      = "interstitial"
element   = "H"
position  = [0.25, 0.25, 0.25]
cartesian = false
```

---

## Common crystal systems

=== "Metals (FCC)"

    ```toml
    [geometry.builder]
    type             = "crystal"
    name             = "Au"
    crystalstructure = "fcc"
    a                = 4.07
    cubic            = true
    supercell        = [2, 2, 2]
    ```

=== "Metals (BCC)"

    ```toml
    [geometry.builder]
    type             = "crystal"
    name             = "W"
    crystalstructure = "bcc"
    a                = 3.16
    cubic            = true
    supercell        = [2, 2, 2]
    ```

=== "Semiconductors"

    ```toml
    [geometry.builder]
    type             = "crystal"
    name             = "Si"
    crystalstructure = "diamond"
    a                = 5.43
    cubic            = true
    supercell        = [2, 2, 2]
    ```

=== "Ionic crystals"

    ```toml
    [geometry.builder]
    type             = "crystal"
    name             = "NaCl"
    crystalstructure = "rocksalt"
    a                = 5.64
    cubic            = true
    supercell        = [2, 2, 2]
    ```

=== "HCP metals"

    ```toml
    [geometry.builder]
    type             = "crystal"
    name             = "Ti"
    crystalstructure = "hcp"
    a                = 2.95
    c                = 4.68
    supercell        = [3, 3, 2]
    ```

---

## Inspecting defect provenance

Defect information is stored in `provenance.extra`:

```python
from traincraft.geometry import build_geometry
from traincraft.config.models import CrystalBuilder, DefectSpec, GeometryConfig

s = build_geometry(GeometryConfig(
    builder=CrystalBuilder(
        name="Cu", crystalstructure="fcc", a=3.6, cubic=True,
        supercell=(2, 2, 2),
        defects=[DefectSpec(kind="vacancy", index=0)],
    )
))
print(s.provenance.extra)
# {'defects': [{'kind': 'vacancy', 'index': 0, 'element': None}]}
```

---

## Summary

| Defect type | TOML `kind` | Required fields |
|---|---|---|
| Remove an atom | `"vacancy"` | `index` |
| Swap an element | `"substitution"` | `index`, `element` |
| Add an atom | `"interstitial"` | `element`, `position` |

**Next:** [Tutorial 5](05-slabs-strain.md) — standalone surface slabs and mechanical deformation.
