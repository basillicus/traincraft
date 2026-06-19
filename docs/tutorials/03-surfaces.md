# Tutorial 3 · Molecules on Surfaces

**What you'll learn:** how to place one or more adsorbate molecules on a
crystalline slab, and how the Monte Carlo sampler explores adsorption sites and
molecular orientations while keeping the substrate fixed.

**Prerequisites:** Tutorial 1. For Packmol (multi-molecule): `pixi install -e science`.

---

## Single adsorbate — `surface_adsorbate`

The simplest surface system: one molecule on a crystalline slab.

```toml title="examples/07_co_on_cu_mc.toml"
[run]
name  = "co_on_cu111"
seed  = 42

[geometry.builder]
type          = "surface_adsorbate"
# --- substrate ---
element       = "Cu"           # (1)
facet         = "fcc111"       # (2)
size          = [3, 3, 4]      # (3)
vacuum        = 12.0           # (4)
# --- adsorbate ---
molecule_name = "CO"           # (5)
site          = "ontop"        # (6)
height        = 2.0            # (7)

[calculator]
type = "emt"

[sampling]
type          = "monte_carlo"  # (8)
steps         = 500
temperature   = 500.0
interval      = 25
p_translate   = 0.6
p_rotate      = 0.4
p_conformer   = 0.0
max_translate = 0.5
max_rotate    = 30.0

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 10

[dataset]
path = "dataset"
```

### Substrate parameters

1. **`element`** — the substrate element (must be in ASE's FCC/BCC/HCP tables).

2. **`facet`** — the crystallographic surface. Available: `fcc111`, `fcc100`,
   `fcc110`, `bcc110`, `bcc100`, `bcc111`, `hcp0001`.

3. **`size`** — `[nx, ny, nz]` repetitions of the surface unit cell.
   `[3,3,4]` → 3×3 surface unit cells, 4 atomic layers.

4. **`vacuum`** — total vacuum thickness in Å. 12 Å is standard; 15–20 Å
   prevents unphysical interaction between periodic images.

### Adsorbate parameters

5. **`molecule_name`** — ASE G2 name. Or use `smiles = "CCO"` for an organic
   adsorbate, or `file = "my_molecule.xyz"` for a pre-optimised geometry.
   Exactly one of these must be set.

6. **`site`** — initial placement site: `ontop`, `bridge`, `hollow`, `fcc`,
   `hcp`.

7. **`height`** — initial height above the topmost substrate atom (Å).

### Fragment tagging

The builder automatically assigns fragment IDs:

```
tc_fragment[i] == -1   → substrate (Cu) atom — never moves in MC
tc_fragment[i] == 0    → adsorbate (CO) atom — sampled by MC
```

You can verify:

```python
from traincraft.geometry.builders.surface import build_surface_adsorbate
from traincraft.config.models import SurfaceAdsorbateBuilder
from traincraft.core.fragments import get_fragments, FRAMEWORK
import numpy as np

cfg = SurfaceAdsorbateBuilder(
    element="Cu", facet="fcc111", size=(3,3,4), molecule_name="CO"
)
s = build_surface_adsorbate(cfg)
frag = get_fragments(s.atoms)
print(np.sum(frag == FRAMEWORK))  # 36 substrate atoms
print(np.sum(frag >= 0))          # 2 CO atoms
```

### Monte Carlo sampler

```toml
[sampling]
type          = "monte_carlo"  # Metropolis MC
steps         = 500            # (1)
temperature   = 500.0          # (2)
interval      = 25             # (3)
p_translate   = 0.6            # (4)
p_rotate      = 0.4            # (4)
p_conformer   = 0.0            # (4)
max_translate = 0.5            # (5)
max_rotate    = 30.0           # (6)
```

1. **`steps`** — total MC moves. Each accepted move is a candidate frame.

2. **`temperature`** — Metropolis acceptance temperature in K. Higher → more
   uphill moves accepted → explores further from the initial geometry.

3. **`interval`** — save a snapshot every N accepted moves.

4. **`p_translate / p_rotate / p_conformer`** — probability of each move type
   (must sum to ≤ 1.0; remaining probability → no-op). Here: 60% translate,
   40% rotate, 0% conformer swap.

5. **`max_translate`** — maximum displacement per translate move (Å).

6. **`max_rotate`** — maximum rotation angle per rotate move (degrees).

!!! tip "Setting move probabilities"
    For a rigid adsorbate like CO: use `p_translate=0.5, p_rotate=0.5, p_conformer=0.0`.
    For a flexible molecule like ethanol: add `p_conformer=0.1` and RDKit
    will propose new conformer geometries via ETKDG.

---

## Multiple adsorbates — `surface_packing`

For realistic surface coverage with multiple molecules, use Packmol to randomly
pack `n` copies of the adsorbate in the region above the slab:

```toml title="examples/09_packing_on_surface.toml"
[run]
name = "water_on_cu111"
seed = 42

[geometry.builder]
type             = "surface_packing"
element          = "Cu"
facet            = "fcc111"
size             = [4, 4, 4]
vacuum           = 20.0
molecule_name    = "H2O"     # (1)
n_molecules      = 6         # (2)
tolerance        = 2.0       # (3)
region_height    = 8.0       # (4)
gap              = 2.0       # (5)

[calculator]
type = "emt"

[sampling]
type        = "monte_carlo"
steps       = 300
temperature = 600.0
interval    = 30
p_translate = 0.5
p_rotate    = 0.5
p_conformer = 0.0

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 8

[dataset]
path = "dataset"
```

1. **`molecule_name`** — each of the 6 copies gets this molecule.
2. **`n_molecules`** — how many molecules to pack above the slab.
3. **`tolerance`** — Packmol's minimum intermolecular distance (Å).
4. **`region_height`** — height of the packing box above the slab (Å).
5. **`gap`** — clearance between the topmost slab atom and the bottom of the packing box.

Each molecule gets a unique fragment ID (0, 1, 2, …, n-1). The MC sampler
moves each fragment independently.

---

## SMILES adsorbates

Use any molecule via SMILES instead of a G2 name:

```toml
[geometry.builder]
type    = "surface_adsorbate"
element = "Cu"
facet   = "fcc111"
size    = [3, 3, 4]
vacuum  = 15.0
smiles  = "CCO"   # ethanol
site    = "ontop"
height  = 2.5
```

The canonical SMILES is stored in `provenance.extra["smiles"]` for traceability.

---

## Conformer move (flexible adsorbates)

For flexible adsorbates like ethanol or butane, enable the conformer move:

```toml
[sampling]
type        = "monte_carlo"
steps       = 1000
temperature = 500.0
interval    = 50
p_translate = 0.4
p_rotate    = 0.4
p_conformer = 0.2    # 20% chance of proposing a new RDKit conformer
max_translate = 0.3
max_rotate    = 20.0
```

The conformer move uses RDKit ETKDG to generate a new geometry for the
adsorbate molecule and proposes it as a Metropolis move.

---

## Summary

| Builder | Adsorbates | Packmol required? |
|---|---|---|
| `surface_adsorbate` | 1 molecule | No |
| `surface_packing` | N molecules | Yes (science env) |

**Next:** [Tutorial 4](04-crystals-defects.md) — bulk crystals and point defects.
