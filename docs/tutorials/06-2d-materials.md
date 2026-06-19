# Tutorial 6 · 2D Materials

**What you'll learn:** how to build graphene, hexagonal boron nitride (hBN),
MX₂ transition-metal dichalcogenides, and vertically stacked bilayers — including
twisted moiré geometries.

**Prerequisites:** Tutorial 4. No extra dependencies.

---

## The `layered` builder

```toml title="examples/14_graphene_bilayer_md.toml"
[run]
name = "graphene_bilayer"
seed = 14

[geometry.builder]
type                = "layered"
material            = "graphene"    # (1)
size                = [3, 3]        # (2)
n_layers            = 2             # (3)
stacking            = "AB"          # (4)
interlayer_spacing  = 3.35          # (5)
vacuum              = 15.0          # (6)

[calculator]
type = "emt"

[sampling]
type        = "md"
temperature = 300.0
steps       = 200
interval    = 20

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 8

[dataset]
path = "dataset"
```

```bash
pixi run example-14
```

---

## Parameters explained

1. **`material`** — the 2D crystal type:

    | Value | Description | Default `a` (Å) |
    |---|---|---|
    | `graphene` | Carbon honeycomb | 2.46 |
    | `hbn` | Hexagonal boron nitride | 2.50 |
    | `mx2` | Transition-metal dichalcogenide | 3.18 |

2. **`size`** — `[nx, ny]` in-plane repetitions of the primitive monolayer cell.
   `[3, 3]` → 9 primitive cells in the layer = 18 C atoms per graphene layer.

3. **`n_layers`** — number of layers to stack.

4. **`stacking`** — how layers are stacked:
   - `"AA"` — atoms sit directly on top of each other.
   - `"AB"` — Bernal stacking (the physical ground state of graphite). Alternate
     layers are shifted by one NN vector, `(a₁ + a₂) / 3`.

5. **`interlayer_spacing`** — vertical distance between consecutive layers (Å).
   Graphite: 3.35 Å. MoS₂ bulk: ~6.15 Å.

6. **`vacuum`** — total vacuum above and below the stack (Å). 15 Å avoids
   interaction with periodic images along z.

---

## Material-specific examples

=== "Graphene monolayer"

    ```toml
    [geometry.builder]
    type     = "layered"
    material = "graphene"
    size     = [4, 4]
    n_layers = 1
    vacuum   = 15.0
    ```

=== "Hexagonal BN bilayer (AA')"

    hBN naturally stacks in AA' order (B over N). Use `AB` stacking to approximate
    the experimental structure:

    ```toml
    [geometry.builder]
    type               = "layered"
    material           = "hbn"
    a                  = 2.504
    size               = [3, 3]
    n_layers           = 2
    stacking           = "AB"
    interlayer_spacing = 3.30
    vacuum             = 15.0
    ```

=== "MoS₂ bilayer"

    Molybdenum disulfide — each layer is 3 atoms thick (S–Mo–S sandwich):

    ```toml
    [geometry.builder]
    type               = "layered"
    material           = "mx2"
    formula            = "MoS2"    # or "WS2", "WSe2", "MoSe2", "MoTe2"
    a                  = 3.18
    size               = [2, 2]
    n_layers           = 2
    interlayer_spacing = 6.15
    vacuum             = 15.0
    ```

=== "WSe₂ trilayer"

    ```toml
    [geometry.builder]
    type               = "layered"
    material           = "mx2"
    formula            = "WSe2"
    a                  = 3.28
    size               = [2, 2]
    n_layers           = 3
    interlayer_spacing = 6.48
    vacuum             = 15.0
    ```

---

## Moiré stacks via twist angle

When you specify a non-zero `twist`, each layer is rotated by an additional
`twist` degrees relative to the previous one. Because moiré supercells are
generally incommensurate with any small periodic cell, TrainCraft returns a
**finite, non-periodic flake** padded with vacuum:

```toml
[geometry.builder]
type    = "layered"
material = "graphene"
size    = [6, 6]        # (1)
n_layers = 2
twist   = 21.787        # (2)
vacuum  = 15.0
```

1. **`size`** controls the flake size. Larger `size` → more atoms → larger moiré period.

2. **`twist`** — the twist angle in degrees. For graphene, magic angles occur near
   1.08°, 1.74°, … but small angles require enormous supercells. Use larger angles
   (e.g. 21.8°, 38.2°) for feasible training-set generation.

```python
from traincraft.config.models import LayeredBuilder, GeometryConfig
from traincraft.geometry import build_geometry

s = build_geometry(GeometryConfig(builder=LayeredBuilder(
    material="graphene", size=(6, 6), n_layers=2, twist=21.787
)))
print(s.atoms.get_pbc())            # [False False False] — non-periodic flake
print(s.provenance.extra["twist"])  # 21.787
print(len(s.atoms))                 # 144 atoms (6×6×2×2 C atoms, before vacuum)
```

!!! warning "Non-periodic structures and EMT/tblite"
    EMT works fine for non-periodic flakes. tblite also handles them. MACE-MP0
    is designed for periodic systems — for large non-periodic flakes, use
    MACE-OFF23 instead.

---

## Controlling the lattice constant

Use `a` to override the default lattice constant:

```toml
[geometry.builder]
type     = "layered"
material = "graphene"
a        = 2.46          # Å — close to DFT-PBE value; ASE default is also 2.46
size     = [4, 4]
n_layers = 2
stacking = "AB"
```

This is useful when you want to match an experimentally measured or DFT-relaxed
lattice constant.

---

## Summary

| `material` | Species | Available formulas |
|---|---|---|
| `graphene` | C | — |
| `hbn` | B, N | — |
| `mx2` | M, X₂ | `MoS2`, `WS2`, `WSe2`, `MoSe2`, `MoTe2`, … |

**Next:** [Tutorial 7](07-mace-mp0.md) — using MACE-MP0 as the energy/force engine.
