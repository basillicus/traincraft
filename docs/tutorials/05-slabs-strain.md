# Tutorial 5 · Slabs & Strain

**What you'll learn:** how to build bare crystalline slabs (without adsorbates)
using named facets or arbitrary Miller indices, and how to apply mechanical
strain and other geometric transforms.

**Prerequisites:** Tutorial 4. No extra dependencies.

---

## The standalone slab builder

The `slab` builder creates a clean crystalline surface model — no adsorbate.
All atoms are tagged as framework (`tc_fragment == -1`), meaning they're
treated as immobile substrate atoms by the MC sampler. You would then use
`strain` or `supercell` transforms, or switch to `surface_adsorbate` to put
something on top.

```toml title="examples/13_slab_strain_md.toml"
[run]
name = "cu111_strained"
seed = 13

[geometry.builder]
type    = "slab"
element = "Cu"
facet   = "fcc111"    # (1)
size    = [3, 3, 4]   # (2)
vacuum  = 12.0        # (3)

[[geometry.transforms]]
type  = "strain"
voigt = [0.02, 0.02, 0.0, 0.0, 0.0, 0.0]  # (4)

[calculator]
type = "emt"

[sampling]
type        = "md"
temperature = 400.0
steps       = 200
interval    = 20

[selection]
steps  = ["physicality", "dedup", "diversity"]
budget = 8

[dataset]
path = "dataset"
```

```bash
pixi run example-13
```

1. **`facet`** — named low-index facets: `fcc111`, `fcc100`, `fcc110`,
   `bcc110`, `bcc100`, `bcc111`, `hcp0001`.

2. **`size`** — `[nx, ny, nz]`. For a slab: `nx` and `ny` repeat the surface
   unit cell in-plane; `nz` (or `layers`) sets the slab thickness.

3. **`vacuum`** — total vacuum gap in Å (added above and below the slab).

4. **Biaxial strain** via the Voigt notation — see below.

---

## Miller index slabs

For any surface orientation, specify Miller indices directly:

```toml
[geometry.builder]
type             = "slab"
element          = "Fe"
miller           = [1, 1, 0]    # (1)
crystalstructure = "bcc"        # (2)
a                = 2.87         # (3)
layers           = 8            # (4)
size             = [2, 2, 0]    # (5)
vacuum           = 15.0
```

1. **`miller`** — the surface Miller indices as a list of 3 integers.
   When `miller` is set, `facet` is ignored.

2. **`crystalstructure`** — needed to construct the bulk cell to cleave.

3. **`a`** — lattice constant. If omitted, ASE uses its tabulated value.

4. **`layers`** — number of atomic layers (used in Miller mode; `nz` from
   `size` is used in named-facet mode).

5. **`size`** — `[nx, ny, 0]` in Miller mode; the third element is ignored.

!!! note "Why does Miller mode give fewer atoms?"
    The surface unit cell for a high-index surface is often smaller than for
    a named facet. Use a larger `size` to get a reasonable supercell.

---

## Transforms

Transforms are composable post-processing steps applied after the builder.
They're specified as an array of `[[geometry.transforms]]` sections.

### `strain` — elastic deformation

Deforms the cell and scales atom positions with it (fractional coordinates
are preserved — this is a proper elastic deformation, not just repositioning).

**Hydrostatic (isotropic) strain:**

```toml
[[geometry.transforms]]
type        = "strain"
hydrostatic = 0.02     # +2% in all three directions
```

**Voigt notation** — `(e_xx, e_yy, e_zz, e_yz, e_xz, e_xy)`:

```toml
[[geometry.transforms]]
type  = "strain"
voigt = [0.02, 0.02, 0.0, 0.0, 0.0, 0.0]   # biaxial in-plane tensile
```

```toml
[[geometry.transforms]]
type  = "strain"
voigt = [0.0, 0.0, -0.01, 0.0, 0.0, 0.0]   # 1% uniaxial compression along z
```

!!! tip "Sign convention"
    Positive strain = expansion. Negative strain = compression.
    A hydrostatic strain of `0.02` increases each lattice parameter by 2%.

### `rotate` — rigid rotation

```toml
[[geometry.transforms]]
type        = "rotate"
angle       = 45.0      # degrees
axis        = "z"       # "x", "y", "z" or a 3-vector [x, y, z]
rotate_cell = false     # (1)
```

1. **`rotate_cell`** — if `true`, also rotates the periodic cell vectors (a
   rigid reorientation of the whole system). If `false`, atoms rotate within
   the existing cell box.

### `set_pbc` — change periodicity flags

```toml
[[geometry.transforms]]
type = "set_pbc"
pbc  = [true, true, true]   # make fully periodic
```

Or to remove all periodicity (cluster/molecule):

```toml
[[geometry.transforms]]
type = "set_pbc"
pbc  = false
```

---

## Combining transforms

Transforms are applied left-to-right (top to bottom in the TOML):

```toml
[[geometry.transforms]]
type   = "supercell"
repeat = [2, 2, 1]        # 1. double in-plane

[[geometry.transforms]]
type        = "strain"
hydrostatic = -0.01       # 2. compress 1% hydrostatically

[[geometry.transforms]]
type   = "perturb"
stddev = 0.03             # 3. add thermal-like noise
```

---

## Generating a strain dataset

A common workflow in computational materials science is to compute energies and
forces at several strain values to fit an equation of state or elastic constants.

```toml
[run]
name = "cu_eos_point"    # change this for each strain value in a scan
seed = 42

[geometry.builder]
type             = "crystal"
name             = "Cu"
crystalstructure = "fcc"
a                = 3.61
cubic            = true
supercell        = [2, 2, 2]

[[geometry.transforms]]
type        = "strain"
hydrostatic = 0.02         # vary this from -0.05 to +0.05

[calculator]
type = "emt"

[sampling]
type        = "rattle"
n_structures = 5
std          = 0.03

[selection]
steps  = ["physicality", "dedup"]
budget = 5

[dataset]
path = "dataset"
```

---

## Summary

| Transform | Effect | Key parameters |
|---|---|---|
| `strain` | Elastic deformation | `hydrostatic` or `voigt` (6-vector) |
| `rotate` | Rigid rotation | `angle`, `axis`, `rotate_cell` |
| `set_pbc` | Change boundary conditions | `pbc` (bool or 3-bool) |
| `supercell` | Tile the cell | `repeat` (3-int) |
| `vacuum` | Add vacuum padding | `amount` (Å) |
| `perturb` | Random displacements | `stddev` (Å) |

**Next:** [Tutorial 6](06-2d-materials.md) — layered 2D materials and moiré stacks.
