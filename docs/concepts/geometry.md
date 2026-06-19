# Geometry System

The geometry subsystem follows a **Source × Builder × Transform** pattern:
choose one way to get a structure (source or builder), then optionally apply
a sequence of transforms.

---

## Sources vs. builders

| | Source | Builder |
|---|---|---|
| **Starting point** | An existing structure | Constructs from scratch |
| **Examples** | `file`, `scratch`, `smiles`, `url` | `crystal`, `slab`, `layered`, `nanotube` |
| **TOML key** | `[geometry.source]` | `[geometry.builder]` |
| **Use when** | You have a geometry already | You want to generate one |

Exactly one of `source` or `builder` must be set per `[geometry]` section.

---

## Sources

### `file`
Reads any [ASE-supported format](https://wiki.fysik.dtu.dk/ase/ase/io/io.html)
(extxyz, POSCAR, CIF, xyz, …):
```toml
[geometry.source]
type = "file"
path = "POSCAR"
```

### `scratch`
Builds a molecule or bulk crystal from ASE's built-in databases:
```toml
[geometry.source]
type     = "scratch"
molecule = "H2O"    # or: bulk = "Cu"
```

### `smiles`
Generates 3D geometry via RDKit ETKDG + MMFF optimisation:
```toml
[geometry.source]
type        = "smiles"
smiles      = "c1ccccc1"
n_conformers = 1
optimize    = true
vacuum      = 8.0
```

### `url`
Downloads a structure file and reads it with ASE:
```toml
[geometry.source]
type   = "url"
url    = "https://example.com/structure.xyz"
format = "xyz"      # optional; inferred from suffix
```

---

## Builders

| Builder | System | Key parameters |
|---|---|---|
| `nanotube` | Carbon nanotube | `n`, `m`, `length`, `bond` |
| `molecule` | ASE g2 molecule | `name` or `smiles` |
| `crystal` | Bulk crystal | `name`, `crystalstructure`, `a`, `supercell`, `defects` |
| `slab` | Crystalline slab | `element`, `facet` or `miller`, `size`, `layers` |
| `layered` | 2D material stack | `material`, `n_layers`, `stacking`, `twist` |
| `surface_adsorbate` | Molecule on slab | `element`, `facet`, `molecule_name`/`smiles`/`file` |
| `surface_packing` | N molecules on slab | same + `n_molecules`, Packmol |

---

## Transforms

Transforms are composable post-operations declared as a TOML array. They are
applied in order — each one receives the output of the previous one.

```toml
[[geometry.transforms]]
type   = "supercell"
repeat = [2, 2, 1]

[[geometry.transforms]]
type        = "strain"
hydrostatic = 0.02

[[geometry.transforms]]
type   = "perturb"
stddev = 0.05
```

All transforms preserve the `Structure` interface — they return a new
`Structure` with updated provenance (the transform name is appended to
`provenance.transforms`).

| Transform | Effect |
|---|---|
| `supercell` | Tile the cell: `atoms.repeat([nx, ny, nz])` |
| `vacuum` | Add vacuum padding: `atoms.center(vacuum=amount/2)` |
| `perturb` | Gaussian random displacements (seeded) |
| `strain` | Elastic deformation: `hydrostatic` or Voigt 6-vector |
| `rotate` | Rigid rotation about an axis |
| `set_pbc` | Change periodic boundary conditions |

---

## The `GeometryConfig` validator

```python
class GeometryConfig(TCModel):
    source: SourceConfig | None = None
    builder: BuilderConfig | None = None
    transforms: list[TransformConfig] = []
```

Pydantic validates that exactly one of `source` or `builder` is set.
Setting both or neither raises a `ValidationError` immediately — before
any structure is built.
