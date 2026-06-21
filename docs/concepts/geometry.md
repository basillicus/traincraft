# Geometry System

The geometry subsystem follows a **Source × Builder × Transform** pattern:
choose one way to get a structure (source or builder), then optionally apply
a sequence of transforms.

---

## Sources vs. builders

| | Source | Builder |
|---|---|---|
| **Starting point** | An existing structure | Constructs from scratch |
| **Examples** | `file`, `scratch`, `smiles`, `url`, `materials_project`, `optimade`, `pubchem` | `crystal`, `slab`, `layered`, `nanotube`, `filled_nanotube`, `liquid`, `intercalation` |
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

### `materials_project`
Fetches a bulk crystal by Materials Project id (needs `mp-api` + `pymatgen`;
API key from config or `$MP_API_KEY`):
```toml
[geometry.source]
type        = "materials_project"
material_id = "mp-149"     # Si
conventional = true
```

### `optimade`
Returns the first hit of an [OPTIMADE](https://www.optimade.org/) filter against
any provider base URL (dependency-free — the JSON is parsed directly):
```toml
[geometry.source]
type     = "optimade"
base_url = "https://optimade.materialsproject.org/v1"
filter   = 'chemical_formula_reduced="SiO2"'
```

### `pubchem`
Downloads a 3D conformer by name / CID / SMILES from the PubChem PUG REST API
(dependency-free — an SDF is fetched and read by ASE):
```toml
[geometry.source]
type = "pubchem"
name = "caffeine"          # or: cid = 2519 / smiles = "..."
```

---

## Builders

| Builder | System | Key parameters |
|---|---|---|
| `nanotube` | Carbon nanotube (empty) | `n`, `m`, `length`, `bond` |
| `filled_nanotube` | CNT filled with molecules ("fillMyTubes") | `n`, `m`, `length`, single molecule **or** `species` mixture, `n_molecules`, Packmol |
| `molecule` | ASE g2 molecule | `name` or `smiles` |
| `crystal` | Bulk crystal (optionally a mixed solid) | `name`, `crystalstructure`, `a`, `supercell`, `defects`, `composition` |
| `slab` | Crystalline slab (optionally a mixed solid) | `element`, `facet` or `miller`, `size`, `layers`, `composition` |
| `layered` | 2D material stack | `material`, `n_layers`, `stacking`, `twist` |
| `surface_adsorbate` | Molecule on slab | `element`, `facet`, `molecule_name`/`smiles`/`file`, `composition` |
| `surface_packing` | Molecule(s) on slab | same + single molecule **or** `species` mixture, `n_molecules`, `composition`, Packmol |
| `liquid` | Liquid / mixture / confined bulk | `species`, `box` *or* `density`, Packmol |
| `intercalation` | Guests in a layered host | `host`, `guest`, `n_per_gallery`, `stage` |

---

## Mixtures — one `species` concept, everywhere

Every builder that **places molecules** — `liquid`, `surface_packing` and
`filled_nanotube` — speaks the same mixture language. A mixture is a list of
**species**, each an *identity* plus an *amount*:

- **identity**: exactly one of `molecule_name` (ASE g2), `smiles` (RDKit), or
  `file` (any ASE-readable structure);
- **amount**: an absolute `count`, **or** a relative `ratio`. Use one style per
  mixture — all counts, or all ratios apportioned from the builder's total
  `n_molecules` (largest-remainder rounding, so the parts always sum exactly).

A single molecule is just the one-species shortcut: `molecule_name`/`smiles`/
`file` directly on the builder (then `n_molecules` is how many copies). So these
two are equivalent fillings of a nanotube:

```toml
# shortcut: 6 waters
[geometry.builder]
type = "filled_nanotube"
molecule_name = "H2O"
n_molecules = 6
```

```toml
# mixture: 3 parts water, 1 part methane, 8 total -> 6 H2O + 2 CH4
[geometry.builder]
type = "filled_nanotube"
n_molecules = 8

[[geometry.builder.species]]
molecule_name = "H2O"
ratio = 3

[[geometry.builder.species]]
molecule_name = "CH4"
ratio = 1
```

The same `[[...species]]` block drops into `surface_packing` (a mixed adlayer)
and `liquid` (a solvent blend) unchanged. Every packed molecule becomes its own
`tc_fragment`, and `provenance.extra["fragment_species"]` records which species
each fragment is — so the MC sampler, constraints and downstream analysis can
tell water from ethanol.

The `liquid` builder packs the mixture into a periodic box — either a literal
`box = [a, b, c]` (Å) or a cubic cell sized from a target mass `density`
(g/cm³):

```toml
[geometry.builder]
type    = "liquid"
density = 1.0        # g/cm^3  (or: box = [15, 15, 15])
tolerance = 2.0

[[geometry.builder.species]]
molecule_name = "H2O"
count = 32
```

---

## Mixed solids — `composition` (random solid solution)

The lattice analogue of a molecular mixture is a **random substitutional alloy**.
Any lattice builder — `crystal`, `slab`, `surface_adsorbate`, `surface_packing` —
takes a `composition` list that swaps a fraction of host sites for other
elements (chosen with the builder's `seed`, ratios apportioned over sites):

```toml
[geometry.builder]
type             = "crystal"
name             = "Cu"
crystalstructure = "fcc"
cubic            = true
supercell        = [3, 3, 3]
seed             = 11

[[geometry.builder.composition]]
element = "Au"
ratio   = 0.25       # 25 % of the Cu sites become Au (remainder stays Cu)
```

Because mixtures and mixed solids compose, the richer combinations fall out for
free: a solvent blend with a dissolved species (`liquid` mixture), a mixed
adlayer on an **alloy** surface (`surface_packing` with both `species` and
`composition`), and so on.

The `intercalation` builder inserts guest atoms on an in-plane grid into each
gallery of a **planar** layered host (graphene / hBN; `mx2` rejected). `stage`
controls electrochemical staging (fill every `stage`-th gallery):

```toml
[geometry.builder]
type          = "intercalation"
guest         = "Li"
n_per_gallery = 2
stage         = 1

[geometry.builder.host]
material = "graphene"
n_layers = 3
```

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
| `constraints` | Freeze atoms (`FixAtoms`) by `indices`/`elements`/`fragments`/`below_z` |

`constraints` selects on the **final** structure (selectors are OR-ed), which is
the correct place to (re)apply constraints after a builder that reorders atoms
(e.g. Packmol) — for example, freezing the bottom of a slab before sampling:

```toml
[[geometry.transforms]]
type    = "constraints"
below_z = 7.0       # freeze every atom with Cartesian z below 7 Å
```

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
