# Architecture

TrainCraft is designed around three principles: **pure functions, a plugin
registry, and config as data**. Understanding these helps you extend the
system confidently.

---

## The pipeline

```
     GEOMETRY                    SPINE                     TRAINING
  ┌──────────┐   candidates  ┌─────────────┐  labeled   ┌──────────┐
  │ sources  │──────────────►│  DATASET  + │───────────►│  MACE    │
  │ builders │               │  SELECTION  │   frames   │  multi-  │
  │ transforms│◄──explore────│   (SPINE)   │◄──train────│  head    │
  └──────────┘  (MD/MC)      └─────────────┘            └──────────┘
       ▲               select │     ▲ label (DFT)              │
       │                      ▼     │                           ▼
       │              ┌──────────────┐                 ┌──────────────┐
       └──────────────│  SAMPLING    │                 │  VALIDATION  │
                      │  MD/MC/rattle│                 │  parity, MD, │
                      └──────────────┘                 │  IR/Raman    │
                             │                         └──────────────┘
                      ┌──────────────┐
                      │ CALCULATORS  │
                      │ potentials + │
                      │ DFT          │
                      └──────────────┘
```

The **dataset + selection layer** is the spine. Everything else feeds into it
or draws from it.

---

## Pure, parameterised functions

Every stage takes typed config objects and returns typed outputs:

```python
build_geometry(GeometryConfig) → Structure
make_calculator(CalculatorConfig) → ASE Calculator
run_sampling(structure, calc, job, SamplingConfig) → list[Structure]
run_funnel(frames, SelectionConfig) → list[Structure]
```

No global state. No `os.chdir`. No import-time side effects. This makes stages:

- **Testable** in isolation
- **Composable** — pass the output of one directly to the next
- **Parallelisable** — the orchestration layer is just ordering and fan-out

---

## The plugin registry

Every builder, calculator, sampler, and selector is registered with a decorator:

```python
@register("builder", "crystal")
def build_crystal(cfg: CrystalBuilder) -> Structure:
    ...
```

The engine resolves them by name:

```python
fn = get("builder", cfg.type)   # looks up "crystal" → build_crystal
structure = fn(cfg)
```

Adding a new plugin is **one new file** — no dispatcher to edit, no `__init__`
to update beyond adding an import.

### See what's registered

```bash
traincraft plugins
```

```
source:    file, scratch, smiles, url
builder:   crystal, layered, molecule, nanotube, slab, ...
transform: perturb, rotate, set_pbc, strain, supercell, vacuum
...
```

---

## Config is data

A single TOML file drives the whole workflow. It's validated by pydantic v2 with
`extra="forbid"`, so typos in field names raise an error immediately.

Plugin types are **discriminated unions** on the `type` field:

```toml
[calculator]
type = "mace"          # → MaceCalc model
model = "mace-mp0"
```

```toml
[calculator]
type = "emt"           # → EmtCalc model (no extra fields)
```

The same TOML structure is what a future node-based workflow editor would emit
and consume — the design keeps that door open.

---

## `Structure` — the unit that flows

Every stage passes `Structure` objects. A `Structure` is:

```python
@dataclass
class Structure:
    atoms: ase.Atoms          # geometry, cell, pbc, per-atom arrays
    properties: dict          # energy, forces, stress, dipole, polarizability
    provenance: Provenance    # origin, source, transforms, calculator, seed, ...
```

The `atoms` object carries per-atom arrays (including `tc_fragment` for
fragment identity) that survive `copy()`, extxyz IO, MD/MC, and `repeat()`.

A content **hash** (`sha1` of species + positions + cell, rounded) provides
stable IDs for deduplication across runs.

---

## Workspace — no `os.chdir`

```python
ws = Workspace(Path("runs") / "my_run")
job = ws.job("sampling")     # creates runs/my_run/sampling/
path = job.path("output.extxyz")
```

`Workspace` and `Job` work with absolute paths. Nothing ever calls `os.chdir`.
ASE calculators receive `directory=` explicitly. This lets multiple pipelines
run in the same process without interfering.

---

## Orchestration

Currently: the local engine (`orchestration/local.py`) runs stages serially.
It's intentionally thin — the science lives in pure functions, so swapping in
QuACC, Parsl, or Prefect means replacing only the ordering layer, not touching
any calculator or sampler code.
