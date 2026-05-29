# TrainCraft — Design

> Status: proposed architecture for the greenfield rewrite. The legacy code
> (git `e39647c`) is the citable baseline for the Raman work and is being
> replaced, not extended. No backward compatibility is maintained.

## 1. What TrainCraft is

TrainCraft builds **training datasets for, and trains, machine-learned
interatomic potentials (MLIPs)** — with a focus on MACE — and closes the loop
with an **active-learning** cycle that only spends expensive DFT on structures
that actually improve the model.

It must learn and validate **energies, forces, dipoles, and polarizabilities**
(the last two enable IR/Raman spectra, the original scientific driver).

## 2. Design principles

1. **Pure, parameterized functions.** Every step takes typed inputs and returns
   typed outputs/artifacts. No import-time side effects, no global config, no
   `os.chdir`. This is what makes the pipeline testable, library-usable, and
   parallelizable.
2. **Plugins behind a registry.** New geometry builders, calculators, samplers,
   selectors, and models register themselves; the engine looks them up by name.
   Adding a capability = adding one file, never editing a dispatcher.
3. **One spine.** The **dataset + selection layer** is the center of the project.
   Geometry generation feeds candidates in; the active-learning loop is just
   *explore → select → label → train → repeat* over that layer.
4. **Orchestration is the last, swappable layer.** A trivial local engine runs
   everything now; QuACC (or Covalent/Parsl/Dask) drops in later without
   touching the science.
5. **Provenance everywhere.** Every structure and dataset frame records how it
   was made (source, builder, transforms, calculator, model version, seed).
6. **Config is data.** Validated pydantic models, serializable to/from TOML —
   the same serialized form a future node editor would emit.

## 3. High-level architecture

```
            ┌─────────────────────────────────────────────────────────┐
            │                    ORCHESTRATION                          │
            │     local engine  →  QuACC / Covalent / Parsl / Dask      │
            └─────────────────────────────────────────────────────────┘
                                       │ (wires pure functions into a DAG)
   GEOMETRY                            ▼                         TRAINING
 ┌───────────┐   candidates   ┌───────────────┐   labeled    ┌───────────┐
 │ sources   │──────────────► │  DATASET +    │ ───────────► │  MACE     │
 │ builders  │                │  SELECTION    │   frames     │  multi-   │
 │ transforms│ ◄───explore─── │   (SPINE)     │ ◄──train──── │  head     │
 └───────────┘   (MD/MC)      └───────────────┘              └───────────┘
       ▲                          │      ▲                        │
       │                    select │      │ label (DFT)            │ metrics
       │                          ▼      │                        ▼
       │                   ┌──────────────┐               ┌──────────────┐
       └───────────────────│  SAMPLING    │               │  VALIDATION  │
                           │  MD / MC /   │               │ parity, MD,  │
                           │  rattle      │               │ IR/Raman     │
                           └──────────────┘               └──────────────┘
                                  │                               ▲
                           ┌──────────────┐                       │
                           │ CALCULATORS  │───────────────────────┘
                           │ MLIP + DFT   │
                           └──────────────┘
```

## 4. Package layout

```
src/traincraft/
  __init__.py            # curated public API
  config/                # pydantic models per stage + TOML loader + validation
  core/
    structure.py         # Structure: ase.Atoms + provenance + properties
    registry.py          # generic plugin registry + @register decorators
    workspace.py         # explicit Run/Job dirs (replaces all os.chdir)
    results.py           # typed results (energy, forces, dipole, polarizability, paths)
    provenance.py        # how-was-this-made records
  geometry/
    sources/             # scratch, smiles, file (any ASE format), url, db/optimade/MP
    builders/            # crystal, defect, surface, layered, intercalation,
                         #   adsorbate, polymer (PySoftK), nanotube, liquid/packmol
    transforms/          # supercell, strain, rotate, perturb, vacuum, pbc, constraints
  calculators/
    base.py              # CalculatorFactory protocol; property capabilities
    mlip.py              # mace (+ registry for MatterSim/Orb/SevenNet/CHGNet); tblite, xtb
    dft.py               # qe, aims; E/F/stress + dipole + polarizability (DFPT)
  sampling/              # md.py, monte_carlo.py, rattle.py  (Sampler plugins)
  selection/             # physicality, dedup, uncertainty, diversity (FPS), budget
  datasets/              # extxyz/db IO, dedup, splitting, provenance, HEALTH tooling
  training/              # MACE fine-tune/train wrapper, multi-head config, metrics
  validation/            # parity, learning curves, MD stability, EOS/phonons, IR/Raman
  active_learning/       # the loop + convergence criteria
  orchestration/         # engine adapters: local (now), quacc/covalent/... (later)
  cli.py                 # thin Typer shell over the public API
tests/
```

## 5. Core abstractions

### 5.1 `Structure`
A light wrapper around `ase.Atoms` that carries:
- the atoms (positions, cell, pbc, constraints)
- computed **properties** (energy, forces, stress, dipole, polarizability) when present
- **provenance** (source/builder/transform chain, seed, parent ids)
- a stable content hash (for dedup and idempotent jobs)

It converts cleanly to/from `ase.Atoms`, `pymatgen.Structure`, and RDKit mols
(via `core`/`converter` helpers) so each builder can use the best library.

### 5.2 Registry
```python
@register("builder", "surface")
def build_surface(cfg: SurfaceConfig) -> Structure: ...
```
Generic registries exist for: `source`, `builder`, `transform`, `calculator`,
`sampler`, `selector`, `model`. Config names the plugin; the engine resolves it.
Capabilities (e.g. "this calculator can produce polarizability") are declared on
registration so the engine can validate a requested workflow up front.

### 5.3 `Workspace` / `Job`
Owns an **absolute** directory. Every step receives its output dir and returns
results + artifact paths. There is **no CWD mutation anywhere**. Jobs are keyed
by a hash of their inputs so reruns skip completed work (idempotency/resume).

### 5.4 `Result`
Typed container for outputs of a calculation: requested properties, their values,
convergence status, wall time, calculator/model identity, and artifact paths.

## 6. Data model & properties

Dataset frames are stored as **extended XYZ** (plus an ASE-db index) with a fixed
schema for per-structure (`energy`, `dipole`, `polarizability`, `stress`) and
per-atom (`forces`) properties, all tagged with provenance and the labeling
method/level of theory. A property is treated consistently across four layers:
**request → label (DFT) → train (MACE head) → validate**.

### Property support matrix

| Property        | DFT labeling                         | MACE training | Validation               |
|-----------------|--------------------------------------|---------------|--------------------------|
| Energy          | SCF (QE/AIMS)                        | yes           | parity, EOS, learning curve |
| Forces          | SCF (QE/AIMS)                        | yes           | parity (per element), MD stability |
| Stress          | SCF (QE/AIMS)                        | yes           | EOS / elastic            |
| Dipole          | SCF (AIMS `output dipole`; QE)       | yes (head)    | parity, **IR spectrum**  |
| Polarizability  | **linear response / DFPT** (AIMS DFPT dielectric; QE `ph.x`) | yes (head) | parity, **Raman spectrum** |

> Polarizability is the expensive one (linear-response DFT). The labeler models
> it as a distinct, costlier task, and selection accounts for that cost.

## 7. Geometry subsystem (the priority) — Source × Builder × Transform

Composable in three layers so any combination works:

- **Sources** (ingest → `Structure`): from scratch (`ase.build`), from **SMILES**
  (RDKit ETKDG + optimize), from **file** (`ase.io.read`, any format), from
  **URL**, and from structured providers (Materials Project, OPTIMADE, PubChem).
- **Builders** (system types): molecules & conformer ensembles; bulk crystals
  (`ase.build.bulk`/spacegroup/pymatgen) + **defects** (vacancy/substitution/
  interstitial); **surfaces/slabs**; **layered/2D** (stacking, spacing, twist);
  **intercalation** (guest between layers); **molecules on surfaces** (adsorption
  sites); **polymers** (PySoftK); **liquids/mixtures/confined** (Packmol — the
  legacy nanotube filling ported here); **nanotubes**.
- **Transforms** (composable post-ops): supercell, strain, rotate/translate,
  perturb, vacuum/cell, set pbc, constraints. Fixes the legacy "constraints lost
  after Packmol" bug via index-based reapplication.

A geometry workflow is a declared `source → builder → [transforms]` pipeline.

## 8. Calculators

A `CalculatorFactory` protocol returns an ASE calculator and declares which
properties it can produce.

- **MLIP** (`mlip.py`): MACE first-class (foundation + fine-tuned models, multi-
  head). `tblite`/`xtb` kept for cheap exploration. A `model` registry makes
  MatterSim/Orb/SevenNet/CHGNet one-file additions. **ANI and NEP removed.**
- **DFT** (`dft.py`): QE and FHI-AIMS, each producing E/F/stress and, on request,
  dipole and polarizability (DFPT). Labeled results are written to extxyz with
  provenance — fixing the legacy gap where the main pipeline never wrote labels.

## 9. Sampling

`md` (Langevin/NVT/NVE), `monte_carlo`, and `rattle` (HiPhive) as `Sampler`
plugins driven by an MLIP. Samplers produce candidate frames; they do **not**
decide what gets labeled — that is the selection layer's job.

## 10. Selection funnel (solves the redundant-DFT problem)

Runs entirely on cheap compute **before any DFT job is dispatched**:

```
raw frames
  → 1. PHYSICALITY  (min interatomic dist; |E|/|F| sanity from MLIP)
  → 2. DEDUP        (near-duplicate removal in descriptor space: SOAP/ACE/MACE latents)
  → 3. UNCERTAINTY  (committee/ensemble disagreement; keep informative frames)
  → 4. DIVERSITY    (farthest-point / max-min sampling to cover the distribution)
  → 5. BUDGET CAP   (top-N per iteration)  → DFT labeling
```

Each stage is a `selector` plugin; the funnel is configurable and reorderable.

## 11. Datasets & health tooling

- **IO:** extxyz + ASE db; dedup; stratified train/val/test split with leakage checks.
- **Health (pre-training):** composition / chemical-space / volume / density
  coverage maps; energy & per-element force distributions with outlier flags;
  redundancy report (shared with selection); **extrapolation grade** (distance of
  each frame from the training distribution).

## 12. Training

A MACE fine-tune/train wrapper (`mace_run_train`, `--foundation_model` for
fine-tuning) with explicit **multi-head** configuration (energy/forces +
optional dipole + polarizability heads), checkpoint/continue support, and
emitted metrics. Built so a new `model` backend implements the same train/eval
interface.

## 13. Validation / potential quality

- Per-property **parity** + RMSE/MAE (E, F, dipole, polarizability), per element.
- **Learning curves** (error vs dataset size) — informs "good enough".
- Physical checks: NVE energy conservation, MD stability, RDF, EOS/phonons.
- **Spectra validation:** reconstruct IR (dipole autocorrelation) and Raman
  (polarizability autocorrelation) from MLIP-driven MD; compare to DFT/experiment.
- These metrics feed the active-learning convergence criterion.

## 14. Active-learning loop

```
seed dataset (few datasets / files)
  → fine-tune MACE foundation model
    → EXPLORE  (MD/MC on target systems; parallel fan-out)
      → SELECT (the §10 funnel)
        → LABEL (DFT on selected; parallel fan-out)
          → APPEND → RETRAIN / continue fine-tune
            → CONVERGENCE? (val force-RMSE & spectral error thresholds)
              ├─ no  → loop
              └─ yes → done
```

Exploration and labeling are embarrassingly parallel; training is one big GPU
job. This is also where "train while generating" naturally falls out.

## 15. Orchestration

Stage functions are pure, so orchestration is pluggable:
- **`local`** (now): plain Python, serial or threadpool — for dev/tests.
- **QuACC** (planned): atomistic-native, engine-agnostic (Covalent/Parsl/Dask/
  Prefect underneath), results store, HPC submission. Chosen later by simplicity.

The adapter only wires existing pure functions into a DAG; no science lives here.

## 16. Configuration

A single validated config (pydantic → TOML) with one section per stage
(`geometry`, `calculator`, `sampling`, `selection`, `training`,
`active_learning`, `orchestration`). Fail-fast validation with clear messages.
The serialized form is exactly what a node editor would produce/consume.

## 17. Library use & public API

`import traincraft` exposes a curated surface: build a `Structure`, attach a
calculator, sample, select, label, train, validate — each usable standalone in a
user's own scripts. The CLI is a thin shell over this API; nothing in the CLI is
privileged.

## 18. Future: node-based workflow editor

Because every step is a registered function with typed input/output configs, a
workflow is a serializable DAG of `(node = registered function, ports = typed
configs)`. A front-end (React-Flow/Rete.js) emits that DAG as the same config the
engine already runs — a clean add-on, not a rewrite. Deferred, but kept feasible.

## 19. Explicitly out / removed vs legacy

- `os.chdir`-driven IO, import-time config singleton, broken package layout.
- ANI (poor quality) and NEP (CUDA-only) calculators.
- Scattered if/elif calculator dispatch (replaced by registry).
- Standalone post-processing scripts (folded into `datasets`).
