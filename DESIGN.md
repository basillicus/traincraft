# TrainCraft — Design

> **Status:** Phase 0 foundation + vertical slice implemented (`d43d3b8`).
> The legacy code (git `e39647c`) is the citable Raman baseline and is no
> longer part of this repository. No backward compatibility is maintained.

---

## 1. What TrainCraft is

TrainCraft builds **training datasets for, and trains, machine-learned
interatomic potentials (MLIPs)** — with a focus on MACE — and closes the loop
with an **active-learning** cycle that only spends expensive DFT on structures
that actually improve the model.

It must learn and validate **energies, forces, dipoles, and polarizabilities**
(the last two enable IR/Raman spectra — the original scientific driver).

---

## 2. Design principles

1. **Pure, parameterized functions.** Every step takes typed inputs and returns
   typed outputs/artifacts. No import-time side effects, no global config, no
   `os.chdir`. This makes the pipeline testable, library-usable, and parallelizable.
2. **Plugins behind a registry.** New geometry builders, calculators, samplers,
   selectors, and models register themselves with `@register`; the engine looks
   them up by name from config. Adding a capability = one new file, no dispatcher
   to edit.
3. **One spine.** The **dataset + selection layer** is the center. Geometry
   generation feeds candidates in; the active-learning loop is
   *explore → select → label → train → repeat* over that layer.
4. **Orchestration is the last, swappable layer.** A local serial engine runs
   everything now; QuACC (or Covalent/Parsl/Dask) drops in later without
   touching the science.
5. **Provenance everywhere.** Every structure and frame records how it was made
   (source, builder, transforms, calculator, model version, seed). An `origin`
   tag (`generated` → `ml_sampled` → `ml_labeled` → `dft_labeled`) keeps the
   expensive DFT dataset cleanly separable.
6. **Config is data.** Validated pydantic models, serializable to/from TOML —
   the same serialized form a future node editor would emit.

---

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
                           │ potentials + │
                           │ DFT          │
                           └──────────────┘
```

---

## 4. Package layout

```
src/traincraft/
  __init__.py            # curated public API
  cli.py                 # Typer: run / validate / new / plugins
  config/
    models.py            # pydantic v2 models (discriminated unions on `type`)
    loader.py            # TOML → validated TrainCraftConfig
  core/
    structure.py         # Structure: ase.Atoms + properties + provenance + hash
    provenance.py        # Provenance record + origin tags
    registry.py          # generic plugin registry + @register
    workspace.py         # Workspace / Job — absolute dirs, no os.chdir
    results.py           # Result: typed calc outputs
    rng.py               # seeded numpy Generator factory
  geometry/
    sources/             # file (any ASE format), scratch (ase.build)
                         #   Phase 2: SMILES (RDKit), URL, OPTIMADE/MP
    builders/            # nanotube, molecule
                         #   Phase 2: crystal+defects, surface, layered,
                         #            intercalation, adsorbate, polymer (PySoftK),
                         #            liquid/confined (Packmol)
    transforms/          # vacuum, supercell, perturb
                         #   Phase 2: strain, rotate, pbc, constraints
  calculators/
    base.py              # CalculatorFactory protocol + capability declarations
    potentials.py        # cheap calcs — these are NOT all MLIPs:
                         #   emt (force field, zero deps)
                         #   tblite / xtb (semiempirical, GFN-xTB)
                         #   mace (MLIP: foundation + fine-tuned models)
    dft.py               # Phase 3: QE + FHI-AIMS, E/F/stress/dipole/polarizability
  sampling/
    base.py              # Sampler protocol (molecule-aware)
    md.py                # Langevin MD (Langevin NVT via ASE)
    rattle.py            # HiPhive standard/MC rattle
    monte_carlo.py       # Phase 2: Metropolis MC + RDKit conformer moves
  selection/
    base.py              # Selector protocol + run_funnel
    physicality.py       # min interatomic distance filter
    dedup.py             # exact dedup by content hash
    diversity.py         # farthest-point sampling (histogram descriptor now;
                         #   SOAP/ACE via dscribe in Phase 2)
                         #   Phase 5: uncertainty / committee selector
  datasets/
    io.py                # extxyz read/write with provenance
    dataset.py           # Dataset: append (hash-dedup) + filter by origin
                         #   Phase 4: health tooling (coverage, distributions,
                         #            extrapolation grade, redundancy report)
  training/              # Phase 4: MACE fine-tune/train, multi-head config
  validation/            # Phase 4: parity, learning curves, MD stability,
                         #          EOS/phonons, IR/Raman spectra reconstruction
  active_learning/       # Phase 5: explore → select → label → retrain → converge
  orchestration/
    local.py             # serial engine (run_pipeline) — now
                         #   Phase 6: QuACC / Covalent / Parsl adapter
tests/
examples/
  01_cnt_emt_md.toml          # walking skeleton (zero deps)
  02_molecule_emt_rattle.toml # rattle sampling
  03_molecule_from_file.toml  # read geometry from file
  04_nanotube_supercell_rattle.toml
  05_selection_funnel_demo.toml
  06_mace_mp0_sampling.toml   # MACE-MP0 as sampler
```

---

## 5. Core abstractions

### 5.1 `Structure`
A dataclass wrapping `ase.Atoms` + a `properties` dict
(`energy`, `forces`, `stress`, `dipole`, `polarizability`) + a `Provenance` +
a stable content `hash` (numbers, positions, cell rounded). Helpers
`to_ase`/`from_ase`; stubs `to_pymatgen`/`to_rdkit` (Phase 2).

### 5.2 Registry
```python
@register("builder", "nanotube")
def build_nanotube(cfg: NanotubeBuilder) -> Structure: ...
```
Kind → name → `{obj, capabilities}`. Built-in modules self-register on import;
the engine resolves by name from config. `capabilities` (set of property names)
lets the engine validate a requested workflow up front.

### 5.3 `Workspace` / `Job`
`Workspace(root)` owns an absolute directory. `Job` owns a sub-directory with a
done-marker for idempotent rerun. ASE calculators receive `directory=` explicitly.
**Nothing calls `os.chdir`.**

### 5.4 `Result`
Typed container: requested properties, values, convergence status, wall time,
calculator identity, artifact paths.

---

## 6. Data model and properties

Frames are stored as **extended XYZ** (plus an ASE-db index) with a fixed schema.
A property is treated consistently across four layers:
**request → label (DFT) → train (MACE head) → validate**.

| Property | DFT labeling | MACE training | Validation |
|----------|-------------|---------------|------------|
| Energy | SCF (QE/AIMS) | yes | parity, EOS, learning curve |
| Forces | SCF (QE/AIMS) | yes | parity per element, MD stability |
| Stress | SCF (QE/AIMS) | yes | EOS / elastic |
| Dipole | SCF (AIMS `output dipole`; QE) | yes (head) | parity, **IR spectrum** |
| Polarizability | **DFPT** (AIMS dielectric; QE `ph.x`) | yes (head) | parity, **Raman spectrum** |

> Polarizability requires linear-response DFT (heavier than SCF). The labeler
> models it as a distinct, costlier task, and selection accounts for the cost.

---

## 7. Geometry subsystem — Source × Builder × Transform

### Sources (→ `Structure`)
`file` (any ASE-readable format), `scratch` (`ase.build` molecule/bulk). Phase 2
adds: `smiles` (RDKit ETKDG), `url`, providers (Materials Project, OPTIMADE, PubChem).

### Builders (system types)
`nanotube`, `molecule`. Phase 2 adds: bulk crystals + defects (vacancy/substitution/
interstitial), surfaces/slabs, layered/2D (stacking, interlayer spacing, twist),
intercalation, molecules on surfaces (adsorption sites), polymers (PySoftK),
liquids/mixtures/confined (Packmol).

### Transforms (composable post-ops)
`vacuum`, `supercell`, `perturb`. Phase 2 adds: strain, rotate, pbc, constraints
(index-based reapplication after Packmol — fixes the legacy bug).

---

## 8. Calculators

A `CalculatorFactory` protocol returns an ASE calculator and declares which
properties it can produce.

**`potentials.py`** — *not all of these are MLIPs*:
- `emt` — ASE force field, zero dependencies; default for tests and CI.
- `tblite` / `xtb` — semiempirical GFN-xTB, good for molecules.
- `mace` — MLIP; foundation models (`mace-mp0`, `mace-off23`) or a local
  fine-tuned model via `model_path`. Plumbing corrected from the legacy code.
  ANI and NEP are intentionally not ported (poor quality / CUDA-only).

**`dft.py`** (Phase 3): QE and FHI-AIMS, each producing E/F/stress, and on
request dipole and polarizability (DFPT). Results written to extxyz with
level-of-theory provenance.

---

## 9. Sampling

`Sampler` protocol: `run(structure, calc, job, cfg) -> list[Structure]`. Designed
molecule-aware so MC can do rigid-body moves and conformer generation.

- **`md`** — Langevin NVT. Writes `.traj` into the job dir, returns frames.
- **`rattle`** — HiPhive standard/MC rattle (optional dep).
- **`monte_carlo`** — Phase 2: Metropolis MC with rigid-body translation/rotation
  of molecules and RDKit conformer generation (ETKDG). Primary tool for complex
  molecules on surfaces.

---

## 10. Selection funnel

Runs entirely on cheap compute **before any DFT job is dispatched**:

```
raw frames
  → 1. physicality  (min interatomic distance — catches clashing atoms)
  → 2. dedup        (exact hash; near-dup via SOAP/dscribe in Phase 2)
  → 3. uncertainty  (Phase 5: committee/ensemble disagreement)
  → 4. diversity    (farthest-point sampling over histogram descriptor now;
                     SOAP/ACE in Phase 2)
  → 5. budget cap   → DFT labeling
```

Each stage is a `selector` plugin; the funnel is configurable and reorderable
from the TOML (`selection.steps = [...]`).

---

## 11. Datasets and health tooling

- **IO:** extxyz + ASE db; hash-dedup; filter by `origin`/level-of-theory.
- **Phase 4 health tooling:** composition/space/volume/density coverage maps;
  energy and per-element force distributions with outlier flags; redundancy
  report; **extrapolation grade** (distance of each frame from the training
  distribution).

---

## 12. Training

A MACE fine-tune/train wrapper (`mace_run_train`, `--foundation_model`) with
explicit **multi-head** configuration (energy/forces + optional dipole +
polarizability heads), checkpoint/continue support, and emitted metrics. The
model interface is pluggable; adding MatterSim/Orb/SevenNet/CHGNet is one file.

---

## 13. Validation / potential quality

- Per-property **parity** + RMSE/MAE (E, F, dipole, polarizability), per element.
- **Learning curves** (error vs dataset size).
- Physical checks: NVE conservation, MD stability, RDF, EOS/phonons.
- **Spectra validation:** IR (dipole autocorrelation) and Raman (polarizability
  autocorrelation) reconstructed from MLIP-driven MD vs DFT/experiment.
- These metrics feed the active-learning convergence criterion.

---

## 14. Active-learning loop

```
seed dataset
  → fine-tune MACE foundation model
    → EXPLORE  (MD/MC; parallel fan-out)
      → SELECT (§10 funnel)
        → LABEL (DFT; parallel fan-out)
          → APPEND → RETRAIN
            → CONVERGED? (val force-RMSE + spectral error thresholds)
              ├─ no  → loop
              └─ yes → done
```

---

## 15. Orchestration

Stage functions are pure; orchestration is pluggable:
- **`local`** (now): serial or threadpool — for dev/tests.
- **QuACC** (Phase 6): atomistic-native, engine-agnostic (Covalent/Parsl/Dask/
  Prefect), results store, HPC submission. Chosen by simplicity once the core
  is proven.

---

## 16. Configuration

One TOML, validated by pydantic, drives the whole workflow. Stage sections are
optional — presence enables the stage. Plugin types are discriminated unions on
a `type` field, so adding a plugin = adding one model + one registry entry. The
serialized form is exactly what a future node editor would produce/consume.

---

## 17. Library API

`import traincraft` exposes a curated surface:

```python
tc.load_config(path)          # → TrainCraftConfig
tc.build_geometry(cfg)        # → Structure
tc.make_calculator(cfg)       # → ASE calculator
tc.run_sampling(s, calc, job, cfg) # → list[Structure]
tc.run_funnel(frames, cfg)    # → list[Structure]
tc.run_pipeline(cfg)          # → summary dict
tc.write_frames(path, frames) # write extxyz
tc.read_frames(path)          # read extxyz → list[Structure]
```

---

## 18. Future: node-based workflow editor

Because every step is a registered function with typed input/output configs, a
workflow is a serializable DAG of `(node = registered function, ports = typed
configs)`. A React-Flow/Rete.js front-end emitting this DAG as TOML is a clean
add-on, not a rewrite. Deferred, but kept feasible by the architecture.

---

## 19. What is out of scope / removed vs. legacy

- `os.chdir`-driven IO, import-time config singleton, broken package layout.
- ANI (poor quality) and NEP (CUDA-only) calculators.
- Scattered `if/elif` calculator dispatch (replaced by registry).
- Standalone post-processing scripts (folded into `datasets/`).
- The name `mlip.py` — the file is `potentials.py` because EMT (force field)
  and tblite/xtb (semiempirical) are **not** MLIPs.
