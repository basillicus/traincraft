# TrainCraft — Roadmap

Phased, sequenced delivery. The **dataset + selection layer is the spine**:
geometry and the active-learning loop both build on it.

See [`DESIGN.md`](DESIGN.md) for the target architecture.

---

## ✅ Phase 0 — Foundation + vertical slice  *(done, `d43d3b8`)*

**Goal:** a clean, installable, testable package with no globals and no `chdir`,
plus a thin end-to-end walking skeleton proving the architecture.

Delivered
- Correct package (`__init__.py`, relative imports, `pip install -e .`, console
  script), pixi env config (`pyproject.toml [tool.pixi.workspace]`), ruff +
  pytest config, GitHub Actions CI.
- `config/` — pydantic v2 models (discriminated unions, `extra=forbid`) + TOML
  loader. **No import-time singleton.** One TOML still drives the whole workflow.
- `core/` — `Structure` (+hash), `registry`, `Workspace`/`Job` (no `os.chdir`),
  `Result`, `provenance`, seeded RNG.
- `geometry/` — Source×Builder×Transform (file/scratch sources; nanotube/molecule
  builders; vacuum/supercell/perturb transforms).
- `calculators/potentials.py` — `emt` (force field), `tblite`/`xtb`
  (semiempirical), `mace` (MLIP). ANI/NEP dropped. MACE plumbing fixed.
- `sampling/` — `md` (Langevin NVT) and `rattle` (HiPhive) ported; molecule-aware
  `monte_carlo` interface stubbed.
- `selection/` — funnel: `physicality` → `dedup` → `diversity` (FPS) → `budget`.
- `datasets/` — extxyz IO with provenance; `Dataset` (hash-dedup + origin filter).
- `orchestration/local.py` — `run_pipeline` serial engine.
- `cli.py` — `run` / `validate` / `new` / `plugins`.
- `examples/` — 6 annotated examples covering different builders, samplers,
  selection configs, and calculators.
- 19 tests pass, ruff clean, CLI verified end-to-end.

---

## Phase 1 — Geometry breadth  *(next)*

**Goal:** generate all the system types that matter scientifically. Delivered in
chunks so each is independently testable and reviewable.

### Chunk 1 — Molecules-on-surfaces + Monte Carlo  *(spec ready)*

Full implementation spec: [`docs/phase1_chunk1_spec.md`](docs/phase1_chunk1_spec.md).
Designs the hard abstraction (**fragment identity**, per-atom `tc_fragment` array +
runtime `infer_fragments` for reactive runs) first, plus the features that depend
on it:
- `core/fragments.py` — fragment-identity layer (the spine of this chunk).
- Source: `smiles` (RDKit ETKDG + optimize).
- Builders: `surface_adsorbate` (single adsorbate, `ase.build.add_adsorbate`,
  zero deps) and `surface_packing` (N-molecule coverage via the Packmol binary).
- `sampling/monte_carlo.py` — Metropolis MC with rigid-body translate/rotate +
  RDKit conformer-swap moves; optional `refresh_fragments` for bond-changing runs.

### Later Phase-1 chunks

- `converter.py` — ASE ↔ pymatgen ↔ RDKit bridge.
- Sources: `url` (download → read), providers (Materials Project, OPTIMADE, PubChem).
- Builders: bulk crystal + defects (vacancy/substitution/interstitial); surface/
  slab; layered/2D (stacking, interlayer spacing, twist); intercalation; polymers
  (PySoftK); liquids/mixtures/confined bulk (Packmol).
- Transforms: strain, rotate, set pbc, constraints (index-based reapplication after
  Packmol — fixes the legacy bug at `gengeom.py:155`).

Acceptance
- Each builder produces a valid `Structure` with correct pbc/cell and provenance.
- Legacy nanotube+CO2 dataset reproducible via the new geometry path.
- `examples/07_co_on_cu_mc.toml` (zero deps), `08_smiles_molecule.toml`, and
  `09_packing_on_surface.toml` work.

---

## Phase 2 — DFT labeling with full property set

**Goal:** label E/F/stress + dipole + polarizability.

Deliverables
- `calculators/dft.py` — QE and FHI-AIMS factories; SCF for E/F/stress; dipole
  output; **polarizability via DFPT** (AIMS `DFPT dielectric`; QE `ph.x`).
- Labeled results written to extxyz with level-of-theory provenance.
- Cost-aware labeling (polarizability flagged as the expensive task; selection
  accounts for the cost).
- `runs/<name>/labeled_dft/` tree + `manifest.json` (level-of-theory, counts).

Acceptance
- A selected frame is labeled with all requested properties and lands in
  `labeled_dft/` with provenance; QE and AIMS paths verified on a tiny system.

---

## Phase 3 — Training + validation (multi-head MACE)

**Goal:** train MACE on all properties and measure quality end-to-end.

Deliverables
- `training/` — MACE fine-tune/train wrapper (`mace_run_train --foundation_model`),
  explicit multi-head config (energy/forces + optional dipole + polarizability),
  checkpoints, metrics. Pluggable model interface for future backends.
- `datasets/` health tooling: composition/space/volume coverage maps; energy and
  per-element force distributions with outlier flags; extrapolation grade;
  redundancy report.
- `validation/` — per-property parity + RMSE/MAE (per element), learning curves,
  NVE/MD stability, EOS/phonons, **IR and Raman spectra reconstruction** from
  MLIP-driven MD vs DFT/experiment.

Acceptance
- A fine-tuned model trains on a seed set and produces a quality report covering
  every property, including reconstructed IR/Raman spectra.

---

## Phase 4 — Active-learning loop

**Goal:** close the loop to convergence automatically.

Deliverables
- `selection/uncertainty.py` — committee/ensemble uncertainty selector (plugs into
  the existing funnel).
- `active_learning/` — explore → select → label → retrain → converge; thresholds
  drawn from `validation` (val force-RMSE, spectral error); resume/idempotency
  across iterations.

Acceptance
- Starting from a seed set, the loop runs ≥2 iterations unattended; the dataset
  grows only with informative + diverse frames; validation error decreases.

---

## Phase 5 — Orchestration

**Goal:** parallelism without touching the science.

Deliverables
- `orchestration/` — `local` engine hardened (threadpool for independent jobs);
  QuACC adapter expressing the AL loop as a DAG with parallel explore/label fan-out.
  Engine chosen by simplicity after the core is proven.

Acceptance
- The same AL workflow runs under `local` and the chosen engine with identical
  results; explore and label stages fan out in parallel.

---

## Phase 6 — Polish & extras

- Full public API docs + library-usage tutorials (including the Raman use case).
- Additional MLIP backends (MatterSim/Orb/SevenNet/CHGNet) via the model registry.
- Deferred, architecture-friendly: **node-based workflow editor** emitting the
  serialized TOML DAG.

---

## Dependency graph

```
Phase 0 (done)
  └─► Phase 1 (geometry breadth)
  └─► Phase 2 (DFT labeling)   ─► Phase 3 (training + validation)
                                         └─► Phase 4 (AL loop)
                                                   └─► Phase 5 (orchestration)
                                                             └─► Phase 6 (polish)
```

Phases 1 and 2 are independent and can progress in parallel.
