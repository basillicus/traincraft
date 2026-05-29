# TrainCraft — Roadmap

Phased, sequenced delivery. Each phase lists deliverables and acceptance
criteria. The **dataset + selection layer is the spine**: geometry and the
active-learning loop both build on it, so it lands right after the foundation.

See `DESIGN.md` for the target architecture.

---

## Phase 0 — Foundation (unblocks everything)

**Goal:** a clean, installable, testable package with no globals and no `chdir`.

Deliverables
- Correct package: `__init__.py`, package-relative imports, `pip install -e .`
  works, console script works.
- Dependencies declared in `pyproject.toml`; `environment.yml` for conda/mamba.
- `config/` — pydantic models + TOML loader + fail-fast validation. **No
  import-time singleton**; config is injected.
- `core/` — `Structure`, `registry`, `Workspace`/`Job` (absolute dirs, no
  `chdir`), `Result`, `provenance`.
- pytest + minimal CI; a tiny end-to-end smoke test of one geometry → MLIP path.
- Fix legacy correctness bugs that survive into the new code (label-writing,
  QE vs AIMS handling, MACE param plumbing, no bare `except`).

Acceptance
- `pip install -e .` + `traincraft --help` work from a clean env.
- Tests green in CI; no module mutates global state on import; no `os.chdir`.

---

## Phase 1 — Dataset + selection spine

**Goal:** the center both tracks depend on.

Deliverables
- `datasets/` — extxyz + ASE-db IO, content-hash dedup, stratified split,
  provenance on every frame.
- `selection/` — funnel plugins: physicality, dedup, uncertainty (committee),
  diversity (FPS), budget cap; configurable order.
- `calculators/mlip.py` — MACE (foundation + fine-tuned, multi-head aware),
  `tblite`/`xtb` for cheap work. Registry for future models.
- `sampling/` — `md`, `monte_carlo`, `rattle` as plugins.

Acceptance
- Given a pile of frames, the funnel returns a physical, deduped, informative,
  diverse, budget-capped subset, with a redundancy report.
- A sampler → selector → dataset round-trip runs locally with provenance intact.

---

## Phase 2 — Geometry subsystem (priority breadth)

**Goal:** generate the system types you care about, from any source.

Deliverables (incremental, each with a test)
- Source × Builder × Transform framework; `converter` (ASE/pymatgen/RDKit).
- Sources: scratch, file (any ASE format), SMILES (RDKit), URL, providers
  (Materials Project / OPTIMADE / PubChem).
- Builders: molecules/conformers → crystals + defects → surfaces → layered/2D →
  intercalation → adsorbates (mol-on-surface) → polymers (PySoftK) → liquids/
  confined (Packmol, legacy nanotube ported) → nanotubes.
- Transforms: supercell, strain, rotate, perturb, vacuum, pbc, constraints
  (index-based reapplication after Packmol).

Acceptance
- Each builder produces a valid `Structure` with correct pbc/cell and provenance.
- Legacy nanotube+CO2 dataset reproducible through the new geometry path.

---

## Phase 3 — DFT labeling with full property set

**Goal:** label E/F/stress + dipole + polarizability.

Deliverables
- `calculators/dft.py` — QE and FHI-AIMS factories; SCF for E/F/stress; dipole
  output; **polarizability via DFPT** (AIMS dielectric; QE `ph.x`).
- Labeled results written to extxyz with level-of-theory provenance.
- Cost-aware labeling (polarizability flagged as the expensive task).

Acceptance
- A selected frame is labeled with all requested properties and lands in the
  dataset with provenance; QE and AIMS paths both verified on a tiny system.

---

## Phase 4 — Training + validation (multi-head)

**Goal:** train MACE on all properties and measure quality.

Deliverables
- `training/` — MACE fine-tune/train wrapper (`--foundation_model`), explicit
  multi-head config (energy/forces + dipole + polarizability), checkpoints,
  metrics. Pluggable model interface.
- `datasets/` health tooling: coverage maps, distributions, outliers,
  extrapolation grade.
- `validation/` — per-property parity + RMSE/MAE, learning curves, MD stability,
  EOS/phonons, **IR/Raman spectra reconstruction** vs DFT/experiment.

Acceptance
- A fine-tuned model trains on a seed set and produces a quality report covering
  every requested property, including a reconstructed Raman/IR spectrum.

---

## Phase 5 — Active-learning loop

**Goal:** close the loop to convergence.

Deliverables
- `active_learning/` — explore → select → label → retrain → converge, with
  thresholds drawn from `validation` (val force-RMSE, spectral error).
- Resume/idempotency across iterations.

Acceptance
- Starting from a seed set, the loop runs ≥2 iterations unattended, the dataset
  grows only with informative+diverse frames, and validation error decreases.

---

## Phase 6 — Orchestration

**Goal:** parallelism without touching the science.

Deliverables
- `orchestration/` — `local` engine hardened; QuACC adapter expressing the AL
  loop as a DAG with parallel explore/label fan-out. (Engine chosen by
  simplicity after the core is proven.)

Acceptance
- The same AL workflow runs under `local` and the chosen engine with identical
  results; explore and label stages fan out in parallel.

---

## Phase 7 — Polish & extras

- Curated public API + library usage docs and examples.
- Documentation site; tutorials reproducing the Raman use case.
- Deferred, architecture-friendly: **node-based workflow editor** emitting the
  serialized config DAG.
- Additional MLIP backends (MatterSim/Orb/SevenNet/CHGNet) via the model registry.

---

## Dependency graph (summary)

```
Phase 0 ─► Phase 1 ─► Phase 2 (geometry breadth)
                 └──► Phase 3 ─► Phase 4 ─► Phase 5 ─► Phase 6 ─► Phase 7
```
Phase 2 (geometry) and Phases 3–5 (label→train→loop) can progress in parallel
once Phase 1 (the spine) exists.
```
