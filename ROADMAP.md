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

### ✅ Chunk 1 — Molecules-on-surfaces + Monte Carlo  *(done)*

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

### ✅ Chunk 2 — Mechanical geometry breadth  *(done)*

Unblocked by the fragment abstraction from Chunk 1:
- `core/converter.py` — ASE ↔ pymatgen ↔ RDKit bridge (periodic→`Structure`,
  molecular→`Molecule`; RDKit bond perception via `DetermineBonds`). Wired into
  `Structure.to_/from_pymatgen` and `.to_/from_rdkit`.
- Source: `url` (download → ASE read; `file://` works for local/offline use).
- Builders: `crystal` (bulk via `ase.build.bulk` + supercell + vacancy/
  substitution/interstitial defects with stable indexing); `slab` (standalone
  surface — named facet *or* Miller indices cleaved from bulk; bare slab is all
  framework); `layered` (graphene/hBN/MX2 stacking with interlayer spacing,
  AA/AB Bernal stacking, and twist — twisted stacks return a non-periodic moiré
  flake to avoid incommensurate-cell artefacts).
- Transforms: `strain` (hydrostatic or Voigt; scales atoms with the cell),
  `rotate` (named axis or vector, optional cell rotation), `set_pbc`.

Acceptance (met)
- Each builder/transform produces a valid `Structure` with correct pbc/cell and
  provenance; full unit coverage (`test_converter`, `test_crystal_builder`,
  `test_slab_builder`, `test_layered_builder`, `test_transforms`,
  `test_url_source`).
- `examples/12_bulk_vacancy_md.toml`, `13_slab_strain_md.toml`, and
  `14_graphene_bilayer_md.toml` run end-to-end (zero extra deps).

### ✅ Chunk 3 — Database providers + liquids + intercalation + constraints  *(done)*

- Sources: `materials_project` (mp-api + pymatgen), `optimade` (dependency-free
  JSON parse against any provider base URL), `pubchem` (dependency-free 3D SDF
  via PUG REST).
- Builders: `liquid` (Packmol multi-species box; explicit cell *or* density-driven
  cubic box; per-molecule fragment tagging) and `intercalation` (guests on an
  in-plane grid in each gallery of a planar layered host, with electrochemical
  `stage` control and `gallery_expansion`).
- Transforms: `constraints` (`FixAtoms` by indices / elements / `tc_fragment` /
  `below_z`, OR-ed). Selects on the *final* structure — the correct place to
  (re)apply constraints after a reordering builder (Packmol), fixing the legacy
  index-misalignment bug.

Acceptance (met)
- Each builder/source/transform yields a valid `Structure` with correct
  pbc/cell and provenance; unit coverage in `test_constraints`,
  `test_intercalation_builder`, `test_liquid_builder`, `test_providers`
  (network/Packmol paths skip cleanly when deps are absent).
- `examples/15_h_graphite_intercalation.toml` and `16_slab_constraints_md.toml`
  run end-to-end (zero deps); `17_water_box_liquid.toml` runs in the science env.

### Remaining Phase-1 item

- Builder: `polymer` (PySoftK). Deferred — the PySoftK dependency is not reliably
  resolvable (see `pyproject.toml`) and the wrapper must be verified against the
  live PySoftK API before shipping rather than guessed.

---

## Phase 2 — DFT labeling with full property set  *(mostly done)*

**Goal:** label E/F/stress + dipole + polarizability with **any QM code** that
yields energies and forces (FHI-aims = the reference/primary engine; QE = the
open-source option; VASP/others = a plugin each). Production runs on **any Slurm
cluster** (see "Packaging & HPC deployment" below).

Deliverables
- ✅ `calculators/dft.py` — QE and FHI-AIMS factories; SCF for E/F/stress; dipole
  output; **polarizability via DFPT** (AIMS `dielectric`/`polarizability`; QE
  needs `ph.x`, raises NotImplementedError for now). Command **injected from the
  environment** (`TRAINCRAFT_AIMS_COMMAND`/`TRAINCRAFT_PW_COMMAND`) so the plugin
  stays container-agnostic (DESIGN §20.3).
- ✅ Labeling stage (`labeling.py` + `[labeling]`): labels the *selected* frames,
  tags them `dft_labeled` with `level_of_theory`, writes `runs/<name>/labeled_dft/`
  (`labeled.extxyz`, `manifest.json`, per-frame work dirs). `examples/18`.
- 🔜 Cost-aware labeling (polarizability flagged as the expensive task; selection
  accounts for the cost).

Acceptance
- ✅ A selected frame is labeled and lands in `labeled_dft/` with provenance +
  manifest (verified with EMT standing in for DFT; mechanics covered by tests).
- 🔜 FHI-aims path verified end-to-end inside `traincraft-dft.sif` (a tiny
  molecule, single node) before scaling.

---

## Cross-cutting — Packaging & HPC deployment (any Slurm cluster)

**Goal:** run the real workflow on **any Slurm cluster** via Apptainer (or the
cluster's own binaries). Gates Phase 2 (DFT labeling) and Phase 3 (MACE training
at scale). Nothing is site-specific in the code — account/partitions/modules/binds
and the `runtime`/`mpi` knobs are all config. Full architecture in
[`DESIGN.md` §20](DESIGN.md); Leonardo and LUMI appear only as worked examples.

Four images, dispatched as Slurm steps by the `core` orchestrator:
- `traincraft-core.sif` (CPU) — package + geometry/selection/datasets/orchestration.
- `traincraft-mlip.sif` (GPU) — PyTorch + CUDA + MACE; sampling + training.
- `traincraft-qe.sif` (CPU) — **Quantum ESPRESSO** (open source); source-built.
- `traincraft-dft.sif` (CPU) — **FHI-aims** (MPI/MKL/ScaLAPACK); private build.

Deliverables
- ✅ `containers/` — four Apptainer `*.def` files + README. DFT images are
  **compiled from source** (UCX+PMIx+OpenMPI, self-contained); build via fakeroot /
  off-cluster, then transfer the `.sif`; FHI-aims license/source kept out of the repo.
- ✅ Resumable per-stage execution (`orchestration/stages.py` + `traincraft stage`)
  and a **portable Slurm executor** (`orchestration/slurm.py` + `[orchestration]`):
  renders dependency-chained sbatch scripts with two cluster-agnostic knobs —
  `runtime` (`apptainer` images | `native` host binaries) and `mpi`
  (`pmix`|`cray_shasta`|`pmi2`|`none`). `label` injects `srun --mpi=<plugin> …
  aims.x`/`pw.x`. `traincraft submit CONFIG [--dry-run]`; `examples/19` (Leonardo,
  apptainer+pmix), `examples/20` (LUMI, native+cray_shasta).
- ✅ Command-injection plumbing so `dft.py`/`mace` are container-agnostic.
- 🔜 DFT multi-node verified on a real cluster (PMIx launch + UCX transport on an
  IB cluster; `native`+cray-mpich on Cray). `TODO(site)` markers in the `.def`s for
  target-arch and the MKL link line.

Acceptance
- ✅ `submit --dry-run` renders the full chained pipeline (geometry→…→dataset) with
  correct images/resources/runtime/mpi/command-injection (covered by
  `test_slurm_executor`, incl. both an IB-PMIx and a Cray-native cluster).
- 🔜 End-to-end on a real cluster: GPU sampling (`mlip.sif`) + multi-node DFT label;
  results land in the dataset, identical in shape to a local run.

---

## Phase 3 — Training + validation (multi-head MACE)  *(in progress)*

**Goal:** train MACE on all properties and measure quality end-to-end. Delivered
in chunks (training first — validation builds on it — then health, then validation).

### ✅ Chunk 1 — Training (`training/`)  *(done)*

- `training/base.py` — `trainer` registry kind + `run_training`; `write_training_xyz`
  re-keys TrainCraft `tc_*` properties onto explicit MACE reference keys
  (`REF_energy`/`REF_forces`/`REF_stress`/`REF_dipole`/`REF_polarizability`) so
  labels are never silently dropped.
- `training/mace.py` — MACE backend: deterministic train/valid split, renders the
  full `mace_run_train` command, shells out, locates the `.model`, writes a
  `manifest.json`. **Container-agnostic**: the command is injected from
  `TRAINCRAFT_MACE_TRAIN_COMMAND` (DESIGN §20.3), defaulting to the bare console
  script; no torch import in `core`.
- **Multi-head mapping** from the requested `heads` to MACE model + loss:
  dipole → `AtomicDipolesMACE`/`EnergyDipolesMACE` (`dipole`/`energy_forces_dipole`),
  polarizability → `AtomicDielectricMACE` (`dipole_polar`). The dielectric model
  types track MACE-MDP / `mace-field`; `[training.extra]` overrides `--model`/`--loss`
  for version drift (the same gating philosophy as QE `ph.x` polarizability).
- **Fine-tuning defaults** follow Tompa, Varga-Umbrich, Batatia, Elena, Bernstein &
  Csányi, *"Fine-tuning MLIP foundation models"* (arXiv:2606.12704): `e0s="foundation"`
  (averaging is 2–3× worse), `strategy="multihead"` replay against catastrophic
  forgetting, `weight_decay=0`, `ema_decay=0.995`, constant λ_E = λ_F = 10.
- `train` stage (`stages.py` + `local.py`): consumes the dataset (or labelled
  frames) → `runs/<name>/model/`. Slurm executor maps it to a GPU (`--nv`)
  `traincraft-mlip.sif` step (DESIGN §20.6). `examples/21`.

Acceptance (met)
- ✅ `train` renders a correct `mace_run_train` invocation (flags, reference keys,
  multi-head model/loss mapping, paper-informed defaults, command injection,
  `extra` passthrough) and runs it as a stage; covered by `test_training` (subprocess
  mocked — MACE training is GPU/hours, not a CI op), mirroring `test_dft`/`test_slurm_executor`.
- 🔜 A real fine-tune on a seed set inside `traincraft-mlip.sif` produces a usable
  `.model` (verified on a GPU node).

### 🔜 Chunk 2 — Dataset health tooling (`datasets/`)

Composition/space/volume coverage maps; energy and per-element force distributions
with outlier flags; extrapolation grade; redundancy report.

### 🔜 Chunk 3 — Validation (`validation/`)

Per-property parity + RMSE/MAE (per element), learning curves, NVE/MD stability,
EOS/phonons, **IR and Raman spectra reconstruction** from MLIP-driven MD vs
DFT/experiment.

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
- Deferred, architecture-friendly: **agent workbench** — a purpose-built web UI
  (served from the VM, no X server needed) that pairs the conversational agent
  of `docs/tutorials/11-ai-agent.md` with tabbed views over one workflow. It is
  a front-end over the existing TOML spine (same configs the CLI/agent run — no
  parallel logic):
  - **Chat** — drive the whole pipeline in plain language; the agent writes,
    validates and runs configs and reports back inline.
  - **Geometry** — interactive WebGL 3D of the structure just built
    (weas-widget / py3Dmol), with natural-language edits round-tripping to the
    agent.
  - **Workflow** — the **node-based editor**: the pipeline DAG (geometry →
    sample → select → label → dataset → train) as nodes, edited visually and
    (de)serialised to/from the TOML the CLI runs.
  - **Dataset** — interactive exploration with **chemiscope**
    (https://github.com/lab-cosmo/chemiscope): structure–property maps linked to
    per-frame structures, descriptors, energies and forces.

  Likely Streamlit/Gradio + stmol/py3Dmol + the chemiscope widget. Details TBD.

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
