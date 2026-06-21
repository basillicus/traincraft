# Roadmap

TrainCraft is developed in phased, independently-testable chunks. The dataset +
selection layer is the spine — everything else builds on it.

---

## ✅ Phase 0 — Foundation

**Done.** A clean, installable, testable package with no globals.

- pydantic v2 config models (discriminated unions, `extra=forbid`)
- `Structure` (with content hash), registry, `Workspace`/`Job`, provenance
- Geometry: file/scratch sources; nanotube/molecule builders; vacuum/supercell/perturb transforms
- Calculators: `emt`, `tblite`/`xtb`, `mace` (MP0 + fine-tuned)
- Sampling: `md` (Langevin NVT), `rattle` (HiPhive)
- Selection funnel: physicality → dedup → diversity (FPS)
- Dataset: extxyz IO with provenance; hash-dedup
- CLI: `run` / `validate` / `new` / `plugins`
- 6 annotated examples, 19 tests, CI on GitHub Actions

---

## ✅ Phase 1, Chunk 1 — Molecules on Surfaces + Monte Carlo

**Done.** Fragment identity + surface adsorbate builders + Metropolis MC.

- `core/fragments.py` — per-atom `tc_fragment` array; `infer_fragments` for reactive runs
- `smiles` source — RDKit ETKDG + MMFF
- `surface_adsorbate` builder — single adsorbate on crystalline slab
- `surface_packing` builder — N-molecule coverage via Packmol
- `monte_carlo` sampler — translate/rotate/conformer-swap with Metropolis acceptance
- Examples 07–11 (CO on Cu, ethanol, butane)

---

## ✅ Phase 1, Chunk 2 — Mechanical Geometry Breadth

**Done.** All common bulk, surface, and 2D structure types.

- `core/converter.py` — ASE ↔ pymatgen ↔ RDKit bridge
- `url` source — download → ASE read
- `crystal` builder — bulk + supercell + vacancy/substitution/interstitial defects
- `slab` builder — named facet *or* arbitrary Miller indices; all-framework
- `layered` builder — graphene/hBN/MX₂; AA/AB stacking; twist (non-periodic moiré flake)
- Transforms: `strain` (hydrostatic/Voigt), `rotate`, `set_pbc`
- Examples 12–14, 36 additional tests

---

## ✅ Phase 1, Chunk 3 — Database Providers, Liquids, Intercalation, Constraints

**Done.** Remaining geometry breadth (except the `polymer` builder).

- Sources: `materials_project` (mp-api), `optimade` and `pubchem` (both dependency-free)
- `liquid` builder — Packmol multi-species box (explicit cell *or* density-driven)
- `intercalation` builder — guests per gallery of a planar layered host, with staging
- `constraints` transform — `FixAtoms` on the final structure (fixes the legacy
  index-misalignment bug after reordering builders)
- Examples 15–17, unit tests (network/Packmol paths skip when deps are absent)
- *Deferred:* `polymer` (PySoftK) — dependency is unreliable; wrapper to be verified
  against the live API rather than guessed

---

## 🟡 Phase 2 — DFT Labeling

Label selected frames with energy, forces, stress, dipole, and polarizability.

- ✅ `calculators/dft.py` — FHI-aims (`fhi_aims`) and Quantum ESPRESSO (`qe`) factories;
  FHI-aims polarizability via DFPT (`dielectric` periodic / `polarizability` molecular,
  auto-selected). Run command **injected from the environment** so the plugins stay
  container-agnostic. QE polarizability raises `NotImplementedError` (needs a `ph.x` run).
- ✅ Labeling stage (`[labeling]`): labels the *selected* frames, tags them
  `dft_labeled` with level of theory, writes `labeled_dft/` (`labeled.extxyz`,
  `manifest.json`, per-frame work dirs). `examples/18`.
- 🔜 Cost-aware labeling: polarizability flagged as the expensive task
- 🔜 Production runs on any Slurm cluster via the DFT container (or `runtime=native`)

---

## 🟡 Cross-cutting — Packaging & HPC Deployment (any Slurm cluster)

Run the real workflow on **any Slurm cluster** via **Apptainer** (or the cluster's
own binaries). Nothing is site-specific in the code. See
[`DESIGN.md` §20](https://github.com/basillicus/traincraft/blob/main/DESIGN.md) and
the [Run on HPC (Slurm + Apptainer)](how-to/hpc.md) guide.

- ✅ Architecture + four Apptainer `*.def` files (`containers/`): `traincraft-core`
  (CPU orchestrator), `traincraft-mlip` (GPU MACE), `traincraft-qe` (QE, open
  source), `traincraft-dft` (FHI-aims — private, licensed). DFT images are
  **compiled from source** (self-contained UCX+PMIx+OpenMPI).
- ✅ Resumable per-stage execution (`traincraft stage`) + a **portable Slurm
  executor** that renders dependency-chained sbatch scripts (`traincraft submit`,
  `[orchestration]` config) with two cluster-agnostic knobs: `runtime`
  (`apptainer` images | `native` host binaries) and `mpi`
  (`pmix`|`cray_shasta`|`pmi2`|`none`). `examples/19` (Leonardo, apptainer+pmix),
  `examples/20` (LUMI, native+cray_shasta).
- 🔜 Build + validate the images on a real cluster (single-node DFT, then multi-node)

---

## Phase 3 — Training + Validation *(in progress)*

Train a multi-head MACE model and measure quality end-to-end. Delivered in
chunks: training first (validation builds on it), then dataset health, then
validation.

**✅ Chunk 1 — Training (`training/`).** MACE fine-tune / train-from-scratch
wrapper over `mace_run_train`, as a pluggable `trainer` registry backend.
Multi-head property targets (energy/forces/stress + dipole + polarizability) map
onto MACE's model types and losses (`AtomicDipolesMACE` / `EnergyDipolesMACE` /
`AtomicDielectricMACE`). The `train` stage consumes the dataset and emits a model
tree (`model/<name>.model` + manifest); on HPC it runs as a GPU (`--nv`) step in
`traincraft-mlip.sif` with the command injected from the environment. Fine-tuning
defaults follow [Tompa et al. (arXiv:2606.12704)](https://arxiv.org/abs/2606.12704):
foundation-consistent E0s, multihead replay against forgetting, `weight_decay=0`,
high EMA, constant energy-prioritised loss weights. See
[Training](concepts/training.md); `examples/21`.

**🔜 Chunk 2 — Dataset health tooling (`datasets/`).** Composition/space/volume
coverage maps, per-element force distributions with outlier flags, extrapolation
grade, redundancy report.

**🔜 Chunk 3 — Validation (`validation/`).** Per-property parity + RMSE/MAE per
element, learning curves, NVE/MD stability, EOS/phonons, and **IR/Raman spectra**
reconstructed from MLIP-driven MD vs DFT/experiment.

---

## 🔜 Phase 4 — Active-Learning Loop

Close the loop: explore → select → label → retrain → converge.

- `selection/uncertainty.py` — committee/ensemble uncertainty selector
- `active_learning/` — full loop with resume/idempotency
- Convergence criteria: val force-RMSE + spectral error thresholds

---

## 🔜 Phase 5 — Orchestration

Parallel execution of the active-learning loop.

- Local engine hardened: threadpool for independent jobs
- QuACC adapter: explore + label stages as a parallel DAG
- Identical science, swappable engine

---

## 🔜 Phase 6 — Polish & Extras

- Full public API docs + library-usage tutorials (including Raman use case)
- Additional MLIP backends: MatterSim, Orb, SevenNet, CHGNet

### Agent workbench — a purpose-built web UI

A single browser app, served from the VM (WebGL rendering needs no X server),
that combines the conversational agent of
[Tutorial 11](tutorials/11-ai-agent.md) with tabbed views over one workflow.
It is a *front-end over the existing TOML spine* — the same configs the CLI and
agent already use, no parallel logic:

- **Chat** — drive the whole pipeline in plain language; the agent writes,
  validates and runs configs and reports back inline.
- **Geometry** — interactive 3D of the structure the agent just built
  (weas-widget / py3Dmol), with natural-language edits round-tripping to the
  agent.
- **Workflow** — the node-based editor: the pipeline DAG (geometry → sample →
  select → label → dataset → train) as nodes, edited visually and
  (de)serialised to/from the TOML the CLI runs.
- **Dataset** — interactive exploration of the generated dataset with
  [chemiscope](https://github.com/lab-cosmo/chemiscope) (structure–property
  maps linked to per-frame structures, descriptors, energies and forces).

Likely Streamlit/Gradio + stmol/py3Dmol + the chemiscope widget; the node editor
emits the serialised TOML DAG. Details TBD.

---

## Dependency graph

```mermaid
graph TD
    P0["✅ Phase 0<br/>Foundation"]
    P1A["✅ Phase 1 Ch.1<br/>Surfaces + MC"]
    P1B["✅ Phase 1 Ch.2<br/>Geometry breadth"]
    P1C["✅ Phase 1 Ch.3<br/>providers, liquid, intercalation"]
    P2["🟡 Phase 2<br/>DFT labeling"]
    HPC["🟡 Containers + HPC<br/>any Slurm cluster"]
    P3["🔜 Phase 3<br/>Training + Validation"]
    P4["🔜 Phase 4<br/>Active Learning"]
    P5["🔜 Phase 5<br/>Orchestration"]
    P6["🔜 Phase 6<br/>Polish"]

    P0 --> P1A --> P1B --> P1C
    P0 --> P2
    P2 --> P3 --> P4 --> P5 --> P6
    HPC -.-> P2
    HPC -.-> P3
```
