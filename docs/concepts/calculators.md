# Calculators & DFT Labeling

A calculator turns a structure into properties (energy, forces, …). TrainCraft
uses calculators in **two distinct roles**:

| Role | TOML | Purpose | Typical choice |
|---|---|---|---|
| **Exploration** | `[calculator]` | The cheap force engine that *drives sampling* (MD/MC) | `mace`, `emt`, `tblite` |
| **Labeling** | `[labeling.calculator]` | The expensive engine that *labels the selected frames* | `fhi_aims`, `qe` |

This separation is the whole point of the spine: explore cheaply, run the
selection funnel, then spend the expensive engine only on the survivors.

**TrainCraft is code-agnostic.** Any calculator that produces energies and forces
can label frames — FHI-aims, Quantum ESPRESSO, VASP, or any other QM code can be
plugged in as a single registered factory (an ASE-style calculator + a config
model + a registry entry; see [Write a Custom Builder](../how-to/custom-builder.md)
for the registry pattern — calculators work the same way). Two DFT engines ship
today:

- **`fhi_aims`** — the reference engine. Its academic license (from MS1P e.V.) is
  free for academic groups with a **voluntary donation**; you obtain the source by
  registration (so its container is built from your own copy — see
  [Run on HPC](../how-to/hpc.md)).
- **`qe` (Quantum ESPRESSO)** — fully open source, no license; runs the whole
  workflow with open tooling only.

Both can compute polarizability via DFPT (see below).

## Available calculators

| `type` | Kind | Properties | Extra deps |
|---|---|---|---|
| `emt` | ASE force field | energy, forces (stress for bulk) | none |
| `tblite` | semiempirical GFN-xTB | energy, forces, stress | `tblite` |
| `xtb` | semiempirical GFN-xTB | energy, forces | `xtb` |
| `mace` | MLIP (foundation / fine-tuned) | energy, forces, stress | `mace-torch`, `torch` |
| `fhi_aims` | DFT (FHI-aims) | energy, forces, stress, dipole, polarizability | FHI-aims binary |
| `qe` | DFT (Quantum ESPRESSO) | energy, forces, stress, dipole | QE binaries |

### `mace`
```toml
[calculator]
type   = "mace"
model  = "mace-mp0"      # or "mace-off23"
device = "cuda"          # "cpu" locally; "cuda" on the GPU image
# model_path = "my-finetuned.model"   # a local checkpoint instead of a foundation model
```

## DFT calculators

Both DFT calculators are ASE file-IO calculators. **Their run command is never
configured in the TOML** — it is injected from the environment so the same config
works locally and inside a container on HPC (see
[Run on HPC](../how-to/hpc.md)):

| Env var | Default | Used by |
|---|---|---|
| `TRAINCRAFT_AIMS_COMMAND` | `aims.x` | `fhi_aims` |
| `TRAINCRAFT_PW_COMMAND` | `pw.x` | `qe` |
| `TRAINCRAFT_AIMS_SPECIES_DIR` / `AIMS_SPECIES_DIR` | — | `fhi_aims` species defaults |
| `TRAINCRAFT_PW_PSEUDO_DIR` / `ESPRESSO_PSEUDO` | — | `qe` pseudopotentials |

On HPC the Slurm executor sets these for you, e.g.
`TRAINCRAFT_PW_COMMAND="srun apptainer exec … traincraft-qe.sif pw.x"` — you never
hard-code the command.

### `fhi_aims` (FHI-aims — the reference engine)
```toml
[labeling.calculator]
type             = "fhi_aims"
xc               = "pbe"
species_defaults = "tight"          # basis level (light/intermediate/tight/…)
kpts             = [4, 4, 4]        # Monkhorst-Pack grid (periodic systems)
relativistic     = "atomic_zora scalar"
properties       = ["dipole", "polarizability"]   # beyond E/F/stress
# extra = { ... }                   # any control.in keyword, passed through verbatim
```

**Polarizability via DFPT** (the IR/Raman driver) is computed in a single
`aims.x` run and selected automatically from the requested `properties` and the
system's periodicity: `dfpt = dielectric` for periodic cells, `dfpt =
polarizability` for molecules.

### `qe` (Quantum ESPRESSO — open source)
```toml
[labeling.calculator]
type    = "qe"
ecutwfc = 60.0
kpts    = [4, 4, 4]
pseudopotentials = { Si = "Si.pbe-n-kjpaw_psl.1.0.0.UPF" }
properties = ["dipole"]
```

Fully open source (conda-forge / the `traincraft-qe` image), so the whole workflow
can run with no licensed software. QE **can** compute polarizability/dielectric
response via DFPT (`ph.x` with `epsil`/`fpol`); because QE is a periodic code,
molecular polarizability is done in a vacuum supercell. That is a second binary
(`ph.x`) after the SCF, which TrainCraft's `qe` plugin does **not wire yet** — so
for polarizability today use `fhi_aims`, or extend `build_qe` with a `ph.x` step
(a tracked extension point, not a QE limitation).

## The labeling stage

When `[labeling]` is present, the pipeline labels the **selected** frames after
the funnel:

```toml
[labeling.calculator]
type = "fhi_aims"
# … settings as above …
```

This produces, under the run workspace:

```
labeled_dft/
  labeled.extxyz        # frames with E/F/stress(+dipole/pol) and origin=dft_labeled
  manifest.json         # level of theory, property set, frame count, wall time
  frame_0000/ …         # per-frame work dirs (so file-IO calcs don't collide)
```

Each labeled frame is tagged `origin="dft_labeled"` with its `level_of_theory`,
so the expensive data stays cleanly separable from the cheap, ML-generated points
(see [Provenance](provenance.md)). The dataset is built from these labeled frames.

> Tip: to *try the mechanics* with zero deps, set `[labeling.calculator] type =
> "emt"` — EMT stands in for DFT and produces E/F/stress. `examples/18` does this.
