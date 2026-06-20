# Containers — Apptainer images for HPC

TrainCraft ships as a few small Apptainer images rather than one monolith. The
split mirrors the plugin seam in the code and three separate concerns: hardware
target, rebuild cadence, and licensing. It runs on **any Slurm cluster** —
nothing here is specific to one machine. Full rationale: [`../DESIGN.md` §20](../DESIGN.md).

| Image | Def file | Target | Contents |
|-------|----------|--------|----------|
| `traincraft-core` | [`traincraft-core.def`](traincraft-core.def) | CPU | `traincraft` + CPU science stack (ASE, pymatgen, RDKit, Packmol, hiphive, tblite, dscribe). **The orchestrator.** |
| `traincraft-mlip` | [`traincraft-mlip.def`](traincraft-mlip.def) | GPU (CUDA) | `traincraft` + PyTorch + CUDA + MACE. Sampling + training. |
| `traincraft-qe`   | [`traincraft-qe.def`](traincraft-qe.def)    | CPU | **Quantum ESPRESSO** (MPI). DFT labeling. **Open source — no license.** |
| `traincraft-dft`  | [`traincraft-dft.def`](traincraft-dft.def)  | CPU | **FHI-aims** (MPI/MKL/ScaLAPACK). DFT labeling + polarizability. **Private — licensed.** |

Build only the DFT image you use. The labeler is engine-agnostic: `qe` (open
source) or `fhi_aims` (licensed), selected by `[labeling.calculator].type`.

## Run model

`traincraft-core` is the orchestrator. It does **not** nest container execs.
Each stage that needs a different environment becomes a Slurm step that
`apptainer exec`s the right image:

```bash
# GPU sampling / training — traincraft runs inside the mlip image
srun --nv apptainer exec --bind "$SCRATCH" traincraft-mlip.sif \
     traincraft stage sample config.toml

# DFT label — srun launches the bare engine MPI binary inside its image
#   (QE shown; FHI-aims is identical with traincraft-dft.sif / aims.x)
srun apptainer exec --bind "$SCRATCH",<host-mpi-binds> traincraft-qe.sif pw.x
```

The container wrapper is **injected by the orchestration/executor layer**, never
hard-coded in a plugin. The DFT plugins read their run command from
`TRAINCRAFT_PW_COMMAND` / `TRAINCRAFT_AIMS_COMMAND` (default `pw.x` / `aims.x`);
the executor sets it to the `srun apptainer exec …` line above. So the DFT image
holds **only** the engine — no Python, no `traincraft`. `mlip` *does* embed
`traincraft` because sampling/training are Python stages run on the GPU.

## Building

Some login nodes disallow rooted `apptainer build`. Options:

```bash
# A) fakeroot on a node that allows it
apptainer build --fakeroot traincraft-core.sif traincraft-core.def

# B) build off-cluster / in CI, then transfer the .sif
apptainer build traincraft-core.sif traincraft-core.def      # build host
rsync -avP traincraft-core.sif mycluster:$WORK/sif/          # transfer
```

`core`, `mlip`, and `qe` build with no extra inputs (`core`/`mlip` from
`pyproject.toml` + `pixi.lock`; `qe` from conda-forge). **QE needs no license.**

### FHI-aims (private, optional)

Only needed if you label with `fhi_aims` (e.g. for polarizability).
`traincraft-dft.def` needs the FHI-aims source, which is **license-gated and is
not committed to this repo**. Supply it at build time:

```bash
apptainer build --fakeroot \
  --build-arg AIMS_SRC=/path/to/fhi-aims.tar.gz \
  traincraft-dft.sif traincraft-dft.def
```

Do not push `traincraft-dft.sif` to any public registry.

## MPI across nodes (the hard part)

Single-node DFT works out of the box. Multi-node uses the **hybrid model**, not a
bundled MPI:

1. The engine is built/installed against an MPI **ABI-compatible with the host**;
   for FHI-aims also MKL + ScaLAPACK tuned for the target CPU.
2. At runtime `srun` (host PMI) launches ranks and the host MPI / `libfabric` /
   UCX are **bind-mounted** so the engine drives the real interconnect.

Fill in your site's module versions and bind paths — for FHI-aims see the `TODO`
markers in `traincraft-dft.def` (QE from conda-forge needs no such tuning).
Validate single-node first, then multi-node, before any production labeling.
