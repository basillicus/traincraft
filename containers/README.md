# Containers — Apptainer images for Leonardo

TrainCraft ships as **three** Apptainer images rather than one monolith. The
split mirrors the plugin seam in the code and three separate concerns: hardware
target, rebuild cadence, and licensing. Full rationale: [`../DESIGN.md` §20](../DESIGN.md).

| Image | Def file | Target | Contents |
|-------|----------|--------|----------|
| `traincraft-core` | [`traincraft-core.def`](traincraft-core.def) | CPU (login + DCGP) | `traincraft` + CPU science stack (ASE, pymatgen, RDKit, Packmol, hiphive, tblite, dscribe). **The orchestrator.** |
| `traincraft-mlip` | [`traincraft-mlip.def`](traincraft-mlip.def) | GPU (Booster, A100/CUDA) | `traincraft` + PyTorch + CUDA + MACE. Sampling + training. |
| `traincraft-dft`  | [`traincraft-dft.def`](traincraft-dft.def)  | CPU (DCGP, Sapphire Rapids) | **FHI-aims** (MPI/MKL/ScaLAPACK) + `species_defaults`. **Private — licensed.** |

## Run model

`traincraft-core` is the orchestrator. It does **not** nest container execs.
Each stage that needs a different environment becomes a Slurm step that
`apptainer exec`s the right image:

```bash
# GPU sampling / training — traincraft runs inside the mlip image
srun --nv apptainer exec --bind "$SCRATCH" traincraft-mlip.sif \
     traincraft <stage> --job-dir "$JOB"

# DFT label — srun launches the bare FHI-aims MPI binary inside the dft image
srun apptainer exec --bind "$SCRATCH",<host-mpi-binds> traincraft-dft.sif aims.x
```

The container wrapper is **injected by the orchestration/executor layer**, never
hard-coded in a plugin. `dft.py` reads its run command from
`TRAINCRAFT_AIMS_COMMAND` (default `aims.x`); on Leonardo the executor sets it to
the `srun apptainer exec … aims.x` line above. So the `dft` image holds **only**
FHI-aims — no Python, no `traincraft`. `mlip` *does* embed `traincraft` because
sampling/training are Python stages run directly on the GPU.

## Building

Leonardo login nodes typically disallow rooted `apptainer build`. Options:

```bash
# A) fakeroot on a node that allows it
apptainer build --fakeroot traincraft-core.sif traincraft-core.def

# B) build off-cluster / in CI, then transfer the .sif
apptainer build traincraft-core.sif traincraft-core.def      # build host
rsync -avP traincraft-core.sif leonardo:$WORK/sif/           # transfer
```

`core` and `mlip` build from `pyproject.toml` + `pixi.lock` (reproducible).

### FHI-aims (private)

`traincraft-dft.def` needs the FHI-aims source, which is **license-gated and is
not committed to this repo**. Supply it at build time:

```bash
apptainer build --fakeroot \
  --build-arg AIMS_SRC=/path/to/fhi-aims.tar.gz \
  traincraft-dft.sif traincraft-dft.def
```

Do not push `traincraft-dft.sif` to any public registry.

## MPI on Leonardo (the hard part)

Multi-node FHI-aims uses the **hybrid model**, not a bundled MPI:

1. FHI-aims is compiled in the image against an MPI **ABI-compatible with the
   host** (Intel MPI / MPICH ABI), with MKL + ScaLAPACK, targeting Sapphire Rapids.
2. At runtime `srun` (host PMI) launches ranks and the host MPI / `libfabric` /
   UCX are **bind-mounted** so FHI-aims drives the real InfiniBand fabric.

Exact module versions and bind paths are filled in against the Leonardo modules
in use — see the `TODO(leonardo)` markers in `traincraft-dft.def`. Validate
single-node first, then multi-node, before any production labeling.
