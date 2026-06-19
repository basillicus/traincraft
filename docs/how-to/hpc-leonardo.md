# Run on HPC (Leonardo) with Apptainer

Production DFT labeling and MLIP training run on **CINECA Leonardo**. TrainCraft
ships as **three Apptainer images** rather than one monolith — the split mirrors
the plugin seam in the code and three separate concerns: hardware target, rebuild
cadence, and licensing.

> Full architecture and rationale: [`DESIGN.md` §20](https://github.com/your-org/traincraft/blob/main/DESIGN.md).
> Build/run details and the definition files live in
> [`containers/`](https://github.com/your-org/traincraft/tree/main/containers).

## The three images

| Image | Target | Contents |
|-------|--------|----------|
| `traincraft-core` | CPU (login + DCGP) | `traincraft` + CPU science stack (ASE, pymatgen, RDKit, Packmol, …). **The orchestrator.** |
| `traincraft-mlip` | GPU (Booster, A100) | `traincraft` + PyTorch + CUDA + MACE. Sampling + training. |
| `traincraft-dft` | CPU (DCGP) | **FHI-aims** (MPI/MKL/ScaLAPACK) + species defaults. **Private — licensed.** |

## Run model — the orchestrator dispatches Slurm steps

`traincraft-core` does **not** nest container execs. Each stage that needs a
different environment becomes a Slurm step that `apptainer exec`s the right image:

```bash
# GPU sampling / training — traincraft runs inside the mlip image
srun --nv apptainer exec --bind "$SCRATCH" traincraft-mlip.sif \
     traincraft <stage> --job-dir "$JOB"

# DFT label — srun launches the bare FHI-aims MPI binary inside the dft image
srun apptainer exec --bind "$SCRATCH",<host-mpi-binds> traincraft-dft.sif aims.x
```

The container/`srun` wrapper is injected by the orchestration layer, **never
hard-coded in a plugin**. The FHI-aims calculator reads its run command from
`$TRAINCRAFT_AIMS_COMMAND` (default `aims.x`); on Leonardo you set it to the
`srun apptainer exec … aims.x` line above and the plugin stays
container-agnostic. So `traincraft-dft` holds only FHI-aims — no Python.

## Building

Leonardo login nodes usually disallow rooted `apptainer build`. Either build with
`--fakeroot` where allowed, or build off-cluster / in CI and transfer the `.sif`:

```bash
apptainer build --fakeroot traincraft-core.sif traincraft-core.def
rsync -avP traincraft-core.sif leonardo:$WORK/sif/
```

`traincraft-core` and `traincraft-mlip` build reproducibly from `pyproject.toml`
+ `pixi.lock`.

### FHI-aims (private)

FHI-aims is **not redistributable**. Its source is *not* in this repo — supply it
at build time and never push the resulting `.sif` to a public registry:

```bash
apptainer build --fakeroot \
  --build-arg AIMS_SRC=/path/to/fhi-aims.tar.gz \
  traincraft-dft.sif traincraft-dft.def
```

## MPI across nodes (the hard part)

Multi-node FHI-aims uses the **hybrid model**, not a bundled MPI:

1. FHI-aims is compiled in the image against an MPI **ABI-compatible with the
   host** (Intel MPI / MPICH ABI), with MKL + ScaLAPACK, targeting Sapphire Rapids.
2. At runtime `srun` (host PMI) launches ranks and the host MPI / `libfabric` /
   UCX are **bind-mounted** so FHI-aims drives the real InfiniBand fabric.

Validate single-node first, then multi-node, before any production labeling.
Module versions and bind paths are filled in against the Leonardo modules in use —
see the `TODO(leonardo)` markers in `containers/traincraft-dft.def`.
