# Run on HPC (Slurm + Apptainer)

TrainCraft dispatches the pipeline to **any Slurm cluster** as dependency-chained
jobs, each running in an Apptainer image. Nothing here is specific to one machine —
the account, partitions, modules, image paths, and bind mounts are all config.
CINECA Leonardo appears only as a worked example at the end.

> Architecture and rationale: [`DESIGN.md` §20](https://github.com/your-org/traincraft/blob/main/DESIGN.md).
> Definition files and build notes: [`containers/`](https://github.com/your-org/traincraft/tree/main/containers).

## The images

TrainCraft ships as a small set of images so each concern (hardware, rebuild
cadence, licensing) is isolated:

| Image | Target | Contents |
|-------|--------|----------|
| `traincraft-core` | CPU | `traincraft` + CPU science stack. **The orchestrator.** |
| `traincraft-mlip` | GPU | `traincraft` + PyTorch + CUDA + MACE. Sampling + training. |
| `traincraft-qe` | CPU | **Quantum ESPRESSO** (open source). DFT labeling. |
| `traincraft-dft` | CPU | **FHI-aims** (licensed). DFT labeling — needed for polarizability. |

You only build the DFT image you actually use. **QE is fully open source and needs
no license**, so the whole workflow can run end-to-end with open tooling only.

## Run model — the orchestrator dispatches Slurm steps

`traincraft-core` does not nest container execs. Each enabled stage becomes a
Slurm job that `apptainer exec`s the right image and runs `traincraft stage <name>`:

```bash
# GPU sampling/training — traincraft runs inside the mlip image
srun --nv apptainer exec --bind "$SCRATCH" traincraft-mlip.sif \
     traincraft stage sample config.toml

# DFT label — traincraft runs in core; the engine binary runs in its DFT image
#   under srun via the injected command (QE shown; FHI-aims is identical):
export TRAINCRAFT_PW_COMMAND="srun apptainer exec --bind $SCRATCH traincraft-qe.sif pw.x"
apptainer exec --bind "$SCRATCH" traincraft-core.sif traincraft stage label config.toml
```

`traincraft submit` generates and chains these for you (`--dependency=afterok`).
The DFT command is **injected from the environment**, never hard-coded in a plugin
(DESIGN §20.3), so the same config runs locally, with QE, with FHI-aims, on any
cluster.

## Building images

```bash
cd containers
apptainer build --fakeroot traincraft-core.sif traincraft-core.def
apptainer build --fakeroot traincraft-mlip.sif traincraft-mlip.def
apptainer build --fakeroot traincraft-qe.sif  traincraft-qe.def       # open source

# FHI-aims is licensed: supply the source at build time; never publish the .sif
apptainer build --fakeroot --build-arg AIMS_SRC=/path/to/fhi-aims.tar.gz \
  traincraft-dft.sif traincraft-dft.def
```

If your login nodes disallow rooted builds, build with `--fakeroot` where allowed
or build off-cluster and copy the `.sif` over.

## Configuring dispatch

```toml
[orchestration]
engine = "slurm"

[orchestration.slurm]
account = "<your-account>"        # your scheduler account/project
sif_dir = "$WORK/sif"             # where the .sif images live
modules = ["apptainer"]           # `module load` lines for your site
binds   = ["$SCRATCH", "$WORK"]   # filesystems to bind into the containers

[orchestration.slurm.stages.sample]
image     = "traincraft-mlip.sif"
partition = "<gpu-partition>"
gpus      = 1

[orchestration.slurm.stages.label]
partition = "<cpu-partition>"
nodes     = 2
ntasks    = 224
# pw_command / aims_command default to the qe / dft images; override if needed.
```

Submit:

```bash
traincraft submit config.toml --dry-run   # render + inspect the sbatch scripts
traincraft submit config.toml             # sbatch, dependency-chained
```

## MPI across nodes

Single-node DFT works out of the box. For **multi-node**, use the hybrid model:
the engine in the container is launched by the host `srun` (PMI), and the host
MPI / `libfabric` / UCX are bind-mounted so the DFT code drives the real
interconnect. A bundled MPI alone will not scale across nodes. Validate single-node
first, then multi-node.

## Worked example: CINECA Leonardo

Leonardo has a GPU **Booster** (A100) and a CPU **DCGP** partition. A typical
mapping:

```toml
[orchestration.slurm]
account = "EUHPC_xxxxxxx"
sif_dir = "$WORK/sif"
modules = ["apptainer"]
binds   = ["$SCRATCH", "$WORK"]

[orchestration.slurm.stages.sample]   # MACE on the Booster
image = "traincraft-mlip.sif"
partition = "boost_usr_prod"
gpus = 1

[orchestration.slurm.stages.label]    # DFT on DCGP
partition = "dcgp_usr_prod"
nodes = 2
ntasks = 224
```

The same config shape works on any other Slurm cluster — change the account,
partitions, modules, and binds to match your site. For Leonardo specifically, fill
the `TODO(leonardo)` markers in `containers/traincraft-dft.def` (compiler/MPI/
module versions) when building the FHI-aims image; the QE image needs no such tuning.
