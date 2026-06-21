# Run on HPC (Slurm + Apptainer)

TrainCraft dispatches the pipeline to **any Slurm cluster** as dependency-chained
jobs, each running in an Apptainer image. Nothing here is specific to one machine —
the account, partitions, modules, image paths, and bind mounts are all config.
CINECA Leonardo appears only as a worked example at the end.

> Architecture and rationale: [`DESIGN.md` §20](https://github.com/basillicus/traincraft/blob/main/DESIGN.md).
> Definition files and build notes: [`containers/`](https://github.com/basillicus/traincraft/tree/main/containers).

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
apptainer exec --nv --bind "$SCRATCH" traincraft-mlip.sif \
     traincraft stage sample config.toml

# DFT label — traincraft runs in core; the engine binary runs in its DFT image
#   under srun via the injected command (QE shown; FHI-aims is identical):
export TRAINCRAFT_PW_COMMAND="srun --mpi=pmix apptainer exec --bind $SCRATCH traincraft-qe.sif pw.x"
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
runtime = "apptainer"            # "apptainer" (our images) or "native" (host binaries)
mpi     = "pmix"                  # Slurm MPI plugin — see "Picking the MPI plugin"
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

## Two knobs that make it portable: `runtime` and `mpi`

There is **no universal MPI setup** — the cluster decides. TrainCraft exposes this
as two independent config switches instead of baking in assumptions:

| Switch | Values | What it controls |
|--------|--------|------------------|
| `runtime` | `apptainer` \| `native` | Reach binaries via **our images**, or via **host binaries** already installed (site modules / conda / EasyBuild). `native` drops the container wrapper entirely. |
| `mpi` | `pmix` \| `cray_shasta` \| `pmi2` \| `none` | The Slurm MPI plugin used to launch the multi-node DFT step (`srun --mpi=<plugin>`). |

Both can be overridden per stage in `[orchestration.slurm.stages.*]`.

### Picking the MPI plugin

Run this on the target cluster — it is the ground truth:

```bash
srun --mpi=list
```

- **`pmix` present** (InfiniBand + Slurm, e.g. Leonardo) → `mpi = "pmix"`. Our
  images carry a self-contained OpenMPI+UCX+PMIx, so Slurm does the wire-up and no
  host MPI is needed.
- **No `pmix`, Cray/Slingshot** (e.g. LUMI shows only `cray_shasta` / `pmi2`) →
  `mpi = "cray_shasta"`. On Cray the path of least resistance is `runtime =
  "native"` with the site's `cray-mpich` and its FHI-aims/QE module, rather than
  fighting ABI translation inside a container.
- **Anything else** → `mpi = "pmi2"` is the portable fallback.

### When to use `runtime = "native"`

Use it when the cluster already provides tuned binaries (a site/EasyBuild FHI-aims
or QE, or your own conda/venv), or when bind-mounting containers is awkward (Cray).
TrainCraft then renders bare `srun --mpi=<plugin> aims.x` and runs `traincraft
stage …` directly — put the needed `module load` / `source activate` lines in
`modules` and `pre_commands`. See the LUMI example below.

## Worked example: CINECA Leonardo (Apptainer + PMIx)

Leonardo has a GPU **Booster** (A100) and a CPU **DCGP** partition, and
`srun --mpi=list` shows `pmix`. So: our images, PMIx launch. Full file:
[`examples/19_hpc_leonardo_label.toml`](https://github.com/basillicus/traincraft/blob/main/examples/19_hpc_leonardo_label.toml).

```toml
[orchestration.slurm]
account = "EUHPC_xxxxxxx"
runtime = "apptainer"
mpi     = "pmix"
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

Fill the `TODO(site)` markers in `containers/traincraft-dft.def` (target arch, MKL
link line) when building the FHI-aims image.

## Worked example: LUMI (native + cray_shasta)

LUMI is a Cray EX: Slingshot interconnect, and `srun --mpi=list` shows **no
pmix** (only `cray_shasta` / `pmi2`). The clean path is `runtime = "native"` using
the site's `cray-mpich` and FHI-aims/QE modules. Full file:
[`examples/20_hpc_lumi_native.toml`](https://github.com/basillicus/traincraft/blob/main/examples/20_hpc_lumi_native.toml).

```toml
[orchestration.slurm]
account = "project_465xxxxxx"
runtime = "native"                 # use host binaries, not our .sif images
mpi     = "cray_shasta"            # no pmix on Cray; this drives Slingshot
binds   = []
modules = ["LUMI/24.03", "partition/C", "cray-mpich"]
pre_commands = ["source $HOME/traincraft-venv/bin/activate"]  # traincraft on PATH

[orchestration.slurm.stages.sample]
partition = "standard-g"
gpus = 1

[orchestration.slurm.stages.label]
partition = "standard"
nodes = 2
ntasks = 256
pre_commands = ["module load fhi-aims/240507"]   # the site's FHI-aims build
```

This renders bare `srun --mpi=cray_shasta aims.x` with no container — the two knobs
(`runtime`, `mpi`) are the only things that change between Leonardo and LUMI; the
science config is identical.
