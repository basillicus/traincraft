# Full Workflow: From Scratch to HPC

This is the **line-by-line** walkthrough: clone the repo, build a dataset locally,
then run the *same* workflow on CINECA Leonardo — exploring with MACE on the GPU
Booster and labeling the selected frames with FHI-aims on the CPU partition.

Follow it top to bottom. Where a value is specific to your allocation (account,
partitions, paths, the FHI-aims source) it is called out explicitly.

!!! note "What actually runs"
    Everything in steps 1–5 runs today with zero or minimal deps. Steps 6–9 need
    Apptainer, a Leonardo allocation, and your **FHI-aims license/source**. The
    `TODO(leonardo)` markers in `containers/traincraft-dft.def` (module versions,
    compile flags, MPI binds) must be filled for your cluster.

---

## 0. Prerequisites

- `git`, a Linux shell.
- For the HPC part: a Leonardo account, `apptainer` on the cluster, and the
  **FHI-aims source tarball** (licensed — not shipped with TrainCraft).

---

## 1. Install

```bash
git clone https://github.com/your-org/traincraft.git
cd traincraft

# install pixi if you don't have it
curl -fsSL https://pixi.sh/install.sh | bash

pixi install            # core environment
pixi install -e dev     # + pytest/ruff (for the checks below)
```

Verify the install:

```bash
pixi run -e dev test            # the test suite should pass
pixi run traincraft plugins     # lists every registered source/builder/calculator/…
```

---

## 2. Your first dataset (local, zero deps)

```bash
pixi run traincraft run examples/01_cnt_emt_md.toml
```

Look at what it produced:

```bash
ls runs/cnt_emt_md/
# structures/initial.extxyz   candidates/   selected/   dataset.extxyz
```

The pipeline ran **geometry → sample → select → dataset**. `dataset.extxyz` is
your training set; every frame carries provenance.

---

## 3. Anatomy of a workflow — add DFT labeling

Open `examples/18_label_after_select.toml`. It is the full spine, one TOML:

```toml
[geometry.builder]      # what to build
type = "crystal"
name = "Cu"
crystalstructure = "fcc"
cubic = true
supercell = [2, 2, 2]

[calculator]            # CHEAP engine that drives sampling
type = "emt"

[sampling]              # explore
type = "md"
temperature = 600.0
steps = 200

[selection]             # the funnel (runs before any labeling)
budget = 5

[labeling.calculator]   # EXPENSIVE engine: label only the survivors
type = "emt"            # EMT stands in for DFT here (zero deps)

[dataset]
path = "dataset"
```

Run it:

```bash
pixi run example-18
ls runs/label_after_select/labeled_dft/    # labeled.extxyz  manifest.json  frame_0000/ …
```

The `[labeling]` section is the key addition: it labels the **selected** frames
and tags them `origin=dft_labeled`. See [Calculators & DFT Labeling](../concepts/calculators.md).

---

## 4. Run it stage by stage (this is what HPC does)

The pipeline is a chain of resumable stages, each reading the previous one's
artifact from the workspace. You can run them one at a time:

```bash
pixi run traincraft stage geometry examples/18_label_after_select.toml
pixi run traincraft stage sample   examples/18_label_after_select.toml
pixi run traincraft stage select   examples/18_label_after_select.toml
pixi run traincraft stage label    examples/18_label_after_select.toml
pixi run traincraft stage dataset  examples/18_label_after_select.toml
```

On HPC, each of these becomes its own Slurm job in the right container — but the
commands are identical.

---

## 5. Switch labeling to real DFT (FHI-aims)

Change the labeling section to FHI-aims and request the properties you need:

```toml
[labeling.calculator]
type             = "fhi_aims"
xc               = "pbe"
species_defaults = "tight"
kpts             = [4, 4, 4]
properties       = ["polarizability"]   # E/F/stress always; + DFPT polarizability
```

The FHI-aims **run command is injected from the environment**, so the same TOML
works locally and in a container. Locally, with FHI-aims installed:

```bash
export TRAINCRAFT_AIMS_COMMAND="mpirun -np 8 aims.x"
export AIMS_SPECIES_DIR=/path/to/species_defaults
pixi run traincraft stage label my_workflow.toml
```

You never put the command in the TOML — on Leonardo the executor sets it for you.

---

## 6. Build the three Apptainer containers

TrainCraft ships as three images (see [Run on HPC](../how-to/hpc-leonardo.md) for
the rationale). Build them from the `containers/` directory:

```bash
cd containers

# CPU orchestrator + science stack
apptainer build --fakeroot traincraft-core.sif traincraft-core.def

# GPU image for MACE (Booster)
apptainer build --fakeroot traincraft-mlip.sif traincraft-mlip.def

# DFT image — supply your FHI-aims source tarball (licensed; never committed/published)
apptainer build --fakeroot \
  --build-arg AIMS_SRC=/path/to/fhi-aims.tar.gz \
  traincraft-dft.sif traincraft-dft.def
```

!!! warning
    Leonardo login nodes may not allow rooted builds. Build with `--fakeroot`
    where permitted, or build off-cluster and `rsync` the `.sif` files over.
    **Never push `traincraft-dft.sif` to a public registry.**

---

## 7. Stage images + config on Leonardo

```bash
# from your build host
rsync -avP containers/*.sif        leonardo:$WORK/sif/
rsync -avP my_workflow.toml        leonardo:$WORK/runs/
```

Make sure the config and the run `outdir` live on a filesystem that you `--bind`
into the containers (e.g. `$WORK` / `$SCRATCH`).

---

## 8. Configure HPC dispatch

Add an `[orchestration]` section (full example: `examples/19_hpc_leonardo_label.toml`):

```toml
[orchestration]
engine = "slurm"

[orchestration.slurm]
account = "EUHPC_xxxxxxx"        # <-- your account
sif_dir = "$WORK/sif"
modules = ["apptainer"]
binds   = ["$SCRATCH", "$WORK"]

[orchestration.slurm.stages.sample]   # MACE on the GPU Booster
image     = "traincraft-mlip.sif"
partition = "boost_usr_prod"
gpus      = 1
time      = "02:00:00"

[orchestration.slurm.stages.label]    # FHI-aims on DCGP
partition = "dcgp_usr_prod"
nodes     = 2
ntasks    = 224
time      = "06:00:00"
```

You do **not** set `TRAINCRAFT_AIMS_COMMAND` here — the executor generates it,
pointing at `traincraft-dft.sif` under `srun`.

---

## 9. Submit

From a Leonardo **login node** (where `sbatch` lives), first inspect the plan:

```bash
apptainer exec $WORK/sif/traincraft-core.sif \
  traincraft submit $WORK/runs/my_workflow.toml --dry-run

cat $WORK/runs/<run-name>/slurm/label.sbatch     # check it looks right
```

Then submit for real — the stages are chained with `--dependency=afterok`:

```bash
apptainer exec $WORK/sif/traincraft-core.sif \
  traincraft submit $WORK/runs/my_workflow.toml
```

Monitor:

```bash
squeue --me
tail -f $WORK/runs/<run-name>/slurm/label-*.out
```

Each stage runs in its image: `sample` on the GPU Booster (MACE), `label` runs
`traincraft` in the core image while FHI-aims runs in the dft image under `srun`.

!!! tip "Validate small first"
    Before a production run, label a *single* small structure on one node to
    confirm the FHI-aims build and MPI binds work, then scale to multiple nodes.

---

## 10. Collect results

```bash
ls $WORK/runs/<run-name>/
#   dataset.extxyz                 -> your training set (origin=dft_labeled)
#   labeled_dft/manifest.json      -> level of theory, properties, counts
```

`dataset.extxyz` is ready for MACE training (Phase 3). You have gone from an
empty checkout to a DFT-labeled dataset produced on Leonardo — exploring cheaply,
labeling selectively, with every frame's provenance recorded.
