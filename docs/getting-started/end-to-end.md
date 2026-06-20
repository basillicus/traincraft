# Full Workflow: From Scratch to HPC

This is the **line-by-line** walkthrough: clone the repo, build a dataset locally,
then run the *same* workflow on **any Slurm HPC cluster** — exploring with MACE on
a GPU partition and labeling the selected frames with DFT on a CPU partition.

The labeling engine is your choice: **Quantum ESPRESSO (`qe`) is fully open source
and is used as the default below**; FHI-aims is an alternative (licensed; needed
for polarizability). Follow it top to bottom; values specific to your cluster
(account, partitions, paths) are called out explicitly.

!!! note "What actually runs"
    Steps 1–5 run today with zero or minimal deps. Steps 6–9 need Apptainer and a
    Slurm allocation. With QE the whole path is open source; FHI-aims additionally
    needs your license/source and the `TODO` markers in
    `containers/traincraft-dft.def` filled in for your cluster's compiler/MPI.

---

## 0. Prerequisites

- `git`, a Linux shell.
- For the HPC part: an account on a Slurm cluster and `apptainer` available there.
  Labeling with QE needs nothing extra; FHI-aims needs its **source tarball**
  (licensed — not shipped with TrainCraft).

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

## 5. Switch labeling to real DFT

Replace the EMT stand-in with a real DFT engine. **Quantum ESPRESSO is open source
and the default** here:

```toml
[labeling.calculator]
type    = "qe"
ecutwfc = 60.0
kpts    = [4, 4, 4]
pseudopotentials = { Cu = "Cu.pbe-dn-kjpaw_psl.1.0.0.UPF" }
```

The **run command is injected from the environment**, so the same TOML works
locally and in a container. Locally, with QE installed:

```bash
export TRAINCRAFT_PW_COMMAND="mpirun -np 8 pw.x"
export ESPRESSO_PSEUDO=/path/to/pseudo
pixi run traincraft stage label my_workflow.toml
```

You never put the command in the TOML — on HPC the executor sets it for you.

??? note "Prefer FHI-aims? (e.g. you need polarizability)"
    Use `type = "fhi_aims"` instead; it adds `polarizability` via DFPT. It is
    licensed, so you build its image from your own source (step 6).
    ```toml
    [labeling.calculator]
    type             = "fhi_aims"
    xc               = "pbe"
    species_defaults = "tight"
    kpts             = [4, 4, 4]
    properties       = ["polarizability"]
    ```
    ```bash
    export TRAINCRAFT_AIMS_COMMAND="mpirun -np 8 aims.x"
    export AIMS_SPECIES_DIR=/path/to/species_defaults
    ```

---

## 6. Build the Apptainer containers

TrainCraft ships as a few images (see [Run on HPC](../how-to/hpc.md) for the
rationale). Build the ones you need from the `containers/` directory:

```bash
cd containers

# CPU orchestrator + science stack
apptainer build --fakeroot traincraft-core.sif traincraft-core.def

# GPU image for MACE sampling
apptainer build --fakeroot traincraft-mlip.sif traincraft-mlip.def

# DFT engine — Quantum ESPRESSO (open source, no license)
apptainer build --fakeroot traincraft-qe.sif traincraft-qe.def
```

??? note "Building the FHI-aims image instead"
    FHI-aims is licensed — supply your source tarball at build time and never
    publish the resulting `.sif`:
    ```bash
    apptainer build --fakeroot \
      --build-arg AIMS_SRC=/path/to/fhi-aims.tar.gz \
      traincraft-dft.sif traincraft-dft.def
    ```

!!! warning
    Some login nodes disallow rooted builds. Build with `--fakeroot` where
    permitted, or build off-cluster and `rsync` the `.sif` files over.

---

## 7. Stage images + config on the cluster

```bash
# from your build host
rsync -avP containers/*.sif        mycluster:$WORK/sif/
rsync -avP my_workflow.toml        mycluster:$WORK/runs/
```

Make sure the config and the run `outdir` live on a filesystem that you `--bind`
into the containers (e.g. `$WORK` / `$SCRATCH`).

---

## 8. Configure HPC dispatch

Add an `[orchestration]` section (full example: `examples/19_hpc_leonardo_label.toml`).
Everything below is your cluster's values — there is nothing site-specific in the
code:

```toml
[orchestration]
engine = "slurm"

[orchestration.slurm]
account = "<your-account>"       # your scheduler account/project
sif_dir = "$WORK/sif"
modules = ["apptainer"]          # `module load` names for your site
binds   = ["$SCRATCH", "$WORK"]

[orchestration.slurm.stages.sample]   # MACE on a GPU partition
image     = "traincraft-mlip.sif"
partition = "<gpu-partition>"
gpus      = 1
time      = "02:00:00"

[orchestration.slurm.stages.label]    # DFT on a CPU partition
partition = "<cpu-partition>"
nodes     = 2
ntasks    = 224
time      = "06:00:00"
```

You do **not** set the DFT command here — the executor generates it, pointing at
`traincraft-qe.sif` (or `traincraft-dft.sif` for FHI-aims) under `srun`.

---

## 9. Submit

From a cluster **login node** (where `sbatch` lives), first inspect the plan:

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

Each stage runs in its image: `sample` on the GPU partition (MACE), `label` runs
`traincraft` in the core image while the DFT engine (QE or FHI-aims) runs in its
own image under `srun`.

!!! tip "Validate small first"
    Before a production run, label a *single* small structure on one node to
    confirm the DFT build and MPI binds work, then scale to multiple nodes.

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
