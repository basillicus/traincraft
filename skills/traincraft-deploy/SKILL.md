---
name: traincraft-deploy
description: >-
  Bring up a working TrainCraft pipeline on an HPC cluster — from a laptop driving
  any cluster over SSH, or running on the cluster itself. Probe the site, pick the
  least-effort runtime/MPI, provision it, write a reusable cluster profile, and
  prove it with end-to-end smoke jobs (DFT label on CPU + MACE train/validate on GPU).
---

# TrainCraft Deploy Skill

**Repository:** https://github.com/basillicus/traincraft

You stand up [TrainCraft](https://github.com/basillicus/traincraft) on an HPC
cluster so the pipeline (geometry → sample → select → **label** → dataset →
**train**) runs as dependency-chained Slurm jobs. You usually run on the user's
**laptop** and drive a remote cluster over SSH; you may also run **on the cluster**
itself. You have `bash`, so you operate through `ssh`/`rsync` and the TrainCraft
CLI — no special tools are needed.

This skill is about **bring-up and configuration**, not science. For building
configs and running workflows, use the `traincraft` skill.

## ⛔ Mandatory rules (read before doing ANYTHING)

1. **ALWAYS read the docs first.** Before any probe, build, or submit, you **must**
   have read `docs/getting-started/installation.md`, `docs/how-to/hpc.md`, and
   `DESIGN.md §20`. They are the single source of truth (runtime/MPI knobs, the
   image model, the env-injected DFT command). Point at them; never restate facts
   from memory. No command before the docs.
2. **Confirm before every remote, expensive, or outward action.** Show the exact
   command, explain what it does, and **ask** before running it — every `ssh`
   exec, image build, `rsync`/`scp` transfer, and `sbatch`. Default to dry-run.
   Only skip the per-action prompt if the user explicitly grants looser permission
   for the session.
3. **Never act on credentials — but ALWAYS guide fully.** You must not generate
   SSH keys, fetch tokens, accept licenses, or log in on the user's behalf. You
   **must** still walk them through exactly how to set each up (commands, links,
   what to click). "I can't touch credentials" is never an excuse to leave them
   stuck — pair every refusal-to-act with a complete how-to. See *Credentials*.
4. **Law of minimal effort — probe before you provision.** Discover what the site
   already provides and walk the runtime ladder (below). Do not blindly build
   containers when a site module already gives you the binary.
5. **Done = both smoke jobs pass.** Nothing is "deployed" until a real DFT label
   (CPU path) **and** a real MACE train+validate (GPU path) have run end-to-end
   through the actual scheduler/runtime/MPI. Stamp the profile only then.
6. **When in doubt, ask.** Missing account, unclear partition, ambiguous module —
   ask rather than guess.

## 🔌 Execution context (detect it, don't ask twice)

- **Remote-from-laptop** (common): every cluster action is `ssh <alias> '<cmd>'`
  and files move with `rsync`/`scp`. Use a `~/.ssh/config` host alias. Offer to
  set up an `ssh -M` **ControlMaster** socket (with confirmation) so one
  authenticated login is reused for all subsequent non-interactive commands —
  this is what makes laptop-driven probing work behind OTP/MFA or a jump host. The
  user can decline and point you at a session they opened themselves.
- **Local-on-cluster**: identical steps, just drop the `ssh <alias>` wrapper.

Detect which you're in (is the scheduler local?) and adapt; don't re-ask.

## 🔭 Phase A — Connect & probe (read-only; safe to batch after confirming)

```bash
ssh <alias> true                                   # auth/reachability works at all
ssh <alias> 'command -v sbatch sinfo squeue'       # is this a Slurm cluster?
ssh <alias> 'srun --mpi=list'                      # GROUND TRUTH for the mpi knob
ssh <alias> "sinfo -o '%P %G %l'"                  # partitions, GPUs (%G), time limits
ssh <alias> 'module -t avail 2>&1 | grep -iE "aims|quantum|espresso|mpich|openmpi|cuda|mkl|oneapi|intel"'
ssh <alias> 'command -v apptainer singularity'     # container runtime present?
ssh <alias> 'grep -c "^$USER:" /etc/subuid /etc/subgid; cat /proc/sys/user/max_user_namespaces 2>/dev/null'
#   ^ can we --fakeroot build HERE? login nodes almost always say no → build on the laptop instead
ssh <alias> 'lscpu | grep "Model name"'            # only if the user later wants to optimize the DFT build (optional)
ssh <alias> 'command -v traincraft pixi uv'        # anything already installed?
ssh <alias> 'echo "$WORK $SCRATCH"; quota 2>/dev/null'  # filesystems for binds/sif_dir
```

Summarise findings back to the user before deciding anything.

## 🪜 Phase B — Decide `runtime` + `mpi`

**`mpi`** comes straight from `srun --mpi=list` (see `hpc.md`): `pmix` (InfiniBand+
Slurm, e.g. Leonardo) → `cray_shasta` (Cray/Slingshot, e.g. LUMI — no pmix) →
`pmi2` as the portable fallback.

**`runtime`** — walk the ladder, stop at the first that works (least effort first):

1. **Native + site binaries** — a site module for QE/FHI-aims **and** an MPI module
   exist → `runtime = "native"`. Least effort; let the site own the heavy builds.
2. **Apptainer images** — build where `--fakeroot` is allowed (usually the laptop)
   and `rsync` the `.sif` to `sif_dir`, or copy a prebuilt one → `runtime = "apptainer"`.
3. **Compile on the HPC** — QE (open source) or FHI-aims (cloned from the user's
   GitLab) built against the site MPI; point `runtime = "native"` at the binary.
4. **Plain conda env** — *last resort, discouraged.* Only if 1–3 are impossible,
   and say so explicitly. (TrainCraft itself always goes in a pixi/uv venv — never
   a global install; never conda if avoidable.)

## 🧰 Phase C — Provision (only what Phase B chose)

- **traincraft** must exist wherever a stage runs: a `traincraft-core`/`-mlip`
  image (apptainer) or a pixi/uv venv on the cluster (native). Confirm
  `traincraft plugins` runs there.
- **Native:** confirm the `module load` lines; stand up the venv if absent.
- **FHI-aims source (needed by both the image and a native compile).** If the user
  confirms they have Git access, you may clone it; otherwise stop and guide them
  (see *Credentials*). Never act on the credentials themselves.
  ```bash
  git clone git@aims-git.rz-berlin.mpg.de:aims/FHIaims.git
  ```
  On an authentication error, **stop immediately** and fix credentials first — do
  not retry blindly.
- **Apptainer (DFT image): build on the laptop by default, then transfer.** Image
  builds need `--fakeroot`, which login nodes almost always forbid — so an agent
  that tries to build on the cluster will just fail. **Default: build the `.sif` on
  the laptop and `rsync` it to `sif_dir`** (always works). Only build on the cluster
  if the Phase A subuid/userns probe shows `--fakeroot` actually works there. **Never
  publish the DFT `.sif`** (licensed). The source goes in as a **tarball** build-arg
  (the def untars with `--strip-components=1`, so it needs one top-level dir):
  - From a clone, make a clean archive, then build:
    ```bash
    git -C ./FHIaims archive --prefix=aims/ --format=tar.gz -o /tmp/fhi-aims.tar.gz HEAD
    apptainer build --fakeroot --build-arg AIMS_SRC=/tmp/fhi-aims.tar.gz \
      traincraft-dft.sif containers/traincraft-dft.def
    ```
    (Submodules? `git archive` misses them — fall back to
    `tar czf /tmp/fhi-aims.tar.gz -C <parent-dir> FHIaims`.)
  - From an existing tarball: `--build-arg AIMS_SRC=/path/to/fhi-aims.tar.gz`.
  - The def **builds and runs with portable defaults — you do not need to edit it**
    to get a working image. The `TODO(site)` markers (target-arch flag, MKL link
    line) are **optional performance tuning only**; skip them unless the user
    explicitly asks to optimize, and even then only *propose* values to confirm —
    never silently edit or guess the arch.
- **Compile on the HPC (native):** QE from source; FHI-aims (the same clone) built
  against the site MPI, then point `runtime = "native"` at the resulting binary.

## 📇 Phase D — Write the cluster profile (the one-line switch)

Save the cluster's `[orchestration.slurm]` block **once** as a named profile, so a
workflow targets it by changing one line. Profiles live in
`~/.traincraft/clusters/<name>.toml` (dir overridable via `TRAINCRAFT_CLUSTERS_DIR`)
and hold exactly the keys that go under `[orchestration.slurm]`:

```toml
# ~/.traincraft/clusters/leonardo.toml
account = "EUHPC_xxxxxxx"
runtime = "apptainer"          # from Phase B
mpi     = "pmix"               # from Phase B
sif_dir = "$WORK/sif"
modules = ["apptainer"]
binds   = ["$SCRATCH", "$WORK"]

[stages.sample]                # MACE on the GPU partition
partition = "boost_usr_prod"
gpus = 1

[stages.label]                 # DFT on the CPU partition
partition = "dcgp_usr_prod"
nodes = 2
ntasks = 224
```

A workflow then selects it with one line; inline keys still **override** the
profile (use that when testing):

```toml
[orchestration]
engine = "slurm"

[orchestration.slurm]
profile = "leonardo"           # ← swap to target a different cluster
```

Track progress in a **status file** kept *beside* the profile,
`~/.traincraft/clusters/<name>.deploy.json` — deliberately **not** inside a
`<name>/` directory, so the config loader (which resolves the flat `<name>.toml`
and only reads `*.toml`) never sees it. Write it after each successful phase and
read it on start to **resume** instead of redoing:

```json
{
  "connected":   "2026-06-23T10:00:00Z",
  "runtime":     "apptainer",
  "provisioned": "2026-06-23T10:40:00Z",
  "dry_run_ok":  "2026-06-23T10:42:00Z",
  "smoke_label": "2026-06-23T11:05:00Z",
  "smoke_mace":  null
}
```

Both `smoke_label` and `smoke_mace` must be non-null before the deployment counts
as done. See the *Cluster profiles* section of `docs/how-to/hpc.md`.

## ✅ Phase E — Validate (render only, no jobs)

```bash
ssh <alias> 'cd <project> && traincraft submit config.toml --dry-run'
# then read back runs/<name>/slurm/*.sbatch — confirm the images/partitions/mpi/
# binds are real and the injected TRAINCRAFT_*_COMMAND lines look right.
```

## 🔥 Phase F — Smoke jobs (MANDATORY — this is what "done" means)

Both must pass. Run them for real (with confirmation), watch via `squeue`, inspect
outputs.

- **F1 — DFT label, CPU path.** A deliberately tiny label (≈2 atoms, QE, minimal
  basis, 1 SCF). Proves scheduler + runtime + `srun --mpi=<plugin>` + the DFT
  binary actually fire together.
- **F2 — MACE train+validate, GPU path.** On the GPU partition with `--nv`, run a
  short `train` of a tiny MACE on the smoke data (a handful of epochs), then
  validate: load the model and run a one-point inference / few-step MD and check
  the forces are finite/sane. Proves the GPU was allocated, CUDA is visible inside
  `traincraft-mlip`, torch sees the device, and MACE both trains **and** infers.

Stamp the profile `smoke ✓` only when **both** succeed.

## 🔐 Credentials (guide fully, act never)

For each item: **ask** "is this already set up?" If yes, move on. If no, give the
complete how-to and stop — do not perform it yourself.

- **SSH access:** `ssh-keygen -t ed25519`, `ssh-copy-id <user>@<host>`, and a
  `~/.ssh/config` alias (`Host`, `HostName`, `User`, `IdentityFile`,
  `ControlMaster auto`/`ControlPath`/`ControlPersist`). Explain how the site's
  MFA/OTP and any jump host (`ProxyJump`) fit in.
- **Allocation/account:** point them to the site's project/account portal; you
  need the account string for `[orchestration.slurm] account`.
- **FHI-aims (licensed):** they must hold a license and have their Git/SSH
  credentials configured. Once confirmed, `git clone
  git@aims-git.rz-berlin.mpg.de:aims/FHIaims.git` "just works" and you can
  build/compile it. If the clone fails on authentication, stop and guide them to
  fix their SSH keys / Git config — you never fetch the license or handle their key.
- **Never** write secrets into anything that gets committed, and never echo them.

## 🧱 Human-only blockers to surface early

Rooted/`--fakeroot` build policy on login nodes, project allocation & quotas, the
FHI-aims license, and SSH key/MFA setup. For each: confirm it's handled or guide
the user — then continue.

## Command cheat-sheet

```
ssh <alias> 'srun --mpi=list'                       # pick the mpi knob
ssh <alias> "sinfo -o '%P %G %l'"                   # partitions / GPUs / limits
ssh <alias> 'command -v sbatch apptainer traincraft'# what's available
traincraft submit config.toml --dry-run             # render sbatch scripts only
traincraft submit config.toml                        # submit, dependency-chained
# profile lives at ~/.traincraft/clusters/<name>.toml; select with profile = "<name>"
```
