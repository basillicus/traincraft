# Training (multi-head MACE)

Once the spine has produced a **labelled dataset**, the `[training]` stage turns
it into a trained interatomic potential. TrainCraft wraps MACE's
`mace_run_train` rather than reimplementing training: `core` prepares the
train/valid files and renders the command, the MACE/torch stack does the work
inside `traincraft-mlip.sif` (or a local `mace` env).

```toml
[training]
type = "mace"
name = "cu_finetune"
foundation_model = "medium"      # small | medium | large | mace-mp0 | <path.model>
strategy = "multihead"           # multihead | naive | scratch
heads = ["energy", "forces", "stress"]
e0s = "foundation"
pt_train_file = "mp"             # replay data for multihead fine-tuning
```

The trainer is a registry plugin (`@register("trainer", …)`), so adding a
backend (MatterSim/Orb/SevenNet/…) is one new file and one registry entry — the
same pattern as calculators and samplers.

## Where training sits in the pipeline

```
geometry → sample → select → label → dataset → TRAIN
```

`train` consumes the dataset artifact (or, if no `[dataset]` section, the
labelled frames) and writes a model tree under `runs/<name>/model/`:
`train.xyz`, `valid.xyz`, `model/<name>.model`, `checkpoints/`, `results/`,
`logs/`, and a `manifest.json` recording the foundation model, heads, E0s, the
exact rendered command, and the split sizes.

On HPC the stage runs as a **GPU step** (`--nv`, `traincraft-mlip.sif`); the
command is injected from `$TRAINCRAFT_MACE_TRAIN_COMMAND` so the science stays
container-agnostic (the same mechanism as the DFT labelers — see
[Run on HPC](../how-to/hpc.md)).

## Multi-head property targets

The property set in `heads` selects the MACE model type and loss. Energy/forces
(and stress on periodic data) are the core task; **dipole** and
**polarizability** are the heads that enable IR and Raman spectra reconstruction
(the original scientific driver — see [Calculators & DFT
Labeling](calculators.md)).

| `heads` include … | MACE `--model` | MACE `--loss` |
|---|---|---|
| energy / forces (+ stress) | *(foundation default)* | *(foundation default)* |
| dipole only | `AtomicDipolesMACE` | `dipole` |
| energy / forces + dipole | `EnergyDipolesMACE` | `energy_forces_dipole` |
| + polarizability | `AtomicDielectricMACE` | `dipole_polar` |

The dipole/polarizability model types are MACE's dielectric-property models (cf.
MACE-MDP and `mace-field`). They move faster than the energy/forces path, so the
exact `--model`/`--loss`/keys are **overridable** via `[training.extra]` — verify
them against your installed MACE version before a production polarizability run,
the same way the QE `ph.x` polarizability path is gated for labeling.

Labels are written with explicit reference keys (`REF_energy`, `REF_forces`,
`REF_stress`, `REF_dipole`, `REF_polarizability`) and passed to MACE via
`--energy_key`/`--forces_key`/… so a frame missing a label for one property is
simply skipped for that head rather than trained on zeros.

## Fine-tuning defaults (and why)

The defaults follow **Tompa, Varga-Umbrich, Batatia, Elena, Bernstein &
Csányi, _"Fine-tuning MLIP foundation models: strategies for accuracy and
transferability"_ ([arXiv:2606.12704](https://arxiv.org/abs/2606.12704))**. The
paper's headline finding is that *foundation-model quality, atomic-reference-energy
consistency, and stable optimization matter more than the choice of fine-tuning
method*. Concretely, TrainCraft defaults to:

| Setting | Default | Rationale (paper) |
|---|---|---|
| `e0s` | `"foundation"` | Reuse the foundation's isolated-atom energies. Averaging from data is **2–3× worse** on forces and can destabilize MD. |
| `strategy` | `"multihead"` | Multihead **replay** is the only method tested that consistently preserves out-of-distribution accuracy and the repulsive wall. Use `pt_train_file = "mp"` (or your own) for replay. |
| `weight_decay` | `0.0` | Essential when fine-tuning — non-zero decay pulls weights away from the pretrained solution. |
| `ema_decay` | `0.995` | A higher EMA decay (> 0.99) stabilizes fine-tuning. |
| `energy_weight`, `forces_weight` | `10`, `10` | Constant, energy-prioritised loss weights; single-stage (no force→energy schedule). |
| `lr` | `1e-3` | Naive fine-tuning tolerates `1e-3`–`1e-4`; multihead converges best nearer `1e-4`. |

> **Cost.** Multihead replay needs ~3–15× the compute of naive fine-tuning
> (lower learning rate + dual data batches). Use `strategy = "naive"` for a
> narrow, single-application model where forgetting is not a concern, and
> `strategy = "scratch"` (with `foundation_model` unset) to train from random
> initialization.

LoRA is discussed in the paper as an intermediate regularizer but is **not**
recommended as the primary method for MACE; if you need it, pass the relevant
flags through `[training.extra]`.

## Running it

```bash
# locally (needs the mace env: torch + mace-torch)
pixi run -e mace example-21

# or as part of a Slurm pipeline (train becomes a --nv GPU step)
traincraft submit examples/21_train_mace_finetune.toml --dry-run
```

See also: [Use Your Own MACE Model](../how-to/own-mace-model.md) for loading the
fine-tuned `.model` back in as an exploration calculator, and the
[Roadmap](../roadmap.md) for dataset health tooling and validation (parity,
learning curves, IR/Raman reconstruction), which build on this stage.
