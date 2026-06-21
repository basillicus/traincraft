# Train & Use Your Own MACE Model

TrainCraft both **trains** a MACE model (the `[training]` stage) and lets you
**plug a trained model back in** as the exploration engine for the next
active-learning iteration. This page is the quick task reference; for a guided,
illustrated walkthrough see [Tutorial 10 · Training](../tutorials/10-training.md).

---

## Train a model (the `[training]` stage)

Add a `[training]` section to any config that produces a dataset. The stage
splits the data, re-keys the labels for MACE, and runs `mace_run_train` for you:

```toml
[training]
type             = "mace"
name             = "my_model"
foundation_model = "medium"          # small | medium | large | mace-mp0 | <path.model>
strategy         = "multihead"       # multihead replay (robust) | naive | scratch
heads            = ["energy", "forces", "stress"]
e0s              = "foundation"
pt_train_file    = "mp"              # replay data for multihead fine-tuning
device           = "cuda"            # "cpu" for tiny demos
```

```bash
pixi run -e mace traincraft run my_run.toml      # whole pipeline, incl. training
# or just the training stage on an existing dataset:
pixi run -e mace traincraft stage train my_run.toml
```

The trained potential lands at `runs/<name>/model/<name>.model`, alongside a
`manifest.json` recording exactly how it was trained. See the
[Config Schema](../reference/config.md#training) for every `[training]` field and
the [Training concept page](../concepts/training.md) for the defaults' rationale.

!!! tip "Multi-head for IR / Raman"
    Add `dipole` (IR) and `polarizability` (Raman) to `heads` — provided your
    dataset was labelled with those properties. This selects MACE's dielectric
    model types automatically.

---

## Use a trained model as the exploration engine

Once you have a `.model`, point the `[calculator]` at it so the *next* round of
sampling is driven by your potential instead of the foundation:

```toml
[calculator]
type          = "mace"
model         = "mace-mp0"                       # (1) provenance label
model_path    = "runs/my_run/model/my_model.model"   # (2)
device        = "cuda"                           # or "cpu"
default_dtype = "float32"
```

1. **`model`** — kept as a provenance label when `model_path` is set.
2. **`model_path`** — absolute or relative path to the `.model` from
   `mace_run_train`. When set, it overrides the foundation-model download.

In Python:

```python
import traincraft as tc
from traincraft.config.models import MaceCalc

calc = tc.make_calculator(MaceCalc(
    model_path="runs/my_run/model/my_model.model",
    device="cuda",
))
```

---

## Iterative active-learning workflow

```
Iteration 1
  traincraft run seed_run.toml         ← generate → select → label → dataset → TRAIN
  → runs/seed_run/model/model_v1.model

Iteration 2
  set [calculator].model_path = .../model_v1.model   ← explore with your model
  traincraft run iter2.toml            ← MACE-driven explore → select → label → TRAIN
  → model_v2.model

…repeat until validation error converges
```

Phase 4 automates this loop (committee uncertainty selection + convergence
criteria); today you drive the iterations by hand. See the
[Roadmap](../roadmap.md).

---

## Under the hood

The `[training]` stage is a thin, reproducible front-end over `mace_run_train`.
It writes the train/valid splits with explicit reference keys (`REF_energy`,
`REF_forces`, `REF_stress`, `REF_dipole`, `REF_polarizability`) and passes them
via `--energy_key`/`--forces_key`/…, so labels are read correctly and a frame
missing one property is skipped for that head rather than trained on zeros. The
exact rendered command is saved in `manifest.json` (and viewable with a
`dry_run=True` call — see [Tutorial 10](../tutorials/10-training.md#inspect-before-you-train)).

Because the command is **injected from `TRAINCRAFT_MACE_TRAIN_COMMAND`**, the same
config trains locally or inside `traincraft-mlip.sif` on a GPU node with no change
— the HPC executor sets the variable for you (DESIGN §20.3).
