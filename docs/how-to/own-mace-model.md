# Use Your Own MACE Model

Once you've fine-tuned MACE on your DFT data, plug the checkpoint into
TrainCraft for further active-learning iterations.

---

## In TOML

```toml
[calculator]
type          = "mace"
model         = "mace-mp0"           # used as a provenance label
model_path    = "/data/models/my_finetuned.model"   # (1)
device        = "cuda"               # or "cpu"
default_dtype = "float32"
```

1. **`model_path`** — absolute path to the `.model` file produced by
   `mace_run_train`. When this is set, it overrides the foundation model
   download.

---

## In Python

```python
from traincraft.calculators.potentials import make_mace_calc
from traincraft.config.models import MaceCalc

cfg = MaceCalc(
    model="mace-mp0",
    model_path="/data/models/my_finetuned.model",
    device="cuda",
)
calc = make_mace_calc(cfg)
```

---

## Fine-tuning with MACE

TrainCraft generates the dataset; fine-tuning is done with `mace_run_train`
(from the `mace-torch` package):

```bash
mace_run_train \
  --name="my_model" \
  --train_file="runs/my_run/dataset.extxyz" \
  --valid_fraction=0.1 \
  --foundation_model="mace-mp0" \
  --max_num_epochs=200 \
  --batch_size=4 \
  --device=cuda \
  --energy_key="tc_energy" \
  --forces_key="tc_forces"
```

The output model (`my_model.model`) can then be used in the next iteration of
active learning via `model_path`.

!!! note "Property key names"
    TrainCraft writes energy as `tc_energy` and forces as `tc_forces` in the
    extxyz `info`/`arrays` fields. Pass these to `mace_run_train` via
    `--energy_key` and `--forces_key`.

---

## Iterative active learning workflow

```
Iteration 1
  traincraft run seed_run.toml        ← generate + select
  DFT label selected frames
  mace_run_train → model_v1.model

Iteration 2
  update config: model_path = model_v1.model
  traincraft run iter2.toml           ← MACE-driven explore + select
  DFT label new frames
  mace_run_train → model_v2.model

...repeat until validation error converges
```
