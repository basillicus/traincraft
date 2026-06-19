# Python API

TrainCraft exposes a curated public API via `import traincraft`. Every stage
that the CLI runs is available as a standalone Python function.

---

## Top-level functions

```python
import traincraft as tc

# Config
cfg  = tc.load_config("my_run.toml")   # → TrainCraftConfig
cfg2 = tc.loads_config(toml_string)    # → TrainCraftConfig

# Pipeline stages (pure functions — same ones the CLI calls)
structure  = tc.build_geometry(cfg.geometry)
calc       = tc.make_calculator(cfg.calculator)
frames     = tc.run_sampling(structure, calc, job, cfg.sampling)
selected   = tc.run_funnel(frames, cfg.selection)
summary    = tc.run_pipeline(cfg)      # the whole pipeline

# Dataset IO
tc.write_frames("out.extxyz", frames)
frames = tc.read_frames("dataset.extxyz")
```

---

## `Structure`

::: traincraft.core.structure.Structure
    options:
      members:
        - hash
        - from_ase
        - to_ase
        - copy
        - fragments
        - n_fragments
        - set_fragments
        - to_pymatgen
        - from_pymatgen
        - to_rdkit
        - from_rdkit

---

## `Provenance`

::: traincraft.core.provenance.Provenance
    options:
      members:
        - to_dict
        - from_dict

---

## `Workspace` and `Job`

::: traincraft.core.workspace.Workspace
    options:
      members:
        - subdir
        - job

::: traincraft.core.workspace.Job
    options:
      members:
        - path

---

## Geometry

::: traincraft.geometry.build_geometry

---

## Converter

::: traincraft.core.converter.ase_to_pymatgen

::: traincraft.core.converter.pymatgen_to_ase

::: traincraft.core.converter.ase_to_rdkit

::: traincraft.core.converter.rdkit_to_ase

---

## Fragment helpers

::: traincraft.core.fragments.get_fragments

::: traincraft.core.fragments.set_fragments

::: traincraft.core.fragments.infer_fragments

::: traincraft.core.fragments.fragment_ids

::: traincraft.core.fragments.fragment_mask

---

## Registry

::: traincraft.core.registry.register

::: traincraft.core.registry.get

::: traincraft.core.registry.available

---

## Dataset

::: traincraft.datasets.dataset.Dataset
    options:
      members:
        - append
        - filter
        - write

::: traincraft.datasets.io.write_frames

::: traincraft.datasets.io.read_frames
