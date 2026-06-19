# Provenance

Every `Structure` in TrainCraft carries a `Provenance` record that answers
the question: *"where did this frame come from, and how was it modified?"*

Provenance is the foundation of reproducibility and dataset management.

---

## The `Provenance` dataclass

```python
@dataclass
class Provenance:
    origin: str                      # (1)
    source: str | None               # (2)
    transforms: list[str]            # (3)
    calculator: str | None           # (4)
    level_of_theory: dict            # (5)
    seed: int | None                 # (6)
    parents: list[str]               # (7)
    extra: dict                      # (8)
```

1. **`origin`** — the most important field. One of four increasing-cost tags:
   - `"generated"` — built from scratch (geometry only, no DFT)
   - `"ml_sampled"` — frame produced by an MLIP-driven MD/MC run
   - `"ml_labeled"` — frame labeled by an MLIP (not ground truth)
   - `"dft_labeled"` — frame labeled by DFT (ground truth)

2. **`source`** — a dotted string identifying what produced this frame.
   Examples: `"builder:crystal:Cu"`, `"source:smiles:CCO"`, `"source:url:..."`.

3. **`transforms`** — ordered list of transform names applied after the
   source/builder. Example: `["supercell:(2,2,2)", "strain:hydro=0.02"]`.

4. **`calculator`** — which calculator produced the `properties` dict.
   Set by the sampling engine after running a calculation.

5. **`level_of_theory`** — DFT-specific metadata (functional, basis set, …).
   Reserved for Phase 2 (DFT labeling).

6. **`seed`** — the RNG seed used for reproducible sampling.

7. **`parents`** — hashes of parent structures. Enables tracking lineage
   across active learning iterations.

8. **`extra`** — builder-specific metadata. SMILES builders store
   `{"smiles": "canonical_smiles", "fragment_smiles": {"0": "..."}}`.
   Defect builders store `{"defects": [...]}`.

---

## How provenance flows through the pipeline

```
build_geometry()
    → Structure(provenance=Provenance(
          origin="generated",
          source="builder:crystal:Cu",
          transforms=[],
      ))

apply_transform("supercell")
    → provenance.transforms.append("supercell:(2,2,2)")

apply_transform("strain")
    → provenance.transforms.append("strain:hydro=0.02")

run_sampling()
    → new Structure per frame, provenance.calculator = "emt"
      provenance.parents = [initial_structure.hash]
```

---

## The `origin` tag and dataset management

The `origin` tag is the key to separating cheap and expensive data:

```python
from traincraft import Dataset

ds = Dataset("runs/my_run/dataset")
all_frames = ds.frames()

# Only DFT-labeled frames — use for training
dft_frames = ds.filter(origin="dft_labeled")

# Generated + sampled frames — candidates for labeling
unlabeled = ds.filter(origin=["generated", "ml_sampled"])
```

This separation is critical in an active-learning loop: you never want to
mix DFT ground truth with ML-predicted labels in the same training batch.

---

## Provenance in extxyz files

Provenance is stored in the `tc_provenance` key of each frame's `info` dict:

```
2
Lattice="..." Properties=... tc_provenance={"origin": "generated", "source": "builder:crystal:Cu", ...}
Cu 0.000 0.000 0.000
Cu 1.805 1.805 0.000
```

You can read it with ASE:

```python
from ase.io import read
import json

frame = read("dataset.extxyz")
prov = json.loads(frame.info["tc_provenance"])
print(prov["source"])     # "builder:crystal:Cu"
print(prov["transforms"]) # ["supercell:(2,2,2)"]
```

---

## Serialisation

```python
d = prov.to_dict()         # → plain Python dict (JSON-serializable)
prov2 = Provenance.from_dict(d)  # → reconstruct from dict
```

Unknown keys in `from_dict` are silently ignored — older provenance records
remain readable as the schema evolves.
