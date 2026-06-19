# Manage Datasets

TrainCraft datasets are extended XYZ files (`extxyz`) with provenance metadata
baked into each frame. This page covers reading, writing, filtering, and merging them.

---

## Reading a dataset

```python
from traincraft import read_frames

frames = read_frames("runs/my_run/dataset.extxyz")
print(f"{len(frames)} frames")
print(frames[0].provenance.source)
print(frames[0].properties)
```

Each frame is a `Structure` object with `atoms`, `properties`, and `provenance`.

With ASE directly:

```python
from ase.io import read

atoms_list = read("dataset.extxyz", index=":")
frame = atoms_list[0]
print(frame.info["tc_provenance"])   # raw dict
print(frame.info["tc_hash"])         # content hash
```

---

## Writing frames

```python
from traincraft import write_frames

write_frames("my_dataset.extxyz", frames)
```

This appends to an existing file if it exists (deduplicating by hash).

---

## The `Dataset` class

`Dataset` provides hash-deduplication and `origin`-based filtering:

```python
from traincraft.datasets import Dataset

ds = Dataset("runs/my_run/dataset")   # (1)

# Append frames (deduped by hash)
ds.append(new_frames)

# Write to disk (returns path)
path = ds.write()   # → "runs/my_run/dataset.extxyz"

# Filter
dft_frames = ds.filter(origin="dft_labeled")
generated   = ds.filter(origin=["generated", "ml_sampled"])

print(f"Total: {len(ds)}, DFT-labeled: {len(dft_frames)}")
```

1. `Dataset` takes the path *without* the `.extxyz` extension — it manages
   the extension internally.

---

## Merging datasets from multiple runs

```python
from traincraft.datasets import Dataset
from traincraft import read_frames

# Collect all candidate frames
all_frames = []
for run_dir in ["runs/run1", "runs/run2", "runs/run3"]:
    all_frames.extend(read_frames(f"{run_dir}/dataset.extxyz"))

# Merge into one deduped dataset
merged = Dataset("combined/dataset")
merged.append(all_frames)
merged.write()
print(f"Merged: {len(merged)} unique frames")
```

---

## Inspecting dataset contents

```python
from traincraft import read_frames
from collections import Counter

frames = read_frames("dataset.extxyz")

# Origin breakdown
origins = Counter(f.provenance.origin for f in frames)
print(dict(origins))
# {'generated': 45, 'dft_labeled': 10}

# Source breakdown
sources = Counter(f.provenance.source for f in frames)
print(dict(sources))

# Composition breakdown
from ase.formula import Formula
formulas = Counter(
    str(Formula(f.atoms.get_chemical_formula())) for f in frames
)
print(dict(formulas))
```

---

## Converting to ASE for external tools

```python
from traincraft import read_frames
from ase.io import write

frames = read_frames("dataset.extxyz")
atoms_list = [f.to_ase() for f in frames]

# Export in any ASE format
write("dataset.db", atoms_list)          # ASE database
write("OUTCAR_style.xyz", atoms_list)    # plain xyz
```
