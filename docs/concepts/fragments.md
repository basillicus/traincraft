# Fragment Identity

Fragment identity is TrainCraft's way of encoding **which atoms move together**
during Monte Carlo sampling. It's the key abstraction that makes rigid-body
MC and conformer-swap moves possible for molecules on surfaces.

---

## The `tc_fragment` array

Every `ase.Atoms` object in TrainCraft can carry a per-atom integer array
called `tc_fragment` (stored in `atoms.arrays["tc_fragment"]`).

The convention:

```
tc_fragment[i] == -1   →  atom i is framework / substrate (never moves in MC)
tc_fragment[i] >= 0   →  atom i belongs to mobile fragment with that ID
```

```
Cu(111) slab with CO adsorbate:

  O  ←  tc_fragment = 0   (mobile fragment 0)
  C  ←  tc_fragment = 0   (mobile fragment 0)
  ─────────────────────────
  Cu ←  tc_fragment = -1  (framework)
  Cu ←  tc_fragment = -1  (framework)
  Cu ←  tc_fragment = -1  (framework)
  ...
```

---

## Why this matters for MC sampling

The `monte_carlo` sampler reads the fragment array to know which atoms to move:

- **Translate move** — picks a random mobile fragment and displaces all its
  atoms by the same random vector (rigid body translation).
- **Rotate move** — picks a random mobile fragment and rotates all its atoms
  about the fragment's centroid.
- **Conformer move** — picks a random mobile fragment, generates a new RDKit
  conformer for its SMILES string, and proposes the new geometry as a Metropolis
  move.

Framework atoms (substrate) are **never** included in any MC move. This is
why CO can freely explore the Cu(111) surface without the slab distorting.

---

## How fragments are set

### Automatically by builders

Every builder that involves multiple components sets fragments automatically:

| Builder | Framework | Mobile |
|---|---|---|
| `surface_adsorbate` | Slab atoms (id=-1) | Adsorbate (id=0) |
| `surface_packing` | Slab atoms (id=-1) | Molecule *i* (id=i) |
| `slab` | All atoms (id=-1) | — (no mobile fragments) |
| `crystal` | — (no array set) | — |

### Manually in Python

```python
import numpy as np
from traincraft.core.fragments import set_fragments, FRAMEWORK

# Mark first 36 atoms as substrate, remaining 2 as one fragment
frag = np.full(38, FRAMEWORK, dtype=int)
frag[36:] = 0
set_fragments(atoms, frag)
```

### Inferring from connectivity (`infer_fragments`)

For reactive MD runs where bond topology changes, use `infer_fragments` to
recompute fragment IDs from the current bond graph:

```python
from traincraft.core.fragments import infer_fragments

# framework_mask: True for atoms that should stay fixed (e.g., slab)
framework_mask = np.zeros(len(atoms), dtype=bool)
framework_mask[:n_slab] = True

frag = infer_fragments(atoms, scale=1.2, framework_mask=framework_mask)
set_fragments(atoms, frag)
```

This uses ASE's `NeighborList` with covalent-radius cutoffs and
`scipy.sparse.csgraph.connected_components` to find bonded clusters.

---

## Persistence

Because `tc_fragment` lives in `atoms.arrays`, it survives:
- `atoms.copy()` — the array is copied
- extxyz read/write — ASE stores all arrays in extxyz
- `atoms.repeat()` — each image gets the same fragment ID
- Langevin MD — per-atom arrays are propagated automatically

!!! warning "Fragment IDs after `repeat()`"
    After a `repeat([2,2,1])`, fragment IDs are copied, so all images have the
    same IDs (0, 1, …). If you then add adsorbates manually, you'll need to
    re-assign fragment IDs to avoid collisions.

---

## Querying fragment information

```python
from traincraft.core.fragments import (
    get_fragments,    # → np.ndarray or None
    fragment_ids,     # → sorted list of mobile fragment IDs
    fragment_mask,    # → bool mask for one fragment
    set_fragments,    # → set the array
    infer_fragments,  # → compute from connectivity
    FRAMEWORK,        # == -1
    FRAGMENT_KEY,     # == "tc_fragment"
)

frag = get_fragments(atoms)
mobile = fragment_ids(atoms)        # e.g. [0, 1, 2]
mask0 = fragment_mask(atoms, fid=0) # bool array
```

Or via the `Structure` interface:

```python
s.fragments          # → np.ndarray or None
s.n_fragments        # → number of mobile fragments (excludes framework)
s.set_fragments(frag) # set or overwrite
```
