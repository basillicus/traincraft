"""Fragment identity: which atoms move together.

Convention (project-wide):
  tc_fragment[i] == -1   -> atom i is framework / substrate (immobile)
  tc_fragment[i] >= 0    -> atom i belongs to mobile fragment with that id

The array lives in atoms.arrays["tc_fragment"] so it survives copy / extxyz /
MD / supercell automatically (ASE propagates per-atom arrays).
"""

from __future__ import annotations

import numpy as np
from ase import Atoms

FRAGMENT_KEY = "tc_fragment"
FRAMEWORK = -1


def get_fragments(atoms: Atoms) -> np.ndarray | None:
    """Return the per-atom fragment array, or None if unset."""
    if FRAGMENT_KEY in atoms.arrays:
        return atoms.arrays[FRAGMENT_KEY].astype(int)
    return None


def set_fragments(atoms: Atoms, frag: np.ndarray | list[int]) -> None:
    """Attach/overwrite the per-atom fragment array (length must equal len(atoms))."""
    frag = np.asarray(frag, dtype=int)
    if frag.shape != (len(atoms),):
        raise ValueError(f"fragment array must have shape ({len(atoms)},), got {frag.shape}")
    atoms.set_array(FRAGMENT_KEY, frag)


def fragment_ids(atoms: Atoms) -> list[int]:
    """Sorted list of mobile fragment ids (excludes FRAMEWORK == -1)."""
    frag = get_fragments(atoms)
    if frag is None:
        return []
    return sorted(int(i) for i in np.unique(frag) if i != FRAMEWORK)


def fragment_mask(atoms: Atoms, fid: int) -> np.ndarray:
    """Boolean mask selecting atoms of fragment `fid`."""
    frag = get_fragments(atoms)
    if frag is None:
        raise ValueError("no fragment array set on these atoms")
    return frag == fid


def infer_fragments(
    atoms: Atoms,
    scale: float = 1.2,
    framework_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Assign fragment ids by connected components of a covalent-radius graph.

    Two atoms bond if distance < scale * (r_cov[i] + r_cov[j]).
    `framework_mask` (optional bool array, length == len(atoms)): atoms marked
    True are forced to FRAMEWORK (-1) and excluded from the connectivity graph.
    Returns the array; does NOT mutate `atoms`.
    """
    from ase.neighborlist import NeighborList, natural_cutoffs
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    n = len(atoms)
    result = np.full(n, FRAMEWORK, dtype=int)

    # Determine which atoms are mobile (not in the framework mask).
    mobile = np.ones(n, dtype=bool)
    if framework_mask is not None:
        framework_mask = np.asarray(framework_mask, dtype=bool)
        if framework_mask.shape != (n,):
            raise ValueError(
                f"framework_mask must have shape ({n},), got {framework_mask.shape}"
            )
        mobile[framework_mask] = False

    mobile_idx = np.where(mobile)[0]
    if len(mobile_idx) == 0:
        return result

    cutoffs = natural_cutoffs(atoms, mult=scale)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    # Build adjacency only among mobile atoms.
    idx_map = {int(i): j for j, i in enumerate(mobile_idx)}
    m = len(mobile_idx)
    rows, cols = [], []
    for local, global_i in enumerate(mobile_idx):
        neighbours, _ = nl.get_neighbors(global_i)
        for global_j in neighbours:
            if global_j in idx_map:
                rows.append(local)
                cols.append(idx_map[global_j])

    if rows:
        data = np.ones(len(rows), dtype=np.int8)
        adj = csr_matrix((data, (rows, cols)), shape=(m, m))
    else:
        adj = csr_matrix((m, m), dtype=np.int8)

    n_components, labels = connected_components(adj, directed=False)
    for local, global_i in enumerate(mobile_idx):
        result[global_i] = int(labels[local])

    return result
