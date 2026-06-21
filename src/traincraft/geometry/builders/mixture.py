"""Shared *mixture* machinery: turn a list of :class:`Species` into packed atoms.

Every builder that *places molecules* — ``liquid``, ``surface_packing`` and
``filled_nanotube`` — speaks the same language here:

  1. :func:`species_counts` turns ``count``/``ratio`` amounts into integer copies
     (ratios are apportioned over a total with the largest-remainder method so
     they always sum *exactly* to the requested number of molecules);
  2. :func:`resolve_mixture` builds an ``Atoms`` template per species;
  3. :func:`run_packmol` packs all species into one region (a Packmol ``inside``
     line chosen by the caller: a box, a surface slab, or a cylinder);
  4. :func:`tag_mixture` writes the per-molecule ``tc_fragment`` ids and records
     a per-fragment species/SMILES map for provenance.

Keeping this in one place is what makes mixtures *interusable*: a solvent blend,
a mixed adlayer and a co-filled nanotube differ only in the region line.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii, vdw_radii
from ase.io import read, write

from ._adsorbate import _canonical_smiles, _resolve_adsorbate


def vdw_radius(symbol: str) -> float:
    """Van der Waals radius (Å) for an element, with a covalent-based fallback.

    ``ase.data.vdw_radii`` has NaN for many elements; there we estimate from the
    covalent radius. Used to size packing regions by real atomic extent rather
    than a guessed margin.
    """
    z = atomic_numbers[symbol]
    r = vdw_radii[z]
    if np.isnan(r):
        return float(covalent_radii[z]) + 0.8
    return float(r)


def apportion(weights: list[float], total: int) -> list[int]:
    """Distribute ``total`` integer units across ``weights`` (largest remainder).

    The result sums *exactly* to ``total``; the leftover after flooring goes to
    the species with the largest fractional parts. A zero total gives all zeros.
    """
    s = float(sum(weights))
    if s <= 0:
        raise ValueError("mixture weights must be positive")
    raw = [w / s * total for w in weights]
    base = [int(np.floor(x)) for x in raw]
    rem = total - sum(base)
    order = sorted(range(len(raw)), key=lambda i: raw[i] - base[i], reverse=True)
    for i in order[:rem]:
        base[i] += 1
    return base


def species_counts(species: list, total: int | None) -> list[int]:
    """Resolve each species' integer copy count.

    Count mode (no ``ratio`` anywhere): each species' ``count`` (default 1).
    Ratio mode (any ``ratio`` set): apportion ``total`` by ratio; a bare species
    counts as an equal share (ratio 1).
    """
    has_ratio = any(s.ratio is not None for s in species)
    if not has_ratio:
        return [s.count if s.count is not None else 1 for s in species]
    if total is None:
        raise ValueError(
            "ratio-based mixture needs a total number of molecules (n_molecules)"
        )
    weights = [s.ratio if s.ratio is not None else 1.0 for s in species]
    return apportion(weights, total)


def resolve_mixture(species: list, total: int | None = None) -> list[tuple]:
    """Return ``[(atoms, count, label, canonical_smiles_or_None), ...]``."""
    counts = species_counts(species, total)
    out: list[tuple] = []
    for s, n in zip(species, counts, strict=True):
        atoms = _resolve_adsorbate(s.molecule_name, s.smiles, s.file)
        canonical = _canonical_smiles(s.smiles) if s.smiles is not None else None
        out.append((atoms, n, s.label, canonical))
    return out


def run_packmol(items: list[tuple[Atoms, int]], region_line: str, tol: float, seed: int,
                *, need: str = "this builder", fixed: Atoms | None = None) -> Atoms:
    """Pack ``items`` (atoms, count) into a single Packmol ``region_line``.

    ``region_line`` is a full Packmol constraint such as
    ``"inside box 0 0 0 12 12 12"`` or ``"inside cylinder ..."``. All species
    share the one region; Packmol returns them in listed order (every copy of
    species 0, then species 1, …), which :func:`tag_mixture` relies on.

    If ``fixed`` is given (e.g. a nanotube wall or a slab), it is written first
    as a Packmol *fixed structure* — kept at its exact coordinates and treated as
    an obstacle, so every packed atom is at least ``tol`` from it. The returned
    Atoms then begins with the fixed structure, followed by the packed species.
    """
    if shutil.which("packmol") is None:
        raise ImportError(f"{need} needs Packmol. Install it with: pixi install -e science")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        out_file = tmp / "packed.xyz"
        blocks = [f"tolerance {tol}", f"seed {seed}", f"output {out_file}", "filetype xyz", ""]
        if fixed is not None:
            fixed_file = tmp / "fixed.xyz"
            write(str(fixed_file), fixed, format="xyz")
            blocks += [
                f"structure {fixed_file}",
                "  number 1",
                "  fixed 0. 0. 0. 0. 0. 0.",
                "end structure",
                "",
            ]
        for i, (atoms, n) in enumerate(items):
            sp_file = tmp / f"species_{i}.xyz"
            write(str(sp_file), atoms, format="xyz")
            blocks += [
                f"structure {sp_file}",
                f"  number {n}",
                f"  {region_line}",
                "end structure",
                "",
            ]
        inp_file = tmp / "pack.inp"
        inp_file.write_text("\n".join(blocks))
        with open(inp_file) as inp_fh:
            result = subprocess.run(
                ["packmol"], stdin=inp_fh, capture_output=True, text=True, timeout=300
            )
        if result.returncode != 0 or not out_file.exists():
            raise RuntimeError(
                f"Packmol failed (exit {result.returncode}); the region may be too "
                f"small for the requested molecules. Output:\n{result.stdout[-2000:]}"
            )
        return read(str(out_file), format="xyz")


def tag_mixture(frag: np.ndarray, resolved: list[tuple], cursor: int) -> tuple:
    """Fill ``frag[cursor:]`` with per-molecule fragment ids for ``resolved``.

    Returns ``(n_fragments, fragment_smiles, fragment_species, end_cursor)``.
    ``frag`` ids restart at 0 for the first guest molecule; framework atoms
    (if any) must already be tagged before ``cursor``.
    """
    fragment_smiles: dict[str, str] = {}
    fragment_species: dict[str, str] = {}
    fid = 0
    for atoms, n, label, canonical in resolved:
        n_at = len(atoms)
        for _ in range(n):
            frag[cursor : cursor + n_at] = fid
            fragment_species[str(fid)] = label
            if canonical is not None:
                fragment_smiles[str(fid)] = canonical
            cursor += n_at
            fid += 1
    return fid, fragment_smiles, fragment_species, cursor
