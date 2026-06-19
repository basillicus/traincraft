"""Builder: liquids, mixtures, and confined bulk via Packmol.

Packs one or more molecular species into a periodic box. The box is either given
explicitly (``box = [a, b, c]`` in Å) or derived from a target mass ``density``
(g/cm³) for a cubic cell. Each packed molecule becomes its own ``tc_fragment``
(numbered globally across species), so the MC sampler and constraints can address
them individually.

Reuses the Packmol/adsorbate machinery from the ``surface_*`` builders.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write

from ...core import Provenance, Structure, register, set_fragments
from .surface import _canonical_smiles, _resolve_adsorbate

_AMU_G = 1.66053906660e-24  # grams per atomic mass unit


def _cubic_edge_from_density(total_mass_amu: float, density: float) -> float:
    """Cubic box edge (Å) holding ``total_mass_amu`` at ``density`` g/cm³."""
    volume_cm3 = (total_mass_amu * _AMU_G) / density
    edge_cm = volume_cm3 ** (1.0 / 3.0)
    return edge_cm * 1.0e8  # cm -> Å


def _run_packmol_multi(
    species: list[tuple[Atoms, int]],
    box: tuple[float, float, float],
    margin: float,
    tol: float,
    seed: int,
) -> Atoms:
    """Pack several species (atoms, count) into ``box`` (Å), inset by ``margin``."""
    if shutil.which("packmol") is None:
        raise ImportError("liquid builder needs Packmol. Install it with: pixi install -e science")

    a, b, c = box
    lo, hi = margin, (a - margin, b - margin, c - margin)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        out_file = tmp / "packed.xyz"
        blocks = [f"tolerance {tol}", f"seed {seed}", f"output {out_file}", "filetype xyz", ""]
        for i, (atoms, n) in enumerate(species):
            sp_file = tmp / f"species_{i}.xyz"
            write(str(sp_file), atoms, format="xyz")
            blocks += [
                f"structure {sp_file}",
                f"  number {n}",
                f"  inside box {lo:.4f} {lo:.4f} {lo:.4f} {hi[0]:.4f} {hi[1]:.4f} {hi[2]:.4f}",
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
                f"Packmol failed (exit {result.returncode}):\n{result.stdout[-2000:]}"
            )
        return read(str(out_file), format="xyz")


@register("builder", "liquid")
def build_liquid(cfg) -> Structure:
    # Resolve each species to an Atoms template + its count.
    resolved: list[tuple[Atoms, int]] = []
    for spec in cfg.species:
        atoms = _resolve_adsorbate(spec.molecule_name, spec.smiles, spec.file)
        resolved.append((atoms, spec.count))

    # Box: explicit, or cubic from target density.
    if cfg.box is not None:
        box = tuple(float(x) for x in cfg.box)
    elif cfg.density is not None:
        total_mass = sum(float(a.get_masses().sum()) * n for a, n in resolved)
        edge = _cubic_edge_from_density(total_mass, cfg.density)
        box = (edge, edge, edge)
    else:
        raise ValueError("liquid builder needs either 'box' or 'density'")

    seed = cfg.seed if cfg.seed is not None else 12345
    # Inset packing by `tolerance` from each wall so molecules don't clash across
    # the periodic boundary (Packmol itself is not PBC-aware).
    packed = _run_packmol_multi(resolved, box, margin=cfg.tolerance, tol=cfg.tolerance, seed=seed)

    packed.set_cell(box)
    packed.set_pbc(cfg.pbc)

    # Fragment tagging: every molecule (across all species) is its own fragment.
    frag = np.empty(len(packed), dtype=int)
    fragment_smiles: dict[str, str] = {}
    cursor = fid = 0
    for (atoms, n), spec in zip(resolved, cfg.species, strict=True):
        n_at = len(atoms)
        canonical = _canonical_smiles(spec.smiles) if spec.smiles is not None else None
        for _ in range(n):
            frag[cursor : cursor + n_at] = fid
            if canonical is not None:
                fragment_smiles[str(fid)] = canonical
            cursor += n_at
            fid += 1
    set_fragments(packed, frag)

    extra: dict = {"n_molecules": int(fid), "box": [round(x, 4) for x in box]}
    if fragment_smiles:
        extra["fragment_smiles"] = fragment_smiles

    return Structure.from_ase(
        packed,
        provenance=Provenance(origin="generated", source="builder:liquid", extra=extra),
    )
