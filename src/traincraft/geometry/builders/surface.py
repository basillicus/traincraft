"""Builders: molecule(s) on a crystalline surface.

Two builders:
  surface_adsorbate  — single adsorbate via ase.build.add_adsorbate (zero extra deps)
  surface_packing    — N-molecule coverage via the Packmol binary (conda-forge)

Fragment convention:
  tc_fragment == -1   substrate atoms (immobile)
  tc_fragment >= 0    each adsorbate molecule is a separate mobile fragment
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import add_adsorbate, molecule
from ase.io import read, write

from ...core import Provenance, Structure, register, set_fragments
from ._common import build_named_slab


def _build_slab(element: str, facet: str, size: tuple, vacuum: float) -> Atoms:
    return build_named_slab(element, facet, size, vacuum)


def _resolve_adsorbate(molecule_name, smiles, file) -> Atoms:
    """Return an ase.Atoms for the adsorbate from whichever field is set."""
    if molecule_name is not None:
        return molecule(molecule_name)
    if smiles is not None:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError as e:
            raise ImportError(
                "SMILES adsorbate needs RDKit. "
                "Install it with: pixi install -e science  (or: uv pip install rdkit)"
            ) from e
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        AllChem.EmbedMolecule(mol, params)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer(0)
        positions = np.array(conf.GetPositions())
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        return Atoms(symbols=symbols, positions=positions)
    if file is not None:
        return read(file)
    raise ValueError("no adsorbate specified")


def _canonical_smiles(smiles: str) -> str:
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


@register("builder", "surface_adsorbate")
def build_surface_adsorbate(cfg) -> Structure:
    slab = _build_slab(cfg.element, cfg.facet, cfg.size, cfg.vacuum)
    n_slab = len(slab)

    ads = _resolve_adsorbate(cfg.molecule_name, cfg.smiles, cfg.file)
    n_ads = len(ads)

    kwargs = {"height": cfg.height, "position": cfg.site}
    if cfg.offset is not None:
        kwargs["offset"] = cfg.offset
    add_adsorbate(slab, ads, **kwargs)
    slab.info.pop("adsorbate_info", None)  # not extxyz-serialisable

    # Fragment tagging: substrate=-1, adsorbate=0
    frag = np.full(n_slab + n_ads, -1, dtype=int)
    frag[n_slab:] = 0
    set_fragments(slab, frag)

    extra: dict = {}
    if cfg.smiles is not None:
        canonical = _canonical_smiles(cfg.smiles)
        extra["smiles"] = canonical
        extra["fragment_smiles"] = {"0": canonical}

    return Structure.from_ase(
        slab,
        provenance=Provenance(
            origin="generated",
            source=f"builder:surface_adsorbate:{cfg.element}-{cfg.facet}",
            extra=extra,
        ),
    )


# ---------------------------------------------------------------------------
# Packmol helper

def _run_packmol(
    adsorbate: Atoms,
    n: int,
    box: tuple[float, float, float, float, float, float],
    tol: float,
    seed: int,
) -> Atoms:
    """Pack `n` copies of `adsorbate` into `box` (xmin ymin zmin xmax ymax zmax).

    Writes a temporary Packmol input deck, runs the binary, reads the result.
    Raises ImportError if the `packmol` binary is not on PATH.
    """
    if shutil.which("packmol") is None:
        raise ImportError(
            "surface_packing needs Packmol. "
            "Install it with: pixi install -e science"
        )

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        ads_file = tmp / "adsorbate.xyz"
        out_file = tmp / "packed.xyz"
        inp_file = tmp / "pack.inp"

        write(str(ads_file), adsorbate, format="xyz")

        deck = (
            f"tolerance {tol}\n"
            f"seed {seed}\n"
            f"output {out_file}\n"
            "filetype xyz\n\n"
            f"structure {ads_file}\n"
            f"  number {n}\n"
            "  inside box "
            f"{box[0]:.4f} {box[1]:.4f} {box[2]:.4f} "
            f"{box[3]:.4f} {box[4]:.4f} {box[5]:.4f}\n"
            "end structure\n"
        )
        inp_file.write_text(deck)

        with open(inp_file) as inp_fh:
            result = subprocess.run(
                ["packmol"],
                stdin=inp_fh,
                capture_output=True,
                text=True,
                timeout=120,
            )
        if result.returncode != 0 or not out_file.exists():
            raise RuntimeError(
                f"Packmol failed (exit {result.returncode}):\n{result.stdout[-2000:]}"
            )

        return read(str(out_file), format="xyz")


@register("builder", "surface_packing")
def build_surface_packing(cfg) -> Structure:
    slab = _build_slab(cfg.element, cfg.facet, cfg.size, cfg.vacuum)
    n_slab = len(slab)

    ads = _resolve_adsorbate(cfg.molecule_name, cfg.smiles, cfg.file)
    n_ads = len(ads)

    # Determine packing region: above slab top, within the xy cell footprint.
    cell = slab.get_cell()
    slab_pos = slab.get_positions()
    slab_top = float(slab_pos[:, 2].max())

    z_lo = slab_top + cfg.gap
    z_hi = z_lo + cfg.region_height
    # xy extent from the cell (assume orthorhombic or use bounding box)
    x_hi = float(np.linalg.norm(cell[0]))
    y_hi = float(np.linalg.norm(cell[1]))
    box = (0.0, 0.0, z_lo, x_hi, y_hi, z_hi)

    seed = cfg.seed if cfg.seed is not None else 12345
    packed = _run_packmol(ads, cfg.n_molecules, box, cfg.tolerance, seed)

    combined = slab + packed

    # Fragment tagging: substrate=-1, each molecule=0,1,2,...
    frag = np.full(n_slab + cfg.n_molecules * n_ads, -1, dtype=int)
    for i in range(cfg.n_molecules):
        start = n_slab + i * n_ads
        frag[start : start + n_ads] = i
    set_fragments(combined, frag)

    extra: dict = {"n_molecules": cfg.n_molecules}
    if cfg.smiles is not None:
        canonical = _canonical_smiles(cfg.smiles)
        extra["smiles"] = canonical
        extra["fragment_smiles"] = {str(i): canonical for i in range(cfg.n_molecules)}

    return Structure.from_ase(
        combined,
        provenance=Provenance(
            origin="generated",
            source=f"builder:surface_packing:{cfg.element}-{cfg.facet}",
            extra=extra,
        ),
    )
