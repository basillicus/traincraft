"""Source: build a molecule from a SMILES string via RDKit ETKDG."""

from __future__ import annotations

import numpy as np

from ...core import Provenance, Structure, register, set_fragments


@register("source", "smiles")
def source_smiles(cfg) -> Structure:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as e:
        raise ImportError(
            "SMILES source needs RDKit. "
            "Install it with: pixi install -e science  (or: uv pip install rdkit)"
        ) from e

    mol = Chem.MolFromSmiles(cfg.smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {cfg.smiles!r}")
    mol = Chem.AddHs(mol)

    seed = cfg.seed if cfg.seed is not None else 0xF00D
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1

    ids = AllChem.EmbedMultipleConfs(mol, numConfs=cfg.n_conformers, params=params)
    if len(ids) == 0:
        raise RuntimeError(
            f"ETKDG conformer embedding failed for SMILES {cfg.smiles!r}. "
            "Try a simpler molecule or check the SMILES string."
        )

    if cfg.optimize:
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1)

    conf = mol.GetConformer(0)
    positions = np.array(conf.GetPositions())
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    from ase import Atoms

    atoms = Atoms(symbols=symbols, positions=positions)
    atoms.center(vacuum=cfg.vacuum / 2)

    # All atoms form a single mobile fragment (id 0).
    set_fragments(atoms, np.zeros(len(atoms), dtype=int))

    canonical = Chem.MolToSmiles(Chem.MolFromSmiles(cfg.smiles))

    return Structure.from_ase(
        atoms,
        provenance=Provenance(
            origin="generated",
            source=f"source:smiles:{cfg.smiles}",
            extra={
                "smiles": canonical,
                "n_conformers": cfg.n_conformers,
                "fragment_smiles": {"0": canonical},
            },
        ),
    )
