"""Low-level molecule resolvers shared by the molecule-placing builders.

Kept dependency-free of the higher-level mixture/surface modules so that
``mixture.py`` (and the alloy helper that imports it) can use these without an
import cycle.
"""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import read


def _resolve_adsorbate(molecule_name, smiles, file) -> Atoms:
    """Return an ``ase.Atoms`` for a species from whichever field is set."""
    if molecule_name is not None:
        return molecule(molecule_name)
    if smiles is not None:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError as e:
            raise ImportError(
                "SMILES species needs RDKit. "
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
    raise ValueError("no species specified")


def _canonical_smiles(smiles: str) -> str:
    from rdkit import Chem

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
