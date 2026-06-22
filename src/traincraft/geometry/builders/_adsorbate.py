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


def _reject_smiles_whitespace(smiles: str) -> None:
    """Reject a SMILES containing whitespace.

    In SMILES, a space separates the structure from an optional molecule *name*,
    so RDKit silently keeps only the text before the first space (e.g. ``"OCCCCC O"``
    parses as pentanol, dropping the second oxygen). That would build the wrong
    molecule with no error, so we refuse it up front.
    """
    if any(ch.isspace() for ch in smiles):
        raise ValueError(
            f"SMILES {smiles!r} contains whitespace. In SMILES a space separates the "
            "structure from an optional molecule name, so RDKit would silently keep "
            "only the part before the space and drop the rest. Remove the space "
            "(e.g. 'OCCCCCO' for pentane-1,5-diol)."
        )


def _resolve_adsorbate(molecule_name, smiles, file) -> Atoms:
    """Return an ``ase.Atoms`` for a species from whichever field is set.

    The SMILES path is defensive: tiny/symmetric molecules (water is the classic
    case) can make RDKit's distance-geometry embedder fail or return a degenerate
    conformer, which would otherwise silently yield a heavy-atom-only structure
    (e.g. water without its hydrogens). We retry, then verify the result, and
    raise a clear error rather than emit a broken geometry into a dataset.
    """
    if molecule_name is not None:
        return molecule(molecule_name)
    if smiles is not None:
        _reject_smiles_whitespace(smiles)
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError as e:
            raise ImportError(
                "SMILES species needs RDKit. "
                "Install it with: pixi install -e science  (or: uv pip install rdkit)"
            ) from e
        parsed = Chem.MolFromSmiles(smiles)
        if parsed is None:
            raise ValueError(f"could not parse SMILES {smiles!r}")
        mol = Chem.AddHs(parsed)
        n_expected = mol.GetNumAtoms()  # includes the explicit hydrogens

        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        if AllChem.EmbedMolecule(mol, params) != 0:
            # ETKDG can fail for small/symmetric molecules; random coords help.
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, params) != 0:
                raise RuntimeError(
                    f"RDKit could not generate 3D coordinates for SMILES {smiles!r}. "
                    "For a small molecule, use molecule_name (ASE g2) instead — "
                    'e.g. molecule_name = "H2O".'
                )
        try:
            AllChem.MMFFOptimizeMolecule(mol)  # best-effort; geometry is valid without it
        except Exception:  # noqa: BLE001 - some molecules lack MMFF params
            pass

        conf = mol.GetConformer()
        positions = np.array(conf.GetPositions())
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        atoms = Atoms(symbols=symbols, positions=positions)

        # Guard against silent atom loss / collapsed embeddings.
        if len(atoms) != n_expected:
            raise RuntimeError(
                f"SMILES {smiles!r}: expected {n_expected} atoms but built "
                f"{len(atoms)} — RDKit dropped atoms. Use molecule_name instead."
            )
        if len(atoms) > 1 and float(np.ptp(positions, axis=0).max()) < 1e-3:
            raise RuntimeError(
                f"SMILES {smiles!r}: RDKit returned a degenerate (collapsed) "
                "geometry. Use molecule_name (ASE g2) for small molecules."
            )
        return atoms
    if file is not None:
        return read(file)
    raise ValueError("no species specified")


def _canonical_smiles(smiles: str) -> str:
    from rdkit import Chem

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
