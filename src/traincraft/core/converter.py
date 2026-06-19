"""Interop bridge: ASE ↔ pymatgen ↔ RDKit.

ASE is the internal lingua franca (everything flows as :class:`ase.Atoms`).
This module converts to/from the two other ecosystems we lean on:

  * **pymatgen** — crystal-structure analysis, symmetry, providers (Materials
    Project / OPTIMADE).  Periodic ``Atoms`` map to a pymatgen ``Structure``;
    non-periodic ``Atoms`` map to a ``Molecule``.
  * **RDKit** — cheminformatics, SMILES, conformers.  RDKit molecules are
    inherently non-periodic; converting a periodic cell raises.

Both libraries are optional dependencies (``pixi install -e science``); the
import error message points there.
"""

from __future__ import annotations

import io
from typing import TYPE_CHECKING

import numpy as np
from ase import Atoms
from ase.io import write

if TYPE_CHECKING:  # pragma: no cover - typing only
    from pymatgen.core import Molecule, Structure
    from rdkit.Chem import Mol


# --------------------------------------------------------------------------- pymatgen
def _require_pymatgen():
    try:
        from pymatgen.io.ase import AseAtomsAdaptor
    except ImportError as e:  # pragma: no cover - exercised only without pymatgen
        raise ImportError(
            "pymatgen interop needs pymatgen. "
            "Install it with: pixi install -e science  (or: uv pip install pymatgen)"
        ) from e
    return AseAtomsAdaptor


def ase_to_pymatgen(atoms: Atoms) -> Structure | Molecule:
    """Convert ASE ``Atoms`` to a pymatgen ``Structure`` (periodic) or ``Molecule``.

    The choice is driven by periodicity: an ``Atoms`` periodic in all three
    directions becomes a ``Structure``; otherwise a ``Molecule`` (the cell is
    dropped, since a partially periodic slab/wire has no pymatgen analogue).
    """
    adaptor = _require_pymatgen()
    if bool(np.all(atoms.get_pbc())):
        return adaptor.get_structure(atoms)
    return adaptor.get_molecule(atoms)


def pymatgen_to_ase(obj: Structure | Molecule) -> Atoms:
    """Convert a pymatgen ``Structure`` or ``Molecule`` to ASE ``Atoms``."""
    adaptor = _require_pymatgen()
    return adaptor.get_atoms(obj)


# ----------------------------------------------------------------------------- rdkit
def _require_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError as e:  # pragma: no cover - exercised only without rdkit
        raise ImportError(
            "RDKit interop needs RDKit. "
            "Install it with: pixi install -e science  (or: uv pip install rdkit)"
        ) from e
    return Chem, rdDetermineBonds


def ase_to_rdkit(atoms: Atoms, charge: int = 0) -> Mol:
    """Convert non-periodic ASE ``Atoms`` to an RDKit ``Mol`` with bonds perceived.

    Bonds are inferred from the 3D geometry by RDKit's ``DetermineBonds`` (the
    xyz2mol algorithm).  Raises if the structure is periodic in any direction.
    """
    if bool(np.any(atoms.get_pbc())):
        raise ValueError(
            "ase_to_rdkit needs a non-periodic structure; got pbc="
            f"{atoms.get_pbc().tolist()}. RDKit molecules are not periodic."
        )
    Chem, rdDetermineBonds = _require_rdkit()

    buf = io.StringIO()
    write(buf, atoms, format="xyz")
    mol = Chem.MolFromXYZBlock(buf.getvalue())
    if mol is None:
        raise ValueError("RDKit could not parse the structure as an XYZ molecule")
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
    except ValueError as e:
        raise ValueError(
            f"RDKit could not perceive bonds (charge={charge}): {e}. "
            "Try passing the correct total charge."
        ) from e
    return mol


def rdkit_to_ase(mol: Mol, conf_id: int = 0) -> Atoms:
    """Convert one conformer of an RDKit ``Mol`` to ASE ``Atoms``.

    ``conf_id`` selects which embedded conformer to read (default: the first).
    """
    if mol.GetNumConformers() == 0:
        raise ValueError(
            "RDKit molecule has no conformers; embed one first "
            "(e.g. AllChem.EmbedMolecule)."
        )
    conf = mol.GetConformer(conf_id)
    positions = np.asarray(conf.GetPositions())
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return Atoms(symbols=symbols, positions=positions)


__all__ = [
    "ase_to_pymatgen",
    "ase_to_rdkit",
    "pymatgen_to_ase",
    "rdkit_to_ase",
]
