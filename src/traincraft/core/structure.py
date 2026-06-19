"""``Structure``: the unit that flows through the whole pipeline.

A light wrapper around :class:`ase.Atoms` that also carries computed properties
and provenance, plus a stable content hash used for dedup and idempotent jobs.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from ase import Atoms

from .provenance import Provenance

# Per-structure / per-atom properties we care about end to end.
PROPERTY_KEYS = ("energy", "forces", "stress", "dipole", "polarizability")


@dataclass
class Structure:
    atoms: Atoms
    properties: dict[str, Any] = field(default_factory=dict)
    provenance: Provenance = field(default_factory=Provenance)

    @property
    def hash(self) -> str:
        """Content hash from composition + geometry (rounded for stability)."""
        a = self.atoms
        payload = {
            "numbers": a.get_atomic_numbers().tolist(),
            "positions": np.round(a.get_positions(), 4).tolist(),
            "cell": np.round(np.asarray(a.get_cell()), 4).tolist(),
            "pbc": a.get_pbc().tolist(),
        }
        blob = json.dumps(payload, sort_keys=True).encode()
        return hashlib.sha1(blob).hexdigest()[:16]

    @classmethod
    def from_ase(cls, atoms: Atoms, **kwargs: Any) -> Structure:
        return cls(atoms=atoms.copy(), **kwargs)

    def to_ase(self, with_properties: bool = True) -> Atoms:
        """Return a copy of the atoms with properties/provenance in ``info``."""
        atoms = self.atoms.copy()
        atoms.info["tc_provenance"] = self.provenance.to_dict()
        atoms.info["tc_hash"] = self.hash
        if with_properties:
            for key, value in self.properties.items():
                if key == "forces" and value is not None:
                    atoms.arrays["tc_forces"] = np.asarray(value)
                else:
                    atoms.info[f"tc_{key}"] = value
        return atoms

    def copy(self) -> Structure:
        return Structure(
            atoms=self.atoms.copy(),
            properties=dict(self.properties),
            provenance=Provenance.from_dict(self.provenance.to_dict()),
        )

    # --- fragment identity helpers ----------------------------------------
    @property
    def fragments(self):
        """Per-atom fragment array, or None if unset."""
        from .fragments import get_fragments
        return get_fragments(self.atoms)

    def set_fragments(self, frag) -> None:
        """Attach/overwrite the per-atom fragment array."""
        from .fragments import set_fragments
        set_fragments(self.atoms, frag)

    @property
    def n_fragments(self) -> int:
        """Number of distinct mobile fragments (excludes framework atoms)."""
        from .fragments import fragment_ids
        return len(fragment_ids(self.atoms))

    # --- interop (see core.converter) -------------------------------------
    def to_pymatgen(self):
        """Return a pymatgen ``Structure`` (periodic) or ``Molecule``."""
        from .converter import ase_to_pymatgen
        return ase_to_pymatgen(self.atoms)

    def to_rdkit(self, charge: int = 0):
        """Return an RDKit ``Mol`` with bonds perceived (non-periodic only)."""
        from .converter import ase_to_rdkit
        return ase_to_rdkit(self.atoms, charge=charge)

    @classmethod
    def from_pymatgen(cls, obj, **kwargs: Any) -> Structure:
        """Build a :class:`Structure` from a pymatgen ``Structure``/``Molecule``."""
        from .converter import pymatgen_to_ase
        return cls.from_ase(pymatgen_to_ase(obj), **kwargs)

    @classmethod
    def from_rdkit(cls, mol, conf_id: int = 0, **kwargs: Any) -> Structure:
        """Build a :class:`Structure` from one conformer of an RDKit ``Mol``."""
        from .converter import rdkit_to_ase
        return cls.from_ase(rdkit_to_ase(mol, conf_id=conf_id), **kwargs)
