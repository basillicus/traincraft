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

    # --- interop stubs (implemented in the geometry chunk) -----------------
    def to_pymatgen(self):  # pragma: no cover - Phase 2
        raise NotImplementedError("pymatgen interop lands with the geometry chunk")

    def to_rdkit(self):  # pragma: no cover - Phase 2
        raise NotImplementedError("rdkit interop lands with the geometry chunk")
