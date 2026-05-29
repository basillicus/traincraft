from __future__ import annotations

from ase.build import molecule

from traincraft import Structure
from traincraft.core import Provenance


def test_hash_stable_across_copies():
    a = molecule("H2O")
    assert Structure.from_ase(a).hash == Structure.from_ase(a.copy()).hash


def test_hash_changes_when_moved():
    a = molecule("H2O")
    b = a.copy()
    b.positions[0] += 1.0
    assert Structure.from_ase(a).hash != Structure.from_ase(b).hash


def test_to_ase_carries_metadata():
    a = molecule("CO")
    s = Structure.from_ase(
        a, properties={"energy": -1.23}, provenance=Provenance(origin="ml_sampled")
    )
    atoms = s.to_ase()
    assert atoms.info["tc_energy"] == -1.23
    assert atoms.info["tc_hash"] == s.hash
