"""Regression tests for two SMILES/conformer bugs.

1. A SMILES containing whitespace silently built the wrong molecule (RDKit treats
   text after a space as a molecule *name*), e.g. "OCCCCC O" -> pentanol.
2. The Monte Carlo conformer-swap move copied a freshly-built conformer's
   coordinates into the fragment index-for-index. When the rebuilt atom order
   differed from the fragment (it recorded the *canonical* SMILES), each element's
   coordinates landed on the wrong atom and silently corrupted the molecule
   (hydrogens "migrating" from the CH3 onto the OH).
"""

from __future__ import annotations

import numpy as np
import pytest

from traincraft.config.models import SmilesSource


def _o_neighbor_elements(atoms) -> list[str]:
    """Element symbols bonded to the (single) oxygen, by natural cutoffs."""
    from ase.neighborlist import NeighborList, natural_cutoffs

    nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True)
    nl.update(atoms)
    syms = atoms.get_chemical_symbols()
    oi = syms.index("O")
    nbrs, _ = nl.get_neighbors(oi)
    return sorted(syms[j] for j in nbrs)


# --- Bug 1: whitespace in SMILES -------------------------------------------

def test_smiles_source_config_rejects_whitespace():
    with pytest.raises(ValueError, match="whitespace"):
        SmilesSource(smiles="OCCCCC O")  # would silently parse as pentanol


def test_resolve_adsorbate_rejects_whitespace():
    pytest.importorskip("rdkit")
    from traincraft.geometry.builders._adsorbate import _resolve_adsorbate

    with pytest.raises(ValueError, match="whitespace"):
        _resolve_adsorbate(None, "OCCCCC O", None)


# --- Bug 2: conformer move must not scramble atoms --------------------------

def test_conformer_move_preserves_chemistry():
    pytest.importorskip("rdkit")
    from traincraft.geometry.sources.smiles import source_smiles
    from traincraft.sampling.monte_carlo import _move_conformer

    s = source_smiles(SmilesSource(smiles="OCCCCC", n_conformers=1, optimize=True))
    assert _o_neighbor_elements(s.atoms) == ["C", "H"]  # a hydroxyl, as built

    atoms = s.atoms.copy()
    frag_smiles = s.provenance.extra["fragment_smiles"]
    applied = _move_conformer(atoms, 0, np.random.default_rng(0), frag_smiles)
    assert applied is True
    # The OH is still an OH — no extra hydrogens migrated onto the oxygen.
    assert _o_neighbor_elements(atoms) == ["C", "H"]


def test_conformer_move_skips_on_reordered_smiles():
    """The old bug: a canonical (reordered) SMILES. The move must now REFUSE."""
    pytest.importorskip("rdkit")
    from rdkit import Chem

    from traincraft.geometry.sources.smiles import source_smiles
    from traincraft.sampling.monte_carlo import _move_conformer

    s = source_smiles(SmilesSource(smiles="OCCCCC"))
    atoms = s.atoms.copy()
    before = atoms.get_positions().copy()

    canonical = Chem.MolToSmiles(Chem.MolFromSmiles("OCCCCC"))  # reorders atoms
    applied = _move_conformer(atoms, 0, np.random.default_rng(0), {"0": canonical})

    assert applied is False  # refused rather than scrambling the molecule
    np.testing.assert_array_equal(atoms.get_positions(), before)  # untouched
