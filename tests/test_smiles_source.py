from __future__ import annotations

import pytest

pytest.importorskip("rdkit", reason="rdkit not installed")

from traincraft.config.models import SmilesSource  # noqa: E402
from traincraft.core.fragments import fragment_ids, get_fragments  # noqa: E402
from traincraft.geometry.sources.smiles import source_smiles  # noqa: E402


def test_water_from_smiles():
    cfg = SmilesSource(smiles="O")
    s = source_smiles(cfg)

    assert len(s.atoms) == 3  # H2O: O + 2H
    assert set(s.atoms.get_chemical_symbols()) == {"H", "O"}


def test_fragments_all_zero():
    cfg = SmilesSource(smiles="O")
    s = source_smiles(cfg)

    frag = get_fragments(s.atoms)
    assert frag is not None
    assert list(frag) == [0] * len(s.atoms)
    assert fragment_ids(s.atoms) == [0]


def test_smiles_in_provenance_extra():
    cfg = SmilesSource(smiles="O")
    s = source_smiles(cfg)

    assert "smiles" in s.provenance.extra
    assert s.provenance.extra["smiles"] == "O"
    assert "fragment_smiles" in s.provenance.extra
    assert s.provenance.extra["fragment_smiles"]["0"] == "O"


def test_invalid_smiles_raises():
    cfg = SmilesSource(smiles="not_a_smiles_!!!")
    with pytest.raises((ValueError, RuntimeError)):
        source_smiles(cfg)


def test_ethanol_atom_count():
    cfg = SmilesSource(smiles="CCO")
    s = source_smiles(cfg)
    # ethanol C2H5OH -> 9 heavy+H atoms
    assert len(s.atoms) == 9


def test_deterministic_with_seed():
    cfg1 = SmilesSource(smiles="CCO", seed=42)
    cfg2 = SmilesSource(smiles="CCO", seed=42)
    s1 = source_smiles(cfg1)
    s2 = source_smiles(cfg2)
    import numpy as np
    np.testing.assert_allclose(
        s1.atoms.get_positions(), s2.atoms.get_positions(), atol=1e-4
    )
