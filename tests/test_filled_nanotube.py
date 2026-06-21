from __future__ import annotations

import shutil

import numpy as np
import pytest

from traincraft.config.models import FilledNanotubeBuilder, GeometryConfig
from traincraft.core import FRAMEWORK
from traincraft.geometry import build_geometry

_HAS_PACKMOL = shutil.which("packmol") is not None


# --- validation (no Packmol needed) ----------------------------------------

def test_needs_exactly_one_guest():
    with pytest.raises(ValueError, match="exactly one"):
        FilledNanotubeBuilder()  # no guest specified
    with pytest.raises(ValueError, match="exactly one"):
        FilledNanotubeBuilder(molecule_name="H2O", smiles="O")


def test_n_molecules_must_be_positive():
    with pytest.raises(ValueError, match="n_molecules"):
        FilledNanotubeBuilder(molecule_name="H2O", n_molecules=0)


def test_registered():
    from traincraft.core import get

    assert get("builder", "filled_nanotube") is not None


# --- end-to-end packing (needs the Packmol binary) -------------------------

@pytest.mark.skipif(not _HAS_PACKMOL, reason="packmol binary not installed")
def test_fills_tube_with_water():
    s = build_geometry(
        GeometryConfig(
            builder=FilledNanotubeBuilder(
                n=10, m=10, length=6, molecule_name="H2O", n_molecules=4, seed=1
            )
        )
    )
    atoms = s.atoms
    frag = atoms.arrays["tc_fragment"]
    # tube carbons are framework; 4 waters are fragments 0..3
    n_carbon = int((atoms.get_atomic_numbers() == 6).sum())
    assert (frag[:n_carbon] == FRAMEWORK).all()
    assert sorted(set(int(f) for f in frag if f >= 0)) == [0, 1, 2, 3]
    assert len(atoms) == n_carbon + 4 * 3  # tube + 4 * H2O
    # periodic along the tube (z) axis only
    assert list(atoms.get_pbc()) == [False, False, True]
    # guests sit inside the tube radius
    cx, cy = atoms.get_positions()[:, 0].mean(), atoms.get_positions()[:, 1].mean()
    guest_xy = atoms.get_positions()[n_carbon:, :2]
    radial = np.hypot(guest_xy[:, 0] - cx, guest_xy[:, 1] - cy)
    assert radial.max() < s.provenance.extra["tube_radius"]
    assert s.provenance.extra["n_molecules"] == 4


@pytest.mark.skipif(not _HAS_PACKMOL, reason="packmol binary not installed")
def test_tube_too_small_raises():
    # a narrow (5,0) tube cannot hold benzene
    with pytest.raises((ValueError, RuntimeError)):
        build_geometry(
            GeometryConfig(
                builder=FilledNanotubeBuilder(
                    n=5, m=0, length=4, molecule_name="C6H6", n_molecules=2, radial_margin=1.8
                )
            )
        )


def test_narrow_tube_refused_by_vdw_before_packmol():
    # An (8,0) tube physically cannot hold water without overlapping the wall;
    # the vdW sizing rejects it up front (no Packmol needed) with a clear message.
    with pytest.raises(ValueError, match="too small to fill without overlap"):
        build_geometry(
            GeometryConfig(
                builder=FilledNanotubeBuilder(
                    n=8, m=0, length=5, molecule_name="H2O", n_molecules=5
                )
            )
        )


@pytest.mark.skipif(not _HAS_PACKMOL, reason="packmol binary not installed")
def test_guests_do_not_overlap_wall():
    s = build_geometry(
        GeometryConfig(
            builder=FilledNanotubeBuilder(
                n=12, m=12, length=6, molecule_name="H2O", n_molecules=6,
                tolerance=2.5, seed=1,
            )
        )
    )
    atoms = s.atoms
    frag = atoms.arrays["tc_fragment"]
    wall = atoms.get_positions()[frag == FRAMEWORK]
    guest = atoms.get_positions()[frag >= 0]
    dmin = np.sqrt(((guest[:, None, :] - wall[None, :, :]) ** 2).sum(-1)).min()
    # Packmol's fixed-obstacle constraint keeps every guest atom >= tolerance away.
    assert dmin >= 2.5 - 1e-3
