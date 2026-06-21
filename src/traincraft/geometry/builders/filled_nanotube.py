"""Builder: a carbon nanotube randomly filled with molecules ("fillMyTubes").

The original driver for this project: take a CNT and stuff N small molecules
inside it, then sample/label to study (e.g.) confined-phase Raman. Packmol packs
the guests into a cylinder coaxial with the tube and just narrower than its wall;
the tube carbons are tagged framework (``tc_fragment == -1``) and every guest
molecule gets its own fragment id, so the MC sampler and the ``constraints``
transform can address tube and guests separately.

Reuses the Packmol/adsorbate machinery from the ``surface_*``/``liquid`` builders;
only the Packmol *region* differs (``inside cylinder`` rather than ``inside box``).
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import nanotube
from ase.io import read, write

from ...core import FRAMEWORK, Provenance, Structure, register, set_fragments
from .surface import _canonical_smiles, _resolve_adsorbate


def _tube_axis_and_radius(cnt: Atoms) -> tuple[float, float, float]:
    """Return (cx, cy, R): the tube's xy axis centre and mean carbon radius."""
    pos = cnt.get_positions()
    cx, cy = pos[:, 0].mean(), pos[:, 1].mean()
    radii = np.hypot(pos[:, 0] - cx, pos[:, 1] - cy)
    return float(cx), float(cy), float(radii.mean())


def _run_packmol_cylinder(
    guest: Atoms,
    n: int,
    *,
    centre: tuple[float, float],
    z_lo: float,
    length: float,
    radius: float,
    tol: float,
    seed: int,
) -> Atoms:
    """Pack ``n`` copies of ``guest`` inside a z-aligned cylinder via Packmol."""
    if shutil.which("packmol") is None:
        raise ImportError(
            "filled_nanotube needs Packmol. Install it with: pixi install -e science"
        )
    cx, cy = centre
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        out_file = tmp / "packed.xyz"
        guest_file = tmp / "guest.xyz"
        write(str(guest_file), guest, format="xyz")
        # inside cylinder a1 a2 a3 d1 d2 d3 r l  — base centre, axis direction, radius, length
        blocks = [
            f"tolerance {tol}",
            f"seed {seed}",
            f"output {out_file}",
            "filetype xyz",
            "",
            f"structure {guest_file}",
            f"  number {n}",
            f"  inside cylinder {cx:.4f} {cy:.4f} {z_lo:.4f} 0. 0. 1. {radius:.4f} {length:.4f}",
            "end structure",
            "",
        ]
        inp_file = tmp / "pack.inp"
        inp_file.write_text("\n".join(blocks))
        with open(inp_file) as inp_fh:
            result = subprocess.run(
                ["packmol"], stdin=inp_fh, capture_output=True, text=True, timeout=300
            )
        if result.returncode != 0 or not out_file.exists():
            raise RuntimeError(
                f"Packmol failed (exit {result.returncode}); the tube may be too "
                f"narrow/short for {n} guest(s). Output:\n{result.stdout[-2000:]}"
            )
        return read(str(out_file), format="xyz")


@register("builder", "filled_nanotube")
def build_filled_nanotube(cfg) -> Structure:
    cnt = nanotube(cfg.n, cfg.m, length=cfg.length, bond=cfg.bond)
    cnt.center(vacuum=cfg.vacuum / 2, axis=(0, 1))

    cx, cy, radius = _tube_axis_and_radius(cnt)
    length_z = float(cnt.get_cell()[2][2])
    pack_radius = radius - cfg.radial_margin
    pack_length = length_z - 2 * cfg.axial_margin
    if pack_radius <= 0 or pack_length <= 0:
        raise ValueError(
            f"tube too small to fill: usable radius {pack_radius:.2f} Å, length "
            f"{pack_length:.2f} Å. Widen the tube (n/m) or lengthen it (length), "
            "or reduce radial_margin/axial_margin."
        )

    guest = _resolve_adsorbate(cfg.molecule_name, cfg.smiles, cfg.file)
    seed = cfg.seed if cfg.seed is not None else 12345
    guests = _run_packmol_cylinder(
        guest,
        cfg.n_molecules,
        centre=(cx, cy),
        z_lo=cfg.axial_margin,
        length=pack_length,
        radius=pack_radius,
        tol=cfg.tolerance,
        seed=seed,
    )

    # Combine: tube first (framework), then the packed guests (one fragment each).
    system = cnt + guests
    system.set_cell(cnt.get_cell())
    system.set_pbc((False, False, cfg.pbc))

    frag = np.empty(len(system), dtype=int)
    frag[: len(cnt)] = FRAMEWORK
    n_at = len(guest)
    cursor = len(cnt)
    for fid in range(cfg.n_molecules):
        frag[cursor : cursor + n_at] = fid
        cursor += n_at
    set_fragments(system, frag)

    extra: dict = {
        "n_molecules": cfg.n_molecules,
        "tube_radius": round(radius, 4),
        "pack_radius": round(pack_radius, 4),
    }
    if cfg.smiles is not None:
        canonical = _canonical_smiles(cfg.smiles)
        extra["fragment_smiles"] = {str(i): canonical for i in range(cfg.n_molecules)}

    return Structure.from_ase(
        system,
        provenance=Provenance(
            origin="generated",
            source=f"builder:filled_nanotube:{cfg.n}-{cfg.m}-l{cfg.length}",
            extra=extra,
        ),
    )
