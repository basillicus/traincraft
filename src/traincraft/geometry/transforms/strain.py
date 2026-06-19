"""Transform: apply strain to the cell (atoms scale with it).

Two ways to specify the strain:
  * ``hydrostatic`` — a single isotropic engineering strain applied to all three
    axes (positive = expansion);
  * ``voigt`` — a 6-component engineering-strain vector
    ``(e_xx, e_yy, e_zz, e_yz, e_xz, e_xy)``.

Both build a deformation gradient ``F = I + e`` and set ``new_cell = cell @ Fᵀ``
with ``scale_atoms=True`` so fractional coordinates are preserved. Strain only
makes sense for a periodic cell, so a non-zero strain on a zero cell raises.
"""

from __future__ import annotations

import numpy as np

from ...core import Structure, register


def _strain_tensor(cfg) -> np.ndarray:
    if cfg.hydrostatic is not None:
        return np.eye(3) * cfg.hydrostatic
    exx, eyy, ezz, eyz, exz, exy = cfg.voigt
    return np.array(
        [
            [exx, exy / 2, exz / 2],
            [exy / 2, eyy, eyz / 2],
            [exz / 2, eyz / 2, ezz],
        ]
    )


@register("transform", "strain")
def transform_strain(structure: Structure, cfg) -> Structure:
    out = structure.copy()
    cell = np.asarray(out.atoms.get_cell())
    if not cell.any():
        raise ValueError("strain transform needs a periodic cell (cell is all zero)")

    f = np.eye(3) + _strain_tensor(cfg)
    out.atoms.set_cell(cell @ f.T, scale_atoms=True)
    out.properties = {}  # geometry changed: stale properties dropped

    label = (
        f"strain:hydro={cfg.hydrostatic}"
        if cfg.hydrostatic is not None
        else f"strain:voigt={tuple(cfg.voigt)}"
    )
    out.provenance.transforms.append(label)
    return out
