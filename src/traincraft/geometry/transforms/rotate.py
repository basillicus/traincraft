"""Transform: rotate the atoms (and optionally the cell) about an axis.

``axis`` is either a named Cartesian axis (``"x"``/``"y"``/``"z"``) or an
explicit 3-vector. ``angle`` is in degrees. Rotation is about the centre of
mass by default. Set ``rotate_cell=True`` to rotate the periodic cell with the
atoms (a rigid reorientation); leave it False to rotate atoms within a fixed
cell (useful for molecules in a vacuum box).
"""

from __future__ import annotations

from ...core import Structure, register

_NAMED_AXES = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}


@register("transform", "rotate")
def transform_rotate(structure: Structure, cfg) -> Structure:
    axis = _NAMED_AXES[cfg.axis] if isinstance(cfg.axis, str) else tuple(cfg.axis)
    out = structure.copy()
    out.atoms.rotate(cfg.angle, axis, center="COM", rotate_cell=cfg.rotate_cell)
    out.properties = {}  # geometry changed: stale properties dropped
    out.provenance.transforms.append(f"rotate:{cfg.angle}@{cfg.axis}")
    return out
