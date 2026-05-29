"""Geometry subsystem: Source x Builder x Transform.

A geometry workflow is a declared ``source|builder`` followed by ordered
``transforms``. Plugins self-register on import below.
"""

from __future__ import annotations

from ..core import Structure, get
from . import builders, sources, transforms  # noqa: F401  (register plugins)


def build_source(cfg) -> Structure:
    return get("source", cfg.type)(cfg)


def build_builder(cfg) -> Structure:
    return get("builder", cfg.type)(cfg)


def apply_transform(structure: Structure, cfg) -> Structure:
    return get("transform", cfg.type)(structure, cfg)


def build_geometry(geom_cfg) -> Structure:
    """Resolve a :class:`GeometryConfig` into a single :class:`Structure`."""
    if geom_cfg.builder is not None:
        structure = build_builder(geom_cfg.builder)
    else:
        structure = build_source(geom_cfg.source)
    for transform_cfg in geom_cfg.transforms:
        structure = apply_transform(structure, transform_cfg)
    return structure


__all__ = ["apply_transform", "build_builder", "build_geometry", "build_source"]
