"""Core abstractions: Structure, provenance, registry, workspace, results, rng."""

from __future__ import annotations

from .fragments import (
    FRAGMENT_KEY,
    FRAMEWORK,
    fragment_ids,
    fragment_mask,
    get_fragments,
    infer_fragments,
    set_fragments,
)
from .provenance import ORIGINS, Provenance
from .registry import RegistryError, available, capabilities, get, register
from .results import Result
from .rng import make_rng
from .structure import PROPERTY_KEYS, Structure
from .workspace import Job, Workspace

__all__ = [
    "FRAGMENT_KEY",
    "FRAMEWORK",
    "ORIGINS",
    "PROPERTY_KEYS",
    "Job",
    "Provenance",
    "RegistryError",
    "Result",
    "Structure",
    "Workspace",
    "available",
    "capabilities",
    "fragment_ids",
    "fragment_mask",
    "get_fragments",
    "get",
    "infer_fragments",
    "make_rng",
    "register",
    "set_fragments",
]
