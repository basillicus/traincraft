"""Provenance records: how a structure/frame came to be.

The ``origin`` tag keeps DFT-labeled data cleanly distinct from ML-generated
points so the expensive dataset stays separable and shareable.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

# Allowed origin tags, cheap (generated) -> expensive (dft_labeled).
ORIGINS = ("generated", "ml_sampled", "ml_labeled", "dft_labeled")


@dataclass
class Provenance:
    origin: str = "generated"
    source: str | None = None  # e.g. "builder:nanotube", "source:file"
    transforms: list[str] = field(default_factory=list)
    calculator: str | None = None  # method that produced ``properties``
    level_of_theory: dict[str, Any] = field(default_factory=dict)
    seed: int | None = None
    parents: list[str] = field(default_factory=list)  # parent structure hashes
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.origin not in ORIGINS:
            raise ValueError(f"origin must be one of {ORIGINS}, got {self.origin!r}")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Provenance:
        known = {f for f in cls.__dataclass_fields__}  # noqa: C416
        return cls(**{k: v for k, v in data.items() if k in known})
