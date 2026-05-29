"""Typed result of a calculation step."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class Result:
    properties: dict[str, Any] = field(default_factory=dict)
    success: bool = True
    calculator: str | None = None
    walltime: float | None = None
    error: str | None = None
    artifacts: dict[str, str] = field(default_factory=dict)
