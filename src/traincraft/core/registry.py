"""A tiny plugin registry.

New geometry builders, calculators, samplers, and selectors register themselves
by decorator; the engine resolves them by name from config. Adding a capability
means adding one file, never editing a dispatcher.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

# kind -> name -> {"obj": ..., "capabilities": set[str]}
_REGISTRY: dict[str, dict[str, dict[str, Any]]] = {}


class RegistryError(KeyError):
    pass


def register(kind: str, name: str, *, capabilities: Iterable[str] | None = None):
    """Decorator: register ``obj`` under ``(kind, name)``."""

    def decorator(obj):
        bucket = _REGISTRY.setdefault(kind, {})
        if name in bucket:
            raise RegistryError(f"{kind} {name!r} is already registered")
        bucket[name] = {"obj": obj, "capabilities": set(capabilities or ())}
        return obj

    return decorator


def get(kind: str, name: str):
    try:
        return _REGISTRY[kind][name]["obj"]
    except KeyError:
        raise RegistryError(
            f"unknown {kind} {name!r}; available: {available(kind)}"
        ) from None


def capabilities(kind: str, name: str) -> set[str]:
    get(kind, name)  # raises if missing
    return set(_REGISTRY[kind][name]["capabilities"])


def available(kind: str) -> list[str]:
    return sorted(_REGISTRY.get(kind, {}))
