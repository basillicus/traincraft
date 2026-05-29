"""Selector: remove exact duplicates by content hash.

Descriptor-based near-duplicate removal lands in a later chunk.
"""

from __future__ import annotations

from ..core import Structure, register


@register("selector", "dedup")
def select_dedup(frames: list[Structure], cfg) -> list[Structure]:
    seen: set[str] = set()
    out: list[Structure] = []
    for s in frames:
        h = s.hash
        if h in seen:
            continue
        seen.add(h)
        out.append(s)
    return out
