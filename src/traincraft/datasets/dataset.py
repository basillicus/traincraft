"""A simple, deduplicated, provenance-aware dataset."""

from __future__ import annotations

from pathlib import Path

from ..core import Structure
from .io import read_frames, write_frames


class Dataset:
    def __init__(self, path: str | Path):
        path = Path(path)
        if path.suffix != ".extxyz":
            path = path.with_suffix(".extxyz")
        self.path = path
        self._frames: list[Structure] = []
        self._hashes: set[str] = set()

    def append(self, structures: list[Structure]) -> int:
        """Add new frames, skipping exact duplicates. Returns count added."""
        added = 0
        for s in structures:
            h = s.hash
            if h in self._hashes:
                continue
            self._hashes.add(h)
            self._frames.append(s)
            added += 1
        return added

    def filter(self, origin: str | None = None) -> list[Structure]:
        if origin is None:
            return list(self._frames)
        return [s for s in self._frames if s.provenance.origin == origin]

    @property
    def frames(self) -> list[Structure]:
        return list(self._frames)

    def __len__(self) -> int:
        return len(self._frames)

    def write(self) -> Path:
        return write_frames(self.path, self._frames)

    @classmethod
    def load(cls, path: str | Path) -> Dataset:
        ds = cls(path)
        ds.append(read_frames(ds.path))
        return ds
