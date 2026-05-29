"""Explicit run/job directories.

Replaces the legacy ``os.chdir`` + relative-filename approach: every step is
handed an absolute directory and writes there. Nothing mutates the process CWD,
so steps are safe to run in parallel.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class Job:
    dir: Path

    @property
    def marker(self) -> Path:
        return self.dir / ".tc_done"

    def done(self) -> bool:
        return self.marker.exists()

    def mark_done(self) -> None:
        self.marker.write_text("ok\n")

    def path(self, *parts: str) -> Path:
        return self.dir.joinpath(*parts)


class Workspace:
    """Owns an absolute run directory and hands out sub-directories/jobs."""

    def __init__(self, root: str | Path):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)

    def subdir(self, *parts: str) -> Path:
        p = self.root.joinpath(*parts)
        p.mkdir(parents=True, exist_ok=True)
        return p

    def job(self, *parts: str) -> Job:
        return Job(dir=self.subdir(*parts))
