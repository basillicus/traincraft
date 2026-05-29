"""Load and validate a TOML config into :class:`TrainCraftConfig`."""

from __future__ import annotations

import json
from pathlib import Path

import tomlkit

from .models import TrainCraftConfig


def load_config(path: str | Path) -> TrainCraftConfig:
    """Parse a TOML file and validate it (fail-fast with clear errors)."""
    text = Path(path).read_text()
    return loads_config(text)


def loads_config(text: str) -> TrainCraftConfig:
    doc = tomlkit.parse(text)
    # Unwrap tomlkit's typed containers into plain python for pydantic.
    plain = json.loads(json.dumps(doc))
    return TrainCraftConfig.model_validate(plain)


def dump_starter_config() -> str:
    """A minimal, valid starter config for ``traincraft new``."""
    return (
        '[run]\n'
        'name = "my_run"\n'
        'outdir = "runs"\n'
        'seed = 42\n'
        '\n'
        '[geometry.builder]\n'
        'type = "molecule"\n'
        'name = "H2O"\n'
        '\n'
        '[calculator]\n'
        'type = "emt"\n'
        '\n'
        '[sampling]\n'
        'type = "md"\n'
        'temperature = 300.0\n'
        'steps = 100\n'
        'interval = 10\n'
        '\n'
        '[selection]\n'
        'steps = ["physicality", "dedup", "diversity"]\n'
        'budget = 5\n'
        '\n'
        '[dataset]\n'
        'path = "dataset"\n'
    )
