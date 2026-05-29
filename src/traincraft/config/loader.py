"""Load and validate a TOML config into :class:`TrainCraftConfig`."""

from __future__ import annotations

import json
from pathlib import Path

import tomlkit

from .models import TrainCraftConfig


def load_config(path: str | Path) -> TrainCraftConfig:
    """Parse a TOML file and validate it (fail-fast with clear errors).

    Relative *input* paths (e.g. ``[geometry.source] path = ...``) are
    resolved relative to the config file's directory so that configs are
    portable regardless of the working directory from which the CLI is run.

    Relative *output* paths (``[run] outdir``) are left as-is so that output
    lands next to where the user runs the command, not next to the config file.
    """
    config_path = Path(path).resolve()
    text = config_path.read_text()
    cfg = loads_config(text)
    return _resolve_input_paths(cfg, config_path.parent)


def loads_config(text: str) -> TrainCraftConfig:
    """Parse a TOML string; does NOT resolve relative paths."""
    doc = tomlkit.parse(text)
    # Unwrap tomlkit's typed containers into plain Python for pydantic.
    plain = json.loads(json.dumps(doc))
    return TrainCraftConfig.model_validate(plain)


def _resolve_input_paths(cfg: TrainCraftConfig, base: Path) -> TrainCraftConfig:
    """Resolve relative input file paths against *base* (the config file's dir)."""
    if cfg.geometry is None:
        return cfg

    updates: dict = {}

    # source.path (e.g. FileSource)
    if cfg.geometry.source is not None and hasattr(cfg.geometry.source, "path"):
        src = cfg.geometry.source
        p = Path(src.path)
        if not p.is_absolute():
            updates["source"] = src.model_copy(update={"path": str(base / p)})

    # builder.file (e.g. SurfaceAdsorbateBuilder, SurfacePackingBuilder)
    if cfg.geometry.builder is not None and hasattr(cfg.geometry.builder, "file"):
        builder = cfg.geometry.builder
        if builder.file is not None:
            p = Path(builder.file)
            if not p.is_absolute():
                updates["builder"] = builder.model_copy(update={"file": str(base / p)})

    if not updates:
        return cfg
    new_geom = cfg.geometry.model_copy(update=updates)
    return cfg.model_copy(update={"geometry": new_geom})


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
