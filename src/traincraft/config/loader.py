"""Load and validate a TOML config into :class:`TrainCraftConfig`."""

from __future__ import annotations

import json
import os
from pathlib import Path

import tomlkit

from .models import TrainCraftConfig

# Where saved cluster profiles live. Overridable (tests, shared sites) via
# ``TRAINCRAFT_CLUSTERS_DIR``.
_DEFAULT_CLUSTERS_DIR = "~/.traincraft/clusters"


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
    plain = _resolve_slurm_profile(plain)
    return TrainCraftConfig.model_validate(plain)


def _clusters_dir() -> Path:
    raw = os.environ.get("TRAINCRAFT_CLUSTERS_DIR", _DEFAULT_CLUSTERS_DIR)
    return Path(os.path.expandvars(os.path.expanduser(raw)))


def _load_profile(name: str) -> dict:
    """Read a saved cluster profile (a serialized ``[orchestration.slurm]`` block)."""
    path = _clusters_dir() / f"{name}.toml"
    if not path.is_file():
        available = sorted(p.stem for p in _clusters_dir().glob("*.toml"))
        hint = ", ".join(available) if available else "none found"
        raise FileNotFoundError(
            f"cluster profile {name!r} not found at {path} "
            f"(available in {_clusters_dir()}: {hint})"
        )
    return json.loads(json.dumps(tomlkit.parse(path.read_text())))


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge *override* onto *base*; override wins for scalars/lists."""
    out = dict(base)
    for key, val in override.items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], val)
        else:
            out[key] = val
    return out


def _resolve_slurm_profile(plain: dict) -> dict:
    """Expand ``[orchestration.slurm] profile = "<name>"`` against the registry.

    The named profile is the base; the inline ``[orchestration.slurm]`` keys are
    layered on top (inline wins), so one workflow targets a different cluster by
    changing only the profile name. No-op when no profile is referenced.
    """
    slurm = plain.get("orchestration", {}).get("slurm")
    if not isinstance(slurm, dict) or not slurm.get("profile"):
        return plain
    merged = _deep_merge(_load_profile(slurm["profile"]), slurm)
    plain["orchestration"]["slurm"] = merged
    return plain


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
