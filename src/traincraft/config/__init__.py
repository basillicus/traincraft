"""Configuration: pydantic models + TOML loader."""

from __future__ import annotations

from .loader import dump_starter_config, load_config, loads_config
from .models import TrainCraftConfig

__all__ = ["TrainCraftConfig", "dump_starter_config", "load_config", "loads_config"]
