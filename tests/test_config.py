from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from traincraft import loads_config
from traincraft.config.loader import load_config

GOOD = """
[run]
name = "t"
[geometry.builder]
type = "molecule"
name = "H2O"
[calculator]
type = "emt"
"""


def test_loads_good():
    cfg = loads_config(GOOD)
    assert cfg.run.name == "t"
    assert cfg.geometry.builder.type == "molecule"
    assert cfg.calculator.type == "emt"


def test_unknown_key_rejected():
    with pytest.raises(ValidationError):
        loads_config('[run]\nnam = "x"\n')


def test_unknown_builder_type_rejected():
    with pytest.raises(ValidationError):
        loads_config('[geometry.builder]\ntype = "wat"\n')


def test_geometry_needs_source_or_builder():
    with pytest.raises(ValidationError):
        loads_config("[geometry]\ntransforms = []\n")


def test_file_source_path_resolved_relative_to_config(tmp_path):
    """Relative paths in [geometry.source] resolve relative to the config file,
    not relative to the working directory from which the CLI is run."""
    asset = tmp_path / "struct.xyz"
    asset.write_text("1\ntest\nH 0 0 0\n")

    config_text = '[geometry.source]\ntype = "file"\npath = "struct.xyz"\n'
    config_file = tmp_path / "run.toml"
    config_file.write_text(config_text)

    cfg = load_config(config_file)
    # Path must now be absolute and point at the asset next to the config file.
    assert Path(cfg.geometry.source.path).is_absolute()
    assert Path(cfg.geometry.source.path) == asset


def test_absolute_file_source_path_unchanged(tmp_path):
    asset = tmp_path / "struct.xyz"
    asset.write_text("1\ntest\nH 0 0 0\n")

    config_text = f'[geometry.source]\ntype = "file"\npath = "{asset}"\n'
    config_file = tmp_path / "run.toml"
    config_file.write_text(config_text)

    cfg = load_config(config_file)
    assert Path(cfg.geometry.source.path) == asset
