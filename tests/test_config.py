from __future__ import annotations

import pytest
from pydantic import ValidationError

from traincraft import loads_config

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
