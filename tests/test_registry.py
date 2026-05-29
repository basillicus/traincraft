from __future__ import annotations

import pytest

from traincraft.core import RegistryError, available, get, register


def test_register_and_get():
    @register("ttest_kind", "foo")
    def foo():
        return 1

    assert get("ttest_kind", "foo")() == 1
    assert "foo" in available("ttest_kind")


def test_duplicate_raises():
    @register("ttest_dup", "bar")
    def bar():
        return None

    with pytest.raises(RegistryError):

        @register("ttest_dup", "bar")
        def bar2():
            return None


def test_unknown_raises():
    with pytest.raises(RegistryError):
        get("ttest_missing", "nope")


def test_builtins_registered():
    assert "nanotube" in available("builder")
    assert "emt" in available("calculator")
    assert "md" in available("sampler")
    assert "physicality" in available("selector")
