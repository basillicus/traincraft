"""Cheap calculators used for exploration.

These are *not* all MLIPs: EMT is a force field, tblite/xtb are semiempirical,
MACE is the ML potential. They share one factory interface. Heavy imports are
deferred so the package imports with no optional deps installed. (ANI and NEP
from the legacy code are intentionally dropped.)
"""

from __future__ import annotations

from ..core import register

_EF = {"energy", "forces"}


@register("calculator", "emt", capabilities=_EF)
def build_emt(cfg):
    """ASE EMT — zero extra deps; default for the walking-skeleton and CI."""
    from ase.calculators.emt import EMT

    return EMT()


@register("calculator", "tblite", capabilities=_EF | {"stress"})
def build_tblite(cfg):
    from tblite.ase import TBLite

    return TBLite(method=cfg.method)


@register("calculator", "xtb", capabilities=_EF)
def build_xtb(cfg):
    from xtb.ase.calculator import XTB

    return XTB(method=cfg.method)


@register("calculator", "mace", capabilities=_EF | {"stress"})
def build_mace(cfg):
    """MACE foundation models or a local fine-tuned model.

    Fixes the legacy plumbing: a local ``model_path`` is honoured explicitly,
    and ``mace-off23`` / ``mace-mp0`` map to the right loaders.
    """
    if cfg.model_path:
        from mace.calculators import MACECalculator

        return MACECalculator(
            model_paths=cfg.model_path,
            device=cfg.device,
            default_dtype=cfg.default_dtype,
        )
    if cfg.model == "mace-off23":
        from mace.calculators import mace_off

        return mace_off(model="large", device=cfg.device, default_dtype=cfg.default_dtype)
    # default: mace-mp0
    from mace.calculators import mace_mp

    return mace_mp(model="large", device=cfg.device, default_dtype=cfg.default_dtype)
