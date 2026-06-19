"""Builder: a bare crystalline slab/surface (no adsorbate).

Two ways to specify the surface:
  * a named low-index facet (``facet="fcc111"``) — quick, element-only;
  * arbitrary Miller indices cleaved from a bulk crystal (``miller=[2,1,1]``)
    with an explicit ``crystalstructure`` and lattice constants.

The slab is periodic in-plane and (by default) non-periodic along the surface
normal, with vacuum padding. All atoms are framework (``tc_fragment == -1``):
a bare slab has no mobile fragments until something is adsorbed onto it.
"""

from __future__ import annotations

import numpy as np

from ...core import FRAMEWORK, Provenance, Structure, register, set_fragments
from ._common import build_miller_slab, build_named_slab


@register("builder", "slab")
def build_slab(cfg) -> Structure:
    if cfg.facet is not None:
        slab = build_named_slab(
            cfg.element, cfg.facet, cfg.size, cfg.vacuum, orthogonal=cfg.orthogonal
        )
        tag = f"{cfg.element}-{cfg.facet}"
    else:
        slab = build_miller_slab(
            cfg.element,
            cfg.crystalstructure,
            tuple(cfg.miller),
            cfg.layers,
            cfg.vacuum,
            a=cfg.a,
            c=cfg.c,
            cubic=cfg.cubic,
            size=(cfg.size[0], cfg.size[1]),
            periodic=cfg.periodic,
        )
        miller_str = "".join(str(i) for i in cfg.miller)
        tag = f"{cfg.element}-{miller_str}"

    # A bare slab is all framework; nothing is mobile yet.
    set_fragments(slab, np.full(len(slab), FRAMEWORK, dtype=int))

    return Structure.from_ase(
        slab,
        provenance=Provenance(origin="generated", source=f"builder:slab:{tag}"),
    )
