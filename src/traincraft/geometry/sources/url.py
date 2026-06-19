"""Source: download a structure file from a URL, then read it with ASE.

Useful for pulling a reference geometry straight from a dataset, a Gist, or a
provider's raw-file endpoint without committing it to the repo. The download
goes to a temp file; nothing is left behind.
"""

from __future__ import annotations

import tempfile
import urllib.request
from pathlib import Path

from ase.io import read

from ...core import Provenance, Structure, register


@register("source", "url")
def source_url(cfg) -> Structure:
    suffix = Path(cfg.url.split("?", 1)[0]).suffix or ".xyz"
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        req = urllib.request.Request(cfg.url, headers={"User-Agent": "traincraft"})
        with urllib.request.urlopen(req, timeout=cfg.timeout) as resp:
            tmp_path.write_bytes(resp.read())
        read_kwargs = {"format": cfg.format} if cfg.format else {}
        atoms = read(str(tmp_path), **read_kwargs)
    finally:
        tmp_path.unlink(missing_ok=True)

    return Structure.from_ase(
        atoms,
        provenance=Provenance(origin="generated", source=f"source:url:{cfg.url}"),
    )
