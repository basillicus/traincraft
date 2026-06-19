"""Sources: pull structures from external databases.

  * ``materials_project`` — a bulk crystal by ``mp-xxxx`` id (needs ``mp-api`` +
    ``pymatgen``; API key from config or ``$MP_API_KEY``);
  * ``optimade``          — the first hit of an OPTIMADE filter against any
    provider base URL (dependency-free: parses the OPTIMADE JSON directly);
  * ``pubchem``           — a 3D conformer by name / CID / SMILES from the
    PubChem PUG REST API (dependency-free: downloads an SDF and reads it).
"""

from __future__ import annotations

import json
import os
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read

from ...core import Provenance, Structure, register


def _get(url: str, timeout: float) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "traincraft"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


@register("source", "materials_project")
def source_materials_project(cfg) -> Structure:
    api_key = cfg.api_key or os.environ.get("MP_API_KEY")
    if not api_key:
        raise ValueError(
            "materials_project source needs an API key (cfg.api_key or $MP_API_KEY)"
        )
    try:
        from mp_api.client import MPRester
        from pymatgen.io.ase import AseAtomsAdaptor
    except ImportError as e:
        raise ImportError(
            "materials_project source needs mp-api and pymatgen. "
            "Install with: pixi install -e science  (plus: uv pip install mp-api)"
        ) from e

    with MPRester(api_key) as mpr:
        pmg = mpr.get_structure_by_material_id(
            cfg.material_id, conventional_unit_cell=cfg.conventional
        )
    atoms = AseAtomsAdaptor.get_atoms(pmg)
    return Structure.from_ase(
        atoms,
        provenance=Provenance(
            origin="generated", source=f"source:materials_project:{cfg.material_id}"
        ),
    )


def _atoms_from_optimade(attrs: dict) -> Atoms:
    """Build an ase.Atoms from an OPTIMADE structure resource's attributes."""
    symbol_of = {
        s["name"]: s["chemical_symbols"][0]
        for s in attrs.get("species", [])
        if s.get("chemical_symbols")
    }
    symbols = [symbol_of.get(name, name) for name in attrs["species_at_sites"]]
    positions = np.asarray(attrs["cartesian_site_positions"], dtype=float)
    cell = attrs.get("lattice_vectors")
    pbc = [bool(d) for d in attrs.get("dimension_types", [0, 0, 0])]
    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=pbc)


@register("source", "optimade")
def source_optimade(cfg) -> Structure:
    query = {"page_limit": 1, "response_format": "json"}
    if cfg.filter:
        query["filter"] = cfg.filter
    url = f"{cfg.base_url.rstrip('/')}/structures?{urllib.parse.urlencode(query)}"
    payload = json.loads(_get(url, cfg.timeout))
    data = payload.get("data") or []
    if not data:
        raise RuntimeError(f"optimade query returned no structures: {cfg.filter!r}")
    atoms = _atoms_from_optimade(data[0]["attributes"])
    return Structure.from_ase(
        atoms,
        provenance=Provenance(origin="generated", source=f"source:optimade:{cfg.base_url}"),
    )


@register("source", "pubchem")
def source_pubchem(cfg) -> Structure:
    namespace = "cid" if cfg.cid is not None else ("smiles" if cfg.smiles else "name")
    ident = cfg.cid if cfg.cid is not None else (cfg.smiles or cfg.name)
    if ident is None:
        raise ValueError("pubchem source needs one of 'cid', 'name', or 'smiles'")

    ident_q = urllib.parse.quote(str(ident), safe="")
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        f"{namespace}/{ident_q}/record/SDF?record_type=3d"
    )
    sdf = _get(url, cfg.timeout)
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tmp:
        tmp_path = Path(tmp.name)
        tmp_path.write_bytes(sdf)
    try:
        atoms = read(str(tmp_path), format="sdf")
    finally:
        tmp_path.unlink(missing_ok=True)

    return Structure.from_ase(
        atoms,
        provenance=Provenance(origin="generated", source=f"source:pubchem:{namespace}:{ident}"),
    )
