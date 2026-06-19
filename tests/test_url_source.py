from __future__ import annotations

from traincraft.config.models import GeometryConfig, UrlSource
from traincraft.geometry import build_geometry


def test_url_source_reads_local_file(tmp_path):
    xyz = tmp_path / "co.xyz"
    xyz.write_text("2\nCO\nC 0.0 0.0 0.0\nO 0.0 0.0 1.13\n")
    url = xyz.as_uri()  # file:// URL — exercises the download path, no network

    s = build_geometry(GeometryConfig(source=UrlSource(url=url)))
    assert len(s.atoms) == 2
    assert s.atoms.get_chemical_symbols() == ["C", "O"]
    assert s.provenance.source.startswith("source:url:")


def test_url_source_format_override(tmp_path):
    data = tmp_path / "co.data"  # non-standard suffix -> needs explicit format
    data.write_text("2\nCO\nC 0.0 0.0 0.0\nO 0.0 0.0 1.13\n")
    url = data.as_uri()

    s = build_geometry(GeometryConfig(source=UrlSource(url=url, format="xyz")))
    assert len(s.atoms) == 2
