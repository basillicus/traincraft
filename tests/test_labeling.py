from __future__ import annotations

import json

from traincraft import label_frames, read_frames, run_pipeline
from traincraft.config.models import (
    CrystalBuilder,
    DatasetConfig,
    EmtCalc,
    GeometryConfig,
    LabelingConfig,
    MdSampling,
    RunConfig,
    SelectionConfig,
    TrainCraftConfig,
)
from traincraft.geometry import build_geometry

# EMT stands in for a DFT calculator here: it produces energy/forces/stress with
# zero deps, so the labeling *mechanics* are exercised without a DFT binary.


def _cu_frames(n=2):
    base = build_geometry(
        GeometryConfig(builder=CrystalBuilder(name="Cu", crystalstructure="fcc", a=3.6, cubic=True))
    )
    frames = []
    for i in range(n):
        s = base.copy()
        s.atoms.rattle(0.01, seed=i)
        frames.append(s)
    return frames


def test_label_frames_sets_properties_and_origin():
    labeled = label_frames(_cu_frames(2), EmtCalc())
    assert len(labeled) == 2
    for s in labeled:
        assert s.provenance.origin == "dft_labeled"
        assert s.provenance.calculator == "emt"
        assert s.provenance.level_of_theory["calculator"] == "emt"
        assert "energy" in s.properties
        assert "forces" in s.properties
        assert "stress" in s.properties  # Cu bulk is periodic


def test_label_frames_writes_manifest_and_frame_dirs(tmp_path):
    label_frames(_cu_frames(2), EmtCalc(), out_dir=tmp_path)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["origin"] == "dft_labeled"
    assert manifest["calculator"] == "emt"
    assert manifest["n_frames"] == 2
    assert "energy" in manifest["properties"]
    assert (tmp_path / "frame_0000").is_dir()


def test_pipeline_labels_selected_frames(tmp_path):
    cfg = TrainCraftConfig(
        run=RunConfig(name="lbl", outdir=str(tmp_path), seed=1),
        geometry=GeometryConfig(
            builder=CrystalBuilder(
                name="Cu", crystalstructure="fcc", a=3.6, cubic=True, supercell=(2, 2, 2)
            )
        ),
        calculator=EmtCalc(),
        sampling=MdSampling(temperature=300.0, steps=50, interval=10),
        selection=SelectionConfig(budget=3),
        labeling=LabelingConfig(calculator=EmtCalc()),
        dataset=DatasetConfig(),
    )
    summary = run_pipeline(cfg)
    assert summary["n_labeled"] > 0
    frames = read_frames(summary["dataset"])
    assert len(frames) > 0
    assert all(f.provenance.origin == "dft_labeled" for f in frames)
    assert all("energy" in f.properties for f in frames)
