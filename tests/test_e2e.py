from __future__ import annotations

from pathlib import Path

from traincraft import load_config, run_pipeline
from traincraft.datasets import read_frames

ROOT = Path(__file__).resolve().parents[1]


def test_walking_skeleton(tmp_path):
    cfg = load_config(ROOT / "examples" / "01_cnt_emt_md.toml")
    cfg.run.outdir = str(tmp_path)
    summary = run_pipeline(cfg)

    assert summary["n_candidates"] >= 1
    assert summary["n_selected"] >= 1

    dataset = Path(summary["dataset"])
    assert dataset.exists()

    frames = read_frames(dataset)
    assert 1 <= len(frames) <= 3
    assert frames[0].provenance.origin in ("ml_sampled", "generated")
    # selected frames came from the sampler and carry forces
    assert "forces" in frames[0].properties
