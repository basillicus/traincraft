from __future__ import annotations

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
from traincraft.orchestration import enabled_stages, run_stage, workspace_for


def _cfg(tmp_path, supercell=(2, 2, 2)):
    return TrainCraftConfig(
        run=RunConfig(name="staged", outdir=str(tmp_path), seed=1),
        geometry=GeometryConfig(
            builder=CrystalBuilder(
                name="Cu", crystalstructure="fcc", a=3.6, cubic=True, supercell=supercell
            )
        ),
        calculator=EmtCalc(),
        sampling=MdSampling(temperature=300.0, steps=50, interval=10),
        selection=SelectionConfig(budget=3),
        labeling=LabelingConfig(calculator=EmtCalc()),
        dataset=DatasetConfig(),
    )


def test_enabled_stages_lists_requested_sections(tmp_path):
    assert enabled_stages(_cfg(tmp_path)) == ["geometry", "sample", "select", "label", "dataset"]


def test_stages_chain_through_artifacts(tmp_path):
    cfg = _cfg(tmp_path)
    ws = workspace_for(cfg)

    assert len(run_stage("geometry", cfg, ws)) == 1
    assert (ws.root / "structures" / "initial.extxyz").exists()

    assert len(run_stage("sample", cfg, ws)) > 1
    assert (ws.root / "candidates" / "candidates.extxyz").exists()

    selected = run_stage("select", cfg, ws)
    assert len(selected) <= 3
    assert (ws.root / "selected" / "selected.extxyz").exists()

    labeled = run_stage("label", cfg, ws)
    assert (ws.root / "labeled_dft" / "labeled.extxyz").exists()
    assert (ws.root / "labeled_dft" / "manifest.json").exists()
    assert all(s.provenance.origin == "dft_labeled" for s in labeled)

    run_stage("dataset", cfg, ws)
    assert (ws.root / "dataset.extxyz").exists()


def test_stage_uses_cache_on_rerun(tmp_path):
    cfg = _cfg(tmp_path)
    ws = workspace_for(cfg)
    run_stage("geometry", cfg, ws)
    # Re-running reads the cached artifact rather than rebuilding; same result.
    again = run_stage("geometry", cfg, ws)
    assert len(again) == 1


def test_unchanged_config_does_not_rebuild(tmp_path, monkeypatch):
    import traincraft.orchestration.stages as st

    cfg = _cfg(tmp_path)
    ws = workspace_for(cfg)
    run_stage("geometry", cfg, ws)

    calls = {"n": 0}
    real = st.build_geometry

    def counting(geometry_cfg):
        calls["n"] += 1
        return real(geometry_cfg)

    monkeypatch.setattr(st, "build_geometry", counting)

    run_stage("geometry", cfg, ws)  # same config -> cached, no rebuild
    assert calls["n"] == 0
    run_stage("geometry", cfg, ws, force=True)  # --force -> rebuild
    assert calls["n"] == 1


def test_changed_geometry_config_invalidates(tmp_path):
    ws = workspace_for(_cfg(tmp_path))
    small = run_stage("geometry", _cfg(tmp_path, supercell=(2, 2, 2)), ws)[0]
    # Editing the builder changes the geometry key, so it rebuilds (bigger cell).
    big = run_stage("geometry", _cfg(tmp_path, supercell=(3, 3, 3)), ws)[0]
    assert len(big.atoms) != len(small.atoms)


def test_geometry_change_cascades_to_sample(tmp_path):
    cfg = _cfg(tmp_path)
    ws = workspace_for(cfg)
    run_stage("geometry", cfg, ws)
    n_small = len(run_stage("sample", cfg, ws)[0].atoms)

    # Change only the geometry; sample's upstream key changes, so it must recompute
    # from the new (bigger) structure even though [sampling] itself is unchanged.
    cfg2 = _cfg(tmp_path, supercell=(3, 3, 3))
    run_stage("geometry", cfg2, ws)
    cand = run_stage("sample", cfg2, ws)
    assert len(cand[0].atoms) != n_small
