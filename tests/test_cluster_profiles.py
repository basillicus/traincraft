"""Cluster-profile registry: ``[orchestration.slurm] profile = "<name>"``.

A saved profile (``~/.traincraft/clusters/<name>.toml``, dir overridable via
``TRAINCRAFT_CLUSTERS_DIR``) is loaded as the base; inline ``[orchestration.slurm]``
keys override it, so a workflow moves between clusters by changing one line.
"""

from __future__ import annotations

import pytest

from traincraft.config import loads_config


@pytest.fixture()
def clusters_dir(tmp_path, monkeypatch):
    d = tmp_path / "clusters"
    d.mkdir()
    monkeypatch.setenv("TRAINCRAFT_CLUSTERS_DIR", str(d))
    return d


_LEONARDO = """\
account = "EUHPC_xxx"
runtime = "apptainer"
mpi = "pmix"
sif_dir = "/work/sif"
modules = ["apptainer"]
binds = ["/scratch", "/work"]

[stages.sample]
partition = "boost_usr_prod"
gpus = 1

[stages.label]
partition = "dcgp_usr_prod"
nodes = 2
ntasks = 224
"""


def test_profile_is_loaded_as_base(clusters_dir):
    (clusters_dir / "leonardo.toml").write_text(_LEONARDO)
    cfg = loads_config(
        '[orchestration]\nengine = "slurm"\n'
        '[orchestration.slurm]\nprofile = "leonardo"\n'
    )
    slurm = cfg.orchestration.slurm
    assert slurm.profile == "leonardo"  # recorded for provenance
    assert slurm.account == "EUHPC_xxx"
    assert slurm.mpi == "pmix"
    assert slurm.stages["sample"].partition == "boost_usr_prod"
    assert slurm.stages["label"].ntasks == 224


def test_inline_keys_override_profile(clusters_dir):
    (clusters_dir / "leonardo.toml").write_text(_LEONARDO)
    cfg = loads_config(
        '[orchestration]\nengine = "slurm"\n'
        '[orchestration.slurm]\n'
        'profile = "leonardo"\n'
        'account = "OTHER_acct"\n'      # scalar override wins
        'mpi = "pmi2"\n'
    )
    slurm = cfg.orchestration.slurm
    assert slurm.account == "OTHER_acct"
    assert slurm.mpi == "pmi2"
    # untouched profile keys survive
    assert slurm.runtime == "apptainer"
    assert slurm.stages["label"].ntasks == 224


def test_inline_stage_deep_merges_over_profile_stage(clusters_dir):
    (clusters_dir / "leonardo.toml").write_text(_LEONARDO)
    cfg = loads_config(
        '[orchestration]\nengine = "slurm"\n'
        '[orchestration.slurm]\nprofile = "leonardo"\n'
        '[orchestration.slurm.stages.label]\nnodes = 4\n'  # override one field only
    )
    label = cfg.orchestration.slurm.stages["label"]
    assert label.nodes == 4                 # inline wins
    assert label.partition == "dcgp_usr_prod"  # profile field preserved
    assert label.ntasks == 224
    # a profile-only stage is still present
    assert cfg.orchestration.slurm.stages["sample"].gpus == 1


def test_missing_profile_raises_with_available_listed(clusters_dir):
    (clusters_dir / "leonardo.toml").write_text(_LEONARDO)
    with pytest.raises(FileNotFoundError, match="lumi.*leonardo"):
        loads_config(
            '[orchestration]\nengine = "slurm"\n'
            '[orchestration.slurm]\nprofile = "lumi"\n'
        )


def test_no_profile_is_a_noop(clusters_dir):
    cfg = loads_config(
        '[orchestration]\nengine = "slurm"\n'
        '[orchestration.slurm]\naccount = "inline_only"\nmpi = "cray_shasta"\n'
    )
    assert cfg.orchestration.slurm.profile is None
    assert cfg.orchestration.slurm.account == "inline_only"
    assert cfg.orchestration.slurm.mpi == "cray_shasta"
