from __future__ import annotations

import os

from traincraft.core import Workspace


def test_absolute_dirs_and_no_chdir(tmp_path):
    cwd = os.getcwd()
    ws = Workspace(tmp_path / "run")
    j1 = ws.job("a")
    j2 = ws.job("b")
    assert j1.dir.is_absolute() and j2.dir.is_absolute()
    assert j1.dir != j2.dir
    assert os.getcwd() == cwd  # nothing changed the process CWD


def test_done_marker(tmp_path):
    ws = Workspace(tmp_path / "run")
    job = ws.job("step")
    assert not job.done()
    job.mark_done()
    assert job.done()
