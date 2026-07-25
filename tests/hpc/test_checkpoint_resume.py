"""Deepen checkpoint-resume end-to-end (Plan 28). A run with a zero time budget
stops cleanly at a checkpoint (exit 75, no result written); a rerun resumes and
completes (exit 0), and the final HH_* equals an uninterrupted run.

Marked deep (mirrors tests/engine/test_deepen.py's bucket): it exercises the
minimal-A^e resolution driver, so it belongs on the heavy CI leg."""
import json

import pytest

from quiverlab.hpc import cli

pytestmark = pytest.mark.deep

_CFG = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["hh_homology:0..6"],
    "artifacts": {"pdf": False, "tikz": False},
}


def _cfg(tmp_path):
    p = tmp_path / "deep.json"
    p.write_text(json.dumps(_CFG), encoding="utf-8")
    return p


def test_checkpoint_stop_then_resume_equals_uninterrupted(tmp_path):
    cfg = _cfg(tmp_path)

    # 1) time-limit 0 -> the deepen path stops at a checkpoint after one degree.
    out = tmp_path / "result.json"
    ck = tmp_path / "ck"
    rc = cli.main(["run", str(cfg), "-o", str(out),
                   "--checkpoint-dir", str(ck), "--time-limit", "0"])
    assert rc == 75
    assert not out.exists(), "no result must be written on a clean checkpoint stop"
    assert (ck / "latest.txt").exists(), "deepen must have written a checkpoint"

    # 2) rerun (no limit) -> resume from the checkpoint to completion.
    rc = cli.main(["run", str(cfg), "-o", str(out), "--checkpoint-dir", str(ck)])
    assert rc == 0
    resumed = json.loads(out.read_text(encoding="utf-8"))["results"]["hh_homology"]["dims"]

    # 3) a fresh, uninterrupted run.
    out2 = tmp_path / "result2.json"
    rc = cli.main(["run", str(cfg), "-o", str(out2), "--checkpoint-dir", str(tmp_path / "ck2")])
    assert rc == 0
    fresh = json.loads(out2.read_text(encoding="utf-8"))["results"]["hh_homology"]["dims"]

    assert resumed == fresh == [3, 2, 2, 2, 2, 2, 2]


def test_deepen_block_carries_engine_and_references(tmp_path):
    out = tmp_path / "r.json"
    rc = cli.main(["run", str(_cfg(tmp_path)), "-o", str(out),
                   "--checkpoint-dir", str(tmp_path / "ck")])
    assert rc == 0
    block = json.loads(out.read_text(encoding="utf-8"))["results"]["hh_homology"]
    assert block["kind"] == "HH_"
    assert "deepen" in block["engine"]
    assert isinstance(block["references"], list)
