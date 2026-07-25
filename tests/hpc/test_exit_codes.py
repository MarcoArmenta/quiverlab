"""The CLI exit-code matrix (BSD sysexits): 0 ok, 64 usage, 65 data/config,
73 cannot-write, 75 checkpoint. (75 is exercised in test_checkpoint_resume.)"""
import json

import pytest

from quiverlab.hpc import cli


def _write(tmp_path, obj, name="c.json"):
    p = tmp_path / name
    p.write_text(json.dumps(obj), encoding="utf-8")
    return p


def test_ok_is_zero(tmp_path):
    cfg = {"schema": 1,
           "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                       "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
           "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}}
    rc = cli.main(["run", str(_write(tmp_path, cfg)), "-o", str(tmp_path / "o.json")])
    assert rc == 0


def test_bad_config_is_65(tmp_path):
    # empty vertices -> a validation refusal
    cfg = {"schema": 1,
           "algebra": {"kind": "quiver", "vertices": [], "arrows": {},
                       "relations": [], "field": {"kind": "GF", "p": 2}},
           "compute": ["cartan"], "artifacts": {}}
    rc = cli.main(["run", str(_write(tmp_path, cfg)), "-o", str(tmp_path / "o.json")])
    assert rc == 65


def test_bad_relation_is_65(tmp_path):
    cfg = {"schema": 1,
           "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                       "relations": ["x*x*nonsense"], "field": {"kind": "GF", "p": 2}},
           "compute": ["cartan"], "artifacts": {}}
    rc = cli.main(["run", str(_write(tmp_path, cfg)), "-o", str(tmp_path / "o.json")])
    assert rc == 65


def test_missing_config_file_is_65(tmp_path):
    rc = cli.main(["run", str(tmp_path / "nope.json"), "-o", str(tmp_path / "o.json")])
    assert rc == 65


def test_no_verb_is_usage_64():
    assert cli.main([]) == 64


def test_argparse_usage_error_is_64():
    # missing required positional -> argparse error -> SystemExit(64)
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["run"])
    assert excinfo.value.code == 64


def test_unwritable_render_output_is_73(tmp_path):
    # produce a valid result, then render into a non-existent parent directory.
    from quiverlab.hpc.spec import run as spec_run, RESULT_SCHEMA
    cfg = {"schema": 1,
           "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                       "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
           "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}}
    rpath = tmp_path / "result.json"
    result = spec_run(cfg, tmp_path, result_schema=RESULT_SCHEMA)
    rpath.write_text(json.dumps(result), encoding="utf-8")
    bad_out = tmp_path / "no_such_dir" / "deeper" / "report.txt"
    rc = cli.main(["render", str(rpath), "-o", str(bad_out), "--format", "txt"])
    assert rc == 73


def test_gui_verb_missing_dependency_is_clean(monkeypatch):
    # Simulate the offline GUI module being unavailable (no [web]) -> the verb
    # must fail cleanly with a message, never a traceback, and never start a server.
    import sys
    monkeypatch.setitem(sys.modules, "webapp.server.offline", None)
    rc = cli.main(["gui", "--no-open"])
    assert rc == 70


def test_future_result_schema_render_is_65(tmp_path):
    p = tmp_path / "future.json"
    p.write_text(json.dumps({"result_schema": 999, "results": {},
                             "quiverlab_version": "9.9.9"}), encoding="utf-8")
    rc = cli.main(["render", str(p), "-o", str(tmp_path / "r.txt"), "--format", "txt"])
    assert rc == 65
