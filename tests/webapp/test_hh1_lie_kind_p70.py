"""Plan 70 / R11 -- the ``hh1_lie`` algebra-only compute kind, served end-to-end.

Cross-runner CONTRACT (unmarked, per the Plan-32 extras-gated-dir ruling: no oracle_*
in tests/webapp): the wheel's hpc.spec core and the Pyodide GUI twin emit the block
byte-for-byte identical, so the browser cannot drift from the server / cluster run."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")


def _kron_body(field):
    return {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1, 2],
                        "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
                        "field": field},
            "compute": ["hh1_lie"],
            "artifacts": {"pdf": False, "tikz": False}}


def _truncpoly_body(a, field):
    return {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1],
                        "arrows": {"x": [1, 1]}, "relations": ["*".join(["x"] * a)],
                        "field": field},
            "compute": ["hh1_lie"],
            "artifacts": {"pdf": False, "tikz": False}}


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["hh1_lie"]


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_hh1lie", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("hh1_lie"))
    assert out["ok"], out
    return out["block"]


def test_kronecker_block_is_sl2(tmp_path):
    b = _server_block(_kron_body({"kind": "QQ"}), tmp_path)
    assert b["dim"] == 3 and b["solvable"] is False and b["perfect"] is True
    assert b["sl2_count"] == 1 and b["levi_type"] == "A1" and b["toral_rank"] == 1
    assert "rss_hh1_lie" in b["references"]


def test_truncpoly_qq_block_solvable(tmp_path):
    b = _server_block(_truncpoly_body(3, {"kind": "QQ"}), tmp_path)
    assert b["solvable"] is True and b["derived_series_dims"] == [2, 1, 0]


def test_default_field_cc_block_computes_with_base_change_note(tmp_path):
    """MAJOR-fix pin at the block layer: the DEFAULT field CC is FORMALLY closed
    (algebraically_closed True) but computes in exact QQ -- the block computes, the
    values equal the QQ run, and the base_change_note is present."""
    b = _server_block(_truncpoly_body(3, {"kind": "CC"}), tmp_path)
    assert b["algebraically_closed"] is True and b["dim"] == 2
    assert b["solvable"] is True and b["base_change_note"]


def test_char_p_block_gates_classification(tmp_path):
    b = _server_block(_truncpoly_body(3, {"kind": "GF", "p": 3, "n": 1}), tmp_path)
    assert b["dim"] == 3 and b["solvable"] is False        # W_1, field-general still computed
    assert b["radical_dim"] is None and b["char0_note"]     # loud char-p gate


def test_block_matches_across_runners_kronecker(tmp_path):
    body = _kron_body({"kind": "QQ"})
    assert (json.dumps(_gui_block(body), sort_keys=True, default=str)
            == json.dumps(_server_block(body, tmp_path), sort_keys=True, default=str))


def test_block_matches_across_runners_char_p(tmp_path):
    body = _truncpoly_body(3, {"kind": "GF", "p": 3, "n": 1})
    assert (json.dumps(_gui_block(body), sort_keys=True, default=str)
            == json.dumps(_server_block(body, tmp_path), sort_keys=True, default=str))


def test_block_matches_across_runners_default_cc(tmp_path):
    body = _truncpoly_body(3, {"kind": "CC"})
    assert (json.dumps(_gui_block(body), sort_keys=True, default=str)
            == json.dumps(_server_block(body, tmp_path), sort_keys=True, default=str))
