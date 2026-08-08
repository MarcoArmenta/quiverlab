"""Plan 70 / R11 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin)
emit the ``hh1_lie`` block byte-for-byte identical, so the browser cannot drift from the
server report / cluster run. Algebra-level DIM-budget kind (no module block); the budget
defaults to 48 and rides in the compute string ('hh1_lie' or 'hh1_lie:48')."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

# kK2 = sl2 (char 0), k[x]/x^3 over GF(3) = W_1 (char-p gate), k[x]/x^3 default CC.
_KRON_QQ = {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1, 2],
                        "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
                        "field": {"kind": "QQ"}},
            "compute": ["hh1_lie"],
            "artifacts": {"pdf": False, "tikz": False}}

_WITT_GF3 = {"schema": 1,
             "algebra": {"kind": "quiver", "vertices": [1],
                         "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                         "field": {"kind": "GF", "p": 3, "n": 1}},
             "compute": ["hh1_lie"],
             "artifacts": {"pdf": False, "tikz": False}}

_CC_DEFAULT = {"schema": 1,
               "algebra": {"kind": "quiver", "vertices": [1],
                           "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                           "field": {"kind": "CC"}},
               "compute": ["hh1_lie"],
               "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_hh1lie_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("hh1_lie"))
    assert out["ok"], out
    return out["block"]


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["hh1_lie"]


def _assert_twins(body, tmp_path):
    gui = _gui_block(body)
    ref = _server_block(body, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
    return ref


def test_hh1_lie_twin_kronecker(tmp_path):
    b = _assert_twins(_KRON_QQ, tmp_path)
    assert b["kind"] == "hh1_lie" and b["dim"] == 3 and b["sl2_count"] == 1


def test_hh1_lie_twin_witt_char_p(tmp_path):
    b = _assert_twins(_WITT_GF3, tmp_path)
    assert b["radical_dim"] is None and b["char0_note"]


def test_hh1_lie_twin_default_cc(tmp_path):
    b = _assert_twins(_CC_DEFAULT, tmp_path)
    assert b["algebraically_closed"] is True and b["base_change_note"]
