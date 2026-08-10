"""Plan 71 / R12 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin)
emit the ``hh_lie_module`` block byte-for-byte identical, so the browser cannot drift
from the server report / cluster run. A top-carrying HH kind (the ``bracket`` precedent);
top defaults to 2 and rides in the compute string ('hh_lie_module' or 'hh_lie_module:0..2')."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

# kK2 = sl2 (char 0, has weights) ; k[x]/x^3 over GF(3) = W_1 (char-p gate -> weights None).
_KRON_QQ = {"schema": 2,
            "algebra": {"kind": "quiver", "vertices": [1, 2],
                        "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
                        "field": {"kind": "QQ"}},
            "compute": ["hh_lie_module:0..2"],
            "artifacts": {"pdf": False, "tikz": False}}

_WITT_GF3 = {"schema": 2,
             "algebra": {"kind": "quiver", "vertices": [1],
                         "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                         "field": {"kind": "GF", "p": 3, "n": 1}},
             "compute": ["hh_lie_module:0..2"],
             "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body, spec):
    s = importlib.util.spec_from_file_location("gui_runner_hhliemodule_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(spec))
    assert out["ok"], out
    return out["block"]


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["hh_lie_module"]


def _assert_twins(body, tmp_path):
    gui = _gui_block(body, body["compute"][0])
    ref = _server_block(body, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
    return ref


def test_hh_lie_module_twin_kronecker(tmp_path):
    b = _assert_twins(_KRON_QQ, tmp_path)
    assert b["kind"] == "hh_lie_module" and b["hh_dims"] == [1, 3, 0] and b["hh1_dim"] == 3
    assert b["module_axiom_ok"] is True and b["inner_acts_zero"] is True
    # degree-1 = one irreducible summand of dim 3 (the sl2 adjoint L(2))
    deg1 = next(e for e in b["summands"] if e["n"] == 1)
    assert len(deg1["parts"]) == 1 and deg1["parts"][0]["dim"] == 3
    # degree-1 weight table is a symmetric sl2-string (three weights, symmetric about 0)
    w1 = next(e for e in b["weights"] if e["n"] == 1)
    vals = sorted(int(pair[0][0]) for pair in w1["weights"])
    assert len(vals) == 3 and vals[0] == -vals[2] and vals[1] == 0
    assert "alsolotar_toupie" in b["references"]


def test_hh_lie_module_twin_char_p_gate(tmp_path):
    b = _assert_twins(_WITT_GF3, tmp_path)
    assert b["weights"] is None and b["char0_note"]        # weight/torus char-0 gated
    assert b["module_axiom_ok"] is True                    # field-general part still computed
    # summands governed INDEPENDENTLY by decompose's guard; char 3 <= d_n=2 -> per-degree note
    assert b["summands"] is not None
