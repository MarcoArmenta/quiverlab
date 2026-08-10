"""Plan 73 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin) emit
the ``han_transport`` block byte-for-byte identical, so the browser cannot drift from
the server report / cluster run. An ALGEBRA-level CERTIFICATE kind carrying the
new-arrow subset ``new_arrows`` (the derived_compare second-input precedent). Unmarked
per the Plan-32 extras-gated ruling (a cross-runner contract, not an oracle)."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

# CLMS Ex. 5.3 (bounded), F = {a}; GF(32003) > dim A/B = 27.
_EX53 = {"schema": 1,
         "algebra": {"kind": "quiver", "vertices": [1, 2, 3, 4, 5],
                     "arrows": {"al": [5, 1], "be": [1, 5], "d": [1, 4],
                                "c": [4, 3], "b": [3, 2], "a": [2, 1]},
                     "relations": ["al*be", "d*c*b*a - be*al"],
                     "field": {"kind": "GF", "p": 32003, "n": 1}},
         "compute": ["han_transport"],
         "new_arrows": ["a"],
         "artifacts": {"pdf": False, "tikz": False}}

# CLMS Ex. 5.5 (NOT bounded), F = {d}.
_EX55 = {"schema": 1,
         "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                     "arrows": {"a": [1, 2], "b": [2, 3], "c": [1, 3], "d": [3, 1]},
                     "relations": ["a*b", "d*a", "b*d", "d*c*d"],
                     "field": {"kind": "GF", "p": 32003, "n": 1}},
         "compute": ["han_transport"],
         "new_arrows": ["d"],
         "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body, compute):
    spec = importlib.util.spec_from_file_location("gui_runner_p73", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(compute))
    assert out["ok"], out
    return out["block"]


def _server_block(body, kind, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"][kind]


def test_han_transport_bounded_matches_across_runners(tmp_path):
    gui = _gui_block(_EX53, "han_transport")
    ref = _server_block(_EX53, "han_transport", tmp_path)
    assert gui["transport"] == "bounded"
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))


def test_han_transport_not_bounded_matches_across_runners(tmp_path):
    gui = _gui_block(_EX55, "han_transport")
    ref = _server_block(_EX55, "han_transport", tmp_path)
    assert gui["transport"] == "not_bounded"
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
