"""The `string_homological` + `toupie` algebra-only kinds (Plan 59): served by
hpc.spec, mirrored by the Pyodide twin (docs/gui/runner.py), byte-identical blocks.
Unmarked (extras-gated dir): a cross-runner CONTRACT test, not an oracle."""
import importlib.util
import json
import pathlib
import tempfile

from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest

ROOT = pathlib.Path(__file__).resolve().parents[2]
_PYO = None

# kD4 (subspace orientation): rep-finite, NOT string -> the 3-summand-middle witness.
_KD4 = {"schema": 1,
        "algebra": {"kind": "quiver", "vertices": [0, 1, 2, 3],
                    "arrows": {"a": [1, 0], "b": [2, 0], "c": [3, 0]},
                    "relations": [], "field": {"kind": "QQ"}},
        "compute": ["string_homological"], "artifacts": {"pdf": False, "tikz": False}}
# the 3-Kronecker = ToupieAlgebra([1,1,1]) as a raw quiver: 3 parallel arrows 1 -> 2.
_KRON3 = {"schema": 1,
          "algebra": {"kind": "quiver", "vertices": [1, 2],
                      "arrows": {"a": [1, 2], "b": [1, 2], "c": [1, 2]},
                      "relations": [], "field": {"kind": "QQ"}},
          "compute": ["toupie"], "artifacts": {"pdf": False, "tikz": False}}


def _pyodide_runner():
    global _PYO
    if _PYO is None:
        path = ROOT / "docs" / "gui" / "runner.py"
        s = importlib.util.spec_from_file_location("_p59_gui_runner_twin", path)
        m = importlib.util.module_from_spec(s)
        s.loader.exec_module(m)
        _PYO = m
    return _PYO


def _server(body, kind):
    with tempfile.TemporaryDirectory() as d:
        res = run_spec(ComputeRequest.model_validate(body), d)
    assert "error" not in res, res.get("error")
    return res["results"][kind]


def _pyo(body, item):
    gui = _pyodide_runner()
    assert json.loads(gui.run_build(json.dumps(body)))["ok"], "Pyodide run_build failed"
    out = json.loads(gui.compute_one(item))
    assert out["ok"], out
    return out["block"]


def test_string_homological_block_kD4():
    block = _server(_KD4, "string_homological")
    assert block["is_string"] is False
    assert block["verdict"] == "not_string"
    assert block["witness"]["summand_count"] >= 3
    # `references` carries the registry keys; `citations` the resolved (bibtex_key, text).
    assert "suarez_alvarez" in block["references"]
    assert block["citations"], "citation pairs should resolve"


def test_toupie_block_a_kronecker():
    block = _server(_KRON3, "toupie")
    assert block["is_toupie"] is True
    assert block["branch_count"] == 3 and block["direct_arrow_count"] == 3
    assert block["hh"][:2] == [1, 8]
    assert block["sl_a"]["dim"] == 8 and block["sl_a"]["char0"] is True
    assert "alsolotar_toupie" in block["references"]
    assert block["citations"], "citation pairs should resolve"


def test_twin_parity():
    for body, kind in ((_KD4, "string_homological"), (_KRON3, "toupie")):
        server = json.dumps(_server(body, kind), sort_keys=True, default=str)
        gui = json.dumps(_pyo(body, kind), sort_keys=True, default=str)
        assert server == gui, f"runner twin drift for {kind}"
