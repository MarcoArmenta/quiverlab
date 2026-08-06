"""The tilting_check module kind (Plan 44 / C7): served by quiverlab.hpc.spec (server /
container / CLI), mirrored byte-for-byte on the math subkeys by the Pyodide twin
docs/gui/runner.py. Schema-2 (the candidate T is the request module block; no second
module). No oracle-class marker (tests/webapp collects only under the [web] extra)."""
import json
import pathlib
import tempfile

import pytest

from quiverlab.hpc import spec

_A2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
       "relations": [], "field": {"kind": "GF", "p": 5}}   # char 5 > dim: decompose rigorous
_MATH_KEYS = ("kind", "is_tilting", "n", "pd", "self_ext_vanishes",
              "num_summands", "num_vertices", "note")


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    import importlib.util
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p44", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


def _req(builtin, compute="tilting_check"):
    kind, vertex = builtin
    return {"schema": 2, "algebra": _A2,
            "module": {"builtin": {"kind": kind, "vertex": vertex}},
            "compute": [compute]}


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_tilting_check_block_shape(run):
    # module = P1 over kA2: only ONE indecomposable summand for TWO vertices -> not tilting.
    b = run(_req(("projective", 1)))["tilting_check"]
    assert b["is_tilting"] is False
    assert b["num_vertices"] == 2 and b["num_summands"] == 1
    assert b["pd"] == 0 and b["self_ext_vanishes"] is True
    assert "citations" in b


def test_twin_parity(tmp_path):
    # the SAME request through both runners agrees on every math subkey (byte-identical
    # via json.dumps(sort_keys=True)); the tie is: the draw page and the CLI/report agree.
    req = _req(("projective", 1))
    s = _server(req)["tilting_check"]
    p = _pyodide(req)["tilting_check"]
    sub_s = {k: s[k] for k in _MATH_KEYS}
    sub_p = {k: p[k] for k in _MATH_KEYS}
    assert json.dumps(sub_s, sort_keys=True) == json.dumps(sub_p, sort_keys=True)
