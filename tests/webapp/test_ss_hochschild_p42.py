"""ss_hochschild scalar/range kind (Plan 42): served by hpc.spec, mirrored
byte-identically by the Pyodide twin (docs/gui/runner.py). Unmarked -- the
extras-gated webapp dir per the Plan-32 ruling."""
import importlib.util
import json
import pathlib
import tempfile

import pytest

from quiverlab import GF, truncated_polynomial
from quiverlab.hpc import spec

# k[x]/(x^2) over GF(5), schema v1 (algebra-only kind).
_DUAL = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
         "relations": ["x*x"], "field": {"kind": "GF", "p": 5}}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p42", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_block_shape(run):
    req = {"schema": 1, "algebra": _DUAL, "compute": ["ss_hochschild:0..4"]}
    b = run(req)["ss_hochschild"]
    assert b["abutment"] == list(
        truncated_polynomial(2, field=GF(5)).cyclic_homology(4).dims)
    assert b["abutment"] == [2, 0, 2, 0, 2]
    assert "grid" in b and "prose" in b
    assert b["collapse"] is True and b["degenerates_at"] == 2
    assert "weibel_homological" in b["references"]


def test_twin_parity():
    req = {"schema": 1, "algebra": _DUAL, "compute": ["ss_hochschild:0..4"]}
    srv, pyo = _server(req), _pyodide(req)
    assert set(srv) == set(pyo)
    for kind in srv:
        assert (json.dumps(srv[kind], sort_keys=True, default=str)
                == json.dumps(pyo[kind], sort_keys=True, default=str)), kind
