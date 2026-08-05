"""derived_fingerprint scalar kind (Plan 43): schema-v1, served by hpc.spec,
mirrored byte-identically by the Pyodide twin (docs/gui/runner.py).

NO oracle-class marker: tests/webapp/ collects only with the [web] extra, and the
Plan-32 audit requires the audited class counts to be environment-independent.
These are contract tests for the two runners' block shapes."""
import json
import tempfile

import pytest

from quiverlab.hpc import spec

_A4 = {"kind": "quiver", "vertices": [1, 2, 3, 4],
       "arrows": {"a": [1, 2], "b": [2, 3], "c": [3, 4]},
       "relations": [], "field": {"kind": "GF", "p": 5}}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    import importlib.util
    import pathlib
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p43", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


def test_derived_fingerprint_block_shape():
    req = {"schema": 1, "algebra": _A4, "compute": ["derived_fingerprint"]}
    block = _server(req)["derived_fingerprint"]
    assert block["kind"] == "derived_fingerprint"
    assert block["fingerprint"]["coxeter_polynomial"] == "t**4 + t**3 + t**2 + t + 1"
    assert block["fingerprint"]["cartan_det"] == 1
    # references are the registry keys; citations resolve them to (bibkey, formatted).
    assert "happel_triangulated" in block["references"]
    assert "Happel1988" in [k for k, _ in block["citations"]]


def test_twin_parity():
    req = {"schema": 1, "algebra": _A4, "compute": ["derived_fingerprint"]}
    srv = _server(req)["derived_fingerprint"]
    pyo = _pyodide(req)["derived_fingerprint"]
    assert (json.dumps(srv, sort_keys=True, default=str)
            == json.dumps(pyo, sort_keys=True, default=str))
