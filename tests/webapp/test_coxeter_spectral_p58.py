"""The coxeter_spectral compute kind across BOTH runners (Plan 58 / R20).

``quiverlab.hpc.spec`` (server / container / CLI) and its Pyodide twin
``docs/gui/runner.py`` must produce byte-identical blocks, or the draw page and the
desktop app would disagree about the same certified spectral report.

NO oracle-class marker: tests/webapp/ collects only with the [web] extra, and the
Plan-32 audit requires the audited class counts to be environment-independent. These
are cross-runner contract tests.
"""
import json
import tempfile

import pytest

from quiverlab.hpc import spec

_KA2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
        "relations": [], "field": {"kind": "QQ"}}
_KRON3 = {"kind": "quiver", "vertices": [1, 2],
          "arrows": {"a0": [1, 2], "a1": [1, 2], "a2": [1, 2]},
          "relations": [], "field": {"kind": "QQ"}}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    import importlib.util
    import pathlib
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p58", path)
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
def test_coxeter_spectral_block_shape(run):
    b = run({"schema": 1, "algebra": _KRON3, "compute": ["coxeter_spectral"]})[
        "coxeter_spectral"]
    assert b["kind"] == "coxeter_spectral"
    assert b["cyclotomic"] is False
    assert b["outside_unit_circle_count"] == 1
    assert b["spectral_radius"]["minpoly"] == [1, -7, 1]
    assert "dlPena2014mahler" in b["references"]


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_ka2_block_is_phi3(run):
    b = run({"schema": 1, "algebra": _KA2, "compute": ["coxeter_spectral"]})[
        "coxeter_spectral"]
    assert b["factorization"][0]["cyclotomic_index"] == 3
    assert b["coxeter_order"] == 3


def test_twin_parity():
    """The SAME request through hpc.spec and docs/gui/runner.py yields a
    byte-identical coxeter_spectral block (sorted-keys JSON)."""
    req = {"schema": 1, "algebra": _KRON3, "compute": ["coxeter_spectral"]}
    sblock = _server(req)["coxeter_spectral"]
    pblock = _pyodide(req)["coxeter_spectral"]
    assert json.dumps(sblock, sort_keys=True) == json.dumps(pblock, sort_keys=True)
