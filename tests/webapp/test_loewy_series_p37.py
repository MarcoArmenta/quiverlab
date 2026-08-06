"""The rad_top_soc block carries the Loewy series; both runners agree (Plan 37 C1).

Follows the ``test_module_blocks_m0729.py`` cross-runner pattern: the server core
``quiverlab.hpc.spec`` and its Pyodide twin ``docs/gui/runner.py`` must emit the
same ``series`` (the Loewy / radical layers, top-to-bottom), or the draw page and
the desktop app disagree about the same module.

NO oracle-class marker: tests/webapp/ collects only with the [web] extra, and the
Plan-32 audit requires the audited class counts to be environment-independent --
these are contract tests for the two runners' block shapes.
"""
import json
import tempfile

import pytest

from quiverlab.hpc import spec

# kA3 / (a*b): P1 has dim 2 (dimvec {1:1, 2:1}) -- a two-layer Loewy series.
_A3 = {"kind": "quiver", "vertices": [1, 2, 3], "arrows": {"a": [1, 2], "b": [2, 3]},
       "relations": ["a*b"], "field": {"kind": "GF", "p": 5}}
_REQ = {"schema": 2, "algebra": _A3,
        "module": {"builtin": {"kind": "projective", "vertex": 1}},
        "compute": ["rad_top_soc"]}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    import importlib.util
    import pathlib
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p37", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


def _expected_layers():
    from quiverlab import GF, Quiver
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))
    return [dict(layer) for layer in A.projective(1).loewy_layers()]


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_rad_top_soc_block_has_series(run):
    b = run(_REQ)["rad_top_soc"]
    series = b["series"]
    assert isinstance(series, list) and series
    assert all(isinstance(layer, dict) for layer in series)
    # all simples are 1-dimensional here, so the total multiplicity == dim P1 == 2
    total = sum(v for layer in series for v in layer.values())
    assert total == 2
    assert series == _expected_layers()          # matches the public Module.loewy_layers()


def test_both_runners_agree_on_series():
    assert _server(_REQ)["rad_top_soc"]["series"] == _pyodide(_REQ)["rad_top_soc"]["series"]
