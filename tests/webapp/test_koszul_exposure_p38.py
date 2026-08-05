"""ext_algebra + recognizers kinds (Plan 38): served by hpc.spec, mirrored by
the Pyodide twin (docs/gui/runner.py), three-valued Koszul verdict rendered
honestly. Both runners must produce byte-identical blocks.

NO oracle-class marker: tests/webapp/ collects only with the [web] extra, so per
the Plan-32 audit these extras-gated dirs carry no oracle markers -- they are
cross-runner contract tests."""
import importlib.util
import json
import pathlib
import tempfile

import pytest

from quiverlab.hpc import spec

# the Koszul poster child k<x,y>/(x^2, y^2, x*y + y*x) over GF(7): its Yoneda
# algebra E(A) is the symmetric algebra k[x, y], graded dims 1, 2, 3, 4, ...
# (a literature-true pin: exterior <-> symmetric are Koszul dual).
_EXTERIOR = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1], "y": [1, 1]},
             "relations": ["x*x", "y*y", "x*y+y*x"], "field": {"kind": "GF", "p": 7}}
# the standard gentle A3: 1 --a--> 2 --b--> 3 with a*b = 0.
_GENTLE_A3 = {"kind": "quiver", "vertices": [1, 2, 3],
              "arrows": {"a": [1, 2], "b": [2, 3]},
              "relations": ["a*b"], "field": {"kind": "GF", "p": 5}}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p38", path)
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
def test_ext_algebra_block_shape(run):
    req = {"schema": 1, "algebra": _EXTERIOR, "compute": ["ext_algebra:0..4"]}
    b = run(req)["ext_algebra"]
    assert b["koszul"] is True
    assert b["graded_dims"][0:3] == [1, 2, 3]     # 1, 2, 3, ... = symmetric algebra
    assert "priddy" in b["references"]
    assert b["obstruction"] is None
    assert b["certified_through_degree"] == 4


@pytest.mark.parametrize("run", [_server, _pyodide])
def test_recognizers_block_flags(run):
    req = {"schema": 1, "algebra": _GENTLE_A3, "compute": ["recognizers"]}
    b = run(req)["recognizers"]
    f = b["flags"]
    assert f["is_gentle"] is True and f["is_string"] is True
    assert f["is_special_biserial"] is True
    assert f["is_hereditary"] is False
    assert b["dynkin_type"] == "A3"
    assert b["form_type"] == "finite"
    assert "assem_book" in b["references"]


def test_twin_parity():
    for req in ({"schema": 1, "algebra": _EXTERIOR, "compute": ["ext_algebra:0..4"]},
                {"schema": 1, "algebra": _GENTLE_A3, "compute": ["recognizers"]}):
        srv, pyo = _server(req), _pyodide(req)
        assert set(srv) == set(pyo)
        for kind in srv:
            assert (json.dumps(srv[kind], sort_keys=True, default=str)
                    == json.dumps(pyo[kind], sort_keys=True, default=str)), kind
