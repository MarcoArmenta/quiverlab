"""Pyodide-twin parity for the orbit_geometry module kind (Plan 49 / C8):
docs/gui/runner.py produces the SAME orbit_geometry block as quiverlab.hpc.spec
(both import the shared invariants.geometry.orbit_geometry_block), so the draw
page and the server/CLI agree. Mirrors tests/gui/test_tilting_check_runner_twin.py."""
import importlib.util
import json
import pathlib
import tempfile

from quiverlab.hpc import spec

_A3 = {"kind": "quiver", "vertices": [1, 2, 3], "arrows": {"a1": [1, 2], "a2": [2, 3]},
       "relations": [], "field": {"kind": "GF", "p": 32003}}


def _pyodide_block(req):
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_twin_p49b", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    r = json.loads(mod.compute_one("orbit_geometry"))
    assert r["ok"], r
    return r["block"]


def _server_block(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]["orbit_geometry"]


def test_orbit_geometry_twin_byte_identical():
    # every builtin exercises the shared builder over a large-char field (canonical
    # decomposition alive: char 32003 > dim on kA3).
    for builtin in ({"kind": "simple", "vertex": 1},
                    {"kind": "projective", "vertex": 1},
                    {"kind": "injective", "vertex": 2}):
        req = {"schema": 2, "algebra": _A3, "module": {"builtin": builtin},
               "compute": ["orbit_geometry"]}
        s = _server_block(req)
        p = _pyodide_block(req)
        assert json.dumps(s, sort_keys=True) == json.dumps(p, sort_keys=True)
