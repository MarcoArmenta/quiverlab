"""Pyodide-twin parity for the tilting_check module kind (Plan 44 / C7): docs/gui/runner.py
produces the SAME tilting_check block (math subkeys) as quiverlab.hpc.spec, so the draw
page and the server/CLI agree. Mirrors tests/gui/test_runner_modules.py's twin style."""
import importlib.util
import json
import pathlib
import tempfile

from quiverlab.hpc import spec

_A2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
       "relations": [], "field": {"kind": "GF", "p": 5}}
_MATH_KEYS = ("kind", "is_tilting", "n", "pd", "self_ext_vanishes",
              "num_summands", "num_vertices", "note")


def _pyodide_block(req, item):
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_twin_p44", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    r = json.loads(mod.compute_one(item))
    assert r["ok"], r
    return r["block"]


def _server_block(req, item):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"][item.split(":")[0]]


def test_tilting_check_twin_math_subkeys():
    # projective P1 (not tilting) and the regular-ish injective I2 exercise both verdicts.
    for builtin in ({"kind": "projective", "vertex": 1},
                    {"kind": "injective", "vertex": 2}):
        req = {"schema": 2, "algebra": _A2, "module": {"builtin": builtin},
               "compute": ["tilting_check"]}
        s = _server_block(req, "tilting_check")
        p = _pyodide_block(req, "tilting_check")
        sub_s = {k: s[k] for k in _MATH_KEYS}
        sub_p = {k: p[k] for k in _MATH_KEYS}
        assert json.dumps(sub_s, sort_keys=True) == json.dumps(sub_p, sort_keys=True)
