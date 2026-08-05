"""The almost_split module kind (Plan 41): schema-2 gated, served by hpc.spec,
mirrored byte-for-byte by the Pyodide twin, honest refusal for projective input.

Both runners are checked (server ``quiverlab.hpc.spec`` and the Pyodide twin
``docs/gui/runner.py``): the draw page and the desktop app must agree.
"""
import json
import pathlib
import tempfile

from quiverlab.hpc import spec

# NO oracle-class marker: tests/webapp/ collects only with the [web] extra and the
# Plan-32 audit requires environment-independent class counts. Contract test.

# kA3 (1 -> 2 -> 3, linear) over GF(5): all indecomposable dims <= 3 < 5, so the
# Krull-Schmidt/rad-End engines are reliable (char > dim).
_KA3 = {"kind": "quiver", "vertices": [1, 2, 3],
        "arrows": {"a": [1, 2], "b": [2, 3]}, "relations": [],
        "field": {"kind": "GF", "p": 5}}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    import importlib.util
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p41", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


def test_almost_split_block_shape():
    req = {"schema": 2, "algebra": _KA3,
           "module": {"builtin": {"kind": "simple", "vertex": 2}},
           "compute": ["almost_split"]}
    block = _server(req)["almost_split"]
    assert block["exists"] is True
    assert set(block["tau"]) >= {"dims", "maps"}          # tau M as a full representation
    assert isinstance(block["middle"]["summands"], list) and block["middle"]["summands"]
    assert block["indecomposable"] is True
    assert block["latex"] == r"0 \to \tau M \to E \to M \to 0"
    # references carry the registry keys; citations carry the resolved bibtex-keyed pairs
    assert "assem_book" in block["references"] and "ars_book" in block["references"]
    assert "ARS1995" in [k for k, _ in block["citations"]]


def test_almost_split_projective_refused():
    req = {"schema": 2, "algebra": _KA3,
           "module": {"builtin": {"kind": "projective", "vertex": 1}},
           "compute": ["almost_split"]}
    block = _server(req)["almost_split"]                  # NO exception escapes
    assert block["exists"] is False
    assert "reason" in block and "projective" in block["reason"]


def test_twin_parity():
    req = {"schema": 2, "algebra": _KA3,
           "module": {"builtin": {"kind": "simple", "vertex": 2}},
           "compute": ["almost_split"]}
    s = _server(req)["almost_split"]
    p = _pyodide(req)["almost_split"]
    assert json.dumps(s, sort_keys=True) == json.dumps(p, sort_keys=True)


def test_twin_parity_on_refusal():
    req = {"schema": 2, "algebra": _KA3,
           "module": {"builtin": {"kind": "projective", "vertex": 1}},
           "compute": ["almost_split"]}
    s = _server(req)["almost_split"]
    p = _pyodide(req)["almost_split"]
    assert json.dumps(s, sort_keys=True) == json.dumps(p, sort_keys=True)
