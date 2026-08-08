"""Plan 51 Task G.1: the CS-native bracket is clickable end-to-end and BYTE-IDENTICAL
across the two runners (quiverlab.hpc.spec server + the docs/gui Pyodide twin).

Unmarked (extras-gated dir -- no oracle markers in tests/webapp, per Plan 32).
Every GUI/webapp algebra is presented (quiver/family), so the CS bracket ALWAYS
serves here; the presentation-less refusal is a library-level guard, covered by
tests/hochschild/test_products_api.py::test_bracket_refuses_presentationless_off_gfp.
"""
import importlib.util
import json
import pathlib
import tempfile

import pytest

from quiverlab.hpc import spec

_RUNNER = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"


def _server_block(body, kind):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(body, d)["results"][kind]


def _gui_block(body, item):
    s = importlib.util.spec_from_file_location("_gui_runner_p51", _RUNNER)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    r = json.loads(mod.compute_one(item))
    assert r["ok"], r
    return r["block"]


_QQ_BODY = {
    "schema": 2,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x"], "field": {"kind": "QQ"}},
    "artifacts": {"pdf": False, "tikz": False},
    "compute": ["bracket:0..2"],
}

_GFP_BODY = {
    "schema": 2,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 7, "n": 1}},
    "artifacts": {"pdf": False, "tikz": False},
    "compute": ["bracket:0..2"],
}


def test_qq_bracket_served_with_cs_provenance():
    b = _server_block(_QQ_BODY, "bracket")
    assert b["kind"] == "bracket"
    assert b["engine"] == "Chouhy-Solotar native diagonal (homotopy lifting)"
    assert b["basis"] == "cs/QQ"
    assert "window" not in b                    # native -> no bar window field
    assert "oke_koszul" in b["references"]


def test_qq_bracket_byte_identical_across_runners():
    server = _server_block(_QQ_BODY, "bracket")
    gui = _gui_block(_QQ_BODY, "bracket:0..2")
    assert (json.dumps(server, sort_keys=True, default=str)
            == json.dumps(gui, sort_keys=True, default=str))


def test_gfp_in_window_bracket_carries_window_and_matches_runners():
    server = _server_block(_GFP_BODY, "bracket")
    assert server.get("window") == 2                    # GF(p) bar route: bounded
    assert server["basis"].startswith("bar/GF(")
    gui = _gui_block(_GFP_BODY, "bracket:0..2")
    assert (json.dumps(server, sort_keys=True, default=str)
            == json.dumps(gui, sort_keys=True, default=str))
