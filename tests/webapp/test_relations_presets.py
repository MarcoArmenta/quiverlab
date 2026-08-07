"""Predetermined-relation presets (GitHub #2, Samuel Leblanc).

The two generators (rad²=0 and commutativity) live between the
``// QLGUI-RELGEN-BEGIN/END`` sentinels in gui.js. Here we extract that block,
run it under node, and check every output against an INDEPENDENT Python
enumeration (a cross-implementation oracle) -- then feed the generated relations
to the real ``quiverlab.hpc.spec.build_algebra`` and pin the resulting algebra
dimensions.
"""
import json
import pathlib
import re
import shutil
import subprocess
import tempfile

import pytest

from quiverlab.hpc.spec import build_algebra

_NODE = shutil.which("node")
_ROOT = pathlib.Path(__file__).resolve().parents[2]
_GUI_JS = _ROOT / "webapp" / "static" / "gui" / "gui.js"

pytestmark = pytest.mark.skipif(_NODE is None, reason="node not installed")


def _relgen_block():
    src = _GUI_JS.read_text(encoding="utf-8")
    m = re.search(r"// QLGUI-RELGEN-BEGIN\n(.*?)\n\s*// QLGUI-RELGEN-END", src, re.S)
    assert m, "the QLGUI-RELGEN sentinel block is missing from gui.js"
    return m.group(1)


def _run_js(fn_name, arrows):
    """Run one generator on one `arrows` object under node; return its result."""
    driver = (_relgen_block() + "\nconsole.log(JSON.stringify(" + fn_name + "("
              + json.dumps(arrows) + ")));\n")
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False) as f:
        f.write(driver)
        path = f.name
    try:
        out = subprocess.check_output([_NODE, path], encoding="utf-8")
    finally:
        pathlib.Path(path).unlink()
    return json.loads(out)


# --------------------------------------------------------------------------- #
# Independent Python oracle (a second implementation, deliberately not sharing
# code with the JS -- if both agree, the logic is right).
# --------------------------------------------------------------------------- #
def _py_rad2(arrows):
    names = sorted(arrows)
    rels = [a + "*" + b for a in names for b in names
            if arrows[a][1] == arrows[b][0]]
    return sorted(rels)


def _py_commutativity(arrows):
    names = sorted(arrows)
    verts = set()
    for a in names:
        verts.update(arrows[a])
    out_arrows = {v: [] for v in verts}
    for a in names:
        out_arrows[arrows[a][0]].append(a)
    # cycle detection
    color = {}

    def visit(v):
        color[v] = 1
        for a in out_arrows[v]:
            w = arrows[a][1]
            if color.get(w) == 1:
                return True
            if not color.get(w) and visit(w):
                return True
        color[v] = 2
        return False
    for v in verts:
        if not color.get(v) and visit(v):
            return {"error": "cycle"}
    n = len(verts)
    classes = {}

    def walk(start, cur, acc):
        if acc and len(acc) <= n:
            classes.setdefault((start, cur), []).append("*".join(acc))
        if len(acc) >= n:
            return
        for a in out_arrows[cur]:
            walk(start, arrows[a][1], acc + [a])
    for v in verts:
        walk(v, v, [])
    rels = []
    for key in sorted(classes, key=lambda k: (str(k[0]), str(k[1]))):
        paths = sorted(classes[key])
        rels.extend(paths[0] + " - " + p for p in paths[1:])
    if len(rels) > 500:
        return {"error": "cap"}
    return {"relations": rels}


# `arrows` uses string vertex keys so JSON/Python compare identically (the JS
# object keys are strings; build_algebra takes int vertices, handled per case).
_QUIVERS = {
    "kA3": {"a": [1, 2], "b": [2, 3]},
    "square": {"a": [1, 2], "b": [1, 3], "c": [2, 4], "d": [3, 4]},
    "three_parallel": {"a": [1, 2], "b": [1, 3], "c": [2, 4], "d": [3, 4],
                       "e": [1, 5], "f": [5, 4]},
    "loop": {"x": [1, 1]},
    "two_loops": {"x": [1, 1], "y": [1, 1]},
    "disconnected": {"a": [1, 2], "b": [3, 4]},
    "diamond_with_direct": {"a": [1, 2], "b": [1, 3], "c": [2, 4], "d": [3, 4],
                            "e": [1, 4]},
}


@pytest.mark.parametrize("name", sorted(_QUIVERS))
def test_rad2_matches_python_oracle(name):
    arrows = _QUIVERS[name]
    assert _run_js("qlguiRelRad2", arrows) == {"relations": _py_rad2(arrows)}


@pytest.mark.parametrize("name", sorted(_QUIVERS))
def test_commutativity_matches_python_oracle(name):
    arrows = _QUIVERS[name]
    assert _run_js("qlguiRelCommutativity", arrows) == _py_commutativity(arrows)


def test_known_small_cases():
    assert _run_js("qlguiRelRad2", {"a": [1, 2], "b": [2, 3]}) == {"relations": ["a*b"]}
    assert _run_js("qlguiRelCommutativity",
                   {"a": [1, 2], "b": [1, 3], "c": [2, 4], "d": [3, 4]}) \
        == {"relations": ["a*c - b*d"]}
    # a loop is a cycle -> commutativity refuses; rad² kills the square
    assert _run_js("qlguiRelCommutativity", {"x": [1, 1]}) == {"error": "cycle"}
    assert _run_js("qlguiRelRad2", {"x": [1, 1]}) == {"relations": ["x*x"]}
    # three parallel 1->4 paths: first-vs-rest = two relations
    assert _run_js("qlguiRelCommutativity", _QUIVERS["three_parallel"]) \
        == {"relations": ["a*c - b*d", "a*c - e*f"]}


def _field():
    return {"kind": "GF", "p": 7, "n": 1}


def _build(vertices, arrows, relations):
    return build_algebra({"kind": "quiver", "vertices": vertices,
                          "arrows": arrows, "relations": relations,
                          "field": _field()})


def test_generated_relations_build_expected_algebras():
    # kA3 with rad²=0 -> the radical-square-zero Nakayama, dim 5.
    rad2 = _run_js("qlguiRelRad2", {"a": [1, 2], "b": [2, 3]})["relations"]
    A = _build([1, 2, 3], {"a": [1, 2], "b": [2, 3]}, rad2)
    assert A.dim == 5

    # the commutative square -> dim 9, equal to the diamond poset's incidence algebra.
    comm = _run_js("qlguiRelCommutativity",
                   {"a": [1, 2], "b": [1, 3], "c": [2, 4], "d": [3, 4]})["relations"]
    B = _build([1, 2, 3, 4],
               {"a": [1, 2], "b": [1, 3], "c": [2, 4], "d": [3, 4]}, comm)
    assert B.dim == 9
    diamond = build_algebra({"kind": "family", "family": "IncidenceAlgebra",
                             "params": {"poset_or_covers": [[1, 2], [1, 3],
                                                            [2, 4], [3, 4]]},
                             "field": _field()})
    assert B.dim == diamond.dim
