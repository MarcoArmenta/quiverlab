"""Behavioral tests for the wave-2 GUI wiring fix round (derived_compare's
algebra B, the potential box, and the algebra_b gate).

Two layers:

* the pure ``dynkinQuiver`` synthesiser (between the ``// QLGUI-DYNKIN`` sentinels
  in gui.js) is node-executed and its output compared EXACTLY to the quiver the
  real ``quiverlab.PathAlgebra`` builds -- a cross-implementation oracle, so the
  client-synthesised algebra B is byte-for-byte the library's Dynkin quiver and
  therefore computes on every tier (including the browser Pyodide runner, which
  refuses ``kind:"family"``); and

* control-flow invariants on the real gui.js source: ``req.algebra_b`` is written
  ONLY under the ``derived_compare`` checkbox guard (so it vanishes when the kind
  is off), and every programmatic relations write clears the potential box +
  refreshes the mutual-exclusion hint (so a stale potential can't ride along).
"""
import json
import pathlib
import re
import shutil
import subprocess
import tempfile

import pytest

import quiverlab as ql

_NODE = shutil.which("node")
_ROOT = pathlib.Path(__file__).resolve().parents[2]
_GUI_JS = _ROOT / "webapp" / "static" / "gui" / "gui.js"

pytestmark = pytest.mark.skipif(_NODE is None, reason="node not installed")


# --------------------------------------------------------------------------- #
# (a) dynkinQuiver === the library's PathAlgebra quiver
# --------------------------------------------------------------------------- #
def _dynkin_block():
    src = _GUI_JS.read_text(encoding="utf-8")
    m = re.search(r"// QLGUI-DYNKIN-BEGIN\n(.*?)\n\s*// QLGUI-DYNKIN-END", src, re.S)
    assert m, "the QLGUI-DYNKIN sentinel block is missing from gui.js"
    return m.group(1)


def _run_dynkin(type_str):
    driver = (_dynkin_block() + "\nconsole.log(JSON.stringify(dynkinQuiver("
              + json.dumps(type_str) + ")));\n")
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False) as f:
        f.write(driver)
        path = f.name
    try:
        out = subprocess.check_output([_NODE, path], encoding="utf-8")
    finally:
        pathlib.Path(path).unlink()
    return json.loads(out)


def _library_quiver(type_str):
    """The library's PathAlgebra quiver as {vertices:[ints], arrows:{name:[s,t]}}
    -- the exact shape dynkinQuiver must reproduce."""
    Q = ql.PathAlgebra(type_str).quiver
    return {"vertices": list(Q.vertices),
            "arrows": {a: [Q.source(a), Q.target(a)] for a in Q.arrows}}


# A1 is excluded: quiverlab.PathAlgebra("A1") itself raises (no oracle to compare
# against), yet a single-vertex quiver is a perfectly good algebra -- so client
# synthesis is strictly MORE capable there (pinned in test_dynkin_A1 below).
@pytest.mark.parametrize("type_str", ["A2", "A3", "A5", "D4", "D5", "E6", "E7", "E8"])
def test_dynkin_quiver_matches_the_library(type_str):
    got = _run_dynkin(type_str)
    assert got is not None, f"dynkinQuiver refused a supported type {type_str}"
    assert got == _library_quiver(type_str), (
        f"{type_str}: client synthesis differs from quiverlab.PathAlgebra")


def test_dynkin_A3_exact_block():
    # The exact JSON the "A3" B-input produces (the critic's pin (a)).
    assert _run_dynkin("A3") == {
        "vertices": [1, 2, 3],
        "arrows": {"e12": [1, 2], "e23": [2, 3]}}


def test_dynkin_A1_is_the_single_vertex_quiver():
    # The library's PathAlgebra("A1") string raises, but the single-vertex quiver
    # block builds fine (dim 1) -- client synthesis works where the string fails.
    from quiverlab.hpc.spec import build_algebra
    assert _run_dynkin("A1") == {"vertices": [1], "arrows": {}}
    A = build_algebra({"kind": "quiver", "vertices": [1], "arrows": {},
                       "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}})
    assert A.dim == 1


@pytest.mark.parametrize("bad", ["A0", "D3", "E5", "E9", "B2", "F4", "", "junk", "A"])
def test_dynkin_quiver_refuses_unsupported(bad):
    assert _run_dynkin(bad) is None


def test_synthesised_block_builds_a_real_algebra():
    # A dynkinQuiver block, wrapped as algebra B and fed to the real builder,
    # must construct -- the whole point of client synthesis (works on every tier).
    from quiverlab.hpc.spec import build_algebra
    q = _run_dynkin("D4")
    A = build_algebra({"kind": "quiver", "vertices": q["vertices"],
                       "arrows": q["arrows"], "relations": [],
                       "field": {"kind": "GF", "p": 5, "n": 1}})
    assert A.dim == ql.PathAlgebra("D4", field=ql.GF(5)).dim


# --------------------------------------------------------------------------- #
# (b)/(c) control-flow invariants on the real gui.js source
# --------------------------------------------------------------------------- #
def _gui_src():
    return _GUI_JS.read_text(encoding="utf-8")


def _function_body(src, header):
    """Brace-match a `function name(...) {...}` (or a labelled handler) body."""
    i = src.index(header)
    b = src.index("{", i)
    depth, j = 0, b
    while j < len(src):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[b:j + 1]
        j += 1
    raise AssertionError(f"unbalanced braces after {header!r}")


def test_algebra_b_only_written_under_the_derived_compare_guard():
    src = _gui_src()
    # every occurrence of the assignment must sit inside the checkbox guard block
    guard = _function_body(src, "if (el.derived_compare.checked) {\n      var algb")
    assert "req.algebra_b = algb" in guard
    assert src.count("req.algebra_b") == guard.count("req.algebra_b"), (
        "req.algebra_b is referenced outside the derived_compare guard -- it "
        "could leak into a request whose derived_compare box is unchecked")


def test_potential_attached_only_when_non_empty():
    src = _gui_src()
    # req.algebra.potential is written once, guarded by a non-empty `pot`
    assert re.search(r"if \(pot\)\s*req\.algebra\.potential = pot;", src), (
        "the potential must be attached only when typed (cache-key discipline)")


def test_programmatic_relation_writes_clear_the_potential():
    src = _gui_src()
    load = _function_body(src, "function loadQuiver(p)")
    assert 'el.potential.value = ""' in load, "loadQuiver leaves a stale potential"
    assert "syncRelPotConflict" in load, "loadQuiver doesn't refresh the hint"
    clear = _function_body(src, "el.clear.addEventListener")
    assert 'el.potential.value = ""' in clear, "Clear leaves a stale potential"
    preset = _function_body(src, "function applyRelPreset(kind)")
    assert "syncRelPotConflict" in preset, (
        "a generated relation set doesn't refresh the potential-vs-relations hint")
