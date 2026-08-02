"""Marco 2026-08-01: the draw-page GUI renders HH product bidegrees as CAYLEY GRIDS.

Static-source checks that both byte-identical ``gui.js`` copies carry the Cayley
builders, plus a node harness that runs the pure JS helpers in isolation and pins
that they agree with ``quiverlab.trace.products`` byte-for-byte: balanced-rep coeff
display, the cell text in the target basis, and the derived structural notes (true
AND sabotaged-false).
"""
import json
import pathlib
import re
import shutil
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
GUI_DOCS = ROOT / "docs" / "gui" / "gui.js"
GUI_WEBAPP = ROOT / "webapp" / "static" / "gui" / "gui.js"


def test_gui_copies_carry_the_cayley_builders():
    for path in (GUI_DOCS, GUI_WEBAPP):
        src = path.read_text(encoding="utf-8")
        for name in ("primeFromBasis", "balancedCoeff", "balancedRepNote",
                     "cellTex", "cayleyStructuralNotes", "cayleyGrid",
                     "cayleyNoteLine", "combinedCayley", "combinedNote",
                     "cayleyBigGrid", "beyondWindowNote"):
            assert ("function %s" % name) in src, (path.name, name)
        assert "qlgui-cayley" in src and "qlgui-degrow" in src, path.name
        assert "Cayley" in src, path.name


def test_gui_copies_byte_identical():
    assert GUI_DOCS.read_bytes() == GUI_WEBAPP.read_bytes()


# --------------------------------------------------------------------------- #
# The pure JS helpers under node.
# --------------------------------------------------------------------------- #
def _grab_fn(src, name):
    start = src.index("function %s" % name)
    depth, started, i = 0, False, start
    while i < len(src):
        c = src[i]
        if c == "{":
            depth += 1
            started = True
        elif c == "}":
            depth -= 1
            if started and depth == 0:
                return src[start:i + 1]
        i += 1
    raise AssertionError("function %s not found / unbalanced" % name)


def _grab_var(src, name):
    m = re.search(r"var\s+%s\s*=\s*\{[^}]*\};" % re.escape(name), src)
    assert m, "var %s not found" % name
    return m.group(0)


def _node(harness):
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not available")
    proc = subprocess.run([node, "-e", harness], capture_output=True, text=True,
                          timeout=30)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


def test_js_cayley_helpers_match_python():
    src = GUI_WEBAPP.read_text(encoding="utf-8")
    pieces = [_grab_var(src, "PRODUCT_CORNER")]
    for fn in ("primeFromBasis", "balancedCoeff", "balancedRepNote", "signedJoinTex",
               "cellTex", "mirrorSign", "isIntStr", "cayleyStructuralNotes",
               "cayleyNoteLine"):
        pieces.append(_grab_fn(src, fn))
    harness = "\n".join(pieces) + r"""
// GF(2) cup(1,1): [[1,0],[0,0]] and [[0,1],[1,0]] over generators 1,2 of HH^2
var K = [[["1","0"],["0","0"]], [["0","1"],["1","0"]]];
function cell(i,j){var c=[];for(var k=0;k<2;k++)c.push(K[k][i][j]);return cellTex("cup",2,c,2);}
// a sabotaged (non-antisymmetric) single-generator table over GF(3)
var Ksab = [[["0","1"],["1","0"]]];
process.stdout.write(JSON.stringify({
  bal6: balancedCoeff("6",7), bal5: balancedCoeff("5",7),
  bal3: balancedCoeff("3",7), bal1_2: balancedCoeff("1",2),
  primeGF: primeFromBasis("bar/GF(7)"), primeQQ: primeFromBasis("cs/QQ"),
  cell00: cell(0,0), cell11: cell(1,1), cell01: cell(0,1),
  combo: cellTex("cup",2,["1","0","5"],7),
  note: cayleyNoteLine("cup",[1,1],[2,2,2],K,2),
  sabNotes: cayleyStructuralNotes("cup",[1,1],[2,2,1],Ksab,3)
}));
"""
    got = _node(harness)
    assert got["bal6"] == "-1" and got["bal5"] == "-2" and got["bal3"] == "3"
    assert got["bal1_2"] == "1"
    assert got["primeGF"] == 7 and got["primeQQ"] is None
    assert got["cell00"] == r"\alpha^{2}_{1}"        # a known nonzero cell
    assert got["cell11"] == "0"                       # a shown zero cell
    assert got["cell01"] == r"\alpha^{2}_{2}"
    assert got["combo"] == r"\alpha^{2}_{1} - 2\,\alpha^{2}_{3}"     # 5 -> balanced -2
    assert got["note"] == "The table is graded-antisymmetric."
    # the sabotaged table keeps "all squares are 0" but is NOT antisymmetric
    assert "all squares are 0" in got["sabNotes"]
    assert not any("antisymmetric" in n for n in got["sabNotes"])


def test_js_combined_cayley_big_table():
    src = GUI_WEBAPP.read_text(encoding="utf-8")
    # CAYLEY_AXIS_CAP / EM_DASH are declared as scalars alongside; inline them.
    pieces = [_grab_var(src, "PRODUCT_CORNER"),
              'var CAYLEY_AXIS_CAP = 50, EM_DASH = "\\u2014";']
    for fn in ("primeFromBasis", "balancedCoeff", "signedJoinTex", "cellTex",
               "mirrorSign", "isIntStr", "combinedOutDegree", "combinedNote",
               "combinedCayley"):
        pieces.append(_grab_fn(src, fn))
    # cup over GF(2), top 2: HH^0 dim 2, HH^1 dim 2, HH^2 dim 2 (k[x]/x^2 style shape).
    harness = "\n".join(pieces) + r"""
var K00 = [[["1","0"],["0","0"]], [["0","1"],["1","0"]]];   // (0,0) -> HH^0
var K01 = [[["1","0"],["0","0"]], [["0","1"],["1","0"]]];   // (0,1) -> HH^1 (identity-ish)
var K10 = [[["1","0"],["0","0"]], [["0","1"],["1","0"]]];   // (1,0) -> HH^1
var K11 = [[["1","0"],["0","0"]], [["0","1"],["1","0"]]];   // (1,1) -> HH^2
var K02 = [[["1","0"],["0","0"]], [["0","1"],["1","0"]]];   // (0,2) -> HH^2
var K20 = [[["1","0"],["0","0"]], [["0","1"],["1","0"]]];   // (2,0) -> HH^2
var T = [
  {degrees:[0,0], out_degree:0, dims:[2,2,2], constants:K00},
  {degrees:[0,1], out_degree:1, dims:[2,2,2], constants:K01},
  {degrees:[1,0], out_degree:1, dims:[2,2,2], constants:K10},
  {degrees:[1,1], out_degree:2, dims:[2,2,2], constants:K11},
  {degrees:[0,2], out_degree:2, dims:[2,2,2], constants:K02},
  {degrees:[2,0], out_degree:2, dims:[2,2,2], constants:K20}
];
var c = combinedCayley("cup", T, 2);
process.stdout.write(JSON.stringify({
  overCap: !!c.overCap, dl: c.dl, dr: c.dr, hasBeyond: c.hasBeyond,
  // (2,1)/(1,2)/(2,2) bidegrees are absent -> those cells must be em dashes
  cellBeyond: c.cells[4][2], degsep: c.rowDegsep
}));
"""
    got = _node(harness)
    assert got["overCap"] is False
    assert got["dl"] == 6 and got["dr"] == 6          # degree-major, all classes
    assert got["hasBeyond"] is True                    # (2,2) etc. beyond the window
    assert got["cellBeyond"] == "—"                    # a beyond-window em dash
    assert got["degsep"] == [False, False, True, False, True, False]  # degree boundaries
