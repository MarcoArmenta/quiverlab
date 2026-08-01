"""Plan 35 wave 3b -- the cyclic-homology (HC) explicit representatives are RENDERED into
the human surfaces (report Computed-results + both gui.js copies), the sibling of the
Ext/Tor per-degree layout. Rendering-only: every number is the capture layer's.

Covered: the per-degree structure (column heading -> ordered Tot_n basis -> classes ->
total differential + verification) with stable anchors, the total-complex column-
structure heading 'Tot_n = C_n (+) C_{n-2} (+) ...', the inline coordinate-vector hand-
check, the D_0 = 0 note, tolerance (a block WITHOUT the reps fields), and the two-copy
gui.js wiring.
"""
import pathlib

import pytest

import quiverlab as ql
from quiverlab.fields import GF
from quiverlab.trace.results_html import results_section


def _hc_block(a=2, top=3, field=None):
    A = ql.truncated_polynomial(a, field=field or GF(7))
    table, reps = A.cyclic_homology(top, with_reps=True)
    return {"kind": table.kind, "top": top, "dims": list(table.dims),
            "engine": table.engine, **reps}


# --------------------------------------------------------------------------- #
# (1) Per-degree structure, anchors, order, column heading.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_results_section_renders_per_degree_sections():
    html = "".join(results_section({"cyclic_homology": _hc_block()}))
    for n in range(4):
        assert "cr-hc-deg-%d" % n in html
    assert "Ordered basis of" in html
    assert "Basis classes" in html
    assert "Verification" in html
    # anchors in increasing degree order
    order = [html.index("cr-hc-deg-%d" % n) for n in range(4)]
    assert order == sorted(order)


@pytest.mark.oracle_selfcert
def test_column_structure_heading_stated():
    """Each degree opens with the total-complex column structure Tot_n = C_n (+) ...
    before the enumeration -- deliverable-3 heading line."""
    html = "".join(results_section({"cyclic_homology": _hc_block()}))
    assert r"\mathrm{Tot}_{2}" in html                 # the LaTeX heading source
    assert r"C_{2} \oplus C_{0}" in html               # Tot_2 = C_2 (+) C_0
    assert "coordinate slices" in html                 # the slice map note


@pytest.mark.oracle_literature
def test_inline_vector_handcheck_dualnumbers():
    """HC_0 of the dual numbers: the two classes are e_1 and x (col C_0), each shown as a
    column-tagged term-sum = coordinate vector; D_0 = 0 states HC_0 = A/[A,A]."""
    html = "".join(results_section({"cyclic_homology": _hc_block()}))
    assert "col C_0: e_1" in html                       # column-tagged term-sum
    assert "= e_1" in html                              # inline coordinate vector
    assert "A/[A,A]" in html                            # the D_0 note


# --------------------------------------------------------------------------- #
# (2) Tolerance: a legacy block (no reps fields) still renders the dims table only.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_legacy_block_tolerated():
    legacy = {"kind": "cyclic_homology", "top": 3, "dims": [2, 0, 2, 0], "engine": "x"}
    html = "".join(results_section({"cyclic_homology": legacy}))
    assert "dim HC_n" in html
    assert "cr-hc-deg-" not in html                     # no reps sections
    assert "Ordered basis of" not in html


# --------------------------------------------------------------------------- #
# (3) Both gui.js copies mirror the rendering (byte-identical, appendCyclicReps).
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_gui_js_wiring_both_copies():
    root = pathlib.Path(__file__).resolve().parents[2]
    webapp = (root / "webapp" / "static" / "gui" / "gui.js").read_text(encoding="utf-8")
    docs = (root / "docs" / "gui" / "gui.js").read_text(encoding="utf-8")
    assert webapp == docs                               # byte-identical
    for needle in ("function appendCyclicReps", "appendCyclicColumnHeading",
                   "gui-cyclic-hc-deg-", r"\\mathrm{Tot}_{"):
        assert needle in webapp
    # wired into the cyclic_homology branch of renderBlock
    assert "appendCyclicReps(div, b);" in webapp
