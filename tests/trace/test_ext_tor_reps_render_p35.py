"""Plan 35 wave 3a -- the module Ext / Tor explicit representatives are RENDERED into
the human surfaces (report Computed-results + the worked-steps module chapter), the
sibling of the HH-product per-degree layout. Rendering-only: every number is the
capture layer's, nothing recomputes.

Covered: the per-degree structure (ordered basis -> classes -> differential +
verification) with stable anchors, the inline coordinate-vector hand-check, the Tor_0
= cokernel note, tolerance (a block / event WITHOUT the reps fields), and the two-copy
gui.js wiring.
"""
import json

import pytest

from quiverlab import GF, Quiver
from quiverlab.trace import modules as tm
from quiverlab.trace.render_html import render_html
from quiverlab.trace.render_json import render_json
from quiverlab.trace.render_text import render_text
from quiverlab.trace.results_html import results_section


def _loop3():
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=["x*x*x"], field=GF(2))


def _loop_block(kind):
    """The ext/tor block a runner would ship for the x^3 loop / GF(2)."""
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    if kind == "ext":
        from quiverlab.modules.ext import ext_dims
        dims, reps = ext_dims(A, M, A.simple(1), 3, with_reps=True)
    else:
        from quiverlab.modules.tor import tor_dims
        dims, reps = tor_dims(A, M, A.simple(1, side="left"), 3, with_reps=True)
    return {"kind": kind, "top": 3, "dims": dims, "target": {"dimvec": {"1": 1}}, **reps}


# --------------------------------------------------------------------------- #
# (1) Computed-results block rendering: structure, order, anchors, hand-check.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("kind", ["ext", "tor"])
def test_results_section_renders_per_degree_sections(kind):
    html = "".join(results_section({kind: _loop_block(kind)}))
    # a per-degree heading + stable anchor per degree, in order
    for n in range(4):
        assert "cr-%s-deg-%d" % (kind, n) in html
    # ordered-basis enumeration + classes + a verification sentence
    assert "Ordered basis of" in html
    assert "Basis classes" in html
    assert "Verification" in html
    # anchors appear in increasing degree order
    order = [html.index("cr-%s-deg-%d" % (kind, n)) for n in range(4)]
    assert order == sorted(order)


@pytest.mark.oracle_literature
def test_ext_results_class_written_over_ordered_basis():
    """The Ext^1 class of the x^3 loop is the single Hom-basis element, written as its
    term-sum [P_1 → n_1,1] over the ordered basis (Marco 2026-07-31: no ``= e_1``
    coordinate tail -- coordinates are in the JSON)."""
    html = "".join(results_section({"ext": _loop_block("ext")}))
    assert "[P_1 → n_1,1]" in html
    assert "[P_1 → n_1,1] = e_1" not in html     # coordinate inline gone
    # the Ext class-label explanation is present (Marco item 6)
    assert "P_v#k" in html and "n_{v,j}" in html


@pytest.mark.oracle_literature
def test_tor0_cokernel_note_rendered():
    html = "".join(results_section({"tor": _loop_block("tor")}))
    assert "cr-tor-deg-0" in html
    assert "coker" in html                       # Tor_0 = M (x)_A N cokernel note
    assert "P_1 ⊗ n_1,1" in html                 # the tensor label


# --------------------------------------------------------------------------- #
# (2) Worked-steps module chapter: the ExtReps event renders per-degree sections.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("kind", ["ext", "tor"])
def test_worked_steps_chapter_renders_reps(kind):
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    if kind == "ext":
        events, _ = tm.trace_ext(A, M, A.simple(1), 3)
    else:
        events, _ = tm.trace_tor(A, M, A.simple(1, side="left"), 3)
    html = render_html(events, title=kind)
    assert "Explicit representatives by degree" in html
    assert "ws-%s-deg-1" % kind in html
    # the surface-distinct prefix keeps the worked-steps anchors separate from cr-
    assert "cr-%s-deg-1" % kind not in html


@pytest.mark.oracle_crossengine
def test_extreps_event_text_and_json_tolerant():
    """render_text ignores the ExtReps event (text goldens unchanged); render_json
    serializes it JSON-safely."""
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    events, _ = tm.trace_ext(A, M, A.simple(1), 2)
    txt = render_text(events, title="ext")            # must not raise
    assert "Ext" in txt
    payload = json.loads(render_json(events, title="ext"))   # must be valid JSON
    ops = [e for e in payload["events"] if e.get("type") == "ExtReps"]
    assert ops and ops[0]["op"] == "ext"


# --------------------------------------------------------------------------- #
# (3) Tolerance: a block / event WITHOUT the reps fields still renders the dims.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_results_section_tolerates_missing_reps():
    legacy = {"kind": "ext", "top": 2, "dims": [1, 1, 1],
              "target": {"dimvec": {"1": 1}}}
    html = "".join(results_section({"ext": legacy}))
    assert "ql-dims" in html                     # the dims table still renders
    assert "cr-ext-deg" not in html              # no reps sections fabricated


@pytest.mark.oracle_selfcert
def test_worked_steps_tolerates_extreps_without_fields():
    from quiverlab.trace.events import ExtReps, StepNote
    events = [StepNote("Ext run", heading=True), ExtReps(op="ext")]
    html = render_html(events, title="ext")       # ExtReps with None fields -> no crash
    assert "Explicit representatives by degree" not in html


# --------------------------------------------------------------------------- #
# (4) Both gui.js copies carry the module-reps wiring and stay byte-identical.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_gui_js_wires_module_reps():
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[2]
    docs = (root / "docs" / "gui" / "gui.js").read_text(encoding="utf-8")
    web = (root / "webapp" / "static" / "gui" / "gui.js").read_text(encoding="utf-8")
    assert docs == web, "the two gui.js copies drifted"
    for token in ("appendModuleReps", "MODULE_REPS", 'gui-" + kind + "-deg-'):
        assert token in docs, token
