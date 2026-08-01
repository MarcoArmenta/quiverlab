"""Plan 35 wave 3c -- the Yoneda exact-sequence interpretation of Ext and the classical
DICTIONARY of every space are RENDERED into the report (Computed-results HTML) and the
GUI. Rendering-only: every number is the capture layer's; the prose is the ONE shared
builder ``trace.interpretations`` both surfaces read.

Covered: the dictionary framing per theory/degree; the Ext Yoneda sequence render
(sequence line + middle module + exactness verified + connecting maps); the HH^1
derivation read-off; tolerance (a block with no interpretation); and the two-copy gui.js
wiring.
"""
import pytest

from quiverlab import GF, Quiver
from quiverlab.trace import interpretations as I
from quiverlab.trace.render_html import ext_interpretation_sections
from quiverlab.trace.results_html import _block_html


def _kA2():
    return Quiver(vertices=[1, 2], arrows={"a": (1, 2)}).algebra(relations=[], field=GF(7))


def _ext_block(A, M, N, top):
    from quiverlab.modules.ext import ext_dims
    dims, reps = ext_dims(A, M, N, top, with_reps=True, interpret=True)
    return {"kind": "ext", "top": top, "dims": [int(d) for d in dims],
            "target": {"dimvec": {"1": 0, "2": 1}, "dim": 1}, **reps}


# --------------------------------------------------------------------------- #
# (1) The classical dictionary framing renders per theory.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("kind,dims,needle", [
    ("hh_cohomology", [1, 1, 1, 1], "CENTRE"),
    ("hh_cohomology", [1, 1, 1, 1], "OUTER DERIVATIONS"),
    ("hh_cohomology", [1, 1, 1, 1], "INFINITESIMAL DEFORMATIONS"),
    ("hh_homology", [1, 0, 1], "COMMUTATOR QUOTIENT"),
    ("cyclic_homology", [2, 0, 2], "TRACE FUNCTIONALS"),
])
def test_dictionary_framing_renders(kind, dims, needle):
    html = "".join(_block_html(kind, {"kind": kind, "dims": dims, "engine": "x"}))
    assert "classical dictionary" in html
    assert needle in html


@pytest.mark.oracle_selfcert
def test_tor_dictionary_framing_renders():
    A = Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(relations=["x*x*x"], field=GF(2))
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    from quiverlab.modules.tor import tor_dims
    dims, reps = tor_dims(A, M, A.simple(1, side="left"), 3, with_reps=True)
    block = {"kind": "tor", "top": 3, "dims": dims, "target": {"dimvec": {"1": 1}}, **reps}
    html = "".join(_block_html("tor", block))
    assert "M ⊗ₐ N" in html or "tensor product" in html
    assert "flatness" in html


# --------------------------------------------------------------------------- #
# (2) The Ext Yoneda interpretation renders: sequence line, middle module, exactness.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
def test_kA2_yoneda_sequence_renders():
    A = _kA2()
    block = _ext_block(A, A.simple(1), A.simple(2), 3)
    html = "".join(_block_html("ext", block))
    # the interpretation heading + anchor
    assert "cr-ext-yoneda-deg-1" in html
    assert "as exact sequences" in html
    # the sequence line 0 -> N -> E -> M -> 0 (x-tex annotation carries the source)
    assert "0 \\to N \\to E \\to M \\to 0" in html
    # E named as the projective cover P_1 (the library's own identify_standard)
    assert "standard indecomposable" in html
    # exactness self-cert stated
    assert "Exactness verified" in html


@pytest.mark.oracle_selfcert
def test_yoneda_sections_structure_and_order():
    A = Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(relations=["x*x"], field=GF(7))
    block = _ext_block(A, A.simple(1), A.simple(1), 3)
    secs = ext_interpretation_sections(block["interpretation"], anchor_prefix="cr")
    joined = "".join(secs)
    # a heading + anchor per degree present, in increasing order
    idxs = [joined.index("cr-ext-yoneda-deg-%d" % n) for n in range(1, 4)]
    assert idxs == sorted(idxs)
    # each degree names its class and states exactness
    assert joined.count("Exactness verified") >= 3


# --------------------------------------------------------------------------- #
# (3) The HH^1 derivation read-off (hand-check on k[x]/(x^2)).
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
def test_derivation_read_off_handcheck():
    """The wave-1 anchor: HH^1 of k[x]/(x^2) has the class [x -> x]. As a derivation that
    is D(x) = x -- read straight off the captured term-sum, inventing nothing."""
    assert I.derivation_values([["1", "x", "x"]]) == ["D(x) = x"]
    # a coefficient and a two-arrow example
    assert I.derivation_values([["2", "a", "b"]]) == ["D(a) = 2 b"]
    assert I.derivation_values([["1", "a", "a"], ["1", "b", "e_1"]]) == \
        ["D(a) = a", "D(b) = e_1"]


@pytest.mark.oracle_literature
def test_dictionary_sentences_are_stable():
    # elementary-degree entries carry their distinctive keyword; the shared higher-degree
    # framing is explicitly labelled "framing".
    assert "CENTRE" in I.hh_cohomology_degree(0)
    assert "Leibniz" in I.hh_cohomology_degree(1)
    assert "framing" in I.hh_cohomology_degree(5)
    assert "Yoneda" in I.ext_degree(3)
    assert I.sentence("HH^", 0) == I.hh_cohomology_degree(0)   # alias resolves


# --------------------------------------------------------------------------- #
# (4) Tolerance: an ext block WITHOUT an interpretation renders (old cache).
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_tolerates_missing_interpretation():
    A = _kA2()
    from quiverlab.modules.ext import ext_dims
    dims, reps = ext_dims(A, A.simple(1), A.simple(2), 3, with_reps=True)   # no interpret
    block = {"kind": "ext", "top": 3, "dims": dims,
             "target": {"dimvec": {"1": 0, "2": 1}}, **reps}
    assert "interpretation" not in block
    html = "".join(_block_html("ext", block))              # must not raise
    assert "dim Ext" in html
    assert ext_interpretation_sections(None, "cr") == []
    assert ext_interpretation_sections({}, "cr") == []


# --------------------------------------------------------------------------- #
# (5) Both gui.js copies mirror the dictionary + Yoneda rendering.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_gui_js_wires_yoneda_and_dictionary():
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[2]
    a = (root / "docs/gui/gui.js").read_text(encoding="utf-8")
    b = (root / "webapp/static/gui/gui.js").read_text(encoding="utf-8")
    assert a == b, "gui.js copies must be byte-identical"
    for needle in ("appendExtInterpretation", "appendDictionaryFraming", "DICTIONARY",
                   "as exact sequences", "OUTER DERIVATIONS", "Exactness verified"):
        assert needle in a
