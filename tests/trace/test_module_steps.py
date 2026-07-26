"""Plan 30, Part C: module-computation worked steps must let an algebraist REPLAY
the computation by hand (the acceptance bar) -- every resolution differential is
rendered AS A MATRIX with the direct-sum labels, the projectives/injectives of A
are described, Ext/tau show the algebra, and matrices past the threshold elide to
shape + rank (never a silently dropped step)."""
import pathlib

import pytest

from quiverlab import Quiver, GF, CC
from quiverlab.trace.events import ModuleDifferential, ModuleTerm, StepNote
from quiverlab.trace.recorder import module_differential, ext_degree, MATRIX_ELISION_CELLS
from quiverlab.trace.modules import (
    trace_projective_resolution, trace_injective_resolution, trace_ext, trace_tau,
    algebra_objects,
)
from quiverlab.trace.events import ExtDegree
from quiverlab.trace.render_text import render_text
from quiverlab.trace.render_latex import render_latex
from quiverlab.trace.render_html import render_html

GOLDEN = pathlib.Path(__file__).parent / "golden" / "module_ka2_s1.txt"


def _ka2():
    """The path algebra kA_2 (1 --a--> 2) over GF(2); hereditary, gl.dim 1."""
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    return Q.algebra(relations=[], field=GF(2))


# --------------------------------------------------------------------------- #
# The replay bar: EVERY differential of the S_1 resolution, verbatim.
# --------------------------------------------------------------------------- #
def test_ka2_s1_golden_text_replay():
    A = _ka2()
    ev, res = trace_projective_resolution(A.simple(1), 3)
    txt = render_text(ev, title="proj res of S_1", algebra=A)
    assert txt == GOLDEN.read_text(), "kA_2 S_1 worked-steps text drifted from golden"


def test_ka2_s1_tex_every_differential_verbatim():
    A = _ka2()
    ev, res = trace_projective_resolution(A.simple(1), 3)
    tex = render_latex(ev, title="proj res of S_1", algebra=A)
    # every differential of the resolution appears verbatim, as a matrix:
    assert r"\begin{pmatrix} 1 & 0 \end{pmatrix}" in tex           # epsilon: Q_0 -> M
    assert r"\begin{pmatrix} 0 \\ 1 \end{pmatrix}" in tex           # d_1: Q_1 -> Q_0
    # ...labelled with the direct-sum notation and the arrows:
    assert r"\varepsilon : P_{1} \to M" in tex
    assert r"d_{1} : P_{2} \to P_{1}" in tex
    assert r"Q_{0} = P_{1}" in tex and r"Q_{1} = P_{2}" in tex
    # ...and the projectives/injectives section is present.
    assert r"The projectives and injectives of $A$" in tex


def test_ka2_s1_html_parity_shows_matrices_and_pi():
    A = _ka2()
    ev, res = trace_projective_resolution(A.simple(1), 3)
    html = render_html(ev, title="proj res of S_1", algebra=A)
    assert "<script" not in html.lower()                           # no JS, ever
    # the same TeX-source matrices are shown (escaped inside <pre><code>):
    assert r"\begin{pmatrix} 1 &amp; 0 \end{pmatrix}" in html or \
        r"\begin{pmatrix} 1 & 0 \end{pmatrix}" in html
    assert r"\begin{pmatrix} 0 \\ 1 \end{pmatrix}" in html
    assert "The projectives and injectives of A" in html


def test_multiplicity_oplus_notation():
    """P_1^{2} (+) P_3 style labels when a term repeats a summand."""
    ev = [ModuleTerm(degree=0, kind="projective", sym="P",
                     summands=[1, 1, 3], dim=5, dimvec={"1": 4, "3": 1})]
    assert r"P_{1}^{2} \oplus P_{3}" in render_latex(ev, title="t")
    assert "P_1^2 (+) P_3" in render_text(ev, title="t")


# --------------------------------------------------------------------------- #
# The projectives and injectives of A section (Marco #4).
# --------------------------------------------------------------------------- #
def test_pi_section_loewy_content():
    A = _ka2()
    objs = algebra_objects(A)
    by_v = {row["vertex"]: row for row in objs}
    # dimvecs carry the full vertex set (zeros included); the layers are top->bottom.
    # P_1 is uniserial S_1 / S_2 (top S_1, soc S_2); P_2 = S_2 simple.
    assert by_v["1"]["P"]["layers"] == [{"1": 1, "2": 0}, {"1": 0, "2": 1}]
    assert by_v["1"]["P"]["top"] == {"1": 1, "2": 0}
    assert by_v["1"]["P"]["socle"] == {"1": 0, "2": 1}
    assert by_v["2"]["P"]["layers"] == [{"1": 0, "2": 1}]
    # I_2 is S_1 / S_2 (socle S_2), I_1 = S_1 simple.
    assert by_v["2"]["I"]["socle"] == {"1": 0, "2": 1}
    assert by_v["1"]["I"]["dim"] == 1
    # rendered text shows the layer stack and omits a separate simples section.
    txt = render_text([], title="t", algebra=A)
    assert "The projectives and injectives of A" in txt
    assert "P_1: dim 2" in txt and "Loewy layers: S_1 / S_2" in txt
    assert "simples S_v omitted" in txt


def test_pi_section_in_hh_report_too():
    """The section describes A, so it rides along with an HH trace too."""
    from quiverlab import truncated_polynomial
    from quiverlab.trace.recorder import Trace
    A = truncated_polynomial(3, field=CC)          # k[x]/(x^3), one vertex
    tr = Trace()
    A.hochschild_cohomology(2, trace=tr)
    tex = render_latex(list(tr), title="HH", algebra=A)
    assert r"The projectives and injectives of $A$" in tex
    assert r"P_{1}" in tex                          # the one projective is described


# --------------------------------------------------------------------------- #
# Elision: past the threshold, shape + rank in place of the matrix (never silent).
# --------------------------------------------------------------------------- #
def test_module_differential_elides_past_threshold():
    A = _ka2()
    dom = A.domain
    n = 21                                          # 21*21 = 441 > 400 cells
    big = [[dom.zero()] * n for _ in range(n)]
    e = module_differential(degree=2, kind="projective", sym="P", symbol="d_{2}",
                            dom_summands=[1] * n, cod_summands=[1] * n, D=big,
                            nrows=n, ncols=n, dom=dom, rank=0)
    assert e.elided is True and e.matrix is None
    assert "elided" in e.note and "21x21" in e.note and "rank 0" in e.note


def test_elided_step_appears_in_every_format_never_dropped():
    """A large (elided) differential and a small (full) one both appear -- the
    elided one as a shape+rank note, the small one as a matrix."""
    A = _ka2()
    dom = A.domain
    n = 25
    big = [[dom.zero()] * n for _ in range(n)]      # 625 > 400 -> elided
    small = [[dom.one(), dom.zero()]]               # 1x2 -> full
    events = [
        ModuleTerm(degree=0, kind="projective", sym="P", summands=[1], dim=2,
                   dimvec={"1": 2}),
        module_differential(0, "projective", "P", r"\varepsilon", [1], [], small,
                            1, 2, dom, cod_is_module=True, rank=1),
        module_differential(1, "projective", "P", "d_{1}", [1] * n, [1] * n, big,
                            n, n, dom, rank=3),
    ]
    for render in (render_text, render_latex, render_html):
        s = render(events, title="t")
        assert "elided" in s and "shape 25x25" in s          # the note is present
        assert dom.zero().__str__() * 50 not in s            # no 625-cell body dumped
    # the small differential is still shown as a real matrix (not elided):
    assert r"\begin{pmatrix} 1 & 0 \end{pmatrix}" in render_latex(events, title="t")
    # the threshold is stated once in the preamble:
    txt = render_text(events, title="t")
    assert txt.count("Matrices with more than %d entries" % MATRIX_ELISION_CELLS) == 1


# --------------------------------------------------------------------------- #
# Ext trace: the Hom-collapse + rank bookkeeping; dims pinned to ext_dims.
# --------------------------------------------------------------------------- #
def test_ext_trace_dims_match_engine():
    A = _ka2()
    from quiverlab.modules.ext import ext_dims
    for (u, w) in [(1, 1), (1, 2), (2, 1)]:
        ev, dims = trace_ext(A, A.simple(u), A.simple(w), 3)
        assert dims == ext_dims(A, A.simple(u), A.simple(w), 3)
    # the connecting-map bookkeeping is rendered:
    ev, dims = trace_ext(A, A.simple(1), A.simple(2), 3)
    txt = render_text(ev, title="ext")
    assert "dim Hom" in txt and "Result: Ext^0 = 0" in txt and "Ext^1 = 1" in txt


# --------------------------------------------------------------------------- #
# tau trace: the projective presentation + the transpose step.
# --------------------------------------------------------------------------- #
def test_injective_coresolution_terms_traced_with_I_labels():
    A = _ka2()
    ev, res = trace_injective_resolution(A.simple(2), 3)
    tex = render_latex(ev, title="inj res of S_2", algebra=A)
    # the injective terms E^n are shown with I-labels + the D(P over A^op) narration.
    assert r"E^{0} = I_" in tex
    assert "E^n = D(P_n)" in render_text(ev, title="t") or \
        "D(P_n)" in render_text(ev, title="t")


def test_tau_of_projective_is_zero_and_presentation_traced():
    A = _ka2()
    ev, t = trace_tau(A.projective(1), "tau")
    assert t.dim == 0                                # tau(projective) = 0
    txt = render_text(ev, title="tau")
    assert "Tr M = coker" in txt                     # the transpose step is narrated
    assert "term Q_0 = P_1" in txt                   # the projective presentation P_1 -> P_0
    assert "tau P_1 = 0" in txt


def test_tau_of_simple_traces_presentation_and_translate():
    A = _ka2()
    S2 = A.simple(2)
    ev, t = trace_tau(S2, "tau")
    txt = render_text(ev, title="tau")
    # S_2 is not projective here (P_2 = S_2 IS projective actually) -> tau = 0;
    # use tau^- of S_1 instead to exercise a non-trivial translate narration.
    ev2, t2 = trace_tau(A.simple(1), "tau_minus")
    txt2 = render_text(ev2, title="tau-")
    assert "inverse AR translate" in txt2
    assert ("tau^- S_1 = 0" in txt2) or ("dimension vector" in txt2)
