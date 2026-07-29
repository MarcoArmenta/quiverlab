"""Plan 34 (post-critique) -- the render_html LaTeX->MathML converter FIDELITY GUARD
(MAJOR-1) and the op-aware Ext/Tor HTML rendering (MAJOR-2).

MAJOR-1: the converter must never leak a backslash-stripped command name into a
typeset atom (the old `\\frac`->"frac 1 2", `\\begin{array}{cc}`->"cc", `\\substack`,
`\\sqrt`/`\\mathbb`/`\\xrightarrow` bugs). Every real module+HH sample is rendered
through render_html and its PRESENTATION MathML (the <annotation> raw source stripped)
is scanned: no <mi>/<mo>/<mn> may carry a backslash, and no typeset atom may be a
leaked command name. The advertised <mtext> fallback is REAL -- an unknown command
routes the whole run to <mtext> of the escaped source (backslashes PRESERVED, so it is
readable source, not garbled math). The leak detector is self-tested so a real
regression is provably caught.

MAJOR-2: an ExtDegree with op="Tor" renders with tensor-dimension wording, d_n, and a
Tor_n SUBSCRIPT (not "dim Hom", delta^n, Ext^n superscript), with the rank arithmetic
spelled out -- so the no-toolchain HTML surface matches the PDF.
"""
import re

import pytest

from quiverlab import Quiver, GF, CC, truncated_polynomial
from quiverlab.trace.recorder import Trace, module_differential, ext_degree
from quiverlab.trace.events import ModuleTerm, StepNote
from tests.trace._result_table import result_dims
from quiverlab.trace.render_html import (
    render_html, _math, _tex_to_mathml_body, _UnknownTeX,
)
from quiverlab.trace.modules import (
    trace_projective_resolution, trace_injective_resolution, trace_ext, trace_tor,
    trace_module_report,
)


# --------------------------------------------------------------------------- #
# Leak detection (self-tested below).
# --------------------------------------------------------------------------- #

# Command NAMES that can only appear in a typeset atom via the backslash-stripping
# leak bug -- they are never legitimate <mi>/<mo>/<mn>/<mtext> content in our grammar
# (unlike dim/Hom/ker/Tor/Ext/im/rank, which ARE legitimate operator identifiers).
_LEAK_NAMES = (
    "frac", "substack", "mathbb", "mathcal", "mathfrak", "mathscr", "sqrt",
    "xrightarrow", "xleftarrow", "resizebox", "array", "cc", "begin", "end",
)
_ATOM = re.compile(r"<(mi|mo|mn)>(.*?)</\1>", re.DOTALL)
_TEXT = re.compile(r"<mtext>(.*?)</mtext>", re.DOTALL)
_ANNOT = re.compile(r"<annotation\b.*?</annotation>", re.DOTALL)


def _presentation(html):
    """The typeset MathML with the x-tex <annotation> (which legitimately holds the raw
    source, backslashes and all) removed -- so the scan sees only typeset atoms."""
    return _ANNOT.sub("", html)


def _leaks(mathml):
    """Return the list of leak offences in a chunk of PRESENTATION MathML:
      * any <mi>/<mo>/<mn> whose content carries a backslash (a raw command that
        escaped into a typeset atom -- e.g. the old `\\substack`->`<mi>\\</mi>`);
      * any typeset atom OR <mtext> equal to a bare leak command name.
    An <mtext> that holds escaped SOURCE (with its backslashes, the honest fallback)
    is NOT a leak -- only a backslash-free bare command name is."""
    offences = []
    for tag, content in _ATOM.findall(mathml):
        if "\\" in content:
            offences.append("backslash in <%s>%s</%s>" % (tag, content, tag))
        if content.strip() in _LEAK_NAMES:
            offences.append("leaked name <%s>%s</%s>" % (tag, content, tag))
    for content in _TEXT.findall(mathml):
        # The honest <mtext> fallback holds escaped LaTeX SOURCE, which always carries a
        # backslash (\\frac{1}{2}, \\begin{array}{cc}) -- exempt it. A backslash-FREE
        # <mtext> that is a bare command name is a genuine typeset-text leak.
        if "\\" in content:
            continue
        for name in _LEAK_NAMES:
            if re.search(r"\b%s\b" % name, content):
                offences.append("leaked name in <mtext>: %r (%s)" % (content, name))
    return offences


# --------------------------------------------------------------------------- #
# Sample algebras / event streams (real builders + hand-built edge events).
# --------------------------------------------------------------------------- #

def _kA2():
    return Quiver(vertices=[1, 2], arrows={"a": (1, 2)}).algebra(relations=[], field=GF(2))


def _square():
    Q = Quiver(vertices=[1, 2, 3, 4],
               arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=GF(2))


def _hh_events():
    A = truncated_polynomial(3, field=CC)          # wide (>10-col) HH^2 differential
    tr = Trace()
    A.hochschild_cohomology(2, trace=tr)
    return A, list(tr)


def _handbuilt_events():
    """Edge events that exercise the full emitted grammar directly (independent of the
    module builders): a numbered heading, both Ext and Tor degrees, a self-map
    differential (dom_is_module), and the exotic symbols rho_M / varphi / d_1^*."""
    dom = GF(2)
    D = [[dom.one(), dom.zero()], [dom.zero(), dom.one()]]
    return [
        StepNote("Opening step", "with a justification.", heading=True),
        ModuleTerm(degree=0, kind="projective", sym="P", summands=[1, 1, 3], dim=5,
                   dimvec={"1": 4, "3": 1}),
        module_differential(0, "module", "Me", r"\rho_M(a)", [1], [1], D, 2, 2, dom,
                            cod_is_module=True, dom_is_module=True, rank=2),
        module_differential(0, "module", "Me", r"\varphi", [1], [1], D, 2, 2, dom,
                            cod_is_module=True, dom_is_module=True, rank=2),
        module_differential(1, "projective", "P", r"d_{1}^{*}", [1], [2], D, 2, 2, dom,
                            cod_is_module=False, rank=2),
        ext_degree(degree=1, op="Ext", space_dim=3, rank_here=1, rank_prev=1,
                   result_dim=1, D=D, nrows=2, ncols=2, dom=dom),
        ext_degree(degree=1, op="Tor", space_dim=3, rank_here=1, rank_prev=1,
                   result_dim=1, D=D, nrows=2, ncols=2, dom=dom),
    ]


def _all_samples():
    A, hh = _hh_events()
    sq = _square()
    M = A  # placeholder to avoid confusion; not used
    return [
        ("hh", hh, A),
        ("proj_res", trace_projective_resolution(_kA2().simple(1), 3)[0], _kA2()),
        ("inj_res", trace_injective_resolution(_kA2().simple(2), 3)[0], _kA2()),
        ("ext", trace_ext(sq, sq.simple(1), sq.simple(4), 3)[0], sq),
        ("tor", trace_tor(sq, sq.simple(1), sq.simple(1, side="left"), 3)[0], sq),
        ("module_report", trace_module_report(
            sq, sq.projective(1).radical(), N=sq.simple(4), top=3), sq),
        ("handbuilt", _handbuilt_events(), None),
    ]


# --------------------------------------------------------------------------- #
# The leak detector itself must catch leaks and pass honest fallbacks.
# --------------------------------------------------------------------------- #

def test_leak_detector_flags_a_leak_but_not_the_honest_fallback():
    # the OLD-bug forms are flagged:
    assert _leaks("<mi>frac</mi>"), "detector missed a leaked command name"
    assert _leaks(r"<mi>\</mi>"), "detector missed a backslash in a typeset atom"
    assert _leaks("<mo>cc</mo>"), "detector missed a leaked array colspec"
    # the honest <mtext> fallback (source WITH its backslash) is NOT a leak:
    assert not _leaks(r"<mtext>\frac{1}{2}</mtext>")
    assert not _leaks(r"<mtext>\begin{array}{cc}</mtext>")
    # legitimate operator identifiers are NOT leaks:
    assert not _leaks("<mi>dim</mi><mi>Hom</mi><mi>Tor</mi><mi>ker</mi>")


# --------------------------------------------------------------------------- #
# MAJOR-1: real samples convert with zero leaks.
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", [s[0] for s in _all_samples()])
def test_every_sample_converts_without_a_leak(name):
    events, algebra = {s[0]: (s[1], s[2]) for s in _all_samples()}[name]
    html = render_html(events, title=name, algebra=algebra)
    offences = _leaks(_presentation(html))
    assert not offences, "render_html(%s) leaked: %s" % (name, offences[:5])


def test_chained_sub_and_superscript_are_both_preserved():
    body = _tex_to_mathml_body(r"a_{i}^{j}")
    assert "<msubsup>" in body, "chained scripts collapsed"
    assert "<mi>i</mi>" in body and "<mi>j</mi>" in body, "a script was dropped"
    # a real emitted form: P_{1}^{2} keeps both the vertex sub and the multiplicity sup.
    p = _tex_to_mathml_body(r"P_{1}^{2}")
    assert "<msubsup>" in p and "<mn>1</mn>" in p and "<mn>2</mn>" in p


def test_unknown_command_routes_the_whole_run_to_the_mtext_fallback():
    # unknown command -> _UnknownTeX -> _math renders the WHOLE run as escaped source
    with pytest.raises(_UnknownTeX):
        _tex_to_mathml_body(r"\wibble x")
    html = _math(r"\frac{1}{2} + \wibble")
    pres = _presentation(html)
    assert "<mtext>" in pres, "the fallback must be a real <mtext>, not a leak"
    # NO backslash-stripped command name typeset as an identifier:
    assert not _leaks(pres), "the fallback itself leaked: %s" % _leaks(pres)
    # the source is still readable (kept verbatim in the annotation):
    assert r"\frac{1}{2}" in html


def test_resizebox_is_unwrapped_not_leaked():
    body = _tex_to_mathml_body(r"\resizebox{5pt}{3pt}{\delta^{2}}")
    assert "<msup>" in body and "&#948;" in body      # the content typesets
    assert not _leaks(body)                            # no "resizebox" leak


# --------------------------------------------------------------------------- #
# MAJOR-2: op-aware Ext vs Tor rendering (wording, subscript vs superscript, arith).
# --------------------------------------------------------------------------- #

def _ext_html():
    dom = GF(2)
    D = [[dom.one(), dom.zero()], [dom.zero(), dom.one()]]
    ev = [ext_degree(degree=1, op="Ext", space_dim=3, rank_here=1, rank_prev=1,
                     result_dim=1, D=D, nrows=2, ncols=2, dom=dom)]
    return render_html(ev, title="ext")


def _tor_html():
    dom = GF(2)
    D = [[dom.one(), dom.zero()], [dom.zero(), dom.one()]]
    ev = [ext_degree(degree=1, op="Tor", space_dim=3, rank_here=1, rank_prev=1,
                     result_dim=1, D=D, nrows=2, ncols=2, dom=dom)]
    return render_html(ev, title="tor")


def test_ext_html_uses_hom_delta_and_ext_superscript():
    html = _ext_html()
    assert r"\operatorname{Ext}^{1}(M,N)" in html
    assert r"\dim\operatorname{Hom}(Q_{1},N)" in html
    assert r"\ker\delta^{1}/\operatorname{im}\delta^{0}" in html
    assert r"\dim = 3 - 1 - 1 = 1" in html            # the rank arithmetic, spelled out
    # the Result footer is a degree table (Marco 2026-07-29) labelled with the Ext
    # SUPERSCRIPT (degrees re-indexed from 0, as in the shared ext_result_dims):
    assert result_dims(html, "dim Ext^n") == [1]
    # NOT Tor wording:
    assert r"\operatorname{Tor}" not in html and r"\otimes_A" not in html


def test_tor_html_uses_tensor_d_n_and_tor_subscript():
    html = _tor_html()
    assert r"\operatorname{Tor}_{1}(M,N)" in html
    assert r"\dim(Q_{1}\otimes_A N)" in html          # tensor-dimension wording
    assert r"d_{2}" in html                            # d_{n+1} is the shown map
    assert r"\ker d_{1}/\operatorname{im}d_{2}" in html
    assert r"\dim = 3 - 1 - 1 = 1" in html
    # the Result footer is a degree table labelled with the Tor SUBSCRIPT, not a
    # superscript (degrees re-indexed from 0, as in the shared ext_result_dims):
    assert result_dims(html, "dim Tor_n") == [1]
    # NOT Ext/Hom wording:
    assert r"\operatorname{Ext}" not in html
    assert "dim Hom" not in html and r"\delta" not in html


def test_numbered_step_headings_render_in_html():
    ev = [StepNote("Radical of M", "By definition rad M = M J.", heading=True),
          StepNote("Top of M", "top M = M / rad M.", heading=True)]
    html = render_html(ev, title="t")
    assert "<b>Step 1. Radical of M</b>" in html
    assert "<b>Step 2. Top of M</b>" in html


def test_self_map_differential_shows_M_on_both_sides():
    dom = GF(2)
    D = [[dom.one()]]
    ev = [module_differential(0, "module", "Me", r"\rho_M(a)", [1], [1], D, 1, 1, dom,
                              cod_is_module=True, dom_is_module=True, rank=1)]
    html = render_html(ev, title="t")
    assert r"\rho_M(a) : M \to M" in html, "dom_is_module not honored in HTML"


def test_zero_dimensional_matrix_renders_zero_in_html():
    """A 0xk / kx0 differential (every zero Ext/Tor connecting map) renders as the
    symbol ``0`` in the shared ``_pmatrix`` -- never an empty ``\\begin{pmatrix}
    \\end{pmatrix}`` that would typeset as a stray ``()`` (Plan 34 BLOCKING-2)."""
    from quiverlab.trace.render_html import _pmatrix

    class _E:                                   # a 0xk / kx0 differential event
        matrix, nrows, ncols, elided, note = [], 0, 3, False, ""
    assert _pmatrix(_E()) == "0"
    _E.matrix, _E.nrows, _E.ncols = [[], [], []], 3, 0
    assert _pmatrix(_E()) == "0"
    # ...and in a real Ext report the empty environment never appears:
    A = _square()
    ev, _ = trace_ext(A, A.simple(1), A.projective(1), 3)
    html = render_html(ev, title="t")
    assert r"\begin{pmatrix} \end{pmatrix}" not in html
    assert r"\begin{pmatrix}\end{pmatrix}" not in html


def test_render_html_module_report_is_byte_deterministic():
    sq = _square()
    M = sq.projective(1).radical()
    M.name = "M"
    ev = trace_module_report(sq, M, N=sq.simple(4), top=3)
    assert render_html(ev, title="t", algebra=sq) == render_html(ev, title="t", algebra=sq)
