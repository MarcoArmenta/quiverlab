"""Plan 34: the worked-steps CONTENT must reach homework standard (Marco's blocking
feedback: "say how every object is computed as I would demand an undergraduate student
in his homework, writing everything and with justifications").

The bar, pinned here on the rendered HTML report of a small module report: every
traced object states WHAT is computed, the DEFINITION used, the ACTUAL matrices, and
WHY the conclusion follows, with a literature justification. We assert stable,
phrase-level markers for each of rad / top / soc / resolution / Ext / tau / decompose,
plus byte-determinism. (PDF/TeX report output has been removed; the HTML report -- the
same event stream, with the math typeset as MathML and the TeX source embedded in an
x-tex annotation -- is the homework document, exported to PDF via the browser.)"""
from quiverlab import Quiver, GF
from quiverlab.trace.modules import (
    trace_module_report, trace_radical, trace_top, trace_socle,
    trace_projective_resolution, trace_ext, trace_tor, trace_tau, trace_decompose,
)
from quiverlab.trace.render_html import render_html
from tests.trace._matrix_grid import grids, has_grid, grid_indices
from quiverlab.trace.provenance import references_for, resolve_references


def _square():
    """The commutative square kQ/(ab - cd) over GF(2) (Assem's worked example)."""
    Q = Quiver(vertices=[1, 2, 3, 4],
               arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=GF(2))


def _square_report_html():
    A = _square()
    M = A.projective(1).radical()          # dim 3, non-projective, rich rad/top/soc
    M.name = "M"
    ev = trace_module_report(A, M, N=A.simple(4), top=3)
    refs = resolve_references(references_for(ev))
    return render_html(ev, title="M = rad P_1 over the square", references=refs, algebra=A)


# --------------------------------------------------------------------------- #
# Homework depth: the definition sentences, matrices, rank arithmetic, justifications.
# The narration prose rides in <p> paragraphs; every math run keeps its TeX source
# verbatim in the x-tex <annotation>, so the LaTeX-source markers below are present.
# --------------------------------------------------------------------------- #

# Phrase markers introduced by Plan 34 -- prose (in <p>) or genuine display-math source
# (in the x-tex annotation). Grouped by the object they certify.
_DEFINITION_MARKERS = [
    "rad M = M J",                                   # radical: rad M = M J
    "top M = M / rad M",                             # top:  M / rad M
    "the largest semisimple submodule",              # socle definition
    "collapses the complex to finite-dimensional",   # Ext: the Hom collapse
    "tau M = D(Tr M)",                               # AR translate definition
    "endomorphism ring",                             # decompose via End_A(M)
]
# MINOR-6: rad M = M J is justified by A/J semisimple (J annihilates every simple),
# NOT by "Nakayama's lemma" (the old, mis-attributed sentence -- removed Plan 34).
_JUSTIFICATION_MARKERS = [
    "annihilates every simple",                      # rad = MJ justification (A/J s.s.)
    "Assem-Simson-Skowronski",                       # ASS2006 module theory
    "Green-Solberg-Zacharia",                        # GSZ2001 resolution/Ext
    "Krull-Schmidt",                                 # decomposition theorem
]


def test_nakayama_misattribution_removed():
    """MINOR-6: the radical justification no longer claims 'Nakayama's lemma' identifies
    rad M with the arrow-image sum (that is a consequence of A/J being semisimple)."""
    html = _square_report_html()
    assert "Nakayama" not in html, "the mis-attributed Nakayama justification is back"
    assert "A/J is semisimple" in html


def test_definition_sentences_present():
    html = _square_report_html()
    for marker in _DEFINITION_MARKERS:
        assert marker in html, "missing definition sentence: %r" % marker


def test_justification_markers_present():
    html = _square_report_html()
    for marker in _JUSTIFICATION_MARKERS:
        assert marker in html, "missing justification/citation marker: %r" % marker


def test_numbered_worked_steps():
    """Each object opens a numbered Step (rad/top/soc/resolution/Ext/tau/decompose =
    at least seven numbered steps)."""
    html = _square_report_html()
    for n in range(1, 8):
        assert "<b>Step %d. " % n in html, "missing numbered Step %d" % n


def test_per_arrow_matrices_shown():
    """rad shows each arrow's action matrix rho_M(a): M -> M as a pmatrix."""
    html = _square_report_html()
    assert r"\rho_M(a)" in html and r"\rho_M(b)" in html
    assert grids(html)                                       # matrices are rendered
    assert r"\rho_M(a) : M \to M" in html                    # self-map declaration


def test_radical_column_reduction_shown():
    """rad = colspace(G): the assembled matrix and its column reduction appear."""
    html = _square_report_html()
    assert "colspace(G)" in html
    assert "Column-reducing" in html                         # the reduction step
    assert r"G :" in html and r"G = " in html                # G shown as a map + matrix


def test_socle_stacked_kernel_shown():
    html = _square_report_html()
    assert "intersection over the arrows" in html            # soc definition
    assert r"K :" in html and r"K = " in html                # the stacked system K


def test_ext_rank_arithmetic_spelled_out():
    """Ext^n = ker delta^n / im delta^{n-1} with dim = space - rank - rank = result."""
    html = _square_report_html()
    assert r"\ker\delta" in html and r"\operatorname{im}\delta" in html
    assert r"\operatorname{rank}\delta" in html              # the rank lines
    assert r"\dim = " in html                                # the dimension count
    # the emitted Ext dims equal the engine's ext_dims (the binding discipline):
    A = _square()
    from quiverlab.modules.ext import ext_dims
    ev, dims = trace_ext(A, A.simple(1), A.simple(4), 3)
    assert dims == ext_dims(A, A.simple(1), A.simple(4), 3)


def test_tau_step_shows_transpose_matrix():
    html = _square_report_html()
    assert "transpose" in html
    assert r"d_{1}^{*}" in html                              # the transposed differential
    assert "IV.2" in html                                    # ASS IV.2 (AR translate)


def test_resolution_syzygy_and_cover_narrated():
    html = _square_report_html()
    assert "iterated projective covers" in html
    assert "Syzygy: Omega" in html                           # the syzygy computation
    assert r"\varepsilon : P" in html                        # the augmentation
    assert "Betti number" in html                            # minimality justification


def test_decompose_certificate_labelled():
    html = _square_report_html()
    assert "indecomposable" in html
    # M = rad P_1 is indecomposable here, certified by dim End = 1 (End = k*id local):
    assert "k*id is a field" in html or "End_A(M) is local" in html


# --------------------------------------------------------------------------- #
# Determinism: two independent builds render byte-for-byte identically.
# --------------------------------------------------------------------------- #

def test_render_is_byte_deterministic():
    a = _square_report_html()
    b = _square_report_html()
    assert a == b, "the homework report render is not byte-deterministic"


def test_tor_and_ext_and_translate_builders_deterministic():
    """Each individual builder (built twice from scratch) renders identically -- catches
    any nondeterminism in the underlying linear algebra (Hom bases, decompose order)."""
    A = _square()
    M = A.projective(1).radical()
    M.name = "M"
    builders = [
        lambda: trace_radical(M)[0],
        lambda: trace_top(M)[0],
        lambda: trace_socle(M)[0],
        lambda: trace_projective_resolution(M, 3)[0],
        lambda: trace_ext(A, M, A.simple(4), 3)[0],
        lambda: trace_tor(A, A.simple(1), A.simple(1, side="left"), 3)[0],
        lambda: trace_tau(M, "tau")[0],
        lambda: trace_decompose(M)[0],
    ]
    for build in builders:
        assert render_html(build(), title="t") == render_html(build(), title="t")


# --------------------------------------------------------------------------- #
# Citations flow through the registry plumbing (only registry keys).
# --------------------------------------------------------------------------- #

def test_module_report_cites_registry_keys():
    A = _square()
    M = A.projective(1).radical()
    M.name = "M"
    ev = trace_module_report(A, M, N=A.simple(4), top=3)
    keys = references_for(ev)
    assert "assem_book" in keys                              # ASS2006 module theory
    assert "minimal_resolution" in keys                     # GSZ2001 resolution
    assert "module_ext" in keys                             # GSZ2001 Ext
    pairs = resolve_references(keys)                         # resolves (no drift)
    assert any("Assem" in formatted for _, formatted in pairs)
    # distinct registry keys sharing GSZ2001 collapse to ONE bibliography entry:
    assert len(pairs) == len(set(pairs))


def test_tor_report_cites_tensor_product():
    A = _square()
    ev, _ = trace_tor(A, A.simple(1), A.simple(1, side="left"), 3)
    assert "tensor_product" in references_for(ev)            # Cartan-Eilenberg (Tor)
