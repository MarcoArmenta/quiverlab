"""Plan 34: the worked-steps CONTENT must reach homework standard (Marco's blocking
feedback: "say how every object is computed as I would demand an undergraduate student
in his homework, writing everything and with justifications").

The bar, pinned here on the rendered LaTeX of a small module report: every traced
object states WHAT is computed, the DEFINITION used, the ACTUAL matrices, and WHY the
conclusion follows, with a literature justification. We assert stable, phrase-level
markers (chosen brace/underscore-free so they survive TeX escaping) for each of
rad / top / soc / resolution / Ext / tau / decompose, plus byte-determinism and a
LaTeX-compile smoke (skipped when no toolchain is on PATH -- so the fast CI matrix,
which has no TeX, skips it and it stays quick locally)."""
import shutil
import subprocess

import pytest

from quiverlab import Quiver, GF
from quiverlab.trace.modules import (
    trace_module_report, trace_radical, trace_top, trace_socle,
    trace_projective_resolution, trace_ext, trace_tor, trace_tau, trace_decompose,
)
from quiverlab.trace.render_latex import render_latex
from quiverlab.trace.provenance import references_for, resolve_references


def _square():
    """The commutative square kQ/(ab - cd) over GF(2) (Assem's worked example)."""
    Q = Quiver(vertices=[1, 2, 3, 4],
               arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=GF(2))


def _kA2():
    return Quiver(vertices=[1, 2], arrows={"a": (1, 2)}).algebra(relations=[], field=GF(2))


def _square_report_tex():
    A = _square()
    M = A.projective(1).radical()          # dim 3, non-projective, rich rad/top/soc
    M.name = "M"
    ev = trace_module_report(A, M, N=A.simple(4), top=3)
    refs = resolve_references(references_for(ev))
    return render_latex(ev, title="M = rad P_1 over the square", references=refs, algebra=A)


# --------------------------------------------------------------------------- #
# Homework depth: the definition sentences, matrices, rank arithmetic, justifications.
# --------------------------------------------------------------------------- #

# Phrase markers introduced by Plan 34, chosen to be STABLE BY DESIGN: each avoids
# TeX-special characters (_, ^, {, }) so it survives _tex_escape verbatim in prose, or
# is genuine display-math source (not escaped). Grouped by the object they certify.
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
    tex = _square_report_tex()
    assert "Nakayama" not in tex, "the mis-attributed Nakayama justification is back"
    assert "A/J is semisimple" in tex


def test_definition_sentences_present():
    tex = _square_report_tex()
    for marker in _DEFINITION_MARKERS:
        assert marker in tex, "missing definition sentence: %r" % marker


def test_justification_markers_present():
    tex = _square_report_tex()
    for marker in _JUSTIFICATION_MARKERS:
        assert marker in tex, "missing justification/citation marker: %r" % marker


def test_numbered_worked_steps():
    """Each object opens a numbered \\paragraph step (rad/top/soc/resolution/Ext/tau/
    decompose = at least seven numbered steps)."""
    tex = _square_report_tex()
    for n in range(1, 8):
        assert r"\paragraph{Step %d." % n in tex, "missing numbered Step %d" % n


def test_per_arrow_matrices_shown():
    """rad shows each arrow's action matrix rho_M(a): M -> M as a pmatrix."""
    tex = _square_report_tex()
    assert r"\rho_M(a)" in tex and r"\rho_M(b)" in tex
    assert r"\begin{pmatrix}" in tex                         # matrices are rendered
    assert r"\rho_M(a) : M \to M" in tex                     # self-map declaration


def test_radical_column_reduction_shown():
    """rad = colspace(G): the assembled matrix and its column reduction appear."""
    tex = _square_report_tex()
    assert "colspace(G)" in tex
    assert "Column-reducing" in tex                          # the reduction step
    assert r"G :" in tex and r"G = " in tex                  # G shown as a map + matrix


def test_socle_stacked_kernel_shown():
    tex = _square_report_tex()
    assert "intersection over the arrows" in tex             # soc definition
    assert r"K :" in tex and r"K = " in tex                  # the stacked system K


def test_ext_rank_arithmetic_spelled_out():
    """Ext^n = ker delta^n / im delta^{n-1} with dim = space - rank - rank = result."""
    tex = _square_report_tex()
    assert r"\ker\delta" in tex and r"\operatorname{im}\delta" in tex
    assert r"\operatorname{rank}\delta" in tex               # the rank lines
    assert r"\dim = " in tex                                 # the dimension count
    # the emitted Ext dims equal the engine's ext_dims (the binding discipline):
    A = _square()
    from quiverlab.modules.ext import ext_dims
    ev, dims = trace_ext(A, A.simple(1), A.simple(4), 3)
    assert dims == ext_dims(A, A.simple(1), A.simple(4), 3)


def test_tau_step_shows_transpose_matrix():
    tex = _square_report_tex()
    assert "transpose" in tex
    assert r"d_{1}^{*}" in tex                               # the transposed differential
    assert "IV.2" in tex                                     # ASS IV.2 (AR translate)


def test_resolution_syzygy_and_cover_narrated():
    tex = _square_report_tex()
    assert "iterated projective covers" in tex
    assert "Syzygy: Omega" in tex                            # the syzygy computation
    assert r"\varepsilon : P" in tex                         # the augmentation
    assert "Betti number" in tex                             # minimality justification


def test_decompose_certificate_labelled():
    tex = _square_report_tex()
    assert "indecomposable" in tex
    # M = rad P_1 is indecomposable here, certified by dim End = 1 (End = k*id local):
    assert "k*id is a field" in tex or "End_A(M) is local" in tex


# --------------------------------------------------------------------------- #
# Determinism: two independent builds render byte-for-byte identically.
# --------------------------------------------------------------------------- #

def test_render_is_byte_deterministic():
    a = _square_report_tex()
    b = _square_report_tex()
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
        assert render_latex(build(), title="t") == render_latex(build(), title="t")


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


# --------------------------------------------------------------------------- #
# LaTeX compile smoke: the rendered .tex compiles cleanly (skipped without a
# toolchain, so the fast CI matrix -- which has no TeX -- skips it; quick locally).
# --------------------------------------------------------------------------- #

def _latex_engine():
    for engine in ("pdflatex", "tectonic"):
        if shutil.which(engine):
            return engine
    return None


@pytest.mark.skipif(_latex_engine() is None, reason="no LaTeX toolchain on PATH")
def test_report_tex_compiles_cleanly(tmp_path):
    # Use the small kA_2 report so the smoke stays quick even in the fast bucket.
    A = _kA2()
    ev = trace_module_report(A, A.simple(1), N=A.simple(2), top=3)
    refs = resolve_references(references_for(ev))
    tex = render_latex(ev, title="S_1 over kA_2", references=refs, algebra=A)
    src = tmp_path / "report.tex"
    src.write_text(tex)
    engine = _latex_engine()
    if engine == "tectonic":
        cmd = ["tectonic", "-o", str(tmp_path), str(src)]
    else:
        cmd = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
               "-output-directory", str(tmp_path), str(src)]
    proc = subprocess.run(cmd, cwd=str(tmp_path), stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, timeout=120)
    assert proc.returncode == 0, (
        "LaTeX did not compile the homework report:\n"
        + proc.stdout.decode("utf-8", "replace")[-2000:])
    assert (tmp_path / "report.pdf").exists()
    # no overfull-box disasters (the report must be typeset cleanly):
    log = (tmp_path / "report.log").read_text(encoding="utf-8", errors="replace")
    assert "Overfull" not in log, "the report has overfull boxes (bad line breaking)"
