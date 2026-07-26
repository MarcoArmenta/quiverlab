"""Plan 34 (post-critique) -- the worked-steps PDF is a PAGE-BOUNDED homework document
that actually COMPILES, and the flagship arithmetic is exercised on non-degenerate data:

  * BLOCKING-1: an 11..25-column matrix (homework-scale, e.g. a dim-12 module's 12x12
    arrow action) is too wide for amsmath's default pmatrix but too small to have elided
    under the old 400-cell rule, so pdflatex aborted ("Extra alignment tab"). The renderer
    now emits \\setcounter{MaxMatrixCols} + a graphicx \\qlmat shrink-to-width wrap, and a
    matrix past 25 rows/cols is STATED-elided (shown in full in the HTML/JSON report).
  * BLOCKING-2: a zero-DIMENSIONAL matrix (every zero Ext/Tor differential) rendered as an
    empty \\begin{pmatrix}\\end{pmatrix} -> a stray "(" in the PDF; it now renders "0".
  * MAJOR-3: the acceptance samples used simple second arguments, so the Ext rank
    arithmetic only ever showed 0-0-0=0; here Ext(S_1, P_1) over the commutative square has
    a genuinely NONZERO differential delta^0 = (1;1) (rank 1), so the arithmetic reads
    1 - 1 - 0 and 2 - 1 - 1.

The compile smokes are skipped without a LaTeX toolchain (so the fast CI matrix skips
them); the byte-determinism + rendered-token checks always run."""
import shutil
import subprocess

import pytest

from quiverlab import Quiver, GF, truncated_polynomial
from quiverlab.trace.events import ModuleDifferential
from quiverlab.trace.recorder import module_differential
from quiverlab.trace.modules import trace_module_report, trace_ext
from quiverlab.trace.render_latex import (
    render_latex, _pmatrix, matrix_preamble_lines, _events_max_cols,
)
from quiverlab.trace.provenance import references_for, resolve_references


def _latex_engine():
    for engine in ("pdflatex", "tectonic"):
        if shutil.which(engine):
            return engine
    return None


def _compile(tex, tmp_path):
    """Compile ``tex`` and return (returncode, overfull_hbox_count)."""
    src = tmp_path / "r.tex"
    src.write_text(tex)
    engine = _latex_engine()
    if engine == "tectonic":
        cmd = ["tectonic", "-o", str(tmp_path), str(src)]
    else:
        cmd = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
               "-output-directory", str(tmp_path), str(src)]
    proc = subprocess.run(cmd, cwd=str(tmp_path), stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, timeout=180)
    log = (tmp_path / "r.log")
    overfull = (log.read_text(encoding="utf-8", errors="replace").count("Overfull \\hbox")
                if log.exists() else 0)
    return proc, overfull


def _square():
    Q = Quiver(vertices=[1, 2, 3, 4],
               arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=GF(2))


# --------------------------------------------------------------------------- #
# BLOCKING-2: zero-dimensional matrices render "0", never an empty pmatrix.
# --------------------------------------------------------------------------- #
def test_zero_dimensional_matrix_renders_zero_not_empty_pmatrix():
    class _E:                                   # a 0xk / kx0 differential event
        matrix, nrows, ncols, elided, note = [], 0, 3, False, ""
    assert _pmatrix(_E()) == "0"
    _E.matrix, _E.nrows, _E.ncols = [[], [], []], 3, 0
    assert _pmatrix(_E()) == "0"
    # ...and in a real Ext report the empty environment never appears:
    A = _square()
    ev, _ = trace_ext(A, A.simple(1), A.projective(1), 3)
    tex = render_latex(ev, title="t")
    assert r"\begin{pmatrix} \end{pmatrix}" not in tex
    assert r"\begin{pmatrix}\end{pmatrix}" not in tex
    assert r"(\text{the zero map})" in tex      # the disambiguating parenthetical


# --------------------------------------------------------------------------- #
# MAJOR-3: the Ext arithmetic runs on a genuinely nonzero differential.
# --------------------------------------------------------------------------- #
def test_ext_report_has_nonzero_differential_arithmetic():
    A = _square()
    ev, dims = trace_ext(A, A.simple(1), A.projective(1), 3)
    tex = render_latex(ev, title="Ext(S_1,P_1)")
    # delta^0 = (1;1): a genuinely nonzero connecting map (rank 1), not the old 0-0-0.
    assert r"\begin{pmatrix} 1 \\ 1 \end{pmatrix}" in tex
    assert r"\operatorname{rank}\delta^{0} = 1" in tex
    assert "1 - 1 - 0" in tex and "2 - 1 - 1" in tex     # nonzero rank arithmetic
    # the degree-2 connecting map is the zero map -> "0", disambiguated:
    assert r"\delta^{2} = 0" in tex


def test_render_is_byte_deterministic_on_nonzero_ext():
    A = _square()

    def build():
        ev, _ = trace_ext(A, A.simple(1), A.projective(1), 3)
        return render_latex(ev, title="t")
    assert build() == build()


# --------------------------------------------------------------------------- #
# BLOCKING-1: the preamble lifts MaxMatrixCols; wide matrices scale, >25 elide.
# --------------------------------------------------------------------------- #
def test_preamble_lifts_maxmatrixcols_for_homework_scale_matrix():
    A = truncated_polynomial(13, field=GF(32003))
    M = A.projective(1).radical()                # dim 12, 12x12 arrow action
    M.name = "M"
    ev = trace_module_report(A, M, N=A.simple(1), top=3, with_tau=True, with_decompose=True)
    tex = render_latex(ev, title="dim 12", algebra=A)
    assert r"\usepackage{graphicx}" in tex
    import re
    m = re.search(r"\\setcounter\{MaxMatrixCols\}\{(\d+)\}", tex)
    assert m and int(m.group(1)) >= 12                     # homework-scale (>=12) col matrix
    assert r"\qlmat" in tex                                 # the shrink-to-width wrap is used


def test_preamble_helper_only_sets_counter_when_needed():
    assert r"\setcounter{MaxMatrixCols}" not in "".join(matrix_preamble_lines(10))
    assert r"\setcounter{MaxMatrixCols}{12}" in "".join(matrix_preamble_lines(12))
    # graphicx + the macro are always present (the P/I Loewy displays use \qlmat too):
    assert r"\usepackage{graphicx}" in matrix_preamble_lines(3)


def test_wide_matrix_over_25_cols_is_stated_elided_not_pmatrix():
    from quiverlab import CC
    dom = CC.make_domain([CC.parse_entry(0), CC.parse_entry(1)])
    D = [[dom.zero()] * 30 for _ in range(3)]              # 3 x 30, recorded in full
    e = module_differential(1, "projective", "P", "d_{1}", [1] * 3, [1] * 30, D,
                            3, 30, dom, rank=1)
    assert e.elided is False and e.matrix is not None       # recorded IN FULL
    tex = render_latex([e], title="t")
    assert "shown in full in the HTML/JSON report" in tex   # STATED page-elision
    assert r"\begin{pmatrix}" not in tex                    # the 30-col body is NOT typeset
    assert _events_max_cols([e]) == 0                       # so MaxMatrixCols stays default


# --------------------------------------------------------------------------- #
# Compile smokes (skipped without a toolchain): the samples build to a clean PDF.
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(_latex_engine() is None, reason="no LaTeX toolchain on PATH")
def test_dim12_module_report_compiles_to_pdf(tmp_path):
    A = truncated_polynomial(13, field=GF(32003))
    M = A.projective(1).radical()
    M.name = "M"
    ev = trace_module_report(A, M, N=A.simple(1), top=3)
    refs = resolve_references(references_for(ev))
    tex = render_latex(ev, title="M = rad P_1 over k[x]/(x^13), dim 12",
                       references=refs, algebra=A)
    proc, overfull = _compile(tex, tmp_path)
    assert proc.returncode == 0, proc.stdout.decode("utf-8", "replace")[-2000:]
    assert (tmp_path / "r.pdf").exists()
    assert (tmp_path / "r.pdf").read_bytes()[:5] == b"%PDF-"
    assert overfull == 0, "the dim-12 homework report has overfull boxes"


@pytest.mark.skipif(_latex_engine() is None, reason="no LaTeX toolchain on PATH")
def test_nonzero_ext_report_compiles_to_pdf(tmp_path):
    A = _square()
    ev, _ = trace_ext(A, A.simple(1), A.projective(1), 3)
    refs = resolve_references(references_for(ev))
    tex = render_latex(ev, title="Ext(S_1, P_1) over the square", references=refs, algebra=A)
    proc, overfull = _compile(tex, tmp_path)
    assert proc.returncode == 0, proc.stdout.decode("utf-8", "replace")[-2000:]
    assert (tmp_path / "r.pdf").read_bytes()[:5] == b"%PDF-"
    assert overfull == 0


@pytest.mark.skipif(_latex_engine() is None, reason="no LaTeX toolchain on PATH")
def test_wide_over_25col_report_compiles_with_stated_notice(tmp_path):
    from quiverlab import CC
    dom = CC.make_domain([CC.parse_entry(0), CC.parse_entry(1)])
    D = [[dom.zero()] * 30 for _ in range(4)]
    e = module_differential(1, "projective", "P", "d_{1}", [1] * 4, [1] * 30, D,
                            4, 30, dom, rank=2)
    tex = render_latex([e], title="wide (30-col) differential")
    assert "shown in full in the HTML/JSON report" in tex
    proc, overfull = _compile(tex, tmp_path)
    assert proc.returncode == 0, proc.stdout.decode("utf-8", "replace")[-2000:]
    assert (tmp_path / "r.pdf").read_bytes()[:5] == b"%PDF-"
    assert overfull == 0
