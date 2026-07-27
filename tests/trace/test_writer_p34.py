"""Plan 34 -- the worked-steps report renders to a REAL PDF, and the no-LaTeX HTML
fallback is print-ready + typeset.

Regression for Marco's "the print report PDF does not give a rendered PDF": a
differential wider than amsmath's default MaxMatrixCols (10) aborted the LaTeX
compile ("Extra alignment tab has been changed to \\cr"), so write_trace silently
degraded to the HTML fallback even with a toolchain present. The writer now lifts
the column ceiling on the persisted .tex AND the compiled PDF.
"""
import pathlib

import pytest

from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace import writer as W
from quiverlab.trace.render_html import render_html


def _wide_events():
    """HH^*(k[x]/(x^3)) -- its degree-2 differential is a 12-column pmatrix (>10),
    exactly the case that used to abort the amsmath compile."""
    A = truncated_polynomial(3, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    return A, list(tr), table


def _small_events():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    A.hochschild_cohomology(2, trace=tr)
    return A, list(tr)


# --------------------------------------------------------------------------- #
# MaxMatrixCols: the compile-abort root cause.
# --------------------------------------------------------------------------- #
def test_needed_matrix_cols_detects_wide_differential():
    _A, ev, _t = _wide_events()
    assert W._needed_matrix_cols(ev) > 10, \
        "the k[x]/(x^3) HH^2 differential should exceed amsmath's default of 10"


def test_widen_matrix_cols_injects_before_document_idempotently():
    tex = (r"\documentclass{article}" "\n" r"\usepackage{amsmath}" "\n"
           r"\begin{document}" "\n" r"X" "\n" r"\end{document}" "\n")
    out = W._widen_matrix_cols(tex, 12)
    assert r"\setcounter{MaxMatrixCols}{12}" in out
    assert out.index("MaxMatrixCols") < out.index(r"\begin{document}")
    assert W._widen_matrix_cols(out, 20) == out          # idempotent: never doubled
    assert W._widen_matrix_cols(tex, 10) == tex          # no-op at/under the default


def test_persisted_tex_lifts_maxmatrixcols_for_wide_matrix(tmp_path, monkeypatch):
    A, ev, table = _wide_events()
    monkeypatch.setattr(W, "have_latex", lambda: None)          # force the HTML branch
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, out_dir=str(tmp_path))
    tex = pathlib.Path(path).with_suffix(".tex").read_text()
    assert r"\setcounter{MaxMatrixCols}" in tex, \
        "the compile-it-yourself .tex must build standalone (wide differential)"


# --------------------------------------------------------------------------- #
# The ladder: PDF when an engine is present; the widened tex reaches the compiler.
# --------------------------------------------------------------------------- #
def test_ladder_picks_pdf_and_widens_the_compiled_tex(tmp_path, monkeypatch):
    A, ev, table = _wide_events()
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")
    seen = {}

    def fake_compile(tex, out_pdf, engine):
        seen["tex"] = tex
        pathlib.Path(out_pdf).write_bytes(b"%PDF-1.5 fake\n%%EOF\n")
        return 1

    monkeypatch.setattr(W, "_compile_pdf", fake_compile)
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, out_dir=str(tmp_path))
    assert path.endswith(".pdf") and pathlib.Path(path).exists()
    assert r"\setcounter{MaxMatrixCols}" in seen["tex"], \
        "the tex handed to the compiler must carry the lifted column ceiling"


@pytest.mark.skipif(W.have_latex() is None, reason="no LaTeX toolchain on PATH")
def test_real_latex_compiles_wide_matrix_to_pdf(tmp_path):
    """The bug, end to end with a REAL engine: a >10-column differential compiles to
    a PDF instead of silently degrading to HTML."""
    A, ev, table = _wide_events()
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, out_dir=str(tmp_path))
    assert path.endswith(".pdf"), "a wide-matrix trace must render a PDF, not HTML"
    assert pathlib.Path(path).read_bytes()[:5] == b"%PDF-", "must be a real PDF"


# --------------------------------------------------------------------------- #
# The HTML fallback is print-ready and TYPESET (no apology, MathML, print CSS).
# --------------------------------------------------------------------------- #
def test_html_fallback_is_print_ready_and_typeset():
    _A, ev = _small_events()
    html = render_html(ev, title="HH")
    # print CSS: page setup + page-break control
    assert "@media print" in html and "@page" in html
    assert "break-inside" in html
    # math is typeset as MathML (not raw source dumped as text)
    assert "<math" in html and "<mtable>" in html
    assert 'encoding="application/x-tex"' in html          # source kept for copy/paste
    # the old apology is gone; an honest print-to-PDF hint replaces it
    assert "Math is shown as TeX source" not in html
    assert "Print" in html and "PDF" in html
    # still fully self-contained (no JS, no external resource, no network)
    assert "<script" not in html.lower()
    assert "http://" not in html and "https://" not in html
    assert "<link" not in html.lower()


def test_render_html_is_byte_deterministic():
    _A, ev = _small_events()
    refs = (("Refkey2020", "A. Author, A Journal 1 (2020), 1-2."),)
    assert render_html(ev, title="HH", references=refs) == \
        render_html(ev, title="HH", references=refs)
