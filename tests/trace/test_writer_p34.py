"""Plan 34 -- the no-code worked-steps HTML report is print-ready + typeset.

The former PDF path (a LaTeX toolchain compiling a page-bounded document) has been
removed; the print-ready HTML report is the sole rendered artifact and is exported
to PDF via the browser's Print -> Save as PDF. What survives here is the HTML
report's print-readiness (print CSS, typeset MathML, self-containment)."""
from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace.render_html import render_html


def _small_events():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    A.hochschild_cohomology(2, trace=tr)
    return A, list(tr)


def _wide_events():
    """HH^*(k[x]/(x^3)) -- its degree-2 differential is a 12-column pmatrix (>10)."""
    A = truncated_polynomial(3, field=CC)
    tr = Trace()
    A.hochschild_cohomology(2, trace=tr)
    return A, list(tr)


# --------------------------------------------------------------------------- #
# The HTML report is print-ready and TYPESET (no apology, MathML, print CSS).
# --------------------------------------------------------------------------- #
def test_html_report_is_print_ready_and_typeset():
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


def test_wide_differential_renders_full_matrix_in_html():
    """A >10-column differential is shown IN FULL in the HTML report (the HTML has
    no page bound, so nothing is elided or scaled away)."""
    _A, ev = _wide_events()
    html = render_html(ev, title="HH")
    assert r"\begin{pmatrix}" in html
    assert "<mtable>" in html


def test_render_html_is_byte_deterministic():
    _A, ev = _small_events()
    refs = (("Refkey2020", "A. Author, A Journal 1 (2020), 1-2."),)
    assert render_html(ev, title="HH", references=refs) == \
        render_html(ev, title="HH", references=refs)
