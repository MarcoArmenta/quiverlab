"""Plan 34 -- the no-code worked-steps HTML report is print-ready + typeset.

The former PDF path (a LaTeX toolchain compiling a page-bounded document) has been
removed; the print-ready HTML report is the sole rendered artifact and is exported
to PDF via the browser's Print -> Save as PDF. What survives here is the HTML
report's print-readiness (print CSS, typeset MathML, self-containment)."""
from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace.render_html import render_html
from tests.trace._matrix_grid import grids, has_grid, grid_indices


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
    assert "<math" in html and grids(html)
    assert 'encoding="application/x-tex"' in html          # source kept for copy/paste
    # the old apology is gone; and so is the print-to-PDF suggestion -- the
    # deliverables are HTML + JSON ONLY (Marco 2026-07-28), stated up front.
    assert "Math is shown as TeX source" not in html
    assert "Save as PDF" not in html
    assert "HTML page and the JSON records" in html
    # the JSON-record note documents the machine deliverable's structure
    assert "quiverlab_trace_schema" in html
    # still fully self-contained (no JS, no external resource, no network)
    assert "<script" not in html.lower()
    assert "http://" not in html and "https://" not in html
    assert "<link" not in html.lower()


def test_wide_differential_over_cap_points_to_json():
    """Marco 2026-08-02: a differential exceeding the 20-line display cap on either
    dimension is stated as a size note pointing at the JSON, not a full grid (the
    complete matrix always lives in trace.json). The HH^*(k[x]/(x^3)) degree-2
    differential is 24x12 (24 rows > 20), so it elides to the note."""
    _A, ev = _wide_events()
    html = render_html(ev, title="HH")
    assert "24 rows and 12 columns exceed the 20-line display cap" in html
    assert "JSON record" in html


def test_render_html_is_byte_deterministic():
    _A, ev = _small_events()
    refs = (("Refkey2020", "A. Author, A Journal 1 (2020), 1-2."),)
    assert render_html(ev, title="HH", references=refs) == \
        render_html(ev, title="HH", references=refs)
