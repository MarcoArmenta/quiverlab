"""Plan 34 -- the hpc result.json -> report path shows wide matrices sensibly and the
HTML report is print-ready:

  * ``_matrix_html``: a matrix up to 25 columns is a rendered HTML table; past 25
    rows/cols it is STATED-elided (shown in full in the JSON result) so an oversized
    module action matrix does not dominate the page.
  * the HTML report shows RENDERED math (tables, sub/superscripts), never LaTeX
    source, and no removed .tex/PDF compile step is advertised.

(PDF/TeX report output has been removed; the report renders to HTML/text/JSON only.)"""
from quiverlab.hpc import report


def _result_with_maps(cols_rad, cols_soc):
    """A minimal result carrying a rad/top/soc block whose radical map is ``cols_rad``
    wide and socle map ``cols_soc`` wide (square). ``algebra={}`` is unbuildable, so the
    descriptive P/I section is skipped -- keeping the fixture about the matrices."""
    def sq(n):
        m = [[0] * n for _ in range(n)]
        if n:
            m[0][0] = 1
        return m
    return {
        "quiverlab_version": "0.1.0-test",
        "result_schema": 1,
        "algebra": {},
        "results": {
            "rad_top_soc": {
                "side": "right",
                "radical": {"dims": {"1": cols_rad}, "maps": {"x": sq(cols_rad)}},
                "top": {"dims": {"1": 1}, "maps": {}},
                "socle": {"dims": {"1": cols_soc}, "maps": {"y": sq(cols_soc)}},
            }
        },
    }


# --------------------------------------------------------------------------- #
# _matrix_html: <=25 -> rendered table, >25 -> stated elision.
# --------------------------------------------------------------------------- #
def test_matrix_html_elides_past_25():
    assert "<table class='matrix'>" in report._matrix_html([[0] * 12 for _ in range(12)])
    wide = report._matrix_html([[0] * 30 for _ in range(3)])
    assert "<table" not in wide
    assert "shown in full in the JSON result" in wide


def test_render_html_shows_tables_not_latex_source():
    html = report.render_html(_result_with_maps(12, 3))
    # matrices are RENDERED, labeled by their arrow, never LaTeX source:
    assert "<table class='matrix'>" in html
    assert "action of arrow <code>x</code>" in html
    assert "pmatrix" not in html and "\\begin{" not in html


def test_render_html_apology_replaced_with_print_hint():
    html = report.render_html(_result_with_maps(12, 3))
    assert "Math is shown as TeX source" not in html        # the apology is gone
    assert "print-ready" in html.lower() and "PDF" in html
    # the removed .tex/PDF-compile step is not advertised:
    assert "pdflatex" not in html and "tectonic" not in html
