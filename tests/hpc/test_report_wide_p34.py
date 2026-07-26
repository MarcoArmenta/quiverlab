"""Plan 34 (post-critique) -- the hpc result.json -> report path has the SAME page-bounded
matrix handling as the trace worked-steps PDF:

  * BLOCKING-1 (hpc twin): render_latex lifted no MaxMatrixCols ceiling, so a dim>10
    module's action matrix (an 11..25-column pmatrix) aborted the compile; it now emits
    \\setcounter{MaxMatrixCols} + graphicx \\qlmat, and a matrix past 25 rows/cols is
    STATED-elided (shown in full in the HTML/JSON report).
  * the render_html apology ("Math is shown as TeX source ...") is replaced by an honest
    print-ready hint."""
import shutil
import subprocess

import pytest

from quiverlab.hpc import report


def _latex_engine():
    for engine in ("pdflatex", "tectonic"):
        if shutil.which(engine):
            return engine
    return None


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
# _pmatrix_latex: <=25 -> pmatrix, >25 -> stated elision.
# --------------------------------------------------------------------------- #
def test_pmatrix_latex_elides_past_25():
    assert r"\begin{pmatrix}" in report._pmatrix_latex([[0] * 12 for _ in range(12)])
    wide = report._pmatrix_latex([[0] * 30 for _ in range(3)])
    assert r"\begin{pmatrix}" not in wide
    assert "shown in full in the HTML/JSON report" in wide


def test_tex_matrix_cols_counts_widest_row():
    assert report._tex_matrix_cols(r"\begin{pmatrix} 1 & 2 & 3 \end{pmatrix}") == 3
    assert report._tex_matrix_cols(r"\begin{bmatrix} 1 \\ 2 \end{bmatrix}") == 1
    assert report._tex_matrix_cols(r"no matrix here") == 0


# --------------------------------------------------------------------------- #
# render_latex: MaxMatrixCols + \qlmat for the 12-col map, elision for the 30-col map.
# --------------------------------------------------------------------------- #
def test_render_latex_lifts_ceiling_and_elides_wide():
    tex = report.render_latex(_result_with_maps(12, 30))
    assert r"\usepackage{graphicx}" in tex
    assert r"\setcounter{MaxMatrixCols}{12}" in tex        # the 12-col map is typeset
    assert r"\qlmat" in tex                                 # ...and scaled to the page
    assert "shown in full in the HTML/JSON report" in tex   # the 30-col map is elided


def test_render_html_apology_replaced_with_print_hint():
    html = report.render_html(_result_with_maps(12, 3))
    assert "Math is shown as TeX source" not in html        # the apology is gone
    assert "print-ready" in html.lower() and "PDF" in html


# --------------------------------------------------------------------------- #
# Compile smoke (skipped without a toolchain).
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(_latex_engine() is None, reason="no LaTeX toolchain on PATH")
def test_hpc_report_with_wide_maps_compiles(tmp_path):
    tex = report.render_latex(_result_with_maps(12, 30))
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
    assert proc.returncode == 0, proc.stdout.decode("utf-8", "replace")[-2000:]
    assert (tmp_path / "r.pdf").read_bytes()[:5] == b"%PDF-"
