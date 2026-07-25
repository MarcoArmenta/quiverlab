"""Renderer goldens: token asserts on the .txt/.html/.tex output, PDF only when a
LaTeX toolchain is on PATH (else cleanly skipped), and the future-``result_schema``
refusal."""
import shutil

import pytest

from quiverlab.hpc import report
from quiverlab.hpc.spec import RESULT_SCHEMA
from quiverlab.hpc.spec import run as spec_run

_CFG = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["hh_cohomology:0..3", "cartan", "dimension"],
    "artifacts": {"pdf": False, "tikz": False},
}


@pytest.fixture()
def result(tmp_path):
    return spec_run(_CFG, tmp_path / "compute", result_schema=RESULT_SCHEMA)


def test_text_tokens(result, tmp_path):
    out, fmt = report.render(result, tmp_path / "r.txt", fmt="txt")
    assert fmt == "txt"
    text = out.read_text(encoding="utf-8")
    assert "quiverlab report" in text
    assert "Hochschild cohomology" in text
    assert "HH^0 = 3" in text
    assert "Cartan matrix" in text
    assert "dim A = 3" in text


def test_html_tokens_and_escaping(result, tmp_path):
    out, fmt = report.render(result, tmp_path / "r.html", fmt="html")
    assert fmt == "html"
    html = out.read_text(encoding="utf-8")
    assert "<!doctype html>" in html
    assert "<h1>quiverlab report</h1>" in html
    assert "Hochschild cohomology" in html
    assert "<script" not in html                  # no JavaScript, ever
    # The pmatrix TeX source is shown escaped, not raw-injected.
    assert "pmatrix" in html


def test_latex_tokens(result):
    tex = report.render_latex(result)
    assert r"\documentclass{article}" in tex
    assert r"\begin{document}" in tex and r"\end{document}" in tex
    assert "Hochschild cohomology" in tex
    assert r"\begin{pmatrix}" in tex


def test_auto_falls_back_to_html_without_toolchain(result, tmp_path, monkeypatch):
    monkeypatch.setattr(report, "have_latex", lambda: None)
    out, fmt = report.render(result, tmp_path / "r", fmt="auto")
    assert fmt == "html"


def test_pdf_requires_toolchain(result, tmp_path, monkeypatch):
    monkeypatch.setattr(report, "have_latex", lambda: None)
    with pytest.raises(report.ReportError):
        report.render(result, tmp_path / "r.pdf", fmt="pdf")


def test_pdf_when_toolchain_present(result, tmp_path):
    if not (shutil.which("tectonic") or shutil.which("pdflatex")):
        pytest.skip("no LaTeX toolchain (tectonic/pdflatex) on PATH")
    out, fmt = report.render(result, tmp_path / "r.pdf", fmt="pdf")
    assert fmt == "pdf"
    assert out.read_bytes()[:4] == b"%PDF"


def test_future_result_schema_refused(result, tmp_path):
    result["result_schema"] = RESULT_SCHEMA + 1
    with pytest.raises(report.ResultSchemaError):
        report.render(result, tmp_path / "r.txt", fmt="txt")


def test_version_skew_warns(result, tmp_path):
    result["quiverlab_version"] = "0.0.1-ancient"
    warnings = []
    report.render(result, tmp_path / "r.txt", fmt="txt", on_warn=warnings.append)
    assert warnings and "0.0.1-ancient" in warnings[0]
