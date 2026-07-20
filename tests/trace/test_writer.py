"""LaTeX/HTML renderers + renderer selection + output-path contract (spec §3.8)."""
import pathlib

import pytest

from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace import writer as W
from quiverlab.trace.render_latex import render_latex
from quiverlab.trace.render_html import render_html


def _events():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    return A, list(tr), table


# A generic (bibtex_key, formatted) fixture; writer/renderer tests do not assert on
# any specific citation, so they stay decoupled from Plan 06's registry.
REFS = (("Refkey2020", "A. Author, A Journal 1 (2020), 1-2."),)


def test_latex_has_matrices_dims_and_bibliography():
    A, ev, table = _events()
    src = render_latex(ev, title="HH", references=REFS)
    assert r"\begin{pmatrix}" in src
    assert "HH^{0} = 2" in src or "HH^0 = 2" in src
    assert r"\begin{thebibliography}" in src and "Refkey2020" in src
    assert src.startswith(r"\documentclass")


def test_html_is_self_contained_no_js_tex_source():
    A, ev, table = _events()
    html = render_html(ev, title="HH", references=REFS)
    assert "<html" in html.lower()
    # Marco's decision: the HTML fallback is JavaScript-free and self-contained.
    assert "<script" not in html.lower() and "mathjax" not in html.lower()
    assert "polyfill.io" not in html and "jsdelivr" not in html
    # ...and, more generally, no network fetch, no external resource, no inline JS:
    assert "http://" not in html and "https://" not in html
    assert "<link" not in html.lower()
    assert "<iframe" not in html.lower()
    # no inline event handlers (onload=, onclick=, onerror=, ... i.e. any " on...=")
    assert "onload=" not in html.lower()
    import re as _re
    assert not _re.search(r"\son[a-z]+\s*=", html.lower())
    # math is shown as readable TeX source (not typeset)
    assert r"\begin{pmatrix}" in html and r"\operatorname{rank}" in html
    assert "HH" in html and "Refkey2020" in html


def test_selection_prefers_pdf_when_toolchain_present(tmp_path, monkeypatch):
    A, ev, table = _events()
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")
    compiled = {}
    def fake_compile(tex, out_pdf, engine):
        pathlib.Path(out_pdf).write_bytes(b"%PDF-1.5 fake\n%%EOF\n")
        compiled["engine"] = engine
        return 1  # page count
    monkeypatch.setattr(W, "_compile_pdf", fake_compile)
    printed = {}
    monkeypatch.setattr("builtins.print", lambda *a, **k: printed.setdefault("line", " ".join(map(str, a))))
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS, out_dir=str(tmp_path))
    assert path.endswith(".pdf") and pathlib.Path(path).exists()
    assert compiled["engine"] == "tectonic"
    assert printed["line"].startswith("Worked steps: ") and ".pdf" in printed["line"]


def test_selection_falls_back_to_html_with_loud_message(tmp_path, monkeypatch):
    A, ev, table = _events()
    monkeypatch.setattr(W, "have_latex", lambda: None)   # no toolchain
    printed = {}
    monkeypatch.setattr("builtins.print", lambda *a, **k: printed.setdefault("line", " ".join(map(str, a))))
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS, out_dir=str(tmp_path))
    assert path.endswith(".html") and pathlib.Path(path).exists()
    line = printed["line"]
    # loud one-line explanation of BOTH facts: no toolchain -> HTML written
    assert line.startswith("Worked steps: ") and ".html" in line
    assert "no LaTeX toolchain found" in line and ("pdflatex" in line or "tectonic" in line)


def test_html_fallback_when_compile_fails(tmp_path, monkeypatch):
    """A toolchain IS present but compilation raises: the message must be honest
    ('compilation failed'), never 'no toolchain found'."""
    A, ev, table = _events()
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")
    def boom(tex, out_pdf, engine):
        raise RuntimeError("tectonic exploded")
    monkeypatch.setattr(W, "_compile_pdf", boom)
    printed = {}
    monkeypatch.setattr("builtins.print",
                        lambda *a, **k: printed.setdefault("line", " ".join(map(str, a))))
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS, out_dir=str(tmp_path))
    assert path.endswith(".html") and pathlib.Path(path).exists()
    line = printed["line"]
    assert "compilation failed" in line
    assert "no LaTeX toolchain found" not in line


def test_output_dir_and_filename_contract(tmp_path, monkeypatch):
    A, ev, table = _events()
    monkeypatch.setattr(W, "have_latex", lambda: None)
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS, out_dir=str(tmp_path))
    p = pathlib.Path(path)
    assert p.parent == tmp_path
    assert p.name.startswith("HH_") and p.suffix == ".html"  # <kind>_<hash>.<ext>
    # hash must stay wide (guards the [:12] widening against a silent [:4] revert)
    assert len(p.stem.split("_")[-1]) >= 8
