"""Worked-steps writer: HTML + JSON output contract + the printed one-liner
(spec §3.8). PDF/TeX report output has been removed; the writer produces a
self-contained, print-ready HTML report and its JSON machine record."""
import pathlib

from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace import writer as W
from quiverlab.trace.render_html import render_html


def _events():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    return A, list(tr), table


# A generic (bibtex_key, formatted) fixture; writer/renderer tests do not assert on
# any specific citation, so they stay decoupled from Plan 06's registry.
REFS = (("Refkey2020", "A. Author, A Journal 1 (2020), 1-2."),)


def test_html_is_self_contained_no_js_tex_source():
    A, ev, table = _events()
    html = render_html(ev, title="HH", references=REFS)
    assert "<html" in html.lower()
    # Marco's decision: the HTML report is JavaScript-free and self-contained.
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
    # math is TYPESET as MathML, with the LaTeX source preserved verbatim in an
    # x-tex <annotation> (self-contained, no JS) -- the display renders/prints
    # typeset AND the source stays available for copy/paste (Plan 34).
    assert "<math" in html and "<mtable>" in html          # typeset display
    assert 'encoding="application/x-tex"' in html           # embedded source
    assert r"\begin{pmatrix}" in html and r"\operatorname{rank}" in html
    assert "HH" in html and "Refkey2020" in html


def test_write_trace_produces_html_and_json(tmp_path, monkeypatch):
    A, ev, table = _events()
    printed = {}
    monkeypatch.setattr("builtins.print",
                        lambda *a, **k: printed.setdefault("line", " ".join(map(str, a))))
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS,
                         out_dir=str(tmp_path))
    p = pathlib.Path(path)
    assert path.endswith(".html") and p.exists()
    # the JSON machine record sits beside the HTML, sharing the stem
    j = p.with_suffix(".json")
    assert j.exists()
    # PDF/TeX report artifacts are never produced
    assert not p.with_suffix(".pdf").exists()
    assert not p.with_suffix(".tex").exists()
    # the loud one-liner points at the HTML report
    assert printed["line"].startswith("Worked steps: ") and ".html" in printed["line"]


def test_output_dir_and_filename_contract(tmp_path, monkeypatch):
    A, ev, table = _events()
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS,
                         out_dir=str(tmp_path))
    p = pathlib.Path(path)
    assert p.parent == tmp_path
    assert p.name.startswith("HH_") and p.suffix == ".html"  # <kind>_<hash>.<ext>
    # hash must stay wide (guards the [:12] widening against a silent [:4] revert)
    assert len(p.stem.split("_")[-1]) >= 8


def test_render_html_is_byte_deterministic():
    A, ev, table = _events()
    assert render_html(ev, title="HH", references=REFS) == \
        render_html(ev, title="HH", references=REFS)
