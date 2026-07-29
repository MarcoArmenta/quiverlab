"""Regression battery for the F2 trace-layer fixes (Plan 34):

  * the CS worked-steps Result shows the ENGINE's dims with the correct HH^/HH_
    variance (the headline bug: CS used to record only ResolutionTerms, so the
    renderers re-derived the cochain-SPACE dims and mislabelled the variance --
    QuantumCI(-1)/GF(5) cohomology rendered ``HH_0=4 HH_1=8 HH_2=12 HH_3=16``
    instead of the engine's ``HH^ = [4,4,5,6]``);
  * the HH drift gate refuses a report whose event-derived dims disagree with
    the engine's;
  * the default GF(p) fast-engine report (empty worked steps) still gets an
    authoritative Result line + an honest note;
  * a GF(p^n) matrix entry renders readably (``x^1``) -- never the raw internal
    coefficient tuple ``(0, 1)``;
  * the renderers refuse a foreign (non-event) object loudly.
"""
import json
import pathlib

import pytest

from quiverlab import QuantumCI, truncated_polynomial, GF, CC
from quiverlab.errors import QuiverlabError
from quiverlab.trace.recorder import Trace, _cell_str, ext_degree
from quiverlab.trace import writer as W
from quiverlab.trace.render_text import render_text, derive_dims, _dims_kind
from quiverlab.trace.render_html import render_html
from quiverlab.trace.render_json import render_json
from quiverlab.trace.events import RankStep, ResolutionTerm, Dispatch, ResultDims


def _write_html(events, table, algebra, kind, monkeypatch, tmp_path):
    """Run the real writer (HTML + JSON only); return (html, json)."""
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    path = W.write_trace(list(events), table, algebra=algebra, kind=kind, top=table.top,
                         out_dir=str(tmp_path))
    p = pathlib.Path(path)
    assert p.suffix == ".html"
    return p.read_text(encoding="utf-8"), json.loads(p.with_suffix(".json").read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
# Fix 1: the CS Result -- correct dims AND correct variance, across engines/sides/fields.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("engine", ["bar", "cs"])
@pytest.mark.parametrize("side", ["coh", "hom"])
@pytest.mark.parametrize("field", [GF(5), CC])
def test_rendered_result_equals_engine_dims_with_correct_variance(
        engine, side, field, tmp_path, monkeypatch):
    A = QuantumCI(-1, field=field)
    top = 3
    tr = Trace()
    if side == "coh":
        table = A.hochschild_cohomology(top, engine=engine, trace=tr)
        want_kind, wrong_kind = "HH^", "HH_"
    else:
        table = A.hochschild_homology(top, engine=engine, trace=tr)
        want_kind, wrong_kind = "HH_", "HH^"

    ev = list(tr)
    # Root-cause fix (renderer-agnostic): the CS trace now carries the per-degree
    # RankSteps, so the event-derived dims and variance match the engine EXACTLY.
    assert derive_dims(ev) == table.dims
    assert _dims_kind(ev) == table.kind == want_kind

    # And the shipped artifact (HTML + JSON) shows the engine's dims with that variance.
    html, obj = _write_html(ev, table, A, want_kind, monkeypatch, tmp_path)
    for i, d in enumerate(table.dims):
        assert r"%s{%d} = %d" % (want_kind, i, d) in html
    assert wrong_kind + "{" not in html                       # never the swapped label
    (rd,) = [e for e in obj["events"] if e["type"] == "ResultDims"]
    assert rd["dims"] == table.dims and rd["kind"] == want_kind


def test_quantumci_gf5_cohomology_headline_value(tmp_path, monkeypatch):
    """The exact bug in the task: QuantumCI(-1)/GF(5) cohomology must render
    HH^ = [4,4,5,6], never the pre-fix HH_ = [4,8,12,16] (the cochain-space dims)."""
    A = QuantumCI(-1, field=GF(5))
    tr = Trace()
    table = A.hochschild_cohomology(3, engine="cs", trace=tr)
    assert table.dims == [4, 4, 5, 6] and table.kind == "HH^"
    html, obj = _write_html(list(tr), table, A, "HH^", monkeypatch, tmp_path)
    for i, d in zip(range(4), [4, 4, 5, 6]):
        assert r"HH^{%d} = %d" % (i, d) in html
    # the pre-fix wrong render must be gone, numbers AND label
    for bad in ("HH_{1} = 8", "HH_{2} = 12", "HH_{3} = 16"):
        assert bad not in html
    txt = render_text(list(tr) + [ResultDims("HH^", [4, 4, 5, 6])], title="t")
    assert "Result: HH^0 = 4   HH^1 = 4   HH^2 = 5   HH^3 = 6" in txt


# --------------------------------------------------------------------------- #
# Fix 2: the HH drift gate refuses a report that misstates the computation.
# --------------------------------------------------------------------------- #
def test_drift_gate_raises_on_dim_mismatch():
    from quiverlab.trace.writer import _authoritative_result
    from quiverlab.hochschild.table import HHTable

    # A minimal cohomology trace whose derived dim is 1 at degree 0 (C_0=1, rank 0).
    ev = [ResolutionTerm(degree=0, n_generators=1, collapsed_dim=1),
          RankStep(degree=0, side="cochain", nrows=0, ncols=1, rank=0, field="GF(5)",
                   matrix=[])]
    assert derive_dims(ev) == [1]
    good = HHTable([1], "HH^", "A")
    # matching table: no raise, ResultDims appended
    out = _authoritative_result(list(ev), good)
    assert any(isinstance(e, ResultDims) for e in out)
    # a LYING table (wrong dims) must be refused loudly
    with pytest.raises(QuiverlabError):
        _authoritative_result(list(ev), HHTable([999], "HH^", "A"))
    # a swapped VARIANCE must be refused loudly
    with pytest.raises(QuiverlabError):
        _authoritative_result(list(ev), HHTable([1], "HH_", "A"))


# --------------------------------------------------------------------------- #
# Fix 3: the default GF(p) fast-engine report is never empty.
# --------------------------------------------------------------------------- #
def test_fast_gfp_report_has_authoritative_result_and_honest_note(tmp_path, monkeypatch):
    A = truncated_polynomial(2, field=GF(2))          # prime field -> fast engine
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    # the fast engine records only a Dispatch (no per-degree worked steps)
    assert [type(e).__name__ for e in tr] == ["Dispatch"]
    html, obj = _write_html(list(tr), table, A, "HH^", monkeypatch, tmp_path)
    # ...yet the report has a Result line with the engine's dims + an HONEST note.
    assert "<h2 id='result'>Result</h2>" in html
    for i, d in enumerate(table.dims):
        assert r"HH^{%d} = %d" % (i, d) in html
    assert "records no per-degree worked steps" in html
    (rd,) = [e for e in obj["events"] if e["type"] == "ResultDims"]
    assert rd["dims"] == table.dims and rd["note"]


# --------------------------------------------------------------------------- #
# Fix 5: a GF(p^n) matrix entry renders readably, never the raw tuple.
# --------------------------------------------------------------------------- #
def test_gf4_cell_str_is_readable_not_a_tuple():
    F = GF(4)
    dom = F.make_domain([F.parse_entry(0), F.parse_entry(1)])
    assert str(dom.one()) == "(1, 0)"                 # the raw internal form leaks a tuple
    assert _cell_str(dom.one(), dom) == "1"           # ...but the recorder renders it readably
    assert _cell_str((0, 1), dom) == "x^1"
    assert _cell_str((1, 1), dom) == "1 + x^1"
    # int/fraction entries stay byte-identical to plain str (no regression for GF(p)/CC).
    assert _cell_str(3, dom) == "3" and _cell_str("-1/2", dom) == "-1/2"


def test_gf4_trace_leaks_no_tuple_in_any_render():
    A = truncated_polynomial(3, field=GF(4))          # GF(p^n): not prime -> bar engine
    tr = Trace()
    table = A.hochschild_cohomology(3, trace=tr)
    ev = list(tr)
    # at least one recorded entry came from a GF(4) tuple (else the test is vacuous):
    cells = [c for e in ev if isinstance(e, RankStep) and e.matrix
             for row in e.matrix for c in row]
    assert cells and all("(" not in c and ")" not in c for c in cells)
    for render in (render_text, render_html, render_json):
        s = render(ev, title="GF4")
        assert "(0, 1)" not in s and "(1, 0)" not in s


# --------------------------------------------------------------------------- #
# Fix 6: the Ext/Tor footer shows ALL runs (an Ext run AND a Tor run), not just the first.
# --------------------------------------------------------------------------- #
def test_ext_then_tor_footer_shows_both_runs():
    dom = GF(2)
    D = [[dom.one()]]
    ev = [ext_degree(0, "Ext", 1, 0, 0, 1, D, 1, 1, dom),
          ext_degree(1, "Ext", 1, 1, 0, 0, D, 1, 1, dom),
          ext_degree(0, "Tor", 2, 0, 0, 2, D, 1, 1, dom),   # a second run
          ext_degree(1, "Tor", 1, 0, 0, 1, D, 1, 1, dom)]
    txt = render_text(ev, title="both")
    assert "Result: Ext^0 = 1   Ext^1 = 0" in txt          # the Ext run (superscript)
    assert "Result: Tor_0 = 2   Tor_1 = 1" in txt          # the Tor run (subscript)
    html = render_html(ev, title="both")
    assert r"\operatorname{Ext}^{0} = 1" in html
    assert r"\operatorname{Tor}_{0} = 2" in html and r"\operatorname{Tor}_{1} = 1" in html


# --------------------------------------------------------------------------- #
# Fix 7: a full event buffer surfaces a "steps elided" note in the rendered report.
# --------------------------------------------------------------------------- #
def test_buffer_full_note_is_rendered():
    tr = Trace(max_events=3)
    for i in range(10):
        tr.append(Dispatch(route=str(i), reason="", n_relations=0))
    txt = render_text(list(tr), title="t")
    html = render_html(list(tr), title="t")
    assert "later steps were elided" in txt
    assert "later steps were elided" in html


# --------------------------------------------------------------------------- #
# Fix 4: the renderers refuse a foreign (non-event) object loudly (ALL_EVENTS gate).
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("render", [render_text, render_html, render_json])
def test_foreign_object_in_stream_is_refused(render):
    good = [Dispatch(route="normalized bar complex", reason="x", n_relations=0)]
    # a stray (events, result) tuple element -- the classic caller bug
    bad = good + [("events", "result")]
    with pytest.raises(QuiverlabError) as ei:
        render(bad, title="t")
    assert "tuple" in str(ei.value)
    # a plain string leaks too
    with pytest.raises(QuiverlabError):
        render(good + ["not an event"], title="t")
    # the clean stream still renders
    render(good, title="t")
