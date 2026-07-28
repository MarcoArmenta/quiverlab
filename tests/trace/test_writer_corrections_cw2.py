"""Correction-worker CW2 follow-ups to the Plan-34 F2 trace-layer fixes.

Two self-contradictions the F2 critic found in the shipped trace layer:

  * Correction 1 -- the HH DRIFT GATE (writer ``_authoritative_result``) compared the
    event-derived dims to the engine's ``table.dims`` and HARD-RAISED on mismatch. Under
    the MAX_EVENTS buffer cap (Fix 7) the surviving event list is TRUNCATED, so
    ``derive_dims`` reconstructs a short/partial dim list -- and the gate raised, defeating
    the cap (whose point is to let an over-long report still ship with an elision note).
    The gate now SKIPS the derive-vs-table cross-check when the buffer elided (detected via
    the recorder's reserved BUFFER_FULL StepNote); the authoritative ResultDims still ships.

  * Correction 2 -- the auto->CS depth fallback appends CS worked steps onto the PARTIAL
    bar events already in the recorder, so one stream carries BOTH the abandoned bar
    differentials AND the CS ones (the drift gate survived only by last-wins luck). The
    writer now drops the abandoned bar worked steps (keeping both Dispatch lines -- the
    routing story) so the report is a single coherent CS worked-steps set, and the gate
    passes by construction.
"""
import io
import json
import contextlib
import pathlib

import pytest

from quiverlab import truncated_polynomial, CC
from quiverlab.errors import QuiverlabError
from quiverlab.hochschild.table import HHTable
from quiverlab.trace.recorder import Trace, BUFFER_FULL_PREFIX
from quiverlab.trace import writer as W
from quiverlab.trace.writer import (
    _drop_superseded, _authoritative_result, _superseding_dispatch_index,
)
from quiverlab.trace.render_text import derive_dims
from quiverlab.trace.events import (
    Dispatch, ResolutionTerm, RankStep, StepNote, ResultDims,
)


def _write_html(events, table, algebra, kind, top, monkeypatch, tmp_path):
    """Run the real writer (HTML + JSON only); return (html, json)."""
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    path = W.write_trace(list(events), table, algebra=algebra, kind=kind, top=top,
                         out_dir=str(tmp_path))
    p = pathlib.Path(path)
    assert p.suffix == ".html"
    return p.read_text(encoding="utf-8"), json.loads(p.with_suffix(".json").read_text(encoding="utf-8"))


def _fallback_stream():
    """A REAL auto->CS depth fallback: a tiny max_cells trips the bar wall at degree 1,
    then auto reroutes to Chouhy-Solotar. The recorder holds the abandoned bar steps
    (6x3, 12x6) AND the CS steps (all 3x3) -- distinguishable by shape."""
    A = truncated_polynomial(3, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(4, max_cells=200, trace=tr)
    return A, list(tr), table


# --------------------------------------------------------------------------- #
# Correction 2: the auto->CS fallback report is a single coherent worked-steps set.
# --------------------------------------------------------------------------- #
def test_fallback_stream_really_carries_both_bar_and_cs_steps():
    """Guard the premise: without the drop, one stream carries BOTH routes' steps."""
    _A, ev, _t = _fallback_stream()
    assert _superseding_dispatch_index(ev) is not None            # two Dispatches
    ranks = [(e.nrows, e.ncols) for e in ev if isinstance(e, RankStep)]
    assert (6, 3) in ranks and (12, 6) in ranks                   # abandoned bar shapes
    assert (3, 3) in ranks                                        # CS shapes


def test_drop_superseded_removes_abandoned_bar_steps():
    _A, ev, _t = _fallback_stream()
    cleaned = _drop_superseded(ev)
    # both Dispatch lines survive (the routing story), abandoned bar steps are gone
    assert sum(isinstance(e, Dispatch) for e in cleaned) == 2
    shapes = [(e.nrows, e.ncols) for e in cleaned if isinstance(e, RankStep)]
    assert shapes and all(s == (3, 3) for s in shapes), \
        "only the CS 3x3 differentials must remain"
    # every CS degree appears exactly once (no bar/CS degree collision)
    degs = [e.degree for e in cleaned if isinstance(e, ResolutionTerm)]
    assert degs == sorted(set(degs)) and degs == [0, 1, 2, 3, 4]


def test_drift_gate_passes_by_construction_on_cleaned_fallback():
    """After the drop the gate derives over the CS steps ALONE -- it matches the engine
    without relying on last-wins luck, and appends the authoritative ResultDims."""
    _A, ev, table = _fallback_stream()
    cleaned = _drop_superseded(ev)
    assert derive_dims(cleaned) == table.dims                     # not by luck
    out = _authoritative_result(cleaned, table)                  # must not raise
    (rd,) = [e for e in out if isinstance(e, ResultDims)]
    assert rd.dims == table.dims and rd.kind == table.kind


def test_fallback_report_ships_single_coherent_set(tmp_path, monkeypatch):
    A, ev, table = _fallback_stream()
    html, obj = _write_html(ev, table, A, "HH^", table.top, monkeypatch, tmp_path)
    # the JSON machine record (lists EVERY event) no longer carries the bar wide matrices
    rank_rows = [e for e in obj["events"] if e["type"] == "RankStep"]
    assert rank_rows and all(e["nrows"] == 3 and e["ncols"] == 3 for e in rank_rows), \
        "abandoned bar differentials must not leak into the artifact"
    # the authoritative Result is correct AND the superseding routing story is present
    for i, d in enumerate(table.dims):
        assert r"HH^{%d} = %d" % (i, d) in html
    assert "chouhy-solotar" in html.lower() or "Chouhy-Solotar" in html
    (rd,) = [e for e in obj["events"] if e["type"] == "ResultDims"]
    assert rd["dims"] == table.dims and rd["kind"] == "HH^"


def test_drop_superseded_is_noop_for_normal_and_module_streams():
    # a normal single-Dispatch HH stream is returned UNCHANGED (byte-stable goldens)
    normal = [Dispatch(route="normalized bar complex", reason="x", n_relations=1),
              ResolutionTerm(degree=0, n_generators=1, collapsed_dim=2),
              RankStep(degree=0, side="cochain", nrows=0, ncols=2, rank=0, field="CC",
                       matrix=[])]
    assert _drop_superseded(normal) == normal
    # a module-style stream (zero Dispatches) is likewise untouched
    module = [ResolutionTerm(degree=0, n_generators=1, collapsed_dim=2),
              StepNote("chose the projective cover")]
    assert _drop_superseded(module) == module


# --------------------------------------------------------------------------- #
# Correction 1: the buffer-elision cap and the drift gate no longer contradict.
# --------------------------------------------------------------------------- #
def _elided_stream(max_events=5):
    """A bar HH^ computation whose FULL dims are [3,2,2,2,2], captured under a buffer cap
    small enough to ELIDE -- so ``derive_dims`` over the survivors is a SHORT prefix."""
    A = truncated_polynomial(3, field=CC)
    tr = Trace(max_events=max_events)
    table = A.hochschild_cohomology(4, trace=tr)
    return A, list(tr), table


def test_elision_present_and_derive_is_truncated():
    """Guard the premise: the buffer elided and the derived dims are a SHORT prefix that
    disagrees with the engine (exactly what would trip the gate)."""
    _A, ev, table = _elided_stream()
    assert any(isinstance(e, StepNote) and BUFFER_FULL_PREFIX in e.text for e in ev)
    assert derive_dims(ev) != table.dims and len(derive_dims(ev)) < len(table.dims)


def test_gate_skips_under_elision_and_ships_correct_result():
    _A, ev, table = _elided_stream()
    out = _authoritative_result(list(ev), table)                 # must NOT raise
    (rd,) = [e for e in out if isinstance(e, ResultDims)]
    assert rd.dims == table.dims == [3, 2, 2, 2, 2] and rd.kind == "HH^"


def test_gate_still_has_teeth_without_the_elision_note():
    """The skip is elision-SPECIFIC: the SAME truncated survivors WITHOUT the reserved
    BUFFER_FULL StepNote still raise, so a genuinely-drifting report is not waved through."""
    _A, ev, table = _elided_stream()
    no_note = [e for e in ev
               if not (isinstance(e, StepNote) and BUFFER_FULL_PREFIX in e.text)]
    with pytest.raises(QuiverlabError, match="drift"):
        _authoritative_result(no_note, table)


def test_elided_report_ships_with_note_and_authoritative_result(tmp_path, monkeypatch):
    A, ev, table = _elided_stream()
    html, obj = _write_html(ev, table, A, "HH^", table.top, monkeypatch, tmp_path)
    # the over-long report SHIPS: the elision note is shown AND the Result is the engine's
    assert "later steps were elided" in html
    for i, d in enumerate(table.dims):
        assert r"HH^{%d} = %d" % (i, d) in html
    (rd,) = [e for e in obj["events"] if e["type"] == "ResultDims"]
    assert rd["dims"] == table.dims
    # under elision the ResultDims must NOT carry the misleading "fast engine" note
    # (per-degree steps DID exist -- they were elided, not absent)
    assert "records no per-degree worked steps" not in (rd["note"] or "")
