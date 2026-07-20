"""Trace recorder: buffering caps, matrix elision, verbose resolution (spec §3.8
'performance guard: verbose must not blow up long computations')."""
from quiverlab.fields import CC
from quiverlab.trace.recorder import (
    Trace, rankstep, resolve_verbose, MAX_EVENTS, MATRIX_ELISION_CELLS,
)
from quiverlab.trace.events import RankStep, Dispatch


def test_trace_is_list_compatible_for_groebner():
    tr = Trace()
    tr.append(Dispatch(route="monomial", reason="x", n_relations=1))  # .append == .record
    assert len(tr) == 1 and list(tr)[0].route == "monomial"


def test_buffer_cap_drops_and_counts_overflow():
    tr = Trace(max_events=3)
    for i in range(10):
        tr.append(Dispatch(route=str(i), reason="", n_relations=0))
    assert len(tr) == 3
    assert tr.elided_events == 7
    assert any("elided" in n.lower() for n in tr.elision_notes)


def test_rankstep_keeps_small_matrix():
    # non-square (2 rows x 3 cols) pins the (nrows, ncols) param order and row-major
    # D[row][col] rendering -- a swap of nrows/ncols would IndexError or transpose here.
    dom = CC.make_domain([CC.parse_entry(0), CC.parse_entry(1), CC.parse_entry(2)])
    D = [[dom.zero(), dom.zero(), dom.zero()], [dom.coerce(2), dom.zero(), dom.zero()]]
    r = rankstep(1, "cochain", D, 2, 3, 1, dom)
    assert isinstance(r, RankStep)
    assert r.elided is False
    assert r.nrows == 2 and r.ncols == 3
    assert r.matrix == [["0", "0", "0"], ["2", "0", "0"]]
    assert r.field == dom.name


def test_rankstep_elides_large_matrix():
    dom = CC.make_domain([CC.parse_entry(0), CC.parse_entry(1)])
    n = 30  # 30*30 = 900 > MATRIX_ELISION_CELLS (400)
    D = [[dom.zero()] * n for _ in range(n)]
    r = rankstep(4, "cochain", D, n, n, 0, dom)
    assert r.elided is True and r.matrix is None
    assert str(MATRIX_ELISION_CELLS) in r.note or "elided" in r.note.lower()
    assert r.nrows == n and r.ncols == n and r.rank == 0  # shape + rank always kept


def test_resolve_verbose_per_call_overrides_global():
    assert resolve_verbose(None, True) is True
    assert resolve_verbose(None, False) is False
    assert resolve_verbose(False, True) is False
    assert resolve_verbose(True, False) is True
