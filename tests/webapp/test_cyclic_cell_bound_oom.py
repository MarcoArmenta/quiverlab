"""Regression: GF(p) cyclic homology must honour ``max_cells`` (the
hotfix-cyclic-cell-guard OOM bug), the sibling of ``test_connes_cell_bound_oom``.

``Algebra.cyclic_homology(top, max_cells=4_000_000)`` promised a cell contract, but the
GF(p) engine branch (``core/algebra.py`` -> ``engine/cyclic.py::cyclic_homology_dims``)
never received or applied it -- honoured on the generic-Domain branch
(``hochschild/cyclic.py``) and by the connes_b guard. On the curated dim-36 multi-vertex
grid3x3 algebra a ``cyclic_homology:0..2`` request SIGKILLed the memory-capped worker
(exit 137): the bar chain bases grow like 36*35^n, so b_2 alone is 1260 x 44100 = 56M
cells (>> the 4M cap) and even MATERIALIZING the degree-3 basis (1.5M tuples) exhausts
memory. The wave-3b ``with_reps`` path rides the same engine call, so the guard covers it.

The fix: (1) the GF(p) engine honours ``max_cells`` and refuses LOUDLY before allocating
(over GF(p) cyclic homology has no Chouhy-Solotar route, so -- like an explicit bar engine
-- it raises ``DepthLimitError`` rather than silently OOMing), computing the matrix sizes
by LENGTH ARITHMETIC (m*(m-1)^n), never by enumerating a basis to count it; (2) the
server runner already degrades per-item on ``DepthLimitError`` into an honest error block
(the bcd90f2 pattern), so the whole request COMPLETES -- every other invariant + the
worked-steps report still ship. Unmarked (contract/infrastructure) and cheap: the guard
trips before any matrix is built, so this never allocates gigabytes in CI.
"""
import json
import pathlib

import pytest

import quiverlab as ql
from quiverlab.errors import DepthLimitError
from quiverlab.fields import GF
from quiverlab.hochschild.cyclic import cyclic_homology_dims as generic_hc
from webapp.server.runner import build_algebra, run_spec
from webapp.server.schema import ComputeRequest

_GRID3X3 = json.loads(
    (pathlib.Path(__file__).resolve().parents[2]
     / "webapp" / "precomputed" / "examples" / "grid3x3" / "request.json")
    .read_text(encoding="utf-8"))


def _req(compute, pdf=False):
    return ComputeRequest.model_validate({
        "schema": 2,
        "algebra": _GRID3X3["algebra"],           # the dim-36 grid3x3 quiver algebra
        "compute": compute,
        "artifacts": {"pdf": pdf, "tikz": False}})


def test_cyclic_over_cap_raises_before_allocating():
    """The library guard trips BEFORE any dense matrix (or its exponential basis) is
    built -- cheap + loud. HC:0..0 is the only feasible degree on this dim-36 algebra
    (b_1 = 45360 cells); HC:0..1 already needs b_2 = 1260 x 44100 = 56M cells > 4M, and
    HC:0..2 more. A feasible top's dims are the generic route's (cross-engine anchor)."""
    A = build_algebra(_req(["dimension"]).algebra)
    assert A.dim == 36
    # feasible top: computes, and agrees with the independent generic (b, B) route.
    assert A.cyclic_homology(0).dims == generic_hc(A, 0).dims == [9]
    for top in (1, 2):                             # over-cap: refuses loudly
        with pytest.raises(DepthLimitError) as ei:
            A.cyclic_homology(top)
        assert "max_cells" in str(ei.value)
    # a raised max_cells lifts the bound and still computes (feature preserved).
    assert A.cyclic_homology(1, max_cells=200_000_000).dims == [9, 0]


def test_feasible_top_matches_generic_cross_engine():
    """A feasible computation is UNCHANGED by the guard: the GF(p) engine and the
    generic-Domain (b, B) route agree degreewise -- the correctness anchor. b,B here is
    tiny (dim 3), far under the cap, so nothing trips."""
    A = ql.truncated_polynomial(3, field=GF(101))
    top = 4
    assert A.cyclic_homology(top).dims == generic_hc(A, top).dims == [3, 0, 3, 0, 3]


def test_run_spec_completes_when_cyclic_is_over_cap(tmp_path):
    """The whole request completes with pdf=true: cyclic_homology becomes an honest
    per-item error block, every other invariant + the worked-steps report still ship."""
    res = run_spec(_req(["cyclic_homology:0..2", "dimension"], pdf=True), tmp_path)
    r = res["results"]
    # (a) completion -- all kinds present, no exception, no OOM.
    assert set(r) == {"cyclic_homology", "dimension"}
    # cyclic_homology is the honest refusal, not a computed table.
    assert r["cyclic_homology"]["error"]["type"] == "DepthLimitError"
    assert "max_cells" in r["cyclic_homology"]["error"]["message"]
    # (b) the other invariant stands, exact.
    assert r["dimension"]["value"] == 36
    # (c) the report bundle is written and states the honest refusal.
    html = (tmp_path / "trace_steps.html").read_text(encoding="utf-8")
    assert "DepthLimitError" in html
    assert (tmp_path / "trace.json").exists()


def test_over_cap_cyclic_never_perturbs_the_other_results(tmp_path):
    """The other results are byte-identical with and without the over-cap cyclic
    computation riding along -- the refusal is isolated to its own block."""
    a = run_spec(_req(["dimension"]), tmp_path / "a")
    b = run_spec(_req(["dimension", "cyclic_homology:0..2"]), tmp_path / "b")
    assert a["results"]["dimension"] == b["results"]["dimension"]
    assert b["results"]["cyclic_homology"]["error"]["type"] == "DepthLimitError"
