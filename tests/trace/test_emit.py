"""Engines emit typed events; claims equal computed values (spec §3.8, §8)."""
import quiverlab
from quiverlab import truncated_polynomial, CC, GF
from quiverlab.trace.events import Dispatch, ResolutionTerm, RankStep
from quiverlab.trace.recorder import Trace
from quiverlab.hochschild.bar import coboundary_matrix, hochschild_cohomology_dims
from quiverlab.fields.linalg import rank


def test_bar_dims_emit_terms_and_ranks():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = hochschild_cohomology_dims(A, 2, trace=tr)
    terms = [e for e in tr if isinstance(e, ResolutionTerm)]
    ranks = [e for e in tr if isinstance(e, RankStep)]
    assert [t.degree for t in terms] == [0, 1, 2]
    assert [t.collapsed_dim for t in terms] == [2, 2, 2]
    assert [r.rank for r in ranks] == [0, 1, 0]           # machine-verified (F5)
    assert ranks[1].matrix == [["0", "0"], ["2", "0"]]     # d^1 over CC
    # claims == computed: recorded rank equals an independent recomputation
    B = A.unit_adapted()
    D, nc, nr = coboundary_matrix(B, 1, 4_000_000)
    assert ranks[1].rank == rank(D, B.domain)
    assert table.dims == [2, 1, 1]


def test_bar_homology_dims_emit_terms_and_ranks():
    # mirror of the cohomology emit test, pinning the homology recording contract
    from quiverlab.hochschild.bar import boundary_matrix, hochschild_homology_dims
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = hochschild_homology_dims(A, 2, trace=tr)
    terms = [e for e in tr if isinstance(e, ResolutionTerm)]
    ranks = [e for e in tr if isinstance(e, RankStep)]
    assert [t.degree for t in terms] == [0, 1, 2]
    assert ranks and all(r.side == "chain" for r in ranks)   # homology side (not "cochain")
    # claims == computed: each recorded rank + shape equals an independent recomputation
    B = A.unit_adapted()
    for r in ranks:
        D, ncols, nrows = boundary_matrix(B, r.degree, 4_000_000)
        assert r.rank == rank(D, B.domain)
        assert (r.nrows, r.ncols) == (nrows, ncols)
    # recording is a pure side-channel: dims identical to a no-trace call
    assert table.dims == hochschild_homology_dims(A, 2).dims


def test_algebra_records_engine_dispatch_via_trace_param():
    A = truncated_polynomial(2, field=CC)
    tr = []
    A.hochschild_cohomology(2, trace=tr)
    disp = [e for e in tr if isinstance(e, Dispatch)]
    assert disp and disp[0].route == "normalized bar complex"
    assert disp[0].n_relations == 1                        # x^2


def test_prime_field_uses_fast_engine_dispatch():
    A = truncated_polynomial(2, field=GF(2))
    tr = []
    A.hochschild_cohomology(2, trace=tr)
    routes = [e.route for e in tr if isinstance(e, Dispatch)]
    assert routes and "fast" in routes[0].lower()


def test_explicit_trace_does_not_flip_global_verbose(tmp_path, monkeypatch):
    # passing trace=[] is programmatic: it must NOT write a file even if verbose is on
    monkeypatch.chdir(tmp_path)
    quiverlab.verbose = True
    try:
        A = truncated_polynomial(2, field=CC)
        A.hochschild_cohomology(2, trace=[])
    finally:
        quiverlab.verbose = False
    assert not (tmp_path / "quiverlab_traces").exists()
