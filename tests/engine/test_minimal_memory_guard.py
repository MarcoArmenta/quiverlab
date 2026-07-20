"""Re-port of the deleted bank test_minimal_memory_guard (ledger). open_33_0 =
k<x,y>/(x^3 - y^2, y^3, y*x + x*y), dim 9, local, non-monomial -- now buildable via the
Plan-03 general kQ/I lowering (reduction_algebra is superseded by Quiver.algebra). The
minimal A^e resolution's radK transient grows degree by degree, exercising the
max_transient_bytes guard: any budget yields a PREFIX; a 1-byte budget yields
truncated_at == 0, hh == []."""
from quiverlab import Quiver, GF
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.resolutions_minimal import minimal_homology_dims, minimal_resolution

P = 32003


def _open_33_0():
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x^3 - y^2", "y^3", "y*x + x*y"], field=GF(P))
    assert A.dim == 9
    return to_engine(A)


def test_default_none_matches_existing_build():
    A = _open_33_0()
    base = minimal_homology_dims(A, 5, primes=(P,))
    same = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=None)
    assert base == same


def test_any_transient_budget_yields_a_prefix():
    A = _open_33_0()
    full = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=10**18)[P]
    for budget in (1, 10**5, 10**6, 10**7, 10**9):
        got = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=budget)[P]
        assert got == full[:len(got)], f"budget={budget}: {got} not a prefix of {full}"


def test_tiny_budget_truncates_strictly_earlier():
    A = _open_33_0()
    full = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=10**18)[P]
    tiny = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=1)[P]
    assert len(tiny) < len(full)


def test_truncated_list_length_equals_truncated_at_and_is_a_prefix():
    A = _open_33_0()
    full = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=10**18)[P]
    budget = 10**7
    hh = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=budget)[P]
    _rks, _cols, _eng, trunc = minimal_resolution(A, 5, P, max_transient_bytes=budget)
    assert trunc is not None
    assert len(hh) == trunc
    assert hh == full[:trunc]


def test_one_byte_budget_stops_before_any_differential():
    A = _open_33_0()
    hh = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=1)[P]
    _rks, _cols, _eng, trunc = minimal_resolution(A, 5, P, max_transient_bytes=1)
    assert trunc == 0
    assert hh == []
