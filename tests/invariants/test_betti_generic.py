"""Generic-Domain Betti numbers of the minimal A^e resolution (Plan 19):
exact GF(p) parity with engine/resolutions_minimal (incl. multi-vertex and
straddling Plan-18 records), char-0/GF(4) closed forms, complexity routing."""
import pytest

from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.families import QuantumCI
from quiverlab.families.zoo import build_from_record, load_catalog
from quiverlab.fields import QQ


def _rec(name):
    return next(r for r in load_catalog() if r.get("name") == name)


def test_betti_matches_engine_rks_over_gfp():
    # depth N is per-algebra: the generic chain count is (dim r)^n on a
    # single vertex, so 5-radical algebras get N=2 (still past the degree
    # where the Plan-12 straddling chains appear) and 3-radical ones N=3
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import minimal_resolution
    from quiverlab.invariants.betti import relative_betti_numbers
    zoo = [(truncated_polynomial(3, field=GF(5)), 5),
           (QuantumCI(2, field=GF(5)), 3),
           (linear_path_algebra(2, field=GF(5)), 5),
           (build_from_record(_rec("straddle_xx_yy_xyx"), field=GF(5)), 2),
           (build_from_record(_rec("comm_square"), field=GF(5)), 5),
           (build_from_record(_rec("cn_3_2"), field=GF(3)), 5)]
    for A, N in zoo:
        p = A.domain.p
        rks, _cols, _e, _trunc = minimal_resolution(to_engine(A), N, p)
        want = [max(0, rks[k]) for k in sorted(rks)]
        assert len(want) == N + 2  # engine convention: rks covers degrees 0..N+1
        got = relative_betti_numbers(A, len(want) - 1)
        assert got == want, f"{getattr(A, 'zoo_name', repr(A))}: {got} != {want}"


def test_betti_char0_closed_forms():
    from quiverlab.invariants.betti import relative_betti_numbers
    # monomial: char-independent closed forms
    assert relative_betti_numbers(truncated_polynomial(3, field=QQ), 5) == [1] * 6
    assert relative_betti_numbers(linear_path_algebra(2, field=QQ), 4) == [2, 1, 0, 0, 0]
    # commutative complete intersection k[x,y]/(x^2, y^2): Betti_n = n + 1
    assert relative_betti_numbers(QuantumCI(-1, field=QQ), 4) == [1, 2, 3, 4, 5]
    # gldim-2 multi-vertex: (vertices, arrows, relations, 0, 0, ...)
    assert relative_betti_numbers(
        build_from_record(_rec("comm_square"), field=QQ), 5) == [4, 4, 1, 0, 0, 0]


def test_betti_dd_zero_over_qq():
    # the middle-face differential is a complex: d o d = 0, asserted exactly
    from quiverlab.invariants.betti import _relative_complex
    A = QuantumCI(-1, field=QQ)
    dom = A.domain
    chains, mats = _relative_complex(A, 4, 4_000_000)
    for n in range(3, 5):
        M1, M2 = mats[n], mats[n - 1]          # d_n, d_{n-1}
        for ci in range(len(chains[n])):
            col = [M1[r][ci] for r in range(len(chains[n - 1]))]
            acc = [dom.zero()] * len(chains[n - 2])
            for j, cv in enumerate(col):
                if dom.is_zero(cv):
                    continue
                for r in range(len(chains[n - 2])):
                    acc[r] = dom.add(acc[r], dom.mul(cv, M2[r][j]))
            assert all(dom.is_zero(x) for x in acc), f"d o d != 0 at {n}"


def test_complexity_field_generality():
    assert truncated_polynomial(2, field=CC).complexity(4) == 1
    assert truncated_polynomial(2, field=GF(4)).complexity(4) == 1
    assert QuantumCI(-1, field=QQ).complexity(3) == 2
    assert linear_path_algebra(2, field=QQ).complexity(4) == 0
