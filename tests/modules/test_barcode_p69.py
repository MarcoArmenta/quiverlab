"""Barcodes = interval decompositions of A_n / zigzag modules (Plan 69 / R33).
Literature/self-cert: the A_5 filtration barcode {[1,5],[2,2],[4,4]} (H_0 of a filtration
with dims (1,2,1,2,1); HAND-COMPUTED + LIVE-VERIFIED, identical over QQ and GF(2)).
Cross-engine: decompose == rank-formula (forward); GF(2) == QQ (field-robust bricks);
barcode == decompose summands."""
import pytest
from quiverlab import Quiver, GF
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.barcode import barcode

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert

FWD = {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 5)}


def _H0(field):
    # H_0 of a filtration: dims (1,2,1,2,1); expected barcode {[1,5],[2,2],[4,4]}.
    # blocks V1={0},V2={1,2},V3={3},V4={4,5},V5={6}; full 7x7 arrow matrices.
    A = Quiver([1, 2, 3, 4, 5], FWD).algebra(relations=[], field=field)
    o, z = field.one(), field.zero()

    def M(pairs):
        m = [[z] * 7 for _ in range(7)]
        for (i, j) in pairs:
            m[i][j] = o
        return m
    return A.module({1: 1, 2: 2, 3: 1, 4: 2, 5: 1},
                    {"a": M([(1, 0)]), "b": M([(3, 1), (3, 2)]),
                     "c": M([(4, 3)]), "d": M([(6, 4), (6, 5)])}, name="H0")


@lit
@pytest.mark.parametrize("field", [QQ, GF(2)])
def test_a5_filtration_barcode(field):
    bc = barcode(_H0(field))
    got = sorted((b.birth, b.death, b.multiplicity) for b in bc.bars)
    assert got == [(1, 5, 1), (2, 2, 1), (4, 4, 1)]     # LIVE-VERIFIED
    assert bc.field_robust is True
    assert any(b.death == 5 and b.essential for b in bc.bars)   # [1,5] essential


@xeng
def test_gf2_equals_qq_field_robust():
    q = sorted((b.birth, b.death, b.multiplicity) for b in barcode(_H0(QQ)).bars)
    g = sorted((b.birth, b.death, b.multiplicity) for b in barcode(_H0(GF(2))).bars)
    assert q == g


@xeng
def test_barcode_equals_rank_formula_forward():
    from quiverlab.modules.barcode import _barcode_by_ranks
    M = _H0(QQ)
    a = sorted((b.birth, b.death, b.multiplicity) for b in barcode(M).bars)
    b = sorted(_barcode_by_ranks(M))                    # independent rank route
    assert a == b


@selfcert
def test_interval_sum_identity_and_support_contiguous():
    M = _H0(QQ)
    bc = barcode(M)
    acc = {v: 0 for v in [1, 2, 3, 4, 5]}
    for bar in bc.bars:
        assert list(range(bar.birth, bar.death + 1)) == sorted(bar.dimvec)  # contiguous
        for v, d in bar.dimvec.items():
            acc[v] += bar.multiplicity * d
    assert acc == M.dimension_vector()                  # Sigma m_i dimvec == M


@selfcert
def test_zigzag_barcode(field=QQ):
    # 1->2<-3->4<-5 : intervals are contiguous in the line order 1-2-3-4-5.
    Z = Quiver([1, 2, 3, 4, 5], {"a": (1, 2), "b": (3, 2), "c": (3, 4), "d": (5, 4)})
    A = Z.algebra(relations=[], field=field)
    from quiverlab.modules.morphism import direct_sum
    M = direct_sum(A.projective(1), A.simple(4), A.injective(2))[0]
    bc = barcode(M)
    assert bc.kind == "zigzag" and bc.field_robust is True
    # M2: the essential flag is forward-only -- every zigzag bar has essential=False,
    # even one whose death reaches the terminal line index (no monotone "top").
    assert all(b.essential is False for b in bc.bars)


@selfcert
def test_a1_single_vertex_handled():
    # M1: a single-vertex A_1 = k has ZERO degree-1 endpoints, so the general line
    # classifier can't see it -- barcode HANDLES it as a special case (one bar [1,1]
    # with multiplicity = dim M), it does NOT refuse.
    A = Quiver([1], {}).algebra(relations=[], field=QQ)
    M = A.module({1: 3}, {}, name="V3")          # a 3-dim vector space at the one vertex
    bc = barcode(M)
    assert bc.kind == "persistence" and bc.n == 1
    got = sorted((b.birth, b.death, b.multiplicity, b.essential) for b in bc.bars)
    assert got == [(1, 1, 3, True)]              # one bar [1,1], mult 3, essential (death==n)


@selfcert
def test_loud_refusals():
    # non-A_n quiver (a vertex of degree 3): not a persistence line.
    kD4 = Quiver([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)}).algebra(
        relations=[], field=QQ)
    with pytest.raises(QuiverlabError):
        barcode(kD4.simple(0))                          # not A_n / not CL
