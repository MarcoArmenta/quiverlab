"""Plan 69 / R33 barcode battery -- a zoo of A_n / zigzag interval sums (n = 2..6) over
QQ AND GF(2), the field-robustness (GF(2) == QQ) and the two-route agreement
(decompose == rank-formula, forward), plus the essential-bar flag correctness. No AR knit
here (the CL diagram tests live in test_barcode_cl_p69.py); this battery is interval-only
and fast."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.fields import QQ
from quiverlab.modules.barcode import barcode, _barcode_by_ranks
from quiverlab.modules.morphism import direct_sum

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert

_FIELDS = [QQ, GF(2)]


def _forward(n, field):
    arrows = {chr(ord("a") + i - 1): (i, i + 1) for i in range(1, n)}
    return Quiver(list(range(1, n + 1)), arrows).algebra(relations=[], field=field)


def _interval_sum(A):
    # a spread of intervals: projectives + simples + injectives (all interval modules).
    verts = A.quiver.vertices
    pieces = [A.projective(verts[0]), A.simple(verts[len(verts) // 2]),
              A.injective(verts[-1])]
    return direct_sum(*pieces)[0]


@xeng
@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
def test_forward_zoo_gf2_equals_qq_and_rank_formula(n):
    key = lambda bc: sorted((b.birth, b.death, b.multiplicity) for b in bc.bars)
    q = key(barcode(_interval_sum(_forward(n, QQ))))
    g = key(barcode(_interval_sum(_forward(n, GF(2)))))
    assert q == g                                      # GF(2) == QQ (bricks)
    r = sorted(_barcode_by_ranks(_interval_sum(_forward(n, QQ))))
    assert q == r                                      # decompose == rank formula (forward)


@selfcert
@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
@pytest.mark.parametrize("field", _FIELDS)
def test_forward_interval_sum_identity_and_essential(n, field):
    A = _forward(n, field)
    M = _interval_sum(A)
    bc = barcode(M)
    assert bc.kind == "persistence" and bc.field_robust is True
    acc = {v: 0 for v in A.quiver.vertices}
    for bar in bc.bars:
        assert list(range(bar.birth, bar.death + 1)) == sorted(bar.dimvec)  # contiguous
        # essential is EXACTLY the death==n forward flag.
        assert bar.essential is (bar.death == n)
        for v, d in bar.dimvec.items():
            acc[v] += bar.multiplicity * d
    assert acc == M.dimension_vector()                 # Sigma m_i dimvec == M


@selfcert
@pytest.mark.parametrize("field", _FIELDS)
def test_zigzag_essential_always_false(field):
    # 1->2<-3->4<-5<-6 (alternating): every bar essential=False (no monotone top, M2).
    Z = Quiver([1, 2, 3, 4, 5, 6],
               {"a": (1, 2), "b": (3, 2), "c": (3, 4), "d": (5, 4), "e": (6, 5)})
    A = Z.algebra(relations=[], field=field)
    M = direct_sum(A.projective(1), A.simple(4), A.injective(2))[0]
    bc = barcode(M)
    assert bc.kind == "zigzag"
    assert all(b.essential is False for b in bc.bars)
    for bar in bc.bars:                                # supports still contiguous in line order
        assert list(range(bar.birth, bar.death + 1)) == sorted(bar.dimvec)
