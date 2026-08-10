"""Plan 74 Task 1: QuiverAutomorphism + GroupAction with a per-instance
automorphism certificate."""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.families.basic import truncated_polynomial
from quiverlab.families.skew_group import GroupAction, QuiverAutomorphism
from quiverlab.fields import GF, QQ

pytestmark = pytest.mark.oracle_selfcert


def _sigma():
    # x |-> -x on k[x]/(x^2): identity vertex-perm, arrow scalar s(x) = -1
    return QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": -1})


def test_sigma_matrix_and_check():
    A = truncated_polynomial(2, field=QQ)
    dom = A.domain
    sig = _sigma()
    M = sig.matrix(A)
    expected = [[dom.one(), dom.zero()], [dom.zero(), dom.coerce(-1)]]
    assert M == expected
    assert sig.check(A) is True


def test_default_scalar_is_one():
    A = truncated_polynomial(2, field=QQ)
    idn = QuiverAutomorphism({1: 1}, {"x": "x"})       # no scalars -> identity
    dom = A.domain
    assert idn.matrix(A) == [[dom.one(), dom.zero()], [dom.zero(), dom.one()]]
    assert idn.check(A) is True


def test_float_scalar_refused():
    with pytest.raises(QuiverlabError):
        QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": 0.5})


def test_non_automorphism_relation_breaker_raises():
    # k<x,y>/(xx, xy, yx, yyy): x is nilpotent of order 2, y of order 3 (ASYMMETRIC).
    # Swapping x<->y sends x^2 |-> y^2 (nonzero, NOT in I) -> not an automorphism.
    Q = Quiver(vertices=[1], arrows={"x": (1, 1), "y": (1, 1)})
    A = Q.algebra(relations=["x*x", "x*y", "y*x", "y*y*y"], field=QQ)
    swap = QuiverAutomorphism({1: 1}, {"x": "y", "y": "x"})
    with pytest.raises(QuiverlabError):
        swap.check(A)


def test_incompatible_perm_raises():
    # kA2 (1 --a--> 2): swapping vertices without reversing the arrow is incompatible.
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    A = Q.algebra(relations=[], field=QQ)
    bad = QuiverAutomorphism({1: 2, 2: 1}, {"a": "a"})
    with pytest.raises(QuiverlabError):
        bad.matrix(A)


def test_cyclic_group_structure():
    ga = GroupAction.cyclic(_sigma(), 2)
    assert ga.order == 2
    assert ga.conjugacy_classes == [[0], [1]]
    assert ga.centralizer(1) == [0, 1]
    A = truncated_polynomial(2, field=QQ)
    assert ga.check(A) is True


def test_cyclic_order3_over_gf7():
    # sigma(x) = 2x on k[x]/(x^3) over GF(7): 2^3 = 8 = 1 mod 7, order 3.
    A = truncated_polynomial(3, field=GF(7))
    sig = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": 2})
    ga = GroupAction.cyclic(sig, 3)
    assert ga.order == 3
    assert ga.check(A) is True
    # the three element matrices are diag(1,1,1), diag(1,2,4), diag(1,4,2)
    dom = A.domain
    mats = [ga.matrix_of(i, A) for i in range(3)]
    assert mats[1][1][1] == dom.coerce(2) and mats[1][2][2] == dom.coerce(4)
    assert mats[2][1][1] == dom.coerce(4) and mats[2][2][2] == dom.coerce(2)


def test_cyclic_wrong_order_over_field_refused():
    # claiming order 2 for sigma(x)=2x over GF(7) is a lie (true order 3): check
    # must FAIL because the recorded elements are not distinct / not closed.
    A = truncated_polynomial(3, field=GF(7))
    sig = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": 2})
    ga = GroupAction.cyclic(sig, 2)
    with pytest.raises(QuiverlabError):
        ga.check(A)


def test_from_generators_s3_non_abelian():
    # S3 permuting the 3 vertices of k^3 (semisimple): order 6, three conj classes.
    A = Quiver(vertices=[1, 2, 3], arrows={}).algebra(relations=[], field=GF(7))
    transp = QuiverAutomorphism({1: 2, 2: 1, 3: 3}, {})
    cyc = QuiverAutomorphism({1: 2, 2: 3, 3: 1}, {})
    ga = GroupAction.from_generators([transp, cyc], A)
    assert ga.order == 6
    assert ga.check(A) is True
    sizes = sorted(len(c) for c in ga.conjugacy_classes)
    assert sizes == [1, 2, 3]
