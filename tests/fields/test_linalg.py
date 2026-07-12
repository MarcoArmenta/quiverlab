from fractions import Fraction as F

from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace, rank, rref, solve


def _m(rows):
    return [[QQ.coerce(x) for x in r] for r in rows]


def test_rank_and_rref():
    A = _m([[1, 2, 3], [2, 4, 6], [1, 0, 1]])
    assert rank(A, QQ) == 2
    R, piv = rref(A, QQ)
    assert piv == [0, 1]
    assert R[0][0] == F(1) and R[1][1] == F(1)


def test_nullspace_is_kernel():
    A = _m([[1, 2, 3], [2, 4, 6], [1, 0, 1]])
    N = nullspace(A, QQ)
    assert len(N) == 1
    v = N[0]
    for row in A:
        s = QQ.zero()
        for a, x in zip(row, v):
            s = QQ.add(s, QQ.mul(a, x))
        assert QQ.is_zero(s)


def test_solve():
    A = _m([[1, 1], [1, -1]])
    b = [QQ.coerce(3), QQ.coerce(1)]
    x = solve(A, b, QQ)
    assert x == [F(2), F(1)]
    A2 = _m([[1, 1], [1, 1]])
    assert solve(A2, [QQ.coerce(0), QQ.coerce(1)], QQ) is None


def test_empty_matrix():
    assert rank([], QQ) == 0
    assert nullspace([], QQ) == []
