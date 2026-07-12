from fractions import Fraction as F

import pytest
from quiverlab.core.algebra import Algebra
from quiverlab.errors import ExactnessError, FieldError, QuiverlabError
from quiverlab.fields import CC, GF


def _dual_numbers(field=CC):
    # basis (1, x), x^2 = 0
    T = [
        [[1, 0], [0, 1]],
        [[0, 1], [0, 0]],
    ]
    return Algebra.from_structure_constants(T, unit=[1, 0], field=field)


def test_dual_numbers_multiply():
    A = _dual_numbers()
    dom = A.domain
    one = [dom.coerce(1), dom.coerce(0)]
    x = [dom.coerce(0), dom.coerce(1)]
    assert A.multiply(x, x) == [dom.zero(), dom.zero()]
    assert A.multiply(one, x) == x
    assert A.is_unit_adapted


def test_associativity_check_catches_garbage():
    T = [
        [[1, 0], [0, 1]],
        [[0, 1], [1, 0]],  # x*x = 1 ... with x*1 = x this IS associative (k[x]/(x^2-1)); tweak:
    ]
    # break unit instead: claim unit = [0, 1]
    with pytest.raises(QuiverlabError):
        Algebra.from_structure_constants(T, unit=[0, 1], field=CC)


def test_nonassociative_rejected():
    # b1*b1 = b1 with b1 the unit, b2*b2 = b1, b2*b1 = b2, but b1*b2 = 0: not unital/associative
    T = [
        [[1, 0], [0, 0]],
        [[0, 1], [1, 0]],
    ]
    with pytest.raises(QuiverlabError):
        Algebra.from_structure_constants(T, unit=[1, 0], field=CC)


def test_floats_rejected_in_T():
    T = [
        [[1, 0], [0, 1]],
        [[0, 1], [0.0, 0]],
    ]
    with pytest.raises(ExactnessError):
        Algebra.from_structure_constants(T, unit=[1, 0], field=CC)


def test_gf_algebra():
    A = _dual_numbers(field=GF(2))
    assert A.domain.characteristic == 2
    assert A.dim == 2


def test_unit_adapted_transform():
    # k x k with basis (e1, e2): unit = e1 + e2 is NOT a basis vector
    T = [
        [[1, 0], [0, 0]],
        [[0, 0], [0, 1]],
    ]
    A = Algebra.from_structure_constants(T, unit=[1, 1], field=CC)
    assert not A.is_unit_adapted
    B = A.unit_adapted()
    assert B.is_unit_adapted
    dom = B.domain
    # in the new basis, b0 is the unit
    e0 = [dom.one()] + [dom.zero()] * (B.dim - 1)
    v = [dom.coerce(3), dom.coerce(-2)]
    assert B.multiply(e0, v) == v and B.multiply(v, e0) == v


def test_singular_change_of_basis_rejected():
    A = _dual_numbers()
    dom = A.domain
    one, zero = dom.one(), dom.zero()
    P = [[one, one], [zero, zero]]  # rank 1: columns do not form a basis
    with pytest.raises(QuiverlabError) as ei:
        A.change_of_basis(P)
    assert "singular" in str(ei.value)


def test_genuinely_nonassociative_with_good_unit_rejected():
    T = [
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        [[0, 1, 0], [0, 0, 1], [0, 0, 0]],
        [[0, 0, 1], [0, 1, 0], [0, 0, 0]],
    ]
    with pytest.raises(QuiverlabError) as ei:
        Algebra.from_structure_constants(T, unit=[1, 0, 0], field=CC)
    assert "associative" in str(ei.value)
