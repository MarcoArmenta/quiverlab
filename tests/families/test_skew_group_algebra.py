"""Plan 74 Task 2: the skew_group_algebra A rtimes G constructor -- dim law,
trivial-G byte-identity, presented_form orbit-quiver view."""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.families.basic import truncated_polynomial
from quiverlab.families.skew_group import (
    GroupAction, QuiverAutomorphism, is_free_action, skew_group_algebra)
from quiverlab.fields import QQ

pytestmark = [pytest.mark.oracle_selfcert, pytest.mark.oracle_crossengine]


def _dual_numbers_z2():
    A = truncated_polynomial(2, field=QQ)
    sig = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": -1})
    return A, GroupAction.cyclic(sig, 2)


def test_dim_law_and_multiplication_table():
    A, ga = _dual_numbers_z2()
    S = skew_group_algebra(A, ga)
    dom = A.domain
    assert S.dim == 2 * A.dim == 4
    # basis index n(a, g) = a*order + g; labels "<a>.g<g>"
    assert S.basis_labels == ["e_1.g0", "e_1.g1", "x.g0", "x.g1"]
    # (x.1)(x.1) = x*1(x).1 = x^2.1 = 0
    x1 = S._basis_vec(2)
    assert S.multiply(x1, x1) == [dom.zero()] * 4
    # (e.sigma)(x.1) = e*sigma(x).sigma = -x.sigma  -> coeff -1 at index n(x,sigma)=3
    esig = S._basis_vec(1)
    got = S.multiply(esig, x1)
    assert got[3] == dom.coerce(-1)
    assert all(dom.is_zero(got[i]) for i in range(4) if i != 3)
    # (e.sigma)^2 = e.1  -> index n(e,1)=0
    got2 = S.multiply(esig, esig)
    assert got2[0] == dom.one()
    assert all(dom.is_zero(got2[i]) for i in range(1, 4))


@pytest.mark.oracle_crossengine
def test_trivial_group_byte_identity():
    A = truncated_polynomial(2, field=QQ)
    idn = QuiverAutomorphism({1: 1}, {"x": "x"})
    ga = GroupAction.cyclic(idn, 1)
    S = skew_group_algebra(A, ga)
    assert S.T == A.T                      # structure constants byte-identical
    assert list(S.unit) == list(A.unit)
    assert S.hochschild_cohomology(3).dims == A.hochschild_cohomology(3).dims
    assert S.hochschild_homology(3).dims == A.hochschild_homology(3).dims


def test_presented_form_recovers_orbit_nakayama():
    A, ga = _dual_numbers_z2()
    S = skew_group_algebra(A, ga)
    pf = S.presented_form()
    assert pf.dim == 4
    assert sorted(pf.quiver.vertices) == [1, 2]
    assert set(pf.quiver.arrows.values()) == {(1, 2), (2, 1)}
    assert sorted(str(r) for r in pf.relations) == ["a1*a2", "a2*a1"]


def test_flagship_is_not_free_action():
    A, ga = _dual_numbers_z2()
    assert is_free_action(A, ga) is False   # sigma fixes the single vertex


def test_delegate():
    A, ga = _dual_numbers_z2()
    S = A.skew_group(ga)
    assert S.dim == 4
