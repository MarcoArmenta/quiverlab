"""Reduced walks, string signs, string enumeration, band detection (Plan 46 / C5).

Self-cert: every enumerated walk is valid and canonicalised; the sign functions
satisfy the Butler-Ringel conditions. Literature: kA_n has n(n+1)/2 strings (the
interval modules) and NO bands (rep-finite); the Kronecker quiver HAS a band.

PLAN CORRECTION (documented): the plan's original band example, the 2-cycle
kQ/(ab,ba), is a self-injective Nakayama algebra -- in/out-degree <=1 at every
vertex makes it uniserial and representation-FINITE, so it has NO band. Any algebra
on the 2-cycle quiver is Nakayama, hence never has a band. The correct minimal
gentle band algebra is the Kronecker quiver (two parallel arrows), whose band is
a * b^{-1}. We pin BOTH: the Kronecker finds its band, and the 2-cycle honestly
reports rep-finite with no band."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.fields import QQ
from quiverlab.strings.walks import (enumerate_strings, find_bands, invert,
                                     is_valid_walk, string_signs)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _kA3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=QQ)


def _gentle_a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)


def _two_cycle():
    # kQ/(ab, ba): a: 1->2, b: 2->1, self-injective gentle Nakayama, dim 4. NO band.
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


def _kronecker(field=QQ):
    # two parallel arrows a, b: 1 -> 2, no relations; gentle, rep-infinite; band a b^{-1}.
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(
        relations=[], field=field)


@selfcert
def test_signs_satisfy_string_conditions():
    A = _gentle_a3()
    sigma, epsilon = string_signs(A)
    assert set(sigma) == set(A.quiver.arrows) == set(epsilon)
    assert all(v in (1, -1) for v in sigma.values())
    assert all(v in (1, -1) for v in epsilon.values())
    # (S3) applies exactly to a nonzero composable pair. In gentle kA3/(ab), a*b IS
    # in I so (S3) is vacuous; test it on kA3 with NO relations where a*b NOT in I:
    B = _kA3()
    sb, eb = string_signs(B)
    assert sb["b"] == -eb["a"]                          # sigma(gamma) = -epsilon(beta)


@selfcert
def test_enumerated_walks_are_valid_and_canonical():
    A = _gentle_a3()
    cen = enumerate_strings(A, max_length=6)
    for w in cen.walks:
        assert is_valid_walk(A, w)
        assert tuple(w) <= tuple(invert(e) for e in reversed(w))  # canonical rep


@lit
@pytest.mark.parametrize("n, count", [(2, 3), (3, 6), (4, 10)])
def test_linear_a_n_string_count_is_interval_count(n, count):
    # kA_n indecomposables = the n(n+1)/2 interval modules, ALL string modules; no
    # bands (rep-finite). Butler-Ringel; the interval classification.
    arrows = {chr(ord("a") + i): (i + 1, i + 2) for i in range(n - 1)}
    A = Quiver(list(range(1, n + 1)), arrows).algebra(relations=[], field=QQ)
    cen = enumerate_strings(A, max_length=2 * n)
    assert cen.status == "complete"
    assert cen.count == count          # n(n+1)/2
    assert find_bands(A, max_length=2 * n) == []


@lit
def test_kronecker_has_a_band_and_is_not_complete():
    A = _kronecker()
    bands = find_bands(A, max_length=6)
    assert bands != []                                   # rep-infinite
    cen = enumerate_strings(A, max_length=6)
    assert cen.status == "budget" and cen.has_bands      # never claim "complete"


@lit
def test_two_cycle_is_repfinite_no_bands():
    # PLAN CORRECTION: the 2-cycle is self-injective Nakayama, rep-finite, no band.
    A = _two_cycle()
    assert find_bands(A, max_length=6) == []
    cen = enumerate_strings(A, max_length=6)
    assert cen.has_bands is False
    # 4 indecomposables: S_1, S_2, and the two length-2 projective-injectives.
    assert cen.status == "complete"
    assert cen.count == 4


@selfcert
def test_non_string_algebra_refused():
    from quiverlab.errors import QuiverlabError
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3), "c": (1, 4)}).algebra(
        relations=[], field=QQ)                          # 3 arrows out of 1: not SB
    with pytest.raises(QuiverlabError):
        string_signs(A)
