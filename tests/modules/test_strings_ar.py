"""String-module tau by hooks/cohooks (Plan 46 / C5), ARBITRATED against the trusted
Plan-41/23 engine tau: on every walk, string_module(string_tau(w)) ~
string_module(w).tau(). Cross-engine: the string census count equals knit_ar_quiver's
vertex count on a rep-finite string algebra."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.ar import knit_ar_quiver
from quiverlab.modules.hom import is_isomorphic
from quiverlab.strings.ar_strings import string_tau, string_tau_minus
from quiverlab.strings.modules import string_module
from quiverlab.strings.walks import enumerate_strings

xeng = pytest.mark.oracle_crossengine


def _gentle_a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)


def _kA3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=QQ)


@xeng
@pytest.mark.parametrize("factory", [_kA3, _gentle_a3])
def test_string_tau_matches_engine_tau(factory):
    A = factory()
    cen = enumerate_strings(A, max_length=6)
    for w in cen.walks:
        M = string_module(A, w)
        if M.tau().dim == 0:                          # projective: no non-proj tau
            assert string_tau(A, w) is None
            continue
        wt = string_tau(A, w)
        assert is_isomorphic(string_module(A, wt), M.tau())


@xeng
@pytest.mark.parametrize("factory", [_kA3, _gentle_a3])
def test_string_tau_minus_matches_engine(factory):
    A = factory()
    for w in enumerate_strings(A, max_length=6).walks:
        M = string_module(A, w)
        if M.tau_minus().dim == 0:                    # injective: no non-inj tau^-
            assert string_tau_minus(A, w) is None
            continue
        wm = string_tau_minus(A, w)
        assert is_isomorphic(string_module(A, wm), M.tau_minus())


@xeng
def test_tau_and_tau_minus_are_inverse_on_non_projective_injective():
    A = _gentle_a3()
    for w in enumerate_strings(A, max_length=6).walks:
        M = string_module(A, w)
        if M.tau().dim == 0 or M.tau_minus().dim == 0:
            continue
        back = string_tau_minus(A, string_tau(A, w))
        assert is_isomorphic(string_module(A, back), M)


@xeng
@pytest.mark.parametrize("n, count", [(2, 3), (3, 6), (4, 10)])
def test_census_count_equals_ar_quiver_vertex_count(n, count):
    arrows = {chr(ord("a") + i): (i + 1, i + 2) for i in range(n - 1)}
    A = Quiver(list(range(1, n + 1)), arrows).algebra(relations=[], field=QQ)
    cen = enumerate_strings(A, max_length=2 * n)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    assert cen.count == len(ar.vertices) == count
