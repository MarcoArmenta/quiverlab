"""Quasi-heredity test + Delta-multiplicities + BGG reciprocity (Plan 47). Literature: kA_n
is quasi-hereditary for BOTH the natural and opposite orders (Dlab-Ringel: hereditary => qh
for every order); k[x]/(x^2) is NOT (infinite gl.dim / End Delta != k). BGG:
(P(i):Delta(j))=[Nabla(j):S(i)]. Char-clean (a GF(2) cell)."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.quasihereditary import (costandard_modules, delta_multiplicities,
                                               is_quasi_hereditary, standard_modules)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _a3(field=QQ):
    return linear_path_algebra(3, field=field)


def _dual_numbers(field=QQ):
    # k[x]/(x^2): one vertex, one loop x, x^2 = 0. Local, self-injective, gl.dim = infinity.
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=field)


@lit
@pytest.mark.parametrize("order", [None, [3, 2, 1]])
def test_ka3_is_quasi_hereditary_for_two_orders(order):
    A = _a3()
    rep = is_quasi_hereditary(A, order)
    assert rep.is_quasi_hereditary is True and bool(rep) is True
    assert rep.gl_dim.exact and int(rep.gl_dim) == 1          # hereditary
    assert all(rep.per_index[v]["brick"] for v in (1, 2, 3))
    assert all(rep.per_index[v]["delta_filters_P"] for v in (1, 2, 3))


@lit
def test_dual_numbers_not_quasi_hereditary():
    A = _dual_numbers()
    rep = is_quasi_hereditary(A)
    assert rep.is_quasi_hereditary is False and bool(rep) is False
    assert rep.note and ("gl.dim" in rep.note or "End" in rep.note or "1" in rep.note)


@selfcert
def test_qh_char_clean_gf2():
    A = _a3(GF(2))
    assert is_quasi_hereditary(A).is_quasi_hereditary is True   # no char caveat on the test


@selfcert
def test_delta_filtration_of_P_is_certified():
    A = _a3()
    D = standard_modules(A)                        # natural order: Delta(i)=S_i
    for v in (1, 2, 3):
        mult, ok = delta_multiplicities(A.projective(v), D)
        assert ok is True
        assert mult == {w: (1 if w >= v else 0) for w in (1, 2, 3)}


@selfcert
def test_non_filtered_module_is_uncertified_loudly():
    # a module whose top is not a Delta-top cannot be Delta-peeled: certified=False.
    A = _a3()
    D = {1: A.simple(1)}                            # deliberately incomplete Delta set
    mult, ok = delta_multiplicities(A.projective(1), D)   # P(1) needs Delta(2),Delta(3) too
    assert ok is False


@xeng
def test_bgg_reciprocity():
    A = _a3()
    for order in (None, [3, 2, 1]):
        D = standard_modules(A, order)
        N = costandard_modules(A, order)
        for i in (1, 2, 3):
            mult, ok = delta_multiplicities(A.projective(i), D, order)
            assert ok
            for j in (1, 2, 3):
                # [Nabla(j):S(i)] = multiplicity of vertex i in Nabla(j)'s composition factors
                comp = N[j].composition_factors().get(str(i), 0)
                assert mult[j] == comp, (order, i, j, mult[j], comp)
