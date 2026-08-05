"""Euler/Tits forms off the Cartan matrix. Literature: Gabriel's theorem
(finite type <=> Dynkin <=> positive definite Tits form), ASS Ch. VII."""
import sympy as sp
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.forms import (euler_form, euler_form_matrix,
                                        form_type, tits_form)

pytestmark = pytest.mark.oracle_literature


def _kronecker(m):
    arrows = {f"a{i}": (1, 2) for i in range(m)}
    return Quiver([1, 2], arrows).algebra()


def test_hereditary_euler_form_is_the_arrow_formula():
    # <d,e> = sum d_v e_v - sum_{a: s->t} d_s e_t   (hereditary)
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra()
    d, e = [1, 2, 1], [3, 1, 2]
    expected = (1*3 + 2*1 + 1*2) - (1*1 + 2*2)
    assert euler_form(A, d, e) == expected


def test_euler_form_computes_hom_minus_ext_dimension():
    # homological meaning, pinned on a hereditary example:
    # <dim M, dim N> = dim Hom(M,N) - dim Ext^1(M,N)
    A = Quiver([1, 2], {"a": (1, 2)}).algebra()
    S1, S2 = A.simple(1), A.simple(2)
    d1 = [1, 0]; d2 = [0, 1]
    assert euler_form(A, d1, d2) == A.hom(S1, S2) - A.ext(S1, S2, 1)
    assert euler_form(A, d2, d1) == A.hom(S2, S1) - A.ext(S2, S1, 1)


def test_form_type_dynkin_euclidean_wild():
    assert form_type(Quiver([1, 2], {"a": (1, 2)}).algebra()) == "finite"   # A2
    assert form_type(_kronecker(2)) == "tame"                               # ~A1
    assert form_type(_kronecker(3)) == "wild"                               # 3-Kronecker


def test_tits_form_values():
    # A2: q(d) = d1^2 + d2^2 - d1 d2; q(1,1) = 1 (a root!)
    A = Quiver([1, 2], {"a": (1, 2)}).algebra()
    assert tits_form(A, [1, 1]) == 1
    # Kronecker: q(1,1) = 0 (the isotropic imaginary root)
    assert tits_form(_kronecker(2), [1, 1]) == 0


def test_singular_cartan_refused():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"])   # k[x]/(x^2)
    with pytest.raises(QuiverlabError, match="[Cc]artan"):
        euler_form_matrix(A)
