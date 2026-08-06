"""Characteristic tilting module + Ringel dual (Plan 47). Literature (Ringel, Math. Z. 1991):
kA_n natural order T = D(A) (injective cogenerator); opposite order T = A (regular). Self-cert:
Ext^1(T,Nabla)=0 AND is_tilting_module(T). Oracle: double Ringel dual has the same Cartan
invariant factors (Smith normal form) as A. Over QQ (decompose/presented_form char-clean at
char 0). ringel_dual(A) must be PRESENTED for the double-dual (its projectives get built)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.ext import ext_dims
from quiverlab.modules.hom import is_isomorphic
from quiverlab.modules.morphism import direct_sum
from quiverlab.modules.quasihereditary import (characteristic_tilting, costandard_modules,
                                               ringel_dual, standard_modules)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _a3():
    return linear_path_algebra(3, field=QQ)


def _square():
    # commutative square 1->2->4, 1->3->4 with ab = cd: a NON-hereditary quasi-hereditary
    # algebra (gl.dim 2, every Delta simple). The discriminating oracle: the classical
    # pd<=1 tilting certificate and the single-pass T(j)-extension both got this wrong.
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b-c*d"], field=QQ)


@lit
def test_natural_order_T_is_D_of_A():
    A = _a3()
    T = characteristic_tilting(A)                       # natural order
    DA, _, _ = direct_sum(*[A.injective(v) for v in (1, 2, 3)])   # D(A) = (+) I(v)
    assert is_isomorphic(T, DA)


@lit
def test_opposite_order_T_is_regular():
    A = _a3()
    T = characteristic_tilting(A, [3, 2, 1])
    reg, _, _ = direct_sum(*[A.projective(v) for v in (1, 2, 3)])  # A_A
    assert is_isomorphic(T, reg)


@selfcert
def test_T_is_tilting_and_ext_perp_nabla():
    from quiverlab.modules.tilting import is_tilting_module     # P44
    A = _a3()
    T = characteristic_tilting(A)
    assert is_tilting_module(T).is_tilting is True
    N = costandard_modules(A)
    for j in (1, 2, 3):
        assert ext_dims(A, T, N[j], 1)[1] == 0                  # Ext^1(T, Nabla(j)) = 0


@lit
def test_commutative_square_T_is_D_of_A():
    # THE discriminating non-hereditary oracle (Dlab-Ringel). All Delta(i) are simple, but
    # T = D(A) has pd 2, so the classical pd<=1 self-cert rejected it and the single-pass
    # universal extension over T(j) built a module outside F(Nabla). The genuine Ringel
    # iteration + the gl.dim-degree certificate must land exactly on D(A).
    A = _square()
    assert all(standard_modules(A)[v].dim == 1 for v in (1, 2, 3, 4))   # all Delta simple
    T = characteristic_tilting(A)                                       # natural order
    DA, _, _ = direct_sum(*[A.injective(v) for v in (1, 2, 3, 4)])       # D(A) = (+) I(v)
    assert T.dim == 9
    assert T.dimension_vector() == {1: 4, 2: 2, 3: 2, 4: 1}
    assert is_isomorphic(T, DA)


@selfcert
def test_commutative_square_T_both_sided_ext_perp():
    # F(Delta) cap F(Nabla): Ext^1(Delta(j), T) = 0 AND Ext^1(T, Nabla(j)) = 0 for all j,
    # and T is 2-tilting (pd = gl.dim = 2, not the classical 1).
    from quiverlab.modules.tilting import is_tilting_module
    A = _square()
    T = characteristic_tilting(A)
    D = standard_modules(A)
    N = costandard_modules(A)
    for j in (1, 2, 3, 4):
        assert ext_dims(A, D[j], T, 1)[1] == 0                          # Ext^1(Delta(j), T) = 0
        assert ext_dims(A, T, N[j], 1)[1] == 0                          # Ext^1(T, Nabla(j)) = 0
    rep = is_tilting_module(T, n=2)
    assert rep.is_tilting is True and rep.pd == 2                       # NOT pd <= 1


@lit
def test_double_ringel_dual_cartan_smith_form():
    import sympy
    from sympy import ZZ
    from sympy.matrices.normalforms import invariant_factors
    A = _a3()
    R2 = ringel_dual(ringel_dual(A))

    def invf(M):
        S = sympy.Matrix([[int(x) for x in row] for row in M.cartan_matrix()])
        return list(invariant_factors(S, domain=ZZ))

    assert invf(R2) == invf(A)                                  # Ringel duality involution


@selfcert
def test_ringel_dual_of_ka3_is_a_path_algebra():
    A = _a3()
    R = ringel_dual(A)
    assert R.quiver is not None                                 # presented (char 0 => kQ/I)
    assert R.dim == A.dim                                       # Morita to a kA3-shaped algebra
