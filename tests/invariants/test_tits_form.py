"""The combinatorial Tits form q_A = sum x_i^2 - sum_arrows + sum r_ij, r_ij the
minimal-relation count (P56's I/(rad.I+I.rad)) cross-checked against dim Ext^2(S_i,S_j).
Refs: Bongartz Math. Ann. 269 (1984); Encyclopedia of Mathematics 'Tits quadratic form'."""
import pytest

from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.coverings import minimal_relation_counts
from quiverlab.invariants.tits import as_unit_form, tits_form_combinatorial as titsc

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _dual_numbers_sc(field=CC):
    """k[x]/(x^2) via raw structure constants: NO quiver."""
    T = [[[1, 0], [0, 1]],
         [[0, 1], [0, 0]]]
    return Algebra.from_structure_constants(T, [1, 0], field=field)


@xeng
def test_r_ij_matches_ext2_of_simples():
    # the named cross-plan contract: combinatorial count (P56) == homological Ext^2.
    # 1->2->3 with the length-2 relation a*b=0 : ONE minimal relation 1->3.
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=GF(7))
    r = minimal_relation_counts(A)
    assert r.get((1, 3), 0) == 1 and sum(r.values()) == 1
    assert A.ext(A.simple(1), A.simple(3), 2) == 1          # independent engine
    assert A.ext(A.simple(1), A.simple(2), 2) == 0


@lit
def test_kronecker_tits_values():
    # m-Kronecker (hereditary): q(x,y) = x^2 + y^2 - m x y. r_ij all 0.
    for m, val11 in [(1, 1), (2, 0), (3, -1), (4, -2)]:
        A = Quiver([1, 2], {f"a{k}": (1, 2) for k in range(m)}).algebra(field=CC)
        assert minimal_relation_counts(A) == {}
        assert titsc(A, [1, 1]) == val11                     # 2 - m


@lit
def test_A2_and_nakayama_relation_form():
    A2 = Quiver([1, 2], {"a": (1, 2)}).algebra(field=CC)
    assert titsc(A2, [1, 1]) == 1                            # x^2+y^2-xy, a root
    # kA3 with a*b=0 : q = x1^2+x2^2+x3^2 - x1x2 - x2x3 + x1x3 (one relation 1->3)
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=CC)
    assert titsc(A, [1, 1, 1]) == 3 - 1 - 1 + 1              # == 2
    # dict form agrees with the vertex-order sequence form
    assert titsc(A, {1: 1, 2: 1, 3: 1}) == 2


@xeng
def test_gldim_le_2_agrees_with_P38_euler_symmetrization():
    # W2: where BOTH are defined (gl.dim <= 2), the combinatorial Tits form == P38's
    # homological Euler symmetrization (E + E^T). The self-cert cross-check.
    import sympy as sp

    from quiverlab.invariants.forms import tits_matrix as euler_sym
    A = linear_path_algebra(4, field=CC)                     # hereditary, gl.dim 1
    G = as_unit_form(A).gram
    assert sp.Matrix([list(row) for row in G]) == euler_sym(A)
    # a gl.dim-2 relation algebra also coincides
    B = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=CC)
    from quiverlab.invariants.tits import tits_matrix_combinatorial
    assert tits_matrix_combinatorial(B) == euler_sym(B)


def test_non_triangular_refused():
    # M2: k[x]/(x^3) has a loop => G_11 == 2 by cancellation (a_11=1, r_11=1), so the
    # G_ii==2 heuristic would WRONGLY accept it. is_acyclic() is the honest gate.
    A = truncated_polynomial(3, field=CC)                    # k[x]/(x^3): a loop
    with pytest.raises(QuiverlabError, match="triangular|acyclic|cycle"):
        as_unit_form(A)


def test_presentationless_refused():
    A = _dual_numbers_sc(CC)                                 # NO quiver
    with pytest.raises(QuiverlabError, match="quiver|presentation"):
        minimal_relation_counts(A)
    with pytest.raises(QuiverlabError, match="quiver|presentation"):
        titsc(A, [1, 1])
    with pytest.raises(QuiverlabError, match="quiver|presentation"):
        as_unit_form(A)
