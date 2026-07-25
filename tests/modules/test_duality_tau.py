"""Auslander-Reiten translates tau / tau^- via A^op + duality D + transpose Tr
(Plan 23, Tier 1b item 1). Theory/literature oracles, no [qpa] extra needed.

Oracles:
  * hereditary Coxeter transformation dim tau M = Phi^{-T} dim M (ASS);
  * explicit kA_2 / kA_3 AR tables (values cross-checked live against QPA DTr/TrD);
  * tau(projective)=0, tau^-(injective)=0;
  * tau^- tau M ~= M for non-projective indecomposables;
  * (A^op)^op == A, D.D == id, D(projective(A^op,v)) ~= injective(A,v).
"""
import sympy
import pytest

from quiverlab import Quiver, CC, GF, linear_path_algebra
from quiverlab.fields import QQ


def _cn(n, field=CC):
    """Cyclic Nakayama kZ_n / rad^2 (self-injective): vertices 1..n in a cycle."""
    verts = list(range(1, n + 1))
    arrows = {f"a{i}": (i, i % n + 1) for i in range(1, n + 1)}
    rels = [f"a{i}*a{i % n + 1}" for i in range(1, n + 1)]
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


# --- opposite algebra / duality involutions --------------------------------
def test_opposite_is_involution():
    A = linear_path_algebra(3)
    Aop = A.opposite()
    App = Aop.opposite()
    assert App is A                        # cached cross-link => exact identity
    # reversed arrows
    assert Aop.quiver.arrows["a1"] == (2, 1)
    assert Aop.quiver.arrows["a2"] == (3, 2)


def test_duality_is_involution():
    A = linear_path_algebra(3)
    for M in (A.simple(2), A.projective(1), A.injective(2)):
        DDM = M.dualize().dualize()
        assert DDM.algebra is A
        assert DDM.dimension_vector() == M.dimension_vector()
        assert M.is_isomorphic(DDM)


def test_D_of_opposite_projective_is_injective():
    # builders.injective's implicit D must agree with the explicit D on A^op.
    for field in (CC, GF(5)):
        A = _square(field=field)
        Aop = A.opposite()
        for v in A.quiver.vertices:
            explicit = Aop.projective(v).dualize()      # D(P_v^{op}) as a right A-module
            implicit = A.injective(v)                   # builders.injective = D(A e_v)
            assert explicit.dimension_vector() == implicit.dimension_vector()
            assert explicit.is_isomorphic(implicit), f"D disagreement at v={v}"


# --- tau tables (explicit AR quiver values) --------------------------------
@pytest.mark.parametrize("field", [CC, GF(2), GF(7)])
def test_ka2_tau_table(field):
    A = linear_path_algebra(2, field=field)
    S1, S2 = A.simple(1), A.simple(2)
    P1 = A.projective(1)
    # kA_2 indecomposables: S_2=P_2 (proj), P_1 (proj), S_1=I_1 (inj, non-proj)
    assert S1.tau().dimension_vector() == {1: 0, 2: 1}     # tau S_1 = S_2
    assert S2.tau().dimension_vector() == {1: 0, 2: 0}     # S_2 projective => 0
    assert P1.tau().dimension_vector() == {1: 0, 2: 0}     # projective => 0
    assert S1.tau_minus().dimension_vector() == {1: 0, 2: 0}   # S_1 injective => 0
    assert S2.tau_minus().dimension_vector() == {1: 1, 2: 0}   # tau^- S_2 = S_1


@pytest.mark.parametrize("field", [CC, GF(3)])
def test_ka3_tau_table(field):
    A = linear_path_algebra(3, field=field)
    S = {v: A.simple(v) for v in (1, 2, 3)}
    # live-QPA-confirmed DTr table for 1->2->3
    assert S[1].tau().dimension_vector() == {1: 0, 2: 1, 3: 0}   # tau S_1 = S_2
    assert S[2].tau().dimension_vector() == {1: 0, 2: 0, 3: 1}   # tau S_2 = S_3
    assert S[3].tau().dimension_vector() == {1: 0, 2: 0, 3: 0}   # S_3 projective
    # TrD table
    assert S[1].tau_minus().dimension_vector() == {1: 0, 2: 0, 3: 0}  # S_1 injective
    assert S[2].tau_minus().dimension_vector() == {1: 1, 2: 0, 3: 0}  # tau^- S_2 = S_1
    assert S[3].tau_minus().dimension_vector() == {1: 0, 2: 1, 3: 0}  # tau^- S_3 = S_2


@pytest.mark.parametrize("n", [2, 3, 4])
def test_hereditary_coxeter_transformation(n):
    # ASS: for hereditary kA_n and M non-projective indecomposable,
    # dim tau M = Phi^{-T} dim M (Phi = quiverlab coxeter_matrix, e_iAe_j convention).
    A = linear_path_algebra(n)
    verts = list(A.quiver.vertices)
    Phi = sympy.Matrix(A.coxeter_matrix())
    cox = Phi.inv().T
    for v in verts:
        S = A.simple(v)
        tau = S.tau()
        if tau.dim == 0:
            continue                        # projective simple
        dv = sympy.Matrix([S.dimension_vector()[w] for w in verts])
        predicted = cox * dv
        got = [tau.dimension_vector()[w] for w in verts]
        assert [int(x) for x in predicted] == got


# --- tau^- tau ~= M for non-projective indecomposables ---------------------
@pytest.mark.parametrize("field", [CC, GF(2)])
def test_tauminus_tau_roundtrip_hereditary(field):
    A = linear_path_algebra(3, field=field)
    for v in (1, 2):                        # S_1, S_2 are non-projective in kA_3
        S = A.simple(v)
        back = S.tau().tau_minus()
        assert back.is_isomorphic(S)


def test_tau_roundtrip_selfinjective_nakayama():
    A = _cn(3)
    # kZ_3/rad^2 self-injective: all simples non-projective; tau permutes them.
    orbit = []
    S = A.simple(1)
    cur = S
    for _ in range(3):
        cur = cur.tau()
        assert cur.dim == 1                 # stays simple
        orbit.append(tuple(sorted(cur.dimension_vector().items())))
    # tau^3 returns to S_1 (cyclic orbit of length 3)
    assert A.simple(1).is_isomorphic(A.simple(1).tau().tau().tau())
    # tau^- tau = id on each simple (non-projective)
    for v in (1, 2, 3):
        assert A.simple(v).tau().tau_minus().is_isomorphic(A.simple(v))


def test_tau_outputs_are_genuine_modules():
    A = _square(field=GF(5))
    for v in A.quiver.vertices:
        for X in (A.simple(v).tau(), A.simple(v).tau_minus(),
                  A.simple(v).dualize(), A.simple(v).transpose()):
            ok, why = X.check_module()
            assert ok, why
