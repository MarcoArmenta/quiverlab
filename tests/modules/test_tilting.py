"""Tilting modules (Plan 44 / C7). Classical n=1 uses the Bongartz COUNT criterion
(#indec summands = #vertices, ASS VI.4); Bongartz completion is self-certified (T (+) E
tilting). Over QQ -- the summand count leans on decompose (char caveat)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.morphism import direct_sum
from quiverlab.modules.tilting import (bongartz_completion, is_cotilting_module,
                                       is_tilting_module)


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@pytest.mark.oracle_selfcert
def test_regular_module_is_always_tilting():
    A = _kA2()
    reg, _, _ = direct_sum(A.projective(1), A.projective(2))
    rep = is_tilting_module(reg)
    assert rep.is_tilting is True and bool(rep) is True
    assert rep.pd == 0 and rep.self_ext_vanishes and rep.num_summands == rep.num_vertices


@pytest.mark.oracle_literature
def test_apr_tilt_of_kA2():
    # T = P1 (+) S1 is a tilting module for kA2 (hereditary => pd S1 = 1, Ext^1(T,T)=0,
    # 2 non-iso summands = 2 vertices). ASS VI worked example (the APR tilt at vertex 1).
    A = _kA2()
    T, _, _ = direct_sum(A.projective(1), A.simple(1))
    assert is_tilting_module(T).is_tilting is True


@pytest.mark.oracle_selfcert
def test_partial_tilting_not_complete_then_bongartz_completes():
    A = _kA2()
    P1 = A.projective(1)
    assert is_tilting_module(P1).is_tilting is False       # only 1 summand, 2 vertices
    E = bongartz_completion(P1)
    T, _, _ = direct_sum(P1, E)
    assert is_tilting_module(T).is_tilting is True          # the certificate IS acceptance


@pytest.mark.oracle_selfcert
def test_bongartz_completes_a_nonprojective_partial_tilting():
    # S1 over kA2 is partial tilting (pd 1, Ext^1(S1,S1)=0) with d = dim Ext^1(S1, A_A) = 1,
    # so bongartz_completion exercises ONE genuine universal extension 0 -> A -> E -> S1 -> 0
    # (the d>0 path, unlike the projective P1 above). T = S1 (+) E must be tilting.
    A = _kA2()
    S1 = A.simple(1)
    assert is_tilting_module(S1).is_tilting is False        # partial (1 summand, 2 vertices)
    E = bongartz_completion(S1)
    T, _, _ = direct_sum(S1, E)
    assert is_tilting_module(T).is_tilting is True


@pytest.mark.oracle_selfcert
def test_cotilting_is_dual():
    A = _kA2()
    reg, _, _ = direct_sum(A.projective(1), A.projective(2))
    # D(A_A) = the injective cogenerator is cotilting.
    assert is_cotilting_module(reg.dualize()).is_tilting is True


@pytest.mark.oracle_selfcert
def test_self_ext_nonvanishing_rejected():
    # Kronecker K = (1 => 2) with two arrows is hereditary (pd <= 1). The regular module
    # R = (1,1) with both arrows the identity has Ext^1(R,R) = 1 != 0: Euler
    # <(1,1),(1,1)> = 1 + 1 - 2 = 0 and dim Hom(R,R) = 1, so dim Ext^1 = 1 - 0 = 1.
    # pd(R) = 1 <= 1, so the pd clause PASSES and the failing clause is the self-Ext one.
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    # full vertex-ordered 2x2 action (basis: index 0 = vertex 1, index 1 = vertex 2);
    # a, b both send the vertex-1 generator to the vertex-2 generator (the identity block).
    zo = [[0, 0], [1, 0]]
    R = K.module({1: 1, 2: 1}, {"a": zo, "b": zo})
    rep = is_tilting_module(R)
    assert rep.is_tilting is False
    assert rep.self_ext_vanishes is False
    assert "Ext" in rep.note


@pytest.mark.oracle_selfcert
def test_algebra_wrappers():
    A = _kA2()
    reg, _, _ = direct_sum(A.projective(1), A.projective(2))
    assert A.is_tilting_module(reg).is_tilting is True
    E = A.bongartz_completion(A.projective(1))
    T, _, _ = direct_sum(A.projective(1), E)
    assert A.is_tilting_module(T).is_tilting is True
