"""Orbit dimension + Voigt rigidity (Plan 49 / C8). Self-cert: the orbit-dim
identity dim O_M = sum d_v^2 - dim End(M); the Voigt codim identity
dim Rep - dim O_M == dim Ext^1(M,M) on hereditary; is_rigid one call.
Literature: every indecomposable over a Dynkin path algebra is rigid (codim 0)."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.invariants.geometry import (group_dim, is_rigid,
                                           orbit_dimension, rigidity_codim,
                                           representation_variety_dim)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@selfcert
def test_orbit_dim_identity():
    A = _kA2()
    for M in (A.projective(1), A.simple(1), A.simple(2)):
        dv = M.dimension_vector()
        assert orbit_dimension(M) == group_dim(dv) - A.hom(M, M)   # end_dim == hom(M,M)


@xeng
def test_voigt_codim_identity_on_hereditary():
    # hereditary => codim O_M in Rep(Q,d) == dim Ext^1(M, M) (Voigt).
    A = linear_path_algebra(3, field=QQ)
    for v in (1, 2, 3):
        for M in (A.simple(v), A.projective(v), A.injective(v)):
            dv = M.dimension_vector()
            codim = representation_variety_dim(A, dv) - orbit_dimension(M)
            assert codim == A.ext(M, M, 1) == rigidity_codim(M)
            # and the P38 tie: dim Rep - dim GL(d) == -<d,d> = -tits_form
            dlist = [dv[w] for w in A.quiver.vertices]
            assert representation_variety_dim(A, dv) - group_dim(dv) == -A.tits_form(dlist)


@lit
def test_every_dynkin_indecomposable_is_rigid():
    A = linear_path_algebra(4, field=QQ)          # kA4, hereditary Dynkin
    ar = A.ar_quiver()
    assert ar.is_complete
    for vtx in ar.vertices:                        # all 10 indecomposables
        M = vtx["module"]
        assert is_rigid(M)                         # real-root modules: Ext^1(M,M)=0
        dv = M.dimension_vector()
        assert representation_variety_dim(A, dv) - orbit_dimension(M) == 0   # codim 0 (open)


@selfcert
def test_non_rigid_has_positive_codim():
    # the 2-Kronecker generic (1,1)-module is NOT rigid: Ext^1(M,M) = 1 (isotropic delta).
    # NB: A.module arrow matrices are full n x n in the vertex-ordered basis (the
    # verified library convention); the arrow a:1->2 sends the v1 slot to the v2 slot.
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    M = A.module({1: 1, 2: 1}, {"a": [[0, 0], [1, 0]], "b": [[0, 0], [1, 0]]})
    assert is_rigid(M) is False
    assert rigidity_codim(M) == 1


@selfcert
def test_orbit_dim_field_parity():
    # orbit dim / rigidity are field-independent for these examples (QQ vs GF(32003)).
    for field in (QQ, GF(32003)):
        A = linear_path_algebra(3, field=field)
        M = A.projective(1)
        assert orbit_dimension(M) == group_dim(M.dimension_vector()) - A.hom(M, M)
        assert is_rigid(M) is True
