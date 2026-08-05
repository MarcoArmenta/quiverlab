"""Flip <-> Fomin-Zelevinsky mutation (Plan 48, deep). Cross-engine: quiver_of(flip(T,a))
matches the exact matrix mutation of quiver_of(T) at vertex a, on every interior arc of
the disc/annulus zoo. Self-cert: flip is an involution."""
import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.flip import (certify_flip_mutation, exchange_matrix, flip,
                                     matrix_mutation)
from quiverlab.surfaces.qp import quiver_of
from quiverlab.surfaces.triangulation import annulus_triangulation, fan_triangulation

xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


@selfcert
@pytest.mark.parametrize("factory", [lambda: fan_triangulation(6),
                                     lambda: fan_triangulation(7),
                                     lambda: annulus_triangulation(2, 2)])
def test_flip_is_an_involution(factory):
    T = factory()
    for a in T.arcs():
        T2 = flip(flip(T, a), a)
        assert set(T2.arcs()) == set(T.arcs())
        # same quiver back (up to arrow renaming): compare exchange matrices
        assert exchange_matrix(quiver_of(T2)) == exchange_matrix(quiver_of(T))


@xeng
@pytest.mark.parametrize("factory", [lambda: fan_triangulation(6),
                                     lambda: fan_triangulation(7),
                                     lambda: annulus_triangulation(2, 1)])
def test_every_flip_certifies_against_matrix_mutation(factory):
    T = factory()
    for a in T.arcs():
        assert certify_flip_mutation(T, a)


@xeng
def test_matrix_mutation_is_an_involution_on_B():
    Q = quiver_of(fan_triangulation(6))
    B = exchange_matrix(Q)
    k = 1
    assert matrix_mutation(matrix_mutation(B, k), k) == B


@selfcert
def test_flip_of_boundary_segment_refused():
    T = fan_triangulation(6)
    with pytest.raises(QuiverlabError):
        flip(T, "b0_0")                                  # boundary segment: not flippable
