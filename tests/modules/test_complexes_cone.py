"""Cones and triangles. Self-certifying: cone d.d=0 under full validation;
LES rank identities; quasi-iso <=> cone acyclic on known cases."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import (ChainComplex, ChainMap,
                                         identity_chain_map)

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def test_cone_of_identity_is_acyclic():
    A = _a2()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=3)
    f = identity_chain_map(X)
    C = f.cone()
    assert C.is_acyclic()
    assert f.is_quasi_iso()


def test_cone_of_zero_map_is_direct_sum_shift():
    A = _a2()
    S1 = A.simple(1)
    X = ChainComplex.stalk(S1, 0)
    Y = ChainComplex.stalk(S1, 0)
    z = ChainMap(X, Y, {0: [[0]]})
    C = z.cone()
    # cone(0) = X[1] (+) Y: homology in degrees 1 and 0
    assert C.homology_dims().get(1) == 1 and C.homology_dims().get(0) == 1


def test_resolution_augmentation_is_quasi_iso():
    # the augmentation Q_* -> M (stalk) is the canonical quasi-iso
    A = _a2()
    S1 = A.simple(1)
    Q = ChainComplex.from_projective_resolution(S1, length=3)
    M = ChainComplex.stalk(S1, 0)
    from quiverlab.modules.resolution import projective_cover
    Q0, d0, _ = projective_cover(S1)
    aug = ChainMap(Q, M, {0: d0})
    assert aug.is_quasi_iso()


def test_bad_square_refused():
    A = _a2()
    P1 = A.projective(1)
    X = ChainComplex.from_projective_resolution(A.simple(1), length=2)
    with pytest.raises(QuiverlabError, match="square|chain map"):
        ChainMap(X, X, {n: [[1] * X.term(n).dim for _ in range(X.term(n).dim)]
                        for n in X.degrees() if X.term(n).dim})


def test_triangle_euler_characteristic():
    # LES consequence: chi(cone) = chi(X[1]) + chi(Y) = -chi(X) + chi(Y)
    A = _a2()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=2)
    Y = ChainComplex.stalk(A.simple(2), 0)
    z = ChainMap(X, Y, {0: [[0] * X.term(0).dim]})
    C = z.cone()
    chi = lambda Z: sum((-1) ** n * d for n, d in Z.homology_dims().items())
    assert chi(C) == -chi(X) + chi(Y)
