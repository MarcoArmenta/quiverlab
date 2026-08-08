"""Hom(pi1, k+) into HH^1 for triangular algebras (Assem-de la Pena; CRS Hurewicz
iso for Schurian). Cross-engine: pi1 via SNF vs HH via the shipped engines."""
import pytest
from quiverlab import GF, Quiver, linear_path_algebra, zoo
from quiverlab.invariants.coverings import fundamental_group

pytestmark = pytest.mark.oracle_crossengine

P = 32003


def _dimhh1(A):
    return A.hochschild_cohomology(1)[1]


def test_square_equality_schurian():
    # commutative square (Schurian, triangular): Hom(pi1,k+) = dim HH^1 (Hurewicz iso)
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
               ).algebra(relations=["a*b - c*d"], field=GF(P))
    g = fundamental_group(A)
    assert g.hom_to_additive_dim(P) == 0
    assert _dimhh1(A) == 0                                 # square with relation: HH^1 = 0


def test_hereditary_square_Z_bound():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
               ).algebra(relations=[], field=GF(P))
    g = fundamental_group(A)
    assert g.hom_to_additive_dim(P) == 1 <= _dimhh1(A)     # Z: Hom-dim 1 <= dim HH^1


def test_bound_holds_on_triangular_zoo():
    for A in zoo(dim_max=9, field=GF(P)):
        if A.quiver is None or not A.quiver.is_acyclic():  # embedding is a theorem for triangular only
            continue
        g = fundamental_group(A)
        assert g.hom_to_additive_dim(P) <= _dimhh1(A)
