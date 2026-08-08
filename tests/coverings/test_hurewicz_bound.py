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


def _kronecker(m):
    return Quiver([1, 2], {("x%d" % i): (1, 2) for i in range(m)}).algebra(field=GF(P))


def test_bound_on_constructed_high_hh1():
    # Widen the crossengine oracle beyond the 1-2 triangular zoo members (adjudicated
    # Fix 3): constructed triangular algebras with dim HH^1 >= 2. The m-Kronecker is
    # hereditary with dim HH^1 = m^2 - 1 (>= 3 for m >= 2); the hereditary double-double
    # (two parallel 1->2 and two parallel 2->3, NO relations) has underlying Betti 2.
    # In each case Hom(pi1, k+) = free rank <= dim HH^1 (Assem-de la Pena embedding).
    cases = [_kronecker(2), _kronecker(3)]                # dim HH^1 = 3, 8
    cases.append(Quiver([1, 2, 3],
                        {"a": (1, 2), "b": (1, 2), "c": (2, 3), "d": (2, 3)}
                        ).algebra(relations=[], field=GF(P)))   # hereditary Betti-2
    high = 0
    for A in cases:
        g = fundamental_group(A)
        d = _dimhh1(A)
        assert g.hom_to_additive_dim(P) <= d              # the Hurewicz bound
        high += 1 if d >= 2 else 0
    assert _dimhh1(_kronecker(2)) >= 2 and _dimhh1(_kronecker(3)) >= 2   # dim HH^1 >= 2
