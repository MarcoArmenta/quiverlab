"""Kac canonical decomposition over hereditary Dynkin (Plan 49 / C8). Literature:
(2,1) over kA2 = P1 (+) S1 (hand-derived). Self-cert: sum m_i root_i == d, each
root a positive root, pairwise ext 0 (=> generic module rigid). Refusals: the
Euclidean 2-Kronecker (deferred) and the non-hereditary k[x]/(x^2)."""
import pytest

from quiverlab import Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.geometry import canonical_decomposition

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@lit
def test_kA2_21_is_P1_plus_S1():
    A = _kA2()
    cd = canonical_decomposition(A, {1: 2, 2: 1})
    # (2,1) = (1,1) + (1,0) = P1 (+) S1 (the rank-1 generic map k^2 -> k)
    got = sorted((c["root"], c["multiplicity"]) for c in cd)
    assert got == [((1, 0), 1), ((1, 1), 1)]
    names = sorted(c["name"] for c in cd)
    assert names == ["P_1", "S_1"]


@selfcert
@pytest.mark.parametrize("d", [{1: 1, 2: 1}, {1: 2, 2: 1}, {1: 2, 2: 2}, {1: 3, 2: 2}])
def test_kA2_decomposition_is_certified(d):
    A = _kA2()
    cd = canonical_decomposition(A, d)
    verts = list(A.quiver.vertices)
    # sum of components == d
    total = {v: 0 for v in verts}
    for c in cd:
        for k, v in enumerate(verts):
            total[v] += c["multiplicity"] * c["root"][k]
    assert total == d
    # each component is a positive root
    roots = {tuple(r) for r in A.positive_roots()}
    assert all(tuple(c["root"]) in roots for c in cd)


@selfcert
def test_euclidean_is_deferred():
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    with pytest.raises(QuiverlabError, match="Dynkin|finite|deferred|Euclidean"):
        canonical_decomposition(A, {1: 1, 2: 1})     # delta = (1,1) isotropic, tame


@selfcert
def test_non_hereditary_refused():
    A = truncated_polynomial(2, field=QQ)            # k[x]/(x^2): not hereditary
    with pytest.raises(QuiverlabError, match="hereditary"):
        canonical_decomposition(A, {1: 1})
