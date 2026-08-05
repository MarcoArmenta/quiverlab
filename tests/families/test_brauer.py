"""Brauer graph algebras (Plan 46 / C5). Literature: dim = sum_v m_v*val(v)^2; BGAs
are symmetric. Cross-engine: the Brauer STAR (one central vertex, n truncated leaves)
is NakayamaAlgebra(n, mn+1, cyclic=True) -- byte-equal Cartan matrices; m=1 recovers
T(kA_n) = kZ_n/J^{n+1}.

PLAN CORRECTION (documented): the plan's _star sketch placed each leaf's ribbon end at
edge index i at vertex i, but edge (0,i) is `edges[i-1]`; the corrected _star places
leaf vertex k+1's single end at edge index k."""
import pytest

from quiverlab.families.brauer import BrauerGraph, BrauerGraphAlgebra
from quiverlab.families.nakayama import NakayamaAlgebra
from quiverlab.fields import QQ

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _star(n):
    # center 0, leaves 1..n; edges[k] = (0, k+1); cyclic order at 0 = edges 0..n-1;
    # each leaf vertex k+1 carries the single end of edge k.
    edges = tuple((0, i) for i in range(1, n + 1))
    cyclic = {0: tuple((k, "at0") for k in range(n))}
    for k in range(n):
        cyclic[k + 1] = ((k, "atleaf"),)
    return BrauerGraph(edges=edges, cyclic_order=cyclic)


@lit
@pytest.mark.parametrize("n, m", [(3, 1), (4, 1), (3, 2)])
def test_star_dimension_formula(n, m):
    G = _star(n)
    mult = {0: m, **{i: 1 for i in range(1, n + 1)}}
    A = BrauerGraphAlgebra(G, mult, field=QQ)
    assert A.dim == n * (m * n + 1)                     # sum_v m_v*val(v)^2 = m*n^2 + n


@xeng
@pytest.mark.parametrize("n, m", [(3, 1), (4, 1), (5, 1), (3, 2)])
def test_brauer_star_is_symmetric_nakayama(n, m):
    G = _star(n)
    mult = {0: m, **{i: 1 for i in range(1, n + 1)}}
    A = BrauerGraphAlgebra(G, mult, field=QQ)
    N = NakayamaAlgebra(n=n, l=m * n + 1, cyclic=True, field=QQ)
    assert A.dim == N.dim
    assert A.cartan_matrix() == N.cartan_matrix()      # byte-equal Cartan (arbiter)


@lit
def test_brauer_graph_algebras_are_symmetric():
    G = _star(3)
    A = BrauerGraphAlgebra(G, {0: 2, 1: 1, 2: 1, 3: 1}, field=QQ)
    assert A.is_symmetric() is True                    # Plan-29 trace-form certifier


@lit
def test_line_graph_dimension_certifies():
    # path graph u=0 -e0- v=1 -e1- w=2 : val(v)=2 (middle), leaves truncated.
    # dim = m_u*1 + m_v*4 + m_w*1 = 1+4+1 = 6 (a Brauer tree = Nakayama kZ_2/J^3).
    edges = ((0, 1), (1, 2))
    cyclic = {0: ((0, "u"),), 1: ((0, "vL"), (1, "vR")), 2: ((1, "w"),)}
    G = BrauerGraph(edges=edges, cyclic_order=cyclic)
    A = BrauerGraphAlgebra(G, {0: 1, 1: 1, 2: 1}, field=QQ)
    assert A.dim == 6
    assert A.is_symmetric() is True


@lit
def test_dimension_certificate_refuses_on_bad_wiring():
    from quiverlab.errors import QuiverlabError
    bad = BrauerGraph(edges=((0, 1),), cyclic_order={0: (), 1: ()})   # no ends ordered
    with pytest.raises(QuiverlabError):
        BrauerGraphAlgebra(bad, {0: 1, 1: 1}, field=QQ)
