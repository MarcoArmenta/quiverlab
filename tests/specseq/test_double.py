"""DoubleComplex: anticommutation gate, total-complex assembly, and both
filtrations. Self-certifying (Tot d.d=0 by construction; filtration is a
subcomplex filtration -- FilteredComplex re-checks it)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.specseq.double import DoubleComplex
from quiverlab.specseq.filtered import FilteredComplex

pytestmark = pytest.mark.oracle_selfcert


def _dom():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(7)).domain


def test_anticommute_gate_refuses_commuting_square():
    dom = _dom()
    # D_{1,1}=D_{0,1}=D_{1,0}=D_{0,0}=1-dim; d_h=d_v=1 everywhere => they COMMUTE
    # (d_h d_v = d_v d_h = 1), so d_h d_v + d_v d_h = 2 != 0 over GF(7): refused.
    terms = {(1, 1): 1, (0, 1): 1, (1, 0): 1, (0, 0): 1}
    d_h = {(1, 1): [[1]], (1, 0): [[1]]}
    d_v = {(1, 1): [[1]], (0, 1): [[1]]}
    with pytest.raises(QuiverlabError, match="anticommut|d_h.*d_v"):
        DoubleComplex(terms, d_h, d_v, dom)


def test_signed_square_totals_correctly():
    dom = _dom()
    # put the (-1)^p sign on d_v so the square anticommutes; Tot of the 2x2 square
    # [P (+) P -> P] has the homology of an exact square (acyclic in the middle).
    neg = dom.neg(dom.one())
    terms = {(1, 1): 1, (0, 1): 1, (1, 0): 1, (0, 0): 1}
    d_h = {(1, 1): [[1]], (1, 0): [[1]]}
    d_v = {(1, 1): [[neg]], (0, 1): [[1]]}     # sign on the p=1 column
    dc = DoubleComplex(terms, d_h, d_v, dom)
    terms_tot, dmats_tot, _ = dc.total()
    assert terms_tot == {2: 1, 1: 2, 0: 1}     # Tot_2=D11, Tot_1=D01(+)D10, Tot_0=D00
    # the total of a commuting (sign-fixed) 2x2 square with all maps iso is acyclic
    F = dc.column_filtration()
    assert F.total_homology_dims() == {2: 0, 1: 0, 0: 0}


def test_both_filtrations_are_valid_subcomplex_filtrations():
    # column_filtration / row_filtration must both pass FilteredComplex's closed +
    # exhaustive gates (they are built to; this pins that they do).
    dom = _dom()
    neg = dom.neg(dom.one())
    terms = {(1, 1): 1, (0, 1): 1, (1, 0): 1, (0, 0): 1}
    d_h = {(1, 1): [[1]], (1, 0): [[1]]}
    d_v = {(1, 1): [[neg]], (0, 1): [[1]]}
    dc = DoubleComplex(terms, d_h, d_v, dom)
    assert isinstance(dc.column_filtration(), FilteredComplex)   # constructs (checks)
    assert isinstance(dc.row_filtration(), FilteredComplex)
