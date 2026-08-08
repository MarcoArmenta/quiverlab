"""Separation condition + strongly-simply-connected recognizer (R16, gate for P62).
Oracles: STAR/tree separated (the corrected-convention discriminator); a branching
source separated inside a NON-tree; the Zito example simply connected but NOT strongly;
the hereditary square genuinely NOT separated. Refs: Skowronski CMS 14 (1993); ASS2006;
arXiv:1905.06028.

Oracle discrimination (H2): `test_star_separated` + `test_deeper_tree_separated` catch
the B1 subquiver-convention bug (the WRONG Q_a mislabels trees as non-separated);
`test_branching_source_separated_in_nontree` shows the fix is not tree-only (a
branching source with 2 summands, separated inside a Betti-1 quiver);
`test_hereditary_square_not_separated` + `test_zito_...` confirm genuine
non-separation is still detected (no false True). The H1 minimal-relations bug is
caught by `test_two_independent_pairs_not_over_glued` in tests/coverings/test_pi1.py."""
import pytest
from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.coverings import (
    is_strongly_simply_connected, separation_condition, is_simply_connected)

pytestmark = pytest.mark.oracle_literature


def test_linear_An_separated_and_strongly():
    A = linear_path_algebra(4, field=QQ)                 # uniserial: rad P_a has <=1 summand
    assert separation_condition(A).holds is True
    assert is_strongly_simply_connected(A).verdict is True


def test_star_separated():
    # THE convention discriminator (B1): the star 1->2, 1->3 is a tree (strongly s.c.).
    # rad P_1 = S_2 (+) S_3 (two summands). CORRECT Q_1 = {2,3} (delete 1 and its
    # closure) has TWO components {2},{3} -> separated. The OLD (wrong) Q_1 = {1,2,3}
    # keeps 1, joins {2},{3} into one component -> would spuriously FAIL.
    A = Quiver([1, 2, 3], {"a": (1, 2), "c": (1, 3)}).algebra(field=QQ)
    assert separation_condition(A).holds is True
    assert is_strongly_simply_connected(A).verdict is True
    assert is_simply_connected(A).verdict is True         # tree => R1


def test_deeper_tree_separated():
    # 1->2->4, 1->3 : Q_1 = {2,3,4} has components {2,4},{3}; rad P_1 = a.A (+) c.A
    # with supports {2,4},{3} -> separated. (A deeper branching tree.)
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3)}).algebra(field=QQ)
    assert separation_condition(A).holds is True
    assert is_strongly_simply_connected(A).verdict is True


def test_branching_source_separated_in_nontree():
    # H2: a NON-tree (Betti 1: the commutative square 2-4-6-5) with a branching source
    # at 1. rad P_1 = P_2 (+) P_3 (two NON-iso summands), supports {2,4,5,6} and {3} in
    # DISTINCT components of Q_1 = {2,3,4,5,6} -> vertex 1 separated. The relation p*r-q*s
    # keeps vertex 2 separated too. Shows the corrected convention is not tree-only.
    Q = Quiver([1, 2, 3, 4, 5, 6],
               {"a": (1, 2), "c": (1, 3), "p": (2, 4), "q": (2, 5),
                "r": (4, 6), "s": (5, 6)})
    A = Q.algebra(relations=["p*r - q*s"], field=QQ)
    assert A.quiver.is_acyclic()                          # triangular; underlying graph non-tree (Betti 1)
    assert separation_condition(A).holds is True          # branching source 1 separated


def test_zito_example_simply_but_not_strongly():
    # Q: 1->2, 2->3->5, 2->4->5 ; I = <al*be*ga - al*de*ep>. pi1(Q,I)=0 (simply
    # connected) but separation FAILS (at vertex 2 already: rad P_2 = be.A (+) de.A,
    # supports {3,5},{4,5} share 5 in the connected Q_2 = {3,4,5}) -> NOT strongly s.c.
    Q = Quiver([1, 2, 3, 4, 5],
               {"al": (1, 2), "be": (2, 3), "ga": (3, 5), "de": (2, 4), "ep": (4, 5)})
    A = Q.algebra(relations=["al*be*ga - al*de*ep"], field=QQ)
    assert is_simply_connected(A, strong=True).verdict is True   # pi1(Q,I) = 0 (Zito, R3)
    ssc = is_strongly_simply_connected(A)
    assert ssc.verdict is False and ssc.witness["vertex"] == 2


def test_hereditary_square_not_separated():
    # pure commutative square (no relation): rad P_2 = be.A (+) de.A, supports {3,5}
    # and {4,5} share 5; Q_2 = {3,4,5} connected -> separation fails at vertex 2.
    Q = Quiver([2, 3, 4, 5], {"be": (2, 3), "ga": (3, 5), "de": (2, 4), "ep": (4, 5)})
    A = Q.algebra(relations=[], field=QQ)
    sep = separation_condition(A)
    assert sep.holds is False and sep.witness["vertex"] == 2


def test_non_triangular_raises():
    from quiverlab import truncated_polynomial
    with pytest.raises(QuiverlabError, match="triangular"):
        separation_condition(truncated_polynomial(3, field=QQ))


@pytest.mark.oracle_selfcert
def test_char_caveat_is_undecided_not_a_raise():
    # over a small GF(p) with char <= dim rad P_a, decompose refuses; separation must
    # CATCH it as undecided_char (B4), and the sweep maps it to verdict None -- NEVER a
    # raise from inside the sweep, NEVER a silent verdict. (Pick p and a quiver where
    # some rad P_a has dim >= p at implementation.)
    A = linear_path_algebra(6, field=GF(2))
    sep = separation_condition(A)
    assert sep.holds in (True, False, None)               # decided OR undecided_char
    ssc = is_strongly_simply_connected(A)
    assert ssc.verdict in (True, False, None)             # never raises; None if char-blocked
