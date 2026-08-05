"""QPA as the oracle for almost-split sequences + AR predecessors (Plan 41).
qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.

Honest scope (documented): QPA's ``AlmostSplitSequence`` runs over QQ, so the
middle-term DIMENSION VECTOR is crosschecked over QQ; QPA's ``DecomposeModule`` /
``PredecessorsOfModule`` require a FINITE field, so the summand-multiset and the
immediate-predecessor multiset are crosschecked over GF(p)."""
import pytest

from quiverlab import GF, NakayamaAlgebra, linear_path_algebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_exposes_almost_split_surface():
    lg = session.libgap_handle()
    for name in ("AlmostSplitSequence", "PredecessorsOfModule"):
        assert bool(lg.eval(f'IsBoundGlobal("{name}")')), name


@pytest.mark.parametrize("A", [linear_path_algebra(3, field=QQ),
                               NakayamaAlgebra(n=3, l=2, cyclic=False, field=QQ)])
def test_almost_split_middle_vs_qpa(A):
    for v in A.quiver.vertices:
        M = A.simple(v)
        if M.tau().dim == 0:                     # projective end: no AR sequence
            continue
        A.crosscheck("almost_split", M).assert_agree()


@pytest.mark.parametrize("A", [linear_path_algebra(3, field=GF(5)),
                               NakayamaAlgebra(n=3, l=2, cyclic=False, field=GF(5))])
def test_almost_split_summands_and_predecessors_vs_qpa(A):
    # over GF(p) QPA's DecomposeModule / PredecessorsOfModule are available, so the
    # summand multiset (in the almost_split crosscheck) and the immediate-predecessor
    # multiset are compared too.
    for v in A.quiver.vertices:
        M = A.simple(v)
        if M.tau().dim == 0:
            continue
        A.crosscheck("almost_split", M).assert_agree()
        A.crosscheck("predecessors", M).assert_agree()


def test_almost_split_non_brick_middle_vs_qpa():
    # Devil's-advocate round (Finding B): the NON-BRICK case the battery lacked. Over
    # k[x]/(x^4) (GF(32003)), M = k[x]/(x^2) has dim Ext^1(M, tau M) = 2, so the socle pick
    # -- not the three sanity guards -- decides the middle term. QPA's AlmostSplitSequence +
    # DecomposeModule are the live oracle: the middle is {1, 3}, never the projective {4}.
    A = truncated_polynomial(4, field=GF(32003))
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")   # k[x]/(x^2)
    A.crosscheck("almost_split", M).assert_agree()
