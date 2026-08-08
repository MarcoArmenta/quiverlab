"""Commutative ladders CL(n) = A_n [] A_2 (Plan 69 / R33, Escolar-Hiraoka 2016).
Literature (oracle_literature): CL rep-finite iff n <= 4 -- the n<=4 knit terminates
complete AND CL(5) is refused (the theorem boundary, literature-confirmed verbatim).
Self-cert (oracle_selfcert): the fully-forward CL(n) = incidence algebra of [n]x[2],
dim = 3*C(n+1,2); the vertex COUNTS CL(2)=11 / CL(3)=29 are SINGLE-ENGINE values asserted
against the same AR knitter that produced them -- they rest on one engine until reconciled
against the Escolar-Hiraoka figure (BLOCKING, Task 6 Step 1a) / QPA (Task 4), so they are
oracle_selfcert NOT oracle_literature (H3); the recognizer accepts constructed ladders and
rejects non-ladders; VERTICES ARE SCALAR (tuple labels crash the knit, H1)."""
import pytest
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.families.commutative_ladder import (CommutativeLadder,
                                                   is_commutative_ladder, persistence_line)
from quiverlab.modules.ar import knit_ar_quiver

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
@pytest.mark.parametrize("n, dim", [(2, 9), (3, 18), (4, 30)])   # 3*C(n+1,2)
def test_forward_ladder_dimension(n, dim):
    A = CommutativeLadder(n, field=QQ)
    assert A.dim == dim
    ok, m = is_commutative_ladder(A)
    assert ok is True and m == n


@lit
@pytest.mark.parametrize("n", [2, 3])
def test_cl_is_rep_finite_for_n_le_4(n):
    # Escolar-Hiraoka: rep-finite for n <= 4 -- the engine corroborates by terminating
    # the knit "complete". (The vertex COUNT is a separate SELF-CERT test below; the
    # literature content here is the finite/complete BOOLEAN, not the number.)
    # n=4 is EXCLUDED from the knit here: the CL(4) (dim 30) knit did not complete within
    # 30+ min in-session (impractical for the deep bucket, whose target is ~19 min total),
    # so CL(4) rep-finiteness rests on the THEOREM (literature-confirmed n<=4) + the
    # build/dim-30 certificate (test_cl4_dim_certificate_no_knit); the knit corroboration
    # is exercised at n=2,3. See the worker report's declared deviation (Task 1).
    A = CommutativeLadder(n, field=QQ)
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)
    assert ar.is_complete and ar.status == "complete"


@lit
def test_cl4_dim_certificate_no_knit():
    # CL(4) is rep-finite by the Escolar-Hiraoka theorem (n <= 4). We certify its
    # construction (dim 30 = 3*C(5,2), scalar vertex names, recognizer) WITHOUT knitting
    # -- the CL(4) AR knit is dim-30-heavy and did not complete in 30+ min in-session, so
    # it is not a standing CI test (declared deviation, Task 1). The theorem is the
    # rep-finiteness oracle; CL(5) refusal (below) is the other half of the boundary.
    A = CommutativeLadder(4, field=QQ)
    assert A.dim == 30 and all(not isinstance(v, tuple) for v in A.quiver.vertices)
    ok, m = is_commutative_ladder(A)
    assert ok is True and m == 4


@lit
@selfcert
@pytest.mark.parametrize("n, count", [(2, 11), (3, 29)])   # FIGURE-CONFIRMED counts
def test_cl_verified_counts(n, count):
    # FIGURE-CONFIRMED literature pin (H3 reconciliation DONE): the Escolar-Hiraoka
    # AR-quiver figures (arXiv:1404.7588, DCG 2016) were rendered (pdftotext -layout +
    # pages 46-47 @150dpi) and counted -- Fig. 13 = CL(f) = this CL(2): 3+5+3 = 11 (all
    # thin); Fig. 14 = CL(ff) = this CL(3): 1+6+11+6+5 = 29, exactly 2 non-thin, matching
    # the engine's 27+2 thin/non-thin split. The equioriented ladder IS the incidence
    # algebra of the [n]x[2] grid poset (iso-invariant), so the orientation identification
    # is airtight. Kept oracle_selfcert too (the same knit reproduces it; Plan-32 overlap).
    A = CommutativeLadder(n, field=QQ)
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)
    assert ar.is_complete and len(ar.vertices) == count


@selfcert
def test_cl_vertices_are_scalar_and_knit_completes():
    # REGRESSION for H1: the natural [n]x[2] grid-poset labels (i,j) are TUPLES, which
    # crash knit_ar_quiver at complex_reps.py:179 (`"e_%s" % v`, TypeError). CommutativeLadder
    # MUST use scalar vertex names, so the knit that produced the 29 count above completes.
    A = CommutativeLadder(3, field=QQ)
    assert all(not isinstance(v, tuple) for v in A.quiver.vertices)   # no tuple labels
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)       # must NOT crash
    assert ar.is_complete and len(ar.vertices) == 29


@lit
def test_cl5_is_refused_rep_infinite():
    with pytest.raises(QuiverlabError):
        CommutativeLadder(5, field=QQ)   # Escolar-Hiraoka: rep-infinite for n >= 5


@selfcert
def test_persistence_line_orientations():
    fwd = persistence_line(5, "forward")
    assert fwd.is_acyclic() and len(fwd.arrows) == 4
    zz = persistence_line(5, "zigzag")     # 1->2<-3->4<-5
    assert zz.is_acyclic() and len(zz.arrows) == 4


@selfcert
def test_recognizer_rejects_non_ladder():
    from quiverlab import Quiver
    line = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=QQ)
    ok, _ = is_commutative_ladder(line)
    assert ok is False
