"""Plan 34 (Marco, 2026-07-26): the auto->CS depth fallback.

Marco hit ``DepthLimitError: bar coboundary d^4 ... deeper engines ... arrive in
later phases`` on a dim-5 algebra over CC -- a stale hint AND a dead end, since
every deeper engine had long shipped.  Amended contract: ``engine="auto"`` falls
back to Chouhy-Solotar at the exact depth where bar/fast would raise, for
quiver-presented algebras, recorded in the dispatch trace.  In-window behaviour
is byte-unchanged; explicit engines keep their honest walls; presentation-less
algebras still refuse.  Oracle: the k[x]/(x^a) closed form HH^0 = a,
HH^n = a - 1 (n >= 1) over a characteristic-0 field [bank/CS battery]."""
import pytest

from quiverlab import CC, GF, Quiver, truncated_polynomial
from quiverlab.errors import DepthLimitError, FieldError, QuiverlabError

pytestmark = [pytest.mark.oracle_crossengine]


def test_marco_case_dim5_cc_computes_past_the_bar_wall():
    A = truncated_polynomial(5, field=CC)
    t = A.hochschild_cohomology(6)                      # died at top=4 before
    assert t.dims == [5, 4, 4, 4, 4, 4, 4]
    assert "Chouhy-Solotar" in (t.engine or "")
    h = A.hochschild_homology(6)
    assert h.dims == [5, 4, 4, 4, 4, 4, 4]


def test_fallback_recorded_in_dispatch_never_silent():
    A = truncated_polynomial(5, field=CC)
    rec = []
    A.hochschild_cohomology(6, trace=rec)
    reasons = " | ".join(getattr(e, "reason", "") for e in rec)
    assert "exceeded max_cells" in reasons and "Chouhy-Solotar" in reasons


def test_fast_gfp_wall_also_falls_back():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x^4"], field=GF(3))
    t = A.hochschild_cohomology(30, max_cells=1000)     # fast wall at tiny budget
    assert len(t.dims) == 31 and "Chouhy-Solotar" in (t.engine or "")


def test_explicit_engines_keep_their_honest_walls():
    A = truncated_polynomial(5, field=CC)
    with pytest.raises(DepthLimitError):
        A.hochschild_cohomology(6, engine="bar")


def test_presentationless_algebra_still_refuses():
    # Raw structure constants of k[x]/(x^5): dim 5, bar basis 5*4^n -- hits the
    # wall, and with NO quiver the fallback must NOT fire (CS needs the
    # presentation): the DepthLimitError stays, honestly.
    from quiverlab.core.algebra import Algebra
    n = 5
    T = [[[1 if t == i + k else 0 for t in range(n)] for k in range(n)]
         for i in range(n)]                              # x^i * x^k = x^(i+k), truncated
    A = Algebra.from_structure_constants(T, [1, 0, 0, 0, 0], field=CC)
    assert A.quiver is None
    with pytest.raises(DepthLimitError):
        A.hochschild_cohomology(12, max_cells=1000)
