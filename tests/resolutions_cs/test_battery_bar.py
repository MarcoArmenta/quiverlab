"""Validation battery (b): CS resolution == normalized-bar oracle, degreewise, on
NON-monomial algebras (spec-§6; Plan-04 Task 9).

Two INDEPENDENT engines are compared, degree by degree, over several fields:

  * CS side:  cs_(co)homology_dims(A, top) -> HH(A) from the Chouhy-Solotar
    resolution built on A's Plan-03 Groebner reduction system, with the binding
    d^2 = 0 and order-condition gates run first (resolutions_cs/homology.py).
  * bar side: hochschild_(co)homology_dims(A, top) -> HH(A) from the exponential
    normalized bar complex (hochschild/bar.py), disjoint code.

The bar dimensions are computed LIVE at assert time (function calls on the same
Algebra A), never hardcoded constants, so agreement is a genuine cross-check of
two disjoint differentials -- and the load-bearing cell of Plan-04's standing
CS<->bar deferral (T7 review).

===========================================================================
TRANSCRIPTION DEVIATIONS FROM p04-task-9-brief.md (flagged per brief's
minimal-fix clause; the brief's transcribed code fails its OWN tests here for
grammar/depth reasons, NOT for any oracle mismatch -- the oracle was never
loosened):

  (1) qci relation string.  The brief writes `f"y*x - ({xi})*x*y"`.  The Plan-03
      relation grammar (quiverlab.combinat.relations) rejects parenthesized
      coefficients -- "y*x - (2)*x*y" raises RelationError "unknown arrow '(2)'".
      Parens removed, matching the documented workaround in this package's
      conftest.qci_rs.  (Grammar extension is Plan-06 territory.)

  (2) qci xi = "i" (the intended Q(i) number-field case) is UNBUILDABLE via the
      public relation surface.  Coefficients flow through parse_rational and are
      stored as fractions.Fraction (relations.Relation.terms); the grammar has no
      hook for field-specific literals, so 'i' parses as an (unknown) arrow name,
      and no relation-derived Algebra can ever carry Q(i)-valued structure
      constants in this phase.  The case is kept VISIBLE as an xfail-style skip
      documenting the intended coverage and its Plan-06 blocker; it is NOT
      silently dropped and NOT faked with a rational stand-in.  (Computing the
      other qci cases over CC still exercises the exact-complex CC domain.)

  (3) qci xi = "-1".  The bare form "y*x - -1*x*y" is a double sign the grammar
      rejects ("malformed term ''").  Rewritten to the identical relation
      "y*x + x*y" (q = -1).

  (4) square depth 4 -> 2.  The square algebra has dim 9, so the bar cochain
      C^n = 9*8^n; the coboundary d^3 is 36864 x 4608 (~1.7e8 entries) and trips
      hochschild/bar.py's max_cells=4_000_000 DepthLimitError guard (and exact CC
      rank at that size is infeasible).  top is reduced to the largest value
      inside the bar window, 2 -- same minimal-fix pattern as Task 8's
      local_radsq depth reduction.  All four fields are retained (CC gives char 0,
      GF(2) char 2), so the "char 0 AND char 2" cross-check still lands.

  (5) kx2 added (char 0 = CC AND char 2 = GF(2)).  The launching controller's
      standing-deferral minimum names kx2 explicitly; the brief omitted it.  Added
      here and flagged.  kx2 = k[x]/(x^2) is monomial (a base sanity anchor); the
      non-monomial content is carried by the square (binomial a*b - c*d) and the
      qci (binomial y*x -/+ x*y).

Net count: 4 (square) + 2 run + 1 skip (qci) + 2 (kx2) = 8 passed, 1 skipped.
===========================================================================
"""
import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
from quiverlab.hochschild.bar import hochschild_cohomology_dims, hochschild_homology_dims
pytest.importorskip("quiverlab.groebner")


def _square(field):
    # Commutative-square quiver with the single BINOMIAL (non-monomial) relation.
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=field)


def _qci(field, rel):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    return Q.algebra(relations=["x*x", "y*y", rel], field=field)


def _kx2(field):
    Q = Quiver([1], {"x": (1, 1)})
    return Q.algebra(relations=["x*x"], field=field)


# Buildable quantum parameters (see deviation notes 2 and 3).
_QCI_REL = {"2": "y*x - 2*x*y", "-1": "y*x + x*y"}


@pytest.mark.parametrize("field", [CC, GF(2), GF(3), GF(5)])
def test_cs_equals_bar_square(field):
    A = _square(field)                       # deviation (4): top 4 -> 2 (bar window)
    assert cs_cohomology_dims(A, 2).dims == hochschild_cohomology_dims(A, 2).dims
    assert cs_homology_dims(A, 2).dims == hochschild_homology_dims(A, 2).dims


@pytest.mark.parametrize("xi", ["2", "i", "-1"])
def test_cs_equals_bar_qci(xi):
    if xi == "i":                            # deviation (2): Q(i) unreachable, Plan-06
        pytest.skip("Q(i) coefficients unreachable via the Plan-03 relation grammar "
                    "(rationals only); number-field literals are Plan-06 territory")
    A = _qci(CC, _QCI_REL[xi])               # CC exact-complex domain exercised
    assert cs_homology_dims(A, 3).dims == hochschild_homology_dims(A, 3).dims
    assert cs_cohomology_dims(A, 3).dims == hochschild_cohomology_dims(A, 3).dims


@pytest.mark.parametrize("field", [CC, GF(2)])
def test_cs_equals_bar_kx2(field):           # deviation (5): standing-deferral anchor
    A = _kx2(field)                          # char 0 (CC) AND char 2 (GF(2))
    assert cs_cohomology_dims(A, 6).dims == hochschild_cohomology_dims(A, 6).dims
    assert cs_homology_dims(A, 6).dims == hochschild_homology_dims(A, 6).dims
