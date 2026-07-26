"""Plan 29 (Part 3) literature-oracle battery -- engine-level HH/HC value pins.

Companion to ``test_engine_validation.py`` / ``test_cyclic_homology.py``.  Pins the
Hochschild (co)homology + cyclic engines against published computations from the
Cibils and Solotar/Taillefer clusters and the Schremmer/Happel canonical-algebra
formula.  Every pin carries its provenance + the verification-status label from
``docs/plans/2026-07-25-literature-oracles-deep-research.md`` (BINDING labels).

Batteries here (plan Part 3, items 3/4/5):

  2a. Cibils (1998) radical-square-zero char discriminator: k[x]/(x^2) doubles in
      char 2; the 2-Kronecker rad^2=0 path algebra HH^* = [1,3,0,0] in every char.
  2b. Taillefer (2001) cyclic homology of the Taft / self-injective cyclic Nakayama
      algebras (char 0 ONLY).
  2c. Happel-via-Schremmer canonical-algebra HH (HH^2 = t-3; HH_0 = #vertices).

===========================================================================
CHAR NOTES:
  * Taillefer requires a field containing Q (char 0).  HH_* is pinned over CC
    (genuine char 0, via the CS engine which stays cheap where the bar blows up).
    Cyclic HC over CC is the exponential (b,B) mixed complex: feasible for the
    dim-4 Taft_2 (to degree 5) but NOT the dim-9 Taft_3.  For Taft_3 the HC pattern
    is pinned over GF(32003) as an explicit CHARACTERISTIC-0 PROXY -- the same
    convention ``test_cyclic_homology.py`` already uses ("characteristic-0 proxy
    (p = 32003)"); p = 32003 divides neither n nor n-1 for n in {2,3}, so the proxy
    is faithful.  Depth is capped at 2 for Taft_3 (HC to degree 3 costs ~4 GB --
    OOM-risky in the deep bucket); [3,2,3] already exhibits the even/odd (n / n-1)
    alternation decisively.
  * Bergh-Erdmann char-p homology branches are NOT pinned (need infinite char-p
    fields) -- noted in the sibling resolutions_cs battery.
===========================================================================
"""
import pytest

from quiverlab import Quiver, CC, GF
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
from quiverlab.families import NakayamaAlgebra

pytest.importorskip("quiverlab.groebner")


# =========================================================================== #
# 2a.  Cibils (1998) radical-square-zero char discriminator                     #
# =========================================================================== #
# Citation: C. Cibils, "Hochschild cohomology algebra of radical square zero
# algebras," in Algebras and Modules II (Geiranger 1996), CMS Conf. Proc. 24,
# AMS, 1998, 93-101.  Framework: for A = kQ/J^2 the reduced cochain complex splits
# into parallel-path terms; the differential dies in char 2 on k[x]/(x^2).
# Verification status: framework fetched-source (via Wang arXiv:1511.08348 and
# Sanchez-Flores arXiv:0711.2810); the loop/Kronecker values derived + hand-verified.
def test_cibils_kx2_char2_doubling():
    """k[x]/(x^2) (one loop, rad^2=0): char 2 => HH^0 = 2 and HH^n = 2 for all n>=1
    (the Cibils differential dies).  A first-class GF(2)-vs-char-0 regression pin."""
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(2))
    assert A.hochschild_cohomology(5, verbose=False).dims == [2, 2, 2, 2, 2, 2]


@pytest.mark.parametrize("fld", [GF(3), GF(32003), CC], ids=["GF3", "GF32003", "CC"])
def test_cibils_kx2_generic_char(fld):
    """k[x]/(x^2), char != 2: HH^0 = 2 (commutative center), HH^n = 1 for n>=1."""
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=fld)
    assert A.hochschild_cohomology(5, verbose=False).dims == [2, 1, 1, 1, 1, 1]


@pytest.mark.parametrize("fld", [GF(2), GF(3), CC], ids=["GF2", "GF3", "CC"])
def test_cibils_kronecker_radsq(fld):
    """The 2-Kronecker path algebra kK_2 (2 vertices, 2 parallel arrows, no
    relations) is ALREADY radical-square-zero (no length-2 paths): HH^0 = 1,
    HH^1 = 3 (= sl_2), HH^{>=2} = 0, in every characteristic.  Euler check
    2 - 4 = -2 (= HH^0 - HH^1).  This is the Cibils-framework instance of the
    m-Kronecker HH^* = [1, m^2-1, 0, 0] value."""
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=fld)
    assert A.hochschild_cohomology(3, verbose=False).dims == [1, 3, 0, 0]


# =========================================================================== #
# 2b.  Taillefer (2001) cyclic homology of Taft / cyclic-Nakayama algebras       #
# =========================================================================== #
# Citation: R. Taillefer, "Cyclic homology of the Taft algebras and of their
# Auslander algebras," arXiv:math/0009214 (cf. K-Theory 24 (2001)); computations
# via Cibils' mixed complex.  Field containing Q (char 0 only).  Thm 2.2 + example,
# Cor 2.8.  Algebra: the Taft algebra Lambda_n = k Z_n / J^n is the self-injective
# cyclic Nakayama algebra (n vertices, Loewy length n), dim n^2.  Verification
# status: fetched-source (ar5iv).
#   HH_0 = k^n,  HH_p = k^{n-1} (p>=1);  HC_{2c} = k^n,  HC_{2c+1} = k^{n-1}.
@pytest.mark.parametrize("n", [2, 3])
def test_taillefer_taft_hochschild_homology(n):
    """HH_*(Lambda_n) = [n, n-1, n-1, ...] over CC (genuine char 0).  Computed via
    the CS engine (monomial; cheap where the bar oracle blows up at dim n^2)."""
    A = NakayamaAlgebra(n=n, l=n, cyclic=True, field=CC)
    assert cs_homology_dims(A, 6).dims == [n] + [n - 1] * 6


def test_taillefer_taft_cyclic_homology_n2_charzero():
    """HC_*(Lambda_2): HC_{2c} = 2, HC_{2c+1} = 1 => [2,1,2,1,2,1] to degree 5 over
    CC (GENUINE char 0), cross-anchored to degree 6 by the char-0 proxy GF(32003)."""
    A_cc = NakayamaAlgebra(n=2, l=2, cyclic=True, field=CC)
    assert A_cc.cyclic_homology(5).dims == [2, 1, 2, 1, 2, 1]
    A_proxy = NakayamaAlgebra(n=2, l=2, cyclic=True, field=GF(32003))
    assert A_proxy.cyclic_homology(6).dims == [2, 1, 2, 1, 2, 1, 2]


def test_taillefer_taft_cyclic_homology_n3_charzero_proxy():
    """HC_*(Lambda_3): HC_{2c} = 3, HC_{2c+1} = 2.  CC is infeasible at dim 9
    (the (b,B) mixed complex exceeds max_cells); pinned over GF(32003) as an
    explicit char-0 proxy (see module CHAR NOTES).  [3, 2, 3] to degree 2 exhibits
    the even/odd (n / n-1) alternation decisively; deeper is OOM-risky."""
    A = NakayamaAlgebra(n=3, l=3, cyclic=True, field=GF(32003))
    assert A.cyclic_homology(2).dims == [3, 2, 3]


# =========================================================================== #
# 2c.  Canonical-algebra HH  (Happel via Schremmer, Prop 4.2.8)                  #
# =========================================================================== #
# Citation: F. Schremmer, "Weighted projective lines and Hochschild cohomology,"
# arXiv:2512.08414 (2025), Prop 4.2.8 / Cor 4.2.9, attributing Happel LNM 1404
# [Hap98 Thm 2.4].  For a canonical algebra C(p_1..p_t), t>=3:
#   dim HH^i = 1 (i=0), 0 (i=1), t-3 (i=2), 0 (i>=3);
#   HH_0 = #vertices = 2 + sum(p_i - 1),  HH_{>=1} = 0 (acyclic).
# Verification status: fetched-source (Schremmer; primary Happel LNM 1404).  The
# formula is char-0-claimed; here CC AND GF(32003) AGREE (no char discrepancy), so
# both are pinned.
def _canonical_equal2(t_arms, weights, field):
    """Ringel canonical algebra C(2,2,...,2) with t arms (all weights p_i = 2):
    source vertex 0, sink vertex 1, one internal vertex per arm; arm i = arrows
    a_i:0->m_i and b_i:m_i->1, arm path pi_i = a_i b_i.  Relations (t-2 of them,
    non-monomial): pi_i - pi_2 + lambda_i pi_1 = 0 for i = 3..t, with distinct
    nonzero weights lambda_i != 1.  #vertices = 2 + t."""
    verts = [0, 1] + [10 + i for i in range(1, t_arms + 1)]
    arrows = {}
    for i in range(1, t_arms + 1):
        mi = 10 + i
        arrows[f"a{i}"] = (0, mi)
        arrows[f"b{i}"] = (mi, 1)
    rels = []
    for idx, i in enumerate(range(3, t_arms + 1)):
        lam = weights[idx]
        rels.append(f"a{i}*b{i} - a2*b2 + {lam}*a1*b1")
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


@pytest.mark.parametrize("fld", [CC, GF(32003)], ids=["CC", "GF32003"])
def test_canonical_C222_t3(fld):
    """C(2,2,2) (t=3): HH^* = [1,0,0,0,...] (HH^2 = t-3 = 0); HH_0 = 5 = #vertices,
    HH_{>=1} = 0.  Non-monomial admissible => engine = CS."""
    A = _canonical_equal2(3, [2], fld)
    assert cs_cohomology_dims(A, 4).dims == [1, 0, 0, 0, 0]
    assert cs_homology_dims(A, 4).dims == [5, 0, 0, 0, 0]


@pytest.mark.parametrize("fld", [CC, GF(32003)], ids=["CC", "GF32003"])
def test_canonical_C2222_t4(fld):
    """C(2,2,2,2) (t=4): HH^0 = 1, HH^1 = 0, HH^2 = t-3 = 1, HH^{>=3} = 0; HH_0 = 6
    = #vertices, HH_{>=1} = 0.  Non-monomial admissible => engine = CS."""
    A = _canonical_equal2(4, [2, 3], fld)
    assert cs_cohomology_dims(A, 4).dims == [1, 0, 1, 0, 0]
    assert cs_homology_dims(A, 4).dims == [6, 0, 0, 0, 0]
