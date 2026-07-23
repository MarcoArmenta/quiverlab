"""Validation battery (d): the READ-ONLY bank's closed-form CS families as
byte-level oracles (spec-§6; Plan-04 Task 11).

Batteries (a) [CS == Bardzell] and (b) [CS == bar] compare two of the library's OWN
engines; battery (c) pins CS HH against classical/literature dim vectors.  Battery (d)
pins the Plan-04 CS resolution against the ORIGINAL HanLab bank resolution
`hanlab/resolutions_cs.py` -- a wholly separate implementation of the Chouhy-Solotar
differentials (the closed forms of arXiv:1406.2300 for the two validated families
k[x]/(x^a) and the quantum CI k<x,y>/(x^2,y^2,yx-xi xy)).  The bank is READ-ONLY
(project law): these tests IMPORT it by path and never regenerate or modify it.

Two independent CS implementations, cross-checked:
  * quiverlab: CSResolution(A) over GF(p) -- the Plan-04 leading map delta_n PLUS the
    order-condition correction gamma (resolutions_cs/resolution.py).
  * bank: ChouhySolotarResolution(rs, alg) -- the hand-derived closed-form CS
    differentials (four-parity-case formulae), raw int64, never reduced mod p.

TWO tiers of comparison:
  1. BINDING -- HH_* dimension vectors must agree, family by family, prime by prime.
     Invariant under the correction's nullspace non-uniqueness (a permutation or a
     nullspace shift of the differential cannot change its rank), yet swap-sensitive
     to any genuine Task-6 differential bug (which would move a rank).  These MUST pass.
  2. STRICT (since Plan 17) -- the collapsed differential matrices agree
     ENTRY-BY-ENTRY mod p, once the two generator orders are aligned (quiverlab
     Chain.word <-> bank chain label ("c",)/("v",) for truncpoly, (s,t) for qci;
     corner index is identity -- both libraries order the A-basis as
     [1, x, y, xy] / [1, x, x^2, ...]).

===========================================================================
HISTORY (the pre-Plan-17 EMPIRICAL FINDING, kept for the record):

  *** The byte-exact pins PASS byte-for-byte -- they XPASS, not xfail. ***

  The brief marked the entrywise pins `xfail(strict=False, reason="canonicalization
  pending")`, expecting the Plan-04 order-condition correction gamma (which is unique
  only MOD THE NULLSPACE for a general admissible algebra) to make the collapsed
  matrices differ from the bank's hand-derived closed forms.  Measured on THIS branch,
  it does not: for BOTH families, EVERY degree n = 1..N, and EVERY prime in
  {32003, 2, 3, 5}, the aligned quiverlab and bank differential matrices are equal
  entry-by-entry mod p (0 mismatches).  But the MECHANISM is NOT the same for the two
  families, and the naive "gamma is trivially pinned everywhere" story is FALSE:

    * truncpoly / monomial family k[x]/(x^a) -- byte identity is FORCED.  The
      order-condition correction gamma is genuinely 0 here (the lower-generator solve
      is empty/degenerate on a monomial algebra), so the Plan-04 leading map delta_n
      already IS the canonical CS closed form.  Nothing to choose; no nullspace.

    * quantum-CI family k<x,y>/(x^2,y^2,yx - xi xy) -- byte identity is a SOLVER-CHOICE
      COINCIDENCE, not uniqueness.  Here the correction system has genuine nullspace
      freedom: its nullity GROWS WITH DEGREE (measured up to 6), gamma is NONZERO, and
      DISTINCT valid solutions yield DISTINCT byte matrices that ALL still satisfy the
      order condition / d^2 = 0 locally.  The pins pass only because quiverlab's
      deterministic particular solution -- RREF with free variables set to 0
      (fields/linalg.py `solve`) -- happens to land on the SAME representative as the
      bank's minimal hand-derived closed form.  This is a canonical-representative
      COINCIDENCE (two independent tie-breakers agreeing), NOT a mathematical
      guarantee of uniqueness.

  RESOLVED BY PLAN 17 (2026-07-23): the tie-breaking coincidence became a
  guarantee.  `_d_general` now reduces gamma through
  `fields.linalg.reduce_mod_nullspace` -- the unique coset representative with
  zero coordinates at every free column of the correction system's RREF, i.e.
  exactly the representative both the solver and the bank's hand derivation were
  already choosing.  Byte identity is therefore pinned BY CONSTRUCTION
  (adversarial-solver gate: tests/resolutions_cs/test_canonicalization.py shifts
  the solve by a nullspace vector and the differential bytes must not move), and
  the entrywise pins below are STRICT plain tests.  Tier 1 (HH-dim equality)
  remains the rank-based, nullspace-invariant load-bearing guarantee.
===========================================================================
"""
import functools
import importlib.util
import pathlib
import sys

import numpy as np
import pytest

from quiverlab import GF, Quiver
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.hh_engine import hochschild_homology_dims
from quiverlab.families import truncated_polynomial
from quiverlab.resolutions_cs.engine_facade import CSResolution

pytest.importorskip("quiverlab.groebner")

BANK = pathlib.Path(
    "/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/HansConjecture")

if not (BANK / "hanlab" / "resolutions_cs.py").exists():
    pytest.skip("bank HansConjecture/hanlab not present", allow_module_level=True)


@functools.lru_cache(maxsize=1)
def _bank():
    """(hanlab module, bank resolutions_cs module).  READ-ONLY bank, imported by
    path (never regenerated, never modified).  `hanlab` supplies the bank Algebra
    builders (truncated_polynomial, quantum_ci) and homology_dims; the closed-form
    `ChouhySolotarResolution` + reduction-system factories are loaded straight from
    the bank file by spec (per the brief's `_bank_cs`)."""
    added = [p for p in (str(BANK), str(BANK / "hanlab")) if p not in sys.path]
    for p in added:
        sys.path.insert(0, p)
    try:
        import hanlab as H
        spec = importlib.util.spec_from_file_location(
            "bank_cs", BANK / "hanlab/resolutions_cs.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
    finally:
        # Pop only the path entries WE added; the loaded modules stay in
        # sys.modules (needed by callers) -- only the sys.path pollution is undone.
        for p in added:
            try:
                sys.path.remove(p)
            except ValueError:
                pass
    return H, m


# Primes: the bank's default battery -- a large char-0 proxy plus the small primes
# that probe characteristic pathology (p | a collapses the truncpoly norm map;
# p | xi^2 - 1 etc. reshape the qci homology).
PRIMES = (32003, 2, 3, 5)


# --- quiverlab-side algebra builders (over GF(p)) --------------------------------
def _ql_truncpoly(a, p):
    return truncated_polynomial(a, field=GF(p))


def _ql_qci(xi, p):
    # bare coefficient: the Plan-03 relation grammar rejects "y*x - (2)*x*y".
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", f"y*x - {xi}*x*y"], field=GF(p))


# --- generator-order alignment: quiverlab (Chain, corner) -> bank ((label), corner)
def _key_truncpoly(g):
    ch, j = g
    return (("v",) if ch.degree == 0 else ("c",), j)


def _key_qci(g):
    ch, j = g
    w = ch.word
    return ((w.count("y"), w.count("x")), j)


# --- tier-1 helper: HH-dim equality (binding, swap-sensitive) --------------------
def _assert_hh_matches_bank(ql_build, bank_alg, bank_rs, N):
    """HH_*(bank_alg) through the bank CS resolution vs the Plan-04 CSResolution must
    agree as full dim vectors, prime by prime.  The bank build is prime-independent
    (raw int64) so one resolution serves every prime; the quiverlab facade computes
    over GF(p), so A is rebuilt per prime and its ranks taken at that same p."""
    H, BCS = _bank()
    res = BCS.ChouhySolotarResolution(bank_rs, bank_alg)
    bank = H.homology_dims(bank_alg, N, primes=PRIMES, resolution=res)
    for p in PRIMES:
        A = ql_build(p)
        E = to_engine(A.unit_adapted())
        ql = hochschild_homology_dims(
            E, N, primes=(p,), resolution=CSResolution(A))[p]
        assert list(ql) == list(bank[p]), (
            f"HH-dim disagreement p={p}: quiverlab {list(ql)} != bank {list(bank[p])}")


# --- tier-2 helper: entrywise-mod-p equality of the collapsed differentials ------
def assert_matrices_equal_mod_p(ql_build, bank_alg, bank_rs, N, p, keyfn):
    """Entry-by-entry equality mod p of the collapsed CS homology differentials
    d_n : C_n -> C_{n-1}, n = 1..N, with the two generator orders aligned.  Reads the
    bank `ChouhySolotarResolution.differential_matrix` (raw int64) and reindexes the
    quiverlab CSResolution matrix into the bank generator order via `keyfn` (quiverlab
    (Chain, corner) -> bank ((label), corner)), then asserts congruence mod p."""
    H, BCS = _bank()
    A = ql_build(p)
    E = to_engine(A.unit_adapted())
    cs = CSResolution(A)
    res = BCS.ChouhySolotarResolution(bank_rs, bank_alg)
    for n in range(1, N + 1):
        ql_n = cs.term_basis(E, n)
        ql_nm1 = cs.term_basis(E, n - 1)
        Mq = cs.differential_matrix(
            E, n, ql_n, {g: i for i, g in enumerate(ql_nm1)})
        bank_n = res.term_basis(bank_alg, n)
        bank_nm1 = res.term_basis(bank_alg, n - 1)
        Mb = res.differential_matrix(
            bank_alg, n, bank_n, {g: i for i, g in enumerate(bank_nm1)})
        # position of each bank generator within the quiverlab ordering
        pos_n = {keyfn(g): i for i, g in enumerate(ql_n)}
        pos_nm1 = {keyfn(g): i for i, g in enumerate(ql_nm1)}
        cols = [pos_n[bg] for bg in bank_n]
        rows = [pos_nm1[bg] for bg in bank_nm1]
        Mq_aligned = Mq[np.ix_(rows, cols)]
        assert Mq_aligned.shape == Mb.shape, (
            f"shape mismatch n={n} p={p}: {Mq_aligned.shape} vs {Mb.shape}")
        assert ((Mq_aligned - Mb) % p == 0).all(), (
            f"byte disagreement n={n} p={p} (aligned generator order)")


# ================================ tier 1 (binding) ===============================
@pytest.mark.parametrize("a", [2, 3, 4])
def test_cs_hh_matches_bank_truncpoly(a):
    # k[x]/(x^a): HH_* dims via Plan-04 CSResolution == via bank ChouhySolotarResolution.
    H, BCS = _bank()
    _assert_hh_matches_bank(
        lambda p: _ql_truncpoly(a, p),
        H.truncated_polynomial(a),
        BCS.truncpoly_reduction_system(a), N=6)


@pytest.mark.parametrize("xi", [1, 2, 3])
def test_cs_hh_matches_bank_quantum_ci(xi):
    # k<x,y>/(x^2,y^2,yx - xi xy): same, HH_* dims equal (both resolutions), all primes.
    H, BCS = _bank()
    _assert_hh_matches_bank(
        lambda p: _ql_qci(xi, p),
        H.quantum_ci(xi),
        BCS.qci_reduction_system(xi), N=8)


# ========================= tier 2 (aspirational byte pins) =======================
@pytest.mark.parametrize("a", [2, 3, 4])
def test_cs_matches_bank_truncpoly_bytes(a):
    H, BCS = _bank()
    for p in PRIMES:
        assert_matrices_equal_mod_p(
            lambda pp: _ql_truncpoly(a, pp),
            H.truncated_polynomial(a),
            BCS.truncpoly_reduction_system(a), N=6, p=p, keyfn=_key_truncpoly)


@pytest.mark.parametrize("xi", [1, 2, 3])
def test_cs_matches_bank_quantum_ci_bytes(xi):
    H, BCS = _bank()
    for p in PRIMES:
        assert_matrices_equal_mod_p(
            lambda pp: _ql_qci(xi, pp),
            H.quantum_ci(xi),
            BCS.qci_reduction_system(xi), N=8, p=p, keyfn=_key_qci)
