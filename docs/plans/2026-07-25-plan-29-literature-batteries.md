# Plan 29 — Literature-oracle batteries + the is_symmetric fix (Tier 1a)

**Date:** 2026-07-25. **Branch:** `plan-29-literature-batteries`. **Backlog:**
Tier-1a "Literature-oracle battery expansion (Marco, 2026-07-25 deep research)".
**Evidence base:** `2026-07-25-literature-oracles-deep-research.md` (the four
cluster reports; verification-status labels there are BINDING — nothing flagged
abstract-only / lossy / figure-unverified gets frozen as a strict pin).

## Part 0 — DO FIRST: the `is_symmetric` bug (QPA-verified, live)

`invariants/frobenius.py::is_symmetric_generic` returns **False** on provably
symmetric algebras:

- multi-vertex symmetric Nakayama `kZ_n/J^L` with `n | (L−1)`, `n ≥ 2`
  (Brauer stars: kZ₂/J³, kZ₃/J⁴, kZ₄/J⁵ — QPA `IsSymmetricAlgebra = true`,
  we say False). Weakly-symmetric detection (Nakayama permutation = id) is
  CORRECT; the downstream "ν is inner" search fails only for multi-vertex.
  Root-cause properly (systematic debugging — understand WHY the inner-ν
  sweep misses; do not patch symptoms). Single-vertex and genuinely
  non-symmetric cases are correct today and must stay correct.
- `TrivialExtension(A)` for every A (classically symmetric — A ≅ D(A) gives
  the symmetrizing form). Root cause differs: T(A) carries no quiver
  presentation, so the path-basis route can't run. Fix EITHER by a
  double-quiver presentation for `TrivialExtension` OR by making the symmetry
  certifier presentation-free (a symmetrizing linear form λ with λ(ab)=λ(ba)
  and no-radical-kernel is checkable from structure constants alone) —
  implementer chooses after diagnosis, documents why.

Regression battery: `tests/invariants/` (fast) — Brauer stars n=2,3,4
symmetric; kZ₃/J³ NOT symmetric (n∤L−1); `is_weakly_symmetric` one-liner added
and tested; `TrivialExtension(kA₂/kA₃/Kronecker/comm-square).is_symmetric()`
True; existing single-vertex pins unchanged. `tests/qpa/` (qpa) —
`crosscheck_symmetric` vs `IsSymmetricAlgebra`/`IsWeaklySymmetricAlgebra`
across the zoo incl. the Brauer stars. Citation: Skowroński–Yamagata
(*Frobenius Algebras I*, EMS 2011) for the n | (L−1) criterion.

## Part 1 — Coxeter / spectral batteries (fast; invariants exist, tests only)

Source: de la Peña/Lenzing/Marcos cluster report (convention MAPPED:
`charpoly(−C^{−T}C) ≡ charpoly(−C^{−1}C^t)` — our `coxeter_polynomial` equals
the papers' exactly; six Nakayama polynomials already recomputed and matched).

1. **Nakayama `N_n(r)` (LMR arXiv:2112.15587 §2.7 + Prop 6.1):** the six exact
   polynomials χ(17,8), χ(16,3), χ(15,6)=χ(15,4), χ(15,5), χ(14,7) + the
   `r ≥ 9 ⇒ χ(r+7,r) = (λ+1)(λ⁶−λ³+1)(λ^r+1)` family (2 instances) + Coxeter
   number = min m with Φ^m = I where stated.
2. **Dynkin table (LdlP 2008, arXiv:0805.1018):** Aₙ = v_{n+1} (n = 2..6),
   Dₙ = Φ₂·(x^{n−1}+1) via the v-factorization for n = 4,5,6 (NOT the flagged
   cyclotomic-condition column), E₆ = Φ₃Φ₁₂, E₇ = Φ₂Φ₁₈, E₈ = Φ₃₀ — and
   orientation-independence (two orientations per type agree).
3. **Affine + canonical (triple-corroborated):** Ã_{p,q} = (x−1)²v_p v_q (two
   (p,q) points, orientation-DEPENDENT — test that too), D̃₄, Ẽ₆; canonical
   `χ = (x−1)²∏v_{p_i}` for (2,2,2), (2,3,5), (3,3,3) built as
   quiver-with-relations.
4. **Lehmer [2,3,7] (dlP arXiv:1310.1910 §2.5):** the wild star's
   `coxeter_polynomial` == Lehmer's polynomial exactly;
   `mahler_measure` = `spectral_radius` = the Salem root; the strict ordering
   ρ[2,4,5] > ρ[2,3,8] > ρ[2,3,7]; plus `mahler_measure == 1` ⇔
   cyclotomic-type on the Dynkin/affine battery members.
5. **B₃ parametric (LdlP 2008):** Cartan [[1,a,b],[0,1,c],[0,0,1]] ⇒
   χ = x³+αx²+αx+1, α = abc−a²−b²−c²+3 — via a raw structure-constant/Cartan
   entry point if buildable, else the (1,1,1) instance as quiver algebra.
6. **Happel Prop 1.5:** χ_A(−1) is a perfect square — asserted across the
   whole battery (cheap sanity sweep).

## Part 2 — Identity oracles (fast where cheap, deep where HH is deep)

1. **Happel-1997 trace identity** (LAA 258; sign PINNED on A₃ first:
   `tr(coxeter_matrix) == −Σ(−1)^i dim HH^i`): hereditary set (2/3-Kronecker,
   A₃, A₄, D₄, acyclic triangle — top=2 suffices) + two finite-gl.dim
   non-hereditary members (comm square; a monomial gl.dim-2 zoo member) with
   honest top = gl.dim.
2. **Derived invariance across orientations** (Keller/Rickard): the two
   acyclic Ã₂-triangle orientations give equal HH^•/HH_•/HC_• (pinned values
   [1,2,0]/[3,0,0]/[3,0,3] from the report, re-verified live); one Dynkin
   graph with two orientations (trivial [1,0,…] but pins the scheme).
3. **Hereditary HH¹ formula** (Happel 1989 via arXiv:2509.05135 Cor 4.10):
   m-Kronecker HH^• = [1, m²−1, 0, 0] for m = 1,2,3 (= a-Kronecker toupie pin,
   Euler-verified); trees → [1,0,0]; acyclic triangle → [1,2,0]; a-Kronecker
   over CC AND GF(p) (char-independent).
4. **Trivial extension HH¹** (CMRS 2003 + CRS 2004): `HH¹(T(A)) ≠ 0` for every
   zoo/battery member (structural, any char); Example 2.20 (Z₅ cycle, two
   length-3 relations — TRANSLATE the paper's composition order to ours and
   note it) ⇒ HH¹(A) = 0 over char 0.
5. **Incidence ≅ nerve** (Cibils 1989 / Redondo 2008): crown/diamond poset
   (nerve ≃ S¹) ⇒ HH^• = [1,1,0,0]; a contractible-nerve poset ⇒ HH^{≥1} = 0.
   Label: derived-from-theorem.
6. **Acyclic vanishing + truncated finiteness boolean** (Cibils 1986;
   Xu–Han–Jiang Thm 3): acyclic members ⇒ HH_{≥1} = 0; truncated kQ/R^N:
   finite total HH^• ⇔ Q acyclic — assert on one acyclic + one cyclic
   truncated pair (boolean only; the lossy c_{N,e,i} counts are NOT pinned).

## Part 3 — HH/HC value batteries (deep bucket; CS/Bardzell/cyclic engines)

1. **Bergh–Erdmann QCI, general (a,b)** (AN&T 2 (2008) Thms 3.1/3.2;
   `BerghErdmann2008` exists — extend): over CC ONLY (q = 2 not a root of
   unity; the GF(p) trap is documented in the research doc and in the test
   docstring): cohomology `[2,2,1,0,…]` and homology `[a+b−1, a+b−2, …]` for
   (a,b) ∈ {(2,3),(3,3)} to degree 6 (a=b=2 is already pinned). `engine="cs"`.
2. **Triangular string Aₙ** (Redondo–Román JPAA 218 (2014), Example 3;
   monomial ⇒ Bardzell + CS; char-independent — GF(2), GF(3), CC):
   A₁ [1,3,0,0], A₂ [1,2,0,0], A₃ [1,3,0,2,0], A₄ [1,4,0,0,0] — **bar
   cross-check the A₃ revival value in the test itself** (agent
   recommendation) before pinning.
3. **Rad²=0 char discriminator** (Cibils 1998 framework, hand-verified):
   `k[x]/(x²)`: GF(2) ⇒ HH^n = 2 ∀n≥1 vs char ≠ 2 ⇒ HH^n = 1 (existing pins
   may partially cover — extend to the GF(2) side explicitly); Kronecker
   HH^• = [1,3,0,0] any char (= Part 2.3, cross-file OK to share one pin).
4. **Taft/Nakayama cyclic homology** (Taillefer arXiv:math/0009214, char 0):
   Λₙ = cyclic Nakayama kZ_n/J^n for n = 2,3: HH_• = [n, n−1, n−1, …],
   HC_{2c} = n, HC_{2c+1} = n−1, to degree 6.
5. **Canonical algebra HH** (Happel via Schremmer arXiv:2512.08414 Prop
   4.2.8): (2,2,2) t=3 ⇒ HH^• = [1,0,0,0,…]; (2,2,2,2) t=4 ⇒ HH² = 1;
   HH_0 = #vertices, HH_{≥1} = 0 — and cross-link: the Part-2.1 trace
   identity on the same algebras (two-invariant consistency).
6. **Cup predicates** (exercises Plans 20/21): triangular-string cup
   triviality (RR 2014 Thm 5.2 — HH^n ∪ HH^m = 0 for n,m > 0 on A₁/A₂);
   the RR-2018 quadratic-string Example 3.1 cup ≠ 0 ONLY IF the bar
   cross-check pins a nonzero product cleanly — otherwise defer with a note
   (the paper's integer vectors are convention-risky per the report).

## Part 4 — module Tor_n^A(M, N) (Marco, 2026-07-25 mid-plan request)

Ext^n(M,N) for arbitrary same-side modules EXISTS (`modules/ext.py`, surfaced
incl. no-code `ext_target`). Tor does not — add it:
`src/quiverlab/modules/tor.py::tor_dims(A, M, N, top)` for a RIGHT module M
and LEFT module N (Plan-24 sides; loud refusal on wrong sides/cross-algebra,
mirroring `hom.py::_assert_comparable`): H_n(P_• ⊗_A N) with P_• the existing
minimal projective resolution of M; `P_i ⊗_A N` collapses summand-wise
(`e_v A ⊗_A N ≅ e_v N`) — exact linear algebra over any Domain. Delegator
`Algebra.tor(M, N, n)` + `tor_dims` surface, Plan-23/24 style.
**Oracles (self-certifying + literature):** (i) the classical duality
`dim Tor_n^A(M, N) = dim Ext_A^n(M, DN)` (D side-aware, Plan 24) asserted
degreewise across the battery — every Tor value cross-checked against the
EXISTING Ext engine; (ii) resolve-M vs resolve-N agreement
(Tor_n(M,N) via P_•(M) ⊗ N ≡ M ⊗ P_•(N) computed through A^op); (iii)
hereditary vanishing Tor_{≥2} = 0; semisimple Tor_{>0} = 0; simples over
kA_n; `Tor_0 = M ⊗_A N` dim checked by hand on kA₂. Cite Cartan–Eilenberg
(existing `tensor_product` key) + assem_book. QPA: probe `TorOverAlgebra`
live in the qpa test file (session-side, do NOT edit crosscheck.py — Part 0
owns it); if absent, the honest-scope note names oracle (i) as coverage.
GUI/webapp `tor` compute kind = follow-up backlog note, NOT this plan.

## Citations (ONE agent owns registry/bib edits — no merge races)

New keys (BibTeX in the research doc): `redondo_roman_2014`,
`redondo_roman_2018` (only if 3.6 lands), `toupie` (ArtensteinLanzilottaSolotar2020),
`ladkani_gentle` (only if used), `cibils_radsq` (Cibils1998),
`taillefer_taft`, `cibils_incidence` (Cibils1989 + Redondo2008),
`cmrs_split` (CMRS2003), `crs_trivial_ext_hh1` (CRS2004), `xhj_truncated`
(XuHanJiang2007), `cibils_acyclic` (Cibils1986), `happel_trace` (Happel1997),
`keller_cyclic_invariance` (Keller1998), `rickard_morita`/`rickard_stable`
(one key suffices for the derived-invariance scheme), `lenzing_meltzer_ruan`,
`lenzing_delapena_spectral`, `delapena_mahler`, `skowronski_yamagata`
(Frobenius Algebras I), `schremmer_wpl` OR cite Happel LNM 1404 via the
existing `happel_question` key where the repo already records it (prefer
existing keys; add only what a test actually cites — the standing
no-guessed-numbers rule applies).

## NOT pinned (honest-scope, listed on the verification page)

Artenstein Ex. 7.4.1 (figure unverified); Xu–Han–Jiang c_{N,e,i} counts
(lossy transcription); Cibils 1990 rigidity count (walled); Barot–de la Peña
D₁₁ (figure); Ladkani φ_Λ hand-computations (prefer the RR integer vectors);
Bergh–Erdmann char-p homology branches (need infinite char-p fields);
RR-2018 exact integer vectors (unless bar-anchored in-test).

## Acceptance

Part-0 fix + regression battery green incl. QPA; every battery cites its
source in the test docstring (registry key where added); `docs/verification.md`
gains the new pins (subsystem rows + Named-literature-pins section +
References) and the not-pinned honest-scope entries; counts resync at merge;
backlog tick; fast/deep/qpa suites green. Merge/push only when Marco asks.
