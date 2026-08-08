# Plan 54 — The BV operator Δ on HH of Frobenius / self-injective algebras (R2)

> **For agentic workers:** REQUIRED SUB-SKILL — use `superpowers:subagent-driven-development`
> (or `superpowers:executing-plans`) to implement this plan task-by-task. Steps
> use checkbox (`- [ ]`) syntax for tracking. TDD is mandatory: write the failing
> test first, watch it fail for the stated reason, implement, watch it pass.

**Date:** 2026-08-07 · **Branch:** `plan-54-bv-operator` (off `dev`) ·
**Status:** spec — awaiting P52 alignment then implementation ·
**Wave:** 1 (tier α, flagship engine upgrades) · **Size:** M ·
**Prerequisite:** **P52** (HH with bimodule coefficients —
`docs/plans/2026-08-07-plan-52-hh-coefficients.md`, being written in parallel;
NOT yet on disk at authoring time — see §3 *P52 alignment* for the exact
consumption contract and the fallback if P52's API lands differently).

**Change-log.** *2026-08-07 adversarial review: 4 majors + 9 minors applied
(BIKLZ pins UN-BLOCKED via the arXiv HTML full text; twist-direction honesty —
the arbiter determines direction, the dim-match is demoted to a coarse
coefficient-sanity check, and direction validation is P52-gated; GF(2) SERVED
with a char-2 certification caveat; a P52 chain-level convention contract added
as a cross-plan alignment block).*

---

## Record (R2, verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R2 — BV operator Δ on HH of Frobenius/self-injective algebras.** Object:
> Δ: HH^n → HH^{n-1} by transporting Connes B through the σ-twisted Frobenius
> duality (design sketch Δ = ∂∘B_σ∘∂⁻¹ — a paraphrase, not a verbatim formula);
> bracket recovered from Δ (7-term BV axioms self-certified against the
> independent bracket). Corrected refs: Lambre–Zhou–Zimmermann, arXiv:1405.5325,
> *J. Algebra* 2016 (hypothesis: SEMISIMPLE Nakayama automorphism); Volkov,
> arXiv:1405.5155 (hypothesis: ord(ν) coprime to char k — a DIFFERENT computable
> condition; implementation must pick deliberately);
> Bian–Itagaki–Kou–Lyu–Zhou arXiv:2603.04834 (2026; self-injective Nakayama,
> semisimplicity hypothesis removed — the explicit kZ₁/J^N Δ/bracket pins MUST be
> transcribed verbatim from the paper body with equation numbers before any test
> uses them; the scout's quoted formulas are unverified). Oracles: Δ²=0,
> bracket-from-Δ ≡ independent bracket, symmetric case = Tradler on
> k[x,y]/(x²,y²) (HH_• = [4,4,5,6] pin exists), QuantumCI ties. Size M. Deps:
> connes_differentials, nakayama_automorphism, is_frobenius, cup/bracket.

---

## 0. Reference re-verification (mandatory, standing rule §1.5 of the metaplan)

Done at spec time via WebSearch/WebFetch on 2026-08-07. Recorded exactly:

| Record ref | Verified | Finding |
|---|---|---|
| arXiv:**1405.5325** | ✅ verified | Lambre, Zhou, Zimmermann, *"The Hochschild cohomology ring of a Frobenius algebra with **semisimple Nakayama automorphism** is a Batalin–Vilkovisky algebra"*, J. Algebra 446 (2016). Hypothesis **confirmed: semisimple ν**. The paper *also supplies a decidable criterion* for when a Frobenius quiver-with-relations algebra has semisimple ν — relevant to our hypothesis gate. |
| arXiv:**1405.5155** | ✅ verified | Volkov, *"BV-differential on Hochschild cohomology of Frobenius algebras"*. Hypothesis **confirmed: ord(ν) NOT divisible by char k** (= coprime to char over a field). Volkov constructs an algebra `HH*(R)^{ν↑}` that is ALWAYS a BV algebra, and is isomorphic to `HH*(R)` exactly when `char k ∤ ord(ν)`. The scout's *swap* of 5325↔5155 is confirmed corrected. |
| arXiv:**2603.04834** (2026) | ✅ verified — **full text available, formulas transcribed** | Bian, Itagaki, Kou, Lyu, Zhou, *"The Hochschild cohomology ring of a self-injective Nakayama algebra is a Batalin–Vilkovisky algebra"*. Title, authorship, and **hypothesis (self-injective Nakayama, semisimplicity removed)** verified. **A clean full-text HTML exists at `https://arxiv.org/html/2603.04834`; §3.2 ("The case e = 1 for truncated basic cycle algebras Λ = 𝕂Z_e/J^N") gives the kZ₁/J^N BV operator and bracket EXPLICITLY.** The dimension summary is now sourced directly from §3.2 (see the transcribed draft below), and the classical truncated-polynomial computation agrees. The earlier "inaccessible / lossy" assessment was WRONG — the HTML renders cleanly and is the transcription target. **The formulas are usable; the implementer re-verifies the transcription below against the HTML at implementation time (this plan does not claim the text is unavailable).** |
| Tradler AIF 2008 | ✅ verified | Tradler, *"The Batalin–Vilkovisky algebra on Hochschild cohomology induced by infinity inner products"*, Ann. Inst. Fourier 58 (2008) no. 7, 2351–2379, DOI 10.5802/aif.2417. BV structure on `HH*` from a **symmetric, invariant, non-degenerate** inner product — the symmetric-algebra anchor (ν inner). |

**BIKLZ §3.2 (e = 1) — critic-transcribed draft (re-verify against the HTML at
implementation time).** For `Λ = 𝕂Z₁/Jᴺ = k[x]/(xᴺ)`, `N ≥ 2`, arXiv:2603.04834
§3.2 records:

* **`char 𝕂 ∤ N`.** Basis generators `x₀, y, z` in degrees `0, 1, 2`
  (relations `x₀ᴺ = 0`, `y·x₀ᴺ⁻¹ = 0`, `x₀ᴺ⁻¹·z = 0`, `y² = 0`).
  `dim HH⁰ = N`, `dim HHⁿ = N − 1` for `n ≥ 1`.
  * **BV:**  `Δ(y·x₀ᵃ·zᵇ) = (bN + N − a − 1)·x₀ᵃ·zᵇ`, for `0 ≤ a ≤ N−2`, `b ≥ 0`.
  * **Bracket:**  `[x₀, y] = −x₀`,  `[y, z] = −N·z`.
* **`char 𝕂 | N`.** Basis generators `x₀, y, w, z'` in degrees `0, 1, 1, 2`
  (`x₀ᴺ = 0`, `y = w·x₀`, `w² = −(N(N−1)/2)·x₀ᴺ⁻²·z'`).
  `dim HH⁰ = N`, `dim HHⁿ = N` for `n ≥ 1`.
  * **BV:**  `Δ(w·z'ᵇ) = 0`  and  `Δ(w·x₀ᵃ·z'ᵇ) = −a·x₀ᵃ⁻¹·z'ᵇ`, for `1 ≤ a ≤ N−1`, `b ≥ 0`.
  * **Bracket:**  `[x₀, w] = −1`,  `[y, w] = −w`.

**Un-blocked.** kZ₁/J^N **= k[x]/(x^N) is symmetric** (commutative Frobenius ⇒
ν = id), so our engine computes its Δ on the **symmetric/Tradler route in every
characteristic**, and the §3.2 values above are the exact-VALUE literature
oracle for that route — a **live `oracle_literature` assert**
(`test_bv_biklz_kz1.py`, §7 / Task F), NOT a fence. The transcription is
re-verified against the HTML during implementation (equation/label identifiers
copied verbatim into the test docstring); we NEVER fabricate a value the HTML
does not print. What remains out of v1 scope is the *general* BIKLZ construction
for `e ≥ 2` self-injective Nakayama with **non-semisimple** ν (§2.4 gate 3) —
that is a different, harder deliverable than validating the e = 1 values, and it
stays a loud refusal.

---

## 1. Motivation and the shape of the deliverable

Plan 35 surfaced the whole Tamarkin–Tsygan calculus as public product objects:
`cup_products`, `cap_products`, `gerstenhaber_brackets`, `connes_differentials`,
each a frozen structure-constant container with one `.blocks()` serialization
(`hochschild/products.py`: `HHProducts`, `ConnesB`, `ProductTable`), the
**basis-dependence doctrine** (constants record WHICH basis; cross-engine gates
compare only basis-independent data), and the induced Connes `B: HH_n → HH_{n+1}`
with `B² = 0` at the induced level.

What is missing — and what a representation theorist working with a symmetric or
self-injective algebra reaches for — is the **BV operator** `Δ: HHⁿ → HHⁿ⁻¹`, the
degree `−1` square-zero operator that makes `HH•` of a Frobenius algebra a
Batalin–Vilkovisky algebra and whose *defect from being a cup-derivation is the
Gerstenhaber bracket*. It is the last structure of the calculus quiverlab does
not expose, and it is the natural companion of `connes_differentials` (its
Frobenius-dual mirror).

The mathematics is a **transport**, not a new engine: `Δ` is Connes `B` carried
across the Frobenius duality `HHⁿ(A) ≅ D(HH_n(A, {}_1A_ν))`. So P54 is a thin
exact-linear-algebra layer over three things quiverlab already has
(`connes_differentials`, `nakayama_automorphism`, `is_frobenius`/`is_symmetric`,
`gerstenhaber_brackets`) plus one thing P52 adds (twisted-coefficient homology).
Every serving path (GUI/webapp/report) already exists for the Plan-35 kinds; P54
adds one more kind of exactly the same shape.

Two non-negotiable house principles drive every decision below:

* **Exact only, loud typed refusals.** The BV structure exists only under a
  hypothesis on ν; we implement the **decidable** hypotheses, route per-instance
  to whichever holds, record which in provenance, and refuse loudly when none
  applies. No float ever enters `src/`.
* **The self-certifications ARBITRATE the conventions.** The twist direction
  (ν vs ν⁻¹), the pairing side, and the BV-relation signs are NOT assumed from a
  paper we cannot fully verify — they are **pinned per instance** by `Δ² = 0` and
  by *the derived bracket equalling the independent Plan-35 bracket in-window*.
  This is the Plan-21 cap sign-arbitration precedent applied verbatim.

---

## 2. The mathematics (design-level, with degree bookkeeping)

### 2.1 The Frobenius duality and the transport

Let `A` be a finite-dimensional Frobenius `k`-algebra with Frobenius form `λ` and
Nakayama automorphism `ν` (the unique-up-to-inner automorphism with
`λ(ab) = λ(b ν(a))`; `nakayama_automorphism()` returns a concrete matrix
`N = G⁻¹Gᵀ`). The Frobenius form induces a **perfect pairing in each degree**

```
⟨ − , − ⟩_n :  HHⁿ(A)  ⊗  HH_n(A, {}_1A_ν)  ⟶  k              (†)
```

between ordinary Hochschild cohomology and the **ν-twisted** Hochschild homology
(coefficients in the bimodule `{}_1A_ν` = `A` with the right action twisted by ν).
At chain level, for a cochain `f: Ā^⊗n → A` and a twisted chain
`z = a_0 ⊗ a_1 ⊗ ⋯ ⊗ a_n`,

```
⟨ f , z ⟩  =  λ( a_0 · f(a_1 ⊗ ⋯ ⊗ a_n) ).
```

The twist by ν on `{}_1A_ν` is *exactly* the correction that makes `⟨δf, z⟩ =
±⟨f, ∂z⟩` (the ν-twisted last face `d_n(z) = ν(a_n)a_0 ⊗ ⋯` uses
`λ(ab)=λ(bν(a))`), so `(†)` descends to (co)homology and is perfect — this is why
**untwisted** homology only works when `ν` is inner. This is the classical
duality of Lambre–Zhou–Zimmermann (following Tradler in the symmetric case).

Connes' operator on the twisted complex, `B_ν : HH_{n-1}(A, {}_1A_ν) →
HH_n(A, {}_1A_ν)`, RAISES homological degree by 1. Define `Δ` as its **adjoint
under the perfect pairing `(†)`**:

```
⟨ Δα , z ⟩_{n-1}  =  ⟨ α , B_ν z ⟩_n        for all z ∈ HH_{n-1}(A, {}_1A_ν).      (‡)
```

Degrees: `α ∈ HHⁿ`, `z ∈ HH_{n-1}` twisted, `B_ν z ∈ HH_n` twisted, so `⟨α, B_ν z⟩`
is defined, and `(‡)` determines `Δα ∈ D(HH_{n-1} twisted) = HHⁿ⁻¹(A)`. Thus
**`Δ: HHⁿ → HHⁿ⁻¹` lowers cohomological degree by 1**, as required. This is the
honest form of the record's paraphrase `Δ = ∂∘B_σ∘∂⁻¹`: `∂` is the duality iso
`(†)` and `B_σ = B_ν` is the twisted Connes operator; the composite is the
transpose `(‡)`. We implement `(‡)` (unambiguous), NOT a chosen composite.

### 2.2 The symmetric special case (ν inner) is P52-free — the anchor

When `A` is **symmetric**, `ν` is inner, so `{}_1A_ν ≅ A` as bimodules and the
twisted homology collapses to ordinary homology: the pairing `(†)` is simply
`HHⁿ(A) ⊗ HH_n(A) → k`, `⟨f,z⟩ = λ(a_0 f(a_1⊗⋯⊗a_n))` with `λ` the **symmetric**
trace form (`is_symmetric`'s certified nondegenerate `λ`), and `B_ν` is the
**ordinary** Connes `B` that `connes_differentials` already computes. This is
Tradler's theorem, and it needs **no P52 dependency**. We therefore make the
symmetric route the `ν = id` specialization that lands first and de-risks the
whole sign/twist arbitration before the twisted machinery is wired.

### 2.3 The recovered bracket and the BV relation

Grading `|α| = n` for `α ∈ HHⁿ`. The BV relation expresses the Gerstenhaber
bracket as the deviation of `Δ` from being a derivation of the cup product:

```
[α, β]  =  ε · ( Δ(α ∪ β)  −  Δ(α) ∪ β  −  (−1)^{|α|} α ∪ Δ(β) )                  (BV)
```

with an overall sign `ε = ε(|α|,|β|)` (sources differ: some carry
`ε = −(−1)^{(|α|−1)|β|}`). **We do NOT assume `ε`.** The plan tests `(BV)` with a
stated candidate convention, and the sign is **arbitrated once** by the in-window
equality of the derived bracket with the independent Plan-35 bracket
(`gerstenhaber_brackets`) over ODD primes (`3, 5, 32003` — where the `(−1)`
factors are visible). The arbitrated `ε` is then a fixed global convention that
the same code path applies over every field, GF(2) included (§Major-3 rationale
in §7).

**The 7-term relation the record names.** A BV algebra is a Gerstenhaber algebra
`(HH•, ∪, [−,−])` with a degree `−1` operator `Δ` such that `Δ² = 0` and `Δ` is a
**differential operator of order ≤ 2** with respect to `∪` — the *seven-term*
(BV / Koszul) relation

```
Δ(a∪b∪c) = Δ(a∪b)∪c + (−1)^{|a|}a∪Δ(b∪c) + (−1)^{(|a|−1)|b|}b∪Δ(a∪c)
           − Δ(a)∪b∪c − (−1)^{|a|}a∪Δ(b)∪c − (−1)^{|a|+|b|}a∪b∪Δ(c).           (7T)
```

Given a Gerstenhaber algebra and `Δ² = 0`, `(7T)` holds **iff** `[−,−]` is the
`Δ`-generated bracket `(BV)` (Getzler 1994; Koszul 1985). So `Δ²=0` + `(7T)` is
the record's "7-term BV axioms", and `Δ²=0` + `(BV)` is the equivalent 2-term
form. **We self-cert BOTH** — `Δ²=0`, `(7T)` directly, and (cross-engine) `(BV)`
against the independent bracket — so the equivalence is checked, not merely
invoked.

### 2.4 Hypothesis gates — the DELIBERATE pick (record §"implementation must pick")

The transport `(‡)` is a theorem only when `ν` is "nice". The record names three
computable conditions; we settle them exactly, and note a decisive containment
that makes the pick clean:

> **Fact (settled at spec time).** Over any field, `ord(ν) < ∞ and char k ∤
> ord(ν)`  ⟺  `ν is semisimple with finite order`. (An element with `ord(ν) = m`
> and `char k ∤ m` is semisimple because `x^m − 1` is separable when `char ∤ m`;
> conversely a finite-order semisimple element over a field of char `p` has
> `p ∤ ord(ν)`, since the only `p`-power root of unity in char `p` is `1`.)
> Volkov's condition is stated as `char k ∤ ord(ν)` (NOT a gcd — in char 0,
> `0 ∤ m` holds for every `m ≥ 1`, so char 0 imposes only `ord(ν) < ∞`). Hence
> **Volkov's condition (1405.5155) ⊆ LZZ's condition (1405.5325)**; the
> difference `{semisimple, infinite order}` occurs only in char 0.

So the **primary decidable gate is LZZ semisimplicity**, which subsumes Volkov
and is decided by ONE exact test:

* **ν semisimple ⟺ minpoly(ν) is squarefree.** `ν` is the matrix from
  `nakayama_automorphism()` (integers over `GF(p)`, Domain elements otherwise).
  Squarefreeness is exact over a **perfect** field (`GF(p)` and `QQ` are both
  perfect) via `gcd(m, m') = 1` for the minimal polynomial `m` (equivalently the
  characteristic polynomial's radical test). No float, conclusive. Over `GF(p)`
  this coincides with Volkov's `p ∤ ord(ν)`; we report BOTH phrasings in
  provenance.
  * **Up-to-inner caveat (stated, not hidden):** `ν` is defined only up to inner
    automorphism and `nakayama_automorphism()` returns one specific representative
    `G⁻¹Gᵀ`. Semisimplicity of that matrix is what we test; LZZ's hypothesis is on
    "the" Nakayama automorphism. For the algebras in scope (symmetric, quantum
    complete intersections with diagonal ν, self-injective Nakayama) the canonical
    representative is the natural one, and the **arbiter (`(BV)` ≡ independent
    bracket) is the per-instance correctness gate** — a wrong representative would
    fail the arbiter and refuse, never return a wrong Δ silently. Precise scope
    (fix-round): the arbiter pins Δ **modulo cup-derivations** (the data the `(BV)`
    relation constrains) together with `Δ²=0`, not every last coordinate. Recorded
    on the verification page.

* **BIKLZ (self-injective Nakayama, semisimplicity REMOVED)** is the only regime
  that adds genuinely new instances beyond LZZ: a **non-semisimple** ν, which
  (by the Fact) forces `char | ord(ν)` — the char-`p`, `p | N` corner. Decidable
  gate: `is_selfinjective()` AND the quiver is Nakayama (linear `A_n` or a single
  cycle `Z_e`, admissible). **In v1 the non-semisimple case is REFUSED loudly**
  (implementing the general BIKLZ construction — the part they add over LZZ — is
  out of v1 scope; this is a scope decision, not a transcription blocker: the
  paper's formulas ARE transcribed, but only the `e = 1` slice is needed for what
  v1 serves). The important special case `kZ₁/J^N = k[x]/(x^N)` is *symmetric*
  (`ν = id`, trivially semisimple) and is fully served by the symmetric route in
  every characteristic — including `char | N` — so the char-sensitivity
  phenomenon (and the BIKLZ §3.2 explicit-Δ values) are in scope without the
  general BIKLZ construction.

**Routing (per instance, in order; provenance records which fired):**

1. `is_symmetric()` → **symmetric / Tradler** route (ν inner, no twist, ordinary
   Connes B, symmetric trace-form pairing). *P52-free.* Provenance
   `"symmetric (Tradler AIF 2008)"`.
2. else `is_frobenius()` and `minpoly(ν)` squarefree → **semisimple-ν / LZZ**
   route (twisted homology + twisted Connes `B_ν` + twisted pairing). Provenance
   `"semisimple Nakayama automorphism (Lambre–Zhou–Zimmermann 2016; Volkov 2016
   over GF(p): char ∤ ord ν = k)"`.
3. else `is_frobenius()` and self-injective Nakayama and **ν semisimple** →
   already caught by (2); the non-semisimple sub-case (necessarily multi-vertex
   `e ≥ 2`, `char | ord ν`) → **loud refusal**: `"BV operator: self-injective
   Nakayama with NON-semisimple Nakayama automorphism (char | ord ν) — the
   general Bian–Itagaki–Kou–Lyu–Zhou construction (arXiv:2603.04834) is NOT
   implemented in v1 (out of scope); note the e = 1 case k[x]/(x^N) is symmetric
   and served by route (1)"`.
4. else → loud refusal: `"BV operator needs a Frobenius algebra whose Nakayama
   automorphism is semisimple (or a symmetric algebra); this algebra is
   {not Frobenius | Frobenius with non-semisimple ν} — no implemented BV
   hypothesis applies"`.

There is **one transport engine**; the gate decides *whether* it is a theorem and
records *why*. `Δ² = 0` (always) and the bracket arbiter (in-window, GF(p))
finalize correctness per instance.

**Weakly-symmetric-but-not-symmetric algebras.** `is_weakly_symmetric_generic`
already exists (`invariants/frobenius.py`): weakly symmetric = Frobenius with the
IDENTITY Nakayama *permutation* but a possibly **non-inner** ν (exterior algebras
and quantum complete intersections are the standard examples — identity
permutation, diagonal-but-non-inner ν). These are NOT symmetric, so gate 1 does
not fire. They are served by **gate 2 (semisimple-ν / LZZ) exactly when their ν
is semisimple** — which for the diagonal ν of exterior/QCI algebras it is
(`minpoly` squarefree), so they route through the twisted machinery — and by a
**loud refusal otherwise** (a weakly-symmetric algebra whose ν has a non-trivial
unipotent part in bad characteristic falls into the BIKLZ-blocked corner). This
is stated so the reader knows weak symmetry is not a separate route: it collapses
to the semisimple-ν gate or a refusal, never a silent third path.

---

## 3. P52 alignment (prerequisite consumption contract)

P52 (`docs/plans/2026-08-07-plan-52-hh-coefficients.md`) is **not on disk** at
authoring time. P54's twisted routes (§2.1, gate 2) consume P52's
**twisted-coefficient Hochschild homology** `HH_•(A, {}_1A_ν)`. The alignment
point, to be reconciled the moment P52 lands (the implementer's Task B gate):

* **What P54 needs from P52:** for a bimodule `M = {}_1A_σ` (right action twisted
  by an algebra automorphism `σ`), over `GF(p)` in the bar window: (a) the
  homology dimensions `dim HH_n(A, M)`; (b) an explicit **cycle basis** for each
  `HH_n(A, M)` in the (twisted) bar chain coordinates, plus the boundary matrices
  — i.e. the `(reps, image_cols)` quotient data the Plan-35 `_class_coords`
  machinery already consumes for ordinary homology
  (`hochschild/products.py::_generic_homology_quotient`, `homology_classes`).
* **Expected P52 surface (to confirm):** a bimodule-coefficient argument on the
  homology entry point — most likely
  `A.hochschild_homology(top, coefficients=M)` or a
  `hochschild/coefficients.py::twisted_homology_classes(A, sigma, n)` helper.
  P54 wraps whichever P52 ships in a thin `bv/twist.py::twisted_homology_quotient`
  adapter so the rest of P54 is P52-API-agnostic.
* **What P54 adds (NOT in P52):** the **twisted Connes operator** `B_ν` on the
  twisted bar complex (P52 gives the homology *spaces*; the twisted *cyclic*
  operator is P54's, `bv/twisted_connes.py`, citing LZZ's construction), and the
  pairing `(†)`.

* **CROSS-PLAN CONTRACT — the twisted chain-level convention (MUST be mirrored
  into P52's plan doc at its own fix round).** `B_ν` provably descends against
  P52's boundary ONLY if the two agree on the twisted complex convention. Both
  plans pin, verbatim, this single convention block:
  1. **Complex.** `C_n(A, {}_1A_σ) = A ⊗ Ā^{⊗n}` on the unit-adapted basis
     `(s, J)`, `s ∈ 0..m−1`, `J ∈ {1..m−1}^n` — **identical shapes to
     `hochschild/bar.py::_cochain_basis`** (basis index 0 = `1_A`; the same
     ordering the Plan-35 homology machinery uses).
  2. **Twist placement — the LAST face map only.** The Hochschild boundary
     `b = Σ_{i=0}^{n} (−1)^i d_i` has `d_i` unchanged for `0 ≤ i ≤ n−1`
     (`d_i(a_0⊗⋯⊗a_n) = a_0⊗⋯⊗a_i a_{i+1}⊗⋯⊗a_n`) and the wrap-around face
     **twisted on the coefficient slot**:
     `d_n(a_0⊗⋯⊗a_n) = σ(a_n)·a_0 ⊗ a_1 ⊗ ⋯ ⊗ a_{n−1}` (σ applied to the last
     factor as it multiplies into the `A`-slot). This is the placement the
     pairing compatibility `λ(ab) = λ(b ν(a))` forces (§2.1); P52's twisted
     boundary MUST use exactly this `d_n`.
  3. **Normalization.** Degenerate chains (any bar slot `= 1_A`, i.e. basis
     index 0 in a `J`-slot) are killed — the normalized complex, matching
     `hochschild/cyclic.py`.
  4. **Twisted rotation / `B_ν`.** `t_σ(a_0⊗⋯⊗a_n) = (−1)^n · σ(a_n) ⊗ a_0 ⊗ ⋯ ⊗
     a_{n−1}`, `s(a_0⊗⋯⊗a_n) = 1 ⊗ a_0 ⊗ ⋯ ⊗ a_n`, `B_ν = (Σ_{i=0}^{n} t_σ^i)∘s`
     on the normalized complex (the ordinary `connes_B_matrix` of
     `hochschild/cyclic.py` with the untwisted rotation replaced by `t_σ`). The
     Koszul sign is inherited from that engine; the arbiter finalizes it.
  5. **Basis ordering.** P52 returns cycle bases as coordinate columns over the
     ordering of (1), so `bv/twist.py` and `_class_coords` read them without
     re-indexing. If P52's ordering differs, `bv/twist.py` permutes ONCE and the
     permutation is pinned by a round-trip test.
  A **self-cert (Task E)** rebuilds P52's twisted `b` from this convention and
  asserts `b∘b = 0` and `b∘B_ν + B_ν∘b = 0` on the twisted complex — the concrete
  proof that `B_ν` descends against P52's boundary. If P52's delivered convention
  differs, THIS is the reconciliation point (adjust `bv/twist.py`, never the math).

* **Fallback if P52 slips or lands with a different API:** the **symmetric route
  (§2.2) has NO P52 dependency** and is a complete, shippable deliverable on its
  own (Tradler anchor + `k[x]/(x^N)` all-char incl. the BIKLZ §3.2 value pin +
  `k[x,y]/(x²,y²)` [4,4,5,6]). If P52 is unavailable at implementation time, P54
  ships the symmetric route and the semisimple-ν route is **shipped DISABLED
  behind a loud typed refusal** (`"twisted-coefficient homology unavailable
  (P52 not present) — the semisimple-ν BV route is disabled"`), with a ledger
  entry P80 reconciles. **Consequence to state plainly (MAJOR 2):** the symmetric
  fallback NEVER exercises the twist machinery, so the twist-direction determiner
  (the QuantumCI arbiter, §5) and the `b∘B_ν+B_ν∘b=0` descent self-cert are
  **P52-gated** — until P52 lands, the twisted route is disabled, not silently
  unvalidated. The QuantumCI LZZ ties then move to the P52-follow-up. **This is
  the explicit alignment flag the task requested.**

---

## 4. Public API (`core/algebra.py`)

One method, mirroring the Plan-35 product dispatch signatures:

```python
A.bv_operator(top, engine="auto", max_cells=4_000_000)
    # The Batalin–Vilkovisky operator Δ: HH^n -> HH^{n-1} for 1 <= n <= top,
    # on the recorded HH basis, exact. Requires a Frobenius algebra whose
    # Nakayama automorphism is semisimple (symmetric algebras included as the
    # ν-inner anchor); loud typed refusal otherwise. engine: 'auto' (GF(p) bar/tt
    # route in v1); 'bar' is the explicit GF(p) route (loud off GF(p)); 'cs' is
    # reserved for the P51 past-window/off-GF(p) enhancer (loud "not available"
    # until P51 lands — never a silent fallback).
```

Return object: a frozen **`BVOperator`** in `hochschild/products.py` beside
`HHProducts`/`ConnesB`, with the same `.blocks()` single-shape serialization
consumed identically by both runners:

```python
class BVOperator:
    top            # int
    hh_dims        # [dim HH^0 .. dim HH^top]
    matrices       # {n: rows-of-str}, n in 1..top; Δ_n is dim HH^{n-1} x dim HH^n
                   #   (rows indexed by the OUTPUT degree HH^{n-1}, mirroring ConnesB)
    ranks          # {n: int}
    hypothesis     # provenance string (§2.4 routing) — WHICH hypothesis certified it
    nakayama       # {"matrix": [[str]], "semisimple": bool, "order": int|None,
                   #  "inner": bool}  — the ν data the gate used
    basis          # "bar/GF(p)" — WHICH basis the constants live in (Plan-35 doctrine)
    window         # int: the certified degree window served (bracket-arbiter bound)
    bracket_check  # {"engine": "...", "window": int, "agrees": bool}
                   #  the in-window derived==independent arbiter result (self-cert flag)
    derived_bracket   # an HHProducts(kind="bracket") built FROM Δ via (BV) — serialized
                      #  identically to gerstenhaber_brackets so the oracle diffs tables
    references     # ["bv_tradler","bv_lzz","bv_volkov","cyclic","bracket", ...]
    basis_classes / chain_basis / differentials   # Plan-35 explicit-reps parity (best-effort)

    def blocks(self): ...   # {"kind":"bv_operator", "top", "hh_dims", "matrices",
                            #  "ranks", "hypothesis", "nakayama", "basis", "window",
                            #  "bracket_check", "derived_bracket": <bracket blocks>,
                            #  "references"}  (+ explicit-reps fields when present)
```

Degenerate degrees (`dim HHⁿ = 0`) are present-but-empty, never omitted (Plan-35
convention: a zero matrix is a statement). `Δ_0 = 0` by degree (no `HH^{-1}`), so
`matrices` starts at `n = 1`.

---

## 5. Engine routing and the refusal matrix

`engine="auto"` resolves per request; **one route serves one call end-to-end** —
a single call's Δ matrices never mix bases across degrees (Plan-35 rule).

| Algebra / field | route | Connes B | pairing | needs P52? |
|---|---|---|---|---|
| symmetric, GF(p), `top` in bar window | Tradler (§2.2) | ordinary `connes_b_tables` | symmetric trace-form `λ` | no |
| Frobenius, semisimple ν, GF(p), in window | LZZ (§2.1) | twisted `B_ν` on `HH_•(A,{}_1A_ν)` | ν-twisted `(†)` | **yes** |
| Frobenius, non-semisimple ν (self-inj. Nakayama, char\|N) | — | — | — | **loud refusal** (BIKLZ blocked) |
| not Frobenius | — | — | — | **loud refusal** |
| any Frobenius+ss ν, off GF(p) OR past window | — | — | — | **loud refusal** in v1 (P51 enhancer, §8) |

* Explicit `engine="bar"` keeps an honest wall (error, never silent fallback), the
  Plan-35 contract. `engine="cs"` raises `"not available until P51"` in v1.
* `max_cells` guards the bar/twisted-bar blow-up exactly as Plan-35 (the twisted
  complex has the same cardinality as the untwisted one; the same cochain/chain
  pair guard applies before any matrix is built).
* **Perfect-pairing self-cert (a COARSE sanity check, NOT the direction
  determiner):** the transport builds `P_n` for `(†)` in each degree and requires
  it **square and invertible** (`dim HHⁿ = dim HH_n(twisted)` and `det P_n ≠ 0`).
  This rejects a grossly wrong coefficient bimodule. **It does NOT determine the
  twist direction:** `{}_1A_ν` and `{}_1A_{ν⁻¹}` give the *same* twisted-homology
  dimensions in general (ν ↦ ν⁻¹ is an isomorphism of the twisted complexes'
  underlying graded spaces), so both directions pass the dim-match. A one-time
  **pin task (Task E, §9)** verifies `dim HH_•(A,{}_1A_ν) = dim HH_•(A,{}_1A_{ν⁻¹})`
  on a discriminating instance so this demotion is a checked fact, not a claim.
* **The bracket ARBITER determines the twist direction.** For each candidate twist
  the transport produces a Δ; the direction is the one whose **derived bracket
  equals the independent Plan-35 bracket in-window** (`bracket_check.agrees`). The
  concrete discriminating instance the plan REQUIRES exercised is
  **`QuantumCI(q=2)`**: its Nakayama automorphism has order 4 with `ν ≠ ν⁻¹`
  (over `GF(5)`, `ν = diag(1, 3, 2, 1)` on the `{1,x,y,xy}`-type basis, `ν² =
  diag(1,4,4,1) ≠ 1`), so the two twist directions give genuinely different Δ and
  only one reproduces the bracket. A symmetric or order-2 (`ν = ν⁻¹`) instance
  can NOT discriminate direction and must not be used to "validate" it.
  * **Direction validation is therefore P52-GATED.** The symmetric-only fallback
    (§3) never builds a twisted complex, so it never exercises the direction
    determiner. If P52 slips, the twisted route ships **DISABLED** behind a loud
    refusal — never silently unvalidated with an unchecked direction.
* Cross-engine agreement of the *derived bracket* with the independent Plan-35
  bracket is checked in-window and RECORDED in `bracket_check`; a mismatch under
  BOTH twist directions is a **loud refusal** (`"BV transport does not reproduce
  the independent Gerstenhaber bracket in-window under either twist — the
  hypothesis/convention does not certify for this instance"`), never a silent
  wrong Δ. This is the load-bearing correctness gate for the non-symmetric routes.

---

## 6. Webapp + GUI (`bv_operator` kind on the hochschild theme)

Mirrors the Plan-35 product-kind wiring precisely (§3 of Plan 35; §Task H of
Plan 40 is the seven-touchpoint checklist).

* **Schema:** `bv_operator:0..n` accepted wherever `bracket:0..n` /
  `connes_b:0..n` are (algebra-level range kind, no module block). It is a member
  of `PRODUCT_KINDS` (`hpc/spec.py:89`) — routes through `_product_object` /
  `_dispatch`'s product branch, not the module dispatch. **Canonical-key
  stability:** the kind string flows through `canonical_key` via the compute list
  exactly like the other product kinds; every existing request keeps its
  byte-identical key (frozen-golden gated). The GUI push order places
  `bv_operator` immediately after `connes_b` (the curated `_GUI_ORDER`, so the
  reachability gate `tests/webapp/test_curated_reachability.py` extends by one).
* **Runners:** dispatch lands once in `hpc/spec.py::_product_object`
  (`method = {... "bv_operator": A.bv_operator}[kind]`) and its Pyodide twin
  `docs/gui/runner.py` (the `name in ("cup","cap","bracket","connes_b")` elif at
  `runner.py:939` gains `"bv_operator"`; the `_snip` recipe at `spec.py`/
  `runner.py:1180` gains `"bv_operator": "A.bv_operator(%d)"`). Both runners
  **byte-identical** on the block (`json.dumps(sort_keys=True)` parity test).
* **Estimator:** `bv_operator` sizes like `hh_cohomology` at the same `top` for
  tier routing (its cost is HH + one Connes + one bracket + linear-algebra
  transposes); GUI wait-estimate entry added (`ETA_MODEL` scalar/curve).
* **GUI:** one checkbox `qlgui-bv_operator` + degree picker in both vendored
  `gui.js` copies (`docs/gui/gui.js` and `webapp/static/gui/gui.js`, kept
  byte-identical), added to `S.ids` and the push list; a block renderer reusing
  `matrixGrid` (indexed grids) and the `matIsZero` one-line zero-matrix
  convention; the block renders the Δ matrices per degree, the `hypothesis`
  provenance line, the ν data, and the `bracket_check` agreement flag. New EN/ES
  i18n keys for every user-facing string (`inv.bv_operator`,
  `block.bv_operator.title`, per-row labels, the hypothesis/provenance gloss —
  `webapp/server/i18n/en.json` + `es.json`).
* **Report (`trace/products.py` + `trace/results_html.py`):** a BV chapter in the
  worked-steps bundle — definitional preamble (`Δ: HH^n → HH^{n-1}`, `Δ²=0`, the
  `(BV)` relation), the recorded HH basis dims, the **hypothesis + ν provenance
  gloss** (which of Tradler / LZZ / Volkov certified it, and the plain-language
  meaning), the pairing-matrix invertibility line (the perfect-pairing
  self-cert), the Δ matrices as indexed grids, and the **arbiter statement** ("no
  external oracle: the bracket recovered from Δ equals the independently computed
  Gerstenhaber bracket in this window"). A **loud drift gate** asserts the
  narrated dims equal the block's. Citations ship from the registry (§Task G).
* **Curated cache examples:** `bv_operator` is added ONLY to Frobenius curated
  examples (the dual-numbers `k[x]/(x²)` and a QuantumCI example), at each
  example's existing HH top; deep non-Frobenius or dim≥220 examples get nothing
  (they already omit products, `webapp/precomputed/manifest.yaml`). Regeneration
  reuses the Plan-35 discipline: every pre-existing block byte-identical, the only
  addition the `bv_operator` block; manifest comment updated. If a curated
  example is not Frobenius the request simply does not carry the kind.

---

## 7. Testing, oracles, verification (Plan-32 classes)

* **`oracle_selfcert`** (`tests/hochschild/test_bv_identities.py`): `Δ² = 0`
  (`Δ_{n-1}∘Δ_n = 0` on HH) over several primes; the **seven-term relation `(7T)`
  of §2.3** (Δ is a differential operator of order ≤ 2 — checked directly, so the
  BV-algebra characterization is verified, not merely invoked); the perfect-pairing
  certificate (`P_n` square + invertible in-window); `bracket_check.agrees is True`
  recorded. All routes. *(No `Δ(1)=0` pin — `Δ_0` is undefined by degree, there is
  no `HH^{−1}`; the family starts at `n = 1`.)*
* **GF(2) IS SERVED** (`tests/hochschild/test_bv_gf2.py`, `oracle_selfcert`).
  **Rationale:** the transpose arrangement and the BV sign `ε` are a *global
  code-path property* fixed ONCE by the odd-prime arbiter (§2.3), not a per-
  instance choice; over GF(2) the identical code path runs and is certified per
  instance by `Δ² = 0`, the **BV relation mod 2** (still a genuine constraint —
  the `(−1)` factors collapse but `Δ(a∪b) − Δ(a)∪b − a∪Δ(b)` must still equal the
  independent GF(2) bracket), and pairing invertibility. The test runs on
  `products_loop_gf2` (the existing GF(2) loop algebra from the Plan-35 corpus)
  and asserts those three certs pass; a documented honest-scope caveat records
  that char-2 certification is `Δ²=0` + BV-relation + pairing-invertibility, not
  sign-determination.
* **`oracle_crossengine`** (`tests/hochschild/test_bv_bracket_arbiter.py`): **the
  derived bracket ≡ `A.gerstenhaber_brackets(top)` in-window**, table-for-table,
  over ODD primes `(3, 5, 32003)` — the sign/twist arbiter (where the `(−1)`
  factors are visible and determine `ε` and the twist direction). QuantumCI
  (semisimple diagonal ν, non-symmetric — `QuantumCI(q=2)` for the direction
  discrimination of §5) is the flagship non-trivial-twist instance here.
* **`oracle_literature`**:
  * `tests/hochschild/test_bv_tradler_symmetric.py` — `k[x,y]/(x²,y²)`: symmetric,
    ν = id, `HH_• = [4,4,5,6]` pin (existing); the BV structure matches Tradler
    (Δ nonzero where the theory predicts; the derived bracket ≡ independent). The
    Tradler anchor.
  * `tests/hochschild/test_bv_kx_xN_charsens.py` — `k[x]/(x^N)`: the char-sensitive
    HH dim pattern (`dim HH⁰ = N`; for `n ≥ 1`, `dim HHⁿ = N−1` if `char ∤ N`,
    `= N` if `char | N`), sourced from **BIKLZ §3.2** and the classical
    truncated-polynomial computation (NOT a vague "folklore" — the §0 transcribed
    draft + the standard `k[x]/(x^N)` HH computation are the two agreeing
    sources); served on the symmetric route in BOTH regimes (`char ∤ N` over
    `GF(32003)`, `char | N` over `GF(N)` with `N` prime), with `Δ² = 0` in both.
  * `tests/hochschild/test_bv_biklz_kz1.py` — **LIVE `oracle_literature` assert**
    (NOT fenced): the explicit BIKLZ §3.2 Δ VALUES for `kZ₁/J^N = k[x]/(x^N)`,
    validating the symmetric-route Δ against the §0 transcribed formulas in BOTH
    char regimes — `char ∤ N`: `Δ(y·x₀ᵃ·zᵇ) = (bN+N−a−1)·x₀ᵃ·zᵇ`; `char | N`:
    `Δ(w·x₀ᵃ·z'ᵇ) = −a·x₀ᵃ⁻¹·z'ᵇ`, `Δ(w·z'ᵇ)=0`. The implementer copies the
    §3.2 section identifier + the printed formulas into the test docstring
    **re-verified against `arxiv.org/html/2603.04834` at implementation time**,
    maps the paper's basis symbols to the engine's HH classes, and asserts
    equality. This is the only exact-Δ-VALUE literature oracle for the symmetric
    route — it MUST be live. Never fabricated (any symbol the HTML does not print
    is not asserted).
  * `tests/hochschild/test_bv_quantumci.py` — **P52-gated** QuantumCI ties: the
    derived bracket agrees with the existing Plan-35 QuantumCI bracket pins
    (`qci_hh_oracle`), semisimple-ν / LZZ route, diagonal ν. Skips with a recorded
    reason until P52 lands (the twisted route is disabled in the symmetric-only
    fallback, §3); lands live in Task E.
* **`qpa`** (`tests/qpa/test_bv_qpa.py`): probe live (the `NamesGVars()` precedent,
  `tests/qpa/test_products_qpa.py`) that QPA 1.37 exposes NO BV / Δ surface
  (expected: none — Plan 35 established QPA has no HH product surface at all). An
  honest **skip that FAILS if that ever changes**; the covering theory oracle is
  named (the Δ²=0 + arbiter + Tradler/LZZ literature).
* **Runner/GUI contract** (`tests/webapp/test_bv_block_p54.py`, unmarked —
  extras-gated dir, the m0729 cross-runner-pair pattern): one NEW golden fixture
  `bv_operator_kxy` in `_runner_goldens.json` (existing goldens untouched,
  documented in `test_runner_delegation.py`'s change-log); both runners
  byte-identical on the block; a report test reading the Δ grids back out via the
  `_matrix_grid` helper.
* `docs/verification.md`: a P54 subsystem→oracles→tests row; the **honest-scope
  entries**: (a) GF(p)-window-bounded in v1 (P51 optional enhancer for
  past-window/off-GF(p)); (b) the **general** BIKLZ construction for `e ≥ 2`
  self-injective Nakayama with non-semisimple ν is refused (out of v1 scope) —
  but the `e = 1` explicit Δ VALUES (BIKLZ §3.2) are a LIVE oracle for the
  symmetric route (no blocked pin remains); (c) QPA has NO BV surface (the theory
  oracle covers it); (d) the ν-up-to-inner caveat and the arbiter that guards it;
  (e) **char 2** — the sign/transpose convention is INHERITED from the odd-prime
  anchor; GF(2) certification is `Δ²=0` + the BV relation + pairing invertibility,
  not sign-determination. Recount the class table (`tests/release/
  test_oracle_classes.py` drives the numbers — collect, paste, re-run to green).

---

## 8. Non-goals (v1)

* **Off-GF(p) and past-window Δ.** The bracket arbiter is the load-bearing oracle
  and it is `GF(p)`-window-bounded (Plan 35). v1 is `GF(p)`-in-window and SAYS so.
  **P51 (native past-window / any-Domain bracket) is an OPTIONAL ENHANCER**, not a
  prerequisite: if P51 has merged, a follow-up wires `engine="cs"` to extend both
  Δ and the arbiter past the window / off GF(p). Until then `engine="cs"` refuses
  loudly (`"not available until P51"`). Recorded in the GUI-deferral ledger for
  P80.
* **The GENERAL BIKLZ non-semisimple construction** (`e ≥ 2` self-injective
  Nakayama with non-semisimple ν, §2.4 gate 3). Specified but refused loudly in
  v1; a successor plan implements the explicit construction. NOTE this is distinct
  from the `e = 1` (`k[x]/(x^N)`) explicit Δ VALUES, which ARE live in v1 via the
  symmetric route + the BIKLZ §3.2 oracle — nothing about the e = 1 values is
  deferred.
* **Class representatives in the public output.** As in Plan 35, the structure
  constants pin Δ up to the recorded basis; representative canonicality is its own
  problem (the explicit-reps capture is best-effort provenance only).
* **A dedicated BV showcase example.** The Frobenius curated examples gain the
  kind; a bespoke showcase is a later curation decision.
* **Per-step trace i18n.** Worked-steps stay EN (existing contract); the
  GUI/webapp block labels ARE bilingual.

---

## 9. TDD task breakdown

Buckets (auto by directory, `tests/conftest.py`): `tests/hochschild/` → **deep**;
`tests/qpa/` → **qpa**; `tests/webapp/` → fast/unmarked. Run new deep tests by
path during development; finish with a `-m deep` spot-run of touched files. Test
cmd throughout:
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q <path>`.
Conventional commits; green at every commit; branch `plan-54-bv-operator` off
`dev`.

---

### Task A — ν hypothesis gate (`bv/hypothesis.py`)

**Files:** create `src/quiverlab/hochschild/bv/__init__.py`,
`src/quiverlab/hochschild/bv/hypothesis.py`; test
`tests/hochschild/test_bv_hypothesis.py`.

**Interfaces:**
- Consumes: `A.is_symmetric()`, `A.is_frobenius()`, `A.is_selfinjective()`,
  `A.nakayama_automorphism()` (matrix); sympy for the minimal-polynomial
  squarefree test over the Domain; `A.domain` for char.
- Produces:
  ```python
  @dataclass(frozen=True)
  class BVHypothesis:
      applies: bool
      route: str          # "symmetric" | "semisimple" | None
      label: str          # human provenance (the §2.4 strings)
      nu_semisimple: bool
      nu_order: int | None   # ord(ν) when finite & cheaply found, else None
      nu_inner: bool         # symmetric => True
      refusal: str | None    # the loud message when applies is False
  def classify_bv(A) -> BVHypothesis
  def nu_is_semisimple(nu_matrix, domain) -> bool   # minpoly squarefree, exact
  ```
  `classify_bv` implements the §2.4 order: symmetric → semisimple → (self-inj.
  Nakayama non-semisimple → refuse with the BIKLZ-blocked message) → refuse.

- [ ] **Step 1 — failing tests.** Symmetric `k[x,y]/(x²,y²)` → `route=="symmetric"`,
  `nu_inner`; QuantumCI (diagonal ν) → `route=="semisimple"`, `nu_semisimple`;
  `k[x]/(x^p)` over `GF(p)` → symmetric (ν=id), `route=="symmetric"` (NOT the
  BIKLZ refusal — it is symmetric); a genuinely non-symmetric self-injective
  Nakayama with `char | ord(ν)` (multi-vertex `kZ_e/J^L`, chosen so ν is a
  non-semisimple permutation-type map in char p) → `applies is False`,
  `"general ... construction is NOT implemented in v1"` in `refusal` (the §2.4
  gate-3 message — the e≥2 non-semisimple case is out of scope, NOT a
  transcription blocker); a non-Frobenius algebra (`kA_2`) → `applies is False`,
  the non-Frobenius refusal.
- [ ] **Step 2 — verify failure** (`ModuleNotFoundError: quiverlab.hochschild.bv`).
- [ ] **Step 3 — implement.** `nu_is_semisimple`: build the sympy matrix over the
  Domain, compute its minimal polynomial, test `gcd(m, m') == 1` (squarefree).
  Over `GF(p)` use sympy `GF(p)` domain; over `QQ` rationals. **Adjust to reality:**
  confirm `nakayama_automorphism()`'s matrix orientation (columns = images) and
  that sympy's minpoly over `GF(p)` is available; if a Domain (e.g. `GF(p^n)`) has
  no direct sympy minpoly, fall back to the squarefree test on the characteristic
  polynomial's radical (equivalent for the semisimplicity verdict via
  `gcd(charpoly, charpoly')` degree bookkeeping) or refuse loudly. `nu_order`:
  cheap bounded power-search (cap ~ a few hundred) recording `None` past the cap —
  it is provenance only, never load-bearing (the squarefree test is the gate).
- [ ] **Step 4 — pass.** `-m pytest tests/hochschild/test_bv_hypothesis.py -q`.
- [ ] **Step 5 — commit** `feat(hochschild): BV hypothesis gate — semisimple-ν
  (minpoly squarefree) + symmetric anchor, loud refusals`.

---

### Task B — the pairing and the symmetric-route transport (`bv/transport.py`)

**Files:** create `src/quiverlab/hochschild/bv/transport.py`; extend
`hochschild/products.py` with `class BVOperator`; test
`tests/hochschild/test_bv_symmetric.py` (the Tradler anchor, **P52-free**).

**Interfaces:**
- Consumes: `engine.tt_calculus.cohomology_classes` / `homology_classes` (GF(p) bar
  basis, already used by `gfp_product_tables`); `connes_b_tables` (ordinary Connes
  B, symmetric route); `invariants.frobenius.frobenius_form_generic` (the trace
  form `λ` and its Gram — for `GF(p)` route through the engine unit-adapted basis
  consistently with `nakayama_automorphism`'s `GF(p)` branch); the shared
  `_class_coords` solver in `products.py`.
- Produces:
  ```python
  def pairing_matrix(A, n, coh_classes, hom_classes, lam, ...) -> matrix   # (†) matrix
  def bv_matrices_symmetric(A, top, max_cells) -> BVOperator               # §2.2
  ```
  `Δ_n` via `(‡)`: with `P_n[i][j] = ⟨α_i, z_j⟩` (dim HHⁿ × dim HH_n) the pairing
  matrices and `B[k][j]` (dim HH_n × dim HH_{n−1}) the ordinary Connes matrix
  `HH_{n-1} → HH_n` (`connes_b_tables` matrices — mind the `ConnesB` row/column
  orientation, `matrices[n]` is `hh_{n+1} × hh_n`), solving
  `P_{n-1}ᵀ · Δ_n = Bᵀ · P_nᵀ` gives

  ```
  Δ_n = (P_{n-1}ᵀ)⁻¹ · Bᵀ · P_nᵀ .
  ```

  This is the literal transpose bookkeeping of `(‡)` (derived, not guessed —
  see the Methodology note); any residual Koszul/pairing sign is **derivation-
  corrected and arbiter-pinned** (Task D), not assumed. **Perfect-pairing
  self-cert:** require every `P_n` square and invertible in-window, else loud
  refusal.
- [ ] **Step 1 — failing tests.** `k[x]/(x^N)` over `GF(32003)` (char ∤ N) and
  over `GF(N)` (char | N): `bv_operator` returns a `BVOperator` with `Δ²=0`,
  `hypothesis` symmetric, `hh_dims` matching the char-sensitive pattern; the
  pairing matrices are invertible; `k[x,y]/(x²,y²)` returns Δ with the [4,4,5,6]
  hh_dims. (No bracket arbiter yet — Task D.)
- [ ] **Step 2 — verify failure.**
- [ ] **Step 3 — implement** `pairing_matrix` and `bv_matrices_symmetric`; add
  `BVOperator` to `products.py` with `.blocks()`. Concretely, **reuse the Plan-35
  evaluation machinery** rather than re-deriving it:
  * A cohomology class `f ∈ HHⁿ` is a cochain `Ā^{⊗n} → A` in the engine's
    representation; **evaluate it on a bar chain tuple with the SAME cochain-on-
    chain application `engine.tt_calculus` already uses for the cup/cap product
    matrices** (`cup_product_matrix`/`cap_product_matrix` evaluate cochains on the
    bar basis — reuse that evaluation path, do not write a second one).
  * The pairing entry is then `⟨f, z⟩ = λ( a_0 · f(a_1⊗⋯⊗a_n) )`: apply `λ` (the
    `frobenius_form_generic` covector; over `GF(p)` obtain it consistently with
    `nakayama_automorphism`'s `GF(p)` engine branch on the unit-adapted basis) to
    the structure-constant product `a_0 · f(⋯)`. Assemble over the ordered class
    bases `coh(n)`/`hom(n)` exactly as `gfp_product_tables` orders them
    (`coh(n).reps`, `hom(n).reps`, `.coords`) so `P_n` is in the SAME basis the
    constants live in.
  * **Descent to classes** reuses `products.py::_class_coords` /
    `hochschild.basis_reps.classes_from_columns` — the identical
    coordinates-in-the-class-basis solve the Plan-35 product tables use — so the
    pairing is a chain-level bilinear form read off in class coordinates, never a
    bespoke re-basing.
  **Adjust to reality:** confirm the `tt_calculus` cochain-evaluation entry point
  name and that `frobenius_form_generic` is reachable on the `GF(p)` route (if the
  engine exposes `λ` only through the unit-adapted structure constants, thread it
  through there); pin the final Koszul/transpose sign by the Task-D arbiter and
  record the choice in a comment.
- [ ] **Step 4 — pass.** `-m pytest tests/hochschild/test_bv_symmetric.py -q`.
- [ ] **Step 5 — commit** `feat(hochschild): symmetric-route BV operator via the
  trace-form pairing + ordinary Connes B (Tradler anchor, P52-free)`.

---

### Task C — the public method + refusals (`core/algebra.py`)

**Files:** modify `src/quiverlab/core/algebra.py`; test
`tests/hochschild/test_bv_api.py`.

- [ ] **Step 1 — failing tests.** `A.bv_operator(top)` exists; unknown `engine`
  raises the standard `QuiverlabError`; `engine="cs"` raises `"not available until
  P51"`; `engine="bar"` off `GF(p)` raises; a non-Frobenius `kA_2` raises the
  non-Frobenius refusal; a non-semisimple self-inj. Nakayama raises the BIKLZ-blocked
  refusal; the returned object's `.blocks()` has the documented shape.
- [ ] **Step 2 — verify failure.**
- [ ] **Step 3 — implement** `Algebra.bv_operator` mirroring
  `_product_dispatch`/`connes_differentials`: validate `engine`, `classify_bv`,
  route symmetric now (Task B) and semisimple to `bv_matrices_semisimple` (Task E,
  a loud `"P52 twisted homology unavailable"` stub until Task E lands), stamp
  `hypothesis`/`nakayama`/`basis`/`window`.
- [ ] **Step 4 — pass.** **Step 5 — commit** `feat(core): Algebra.bv_operator
  public method + engine routing + loud typed refusals`.

---

### Task D — the derived bracket + the arbiter (`bv/bracket.py`)

**Files:** create `src/quiverlab/hochschild/bv/bracket.py`; wire
`derived_bracket` + `bracket_check` into `BVOperator`/`bv_matrices_*`; tests
`tests/hochschild/test_bv_identities.py` (`oracle_selfcert`),
`tests/hochschild/test_bv_bracket_arbiter.py` (`oracle_crossengine`).

**Interfaces:**
- Consumes: the Δ matrices, the Plan-35 `cup_products` and `gerstenhaber_brackets`
  tables (same call's engine/basis), the (BV) relation `(BV)`.
- Produces: `derived_bracket_tables(A, bv, cup, top) -> HHProducts(kind="bracket")`
  (build `[α,β]` from Δ and ∪ via `(BV)` with the candidate sign `ε`), and
  `bracket_check(derived, independent, window) -> {"agrees", "window", ...}`.
- [ ] **Step 1 — failing tests (SYMMETRIC instances only — the twisted QuantumCI
  ties live in Task E, P52-gated).** `Δ²=0` and the seven-term relation `(7T)`
  over `(3,5,32003)`; the derived bracket ≡ `gerstenhaber_brackets` in-window over
  ODD primes on `k[x,y]/(x²,y²)` and `k[x]/(x^N)`; the perfect-pairing certificate
  asserted; **GF(2) served** — on `products_loop_gf2`, assert `Δ²=0`, the BV
  relation mod 2 vs the independent GF(2) bracket, and pairing invertibility all
  pass (the same code path, convention inherited from the odd-prime anchor).
- [ ] **Step 2 — verify failure.**
- [ ] **Step 3 — implement.** Pin the sign `ε` ONCE by making the odd-prime
  arbiter pass (Plan-21 precedent — the in-window derived≡independent equality IS
  the sign definition; document the chosen `ε` in the docstring as arbitrated, not
  assumed). Implement `(7T)` directly. Wire `bracket_check` into the transport; a
  mismatch is a **loud refusal**.
- [ ] **Step 4 — pass.** **Step 5 — commit** `feat(hochschild): bracket recovered
  from Δ via the BV relation + in-window arbiter vs the independent Gerstenhaber
  bracket (sign convention arbitrated, Plan-21 precedent)`.

---

### Task E — the twisted route (`bv/twisted_connes.py`, `bv/twist.py`) — **P52-gated**

**Files:** create `src/quiverlab/hochschild/bv/twisted_connes.py`,
`src/quiverlab/hochschild/bv/twist.py` (the P52 adapter, §3); extend
`bv/transport.py` with `bv_matrices_semisimple`; test
`tests/hochschild/test_bv_twisted.py`.

**GATE:** begin ONLY after P52 is merged to `dev`. First action: read P52's
delivered API and reconcile `bv/twist.py::twisted_homology_quotient` to it
(§3 alignment); if P52 is absent, this task is DEFERRED with a ledger entry and
`bv_matrices_semisimple` stays the loud stub from Task C (v1 ships symmetric-only —
the §3 fallback).

**Interfaces:**
- Consumes (P52): twisted homology `HH_•(A, {}_1A_σ)` cycle bases + boundaries,
  built to the §3 CROSS-PLAN CONVENTION block (twist on the last face map,
  normalized, basis ordering matched via `bv/twist.py`).
- Produces: `twisted_connes_matrix(A, sigma, n, max_cells)` — `B_σ = (Σ_{i=0}^{n}
  t_σ^i)∘s` with `t_σ(a_0⊗⋯⊗a_n) = (−1)^n σ(a_n)⊗a_0⊗⋯⊗a_{n−1}` on the normalized
  twisted complex (§3 convention, citing LZZ); `bv_matrices_semisimple(A, top,
  max_cells)` — the §2.1 transport, twist direction determined by the ARBITER
  (§5), NOT by dim-match.
- [ ] **Step 1 — failing tests (P52-gated; skip-with-reason if P52 absent):**
  * **Descent self-cert (MAJOR 4):** rebuild P52's twisted boundary `b` from the
    §3 convention and assert `b∘b = 0` and `b∘B_σ + B_σ∘b = 0` on the twisted
    complex — the concrete proof `B_σ` descends against P52's boundary.
  * **Twist-dim pin (MAJOR 2):** `dim HH_•(A, {}_1A_ν) = dim HH_•(A, {}_1A_{ν⁻¹})`
    in-window on `QuantumCI(q=2)` — the checked fact that dim-match does NOT
    determine direction.
  * **Direction discrimination:** on `QuantumCI(q=2)` (`ν` of order 4, `ν ≠ ν⁻¹`;
    over `GF(5)`, `ν = diag(1,3,2,1)`), the two twist directions give different Δ;
    assert exactly ONE has `bracket_check.agrees is True`, and that `bv_operator`
    auto-selects it. `Δ²=0` and `(7T)` on the selected Δ.
  * **QuantumCI ties** (`test_bv_quantumci.py`, moved here from Task D): the
    derived bracket ≡ the Plan-35 QuantumCI bracket pins (`qci_hh_oracle`),
    semisimple/LZZ route.
- [ ] **Steps 2–5** as the TDD ritual; commit `feat(hochschild): twisted-coefficient
  BV route (semisimple ν) — twisted Connes B_σ + arbiter-determined twist direction
  + descent self-cert, over P52 twisted homology`.

---

### Task F — symmetric-route char-sensitivity + BIKLZ value literature batteries

**Files:** `tests/hochschild/test_bv_tradler_symmetric.py`,
`test_bv_kx_xN_charsens.py`, `test_bv_biklz_kz1.py`. All `oracle_literature`.
Concrete pins per §7. **`test_bv_biklz_kz1.py` is a LIVE assert** — it copies the
BIKLZ §3.2 section identifier + printed formulas (§0 draft) into the docstring
**re-verified against `arxiv.org/html/2603.04834` at implementation time**, maps
the paper's basis symbols (`x₀, y, z` / `x₀, y, w, z'`) to the engine's HH
classes, and asserts the symmetric-route Δ VALUES match in both char regimes.
(`test_bv_quantumci.py` is NOT here — it is the P52-gated twisted tie in Task E.)

- [ ] Steps per §7; commit `test(hochschild): BV symmetric-route literature
  batteries — Tradler k[x,y]/(x²,y²), k[x]/(x^N) char-sensitivity, BIKLZ §3.2
  kZ₁/J^N explicit Δ values (live)`.

---

### Task G — citations

**Files:** `src/quiverlab/citations/registry.py` (+ `references.bib`).

Add (verified BibTeX only, the `_r(key, bibtex_key, kind, title, annotation,
*tags)` precedent at `registry.py`):
- `bv_tradler` → Tradler, Ann. Inst. Fourier 58 (2008) 2351–2379, DOI
  10.5802/aif.2417 (verified).
- `bv_lzz` → Lambre–Zhou–Zimmermann, J. Algebra 446 (2016) — semisimple ν
  (verified arXiv:1405.5325; confirm volume/pages at implementation).
- `bv_volkov` → Volkov, arXiv:1405.5155 (J. Pure Appl. Algebra 220 (2016) —
  **confirm the exact journal/volume at implementation before minting BibTeX**;
  the arXiv id + hypothesis are verified, the journal was not firmly pinned at
  spec time).
- `bv_biklz` → Bian–Itagaki–Kou–Lyu–Zhou, arXiv:2603.04834 (2026) — backs the
  LIVE `k[x]/(x^N)` §3.2 value oracle AND the char-sensitivity dim source; a
  verified `@misc` arXiv entry (add the journal ref if/when it is published).
- Reuse `cyclic` (Connes' B), `bracket`/`gerstenhaber`/`cup` (already present),
  `quantum_ci`/`qci_hh_oracle` for the QuantumCI ties.
- [ ] Commit `docs(citations): BV operator references (Tradler, LZZ, Volkov,
  BIKLZ)`.

---

### Task H — GUI/webapp wiring + verification page + suite gate

**Files:** `hpc/spec.py` (`PRODUCT_KINDS`, `_product_object`, `_snip`,
`_write_product_worked_steps` reachability), `docs/gui/runner.py` (twin +
`_snip` + ETA), `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox,
`S.ids`, push list, renderer), `webapp/static/app.js` (webapp renderer),
`webapp/server/i18n/{en,es}.json`, `trace/products.py` + `trace/results_html.py`
(BV chapter), `tests/webapp/_runner_goldens.json` +
`test_runner_delegation.py` (ONE new fixture, documented),
`webapp/precomputed/*` (Frobenius examples gain the kind), `docs/verification.md`,
`README.md`; tests `tests/webapp/test_bv_block_p54.py`,
`tests/qpa/test_bv_qpa.py`.

- [ ] Follow the Plan-35 §3/§6 + Plan-40 Task-H seven-touchpoint checklist:
  runner parity golden, cross-runner block test, QPA no-surface probe (skip that
  FAILS on change), curated Frobenius examples regenerated (pre-existing blocks
  byte-identical), verification-page P54 row + honest-scope entries (a–d of §7) +
  recounted class table (`tests/release/test_oracle_classes.py` green), README
  line, `tests/webapp/test_js_parses.py` + `node --check` on the touched `.js`.
- [ ] Full gate: `-m pytest tests/hochschild -q` (deep, touched);
  `-m pytest -q -m fast`; `-m pytest tests/qpa -q -m qpa`;
  `-m pytest tests/release -q` — all green; end the chain with `; echo EXIT=$?`.
- [ ] Commit `feat(gui,webapp,hpc): bv_operator compute kind end-to-end +
  verification page + recounted oracle classes`.

---

## 10. Acceptance criteria (definition of done)

1. `A.bv_operator(top, engine=...)` exists, documented, exact, GF(p)-window-routed
   per §5, loud on every refusal path (non-Frobenius, non-semisimple/BIKLZ-blocked,
   off-GF(p), past-window, `engine="cs"`-until-P51).
2. The symmetric/Tradler route is complete and **P52-free**; the semisimple-ν/LZZ
   route lands over P52 twisted homology (or is deferred behind a loud typed
   refusal + ledger entry per §3 if P52 slips), with the twist direction
   **empirically determined and certified** by the perfect-pairing + arbiter
   self-certs.
3. `Δ² = 0`, the seven-term relation `(7T)` (self-cert), and **derived bracket ≡
   independent Gerstenhaber bracket in-window over odd primes** (cross-engine
   arbiter) hold on every served instance; a mismatch (under both twist directions
   on the twisted route) refuses loudly, never a silent wrong Δ. **GF(2) is SERVED**
   (same code path, convention inherited from the odd-prime anchor; certified by
   `Δ²=0` + BV-relation mod 2 + pairing invertibility) with the char-2 honest-scope
   caveat recorded.
4. Literature pins: Tradler `k[x,y]/(x²,y²)` ([4,4,5,6]); `k[x]/(x^N)`
   char-sensitivity in both regimes; **the BIKLZ §3.2 `kZ₁/J^N` explicit-Δ VALUES
   as a LIVE `oracle_literature` assert** (transcribed from `arxiv.org/html/
   2603.04834`, re-verified at implementation time — never fenced, never
   fabricated); QuantumCI ties to the Plan-35 pins (P52-gated, Task E).
5. `bv_operator` clickable end-to-end (GUI canvas → block → report) in EN+ES, both
   runners byte-identical, ONE golden added with a documented change-log entry,
   and **the canonical key of every request that does NOT gain the kind stays
   byte-identical**; the curated Frobenius examples that DO gain `bv_operator`
   change their keys BY DESIGN (Plan-25/Plan-35 precedent — regenerated with every
   pre-existing block byte-identical, the manifest noting the addition).
6. QPA no-BV-surface honest-scope entry with the FAIL-on-change probe; the theory
   oracle covering the gap named.
7. `docs/verification.md` recounted (new oracle rows + honest-scope a–d +
   `test_oracle_classes.py` green); README line added; deep (touched) + fast + qpa
   + release suites green on the merge commit; no golden / cache key / frozen pin
   moved except the documented runner-goldens addition and the Frobenius curated
   regeneration.

---

## Methodology & assumptions

**Approach.** I read Plan 35 end-to-end (the surface P54 extends — `HHProducts`/
`ConnesB`/`ProductTable`, the `.blocks()` single-shape serialization, the
basis-dependence doctrine, the GF(p) bar/tt vs CS routing, the sign-arbitration
precedent) and skimmed Plan 40 for TDD task granularity (failing test → verify
failure → implement with an "adjust to reality" note → pass → commit, per task).
I studied the live code the plan builds on: `hochschild/products.py`
(`connes_b_tables`, `gfp_product_tables`, `_class_coords`), `hochschild/cyclic.py`
(the generic `(b,B)` Connes B), `invariants/frobenius.py`
(`nakayama_automorphism_generic = G⁻¹Gᵀ`, the trace-form symmetry certifier,
socle criterion), `core/algebra.py` (the product/HH dispatch signatures and the
`nakayama_automorphism`/`is_frobenius`/`is_symmetric` public methods), and the
GUI/webapp product wiring (`hpc/spec.py::PRODUCT_KINDS`/`_product_object`,
`docs/gui/runner.py`, `gui.js`). I re-verified all four references via
WebSearch/WebFetch and recorded the outcome in §0. I worked out the transport
mathematics (the perfect pairing `(†)`, the adjoint definition `(‡)`, the degree
bookkeeping giving `Δ: HHⁿ→HHⁿ⁻¹`, the symmetric collapse, the BV relation) and
the decidability of each hypothesis, deriving the containment
`Volkov ⟺ (semisimple + finite order) ⊆ LZZ (semisimple)` that makes the
"deliberate pick" clean.

**Assumptions (explicit).**
1. **P52's twisted-homology API is not yet fixed.** I designed the twisted route
   against an *adapter* (`bv/twist.py`) and flagged the exact consumption contract
   + a complete symmetric-only fallback (§3) so P54 is not blocked by P52's final
   shape. This is the alignment point the task asked me to flag; P52 does not exist
   on disk yet.
2. **The perfect pairing is `HHⁿ(A) ≅ D(HH_n(A,{}_1A_ν))`** (LZZ/Tradler duality).
   I assumed the standard degree-matching Frobenius duality. The adjoint
   bookkeeping of `(‡)` gives `Δ_n = (P_{n-1}ᵀ)⁻¹·Bᵀ·P_nᵀ` (derived in §Task B,
   not guessed); any residual Koszul sign is arbiter-pinned. The twist *direction*
   (ν vs ν⁻¹) is **determined by the bracket ARBITER**, NOT by the dim-match:
   post-review I corrected the earlier claim — ν and ν⁻¹ give equal twisted-
   homology dims in general, so the dim-match is only a coarse coefficient-sanity
   check (a one-time pin verifies this), and direction validation is P52-gated
   (the symmetric fallback never exercises it). This is the Plan-21 arbitration
   precedent applied honestly.
3. **`nu_is_semisimple` = minpoly squarefree** is exact over perfect fields
   (`GF(p)`, `QQ`); I assumed the Domains in scope are perfect (they are) and noted
   the loud-refuse fallback for any Domain where sympy's minpoly is unavailable.
4. **`kZ₁/J^N = k[x]/(x^N)` is symmetric in every characteristic** (commutative
   Frobenius ⇒ ν = id), so its char-sensitivity is served by the symmetric route
   and does NOT depend on the blocked BIKLZ pins — a load-bearing simplification I
   verified from the definition of the Nakayama automorphism.
5. **GF(2) cannot arbitrate the BV signs** (characteristic-2 collapses `(−1)`); I
   assumed the arbiter must run over odd primes, and mandated a recorded,
   test-enforced GF(2) exclusion.

**What I deliberately did NOT check.**
- **BIKLZ (arXiv:2603.04834) — corrected post-review.** My first pass wrongly
  called the full text inaccessible/lossy and BLOCKED the pins. The adversarial
  critic (and a re-fetch I ran) confirmed a clean HTML full text at
  `arxiv.org/html/2603.04834`, §3.2, with the e = 1 BV/bracket formulas printed
  explicitly. They are now transcribed into §0 as a draft and drive a **LIVE**
  `oracle_literature` assert on the symmetric route; the plan no longer claims the
  text is unavailable. I did **not** independently re-derive the formulas from
  first principles — the implementer re-verifies the transcription against the
  HTML at implementation time and copies the section/label identifiers into the
  test docstring; no value the HTML does not print is ever asserted.
- I did **not** pin the exact journal/volume for Volkov (1405.5155) or re-confirm
  LZZ's volume/pages beyond "J. Algebra 2016" — the arXiv ids + hypotheses are
  firmly verified; I flagged the BibTeX metadata for confirmation at citation
  time (Task G) rather than assert an unverified volume.
- I did **not** run any code or write to `src/`/`tests/` (deliverable is the plan
  only, per the task); the interface signatures and file line-references are from
  reading, and each task carries an "adjust to reality" step so the implementer
  reconciles against the tree.
- I did **not** design the off-GF(p)/past-window route (P51 enhancer) beyond
  scoping it as an optional non-goal, since the arbiter that certifies correctness
  is itself GF(p)-window-bounded in the current tree.

**Why I believe the result is correct.** The design reduces the BV operator to a
*transport already grounded in the codebase*: the pairing is the classical
Frobenius duality (Tradler/LZZ), the Connes B is the existing engine, and the ν
gate is an exact squarefree test over perfect fields. The two places where the
literature could mislead (the twist direction and the BV sign) are not trusted —
they are pinned per instance by self-certifications (perfect-pairing invertibility,
`Δ²=0`, the seven-term relation) and, decisively, by the *independent*
Gerstenhaber bracket that quiverlab already computes with no Frobenius input;
matching it in-window over odd primes both fixes the sign convention (applied
uniformly, GF(2) included) and, on `QuantumCI(q=2)`, discriminates the twist
direction — a genuine correctness proof, not an internal consistency check. The
scope is honestly bounded (GF(p), in-window, semisimple ν or symmetric), every
boundary is a loud typed refusal, and the one exact-Δ-VALUE literature oracle (the
BIKLZ §3.2 `k[x]/(x^N)` pin) is transcribed from the verified HTML full text and
re-verified at implementation time — never fabricated, and no longer wrongly
called unavailable. Post-adversarial-review, 4 majors + 9 minors are applied
(header change-log).
