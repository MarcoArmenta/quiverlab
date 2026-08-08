# Plan 53: Homological Invariants II — φdim/ψdim as algebra invariants + fractional Calabi–Yau dimension (R23 + R24)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Lift the module-level Igusa–Todorov surface that Plan 40 already ships to
the **algebra level**, and add the **stable-category fractional Calabi–Yau
dimension** of self-injective algebras — both exact, both honesty-gated. Concretely:
(a) `phi_dim(A)` / `psi_dim(A)` as algebra invariants (a certified exact value for
representation-finite `A` via the ⊕-of-all-indecomposables theorem; a certified
LOWER bound otherwise — never a claimed sup); (b) the **φ-spectrum** `{φ(M) : M
indecomposable}` with **gap** reporting (Barrios–Mata–Rama); (c) the standing chain
self-certificate `findim ≤ φdim ≤ ψdim ≤ gldim`, checked only on the terms the
bounded engine actually resolved; (d) **Lat-Igusa-Todorov (LIT) finitistic
certificates** — four decidable families each emitting a *proof-carrying* finite
`findim` upper bound (self-injective, Iwanaga–Gorenstein, finite φdim, finite
`id(A_A)`), which flip Plan 40's honest `None` upper bound exactly where a certified
bound exists; and (e) the **fractional Calabi–Yau dimension** `(m, ℓ)` of a
self-injective algebra `A` in its stable module category `mod̲ A`, certified by the
Ivanov–Volkov / Erdmann–Skowroński criterion `Ω^{n+1} ≅ ν^{-1}` (and its `ℓ ≥ 2`
generalisation `S^ℓ ≅ Σ^m`), object-wise on a generating set with a bounded search +
loud budget refusal. Two no-code GUI touches: the existing `homological_profile`
block gains φdim/ψdim/spectrum/LIT rows, and a new `fractional_cy` compute kind puts
the stable CY dimension one click away.

**This plan EXTENDS Plan 40** (`docs/plans/2026-08-05-plan-40-homdims.md`,
`src/quiverlab/modules/homdims.py`). Everything Plan 40 shipped is a hard dependency
and is reused byte-unchanged: `igusa_todorov_phi`/`_psi` (last-strict-drop φ with the
Fitting-closure termination window; ψ = φ + fpd), the `φ = ψ = pd` identity for
finite projective dimension, the self-injective `⇒ φ ≡ 0` short-circuit (returns `0`
instantly — this plan's φdim self-injective route is `0` by the same theorem), and
`finitistic_dimension_bounds`'s honest `upper=None` degrade (this plan's LIT
certificate is what supplies a certified `upper` where one exists). Nothing in
`homdims.py`'s existing public functions changes signature; the new work is additive.

**Architecture:** Two thin exact layers over primitives that already exist:

- **Extend `src/quiverlab/modules/homdims.py`** with the algebra-level invariants
  (`phi_dim`, `psi_dim`, `phi_spectrum`, `lit_finitistic_certificate`) and their
  honesty dataclasses (`PhiDim`/`PsiDim`/`PhiSpectrum`/`LITCertificate`), mirroring
  the `GlobalDimension`/`DominantDimension` pattern already in the file. The
  rep-finite exhaustive route consumes P41 `knit_ar_quiver` (the finite
  indecomposable universe, `.is_complete`/`.status` budget contract, each vertex's
  `["module"]`), exactly as Plan 49's `degeneration_order` does; add-monotonicity of
  φ and ψ makes φdim/ψdim a **single** `igusa_todorov_phi`/`_psi` call on the ⊕ of
  all indecomposables. The chain self-cert and the LIT families compose the shipped
  `global_dimension`, `gorenstein_dimension`, `is_selfinjective`, and
  `finitistic_dimension_bounds`.
- **Create `src/quiverlab/modules/fractional_cy.py`** — the stable-category
  fractional CY certificate `fractional_calabi_yau(A, ...)` + the `FractionalCY`
  dataclass. It leans ENTIRELY on module operations that already ship: the Nakayama
  functor `ν` = `modules/ar.py::nakayama_functor` / `nakayama_functor_minus`
  (`Module.nakayama()`/`.nakayama_minus()`), the syzygy/cosyzygy `Ω`/`Ω^{-1}` =
  `modules/resolution.py::syzygy`/`cosyzygy` (`Module.syzygy()`/`.cosyzygy()`), the
  stable normal form via `modules/decompose.py::decompose` (drop projective
  summands), and `modules/hom.py::is_isomorphic` (the loud-when-undecidable
  certificate — the SAME primitive Plan 40's `omega_periodicity` uses). No new math
  engine; no Hochschild `A^e` resolution is needed for the module-level certificate
  (the bimodule-syzygy functorial upgrade is a named successor — see the honest-scope
  note in Task D).
- Thin `Algebra` delegates in `src/quiverlab/core/algebra.py` beside
  `global_dimension`/`is_selfinjective`: `A.phi_dim()`, `A.psi_dim()`,
  `A.phi_spectrum()`, `A.finitistic_certificate()`,
  `A.fractional_calabi_yau_dimension()` / `A.is_fractionally_calabi_yau()`.

**Tech Stack:** pure exact linear algebra over `Domain` (`modules/hom.py`,
`modules/linalg_mod`); **sympy integer rank** for the Igusa–Todorov K₀ (inherited
from Plan 40, a ℤ-computation independent of the field); P41 `knit_ar_quiver` for the
rep-finite indecomposable universe. No floats in `src/` (AST-gated by
`tests/test_no_floats.py`).

---

## The Records (verbatim — the adjudicated spec layer, copied from
`docs/plans/2026-08-06-computability-expansion-deep-research.md`)

**R23 — φdim/ψdim as ALGEBRA invariants + the φ-spectrum + LIT certificates.**
[B-scout P1+P2, MAJOR re-scope: quiverlab ALREADY SHIPS φ(M)/ψ(M), the truncated
closed forms, φ=pd, and the self-injective ⇔ φ≡0 oracles (Plan 40)] Surviving new
objects: (a) φdim(A)/ψdim(A) — rep-finite exhaustive sup over the knitted
indecomposables (monotonicity on add makes ⊕-of-all sufficient); certified lower
bounds otherwise; the chain findim ≤ φdim ≤ ψdim ≤ gldim as a standing self-cert;
(b) the φ-spectrum/gaps (Barrios–Mata–Rama arXiv:1810.12112); (c) LIT-algebra
certificates: no known decision procedure in general (NOT 'proven undecidable' —
reworded per critic); the decidable certificate families (φdim ≤ 1, Gorenstein with
𝒟 = Gproj, self-injective, id(A_A) < ∞) emit a proof-carrying findim bound
ψ_𝒟(V) + n + 1. Corrected refs: Fernandes–Lanzilotta–Mendoza arXiv:1304.0754 (φdim
paper — NOT Huard); Barrios–Lanzilotta–Mata survey arXiv:2310.09283; LIT definitions
in Bravo–Lanzilotta–Mendoza–Vivero arXiv:2002.07866 (JPAA) — arXiv:2105.06273 is
Barrios–Mata 'On Lat-Igusa-Todorov algebras'; arXiv:2103.12120, 2311.06148. Size S–M
on top of Plan 40.

**R24 — Fractional Calabi–Yau dimension of self-injective algebras.** [D-critic
omission-find; new record] Object: the (m, n) with Ω^{m}-shifted Nakayama-twist
periodicity certifying A fractionally CY of dimension m/n in the stable category;
every prerequisite shipped (syzygy periods, ν, Frobenius certifiers, Π(A_n) known
fractional CY values as oracles). Refs to be pinned at plan time from the
Keller-school literature (e.g. the standard fractionally-CY computations for
preprojective and Nakayama algebras). Size S–M.

---

## Reference re-verification (mandatory, done at spec time — 2026-08-07)

All citations below were re-verified via web search/fetch during authoring. Findings:

1. **arXiv:1304.0754 — CONFIRMED authorship correction.** "The Φ-dimension: A new
   homological measure," by **Sonia Fernandes, Marcelo Lanzilotta, Octavio Mendoza**
   (v1 2013, v2 2014; published Algebr. Represent. Theory 18 (2015) 463–476,
   `10.1007/s10468-014-9504-9`). **NOT Huard** — the critic was right. It defines
   `φdim(A) = sup{φ(M)}`, characterises it via the Ext/Tor bifunctors, and proves
   finiteness of φdim is invariant under derived equivalence.
2. **arXiv:1810.12112 — CONFIRMED.** "Gaps for the Igusa-Todorov function," by
   **Marcos Barrios, Gustavo Mata, Gustavo Rama** (2018). Theorem: if
   `0 < φdim(A) = m < ∞` then there exist modules `M, N` with `φ(M) = m − 1` and
   `φ(N) = 1` (so `1` and `m − 1` are always in the φ-spectrum). A **gap** is an
   integer in `[1, φdim]` not attained by φ; algebras with gaps satisfy the
   finitistic dimension conjecture. (Already in the registry as `barrios_mata` →
   `BarriosMataRama2020`.)
3. **arXiv:2002.07866 — CONFIRMED (LIT definitions).** "Generalised Igusa-Todorov
   functions and Lat-Igusa-Todorov algebras," by **Diego Bravo, Marcelo Lanzilotta,
   Octavio Mendoza, José Vivero** (2020; J. Algebra 2021). Introduces LIT algebras — a
   class satisfying the finitistic dimension conjecture that INCLUDES both the
   Igusa–Todorov algebras (Wei) AND the self-injective algebras. The proof-carrying
   findim bound is the record's `ψ_𝒟(V) + n + 1` form (copied verbatim into Task C;
   the exact statement of the relative function `ψ_𝒟` and the four families' data is
   pinned from this paper at implementation — see Task C's `# PIN` note).
4. **arXiv:2310.09283 — CONFIRMED (survey).** "A survey on Igusa-Todorov functions"
   (Barrios–Lanzilotta–Mata). Explicitly states the **chain
   `findim(A) ≤ φdim(A) ≤ ψdim(A) ≤ gldim(A)`**, the **add-monotonicity of BOTH φ and
   ψ** (`add M ⊆ add N ⇒ φ(M) ≤ φ(N)` and `ψ(M) ≤ ψ(N)`), and that for a
   representation-finite algebra **`φdim(A) = φ(M₀)` with `M₀` = the direct sum of
   all indecomposables** — the three facts Task A is built on.
5. **arXiv:2105.06273 / 2103.12120 / 2311.06148 — CONFIRMED** as Barrios–Mata "On
   Lat-Igusa-Todorov algebras," "Triangular Lat-Igusa-Todorov algebras," and
   "Generalised Lat-Igusa-Todorov Algebras and Morita Contexts" respectively — the
   LIT extension papers named in R23.

6. **Fractional CY — the criterion is PINNED, not left open.** The authoritative
   source is **Ivanov–Volkov, arXiv:1212.2619, "Stable Calabi–Yau dimension of
   self-injective algebras of finite type"** (§1.3, transcribed verbatim below), which
   completes the program of **Erdmann–Skowroński** (who introduced the stable CY
   dimension of a self-injective algebra as the weak CY dimension of `mod̲ A`):
   - `mod̲ A` (stable module category of a self-injective `A`) is triangulated with
     **suspension `Σ = Ω^{-1}`** (Heller's cosyzygy) and **Serre functor `S = Ω ν`**
     (`ν` = the Nakayama functor `D Hom_A(−, A)`).
   - `ν ≅ id` iff `A` is symmetric; `ν` is an autoequivalence iff `A` is
     self-injective (so `ν^{-1}` = `nakayama_functor_minus` exists exactly in scope).
   - **THREE TIERS of the Calabi–Yau condition — kept distinct EVERYWHERE below (the
     definitional correction from adversarial review):**
     1. **strong `n`-CY** — `S ≅ Σ^n` **as functors** (a natural isomorphism);
        equivalently `Ω^{n+1} ≅ ν^{-1}` **as functors**. This is Keller's tier and the
        phrase "as functors" is RESERVED for it.
     2. **weak `n`-CY** — `S X ≅ Σ^n X` for **all** objects `X` (object-wise on the
        whole category, no naturality); equivalently `Ω^{n+1} X ≅ ν^{-1} X` for all
        `X`. This is Erdmann–Skowroński's / Ivanov–Volkov's DEFINITION of "weakly
        `n`-CY" (their eq. 0.1 is object-wise) — NOT a functor isomorphism. The
        **stable CY dimension** is the least `n ≥ 0` for which `mod̲ A` is weakly
        `n`-CY, or the algebra is only *fractionally* CY (least `ℓ ≥ 1` with
        `S^ℓ X ≅ Σ^m X` for all `X`, some `m`, dimension `m/ℓ`). The "Ω^m-shifted
        Nakayama-twist periodicity" of R24 is exactly `ν^ℓ X ≅ Ω^{-(m+ℓ)} X`
        (⟺ `Ω^{m+ℓ} ν^ℓ X ≅ X`).
     3. **weak-on-generators `(m,ℓ)`-CY (THE TIER THIS PLAN CERTIFIES)** — `S^ℓ X ≅
        Σ^m X` for all `X` in a fixed generating set (the simples + the `Ω^i ν^j`-orbit
        reps met in the search), object-wise. This is a NECESSARY condition for tier 2
        (and hence tier 1); it is not sufficient for either (see Task D honest-scope).
        Every payload field, docstring, and verification-page row labels the shipped
        result as this third tier.
   - Ivanov–Volkov **Theorem 1.8** gives the tier-1 (strong, functorial) bimodule
     upgrade `Ω^{n+1}_{Aᵉ}(A) ≅ (A^∨)_φ` for a stably-inner `φ` — recorded as the named
     successor certificate in Task D (uses quiverlab's minimal `A^e` engine).
7. **Fractional CY oracle VALUES (pinned):**
   - **`k[x]/(x^a)` (symmetric local Nakayama), `a ≥ 3`: stable CY dimension `1`,
     `(m, ℓ) = (1, 1)`.** DERIVED from the criterion: symmetric ⇒ `ν = id` ⇒
     `ν^{-1} = id`; `Ω(k[x]/(x^i)) = k[x]/(x^{a−i})` so `Ω² ≅ id` (a ≥ 3); least
     `n ≥ 0` with `Ω^{n+1} ≅ id` is `n = 1`. In-engine verifiable
     (`Ω²(S) ≅ S`, `ν = id`). For `a = 2` (dual numbers) `Ω ≅ id`, the shift is
     trivial, `(m, ℓ) = (0, 1)` (degenerate — reported with a note). This is a THEORY
     oracle (self-derived from the pinned Serre-functor formula).
   - **Preprojective algebra `Π(Δ)` of Dynkin type (self-injective): stable category
     is 2-Calabi-Yau, `(m, ℓ) = (2, 1)`** (Geiß–Leclerc–Schröer; Erdmann–Snashall
     periodicity). Pinned as a LITERATURE oracle with an in-engine cross-check
     (`Ω³(S_v) ≅ ν^{-1}(S_v)`, i.e. `ν(Ω³ S_v) ≅ S_v`) on `Π("A3")` / `Π("A4")` (the
     rep-finite, hence periodic, small preprojectives — `Π(A_n)` is rep-finite iff
     `n ≤ 4`, so periodic there; `PreprojectiveAlgebra("A4")` builds, dim 20, and is
     self-injective per Plan 33). Implementation MUST verify in-engine and pick the
     smallest instance that is rep-finite AND whose simples are `S`-periodic; if a
     chosen `n` proves non-periodic the search caps loudly (never a wrong value).
   - **`BLOCKED-until-transcribed`:** Ivanov–Volkov **Table 1** gives the full stable
     CY dimensions of standard/nonstandard self-injective algebras of finite type by
     Asashiba type `(Δ, f, t)` (e.g. `(A_n, r/n, 1)`, `(D_{3n}, r/3, 1)`,
     `(A_{2n+1}, r, 2)`, and `∞` rows where `S` is not itself a shift). These are NOT
     wired as oracles in this plan: mapping an Asashiba type triple to a quiverlab
     presentation is non-trivial and is deferred. The plan's fractional-CY acceptance
     rests on the two pinned values above + the internal-consistency oracle (the SAME
     `(m, ℓ)` certifies every simple simultaneously) + the ℓ=1 ≡ Ivanov–Volkov
     integer-criterion cross-check. Transcribing a Table-1 row is a good successor
     oracle and is listed on the verification page as honest scope.

**Sources.** [1304.0754](https://arxiv.org/abs/1304.0754) ·
[1810.12112](https://arxiv.org/abs/1810.12112) ·
[2002.07866](https://arxiv.org/abs/2002.07866) ·
[2310.09283](https://arxiv.org/abs/2310.09283) ·
[2105.06273](https://arxiv.org/abs/2105.06273) ·
[1212.2619](https://arxiv.org/abs/1212.2619) ·
Erdmann–Skowroński (stable CY dimension of tame symmetric algebras, J. Math. Soc.
Japan 58 (2006)) · Geiß–Leclerc–Schröer, "Rigid modules over preprojective
algebras," Invent. Math. 165 (2006).

---

## Mathematical foundation (the plan's ground truth)

### φdim / ψdim as algebra invariants

Per the survey (2310.09283) and Fernandes–Lanzilotta–Mendoza (1304.0754):

- `φdim(A) = sup{ φ(M) : M ∈ mod A }`, `ψdim(A) = sup{ ψ(M) : M ∈ mod A }`.
- **Add-monotonicity:** `add M ⊆ add N ⇒ φ(M) ≤ φ(N)` and `ψ(M) ≤ ψ(N)`.
- **Rep-finite reduction:** if `A` is representation-finite with `M₀ = ⊕_i X_i` the
  direct sum of a representative of every indecomposable, then every `M` has
  `add M ⊆ add M₀`, so `φ(M) ≤ φ(M₀)` and `ψ(M) ≤ ψ(M₀)`; hence
  `φdim(A) = φ(M₀)` and `ψdim(A) = ψ(M₀)` — a SINGLE φ/ψ evaluation.
- **Chain:** `findim(A) ≤ φdim(A) ≤ ψdim(A) ≤ gldim(A)`. When `gldim(A) < ∞` all four
  are EQUAL (every module has finite pd, so `φ = ψ = pd`). When `A` is self-injective
  `φdim = ψdim = 0` while `gldim = ∞` (Plan 40's φ ≡ 0 theorem).
- **Honest computation.** `φdim`/`ψdim` are EXACT only when the indecomposable
  universe is finite and closed (`knit_ar_quiver` returns `is_complete`). Otherwise
  they are a **certified lower bound** `≥ φ(V)` / `≥ ψ(V)` for any finite subfamily
  `V ⊆ mod A` (a rigorous lower bound because `φ`/`ψ` are add-monotone) — NEVER a
  claimed sup. Self-injective short-circuits to the exact `0` (no knit needed — and
  `knit_ar_quiver` REFUSES self-injective anyway, `status="unsupported"`).

### The φ-spectrum and gaps (1810.12112)

For representation-finite `A`, the **φ-spectrum** is `Spec_φ(A) = { φ(X) : X
indecomposable }` (a subset of `{0, …, φdim(A)}`). A **gap** is an integer
`g ∈ [1, φdim(A)]` with `g ∉ Spec_φ(A)`. Barrios–Mata–Rama: if `0 < φdim(A) = m < ∞`
then `1 ∈ Spec_φ(A)` and `m − 1 ∈ Spec_φ(A)` (self-cert oracle). Rep-infinite input
gives an honest partial spectrum labelled a lower bound (never a claimed complete
spectrum).

### LIT finitistic certificate (2002.07866; the record's four families)

An algebra is **Lat-Igusa-Todorov** if it carries data `(𝒟, n)` (a summand-closed,
syzygy-behaved subcategory `𝒟` and `n ≥ 0`) yielding the proof-carrying bound
`findim(A) ≤ ψ_𝒟(V) + n + 1` (record R23c; `V` = generators of `𝒟`, `ψ_𝒟` = the
generalised IT function relative to `𝒟`). This plan does not decide LIT-ness in
general (**"no known decision procedure in general"** — explicitly NOT "undecidable").
It ships the FOUR decidable families, each producing a CERTIFIED finite findim upper
bound from data quiverlab already computes:

1. **Self-injective** (`is_selfinjective`) — EXECUTING TEST: `𝒟 = mod A`, `n = 0`;
   only projectives have finite pd, so **`findim(A) = 0`**. Tested on
   `truncated_polynomial(4)`.
2. **Iwanaga–Gorenstein** (Plan 40 `gorenstein_dimension` both-sided finite) —
   EXECUTING TEST: `𝒟 = Gproj(A)`; **`findim(A) = id(A_A) = id(_AA)` = the Gorenstein
   dimension** (a theorem for Gorenstein algebras — the sup of finite pd's is the
   common injective dimension of the regular module). The bound is
   `max(gd.right_id, gd.left_id)`. **Test algebra: `kA₃/J²` =
   `TruncatedPathAlgebra("A3", 2)`** — NOT self-injective (so family 1 does not
   preempt), and **Iwanaga–Gorenstein because it has finite global dimension** (a
   directed/acyclic quiver ⇒ `gl.dim < ∞` ⇒ every module has finite injective
   dimension ⇒ Iwanaga–Gorenstein). For finite global dimension the three quantities
   coincide, `findim = gl.dim = id(A_A) = Gorenstein dimension`, so the family-2 bound
   equals `int(A.global_dimension())` (`= 2` for `kA₃/J²`) — the test asserts
   `c.family == "gorenstein"` and `c.findim_upper == int(A.global_dimension())`, which
   both confirms the branch executed AND is the decisive family (returned before
   family 3). A stronger *infinite*-gl.dim Gorenstein instance (to isolate family 2
   from family 3 rather than merely order-preempt it) needs a Gorenstein,
   non-self-injective, `gl.dim = ∞` presentation, `# PIN`'d as a successor oracle
   (`kQ/I` Nakayama/one-point-extension examples — not built from a shipped family
   yet).
3. **Finite φdim** (Task A exact, i.e. rep-finite): the chain gives
   **`findim(A) ≤ φdim(A)`** — a certified finite bound (this is the `φdim ≤ 1`
   family generalised to any finite computed φdim). **Coverage note:** on any
   *finite*-gl.dim rep-finite algebra family 2 preempts family 3 (both apply, order
   1→2→3), so family 3 is *decisive* only on a rep-finite, **infinite**-gl.dim,
   non-self-injective, non-Gorenstein input — which no shipped builder cleanly
   certifies (Nakayama algebras are either cyclic ⇒ self-injective, or linear ⇒ finite
   gl.dim). Family 3's inequality `findim ≤ φdim` is therefore covered by the Task A
   **chain self-cert** (it IS family 3's bound) plus a LIT bound-consistency test
   (`findim_lower ≤ c.findim_upper ≤ gl.dim`); a *decisive* family-3 LIT test is
   `# PIN`'d pending such an instance.
4. **Finite `id(A_A)` one-sided (DEMOTED — machinery + `# PIN`, NO executing test)** —
   the LIT bound with `n = id(A_A)` is **`findim upper = the record's ψ_𝒟(V) + n + 1`**
   with `𝒟` the appropriate subcategory; `# PIN` the exact `ψ_𝒟`/`V` from 2002.07866
   at implementation and ARBITRATE the `+ n + 1` constant against the family-1/2 exact
   values (the bound must not undershoot them); HONEST-DEGRADE to `None` if the exact
   statement cannot be certified in-session (never a folklore number — the Plan 40
   rule). **Why no executing test:** family 4 fires only when families 1–3 all fail
   yet `id(A_A) < ∞` — i.e. `id(A_A)` finite but `id(_AA)` **infinite** (else family 2
   fires). But the **bounded engine never PROVES an injective dimension infinite** (it
   only ever returns a lower bound), so it cannot certify the "one-sided-finite,
   not Gorenstein" precondition; `gorenstein_dimension` returns `None` (undecided), not
   a decidable "not Gorenstein". Family 4's activation condition is thus not decidably
   reachable with shipped tools. The `ψ_𝒟(V)+n+1` machinery ships and is unit-covered
   by a direct call on constructed `(𝒟, n)` data, but there is NO end-to-end executing
   `lit_finitistic_certificate` test for family 4 (the acceptance criterion and the
   test file reflect exactly this).

`lit_finitistic_certificate(A)` tries families 1→2→3→4 and returns the FIRST that
applies (the sharpest is usually earliest), with a `family` label and the
proof-carrying bound; if none applies it returns `status="no known decision procedure
in general"` and `upper=None` — the Plan 40 degrade, preserved. This function is what
`finitistic_dimension_bounds`'s `upper` consults BEFORE degrading (see Task C wiring):
Plan 40's honest `None` becomes a certified number exactly in these families.

### Fractional Calabi–Yau dimension (stable category, self-injective) — R24

Scope: `A` **self-injective** (loud `QuiverlabError` otherwise — reuse
`is_selfinjective`; off scope `ν` is not an autoequivalence and `mod̲ A` is not
triangulated). In `mod̲ A`: suspension `Σ = Ω^{-1}` (cosyzygy), Serre functor
`S = Ω ν` (Ivanov–Volkov §1.3). The **fractionally-CY notion** is the standard one of
**Herschend–Iyama–Oppermann** / **Keller** ("On triangulated orbit categories"): a
Hom-finite triangulated category `T` with Serre functor `S` is **fractionally CY of
dimension `m/ℓ`** iff `S^ℓ ≅ Σ^m`. For `mod̲ A` this reads `ℓ` = the least positive
integer for which `S^ℓ X ≅ Σ^m X` holds for all objects `X` (some `m ∈ ℤ`). Because
`ν` and `Ω` commute on `mod̲ A`, `S^ℓ = Ω^ℓ ν^ℓ`, so — fixing the sign to match the
`ℓ=1 ≡` Ivanov–Volkov form below — the **certified identity** is, per generator `X`
(stable iso, i.e. iso after dropping projective summands):

```
S^ℓ(X) ≅ Σ^{+m}(X)   ⟺   Ω^ℓ ν^ℓ(X) ≅ Ω^{-m}(X)   ⟺   Ω^{m+ℓ}( ν^ℓ(X) ) ≅ X.
```

The `ℓ = 1` case is exactly Ivanov–Volkov's integer stable CY dimension: least
`n ≥ 0` with `Ω^{n+1} ≅ ν^{-1}` (set `m = n`, `ℓ = 1`: `Ω^{n+1} ν X ≅ X ⟺ Ω^{n+1} X ≅
ν^{-1} X`, object-wise) — so the derived sign is fixed, NOT free (the arbiter test in
Task D CONFIRMS this sign, it does not resolve an open choice).

**The `m` is only defined modulo the Σ-period.** On a periodic stable category `Σ` has
finite order `p` on the generators (`Σ^p X ≅ X`), so `S^ℓ X ≅ Σ^m X` implies
`S^ℓ X ≅ Σ^{m + kp} X` for all `k`: the numerator `m` is a residue mod `p`, a genuine
ambiguity, not a search tie. The plan fixes the **canonical representative** by
"smallest `ℓ ≥ 1`, then smallest `m ≥ 0`" and reports `p` (the Σ-period on the
generators, when the search finds it) in a payload `sigma_period` field so the
normalization is auditable. Reduced rational `m/ℓ` is the CY dimension; when `ℓ = 1`
also report the Ivanov–Volkov integer `n = m`.

**Generating set + honest scope (tier 3 of the criterion pin).** The certificate is
checked object-wise on the **simple modules** (they generate `mod̲ A` as a thick
triangulated subcategory — every module has a finite radical filtration by simples,
giving triangles) together with the `Ω^i ν^j`-orbit representatives encountered during
the search (genuine indecomposables). This is **weak-on-generators `(m,ℓ)`-CY** — a
NECESSARY condition for the tier-2 weak CY property (and hence for tier-1 strong CY),
the SAME certificate style Plan 40's periodicity and the AR/almost-split engine use; it
is ANCHORED by (i) the pinned literature values and (ii) internal consistency (one
`(m, ℓ)` for all generators). **Two ways it can be too optimistic, stated plainly:**
the object-wise-on-generators check can (a) mis-report the numerator `m` and, equally,
(b) **under-report the denominator `ℓ`** — a spurious small `ℓ` may pass on the simples
while failing on a non-generator object, since object-wise agreement on a generating
set need not propagate to the whole category. The FULLY functorial (tier-1) certificate
— `Ω^{n+1}_{Aᵉ}(A) ≅ (A^∨)_φ` (Ivanov–Volkov Thm 1.8, via the minimal `A^e` engine
`engine/resolutions_minimal.py`) — is the named successor recorded on the verification
page; this plan does NOT claim a functor isomorphism, only the certified
weak-on-generators witness. **Completeness:** a finite stable CY dimension is found
exactly when the
generators are `S`-periodic, which holds for **finite representation type**
self-injective algebras (every f.d. self-injective of finite type is periodic); a
rep-infinite self-injective input caps at the budget and returns `status="budget"`
(honest — never "not fractionally CY"). The decompose char caveat applies (stable
normal form uses `decompose`): batteries run over `QQ` / `GF(32003)`.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **Plan 40 (C6 homdims), P37 (`direct_sum`), P41 (`knit_ar_quiver`), P30
  (`decompose`), P23/P24 (`Algebra.opposite`, `duality.tau`), and the module-level
  Nakayama functor (`modules/ar.py::nakayama_functor`/`nakayama_functor_minus`,
  `Module.nakayama()`/`.nakayama_minus()`) and `syzygy`/`cosyzygy`
  (`modules/resolution.py`, `Module.syzygy()`/`.cosyzygy()`) are all MERGED to `dev`
  and are hard prerequisites** — verified at authoring. Branch `plan-53-phidim-fcy`
  off `dev` (metaplan Wave 1; R23+R24 are `[independent]`).
- **The decompose char caveat is load-bearing** (the Plan 40 precedent): `decompose`
  / `is_isomorphic` RAISE loudly over `GF(p)` when `char ≤ dim M`. Both the φdim K₀
  bookkeeping AND the fractional-CY stable normal form inherit this loud refusal.
  Every φdim/ψdim/spectrum/fractional-CY battery runs over **`QQ`** or **`GF(32003)`**
  so both routines decide; a small-`GF(p)` test with `char ≤ dim` is a test bug, not
  an engine bug.
- **The honesty pattern is mandatory** (`GlobalDimension`/`DominantDimension`/
  `FinitisticBounds` precedents in `homdims.py`): φdim/ψdim are an exact value OR a
  labelled certified LOWER bound (never a claimed sup); the LIT findim upper is a
  certified number OR an honest `None` (never folklore); the fractional CY is a
  certified `(m, ℓ)` OR `status="budget"`/`"shift_trivial"` (never a fabricated
  "not CY"); off-scope inputs RAISE a typed `QuiverlabError`, never a silent wrong
  answer.
- **`knit_ar_quiver` REFUSES self-injective input** (`status="unsupported"`). The
  φdim rep-finite route inherits that — but self-injective is short-circuited to the
  exact `φdim = ψdim = 0` BEFORE any knit (Plan 40's φ ≡ 0 theorem), so the
  unsupported status is never hit for a legitimate φdim.
- Plan-32 markers: the chain `findim ≤ φdim ≤ ψdim ≤ gldim` (on computed terms), the
  `φdim = ψdim = gldim = findim` collapse when gldim is exact-finite, the
  Barrios–Mata–Rama `{1, m−1} ⊆ Spec_φ` fact, the LIT proof-carrying bounds' internal
  consistency (bound ≥ the exact `findim` where known), and the fractional-CY
  `Ω^{m+ℓ} ν^ℓ ≅ id` self-identity + `ℓ=1 ≡` Ivanov–Volkov integer form are
  `oracle_selfcert`; the rep-finite `φdim` closed values (kA_n, Nakayama), the
  `k[x]/(x^a)` and `Π(Δ)` fractional-CY values are `oracle_literature`; two-route
  agreement (φdim via ⊕-of-all ≡ `max` of φ over the knit's indecomposables;
  fractional-CY ℓ=1 via `Ω^{n+1} ≅ ν^{-1}` ≡ via `S ≅ Σ^n`) is `oracle_crossengine`;
  QPA lives in `tests/qpa/` (bucket = the class, never double-marked).
- Buckets (auto by directory, `tests/conftest.py`): `tests/modules/` → **deep**;
  `tests/invariants/` → **fast**; `tests/qpa/` → **qpa**; `tests/webapp/` extras-gated
  cross-runner tests unmarked (Plan-32 ruling). The φdim/spectrum/LIT/fractional-CY
  batteries go in `tests/modules/` (they consume `knit_ar_quiver`/`decompose` and are
  naturally heavy) = **deep**.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class
  table, `tests/release/test_oracle_classes.py` green) and adds its citations to
  `citations/references.bib` + `registry.py` (the `bibtex()` helper hard-fails if the
  two are out of sync). Conventional commits; green at every commit.
- GUI is in-plan (metaplan §1.2): the `homological_profile` block extension + the new
  `fractional_cy` kind ship in this plan (both runners, schema guard, canonical-key
  stability, i18n ×4, report rendering, citations) — no ledger deferral.

---

### Task A: φdim / ψdim as algebra invariants + the standing chain

**Files:**
- Modify: `src/quiverlab/modules/homdims.py` (add `PhiDim`/`PsiDim` dataclasses,
  `phi_dim`, `psi_dim`, `_all_indecomposables_or_bound`, `_chain_selfcheck`)
- Modify: `src/quiverlab/core/algebra.py` (`Algebra.phi_dim`, `Algebra.psi_dim`)
- Test: `tests/modules/test_phidim.py`

**Interfaces:**
- Consumes: `igusa_todorov_phi`/`igusa_todorov_psi` (Plan 40, byte-unchanged),
  `is_selfinjective` (`modules/ext.py`), `knit_ar_quiver`
  (`modules/ar.py`, `.is_complete`/`.status`/`.vertices[i]["module"]`),
  `modules/morphism.py::direct_sum`, `modules/ext.py::global_dimension`,
  `finitistic_dimension_bounds` (Plan 40).
- Produces:
  ```python
  @dataclass
  class PhiDim:
      value: int            # the exact φdim (rep-finite) OR the certified lower bound
      exact: bool           # True iff the indecomposable universe closed (knit complete)
      status: str           # "complete" | "budget" | "self-injective"
      # __int__/__eq__(int) mirroring GlobalDimension; __repr__:
      #   exact -> "φdim = {value}"; else -> ">= {value} (certified lower bound; the
      #   indecomposable universe did not close -- knit status {status})"
  # PsiDim identical shape (ψdim).
  def phi_dim(A, *, budget_modules=256, phi_budget=512, phi_bound=64) -> PhiDim
  def psi_dim(A, *, ...) -> PsiDim
  ```
  `phi_dim`: self-injective ⇒ `PhiDim(0, exact=True, status="self-injective")`
  (Plan 40 φ ≡ 0). Else `knit_ar_quiver(A, budget_modules=...)`; branch on the knit
  outcome — **`status="error"` and `status="budget"` are NOT the same degrade
  (M3, adversarial review):**
  - `is_complete` ⇒ `M0 = direct_sum(*[v["module"] for v in ar.vertices])`,
    `PhiDim(igusa_todorov_phi(M0), exact=True, status="complete")`.
  - `status == "budget"` ⇒ SOFT degrade to the certified LOWER bound
    `igusa_todorov_phi(⊕ discovered-prefix ∪ ⊕_v S_v)` with
    `exact=False, status="budget"` (the discovered indecomposables + the simples are a
    genuine finite subfamily; add-monotonicity ⇒ their φ is `≤ φdim`, a valid lower
    bound).
  - `status == "error"` ⇒ **RAISE loudly** — a genuine AR-knitting failure is a bug
    or an unhandled input, NOT an honest partial answer; re-raise as `QuiverlabError`
    carrying `ar.note` (the offending dimension vector). Never return a `PhiDim` with
    `status="error"` (a lower bound built on a broken knit is not certified). The
    `status="unsupported"` case is unreachable (self-injective short-circuited above);
    if ever hit, treat it as `error` (loud).

  `psi_dim` identical with `igusa_todorov_psi`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_phidim.py
"""phidim/psidim as algebra invariants (Plan 53 / R23a). Literature: gldim finite =>
findim = phidim = psidim = gldim (survey 2310.09283 chain collapse); self-injective =>
phidim = psidim = 0 (Plan 40 phi==0). Self-cert: the standing chain
findim <= phidim <= psidim <= gldim on computed terms; the rep-finite
phidim = phi(+ all indec) route == max phi over indecomposables. Over QQ (decompose
char caveat)."""
import pytest

from quiverlab import GF, Quiver, TruncatedPathAlgebra, linear_path_algebra, \
    truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import phi_dim, psi_dim

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
def test_finite_gldim_collapse():
    # kA3/rad^2: gl.dim 2 => findim = phidim = psidim = gldim = 2 (chain collapse).
    A = TruncatedPathAlgebra("A3", 2, field=QQ)
    assert phi_dim(A) == 2 and phi_dim(A).exact is True
    assert psi_dim(A) == 2

@lit
def test_hereditary_kA3_phidim_is_gldim():
    A = linear_path_algebra(3, field=QQ)                # hereditary, gl.dim 1
    assert phi_dim(A) == 1 and psi_dim(A) == 1

@lit
def test_self_injective_phidim_zero():
    A = truncated_polynomial(4, field=QQ)               # k[x]/(x^4): self-injective
    pd = phi_dim(A)
    assert pd == 0 and pd.exact is True and pd.status == "self-injective"
    assert psi_dim(A) == 0

@xeng
def test_phidim_equals_max_phi_over_indecomposables():
    # the +-of-all route (add-monotonicity) == max phi over the knit indecomposables.
    from quiverlab.modules.ar import knit_ar_quiver
    from quiverlab.modules.homdims import igusa_todorov_phi
    A = TruncatedPathAlgebra("A3", 2, field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    manual = max(igusa_todorov_phi(v["module"]) for v in ar.vertices)
    assert phi_dim(A).value == manual

@selfcert
def test_standing_chain_holds():
    from quiverlab.modules.homdims import _chain_selfcheck    # returns dict of checks
    for A in (linear_path_algebra(3, field=QQ),
              TruncatedPathAlgebra("A3", 2, field=QQ),
              truncated_polynomial(4, field=QQ)):
        chk = _chain_selfcheck(A)
        assert chk["ok"]                                # findim<=phidim<=psidim<=gldim

@selfcert
def test_char_caveat_refuses_loudly():
    from quiverlab.errors import QuiverlabError
    A = TruncatedPathAlgebra("A3", 2, field=GF(2))      # char 2 <= module dims
    with pytest.raises(QuiverlabError):
        phi_dim(A)
```

- [ ] **Step 2: Run to verify failure** — `ImportError: phi_dim`
- [ ] **Step 3: Implement.** Add `PhiDim`/`PsiDim` beside `DominantDimension`; a
  shared `_phipsi_dim(A, fn, ...)` (fn = `igusa_todorov_phi`/`_psi`) with the
  self-injective / knit-complete / lower-bound branches above; `_chain_selfcheck(A)`
  computes `findim` (lower from `finitistic_dimension_bounds`), `phidim`, `psidim`,
  `gldim` (from `global_dimension`) and asserts every inequality **whose both terms
  are exact/computed** (findim's `lower` ≤ φdim.value ≤ ψdim.value; ψdim.value ≤
  gldim.value when `gldim.exact`; and when `gldim.exact` all four are EQUAL). Wire the
  two `Algebra` delegates.
- **Adjust to reality:** confirm `direct_sum(*mods)` over the whole indecomposable
  set is within budget on the oracle inputs (kA_n small); the lower-bound branch must
  DEDUPLICATE the discovered prefix ∪ simples by `is_isomorphic` before `direct_sum`
  (a repeated summand does not change φ but wastes the decompose). `knit_ar_quiver`
  self-injective returns `status="unsupported"` — unreachable here (short-circuited).
- **Cost note (W4, adversarial review).** `igusa_todorov_phi(M0)` opens with
  `decompose(M0)`, and `M0` is the ⊕ of every indecomposable — so its dimension can
  reach the AR-knit's `budget_dim` (4096) across up to `budget_modules` (~256)
  summands; `decompose` on a module that large is the dominant cost. **The knit budget
  IS the guard** (a bigger indecomposable universe trips `status="budget"` and never
  reaches the `decompose`), so no separate cap is added — but the oracle inputs must
  stay small. Add a **wall-clock smoke bound** on the largest cited φdim oracle
  (`TruncatedPathAlgebra("A4", 2)` in Task B, ~10 indecomposables, `M0` dim ≈ 20): the
  `deep` test must complete `phi_dim` on it in `< 30 s` on the CI cell (assert nothing
  — a `pytest` `-k` timing note in the module docstring; if it exceeds, shrink the
  oracle to `"A3"`). Do NOT cite a φdim oracle whose `M0` approaches `budget_dim`.
- [ ] **Step 4: Run tests** — PASS.
- [ ] **Step 5: Commit** — `feat(modules): phidim/psidim as algebra invariants (rep-finite +-of-all via add-monotonicity; certified lower bounds; standing findim<=phidim<=psidim<=gldim chain)`

---

### Task B: the φ-spectrum + gaps

**Files:**
- Modify: `src/quiverlab/modules/homdims.py` (`PhiSpectrum` dataclass, `phi_spectrum`)
- Modify: `src/quiverlab/core/algebra.py` (`Algebra.phi_spectrum`)
- Test: `tests/modules/test_phi_spectrum.py`

**Interfaces:**
- Consumes: `knit_ar_quiver`, `igusa_todorov_phi`, `phi_dim` (Task A).
- Produces:
  ```python
  @dataclass
  class PhiSpectrum:
      values: list[int]     # sorted distinct phi over the (knit) indecomposables
      gaps: list[int]       # integers in [1, phidim] not in `values` (rep-finite only)
      phidim: int
      complete: bool        # True iff the knit closed (a genuine spectrum)
      status: str
  def phi_spectrum(A, *, budget_modules=256) -> PhiSpectrum
  ```
  Rep-finite: `values = sorted({ igusa_todorov_phi(X) : X ∈ knit vertices })`,
  `phidim = max(values)`, `gaps = [g for g in range(1, phidim) if g not in values]`
  (`0` and `phidim` are always attained; Barrios–Mata–Rama guarantee `1` and
  `phidim−1` are too). Not closed ⇒ `complete=False`, `gaps=[]` (a partial spectrum
  claims no gaps — a "gap" needs a proven complete spectrum).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_phi_spectrum.py
"""phi-spectrum + gaps (Plan 53 / R23b; Barrios-Mata-Rama 1810.12112). Self-cert:
0 and phidim are attained; when 0 < phidim < infinity, 1 and phidim-1 are attained
(BMR theorem). Rep-finite only; a partial spectrum claims no gaps. Over QQ."""
import pytest

from quiverlab import TruncatedPathAlgebra, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.homdims import phi_spectrum

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
def test_bmr_endpoints_attained():
    A = TruncatedPathAlgebra("A4", 2, field=QQ)          # phidim = 3
    S = phi_spectrum(A)
    assert S.complete
    assert 0 in S.values and S.phidim in S.values
    if 0 < S.phidim < 10**6:
        assert 1 in S.values and (S.phidim - 1) in S.values   # BMR

@lit
def test_hereditary_spectrum_is_zero_one():
    A = linear_path_algebra(3, field=QQ)                 # hereditary: phi = pd in {0,1}
    S = phi_spectrum(A)
    assert S.values == [0, 1] and S.gaps == [] and S.phidim == 1
```

- [ ] **Step 2: Run to verify failure**
- [ ] **Step 3: Implement** `phi_spectrum` + `PhiSpectrum`; wire the delegate.
  **Adjust to reality:** the knit vertices are ALREADY one-per-indecomposable, so no
  dedup needed for `values`; the BMR endpoints are an assertion in the TEST (the
  self-cert), not enforced in the library (the library just reports what it found — if
  BMR failed the test would catch an engine bug).
- [ ] **Step 4: Run tests** — PASS.
- [ ] **Step 5: Commit** — `feat(modules): phi_spectrum + gaps (Barrios-Mata-Rama; rep-finite, partial spectrum claims no gaps)`

---

### Task C: LIT finitistic certificate (four families) — flips Plan 40's honest None

**Files:**
- Modify: `src/quiverlab/modules/homdims.py` (`LITCertificate` dataclass,
  `lit_finitistic_certificate`, and the `finitistic_dimension_bounds` wiring)
- Modify: `src/quiverlab/core/algebra.py` (`Algebra.finitistic_certificate`)
- Test: `tests/modules/test_lit_certificate.py`

**Interfaces:**
- Consumes: `is_selfinjective`, `gorenstein_dimension` (Plan 40), `phi_dim` (Task A),
  `injective_dimension` of the regular module (Plan 40 `_regular_injective_dimension`
  — refactor it to a module-private helper reused here), `global_dimension`.
- Produces:
  ```python
  @dataclass
  class LITCertificate:
      findim_upper: int | None   # certified finite findim upper bound, or honest None
      family: str | None         # "self-injective"|"gorenstein"|"finite-phidim"|"finite-id-AA"|None
      proof: str                 # the proof-carrying justification (the theorem + numbers)
      # __repr__: family -> "findim(A) <= {findim_upper}  [{family}: {proof}]";
      #           None    -> "no known decision procedure for a finite findim bound
      #                       here (A is not in a shipped LIT-decidable family)"
  def lit_finitistic_certificate(A, bound=32) -> LITCertificate
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_lit_certificate.py
"""Lat-Igusa-Todorov finitistic certificates (Plan 53 / R23c; Bravo-Lanzilotta-
Mendoza-Vivero 2002.07866). Each decidable family emits a proof-carrying finite
findim upper bound; the bound is never below the known exact findim. 'No known
decision procedure in general' -- NOT undecidable. Over QQ."""
import pytest

from quiverlab import Quiver, TruncatedPathAlgebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import (finitistic_dimension_bounds,
                                       lit_finitistic_certificate)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
def test_self_injective_findim_zero():
    A = truncated_polynomial(4, field=QQ)
    c = lit_finitistic_certificate(A)
    assert c.family == "self-injective" and c.findim_upper == 0

@lit
def test_gorenstein_nonselfinjective_findim_is_gorenstein_dim():
    # FAMILY 2 executing test (adversarial-review fix): kA3/J^2 is NOT self-injective
    # (family 1 does not preempt) and IS Iwanaga-Gorenstein because it has finite
    # global dimension (acyclic quiver => gl.dim < infinity => finite injective dim
    # both sides). For finite gl.dim, findim = gl.dim = id(A_A) = Gorenstein dim, so
    # the family-2 bound == int(A.global_dimension()). Asserting family == "gorenstein"
    # also confirms family 2 is the DECISIVE family (returned before family 3).
    A = TruncatedPathAlgebra("A3", 2, field=QQ)          # gl.dim 2, not self-injective
    c = lit_finitistic_certificate(A)
    assert c.family == "gorenstein"
    assert c.findim_upper == int(A.global_dimension())   # == 2 = Gorenstein dim

@selfcert
def test_lit_bound_respects_chain():
    # the LIT bound sits between the rigorous findim lower bound and gl.dim (family 3's
    # inequality findim <= phidim is the Task-A chain self-cert; here we tie the LIT
    # bound to [findim_lower, gl.dim]). Whichever family fires, the bound is valid.
    from quiverlab.modules.homdims import finitistic_dimension_bounds
    A = TruncatedPathAlgebra("A3", 2, field=QQ)
    c = lit_finitistic_certificate(A)
    lo = finitistic_dimension_bounds(A).lower
    assert c.findim_upper is not None
    assert lo <= c.findim_upper <= int(A.global_dimension())

@selfcert
def test_finitistic_upper_now_certified_where_it_was_none():
    # Plan 40 returned upper=None for self-injective (gl.dim infinite); the LIT
    # certificate now supplies the certified upper = 0 (family 1).
    A = truncated_polynomial(4, field=QQ)
    fb = finitistic_dimension_bounds(A)
    assert fb.lower == 0
    assert fb.upper == 0 and fb.note.startswith("LIT")       # was None pre-Plan-53

# NOTE (family 4 = finite one-sided id(A_A)): DEMOTED -- no end-to-end executing test.
# It fires only when families 1-3 fail yet id(A_A) < infinity, i.e. id(A_A) finite but
# id(_AA) infinite; the bounded engine never PROVES id = infinity, so this precondition
# is not decidably reachable. The psi_D(V)+n+1 machinery is unit-covered by a direct
# call on constructed (D, n) data (test_lit_family4_machinery), NOT through
# lit_finitistic_certificate.
```

- [ ] **Step 2: Run to verify failure**
- [ ] **Step 3: Implement** `lit_finitistic_certificate` (families 1→2→3→4, first
  applicable wins; families 1–3 have executing tests, family 4 is machinery + `# PIN`
  per the note above). **`# PIN`:** transcribe the exact `ψ_𝒟(V) + n + 1` statement
  from 2002.07866 (definition of `ψ_𝒟`, the `𝒟`, and `V`) and ARBITRATE the `+ n + 1`
  constant against families 1–2 exact values (a bound below the known `findim` is a
  bug — raise, do not clamp); HONEST-DEGRADE family 4 to `None` if the statement cannot
  be certified in-session (Plan 40 rule). Unit-cover the machinery by a direct call on
  constructed `(𝒟, n)` data (`test_lit_family4_machinery`). Then **wire
  `finitistic_dimension_bounds`**: after the `gl.dim` exact-finite branch (unchanged),
  before the `_igusa_todorov_finitistic_upper` degrade, consult
  `lit_finitistic_certificate(A)`; if it returns a finite `findim_upper`, set
  `upper = findim_upper`, `note = "LIT [{family}]: {proof}"`, keep the mandatory
  `upper >= lower` sanity gate (raise on violation). This is a BEHAVIOUR CHANGE to
  `finitistic_dimension_bounds`'s `upper`/`note` for LIT-family inputs — Plan 40's
  `test_self_injective_lower_zero_upper_valid` asserted `upper is None OR upper >=
  lower`, which STILL holds (0 ≥ 0), so no Plan-40 test breaks. **Which goldens change
  (W2, adversarial review):** the flip `None → certified` affects only inputs whose
  Plan-40 finitistic `upper` was `None` — i.e. `gl.dim = ∞` LIT-family inputs
  (self-injective / Gorenstein-with-infinite-gl.dim). The kA₂ `homological_profile`
  golden is HEREDITARY (`gl.dim 1`, finitistic already exact) so its `finitistic` entry
  does NOT change from Task C — kA₂ changes SOLELY from Task E's four new additive
  keys. Task C's cross-runner coverage of the `None → 0` flip is the NEW self-injective
  golden (`truncated_polynomial(3)`) added in Task E.
- [ ] **Step 4: Run tests** — PASS; run `tests/modules/test_finitistic.py` (Plan 40)
  to confirm no regression.
- [ ] **Step 5: Commit** — `feat(modules): LIT finitistic certificate (self-injective/Gorenstein/finite-phidim/finite-id families) -- certified findim upper flips Plan-40 None`

---

### Task D: fractional Calabi–Yau dimension (stable category, self-injective)

**Files:**
- Create: `src/quiverlab/modules/fractional_cy.py`
- Modify: `src/quiverlab/core/algebra.py`
  (`Algebra.fractional_calabi_yau_dimension`, `Algebra.is_fractionally_calabi_yau`)
- Test: `tests/modules/test_fractional_cy.py`

**Interfaces:**
- Consumes: `is_selfinjective`; `Module.nakayama()`/`.nakayama_minus()`
  (`modules/ar.py::nakayama_functor`/`_minus`); `Module.syzygy()`/`.cosyzygy()`
  (`modules/resolution.py`); `modules/decompose.py::decompose`;
  `modules/hom.py::is_isomorphic`; `modules/morphism.py::direct_sum`.
- Produces:
  ```python
  @dataclass
  class FractionalCY:
      m: int | None
      ell: int | None            # the denominator (least positive with S^ell a shift)
      cy_dimension: str | None   # "m/ell" reduced, or None
      weakly_n_cy: int | None    # the integer n = m when ell == 1 (Ivanov-Volkov), else None
      sigma_period: int | None   # order of Sigma on the generators (H2: m is a residue
                                 #  mod this; None when the search did not determine it)
      status: str                # "certified" | "budget" | "shift-trivial"
      tier: str                  # ALWAYS "weak-on-generators" -- the tier-3 certificate
                                 #  (object-wise on the simples + orbit reps; a NECESSARY
                                 #  condition for weak/strong CY, never a functor iso)
      checked_on: str            # "simples + Omega/nu-orbit reps (object-wise on a
                                 #  generating set; anchored by literature + one (m,ell)
                                 #  for all generators; can under-report m AND ell)"
      certificate: list          # per-generator [(name, "Omega^{m+ell} nu^ell ~ id")]
      # __repr__: certified -> "stable CY dimension {cy_dimension} (weak-on-generators;
      #                         weakly {weakly_n_cy}-CY)" [drop the last clause if ell!=1];
      #           budget    -> "not certified within the (m,ell) search window
      #                         (self-injective but the generators were not S-periodic in
      #                          budget -- likely representation-infinite)";
      #           shift-trivial -> "the suspension Omega^{-1} is trivial (A is radical-
      #                             square-zero self-injective local); CY dimension
      #                             degenerate, reported (0, 1)"
  def fractional_calabi_yau(A, *, ell_max=4, m_window=8, dim_budget=4096) -> FractionalCY
  def is_fractionally_calabi_yau(A, **kw) -> bool     # status == "certified"
  ```

**Algorithm (the pinned criterion):**
1. `is_selfinjective(A)` else RAISE `QuiverlabError` ("fractional CY of the stable
   category is defined for self-injective A only; ν is not an autoequivalence and
   `mod̲ A` is not triangulated otherwise", hint pointing at derived-category CY for
   finite gl.dim — the named successor).
2. Generators `G = [A.simple(v) for v in A.quiver.vertices]` (non-projective ones;
   drop any projective simple). Stable normal form `nf(X)` = `direct_sum` of the
   non-projective indecomposable summands of `decompose(X)` (loud on char caveat).
3. `Sfun(X) = nf(X.nakayama().syzygy())` (`S = Ω ν`). Detect **shift-trivial**: if
   `nf(X.cosyzygy()) ≅ nf(X)` for every generator (`Σ ≅ id`), return
   `FractionalCY(m=0, ell=1, cy_dimension="0/1", weakly_n_cy=0, sigma_period=1,
   status="shift-trivial", tier="weak-on-generators", …)`.
4. For `ell = 1 … ell_max`: compute `T_v = Sfun^ell(G_v)` (apply `Sfun` `ell` times).
   For each `m = 0, 1, …, m_window`: the pair `(m, ell)` certifies iff for EVERY `v`
   `nf(Ω^{m}(T_v)) ≅ nf(G_v)`. This is the DERIVED identity `S^ell(G_v) ≅ Σ^{+m}(G_v)`
   — apply `Ω^m` to `S^ℓ G_v ≅ Σ^m G_v = Ω^{-m} G_v` to get `Ω^m(T_v) ≅ G_v`; the sign
   is `Σ^{+m}` (matching the Math-foundation line), NOT free. Only `m ≥ 0` is searched:
   `Σ` has finite period `p` on the generators (periodic stable category), so every
   residue is hit by `m ∈ {0, …, p−1}` via `syzygy^m` (`m_window ≥ p`); `cosyzygy` is
   not needed for the canonical `m ≥ 0`. The `ell=1 ≡ Ivanov–Volkov` cross-engine test
   (`Ω^{n+1} ≅ ν^{-1}`) is a CONFIRMATION of this derived sign, not a resolution of an
   open choice (contrast P20/P39, where the sign genuinely had to be arbitrated —
   here it is derived). Return the FIRST `(m, ell)` (smallest `ell`, then smallest
   `m ≥ 0` — the canonical representative for the mod-`p` ambiguity), reducing `m/ell`,
   setting `weakly_n_cy = m` iff `ell == 1`, and recording `sigma_period = p` when the
   search determined it.
5. Budget guard: if any `nf`/`syzygy`/`nakayama` iterate exceeds `dim_budget`, or no
   `(m, ell)` is found within `ell_max`/`m_window`, return `status="budget"` (honest —
   likely representation-infinite self-injective, whose modules need not be
   `S`-periodic). `is_isomorphic`'s loud refusal propagates unchanged (never a silent
   `budget` where an iso was undecidable).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_fractional_cy.py
"""Fractional Calabi-Yau dimension of self-injective algebras in the stable category
(Plan 53 / R24). Serre functor S = Omega.nu, suspension Sigma = Omega^{-1}
(Ivanov-Volkov 1212.2619 sec 1.3; Erdmann-Skowronski). Criterion: least ell with
S^ell ~ Sigma^m. Over QQ (decompose char caveat)."""
import pytest

from quiverlab import PreprojectiveAlgebra, Quiver, linear_path_algebra, \
    truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.fractional_cy import (fractional_calabi_yau,
                                             is_fractionally_calabi_yau)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
def test_truncated_polynomial_is_1_CY():
    # k[x]/(x^a), a>=3: symmetric (nu = id), Omega^2 ~ id => weakly 1-CY, (m,ell)=(1,1).
    A = truncated_polynomial(3, field=QQ)
    fcy = fractional_calabi_yau(A)
    assert fcy.status == "certified"
    assert (fcy.m, fcy.ell) == (1, 1) and fcy.cy_dimension == "1/1"
    assert fcy.weakly_n_cy == 1
    A5 = truncated_polynomial(5, field=QQ)
    assert (fractional_calabi_yau(A5).m, fractional_calabi_yau(A5).ell) == (1, 1)

@lit
def test_dual_numbers_shift_trivial():
    A = truncated_polynomial(2, field=QQ)                # k[x]/(x^2): Omega ~ id
    fcy = fractional_calabi_yau(A)
    assert fcy.status == "shift-trivial" and (fcy.m, fcy.ell) == (0, 1)

@lit
def test_preprojective_Dynkin_is_2_CY():
    # stab Pi(Delta) is 2-Calabi-Yau (Geiss-Leclerc-Schroer): (m,ell)=(2,1).
    A = PreprojectiveAlgebra("A3", field=QQ)             # rep-finite (n<=4) => periodic
    fcy = fractional_calabi_yau(A)
    assert fcy.status == "certified"
    assert (fcy.m, fcy.ell) == (2, 1) and fcy.weakly_n_cy == 2

@xeng
def test_ell1_equals_ivanov_volkov_integer_form():
    # ell=1 certificate S ~ Sigma^n  <=>  Omega^{n+1} ~ nu^{-1} on the simples.
    A = truncated_polynomial(3, field=QQ)
    fcy = fractional_calabi_yau(A)
    assert fcy.ell == 1
    n = fcy.m
    for v in A.quiver.vertices:
        S = A.simple(v)
        lhs = S
        for _ in range(n + 1):
            lhs = lhs.syzygy()                           # Omega^{n+1} S
        from quiverlab.modules.hom import is_isomorphic
        # Omega^{n+1} S ~ nu^{-1} S  (stably; both non-projective here)
        assert is_isomorphic(_nf(lhs), _nf(S.nakayama_minus()))

@selfcert
def test_not_self_injective_refused():
    A = linear_path_algebra(3, field=QQ)                 # hereditary: not self-injective
    with pytest.raises(QuiverlabError, match="self-injective"):
        fractional_calabi_yau(A)
```

(`_nf` is a one-line test helper = the module's stable normal form; import the private
`_nf`/`_stable_nf` from `fractional_cy`. If `PreprojectiveAlgebra("A3")` proves
non-periodic in budget at implementation, switch the pin to the smallest rep-finite
`Π(A_n)`/`Π(D_4)` that IS `S`-periodic and record which in the test docstring — the
2-CY value is type-independent for Dynkin.)

- [ ] **Step 2: Run to verify failure**
- [ ] **Step 3: Implement `src/quiverlab/modules/fractional_cy.py`** per the algorithm
  above. **Adjust to reality:**
  - `nakayama_functor` self-certifies `ker g ≅ τ M` internally and returns a genuine
    `Module`; on a self-injective algebra `ν` and `ν^{-1}` are autoequivalences, so
    every iterate is a genuine module — no special-casing.
  - the SIGN of `m` in step 4 is **DERIVED**, not arbitrated: `S^ℓ ≅ Σ^{+m}` ⟺
    `Ω^m(T_v) ≅ G_v`, searched with `m ≥ 0` via `syzygy^m` only (the canonical
    representative; `Σ`-periodicity makes negative `m` unnecessary — see step 4).
    `test_ell1_equals_ivanov_volkov_integer_form` (`Ω^{n+1} ≅ ν^{-1}`) is a
    CONFIRMATION of the derived sign, NOT the arbiter of an open choice (unlike P20/P39,
    whose signs genuinely had to be arbitrated).
  - the two commuting facts (`ν Ω ≅ Ω ν`) are NOT relied on in code (compute `S = Ων`
    literally); they are only used in the docstring derivation.
  - **Honest-scope docstring (three-tier language):** state plainly that this certifies
    the fractional CY dimension at the **weak-on-generators tier** — object-wise on the
    simples + orbit reps, a NECESSARY condition for weak `(m,ℓ)`-CY (and hence strong),
    never a functor isomorphism — anchored by the pinned values + one `(m, ℓ)` for all
    generators, and that it can under-report BOTH `m` and `ℓ`. The tier-1 functorial
    certificate is the successor `Ω^{n+1}_{Aᵉ}(A) ≅ (A^∨)_φ` (Ivanov–Volkov Thm 1.8, via
    `engine/resolutions_minimal.py`).
  - wire the two `Algebra` delegates.
- [ ] **Step 4: Run tests** — PASS.
- [ ] **Step 5: Commit** — `feat(modules): stable-category fractional Calabi-Yau dimension of self-injective algebras (S=Omega.nu, Sigma=Omega^{-1}; Ivanov-Volkov criterion; bounded search + loud budget)`

---

### Task E: GUI — extend `homological_profile` + new `fractional_cy` kind

**Files:**
- Modify: `src/quiverlab/modules/homdims.py` (`homological_profile` block builder:
  add `phidim`/`psidim`/`spectrum`/`lit` entries) + a new `fractional_cy_block(A)`
  shared builder in `fractional_cy.py`
- Modify: `src/quiverlab/hpc/spec.py` (extend the `homological_profile` branch at
  `spec.py:1462`; add a `fractional_cy` branch; `_snip` recipe ~`spec.py:2292`)
- Modify: `docs/gui/runner.py` (twin: extend the `homological_profile` elif at
  `runner.py:863`; add `fractional_cy`; ETA `"scalars"` entry ~`runner.py:1259`)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox + `S.ids` +
  push-list + block renderer for the new rows and the new kind)
- Modify: `webapp/static/app.js` (webapp block renderer)
- Modify: `webapp/server/i18n/en.json` + `es.json`
- Modify: `src/quiverlab/trace/results_html.py` (report renderer branches)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (three golden changes, all documented in the change-log docstring):
  (1) re-freeze `homological_profile_kA2` — kA₂ is HEREDITARY (`gl.dim 1`), so its
  `finitistic` entry is UNCHANGED by Task C; it changes SOLELY from Task E's four new
  additive keys (`phidim`/`psidim`/`phi_spectrum`/`lit`);
  (2) ADD a NEW self-injective `homological_profile_kxx3` golden
  (`truncated_polynomial(3)`) — this is where Task C's `None → 0` `finitistic` flip is
  covered CROSS-RUNNER (its `finitistic.upper` is `0`, `note` starts `"LIT"`, and its
  `lit.family == "self-injective"`);
  (3) ADD ONE new `fractional_cy_kxx3` golden.
- Test: `tests/webapp/test_phidim_fcy_gui_p53.py`

**Interfaces & produces:**
- `homological_profile(A)` gains (schema-stable — same block name, additive keys):
  ```python
  "phidim":   {"value": int, "exact": bool, "status": str, "text": str},
  "psidim":   {"value": int, "exact": bool, "status": str, "text": str},
  "phi_spectrum": {"values": [int], "gaps": [int], "complete": bool} | {"error": ...},
  "lit": {"findim_upper": int | None, "family": str | None, "proof": str},
  ```
  each honest-noted; the φdim/spectrum entries catch the decompose char-caveat refusal
  into `{"error": ...}` (the Plan-40 IT-entry precedent), never a silent omission.
  `"references"` gains `fernandes_lanzilotta_mendoza`, `barrios_mata`,
  `bravo_lanzilotta_mendoza_vivero`.
- **Additive-key + renderer-tolerance contract (W3, adversarial review — the
  Plan-25 result cache).** The new `homological_profile` keys are strictly ADDITIVE
  (the block name is unchanged), so canonical keys are unaffected and a bump of the
  library `__version__` is NOT required for this change. Pre-P53 cached blocks replay
  with the OLD shape (missing the four new keys). The chosen route is **renderer
  tolerance** (the established additive-golden pattern): both GUI renderers and the
  `results_html.py` report branch MUST guard every new key
  (`block.get("phidim")`, `block.phi_spectrum ?? null`, etc.) and render nothing for a
  missing key — a replayed pre-P53 block shows the old four dimensions, never a
  crash. State this explicitly in both renderers' comments. (The `finitistic` value
  itself changing under Task C for self-injective inputs is fine for the cache: a
  cached pre-P53 self-injective block simply carries the old `upper=None` text and
  replays as-is; new computes get the certified `0`.)
- New `fractional_cy` block (algebra-level scalar kind — routes through `_dispatch`,
  NOT `_dispatch_module`; schema stays v1, no request block):
  ```python
  {"kind": "fractional_cy",
   "cy_dimension": str | None, "m": int | None, "ell": int | None,
   "weakly_n_cy": int | None, "sigma_period": int | None,
   "status": str, "tier": "weak-on-generators", "checked_on": str,
   "certificate": [...], "text": str,
   "references": ["ivanov_volkov", "erdmann_skowronski", "assem_book"]}
   | {"kind": "fractional_cy", "error": "<the loud not-self-injective message>"}
  ```
  A non-self-injective input is caught into `{"error": ...}` (a clean typed 4xx in the
  webapp, never a 500 — the Plan-26 relation-violation precedent). The `tier` field is
  ALWAYS `"weak-on-generators"` and the renderer prints it beside the value so the
  displayed result never overstates the certificate.

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir;
  copy the runner-pair fixture from `tests/webapp/test_module_blocks_m0729.py`):
  assert the extended `homological_profile` block shape (phidim/psidim/spectrum/lit
  present with values for kA3/rad² over GF(32003)); assert the `fractional_cy` block
  for `k[x]/(x^3)` (`cy_dimension == "1/1"`); assert BOTH runners byte-identical
  (`json.dumps(sort_keys=True)`) on both blocks.
- [ ] **Step 2: Implement the builders + the seven-touchpoint GUI wiring** (P38 Task 6
  checklist): checkbox id `qlgui-fractional_cy` + `S.ids` + push-list in both
  `gui.js`; block renderers in both `gui.js` + `app.js` — the `homological_profile`
  renderer adds the φdim/ψdim rows, the spectrum as a value list with gaps
  highlighted, and the LIT proof line, EACH BEHIND A MISSING-KEY GUARD (W3:
  `block.get(...)` / `?? null`, render nothing when absent, so a replayed pre-P53
  cached block shows the old four dimensions and never crashes); the CY block renders
  `cy_dimension`, `weakly_n_cy`, `sigma_period`, the `tier` = `weak-on-generators`
  label beside the value (so the display never overstates the certificate), the
  `checked_on` honesty line, and the per-generator certificate. ETA scalars
  (`"fractional_cy": 3.0` — it iterates ν/Ω); i18n keys (`inv.fractional_cy`,
  `block.fractional_cy.title` + per-row labels incl. `tier`/`sigma_period`, and the new
  `homological_profile` sub-labels `phidim`/`psidim`/`spectrum`/`gaps`/`lit`, EN + ES);
  `_snip` recipe (`"fractional_cy": lambda it: "A.fractional_calabi_yau_dimension()"`);
  `results_html.py` report branches (the CY block states `S = Ω∘ν`, `Σ = Ω⁻¹`, the
  criterion, and the honest weak-on-generators tier-3 scope) — also guarded for the
  new `homological_profile` keys.
- [ ] **Step 3: Re-freeze goldens (W2 narrative).** Run `test_runner_delegation.py`
  FIRST to see the diff: `homological_profile_kA2` changes ONLY by the four new
  additive keys (kA₂ is hereditary — its `finitistic` entry is untouched by Task C; do
  NOT add a "confirm only finitistic changed" check, that claim was wrong). Merge
  `_runner_goldens.json` SEMANTICALLY (dict-union of both parents, never textually),
  re-freeze `homological_profile_kA2`, and ADD the two new goldens
  `homological_profile_kxx3` (self-injective — the cross-runner home of Task C's
  `None → 0` flip) and `fractional_cy_kxx3`. Document ALL THREE in the delegation
  test's change-log docstring.
- [ ] **Step 4: Run the gates** —
  `... -m pytest tests/webapp/test_phidim_fcy_gui_p53.py tests/webapp/test_runner_delegation.py tests/hpc -q`
  then `.venv/bin/node --check` on both `gui.js` + `app.js` and
  `tests/webapp/test_js_parses.py`. PASS.
- [ ] **Step 5: Commit** — `feat(gui,webapp,hpc): homological_profile phidim/psidim/spectrum/LIT rows + fractional_cy compute kind (both runners byte-identical, EN/ES)`

---

### Task F: QPA crosschecks + honest scope

**Files:**
- Create: `tests/qpa/test_phidim_fcy_qpa.py`
- Modify (if a live verb is found): `src/quiverlab/qpa/scripts.py`,
  `src/quiverlab/qpa/crosscheck.py`

**Interfaces:** probe live via the `NamesGVars()` precedent
(`tests/qpa/test_products_qpa.py`). **QPA has NO φdim/ψdim, NO Igusa–Todorov
φ-dimension, and NO stable-CY-dimension surface** (verified expectation — the probe
SKIPS honestly for these and FAILS only if QPA ever ships them). What QPA CAN
crosscheck, and this battery does, are the PREREQUISITES the values rest on:
- `IsSelfinjectiveAlgebra` ↔ `A.is_selfinjective()` (the fractional-CY scope gate) on
  `k[x]/(x^a)`, `Π(A_n)`, and a non-self-injective control;
- `NakayamaAutomorphism`/`IsSymmetricAlgebra` ↔ `A.nakayama_automorphism()` /
  `A.is_symmetric()` (ν = id ⇔ symmetric, the `k[x]/(x^a)` fractional-CY derivation);
- the syzygy period underpinning the CY value: QPA `NthSyzygy`/`1st_syzygy` of the
  simples on `k[x]/(x^a)` ↔ `Module.syzygy()` orbit (the `Ω² ≅ id` fact), where a QPA
  syzygy verb exists.

- [ ] **Step 1: Probe + write the battery** (concrete tests in the
  `tests/qpa/test_tor_qpa.py` style; header
  `pytestmark = pytest.mark.skipif(session.should_skip_qpa(), ...)`). The φdim/CY
  probe SKIPS honestly (no QPA surface); the self-injective/symmetric/syzygy
  crosschecks run live.
- [ ] **Step 2: Run** `... -m pytest tests/qpa/test_phidim_fcy_qpa.py -q -m qpa`
  (venv has `[qpa]`). PASS; the φdim/CY probe SKIPS, the prerequisite crosschecks
  agree.
- [ ] **Step 3: Commit** — `test(qpa): self-injective/symmetric/syzygy prerequisites for phidim + fractional-CY (QPA has no phidim/CY surface -- honest skip)`

---

### Task G: citations, verification page, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/registry.py` + `references.bib`
- Modify: `docs/verification.md`, `README.md`
- Test: existing release gates

- [ ] **Step 1: Citations** (verified BibTeX; `_r(...)` registry precedent; the
  `bibtex()` helper hard-fails if `registry.py`/`references.bib` disagree). Add:
  - `fernandes_lanzilotta_mendoza` — Fernandes, Lanzilotta, Mendoza, "The Φ-dimension:
    A new homological measure," Algebr. Represent. Theory 18 (2015) 463–476
    (arXiv:1304.0754). [φdim definition + derived-invariance.]
  - `bravo_lanzilotta_mendoza_vivero` — Bravo, Lanzilotta, Mendoza, Vivero,
    "Generalised Igusa-Todorov functions and Lat-Igusa-Todorov algebras," J. Algebra
    (2021) (arXiv:2002.07866). [LIT definition + the `ψ_𝒟(V)+n+1` findim bound.]
  - `ivanov_volkov` — Ivanov, Volkov, "Stable Calabi–Yau dimension of self-injective
    algebras of finite type" (arXiv:1212.2619). [`S = Ω ν`, `Σ = Ω⁻¹`, the
    `Ω^{n+1} ≅ ν⁻¹` criterion + Table 1.]
  - `erdmann_skowronski_scy` — Erdmann, Skowroński, "The stable Calabi–Yau dimension
    of tame symmetric algebras," J. Math. Soc. Japan 58 (2006). [Introduced the stable
    CY dimension of a self-injective algebra.]
  - `geiss_leclerc_schroer` — Geiß, Leclerc, Schröer, "Rigid modules over preprojective
    algebras," Invent. Math. 165 (2006). [`stab Π(Δ)` is 2-CY.]
  Reuse `igusa_todorov`, `barrios_mata`, `preprojective`, `assem_book` (present).
- [ ] **Step 2: Verification page** — add the Plan-53 subsystem rows:
  - `modules/homdims.py` (φdim/ψdim/spectrum/LIT): `oracle_selfcert` — the standing
    chain `findim ≤ φdim ≤ ψdim ≤ gldim` + the `gldim`-finite collapse + BMR
    `{1, m−1} ⊆ Spec_φ` + LIT bound ≥ known-findim; `oracle_literature` — rep-finite
    φdim closed values (kA_n, Nakayama), self-injective φdim = 0; `oracle_crossengine`
    — ⊕-of-all ≡ max-over-indecomposables; `qpa` — self-injective/symmetric
    prerequisites.
  - `modules/fractional_cy.py`: `oracle_literature` — `k[x]/(x^a)` (1,1), `Π(Δ)` (2,1),
    dual-numbers shift-trivial; `oracle_selfcert` — `Ω^{m+ℓ} ν^ℓ ≅ id` identity;
    `oracle_crossengine` — `ℓ=1 ≡` Ivanov–Volkov `Ω^{n+1} ≅ ν⁻¹`.
  - **Honest-scope entries:** (a) φdim/ψdim are EXACT only for representation-finite
    input; otherwise a certified LOWER bound (never a claimed sup) — the ⊕-of-all
    theorem needs the closed AR knit; a `status="error"` knit RAISES (only
    `status="budget"` degrades softly). (b) LIT: **no known decision procedure in
    general** (NOT undecidable); families 1–3 emit certified executing bounds, family 4
    (one-sided finite `id`) ships machinery + `# PIN` but has NO end-to-end executing
    test (its precondition is not certifiable by the bounded engine, which never proves
    an injective dimension infinite). (c) Fractional CY is certified at the
    **weak-on-generators tier** (tier 3 of three: strong = functor iso, weak =
    object-wise on all objects, weak-on-generators = object-wise on the simples + orbit
    reps). This is a NECESSARY condition for the weak (hence strong) CY property; it can
    under-report BOTH the numerator `m` (defined mod the Σ-period) AND the denominator
    `ℓ` (a spurious small `ℓ` may pass on the generators while failing on a
    non-generator object). The tier-1 functorial certificate `Ω^{n+1}_{Aᵉ}(A) ≅
    (A^∨)_φ` (Ivanov–Volkov Thm 1.8) is the named successor. Complete only for
    finite-type self-injective (periodicity); rep-infinite caps at budget. (d)
    Ivanov–Volkov Table 1 rows are BLOCKED-until-transcribed (Asashiba type →
    presentation mapping deferred). (e) QPA has NO φdim/ψdim/stable-CY surface — the
    literature battery is their oracle.
  Recount the class table (`tests/release/test_oracle_classes.py` drives the numbers —
  run collection, paste, re-run to green).
- [ ] **Step 3: README** — one features line: "φdim / ψdim as algebra invariants, the
  φ-spectrum and its gaps, Lat-Igusa-Todorov finitistic certificates, and the stable
  fractional Calabi–Yau dimension of self-injective algebras — clickable via
  `homological_profile` + `fractional_cy`."
- [ ] **Step 4: Full gate** —
  `... -m pytest tests/modules -q` (deep, touched files),
  `... -m pytest -q -m fast`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release -q`,
  `... -m pytest tests/citations/test_bib_structure.py tests/webapp/test_js_parses.py tests/webapp/test_runner_delegation.py -q` — all green.
- [ ] **Step 5: Commit** — `docs(verification): Plan-53 oracle rows + honest scope (phidim rep-finite, LIT no-general-procedure, fractional-CY object-wise + IV Table-1 deferred) + recounted classes; citations`

---

## Acceptance (Plan-53 definition of done)

1. `phi_dim`/`psi_dim` public on `Algebra`/library, EXACT for representation-finite
   input (⊕-of-all-indecomposables via add-monotonicity), a labelled certified LOWER
   bound otherwise (never a claimed sup), self-injective short-circuiting to the exact
   `0`; the standing chain `findim ≤ φdim ≤ ψdim ≤ gldim` self-certified on the
   computed terms (all equal when gldim is exact-finite). Decompose char caveat
   inherited and tested (loud over `char ≤ dim`).
2. `phi_spectrum` public: the rep-finite φ-value set + gaps, with the Barrios–Mata–Rama
   `{1, φdim−1} ⊆ Spec_φ` fact pinned; a partial spectrum claims no gaps.
3. `lit_finitistic_certificate` public: **families 1–3 ship certified EXECUTING
   bounds** (self-injective → 0, tested on `k[x]/(x⁴)`; Iwanaga–Gorenstein →
   Gorenstein dim, tested on the non-self-injective `kA₃/J²`; finite φdim → φdim,
   covered by the Task-A chain self-cert + the LIT bound-consistency test), each never
   below a known exact findim; **family 4 (one-sided finite `id`) ships the
   `ψ_𝒟(V)+n+1` machinery + `# PIN`, unit-covered on constructed `(𝒟,n)` data but with
   NO end-to-end executing test** (its "finite one-sided but not Gorenstein"
   precondition is not certifiable by the bounded engine, which never proves an
   injective dimension infinite). Wired into `finitistic_dimension_bounds` so Plan 40's
   honest `None` upper becomes a certified number exactly in the applicable families;
   "no known decision procedure in general" (NOT undecidable) elsewhere.
4. `fractional_calabi_yau` / `is_fractionally_calabi_yau` public: the stable-category
   `(m, ℓ)` of a self-injective algebra via `S = Ω ν`, `Σ = Ω⁻¹`, certified at the
   **weak-on-generators tier** (object-wise on the simples + orbit reps — a necessary
   condition, labelled as such in the payload `tier` field, honestly noted to
   under-report both `m` and `ℓ`) with a bounded `(m, ℓ)` search + `sigma_period`
   reported + loud budget refusal + loud non-self-injective refusal; `k[x]/(x^a)` =
   (1,1), dual numbers = shift-trivial (0,1), `Π(Δ)` = (2,1) pinned, and `ℓ=1 ≡`
   Ivanov–Volkov `Ω^{n+1} ≅ ν⁻¹` cross-checked (the derived sign confirmed).
5. GUI: `homological_profile` gains φdim/ψdim/spectrum/LIT rows (strictly additive,
   renderers guard missing keys so pre-P53 cached blocks replay) and a new
   `fractional_cy` compute kind is clickable end-to-end (canvas → block → report) in
   EN+ES, both runners byte-identical, schema stable (v1), the THREE touched/added
   goldens documented (`homological_profile_kA2` re-frozen for the new keys, new
   self-injective `homological_profile_kxx3`, new `fractional_cy_kxx3`); a
   non-self-injective `fractional_cy` request returns a clean typed error entry, never
   a 500.
6. QPA battery green live (`-m qpa`): self-injective/symmetric/syzygy prerequisites
   agree; the φdim/CY probe SKIPS honestly and FAILS if QPA ever ships such a surface.
7. `docs/verification.md` recounted with the five honest-scope entries (rep-finite
   φdim, LIT no-general-procedure, fractional-CY object-wise-on-generators, IV Table-1
   deferred, QPA no-surface); README line added; citations added and bib-synced; deep
   (touched dirs) + fast + qpa + release + js-parse + bib-structure suites green.

---

## Methodology & assumptions

**Approach.** I read Plan 40 end-to-end (the plan this extends) and the shipped
`src/quiverlab/modules/homdims.py` to establish EXACTLY what already exists (φ/ψ with
last-strict-drop + Fitting-closure termination, the self-injective φ≡0 short-circuit,
`finitistic_dimension_bounds`'s honest `None` upper, `DominantDimension`/
`GorensteinDimension` honesty dataclasses); I skimmed Plan 49 for the AR-knit budget
pattern (`knit_ar_quiver` `.is_complete`/`.status`, and its REFUSAL of self-injective
input) which the φdim rep-finite route reuses. I grepped the codebase and confirmed the
fractional-CY prerequisites already ship at the MODULE level: `nakayama_functor`/
`nakayama_functor_minus` (ν, ν⁻¹) and `stable_hom_dim` in `modules/ar.py`,
`syzygy`/`cosyzygy` (Ω, Ω⁻¹) in `modules/resolution.py`, `is_selfinjective` in
`modules/ext.py`, `decompose`/`is_isomorphic` in `modules/{decompose,hom}.py`, and the
`Module.nakayama()/.syzygy()/.cosyzygy()` methods — so the fractional-CY certificate
needs NO new engine. I re-verified all six citation clusters via WebSearch/WebFetch and
transcribed the authoritative Serre-functor criterion and table page-image from
Ivanov–Volkov 1212.2619.

**Assumptions (each load-bearing).**
1. **φdim = φ(⊕ all indec) for rep-finite, and both φ, ψ are add-monotone** — taken
   from the survey 2310.09283 (WebFetch, explicit) and 1304.0754. If the ⊕-of-all
   route ever disagreed with `max φ over indecomposables`, Task A's `oracle_crossengine`
   test catches it (I designed that test precisely as the guard).
2. **The chain is `findim ≤ φdim ≤ ψdim ≤ gldim`** with the `≤ gldim` end — confirmed
   verbatim by the survey WebFetch. The self-cert checks only computed/exact terms
   because `findim` and `gldim` are generally lower bounds only.
3. **Serre functor `S = Ω ν`, suspension `Σ = Ω⁻¹`, criterion `Ω^{n+1} ≅ ν⁻¹`** for
   self-injective `A` — transcribed from Ivanov–Volkov 1212.2619 §1.3 (I read the page
   images: eq. 0.1, the `mod̲ A weakly n-CY ⟺ Ω^{n+1} ≅ ν⁻¹` line, and `ν ≅ id ⟺
   symmetric`). The `ℓ ≥ 2` fractional generalisation `S^ℓ ≅ Σ^m` is the standard
   fractional-CY definition (also confirmed by the WebSearch result quoting the
   generic `ν^ℓ ≅ [m]` definition). The `k[x]/(x^a) = (1,1)` value is my derivation
   FROM this criterion (symmetric ⇒ ν=id; `Ω(k[x]/(x^i)) = k[x]/(x^{a−i})` ⇒ Ω²≅id;
   least `n` with Ω^{n+1}≅id is 1), and it is in-engine verifiable, so the plan's test
   is self-checking, not memory-dependent.
4. **`stab Π(Δ)` is 2-CY** — literature (GLS; the WebSearch results state it directly).
   The value is Dynkin-type-independent; the plan pins `Π("A3")`/`Π("A4")` (rep-finite,
   hence periodic) with an in-engine cross-check and an explicit implementation-time
   fallback if the chosen instance is non-periodic. I did NOT hand-run the engine on
   `Π(A3)` to confirm the exact numbers (I am the plan writer, not the implementer);
   the test is written to catch a wrong pin.

**What I deliberately did NOT do, and why.**
- I did NOT run any code, modify `src/` or `tests/`, or commit — the deliverable is the
  plan document only (task constraint).
- I did NOT transcribe Ivanov–Volkov Table 1 into concrete oracle values: the table is
  keyed by Asashiba type `(Δ, f, t)`, and mapping a type triple to a quiverlab
  presentation is itself non-trivial. I marked these `BLOCKED-until-transcribed` and
  rested fractional-CY acceptance on the two independently-solid pins (`k[x]/(x^a)`
  self-derived; `Π(Δ)` bedrock literature) + the internal-consistency + Ivanov–Volkov
  integer-form cross-engine oracle. This is the honest move the task explicitly
  requested.
- I did NOT design a FUNCTORIAL fractional-CY certificate. Object-wise iso on a
  generating set is a necessary condition, not sufficient for a functor isomorphism; I
  named the rigorous upgrade (Ivanov–Volkov Thm 1.8, `Ω^{n+1}_{Aᵉ}(A) ≅ (A^∨)_φ`, via
  the existing minimal `A^e` engine) as a successor and put the limitation on the
  verification page — consistent with how the shipped periodicity/AR/almost-split
  engines already certify (object iso + literature anchor).
- I could not extract the LIT `ψ_𝒟(V)+n+1` internals from the 2002.07866 PDF (compressed
  streams defeated WebFetch), so I kept the record's verbatim formula and marked the
  exact `ψ_𝒟`/`V`/constant as a `# PIN` for the implementation session, with the
  mandatory arbitration-against-known-findim safety gate — the same honest-degrade
  discipline Plan 40 uses for its own IT finitistic bound.

**Why I believe the result is correct.** The φdim/ψdim/chain/monotonicity facts are
quoted verbatim from the survey and the primary φdim paper; the fractional-CY criterion
is transcribed from the authoritative Ivanov–Volkov page images (not memory); every new
verb reuses shipped, tested primitives (φ/ψ, ν, Ω, decompose, is_isomorphic, the AR
knit) with no new math engine; every honesty boundary (rep-finite exactness, LIT
decidable-families-only, object-wise-on-generators, budget/scope refusals) is a loud
typed refusal or a labelled marker matching existing Plan-40/49 precedents; and every
pinned number is either self-derivable in-engine (so the tests are self-checking) or a
bedrock literature value with an in-engine cross-check and an implementation-time
fallback. The two places I could not fully verify at spec time (LIT internals; IV
Table 1) are explicitly flagged BLOCKED/PIN with safety gates, per the house honesty
rule.

---

## Change log

- **2026-08-07 — plan written** (spec P53, R23 + R24; all six citation clusters
  re-verified via WebSearch/WebFetch; Ivanov–Volkov 1212.2619 §1.3 criterion + table
  transcribed from page images).
- **2026-08-07 adversarial review: 2 majors + 7 minors applied** (weak-CY criterion
  relabeled object-wise, three-tier language unified; LIT families 2/4 given executing
  tests or honestly demoted; Σ-period normalization stated; knit error/budget split;
  sign gloss fixed; golden narrative corrected + self-injective cross-runner golden;
  cache-shape tolerance; decompose cost note).
