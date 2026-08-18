# Plan 71: HH• as a graded Lie module over HH¹ (P71 / R12)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking. Do the reference re-verification in the
> Record/Reference sections BEFORE writing any oracle pin — the `# PIN` fields below must be
> resolved against the sources, not this plan's prose. **This plan BUILDS ON Plan 70** (the
> merged doc `docs/plans/2026-08-08-plan-70-hh1-lie.md` on `dev` — commit `d5a26c7` at authoring,
> but `dev` advances, so refer to the DOC, not any tip hash): it consumes
> P70's `invariants/hh1_lie.py::derivations(A)` / `inner_derivations(A)` (the Der/Inn
> representation), its char-0 classification (`hh1_lie_structure`), its `DEFAULT_MAXDIM = 48`
> budget and its arithmetic-field base-change-note gating — **honor those contracts, do not
> respec them.** P70's critic explicitly verified the P71 seam: the degree-1 Lie derivative
> `L_D f = D∘f − Σ f(…,D a_i,…)` uses `D` only as a map `A → A`, exactly what `derivations(A)`
> returns — **no bracket-engine / basis reconciliation is needed on the primary route.**
>
> **Two card corrections are binding (evidence in the Reference and Live-verified sections):**
> (1) the reference the metaplan card wrote as "Assem–Lanzilotta–Schroll(?) 1803.10310" is in
> fact **Artenstein–Lanzilotta–Solotar, "Gerstenhaber structure on Hochschild cohomology of
> toupie algebras"** (arXiv abstract fetched this round — the card's guess was wrong, R12's
> attribution is right and now confirmed); (2) the metaplan sl₂ example inherited from P70 is the
> **Kronecker path algebra `kK₂`** (`HH¹ ≅ sl₂`, `HH¹ = HH^1` = the *irreducible* 3-dimensional
> adjoint module `L(2)`), NOT its trivial extension.

**Goal.** Make the whole Hochschild cohomology `HH•(A) = ⊕ₙ HH^n(A)` a computable **graded Lie
module over the Lie algebra `HH¹(A)`**, field-general on the primary route, with the char-0
weight/torus decomposition and the indecomposable-summand decomposition. The action is the
**degree-1 component of the Gerstenhaber bracket** `[·,·]: HH¹ ⊗ HH^n → HH^n`, which on cochains
is the **Lie derivative**

```
   (L_D f)(a_1, …, a_n) = D(f(a_1, …, a_n)) − Σ_{i=1}^{n} f(a_1, …, D(a_i), …, a_n),
```

for a derivation `D` representing a class in `HH¹ = Der(A)/Inn(A)` (P70) and an `n`-cochain `f`.

- **Any exact field** (`CC`/`QQ`, `GF(p)`, `GF(pⁿ)`): the action matrices `ρ_n(D) ∈ End(HH^n)`
  for `D` in a `Der/Inn` basis of `HH¹` (P70), the graded-module structure `⊕ₙ HH^n`, and the
  self-certified module axiom `ρ_n([D,E]) = [ρ_n(D), ρ_n(E)]`. **`L_D` needs only `D` as a map
  `A → A` and the normalized bar cochain complex over `A.domain`** — it does NOT need the
  window-bounded GF(p) Gerstenhaber-bracket engine (Plan 35 / P51), which enters ONLY as an
  in-window cross-engine SIGN ARBITER (Task 4).
- **Characteristic 0 (a HARD, loud gate, inherited from P70):** the **weight/torus
  decomposition** of each `HH^n` with respect to a maximal torus of the classified `HH¹`
  (P70's char-0 classification), and the **indecomposable Lie-module summand decomposition**
  `HH^n = ⊕ M_i` via the `modules/decompose.py` Fitting machinery (char-0 cleanest; the
  Dickson/Cohen–Ivanyos–Wales `char > dim` allowance carries over).

**GUI.** ONE new **algebra-only** budget-carrying compute kind `hh_lie_module` (a `top`-carrying
HH kind, the `gerstenhaber_brackets`/`hochschild` precedent), served by all three tiers
(`hpc/spec.py` + the Pyodide twin `docs/gui/runner.py` byte-identical, plus the
`webapp/server/schema.py` grammar site), the theme's picker row, all four locales
(en/es/fr/zh) with exact key parity, canonical keys stable, a **decomposition-table report
block**, and citations.

**Architecture.** ONE new source module + one HH kind; everything else is a thin exact layer
over primitives already on `dev` (Plan-35 TT-calculus siblings) plus P70's `Der/Inn` surface:

- **`src/quiverlab/hochschild/lie_module.py`** (new) — the graded-Lie-module surface, a sibling
  of `hochschild/products.py` (the Plan-35 TT-calculus cup/cap/bracket surface). It is the
  **Lie-MODULE side of the Tamarkin–Tsygan calculus** (CSSS 2311.08003 frames exactly this:
  `HH•` as a module over the graded Lie algebra `HH¹`). Public:
  - `lie_module_action(A, top, *, budget=DEFAULT_MAXDIM) -> HHLieModule` — the full report: per
    degree `n ∈ 0..top`, the action matrices `ρ_n(D)` (structure constants, basis-provenance
    tagged), the module-axiom self-cert, the char-0 weight decomposition, and the indecomposable
    summands.
  - `lie_derivative_on_hh(A, D, n) -> matrix` — the single-derivation action `ρ_n(D)` on
    `HH^n` (the field-general primitive; the whole report is built from it).
  - `hh_lie_module_block(A, top, *, budget=...) -> dict` — the JSON block for the runners.
- The `hh_lie_module` kind is wired into `hpc/spec.py::_dispatch` + `docs/gui/runner.py`
  (byte-identical twins) + `webapp/server/schema.py` (the third grammar-parse site, mirroring the
  Plan-35 `gerstenhaber_brackets`/`cup_products` `top`-suffix branch), the GUI touchpoints
  (checkbox / `S.ids` / push-list / `renderBlock` / `scheduleProbe`), the i18n chains, and
  `trace/results_html.py` (the decomposition-table render branch).
- **`modules/decompose.py`** gains a thin **matrix-level** entry
  `decompose_representation(gens, domain, *, budget) -> [Summand]` (see §4): the existing
  `decompose(M)` is quiver-coupled (`hom_space` enumerates `M.algebra.quiver`), so P71 factors
  the Fitting/min-poly/trace-form CORE out to drive it on the Lie-module representation
  `{ρ_n(D)}` directly. `decompose(M)` becomes a thin wrapper over it — byte-stable behaviour.

**Tech stack.** Exact `Domain` linear algebra via `quiverlab.fields.linalg` (`nullspace`,
`rank`, `rref`, `solve`) on the **normalized bar cochain complex** over any Domain
(`hochschild/bar.py::coboundary_matrix` / `_cochain_basis`, the field-general substrate P19/P35
already compute HH dims on) plus P70's `derivations(A)`. **No floats in `src/`** (AST-gated by
`tests/test_no_floats.py`): action structure constants are exact-`Domain` strings, dims/weights
are `int`, verdicts `bool`/`str`. No resolution engine on the primary path.

---

## Record (verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`, R12)

> **R12 — HH• as a graded Lie module over HH¹.** [A-scout P2; keep-with-corrections]
> Object: the bracket action [HH¹, HH^n] (degree-1 case = Lie derivative of a derivation —
> implementable field-generally WITHOUT the full bracket engine; state which method),
> weight/torus decomposition, indecomposable Lie-module summands (reuse modules/decompose.py;
> char-0 cleanest). Refs: Artenstein–Lanzilotta–Solotar arXiv:1803.10310;
> Meinel–Nguyen–Pauwels–Redondo–Solotar arXiv:1803.10909 (J. Algebra 580 (2021);
> Virasoro-subquotient statement verified verbatim);
> Chaparro–Schroll–Solotar–Suárez-Álvarez arXiv:2311.08003 (J. Algebra 708 (2026) 138–231 —
> journal ref independently confirmed). Size M–L. Deps: R11.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P71 (Wave 4, tier δ, **needs P70**;
builds the Lie-module action `[HH¹, HH^n]` on P70's `Der/Inn` surface). Size M–L.

---

## Reference re-verification (done at authoring; findings binding, `# PIN` = resolve at citation-add)

The metaplan standing rule (§1.5) requires the plan writer to re-verify its record's citations.
Findings (all three arXiv abstracts fetched this authoring round; the MNPRS journal vol/pages
matched a ScienceDirect search hit but is `# PIN`'d for a landing-page re-confirm — see finding 2):

1. **Artenstein, D., Lanzilotta, M. & Solotar, A. — "Gerstenhaber structure on Hochschild
   cohomology of toupie algebras."** arXiv:1803.10310 [math.RT]. **Abstract verified verbatim
   this round:** *"We study homological properties of a family of algebras called toupie
   algebras. Our main objective is to obtain the Gerstenhaber structure of their Hochschild
   cohomology, with the purpose of describing the Lie algebra structure of the first Hochschild
   cohomology space, together with **the Lie module structure of the whole Hochschild
   cohomology**."* This is precisely P71's deliverable for the *toupie* family (the source of
   the "Lie-module of `HH•`" surface). **CARD CORRECTION (binding):** the metaplan §5 P71 card
   wrote the authors as *"Assem–Lanzilotta–Schroll(?)"* — that guess is WRONG; the actual
   authors are **Artenstein–Lanzilotta–Solotar** (as R12 already recorded), now confirmed
   against the arXiv abstract. `# PIN`: journal ref (the arXiv page shows none this round; ship
   `@misc` with `note = {arXiv:1803.10310}` unless a journal ref is confirmed at citation-add —
   never fabricate a venue).

2. **Meinel, J., Nguyen, V. C., Pauwels, B., Redondo, M. J. & Solotar, A. — "The Gerstenhaber
   structure on the Hochschild cohomology of a class of special biserial algebras."**
   arXiv:1803.10909 [math.RA]; **J. Algebra 580 (2021) 264–298** (independently confirmed this
   round via ScienceDirect `S0021869321002015`). **Abstract verified verbatim; the Virasoro
   sentence quoted exactly:** *"In degree one, we show that the cohomology is isomorphic, as a
   Lie algebra, to a direct sum of copies of **a subquotient of the Virasoro algebra**. These
   copies share Virasoro degree 0 and commute otherwise. Finally, **we describe the cohomology
   in degree `n` as a module over this Lie algebra by providing its decomposition as a direct
   sum of indecomposable modules**."* This is the literal statement of P71's char-0 indecomposable-
   summand deliverable and the source of the **Virasoro-subquotient** structure (a truncated
   Witt algebra — see §5, where `HH¹(k[x]/(x^n))` is a positive-Witt truncation, an *analogue* of
   a Virasoro subquotient). `# PIN`: the **volume/pages (J. Algebra 580 (2021) 264–298)** matched a
   ScienceDirect search hit (`S0021869321002015`) this round but were NOT confirmed against the
   publisher's article landing page — re-verify at citation-add (the critic could not
   independently re-confirm; the arXiv id + authors + Virasoro abstract ARE confirmed).

3. **Chaparro, C., Schroll, S., Solotar, A. & Suárez-Álvarez, M. — "The Hochschild cohomology
   and the Tamarkin–Tsygan calculus of gentle algebras."** arXiv:2311.08003 [math.RT]; J.
   Algebra 708 (2026) 138–231 (per R12, independently confirmed). **Abstract verified verbatim
   this round:** the paper computes *"the whole of the Tamarkin–Tsygan calculus for the class of
   gentle algebras … the Hochschild cohomology, its structure both as a graded commutative
   algebra under the cup product and as a graded Lie algebra under the Gerstenhaber bracket,
   together with the Hochschild homology and **its module structure over the Hochschild
   cohomology** given by the cap product."* The **module-over-`HH¹`** structure this plan
   computes is one facet of that TT calculus; CSSS is the gentle-family literature anchor.
   `# PIN`: re-confirm J. Algebra 708 (2026) 138–231 (authors/title/abstract verified; the
   volume/pages I take from R12's "independently confirmed" note and re-verify at citation-add).

**Reused / algorithmic references (verify live in `registry.py` at implementation):**

4. **Gerstenhaber, M. — "The cohomology structure of an associative ring."** Ann. of Math. (2)
   78 (1963), 267–288 (P70 already added `gerstenhaber1963`). The Gerstenhaber bracket; **in
   degree `(1,n)` it is exactly the Lie derivative `L_D f`** — the theorem that makes the
   primary route the genuine Gerstenhaber Lie-module action (see §1 for the sign proof).

5. **de Graaf, W. A. — "Lie Algebras: Theory and Algorithms."** North-Holland Math. Library 56,
   Elsevier, 2000 (P70 already added `degraaf_lie`). The Cartan-subalgebra algorithm (the
   maximal-torus computation, §2) and the module-decomposition context.

**Already shipped and reused:** P70's `rss_hh1_lie`, `eisele_raedschelders`, `strametz_hh1_lie`,
`gerstenhaber1963`, `degraaf_lie`; the Plan-35 product citations (`gerstenhaber`, `bracket`) for
the in-window cross-engine oracle; `modules/decompose.py`'s Dickson/Cohen–Ivanyos–Wales locality
references (Plan 30). **New keys to add:** `als_toupie_hh_lie` (arXiv:1803.10310),
`mnprs_special_biserial` (arXiv:1803.10909, J. Algebra 580 (2021) 264–298), `csss_gentle_tt`
(arXiv:2311.08003).

---

## Mathematical foundation (the definitions, exactly — the plan's ground truth)

Throughout `A` is a finite-dimensional associative unital `k`-algebra with structure-constant
tensor `T` (`A.T`), `d = A.dim`, `k = A.domain` an exact `Domain`. Composition is left-to-right.
`C^n = Hom_k(A^{⊗n}, A)` is the Hochschild cochain space, `δ` the Hochschild differential

```
   (δf)(a_0,…,a_n) = a_0·f(a_1,…,a_n)
                     + Σ_{i=1}^{n} (−1)^i f(a_0,…,a_{i-1}a_i,…,a_n)
                     + (−1)^{n+1} f(a_0,…,a_{n-1})·a_n,
```

`HH^n = ker δ^n / im δ^{n-1}`. `HH^0 = Z(A)`, `HH^1 = Der(A)/Inn(A)` (P70). The shipped
field-general substrate is the **normalized** bar cochain complex `Hom_k(Ā^{⊗n}, A)`,
`Ā = A/k·1` (`hochschild/bar.py::_cochain_basis` / `coboundary_matrix`, over any Domain — the
same complex `A.hochschild_cohomology` reports dims on); it computes the same `HH^n` and is
smaller. Everything below is representative-independent on `HH^n` and holds over the normalized
complex verbatim, with the one caveat that inserted derivation values are reduced mod `k·1`
(exactly the shipped circle product's normalization; §1).

### 1. The action = the degree-1 Gerstenhaber bracket = the Lie derivative (field-general)

The Gerstenhaber bracket `[·,·]: C^p ⊗ C^q → C^{p+q-1}` is the graded commutator of the circle
product `f ∘ g`. **The shipped in-house convention** (`engine/tt_calculus.py::circle_cochain`,
verified against source this round) is

```
   (f ∘ g)(a_1,…,a_{p+q-1}) = Σ_{i=0}^{p-1} (−1)^{(q-1) i}
        f(a_1,…,a_i, g(a_{i+1},…,a_{i+q}), a_{i+q+1},…,a_{p+q-1}),
   [f, g] = f ∘ g − (−1)^{(p-1)(q-1)} g ∘ f,
```

with `g`'s value reduced mod `k·1` before it is fed to `f` (`f` normalized).
**Specialize to `p = 1` (`f = D`, a derivation / 1-cochain) and `q = n` (`g` an `n`-cochain):**

- `D ∘ g` has a single term (`i = 0`, sign `(−1)^{(n-1)·0} = 1`): `(D∘g)(a₁…aₙ) = D(g(a₁…aₙ))`.
- `g ∘ D` has terms `i = 0…n−1`, **every sign `(−1)^{(1-1)i} = 1`**:
  `(g∘D)(a₁…aₙ) = Σⱼ g(a₁,…,D(aⱼ),…,aₙ)`.
- the bracket sign is `(−1)^{(1-1)(n-1)} = 1`, so
  **`[D,g] = D∘g − g∘D`**, i.e.
  ```
     (L_D g)(a_1,…,a_n) = D(g(a_1,…,a_n)) − Σ_{i=1}^{n} g(a_1,…,D(a_i),…,a_n).
  ```

**All Koszul signs collapse to `+1` because `p = 1`.** So the field-general Lie derivative `L_D`
is **byte-identical to the shipped `gerstenhaber_bracket_cochain(alg, 1, n, D, g)`** — the #1 sign
risk the card flags is discharged by this source-level identity (and re-confirmed in-window by
the Task-4 arbiter). `L_D` needs `D` only as a `k`-linear map `A → A`, exactly P70's
`derivations(A)`; **no circle-product/bracket engine and no window are needed on this route.**

**Well-definedness on `HH^n`.** `L_D = [D, −]` is a chain map of degree 0 on `(C^•, δ)` (the
Gerstenhaber bracket is compatible with `δ`): `δ(L_D g) = L_D(δ g)`, so `L_D` maps cocycles to
cocycles and coboundaries to coboundaries, and the induced `ρ_n(D): HH^n → HH^n` is independent
of the cochain representative of both `[D] ∈ HH¹` and `[g] ∈ HH^n`. In particular:

- **Inner derivations act as zero.** `[ad_x] = 0` in `HH¹` (`ad_x = −δ^0(x)` is a Hochschild
  coboundary), and the bracket of a coboundary with a cocycle is a coboundary, so `ρ_n(ad_x) = 0`
  on `HH^n`. **The action factors through `HH¹ = Der/Inn`** — pinned as a self-cert.
- **Lie-module axiom.** `ρ_n: HH¹ → gl(HH^n)` is a Lie-algebra homomorphism:
  `ρ_n([D,E]) = ρ_n(D)ρ_n(E) − ρ_n(E)ρ_n(D)`, where `[D,E] = D∘E − E∘D` is P70's commutator of
  derivations. **This makes `HH^n` a module over the Lie algebra `HH¹`**, and
  `HH• = ⊕ₙ HH^n` a **graded** module (the action preserves `n`: degree `1 + n − 1 = n`).

**The primary computation** (per `(A, n)`, any Domain): build `δ^{n-1}, δ^n` via
`bar.coboundary_matrix`; `Z^n = nullspace(δ^n)`, `B^n = columnspace(δ^{n-1})`, and a complement
basis `reps` of `B^n` inside `Z^n` (the `HH^n` representatives, `dim HH^n = len(reps)`). For each
`D` in a `Der/Inn` basis of `HH¹` and each rep `r_j`, `L_D r_j ∈ Z^n` (assert — the descent
check), and its class coordinates give the `j`-th column of `ρ_n(D)`. Reduce inserted derivation
values mod `k·1` on the normalized complex (the shipped convention).

### 2. Weight / torus decomposition (char-0 gate) — NET-NEW code, P70 does not provide it

**This is the plan's hardest deliverable and it is net-new: P70 provides NONE of the machinery
it needs.** P70's `HH1LieStructure` exposes only integers/bools (`radical_dim`, `levi_dim`,
`sl2_count`, `toral_rank`, `levi_type`) — **no Cartan basis, no ad-semisimple elements, no torus
matrices** — and P70 *explicitly scopes OUT* the very torus P71 needs: its honest-scope entry (c)
reads "the broader maximal torus of `HH¹` including ad-semisimple elements of `rad(L)` is out of
scope — the grading element of `k[x]/(x^n)` is such a toral element in the radical." So P71 must
compute the maximal torus **from scratch** (the algorithm is spelled out in Task 3). It reuses
only P70's `char == 0` gate (as a *predicate*) and P70's `derivations`/`inner_derivations`
representation; the torus, its ad-semisimple part, the radical-toral extension, and the
simultaneous diagonalization are all new.

Over `char k = 0`, compute a **maximal torus** `t ⊆ HH¹` — a maximal ad-diagonalizable abelian
subalgebra (Task 3's from-scratch algorithm: a Cartan subalgebra by de Graaf split-element/Fitting,
then its ad-semisimple part, then the radical-toral extension). Because `t` is abelian and `ρ_n` a
Lie hom, `{ρ_n(h) : h ∈ t}` is a **commuting family**; when each `ρ_n(h)` is diagonalizable over
`k` they are **simultaneously diagonalizable**, giving the **weight-space decomposition**

```
   HH^n = ⊕_{λ ∈ t*} HH^n_λ,   HH^n_λ = { v : ρ_n(h) v = λ(h) v  ∀ h ∈ t },
```

with weights `λ` recorded as integer/rational tuples over a chosen basis of `t` and per-weight
dims. `Σ_λ dim HH^n_λ = dim HH^n` (self-cert).

- **Torus provenance (Plan-35 precedent, binding).** The torus `t` is **basis-dependent** (the
  null-space / Cartan pivots, and the normalization of each `h`, are choices); the report tags
  the torus with its provenance and the normalization used. The weight **multiset up to the
  torus normalization** (e.g. sl₂ roots normalized to `±2`; a grading normalized so the
  degree-1 generator has weight `1`) is the reproducible content — the exact integer labels are
  reported but **NON-NORMATIVE** for tests (the sl₂-triple precedent of P70) — a scalar-doubled
  torus generator `2h` turns `{−2,0,2}` into `{−4,0,4}`, so the literal integers are provenance,
  never an assertion; the **normative** weight facts are basis-independent: the number of distinct
  weights, the per-weight dims, and the weight *pattern* (e.g. "a symmetric sl₂-string `{−c,0,c}`,
  `c ≠ 0`", "an equal-gap arithmetic progression").
- **SCOPE — the torus is a DIFFERENT object from P70's `toral_rank` (net-new, not an extension of
  P70 scaffolding).** P70's `toral_rank` is the rank of the **semisimple Levi factor** only, and
  P70 explicitly scopes OUT "the maximal torus of `HH¹` including ad-semisimple elements of the
  radical." P71's weight torus is exactly that scoped-out object — the maximal ad-diagonalizable
  abelian subalgebra of **all** of `HH¹` (the grading / Euler operator `x∂` of `k[x]/(x^n)` is a
  semisimple element of the *solvable radical*, and MNPRS's "Virasoro degree" is precisely such a
  grading). P71 computes it from scratch (Task 3), reads NOTHING from P70's `toral_rank`, and never
  re-derives or contradicts that field.
- **Honest three-valued (weights only).** (i) **no nontrivial torus** (`HH¹` nilpotent, or `t =
  0`): weights are all `0` and the decomposition is trivial — reported as `weights = None` with a
  note. (ii) **`t` does not act semisimply over `k`** (anisotropic: an eigenvalue lies outside
  `k`): weights are `"unavailable over k — base change to k̄ needed"` (the P70 base-change-note
  precedent), never a fabricated split. (iii) **char > 0**: the **weight/torus** block is a loud
  `QuiverlabError` behind the P70 char-0 gate (no Cartan/Weyl theory); the field-general action
  matrices + module axiom are still returned, and the **indecomposable-summand decomposition is
  governed independently** by `decompose`'s own char guard (`char 0 or char > d_n`; §3–§4), NOT by
  this weight gate.

### 3. `HH^n` as a module over `HH¹` — the associative envelope (the decompose input)

The Lie-module structure of `HH^n` is the same data as the module structure over the **associative
envelope**

```
   B_n := ⟨ ρ_n(D) : D ∈ HH¹ ⟩ + k·I   ⊆   M_{d_n}(k),      d_n = dim HH^n,
```

the unital associative subalgebra of `M_{d_n}` generated by the action matrices — equivalently the
image of the universal enveloping algebra `U(HH¹) → End(HH^n)`. **Why this is the right object
(justification, live-verified in §5):**

- A subspace `W ⊆ HH^n` is `HH¹`-Lie-invariant **iff** it is `B_n`-invariant (both mean stable
  under every `ρ_n(D)`), so `HH^n` has the **same invariant-subspace lattice** as a Lie-module and
  as a `B_n`-module.
- `End_{B_n}(HH^n) = { φ : φ ρ_n(D) = ρ_n(D) φ  ∀ D } = ` the **commutant** of `{ρ_n(D)}` `=
  End_{HH¹}(HH^n)`. The endomorphism ring is the same object.
- Therefore the **Krull–Schmidt decomposition of `HH^n` as a `B_n`-module coincides with its
  decomposition as an `HH¹`-Lie-module** (Krull–Schmidt is governed by the local structure of
  `End`). `B_n` is finite-dimensional (`⊆ M_{d_n}`) and exactly constructible (matrix closure).

This is the concrete answer to the card's "state how `HH^n` becomes a Module over what Algebra":
**over `B_n`, the associative subalgebra generated by the action matrices** (NOT an abstract
`U(g)`-truncation — `B_n` IS the finite image of `U(HH¹)` and is directly constructible).

### 4. Reusing `modules/decompose.py` on the representation (mechanics, exactly)

`modules/decompose.py::decompose(M)` (Plan 30) splits a module by Fitting's method: it computes
`End_A(M) = hom_space(M, M)`, scans candidate endomorphisms for a coprime min-poly factorization
`M = ker f(φ) ⊕ ker g(φ)`, and, when no split is found, certifies indecomposability via
`End_A(M)` local (`dim End = 1 ⟹ local`, rigorous in **every** characteristic; else the
Dickson/Cohen–Ivanyos–Wales trace-form rank under the `char 0 or char > dim M` bound, loud
refusal at `char p ≤ dim`). **All of this machinery is matrix-level** — but the public
`decompose(M)` reads `M.algebra.quiver.vertices/arrows` (via `hom_space`), so it cannot be fed a
bare matrix algebra `B_n` (which has no quiver). **P71 therefore factors the Fitting core out:**

- **Add** `modules/decompose.py::decompose_representation(gens, domain, *, budget=...) ->
  [Summand]` — the min-poly split + local-endomorphism certificate driven by `End = ` the
  **directly computed commutant** of `gens` (exact `fields.linalg` nullspace of the simultaneous
  commutation system `X gᵢ − gᵢ X = 0`), returning certified-indecomposable summands as
  `(subspace basis, restricted action matrices)`. The existing `decompose(M)` becomes a thin
  wrapper (module → its generator action matrices → `decompose_representation`), **byte-stable**
  (a golden-pinned regression: existing decompositions unchanged).
- P71 calls `decompose_representation([ρ_n(D₁),…,ρ_n(D_r)], A.domain)` (the `Der/Inn` basis
  action matrices). Because `End = ` commutant `= End_{HH¹}(HH^n)`, the certified summands are the
  **indecomposable Lie-module summands** of `HH^n` (§3). The **same char guard** applies (char-0
  cleanest; `char > d_n` OK; loud refusal at `char p ≤ d_n`).

**Live-verified that decompose's certificate fires on the representation (§5):** `kK₂`'s
`HH¹ = sl₂` acting on `HH^1` (the adjoint) has `dim End = 1` → the "`dim End = 1 ⟹
indecomposable`" certificate (rigorous every char) returns **one** irreducible summand;
`k[x,y]/(x,y)²`'s `HH^n` have `dim End = 2, 2, 6 > 1` → the multi-summand trace-form path fires.

### 5. The pinned examples, exactly (all recomputed in the venv, §Live-verified)

- **`kK₂`** (Kronecker, `1 ⇉ 2`, hereditary, `dim 4`): `HH• = [1, 3, 0, 0, …]`. `HH¹ ≅ sl₂`
  (P70) `= HH^1` = the **irreducible 3-dimensional adjoint module `L(2)`**, torus weights
  a symmetric sl₂-string `{−c, 0, c}` (`c ≠ 0`; `c = 2` under the roots-normalized Cartan);
  `HH^0 = k` = the trivial module `L(0)` (weight `{0}`); `HH^{≥2} = 0`. The sl₂ case, exactly
  what the card asks ("which irreducibles appear? pin dims + highest weights": `L(0)` in degree 0,
  `L(2)` in degree 1). **`kK₂` IS a toupie algebra** (a single source `1`, single sink `2`, two
  parallel paths — the simplest toupie), so its `HH•` Lie-module structure is a member of the ALS
  toupie family (`als_toupie_hh_lie`) — a REAL literature anchor, not just an analogy.
- **`k[x]/(x^n)`** (`char ∤ n`): `HH• = [n, n−1, n−1, …]`; `HH¹ = ` the **positive-Witt
  truncation** `⟨x^{i+1}∂ : 0 ≤ i ≤ n−2⟩`, `[x^{i+1}∂, x^{j+1}∂] = (j−i) x^{i+j+1}∂` — a
  **subquotient of the Witt (centreless Virasoro) algebra**, which is an ANALOGUE of the MNPRS
  Virasoro-subquotient structure (MNPRS's theorem proper is about their *special-biserial* family,
  §Reference 2 — `k[x]/(x^n)` is a self-injective Nakayama, so the analogy is close but the pin is
  self-contained, not a citation of MNPRS's result). Solvable for `n ≥ 3`. Torus `t = k·(x∂)` (the
  grading Euler operator — a *radical* toral element, §2 scope note). Each `HH^n` (`n ≥ 1`, dim
  `n−1`) is a **single indecomposable** `HH¹`-module (`dim End = 1`), with two distinct grading
  weights an equal gap apart (`k[x]/x³` under the `x`-degree normalization: `{0,1}, {−2,−3},
  {−2,−3}, {−5,−6}` in degrees `1,2,3,4` — NON-NORMATIVE labels, §2).
- **`k[x,y]/(x,y)²`** (2-loops, `rad² = 0`, non-hereditary, `dim 3`): `HH• = [3, 4, 6, 12, …]`,
  nonzero in every degree. Commutative, so `Inn = 0` and `HH¹ = Der = End(rad) ≅ gl₂` (dim 4,
  reductive). Torus `t = ` diagonal of `gl₂` (2-dim). The `HH^n` are `gl₂`-modules that
  **decompose into ≥ 2 indecomposables** (`dim End = 2, 2, 6` in degrees `1, 2, 3`) — the
  multi-summand example.
- **`k[x]/(x²)`** (dual numbers): `HH• = [2, 1, 1, …]`, `HH¹ = k·(x∂)` (1-dim abelian). Each
  `HH^n` (`n ≥ 1`) is a 1-dim weight space, weight ladder `{0}, {−2}, {−2}, {−4}` (degrees
  `1..4`) — the cleanest weight illustration.

---

## Live-verified facts (all recomputed in the venv at spec time)

Two independent computations agree on every row: **(A)** a **standalone exact Hochschild cochain
complex** over `QQ` (sympy `Rational`) — the unnormalized `Hom(A^{⊗n}, A)` with the textbook
differential, cocycles/coboundaries by exact nullspace/columnspace, `Der(A)`/`Inn(A)` by the
Leibniz null space, the Lie derivative `L_D` per §1, the induced action `ρ_n(D)`, the module
axiom, the inner-acts-zero check, the torus weights, the associative envelope `B_n` and the
commutant `End_{HH¹}(HH^n)`. **(B)** `A.hochschild_cohomology(top=…).dims` (the shipped bar/CS
engine) for the `HH^n` dims, built in quiverlab (`Quiver(...).algebra`, the 2-loop `rad²=0`
quotient). The shipped `gerstenhaber_brackets` degree-`(1,n)` surface was confirmed present.

**Oracle 1 — `HH^n` dims (cross-engine anchor: standalone probe ≡ `A.hochschild_cohomology`).**

| algebra | field | `HH•` dims (probe) | `A.hochschild_cohomology` | agree |
|---|---|---|---|---|
| `k[x]/(x²)` | `QQ` | `[2,1,1,1,1]` | `[2,1,1,1,1]` | ✓ |
| `k[x]/(x³)` | `QQ` | `[3,2,2,2,2]` | `[3,2,2,2,2]` | ✓ |
| `k[x,y]/(x,y)²` | `QQ` | `[3,4,6,12]` | `[3,4,6,12]` | ✓ |
| `kK₂` | `QQ` | `[1,3,0,0]` | `[1,3,0,0,0]` | ✓ |

**Oracle 2 — the Lie-module action `ρ_n` (self-cert: axiom + inner-zero, every row).**
For each algebra, at every degree `n` with `HH^n ≠ 0`: the module axiom
`ρ_n([D,E]) = [ρ_n(D), ρ_n(E)]` holds on all `Der/Inn` basis pairs, and every **inner** derivation
acts as the **zero** matrix on `HH^n` (the action factors through `HH¹ = Der/Inn`). Verified for
`k[x]/(x²)` (`n=1..4`), `k[x]/(x³)` (`n=1..4`), `k[x,y]/(x,y)²` (`n=1..3`), `kK₂` (`n=1`). The
`L_D` used is the §1 formula = the shipped `(1,n)` bracket (source-identical signs).

**Oracle 3 — weights under an explicit maximal torus (char-0; NON-NORMATIVE exact labels).**

| algebra | torus `h` (normalization) | `HH^0` | `HH^1` | `HH^2` | `HH^3` | `HH^4` |
|---|---|---|---|---|---|---|
| `k[x]/(x²)` | `x∂` (grading) | `{0,1}` | `{0}` | `{−2}` | `{−2}` | `{−4}` |
| `k[x]/(x³)` | `x∂` (grading) | `{0,1,2}` | `{0,1}` | `{−2,−3}` | `{−2,−3}` | `{−5,−6}` |
| `k[x,y]/(x,y)²` | `(x∂, y∂)` (diag `gl₂`) | `(0,0)²,(1,1)`* | 4 wts* | 6 wts* | 12 wts* | — |
| `kK₂` | sl₂-Cartan `h` (roots `±2`) | `{0}` | `{−2, 0, 2}` | 0 | 0 | — |

\* `k[x,y]/(x,y)²` weights are 2-tuples; the marginal `x∂`-weights are, per degree: `HH^0:
{0²,1}`; `HH^1: {−1,0²,1}` (= `gl₂` adjoint); `HH^2: {−1²,−2,0²,1}`; `HH^3:
{−1⁴,−2³,−3,0³,1}` (`y∂` marginals identical by symmetry). The integers above are what these
SPECIFIC torus normalizations (roots-`±2` sl₂-Cartan; the `x`-degree grading) produced this round.
**Normative (test-asserted) content is the basis-independent PATTERN, NOT these integers** (a
scalar-doubled generator sends `{−2,0,2}` → `{−4,0,4}`): `kK₂ HH^1` is a symmetric **sl₂-string
`{−c,0,c}`** (`c ≠ 0`, the adjoint `L(2)`); `k[x]/(x^n) HH^n` has an **equal-gap** grading-weight
progression (gap `1` only under the `x`-degree normalization); per-weight dims and distinct-weight
counts are normative (the P70 sl₂-triple NON-NORMATIVE precedent).

**Oracle 4 — the associative envelope `B_n` and the decompose certificate (the reuse validation).**

| algebra | `n` | `dim HH^n` | `dim B_n = ⟨ρ_n(HH¹),I⟩` | `dim End_{HH¹}(HH^n)` (commutant) | decompose verdict |
|---|---|---|---|---|---|
| `kK₂` | 1 | 3 | **9 = full `M₃`** | **1** | `dim End=1` ⟹ **1 irreducible** (adjoint `L(2)`) |
| `k[x]/(x³)` | 1..4 | 2 | 3 (triangular ⊂ `M₂`) | **1** | `dim End=1` ⟹ **1 indecomposable** (uniserial) |
| `k[x]/(x²)` | 1..4 | 1 | 1 | **1** | `dim End=1` ⟹ **1** (the line) |
| `k[x,y]/(x,y)²` | 1 | 4 | 10 | **2** | `End>1` ⟹ trace-form/Fitting → **multi-summand** |
| `k[x,y]/(x,y)²` | 2 | 6 | 20 | **2** | multi-summand |
| `k[x,y]/(x,y)²` | 3 | 12 | 35 | **6** | multi-summand |

The `kK₂` row is the card's requested validation: `B_1 = M₃` (Burnside — the adjoint is
absolutely irreducible in char ≠ 2), `End = k·I` (dim 1), so `decompose`'s rigorous
`dim End = 1 ⟹ indecomposable` certificate returns the single irreducible summand. `k[x,y]/(x,y)²`
exercises the `End > 1` multi-summand path (char-0 trace-form).

**Sign-convention arbiter (source-level, this round).** `engine/tt_calculus.py::circle_cochain`
+ `gerstenhaber_bracket_cochain` were read: for `p = 1` the `(−1)^{(q-1)i}` and `(−1)^{(p-1)(q-1)}`
signs are both `+1`, so `[D, g] = D∘g − g∘D` with the §1 formula — **identical** to the plan's
field-general `L_D`. `kK₂.gerstenhaber_brackets(top=2)` over `GF(32003)` returned tables for
degrees `(1,1),(1,2),(2,1)` (the `(1,n)` component exists) — the in-window arbiter surface is live.

**What I could NOT live-verify (labelled to-be-verified-at-implementation):**
- **The `decompose_representation` refactor as shipped** — I verified the *equivalence* (`B_n` /
  commutant / `dim End` above) that decides the certificate, and that `decompose(M)` is
  quiver-coupled (so the refactor is required), but not the shipped code path. The pinned counts
  (kK₂ HH¹ = 1 irreducible; `k[x,y]/(x,y)²` multi-summand) are the standing oracles.
- **The general maximal-torus (Cartan-subalgebra) computation** — I used *explicit* torus
  elements (`x∂`, the sl₂-Cartan, the `gl₂` diagonal) hand-identified per example; the general de
  Graaf Cartan computation over `QQ` is the worker's Task-3 work, with the §2 three-valued
  fallbacks (trivial torus / anisotropic / char-p) keeping it honest.
- **The MNPRS / ALS literature decomposition tables at scale** — the pinned small cases
  (`kK₂`, `k[x]/(x^n)`, `k[x,y]/(x,y)²`) reproduce the *structure* (Virasoro-subquotient `HH¹`,
  indecomposable `HH^n`), but the exact MNPRS special-biserial / ALS toupie decomposition
  **tables** for a specific published member are `# PIN`'d to transcribe verbatim at
  implementation (Task 5), not invented here.
- **GAP module cross-check** — QPA has no Hochschild Lie-module surface (Plan-35 precedent); GAP's
  core Lie library can crosscheck the module decomposition. **Honest caveat:** the sl₂ = `HH¹(kK₂)`
  case is a *trivial* cross-check (the adjoint is irreducible → 1 summand — GAP agreeing there
  proves little); the oracle's REAL value is the **multi-summand** `k[x,y]/(x,y)²` case (the `gl₂`
  action on `HH^n` with `dim End = 2, 2, 6` → ≥ 2 indecomposables), where GAP's
  `DirectSumDecomposition` on the fed representation is a genuine independent count. I confirmed the
  access path but did NOT run GAP.

---

## Scope gates (contractual; every boundary is a loud typed refusal)

| surface | scope | refusal |
|---|---|---|
| `lie_derivative_on_hh` / `lie_module_action` (dims, `ρ_n`, module axiom, inner-zero) | **any exact field**; `A.dim ≤ budget`; `top` within the bar cell cap (`max_cells`) | oversized `A.dim > budget` → `QuiverlabError` (Der size + budget); a `top` whose bar cochain basis `d·(d−1)^{top+1}` exceeds `max_cells` → the same loud bar-cell refusal `A.hochschild_cohomology` already raises (reused, not reinvented) |
| **weight / torus decomposition ONLY** | **`char == 0` only** (Cartan/Weyl theory is char-0; no positive-char analogue) — the P70 char-0 predicate | `char > 0` → loud `QuiverlabError` (`require_char0=True`) or `weights = None` + `char0_note`; the field-general `ρ_n` / axiom block AND the summand decomposition (governed by the row below, INDEPENDENTLY) are still returned |
| **indecomposable-summand decomposition** (INDEPENDENT of the weight gate) | **`char == 0 or char > d_n = dim HH^n`** — `decompose`'s OWN Dickson/CIW guard, applied per degree; NOT tied to the weight/torus char-0 gate | `char p ≤ d_n` and no `dim End = 1` certificate → `decompose_representation` raises loudly (the Plan-30 refusal), surfaced as a per-degree note (never a crash, never a guessed decomposition). So at `char p > d_n` summands ARE returned even though weights are `None` |
| weights over a torus that does not split over `k` (anisotropic) | reported with **base-change provenance** (gated on the arithmetic field, the P70 precedent) | never refused; weights are `"unavailable over k — base change to k̄"`, no fabricated split |
| all | `A.dim ≤ DEFAULT_MAXDIM = 48` (inherited from P70) — but for `top ≥ 2` the **real limiter is `max_cells`**: the bar cochain basis `d·(d−1)^{top+1}` blows up exponentially in `top` while the Der solve is fixed by `A.dim`, so `max_cells` (not `DEFAULT_MAXDIM`) is what refuses first as `top` grows | `A.dim > 48` OR an over-`max_cells` `top` → honest `status="budget"` refusal (the `A.hochschild_cohomology` bar-cell wall, reused), never a partial guess; the dim-220 Nakayama webapp examples carry NO `hh_lie_module` (the Plan-35 products-omission precedent on the same examples) |

Action **structure constants** always carry `basis="der_inn/bar"` provenance; cross-engine and
cross-run comparisons use only **basis-independent** data (dims, weight multisets, the induced-map
invariants char-poly/rank/trace, `dim End`, summand dims) — the Plan-35 rule (never compare raw
constants across bases). The **torus choice** is basis-dependent provenance (§2).

---

## API surface (public via `import quiverlab`; exact only)

```python
# src/quiverlab/hochschild/lie_module.py
def lie_derivative_on_hh(A, D, n, *, max_cells=4_000_000) -> list   # rho_n(D): d_n x d_n over A.domain
def lie_module_action(A, top, *, budget=DEFAULT_MAXDIM, max_cells=4_000_000) -> HHLieModule
def hh_lie_module_block(A, top, *, budget=...) -> dict              # the runner block

@dataclass(frozen=True)
class HHLieModule:
    top: int
    hh_dims: list                # [dim HH^0, ..., dim HH^top]  (== A.hochschild_cohomology)
    hh1_dim: int                 # dim HH^1 (the acting Lie algebra; P70)
    basis: str                   # provenance: "der_inn/bar"
    action: list                 # per n: {"n": n, "gens": [flattened rho_n(D) as EXACT strings], ...}
    module_axiom_ok: bool        # rho_n([D,E]) == [rho_n(D), rho_n(E)] on all pairs (self-cert)
    inner_acts_zero: bool        # rho_n(ad_x) == 0 (action factors through Der/Inn) (self-cert)
    characteristic: int
    # WEIGHTS: char-0 ONLY (None over char p, with char0_note):
    weights: list | None         # per n: {"n": n, "torus_rank": int, "weights": [[lam..],dim], ...}
    torus_provenance: str | None # which maximal torus + its normalization (basis-dependent)
    # SUMMANDS: governed INDEPENDENTLY by decompose's own guard (char 0 OR char > d_n), NOT the
    # weight gate -- so `summands` can be non-None at char p > d_n while `weights` is None; a
    # per-degree entry is None+note only where decompose refuses (char p <= d_n, no dim-End=1 cert):
    summands: list | None        # per n: [{"dim": int, "label": str|None, "multiplicity": int}] | {"error": note}
    weight_base_change_note: str | None   # anisotropic-torus / arithmetic-field note (P70 precedent)
    char0_note: str | None       # the WEIGHT/torus loud-gate explanation when char > 0
    status: str                  # "complete" | "budget" | "unsupported"
    note: str
    references: list

# modules/decompose.py (new matrix-level entry; decompose(M) delegates to it — byte-stable):
def decompose_representation(gens, domain, *, budget=_DEFAULT_BUDGET) -> list   # certified summands

# Algebra delegate (core/algebra.py, thin lazy-import beside gerstenhaber_brackets):
Algebra.hh_lie_module(top, budget=..., max_cells=...) -> HHLieModule
```

`HHLieModule` is a frozen data report (like Plan-35's `HHProducts`); `__bool__` is NOT defined.

---

## Global Constraints

- Python always `.venv/bin/python`; tests `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2
  .venv/bin/python -m pytest -q ...`.
- **Reuse, do not reinvent (verified live in `dev` at spec time):**
  - **P70 `invariants/hh1_lie.py`**: `derivations(A)` (a `Der(A)` basis as `d×d` matrices over
    `A.domain`), `inner_derivations(A)` (the `Inn` ideal), `hh1_lie_structure(A)` (the char-0
    classification — radical/Levi/`sl2_count` as ints/bools). Honor its `DEFAULT_MAXDIM = 48`, its
    `char == 0` gate (as a *predicate*), and its arithmetic-field base-change gating (do NOT
    re-derive them). **P70 provides NO torus matrices / Cartan basis** — the maximal torus is
    NET-NEW code in Task 3 (§2); P71 reads only the reps + the char-0 predicate from P70.
  - **`hochschild/bar.py`**: `coboundary_matrix(A, n, max_cells)` (the normalized Hochschild
    cochain differential over any Domain — the field-general substrate), `_cochain_basis`,
    `hochschild_cohomology_dims` (the dim anchor). `cyclic.py` is the generic-`(b,B)` precedent.
  - **`quiverlab.fields.linalg`**: `nullspace`, `rank`, `rref`, `solve` (exact over any Domain).
  - **`A.hochschild_cohomology(top)`** — the cross-engine `HH^n` dim anchor.
  - **`A.gerstenhaber_brackets(top, engine=...)`** (Plan 35 / P51) — the degree-`(1,n)` in-window
    SIGN ARBITER (GF(p) in-window; NOT the primary route).
  - **`modules/decompose.py`** — the Fitting-split + local-endomorphism (`dim End = 1`;
    Dickson/CIW trace-form, `char 0 or char > dim`) machinery; P71 exposes its matrix-level core.
- **The char-0 gate is P70's `char == 0`** (weight/torus/Levi/decomposition-cleanest); the
  field-general `ρ_n` / module-axiom / graded-module surface stays live in **every**
  characteristic. Document at the gate; do NOT re-argue P70's rationale (no positive-char
  Cartan/Weyl), reference it.
- **No floats in `src/`.** Action structure constants are exact-`Domain` strings; dims/weights
  `int`; verdicts `bool`/`str`. Client-side numeric conversion only (`gui.js`, exempt).
- **Composition is left-to-right**; the Lie derivative and `δ` respect `A.multiply`. All refusals
  are `QuiverlabError`.
- **Plan-32 markers** (orthogonal; `oracle_*` FORBIDDEN in `tests/{webapp,gui,hpc}`):
  - `oracle_literature`: `kK₂ HH^1 = ` irreducible sl₂-adjoint `L(2)`, weights a symmetric
    sl₂-string `{−c,0,c}` (`c ≠ 0`; the PATTERN, not the literal integers — §2 NON-NORMATIVE) with
    all per-weight dims 1 (+ `HH^0 = L(0)`); `k[x]/(x^n) HH^n` a single indecomposable over the
    truncated-Witt `HH¹`, weights an equal-gap progression; the `k[x]/(x²)` weight ladder;
    `k[x,y]/(x,y)²` `HH^n` multi-summand `gl₂`-modules; the MNPRS special-biserial / ALS toupie
    decomposition anchor (`# PIN` verbatim table for the pinned small member, Task 5 — plus `kK₂` as
    a live ALS-toupie member and `k[x]/(x^n)` as the self-contained truncated-Witt anchor).
  - `oracle_crossengine`: `hh_dims == A.hochschild_cohomology` on the zoo; the **STRONG entry-wise
    sign check** — the `L_D` operator built on the engine's degree-`n` cochain basis equals the
    shipped `engine.tt_calculus.gerstenhaber_bracket_cochain(alg,1,n,D,·)` ENTRY-WISE on the whole
    cochain space over GF(p) (strictly stronger than an induced-map invariant; the critic ran this
    and it passes); AND the weaker in-window induced-map arbiter (`ρ_n` vs `gerstenhaber_brackets`
    degree-`(1,n)` on basis-independent char-poly / rank / trace) for the field-general/basis-
    independent comparison; the normalized-bar route ≡ the unnormalized-probe structure degreewise.
  - `oracle_selfcert`: module axiom `ρ_n([D,E]) = [ρ_n(D),ρ_n(E)]`; `ρ_n(ad_x) = 0` (factors
    through `Der/Inn`); `ρ_n` maps cocycles to cocycles (representative-independent); torus
    elements commute (`[ρ_n(h),ρ_n(h')] = 0`); `Σ_λ dim HH^n_λ = dim HH^n`; each certified summand
    reassembles `HH^n`; the char-p loud gate on **weights/torus** (the summand decomposition is
    governed independently, §3–§4); the oversize/`max_cells` loud refusal.
  - `qpa`: QPA has NO Hochschild Lie-module surface (an honest skip that FAILS if it appears, the
    Plan-35 precedent); the substantive GAP oracle is the **multi-summand** `k[x,y]/(x,y)²` case
    (`gl₂` acting, `dim End = 2,2,6` → ≥ 2 indecomposables), where GAP's `DirectSumDecomposition`
    on the fed representation is a genuine independent count (the sl₂ = `kK₂` case is a trivial
    irreducible check, kept only as a sanity sentinel).
- **Mid-merge-train counts.** The verification task recounts the oracle-class table at merge time
  (`tests/release/test_oracle_classes.py`; paste live numbers, claim only this plan's deltas).
- Every merge updates `docs/verification.md` (new oracle rows + recounted class table) and adds
  citations to `references.bib` + `registry.py` (`bibtex()` hard-fails if the two disagree).
  Conventional commits; green at every commit.

**No-red-commits discipline.** Each Task ships a self-contained green slice: the action primitive
before the report before the char-0 weight/decomposition before the GUI. Char-0 fields default to
`None` so intermediate task boundaries are green (no forward reference, no `xfail`).

---

## Task 0 — citations (do first; later tasks reference the keys)

**Files:** `src/quiverlab/citations/references.bib` + `registry.py`.

- [ ] **Step 1: Add BibTeX (verified fields only; `# PIN` the rest, never guess).**
```bibtex
@misc{als_toupie_hh_lie,
  author = {Artenstein, Dalia and Lanzilotta, Marcelo and Solotar, Andrea},
  title  = {Gerstenhaber structure on {H}ochschild cohomology of toupie algebras},
  year   = {2018}, note = {arXiv:1803.10310}}                       % journal ref # PIN if any
@article{mnprs_special_biserial,
  author  = {Meinel, Joanna and Nguyen, Van C. and Pauwels, Bregje and Redondo, Mar\'ia Julia and
             Solotar, Andrea},
  title   = {The {G}erstenhaber structure on the {H}ochschild cohomology of a class of special
             biserial algebras},
  journal = {J. Algebra}, volume = {580}, pages = {264--298}, year = {2021}}   % vol/pages ScienceDirect search hit, re-confirm on landing page # PIN
@article{csss_gentle_tt,
  author  = {Chaparro, Cristian and Schroll, Sibylle and Solotar, Andrea and
             Su\'arez-\'Alvarez, Mariano},
  title   = {The {H}ochschild cohomology and the {T}amarkin--{T}sygan calculus of gentle algebras},
  journal = {J. Algebra}, volume = {708}, pages = {138--231}, year = {2026},
  note    = {arXiv:2311.08003}}                                     % vol/pages R12-confirmed # PIN
```
- [ ] **Step 2:** add the `registry.py` entries (`als_toupie_hh_lie` = the Lie-module of `HH•` for
  toupie algebras; `mnprs_special_biserial` = the Virasoro-subquotient `HH¹` and the degree-`n`
  indecomposable decomposition; `csss_gentle_tt` = the gentle-algebra TT calculus, the module-
  over-`HH¹` facet) with annotations naming what each grounds. Verify
  `tests/citations/test_bib_structure.py` green.
- [ ] **Step 3:** commit `docs(citations): P71 HH-Lie-module references (Artenstein-Lanzilotta-Solotar toupie, MNPRS special biserial/Virasoro, CSSS gentle TT-calculus)`.

---

## Task 1 — the action primitive `ρ_n(D)` + graded module (any exact field)

**Files:** create `src/quiverlab/hochschild/lie_module.py`; test `tests/hochschild/test_lie_module.py`.
**`tests/hochschild/` is the FAST bucket** by `tests/conftest.py` (only
`tests/{engine,resolutions_cs,modules,families,batch}/` auto-assign to `deep`) — so keep every test
here small (`top ≤ 3`, `A.dim ≤ 4`, which all the pins are). Any heavier probe (large `top`, dim
> ~8) must go in a **deep-bucketed dir** (e.g. `tests/engine/` or `tests/modules/`) or be honestly
marked — never a slow test masquerading as fast in `tests/hochschild/`.

- [ ] **Step 1: Failing tests.**
```python
# tests/hochschild/test_lie_module.py
"""HH^* as a graded Lie module over HH^1 (Plan 71 / R12; Gerstenhaber deg-1 = Lie derivative).
Field-general on the normalized bar cochain complex; char-0 weights/decomposition behind P70's gate."""
import pytest
from quiverlab.fields import QQ, GF
from quiverlab.combinat.quiver import Quiver
from quiverlab.families import truncated_polynomial
from quiverlab.hochschild.lie_module import lie_module_action, lie_derivative_on_hh
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert

@xeng
@pytest.mark.parametrize("n,dims", [(2,[2,1,1]),(3,[3,2,2])])   # k[x]/x^n : HH^0=n, HH^k=n-1
def test_hh_dims_match_bar_engine(n, dims):
    A = truncated_polynomial(n, field=QQ)
    L = lie_module_action(A, top=2)
    assert L.hh_dims == dims == A.hochschild_cohomology(top=2).dims

@xeng
def test_kronecker_hh_dims():
    A = Quiver([1,2], {"a":(1,2),"b":(1,2)}).algebra(field=QQ)    # kK2
    L = lie_module_action(A, top=2)
    assert L.hh_dims == [1,3,0] and L.hh1_dim == 3

@selfcert
@pytest.mark.parametrize("build", ["trunc3", "rad2loops", "kron"])
def test_module_axiom_and_inner_zero(build):
    A = ({"trunc3": lambda: truncated_polynomial(3, field=QQ),
          "rad2loops": lambda: Quiver([1], {"x":(1,1),"y":(1,1)}).algebra(
                          relations=["x*x","x*y","y*x","y*y"], field=QQ),
          "kron": lambda: Quiver([1,2], {"a":(1,2),"b":(1,2)}).algebra(field=QQ)}[build])()
    L = lie_module_action(A, top=2)
    assert L.module_axiom_ok is True and L.inner_acts_zero is True
```
- [ ] **Step 2:** confirm failure (`ModuleNotFoundError`).
- [ ] **Step 3: Implement** `lie_derivative_on_hh(A, D, n)` (build `δ^{n-1}, δ^n` via
  `bar.coboundary_matrix`; `Z^n = nullspace`, `B^n = columnspace`, `reps` = complement; `L_D` per
  §1 on the normalized cochain basis, inserted derivation values reduced mod `k·1`; assert
  `L_D r ∈ Z^n`; return `ρ_n(D)` in the class basis) and the field-general part of
  `lie_module_action` (`hh_dims`, `hh1_dim` from P70's `derivations`/`inner_derivations`, the
  per-degree `ρ_n(D)` structure constants tagged `basis="der_inn/bar"`, `module_axiom_ok`,
  `inner_acts_zero`). `budget` guards `A.dim` (P70's `DEFAULT_MAXDIM`); `max_cells` guards `top`.
- [ ] **Step 4:** run; **Step 5:** commit `feat(hochschild): HH^* as a graded Lie module over HH^1 -- field-general rho_n(D) = Lie derivative (Plan 71/R12); dims cross-checked vs the bar engine, module axiom + inner-zero self-certified`.

---

## Task 2 — the Gerstenhaber-bracket in-window SIGN ARBITER (cross-engine)

**Files:** modify `lie_module.py`; extend the test (an `oracle_crossengine` battery).

- [ ] **Step 1: Failing tests.**
```python
@xeng
@pytest.mark.parametrize("build,p", [("trunc4", 5), ("kron", 3)])
def test_sign_arbiter_in_window(build, p):
    """The field-general L_D route agrees IN-WINDOW over GF(p) with gerstenhaber_brackets
    degree-(1,n) on BASIS-INDEPENDENT data (induced-map char-poly/rank/trace). The #1 sign
    risk: p=1 collapses all Koszul signs, so L_D == the shipped (1,n) bracket (source-identical).
    Structure constants are NOT compared (basis-dependent, Plan-35 rule)."""
    A = (truncated_polynomial(4, field=GF(p)) if build=="trunc4"
         else Quiver([1,2], {"a":(1,2),"b":(1,2)}).algebra(field=GF(p)))
    L = lie_module_action(A, top=2)
    g = A.gerstenhaber_brackets(top=2)                     # in-window (1,1),(1,2)
    inv = _bracket_action_invariants(g)                    # per (1,n): rank/char-poly of the induced map
    assert _action_invariants(L) == inv                    # basis-independent agreement

@xeng
@pytest.mark.parametrize("build,p,n", [("trunc3", 7, 2), ("trunc3", 7, 3), ("kron", 5, 1)])
def test_L_D_equals_circle_cochain_ENTRY_WISE(build, p, n):
    """STRONG form (strictly stronger than the induced-map arbiter): the L_D operator built on
    the engine's degree-n COCHAIN basis equals the shipped Gerstenhaber bracket cochain-for-cochain
    over GF(p). The critic ran this on k[x]/x^3 over GF(7) and it passes -- the p=1 sign collapse
    is exact, not just up to invariants. (Reaches into engine.tt_calculus, as tests may -- the
    no-reach-into-engine rule is for app/GUI code, not tests: test_gerstenhaber.py already does.)"""
    from quiverlab.engine import tt_calculus as TT
    from quiverlab.engine.adapter import to_engine
    A = (truncated_polynomial(3, field=GF(p)) if build=="trunc3"
         else Quiver([1,2], {"a":(1,2),"b":(1,2)}).algebra(field=GF(p)))
    E = to_engine(A.unit_adapted())
    for D in _der_reps_on_engine_basis(A, E):              # each Der/Inn rep as a 1-cochain over E
        M = _L_D_matrix_on_cochains(E, D, n)               # our L_D on the deg-n cochain basis
        for f in _cochain_basis_vectors(E, n):             # every basis cochain (or a random sample)
            assert (M @ f == TT.gerstenhaber_bracket_cochain(E, 1, n, D, f)).all()
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement** (i) `_bracket_action_invariants` (read the degree-`(1,n)` table off the
  Plan-35 `HHProducts` and form the induced maps `HH^n → HH^n`, `D ↦ [D,·]`, on the GF(p) class
  basis) and `_action_invariants` (the same invariants off `ρ_n`): compare **basis-independent**
  data only (rank, char-poly, trace, weight multiset). (ii) the entry-wise helpers
  `_L_D_matrix_on_cochains` (our `L_D` per §1 built on the engine's degree-`n` cochain basis) +
  `_der_reps_on_engine_basis` / `_cochain_basis_vectors` for the STRONG check
  (`_L_D_matrix_on_cochains ≡ engine.tt_calculus.gerstenhaber_bracket_cochain(·,1,n,·)` entry-wise
  over GF(p)). Record the arbiter result in the report `note`. **The primary `ρ_n` is unchanged** —
  these are checks, not routes.
- [ ] **Step 4:** run; **Step 5:** commit `feat(hochschild): in-window Gerstenhaber-bracket sign arbiter (induced-map invariants) + STRONG entry-wise L_D == circle_cochain check for the HH-Lie-module action (Plan 71; deg-(1,n) == L_D, GF(p))`.

---

## Task 3 — char-0 weight / torus decomposition (the plan's HARDEST deliverable; NET-NEW code)

**This is the hardest task and NONE of its machinery comes from P70** (§2): P70 exposes only
integers/bools and *explicitly scopes out* the radical-inclusive torus P71 needs. The maximal torus
is computed from scratch here.

**Files:** modify `lie_module.py`; extend the test.

- [ ] **Step 1: Failing tests (assert the basis-independent PATTERN, never the literal integers —
  §2 NON-NORMATIVE; a scalar-doubled generator sends `{−2,0,2}` to `{−4,0,4}`).**
```python
lit = pytest.mark.oracle_literature
from quiverlab.errors import QuiverlabError

def _is_sl2_string(ws):                       # symmetric {-c,0,c}, c != 0, each mult 1
    s = sorted(ws); c = -s[0]
    return len(s) == 3 and c != 0 and s == [-c, 0, c]

def _is_equal_gap(ws):                        # arithmetic progression, common difference != 0
    s = sorted(set(ws)); diffs = {b - a for a, b in zip(s, s[1:])}
    return len(s) < 2 or (len(diffs) == 1 and 0 not in diffs)

@lit
def test_kronecker_hh1_is_sl2_string():
    A = Quiver([1,2], {"a":(1,2),"b":(1,2)}).algebra(field=QQ)     # HH^1 = adjoint sl2
    L = lie_module_action(A, top=1)
    assert _is_sl2_string(_weights_at(L, 1))                       # PATTERN, normalization-free
    assert set(_weights_at(L, 0)) == {0}                           # HH^0 trivial L(0), one weight 0

@lit
@pytest.mark.parametrize("n", [2, 3])
def test_truncpoly_equal_gap_weights_indecomposable(n):
    """k[x]/x^n: each HH^m has an equal-gap weight progression (gap is normalization-dependent, so
    assert equal-gap, NOT gap==1) and is a SINGLE indecomposable over the truncated-Witt HH^1."""
    A = truncated_polynomial(n, field=QQ)
    L = lie_module_action(A, top=2)
    for m in (1, 2):
        assert _is_equal_gap(_weights_at(L, m))                    # PATTERN, not literal gap
        assert len(L.summands[m]) == 1                             # single indecomposable

@selfcert
def test_weights_gated_char_p_but_summands_independent():
    """WEIGHTS are char-0-gated (loud); the SUMMAND decomposition is governed independently by
    decompose's own char guard (char 0 OR char > d_n) -- NOT by the weight gate (MAJOR-fix b)."""
    A = truncated_polynomial(3, field=GF(3))                       # W_1 in char 3
    L = lie_module_action(A, top=2)
    assert L.module_axiom_ok is True                               # field-general part still computed
    assert L.weights is None and L.char0_note                      # WEIGHTS gated + explained
    with pytest.raises(QuiverlabError):
        lie_module_action(A, top=2, require_char0=True)            # the weight/torus accessor raises
    # NB: L.summands is NOT asserted None -- it follows decompose's guard per degree (present where
    # char > d_n, a per-degree {"error": note} where char p <= d_n). See Task 4.

@selfcert
def test_weight_spaces_sum_to_dim():
    A = truncated_polynomial(3, field=QQ)
    L = lie_module_action(A, top=3)
    for n in range(4):
        assert sum(d for _, d in _weight_dims(L, n)) == L.hh_dims[n]
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement the from-scratch maximal-torus algorithm** (char-0 gate first: else
  `weights = None` + `char0_note`, and `require_char0=True` raises — the SUMMANDS are computed
  separately in Task 4, independent of this gate). The algorithm, each step with its own self-cert:
  1. **Cartan subalgebra of `HH¹`** — the de Graaf split-element / Fitting-null-space method on the
     `Der/Inn` bracket (a non-nilpotent element `x`, the Fitting-null component of `ad x`, iterate
     until nilpotent and self-normalizing). *Self-cert:* the result `h` is nilpotent-acting on
     itself and equals its own normalizer in `HH¹`.
  2. **ad-semisimple part** — within the Cartan, take the elements whose `ad`-action on `HH¹` is
     diagonalizable over `k` (the semisimple part of the abstract Jordan decomposition, exact over
     `QQ`). *Self-cert:* each extracted `h` has `ad h` with a squarefree minimal polynomial
     splitting over `k` (semisimple), and they pairwise commute.
  3. **Radical-toral extension** — extend that toral set to a MAXIMAL ad-diagonalizable abelian
     subalgebra of ALL of `HH¹` (greedily add ad-semisimple elements — including ones in the
     solvable radical, e.g. the `k[x]/x^n` grading `x∂` — that commute with the current torus and
     keep it ad-diagonalizable). This is the object P70 scoped out; record `torus_provenance`
     (basis-dependent — which subalgebra, and each `h`'s normalization). *Self-cert:* the final `t`
     is abelian, every `h ∈ t` is ad-semisimple, and no ad-semisimple element outside `t` commutes
     with all of `t` (maximality).
  4. **Simultaneous rational diagonalization** — the `{ρ_n(h) : h ∈ t}` commute (self-cert, from
     the module axiom); simultaneously diagonalize them over `k` and read the joint weights
     (`weights`, per-weight dims, `Σ = dim HH^n`). **Anisotropic three-valued fallback:** if some
     `ρ_n(h)` has an eigenvalue outside `k`, set `weight_base_change_note` (P70 precedent) and do
     NOT fabricate a split; if `t = 0` (nilpotent `HH¹`), `weights = None` + note.
  The literal weight labels are recorded as provenance only; the plan asserts the PATTERN (§2).
- [ ] **Step 4:** run; **Step 5:** commit `feat(hochschild): from-scratch maximal-torus weight decomposition of HH^n over HH^1 (Plan 71; Cartan + ad-semisimple + radical-toral extension + rational simultaneous diagonalization; char-0-gated, honest three-valued) -- NET-NEW, not in P70`.

---

## Task 4 — indecomposable Lie-module summands via `modules/decompose.py`

**Files:** add `modules/decompose.py::decompose_representation` (refactor: extract the matrix-level
Fitting core; `decompose(M)` delegates — byte-stable); modify `lie_module.py`; extend the tests.

- [ ] **Step 1: Failing tests.**
```python
from quiverlab.modules.decompose import decompose_representation

@selfcert
def test_decompose_representation_matches_module_decompose():
    """The refactor is byte-stable: decompose(M) still returns the same summands (it now
    delegates to decompose_representation on M's action matrices)."""
    A = truncated_polynomial(3, field=QQ)                          # a small quiver-module regression
    # (worker: assert decompose(some standard module) unchanged vs the pre-refactor golden)

@lit
def test_kronecker_hh1_single_irreducible():
    A = Quiver([1,2], {"a":(1,2),"b":(1,2)}).algebra(field=QQ)
    L = lie_module_action(A, top=1)
    s = L.summands[1]                                              # HH^1 summands
    assert len(s) == 1 and s[0]["dim"] == 3                        # dim End=1 => irreducible adjoint L(2)

@lit
@pytest.mark.parametrize("n", [2,3])
def test_truncpoly_hh_n_indecomposable(n):
    A = truncated_polynomial(3, field=QQ)
    L = lie_module_action(A, top=3)
    assert len(L.summands[n]) == 1                                 # single uniserial indecomposable

@lit
def test_rad2_loops_multisummand():
    """k[x,y]/(x,y)^2: HH^n (gl2-modules) split into >= 2 indecomposables (dim End = 2,2,6)."""
    A = Quiver([1], {"x":(1,1),"y":(1,1)}).algebra(relations=["x*x","x*y","y*x","y*y"], field=QQ)
    L = lie_module_action(A, top=2)
    assert len(L.summands[1]) >= 2 and len(L.summands[2]) >= 2
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement** `decompose_representation(gens, domain, budget)`: `End = ` commutant of
  `gens` (exact `fields.linalg.nullspace` of `X gᵢ − gᵢ X = 0`); drive the SAME Fitting min-poly
  split + local-endomorphism certificate (`dim End = 1 ⟹ indecomposable`, every char; else the
  Dickson/CIW trace-form under `char 0 or char > dim`, loud refusal otherwise) that
  `decompose(M)` uses; return certified summands `(subspace basis, restricted action)`. Refactor
  `decompose(M)` to delegate (module → action matrices → `decompose_representation`), **byte-stable**
  (pin a pre-refactor golden). In `lie_module_action`, call `decompose_representation([ρ_n(D)…],
  A.domain)` per degree **INDEPENDENTLY of the char-0 weight gate** — whenever decompose's own guard
  holds (`char 0 or char > d_n`), so summands appear even at `char p > d_n` where `weights` is
  `None` (MAJOR-fix b). `summands[n]` = dims + `label` when a summand is a recognizable irreducible
  (e.g. an sl₂ `L(k)` by its weight string — else `None`, never guessed); where decompose refuses
  (`char p ≤ d_n`, no `dim End = 1` cert) the per-degree entry is `{"error": note}` (the loud
  refusal surfaced as a note, not a crash, not a guess).
- [ ] **Step 4:** run; **Step 5:** commit `feat(modules,hochschild): decompose_representation (matrix-level Fitting core) + indecomposable HH-Lie-module summands (Plan 71; kK2 HH^1 irreducible, k[x,y]/(x,y)^2 multi-summand); decompose(M) delegates, byte-stable`.

---

## Task 5 — literature decomposition anchors (MNPRS / ALS) + GAP module cross-check (`-m qpa`)

**Files:** extend `tests/hochschild/test_lie_module.py` (lit) + `src/quiverlab/qpa/crosscheck.py`
(add `crosscheck_hh_lie_module` where feasible); test `tests/qpa/test_hh_lie_module_qpa.py`.

- [ ] **Step 1: Failing tests.**
  - **`kK₂` as a live ALS-toupie member** (self-contained, no `# PIN`): `kK₂` IS a toupie algebra
    (single source, single sink, parallel paths — §5), so `als_toupie_hh_lie` is a REAL citation
    for its `HH•` Lie-module structure (`HH^1 = ` irreducible `L(2)`), not just an analogy.
  - **`k[x]/(x^n)` as the self-contained truncated-Witt anchor** (no `# PIN`): `HH¹` = the
    positive-Witt truncation (an *analogue* of a Virasoro subquotient — NOT a citation of MNPRS's
    theorem, whose class is special-biserial), each `HH^n` a single indecomposable — the §5 pin.
  - **MNPRS Virasoro / larger ALS toupie table** (`# PIN`, deepest anchor): transcribe verbatim
    (from the paper's decomposition table, BEFORE the test) the degree-`n` indecomposable-summand
    dims for the SMALLEST tractable published member (a small self-injective special-biserial
    `kZ_m/J^ℓ` from MNPRS, or a larger ALS toupie), and assert `lie_module_action` reproduces them.
    The worker chooses the smallest that fits the budget and matches a published table — do NOT
    invent numbers.
```python
# tests/qpa/test_hh_lie_module_qpa.py
import pytest
pytestmark = pytest.mark.qpa
from quiverlab.qpa.session import gap_available

@pytest.mark.skipif(not gap_available(), reason="[qpa] libgap backend not installed")
def test_hh_lie_module_gap_multisummand_decomposition():
    """QPA has NO Hochschild Lie-module surface (honest skip that FAILS if it appears). The
    SUBSTANTIVE GAP oracle is the MULTI-summand case: k[x,y]/(x,y)^2, gl2 acting on HH^n with
    dim End = 2,2,6 -> >= 2 indecomposables. GAP's DirectSumDecomposition on the fed representation
    is a genuine independent count (the sl2 = kK2 case is only a trivial irreducible sentinel)."""
    from quiverlab.combinat.quiver import Quiver
    from quiverlab.fields import QQ
    from quiverlab.qpa.crosscheck import crosscheck_hh_lie_module
    A = Quiver([1], {"x":(1,1),"y":(1,1)}).algebra(relations=["x*x","x*y","y*x","y*y"], field=QQ)
    r = crosscheck_hh_lie_module(A, top=2)               # gl2 modules, multi-summand
    r.assert_agree()                                     # GAP summand dims == ours, per degree
    # trivial sentinel: sl2 = HH^1(kK2) acts irreducibly (1 summand) -- kept, but proves little
    r2 = crosscheck_hh_lie_module(Quiver([1,2],{"a":(1,2),"b":(1,2)}).algebra(field=QQ), top=1)
    r2.assert_agree()
```
- [ ] **Step 2:** confirm failure/skip.
- [ ] **Step 3: Implement** the lit anchors (transcribed tables) and `crosscheck_hh_lie_module(A,
  top)`: probe whether QPA exposes any Hochschild-Lie-module verb (`NamesGVars()` sweep; FAIL the
  skip if one appears — the Plan-35 precedent); feed the `HH¹` structure constants + the `ρ_n`
  matrices to GAP's Lie-module tools (`LieAlgebraByStructureConstants` + a matrix representation +
  `DirectSumDecomposition` / module irreducibility) and compare **per-degree summand dims** — the
  `gl₂` multi-summand case is the real cross-check, the sl₂ case the trivial sentinel. `# PIN`:
  confirm the exact GAP Lie-module function names against the installed GAP.
- [ ] **Step 4:** run `-m qpa`; **Step 5:** commit `test(qpa,hochschild): MNPRS/ALS decomposition anchors + GAP multi-summand (gl2) module cross-check for HH-Lie-module (Plan 71)`.

---

## Task 6 — the `hh_lie_module` HH kind (all three tiers, i18n ×4, decomposition-table block)

A `top`-carrying HH kind (the `gerstenhaber_brackets`/`cup_products` precedent — the budget caps
`A.dim` for the Der solve, `top` caps the bar degree). GUI touchpoints follow
`gerstenhaber_brackets`.

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (grammar parse — the Plan-35 `gerstenhaber_brackets`
  `top`-suffix branch: `hh_lie_module` or `hh_lie_module:3` — + `_dispatch` branch + `_snip`),
  `docs/gui/runner.py` (the Pyodide twin — matching grammar + dispatch + `_snip` + ETA,
  shape-identical), **`webapp/server/schema.py`** (the THIRD grammar-parse site — mirror the
  `gerstenhaber_brackets` branch), `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox
  `qlgui-hh-lie-module`, `S.ids`, push-list, `renderBlock` — the **decomposition table**,
  `scheduleProbe`), `webapp/templates/index.html` (one checkbox), `webapp/server/i18n/{en,es,fr,
  zh}.json` (all four locales), `src/quiverlab/trace/results_html.py` (`_HEADINGS` + a
  decomposition-table render branch), `tests/webapp/_runner_goldens.json` +
  `test_runner_delegation.py` (ONE golden).
- Test: `tests/webapp/test_hh_lie_module_kind_p71.py`, `tests/gui/test_hh_lie_module_twin_p71.py`.

**Block shape** (one shared `hh_lie_module_block(A, top, budget)`):
```python
{"kind": "hh_lie_module", "top": int, "hh_dims": [int,...], "hh1_dim": int,
 "module_axiom_ok": bool, "inner_acts_zero": bool, "characteristic": int,
 # WEIGHTS: char-0 ONLY (None over char p, char0_note explains):
 "weights": [{"n": int, "torus_rank": int, "weights": [[ [lam...], dim ], ...]}] | None,
 "torus_provenance": str|None,
 # SUMMANDS: INDEPENDENT of the weight gate -- present whenever decompose's guard holds
 # (char 0 OR char > d_n), a per-degree {"error": note} where char p <= d_n:
 "summands": [{"n": int, "parts": [{"dim": int, "label": str|None, "mult": int}] | {"error": str}}] | None,
 "weight_base_change_note": str|None, "char0_note": str|None,
 "status": str, "note": str|None,
 "references": ["als_toupie_hh_lie","mnprs_special_biserial","csss_gentle_tt","gerstenhaber1963","rss_hh1_lie"]}
# refusal (oversize A.dim / over-max_cells top) -> {"error": msg, "references":[...]}, never a 500.
```
The action **structure constants are NOT shipped in the block** (they explode and are
basis-dependent); the block reports `hh_dims`, the verdicts, the char-0 **weight table** and the
**decomposition table** (per degree `n`: `HH^n = ⊕ M_i^{mult}`, dims + labels).

- [ ] **Step 1: Failing cross-runner tests** (unmarked — extras-gated dir): `kK₂` block has
  `hh_dims == [1,3,0]`, `hh1_dim == 3`, a `summands` entry with one degree-1 part of `dim 3`, a
  degree-1 weight table that is a symmetric sl₂-string `{−c,0,c}` (the PATTERN — assert three
  weights symmetric about 0, each dim 1, NOT the literal `[-2,0,2]`; §2 NON-NORMATIVE),
  `"als_toupie_hh_lie"` in the citation keys; a `GF(3)` `k[x]/(x³)` block has `weights is None` + a
  `char0_note` (summands NOT asserted None — decompose's guard governs them independently); twin
  parity (`json.dumps(sort_keys=True)` across `hpc/spec.py` and `docs/gui/runner.py`).
- [ ] **Step 2: Implement** the grammar parse in **all three** sites (`top` suffix
  `hh_lie_module` or `hh_lie_module:3`, default `top`, budget `DEFAULT_MAXDIM = 48`; copy the
  `gerstenhaber_brackets` branch), the `_dispatch` branch
  (`from quiverlab.hochschild.lie_module import hh_lie_module_block`), the twin, the `gui.js`
  checkbox + decomposition-table `renderBlock` + `scheduleProbe` ETA (**estimator note**: the
  action matrices scale with `dim HH^n × dim HH^1`, and the bar cochain complex is exponential in
  `top` — so `A.dim` (`DEFAULT_MAXDIM`) bounds the Der solve but **`max_cells` is the binding
  limiter once `top ≥ 2`** (the `d·(d−1)^{top+1}` bar-cell wall `A.hochschild_cohomology` already
  enforces). Size `hh_lie_module` ABOVE `gerstenhaber_brackets` on `A.dim` AND `top`, reuse the
  Plan-35 sizing bucket so anything past the instant box routes OFF instant; the dim-220 Nakayama
  examples get an honest `status="budget"`, NO `hh_lie_module` — the Plan-35 omission precedent),
  the `_snip` recipe (`"hh_lie_module": "A.hh_lie_module(top=2)"`),
  `_HEADINGS["hh_lie_module"] = "HH• as a Lie module over HH¹"` + the decomposition-table render
  branch.
- [ ] **Step 3: Add ONE golden** (`hh_lie_module_kronecker`) to `_runner_goldens.json`; note it in
  `test_runner_delegation.py`'s change-log docstring; confirm existing goldens byte-identical
  first. **Canonical-key stability:** schema-v1 algebra-only, no `module` block → canonicalizes
  unchanged (confirm an existing family request's key is byte-stable).
- [ ] **Step 4:** run `tests/webapp/test_hh_lie_module_kind_p71.py tests/webapp/test_runner_
  delegation.py tests/gui/test_hh_lie_module_twin_p71.py tests/webapp/test_i18n.py tests/hpc -q`.
- [ ] **Step 5:** commit `feat(gui,webapp,hpc): hh_lie_module HH kind (Plan 71) -- HH^n weights + indecomposable-summand decomposition table, both runners byte-identical, grammar in 3 sites, i18n x4, one golden`.

**i18n keys — ALL FOUR LOCALES** (mirror `pick.kind.gerstenhaber_brackets`):

| key | en | es | fr | zh |
|---|---|---|---|---|
| `pick.kind.hh_lie_module` | HH• as a Lie module over HH¹ | HH• como módulo de Lie sobre HH¹ | HH• comme module de Lie sur HH¹ | HH• 作为 HH¹ 上的李模 |
| `inv.hh_lie_module` | HH• Lie-module structure (weights, indecomposable summands) | Estructura de módulo de Lie de HH• (pesos, sumandos indescomponibles) | Structure de module de Lie de HH• (poids, facteurs indécomposables) | HH• 李模结构（权、不可分解直和项） |
| `block.hh_lie_module.title` | HH• as a Lie module over HH¹ | HH• como módulo de Lie sobre HH¹ | HH• comme module de Lie sur HH¹ | HH• 作为 HH¹ 上的李模 |
| `block.hh_lie_module.weights` | Weights (maximal torus) | Pesos (toro maximal) | Poids (tore maximal) | 权（极大环面） |
| `block.hh_lie_module.decomp` | Indecomposable summands | Sumandos indescomponibles | Facteurs indécomposables | 不可分解直和项 |
| `block.hh_lie_module.charp` | Weights/decomposition need characteristic 0 | Pesos/descomposición requieren característica 0 | Poids/décomposition exigent la caractéristique 0 | 权/分解需要特征 0 |

---

## Task 7 — verification page, README, suite gate

**Files:** `docs/verification.md`; `README.md`; `docs/plans/2026-08-07-metaplan-v1.0.0.md`
(tick P71); existing release gates.

- [ ] **Step 1: Verification page.** Add the P71 subsystem rows for
  `hochschild/lie_module.py` + `modules/decompose.py::decompose_representation`:
  - `oracle_literature`: `kK₂ HH^1 = ` irreducible sl₂-adjoint `L(2)`, weights a symmetric
    sl₂-string `{−c,0,c}` (the PATTERN, `c ≠ 0` — NOT the literal integers, §2 NON-NORMATIVE)
    (`HH^0 = L(0)`); `k[x]/(x^n) HH^n` a single indecomposable over the truncated-Witt `HH¹` (a
    Virasoro-subquotient *analogue*); `k[x]/(x²)` weight ladder; `k[x,y]/(x,y)²` `HH^n`
    multi-summand `gl₂`-modules; `kK₂` a live ALS-toupie member + the MNPRS/ALS transcribed
    decomposition anchor.
  - `oracle_crossengine`: `hh_dims == A.hochschild_cohomology`; the **STRONG entry-wise** check
    (`L_D` on the engine cochain basis ≡ shipped `gerstenhaber_bracket_cochain(·,1,n,·)` over GF(p))
    AND the weaker in-window `gerstenhaber_brackets` degree-`(1,n)` induced-map arbiter (basis-
    independent invariants); normalized-bar ≡ unnormalized-probe structure.
  - `oracle_selfcert`: module axiom; `ρ_n(ad_x) = 0`; representative-independence; torus commuting;
    weight-dims sum to `dim HH^n`; summands reassemble `HH^n`; the char-p loud gate **on
    weights/torus** (summands governed independently); the oversize/`max_cells` refusal.
  - `qpa`: no QPA Lie-module surface (honest skip that FAILS if it appears); the substantive GAP
    oracle is the MULTI-summand `k[x,y]/(x,y)²` (`gl₂`) decomposition (`dim End = 2,2,6`); the sl₂
    = `kK₂` case is a trivial irreducible sentinel.
  - **Honest-scope entries (binding):**
    (a) **The weight/torus decomposition is char-0 ONLY** (P70's gate — no positive-char
    Cartan/Weyl; reference P70's rationale, do not re-argue). The **indecomposable-summand
    decomposition is NOT tied to that gate** — it follows `decompose`'s own char guard
    independently (entry (b)). The field-general `ρ_n` / module-axiom / graded-module surface stays
    live in every characteristic.
    (b) **The summand decomposition inherits `decompose.py`'s char guard** (`char 0 or char > d_n`;
    loud refusal at `char p ≤ d_n` with no `dim End = 1` certificate) — never a guessed
    decomposition.
    (c) **The weight torus is a DIFFERENT object from P70's `toral_rank` — NET-NEW code P70 scoped
    OUT** (not a reuse/extension of P70 scaffolding): a maximal ad-diagonalizable abelian subalgebra
    of ALL of `HH¹`, including radical toral elements like the `k[x]/x^n` grading (MNPRS's "Virasoro
    degree"). P70 exposes only integers/bools (no Cartan/torus matrices) and explicitly scopes this
    torus out, so P71 computes it from scratch (Task 3 — the plan's hardest deliverable). The exact
    weight labels are **basis-dependent / NON-NORMATIVE** (scalar-doubling a generator doubles every
    weight; the torus-normalization provenance is recorded — the P70 sl₂-triple precedent);
    normative = distinct-weight count, per-weight dims, weight PATTERN (sl₂-string / equal-gap).
    (d) **The action structure constants are basis-dependent** (`basis="der_inn/bar"`);
    cross-engine/cross-run comparison uses only basis-independent data — the Plan-35 rule.
    (e) **The sign arbiter is GF(p) in-window** (Plan 35 is GF(p); P51 extends to any Domain
    in-window; degree `(1,n)` is in-window). **The field-general `L_D` = the shipped `(1,n)`
    bracket by a source-level sign identity** (`p=1` collapses all Koszul signs) — the primary
    route needs no bracket engine (the whole point of R12).
    (f) **QPA has no Hochschild Lie-module surface** (the Plan-35 products precedent); the oracle is
    GAP's core Lie library where the acting algebra is GAP-native, plus the internal identity
    batteries.
    (g) **Cost/budget** — the action is the P70 Der solve (`≈ d^5.4` over QQ, bounded by
    `DEFAULT_MAXDIM = 48` on `A.dim`) PLUS the exponential bar cochain complex in `top`; **for
    `top ≥ 2` the binding limiter is `max_cells`, not `DEFAULT_MAXDIM`** (the `d·(d−1)^{top+1}`
    bar-cell wall `A.hochschild_cohomology` already enforces, reused). Every oracle (max `dim` 4,
    `top ≤ 4`) fits trivially; the dim-220 Nakayama webapp examples carry NO `hh_lie_module` (the
    Plan-35 omission precedent).
    (h) **The MNPRS / ALS decomposition tables at scale are `# PIN`'d, not invented** — the pinned
    small cases reproduce the *structure*; the exact published tables are transcribed verbatim at
    implementation for one small member.
  - **Recount the class table** (`tests/release/test_oracle_classes.py` drives the numbers — run
    collection, paste LIVE counts, re-run to green; no guessed at-authoring number).
- [ ] **Step 2: README.** One features line: "HH• as a graded Lie module over HH¹ — the
  Gerstenhaber degree-1 action (the field-general Lie derivative, over any exact field), its
  weight/torus decomposition over characteristic 0 and indecomposable Lie-module summands; the
  Kronecker `HH¹(kK₂) ≅ sl₂` acting irreducibly on `HH^1`, the `k[x]/(x^n)` truncated-Witt grading
  (a Virasoro-subquotient analogue) — R12."
- [ ] **Step 3: Full gate:** `... tests/hochschild -q -m fast`, `... tests/webapp tests/gui
  tests/hpc -q -m fast`, `... tests/qpa -q -m qpa`, `... tests/modules -q -m deep` (the
  `decompose` refactor regression), `... tests/release tests/citations -q`, plus a
  `QUIVERLAB_NO_NUMBA=1` parity spot-run of `tests/hochschild/test_lie_module.py` — all green.
- [ ] **Step 4:** commit `docs(verification): P71 HH-Lie-module oracle rows + honest scope (char-0 weights/decomposition gate, broader torus, sign arbiter, MNPRS/ALS anchors) + recounted classes; README line; tick P71`.

---

## Acceptance (Plan-71 definition of done)

1. **`quiverlab.hochschild.lie_module` public:** `lie_derivative_on_hh`, `lie_module_action`,
   `hh_lie_module_block`; `Algebra.hh_lie_module` delegate; `modules/decompose.py::
   decompose_representation` (with `decompose(M)` delegating, byte-stable). Computes over **any
   exact field** on the normalized bar cochain complex + P70's `Der/Inn` (no resolution engine).
2. **Field-general action green (any characteristic + field kind):** `hh_dims ==
   A.hochschild_cohomology` on the zoo; the module axiom `ρ_n([D,E]) = [ρ_n(D),ρ_n(E)]` and
   `ρ_n(ad_x) = 0` (factors through `Der/Inn`) self-certified on `k[x]/(x^n)`, `k[x,y]/(x,y)²`,
   `kK₂`; the **STRONG entry-wise** `L_D ≡ gerstenhaber_bracket_cochain(·,1,n,·)` check over GF(p)
   AND the in-window `gerstenhaber_brackets` degree-`(1,n)` induced-map arbiter agree.
3. **char-0 weights green + gated (WEIGHTS only):** the from-scratch maximal-torus algorithm gives
   `kK₂ HH^1` weights a symmetric sl₂-string `{−c,0,c}` (the PATTERN — NOT the literal integers,
   §2 NON-NORMATIVE; `HH^0 = {0}`); `k[x]/(x^n) HH^n` an equal-gap weight progression; `k[x]/(x²)`
   weight ladder; weight-dims sum to `dim HH^n`; over `GF(p)` the **weight/torus** block is a
   **loud `QuiverlabError`** (`require_char0=True`) or `None` + a `char0_note`; anisotropic torus →
   `weight_base_change_note`, no fabricated split.
4. **Indecomposable summands green (INDEPENDENT of the weight gate):** `decompose_representation`
   reproduces `decompose(M)` on the module regression (byte-stable); `kK₂ HH^1` = **one
   irreducible** (`dim 3`, adjoint `L(2)`, `dim End = 1`); `k[x]/(x^n) HH^n` = one indecomposable;
   `k[x,y]/(x,y)²` `HH^n` multi-summand (`dim End = 2,2,6`); summands follow `decompose`'s own guard
   (`char 0 or char > d_n`) so they appear even at `char p > d_n` while `weights` is `None`; the
   char-`p ≤ d_n` loud refusal is surfaced as a per-degree note, honored.
5. **Literature + QPA oracles:** `kK₂` pinned as a live ALS-toupie member + the `k[x]/(x^n)`
   truncated-Witt anchor (self-contained); the MNPRS special-biserial / larger-ALS-toupie
   decomposition table transcribed for one small member (`# PIN`); GAP crosschecks the **multi-
   summand** `gl₂` decomposition (the substantive oracle) with the sl₂ case as a trivial sentinel
   (QPA has no such surface — honest skip).
6. **One HH kind `hh_lie_module`** clickable end-to-end (GUI canvas → block → **decomposition
   table**) in **all four locales** with `pick.kind.*` keys and exact key parity, schema v1, both
   runners byte-identical, grammar parsed in **all three sites**, ONE golden with a documented
   change-log entry, canonical keys byte-stable (no `module` block), the estimator sizing on
   `A.dim × top`; the weights + summand tables render.
7. `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the eight
   honest-scope entries (a)–(h); citations added and BibTeX-verified (Artenstein–Lanzilotta–Solotar
   authors CORRECTED from the card's guess; MNPRS `J. Algebra 580 (2021) 264–298` `# PIN`'d for a
   landing-page re-confirm; the Virasoro-subquotient claim verbatim-verified); README line; metaplan
   P71 ticked; fast + qpa +
   release + citations + the `tests/modules` decompose regression green; a `QUIVERLAB_NO_NUMBA=1`
   parity spot-run of the new suite green. **Builds on P70** (`derivations`/`inner_derivations` +
   the char-0 classification) with no respec of P70's contracts; **feeds P78** (the L∞ /
   Maurer–Cartan deformation plan consumes this Lie-module surface).

---

## Methodology & assumptions

**Approach.** I (1) read the metaplan §5 P71 card + §§1–4 conventions and the R12 record
verbatim; (2) read the **committed P70 doc** (`docs/plans/2026-08-08-plan-70-hh1-lie.md` on `dev`;
`d5a26c7` at authoring, but `dev` advances — I refer to the doc, not a tip hash) end-to-end — its
`Der/Inn` surface (`derivations`/`inner_derivations`), char-0 classification, `DEFAULT_MAXDIM = 48`,
and arithmetic-field base-change gating are the contracts P71 builds on — plus plan-65/plan-63 as
house-style templates; (3) discovered my worktree was a STALE checkout (`b7fa566`, v0.3.0) behind
`dev` — so I read those inputs from `dev` (`git show dev:…`) and ran the authoring-round live
verification against the importable `quiverlab` (which ships the Plan-35 bracket surface and the
`hochschild/bar.py` field-general substrate P71 needs; the fix round re-verified against current
`dev`, P67 merged — see the Fix-round note below); (4) **read the shipped Gerstenhaber
bracket source** (`engine/tt_calculus.py::circle_cochain` / `gerstenhaber_bracket_cochain`) and
PROVED the degree-`(1,n)` bracket equals the plan's field-general Lie derivative sign-for-sign
(`p = 1` collapses every Koszul sign) — discharging the card's #1 sign risk at the source level;
(5) confirmed the build-on surfaces live (`derivations` semantics from P70; `bar.coboundary_matrix`
+ `_cochain_basis` over any Domain; `A.hochschild_cohomology`; `A.gerstenhaber_brackets(top=2)`
returning `(1,1)/(1,2)` tables; `modules/decompose.py` is quiver-coupled via `hom_space`); (6)
**live-verified every numeric/structural pin two independent ways** — a standalone exact
Hochschild-cochain + Lie-derivative + envelope/commutant probe over `QQ` (sympy), AND
`A.hochschild_cohomology` for the `HH^n` dims (agreeing on `kK₂ [1,3,0]`, `k[x]/x² [2,1,1,1,1]`,
`k[x]/x³ [3,2,2,2,2]`, `k[x,y]/(x,y)² [3,4,6,12]`); (7) fetched all three R12 arXiv abstracts,
CORRECTED the card's 1803.10310 authorship (Artenstein–Lanzilotta–Solotar, "toupie algebras"),
verbatim-verified the MNPRS Virasoro-subquotient claim, and matched the MNPRS journal ref (J.
Algebra 580 (2021) 264–298) to a ScienceDirect search hit (`# PIN`'d for a landing-page re-confirm,
per the fix round).

**Fix round (2026-08-08, adversarial review NEEDS WORK → all findings VALID; document-only, dev
re-verified).** The critic confirmed the mathematical core unusually thoroughly (the `p=1` sign
collapse verified ENTRY-WISE against the shipped `circle_cochain` on `k[x]/x³/GF(7)`, every pin
reproduced on current `dev`, the envelope/commutant design validated, both citation corrections
confirmed). I re-verified against **current `dev`** (the plain venv `import quiverlab` resolves to
the main-repo `dev` src — `kK₂ [1,3,0,0]`, `k[x]/x³ [3,2,2,2]` reproduced with P67 merged) and
confirmed `engine.tt_calculus.circle_cochain` / `gerstenhaber_bracket_cochain` exist for the
entry-wise self-cert. Edits made: **(MAJOR a)** Task 3 rewritten as a from-scratch maximal-torus
algorithm (Cartan → ad-semisimple → radical-toral extension → rational simultaneous diagonalization,
each with its own self-cert), all "reuse/extends P70 scaffolding" language dropped and replaced by
"NET-NEW code P70 scoped out — the plan's hardest deliverable"; **(MAJOR b)** the weight/torus gate
(char-0) and the summand gate (`decompose`'s `char 0 or char > d_n`) made consistent and
INDEPENDENT everywhere (§2(iii), scope table, API, block shape, Task 3/4, honest-scope (a), API
`summands` doc, Acceptance 3/4); **(MAJOR c)** every weight test now asserts the basis-independent
PATTERN (`_is_sl2_string`, `_is_equal_gap`), never the literal integers (a scalar-doubled torus
gives `{−4,0,4}`), and the STRONG entry-wise `L_D ≡ circle_cochain` cross-engine check was added to
Task 2 + the marker list; **(MINORS)** dev-tip references made tip-agnostic, the MNPRS vol/pages
`# PIN`'d, the GAP oracle re-pointed at the substantive multi-summand `gl₂` case (sl₂ demoted to a
sentinel), the `k[x]/x^n` "Virasoro-subquotient" softened to an analogue + `kK₂` added as a real
ALS-toupie anchor, the `DEFAULT_MAXDIM`-vs-`max_cells` framing corrected (`max_cells` binds for
`top ≥ 2`), and the Task-1 test-bucket note fixed (`tests/hochschild/` is FAST by conftest; heavy
probes go to a deep-bucketed dir).

**Assumptions.** (1) My worktree's `quiverlab 0.3.0` engine is a faithful verification target for
the field-general mathematics (the Lie derivative, HH cochain complex, decompose char guards) —
these surfaces predate P51 and are stable; P70's `derivations` is a *planned* surface (unimplemented
anywhere), so I hand-rolled it (the Leibniz null space) exactly as P70 specs, matching P70's own
live-verified `HH¹` dims (`kK₂ = 3`, `k[x]/x³ = 2`). (2) The degree-`(1,n)` Gerstenhaber bracket
equals the commutator/Lie-derivative action (Gerstenhaber 1963) — I verified the *sign identity at
the source level* and the *module axiom + inner-zero + dim agreement* computationally, and rely on
the classical theorem for the cohomology-level identification, with the shipped `(1,n)` bracket as
the in-window arbiter. (3) The associative envelope `B_n` / commutant is the correct decompose
input — I verified `End_{B_n} = ` commutant `= End_{HH¹}` and the certificate outcomes (kK₂ `dim
End = 1`; `k[x,y]/(x,y)²` `dim End = 2,2,6`). (4) GAP's core Lie-module tools are reachable via the
`[qpa]` libgap backend (standard GAP) — I confirmed P70's access-path finding but did NOT run GAP.
(5) Importing `quiverlab` in a plain computation is side-effect-free w.r.t. the webapp job queue
(the data-dir-collision note is about the desktop worker, not library imports) — though I did note
and clean a stray `quiverlab_traces/` dir the `hochschild_cohomology` call wrote to CWD.

**Deliberately NOT checked, and why.** (1) I did not run any pytest suite (other agents run
`-m qpa`/`-m deep` concurrently; I kept probes short, single-process). (2) I did not run GAP (same
reason + heaviness); the GAP Lie-module function names are `# PIN`'d for the worker. (3) I did not
compute the general maximal-torus (Cartan-subalgebra) algorithm — I used explicit, hand-identified
torus elements (`x∂`, the sl₂-Cartan, the `gl₂` diagonal), which is what the weight pins assert;
the general de Graaf Cartan computation is the worker's Task-3 work with the §2 three-valued
fallbacks. (4) I did not transcribe the MNPRS/ALS decomposition tables (their algebras need a
per-paper read; the tables are `# PIN`'d for the worker, and the `k[x]/(x^n)` Virasoro-subquotient
instance is the self-contained small anchor I DID verify). (5) I did not implement
`decompose_representation` — I verified the equivalence (`B_n`/commutant/`dim End`) that decides its
certificate. (6) I did not independently re-confirm the CSSS J. Algebra 708 (2026) 138–231
volume/pages (authors/title/abstract verified; volume/pages taken from R12's "independently
confirmed" note, `# PIN`'d). (7) I did not mutate my worktree's git state — only the plan doc is
written, per the task.

**Why I believe the result is correct.** Every numeric pin is verified two independent ways and
they agree everywhere: the standalone Hochschild-cochain + Lie-derivative probe reproduces
`A.hochschild_cohomology`'s `HH^n` dims on all four algebras; the module axiom `ρ_n([D,E]) =
[ρ_n(D),ρ_n(E)]` and `ρ_n(ad_x) = 0` hold at every computed degree (an independent confirmation
that the Lie-derivative implementation IS a Lie-module action factoring through `HH¹ = Der/Inn`);
the sl₂-string weights on `kK₂ HH^1` (a symmetric `{−c,0,c}` — the PATTERN, not literal integers),
the truncated-Witt equal-gap weights on `k[x]/(x^n)`, and the `gl₂`-weight pairs on
`k[x,y]/(x,y)²` all reproduce the expected representation theory; and the envelope/commutant data
(`kK₂` `B_1 = M₃`, `dim End = 1`; `k[x,y]/(x,y)²` `dim End = 2,2,6`) proves decompose's certificate
fires exactly as the design predicts. The single highest-risk piece — the sign convention — is
discharged by a source-level identity (`p = 1` kills every Koszul sign), reinforced by the strong
entry-wise + in-window arbiter oracles. Every scope boundary R12 names (field-general action
WITHOUT the bracket engine; char-0 for the weight/torus decomposition; `decompose`'s own char guard
for the summands, independently)
becomes a loud typed refusal or an honest provenance note, and every design deviation from the
literal card (the corrected 1803.10310 authorship, the torus broader than P70's `toral_rank`, the
`decompose_representation` refactor rather than a quiver wrapper) is stated with its evidence.

## Open design risks (what a critic will likely attack)

1. **The general maximal-torus / weight computation.** The least-verified piece: I used explicit
   torus elements per example, not a general Cartan-subalgebra algorithm. *Mitigation:* the de
   Graaf Cartan computation over char 0 (Task 3) with the §2 three-valued fallbacks (trivial torus /
   anisotropic → base-change note / char-p loud gate); the exact weight labels are NON-NORMATIVE
   (the P70 sl₂-triple precedent), only the pattern/multiset is pinned.
2. **The `decompose(M)` refactor byte-stability.** Extracting the matrix-level core and delegating
   risks perturbing the existing (Plan-30) decompositions. *Mitigation:* a pre-refactor golden
   pins `decompose(M)` on the standard module zoo; the refactor is purely a re-plumb (the Fitting
   math is unchanged), and `decompose_representation` reuses the identical char guard.
3. **The sign convention.** The card's #1 risk. *Mitigation:* discharged by the source-level
   identity (`circle_cochain` with `p = 1` has all-`+1` signs, so `[D,g] = D∘g − g∘D` = the §1
   `L_D`) AND the in-window GF(p) arbiter oracle (Task 2), which compares basis-independent
   induced-map invariants against the shipped `gerstenhaber_brackets`.
4. **Cost.** The Der solve (`≈ d^5.4`, P70) PLUS the exponential bar cochain complex in `top`.
   *Mitigation:* `DEFAULT_MAXDIM = 48` + the reused `max_cells` bar-cell wall + honest
   `status="budget"`; the estimator sizes on `A.dim × top`; the dim-220 Nakayama examples carry no
   `hh_lie_module` (the Plan-35 omission precedent).
5. **MNPRS/ALS literature scope.** The Virasoro-subquotient and toupie decomposition tables are the
   deep literature anchors, but their algebras may exceed a tractable size or need careful
   transcription. *Mitigation:* the `k[x]/(x^n)` truncated-Witt (Virasoro-subquotient) instance is
   the self-contained small anchor I verified; the published tables are `# PIN`'d to transcribe
   verbatim for one small member, never invented.
6. **The normalized-bar reduction mod `k·1`.** On `Hom(Ā^{⊗n}, A)` the Lie derivative feeds
   `D(a_i)` (possibly with a unit component) into a bar slot living in `Ā`. *Mitigation:* reduce
   inserted values mod `k·1` — exactly the shipped `circle_cochain` normalization — so the
   normalized route and the in-window bracket agree; my probe (unnormalized) cross-validates the
   dims and structure.

## Change log

- **2026-08-08 authoring.** Initial plan (R12 HH• as a graded Lie module over HH¹ via the
  field-general Lie derivative `L_D = D∘f − Σ f(…,Da_i,…)`, char-0 weight/torus decomposition, and
  indecomposable summands via a `decompose_representation` refactor of `modules/decompose.py`).
  Builds on the committed P70 `Der/Inn` surface. **All three R12 references re-verified against
  their arXiv abstracts**; the metaplan card's 1803.10310 authorship (guessed
  "Assem–Lanzilotta–Schroll(?)") CORRECTED to **Artenstein–Lanzilotta–Solotar** (toupie algebras);
  the MNPRS Virasoro-subquotient claim verbatim-verified and its journal ref (J. Algebra 580 (2021)
  264–298) confirmed via ScienceDirect. **The #1 sign risk discharged at the source level**: the
  shipped `circle_cochain`/`gerstenhaber_bracket_cochain` in degree `(1,n)` (`p = 1`) has all-`+1`
  Koszul signs, so `[D,g] = D∘g − g∘D` is byte-identical to `L_D`. **Every numeric/structural pin
  live-verified two ways** (a standalone exact Hochschild-cochain + Lie-derivative + envelope
  probe, and `A.hochschild_cohomology`): `kK₂ HH• = [1,3,0]` with `HH^1 = ` irreducible sl₂-adjoint
  `L(2)` (weights `{−2,0,2}`, `dim End = 1`); `k[x]/(x^n)` HH^n indecomposable over the
  truncated-Witt (Virasoro-subquotient) `HH¹`; `k[x,y]/(x,y)² HH• = [3,4,6,12]` multi-summand
  `gl₂`-modules (`dim End = 2,2,6`); the module axiom and inner-acts-zero self-cert at every degree.
  **Three design decisions recorded with evidence:** (i) `HH^n` becomes a module over the
  associative envelope `B_n = ⟨ρ_n(HH¹),I⟩` (its Krull–Schmidt = the Lie-module decomposition,
  since `End_{B_n} = ` commutant `= End_{HH¹}`); (ii) `decompose(M)` is quiver-coupled, so P71 adds
  a matrix-level `decompose_representation` driven by the commutant rather than wrapping `B_n` as a
  quiver algebra; (iii) the weight torus is a maximal ad-semisimple abelian subalgebra of ALL of
  `HH¹` (including radical toral elements like the `k[x]/x^n` grading), broader than P70's
  Levi-only `toral_rank` — a stated extension, with basis-dependent torus provenance.
- **2026-08-08 fix round (adversarial review NEEDS WORK → all findings adjudicated VALID;
  document-only, current-`dev` re-verified with P67 merged).** The critic confirmed the math core
  entry-wise (the `p=1` sign collapse vs the shipped `circle_cochain` on `k[x]/x³/GF(7)`), reproduced
  every pin on `dev`, validated the envelope/commutant design, and confirmed both citation
  corrections. Edits: **(MAJOR a)** Task 3 rewritten as a from-scratch maximal-torus algorithm
  (Cartan → ad-semisimple → radical-toral extension → rational simultaneous diagonalization, each
  self-certified); all "reuse/extends P70 scaffolding" language replaced by "NET-NEW code P70 scoped
  out — the plan's hardest deliverable" (P70's `HH1LieStructure` exposes only ints/bools, no torus
  matrices, and *explicitly* scopes out the radical-inclusive torus). **(MAJOR b)** the weight/torus
  char-0 gate and the summand `char 0 or char > d_n` guard made consistent and INDEPENDENT
  everywhere (no unreachable gate rows): weights are char-0-only, summands follow `decompose`'s own
  guard and can appear at `char p > d_n` while `weights` is `None`. **(MAJOR c)** every weight test
  asserts the basis-independent PATTERN (symmetric sl₂-string `{−c,0,c}`; equal-gap progression),
  never the literal integers (a scalar-doubled torus generator gives `{−4,0,4}`); the STRONG
  entry-wise `L_D ≡ gerstenhaber_bracket_cochain(·,1,n,·)` cross-engine check (the critic ran it,
  passes) added to Task 2 + the marker list. **(MINORS)** tip-agnostic `dev` references; MNPRS
  vol/pages `# PIN`'d (ScienceDirect search hit, not landing-page-confirmed); the GAP oracle
  re-pointed at the substantive multi-summand `gl₂` case (sl₂ demoted to a sentinel); `k[x]/x^n`
  "Virasoro-subquotient" softened to an *analogue* and `kK₂` added as a REAL ALS-toupie anchor; the
  `DEFAULT_MAXDIM`-vs-`max_cells` framing corrected (`max_cells` binds for `top ≥ 2`); the Task-1
  test-bucket note fixed (`tests/hochschild/` is FAST by `conftest`; heavy probes → a deep-bucketed
  dir).
