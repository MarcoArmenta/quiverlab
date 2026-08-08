# Plan 67: Silting — Verifier, Single Mutation, Bounded Exploration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The silting-theory surface on top of P43's shipped derived stack: a
**silting-object VERIFIER** in `K^b(proj A)` (presilting `Hom_{D^b}(T, T[n]) = 0`
for **`n > 0` only** — the silting vanishing, weaker than tilting's `n != 0` —
plus generation `thick(T) = K^b(proj A)`, certified honestly and **three-valued**
where K₀ alone cannot decide), a **single silting MUTATION** `μ_X^±(T)` via **one
approximation triangle** (Aihara–Iyama Def 2.30/2.34, built from a complex-level
minimal `add(T/X)`-approximation + `ChainMap.cone`), a **bounded-radius
exploration** with **loud truncation** (NO general BFS — the silting quiver can be
infinite and mutation-transitivity is proven only for special classes), and the
**co-t-structure dictionary** (silting `T` ↦ the bounded co-t-structure with
coheart `add(T)`, AI Prop 2.23(b) / Jørgensen). The 2-term slice is a **cross-check
against P45's `τ`-tilting engine, not a deliverable**. This is the metaplan's P67
card (record R30, tier γ, size L).

**Architecture:** One new module `src/quiverlab/derived/silting.py` — a thin exact
layer over primitives that already ship: P43's `is_tilting_complex` / `TiltingReport`
/ `hyper_hom_basis` / `end_algebra_of_complex` / the private
`_direct_sum_complex` / `_span` / `_cochain_vec` (`derived/tilting.py`), **plus a new
shared `g_proj` K₀-basis helper introduced in Task 0** (the g-matrix must live in
`K₀(K^b proj A) = ⊕_v Z[P_v]`, the **projective** basis — NOT `_chi`'s
composition-factor basis; see the Task-0 bug fix and the sufficiency ruling),
P39's `ChainComplex` / `ChainMap.cone` / `ChainMap.triangle` / `ChainMap.then` /
`hyper_hom_dims` / `is_perfect` / `from_projective_resolution` (`modules/complexes.py`),
P44's retract-pruning **algorithm** (`modules/approximations.py`, lifted here from
modules to complexes), and P45's `two_term_silting` / `exchange_graph`
(`tautilting/`) for the cross-check. **No new math engine.** Every constructed
object self-certifies (the mutant is re-verified silting; the triangle's `d²=0` is
`ChainComplex(check=True)`; the neighbour shares `n-1` summands) and is pinned
against Aihara–Iyama's own worked examples and Oppermann's quiver rule —
correctness never rests on trusting a construction. Plus a `silting_block` in
`derived/block.py` and the `silting` compute kind (both runners + i18n ×4).

**Tech Stack:** `modules/complexes` (`ChainComplex`, `ChainMap`, `hyper_hom_dims`,
`identity_chain_map`, `projective_model`), `derived/tilting`
(`is_tilting_complex`, `TiltingReport`, `_direct_sum_complex`, `_span`, `g_proj`
(NEW, Task 0 — K₀ in the projective basis), `_cochain_vec`, `end_algebra_of_complex`),
`derived/homs` (`hyper_hom_basis`),
`modules/linalg_mod` + `fields.linalg` (`solve`, `reduce_mod_nullspace`, `rref`,
`mat_rank`, `cols_to_matrix`), `sympy` (integer `det` of the **g_proj** matrix, i.e.
the K₀ classes in the projective basis). No floats in `src/`.

## Record (verbatim — R30, from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R30 — Silting: verifier + single mutation + bounded exploration.** [D-scout P4,
> MAJOR re-scope per critic] Object: silting-object VERIFIER in K^b(proj)
> (Hom(T,T[n≠0]) = 0 + generation — on the shipped derived stack), silting MUTATION
> (one approximation triangle), and the co-t-structure dictionary. HONEST SCOPE: the
> silting quiver can be infinite and mutation-transitivity is known only for special
> classes (Aihara–Iyama's own statement) — no general BFS enumeration; 2-term slice =
> existing τ-tilting (cross-check only, not a deliverable); bounded-radius exploration
> with loud truncation; certified-complete only where finiteness is proven
> (derived-discrete, local, τ-tilting-finite 2-term). CORRECTED refs: Aihara–Iyama
> arXiv:1009.3370 (J. LMS 85 (2012)); Oppermann arXiv:1504.02617 'Quivers for silting
> mutation'; Jørgensen arXiv:1603.09379 'Co-t-structures: the first decade'. Oracles:
> k[x]/(x²) silting = shifts; kA₂ silting quiver (Aihara–Iyama); End(μT) quiver vs
> Oppermann's rule. Size L.

**Record correction folded back (dated 2026-08-07, per metaplan §2):** the record's
one-line summary writes the verifier vanishing condition as `Hom(T,T[n≠0]) = 0`.
That is the **tilting** condition. The **silting** condition — the object of this
plan — is `Hom(T,T[n]) = 0` for **`n > 0` only** (presilting), plus generation
(AI Def 2.1(a); tilting is Def 2.1(b) with `n != 0`). The verifier below scans the
**positive** window; the `n < 0` direction is what separates tilting ⊂ silting.
The research doc gets this right in its body; the P67 card's shorthand is corrected
here.

## Reference re-verification (mandatory, metaplan §5) — findings

Verified against the primary sources (arXiv PDFs read directly this session):

- **Aihara–Iyama, "Silting mutation in triangulated categories", arXiv:1009.3370,
  J. London Math. Soc. (2) 85 (2012), no. 3, 633–668.** CONFIRMED verbatim:
  - **Def 2.1.** `M ⊆ T` is **silting** if `Hom_T(M, M[>0]) = 0` and `T = thick M`;
    **tilting** if `Hom_T(M, M[≠0]) = 0` and `T = thick M`. (So silting checks only
    positive shifts; this is the plan's verifier.)
  - **Def 1.4.** `(X, Y)` is a **co-t-structure** if `(X[-1], Y)` is a torsion pair
    and `X ⊆ X[1]`; the **coheart** is `X ∩ Y`.
  - **Example 2.2(a).** `A` (as a stalk complex) is a tilting object of `K^b(proj A)`;
    every tilting complex is a tilting object. **Example 2.5(a).** `D^b(mod A)` has
    silting subcategories iff `gl.dim A < ∞` (but `K^b(proj A)` always has one: `A`).
  - **Thm 2.26.** If `T` has an **indecomposable** silting object `M`, then
    `silt T = {M[i] : i ∈ Z}`. (→ the `k[x]/(x²)` = shifts oracle, `n = 1`.)
  - **Thm 2.27 / Cor 2.28.** For a silting subcategory `M`, `ind M` is a **Z-basis of
    `K₀(T)`**; any two silting objects have the same number of indecomposable summands
    (`= n = rk K₀`). This is the **necessary** K₀ condition (det ±1); see the
    sufficiency ruling below.
  - **Def 2.30 (mutation).** For `M ∈ silt T` and a covariantly-finite `D ⊆ M`, take
    a left `D`-approximation `f: M → D` and the triangle
    `M -f-> D -g-> N_M -> M[1]` (triangle (4)); put
    `μ^+(M;D) := add(D ∪ {N_M | M ∈ M})` — the **left mutation**. Dually `μ^-(M;D)`
    (right mutation) for contravariantly-finite `D`.
  - **Def 2.34 (irreducible).** `μ_X^±(M) := μ^±(M; M_X)` where `ind M_X = ind M \ {X}`
    — mutation at a **single indecomposable summand** `X` (the "one approximation
    triangle"). AI's `μ^+` is LEFT; `M ≥ μ^+(M;D)` (Prop 2.33 — left mutation goes
    DOWN the silting order).
  - **Thm 2.31.** Any mutation of a silting subcategory is again silting. (→ the
    mutant is silting **by construction**; the plan re-certifies it as a self-cert,
    not to establish it.)
  - **Thm 2.35 + Def 2.41.** Under condition (F) (`T` Krull–Schmidt Hom-finite with a
    silting object — satisfied by `K^b(proj A)`, Prop 2.20), `N` is an irreducible
    left mutation of `M` iff `M > N` with no `L` strictly between; the **silting
    quiver** (vertices `silt T`, arrows = irreducible left mutations) is the Hasse
    quiver of the poset.
  - **Cor 2.43.** If `T` has an indecomposable silting object, iterated irreducible
    mutation is transitive (→ **local** algebras). **Thm 1.2.** (T) transitivity holds
    for **local, hereditary, or canonical** algebras; there is a symmetric algebra
    where it fails ([AGI]) — so **no general transitivity/BFS**.
  - **Prop 2.23(b).** For `M ∈ silt T`, `(⊥(T_M^{≤0}), T_M^{≤0})` is a co-t-structure
    with **coheart `M`** — the co-t-structure dictionary payload.
  - **Example 2.45 (kA₂, `1 → 2`).** The AR-quiver of `K^b(proj A)` is the `ZA₂`
    zigzag `… X_{-1}, X_0, X_1, X_2 …`; the silting objects are the pairs `X_a ⊕ X_b`,
    and the **silting quiver is drawn explicitly** (mod shift: three infinite `[1]`-
    chains `X_0⊕X_1 → X_0⊕X_4 → X_0⊕X_7 → …`, `X_1⊕X_2 → X_1⊕X_5 → …`,
    `X_2⊕X_3 → X_2⊕X_6 → …` with cross arrows). **kA₂ silting is INFINITE** → the
    exploration truncates loudly (honest), transitivity is by Thm 1.2 (hereditary).
  - **Example 2.47 (`1 ⇄ 2` with `ab = 0 = ba`).** Concrete mutation: from
    `P₁ ⊕ P₂`, the two irreducible left mutations are `X ⊕ P₂` and `Y ⊕ P₁` with
    `X := cone(P₁ → P₂)`, `Y := cone(P₂ → P₁)` — a **direct cone-computation oracle**
    for the mutation engine.
- **Oppermann, "Quivers for silting mutation", arXiv:1504.02617, Adv. Math. 307
  (2017) 684–714.** CONFIRMED. **Theorem 1.1:** for a dg quiver algebra `(kQ, d)`
  concentrated in non-negative degree and a vertex `i` with no degree-0 loops, the
  derived endomorphism ring of the **left** AI silting mutation of `kQ` at `i` is the
  **mutated quiver `(kM, ∂)`**, given by three steps: **(1) Rotation** — raise degrees
  of arrows *into* `i` by 1, lower degrees of arrows *out of* `i` by 1, and reverse
  each outgoing degree-0 arrow `α ↦ α*` (degree 0); **(2) Composition arrows** `αφ`
  for `φ` ending at `i`, `α` outgoing degree-0 from `i`, degree `‖φ‖ − 1`; **(3)
  Anti-composition arrows** `φα⁻¹`. Dual for right mutation. Hypotheses: `Λ ≅ T_R(A)`
  non-negatively graded, degree-wise finite over `R = k^n`, mutation vertex has no
  degree-0 self-loop (covers hereditary; general f.d. via Construction 2.2).
  Example 2.3 works `k[1→2→3]/(ψφ)`. **# PIN** the full graded dg rule (degree
  bookkeeping is out of the engine's scope — see Task 4); the plan verifies the
  **underlying degree-0 quiver** of `End(μT)` on worked instances.
- **Jørgensen, "Co-t-structures: the first decade", arXiv:1603.09379.** CONFIRMED a
  survey; co-t-structure ≡ weight structure (Bondarko), the bijection
  {bounded co-t-structures} ↔ {silting subcategories} via **coheart = add(silting)**
  (Bondarko; Mendoza–Sáenz–Santiago–Souto Salorio; Keller–Nicolás). Used as the
  documentation reference for the co-t-structure dictionary. **BibTeX-verify the
  publication venue at merge** (Abel Symposia proceedings, likely vol. 17, Springer
  2024; the arXiv is the durable reference).
- **Keller–Vossieck, "Aisles in derived categories", Bull. Soc. Math. Belg. Sér. A
  40 (1988) 239–253** — the origin of silting objects (AI Def 2.1 attribution).
  Add if BibTeX-verifiable at merge; else cite AI as the primary.

**Transcription caveat (Ruling 8, honest).** The theorem/definition NUMBERS cited
throughout this plan (AI Def 2.1 / 2.30 / 2.34, Thm 2.26 / 2.27 / 2.31 / 2.35, Cor
2.28 / 2.43, Prop 2.23(b) / 2.33, Thm 1.2, Examples 2.45 / 2.47; Oppermann Thm 1.1 /
Example 2.3; Jørgensen) are transcribed from the primary sources read this session but
**remain to-be-spot-checked at implementation** against the published (not just arXiv)
numbering — arXiv and journal versions occasionally renumber. The MATH each label names
is the load-bearing content and is what the oracles pin; a mis-transcribed number is a
citation-hygiene fix, not a design change. BibTeX venues carry their own
verify-at-merge flags (below and Task 6).

**The K₀ basis (BLOCKING correction — the g-matrix must live in the projective basis).**
`K₀(K^b(proj A)) = ⊕_v Z[P_v]` — the free abelian group on the **indecomposable
projectives**. A silting object's g-matrix is the matrix of its summands' classes in
THIS basis; Thm 2.27's "Z-basis of K₀" and the unimodularity `det = ±1` are statements
about the **projective** basis. **`derived/tilting.py::_chi` computes the wrong thing:**
`_chi(cx)` sums `(-1)^n · dimension_vector(cx_n)`, and a module's `dimension_vector` is
its class in the **composition-factor** basis `⊕_v Z[S_v]`. The change-of-basis from the
projective to the composition-factor basis is the **Cartan matrix** `C` (each `[P_v]`
maps to the `v`-th Cartan column), so `_chi`'s g-matrix equals `C · g_proj` and
`det(_chi) = det(C) · det(g_proj)`. For a **non-unimodular-Cartan** algebra (every
self-injective / symmetric algebra, and the plan's headline Example 2.47 where
`det C = 0`) this makes `det(_chi) ≠ ±1` **even for the regular object `A = ⊕ P_v`**,
whose true g-matrix is the identity. **Verified live (this session):** over Example 2.47
(`1 ⇄ 2`, `ab = ba = 0`), `_chi` gives the singular Cartan `[[1,1],[1,1]]`, `det = 0`,
and `is_tilting_complex([P₁,P₂]).generates = False` — a systematic false-negative
(likewise `NakayamaAlgebra([2,2,2], cyclic)`, `det C = 2`, and every `kZ_n/J^L`).
- **Fix (Task 0): compute the g-matrix with a new `g_proj(cx, verts)`** — the K₀ class
  in the **projective** basis: `Σ_n (-1)^n · (multiplicity of P_v in cx_n)`. Read the
  per-degree projective-summand multiplicities either from the provenance
  `_proj_vertices` when present (`from_projective_resolution` / `_direct_sum_complex`
  carry it — tilting.py:76–80) **or, for provenance-free summands such as
  `ChainComplex.stalk(A.projective(v))`, from the top of each (projective) term**
  (`top(P_v^{m_v}) = S_v^{m_v}`, so `dim-vec(top(cx_n))_v = mult of P_v`). Both routes
  give the same integer multiset; the top route is the robust primary since the plan's
  test summands are stalks that do **not** carry `_proj_vertices` (confirmed live).
- **Necessary (Thm 2.27):** a silting object's `n` summands are a Z-basis of `K₀`
  (g_proj square + `det(g_proj) = ±1`). The plan computes this exactly via `g_proj` +
  sympy. **Verified live:** `det(g_proj) = 1` for the regular object over both Example
  2.47 and `NakayamaAlgebra([2,2,2], cyclic)` — the identity matrix, as it must be for
  `A = ⊕ P_v` over EVERY algebra (Aihara–Iyama Thm 2.27's actual invariant).
- **Not sufficient:** thick subcategories of `K^b(proj A)` are **not** classified by
  `K₀` — J. Krah's example gives a f.d. algebra whose `K^b(proj)` has a **phantom**
  (a nonzero object with zero `K₀`-class), and there exist **presilting objects that
  are not completable to silting** (arXiv:2302.12502, arXiv:2304.08417). So
  `presilting + det ±1` alone does **not** imply `thick(T) = K^b(proj)`.
- **Where it IS decidable (certified TRUE):** (a) `T` is a **tilting** complex
  (reuse `is_tilting_complex`; tilting ⟹ silting — Example 2.2 covers `T = A`);
  (b) `T` is **2-term** presilting with `n` summands (AIR / Iyama–Jørgensen–Yang:
  every 2-term presilting completes to 2-term silting — cross-checked vs P45); (c)
  `A` is **local** and `T` indecomposable (Thm 2.26 → `T = A[i]`; covers `k[x]/(x^n)`).
- **Otherwise: `is_silting = "unknown"`** — the honest three-valued verdict:
  "presilting + K₀-basis (necessary conditions for silting) hold; thick-generation is
  not independently certified for this algebra class (K₀ does not detect thick
  subcategories in general — Krah phantom)." This is the metaplan-§1.3 loud-scope
  contract applied to a genuinely-undecided-by-K₀ question — never a silent `True`.
  (Gentle algebras also satisfy the sufficiency — presilting with `r = rk K₀` summands
  is silting — but a shipped gentle recognizer is P68's; not relied on here.)

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **Prerequisites: P43 (derived stack) AND P45 (τ-tilting) MERGED to `dev`** — both
  are confirmed merged (this is a Wave-3 plan; `src/quiverlab/derived/` and
  `src/quiverlab/tautilting/` are present). Branch `plan-67-silting` off `dev`.
  P67 is **[independent]** in the wave map (P79 later depends on P67).
- **Test bucket (check `tests/conftest.py::_DEEP_DIRS`).** `_DEEP_DIRS = (engine,
  resolutions_cs, modules, families, batch, specseq)` — **there is NO `tests/derived`
  dir and it is NOT a deep dir.** Follow the P43 precedent exactly: **all P67 engine
  tests live in `tests/modules/` as `tests/modules/test_silting_*.py`** (→ **deep**
  bucket). `tests/webapp/` + `tests/gui/` are extras-gated (unmarked); `tests/qpa/`
  → qpa. Do **not** create a new top-level tests dir. Run new tests by path during
  development; finish each task with a `-m deep` (or the extras-gated dir) spot-run.
- **The `decompose` char caveat is load-bearing.** The complex-level approximation
  and any indecomposability check inherit `decompose`'s scope (rigorous over char 0
  or char > dim; loud `QuiverlabError` otherwise — `modules/decompose.py`). **Every
  battery runs over QQ or GF(32003)** so both the split search and the locality
  certificate decide; state this per battery. Silting summands are **supplied
  indecomposable** (the verifier does not re-split them), but the approximation's
  Hom-space coordinatisation and the mutant's re-verification stay in scope.
- House conventions: homological complexes (`d_n: C_n → C_{n-1}`, rows = target —
  P39); a chain map's component is a `tgt.dim × src.dim` matrix over the shared
  Domain; composition is left-to-right (`f.then(g)`). `ChainMap.then` requires the
  middle complex be the **same Python object** (`g.src is self.tgt`) — pass the SAME
  summand objects into the approximation so `hyper_hom_basis`-produced maps compose.
  All refusals are `QuiverlabError`; `check=False` only for internally-constructed
  data whose certificate is asserted separately.
- **Honest scope (metaplan §6; verification page at merge).** (i) The silting quiver
  can be **infinite** (kA₂ already is) and mutation-transitivity is proven only for
  **local / hereditary / canonical** (AI Thm 1.2, fails for a symmetric algebra) —
  **no general BFS enumeration claim**; the exploration is bounded-radius with a loud
  `status`, certified **complete only for local** (`finite_class ∈ {"local", None}` — no
  `"two_term"`; the 2-term slice is a direct P45 cross-check, not a walk mode, Ruling 7).
  (ii) Generation is three-valued — certified `True` only on tilting /
  2-term / local, else `"unknown"` (K₀ **projective-basis** det ±1 necessary,
  Krah-phantom scope). (iii) The
  co-t-structure is a **documentation record** (coheart + aisle descriptors +
  references), NOT a computed subcategory (the aisles are infinite). (iv) `End(μT)`
  vs Oppermann is verified at the **underlying-quiver** level; the full graded dg
  rule is # PIN'd. (v) **QPA has no silting surface** — honest guard, no `qpa`
  battery (probe-first, FAIL if that ever changes).
- Plan-32 markers: presilting window scan / mutant re-verification / triangle
  `d²=0` / shares-`n-1` / involution `μ^-∘μ^+ = id` = `oracle_selfcert`;
  2-term silting ≡ P45 τ-tilting, and `End(μT)` quiver ≡ Oppermann rule =
  `oracle_crossengine`; AI worked values (k[x]/(x²) = shifts, kA₂ Example 2.45
  quiver, Example 2.47 cones, transitivity classes) = `oracle_literature`; **no
  `qpa` bucket** (QPA cannot compare — stated on the verification page). Every merge
  updates `docs/verification.md` (new oracle rows + recounted class table green,
  `tests/release/test_oracle_classes.py`) and adds citations to `references.bib`
  (BibTeX-verified) + `registry.py`.
- **Merge-train count note:** v1.0.0 runs Plans 51–79 in waves; suite/oracle-class
  counts drift as siblings merge. Recount from a live collection at *this* plan's
  merge (Task 6); note "counts as of the P67 merge; sibling plans in flight may shift
  them", matching the P43/P45 recount discipline.
- Conventional commits; green tests at every commit.

---

### Task 0: fix the inherited K₀-basis bug in `derived/tilting.py` (P43 g-matrix on the wrong basis)

**Why FIRST (BLOCKING).** P67's verifier reuses P43's generation machinery, and that
machinery is wrong for any non-unimodular-Cartan algebra (see the sufficiency ruling
above). `is_tilting_complex` computes its g-matrix with `_chi` (the composition-factor
basis), so `generates` is `det(C·g_proj) = det C · det g_proj`, not `det g_proj`. **This
is a live bug in shipped P43 code**, and it makes P67's tilting sufficiency rung
(`is_silting_object` step 1) reject the regular object `A = ⊕ P_v` over exactly the
algebras P67 must handle (self-injective/symmetric, Example 2.47). It is fixed once,
here, in `derived/tilting.py`, before any P67 task consumes it.

**Verified live (this session), pins to encode:**
- Example 2.47 (`1 ⇄ 2`, `ab = ba = 0`): `cartan_matrix = [[1,1],[1,1]]`, `det C = 0`;
  OLD `is_tilting_complex([P₁,P₂])` → `is_tilting=False, rigid=True, generates=False,
  det=0`. With `g_proj`: g-matrix `[[1,0],[0,1]]`, `det=1`, `generates=True`,
  `is_tilting=True`.
- `NakayamaAlgebra([2,2,2], cyclic=True)` (self-injective, 3 simples, `det C = 2`): OLD
  `generates=False, det=2`; `g_proj` → `det=1, generates=True`. (Also `[3,3,3]`:
  `det C = 0` → `det g_proj = 1`; `[3,3]`: `det C = 3` → `1`.)

**Files:**
- Modify: `src/quiverlab/derived/tilting.py` (add `g_proj`; route `is_tilting_complex`'s
  generation through it; respec `TiltingReport.g_matrix`/`det` semantics onto the
  projective basis).
- Test: `tests/modules/test_derived_tilting.py` (add the non-unimodular-Cartan
  regression oracle; keep every existing test green).

- [ ] **Step 1: Add `g_proj` and route `is_tilting_complex` through it.**

```python
def g_proj(cx, verts):
    """K0 class of a perfect complex in the PROJECTIVE basis of K0(K^b proj A) =
    (+)_v Z[P_v]: g_proj(cx)[v] = sum_n (-1)^n * (multiplicity of P_v in cx.term(n)).
    The multiplicity of P_v in a projective Q is dim_k (top Q)_v (top(P_v^{m}) = S_v^{m}).
    Uses the `_proj_vertices` provenance when present (from_projective_resolution /
    _direct_sum_complex); else reads it from the top of each term -- so it is correct for
    a bare ``ChainComplex.stalk(A.projective(v))`` (which carries NO provenance)."""
    idx = {v: i for i, v in enumerate(verts)}
    out = [0] * len(verts)
    prov = getattr(cx, "_proj_vertices", None)
    for n in cx.degrees():
        if prov is not None and n in prov:
            mult = {}
            for v in prov[n]:
                mult[v] = mult.get(v, 0) + 1
        else:
            mult = cx.term(n).top().dimension_vector()
        for v, m in mult.items():
            if v in idx:
                out[idx[v]] += (-1) ** n * m
    return out
```

  In `is_tilting_complex`, replace the generation leg
  `g = [_chi(T, verts) for T in summands]` with `g = [g_proj(T, verts) for T in
  summands]` (det + `generates` unchanged in shape — `len(g) == len(verts)` and
  `det ∈ {±1}`). Update the `TiltingReport` field docstrings: `g_matrix` rows are the
  summand K₀ classes in the **projective** basis; `det` is `det(g_proj)`. `_chi` (the
  composition-factor χ) then has no remaining consumer in `src/` (verified: only
  `is_tilting_complex` used it) — remove it, or retain it with a docstring stating it is
  the composition-factor class and is NOT the K₀-basis used for generation. Do **not**
  leave two silently-different g-matrices in play.

- [ ] **Step 2: Regression oracle + no-regression gate.**

```python
# tests/modules/test_derived_tilting.py  (append)
@pytest.mark.oracle_selfcert
def test_regular_object_is_tilting_over_nonunimodular_cartan():
    # The g-matrix must live in K0(K^b proj) = (+)_v Z[P_v], NOT the composition-factor
    # basis: the regular object A = (+) P_v is a tilting complex over EVERY algebra
    # (AI Ex 2.2). The old _chi g-matrix gave det(Cartan) and FAILED here (det C = 0/2).
    import sympy as sp
    from quiverlab import GF, NakayamaAlgebra, Quiver
    from quiverlab.modules.complexes import ChainComplex
    from quiverlab.derived.tilting import is_tilting_complex
    # Example 2.47: det Cartan = 0
    A = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)
    assert int(sp.Matrix(A.cartan_matrix()).det()) == 0        # non-unimodular
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    rep = is_tilting_complex(T)
    assert rep.is_tilting is True and rep.generates is True and rep.det in (1, -1)
    # self-injective Nakayama kZ3/J2: det Cartan = 2
    B = NakayamaAlgebra([2, 2, 2], cyclic=True, field=GF(32003))
    assert int(sp.Matrix(B.cartan_matrix()).det()) == 2
    TB = [ChainComplex.stalk(B.projective(v), 0) for v in B.quiver.vertices]
    repB = is_tilting_complex(TB)
    assert repB.is_tilting is True and repB.det in (1, -1)
```

  **No-regression:** the existing `tests/modules/test_derived_tilting.py` +
  `test_derived_tau.py` must stay byte-for-byte green. They only exercise **kA₂ and its
  APR tilts** — hereditary, so the Cartan is unimodular (`det = 1`) and `_chi` and
  `g_proj` COINCIDE (`det _chi = det C · det g_proj = 1 · det g_proj`). **That is exactly
  why the bug survived P43's green suite** — state this in the test docstring. Run
  `... -m pytest tests/modules/test_derived_tilting.py tests/modules/test_derived_tau.py -q`.

- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/derived/tilting.py tests/modules/test_derived_tilting.py
git commit -m "fix(derived): tilting-complex g-matrix on the projective K0 basis (g_proj), not the composition-factor basis -- P43 was a systematic false-negative on non-unimodular Cartan (self-injective/symmetric); regression pinned over Example 2.47 + kZ3/J2"
```

---

### Task 1: `silting.py` — the silting-object verifier + the co-t-structure dictionary

**Files:**
- Create: `src/quiverlab/derived/silting.py`
- Modify: `src/quiverlab/derived/__init__.py` (export the new public surface)
- Test: `tests/modules/test_silting_verifier.py`

**Interfaces:**
- Consumes: `derived/tilting.py` — `is_tilting_complex`, `TiltingReport`, `g_proj`
  (the Task-0 K₀-projective-basis helper), and the private helpers `_direct_sum_complex`,
  `_span` (reuse; do NOT re-derive the K₀ / window machinery, and do NOT use `_chi` for
  generation — it is the wrong basis); `modules/complexes.py` — `hyper_hom_dims`,
  `ChainComplex.is_perfect`/`degrees`/`term`; `sympy` for the integer determinant.
- Produces:
  ```python
  @dataclass
  class SiltingReport:
      is_silting: object          # True | False | "unknown"  (three-valued; see ruling)
      is_presilting: bool         # Hom_{D^b}(T, T[n]) == 0 for all n > 0 in window
      k0_basis: bool              # g_proj square (#summands == #simples) AND det(g_proj) +-1
      window: tuple               # (1, n_max): the EXACT positive rigidity-check range
      g_matrix: list              # rows = summand K0 classes in the PROJECTIVE basis
                                  # (+)_v Z[P_v] via g_proj -- NOT the composition-factor
                                  # basis (see Task 0 / the sufficiency ruling)
      det: int                    # det(g_proj) -- unimodularity in the projective basis
      generation_certified_by: str  # "tilting (AI Ex 2.2)" | "2-term (AIR/IJY)" |
                                     # "local (AI Thm 2.26)" | "K0-basis only (necessary;
                                     # not sufficient -- Krah phantom)"

  def is_silting_object(summands) -> SiltingReport
      # summands: a list of INDECOMPOSABLE perfect complexes (the intended silting
      # summands). Presilting is DECIDED on the exact positive window (perfect =>
      # bounded => finite window, outside it hyper-Hom is provably 0). Generation is
      # three-valued per the sufficiency ladder. Loud if any summand is not perfect.

  def co_t_structure_of(summands) -> dict
      # The co-t-structure dictionary (AI Prop 2.23(b) / Jorgensen). Documentation
      # record ONLY: {"coheart": [<summand labels/dimvecs>], "aisle": <descriptor
      # string>, "coaisle": <descriptor string>, "bounded": True, "references": [...]}.
      # NOT a computed subcategory (the aisles are infinite) -- refuses loudly if
      # `summands` is not a certified silting object (is_silting in (True,"unknown")).
  ```

**The positive rigidity window (decided, not semi-decided).** For perfect summands
each spanning degrees `[lo_i, hi_i]`, `Hom^n(T_i, T_j)` is the zero cochain group
unless some `p` has both `T_i` at `p` and `T_j` at `p−n` nonzero, i.e.
`n ∈ [lo_i − hi_j, hi_i − lo_j]`. The honest **positive** window is
`n_max = max_i hi_i − min_j lo_j`; for `n > n_max` hyper-Hom is provably 0, so
presilting is fully decided by scanning `n ∈ [1, n_max]`. Report `window = (1, n_max)`
(metaplan honesty: name what was checked). Reuse `_direct_sum_complex` + `_span` from
`derived/tilting.py`; scan `hyper_hom_dims(Tsum, Tsum, 1, n_max)`.

**The sufficiency ladder (`generation_certified_by`).** In order:
1. `is_tilting_complex(summands).is_tilting` → `is_silting = True`,
   `generation_certified_by = "tilting (AI Ex 2.2)"` (tilting ⟹ silting; covers
   `T = A`, the regular silting object).
2. else if presilting ∧ **2-term** (every summand's degrees ⊆ two consecutive integers,
   up to a common shift) ∧ `k0_basis` → `True`, `"2-term (AIR/IJY)"`.
   **This is the genuine IJY completion criterion** (Iyama–Jørgensen–Yang / AIR): a
   2-term **presilting** object whose `#(distinct indec) summands == #simples (= rk K₀)`
   IS 2-term silting. The rung condition is `presilting ∧ 2-term ∧ k0_basis`, and
   `k0_basis = (#summands == #simples) ∧ det(g_proj) = ±1`: the `det = ±1` half is
   automatic for a genuine (basic) 2-term silting object, and it is **load-bearing** for
   the guard — it correctly rejects degenerate NON-basic inputs (a repeated summand
   `[P, P]` has `#summands == #simples` but `det(g_proj) = 0`, so it falls through to
   `False`, never a false `True` or a crash). Do NOT drop `det` from the condition.
3. else if presilting ∧ **`A` local** (one simple) ∧ `#summands == 1` → `True`,
   `"local (AI Thm 2.26)"`.
4. else if presilting ∧ `k0_basis` → `"unknown"`,
   `"K0-basis only (necessary; not sufficient -- Krah phantom)"`.
5. else (`not presilting`, or `#summands != n`, or `det ∉ {±1}`) → `False`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_silting_verifier.py
"""Silting-object verifier (Plan 67 / AI Def 2.1). Self-cert: T = A is silting (in fact
tilting) with the positive window (1,0); T = P2 (+) P1[1] over kA2 is a GENUINE
silting-but-not-tilting object (presilting via the clean positive window, but the P43
tilting verifier fails rigidity at n = -1) -- the silting-vs-tilting separation; a
single projective P1 over kA2 is presilting but NOT silting (#summands != #simples) --
the presilting-vs-silting (generation) separation; a single projective stalk over a
local algebra is silting (Thm 2.26). SMOKE (would have caught the Task-0 K0-basis bug):
the regular object over the NON-unimodular-Cartan Example 2.47 and a self-injective
Nakayama both verify silting=True. Literature: k[x]/(x^2) silting = shifts. Honest
three-valued: the K0-basis-only case reports 'unknown', never a silent True.

WHY the kA2/local oracles ALONE mask the K0-basis bug (Task 0): kA2 is hereditary so its
Cartan is unimodular (det = 1), where g_proj and the old composition-factor _chi
coincide; a LOCAL algebra hits the #summands==1 'local' branch that never touches the
K0 g-matrix. Only a non-unimodular-Cartan, non-local input (Example 2.47, kZ3/J2)
exercises the projective-basis g_proj -- hence the two smoke tests below."""
import pytest

from quiverlab import GF, NakayamaAlgebra, Quiver, linear_path_algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.silting import (is_silting_object, co_t_structure_of)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _a2():
    return linear_path_algebra(2, field=QQ)          # kA2, 1->2, hereditary


def _comm_square():
    # 1 <=> 2 with ab = ba = 0 (AI Example 2.47) -- det Cartan = 0 (NON-unimodular).
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


@selfcert
def test_regular_module_is_silting_and_tilting():
    A = _a2()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    rep = is_silting_object(T)
    assert rep.is_presilting and rep.k0_basis
    assert rep.is_silting is True
    assert rep.generation_certified_by.startswith("tilting")
    assert rep.window == (1, 0)                       # no positive shifts to scan
    assert rep.det in (1, -1)                         # g_proj identity (projective basis)


@selfcert
def test_silting_not_tilting_separation():
    # T = P2 (+) P1[1] over kA2 is a GENUINE silting object that is NOT tilting.
    # Presilting: the positive window (1,1) is clean (Hom_{D^b}(T, T[1]) = 0). It is
    # silting via the 2-term (IJY) rung. But it is NOT tilting: the P43 tilting verifier
    # fails rigidity at n = -1 (Hom_{D^b}(T, T[-1]) = End(P1)-arrow != 0). This is the
    # vanishing-condition separation: silting checks only n > 0, tilting checks n != 0.
    # (Live-verified: full window {-1:1, 0:2, 1:0}; is_presilting True; det g_proj = 1.)
    A = _a2()
    P1 = ChainComplex.stalk(A.projective(1), 0)
    P2 = ChainComplex.stalk(A.projective(2), 0)
    T = [P2, P1.shift(1)]
    rep = is_silting_object(T)
    assert rep.is_presilting is True                  # positive window is clean
    assert rep.window == (1, 1)
    assert rep.is_silting is True                     # IJY 2-term rung
    assert rep.generation_certified_by.startswith("2-term")
    # but it is NOT tilting -- the negative window is nonzero:
    from quiverlab.derived.tilting import is_tilting_complex
    assert is_tilting_complex(T).rigid is False       # fails rigidity at n = -1


@selfcert
def test_single_projective_is_presilting_not_silting():
    # A single indecomposable P1 over kA2: PRESILTING (positive window (1,0) is
    # vacuously clean) but NOT silting -- #summands (1) != #simples (2), so the classes
    # cannot be a K0-basis (k0_basis False). The presilting-vs-silting separation.
    A = _a2()
    rep = is_silting_object([ChainComplex.stalk(A.projective(1), 0)])
    assert rep.is_presilting is True
    assert rep.window == (1, 0)
    assert rep.k0_basis is False and rep.is_silting is False


@selfcert
def test_verifier_smoke_nonunimodular_cartan():
    # SMOKE (Ruling 6a): the ONE probe that would have caught the Task-0 K0-basis bug.
    # Over Example 2.47 (det Cartan = 0) the regular object A = P1 (+) P2 IS silting
    # (in fact tilting, AI Ex 2.2). The old _chi g-matrix gave det = det Cartan = 0 and
    # a FALSE is_silting=False. (Live-verified post-fix: is_silting True, det g_proj 1.)
    A = _comm_square()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    rep = is_silting_object(T)
    assert rep.is_silting is True and rep.det in (1, -1)
    assert rep.generation_certified_by.startswith("tilting")


@selfcert
def test_selfinjective_nakayama_regular_is_silting():
    # SMOKE (Ruling 6b): a NON-local, NON-unimodular-Cartan self-injective algebra.
    # kZ3/J2 = NakayamaAlgebra([2,2,2], cyclic): 3 simples, det Cartan = 2, self-injective.
    # The regular object is silting (tilting). Old _chi gave det = 2 -> false negative.
    B = NakayamaAlgebra([2, 2, 2], cyclic=True, field=GF(32003))
    T = [ChainComplex.stalk(B.projective(v), 0) for v in B.quiver.vertices]
    rep = is_silting_object(T)
    assert rep.is_silting is True and rep.det in (1, -1)


@lit
def test_local_algebra_silting_is_shifts():
    # k[x]/(x^2): local, indecomposable silting object A -> silt = {A[i]} (Thm 2.26).
    from quiverlab.families import truncated_polynomial
    A = truncated_polynomial(2, field=GF(32003))
    Astalk = ChainComplex.stalk(A.projective(1), 0)
    rep = is_silting_object([Astalk])
    assert rep.is_silting is True
    assert rep.generation_certified_by.startswith(("tilting", "local"))
    # a shift A[3] is also silting:
    assert is_silting_object([Astalk.shift(3)]).is_silting is True


@selfcert
def test_missing_summand_is_not_silting():
    A = _a2()
    T = [ChainComplex.stalk(A.projective(1), 0)]      # one summand, two simples
    rep = is_silting_object(T)
    assert rep.k0_basis is False and rep.is_silting is False


@selfcert
def test_nonperfect_summand_refused():
    A = _a2()
    T = [ChainComplex.stalk(A.simple(1), 0)]          # simple: not projective
    with pytest.raises(QuiverlabError, match="perfect"):
        is_silting_object(T)


@selfcert
def test_co_t_structure_record_and_refusal():
    A = _a2()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    ct = co_t_structure_of(T)
    assert ct["bounded"] is True and len(ct["coheart"]) == 2
    assert "aihara" in " ".join(ct["references"]).lower() or ct["references"]
    # refuses on a non-silting input (does not fabricate a co-t-structure):
    with pytest.raises(QuiverlabError):
        co_t_structure_of([ChainComplex.stalk(A.projective(1), 0)])
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.derived.silting`.

- [ ] **Step 3: Implement `src/quiverlab/derived/silting.py`** (verifier + co-t-structure).

```python
"""Silting objects in K^b(proj A): verifier, single mutation, bounded exploration,
and the co-t-structure dictionary (Plan 67 / Aihara-Iyama arXiv:1009.3370).

A silting object T satisfies Hom_{D^b}(T, T[n]) = 0 for all n > 0 (PRESILTING) and
thick(T) = K^b(proj A) (GENERATION). This is WEAKER than tilting (which also needs
n < 0 vanishing); tilting => silting. Presilting is DECIDED on the exact positive
window (perfect => bounded => finite). Generation is three-valued: certified True on
the tilting / 2-term / local classes (where a completion theorem applies), else
'unknown' -- det(g_proj) +-1 is NECESSARY (Thm 2.27) but not sufficient in general
(thick subcategories are not K0-classified; Krah phantom). The g-matrix is computed in
the PROJECTIVE basis K0(K^b proj A) = (+)_v Z[P_v] via `g_proj` (Task 0), NOT the
composition-factor basis -- `_chi` gives det(Cartan.g) and is wrong on non-unimodular
Cartan (self-injective/symmetric; Example 2.47). Float-free; every mutant is re-verified
silting."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import hyper_hom_dims
from quiverlab.derived.tilting import (is_tilting_complex, _direct_sum_complex, _span,
                                       g_proj)


@dataclass
class SiltingReport:
    is_silting: object
    is_presilting: bool
    k0_basis: bool
    window: tuple
    g_matrix: list
    det: int
    generation_certified_by: str


def _positive_window(summands):
    spans = [_span(T) for T in summands]
    return 1, max(hi for _, hi in spans) - min(lo for lo, _ in spans)


def _is_two_term(summands):
    """Every summand concentrated in two consecutive degrees, all up to ONE common
    shift (2-term silting is defined up to shift). True for the empty/stalk cases."""
    spans = [_span(T) for T in summands if T.degrees()]
    if not spans:
        return True
    widths = {hi - lo for lo, hi in spans}
    if not (widths <= {0, 1}):
        return False
    tops = {hi for _, hi in spans}
    return max(tops) - min(tops) <= 1          # a common 2-window covers them


def is_silting_object(summands):
    import sympy as sp
    if not summands:
        raise QuiverlabError("is_silting_object: need at least one summand")
    for T in summands:
        if not T.is_perfect():
            raise QuiverlabError("is_silting_object: every summand must be a "
                                 "certified perfect complex")
    A = summands[0].algebra
    verts = list(A.quiver.vertices)
    n_lo, n_max = _positive_window(summands)
    Tsum = _direct_sum_complex(summands)
    hh = hyper_hom_dims(Tsum, Tsum, n_lo, n_max) if n_max >= n_lo else {}
    presilting = all(hh.get(n, 0) == 0 for n in range(n_lo, n_max + 1))
    g = [g_proj(T, verts) for T in summands]       # K0 in the PROJECTIVE basis (Task 0)
    det = int(sp.Matrix(g).det()) if len(g) == len(verts) else 0
    k0_basis = (len(g) == len(verts)) and det in (1, -1)
    # sufficiency ladder (see the module/plan sufficiency ruling)
    is_silting, why = False, ""
    if is_tilting_complex(summands).is_tilting:     # is_tilting_complex now g_proj-routed
        is_silting, why = True, "tilting (AI Ex 2.2)"
    elif presilting and _is_two_term(summands) and k0_basis:
        # IJY: a 2-term PRESILTING object whose #(distinct indec) summands == #simples is
        # silting. k0_basis = (#summands == #simples) AND det(g_proj) in {+-1}: the det
        # condition is automatic for a genuine (basic) 2-term silting object AND it
        # correctly REJECTS degenerate non-basic inputs (e.g. a repeated summand [P,P]
        # has #summands == #simples but det(g_proj) = 0 -> falls through to False, never
        # a crash). Do NOT drop det from the condition.
        is_silting, why = True, "2-term (AIR/IJY completion)"
    elif presilting and len(verts) == 1 and len(summands) == 1:
        is_silting, why = True, "local (AI Thm 2.26)"
    elif presilting and k0_basis:
        is_silting = "unknown"
        why = "K0-basis only (necessary; not sufficient -- Krah phantom)"
    return SiltingReport(is_silting=is_silting, is_presilting=presilting,
                         k0_basis=k0_basis, window=(n_lo, n_max), g_matrix=g,
                         det=det, generation_certified_by=why)


_CT_REFS = ["aihara_iyama_silting", "jorgensen_cotstructures"]


def co_t_structure_of(summands):
    """The bounded co-t-structure with coheart add(T) (AI Prop 2.23(b)). A
    documentation record -- NOT a computed subcategory (the aisles are infinite)."""
    rep = is_silting_object(summands)
    if rep.is_silting not in (True, "unknown"):
        raise QuiverlabError(
            "co_t_structure_of: input is not a (certified or candidate) silting "
            "object; no co-t-structure to report", hint=str(rep.generation_certified_by))
    coheart = [{"dimvecs": {n: T.term(n).dimension_vector() for n in T.degrees()}}
               for T in summands]
    return {"coheart": coheart, "bounded": True,
            "aisle": "T_M^{<=0} = { X : Hom_{D^b}(T, X[>0]) = 0 }  (AI Def 2.12)",
            "coaisle": "^perp(T_M^{<=0})  (AI Prop 2.23(b) torsion pair)",
            "coheart_is_addT": True, "references": list(_CT_REFS)}
```

**Adjust to reality (Task 1):**
- `_direct_sum_complex`, `_span` are private in `derived/tilting.py` (confirmed present:
  `derived/tilting.py:45-103`). Import them; do **not** copy. `g_proj` is the **new**
  Task-0 helper in the same module — import it too; do NOT re-derive a g-matrix here.
- **`g_proj` must handle stalk summands (no `_proj_vertices`).** The Task-1 tests build
  summands via `ChainComplex.stalk(A.projective(v), 0)`, which carries **no**
  `_proj_vertices` provenance (confirmed live). `g_proj` therefore reads projective
  multiplicities from the **top of each term** in that case (`top(P_v^{m}) = S_v^{m}`),
  and from `_proj_vertices` only as a fast path when present. Do NOT rely on
  `_direct_sum_complex` carrying provenance for stalk inputs — it will not.
- `is_tilting_complex` (now Task-0 g_proj-routed) internally scans `n != 0` and builds
  its own `_direct_sum_complex` + `hyper_hom_dims` — calling it in the ladder
  double-computes the direct sum once. That is acceptable (path-algebra hyper-Hom is
  cheap); if profiling shows it, factor a private `_rigidity(summands, lo, hi)` shared by
  both. Do **not** create a second independent hyper-Hom route.
- The `window == (1, 0)` for `T = A` (all stalks in degree 0): `n_max = 0 - 0 = 0 < 1`,
  so the scan range is empty and `presilting = True` vacuously — correct
  (`Ext^{>0}(P, P) = 0`).
- `co_t_structure_of` ships the coheart's per-summand dim-vectors (the durable,
  finite data); the aisle/coaisle are **descriptor strings** with the AI reference,
  never enumerated. This is the honest v1 (a fuller weight-filtration is a strict
  superset for a later plan).

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/derived/silting.py src/quiverlab/derived/__init__.py \
        tests/modules/test_silting_verifier.py
git commit -m "feat(derived): silting-object verifier (presilting positive window + three-valued generation) + co-t-structure dictionary"
```

---

### Task 2: `silting.py` (mutation) — single silting mutation via one approximation triangle

**Files:**
- Modify: `src/quiverlab/derived/silting.py` (add the mutation + neighbours)
- Test: `tests/modules/test_silting_mutation.py`

**Precondition (depends on Task 0).** The whole mutation battery seeds from the regular
object over Example 2.47 (`_comm_square`) and asserts `is_silting_object([P₁,P₂]).is_silting
is True`. That assertion **only holds once Task 0's g_proj fix is in** — with the old
`_chi` g-matrix, `det = det Cartan = 0` and the seed verified `False`, so `silting_mutate`
would refuse its own input at the guard and the battery could not run. Task 0 is a hard
prerequisite of Task 2 (live-verified: post-fix `is_silting_object([P₁,P₂]).is_silting is
True`, `generation_certified_by = "tilting (AI Ex 2.2)"`). The mutant re-verification
(`_assert_neighbour` below) likewise leans on the repaired certifier.

**The math (AI Def 2.30 + 2.34, verified).** Irreducible **left** mutation at an
indecomposable summand `X` of a silting object `T = X ⊕ rest`: take the minimal
**left `add(rest)`-approximation** `f: X → D` (`D ∈ add(rest)`), form the triangle
`X -f-> D -g-> N_X -> X[1]` (so `N_X = cone(f)`), and set
`μ_X^+(T) = rest ⊕ {N_X}`. By Thm 2.31 this is again silting. Right mutation `μ_X^-`
is dual: the minimal **right `add(rest)`-approximation** `f': D' → X`, cocone triangle
`N'_X -> D' -f'-> X -> N'_X[1]` (so `N'_X = cone(f')[-1]`), `μ_X^-(T) = rest ⊕ {N'_X}`.
AI's sign convention: **`μ^+` = LEFT** and `T ≥ μ^+(T)` (Prop 2.33). **Primary =
left**; both ship.

**The complex-level minimal `add(rest)`-approximation (the new engine piece).** Lift
P44's retract-pruning **algorithm** (`modules/approximations.py::_left_approx_blocks`
/ `_prune`) from modules to `K^b(proj A)`:
- `types = rest` (the silting summands `≠ X`; already indecomposable — no `decompose`);
- `from_X[t] = hyper_hom_basis(X, rest[t], 0)` — degree-0 chain maps `X → rest_t`;
- `between[s][t] = hyper_hom_basis(rest[s], rest[t], 0)`; compose with `ChainMap.then`
  (**pass the SAME `rest` objects** so `then`'s `g.src is self.tgt` identity holds);
- `blocks = [(t, g) for t in range(len(rest)) for g in from_X[t]]` (the full,
  trivially-surjective start); `is_approx(bl)`: for every target `rest_s`, the image
  span of `{ g.then(chi) : (t, g) in bl, chi in between[t][s] }` fills
  `Hom_{D^b}(X, rest_s)` (= `hyper_hom_basis(X, rest_s, 0)`). "Fills" is decided in
  **homotopy coordinates**: coordinatise each `ChainMap` in the degree-0 hyper-Hom
  basis via `derived/tilting._cochain_vec` (which already handles the coset/coboundary
  reduction), then `mat_rank == dim Hom_{D^b}(X, rest_s)` (reuse the P44 `_surjects`).
- `_prune` in place → the minimal (right/left-minimal ⟹ unique-up-to-iso, so
  `N_X` indecomposable — AI Thm 2.35 proof); assemble `f: X → D = (+) rest_t` (over
  surviving blocks) as a `ChainMap`; `N_X = f.cone()`.

**Interfaces:**
- Consumes: `derived/homs.hyper_hom_basis`, `ChainMap.then`/`.cone`/`.triangle`
  (`modules/complexes.py:448-519`), `derived/tilting._cochain_vec`,
  `modules/linalg_mod` (`mat_rank`, `cols_to_matrix`), `is_silting_object` (Task 1).
- Produces:
  ```python
  def silting_mutate(summands, i, direction="left") -> list   # the mutated summand list
      # Irreducible mutation at summand index i (AI Def 2.34). direction in {"left","right"}.
      # SELF-CERTIFIED: input is silting (or presilting+k0_basis; loud otherwise); the
      # result RE-VERIFIES as silting, shares exactly n-1 summands with the input
      # (Hasse-neighbour), and differs from it. Returns the new list [rest..., N_X].
  def silting_neighbors(summands, direction="left") -> list   # [ (i, mutant) for i in range(n) ]
      # the (up to) n irreducible mutations -- the single-step neighbourhood. A summand
      # whose approximation degenerates (mutation not defined at that summand) is
      # reported with mutant=None + a reason (never a silent skip).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_silting_mutation.py
"""Single silting mutation via one approximation triangle (Plan 67 / AI Def 2.30/2.34).
Literature: the 1<=>2 algebra with ab=ba=0 (AI Example 2.47): mutating P1 (+) P2 at P1
gives P2 (+) cone(P1->P2), at P2 gives P1 (+) cone(P2->P1). Self-cert: every mutant is
silting, shares n-1 summands, differs from the input; left-then-right is the identity
(involution)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.silting import (is_silting_object, silting_mutate,
                                       silting_neighbors)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _comm_square():
    # 1 <=> 2 with ab = ba = 0 (AI Example 2.47): a: 1->2, b: 2->1.
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


def _key(cx):
    # a shift-insensitive fingerprint of a summand: its K0 class up to sign, plus the
    # sorted per-degree dim-vectors -- enough to recognise cone(P1->P2) etc.
    return (tuple(sorted((n, tuple(sorted(cx.term(n).dimension_vector().items())))
                         for n in cx.degrees())),)


@selfcert
def test_mutant_is_silting_and_shares_all_but_one():
    A = _comm_square()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    assert is_silting_object(T).is_silting is True
    for i in range(2):
        mut = silting_mutate(T, i, direction="left")
        assert is_silting_object(mut).is_silting in (True, "unknown")
        shared = {_key(c) for c in T} & {_key(c) for c in mut}
        assert len(shared) == 1                         # n - 1 == 1 summand shared
        assert {_key(c) for c in mut} != {_key(c) for c in T}


@lit
def test_example_2_47_cones():
    # AI Example 2.47: mu_{P1}^+(P1 (+) P2) = P2 (+) cone(P1->P2).
    A = _comm_square()
    P1 = ChainComplex.stalk(A.projective(1), 0)
    P2 = ChainComplex.stalk(A.projective(2), 0)
    mut = silting_mutate([P1, P2], 0, direction="left")   # mutate at P1 (index 0)
    # the shared summand is P2; the other is cone(P1 -> P2) (degrees {0,-1}) -- verify
    # via homology: cone of the P1->P2 approximation is a genuine 2-term complex, not a
    # stalk, and P2 survives.
    keys = {_key(c) for c in mut}
    assert _key(P2) in keys
    other = [c for c in mut if _key(c) != _key(P2)][0]
    assert set(other.degrees()) != {0}                    # a cone, not a stalk


@selfcert
def test_left_then_right_is_identity():
    # mu^-_?(mu^+_X(T)) recovers T at the matching summand (AI Prop 2.33 involution).
    A = _comm_square()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    mut = silting_mutate(T, 0, direction="left")
    # the new summand N_X is the one NOT shared with T; mutate it back RIGHT:
    shared = {_key(c) for c in T}
    j = [k for k, c in enumerate(mut) if _key(c) not in shared][0]
    back = silting_mutate(mut, j, direction="right")
    assert {_key(c) for c in back} == {_key(c) for c in T}


@selfcert
def test_neighbors_reports_all_summands():
    A = _comm_square()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    nb = silting_neighbors(T, direction="left")
    assert len(nb) == 2
    for (i, mutant) in nb:
        assert mutant is not None and is_silting_object(mutant).is_silting in (True, "unknown")
```

- [ ] **Step 2: Run to verify failure** — `ImportError: silting_mutate`.

- [ ] **Step 3: Implement** the complex-level approximation + `silting_mutate` /
  `silting_neighbors`, appended to `derived/silting.py`.

```python
def _left_add_approx_complex(rest, X):
    """Minimal left add(rest)-approximation f: X -> D (D in add rest), as a ChainMap,
    by the P44 retract-prune lifted to K^b(proj) over degree-0 hyper-Hom + ChainMap.then.
    `rest` are the (indecomposable, SAME-object) silting summands != X."""
    from quiverlab.derived.homs import hyper_hom_basis
    from quiverlab.derived.tilting import _cochain_vec
    from quiverlab.modules import linalg_mod as lm
    dom = X.domain
    from_X = [hyper_hom_basis(X, Rt, 0) for Rt in rest]           # X -> rest_t
    between = [[hyper_hom_basis(rest[s], rest[t], 0) for t in range(len(rest))]
              for s in range(len(rest))]
    target = [len(hyper_hom_basis(X, Rs, 0)) for Rs in rest]      # dim Hom(X, rest_s)
    blocks = [(t, g) for t in range(len(rest)) for g in from_X[t]]

    def is_approx(bl):
        for s in range(len(rest)):
            if target[s] == 0:
                continue
            cols = [_cochain_vec(g.then(chi), X, rest[s], dom)
                    for (t, g) in bl for chi in between[t][s]]
            if not cols or lm.mat_rank(lm.cols_to_matrix(cols), dom) != target[s]:
                return False
        return True

    # retract-prune (P44 _prune, inlined)
    changed = True
    while changed:
        changed = False
        for idx in range(len(blocks)):
            if is_approx(blocks[:idx] + blocks[idx + 1:]):
                del blocks[idx]; changed = True; break
    return _assemble_left_complex(rest, blocks, X, dom)            # ChainMap X -> (+)rest_t


def silting_mutate(summands, i, direction="left"):
    rep = is_silting_object(summands)
    if rep.is_silting not in (True, "unknown"):
        raise QuiverlabError("silting_mutate: input is not a (candidate) silting object",
                             hint=rep.generation_certified_by)
    n = len(summands)
    if not (0 <= i < n):
        raise QuiverlabError(f"silting_mutate: summand index {i} out of range 0..{n-1}")
    X = summands[i]
    rest = summands[:i] + summands[i + 1:]
    if direction == "left":
        f = _left_add_approx_complex(rest, X)                    # X -> D
        N = f.cone()                                             # N_X, triangle (4)
    elif direction == "right":
        fp = _right_add_approx_complex(rest, X)                  # D' -> X
        N = fp.cone().shift(-1)                                  # N'_X = cone(f')[-1]
    else:
        raise QuiverlabError("silting_mutate: direction must be 'left' or 'right'")
    N._perfect = True                                            # cone of perfects
    mutant = rest + [N]
    _assert_neighbour(summands, mutant, n)                       # self-cert (below)
    return mutant
```

**Adjust to reality (Task 2) — the certificate-arbitrated details:**
- **`_assemble_left_complex(rest, blocks, X, dom)`** mirrors P44's `_assemble_left`
  but over complexes: `D = _direct_sum_complex([rest[t] for (t,_g) in blocks])`
  (or the zero complex when `blocks == []`), and `f: X → D` has, per degree `p`, the
  stacked component matrices of the surviving `g`'s. Reuse the block-diagonal degree
  layout `_direct_sum_complex` produces so the target grading matches; then
  `ChainMap(X, D, comps, check=True)` re-certifies `f` is a genuine chain map. The
  empty-`blocks` case (`Hom(X, rest) = 0`) gives `D = 0`, `f: X → 0`, and
  `N_X = cone(f) = X[1]` — a shift, the correct degenerate mutation.
- **`_assert_neighbour(T, mutant, n)`**: (a) `len(mutant) == n`; (b)
  `is_silting_object(mutant).is_presilting` (Thm 2.31 guarantees it — this is the
  self-cert, not the proof); (c) the shift-insensitive summand fingerprints share
  exactly `n − 1` with `T` and the mutant `≠ T`. If (b)/(c) fail, the left/right
  branch or the approximation is wrong — **the involution + Example-2.47 oracles are
  the hard arbiters** (the P41/P44 "the certificate decides the convention" pattern);
  fix once and pin with a derivation comment. Do **not** special-case per summand.
- **`ChainMap.then` object-identity.** `hyper_hom_basis(X, rest[t], 0)` returns maps
  whose `tgt` **is** the passed `rest[t]` object (n=0 uses `Yn = Y`, confirmed
  `derived/homs.py:92`); `between[t][s] = hyper_hom_basis(rest[t], rest[s], 0)` has
  `src` **is** `rest[t]`. So `g.then(chi)` composes iff the **same** `rest[t]` Python
  object is used in both — pass `rest` (the caller's objects) straight through, never
  copies.
- **Right mutation** reuses the identical prune over the DUAL Hom direction
  (`_right_add_approx_complex(rest, X)`: `from_rest[t] = hyper_hom_basis(rest[t], X, 0)`,
  test `Hom(rest_s, X)`-surjectivity), assembling `f': D' → X`; `N'_X = f'.cone().shift(-1)`
  (the cocone). The `test_left_then_right_is_identity` involution pins the `shift(-1)`.
- **`silting_neighbors`** wraps `silting_mutate` over `i in range(n)` in a `try/except
  QuiverlabError`, recording `(i, None, reason)` on a summand where the mutation
  degenerates or refuses — never a silent skip.

- [ ] **Step 4: Run tests** — Expected: PASS (expect 1–2 iterations to lock the
  left/right + `shift(-1)` conventions against the Example-2.47 + involution oracles).
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/derived/silting.py tests/modules/test_silting_mutation.py
git commit -m "feat(derived): single silting mutation mu_X^± via one approximation triangle (complex-level add-approx + cone), AI Example-2.47 pinned"
```

---

### Task 3: `silting.py` (exploration) — bounded-radius exploration with loud truncation

**Files:**
- Modify: `src/quiverlab/derived/silting.py` (add the bounded exploration)
- Modify: `src/quiverlab/core/algebra.py` (thin `Algebra.silting_report` /
  `Algebra.silting_exploration` delegates — lazy-import, beside `global_dimension`;
  read the delegation idiom)
- Test: `tests/modules/test_silting_exploration.py`

**The honest contract (metaplan §6; the record's MAJOR re-scope).** There is **NO
general BFS** and **NO enumeration claim**: the silting quiver can be infinite (kA₂
already is — Example 2.45) and transitivity is proven only for local / hereditary /
canonical (AI Thm 1.2). `bounded_silting_exploration(A_or_T, radius, budget)` does a
**bounded-radius** BFS from the regular silting object `A` (or a supplied silting `T`)
through irreducible left+right mutations, dedup by a **shift-insensitive silting key**,
and STOPS with a loud `status`:
- `"complete"` — the frontier closed **and** the algebra is in a **proven-finite
  class** (see below), so the discovered set is certified to be all silting objects
  (up to shift);
- `"radius"` — the radius bound was reached with the frontier still open (honest
  truncation — the reported set is a ball, not the whole quiver);
- `"budget"` — the vertex budget was hit first (loud cap).

**Certified-complete only where finiteness is proven — LOCAL ONLY (Ruling 7):**
- **local** (`#simples == 1`): `silt = {A[i]}` (Thm 2.26) — mod shift a single vertex,
  so radius-0 already closes → `"complete"`. This is the **only** class the bounded
  exploration certifies complete.
- **The 2-term slice is NOT an exploration mode.** `bounded_silting_exploration` walks
  the FULL silting quiver via `silting_neighbors` (general irreducible mutations), and a
  general silting mutation **leaves the 2-term slice** — so there is no way to restrict
  this walk to 2-term objects by mutate+filter (a filtered walk would miss 2-term
  objects reached only through non-2-term intermediates, and there is no `scope="two_term"`
  parameter). The finite 2-term slice is instead a **direct cross-check against P45's
  `exchange_graph`** (Task 4a) — enumerate the support τ-tilting pairs, translate each to
  its 2-term silting object, verify silting per object, and tie the **count** to
  `len(exchange_graph.vertices)`. That crosscheck lives entirely in Task 4, NOT in this
  exploration, and matches the record's "2-term slice = cross-check only, not a
  deliverable". So `finite_class` is `"local"` or `None` — there is no `"two_term"`.
- **derived-discrete**: **NO recognizer shipped** → **out of scope**, stated. The
  exploration on such an algebra returns `"radius"`/`"budget"` honestly.
- everything else (general, incl. kA₂ which is hereditary-but-infinite): the full
  silting quiver is infinite or unproven-finite → `"radius"`/`"budget"`.

**Interfaces:**
- Produces:
  ```python
  @dataclass
  class SiltingExploration:
      vertices    # list of dicts: {"summands": [<ChainComplex>], "key": <silting key>,
                  #                  "is_initial": bool, "silting": <True|"unknown">}
      arrows      # dict {(i, j): {"direction": "left"|"right", "summand": <index>}}
      status      # "complete" | "radius" | "budget"
      radius      # the radius actually reached
      finite_class  # "local" | None  -- why (or why not) it is complete (local is the
                    # ONLY certified-complete class; there is NO "two_term" -- see Ruling 7)
  def bounded_silting_exploration(start, radius=3, budget=256) -> SiltingExploration
      # `start` is an Algebra (=> begin at the regular silting A) or a silting summand
      # list. BFS via silting_neighbors (left+right), dedup by the shift-insensitive
      # key, LOUD status. NEVER claims completeness outside a proven-finite class.
  ```
  `Algebra.silting_report()` = `is_silting_object([stalk(P_v) for v])`;
  `Algebra.silting_exploration(radius, budget)` = `bounded_silting_exploration(A, ...)`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_silting_exploration.py
"""Bounded-radius silting exploration with LOUD truncation (Plan 67 / AI Thm 1.2 --
no general BFS). Literature: k[x]/(x^2) is local -> complete at radius 0 (silt = {A[i]},
mod shift a single vertex). Honest: kA_2 is INFINITE -> the exploration truncates with
status 'radius'/'budget', never 'complete'. Self-cert: every discovered vertex is a
silting object and every arrow is a single mutation."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.derived.silting import bounded_silting_exploration

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


@lit
def test_local_is_complete_mod_shift():
    from quiverlab.families import truncated_polynomial
    A = truncated_polynomial(2, field=GF(32003))     # local: silt = {A[i]}
    ex = bounded_silting_exploration(A, radius=2)
    assert ex.status == "complete" and ex.finite_class == "local"
    assert len(ex.vertices) == 1                      # mod shift: one silting object


@lit
def test_hereditary_kA2_is_infinite_truncates_loudly():
    A = linear_path_algebra(2, field=QQ)             # kA2: silting quiver INFINITE
    ex = bounded_silting_exploration(A, radius=2, budget=64)
    assert ex.status in ("radius", "budget")         # never "complete"
    assert ex.finite_class is None
    # the discovered ball is genuine: A itself is present, every vertex is silting
    assert any(v["is_initial"] for v in ex.vertices)
    from quiverlab.derived.silting import is_silting_object
    for v in ex.vertices:
        assert is_silting_object(v["summands"]).is_silting in (True, "unknown")


@selfcert
def test_budget_trips_before_radius():
    A = linear_path_algebra(3, field=QQ)
    ex = bounded_silting_exploration(A, radius=99, budget=8)
    assert ex.status == "budget" and len(ex.vertices) <= 8


@selfcert
def test_every_arrow_is_a_single_mutation():
    A = linear_path_algebra(2, field=QQ)
    ex = bounded_silting_exploration(A, radius=2, budget=64)
    for (i, j), lab in ex.arrows.items():
        assert lab["direction"] in ("left", "right")
        assert ex.vertices[i]["key"] != ex.vertices[j]["key"]
```

- [ ] **Step 2: Run to verify failure** — `ImportError: bounded_silting_exploration`.

- [ ] **Step 3: Implement** the bounded BFS + the `finite_class` decision + the
  delegates.

**Adjust to reality (Task 3):**
- **The shift-insensitive silting key.** Two silting objects `T`, `T[i]` are the same
  vertex of the reduced silting quiver (AI's "identifying `T` with `T[i]`", Example
  2.45). Key each summand-list by: sort the per-summand K₀ classes, then **normalise
  the common shift** (subtract the minimal degree across all summands so the object is
  anchored) and take the `frozenset` of per-summand (degree-shifted) dim-vector
  tuples. This dedups `A` and `A[3]`. Document that the key is shift-normalised.
- **`finite_class`:** `"local"` when `#simples == 1` (Thm 2.26 → return radius-0 with
  the single vertex, `status="complete"`) — **the only certified-complete class**.
  Otherwise `None` and the status is `"radius"`/`"budget"`. **Never** set `"complete"`
  outside the local class. (There is deliberately **no** `"two_term"` finite_class: the
  general mutation walk cannot be restricted to the 2-term slice — Ruling 7. The 2-term
  finiteness is P45's, cross-checked directly in Task 4a, not certified by this walk.)
- **Budget vs radius.** Trip `"budget"` on `len(vertices) >= budget` BEFORE appending
  the (budget+1)-th vertex (the P41 `ARQuiver` loud-cap contract); trip `"radius"` when
  the frontier still has unexpanded vertices at depth `radius`. If BOTH the frontier
  empties AND no finite-class certificate applies, the status is **still `"radius"`**
  (we reached a closed ball within radius but cannot certify it is the WHOLE quiver) —
  the honest wording: a closed ball is not a proof of completeness without a
  finiteness theorem. Only `finite_class == "local"` upgrades to `"complete"`.
- **Dedup + left+right.** Explore both left and right irreducible mutations
  (`silting_neighbors(..., "left")` and `"right"`); an undirected edge may be found
  from both ends — key the arrow on `(min, max)` vertex indices to avoid double count,
  storing the direction that produced it.
- `Algebra.silting_report` / `Algebra.silting_exploration` are lazy-import delegates
  (avoid a `derived → core` import cycle), mirroring `Algebra.exchange_graph` (P45).

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/derived/silting.py src/quiverlab/core/algebra.py \
        tests/modules/test_silting_exploration.py
git commit -m "feat(derived): bounded-radius silting exploration -- loud truncation, certified-complete only local (NO general BFS; 2-term is a direct P45 cross-check, not a walk mode)"
```

---

### Task 4: oracle batteries — 2-term ≡ P45 τ-tilting, End(μT) ≡ Oppermann, kA₂ quiver, QPA guard

**Files:**
- Test: `tests/modules/test_silting_oracles.py` (deep),
  `tests/qpa/test_silting_qpa.py` (qpa — the honest no-surface guard)
- No `src/` change beyond what Tasks 1–3 shipped (this task is oracles).

**(a) The 2-term slice ≡ P45 τ-tilting (cross-check ONLY, not a deliverable — record
mandate; done DIRECTLY against P45, NOT via the silting exploration — Rulings 5 & 7).**
P45 already builds 2-term silting objects from support τ-tilting pairs
(`tautilting.silting.two_term_silting(pair)`, returns `{"summands": [{"P1":[..],"P0":[..],
"d1":[[..]]}, ...], "complex": ChainComplex|None}`, AIR Thm 3.2). This battery certifies
the two surfaces AGREE, **per object over the SUMMAND LIST** (not the direct-summed
`["complex"]` blob):
- enumerate the support τ-tilting pairs as objects: `pairs = [rec["pair"] for rec in
  A.exchange_graph(budget_pairs=64).vertices]` (`exchange_graph.is_complete` is asserted
  for kA₂ — τ-tilting-finite);
- for each `pair`, build the **LIST** of per-summand perfect complexes of its 2-term
  silting object `T(M,P) = M ⊕ P[1]` — each module summand `M_i` as its 2-term
  presentation `[P₁ → P₀]` (via `two_term_silting_from_presentation(M_i)` / a
  `_pair_to_summand_complexes(pair)` translation over `pair.summands` +
  `pair.support`, the killed projectives `P_v[1]` as `stalk(P_v).shift(1)`), then call
  `is_silting_object(summand_list)` — and assert `is_silting is True` with
  `generation_certified_by.startswith(("tilting", "2-term"))` (the regular pair `(A,0)`
  is `A` itself → labeled `"tilting"`, a special 2-term silting; every non-regular pair
  → `"2-term"`). **Do NOT** pass `["complex"]` (a single object) to `is_silting_object`
  — with one summand the g-matrix is `1×n`, non-square, so `k0_basis` is vacuously
  `False` and the check is a no-op that asserts nothing (the original sketch's defect);
- the **count tie**: the number of these objects `== len(A.exchange_graph(...).vertices)`
  — for kA₂, `5` (P45's support-τ-tilting count; the AIR four-way identity's silting
  leg). This count comes straight from `exchange_graph`, NOT from a 2-term-restricted
  silting walk (which cannot exist — Ruling 7). (`oracle_crossengine`.)

**(b) `End(μT)` quiver ≡ Oppermann's rule (Thm 1.1).** For a concrete mutation,
compute `End_{D^b}(μ_X T)` as an `Algebra` via P43's `end_algebra_of_complex`, take its
**Ext-quiver** (P27 `Algebra.ext_algebra` / the Gabriel quiver of `End`), and compare
the **underlying quiver** to Oppermann's mutated quiver. Worked instance: the
`1 ⇄ 2, ab=ba=0` algebra (Example 2.47) — `End(P₁ ⊕ P₂) ≅ A`; after left mutation at
`P₁`, `End(P₂ ⊕ cone(P₁→P₂))` has the Oppermann-mutated quiver at vertex 1. **# PIN**
the full graded dg rule (degree bookkeeping out of scope); the battery verifies the
**degree-0 underlying quiver** matches (arrows reversed at the mutated vertex +
composites), `oracle_literature` + `oracle_crossengine`.

**(c) kA₂ silting quiver (Example 2.45).** From `A = P₁ ⊕ P₂`, verify the first ring of
irreducible mutations reproduces the AI-drawn neighbours (the two left mutations and
their re-verification as silting), `oracle_literature`. Since kA₂ silting is infinite,
pin only the **local** structure AI draws (n=2 neighbours per vertex, the concrete
cones), not a global count.

**(d) QPA honest-scope guard.** QPA 1.37 has **no silting surface** (probe live with
`NamesGVars()` / `IsBoundGlobal` for `SiltingMutation`/`SiltingObjects`/`SiltingQuiver`
— the Plan-35 precedent). `tests/qpa/test_silting_qpa.py` is an honest **skip that
FAILS if a silting surface ever appears** in QPA (so the honest-scope entry stays
truthful).

- [ ] **Step 1: Write the batteries**

```python
# tests/modules/test_silting_oracles.py
"""Silting oracles (Plan 67). Cross-engine: the 2-term slice == P45 tau-tilting (count
+ per-object silting verdict); End(mu T) underlying quiver == Oppermann's rule.
Literature: kA2 first mutation ring (AI Example 2.45); Example 2.47 cones."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.derived.silting import is_silting_object, silting_mutate

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _pair_to_summand_complexes(pair):
    """The 2-term silting object T(M,P) = M (+) P[1] of a support tau-tilting pair, as a
    LIST of per-summand perfect complexes (NOT the direct-summed blob): each module
    summand M_i -> its 2-term presentation [P1 -> P0] (degrees 1,0); each killed
    projective P_v -> the stalk P_v[1]. This is the list is_silting_object consumes."""
    from quiverlab.modules.complexes import ChainComplex
    from quiverlab.derived.tilting import two_term_silting_from_presentation
    A = pair.algebra
    out = []
    for Mi in pair.summands:
        cx, _rep = two_term_silting_from_presentation(Mi)   # [P1 -> P0], perfect
        out.append(cx)
    for v in pair.support:                                   # killed projective -> P_v[1]
        out.append(ChainComplex.stalk(A.projective(v), 0).shift(1))
    return out


@xeng
def test_two_term_silting_matches_p45():
    # Per-object, over the SUMMAND LIST (Ruling 5): every support tau-tilting pair's
    # 2-term silting object verifies silting via is_silting_object; the count ties to
    # #support tau-tilting pairs (exchange_graph) -- for kA2, 5 (AIR four-way identity).
    A = linear_path_algebra(2, field=QQ)
    eg = A.exchange_graph(budget_pairs=64)
    assert eg.is_complete
    pairs = [rec["pair"] for rec in eg.vertices]
    for pair in pairs:
        summand_list = _pair_to_summand_complexes(pair)      # a LIST, not ["complex"]
        assert len(summand_list) == len(list(A.quiver.vertices))   # n summands = n simples
        rep = is_silting_object(summand_list)
        assert rep.is_silting is True                        # NOT a vacuous no-op
        assert rep.generation_certified_by.startswith(("tilting", "2-term"))
    assert len(pairs) == 5                                    # kA2 count tie (P45)


@xeng
def test_end_of_mutant_quiver_matches_oppermann():
    # 1 <=> 2, ab = ba = 0 (Example 2.47): End(P2 (+) cone(P1->P2)) has the vertex-1
    # Oppermann-mutated quiver. Compare the UNDERLYING quiver (degree-0), # PIN the
    # graded rule. Uses P43 end_algebra_of_complex + P27 ext-quiver.
    from quiverlab.modules.complexes import ChainComplex
    from quiverlab.derived.tilting import end_algebra_of_complex
    A = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    mut = silting_mutate(T, 0, direction="left")
    E = end_algebra_of_complex(mut)                   # the derived-equivalent algebra
    # its Gabriel/Ext quiver arrow-set matches the Oppermann prediction (transcribed
    # in the test as the expected (source,target) multiset).
    ...


@lit
def test_kA2_first_mutation_ring():
    A = linear_path_algebra(2, field=QQ)
    from quiverlab.modules.complexes import ChainComplex
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    nbrs = [silting_mutate(T, i, "left") for i in range(2)]
    for mut in nbrs:
        assert is_silting_object(mut).is_silting in (True, "unknown")
    assert {tuple(sorted(str(c.degrees()) for c in m)) for m in nbrs}  # two distinct
```

```python
# tests/qpa/test_silting_qpa.py
"""QPA has NO silting surface (Plan 67 honest scope). This test SKIPS, and FAILS if a
silting surface ever appears in QPA (so the verification-page honest-scope entry stays
truthful) -- the Plan-35 no-surface precedent."""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_has_no_silting_surface():
    lg = session.libgap_handle()
    names = ("SiltingObjects", "SiltingMutation", "SiltingQuiver", "SiltingComplexes")
    present = {nm: bool(lg.eval(f'IsBoundGlobal("{nm}")')) for nm in names}
    assert not any(present.values()), (
        f"QPA now exposes a silting surface {present} -- wire a live crosscheck and "
        "update the verification page's honest-scope entry")
```

- [ ] **Step 2: Fill the `...`** — transcribe the Oppermann-predicted quiver arrow
  multiset for the Example-2.47 mutation (derive by hand from Thm 1.1's rule at vertex
  1; **# PIN** the graded degrees, verify the underlying `(source, target)` multiset).
  (The 2-term count tie in `test_two_term_silting_matches_p45` comes DIRECTLY from
  `exchange_graph` — there is no `two_term`-scoped exploration mode to wire, Ruling 7.)
- [ ] **Step 3: Run live** `... -m pytest tests/modules/test_silting_oracles.py -q`
  (deep) and `... -m pytest tests/qpa/test_silting_qpa.py -q -m qpa` (the venv has the
  `[qpa]` extra). Expected: PASS (QPA test = honest skip that would FAIL on a new
  surface).
- [ ] **Step 4: Commit**

```bash
git add tests/modules/test_silting_oracles.py tests/qpa/test_silting_qpa.py
git commit -m "test(derived,qpa): silting oracles -- 2-term == P45, End(muT) == Oppermann, kA2 ring; QPA no-surface guard"
```

---

### Task 5: GUI/webapp/HPC — the `silting` compute kind

One **algebra-level** compute kind (schema v1 — no new request block; sizes on
`A.dim` like the HH/`tau_tilting`/`derived_fingerprint` kinds), wired through the
standard **seven-touchpoint** pattern (the P43 `derived_fingerprint` and P45
`tau_tilting` precedents).

**Files:**
- Create/modify: `src/quiverlab/derived/block.py` (add `silting_block(A, radius,
  budget)`)
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch` scalar branch + `_snip` recipe)
- Modify: `docs/gui/runner.py` (the byte-identical Pyodide twin handler + ETA)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox id
  `qlgui-silting` + a radius/budget picker, `S.ids`, push list, block renderer;
  layout-picker row for the exploration graph — the P45 `tau_tilting` layout precedent)
- Modify: `webapp/static/app.js` (webapp block renderer)
- Modify: `webapp/server/i18n/en.json`, `es.json` (labels — EN + ES; ×2 renderers ⇒
  ×4 surfaces)
- Modify: `src/quiverlab/trace/results_html.py` (report renderer branch)
- Modify: `tests/webapp/_runner_goldens.json` + `test_runner_delegation.py` (ONE new
  fixture, existing entries byte-identical)
- Test: `tests/webapp/test_silting_p67.py`, `tests/gui/test_silting_runner_twin.py`

**The block shape:**
```python
{"kind": "silting",
 "n": int,                                # |Q_0|
 "regular": {"is_silting": True, "generation_certified_by": "tilting (AI Ex 2.2)",
             "window": [1, 0], "g_matrix": [[...]], "det": 1},   # the verifier on A
 "neighbors": [ {"summand": i, "direction": "left",
                 "is_silting": True|"unknown", "summand_dimvecs": [...]}, ... ],
 "exploration": {"status": "complete"|"radius"|"budget", "radius": r,
                 "finite_class": "local"|None,   # local is the only complete class
                 "vertices": [ {"key": "...", "is_initial": bool}, ... ],
                 "arrows": [ {"from": i, "to": j, "direction": "left"|"right"}, ... ]},
 "co_t_structure": {<co_t_structure_of(A) record>},
 "scope": "<the honest three-valued + no-general-BFS + infinite-quiver note>",
 "references": ["aihara_iyama_silting", "oppermann_silting_quivers",
                "jorgensen_cotstructures"], "citations": [...]}
```

- [ ] **Step 1: Write the failing cross-runner test** (unmarked, extras-gated dir —
  copy `tests/webapp/test_tau_tilting_p45.py`'s runner-pair fixture)

```python
# tests/webapp/test_silting_p67.py
"""The silting algebra-level compute kind (Plan 67): served by hpc.spec, mirrored
byte-identically by the Pyodide twin. Local k[x]/(x^2) -> complete; kA2 -> honest
truncation; the regular verdict is silting."""
import json


def test_silting_block_local_complete(tmp_path):
    # request: k[x]/(x^2) via the standard algebra request dict; compute ["silting"].
    # Assert: block["regular"]["is_silting"] is True;
    #         block["exploration"]["status"] == "complete";
    #         block["exploration"]["finite_class"] == "local";
    #         "aihara_iyama_silting" in [k for k,_ in block["citations"]].
    ...


def test_silting_block_kA2_truncates(tmp_path):
    # kA2 -> block["exploration"]["status"] in ("radius","budget"),
    #        block["exploration"]["finite_class"] is None; no crash.
    ...


def test_twin_parity(tmp_path):
    # run the SAME request through docs/gui/runner.py; assert json.dumps(sort_keys=True)
    # equality on the silting block (both runners byte-identical), the way
    # tests/gui/test_tau_tilting_runner_twin.py does.
    ...
```

- [ ] **Step 2: Implement** `silting_block(A, radius=3, budget=64)` in
  `derived/block.py` (verifier on the regular `A` + `silting_neighbors` + a
  `bounded_silting_exploration` + `co_t_structure_of`, with the honest `scope`
  string), then the `spec.py::_dispatch` scalar branch (parse `silting:radius,budget`
  like `tau_tilting:budget` parses its suffix; catch the `decompose` char-caveat
  `QuiverlabError` into `{"error": <loud message>}` per the Plan-30 honest-per-entry
  precedent — never a 500), the byte-identical `docs/gui/runner.py` twin, the two
  `gui.js` renderers + `app.js` (checkbox + radius/budget picker + the exploration
  graph fieldset with a layout picker, mirroring `tau_tilting`), the ETA entry
  (`"silting": ~1.5` — the bounded BFS of hyper-Hom dominates), i18n keys
  (`inv.silting`, `block.silting.title`, `block.silting.scope`,
  `block.silting.status_*`, field labels — EN and ES), the `_snip` recipe
  (`"silting": lambda it: f"A.silting_exploration(radius={it.radius}, budget={it.budget})"`),
  and the `results_html.py` branch (verifier verdict + neighbours table + exploration
  edges + co-t-structure note).

- [ ] **Step 3: Add ONE golden fixture** (`silting_local` — the `k[x]/(x²)` complete
  run) to `_runner_goldens.json`; note it in `test_runner_delegation.py`'s docstring
  change-log (the `products_loop_gf2` precedent — new fixture, existing entries
  byte-identical). Run the delegation test BEFORE adding to confirm existing entries
  are untouched. **Canonical-key note:** `silting` is a schema-v1 scalar kind with a
  `radius,budget` suffix — it keys through the frozen `_V` constant exactly like
  `tau_tilting:budget`; no re-keying (metaplan §3).

- [ ] **Step 4: Run the gates**

Run: `... -m pytest tests/webapp/test_silting_p67.py tests/webapp/test_runner_delegation.py tests/gui/test_silting_runner_twin.py tests/hpc -q`
Expected: PASS (both runners byte-identical; the kA₂ case returns a truncated status
cleanly; the local case is `"complete"`).

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc,trace): silting compute kind -- verifier verdict + mutation neighbours + bounded exploration + co-t-structure, both runners + i18n x4 + one golden"
```

---

### Task 6: citations, verification page, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick the P67 card in §6)
- Test: existing release gates (`tests/release/test_oracle_classes.py`,
  `tests/citations/test_bib_structure.py`) + one citation-presence assertion.

- [ ] **Step 1: Citations** (VERIFIED BibTeX only — the `_r(...)` registry precedent).
  Add to `references.bib`:

```bibtex
@article{AiharaIyama2012,
  author  = {Aihara, Takuma and Iyama, Osamu},
  title   = {Silting mutation in triangulated categories},
  journal = {Journal of the London Mathematical Society},
  volume  = {85},
  number  = {3},
  pages   = {633--668},
  year    = {2012},
}
@article{Oppermann2017,
  author  = {Oppermann, Steffen},
  title   = {Quivers for silting mutation},
  journal = {Advances in Mathematics},
  volume  = {307},
  pages   = {684--714},
  year    = {2017},
}
@misc{Jorgensen2016cotstructures,
  author  = {J{\o}rgensen, Peter},
  title   = {Co-t-structures: the first decade},
  note    = {arXiv:1603.09379},
  year    = {2016},
}
```

  and in `registry.py` (mirroring `_r("air_tau_tilting", "AIR2014", ...)`):

```python
_r("aihara_iyama_silting", "AiharaIyama2012", "foundation",
   "Silting mutation in triangulated categories",
   "Aihara-Iyama: silting/presilting objects (Hom(T,T[>0])=0 + generation), silting "
   "mutation via one approximation triangle, the silting quiver = Hasse quiver, and "
   "transitivity for local/hereditary/canonical -- the ground truth for Plan 67.",
   "article"),
_r("oppermann_silting_quivers", "Oppermann2017", "foundation",
   "Quivers for silting mutation",
   "Oppermann: the quiver of the derived endomorphism ring of a left/right silting "
   "mutation -- the End(muT) quiver-mutation rule (Plan 67 oracle).", "article"),
_r("jorgensen_cotstructures", "Jorgensen2016cotstructures", "foundation",
   "Co-t-structures: the first decade",
   "Jorgensen's survey: bounded co-t-structures <-> silting subcategories via "
   "coheart = add(silting) -- the co-t-structure dictionary reference (Plan 67).",
   "misc"),
```

  **BibTeX-verify at merge:** the Jørgensen venue (arXiv `@misc` is the durable
  reference; upgrade to the Abel Symposia proceedings entry if verifiable —
  `tests/citations/test_bib_structure.py` must stay green). Optionally add
  `keller_vossieck_aisles → KellerVossieck1988` (silting origin) if BibTeX-verifiable;
  else AI is the primary attribution.

- [ ] **Step 2: Verification page.** Add the Plan-67 subsystem row (`derived/silting.py`)
  — AND note the Task-0 fix to `derived/tilting.py` (the P43 g-matrix basis bug) with
  its regression oracle — with the oracles and the honest-scope entries:
  - `oracle_selfcert`: g_proj K₀-basis (`det = 1` for the regular object over EVERY
    algebra, incl. non-unimodular-Cartan Example 2.47 / kZ3/J2 — the Task-0 regression
    oracle); presilting positive-window scan; the silting-vs-tilting separation
    (`P₂ ⊕ P₁[1]` presilting + silting but NOT tilting; `P₁` presilting but not silting);
    every mutant re-verified silting; the triangle `d²=0`; shares-`n−1`; `μ^-∘μ^+ = id`
    involution; `co_t_structure_of` refuses on non-silting input.
  - `oracle_crossengine`: 2-term slice ≡ P45 τ-tilting — **per support τ-tilting pair,
    the 2-term silting object (built from `two_term_silting` summands as a LIST) verifies
    silting via `is_silting_object`, and the count `== len(exchange_graph.vertices)`**
    (kA₂ = 5); `End(μT)` underlying quiver ≡ Oppermann's rule.
  - `oracle_literature`: `k[x]/(x²)` silting = shifts (Thm 2.26); kA₂ first mutation
    ring (Example 2.45); Example 2.47 cones; the local-complete exploration.
  - **honest-scope entries (binding):** (i) the silting quiver can be **infinite**
    (kA₂) and transitivity is proven only for local/hereditary/canonical (AI Thm 1.2)
    — **no general BFS/enumeration**; bounded-radius with loud `status`; the exploration
    certifies **complete only for local** (`finite_class ∈ {"local", None}`; there is no
    `"two_term"` — a general mutation walk cannot be restricted to the 2-term slice, so
    the 2-term finiteness is P45's, cross-checked directly, not claimed by the walk). (ii)
    generation is **three-valued** — certified `True` only on tilting/2-term/local,
    else `"unknown"` (`det(g_proj) = ±1` in the **projective** K₀ basis is necessary but
    not sufficient — thick subcategories are not K₀-classified; Krah phantom). (iii) the
    co-t-structure is a **documentation record**, not a computed subcategory. (iv)
    `End(μT)` vs Oppermann is at the **underlying-quiver** level; the full graded dg rule
    is `# PIN`'d. (v) **QPA has no silting surface** — the guard test FAILS if that
    changes. (vi) **derived-discrete recognition is out of scope** (no shipped
    recognizer). Recount the class table
    (`tests/release/test_oracle_classes.py` drives the numbers — run collection, paste
    the LIVE counts, re-run to green; note "counts as of the P67 merge; sibling plans
    in flight may shift them").

- [ ] **Step 3: README.** One features line: "silting theory (Aihara–Iyama): a
  silting-object verifier in `K^b(proj A)` (presilting + honest three-valued
  generation), single silting mutation via one approximation triangle, a
  bounded-radius exploration with loud truncation, and the co-t-structure dictionary —
  no-code in the browser."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/modules/test_silting*.py -q` (deep, the touched files),
  `... -m pytest -q -m fast`, `... -m pytest tests/webapp tests/gui -q`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q`,
  and a citation-presence check (`aihara_iyama_silting` / `oppermann_silting_quivers`
  / `jorgensen_cotstructures` resolve; the `silting` block carries them) — all green.
  Also `node --check` the touched JS + `tests/webapp/test_js_parses.py` (merge-train
  lesson).

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-67 silting oracle rows + AI/Oppermann/Jorgensen citations + honest scope (no general BFS, three-valued generation, no-QPA) + recounted classes"
```

---

## Acceptance (Plan-67 definition of done)

1. `is_silting_object`/`SiltingReport`, `silting_mutate`/`silting_neighbors`,
   `bounded_silting_exploration`/`SiltingExploration`, `co_t_structure_of`, and
   `silting_block` all public in `src/quiverlab/derived/silting.py` (+ the
   `derived/__init__.py` exports and the `Algebra.silting_report`/
   `Algebra.silting_exploration` delegates), each certified per instance or loudly
   refusing.
2. The verifier DECIDES **presilting** on the exact positive window
   (`Hom_{D^b}(T,T[n]) = 0`, `n > 0`), separating silting from tilting (the genuine
   `P₂ ⊕ P₁[1]` over kA₂: presilting + silting via the IJY 2-term rung, but NOT tilting
   — nonzero `Hom_{D^b}(T,T[-1])`) and separating silting from mere presilting (`P₁`
   alone: presilting but `#summands ≠ #simples`, so not silting). **Generation** is
   honest **three-valued** and computed in the **projective** K₀ basis via `g_proj`
   (Task 0 — det = 1 for the regular object over EVERY algebra, non-unimodular Cartan
   included) — `True` on tilting/2-term/local (`T = A`, `k[x]/(x²)`, 2-term ≡ P45),
   `"unknown"` on K₀-basis-only (Krah-phantom scope), `False` otherwise.
3. **Single silting mutation** `μ_X^±(T)` is built from **one approximation triangle**
   (complex-level minimal `add(T/X)`-approximation + `ChainMap.cone`); the mutant
   re-verifies as silting, shares exactly `n−1` summands, and `μ^-∘μ^+ = id`; **AI
   Example 2.47** (`X = cone(P₁→P₂)`, `Y = cone(P₂→P₁)`) is reproduced concretely.
4. **NO general BFS**: `bounded_silting_exploration` is bounded-radius with a loud
   `status` (`complete`/`radius`/`budget`) and is certified **complete only for local**
   (Thm 2.26; `finite_class ∈ {"local", None}`); kA₂ truncates loudly (infinite silting
   quiver), never `"complete"`. The 2-term τ-tilting-finite slice is a **direct P45
   cross-check** (Task 4a), NOT an exploration mode — a general mutation walk leaves the
   2-term slice, so there is no `"two_term"` finite_class and no 2-term-scoped walk.
5. The **co-t-structure dictionary** is a documentation record (coheart `add(T)` + AI
   Prop 2.23(b)/Jørgensen references), refusing on non-silting input.
6. Oracles green: the Task-0 g_proj regression (regular object is tilting/silting over
   the non-unimodular-Cartan Example 2.47 + kZ3/J2) and the verifier smoke tests;
   `k[x]/(x²)` = shifts, kA₂ first ring (Example 2.45), Example 2.47 cones (literature);
   2-term ≡ P45 τ-tilting — **per-object over the `two_term_silting` summand LIST plus
   the count `== len(exchange_graph.vertices)` (kA₂ = 5)** — and `End(μT)` underlying
   quiver ≡ Oppermann Thm 1.1 (cross-engine, graded rule `# PIN`'d); **QPA no-surface
   guard** (fails if QPA ever adds silting).
7. The `silting` compute kind clickable end-to-end (GUI canvas → block → report) in
   EN+ES, schema-v1 algebra-level scalar kind with a `radius,budget` suffix, both
   runners byte-identical, ONE golden added with a documented change-log entry.
8. **Task 0 landed FIRST**: `derived/tilting.py::g_proj` is the K₀-basis helper, both
   `is_tilting_complex` and `is_silting_object` compute generation in the projective
   basis, the non-unimodular-Cartan regression oracle is green
   (`tests/modules/test_derived_tilting.py`), and the existing P43 tilting/AR-translate
   suites (`test_derived_tilting.py` + `test_derived_tau.py`) stay green (no regression).
9. `docs/verification.md` updated (new oracle rows incl. the Task-0 g_proj fix,
   recounted classes with the merge-train note, the honest-scope entries);
   `AiharaIyama2012`/`Oppermann2017`/`Jorgensen2016cotstructures` registry keys added
   and BibTeX-verified; `tests/release/` + `tests/citations/` green; deep + qpa + fast +
   webapp/gui buckets green on the touched surface; the P67 card ticked in the
   metaplan §6 ledger.
