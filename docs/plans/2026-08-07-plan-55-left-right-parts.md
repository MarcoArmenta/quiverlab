# Plan 55: left/right parts L_A, R_A + support algebras (P55 / R15) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The tier-β **substrate** for the Assem-school recognizer ladder. Make the
**left/right parts** `L_A`, `R_A` of the module category no-code and computable:
`L_A` = the indecomposables all of whose predecessors have projective dimension ≤ 1,
`R_A` dually (successors, injective dimension ≤ 1); the finite **complement**
`ind A ∖ (L_A ∪ R_A)` (laura's cofiniteness datum — which need NOT be empty even for an ada
algebra); the **Ext-injectives of add L_A** (the τ⁻¹-criterion — P60's left-section hook) and
their dual Ext-projectives of `add R_A`; and the **left/right support algebras**
`A_λ = End(⊕ P_x : P_x ∈ L_A)` and `A_ρ` dually, each reported as a genuine presented
`Algebra` with its **connected-component factors** (the "product of tilted / quasi-tilted
algebras" — verified downstream in P60). Every part is computed by a **closed-under-
predecessors pd/id sweep on the knitted AR quiver** — **representation-finite scope only** (the
knit refuses self-injective and rep-infinite input; we inherit that loud refusal). One no-code
GUI compute kind — `left_right_parts` — puts the whole **module-category atlas** (both parts
named with `S_v/P_v/I_v`, the complement, the Ext-injectives, the two support algebras) one
click away. Every result self-certifies (predecessor closure, D-duality `D R_A = L_{A^op}`,
support-dim = `dim End`, gl.dim ≤ 1 ⇒ both parts total) or refuses **loudly**; exact only, no
floats.

**Architecture:** ONE new module, a thin exact layer over primitives that already exist
(P41 AR knitting, P05/P23 pd/id, P37 hom/morphism/End-algebra, P23/24 opposite + duality):

- **`src/quiverlab/modules/left_right.py`** (new) — the whole surface. This verb ENUMERATES
  and COMPARES whole iso-classes of modules (it consumes the P41 `knit_ar_quiver`
  indecomposable universe + the `hom_dim` predecessor relation + per-module pd/id probes),
  so it is a `modules/` citizen exactly like `modules/degeneration.py` (P49), and its
  battery lives in `tests/modules/` (the **deep** bucket, `tests/conftest.py`). Public
  surface:
  - `left_right_parts(A, *, budget=256) -> LeftRightAtlas` — the primary compute: knit,
    predecessor sweep, pd/id ≤ 1, both parts + intersection + complement + Ext-injectives +
    the two support algebras, in **one pass** over the shared knit/Hom data.
  - `@dataclass LeftRightAtlas` — mirrors P41 `ARQuiver` / P49 `DegenerationPoset` shape and
    the honest `is_complete`/`status` contract; carries the actual `Module`s (for the P60/P61
    consumers in-process) plus JSON-ready summaries. **All its fields are assembled in Task A**
    (ext-injective and support fields default to `[]`/`None` there); Tasks B and C populate
    them — so every task boundary is green (see M3 resolution below).
  - `@dataclass SupportAlgebra` — `{vertices, dim, algebra, components}`; `algebra` a genuine
    presented `Algebra` (the induced convex full subquiver), `components` its
    connected-component factors.
  - `left_right_parts_block(A, *, budget=256) -> dict` — the JSON block the two runners share
    (strips `Module`s to names + dim-vectors).
  - Thin `Algebra` delegates in `core/algebra.py` beside `ar_quiver` / `degeneration_order`:
    `Algebra.left_right_parts(budget=256)`, `Algebra.left_part(budget=256)`,
    `Algebra.right_part(budget=256)`, `Algebra.support_algebras(budget=256)`.

No new math engines; no new AR machinery. Everything is a bounded sweep over the finite
indecomposable universe the knit already produces, plus the shipped exact `hom_dim`,
`projective_resolution`, `injective_resolution`, `tau_minus`, `end_algebra`, and
`Quiver.algebra` constructors.

**Tech Stack:** exact `Domain` linear algebra throughout — `modules/hom.py::hom_dim`/`end_dim`
(exact in **every** characteristic — the predecessor relation carries no char caveat),
`modules/resolution.py::ProjectiveResolution.betti` (pd ≤ 1 ⟺ `projective_resolution(2)
.betti(2) == 0`), `modules/injective.py::injective_resolution` (id ≤ 1 ⟺
`injective_resolution(M,2).betti(2) == 0`), `modules/duality.py::tau_minus`/`dualize`,
`modules/endomorphism.py::end_algebra`, `modules/morphism.py::direct_sum`,
`combinat/quiver.py::Quiver.algebra` (the induced subquiver build),
`modules/hom.py::is_isomorphic`/`identify_standard` (**QQ / char-scope-gated** — see M1).
No floats in `src/` (AST-gated by `tests/test_no_floats.py`); all dim-vectors are
`dict[vertex,int]`, all dims `int`.

---

## Record (R15, verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R15 — Left/right parts L_A, R_A + support algebras (the substrate).**
> [C-scout P2; keep — prerequisite for R17/R18/R21]
> Object: L_A/R_A by the closed-under-predecessors pd/id sweep on the knitted AR
> quiver (rep-finite scope, loud otherwise); Ext-injectives of add L_A; left/right
> support algebras A_λ, A_ρ (theorem: products of tilted algebras — verified by R17).
> Refs: Assem–Coelho–Trepode J. Algebra 281 (2004) 518–534;
> Assem–Castonguay–Lanzilotta–Vargas arXiv:1102.1188; the "Organising the module
> category" survey (São Paulo J. Math. Sci. 2021). Size M.

Tier β, C-cluster spine; dependency chain (research doc §9): **R15 → R16 → {R17, R18,
R21, R37}**. P55 is the topmost, independent substrate of Wave 2 (metaplan §4).

---

## Reference re-verification (mandatory standing rule; done at spec time)

Re-verified via web + the arXiv PDF of 1102.1188 (read pages 1–11). **Findings:**

1. **Assem–Coelho–Trepode, J. Algebra 281 (2004) 518–534** — title **"The left and the
   right parts of a module category"**, no. 2, pp. 518–534. Vol/issue/pages/year verified.
   `# PIN`: the exact title wording ("part" vs "the right parts") differs across secondary
   listings — the worker confirms it against the J. Algebra 281(2) table of contents at
   citation-add time. This is reference `[5]` inside ACLV and the ground truth for `L_A/R_A`,
   the Ext-injective criterion, and the support algebra `A_λ`.

2. **Assem–Castonguay–Lanzilotta–Vargas, "Algebras determined by their supports",
   arXiv:1102.1188** — authors and title verified verbatim from the PDF ("IBRAHIM ASSEM,
   DIANE CASTONGUAY, MARCELO LANZILOTTA, AND ROSANA R. S. VARGAS"). Ships as an **arXiv
   preprint** BibTeX entry (every field verified from the PDF); `# PIN`: the published venue
   volume/pages is *not* verified here — the worker either verifies and upgrades the entry to
   `@article` or leaves the verified `@misc` arXiv entry (house rule: no fabricated
   bibliographic fields).

3. **"Organising the module category" survey** — authors **Alvares, E. R.; Assem, I.;
   Castonguay, D.; Vargas, R. R. S.**, São Paulo J. Math. Sci. **(2022) 16:62–82**, DOI
   `10.1007/s40863-021-00241-4`. `# PIN`: authors/pages confirmed at citation-add time.

**Definitions extracted verbatim (ACLV §1.2, §1.1, §2, and the Introduction):**

- **Paths / predecessors (§1.2).** For `M, N ∈ ind A`, a *path* `M ⤳ N` is a sequence of
  **non-zero morphisms** `M = X_0 →^{f₁} X_1 → … →^{fₜ} X_t = N` (t ≥ 1, `X_i ∈ ind A`);
  then `M` is a **predecessor** of `N`, `N` a **successor** of `M` (`M ≤ N`).
- **Left part (§1.2).** `𝓛_A = { M ∈ ind A | for any L ⤳ M, we have pd L ≤ 1 }`. "Note that
  `𝓛_A` is closed under predecessors." **Right part `𝓡_A`** dually (successors, `id ≤ 1`),
  closed under successors. (Consequently `M ∈ 𝓛_A ⟹ pd M ≤ 1`; take the trivial predecessor.)
- **Ext-injective (§1.2).** For `𝓒 = add 𝓛_A`, an indecomposable `M ∈ 𝓒` is *Ext-injective*
  in `𝓒` iff `Ext¹_A(−, M)|_𝓒 = 0`. **"M is Ext-injective in add 𝓛_A if and only if
  τ⁻¹_A M ∉ 𝓛_A"** (ACLV cites [13](3.4) = ACT 2004). Dually **`M` is Ext-projective in
  `add 𝓡_A` iff `τ_A M ∉ 𝓡_A`**. ✅ confirms the task's design question (3) exactly.
- **Support algebra (§2, and the Introduction).** "the left support `A_λ` of an artin
  algebra is the **endomorphism ring of the direct sum of all the indecomposable projective
  modules lying in the left part** of mod A, and the right support `A_ρ` is defined dually."
  ACLV p. 4: **"Let P denote the direct sum of a complete set of representatives of the
  isomorphism classes of indecomposable projective A-modules lying in `𝓛_A`. Then the
  algebra `A_λ = End P` is called the left support of A."** `A_λ` is a **full convex
  subcategory** of A (= `e_λ A e_λ`, `e_λ = Σ_{P_x ∈ 𝓛_A} e_x`), `𝓛_A ⊆ ind A_λ`, and —
  ACT [5](2.3) — **`A_λ` is a direct product of quasi-tilted algebras**; for **ada** algebras
  (ACLV Thm A) it is a direct product of **tilted** algebras. ⟵ the record's "products of
  tilted algebras — verified by R17/P60" is precisely this, with the honest general
  statement being *quasi-tilted*.
- **Duality (§2, verbatim).** `D 𝓛_A = 𝓡_{A^op}` and `D 𝓡_A = 𝓛_{A^op}`. ⟵ the self-cert.
- **Laura / complement (§3, Cor. 3.10 context, verbatim).** "an artin algebra A is *laura*
  if the class `ind A ∖ (𝓛_A ∪ 𝓡_A)` contains only finitely many objects." ⟵ P61's datum.

**Worked example (ACLV Example 2.2(b) — re-verified by hand AND computed live in the venv).**
`A = kQ/rad²` for the **linear Nakayama** `1 ← 2 ← 3 ← 4 ← 5` (arrows `α_i : i+1 → i`,
`rad²A = 0`). This is a **representation-finite ada** algebra with **9 indecomposables**,
`gl.dim A = 4` (ACLV Remark 2.7(a): the bound `gl.dim ≤ 4` is sharp here). The full atlas
(pd/id per indecomposable, recomputed with the shipped engine):

| indec | dim-vector | pd | id | part |
|-------|-----------|----|----|------|
| `S₁ = P₁`   | (1,0,0,0,0) | 0 | 4 | `𝓛_A` |
| `S₂`        | (0,1,0,0,0) | 1 | 3 | `𝓛_A` |
| `S₃`        | (0,0,1,0,0) | **2** | **2** | **complement (neither part)** |
| `S₄`        | (0,0,0,1,0) | 3 | 1 | `𝓡_A` |
| `S₅ = I₅`   | (0,0,0,0,1) | 4 | 0 | `𝓡_A` |
| `P₂ = I₁`   | (1,1,0,0,0) | 0 | 0 | `𝓛_A` |
| `P₃ = I₂`   | (0,1,1,0,0) | 0 | 0 | `𝓛_A` |
| `P₄ = I₃`   | (0,0,1,1,0) | 0 | 0 | `𝓡_A` |
| `P₅ = I₄`   | (0,0,0,1,1) | 0 | 0 | `𝓡_A` |

So `𝓛_A = {S₁, S₂, P₂, P₃}`, `𝓡_A = {S₄, S₅, P₄, P₅}`, **`𝓛_A ∩ 𝓡_A = ∅`**, and
**`complement = {S₃}` (dim-vector `{3:1}`, pd 2, id 2) — NOT empty.** The support-vertex sets
are `e_λ = {1,2,3}` (`P₁, P₂, P₃ ∈ 𝓛_A`) and `e_ρ = {3,4,5}` (`I₃=P₄, I₄=P₅, I₅=S₅ ∈ 𝓡_A`),
overlapping at vertex 3. **`A` is ada** (every `P_x` and every `I_x` lies in `𝓛_A ∪ 𝓡_A`) yet
its complement is non-empty — ada does **NOT** imply `𝓛_A ∪ 𝓡_A = ind A`. This is exactly the
interesting P61 *laura* datum (a rep-finite ada algebra with a non-empty finite complement),
and it is corroborated by our own quasi-tilted criterion: quasi-tilted would force `gl.dim ≤ 2`,
but `gl.dim = pd S₅ = 4`, so the complement CANNOT be empty. (An earlier spec draft's claim
"ada ⇒ 𝓛∪𝓡 = ind A" was a false implication — corrected here per adjudicated review H1.)

**Rep-infinite ada (ACLV Example 2.2(c), transcribed).** `1 ⇉ 2 ⇉ 3 ⇉ 4` bound by
`rad²A = 0` is a **representation-infinite** ada algebra — "an ada algebra may have infinitely
many indecomposables which are not in `𝓛_A ∪ 𝓡_A`." ⟵ the honest rep-infinite refusal oracle
(mathematically ada, but our engine refuses because rep-infinite).

---

## Consumers (design the API so P60 / P61 are served — name their points)

This plan is the **substrate** (metaplan §5 P55). Its API MUST expose, as first-class atlas
fields, exactly what the two downstream recognizers consume:

- **P60 (R17, tilted-algebra recognizer, `needs P55`)** consumes:
  1. `atlas.ext_injectives_left` — the indecomposable **Ext-injectives of add `𝓛_A`** (the
     `τ⁻¹ ∉ 𝓛_A` set). ACLV Thm 3.1 identifies these as the **left sections** `Σ′_j`; P60's
     faithful-section search seeds from them. Dual `atlas.ext_projectives_right` = the
     **right sections** `Σ` (Ext-projectives of `add 𝓡_A`, ACLV §1.3/§3.1).
  2. `atlas.left_support` / `atlas.right_support` — the support algebras **as presented
     `Algebra`s** with their component factors; P60 verifies each factor is *tilted* (the
     `# PIN` cross-check this plan leaves open).
- **P61 (R18, quasi-tilted/shod/laura/ada ladder, `needs P55`)** consumes:
  1. `atlas.complement` — `ind A ∖ (𝓛_A ∪ 𝓡_A)` (the **finite complement**); *laura* ⟺ this
     is finite (always true in our rep-finite scope; P61 reports its size), *quasi-tilted* ⟺
     it is empty and every indecomposable has pd ≤ 1 or id ≤ 1 (⟹ `gl.dim ≤ 2`). The ACLV
     2.2(b) example (complement `{S₃}`) is the ada-but-not-quasi-tilted witness.
  2. `atlas.projective_placement` / `atlas.injective_placement` — per-vertex membership maps
     `{v → "L"|"R"|"both"|"neither"}`; *ada* ⟺ every `P_v` and every `I_v` is in `𝓛_A ∪ 𝓡_A`
     (ACLV Def. 2.1) — computable from these maps **even when the complement is non-empty**.
  3. `atlas.left` / `atlas.right` union membership (`atlas.in_union(index)`).

Both consumers also read `is_complete`/`status` so their own rep-finite scope refuses in
lockstep.

---

## Mathematical foundation (derived — the plan's ground truth)

*The sweep.* Knit the AR quiver → the finite universe `U = [X₀ … X_{N−1}]` of
indecomposables (each an actual `Module`) with names `S_v/P_v/I_v|None`. Build the
**predecessor relation** `≤` on `U` and the two **Boolean sweeps** pd ≤ 1, id ≤ 1. Then:

- `bad_pd = { i : pd Xᵢ ≥ 2 }`, `bad_id = { i : id Xᵢ ≥ 2 }`.
- **`𝓛_A = { j : ∀ i with Xᵢ ≤ Xⱼ (incl. i = j), i ∉ bad_pd }`** = `U ∖` (successor-closure of
  `bad_pd`). Equivalently, `Xⱼ ∉ 𝓛_A` ⟺ some predecessor of `Xⱼ` (or `Xⱼ` itself) has pd ≥ 2.
- **`𝓡_A = { j : ∀ i with Xⱼ ≤ Xᵢ (incl. i = j), i ∉ bad_id }`** = `U ∖` (predecessor-closure
  of `bad_id`). These are automatically closed under predecessors / successors resp.
- `intersection = 𝓛_A ∩ 𝓡_A`, `complement = U ∖ (𝓛_A ∪ 𝓡_A)`.

*The predecessor relation — the decisive design choice.* The definition is the **transitive
closure of "there is a non-zero morphism between the indecomposables"**. A non-zero morphism
`Xᵢ → Xⱼ` with `i ≠ j` exists ⟺ `hom_dim(Xᵢ, Xⱼ) > 0` (and `Xᵢ ≇ Xⱼ` ⟹ any such map is a
non-iso), so `≤` is **exactly** the reflexive-transitive closure of `{(i,j) : hom_dim(Xᵢ,Xⱼ)
> 0}`. This is the **shipped default**: it is *definition-faithful* (not an AR-graph proxy),
**exact in every characteristic** (`hom_dim` = `dim ker` of the intertwiner, no char caveat,
no `is_isomorphic` search), and cannot silently under- or over-report predecessors. Cost:
the `N×N` `hom_dim` matrix over the finite universe — the same order as P49's
`degeneration_order` Hom matrix, deep-bucket and budget-bounded.

*AR-reachability as a cross-check, not a proxy (closure statement — NOT a per-pair
biconditional).* On the knitted AR quiver the arrow set (the actual irreducible maps, computed
by P41 only when the knit is complete) gives a cheap reachability relation. Every AR-arrow is a
non-zero (irreducible) map, hence a Hom-edge, so **AR-reachability ⊆ Hom-closure**. Conversely,
for representation-finite algebras `rad^∞(mod A) = 0` (Auslander), so any non-zero map
`X → Y` lies in some `rad^n(X,Y)` and therefore factors through a chain of irreducible maps
`X ⤳ Y` — hence **Hom-closure ⊆ AR-reachability**. The two **closures coincide** for rep-finite
(this is a statement about the transitive *closures*; it is NOT the false per-pair claim
"`hom(X,Y) ≠ 0 ⟺ a single AR-path with non-zero composite" — mesh-relation composites can
vanish). The agreement holds *independently of standardness* (we use the knit's actual
irreducible-map arrows, never the mesh-category combinatorics, so the non-standard char-2
categories are not a hazard). We ship this closure agreement as an `oracle_crossengine` tie (it
also validates the knit's arrow completeness); the **default relation stays the Hom closure**
(belt-and-suspenders — it never depends on the arrow set being complete).

*pd ≤ 1 and id ≤ 1 via `betti(2)`.* We need only the predicate `≤ 1`, and `pd X ≤ 1 ⟺ P₂ = 0`
in the minimal projective presentation. `ProjectiveResolution.term(n)` returns the LIST of
degree-`n` summand vertices and `betti(n) = len(term(n))` (`resolution.py:181-196`), so
`P₂ = 0 ⟺ betti(2) == 0`:
`_pd_le_1(X) := (X.projective_resolution(2).betti(2) == 0)`. Dually `id X ≤ 1 ⟺ E² = 0` in
the minimal injective coresolution, so `_id_le_1(X) := (injective_resolution(X, 2).betti(2)
== 0)` (`modules/injective.py:37,58`, `id_A M = pd_{A^op}(DM)` under the hood). At depth 2 the
resolution's own `pd()`/`injective_dimension()` are equally decisive for the `≤ 1` question;
`betti(2) == 0` is simply cleaner (a direct summand-count test, no `None` bookkeeping). A
`DepthLimitError` (a syzygy exceeding `max_term_dim`) is a **loud refusal of the whole
compute** (never a silent wrong flag) — but on rep-finite indecomposables it does not fire.

*Ext-injectives of add `𝓛_A` (ACT/[13](3.4)).* `X ∈ 𝓛_A` is Ext-injective in `add 𝓛_A` ⟺
`τ⁻¹X ∉ 𝓛_A`. Compute `Y = X.tau_minus()`; if `Y.dim == 0` (X injective) then `X` is
Ext-injective; else locate `Y` in `U` (by `is_isomorphic`, QQ-scope — see M1) and `X` is
Ext-injective ⟺ that index ∉ `𝓛_A`. Dually the Ext-projectives of `add 𝓡_A` use `τ` and `𝓡_A`.

*Support algebras.* `e_λ = { x ∈ Q₀ : P_x ∈ 𝓛_A }`, `e_ρ = { x : I_x ∈ 𝓡_A }` (locate
`A.projective(x)` / `A.injective(x)` in `U`). `A_λ = End_A(⊕_{x∈e_λ} P_x) ≅ e_λ A e_λ`, a
**full convex subcategory**; since `e_λ` is convex its presentation is the **induced full
subquiver** on `e_λ` with the restricted relations (convexity ⟹ no path between chosen
vertices leaves `e_λ`, so no relation is lost or created). Build it with `Quiver(e_λ, arrows
∩ e_λ).algebra(relations|_{e_λ})`. The value `dim_k A_λ = dim_k End_A(⊕ P_x) =
Σ_{x,y∈e_λ} hom_dim(P_x, P_y)` (additivity of `Hom` on a direct sum — this is *one* quantity,
not two independent computations). The genuine **independent certificate** is that the
**presented induced-subquiver algebra reproduces that dimension**: `B.dim == Σ hom_dim(P_x,
P_y)` (a wrong convexity/relation handling would break it), cross-checked against the shipped
`end_algebra(⊕ P_x).dim` (a second, structure-constant, implementation of `dim End`); and
`e_λ` is certified convex. The **connected components** of `A_λ`'s quiver are the "**product of
(quasi-)tilted algebras**" factors — reported here, each factor's *tiltedness* `# PIN`'d for
P60. `A_ρ` dually (`e_ρ`, injectives).

*Honest scope.* `left_right_parts` is complete **iff** `A` is representation-finite **and not
self-injective** (the P41 `knit_ar_quiver` refuses self-injective with `status="unsupported"`
and rep-infinite with `status="budget"`); we surface that status and never a partial atlas.
The support-algebra build additionally needs `A` to be **quiver-presented** (`A.quiver is not
None` and `A.projective`/`A.injective` available — the GUI/webapp always are); a
structure-constants-only `A` cannot present the induced subquiver, so the support build refuses
**loudly** (`QuiverlabError`), never fabricating a quiver.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P41 (AR knitting), P37 (categorical glue: `hom`/`morphism`/`end_algebra`), P05/P23
  (projective/injective resolutions + dimensions), P23/24 (opposite + duality) are the hard
  prerequisites** — all shipped on `dev`. This plan consumes, at the signatures verified in
  `dev`:
  - P41: `Algebra.ar_quiver(budget_modules=256, budget_dim=4096) -> ARQuiver`
    (`core/algebra.py:421` → `modules/ar.py:599`), `ARQuiver.vertices` (list of
    `{"name": str|None, "dimvec": dict, "module": Module}`), `.arrows` (`{(i,j): mult}`,
    trustworthy only when `is_complete`), `.is_complete`, `.status ∈
    {"complete","budget","unsupported","error"}`, `.note`. **Self-injective ⇒
    `status="unsupported"`, empty vertices, NO raise** (`ar.py:632`); rep-infinite ⇒
    `"budget"`.
  - P05/P23: `Module.projective_resolution(length) -> ProjectiveResolution` with
    `betti(n) = len(term(n))` and `term(n)` = the degree-`n` summand-vertex **list**
    (`resolution.py:181-196` — a list has no `.dim`, so the `≤ 1` test is `betti(2) == 0`);
    `modules/injective.py::injective_resolution(M, length) -> InjectiveResolution` with the
    same `betti(n)` (`injective.py:37,58`). **There is no `projective_dimension` method** — pd
    is `projective_resolution(k).pd()`; we use the `betti(2) == 0` predicate for the `≤ 1`
    question.
  - P37: `modules/hom.py::hom_dim(M,N)`/`end_dim(M)`/`is_isomorphic`/`identify_standard`
    (`("simple"|"projective"|"injective", v)` or `None`)/`_assert_comparable`;
    `modules/morphism.py::direct_sum(*mods) -> (D, incls, projs)`;
    `modules/endomorphism.py::end_algebra(M) -> Algebra` (structure-constant, presentation-less).
  - P23/24: `Algebra.opposite()` (`modules/opposite.py:40`, needs the quiver, raises
    otherwise), `modules/duality.py::dualize(M)` (side-aware `D`, preserves the dimension
    vector), `tau_minus(M)` / `Module.tau_minus()` (right `A`-module; `= 0` on injectives),
    `Module.tau()`.
  - Core: `combinat/quiver.py::Quiver(vertices, arrows).algebra(relations=(), field=None)`;
    `families/radical_square_zero.py::RadicalSquareZero(Q, field)` (= `kQ/rad²` — the ACLV
    2.2(b) builder); `Algebra.global_dimension() -> GlobalDimension` (`.value`,`.exact`);
    `invariants/recognizers.py::is_hereditary`.
  Branch `plan-55-left-right-parts` off `dev` (metaplan Wave 2; **independent** — no edge into
  P56–P59; P60/P61 branch after this merges).
- **Char scope is load-bearing (M1 — the identification refusal hazard).** The predecessor
  relation (`hom_dim`) and the pd/id ≤ 1 probes are exact over **every** `Domain`. But every
  step that *identifies* a module in `U` — `_index_in_U` locating `τ⁻¹X` / `P_x` / `I_x` — calls
  `is_isomorphic`, which over **char 0** is decisive (generic-rank / Noether–Deuring,
  `hom.py:171-176`) but over **large GF(p) / GF(p^n)** is **positive-only and RAISES** when it
  cannot exhibit an isomorphism (`hom.py:177-183`). So on an in-scope algebra with two
  non-isomorphic indecomposables sharing a dimension vector, `_index_in_U` would propagate a
  **loud whole-compute refusal** over large GF(p). Therefore **all identification-touching
  batteries run over QQ** (decisive negatives). The `GF(32003)` byte-parity spot-check is kept
  **only** for the `kA_n` family, where each indecomposable has a *distinct* dimension vector,
  the dim-vector prefilter separates them, and `is_isomorphic` never enters the positive-only
  branch — the plan says exactly that where it uses it. Over `char ≤ dim` the engine inherits
  the loud `QuiverlabError` — never a silent wrong part or support.
- **Honest semi-decision contract (metaplan §1.3; P41 `ARQuiver` loud-cap).** `left_right_parts`
  is complete **iff** `A` is rep-finite and not self-injective. A rep-infinite, self-injective,
  or knit-erroring `A` yields `is_complete=False` with `status ∈
  {"budget","unsupported","error"}` — never a partial atlas. `LeftRightAtlas.is_complete`/
  `status` mirror `ARQuiver` exactly.
- **No floats in `src/`.** Dim-vectors `dict[vertex,int]`; all dims/multiplicities `int`;
  memberships and closures are Boolean/index sets. The only float conversion is client-side
  (`docs/gui/gui.js`, exempt).
- **Composition is left-to-right** (`a*b` = first `a` then `b`). `_assert_comparable` guards
  every cross-module/side/algebra `hom_dim`. All refusals are `QuiverlabError`.
- **Plan-32 markers.** Predecessor-closure of the parts, intersection/complement consistency,
  convexity of `e_λ`/`e_ρ`, support presented-dim = `dim End`, Ext-injective self-consistency
  (injective ⇒ Ext-injective), gl.dim ≤ 1 ⇒ both parts total = `oracle_selfcert`; the ACLV
  2.2(b) placement + complement `{S₃}`, hereditary ⇒ parts total = `oracle_literature`; the
  D-duality `D 𝓡_A = 𝓛_{A^op}`, the Hom-closure ≡ AR-reachability predecessor tie, and the
  presented-`A_λ` dim ≡ `end_algebra` dim = `oracle_crossengine`; the QPA pointwise pd/id/τ⁻
  crosscheck = the `qpa` bucket. `oracle_*` markers are FORBIDDEN in `tests/{webapp,gui,hpc}`
  (`test_oracle_classes.py`) — engine batteries live in `tests/modules/`, the cross-runner
  twins in `tests/webapp/` unmarked.
- **Mid-merge-train counts.** v1.0.0 lands ~29 subplans in waves; absolute suite counts drift
  between authoring and merge. **Task G recounts the oracle-class table at merge time by
  running `tests/release/test_oracle_classes.py`** (paste the live numbers, never a
  guessed-at-authoring count) and claims only this plan's deltas.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table
  green) and adds its citations to `citations/references.bib` + `registry.py` (`bibtex()`
  hard-fails if the two disagree). Conventional commits; green at every commit.

**M3 resolution (no red commits).** The `LeftRightAtlas` dataclass, `left_right_parts`, and the
`Algebra` delegates are assembled in **Task A**, with `ext_injectives_left` /
`ext_projectives_right` defaulting to `[]` and `left_support` / `right_support` to `None`.
**Task B** populates the two Ext fields; **Task C** populates the two support fields. Because
the fields already exist (as safe defaults) from Task A, the Ext tests (Task B) and support
tests (Task C) reference real attributes — **every task boundary is green**, no forward
reference, no `xfail` fence.

---

### Task A: `left_right.py` — sweep + `LeftRightAtlas` assembly + delegates

**Files:**
- Create: `src/quiverlab/modules/left_right.py`
- Modify: `src/quiverlab/core/algebra.py` (thin delegates beside `ar_quiver`,
  `core/algebra.py:421`, and `degeneration_order`)
- Test: `tests/modules/test_left_right_parts.py`

**Interfaces:**
- Consumes: `Algebra.ar_quiver(budget_modules=...)`, `modules/hom.py::hom_dim`,
  `Module.projective_resolution`, `modules/injective.py::injective_resolution`.
- Produces:
  ```python
  @dataclass
  class LeftRightAtlas:
      algebra: object
      left: list           # [ {"index","name","dimvec"} , ... ]  X in L_A
      right: list
      intersection: list
      complement: list     # ind A \ (L_A u R_A)          -- P61 laura datum (need NOT be empty)
      ext_injectives_left: list = ()   # populated in Task B (default empty)
      ext_projectives_right: list = ()  # populated in Task B
      left_support: object = None       # SupportAlgebra, populated in Task C
      right_support: object = None      # populated in Task C
      projective_placement: dict = None # {v: "L"|"R"|"both"|"neither"}   -- P61 ada check
      injective_placement: dict = None
      universe_size: int = 0
      is_complete: bool = False
      status: str = "error"
      note: str = ""
      _modules: tuple = ()   # actual Module universe in index order (P60/P61 in-process)
      _leq: tuple = ()       # reachability matrix (self-cert hook)
      def in_union(self, index) -> bool: ...
  def _universe(A, budget) -> tuple      # (ar, U, names, recs); refusal via ar.is_complete
  def _leq_matrix(U) -> list[list[bool]] # DEFAULT: Hom-nonzero transitive closure
  def _ar_reachability(ar, N) -> list[list[bool]]   # the CROSS-CHECK closure
  def _pd_le_1(M) -> bool     # M.projective_resolution(2).betti(2) == 0
  def _id_le_1(M) -> bool     # injective_resolution(M, 2).betti(2) == 0
  def _left_indices(leq, pd_ok) -> set[int]
  def _right_indices(leq, id_ok) -> set[int]
  def left_right_parts(A, *, budget=256) -> LeftRightAtlas
  ```
  `Algebra.left_right_parts(budget=256)` / `.left_part(...)` / `.right_part(...)` /
  `.support_algebras(...)` delegate (lazy-import to avoid `modules → core` cycles).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_left_right_parts.py
"""Left/right parts L_A, R_A via the closed-under-predecessors pd/id sweep (Plan 55 / R15).
Literature: over a hereditary algebra (pd/id <= 1 everywhere) L_A = R_A = ind A, empty
complement; the ACLV Example 2.2(b) rad^2=0 linear Nakayama splits as L_A={S1,S2,P2,P3},
R_A={S4,S5,P4,P5}, complement={S3} (pd 2, id 2) -- ada with NON-empty complement. Self-cert:
L_A closed under predecessors, intersection/complement consistent. Cross-engine: the Hom
closure and AR-quiver reachability predecessor relations coincide (rad^infty = 0)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _radsq_nakayama_a5():
    # ACLV Example 2.2(b): 1 <- 2 <- 3 <- 4 <- 5 bound by rad^2 = 0.
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


@lit
@pytest.mark.parametrize("n, count", [(2, 3), (3, 6), (4, 10)])
def test_hereditary_both_parts_are_total(n, count):
    A = linear_path_algebra(n, field=QQ)          # kA_n hereditary, pd/id <= 1 everywhere
    atlas = left_right_parts(A)
    assert atlas.is_complete and atlas.status == "complete"
    assert atlas.universe_size == count == n * (n + 1) // 2
    assert len(atlas.left) == len(atlas.right) == count      # L_A = R_A = ind A
    assert atlas.complement == []                            # hereditary => empty complement


@lit
def test_aclv_example_22b_complement_is_S3():
    A = _radsq_nakayama_a5()
    atlas = left_right_parts(A)
    assert atlas.is_complete and atlas.universe_size == 9
    assert A.global_dimension().value == 4                   # sharp (ACLV Rem 2.7a)
    # S3 (dim-vector {3:1}) is in NEITHER part -- the ada-with-non-empty-complement datum.
    comp = [r["dimvec"] for r in atlas.complement]
    assert comp == [{1: 0, 2: 0, 3: 1, 4: 0, 5: 0}]          # exactly S3
    assert {r["name"] for r in atlas.complement} == {"S_3"}
    left_names = {r["name"] for r in atlas.left}
    right_names = {r["name"] for r in atlas.right}
    assert "S_3" not in left_names and "S_3" not in right_names   # positive: S3 in neither
    assert len(atlas.left) == 4 and len(atlas.right) == 4        # proper subsets (of 9)
    # placement (up to the projective/injective side convention -- see note):
    assert atlas.projective_placement[1] in ("L", "both")
    assert atlas.projective_placement[3] in ("L", "both")
    assert atlas.projective_placement[4] in ("R", "both")
    assert atlas.projective_placement[5] in ("R", "both")


@selfcert
def test_intersection_and_complement_consistent():
    A = _radsq_nakayama_a5()
    atlas = left_right_parts(A)
    L = {r["index"] for r in atlas.left}
    R = {r["index"] for r in atlas.right}
    assert {r["index"] for r in atlas.intersection} == (L & R)   # here == empty set
    assert {r["index"] for r in atlas.complement} == set(range(atlas.universe_size)) - (L | R)


@selfcert
def test_left_part_closed_under_predecessors():
    A = linear_path_algebra(4, field=QQ)
    atlas = left_right_parts(A)
    leq = atlas._leq
    left = {r["index"] for r in atlas.left}
    for j in left:
        for i in range(atlas.universe_size):
            if leq[i][j]:
                assert i in left


@xeng
def test_predecessor_relation_two_routes_agree():
    # Hom-nonzero transitive closure == AR-quiver reachability (rad^infty = 0, rep-finite).
    from quiverlab.modules.left_right import _ar_reachability, _leq_matrix, _universe
    A = _radsq_nakayama_a5()
    ar, U, _names, _recs = _universe(A, 256)
    assert ar.is_complete
    assert _leq_matrix(U) == _ar_reachability(ar, len(U))


@selfcert
def test_selfinjective_refused_loudly():
    from quiverlab import truncated_polynomial
    A = truncated_polynomial(3, field=QQ)                    # k[x]/(x^3), self-injective
    atlas = left_right_parts(A)
    assert atlas.is_complete is False and atlas.status == "unsupported"
    assert atlas.left == [] and atlas.right == []            # never a partial atlas


@selfcert
def test_rep_infinite_refused_loudly():
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    atlas = left_right_parts(A, budget=40)                   # 2-Kronecker, rep-infinite
    assert atlas.is_complete is False and atlas.status in ("budget", "error", "unsupported")


@selfcert
def test_delegate_matches_free_function():
    A = linear_path_algebra(2, field=QQ)
    assert [r["index"] for r in A.left_right_parts().left] == \
           [r["index"] for r in left_right_parts(A).left]
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.modules.left_right`.

  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_left_right_parts.py -v`

- [ ] **Step 3: Implement the sweep core + atlas assembly** in `left_right.py`.

```python
"""Left and right parts of the module category, and the support algebras (Plan 55 / R15,
Assem-Coelho-Trepode J. Algebra 281 (2004); Assem-Castonguay-Lanzilotta-Vargas arXiv:1102.1188).

L_A = { M in ind A : pd L <= 1 for every predecessor L of M }, closed under predecessors;
R_A dually (successors, id <= 1). A predecessor of M is the source of a path of non-zero
morphisms between indecomposables ending at M -- the reflexive-transitive closure of
"hom_dim(X, Y) > 0" over the finite (representation-finite) indecomposable universe knitted
by P41. Exact over every Domain; representation-finite scope only (the knit refuses
self-injective and rep-infinite input -- surfaced, never a partial atlas)."""
from __future__ import annotations

from dataclasses import dataclass, field

from quiverlab.errors import QuiverlabError
from quiverlab.modules.hom import hom_dim
from quiverlab.modules.injective import injective_resolution


def _universe(A, budget):
    ar = A.ar_quiver(budget_modules=budget)
    U = [v["module"] for v in ar.vertices]
    names = [v["name"] for v in ar.vertices]
    recs = [{"index": i, "name": v["name"], "dimvec": v["dimvec"]}
            for i, v in enumerate(ar.vertices)]
    return ar, U, names, recs


def _leq_matrix(U):
    n = len(U)
    leq = [[i == j for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and hom_dim(U[i], U[j]) > 0:
                leq[i][j] = True
    _transitive_closure(leq)
    return leq


def _ar_reachability(ar, n):
    leq = [[i == j for j in range(n)] for i in range(n)]
    for (i, j) in ar.arrows:
        leq[i][j] = True
    _transitive_closure(leq)
    return leq


def _transitive_closure(leq):
    n = len(leq)
    for k in range(n):
        rk = leq[k]
        for i in range(n):
            if leq[i][k]:
                ri = leq[i]
                for j in range(n):
                    if rk[j]:
                        ri[j] = True


def _pd_le_1(M):
    return M.projective_resolution(2).betti(2) == 0            # P_2 = 0  <=>  pd <= 1


def _id_le_1(M):
    return injective_resolution(M, 2).betti(2) == 0            # E^2 = 0  <=>  id <= 1


def _left_indices(leq, pd_ok):
    n = len(leq)
    return {j for j in range(n)
            if all((not leq[i][j]) or pd_ok[i] for i in range(n))}


def _right_indices(leq, id_ok):
    n = len(leq)
    return {j for j in range(n)
            if all((not leq[j][i]) or id_ok[i] for i in range(n))}


def left_right_parts(A, *, budget=256):
    ar, U, names, recs = _universe(A, budget)
    if not ar.is_complete:
        return LeftRightAtlas(A, [], [], [], [], is_complete=False,
                              status=ar.status, note=ar.note or "")
    leq = _leq_matrix(U)                                        # DEFAULT: Hom closure
    pd_ok = [_pd_le_1(M) for M in U]
    id_ok = [_id_le_1(M) for M in U]
    L = _left_indices(leq, pd_ok)
    R = _right_indices(leq, id_ok)
    sel = lambda S: [recs[i] for i in sorted(S)]
    proj_place = _placement(A, U, L, R, kind="projective")     # Task-A helper (below)
    inj_place = _placement(A, U, L, R, kind="injective")
    return LeftRightAtlas(
        A, sel(L), sel(R), sel(L & R), sel(set(range(len(U))) - (L | R)),
        projective_placement=proj_place, injective_placement=inj_place,
        universe_size=len(U), is_complete=True, status="complete",
        _modules=tuple(U), _leq=tuple(tuple(r) for r in leq))
    # ext_injectives_left / ext_projectives_right and left_support / right_support keep their
    # []/None defaults here; Task B and Task C fill them (M3: no forward reference).
```

`_placement(A, U, L, R, kind)` builds `{v: "both"|"L"|"R"|"neither"}` from the membership of
`_index_in_U(U, A.projective(v))` (resp. `A.injective(v)`) in `L`/`R`. `_index_in_U` is added
in Task B (used first there); in Task A, compute `_placement` from the projective/injective
membership via the same dim-vector-prefiltered `is_isomorphic` lookup (QQ-scope; see M1). The
`Algebra` delegates are one-liners lazy-importing `quiverlab.modules.left_right`.

**Adjust to reality (Task A):**
- **The refusal is raise-free at the knit** (`ar_quiver` returns an `ARQuiver` with
  `status ∈ {"budget","unsupported","error"}` and empty `vertices`, never raising). The
  `not ar.is_complete` branch returns the empty-list atlas with that `status`/`note`.
- **`_pd_le_1` / `_id_le_1`** — confirm `projective_resolution(2)` / `injective_resolution(M,
  2)` populate `betti(2)` (they do; `betti(n) = len(term(n))`, `resolution.py:184`,
  `injective.py:37`). A `pd = 1` module has `betti(2) == 0`; `pd >= 2` has `betti(2) > 0`.
  Let any `DepthLimitError` propagate as a loud refusal — it does not fire on rep-finite
  indecomposables.
- **The ACLV 2.2(b) side convention.** quiverlab's `projective`/`injective` and the paper's
  right-module convention may exchange `L ↔ R` / `λ ↔ ρ` (an op-duality). The test asserts the
  convention-robust facts (complement `= {S₃}`, `|𝓛_A| = |𝓡_A| = 4` of 9, the placement
  pattern); if the built orientation flips the roles, swap the `L`/`R` assertions — the paper's
  statement is convention-free up to `D`. (Verified live in the venv with arrows
  `α_i : i+1 → i`: complement `= {S₃}`, `e_λ = {1,2,3}`, `e_ρ = {3,4,5}`.)
- **`_placement`/`_index_in_U` over QQ only** (M1): `is_isomorphic` is decisive over char 0;
  the batteries use QQ. Do not call these over large GF(p) with repeated dim-vectors.

- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/left_right.py src/quiverlab/core/algebra.py \
        tests/modules/test_left_right_parts.py
git commit -m "feat(modules): left/right parts L_A, R_A -- Hom-closure predecessor sweep + betti(2)<=1 pd/id, rep-finite scope, honest refusal; ACLV 2.2(b) complement={S3}"
```

---

### Task B: Ext-injectives of `add 𝓛_A` (+ dual Ext-projectives of `add 𝓡_A`)

**Files:**
- Modify: `src/quiverlab/modules/left_right.py` (add the helpers; populate the two Ext fields
  in `left_right_parts`)
- Test: `tests/modules/test_left_right_ext.py`

**Interfaces:**
- Consumes: `Module.tau_minus()` / `Module.tau()` (right `A`-modules; `= 0` on injectives /
  projectives), `modules/hom.py::is_isomorphic` (locate `τ⁻¹X` / `τX` in `U` — QQ-scope, M1).
- Produces:
  ```python
  def _index_in_U(U, M) -> int | None    # first j with dimvec(U[j])==dimvec(M) and
                                          # is_isomorphic(U[j], M); None if absent. QQ-scope.
  def _ext_injectives_left(U, left_idx) -> set[int]
      # { i in left_idx : tau_minus(U[i]).dim == 0  OR  _index_in_U(tau_minus(U[i])) not in left_idx }
      #   (ACT/[13](3.4): X Ext-injective in add L_A  <=>  tau^{-1} X notin L_A).
  def _ext_projectives_right(U, right_idx) -> set[int]
      # dual: tau in place of tau_minus, right_idx in place of left_idx.
  ```
  `left_right_parts` now fills `ext_injectives_left` / `ext_projectives_right` (the `[]`
  defaults from Task A become populated lists).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_left_right_ext.py
"""Ext-injectives of add L_A via the tau^{-1} criterion (Plan 55 / R15, ACT [13](3.4)).
Self-cert: an injective module in L_A is always Ext-injective (tau^{-1} I = 0); every
Ext-injective lies in L_A. Literature: over kA_n, L_A = ind A, so X is Ext-injective in
add L_A iff tau^{-1}X = 0 iff X is injective (there are exactly n such). QQ-scope (M1)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
def test_ext_injectives_are_in_left():
    A = linear_path_algebra(3, field=QQ)
    atlas = left_right_parts(A)
    left_ix = {r["index"] for r in atlas.left}
    assert all(r["index"] in left_ix for r in atlas.ext_injectives_left)    # subset of L_A


@lit
def test_kA3_ext_injectives_count_is_n():
    # kA_n hereditary => L_A = ind A, so Ext-injective in add L_A iff tau^{-1}X = 0 iff X
    # injective. kA_3 has exactly 3 indecomposable injectives.
    A = linear_path_algebra(3, field=QQ)
    assert len(left_right_parts(A).ext_injectives_left) == 3


@selfcert
def test_dual_ext_projectives_are_in_right():
    A = linear_path_algebra(3, field=QQ)
    atlas = left_right_parts(A)
    right_ix = {r["index"] for r in atlas.right}
    assert all(r["index"] in right_ix for r in atlas.ext_projectives_right)
```

- [ ] **Step 2: Run to verify failure** — the fields exist as `[]` (default from Task A), so
  `test_kA3_ext_injectives_count_is_n` fails (`0 != 3`) — the red that drives Step 3. (No
  forward reference: the attribute is real, just empty — M3.)

- [ ] **Step 3: Implement** `_index_in_U`, `_ext_injectives_left`, `_ext_projectives_right`,
  and wire them into `left_right_parts` (fill the two fields with `sel(...)`). `tau_minus`/`tau`
  return a right `A`-module or the zero module (`.dim == 0`); the zero case IS the
  injective/projective edge and marks Ext-injectivity/-projectivity directly. For a non-zero
  `τ⁻¹X`, `_index_in_U` locates it (rep-finite ⟹ present); if `is_isomorphic` raises off
  char-scope, propagate the loud refusal.

**Adjust to reality (Task B):**
- **`_index_in_U` cost + M1 scope.** Filter candidates by dim-vector first (like
  `identify_standard`, `hom.py:203-207`) before the `is_isomorphic` certificate. Run only over
  QQ / char 0 (decisive); over large GF(p) with a repeated dim-vector `is_isomorphic` may raise
  — that is the honest loud refusal, not a silent skip.
- **The injective count** (`test_kA3_...`) is the hard pin; `identify_standard` may name some
  injectives `S_`/`P_` (order simple → projective → injective), so assert the *count* (3), not
  the `I_` prefix.

- [ ] **Step 4: Run Tasks A+B tests, verify pass.**
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/left_right.py tests/modules/test_left_right_ext.py
git commit -m "feat(modules): Ext-injectives of add L_A (tau^{-1} criterion) + dual Ext-projectives of add R_A -- P60 section hooks, QQ-scoped identification"
```

---

### Task C: support algebras `A_λ`, `A_ρ` — induced convex subquiver + End certificate + factors

**Files:**
- Modify: `src/quiverlab/modules/left_right.py` (add `SupportAlgebra` + `_support_algebra`;
  populate `left_support` / `right_support`)
- Test: `tests/modules/test_left_right_support.py`

**Interfaces:**
- Consumes: `modules/morphism.py::direct_sum`, `modules/hom.py::hom_dim`,
  `modules/endomorphism.py::end_algebra`, `combinat/quiver.py::Quiver.algebra`,
  `Algebra.quiver`/`Algebra.relations`, `Algebra.projective`/`Algebra.injective`.
- Produces:
  ```python
  @dataclass(frozen=True)
  class SupportAlgebra:
      vertices: tuple            # e_lambda (or e_rho) -- the chosen quiver vertices
      dim: int                   # dim_k A_lambda = sum_{x,y in e} hom_dim(P_x, P_y) = dim End(+ P_x)
      algebra: object            # the presented induced-subquiver Algebra
      components: tuple          # tuple of {"vertices": (...), "algebra": Algebra} factors
  def _support_algebra(A, U, verts, projectives=True) -> SupportAlgebra
      # A_lambda = induced full subquiver on verts (convex) with restricted relations.
      # CERTIFICATE (independent): B.dim == sum_{x,y} hom_dim(gen_x, gen_y) == end_algebra(+ gen_x).dim,
      # and verts is convex.  REFUSES loudly if A has no quiver presentation.
  def _quiver_components(verts, arrows) -> list[tuple]    # connected components (undirected BFS).
  ```
  `left_right_parts` now fills `left_support` (`e_λ`, projectives) and `right_support`
  (`e_ρ`, injectives) — the `None` defaults from Task A become `SupportAlgebra`s.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_left_right_support.py
"""Support algebras A_lambda = End(+ P_x : P_x in L_A), A_rho dually (Plan 55 / R15).
Self-cert: e_lambda is convex; dim A_lambda = dim End(+ P_x) (by additivity), and the presented
induced-subquiver algebra REPRODUCES that dimension (the independent check) -- cross-checked
against end_algebra. Literature: hereditary => A_lambda = A_rho = A (connected); ACLV 2.2(b)
=> e_lambda = {1,2,3}, e_rho = {3,4,5}. The tiltedness of each factor is PIN'd for P60."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _radsq_nakayama_a5():
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


@lit
def test_hereditary_support_is_the_whole_algebra():
    A = linear_path_algebra(4, field=QQ)               # L_A = R_A = ind A
    atlas = left_right_parts(A)
    assert set(atlas.left_support.vertices) == {1, 2, 3, 4}
    assert set(atlas.right_support.vertices) == {1, 2, 3, 4}
    assert atlas.left_support.dim == A.dim             # A_lambda = A
    assert len(atlas.left_support.components) == 1      # connected


@lit
def test_aclv_22b_support_vertices():
    atlas = left_right_parts(_radsq_nakayama_a5())
    assert set(atlas.left_support.vertices) == {1, 2, 3}       # e_lambda
    assert set(atlas.right_support.vertices) == {3, 4, 5}      # e_rho (overlap at 3)


@xeng
def test_presented_support_dim_equals_end_dim():
    # the INDEPENDENT check: the presented induced-subquiver algebra reproduces dim End.
    from quiverlab.modules.endomorphism import end_algebra
    from quiverlab.modules.morphism import direct_sum
    A = linear_path_algebra(3, field=QQ)
    atlas = left_right_parts(A)
    verts = list(atlas.left_support.vertices)
    projs = [A.projective(x) for x in verts]
    D = direct_sum(*projs)[0] if len(projs) > 1 else projs[0]
    assert atlas.left_support.algebra.dim == end_algebra(D).dim   # presented build == End
    assert atlas.left_support.dim == end_algebra(D).dim


@pytest.mark.skip(reason="PIN: each support component is tilted -- VERIFIED in P60 "
                         "(tilted recognizer). Auto-flips to a real assert when P60 lands.")
def test_support_components_are_tilted_PIN():
    atlas = left_right_parts(linear_path_algebra(3, field=QQ))
    for comp in atlas.left_support.components:
        assert comp["algebra"].is_tilted()             # P60 API
```

- [ ] **Step 2: Run to verify failure** — `AttributeError: 'NoneType' ... vertices`
  (`left_support` is `None` until this task populates it).

- [ ] **Step 3: Implement `_support_algebra`.** Refuse loudly if `A.quiver is None` (checked
  BEFORE reading vertices — no dead fallback). `verts` = sorted `e_λ` (resp. `e_ρ`); generators
  `gen_x = A.projective(x)` (resp. `A.injective(x)`). Build the induced full subquiver:
  `new_arrows = {name: (s, t) for name, (s, t) in A.quiver.arrows.items() if s in verts and t
  in verts}`, `new_rels = [r for r in A.relations if _support(r) <= set(verts)]`, then
  `B = Quiver(verts, new_arrows).algebra(relations=new_rels, field=<A's field>)`. **Certify per
  instance**: `end = Σ_{x,y∈verts} hom_dim(gen_x, gen_y)` (= `dim End`); assert
  `B.dim == end == end_algebra(⊕ gen_x).dim` and `verts` is convex (`_is_convex(A, verts)` —
  every path between two chosen vertices stays inside `verts`); raise `QuiverlabError` on any
  mismatch (a convexity or presentation bug, never silent). `components =
  _quiver_components(verts, new_arrows)`, each rebuilt as its own induced-subquiver `Algebra`.
  Wire `left_support = _support_algebra(A, U, e_λ, projectives=True)` and `right_support = ...
  e_ρ ... projectives=False` into `left_right_parts`.

**Adjust to reality (Task C):**
- **Relation filtering.** `A.relations` are relation objects (paths / linear combinations);
  `_support(r)` = the set of vertices its paths touch. By convexity of `e_λ`, a relation with
  both endpoints in `e_λ` has all intermediate vertices in `e_λ`, so the endpoint test
  suffices — but compute the full support to be safe (grep the relation representation used by
  `groebner`/`core.monomial`; adapt `_support(r)` to it). Match `Quiver.algebra`'s `field=`
  signature to `A`'s field.
- **`dim End` is ONE quantity, checked two ways** (W2): `Σ_{x,y} hom_dim(gen_x, gen_y)` IS
  `dim End(⊕ gen_x)` by additivity of `Hom` — not an independent computation. The genuine
  independent certificate is `B.dim == dim End` (the presented induced-subquiver build vs the
  `End` definition); `end_algebra(⊕ gen_x).dim` is a second, structure-constant implementation
  of `dim End` used as a cheap sanity cross-check.
- **No presentation-less fallback** (M2): the build reads `A.quiver`/`A.relations` and
  `A.projective`/`A.injective`; a structure-constants-only `A` (`A.quiver is None`) cannot
  present the induced subquiver and cannot build `P_x` by vertex, so `_support_algebra` raises
  `QuiverlabError` up front (never a fabricated quiver, never a dead code path). The GUI/webapp
  always feed a quiver-presented `A`.
- **The tiltedness `# PIN`.** P55 reports the components and the *property* ("`A_λ` is a
  product of quasi-tilted algebras in general, tilted for ada" — ACT [5](2.3) / ACLV Thm A). It
  does NOT certify tiltedness (that is P60's `tilted_check`). Ship the skipped
  `test_support_components_are_tilted_PIN` as the auto-flipping fence (the Plan-29/31 xfail-
  fence precedent — a real assert when P60's `is_tilted` lands).

- [ ] **Step 4: Run Tasks A+B+C tests together, verify pass.**

  Run: `... -m pytest tests/modules/test_left_right_parts.py tests/modules/test_left_right_ext.py tests/modules/test_left_right_support.py -v`

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/left_right.py tests/modules/test_left_right_support.py
git commit -m "feat(modules): support algebras A_lambda/A_rho as presented induced-convex-subquiver algebras + End certificate + component factors (product-of-tilted, PIN P60)"
```

---

### Task D: op/dual self-cert + gl.dim consistency + literature oracles

**Files:**
- Test: `tests/modules/test_left_right_oracles.py` (no `src/` change — pure oracle battery)

**Interfaces:**
- Consumes: `Algebra.opposite()`, `modules/duality.py::dualize`, `left_right_parts`,
  `Algebra.global_dimension()`.

- [ ] **Step 1: Write the oracle battery**

```python
# tests/modules/test_left_right_oracles.py
"""Cross-engine + literature oracles for the left/right parts (Plan 55 / R15). Cross-engine:
D R_A = L_{A^op} (ACLV, verbatim) via the shipped opposite + duality, run on the ACLV 2.2(b)
example where R_A is a PROPER subset (4 of 9) so the duality is NON-vacuous; the hereditary
case is a smoke test only (both parts total, match trivial). Literature: gl.dim <= 1 => both
parts total. Honest scope: self-injective and rep-infinite refuse loudly. QQ-scope (M1)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


def _radsq_nakayama_a5():
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


def _dvms(records):
    return sorted(tuple(sorted(r["dimvec"].items())) for r in records)


@xeng
def test_D_of_right_equals_left_of_opposite_on_proper_subset():
    # ACLV: D R_A = L_{A^op}. On the rad^2=0 A5, R_A is a PROPER subset (4 of 9), so the
    # dim-vector multisets of R_A (over A) and L_{A^op} (over A^op) coincide non-trivially.
    A = _radsq_nakayama_a5()
    atlas_A = left_right_parts(A)
    atlas_op = left_right_parts(A.opposite())
    assert len(atlas_A.right) == 4 and len(atlas_A.right) < atlas_A.universe_size  # non-vacuous
    assert _dvms(atlas_A.right) == _dvms(atlas_op.left)
    assert _dvms(atlas_A.left) == _dvms(atlas_op.right)          # D L_A = R_{A^op}


@xeng
def test_D_duality_smoke_hereditary():
    # hereditary: both parts total; the duality holds trivially -- a smoke case only.
    A = linear_path_algebra(3, field=QQ)
    assert _dvms(left_right_parts(A).right) == _dvms(left_right_parts(A.opposite()).left)


@lit
def test_gldim_le_1_forces_total_parts():
    A = linear_path_algebra(3, field=QQ)            # gl.dim = 1
    assert A.global_dimension().value <= 1
    atlas = left_right_parts(A)
    assert len(atlas.left) == len(atlas.right) == atlas.universe_size
    assert atlas.complement == []


@selfcert
def test_selfinjective_and_rep_infinite_refuse():
    from quiverlab import truncated_polynomial
    assert left_right_parts(truncated_polynomial(4, field=QQ)).status == "unsupported"
    K2 = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    assert left_right_parts(K2, budget=40).is_complete is False
```

- [ ] **Step 2: Run** — the D-duality (proper-subset) and gl.dim ties pass (verified live in
  the venv: `D R_A == L_{A^op}` and `D L_A == R_{A^op}` on the rad²=0 A₅, `|R_A| = 4` of 9).
  Add the ACT-2004 worked examples as a `# PIN` (the worker transcribes ACT 281(2)'s explicit
  `L_A`/`R_A` examples at implementation, or leaves the ACLV 2.2(b) transcription — already the
  primary literature oracle in Task A — as the anchor; both are Assem-school).

- [ ] **Step 3: Commit**

```bash
git add tests/modules/test_left_right_oracles.py
git commit -m "test(modules): left/right parts oracles -- D R_A = L_{A^op} on the proper-subset rad^2=0 A5, gl.dim consistency, honest self-injective/rep-infinite refusals"
```

---

### Task E: GUI/webapp/report — the `left_right_parts` module-category atlas kind

`left_right_parts` is an **ALGEBRA-level compute kind** (computed on the drawn `A`, like
`ar_quiver` / `tau_tilting` — routed through `_dispatch`, NOT `_dispatch_module`; **schema
v1**, no module block; it sizes on `A.dim` like `ar_quiver`). It carries a budget:
`compute.push("left_right_parts:256")`.

**Files:**
- Modify: `src/quiverlab/modules/left_right.py` (`left_right_parts_block(A, budget)`)
- Modify: `src/quiverlab/hpc/spec.py` (`parse_compute_item` budget interceptor; `_dispatch`
  branch; `_snip` recipe)
- Modify: `webapp/server/schema.py` (the parse twin — budget interceptor) +
  `webapp/server/estimator.py` (`_max_degree` exclusion tuple)
- Modify: `docs/gui/runner.py` (the Pyodide twin: `_parse_compute` interceptor, `compute_one`
  branch, `calls` snippet, `ETA_MODEL["scalars"]`)
- Modify: `docs/gui/gui.js` (+ re-copy byte-identical to `webapp/static/gui/gui.js`): the
  checkbox + budget picker, the `el` id list, the `buildRequest` push, a `renderBlock` branch
  `renderLeftRightParts(div, b)` (two named lists + the complement + the Ext-injectives + the
  two support-algebra summaries), `KIND_CTRL`, `THEMES` (`"structure"`), `SEARCH_INDEX`
- Modify: `webapp/server/i18n/{en,es,fr,zh}.json` (the `inv.left_right_parts` label ×4 + any
  `block.left_right_parts.*` strings the renderer reads)
- Modify: `src/quiverlab/trace/results_html.py` (`_HEADINGS["left_right_parts"] = "Left / right
  parts"` + a `_block_html` branch `_left_right_parts_html(b)`)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE new fixture `left_right_parts_kA3`, existing entries byte-identical; docstring bullet)
- Test: `tests/webapp/test_left_right_parts_p55.py`, `tests/gui/test_left_right_runner_twin.py`

**Block shape (returned identically by both runners — the byte-parity contract):**
```python
{"kind": "left_right_parts",
 "n": int,                                  # |Q_0|
 "complete": bool, "status": "complete"|"budget"|"unsupported"|"error",
 "universe_size": int | None,
 "left":  [ {"name": "S_1"|..., "dimvec": {...}} , ... ],
 "right": [ ... ],
 "intersection": [ ... ],
 "complement": [ ... ],                     # ind A \ (L_A u R_A)   (P61 laura datum; may be non-empty)
 "ext_injectives_left": [ ... ],            # P60 left-section hook
 "ext_projectives_right": [ ... ],          # P60 right-section hook
 "left_support":  {"vertices": [...], "dim": int, "components": [{"vertices":[...]}, ...]},
 "right_support": {...},
 "projective_placement": {str(v): "L"|"R"|"both"|"neither"},
 "injective_placement":  {str(v): "..."},
 "note": str | None,
 "references": ["act_left_right", "aclv_supports", "organising_module_category"],
 "citations": [...]}
```

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir; copy the
  `tests/webapp/test_tau_tilting_p45.py` runner-pair fixture):

```python
# tests/webapp/test_left_right_parts_p55.py
"""The left_right_parts algebra-level compute kind: served by hpc.spec, mirrored byte-for-byte
by the Pyodide twin. kA3 -> L_A = R_A = ind A (6 indecs), empty complement, both supports = A;
the rad^2=0 A5 -> complement = [S3] (the non-empty ada laura datum)."""


def test_left_right_parts_block_shape(tmp_path):
    from quiverlab.hpc.spec import ComputeRequest, run
    req = _lr_request(quiver=("kA3",), budget=256)          # helper: algebra + compute
    out = run(ComputeRequest.model_validate(req), tmp_path)
    b = out["results"]["left_right_parts"]
    assert b["complete"] and b["n"] == 3 and b["universe_size"] == 6
    assert len(b["left"]) == len(b["right"]) == 6 and b["complement"] == []
    assert set(b["left_support"]["vertices"]) == {1, 2, 3}


def test_nonempty_complement_block(tmp_path):
    # the rad^2=0 A5: complement carries exactly S3 (dim-vector {3:1}) -- the ada-with-non-empty
    # complement datum P61 consumes.
    ...


def test_twin_parity(tmp_path):
    # same request through docs/gui/runner.py; json.dumps(sort_keys=True) equality on the block
    # (both runners call the shared left_right_parts_block builder).
    ...


def test_selfinjective_status(tmp_path):
    # k[x]/(x^3) -> status "unsupported", complete False, no crash, empty lists.
    ...
```

- [ ] **Step 2: Implement `left_right_parts_block(A, budget)`** — call `left_right_parts(A,
  budget=budget)`, strip `_modules`/`_leq`, serialize the `SupportAlgebra`s to
  `{vertices, dim, components:[{vertices}]}`, stamp `references = ["act_left_right",
  "aclv_supports", "organising_module_category"]`. On the `is_complete=False` path return the
  same shape with empty lists + `status`/`note` (the honest refusal, like `ar_quiver_block`).

- [ ] **Step 3: Wire the compute kind** (the algebra-level, `ar_quiver`-style touchpoints):
  - `spec.py::parse_compute_item`: add a budget interceptor mirroring `tau_tilting`/`ar_quiver`
    — `if s == "left_right_parts" or s.startswith("left_right_parts:"): ... return
    ComputeItem("left_right_parts", None, int(b) if b else None)`.
  - `spec.py::_dispatch` (before the catch-all): `if kind == "left_right_parts": from
    quiverlab.modules.left_right import left_right_parts_block; block =
    left_right_parts_block(A, budget=item.hi if item.hi is not None else 256);
    block["citations"] = _citation_pairs(block["references"]); return block, None`. Catch the
    `is_isomorphic` / `DepthLimitError` char-scope refusal into `{"kind": "left_right_parts",
    "error": "<loud message>"}` (the Plan-30 per-entry precedent — never a 500).
  - `spec.py::_snip`: `"left_right_parts": lambda it: ("A.left_right_parts("
    f"budget={it.hi if it.hi is not None else 256})")`.
  - `webapp/server/schema.py::parse_compute_item`: the same budget interceptor (parse twin #2).
    `webapp/server/estimator.py::_max_degree`: add `"left_right_parts"` to the
    `("tau_tilting","ar_quiver")` skip tuple (its budget is not a homological degree).
  - `docs/gui/runner.py`: `_parse_compute` interceptor, `compute_one` branch mirroring
    `tau_tilting` calling the **same** `left_right_parts_block`, `calls` snippet,
    `ETA_MODEL["scalars"]["left_right_parts"]`.
  - GUI (`docs/gui/gui.js`, then re-copy to `webapp/static/gui/gui.js` — the byte-identical
    gate is `test_draw_page.py`): the checkbox + budget picker (mirror the `ar_quiver` row),
    the `el` id list entry `"left_right_parts", "left_right_parts-budget"`, the `buildRequest`
    push `if (el.left_right_parts.checked) compute.push("left_right_parts:" +
    el["left_right_parts-budget"].value)`, a `renderBlock` branch `else if (name ===
    "left_right_parts") renderLeftRightParts(div, b)`, a `renderLeftRightParts(div, b)`
    fieldset (two lists of `S_v/P_v/I_v` chips, the complement, the Ext-injective list, the two
    support-algebra summaries `A_λ: vertices {…}, k factors`; `citesLine(b)` appended
    automatically), `KIND_CTRL` `{cb, budget:true}`, `THEMES` add to `"structure"`,
    `SEARCH_INDEX` entry with multilingual keywords + an `example`.
  - i18n ×4 (`en/es/fr/zh.json`): `"inv.left_right_parts"` in ALL FOUR (the `test_i18n.py`
    key-parity gate + the no-placeholder gate) + any `block.left_right_parts.*` labels the
    renderer reads.
  - `results_html.py`: `_HEADINGS["left_right_parts"] = "Left / right parts"` + a `_block_html`
    branch `if kind == "left_right_parts": return _left_right_parts_html(b)` building two named
    lists + the complement + the support-algebra rows (reuse the `_derived_fingerprint_html`
    table idiom).

- [ ] **Step 4: Add ONE golden** `left_right_parts_kA3` to `_runner_goldens.json` (full kA3
  run: 6 left, 6 right, empty complement, both supports = all vertices). Verify the existing
  entries are byte-identical BEFORE appending (a schema-v1 budget kind adds no request field, so
  their `canonical_key`s and `result_json`s are untouched); add a docstring bullet to
  `test_runner_delegation.py` documenting the pure addition.

- [ ] **Step 5: Run the gates**

  Run: `... -m pytest tests/webapp/test_left_right_parts_p55.py tests/webapp/test_runner_delegation.py tests/gui/test_left_right_runner_twin.py tests/webapp/test_i18n.py tests/webapp/test_draw_page.py tests/hpc -q`
  Expected: PASS (both runners byte-identical; self-injective returns `status="unsupported"`
  cleanly; gui.js twins byte-identical; i18n key-parity green).

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc,trace): left_right_parts compute kind -- module-category atlas (both parts + complement + ext-injectives + supports), both runners, i18n x4, one golden"
```

---

### Task F: QPA pointwise cross-oracle (probe-first) + honest fallback

**Files:**
- Test: `tests/qpa/test_left_right_qpa.py`

**Interfaces:**
- Consumes: the QPA session (`session.should_skip_qpa()`, `session.libgap_handle()`) — mirror
  the Plan-35 products probe. **QPA has no `L_A`/`R_A`/support-algebra surface** — but its
  `ProjDimensionOfModule` / `InjDimensionOfModule` / `DTr` (used by P05/P23/P41) corroborate
  the **pointwise** ingredients (the pd ≤ 1 / id ≤ 1 flags and the `τ⁻¹` membership).

- [ ] **Step 1: Probe live QPA** (the Plan-35 fail-if-appears pattern):

```python
# tests/qpa/test_left_right_qpa.py
"""QPA cross-check for the left/right parts (Plan 55). qpa-marked: skips locally, mandatory
under QUIVERLAB_REQUIRE_QPA=1. QPA has NO left/right-part or support-algebra verb -- this test
documents that (fail-if-appears) and corroborates the POINTWISE ingredients (pd<=1, id<=1) via
QPA's ProjDimensionOfModule / InjDimensionOfModule on the same kA3 indecomposables."""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_has_no_left_right_surface():
    lg = session.libgap_handle()
    for name in ("LeftPartOfModuleCategory", "RightPartOfModuleCategory",
                 "SupportAlgebra", "LeftSupport", "RightSupport"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), (
            f"QPA now exposes {name} -- wire a real crosscheck (this assert is the trip-wire "
            "that Plan 55's QPA scope note is stale)")


def test_pd_id_flags_match_qpa_pointwise():
    # the pd<=1 / id<=1 flags that DEFINE L_A/R_A, corroborated per indecomposable by QPA.
    from quiverlab import linear_path_algebra
    from quiverlab.fields import QQ
    A = linear_path_algebra(3, field=QQ)
    atlas = A.left_right_parts()
    for M in atlas._modules:
        A.crosscheck("proj_dim", M).assert_agree()          # existing P05/P23 pd crosscheck verb
```

- [ ] **Step 2:** Confirm the crosscheck verb name/signature in `qpa/crosscheck.py` (the pd/id
  crosscheck may be `"proj_dim"`/`"inj_dim"` or ride a broader profile verb — grep before
  wiring). If no clean two-argument pd crosscheck exists, the probe
  (`test_qpa_has_no_left_right_surface`) is the deliverable on its own and the pointwise pd/id
  corroboration falls back to the P05/P23 QPA batteries — honest, never a silent skip.

- [ ] **Step 3: Run live** `... -m pytest tests/qpa/test_left_right_qpa.py -v` (the venv has
  `[qpa]`). Expected: PASS live (probe absent; pd/id agree pointwise).

- [ ] **Step 4: Commit**

```bash
git add tests/qpa/test_left_right_qpa.py
git commit -m "test(qpa): left/right parts probe -- QPA has no L_A/R_A/support verb (fail-if-appears); pd<=1/id<=1 flags corroborated pointwise"
```

---

### Task G: verification page, citations, README, metaplan, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick the P55 progress ledger)
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`)

- [ ] **Step 1: Citations** (VERIFIED BibTeX only — `_r(key, bibtex_key, kind, title,
  annotation, *tags)`, `registry.py:24`; `bibtex()` hard-fails if `.bib` and registry
  disagree). Add to `references.bib`:

```bibtex
@article{ACT2004,
  author  = {Assem, Ibrahim and Coelho, Fl{\'a}vio U. and Trepode, Sonia},
  title   = {The left and the right parts of a module category},
  journal = {Journal of Algebra},
  volume  = {281},
  number  = {2},
  pages   = {518--534},
  year    = {2004},
}
@misc{ACLV2011,
  author        = {Assem, Ibrahim and Castonguay, Diane and Lanzilotta, Marcelo and Vargas, Rosana R. S.},
  title         = {Algebras determined by their supports},
  year          = {2011},
  eprint        = {1102.1188},
  archivePrefix = {arXiv},
  primaryClass  = {math.RT},
}
@article{AACV2021,
  author  = {Alvares, Edson Ribeiro and Assem, Ibrahim and Castonguay, Diane and Vargas, Rosana R. S.},
  title   = {Organising the module category},
  journal = {S{\~a}o Paulo Journal of Mathematical Sciences},
  volume  = {16},
  number  = {1},
  pages   = {62--82},
  year    = {2022},
  doi     = {10.1007/s40863-021-00241-4},
}
```

  and in `registry.py`:

```python
_r("act_left_right", "ACT2004", "foundation",
   "The left and the right parts of a module category",
   "Assem-Coelho-Trepode: L_A / R_A via predecessor/successor closure of pd<=1 / id<=1, the "
   "Ext-injective criterion tau^{-1}X notin L_A, and the support algebra A_lambda = End of "
   "the projectives in L_A -- the ground truth for Plan 55.", "recognizer"),
_r("aclv_supports", "ACLV2011", "foundation",
   "Algebras determined by their supports",
   "Assem-Castonguay-Lanzilotta-Vargas: A_lambda / A_rho are products of tilted algebras for "
   "ada algebras (quasi-tilted in general), D L_A = R_{A^op}, and the laura complement "
   "ind A minus (L_A u R_A) -- non-empty even for ada; feeds P60 (tilted) and P61 (laura/ada).",
   "recognizer"),
_r("organising_module_category", "AACV2021", "foundation",
   "Organising the module category",
   "Alvares-Assem-Castonguay-Vargas survey of the left/right parts, supports, and the "
   "quasi-tilted/laura/ada organisation of mod A -- Plan 55's secondary reference.", "survey"),
```

  **Spec-ambiguity resolution (recorded):** BibTeX entry keys `ACT2004`/`ACLV2011`/`AACV2021`;
  snake-case registry keys `act_left_right`/`aclv_supports`/`organising_module_category`
  (house convention, `assem_book → ASS2006`). ACLV ships as the **verified arXiv `@misc`**
  (published-venue fields not verified here — worker upgrades to `@article` iff a BibTeX-
  verifiable record is confirmed at merge; no fabricated fields). ACT title and survey
  authors/pages `# PIN`'d for the merge-time confirmation.

- [ ] **Step 2: Verification page.** Add the Plan-55 subsystem rows (all NON-`qpa` for the
  engine; the pointwise pd/id crosscheck is the one `qpa` row):
  - `modules/left_right.py` (parts) — `oracle_selfcert` (predecessor-closure of `𝓛_A`/`𝓡_A`,
    intersection/complement consistency, gl.dim ≤ 1 ⇒ total); `oracle_literature` (hereditary
    parts total; **ACLV Example 2.2(b): `𝓛_A={S₁,S₂,P₂,P₃}`, `𝓡_A={S₄,S₅,P₄,P₅}`,
    complement `{S₃}`, gl.dim 4 — ada with non-empty complement**); `oracle_crossengine`
    (Hom-closure ≡ AR-reachability predecessor relation; `D 𝓡_A = 𝓛_{A^op}` on the
    proper-subset A₅).
  - `modules/left_right.py` (Ext-injectives) — `oracle_selfcert` (Ext-injectives ⊆ `𝓛_A`,
    injective ⇒ Ext-injective); `oracle_literature` (kA_n Ext-injectives count = #injectives).
  - `modules/left_right.py` (support algebras) — `oracle_crossengine` (presented-`A_λ` dim ≡
    `end_algebra` dim); `oracle_selfcert` (`e_λ`/`e_ρ` convex, presented dim = `dim End`,
    hereditary ⇒ `A_λ = A_ρ = A`).
  - Add the **honest-scope entries**: (a) **representation-finite and non-self-injective only**
    — the P41 knit refuses self-injective (`status="unsupported"`) and rep-infinite
    (`status="budget"`); ACLV Example 2.2(c) is a *mathematically ada* algebra we refuse
    because rep-infinite (the honest boundary); (b) the support-algebra build needs a **quiver
    presentation** (structure-constants-only ⇒ loud refusal, no fabricated quiver); (c)
    **identification (`is_isomorphic`) is QQ / char-0 decisive; over large GF(p)/GF(p^n) it is
    positive-only and RAISES, so `_index_in_U` may propagate a loud whole-compute refusal on
    an algebra with two non-iso indecomposables sharing a dim-vector** — batteries run over QQ;
    the `GF(32003)` parity check is used only for `kA_n` (distinct dim-vectors, no raise); (d)
    the **"product of tilted algebras"** property is reported, its per-factor tiltedness
    `# PIN`'d for **P60** (skipped auto-flipping test); (e) **QPA cannot compare** the
    left/right-part surface (fail-if-appears probe; the defining pd ≤ 1 / id ≤ 1 flags ARE
    QPA-checked pointwise). Recount the class table (`tests/release/test_oracle_classes.py`
    drives the numbers — run collection for each of `oracle_literature`/`oracle_crossengine`/
    `oracle_selfcert`/`qpa`/`m2` + the union, paste the LIVE counts into the class-count table,
    re-run to green; mid-merge-train honest — as of authoring the page reads lit 882 / xeng
    528 / selfcert 1084 / qpa 189 / m2 11 / union 2122, which WILL have drifted).

- [ ] **Step 3: README.** One features line: "left/right parts of the module category
  (Assem–Coelho–Trepode): `L_A`, `R_A` via the closed-under-predecessors pd/id sweep, the
  finite complement, the Ext-injectives of `add L_A`, and the left/right support algebras
  `A_λ`, `A_ρ` (products of tilted algebras) — the recognizer-ladder substrate, no-code in the
  browser."

- [ ] **Step 4: Metaplan.** Tick `P55` in the §6 progress ledger.

- [ ] **Step 5: Full gate:**
  `... -m pytest tests/modules/test_left_right_*.py -q` (deep — the touched files),
  `... -m pytest -q -m fast`, `... -m pytest tests/webapp tests/gui -q`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q` — all green.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-55 left/right-parts oracle rows + ACT/ACLV/survey citations + honest scope (rep-finite, non-self-injective, QQ identification, tilted PIN, no-QPA) + recounted classes"
```

---

## Acceptance (Plan-55 definition of done)

1. `left_right_parts`, `LeftRightAtlas`, `SupportAlgebra`, `left_right_parts_block` public in
   `src/quiverlab/modules/left_right.py`; the `Algebra` delegates
   (`left_right_parts`/`left_part`/`right_part`/`support_algebras`) named; every result
   self-certified or loudly refusing.
2. `𝓛_A` / `𝓡_A` are the closed-under-predecessors pd ≤ 1 / closed-under-successors id ≤ 1
   sets computed via the **Hom-nonzero transitive-closure** predecessor relation (exact in
   every characteristic; pd/id ≤ 1 via `projective_resolution(2).betti(2) == 0` /
   `injective_resolution(M,2).betti(2) == 0`), with the AR-reachability **closure** agreement
   pinned cross-engine; `𝓛_A` closed under predecessors, `𝓡_A` under successors,
   `intersection = 𝓛_A ∩ 𝓡_A`, `complement = ind A ∖ (𝓛_A ∪ 𝓡_A)`.
3. Hereditary ⇒ `𝓛_A = 𝓡_A = ind A`, empty complement (kA_n; QQ, plus GF(32003) parity for
   the distinct-dim-vector kA_n family only). **ACLV Example 2.2(b)** pinned:
   `𝓛_A = {S₁,S₂,P₂,P₃}`, `𝓡_A = {S₄,S₅,P₄,P₅}`, `𝓛_A ∩ 𝓡_A = ∅`, **`complement = {S₃}`**
   (pd 2, id 2 — the ada-with-non-empty-complement datum), `e_λ = {1,2,3}`, `e_ρ = {3,4,5}`,
   `gl.dim = 4` (up to the side convention); gl.dim ≤ 1 ⇒ both parts total.
4. Ext-injectives of `add 𝓛_A` via `τ⁻¹X ∉ 𝓛_A` (⊆ `𝓛_A`; injective ⇒ Ext-injective; kA_n ⇒
   count = #injectives), and the dual Ext-projectives of `add 𝓡_A` — the P60 section hooks;
   identification QQ-scoped.
5. Support algebras `A_λ`, `A_ρ` shipped as **presented induced-convex-subquiver `Algebra`s**,
   certified per instance (the presented build reproduces `dim End P`, cross-checked against
   `end_algebra`; `e_λ`/`e_ρ` convex; a presentation-less `A` refuses loudly), with their
   **connected-component factors** reported; the "product of tilted algebras" property
   documented and its per-factor tiltedness `# PIN`'d for P60 (auto-flipping skipped test).
6. The **honest semi-decision contract**: complete **iff** `A` is representation-finite and
   non-self-injective; self-injective (`k[x]/(xⁿ)`) ⇒ `status="unsupported"`, rep-infinite
   (2-Kronecker / ACLV 2.2(c)) ⇒ `status="budget"` — never a partial atlas; identical to P41's
   `ARQuiver` loud cap.
7. The `D 𝓡_A = 𝓛_{A^op}` self-cert holds on the **proper-subset** rad²=0 A₅ (non-vacuous:
   `|𝓡_A| = 4` of 9), with hereditary kA₃ kept only as a smoke case.
8. `left_right_parts` clickable end-to-end (GUI canvas → module-category atlas block →
   report), algebra-level compute kind (schema v1, NO module block), both runners
   byte-identical via the shared `left_right_parts_block`, EN+ES+FR+ZH i18n, one golden with a
   documented change-log entry.
9. QPA probe green (`-m qpa`): no `L_A`/`R_A`/support verb (fail-if-appears); the defining
   pd ≤ 1 / id ≤ 1 flags corroborated pointwise via QPA's proj/inj-dimension.
10. **Consumers served:** `atlas.complement` (may be non-empty for ada) +
    `atlas.projective_placement`/`injective_placement` (P61 laura/ada) and
    `atlas.ext_injectives_left`/`ext_projectives_right` + `atlas.left_support`/`right_support`
    (P60 tilted recognizer) are first-class, documented atlas fields.
11. `docs/verification.md` recounted (live numbers, mid-merge-train honest); the three
    citations added and BibTeX-verified (ACLV as the verified arXiv entry; ACT title + survey
    authors `# PIN`'d); README line + metaplan P55 tick; deep (`test_left_right_*`) + fast +
    webapp/gui + qpa + release + citations suites green. Honest scope recorded: rep-finite +
    non-self-injective only; presented supports need a quiver; QQ / char-0 for the
    identification certificates (large-GF(p) positive-only refusal named); tilted-factor claim
    `# PIN`'d for P60; QPA cannot compare the parts.

---

## Change-log

- **2026-08-07** — plan authored (P55 / R15).
- **2026-08-07 adversarial review: 1 blocking + 3 majors + 4 minors applied** (ACLV 2.2(b) pin
  corrected — complement = {S₃}, the ada-with-nonempty-complement datum; betti(2) predicate;
  QQ-scoped identification with large-GF(p) honesty; duality oracle on the proper-subset
  example; B/D task merge). Verified live in the venv: rad²=0 A₅ has 9 indecomposables,
  `𝓛_A={S₁,S₂,P₂,P₃}`, `𝓡_A={S₄,S₅,P₄,P₅}`, complement `{S₃}` (pd 2, id 2), `e_λ={1,2,3}`,
  `e_ρ={3,4,5}`, gl.dim 4, and `D 𝓡_A = 𝓛_{A^op}` (with `|𝓡_A|=4` of 9).

---

## Addendum (2026-08-07, P61 contract): atlas `pd_le_1` / `id_le_1` vectors

P61 (the recognizer ladder, R18 — `docs/plans/2026-08-07-plan-61-recognizer-ladder.md`) needs,
per indecomposable, the two booleans this plan **already computes** for its parts sweep —
`pd X ≤ 1` and `id X ≤ 1` (the `_pd_le_1` / `_id_le_1` helpers, `X.projective_resolution(2)
.betti(2) == 0` and `injective_resolution(X,2).betti(2) == 0`). Expose them as two **additive**
`LeftRightAtlas` fields, index-aligned with `_modules` / the `left`/`right` records:

- `pd_le_1: tuple[bool, ...]` — `pd_le_1[i]` is `True` iff `pd(U[i]) ≤ 1`.
- `id_le_1: tuple[bool, ...]` — `id_le_1[i]` is `True` iff `id(U[i]) ≤ 1`.

**Zero new mathematics** — the sweep already builds `pd_ok` / `id_ok` (Task A Step 3); this is a
payload addition (assign the existing vectors into the dataclass, alongside `_leq`). It lets P61
compute the **shod QT2 route** (`A` is shod ⟺ every indec has `pd ≤ 1 or id ≤ 1`) as an
INDEPENDENT cross-check of the complement-empty route (survey Thm 4.1 (a)⟺(b)), and gives the
shod/ada witnesses honest per-module `{pd_le_1, id_le_1}` flags (rather than exact pd/id
numbers, which this plan does not compute past 2). Populated in **Task A** (the fields are
assembled with the other atlas fields — no task-boundary regression); keep them **engine-internal**
(P61 reads the dataclass in-process) so **no golden / canonical-key change** (`left_right_parts_block`
need not serialize them — if it later does, add them under the byte gate). Requested at P61's
adversarial-review round 2 (2026-08-07); this addendum is committed together with the P61 plan.
