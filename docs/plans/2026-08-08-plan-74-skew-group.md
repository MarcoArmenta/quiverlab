# Plan 74: Skew group algebras A⋊G — constructor + Ştefan conjugacy-class HH decomposition (P74 / R8)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking. Do the reference re-verification in
> the Record/Reference sections BEFORE writing any oracle pin — the flagged attributions
> below (`# PIN` / the withdrawn-preprint / the mis-attributed "CMRS" defect) must be
> resolved against the sources, not this plan's prose.

**Goal.** Make the **skew group algebra** `A⋊G` (the smash product `A#kG`) a first-class
quiverlab **input**, and compute its Hochschild (co)homology two ways that cross-check:

- **(1 — the constructor, the headline).** From a base algebra `A = kQ/I` and an **explicit**
  finite group action `G → Aut(A)` (never inferred), build `A⋊G` as a first-class
  structure-constant `Algebra` of dimension `|G|·dim A`, basis `{a·g : a ∈ A-basis, g ∈ G}`,
  product `(a·g)(b·h) = a·g(b)·(gh)`. **Characteristic-agnostic** (the smash product is
  defined in every characteristic). A genuine `kQ'/I'` presentation — the
  Reiten–Riedtmann/Demonet quiver — is recovered **best-effort** via the shipped, per-instance
  certified `Algebra.presented_form()` when `A⋊G` is basic and char-scope holds (no bespoke
  presented builder — live-verified below that `presented_form` recovers the orbit quiver
  automatically).

- **(2 — the Ştefan decomposition, the value-add).** Over `char k ∤ |G|` the Hochschild
  cohomology decomposes along **conjugacy classes**:
  `HH^n(A⋊G) ≅ (⊕_{[g]} HH^n(A, {}_gA))^{Z(g)}` (Ştefan 1995; Shepler–Witherspoon). Each
  summand `HH^n(A, {}_gA)` is served **today** by P52's `Bimodule.twisted(A, phi=g_matrix)` +
  `Algebra.hochschild_cohomology(coefficients=)` (P52's own docstring names P74 as this
  consumer). The **new engine work** is the `Z(g)`-action on the twisted HH classes + the
  **Reynolds invariants** functor (a bridge P52 does *not* ship — see the design finding).

- **The oracle spine (the two routes meet).** The **direct** computation of `HH^•(A⋊G)` by the
  shipped engine on the constructed algebra is the **trustworthy ground truth** (it needs no
  new math and works in any characteristic); the Ştefan **decomposition** is cross-checked
  **degreewise against it**. A mismatch is a loud bug in the decomposition, never a silent
  guess. A third, invariants-free oracle: when the action makes `A⋊G` an **orbit algebra**
  (e.g. a Nakayama algebra), `HH^•(A⋊G) ≡ HH^•` of that independently-built family.

**GUI.** ONE new **construction family** `SkewGroupAlgebra` (the `TrivialExtension`/`OppositeAlgebra`
no-code-INPUT precedent — a base + an explicit action; once built it flows through every
existing compute kind), PLUS one **algebra-only** compute kind `skew_group_hh` that renders the
conjugacy-class decomposition table with its honest char-scope. Served by **all three tiers**,
both runners byte-identical, **all four locales (en/es/fr/zh)** with exact key parity, canonical
keys stable, report rendering, and citations.

**Architecture.** TWO new source modules + one construction family + one compute kind; the HH
summands ride entirely on P52 (no new HH engine — only the group-action transport + invariants
on top of P52's chain-level accessors):

- **`src/quiverlab/families/skew_group.py`** (new — the constructor). Public:
  - `class QuiverAutomorphism` — a quiver automorphism `(vertex_perm, arrow_perm, arrow_scalars)`
    inducing an algebra automorphism; `.matrix(A)` → the `m×m` Domain automorphism matrix in
    `A`'s own basis (the form `Bimodule.twisted` consumes); `.check(A)` **certifies it is an
    algebra automorphism** (bijective, unital, `φ(uv)=φ(u)φ(v)` on structure constants;
    equivalently preserves `I`) — a per-instance certificate.
  - `class GroupAction` — a finite group given by **generators** (each a `QuiverAutomorphism`),
    closed under composition by a finite BFS (budget-capped, loud on non-closure); records
    `order` = `|G|`, the Cayley table, `conjugacy_classes`, `centralizer(g)`; `.check(A)`
    certifies every element is an algebra automorphism. The **matrix route**
    (`GroupAction.from_matrices`) is the power-user twin (explicit `m×m` automorphism matrices)
    — the quiver route is primary and the only one the no-code GUI exposes.
  - `skew_group_algebra(A, action, *, field=None) -> Algebra` — `A⋊G` in structure constants,
    `dim = |G|·dim A`, `basis_labels` `"<a>.<g>"`; `_validate` + the **dim-law certificate**.
  - `Algebra.skew_group(action)` delegate; `is_free_action(A, action) -> bool` (the finite
    idempotent-orbit check enabling the Galois-covering reading).
- **`src/quiverlab/hochschild/skew_group.py`** (new — the decomposition). Public:
  - `stefan_decomposition(A, action, top, *, side="coh", max_cells=4_000_000) -> SkewGroupHH` —
    per conjugacy-class rep `g`: `{}_gA = Bimodule.twisted(A, phi=g.matrix(A))`,
    `twisted_cohomology_classes(A, {}_gA, top)` (P52), the **`Z(g)`-action** on those classes
    (the new transport), **Reynolds invariants** (char `∤ |G|`), assemble `⊕`. `side="hom"`
    twin via `twisted_homology_classes`. Char-gated; cross-checked against the direct engine
    when `verify_direct=True`.
  - internal `_autom_action_on_classes(...)` — the chain map on the twisted bar cochain
    complex induced by an algebra automorphism `h`, read on cohomology classes.
- The `SkewGroupAlgebra` construction family is wired into `hpc/spec.py::_build_synthetic` +
  its Pyodide twin `docs/gui/runner.py` + `webapp/server/catalog.py::_SYNTHETIC_FAMILIES` (the
  `TrivialExtension`/`OppositeAlgebra`/`SkewGentleAlgebra` precedent — a base + flattened action
  params). The `skew_group_hh` kind is wired into `hpc/spec.py::_dispatch` + the twin +
  `webapp/server/schema.py` (grammar parse, a THIRD site) + `gui.js` ×2 + i18n ×4 +
  `trace/results_html.py` + one golden.

**Tech Stack.** Pure exact `Domain` linear algebra + finite group combinatorics. **No floats in
`src/`** (AST-gated by `tests/test_no_floats.py`): structure constants are exact Domain vectors,
automorphism matrices exact, dims/counts `int`, verdicts `bool`/`str`. All HH exact via P52
(`Bimodule.twisted`, `hochschild_cohomology(coefficients=)`, `twisted_cohomology_classes`); the
constructor via `Algebra.from_structure_constants`; the presentation via the shipped
`presented_form`.

---

## Records (verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R8 — Skew group / smash algebras A⋊G as inputs, with HH transfer.** [merged B-P8 + D-P2;
> keep] Object: constructor (A, G-action) → A⋊G (dim |G|·dim A) as a first-class Algebra; HH via
> the Ştefan conjugacy-class decomposition (⊕_{[g]} HH^n(A, {}_gA))^{Z(g)} — needs R4. HARD
> SCOPE: char k ∤ |G| (Reynolds averaging); the action is explicit input (never inferred);
> free-action detection on the quiver is a finite check enabling the Galois-covering reduction.
> Refs: Shepler–Witherspoon (Adv. Math., people.tamu.edu/~sjw/pub/ring.pdf); arXiv:0911.0938;
> Cibils–Marcos arXiv:math/0312214 (Proc. AMS 134 (2006) 39–50); CMRS arXiv:1804.02223. Oracles:
> trivial G byte-identity; dim law; Z/2 on dual numbers direct-vs-decomposition; Nakayama as
> orbit algebras. Size L.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P74 (Wave 4, tier δ, **needs P52**):

> **P74 — skew group algebras (R8, size L; needs P52).** Constructor (A, explicit G-action) →
> A⋊G as first-class Algebra; HH via Ştefan conjugacy-class decomposition (needs twisted
> coefficients); HARD scope char k ∤ |G|; free-action detection enabling the Galois-covering
> reduction. Refs: Shepler–Witherspoon; 0911.0938; Cibils–Marcos math/0312214; CMRS 1804.02223.
> Oracles: trivial-G byte-identity; dim law; Z/2 on dual numbers direct-vs-decomposition;
> Nakayama as orbit algebras. GUI: skew-group constructor preset.

`R4` = P52 (HH with bimodule coefficients), **merged on `dev`** (`77f6cc9`, Wave 1 complete) —
verified by live import at spec time (`quiverlab.hochschild.coefficients.Bimodule.twisted`).

---

## Reference re-verification (done at authoring; findings binding, `# PIN` = resolve at citation-add)

The metaplan standing rule (§1.5) requires the plan writer to re-verify citations. Findings — the
card's reference list has a **real defect** and a **material omission**, both resolved here:

1. **Ştefan, D. — "Hochschild cohomology of Hopf Galois extensions", J. Pure Appl. Algebra 103
   (1995), no. 3, 221–233. [MISSING FROM THE CARD — ADD.]** Verified this round (web): this is the
   ORIGINAL source of the **spectral sequence** for an H-Galois extension `A ⊇ A^{coH}`; the web
   search confirms that **for `H = kG` a group algebra the spectral sequence decomposes as a
   direct sum indexed by the conjugacy classes of `G`** — the origin of the R8 formula
   `(⊕_{[g]} HH^n(A, {}_gA))^{Z(g)}`. **Scope of the anchor (finding 3a — honest split):** I could
   verify (web abstract + secondary sources) that Ştefan 1995 constructs the spectral sequence and
   that the group-algebra case decomposes along conjugacy classes, but I could NOT fetch the paper
   body to transcribe the **additive `E₂` formula verbatim**. Therefore Ştefan 1995 is cited for
   the **spectral-sequence origin**, and the **explicit additive `(⊕_{[g]} …)^{Z(g)}` statement is
   anchored primarily on Shepler–Witherspoon (item 2)**, which states the skew-group-algebra HH
   decomposition directly. `# PIN`: read the Ştefan PDF at citation-add and, if it states the
   additive degeneration explicitly, promote it back to co-primary; the withdrawn `1804.02223`
   (item 4) may be cited AS withdrawn-context for the same statement, never load-bearing. (For
   `|G|` invertible the SS **degenerates at E₂** — the direct sum, no higher differentials — which
   is what this plan computes; the full multi-page SS is out of scope, honest note.)

2. **Shepler, A. V. & Witherspoon, S. — "Group actions on algebras and the graded Lie structure
   of Hochschild cohomology", J. Algebra 351 (2012), 350–381 (arXiv:0911.0938).** **Verified this
   round** (arXiv abstract): title and authors exact; studies `HH^•` of the skew group algebra of
   a finite group acting by automorphisms, the Gerstenhaber bracket, and the abelian-group case.
   This is THE **load-bearing anchor for the additive conjugacy-class decomposition + the
   `Z(g)`-invariance** (it does the skew-group-algebra HH decomposition explicitly; Ştefan item 1
   supplies the spectral-sequence origin). **Venue fix (this fix round, finding 3b):** the earlier
   stub said "Adv. Math."; the paper is **J. Algebra 351 (2012), 350–381** — corrected. `# PIN`:
   the **left-vs-right twist convention** for `{}_gA` and the precise `Z(g)`-action formula
   transcribed verbatim before `_autom_action_on_classes` is coded. **Correction (finding 2a — do
   NOT claim the direct oracle arbitrates the convention at summand level):** the summand
   dimensions `dim HH^•(A, {}_gA)` are provably **insensitive** to the left-vs-right twist (the
   bimodule iso `{}_φA_ψ ≅ {}_1A_{φ⁻¹ψ}` makes them twist-symmetric — live-verified below on
   order-2 AND order-3 examples). The convention is therefore **source-pinned** (from this paper),
   NOT arbitrated by any summand; the only discriminator is the **assembled decomposition + the
   per-summand split** (which the correctly-signed `Z(g)`-transport determines) checked against
   the direct engine — see the twist-convention subsection and the live-facts table.

3. **Cibils, C. & Marcos, E. N. — "Skew category, Galois covering and smash product of a
   `k`-category", Proc. Amer. Math. Soc. 134 (2006), no. 1, 39–50 (arXiv:math/0312214).**
   **Verified this round** (arXiv abstract): title, authors, PUBLISHED venue confirmed. This is
   the anchor for the **smash-product / skew-category construction, the free action, and the
   Galois covering** — NOT for the HH conjugacy decomposition (its abstract does not state that;
   it is the categorical construction + duality framework). Cite it for the constructor +
   free-action + covering, not for the decomposition. `# PIN`: issue/pages at citation-add.

4. **CARD DEFECT — "CMRS 1804.02223" is WRONG on two counts. [RESOLVED HERE.]** Verified this
   round (arXiv abstract): **arXiv:1804.02223 is Cibils–Marcos** (TWO authors — NOT "CMRS" =
   Cibils–Marcos–Redondo–Solotar, four authors), "Hochschild–Mitchell (co)homology of skew
   categories and of Galois coverings", **and it has been WITHDRAWN by the author** (a revised
   version with a "resolving category" viewpoint was promised). It *does* give the conjugacy-class
   decomposition (stated to match M. Lorenz's), but **a withdrawn preprint must not be a
   load-bearing anchor**. **Resolution:** drop 1804.02223 as an anchor; the decomposition is
   anchored on Ştefan (item 1) + Shepler–Witherspoon (item 2). Optionally cite 1804.02223 in an
   annotation ONLY as "withdrawn; superseded" — never as a pin. The worker must not transcribe
   any number from it.

5. **Marcos, E. N., Martínez-Villa, R. & Martins, M. I. R. — "Hochschild cohomology of skew
   group rings and invariants", arXiv:math/0312460.** **Verified this round** (arXiv abstract):
   establishes the ring monomorphism `HH^•(A)^G ↪ HH^•(A⋊G)` (and the Galois-covering
   monomorphism). A useful **invariants-side** anchor (the identity summand `HH^•(A)^G` is a
   genuine subalgebra of `HH^•(A⋊G)` — a self-cert the plan can check). Add as a supporting ref.

6. **Shepler–Witherspoon survey (`people.tamu.edu/~sjw/pub/ring.pdf`).** The card's survey URL:
   `# PIN` the exact title/venue at citation-add (I could not uniquely resolve `ring.pdf` to a
   published survey this round). **Already shipped and reusable:** `witherspoon_gsm204` (S.
   Witherspoon, *Hochschild Cohomology for Algebras*, GSM 204, AMS 2019) IS the canonical survey
   textbook — verified live in `registry.py` — use it as the survey anchor; only add a distinct
   `ring.pdf` entry if its bibliographic identity is confirmed.

7. **Positive-characteristic honesty.** Shepler–Witherspoon, "Gerstenhaber brackets for skew
   group algebras in positive characteristic" (arXiv:1905.09613) confirms the MODULAR case
   (char `| |G|`) is genuinely different — the clean conjugacy-class decomposition **fails** and
   the bracket needs a twisted product. Cited (annotation) to justify the loud modular refusal of
   the *decomposition* route (the constructor + direct HH stay available). `# PIN`: venue.

**Already shipped and reused (verified live in `registry.py`):** `witherspoon_gsm204`
(WitherspoonGSM204), the `cibils_*` cluster, `assem_book` (ASS2006), `crs_trivial_ext_hh1`.
**New keys to add:** `stefan_hopf_galois` (JPAA 103, 1995), `shepler_witherspoon_group_actions`
(arXiv:0911.0938), `cibils_marcos_smash` (Proc. AMS 134, arXiv:math/0312214),
`marcos_mv_invariants` (arXiv:math/0312460). (1804.02223 is NOT added as an anchor.)

---

## Mathematical foundation (the definitions, exactly — the plan's ground truth)

Throughout `A = kQ/I` is basic, connected, finite-dimensional over an exact `Domain`, with the
quiverlab conventions (right modules, composition left-to-right `a*b` = first `a` then `b`).
`G` is a **finite** group acting on `A` by **algebra automorphisms**, given **explicitly** (never
inferred). `m = dim_k A`.

### The skew group algebra (smash product) `A⋊G`

- **As a vector space** `A⋊G = A ⊗_k kG`, with basis `{a·g : a ∈ (A-basis), g ∈ G}`,
  `dim_k(A⋊G) = |G|·dim_k A` (**the dim law**).
- **Multiplication** `(a·g)(b·h) = (a · g(b)) · (gh)`, where `g(b)` is the automorphism `g`
  applied to `b ∈ A`. `1_A·1_G` is the unit. Associativity is a consequence of `g` being an
  automorphism and the group law (certified per instance by `_validate`).
- **Characteristic-agnostic.** The structure constants are exact in ANY characteristic; the
  constructor and the DIRECT `HH^•(A⋊G)` never need `|G|` invertible. Only the **decomposition**
  and the **covering reduction** need `char k ∤ |G|` (Maschke/Reynolds — item below).
- **Presentation (best-effort, certified).** When `A⋊G` is basic and char-scope holds, its
  Reiten–Riedtmann/Demonet quiver `kQ'/I'` is recovered by the shipped `presented_form()`
  (`dim(kQ'/I') = dim(basic(A⋊G))` + a multiplicativity certificate). **Live-verified**: for the
  `Z/2`-on-dual-numbers example, `presented_form` on the structure-constant `A⋊G` returns the
  2-vertex quiver `1 ⇄ 2` (arrows `a1:1→2`, `a2:2→1`) with `I' = ⟨a1a2, a2a1⟩` — the orbit
  Nakayama algebra, with NO bespoke code. When `A⋊G` is **non-basic** (e.g. `kG` contributes
  matrix blocks — non-abelian `G` with higher-dimensional irreducibles), `presented_form`
  returns the presentation of the **basic algebra** of `A⋊G` (Morita-equivalent; `HH` is
  Morita-invariant, so the HH numbers are unaffected, but the **dim law is a statement about the
  full `A⋊G`**, not its basic algebra). The constructor therefore ALWAYS returns the full
  structure-constant `A⋊G`; the presentation is an optional certified view.

### The G-action, exactly (how a user specifies it — design decision (1))

A `G`-action is data on the **quiver**, lifted to `A`, then verified:

- **Quiver automorphism (primary, no-code).** A triple `(π, ρ, s)`: `π` a permutation of `Q_0`
  (vertices), `ρ` a permutation of `Q_1` (arrows) compatible with `π`
  (`source(ρα)=π(source α)`, `target(ρα)=π(target α)`), and `s: Q_1 → k^×` **arrow scalars**
  (nonzero exact entries). It induces the algebra automorphism `φ` on paths by
  `φ(e_v) = e_{π(v)}`, `φ(α) = s(α)·ρ(α)`, extended multiplicatively. **The arrow scalars are
  essential even for the flagship example** — `σ: x ↦ −x` on `k[x]/(x²)` is the identity
  permutation with the scalar `s(x) = −1`; a pure permutation cannot express it. `φ` is an
  algebra automorphism **iff** it preserves `I` (checked on the relation generators).
- **Matrix route (power-user).** An explicit `m×m` Domain automorphism matrix (columns = images
  of the `A`-basis) — the form `Bimodule.twisted` already consumes. `GroupAction.from_matrices`.
  Not exposed in the no-code GUI (matrices of the whole basis are not homework-grade input).
- **The group.** `G` is generated by a finite list of quiver automorphisms; the constructor
  **closes under composition** by BFS (budget-capped; loud `QuiverlabError` if it does not close
  within the budget or is infinite), builds the Cayley table, and derives the **conjugacy
  classes** and **centralizers** `Z(g)` combinatorially. Cyclic `G = ⟨σ⟩` needs only `σ`.
- **Certificate.** `GroupAction.check(A)` verifies (per instance, loud) that (i) each generator
  is a bijective unital algebra automorphism preserving `I`, and (ii) the closed set is a group
  of the claimed order. **Verification that the action is by algebra automorphisms is a
  per-instance certificate**, not assumed.

### The Ştefan conjugacy-class decomposition (design decisions (2)–(4))

For `char k ∤ |G|` (Ştefan 1995; Shepler–Witherspoon 0911.0938):
```
HH^n(A⋊G) ≅ ( ⊕_{g ∈ G} HH^n(A, {}_gA) )^G  =  ⊕_{[g] conj. class} HH^n(A, {}_gA)^{Z(g)} .
```
- **The summand coefficient `{}_gA`** is `A` with the bimodule structure twisted by `g` on one
  side — `a · x · b = g(a) x b` (left twist) — realized by **P52's `Bimodule.twisted(A,
  phi=g.matrix(A))`** (the `phi=g_matrix` path P52's docstring names for P74).
  - **Twist convention — NOT summand-arbitrable (finding 2a, live-verified; corrects the earlier
    "direct oracle arbitrates loudly" claim).** The bimodule iso `{}_φA_ψ ≅ {}_1A_{φ⁻¹ψ}` (apply
    `φ⁻¹` to the underlying space) means the **summand dimensions** `dim HH^•(A, {}_gA)` are
    **independent of the left-vs-right twist choice**: e.g. `HH^•(A, phi=g) = HH^•(A, psi=g⁻¹)`.
    For an **involution** (`g²=1`, `g⁻¹=g`) this makes left-`g` and right-`g` **provably
    indistinguishable** at the summand level — every planned `σ` was order 2, so nothing there
    could witness a wrong convention. I therefore tested the **order-3** case `Z/3` on `k[x]/(x³)`
    over `GF(7)` (`σ(x)=2x`, `2³≡1`): `HH^•(A, phi=σ)`, `HH^•(A, psi=σ)`, AND `HH^•(A, phi=σ²)`
    are **all `[1,1,1,1,1]`** — still indistinguishable, because `HH^•(A, A_σ) = HH^•(A, A_{σ²})`
    for this family. **Conclusion:** the convention is **source-pinned** from Shepler–Witherspoon,
    NOT arbitrated by any summand dimension; the only discriminator is the **assembled
    decomposition + the per-summand split** (below), which depends on the correctly-signed
    `Z(g)`-transport and is checked against the direct engine. The `# PIN` stands, but its
    validation is the assembled total, not a summand comparison.
- **The `G`-action on the sum** permutes summands by conjugation `g ↦ hgh^{-1}` (for **abelian**
  `G`, every class is a singleton and `Z(g)=G`, so the sum is `⊕_{g} HH^n(A, {}_gA)^G`) and acts
  **within** each summand by the automorphism transport `h·(m⊗a_1⊗…⊗a_n) = h(m)⊗h(a_1)⊗…⊗h(a_n)`
  **together with the induced sign/scalar on the resolution generators** (the crux — a
  live-verified degree-1 witness below shows the resolution sign decides which summand carries a
  class). Taking invariants (`(-)^G`) needs the **Reynolds idempotent** `e = (1/|G|)Σ_h h` — hence
  **char `∤ |G|`**.
- **Design decision (3): the DIRECT computation is the trustworthy primary; the DECOMPOSITION
  is the value-add, cross-checked against it.** `HH^•(A⋊G)` computed directly by the shipped
  engine on the constructed algebra is the ground truth (no new math, any characteristic). The
  decomposition is the "interesting" content and is **verified degreewise against the direct
  result** — a mismatch is a loud bug in the decomposition, never a silent number. This inverts
  the naive "spectral-sequence-first" reading in favour of the honest engineering ground truth.
- **Design decision (4): char `| |G|` (modular) refuses the DECOMPOSITION + covering loudly**,
  with the modular-invariant-theory honesty: Maschke fails, the Reynolds average is unavailable,
  and `HH^•(A⋊G) ≇ (⊕HH^•(A,{}_gA))^G` in general (Shepler–Witherspoon 1905.09613 — the modular
  bracket needs a twisted product). **The constructor and the direct `HH^•(A⋊G)` remain
  available in every characteristic** (a card refinement: the char scope is on the decomposition
  route, not the constructor — evidenced by the char-free structure-constant build).
- **P70/P71 seam (this plan ships only DIMENSIONS, not the graded Lie structure).** Shepler–
  Witherspoon 0911.0938 is fundamentally about the **graded Lie / Gerstenhaber-bracket** structure
  of `HH^•(A⋊G)`; this plan computes only the **decomposition of dimensions** (and the direct HH),
  NOT the bracket. The `HH¹`-Lie-algebra classification of `A⋊G` and the bracket transport are
  **downstream of P70/P71** (`hh1_lie` + `HH•`-as-Lie-module) and P51 (bracket beyond the window),
  which run on `A⋊G` as an ordinary input once this constructor exists. The plan states this seam
  explicitly (verification honest-scope entry) so no reviewer mistakes the dimension decomposition
  for the Lie/bracket content — that is a later plan's job, enabled by this one's constructor.

### The design finding — what P52 ships vs. what this plan builds

**P52 ships (verified live, all reusable today):**
- `Bimodule.twisted(A, phi, psi)` → the summand coefficient `{}_gA` (`.check()` certifies the
  `A^e`-axioms). P52's docstring: *"P74 (skew group algebras) — `twisted` with `phi=g_matrix`."*
- `Algebra.hochschild_cohomology(top, coefficients=M)` / `hochschild_homology(...)` → the summand
  dims `HH^•(A, {}_gA)`.
- `twisted_cohomology_classes(A, M, top)` / `twisted_homology_classes(...)` → the **chain-level**
  data per degree: `basis`, `coboundary`/`boundary`, and a basis of **class representatives** in
  the twisted bar (co)chain complex — exactly the handles a `G`-action needs.

**P52 does NOT ship (the bridge this plan builds):**
1. the **`Z(g)`-action on `HH^•(A, {}_gA)`** — the chain map on the twisted bar (co)chain complex
   induced by an algebra automorphism `h ∈ Z(g)`, and its induced map on classes; and
2. the **invariants functor** `(-)^{Z(g)}` (Reynolds average, char `∤ |G|`) + the assembly `⊕`.

So the conjugacy-class **summands are computable today**; the **`Z(g)`-action + invariants +
assembly** is the genuinely new engine slice, feasible on top of the P52 chain-level accessors.
This is an honest "what exists vs what must be built" statement, not a blocker.

### The bridge is genuinely non-trivial — a live-verified per-class witness (finding 1a/2b)

A critic's decisive worry: if `g` acts as `+1` on every `HH^n(A,A)`, the **identity summand alone
carries the total** and the twisted invariants are `0` everywhere — nothing would distinguish a
correct bridge from a **dead** one (twisted `≡ 0`). The flagship COHOMOLOGY exhibits exactly this
trap (`σ` acts as `+1` on `HH^1(A,A)=Der`, so the twisted summand is `0` there). The fix is a
degree where a **twisted summand is provably nonzero**:

- **Flagship, HOMOLOGY, degree `n=1` (fully hand-derived, no bridge needed, confirmed two ways).**
  `HH_1(A,A) ≅ Ω¹_{A/k}` for commutative `A` (standard). For `A=k[x]/(x²)`,
  `Ω¹ = A\,dx / (d(x²)) = A\,dx/(2x\,dx) = k\,dx` (dim 1, `x\,dx = 0`). `σ(dx) = d(σx) = d(−x) =
  −dx`, so **`σ` acts as `−1`** on `HH_1(A,A)`. Hence the identity coinvariant
  `HH_1(A,A)_G = Ω¹/(σ−1)Ω¹ = k\,dx / (−2\,k\,dx) = 0` (`−2` invertible). But **`HH_1(A⋊G)`
  direct `= 1`** (live: the flagship `HH_• = [2,1,1,1]`). Therefore the **twisted summand
  `HH_1(A,{}_σA)^G = 1 − 0 = 1` — NONZERO**, and it **alone** carries `HH_1(A⋊G)`. A dead bridge
  gives `0 ≠ 1`; a **sign-wrong** bridge gives the split `(1,0)` instead of the correct `(0,1)`
  (same total `1`, wrong per-summand). Independently: the minimal 2-periodic resolution forces the
  generator sign `σ_P(ξ_1) = −ξ_1`, giving `σ = −1` on `HH_1(A,A)` and `σ = +1` on the twisted
  class — same conclusion.
- **Homology degree 0 (finding 2e, hand-derived).** Flagship: `HH_0(A,A)_G = (A)_G = 1` (`σ`=
  `diag(1,−1)`, `(σ−1)A = span{x}`), twisted `HH_0(A,{}_σA)_G = 1`, total `2 =` direct `HH_0(A⋊G)`
  — again **twisted nonzero**. `Z/3` on `k[x]/(x³)` (`GF(7)`): `HH_0(A⋊G)` direct `= 3` splits as
  `1 + 1 + 1` (identity `+` `σ`-twisted `+` `σ²`-twisted, each `1`).
- **The per-summand oracle (beyond total == direct).** `stefan_decomposition` must report
  `summands` per class, and the test pins the **split**, not only the sum: at flagship homology
  `n=1`, `summands = [{g=e: 0}, {g=σ: 1}]` (NOT `[1,0]`); at `n=0`, `[1,1]`. This is the check a
  dead-or-sign-wrong bridge fails. Additionally the **identity-summand monomorphism**
  `HH^n(A)^G ↪ HH^n(A⋊G)` (Marcos–Martínez-Villa, cohomology) is a degreewise dim `≤` self-cert,
  and a **QPA input-level crosscheck** feeds QPA the presented `A⋊G` (`= kQ'/I'`) and compares its
  `HH` dims to the assembled total (QPA has no isotypic/decomposition surface — honest scope).

### The bridge is exercised for non-abelian `G` too — the smallest pin (finding 2c)

Non-abelian `G` is **in scope** (not downgraded), pinned by **`S₃` acting on `k³ = k×k×k`**
permuting the three factors (live-verified over `GF(7)`): `dim(k³⋊S₃) = 6·3 = 18` (dim law) and
`HH^0(k³⋊S₃) = 2` (`k³⋊S₃` is semisimple over `char ∤ 6`, Morita-equivalent to `k[S₂] ≅ k×k` via
the point-stabilizer `S₂`, so `2` simple factors ⇒ `HH^{≥1}=0`, `HH^0 = 2`). The **conjugacy-class
shape drives the decomposition**: `S₃` has classes `{e}` (`Z(e)=S₃`, order 6), `{(12),(13),(23)}`
(`Z=⟨\text{transp}⟩≅Z/2`), `{(123),(132)}` (`Z=⟨\text{3-cycle}⟩≅Z/3`) — so `Z(g)⊊G` genuinely
occurs and the per-class `(-)^{Z(g)}` (a **proper centralizer**, not all of `G`) is exercised. The
direct `HH^{≥1}` of the dim-18 algebra exceeds the bar-engine cap (`88434×5202 > max_cells` at
`d^2`) — an honest `status="budget"` example; the dim law + `HH^0` + the class shape are the pins.

### Free action + Galois covering (honest slice)

`is_free_action(A, action)` is the finite check: no non-identity `g ∈ G` fixes a primitive
idempotent (vertex). When the action is free (on `Q_0`), `A⋊G` is the total algebra of a
**Galois `G`-covering** of the orbit algebra `A/G` (Cibils–Marcos math/0312214), and (over
`char ∤ |G|`) `A⋊G` is Morita-equivalent to the smash — the **orbit-algebra** picture. **v1
scope:** ship the detection + the **genuine free-orbit oracle** (below); the general
covering-reduction HH transport is honest-scope (a ledger item if it proves heavy), because the
Ştefan decomposition already gives the HH answer within scope.

- **The flagship is NOT a free action (finding 2d).** `Z/2` on `k[x]/(x²)` fixes the single vertex
  (`σ` acts trivially on `Q_0`); `is_free_action` returns `False`. Its `A⋊G = ` 2-vertex Nakayama
  arises from the **`kG`-idempotent splitting** `k[Z/2] ≅ k×k` (char `≠ 2`), NOT from a free vertex
  orbit — the plan labels it as such, correcting the card's "orbit algebra" framing for this case.
- **A genuinely FREE orbit oracle (live-verified, QQ).** `Z/2` acting on the 2-cycle Nakayama
  `A = k(1⇄2)/\mathrm{rad}^2` (dim 4) by **swapping the two vertices** (`1↔2`, `a↔b`) is **free on
  `Q_0`** (vertex orbit `{1,2}` merges). `A⋊(Z/2)` has dim `8`, and **`presented_form` recovers the
  orbit quiver `k[x]/(x²)`** — 1 vertex, 1 loop `a1`, relation `a1²=0` — with
  `HH^•(A⋊G) = [2,1,1] = HH^•(k[x]/(x²))` and `is_selfinjective True`. This is the honest
  "**Nakayama/local as orbit algebra**" oracle: free action ⇒ `presented_form = ` orbit quiver,
  HH matches the orbit algebra. **Elegant duality (worth a note in the report):** this free example
  and the flagship are **inverse constructions** — the flagship's `A⋊G` (2-cycle Nakayama) is this
  example's base `A`, and this example's `A⋊G` (`~ k[x]/(x²)`) is the flagship's base.
  (Over `GF(7)` this same free build hits the `char 7 ≤ dim 8` `presented_form` refusal — the
  char-`>`-dim caveat, finding 3c — so the orbit-quiver leg runs over QQ.)

---

## Live-verified facts (all recomputed in the venv at spec time; QQ unless noted, plus GF(7))

Flagship example: `A = k[x]/(x²)` (dual numbers, `truncated_polynomial(2)`, basis `[e_1, x]`) over
`QQ`; `G = Z/2 = ⟨σ⟩`, `σ(x) = −x` (identity vertex-perm, arrow scalar `s(x) = −1`); char `0 ∤ 2`.
The `σ`-matrix in `A`'s basis is `[[1,0],[0,−1]]`.

| object | value (degrees `n = 0..3`) | oracle class |
|---|---|---|
| `dim(A⋊G) = |G|·dim A` | `2·2 = 4` ✓ | selfcert (dim law) |
| `HH^•(A, A)` (ordinary) | `[2, 1, 1, 1]` | literature (dual numbers) |
| `HH^•(A, A)` via `Bimodule.regular(A)` | `[2, 1, 1, 1]` (`≡` ordinary — **P52 `M=A` oracle**) | crossengine |
| `HH^•(A, {}_σA)` (`Bimodule.twisted(A, phi=σ)`) | `[1, 1, 1, 1]` | crossengine (summand) |
| `HH_•(A, {}_σA)` | `[1, 1, 1, 1]` | crossengine (summand) |
| raw `⊕` pre-invariants `HH^•(A,A) ⊕ HH^•(A,{}_σA)` | `[3, 2, 2, 2]` | (intermediate) |
| **`HH^•(A⋊G)` DIRECT** (shipped engine) | **`[1, 1, 1, 1]`** | crossengine (ground truth) |
| **`HH_•(A⋊G)` DIRECT** | **`[2, 1, 1, 1]`** | crossengine (ground truth) |
| **decomposition `(⊕ HH^•)^G`** | **`[1, 1, 1, 1]`** — matches DIRECT ✓ | crossengine (the R8 oracle) |
| orbit Nakayama `1⇄2 / rad²` `HH^•` | `[1, 1, 1, 1]`; `HH_•` `[2, 1, 1, 1]`; self-injective | crossengine (orbit) |
| `presented_form(A⋊G)` | `1⇄2`, `I'=⟨a1a2, a2a1⟩` (auto-recovered) | selfcert (dim + multiplicativity) |
| **flagship HOM `n=1` per-summand split** (bridge witness, finding 1a/2b) | identity `HH_1(A,A)^G = 0`; **twisted `HH_1(A,{}_σA)^G = 1`**; total `= 1 =` DIRECT | crossengine + hand |
| **`Z/3` on `k[x]/(x³)` / `GF(7)`** twist arbiter (finding 2a): `HH^•(A,phi=σ)` vs `HH^•(A,psi=σ)` vs `HH^•(A,phi=σ²)` | **all `[1,1,1,1,1]`** — convention NOT summand-arbitrable | selfcert (duality) |
| `Z/3` on `k[x]/(x³)`: `HH^•(A⋊G)` / `HH_•(A⋊G)` DIRECT `[0..2]` | `[1,1,1]` / `[3,2,2]`; `HH_0` splits `1+1+1` | crossengine |
| **`S₃` on `k³`** (non-abelian, finding 2c): dim law / `HH^0` | `dim = 18` ✓ / `HH^0 = 2` (Morita `k×k`) | crossengine + literature |
| **FREE `Z/2`-swap on `k(1⇄2)/rad²`** (finding 2d): dim / `presented_form` / `HH^•` | `dim 8` / orbit `= k[x]/(x²)` / `[2,1,1]`, self-inj | crossengine (free orbit) |

**The decomposition oracle, degree 0 hand-verified (the piece the bridge must reproduce):**
- COHOMOLOGY: `HH^0(A,A) = Z(A) = A = {1, x}`; `σ` acts `diag(1, −1)`; `(A)^G = span{1}`, **dim 1**.
  `HH^0(A,{}_σA) = {a : σ(b)a = ab ∀ b} = span{x}` (**dim 1**); `σ` acts `−1`, so `(span{x})^G = 0`.
  `HH^0(A⋊G) = 1 + 0 = 1` = DIRECT. ✓ (Here the identity summand carries it; twisted is 0.)
- HOMOLOGY (finding 2e): `HH_0(A,A)_G = (A)_G = 1` (`(σ−1)A = span{x}`); `HH_0(A,{}_σA)_G = 1`
  (`M/[A,M] = span{1̄}`, `σ`-coinv `1`); `HH_0(A⋊G) = 1 + 1 = 2` = DIRECT. ✓ (**Here BOTH summands
  contribute — the twisted one is nonzero.**)
- **HOMOLOGY `n=1` (the decisive n≥1 witness, finding 1a/2b):** `HH_1(A,A) ≅ Ω¹ = k\,dx` with
  `σ(dx) = −dx`, so identity coinvariant `= Ω¹/(−2)Ω¹ = 0`; DIRECT `HH_1(A⋊G) = 1`; therefore
  **twisted `HH_1(A,{}_σA)^G = 1`, nonzero, carrying the whole degree.** A dead bridge gives
  `0 ≠ 1`; a sign-wrong bridge gives the split `(1,0)` not `(0,1)` (see the design-finding section).

**Two invariants-free oracles (need no bridge — robust ground truth):**
- **DIRECT ≡ orbit Nakayama.** `A⋊G` (`Z/2` on dual numbers) equals the 2-vertex `rad²`-zero
  self-injective Nakayama algebra as an algebra (`presented_form` recovers `1⇄2 / rad²`), and its
  HH `[1,1,1,1]/[2,1,1,1]` matches the DIRECT `A⋊G` HH — **"Nakayama as orbit algebras"** (the
  R8 card oracle), verified two ways (structure-constant build + independent family).
- **Trivial-`G` byte-identity.** For `G = {1}`, `A⋊{1}` has basis `{a·1}` and product
  `(a·1)(b·1) = ab·1`, so its structure-constant table is **byte-identical** to `A.T` and its HH
  is byte-identical to `A`'s — a construction invariant (test builds `A⋊{trivial}` and asserts
  `T`-equality + HH-equality).

**What I could NOT fully verify at spec time (labelled to-be-verified-at-implementation):**
- The **general `Z(g)`-action machinery** `_autom_action_on_classes` — I hand-derived the flagship
  HOM `n=0`/`n=1` per-summand splits (two ways: Kähler differentials + the minimal-resolution
  generator sign) and the `Z/3` HOM `n=0` split (`1+1+1`), which pin the bridge's non-triviality;
  the general implementation across all degrees/classes IS the plan's Task 3, TDD-verified against
  the DIRECT oracle + the per-summand pins.
- **The left-vs-right twist convention** for `{}_gA` — `# PIN` to Shepler–Witherspoon. **NOT
  summand-arbitrable** (finding 2a, live-verified: order-2 AND order-3 `k[x]/(x^n)` give
  twist-symmetric summand dims via `{}_φA_ψ ≅ {}_1A_{φ⁻¹ψ}`); the convention is source-pinned and
  validated only through the assembled total + per-summand split against DIRECT.
- **Deeper degrees / larger algebras** — the `A⋊G` bar complex grows as `|G|·dim A`; the dim-9
  (`Z/3`) direct capped at degree 2 and the dim-18 (`S₃`) at degree 0 (`max_cells`), honest
  `status="budget"`. Deep degrees are TDD-verified against the DIRECT oracle at implementation
  (or via the cheaper base-`A` summand route, which caps later).

---

## Scope gates (contractual; every boundary is a loud typed refusal)

| surface | scope | refusal |
|---|---|---|
| `QuiverAutomorphism.check` / `GroupAction.check` | the triple is a genuine algebra automorphism of `A` (preserves `I`), and the generators close to a finite group within budget | not-an-automorphism / non-closing / infinite → `QuiverlabError` (loud, with the failing relation or the budget) |
| `skew_group_algebra(A, action)` | quiver-presented `A` (`A.quiver is not None`), explicit `action`; **any characteristic** | presentation-less base → loud (the action is quiver data); dim-law certificate failure → loud (internal invariant) |
| `Algebra.presented_form()` on `A⋊G` (the kQ'/I' view) | `A⋊G` basic **and** char 0 or char > dim | non-basic → returns the **basic algebra**'s presentation (documented; HH unaffected); off char-scope → the shipped `presented_form` loud refusal |
| `stefan_decomposition(A, action, top)` (the DECOMPOSITION + the `Z(g)`-invariants) | **char k ∤ |G|** (Reynolds) | **char `| |G|` → loud `QuiverlabError`** naming the modular-invariant-theory obstruction (the clean decomposition fails; use the DIRECT route) |
| the DIRECT `HH^•(A⋊G)` (via the shipped engine on the built algebra) | any characteristic, any admissible `A⋊G` | inherits the shipped engine's honest caps (`max_cells`, CS depth) — never a guessed number |
| `is_free_action` / covering reduction | free action on `Q_0`; covering HH transport is **honest-scope v1** (detection + orbit oracle only) | non-free → `is_free_action` returns `False` (not an error); the general covering HH reduction is ledger-scoped |
| over `GF(p)` | the twisted-coefficient / `presented_form` char caveats (`char > dim` where an iso/decompose is invoked) | `char ≤ dim` → the shipped primitives raise loudly (never a silent wrong verdict) |

All batteries run over **QQ**; a `GF(32003)` parity spot-check is kept for the dim-law +
direct-HH legs (distinct enough to avoid the `char ≤ dim` raises).

---

## API surface (public via `import quiverlab`; exact only)

```python
# --- the constructor (src/quiverlab/families/skew_group.py) ---
class QuiverAutomorphism:
    def __init__(self, vertex_perm: dict, arrow_perm: dict, arrow_scalars: dict | None = None)
    def matrix(self, A) -> list[list]           # m x m Domain automorphism matrix (A's basis)
    def check(self, A) -> bool                  # certifies: algebra automorphism preserving I

class GroupAction:
    generators: list                            # [QuiverAutomorphism]
    order: int                                  # |G|
    elements: list                              # closed set (as automorphisms), Cayley-indexed
    conjugacy_classes: list                     # [[element indices], ...]
    def centralizer(self, g_index) -> list      # Z(g) element indices
    def check(self, A) -> bool                  # every element an algebra automorphism; closed
    @classmethod
    def cyclic(cls, generator, order) -> "GroupAction"
    @classmethod
    def from_matrices(cls, matrices) -> "GroupAction"      # power-user (explicit m x m)

def skew_group_algebra(A, action, *, field=None) -> Algebra      # A⋊G, dim |G|*dim A
def is_free_action(A, action) -> bool

# Algebra delegate (core/algebra.py, thin lazy-import):
Algebra.skew_group(action) -> Algebra

# --- the Ştefan decomposition (src/quiverlab/hochschild/skew_group.py) ---
def stefan_decomposition(A, action, top, *, side="coh",
                         max_cells=4_000_000, verify_direct=False) -> SkewGroupHH
@dataclass
class SkewGroupHH:
    algebra: object                 # the base A
    action: object                  # the GroupAction
    side: str                       # "coh" | "hom"
    top: int
    dims: list[int]                 # assembled (⊕ HH^n(A,{}_gA))^{Z(g)}, degrees 0..top
    summands: list                  # per class rep: {"g": label, "hh": [...], "inv": [...]}
    direct_dims: list[int] | None   # HH^•(A⋊G) by the shipped engine (when verify_direct)
    agrees: bool | None             # dims == direct_dims (the R8 oracle) when verify_direct
    status: str                     # "complete" | "budget" | "modular" | "unsupported"
    note: str
```

`SkewGroupHH.__bool__` is NOT defined (a data report, not a predicate) — mirroring `HHTable`.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P52 (HH with bimodule coefficients) is merged on `dev`** (`77f6cc9`) — verified by live
  import at spec time. This plan consumes, with the exact signatures verified in `dev`:
  - `quiverlab.hochschild.coefficients.Bimodule`: `twisted(A, phi=, psi=, name=)`, `regular(A)`,
    `dual(A)`, `check()`, `describe()`, `invariants_dim()`/`coinvariants_dim()`;
    `twisted_cohomology_classes(A, M, top)` / `twisted_homology_classes(A, M, top)` (per-degree
    `basis` / `coboundary`|`boundary` / `classes`).
  - `Algebra.hochschild_cohomology(top, coefficients=M)` / `hochschild_homology(...)` (the
    `coefficients=` hook + provenance stamping).
  - `Algebra.from_structure_constants(T, unit, field=, basis_labels=)`; `Algebra._validate`;
    `Algebra.presented_form()` / `gabriel_quiver()` (`core/basic.py`); `Quiver(...).algebra(...)`.
  - `families.basic.truncated_polynomial`; `Algebra.is_selfinjective` (oracle corroboration).
  Branch `plan-74-skew-group` off `dev`.
- **Honest semi-decision contract (metaplan §1.3).** Every report is **complete** iff its scope
  gate holds and its budget did not trip; otherwise `status ∈ {"budget","modular","unsupported"}`
  with a note — **never a guessed number/decomposition**.
- **Char scope is load-bearing.** The **constructor + DIRECT HH** are char-agnostic. The
  **decomposition + covering** are char `∤ |G|` only (loud `"modular"` refusal otherwise). The
  `presented_form` view inherits the char `> dim` caveat. Batteries over **QQ**; `GF(32003)`
  parity on the char-free legs only.
- **No floats in `src/`.** Structure constants / automorphism matrices are exact Domain entries;
  arrow scalars exact (`int`/`"-1"`/`"1/2"` strings, never floats — `reject_inexact`); dims/orders
  `int`; verdicts `bool`/`str`. Client-side float only (`gui.js`, exempt).
- **Composition is left-to-right** (`a*b` = first `a` then `b`); the smash product
  `(a·g)(b·h) = a·g(b)·gh` respects it. All refusals are `QuiverlabError`.
- **Plan-32 markers** (orthogonal; `oracle_*` FORBIDDEN in `tests/{webapp,gui,hpc}` per
  `test_oracle_classes.py`):
  - `oracle_literature`: `HH^•(k[x]/(x²)) = [2,1,1,1]`; the dim law `dim(A⋊G) = |G|·dim A`;
    the orbit-Nakayama HH `[1,1,1,1]`/`[2,1,1,1]` (self-injective Nakayama literature); the
    `S₃`-on-`k³` `HH^0 = 2` (semisimple Morita `k×k`).
  - `oracle_crossengine`: **the R8 oracle** `stefan_decomposition.dims == HH^•(A⋊G) direct`
    (`Z/2` on dual numbers, `[1,1,1,1]`); `HH^•(A⋊G) direct == HH^•(orbit Nakayama)`; the FREE
    `Z/2`-swap orbit `HH^•(A⋊G) == HH^•(k[x]/(x²)) == [2,1,1]`; `Bimodule.regular ≡ ordinary`
    (P52 `M=A`); `HH^•(A,{}_σA)` twisted summand; the HOM `n=1` **per-summand split**
    `[{e:0},{σ:1}]` (nonzero twisted, findings 1a/2b); the `Z/3`/`GF(7)` non-involution direct.
  - `oracle_selfcert`: the dim-law certificate; `GroupAction.check` (each element an algebra
    automorphism; closure); trivial-`G` `T`-byte-identity; `presented_form` dim +
    multiplicativity certificate; the Reynolds idempotent `e² = e`; the identity-summand
    monomorphism `HH^•(A)^G ↪ HH^•(A⋊G)` (Marcos–MV, degreewise dim `≤`); the degree-0
    hand-pin `HH^0(A⋊G) = 1`; the **twist-symmetry** `HH^•(A,phi=g) == HH^•(A,psi=g⁻¹)` (2a);
    the canonical-key generator-order normalization.
  - `qpa`: QPA 1.37 has **no** skew-group-algebra HH-decomposition surface (live `NamesGVars()`
    guard that FAILS if that ever changes — the P35-products precedent). QPA *can* build a skew
    group algebra by hand and compute its `HH` dims — used as an **input-level** crosscheck of
    the DIRECT `HH^•(A⋊G)` on the small `Z/2`/`Z/3` examples (a construction anchor, not a
    decomposition compare).
- **Mid-merge-train counts.** Absolute suite counts drift; the verification task **recounts the
  oracle-class table at merge time** (`tests/release/test_oracle_classes.py`; paste live numbers,
  claim only this plan's deltas).
- Every merge updates `docs/verification.md` (new oracle rows + recounted class table green) and
  adds citations to `references.bib` + `registry.py` (`bibtex()` hard-fails if the two disagree).
  Conventional commits; green at every commit.

**No-red-commits discipline.** Each Task ships a self-contained green slice: the constructor
(Task 2) is usable before the decomposition (Tasks 3–4) exists; the GUI (Task G1) rides the
constructor.

---

## Tasks (TDD-shaped; write the failing test first)

### Task 0: citations (do first — later tasks reference the keys)

- [ ] Add to `references.bib` (BibTeX-VERIFIED only; delete any field that cannot be confirmed
      rather than guess — house rule):
```bibtex
@article{StefanHopfGalois1995,
  author = {{\c{S}}tefan, Drago{\c{s}}},
  title  = {Hochschild cohomology on {H}opf {G}alois extensions},
  journal = {J. Pure Appl. Algebra}, volume = {103}, number = {3},
  pages = {221--233}, year = {1995}}                       % # PIN vol/pages/DOI at citation-add
@article{SheplerWitherspoonGroupActions,
  author = {Shepler, Anne V. and Witherspoon, Sarah},
  title  = {Group actions on algebras and the graded {L}ie structure of {H}ochschild cohomology},
  journal = {J. Algebra}, volume = {351}, pages = {350--381}, year = {2012},
  note   = {arXiv:0911.0938}}      % venue FIXED this fix round: J. Algebra 351 (2012), NOT Adv. Math. (finding 3b)
@article{CibilsMarcosSmash2006,
  author = {Cibils, Claude and Marcos, Eduardo N.},
  title  = {Skew category, {G}alois covering and smash product of a $k$-category},
  journal = {Proc. Amer. Math. Soc.}, volume = {134}, number = {1},
  pages = {39--50}, year = {2006}, note = {arXiv:math/0312214}}   % # PIN pages at citation-add
@article{MarcosMartinezVillaMartinsInvariants,
  author = {Marcos, Eduardo N. and Mart{\'\i}nez-Villa, Roberto and Martins, Maria Izabel Ramalho},
  title  = {Hochschild cohomology of skew group rings and invariants},
  note   = {arXiv:math/0312460}, year = {2004}}            % # PIN journal if published
% NOT added as an anchor: arXiv:1804.02223 (Cibils-Marcos, WITHDRAWN) -- annotation only if at all.
```
      and in `registry.py` (`_r(key, bibtex_key, kind, title, annotation, *tags)`):
      `stefan_hopf_galois`, `shepler_witherspoon_group_actions`, `cibils_marcos_smash`,
      `marcos_mv_invariants`; reuse `witherspoon_gsm204` (survey) + `assem_book`.
- [ ] `.venv/bin/python -m pytest tests/citations/test_bib_structure.py -q` green.
- [ ] Commit `docs(citations): Plan 74 skew-group refs -- Stefan 1995 (added), Shepler-Witherspoon 0911.0938, Cibils-Marcos math/0312214; 1804.02223 flagged withdrawn/mis-attributed`.

### Task 1: `QuiverAutomorphism` + `GroupAction` (the action, certified)

- [ ] **Failing test** `tests/families/test_skew_group_action.py` (`oracle_selfcert`): the
      `σ: x ↦ −x` triple `(π=id, ρ=id, s={"x": -1})` on `truncated_polynomial(2)` has
      `.matrix(A) == [[1,0],[0,-1]]`; `.check(A) is True`; a **non**-automorphism (e.g.
      `s={"x": 2}` on `k[x]/(x²)`? — still an automorphism; use a genuine relation-breaker like a
      vertex-perm that does not preserve `I` on a larger example) raises loudly; `GroupAction.cyclic(σ, 2)`
      has `order == 2`, `conjugacy_classes == [[0],[1]]`, `centralizer(1) == [0,1]`.
- [ ] **Implement** the induced-automorphism matrix (paths → `φ(e_v)=e_{π v}`, `φ(α)=s(α)ρ(α)`,
      extended by the structure constants), the `check` (bijective, `φ(1)=1`,
      `φ(uv)=φ(u)φ(v)` on all `A.T`, equivalently preserves `I`), and the BFS closure
      (budget-capped, Cayley table, conjugacy classes, centralizers). Refuse floats in scalars.
- [ ] Green; commit `feat(families): QuiverAutomorphism + GroupAction with per-instance automorphism certificate (Plan 74)`.

### Task 2: `skew_group_algebra` constructor + dim law + trivial-G + presented view

- [ ] **Failing test** `tests/families/test_skew_group_algebra.py`:
      - (`oracle_selfcert`) `dim(skew_group_algebra(A, cyclic(σ,2))) == 2·dim A == 4`; `_validate`
        passes; the dim-law certificate holds; the multiplication table matches the hand table
        (`(x·1)(x·1)=0`, `(e·σ)(x·1) = −x·σ`, `(e·σ)²=e·1`).
      - (`oracle_crossengine`) **trivial-G byte-identity**: `skew_group_algebra(A,
        GroupAction.cyclic(id, 1)).T == A.T` and `HH^•` byte-identical.
      - (`oracle_selfcert`) `presented_form()` of the `Z/2`-dual-numbers `A⋊G` returns the
        2-vertex quiver with `relations == ["a1*a2", "a2*a1"]` (auto-recovered orbit Nakayama);
        `dim` certificate holds.
- [ ] **Implement** `skew_group_algebra` (structure constants from `(a·g)(b·h)=a·g(b)·gh`, basis
      labels `"<a>.<g>"`, `from_structure_constants(check=True)` + the dim-law assert),
      `Algebra.skew_group` delegate, `is_free_action`.
- [ ] Green; commit `feat(families): skew_group_algebra A⋊G constructor -- dim law, trivial-G identity, presented_form orbit-quiver view (Plan 74)`.

### Task 3: the twisted summand + `Z(g)`-action transport + invariants (the P52 bridge)

- [ ] **Failing test** `tests/hochschild/test_skew_group_summands.py`:
      - (`oracle_crossengine`) `Bimodule.twisted(A, phi=σ)` `.check()` passes;
        `HH^•(A, {}_σA) == [1,1,1,1]`; `Bimodule.regular(A) ≡ ordinary` `[2,1,1,1]` (P52 `M=A`).
      - (`oracle_selfcert`) the Reynolds idempotent `e = (1/|G|)Σ h` on `HH^0` satisfies `e²=e`;
        the **cohomology degree-0 hand-pins**: `HH^0(A,A)^G` dim `1`, `HH^0(A,{}_σA)^G` dim `0`.
      - (`oracle_crossengine` — **the bridge-non-triviality pins, findings 1a/2b/2e**) the
        **HOMOLOGY** per-summand splits: `HH_0` → identity coinv `1` **and** twisted coinv `1`;
        **`HH_1` → identity coinv `0` (σ = −1 on `HH_1 ≅ Ω¹ = k\,dx`) and twisted coinv `1`
        (nonzero, carries the degree)** — a dead/sign-wrong bridge fails this even though the total
        (`= 1`) is unchanged. Cross-anchored to the DIRECT `HH_•(A⋊G) = [2,1,1,1]`.
- [ ] **Implement** `_autom_action_on_classes(A, M, h, classes_data)` — the chain map on the
      twisted bar (co)chain complex induced by an automorphism `h` (transport
      `(h·f)(a_1⊗…⊗a_n) = h(f(h^{-1}a_1 ⊗ … ⊗ h^{-1}a_n))`, twisted per the `# PIN`ed
      convention), expressed on the P52 `classes` representatives via the shipped `coboundary`
      and `reduce`-to-class linear solves; and the Reynolds invariants `(-)^{Z(g)}` (average +
      rank). **`# PIN` the exact transport/twist against Shepler–Witherspoon before coding.**
- [ ] Green; commit `feat(hochschild): Z(g)-action transport on P52 twisted HH classes + Reynolds invariants (Plan 74 bridge)`.

### Task 4: `stefan_decomposition` assembly + the DIRECT cross-check + orbit oracle

- [ ] **Failing test** `tests/hochschild/test_skew_group_decomposition.py`:
      - (`oracle_crossengine`) **the R8 oracle**: `stefan_decomposition(A, cyclic(σ,2), 3,
        verify_direct=True)` has `dims == [1,1,1,1] == direct_dims` and `agrees is True`;
        `side="hom"` twin `dims == [2,1,1,1] == direct` with the **per-summand split pinned**
        (`n=0`: `[{e:1},{σ:1}]`; `n=1`: `[{e:0},{σ:1}]` — findings 1a/2b, NOT `[1,0]`).
      - (`oracle_crossengine`) **free-orbit oracle (finding 2d)**: the FREE `Z/2`-swap on
        `k(1⇄2)/rad²` — `is_free_action True`, `dim(A⋊G) = 8`, `presented_form == k[x]/(x²)`,
        `HH^•(A⋊G) == [2,1,1] == HH^•(k[x]/(x²))`, self-injective (over **QQ**; `GF(7)` shows the
        `char ≤ dim` `presented_form` refusal). The NON-free flagship is labelled `kG`-splitting.
      - (`oracle_crossengine`) **`Z/3` on `k[x]/(x³)` / `GF(7)`** (non-involution, finding 2a): the
        three twisted summands + identity; `HH^•(A⋊G)` DIRECT `[0..2] == [1,1,1]`, `HH_•` `[3,2,2]`
        with the `HH_0` split `1+1+1`; budget-capped past degree 2 (honest `status`).
      - (`oracle_crossengine` + `oracle_literature`) **`S₃` on `k³`** (non-abelian, finding 2c):
        `dim == 18` (dim law), `HH^0 == 2` (Morita `k×k`); the report states the 3-class shape
        (`Z(g) = S₃ / Z/2 / Z/3`); `HH^{≥1}` budget-capped (honest `status`).
      - (`oracle_selfcert`) the identity-summand monomorphism `dim HH^n(A)^G ≤ dim HH^n(A⋊G)`
        degreewise (Marcos–Martínez-Villa); the `"modular"` loud refusal on `GF(2)` (char `| |G|`).
- [ ] **Implement** the per-class assembly (rep per conjugacy class, `{}_gA`, `twisted_*_classes`,
      `_autom_action_on_classes` over `Z(g)`, invariants, `⊕`, **the per-class `summands` payload**),
      the char `∤ |G|` gate (loud `"modular"` otherwise), and the optional `verify_direct` DIRECT
      comparison.
- [ ] Green; commit `feat(hochschild): Stefan conjugacy-class HH decomposition of A⋊G, cross-checked vs direct + orbit-algebra oracle (Plan 74)`.

### Task 5: free-action detection + the honest covering slice

- [ ] **Failing test** (`oracle_selfcert`): `is_free_action` is `False` for the `Z/2`-dual-numbers
      action (σ fixes the single vertex) yet `A⋊G` is still an orbit Nakayama algebra (the
      splitting comes from `kG`, not from a free vertex action — a documented subtlety); `is_free_action`
      is `True` for a `Z/2` that **swaps** two vertices of a symmetric base (e.g. `Z/2` on
      `kA_2`-with-a-symmetric-orientation, picked at implementation). The covering HH transport is
      **NOT** claimed in v1 — assert the honest-scope note is present.
- [ ] **Implement** `is_free_action` (finite vertex-orbit check); record the covering-reduction
      HH transport as a **ledger deferral** (`DEEPER-ENGINES-BACKLOG.md`) — the Ştefan route
      already delivers HH within scope.
- [ ] Green; commit `feat(families): is_free_action detection + honest covering-slice scope note (Plan 74)`.

### Task G1: the `SkewGroupAlgebra` construction family + `skew_group_hh` kind (3 tiers, i18n ×4)

The construction family follows `TrivialExtension`/`OppositeAlgebra`/`SkewGentleAlgebra`
(no-code INPUT, `_build_synthetic`); the compute kind follows `tau_tilting`/`ar_quiver`
(algebra-only budget scalar). Once built, `A⋊G` flows through EVERY existing compute kind
(estimator sizes on the built algebra's dim `= |G|·dim A` automatically).

**Files:**
- `src/quiverlab/hpc/spec.py`: add `"SkewGroupAlgebra"` to `_SYNTHETIC_FAMILY_PARAMS`
  (`{"base", "generators", "orders"}` — `base` a Dynkin/family string, `generators` a list of
  quiver-automorphism dicts `{"vertex_perm", "arrow_perm", "arrow_scalars"}`, `orders` the
  generator orders), a `_build_skew_group(params, field)` builder, AND the `skew_group_hh` kind
  in `_dispatch` + grammar (budget suffix) + `_snip`.
- `docs/gui/runner.py`: the Pyodide twin — SAME builder + dispatch + `_snip` + ETA (byte-identical).
- `webapp/server/catalog.py`: mirror `SkewGroupAlgebra` in `_SYNTHETIC_FAMILIES` (flattened
  params + a prefill example: the `Z/2`-dual-numbers action); `webapp/server/schema.py`: the
  `skew_group_hh` grammar parse (THIRD site).
- `docs/gui/gui.js` + `webapp/static/gui/gui.js`: the construction-family panel entry (base
  picker + a small "group generator" editor: vertex/arrow permutation + arrow-scalar cells) and
  the `skew_group_hh` checkbox / `S.ids` / push-list / `renderBlock` / `scheduleProbe`.
- `webapp/templates/index.html`: one checkbox + the constructor panel row.
- `webapp/server/i18n/{en,es,fr,zh}.json`: **all four locales** — the family label + params + the
  `skew_group_hh` block key chain + `pick.kind.skew_group_hh` (exact key parity, gated by
  `tests/webapp/test_i18n.py`).
- `src/quiverlab/trace/results_html.py`: `_HEADINGS["skew_group_hh"]` + a render branch (the
  per-class decomposition table: class rep, `HH^n(A,{}_gA)` dims, invariant dims, the assembled
  total, and the DIRECT cross-check row with its `agrees` flag; honest `status`/`note`).
- `tests/webapp/_runner_goldens.json` + `test_runner_delegation.py`: ONE golden
  (`skew_group_hh_z2dual`); existing goldens byte-identical first.
- Tests: `tests/webapp/test_skew_group_kind_p74.py`, `tests/gui/test_skew_group_runner_twin_p74.py`,
  `tests/webapp/test_construction_families_input.py` (extend for `SkewGroupAlgebra`).

**Block shape** (`skew_group_hh_block(A_or_spec, budget)`):
```python
{"kind": "skew_group_hh", "base_dim": int, "group_order": int, "dim": int,   # = |G|*base_dim
 "char_ok": bool,                                                            # char ∤ |G|
 "decomposition": [{"class": str, "hh": [int,...], "inv": [int,...]}, ...] | None,
 "dims": [int,...], "direct_dims": [int,...] | None, "agrees": bool | None,
 "status": str, "note": str|None,
 "references": ["stefan_hopf_galois", "shepler_witherspoon_group_actions",
                "cibils_marcos_smash", "witherspoon_gsm204"]}
# modular (char | |G|) -> decomposition=None, status="modular", dims=direct only (constructor still builds)
# refusal (bad action) -> {"error": msg, "references": [...]}, never a 500.
```

- [ ] **Step 1** failing cross-runner tests (unmarked — extras-gated dir): the `SkewGroupAlgebra`
      prefill builds a dim-4 algebra; `skew_group_hh` block has `dims == [1,1,1,1]`,
      `direct_dims == [1,1,1,1]`, `agrees is True`, `"stefan_hopf_galois"` in citations; twin
      parity (`json.dumps(sort_keys=True)` across `hpc/spec.py` and `docs/gui/runner.py`).
- [ ] **Step 2** implement all three grammar sites + both builders + `gui.js` ×2 + `_snip` +
      `_HEADINGS` + render branch; canonical-key stability (family input canonicalizes unchanged;
      confirm an existing family request's key is byte-stable).
      - **Canonical-key generator-ORDER normalization (critic ask).** Mathematically-identical
        actions must not key differently. The `SkewGroupAlgebra` params carry a `generators` LIST,
        but the list order is presentation, not mathematics: `_build_skew_group` **canonicalizes
        before the key is taken** — each generator normalized to a stable form (sorted
        `vertex_perm`/`arrow_perm` maps + sorted `arrow_scalars`), then the generator list itself
        **sorted by that canonical form** (and de-duplicated). So `[σ, τ]` and `[τ, σ]` (same
        group, same action) collide to ONE canonical key. Add a test:
        two generator orderings of the same action produce the **same** `canonical_key` and the
        same golden. (Full up-to-iso group normalization — e.g. a different generating SET of the
        same group — is out of v1 scope; documented honestly, ledgered.)
- [ ] **Step 3** add ONE golden; note it in `test_runner_delegation.py`'s change-log docstring.
- [ ] **Step 4** run `tests/webapp/test_skew_group_kind_p74.py test_runner_delegation.py
      tests/gui/test_skew_group_runner_twin_p74.py tests/webapp/test_i18n.py
      tests/webapp/test_catalog.py tests/hpc -q`.
- [ ] **Step 5** commit `feat(gui,webapp,hpc): SkewGroupAlgebra construction family + skew_group_hh kind (Plan 74) -- both runners byte-identical, grammar in 3 sites, i18n x4, one golden`.

**i18n keys — ALL FOUR LOCALES** (mirror `pick.kind.tau_tilting` + the family label pattern):

| key | en | es | fr | zh |
|---|---|---|---|---|
| `family.SkewGroupAlgebra` | Skew group algebra A⋊G | Álgebra de grupo torcido A⋊G | Algèbre de groupe gauche A⋊G | 斜群代数 A⋊G |
| `pick.kind.skew_group_hh` | Skew-group HH decomposition | Descomposición HH de grupo torcido | Décomposition HH de groupe gauche | 斜群 HH 分解 |
| `block.skew_group_hh.title` | Hochschild HH of A⋊G (Ştefan decomposition) | HH de A⋊G (descomposición de Ştefan) | HH de A⋊G (décomposition de Ştefan) | A⋊G 的 HH（Ştefan 分解） |
| `block.skew_group_hh.class` | Conjugacy class | Clase de conjugación | Classe de conjugaison | 共轭类 |
| `block.skew_group_hh.invariants` | Z(g)-invariants | Invariantes de Z(g) | Invariants de Z(g) | Z(g) 不变量 |
| `block.skew_group_hh.direct` | Direct (cross-check) | Directo (verificación) | Direct (recoupement) | 直接（交叉验证） |
| `block.skew_group_hh.modular` | Modular case (char divides |G|): direct only | Caso modular (car divide |G|): solo directo | Cas modulaire (car divise |G|) : direct seul | 模情形（特征整除 |G|）：仅直接 |

### Task G2: verification page, README, suite gate, metaplan tick

**Files:** `docs/verification.md`; `README.md`; `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick
P74); existing release gates.

- [ ] **Verification rows** (new oracles): the dim law; trivial-G byte-identity; the R8
      direct-vs-decomposition oracle (`Z/2` on dual numbers `[1,1,1,1]`); the **HOM `n=1`
      per-summand split** `[{e:0},{σ:1}]` (bridge-non-triviality, findings 1a/2b); the FREE-orbit
      oracle (`Z/2`-swap → `k[x]/(x²)`, finding 2d); the non-abelian `S₃`-on-`k³` dim-law + `HH^0=2`
      (finding 2c); the orbit-Nakayama cross-check; the `Bimodule.regular ≡ ordinary` reuse (P52);
      the identity-summand monomorphism; the modular loud refusal. **Honest-scope entries
      (binding):** (a) the DIRECT computation is the trustworthy primary, the decomposition is
      cross-checked against it; (b) char `∤ |G|` for the decomposition/covering ONLY (constructor +
      direct are char-agnostic); (c) the presented `kQ'/I'` view is the **basic** algebra of `A⋊G`
      (Morita-equivalent; HH-faithful, dim-law is on the full algebra); (d) the covering-reduction
      HH transport is v1-deferred (detection + FREE-orbit oracle shipped; the flagship is
      `kG`-idempotent-splitting, NOT a free orbit — finding 2d); (e) QPA has no
      skew-group-HH-decomposition surface (input-level DIRECT crosscheck only); (f) **the twist
      convention is source-pinned from Shepler–Witherspoon and NOT summand-arbitrable** (order-2 AND
      order-3 give twist-symmetric summand dims — finding 2a; it is validated only through the
      assembled total + per-summand split); (g) deeper degrees / larger `|G|` are budget-capped,
      direct-oracle-pinned (dim-9 `Z/3` caps at degree 2, dim-18 `S₃` at degree 0); non-abelian `G`
      is IN scope (finding 2c) with proper-centralizer `(-)^{Z(g)}`; (h) `1804.02223` is withdrawn
      and mis-attributed in the source card (Cibils–Marcos, not CMRS) — noted, not anchored, Ştefan
      1995 / Shepler–Witherspoon 0911.0938 anchor instead; (i) **char `| |G|` (modular): the
      DECOMPOSITION refuses, and even the presented `A⋊G` is a GENUINELY DIFFERENT (modular) algebra
      than the semisimple-`kG` reading — the direct `A⋊G` is still computable but its
      `presented_form`/`is_selfinjective` describe the modular algebra** (finding 3c; distinct from
      the `char ≤ dim` `presented_form` refusal seen live on the dim-8 free example over `GF(7)`);
      (j) **the graded Lie / Gerstenhaber-bracket structure of `HH^•(A⋊G)` is OUT of scope —
      downstream of P70/P71/P51** on `A⋊G` as an input (the P70 seam; S–W 0911.0938 is about the
      bracket, this plan ships only dimensions).
- [ ] **Recount** the oracle-class table live at merge (`tests/release/test_oracle_classes.py`);
      paste live numbers, claim only this plan's deltas.
- [ ] **README** one line (representation-theory-first): "skew group algebras `A⋊G` as no-code
      inputs, with the Ştefan conjugacy-class Hochschild decomposition cross-checked against the
      direct engine".
- [ ] Gate chain (end with `; echo EXIT=$?`): `deep` (both kernel paths) + `fast` + `qpa` +
      `test_bib_structure.py` + `test_js_parses.py` + `test_runner_delegation.py` +
      `test_oracle_classes.py` + `mkdocs build --strict` green. Tick P74 in the metaplan ledger.
- [ ] Commit `docs(verification,readme): Plan 74 skew-group oracles + honest scope + recount; tick P74`.

---

## Acceptance (Plan-74 definition of done)

1. **Constructor `quiverlab.families.skew_group`:** `QuiverAutomorphism`, `GroupAction`,
   `skew_group_algebra`, `is_free_action`, `Algebra.skew_group` public. `GroupAction.check`
   certifies each element is an algebra automorphism preserving `I` and the generators close to
   the claimed finite group (loud otherwise). **HARD pins green:** dim law
   `dim(A⋊G) = |G|·dim A`; **trivial-G byte-identity** (`A⋊{1}.T == A.T`, HH byte-identical);
   `presented_form(A⋊G)` recovers the orbit `kQ'/I'` where basic + char-scoped
   (`Z/2`-dual-numbers → `1⇄2 / rad²`).
2. **Decomposition `quiverlab.hochschild.skew_group`:** `stefan_decomposition` public; char
   `∤ |G|` gate (loud `"modular"` otherwise). **The R8 oracle green:**
   `stefan_decomposition.dims == HH^•(A⋊G) direct == [1,1,1,1]` for `Z/2` on the dual numbers
   (both coh and hom sides). Each summand `HH^n(A,{}_gA)` comes from P52 (`Bimodule.twisted`); the
   `Z(g)`-action transport + Reynolds invariants are the new, self-certified bridge (Reynolds
   `e²=e`; identity-summand monomorphism `HH^n(A)^G ↪ HH^n(A⋊G)`). **The bridge-non-triviality
   pins green (findings 1a/2b):** the HOMOLOGY **per-summand split** `n=0 → [{e:1},{σ:1}]`,
   `n=1 → [{e:0},{σ:1}]` (a nonzero twisted invariant at `n≥1`, NOT `[1,0]`). **`Z/3` on
   `k[x]/(x³)`** (non-involution, `GF(7)`) and the **non-abelian `S₃` on `k³`** (dim `18`,
   `HH^0=2`, proper-centralizer `Z(g)`) both `verify_direct` where the budget allows, honest
   `status` where capped.
3. **Two invariants-free oracles green:** DIRECT `HH^•(A⋊G) == HH^•(orbit Nakayama)`
   (`[1,1,1,1]`/`[2,1,1,1]`, self-injective); `Bimodule.regular ≡ ordinary` (P52 `M=A`). **Plus the
   GENUINELY FREE orbit oracle (finding 2d):** the `Z/2`-swap on `k(1⇄2)/rad²` — `is_free_action
   True`, `dim 8`, `presented_form == k[x]/(x²)`, `HH^• == [2,1,1]`, self-injective (QQ) — with the
   flagship correctly labelled `kG`-idempotent-splitting (not a free orbit).
4. **Honest scope on the verification page (binding):** the ten entries (a)–(j) above, including
   the char-scope-is-on-the-decomposition refinement, the twist convention being source-pinned
   (NOT summand-arbitrable, 2a), non-abelian in scope (2c), the modular presented-algebra
   difference (3c), the P70/P71 bracket seam (j), and the `1804.02223`-withdrawn/mis-attributed
   finding.
5. **One construction family `SkewGroupAlgebra` + one kind `skew_group_hh`** clickable
   end-to-end (GUI canvas → build `A⋊G` → decomposition block → report) in **all four locales
   (en/es/fr/zh)** with `pick.kind.*` + family keys and exact key parity, both runners
   byte-identical, the grammar parsed in **all three sites** (`hpc/spec.py`, `docs/gui/runner.py`,
   `webapp/server/schema.py`), the family mirrored in the catalog, ONE golden with a documented
   change-log entry, canonical keys byte-stable (family input canonicalizes unchanged); the
   decomposition table + direct cross-check row render.
6. **QPA battery green (`-m qpa`):** the `NamesGVars()` guard (FAILS if QPA ships a skew-group-HH
   surface) + the input-level DIRECT `HH^•(A⋊G)` crosscheck on `Z/2`/`Z/3` small examples (build
   `A⋊G` by hand in QPA, compare dims).
7. `docs/verification.md` recounted (live numbers, mid-merge-train honest); citations added and
   BibTeX-verified — `stefan_hopf_galois` (ADDED, the missing anchor),
   `shepler_witherspoon_group_actions` (0911.0938), `cibils_marcos_smash` (Proc. AMS 134,
   math/0312214), `marcos_mv_invariants` (math/0312460), `witherspoon_gsm204` reused;
   `1804.02223` NOT anchored. README line; deep + fast + qpa + release + citations + mkdocs
   `--strict` green. Depends only on P52 (merged) — independent merge to `dev`.

---

## Methodology & assumptions

**Approach.** I (1) discovered the worktree was STALE (checked out at `b7fa566`/v0.3.0 — plans
only through P50), while the real base is `dev` at `26a8412` (P51–P68 merged, P52 among them);
git objects are shared, so I materialized the `dev` `src` tree into scratchpad
(`git archive dev src | tar`) and confirmed `quiverlab` imports from it with
`hochschild/coefficients.py` present. (2) Read the metaplan P74 card + §§1–4 and R8 verbatim, and
the plan-65/plan-63 house templates. (3) Read P52's `coefficients.py` in full and `core/algebra.py`
/`core/basic.py`/`families/basic.py`/`hpc/spec.py` for the exact live signatures. (4)
**Live-verified every numeric pin in the venv over QQ**, including building the `Z/2`-on-dual-numbers
`A⋊G` by hand as structure constants and computing BOTH sides of the decomposition oracle. (5)
Re-verified the citations via the web (arXiv abstracts), finding a real card defect. (6) Mirrored
the house form (records verbatim, reference re-verification with `# PIN`, live-facts table, scope
gates, Plan-32 markers, three-tier GUI, verification rows, TDD tasks).

**Assumptions (stated).** (a) The `dev` state I read is the true P74 base — verified by
`git ls-tree dev` and the P52 import; the deliverable is a new file that applies cleanly on top of
`dev` regardless of my stale worktree HEAD. (b) P52's `Bimodule.twisted(A, phi=g_matrix)` realizes
the Ştefan summand coefficient `{}_gA` — grounded by P52's own docstring naming P74 and by the
degree-0 hand-check matching. (c) The `Z(g)`-action-on-classes + Reynolds invariants is genuinely
absent from `dev` — grep found no skew-group/smash/invariants-of-HH code, and P52 exposes only
chain-level classes, not a group action on them. (d) The left-vs-right twist convention is pinned
from Shepler–Witherspoon at implementation — and it is **NOT summand-arbitrable** (fix round,
finding 2a: verified twist-symmetric summand dims on order-2 AND order-3 examples), so its
validation is the assembled total + per-summand split, not a summand comparison. (e) The GUI
construction-family pattern (`_SYNTHETIC_FAMILY_PARAMS`/`_build_synthetic` + catalog + twin) is the
right home for the constructor — evidenced by `TrivialExtension`/`OppositeAlgebra`/
`SkewGentleAlgebra` living exactly there.

**Card deviations (evidence-backed).** (1) **Char scope is on the DECOMPOSITION + covering, not
the constructor.** The card says "HARD scope char k ∤ |G|" flatly; but the smash-product structure
constants and the DIRECT `HH^•(A⋊G)` are char-agnostic (the char-free build `_validate`d over QQ,
and the smash formula never divides by `|G|`). Only Reynolds/Maschke — the decomposition and the
covering reduction — need `|G|` invertible. The refusal is scoped accordingly. (2) **The presented
`kQ'/I'` is reachable via the SHIPPED `presented_form()`, no bespoke Reiten–Riedtmann builder** —
live-verified it auto-recovers the orbit quiver `1⇄2 / rad²` from the structure-constant `A⋊G`
(caveat: it presents the *basic* algebra of `A⋊G`; the constructor still returns the full
structure-constant algebra so the dim law holds). (3) **The DIRECT computation is the primary
trustworthy route; the Ştefan decomposition is the value-add cross-checked against it** — a
deliberate inversion of a "spectral-sequence-first" reading, for honest ground truth (the SS
degenerates at E₂ for `|G|` invertible, so we compute the direct sum, not multi-page pages). (4)
**Citation defect fixed:** the card's "CMRS 1804.02223" is Cibils–Marcos (two authors) AND
withdrawn — replaced by the missing-but-load-bearing Ştefan 1995 + Shepler–Witherspoon 0911.0938 +
the published Cibils–Marcos math/0312214.

**Live-verified (venv, QQ) — observed values.** `A = k[x]/(x²)`, `G = Z/2`, `σ(x)=−x`:
`HH^•(A,A) = [2,1,1,1]`; `HH^•(A,A)` via `Bimodule.regular` `= [2,1,1,1]` (P52 `M=A`);
`HH^•(A,{}_σA) = [1,1,1,1]` and `HH_•(A,{}_σA) = [1,1,1,1]` (`Bimodule.twisted(A, phi=[[1,0],[0,-1]])`,
`.check()` passing); `HH^•(A⋊G)` DIRECT `= [1,1,1,1]`, `HH_•(A⋊G)` DIRECT `= [2,1,1,1]`;
`dim(A⋊G) = 4 = 2·2`; the orbit Nakayama `1⇄2 / rad²` `HH^• = [1,1,1,1]`, `HH_• = [2,1,1,1]`,
`is_selfinjective True`; `presented_form(A⋊G)` `→ 1⇄2, I'=⟨a1a2, a2a1⟩`. Degree-0 invariants
hand-verified: `HH^0(A,A)^G = 1`, `HH^0(A,{}_σA)^G = 0`, total `1 =` DIRECT.

**Fix-round live-verifications (venv; the critic's four majors + minors).** (2a, `GF(7)`) `Z/3` on
`k[x]/(x³)` (`σ(x)=2x`, `2³≡1`): `HH^•(A,phi=σ) = HH^•(A,psi=σ) = HH^•(A,phi=σ²) = [1,1,1,1,1]` —
the twist convention is NOT summand-arbitrable (nor is the involution case); direct
`HH^•(A⋊G) = [1,1,1]`, `HH_• = [3,2,2]` (top 2; `HH_0` splits `1+1+1`). (1a/2b) the decisive n≥1
witness: flagship HOM `n=1`, `HH_1(A,A) ≅ Ω¹ = k\,dx`, `σ(dx) = −dx` ⇒ identity coinv `0`, direct
`HH_1(A⋊G) = 1` ⇒ **twisted coinv `1`, nonzero** (confirmed by the minimal-resolution generator
sign `σ_P(ξ_1)=−ξ_1` independently). (2c, `GF(7)`) `S₃` on `k³`: `dim = 18` ✓, `HH^0 = 2`
(semisimple, Morita `k×k`), classes `{e}/3\text{ transp}/2\text{ 3-cyc}` with `Z(g)=S₃/Z/2/Z/3`.
(2d, QQ) FREE `Z/2`-swap on `k(1⇄2)/rad²`: `dim(A⋊G)=8`, `presented_form → k[x]/(x²)` (loop, `x²=0`),
`HH^• = [2,1,1]`, self-injective (over `GF(7)` the same build hits the `char ≤ dim` refusal — 3c).
(2e) flagship HOM `n=0`: identity coinv `1` + twisted coinv `1` = `2` = direct.

**Could NOT verify (labelled to-be-verified-at-implementation).** (1) The **general
`_autom_action_on_classes` machinery across all degrees/classes** — that IS the plan's Task 3; but
the fix round CLOSED the "dead bridge" gap: the flagship HOM `n=0`/`n=1` per-summand splits are
hand-derived two ways (Kähler + minimal-resolution sign) and pin a NONZERO twisted invariant at
`n=1`. (2) The **twist convention** (left vs right) — `# PIN`ed to Shepler–Witherspoon; **NOT
summand-arbitrable** (finding 2a), validated by the assembled split. (3) **Deeper degrees / larger
algebras** — dim-9 `Z/3` capped at degree 2, dim-18 `S₃` at degree 0 (`max_cells`); TDD-verified
against the direct oracle at implementation, budget-capped. (4) The **modular refusal wording** —
the mechanism (Maschke/Reynolds failure) is certain; the exact `QuiverlabError` string is written
at implementation. (5) The **Ştefan JPAA vol/pages** are `# PIN`ed against the source PDF; the
**Shepler–Witherspoon venue is FIXED this round to J. Algebra 351 (2012), 350–381** (was wrongly
stubbed "Adv. Math.", finding 3b).

**Deliberately not checked, and why.** I did not run any pytest suite (other agents run tests
concurrently in shared worktrees — I kept probes light; the flagship HH pins were computed
directly). I did not implement `_autom_action_on_classes` (that is the plan's Task 3, not a
spec-time probe) — instead I proved the oracle is satisfiable by hand at degree 0 and by the
DIRECT engine at all degrees. I did not resolve the `people.tamu.edu/~sjw/pub/ring.pdf` survey to
a published title (`# PIN`; `witherspoon_gsm204` already ships as the survey anchor). I did not
push deeper than degree 3 on `A⋊G` (the 4-dim bar complex at degree 5 exceeded a 120 s probe box —
the deep-degree behaviour is TDD-verified against the direct oracle, not at spec time).

**Why I believe the result is correct.** The R8 oracle is satisfied by an independent, robust
ground truth (the DIRECT engine on the constructed algebra) with a numeric match `[1,1,1,1]`, and
is corroborated two further ways that need no new machinery: the orbit-Nakayama cross-check
(same HH, independently built) and P52's `M=A` reuse. The summand machinery already exists and is
named for P74 by P52's own author; the only genuinely new slice (the `Z(g)`-action + Reynolds
invariants) is isolated in Task 3, is **hand-verified at homology degrees 0 AND 1 (a nonzero
twisted summand at `n=1`, two independent derivations — closing the fix round's "dead bridge"
worry)**, and is pinned at every degree by the direct oracle + the per-summand split. Every scope
boundary the record names becomes a loud typed refusal, and the card deviations (char scope on the
decomposition; twist convention source-pinned not summand-arbitrable; non-abelian in scope) are
each evidenced by a live computation.

## Open design risks (what a critic will likely attack)

1. **The `Z(g)`-action transport formula (Task 3), incl. the resolution SIGN.** The induced chain
   map on the twisted bar (co)chain complex — and crucially the **sign/scalar on the resolution
   generators** — decides which summand carries a class; a sign error gives the right TOTAL but the
   wrong per-summand split. *Mitigation:* the fix round pins the **per-summand split** as an oracle
   (flagship HOM `n=1` `= [{e:0},{σ:1}]`, hand-derived two ways), not only the total; the Reynolds
   `e²=e` + the identity-summand monomorphism bracket it; the DIRECT oracle pins totals everywhere.
2. **Non-abelian `G` / non-singleton conjugacy classes.** *Resolved (fix round):* the `S₃`-on-`k³`
   pin (dim law `18`, `HH^0=2`, classes `{e}/3\text{ transp}/2\text{ 3-cyc}` with `Z(g)=S₃/Z/2/Z/3`)
   is live-verified and non-abelian is IN scope; the proper-centralizer `(-)^{Z(g)}` path is
   exercised. *Residual:* the dim-18 direct `HH^{≥1}` is `max_cells`-capped — honest `status`; the
   deep non-abelian HH is direct-oracle-pinned only where the budget allows (else the base-`A`
   summand route, cheaper).
3. **Non-basic `A⋊G` and the presented view.** `presented_form` returns the *basic* algebra's
   quiver; a critic may conflate that with `A⋊G` itself. *Mitigation:* the constructor always
   returns the full structure-constant algebra (dim law), and the honest-scope entry (c) states
   the presented view is Morita-equivalent (HH-faithful) but not dimension-faithful.
4. **Char scope wording.** Claiming the constructor is char-agnostic while the decomposition is
   gated must be air-tight. *Mitigation:* the structure-constant build was `_validate`d over QQ
   with no `|G|`-division anywhere; the `"modular"` refusal is only on `stefan_decomposition`.
5. **The withdrawn `1804.02223`.** A reviewer re-checking the card will hit the same withdrawn
   preprint. *Mitigation:* the reference re-verification section documents it explicitly and
   anchors the math on Ştefan 1995 + Shepler–Witherspoon instead — nothing is transcribed from
   `1804.02223`.
6. **Cost / budget.** `A⋊G` bar HH grows with `|G|·dim A`; deep degrees or `|G|≥3` on a
   non-trivial base can be heavy. *Mitigation:* the decomposition rides the base-`A` bar complex
   (cheaper than the `|G|·dim A` full algebra) for the summands, and `verify_direct` is opt-in;
   honest `status="budget"` on the cap; the GUI kind is budget-carrying like `tau_tilting`.

## Change log

- **2026-08-08 authoring.** Initial plan (R8 skew group algebras `A⋊G`). Discovered the stale
  worktree, based the design on the `dev` P52 state (materialized for live-verification).
  Constructor (quiver-automorphism action, per-instance certificate, dim law, trivial-G identity,
  `presented_form` orbit-quiver view) + the Ştefan conjugacy-class HH decomposition (P52 twisted
  summands + a new `Z(g)`-action/Reynolds-invariants bridge) cross-checked against the DIRECT
  engine. Flagship `Z/2`-on-dual-numbers oracle live-verified BOTH sides (`[1,1,1,1]`), degree-0
  invariants hand-verified, orbit-Nakayama + `M=A` corroborations live. Citations re-verified:
  **Ştefan 1995 added** (the missing anchor); **card defect fixed** — `1804.02223` is
  Cibils–Marcos (two authors) and WITHDRAWN, replaced by the published anchors. Card deviations
  recorded (char scope on the decomposition not the constructor; direct-as-primary; presented via
  the shipped `presented_form`). GUI: `SkewGroupAlgebra` construction family + `skew_group_hh`
  kind, three tiers, i18n ×4, one golden.
- **2026-08-08 fix round (adversarial review NEEDS WORK → all findings adjudicated valid, applied).**
  (2a) The twist convention is **NOT summand-arbitrable** — live-verified order-2 AND order-3
  (`Z/3` on `k[x]/(x³)` over `GF(7)`) both give twist-symmetric summand dims (`{}_φA_ψ ≅
  {}_1A_{φ⁻¹ψ}`); removed the false "direct oracle arbitrates a wrong convention at summand level"
  claim, re-scoped to source-pinned + assembled-split validation. (1a/2b) Added the decisive
  **nonzero twisted-summand witness at `n≥1`**: flagship HOMOLOGY `n=1`, `HH_1(A,A)≅Ω¹` with
  `σ=−1`, so identity coinv `0`, twisted coinv `1` (carries the degree) — hand-derived two ways
  (Kähler + minimal-resolution sign); added the **per-summand split** oracle beyond `total==direct`
  + the identity-summand monomorphism + QPA input-level anchor. (2c) Added the **non-abelian `S₃`
  on `k³`** pin (dim `18`, `HH^0=2`, 3-class shape) — non-abelian stays in scope. (3a) Anchored the
  additive formula primarily on Shepler–Witherspoon (Ştefan for the SS origin; body not fetched).
  (2d) Added the **genuinely FREE** `Z/2`-swap orbit oracle (`→ k[x]/(x²)`), relabelled the
  flagship "`kG`-idempotent splitting, not a free orbit". (2e) Hand-verified the homology degree-0
  split. (3b) **Fixed S–W venue: J. Algebra 351 (2012), 350–381** (was "Adv. Math."). (3c)
  Foregrounded the char-`|`-`|G|` modular presented-algebra difference + the char-`≤`-dim
  `presented_form` refusal (seen live). Added the canonical-key **generator-order normalization**
  and the **P70/P71 seam** (bracket/Lie downstream). All new pins live-verified in the venv.
