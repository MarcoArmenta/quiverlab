# Plan 72: Split-extension Hochschild LES + certified arrow removal/addition (P72 / R5 + R6)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps
> use checkbox (`- [ ]`) syntax for tracking. Do the reference re-verification in the
> Record / Reference sections BEFORE writing any oracle pin — every `# PIN` below must be
> resolved against the sources (not this plan's prose) at citation-add / test-freeze time.
> The two load-bearing theorems (CMRS Thm 4.1 connecting-map cup formula; CLMS Thm 3.2 / 4.1
> arrow-removal isomorphism + correction) must be transcribed **verbatim with equation
> numbers** from the PDFs before the code that depends on them is written.

**Goal.** Two self-verifying Hochschild reductions, each grounded on machinery quiverlab
already ships (Plan 52 HH-with-coefficients merged `77f6cc9`; Plan 31 `TrivialExtension`
presented build; Plan 35 products for the connecting-map arbiter):

- **(a — R5) The split-extension / trivial-extension Hochschild long exact sequence
  (Cibils–Marcos–Redondo–Solotar, `math/0102194`).** For the split (square-zero) extension
  `L = B ⋉ M` (`B` a subalgebra, `M` a two-sided ideal with `M·M = 0`, `L = B ⊕ M` as
  `k`-spaces), the short exact sequence of `L`-bimodules `0 → M → L → B → 0` induces the
  **long exact sequence in `HH^•(L, −)`**
  ```
  ··· → HH^n(L,M) →^{ι} HH^n(L,L) →^{π} HH^n(L,B) →^{δ^n} HH^{n+1}(L,M) → ···
  ```
  whose middle term `HH^n(L,L) = HH^n(L)` is the ordinary Hochschild cohomology of the
  extension. **Degree-range boundary rule (MAJOR, live-verified below):** to assemble
  `HH^n(L)` for `n = 0..top` the LES needs `δ^n : HH^n(L,B) → HH^{n+1}(L,M)` for **every**
  `n ≤ top`, so the **`M`-flank must be computed to degree `top+1`** (the `B`-flank to
  `top`); truncating the `M`-flank at `top` silently drops the `HH^{n+1}(L,M)` term into
  which `δ^top` lands (it happened to be `0` for `T(kA₂)`, but `HH^5(T(kD₄),M) = 4 ≠ 0` —
  see the flank tables). The **connecting map `δ`** is CMRS **Theorem 4.1**'s cup-product
  with the identity `1_M`; the two flanks decompose (CMRS **Corollary 3.2**, the tensor-power
  grading) `HH^n(L,X) = ⊕_{p+q=n} Ext^q_{B^e}(M^{⊗_B p}, X)` **when `M` is one-sided
  projective as a `B`-module** (Cor 3.2's hypothesis — see the honest-scope note below),
  whose **`p = 0` leading piece is `HH^n(B)` (for `X = B`) / `H^n(B, M)` (for `X = M`)** —
  the card's "from `HH^*(B) + H^*(B,M)`", surfaced as the `n = 0` identity + the summand
  inequality (not a full splitting for `T(kA_n)`/`T(kD₄)`, whose `M = D(B)` is NOT
  bimodule-projective). The **flagship instance is the trivial extension**
  `L = T(B) = B ⋉ D(B)` (`M = D(B) = Hom_k(B,k)`), the `TrivialExtension(B)` quiverlab
  already builds (Plan 31). R5's stated value: *"turns the existing brute `T(B)` computation
  into a structured, self-verifying decomposition."* Plus the **theorem `HH^1(B ⋉ M) ≠ 0`
  for `M ≠ 0`** (CMRS **Theorem 5.5**): the **grading derivation** of `L = B ⊕ M` (`M` in
  degree 1; CMRS **Def. 5.10**), `E(b+m) = m`, is an **outer** derivation, so `k ↪ HH^1(L)`
  always — realized as the grading-derivation witness class; for a **directed**
  (acyclic-quiver) `B` the sharper `HH^1(T(B)) = k ⊕ HH^1(B)` holds (CMRS's "one-way" case).

- **(b — R6) Certified arrow removal / addition Hochschild reductions (Cibils–Lanzilotta–
  Marcos–Solotar, `1812.07655`, Proc. AMS 148 (2020) 2421–2432).** An **inert arrow** (CLMS
  **Def. 3.1**) is an arrow of `Q` that appears in **no** minimal relation of `I`. Deleting a
  set `D` of inert arrows produces `B = A ∖ D = kQ'/I'` (`Q' = Q ∖ D`, `I' = I` verbatim —
  no relation mentions `D`), and CLMS give a **clean homology isomorphism**
  `HH_n(A) ≅ HH_n(B)` for **`n ≥ 2`** (**Thm 3.2**) — the certificate asserts this range
  ONLY; `HH_0`, `HH_1` are reported informationally (in every battery instance below they
  ALSO agree, and `HH_0 = A/[A,A]` is **provably** unchanged — an inert arrow never lies on
  an oriented cycle of a finite-dimensional `kQ/I`, so no cyclic-path class moves). By
  contrast **cohomology carries a `Ext` correction** `dim HH^n(B) = dim HH^n(A) − dim Ext^n_A(⋯)`
  (**Thm 4.2**, degree range and exact `Ext` term `# PIN` from the PROC PDF) that can be
  **nonzero for `n ≥ 2`** (live: `HH^2` differs on pair P2), while the `n = 0, 1` cohomology
  deltas are a **separate center/disconnection effect** (`B` may be disconnected — `HH^0 =
  Z` grows — a phenomenon distinct from the `Ext` term); the plan reports the `n ≥ 2`
  correction and the `n = 0, 1` deltas under **separate honest labels**. The asymmetry
  (homology clean `n ≥ 2`, cohomology corrected) is contractual. The dual **arrow addition**
  (**Thm 3.5 / 3.6**) builds `A = B_F` as the tensor algebra `T_B(N)`, finite-dimensional
  **iff** adding `F` creates **no relative cycle**, with the same `HH_{≥2}` isomorphism.
  P72 owns the **reduction primitive**; **P73 (Han recognizer + Jacobi–Zariski, R7) owns the
  recognizer / transport** built on top — the seam is `inert_arrows` / `remove_arrows` /
  `add_arrows` as reusable callables, no Han logic in P72. **Shared-substrate note (P73
  friction):** P73 specs a `arrow_removal_subalgebra(A, F)` with *different, more general
  semantics* — arbitrary arrows `F` with the reduced ideal `I_B = ker(π)` (the induced
  presentation), NOT the inert-only `I' = I`-verbatim case P72 needs. The two constructors
  are the **same substrate**: whichever plan lands first **owns the module**; the second
  **consumes** it (and P72's `remove_arrows` is the inert special case
  `arrow_removal_subalgebra(A, inert F)` where `I_B = I`). A **cross-consistency test**
  (`remove_arrows ≡ arrow_removal_subalgebra` on inert `F`) is a task in whichever merges
  second; the shared-file merge friction is flagged here and in the P73 doc.

- **The seam between (a) and (b).** Both replace a "big" Hochschild computation by a
  structured one certified against the direct answer. They share the P52 coefficient engine
  (a: `HH^•(L, M)`; b: the `Ext` correction) and the same honest contract: the reduction is
  **complete iff its scope gate holds**, else a loud typed refusal — never a guessed value.

**GUI.** TWO new **algebra-only, budget-carrying** compute kinds (the `tau_tilting` /
`ar_quiver` precedent — expensive scalar kinds, NOT the fast `recognizers` block), each
served by **all three tiers** (`hpc/spec.py::_dispatch` + the byte-identical Pyodide twin
`docs/gui/runner.py` + the `webapp/server/schema.py` grammar), with a budget
(`split_extension:6` / `arrow_removal`), the recognizer-panel rows, **all four locales
(en / es / fr / zh) with exact key parity**, canonical-key stability, report rendering
(`trace/results_html.py`), and citations:

- **`split_extension`** — on the drawn algebra `A` (interpreted as `B`), form
  `L = T(B) = B ⋉ D(B)`, render the **structured LES decomposition**: the flanks
  `HH^•(L, D(B))` and `HH^•(L, B)` with their `p = 0` leading pieces `HH^•(B)` / `H^•(B, D(B))`
  named, the connecting-map ranks, the assembled `HH^•(L)`, the direct `HH^•(L)` cross-check
  row, and the `HH^1(L) ≠ 0` grading-derivation witness. (The card's *"structured
  decomposition rendering on trivial-extension inputs."*)
- **`arrow_removal`** — auto-detect the inert arrows of `A`, build `B = A ∖ (inert)`, render
  the certificate: the inert-arrow list with the "appears in no relation" witness, the
  `HH_n(A) = HH_n(B)` (`n ≥ 2`) homology-isomorphism table, and the cohomology `Ext`-correction
  row.

**Architecture.** TWO new source modules + one small P52 constructor + two algebra-only kinds;
everything else is a thin exact layer over primitives already on `dev`:

- **`src/quiverlab/hochschild/split_extension.py`** (new — part a). Public:
  - `split_extension(B, M=None, *, name=None) -> Algebra` — the extension `L = B ⋉ M`. With
    `M = None` (default) returns `TrivialExtension(B)` (`M = D(B)`, presented, Plan 31). A
    general `M` (a P52 `Bimodule` over `B`) is accepted **only when the resulting `L` can be
    presented** (has a quiver) — otherwise a loud refusal (the inflation of coefficients
    needs `L.quiver`; DD-A1). Self-cert: `dim L = dim B + dim M`; for `M = D(B)`, `L`
    byte-identical to `TrivialExtension(B)`.
  - `split_extension_cohomology(B, top, *, M=None, max_cells=4_000_000) -> SplitExtReport` —
    the LES assembly + the direct cross-check.
  - `split_extension_homology(B, top, *, M=None, max_cells=...) -> SplitExtReport` — the
    homology twin (the dual SES `0 → M → L → B → 0` in `HH_•(L, −)`).
  - `hh1_grading_witness(B, *, M=None) -> dict` — the `HH^1(L) ≠ 0` witness (the grading
    derivation as an explicit cocycle + its nonzero class), and the directed-`B` sharpening
    `HH^1(L) = k ⊕ HH^1(B)`.
- **`src/quiverlab/hochschild/arrow_removal.py`** (new — part b). Public:
  - `inert_arrows(A) -> list[str]` — the arrows of `Q` in no minimal relation of `I` (Def
    3.1), read off `A.relations` (each `Relation.terms` = `(coeff, word-of-arrow-names)`).
  - `remove_arrows(A, arrows) -> Algebra` — `B = A ∖ arrows` (loud unless every named arrow
    is inert); `Q' = Q ∖ arrows`, `I' = I` verbatim; certified `dim` drop.
  - `add_arrows(A, new_arrows) -> Algebra` — the dual `A_F = T_A(N)` (Thm 3.5/3.6); loud
    refusal when `F` creates a relative cycle (would be infinite-dimensional).
  - `arrow_removal(A, arrows=None, top=..., *, side="both") -> ArrowRemovalReport` — the
    certified HH reduction: build `B`, the `HH_{≥2}` homology isomorphism (Thm 3.2), the
    cohomology `Ext`-correction (Thm 4.1/4.2), the direct cross-check. `arrows=None`
    auto-detects and removes ALL inert arrows (the maximal reduction).
- **`src/quiverlab/hochschild/coefficients.py`** (extend P52 — ONE constructor):
  `Bimodule.inflate(L, B, coeff_over_B, *, vanishing_arrows) -> Bimodule` — the `L`-bimodule
  obtained from a `B`-bimodule `coeff_over_B` by **inflation along the split projection**
  `π: L ↠ B`: the generators of `L` that come from `B` (vertices + `B`'s arrows) act as in
  `coeff_over_B`; the ideal generators `vanishing_arrows` (the new `M`-part arrows) act as
  **zero** (`M·M = 0`). Implemented over `Bimodule.from_actions` (folds to the full basis,
  `.check()` certifies the `L^e`-axioms). Reused by P74 (skew group) — noted for the record.
- **`Algebra` delegates** (`core/algebra.py`, thin lazy-import beside the P52
  `hochschild_cohomology(..., coefficients=)`):
  `A.split_extension_cohomology(top, M=None)`, `A.split_extension_homology(top, M=None)`,
  `A.arrow_removal(arrows=None, top=6)`, `A.inert_arrows()`.
- The two kinds wire into `hpc/spec.py::_dispatch` (the `if kind == …` chain beside
  `tau_tilting` / `ar_quiver` at `spec.py:1567`+) + `docs/gui/runner.py` (byte-identical twin)
  + `webapp/server/schema.py` (the grammar parse — a THIRD site) + the GUI touchpoints
  (checkbox / `S.ids` / push-list / `renderBlock` / `scheduleProbe`) + the i18n chains ×4 +
  `trace/results_html.py` — exactly following the `ar_quiver` / `tau_tilting` algebra-only
  kinds.

**Tech Stack.** Pure exact-`Domain` linear algebra over the P52 bar cochain complexes +
exact combinatorics on `A.relations`. **No floats in `src/`** (AST-gated by
`tests/test_no_floats.py`): dims/ranks are `int`, verdicts `bool`/`str`, connecting maps are
`Domain` matrices, the assembled dim sequences are `list[int]`. All homology exact
(`Algebra.hochschild_cohomology`/`_homology` with `coefficients=`, `nullspace`/`rank`/`solve`
from `fields.linalg`). Composition left-to-right (`a*b` = first `a` then `b`).

---

## Records (verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R5 — Split-algebra / trivial-extension HH long exact sequence.** [B-scout P7; keep]
> Object: HH^*(B⋉M) assembled from HH^*(B) + H^*(B,M) via the LES with cup-product connecting
> map; theorem HH¹(B⋉M) ≠ 0 for M ≠ 0 (grading-derivation witness; state M ≠ 0). Ref:
> Cibils–Marcos–Redondo–Solotar arXiv:math/0102194. Turns the existing brute T(B) computation
> into a structured, self-verifying decomposition. Oracles: LES-vs-direct degreewise on
> T(kA_n)/T(kD₄) (in suite); HH¹≠0 across the zoo. Size M. Deps: R4.

> **R6 — Arrow removal/addition HH reductions (Han accelerator).** [A-critic promotion of a
> demoted item] Object: certified finite HH reductions along deleting/adding arrows of kQ/I.
> Ref: Cibils–Lanzilotta–Marcos–Solotar arXiv:1812.07655. Pairs with R7's Jacobi–Zariski
> machinery; feeds the backlogged Han's campaigns. Size M.

Metaplan card (`docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P72, Wave 4, tier δ, **needs
P52**):

> **P72 — split-extension LES + arrow removal (R5+R6, size M+M; needs P52).** (a) HH^*(B⋉M)
> from HH^*(B) + H^*(B,M) via the LES with cup connecting map; HH¹(B⋉M) ≠ 0 for M ≠ 0
> (grading-derivation witness). CMRS math/0102194. (b) Certified arrow removal/addition HH
> reductions (CLMS 1812.07655), feeding Han's-conjecture campaigns. Oracles: LES-vs-direct on
> T(kA_n)/T(kD₄); HH¹ ≠ 0 zoo sweep; reduction ≡ direct on removed-arrow pairs. GUI: structured
> decomposition rendering on trivial-extension inputs.

---

## Reference re-verification (done at authoring; findings binding, `# PIN` = resolve at citation-add)

The metaplan standing rule (§1.5) requires the plan writer to re-verify citations. Findings:

1. **Cibils, C., Marcos, E., Redondo, M. J. & Solotar, A. — "Cohomology of split algebras and
   of trivial extensions".** arXiv:`math/0102194`. **ALREADY REGISTERED** as `cmrs_split`
   (bibtex `cibilsmarcosredondosolotar2003`, annotated "HH^1(T(A)) is never zero … the
   trivial-extension HH^1 oracle") — **reuse, no new key.** `# PIN`: the abstract (fetched
   this authoring round) confirms the paper *"describe[s] a long exact sequence computing the
   Hochschild cohomology of L … study the connecting homomorphism using the cup-product …
   the first Hochschild cohomology group of a trivial extension never vanishes"* — the exact
   **LES**, **Theorem 4.1** (connecting map `δ^{p,q}φ = 1_M ⌣ φ + (−1)^{p+q+1} φ ⌣ 1_M`),
   the tensor-power **decomposition Corollary 3.2** `H^n(L,X) = ⊕_{p+q=n} Ext^q_{B^e}(M^{⊗_B p}, X)`
   (its **hypothesis** — `M` one-sided `B`-projective — recorded: it FAILS for `M = D(B)`
   over a non-self-injective `B`, so for `T(kA_n)`/`T(kD₄)` only the `p = 0` leading identity
   + summand inequality survive, NOT the full splitting), the **Theorem 5.5 `HH^1(TA) ≠ 0`**
   statement (with the `HH^1(TA) = k ⊕ HH^1(A)` one-way case), and the grading-derivation
   **Definition 5.10** are transcribed **verbatim with theorem numbers from the published PDF**
   BEFORE the assembly / witness code uses them (the LES sign convention and the direction of
   `δ` are load-bearing).
   Confirm the venue **Glasgow Math. J. 45 (2003)** and page range against the journal
   landing page at merge (the bibtex is present; the venue string is `# PIN`).
2. **Cibils, C., Redondo, M. J. & Saorín, M. — "The first cohomology group of the trivial
   extension of a monomial algebra".** **ALREADY REGISTERED** as `crs_trivial_ext_hh1`
   (bibtex `cibilsredondosaorin2004`, "the trivial-extension first-cohomology oracle").
   **Reuse** as the corroborating reference for the `HH^1(T(B))` decomposition / the directed
   sharpening. `# PIN`: verify it states the `HH^1(T(A))` decomposition the plan pins as the
   "one-way" formula (or scope the pin to CMRS Thm 5.5 alone if CRS-2004's hypotheses differ).
3. **Cibils, C., Lanzilotta, M., Marcos, E. N. & Solotar, A. — "Deleting or adding arrows of
   a bound quiver algebra and Hochschild (co)homology".** arXiv:`1812.07655`; **published
   Proc. Amer. Math. Soc. 148 (2020), no. 6, 2421–2432** (verified this authoring round via
   the AMS / CONICET / HAL listings). **NEW key required** — `clms_arrow_removal` (bibtex
   `cibilslanzilottamarcossolotar2020`). `# PIN`: **Def. 3.1** (inert arrow), **Thm 3.2**
   (`HH_∗(A ∖ D) ≅ HH_∗(A)` for `∗ ≥ 2`), **Thm 4.1 / 4.2** (the cohomology `Ext`-correction
   formula and its exact degree range), **Thm 3.5 / 3.6** (arrow addition = `T_B(N)`, the
   relative-cycle finiteness criterion) transcribed **verbatim with theorem numbers from the
   PROC PDF** BEFORE any oracle pin — the fetched paraphrase of the correction formula
   `dim HH^∗(B) = dim HH^∗(A) − dim Ext^∗_A((fA)', Ae)` is **NOT trustworthy for the exact
   formula/degree range** (the plan's live data below fixes the *phenomenology*, the PDF fixes
   the *statement*). The paper's own tools — relative `HH` (P52 territory), **Kaygun's
   Jacobi–Zariski LES**, and the one-step relative projective resolution of a tensor algebra —
   are the exact machinery P73 will consume; name them but do not implement JZ here.
4. **Han, Y. — "Hochschild homology and Han's conjecture".** **ALREADY REGISTERED** as
   `han_conjecture` (bibtex `Han2006`). **Reuse** only in the annotation naming R6's downstream
   purpose (the seam to P73); no HH-of-Han logic in P72.
5. **Green, E. L., Psaroudakis, C. & Schroll, S. — "arrow removal" (singular equivalence).**
   OPTIONAL corroboration for the inert-arrow/`A ∖ ⟨a⟩` construction. `# PIN`: add an `@misc`/
   `@article` ONLY if cited in an annotation; not load-bearing for any pin (CLMS is the HH
   anchor). Verify authorship/venue before adding, or drop (never guess).

**Already shipped and reused (verified live at spec time):** P52 `hochschild/coefficients.py`
(`Bimodule.regular/dual/twisted/twisted_by_nakayama/mod_socle/from_actions`,
`Algebra.hochschild_cohomology(top, coefficients=M)` / `hochschild_homology(...)`, the bar
route forced by a non-`None` coefficient); Plan 31 `families.TrivialExtension` (presented
`kQ_T/I_T`, `dim = 2·dim B`); Plan 35 `hochschild/products.py` (cup — the connecting-map
theoretical arbiter, see the honest-scope note); `combinat.Quiver(...).algebra(relations=…)`;
`Algebra.relations` (`list[Relation]`, `Relation.terms = ((coeff, (arrow, …)), …)`).
**New keys to add:** `clms_arrow_removal` (arXiv:1812.07655 / Proc. AMS 148 (2020) 2421–2432).

---

## Mathematical foundation (the definitions, exactly — the plan's ground truth)

Throughout `B = kQ/I` is basic, connected, finite-dimensional over an exact `Domain`; `L`
denotes a split square-zero extension `L = B ⋉ M` (`M` a `B`-bimodule, `M·M = 0`,
`L = B ⊕ M` as `k`-spaces, multiplication `(b+m)(b'+m') = bb' + (bm' + mb')`); the flagship is
`M = D(B) = Hom_k(B,k)` with `L = T(B) = TrivialExtension(B)`. `HH^n(L, X) = Ext^n_{L^e}(L, X)`
is Hochschild cohomology with bimodule coefficients (P52); `HH^n(L) = HH^n(L, L)`. Composition
is left-to-right.

### (a) The CMRS long exact sequence — what is COMPUTED vs. dimension-bookkeeping

- **The tool.** `0 → M → L → B → 0` is a short exact sequence of `L`-bimodules (`M` an ideal,
  `B = L/M`). The bar resolution `bar_•(L)` is `L^e`-free, so `Hom_{L^e}(bar_•(L), −)` is
  **exact**: applying it gives a **short exact sequence of cochain complexes**
  `0 → C^•(L, M) → C^•(L, L) → C^•(L, B) → 0`, hence the long exact sequence
  ```
  ··· → HH^n(L,M) →^{ι^n} HH^n(L,L) →^{π^n} HH^n(L,B) →^{δ^n} HH^{n+1}(L,M) → ···   (LES)
  ```
- **DEGREE-RANGE BOUNDARY RULE (MAJOR — the flank contract).** Assembling `HH^n(L)` for
  `n = 0..top` reads `δ^{n-1}` and `δ^n`, so `δ^{top} : HH^{top}(L,B) → HH^{top+1}(L,M)`
  is required — the **`M`-flank runs to degree `top+1`, the `B`-flank to `top`, the direct
  middle to `top`**. Truncating the `M`-flank at `top` drops the `HH^{top+1}(L,M)` codomain
  of `δ^{top}` and makes the top-degree assembly **unverifiable** (it silently worked on
  `T(kA₂)` only because `HH^5(L,M) = 0`; for `T(kD₄)` at `top = 4`,
  `HH^4(L,B) = 4 ≠ 0` maps into `HH^5(L,M) = 4` — see the flank tables). `Bimodule.inflate`
  is called once; `flank_M = HH^{0..top+1}(L, M_L)`, `flank_B = HH^{0..top}(L, B_L)`.
- **Which coefficients are COMPUTED.** The **flanks** `HH^n(L, M)` and `HH^n(L, B)` are
  computed **directly by the P52 engine** — `M` and `B` are built as `L`-bimodules by
  **inflation along `π: L ↠ B`** (`Bimodule.inflate`): the `B`-generators of `L` act as `M`
  resp. the regular `B`-action, and the new `M`-part arrows act as **zero** (because
  `M·M = 0`). The middle term `HH^n(L, L) = HH^n(L)` is the ordinary P52/bar computation
  (`coefficients=None`). All three are genuine cocycle-space dimensions — **COMPUTED**.
- **The connecting map `δ` — COMPUTED via the snake, IDENTIFIED with the cup (Thm 4.1).**
  `δ^n : HH^n(L,B) → HH^{n+1}(L,M)` is the **snake-lemma connecting homomorphism** of (LES):
  lift a cocycle `φ ∈ Z^n(L,B)` to `φ̃ ∈ C^n(L,L)` (any preimage under the surjection
  `C^n(L,L) ↠ C^n(L,B)`), apply the `L`-coboundary `δ_L φ̃ ∈ C^{n+1}(L,L)`, which lies in the
  image of `C^{n+1}(L,M)` (since `π(δ_L φ̃) = δ_B φ = 0`); its `M`-preimage's class is `δ^n[φ]`.
  This is **pure exact linear algebra over the P52 bar coboundary matrices** — no new engine.
  CMRS **Theorem 4.1** identifies this same connecting map with the **cup product**
  `δ^{p,q}φ = 1_M ⌣ φ + (−1)^{p+q+1} φ ⌣ 1_M` — the card's "cup connecting map"; the snake
  realization IS that cup (Yoneda product with the extension class), and the plan **records
  Thm 4.1 as the theoretical identity** (`# PIN` verbatim). **Honest-scope (DD-A2):** a
  *literal* independent recomputation of `δ` through the Plan-35 cup surface is **out of
  scope** — Plan 35's cup is coefficient-`A` only (`HH^•(A) ⊗ HH^•(A) → HH^•(A)`), whereas
  Thm 4.1's cup pairs *different* coefficient bimodules (`Hom(M,M)`, `B`, `M`); building that
  coefficient-cup is not this plan's work. The cup identity is therefore a **cited theorem +
  a documented deferral**, and the connecting map is validated **by the exactness self-cert
  and the assembled-≡-direct cross-check** (below), not by a second cup computation.
- **How `HH^n(L)` is ASSEMBLED (and why assembled-≡-direct is a genuine oracle, not a
  tautology).** From (LES) exactness,
  ```
  dim HH^n(L) = [dim HH^n(L,M) − rank δ^{n-1}] + [dim HH^n(L,B) − rank δ^n].
  ```
  The **left-hand `HH^n(L)`** is computed independently by the ordinary bar engine (the
  MIDDLE complex `C^•(L,L)`); the **right-hand side** is assembled from the SUB and QUOTIENT
  complexes `C^•(L,M)`, `C^•(L,B)` and the snake `δ` (which lifts through the middle complex).
  Bug-free, they must agree by exactness — so the check is an **integration `oracle_selfcert`
  (exactness holds) + `oracle_crossengine`** (the sub/quotient-complex route reproduces the
  standalone middle-complex HH): a defect in the inflation, the snake lift, the rank, or the
  degreewise SES wiring breaks it loudly. **This is R5's "self-verifying decomposition."**
- **The card's "`HH^*(B) + H^*(B,M)`" — the `p = 0` leading piece (identity + INEQUALITY,
  NOT a splitting).** CMRS's decomposition `HH^n(L,X) = ⊕_{p+q=n} Ext^q_{B^e}(M^{⊗_B p}, X)`
  (**Cor. 3.2**) is a genuine direct-sum splitting only under its **hypothesis** that `M` is
  **one-sided projective as a `B`-module**. For the flagship `M = D(B)`, `D(B)` is
  bimodule-projective **iff `B` is self-injective** — which `kA_n` / `kD₄` are **NOT** — so
  the full `p`-graded splitting **does not hold** for `T(kA_n)`/`T(kD₄)`, and pinning
  `flank_B == HH^•(B)` degreewise would be **FALSE** (and, as `leading_B` is defined as
  `HH^•(B)`, self-comparing — a tautology). What the plan pins instead, both true and
  non-vacuous:
  - **the `n = 0` identity** `flank_B[0] == HH^0(B)` and `flank_M[0] == H^0(B,M)` (the `p = 0`
    term is all of degree 0; live-verified — `T(kA₂)`: `flank_B[0]=1=HH^0(kA₂)`,
    `flank_M[0]=2=H^0(kA₂,D(kA₂))`);
  - **the summand inequality** `flank_B[n] ≥ HH^n(B)` and `flank_M[n] ≥ H^n(B,M)` for all `n`
    (the `p = 0` term is always a subquotient contributing at least its dimension — live-
    verified, and **tight and non-vacuous** on `T(2`-Kronecker`)`: `HH^1(B) = 3 = flank_B[1]`,
    while `T(kA₂)` is strict, `flank_B = [1,0,1,2,1,0] ≥ HH^•(kA₂) = [1,0,0,0,0,0]`).
  The plan **computes `HH^•(B)` and `H^•(B, M)` independently over `B`** (P52 over the smaller
  algebra) and surfaces them as informational leading terms; the **`p ≥ 1` tensor-power terms
  are captured automatically inside the `L`-flanks** `HH^•(L, ±)` (the engine sees all of `L`);
  a **standalone `⊗_B`-tensor-power + bimodule-`Ext` engine realizing the graded decomposition
  term-by-term is an explicit honest-scope deferral** (DD-A3) — the LES-over-`L` route needs
  none of it and is card-complete.
- **`HH^1(L) ≠ 0` (CMRS Thm 5.5) — the grading-derivation witness.** `L = B ⊕ M` is `ℤ`-graded
  (`B` in degree 0, `M` in degree 1); the **grading derivation** `E: L → L`, `E(b+m) = m`
  (`E|_B = 0`, `E|_M = id`), is a `k`-derivation. It is **outer** (`M ≠ 0` ⇒ `E` is not
  inner, since an inner derivation `[l,−]` restricted to the grading is trivial on the center
  contribution), so `0 ≠ [E] ∈ HH^1(L)` — hence `HH^1(B ⋉ M) ≠ 0` whenever `M ≠ 0`. The plan
  materializes `E` as an explicit degree-1 bar cocycle, certifies `[E] ≠ 0` (its class is
  nonzero in `HH^1(L)`), and for **directed `B`** (acyclic `Q`) pins the sharper
  `HH^1(T(B)) = k ⊕ HH^1(B)` (CMRS one-way case) — `[E]` spans the `k`.

### (b) Arrow removal / addition (CLMS `1812.07655`)

- **Inert arrow (Def. 3.1).** An arrow `a ∈ Q_1` is **inert** iff it appears in **no** word
  of any minimal relation generating `I = ⟨R⟩`. Decidable from the presentation:
  ```
  inert_arrows(A) = { a ∈ Q_1 : a ∉ ⋃_{rel ∈ A.relations} ⋃_{(c,w) ∈ rel.terms} set(w) }.
  ```
  (`A.relations` is the reduced/Gröbner presentation of `I`; `# PIN` — confirm minimality of
  that generating set for non-monomial `I` at implementation, else compute the minimal
  generators from the Gröbner basis.)
- **Deletion `B = A ∖ D` (Thm 3.2).** For `D` a set of inert arrows, `B = kQ'/I'` with
  `Q' = Q ∖ D` and `I' = I` (every relation survives verbatim — none mentions `D`). Then
  **`HH_n(A) ≅ HH_n(B)` for all `n ≥ 2`** (homology; **clean**). The certificate `hom_agrees`
  asserts this range **only**. The `n = 0, 1` homology is **reported informationally** (never
  asserted), because Thm 3.2 does not cover it — but:
  - **`HH_0` is PROVABLY unchanged.** `HH_0(A) = A/[A,A]` has a basis of the vertices plus the
    cyclic (loop-based) paths mod commutators. An inert arrow `a: i → j` in a **finite-
    dimensional** `kQ/I` **cannot lie on an oriented cycle** — a cycle through `a` forces its
    cyclic powers to vanish, i.e. to enter `I`, but every such relation would contain a letter
    of the cycle including `a`, contradicting inertness (or, with `i ≠ j`, `a = [e_i, a] ∈
    [A,A]` already). So no cyclic-path class involves `a`; deleting `a` (and the paths through
    it, all non-cyclic) leaves `A/[A,A]` fixed. Hence **inert removal never changes `HH_0`**
    (live: `HH_0` agrees on every battery pair, incl. the disconnecting P1/P3).
  - **`HH_1`** may in principle differ; in **every** battery pair below it also agrees. The
    plan therefore states the contract honestly (Thm 3.2 asserts `n ≥ 2`; `HH_0` provably
    fixed; `HH_1` reported, observed-equal on the batteries) and does NOT manufacture a
    convenience claim — a genuine `HH_{0,1}` homology divergence, if it exists, is outside the
    inert-removal regime this plan constructs, and constructing one is flagged as an open
    boundary probe (Task B2 Step 3).
- **Cohomology correction (Thm 4.2, `n ≥ 2`) vs. the low-degree center/disconnection deltas
  (`n = 0, 1`).** For `n ≥ 2`, `dim HH^n(B) = dim HH^n(A) − dim Ext^n_A(⋯)` (**Thm 4.2**; the
  exact `Ext` term and its degree range `# PIN` from the PDF), so `coh_correction[n] :=
  coh_A[n] − coh_B[n] = dim Ext^n ≥ 0` for `n ≥ 2` — **nonzero on P2 (`n = 2`: `1`)**, so
  **cohomology arrow-removal is NOT a clean isomorphism**. The `n = 0, 1` cohomology deltas
  are a **DIFFERENT phenomenon** — a center/disconnection effect (`B` may be disconnected, so
  `HH^0(B) = Z(B)` grows; e.g. P1 `HH^0(A) = 1`, `HH^0(B) = 3`, delta `−2` from the isolated
  block, NOT an `Ext` term). The report carries `coh_correction` (`n ≥ 2`, the `Ext` term)
  and `coh_low_delta` (`n = 0, 1`, labelled "center/disconnection") as **separate fields**,
  and never conflates them.
- **Addition `A = B_F` (Thm 3.5 / 3.6).** Adding arrows `F` gives the tensor algebra
  `A = T_B(N)` (`N` = the `B`-bimodule spanned by `F`); **finite-dimensional iff `F` creates
  no relative cycle** (`s(a_2) B t(a_1) ≠ 0` links `(a_2, a_1)`; a cycle of links ⇒ infinite
  dim). Same `HH_{n ≥ 2}` isomorphism. Loud refusal on a relative cycle.
- **The homology/cohomology asymmetry is design-load-bearing.** `arrow_removal` reports BOTH:
  the clean `HH_{≥2}` homology isomorphism (the primary certificate) AND the cohomology
  `Ext`-correction row (honest — the direct `HH^•(A)` recovered from `HH^•(B)` + correction).

### The two consequences the plan leans on

1. **Self-verification, not acceleration (R5).** The split-extension LES over `L` costs the
   same order as the direct `T(B)` computation (both run bar over `L`); its value is
   **structure + a self-cert against the direct answer** (R5's own framing). The card's
   acceleration reading ("from the smaller `HH^*(B)`") is realized only as the surfaced
   `p = 0` leading pieces; the genuine `⊗_B` acceleration is deferred (DD-A3). Arrow removal
   (R6) **is** a genuine reduction (`B` smaller than `A`).
2. **Scope gates make refusals instant and loud.** Presentation-less `L`/`A` ⇒ loud (the
   inflation and inert-detection need the quiver). Non-inert arrow ⇒ loud with the witnessing
   relation. Relative cycle on addition ⇒ loud. `char ≤ dim` over `GF(p)` ⇒ the shipped P52 /
   `TrivialExtension` primitives raise (the socle/dual/presented builds inherit their caveats).

---

## Live-verified facts (all recomputed in the venv at spec time, QQ, `quiverlab.__file__` =
## the `dev` editable install; top = 4 unless noted)

**HH pins — the LES flanks + trivial extension (part a).**

| `B` | `dim B` | `HH^•(B)` | `H^•(B, D(B))` | `dim T(B)` | `HH^•(T(B))` | `HH^1(T(B))` |
|---|---|---|---|---|---|---|
| `kA₂` (`1→2`) | 3 | `[1,0,0,0,0]` | `[2,0,0,0,0]` | 6 | `[3,1,1,1,1]` | **1** |
| `kA₃` (`1→2→3`) | 6 | `[1,0,0,0,0]` | `[3,0,0,0,0]` | 12 | `[4,1,1,1,1]` | **1** |
| `kD₄` (subspace `2,3,4→1`) | 7 | `[1,0,0,0,0]` | — | 14 | `[5,1,0,0,2]` | **1** |

Homology twin (`kA₂`): `HH_•(B) = [2,0,0,0,0]`, `H_•(B,D(B)) = [1,0,0,0,0]`,
`HH_•(T(B)) = [3,1,1,1,1]`. (`kA₃`): `HH_•(B) = [3,0,0,0,0]`, `HH_•(T(B)) = [4,1,1,1,1]`.

**The `L`-coefficient LES — flanks (`M` to `top+1`, `B` to `top`, per the boundary rule) +
dimensional consistency (part a); inflated `L`-bimodules `M_L = D(B)↑`, `B_L = B↑`, both
`.check()` = `True`. Timings measured this fix round (flanks FORCE the bar route; the direct
middle over a non-monomial presented `T(B)` routes via the auto/CS engine — a DIFFERENT engine,
so the direct timing does not bound the flank timing):**

`L = T(kA₂)`, `top = 5` (M-flank to 6):

| `n` | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| `HH^n(L,L)` **direct** | 3 | 1 | 1 | 1 | 1 | 1 | — |
| `HH^n(L, M_L)` (M-flank → `top+1`) | 2 | 1 | 0 | 1 | 2 | 1 | 0 |
| `HH^n(L, B_L)` (B-flank → `top`) | 1 | 0 | 1 | 2 | 1 | 0 | — |
| implied `rank δ^n` | 0 | 0 | 0 | 2 | 0 | 0 | — |

`L = T(kA₃)`, `top = 4` (build 0.12 s, flanks 0.48 s):

| `n` | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| `HH^n(L,L)` **direct** | 4 | 1 | 1 | 1 | 1 | — |
| `HH^n(L, M_L)` (→ `top+1=5`) | 3 | 1 | 0 | 1 | 0 | 1 |
| `HH^n(L, B_L)` (→ `top=4`) | 1 | 0 | 1 | 0 | 1 | — |
| implied `rank δ^n` | 0 | 0 | 0 | 0 | 0 | — |

`L = T(kD₄)`, `top = 4` (build 0.19 s, flanks ~2.9–3.0 s) — **the pin that FAILS on truncated
flanks:** `HH^4(L,B_L) = 4 ≠ 0` maps by `δ^4` into `HH^5(L,M_L) = 4`, so assembling `HH^4(L)=2`
is UNVERIFIABLE without the `top+1` M-flank:

| `n` | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| `HH^n(L,L)` **direct** | 5 | 1 | 0 | 0 | 2 | — |
| `HH^n(L, M_L)` (→ `top+1=5`) | 4 | 1 | 0 | 0 | 1 | **4** |
| `HH^n(L, B_L)` (→ `top=4`) | 1 | 0 | 0 | 1 | **4** | — |
| implied `rank δ^n` | 0 | 0 | 0 | 1 | 2 | — |

- The implied connecting-map ranks are **all within `0 ≤ rank δ^n ≤ min(HH^n(L,B), HH^{n+1}(L,M))`**
  on ALL three (checked with the `top+1` M-flank) — the LES is dimensionally consistent, a
  strong de-risking of the snake-`δ` + assembly. `# PIN` — the implementation replaces the
  *implied* ranks with the *computed* snake ranks and asserts the identity degreewise.
- **`p = 0` leading piece — the identity + the (non-vacuous) inequality:** `flank_B[0] = HH^0(B)`
  and `flank_M[0] = H^0(B,M)` on all three (`T(kA₂)`: `1 = HH^0(kA₂)`, `2 = H^0(kA₂,D(kA₂))`);
  and `flank_B[n] ≥ HH^n(B)` for all `n` — **tight** on `T(2`-Kronecker`)`:
  `HH^•(B) = [1,3,0,0]`, `flank_B = [1,3,3,5]`, so `flank_B[1] = 3 = HH^1(B)` (the summand
  inequality is a genuine, non-tautological constraint), and strict on `T(kA₂)`
  (`flank_B = [1,0,1,2,1,0] > HH^•(kA₂) = [1,0,0,0,0,0]` at `n ≥ 2`).

**`HH^1(T(B)) ≠ 0` sweep (part a), and the directed sharpening `HH^1(T(B)) = 1 + HH^1(B)`:**

| `B` | directed? | `HH^1(B)` | `HH^1(T(B))` | `1 + HH^1(B)` | formula holds? | `≠ 0`? |
|---|---|---|---|---|---|---|
| `kA₂` | yes | 0 | 1 | 1 | ✓ | ✓ |
| `kA₃` | yes | 0 | 1 | 1 | ✓ | ✓ |
| `kD₄` | yes | 0 | 1 | 1 | ✓ | ✓ |
| `2`-Kronecker (`1⇉2`) | yes | 3 | 4 | 4 | ✓ | ✓ |
| `k[x]/(x²)` | **no** (loop) | 1 | 4 | 2 | N/A (not directed) | ✓ |
| `k[x]/(x³)` | **no** (loop) | 2 | 7 | 3 | N/A (not directed) | ✓ |

`HH^1(T(B)) > 0` in **every** case (the grading-derivation `k` always injects); the sharper
`= k ⊕ HH^1(B)` holds exactly on the **directed** (`m`-Kronecker included) `B` and is
honestly N/A for the non-directed `k[x]/(xᵃ)`.

**Arrow removal (part b) — `reduction ≡ direct` on removed-arrow pairs, QQ, top = 5.**

| pair | `A = kQ/I` (inert arrow) | `dim A` | `B = A ∖ D` | `dim B` | `inert_arrows(A)` |
|---|---|---|---|---|---|
| P1 | loop `x` at 1 (`x²=0`) + `c:1→2` | 5 | `k[x]/(x²) × k` | 3 | `['c']` |
| P2 | `1→2→3` (`ab=0`) + `c:1→3` | 6 | `kA₃/(ab)` | 5 | `['c']` |
| P3 | loops `x` at 1, `y` at 2 (`x²=y²=0`) + `c:1→2` | 8 | `k[x]/(x²) × k[y]/(y²)` | 4 | `['c']` |
| P4 (binomial, non-monomial) | comm. square `a c − b d` + inert `e:1→4` | 10 | comm. square (dim 9) | 9 | `['e']` |

- **P1:** `HH_•(A) = HH_•(B) = [3,1,1,1,1,1]` — **homology iso in all degrees** (`n ≥ 2`
  **Thm 3.2 ✓**; `HH_0 = HH_1` agree — `HH_0` provably fixed). Cohomology `HH^•(A) =
  [1,1,1,1,1,1]`, `HH^•(B) = [3,1,1,1,1,1]`: differ **only at `n = 0`** (`coh_low_delta[0] = −2`,
  the **center/disconnection** effect — `B` is disconnected, `Z(B) = 3`); the `n ≥ 2` `Ext`
  correction is `0`.
- **P2:** `HH_•(A) = HH_•(B) = [3,0,0,0,0,0]` (**Thm 3.2 ✓**; `HH_{0,1}` agree). Cohomology
  `HH^•(A) = [1,1,1,0,0,0]`, `HH^•(B) = [1,0,0,0,0,0]`: `coh_low_delta[1] = 1` (`n = 1`) AND
  `coh_correction[2] = 1` (`n = 2`, the **Thm 4.2 `Ext` term, nonzero for `n ≥ 2`**) — the
  homology/cohomology asymmetry made concrete.
- **P3 (the cohomology low-degree boundary):** `HH_•(A) = HH_•(B) = [4,2,2,2,2,2]` — homology
  iso in **all** degrees (`HH_0 = 4` fixed despite `B` disconnecting into two loop-blocks).
  Cohomology `HH^•(A) = [1,3,2,2,2,2]`, `HH^•(B) = [4,2,2,2,2,2]`: differ at `n = 0`
  (`coh_low_delta[0] = −3`, disconnection) AND `n = 1` (`coh_low_delta[1] = +1`, the bridge
  `c` contributes an extra outer derivation), while `n ≥ 2` **agree** (`coh_correction ≡ 0`) —
  a clean instance of the `n = 0,1` low-degree deltas being distinct from the `n ≥ 2` `Ext`
  term.
- **P4 (non-monomial inert detection):** the binomial relation `a c − b d` (a `Relation` with
  two terms) is scanned correctly — `inert_arrows(A) = ['e']` (`a,b,c,d` all appear;
  `e` does not), `dim A = 10`, `dim B = 9`, `HH_•(A) = HH_•(B) = [4,0,0,0,0,0]` (hereditary-
  directed, `n ≥ 2` trivially ✓). This pins the non-monomial `# PIN` (Task B1 Step 1).

**What I could NOT live-verify (labelled to-be-verified-at-implementation):**
- The **snake connecting map `δ`** and its equality with CMRS Thm 4.1's cup — I verified the
  LES is *dimensionally* consistent (the implied ranks fit, on `T(kA₂)`/`T(kA₃)`/`T(kD₄)` with
  the `top+1` M-flank), not the map itself. Building the snake + asserting the degreewise rank
  identity is the plan's work (Task A2).
- **`L`-coefficient flanks re-measured (this fix round) — feasibility now empirical, not
  argued:** the general `inflate` helper ran the flanks (bar route, forced by coefficients)
  for `T(kA₂)`/`T(kA₃)`/`T(kD₄)`; timings build+flanks ≈ `T(kA₃)` 0.6 s, `T(kD₄)` ~3.0 s
  (dominated by the 5-degree M-flank). NOTE the earlier "direct `T(kD₄)` ≈ 1 s ⇒ flanks
  feasible" reasoning was a **non-sequitur** (the direct route is CS/auto over the non-monomial
  presented `T(B)`; the flanks force the bar route — different engines); the measured flank
  timings above are the real feasibility evidence.
- A genuine **`HH_{0,1}` HOMOLOGY divergence** under inert removal — NOT found (four batteries
  P1–P4 all agree in homology at every degree; `HH_0` is provably fixed). The plan states the
  contract from the theorem (`n ≥ 2`, Thm 3.2) and reports `HH_{0,1}` informationally; the
  open boundary probe is Task B2 Step 3.
- The **exact CLMS cohomology-correction `Ext` formula / degree range** — pinned by the
  *phenomenology* (P2 `n = 2`; P3 `n = 1` low-degree) but the *statement* is `# PIN` (Thm 4.2)
  from the PROC PDF.

---

## Scope gates (contractual; every boundary is a loud typed refusal, `QuiverlabError`)

| surface | scope | refusal |
|---|---|---|
| `split_extension(B, M)` / `split_extension_cohomology` / `_homology` | `B` **quiver-presented** and the extension `L` **presentable** (`M = D(B)` always is, via `TrivialExtension`; a general `M` only when a presented `L` is supplied/derivable) | presentation-less `B`, or a general `M` whose `L` cannot be presented ⇒ `QuiverlabError` ("split-extension LES needs a quiver presentation of B ⋉ M to inflate coefficients") |
| the inflated flanks `HH^•(L, M_L)`, `HH^•(L, B_L)` | P52 coefficient scope (any exact Domain via the bar route) | inherits P52's refusals (e.g. structure-constant-only coefficient over GF(p) for path-basis needs) |
| `hh1_grading_witness` directed sharpening `HH^1(T(B)) = k ⊕ HH^1(B)` | `B` **directed** (acyclic `Q`) | non-directed `B` ⇒ the sharper formula is **not asserted** (only `HH^1(L) ≠ 0` is), honestly labelled — never a wrong `k ⊕ HH^1(B)` claim |
| `inert_arrows` / `remove_arrows` / `arrow_removal` | `A` **quiver-presented** | structure-constant-only ⇒ loud (no relations to read Def 3.1) |
| `remove_arrows(A, D)` | every arrow in `D` **inert** (Def 3.1) | a non-inert arrow ⇒ loud with the **witnessing relation** ("arrow `a` appears in relation `r` — not inert; deletion would change the ideal") |
| `arrow_removal` cohomology | reports the **`Ext`-correction** (Thm 4.1/4.2), never a clean iso | — (the homology `HH_{≥2}` iso is the only iso claim) |
| `add_arrows(A, F)` | `F` creates **no relative cycle** (Thm 3.6 finiteness) | a relative cycle ⇒ loud ("adding `F` creates a relative cycle — `B_F` is infinite-dimensional") |
| both, over `GF(p)` | **char 0 or char > dim** where the socle/dual/presented `TrivialExtension` build needs it | `char ≤ dim` ⇒ loud from the shipped P31/P52 primitives (never a silent wrong verdict) |

All batteries run over **QQ**; a `GF(32003)` parity spot-check is kept for the distinct-dim
`kA_n` family (large prime ⇒ char > dim, no `TrivialExtension`/socle raise).

---

## API surface (public via `import quiverlab`; exact only)

```python
# --- part (a): split-extension LES (src/quiverlab/hochschild/split_extension.py) ---
def split_extension(B, M=None, *, name=None)              -> Algebra     # L = B (x) M; M=None -> T(B)
def split_extension_cohomology(B, top, *, M=None, max_cells=4_000_000) -> SplitExtReport
def split_extension_homology(B, top, *, M=None, max_cells=4_000_000)   -> SplitExtReport
def hh1_grading_witness(B, *, M=None)                     -> dict        # the HH^1(L)!=0 witness

@dataclass
class SplitExtReport:
    base: object              # B
    extension: object         # L = B (x) M
    coeff_name: str           # "D(B)" | the Bimodule.describe()
    side: str                 # "cohomology" | "homology"
    n_range: int              # top
    flank_M: list[int]        # dim HH^n(L, M)   (n = 0..top+1  -- BOUNDARY RULE: to top+1)
    flank_B: list[int]        # dim HH^n(L, B)   (n = 0..top)
    delta_ranks: list[int]    # rank delta^n     (computed snake ranks, n = 0..top)
    assembled: list[int]      # dim HH^n(L) from the LES  (n = 0..top)
    direct: list[int] | None  # dim HH^n(L) directly (the cross-check; None if skipped)
    leading_B: list[int]      # HH^n(B)  (informational; the p=0 term -- oracle is the
                              #   flank_B[0]==HH^0(B) identity + flank_B[n] >= HH^n(B) inequality,
                              #   NOT equality: D(B) is bimodule-projective only iff B self-injective)
    leading_M: list[int]      # H^n(B,M) (informational; same identity + inequality contract)
    hh1_nonzero: bool         # HH^1(L) != 0 (the grading-derivation witness)
    hh1_directed_formula: bool | None  # HH^1(L) == 1 + HH^1(B) (directed B; None if not directed)
    exact: bool               # the LES exactness self-cert (rank alternating-sum) holds
    agrees: bool | None       # assembled == direct
    status: str               # "complete" | "unsupported" | "budget" | "error"
    note: str

# --- part (b): arrow removal / addition (src/quiverlab/hochschild/arrow_removal.py) ---
def inert_arrows(A)                                       -> list[str]
def remove_arrows(A, arrows)                              -> Algebra     # B = A \ arrows
def add_arrows(A, new_arrows)                             -> Algebra     # A_F = T_A(N)
def arrow_removal(A, arrows=None, top=6, *, side="both")  -> ArrowRemovalReport

@dataclass
class ArrowRemovalReport:
    algebra: object           # A
    reduced: object           # B = A \ D
    removed: list[str]        # D (the inert arrows deleted; = inert_arrows(A) if arrows=None)
    inert_certificate: dict   # {arrow: "appears in no relation"} witnesses (Def 3.1)
    hom_A: list[int]          # HH_n(A)
    hom_B: list[int]          # HH_n(B)
    hom_iso_from: int         # 2 (HH_n iso for n >= 2, Thm 3.2 -- the ONLY asserted range)
    hom_agrees: bool          # HH_n(A) == HH_n(B) for n >= 2   (n=0,1 NOT asserted)
    hom_low_agrees: list[bool]  # [HH_0 same, HH_1 same]  (informational; HH_0 provably True)
    coh_A: list[int] | None   # HH^n(A)
    coh_B: list[int] | None   # HH^n(B)
    coh_correction: list[int] | None  # (HH^n(A) - HH^n(B)) for n >= 2  (the Thm 4.2 Ext term)
    coh_low_delta: list[int] | None   # (HH^n(A) - HH^n(B)) for n = 0,1  (center/disconnection)
    status: str               # "complete" | "unsupported" | "error"
    note: str

# --- Bimodule constructor (extend src/quiverlab/hochschild/coefficients.py) ---
@classmethod
def Bimodule.inflate(cls, L, B, coeff_over_B, *, vanishing_arrows) -> "Bimodule"

# Algebra delegates (core/algebra.py, thin lazy-import beside coefficients-HH):
Algebra.split_extension_cohomology(top, M=None)   -> SplitExtReport
Algebra.split_extension_homology(top, M=None)     -> SplitExtReport
Algebra.arrow_removal(arrows=None, top=6)         -> ArrowRemovalReport
Algebra.inert_arrows()                            -> list[str]
```

`SplitExtReport.__bool__` / `ArrowRemovalReport.__bool__` are NOT defined (data reports, not
predicates) — mirroring `HHProducts` / `ARQuiver`, whose truthiness is never relied on.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **P52 (HH with coefficients, merged `77f6cc9`), P31 (`TrivialExtension` presented), P35
  (products — the connecting-map arbiter, honest-scope-deferred), P01 (`Quiver.algebra`) are
  merged on `dev` — verified by live import at spec time.** This plan consumes, with the exact
  signatures verified in `dev`:
  - P52 `hochschild/coefficients.py`: `Bimodule(algebra, dim_M, Lact, Ract, name)`,
    `.regular(B)`, `.dual(B)`, `.from_actions(A, dim_M, left_maps, right_maps, name)` (folds
    generator-indexed action maps to the full basis, `.check()` certifies the `A^e`-axioms),
    `.invariants_dim()` / `.coinvariants_dim()` (`HH^0`/`HH_0` closed forms),
    `.change_of_basis(P, newA)`; `Algebra.hochschild_cohomology(top, coefficients=M)` /
    `hochschild_homology(...)` (a non-`None` coefficient forces the exact **bar** route on
    any Domain — verified: `_use_fast_engine` returns `False`).
  - P31 `families.TrivialExtension(B)` → presented `kQ_T/I_T` Algebra, `dim = 2·dim B`,
    `.quiver` present (verified: `T(kA₂)` has arrows `{a1, te0}`, dual arrow `te0` per socle
    generator).
  - P01 `combinat.Quiver(vertices, arrows).algebra(relations=…, field=…)`; `Algebra.relations`
    (`list[Relation]`, `Relation.terms = ((coeff, (arrow, …)), …)`, `Relation.is_monomial`);
    `Algebra.quiver` (`.arrows` dict, `.vertices`), `Algebra.basis_labels`.
  - P35 `hochschild/products.py` — cited as the Thm-4.1 cup arbiter; **not called** in `src/`
    (coefficient-cup is out of scope, DD-A2), only named on the verification page.
  Branch `plan-72-split-extension-arrow-removal` off `dev`.
- **Merge context (metaplan §3).** `dev` tip at authoring is `d617108` (P62 merged). Before
  this plan's implementation, several branches (P54 BV, P63 wall-chamber, P65 exceptional,
  P67 silting, P68 skew-gentle) merge to `dev` sequentially — **do not assume their internals**;
  consume only merged public surfaces. Re-freeze goldens / recount only under the byte gates.
- **Honest semi-decision contract (metaplan §1.3).** Every report is **complete** iff its
  scope gate holds and the budget did not trip; otherwise `status ∈ {"unsupported","budget",
  "error"}` with a note — **never a guessed dim/certificate**.
- **Char scope is load-bearing (inherited P31/P52).** Batteries run over **QQ**; the
  `GF(32003)` parity check only on the distinct-dim `kA_n` family. Over `char ≤ dim` the
  shipped socle/dual/presented builds raise loudly.
- **No floats in `src/`.** dims/ranks `int`; verdicts `bool`/`str`; connecting maps `Domain`
  matrices; sequences `list[int]`. Client-side float conversion only (`gui.js`, exempt).
- **Composition is left-to-right** (`a*b` = first `a` then `b`); all refusals `QuiverlabError`.
- **Plan-32 markers** (orthogonal; `oracle_*` FORBIDDEN in `tests/{webapp,gui,hpc}` per
  `test_oracle_classes.py`):
  - `oracle_literature`: `HH^•(T(kA_n))` / `HH^•(T(kD₄))` frozen dims (`n = 2,3`; `kD₄`);
    `HH^1(T(B)) = k` on directed `B`; the CLMS `HH_{≥2}` homology iso on the P1–P4
    removed-arrow pairs; the P2 `n = 2` `Ext` correction `= 1`.
  - `oracle_crossengine`: **assembled `HH^•(L)` == direct `HH^•(L)`** (the LES route vs. the
    standalone bar engine, `T(kA₂)`/`T(kA₃)`/`T(kD₄)` — the `T(kD₄)` case is the one that
    EXERCISES the `top+1` M-flank boundary rule); **`arrow_removal` `HH_n(A) == HH_n(B)` for
    `n ≥ 2` == direct** on P1–P4; the leading-piece **`flank_B[0] == HH^0(B)` identity**
    (computed independently over `B`) and the **summand inequality `flank_B[n] ≥ HH^n(B)`**,
    tight/non-vacuous on `T(2`-Kronecker`)` (`flank_B[1] = 3 = HH^1(B)`) — NOT an equality pin
    (would be tautological; the full splitting fails since `D(B)` is bimodule-projective only
    iff `B` self-injective).
  - `oracle_selfcert`: the **LES exactness** (rank alternating-sum at each spot, with the
    `M`-flank to `top+1` so `δ^top`'s codomain is real); the SES of cochain complexes is
    degreewise exact (`dim C^n(L,M) + dim C^n(L,B) = dim C^n(L,L)`); the snake `δ` well-defined
    (independent of the lift chosen — an adversarial-lift test); `Bimodule.inflate(...).check()`;
    the `HH^1` grading cocycle is a cocycle with nonzero class; `dim T(B) = 2·dim B`;
    `remove_arrows` `dim` drop; `inert_arrows` witnesses (incl. the P4 binomial relation);
    `HH_0` unchanged under inert removal on P1–P4 (the provable-invariance check).
  - `qpa`: QPA 1.37 has **no** split-extension-LES or arrow-removal HH surface (live
    `NamesGVars()` guard that FAILS if that ever changes — the P35 precedent); the **inputs**
    are crosscheckable — `TrivialExtensionOfQuiverAlgebra` (dim + `IsSymmetric`, the P31
    oracle, reused) and `HochschildCohomologyHomology` of `A`/`B` on the removed-arrow pairs
    (an input-level anchor for the `HH_{≥2}` iso, not a verdict compare of the LES itself).
- **Mid-merge-train counts.** Absolute suite counts drift; the verification task **recounts
  the oracle-class table at merge time** (`tests/release/test_oracle_classes.py`; paste live
  numbers, claim only this plan's deltas).
- Every merge updates `docs/verification.md` (new oracle rows + recounted class table green)
  and adds `clms_arrow_removal` to `references.bib` + `registry.py` (`bibtex()` hard-fails if
  the two disagree). Conventional commits; green at every commit.

**No-red-commits discipline.** Each Task ships a self-contained green slice: the inflate
primitive + LES machinery before the report before the GUI; recognizers before the reduction.
Fields default to `None`/empty so intermediate task boundaries are green (no forward
reference, no `xfail`).

---

# Task group I — part (a): the split-extension LES (R5)

### Task A1: `Bimodule.inflate` + `split_extension` constructor + the flanks

**Files:**
- Modify: `src/quiverlab/hochschild/coefficients.py` (add `Bimodule.inflate`).
- Create: `src/quiverlab/hochschild/split_extension.py` (`split_extension`, the flank helpers).
- Create: `tests/families/test_split_extension_les_p72.py` (deep bucket).

**Mathematical content.** `Bimodule.inflate(L, B, coeff_over_B, vanishing_arrows)` builds the
`L`-bimodule from a `B`-bimodule by inflation along `π: L ↠ B`: for each generator label of
`L` that is a **vertex `e_v`** or a **`B`-arrow**, copy `coeff_over_B`'s `Lact`/`Ract` for that
label; for each label in `vanishing_arrows` (the `M`-part / dual arrows), use the zero matrix;
then `from_actions` folds to the full `L`-basis and `.check()` certifies. (The general helper
maps EVERY vertex + every `B`-arrow — the multi-vertex `kD₄` case exercises this, not just the
single-arrow `kA₂`.) `split_extension(B, M)` returns `TrivialExtension(B)` when `M is None`
(self-cert `dim = 2·dim B`, `.quiver` present); a general `Bimodule M` is accepted only when a
presented `L` is available (else loud).

- [ ] **Step 1 (RED):** test `Bimodule.inflate(T(kA₂), kA₂, Bimodule.dual(kA₂),
  vanishing_arrows=["te0"]).check()` is `True` and `dim_M = 3`; likewise for
  `Bimodule.regular(kA₂)`. Test `split_extension(kA₂)` is byte-identical to
  `TrivialExtension(kA₂)` (`dim 6`, same quiver arrows `{a1, te0}`). Test the **multi-vertex**
  `T(kD₄)` inflation `.check()` (vanishing arrows = the 3 dual arrows).
- [ ] **Step 2 (GREEN):** implement `Bimodule.inflate` (the label-indexed action copy over ALL
  vertices + `B`-arrows, zero on `vanishing_arrows`, `from_actions` fold, `.check()`);
  `split_extension` (delegate to `TrivialExtension` for `M=None`; the vanishing arrows = the
  dual `te*` arrows = `L.quiver.arrows ∖ B.quiver.arrows`). Refuse loudly on presentation-less
  input.
- [ ] **Step 3 (GREEN):** the flank helpers `_flank_M(L, B, M, top)` = `HH^{0..top+1}(L, M_L)`
  (**to `top+1`, the BOUNDARY RULE**), `_flank_B(L, B, top)` = `HH^{0..top}(L, B_L)`. Assert
  the `T(kA₂)` flank pins (top = 5): `flank_M = [2,1,0,1,2,1,0]` (0..6), `flank_B =
  [1,0,1,2,1,0]` (0..5); and the `T(kD₄)` pins (top = 4): `flank_M = [4,1,0,0,1,4]` (0..5),
  `flank_B = [1,0,0,1,4]` (0..4) — the case with `HH^4(L,B) = 4 ≠ 0`. `.check()` on both
  inflated bimodules (`oracle_selfcert`).
- [ ] **Step 4:** commit `feat(hochschild): Bimodule.inflate + split_extension(B, M) + LES flanks`.

### Task A2: the LES machinery — SES of complexes, snake `δ`, assembly, exactness

**Files:**
- Modify: `src/quiverlab/hochschild/split_extension.py`.
- Modify: `tests/families/test_split_extension_les_p72.py`.

**Mathematical content.** Build the degreewise SES of P52 bar cochain complexes
`0 → C^•(L, M_L) → C^•(L, L) → C^•(L, B_L) → 0` (the coboundary matrices come from
`hochschild.bar.coboundary_matrix(L, n, coefficients=…)`); the maps `ι`/`π` are the coordinate
inclusion/projection induced by `M ⊕ B = L` (the `M`-block / `B`-block of `L`'s basis). Compute
the **snake `δ^n`**: for each cocycle representative `φ ∈ Z^n(L, B_L)`, lift to
`φ̃ ∈ C^n(L,L)`, apply `δ_L`, pull the result back to `C^{n+1}(L, M_L)`, take its class in
`HH^{n+1}(L, M_L)`; `rank δ^n` on classes. Assemble
`dim HH^n(L) = (dim HH^n(L,M) − rank δ^{n-1}) + (dim HH^n(L,B) − rank δ^n)`. Cross-check
against the direct `HH^n(L)`.

- [ ] **Step 1 (RED):** on `T(kA₂)` (top = 5), assert the **assembled** `HH^•(L)` equals the
  **direct** `[3,1,1,1,1,1]`, and the computed snake `delta_ranks` match the dimensionally-
  implied `[0,0,0,2,0,0]`; AND on `T(kD₄)` (top = 4) the assembled `HH^•(L) = [5,1,0,0,2]`
  equals direct with `delta_ranks = [0,0,0,1,2]` — **the case that would FAIL with a `top`-only
  M-flank** (`δ^4` lands in `HH^5(L,M) = 4`) (`oracle_crossengine` + `oracle_selfcert`).
- [ ] **Step 2 (GREEN):** implement the SES-of-complexes builder (assert the degreewise
  `dim C^n(L,M) + dim C^n(L,B) = dim C^n(L,L)` self-cert), the snake `δ` (via `solve` for the
  lift + `reduce_mod_nullspace`/`nullspace` for the class, canonical), the assembly, and the
  exactness gate (the rank alternating-sum holds at every spot).
- [ ] **Step 3 (GREEN):** the **adversarial-lift** self-cert — shifting `φ̃` by any element of
  `C^n(L, M_L)` (the kernel of `π`) leaves `δ^n[φ]` unchanged (well-definedness). And the
  **homology twin** `split_extension_homology` on the dual SES (bar boundary; `HH_•(T(kA₂))`
  assembled == direct `[3,1,1,1,1,1]`).
- [ ] **Step 4:** commit `feat(hochschild): split-extension LES snake connecting map + assembly`.

### Task A3: `SplitExtReport` + the `HH^1 != 0` witness + `Algebra` delegates + scale battery

**Files:**
- Modify: `src/quiverlab/hochschild/split_extension.py`, `src/quiverlab/core/algebra.py`.
- Modify: `tests/families/test_split_extension_les_p72.py`.

**Mathematical content.** `hh1_grading_witness(B, M)`: the grading derivation `E` as an explicit
degree-1 bar cocycle of `L` (`E|_M = id`, `E|_B = 0`), assert `[E] ≠ 0` in `HH^1(L)`; the
directed sharpening `HH^1(L) = 1 + HH^1(B)` when `Q(B)` is acyclic. `SplitExtReport` bundles the
flanks, the `p = 0` leading pieces (`HH^•(B)`, `H^•(B, D(B))` computed over `B`, informational),
the assembly, the direct cross-check, the witness, and the exactness/agreement flags.

- [ ] **Step 1 (RED):** `HH^1(T(B)) != 0` sweep over the zoo — assert `hh1_nonzero` on
  `kA₂/kA₃/kD₄/2-Kronecker/k[x]/(x²)/k[x]/(x³)` and `hh1_directed_formula` True exactly on the
  directed ones (`HH^1(T(B)) = 1 + HH^1(B)`: `kA₂/kA₃/kD₄ → 1`, `2-Kronecker → 4`), N/A
  (None) on the two loop algebras (`oracle_literature`). **Leading-piece pins (NOT equality —
  the tautology fix):** the identity `flank_B[0] == HH^0(B)` and `flank_M[0] == H^0(B,M)`, and
  the inequality `flank_B[n] ≥ HH^n(B)` for all `n`, with the **tight/non-vacuous** case
  `T(2`-Kronecker`)`: `flank_B[1] == HH^1(B) == 3` (`oracle_crossengine`).
- [ ] **Step 2 (GREEN):** implement `hh1_grading_witness`, `SplitExtReport`, and the four
  `Algebra` delegates (lazy import). Directed detection via `A.quiver` acyclicity.
- [ ] **Step 3 (GREEN):** the **scale battery** — assembled == direct on `T(kA₃)` (`[4,1,1,1,1]`)
  and `T(kD₄)` (`[5,1,0,0,2]`), and a `GF(32003)` parity spot-check on `T(kA₂)` (char > dim).
  `oracle_literature` (the frozen `T`-dims) + `oracle_crossengine` (assembled == direct).
- [ ] **Step 4:** commit `feat(algebra): split-extension report + HH^1 grading witness + delegates`.

---

# Task group II — part (b): certified arrow removal / addition (R6)

### Task B1: `inert_arrows` + `remove_arrows` + `add_arrows`

**Files:**
- Create: `src/quiverlab/hochschild/arrow_removal.py`.
- Modify: `src/quiverlab/core/algebra.py` (`inert_arrows` delegate).
- Create: `tests/families/test_arrow_removal_p72.py` (deep bucket).

**Mathematical content.** `inert_arrows(A)` = arrows in no relation word (Def 3.1, read off
`A.relations`, scanning EVERY term of every relation — monomial AND binomial/non-monomial).
`remove_arrows(A, D)` = `Quiver(Q ∖ D).algebra(relations = A.relations verbatim)` after
asserting every `d ∈ D` is inert (else loud with the witnessing relation). `add_arrows(A, F)`
= the tensor-algebra `A_F`, with the relative-cycle finiteness gate (Thm 3.6): a link
`(a₂, a₁)` exists when `s(a₂) A t(a₁) ≠ 0`; a directed cycle of links ⇒ loud infinite-dim refusal.

- [ ] **Step 1 (RED):** `inert_arrows(P1_algebra) == ['c']`, `inert_arrows(P2_algebra) == ['c']`,
  and `x/a/b` NOT inert (they appear in relations); `remove_arrows(A, ['x'])` raises loudly
  (`x` in `x*x`); `remove_arrows(P1, ['c'])` has `dim 3` (`k[x]/(x²) × k`);
  `remove_arrows(P2, ['c'])` has `dim 5`. **The non-monomial `# PIN` (minor h):** on **P4** —
  the commutative square with the **binomial** relation `a*c - b*d` (a `Relation` with two
  terms) plus inert `e:1→4` — `inert_arrows(P4) == ['e']` (`a,b,c,d` all appear across the two
  terms; `e` does not), `remove_arrows(P4, ['e'])` has `dim 9`. `oracle_selfcert` + a
  loud-refusal test.
- [ ] **Step 2 (GREEN):** implement `inert_arrows` (the relation-word scan), `remove_arrows`
  (inert gate + verbatim-relations rebuild), `add_arrows` (tensor algebra + relative-cycle
  gate), and the `Algebra.inert_arrows()` delegate.
- [ ] **Step 3 (GREEN):** the **addition round-trip** self-cert — `add_arrows(remove_arrows(A,
  D), D_as_new) ≅` `A` on the P1/P2 pairs (dim + HH agreement), and a relative-cycle refusal
  (adding a loop to `k[x]/(x²)` ⇒ loud).
- [ ] **Step 4:** commit `feat(hochschild): inert arrows + certified arrow removal/addition`.

### Task B2: `arrow_removal` — the certified HH reduction (Thm 3.2 homology / Thm 4.2 cohomology)

**Files:**
- Modify: `src/quiverlab/hochschild/arrow_removal.py`, `src/quiverlab/core/algebra.py`.
- Modify: `tests/families/test_arrow_removal_p72.py`.

**Mathematical content.** `arrow_removal(A, arrows=None, top, side)`: build `B = remove_arrows(A,
arrows or inert_arrows(A))`; compute `HH_•(A)`, `HH_•(B)` and **assert `HH_n(A) = HH_n(B)` for
`n ≥ 2` ONLY** (Thm 3.2 — the clean homology certificate `hom_agrees`), recording `HH_{0,1}`
agreement informationally in `hom_low_agrees` (`HH_0` provably True); compute `HH^•(A)`,
`HH^•(B)` and split the difference into `coh_correction` (`n ≥ 2`, the Thm 4.2 `Ext` term) and
`coh_low_delta` (`n = 0,1`, the center/disconnection effect) — never a clean-cohomology-iso
claim. Bundle into `ArrowRemovalReport`.

- [ ] **Step 1 (RED):** on **P1**, `hom_A == hom_B == [3,1,1,1,1,1]` (`hom_agrees` `n ≥ 2`;
  `hom_low_agrees == [True, True]`), `coh_A == [1,1,1,1,1,1]`, `coh_B == [3,1,1,1,1,1]`,
  `coh_low_delta[0] == -2` (disconnection), `coh_correction == [0,0,0,0]` (`n ≥ 2`); on **P2**,
  `hom_A == hom_B == [3,0,0,0,0,0]`, `coh_A == [1,1,1,0,0,0]`, `coh_B == [1,0,0,0,0,0]`,
  `coh_low_delta[1] == 1` AND `coh_correction[2] == 1` (the nonzero `n ≥ 2` `Ext` term); on
  **P3** (two loops bridged by inert `c`), `hom_A == hom_B == [4,2,2,2,2,2]`, `coh_A ==
  [1,3,2,2,2,2]`, `coh_B == [4,2,2,2,2,2]`, `coh_low_delta == [-3, 1]`, `coh_correction ==
  [0,0,0,0]` (the low-degree-vs-`Ext` distinction). `oracle_literature` (frozen dims) +
  `oracle_crossengine` (`HH_{≥2}` iso == direct).
- [ ] **Step 2 (GREEN):** implement `arrow_removal` + `ArrowRemovalReport` + the
  `Algebra.arrow_removal(...)` delegate; the inert certificate carries the per-arrow witness;
  `HH_0`-unchanged is asserted on every pair (the provable-invariance self-cert). Honest
  `status`/`note` on any refusal path (non-inert subset, presentation-less).
- [ ] **Step 3 (informational — the open `HH_{0,1}` boundary probe):** record that P1–P4 all
  agree in homology at EVERY degree (so the `n ≥ 2` contract is not vacuously exercised at the
  boundary); document `HH_0`-invariance as **proven** (no cyclic path through an inert arrow),
  and leave the search for a genuine `HH_{0,1}` homology divergence as a NOTE (decided from
  Thm 3.2, not fabricated). This step ships NO false claim — it states the honest contract.
- [ ] **Step 4 (GREEN):** the **trivial-extension crosscheck battery** — `arrow_removal` on a
  `T(B)`-with-an-added-inert-arrow reduces back to `T(B)` with `HH_{≥2}` preserved (ties R6 to
  the R5 flagship + the Han seam). `GF(32003)` parity on a distinct-dim pair.
- [ ] **Step 5 (P73 seam — cross-consistency, whichever plan lands second):** if P73's
  `arrow_removal_subalgebra` is already merged, assert `remove_arrows(A, F) ==
  arrow_removal_subalgebra(A, F)` on inert `F` (same algebra, byte-level dim + HH); if P72
  lands first, this step is a `# NOTE` deferred to P73's merge. Shared-file friction flagged
  in DD-B2.
- [ ] **Step 6:** commit `feat(algebra): certified arrow-removal HH reduction (CLMS Thm 3.2/4.2)`.

---

# Task group III — GUI / webapp + verification

### Task G1: the `split_extension` + `arrow_removal` compute kinds (all three tiers, i18n ×4)

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch` — two `if kind == …` branches beside
  `tau_tilting`/`ar_quiver`).
- Modify: `docs/gui/runner.py` (byte-identical Pyodide twin) + `webapp/static/gui/gui.js`
  (vendored twin) — checkbox / `S.ids` / push-list / `renderBlock` / `scheduleProbe`.
- Modify: `webapp/server/schema.py` (the grammar parse — the third site; algebra-only,
  budget-carrying like `tau_tilting`).
- Modify: `webapp/server/estimator.py` (sizing note — below).
- Modify: `webapp/server/i18n/{en,es,fr,zh}.json` (exact key parity ×4).
- Modify: `webapp/server/runner.py` if it maintains its own kind list.
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE new golden per kind).
- Create: `tests/hpc/test_split_arrow_kinds.py`, `tests/gui/test_split_arrow_runner_twin.py`
  (cross-runner shape parity; unmarked per the Plan-32 extras-gated ruling).

**Mathematical content.** `split_extension` interprets the drawn algebra as `B`, runs
`split_extension_cohomology` (default `top`), and returns the structured block (flanks, leading
pieces, assembled vs. direct, `HH^1` witness). `arrow_removal` runs the auto-inert reduction.
Both are algebra-only, budget-carrying kinds (`split_extension:6` caps `top`; `arrow_removal`
caps `top`). **Estimator sizing (DD-G1):** `split_extension` runs the LES over `L = T(B)`
(`dim = 2·dim B`) with coefficient bar, so `sizing_dim` keys off `2·algebra_dim` (it routes off
"instant" like an oversized family); `arrow_removal` runs HH of `A` and `B` (both `≤ dim A`),
size on `algebra_dim`. Document both in `estimator.py`.

- [ ] **Step 1 (RED):** `tests/hpc/test_split_arrow_kinds.py` — the two kinds through
  `hpc/spec.run` produce the `SplitExtReport`/`ArrowRemovalReport` blocks with the pinned dims
  on `kA₂`/P1; a relation-violating / presentation-less input returns a clean typed 4xx-shaped
  error, never a 500.
- [ ] **Step 2 (GREEN):** wire `_dispatch` (both kinds), the schema grammar (budget parse), the
  runner twin (byte-identical — the cross-runner contract test asserts key-for-key equality),
  the GUI touchpoints, `results_html.py` rendering (the structured-decomposition table for
  `split_extension`; the reduction certificate for `arrow_removal`), and the i18n ×4 keys.
- [ ] **Step 3 (GREEN):** freeze ONE runner golden per kind in `_runner_goldens.json` (dump
  `webapp.server.runner.run_spec` once; dated docstring entry per `test_runner_delegation.py`);
  `tests/webapp/test_js_parses.py` + `test_i18n_parity` green (key parity ×4).
- [ ] **Step 4:** commit `feat(hpc,webapp,gui): split_extension + arrow_removal compute kinds`.

### Task G2: verification page, citations, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/registry.py` (`_r("clms_arrow_removal",
  "cibilslanzilottamarcossolotar2020", "foundation", …)`) + `src/quiverlab/citations/references.bib`
  (the matching entry — `bibtex()` hard-fails on disagreement). Reuse `cmrs_split`,
  `crs_trivial_ext_hh1`, `han_conjecture` (already registered).
- Modify: `docs/verification.md` (new oracle rows: the split-extension LES-vs-direct oracle,
  the `HH^1 ≠ 0` sweep, the arrow-removal `HH_{≥2}` iso + cohomology correction; the honest-scope
  entries: the coefficient-cup deferral DD-A2, the `⊗_B` tensor-power deferral DD-A3, the
  cohomology-correction non-iso). Recount the Plan-32 class table (live numbers).
- Modify: `README.md` (one feature line under the R5/R6 heading).
- Modify: `docs/plans/DEEPER-ENGINES-BACKLOG.md` (record DD-A2 / DD-A3 deferrals for P80).

- [ ] **Step 1 (RED):** `tests/citations/test_bib_structure.py` green with the new key;
  `tests/release/test_oracle_classes.py` recount matches the verification-page table.
- [ ] **Step 2 (GREEN):** add the citation (verify `Proc. AMS 148 (2020) 2421–2432` against
  the AMS landing page at add-time), the verification rows, the README line, the backlog
  deferral entries.
- [ ] **Step 3 (GREEN):** full local gate — `-m fast` + `-m deep` on BOTH kernel paths
  (`QUIVERLAB_NO_NUMBA=1` and numba), `-m qpa` (the `NamesGVars` guard), `test_js_parses`,
  `test_bib_structure`, `test_runner_delegation`; `mkdocs build --strict` if docs changed.
- [ ] **Step 4:** commit `docs(verification,citations): P72 split-extension + arrow-removal oracles`.

---

## Acceptance (Plan-72 definition of done)

1. **Part (a) — the split-extension LES.** `Algebra.split_extension_cohomology(top)` /
   `_homology(top)` return a `SplitExtReport` whose **assembled `HH^•(L)` equals the direct
   `HH^•(L)`** on `T(kA₂)` (`[3,1,1,1,1,1]`), `T(kA₃)` (`[4,1,1,1,1]`), `T(kD₄)`
   (`[5,1,0,0,2]`) — with the flanks computed to the **boundary rule** (`M`-flank to `top+1`,
   the `T(kD₄)` case exercising `HH^5(L,M) = 4 ≠ 0`), the **exactness self-cert** (rank
   alternating-sum) green, the **snake `δ`** well-defined (adversarial-lift test), the `p = 0`
   **leading piece as the `flank_B[0] == HH^0(B)` identity + the `flank_B[n] ≥ HH^n(B)`
   inequality** (tight/non-vacuous on `T(2`-Kronecker`)`, NOT an equality tautology), and the
   **`HH^1(L) ≠ 0` grading-derivation witness** on the whole zoo (with the directed
   `= k ⊕ HH^1(B)` sharpening on `kA_n`/`kD₄`/`m`-Kronecker). CMRS Thm 4.1 (cup) / Cor 3.2
   (decomposition) / Thm 5.5 (`HH^1≠0`) / Def 5.10 (grading derivation) transcribed + cited;
   the coefficient-cup validation and the `⊗_B` tensor-power decomposition are documented
   honest-scope deferrals (DD-A2 / DD-A3).
2. **Part (b) — certified arrow removal/addition.** `Algebra.arrow_removal()` returns an
   `ArrowRemovalReport` with the inert certificate (Def 3.1 witnesses, incl. the P4 binomial
   relation), the clean **`HH_n(A) = HH_n(B)` for `n ≥ 2` ONLY** (Thm 3.2) on the P1–P4 pairs
   `== direct` (with `HH_0` provably invariant, `HH_{0,1}` reported), and the honest cohomology
   split — **`coh_correction` (`n ≥ 2`, Thm 4.2 `Ext` term, nonzero `HH^2` on P2)** separated
   from **`coh_low_delta` (`n = 0,1`, center/disconnection)**. Loud refusals: non-inert arrow
   (with witness), presentation-less, relative-cycle addition. The `inert_arrows` /
   `remove_arrows` / `add_arrows` primitives are exposed for **P73** (the Han seam — no Han
   logic in P72; shared-substrate cross-consistency with `arrow_removal_subalgebra` per DD-B2).
3. **GUI / three tiers.** `split_extension` + `arrow_removal` algebra-only budget kinds served
   by `hpc/spec.py` + the byte-identical `docs/gui/runner.py` twin + `webapp/server/schema.py`,
   with the recognizer-panel rows, report rendering, i18n ×4 exact parity, canonical-key
   stability, one frozen runner golden per kind, and the estimator sizing note (DD-G1).
4. **Oracles + verification.** All Plan-32-marked (`oracle_literature` / `oracle_crossengine` /
   `oracle_selfcert` / `qpa`); `docs/verification.md` carries the new rows + the honest-scope
   entries and the recounted class table; `clms_arrow_removal` added to `references.bib` +
   `registry.py` (BibTeX-verified `Proc. AMS 148 (2020) 2421–2432`); `cmrs_split` /
   `crs_trivial_ext_hh1` / `han_conjecture` reused.
5. **Green everywhere.** `-m fast` + `-m deep` on both kernel paths, `-m qpa`, `test_js_parses`,
   `test_bib_structure`, `test_runner_delegation`, `test_oracle_classes` all green;
   conventional commits, green at every commit; branch `plan-72-split-extension-arrow-removal`
   merged to `dev` with a recount.

---

## Design decisions (recorded for the critic + P80)

- **DD-A1 (the flagship is the trivial extension).** `split_extension` centers on
  `M = D(B)` (`L = TrivialExtension(B)`, presented, Plan 31) because the coefficient inflation
  needs `L.quiver`. General `B ⋉ M` (arbitrary `Bimodule M`) is supported only with a presented
  `L`; else loud. This matches the card ("trivial-extension inputs") and R5's oracles
  (`T(kA_n)`/`T(kD₄)`).
- **DD-A2 (the connecting map is the snake, cited as the cup).** `δ` is COMPUTED as the
  snake-lemma connecting homomorphism of the P52 bar SES (linear algebra, no new engine); CMRS
  Thm 4.1 identifies it with the cup-with-`1_M` and is cited as the theoretical identity. A
  literal coefficient-cup recomputation via Plan-35 is out of scope (Plan-35 cup is
  coefficient-`A`); `δ` is validated by exactness + assembled-≡-direct. Honest-scope on the
  verification page.
- **DD-A3 (the `⊗_B` tensor-power decomposition — Cor 3.2 — is deferred, and does NOT fully
  split here).** The card's "from `HH^*(B) + H^*(B,M)`" is realized as the `p = 0` **leading
  piece** surfaced from the `L`-flanks (which already contain all `p`), pinned as the `n = 0`
  identity + the summand inequality — **NOT a degreewise equality**, because Cor 3.2's full
  splitting requires `M` one-sided `B`-projective, which `D(B)` is **iff `B` is self-injective**
  (false for `kA_n`/`kD₄`). A standalone `M^{⊗_B p}` + bimodule-`Ext` engine (the genuine
  acceleration) is a ledger deferral for a later plan. The LES-over-`L` route is card-complete
  without it.
- **DD-B1 (homology clean `n ≥ 2`, cohomology corrected + low-degree separated).** Arrow
  removal reports (i) the clean `HH_{n≥2}` homology isomorphism (Thm 3.2, the only asserted
  range; `HH_0` provably invariant, `HH_{0,1}` reported), and (ii) the cohomology difference
  split into `coh_correction` (`n ≥ 2`, Thm 4.2 `Ext` term — nonzero `HH^2` on P2) and
  `coh_low_delta` (`n = 0,1`, center/disconnection — e.g. P1 `−2`, P3 `[−3, +1]`), never a
  false clean-cohomology-iso claim and never conflating the two phenomena.
- **DD-B2 (the P73 seam — shared substrate).** P72 owns the reduction primitives
  (`inert_arrows`, `remove_arrows`, `add_arrows`, `arrow_removal`); P73 (R7, Han +
  Jacobi–Zariski) owns the recognizer/transport built on them — no Han decision, no JZ sequence
  in P72. P73 specs `arrow_removal_subalgebra(A, F)` with **more general semantics** (arbitrary
  `F`, induced ideal `I_B = ker π`), of which P72's `remove_arrows` is the **inert special case**
  (`I_B = I` verbatim). **Whichever plan's module lands on `dev` first OWNS the constructor;
  the second CONSUMES it** and adds the cross-consistency test (`remove_arrows ≡
  arrow_removal_subalgebra` on inert `F`, Task B2 Step 5). The shared-file merge friction is
  flagged in both plan docs.
```
