# Plan 59: Recognizer Batteries — the homological string test + toupie algebras (R34 + R35)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking. Do the reference re-verification in
> the Record section BEFORE writing any oracle pin — two formulas below are `# PIN`/
> `# BLOCK` and must be resolved against the PDF, not from this plan's prose.

**Goal.** Two independent recognizer/oracle batteries, each grounding a *research
characterization* on the machinery quiverlab already ships:

- **(R34) The homological string-algebra test.** Suárez-Álvarez (2023): among
  **representation-finite** algebras, `A` is a **string** algebra **iff the middle term
  of every extension of indecomposable modules has at most two indecomposable
  summands** — the quantifier is over **all** `Ext¹(X,Y)` classes, not merely
  Auslander–Reiten sequences. We deliver `homological_string_test(A)` — an **honest
  three-valued semi-decision** (definitive `True`/`False` where the quantifier is
  finitely checkable, loud `inconclusive` otherwise) built on `knit_ar_quiver` (the
  finite indecomposable universe), `baer_extension` (realizing the middle term of a
  cocycle), and `decompose`/Fitting-splits (counting/witnessing summands). The
  headline oracle is the **discriminating battery vs `is_string`** (the P38/P46
  syntactic recognizer): wherever the homological test is definitive it must equal
  `is_string`, and a definitive disagreement is a **loud error** — one of the two is a
  bug.

- **(R35) Toupie algebras.** Artenstein–Lanzilotta–Solotar (2020): a **toupie** is a
  `kQ/I` whose quiver has a **unique source, a unique sink, and every other vertex has
  in-degree 1 and out-degree 1** (a vertex-disjoint bundle of `a` directed *branches*
  from source to sink), `I` an admissible ideal (monomial branch-truncations +
  non-monomial cross-branch linear relations). We deliver a **`ToupieAlgebra`
  constructor** + an **`is_toupie` graph-shape recognizer**, a **closed-form HH oracle
  family** (the `a`-Kronecker `HH^• = [1, a²−1, 0, …]` — reusing the existing P29 pin —
  and cross-engine bar≡CS on subdivided toupies), and the **`sl_a ⊆ HH¹` inclusion**
  identified over **char 0 only** (`a` = number of branches = out-deg(source) =
  in-deg(sink); `dim HH¹ ≥ a²−1`, with equality on the `a`-Kronecker).

**GUI.** Two new **algebra-only** compute kinds (the expensive-but-scalar `ar_quiver`
precedent, NOT the fast `recognizers` block): `string_homological` and `toupie`, each a
recognizer-panel row, both runners byte-identical, **all four locales (en/es/fr/zh —
`webapp/server/i18n/*.json`, exact key parity enforced by `tests/webapp/test_i18n.py`)**,
canonical keys stable, with report rendering + citations. A `ToupieAlgebra` preset in the
families input surface.

**Architecture.** Two new source modules + two new algebra-only kinds; everything else
is a thin exact layer over primitives that already exist on `dev` (no new math engine):

- `src/quiverlab/modules/string_homological.py` — `homological_string_test(A, …) ->
  StringHomologyVerdict` and `string_homological_block(A)`. Consumes
  `modules/ar.py::knit_ar_quiver` (finite indecomposables), `modules/ar.py::_ext1_data`
  + `_combine_matrices` (cocycle for a class over the `Ext¹` basis — the
  `almost_split_sequence` idiom), `modules/yoneda.py::baer_extension` (the middle
  term), `modules/ses.py::ShortExactSequence` (exactness self-cert),
  `modules/decompose.py::decompose`/`is_indecomposable` (summand count / ≥3 witness),
  `invariants/recognizers.py::is_string` (the arbiter), `modules/ext.py::ext_dims`
  (`with_reps` route for the basis classes).
- `src/quiverlab/families/toupie.py` — `ToupieAlgebra(branches, relations=None,
  field=None) -> Algebra`, `is_toupie(A) -> bool`, `toupie_branch_count(A) -> int`,
  `toupie_block(A)`. Consumes `combinat/quiver.py::Quiver` + `Quiver.algebra` (the
  `QuantumCI` non-monomial template), `Algebra.hochschild_cohomology` (the closed
  forms + `sl_a`), and mirrors the `families/nakayama.py` loud-hint validation /
  `families/trivial_extension.py` dimension-certificate idiom.
- The two algebra-only kinds are wired into `hpc/spec.py::_dispatch` +
  `docs/gui/runner.py` (byte-identical twins), the GUI touchpoints (checkbox / `S.ids`
  / push-list / `renderBlock` / `scheduleProbe`), the i18n chains, and
  `trace/results_html.py`, exactly following the `ar_quiver` kind (also an expensive
  algebra-only kind seeded from the AR machinery).

**Tech Stack.** Pure exact linear algebra over `Domain` (`modules/linalg_mod`,
`fields/linalg`) + exact combinatorics on the quiver; HH via the shipped engines
(`engine`/`resolutions_cs`, dispatched through `Algebra.hochschild_cohomology`). **No
floats in `src/`** (AST-gated by `tests/test_no_floats.py`): toupie relations are exact
grammar tokens (the `QuantumCI` `_q_token` precedent), all matrices `0`/`1`/int-mod-`p`/
`Fraction`. The `sl_a` inclusion is char-0 (`QQ`); the HH-dimension pins are
char-independent.

---

## Records (verbatim)

> **R34 — Homological rep-finite string-algebra test.** [A-scout P6, CORRECTED
> statement] Object: among rep-finite algebras, A is string ⇔ the middle term of every
> extension of indecomposables (arbitrary Ext¹ classes — NOT just AR sequences) has ≤ 2
> indecomposable summands; finite check over all Ext¹(X,Y)-middle-terms on the knitted
> category. Ref: Suárez-Álvarez arXiv:2105.02948 (Alg. Rep. Theory 26 (2023)
> 1759–1772); the Huisgen-Zimmermann–Smalø citation dropped (wrong content — critic).
> Oracle: discriminating battery vs is_string across the zoo. Size M.

> **R35 — Toupie algebras: recognizer + closed-form HH + sl_a structure.** [A-scout P5;
> keep-with-corrections] Object: trivial graph recognizer (unique source/sink, parallel
> paths); closed HH forms as an oracle family; HH¹ ⊇ sl_a identified over ℂ/char 0 ONLY
> (scope per critic). Refs: Artenstein–Lanzilotta–Solotar arXiv:1803.10310 (Alg. Rep.
> Theory 2020 — cite the journal version); Artenstein thesis (Colibri UdelaR). Size S–M.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P59 (Wave 2, tier β,
independent).

---

## Reference re-verification (done at authoring; findings binding)

**R34 — Suárez-Álvarez, arXiv:2105.02948** (abstract + ar5iv full text fetched):
- **Theorem A (verbatim intent):** *"A finite dimensional algebra of finite
  representation type is a string algebra if and only if the middle term of every
  extension between its indecomposable modules has at most two direct factors."*
- **Quantifier:** over **all** extensions of indecomposables, **not** AR sequences —
  the introduction is explicit: *"It is well-known that the middle term of every almost
  split short exact sequence … has at most two direct factors. In this note I will show
  that, in fact, this is true of all extensions of indecomposable modules."* This is
  the load-bearing correction: the AR sequences alone do **not** decide the property.
- **Corollary 4 (`|E| ≤ |M| + |N|`):** this is the **string-algebra direction**
  (*only if*), NOT an unconditional fact. Sanity check that pins the sign of the whole
  battery: `kD₄` (path algebra, subspace orientation, three arrows into the centre) is
  rep-finite and **not** special-biserial; its central AR mesh `0 → τM → E → M → 0` has
  a **3-summand** middle with `M`, `τM` indecomposable, so `|E| = 3 > 2 = |M| + |N|` —
  which is exactly why `kD₄` is **not** a string algebra. If Corollary 4 held
  unconditionally the theorem would be vacuous; it is not. **`# PIN`** the exact
  statement/number of Corollary 4 and the definition of `|·|` (= number of
  indecomposable direct summands) against the PDF before quoting it in a docstring.
- **Field:** *"Throughout this paper 𝕜 denotes an algebraically closed field."* No
  general characteristic restriction. **Scope consequence (honest):** our exact test
  runs over the **ground field** (`GF(p)`/`GF(pⁿ)`/`QQ`), not `k̄`; the theorem's *iff*
  is stated over `k̄`. The specific battery examples are chosen field-robust (see Task 2)
  and the `k̄` gap is documented on the verification page.
- **Decidability:** the paper gives the *characterization*, **no algorithm**. Design
  question 1 (below) is therefore ours to settle honestly.

**R35 — Artenstein–Lanzilotta–Solotar, arXiv:1803.10310 = ARTh 23 (2020) 421–456,
DOI 10.1007/s10468-019-09854-y** (ar5iv full text + the P29 deep-research entry
`docs/plans/2026-07-25-literature-oracles-deep-research.md` §`[alsolotar_toupie]`):
- **Definition (Def. 1):** unique source `0`, unique sink `ω`, *"any other vertex is
  the source of exactly one arrow and the target of exactly one arrow"* → a
  vertex-disjoint bundle of `a` directed **branches** (parallel paths) from `0` to `ω`;
  `I` an **admissible** ideal of **monomial** (branch-internal) + **non-monomial**
  (cross-branch linear combinations of the full `0→ω` branch-paths) relations. Same
  ordinary quiver as the canonical algebras.
- **`a`-Kronecker `Q_a`** (2 vertices `0 ⇒ ω`, `a` parallel arrows, **no relations** —
  the length-1 toupie): **`HH^• = [1, a²−1, 0, …]`** (hereditary ⇒ `HHⁱ = 0`, `i ≥ 2`;
  `dim HH¹ = (a−1) + a(a−1) = a²−1`; **char-independent**, Euler-verified `χ = 2 − a²`).
  `a=2 → [1,3,0,…]`, `a=3 → [1,8,0,…]`, `a=4 → [1,15,0,…]`. **This pin already ships**
  as the `m`-Kronecker literature oracle at `tests/engine/test_literature_p29.py:214`
  (`assert hh == [1, m*m - 1, 0, 0]`) — P59 reuses it and re-cites it under
  `alsolotar_toupie`.
- **`sl_a` (Thm 6.5):** for **`k = ℂ`** (stated *"Let 𝕜 be a field of characteristic
  zero"*), `HH¹(A)` contains a Lie subalgebra `≅ sl_a(ℂ)`, where **`a` = the number of
  arrows/branches from `0` to `ω`**. On the `a`-Kronecker, `dim HH¹ = a²−1 = dim sl_a`
  exactly (the whole `HH¹` is `sl_a`). **`# PIN`** the exact wording of "the number of
  arrows from `0` to `ω`" — we implement `a = out_degree(source) = in_degree(sink)`
  (asserting the equality) and pin `dim sl_a = a²−1`; confirm against Thm 6.5 that `a`
  counts **branches** (not e.g. only length-1 direct arrows) before the test asserts.
- **`# BLOCK` — Example 7.4.1** (13 vertices / 15 arrows, `HH^• = [1,10,3,0,4,0,…]`):
  the quiver is given **only as a figure**; the P29 deep-research reconstruction did
  **not** reconcile with the stated vertex/arrow counts. **Do NOT pin Example 7.4.1
  without consulting the PDF figure.** The `a`-Kronecker is the safe literature pin; all
  other toupie HH values are **cross-engine** (bar ≡ CS on constructed toupies), which
  is the honest oracle for the general closed forms (ALS defers the *dimension* formulas
  to prior work `[14]` and gives *bases*/Gerstenhaber structure — so there is no single
  transcribable closed-form-dimension theorem to pin beyond the `a`-Kronecker).

**Citations to add** (both BibTeX-verified; `alsolotar_toupie` DOI confirmed above):
`suarez_alvarez` (2105.02948 / ARTh 26 (2023) 1759–1772) and `alsolotar_toupie`
(`ArtensteinLanzilottaSolotar2020`). `butler_ringel` (Butler–Ringel 1987) already ships
from P46. The Artenstein thesis (Colibri UdelaR) is a `note=` on the same entry (not a
separate BibTeX key — thesis metadata not journal-verifiable).

---

## Design question 1 — the quantifier, and what is decidable (R34)

The record says *"finite check over all Ext¹(X,Y)-middle-terms on the knitted
category."* The subtlety the reviewer flagged is real and must be handled honestly:

1. **The indecomposable universe is finite.** `A` rep-finite ⇒ `knit_ar_quiver(A)`
   returns `status="complete"` with finitely many indecomposables `{X_1,…,X_n}` (their
   `Module`s). We only run when the knit is complete; otherwise **loud refusal** (out of
   scope — the theorem is stated for rep-finite `A`). Self-injective input is refused by
   the knitter already (`status="unsupported"`) → refuse loudly there too.

2. **The middle term varies across `Ext¹(X,Y)`, and it is NOT constant.** The zero class
   gives `X ⊕ Y` (exactly 2 summands — always fine). Nonzero classes give middle terms
   that depend on the class (up to scalar and up to the `End(X) × End(Y)` action), and
   over an **infinite** field there are infinitely many classes. The number of summands
   is **not** monotone/semicontinuous in a way that lets a generic sample certify the
   `≤ 2` bound — a *special* (closed-locus) class can carry **more** summands than the
   generic one, so no finite sample proves the universal `≤ 2`.

3. **Two directions, asymmetric difficulty (ORCHESTRATOR RULING, 2026-08-07 review):**
   - **Refute (definitive `not_string`, `k̄`-sound)** — *find one class with a `≥ 3`
     middle.* "String" is a property of the **presentation** `(Q, I)`, field-independent;
     a `≥ 3`-summand middle over **any** ground field certifies non-string (a direct-sum
     `E = A⊕B⊕C` over `k` base-changes to `E⊗k̄ = A⊗k̄ ⊕ B⊗k̄ ⊕ C⊗k̄`, still `≥ 3`
     summands — so a ground-field witness **base-changes up** to `k̄`). Only the
     **terminal split step** is char-robust (a found idempotent is exact); the **refute
     PIPELINE** feeding it (`knit_ar_quiver` / `almost_split_sequence` / `decompose`)
     needs `char 0` or `char > dim` — so refute runs over **QQ** in practice (and over
     small primes only on tiny algebras where `char > dim` holds throughout the
     pipeline). Checked classes: the split (trivially 2), **every AR-sequence** middle
     (`almost_split_sequence`), and a **basis** of each `Ext¹(X,Y)` (via
     `ext_dims(..., with_reps=True)` → the `interpretation` middle modules, or the
     `_ext1_data` + `_combine_matrices` idiom on unit coordinate vectors). Any `≥ 3` ⇒
     definitive `not_string` **with a witness**.
   - **Confirm (NON-definitive w.r.t. the `k̄` theorem) — `verdict="string_over_ground_field"`.**
     *Every class over the ground field has `≤ 2`.* Only finitely checkable over a
     **finite ground field**: enumerate `Ext¹(X,Y)` up to the `End(X)×End(Y)`-orbit +
     scalar (`end_action_on_ext1`, budget-capped), build each middle via `baer_extension`,
     certify `≤ 2` via `decompose` (needs `char > dim E` OR `dim End = 1` OR a found ≤2
     split — else that pair is char-undecidable). This is **INCONCLUSIVE** with respect
     to the Suárez-Álvarez theorem: over `k̄` there may be **more** indecomposables and
     **non-rational** `Ext¹` classes that a finite ground field cannot see, so "all `≤ 2`
     over `GF(q)`" does **not** prove "all `≤ 2` over `k̄`". The **definitive string
     verdict is carried by the syntactic `is_string`** (below), not by this pass.

4. **The verdict enum is `{"not_string", "string_over_ground_field", "inconclusive"}`**
   (the house semi-decision contract, cf. P46 `enumerate_strings`). There is **no** bare
   `"string"` verdict — a definitive `k̄` "string" answer is not producible by the
   homological test alone (it is the syntactic `is_string`'s job).
   - `verdict="not_string"` — a `≥ 3` middle witnessed; `k̄`-sound and field-independent.
   - `verdict="string_over_ground_field"` — finite-ground-field exhaustive confirm closed
     (every middle `≤ 2` over `GF(q)`); inconclusive w.r.t. the `k̄` theorem; defers to
     `is_string`.
   - `verdict="inconclusive"` — infinite ground field with no `≥ 3` witness in the checked
     sample, OR a char-undecidable middle, OR budget/knit tripped. The payload names the
     reason and the checked-class set (`split`, `ar`, `basis`, `exhaustive`). **Never** a
     silent string claim.

5. **The contract is ONE thing: refute-side disagreement RAISES; everything else RETURNS.**
   `is_string` (P38/P46) is a decidable syntactic recognizer (special-biserial + monomial
   ideal), always definitive — it carries the definitive `string` answer.
   - The **only** loud disagreement: if a `≥ 3` witness is found (⇒ genuinely non-string,
     `k̄`-sound) **and** `is_string(A) is True`, raise `QuiverlabError` — one engine is a
     bug. This is the discriminating oracle's teeth. In all other cases the function
     **returns** the structured `StringHomologyVerdict` (no raise).
   - A `string_over_ground_field` verdict on an algebra with `is_string(A) is False` is a
     **`k̄`-gap**, recorded in `kbar_gap_note` and surfaced in the block/report — it is
     **never raised** (the finite ground field simply did not witness the non-stringness
     the `k̄` theorem guarantees). `inconclusive` never asserts agreement either.

6. **Deliberately NOT implemented (documented non-goal).** The alternative
   "middle-term-candidate" route — enumerate candidate modules `M'` from the knit with
   `dim M' = dim X + dim Y` and decide *"∃ SES `0→Y→M'→X→0`"* — is **not** a shortcut to
   a field-general *confirm*: realizability of a fixed `M'` as an extension of `X` by `Y`
   is itself an existence question over an infinite space (an epi `M' ↠ X` with kernel
   `≅ Y`), decidable in principle only with algebraic-geometry machinery quiverlab does
   not have. v1 uses the cocycle route above; the candidate route is stated as a scoped
   non-goal on the verification page, not silently skipped.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **Prerequisites present on `dev` (verify before branching).** P41 (AR completion —
  `modules/ar.py`: `knit_ar_quiver`, `almost_split_sequence`, `end_action_on_ext1`,
  `_ext1_data`, `_combine_matrices`, `ARQuiver`), P37 (`modules/ses.py`,
  `modules/yoneda.py::baer_extension`), P38 (`invariants/recognizers.py::is_string` /
  `is_special_biserial`), P30 (`modules/decompose.py::decompose`/`is_indecomposable`),
  P23/24 (`modules/duality.py::tau`, `modules/hom.py::is_isomorphic`/`identify_standard`),
  P05/40 (`modules/ext.py::ext_dims`, `modules/complex_reps.py::ext_reps`). Verify:
  `python -c "import quiverlab.modules.ar, quiverlab.modules.yoneda,
  quiverlab.modules.decompose, quiverlab.invariants.recognizers"`. **No dependency on
  the other Wave-2 plans** (P55–P58, P60–P62) — P59 merges independently. If P46's
  `strings/` package is on `dev` it is convenient (a source of string zoo examples) but
  **not required**: the R34 arbiter is `is_string` (P38), which does not need P46.
- **Test buckets are auto-assigned by directory** (`tests/conftest.py`):
  `tests/modules/` and `tests/families/` → **deep**; `tests/webapp/`, `tests/gui/`,
  `tests/hpc/`, `tests/invariants/` → **fast**; `tests/qpa/` → **qpa**. R34 batteries
  live in `tests/modules/test_string_homological_p59.py` (deep — they knit, realize
  extensions, decompose); R35 in `tests/families/test_toupie_p59.py` (deep — it
  constructs algebras and computes HH). GUI/runner cross-tests in `tests/webapp/` +
  `tests/gui/` (fast); QPA in `tests/qpa/` (qpa). Run by path during development; finish
  each task with a `-m deep` / `-m fast` / `-m qpa` spot-run of touched files.
- **The `decompose` char-caveat is load-bearing (R34).** `decompose`/`is_indecomposable`
  refuse when `char ≤ dim M` and neither `dim End = 1` nor a Fitting split decides
  (`modules/decompose.py::_certify_local`). Consequence for R34: the **refute** side
  (finding a `≥ 3` split) is char-robust; the **confirm** side (certifying `≤ 2`) needs
  `char > dim E` or `char 0`. The batteries therefore run the **witness** examples
  (`kD₄`, …) over **QQ** (char 0 — `decompose` always certifies), and the **exhaustive
  confirm** examples over **small primes on small algebras** (`kA₂`/`kA₃`, middles
  `dim ≤ 4`) where `char > dim` holds or the middle splits into `≤ 2` pieces outright. A
  char-undecidable middle makes that pair `inconclusive`, never a guess.
- **Loud refusals, honest three-valued verdicts** (the `GlobalDimension` /
  `enumerate_strings` precedent): `homological_string_test` refuses loudly on
  rep-infinite / self-injective / presentation-less / non-complete-knit input;
  `ToupieAlgebra` refuses loudly on a non-toupie shape or a failed dimension
  certificate; `is_toupie` refuses loudly on a presentation-less (structure-constant)
  algebra; the `sl_a` block refuses loudly off char 0. Every refusal is `QuiverlabError`
  with a `hint=`.
- **Convention/scope choices are ARBITRATED or GATED, not assumed** (house
  cup-sign/composition-order precedent): the `is_string`↔homological agreement is the
  arbiter for R34 (a disagreement is loud); the toupie `a`-count is arbitrated by
  `out_deg(source) == in_deg(sink)` + the `a`-Kronecker `dim HH¹ = a²−1` pin; the
  `sl_a` inclusion is **hard-gated** to char 0.
- **Plan-32 markers** (orthogonal): `oracle_selfcert` = every realized middle
  self-certifies (`baer_extension` cocycle check + `ShortExactSequence` exactness), the
  AR-mesh `≤ 2`, the toupie dimension certificate + recognizer partition, `Δ`-free
  identities; `oracle_crossengine` = the **discriminating battery** homological ≡
  `is_string`, toupie HH bar ≡ CS degreewise, `is_toupie` ≡ hand quiver shape;
  `oracle_literature` = `kD₄` 3-middle (Butler–Ringel/ARS), `kAₙ` all-`≤2`, the
  `a`-Kronecker `[1,a²−1,0,…]` + `dim sl_a = a²−1`; QPA lives in `tests/qpa/` (bucket =
  the class, never double-marked).
- **Mid-merge-train counts.** v1.0.0 lands many subplans in overlapping waves; absolute
  suite counts drift. **Task 6 recounts the oracle-class table at merge time** by
  running `tests/release/test_oracle_classes.py` (paste the LIVE numbers, never a
  guessed-at-authoring count) and claims only the deltas this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class
  table green) and adds its citations to `citations/references.bib` + `registry.py`.
  Conventional commits; green tests at every commit; branch `plan-59-recognizer-batteries`
  off `dev` (do not commit/push until asked).

---

### Task 1: `homological_string_test` — the extension/middle-term semi-decision (R34)

The core. Build the finite indecomposable universe, walk ordered pairs, realize
middle terms of extensions, and return the three-valued verdict + witness.

**Files:**
- Create: `src/quiverlab/modules/string_homological.py`
- Test: `tests/modules/test_string_homological_p59.py`

**Interfaces:**
- Consumes: `modules/ar.py::knit_ar_quiver(A) -> ARQuiver` (`.vertices` = list of
  `{"name","dimvec","module"}`, `.is_complete`, `.status`),
  `modules/ar.py::almost_split_sequence(M) -> ShortExactSequence`,
  `modules/ar.py::_ext1_data(A, X, Y) -> (cocycle_mats, cob, terms, dmats)` and
  `_combine_matrices(cocycle_mats, coeffs, rows, cols, dom)` (the class→cochain idiom
  from `almost_split_sequence`; `rows = Y.dim`, `cols = width` of `dmats[1]`),
  `modules/ext.py::ext_dims(A, X, Y, 1, with_reps=True, interpret=True)` (basis classes
  + Yoneda `interpretation` middle modules), `modules/yoneda.py::baer_extension(X, Y, f,
  terms, dmats) -> YonedaSequence` (`.modules[1]` is `E`; refuses if `f` is not a
  cocycle), `modules/ses.py::ShortExactSequence`, `modules/decompose.py::decompose` /
  `is_indecomposable`, `invariants/recognizers.py::is_string`,
  `modules/ar.py::end_action_on_ext1(X, Y)` (orbit reduction for the finite-field
  exhaustive route).
- Produces:
  ```python
  @dataclass(frozen=True)
  class StringHomologyVerdict:
      verdict: str            # "not_string" | "string_over_ground_field" | "inconclusive"
      is_string: bool         # the P38 syntactic arbiter (carries the definitive k-bar
                              #   "string" answer)
      witness: dict | None    # on "not_string" -- unified schema across both refute
                              #   sub-routes (see _mk_witness): {"X","Y","class","coeffs",
                              #   "E_dimvec","summand_count","summand_dimvecs"};
                              #   class in {"ar_socle","basis"}; coeffs is None on the AR route
      kbar_gap_note: str | None  # set iff verdict=="string_over_ground_field" AND is_string
                              #   is False: the finite ground field did not witness the
                              #   non-stringness the k-bar theorem guarantees (NOT an error)
      checked: dict           # {"classes": int, "route": [..], "field": str, "n_indec": int}
      reason: str | None      # why inconclusive (field-infinite / char-undecidable /
                              #   budget / knit-incomplete)
      references: tuple       # ("suarez_alvarez", "butler_ringel", "assem_book")

  def homological_string_test(A, budget_classes=4096, budget_pairs=4096)
          -> StringHomologyVerdict:
      # rep-finite ONLY (knit complete). Refuses loudly (QuiverlabError) on out-of-scope
      # input. RAISES iff a >= 3 witness is found AND is_string(A) is True (the one loud
      # disagreement -- k-bar-sound). Otherwise RETURNS the structured verdict.
  def string_homological_block(A) -> dict          # the algebra-only compute kind (Task 5)
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_string_homological_p59.py
"""The homological string-algebra test (Plan 59 / R34, Suarez-Alvarez 2023).
Discriminating oracle: a >= 3-summand middle (verdict "not_string") is k-bar-sound and,
if is_string(A) is True, RAISES (one engine is a bug). The finite-field confirm yields
"string_over_ground_field" (inconclusive w.r.t. the k-bar theorem) and defers to the
syntactic is_string; a mismatch there is a recorded k-bar-gap, never raised. Refute
runs over QQ (the pipeline needs char 0 / char > dim); the finite-field confirm runs
over small primes on small algebras (char > dim of the middles)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import is_string
from quiverlab.modules.string_homological import homological_string_test

xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _kA(n, field=QQ):
    arrows = {chr(ord("a") + i): (i + 1, i + 2) for i in range(n - 1)}
    return Quiver(list(range(1, n + 1)), arrows).algebra(relations=[], field=field)


def _kD4(field=QQ):
    # subspace orientation: three arrows into the centre 0. Rep-finite (Dynkin D4),
    # NOT special-biserial (out/in-degree 3), so NOT a string algebra.
    Q = Quiver([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)})
    return Q.algebra(relations=[], field=field)


def _gentle_a3(field=QQ):
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=field)


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_linear_A_over_QQ_never_refutes_string(n):
    # kA_n: every extension of indecomposables has <= 2 middle summands (string). Over
    # QQ (infinite field) the confirm side cannot exhaust, so the honest verdict is
    # "inconclusive" -- and NEVER "not_string" (no >= 3 witness exists).
    v = homological_string_test(_kA(n))
    assert v.is_string is True
    assert v.verdict == "inconclusive"                 # infinite field: no k-bar confirm
    assert v.witness is None


@lit
def test_kD4_is_not_string_with_the_verified_three_middle_witness():
    # THE discriminating witness (VERIFIED LIVE 2026-08-07, dev, QUIVERLAB_NO_NUMBA=1):
    # knit_ar_quiver(kD4) is complete with 12 indecomposables; the AR sequence ending at
    # the module M with dim vector {0:2,1:1,2:1,3:1} has middle E = P1 (+) P2 (+) P3, the
    # three arm-projectives {0:1,1:1} / {0:1,2:1} / {0:1,3:1} -- EXACTLY 3 summands.
    v = homological_string_test(_kD4())
    assert v.is_string is False                        # kD4 is not special-biserial
    assert v.verdict == "not_string"
    assert v.witness is not None and v.witness["summand_count"] == 3
    assert v.witness["E_dimvec"] == {0: 3, 1: 1, 2: 1, 3: 1}
    # this verdict is k-bar-sound and is_string is False, so the function RETURNS (no
    # raise); a raise would occur only if is_string had (wrongly) been True.


@xeng
@pytest.mark.parametrize("factory", [_kA, _gentle_a3])
def test_returned_verdict_is_consistent_with_is_string(factory):
    A = factory(3) if factory is _kA else factory()
    v = homological_string_test(A)
    # a returned "not_string" is only produced when is_string is False (else it raises);
    # "string_over_ground_field"/"inconclusive" defer to is_string (no assertion of equality).
    if v.verdict == "not_string":
        assert v.is_string is False


@xeng
def test_finite_field_exhaustive_confirm_on_kA3():
    # small prime, small middles (dim E <= 4 < char 5 so decompose certifies): the confirm
    # side closes -> "string_over_ground_field" (NOT a definitive k-bar "string"; that is
    # is_string's job). No k-bar gap here since is_string is True.
    v = homological_string_test(_kA(3, field=GF(5)))
    assert v.verdict == "string_over_ground_field"
    assert v.is_string is True and v.kbar_gap_note is None


@selfcert
def test_witness_middle_self_certifies():
    v = homological_string_test(_kD4())
    # the recorded witness is a genuine >= 3 decomposition of a genuine extension
    assert sum(1 for _ in v.witness["summand_dimvecs"]) == v.witness["summand_count"]
    assert v.witness["class"] in ("ar_socle", "basis")


@selfcert
def test_rep_infinite_and_selfinjective_and_presentationless_refused():
    from quiverlab.families import NakayamaAlgebra
    kron = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    with pytest.raises(QuiverlabError):
        homological_string_test(kron)                  # 2-Kronecker: rep-infinite (knit incomplete)
    sinj = NakayamaAlgebra(n=3, l=4, cyclic=True, field=QQ)   # self-injective
    with pytest.raises(QuiverlabError):
        homological_string_test(sinj)                  # knitter unsupported
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError:
  quiverlab.modules.string_homological`

- [ ] **Step 3: Implement**

```python
# src/quiverlab/modules/string_homological.py
"""The homological string-algebra test (Plan 59 / R34).

Suarez-Alvarez (Alg. Rep. Theory 26 (2023) 1759-1772, arXiv:2105.02948): among
REPRESENTATION-FINITE algebras, A is a STRING algebra iff the middle term of EVERY
extension of indecomposable modules has at most two indecomposable summands. The
quantifier is over all Ext^1(X, Y) classes, NOT just AR sequences.

Honest semi-decision (see the Plan-59 design note); verdict in
{"not_string", "string_over_ground_field", "inconclusive"}:
  * REFUTE ("not_string", k-bar-sound): a >= 3-summand middle E. "String" is a property
    of (Q, I) and is field-independent, and a >= 3 decomposition base-changes up to
    k-bar, so a ground-field witness certifies non-string. Only the terminal split step
    is char-robust; the REFUTE PIPELINE (knit / almost_split / decompose) needs char 0
    or char > dim, so refute runs over QQ (small primes only on tiny algebras). Checked
    classes: the split (2), every AR sequence, and a basis of each Ext^1.
  * CONFIRM ("string_over_ground_field", INCONCLUSIVE w.r.t. the k-bar theorem): finite
    ground field only -- exhaustive End(X) x End(Y)-orbit enumeration of every
    Ext^1(X, Y), each middle certified <= 2 (needs char > dim E or char 0). Over k-bar
    there may be more indecomposables / non-rational classes, so this does NOT prove the
    theorem's "string"; the definitive string answer is the syntactic is_string.
  * else INCONCLUSIVE, loud reason.

CONTRACT (one thing): RAISES iff a >= 3 witness is found AND is_string(A) is True (one
engine is a bug -- k-bar-sound). Otherwise RETURNS the verdict; a
"string_over_ground_field" with is_string False records a kbar_gap_note (never raised).
Refuses loudly for rep-infinite / self-injective / presentation-less A (out of scope --
the theorem is stated for rep-finite algebras)."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.recognizers import is_string
from quiverlab.modules.ar import (_combine_matrices, _ext1_data,
                                   almost_split_sequence, knit_ar_quiver)
from quiverlab.modules.decompose import decompose
from quiverlab.modules.ses import ShortExactSequence
from quiverlab.modules.morphism import ModuleHom
from quiverlab.modules.yoneda import baer_extension

_REFS = ("suarez_alvarez", "butler_ringel", "assem_book")


@dataclass(frozen=True)
class StringHomologyVerdict:
    verdict: str                 # not_string | string_over_ground_field | inconclusive
    is_string: bool              # the P38 syntactic arbiter (definitive k-bar "string")
    witness: "dict | None"       # unified schema (see _mk_witness)
    kbar_gap_note: "str | None"  # confirm-vs-is_string mismatch (never raised)
    checked: dict
    reason: "str | None"
    references: tuple = _REFS


def _require_in_scope(A):
    if A.quiver is None:
        raise QuiverlabError(
            "homological_string_test: needs a quiver-presented algebra",
            hint="the extension/middle-term test walks the knitted indecomposables")
    ar = knit_ar_quiver(A)
    if not ar.is_complete or ar.status != "complete":
        raise QuiverlabError(
            "homological_string_test: A is not representation-finite in scope "
            f"(AR knit status={ar.status!r}); the Suarez-Alvarez characterization is "
            "stated for rep-finite algebras",
            hint="self-injective and rep-infinite inputs are refused (knit "
                 "unsupported/budget) -- this test does not extend past them")
    return ar


def _middle_of_class(A, X, Y, coeffs):
    """The middle E of the Ext^1(X, Y) class with coordinates `coeffs` over the cocycle
    basis (the almost_split_sequence idiom). Returns (E, ses) with the exact sequence
    0 -> Y -> E -> X -> 0; baer_extension self-certifies `coeffs` is a genuine class."""
    dom = A.domain
    cocycle_mats, _cob, terms, dmats = _ext1_data(A, X, Y)
    width = len(dmats[1][0]) if (len(dmats) > 1 and dmats[1] and dmats[1][0]) else 0
    f = _combine_matrices(cocycle_mats, coeffs, Y.dim, width, dom)
    seq = baer_extension(X, Y, f, terms, dmats)          # 0 -> Y -> E -> X -> 0
    seq.assert_exact()
    E = seq.modules[1]
    ses = ShortExactSequence(
        ModuleHom(seq.modules[0], E, seq.maps[0], check=False),
        ModuleHom(E, seq.modules[2], seq.maps[1], check=False))
    return E, ses


def _summands_at_least(E, k):
    """(count, dimvecs) via decompose; a >= k split is char-robust (a found idempotent
    certifies), so surface `>= k` even if the LEAF-locality certificate is char-blocked:
    catch the decompose refusal and fall back to a bounded Fitting-split search that
    reports the number of pieces FOUND (a lower bound). See Adjust-to-reality."""
    ...


def _mk_witness(X, Y, klass, coeffs, E, cnt, dvs):
    """Unified witness schema across BOTH refute sub-routes. `klass` in
    {"ar_socle","basis"}; `coeffs` is None on the AR route (the almost-split socle class
    is not a coordinate vector), the unit/combination vector on the basis route."""
    return {"X": X.dimension_vector(), "Y": Y.dimension_vector(),
            "class": klass, "coeffs": (None if coeffs is None else _fmt(coeffs)),
            "E_dimvec": E.dimension_vector(),
            "summand_count": cnt, "summand_dimvecs": dvs}


def homological_string_test(A, budget_classes=4096, budget_pairs=4096):
    ar = _require_in_scope(A)
    syn = is_string(A)
    dom = A.domain
    q = getattr(dom, "order", None)          # finite field size, or None (infinite)
    indec = [v["module"] for v in ar.vertices]
    route, n_classes = [], 0
    witness = None

    # -- REFUTE pass: split classes are trivially 2; check every AR sequence middle and a
    #    basis of each Ext^1. Any >= 3 => not_string (k-bar-sound). --
    # (a) AR sequences (end Y; ses.L = tau Y, ses.M = E, ses.N = Y)
    route.append("ar")
    for Y in indec:
        try:
            ses = almost_split_sequence(Y)   # 0 -> tauY -> E -> Y -> 0
        except QuiverlabError:
            continue                         # projective end / char-undecidable -> skip
        cnt, dvs = _summands_at_least(ses.M, 3)
        if cnt >= 3:
            witness = _mk_witness(Y, ses.L, "ar_socle", None, ses.M, cnt, dvs); break
    # (b) basis classes of every ordered pair 0 -> Y -> E -> X -> 0
    if witness is None:
        route.append("basis")
        for X in indec:
            for Y in indec:
                d = _ext1_dim(A, X, Y)       # dim Ext^1(X, Y) via _ext1_data / ext_dims
                for e in _unit_vectors(d, dom):
                    n_classes += 1
                    E, _ = _middle_of_class(A, X, Y, e)
                    cnt, dvs = _summands_at_least(E, 3)
                    if cnt >= 3:
                        witness = _mk_witness(X, Y, "basis", e, E, cnt, dvs); break
                if witness or n_classes > budget_classes:
                    break
            if witness or n_classes > budget_classes:
                break

    checked = {"n_indec": len(indec), "route": route, "field": str(dom),
               "classes": n_classes}
    if witness is not None:
        # THE one loud disagreement: a k-bar-sound non-string witness vs is_string True.
        if syn is True:
            raise QuiverlabError(
                "homological_string_test: found an extension of indecomposables with a "
                ">= 3-summand middle, so A is NOT a string algebra, yet is_string(A) is "
                "True -- one engine is a bug",
                hint=f"witness middle dim vector {witness['E_dimvec']} "
                     f"({witness['summand_count']} summands)")
        return StringHomologyVerdict("not_string", syn, witness, None, checked, None)

    # -- CONFIRM pass (finite ground field only): exhaustive orbit enumeration. --
    if q is None:
        return StringHomologyVerdict(
            "inconclusive", syn, None, None, checked,
            "infinite ground field: the quantifier over ALL Ext^1 classes cannot be "
            "exhausted; checked split/AR/basis only (no >= 3 witness). Defers to is_string.")
    ok, reason = _exhaustive_confirm(A, indec, budget_classes)   # every middle <= 2 ?
    if ok:
        route.append("exhaustive")
        checked["classes"] = q
        gap = (None if syn else
               "string_over_ground_field but is_string(A) is False: the finite ground "
               f"field GF({q}) did not witness the non-stringness the k-bar theorem "
               "guarantees (more indecomposables / non-rational classes over k-bar)")
        return StringHomologyVerdict("string_over_ground_field", syn, None, gap, checked, None)
    return StringHomologyVerdict("inconclusive", syn, None, None, checked, reason)


def string_homological_block(A):
    """The `string_homological` algebra-only compute kind (Task 5)."""
    try:
        v = homological_string_test(A)
    except QuiverlabError as exc:
        return {"error": str(exc), "references": list(_REFS)}
    return {"verdict": v.verdict, "is_string": v.is_string, "witness": v.witness,
            "kbar_gap_note": v.kbar_gap_note, "checked": v.checked, "reason": v.reason,
            "references": list(v.references)}
```

**Adjust to reality (Task 1):**
- **`_ext1_data` / `_combine_matrices` are `modules/ar.py` internals** — read
  `almost_split_sequence` (lines ~479–482) for the exact call shape and reuse it
  verbatim; `modules/` code is allowed to use `modules/` internals (this is not
  reaching into `engine.*`). Prefer the public `ext_dims(A, X, Y, 1, with_reps=True,
  interpret=True)` for the **basis** middle terms (its `interpretation` gives the Yoneda
  `middle` modules directly) and fall back to `_ext1_data` + `_combine_matrices` only
  for **non-basis** coordinate vectors in the finite-field exhaustive pass. Whichever you
  use, the middle `E` must be `is_isomorphic`-comparable and `decompose`-able.
- **`_summands_at_least(E, k)` is the char scope crux (W5).** The refute PIPELINE is
  **not** char-robust end-to-end: `knit_ar_quiver` / `almost_split_sequence` /
  `decompose` all lean on `char 0` or `char > dim` (the trace-form radical + Fitting
  locality). Only the **terminal split step** is char-robust (a found idempotent is exact
  linear algebra). So `_summands_at_least` catches a `decompose(E)` refusal and falls
  back to a **bounded Fitting-split search** (the `_try_split` primitive `decompose`
  uses) that reports the number of pieces *found* — a `≥ k` split is a valid witness in
  any characteristic. But because the pipeline that PRODUCES the middles is char-gated,
  **refute runs over QQ in practice** (`char 0`); over small primes the pipeline is
  blocked upstream except on tiny algebras where `char > dim` holds throughout. State
  this scope in the docstring; do not claim "char-robust in any characteristic" for the
  pipeline — only for the final split.
- **`_exhaustive_confirm` orbit reduction.** Enumerating all `qᵈ` classes is only feasible
  for small `q` and small `d`; use `end_action_on_ext1(X, Y)` to enumerate
  `End(X)×End(Y)`-orbit representatives (+ scalar) and cap at `budget_classes`. On budget
  overflow return `inconclusive` (reason `"budget"`), never a partial confirm. If ANY
  middle in the confirm pass is char-undecidable (`decompose` refuses and no ≤2 split is
  found) return `inconclusive` (reason `"char-undecidable"`), naming the pair.
- **`kD₄` witness (VERIFIED LIVE 2026-08-07 on `dev`, `QUIVERLAB_NO_NUMBA=1`).** The
  explicit subspace-orientation `Quiver([0,1,2,3], {"a":(1,0),"b":(2,0),"c":(3,0)})` over
  QQ has `dim = 7`, `is_string = False`, `is_special_biserial = False`;
  `knit_ar_quiver` is `status="complete"` with **12** indecomposables. TWO AR sequences
  carry a 3-summand middle; the headline pin is the one ending at the module `M` with
  dim vector `{0:2, 1:1, 2:1, 3:1}`, whose AR middle `E` has dim vector
  `{0:3, 1:1, 2:1, 3:1}` and decomposes as **`P₁ ⊕ P₂ ⊕ P₃`** (the three arm-projectives
  `{0:1,1:1}`, `{0:1,2:1}`, `{0:1,3:1}`) — exactly 3 summands. The second (unnamed)
  3-middle ends at `I₀ = {0:1,1:1,2:1,3:1}`, `E = {0:3,2,2,2}`, three summands each with
  one arm zeroed. The `@lit` test pins the `E` dim vector `{0:3,1:1,2:1,3:1}` and
  `summand_count == 3`; if `almost_split_sequence` surfaces the `I₀` mesh first that is
  also `≥ 3` and still passes — but do NOT weaken the pin, fix the summand counter if it
  under-reports.
- **`inconclusive`/`string_over_ground_field` carry no agreement claim.** Only a returned
  `not_string` implies `is_string is False` (else it would have raised). On `kA_n` over
  `GF(5)` the confirm pass closes to `verdict="string_over_ground_field"` (NOT a bare
  `"string"`); on `kA_n` over `QQ` it is `inconclusive` (infinite field). **Never**
  silently base-change QQ→finite to manufacture a confirm — the honest QQ answer is
  `inconclusive`, deferring to `is_string`.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/string_homological.py \
        tests/modules/test_string_homological_p59.py
git commit -m "feat(modules): homological_string_test -- Suarez-Alvarez extension/middle-term semi-decision (three-valued; char-robust >=3 witness, finite-field exhaustive confirm), kD4 3-middle discriminator"
```

---

### Task 2: the discriminating battery + QPA (R34)

The headline oracle: over the rep-finite zoo, a returned `not_string` (a `≥ 3` witness,
`k̄`-sound) always has `is_string False` — and if `is_string` were `True` the function
would already have raised. `string_over_ground_field`/`inconclusive` defer to `is_string`;
a `k̄`-gap (confirm-vs-`is_string` mismatch) is recorded, never raised.

**Files:**
- Test: `tests/modules/test_string_homological_battery_p59.py` (deep),
  `tests/qpa/test_string_homological_qpa.py` (qpa)
- **No `src/` change** — the QPA crosscheck is a **direct session parity call** (below);
  it does NOT go through `qpa/crosscheck.py::crosscheck` (which has no recognizer verb —
  an unknown `what` raises `QuiverlabError` at its dispatch tail; verified on `dev`).

**Interfaces:**
- Consumes: `homological_string_test`, `invariants/recognizers.py::is_string` /
  `is_special_biserial`, a spread of rep-finite algebras (string AND non-string),
  `qpa/session.py` (`should_skip_qpa`, `require_gap`, `run`, `libgap_handle`),
  `qpa/scripts.py::quiver_and_algebra_script(A)` (the same builder
  `crosscheck_symmetric` uses to hand a presented `kQ/I` to QPA).
- Produces: the discriminating battery + honest QPA scope (recognizer parity via a
  direct `IsSpecialBiserialAlgebra` call + a no-homological-surface guard).

- [ ] **Step 1: Write the failing battery**

```python
# tests/modules/test_string_homological_battery_p59.py
"""Discriminating battery (Plan 59 / R34). A returned "not_string" (a >= 3-summand
middle) is k-bar-sound and only occurs when is_string is False (else the function
raises). "string_over_ground_field"/"inconclusive" defer to is_string; a k-bar-gap is
recorded, never raised. String side: kA_n / gentle kA_n(ab) -> no >= 3 witness.
Non-string rep-finite side: Dynkin D/E path algebras (a valency->=3 vertex -> a >= 3 AR
mesh) -> not_string with a witness."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import is_string
from quiverlab.modules.string_homological import homological_string_test

xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


def _path(vertices, arrows):
    return Quiver(vertices, arrows).algebra(relations=[], field=QQ)


STRING = [                                   # rep-finite string algebras
    _path([1, 2], {"a": (1, 2)}),                          # kA2
    _path([1, 2, 3], {"a": (1, 2), "b": (2, 3)}),          # kA3
    Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ),
]
# rep-finite, NOT string (a valency->=3 vertex). VERIFIED LIVE: kD4 (subspace) knits
# complete (12 indec) with a 3-summand AR middle. Build D5/E6 via families.dynkin.
from quiverlab.families.dynkin import dynkin_quiver
NOT_STRING = [
    _path([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)}),          # kD4 (subspace)
    dynkin_quiver("D5").algebra(relations=[], field=QQ),                   # kD5
    dynkin_quiver("E6").algebra(relations=[], field=QQ),                   # kE6
]


@xeng
@pytest.mark.parametrize("A", STRING)
def test_string_side_no_false_refutation(A):
    assert is_string(A) is True
    v = homological_string_test(A)               # cannot raise: no >=3 witness on strings
    assert v.verdict != "not_string"             # never a false >= 3 witness


@xeng
@pytest.mark.parametrize("A", NOT_STRING)
def test_not_string_side_finds_a_witness(A):
    assert is_string(A) is False                 # a valency->=3 vertex: not special-biserial
    v = homological_string_test(A)
    assert v.verdict == "not_string" and v.witness["summand_count"] >= 3


@lit
def test_refute_side_raises_on_a_contradiction_never_returns_it():
    # The one loud disagreement: a >=3 witness with is_string True would RAISE. Since our
    # zoo has none, every returned verdict is self-consistent: not_string => is_string False.
    from quiverlab.errors import QuiverlabError  # noqa: F401 (documents the raise contract)
    for A in STRING + NOT_STRING:
        v = homological_string_test(A)           # returns (no contradiction in the zoo)
        if v.verdict == "not_string":
            assert v.is_string is False
```

- [ ] **Step 2: Non-string set is LIVE-VERIFIED.** `kD₄` (subspace orientation) knits
  complete (12 indecomposables, verified) with a 3-summand AR middle; `kD₅`/`kE₆` via
  `families.dynkin.dynkin_quiver` are rep-finite Dynkin with a valency-3 vertex (same
  mesh mechanism). If a rep-finite special-biserial **non-monomial** algebra (the subtle
  discriminator) is added, keep it only where `homological_string_test` **returns
  definitive**; where it is `inconclusive`, assert it **defers to** `is_string`
  (`verdict == "inconclusive"`, no agreement asserted) and record it as an honest-scope
  example. Never assert a `not_string` the engine did not produce. NOTE: the commutative
  `[2,2]` toupie (Task 3) is special-biserial-non-monomial and — verified live — knits...
  (leave rep-finiteness of that specific example to the implementer; only use it in the
  battery if `knit_ar_quiver` returns `status="complete"`).

- [ ] **Step 3: QPA scope (W2 — the real crosscheck).** `A.crosscheck("is_special_biserial")`
  does NOT exist — `qpa/crosscheck.py::crosscheck` has no recognizer verb (unknown `what`
  raises). Use a **direct session parity call** (no `src/` change), mirroring
  `crosscheck_symmetric`'s internals: build the QPA algebra with
  `scripts.quiver_and_algebra_script(A)` and run `IsSpecialBiserialAlgebra(A)`, comparing
  to our `is_special_biserial(A)`. Plus the standing no-surface guard (P35 precedent):

```python
# tests/qpa/test_string_homological_qpa.py
"""QPA scope for R34 (Plan 59). QPA recognizes special-biserial/gentle but has NO
homological string test (no extension/middle-term verb). Crosschecks: a DIRECT-session
IsSpecialBiserialAlgebra parity call (crosscheck.py has no recognizer verb) + a standing
guard that FAILS if such a surface ever appears."""
import pytest
from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import is_special_biserial
from quiverlab.qpa import scripts, session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_has_no_homological_string_surface():
    lg = session.libgap_handle()
    for name in ("MiddleTermsOfExtensions", "HomologicalStringTest"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), \
            f"QPA now ships {name} -- add a real crosscheck (honest scope changed)"


def test_special_biserial_parity_direct_session():
    session.require_gap()
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ)
    base = scripts.quiver_and_algebra_script(A)            # binds `A` in the GAP session
    qpa = bool(session.run(base + "\nIsSpecialBiserialAlgebra(A);"))
    assert qpa == is_special_biserial(A)                   # gentle kA3/(ab): both True
```

- [ ] **Step 4: Run** the deep battery + `-m qpa` live. Expected: PASS (the `[qpa]` extra
  is installed in the venv). Confirm `IsSpecialBiserialAlgebra` exists in this QPA
  version first (`NamesGVars()` probe); if a QPA version lacks it, guard that one call
  with a documented `pytest.skip` naming the missing verb (never a silent pass).
- [ ] **Step 5: Commit**

```bash
git add tests/modules/test_string_homological_battery_p59.py \
        tests/qpa/test_string_homological_qpa.py
git commit -m "test(modules,qpa): R34 discriminating battery (kD4/kD5/kE6 live-verified 3-middle witnesses; refute-side raises on contradiction) + direct-session IsSpecialBiserialAlgebra parity + no-homological-surface guard"
```

---

### Task 3: `ToupieAlgebra` constructor + `is_toupie` recognizer (R35)

**Files:**
- Create: `src/quiverlab/families/toupie.py`
- Modify (exports): `src/quiverlab/families/__init__.py` (export `ToupieAlgebra`,
  `is_toupie`) **and** `src/quiverlab/__init__.py` (top-level export + `__all__`, beside
  the `BrauerGraph`/`BrauerGraphAlgebra` line at ~26/55 — the webapp discovery reads the
  top-level export).
- Modify (CATALOG): `src/quiverlab/families/discover.py` — add a `FamilyInfo` entry.
- Modify (BOTH discovery skip-sets — verified on `dev`; `discover.py` has NEITHER):
  `src/quiverlab/hpc/spec.py:328` — the inline skip tuple `("zoo",
  "BrauerGraphAlgebra")` in `_iter_families()` gains `"ToupieAlgebra"`; **and**
  `webapp/server/catalog.py:27` — the `_NON_FORM_FAMILIES` frozenset gains
  `"ToupieAlgebra"`. Both keep the scalar form-builder from introspecting the list-arg
  constructor; it is surfaced as a GUI preset (Task 5) instead.
- Test: `tests/families/test_toupie_p59.py`

**Interfaces:**
- Consumes: `combinat/quiver.py::Quiver` + `Quiver.algebra(relations, field)` (the
  `QuantumCI` non-monomial template — relations are exact grammar tokens over the
  branch-path names), `families/nakayama.py` loud-hint validation template,
  `families/trivial_extension.py` dimension-certificate idiom,
  `invariants/recognizers.py::_in_degree`/`_out_degree`.
- Produces:
  ```python
  def ToupieAlgebra(branches, relations=None, field=None) -> Algebra
      # `branches` = list of branch lengths [l_1,..,l_a] (each l_i >= 1). Builds the
      # toupie quiver: source `s`, sink `t`, branch i a fresh directed path
      # s -> v_{i,1} -> ... -> v_{i,l_i-1} -> t of length l_i (a length-1 branch is a
      # single arrow s->t). `relations` (optional) = grammar-token strings over the
      # branch-path tokens p0..p_{a-1} (monomial truncations and/or NON-MONOMIAL
      # cross-branch linear combos). Presented via Quiver.algebra.
      # Per-instance dimension certificate for the RELATION-FREE case (VERIFIED LIVE):
      #   dim = 2 + sum_i (l_i - 1)  +  sum_i l_i*(l_i+1)/2      (|Q_0| + #nonzero paths)
      #   [1,1]->4  [1,1,1]->5  [2]->6  [1,2]->7  [2,2]->10   (all confirmed on dev)
      # (a-Kronecker: l_i=1 all -> dim = 2 + a; single branch l=n -> dim = kA_{n+1}).
      # With relations: certify finite-dimensionality (Quiver.algebra raises otherwise)
      # + an OPTIONAL user `expected_dim` (loud on mismatch). a-Kronecker = ToupieAlgebra
      # ([1]*a). Loud on a < 1, a bad branch length, a length-1 cross-branch relation
      # (not in rad^2), or a failed certificate.
  def is_toupie(A) -> bool
      # graph-shape recognizer (the "trivial graph" of R35): the quiver is CONNECTED and
      # ACYCLIC, has a UNIQUE source (in-deg 0), a UNIQUE sink (out-deg 0), and EVERY
      # other vertex has in-deg 1 AND out-deg 1. CONNECTED + ACYCLIC are LOAD-BEARING:
      # without them a (path U oriented-cycle) disjoint quiver passes the degree checks
      # (unique source/sink, others in/out-deg 1 -- VERIFIED) yet is disconnected with a
      # cycle => kQ infinite-dimensional. Loud on a presentation-less algebra.
  def toupie_branch_count(A) -> int
      # a_branch = out_deg(source) == in_deg(sink) (# branches); asserts the equality.
      # This is NOT the sl_a `a` -- see toupie_direct_arrow_count (Task 4).
  def toupie_direct_arrow_count(A) -> int
      # the sl_a `a`: the number of arrows DIRECTLY source->sink (length-1 branches).
      # a_direct <= a_branch, with equality iff every branch has length 1 (a-Kronecker).
  def toupie_block(A) -> dict                    # the `toupie` algebra-only compute kind (Task 5)
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_toupie_p59.py
"""Toupie algebras (Plan 59 / R35, Artenstein-Lanzilotta-Solotar 2020). Self-cert:
dimension certificate; recognizer accepts every constructed toupie and rejects
non-toupie shapes. Literature: the a-Kronecker HH = [1, a^2-1, 0, ..]. sl_a: char 0."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.families.toupie import (ToupieAlgebra, is_toupie, toupie_branch_count)

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@selfcert
@pytest.mark.parametrize("branches, dim", [
    ([1, 1], 4),            # 2-Kronecker: 2 + 2
    ([1, 1, 1], 5),         # 3-Kronecker: 2 + 3
    ([2], 6),               # single length-2 branch = kA3
    ([2, 2], 10),           # two length-2 branches, no rel: 4 verts + 6 paths
    ([1, 2], 7),            # mixed lengths: 3 verts + 4 paths
])
def test_relation_free_dimension_certificate(branches, dim):
    A = ToupieAlgebra(branches, field=QQ)
    assert A.dim == dim
    assert is_toupie(A) is True
    assert toupie_branch_count(A) == len(branches)


@selfcert
def test_recognizer_rejects_non_toupie():
    kD4 = Quiver([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)}).algebra(
        relations=[], field=QQ)                            # centre has in-deg 3 -> not a toupie
    assert is_toupie(kD4) is False
    line = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=QQ)
    assert is_toupie(line) is True                         # single branch = a toupie (a=1)


@selfcert
def test_recognizer_rejects_disconnected_or_cyclic_trap():
    # VERIFIED LIVE: (path U oriented-cycle) passes the DEGREE checks -- unique source
    # [0], unique sink [1], all other vertices in/out-deg 1 -- yet is disconnected AND has
    # a cycle, so kQ is infinite-dimensional. is_toupie MUST reject it (connected+acyclic).
    trap = Quiver([0, 1, "c0", "c1"],
                  {"p": (0, 1), "e": ("c0", "c1"), "f": ("c1", "c0")}).algebra(
        relations=["e*f", "f*e"], field=QQ)                # kill the cycle to stay f.d.
    assert is_toupie(trap) is False                        # disconnected + (pre-truncation) cyclic


@selfcert
def test_non_monomial_cross_branch_relation():
    # two length-2 branches p1 = x1*x2, p2 = y1*y2; the commutative toupie p1 - p2 = 0.
    A = ToupieAlgebra([2, 2], relations=["p0 - p1"], field=QQ)   # branch-path token names
    assert is_toupie(A) is True
    assert A.dim == 9                                      # 10 (free) - 1 (identifies the two 0->w paths)


@selfcert
def test_loud_refusals():
    with pytest.raises(QuiverlabError):
        ToupieAlgebra([], field=QQ)                        # no branches
    with pytest.raises(QuiverlabError):
        ToupieAlgebra([1, 0], field=QQ)                    # branch length 0
    from quiverlab.families import zoo
    # a structure-constant (presentation-less) algebra: is_toupie refuses loudly
    # (choose one from the zoo that carries no quiver) -- see Adjust to reality.
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError:
  quiverlab.families.toupie`

- [ ] **Step 3: Implement.** Build the quiver (fresh internal vertices per branch;
  a length-1 branch is a single `s->t` arrow), name branch arrows deterministically,
  expose the full branch-paths `p_0..p_{a-1}` as grammar tokens the `relations=` strings
  can reference (the `QuantumCI` token idiom), present via `Quiver.algebra`, and certify
  the dimension. The recognizer is `_in_degree`/`_out_degree` bookkeeping PLUS a
  connectivity (union-find on the underlying graph) and acyclicity (topological-sort /
  no directed cycle) check.

**Adjust to reality (Task 3):**
- **Branch-path token names.** `relations=` strings must reference the full `0→ω`
  branch-paths. Decide a stable spelling — either explicit arrow products
  (`"x1*x2 - y1*y2"`) or convenience path tokens `p0,p1,…` the constructor expands to
  the arrow product before handing to `Quiver.algebra`. Pick one, document it, and make
  `test_non_monomial_cross_branch_relation` use exactly that spelling. A cross-branch
  linear relation is admissible only when the involved branches have length `≥ 2`
  (a relation must lie in `rad²`; parallel length-1 arrows sit in `rad`, not `rad²`) —
  **validate this loudly** (a length-1 cross relation is a user error, not a toupie).
- **Dimension certificate (W7).** Relation-free: the closed form
  `2 + Σ(lᵢ−1) + Σ lᵢ(lᵢ+1)/2`, **all verified live on `dev`**: `[1,1]→4`,
  `[1,1,1]→5`, `[2]→6`, `[1,2]→7`, `[2,2]→10` (NOT 8 — the old sketch was wrong). With
  relations: `Quiver.algebra` computes `dim` (raises `NotFiniteDimensionalError` if not
  finite-dim — a wiring bug, since a genuine toupie is finite-dim); accept an optional
  `expected_dim=` and raise loudly on mismatch (the `TrivialExtension`/`Brauer`
  precedent; verified: the commutative `[2,2]` with `p0−p1` is `dim 9`). Do **not**
  invent a general with-relations closed form.
- **`is_toupie` must check CONNECTED + ACYCLIC (M1), not just degrees.** Compute
  `_in_degree`/`_out_degree`: unique in-deg-0 source, unique out-deg-0 sink, every other
  vertex in/out-deg 1 — AND the underlying graph connected AND the quiver acyclic (no
  oriented cycle). The (path ⊔ oriented-cycle) trap passes the degree checks (VERIFIED)
  but must be rejected. On a presentation-less algebra refuse loudly (`QuiverlabError`,
  `hint=`; mirror `_require_quiver`) — never a silent `False`.
- **`toupie_branch_count` vs `toupie_direct_arrow_count`.** `toupie_branch_count` =
  `out_deg(source)` (assert `== in_deg(sink)`, raise on inequality). `toupie_direct_arrow_count`
  = the number of arrows whose source IS the source and target IS the sink (length-1
  branches) — **this** is the `sl_a` `a` (Task 4), NOT the branch count.
- **CATALOG + BOTH discovery skip-sets (W3).** Add `FamilyInfo("ToupieAlgebra",
  "ToupieAlgebra([l1,..,la], relations=...)", "general", ("alsolotar_toupie",
  "assem_book"), "Toupie: unique source/sink, a parallel branches; admissible I.")` to
  `families/discover.py::CATALOG`. `ToupieAlgebra` takes a **list** (non-scalar), so add
  `"ToupieAlgebra"` to **both** real skip sites (there is NO `_iter_families` in
  `discover.py`): the inline tuple `("zoo", "BrauerGraphAlgebra")` in
  `hpc/spec.py:328::_iter_families`, and the `_NON_FORM_FAMILIES` frozenset in
  `webapp/server/catalog.py:27`. Surface it as a GUI **preset** (Task 5).

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/families/toupie.py src/quiverlab/families/__init__.py \
        src/quiverlab/__init__.py src/quiverlab/families/discover.py \
        src/quiverlab/hpc/spec.py webapp/server/catalog.py \
        tests/families/test_toupie_p59.py
git commit -m "feat(families): ToupieAlgebra(branches, relations) presented kQ/I + is_toupie (connected+acyclic graph-shape) + branch/direct-arrow counts; dim certificate, non-monomial cross-branch relations, both discovery skip-sets, loud refusals"
```

---

### Task 4: toupie HH closed-form oracle family + `sl_a ⊆ HH¹` (R35)

**Files:**
- Modify: `src/quiverlab/families/toupie.py` (fold the HH/`sl_a` logic into
  `toupie_block`; add a small `toupie_sl_a_lower_bound(A)` helper)
- Test: `tests/families/test_toupie_hh_p59.py` (deep)

**Interfaces:**
- Consumes: `Algebra.hochschild_cohomology(top, engine=…, auto_cs=…)` (the shipped bar/CS
  dispatch), `toupie_direct_arrow_count`, `is_toupie`. Char-0 (`QQ`) for the `sl_a`
  inclusion.
- Produces: the HH oracle family + the char-0 `sl_a` lower bound (on the DIRECT-arrow `a`).

**sl_a correction (found via live verification, 2026-08-07 — supersedes the review's M2
premise).** `a` in ALS Thm 6.5 is the number of **direct arrows `source→sink`** (length-1
branches), NOT the branch count. Live evidence: the commutative `[2,2]` toupie (2
branches, 0 direct arrows, relation `a1a2 − b1b2`) has `HH^• = [1,0,0,0,0]` — so
`dim HH¹ = 0`, which would VIOLATE a `branch-count`-based `a²−1 = 3` bound but is exactly
right for `a_direct = 0` (no `sl_a`). The relation-free mixed toupie (2 direct arrows + 1
length-2 branch) has `HH^• = [1,6,0,0]`, and `6 ≥ a_direct²−1 = 3` ✓ (strict). So the
`sl_a` lower bound uses `a_direct = toupie_direct_arrow_count`, and only claims content
when `a_direct ≥ 2`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_toupie_hh_p59.py
"""Toupie Hochschild cohomology (Plan 59 / R35). Literature (SAFE pin): the a-Kronecker
HH^* = [1, a^2-1, 0, ...] (ALS 2020; the SAME VALUE the m-Kronecker engine test pins,
recomputed here independently over QQ through Algebra.hochschild_cohomology at top=4).
Cross-engine: bar == CS degreewise where bar survives (non-hereditary toupies blow the
bar window past low degree). sl_a (char 0 only): a = # DIRECT source->sink arrows;
dim HH^1 >= a^2-1, equality on the a-Kronecker. ALL VALUES BELOW VERIFIED LIVE on dev."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.toupie import (ToupieAlgebra, toupie_direct_arrow_count,
                                       toupie_sl_a_lower_bound)

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("a", [2, 3, 4])
def test_a_kronecker_hh_is_closed_form(a):
    A = ToupieAlgebra([1] * a, field=QQ)                   # a-Kronecker, hereditary
    hh = A.hochschild_cohomology(4).dims                   # VERIFIED: a=2->[1,3,..], a=3->[1,8,..]
    assert hh == [1, a * a - 1, 0, 0, 0]                   # ALS 2020: [1, a^2-1, 0, ..]


@lit
@pytest.mark.parametrize("a", [2, 3, 4])
def test_a_kronecker_hh1_is_dim_sl_a(a):
    # char 0: HH^1 = sl_a exactly on the a-Kronecker (all a arrows are direct).
    A = ToupieAlgebra([1] * a, field=QQ)
    hh1 = A.hochschild_cohomology(1).dims[1]
    assert toupie_direct_arrow_count(A) == a               # all branches length 1
    assert hh1 == a * a - 1 == toupie_sl_a_lower_bound(A)  # dim sl_a = a^2 - 1


@xeng
def test_bar_equals_cs_on_a_kronecker_and_empty_ideal():
    # CS handles the EMPTY ideal (VERIFIED): a-Kronecker cs == cs(auto_cs) == bar == auto.
    A = ToupieAlgebra([1, 1, 1], field=QQ)                 # dim 5, hereditary
    bar = A.hochschild_cohomology(4, engine="bar").dims
    cs = A.hochschild_cohomology(4, engine="cs", auto_cs=True).dims
    cs_plain = A.hochschild_cohomology(4, engine="cs").dims
    assert bar == cs == cs_plain == [1, 8, 0, 0, 0]        # VERIFIED live


@xeng
def test_bar_equals_cs_on_a_NON_hereditary_toupie():
    # M2: a genuine NON-monomial cross-branch toupie. Bar blows up past deg 2 (36864x4608
    # > max_cells at d^3), so compare where bar SURVIVES: top=2. VERIFIED: both [1,0,0].
    A = ToupieAlgebra([2, 2], relations=["p0 - p1"], field=QQ)   # commutative square, dim 9
    bar = A.hochschild_cohomology(2, engine="bar").dims
    cs = A.hochschild_cohomology(2, engine="cs", auto_cs=True).dims
    assert bar == cs == [1, 0, 0]                          # VERIFIED live (a_direct=0 -> no sl_a)
    full = A.hochschild_cohomology(4, engine="cs", auto_cs=True).dims
    assert full == [1, 0, 0, 0, 0]                         # CS past the bar window (VERIFIED)


@xeng
def test_sl_a_lower_bound_with_direct_arrows():
    # char 0: dim HH^1 >= dim sl_a = a_direct^2 - 1, a_direct = # direct source->sink arrows.
    # mixed toupie: 2 DIRECT arrows + one length-2 branch (relation-free). VERIFIED [1,6,..].
    A = Quiver([0, "v", 1],
               {"d0": (0, 1), "d1": (0, 1), "b1": (0, "v"), "b2": ("v", 1)}
               ).algebra(relations=[], field=QQ)
    assert toupie_direct_arrow_count(A) == 2
    hh1 = A.hochschild_cohomology(3, engine="cs", auto_cs=True).dims[1]
    assert hh1 == 6 >= toupie_sl_a_lower_bound(A) == 2 * 2 - 1   # 6 >= 3 (strict inclusion)


def test_sl_a_refused_off_char_zero():
    from quiverlab import GF
    from quiverlab.errors import QuiverlabError
    A = ToupieAlgebra([1, 1], field=GF(7))
    with pytest.raises(QuiverlabError):
        toupie_sl_a_lower_bound(A)                         # sl_a inclusion is char-0 only
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** `toupie_sl_a_lower_bound(A)` = `a²−1` where
  `a = toupie_direct_arrow_count(A)` (the # of direct source→sink arrows — NOT the branch
  count), **hard-gated to char 0** (raise `QuiverlabError` off char 0 — the ALS Thm 6.5
  hypothesis is `k = ℂ`/char 0; the critic scoped this explicitly). Fold the HH forms +
  the `sl_a` line into `toupie_block` (Task 5).

**Adjust to reality (Task 4):**
- **`a` is the DIRECT-arrow count — `# PIN` against Thm 6.5.** ALS phrase it "the number
  of arrows from `0` to `ω`"; the live evidence (commutative `[2,2]`: 2 branches, 0 direct
  arrows, `HH¹ = 0`) forces `a = a_direct`, not the branch count. Confirm the exact wording
  in Thm 6.5 before the `sl_a` test asserts; the a-Kronecker (all branches length 1,
  `a_direct = a_branch`) is unaffected either way and is the clean pin.
- **The `a`-Kronecker `[1, a²−1, 0, …]` is the ONLY transcribable closed-form pin, and it
  is a NEW QQ recomputation (W6).** Do NOT claim to "reuse" the P29 `test_literature_p29`
  machinery — that test pins the SAME VALUE for the `m`-Kronecker via the engine path;
  here we independently recompute over QQ through `Algebra.hochschild_cohomology` at
  `top=4` (verified: `a=2 → [1,3,0,0,0]`, `a=3 → [1,8,0,0,0]`). Cross-reference the P29
  value in a comment; the citation edit is `alsolotar_toupie` (Task 6 Step 1), not a code
  reuse. ALS 2020 defers general HH *dimensions* to prior work and gives *bases* +
  Gerstenhaber structure, so the general toupie HH oracle is **cross-engine** (bar ≡ CS
  where bar survives). The a-Kronecker is char-independent (Euler `χ = 2 − a²`).
- **CS handles the EMPTY ideal (VERIFIED).** For the relation-free a-Kronecker,
  `engine="cs"` and `engine="cs", auto_cs=True` both return the right value
  (`[1,3,0,0,0]` for a=2), matching bar/auto — CS does NOT refuse on an empty ideal. So
  the bar≡CS battery is valid on relation-free toupies (they are just hereditary, hence
  HH^{≥2}=0 trivially — that is why M2 additionally requires a non-hereditary example).
- **Non-hereditary bar≡CS is checkable only where bar survives (M2, VERIFIED).** On the
  commutative `[2,2]` (non-monomial, dim 9) the bar coboundary `d^3` is `36864×4608 >
  max_cells` → `DepthLimitError`; bar succeeds to `top=2` (`[1,0,0]`), CS agrees
  (`[1,0,0]`) and continues past the window (`[1,0,0,0,0]`). So the M2 non-hereditary
  cross-engine check pins bar≡CS at `top=2` AND records the full CS value; do NOT ask
  bar for `top ≥ 3` on this input (it will raise, correctly).
- **`# BLOCK` Example 7.4.1.** Do NOT add the `[1,10,3,0,4,0,…]` pin — figure-only,
  unverified (see the reference section). Deferred item on the verification page.
- **`sl_a` is a LOWER BOUND, not equality, in general.** Thm 6.5 gives `HH¹ ⊇ sl_a`
  (`a = a_direct`); equality on the a-Kronecker (`HH¹ = sl_a`, `dim = a²−1`), strict on a
  mixed toupie with extra branches (verified: 2 direct + one length-2 branch → `HH¹ = 6 >
  3`). Test `dim HH¹ == a²−1` on the a-Kronecker and `dim HH¹ ≥ a²−1` on a mixed toupie
  with `a_direct ≥ 2`. Do NOT assert a `sl_a`-summand decomposition of `HH¹` (P70/R11).
- **`hochschild_cohomology` engine choice.** The a-Kronecker (`a ≥ 3`) is wild — fine for
  HH. Non-hereditary toupies blow the bar window fast; use `engine="cs", auto_cs=True`
  for the full profile and only ask bar where the window holds. Cap `top` at 4.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/families/toupie.py tests/families/test_toupie_hh_p59.py
git commit -m "feat(families): toupie HH oracle family -- a-Kronecker [1,a^2-1,0,..] pin (ALS 2020) + bar==CS on subdivided toupies + char-0 sl_a lower bound dim HH^1 >= a^2-1"
```

---

### Task 5: two algebra-only GUI kinds — `string_homological` + `toupie` (R34 + R35)

Both are expensive **algebra-only** scalar kinds (the `ar_quiver` precedent, NOT the
fast `recognizers` block), so no estimator edit beyond the algebra-dim sizing already in
place. The GUI touchpoints follow `ar_quiver`/`recognizers` verbatim.

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (two `_dispatch` branches + two `_snip` recipes),
  `docs/gui/runner.py` (the Pyodide twin — matching dispatch + `_snip` + ETA, kept
  shape-identical), `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkboxes,
  `S.ids`, push-list, `renderBlock`, `scheduleProbe`), `webapp/templates/index.html`
  (two checkboxes + the `ToupieAlgebra` preset button),
  `webapp/server/i18n/en.json` + `es.json` + `fr.json` + `zh.json` (**all four
  locales** — the two block key chains + the two `pick.kind.*` keys; exact key parity is
  gated by `tests/webapp/test_i18n.py::test_key_parity_with_english[es|fr|zh]`, so a
  missing fr/zh key lands the fast matrix RED),
  `src/quiverlab/trace/results_html.py` (`_HEADINGS` + two `_block_html` branches),
  `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (TWO goldens; existing ones byte-identical)
- Test: `tests/webapp/test_recognizer_kinds_p59.py`,
  `tests/gui/test_recognizer_runner_twin_p59.py`

**Interfaces (the `ar_quiver` kind is the exact template — algebra-only, schema v1,
shared block, per-runner `citations`):**
- `string_homological` → `string_homological_block(A)` (Task 1); block shape:
  ```python
  {"verdict": "not_string"|"string_over_ground_field"|"inconclusive", "is_string": bool,
   "witness": {...}|None, "kbar_gap_note": str|None, "reason": str|None,
   "checked": {"n_indec": int, "route": [...], "field": str, "classes": int},
   "references": ["suarez_alvarez", "butler_ringel", "assem_book"]}
  # a refusal (rep-infinite / self-injective / presentation-less) -> {"error": msg,
  #   "references": [...]}, never a 500 (the recognizers-block per-flag precedent).
  ```
- `toupie` → `toupie_block(A)`; block shape:
  ```python
  {"is_toupie": bool, "branch_count": int|None,          # out_deg(source) == in_deg(sink)
   "direct_arrow_count": int|None,                       # the sl_a `a` (length-1 branches)
   "hh": [int, ...]|None, "hh_top": int,                 # HH^0..HH^top dims (engine="cs",auto_cs)
   "sl_a": {"a": int, "dim": int, "char0": bool}|None,   # dim sl_a = a^2-1 (char 0); a=direct_arrow_count
   "note": str|None, "references": ["alsolotar_toupie", "assem_book"]}
  # is_toupie False -> {"is_toupie": False, ... nulls ..., "note": "not a toupie"}.
  # a_direct < 2 -> sl_a is {"a": a, "dim": max(a*a-1,0), ...} with a note (no content when a<2).
  # off char 0 -> sl_a: {"a": .., "dim": .., "char0": False} + a note (no inclusion claim).
  ```
- Both are **algebra-only scalar kinds** (schema v1, NO module block), routed by
  `spec.py::_dispatch` (NOT `_dispatch_module`); `estimator.sizing_dim` is algebra-dim
  based → **no estimator edit** (`ar_quiver`/`recognizers` precedent). `string_homological`
  is the pricier of the two (it knits + realizes extensions) — its `_snip`/ETA entry
  sizes like `ar_quiver` (a module budget), the `toupie` kind like a small HH.

- [ ] **Step 1: Write the failing cross-runner tests** (unmarked — extras-gated dir):

```python
# tests/webapp/test_recognizer_kinds_p59.py
"""The `string_homological` + `toupie` algebra-only kinds: served by hpc.spec, mirrored
by the Pyodide twin, byte-identical blocks."""


def test_string_homological_block_kD4(tmp_path):
    # schema-1 request: kD4 over QQ, compute ["string_homological"]. Assert:
    #   block["is_string"] is False; block["verdict"] == "not_string"
    #   block["witness"]["summand_count"] >= 3
    #   "suarez_alvarez" in [k for k, _ in block["citations"]]
    ...

def test_toupie_block_a_kronecker(tmp_path):
    # ToupieAlgebra([1,1,1]) over QQ, compute ["toupie"]. Assert:
    #   block["is_toupie"] is True; block["branch_count"] == 3; block["direct_arrow_count"] == 3
    #   block["hh"][:2] == [1, 8]; block["sl_a"]["dim"] == 8 and block["sl_a"]["char0"] is True
    #   "alsolotar_toupie" in [k for k, _ in block["citations"]]
    ...

def test_twin_parity(tmp_path):
    # run both requests through docs/gui/runner.py; json.dumps(sort_keys=True) equality
    # on each block (both runners byte-identical).
    ...
```

- [ ] **Step 2: Implement** both blocks' dispatch + the `docs/gui/runner.py` twin
  (shape-identical), the two `gui.js` checkboxes (`qlgui-string-homological`,
  `qlgui-toupie`) + `index.html` + the `ToupieAlgebra` **preset** button (a small
  `[1,1,1]` `a`-Kronecker seed drawn onto the canvas, mirroring the Brauer-star preset
  hook if present), the ETA entries, the `_snip` recipes
  (`"string_homological": "A.homological_string_test()"`,
  `"toupie": "A.is_toupie(), A.hochschild_cohomology(4)"` or the public spellings you
  expose), and `results_html.py` `_HEADINGS["string_homological"] = "Homological string
  test"` / `_HEADINGS["toupie"] = "Toupie structure"` + the two render branches (the
  witness rendered as a named `≥ 3` decomposition; the HH row as a degree table like the
  existing HH kinds; the `sl_a` line stating the char-0 scope + the `k̄`-gap note when
  present).

  **i18n keys — ALL FOUR LOCALES (`en`/`es`/`fr`/`zh`), H1.** The `ar_quiver` kind ships
  `pick.kind.ar_quiver` in every locale (line 246 of each file); mirror that. Add
  `pick.kind.string_homological`, `inv.string_homological`,
  `block.string_homological.{title,verdict,witness,arbiter}` and `pick.kind.toupie`,
  `inv.toupie`, `block.toupie.{title,recognizer,hh,sla}` to **all four** JSONs with
  exact key parity (else `tests/webapp/test_i18n.py::test_key_parity_with_english`
  fails). Suggested strings (the implementer may refine, but every key MUST be present
  in every locale):

  | key | en | es | fr | zh |
  |---|---|---|---|---|
  | `pick.kind.string_homological` | Homological string test | Prueba homológica de cuerdas | Test homologique de cordes | 同调弦代数判定 |
  | `inv.string_homological` | Homological string-algebra test | Prueba homológica de álgebra de cuerdas | Test homologique d'algèbre à cordes | 同调弦代数判定 |
  | `block.string_homological.title` | Homological string test (Suárez-Álvarez) | Prueba homológica de cuerdas (Suárez-Álvarez) | Test homologique de cordes (Suárez-Álvarez) | 同调弦判定 (Suárez-Álvarez) |
  | `block.string_homological.verdict` | Verdict | Veredicto | Verdict | 结论 |
  | `block.string_homological.witness` | Witness (≥3-summand middle) | Testigo (término medio con ≥3 sumandos) | Témoin (terme du milieu à ≥3 facteurs) | 见证 (中间项含≥3个直和项) |
  | `block.string_homological.arbiter` | Syntactic is_string | is_string sintáctico | is_string syntaxique | 句法 is_string |
  | `pick.kind.toupie` | Toupie structure | Estructura toupie | Structure toupie | 陀螺代数结构 |
  | `inv.toupie` | Toupie recognizer + HH | Reconocedor toupie + HH | Reconnaisseur toupie + HH | 陀螺代数识别 + HH |
  | `block.toupie.title` | Toupie structure | Estructura toupie | Structure toupie | 陀螺代数结构 |
  | `block.toupie.recognizer` | Is toupie | Es toupie | Est toupie | 是否陀螺代数 |
  | `block.toupie.hh` | Hochschild cohomology | Cohomología de Hochschild | Cohomologie de Hochschild | Hochschild 上同调 |
  | `block.toupie.sla` | sl_a ⊆ HH¹ (char 0) | sl_a ⊆ HH¹ (car. 0) | sl_a ⊆ HH¹ (car. 0) | sl_a ⊆ HH¹ (特征 0) |

- [ ] **Step 3: Add TWO goldens** (`string_homological_kD4`, `toupie_a_kronecker`) to
  `_runner_goldens.json`; note them in `test_runner_delegation.py`'s change-log
  docstring. Run the delegation test BEFORE adding to confirm existing goldens stay
  byte-identical. **Canonical-key stability:** both kinds are schema-v1 algebra-only, so
  they canonicalize through the Plan-25 `canonical_key` unchanged (no `module` block);
  confirm an existing family/quiver request's key is byte-identical after the wiring.

- [ ] **Step 4: Run the gates**
  `... -m pytest tests/webapp/test_recognizer_kinds_p59.py
  tests/webapp/test_runner_delegation.py tests/gui/test_recognizer_runner_twin_p59.py
  tests/hpc -q` — Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): string_homological + toupie algebra-only scalar kinds -- recognizer-panel rows, ToupieAlgebra preset, both runners byte-identical, i18n x4 (en/es/fr/zh) + pick.kind keys, two goldens"
```

---

### Task 6: verification page, citations, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick the P59 card)
- Test: existing release gates (`tests/release/test_oracle_classes.py`,
  `tests/citations/`)

- [ ] **Step 1: Citations** (BibTeX-VERIFIED only; the `_r(key, bibtex_key, kind, title,
  annotation, *tags)` registry precedent). Add to `references.bib`:

```bibtex
@article{SuarezAlvarez2023,
  author  = {Su\'arez-\'Alvarez, Mariano},
  title   = {A Simple Homological Characterization of String Algebras of Finite
             Representation Type},
  journal = {Algebras and Representation Theory},
  volume  = {26}, number = {5}, pages = {1759--1772}, year = {2023},
  doi     = {10.1007/s10468-022-10154-1}, note = {arXiv:2105.02948},
}
@article{ArtensteinLanzilottaSolotar2020,
  author  = {Artenstein, Dalia and Lanzilotta, Marcelo and Solotar, Andrea},
  title   = {Gerstenhaber structure on {H}ochschild cohomology of toupie algebras},
  journal = {Algebras and Representation Theory},
  volume  = {23}, number = {2}, pages = {421--456}, year = {2020},
  doi     = {10.1007/s10468-019-09854-y},
  note    = {arXiv:1803.10310; dimension bases in the Artenstein thesis, Colibri UdelaR},
}
```

  and in `registry.py` (mirror `_r("butler_ringel", "ButlerRingel1987", ...)`):

```python
_r("suarez_alvarez", "SuarezAlvarez2023", "algorithm",
   "A simple homological characterization of string algebras of finite rep. type",
   "Suarez-Alvarez: among rep-finite algebras, string <=> the middle term of EVERY "
   "extension of indecomposables has <= 2 summands (all Ext^1 classes, not just AR "
   "sequences) -- the homological string test.", "modules"),
_r("alsolotar_toupie", "ArtensteinLanzilottaSolotar2020", "family",
   "Hochschild cohomology of toupie algebras",
   "Artenstein-Lanzilotta-Solotar: toupie = unique source/sink + a parallel branches; "
   "a-Kronecker HH^* = [1, a^2-1, 0, ..]; HH^1 contains sl_a (char 0), a = #branches.",
   "families", "hochschild"),
```

  Verify BibTeX presence (`tests/citations/`) and that both keys resolve.

- [ ] **Step 2: Verification page.** Add the P59 subsystem rows:
  - `modules/string_homological.py` — `oracle_literature` (`kD₄`/`kD₅`/`kE₆` 3-middle
    witnesses; Suárez-Álvarez 2023); `oracle_crossengine` (the **discriminating battery**:
    a returned `not_string` always has `is_string False`, and the refute side RAISES on a
    contradiction); `oracle_selfcert` (every realized middle self-certifies —
    `baer_extension` cocycle + `ShortExactSequence` exactness; witness = genuine `≥ 3`
    decomposition).
  - `families/toupie.py` — `oracle_literature` (`a`-Kronecker `[1,a²−1,0,…]`, recomputed
    over QQ; `dim sl_a = a²−1` with `a` = # direct arrows, char 0, ALS 2020);
    `oracle_crossengine` (toupie HH bar ≡ CS where bar survives, incl. the non-hereditary
    commutative `[2,2]` at `top=2`; `is_toupie` ≡ hand quiver shape); `oracle_selfcert`
    (dimension certificate; connected+acyclic recognizer; `branch_count =
    out_deg(source) == in_deg(sink)`).
  - `qpa` — direct-session `IsSpecialBiserialAlgebra` parity (no `crosscheck.py` verb) on
    the R34 zoo; standing `NamesGVars()` guard that FAILS if QPA ever ships a
    homological-string / toupie surface.
  - **Honest-scope entries (binding):**
    (a) R34 is a **semi-decision**: a definitive `not_string` (a `≥ 3` witness,
    `k̄`-sound, field-independent) OR `string_over_ground_field` (a **finite ground field**
    exhaustive confirm that is **INCONCLUSIVE w.r.t. the Suárez-Álvarez `k̄` theorem** —
    over `k̄` there may be more indecomposables / non-rational classes) OR `inconclusive`.
    There is **no** bare `string` verdict; the definitive `string` answer is the syntactic
    `is_string`. The refute side RAISES iff a `≥ 3` witness meets `is_string True`
    (one engine is a bug); a `string_over_ground_field`-vs-`is_string` mismatch is a
    recorded **`k̄`-gap**, never raised. The candidate-realizability confirm route is a
    **scoped non-goal** (undecidable without algebraic-geometry machinery quiverlab lacks).
    (b) The refute PIPELINE (`knit`/`almost_split`/`decompose`) needs `char 0` or
    `char > dim` — only the terminal Fitting-split step is char-robust — so refute runs
    over **QQ** (small primes only on tiny algebras); the finite-field confirm runs over
    small primes with `char > dim` of the middles.
    (c) R35 pins **only** the `a`-Kronecker HH closed form from the literature; all other
    toupie HH values are **cross-engine** (ALS defers the dimension formulas to prior
    work). **Example 7.4.1 is deferred** (`# BLOCK`: quiver figure-only, unverified).
    (d) The `sl_a ⊆ HH¹` claim is **char 0 only** (ALS Thm 6.5, `k = ℂ`), with `a` = the
    number of **direct source→sink arrows** (NOT the branch count — verified: the
    commutative `[2,2]` has 2 branches but `HH¹ = 0`), delivered at the **dimension
    level** (`dim HH¹ ≥ a²−1`), not as a Lie-algebra decomposition (that is P70/R11).
  - **Recount the class table** (`tests/release/test_oracle_classes.py` drives the
    numbers — run collection, paste the LIVE counts, re-run to green; do NOT guess an
    at-authoring number given the mid-merge-train drift).

- [ ] **Step 3: README.** One features line: "the homological string-algebra test
  (Suárez-Álvarez: string ⇔ every extension of indecomposables has a ≤ 2-summand middle
  — a discriminating oracle against the syntactic recognizer), and toupie algebras
  (constructor + graph-shape recognizer + the `a`-Kronecker HH closed form + the char-0
  `sl_a ⊆ HH¹` inclusion) — the R34+R35 recognizer batteries."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/modules tests/families -q -m deep` (touched deep dirs),
  `... -m pytest tests/webapp tests/gui tests/hpc tests/invariants -q -m fast`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q`,
  and a citation-presence check (`suarez_alvarez`/`alsolotar_toupie` resolve; the
  `string_homological` block carries `suarez_alvarez`, the `toupie` block carries
  `alsolotar_toupie`) — all green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): P59 R34+R35 oracle rows + honest scope (three-valued semi-decision, k-bar caveat, a-Kronecker-only HH pin, char-0 sl_a) + suarez_alvarez/alsolotar_toupie citations + recounted classes"
```

---

## Acceptance (Plan-59 definition of done)

1. `homological_string_test(A) -> StringHomologyVerdict` public in
   `quiverlab.modules.string_homological`, rep-finite-only (loud refusal on
   rep-infinite / self-injective / presentation-less / non-complete-knit input); verdict
   in `{"not_string", "string_over_ground_field", "inconclusive"}` (no bare `"string"`);
   a self-certified witness (`baer_extension` cocycle + `ShortExactSequence` exactness,
   unified `_mk_witness` schema) on the `not_string` branch.
2. **The discriminating oracle holds:** the function RAISES iff a `≥ 3` witness meets
   `is_string(A) is True` (one engine is a bug); otherwise it RETURNS, and a returned
   `not_string` always has `is_string False`. `kD₄`/`kD₅`/`kE₆` return `not_string` with
   a `≥ 3`-summand witness (`kD₄` VERIFIED: `E = P₁⊕P₂⊕P₃`, dim vector `{0:3,1:1,2:1,3:1}`);
   `kAₙ`/gentle never return `not_string`. `string_over_ground_field`-vs-`is_string`
   mismatch is recorded as a `k̄`-gap, never raised.
3. **The quantifier is honest:** the refute side checks split + every AR sequence + a
   basis of each `Ext¹` (over QQ — the pipeline needs char 0 / char > dim; only the
   terminal split is char-robust); the confirm side is finite-ground-field exhaustive
   (orbit-reduced, budget-capped, `char > dim`) and yields `string_over_ground_field`
   (inconclusive w.r.t. the `k̄` theorem); infinite-field / char-undecidable / budget →
   `inconclusive`. The `k̄` scope and the candidate-realizability non-goal are on the
   verification page.
4. `ToupieAlgebra(branches, relations=None, field=None)` + `is_toupie(A)` +
   `toupie_branch_count(A)` + `toupie_direct_arrow_count(A)` public in
   `quiverlab.families` (and top-level `quiverlab`): the relation-free dimension
   certificate `dim = 2 + Σ(lᵢ−1) + Σ lᵢ(lᵢ+1)/2` holds (VERIFIED `[1,1]→4`,
   `[1,1,1]→5`, `[2]→6`, `[1,2]→7`, `[2,2]→10`; `a`-Kronecker `= ToupieAlgebra([1]*a)`,
   `dim = 2+a`); non-monomial cross-branch relations present (commutative `[2,2]` dim 9);
   the recognizer is **connected + acyclic** + degree-shape (rejects `kD₄`, and the
   path ⊔ oriented-cycle trap), and refuses a presentation-less algebra loudly; both
   discovery skip-sets (`hpc/spec.py`, `webapp/server/catalog.py`) carry `ToupieAlgebra`.
5. **Toupie HH oracle family:** the `a`-Kronecker `HH^• = [1, a²−1, 0, …]` recomputed over
   QQ (VERIFIED `a=2→[1,3,..]`, `a=3→[1,8,..]`; ALS 2020, char-independent) for
   `a ∈ {2,3,4}`; bar ≡ CS where bar survives, incl. the **non-hereditary** commutative
   `[2,2]` at `top=2` (`[1,0,0]==[1,0,0]`, full CS `[1,0,0,0,0]`); CS handles the empty
   ideal; Example 7.4.1 is **not** pinned (`# BLOCK`, figure-only).
6. **`sl_a` inclusion, char 0:** `toupie_sl_a_lower_bound(A) = a²−1` with
   `a = toupie_direct_arrow_count(A)` (# direct source→sink arrows, NOT branch count),
   hard-gated to char 0 (loud off it); `dim HH¹ = a²−1` on the `a`-Kronecker (equality)
   and `dim HH¹ ≥ a²−1` on a mixed toupie with `a ≥ 2` (VERIFIED: 2 direct + one length-2
   branch → `HH¹ = 6 ≥ 3`); ALS Thm 6.5.
7. Two **algebra-only** scalar kinds `string_homological` + `toupie` clickable
   end-to-end (GUI canvas → block → report) in **all four locales (en/es/fr/zh)** with
   `pick.kind.*` keys and exact key parity, schema v1, both runners
   byte-identical, TWO goldens added with a documented change-log entry; the witness /
   HH table / `sl_a` line render; the `ToupieAlgebra` preset seeds the canvas; canonical
   keys byte-stable (no `module` block).
8. Live QPA battery green (`-m qpa`): recognizer parity on the R34 zoo + the standing
   `NamesGVars()` guard that FAILS if QPA ever ships a homological-string / toupie
   surface (honest scope: QPA has neither today).
9. `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the four
   honest-scope entries (three-valued semi-decision + `k̄` gap + candidate-route non-goal;
   `decompose` char caveat; `a`-Kronecker-only HH pin + Example-7.4.1 deferral; char-0
   `sl_a`); `suarez_alvarez` + `alsolotar_toupie` citations added and BibTeX-verified;
   README line added; deep (touched dirs) + fast + qpa + release + citations suites
   green. No dependency taken on the other Wave-2 plans (merged independently to `dev`).

---

## Change log

- **2026-08-07 authoring.** Initial plan (R34 homological string test + R35 toupie
  algebras); references re-verified (ar5iv); toupie dimension formula hand-derived +
  arithmetic-checked.
- **2026-08-07 adversarial review: 1 blocking + 5 majors + 4 minors applied** (i18n ×4
  en/es/fr/zh + `pick.kind.*` keys; confirm demoted to `string_over_ground_field` with
  refute-only raising and a `k̄`-gap note; real QPA crosscheck = direct-session
  `IsSpecialBiserialAlgebra` parity, not a nonexistent `crosscheck` verb; both discovery
  skip-sets `hpc/spec.py`+`webapp/server/catalog.py` (NOT `discover.py`); non-hereditary
  relation-toupie in the HH battery, bar≡CS at `top=2`; `kD₄` witness live-verified
  `E = P₁⊕P₂⊕P₃`; char-scope honesty — refute pipeline is char-gated, only the terminal
  split is char-robust). **Plus a self-caught correction exposed by the M2 example:** the
  `sl_a` `a` is the number of **direct source→sink arrows**, NOT the branch count
  (commutative `[2,2]`: 2 branches, `HH¹ = 0`); `toupie_direct_arrow_count` added, the
  `[2,2,2]` sl_a test replaced by a live-verified mixed toupie (`HH¹ = 6 ≥ 3`); `is_toupie`
  gains connected+acyclic (the path ⊔ cycle trap); the `[2,2]→8` dimension typo fixed to
  `→10`; the "reuse P29 pin" claim replaced by an honest new-QQ-recomputation note.
