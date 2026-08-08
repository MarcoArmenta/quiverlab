# Plan 77: N-Koszul / K₂ / multi-Koszul recognizers (R36) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote the **Koszulity ladder beyond the quadratic case** to a first-class no-code
surface. For a finite-dimensional `A = kQ/I` over any exact Domain the plan ships four
generalized recognizers, all read off the **shipped Plan-05 minimal projective resolutions of
the simples** (= the graded Betti numbers of `Ext•_A(A/J, A/J)`, the Plan-27 Yoneda engine):

1. **Generation degrees of `Ext•(k,k)`** — the headline primitive: for every simple `S_i`
   the **internal (path-length) generation degrees** `ℓ_i(0), ℓ_i(1), ℓ_i(2), …` of its
   minimal resolution, recovered from the **homogeneous differentials** (each `d_n` entry is a
   path of a known length). Over a length-graded `A` these are exact; the extraction refuses
   loudly on a non-length-graded (inhomogeneous) presentation.
2. **N-Koszul (Berger)** — for an **N-homogeneous** `A` (every defining relation homogeneous of
   one length `N ≥ 2`), the **2-N alternation certificate**: the resolution of every simple is
   *pure*, with `ℓ_i(n)` equal to Berger's degree `δ(n) = (N/2)·n` (`n` even) / `(N/2)(n−1)+1`
   (`n` odd) through the certified window — the jumps alternate `1, N−1, 1, N−1, …`.
   Cross-checked against Berger's Ext-algebra characterization (for `N ≥ 3`, N-Koszul ⟺ the
   Yoneda algebra is generated in degrees `0, 1, 2`).
3. **K₂ (Cassidy–Shelton)** — the Yoneda algebra `E(A) = Ext•(k,k)` is **generated in
   cohomological degrees 1 and 2** as an algebra, decided from the shipped
   `generators_by_degree` through an **explicit certified window** `W`: complete when
   `gl.dim A` is finite (exact), **honestly inconclusive beyond `W`** otherwise; a genuine
   generator in degree `≥ 3` is a definitive **False**.
4. **(p,q)-almost-Koszul (Brenner–Butler–King)** — the classifier for the algebras P75's
   GHMS-Koszul route *refuses*: `A` is concentrated in degrees `0..p`, its trivial-module
   resolution is linear then breaks **once** with the error at internal degree `p + q`. The
   engine computes `p = ` top degree of `A` and `q = e − p` where `e` is the first
   internal-degree jump. This **classifies** the Dynkin preprojective algebras (`Π(A₃)`
   is (2,2)-, `Π(A₄)` is (3,2)-almost-Koszul, matching BBK's `(h−2, 2)`).

Plus a **multi-Koszul (Herscovich)** recognizer **scoped to connected-graded (local,
`A₀ = k`) input** — settled at spec time (§ "Multi-Koszul settlement"): Herscovich's notion is
a *connected-algebras* notion; the **transfer to multi-vertex `kQ/I` is via K₂**
(multi-Koszul `⟹` K₂ for finitely generated algebras with finite-dimensional relation space).
For multi-vertex input the recognizer **refuses loudly and points at K₂**.

One algebra-level compute kind — **`koszul`** — exposes the whole profile (the Plan-27
quadratic verdict + the four generalized fields + the generation-degree table) end-to-end in
all four locales. The metaplan's mandate ("koszul panel extension") is realized by growing the
existing Koszulity verdict block with the `n_koszul` / `k2` / `almost_koszul` /
`generation_degrees` fields.

**Quadratic overlap with Plan 27 (the scope boundary — read this first).** **Plan 27 OWNS
quadratic Koszulity.** It ships `modules/koszul.py` (`is_quadratic`, `g_quadratic_certificate`
= Priddy PBW, `quadratic_dual` = `A^!`, `froberg_obstruction`) and the three-valued
`YonedaPresentation.koszul` verdict in `modules/ext_algebra.py`. **Plan 77 OWNS the
generalized recognizers** (N-Koszul, K₂, almost-Koszul, multi-Koszul, generation degrees) and
does **not** re-implement or modify any quadratic surface. Concretely:

- The **N = 2 case defers to Plan 27**: `n_koszul_certificate` on a quadratic `A` reports
  `N = 2` and returns Plan-27's `koszul` verdict verbatim (no second opinion). `koszul.py`,
  `ext_algebra.py` are consumed **read-only**.
- Plan 77's new content is exactly the part Plan 27 does **not** have: the **internal
  (path-length) generation degrees**, the **N ≥ 3** alternation certificate, the **K₂**
  window-verdict (Plan 27 has no K₂ notion), the **almost-Koszul** classifier (Plan 27's
  `koszul` merely returns `False` with an obstruction — Plan 77 *classifies* the obstruction),
  and the scoped **multi-Koszul** recognizer.

**Architecture.** One new module `src/quiverlab/modules/nkoszul.py`, a thin recognizer layer
over the shipped Plan-05 minimal resolutions + the Plan-27 Yoneda engine (no new homological
engines):

- `nkoszul.py::generation_degrees(A, i, length)` — the internal generation degrees `ℓ_i(n)` of
  the minimal resolution of `S_i` (via the homogeneous differentials). Loud `QuiverlabError`
  on a non-length-graded presentation (`koszul._is_length_graded` gate) and on an inhomogeneous
  differential column (`_NonPure`).
- `nkoszul.py::n_homogeneous_degree(A)` — the single relation degree `N ≥ 2` when every
  defining relation is homogeneous of the same length, else `None`.
- `nkoszul.py::n_koszul_certificate(A, top)` — the three-valued N-Koszul verdict + `N` +
  window; the Berger-pattern check on every simple, cross-checked against K₂ + N-homogeneity.
- `nkoszul.py::k2_certificate(A, top)` — the three-valued K₂ verdict + the explicit window.
- `nkoszul.py::almost_koszul_certificate(A, top)` — the `(p, q)` classifier + the BBK data.
- `nkoszul.py::multi_koszul_certificate(A, top)` — connected-graded-only; refuses multi-vertex
  loudly, reporting K₂ as the transferable property.
- `nkoszul.py::koszul_profile(A, top)` and `koszul_profile_block(A, top)` — the assembled
  record + the no-code block (shared by both runners, byte-identical).
- `Algebra.koszul_profile(top=8)` — thin lazy-import delegate (`core/algebra.py`), beside
  `Algebra.ext_algebra`.

**Plan 27 is left byte-unchanged.** `koszul.is_quadratic`/`g_quadratic_certificate`/
`_is_length_graded`, `ext_algebra.ext_algebra`/`YonedaPresentation` are consumed read-only. No
Plan-27 test moves.

**Tech Stack.** Pure exact arithmetic over `A.domain` (the Plan-05/27 precedent — the module
engine runs over any exact Domain). Internal degrees are **integer path lengths** (field-free
combinatorics, exactly like `koszul._graded_count_matrices`). No floats in `src/` (AST-gated by
`tests/test_no_floats.py`). The only linear algebra is inherited (the Plan-05 covers, the
Plan-27 Yoneda lifts).

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`. The
  recognizer path is pure-Python (module engine); `QUIVERLAB_NO_NUMBA=1` changes nothing.
- **Plan 27 (Yoneda / Ext-algebra + Koszulity) is a HARD prerequisite, MERGED to `main` since
  2026-07-25** (byte-identical on `dev` — verified: `git diff v0.3.0 dev` touches neither
  `modules/ext_algebra.py` nor `modules/koszul.py` nor `modules/resolution.py`). This plan
  consumes, read-only: `modules.ext_algebra.ext_algebra`/`YonedaPresentation`
  (`generators_by_degree`, `graded_dims_through`, `hilbert_matrix_through`,
  `certified_through_degree`, `koszul`/`koszul_obstruction`, `is_finite_dimensional`),
  `modules.ext_algebra.ext_algebra_block`, `modules.koszul.is_quadratic`/
  `g_quadratic_certificate`/`_is_length_graded`, `modules.resolution.minimal_resolution`,
  `modules.builders.projective`. **If Plan 27's surface has shifted when a worker picks this
  up, STOP and re-grep before forking anything** — do not re-implement the Yoneda engine.
- **Wave-4 independence.** P77 is `[independent]` in the metaplan wave map. It does **not**
  depend on P75 (incidence + GHMS-Koszul); the P75 **seam is conceptual** (both read the same
  `ext_algebra.koszul` obstruction). P75 may or may not be merged when P77 is picked up — P77
  neither imports P75 code nor moves any P75 test. Wave-1/2/3 merges (P51–P69) are irrelevant
  to P77 (no shared surface).
- **Length-graded gate.** The internal-degree extraction and the N-Koszul / almost-Koszul
  recognizers require `A` **length-graded** (`I` homogeneous — every relation length-pure), so
  the minimal resolution is graded and the differential path-lengths are well-defined. The gate
  is the shipped `koszul._is_length_graded(A)`. On a non-length-graded `A` the recognizers
  return `None` with an honest reason (`"not length-graded: internal degrees undefined"`) — a
  clean refusal, never a wrong pattern. **K₂ needs only the shipped `generators_by_degree`
  (homological), so it is defined for any length-graded `A` too**; an inhomogeneous ideal makes
  the whole graded story vacuous and every recognizer reports `None`.
- **The certified window is contractual (metaplan §3 honest-scope).** Every verdict that reads
  the Yoneda engine is bounded by `YonedaPresentation.certified_through_degree`:
  - `gl.dim A` **finite (exact)** ⇒ `E(A)` is finite-dimensional, the window is **complete**,
    the verdicts are decisive (`complete=True`).
  - `gl.dim A` **not exact-finite** (self-injective, radical-square-zero on several
    loops, …) ⇒ verdicts are certified **through the window `W = certified_through_degree`**
    only, with `complete=False` and an honest note; a `True` becomes "holds through `W`", a
    definitive `False` (a generator/pattern-break **inside** `W`) stays `False`. **Never claim
    a window-bounded `True` as unconditional.**
- **Composition is left-to-right** (ASS); all module maps inherit the Plan-05/27 conventions.
- All refusals are `QuiverlabError(message, hint=...)` (the two-space `[hint: ...]` convention,
  `errors.py`). Presentation-less (structure-constant) algebras carry no path basis; the Yoneda
  engine + `koszul` routines already refuse loudly — Plan 77 surfaces that as a clean typed 4xx,
  never a 500.
- **Every citation key is BibTeX-verified before use** (Plan-29 rule); `_r(...)` registry entry
  + `references.bib` entry both land, gated by `tests/citations/test_bib_structure.py`. **BBK
  is NOT on arXiv** — the `references.bib` entry is a plain `@article` (Algebras and
  Representation Theory), no `eprint`.
- **`koszul` is an ALGEBRA-LEVEL compute kind — it adds NO request-schema field** (the
  `ext_algebra`/`derived_fingerprint` precedent: a new kind is a new *compute-list string*, not
  a new request key; the request schema — which the live tree already speaks at v1/v2/v3 — is
  untouched). **Canonical keys are request-derived** (`cache.py::canonical_key`): a request that
  does not use the new `koszul` compute string keys byte-identically to before, so existing
  goldens are untouched. Verify with `test_runner_delegation.py` BEFORE and AFTER adding the
  golden.
- **Plan-32 markers:** the Berger-pattern pins (`k[x]/x^N` internal degrees), the K₂
  Cassidy–Shelton definition on `k[x]/x^N`, the BBK preprojective classification
  (`Π(A₃)`=(2,2), `Π(A₄)`=(3,2)) = `oracle_literature`; the **two-certificate agreement**
  (internal-degree purity ≡ K₂ + N-homogeneous), `generation_degrees ≡ Berger`, and the
  P77-Yoneda-gens ≡ Plan-27 `generators_by_degree` = `oracle_crossengine`; the extraction
  self-consistency (homogeneity), the Berger closed form, `q = e − p`, and the **P75-seam
  agreement** (the `ext_algebra.koszul` obstruction degree = the almost-Koszul break degree) =
  `oracle_selfcert`; **QPA has NO Koszul surface of any kind** (Plan-27 finding — "QPA has no
  IsKoszul") — the `qpa` bucket is an honest-scope probe that FAILS if QPA ever ships one.
- **Mid-merge-train counts drift.** Task 8 recounts the oracle-class table at merge time by
  running `tests/release/test_oracle_classes.py` (paste-the-live-numbers, never a guess) and
  claims only this plan's deltas.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table
  green) and the README line. Conventional commits; green at every commit; branch
  `plan-77-nkoszul` off `dev`.

---

## Record (verbatim — the mandate for this plan)

> **R36 — N-Koszul / K₂ / multi-Koszul certifiers.** [A-scout P7; keep-with-corrections]
> Object: generation degrees of Ext•(k,k) off the shipped minimal resolutions; N-Koszul
> (single-degree relations + the 2-N alternation pattern — checkable), K₂ (Cassidy–Shelton,
> generated in degrees 1,2 — needs an explicit certified WINDOW, honest inconclusive beyond),
> multi-Koszul (Herscovich — stated for connected graded; the f.d. kQ/I transfer needs care at
> spec time). Quadratic case = existing Plan-27 (overlap named). Refs: Herscovich
> arXiv:1305.1678; Herscovich JPAA 223 (2019) 1054–1072 (year corrected); Chouhy
> arXiv:1708.02933 (N-Koszul degeneration stability). Size M.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P77.

**Card additions binding on this plan (metaplan §5 P77 + the writer's dispatch):**
- N-Koszul recognition = the **2-N alternation check** (pd pattern of the linear/alternating
  resolution degrees per Berger).
- K₂ with an **explicit certified window** (honest inconclusive beyond).
- **multi-Koszul with the f.d. `kQ/I` transfer SETTLED at spec time** (done — § below).
- Extra ref binding: **Brenner–Butler–King, "Periodic algebras which are almost Koszul"** —
  the (p,q)-almost-Koszul classifier for the P75-refused preprojectives.
- The quadratic overlap with Plan 27 must be **NAMED** (done — the scope boundary above).

---

## Multi-Koszul settlement (spec-time decision — the card's mandate)

**The question (card):** Herscovich's multi-Koszul is *"stated for connected graded"*; does it
transfer to finite-dimensional `kQ/I`?

**What the sources say (re-read at authoring — see Reference re-verification):**

- Herscovich, *On the multi-Koszul property for connected algebras* (arXiv:1305.1678, 2013)
  defines multi-Koszul for a **"locally finite dimensional nonnegatively graded connected
  algebra"** — i.e. `A₀ = k` (connected / **local**), `A = ⊕_{n≥0} A_n`, generated in
  degree 1, `dim A_n < ∞`. The definition generalizes Berger's N-Koszul to algebras whose
  relations live in **several** degrees (not one degree `N`).
- Herscovich, *Applications of one-point extensions … A∞-(co)module structure …*, **J. Pure
  Appl. Algebra 223 (2019), no. 3, 1054–1072** (the card's "JPAA 223") gives the A∞-structure
  of the Yoneda algebra of a multi-Koszul algebra — same connected-graded setting.
- The **transfer theorem** (Herscovich 1305.1678): *a finitely generated multi-Koszul algebra
  with a finite-dimensional space of relations is a **K₂ algebra** in the sense of
  Cassidy–Shelton.* So **multi-Koszul `⟹` K₂**.

**The settlement (decided here, not deferred):**

1. **Multi-Koszul is a CONNECTED-graded (`A₀ = k`) notion.** In quiverlab's world this is
   exactly the **LOCAL** finite-dimensional algebras: a single-vertex `kQ/I` (one vertex, any
   number of loops) with a length-homogeneous ideal — `A₀ = k`. `k[x]/(x^N)`, `k⟨x,y⟩/I`
   (local) qualify. **Multi-vertex `kQ/I` has `A₀ = k^{Q₀}` semisimple, NOT connected — outside
   Herscovich's stated scope.**
2. **The transfer to multi-vertex `kQ/I` is K₂.** Since multi-Koszul `⟹` K₂ (for f.g.
   finite-relation-space algebras) and K₂ is defined over any semisimple `A₀`
   (Cassidy–Shelton work graded over a field, but the Yoneda-generation criterion is
   basis-independent and quiverlab's Yoneda engine is already multi-vertex), **K₂ is the
   property quiverlab offers for every `kQ/I`**; multi-Koszul is offered **only for connected
   (local) input**, refusing multi-vertex loudly with a pointer to K₂.
3. **Honest scope (binding on the verification page):** Herscovich's headline multi-Koszul
   examples — the **Yang–Mills and super-Yang–Mills algebras** — are **infinite-dimensional**
   connected graded algebras, *outside* quiverlab's finite-dimensional engine. So the
   f.d.-local multi-Koszul examples are genuinely sparse: `k[x]/(x^N)` (a single-relation-degree
   multi-Koszul = N-Koszul) and constructed two-relation-degree local algebras. This is stated,
   not hidden.
4. **The transfer is PINNED NOW as Herscovich Prop. 3.30; the multi-Koszul DECISION is scoped
   OUT (not merely un-transcribed).** Re-read from the **ar5iv HTML**
   (`https://ar5iv.labs.arxiv.org/abs/1305.1678`) at the fix round: **Proposition 3.30** —
   *"the Yoneda algebra of a finitely generated multi-Koszul algebra with a finite dimensional
   space of relations is generated in degrees 1 and 2, so a `𝒦₂` algebra"* (multi-Koszul ⟹ K₂,
   quoted). **Multi-Koszul itself is NOT a single `Definition` environment** — it is built up in
   **§3.2 ("The definition of multi-Koszul algebras")** as a construction requiring specific
   vanishing/finiteness of the **Tor / Ext groups of the minimal graded projective resolution of
   `A` as a bimodule**. **That bimodule-Tor decision machinery is genuinely OUT of this plan's
   scope** (it is a homological engine of its own, not a recognizer over the shipped surfaces) —
   this is the honest deferral, and it is stated as scope, not hidden as a missing transcription.
5. **What v1 ships for multi-Koszul (binding):** on connected (local, single-vertex, `A₀ = k`)
   length-graded input, the recognizer reports the **generation degrees of `Ext•(k,k)`** (the
   headline primitive) + the **K₂ verdict** (which multi-Koszul implies, Prop. 3.30) + an honest
   `status` naming the deferred decision (`"the multi-Koszul DECISION requires bimodule Tor/Ext
   vanishing per Herscovich §3.2 -- out of this plan's scope; K2 (Prop 3.30) and the generation
   degrees are reported"`). **`verdict` is `None` (never a fabricated `True`/`False`).** The
   full multi-Koszul decision is a named entry in `docs/plans/DEEPER-ENGINES-BACKLOG.md`
   (metaplan §1.2 GUI-deferral-ledger discipline). On **multi-vertex** input the recognizer
   refuses loudly and points at K₂ (item 2).

This satisfies the card both ways: the transfer **is** settled (multi-Koszul is connected-only;
the multi-vertex bridge is K₂, Prop. 3.30 pinned), and the one piece not shipped — the
multi-Koszul *decision* — is an honestly-scoped deferral (the §3.2 bimodule-Tor machinery), not
a silent gap or a guessed verdict.

---

## Reference re-verification (mandatory; findings recorded)

All statements below were re-read from the primary sources during authoring (WebFetch /
WebSearch on the arXiv abstracts, the Bath research-portal record for BBK, and the published
venues). The **empirical ground truth** (every pin) was independently produced by running the
shipped `ext_algebra` + a prototype internal-degree extractor on the `.venv` (v0.3.0 tree,
whose `ext_algebra.py`/`koszul.py`/`resolution.py` are byte-identical to `dev` — verified) and
is transcribed in the Live-verified pins section.

- **Berger, *Koszulity for nonquadratic algebras*, J. Algebra 239 (2001), no. 2, 705–734 —
  VERIFIED.** N-Koszul generalizes Koszulity from quadratic to `N`-homogeneous algebras
  (relations homogeneous of one degree `N > 2`) via a *pure* minimal resolution of the trivial
  module: `P_n` generated in the single internal degree `δ(n) = (N/2)·n` (`n` even),
  `(N/2)(n−1) + 1` (`n` odd). **Ext-algebra characterization (`N ≥ 3`):** the algebra is
  N-Koszul **iff the Yoneda algebra `Ext•(k,k)` is generated in degrees `0, 1, 2`** (WebSearch,
  Berger's classification; also Green–Marcos–Martínez-Villa–Zhang for the multi-vertex/
  semisimple-base version). This is the double certificate Plan 77 cross-checks: **pure
  resolution (internal degrees = `δ(n)`) ≡ (`N`-homogeneous ∧ K₂)**.
- **Cassidy, Shelton, *Generalizing the notion of Koszul algebra*, Math. Z. 260 (2008), no. 1,
  93–114 (arXiv:0704.3752) — VERIFIED VERBATIM (definition).** *"A graded `k`-algebra `A` is
  **K₂** if its Yoneda algebra `Ext•_A(k,k)` is generated as an algebra in cohomological
  degrees 1 and 2."* K₂ generalizes both Koszul and N-Koszul: **Koszul `⟹` K₂**,
  **N-Koszul `⟹` K₂**, and K₂ *"may have its ideal of relations generated in different
  degrees"* — strictly more general. This is the exact criterion `k2_certificate` decides
  (through the window).
- **Brenner, Butler, King, *Periodic algebras which are almost Koszul*, Algebras and
  Representation Theory 5 (2002), no. 4, 331–368 — VERIFIED (definition + preprojective
  classification), Bath research portal + Springer.** *"An algebra `A` is almost Koszul, or
  **(p,q)-Koszul**, if `A` is concentrated in degrees 0 to `p` and there is a linear complex of
  projective modules which resolves `S` up to an error given by the degree `p+q` part of the
  resolution."* **The classical preprojective algebras of simply-laced Dynkin type are
  `(h−2, 2)`-Koszul, where `h` is the Coxeter number** (type `Aₛ`: `h = s+1`), and the algebra
  is periodic of period `2(h−1)`. **CRITICAL READING (resolved at authoring — see Methodology):
  `p` is the algebra's TOP DEGREE (concentration), NOT the homological length of the linear
  strand of the trivial-module minimal resolution.** The engine's minimal resolution of the
  simple has its **first** non-linear jump at a *fixed* homological degree (empirically 3 for
  every Dynkin preprojective tried), to internal degree `e = p + q`; recovering `q = e − p` with
  `p = ` top degree reproduces BBK's `(h−2, 2)` on the nose for `Π(A₃)`, `Π(A₄)`, `Π(A₅)`,
  `Π(D₄)`. **BBK is not on arXiv** — plain `@article`, no `eprint`. This is the metaplan's
  named extra reference.
- **Herscovich, *On the multi-Koszul property for connected algebras*, arXiv:1305.1678 (2013) —
  VERIFIED (scope + transfer + structure), via the ar5iv HTML.** Hypotheses (§2): *"locally
  finite dimensional nonnegatively graded connected algebra"*, `A = ⊕_{n≥0} A_n` with
  `A₀ = k` — connected/**local**, `dim A_n < ∞`. **Transfer — Proposition 3.30 (quoted
  verbatim):** *"the Yoneda algebra of a finitely generated multi-Koszul algebra with a finite
  dimensional space of relations is generated in degrees 1 and 2, so a `𝒦₂` algebra."* So
  **multi-Koszul ⟹ K₂**, PINNED (no `# PIN` left on the transfer). **Multi-Koszul is NOT a
  single `Definition` environment** — it is constructed in **§3.2** via Tor/Ext-vanishing of the
  minimal graded bimodule resolution; **that decision machinery is scoped OUT of this plan**
  (bimodule-Tor engine — Multi-Koszul settlement items 4–5). The paper does **not** treat
  semisimple `A₀` (several vertices), confirming the connected-only scope.
- **Herscovich, J. Pure Appl. Algebra 223 (2019), no. 3, 1054–1072 — IDENTIFIED (title
  corrected).** Exact title *"Applications of one-point extensions to compute the
  A∞-(co)module structure of several Ext (resp., Tor) groups"* — the A∞-structure of the
  Yoneda algebra of a multi-Koszul algebra (the card's "JPAA 223 … multi-Koszul", year
  corrected to 2019). Cite as **context** (the A∞ refinement), not a computed oracle.
- **Chouhy, *On geometric degenerations and Gerstenhaber formal deformations*, arXiv:1708.02933,
  Bull. London Math. Soc. 51 (2019) — VERIFIED.** *"For finite dimensional associative
  algebras, the N-Koszul property is preserved under the degeneration relation for all
  `N ≥ 2`."* This is the card's "N-Koszul degeneration stability". Cited as **context / a
  theorem oracle** (computing degenerations is out of scope; the stability statement frames the
  robustness of the N-Koszul verdict — no degeneration engine is built).
- **Green, Marcos, Martínez-Villa, Zhang, *D-Koszul algebras*, J. Pure Appl. Algebra 193
  (2004), no. 1–3, 141–162 — CITED (multi-vertex/semisimple-base foundation).** The N-Koszul /
  δ-Koszul theory **over a semisimple base `A₀ = k^{Q₀}`** (several vertices) — the foundation
  that legitimizes Plan 77's N-Koszul recognizer on multi-vertex `kQ/I`. **`# verify` exact
  volume/pages** at implementation before committing the key.

---

## Live-verified pins (produced on the `.venv` at authoring; transcribe into tests as-is)

Every table below is a LIVE run of the shipped `modules.ext_algebra.ext_algebra_block` and a
prototype of `generation_degrees` (the internal-degree extractor). The prototype was **validated
against Berger's closed form on four `k[x]/(x^N)`** before any preprojective pin was trusted.

**Constructors:** `truncated_polynomial(N)` = `k[x]/(x^N)`; `TruncatedPathAlgebra("A4", 3)` =
`kA₄/J³`; `PreprojectiveAlgebra("A3")` etc. All over `QQ`.

### A. `k[x]/(x^N)` — the N-Koszul spine (Berger)

| `A` | `koszul` (P27) | obstruction | `gens_by_deg` (homological) | internal degrees `ℓ(n)` | Berger `δ(n)` |
|---|---|---|---|---|---|
| `k[x]/(x²)` | **True** (G-quadratic) | — | `{1:1}` | `0,1,2,3,4,5,6,7,8` | `0,1,2,3,…` ✓ |
| `k[x]/(x³)` | **False** | `(2,'non-quadratic: a defining relation is not length 2')` | `{1:1, 2:1}` | `0,1,3,4,6,7,9,10,12` | ✓ (jumps 1,2,1,2…) |
| `k[x]/(x⁴)` | **False** | same | `{1:1, 2:1}` | `0,1,4,5,8,9,12,13,16` | ✓ (jumps 1,3,1,3…) |
| `k[x]/(x⁵)` | **False** | same | `{1:1, 2:1}` | `0,1,5,6,10,11,15,16,20` | ✓ (jumps 1,4,1,4…) |

**Reading:** `gens_by_deg` is IDENTICAL `{1:1, 2:1}` for `N = 3, 4, 5` — **`N` is NOT
recoverable from homological degrees; the internal degrees are essential** (this is why Plan 77
needs the extraction primitive, not just Plan 27's Yoneda gens). Every `N ≥ 3` case is **K₂**
(gens ⊆ {1,2}) AND N-Koszul (internal degrees = `δ(n)`), the two certificates agreeing.
`E(A)`'s graded dims are all `1` (`Ext^n(k,k) = k`).

### B. `kA₄/J³` — a finite-gl.dim cubic monomial (multi-vertex N-Koszul)

`TruncatedPathAlgebra("A4", 3)` = `kA₄` mod the single length-3 path `e12·e23·e34`. dim 9,
`gl.dim = 2` (finite): `koszul=False`, obstruction `(2,'non-quadratic …')`, `gens_by_deg =
{1:3, 2:1}` (K₂ ✓), `E` graded dims `[4,3,1]`. Internal degrees `S₁: [0,1,3]` = Berger `δ`
for `N=3` (the resolution is finite, pd `S₁ = 2`). **⟹ 3-Koszul, multi-vertex, both
certificates agree.**

### C. Preprojective `Π(Δ)` — the almost-Koszul boundary (BBK), the P75 seam

| `A` | dim | `koszul` | obstruction | `gens_by_deg` | internal degrees `ℓ(n)` | top deg `p` | first jump `e` | `q=e−p` | BBK `(h−2,2)` |
|---|---|---|---|---|---|---|---|---|---|
| `Π(A₂)` | 4 | **True** | — | `{1:2}` | `0,1,2,3,4,5,6,7,8` | 1 | — (Koszul) | — | Koszul boundary (`rad²=0`) |
| `Π(A₃)` | 10 | **False** | `(3,'a new Ext-algebra generator appears in degree 3')` | `{1:4, 3:3}` | `0,1,2,4,5,6,8` | 2 | 4 | **2** | **(2,2)**, `h=4` ✓ |
| `Π(A₄)` | 20 | **False** | `(3,…degree 3)` | `{1:6, 3:4}` | `0,1,2,5,6,7,10` | 3 | 5 | **2** | **(3,2)**, `h=5` ✓ |
| `Π(A₅)` | 35 | **False** | `(3,…degree 3)` | — | `0,1,2,6,7,8` | 4 | 6 | **2** | **(4,2)**, `h=6` ✓ |
| `Π(D₄)` | 28 | **False** | `(3,…degree 3)` | — | `0,1,2,6,7,8` | 4 | 6 | **2** | **(4,2)**, `h=6` ✓ |

**Readings:**
- `Π(A₂)` is genuinely **Koszul** (`rad²=0`; linear resolution forever): the boundary — a
  preprojective that IS Koszul. Not almost-Koszul (no jump).
- `Π(A_{n≥3})` is **(p,q)-almost-Koszul** with `p = ` top degree `= h−2`, `q = e − p = 2` — the
  engine reproduces BBK's `(h−2, 2)` exactly via `q = e − p`. **All break at homological degree
  3** with the shipped obstruction `(3, 'a new Ext-algebra generator appears in degree 3')`.
- `Π(A_{n≥3})` is **NOT K₂** (a genuine Yoneda generator at cohomological degree 3 —
  `gens_by_deg` has key `3`): a clean discriminating K₂-False family.
- **The P75 seam (the card's mandate):** P75's GHMS-Koszul route REFUSES `Π(A₃)` with exactly
  the `ext_algebra.koszul` obstruction `(3, …)`; **Plan 77's `almost_koszul_certificate`
  CLASSIFIES the same input as `(2,2)`** — the seam is that the obstruction *degree* equals the
  almost-Koszul *break degree* (a `oracle_selfcert` pin).

---

### Task 0: citations (do first — later tasks reference the keys)

**Files:**
- Modify: `src/quiverlab/citations/references.bib`, `src/quiverlab/citations/registry.py`
- Test: `tests/citations/test_bib_structure.py` (existing gate)

**Interfaces:** the `_r(key, bibtex_key, kind, title, annotation, *tags)` registry helper.
`priddy`, `froberg_koszul`, `polishchuk_positselski` **already exist** (Plan 27) — reuse for
the quadratic overlap, do NOT re-add. `_citation_pairs` (`hpc/spec.py`) resolves keys →
`[key, formatted]` pairs for payloads.

- [ ] **Step 1: add the keys** (BibTeX-verify each before committing — Plan-29 rule; **BBK has
  no arXiv eprint**):
  - `berger_nonquadratic` → Berger, *Koszulity for nonquadratic algebras*, J. Algebra **239**
    (2001), no. 2, 705–734, doi 10.1006/jabr.2000.8703. (N-Koszul; pure resolution `δ(n)`;
    `N≥3` ⟺ Ext-algebra generated in degrees 0,1,2.)
  - `cassidy_shelton` → Cassidy, Shelton, *Generalizing the notion of Koszul algebra*, Math. Z.
    **260** (2008), no. 1, 93–114, arXiv:0704.3752, doi 10.1007/s00209-007-0266-5. (K₂ = Yoneda
    generated in degrees 1,2.)
  - `brenner_butler_king` → Brenner, Butler, King, *Periodic algebras which are almost Koszul*,
    Algebr. Represent. Theory **5** (2002), no. 4, 331–368, doi 10.1023/A:1020146502185. (**No
    eprint.** (p,q)-almost-Koszul; preprojective `(h−2,2)`; period `2(h−1)`.)
  - `herscovich_multikoszul` → Herscovich, *On the multi-Koszul property for connected
    algebras*, arXiv:1305.1678 (2013). (multi-Koszul for connected graded; multi-Koszul ⟹ K₂.)
  - `herscovich_ainfty_ext` → Herscovich, *Applications of one-point extensions to compute the
    A∞-(co)module structure of several Ext (resp., Tor) groups*, J. Pure Appl. Algebra **223**
    (2019), no. 3, 1054–1072. (the card's "JPAA 223"; A∞ Yoneda of a multi-Koszul algebra —
    context.) **Title + volume/no./pages CONFIRMED** (WebSearch, fix round); `# verify` the doi
    (`10.1016/j.jpaa.2018.05.015` plausible, not independently resolved).
  - `chouhy_degenerations` → Chouhy, *On geometric degenerations and Gerstenhaber formal
    deformations*, Bull. Lond. Math. Soc. **51** (2019), arXiv:1708.02933,
    doi **10.1112/blms.12277**. (N-Koszul preserved under degeneration for all `N≥2` — context.)
    **doi CONFIRMED** (Wiley record, fix round); `# verify` the exact issue/pages
    (no. 5, 890–908 — from the arXiv listing, not the Wiley front matter).
  - `green_marcos_martinezvilla_zhang` → Green, Marcos, Martínez-Villa, Zhang, *D-Koszul
    algebras*, J. Pure Appl. Algebra **193** (2004), no. 1–3, 141–162. (N-Koszul over a
    semisimple base — the multi-vertex foundation.) `# verify` volume/pages.

  Registry entries mirror the Plan-27 Koszul keys' shape. **Annotations must name the corrected
  attributions:** Cassidy–Shelton = K₂ (degrees 1,2); BBK = almost-Koszul `(p,q)`, `p`= top
  degree; Herscovich 1305.1678 = connected-graded multi-Koszul, transfer = K₂; Chouhy/JPAA =
  context only.
- [ ] **Step 2:** `... -m pytest tests/citations/test_bib_structure.py -q` green.
- [ ] **Step 3: Commit** — `docs(citations): N-Koszul/K2/almost-Koszul/multi-Koszul refs (Berger, Cassidy-Shelton, Brenner-Butler-King, Herscovich, Chouhy, GMMVZ)`

---

### Task 1: `nkoszul.py` — generation degrees of `Ext•(k,k)` (THE headline primitive)

**Files:**
- Create: `src/quiverlab/modules/nkoszul.py`
- Test: `tests/modules/test_nkoszul_generation_degrees.py`

**Interfaces:**
- Consumes: `modules.resolution.minimal_resolution(M, length)` (Plan-05 terms + differentials),
  `modules.builders.projective(A, v)` (the `_pv_basis_labels` path labels),
  `modules.koszul._is_length_graded(A)` (the homogeneity gate).
- Produces:
  ```python
  def n_homogeneous_degree(A) -> int | None:
      # The single relation length N>=2 iff every defining relation is homogeneous of the SAME
      # length (rel.min_length == rel.max_length == N for all).
      #   * HEREDITARY kQ (relations == []) -> returns 2 (quadratic-TRIVIALLY, matching
      #     koszul.is_quadratic's vacuous-True convention: kQ is the quadratic algebra with an
      #     empty relation space). So callers treat hereditary as the N=2 (defer-to-Plan-27) case
      #     -- no special sentinel branch. This is what makes the kA2 test pin n_homogeneous==2.
      #   * None otherwise (mixed / inhomogeneous relation degrees).

  def generation_degrees(A, i, length, *, max_term_dim=200000) -> list[list[int]]:
      # per homological degree n = 0..: the SORTED SET of internal (path-length) generation
      # degrees of the summand generators of P_n in the minimal resolution of S_i. Requires
      # koszul._is_length_graded(A) (else QuiverlabError). Raises _NonPure(n) if a generator
      # column mixes internal degrees (the resolution is not graded-pure -> caller reports None).
  ```

**The extraction (validated against Berger — Methodology).** For a length-graded `A` the
minimal resolution of `S_i` is graded; each `d_n : P_n -> P_{n-1}` is homogeneous, so the
internal degree of a `P_n` summand generator is well-defined:
- `P_0 = P_i` covers `S_i` in internal degree 0 ⇒ every `P_0` generator has degree 0.
- For `n ≥ 1`, summand `s` of `P_n` has generator column `off_s`. For each nonzero entry
  `d_n[r][off_s]`, locate the summand `s'` of `P_{n-1}` containing row `r` and the **path length
  `ℓ_path`** of that basis vector; the target internal degree is `ℓ_{n-1}(s') + ℓ_path`. **All
  nonzero entries of the column must agree** (homogeneity). That common value is `ℓ_n(s)`.
- **`_NonPure` is a DEFENSIVE assert, not the validation (`H3`).** On a length-graded `A` the
  resolution is graded, so the column-agreement holds by construction — `_NonPure` fires only on
  a bug or a mis-gated non-graded input. **The REAL correctness anchor is Berger's CLOSED FORM**
  (`δ(n)`), which the extractor reproduces exactly on `k[x]/(x^N)` (`N=2..5`) — an external
  ground truth, not a self-consistency check. **Any independent re-derivation of the internal
  degrees uses the SAME method** (internal degree of a `P_n` generator = predecessor generator's
  degree + differential path length — there is no second route), so the anchor MUST be Berger's
  formula, not "another extractor agrees". Keep `_NonPure` as a raising assert with an honest
  comment saying so.
- **Path length of a `P_v` basis label — ROBUST, not string-counting.** The label is `e_v`
  (length 0) or a `*`-joined arrow word. Compute the length via the **quiver arrow set** (count
  the label's arrow letters against `A.quiver.arrows`), NOT a bare `label.count("*")+1` (which
  is fragile if an arrow name ever contains `*`). Add a module-internal
  `_projective_basis_lengths(A, v)` that returns the length of each `projective(A, v)`
  basis label, memoized per `(A, v)`. (Adjust-to-reality: if `projective` grows a first-class
  per-label length accessor, use it and delete the helper.)

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_nkoszul_generation_degrees.py
"""Internal (path-length) generation degrees of Ext(k,k) off the shipped minimal resolutions
(Plan 77 / R36). Validated against Berger's closed form on k[x]/(x^N): delta(n) = (N/2) n
(n even) / (N/2)(n-1)+1 (n odd)."""
import pytest
from quiverlab import truncated_polynomial, TruncatedPathAlgebra
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import generation_degrees, n_homogeneous_degree

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _berger(n, N):
    return (n // 2) * N if n % 2 == 0 else (n // 2) * N + 1


@lit
@pytest.mark.parametrize("N,length,expected", [
    (2, 8, [0, 1, 2, 3, 4, 5, 6, 7, 8]),
    (3, 8, [0, 1, 3, 4, 6, 7, 9, 10, 12]),
    (4, 8, [0, 1, 4, 5, 8, 9, 12, 13, 16]),
    (5, 8, [0, 1, 5, 6, 10, 11, 15, 16, 20]),
])
def test_kxN_internal_degrees_match_berger(N, length, expected):
    A = truncated_polynomial(N, field=QQ)
    gd = generation_degrees(A, 1, length)         # single vertex -> S_1
    flat = [s[0] for s in gd]                      # one summand each degree
    assert flat == expected[:len(flat)]
    assert flat == [_berger(n, N) for n in range(len(flat))]
    assert n_homogeneous_degree(A) == N


@lit
def test_cubic_monomial_multivertex_is_berger_N3():
    A = TruncatedPathAlgebra("A4", 3, field=QQ)    # single relation e12*e23*e34, N=3
    assert n_homogeneous_degree(A) == 3
    gd = generation_degrees(A, 1, 4)               # S_1 pd = 2
    assert [s[0] for s in gd if s] == [0, 1, 3]    # Berger delta for N=3, finite


@selfcert
def test_generation_degrees_single_element_on_pure_resolution():
    # DEFENSIVE (H3): on k[x]/x^3 every P_n is generated in ONE internal degree, so each returned
    # set is a singleton. This is NOT the correctness anchor (that is Berger's closed form, above)
    # -- it guards the _NonPure branch stays dormant on a genuinely graded-pure input.
    A = truncated_polynomial(3, field=QQ)
    for s in generation_degrees(A, 1, 6):
        assert len(s) == 1


@selfcert
def test_refuses_non_length_graded():
    from quiverlab import Quiver
    from quiverlab.errors import QuiverlabError
    # x^3 - x^2 is ADMISSIBLE (ideal inside rad^2) but genuinely non-length-graded (a length-3
    # and a length-2 term) -- verified at authoring it BUILDS and _is_length_graded==False.
    # (x^2 - x, the naive choice, raises AdmissibilityError AT CONSTRUCTION -- ideal not in
    # rad^2 -- so it would never exercise the _is_length_graded gate: a false green.)
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x*x - x*x"], field=QQ)
    with pytest.raises(QuiverlabError, match="length-graded"):   # the GATE, scoped to the call
        generation_degrees(A, 1, 4)
```

- [ ] **Step 2: Run to verify failure** (`ModuleNotFoundError: quiverlab.modules.nkoszul`).
- [ ] **Step 3: Implement** `n_homogeneous_degree`, `_projective_basis_lengths`,
  `generation_degrees`, `_NonPure`. Gate on `koszul._is_length_graded`; propagate
  `DepthLimitError` from `minimal_resolution` honestly (a syzygy that overshoots `max_term_dim`
  ⇒ the window is capped, reported by the caller).
- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(modules): nkoszul.generation_degrees -- internal generation degrees of Ext(k,k) off the minimal resolution, Berger-validated on k[x]/x^N`

**Adjust to reality (Task 1):** confirm `_pv_basis_labels` are `*`-joined arrow words with the
trivial path labelled `e_<v>` (verified at authoring: `k[x]/x^N` → `['e_1','x','x*x',…]`,
`Π(A4)` `P_1` → `['e_1','e12','e12*e23','e12*e23*e34']`). If a later refactor changes the label
grammar, adjust `_projective_basis_lengths` (the Berger pins are the arbiter).

---

### Task 2: N-Koszul recognizer (Berger) + the two-certificate cross-check

**Files:**
- Modify: `src/quiverlab/modules/nkoszul.py` (add `n_koszul_certificate`)
- Test: `tests/modules/test_nkoszul_ncertificate.py`

**Interfaces:**
```python
def n_koszul_certificate(A, top=8) -> dict:
    # {"n_homogeneous": N|None,       # single relation degree, or None (mixed/inhomogeneous)
    #  "verdict": True|False|None,    # N-Koszul: True (certified), False (pattern break),
    #                                 #   None (inconclusive beyond window / not N-homogeneous)
    #  "window": W,                   # certified_through_degree of the Yoneda engine
    #  "complete": bool,              # gl.dim finite (exact) -> window complete
    #  "berger_expected": [..],       # delta(n) for n=0..W
    #  "internal_degrees": {i: [..]}, # per-simple generation degrees through W
    #  "k2_agrees": bool|None,        # cross-check: (N-homogeneous AND K2) == purity verdict
    #  "reason": str}
```

**Logic.**
- `N = n_homogeneous_degree(A)`. If `N is None` (mixed relation degrees) ⇒ **not N-homogeneous**;
  `verdict=None`, `reason="relations are not homogeneous of a single degree; try K2/multi-Koszul"`.
- **`N == 2` DEFERS TO PLAN 27** (the overlap). **This branch also absorbs the HEREDITARY case**
  (`relations == []` ⇒ `n_homogeneous_degree` returns `2` per its docstring, quadratic-trivially),
  so there is no separate empty-relations fall-through: `verdict = A.ext_algebra(top).koszul`
  verbatim, `reason="N=2: quadratic Koszulity is Plan 27's g_quadratic_certificate/Froberg
  verdict (hereditary kQ is quadratic-trivially)"`. Plan 77 adds no second opinion in the
  quadratic case. (kA₂: `n_homogeneous == 2`, `verdict == True` = Plan-27's Koszul.)
- `N ≥ 3`: for every simple `S_i`, compute `generation_degrees(A, i, W)`. **Berger purity
  requires each `ℓ_i(n)` to be the SINGLETON `{δ(n)}`** (`P_n` generated in one internal
  degree). **`M4` — a MULTI-ELEMENT `ℓ_i(n)` set (a `P_n` with summand generators in different
  internal degrees) is a purity FAILURE ⇒ `verdict=False`**, `reason="P_n of S_i is generated in
  multiple internal degrees {..} != the single Berger degree δ(n): the resolution is not pure"`.
  Otherwise compare the singleton to Berger `δ(n)`: **`verdict=True`** iff every
  `ℓ_i(n) == {δ(n)}` through `W` **and** (`complete` — the finite-gl.dim resolution terminates
  within the pure pattern — **or** K₂ holds AND N-homogeneous, Berger's `N≥3` characterization).
  **`verdict=False`** at the first `ℓ_i(n) ≠ {δ(n)}` (a genuine break — e.g. `Π(A₃)` at `n=3`;
  but `Π(A₃)` is `N=2`, so it defers to Plan 27 — a true `N≥3` break example is a cubic monomial
  with a non-Berger overlap). **`verdict=None`** when the pattern holds through `W` but `gl.dim`
  is not exact-finite and K₂ is only window-certified (honest inconclusive beyond).
- **`k2_agrees`** = the `oracle_crossengine` self-check: `(N-homogeneous ∧ k2_certificate.verdict)`
  matches the purity verdict on the same window (Berger's two certificates agree).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_nkoszul_ncertificate.py
"""N-Koszul recognizer (Berger) + the Ext-generation cross-check (Plan 77 / R36).
k[x]/(x^N), N>=3 -> N-Koszul (pure resolution, delta(n)); quadratic defers to Plan 27; the
preprojectives are NOT N-Koszul (pattern break at n=3)."""
import pytest
from quiverlab import truncated_polynomial, TruncatedPathAlgebra, PreprojectiveAlgebra
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import n_koszul_certificate

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("N", [3, 4, 5])
def test_kxN_is_N_koszul(N):
    c = n_koszul_certificate(truncated_polynomial(N, field=QQ), top=8)
    assert c["n_homogeneous"] == N
    # k[x]/x^N is self-injective (gl.dim infinite): verdict True-through-window, backed by K2
    assert c["verdict"] in (True, None)            # certified through window (K2 holds)
    assert c["berger_expected"][:4] == [0, 1, N, N + 1]
    assert c["internal_degrees"][1][:4] == [0, 1, N, N + 1]


@lit
def test_cubic_monomial_finite_gldim_is_3koszul():
    c = n_koszul_certificate(TruncatedPathAlgebra("A4", 3, field=QQ), top=6)
    assert c["n_homogeneous"] == 3 and c["verdict"] is True and c["complete"] is True


@xeng
def test_quadratic_defers_to_plan27():
    # kA2 (hereditary) is quadratic -> N=2, verdict == ext_algebra.koszul
    from quiverlab import linear_path_algebra
    A = linear_path_algebra(2, field=QQ)
    c = n_koszul_certificate(A, top=6)
    assert c["n_homogeneous"] == 2
    assert c["verdict"] == A.ext_algebra(6).koszul


@lit
def test_preprojective_A3_is_not_N_koszul():
    # Pi(A3) is quadratic (N=2) but NOT Koszul: defers to Plan 27's False.
    c = n_koszul_certificate(PreprojectiveAlgebra("A3", field=QQ), top=6)
    assert c["n_homogeneous"] == 2 and c["verdict"] is False


@xeng
@pytest.mark.parametrize("N", [3, 4, 5])
def test_two_certificates_agree(N):
    # oracle_crossengine: purity verdict == (N-homogeneous AND K2). Berger's N>=3 theorem.
    c = n_koszul_certificate(truncated_polynomial(N, field=QQ), top=8)
    assert c["k2_agrees"] is True
```

- [ ] **Step 2–4:** run-fail, implement, run-pass.
- [ ] **Step 5: Commit** — `feat(modules): nkoszul.n_koszul_certificate -- Berger 2-N alternation + Ext-generation cross-check; N=2 defers to Plan 27`

**Adjust to reality (Task 2):** the `verdict True vs None` split for the self-injective
`k[x]/x^N` (infinite gl.dim) hinges on whether K₂ (the Berger `N≥3` characterization) certifies
completeness through the window — assert `verdict in (True, None)` and pin the internal degrees
(the decisive, window-free data) rather than over-committing the three-valued flag. The finite
`gl.dim` cases (`kA₄/J³`) are unambiguously `True`.

---

### Task 3: K₂ recognizer (Cassidy–Shelton) + the explicit certified window

**Files:**
- Modify: `src/quiverlab/modules/nkoszul.py` (add `k2_certificate`)
- Test: `tests/modules/test_nkoszul_k2.py`

**Interfaces:**
```python
def k2_certificate(A, top=8) -> dict:
    # {"verdict": True|False|None, "window": W, "complete": bool,
    #  "generator_degrees": [d,...],   # the homological degrees carrying Yoneda generators
    #  "reason": str}
    # K2 (Cassidy-Shelton): E(A)=Ext(k,k) generated in cohomological degrees 1,2.
    # Reads YonedaPresentation.generators_by_degree (Plan 27) through certified_through_degree.
```

**Logic.** `Y = A.ext_algebra(top)`; `degs = sorted(Y.generators_by_degree)`; `W =
Y.certified_through_degree`; `complete = Y.is_finite_dimensional is True` (gl.dim finite exact).
- **`verdict=False`** iff some generator sits in degree `≥ 3` **within `W`** (a genuine
  obstruction — decisive regardless of completeness; e.g. `Π(A₃)` has a degree-3 generator).
- **`verdict=True`** iff `degs ⊆ {1, 2}` **and** `complete` (the whole finite `E(A)` is
  generated in degrees 1,2).
- **`verdict=None`** iff `degs ⊆ {1, 2}` through `W` but **not** `complete` — the honest
  certified window: "K₂ holds through degree `W`; a generator could still appear beyond it since
  `gl.dim A` is not finite" (`reason` states `W`).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_nkoszul_k2.py
"""K2 (Cassidy-Shelton): E(A)=Ext(k,k) generated in cohomological degrees 1,2, decided through
an explicit certified window (Plan 77 / R36)."""
import pytest
from quiverlab import truncated_polynomial, TruncatedPathAlgebra, PreprojectiveAlgebra
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import k2_certificate

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("N", [3, 4, 5])
def test_kxN_is_K2_through_window(N):
    c = k2_certificate(truncated_polynomial(N, field=QQ), top=8)
    assert c["generator_degrees"] == [1, 2]        # K2: gens in degrees 1,2
    assert c["verdict"] in (True, None)            # window-certified (self-injective)
    assert c["window"] >= 8


@lit
def test_finite_gldim_K2_is_complete_true():
    c = k2_certificate(TruncatedPathAlgebra("A4", 3, field=QQ), top=6)
    assert c["generator_degrees"] == [1, 2]
    assert c["verdict"] is True and c["complete"] is True


@lit
def test_preprojective_A3_is_NOT_K2():
    # THE discriminator: Pi(A3) has a genuine Yoneda generator in cohomological degree 3.
    c = k2_certificate(PreprojectiveAlgebra("A3", field=QQ), top=6)
    assert 3 in c["generator_degrees"]
    assert c["verdict"] is False                   # decisive: a degree-3 generator


@selfcert
def test_koszul_implies_K2():
    # Koszul (degree-1 gens) => K2. kA2 hereditary.
    from quiverlab import linear_path_algebra
    c = k2_certificate(linear_path_algebra(2, field=QQ), top=6)
    assert c["generator_degrees"] == [1] and c["verdict"] is True
```

- [ ] **Step 2–4:** run-fail, implement, run-pass.
- [ ] **Step 5: Commit** — `feat(modules): nkoszul.k2_certificate -- Cassidy-Shelton K2 with an explicit certified window (three-valued, honest beyond)`

**Adjust to reality (Task 3):** an empty `generators_by_degree` (semisimple `A`, no arrows) has
`degs = []` ⊆ {1,2} — Koszul/K₂ vacuously; the `verdict` follows the `complete` flag. Confirm
`Y.is_finite_dimensional` is the right completeness signal (`True` when `gd.exact`, else
`None` — the Plan-27 field), not a recompute.

---

### Task 4: (p,q)-almost-Koszul classifier (BBK) + the P75 seam

**Files:**
- Modify: `src/quiverlab/modules/nkoszul.py` (add `almost_koszul_certificate`, `_top_degree`)
- Test: `tests/modules/test_nkoszul_almost.py`

**Interfaces:**
```python
def _top_degree(A) -> int:
    # p = top graded degree of A (concentration) = the highest path length over A's OWN basis:
    #     max(_projective_basis_lengths-style path length of `label` for label in A.basis_labels).
    # A.basis_labels IS the irreducible-path basis of A (e_v = length 0, "*"-joined arrow words
    # otherwise) -- the length-graded concentration is exactly its max length (= Loewy length -1).
    # NO Groebner completion involved. Requires length-graded A.
    #
    # WHY NOT koszul._algebra_graded_matrices: it re-completes the reduction system under the
    # certificate's default degree_bound=8, which is INSUFFICIENT for Pi(A5)/Pi(D4) (their tips
    # reach length >= 9) -- it raises AdmissibilityError live (verified at the fix round). The
    # basis-label route never completes anything, so it is total on every finite-dim A.

def almost_koszul_certificate(A, top=8) -> dict:
    # {"verdict": True|False|None,      # (p,q)-almost-Koszul recognized
    #  "p": p, "q": q,                  # p = top degree; q = e - p (first-jump internal - p)
    #  "break_hom_degree": n_star,      # first homological n with internal(n) > n
    #  "break_internal_degree": e,      # internal degree there (= p + q)
    #  "period": int|None,              # engine-observed period of the jump pattern (window)
    #  "window": W, "complete": bool,
    #  "seam_obstruction_degree": d,    # ext_algebra.koszul obstruction degree (the P75 refusal)
    #  "reason": str}
```

**Logic — the SIGNATURE the recognizer checks (NOT BBK's full finite-complex condition).**
Requires length-graded `A`, `A₀` semisimple, generated in degree 1 (loud refusal otherwise).
`p = _top_degree(A)`. Walk `generation_degrees(A, i, W)` for each simple. Because every
differential has length `≥ 1`, `ℓ_i(n) ≥ n` always, with `ℓ_i(n) = n` exactly on a **linear
(Koszul) strand**. Define:
- the **linear prefix** = the initial run `ℓ_i(n) = n` for `n < n★`;
- the **first jump** `n★` = the least `n` with `ℓ_i(n★) > n★`; `e = ℓ_i(n★)` (must agree across
  simples — else `verdict=None`, a non-uniform algebra);
- `q = e − p`.

**The recognizer decides the (p,q)-almost-Koszul SIGNATURE:** a linear prefix, then a first jump
to internal degree `e = p + q`. **It does NOT verify BBK's full definition** (the linear complex
of projectives resolving `S` up to a single error, and the `2(h−1)` periodicity) — that heavier
verifier is **explicitly scoped OUT** and stated as such; the signature is what BBK's definition
implies pointwise and is exactly what the engine can certify from the shipped resolution.

- **`verdict=True`** iff a linear prefix exists, a first jump exists with `e` agreeing across
  simples, **`q ≥ 2`**, and `A` is concentrated in `0..p`. The engine-observed `period` (distance
  between successive jumps through `W`) is recorded with an honest note it is **window-observed**,
  not proven (BBK's `2(h−1)` is cited, not asserted from the window).
- **`verdict=None`** iff no first jump through `W` (a genuinely Koszul linear resolution — e.g.
  `Π(A₂)`, `rad²=0`), **or `q = 1`** (see the negative pin), or `e` disagrees across simples.
- **The `q ≥ 2` gate is BBK's proper almost-Koszul boundary — the negative pin (`M3`).** For
  `k[x]/(x^N)` the top degree is `p = N−1`, the first jump lands at `e = N`, so **`q = e − p =
  1` for EVERY `N`** (engine-verified `N=2,3,4`). `q = 1` means the "error" sits at degree
  `p+1` — the *N-Koszul / Koszul* regime, **not** almost-Koszul. So `almost_koszul_certificate`
  **returns `verdict=None` for `k[x]/(x^N)`** (it is caught by `n_koszul_certificate` instead),
  and `Π(A₃)` (`q = 4−2 = 2`), `Π(A₄)` (`q = 5−3 = 2`), … are the `q = 2` almost-Koszul cases.
  This matches BBK: `(p,1)` is the Koszul-type degenerate boundary, `q ≥ 2` is "genuinely almost
  Koszul".
- **`Π(A₂)` vs BBK's formal `(h−2,2) = (1,2)` — reconciled, not asserted.** `Π(A₂) = kZ₂/rad²`
  is **radical-square-zero**, hence a **quadratic monomial (Priddy G-quadratic) algebra — it is
  genuinely Koszul** (`ext_algebra.koszul == True`, engine-verified; linear resolution
  `ℓ(n) = n` forever, no first jump). BBK's periodic-preprojective theorem is a statement about
  the *non-Koszul* Dynkin range (type `A_n`, `n ≥ 3`); A₂'s nominal `(1,2)` label is **subsumed
  by genuine Koszulity** (a Koszul algebra needs no "error"). The recognizer therefore reports
  `verdict=None` for `Π(A₂)` (no almost-Koszul break) while the **Koszul/K₂ fields report it
  Koszul** — the honest split, and the reason string says exactly this.
- **`seam_obstruction_degree`** = `A.ext_algebra(top).koszul_obstruction[0]` (present iff the
  Plan-27 verdict is `False`). **`W3`: this is a STRUCTURAL-CONSISTENCY check, not two
  independent computations** — the `ext_algebra.koszul` obstruction degree and the
  almost-Koszul `break_hom_degree` are *the same first-failure-of-linearity degree* read two
  ways (the Yoneda engine's "a new generator appears" ⟺ the resolution's first jump above the
  diagonal). The `oracle_selfcert` pin `seam_obstruction_degree == break_hom_degree` names the
  P75 seam (P75 REFUSES on that degree; P77 CLASSIFIES from it) as a consistency invariant, not
  a cross-engine oracle.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_nkoszul_almost.py
"""(p,q)-almost-Koszul (Brenner-Butler-King): p = top degree, q = e - p (first internal jump
minus p). Reproduces BBK's (h-2, 2) on the Dynkin preprojectives. Classifies exactly the
algebras P75's GHMS-Koszul route refuses (Plan 77 / R36)."""
import pytest
from quiverlab import PreprojectiveAlgebra
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import almost_koszul_certificate

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("typ,p,h", [("A3", 2, 4), ("A4", 3, 5), ("A5", 4, 6), ("D4", 4, 6)])
def test_preprojective_is_p_2_almost_koszul(typ, p, h):
    c = almost_koszul_certificate(PreprojectiveAlgebra(typ, field=QQ), top=6)
    assert c["verdict"] is True
    assert c["p"] == p == h - 2                     # BBK: p = h-2 (top degree)
    assert c["q"] == 2                              # BBK: q = 2
    assert c["break_hom_degree"] == 3               # first non-linear jump (engine-observed)
    assert c["break_internal_degree"] == p + 2      # e = p + q


@lit
def test_preprojective_A2_is_koszul_not_almost():
    # Pi(A2) = kZ2/rad^2 is genuinely Koszul (linear forever): no break -> verdict None.
    # (BBK's nominal (1,2) label is subsumed by genuine Koszulity -- see Logic reconciliation.)
    from quiverlab import PreprojectiveAlgebra
    A = PreprojectiveAlgebra("A2", field=QQ)
    c = almost_koszul_certificate(A, top=8)
    assert c["verdict"] is None                     # no almost-Koszul break (it IS Koszul)
    assert A.ext_algebra(8).koszul is True          # the Koszul field reports it Koszul


@lit
@pytest.mark.parametrize("N", [2, 3, 4])
def test_kxN_is_not_almost_koszul_q_equals_1(N):
    # THE negative pin (q>=2 gate). k[x]/x^N has top degree p=N-1 and first jump e=N, so
    # q = e - p = 1 for EVERY N -> NOT almost-Koszul (it is N-Koszul; caught by
    # n_koszul_certificate). q=1 is BBK's Koszul-type degenerate boundary, not "almost".
    from quiverlab import truncated_polynomial
    c = almost_koszul_certificate(truncated_polynomial(N, field=QQ), top=8)
    assert c["p"] == N - 1 and c["q"] == 1
    assert c["verdict"] is None                     # q=1 refused: N-Koszul, not almost-Koszul


@selfcert
@pytest.mark.parametrize("typ", ["A3", "A4"])
def test_p75_seam_obstruction_equals_break(typ):
    # THE P75 seam (card): the ext_algebra.koszul obstruction degree that P75's GHMS route
    # REFUSES on equals the almost-Koszul break degree that P77 CLASSIFIES.
    c = almost_koszul_certificate(PreprojectiveAlgebra(typ, field=QQ), top=6)
    assert c["seam_obstruction_degree"] == c["break_hom_degree"] == 3
```

- [ ] **Step 2–4:** run-fail, implement, run-pass.
- [ ] **Step 5: Commit** — `feat(modules): nkoszul.almost_koszul_certificate -- BBK (p,q) via p=top-degree, q=e-p; reproduces (h-2,2) on Dynkin preprojectives; the P75 classification seam`

**Adjust to reality (Task 4):** `q = e − p` reproduces BBK `(h−2, 2)` for all preprojectives
tried (Live-verified pins C). If a reviewer wants the FULL BBK certificate (the linear-complex
error + the `2(h−1)` periodicity proven, not window-observed), that is a heavier verifier —
scope it OUT with an honest note (`period` is window-observed; the classification is the
`(p, q)` label + the single-break structure, which is what BBK's definition checks pointwise).
`_top_degree` is `max` path length over `A.basis_labels` (NOT `koszul._algebra_graded_matrices`,
which re-completes under `degree_bound=8` and raises `AdmissibilityError` on `Π(A₅)`/`Π(D₄)` —
verified at the fix round; the basis-label route is total, giving `Π(A₂..A₆)` ⇒ `1,2,3,4,5` and
`Π(D₄)` ⇒ `4`, all `= h−2`).

---

### Task 5: multi-Koszul recognizer (scoped, connected-graded) + the K₂ transfer

**Files:**
- Modify: `src/quiverlab/modules/nkoszul.py` (add `multi_koszul_certificate`, `_is_connected_graded`)
- Test: `tests/modules/test_nkoszul_multi.py`

**Interfaces:**
```python
def multi_koszul_certificate(A, top=8) -> dict:
    # {"applicable": bool,       # connected graded (single vertex, A_0 = k)?
    #  "verdict": True|False|None,
    #  "k2": {...},              # the transferable property (multi-Koszul => K2)
    #  "generation_degrees": {1: [..]},   # the headline primitive on the local resolution
    #  "status": str,            # e.g. "requires the verbatim Herscovich criterion (# PIN)"
    #  "reason": str}
```

**Logic (the settlement, § Multi-Koszul settlement).**
- `_is_connected_graded(A)` = single vertex (`len(A.quiver.vertices) == 1`) ∧ length-graded
  (`A₀ = k`). **Multi-vertex ⇒ `applicable=False`**, `verdict=None`,
  `reason="multi-Koszul is a connected-graded (A_0=k) notion (Herscovich 1305.1678); for
  multi-vertex kQ/I the transferable property is K2 -- see the k2 field"`, and `k2 =
  k2_certificate(A, top)`. **No fabricated verdict.**
- Connected (local) input: report `generation_degrees(A, 1, W)` (the headline) + `k2` (which
  multi-Koszul implies, **Prop. 3.30**). **`verdict` is `None`** — the full multi-Koszul
  *decision* (Herscovich §3.2: Tor/Ext-vanishing of the minimal graded **bimodule** resolution)
  is **scoped OUT** of this plan (a bimodule-Tor engine, not a recognizer over shipped surfaces),
  with `status="the multi-Koszul DECISION requires bimodule Tor/Ext vanishing per Herscovich
  §3.2 -- out of this plan's scope; K2 (Prop 3.30) and the generation degrees ARE reported"`.
  **No fabricated verdict, ever.** The N-Koszul case (`k[x]/x^N`, a single-relation-degree
  multi-Koszul) is the sanity anchor: its `k2` is True-through-window and its generation degrees
  are Berger's `δ(n)`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_nkoszul_multi.py
"""Multi-Koszul (Herscovich): a CONNECTED-graded (A_0=k) notion; multi-vertex kQ/I refuses and
reports K2 (the settled transfer). Plan 77 / R36."""
import pytest
from quiverlab import truncated_polynomial, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import multi_koszul_certificate

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
def test_multivertex_refuses_and_points_at_K2():
    A = linear_path_algebra(3, field=QQ)           # 3 vertices -> A_0 = k^3, not connected
    c = multi_koszul_certificate(A, top=6)
    assert c["applicable"] is False and c["verdict"] is None
    assert "K2" in c["reason"] and c["k2"]["verdict"] is not None or c["k2"] is not None


@lit
def test_local_kxN_is_connected_and_reports_generation_and_K2():
    A = truncated_polynomial(3, field=QQ)          # single vertex -> connected graded
    c = multi_koszul_certificate(A, top=8)
    assert c["applicable"] is True
    assert c["generation_degrees"][1][:4] == [0, 1, 3, 4]   # Berger N=3
    assert c["k2"]["generator_degrees"] == [1, 2]
    # verdict is None: the multi-Koszul DECISION (Herscovich §3.2 bimodule Tor/Ext) is scoped
    # out; K2 (Prop 3.30) + generation degrees are what v1 reports.
    assert c["verdict"] is None and "3.2" in c["status"] and "scope" in c["status"].lower()
```

- [ ] **Step 2–4:** run-fail, implement, run-pass. The transfer is **Prop. 3.30** (multi-Koszul
  ⟹ K₂, quoted verbatim in the docstring from the ar5iv HTML — already re-read at the fix round,
  not a pending fetch). The multi-Koszul **decision** stays scoped out (`verdict=None` + the
  §3.2 status), so there is no verdict to "flip".
- [ ] **Step 5 (BLOCKING before merge): record the scoped-out decision in the ledger.** Add a
  named entry to `docs/plans/DEEPER-ENGINES-BACKLOG.md` — *"multi-Koszul DECISION (Herscovich
  §3.2): a bimodule Tor/Ext-vanishing engine; P77 ships K₂ (Prop 3.30) + generation degrees +
  the connected-only scope refusal, `verdict=None` for the decision itself"* (metaplan §1.2
  GUI-deferral-ledger discipline — a named deferral, never silent). **No guessed verdict is ever
  emitted.**
- [ ] **Step 6: Commit** — `feat(modules): nkoszul.multi_koszul_certificate -- connected-graded scope, K2 transfer (Herscovich Prop 3.30) for multi-vertex, decision scoped out (§3.2 bimodule Tor)`

---

### Task 6: the assembled profile + QPA honest-scope probe

**Files:**
- Modify: `src/quiverlab/modules/nkoszul.py` (`koszul_profile`, `koszul_profile_block`),
  `src/quiverlab/core/algebra.py` (`Algebra.koszul_profile(top=8)` lazy delegate)
- Create: `tests/qpa/test_nkoszul_qpa.py`
- Test: `tests/modules/test_nkoszul_profile.py`

**Interfaces:**
```python
def koszul_profile(A, top=8) -> dict          # the full record (dataclass-like dict)
def koszul_profile_block(A, top=8) -> dict:
    # the no-code block SHARED by both runners (byte-identical); each runner adds `citations`.
    # {"kind": "koszul", "top": top,
    #  "quadratic_koszul": bool|None,          # Plan 27's ext_algebra.koszul (the overlap)
    #  "quadratic_reason": str, "quadratic_obstruction": [d, msg]|None,
    #  "n_homogeneous": N|None, "n_koszul": {...},
    #  "k2": {...}, "almost_koszul": {...}, "multi_koszul": {...},
    #  "generation_degrees": {str(i): [..]},   # per-simple internal degrees through the window
    #  "certified_through_degree": W, "complete": bool,
    #  "latex": "...",                          # a compact "A is N-Koszul / K2 / (p,q)-Koszul" line
    #  "references": [keys...]}
# Algebra.koszul_profile(top=8) delegates to koszul_profile.
```

**Logic.** `koszul_profile_block` calls `A.ext_algebra(top)` ONCE, threads the
`YonedaPresentation` into every sub-certificate (avoid recomputing the resolution), assembles
the four verdicts + the generation-degree table, and builds the `latex` summary line. The
`quadratic_*` fields are Plan-27's verdict verbatim (the NAMED overlap). `references =
["berger_nonquadratic", "cassidy_shelton", "brenner_butler_king", "herscovich_multikoszul",
"chouhy_degenerations", "priddy", "froberg_koszul"]` (Plan-27 keys included for the quadratic
part).

**QPA scope.** QPA 1.37 has **NO Koszul surface of any kind** (Plan-27 finding: "QPA has no
IsKoszul"; no `IsNKoszul`, `IsK2Algebra`, `IsAlmostKoszul`, `MultiKoszul`). Probe live via
`NamesGVars()` (the Plan-27/45 precedent): SKIP honestly, **FAIL if QPA ever ships** any such
verb, so the honest-scope claim cannot silently rot.

- [ ] **Step 1: Write the failing tests** (`test_nkoszul_profile.py`: the block shape on
  `k[x]/x³` (N-Koszul, K₂), `Π(A₃)` (almost-Koszul, not-K₂), `kA₂` (quadratic-Koszul defers to
  P27); references present; `Algebra.koszul_profile` delegates). QPA probe: skip + fail-if-shipped.

```python
# tests/modules/test_nkoszul_profile.py (excerpt)
import pytest
from quiverlab import truncated_polynomial, PreprojectiveAlgebra
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import koszul_profile_block

lit = pytest.mark.oracle_literature


@lit
def test_profile_kx3():
    b = koszul_profile_block(truncated_polynomial(3, field=QQ), top=8)
    assert b["kind"] == "koszul"
    assert b["quadratic_koszul"] is False                # Plan 27's verdict (N=2 test fails)
    assert b["n_homogeneous"] == 3
    assert b["k2"]["generator_degrees"] == [1, 2]
    assert b["generation_degrees"]["1"][:4] == [0, 1, 3, 4]
    assert "berger_nonquadratic" in b["references"]


@lit
def test_profile_preprojective_A3_almost():
    b = koszul_profile_block(PreprojectiveAlgebra("A3", field=QQ), top=6)
    assert b["quadratic_koszul"] is False
    assert b["almost_koszul"]["verdict"] is True
    assert b["almost_koszul"]["p"] == 2 and b["almost_koszul"]["q"] == 2
    assert b["k2"]["verdict"] is False                   # a degree-3 generator
```

- [ ] **Step 2–4:** run-fail, implement `koszul_profile`/`koszul_profile_block` + the
  `Algebra.koszul_profile` delegate + the QPA probe, run-pass.
- [ ] **Step 5: Commit** — `feat(modules,core): nkoszul.koszul_profile[_block] + Algebra.koszul_profile delegate; QPA honest no-Koszul-surface probe`

---

### Task 7: no-code exposure — the `koszul` compute kind (GUI panel extension)

The organizing invariant is **byte-identity across the two runners**, achieved by routing BOTH
through the single shared builder `nkoszul.koszul_profile_block` (the `ext_algebra` precedent —
each runner only adds `citations` locally). Anchors below were grepped against the live tree at
authoring; still "adjust to reality" if a later merge shifts them.

**`koszul` is an ALGEBRA-level compute kind** carrying an optional TOP DEGREE (default 8),
parsed exactly like `ext_algebra`: `koszul` or `koszul:10`. NOT a `MODULE_KINDS` entry; sized on
`A.dim` like `ext_algebra`. **No `koszul` kind exists yet** (grepped — free name).

**Files (the seven-touchpoint checklist, mirroring `ext_algebra`):**
1. **Parse grammar** — `hpc/spec.py::parse_compute_item` (the `ext_algebra` branch): add
   `koszul` / `koszul:<n>` → `ComputeItem(kind="koszul", lo=None, hi=top)`. Mirror in
   `docs/gui/runner.py` parse.
2. **Server/HPC dispatch** — `spec.py::_dispatch` (beside the `ext_algebra` branch, `spec.py:1789`):
   `if kind == "koszul": top = item.hi or 8; from quiverlab.modules.nkoszul import
   koszul_profile_block; block = koszul_profile_block(A, top); block["citations"] =
   _citation_pairs(block["references"]); return block, None`. Add the `_snippet` entry
   (`spec.py:2706`, beside `ext_algebra`): `"koszul": lambda it:
   f"A.koszul_profile({it.hi if it.hi is not None else 8})"`. Catch the presentation/length-graded
   `QuiverlabError` into `{"error": "<loud message>"}` (never a 500).
3. **Pyodide twin** — `docs/gui/runner.py::compute_one` (beside the `ext_algebra` branch,
   `runner.py:1084`): the byte-identical `elif name == "koszul":` branch calling
   `koszul_profile_block`; one `calls` reproduce entry (`runner.py:1456`): `"koszul":
   "A.koszul_profile()"`; one `ETA_MODEL["scalars"]` cost (`runner.py:1593`, beside
   `"ext_algebra": 2.0`) — `"koszul": 2.2` (ext_algebra + the four recognizers, a little above).
4. **GUI JS** — `docs/gui/gui.js` AND `webapp/static/gui/gui.js` (diff-identical — apply every
   edit to both): a checkbox + top-degree picker in the SAME structure panel that hosts
   `ext_algebra` (the "koszul panel extension" of the card — put `koszul` right beside
   `ext_algebra`); the id-registry entries; the request push-list entry
   (`if (el.koszul.checked) compute.push("koszul:" + el["koszul-top"].value)`); a `renderBlock`
   `else if (name === "koszul")` branch calling a new `renderKoszul(div, b)`; the `THEMES`
   structure `kinds` array entry `"koszul"` inside the `QLGUI-THEMES-BEGIN/END` sentinels; the
   `scheduleProbe`/`KIND_CTRL` entries. **`renderKoszul`** renders (a) the headline verdict line
   (`A is Koszul / N-Koszul (N=..) / K₂ / (p,q)-almost-Koszul / not …`), (b) the
   **generation-degree table** (per simple, `ℓ(n)` vs the homological degree — the visual of the
   2-N alternation), (c) the four sub-verdicts with their windows + honest "through degree W"
   labels, and (d) the Plan-27 quadratic obstruction when present. MathJax for the `δ(n)` line.
5. **i18n (FOUR locales) + picker template** — `webapp/server/i18n/{en,es,fr,zh}.json` (the
   shipped `LANGS`): add `pick.kind.koszul` (REQUIRED — gated by `test_layout_picker.py`),
   `struct.koszul` (checkbox label), and `kz.*` block strings (`kz.title`, `kz.koszul`,
   `kz.nkoszul`, `kz.k2`, `kz.almost`, `kz.multi`, `kz.gen_degrees`, `kz.window`,
   `kz.through_w`, `kz.top`, `kz.not_applicable`). Add one
   `data-pick-kind-koszul="{{ t('pick.kind.koszul') }}"` line in `webapp/templates/draw.html`.
6. **Report renderer** — `src/quiverlab/trace/results_html.py`: `_HEADINGS["koszul"]` + a
   `_koszul_html(b)` branch in `_block_html` (beside the `ext_algebra` branch) rendering the
   headline verdict, the generation-degree table (an indexed grid — the Plan-34 `matrix_grid`
   idiom), the four sub-verdicts with windows, and the quadratic obstruction. A TikZ twin
   `viz/tikz.py::tikz_koszul(b)` (the 2-N alternation staircase for `n ≤ W`; an honest "window
   W" node) mirrors the fan/tikz precedent.
7. **Tests / gates** — add `koszul` to `ALL_KINDS` in `tests/webapp/test_layout_picker.py`
   (THEMES == `ALL_KINDS`, each once; the four-locale `pick.kind.*` gate); add ONE golden
   `koszul_kx3` to `tests/webapp/_runner_goldens.json` with a dated change-log bullet in
   `tests/webapp/test_runner_delegation.py` (verify existing entries byte-identical FIRST; on
   rebase merge the JSON semantically, never textually); add a twin-parity test.
- Test files: `tests/webapp/test_koszul_p77.py`, `tests/gui/test_koszul_runner_twin.py`.

```python
# tests/webapp/test_koszul_p77.py (excerpt)
def test_koszul_block_shape(tmp_path):
    from quiverlab.hpc.spec import ComputeRequest, run
    req = _koszul_request(quiver=("kx3",), top=8)     # k[x]/x^3
    out = run(ComputeRequest.model_validate(req), tmp_path)
    b = out["results"]["koszul"]
    assert b["n_homogeneous"] == 3 and b["k2"]["generator_degrees"] == [1, 2]
    assert "berger_nonquadratic" in b["references"]

def test_twin_parity(tmp_path):
    # json.dumps(block, sort_keys=True) equality: docs/gui/runner.py vs hpc.spec on `koszul`.
    ...
```

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir).
- [ ] **Step 2: Implement** the parse grammar (both runners), server + twin dispatch, the GUI JS
  in both `gui.js` files, four-locale i18n + `draw.html`, `results_html.py` + `viz/tikz.py`.
- [ ] **Step 3: Gate layout/i18n first** —
  `... -m pytest tests/webapp/test_layout_picker.py tests/webapp/test_js_parses.py -q`.
- [ ] **Step 4: Add the golden** (`koszul_kx3`); note it in `test_runner_delegation.py`'s
  change-log. Run the delegation test BEFORE adding to confirm existing entries byte-identical.
- [ ] **Step 5: Run the gates** —
  `... -m pytest tests/webapp/test_koszul_p77.py tests/webapp/test_runner_delegation.py tests/gui/test_koszul_runner_twin.py tests/hpc -q`.
- [ ] **Step 6: Commit** — `feat(gui,webapp,hpc,trace): koszul compute kind -- generalized Koszulity profile (N-Koszul/K2/almost-Koszul + generation degrees), 4 locales, v1 schema, one golden`

**Adjust to reality (Task 7):** if `estimator.sizing_dim` needs a `koszul` entry, size on
`A.dim` exactly as `ext_algebra` (algebra-level). Confirm the THEMES sentinel block in BOTH
`gui.js` files gets the new kind. The block already carries `"kind": "koszul"` (set in the
builder) so both runners agree.

---

### Task 8: verification page, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`,
  `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick P77)
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`)

- [ ] **Step 1: Verification page.** Add the Plan-77 subsystem row (`modules/nkoszul.py`):
  - `oracle_literature` — `k[x]/x^N` internal degrees = Berger `δ(n)` (`N=2,3,4,5`); the
    Cassidy–Shelton K₂ verdict on `k[x]/x^N` (gens degrees {1,2}); the **BBK preprojective
    classification** (`Π(A₃)`=(2,2), `Π(A₄)`=(3,2), `Π(A₅)`/`Π(D₄)`=(4,2) — engine-verified
    against BBK `(h−2,2)` via `q=e−p`); `Π(A₂)` is the Koszul boundary; `Π(A_{n≥3})` is
    K₂-False (a degree-3 Yoneda generator); `kA₄/J³` is finite-gl.dim 3-Koszul.
  - `oracle_crossengine` — the **two N-Koszul certificates agree** (internal-degree purity ≡
    N-homogeneous ∧ K₂ — Berger's `N≥3` theorem); `generation_degrees ≡ Berger δ(n)`; P77's
    Yoneda generation degrees ≡ Plan-27 `generators_by_degree`; the quadratic case ≡ Plan-27's
    `koszul`.
  - `oracle_selfcert` — the extraction homogeneity (singleton internal-degree per generator);
    the Berger closed form; `q = e − p`; the **P75 seam** (`ext_algebra.koszul` obstruction
    degree == almost-Koszul break degree); K₂ vacuity on semisimple/hereditary; the length-graded
    refusal.
  - the **honest-scope entries**: (a) recognizers require **length-graded** `A` (loud refusal
    off it — internal degrees undefined); (b) **certified window** — K₂/N-Koszul `True` beyond a
    finite `gl.dim` is "holds through degree W", never unconditional; (c) **multi-Koszul is
    connected-graded only** (multi-vertex refuses → K₂; the canonical Yang–Mills examples are
    infinite-dimensional, out of the f.d. engine); (d) the almost-Koszul **period is
    window-observed**, not proven (BBK's `2(h−1)` is cited); (e) the exact **Herscovich
    multi-Koszul Definition is a `# PIN`** (transcribe at implementation) — if kept a pin, a
    named GUI-deferral-ledger entry; (f) **QPA cannot compare** — no Koszul surface (honest
    no-surface probe, fails if QPA ships one).
  Recount the class table (`tests/release/test_oracle_classes.py` — collect, paste LIVE counts,
  re-run to green; mid-merge-train honest).
- [ ] **Step 2: README.** One features line: "generalized Koszulity (Berger N-Koszul,
  Cassidy–Shelton K₂, Brenner–Butler–King (p,q)-almost-Koszul): the generation degrees of
  `Ext•(k,k)` off the shipped minimal resolutions, the 2-N alternation certificate, K₂ with an
  explicit certified window, and the almost-Koszul classifier for the Dynkin preprojectives —
  clickable via `koszul`."
- [ ] **Step 3: Full gate:**
  `... -m pytest tests/modules/test_nkoszul*.py -q` (deep, the touched files),
  `... -m pytest -q -m fast`, `... -m pytest tests/webapp tests/gui -q`,
  `... -m pytest tests/qpa -q -m qpa`, `... -m pytest tests/release tests/citations -q` — all
  green; both runners byte-identical.
- [ ] **Step 4: Commit** — `docs(verification): Plan-77 N-Koszul/K2/almost-Koszul oracle rows + Berger/Cassidy-Shelton/BBK/Herscovich/Chouhy citations + honest scope (length-graded, window, multi-Koszul connected-only, no-QPA) + recounted classes`

---

## Acceptance (Plan-77 definition of done)

1. `generation_degrees(A, i, length)` returns the EXACT internal generation degrees of the
   minimal resolution of `S_i`, validated against Berger `δ(n)` on `k[x]/x^N` (`N=2,3,4,5`) and
   the multi-vertex `kA₄/J³`; loud refusal on a non-length-graded presentation.
   `Algebra.koszul_profile(top)` returns the full profile — float-free, over any exact Domain.
2. **N-Koszul (Berger):** `n_koszul_certificate` reproduces the 2-N alternation `δ(n)` on
   `k[x]/x^N` (`N≥3`) and `kA₄/J³` (the internal-degree pins are exact and window-free). The
   **verdict is window-honest**: `kA₄/J³` (finite `gl.dim`) is a clean `True`; `k[x]/x^N`
   (self-injective, infinite `gl.dim`) is `True`-**through-the-window** (or `None`), never an
   unconditional `True` — the decisive claim is the internal-degree table, anchored by Berger's
   theorem (a `oracle_literature` fact). It **defers the `N=2` case to Plan 27 verbatim** (the
   named overlap, hereditary absorbed), and the two certificates (internal-degree purity ≡
   N-homogeneous ∧ K₂) agree cross-engine.
3. **K₂ (Cassidy–Shelton):** `k2_certificate` decides `E(A)` generated in degrees 1,2 through an
   **explicit certified window** (three-valued: complete when `gl.dim` finite, honest
   inconclusive beyond, decisive False on a degree-≥3 generator). `Π(A_{n≥3})` is K₂-False;
   `k[x]/x^N` is K₂-through-window.
4. **(p,q)-almost-Koszul (BBK):** `almost_koszul_certificate` classifies the Dynkin
   preprojectives via `p = ` top degree, `q = e − p`, reproducing BBK `(h−2, 2)` for
   `Π(A₃)/A₄/A₅/D₄`; `Π(A₂)` is the Koszul boundary. **The P75 seam is pinned**: the
   `ext_algebra.koszul` obstruction degree that P75's GHMS route refuses == the almost-Koszul
   break degree P77 classifies.
5. **Multi-Koszul settled:** connected-graded (local) only; multi-vertex refuses loudly and
   reports K₂ (the transfer); the exact Herscovich Definition transcribed (# PIN resolved at
   implementation) or a named ledger deferral — **no fabricated verdict**.
6. `koszul` clickable end-to-end (GUI canvas → block → generation-degree table + four verdicts →
   report TikZ) in all four locales (en/es/fr/zh), both runners byte-identical via the shared
   builder, ONE golden added with a documented change-log entry, `test_layout_picker.py`
   `ALL_KINDS`/THEMES gate green, schema still v1, canonical keys unchanged.
7. QPA honest-scope probe skips (fails if QPA ever ships any Koszul surface);
   `docs/verification.md` recounted (live numbers) with the honest scope; README line added;
   deep (nkoszul files) + fast + webapp/gui + qpa + release + citations suites green. **Plan 27
   left byte-unchanged.** Honest scope recorded: length-graded requirement; certified window;
   multi-Koszul connected-only + infinite-dim canonical examples out of scope; almost-Koszul
   period window-observed; Herscovich Definition `# PIN`; QPA cannot compare.

---

## Methodology & assumptions

**Approach.** My worktree was stale (v0.3.0 vs `dev` at Plan 74), so I read every input via
`git show dev:<path>` (the metaplan §5 P77 + §§1–4, the R36 record, the P63/P65 style templates)
and confirmed the code surfaces on the `.venv`, which resolves the **main-repo v0.3.0 tree**
(`quiverlab.__file__` = `…/quiverlab/src/quiverlab/__init__.py`, version 0.3.0). I verified that
`modules/ext_algebra.py`, `modules/koszul.py`, `modules/resolution.py` are **byte-identical
between v0.3.0 and `dev`** (`git diff v0.3.0 dev` touches none of them — only
`engine/resolutions_minimal.py` grew, additively, with P52's `_coeff_*` bimodule-HH functions,
which the module-Ext path does not use). **So every live run below is faithful to `dev`.** I
read `ext_algebra.py` and `koszul.py` end-to-end (the Yoneda engine, `generators_by_degree`, the
three-valued `koszul` verdict; `is_quadratic`, `g_quadratic_certificate`, `_is_length_graded`,
`froberg_obstruction`) and `modules/resolution.py` (the Plan-05 minimal resolution — terms +
differentials), and mapped the `ext_algebra` compute-kind wiring in `hpc/spec.py` (parse,
`_dispatch`, `_snippet`) and `docs/gui/runner.py` (`compute_one`, `calls`, `ETA_MODEL`).

**Empirical ground truth (run this authoring).** I ran the shipped `ext_algebra_block` and a
prototype internal-degree extractor on: `k[x]/(x^N)` (`N=2..5`), `kA₄/J³`, and the Dynkin
preprojectives `Π(A₂), Π(A₃), Π(A₄), Π(A₅), Π(D₄)`. **The extractor was validated against
Berger's closed form on all four `k[x]/x^N` BEFORE any preprojective pin was trusted** — it
reproduces `δ(n)` exactly (the Live-verified pins tables are verbatim engine output). The
homological Yoneda-generation degrees are the shipped `generators_by_degree`; the internal
(path-length) degrees I extract from the homogeneous differentials of the shipped minimal
resolution (each `d_n` entry a path of a known length, propagated from `P_0` in degree 0). The
extraction is self-validating: it raises on any inhomogeneous column, and no such column arose
on any length-graded input.

**The BBK reconciliation (the one non-obvious call).** My first reading predicted the
preprojective resolution's linear strand would run through homological degree `p = h−2`; the
engine instead breaks at homological degree **3 for every Dynkin preprojective** (A₃,A₄,A₅,D₄).
I resolved the apparent conflict by re-reading BBK's definition (Bath research portal): **`p` is
the algebra's TOP DEGREE (concentration in degrees 0..p), not the homological length of the
linear strand.** Computing `p = ` top degree and `q = (first-jump internal degree) − p` then
reproduces BBK's `(h−2, 2)` **exactly** on all four cases (Live-verified pins C — e.g. `Π(A₄)`:
top degree 3, first jump to internal 5, `q = 5 − 3 = 2` ⇒ `(3,2)`, `h=5`). This is the recognizer
Plan 77 ships. I traced the extraction per-summand on `Π(A₃)` and `Π(A₄)` to confirm the
single-jump structure is real (not an artifact) before committing the `q = e − p` design.

**Multi-Koszul settlement.** I settled it per the card: Herscovich's multi-Koszul is a
**connected-graded (`A₀ = k`) notion** (§2 hypotheses). The transfer **is pinned NOW** — the
fix-round re-read of the **ar5iv HTML** (my authoring WebFetch of the raw PDF had returned the
compressed stream; ar5iv is the readable route) gave **Proposition 3.30** verbatim: the Yoneda
algebra of a finitely generated multi-Koszul algebra with a finite-dimensional relation space is
generated in degrees 1 and 2, i.e. **multi-Koszul ⟹ K₂**. So the multi-vertex bridge is **K₂**,
no `# PIN` on the transfer. Multi-Koszul is **not** a single Definition environment — it is
constructed in **§3.2** via Tor/Ext-vanishing of the minimal graded **bimodule** resolution, and
I **scope that decision engine OUT** of this plan (it is not a recognizer over the shipped
surfaces): v1 ships K₂ (Prop. 3.30) + the generation-degree table + the connected-only refusal,
`verdict=None`, with the decision recorded as a named backlog deferral — never a guessed verdict.
I flagged that Herscovich's canonical examples (Yang–Mills, super-Yang–Mills) are
**infinite-dimensional**, hence outside quiverlab's f.d. engine — an honest-scope note.

**Reference verification.** I re-verified every citation at authoring: Berger (J. Algebra 239,
2001 — N-Koszul, `N≥3` ⟺ Ext gen degrees 0,1,2); Cassidy–Shelton (Math. Z. 260, 2008,
arXiv:0704.3752 — K₂ = Yoneda gen degrees 1,2, quoted verbatim); BBK (Algebras Represent.
Theory 5, 2002 — `(p,q)` def + preprojective `(h−2,2)`, **not on arXiv**); Herscovich 1305.1678
+ the JPAA 223 (2019) 1054–1072 title corrected to *"Applications of one-point extensions …"*;
Chouhy 1708.02933 (BLMS 2019 — N-Koszul degeneration-stable). I added GMMVZ *D-Koszul algebras*
(JPAA 193, 2004) as the multi-vertex/semisimple-base N-Koszul foundation. The three `# verify`
tags (JPAA/BLMS/JPAA volume-pages) are flagged for the implementer's BibTeX pass.

**Deliberately NOT checked / assumptions.** (1) I did not build the multi-Koszul **decision**
engine (Herscovich §3.2 bimodule Tor/Ext vanishing) — it is honestly **scoped out** with a named
backlog deferral; v1 ships the pinned transfer (K₂, Prop. 3.30) + generation degrees + the
connected-only refusal, `verdict=None` (never guessed). (2) I did not prove the almost-Koszul
`2(h−1)` periodicity from the window — the recognizer records the **window-observed** period and
cites BBK; the classification claim is the `(p,q)` label + the single-break structure (what BBK's
definition checks pointwise), which the engine verifies. (3) I did not deep-run `Π(A₅)`/`Π(D₄)`
past length 5 or `Π(E₆)` (dim 156 — out of a fast probe box); the pins use the depths reached
(the `(4,2)` label needs only the first jump, which lands well inside). (4) I did not build a
degeneration engine for Chouhy's stability — it is cited as context, not computed. (5) I did not
verify `estimator.sizing_dim` needs a `koszul` entry (likely sized generically like
`ext_algebra`) — an adjust-to-reality in Task 7. (6) The GUI/i18n/tikz anchors are from the
v0.3.0 tree; a mid-merge-train shift is handled by the "adjust to reality" hedges (the
`ext_algebra` sibling is the live template to mirror).

**Fix-round re-verification (2026-08-08, all findings VALID).** Every touched number was re-run
on the `.venv`: (H1) `koszul._algebra_graded_matrices(Π(A₅))` / `Π(D₄)` **raise
`AdmissibilityError`** live (`degree_bound=8` < the length-9 tips) — respec'd `_top_degree` to
`max` path length over `A.basis_labels`, verified total and `= h−2` for `Π(A₂..A₆), Π(D₄)` (⇒
`1,2,3,4,5,4`). (M1) `x*x - x` **raises at construction** (ideal ⊄ rad²) — swapped the
non-length-graded pin to `x*x*x - x*x` (verified it BUILDS, `_is_length_graded == False`), scoped
`pytest.raises` to the `generation_degrees` call with a `match="length-graded"`. (W1) restated
the almost-Koszul predicate as the **signature** (linear prefix + first jump to `e = p+q`,
`q ≥ 2` gate), explicitly scoping OUT BBK's full finite-linear-complex + `2(h−1)` condition; added
the **negative pin** `k[x]/x^N` (`p = N−1`, first jump `e = N`, `q = 1` for all `N` — verified
`N=2,3,4` — so refused, it is N-Koszul); reconciled `Π(A₂)` (`rad²=0 ⇒` genuinely Koszul, BBK's
nominal `(1,2)` subsumed). (H2) re-read the **ar5iv HTML**: **Prop. 3.30** (multi-Koszul ⟹ K₂)
pinned; multi-Koszul lives in **§3.2** as a bimodule-Tor construction (no single Definition),
scoped OUT with `verdict=None` + a backlog deferral. (M2) hereditary `relations==[]` ⇒
`n_homogeneous_degree` returns **2** (quadratic-trivially), absorbed by the N=2/Plan-27 branch —
kA₂ pins `n_homogeneous==2`. (M3) added the `q=1` negative almost-Koszul test. (M4) a
multi-element `ℓ_i(n)` set ⇒ N-Koszul `verdict=False` (purity failure) stated. (W2) Acceptance #2
made window-honest. (H3) noted the extractor shares its single method with any re-derivation —
the anchor is **Berger's closed form**, not self-agreement; `_NonPure` demoted to a defensive
assert. (W3) the seam pin demoted to a structural-consistency invariant (one first-failure
degree, two readings). (W4) `koszul` adds no schema field (algebra-level compute string);
Chouhy doi (`10.1112/blms.12277`) and the Herscovich-JPAA title/volume/pages **confirmed**,
residual `# verify` tags narrowed to the un-resolved dois/issue-pages with reasons.

**Why this is correct.** The N-Koszul, K₂, and almost-Koszul definitions are verbatim from
Berger, Cassidy–Shelton, and BBK; every pin is independently reproduced by the shipped Plan-27
Yoneda engine + a Berger-validated internal-degree extractor on the exact algebras; the two
N-Koszul certificates (pure resolution ≡ N-homogeneous ∧ K₂) cross-validate as an
`oracle_crossengine` pin; the BBK `(h−2,2)` classification is reproduced on four preprojectives
via `q = e − p` after the `p = top-degree` reconciliation; the quadratic case defers to Plan 27
(the named overlap, no fork); the multi-Koszul transfer is settled (connected-only → K₂) with the
one un-transcribable piece an explicit `# PIN`; and every window-bounded verdict is three-valued
with an honest "through degree W" so no infinite-gl.dim `True` is ever over-claimed.

## Change log

- **2026-08-08 authoring** — plan drafted against `dev` (Plan 74 tip) read via `git show`;
  code surfaces verified byte-identical v0.3.0↔dev (`ext_algebra`/`koszul`/`resolution`). Scope
  boundary vs Plan 27 fixed (P27 owns quadratic Koszulity; P77 owns N-Koszul/K₂/almost-Koszul/
  multi-Koszul + generation degrees; N=2 defers to P27). Internal-degree extractor prototyped
  and **Berger-validated on `k[x]/x^N` (N=2..5)**; BBK `(p,q)` reconciled (`p`=top degree,
  `q=e−p`) after the homological-break discrepancy — reproduces `(h−2,2)` on Π(A₃/A₄/A₅/D₄).
  Multi-Koszul settled (connected-graded only; multi-vertex→K₂; Herscovich Definition a `# PIN`).
  All Live-verified pins are verbatim engine output. Citations verified (Berger, Cassidy–Shelton,
  BBK, Herscovich ×2, Chouhy, GMMVZ); three `# verify` volume/page tags flagged.
- **2026-08-08 critic fix round (NEEDS WORK → adjudicated all-valid; the BBK correction was
  INDEPENDENTLY CONFIRMED)** — twelve findings applied, DOCUMENT-ONLY, every touched number
  re-run on the `.venv`:
  - **H1 (major):** `_top_degree` respec'd from `koszul._algebra_graded_matrices` (which RAISES
    `AdmissibilityError` on `Π(A₅)`/`Π(D₄)`, `degree_bound=8` < length-9 tips) to `max` path
    length over `A.basis_labels` (total; `= h−2`, verified `Π(A₂..A₆), Π(D₄)`).
  - **M1 (major):** the non-length-graded pin swapped `x*x - x` (raises AT CONSTRUCTION — false
    green) → `x*x*x - x*x` (admissible, `_is_length_graded==False`); `pytest.raises` scoped to
    the call with `match="length-graded"`.
  - **W1 (major):** the almost-Koszul predicate restated as the **signature** (linear prefix +
    single first jump to `e=p+q`, `q≥2` gate); BBK's full finite-linear-complex + `2(h−1)`
    condition **scoped OUT** (stated); negative pin `k[x]/x^N` (`q=1`, refused); `Π(A₂)`
    Koszul-vs-`(1,2)` reconciled.
  - **H2 (major):** ar5iv re-read — **Prop. 3.30** (multi-Koszul ⟹ K₂) pinned NOW; multi-Koszul
    = §3.2 bimodule-Tor construction (no single Definition), the *decision* scoped OUT
    (`verdict=None` + backlog deferral), K₂ + generation degrees are what v1 ships.
  - **M2/M3/M4 (moderate):** hereditary → `n_homogeneous_degree==2` (Plan-27 branch); the `q=1`
    negative almost-Koszul test added; multi-element `ℓ_i(n)` ⇒ N-Koszul `False` stated.
  - **W2/H3/W3/W4 (moderate/minor):** Acceptance #2 window-honest; extractor's single method +
    Berger closed-form anchor noted, `_NonPure` demoted to defensive; seam pin → structural
    consistency; `koszul` adds no schema field; Chouhy doi + Herscovich-JPAA title/vol/pages
    confirmed, residual `# verify` narrowed.
