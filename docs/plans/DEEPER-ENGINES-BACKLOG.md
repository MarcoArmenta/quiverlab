# Deeper-engines backlog (standing work list)

**Protocol (for future sessions):** when Marco says "continue" (deeper engines), open
this file, take the **topmost unchecked item**, and run the repo's standard flow:
research → dated plan doc in `docs/plans/` → branch `plan-NN-<slug>` → TDD with oracle
validation → deep suites green → docs/CLAUDE.md/ROADMAP updates → merge+push only when
Marco asks. Mark the item here (checkbox + plan number + date) as part of the branch.
Add newly discovered debts to the right tier instead of doing them ad hoc.

Assessment date: 2026-07-22 (post Plan 12 + Plan 13). Grounded in the spec
(`docs/specs/2026-07-12-quiverlab-design.md`, §12 non-goals), ROADMAP, and the refusal
markers in `src/` (`NotImplementedError` / "later phase" strings).

## Tier 1 — debts that limit what already ships (do these first, in order)

- [x] **1. Operations on the newly admitted algebras** — DONE, Plan 14
  (`2026-07-23-plan-14-cs-operations.md`, branch `plan-14-cs-operations`). Found and
  fixed a third uniform-zoo latent bug on the way: the closed-form block map was not
  a chain map for any tip of length ≥ 3 (even monomial). Native deep-degree CS cup
  (lifted diagonal) added to Tier 2. Original item:
  Sköldberg homotopy expansion for the CS↔bar comparison maps Φ/Ψ —
  `resolutions_cs/comparison.py` refuses beyond the quadratic/degree-2 window in three
  places ("Skoldberg homotopy expansion (a later phase)", `_WINDOW_MSG`, "native CS cup
  is a later phase") — and/or a **native CS cup product**; add the missing **cap**
  transport wrapper (internals 09: "no CS transport wrapper yet"). Payoff: Plan 12
  gave deep *dimensions* for every admissible presentation; this gives the full
  calculus (cup/cap/bracket) there too. The paper's contracting homotopy S_n is
  already transcribed in CS §4 (arXiv:1406.2300 TeX ~880–905).
- [x] **2. Corner-mode checkpoint format for `deepen`** — DONE, Plan 15
  (`2026-07-23-plan-15-deepen-corner-checkpoint.md`, branch
  `plan-15-deepen-corner-checkpoint`). Original item: `engine/deepen.py` refuses
  multi-vertex (Plan-13 boundary). Corner data is deterministic from A; persist only
  `cur`/`cur_r`/`rks`/`tags`. Unlocks cluster-scale multi-vertex scans.
- [x] **3. Cohomology from the minimal/corner resolution** — DONE, Plan 16
  (`2026-07-23-plan-16-minimal-cohomology.md`, branch `plan-16-minimal-cohomology`).
  Original item: the minimal engine is homology-only; the Hom-collapse (`e_o A e_t`
  corner, like CS's `side="coh"`) on the same resolution gives deep HH^• for any
  f.d. algebra + a second deep oracle vs CS.
- [x] **4. CS canonicalization** — DONE, Plan 17
  (`2026-07-23-plan-17-cs-canonicalization.md`, branch `plan-17-cs-canonicalization`).
  Original item: reduce the correction-solve solution modulo its nullspace to a
  normal form; flips the 7 `xfail(strict=False)` coefficient pins strict (Plan-04
  stretch item E2). Makes CS differentials byte-reproducible.
- [x] **5. Battery diversity audit** — DONE, Plan 18
  (`2026-07-23-plan-18-zoo-diversity.md`, branch `plan-18-zoo-diversity`).
  Original item: add mixed-length-tip (straddling) and multi-vertex presentations
  to the *standing* zoo (`families/`, batch scans), not just the Plan-12/13 test
  files. Uniform-length zoos hid both 2026-07-22 bugs.
- [x] **6. Field generality of engine-backed invariants** — DONE, Plan 19
  (`2026-07-23-plan-19-field-generality.md`, branch `plan-19-field-generality`).
  Delivered the generic-Domain path for all five (generic (b,B) mixed complex;
  relative-Tor Betti complex ≡ engine rks over every field; socle-criterion
  Frobenius + inner-ν symmetry) AND the reword (the residual refusal —
  structure-constants algebras off GF(p) needing a path-type basis — no longer
  promises a "later phase"). Original item: `complexity`, cyclic homology etc.
  are GF(p)-only; `_require_prime_field`'s hint promises a "later phase that
  generalizes this invariant". Deliver a generic-Domain path or reword.
- [x] **7. Plan 09 — the server tier** — DONE, 2026-07-24, branch `plan-09-web`
  (executed subagent-driven with adversarial critics; plan
  `docs/plans/2026-07-18-plan-09-web.md` + its 2026-07-24 interface-drift
  amendment). The `webapp/` tier: FastAPI app (`create_app`) + a spawn-child
  resource-capped worker fleet over a single SQLite/WAL queue, two-tier compute
  (instant sync / queued jobs) + the email magic-link big-job tier, bilingual
  EN/ES pages with vendored KaTeX, `/literature`, feedback+admin, deploy assets
  (Dockerfile/compose/Caddy/PROVISIONING). All algebra delegated to the library
  (no user-code exec; engine internals never imported). Whole-branch adversarial
  review's four cross-layer majors fixed (async-error genericization, header
  admin token, last-hop XFF trust, lean base wheel). Post-merge backlog captured
  (error-envelope unification, app.js /es dynamic-label localization,
  tier-ordered claim, real-SMTP/TLS/concurrency acceptance gaps). Merged to main 2026-07-25.

## Tier 1a — verification transparency (Marco, 2026-07-25)

- [x] **Verification transparency** — DONE, Plan 22, 2026-07-25, branch
  `plan-22-verification-transparency`
  (`2026-07-25-plan-22-verification-transparency.md`). Display in the repo that
  everything is unit tested and say HOW: both against QPA and against theory from
  the literature (constructed examples the literature/known theorems already
  resolve). Delivered: `docs/verification.md` (the two oracle classes in Marco's
  framing; a subsystem → oracles → test-file table; the marker/bucket scheme with
  audited counts — 1377 tests, fast 504 / deep 868 / qpa 5; the CI matrix; an
  honest-scope section naming where QPA cannot compare and which theory oracle
  covers that ground) wired into the mkdocs nav after "Under the hood"; a concise
  "How quiverlab is verified" README section; the **standing rule** that every
  future plan adds its new oracles to the verification page as part of acceptance
  (mirrored in ROADMAP). Audit-first, no `src/` change; no sentence written that
  was not traced to an actual test file.
- [x] **Literature citations for oracle tests** — DONE, Plan 22 follow-up (Marco,
  2026-07-25), branch `plan-22-verification-transparency`. "For unit testing
  against literature we should cite the literature we use to test it." Every
  literature/theory pin on `docs/verification.md` now names its precise source
  (author, year, venue; a theorem/example number only where actually recorded —
  Loday *Cyclic Homology* Thm 2.1.5 is the one such number; Happel's is left as
  "theorem" because no number is in the repo), reusing the `references.bib`
  registry keys (`bar`, `bardzell`, `chouhy_solotar`, `happel_question`,
  `quantum_ci`, `qci_hh_oracle`, `tensor_product`, `cyclic`, `minimal_resolution`,
  `assem_book`/`nakayama`, `han_conjecture`, `qpa`) where they exist and naming
  Loday / "BACH" / Assem–Skowroński (1987) at test-verified precision where they
  do not. Added a full **References** section to the page (consistent with
  `quiverlab.bibliography(...)`); added the citation obligation to the standing
  rule (page + ROADMAP constraint); and added the missing citations to the
  value-pin test docstrings (`test_tensor.py`, `test_minimal_cohomology.py`,
  `test_ext.py`, `test_engine_validation.py`, `test_battery_bardzell.py` —
  docstring/comment-only, no logic change). No theorem number guessed.

- [x] **TrivialExtension double-quiver presentation** — DONE, Plan 31,
  2026-07-26, branch `plan-31-trivial-extension`
  (`2026-07-26-plan-31-trivial-extension.md`). `TrivialExtension(A)` of a
  presented `A` over QQ/GF(p) now returns a genuine `kQ_T/I_T` (quiver of `A`
  plus one arrow dual to each corner-homogeneous bimodule-socle basis element,
  direction reversed; relations extracted from the ⋉ structure by a length-lex
  kernel enumeration), certified per instance by `dim = 2·dim A` (a loud
  `QuiverlabError` otherwise) and QPA-oracled against the native
  `TrivialExtensionOfQuiverAlgebra`. `is_symmetric`/`is_weakly_symmetric`/
  `is_frobenius`/`is_selfinjective` return True on every `T(A)` via the Plan-29
  trace-form certifier; the four `xfail(strict=False)` fences in
  `tests/invariants/test_symmetric_regression.py` are now real asserts.
  Presentation-less bases keep the unchanged ⋉ structure-constants build (honest
  refusals preserved, doubling as the iso-invariance oracle). Oracles:
  `T(kA_n) ≅ kZ_n/J^{n+1}` (n=2,3,4), `T(k[x]/(x^a)) = k⟨x,y⟩/(x^a,y²,xy−yx)`,
  Cartan `C_T = C_A + C_Aᵀ`, presented ≡ ⋉ iso-invariance, CS ≡ bar on presented
  `T`; new `tests/families/test_trivial_extension_presented.py` (deep) +
  `tests/qpa/test_trivial_extension_qpa.py` (qpa). No Fernández–Platzeck
  citation (metadata not BibTeX-verifiable; per-instance-certified + QPA-oracled).
  Original item (found by Plan 29's is_symmetric fix): give
  `families/trivial_extension.py` a genuine quiver+relations presentation (quiver
  of `T(A)` = Q plus the dual arrows; relations from the ⋉ structure) so the
  path-basis invariants serve it — `is_symmetric(T(A))` then returns True via the
  trace-form certifier (verified: every tested `T(A)` has a nondegenerate
  symmetric trace form); the `xfail(strict=False)` fences in
  `tests/invariants/test_symmetric_regression.py` auto-flip to xpass.
- [x] **Oracle-class test markers (v0.1.0 release gate — Marco, 2026-07-26):**
  DONE, Plan 32, 2026-07-26, branch `plan-32-oracle-markers`
  (`2026-07-26-plan-32-oracle-markers.md`). Three orthogonal markers registered
  (`oracle_literature` 670 / `oracle_crossengine` 396 / `oracle_selfcert` 604;
  the existing `qpa` = 112 IS the fourth class, no new marker); the sweep marks
  every oracle test across engine/resolutions_cs/hochschild/modules/invariants/
  families/batch (module-level `pytestmark`, per-test splits in mixed files),
  contract/infra left UNMARKED; audited class×count table on
  `docs/verification.md` gated by `tests/release/test_oracle_classes.py`
  (badge==page, subprocess collection). The marker sweep is byte-identical
  (fast/deep/qpa still 804/1180/112 excluding the new gate); the gate adds 3
  deep tests -> 1183/2099 (page + README badge bumped). Zero file moves.
  Original item text follows. —
  the four-way test taxonomy as ORTHOGONAL pytest markers, not a file
  reorganization (tests overlap classes — batteries pin a literature value AND
  assert cross-engine agreement in one test; directories stay the runtime-bucket
  dimension): `oracle_qpa` (= the existing `-m qpa`), `oracle_literature`
  (literature/theory value pins), `oracle_crossengine` (independent
  implementations agreeing: CS≡bar≡Bardzell≡minimal, pure/numba parity,
  presented≡⋉ iso-invariance), `oracle_selfcert` (internal certificates:
  d∘d=0, order gates, dimension certificates, byte-reproducibility pins);
  unmarked = contract & infrastructure (refusal honesty, no-floats gate,
  freshness gates, webapp/GUI/release). Deliverables: markers registered in
  `pyproject.toml`, the assignment sweep, an audited oracle-class × count
  table on `docs/verification.md` gated by a collection test (the badge==page
  doctrine), and the forward-going battery file-naming convention recorded.
  Zero file moves; fast/deep/qpa buckets byte-identical. SHIPS ON MAIN BEFORE
  THE v0.1.0 TAG (right after Plan 31 merges; a JOSS reviewer can then run
  each oracle class as a one-liner).
- [x] **Nontrivial literature examples at scale** — DONE, Plan 33, 2026-07-26,
  branch `plan-33-nontrivial-examples`
  (`2026-07-26-plan-33-nontrivial-examples.md`). Marco's feedback: the test/paper
  examples were too small (kA₂/kA₃/single loops); scale them to the
  books-and-literature cases and push the quantum ones to higher degree.
  Delivered as SCALE, not new theorems (Plan 29 already pins the small
  directions). Scale batteries: the generalized quantum CI
  `k⟨x,y⟩/(x^a,y^b,yx−q·xy)` for `(a,b) ∈ {(2,4),(3,4),(4,4),(2,5),(5,5)}` with
  Bergh–Erdmann cohomology `[2,2,1,0,…]` pushed past degree 8 and homology
  `[a+b−1,…]` (`qci_hh_oracle`); the preprojective family Π(A₄/A₅/D₄/D₅) (dims
  20/35/28/60) with dim / `is_selfinjective` / Loewy = h−1 structural pins
  (`preprojective`, `assem_book`); the Bardzell depth showcase kZ₂₀/J¹¹ (dim 220)
  reaching Hochschild degree 300, guarded by a context-managed recursion-limit
  raise around the Bardzell walk (`bardzell`); the symmetric Brauer stars
  kZ₄/J⁹, kZ₅/J¹¹ (`skowronski_yamagata`); Taft Λ₅/Λ₆ HH + cyclic homology
  (`taillefer_taft`); the canonical algebra C(2,2,2,2,2) with `HH² = t−3 = 2`, the
  first ≥2 case (`schremmer_wpl`); the Boolean-lattice B₃ incidence algebra (dim
  27, nerve vanishing `HH^{≥1}=0`, `cibils_incidence`); presented trivial
  extensions T(kD₄)/T(kA₅)/T(kA₆) (`cmrs_split`); the wild m-Kronecker
  `HH¹ = m²−1` with Coxeter `t²−(m²−2)t+1` (`happel_question`,
  `lenzing_delapena_spectral`); the exterior algebras Λ(k³)/Λ(k⁴) Koszul via
  `g_quadratic_certificate` (`priddy`, `froberg_koszul`). Builders (C2 src):
  `QuantumCI(q, a, b)` generalization (`QuantumCI(q)` byte-identical), the
  preprojective auto degree-bound table, and the Bardzell recursion-limit guard.
  Papers (C3): the JSC manuscript's worked-examples section rebuilt around the
  research top-10 (Π(D₅) dim 60; QCI dim 16 cohomology-dies / homology-persists;
  Bardzell deg-300 at dim 220; Brauer star dim 55; canonical t=5; B₃ incidence
  dim 27; T(kD₄); 3-Kronecker; Taft Λ₅; Λ(k⁴)) with every printed number
  recomputed by replayable scripts in `paper-jsc/computations/`, plus the
  representation-theory-first interior pass; the JOSS Research-impact folds in the
  QCI-(a,b), m-Kronecker, and Bardzell-depth one-liners in band. Honest-scope
  labels binding: preprojective and exterior HH values are cross-engine-only (no
  published table); Π(D₅) HH at scale, Λ(kⁿ≥4) depth, decompose ≳ dim 50,
  dim-30+ non-monomial HH past ~degree 10, and the Π(E₆)/Π(D₆) builds are
  deferred to the SUBMISSION step-4 cluster list.
- [x] **Literature-oracle battery expansion** — DONE, Plan 29, delivered on
  branch 2026-07-25 (merge pre-authorized), branch
  `plan-29-literature-batteries` (`2026-07-25-plan-29-literature-batteries.md`).
  Delivered: the is_symmetric fix (trace-form certifier; the GF(p)
  ν==identity branch was sufficient-not-necessary — silent False on Brauer
  stars, QPA-confirmed) + is_weakly_symmetric + QPA symmetry crosschecks;
  the Coxeter/spectral, identity-oracle, and HH/HC value batteries (79 new
  tests across invariants/engine/resolutions_cs/families, zero
  literature-vs-engine mismatches EXCEPT the documented CRS-2004 Ex. 2.20
  discrepancy where the verified engine value is pinned instead); module
  **Tor_n(M,N)** (Marco request — duality-anchored, QPA-bridged); 17 new
  citation registry keys. Honest-scope entries on the verification page for
  Ex. 2.20, no-native-QPA-Tor, the T(A) refusal, and the RR-2018 deferral.
  Original item text follows. —
  implement the test batteries from
  `docs/plans/2026-07-25-literature-oracles-deep-research.md` (four cluster
  reports — Cibils, Solotar, Happel/Keller/Rickard, de la Peña/Lenzing/Marcos —
  with per-cluster "best 5" rankings, verification-status labels, and BibTeX
  ready for `references.bib`; the consolidated priority list is at the end of
  that doc). Headliners: Nakayama Coxeter polynomials (LMR 2022, all six
  already recomputed under our convention), Dynkin/affine/canonical Coxeter
  tables + Lehmer [2,3,7] spectral pin, the Happel-1997 trace identity
  `tr Φ = −Σ(−1)^i dim HH^i` (sign pinned on A₃), derived-invariance
  orientation pairs (HH*/HH_*/HC_*), Bergh–Erdmann QCI cohomology
  [2,2,1,0,…] for all a,b (char 0 only — GF(p) root-of-unity trap documented),
  triangular-string revival family, rad²=0 char-2 doubling, a-Kronecker
  [1,a²−1,0,…], Taft cyclic homology (first strong HC family), canonical
  HH²=t−3, incidence ≅ nerve, trivial-extension HH¹, and the Plan-27 feeders
  (Ext k[u,v]/(u²) bidegrees, Cassidy quadratic-non-Koszul witnesses,
  (D,A)-stacked Example 1.2, preprojective Koszul/self-injective facts).
  Every battery cites its source in the test docstring and lands on
  `docs/verification.md` per the standing rule; honest-scope flags (paywalled
  primaries, unverified figures, lossy transcriptions) are in the research doc
  and must be respected — nothing flagged gets frozen as a strict pin without
  the primary. **INCLUDES A LIVE, QPA-VERIFIED BUG FIX (do first):**
  `is_symmetric` returns False on multi-vertex symmetric Nakayama kZ_n/J^L
  with n | (L−1) (Brauer stars; QPA `IsSymmetricAlgebra = true`) and on
  `TrivialExtension(A)` (always symmetric; fails because it carries no quiver
  presentation) — fix the ν-inner sweep in `invariants/frobenius.py` for the
  multi-vertex weakly-symmetric case, consider a double-quiver presentation
  for `TrivialExtension`, and add the QPA-crosschecked symmetry regression
  battery (every currently passing `is_symmetric is True` test is
  single-vertex).

## Tier 1b — the module-theoretic surface (Marco, 2026-07-24)

**Vision:** every representation theorist can use this tool — specify a module in
the GUI without writing code, and read off the classical module-level invariants.
Ordering: these items come right after the in-flight native CS cup/cap (Plan 20);
they share one engine (the opposite-algebra functor + the duality D) and should be
planned together even if delivered in slices.

- [x] **AR translates τ / τ⁻** — DONE, Plan 23, 2026-07-25
  (`2026-07-25-plan-23-module-surface.md`, branch `plan-23-module-surface`).
  Native Auslander–Reiten translate of a f.d. right
  module: minimal projective presentation `P₁ → P₀ → M → 0` (exists —
  `modules/resolution.py::projective_cover` / `minimal_resolution`, Plan 05),
  transpose `Tr M = coker(Hom_A(P₀,A) → Hom_A(P₁,A))` as a right `A^op`-module,
  then `τM = D Tr M`, `τ⁻M = Tr D M`. New machinery: the opposite algebra `A^op`
  as a first-class Algebra (reversed quiver, transposed structure constants —
  `modules/opposite.py`), the duality functor `D` on arbitrary Modules (was
  implicit only in `builders.injective`: `I_v = D(Ae_v)`; explicit
  `D(P_v^{op}) ≅ injective(A,v)` tested), and the corner-transpose `Tr`
  (`modules/duality.py`). Surface: `M.tau()`, `M.tau_minus()`,
  `M.dualize()`, `M.transpose()`, `Algebra.opposite()`; honest gates:
  `τ⁻τM ≅ M` for non-projective indecomposables via `modules/hom.py::is_isomorphic`
  (exact invertible-hom certificate; decomposability caveat documented).
- [x] **Injective resolutions + injective dimension** — DONE, Plan 23, 2026-07-25.
  The dual of Plan 05's `ProjectiveResolution` on the same op+D engine
  (`modules/injective.py`): `E(M) = D(projective cover of DM over A^op)`, cosyzygy
  iteration, `injective_resolution(M, length)`, `injective_dimension(M)` =
  `pd_{A^op}(DM)` (int, or `None` = infinite).
- [x] **Left modules alongside right, right as default** — DONE, Plan 24,
  2026-07-25 (`2026-07-25-plan-24-left-modules.md`, branch `plan-24-left-modules`).
  The module surface accepts `side="right" | "left"` (right stays the default and
  byte-unchanged — right-module repr is byte-identical, whole existing suite green).
  Internally a left A-module IS a right A^op-module; the Plan-23 op+D engine makes
  the wrapper a presentation-only tag: constructors (`Algebra.module(...)` +
  S/P/I builders), Hom/End/Ext, projective AND injective resolutions, τ/τ⁻,
  dimensions all read only `(algebra, action)` and are blind to the side — no math
  forked. `D` (now side-aware) exchanges the two sides over the SAME base algebra
  (its classical contravariant form); `Tr` likewise. Two Plan-23 tests that spelled
  `I_v = D(A e_v)` via `A.opposite().projective(v).dualize()` moved to the honest
  side-aware form `A.projective(v, side="left").dualize()` (documented). QPA
  crosschecks left-side τ/τ⁻/inj.dim/resolutions by feeding QPA the opposite
  algebra. The no-code GUI/webapp item below inherits the side picker in its schema.
- [x] **No-code module input (GUI + webapp)** — DONE, Plan 26, 2026-07-25
  (`2026-07-25-plan-26-no-code-modules.md`, branch `plan-26-no-code-modules`).
  The constructor already existed (`modules/module.py`: dimension vector + one
  matrix per arrow, exact entries, loud `RelationError` when the matrices violate
  the relations); now exposed with zero code. Webapp request **schema v2** gains a
  `module` block (`{dims: {v: n_v}, maps: {arrow: [[…]]}, side}` OR
  `{builtin: {kind, vertex, side}}`) + an `ext_target` (the N in Ext), guarded to
  schema 2, canonicalizing through Plan-25 `canonical_key` (side default explicit;
  a non-module request's key is byte-unchanged). Module compute kinds
  (dimension_vector, rad_top_soc from `radtopsoc.py`, Ext from `modules/ext.py`,
  τ/τ⁻, projective/injective resolution + projective/injective dimension) served
  by ALL THREE tiers (instant/queued/big) with references/citations; a
  relation-violating module surfaces as a clean typed 4xx (never a 500).
  Module-aware tier sizing routes big modules off instant like oversized families.
  GUI: a no-code module panel on the Pyodide canvas — per-vertex dimension picker,
  per-arrow exact-entry matrix grid (dims follow source/target live), a right/left
  side toggle, and the S(v)/P(v)/I(v) pick-lists; the client runner
  (`docs/gui/runner.py`) carries the twin dispatch. Versioned schema bump, served
  by BOTH tiers (Pyodide GUI + Plan-09 server).
- [x] **QPA (GAP) as the oracle** — DONE, Plan 23, 2026-07-25. `-m qpa`
  crosschecks (`tests/qpa/test_module_ar_crosscheck.py`,
  `qpa/crosscheck.py::crosscheck_tau`/`crosscheck_proj_resolution`/
  `crosscheck_inj_resolution`/`crosscheck_inj_dimension`) pin τ/τ⁻ (dimension
  vectors + iso class via `IsomorphicModules` on a translated module —
  `modules/qpa_module.py::graded_form`, `qpa/scripts.py::module_decl`),
  projective/injective resolution term dim-vectors (`ProjectiveResolution` /
  `DualOfModule`), and injective dimension (`InjDimensionOfModule`, `false ↔
  None`) across the zoo incl. the Plan-18 multi-vertex `line_abc_cde`. Theory
  oracles (Coxeter transformation, kA_n tables, self-injective/Nakayama,
  inj.dim vs gl.dim) run without the extra.

## Tier 1c — container & HPC surface (Marco, 2026-07-25)

- [x] **Plan 28 — quiverlab as a container: HPC batch tier + offline laptop app**
  — DONE, Plan 28, merged 2026-07-25, branch `plan-28-hpc-container`
  (`2026-07-25-plan-28-hpc-container.md`). Marco's workflow: download an
  Apptainer-ready image; send to a SLURM cluster with extremely simple
  instructions; edit a sample config YAML (or design in the GUI, which prints
  the config); sbatch → `result.json` on disk; download; render a PDF with the
  same container locally. PLUS (Marco, same day): the image doubles as a fully
  OFFLINE laptop app — `quiverlab-hpc gui` serves the Plan-09 webapp locally
  (no internet, zero code), showing memory limits and time estimates, with a
  build-time-seeded Plan-25 result cache of precomputed examples (curated list
  = Marco's open decision; mechanism + placeholder manifest ship). Wheel-level:
  `[hpc]` extra, `quiverlab-hpc` CLI, `quiverlab/hpc/` spec core (runner
  delegation, byte-stable cache keys), `engine/deepen.py` finally wired
  (checkpoint/resume, exit-75 requeue).

## Tier 1d — GUI/report completeness + module decomposition (Marco, 2026-07-25)

- [x] **Plan 30 — Krull–Schmidt decomposition + GUI/report completeness** —
  DONE, delivered on branch 2026-07-25 (merge pre-authorized), branch
  `plan-30-decompose-reporting`
  (`2026-07-25-plan-30-decompose-reporting.md`). Marco's GUI/PDF feedback:
  τ must certify indecomposability (with `decompose()` into indecomposables —
  coverage-program C1 flagship pulled forward, End-idempotent splitting, QPA
  `DecomposeModule` oracle); Ext/Tor accept ANY module both-args in the GUI;
  resolution tables per-term as `P₁²⊕P₃` (drop #summands + dim-vector
  columns); the pdf/tex bundle exhaustive (P_v/I_v description section,
  rendered differential matrices in worked steps, `.tex` downloadable
  everywhere). Merges AFTER Plan 29 (tor-GUI wiring skip-guarded).

- [x] **Plan 34 — homework-grade reports, honest artifacts, Ext/Tor test gaps,
  the auto→CS depth fallback** — DONE, 2026-07-26, branch
  `plan-34-homework-reports` (`2026-07-26-plan-34-homework-reports.md`),
  RELEASE-GATING (Marco: "I don't think it is ready"). /spawn /devils: 4
  workers + 4 paired adversarial critics (verdicts: NEEDS WORK / HOLDS UP /
  NEEDS WORK / REJECTED — all findings reproduced, adjudicated, fixed) + a
  3-agent fix round. Delivered: rad/top/soc as full representations in
  GUI/webapp/CLI (shared `module_blocks`, dim column dropped); Ext/Tor
  batteries over certified non-projective/non-injective/non-simple modules
  (kA₅ interiors, Euler-form arbiter, Tor≡Ext∘D + balance + QPA anchors);
  the PDF+HTML+JSON artifact contract, all cached (HTML complete/no-elision,
  trace.json exact versioned event stream, PDF page-bounded with
  event-computed MaxMatrixCols + resizebox + stated elision; recorder keeps
  full matrices); homework-depth worked steps (definitions, matrices,
  pivot/rank lines, justifications, side-aware, drift-gated); loud renderers
  (ALL_EVENTS gate), segmented Result footers, named self-map labels;
  worker-recorded why-no-PDF; MathML fallback + leak guard. PLUS the
  dispatch amendment: engine="auto" falls back to Chouhy-Solotar exactly
  where bar/fast raised DepthLimitError (presented algebras only, recorded
  in the trace, in-window byte-unchanged; stale "later phases" hint
  rewritten; Pillar-4 pin updated).

## Tier 2 — natural extensions (v1 non-goals worth revisiting, roughly ordered)

- [x] **Native deep-degree CS cup/cap** (added by Plan 14) — **DONE (cup),
  Plan 20, 2026-07-24, branch `plan-20-native-cs-cup`**
  (`2026-07-24-plan-20-native-cs-cup.md`). A comparison-lifted diagonal
  `Δ: P → P ⊗_A P` built degreewise by a per-degree lift-solve on the CS small
  model (`fields.linalg.solve` + `reduce_mod_nullspace`, canonical/byte-reproducible;
  the same loud `NotImplementedError` scope edge as `_d_general`), giving the CUP
  PAST the bar window (`resolutions_cs/diagonal.py`, `cup.py`;
  `Comparison.cup_of_cs_classes(engine="auto"/"native"/"transport")`). Leibniz is the
  sign arbiter; the transported cup is the in-window anchor. The cap is the follow-up
  item below; the bracket stays transported/window-bounded by design.
- [x] **Plan 21: native CS cap via the Plan-20 diagonal** — DONE, Plan 21,
  2026-07-25 (`2026-07-25-plan-21-native-cs-cap.md`, branch `plan-21-native-cs-cap`,
  merged 2026-07-25). The homology-side sign-free `b·w·a` collapse of the SAME lifted diagonal Δ
  (`resolutions_cs/cap.py::native_cap`; `Comparison.cap_of_cs_classes(engine="auto"/
  "native"/"transport")` routes native past-window, transported in-window byte-unchanged).
  `f ∩ z = Σ_Δ coeff · b_c·x·b_a·f(τ)·b_mid` — f eats the degree-p first factor τ, ρ
  survives; the sign convention is arbitrated (not assumed) by the in-window native ≡
  transported anchor (the non-commutative quantum CI distinguishes `b·w·a` from `a·w·b`),
  the exact unit cap, the exact cap-Leibniz `b(f∩z)=(-1)^{p+1}(δf∩z)+(-1)^p(f∩bz)`, and
  the module identity `(z∩f)∩g ~ z∩(f∪g)` via the native cup — all holding simultaneously.
  Degree edges: n=p → C_0; p>n raises. The last window-bounded operation to gain a native
  route (the bracket stays transported by design). Folded-in Plan-20 review backlog all
  done: hardened the QCI skip-guards (a dims regression now FAILs, not skips), tightened
  the native-cup bridge to element-wise CS-basis equality, a session-scoped Δ_4 fixture
  (cup+cap deep pins share one build), routed `native_cup`'s chain resolution through the
  TensorComplex cache, and a constructed refusal test for the diagonal's inconsistent-lift
  `NotImplementedError` scope edge. NOTE: `docs/verification.md` (concurrent
  `plan-22-verification-transparency` branch) must gain Plan 21's cap oracles at merge time.

- [x] **Ext-algebra / Yoneda-ring presentations** — DONE, Plan 27, merged 2026-07-25 (`2026-07-25-plan-27-ext-algebra.md`, branch
  `plan-27-ext-algebra`). The Yoneda algebra `E(A) = Ext•_A(A/J, A/J)` for any
  admissible kQ/I over any exact Domain: Ext dims = graded Betti numbers off the
  Plan-05 minimal resolutions (minimality ⇒ δ≡0, cocycles = summand
  projections), the product by canonical chain-map lifting
  (`reduce_mod_nullspace`, byte-reproducible), degreewise minimal
  generators/relations over R = k^{Q_0}, `YonedaPresentation` with honest
  certification (`as_algebra()` iff gl.dim exact-finite), and the locked
  convention Ext-quiver = Q / `E(A) ≅ (A^!)^op` (arbitrated by worked anchors +
  the QPA ext-quiver pin). Koszulity three-valued: G-quadratic certifier
  (Priddy PBW), non-quadratic/new-generator/Fröberg falsifiers, quadratic dual
  on Q^op (`modules/koszul.py`). Surface: `Algebra.ext_algebra(top)`. Oracles:
  the 7-battery (k[x]/xⁿ char-independent, hereditary round-trip, rad²=0 path
  counts, QCI n+1 = CS chain count, commutative square self-hosting) + monomial
  Anick gate + 43 live QPA tests (`ExtAlgebraGenerators` dims/generator-degrees,
  `ExtOverAlgebra` ext-quiver, `IsQuadraticIdeal`; QPA has NO IsKoszul — honest
  scope on the verification page). Citations added: `priddy`, `froberg_koszul`,
  `polishchuk_positselski`. Deferred: CS tensored-down deep-Ext accelerator,
  native CS Yoneda coproduct, N-Koszul certifier, Ext(M,M) for arbitrary M.
  Original item: generators/relations of `Ext_A(⊕S, ⊕S)` from Plan-05 module
  resolutions + deep CS; Koszulity checks.
- [ ] **HH cohomology ring structure + support varieties**: after Tier-1 item 1,
  finite generation over the even part; support varieties per module.
- [ ] **BV structure** for symmetric/Frobenius algebras: Connes B is ported
  (`engine/cyclic`); add the Tradler-style BV operator + Δ-bracket compatibility.
- [ ] **Periodicity detection with certificates**: general "syzygy ≅ shifted syzygy"
  detection on the corner engine (beyond the two wrapped families), exact isomorphism
  as certificate.
- [ ] **Single-degree HH mode** (Marco, 2026-07-23): `degree=n` on
  `minimal_homology_dims` / `cs_homology_dims` / `deepen`, surfaced as
  `A.hochschild_homology(n, single_degree=True)`. dim HH_n needs only d_n, d_{n+1}
  + two ranks. Honest expectations per engine: Bardzell/monomial ≈ n-fold win
  (closed-form differentials per degree); minimal/CS: the sequential resolution
  build to n+1 is an unavoidable floor — skip the other degrees' collapse+rank
  (20–50% compute) and roll memory to two consecutive differentials (deepen
  already rolls; add a flag skipping non-target finalizations); bar: no win.
  Composes with periodicity certificates (item above): a certified period turns
  deep single-degree into a lookup — implement both together if possible.
- [ ] **Han's-conjecture batch campaigns**: sweep the open zone with the now-correct
  multi-vertex engines through `quiverlab.batch` (Plan 18 opened the scan surface:
  specs carry quiver data, `_analyze_open` serves multi-vertex).
- [ ] **Open-zone scan cohomology** (found by Plan 18): `_analyze_open` still records
  `cx_cohomology = None` ("no cohomology on this engine") — stale since Plan 16's
  `minimal_cohomology_dims`; wire it in and drop the None-side special-casing.
- [ ] **A∞-structure (Kadeishvili) on Ext** — ambitious flagship; CS small models make
  it feasible.
- [ ] **Propagate `degree_bound` into the CS reduction system** (found by
  Plan 33's paper agent): `PreprojectiveAlgebra("D5")` builds fine with its
  auto bound, but `engine="cs"` on the result rebuilds a reduction system via
  `resolutions_cs/build.py::reduction_system_of(A)` WITHOUT the constructor's
  bound, so deep HH on Π(D₅)-scale algebras refuses/stalls where the build
  succeeded. Thread the certified bound (store it on the Algebra at build
  time, read it in `reduction_system_of`) so CS serves what the constructor
  certifies. Unblocks preprojective deep HH (currently cluster-deferred).
- [ ] **Performance**: numba kernels for the Plan-13 corner path (pure Python today);
  GF(p^n) fast-engine acceleration (int64 stack is GF(p)-only). GPU exploration
  (Marco, 2026-07-25): exact GF(p) rank/elimination on GPUs (FFLAS-style
  delayed-reduction in doubles is provably exact for small p, or int32/int64
  kernels via numba-CUDA) — research-tier: big engineering + certification cost,
  uncertain gains on memory-bound sparse ranks; benchmark before committing.
  Until then every surface says honestly: CPU cores + RAM help, GPUs idle.
- [ ] **GUI**: surface deeper engines (CS depth, Betti sequences) in the Pyodide
  landing-page GUI.
- [ ] **Native AR-quiver** (v1 non-goal; `[qpa]` extra covers it today).

## Tier 3 — quiverlab-web (`plan-09-web`) post-merge polish (from the 2026-07-24 whole-branch review)

- [x] **Result cache — never recompute a known example** (Marco, 2026-07-25) —
  DONE, Plan 25, 2026-07-25 (`2026-07-25-plan-25-webapp-result-cache.md`, branch
  `plan-25-webapp-result-cache`). Original item: cache finished results keyed by
  the CANONICALIZED request (family/quiver/module spec + field + invariant +
  parameters + library version — exact results are deterministic, so a key hit is
  a correctness-safe replay). Identical requests are served from the cache
  instantly, including across users. Big-job tier interaction: email verification
  gates the COST of computing, not access to mathematics — so a big-job request
  whose canonical key is already cached is served immediately WITHOUT email
  verification (no token minted, no email stored); only genuinely new big examples
  go through the magic-link flow, and their results enter the cache for everyone
  after. Invalidate on library version bump; nothing user-identifying in the key
  or the cached record; retention/size cap with LRU sweep alongside the existing
  retention sweep. Delivered: `webapp/server/cache.py` canonicalizer + a
  `result_cache` table (new table, not a jobs column) that "pins" its finished job
  against retention; cache-first checks on `/api/compute`, `/api/jobs`, AND
  `/api/jobs/big` (big check is the first statement, before `big_jobs_enabled`, so
  a cached big example is served even with SMTP off, no token/email/pending row);
  worker records on success only; `sweep_cache_once` (version purge + LRU) beside
  `sweep_once`; bilingual "previously computed" note. In-flight dedup documented as
  benign and skipped (would need a jobs-schema change — not "within the existing
  schema"); `cache_put` is idempotent so concurrent identical completions never
  crash.
- [ ] **Unify client error envelopes**: the API mixes `{error_type,message}`
  (compute/feedback/bigjobs), `{detail}` (FastAPI 422 + the big-job 502
  `HTTPException`), and `{message}` (404s); `app.js` reads only the first, so an
  empty-`compute` submit renders "undefined: undefined". Map all to one shape or
  make `app.js` tolerate `{detail}`.
- [ ] **Localize `app.js` dynamic labels** (`Result`, `Reproduce locally`, and the
  polled job status overwriting the server-localized text) via `data-*` so `/es`
  has no English leak in JS-rendered output.
- [ ] **verify page transient-vs-terminal**: `bigjobs` renders `big.link_used` when
  the big queue is full at verify time even though the link is still valid — a
  retry mints a needless pending row + email. Distinguish queue-full (retry) from
  consumed/expired.
- [ ] **Tier-ordered `claim_next`**: shared FIFO can let up to `big_queue_max` 4h
  big jobs sit ahead of an anonymous queued job (instant tier unaffected). Add a
  tier priority or a dedicated anonymous worker.
- [ ] **Acceptance coverage gaps**: real SMTP relay (STARTTLS/auth in `mail.py`),
  TLS/Caddy end-to-end, and multi-worker concurrency (`claim_next` double-claim,
  `requeue_stale_running`, graceful stop) are asserted by code-reading only — the
  smoke ran single-worker with a fake mailer over bare HTTP; the prod HTML
  worked-steps fallback (no TeX in the image) is unit-tested but not in the smoke.
- [ ] **De-flake `test_instant_compute`**: timing-sensitive tier assertion flaked
  once under a full-dir run (passed isolated + rerun); make it deterministic.

## Done (this backlog's history)

- [x] Plan 12 (2026-07-22): straddling ambiguities + `right_decomposition` + CS
  non-quadratic non-monomial scope lift. Merged.
- [x] Plan 13 (2026-07-22): minimal A^e engine multi-vertex (corner-typed projective
  resolution); loud guard for non-path-type bases. Merged.
- [x] Plan 22 (2026-07-25): verification transparency — `docs/verification.md`
  (two oracle classes, subsystem→oracles→tests table, CI matrix, honest scope) +
  README section + the standing "every plan adds its oracles to the verification
  page" rule. Docs+audit only; no `src/` change. Branch
  `plan-22-verification-transparency`, merged 2026-07-25.
