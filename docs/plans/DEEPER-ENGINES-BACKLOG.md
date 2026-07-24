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
  tier-ordered claim, real-SMTP/TLS/concurrency acceptance gaps). UNMERGED —
  merge/push only when Marco asks.

## Tier 2 — natural extensions (v1 non-goals worth revisiting, roughly ordered)

- [ ] **Native deep-degree CS cup/cap** (added by Plan 14): a comparison-lifted
  diagonal `P → P ⊗_A P` computed degreewise like Φ, giving cup/cap PAST the bar
  window (the transported operations are window-bounded by construction).

- [ ] **Ext-algebra / Yoneda-ring presentations** (v1 non-goal): generators/relations
  of `Ext_A(⊕S, ⊕S)` from Plan-05 module resolutions + deep CS; Koszulity checks.
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
- [ ] **Performance**: numba kernels for the Plan-13 corner path (pure Python today);
  GF(p^n) fast-engine acceleration (int64 stack is GF(p)-only).
- [ ] **GUI**: surface deeper engines (CS depth, Betti sequences) in the Pyodide
  landing-page GUI.
- [ ] **Native AR-quiver** (v1 non-goal; `[qpa]` extra covers it today).

## Tier 3 — quiverlab-web (`plan-09-web`) post-merge polish (from the 2026-07-24 whole-branch review)

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
