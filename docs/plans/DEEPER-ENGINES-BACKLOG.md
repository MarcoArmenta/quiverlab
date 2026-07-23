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
- [ ] **2. Corner-mode checkpoint format for `deepen`** — `engine/deepen.py` refuses
  multi-vertex (Plan-13 boundary). Corner data is deterministic from A; persist only
  `cur`/`cur_r`/`rks`/`tags`. Unlocks cluster-scale multi-vertex scans.
- [ ] **3. Cohomology from the minimal/corner resolution** — the minimal engine is
  homology-only; the Hom-collapse (`e_o A e_t` corner, like CS's `side="coh"`) on the
  same resolution gives deep HH^• for any f.d. algebra + a second deep oracle vs CS.
- [ ] **4. CS canonicalization** — reduce the correction-solve solution modulo its
  nullspace to a normal form; flips the 7 `xfail(strict=False)` coefficient pins
  strict (Plan-04 stretch item E2). Makes CS differentials byte-reproducible.
- [ ] **5. Battery diversity audit** — add mixed-length-tip (straddling) and
  multi-vertex presentations to the *standing* zoo (`families/`, batch scans), not
  just the Plan-12/13 test files. Uniform-length zoos hid both 2026-07-22 bugs.
- [ ] **6. Field generality of engine-backed invariants** — `complexity`, cyclic
  homology etc. are GF(p)-only; `_require_prime_field`'s hint promises a "later phase
  that generalizes this invariant". Deliver a generic-Domain path or reword.
- [ ] **7. Plan 09 — the server tier** — spec exists
  (`docs/specs/2026-07-18-quiverlab-web-design.md`); citations/trace already carry its
  hooks; the only planned-but-unbuilt tier.

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
- [ ] **Han's-conjecture batch campaigns**: sweep the open zone with the now-correct
  multi-vertex engines through `quiverlab.batch`.
- [ ] **A∞-structure (Kadeishvili) on Ext** — ambitious flagship; CS small models make
  it feasible.
- [ ] **Performance**: numba kernels for the Plan-13 corner path (pure Python today);
  GF(p^n) fast-engine acceleration (int64 stack is GF(p)-only).
- [ ] **GUI**: surface deeper engines (CS depth, Betti sequences) in the Pyodide
  landing-page GUI.
- [ ] **Native AR-quiver** (v1 non-goal; `[qpa]` extra covers it today).

## Done (this backlog's history)

- [x] Plan 12 (2026-07-22): straddling ambiguities + `right_decomposition` + CS
  non-quadratic non-monomial scope lift. Merged.
- [x] Plan 13 (2026-07-22): minimal A^e engine multi-vertex (corner-typed projective
  resolution); loud guard for non-path-type bases. Merged.
