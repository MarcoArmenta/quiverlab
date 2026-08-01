# Plan 35 — The HH product surface: cup, cap, Gerstenhaber bracket, Connes B (v0.1.0 release gate)

**Date:** 2026-08-01 · **Branch:** `plan-35-hh-products` · **Status:** spec approved, awaiting implementation plan
**Release-gating:** the `v0.1.0` tag waits for this plan to merge green. Main already
carries version 0.1.0 (c55687b) untagged.

## Motivation

Marco (2026-08-01): *"I don't see the computation of the products: cup, cap,
Gerstenhaber and Connes differentials. Implement them if they are not, this
should go with the 0.1.0 version."*

The mathematics is already implemented and test-gated — what is missing is
reachability. Audit result:

| Structure | Implementation | Public API | GUI/webapp |
|---|---|---|---|
| Cup | `engine/tt_calculus.py::cup_product_matrix` (GF(p), bar basis); CS native past-window `resolutions_cs/cup.py` + `Comparison.cup_of_cs_classes` (Plan 20) | none | none |
| Cap | `tt_calculus::cap_product_matrix`; CS native `resolutions_cs/cap.py` + `Comparison.cap_of_cs_classes` (Plan 21) | none | none |
| Gerstenhaber bracket | `tt_calculus::gerstenhaber_bracket_matrix` (circle product, descent test-gated); windowed CS transport `Comparison.bracket_of_cs_classes` (Plan 14) | none | none |
| Connes B | `engine/cyclic.py::connes_B_matrix` (GF(p)); `hochschild/cyclic.py::connes_B_matrix` (any exact Domain, Plan 19) | only inside `cyclic_homology` (dims) | none |

`HH^*` is a Gerstenhaber algebra internally; a user can only see dimensions.
This plan is a surfacing layer plus its full no-code treatment — no new
mathematics is invented; every serving path already exists and is oracle-gated.

## 1. Public API (`core/algebra.py`)

Four methods, mirroring the `hochschild_*` dispatch signatures:

```python
A.cup_products(top, engine="auto", max_cells=4_000_000)
    # for every (p, q) with p+q <= top: the table e_i ∪ e_j expressed in the
    # recorded HH basis of degree p+q
A.cap_products(top, engine="auto", max_cells=4_000_000)
    # tables HH^p ⊗ HH_n -> HH_{n-p} for 0 <= p <= n <= top
A.gerstenhaber_brackets(top, engine="auto", max_cells=4_000_000)
    # tables [·,·]: HH^p ⊗ HH^q -> HH^{p+q-1} for p+q-1 <= top
    # (window-bounded past the bar window; the result SAYS so)
    # tables cover pairs p, q >= 1 (the degree-0 insertion action is out of
    # scope, stated in the block)
A.connes_differentials(top, max_cells=4_000_000)
    # induced B: HH_n -> HH_{n+1} matrices + ranks for 0 <= n < top
    # (no engine= — the two paths are a field property, not a user choice)
```

Return objects: frozen `HHProducts` / `HHBrackets` / `ConnesB` holding

* per-degree dims,
* the exact structure-constant tensors (`dict[(p, q)] -> nested Domain-exact
  entries`; ints over GF(p), Domain elements otherwise — **no floats**, gated
  by the standing AST float gate),
* engine + basis provenance (`"bar/GF(p)"` vs `"cs"` + the CS
  reduction-system fingerprint — constants are basis-dependent and the object
  says which basis),
* for the bracket: the certified degree window actually served,
* `.blocks()` — the canonical serialization both runners consume (the same
  single-shape pattern every existing invariant uses).

Degenerate degrees (dim 0) are present-but-empty, never omitted: a zero table
is a statement, not an absence. Naming is plural (`cup_products`) because the
deliverable is the family of tables up to `top`.

## 2. Engine routing and the refusal matrix

`engine="auto"` resolves per request; **one engine serves one call end-to-end**
— a single call's tables never mix bar-basis and CS-basis classes across
degrees.

| Algebra / field | cup & cap | bracket | Connes B |
|---|---|---|---|
| GF(p), `top` within bar window | bar/`tt_calculus` (fast rank engine, bar HH basis) | bar/`tt_calculus` | `engine/cyclic.py`, induced on the bar HH basis |
| Quiver-presented, any exact Domain (incl. GF(p) past window, QQ, CC) | CS native (Plan-20/21 diagonal, CS basis via `cs_hh_basis`) | no route — loud refusal off GF(p) (the transport itself is the GF(p) tt facade, `comparison.py:66`) | generic `(b,B)` (`hochschild/cyclic.py`), induced on HH over the Domain |
| Structure-constants only, off GF(p) | loud refusal (Plan-19 wording: path basis needed) | loud refusal | generic `(b,B)` — works, no quiver needed |

* Explicit `engine="bar"` / `engine="cs"` keep honest walls — error, never a
  silent fallback — the same contract as `hochschild_cohomology` after the
  Plan-34 dispatch amendment, and the auto decision is recorded in the trace.
* `max_cells` guards the bar side exactly as today.
* Cross-engine agreement is NOT re-proven per call: the Plan-20/21 in-window
  batteries gate the class level; this plan adds the same gate at table level
  (§6) and then trusts it.
* Known limitation, stated not hidden: the bracket has no native past-window
  path (Plan 14's design decision — no CS brace machinery). Past the window
  the bracket tables stop at the window edge and the block says exactly that.

## 3. Webapp + GUI

* **Schema v2**: four algebra-level range kinds — `cup:0..n`, `cap:0..n`,
  `bracket:0..n`, `connes_b:0..n` — accepted wherever `hh_cohomology` is.
  No module block required. Cache keys need zero new logic (kind strings flow
  through `canonical_key` via the compute list); every existing request keeps
  its byte-identical key, frozen-golden gated.
* **Runners**: dispatch lands once in `quiverlab.hpc.spec` (shared by the
  webapp runner since Plan 28) and its Pyodide twin `docs/gui/runner.py`,
  shape-identical. Each kind emits one block: dims per degree, equation-style
  tables (`e_1 ∪ e_2 = 3·f_1`), provenance, the bracket's served window, B's
  matrices + ranks. Library refusals surface as the loud typed error → clean
  4xx, never a 500.
* **Estimator**: the four kinds size like `hh_cohomology` at the same `top`
  for tier routing (their cost is HH + multiplication); GUI wait estimates
  gain entries.
* **GUI**: four checkboxes in the compute pick-list; both vendored `gui.js`
  copies stay byte-identical; rendering reuses `matrixGrid` (indexed grids)
  and the `matIsZero` one-line convention for zero tables. New EN/ES i18n
  keys for every user-facing string.

## 4. Reports — full worked-steps treatment (Marco: all of it in 0.1.0)

The four kinds are fully worked-steps traced, Plan-34 conventions throughout.
The report gains a **products chapter**:

* definitional preamble per structure (what the cup product IS; left-to-right
  path conventions restated),
* HH basis dims per degree,
* the chain/cochain-level computation with rendered matrices — full matrices
  in the recorder, stated elision past the width threshold,
* the reduction to the recorded basis with pivot/rank lines,
* the final tables, with **loud drift gates** asserting the narrated dims
  equal the result block's,
* bracket steps show both circle products and the graded commutator; B shows
  the chain-level matrix and the induced-on-HH quotient step.

`spec.run`'s traced-computation selection is extended so a products request
backs the worked-steps bundle (today only HH and module kinds can). Result
blocks flow into the "Computed results" section automatically (`results_html`
renders every block). Citations ship from the existing registry: `cup`/`bracket`
(Gerstenhaber1963), `cyclic` (Connes1985), the Suárez-Álvarez lifting entry
for the bracket transport.

## 5. Curated cache examples — products added (Marco's amendment)

The six stored bundles (`webapp/precomputed/examples/*`) **gain the product
kinds in their `request.json`** and are recomputed:

* Degree rule: each example's product tops equal its existing HH range top
  (e.g. `qci-q2` computes HH to 8 → `cup:0..8`, `cap:0..8`, `bracket:0..8`,
  `connes_b:0..8`). If a deep example's product computation proves infeasible
  at that top (the dim-312 `nakayama-kz24-deep` is the risk), the spec allows
  trimming that example's product tops with an honest note in the manifest —
  never a silent partial result.
* Changing `request.json` **changes the Plan-25 canonical keys** — accepted
  and intended. The seed key (`container/seed_cache.py`) is
  `canonical_key(model_validate(request.json).model_dump())`; sorted-key hashing
  makes dict order irrelevant, but the compute list is ORDER-SENSITIVE. The
  reachability gate (`tests/webapp/test_curated_reachability.py`, §6) pins TWO
  things: (a) the seed key is a stable `validate→dump` fixed point (a live request
  re-travels the same path, so it shares the key), and (b) the Plan-35 PRODUCTS
  sit in exact gui.js relative order right after `hh_homology`, plus every curated
  request is internally consistent with the one curated-seed-convention order
  (`_GUI_ORDER`). **Full GUI-composability of the curated requests is a
  pre-existing NON-GOAL** — they carry kinds the GUI panel does not offer
  (a bare `dimension`) and order the resolution/dimension/decompose block
  differently from gui.js (which pushes `projective_dimension, injective_dimension,
  decompose` before the resolutions). So `_GUI_ORDER` is the seed convention, not a
  byte-for-byte gui.js mirror, and the test is not a GUI-composability proof.
  Honest consequence, unchanged: a user who DOES tick the products (the kinds the
  GUI offers) computes fresh under a different key — seeded hits stay rare, Marco's
  accepted trade — and the seeded bundles exist to ship the worked examples, not to
  be reconstructed box-by-box.
* The bundle regeneration reuses the Plan-34-era gate discipline from the
  2026-07-31 release prep: recompute through the real runner, assert every
  pre-existing block byte-identical, the only additions being the four new
  product blocks; `tikz.tex` byte-identical; manifest comment updated.
* Seeding (`container/seed_cache.py`) needs no change — it keys fresh from
  `request.json` at build time.

## 6. Testing, oracles, verification (Plan-32 classes)

* **`oracle_selfcert`**: cup unit `1∪f = f`; graded commutativity
  `f∪g = (−1)^{pq} g∪f`; associativity on in-window triples; bracket graded
  antisymmetry + graded Jacobi + cup-Leibniz (Gerstenhaber compatibility);
  cap module identity `(z∩f)∩g = z∩(f∪g)` (the Plan-21 anchor, now at table
  level); `B² = 0`; `Bb + bB = 0`; B-vs-`cyclic_homology` SBI consistency.
* **`oracle_crossengine`**: bar vs CS product tables degreewise on the shared
  window over primes `(32003, 2, 3, 5)`; pure/numba parity rides the existing
  engine gates.
* **`oracle_literature`**: `k[x]/(x²)` — the classic complete Gerstenhaber
  structure; BGMS / Bergh–Erdmann QuantumCI cup structure (`q=1`, GF(2) ties
  into the standing `[4,8,12]` pin); truncated `k[x]/(x^a)` cup ring vs the
  bank closed-form; Connes B ranks vs known HC of these families.
* **`qpa`**: probe live whether QPA 1.37 exposes any comparable product
  (expected: none — then the honest-scope entry names the theory oracle that
  covers the gap, per the Plan-22 standing rule).
* `docs/verification.md`: products row in the subsystem→oracles→tests table;
  audited counts re-counted; `tests/release/test_oracle_classes.py` re-audited.
* Webapp/GUI: one NEW products entry in the runner goldens corpus (existing
  six goldens untouched); cross-runner contract tests in the
  `test_module_blocks` style; a report test reading the product tables back
  out of rendered HTML via the `_result_table`/`_matrix_grid` helpers; and
  the curated-reachability gate — for each curated example, compose the
  request the way the GUI does and assert its canonical key equals the
  seeded row's (the §5 guarantee, kept honest by a test).

## 7. Sequencing and release mechanics

1. Implementation on `plan-35-hh-products`, conventional commits, green at
   every commit; executed subagent-driven with this doc as contract
   (Plans 29–31 pattern).
2. Curated bundle recomputation LAST (after the kinds are final), deep ones
   under nohup + Monitor (the kz24 run is ~90+ min before products).
3. Merge to main → full CI green on the merge commit.
4. Marco tags `v0.1.0` → `release.yml` publishes to PyPI, `container.yml`
   builds GHCR + SIF.
5. Dispatch `desktop.yml` so the desktop app ships the products GUI at 0.1.0.

## 8. Non-goals (0.1.0)

* Class representatives in the public output (`basis=True`) — the structure
  constants pin the algebra structure up to the recorded basis; representative
  rendering/canonicality is its own design problem.
* A native past-window bracket (CS brace machinery) — stays window-bounded by
  design, honestly labeled.
* A new curated products showcase example — the six existing examples now
  carry the products; a dedicated showcase is a later curation decision.
* Per-step trace i18n — the worked-steps report stays EN (existing contract);
  the GUI/webapp block labels ARE bilingual.

## 9. Acceptance criteria

1. `A.cup_products / cap_products / gerstenhaber_brackets / connes_differentials`
   exist, documented, exact, engine-routed per §2, loud on every refusal path.
2. The four compute kinds serve on all three tiers + Pyodide GUI with
   rendered tables, citations, EN/ES labels.
3. The worked-steps report carries the products chapter with drift gates.
4. The six curated bundles carry product blocks; pre-existing blocks
   byte-identical; seeded replay hits by construction.
5. All §6 batteries green; verification page updated with audited counts.
6. Whole suite green on the merge commit; no existing golden, cache key, or
   frozen pin moved except the documented runner-goldens addition and the
   curated `request.json`/`result.json` regeneration.
