# Plan 20 — Native deep-degree CS cup via a lifted diagonal Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** the Hochschild cup product on the Chouhy–Solotar resolution computed
NATIVELY — a diagonal approximation `Δ: P → P ⊗_A P` built degreewise on the CS
small model — so `Comparison` cup works PAST the bar-comparison window (today a
degree-n cup needs bar cochains to degree 2n+1 and `_check_window` refuses;
backlog Tier-2 item 1, filed by Plan 14). Cap via the same Δ is **Plan 21**;
the Gerstenhaber **bracket stays transported/window-bounded** (it needs the CS
brace/circle machinery — out of scope, say so in docs, it is not a gap).

**Architecture (research-verified 2026-07-24; construction feasibility, the
missing-homotopy question, and the oracle inventory were confirmed against the
actual code by a dedicated research pass):**
(1) `P_p ⊗_A P_q = ⊕_{τ∈S_p, ρ∈S_q} A e_{o(τ)} ⊗ (e_{t(τ)}Ae_{o(ρ)}) ⊗ e_{t(ρ)}A`
is again a projective A^e-resolution of A, so Δ exists and is unique up to
homotopy. Represent its elements as a **double-PELT** dict
`{(a_idx, τ.word, mid_idx, ρ.word, c_idx): coeff}` extending
`resolutions_cs/pelt.py`'s single-PELT; the tensor differential is
`d_p ⊗ 1 + (−1)^p 1 ⊗ d_q` (Koszul sign on the second summand), applied via
`d_terms` + `ar.mul` gluing exactly like `pelt.apply_lower`.
(2) **No contracting homotopy is needed** (the CS resolution exposes none, and
the Ψ-synthesis routes are blocked by the unbuilt chain-level Ψ): build Δ
degreewise by solving the lifting equation
`d^{P⊗P}_n · Δ_n(σ) = Δ_{n−1}(d_n σ)` per generator σ∈S_n over the domain —
the exact structural clone of `_d_general`'s correction solve
(`resolution.py:114-147`): `fields.linalg.solve` + `reduce_mod_nullspace`
(canonical free-variables-zero representative, Plan-17 pattern, so Δ is
byte-reproducible). Base case `Δ_0(σ_v) = e_v ⊗ σ_v ⊗ e_v ⊗ σ_v ⊗ e_v`.
An inconsistent solve = the identical scope edge `_d_general` already flags —
raise the same loud `NotImplementedError` ("higher CS homotopy"), never silent.
(3) **Cup**: for CS cochains f∈C^p, g∈C^q (coordinates over `_basis(·,"coh")`),
`(f∪g)(σ) = Σ_Δ coeff · ar.mul(b_a, ar.mul(f(τ), ar.mul(b_mid, ar.mul(g(ρ), b_c))))`
— the a·w·b cohomology collapse (`resolution.py:201`) generalized to two
interior values; sign-free, matching the engine's bar cup convention
(`tt_calculus.py:20-25`). Everything stays CS-sized (|S_p|·|S_q|·corners);
no bar object is ever built.

**Tech Stack:** `resolutions_cs/` (AArith corners/mul, SSequence, d_terms,
pelt), `fields.linalg` (solve, reduce_mod_nullspace, nullspace), the existing
`Comparison` transport as the in-window anchor oracle, GF(p) contexts
(Comparison is GF(p)-gated; the diagonal itself is Domain-generic — build it
domain-generic, test it under GF(p) via Comparison and directly over QQ on one
small algebra). Tests in `tests/resolutions_cs/` → deep bucket.

## Global Constraints

- Exact only; no floats in `src/` (AST gate). All arithmetic through the
  algebra's Domain / `fields.linalg`.
- **Canonicalization is mandatory:** every lift solve is reduced with
  `reduce_mod_nullspace` — Δ, and hence native cup vectors, must be
  byte-reproducible run-to-run (Plan-17 law).
- **Loud scope edge, never silent:** an inconsistent lift solve raises
  `NotImplementedError` with the same "higher CS homotopy correction …
  spec §6 risk register" framing as `resolution.py:131-135`. No fallback.
- **The transported path stays intact:** `cup_of_cs_classes` /
  `transport_then_bar_cup` / `bracket_of_cs_classes` / `_check_window` keep
  their current semantics for bracket; cup routing changes are additive
  (Task 4). In-window, native and transported must agree mod coboundary —
  that agreement is a permanent test, not a migration step.
- Signs: Koszul sign lives in the tensor differential's second summand ONLY;
  the cup evaluation is sign-free (bar convention). The Leibniz gate
  (Task 3) is the arbiter — if it fails, the sign placement is wrong; fix
  the code, never weaken the gate.
- Oracles live, never hardcoded: transported cup (in-window anchor),
  `(d^{P⊗P})²=0`, chain-map identity, Leibniz, graded-commutativity /
  associativity mod coboundary (`same_cohomology_class`), and one
  bumped-`max_cells` longer-bar bridge just past the default window.
- Conventional commits; green at every commit; merge/push only when Marco
  asks. Branch `plan-20-native-cs-cup` stacks on `plan-09-web`
  (merge train 18 → 19 → 09 → 20).

---

### Task 1: Double-PELT + tensor differential (`resolutions_cs/diagonal.py`)

**Files:**
- Create: `src/quiverlab/resolutions_cs/diagonal.py`
- Test: `tests/resolutions_cs/test_diagonal.py` (new)

**Interfaces:**
- Consumes: `ChouhySolotarResolution` (`_res.ss.S(n)`, `d_terms`, `ar.corner`,
  `ar.mul`, `ar.path_vec`, `to_int`, `_chain`), `pelt.py` conventions.
- Produces: `class TensorComplex` (or module functions) with
  `basis(p, q, o, t)` → the (o,t)-corner basis of `P_p⊗_A P_q` as
  `(a_idx, τ, mid_idx, ρ, c_idx)` tuples;
  `tensor_matrix(n, o, t)` → matrix of `d^{P⊗P}_n` on the (o,t)-corner of
  `⊕_{p+q=n}`; `apply_tensor_d(n, dpelt)` → double-PELT of the differential of
  a double-PELT element. Task 2 consumes all three.

- [ ] **Step 1: failing test** — `(d^{P⊗P})² = 0` as an exact matrix identity
  on the zoo `["x*x"]` (k[x]/x², quadratic), `["x*x","y*y","x*y*x"]`
  (straddle monomial), `["x*x","y*y","y*x-2*x*y"]` (quantum CI, non-monomial
  correction path) over GF(5), degrees n = 1..4: `tensor_matrix(n−1)·
  tensor_matrix(n) == 0` for every corner. Also a Koszul-sign unit test:
  for p=q=1, the (1,1)-block of `d^{P⊗P}_2` restricted to `1⊗d` carries the
  −1 (compare one hand-computed entry on k[x]/x³).
- [ ] **Step 2: run, verify FAIL** (module missing).
- [ ] **Step 3: implement** — extend the pelt idiom: gluing a `d_terms`
  output `(coeff, a', τ', c')` into slot 1 multiplies `a'` into `a_idx` from
  the right and `c'` into `mid_idx` from the left; into slot 2 multiplies
  `a'` into `mid_idx` from the right and `c'` into `c_idx` from the left,
  with the extra `(−1)^p`. Expand A-coefficients over corner bases exactly
  as `terms_to_pelt` does. Cache per (n, corner).
- [ ] **Step 4: run, verify PASS** (deep marker: run the file directly).
- [ ] **Step 5: commit** `feat(cs): tensor complex P (x)_A P with Koszul-signed differential (Plan 20)`.

### Task 2: The lifted diagonal Δ

**Files:**
- Modify: `src/quiverlab/resolutions_cs/diagonal.py`
- Test: `tests/resolutions_cs/test_diagonal.py` (append)

**Interfaces:**
- Produces: `diagonal(res, n) -> {σ.word: double-PELT}` for σ∈S_n (cached,
  recursive); raises the loud `NotImplementedError` on an inconsistent solve.
  Task 3 consumes it.

- [ ] **Step 1: failing tests** — (a) chain-map identity
  `apply_tensor_d(n, Δ_n(σ)) == Δ_{n−1}(d_n σ)` for every σ∈S_n, n = 1..4,
  on the three Task-1 algebras (exact dict equality after canonicalization);
  (b) byte-reproducibility: two fresh resolutions give identical Δ dicts;
  (c) `Δ_0` base shape.
- [ ] **Step 2: run, verify FAIL.**
- [ ] **Step 3: implement** — per σ: RHS = glue `d_terms(n,σ)` into cached
  `Δ_{n−1}` (double-PELT); assemble the corner-restricted linear system from
  `tensor_matrix(n, o(σ), t(σ))`; `solve` + `reduce_mod_nullspace`; loud
  `NotImplementedError` on None, mirroring `resolution.py:131-135` verbatim
  in tone.
- [ ] **Step 4: run, verify PASS.**
- [ ] **Step 5: commit** `feat(cs): comparison-lifted diagonal on the CS resolution -- per-degree lift-solve, canonical (Plan 20)`.

### Task 3: Native cup + the algebraic gates

**Files:**
- Create: `src/quiverlab/resolutions_cs/cup.py`
- Test: `tests/resolutions_cs/test_native_cup.py` (new)

**Interfaces:**
- Consumes: `diagonal`, `ChouhySolotarResolution._basis(n, "coh")`, AArith.
- Produces: `native_cup(res, f_vec, p, g_vec, q) -> vec` (coordinates over
  `_basis(p+q, "coh")`), pure cochain-level; `native_cup_matrix(res, p, q)`
  optional. Task 4 wires it into `Comparison`.

- [ ] **Step 1: failing tests** —
  (a) **Leibniz** (the sign arbiter): for basis cochains f∈C^p, g∈C^q on the
  three zoo algebras, p,q ≤ 2:
  `dcs(f∪g) == (dcs f)∪g + (−1)^p f∪(dcs g)` exactly (matrix identity over
  GF(5), using `res.matrix(·,"coh")`);
  (b) unit: the augmentation-degree-0 identity cochain is a two-sided unit;
  (c) **in-window anchor**: via `Comparison` on `_kx2_gf()`-style fixtures,
  `native_cup` of two HH-class representatives is cohomologous
  (`same_cohomology_class`) to `cup_of_cs_classes` for all nonzero
  (p,q) with p+q ≤ window — the permanent overlap oracle;
  (d) graded-commutativity and associativity **mod coboundary** on HH
  representatives (pattern: `test_operations_deep.py:45-59`).
- [ ] **Step 2: run, verify FAIL.**
- [ ] **Step 3: implement** the collapse formula from the Architecture block
  (two interior values glued by `b_mid`; read off against
  `corner(o(σ), t(σ), "coh")`).
- [ ] **Step 4: run, verify PASS.**
- [ ] **Step 5: commit** `feat(cs): native cup on the CS resolution -- Leibniz-gated, anchored to transport in-window (Plan 20)`.

### Task 4: Past-window delivery + `Comparison` wiring

**Files:**
- Modify: `src/quiverlab/resolutions_cs/comparison.py` (cup routing only)
- Test: `tests/resolutions_cs/test_native_cup.py` (append),
  `tests/resolutions_cs/test_comparison.py` (adjust ONE window test)

**Interfaces:**
- Produces: `Comparison.cup_of_cs_classes(u, v, engine="auto")`:
  `"auto"` = transported in-window (unchanged behavior + its cross-checks),
  native past-window (no more `NotImplementedError` for cup);
  `"native"`/`"transport"` force a route (transport keeps its window refusal).
  `bracket_of_cs_classes` untouched (still window-bounded; docstring points
  at Plan-20's scope note). `_WINDOW_MSG` reworded honestly: bracket-only.

- [ ] **Step 1: failing tests** — (a) a cup at degree window+1 on
  `truncated_polynomial(2, GF(5))` COMPUTES natively (today raises) and its
  class is nonzero exactly where the closed-form HH ring of k[x]/x² over
  GF(5)≠char-2 says (α∪α at even degrees — pin dims, not bytes);
  (b) **bridge**: on a small fixture, pick n = window+1, recompute the
  transported cup with a second `Comparison(max_cells=BUMPED)` whose window
  now covers n, and assert native ≡ that longer transport mod coboundary;
  (c) `test_operation_window_boundary` in `test_comparison.py` flips: cup
  no longer raises outside the window, bracket still does.
- [ ] **Step 2: run, verify FAIL.**
- [ ] **Step 3: implement routing + reword `_WINDOW_MSG`** (no "later phase"
  for cup; bracket's message names the brace machinery as the reason).
- [ ] **Step 4: run the whole `tests/resolutions_cs/` battery.**
- [ ] **Step 5: commit** `feat(cs): cup past the bar window -- native route wired into Comparison (Plan 20)`.

### Task 5: Deep pins, docs, suites

**Files:**
- Test: `tests/resolutions_cs/test_native_cup.py` (append deep pins)
- Modify: `docs/internals/09-chouhy-solotar.md` (new section: the diagonal,
  the lift-solve, why no homotopy is needed, the scope edge, bracket note),
  `docs/plans/DEEPER-ENGINES-BACKLOG.md` (tick Tier-2 item 1 for cup; file
  "Plan 21: native cap" as the follow-up item), `docs/plans/ROADMAP.md`
  (row 20), `CLAUDE.md` (status).

- [ ] **Step 1: deep pins** — native cup on `QCI32 = ["x*x*x","y*y","y*x-2*x*y"]`
  and one Plan-18 multi-vertex record (`comm_square`) at the deepest degree
  the deep-bucket budget allows (target: 2–3 degrees past each fixture's
  default window; assert Leibniz + one nonzero product class per fixture);
  one QQ smoke of the domain-generic diagonal (k[x]/x², Δ over QQ, chain-map
  identity — no Comparison involved).
- [ ] **Step 2: docs edits** (facts per Architecture; the honest scope table:
  cup native ✓, cap Plan 21, bracket transported/window-bounded by design).
- [ ] **Step 3: suites** — full `-m deep` (detached + Monitor per the
  long-suite protocol) and `-m fast`; `DISABLE_MKDOCS_2_WARNING=true
  .venv/bin/mkdocs build --strict` (docs changed).
- [ ] **Step 4: commit** `docs: Plan-20 status -- native CS cup delivered; cap is Plan 21`.

## Validation matrix

1. `(d^{P⊗P})² = 0` + Koszul-sign unit pin (3 zoo algebras, GF(5)).
2. Δ chain-map identity per generator to degree 4; byte-reproducible
   (two fresh builds identical); loud NotImplementedError path reachable
   only via the `_d_general` scope edge (no new silent failure).
3. Leibniz exactly; unit; in-window native ≡ transported mod coboundary
   (the permanent anchor); graded-commutativity + associativity mod
   coboundary on HH representatives.
4. Past-window: computes where transport refuses; bridge vs bumped-window
   transport; boundary test flipped for cup, preserved for bracket.
5. Deep pins on QCI + multi-vertex zoo + QQ domain-generic smoke; full
   deep + fast suites green; strict docs build exit 0.

## Status

- [x] **DELIVERED** (2026-07-24, branch `plan-20-native-cs-cup`, all five tasks;
  UNMERGED). Tasks 1–4: double-PELT + Koszul-signed tensor differential, the
  lifted diagonal (per-degree lift-solve, canonical), the native cup + algebraic
  gates, and the `Comparison` past-window routing — all committed green with the
  full oracle set (`(d^{P⊗P})²=0`, Koszul-sign pin, chain-map identity,
  byte-reproducibility, Leibniz, unit, in-window anchor, graded-commutativity/
  associativity, past-window compute, bumped-window bridge, engine selector).
- Task 5 deep pins + docs + suites. **Finding — the small-window pattern replaced
  the infeasible past-default-window pins.** The plan's "2–3 degrees past each
  fixture's DEFAULT window" is infeasible: default windows are ~9–10 while Δ(4)
  costs ~28s on the QCI and is out of budget on the straddle. Per the controller
  adjudication, the deep pins use the Task-4 SMALL/ZERO-`max_cells` (or explicit
  `window=`) Comparison so "past-window" lands at LOW absolute degree:
  * **QCI** (`x²,y²,yx−2xy`, GF(5), `window=0`): Leibniz at (2,1) and (1,2) — the
    pairs Task 3 budgeted out, needing Δ(4) (~28s) — plus one nonzero product
    class. Scope fact discovered: `HH^•(QCI/GF5) = (2,2,1,0,2,4,…)`, so `HH^3 = 0`;
    the nonzero product class is the off-diagonal `cup(HH¹[0],HH¹[1]) → HH²` (a
    degree-3 native cup still COMPUTES past the window, class forced to zero).
  * **`comm_square`** (Plan-18 multi-vertex zoo, GF(5)): the FIRST genuine
    multi-vertex diagonal (empty-corner branch). Δ exists and chain-maps at every
    built degree (gldim 2 ⟹ `S(3)=∅`); native cup **== transported cup exactly as
    cochains** (`HH²=0` makes the class test vacuous, so exact match + a nonzero
    native cochain are the substantive checks). No empty-corner defect.
  * **QQ smoke**: the domain-generic diagonal on `k[x]/x²` over QQ — chain-map
    identity + ζ-cycle exact, no Comparison (it is GF(p)-gated).
- **Battery time delta:** the three added deep pins run in **~30s** (numba;
  QCI Δ(4) dominates), well inside the ~180s budget.
