# Plan 14: CS operations beyond the quadratic window — homotopy-lifted Φ + cap transport

**Goal:** lift the three `NotImplementedError` boundaries in
`resolutions_cs/comparison.py` (the "Sköldberg homotopy expansion (a later phase)" /
window refusals) so cup/bracket transport works for **every admissible presentation**
inside the bar window — in particular the algebras Plan 12 just admitted — and add the
missing **cap** transport wrapper (internals 09: "no CS transport wrapper yet").
Backlog: Tier-1 item 1 of `docs/plans/DEEPER-ENGINES-BACKLOG.md`.

## What exists (verified 2026-07-23)

`Comparison` already has: Φ# matrices from `_expansion`, cochain transport both ways,
`cup_of_cs_classes` / `bracket_of_cs_classes` / `transport_then_bar_cup` via
`engine.tt_calculus`, and two executable gates (`assert_chain_map`,
`assert_transport_roundtrip_identity`). The ONLY gaps:

1. `_expansion(n, σ)` is a chain map only for (a) monomial presentations (pure block
   map, any n) and (b) any presentation at n ≤ 2 (block map − β(tip) correction).
   For any non-monomial presentation at n ≥ 3 (including the *quadratic* QCI!) it
   raises. So all transport is capped at degree 2 off the monomial case.
2. No homology-side comparison ⇒ no cap transport.

## Mathematical design

**Φ by comparison-theorem lifting (Phase A).** The normalized bar resolution has the
canonical contracting homotopy (extra degeneracy, right-A-linear)
`h(a₀ ⊗ ā ⊗ c) = 1 ⊗ [a₀] ⊗ ā ⊗ c` (`[a₀]` = class of a₀ in Ā = A/k·1; kills the unit
component). Given the chain map Φ in degrees < n, define on generators
`Φ_n(1⊗σ⊗1) = h( Φ_{n-1}( d_n^{CS}(1⊗σ⊗1) ) )` and extend A^e-linearly — a chain map
by the standard telescope (`dΦ_n = dhΦ_{n-1}d = (id − hd)Φ_{n-1}d = Φ_{n-1}d` since
`dΦ_{n-1}d = Φ_{n-2}dd = 0`). Because h always outputs left outer factor 1, every
Φ-image term has the shape `1 ⊗ w ⊗ c`:

- **representation:** `{interior word (tuple of reduced unit-adapted indices) : c-vector
  (length-m int64, unit-adapted coords)}` — the scalar coefficient is absorbed into the
  c-vector. The legacy `_expansion` output `{w: coeff}` embeds as `c = coeff·1`.
- **recursion step:** for each CS term `(coeff, a_word, τ, c_word)` of `d_terms(n, σ)`:
  take the cached expansion of τ, right-multiply every c-vector by `c_word` (A-mult in
  unit-adapted coords), left-multiply by `a_word` and apply h — i.e. reduce `a_word` to
  an A-vector, drop its unit component into nothing (h(1⊗…)=0 telescopes) — no: the
  unit component of a contributes `h(1⊗w⊗c)`, which is 0 in the NORMALIZED bar only
  when the first interior leg is degenerate; h(1⊗w⊗c) = 1⊗[1]⊗w⊗c = 0 ✓ (the [1] leg
  is degenerate). Each reduced component r of a with coefficient λ prepends the leg:
  `w ↦ (r,)+w`, coefficient λ·coeff.
- **signs:** fixed empirically against the executable gate `assert_chain_map` (the
  matrix identity δ_cs Φ# = Φ# δ_bar is the ground truth; if the homotopy needs a
  (−1)^n the gate pins it). Record the final convention in the code comment.
- **scope switch (zero blast radius):** keep the closed forms where they are proven
  and pinned — monomial (any n) and any presentation at n ≤ 2; recurse only for
  non-monomial n ≥ 3, which today RAISES (so no existing pin can change).
- **window unchanged:** `_default_window` logic untouched; the recursion only fills
  degrees the window already allows.

**Cap transport (Phase B).** The covariant collapse `A ⊗_{A^e} Φ` gives the
homology-side comparison directly: for a CS chain-basis element `(σ, x)` (loop x in
`e_t A e_o`), `PhiHom_n` sends it to `Σ (w_i ; c_i·x)` in the engine bar chain basis
`cn_basis` (homology collapse convention b·x·a with a = 1). Then:

- `transport_cycle_cs_to_bar(z, n)` = matrix-apply `PhiHom(n)` (direct, covariant);
- `transport_class_bar_to_cs_hom(z̄, n)` = solve `z̄ ≡ PhiHom(cycles) mod bar
  boundaries` (mirror of `transport_cocycle_cs_to_bar`'s solve, homology side);
- `cap_of_cs_classes(z, f)` (z a CS homology class deg n, f a CS cohomology class
  deg m): transport both to bar, `tt_calculus.cap_cochain`, pull the degree-(n−m)
  chain class back. Window check: max(n, m) ≤ window (bar chains to degree n needed).
- gates: unit-cap identity (`z ∩ 1_{HH^0} = z`) and the module identity
  `(z ∩ f) ∩ g ~ z ∩ (f ∪ g)` inside the window.

## Validation battery (tests/resolutions_cs/test_operations_deep.py)

On the Plan-12 algebras — cubic-tail `k⟨x,y⟩/(x², y², xyx−yxy)`, QCI(3,2)
`k⟨x,y⟩/(x³, y², yx−2xy)`, straddle-monomial `k⟨x,y⟩/(xx, yy, xyx)` — over GF(32003)
and GF(2)/GF(3) where feasible:

1. `assert_chain_map(upto=3)` passes (currently raises NotImplementedError).
2. `assert_transport_roundtrip_identity(upto=3)`.
3. Cup: graded commutativity of transported cup on HH classes
   (`u∪v ~ (−1)^{|u||v|} v∪u`), agreement of `cup_of_cs_classes` with
   `transport_then_bar_cup` (both already implemented — cross-gate).
4. Bracket: `[u, u] ~ 0` for odd-degree u; Leibniz spot-check where dims permit.
5. Cap: unit-cap identity and `(z ∩ f) ∩ g ~ z ∩ (f ∪ g)`.
6. Regression: existing `test_comparison.py` unchanged (monomial + quadratic ≤ 2
   closed forms untouched); `test_operation_window_boundary` still raises OUTSIDE
   the window (the window refusal survives; only the presentation refusal dies).

## Execution log

- [ ] Phase A: recursive Φ expansion (`_phi_expansion_general`), general Phi() pairing
      (c-vector slot), scope switch, chain-map gates green on the three algebras.
- [ ] Phase B: `PhiHom`, homology-side transports, `cap_of_cs_classes`, cap gates.
- [ ] Battery + docs (internals 09 operations paragraph, CLAUDE.md, ROADMAP, backlog
      checkbox for Tier-1 item 1).

Branch: `plan-14-cs-operations`. Stretch (NOT this plan): native CS cup via a
homotopy-constructed diagonal `P → P ⊗_A P` (deep-degree operations past the bar
window) — add to the backlog Tier 2 when this plan lands.
