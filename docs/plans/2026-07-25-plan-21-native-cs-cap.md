# Plan 21 — Native deep-degree CS cap via the Plan-20 diagonal

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** the Hochschild **cap** product on the Chouhy–Solotar resolution computed
NATIVELY, PAST the bar-comparison window, by reusing the SAME comparison-lifted
diagonal `Δ: P → P ⊗_A P` Plan 20 built for the cup. The cap is the follow-up named in
backlog Tier-2 item 1 ("Plan 21: native CS cap via the Plan-20 diagonal"). No new
construction: Δ and its Koszul-signed tensor differential are reused verbatim; the cap
is just the HOMOLOGY-side collapse of the same object. The Gerstenhaber **bracket**
stays transported/window-bounded by design (it needs the CS brace/circle machinery,
out of scope) — this is the last window-bounded operation to gain a native route.

## The mathematics

**The covariance flip is the heart of it.** Exactly as the minimal engine's Plan-16
note (`engine/resolutions_minimal.py`): cohomology collapses `a·w·b` while HOMOLOGY
collapses `b·w·a` on SWAPPED corner tags. The cup reads Δ the cohomology way; the cap
reads the SAME Δ the homology way.

For a double-PELT key `(a_idx, τ.word, mid_idx, ρ.word, c_idx)` of the
`(p, n-p)`-component of `Δ_n(σ)` (`τ.degree = p`, `ρ.degree = n-p`), with

    b_a  = e_{o(σ)}A e_{o(τ)},   b_mid = e_{t(τ)}A e_{o(ρ)},   b_c = e_{t(ρ)}A e_{t(σ)},

a cochain `f ∈ C^p` capped with a chain `z ∈ C_n` (value `x = z(σ) ∈ e_{t(σ)}A e_{o(σ)}`,
its `hom`-corner coordinates) lands in `C_{n-p}`, supported on `ρ`, with value

    (f ∩ z)(ρ)  =  Σ_Δ  coeff · b_c · x · b_a · f(τ) · b_mid            (SIGN-FREE).

`f` eats the degree-`p` **first** factor `τ` (`f(τ) ∈ e_{o(τ)}A e_{t(τ)}`); the
degree-`(n-p)` factor `ρ` survives; and the `x ⊗_{A^e} P_{n-p} → C_{n-p}` collapse is
the same op-twisted `R·x·L` (`L = b_a·f(τ)·b_mid`, `R = b_c`) that
`resolution.py`'s homology differential uses (`resolution.py:187`, `c·w·a`). The result
lands in `e_{t(ρ)}A e_{o(ρ)}` = ρ's `hom` corner and is read against
`corner(o(ρ), t(ρ), "hom")`. No bar object is ever built, so it works at ANY degree.

**The sign convention, and how the gates pin it.** The collapse is **sign-free** — all
Koszul sign already lives inside Δ's ONE `(−1)^p` on `1⊗d_q`; the collapse adds none
(the bar cap is likewise sign-free, `engine/tt_calculus.py:29-33`). This was not
assumed but PINNED by the gates, exactly as Leibniz arbitrated the cup's signs:

- **In-window anchor** (the arbiter): native ≡ transported `cap_of_cs_classes` mod
  boundary on every nonzero HH^p × HH_n class pair with `max(p,n) ≤ window`. The
  non-commutative quantum CI (`yx = 2xy`) DISTINGUISHES the op-twisted `b·w·a` from the
  `a·w·b` cohomology order — the two disagree there at `(p,n)=(1,3)`, and only `b·w·a`
  matches transport. (Empirically confirmed before writing the module.)
- **Unit cap** `1 ∩ z = z`: exact on chains (Δ's `(0,n)`-component is the standard AW
  term `e⊗σ_v⊗e⊗σ⊗e`, so the identity is exact, not merely mod boundary).
- **Cap-Leibniz** `b(f ∩ z) = (−1)^{p+1}(δf ∩ z) + (−1)^p(f ∩ bz)`: exact over GF(5) on
  the full basis grid. A programmatic sign search over `{±1}²` returned exactly
  `s1 = (−1)^{p+1}, s2 = (−1)^p` for every `(p,n)` tested — matching
  `tt_calculus` `test_cap_leibniz`.
- **Module identity** `(z∩f)∩g ~ z∩(f∪g)`: holds with the NATIVE cup for `f∪g`. In
  `cap_of_cs_classes(f,z)` notation this is
  `native_cap(g, native_cap(f, z)) ~ native_cap(native_cup(f,g), z)` mod boundary — the
  CS transcription of `tt_calculus`'s `f ∩ (g ∩ z) = (g ⌣ f) ∩ z`.

All four hold **simultaneously** with the sign-free `b·w·a` collapse; no convention
that passes some and skips others was shipped.

**Degree edges.** `n = p` lands in `C_0` (a genuine degree-0 chain, never refused).
`p > n` (a cap into negative degree) **raises `ValueError`**, matching the bar
convention (`tt_calculus.cap_cochain` / `test_cap_degree_guard`) — refused, not
silently zero.

## What shipped

- `src/quiverlab/resolutions_cs/cap.py` — `native_cap(res, f_vec, p, z_vec, n)`,
  domain-generic, sign-free `b·w·a` collapse of `diagonal(res, n)`; reuses the Plan-20
  `TensorComplex` chain cache; raises `ValueError` for `p > n`.
- `Comparison.cap_of_cs_classes(f, z, engine="auto"|"native"|"transport")` — mirrors the
  cup routing byte-compatibly: transported in-window (`engine="transport"` byte-unchanged,
  keeps its window refusal at any degree), native past the window (`"auto"`), native
  forced at any degree (`"native"`). Reworded `_WINDOW_MSG` and module/bracket docstrings:
  the cap now has a native route; only the bracket stays window-bounded.
- Folded-in Plan-20 final-review backlog: hardened the QCI cup/cap skip-guards (a dims
  regression now FAILS, not skips); tightened the native-cup bridge test to element-wise
  CS-basis equality; a session-scoped `qci_gf5_diag4` Δ_4 fixture (recovers the per-test
  Δ rebuild cost — the cup and cap deep pins share one Δ_4 build); routed `native_cup`'s
  chain resolution through the `TensorComplex` cache; and a constructed refusal test for
  the diagonal's inconsistent-lift `NotImplementedError` scope edge.

## Gates / oracles (`tests/resolutions_cs/test_native_cap.py`, both numba + pure paths)

Unit cap (exact); in-window anchor native ≡ transported (kx2/GF(32003), non-commutative
QCI/GF(5) — the distinguisher); module identity via native cup+cap; exact cap-Leibniz
(kx2/straddle/QCI); degree edges (`n=p → C_0`; `p>n` raises); past-window compute
(transport would raise); bumped-window transport bridge (element-wise identical CS
bases); engine selector; `auto` in-window byte-identical to transport; multi-vertex
(comm_square unit cap — acyclic, `HH_{≥1}=0`, so the degree-0 four-corner cap is the
substantive check; cn_3_2 cyclic-Nakayama `HH_=(3,0,1,1)` cap-Leibniz where native caps
are nonzero); QQ domain-generic smoke (unit cap + cap-Leibniz exact over the QQ domain);
deep QCI past-window pins reusing the shared Δ_4.

Sources cited in test docstrings: Tamarkin–Tsygan / `engine.tt_calculus` conventions
(cup `⌣`, cap module law `f ∩ (g ∩ z) = (g ⌣ f) ∩ z`, unit, degree guard); no invented
theorem numbers.

## Note on `docs/verification.md`

`docs/verification.md` does NOT exist on `main` (it lives on the concurrent
`plan-22-verification-transparency` branch). **At merge time, Plan 21's cap oracles**
(unit cap, in-window anchor, module identity, cap-Leibniz, bumped-window bridge, the
multi-vertex and QQ cases) **must be added to that verification page** alongside the
Plan-20 cup oracles.

## Status

- [x] **DELIVERED** (2026-07-25, branch `plan-21-native-cs-cap`, UNMERGED). The native
  cap (`cap.py`), the `Comparison` routing, the full oracle set, the five folded-in
  Plan-20 review items, deep pins (QCI past-window + multi-vertex + QQ smoke), and docs.
  Merge/push only when Marco asks.
