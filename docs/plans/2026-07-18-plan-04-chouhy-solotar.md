# quiverlab Plan 04 — Chouhy–Solotar Resolution + Operation Transport (flagship)

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax. Do every TDD cycle in order: write the failing test, run it, see the stated FAIL, write the implementation, run it, see the stated PASS, commit.

**Goal:** quiverlab gains the first full implementation *anywhere* of the Chouhy–Solotar (CS) projective bimodule resolution for admissible `kQ/I` (arXiv:1406.2300 = J. Algebra 432 (2015) 22–61), specializing exactly to Bardzell's minimal resolution in the monomial case. **What Plan 04 certifies, precisely** (the honest scope — see the certification note in Global Constraints): (1) **deep Hochschild dimensions** — `HH_•`/`HH^•` to large degree — are certified for the two bank closed-form families (`k[x]/(x^a)` and the quantum complete intersection `k⟨x,y⟩/(x²,y²,yx−ξxy)`) under the byte-level bank oracle plus the CS order-gate; (2) for **general admissible `kQ/I`**, dimensions are *computed* via the CS construction and certified **per instance** by three independent gates — `d∘d = 0`, the CS Theorem 4.1 order condition, and degreewise agreement with the bar/minimal/Bardzell engines inside the window those reach; (3) **CS ≡ Bardzell exactly** on every monomial algebra in the zoo; (4) **cup/cap/bracket** of classes are transported through CS↔bar comparison maps and are certified only inside the bar-buildable degree window (dimensions go deep, operations do not). **Scope restriction (binding, RESTRICT option):** the uncollapsed leading map is guaranteed equal to CS's `δ_CS` for **quadratic tips** (`S ⊆ Q_2`; CS Proposition, TeX line 772, where `v_{n-i} = u_i` so the right factorization is the reverse of the left) and — because the tail correction vanishes — for **all monomial presentations** (the collapsed map matches the bar-oracle-validated Bardzell engine at any tip length). A **non-quadratic non-monomial** presentation (a tip of length ≥ 3 together with a nonzero tail) raises `NotImplementedError` at the exact boundary (spec §6 risk register): there `v_n ≠ u_0` and the uncollapsed `δ` feeding the correction solve is not provably `δ_CS`. Lifting this via a `right_decomposition` upgrade (computing the right block factorization so the even `f_n` first term `v_n ⊗ (v_{n-1}⋯v_0) ⊗ 1` is exact for all tips) is an explicit Plan-04-execution **stretch item**. Deliverables: the ambiguity S-sequence, the resolution terms `A ⊗_E kS_n ⊗_E A` over any Plan-01 `Domain`, the CS differential with its order-condition-pinned correction, degreewise-deepening with guards, the CS↔bar comparison maps, the full validation battery, and the internals chapter `docs/internals/08-chouhy-solotar.md`.

**Architecture:** CS lives in a **new top-level package `src/quiverlab/resolutions_cs/`**, engine-adjacent but *not* under `engine/`. Justification (reviewer-verified): `engine/` is the ported hanlab **int64 GF(p) accelerator** (`engine.hh_engine.Algebra`, numpy `int64`, ranks reduced mod p at the end); CS is required by spec §6 to run "over the domain-generic exact arithmetic," i.e. over any Plan-01 `Domain` (ℚ, ℚ(α), GF(p), GF(p^n)), and consumes the Plan-03 Domain-generic `groebner.ReductionSystem` and the Plan-01 `core.Algebra`. Placing it beside the engine keeps the exact/generic contract clean and avoids int64 coefficient blow-up on ℚ. Package modules: `ambiguities.py` (the S-sequence — built by reusing the **ported, validated** `engine.resolutions_bardzell.MonomialPresentation.associated_paths` / `left_decomposition` on the tip monomial algebra `A_S = kQ/⟨S⟩`), `aarith.py` (the `A`-arithmetic layer: paths→normal-form `A`-vectors via `core.Algebra.multiply`, collapse corners), `resolution.py` (`ChouhySolotarResolution`: the leading Bardzell map `δ_n`, the order-condition-pinned correction, the collapsed matrices), `pelt.py` (uncollapsed P-element arithmetic for the correction solve), `homology.py` (HH•/HH• dims + bases via `fields.linalg`, returning `HHTable`), `comparison.py` (CS↔bar comparison maps + class transport), `trace.py` (plain trace dataclasses). A thin `engine_facade.CSResolution(engine.resolutions.Resolution)` re-exposes the same construction over `int64`/GF(p) so CS plugs into the existing resolution dispatch and cross-checks Bardzell/bar on the int64 engine with **one** differential implementation.

**Tech Stack:** Python ≥ 3.10; deps unchanged from Plan 03 (`numpy>=1.21`, `sympy>=1.12`; `[fast]`=numba optional; `[dev]`=pytest). No new hard dependency. Exact arithmetic only.

## Global Constraints

- **Venv is ALWAYS** `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python` (Python 3.12; system python is older). This matches the committed Plan 03 header. Repo root: `/Users/marco/Desktop/HomologicalNetworks/quiverlab`.
- **Every test command** exports `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2` and runs `-m pytest -q` from the repo root. **All test commands run in the FOREGROUND — never backgrounded.** ≤2 parallel processes; no multiprocessing in core.
- **Exact arithmetic only. The float-ban AST gate must stay green** at every commit (`tests/test_no_floats.py`, the Plan-01/03 gate). No float/complex literals, no `float()`, in any `src/` file. Run it in every task's suite.
- **Left-to-right path composition everywhere** (Assem–Simson–Skowroński; the repo standing convention, and the convention Plan 03's `groebner` is written in): `a*b` = "first `a`, then `b`", requiring `target(a)=source(b)`; for an arrow `α: i→j`, `e_i·α·e_j = α`. **CS composes right-to-left** in several places (bar tensors, the `φ_0`/`c=c_n⋯c_1` indexing, the right-ambiguity factorization). Every place this plan transcribes a CS formula it states the translation explicitly (see the Mathematical Preamble "Convention translation", and the per-formula notes).
- **Certification honesty (binding).** "Certified" in this plan means exactly one of: (a) *byte-equal to the bank closed-form oracle* `hanlab/resolutions_cs.py` (the two families, any depth); (b) *equal to a cited literature value*; (c) *equal to the bar/minimal/Bardzell engine* within their feasible window; or (d) *passes the two executable per-degree gates* `assert_dd_zero` and `assert_order_condition` (CS Theorem 4.1 condition (2)). Deep dimensions for **general** admissible `kQ/I` outside the two bank families are labelled "computed (gates d, and c within window; not proved minimal to all degrees)". `HHTable.engine` reads `"Chouhy-Solotar"`; the docstring of every public entry point states which of (a)–(d) back it.
- **Guards fail loudly with the certified range.** Term growth is bounded by a `max_cells`-style guard mirroring `hochschild/bar.py` and `engine/adapter.py`; on hit, `DepthLimitError(msg, hint=…)` states the last certified degree and the resume call.
- **Admissibility gate + NotImplemented boundary.** CS runs only on a certified-admissible `ReductionSystem` (Plan-03 `build_reduction_system` succeeded: `is_confluent` and a finite `irreducibles` basis). Anything where the correction linear system (Task 6) is inconsistent — a non-graded pathology needing the higher CS homotopy correction outside v1's construction — raises `NotImplementedError` **at the exact degree/chain** (spec §6 risk register), never a wrong answer.
- **Trace hooks are plain dataclasses** (`AmbiguityEvent`, `ResolutionTerm`, `DifferentialEvent`, `LiftStep`); **formal rendering (PDF/HTML/text) is Plan 07.** Plan 04 only *populates* them and asserts their claims equal computed values.
- **The bank is READ-ONLY.** `$BANK = /Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/HansConjecture`. CS uses two bank files as **oracles only** — `hanlab/resolutions_cs.py` (closed-form CS differentials; homology only — its `cochain_basis` raises `NotImplementedError`) and the already-ported `engine/resolutions_bardzell.py`. Never write, move, or delete anything under `$BANK`.
- **`reduction_algebra.py` decision — SUPERSEDED, not ported.** Its job (materialize `A` from a reduction system) is Plan-03's `Quiver.algebra()` (general `kQ/I` via `groebner`) + `groebner.ReductionSystem`. Plan 04 consumes those directly (spec §13 "reduction_algebra → 4/9"; Plan-02 deferral "Plans 04/06"): superseded by Plan 03.
- **Conventional commits.** Every commit command below inlines the trailer:
  ```
  Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
  ```
- **Full test suite green at every commit.** Baseline is whatever Plan 03 leaves green; each task states the delta it adds.

---

## Mathematical Preamble (the plan carries the mathematics completely)

The implementer knows Python, not homological algebra. This section is binding; every formula is what the code computes. Citations: Chouhy–Solotar, arXiv:1406.2300v2 = J. Algebra 432 (2015) ("CS"); Bardzell (1997); Gerstenhaber–Schack; Buchweitz–Green–Madsen–Solberg ("BGMS"). Formulas transcribed from the CS TeX source (`§2` Preliminaries, `§3` Ambiguities, `§4` The resolution, `§6` Morphisms in low degrees).

### P0. Setup and the reduction system (CS §2)

`A = kQ/I`, `E = kQ_0 = ∏_v k e_v`. Plan 03 supplies a **reduction system** `R = {(s, f_s)}` satisfying CS's **condition (◆)** (CS Definition, §2): `I = ⟨s − f_s⟩`, every path is reduction-unique, each `f_s` is irreducible. Then (CS Lemma "bases", the Diamond-Lemma restatement) the **irreducible paths `B`** form a `k`-basis of `A`. `S = {s : (s,f_s) ∈ R}` is the set of **tips** (leading words); by CS Remark (§2), `S = Mintip(I)` is the set of minimal non-irreducible paths, and CS assume w.l.o.g. `S ⊆ Q_{≥2}` (admissibility). `β(p)` denotes the normal form of a path `p` (CS Remark: `β = i∘π`).

**The associated monomial algebra** `A_S := kQ/⟨S⟩` (CS §3, opening): its ambiguities depend **only on the tips `S`, never on the tails `f_s`**. This is the structural fact that makes the S-sequence a purely monomial computation on `A_S`, so quiverlab reuses the ported, validated Bardzell machinery for it and CS ≡ Bardzell holds term-for-term in the monomial case.

### P1. The S-sequence (CS §3 "Ambiguities"; = Bardzell `AP^n` = quiverlab `S_n`)

CS index their ambiguity `E`-bimodules `k𝒜_i` for `i ≥ −1` (CS §3): **`𝒜_{-1} = Q_0`** (vertices), **`𝒜_0 = Q_1`** (arrows), **`𝒜_1 = S`** (tips), and for `n ≥ 2`, **`𝒜_n =` the set of `n`-ambiguities**. quiverlab shifts by one so `P_n` is indexed by `𝒜_{n-1}` (bank `resolutions_cs.py` convention: "`P_n` indexed by CS `𝒜_{n-1}`; engine `d_n` = CS `d_{n-1}`"). **Adopted indexing (binding, reviewer-verified):**

```
S_n := 𝒜_{n-1}    for n ≥ 0.
S_0 = Q_0 (vertices),  S_1 = Q_1 (arrows),  S_2 = S (tips),  S_n (n≥3) = (n-1)-ambiguities.
P_n = ⨁_{σ ∈ S_n} A e_{o(σ)} ⊗_k e_{t(σ)} A ;  d_n : P_n → P_{n-1}.
```

**`n`-ambiguity, verbatim (CS §3 Definition).** For `n ≥ 2` and `p ∈ Q_{≥0}`, `p` is a **left `n`-ambiguity** if there exist `u_0 ∈ Q_1` and irreducible paths `u_1,…,u_n` with
```
(i)  p = u_0 u_1 ⋯ u_n     (n+1 blocks),
(ii) for all i, u_i u_{i+1} is reducible but u_i d is irreducible for every proper left divisor d of u_{i+1}.
```
(each consecutive pair `u_i u_{i+1}` is exactly a tip in `S`, minimally). CS Proposition (uniqueness) proves the factorization is unique, and the sets of left and right `n`-ambiguities coincide. CS Proposition (quadratic case, `S ⊆ Q_2`): `𝒜_n = {α_0⋯α_n ∈ Q_{n+1} : α_i ∈ Q_1, α_{i-1}α_i ∈ S}` (every length-2 sub-path is a tip) — the clean case for the quantum CI.

**Implementation (Pillar-1, reviewer-mandated re-base).** `S_n` for `A` = the Bardzell associated paths of the **tip monomial algebra `A_S`**. The ported `engine.resolutions_bardzell.MonomialPresentation` implements CS §3 exactly: `associated_paths(n, maxlen)` returns `AP^n` (with `AP^1 = Q_1`, `AP^2 = S`, `AP^n` the `(n-1)`-ambiguities), and `left_decomposition(p, n)` the unique block decomposition `[u_0,…,u_{n-1}]` (`n` blocks). Since **Bardzell `AP^n` = CS `𝒜_{n-1}` = quiverlab `S_n`**, quiverlab `S_n = pres_S.associated_paths(n, maxlen)` where `pres_S = MonomialPresentation(Q_0, Q_1, relations = leading_words)`. `S_0 = Q_0` handled separately. This reuses code already cross-checked against the bar oracle and makes CS ≡ Bardzell an *identity of the same generator sets*. (My earlier hand-rolled `_extend` was wrong — it stored the stick-out as a block, gave wrong `S_3` blocks, and made `S_4 = ∅` for `k[x]/x²`; this re-base eliminates it.)

**Monomial degeneracy (M2 fix).** For `k[x]/(x^a)` the single tip is `x^a`; each `S_n` is one chain, the chain monomial `x^{d}` with `d = _tp_word_degree(N)` **in CS degree `N = n − 1`** (the bank formula is indexed by CS degree, not quiverlab degree): `_tp_word_degree(N)`: `N<0 → 0` (vertex), even `N=2m → m·a + 1`, odd `N=2m+1 → (m+1)·a`. Thus quiverlab `S_n` carries `x^{_tp_word_degree(n-1)}`; for `a = 2`: `n=1→x, n=2→x², n=3→x³, … , x^n`. For the commutative square (Fixture B) the single tip `cd` has suffix `d ≠` prefix `c`, no self-overlap ⟹ `S_n = ∅ (n ≥ 3)`, resolution length 2 (CS Example: `𝒜^2` empty ⟹ finite resolution).

### P2. The resolution terms and the two collapse maps (unchanged; reviewer-verified correct)

`A ⊗_{A^e}(A e_o ⊗ e_t A) ≅ e_t A e_o` and `Hom_{A^e}(A e_o ⊗ e_t A, A) ≅ e_o A e_t`. If `d_n(1⊗σ⊗1) = Σ_i c_i·(a_i ⊗ τ_i ⊗ b_i)` (a `P_{n-1}` element; `a_i,b_i ∈ A`, `τ_i ∈ S_{n-1}`), then:
- **Homology** boundary on `C_n = ⨁_{σ} e_{t(σ)}Ae_{o(σ)}`: `(σ, w) ↦ Σ_i c_i·(τ_i, b_i·w·a_i)` — **right factor `b` on the left, left factor `a` on the right** (bank `resolutions_cs.py`: `val = b·(e_j·a)`). The `op`-twist of `A^e = A ⊗ A^{op}` under left-to-right composition is exactly why homology puts `b` left.
- **Cohomology** coboundary on `C^n = ⨁_{σ} e_{o(σ)}Ae_{t(σ)}`: cochain `f`, `(δ^n f)(τ) = Σ_i c_i·a_i·f(τ_i)·b_i` for `τ ∈ S_{n+1}` (left `a`, right `b`). Same term list, opposite pairing. The bank supplies only homology; quiverlab adds cohomology via this second collapse.

### P3. The differentials (CS §4 verbatim + §6 explicit; the swap fixed)

**The Bardzell leading maps `f_n` (CS §4, transcribed verbatim).** On `kQ ⊗_E k𝒜_n ⊗_E kQ`:
```
f_{-1}(a ⊗ b) = ab.
n even, q ∈ 𝒜_n, q = u_0⋯u_n = v_n⋯v_0 (left & right factorizations):
    f_n(1 ⊗ q ⊗ 1) = v_n ⊗ (v_{n-1}⋯v_0) ⊗ 1  −  1 ⊗ (u_0⋯u_{n-1}) ⊗ u_n.
n odd, q ∈ 𝒜_n:
    f_n(1 ⊗ q ⊗ 1) = Σ_{a p c = q, p ∈ 𝒜_{n-1}}  a ⊗ p ⊗ c.
```
These induce `δ_n = π_{n-1} ∘ f_n ∘ i_n` (`i_n` embeds `1⊗q⊗1` into `kQ`; `π_{n-1}` applies the normal form `β` to the outer `a,c` factors). `δ_{-1} = μ`, `δ_0` as in §6 below.

**Quiverlab index (binding).** quiverlab `d_n = δ_{n-1}^{CS}` (bank), so CS `𝒜_n` ↔ quiverlab `S_{n+1}` and the parities flip: **quiverlab `d_n` uses CS `f_{n-1}`**.
- **`d_n`, `n` ODD ⟺ CS `f_{n-1}` even ⟹ the 2-term map.** For `σ ∈ S_n` with blocks `u_0…u_{n-1}` (via `left_decomposition(σ, n)`), the term list (format `(coeff, a_word, target_chain_word, c_word)` meaning `coeff·(a ⊗ target ⊗ c)`, where `v_n = u_0`, `v_0 = u_{n-1}`):
  ```
  (+1,  u_0,   [u_1⋯u_{n-1}],  ()),          # v_n⊗(rest)⊗1 = first block u_0 on the LEFT
  (−1,  (),    [u_0⋯u_{n-2}],  u_{n-1}).      # 1⊗(rest)⊗u_n = last block u_{n-1} on the RIGHT
  ```
- **`d_n`, `n` EVEN ⟺ CS `f_{n-1}` odd ⟹ the big sum.** For each factorization `σ = a·p·c` with `p ∈ S_{n-1}` (`a,c` paths): `(+1, a, [p], c)`.

**THE SWAP BUG (fixed).** My earlier draft returned the odd/2-term case as `[(+1, (), P, u0), (−1, ulast, Q, ())]` — `u_0` on the RIGHT and `u_{n-1}` on the LEFT, the mirror image; and the even case with `a,c` swapped. On commutative fixtures this is invisible; on the non-commutative quantum CI it is wrong. The **correct** collapse of the fixed odd map: term 1 → `()·w·u_0 = w·u_0`, term 2 → `u_{n-1}·w` — matching the ported Bardzell "small" formula `(P; w·u_0) − (Q; u_{n-1}·w)` and CS §6 (below). The **swap-catching test** (Task 5) pins the quantum-CI term lists against CS §6 verbatim, where `a` vs `c` carry the `ξ`-asymmetry.

**Low degrees, explicit (CS §6, verbatim).** `φ_0` is the Fox derivative (read **left-to-right**: for a path `c = c_1 c_2 ⋯ c_ℓ` with `c_k ∈ Q_1`, split at each arrow):
```
φ_0(c) = Σ_{k=1}^{ℓ}  (c_1⋯c_{k-1}) ⊗ c_k ⊗ (c_{k+1}⋯c_ℓ)    ∈ A ⊗_E kQ_1 ⊗_E A.
```
(CS write `c = c_n⋯c_1` right-to-left; reading the *path* left-to-right and splitting at each arrow gives this — verified below against `d_1(yx)`.) Then (CS §6):
```
δ_0(π(a) ⊗ α ⊗ π(c)) = π(aα) ⊗ π(c) − π(a) ⊗ π(αc)        [α ∈ Q_1 = S_1]
d_1^{CS}(1 ⊗ s ⊗ 1) = φ_0(s) − φ_0(β(s))     for s ∈ 𝒜_1 = S = quiverlab S_2.
```
**quiverlab `d_2`** (= CS `d_1`) on a tip `s` with rule `s → f_s`: `d_2(1⊗s⊗1) = φ_0(s) − φ_0(β(s))`, where **`β(s) = rs.normal_form(s)`** — the *fully reduced* normal form, **not** the raw stored tail (Plan-03 tails are minimal, not inter-reduced; this normal-form call is the Pillar-4 fix for the "`d_2` on unreduced tails" bug). `φ_0` extends `k`-linearly to `β(s) = Σ λ_i b_i`.

**The correction, pinned by `d∘d = 0` (CS §4 Theorems 4.1/4.2, made concrete via §6).** The tails make the naive leading map fail `d²=0` on non-monomial algebras; CS repair it with a correction in the `≺`-lower span. CS **Theorem 4.1**: if `d_i` satisfy **(1)** `d_{i-1}d_i = 0` and **(2)** `(d_i − δ_i)(1⊗q⊗1) ∈ ⟨\overline{𝓛}_{i-1}^{≺}(q)⟩_k` for all `q ∈ 𝒜_i`, then the complex is exact. CS **Theorem 4.2**: such `d_i` exist. Here `\overline{𝓛}_{i-1}^{≺}(q) = {λ·π(b)⊗p⊗π(b') : b,b'∈B, p∈𝒜_{i-1}, λ·b p b' ≺ q}` — the parallel `P_{i-1}` basis terms **strictly `≺`-below `q`** (`≺` = the reduction order, CS §2, computed from the Plan-03 `order`). CS's own §6 computation for the quantum CI **fixes the correction coefficients by imposing `d²=0`** (their words: "The equality `d_{n-1}∘d_n=0` shows that making the choice `ε=(−1)^s` does the job"). quiverlab does exactly this, generally and deterministically:
```
d_n(1⊗σ⊗1) = δ_n(1⊗σ⊗1) + Σ_i γ_i · gen_i ,
   gen_i ranging over the finite basis of ⟨\overline{𝓛}_{n-1}^{≺}(σ)⟩ (all P_{n-1} terms (b, p, b') with
   p∈S_{n-1}, b,b'∈B, and the path b·p·b' strictly ≺ σ), and γ_i ∈ domain solving the LINEAR system
        d_{n-1}( δ_n(1⊗σ⊗1) + Σ_i γ_i gen_i ) = 0.
```
Existence is Theorem 4.2; any solution gives an exact complex (Theorem 4.1). Solve over the `Domain` with an exact solver; if the system is **inconsistent**, the algebra needs the higher CS homotopy correction outside v1's construction — raise `NotImplementedError` at that exact `(n, σ)`. The two executable gates certify each instance: `assert_dd_zero` (condition 1) and `assert_order_condition` (condition 2: every correction term is `≺`-below `σ`).

**Convention translation (binding, per formula).** (1) `φ_0` and block products are read **left-to-right** as written; (2) CS's *right* `n`-ambiguity `v_n⋯v_0` equals the left one read left-to-right (`v_{n-i} = u_i` in the quadratic case, CS Prop), so `v_n = u_0` (first block), `v_0 = u_n` (last block) — used in the `n`-odd map; (3) homology collapse puts `b` (right factor) on the left; (4) `≺` is the Plan-03 admissible order restricted to parallel paths. Each is stated where used.

### P4. Comparison morphisms and operation transport (CS↔bar; spec §5 comp. 6) — unchanged posture

Operations (cup `⌣`, cap `∩`, bracket `[,]`) are **computed on bar cochains** by the ported `engine/tt_calculus.py` (Plan 02); Plan 04 delivers only the **comparison maps** connecting CS classes to bar representatives (Negron–Witherspoon arXiv:1406.0036 / Volkov arXiv:1610.05741). `Φ: P^{CS} → Bar` (`= i_• ∘ ι` on blocks) and `Ψ: Bar → P^{CS}` (the CS projection), chain maps, mutually quasi-inverse, inducing `HH^•(CS) ≅ HH^•(bar)`. Deliverables: `Φ`, `Ψ`, `transport_cocycle_cs_to_bar`, `transport_class_bar_to_cs`, and `cup_of_cs_classes`/`bracket_of_cs_classes` (transport → `tt_calculus` → optional transport back). **Boundary:** valid only inside the bar-buildable degree window; past it, `NotImplementedError` ("native CS cup is a later phase"). Dimensions go deep; operations are window-bounded.

---

## Hand-verified fixtures (full matrices)

**Fixture A — `A = k[x]/(x²)`** (bank-confirmed: `HH_• = [2,1,1,1,1,1,1]` char 0, `[2,2,…]` char 2). One vertex `v`, loop `x`; `B = {1,x}`, `dim A = 2`; rule `x² → 0` (tip `xx`, `f_s = 0`).
`S_n = {x^n}` (Pillar-1 re-base: `pres_S = MonomialPresentation.truncated_polynomial(2)`, `associated_paths(n)` = the single chain `x^{_tp_word_degree(n-1)} = x^n`). `P_n ≅ A^e`; collapsed `C_n = A`, `dim = 2`.
Differentials (CS specializes to Bardzell; monomial ⟹ no correction): with `u = x⊗1 − 1⊗x`, `v = 1⊗x + x⊗1`, `d_n = u·` (`n` odd), `d_n = v·` (`n` even). Homology collapse on `C_n = ⟨1,x⟩` (char 0): `n` odd → `x·(−)−(−)·x = 0`; `n` even → `2x·(−)`, matrix `[[0,0],[2,0]]`, rank 1. `HH_• = [2,1,1,…]`; char 2 (`2 = 0`): all `d = 0`, `HH_• = [2,2,…]`. Cohomology identical (`A` commutative): `HH^0 = 2`, `HH^n = 1 (n≥1)` char 0. These are the Task-5/6/7 assertions.

**Fixture B — the commutative square `A = kQ/(ab − cd)`** (non-monomial; `= kA_2 ⊗ kA_2`, the diamond-poset incidence algebra). `Q`: `1,2,3,4`; arrows `a:1→2, b:2→4, c:1→3, d:3→4`; arrow order `a<b<c<d` ⟹ tip `cd`, rule `cd → ab`. `B = {e_1..e_4, a,b,c,d, ab}`, `dim A = 9`. `S_2 = {cd}`, `S_n = ∅ (n≥3)`, length 2.
`d_2(1⊗cd⊗1) = φ_0(cd) − φ_0(β(cd)) = φ_0(cd) − φ_0(ab) = (1⊗c⊗d + c⊗d⊗1) − (1⊗a⊗b + a⊗b⊗1)`. `d_1(1⊗α⊗1) = α⊗1 − 1⊗α`.
*Cohomology complex* `C^n = ⨁_{σ∈S_n} e_{o}Ae_{t}`: `C^0 = ⟨ê_1..ê_4⟩` (4), `C^1 = ⟨â,b̂,ĉ,d̂⟩` (4), `C^2 = e_1Ae_4 = ⟨âb⟩` (1).
`δ^0`: `f=(λ_1..λ_4) ↦ (δ^0f)(α:i→j) = (λ_j − λ_i)α`; matrix (rows `a,b,c,d`; cols `e_1..e_4`) `[[-1,1,0,0],[0,-1,0,1],[-1,0,1,0],[0,0,-1,1]]`, rank 3, `ker = ⟨(1,1,1,1)⟩ =` center. `δ^1`: `g(a)=αa,…,g(d)=δd ↦ (δ^1g)(cd) = c·g(d)+g(c)·d−a·g(b)−g(a)·b = (δ+γ−β−α)ab`; matrix (row `cd`; cols `a,b,c,d`) `[−1,−1,1,1]`, rank 1.
`HH^0 = 4−3 = 1`, `HH^2 = 1−1 = 0`, `HH^1 = (4−1)−3 = 0`. **`HH^• = (1,0,0)`.** *Homology complex* `C_n = ⨁ e_tAe_o`: `C_0 = k^4`, `C_1 = 0`, `C_2 = 0` (acyclic). **`HH_• = (4,0,0)`.**
*Three cross-checks:* Euler `4−4+1 = 1 = HH^0`; **Gerstenhaber–Schack** (`HH^•(kP) ≅ H^•(|order complex|)`, diamond order complex = two triangles glued on `{1,4}` = contractible) ⟹ `HH^0=1, HH^{>0}=0`; **Künneth** (Plan-03 Fixture-2 note: `A = kA_2 ⊗ kA_2`, `kA_2` hereditary tree ⟹ `HH^0=k, HH^{≥1}=0`) ⟹ `HH^•(A) = (1,0,0)` over any field. `A` not symmetric ⟹ `HH_0 = 4 ≠ 1 = HH^0` (correct).

**Fixture C — the quantum CI `A = k⟨x,y⟩/(x²,y²,yx−ξxy)`** = CS §6 "the algebra counterexample to Happel's question" (the **non-commutative swap-catcher** and the byte-oracle family). `B = {1,x,y,xy}`, `dim A = 4`; tips `S = {x², y², yx}`, rule `yx → ξxy`. `𝒜_n = {y^s x^t : s+t = n+1}`, quiverlab `S_n = {y^s x^t : s+t = n}`, `|S_n| = n+1`. **CS §6 gives the resolution explicitly** — Task-5/6/11 assertions, verbatim:
```
d_1^{CS}(1⊗x²⊗1) = x⊗x⊗1 + 1⊗x⊗x
d_1^{CS}(1⊗y²⊗1) = y⊗y⊗1 + 1⊗y⊗y
d_1^{CS}(1⊗yx⊗1) = y⊗x⊗1 + 1⊗y⊗x − ξ·x⊗y⊗1 − ξ·1⊗x⊗y      [swap-catcher: a=y LEFT in y⊗x⊗1; c=x RIGHT in 1⊗y⊗x]
d_2^{CS}(1⊗x³⊗1)  = x⊗x²⊗1 − 1⊗x²⊗x
d_2^{CS}(1⊗y²x⊗1) = y⊗yx⊗1 + ξ·1⊗yx⊗y + ξ²·x⊗y²⊗1 − 1⊗y²⊗x
d_2^{CS}(1⊗yx²⊗1) = y⊗x²⊗1 − 1⊗yx⊗x − ξ·x⊗yx⊗1 − ξ²·1⊗x²⊗y
d_2^{CS}(1⊗y³⊗1)  = y⊗y²⊗1 − 1⊗y²⊗y
```
General closed form (CS §6, `s>0, t>0`, CS index `m`, `s+t = m+1`):
```
d_m^{CS}(1⊗y^s x^t⊗1) = y⊗y^{s-1}x^t⊗1 + (−1)^{m+1}·1⊗y^s x^{t-1}⊗x
                        + (−1)^s ξ^s·x⊗y^s x^{t-1}⊗1 + (−1)^s ξ^t·1⊗y^{s-1}x^t⊗y
```
(pure powers `y^{m+1}`, `x^{m+1}` drop the last two terms). The first two terms are `δ_m` (leading); the last two are the `≺`-lower correction with `ε = (−1)^s` **pinned by `d²=0`** — exactly what Task 6's linear solve reproduces. **Bank-confirmed `HH_•` (ξ=2, degree 0..12), from running the read-only oracle:** char 0 `[3,2,2,2,2,2,2,2,2,2,2,2,2]`; p=2 `[3,4,4,…]`; p=3 `[3,4,6,8,10,12,14,16,18,20,22,24,26]` (linear); p=5 `[3,2,3,4,3,2,4,6,4,2,5,8,5]`. These are the Task-10/11/14 assertions.

---

### Task 1: Interface freshness gate against the REAL committed Plan-03 reduction system

> Plan 04 consumes Plan 03's `groebner` package, **committed on `origin/main`** (`docs/plans/2026-07-18-plan-03-groebner.md`). This gate encodes the **real** frozen contract: words are arrow-**name** tuples; `ReductionRule.tail` is a tuple of `(coeff, word)` pairs (not a dict); `Ambiguity(kind, word, left, right, a, c)`; there is **no** `A.reduction_system()`. If `groebner` is missing or has drifted, the gate FAILS — STOP and reconcile before any further Plan-04 work.

**Files:** Create branch `plan-04-chouhy-solotar`; `tests/resolutions_cs/__init__.py` (empty); `tests/resolutions_cs/conftest.py`; `tests/resolutions_cs/test_interface_gate.py`.

**Interfaces:** Consumes `quiverlab.groebner.{ReductionSystem, ReductionRule, Ambiguity, PathOrder, build_reduction_system}`. Produces the freshness gate + shared fixtures.

- [ ] **Step 1: Branch.** Run `cd /Users/marco/Desktop/HomologicalNetworks/quiverlab && git checkout -b plan-04-chouhy-solotar`. Expected `Switched to a new branch 'plan-04-chouhy-solotar'`.

- [ ] **Step 2: Shared fixtures** `tests/resolutions_cs/conftest.py`:
```python
import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.groebner import build_reduction_system

# Reduction systems are the CS entry currency (there is NO A.reduction_system()).
@pytest.fixture
def kx2_rs():
    Q = Quiver([1], {"x": (1, 1)})
    return build_reduction_system(Q, ["x*x"], CC)

@pytest.fixture
def square_rs():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return build_reduction_system(Q, ["a*b - c*d"], CC)

@pytest.fixture
def qci_rs():
    def build(xi="2", field=CC):
        Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
        return build_reduction_system(Q, ["x*x", "y*y", f"y*x - ({xi})*x*y"], field)
    return build
```

- [ ] **Step 3: Freshness gate** `tests/resolutions_cs/test_interface_gate.py`:
```python
"""Plan-04 interface freshness gate: the REAL frozen Plan-03 reduction-system API
(docs/plans/2026-07-18-plan-03-groebner.md). STOP on any drift."""
from dataclasses import fields as dc_fields
import pytest

groebner = pytest.importorskip("quiverlab.groebner",
                               reason="Plan 03 (groebner) must be committed before Plan 04.")


def test_frozen_symbols():
    for n in ("ReductionSystem", "ReductionRule", "Ambiguity", "PathOrder",
              "build_reduction_system"):
        assert hasattr(groebner, n), f"groebner.{n} missing (Plan-03 drift)"


def test_reduction_system_shape():
    names = {f.name for f in dc_fields(groebner.ReductionSystem)}
    assert {"quiver", "domain", "order", "rules", "irreducibles",
            "degree_bound", "is_confluent"} <= names
    for m in ("leading_words", "reduce", "normal_form", "ambiguities"):
        assert callable(getattr(groebner.ReductionSystem, m))


def test_rule_and_ambiguity_shape(square_rs):
    rule = square_rs.rules[0]
    assert rule.lead == ("c", "d")                                   # words are ARROW-NAME tuples
    assert rule.tail == ((square_rs.domain.coerce(1), ("a", "b")),)  # tuple of (coeff, word), NOT a dict
    assert (rule.source, rule.target) == (1, 4)
    for amb in square_rs.ambiguities():
        assert amb.kind in ("overlap", "inclusion")
        assert isinstance(amb.word, tuple) and amb.left is not None and amb.right is not None


def test_normal_form_reduces_tip(square_rs):
    one = square_rs.domain.coerce(1)
    assert square_rs.normal_form(("c", "d")) == {("a", "b"): one}    # cd -> ab (a full dict)


def test_no_algebra_reduction_system_accessor():
    from quiverlab import Quiver, CC
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=CC)
    assert not hasattr(A, "reduction_system")                        # CS gets the RS via build_reduction_system
```

- [ ] **Step 4: Run.** `cd /Users/marco/Desktop/HomologicalNetworks/quiverlab && NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/resolutions_cs/test_interface_gate.py -q`
Expected `5 passed`. If any assertion FAILS: STOP; reconcile the rename with the Plan-03 author and update *every* Task 2–14 call site + this gate in one commit, recording the drift. (`importorskip` fires only if `groebner` is absent — itself the STOP signal.)

- [ ] **Step 5: Field-object contract probe (STOP-and-reconcile, for the Task-7 `_DomainField` shim).** Task 7's `reduction_system_of(A)` calls `build_reduction_system(A.quiver, A.relations, field_for_domain(A.domain))` with a shim standing in for the `field` object; the shim must implement **exactly** the methods `build_reduction_system` calls on `field`. Determine them from the committed source:
```bash
cd /Users/marco/Desktop/HomologicalNetworks/quiverlab && \
  grep -nE 'field\.[A-Za-z_]+' src/quiverlab/groebner/system.py
```
Record every `field.<method>` hit. The committed Plan-03 `build_reduction_system` uses `field.parse_entry(x)` and `field.make_domain(entries)`; **crucially, relation coefficients arrive as `fractions.Fraction`**, and Plan-03 later does `dom.coerce(field.parse_entry(c))`, so the shim's `parse_entry` must return something `A.domain.coerce` accepts (pass the `Fraction` straight through; `make_domain` returns the existing `A.domain`). If the grep shows **any** additional method (e.g. `field.name`, `field.characteristic`), STOP and add exactly those to `_DomainField` in Task 7 (do not guess the surface). Write the recorded contract as a comment at the top of `tests/resolutions_cs/test_interface_gate.py` so Task 7 consumes the verified list, not an assumption.

- [ ] **Step 6: Commit**
```bash
git add tests/resolutions_cs/__init__.py tests/resolutions_cs/conftest.py tests/resolutions_cs/test_interface_gate.py
git commit -m "test(cs): freshness gate for the committed Plan-03 reduction-system interface

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Package skeleton, ambiguity-chain type, trace dataclasses

**Files:** Create `src/quiverlab/resolutions_cs/__init__.py`, `terms.py`, `trace.py`; `tests/resolutions_cs/test_types.py`.

**Interfaces:** Produces `Chain(word, blocks, o, t, degree)` and the four trace dataclasses. No Plan-03 dependency yet.

- [ ] **Step 1: Failing test** `tests/resolutions_cs/test_types.py`:
```python
from quiverlab.resolutions_cs.terms import Chain
from quiverlab.resolutions_cs import trace


def test_chain_records_blocks_and_endpoints():
    ch = Chain(word=("x", "x", "x"), blocks=(("x",), ("x",), ("x",)), o=1, t=1, degree=3)
    assert ch.n_blocks == 3 and ch.degree == 3 and ch.word == ("x", "x", "x")


def test_chain_hashable_and_equal():
    a = Chain(("x", "x"), (("x",), ("x",)), 1, 1, 2)
    b = Chain(("x", "x"), (("x",), ("x",)), 1, 1, 2)
    assert a == b and hash(a) == hash(b) and {a, b} == {a}


def test_trace_dataclasses_are_inert():
    ev = trace.DifferentialEvent(degree=2, chain=("c", "d"), terms=[(1, (), ("c",), ("d",))])
    assert ev.degree == 2 and not hasattr(ev, "render")
```

- [ ] **Step 2: Run — expect FAIL** `ModuleNotFoundError: quiverlab.resolutions_cs`.

- [ ] **Step 3: Implement.** `src/quiverlab/resolutions_cs/__init__.py`:
```python
"""quiverlab.resolutions_cs: the Chouhy-Solotar general bimodule resolution
(arXiv:1406.2300 = J. Algebra 432 (2015)). Domain-generic; consumes the Plan-03
groebner.ReductionSystem and the Plan-01 core.Algebra. See the plan and
docs/internals/08-chouhy-solotar.md for the mathematics (S-sequence = Bardzell
associated paths of the tip monomial algebra; differentials = CS f_n leading map +
order-condition-pinned correction; two collapse maps for HH_• and HH^•)."""
```
`src/quiverlab/resolutions_cs/terms.py`:
```python
"""Ambiguity-chain type for the CS S-sequence. Pure combinatorics; words are tuples of
arrow NAMES read left-to-right (matching the Plan-03 reduction system)."""
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class Chain:
    """sigma in S_n = 𝒜_{n-1}. `word`: the underlying path (arrow names). `blocks`: the
    unique CS left decomposition u_0|...|u_{n-1} (n blocks). `o`,`t`: source/target. `degree`: n."""
    word: tuple
    blocks: tuple
    o: object
    t: object
    degree: int

    @property
    def n_blocks(self) -> int:
        return len(self.blocks)
```
`src/quiverlab/resolutions_cs/trace.py`:
```python
"""Plain trace dataclasses. Plan 07 renders these; Plan 04 only populates them and
asserts their claims equal computed values."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any


@dataclass
class AmbiguityEvent:
    degree: int
    chain_words: list


@dataclass
class ResolutionTerm:
    degree: int
    n_generators: int          # |S_n|
    collapsed_dim: int         # dim C_n / dim C^n


@dataclass
class DifferentialEvent:
    degree: int
    chain: Any                 # source chain word
    terms: list                # [(coeff, a_word, target_word, c_word), ...]


@dataclass
class LiftStep:
    degree: int
    kind: str                  # "delta" | "correction-solve" | "dd-check" | "order-check"
    detail: Any = None
```

- [ ] **Step 4: Run — expect PASS** `3 passed`.

- [ ] **Step 5: Commit**
```bash
git add src/quiverlab/resolutions_cs/__init__.py src/quiverlab/resolutions_cs/terms.py src/quiverlab/resolutions_cs/trace.py tests/resolutions_cs/test_types.py
git commit -m "feat(cs): package skeleton, ambiguity-chain type, trace dataclasses

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: The S-sequence via the ported Bardzell associated-paths (Pillar-1 re-base)

**Files:** Create `src/quiverlab/resolutions_cs/ambiguities.py`; `tests/resolutions_cs/test_ambiguities.py`.

**Interfaces:** Consumes `groebner.ReductionSystem` (`.quiver`, `.leading_words()`, `.order`), `engine.resolutions_bardzell.MonomialPresentation`, `terms.Chain`, `DepthLimitError`. Produces `SSequence(rs, max_degree, max_cells)` with `.S(n) -> list[Chain]` (memoized, guarded), `.tip_presentation()`.

- [ ] **Step 1: Failing test** `tests/resolutions_cs/test_ambiguities.py`:
```python
import pytest
from quiverlab.errors import DepthLimitError
from quiverlab.resolutions_cs.ambiguities import SSequence
pytest.importorskip("quiverlab.groebner")


def test_kx2_single_chain_per_degree(kx2_rs):
    ss = SSequence(kx2_rs, max_degree=6)
    assert [len(ss.S(n)) for n in range(7)] == [1, 1, 1, 1, 1, 1, 1]
    assert ss.S(2)[0].word == ("x", "x")
    assert ss.S(3)[0].word == ("x", "x", "x")
    assert ss.S(3)[0].blocks == (("x",), ("x",), ("x",))          # x|x|x, three blocks


def test_commutative_square_length_two(square_rs):
    ss = SSequence(square_rs, max_degree=5)
    assert len(ss.S(0)) == 4 and len(ss.S(1)) == 4 and len(ss.S(2)) == 1
    assert ss.S(2)[0].word == ("c", "d") and ss.S(2)[0].o == 1 and ss.S(2)[0].t == 4
    assert ss.S(3) == [] and ss.S(4) == []


def test_quantum_ci_chain_counts(qci_rs):
    ss = SSequence(qci_rs(xi="2"), max_degree=6)
    assert [len(ss.S(n)) for n in range(7)] == [1, 2, 3, 4, 5, 6, 7]  # |S_n| = n+1
    assert {c.word for c in ss.S(2)} == {("x", "x"), ("y", "y"), ("y", "x")}
    assert {c.word for c in ss.S(3)} == {("x",)*3, ("y", "x", "x"), ("y", "y", "x"), ("y",)*3}


def test_depth_guard_loud(qci_rs):
    ss = SSequence(qci_rs(xi="2"), max_degree=60, max_cells=5)
    with pytest.raises(DepthLimitError) as e:
        ss.S(40)
    assert "certified" in str(e.value).lower()
```

- [ ] **Step 2: Run — expect FAIL** `ModuleNotFoundError: …ambiguities`.

- [ ] **Step 3: Implement** `src/quiverlab/resolutions_cs/ambiguities.py`:
```python
"""The CS ambiguity S-sequence (arXiv:1406.2300 §3). S_n = 𝒜_{n-1}: S_0=vertices,
S_1=arrows, S_2=tips, S_n(n≥3)=(n-1)-ambiguities. Since ambiguities depend ONLY on the
tips, S_n for A equals Bardzell's AP^n for the tip monomial algebra A_S = kQ/⟨S⟩. We
reuse the ported, bar-oracle-validated
engine.resolutions_bardzell.MonomialPresentation.{associated_paths, left_decomposition}
(Bardzell AP^n = CS 𝒜_{n-1} = quiverlab S_n)."""
from quiverlab.errors import DepthLimitError
from quiverlab.engine.resolutions_bardzell import MonomialPresentation
from quiverlab.resolutions_cs.terms import Chain

_HINT = ("CS ambiguity growth guard; raise max_cells only if the certified range is "
         "genuinely insufficient (SSequence(..., max_degree=D) reuses cached lower degrees)")


class SSequence:
    def __init__(self, rs, max_degree, max_cells=4_000_000):
        self.rs = rs
        self.max_degree = max_degree
        self.max_cells = max_cells
        self._cache = {}
        Q = rs.quiver
        self._vertices = list(Q.vertices)
        arrows = [(name, Q.source(name), Q.target(name)) for name in Q.arrows]  # arrow NAME as id
        self._pres = MonomialPresentation(self._vertices, arrows, list(rs.leading_words()))
        self._maxlen = max(2, max_degree) * max(self._pres.maxrel, 2) + 2       # n-ambiguity length cap

    def tip_presentation(self):
        return self._pres

    def _guard(self, n, count):
        if count > self.max_cells:
            last = max((k for k in self._cache if self._cache.get(k)), default=n - 1)
            raise DepthLimitError(
                f"CS S-sequence at degree {n}: {count} chains exceed max_cells="
                f"{self.max_cells}; certified through degree {last}", hint=_HINT)

    def S(self, n):
        if n in self._cache:
            return self._cache[n]
        if n < 0 or n > self.max_degree:
            return []
        if n == 0:
            out = [Chain((), (), v, v, 0) for v in self._vertices]
        else:
            words = self._pres.associated_paths(n, self._maxlen)   # Bardzell AP^n = S_n
            out = []
            for w in words:
                blocks = tuple(tuple(b) for b in self._pres.left_decomposition(w, n))
                out.append(Chain(tuple(w), blocks,
                                 self._pres.path_src(w), self._pres.path_tgt(w), n))
        self._guard(n, len(out))
        out.sort(key=lambda ch: (self.rs.order.key(ch.word), ch.word))
        self._cache[n] = out
        return out
```
> `MonomialPresentation(vertices, arrows, relations)`, `.associated_paths`, `.left_decomposition`, `.path_src`, `.path_tgt`, `.maxrel` are the ported bank API (verified `src/quiverlab/engine/resolutions_bardzell.py` lines 111–241). `Q.vertices`, `Q.arrows`, `Q.source(name)`, `Q.target(name)`, `rs.order.key` are frozen Plan-01/03.

- [ ] **Step 4: Run — expect PASS** `4 passed`. (Wrong quantum-CI counts ⟹ the tips were not passed as `relations` verbatim.)

- [ ] **Step 5: Commit**
```bash
git add src/quiverlab/resolutions_cs/ambiguities.py tests/resolutions_cs/test_ambiguities.py
git commit -m "feat(cs): S-sequence via ported Bardzell associated-paths on the tip algebra

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: The `A`-arithmetic layer and the two collapse corners

**Files:** Create `src/quiverlab/resolutions_cs/aarith.py`; `tests/resolutions_cs/test_aarith.py`.

**Interfaces:** Consumes `core.Algebra` (`.domain, .dim, .multiply, ._basis_vec`), `groebner.ReductionSystem`. Produces `AArith(A, rs)` with `.path_vec(word)`, `.mul(u,v)`, `.vertex_vec(v)`, `.corner(o,t,side) -> list[int]`, `.basis_word(idx)`.

- [ ] **Step 1: Failing test** `tests/resolutions_cs/test_aarith.py`:
```python
import pytest
from quiverlab import Quiver, CC
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.aarith import AArith
from quiverlab.resolutions_cs.ambiguities import SSequence
pytest.importorskip("quiverlab.groebner")


def _square():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    return A, build_reduction_system(Q, ["a*b - c*d"], CC)


def test_path_vec_reduces_cd_to_ab():
    A, rs = _square()
    ar = AArith(A, rs)
    assert ar.path_vec(("c", "d")) == ar.path_vec(("a", "b"))       # cd normal-forms to ab
    assert ar.path_vec(("a", "b")) != [ar.dom.zero()] * A.dim


def test_square_cohomology_corner_dims():
    A, rs = _square()
    ar, ss = AArith(A, rs), SSequence(rs, 3)
    assert [sum(len(ar.corner(c.o, c.t, "coh")) for c in ss.S(n)) for n in range(3)] == [4, 4, 1]


def test_square_homology_corner_dims():
    A, rs = _square()
    ar, ss = AArith(A, rs), SSequence(rs, 3)
    assert [sum(len(ar.corner(c.o, c.t, "hom")) for c in ss.S(n)) for n in range(3)] == [4, 0, 0]


def test_kx2_corner_is_full_algebra():
    Q = Quiver([1], {"x": (1, 1)})
    A = Q.algebra(relations=["x*x"], field=CC)
    ar = AArith(A, build_reduction_system(Q, ["x*x"], CC))
    assert len(ar.corner(1, 1, "hom")) == A.dim                      # e_v A e_v = A, dim 2
```

- [ ] **Step 2: Run — expect FAIL** `ModuleNotFoundError: …aarith`.

- [ ] **Step 3: Implement** `src/quiverlab/resolutions_cs/aarith.py`:
```python
"""A-arithmetic for the CS resolution over any Domain. A path word maps to its normal-form
A-vector by FOLDING core.Algebra.multiply over the arrow generators (A's structure constants
already encode reduction mod I, so this IS the map β/π). The vertex/irreducible basis order
matches groebner_algebra: vertices first, then rs.irreducibles."""


class AArith:
    def __init__(self, A, rs):
        self.A = A
        self.rs = rs
        self.dom = A.domain
        Q = rs.quiver
        self._vertices = list(Q.vertices)
        self._basis_words = [("v", v) for v in self._vertices] + [("p", w) for w in rs.irreducibles]
        self._word_index = {w: i for i, w in enumerate(self._basis_words)}
        self._arrow_vec = {name: self._lookup_vec(("p", (name,))) for name in Q.arrows}
        self._vertex_idem = {v: self._lookup_vec(("v", v)) for v in self._vertices}

    def _lookup_vec(self, tagged_word):
        v = [self.dom.zero()] * self.A.dim
        v[self._word_index[tagged_word]] = self.dom.one()
        return v

    def basis_word(self, idx):
        return self._basis_words[idx]

    def mul(self, u, v):
        return self.A.multiply(u, v)

    def vertex_vec(self, vtx):
        return list(self._vertex_idem[vtx])

    def path_vec(self, word):
        """Normal-form A-vector of a NONEMPTY path (arrow-name tuple), left-to-right.
        Callers use vertex_vec for the empty path e_v."""
        if len(word) == 0:
            raise ValueError("path_vec needs a nonempty path; use vertex_vec for e_v")
        acc = list(self._arrow_vec[word[0]])
        for name in word[1:]:
            acc = self.A.multiply(acc, self._arrow_vec[name])
        return acc

    def corner(self, o, t, side):
        """Basis indices j with b_j in the corner: side="hom" -> e_t A e_o (paths t->o);
        side="coh" -> e_o A e_t (paths o->t)."""
        left, right = (t, o) if side == "hom" else (o, t)
        el, er = self._vertex_idem[left], self._vertex_idem[right]
        return [j for j in range(self.A.dim)
                if self.A.multiply(el, self.A.multiply(self.A._basis_vec(j), er)) == self.A._basis_vec(j)]
```
> The `("v", v)` / `("p", word)` tagging must match `groebner_algebra`'s basis (vertices first, then `rs.irreducibles` in sorted order). If Plan-03 labels differently, adapt `_basis_words` only; `test_path_vec_reduces_cd_to_ab` catches a mismatch. `A._basis_vec` is the Plan-01 helper.

- [ ] **Step 4: Run — expect PASS** `4 passed`.

- [ ] **Step 5: Commit**
```bash
git add src/quiverlab/resolutions_cs/aarith.py tests/resolutions_cs/test_aarith.py
git commit -m "feat(cs): A-arithmetic layer (path normal forms) and HH/HH collapse corners

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: The leading map `δ_n` (CS `f_n`, corrected) + low-degree `d_1, d_2` + the swap-catching test

**Files:** Create `src/quiverlab/resolutions_cs/resolution.py` (partial); `tests/resolutions_cs/test_leading.py`.

**Interfaces:** Consumes `SSequence`, `AArith`, `groebner.ReductionSystem` (`.normal_form`). Produces `ChouhySolotarResolution(A, rs, max_degree)` with `.delta_terms(n, chain)` (CS leading map) and `.d_terms(n, chain)` for `n ∈ {1,2}`.

- [ ] **Step 1: Failing test** — the reviewer-mandated **non-commutative swap-catcher** pins the quantum-CI term lists against CS §6 verbatim.
```python
import pytest
from quiverlab import Quiver, CC
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
pytest.importorskip("quiverlab.groebner")


def _res(field=CC, xi="2"):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", f"y*x - ({xi})*x*y"]
    A = Q.algebra(relations=rels, field=field)
    return ChouhySolotarResolution(A, build_reduction_system(Q, rels, field), max_degree=6)


def _norm(res, terms):
    return {(res.to_int(c), a, t, cc) for (c, a, t, cc) in terms}


def test_d2_quantum_ci_yx_matches_paper_no_swap():
    """CS §6: d_1^{CS}(1⊗yx⊗1) = y⊗x⊗1 + 1⊗y⊗x − ξ x⊗y⊗1 − ξ 1⊗x⊗y  (ξ=2). a=y is
    LEFT of the first term, c=x is RIGHT of the second: the swap-catcher."""
    res = _res()
    yx = next(c for c in res.ss.S(2) if c.word == ("y", "x"))
    assert _norm(res, res.d_terms(2, yx)) == {
        ( 1, ("y",), ("x",), ()),          #  y ⊗ x ⊗ 1
        ( 1, (),     ("y",), ("x",)),       #  1 ⊗ y ⊗ x
        (-2, ("x",), ("y",), ()),           # −ξ x ⊗ y ⊗ 1
        (-2, (),     ("x",), ("y",)),       # −ξ 1 ⊗ x ⊗ y
    }


def test_delta_leading_x3_not_mirror_reversed():
    """CS §6: δ for x³ = x⊗x²⊗1 − 1⊗x²⊗x. a=x LEFT of the +term, c=x RIGHT of the −term
    (the exact positions my earlier draft had reversed)."""
    res = _res()
    x3 = next(c for c in res.ss.S(3) if c.word == ("x", "x", "x"))
    d = _norm(res, res.delta_terms(3, x3))
    assert d == {(1, ("x",), ("x", "x"), ()), (-1, (), ("x", "x"), ("x",))}


def test_square_d2_fox_derivative():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    res = ChouhySolotarResolution(A, build_reduction_system(Q, ["a*b - c*d"], CC), max_degree=3)
    cd = res.ss.S(2)[0]
    assert _norm(res, res.d_terms(2, cd)) == {
        ( 1, (),     ("c",), ("d",)),       #  1 ⊗ c ⊗ d
        ( 1, ("c",), ("d",), ()),           #  c ⊗ d ⊗ 1
        (-1, (),     ("a",), ("b",)),       # −1 ⊗ a ⊗ b
        (-1, ("a",), ("b",), ()),           # −a ⊗ b ⊗ 1
    }
```

- [ ] **Step 2: Run — expect FAIL** `ModuleNotFoundError: …resolution`.

- [ ] **Step 3: Implement** the leading map + low degrees in `src/quiverlab/resolutions_cs/resolution.py`:
```python
"""The Chouhy-Solotar differential (arXiv:1406.2300 §4, §6), Domain-generic. A TERM is
(coeff, a_word, target_word, c_word) meaning coeff·(a ⊗ target ⊗ c) in P_{n-1}, where
a_word,c_word are PATHS (arrow-name tuples; () = the unit) reduced to normal form only at
collapse time, and target_word identifies a chain in S_{n-1}. Composition is LEFT TO RIGHT."""
from quiverlab.resolutions_cs.ambiguities import SSequence
from quiverlab.resolutions_cs.aarith import AArith


class ChouhySolotarResolution:
    def __init__(self, A, rs, max_degree, max_cells=4_000_000):
        self.A = A
        self.rs = rs
        self.dom = A.domain
        self.ss = SSequence(rs, max_degree, max_cells)
        self.ar = AArith(A, rs)
        self._chain_index = {}
        self._d_cache = {}

    def to_int(self, c):
        return self.dom.to_int(c) if hasattr(self.dom, "to_int") else int(c)

    def _one(self):
        return self.dom.one()

    def _neg(self, c):
        return self.dom.neg(c)

    def _chain(self, degree, word):
        idx = self._chain_index.get(degree)
        if idx is None:
            idx = {c.word: c for c in self.ss.S(degree)}
            self._chain_index[degree] = idx
        return idx.get(tuple(word))

    # -- Fox derivative φ_0 (CS §6), left-to-right --------------------------
    def _fox(self, word, coeff):
        return [(coeff, word[:k], (word[k],), word[k + 1:]) for k in range(len(word))]

    # -- leading map δ_n (CS f_{n-1}), quiverlab index ----------------------
    def delta_terms(self, n, chain):
        if n == 1:
            return self._d1_terms(chain)
        one, none = self._one(), self._neg(self._one())
        blocks = chain.blocks
        if n % 2 == 1:                                            # CS f_{n-1} EVEN: 2-term map
            u0, ulast = blocks[0], blocks[-1]
            P = tuple(x for blk in blocks[1:] for x in blk)      # u_1..u_{n-1} in S_{n-1}
            Q = tuple(x for blk in blocks[:-1] for x in blk)     # u_0..u_{n-2} in S_{n-1}
            return [(one, u0, P, ()), (none, (), Q, ulast)]
        prev = {c.word for c in self.ss.S(n - 1)}                # n even: CS f_{n-1} ODD, big sum
        w, out = chain.word, []
        for i in range(len(w)):
            for j in range(i + 1, len(w) + 1):
                if w[i:j] in prev:
                    out.append((one, w[:i], w[i:j], w[j:]))
        return out

    # -- d_1 (arrows → vertices), CS §6 -------------------------------------
    def _d1_terms(self, chain):
        one, none = self._one(), self._neg(self._one())
        name = chain.word[0]
        return [(one, (name,), ("__v__", chain.t), ()),         # α ⊗ e_t ⊗ 1
                (none, (), ("__v__", chain.o), (name,))]         # 1 ⊗ e_o ⊗ α

    # -- d_2 = φ_0(s) − φ_0(β(s)), CS §6 (β = FULL normal form, Pillar-4 fix) --
    def _d2_terms(self, chain):
        s = chain.word
        terms = list(self._fox(s, self._one()))                  # φ_0(s)
        for word, coeff in self.rs.normal_form(s).items():       # β(s) = Σ λ_i b_i (fully reduced dict)
            for (c, a, tw, cc) in self._fox(word, coeff):
                terms.append((self._neg(c), a, tw, cc))          # − φ_0(β(s))
        return terms

    def d_terms(self, n, chain):
        if n == 1:
            return self._d1_terms(chain)
        if n == 2:
            return self._d2_terms(chain)
        return self._d_general(n, chain)                         # Task 6

    def _d_general(self, n, chain):
        raise NotImplementedError("filled in Task 6 (order-condition-pinned correction)")
```
> Notes: (1) `d_2` targets `S_1` arrows (`_fox` returns `(arrow,)` chain words); `δ_n (n≥3)` targets `S_{n-1}` words. (2) `d_1` targets the two vertex chains, encoded `("__v__", vertex)`; the collapse (Task 6/7) maps these to the `S_0` chains. (3) `to_int` is test-only; if `Domain` lacks it, compare via `dom.eq(c, dom.coerce(k))` in the test helper.

- [ ] **Step 4: Run — expect PASS** `3 passed`. **`test_d2_quantum_ci_yx_matches_paper_no_swap` is load-bearing**: it fails loudly if `a`/`c` are transposed.

- [ ] **Step 5: Commit**
```bash
git add src/quiverlab/resolutions_cs/resolution.py tests/resolutions_cs/test_leading.py
git commit -m "feat(cs): CS leading map δ_n (f_n, factors fixed) + low-degree d_1,d_2; swap-catcher test

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: The general differential — order-condition-pinned correction, collapsed matrices, the two gates

**Files:** Modify `resolution.py` (complete `_d_general`, add `matrix`, `dim_C`, `assert_dd_zero`, `assert_order_condition`, `_lower_generators`, `_idx_word`, `_solve`, `_matmul`); create `src/quiverlab/resolutions_cs/pelt.py`; `tests/resolutions_cs/test_differential.py`.

**Interfaces:** Consumes Task-5 leading map, `groebner.ReductionSystem` (`.order`), `fields.linalg.{rank, solve, nullspace}`. Produces `res.d_terms(n, chain)` for all `n` (CS Theorem 4.2: `δ_n + Σγ_i·gen_i`, `γ` solving `d_{n-1}∘d_n = 0`); `res.matrix(n, side)`; the two gates.

- [ ] **Step 1: Failing test** — Fixture A tower, Fixture C `d²=0` + byte-match to CS §6, both gates.
```python
import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.fields.linalg import rank
pytest.importorskip("quiverlab.groebner")


def _kx2(field=CC):
    Q = Quiver([1], {"x": (1, 1)})
    return ChouhySolotarResolution(Q.algebra(relations=["x*x"], field=field),
                                   build_reduction_system(Q, ["x*x"], field), max_degree=6)


def _qci(field=CC, xi="2"):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", f"y*x - ({xi})*x*y"]
    return ChouhySolotarResolution(Q.algebra(relations=rels, field=field),
                                   build_reduction_system(Q, rels, field), max_degree=8)


def test_kx2_boundary_ranks_char0():
    res = _kx2()
    assert [rank(res.matrix(n, "hom"), res.dom) for n in range(1, 7)] == [0, 1, 0, 1, 0, 1]


def test_kx2_dd_zero_both_sides():
    res = _kx2()
    res.assert_dd_zero(upto=6, side="hom"); res.assert_dd_zero(upto=6, side="coh")


def test_square_dd_zero():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    res = ChouhySolotarResolution(A, build_reduction_system(Q, ["a*b - c*d"], CC), max_degree=3)
    res.assert_dd_zero(upto=3, side="hom"); res.assert_dd_zero(upto=3, side="coh")


def test_qci_dd_zero_and_order_condition():
    for field in (CC, GF(2)):
        res = _qci(field=field)
        res.assert_dd_zero(upto=8, side="hom")               # non-commutative recursion closes
        res.assert_order_condition(upto=8)                   # CS Theorem 4.1 condition (2)


@pytest.mark.xfail(strict=False, reason="canonicalization pending: the d²=0 solve is unique only "
                   "up to the constraint nullspace, so a correct-but-noncanonical correction is "
                   "allowed; binding criteria (d²=0 + order gate + HH dims) are checked separately")
def test_qci_d3_correction_matches_paper():
    """CS §6 (verbatim, kept as an aspirational canonical-form pin): d_2^{CS}(1⊗y²x⊗1) =
    y⊗yx⊗1 + ξ 1⊗yx⊗y + ξ² x⊗y²⊗1 − 1⊗y²⊗x  (ξ=2). Edit #2: the BINDING criteria are
    test_qci_dd_zero_and_order_condition (d²=0 + order) and test_qci_homology_matches_bank_vector
    (HH dims, Task 7) — this exact-coefficient check is xfail-tolerant until the canonicalization
    stretch item (normal-form modulo the solve nullspace) flips it strict."""
    res = _qci()
    q = next(c for c in res.ss.S(3) if c.word == ("y", "y", "x"))
    got = {(res.to_int(c), a, t, cc) for (c, a, t, cc) in res.d_terms(3, q)}
    assert got == {
        ( 1, ("y",), ("y", "x"), ()),       #  y ⊗ yx ⊗ 1
        ( 2, (),     ("y", "x"), ("y",)),    #  ξ 1 ⊗ yx ⊗ y
        ( 4, ("x",), ("y", "y"), ()),        #  ξ² x ⊗ y² ⊗ 1
        (-1, (),     ("y", "y"), ("x",)),    # −1 ⊗ y² ⊗ x
    }


def test_cubic_tip_nonmonomial_raises_notimplemented():
    """RESTRICT boundary (edit #1): a non-quadratic (cubic tip) NON-monomial presentation
    raises NotImplementedError at the exact degree-≥3 differential. A = k<x,y>/(x²,y²,xyx−yxy)
    is finite-dimensional (basis {1,x,y,xy,yx,xyx}, dim 6) with a cubic tip and a nonzero tail.
    (If completion changes the basis, keep any admissible f.d. cubic-tip non-monomial algebra.)"""
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", "x*y*x - y*x*y"]
    A = Q.algebra(relations=rels, field=CC)
    res = ChouhySolotarResolution(A, build_reduction_system(Q, rels, CC), max_degree=4)
    with pytest.raises(NotImplementedError) as e:
        res.d_terms(3, res.ss.S(3)[0])                   # first degree-3 differential trips the guard
    assert "quadratic" in str(e.value).lower() and "right_decomposition" in str(e.value)
```

- [ ] **Step 2: Run — expect FAIL** `NotImplementedError: filled in Task 6`, then `AttributeError: matrix`.

- [ ] **Step 3: Implement.** First `src/quiverlab/resolutions_cs/pelt.py`:
```python
"""Uncollapsed P-element arithmetic for the correction linear-solve and the gates.
A P-ELEMENT (chain-degree n) is a dict (a_idx, chain_word, c_idx) -> coeff, where a_idx,c_idx
are A-basis indices (b,b' ∈ B), i.e. Σ coeff·(b ⊗ chain ⊗ b'). Built from a term list
(coeff, a_word, chain_word, c_word) by reducing a_word,c_word to normal-form A-vectors (π=β)
and expanding over B."""


def _resolve_chain(res, chain_word):
    if isinstance(chain_word, tuple) and len(chain_word) == 2 and chain_word[0] == "__v__":
        for c in res.ss.S(0):
            if c.o == chain_word[1]:
                return c
        raise KeyError(("vertex chain", chain_word))
    for n in range(res.ss.max_degree + 1):
        c = res._chain(n, chain_word)
        if c is not None:
            return c
    raise KeyError(("chain not found", chain_word))


def _accum(out, key, val, dom):
    cur = out.get(key)
    tot = val if cur is None else dom.add(cur, val)
    if dom.is_zero(tot):
        out.pop(key, None)
    else:
        out[key] = tot


def _vecs(res, ch, a_word, c_word):
    a = res.ar.vertex_vec(ch.o) if len(a_word) == 0 else res.ar.path_vec(a_word)
    c = res.ar.vertex_vec(ch.t) if len(c_word) == 0 else res.ar.path_vec(c_word)
    return a, c


def terms_to_pelt(res, term_list):
    dom, out = res.dom, {}
    for (coeff, a_word, chain_word, c_word) in term_list:
        ch = _resolve_chain(res, chain_word)
        a_vec, c_vec = _vecs(res, ch, a_word, c_word)
        for ai, av in enumerate(a_vec):
            if dom.is_zero(av):
                continue
            for ci, cv in enumerate(c_vec):
                if dom.is_zero(cv):
                    continue
                _accum(out, (ai, ch.word, ci), dom.mul(coeff, dom.mul(av, cv)), dom)
    return out


def apply_lower(res, n, pelt):
    """d_{n-1} applied to a P-element at chain-degree n-1, returning chain-degree n-2. For a
    key (b, σ, b') with d_{n-1}(1⊗σ⊗1) = Σ (c, a, τ, c'): image += c·(b·a ⊗ τ ⊗ c'·b')."""
    ar, dom, out = res.ar, res.dom, {}
    for (bi, chain_word, ci), coeff in pelt.items():
        ch = _resolve_chain(res, chain_word)
        bvec, cvec = ar.A._basis_vec(bi), ar.A._basis_vec(ci)
        for (c, a_word, tw, c_word) in res.d_terms(n - 1, ch):
            a_vec, c_vec = _vecs(res, ch, a_word, c_word)
            left, right = ar.mul(bvec, a_vec), ar.mul(c_vec, cvec)
            tch, base = _resolve_chain(res, tw), dom.mul(coeff, c)
            for ai, av in enumerate(left):
                if dom.is_zero(av):
                    continue
                for cj, cv in enumerate(right):
                    if dom.is_zero(cv):
                        continue
                    _accum(out, (ai, tch.word, cj), dom.mul(base, dom.mul(av, cv)), dom)
    return out
```
Then append to `resolution.py`:
```python
    def _idx_word(self, idx):
        tag, val = self.ar.basis_word(idx)
        return () if tag == "v" else val

    def _lower_generators(self, n, chain):
        """Basis of ⟨\\overline{𝓛}_{n-1}^{≺}(σ)⟩: single-term lists (1, b, p.word, b') with
        p∈S_{n-1}, b,b'∈B, and the loop path b·p·b' strictly ≺ σ (parallel to σ). Finite by DCC."""
        order, skey, one = self.rs.order, self.rs.order.key(chain.word), self._one()
        gens = []
        for p in self.ss.S(n - 1):
            for bi in self.ar.corner(chain.o, p.o, "coh"):        # b : o(σ)->o(p)
                bw = self._idx_word(bi)
                for ci in self.ar.corner(p.t, chain.t, "coh"):    # b' : t(p)->t(σ)
                    cw = self._idx_word(ci)
                    if order.key(bw + p.word + cw) < skey:
                        gens.append((one, bw, p.word, cw))
        return gens

    def _require_in_scope(self):
        """RESTRICT (edit #1): the uncollapsed leading map is provably CS's δ only for
        quadratic tips (v_n = u_0; CS Prop, TeX 772) or monomial presentations (no
        correction; collapsed map = validated Bardzell). A non-quadratic non-monomial
        presentation raises NotImplementedError at this exact boundary (spec §6)."""
        quadratic = all(len(w) <= 2 for w in self.rs.leading_words())
        monomial = all(len(r.tail) == 0 for r in self.rs.rules)
        if not quadratic and not monomial:
            raise NotImplementedError(
                "CS is certified for quadratic tips (S ⊆ Q_2; CS Prop, arXiv:1406.2300 TeX "
                "line 772) or monomial presentations; this presentation has a tip of length "
                ">= 3 AND a nonzero tail (non-quadratic non-monomial). There v_n ≠ u_0 and the "
                "uncollapsed leading map feeding the correction solve is not provably δ_CS. "
                "Lifting this needs the right_decomposition upgrade (a Plan-04 stretch item); "
                "boundary per spec §6 risk register.")

    def _d_general(self, n, chain):
        from quiverlab.resolutions_cs.pelt import terms_to_pelt, apply_lower
        self._require_in_scope()                         # NotImplementedError at the exact boundary
        key = (n, chain.word)
        if key in self._d_cache:
            return self._d_cache[key]
        delta = self.delta_terms(n, chain)
        gens = self._lower_generators(n, chain)
        if not gens:
            self._d_cache[key] = delta
            return delta
        dom = self.dom
        rhs_pe = apply_lower(self, n, terms_to_pelt(self, delta))
        cols = [apply_lower(self, n, terms_to_pelt(self, [g])) for g in gens]
        keys = sorted(set(rhs_pe) | {k for col in cols for k in col})
        M = [[col.get(k, dom.zero()) for col in cols] for k in keys]
        rhs = [dom.neg(rhs_pe.get(k, dom.zero())) for k in keys]
        gamma = self._solve(M, rhs, len(gens))
        if gamma is None:
            raise NotImplementedError(
                f"CS correction linear system is inconsistent at degree {n}, chain "
                f"{chain.word}: this admissible algebra needs the higher CS homotopy "
                f"correction, outside quiverlab v1's construction (spec §6 risk register)")
        terms = list(delta)
        for coeff, (c1, a, tw, cc) in zip(gamma, gens):
            if not dom.is_zero(coeff):
                terms.append((dom.mul(coeff, c1), a, tw, cc))
        self._d_cache[key] = terms
        return terms

    def _solve(self, M, rhs, ncols):
        """Return γ (any solution of M·γ = rhs over the domain) or None if inconsistent.
        Uses fields.linalg.solve; if that requires square/consistent input, fall back to a
        particular solution (augmented RREF) + nullspace freedom (exact, over the domain)."""
        from quiverlab.fields.linalg import solve
        if not M:                                        # no equations -> any γ; choose zeros
            return [self.dom.zero()] * ncols
        return solve(M, rhs, self.dom)                   # returns a solution or None (inconsistent)

    def _basis(self, n, side):
        return [(ch, j) for ch in self.ss.S(n) for j in self.ar.corner(ch.o, ch.t, side)]

    def dim_C(self, n, side):
        return len(self._basis(n, side))

    def matrix(self, n, side):
        from quiverlab.resolutions_cs.pelt import _resolve_chain, _vecs
        dom = self.dom
        if side == "hom":
            rows, cols = self._basis(n - 1, "hom"), self._basis(n, "hom")
        else:
            rows, cols = self._basis(n + 1, "coh"), self._basis(n, "coh")
        ridx = {(ch.word, j): i for i, (ch, j) in enumerate(rows)}
        M = [[dom.zero()] * len(cols) for _ in range(len(rows))]
        if side == "hom":
            for cj, (sigma, j) in enumerate(cols):
                ej = self.ar.A._basis_vec(j)
                for (coeff, a_word, tw, c_word) in self.d_terms(n, sigma):
                    a_vec, c_vec = _vecs(self, sigma, a_word, c_word)
                    val = self.ar.mul(c_vec, self.ar.mul(ej, a_vec))     # b·w·a  (homology collapse)
                    tw_word = _resolve_chain(self, tw).word
                    for p, vp in enumerate(val):
                        if not dom.is_zero(vp) and (tw_word, p) in ridx:
                            r = ridx[(tw_word, p)]
                            M[r][cj] = dom.add(M[r][cj], dom.mul(coeff, vp))
        else:
            for cj, (sigma, j) in enumerate(cols):                       # δ^n: C^n -> C^{n+1}
                ej = self.ar.A._basis_vec(j)
                for tau in self.ss.S(n + 1):
                    for (coeff, a_word, tw, c_word) in self.d_terms(n + 1, tau):
                        if _resolve_chain(self, tw).word != sigma.word:
                            continue
                        a_vec, c_vec = _vecs(self, tau, a_word, c_word)
                        val = self.ar.mul(a_vec, self.ar.mul(ej, c_vec)) # a·w·b  (cohomology collapse)
                        for p, vp in enumerate(val):
                            if not dom.is_zero(vp) and (tau.word, p) in ridx:
                                r = ridx[(tau.word, p)]
                                M[r][cj] = dom.add(M[r][cj], dom.mul(coeff, vp))
        return M

    def _matmul(self, A_, B_):
        dom = self.dom
        if not A_ or not B_:
            return []
        inner, rows_out, cols_out = len(B_), len(A_), len(B_[0])
        out = [[dom.zero()] * cols_out for _ in range(rows_out)]
        for i in range(rows_out):
            for k in range(inner):
                aik = A_[i][k]
                if dom.is_zero(aik):
                    continue
                for j in range(cols_out):
                    out[i][j] = dom.add(out[i][j], dom.mul(aik, B_[k][j]))
        return out

    def assert_dd_zero(self, upto, side):
        dom = self.dom
        for n in range(2, upto + 1):
            if side == "hom":
                prod = self._matmul(self.matrix(n - 1, "hom"), self.matrix(n, "hom"))
            else:
                prod = self._matmul(self.matrix(n, "coh"), self.matrix(n - 1, "coh"))
            if any(not dom.is_zero(x) for row in prod for x in row):
                raise AssertionError(
                    f"CS d²≠0 at degree {n} (side={side}); the correction solve failed to "
                    "close — a bug, never an approximation")

    def assert_order_condition(self, upto):
        """CS Theorem 4.1 condition (2): every correction term's loop path b·p·b' is strictly ≺ σ."""
        order = self.rs.order
        for n in range(3, upto + 1):
            for sigma in self.ss.S(n):
                lead = {(a, tw, cc): c for (c, a, tw, cc) in self.delta_terms(n, sigma)}
                for (c, a, tw, cc) in self.d_terms(n, sigma):
                    if lead.get((a, tw, cc)) == c:
                        continue                                 # a leading term
                    if not order.key(a + tw + cc) < order.key(sigma.word):
                        raise AssertionError(
                            f"order condition violated at degree {n}, chain {sigma.word}: "
                            f"correction {(a, tw, cc)} is not ≺ σ")
```

- [ ] **Step 4: Run — expect `5 passed, 1 xfailed`.** The **binding** criteria are `test_qci_dd_zero_and_order_condition` (both gates) and the HH-dim tests (Task 7/9/11); `test_qci_d3_correction_matches_paper` is `xfail(strict=False)` because the `d²=0` solve is unique only modulo the constraint nullspace, so a correct-but-noncanonical correction must not fail the build (edit #2). `test_cubic_tip_nonmonomial_raises_notimplemented` pins the RESTRICT boundary (edit #1). If `test_kx2_boundary_ranks_char0` fails ⟹ leading `δ` parity/sign wrong (Task 5). If Plan-01 `solve` cannot return non-square/underdetermined solutions or signal inconsistency, add the exact augmented-RREF particular-solution helper (stated in `_solve`) in this task.
  > **Stretch item (execution-time): canonicalization.** Add a deterministic normal form for the correction — reduce `γ` modulo the nullspace of the `d²=0` constraint matrix (e.g. RREF-pivot the augmented system so free variables are set to a canonical value), yielding CS §6's verbatim coefficients. When implemented, flip `test_qci_d3_correction_matches_paper` and the Task-11 byte-oracle `xfail` markers to `strict=True`.

- [ ] **Step 5: Commit**
```bash
git add src/quiverlab/resolutions_cs/resolution.py src/quiverlab/resolutions_cs/pelt.py tests/resolutions_cs/test_differential.py
git commit -m "feat(cs): general differential (order-condition-pinned correction) + d²=0 & order gates

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: HH•/HH• dims + bases over any Domain; Resolution facade; admissibility boundary

**Files:** Create `src/quiverlab/resolutions_cs/homology.py`, `engine_facade.py`, `build.py`, `_fieldshim.py`; `tests/resolutions_cs/test_homology.py`.

**Interfaces:** Consumes `res.matrix`, `fields.linalg.rank`, `hochschild.table.HHTable`, `engine.resolutions.Resolution`, `groebner.build_reduction_system`. Produces `reduction_system_of(A)`, `cs_cohomology_dims/cs_homology_dims -> HHTable`, `cs_hh_basis`, `engine_facade.CSResolution(Resolution)`, the admissibility gate.

- [ ] **Step 1: Failing test** `tests/resolutions_cs/test_homology.py`:
```python
import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
pytest.importorskip("quiverlab.groebner")


def _A(field=CC, rels=("x*x",), arrows=None, verts=(1,)):
    return Quiver(list(verts), arrows or {"x": (1, 1)}).algebra(relations=list(rels), field=field)


def test_kx2_dims_char0_and_char2():
    assert cs_cohomology_dims(_A(), 6).dims == [2, 1, 1, 1, 1, 1, 1]
    assert cs_homology_dims(_A(), 6).dims == [2, 1, 1, 1, 1, 1, 1]
    assert cs_cohomology_dims(_A(field=GF(2)), 5).dims == [2, 2, 2, 2, 2, 2]


def test_square_dims():
    A = _A(rels=["a*b - c*d"], verts=(1, 2, 3, 4),
           arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    assert cs_cohomology_dims(A, 4).dims == [1, 0, 0, 0, 0]
    assert cs_homology_dims(A, 4).dims == [4, 0, 0, 0, 0]


def test_qci_homology_matches_bank_vector():
    A = _A(rels=["x*x", "y*y", "y*x - 2*x*y"], arrows={"x": (1, 1), "y": (1, 1)})
    assert cs_homology_dims(A, 12).dims == [3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]


def test_engine_facade_is_resolution_protocol():
    from quiverlab.resolutions_cs.engine_facade import CSResolution
    from quiverlab.engine.resolutions import Resolution
    from quiverlab.engine.adapter import to_engine
    Ap = _A(field=GF(32003), rels=["a*b - c*d"], verts=(1, 2, 3, 4),
            arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    R = CSResolution(Ap)
    assert isinstance(R, Resolution)
    E = to_engine(Ap.unit_adapted())
    b2 = R.term_basis(E, 2)
    assert R.differential_matrix(E, 2, b2, {g: i for i, g in enumerate(R.term_basis(E, 1))}).shape[1] == len(b2)
```

- [ ] **Step 2: Run — expect FAIL** `ModuleNotFoundError: …homology`.

- [ ] **Step 3: Implement.** `build.py`/`_fieldshim.py` (`reduction_system_of`), `homology.py` (dims via `dim − rank d_n − rank d_{n±1}`, gates wired in), `engine_facade.py` (int64 wrapper delegating to `ChouhySolotarResolution` over `PrimeField`, packing `numpy.int64` matrices, un-reduced).
```python
# homology.py
def _require_admissible(rs):
    from quiverlab.errors import AdmissibilityError
    if not rs.is_confluent or not rs.irreducibles:
        raise AdmissibilityError("CS runs only on a certified-admissible reduction system",
                                 hint="Groebner completion did not certify confluence / a finite basis")


def cs_cohomology_dims(A, top, max_cells=4_000_000):
    from quiverlab.fields.linalg import rank
    from quiverlab.hochschild.table import HHTable
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    rs = reduction_system_of(A); _require_admissible(rs)
    res = ChouhySolotarResolution(A, rs, max_degree=top + 1, max_cells=max_cells)
    res.assert_dd_zero(upto=top + 1, side="coh"); res.assert_order_condition(upto=top + 1)
    dom = A.domain
    r = [rank(res.matrix(n, "coh"), dom) for n in range(top + 1)]
    dims = [res.dim_C(n, "coh") - r[n] - (r[n - 1] if n else 0) for n in range(top + 1)]
    return HHTable(dims, "HH^", repr(A).splitlines()[0], engine="Chouhy-Solotar")
```
`cs_homology_dims` mirrors with `side="hom"` and `dim − rank b_n − rank b_{n+1}` (`b_n = matrix(n,"hom")`, `b_{n+1} = matrix(n+1,"hom")`; `b_0 = 0`). `build.reduction_system_of`:
```python
def reduction_system_of(A):
    from quiverlab.groebner import build_reduction_system
    from quiverlab.resolutions_cs._fieldshim import field_for_domain
    if A.quiver is None or A.relations is None:
        raise ValueError("CS needs an algebra built by Quiver.algebra (quiver+relations present)")
    return build_reduction_system(A.quiver, list(A.relations), field_for_domain(A.domain))
```
`_fieldshim.py`:
```python
class _DomainField:
    def __init__(self, dom):
        self._dom = dom
    def parse_entry(self, x):
        return x
    def make_domain(self, entries):
        return self._dom

def field_for_domain(dom):
    return _DomainField(dom)
```
> `reduction_system_of` re-runs Plan-03's public `build_reduction_system` over `A.domain` (there is no `A.reduction_system()` — Pillar-4). The `_DomainField` shim passes the stored relation coefficients straight into the existing domain (`build_reduction_system` calls `field.parse_entry`/`make_domain`; both become identity/`return dom`). If `build_reduction_system`'s `field` contract needs more, add exactly the methods it calls (grep the committed Plan-03 `system.py`). `cs_hh_basis(A, n, side)` returns representative (co)cycles via `fields.linalg.nullspace` of the relevant matrix modulo image.

- [ ] **Step 4: Run — expect PASS** `4 passed`. `test_qci_homology_matches_bank_vector` pins `[3,2,…,2]` over ℚ to degree 12 (deep recursion).

- [ ] **Step 5: Commit**
```bash
git add src/quiverlab/resolutions_cs/homology.py src/quiverlab/resolutions_cs/engine_facade.py src/quiverlab/resolutions_cs/build.py src/quiverlab/resolutions_cs/_fieldshim.py tests/resolutions_cs/test_homology.py
git commit -m "feat(cs): HH•/HH• over any domain (gate-enforced), Resolution facade, admissibility gate

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Validation battery (a) — CS ≡ Bardzell exactly on the monomial zoo

**Files:** Create `tests/resolutions_cs/test_battery_bardzell.py`.

**Interfaces:** Consumes `CSResolution` (facade), `engine.resolutions_bardzell.BardzellResolution`, the monomial builders.

- [ ] **Step 1: Failing test** — term counts, ranks, HH dims identical (GF(32003)); include `cyclic_nakayama(5, 3)` (genuine overlaps).
```python
import pytest, numpy as np
from quiverlab.engine.resolutions_bardzell import BardzellResolution, MonomialPresentation
from quiverlab.engine.hh_engine import hochschild_homology_dims
from quiverlab.engine.adapter import to_engine
from quiverlab.resolutions_cs.engine_facade import CSResolution
pytest.importorskip("quiverlab.groebner")

CASES = [("kx", 3, 12), ("cyclic_nakayama", (5, 3), 10), ("local_radsq", 3, 10)]


def _build(name, param):
    ...  # (core.Algebra over GF(32003), matching MonomialPresentation) per case


@pytest.mark.parametrize("name,param,N", CASES)
def test_cs_equals_bardzell_terms_and_ranks(name, param, N):
    A, pres = _build(name, param)
    E = to_engine(A.unit_adapted())
    cs, bd = CSResolution(A), BardzellResolution(pres)
    for n in range(N + 1):
        assert len(cs.term_basis(E, n)) == len(bd.term_basis(E, n))
    for n in range(1, N + 1):
        Mc = cs.differential_matrix(E, n, cs.term_basis(E, n),
                                    {g: i for i, g in enumerate(cs.term_basis(E, n - 1))})
        Mb = bd.differential_matrix(E, n, bd.term_basis(E, n),
                                    {g: i for i, g in enumerate(bd.term_basis(E, n - 1))})
        assert np.linalg.matrix_rank(Mc % 32003) == np.linalg.matrix_rank(Mb % 32003)


@pytest.mark.parametrize("name,param,N", CASES)
def test_cs_equals_bardzell_hh(name, param, N):
    A, pres = _build(name, param)
    E = to_engine(A.unit_adapted())
    hb = hochschild_homology_dims(E, N, primes=(32003,), resolution=BardzellResolution(pres))[32003]
    hc = hochschild_homology_dims(E, N, primes=(32003,), resolution=CSResolution(A))[32003]
    assert list(hc) == list(hb)
```

- [ ] **Step 2: Run — expect FAIL** (`_build` undefined). Any real CS≠Bardzell discrepancy is a **Task-5/6 bug** (fix `resolution.py`, re-run 5–8), never a test tweak.

- [ ] **Step 3: Implement** `_build`.

- [ ] **Step 4: Run — expect PASS** `6 passed` (spec-§6 "CS specializes to Bardzell").

- [ ] **Step 5: Commit**
```bash
git add tests/resolutions_cs/test_battery_bardzell.py
git commit -m "test(cs): battery (a) — CS ≡ Bardzell exactly (terms, ranks, HH) on the monomial zoo

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: Validation battery (b) — CS ≡ bar ≡ minimal degreewise on non-monomial algebras

**Files:** Create `tests/resolutions_cs/test_battery_bar.py`.

- [ ] **Step 1: Failing test** — square + quantum CI, over ℚ, ℚ(i), GF(p), within the bar window.
```python
import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
from quiverlab.hochschild.bar import hochschild_cohomology_dims, hochschild_homology_dims
pytest.importorskip("quiverlab.groebner")


def _square(field):
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=field)


def _qci(field, xi="2"):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    return Q.algebra(relations=["x*x", "y*y", f"y*x - ({xi})*x*y"], field=field)


@pytest.mark.parametrize("field", [CC, GF(2), GF(3), GF(5)])
def test_cs_equals_bar_square(field):
    A = _square(field)
    assert cs_cohomology_dims(A, 4).dims == hochschild_cohomology_dims(A, 4).dims
    assert cs_homology_dims(A, 4).dims == hochschild_homology_dims(A, 4).dims


@pytest.mark.parametrize("xi", ["2", "i", "-1"])
def test_cs_equals_bar_qci(xi):
    A = _qci(CC, xi)                      # xi="i" exercises the ℚ(i) number field
    assert cs_homology_dims(A, 3).dims == hochschild_homology_dims(A, 3).dims
    assert cs_cohomology_dims(A, 3).dims == hochschild_cohomology_dims(A, 3).dims
```

- [ ] **Step 2: Run — expect FAIL** initially. Mismatches ⟹ Task-6 bug (never loosen the oracle).

- [ ] **Step 3: Implement** — no `src/` change.

- [ ] **Step 4: Run — expect PASS** (`7 passed`).

- [ ] **Step 5: Commit**
```bash
git add tests/resolutions_cs/test_battery_bar.py
git commit -m "test(cs): battery (b) — CS ≡ bar degreewise on non-monomial algebras (ℚ, ℚ(i), GF(p))

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 10: Validation battery (c) — literature oracles (real values)

**Files:** Create `tests/resolutions_cs/test_battery_literature.py`.

- [ ] **Step 1: Failing test** — cited values; BGMS = the bank-derived `[3,2,2,…]`.
```python
import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
from quiverlab.hochschild.bar import hochschild_cohomology_dims
pytest.importorskip("quiverlab.groebner")


def test_kxn_across_characteristics():
    for n in (2, 3, 4, 5):                                # HH_0 = n; HH_i = 1 (i≥1) char 0 (p∤n)
        A = Quiver([1], {"x": (1, 1)}).algebra(relations=[f"x^{n}"], field=CC)
        d = cs_homology_dims(A, 6).dims
        assert d[0] == n and all(v == 1 for v in d[1:])
    A5 = Quiver([1], {"x": (1, 1)}).algebra(relations=["x^5"], field=GF(5))
    assert cs_homology_dims(A5, 6).dims == [5, 5, 5, 5, 5, 5, 5]     # p | n pathology


def test_bgms_quantum_ci_homology():
    # BGMS quantum CI k<x,y>/(x²,y²,yx−ξxy), ξ=2 (char 0) -> HH_• = [3,2,2,...] (bank-derived).
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - 2*x*y"], field=CC)
    assert cs_homology_dims(A, 8).dims == [3, 2, 2, 2, 2, 2, 2, 2, 2]


def test_hereditary_dynkin_happel():
    # Happel: hereditary ⟹ HH^i = 0 (i≥2). Linear A_3 (no relations).
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=CC)
    d = cs_cohomology_dims(A, 4).dims
    assert d[0] == 1 and d[2:] == [0, 0, 0]


def test_gentle_case_matches_bar():
    A = _build_gentle()                                  # small gentle algebra (cite source in-file)
    assert cs_cohomology_dims(A, 3).dims == hochschild_cohomology_dims(A, 3).dims
```

- [ ] **Step 2: Run — expect FAIL** (`_build_gentle` undefined). Fix the gentle builder; keep the citation beside each value. **Do not invent values**; where no published vector is at hand, cross-check bar and cite "bar cross-check".

- [ ] **Step 3: Implement** the builders/constants with citations.

- [ ] **Step 4: Run — expect PASS** (`4 passed`).

- [ ] **Step 5: Commit**
```bash
git add tests/resolutions_cs/test_battery_literature.py
git commit -m "test(cs): battery (c) — literature oracles (k[x]/x^n, BGMS [3,2,…], Happel, gentle)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 11: Validation battery (d) — the bank closed-form CS families as byte-level oracles

**Files:** Create `tests/resolutions_cs/test_battery_bank_oracle.py`.

**Interfaces:** Consumes the READ-ONLY bank `hanlab/resolutions_cs.py` (imported by path), the `CSResolution` facade.

- [ ] **Step 1: Failing test** — entry-by-entry equality of the collapsed **homology** matrices (bank = homology only) with the bank closed forms, both families, several primes.
```python
import importlib.util, pathlib, sys
import numpy as np, pytest
from quiverlab import Quiver, GF
pytest.importorskip("quiverlab.groebner")
BANK = pathlib.Path("/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/HansConjecture")


def _bank_cs():
    sys.path.insert(0, str(BANK / "hanlab"))
    try:
        spec = importlib.util.spec_from_file_location("bank_cs", BANK / "hanlab/resolutions_cs.py")
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
    finally:
        sys.path.pop(0)


# BINDING criterion: HH_* dims computed on the SAME algebra through the bank CS resolution vs the
# Plan-04 CSResolution must agree (invariant under the correction's nullspace non-uniqueness).
@pytest.mark.parametrize("a", [2, 3, 4])
def test_cs_hh_matches_bank_truncpoly(a):
    ...  # k[x]/(x^a): hochschild_homology_dims(E, N, resolution=CSResolution) == via bank ChouhySolotarResolution


@pytest.mark.parametrize("xi", [1, 2, 3])
def test_cs_hh_matches_bank_quantum_ci(xi):
    ...  # k<x,y>/(x²,y²,yx−ξxy): same, HH_* dims equal (both resolutions)


# ASPIRATIONAL byte-exact pins (edit #2): the collapsed matrices agree ENTRY-BY-ENTRY once the
# correction is canonicalized. xfail-tolerant until the canonicalization stretch item lands.
@pytest.mark.xfail(strict=False, reason="canonicalization pending (correction unique mod nullspace)")
@pytest.mark.parametrize("a", [2, 3, 4])
def test_cs_matches_bank_truncpoly_bytes(a):
    ...  # assert ((M_cs - M_bank) % p == 0).all() for n=1..N (aligned generator order)


@pytest.mark.xfail(strict=False, reason="canonicalization pending (correction unique mod nullspace)")
@pytest.mark.parametrize("xi", [1, 2, 3])
def test_cs_matches_bank_quantum_ci_bytes(xi):
    ...  # align (s,t)↔y^s x^t; assert entrywise equality mod p
```
Helper `assert_matrices_equal_mod_p(build_ql, bank_rs, N, p)` reads the bank `ChouhySolotarResolution.differential_matrix(alg, n, basis_n, index_nm1)` (n=1..N), builds the Plan-04 `CSResolution` collapsed homology matrix in the **same generator order** (align `Chain.word` ↔ bank labels: `("c",)/("v",)` for truncpoly, `(s,t)` for qci), and asserts equality mod `p`. The HH-dim helper compares `hochschild_homology_dims(E, N, resolution=…)[p]` for both resolutions.

- [ ] **Step 2: Run — expect FAIL** (helpers undefined; then a real **HH-dim** disagreement ⟹ Task-6 bug — the byte-exact pins may `xfail` under non-canonical corrections but the HH-dim equality is binding and swap-sensitive).

- [ ] **Step 3: Implement** the HH-dim comparison + the byte comparison with generator-order alignment.

- [ ] **Step 4: Run — expect `6 passed, 6 xfailed`** (HH-dim binding tests pass; the entrywise pins are xfail-tolerant until canonicalization flips them strict — edit #2).

- [ ] **Step 5: Commit**
```bash
git add tests/resolutions_cs/test_battery_bank_oracle.py
git commit -m "test(cs): battery (d) — bank resolutions_cs closed forms as byte-level oracles

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 12: Comparison morphisms CS↔bar and transport of cup/cap/bracket classes

**Files:** Create `src/quiverlab/resolutions_cs/comparison.py`; `tests/resolutions_cs/test_comparison.py`.

**Interfaces:** Consumes `ChouhySolotarResolution`, `hochschild.bar`, `engine.tt_calculus` (via the GF(p) facade), `fields.linalg.solve`. Produces `Comparison(A)` with `.Phi(n)`, `.Psi(n)`, `transport_cocycle_cs_to_bar`, `transport_class_bar_to_cs`, `cup_of_cs_classes`, `bracket_of_cs_classes`, `same_cohomology_class`, `assert_chain_map`, `assert_transport_roundtrip_identity`; window-bounded (`NotImplementedError` past the bar window).

- [ ] **Step 1: Failing test** — chain-map law, round-trip identity on `HH`, transported-cup consistency, window boundary.
```python
import pytest
from quiverlab import Quiver, GF
from quiverlab.resolutions_cs.comparison import Comparison
pytest.importorskip("quiverlab.groebner")


def _kx2_gf():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(32003))


def _square_gf():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=GF(32003))


def test_phi_is_chain_map():
    Comparison(_square_gf()).assert_chain_map(upto=2)                 # b_bar ∘ Φ = Φ ∘ d_cs


def test_transport_roundtrip_identity_on_cohomology():
    Comparison(_kx2_gf()).assert_transport_roundtrip_identity(upto=3)  # Ψ*Φ* = id on HH^n


def test_transported_cup_consistent():
    comp = Comparison(_kx2_gf())
    u = comp.hh_class_cs(1, 0)
    assert comp.same_cohomology_class(comp.cup_of_cs_classes(u, u),
                                      comp.transport_then_bar_cup(u, u), degree=2)


def test_operation_window_boundary():
    comp = Comparison(_kx2_gf())
    with pytest.raises(NotImplementedError):
        comp.cup_of_cs_classes(comp.hh_class_cs(20, 0), comp.hh_class_cs(20, 0))
```

- [ ] **Step 2: Run — expect FAIL** `ModuleNotFoundError: …comparison`.

- [ ] **Step 3: Implement.** `Φ_n(1⊗σ⊗1) = 1 ⊗ u_0 ⊗ … ⊗ u_{n-1} ⊗ 1` on the block decomposition (CS `i_n∘ι` into the reduced bar); `Ψ_n` the CS projection. `assert_chain_map` checks `b^{bar}Φ = Φ d^{cs}` and `d^{cs}Ψ = Ψ b^{bar}` as matrix identities. `transport_cocycle_cs_to_bar` = pull through the transpose; `same_cohomology_class` = equal modulo the coboundary image (`fields.linalg.solve`). `cup_of_cs_classes(u,v) = Ψ*( Φ*(u) ⌣_{bar} Φ*(v) )` via `tt_calculus`; **boundary** — if either degree exceeds the bar-buildable window (bar `max_cells`), raise `NotImplementedError("cup/bracket transport is delivered only within the bar-comparison window; native CS cup is a later phase")`.

- [ ] **Step 4: Run — expect PASS** (`4 passed`).

- [ ] **Step 5: Commit**
```bash
git add src/quiverlab/resolutions_cs/comparison.py tests/resolutions_cs/test_comparison.py
git commit -m "feat(cs): CS↔bar comparison maps + windowed transport of cup/cap/bracket classes

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 13: Dispatch — explicit `engine="cs"` + opt-in `auto_cs` (shipped `auto` UNCHANGED)

**Files:** Modify `src/quiverlab/core/algebra.py`; `tests/resolutions_cs/test_dispatch.py`.

**Interfaces:** Produces `A.hochschild_cohomology(n, engine="cs")` / `hochschild_homology(n, engine="cs")`; an **opt-in** `auto_cs=False` keyword so `engine="auto", auto_cs=True` routes non-monomial admissible algebras to CS, while default `engine="auto"` keeps its Plan-02 behaviour **exactly** (Pillar-4: no silent redefinition). Trace list populated; `HHTable.engine == "Chouhy-Solotar"`.

- [ ] **Step 1: Failing test** `tests/resolutions_cs/test_dispatch.py`:
```python
import pytest
from quiverlab import Quiver, CC
pytest.importorskip("quiverlab.groebner")


def _square():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=CC)


def test_engine_cs_selectable():
    t = _square().hochschild_cohomology(4, engine="cs")
    assert t.dims == [1, 0, 0, 0, 0] and "Chouhy-Solotar" in (t.engine or "")


def test_default_auto_is_unchanged():
    """Compatibility: default engine='auto' gives the SAME result & label as before Plan 04."""
    t = _square().hochschild_cohomology(4)               # engine="auto", auto_cs default False
    assert "Chouhy-Solotar" not in (t.engine or "") and t.dims == [1, 0, 0, 0, 0]


def test_opt_in_auto_cs_routes_to_cs():
    t = _square().hochschild_cohomology(4, auto_cs=True)
    assert "Chouhy-Solotar" in (t.engine or "")


def test_trace_claims_equal_values():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=CC)
    A.hochschild_cohomology(4, engine="cs", trace=(tr := []))
    dim_events = [e for e in tr if type(e).__name__ == "ResolutionTerm"]
    assert [e.collapsed_dim for e in dim_events][:3] == [2, 2, 2]
```

- [ ] **Step 2: Run — expect FAIL** `engine 'cs'` unknown (`QuiverlabError`).

- [ ] **Step 3: Implement.** Extend the `engine` set to `("auto","bar","fast","cs")`; add the CS branch (`cs_cohomology_dims`/`cs_homology_dims`); add `auto_cs=False`, `trace=None`. Routing: `engine="cs"` → CS; `engine="auto"` → the **existing Plan-02 dispatch verbatim** unless `auto_cs=True`, in which case a non-monomial admissible algebra routes to CS. Populate the passed `trace` list from resolution events. Keep internals out of `__all__`.

- [ ] **Step 4: Run — expect PASS** `4 passed`, then full suite + float gate green:
`cd /Users/marco/Desktop/HomologicalNetworks/quiverlab && NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q && NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/test_no_floats.py -q`

- [ ] **Step 5: Commit**
```bash
git add src/quiverlab/core/algebra.py tests/resolutions_cs/test_dispatch.py
git commit -m "feat(cs): engine='cs' path + opt-in auto_cs flag (shipped 'auto' unchanged), trace hooks

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 14: Acceptance — deep non-monomial HH impossible for bar + internals chapter 08 + Plan-05 freeze

**Files:** Create `tests/resolutions_cs/test_acceptance.py`, `docs/internals/08-chouhy-solotar.md`.

- [ ] **Step 1: Acceptance test.** **Exact algebra:** quantum CI `A = k⟨x,y⟩/(x², y², yx − 2·xy)` over **ℚ** (CS §6 Happel-counterexample family). **Depth:** `HH_0…HH_{12}`. **Expected (provenance):** bank-derived `[3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]` (Task-11 byte-oracle confirms the differentials; Task-9 confirms bar-agreement in the low window). **Impossible for bar:** bar `dim C_n = 4·3^n`, so the degree-12 boundary matrix has `~4·3^{11} × ~4·3^{12} ≈ 5·10^{11}` entries `≫ max_cells` and bar raises `DepthLimitError`; CS has `dim C_n = n+1`.
```python
import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.errors import DepthLimitError
pytest.importorskip("quiverlab.groebner")

BGMS_QCI_XI2 = [3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]     # bank-derived, degree 0..12


def _qci(field):
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - 2*x*y"], field=field)


def test_deep_nonmonomial_hh_impossible_for_bar():
    assert _qci(CC).hochschild_homology(12, engine="cs").dims == BGMS_QCI_XI2
    with pytest.raises(DepthLimitError):
        _qci(CC).hochschild_homology(12, engine="bar")     # exponential blow-up, loud


def test_characteristic_sweep_shows_pathologies():
    # Bank-derived char-dependence (a headline for the maths paper):
    assert _qci(GF(3)).hochschild_homology(8, engine="cs").dims == [3, 4, 6, 8, 10, 12, 14, 16, 18]
    assert _qci(GF(2)).hochschild_homology(8, engine="cs").dims == [3, 4, 4, 4, 4, 4, 4, 4, 4]
```

- [ ] **Step 2: Run — expect FAIL then PASS.** Every value is independently pinned (Task 11 HH-dim binding oracle per degree; Task 9 bar window; the char-sweep vectors are the bank-run values recorded here). No "fix the constant" placeholder.
  > **Execution-time note (edit #4):** before finalising, re-run the read-only bank oracle to reconfirm the char-sweep vectors on the machine of record — in particular the `p=5` vector `[3,2,3,4,3,2,4,6,4,2,5,8,5]` (recorded in Fixture C; a period-4 pathology because `2` has multiplicative order 4 mod 5). Add `assert _qci(GF(5)).hochschild_homology(12, engine="cs").dims == [...]` once reconfirmed. Command: `cd $BANK/hanlab && NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 <venv>/python -c "from resolutions_cs import qci_reduction_system, ChouhySolotarResolution; from reduction_algebra import algebra_from_reduction_system; from hh_engine import hochschild_homology_dims; rs=qci_reduction_system(2); A=algebra_from_reduction_system(rs,maxlen=6); print(hochschild_homology_dims(A,12,primes=(5,),resolution=ChouhySolotarResolution(rs,A))[5])"`.

- [ ] **Step 3: Write the internals chapter** `docs/internals/08-chouhy-solotar.md` (Marco's standing rule; the first `docs/internals/` chapter, establishing the format for a low-coding algebraist; ~120–180 lines). Complete content:
```markdown
# 08 — The Chouhy–Solotar resolution

## What this computes
Given `A = kQ/I` (any admissible presentation, any exact field), the CS engine builds a
small projective resolution of `A` as an `A`-bimodule and reads Hochschild (co)homology
off it — reaching degrees the bar complex never can.

## The objects (how they are represented in code)
- **Reduction system `rs`** (`quiverlab.groebner.build_reduction_system`): the confluent
  rewriting rules `s → f_s`. `rs.leading_words()` are the tips `S`; `rs.irreducibles` a
  basis of `A`. Words are tuples of arrow NAMES, left to right.
- **Ambiguity chain `Chain`** (`resolutions_cs.terms`): an element of `S_n`. `word` is the
  path; `blocks = (u_0,…,u_{n-1})` its unique CS decomposition; `o,t` its endpoints.
  `S_0`=vertices, `S_1`=arrows, `S_2`=tips, `S_n (n≥3)` = `(n-1)`-fold overlaps of tips.
  Computed by reusing Bardzell's `associated_paths` on the tip monomial algebra `A_S`.
- **Resolution term `P_n = ⨁_{σ∈S_n} A e_{o(σ)} ⊗ e_{t(σ)} A`.** Tensoring/homming down,
  each `σ` contributes the corner `e_tAe_o` (homology) or `e_oAe_t` (cohomology).
- **A term** `(coeff, a, τ, c)` means `coeff · (a ⊗ τ ⊗ c)`; `a,c` are paths, `τ ∈ S_{n-1}`.

## The computation, step by step (on k[x]/(x²))
1. Tips `S = {xx}`. `S_n = {x^n}`; `P_n ≅ A^e`; collapsed `C_n = A = ⟨1,x⟩`.
2. Leading map: `d_n` odd = `x⊗1 − 1⊗x`; even = `1⊗x + x⊗1`. No correction (monomial).
3. Collapse to `C_n`: odd `↦ 0`; even `↦` multiply-by-`2x` = `[[0,0],[2,0]]`, rank 1.
4. `HH_n = dim C_n − rank d_n − rank d_{n+1} = [2,1,1,1,…]` (char 0); `[2,2,…]` (char 2).

## The differential in general (the one subtle step)
The leading map `δ_n` is Bardzell's (depends only on the tips). The tails add a correction
strictly below `σ` in the reduction order; its coefficients are the solution of the linear
system `d_{n-1}∘d_n = 0` (CS Theorems 4.1/4.2 — the same trick CS use in §6). Two gates
certify every run: `assert_dd_zero` and `assert_order_condition`.

## Worked non-monomial example — the commutative square (HH^• = (1,0,0))
`A = kQ/(ab − cd)`, `Q`: `1→2→4`, `1→3→4` (arrows `a,b,c,d`), tip `cd`, `dim A = 9`.
`S_0 = {e_1,e_2,e_3,e_4}`, `S_1 = {a,b,c,d}`, `S_2 = {cd}`, `S_n = ∅ (n≥3)`.
Cohomology terms: `C^0 = ⟨ê_1..ê_4⟩` (dim 4), `C^1 = ⟨â,b̂,ĉ,d̂⟩` (dim 4), `C^2 = ⟨âb⟩` (dim 1).

`δ^0` sends `(λ_1,λ_2,λ_3,λ_4) ↦ ((λ_2−λ_1)a, (λ_4−λ_2)b, (λ_3−λ_1)c, (λ_4−λ_3)d)`:

```
        e1  e2  e3  e4
   a  [ -1   1   0   0 ]
   b  [  0  -1   0   1 ]        rank 3,  ker = <(1,1,1,1)> = the centre  ⟹  HH^0 = 1
   c  [ -1   0   1   0 ]
   d  [  0   0  -1   1 ]
```

`δ^1` sends a 1-cochain `g(a)=αa, g(b)=βb, g(c)=γc, g(d)=δd` to `(δ+γ−β−α)·ab`:

```
        a   b   c   d
  cd  [ -1  -1   1   1 ]        rank 1  ⟹  HH^2 = 1 − 1 = 0
```

`HH^0 = 4 − rank δ^0 = 1`,  `HH^1 = (4 − rank δ^1) − rank δ^0 = 3 − 3 = 0`,  `HH^2 = 1 − rank δ^1 = 0`.
Three cross-checks agree: Euler `4 − 4 + 1 = 1`; Gerstenhaber–Schack (order complex of the diamond
poset = two triangles glued on `{1,4}` = contractible ⟹ `HH^{>0}=0`); Künneth (`A = kA_2 ⊗ kA_2`,
`kA_2` a hereditary tree ⟹ `HH^0=k`, `HH^{≥1}=0`). Homology instead gives `HH_• = (4,0,0)` (the
quiver is acyclic, so `C_1 = C_2 = 0` and `HH_0 = A/[A,A] = k^{|Q_0|} = k^4`); `A` is not symmetric,
so `HH_0 = 4 ≠ 1 = HH^0`.

## What is certified, and what is not
- Deep dims: certified for `k[x]/(x^a)` and the quantum CI (byte-oracle) to any depth.
- General `kQ/I`: computed, certified per instance by `d²=0` + order gate + bar window,
  **restricted to quadratic tips or monomial presentations**; a non-quadratic non-monomial
  presentation raises `NotImplementedError` (the `right_decomposition` stretch item lifts this).
- Operations (cup/cap/bracket): transported to bar; certified only in the bar window.
```

- [ ] **Step 4: Full green + float gate green.** `cd /Users/marco/Desktop/HomologicalNetworks/quiverlab && NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q && NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/test_no_floats.py -q`

- [ ] **Step 5: Commit**
```bash
git add tests/resolutions_cs/test_acceptance.py docs/internals/08-chouhy-solotar.md
git commit -m "feat(cs): acceptance — deep quantum-CI HH_0..12 (ℚ) impossible for bar; internals ch.08

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Frozen interface for Plan 05 (Modules + invariants)

Plan 05 and later consume **exactly** the following from `quiverlab.resolutions_cs`, stable after Plan 04:

1. **`ChouhySolotarResolution(A, rs, max_degree, max_cells=4_000_000)`** — Domain-generic. `.d_terms(n, chain) -> list[(coeff∈domain, a_word, target_word, c_word)]` (term = `coeff·(a ⊗ target ⊗ c)` in `P_{n-1}`); `.delta_terms(n, chain)`; `.matrix(n, side)` (`side="hom"` boundary `d_n: C_n→C_{n-1}`, shape `dimC_{n-1}×dimC_n`; `side="coh"` coboundary `δ^n: C^n→C^{n+1}`, shape `dimC^{n+1}×dimC^n`); `.ss` (`SSequence`); `.dim_C(n, side)`; `.assert_dd_zero(upto, side)`; `.assert_order_condition(upto)`.
2. **`CSResolution(A)`** (`engine_facade`) — **`engine.resolutions.Resolution`-protocol compliant** over int64/GF(p): `term_basis`, `differential_matrix` (shape `(dim_{n-1}, dim_n)`, int64, un-reduced), `cochain_basis`, `coboundary_matrix` (shape `(dim_{n+1}, dim_n)`). Drop-in wherever `BardzellResolution`/`BarResolution` are accepted.
3. **`cs_cohomology_dims(A, top, max_cells=…)` / `cs_homology_dims(...)`** → `HHTable(engine="Chouhy-Solotar")`; **`cs_hh_basis(A, n, side)`** → representative (co)cycles; **`reduction_system_of(A)`** → the `ReductionSystem` for an `A` built by `Quiver.algebra`.
4. **Comparison-map API** (`comparison.Comparison(A)`): `.Phi(n)`, `.Psi(n)`, `transport_cocycle_cs_to_bar`, `transport_class_bar_to_cs`, `cup_of_cs_classes`/`cap_of_cs_classes`/`bracket_of_cs_classes` (valid inside the bar window; `NotImplementedError` past it).
5. **`ambiguities.SSequence`**: `.S(n) -> list[Chain]`, `Chain(word, blocks, o, t, degree)`. gl.dim upper bound `= max{n : S_n ≠ ∅}` when finite.
6. **Guards & errors:** `DepthLimitError` (certified range + resume), `AdmissibilityError` (non-confluent/infinite RS), `NotImplementedError` at the exact `(n, σ)` when the CS correction is inconsistent. The two invariants later plans rely on are **`Resolution`-protocol compliance** and the **comparison-map API**; the correction internals (`pelt.py`, `_d_general`, `_lower_generators`, `_solve`) are private and may change so long as (1)–(6), `assert_dd_zero`, and `assert_order_condition` hold.

---

## Self-review

**Spec §6/§8 coverage map.** §6 step 1 → Task 3 (S-sequence re-based on ported Bardzell `associated_paths`). §6 step 2 → Tasks 4 (corners), 5 (`δ_n`, `d_1`, `d_2`, swap-catcher), 6 (`_d_general = δ + Σγ·gen`, `γ` from CS Thm 4.2's `d²=0`, + the two gates). §6 step 3 → Task 12 (classes transported; operations stay on bar; window boundary). Battery (a)→8, (b)→9, (c)→10 (real BGMS `[3,2,…]`), (d)→11 (bank byte oracle). Risk register: admissibility gate → Task 7; `NotImplementedError` at exact `(n,σ)` → Task 6; loud guards + certified range → Tasks 3, 7, 14. §8 rings/batteries + trace dataclasses (Plan-07 rendering boundary) → Tasks 2, 13. Deepening/checkpoint UX → Tasks 3, 7.

**Four-pillar rework applied.** (1) **S-sequence** re-based on `MonomialPresentation.associated_paths`/`left_decomposition` (Task 3); the buggy hand-rolled `_extend`/`_proper_overlaps` is gone; `[1,1,1,1,1,1,1]` for `k[x]/x²` now holds because `associated_paths` is the validated CS §3 machinery. (2) **Differential** transcribed verbatim from CS §4 (`f_n` even/odd) with **tensor factors fixed** (odd/2-term = `u_0⊗rest⊗1 − 1⊗rest⊗u_last`; even/big-sum over sub-ambiguities), correction pinned by the CS `d²=0` **linear solve** (their §6 method), `δ_2 = φ_0(s) − φ_0(β(s))` with `β = rs.normal_form` (Pillar-4 unreduced-tail fix); a **non-commutative swap-catching test** added (Task 5 `test_d2_quantum_ci_yx_matches_paper_no_swap`, Task 6 `test_qci_d3_correction_matches_paper`, Task 11 byte oracle). All named helpers have complete bodies (`terms_to_pelt`, `apply_lower`, `_resolve_chain`, `_vecs`, `_accum`, `_lower_generators`, `_idx_word`, `_solve`, `_d_general`, `matrix`, `dim_C`, `_basis`, `_matmul`, `assert_dd_zero`, `assert_order_condition`, `_fox`, `delta_terms`, `_d1_terms`, `_d2_terms`); the fictional `tail_correction/_apply_lower/_lift_defect/pi/_is_zero_termlist` names are gone. (3) **Certification** adds `assert_order_condition` (CS Thm 4.1 cond 2) alongside `assert_dd_zero`, both wired into `cs_*_dims` so every reported dim is gated; Goal/acceptance retitled to what is certified (two bank families deep; general per-instance; operations window-bounded); the BGMS `[3,2,…,2]` and char-sweep vectors are the **bank-run** values (recorded from executing the read-only oracle, shown in the fixtures). (4) **Plan-03 consumption** rewritten to the committed interface: `build_reduction_system(quiver, relations, field)`; `ReductionSystem.{quiver,domain,order,rules,irreducibles,degree_bound,is_confluent,leading_words,reduce,normal_form,ambiguities}`; `ReductionRule(lead,tail,source,target)` with `tail` a tuple of `(coeff,word)` (Task-5 `_d2_terms` iterates `rs.normal_form(s).items()`, the **normal-form dict**, not the raw tuple tail); `Ambiguity(kind,word,left,right,a,c)`; words are arrow-name tuples; no `A.reduction_system()` (Task-7 `reduction_system_of(A)` reconstructs via the public builder + a domain shim); Task-1 gate asserts THIS contract and drops the `TAILS_INTERREDUCED` flag (encoded as the normal-form call + comment); Task 13 keeps shipped `engine="auto"` unchanged and adds an opt-in `auto_cs` flag with a compatibility test. **M2 fix:** the monomial-degeneracy formula is restated in CS degree `N = n − 1` (`_tp_word_degree(n-1)`).

**Placeholder scan (HONEST).** Every `src/` helper named in Tasks 2–7 and 12 has a complete body written in the plan (enumerated in the pillar-2 paragraph above). The only deliberate red state is Task-5's `_d_general` `NotImplementedError`, resolved in Task 6. The `...`-marked items are **test-side scaffolding the implementer completes in-task** (`_build`, `_build_gentle`, `assert_matrices_equal_mod_p`, `_res`/comparison helpers in tests), each explicitly flagged — never shipped `src/` placeholders. Two residual risks are stated, not hidden: (i) `fields.linalg.solve` must return any solution of a possibly-underdetermined consistent system and `None`/inconsistent-signal otherwise — if Plan-01's `solve` only handles square/consistent systems, Task 6 `_solve` adds the exact augmented-RREF particular-solution helper over the domain (this is called out inline in `_solve` and in Task-6 Step 4); (ii) the `_DomainField` shim assumes `build_reduction_system` calls only `field.parse_entry`/`make_domain` (grep Plan-03 `system.py`; add exactly the methods it calls if more).

**Signature consistency.** Term format `(coeff, a_word, target_word, c_word)` identical across `delta_terms`, `d_terms`, `pelt`, `matrix`, Task 12. `matrix(n, side)` shapes stated once (P2) and reused (Tasks 6, 7, facade, frozen interface). vs **frozen Plan-03** (verified against the committed doc): uses only `build_reduction_system`, `ReductionSystem.{quiver,domain,order,rules,irreducibles,leading_words,reduce,normal_form,ambiguities,is_confluent,degree_bound}`, `ReductionRule.{lead,tail,source,target}`, `Ambiguity.{kind,word,left,right,a,c}`, `PathOrder.key`, `Quiver.{vertices,arrows,source,target}`, `core.Algebra.{domain,dim,multiply,_basis_vec,unit_adapted,quiver,relations}` — all asserted by the Task-1 gate. vs **actual engine sources** (verified): `engine.resolutions.Resolution` four-method contract + shapes (`resolutions.py`); `to_engine`/`PrimeField` (`adapter.py`); `MonomialPresentation.{associated_paths,left_decomposition,path_src,path_tgt,maxrel}` + `BardzellResolution` (`resolutions_bardzell.py` 111–335); `hochschild_homology_dims(..., resolution=)` (`hh_engine.py:178`), `hochschild_cohomology_dims(..., resolution=)` (`scan3.py:75`); `HHTable(dims, prefix, name, engine=)` (`table.py`); homology-dim identity (`bar.py`); collapse `b·w·a` (bank `resolutions_cs.py`); `DepthLimitError(msg, hint=)`/`AdmissibilityError` (`errors.py`).

**Convention integrity (composition direction, every formula).** Left-to-right in: `φ_0(c_1⋯c_ℓ) = Σ_k (c_1..c_{k-1})⊗c_k⊗(c_{k+1}..c_ℓ)` (verified against CS §6 `d_1(yx) = y⊗x⊗1 + 1⊗y⊗x − ξx⊗y⊗1 − ξ1⊗x⊗y`); the `n`-odd 2-term map `(+1,u_0,[u_1..u_{n-1}],()) , (−1,(),[u_0..u_{n-2}],u_{n-1})` (verified against CS §6 `d_2(x³) = x⊗x²⊗1 − 1⊗x²⊗x`); the `n`-even big sum over sub-ambiguities `σ=a·p·c`; homology collapse `b·w·a`, cohomology `a·w·b` (P2); `≺` = the Plan-03 admissible order on parallel paths (order gate). The "Convention translation" note (P3) lists the four translated items and the `v_{n-i}=u_i` factorization identity. Fixtures computed entirely left-to-right (`δ^0(α:i→j)=(λ_j−λ_i)α`, `δ^1(cd)=(δ+γ−β−α)ab`).

**Fixes applied inline during this rework:** (1) S-sequence re-based on ported Bardzell (Pillar 1); (2) differential transcribed verbatim with fixed tensor factors + swap-catcher + `β` normal-form on tails (Pillar 2); (3) order-condition gate + honest certification retitling + real bank-run vectors (Pillar 3); (4) every Plan-03 call site rewritten to the committed interface, Task-1 gate rewritten, `auto` untouched with opt-in `auto_cs` (Pillar 4); (5) M2 degree-variable fix; (6) added `docs/internals/08-chouhy-solotar.md`; every named `src/` helper has a complete body. **Four post-COMMIT final edits:** (E1) RESTRICT scope — certified for quadratic tips (CS Prop, TeX 772) or monomial presentations; `_require_in_scope` raises `NotImplementedError` on non-quadratic non-monomial at the exact degree-≥3 boundary, with `test_cubic_tip_nonmonomial_raises_notimplemented`; `right_decomposition` named as the lifting stretch item (Goal + Task 6). (E2) solve non-uniqueness — `test_qci_d3_correction_matches_paper` and the Task-11 entrywise byte pins marked `xfail(strict=False)`, binding criteria are `d²=0` + order gate + HH-dim equality; canonicalization (normal form mod solve-nullspace) is the stretch item that flips them strict. (E3) `_DomainField` shim — Task-1 Step 5 greps `groebner/system.py` for every `field.<method>` call (STOP-and-reconcile; `parse_entry` must pass `Fraction` through for `dom.coerce`). (E4) docs/internals/08 Fixture-B matrices filled in (δ^0 incidence rank 3, δ^1 = [−1,−1,1,1] rank 1, three cross-checks) + an execution-time note to reconfirm the `p=5` acceptance vector against the bank oracle.
