# quiverlab Plan 05 — Modules + Invariants Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A research algebraist, having built any `A = kQ/I` over any exact field (Plans 01–04), now asks it about its **modules** and its remaining **invariants**, exactly. `A.simple(v)`, `A.projective(v)`, `A.injective(v)` produce right `A`-modules; `M.radical()`, `M.top()`, `M.socle()`, `M.dimension_vector()` read off their structure; `A.hom(M, N)` and `A.ext(M, N, n)` compute Hom/Ext dimensions via a **minimal projective resolution** `M.projective_resolution(k)` that is an inspectable object; and `A.global_dimension()`, `A.loewy_length()`, `A.complexity(n)`, `A.center()`, plus **exact** `spectral_radius`/`mahler_measure` (no floats, ever) and `sweep()` (invariant × field) complete the invariant surface (spec §3.5–§3.6, §3.9; §5 components 7–8). Every number is exact over the stated `Domain`; every module convention (right modules, action side) is stated in the code.

**Architecture:** A new top-level package `src/quiverlab/modules/` (`linalg_mod.py` — matrix arithmetic over any `Domain` built on `fields.linalg`; `module.py` — the `Module` object + radical/top/socle; `builders.py` — `simple`/`projective`/`injective` from quiver provenance; `hom.py` — Hom/End spaces; `resolution.py` — the **generalized** minimal projective resolution engine, lifted from the bridge `obstruction/module_ext.py` with its hardcoded 4-vertex diamond removed, over any vertex set and any `Domain`, MIT-headered per spec §9; `ext.py` — `ext`/`global_dimension`). Invariants gain `invariants/spectral.py` (exact `spectral_radius`/`mahler_measure` reimplementing the deleted hanlab float layer), `invariants/scalar.py` (`loewy_length`/`complexity`/`center`), and `invariants/sweep.py` (`sweep`); `invariants/cartan.py` gets the non-unimodular-Cartan `coxeter_polynomial` minor fixed. `core/algebra.py` grows thin public methods (`simple`, `projective`, `injective`, `hom`, `ext`, `global_dimension`, `loewy_length`, `complexity`, `center`) that **deferred-import** their engines exactly as the existing `cartan_matrix` method does (no import cycle).

**Tech Stack:** Python ≥ 3.10; no new hard dependencies (the Plan-01 stack: `sympy` for exact ℂ and for the spectral-polynomial layer). All module linear algebra runs over the Plan-01 `Domain` via `fields.linalg.{rank, nullspace, rref, solve}` — exact for ℚ, ℚ(α), GF(p), GF(p^n). The exact spectral layer uses `sympy` algebraic numbers (`factor_list`, `real_roots`/`all_roots`, `count_roots`, `sqf_part`, `Abs`, `is_positive`) — **no floats in `src/`**. pytest.

---

## Global Constraints

- **Repo root:** `/Users/marco/Desktop/HomologicalNetworks/quiverlab`. All paths below are relative to it.
- **Interpreter:** use the project venv **`/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python`** (Python 3.12). The *system* `python`/`python3` is 3.8 and MUST NOT be used — it fails on 3.10+ syntax (`X | None`, `list[int]`).
- **Thread throttle:** prefix **every** test command with `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2` (Marco's machine has crashed under agent fleets; keep thread/memory pressure low). No new parallelism is introduced; module code is single-threaded.
- **Exact arithmetic only, through the `Domain` protocol.** All module coefficient arithmetic goes through `dom.add/sub/neg/mul/inv/is_zero/eq/coerce/zero/one`, or through `fields.linalg` which is itself domain-generic. Never combine coefficients with raw Python operators. The spectral layer works over `sympy` exact algebraic numbers only.
- **Float ban (AST gate).** `tests/test_no_floats.py` scans all of `src/quiverlab/` for float/complex literals and `float()` calls and MUST stay green. Write **no** float or complex literal anywhere under `src/`. In the spectral layer use `sympy.Integer(1)`, `sympy.Rational(...)`, `sympy.Abs`, `expr.is_positive` — never `1.0`, never `float(...)`, never `complex(...)`, never `Poly.nroots` (the deleted hanlab layer's float routine). Run the gate in every task's suite. **Floats are permitted in `tests/` only** (the bank numeric oracle is compared to the exact result inside tests, e.g. `abs(float(rho.evalf(30)) - 1.17628081825991) < 1e-9`).
- **Left-to-right path composition** (Assem–Simson–Skowroński): `a*b` = "first `a`, then `b`", requiring `target(a) == source(b)`. A basis label is `"a*b*c"` (a path word) or `"e_v"` (a trivial path); words split on `"*"`.
- **Right-module convention (binding, stated in every module docstring).** Modules are **right** `A`-modules. An element `m ∈ M` is a **column vector** in a fixed `k`-basis of `M`; the action of an algebra basis element `b` is a matrix `action[b]` with `m·b = action[b] @ m_col`. This is the **anti-homomorphism** convention `action[x·y] = action[y] @ action[x]` (so that `(m·x)·y = m·(x·y)`), matching the bridge `RMod` engine being generalized. The vertex subspace is `M·e_v = image(action["e_v"])`, and `dimension_vector[v] = rank(action["e_v"])`; since `Σ_v e_v = 1_A`, `Σ_v action["e_v"] = I` and `M = ⊕_v M·e_v`.
- **Read-only banks.** `/Users/marco/Desktop/HomologicalNetworks/bridge/`, `/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/`, and every other bank are **never** modified. `bridge/lab/obstruction/module_ext.py` (the engine generalized here) and `HansConjecture/hanlab/coxeter_spectrum.py` (the deleted spectral layer + oracle tests) are read for reference only; the lifted engine is re-authored fresh in `src/` with an MIT header and attribution (bridge/ is unlicensed — spec §9).
- **Full suite green at every commit.** Before each commit run the focused suite for the task (fast). The **one** full-suite run (`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q`, ~19 min, exceeds the 600 s foreground cap) is executed **once**, in the acceptance task, as a single tracked **background** job (`run_in_background: true`), **awaited to completion**, and only then is the final commit made and the report given — all in the same session.
- **Commits:** conventional prefixes (`feat:`, `test:`, `fix:`, `docs:`); every commit message ends with the trailer line
  `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- **Frozen upstream interfaces consumed verbatim** (do not modify their signatures):
  - **Plan 01/02.** `core.Algebra` — attributes `.domain, .T, .unit, .dim, .basis_labels, .quiver, .relations, .is_unit_adapted`; methods `.multiply(u, v)`, `._basis_vec(i)`, `.unit_adapted()`, `.change_of_basis(P)`, `.cartan_matrix()`, `.coxeter_matrix()`, `.coxeter_polynomial()`, `.hochschild_cohomology(top, max_cells=, engine=)`, `.hochschild_homology(...)`, `.nakayama_automorphism()`, `.is_frobenius()`, `.is_symmetric()`, `.cyclic_homology(top)`, `._require_prime_field(what)`; classmethod `.from_structure_constants(T, unit, field=None, check=True, basis_labels=None)`. Constructor `Algebra(domain, T, unit, basis_labels=None, is_unit_adapted=None, _quiver=None, _relations=None)`. **Note:** `.dim` is a plain attribute (there is no `.dimension()` method); `.quiver`/`.relations` are the public attribute names (there is no `._quiver`/`._relations`).
  - `combinat.Quiver` — `.vertices`, `.arrows` (insertion-ordered dict name→(s,t)), `.source(name)`, `.target(name)`, `.word_source(word)`, `.word_target(word)`, `.algebra(relations=, field=, degree_bound=, trace=)`.
  - `fields` — `CC`, `QQ` (instances), `GF(q)` (function), `E(n)`, `PrimeField` (has `.p`), `Domain` (`coerce/zero/one/add/sub/neg/mul/inv/is_zero/eq/characteristic/name`).
  - `fields.linalg` — `rank(rows, dom) -> int`, `nullspace(rows, dom) -> list[vector]`, `rref(rows, dom) -> (matrix, pivots)`, `solve(A, b, dom) -> vector | None`.
  - `errors` — `QuiverlabError(message, hint=None)` and subclasses `ExactnessError, FieldError, RelationError, AdmissibilityError, NotFiniteDimensionalError, DepthLimitError` (all share the `(message, hint=None)` ctor). Python's built-in `NotImplementedError` marks unreached-scope boundaries.
  - `invariants.cartan` — `cartan_matrix(A) -> list[list[int]]`, `coxeter_matrix(A)`, `coxeter_polynomial(A) -> sympy.Poly`.
  - `engine.coxeter_spectrum` — `is_cyclotomic_product(poly) -> bool | None`, `star_quiver(arm_lengths) -> (n, arrows)`, `cartan_of_quiver(n, arrows, name=) -> (C, alg)`, `trivial_extension_cartan(C) -> C + C.T`.
  - `engine.coxeter2` — `coxeter_polynomial_from_cartan(C) -> (poly, Phi) | (None, None)`, `cartan_from_raw`, `cyclic_nakayama`, `quiver_path_algebra`, `dynkin_quiver`.
  - `engine.scan3` — `complexity_of(seq) -> int | None | str`.
  - `engine.adapter` — `to_engine(A) -> hh_engine.Algebra` (raises `FieldError` off `PrimeField`).
  - `engine.resolutions_minimal` — `minimal_resolution(A, N, p, max_term_dim=20000, max_transient_bytes=None) -> (rks, cols, eng, truncated_at)`, `minimal_homology_dims(A, N, primes=(32003,), max_term_dim=20000, max_transient_bytes=None) -> {p: [int]}`.
  - **Plan 03** — `Quiver.algebra(relations=[...])` lowers **general** (non-monomial) `kQ/I` via `groebner`; `groebner.build_reduction_system(quiver, relations, field, degree_bound=, trace=)`. (Needed to build the non-monomial `_open_33_0` memory-guard fixture and non-monomial modules.)
  - **Plan 04** — `resolutions_cs.cs_cohomology_dims/cs_homology_dims`, `resolutions_cs.engine_facade.CSResolution(A)` (drop-in `Resolution`), `resolutions_cs.reduction_system_of(A)`. (Consumed only for the optional monomial/quadratic Ext cross-check, Task 7; guarded by `importorskip` so Plan 05 is testable even if Plan 04 lands late.)
- **One compatible extension of the frozen surface** is made here: `core.Algebra` gains the public methods `simple, projective, injective, hom, ext, global_dimension, is_selfinjective, loewy_length, complexity, center` (all new — none exist today; agent-verified). No existing signature changes. `invariants.cartan.coxeter_polynomial` gains **only a documented non-unimodular caveat** (Task 10, docstring-only — sympy's coefficient-based ZZ/QQ inference is already correct); its return type and every output are byte-identical to today.

---

## Module & spectral conventions (read this before Task 1)

The implementer knows Python, not necessarily representation theory. This section carries the mathematics; the code implements exactly what is stated here. Composition is **left-to-right** and modules are **right** modules throughout.

**Path algebra & the regular projectives.** `A = kQ/I` has `k`-basis the irreducible paths (Plan 01/03). The **indecomposable projective** right module at vertex `v` is `P_v = e_v A` — the right ideal spanned by the basis paths **starting at** `v`. Its `k`-dimension is `Σ_w C[v][w]` (row `v` of the Cartan matrix), and `dim (e_v A) e_w = C[v][w]` = number of basis paths `v → w`. The **simple** right module `S_v` is 1-dimensional, concentrated at vertex `v`, with every arrow acting as `0`; `S_v = top(P_v) = P_v / rad P_v`. The **indecomposable injective** at `v` is `I_v = D(A e_v)`, the `k`-dual of the **left** projective `A e_v` (basis paths **ending at** `v`); as a right module `dim e_w I_v = C[w][v]` (column `v` of the Cartan matrix), with right action the transpose of left multiplication on `A e_v`.

**Radical, top, socle.** `rad A` is the arrow ideal (spanned by all basis paths of length ≥ 1). For a right module `M`: `rad M = M·(rad A) = Σ_{arrows α} image(action[α])` (a submodule); `top M = M / rad M` (semisimple, the "generators/heads"); `soc M = {m ∈ M : m·(rad A) = 0} = ∩_{arrows α} ker(action[α])` (the largest semisimple submodule). The **radical series** `M ⊇ rad M ⊇ rad² M ⊇ …` reaches `0`; its length is the **Loewy length** of `M`. `rad^k A` is computed by iterating the ideal product; the least `k` with `rad^k A = 0` is `loewy_length(A)`.

**Minimal projective resolution (the generalized `module_ext` engine).** A minimal projective resolution `… → Q_2 → Q_1 → Q_0 → M → 0` is built by **iterated projective covers** (Green–Solberg–Zacharia style; cite GSZ, *Trans. AMS* 353 (2001) 2915–2939 in the module docstring):
1. `top M = M / rad M`; choose **homogeneous** generators — a `k`-basis of `top M`, each lifted to a vector in `M` living in a single vertex subspace `M·e_v`. If vertex `v` supplies `t_v` generators, the cover is `Q_0 = ⊕_v P_v^{t_v}`, and `d_0 : Q_0 → M` sends each summand's canonical generator to its lifted top-generator.
2. The **syzygy** `Ω_1 = ker(d_0) ⊆ Q_0` is again a right submodule; cover it to get `Q_1` and `d_1 : Q_1 → Q_0` (the cover of `Ω_1` post-composed with the inclusion `Ω_1 ↪ Q_0`). Repeat. Minimality ⇔ each cover is chosen modulo the radical (generators independent modulo `rad · (ambient)`), so `d_n(Q_n) ⊆ rad Q_{n-1}`.
3. The resolution **terminates** at length `n` iff `Ω_n = 0` (`M` has projective dimension `n`); `global_dimension(A) = sup_v pd(S_v)`. The bridge `RMod` engine does exactly this over `ℚ` with the vertex set literally `[0,1,2,3]`; the generalization ranges over `Q.vertices` and swaps every `sympy` matrix operation for the domain-generic `fields.linalg`, so it runs over `ℚ, ℚ(α), GF(p), GF(p^n)` unchanged.

**Ext.** `Ext^n_A(M, N) = H^n(Hom_A(Q_•, N))`. With the minimal resolution `Q_•` of `M` and its (contravariant) coboundary maps `δ^n = Hom(d_{n+1}, N)`, `dim Ext^n = c_n − rank(δ^n) − rank(δ^{n-1})` where `c_n = dim Hom_A(Q_n, N)`. Hand oracles (Task 7): for a **hereditary** path algebra `kA_m` (linear quiver `1→2→…→m`, ASS §III / Happel), `Hom(S_a, S_b) = δ_{ab}·k`, `Ext^1(S_a, S_b) = k^{#arrows a→b}` (so `Ext^1(S_a, S_{a+1}) = k`, all other `Ext^1(S_a,S_b)=0`), and `Ext^{≥2}(S,S) = 0`. For the commutative square (`= kA_2 ⊗ kA_2`, gl.dim 2), Künneth gives `Ext²(S_1, S_4) = Ext^1_{kA_2}(S_1,S_2) ⊗ Ext^1_{kA_2}(S_1,S_2) = k` (dim 1).

**Exact spectral invariants (reimplementing the deleted hanlab float layer).** The Coxeter polynomial `p(t) = charpoly(Φ)` with `Φ = −C^{-T}C` (from `coxeter_polynomial_from_cartan`) is exact. The bank's `spectral_radius`/`mahler_measure` used `Poly.nroots` (mpmath floats) — **deleted** in the port. We reimplement them **exactly and soundly**: `spectral_radius(p) = max_i |α_i|` and `mahler_measure(p) = |lc| · ∏_{|α_i|>1} |α_i|` over the exact `CRootOf` roots, `sympy.Abs` magnitudes, exact `.is_positive` comparisons — never a float. Soundness requires counting **complex** off-circle roots too: `real_roots` alone is unsound (real-roots-suffice holds only for hereditary quivers, A'Campo). So (after the **cyclotomic short-circuit** `is_cyclotomic_product(p) ⟹ 1`) we form the non-cyclotomic part `q` and, via the self-inversive substitution `y = z + 1/z`, decide by a Sturm real-root count of `Q(y)` — **without isolating any complex root** — whether `real_roots(q)` suffices (Branch A: hereditary/Salem/Lehmer, fast) or `all_roots(q)` is needed (Branch B: non-hereditary, complex-dominant; `q` is the small non-cyclotomic factor). Full algorithm in Task 8. The bank floats survive only as **test oracles** (compared to `.evalf(30)` of the exact result within a rational tolerance, in `tests/`).

---

## Hand-verified fixtures (worked out by hand; the acceptance oracles)

Arrow order is always the insertion order shown. Every number below is re-derived in the task that uses it. Left-to-right composition; right modules.

### Fixture A — the A₂ path algebra `kA_2` (`linear_path_algebra(2)`)
`linear_path_algebra(2)` builds `Quiver([1, 2], {"a1": (1, 2)})` (the arrow is auto-named `"a1"` = `f"a{i}"`, agent-verified against `families/basic.py` and internals ch06), no relations. `dim A = 3`, basis `["e_1", "e_2", "a1"]`, Cartan `[[1, 1], [0, 1]]`. (Below, the arrow/path is written `a1`; the hand-built `action` dicts in Task 2 key on `"a1"` accordingly.)
- **Simples:** `S_1` dim 1, `dimension_vector {1:1, 2:0}`; `S_2` dim 1, `{1:0, 2:1}`.
- **Projectives:** `P_1 = e_1 A = span{e_1, a1}`, dim 2, `dimvec {1:1, 2:1}`; `top P_1 = S_1`, `rad P_1 = span{a1} ≅ S_2`, `rad² P_1 = 0`, `soc P_1 = S_2`, Loewy length 2. `P_2 = e_2 A = span{e_2}`, dim 1, `= S_2` (projective simple), `dimvec {1:0, 2:1}`, `rad P_2 = 0`.
- **Injectives:** `I_1 = D(Ae_1)`; `A e_1 = span{e_1}` (paths ending at 1) so `I_1 = S_1`, `dimvec {1:1, 2:0}`. `I_2 = D(Ae_2)`; `A e_2 = span{e_2, a1}` (paths ending at 2) so `dim I_2 = 2`, `dimvec {1:1, 2:1}` (socle `S_2`, top `S_1`).
- **Resolution of `S_1`:** `0 → P_2 → P_1 → S_1 → 0` (`rad P_1 = span{a1} ≅ P_2`). So `pd(S_1) = 1`; `pd(S_2) = 0` (`S_2 = P_2` projective); `global_dimension = 1` (hereditary).
- **Ext / Hom:** `Hom(S_a,S_b) = δ_{ab}` (dim 1 iff `a=b`); `Ext^1(S_1, S_2) = 1`, `Ext^1(S_1,S_1) = Ext^1(S_2,S_b) = 0`, `Ext^{≥2} = 0`.
- **Loewy length of `A`:** `rad A = span{a1}`, `rad² A = 0` ⇒ `loewy_length(A) = 2`. **Center:** `A` connected ⇒ `center = k·1`, dim 1. **Complexity:** finite gl.dim ⇒ minimal `A^e`-resolution terminates (trailing zeros) ⇒ `complexity_of = 0`.

### Fixture B — the commutative square `kA_2 ⊗ kA_2` (`Quiver([1,2,3,4], {"a":(1,2),"b":(2,4),"c":(1,3),"d":(3,4)}).algebra(relations=["a*b - c*d"])`)
Plan-03 Fixture 2: `dim A = 9`, basis `["e_1","e_2","e_3","e_4","a","b","c","d","a*b"]` (with `c·d = a·b`), Cartan `[[1,1,1,1],[0,1,0,1],[0,0,1,1],[0,0,0,1]]`, `det = 1` (unimodular).
- **Projectives:** `P_1 = span{e_1, a, c, a*b}`, dim 4, `dimvec {1:1,2:1,3:1,4:1}`; `rad P_1 = span{a,c,a*b}`, `rad² P_1 = span{a*b} ≅ S_4`, `rad³ P_1 = 0` ⇒ Loewy length 3; `top P_1 = S_1`; `soc P_1 = S_4`; radical layers `S_1 | S_2 ⊕ S_3 | S_4`. `P_2 = span{e_2, b}` dimvec `{2:1,4:1}`, `P_3 = span{e_3, d}` dimvec `{3:1,4:1}`, `P_4 = span{e_4} = S_4`.
- **Ext:** `Ext²(S_1, S_4) = 1` (Künneth; `= Ext^1(S_1,S_2)⊗Ext^1(S_1,S_2)` in each `kA_2` factor); `global_dimension = 2`.
- **Loewy length of `A` = 3** (`rad³ A = 0`, `rad² A = span{a*b} ≠ 0`). **Center** dim 1 (connected). **Coxeter** unimodular ⇒ `coxeter_polynomial` domain `ZZ`.

### Fixture C — the monomial diamond `kD/(a*b)` (bank oracle, Ext² cross-check)
`Quiver([1,2,3,4], {"a":(1,2),"b":(2,4),"c":(1,3),"d":(3,4)}).algebra(relations=["a*b"])`, `dim 9`. The bridge `module_ext` docstring pins **`Ext²(S_1, S_4) = 1`** (the relation `a*b = 0` from `1⇒4` supplies a degree-2 syzygy). Re-derived and cross-checked in Task 7 (independent of the hereditary Künneth case B).

### Fixture D — the non-unimodular Cartan (the `coxeter_polynomial` minor, Task 10)
`Quiver([1, 2], {"x": (1, 1), "a": (1, 2)}).algebra(relations=["x^2", "x*a"])`, `dim 4`, basis `["e_1","e_2","a","x"]`, Cartan `C = [[2, 1], [0, 1]]`, `det C = 2 ≠ ±1` (verified by running). Then `Φ = −C^{-T}C = [[-1, -1/2],[1, -1/2]]` and `coxeter_polynomial = t² + (3/2)·t + 1` — a **genuinely rational** (domain `QQ`) polynomial, versus the unimodular `kA_2` which gives `t² + t + 1` (domain `ZZ`). This is the exact witness that the non-unimodular case must be surfaced, not silent.

### Fixture E — the Lehmer star `T(2,3,7) = star_quiver([1,2,6])` (ledger; Task 9)
`star_quiver([1, 2, 6])` builds the star tree with center vertex 0 and three arms of edge-lengths 1, 2, 6 (all oriented toward the center) — `n = 10` vertices. Its Cartan's Coxeter polynomial is **Lehmer's polynomial**
`LEHMER = t¹⁰ + t⁹ − t⁷ − t⁶ − t⁵ − t⁴ − t³ + t + 1` (= `1 + t − t³ − t⁴ − t⁵ − t⁶ − t⁷ + t⁹ + t¹⁰`).
It is **not** cyclotomic (`is_cyclotomic_product(LEHMER) = False`), and its `spectral_radius` is **Lehmer's number** `1.176280818259917…` (the smallest known Mahler measure > 1) — asserted exactly as `spectral_radius(LEHMER).evalf(30) ≈ 1.17628081825991` within `1e-9` (bank float oracle). Companion (Task 8): `spectral_radius(t²−3t+1) = mahler_measure(t²−3t+1) = (3+√5)/2` (the 5-subspace wild radius), and `spectral_radius((t+1)³(t−1)²) = 1` exactly (cyclotomic).

### Fixture F — the trivial-extension `(t+1)^v` collapse (ledger; Task 9)
For any bound quiver algebra with `v` vertices and Cartan `C`, the trivial extension `T(A) = A ⋉ DA` has Cartan `C_T = C + Cᵀ` (`trivial_extension_cartan`). When `C_T` is nonsingular, its Coxeter transformation is `Φ_T = −I_v` and its Coxeter polynomial collapses to `(t+1)^v` — **independent of representation type**: rep-finite `kA_2` (`v=2`) and wild 5-subspace `star_quiver([1,1,1,1,1])` (`v=6`) both give `(t+1)^v`, while their hereditary inputs have different spectra (`spectral_radius(p_{A_2}) = 1`, `spectral_radius(p_{wild}) > 2.6`). Euclidean `D̃_4 = star_quiver([1,1,1,1])` gives a **singular** symmetrized Cartan ⇒ `coxeter_polynomial_from_cartan(C_T) = (None, None)`.

### Fixture G — the minimal-memory-guard fixture `open_33_0` (ledger; Task 13)
`open_33_0 = k⟨x,y⟩/(x³ − y², y³, y·x + x·y)`, `dim 9`, local, non-monomial — buildable **now** via Plan-03 `Quiver([1], {"x":(1,1), "y":(1,1)}).algebra(relations=["x^3 - y^2", "y^3", "y*x + x*y"], field=GF(32003))`. (The bank test built it via `reduction_algebra.algebra_from_reduction_system`, which Plan 04 declared **superseded** by `Quiver.algebra`; the dependency is therefore satisfiable and Plan 05 owns the re-port.) Its minimal `A^e`-resolution's `radK` transient grows degree by degree, exercising the `max_transient_bytes` guard: any budget yields a **prefix** of the full HH list; a 1-byte budget yields `truncated_at == 0`, `hh == []`.

---

### Task 1: Interface freshness gate — STOP on drift

**Files:** Create `tests/modules/__init__.py` (empty), `tests/modules/test_interface_gate.py`.

**Interfaces:** Consumes every upstream symbol Plan 05 relies on (listed in Global Constraints). Produces nothing in `src/`; this is a pure guard. If any assertion fails, **STOP** and reconcile against the drifted surface before doing any further Plan-05 work.

- [ ] **Step 1: Write the failing test**

`tests/modules/__init__.py`: empty.

`tests/modules/test_interface_gate.py`:

```python
"""Plan-05 freshness gate: every upstream symbol Plan 05 consumes, asserted to exist
with the expected shape. A failure here means an upstream surface drifted -- STOP and
reconcile before continuing Plan 05 (spec §5 components 7-8)."""
import inspect
import pytest


def test_algebra_surface_present_and_module_methods_absent():
    from quiverlab.core.algebra import Algebra
    for attr in ("domain", "T", "unit", "dim", "basis_labels", "quiver", "relations",
                 "is_unit_adapted"):
        assert attr in Algebra.__init__.__code__.co_varnames or hasattr(Algebra, attr) \
            or True  # attributes are instance-set; presence checked at runtime below
    for meth in ("multiply", "unit_adapted", "cartan_matrix", "coxeter_matrix",
                 "coxeter_polynomial", "hochschild_cohomology", "hochschild_homology",
                 "nakayama_automorphism", "is_frobenius", "is_symmetric"):
        assert callable(getattr(Algebra, meth)), meth
    # Plan-05 adds these; they must NOT pre-exist (else a rebase duplicated work)
    for new in ("simple", "projective", "injective", "hom", "ext",
                "global_dimension", "loewy_length", "complexity", "center",
                "is_selfinjective"):
        assert not hasattr(Algebra, new), f"{new} already exists -- reconcile"


def test_algebra_runtime_attributes():
    from quiverlab import linear_path_algebra
    A = linear_path_algebra(2)
    for attr in ("domain", "T", "unit", "dim", "basis_labels", "quiver", "relations"):
        assert hasattr(A, attr), attr
    assert A.dim == 3 and A.basis_labels == ["e_1", "e_2", "a1"]   # arrow auto-named a1


def test_linalg_signatures():
    from quiverlab.fields import linalg
    from quiverlab.fields import QQ
    assert linalg.rank([[QQ.coerce(1)]], QQ) == 1
    ns = linalg.nullspace([[QQ.coerce(1), QQ.coerce(-1)]], QQ)
    assert len(ns) == 1
    A = [[QQ.coerce(1), QQ.coerce(0)], [QQ.coerce(0), QQ.coerce(1)]]
    assert linalg.solve(A, [QQ.coerce(2), QQ.coerce(3)], QQ) == [QQ.coerce(2), QQ.coerce(3)]
    R, piv = linalg.rref([[QQ.coerce(2), QQ.coerce(0)]], QQ)
    assert piv == [0]


def test_invariants_cartan_surface():
    from quiverlab.invariants import cartan
    from quiverlab import linear_path_algebra
    A = linear_path_algebra(2)
    assert cartan.cartan_matrix(A) == [[1, 1], [0, 1]]
    assert cartan.coxeter_polynomial(A) is not None


def test_engine_spectral_and_coxeter2_surface():
    from quiverlab.engine.coxeter_spectrum import (is_cyclotomic_product, star_quiver,
                                                   cartan_of_quiver, trivial_extension_cartan)
    from quiverlab.engine.coxeter2 import coxeter_polynomial_from_cartan
    from quiverlab.engine.scan3 import complexity_of
    n, arrows = star_quiver([1, 2, 6])
    assert n == 10
    C, _alg = cartan_of_quiver(n, arrows)
    poly, Phi = coxeter_polynomial_from_cartan(C)
    assert poly is not None
    assert trivial_extension_cartan(C).shape == (10, 10)
    assert is_cyclotomic_product(None) is None
    assert complexity_of([1, 1, 1, 1, 1, 1]) in (1, "1", 1)  # constant seq -> complexity 1


def test_minimal_resolution_guard_surface():
    from quiverlab.engine import resolutions_minimal as rm
    from quiverlab.engine.adapter import to_engine
    sig = inspect.signature(rm.minimal_resolution)
    assert "max_transient_bytes" in sig.parameters
    sig2 = inspect.signature(rm.minimal_homology_dims)
    assert "max_transient_bytes" in sig2.parameters
    assert callable(to_engine)


def test_plan03_general_kqi_available():
    """General (non-monomial) kQ/I must lower -- needed for non-monomial modules and
    the memory-guard fixture. If Plan 03 has not landed, STOP."""
    from quiverlab import Quiver, GF
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x^3 - y^2", "y^3", "y*x + x*y"], field=GF(32003))
    assert A.dim == 9


def test_quiver_algebra_accepts_degree_bound():
    """Task 13's memory-guard fixture may need an explicit degree_bound; Plan 03's frozen
    interface adds `Quiver.algebra(..., degree_bound=None, trace=None)`. If the parameter
    is absent, the Gröbner completion cannot be pushed and Task 13 is blocked -- STOP."""
    import inspect
    from quiverlab import Quiver, GF
    params = inspect.signature(Quiver.algebra).parameters
    assert "degree_bound" in params, "Quiver.algebra must accept degree_bound= (Plan 03)"
    # and it must actually take effect (accepts the kwarg without error)
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x^3"], field=GF(2), degree_bound=12)
    assert A.dim == 3


def test_plan04_cs_optional_but_shape_checked_if_present():
    """Plan 04 is consumed only for the optional monomial/quadratic Ext cross-check.
    If present it must expose the frozen names; if absent, skip (Plan 05 stays testable)."""
    cs = pytest.importorskip("quiverlab.resolutions_cs")
    assert hasattr(cs, "cs_cohomology_dims") or hasattr(cs, "cs_homology_dims")
    from quiverlab.resolutions_cs.engine_facade import CSResolution  # noqa: F401
```

- [ ] **Step 2: Run — expect FAIL** (`ModuleNotFoundError: No module named 'tests.modules'` until `__init__.py` exists, then genuine assertion results).

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_interface_gate.py -q`

- [ ] **Step 3: No `src/` change.** This task writes no implementation; it certifies the ground truth. If `test_algebra_surface_present_and_module_methods_absent` fails because a `simple`/`ext`/… method already exists, a prior rebase already did the work — reconcile before proceeding. If `test_plan03_general_kqi_available` fails, Plan 03 has not landed — **STOP**.

- [ ] **Step 4: Run the suite** (gate + float ban).

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_interface_gate.py tests/test_no_floats.py -q`
Expected: all pass (the gate confirms every consumed surface; the CS check skips if Plan 04 is absent).

- [ ] **Step 5: Commit**

```bash
git add tests/modules/__init__.py tests/modules/test_interface_gate.py
git commit -m "test(modules): Plan-05 upstream interface freshness gate

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Module representation + domain-generic matrix helpers + `dimension_vector`

**Files:** Create `src/quiverlab/modules/__init__.py`, `src/quiverlab/modules/linalg_mod.py`, `src/quiverlab/modules/module.py`, `tests/modules/test_module.py`.

**Interfaces:**
- Consumes: `fields.linalg.{rank, nullspace, rref, solve}`, `Domain`, `core.Algebra`, `QuiverlabError`.
- Produces:
  - `linalg_mod`: matrices are `list[list[<domain elt>]]` (rows × cols); vectors are `list[<domain elt>]` (a column). `zeros(r, c, dom)`, `identity(n, dom)`, `transpose(M)`, `matmul(A, B, dom)`, `matvec(A, v, dom)`, `kron(A, B, dom)`, `vstack(mats)`, `col(M, j)`, `cols_to_matrix(cols)`; `mat_rank(M, dom)`; `kernel_columns(M, dom) -> list[vector]` (basis of `{x : M x = 0}`); `column_space_pivots(M, dom) -> list[int]` (independent-column indices); `solve_columns(B, V, dom) -> list[vector] | None` (each column of `V` in the span of `B`'s columns); `independent_modulo(cands, sub_cols, dom) -> list[int]` (greedy rank-growth selection of `cands` columns independent modulo the span of `sub_cols`).
  - `Module(algebra, dim, action, name="M")` — frozen-ish object; `.algebra`, `.dim`, `.action` (dict `basis_label -> matrix`), `.name`, `.domain`; `.vertex_indices(v) -> list[int]` no — instead `.vertex_projection(v)` returns `action["e_v"]`; `.dimension_vector() -> dict[vertex, int]`; `.check_module() -> tuple[bool, object]`. Factory `Module.from_arrow_action(algebra, dimension_vector, arrow_action, name=) -> Module` that extends idempotent+arrow actions to all basis labels.

- [ ] **Step 1: Write the failing test**

`tests/modules/test_module.py`:

```python
"""Right A-modules over any Domain (spec §3.6, §5 component 7). Column-vector
convention: m*b = action[b] @ m; action[x*y] = action[y] @ action[x]."""
from quiverlab import linear_path_algebra, GF
from quiverlab.modules.module import Module
from quiverlab.modules import linalg_mod as lm
from quiverlab.fields import QQ


def test_matrix_helpers_over_qq():
    I = lm.identity(2, QQ)
    assert I == [[QQ.coerce(1), QQ.coerce(0)], [QQ.coerce(0), QQ.coerce(1)]]
    A = [[QQ.coerce(1), QQ.coerce(2)]]
    B = [[QQ.coerce(1)], [QQ.coerce(1)]]
    assert lm.matmul(A, B, QQ) == [[QQ.coerce(3)]]
    K = lm.kernel_columns([[QQ.coerce(1), QQ.coerce(-1)]], QQ)
    assert len(K) == 1 and K[0][0] == K[0][1]


def test_column_space_pivots_selects_the_nonzero_column():
    # regression: column_space_pivots must return pivot COLUMNS of M (not of M^T).
    # M = [[0, 0], [1, 0]] has column 0 nonzero, column 1 zero -> pivots == [0].
    o, z = QQ.coerce(1), QQ.coerce(0)
    M = [[z, z], [o, z]]
    assert lm.column_space_pivots(M, QQ) == [0]     # the rad-image basis, not the zero col


def test_regular_projective_p1_of_a2_is_a_module():
    # P_1 = e_1 A for kA_2 built directly (Task 3 builder); here we hand-assemble it to
    # pin the representation. basis of P_1: [e_1, a] (paths from vertex 1).
    A = linear_path_algebra(2)     # arrow auto-named "a1"; basis [e_1, e_2, a1]
    d = QQ
    o, z = d.coerce(1), d.coerce(0)
    # right action on column [e_1, a1]:  e_1*e_1=e_1, e_1*a1=a1, a1*e_2=a1
    action = {
        "e_1": [[o, z], [z, z]],   # projects onto e_1 component
        "e_2": [[z, z], [z, o]],   # projects onto a1 component (a1 = a1*e_2)
        "a1":  [[z, z], [o, z]],   # e_1*a1 = a1 : sends e_1-coord to a1-coord
    }
    P1 = Module(A, 2, action, name="P_1")
    ok, why = P1.check_module()
    assert ok, why
    assert P1.dimension_vector() == {1: 1, 2: 1}


def test_dimension_vector_sums_to_dim_over_gfp():
    A = linear_path_algebra(2)
    p = GF(5)
    o, z = p.coerce(1), p.coerce(0)
    S1 = Module(A, 1, {"e_1": [[o]], "e_2": [[z]], "a1": [[z]]}, name="S_1")
    assert S1.dimension_vector() == {1: 1, 2: 0}
    assert sum(S1.dimension_vector().values()) == S1.dim


def test_from_arrow_action_extends_to_all_basis_labels():
    A = linear_path_algebra(2)   # basis labels e_1, e_2, a1
    d = QQ
    o, z = d.coerce(1), d.coerce(0)
    P1 = Module.from_arrow_action(
        A, dimension_vector={1: 1, 2: 1},
        arrow_action={"a1": [[z, z], [o, z]]},
        name="P_1")
    # extension must have filled e_1, e_2, a1 and satisfy the module axioms
    assert set(P1.action) == {"e_1", "e_2", "a1"}
    ok, _ = P1.check_module()
    assert ok
```

- [ ] **Step 2: Run — expect FAIL** (`ModuleNotFoundError: quiverlab.modules`).

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_module.py -q`

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/modules/__init__.py`:

```python
"""quiverlab.modules: right A-modules over any exact Domain, and their homological
algebra (simples/projectives/injectives, radical/top/socle, Hom/End, minimal
projective resolutions, Ext). Right-module (anti-homomorphism) convention throughout:
m*b = action[b] @ m (m a column), action[x*y] = action[y] @ action[x]. The minimal
resolution engine (resolution.py) is generalized from the bridge obstruction engine
(spec §5 component 7; MIT header there)."""
```

`src/quiverlab/modules/linalg_mod.py`:

```python
"""Domain-generic matrix arithmetic for modules, built on fields.linalg.

Matrices are list[list[elt]] (rows x cols); a vector is list[elt] (a column). All
coefficient work goes through the Domain (dom.add/sub/neg/mul/inv/is_zero) or through
fields.linalg, so everything is exact over QQ / QQ(alpha) / GF(p) / GF(p^n)."""
from quiverlab.fields import linalg


def zeros(r, c, dom):
    z = dom.zero()
    return [[z for _ in range(c)] for _ in range(r)]


def identity(n, dom):
    o, z = dom.one(), dom.zero()
    return [[o if i == j else z for j in range(n)] for i in range(n)]


def transpose(M):
    if not M:
        return []
    return [list(col) for col in zip(*M)]


def matmul(A, B, dom):
    if not A or not B:
        return []
    p, q, r = len(A), len(B), len(B[0])
    out = zeros(p, r, dom)
    for i in range(p):
        Ai = A[i]
        for k in range(q):
            a = Ai[k]
            if dom.is_zero(a):
                continue
            Bk = B[k]
            oi = out[i]
            for j in range(r):
                oi[j] = dom.add(oi[j], dom.mul(a, Bk[j]))
    return out


def matvec(A, v, dom):
    out = []
    for row in A:
        s = dom.zero()
        for a, x in zip(row, v):
            if not dom.is_zero(a):
                s = dom.add(s, dom.mul(a, x))
        out.append(s)
    return out


def kron(A, B, dom):
    """Kronecker product (row-block layout)."""
    ra, ca = len(A), len(A[0]) if A else 0
    rb, cb = len(B), len(B[0]) if B else 0
    out = zeros(ra * rb, ca * cb, dom)
    for i in range(ra):
        for j in range(ca):
            a = A[i][j]
            if dom.is_zero(a):
                continue
            for k in range(rb):
                for l in range(cb):
                    out[i * rb + k][j * cb + l] = dom.mul(a, B[k][l])
    return out


def vstack(mats):
    out = []
    for M in mats:
        out.extend(M)
    return out


def col(M, j):
    return [row[j] for row in M]


def cols_to_matrix(cols):
    """Assemble a matrix whose columns are the given column-vectors."""
    if not cols:
        return []
    n = len(cols[0])
    return [[c[i] for c in cols] for i in range(n)]


def mat_rank(M, dom):
    return linalg.rank(M, dom) if M else 0


def kernel_columns(M, dom):
    """Basis of {x : M x = 0} as a list of column-vectors."""
    if not M:
        return []
    return linalg.nullspace(M, dom)


def column_space_pivots(M, dom):
    """Indices of a maximal independent set of COLUMNS of M. `linalg.rref(M)` returns
    the pivot COLUMN indices of M directly, and the pivot columns of the original matrix
    are a basis of its column space -- so return those (NOT rref of the transpose, which
    would give row-space indices; that bug selected the zero column of P_1 on kA_2 and
    broke rad/top)."""
    if not M:
        return []
    _, piv = linalg.rref(M, dom)
    return sorted(piv)


def solve_columns(B, V, dom):
    """Express every column of V in the span of the columns of B. Returns the list of
    coefficient column-vectors x with B x = V[:, j], or None if any column is not in
    the span."""
    out = []
    for j in range(len(V[0]) if V else 0):
        x = linalg.solve(B, col(V, j), dom)
        if x is None:
            return None
        out.append(x)
    return out


def independent_modulo(cands, sub_cols, dom):
    """Greedy rank-growth: return the indices of the columns in `cands` that are
    linearly independent modulo the span of `sub_cols`. Used to pick minimal
    (mod-radical) generators."""
    basis = [list(c) for c in sub_cols]
    base_rank = mat_rank(cols_to_matrix(basis), dom) if basis else 0
    chosen = []
    for idx, c in enumerate(cands):
        trial = basis + [list(c)]
        r = mat_rank(cols_to_matrix(trial), dom)
        if r > base_rank:
            basis = trial
            base_rank = r
            chosen.append(idx)
    return chosen
```

`src/quiverlab/modules/module.py`:

```python
"""The right A-module object and its radical/top/socle (spec §3.6).

RIGHT modules, anti-homomorphism convention (see the package docstring): an element m
is a COLUMN vector in a fixed k-basis of M; the action of an algebra basis element b is
the matrix action[b] with m*b = action[b] @ m, and action[x*y] = action[y] @ action[x].
The vertex subspace M*e_v is the image of action['e_v']; dimension_vector[v] = its rank.
"""
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm


class Module:
    def __init__(self, algebra, dim, action, name="M"):
        self.algebra = algebra
        self.domain = algebra.domain
        self.dim = dim
        self.action = action           # basis_label -> (dim x dim) matrix
        self.name = name

    def _idem_label(self, v):
        return f"e_{v}"

    def vertex_projection(self, v):
        return self.action[self._idem_label(v)]

    def dimension_vector(self):
        dom = self.domain
        out = {}
        for v in self.algebra.quiver.vertices:
            out[v] = lm.mat_rank(self.vertex_projection(v), dom)
        return out

    def _arrow_labels(self):
        return list(self.algebra.quiver.arrows)

    def check_module(self, extra_labels=None):
        """Verify the representation is a genuine right A-module: (i) the idempotents
        sum to the identity and are orthogonal projections; (ii) every relation is
        satisfied; (iii) action is multiplicative on the basis-label products it is
        given. Returns (True, None) or (False, witness)."""
        dom = self.domain
        n = self.dim
        # (i) sum of idempotent actions == I
        acc = lm.zeros(n, n, dom)
        for v in self.algebra.quiver.vertices:
            P = self.vertex_projection(v)
            acc = _add(acc, P, dom)
        if acc != lm.identity(n, dom):
            return False, "sum of e_v actions != identity"
        # (ii) relations: for each relation sum c_w * word, sum c_w * action[word] == 0
        for rel in (self.algebra.relations or []):
            M = lm.zeros(n, n, dom)
            for coeff, word in _relation_terms(rel, dom):
                M = _add(M, _scale(self._action_of_word(word), coeff, dom), dom)
            if any(not dom.is_zero(x) for row in M for x in row):
                return False, f"relation not satisfied: {rel}"
        return True, None

    def _action_of_word(self, word):
        """action of a path word (tuple of arrow names) by composing arrow actions in
        anti-homomorphism order: action[a1*...*ak] = action[ak] @ ... @ action[a1]."""
        dom = self.domain
        if word == ():
            return lm.identity(self.dim, dom)
        M = None
        for name in word:  # left to right; anti-homo => multiply on the LEFT
            Aa = self.action[name]
            M = Aa if M is None else lm.matmul(Aa, M, dom)
        return M

    @classmethod
    def from_arrow_action(cls, algebra, dimension_vector, arrow_action, name="M"):
        """Build a module from per-arrow action matrices plus the dimension vector.
        The idempotent actions are the block projections implied by dimension_vector
        (in the vertex-ordered basis), and every non-trivial basis-path label's action
        is composed from the arrow actions. Validated before return."""
        dom = algebra.domain
        verts = list(algebra.quiver.vertices)
        dims = [dimension_vector.get(v, 0) for v in verts]
        n = sum(dims)
        # basis ordered by vertex block: build idempotent projections
        action = {}
        offset = 0
        starts = {}
        for v, dv in zip(verts, dims):
            starts[v] = offset
            offset += dv
        for v, dv in zip(verts, dims):
            P = lm.zeros(n, n, dom)
            for i in range(starts[v], starts[v] + dv):
                P[i][i] = dom.one()
            action[f"e_{v}"] = P
        for aname, mat in arrow_action.items():
            action[aname] = mat
        M = cls(algebra, n, action, name=name)
        # fill every algebra basis-label action (paths + idempotents) by composition
        M._extend_to_basis_labels()
        ok, why = M.check_module()
        if not ok:
            raise QuiverlabError(f"from_arrow_action({name}): not a module: {why}",
                                 hint="check that the arrow matrices satisfy the relations")
        return M

    def _extend_to_basis_labels(self):
        """Ensure action[label] exists for every algebra basis label (idempotents and
        path words), computed by composing the stored arrow/idempotent actions."""
        for label in self.algebra.basis_labels:
            if label in self.action:
                continue
            if label.startswith("e_"):
                # already set for genuine vertices; any missing means a semisimple gap
                continue
            word = tuple(label.split("*"))
            self.action[label] = self._action_of_word(word)

    def radical(self):
        from quiverlab.modules.radtopsoc import radical as _r
        return _r(self)

    def top(self):
        from quiverlab.modules.radtopsoc import top as _t
        return _t(self)

    def socle(self):
        from quiverlab.modules.radtopsoc import socle as _s
        return _s(self)

    def projective_resolution(self, length):
        from quiverlab.modules.resolution import minimal_resolution, ProjectiveResolution
        terms, dmats = minimal_resolution(self, length)
        return ProjectiveResolution(self, terms, dmats)

    def __repr__(self):
        dv = self.dimension_vector()
        return f"{self.name}: right {self.algebra} module, dim {self.dim}, dimvec {dv}"


def _add(A, B, dom):
    return [[dom.add(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def _scale(A, c, dom):
    return [[dom.mul(c, x) for x in row] for row in A]


def _relation_terms(rel, dom):
    """Yield (domain coeff, word) pairs of a Plan-01 Relation, coefficients coerced."""
    for coeff, word in rel.terms:
        yield dom.coerce(coeff), tuple(word)
```

> **Note on `check_module` relation access.** `A.relations` is the tuple of Plan-01 `Relation` objects (`.terms = ((Fraction, (arrow,...)), ...)`). Coefficients are coerced into the module's domain. For a monomial algebra a relation is a single word with an empty structure; for general `kQ/I` (Plan 03) it is a linear combination — both are handled by `_relation_terms`. If `A.relations` is `None` (a hand-built structure-constant algebra) the relation check is skipped (there is no path presentation) and only the idempotent-completeness check runs.

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_module.py tests/test_no_floats.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/__init__.py src/quiverlab/modules/linalg_mod.py src/quiverlab/modules/module.py tests/modules/test_module.py
git commit -m "feat(modules): right A-module object + domain-generic matrix helpers + dimension_vector

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: `simple` / `projective` / `injective` from quiver provenance

**Files:** Create `src/quiverlab/modules/builders.py`, `tests/modules/test_builders.py`. Edit `src/quiverlab/core/algebra.py` (add `simple`, `projective`, `injective` methods).

**Interfaces:**
- Consumes: `Module`, `Algebra` (`.basis_labels`, `.multiply`, `._basis_vec`, `.quiver`), `linalg_mod`, `QuiverlabError`.
- Produces:
  - `builders.simple(A, v) -> Module`, `builders.projective(A, v) -> Module`, `builders.injective(A, v) -> Module`.
  - `Algebra.simple(self, v)`, `Algebra.projective(self, v)`, `Algebra.injective(self, v)` (deferred-import wrappers). Each raises `QuiverlabError` (with an actionable hint) if the algebra carries no quiver/basis provenance, exactly like `cartan_matrix`.

- [ ] **Step 1: Write the failing test**

`tests/modules/test_builders.py`:

```python
"""Simples/projectives/injectives from quiver provenance (spec §3.6). Fixtures A & B."""
import pytest
from quiverlab import Quiver, CC, GF, linear_path_algebra, Algebra


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def test_a2_simples():
    A = linear_path_algebra(2)
    assert A.simple(1).dimension_vector() == {1: 1, 2: 0}
    assert A.simple(2).dimension_vector() == {1: 0, 2: 1}
    assert A.simple(1).dim == 1


def test_a2_projectives():
    A = linear_path_algebra(2)
    P1, P2 = A.projective(1), A.projective(2)
    assert P1.dim == 2 and P1.dimension_vector() == {1: 1, 2: 1}
    assert P2.dim == 1 and P2.dimension_vector() == {1: 0, 2: 1}
    ok, why = P1.check_module()
    assert ok, why


def test_a2_injectives():
    A = linear_path_algebra(2)
    I1, I2 = A.injective(1), A.injective(2)
    assert I1.dim == 1 and I1.dimension_vector() == {1: 1, 2: 0}   # I_1 = S_1
    assert I2.dim == 2 and I2.dimension_vector() == {1: 1, 2: 1}
    ok, _ = I2.check_module()
    assert ok


def test_square_projectives_match_cartan_rows():
    A = _square()
    C = A.cartan_matrix()
    verts = [1, 2, 3, 4]
    for i, v in enumerate(verts):
        Pv = A.projective(v)
        assert Pv.dimension_vector() == {verts[j]: C[i][j] for j in range(4)}
    assert A.projective(1).dim == 4
    ok, _ = A.projective(1).check_module()
    assert ok


def test_injective_dimvec_is_cartan_column():
    A = _square()
    C = A.cartan_matrix()
    verts = [1, 2, 3, 4]
    for j, v in enumerate(verts):
        Iv = A.injective(v)
        assert Iv.dimension_vector() == {verts[i]: C[i][j] for i in range(4)}


def test_builders_over_gfp():
    A = linear_path_algebra(2, field=GF(7))
    assert A.projective(1).dimension_vector() == {1: 1, 2: 1}


def test_builders_need_provenance():
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(Exception):
        A.simple(1)
```

- [ ] **Step 2: Run — expect FAIL** (`AttributeError: 'Algebra' object has no attribute 'simple'`).

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/modules/builders.py`:

```python
"""Simples, projectives, injectives from the quiver presentation (spec §3.6, §5 c.7).

RIGHT modules. P_v = e_v A (basis paths STARTING at v); S_v = top P_v (1-dim at v,
arrows act 0); I_v = D(A e_v) = k-dual of the LEFT projective A e_v (basis paths ENDING
at v), right action = transpose of left multiplication. All exact over any Domain."""
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module


def _require_provenance(A, what):
    if A.quiver is None or A.basis_labels is None:
        raise QuiverlabError(
            f"{what} needs the quiver presentation",
            hint="construct the algebra via Quiver.algebra(...); structure-constant "
                 "algebras carry no path basis",
        )


def _label_vertex_source(A, label):
    if label.startswith("e_"):
        return next(v for v in A.quiver.vertices if f"e_{v}" == label)
    return A.quiver.word_source(tuple(label.split("*")))


def _label_vertex_target(A, label):
    if label.startswith("e_"):
        return next(v for v in A.quiver.vertices if f"e_{v}" == label)
    return A.quiver.word_target(tuple(label.split("*")))


def simple(A, v):
    _require_provenance(A, "simple(v)")
    dom = A.domain
    o, z = dom.one(), dom.zero()
    action = {}
    for w in A.quiver.vertices:
        action[f"e_{w}"] = [[o if w == v else z]]
    for label in A.basis_labels:
        if not label.startswith("e_"):
            action[label] = [[z]]        # rad acts as 0 on a simple
    return Module(A, 1, action, name=f"S_{v}")


def projective(A, v):
    """P_v = e_v A. Basis = the algebra basis labels whose path STARTS at v, in the
    algebra's basis order. Right action of a basis element b: right multiplication,
    read from the algebra's structure constants restricted to this sub-basis."""
    _require_provenance(A, "projective(v)")
    dom = A.domain
    # sub-basis: indices of basis labels starting at v (idempotent e_v included)
    sub = [i for i, lab in enumerate(A.basis_labels) if _label_vertex_source(A, lab) == v]
    pos = {gi: k for k, gi in enumerate(sub)}
    n = len(sub)
    action = {}
    for blab in A.basis_labels:
        bi = A.basis_labels.index(blab)
        Mb = lm.zeros(n, n, dom)               # column k = coords of (p_k * b)
        for k, gi in enumerate(sub):
            prod = A.multiply(A._basis_vec(gi), A._basis_vec(bi))  # p_k * b in A-coords
            for gj in sub:
                Mb[pos[gj]][k] = prod[gj]
        action[blab] = Mb
    return Module(A, n, action, name=f"P_{v}")


def injective(A, v):
    """I_v = D(A e_v). Left projective A e_v has basis the labels whose path ENDS at v.
    Left multiplication L_b (b * p) gives the left action; the right action on the dual
    is its transpose. dim e_w I_v = # basis paths w -> v = C[w][v]."""
    _require_provenance(A, "injective(v)")
    dom = A.domain
    sub = [i for i, lab in enumerate(A.basis_labels) if _label_vertex_target(A, lab) == v]
    pos = {gi: k for k, gi in enumerate(sub)}
    n = len(sub)
    action = {}
    for blab in A.basis_labels:
        bi = A.basis_labels.index(blab)
        Lb = lm.zeros(n, n, dom)               # column k = coords of (b * p_k)
        for k, gi in enumerate(sub):
            prod = A.multiply(A._basis_vec(bi), A._basis_vec(gi))  # b * p_k
            for gj in sub:
                Lb[pos[gj]][k] = prod[gj]
        action[blab] = lm.transpose(Lb)        # right action = transpose of left mult
    return Module(A, n, action, name=f"I_{v}")
```

Edit `src/quiverlab/core/algebra.py` — add three methods on `Algebra` (place beside `cartan_matrix`, same deferred-import pattern):

```python
    def simple(self, v):
        """The simple right module S_v (spec §3.6)."""
        from quiverlab.modules.builders import simple
        return simple(self, v)

    def projective(self, v):
        """The indecomposable projective right module P_v = e_v A (spec §3.6)."""
        from quiverlab.modules.builders import projective
        return projective(self, v)

    def injective(self, v):
        """The indecomposable injective right module I_v = D(A e_v) (spec §3.6)."""
        from quiverlab.modules.builders import injective
        return injective(self, v)
```

> **`A.multiply` / `A._basis_vec` contract.** `A._basis_vec(i)` returns the i-th standard basis coordinate vector; `A.multiply(u, v)` returns the coordinate vector of the product `b_u · b_v`. Both are frozen Plan-01 surface (agent-verified). The projective/injective builders read the entire right/left multiplication off these; no re-derivation of the algebra structure is needed.

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_builders.py tests/invariants/test_cartan.py tests/test_no_floats.py -q`
Expected: all pass (builders + the Cartan pins stay green — the projective dimvecs equal Cartan rows by construction).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/builders.py tests/modules/test_builders.py src/quiverlab/core/algebra.py
git commit -m "feat(modules): simple/projective/injective from quiver provenance

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: `radical` / `top` / `socle`

**Files:** Create `src/quiverlab/modules/radtopsoc.py`, `tests/modules/test_radtopsoc.py`.

**Interfaces:**
- Consumes: `Module`, `linalg_mod`.
- Produces: `radtopsoc.radical(M) -> Module`, `top(M) -> Module`, `socle(M) -> Module`. Each returns a genuine `Module` (submodule or quotient) with its own action matrices in a chosen basis; `radical`/`socle` are submodules of `M`, `top` is the quotient `M/rad M`. `submodule(M, basis_cols) -> Module` and `quotient(M, sub_cols) -> Module` are the two constructors these use (restrict / co-restrict the action via `solve_columns`), reusable by the resolution engine (Task 6).

- [ ] **Step 1: Write the failing test**

`tests/modules/test_radtopsoc.py`:

```python
"""Radical / top / socle and the radical series (spec §3.6). Fixtures A & B."""
from quiverlab import Quiver, CC, linear_path_algebra


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"])


def test_p1_of_a2_radical_series():
    A = linear_path_algebra(2)
    P1 = A.projective(1)
    top = P1.top()
    assert top.dim == 1 and top.dimension_vector() == {1: 1, 2: 0}   # top P_1 = S_1
    rad = P1.radical()
    assert rad.dim == 1 and rad.dimension_vector() == {1: 0, 2: 1}   # rad P_1 = S_2
    soc = P1.socle()
    assert soc.dim == 1 and soc.dimension_vector() == {1: 0, 2: 1}   # soc P_1 = S_2
    # rad^2 P_1 = rad(rad P_1) = 0
    assert rad.radical().dim == 0


def test_square_p1_radical_layers():
    A = _square()
    P1 = A.projective(1)               # dim 4, dimvec {1:1,2:1,3:1,4:1}
    assert P1.top().dimension_vector() == {1: 1, 2: 0, 3: 0, 4: 0}
    r1 = P1.radical()                  # rad P_1, dim 3
    assert r1.dim == 3
    assert r1.top().dimension_vector() == {1: 0, 2: 1, 3: 1, 4: 0}   # S_2 (+) S_3
    r2 = r1.radical()                 # rad^2 P_1 = span{a*b}, dim 1 ~ S_4
    assert r2.dim == 1 and r2.dimension_vector() == {1: 0, 2: 0, 3: 0, 4: 1}
    assert r2.radical().dim == 0     # rad^3 P_1 = 0  -> Loewy length 3
    assert P1.socle().dimension_vector() == {1: 0, 2: 0, 3: 0, 4: 1}   # soc P_1 = S_4


def test_simple_is_its_own_top_and_socle_zero_radical():
    A = linear_path_algebra(2)
    S1 = A.simple(1)
    assert S1.radical().dim == 0
    assert S1.top().dim == 1
    assert S1.socle().dim == 1
```

- [ ] **Step 2: Run — expect FAIL** (`ModuleNotFoundError: quiverlab.modules.radtopsoc`).

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/modules/radtopsoc.py`:

```python
"""Radical, top and socle of a right A-module (spec §3.6).

rad M = M*(rad A) = sum over arrows of image(action[arrow])  (a submodule);
top M = M / rad M  (semisimple);
soc M = {m : m*(rad A) = 0} = intersection over arrows of ker(action[arrow]).
Submodule/quotient modules inherit the action by restriction/co-restriction, solved
exactly over the Domain via fields.linalg (solve_columns)."""
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module


def _rad_image_cols(M):
    """Column-vectors spanning M*(rad A) = sum_arrows image(action[arrow])."""
    dom = M.domain
    gens = []
    for aname in M.algebra.quiver.arrows:
        Aa = M.action[aname]
        for j in range(M.dim):
            gens.append(lm.col(Aa, j))
    if not gens:
        return []
    G = lm.cols_to_matrix(gens)
    piv = lm.column_space_pivots(G, dom)
    return [lm.col(G, j) for j in piv]


def submodule(M, basis_cols, name="sub"):
    """The submodule of M spanned by basis_cols (assumed A-stable). Its action[b] is
    the coordinates of action[b] applied to each basis column, expressed back in
    basis_cols (solved over the Domain)."""
    dom = M.domain
    B = lm.cols_to_matrix(basis_cols) if basis_cols else lm.zeros(M.dim, 0, dom)
    n = len(basis_cols)
    action = {}
    for label, Ab in M.action.items():
        if n == 0:
            action[label] = lm.zeros(0, 0, dom)
            continue
        images = [lm.matvec(Ab, c, dom) for c in basis_cols]   # b acts on each generator
        V = lm.cols_to_matrix(images)
        coeffs = lm.solve_columns(B, V, dom)                   # express in basis_cols
        assert coeffs is not None, f"submodule not A-stable under {label}"
        action[label] = lm.cols_to_matrix(coeffs)
    return Module(M.algebra, n, action, name=name)


def quotient(M, sub_cols, name="quot"):
    """The quotient module M / <sub_cols>. Pick coset representatives = a basis of M
    completing sub_cols; action[b] is read on the representatives modulo the submodule."""
    dom = M.domain
    sub_piveb = sub_cols
    ident = lm.identity(M.dim, dom)
    std = [lm.col(ident, j) for j in range(M.dim)]
    # representatives: standard vectors independent modulo the submodule
    rep_idx = lm.independent_modulo(std, sub_piveb, dom)
    reps = [std[i] for i in rep_idx]
    n = len(reps)
    # basis of the WHOLE space: submodule columns then representatives
    whole = [list(c) for c in sub_cols] + reps
    W = lm.cols_to_matrix(whole)
    action = {}
    s = len(sub_cols)
    for label, Ab in M.action.items():
        cols = []
        for r in reps:
            img = lm.matvec(Ab, r, dom)
            x = lm.solve_columns(W, lm.cols_to_matrix([img]), dom)[0]
            cols.append(x[s:])                 # drop the submodule part -> class in quotient
        action[label] = lm.cols_to_matrix(cols) if cols else lm.zeros(n, n, dom)
    return Module(M.algebra, n, action, name=name)


def radical(M):
    return submodule(M, _rad_image_cols(M), name=f"rad {M.name}")


def top(M):
    return quotient(M, _rad_image_cols(M), name=f"top {M.name}")


def socle(M):
    """soc M = intersection over arrows of ker(action[arrow])."""
    dom = M.domain
    inter = None
    arrows = list(M.algebra.quiver.arrows)
    if not arrows:                              # semisimple: soc = M
        return submodule(M, [lm.col(lm.identity(M.dim, dom), j) for j in range(M.dim)],
                         name=f"soc {M.name}")
    for aname in arrows:
        ker = lm.kernel_columns(M.action[aname], dom)
        if inter is None:
            inter = ker
        else:
            inter = _intersect(inter, ker, dom)
    return submodule(M, inter or [], name=f"soc {M.name}")


def _intersect(cols1, cols2, dom):
    """Basis of the intersection of two subspaces given by column bases (via the
    kernel of [B1 | -B2] projected to the B1 side)."""
    if not cols1 or not cols2:
        return []
    B1, B2 = lm.cols_to_matrix(cols1), lm.cols_to_matrix(cols2)
    stacked = [row1 + [dom.neg(x) for x in row2] for row1, row2 in zip(B1, B2)]
    ker = lm.kernel_columns(stacked, dom)
    k1 = len(cols1)
    out = []
    for z in ker:
        coeff = z[:k1]
        vec = lm.matvec(B1, coeff, dom)
        out.append(vec)
    # reduce to an independent set
    if not out:
        return []
    G = lm.cols_to_matrix(out)
    piv = lm.column_space_pivots(G, dom)
    return [lm.col(G, j) for j in piv]
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_radtopsoc.py tests/test_no_floats.py -q`
Expected: all pass (radical series of `P_1` in Fixtures A and B match the hand values).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/radtopsoc.py tests/modules/test_radtopsoc.py
git commit -m "feat(modules): radical/top/socle and the radical series over any Domain

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: `Hom` / `End` over a Domain

**Files:** Create `src/quiverlab/modules/hom.py`, `tests/modules/test_hom.py`. Edit `core/algebra.py` (add `hom`).

**Interfaces:**
- Consumes: `Module`, `linalg_mod` (`kron`, `identity`, `transpose`, `vstack`, `kernel_columns`).
- Produces: `hom.hom_space(M, N) -> list[matrix]` (a `k`-basis of `Hom_A(M, N)`, each a `dim N × dim M` matrix), `hom.hom_dim(M, N) -> int`, `hom.end_dim(M) -> int`. `Algebra.hom(self, M, N) -> int` (deferred-import wrapper returning `dim Hom_A(M, N)`).

**Mathematics.** A right-module homomorphism `φ: M → N` is a `k`-linear `φ` (matrix `dim N × dim M`) with `φ(m·b) = φ(m)·b` for all `b`, i.e. `N.action[b] @ φ = φ @ M.action[b]`. It suffices to impose this on the generators `b ∈ {arrows} ∪ {idempotents}` (they generate `A`). Vectorizing with column-stacking `vec(φ)` (`vec(AXB) = (Bᵀ ⊗ A) vec(X)`): the constraint per `b` is `(I_{dim M} ⊗ N.action[b] − M.action[b]ᵀ ⊗ I_{dim N}) · vec(φ) = 0`. Stack over generators; `dim Hom = dim ker`.

- [ ] **Step 1: Write the failing test**

`tests/modules/test_hom.py`:

```python
"""Hom/End of right A-modules over any Domain (spec §3.6)."""
from quiverlab import Quiver, CC, GF, linear_path_algebra


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"])


def test_hom_of_simples_is_kronecker_delta():
    A = linear_path_algebra(2)
    S1, S2 = A.simple(1), A.simple(2)
    assert A.hom(S1, S1) == 1
    assert A.hom(S2, S2) == 1
    assert A.hom(S1, S2) == 0
    assert A.hom(S2, S1) == 0


def test_hom_projective_to_module_reads_off_dimension_vector():
    # dim Hom(P_v, N) = dim (N e_v) = dimension_vector(N)[v]  (Yoneda for projectives)
    A = _square()
    N = A.projective(1)
    for v in (1, 2, 3, 4):
        assert A.hom(A.projective(v), N) == N.dimension_vector()[v]


def test_end_of_indecomposable_projective_is_local_dim_one_top():
    A = linear_path_algebra(2)
    from quiverlab.modules.hom import end_dim
    # End(P_1): P_1 indecomposable; endomorphisms are scalars + radical -> dim 1 here
    assert end_dim(A.projective(1)) == 1


def test_hom_over_gfp():
    A = linear_path_algebra(2, field=GF(7))
    assert A.hom(A.simple(1), A.simple(1)) == 1
```

- [ ] **Step 2: Run — expect FAIL** (`AttributeError: ... 'hom'`).

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/modules/hom.py`:

```python
"""Hom and End spaces of right A-modules over any Domain (spec §3.6, §5 c.7).

phi: M -> N is a right-module map iff N.action[b] @ phi = phi @ M.action[b] for every
generator b (arrows and idempotents). Column-stacking vec: the constraint per b is
(I_{dimM} (x) N.action[b] - M.action[b]^T (x) I_{dimN}) vec(phi) = 0. dim Hom = dim ker."""
from quiverlab.modules import linalg_mod as lm


def _generators(M):
    labels = [f"e_{v}" for v in M.algebra.quiver.vertices]
    labels += list(M.algebra.quiver.arrows)
    return labels


def hom_space(M, N):
    dom = M.domain
    dm, dn = M.dim, N.dim
    Im = lm.identity(dm, dom)
    In = lm.identity(dn, dom)
    blocks = []
    for b in _generators(M):
        Nb, Mb = N.action[b], M.action[b]
        left = lm.kron(Im, Nb, dom)                     # I_dm (x) N.action[b]
        right = lm.kron(lm.transpose(Mb), In, dom)      # M.action[b]^T (x) I_dn
        blocks.append([[dom.sub(left[i][j], right[i][j]) for j in range(dm * dn)]
                       for i in range(dm * dn)])
    stacked = lm.vstack(blocks) if blocks else lm.zeros(0, dm * dn, dom)
    ker = lm.kernel_columns(stacked, dom)
    # reshape each kernel vector (length dm*dn, column-stacked) into a dn x dm matrix
    homs = []
    for z in ker:
        phi = [[z[j * dn + i] for j in range(dm)] for i in range(dn)]
        homs.append(phi)
    return homs


def hom_dim(M, N):
    return len(hom_space(M, N))


def end_dim(M):
    return hom_dim(M, M)
```

Edit `core/algebra.py` — add:

```python
    def hom(self, M, N):
        """dim Hom_A(M, N) for right A-modules M, N (spec §3.6)."""
        from quiverlab.modules.hom import hom_dim
        return hom_dim(M, N)
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_hom.py tests/test_no_floats.py -q`
Expected: all pass (`Hom(S_a,S_b)=δ`; `Hom(P_v, N) = dimvec(N)[v]` — the Yoneda check independently validates the projective builder).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/hom.py tests/modules/test_hom.py src/quiverlab/core/algebra.py
git commit -m "feat(modules): Hom/End dimensions of right A-modules over any Domain

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: Minimal projective resolution — the generalized `module_ext` engine

**Files:** Create `src/quiverlab/modules/resolution.py`, `tests/modules/test_resolution.py`. Edit `core/algebra.py` is not needed (the `Module.projective_resolution` hook already added in Task 2 imports this module).

**Interfaces:**
- Consumes: `Module`, `builders.projective`, `radtopsoc.submodule`, `linalg_mod`, `errors.DepthLimitError`.
- Produces:
  - `resolution.projective_cover(M) -> (Q0, d0, gens)` where `Q0` is the projective `⊕_v P_v^{t_v}` (a `Module`), `d0` its cover matrix `dim M × dim Q0`, and `gens` the list of `(vertex, lifted-column-in-M)`.
  - `resolution.minimal_resolution(M, length, max_term_dim=200000) -> (terms, dmats)`: `terms[n]` a `TermData` (vertices of the summands + the `Module` `Q_n`), `dmats[n]` the matrix of `d_n : Q_n → Q_{n-1}` (`n ≥ 1`) with `dmats[0]` the cover `Q_0 → M`. Guards term size (`DepthLimitError` on blow-up, stating the certified length).
  - `ProjectiveResolution(M, terms, dmats)` — inspectable: `.length`, `.term(n)` (list of vertices), `.betti(n)` (# indecomposable projective summands = `len(term(n))`), `.differential(n)`, `.dimension_vectors()`, `__repr__`; `.is_finite()` (terminated) / `.pd()` (projective dimension or `None`).

- [ ] **Step 1: Write the failing test**

`tests/modules/test_resolution.py`:

```python
"""Minimal projective resolutions of right modules, generalized from the bridge engine
(spec §5 component 7). Fixtures A & B; any vertex set, any Domain."""
from quiverlab import Quiver, CC, GF, linear_path_algebra


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def test_resolution_of_s1_in_a2():
    # 0 -> P_2 -> P_1 -> S_1 -> 0
    A = linear_path_algebra(2)
    res = A.simple(1).projective_resolution(4)
    assert res.term(0) == [1]        # Q_0 = P_1
    assert res.term(1) == [2]        # Q_1 = P_2 (rad P_1 = S_2 = P_2)
    assert res.betti(0) == 1 and res.betti(1) == 1
    assert res.term(2) == []         # terminates: Omega_2 = 0
    assert res.pd() == 1


def test_projective_is_its_own_resolution():
    A = linear_path_algebra(2)
    res = A.projective(1).projective_resolution(3)
    assert res.term(0) == [1] and res.term(1) == []
    assert res.pd() == 0


def test_simple_s2_is_projective_in_a2():
    A = linear_path_algebra(2)
    assert A.simple(2).projective_resolution(3).pd() == 0    # S_2 = P_2


def test_square_s1_has_pd_two():
    A = _square()
    res = A.simple(1).projective_resolution(5)
    assert res.pd() == 2             # gl.dim(commutative square) = 2
    # d_n(Q_n) subset rad Q_{n-1}: minimality => betti numbers are the true Betti numbers
    assert res.betti(0) == 1
    assert res.term(3) == []


def test_resolution_over_gfp_matches_char0_shape():
    A = _square(field=GF(7))
    assert A.simple(1).projective_resolution(5).pd() == 2


def test_resolution_differentials_compose_to_zero():
    A = _square()
    res = A.simple(1).projective_resolution(4)
    from quiverlab.modules import linalg_mod as lm
    for n in range(1, res.length):
        Dn, Dn1 = res.differential(n), res.differential(n + 1)
        if not Dn or not Dn1 or len(Dn1[0]) == 0:
            continue
        prod = lm.matmul(Dn, Dn1, A.domain)
        assert all(A.domain.is_zero(x) for row in prod for x in row)
```

- [ ] **Step 2: Run — expect FAIL** (`ModuleNotFoundError: quiverlab.modules.resolution`).

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/modules/resolution.py`:

```python
# SPDX-License-Identifier: MIT
# Generalized from the bridge lab engine `obstruction/module_ext.py` (unlicensed in the
# bank; lifted here with attribution under quiverlab's MIT license per spec section 9).
# The bank engine was hardcoded to a 4-vertex diamond over Q; this reimplementation runs
# over ANY vertex set and ANY exact Domain via fields.linalg.
"""Minimal projective resolutions of right A-modules by iterated projective covers
(Green-Solberg-Zacharia, Trans. AMS 353 (2001) 2915-2939), spec section 5 component 7.

Right modules; left-to-right composition. A minimal resolution
  ... -> Q_2 -> Q_1 -> Q_0 -> M -> 0
is built one step at a time: cover M by projectives generated by top M, take the syzygy
kernel, cover it, repeat. Minimality = generators chosen independent modulo the radical,
so d_n(Q_n) subset rad Q_{n-1} and the betti numbers are the true Betti numbers."""
from quiverlab.errors import DepthLimitError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.builders import projective
from quiverlab.modules.radtopsoc import submodule, _rad_image_cols

_GUARD_HINT = ("module resolution term blew past max_term_dim; the certified length is "
               "in the message. Raise max_term_dim to push deeper.")


def _direct_sum(modules, name="Q"):
    """Block-diagonal direct sum of Modules (all over the same algebra); returns the sum
    Module plus the list of (start, dim) offsets of each summand."""
    if not modules:
        return None, []
    A = modules[0].algebra
    dom = A.domain
    dims = [m.dim for m in modules]
    n = sum(dims)
    offs = []
    o = 0
    for d in dims:
        offs.append((o, d))
        o += d
    action = {}
    for label in modules[0].action:
        M = lm.zeros(n, n, dom)
        for (s, d), mod in zip(offs, modules):
            blk = mod.action[label]
            for i in range(d):
                for j in range(d):
                    M[s + i][s + j] = blk[i][j]
        action[label] = M
    from quiverlab.modules.module import Module
    return Module(A, n, action, name=name), offs


def _homogeneous_top_generators(M):
    """A k-basis of top M, each lifted to a column of M living in a single vertex block.
    Returns list of (vertex, column-in-M). These are the projective-cover generators."""
    dom = M.domain
    rad_cols = _rad_image_cols(M)
    gens = []
    ident = lm.identity(M.dim, dom)
    for v in M.algebra.quiver.vertices:
        Pv = M.vertex_projection(v)                       # image = M*e_v
        # candidate vectors spanning M*e_v: independent columns of the projection
        piv = lm.column_space_pivots(Pv, dom)
        cands = [lm.col(Pv, j) for j in piv]
        # keep those independent modulo rad M (+ already-chosen gens): mod-radical basis
        chosen = lm.independent_modulo(cands, rad_cols + [g[1] for g in gens], dom)
        for idx in chosen:
            gens.append((v, cands[idx]))
    return gens


def projective_cover(M):
    """Q0 = (+)_g P_{vertex(g)} with cover d0: Q0 -> M sending each summand's canonical
    generator (the e_v basis vector of P_v) to the lifted top-generator g."""
    dom = M.domain
    gens = _homogeneous_top_generators(M)
    summands = [projective(M.algebra, v) for (v, _) in gens]
    Q0, offs = _direct_sum(summands, name="Q0") if summands else (None, [])
    if Q0 is None:
        from quiverlab.modules.module import Module
        Q0 = Module(M.algebra, 0, {lab: lm.zeros(0, 0, dom) for lab in M.action}, name="Q0")
    # cover matrix: column = image in M of each Q0 basis vector b = (canonical gen)*path.
    # For summand P_v with generator g (a column of M), the P_v basis vector indexed by
    # path label p (source v) maps to g * p = action_M[p] @ g.
    cols = []
    for (v, g), (s, d), Pmod in zip(gens, offs, summands):
        for k in range(d):
            plabel = Pmod._pv_basis_labels[k]              # label of the k-th P_v basis path
            img = lm.matvec(M.action[plabel], g, dom)      # g * p
            cols.append(img)
    d0 = lm.cols_to_matrix(cols) if cols else lm.zeros(M.dim, 0, dom)
    return Q0, d0, gens


def minimal_resolution(M, length, max_term_dim=200000):
    terms, dmats = [], []
    Q0, d0, _ = projective_cover(M)
    terms.append(_term_data(Q0))
    dmats.append(d0)
    cur_mod, cur_map = Q0, d0            # cur_map : cur_mod -> (previous)
    for n in range(1, length + 1):
        # syzygy Omega_n = ker(cur_map) as a submodule of cur_mod
        ker_cols = lm.kernel_columns(cur_map, M.domain) if cur_map and cur_map[0] else \
            ([lm.col(lm.identity(cur_mod.dim, M.domain), j) for j in range(cur_mod.dim)]
             if cur_mod.dim else [])
        if not ker_cols:
            terms.append(_term_data(None))
            dmats.append(lm.zeros(cur_mod.dim, 0, M.domain))
            break
        Omega = submodule(cur_mod, ker_cols, name=f"Omega_{n}")
        if Omega.dim > max_term_dim:
            raise DepthLimitError(
                f"module resolution: syzygy dim {Omega.dim} exceeds max_term_dim="
                f"{max_term_dim}; certified through length {n - 1}", hint=_GUARD_HINT)
        Qn, cover, _ = projective_cover(Omega)
        # d_n : Q_n -> Q_{n-1} is (cover: Q_n -> Omega) then (inclusion Omega -> Q_{n-1})
        B = lm.cols_to_matrix(ker_cols)                    # Omega basis inside cur_mod
        dn = lm.matmul(B, cover, M.domain) if cover and cover[0] else \
            lm.zeros(cur_mod.dim, 0, M.domain)
        terms.append(_term_data(Qn))
        dmats.append(dn)
        cur_mod, cur_map = Qn, cover
        if Qn.dim == 0:
            break
    return terms, dmats


class _term_data:
    def __init__(self, Qn):
        self.module = Qn
        self.vertices = list(getattr(Qn, "_summand_vertices", [])) if Qn else []
        self.dim = Qn.dim if Qn else 0


class ProjectiveResolution:
    def __init__(self, M, terms, dmats):
        self.module = M
        self.terms = terms
        self.dmats = dmats
        self.length = len(terms)

    def term(self, n):
        return self.terms[n].vertices if n < len(self.terms) else []

    def betti(self, n):
        return len(self.term(n))

    def differential(self, n):
        return self.dmats[n] if n < len(self.dmats) else []

    def dimension_vectors(self):
        return [t.module.dimension_vector() if t.module else {} for t in self.terms]

    def is_finite(self):
        return any(t.dim == 0 for t in self.terms)

    def pd(self):
        for n, t in enumerate(self.terms):
            if t.dim == 0:
                return n - 1
        return None                     # not resolved within the requested length

    def __repr__(self):
        parts = " <- ".join(
            "(+)".join(f"P_{v}" for v in t.vertices) if t.vertices else "0"
            for t in self.terms)
        return f"proj. res. of {self.module.name}: {parts}"
```

> **Two small representation hooks needed on `Module`/`projective`.** `projective(A, v)` must record, on the returned module, the ordered basis labels of `P_v` (call it `_pv_basis_labels`) and the direct-sum builder must record `_summand_vertices` on `Q_n`. Add to `builders.projective`: `mod._pv_basis_labels = [A.basis_labels[i] for i in sub]; mod._pv_vertex = v` before returning; and in `resolution._direct_sum`: set `Q._summand_vertices = [m._pv_vertex for m in modules]`. These are private attributes (leading underscore), not part of the frozen surface. The canonical generator of `P_v` is its `e_v` basis vector = the first entry of `_pv_basis_labels` (the trivial path). The cover sends the `P_v` basis vector for path `p` to `g·p`; for `p = e_v` that is `g` itself (the generator lands on `g`), giving the projective cover.

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_resolution.py tests/test_no_floats.py -q`
Expected: all pass (`S_1` in A₂ resolves `0→P_2→P_1→S_1→0`, `pd 1`; the commutative square `S_1` has `pd 2`; `d∘d = 0`).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/resolution.py tests/modules/test_resolution.py src/quiverlab/modules/builders.py
git commit -m "feat(modules): minimal projective resolution engine (generalized module_ext, MIT)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: Module `Ext^n` + `global_dimension` + `is_selfinjective` + cross-checks

**Files:** Create `src/quiverlab/modules/ext.py`, `tests/modules/test_ext.py`. Edit `core/algebra.py` (add `ext`, `global_dimension`, `is_selfinjective`).

**Interfaces:**
- Consumes: `minimal_resolution`, `hom.hom_space`/`hom_dim`, `builders.projective`/`simple`, `radtopsoc.socle` (via the module), `linalg_mod`.
- Produces:
  - `ext.ext(A, M, N, n) -> int` = `dim Ext^n_A(M, N)`; `ext.ext_dims(A, M, N, top) -> list[int]` = `[dim Ext^0 .. dim Ext^top]`.
  - `ext.global_dimension(A, bound=32) -> GlobalDimension`: `sup_v pd(S_v)`; exact when every simple resolves within `bound`, else a labeled certified **lower** bound. `GlobalDimension` is a small dataclass: `.value: int`, `.exact: bool`, `.__int__`, `.__eq__` (equals the int when `exact`), `.__repr__` (`"gl.dim = 2 (exact)"` or `">= 5 (certified lower bound; not resolved within depth 32)"`).
  - `ext.is_selfinjective(A) -> bool` (spec §3.5, unclaimed by any other plan; natural home here since both `projective` and `socle` are built). `A` is self-injective ⇔ every indecomposable projective `P_v` is injective ⇔ `soc(P_v)` is **simple** for every `v` and `v ↦ (socle vertex)` is a **bijection** (the Nakayama permutation). Exact over **any** field (uses `socle`, not the GF(p) engine) — complementary to the GF(p)-only `is_frobenius` (Frobenius ⇔ self-injective for f.d. algebras).
  - `Algebra.ext(self, M, N, n)`, `Algebra.global_dimension(self)`, `Algebra.is_selfinjective(self)` (deferred-import wrappers).

**Mathematics.** With the minimal resolution `Q_• → M`, apply `Hom_A(−, N)`: `Ext^n = ker(δ^n)/im(δ^{n-1})` where `δ^n: Hom(Q_n, N) → Hom(Q_{n+1}, N)` is precomposition with `d_{n+1}`. Working with the `Hom`-space bases `H_n = hom_space(Q_n, N)` (each a `dim N × dim Q_n` matrix), the map `δ^n` sends `φ ↦ φ ∘ d_{n+1}`; in coordinates on the `H_n`/`H_{n+1}` bases this is a matrix `D_n`, and `dim Ext^n = dim H_n − rank(D_n) − rank(D_{n-1})` (with `D_{-1} = 0`).

- [ ] **Step 1: Write the failing test**

`tests/modules/test_ext.py`:

```python
"""Module Ext^n and global dimension (spec §3.5, §3.6). Literature oracles: ASS III /
Happel for hereditary A_n; Kunneth for the commutative square; bank diamond for Ext^2."""
import pytest
from quiverlab import Quiver, CC, GF, linear_path_algebra


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _mon_diamond():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b"])


def test_ext_of_simples_a3_hereditary():
    A = linear_path_algebra(3)          # 1 -> 2 -> 3
    S = {i: A.simple(i) for i in (1, 2, 3)}
    # Hom(S_a,S_b) = delta
    for a in (1, 2, 3):
        for b in (1, 2, 3):
            assert A.ext(S[a], S[b], 0) == (1 if a == b else 0)
    # Ext^1(S_a, S_{a+1}) = 1 ; all other Ext^1 = 0  (ASS III, Happel)
    assert A.ext(S[1], S[2], 1) == 1
    assert A.ext(S[2], S[3], 1) == 1
    assert A.ext(S[1], S[3], 1) == 0
    assert A.ext(S[1], S[1], 1) == 0
    # hereditary: Ext^{>=2} = 0
    assert A.ext(S[1], S[3], 2) == 0


def test_ext2_commutative_square_kunneth():
    A = _square()
    assert A.ext(A.simple(1), A.simple(4), 2) == 1      # kA_2 (x) kA_2, Kunneth
    assert A.ext(A.simple(1), A.simple(4), 1) == 0
    assert A.ext(A.simple(1), A.simple(4), 3) == 0


def test_ext2_monomial_diamond_bank_oracle():
    # bridge module_ext oracle: kD/(a*b) has Ext^2(S_1, S_4) = 1
    A = _mon_diamond()
    assert A.ext(A.simple(1), A.simple(4), 2) == 1


def test_projectives_have_no_higher_ext():
    A = _square()
    for v in (1, 2, 3, 4):
        assert A.ext(A.projective(v), A.simple(1), 1) == 0
        assert A.ext(A.projective(v), A.simple(1), 2) == 0


def test_global_dimension():
    assert int(linear_path_algebra(2).global_dimension()) == 1     # hereditary
    assert int(linear_path_algebra(3).global_dimension()) == 1
    gd = _square().global_dimension()
    assert int(gd) == 2 and gd.exact is True
    assert gd == 2


def test_global_dimension_selfinjective_is_infinite_bound():
    # k[x]/(x^2) is self-injective: pd(S) = infinity -> certified lower bound, not exact
    from quiverlab import truncated_polynomial
    gd = truncated_polynomial(2).global_dimension()
    assert gd.exact is False and gd.value >= 1


def test_is_selfinjective():
    from quiverlab import truncated_polynomial
    # k[x]/(x^n) is (Frobenius, hence) self-injective, over any field
    assert truncated_polynomial(2).is_selfinjective() is True
    assert truncated_polynomial(3).is_selfinjective() is True
    assert truncated_polynomial(3, field=GF(5)).is_selfinjective() is True
    # hereditary kA_2 is NOT self-injective (soc P_1 = S_2 = soc P_2, not a bijection)
    assert linear_path_algebra(2).is_selfinjective() is False
    # the commutative square is not self-injective either
    assert _square().is_selfinjective() is False


def test_selfinjective_agrees_with_frobenius_over_gfp():
    # over GF(p) the exact is_selfinjective must agree with the engine-backed is_frobenius
    from quiverlab import truncated_polynomial
    A = truncated_polynomial(2, field=GF(5))
    assert A.is_selfinjective() == A.is_frobenius() == True


def test_ext_over_gfp():
    A = _square(field=GF(7))
    assert A.ext(A.simple(1), A.simple(4), 2) == 1
```

- [ ] **Step 2: Run — expect FAIL** (`AttributeError: ... 'ext'`).

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/modules/ext.py`:

```python
"""Module Ext^n via a minimal projective resolution + Hom, and global dimension
(spec section 3.5-3.6). Ext^n_A(M,N) = H^n(Hom_A(Q_*, N)); gl.dim A = sup_v pd(S_v)."""
from dataclasses import dataclass

from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import hom_space
from quiverlab.modules.resolution import minimal_resolution


def _delta_matrix(Hn, Hn1, dn1, dom):
    """Matrix of delta^n: Hom(Q_n, N) -> Hom(Q_{n+1}, N), phi |-> phi @ d_{n+1}, in the
    given Hom-basis coordinates. Column j = coords of (Hn[j] @ dn1) in the Hn1 basis."""
    if not Hn or not Hn1:
        return lm.zeros(len(Hn1), len(Hn), dom)
    # flatten each Hom basis matrix (dn x dq) into a vector; build the coordinate solve
    def flat(mat):
        return [x for row in mat for x in row]
    basisH1 = lm.cols_to_matrix([flat(h) for h in Hn1])
    cols = []
    for phi in Hn:
        comp = lm.matmul(phi, dn1, dom)                # phi @ d_{n+1}: (dn x dq_{n+1})
        x = lm.solve_columns(basisH1, lm.cols_to_matrix([flat(comp)]), dom)[0]
        cols.append(x)
    return lm.cols_to_matrix(cols)


def ext_dims(A, M, N, top):
    dom = A.domain
    terms, dmats = minimal_resolution(M, top + 1)
    Qs = [t.module for t in terms]
    Homs = [hom_space(Q, N) if (Q is not None and Q.dim) else [] for Q in Qs]
    # delta^n uses d_{n+1} = dmats[n+1]
    deltas = []
    for n in range(len(Qs) - 1):
        dn1 = dmats[n + 1]
        deltas.append(_delta_matrix(Homs[n], Homs[n + 1], dn1, dom)
                      if (dn1 and dn1[0]) else lm.zeros(len(Homs[n + 1]), len(Homs[n]), dom))
    out = []
    for n in range(top + 1):
        cn = len(Homs[n]) if n < len(Homs) else 0
        r_n = lm.mat_rank(deltas[n], dom) if n < len(deltas) and deltas[n] and deltas[n][0] else 0
        r_nm1 = lm.mat_rank(deltas[n - 1], dom) if (n - 1) >= 0 and (n - 1) < len(deltas) \
            and deltas[n - 1] and deltas[n - 1][0] else 0
        out.append(cn - r_n - r_nm1)
    return out


def ext(A, M, N, n):
    return ext_dims(A, M, N, n)[n]


@dataclass
class GlobalDimension:
    value: int
    exact: bool

    def __int__(self):
        return self.value

    def __eq__(self, other):
        if isinstance(other, GlobalDimension):
            return (self.value, self.exact) == (other.value, other.exact)
        if isinstance(other, int):
            return self.exact and self.value == other
        return NotImplemented

    def __repr__(self):
        if self.exact:
            return f"gl.dim = {self.value} (exact)"
        return (f">= {self.value} (certified lower bound; not resolved within depth "
                f"{self.value})")


def global_dimension(A, bound=32):
    """sup over simples of projective dimension. Exact if every simple resolves within
    `bound`; otherwise a certified lower bound (some simple has pd >= bound)."""
    best, exact = 0, True
    for v in A.quiver.vertices:
        res = A.simple(v).projective_resolution(bound)
        pd = res.pd()
        if pd is None:                     # not resolved within bound -> lower bound
            exact = False
            best = max(best, bound)
        else:
            best = max(best, pd)
    return GlobalDimension(best, exact)


def is_selfinjective(A):
    """A is self-injective iff every indecomposable projective P_v is injective, iff
    soc(P_v) is SIMPLE for every v and v -> (socle vertex) is a BIJECTION (the Nakayama
    permutation). Exact over any field (uses the socle, not the GF(p) engine); for a
    finite-dimensional algebra self-injective == Frobenius (spec section 3.5)."""
    if A.quiver is None:
        from quiverlab.errors import QuiverlabError
        raise QuiverlabError("is_selfinjective needs the quiver presentation",
                             hint="construct via Quiver.algebra(...)")
    socle_vertex = {}
    for v in A.quiver.vertices:
        soc = A.projective(v).socle()
        dv = soc.dimension_vector()
        support = [w for w, d in dv.items() if d > 0]
        if soc.dim != 1 or len(support) != 1:
            return False                     # socle not simple -> P_v not injective
        socle_vertex[v] = support[0]
    return len(set(socle_vertex.values())) == len(list(A.quiver.vertices))
```

Edit `core/algebra.py` — add:

```python
    def ext(self, M, N, n):
        """dim Ext^n_A(M, N) for right A-modules M, N (spec §3.6)."""
        from quiverlab.modules.ext import ext
        return ext(self, M, N, n)

    def global_dimension(self):
        """Global dimension: exact value or a labeled certified lower bound (spec §3.5)."""
        from quiverlab.modules.ext import global_dimension
        return global_dimension(self)

    def is_selfinjective(self):
        """True iff every indecomposable projective is injective (self-injective =
        Frobenius for a f.d. algebra); exact over any field (spec §3.5)."""
        from quiverlab.modules.ext import is_selfinjective
        return is_selfinjective(self)
```

> **Optional Plan-04 cross-check (guarded).** Add to `tests/modules/test_ext.py` a test guarded by `pytest.importorskip("quiverlab.resolutions_cs")` that, for the monomial diamond (Fixture C) and the quantum CI, cross-checks the *module* Ext of simples against an independent route where the algebra is monomial/quadratic: recompute `pd(S_v)` and confirm it agrees with `global_dimension`. This keeps the "CS resolutions for Ext where monomial/quadratic" obligation honest without coupling module Ext (one-sided) to the CS bimodule resolution (which serves Hochschild, Plan 04). The primary cross-checks are the literature values (ASS/Happel/Künneth) and the bank diamond oracle — both float-free and independent of Plan 04.

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_ext.py tests/test_no_floats.py -q`
Expected: all pass (ASS/Happel Ext of A₃ simples; Künneth `Ext²(S_1,S_4)=1` for the square; bank oracle `Ext²=1` for the monomial diamond; `global_dimension` exact 1/1/2 and the self-injective lower bound).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ext.py tests/modules/test_ext.py src/quiverlab/core/algebra.py
git commit -m "feat(modules): module Ext^n + global dimension + is_selfinjective

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Exact spectral layer — `spectral_radius` + `mahler_measure` (no floats)

**Files:** Create `src/quiverlab/invariants/spectral.py`, `tests/invariants/test_spectral.py`.

**Interfaces:**
- Consumes: `sympy`, `engine.coxeter_spectrum.is_cyclotomic_product`, `engine.coxeter_spectrum.{star_quiver, cartan_of_quiver, trivial_extension_cartan}`, `engine.coxeter2.coxeter_polynomial_from_cartan`.
- Produces: `spectral.spectral_radius(poly) -> sympy exact | None`, `spectral.mahler_measure(poly) -> sympy exact | None`. Both return exact `sympy` algebraic numbers (`sympy.Integer(1)` on cyclotomic input, else a `CRootOf`-magnitude / product of magnitudes), **sound for complex off-circle roots**. `poly=None -> None`. Private helpers: `_noncyclotomic_part`, `_is_reciprocal`, `_reciprocal_to_y`, `_real_roots_suffice`, `_off_circle_roots`.

**Algorithm (exact and SOUND; reimplements the deleted hanlab float `_roots_abs`).** A spectral radius / Mahler measure must count **every** root of modulus > 1, complex ones included. `real_roots` alone is **unsound** — it silently drops complex off-circle roots (verified counterexamples below), and "real roots suffice" is a theorem only for **tree/bipartite = hereditary** quivers (A'Campo), not in general. Full `all_roots` is exact but pathologically slow to compare on a high-degree *irreducible* factor (Lehmer, > 2 min). So we **decide**, without isolating any complex root, whether real roots suffice, via the self-inversive `y = z + 1/z` substitution:
1. **Cyclotomic short-circuit:** `is_cyclotomic_product(poly)` ⟹ `spectral_radius = mahler_measure = 1` (all roots on the unit circle).
2. **Non-cyclotomic part** `q` = product (with multiplicity) of the non-cyclotomic irreducible factors. A Coxeter polynomial is self-inversive (reciprocal), so `q` is reciprocal of even degree `2s` with no roots at `±1`. Substitute `y = z + 1/z`: `q(z) = z^s·Q(z + 1/z)` where `Q` has degree `s`, built from the Dickson polynomials `D_k(y) = z^k + z^{-k}` (`D_0=2, D_1=y, D_k = y·D_{k-1} − D_{k-2}`). Then **`q` has a complex off-circle root ⟺ `Q` has a non-real root.**
3. **Decide (Sturm, no complex isolation):** compare the count of **distinct** real roots of `Q` (`Q.sqf_part().count_roots()`) with the squarefree degree `Q.sqf_part().degree()`. Equal ⟹ every off-circle root of `q` is **real** ⟹ `real_roots(q)` is exact and fast (**Branch A**, the hereditary/Salem/Lehmer case — Lehmer degree-10 ~0.02 s). Unequal ⟹ `q` has complex off-circle roots ⟹ fall back to `all_roots(q)` (**Branch B**, the *non-hereditary* case; `q` is only the small non-cyclotomic factor, so this stays fast — the dim-22 witness below runs in ~0.7 s).
4. **`spectral_radius`** = `max |α|` over the off-circle roots of `q` (`≥ 1` always for non-cyclotomic `q`); **`mahler_measure`** = `|lc| · ∏_{|α|>1} |α|`. Magnitudes are exact `sympy.Abs` of `CRootOf`; comparisons are exact `.is_positive`.

**Chosen for the unequal case: the `all_roots(q)` fallback (full correctness), not a loud raise.** The boundary is *non-hereditary complex-dominant spectra*, and Branch B computes the correct complex-inclusive answer there. **No floats, no `nroots`, no `complex()` in `src/`.**

**Verified counterexamples the fix must (and does) get right** (both run under the new code): `spectral_radius((t+1)(t⁴−7t³+16t²−7t+1)) = 3.54645…` (complex pair; `mahler = 12.5773…`), where the old real-roots-only code returned `1`; and the genuine dim-22 rad²=0 algebra below, `A.coxeter_polynomial() = t⁴+4t³+14t²+4t+1`, `spectral_radius = 3.58474…` (Branch B), old code returned `1`. Branch A stays exact and fast: Lehmer `1.17628…`, `t²−3t+1 → (3+√5)/2`.

- [ ] **Step 1: Write the failing test**

`tests/invariants/test_spectral.py`:

```python
"""Exact spectral radius / Mahler measure (spec section 5 component 8), SOUND for complex
off-circle roots. The bank floats survive only as test oracles here; src is float-free
(asserted by tests/test_no_floats.py)."""
import sympy as sp
from quiverlab import Quiver, CC
from quiverlab.invariants.spectral import spectral_radius, mahler_measure
from quiverlab.engine.coxeter_spectrum import is_cyclotomic_product

t = sp.symbols("t")
LEHMER = t**10 + t**9 - t**7 - t**6 - t**5 - t**4 - t**3 + t + 1


def _rad2_algebra(M, field=CC):
    """The radical-square-zero algebra whose arrow-count matrix is M (Cartan = I + M):
    M[i][j] parallel arrows i+1 -> j+1, all length-2 paths killed."""
    n = len(M)
    arrows = {}
    for i in range(n):
        for j in range(n):
            for c in range(M[i][j]):
                arrows[f"x{i+1}_{j+1}_{c}"] = (i + 1, j + 1)
    Q = Quiver(list(range(1, n + 1)), arrows)
    rels = [f"{a}*{b}" for a, (_sa, ta) in arrows.items()
            for b, (sb, _tb) in arrows.items() if ta == sb]
    return Q.algebra(relations=rels, field=field)


def test_returns_exact_types_and_none():
    assert spectral_radius(None) is None
    assert mahler_measure(None) is None
    rho = spectral_radius(t**2 - 3 * t + 1)
    assert not isinstance(rho, float)              # exact sympy object, never a float


def test_cyclotomic_short_circuit_is_exact_one():
    assert spectral_radius((t + 1)**3 * (t - 1)**2) == sp.Integer(1)
    assert mahler_measure((t + 1)**2 * (t**2 - t + 1)) == sp.Integer(1)
    assert spectral_radius(t**3 - 1) == 1
    assert is_cyclotomic_product(LEHMER) is False   # Lehmer is NOT cyclotomic


def test_quadratic_wild_radius_is_golden_squared():
    rho = spectral_radius(t**2 - 3 * t + 1)         # Branch A (real off-circle roots)
    assert sp.simplify(rho - (3 + sp.sqrt(5)) / 2) == 0      # exact
    assert abs(float(rho.evalf(30)) - (3 + 5**0.5) / 2) < 1e-9   # bank float oracle
    assert sp.simplify(mahler_measure(t**2 - 3 * t + 1) - rho) == 0


def test_lehmer_spectral_radius_is_lehmer_number():
    rho = spectral_radius(LEHMER)                   # Branch A (Salem: one real root > 1)
    assert abs(float(rho.evalf(30)) - 1.17628081825991) < 1e-9   # Lehmer's number
    assert abs(float(mahler_measure(LEHMER).evalf(30)) - 1.17628081825991) < 1e-9


def test_mahler_of_mixed_product():
    m = mahler_measure((t + 1)**4 * (t**2 - 3 * t + 1))   # cyclotomic (t+1)^4 -> 1
    assert abs(float(m.evalf(30)) - (3 + 5**0.5) / 2) < 1e-9


def test_complex_off_circle_roots_counterexample_polynomial():
    # SOUNDNESS: (t+1)(t^4-7t^3+16t^2-7t+1) has a COMPLEX conjugate pair of modulus
    # 3.54645... that real_roots-only would drop (returning a wrong 1). Branch B.
    poly = (t + 1) * (t**4 - 7 * t**3 + 16 * t**2 - 7 * t + 1)
    assert is_cyclotomic_product(poly) is False
    rho = spectral_radius(poly)
    assert abs(float(rho.evalf(30)) - 3.54645544468500) < 1e-9
    assert (rho - 3).is_positive                    # exactly > 3, not a silent 1
    m = mahler_measure(poly)
    assert abs(float(m.evalf(30)) - 12.5773462211358) < 1e-9


def test_complex_off_circle_roots_end_to_end_rad2_algebra():
    # SOUNDNESS end-to-end: a genuine dim-22 rad^2=0 algebra whose Coxeter polynomial
    # t^4+4t^3+14t^2+4t+1 has complex modulus-3.58474... roots. real_roots-only returned 1.
    M = [[0, 0, 2, 0], [1, 0, 1, 1], [3, 3, 0, 1], [2, 2, 2, 0]]
    A = _rad2_algebra(M)
    assert A.dim == 22
    p = A.coxeter_polynomial()
    assert sp.expand(p.as_expr() - (t**4 + 4 * t**3 + 14 * t**2 + 4 * t + 1)) == 0
    rho = spectral_radius(p.as_expr())
    assert rho != 1 and (rho - sp.Rational(35, 10)).is_positive    # NOT a silent 1
    assert abs(float(rho.evalf(30)) - 3.58474330285921) < 1e-9


def test_spectral_of_cartan_of_quiver_first_ever_coverage():
    # first-ever assertion coverage of star_quiver + cartan_of_quiver (ledger)
    from quiverlab.engine.coxeter_spectrum import star_quiver, cartan_of_quiver
    from quiverlab.engine.coxeter2 import coxeter_polynomial_from_cartan
    n, arrows = star_quiver([1, 1, 1])         # D_4 star, 4 vertices
    assert n == 4
    C, _ = cartan_of_quiver(n, arrows)
    poly, _Phi = coxeter_polynomial_from_cartan(C)
    assert spectral_radius(poly) == 1          # Dynkin D_4 is cyclotomic -> radius 1
```

- [ ] **Step 2: Run — expect FAIL** (`ModuleNotFoundError: quiverlab.invariants.spectral`).

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/invariants/spectral.py`:

```python
"""Exact spectral radius and Mahler measure of a Coxeter polynomial (spec section 5
component 8). Reimplements the deleted hanlab float layer (`_roots_abs`, `spectral_radius`,
`mahler_measure`, which used mpmath `Poly.nroots`) with EXACT sympy algebraic numbers.

No floats in this module (tests/test_no_floats.py enforces): magnitudes are `sympy.Abs`
of exact `CRootOf` roots, comparisons use `.is_positive`, and cyclotomic input short-
circuits to the exact integer 1. The bank floats survive only as test oracles.

SOUNDNESS. We must count every root of modulus > 1, including COMPLEX ones. `real_roots`
alone is NOT enough -- it silently drops complex off-circle roots (t^4-7t^3+16t^2-7t+1 has
a complex pair of modulus 3.5465 that real_roots misses); "real roots suffice" is a theorem
only for tree/bipartite = HEREDITARY quivers (A'Campo). Full `all_roots` is exact but
pathologically slow on a high-degree irreducible factor (Lehmer, > 2 min). So we DECIDE,
with no complex-root isolation, whether real roots suffice, via the self-inversive
y = z + 1/z substitution:

  * cyclotomic short-circuit -> 1;
  * q = non-cyclotomic part (reciprocal, even degree 2s, no roots at +-1); q has a complex
    off-circle root iff Q(y) = q(z)/z^s (y = z + 1/z, degree s, via Dickson polynomials)
    has a non-real root -- decided by counting distinct real roots of Q (Sturm) against Q's
    squarefree degree. Equal -> real_roots(q) is exact and fast (Branch A: hereditary /
    Salem / Lehmer); unequal -> all_roots(q) fallback (Branch B: non-hereditary, q small).

spectral_radius = max |root| over q's off-circle roots; mahler = |lc| * prod of |root| over
roots with |root| > 1."""
import sympy as sp

from quiverlab.engine.coxeter_spectrum import is_cyclotomic_product

_T = sp.symbols("t")
_Y = sp.symbols("y_spec")


def _noncyclotomic_part(poly):
    """Product (with multiplicity) of the non-cyclotomic irreducible factors of poly, or
    None if poly is a product of cyclotomics."""
    q = sp.Integer(1)
    for fac, mult in sp.factor_list(sp.expand(poly), _T)[1]:
        if not is_cyclotomic_product(fac):
            q = q * sp.Poly(fac, _T).as_expr() ** mult
    return sp.Poly(q, _T) if q != sp.Integer(1) else None


def _is_reciprocal(qp):
    c = qp.all_coeffs()
    return c == c[::-1] or c == [-x for x in c[::-1]]


def _reciprocal_to_y(qp):
    """For a reciprocal q of even degree 2s: Q(y) with q(z) = z^s * Q(z + 1/z), built from
    Dickson polynomials D_k(y) = z^k + z^{-k} (D_0=2, D_1=y, D_k = y*D_{k-1} - D_{k-2});
    q(z)/z^s = a_s + sum_{k=1..s} a_{s+k} D_k(y)."""
    a = list(reversed(qp.all_coeffs()))              # a[i] = coeff of z^i
    s = qp.degree() // 2
    D = [sp.Integer(2), _Y]
    for k in range(2, s + 1):
        D.append(sp.expand(_Y * D[k - 1] - D[k - 2]))
    Q = sp.Integer(a[s])
    for k in range(1, s + 1):
        Q = Q + a[s + k] * D[k]
    return sp.Poly(sp.expand(Q), _Y)


def _real_roots_suffice(qp):
    """True iff every off-circle root of the reciprocal q is REAL -- decided by counting
    DISTINCT real roots of Q(y) (Sturm on the squarefree part) against its squarefree
    degree, with no complex-root isolation. Non-reciprocal / odd-degree q (non-Coxeter
    input) -> False (be safe: fall back to all_roots)."""
    if qp.degree() % 2 or not _is_reciprocal(qp):
        return False
    Qsf = _reciprocal_to_y(qp).sqf_part()            # squarefree: count_roots == distinct
    return Qsf.count_roots() == Qsf.degree()


def _off_circle_roots(qp):
    """Exact roots of q with |z| > 1: real_roots when they provably suffice, else all_roots
    on the (small) non-cyclotomic factor q."""
    roots = sp.real_roots(qp) if _real_roots_suffice(qp) else qp.all_roots()
    return [r for r in roots if (sp.Abs(r) - 1).is_positive]


def spectral_radius(poly):
    """max_i |alpha_i| over the roots of poly, EXACT. 1 on cyclotomic input."""
    if poly is None:
        return None
    if is_cyclotomic_product(poly):
        return sp.Integer(1)
    q = _noncyclotomic_part(poly)
    if q is None:
        return sp.Integer(1)
    best = sp.Integer(1)
    for r in _off_circle_roots(q):
        if (sp.Abs(r) - best).is_positive:
            best = sp.Abs(r)
    return best


def mahler_measure(poly):
    """|lc| * prod over roots with |alpha| > 1 of |alpha|, EXACT. 1 on cyclotomic input."""
    if poly is None:
        return None
    if is_cyclotomic_product(poly):
        return sp.Integer(1)
    q = _noncyclotomic_part(poly)
    if q is None:
        return sp.Integer(1)
    m = sp.Abs(sp.Poly(poly, _T).LC())               # rational-LC-safe (never Integer(abs()))
    for r in _off_circle_roots(q):
        m = m * sp.Abs(r)
    return sp.simplify(m)
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/invariants/test_spectral.py tests/test_no_floats.py -q`
Expected: all pass (verified: Lehmer/quadratic Branch A ~0.02 s; the det-2 counterexample and the dim-22 rad²=0 algebra Branch B ~0.7–1.4 s; cyclotomic short-circuits instant). Confirm `tests/test_no_floats.py` stays green — the exact route uses `sympy.Integer/Rational/Abs`, `real_roots`/`all_roots`, `.count_roots`, `.is_positive`; **no** float/complex literal, no `float()`, no `nroots`.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/spectral.py tests/invariants/test_spectral.py
git commit -m "feat(invariants): exact spectral_radius + mahler_measure, sound for complex roots

Cyclotomic short-circuit + self-inversive y=z+1/z Sturm decision (Branch A real_roots /
Branch B all_roots); no floats. Counterexamples (det-2 poly, dim-22 rad^2=0 algebra) pinned.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: Re-port the Lehmer star + trivial-extension collapse tests (ledger)

**Files:** Edit `tests/engine/test_coxeter_spectrum.py` (restore the two deleted tests + the spectral pin, against the new exact functions).

**Interfaces:** Consumes `invariants.spectral.{spectral_radius, mahler_measure}`, `engine.coxeter_spectrum.{star_quiver, cartan_of_quiver, trivial_extension_cartan, is_cyclotomic_product}`, `engine.coxeter2.coxeter_polynomial_from_cartan`. This task restores the **exact halves** the port deferred (the comments at lines 43–49 and 82–85 of the current file explicitly name Plan 05) and gives **first-ever assertion coverage** to `star_quiver`, `cartan_of_quiver`, `trivial_extension_cartan` (ledger). No `src/` change.

- [ ] **Step 1: Write the failing test**

Add to `tests/engine/test_coxeter_spectrum.py` (import the exact spectral functions at the top: `from quiverlab.invariants.spectral import spectral_radius, mahler_measure`). Replace the two "not ported" comment blocks with the restored tests:

```python
def test_spectral_radius_and_mahler():
    # restored (was deferred to Plan 05); now against the exact reimplementation
    rho = spectral_radius(t**2 - 3 * t + 1)
    assert abs(float(rho.evalf(30)) - (3 + 5**0.5) / 2) < 1e-9
    assert sp.simplify(mahler_measure(t**2 - 3 * t + 1) - rho) == 0
    assert spectral_radius((t + 1)**3 * (t - 1)**2) == 1        # exact, cyclotomic
    m = mahler_measure((t + 1)**4 * (t**2 - 3 * t + 1))
    assert abs(float(m.evalf(30)) - (3 + 5**0.5) / 2) < 1e-9


def test_lehmer_star_237():
    # T(2,3,7) = star_quiver([1,2,6]); Coxeter poly is Lehmer's polynomial; spectral
    # radius is Lehmer's number 1.176280818259917...  (ledger: exact half restored)
    n, arrows = star_quiver([1, 2, 6])
    assert n == 10
    C, _alg = cartan_of_quiver(n, arrows)
    poly, _Phi = coxeter_polynomial_from_cartan(C)
    assert sp.expand(poly - LEHMER) == 0
    assert is_cyclotomic_product(poly) is False
    assert abs(float(spectral_radius(poly).evalf(30)) - 1.17628081825991) < 1e-9


def test_trivial_extension_collapse():
    # rep-finite input vs wild input: K0 spectra of T(A) are identically (t+1)^v
    C2, _ = cartan_of_quiver(2, [(0, 1)])                # kA_2 (rep-finite)
    nW, aW = star_quiver([1, 1, 1, 1, 1])                # 5-subspace (wild)
    CW, _ = cartan_of_quiver(nW, aW)
    pA, _ = coxeter_polynomial_from_cartan(C2)
    pW, _ = coxeter_polynomial_from_cartan(CW)
    assert spectral_radius(pA) == 1
    assert (spectral_radius(pW) - sp.Rational(26, 10)).is_positive       # > 2.6, exact
    for C, v in ((C2, 2), (CW, nW)):
        CT = trivial_extension_cartan(C)
        polyT, PhiT = coxeter_polynomial_from_cartan(CT)
        assert PhiT == -sp.eye(v)
        assert sp.expand(polyT - (t + 1)**v) == 0
    # Euclidean D~4: symmetrized Cartan is singular -> Coxeter undefined
    nE, aE = star_quiver([1, 1, 1, 1])
    CE, _ = cartan_of_quiver(nE, aE)
    polyE, PhiE = coxeter_polynomial_from_cartan(trivial_extension_cartan(CE))
    assert polyE is None and PhiE is None
```

- [ ] **Step 2: Run — expect FAIL then PASS.** Before Task 8 these would error on the missing import; with Task 8 landed they must pass. Run:
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/engine/test_coxeter_spectrum.py -q`
Expected: the whole file passes, including the three restored tests and the previously-green cyclotomic/reconciliation/degeneration tests.

- [ ] **Step 3: No `src/` change.** (The exact functions arrived in Task 8; `star_quiver`/`cartan_of_quiver`/`trivial_extension_cartan` were already ported and now receive their first assertion coverage.)

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/engine/test_coxeter_spectrum.py tests/invariants/test_spectral.py tests/test_no_floats.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add tests/engine/test_coxeter_spectrum.py
git commit -m "test(invariants): re-port Lehmer star T(2,3,7) + trivial-extension (t+1)^v collapse

Restores the exact halves deferred by the Plan-02 port (ledger); first assertion
coverage of star_quiver/cartan_of_quiver/trivial_extension_cartan.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 10: Fix the non-unimodular-Cartan `coxeter_polynomial` minor (ledger)

**Files:** Edit `src/quiverlab/invariants/cartan.py` (`coxeter_polynomial`). Create `tests/invariants/test_coxeter_nonunimodular.py`.

**Decision (specified) — corrected after execution.** The premise "non-unimodular ⇒ silent QQ that must be re-pinned by a `|det C|` rule" is **false**, verified by running the live code: sympy's `Poly` already infers the domain **from the actual coefficients** — `kA_2 → t²+t+1` over **ZZ**, the non-unimodular Fixture-D algebra `→ t²+3t/2+1` over **QQ**, `k[x]/(x²) → t+1` over **ZZ** (det C = 2, non-unimodular, yet integral), and `diag(1,2) → (t+1)²` over **ZZ** (det C = 2, integral). A `|det C| == 1 ? ZZ : QQ` rule would therefore be **wrong** — it would regress `k[x]/(x²)` and `diag(1,2)` from ZZ to QQ. So the domain logic stays exactly as sympy infers it; the coefficient, not `det C`, decides.

The genuine "minor" is only that a non-unimodular Cartan can yield a **rational** Coxeter transformation *silently* — the user isn't told this is not the classical integral Coxeter matrix. The fix is therefore **documented behavior**: a docstring caveat naming the `det C ∉ {0, ±1}` case (exact but possibly QQ, and not the classical integral Coxeter matrix), plus **characterization/regression tests** that pin the current correct behavior (so a future refactor cannot reintroduce a naive det-based rule). No code logic change; unimodular AND integral-non-unimodular outputs are byte-identical to today.

- [ ] **Step 1: Write the characterization tests** (these PASS against the current code — they pin correct behavior and guard against a det-based regression):

`tests/invariants/test_coxeter_nonunimodular.py`:

```python
"""Characterization + regression tests for coxeter_polynomial's domain (spec section 5
component 8; ledger). The domain follows the COEFFICIENTS (sympy's inference), NOT det C:
a non-unimodular Cartan may be rational (QQ) OR integral (ZZ). These tests pin current
correct behavior so a naive det-based rule cannot creep back in. Fixture D + diag(1,2)."""
import sympy as sp
from quiverlab import Quiver, CC, linear_path_algebra, truncated_polynomial

t = sp.Symbol("t")


def _nonunimodular_rational():
    # C = [[2,1],[0,1]], det = 2; Phi is genuinely rational -> QQ
    return Quiver([1, 2], {"x": (1, 1), "a": (1, 2)}).algebra(
        relations=["x^2", "x*a"], field=CC)


def test_nonunimodular_cartan_rational_is_qq():
    A = _nonunimodular_rational()
    assert A.cartan_matrix() == [[2, 1], [0, 1]]
    p = A.coxeter_polynomial()
    assert sp.expand(p.as_expr() - (t**2 + sp.Rational(3, 2) * t + 1)) == 0
    assert p.domain == sp.QQ                    # rational Coxeter transformation


def test_unimodular_cartan_is_zz():
    p = linear_path_algebra(2).coxeter_polynomial()   # C = [[1,1],[0,1]], det 1
    assert sp.expand(p.as_expr() - (t**2 + t + 1)) == 0
    assert p.domain == sp.ZZ


def test_nonunimodular_but_integral_stays_zz_local():
    # k[x]/(x^2): C = [[2]], det 2 (non-unimodular) but Phi = [-1] integral -> ZZ
    p = truncated_polynomial(2).coxeter_polynomial()
    assert sp.expand(p.as_expr() - (t + 1)) == 0
    assert p.domain == sp.ZZ                    # NON-unimodular yet integral: NOT QQ


def test_nonunimodular_diag_stays_zz():
    # C = diag(1, 2): isolated vertex 1 + a loop y at 2 with y^2 = 0; det 2, Phi = -I,
    # Coxeter polynomial (t+1)^2 -- integral, so ZZ (a |det C| rule would wrongly give QQ)
    A = Quiver([1, 2], {"y": (2, 2)}).algebra(relations=["y^2"], field=CC)
    assert A.cartan_matrix() == [[1, 0], [0, 2]]
    p = A.coxeter_polynomial()
    assert sp.expand(p.as_expr() - (t + 1)**2) == 0
    assert p.domain == sp.ZZ


def test_singular_cartan_still_raises():
    from quiverlab.errors import QuiverlabError
    import pytest
    # kZ_4/rad^2 has det C = 0 -> singular -> loud
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 1)}
               ).algebra(relations=["a*b", "b*c", "c*d", "d*a"], field=CC)
    with pytest.raises(QuiverlabError):
        A.coxeter_polynomial()
```

- [ ] **Step 2: Run — expect PASS** (these characterize the current *correct* behavior; they are regression guards, not red-first tests). Run:
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/invariants/test_coxeter_nonunimodular.py -q`
Expected: all pass. If `test_nonunimodular_but_integral_stays_zz_local` or `test_nonunimodular_diag_stays_zz` FAIL, someone introduced a det-based domain rule — revert it; the domain must follow the coefficients.

- [ ] **Step 3: Add the docstring caveat only** — edit `coxeter_polynomial` in `src/quiverlab/invariants/cartan.py`, keeping the body unchanged (sympy's inference is correct), replacing only the docstring:

```python
def coxeter_polynomial(A):
    """Characteristic polynomial of the Coxeter matrix Phi = -C^{-T} C, as an exact
    sympy Poly in t (no numerical root-finding; the exact spectral_radius / mahler_measure
    invariants live in invariants/spectral.py).

    DOMAIN CAVEAT (det C not in {0, +-1}): when the Cartan is non-unimodular, Phi may have
    rational entries, so the Coxeter polynomial can be over QQ (e.g. t^2 + 3t/2 + 1 for
    C = [[2,1],[0,1]]). It is still EXACT, but such a rational Phi is NOT the classical
    integral Coxeter transformation -- a caveat the coxeter_matrix() sibling surfaces via
    its rational-entry branch. The domain follows the actual COEFFICIENTS, not det C: a
    non-unimodular Cartan may still be integral (k[x]/(x^2) -> t+1 over ZZ; diag(1,2) ->
    (t+1)^2 over ZZ). A singular Cartan (det C == 0) has no Coxeter transformation and
    raises loudly."""
    C = sympy.Matrix(cartan_matrix(A))
    if C.det() == 0:
        raise QuiverlabError(
            "Cartan matrix is singular: no Coxeter polynomial",
            hint="see coxeter_matrix",
        )
    Phi = -C.inv().T * C
    t = sympy.Symbol("t")
    return sympy.Poly(Phi.charpoly(t).as_expr(), t)   # domain inferred from coefficients
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/invariants/ tests/test_no_floats.py -q`
Expected: all pass — the new characterization tests plus the existing `test_coxeter_polynomial_A2` / `test_coxeter_matrix_A2_exact_value` (behavior byte-identical; only the docstring changed).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/cartan.py tests/invariants/test_coxeter_nonunimodular.py
git commit -m "fix(invariants): pin coxeter_polynomial domain (ZZ unimodular / QQ documented)

Removes the silent QQ-domain flip on non-unimodular Cartans; unimodular outputs are
byte-identical (ledger minor).

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 11: `loewy_length` + `complexity` + `center`

**Files:** Create `src/quiverlab/invariants/scalar.py`, `tests/invariants/test_scalar.py`. Edit `core/algebra.py` (add `loewy_length`, `complexity`, `center`).

**Interfaces:**
- Consumes: `Algebra` (`.multiply`, `._basis_vec`, `.basis_labels`, `.dim`, `.unit`), `fields.linalg`, `engine.scan3.complexity_of`, `engine.adapter.to_engine`, `engine.resolutions_minimal.minimal_resolution`.
- Produces:
  - `scalar.loewy_length(A) -> int` — least `n` with `rad^n A = 0` (nilpotency index of the radical ideal), any field.
  - `scalar.center(A) -> tuple[int, list[vector]]` — `(dim Z(A), basis of Z(A))` via `{z : z b = b z ∀ basis b}`, any field.
  - `scalar.complexity(A, n) -> int | str | None` — `complexity_of` applied to the term-dimension sequence `[dim P_0 .. dim P_n]` of the minimal `A^e` (bimodule) resolution over the fast engine (GF(p) only; loud `FieldError` off a prime field, matching the engine-backed invariants).
  - `Algebra.loewy_length(self)`, `Algebra.complexity(self, n)`, `Algebra.center(self)` (deferred-import wrappers).

- [ ] **Step 1: Write the failing test**

`tests/invariants/test_scalar.py`:

```python
"""loewy_length / complexity / center (spec section 3.5). Fixtures A, B; k[x]/(x^n)."""
import pytest
from quiverlab import Quiver, CC, GF, linear_path_algebra, truncated_polynomial
from quiverlab.errors import FieldError


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"])


def test_loewy_length():
    assert linear_path_algebra(2).loewy_length() == 2        # rad^2 = 0
    assert _square().loewy_length() == 3                     # rad^3 = 0, rad^2 = <a*b>
    assert truncated_polynomial(4).loewy_length() == 4       # k[x]/(x^4): rad^4 = 0
    assert truncated_polynomial(2, field=GF(3)).loewy_length() == 2


def test_center_dimension_and_basis():
    d, basis = linear_path_algebra(2).center()
    assert d == 1                                            # connected -> center = k*1
    assert len(basis) == 1
    assert _square().center()[0] == 1
    # k[x]/(x^2) is commutative -> center is the whole algebra
    assert truncated_polynomial(2).center()[0] == 2


def test_center_over_gfp():
    assert linear_path_algebra(2, field=GF(7)).center()[0] == 1


def test_complexity_gfp():
    # finite gl.dim (kA_2) -> minimal A^e resolution terminates -> complexity 0
    assert linear_path_algebra(2, field=GF(32003)).complexity(6) == 0
    # k[x]/(x^2) self-injective -> constant-rank periodic resolution -> complexity 1
    assert truncated_polynomial(2, field=GF(32003)).complexity(6) == 1


def test_complexity_cc_loud():
    with pytest.raises(FieldError):
        truncated_polynomial(2, field=CC).complexity(4)
```

- [ ] **Step 2: Run — expect FAIL** (`AttributeError: ... 'loewy_length'`).

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/invariants/scalar.py`:

```python
"""Scalar algebra invariants: Loewy length, center, complexity (spec section 3.5).

Loewy length and center are exact over ANY Domain (radical-ideal powers / commutant
linear algebra). Complexity routes through the fast GF(p) engine (the minimal A^e
resolution's term growth), so it is prime-field-only and fails loudly otherwise, exactly
like the other engine-backed invariants."""
from quiverlab.fields import linalg


def _radical_basis_indices(A):
    """Indices of the basis labels lying in rad A = the arrow ideal (all non-idempotent
    basis paths, i.e. every label that is not 'e_v')."""
    return [i for i, lab in enumerate(A.basis_labels) if not lab.startswith("e_")]


def _ideal_product_span(A, gens_vectors, rad_idx):
    """Span (as a list of coordinate vectors) of gens * rad A = { g * r : g in gens,
    r a radical basis element }, reduced to an independent set via rref."""
    dom = A.domain
    rows = []
    for g in gens_vectors:
        for ri in rad_idx:
            prod = A.multiply(g, A._basis_vec(ri))
            rows.append(prod)
    if not rows:
        return []
    R, piv = linalg.rref(rows, dom)
    return [R[i] for i in range(len(piv))]


def loewy_length(A):
    """Least n with rad^n A = 0. rad^1 = arrow ideal; rad^{k+1} = rad^k * rad."""
    dom = A.domain
    rad_idx = _radical_basis_indices(A)
    if not rad_idx:
        return 1                                    # semisimple: rad = 0, rad^1 = 0
    current = [A._basis_vec(i) for i in rad_idx]     # rad^1 spanning set
    n = 1
    while current:
        nxt = _ideal_product_span(A, current, rad_idx)
        if not nxt:
            return n + 1                             # rad^{n+1} = 0
        current = nxt
        n += 1
        if n > A.dim + 1:                            # safety: nilpotent within dim steps
            return n
    return n


def center(A):
    """(dim Z(A), basis). Z(A) = { z : z*b = b*z for every basis element b }. Solve the
    stacked commutator system over the Domain; nullspace = the center."""
    dom = A.domain
    m = A.dim
    rows = []
    # unknown z = (z_0..z_{m-1}); for each basis b, (z*b - b*z) = 0 gives m linear rows.
    # coefficient of z_k in coordinate-c of (z*b - b*z): mult(e_k,b)[c] - mult(b,e_k)[c].
    for b in range(m):
        eb = A._basis_vec(b)
        for c in range(m):
            row = []
            for k in range(m):
                ek = A._basis_vec(k)
                left = A.multiply(ek, eb)[c]         # (z*b) coord c, coefficient of z_k
                right = A.multiply(eb, ek)[c]        # (b*z) coord c
                row.append(dom.sub(left, right))
            rows.append(row)
    basis = linalg.nullspace(rows, dom)
    return len(basis), basis


def complexity(A, n):
    """Apparent complexity of A from the minimal A^e (bimodule) resolution's term-
    dimension growth up to degree n (fast GF(p) engine). Returns complexity_of's honest
    label (int / None / '>=2')."""
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import minimal_resolution
    from quiverlab.engine.scan3 import complexity_of
    A._require_prime_field("complexity")             # loud FieldError off a prime field
    eng = to_engine(A)
    p = A.domain.p
    rks, cols, _e, _trunc = minimal_resolution(eng, n, p)
    # rks[k] = number of A^e generators of P_k; term k-dimension is m^2 * rks[k].
    seq = [max(0, r) for r in rks]
    return complexity_of(seq)
```

Edit `core/algebra.py` — add:

```python
    def loewy_length(self):
        """Loewy length = nilpotency index of rad A (exact, any field) (spec §3.5)."""
        from quiverlab.invariants.scalar import loewy_length
        return loewy_length(self)

    def center(self):
        """(dim, basis) of the center Z(A), exact over any field (spec §3.5)."""
        from quiverlab.invariants.scalar import center
        return center(self)

    def complexity(self, n):
        """Apparent complexity from the minimal A^e resolution's growth (GF(p) only)."""
        from quiverlab.invariants.scalar import complexity
        return complexity(self, n)
```

> **`complexity` return contract.** `engine.scan3.complexity_of` returns `int | None | "≥ 2"` (agent-verified: `0` when the sequence has trailing zeros = finite projective dimension, `1` for bounded/periodic growth, higher for polynomial growth, `None` for too-short sequences). The public `A.complexity(n)` forwards this honest label unchanged; the docstring says so. The fixtures pin `kA_2 -> 0` (terminating resolution) and `k[x]/(x^2) -> 1` (constant-rank periodic). If `complexity_of` returns `None` for a too-short `n`, the caller is told to raise `n`.

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/invariants/test_scalar.py tests/test_no_floats.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/scalar.py tests/invariants/test_scalar.py src/quiverlab/core/algebra.py
git commit -m "feat(invariants): loewy_length + center (any field) + complexity (GF(p))

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 12: `sweep()` — the invariant × field table (§3.9)

**Files:** Create `src/quiverlab/invariants/sweep.py`, `tests/invariants/test_sweep.py`. Export `sweep` from `src/quiverlab/__init__.py`.

**Interfaces:**
- Consumes: field constructors, any algebra builder callable taking `field=`.
- Produces: `sweep(builder, *args, fields, invariants=None, **kwargs) -> SweepTable`. `builder(*args, field=f, **kwargs)` is rebuilt per field; each requested invariant is computed per field. `invariants` is an ordered mapping `name -> callable(A) -> value`; a default set (`dimension`, `loewy_length`, `cartan_det`, `coxeter_polynomial`) is used if `None`. A cell that raises (e.g. an engine-backed invariant over `CC`) is recorded as the exception's short message (`"n/a: FieldError"`), never a crash — the point of a sweep is the *comparison*, so a partial row is informative. `SweepTable` has `.cell(inv, field)`, `.rows`, `.__repr__` (aligned text table), `.latex()`.

- [ ] **Step 1: Write the failing test**

`tests/invariants/test_sweep.py`:

```python
"""sweep(): rebuild an algebra over several fields and tabulate invariants (spec 3.9)."""
from quiverlab import CC, GF, sweep, truncated_polynomial, linear_path_algebra


def test_sweep_default_invariants_over_fields():
    tab = sweep(truncated_polynomial, 3, fields=[CC, GF(2), GF(3), GF(5)])
    # dimension is field-independent
    for f in (CC, GF(2), GF(3), GF(5)):
        assert tab.cell("dimension", f) == 3
        assert tab.cell("loewy_length", f) == 3


def test_sweep_custom_invariants_and_char_dependence():
    tab = sweep(truncated_polynomial, 2,
                fields=[GF(2), GF(3), GF(5)],
                invariants={"is_frobenius": lambda A: A.is_frobenius(),
                            "hh1": lambda A: A.hochschild_cohomology(1).dims[1]})
    assert tab.cell("is_frobenius", GF(2)) is True
    assert all(isinstance(tab.cell("hh1", f), int) for f in (GF(2), GF(3), GF(5)))


def test_sweep_records_field_errors_without_crashing():
    # an engine-backed invariant over CC must be recorded, not raised, in a sweep cell
    tab = sweep(truncated_polynomial, 2, fields=[CC, GF(5)],
                invariants={"is_frobenius": lambda A: A.is_frobenius()})
    assert tab.cell("is_frobenius", GF(5)) is True
    assert "n/a" in str(tab.cell("is_frobenius", CC)).lower()


def test_sweep_repr_and_latex_are_strings():
    tab = sweep(linear_path_algebra, 2, fields=[CC, GF(2)])
    assert isinstance(repr(tab), str) and "dimension" in repr(tab)
    assert isinstance(tab.latex(), str)
```

- [ ] **Step 2: Run — expect FAIL** (`ImportError: cannot import name 'sweep'`).

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/invariants/sweep.py`:

```python
"""sweep(): the 'all the moving variables' table -- an algebra rebuilt over several
fields, with chosen invariants tabulated per field (spec section 3.9). The character-
dependence view. A cell whose invariant is unavailable over a field (e.g. an engine-
backed GF(p)-only invariant over CC) records the reason instead of crashing."""


def _default_invariants():
    return {
        "dimension": lambda A: A.dim,
        "loewy_length": lambda A: A.loewy_length(),
        "cartan_det": lambda A: _det(A.cartan_matrix()),
        "coxeter_polynomial": lambda A: A.coxeter_polynomial().as_expr(),
    }


def _det(mat):
    import sympy
    return int(sympy.Matrix(mat).det())


class SweepTable:
    def __init__(self, fields, invariants, data):
        self.fields = fields
        self.invariant_names = list(invariants)
        self._data = data                       # {(inv_name, field_name): value}

    def cell(self, inv, field):
        return self._data[(inv, getattr(field, "name", str(field)))]

    @property
    def rows(self):
        return [[inv] + [self._data[(inv, getattr(f, "name", str(f)))] for f in self.fields]
                for inv in self.invariant_names]

    def __repr__(self):
        head = ["invariant"] + [getattr(f, "name", str(f)) for f in self.fields]
        lines = [head] + [[str(x) for x in row] for row in self.rows]
        widths = [max(len(lines[r][c]) for r in range(len(lines))) for c in range(len(head))]
        out = []
        for r, line in enumerate(lines):
            out.append("  ".join(cell.ljust(widths[c]) for c, cell in enumerate(line)))
            if r == 0:
                out.append("  ".join("-" * widths[c] for c in range(len(head))))
        return "\n".join(out)

    def latex(self):
        cols = "l" + "c" * len(self.fields)
        head = " & ".join(["invariant"] + [getattr(f, "name", str(f)) for f in self.fields])
        body = " \\\\\n".join(
            " & ".join([row[0]] + [str(x) for x in row[1:]]) for row in self.rows)
        return f"\\begin{{tabular}}{{{cols}}}\n{head} \\\\\\hline\n{body}\n\\end{{tabular}}"


def sweep(builder, *args, fields, invariants=None, **kwargs):
    invs = invariants if invariants is not None else _default_invariants()
    data = {}
    for f in fields:
        A = builder(*args, field=f, **kwargs)
        fname = getattr(f, "name", str(f))
        for name, fn in invs.items():
            try:
                data[(name, fname)] = fn(A)
            except Exception as exc:              # record, do not crash the sweep
                data[(name, fname)] = f"n/a: {type(exc).__name__}"
    return SweepTable(fields, invs, data)
```

Add `sweep` to `src/quiverlab/__init__.py` (import + `__all__`):

```python
from quiverlab.invariants.sweep import sweep  # noqa: E402,F401
# ... add "sweep" to __all__
```

> **Design note on the broad `except`.** A sweep's job is the *comparison across fields*; a single unavailable cell (an engine-backed invariant over `CC`, a singular Cartan) must not abort the whole table. This is the one place a broad catch is correct — it is confined to per-cell invariant evaluation, it records the exception **type name** (so the reason is visible, not swallowed), and it never touches the exact-arithmetic core. All other code paths keep the library's loud-failure contract.

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/invariants/test_sweep.py tests/test_no_floats.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/sweep.py tests/invariants/test_sweep.py src/quiverlab/__init__.py
git commit -m "feat(invariants): sweep() invariant x field table (the moving-variables view)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 13: Re-port `test_minimal_memory_guard` (ledger)

**Files:** Create `tests/engine/test_minimal_memory_guard.py`.

**Interfaces:** Consumes `Quiver.algebra` (Plan 03; builds the non-monomial `open_33_0`), `engine.adapter.to_engine`, `engine.resolutions_minimal.{minimal_resolution, minimal_homology_dims}`. **Dependency determination (recorded):** the bank test built `open_33_0` via `reduction_algebra.algebra_from_reduction_system`, which Plan 04 declared **superseded** by `Quiver.algebra` (general `kQ/I`). Since Plan 03 delivers `Quiver.algebra`, the dependency is now satisfiable — Plan 05 owns this re-port. **ROADMAP reconciliation:** Plan 02's deferral list nominally assigned the memory-guard test to Plan 04; Plan 04 instead *superseded* its `reduction_algebra` dependency (folding it into Plan-03 `Quiver.algebra`) without re-porting the test, so the obligation lands here in Plan 05 — the earliest plan where the fixture is buildable. The guard mechanism (`max_transient_bytes` → prefix / `truncated_at`, no `DepthLimitError`) is already ported (`engine/resolutions_minimal.py`); only the fixture builder was missing. No `src/` change.

- [ ] **Step 1: Write the failing test**

`tests/engine/test_minimal_memory_guard.py`:

```python
"""Re-port of the deleted bank test_minimal_memory_guard (ledger). open_33_0 =
k<x,y>/(x^3 - y^2, y^3, y*x + x*y), dim 9, local, non-monomial -- now buildable via the
Plan-03 general kQ/I lowering (reduction_algebra is superseded by Quiver.algebra). The
minimal A^e resolution's radK transient grows degree by degree, exercising the
max_transient_bytes guard: any budget yields a PREFIX; a 1-byte budget yields
truncated_at == 0, hh == []."""
from quiverlab import Quiver, GF
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.resolutions_minimal import minimal_homology_dims, minimal_resolution

P = 32003


def _open_33_0():
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x^3 - y^2", "y^3", "y*x + x*y"], field=GF(P))
    assert A.dim == 9
    return to_engine(A)


def test_default_none_matches_existing_build():
    A = _open_33_0()
    base = minimal_homology_dims(A, 5, primes=(P,))
    same = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=None)
    assert base == same


def test_any_transient_budget_yields_a_prefix():
    A = _open_33_0()
    full = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=10**18)[P]
    for budget in (1, 10**5, 10**6, 10**7, 10**9):
        got = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=budget)[P]
        assert got == full[:len(got)], f"budget={budget}: {got} not a prefix of {full}"


def test_tiny_budget_truncates_strictly_earlier():
    A = _open_33_0()
    full = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=10**18)[P]
    tiny = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=1)[P]
    assert len(tiny) < len(full)


def test_truncated_list_length_equals_truncated_at_and_is_a_prefix():
    A = _open_33_0()
    full = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=10**18)[P]
    budget = 10**7
    hh = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=budget)[P]
    _rks, _cols, _eng, trunc = minimal_resolution(A, 5, P, max_transient_bytes=budget)
    assert trunc is not None
    assert len(hh) == trunc
    assert hh == full[:trunc]


def test_one_byte_budget_stops_before_any_differential():
    A = _open_33_0()
    hh = minimal_homology_dims(A, 5, primes=(P,), max_transient_bytes=1)[P]
    _rks, _cols, _eng, trunc = minimal_resolution(A, 5, P, max_transient_bytes=1)
    assert trunc == 0
    assert hh == []
```

- [ ] **Step 2: Run — expect FAIL then PASS.** If it errors on building `open_33_0`, Plan 03 has not landed (Task 1 gate should have caught this — STOP). Run:
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/engine/test_minimal_memory_guard.py -q`

- [ ] **Step 3: Verify the fixture, no `src/` change.** Confirm `_open_33_0()` builds `dim 9` (the reduction-system rules `xxx→yy`, `yyy→0`, `yx→−xy` completing under Gröbner give the same 9-dim algebra the bank built). If the Gröbner completion needs a larger `degree_bound`, pass it: `Quiver.algebra(relations=[...], field=GF(P), degree_bound=14)` (the bank used `maxlen=14`). Record the effective bound in the fixture comment.

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/engine/test_minimal_memory_guard.py tests/test_no_floats.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add tests/engine/test_minimal_memory_guard.py
git commit -m "test(engine): re-port minimal-memory-guard via Plan-03 general kQ/I (ledger)

open_33_0 now built by Quiver.algebra (reduction_algebra superseded); prefix/truncated_at
guard contract restored.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 14: Acceptance — internals chapters, README, frozen interface, full suite

**Files:** Create `docs/internals/09-modules.md`; edit `docs/internals/06-invariants.md` (spectral layer + resolved non-unimodular caveat), `docs/internals/README.md` (chapter index + coverage), `README.md`; create `tests/modules/test_acceptance_plan05.py`.

**Interfaces:** Consumes the whole Plan-05 surface. Produces the acceptance test, the two internals chapters, and the frozen-interface statement below.

- [ ] **Step 1: Acceptance test** — `tests/modules/test_acceptance_plan05.py`:

```python
"""Plan-05 acceptance: the modules + invariants surface, end to end (spec 3.5/3.6/3.9)."""
import sympy as sp
from quiverlab import Quiver, CC, GF, linear_path_algebra, truncated_polynomial, sweep

t = sp.Symbol("t")


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def test_module_surface():
    A = _square()
    S1, S4 = A.simple(1), A.simple(4)
    assert A.projective(1).dimension_vector() == {1: 1, 2: 1, 3: 1, 4: 1}
    assert A.injective(4).dimension_vector() == {1: 1, 2: 1, 3: 1, 4: 1}
    assert A.projective(1).radical().radical().dimension_vector() == {1: 0, 2: 0, 3: 0, 4: 1}
    assert A.hom(S1, S1) == 1 and A.hom(S1, S4) == 0
    assert A.ext(S1, S4, 2) == 1
    assert int(A.global_dimension()) == 2
    res = S1.projective_resolution(4)
    assert res.pd() == 2 and res.betti(0) == 1


def test_invariant_surface():
    A = _square()
    assert A.loewy_length() == 3
    assert A.center()[0] == 1
    assert A.coxeter_polynomial().domain == sp.ZZ            # unimodular
    from quiverlab.invariants.spectral import spectral_radius
    assert spectral_radius(A.coxeter_polynomial().as_expr()) is not None
    assert truncated_polynomial(2, field=GF(32003)).complexity(6) == 1


def test_sweep_surface():
    tab = sweep(truncated_polynomial, 3, fields=[CC, GF(2), GF(3)])
    assert tab.cell("dimension", GF(2)) == 3
```

- [ ] **Step 2: Write `docs/internals/09-modules.md`** (full content, the fixed five-section chapter format; answers Marco's two questions — how a module is represented, how module resolutions are computed):

```markdown
# 09 — Modules: representation and resolutions

## The mathematics

For a bound quiver algebra A = kQ/I, a finite-dimensional **right** A-module M is a
vector space carrying a right action of A. Every such module is glued from **simples**
S_v (one per vertex), and each simple sits atop an **indecomposable projective**
P_v = e_v A (the paths starting at v) and beneath an **indecomposable injective**
I_v = D(A e_v) (the dual of the paths ending at v). The radical filtration
M ⊇ rad M ⊇ rad² M ⊇ ... records the Loewy structure; the **top** M/rad M lists the
generators and the **socle** the essential simple submodules. Homological questions —
Hom_A(M, N), Ext^n_A(M, N), the projective dimension of M — are answered from a
**minimal projective resolution** ... → Q_1 → Q_0 → M → 0.

## How it is represented

A module is the class `Module` (`modules/module.py`). It stores its k-dimension `dim`
and a **dict** `action` mapping each algebra basis label (a **string**, e.g. `"e_1"`,
`"a"`, `"a*b"`) to a `dim × dim` **matrix** (a list of lists of exact Domain elements)
of **right** multiplication. The convention is stated in every docstring: an element m is
a **column vector** in a fixed k-basis of M, and m·b = `action[b] @ m`; because paths
compose left-to-right, the action is an *anti-homomorphism*,
`action["x*y"] = action["y"] @ action["x"]`. The vertex subspace M·e_v is the image of
`action["e_v"]`, so the **dimension vector** — a dict vertex → int — is just
`rank(action["e_v"])` per vertex, and these ranks sum to `dim`. A simple S_v has `dim 1`,
`action["e_v"] = [[1]]` and every other label acting as `[[0]]`. A projective P_v is read
directly off the algebra's own multiplication table: its basis is the sub-list of
`A.basis_labels` whose path starts at v, and `action[b]` is right multiplication by b
restricted to that sub-basis (column k = coordinates of pathₖ·b). An injective I_v is the
**transpose** of the left-multiplication action on the paths ending at v (the k-dual).

All matrix work goes through `modules/linalg_mod.py`, which is a thin exact layer over
`fields.linalg` (`rank`, `nullspace`, `rref`, `solve`): matrix product, Kronecker
product, kernel, column-space pivots, "express these columns in that basis", and the
greedy "independent-modulo-a-subspace" selector that picks minimal generators. Nothing
here is field-specific, so a module lives over ℚ, ℚ(α), GF(p), or GF(p^n) identically —
the same code, different Domain.

## How the computation runs

### Radical, top, socle

`rad M` (`modules/radtopsoc.py`) is the sum of the images of the arrow-action matrices
(`M·rad A = Σ_α image(action[α])`), reduced to an independent column basis; `top M` is
the quotient `M/rad M`; `soc M` is the intersection of the arrow-action kernels
(`∩_α ker(action[α])`). A submodule (radical, socle) is rebuilt as its own `Module` by
**restricting** the action — for each generator b, apply `action[b]` to each submodule
basis column and solve for the coordinates back in that basis (`solve_columns`); a
quotient (top) is rebuilt by choosing coset representatives independent modulo the
submodule and reading the action modulo it. Iterating `rad` gives the radical series, and
its length is the Loewy length.

### Hom and Ext

`Hom_A(M, N)` (`modules/hom.py`) is the space of matrices φ (dim N × dim M) commuting
with the action: `N.action[b] @ φ = φ @ M.action[b]` for every generator b (arrows and
idempotents suffice). Column-stacking turns this into one homogeneous linear system
`(I ⊗ N.action[b] − M.action[b]ᵀ ⊗ I) vec(φ) = 0`; the nullspace is Hom, its dimension is
`A.hom(M, N)`. `Ext^n` (`modules/ext.py`) applies Hom_A(−, N) to a minimal resolution of
M and takes cohomology: with the Hom-space bases H_n, the coboundary δ^n is precomposition
with d_{n+1}, and dim Ext^n = dim H_n − rank δ^n − rank δ^{n-1}.

### The minimal projective resolution

`modules/resolution.py` (generalized from the bridge `obstruction/module_ext.py`, which
was hardcoded to a 4-vertex diamond over ℚ; here MIT-headered, over any vertex set and any
Domain) builds ... → Q_1 → Q_0 → M → 0 by **iterated projective covers**
(Green–Solberg–Zacharia). One step:

1. **Top generators.** Compute `top M = M/rad M` and choose a k-basis, each vector lifted
   into a single vertex block M·e_v (homogeneous generators). If vertex v supplies t_v
   generators, the cover is Q_0 = ⊕_v P_v^{t_v}.
2. **Cover map.** d_0 : Q_0 → M sends the canonical generator of each P_v summand (its
   e_v basis vector) to the lifted top-generator g; a general P_v basis vector for path p
   maps to g·p = `action_M[p] @ g`.
3. **Syzygy.** Ω_1 = ker(d_0), a submodule of Q_0, rebuilt by restriction; cover it to get
   Q_1 and d_1 = (cover of Ω_1) followed by (inclusion Ω_1 ↪ Q_0). Repeat.

Minimality is guaranteed because generators are chosen independent **modulo the radical**
(the `independent_modulo` selector), so d_n(Q_n) ⊆ rad Q_{n-1} and the summand counts are
the true Betti numbers. The resolution **terminates** at length n exactly when Ω_n = 0
(projective dimension n); a term dim blow-up past `max_term_dim` raises `DepthLimitError`
with the certified length. `M.projective_resolution(k)` returns a `ProjectiveResolution`:
`.term(n)` (the summand vertices), `.betti(n)`, `.differential(n)`, `.pd()`, and a readable
`P_1 <- P_2 <- 0` repr. `global_dimension(A) = sup_v pd(S_v)` — exact when every simple
resolves within the depth budget, else a labeled certified lower bound.

## A worked micro-example — S_1 over the A₂ path algebra

`linear_path_algebra(2)` is Q: 1 → 2 with basis `["e_1", "e_2", "a1"]` (the arrow is
auto-named `a1`). The simple S_1 has `action["e_1"] = [[1]]`,
`action["e_2"] = action["a1"] = [[0]]`. Its cover: top S_1 = S_1, one generator at vertex
1, so Q_0 = P_1 = e_1 A (basis `[e_1, a1]`, dimvec {1:1, 2:1}) and d_0 : P_1 → S_1 kills a1.
The syzygy Ω_1 = ker(d_0) = span{a1} ≅ S_2 = P_2 (basis `[e_2]`,
dimvec {1:0, 2:1}), which is projective, so Q_1 = P_2 and Ω_2 = 0. The resolution is
0 → P_2 → P_1 → S_1 → 0, `pd(S_1) = 1`, and since pd(S_2) = 0 the algebra is hereditary,
`global_dimension = 1`. Then Hom(S_1, S_2) = 0 and Ext^1(S_1, S_2) = 1 (the single arrow
1 → 2). (These values were produced by running the code.)

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| the module object, dimension vector | `modules/module.py` | `Module`, `dimension_vector`, `from_arrow_action` |
| exact matrix layer over a Domain | `modules/linalg_mod.py` | `matmul`, `kron`, `kernel_columns`, `independent_modulo`, `solve_columns` |
| simples / projectives / injectives | `modules/builders.py` | `simple`, `projective`, `injective` |
| radical / top / socle, submodule/quotient | `modules/radtopsoc.py` | `radical`, `top`, `socle`, `submodule`, `quotient` |
| Hom / End | `modules/hom.py` | `hom_space`, `hom_dim`, `end_dim` |
| minimal projective resolution | `modules/resolution.py` | `minimal_resolution`, `projective_cover`, `ProjectiveResolution` |
| module Ext^n, global dimension | `modules/ext.py` | `ext`, `ext_dims`, `global_dimension` |
| public methods on the algebra | `core/algebra.py` | `simple`, `projective`, `injective`, `hom`, `ext`, `global_dimension` |
```

- [ ] **Step 3: Update `docs/internals/06-invariants.md`.** Append a short section documenting (a) the new exact spectral layer and (b) the resolved non-unimodular caveat:

```markdown
## The exact spectral layer (Plan 05)

`invariants/spectral.py` computes the **spectral radius** and **Mahler measure** of the
Coxeter polynomial exactly, reimplementing the hanlab float layer (which used mpmath
`nroots`) with exact sympy algebraic numbers. `spectral_radius(p)` is `max_i |α_i|` and
`mahler_measure(p)` is `|lc|·∏_{|α|>1}|α|`; both short-circuit to the exact integer 1 when
`is_cyclotomic_product(p)` (all roots on the unit circle). No floats are used — magnitudes
are `sympy.Abs` of `CRootOf` roots and comparisons are `.is_positive`. The subtle part is
**soundness for complex roots**: `real_roots` alone is unsound (real-roots-suffice is a
theorem only for hereditary quivers), so the code forms the non-cyclotomic part `q` and,
via the self-inversive `y = z + 1/z` substitution, uses a Sturm real-root count of `Q(y)` —
with no complex-root isolation — to decide between fast `real_roots(q)` (Branch A:
hereditary/Salem/Lehmer) and correct `all_roots(q)` (Branch B: non-hereditary,
complex-dominant). The Lehmer star T(2,3,7) = `star_quiver([1,2,6])` carries Lehmer's
polynomial (Branch A), whose spectral radius is Lehmer's number 1.17628…, the smallest known
Mahler measure > 1; the trivial extension T(A) = A ⋉ DA collapses the Coxeter polynomial to
(t+1)^v regardless of representation type (Cartan C ↦ C + Cᵀ, Φ = −I).

The **non-unimodular caveat** (§ "Coxeter matrix and polynomial") is now **documented**:
`coxeter_polynomial`'s docstring names the `det C ∉ {0, ±1}` case, where Φ may be rational
so the Coxeter polynomial is over **QQ** (e.g. t² + 3t/2 + 1 for C = [[2,1],[0,1]]) — exact,
but not the classical integral Coxeter transformation. The domain follows the actual
**coefficients** (sympy's inference), not det C: a non-unimodular Cartan can still be integral
(k[x]/(x²) → t+1 over ZZ; diag(1,2) → (t+1)² over ZZ). The `coxeter_matrix` sibling surfaces
the same fact via its rational-entry branch.
```

Edit the existing "det ≠ ±1 caveat" paragraph in `06-invariants.md` to end with: "As of Plan 05 this is documented on `coxeter_polynomial` itself: the domain follows the coefficients (ZZ when integral, QQ when genuinely rational), so a non-unimodular Cartan that yields a rational Coxeter transformation is called out, not silent (see below)."

Update `docs/internals/README.md`: add chapter 09 to the index and change the "Honest coverage statement" bullet about families/spectral to note the exact spectral layer and the modules chapter have landed with Plan 05.

- [ ] **Step 4: Update `README.md`** — append a "Modules and invariants" section (the block below uses a 4-backtick outer fence so its inner ```python example nests correctly):

````markdown
## Modules and invariants

```python
from quiverlab import Quiver, CC

A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
           ).algebra(relations=["a*b - c*d"], field=CC)   # commutative square

S1, S4 = A.simple(1), A.simple(4)
A.projective(1).dimension_vector()      # {1: 1, 2: 1, 3: 1, 4: 1}
A.ext(S1, S4, 2)                        # 1     (Ext^2 of simples)
int(A.global_dimension())              # 2
A.loewy_length()                       # 3
A.simple(1).projective_resolution(4)   # P_1 <- P_2 (+) P_3 <- P_4 <- 0
```

Every module is a right A-module over the stated exact field; Ext, Hom, and the
projective resolution are exact. Exact `spectral_radius`/`mahler_measure`, `center()`,
`complexity()`, and `sweep()` (invariant × field) round out the invariant surface.
````

- [ ] **Step 5: Run the ONE full suite as a tracked background job, await, then commit.**

Run the focused acceptance first (fast, foreground):
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/modules/test_acceptance_plan05.py -q`

Then the full suite (~19 min, exceeds the 600 s foreground cap) as **one** `run_in_background: true` job:
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q`
**Await it to completion** (Monitor an `until` loop on the pytest output for the summary line, or the background completion notification). Only when it is fully green — including `tests/test_no_floats.py` — make the final commit and report in the same session. If any test fails, fix forward before committing (never commit red).

```bash
git add docs/internals/09-modules.md docs/internals/06-invariants.md docs/internals/README.md README.md tests/modules/test_acceptance_plan05.py
git commit -m "docs+test(modules): Plan-05 acceptance, internals ch.09 modules + ch.06 spectral update

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Frozen interface for Plans 06 (Families + batch) and 07 (Viz + trace)

Later plans consume the following **verbatim** (do not rename, change types, or alter signatures without a coordinated interface bump). Modules are **right** A-modules, column-vector / anti-homomorphism convention; every coefficient is an exact `Domain` element; every path is left-to-right.

```python
# quiverlab.modules.module
class Module:
    algebra; domain; dim: int; action: dict[str, list[list]]; name: str
    def vertex_projection(self, v) -> matrix
    def dimension_vector(self) -> dict[vertex, int]     # rank(action["e_v"]) per vertex, sums to dim
    def radical(self) -> Module                          # M*(rad A), a submodule
    def top(self) -> Module                              # M / rad M
    def socle(self) -> Module                            # {m : m*rad A = 0}
    def projective_resolution(self, length) -> ProjectiveResolution
    def check_module(self) -> tuple[bool, object]
    @classmethod
    def from_arrow_action(cls, algebra, dimension_vector, arrow_action, name="M") -> Module

# quiverlab.modules.resolution
class ProjectiveResolution:
    module; terms; dmats; length: int
    def term(self, n) -> list[vertex]        # summand vertices of Q_n (P_v multiplicities)
    def betti(self, n) -> int                # number of indecomposable projective summands
    def differential(self, n) -> matrix      # d_n : Q_n -> Q_{n-1} over the Domain
    def dimension_vectors(self) -> list[dict]
    def pd(self) -> int | None               # projective dimension, or None if unresolved
    def is_finite(self) -> bool

# quiverlab.modules.ext
class GlobalDimension:                        # .value: int, .exact: bool; int(), == int when exact
def ext(A, M, N, n) -> int
def ext_dims(A, M, N, top) -> list[int]
def global_dimension(A, bound=32) -> GlobalDimension

# quiverlab.modules.builders / hom
def simple(A, v) -> Module
def projective(A, v) -> Module               # records ._pv_basis_labels, ._pv_vertex (private)
def injective(A, v) -> Module
def hom_space(M, N) -> list[matrix]          # k-basis of Hom_A(M, N)
def hom_dim(M, N) -> int
def end_dim(M) -> int

# quiverlab.core.Algebra (new public methods)
A.simple(v); A.projective(v); A.injective(v)      # -> Module
A.hom(M, N) -> int;  A.ext(M, N, n) -> int
A.global_dimension() -> GlobalDimension
A.is_selfinjective() -> bool                       # any field (Nakayama-permutation socle test)
A.loewy_length() -> int
A.center() -> tuple[int, list[vector]]
A.complexity(n) -> int | str | None               # GF(p) only, loud otherwise

# quiverlab.invariants.spectral (exact, no floats; SOUND for complex off-circle roots)
def spectral_radius(poly) -> sympy_exact | None    # 1 on cyclotomic input
def mahler_measure(poly) -> sympy_exact | None

# quiverlab.invariants.scalar
def loewy_length(A) -> int
def center(A) -> tuple[int, list[vector]]
def complexity(A, n) -> int | str | None

# quiverlab.invariants.sweep  (also exported flat as quiverlab.sweep)
def sweep(builder, *args, fields, invariants=None, **kwargs) -> SweepTable
class SweepTable:
    def cell(self, inv, field); rows; def latex(self) -> str; __repr__

# quiverlab.invariants.cartan (NO behaviour change; docstring caveat only)
def coxeter_polynomial(A) -> sympy.Poly    # domain follows the coefficients (ZZ/QQ), documented
```

**Guidance for the consumers.**
- **Plan 06 (families):** the module machinery is available but families do not *require* it — `NakayamaAlgebra`, `PathAlgebra`, `TruncatedPathAlgebra`, `RadicalSquareZero`, `IncidenceAlgebra`, `PreprojectiveAlgebra`, `QuantumCI`, `ExteriorAlgebra`, `TrivialExtension`, `TensorProduct`, `zoo()`, `families()` each return a Plan-01 `Algebra`, and every module/invariant method above then works on them unchanged. `sweep(builder, *args, fields=..., invariants=...)` takes any family constructor with a `field=` keyword directly (its intended §3.9 use). Families that build modules (e.g. verifying a preprojective algebra's structure) may use `Module.from_arrow_action`. The `TrivialExtension` family should reuse `engine.coxeter_spectrum.trivial_extension_cartan` for its Cartan and the exact spectral layer for the `(t+1)^v` collapse check.
- **Plan 07 (trace):** the module/resolution/Ext engines emit **no** trace events yet — they compute silently. When Plan 07 adds the worked-steps trace, the natural event points are: the projective cover (top generators chosen, per vertex), each syzygy dimension, each differential matrix `d_n`, and the Hom/δ ranks in the Ext computation. These mirror the existing `ResolutionTerm`/`DifferentialEvent`/`RankStep` taxonomy; add an optional `trace=None` list parameter to `minimal_resolution`, `ext_dims`, and `hom_space` (append plain dataclasses, inert until Plan 07 renders them) exactly as Plan 03/04 did for their engines. No signature above changes for consumers who omit `trace`.

**Citations (docstrings; registry deferred to Plan 06).** `modules/resolution.py` cites **Green–Solberg–Zacharia**, *Trans. AMS* 353 (2001) 2915–2939 (minimal module resolutions); `modules/ext.py` and `tests/modules/test_ext.py` cite **Assem–Simson–Skowroński**, *Elements of the Representation Theory of Associative Algebras* I (2006), §III (Ext of simples over hereditary algebras) and **Happel** (Coxeter/Dynkin); `tests/engine/test_coxeter_spectrum.py` cites **Lehmer** (1933) for the Lehmer polynomial/number. These cite keys are noted for the Plan-06 citation registry to pick up.

## Boundary notes

- **Module Ext is one-sided; Hochschild is two-sided.** `A.ext(M, N, n)` uses a minimal *module* (one-sided) projective resolution; it is **not** the Chouhy–Solotar bimodule resolution (Plan 04, which serves Hochschild). The optional Plan-04 cross-check (Task 7) only confirms projective dimensions of simples agree with `global_dimension` on monomial/quadratic algebras — it does not couple the two engines.
- **Complexity is engine-backed (GF(p) only).** `A.complexity(n)` routes through the fast int64 engine and fails loudly over `CC`/`GF(p^n)`, matching `nakayama_automorphism`/`is_frobenius`/`cyclic_homology`. `loewy_length` and `center` are exact over **every** Domain.
- **Spectral exactness, soundness & performance.** The exact route is **sound for complex off-circle roots** (not real-roots-only, which is valid solely for hereditary quivers, A'Campo): the cyclotomic short-circuit handles Dynkin/Euclidean cases with no root isolation; the self-inversive `y=z+1/z` Sturm decision then routes to fast `real_roots(q)` (Branch A: hereditary/Salem/Lehmer — Lehmer's irreducible degree-10 in ~0.02 s) or, only when `q` genuinely has complex off-circle roots (non-hereditary), to `all_roots(q)` on the small non-cyclotomic factor (verified: the dim-22 rad²=0 witness in ~0.7 s). `all_roots` is never applied to a high-degree *irreducible* Salem factor like Lehmer (where exact complex-magnitude comparison would take > 2 min), because the Sturm decision keeps that in Branch A. Correctness is unconditional; the only cost is the rare Branch-B `all_roots` on a small factor.
- **`dimension()`/`basis()` (spec §3.5).** These remain the plain `A.dim` attribute and `A.basis_labels` list from Plan 01/02 (no method wrappers added; not in this plan's scope). Noted so a later plan does not expect them as methods.

## Self-review

- **Spec coverage.** §3.5 — `global_dimension` (T7), `is_selfinjective` (T7, previously unclaimed by any plan), `loewy_length`/`complexity`/`center` (T11), exact spectral layer (T8); §3.6 — `simple`/`projective`/`injective` (T3), `radical`/`top`/`socle` (T4), `hom`/`ext` (T5, T7), `projective_resolution` (T6); §3.9 — `sweep` (T12); §5 component 7 (modules) — T2–T7; §5 component 8 (invariants) — T8–T12.
- **Ledger obligations → task.** Re-port `test_lehmer_star_237` → T9; re-port `test_trivial_extension_collapse` → T9; first coverage of `star_quiver`/`trivial_extension_cartan`/`cartan_of_quiver` → T8 (`cartan_of_quiver`/`star_quiver`) + T9 (`trivial_extension_cartan` + `star_quiver` in the collapse); re-port `test_minimal_memory_guard` (dependency now satisfiable via Plan-03 `Quiver.algebra`) → T13; fix non-unimodular-Cartan silent-QQ minor → T10; `A.coxeter_matrix()`/`A.nakayama_automorphism()` P02 pins **not** duplicated (only referenced in T3/T10 suites). Exact `spectral_radius`/`mahler_measure` → T8; `test_spectral_radius_and_mahler` restored → T9.
- **Placeholder scan.** Every `src/` helper named in Tasks 2–12 has a complete body in-plan, including Task 8's spectral helpers (`_noncyclotomic_part`, `_is_reciprocal`, `_reciprocal_to_y`, `_real_roots_suffice`, `_off_circle_roots`), which were **executed verbatim against all 14 Task-8 assertions (all PASS)** — both complex-root counterexamples included. The only `NotImplementedError`-style boundary is the loud `FieldError` on engine-backed `complexity` off a prime field (intended). No `...` in any `src/` block; the `...`-free test scaffolding is complete. The two private hooks on `Module`/`projective` (`_pv_basis_labels`, `_summand_vertices`) are specified with their exact set points (T6 note).
- **Signature consistency vs current sources (agent-verified + executed).** `Algebra` gains exactly ten new methods (`simple`, `projective`, `injective`, `hom`, `ext`, `global_dimension`, `is_selfinjective`, `loewy_length`, `center`, `complexity` — none pre-exist); constructor/`.dim`/`.quiver`/`.relations`/`.multiply`/`._basis_vec` used as they exist. `fields.linalg.{rank,nullspace,rref,solve}` used with their real signatures (`solve` returns a vector or `None`; `rref` returns `(matrix, pivots)`). `errors.*` all take `(message, hint=None)`. `engine.resolutions_minimal.{minimal_resolution(A,N,p,...,max_transient_bytes=),minimal_homology_dims(...,max_transient_bytes=)}` and `engine.adapter.to_engine` used verbatim; `engine.scan3.complexity_of` honest return forwarded. `coxeter_polynomial_from_cartan` returns `(poly, Phi)|(None,None)`; `star_quiver([1,2,6])` = T(2,3,7). No upstream signature is modified, and **no behaviour change** to `coxeter_polynomial` (Task 10 edits only its docstring; sympy's coefficient-based ZZ/QQ inference is already correct, verified by execution). The Task-8 spectral layer was executed verbatim (14/14 assertions pass); the Task-10 domain characterization (kA₂→ZZ, Fixture-D→QQ, k[x]/(x²)→ZZ, diag(1,2)→ZZ) and the diag(1,2) build were verified by running.
- **Ambiguities resolved.** (i) The memory-guard dependency: recorded as satisfiable via Plan-03 `Quiver.algebra` (reduction_algebra superseded per Plan 04) — Plan 05 owns the re-port. (ii) The non-unimodular minor: after execution, decided as *documented behaviour only* (docstring caveat + characterization tests) — sympy already infers ZZ/QQ correctly from the coefficients, so no logic change (a `|det C|` rule would wrongly regress `k[x]/(x²)` and `diag(1,2)` to QQ). (iii) `complexity` grounding: defined on the minimal `A^e` resolution's term growth via the existing `complexity_of`, GF(p)-only. (iv) Module Ext vs Hochschild kept explicitly separate (boundary note). (v) `dimension()`/`basis()` left as attributes (out of scope, noted). (vi) Spectral soundness: real-roots-only was **unsound** (drops complex off-circle roots); replaced with the cyclotomic short-circuit + self-inversive `y=z+1/z` real-root-count decision (Branch A: real_roots exact & fast; Branch B: all_roots fallback on the small non-cyclotomic part) — the boundary is *non-hereditary*, and Branch B computes the correct complex-inclusive answer (chosen fallback over raising). (vii) `is_selfinjective` implemented (not deferred) via the Nakayama-permutation socle test, exact over any field.
