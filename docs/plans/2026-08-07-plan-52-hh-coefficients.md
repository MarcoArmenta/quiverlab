# Plan 52 — Hochschild (co)homology with bimodule coefficients + relative HH (R4)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task.
> Every step uses checkbox (`- [ ]`) syntax. Do each TDD cycle in order: write the
> failing test, run it, see the stated FAIL, implement, run, see the stated PASS,
> commit (conventional commit, green suite).

**Goal:** quiverlab computes `HH^•(A, M)` and `HH_•(A, M)` for an **arbitrary
finite-dimensional `A`-bimodule `M`** — the named coefficient bimodules (the
regular `A`, the dual `D(A) = Hom_k(A,k)`, a twisted `{}_φA_ψ` including
`{}_1A_ν` with `ν` the Nakayama automorphism, the quotient `A/soc_{A^e}A`) and a
generic no-code `(dims, left maps, right maps)` bimodule — over **every exact
`Domain`** (ℚ, ℚ(α), GF(p), GF(pⁿ)), exact-only, loud on the scope boundaries.
Plus **relative Hochschild (co)homology `HH_•(A|B)`** for `B = kQ₀` (the
separable vertex subalgebra) via the `E`-relative reduced bar resolution, with a
loud typed refusal for a general subalgebra `B`. The coefficient enters **only at
collapse time** (`Hom_{A^e}(P_•, −)` / `− ⊗_{A^e} P_•`); the resolutions of `A`
over `A^e` are coefficient-independent and untouched.

**Architecture.** A new **exact-`Domain`** bimodule container
`src/quiverlab/hochschild/coefficients.py::Bimodule(A, dim_M, Lact, Ract)`
(left- and right-action tables over `A`, compatibility-checked, exact entries) —
the exact twin of the already-ported GF(p)/int64
`src/quiverlab/engine/bimodule.py` (`Bimodule(alg, mu, Lact, Ract)`,
`regular_bimodule`/`dual_bimodule`/`twisted_bimodule`,
`hochschild_{co}homology_with_coefficients`) — which has a live consumer
(`engine/coxeter_spectrum.py:57`) but is **not wired into the public HH
dispatch**, and becomes a **read-only GF(p) reference oracle** here. The coefficient is threaded through the existing public
dispatch `core/algebra.py::hochschild_cohomology`/`hochschild_homology` by a new
`coefficients=None` kwarg (None ≡ regular bimodule, **byte-identical**
short-circuit — a hard gate). Three engine collapses gain the coefficient swap:
**bar** (`hochschild/bar.py`, the reference, any Domain), **Chouhy–Solotar**
(`resolutions_cs/homology.py` + `resolution.py::matrix`, the deep any-Domain
route), and the internal **minimal `A^e`** collapse (`engine/resolutions_minimal.py`,
GF(p), the Plan-16 covariance generalized). The **fast** GF(p) bar-basis route is
hard-wired to the regular bimodule and **refuses coefficients loudly**. Relative
`HH_•(A|kQ₀)` runs the `E`-relative reduced bar complex, sharing the Task-3
coefficient collapse.

**Tech stack:** Python ≥ 3.10; deps unchanged (`numpy`, `sympy`); exact
arithmetic only, no floats in `src/`. No new hard dependency.

**Status / change-log.**
- 2026-08-07 — plan authored (spec layer = R4; house-style TDD).
- 2026-08-07 adversarial review: 8 majors + 6 minors applied (API grounded to
  the real surface — char-0 field is `ql.CC` (the nonexistent `QQ` symbol
  removed), the fictitious engine-ν helper deleted, and the `Algebra` ctor
  corrected to `Algebra(domain, T, unit)`; semantic freshness gate replacing the
  substring match; ν-basis
  reconciliation DD1b + multi-vertex GF(p) discriminator; nontrivial-ν coverage
  via `NakayamaAlgebra(3,3,cyclic)` in the twisted/CS/minimal batteries;
  non-tautological duality oracle reclassified `oracle_selfcert`; estimator
  `sizing_dim` coefficient term + schema-v3 reconciliation with conditional
  versioning; chain-level `twisted_homology_classes` P54 cross-plan contract;
  reference/framing corrections — `engine/bimodule.py` has a live consumer,
  "degreewise" not "byte", 2411.03080 is the in-scope vertex-relative setting).

---

## Record (R4, verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R4 — HH with arbitrary bimodule coefficients + relative HH•(A|B).** [A-scout
> P8; keep — verified genuinely absent from the public surface]
> Object: HH•(A,M), HH_•(A,M) for any f.d. bimodule M (twisted _1A_σ, D(A),
> A/soc…); relative HH over a subalgebra B via the relative bar resolution. Refs:
> Chaparro–Schroll–Solotar arXiv:1811.02211 (gentle HH¹ with coefficients);
> Lindell–Rubio y Degrassi arXiv:2411.03080. Shallow engine extension (swap the
> coefficient in Hom_{A^e}(P•, −)); prerequisite for R2/R8/R9. Size M.

**Wave/consumer context (metaplan §4–5, P52 card).** P52 is Wave-1 tier-α,
independent, and is the **prerequisite of three downstream plans** whose
coefficient needs this plan's API must serve by name:

- **P54 (R2, BV operator)** — needs the **twisted bimodule `{}_1A_ν`** with `ν`
  the Nakayama automorphism (`D(A) ≅ {}_1A_ν` for self-injective `A`); the BV
  operator transports Connes `B` through `σ`-twisted duality. Served by
  `Bimodule.twisted(A, psi=ν)` / `Bimodule.twisted_by_nakayama(A)` and
  `HH^n(A, {}_1A_ν)`. **NOT served by the dims-only surface alone** — the BV
  transport needs **chain-level** twisted (co)homology-class representatives, so
  P52 additionally ships the accessor `twisted_homology_classes` (Task 2, the
  cross-plan contract W2 below).
- **P72 (R5+R6, split-extension LES + arrow removal)** — needs **`H^*(B, M)`**,
  ordinary Hochschild cohomology of an algebra `B` with a `B`-bimodule
  coefficient `M`. Served directly by `B.hochschild_cohomology(top,
  coefficients=M)`.
- **P74 (R8, skew group algebras `A⋊G`)** — needs **`HH^n(A, {}_gA)`**, the
  `g`-twisted bimodule for `g ∈ G` (Ştefan conjugacy-class decomposition).
  Served by `Bimodule.twisted(A, phi=g_matrix)`.
- Also flagged prerequisite of the R9 incidence route (P75).

The coefficient API is designed so **P72/P74/P75 are served by the dims + block
surface without change**; **P54 additionally requires a chain-level accessor**
(next). Each consumer is named again at the interface where it lands (Task 2).

**Cross-plan contract with P54 (chain-level twisted (co)homology data).** BV needs
more than dimensions: it needs the actual class representatives and boundary maps
in twisted bar-chain coordinates to transport Connes `B`. P52 therefore ships a
public accessor `hochschild/coefficients.py::twisted_homology_classes(A, M, top)`
returning, per degree `n ≤ top`: the twisted bar-chain **basis** (the `(s, J)`
enumeration of `M ⊗ \bar A^{⊗n}`), the **boundary matrices** `b_n` (as exact
Domain matrices), and a **basis of homology-class representatives** (a chosen
lift of `ker b_n / im b_{n+1}`). It carries a **CONVENTION BLOCK, to be mirrored
VERBATIM from P54's plan §3** (P54 is amending §3 in its own fix round to contract
exactly this — do not diverge):

> *Twist placement:* the coefficient `M` occupies the leading tensor factor
> `M ⊗ \bar A^{⊗n}`; the twist (right action `m·a = m·ψ(a)`) enters the boundary
> **only** on the **LAST face map** `d_n` (`(−1)^n · (ψ-twisted) a_n · m`),
> matching P1; the first face uses the untwisted left action. *Normalization:*
> reduced bar (`\bar A = A/k·1`), so index 0 is dropped from every bar slot.
> *Basis ordering:* `M`-basis index `s` outermost (slowest), then the bar
> multi-index `J ∈ {1..m-1}^n` in lexicographic order (fastest last). *Cohomology
> twin:* `twisted_cohomology_classes` mirrors with `Hom_{A^e}` collapse and the
> transposed convention.

Marked "**cross-plan contract with P54**"; a divergence between this block and
P54 §3 is a review defect on whichever plan merges second.

---

## Reference re-verification (mandatory; findings recorded 2026-08-07)

Re-verified via web search at spec time (standing rule, metaplan §1.5).

1. **arXiv:1811.02211 — Chaparro, Schroll, Solotar,** *"On the Lie algebra
   structure of the first Hochschild cohomology of gentle algebras and Brauer
   graph algebras,"* **J. Algebra 558 (2020) 293–326.** CONFIRMED authorship,
   title, venue. The abstract confirms the paper **"determine[s] the first
   Hochschild homology and cohomology with different coefficients for gentle
   algebras"** and gives a ribbon-graph/Brauer-graph interpretation — so the
   record's "gentle HH¹ with coefficients" is accurate and on-topic. **BUT** the
   explicit `HH¹(A, M)` values are gentle-family formulas stated through the
   ribbon-graph combinatorics; a clean transcribable numeric pin is **not
   extractable from the abstract**, and a direct PDF fetch returned only the
   compressed byte stream. **Ruling: the 1811.02211 literature pins are
   `BLOCKED-until-transcribed`** — Task 8 transcribes one concrete small worked
   example (quiver + relations + coefficient + the integer `dim HH¹`) verbatim
   from the paper's PDF/TeX *before* any test asserts it, with the equation
   number; until then the internal oracles (Task 6) carry certification. This
   matches the repo's standing posture for hard-to-transcribe references (several
   `docs/verification.md` honest-scope entries are cross-engine-only for exactly
   this reason).

2. **arXiv:2411.03080 — Jonathan Lindell & Lleonard Rubio y Degrassi,** *"On the
   first relative Hochschild cohomology and contracted fundamental group"*
   (28 pp., submitted 2024-11-05, math.RT/math.KT). CONFIRMED authorship and
   title. Content: the Lie-algebra structure of the **first relative Hochschild
   cohomology** and its relation to a relative fundamental group; explicit
   computations for **radical-square-zero** algebras and **dual extension
   algebras of directed monomial** algebras. This is the reference for the
   **relative** part (Task 7); its `HH¹`-relative values are likewise
   `BLOCKED-until-transcribed` (Task 7/8 transcribes the radical-square-zero
   example if used). The paper also motivates the general-`B` follow-up.

Both keys are added to `references.bib` + the citation registry (Task 10) and
surfaced in results regardless of whether a numeric pin lands (the *provenance*
is real even where the *value* stays cross-engine-certified).

---

## Design decisions (with rationale)

### DD1 — What a bimodule IS in quiverlab: a new light exact-`Domain` `Bimodule` container (NOT a module over `A^e`)

**Decision.** Introduce
`src/quiverlab/hochschild/coefficients.py::Bimodule`, an exact-`Domain` container
storing, for a `core.Algebra A` of dimension `m`:

- `dim_M` — `dim_k M`;
- `Lact` — a length-`m` list of `dim_M × dim_M` Domain matrices; `Lact[j] @ v`
  is `b_j · v` (left action by the `j`-th `A`-basis element);
- `Ract` — a length-`m` list of `dim_M × dim_M` Domain matrices; `Ract[j] @ v`
  is `v · b_j` (right action);
- `domain`, `name`, and a cached vertex-grading `corner(v, w)` → the basis
  indices of `e_v M e_w`.

**Rationale.** The two candidate representations:

- **(a) module over the enveloping algebra `A^e = A ⊗ A^op`.** Gets
  Krull–Schmidt / Hom / decompose "for free" from the module surface — *but*
  `A^e` is a **quadratic blow-up** (`dim A^e = (dim A)²`: the zoo's dim-9 square
  → 81; the minimal engine's implicit `A^e` arithmetic tensor is `(m,m,m²,m²)`,
  i.e. `m⁶` entries, dim-9 → ~531 k). Crucially **quiverlab never materializes
  `A^e` as an `Algebra`** — `TensorProduct(A, A.opposite())` exists
  (`families/tensor.py`) but is never called for the enveloping algebra, and
  every engine works **corner-wise** (`e_v A e_w`) or with **`(Lact, Ract)`
  pairs**, precisely to avoid holding a dim-`m²` object.
- **(b) a light `(Lact, Ract)`-over-`A` container.** This is the representation
  the codebase already commits to. A **ported GF(p)/int64 twin already exists**
  at `src/quiverlab/engine/bimodule.py` (`Bimodule(alg, mu, Lact, Ract)` with
  `regular_bimodule`, `dual_bimodule`, `twisted_bimodule(alg, phi, psi)`,
  `check_bimodule`, and full `hochschild_{co}homology_with_coefficients` over the
  bar complex `C_n = M ⊗ (A/k1)^{⊗n}`). It is **not wired into the public HH
  dispatch** — but it is read-only and has a live consumer
  (`engine/coxeter_spectrum.py:57` imports it), so P52 treats it as a reference,
  not dead code. It matches exactly what the record calls "swap the coefficient
  in `Hom_{A^e}(P_•, −)`."

We take **(b)**, exact-`Domain`. The engine (co)homology collapses need only two
things from the coefficient — its **left action** and its **right action** of
`A`-elements, plus the **vertex grading** `e_v M e_w` — and `(Lact, Ract)`
provides exactly those with `O(m·dim_M²)` storage, no `A^e`. The existing
`engine/bimodule.py` becomes a **read-only GF(p) reference oracle** (Task 3:
exact path over GF(p) ≡ ported hanlab path, **degreewise** — the assertions
compare dim vectors, not raw matrices). **Dims bound:**
storage `m · dim_M²` Domain entries; for the named bimodules `dim_M = dim A = m`
(regular / dual / twisted) or `dim_M ≤ m` (`A/soc`), so storage is `O(m³)` —
same order as `A.T`, no blow-up. The generic no-code form is capped by a
`_MAX_COEFF_DIM` allocation guard (mirroring `_MAX_MODULE_DIM`).

**Constructors** (all exact, all validated by `Bimodule.check()`):
`Bimodule.regular(A)`, `Bimodule.dual(A)` = `D(A)`, `Bimodule.twisted(A,
phi=None, psi=None)` = `{}_φA_ψ` (convention `a·x·b = φ(a)·x·ψ(b)`, matching
`engine/bimodule.py`), `Bimodule.twisted_by_nakayama(A)` = `{}_1A_ν`,
`Bimodule.mod_socle(A)` = `A/soc_{A^e}A` (reusing
`families/trivial_extension.py::_bimodule_socle`), and
`Bimodule.from_actions(A, dim_M, left_maps, right_maps)` (the generic no-code
form). `phi`/`psi` are algebra automorphisms given as `m×m` Domain matrices
(columns = images of the basis) **in `A`'s OWN basis** (the basis `A.T` is
written in) — see DD1b, which is load-bearing for `twisted_by_nakayama`.

### DD1b — Nakayama-automorphism BASIS reconciliation (P54-critical, latent-bug fix)

**The trap (verified against the source).** `Algebra.nakayama_automorphism()`
is **basis-inconsistent across fields**: over GF(p) it returns `ν` in the
**engine unit-adapted basis** (`core/algebra.py:896-897`:
`nakayama_automorphism(to_engine(self.unit_adapted()), p)`), whereas off GF(p) it
returns `ν` in `A`'s own basis (`nakayama_automorphism_generic(self)`,
`core/algebra.py:900`). `Bimodule.twisted` builds `Lact`/`Ract` against `A.T`,
which is in `A`'s **own** basis. Feeding the GF(p) engine-basis matrix into
`twisted` is a **silent correctness bug** whenever `unit_adapted()` is a genuine
change of basis (multi-vertex, or any algebra whose stored basis is not already
unit-adapted).

**Decision (single, binding).** `Bimodule.twisted_by_nakayama(A)` obtains `ν`
via the **field-general trace-form route in `A`'s own basis** —
`invariants/frobenius.py::nakayama_automorphism_generic(A)` (the exact
`G⁻¹Gᵀ`-style certifier from Plans 19/29; returns an `m×m` Domain matrix, columns
= images of `A`'s basis elements, in `A`'s basis) — **not**
`Algebra.nakayama_automorphism()`. This is field-uniform (GF(p) and CC take the
same code path) and needs no change-of-basis. *Documented alternative, for the
record only:* had we used the engine matrix `ν_eng`, the correct conversion is
`ν = P⁻¹ · ν_eng · P` where `P` is the `unit_adapted()` change-of-basis matrix
(`A.unit_adapted()` exposes it); we do NOT take this path (a needless extra
solve). A **multi-vertex GF(p) discriminating test** where `unit_adapted()` is a
genuine change of basis (Task 2) fails loudly if the engine-basis matrix is ever
used by mistake.

### DD2 — Engines: bar (reference) + CS (deep) required; minimal `A^e` internal cross-check; fast refuses

- **bar** (`hochschild/bar.py`) — the **reference implementation**, any Domain,
  presentation-free. Primary deliverable. Injection points (verified): in
  `coboundary_matrix` the coefficient enters at `A.T[K[0]][s]` (left, term 0) and
  `A.T[s][K[n]]` (right, last term); in `boundary_matrix` at `A.T[s][J[0]]`
  (right, term 0) and `A.T[J[n-1]][s]` (left, last term); the interior
  `A.T[K[i-1]][K[i]]` is the **algebra multiplication on bar factors and stays**;
  `_cochain_basis` / `cn = m·(m-1)ⁿ` decouple `dim_M` from `m`.
- **Chouhy–Solotar** (`resolutions_cs/`) — the **deep any-Domain** route (needs an
  admissible presentation). Injection: `resolution.py::matrix` collapses through
  `self.ar.mul(c_vec, self.ar.mul(ej, a_vec))` = `b·w·a` (homology) and
  `self.ar.mul(a_vec, self.ar.mul(ej, c_vec))` = `a·w·b` (cohomology), with
  `ej = A._basis_vec(j)`, `j` over the corner `e_t A e_o` / `e_o A e_t`; for `M`,
  `j` ranges over the corner `e_t M e_o` / `e_o M e_t` and the two outer `mul`s
  become `M`'s left/right actions. The differential term-list `d_terms` is
  coefficient-independent and untouched.
- **minimal `A^e`** (`engine/resolutions_minimal.py`) — **GF(p), engine-internal**
  (NOT reachable via public `engine=`, per the standing convention). Delivered as
  the internal `minimal_{co}homology_dims(..., coefficients=M)` for the
  cross-engine oracle and as P54's future GF(p) tool. The Plan-16 covariance
  generalizes cleanly (DD-math below): homology `b·w·a` on corner `e_w A e_v`
  (tag `(v,w)`), cohomology `a·w·b` on the **SWAPPED** corner `e_v A e_w`; for `M`
  replace `cornerA[(i,j)] = e_j A e_i` with `e_j M e_i` and the two structure-
  constant multiplications with `M.Lact`/`M.Ract`. Resolutions unchanged.
- **fast** (`engine/adapter.py` + `scan3`/`hh_engine`) — hard-wired to the
  regular bimodule: `to_engine(A)` carries only `(m, T, unit)`, the coefficient
  basis is literally `range(alg.m)`, and left/right/interior all use the single
  `mult_full == alg.T`. **Refuses coefficients loudly** (`engine="fast"` +
  `coefficients != None` → `QuiverlabError` with a hint pointing at `bar`/`cs`).
  Re-plumbing the fast route for general `M` is out of v1 scope (honest-scope
  entry).

**Cross-engine oracle web (the standing agreement, `oracle_crossengine`):**
`bar ≡ CS` degreewise over **CC and GF(p)** for every coefficient with a
presentation; `bar ≡ minimal` over GF(p); `bar(exact) ≡ engine/bimodule.py`
(ported hanlab) over GF(p) **degreewise** (dim vectors) — three independent
implementations.

### DD3 — Public API: a `coefficients=None` kwarg; None is a byte-identical short-circuit (hard gate)

**Decision.** Extend the two existing public methods (signatures verified):

```python
Algebra.hochschild_cohomology(self, top, max_cells=4_000_000, engine="auto",
                              auto_cs=False, coefficients=None, verbose=None, trace=None)
Algebra.hochschild_homology (self, top, max_cells=4_000_000, engine="auto",
                              auto_cs=False, coefficients=None, verbose=None, trace=None)
```

- **`coefficients=None`** ⇒ the **existing dispatch runs verbatim, byte-identical**
  — the new kwarg threads through but every current branch (`_route_to_cs`,
  `_use_fast_engine`, the Plan-34 depth fallback, the `HHTable` construction, the
  `table.references = self.citations()` frozen contract) is untouched. This is a
  **hard acceptance gate**: every existing golden, every existing HH test, and
  every `_runner_goldens.json` entry stays green with zero re-freezing.
- **`coefficients=M`** (a `Bimodule`) ⇒ route to the coefficient collapse: `bar`
  (any Domain), `cs` (presented), refuse on `fast`. `engine="auto"` with a
  coefficient **never routes fast** (`_use_fast_engine` returns `False` when
  `coefficients is not None`); it routes `bar` in-window and, for a
  quiver-presented algebra, falls back to CS on `DepthLimitError` exactly like
  the regular path (Plan-34 amendment, coefficient-aware). Explicit `engine=`
  keeps honest walls.
- **`HHTable`** gains an optional `coefficients` provenance field (a short
  descriptor, e.g. `"D(A)"`, `"twisted _1A_ν"`, `"A/soc"`, or `None`). Default
  `None` **serializes identically** to today (the field is only surfaced when
  non-None — the Plan-26 "drop-absent" discipline). `HHTable.engine` gains a
  coefficient tag only in the non-regular case.
- A coefficient's algebra **must be `A`** (`M.algebra is A` / equal presentation)
  — a mismatch raises `QuiverlabError` before any compute (mirrors
  `modules/hom.py::_assert_comparable`).

### DD4 — Relative HH: v1 = `B = kQ₀` (separable) only; general `B` refused loudly

**Decision.** Add `relative_to=None` to both methods. `relative_to="vertices"`
computes `HH_•(A|kQ₀, M)` / `HH^•(A|kQ₀, M)` via the **`E`-relative reduced bar
complex** `Bar_n^E = A ⊗_E \bar A_E^{⊗_E n} ⊗_E A` (`\bar A_E = A/E`, `E = kQ₀`) —
Cibils' reduced complex. `relative_to=None` (default) = absolute, unchanged. Any
other `B` (a general subalgebra, however specified) raises
`QuiverlabError("relative HH over B is implemented only for B = kQ_0 in v1", hint=…)`.

**Rationale + verified oracle.** `E = kQ₀ = k^{|Q₀|}` is a **separable**
`k`-algebra (separability idempotent `Σ_v e_v ⊗ e_v ∈ E ⊗_k E`), for **any**
base field and characteristic. Hence each `Bar_n^E` is a **projective
`A^e`-bimodule** and the `E`-relative reduced bar complex is an `A^e`-projective
resolution of `A` — so it computes the **absolute** Hochschild (co)homology:
`HH_•^E(A, M) ≅ HH_•(A, M)` for every `M`. **This is a proven free strong
oracle:** `HH_•(A|kQ₀, M)` must agree **degreewise** with the absolute route
(Task 3). *Nuance verified:* quiverlab's existing normalized bar (`hochschild/bar.py`)
reduces modulo `k·1` only (`\bar A = A/k·1`, index 0 dropped), i.e. it is the
**absolute** normalized complex, NOT the `E`-relative one — so for a **multi-vertex**
algebra the `E`-relative complex is a genuinely **smaller, distinct** complex
that agrees only degreewise (by separability); for a **single-vertex** algebra
`E = k` and the two complexes coincide byte-identically. Task 1's freshness gate
pins `bar.py`'s reduction convention so the implementer builds the right
`\bar A_E`. **Reference placement (corrected):** Lindell–Rubio y Degrassi
(arXiv:2411.03080) work in the **vertex-relative (`E = kQ₀`, separable)** setting
— i.e. **exactly what P52 ships in v1**, not a future general-`B`. Their
radical-square-zero relative-`HH¹` values are therefore an **in-scope** Task-8
oracle (Task 7 substrate), not motivation for the deferred case. The genuinely
deferred follow-up is a **non-separable / non-vertex `B`** (where relative ≠
absolute) — the substrate a later P72 `H^*(B, M)` generalization would use; v1
refuses it loudly.

### DD5 — GUI: builtin coefficient picker delivered; explicit-matrix bimodule editor ledger-deferred

**Decision.** A **new top-level request block `coefficients`** (a `CoefficientSpec`),
**not** smuggled through `HpcConfig` (coefficients are *mathematics*, not resource
tuning, and must enter the cache key). It reuses the Plan-26 `ModuleSpec`
side/entry discipline:

- **builtin (zero-typing, delivered):**
  `{"builtin": {"kind": "regular"|"dual"|"twisted_nakayama"|"quotient_socle"}}` —
  a picklist on the Hochschild kinds. Covers every named bimodule and P54's
  `{}_1A_ν`.
- **explicit (`{dim, left_maps, right_maps}`, one exact-entry matrix per arrow
  per side) — GUI editor LEDGER-DEFERRED** to P80 via an explicit entry in
  `docs/plans/DEEPER-ENGINES-BACKLOG.md` (metaplan §1.2 forbids *silent*
  deferral). The library + server accept the explicit form (Task 2/3); only the
  canvas *matrix editor* for a two-sided bimodule is deferred (the Plan-26 editor
  is one-sided). This keeps P52 sized-M while shipping the full no-code named
  bimodules end-to-end.

**Schema reconciliation (C5 — the current `== 2` gate makes a combined request
unsatisfiable).** `webapp/server/schema.py:282-283` currently requires a
`module`/`ext_target`/`tor_target` block to have `schema_version` **exactly 2**;
a coefficients block needs `>= 3`; a request carrying **both** a module and a
coefficient (e.g. a module-Ext computation with a twisted coefficient — not v1,
but the schema must not forbid it) would be impossible. Fix, precisely:

- change the module gate from `self.schema_version != 2` to
  `self.schema_version < 2` (module blocks accept `schema_version >= 2`);
- the `coefficients` block requires `schema_version >= 3`;
- `_schema_known` / the "known schema" set and its stale message speak **v1 / v2
  / v3** (add 3; update the message text).

**CANONICAL-KEY RULE (verbatim, binding — conditional versioning, the Plan-26
precedent):** *the GUI and every client send `schema_version = 3` **only when a
`coefficients` block is present**; a request without coefficients keeps sending
`schema_version = 2` (or 1), so its `model_dump` and its canonical key are
byte-identical to today.* The `coefficients` field is added to
`ComputeRequest.model_dump`'s absent-block pop-list, so even a v3 request drops
the key when the block is absent — but clients must not gratuitously bump the
version, or an existing example would re-key. Both runners (`hpc/spec.py` +
`docs/gui/runner.py`) carry the same dispatch; i18n ×4; report renders a
**"coefficients: …"** provenance line.

### DD6 — Oracles + Plan-32 markers + QPA

- **M = A ≡ ordinary HH degreewise, every engine** — `Bimodule.regular(A)` through
  the *general* coefficient path must reproduce the `coefficients=None`
  short-circuit dims degreewise (`oracle_crossengine`). This validates that the
  general path specializes correctly.
- **`HH^0(A, M) = M^A` and `HH_0(A, M) = M/[A,M]`** — closed forms computed
  independently from the actions (`oracle_selfcert`).
- **`dim HH^n(A, D(A)) = dim HH_n(A, A)`** — the duality identity, **proven
  exactly** (preamble): `Hom_{A^e}(P_•, DA) ≅ (A ⊗_{A^e} P_•)^*`, cohomology of
  the dual complex = dual of homology over a field, so
  `HH^n(A, DA) ≅ HH_n(A, A)^*`. **`oracle_selfcert`** (Plan-32: a *theory
  identity* the one implementation must satisfy — not two independent
  implementations of the same number; both sides are the same engine's own
  output checked against a proven equation). Made **non-tautological** by driving
  it on a **non-symmetric** self-injective algebra (`NakayamaAlgebra(3, 3,
  cyclic)`, where `D(A) ≇ A`), with `truncated_polynomial(3)` kept only as a
  smoke case (Task 6, C3).
- **Symmetric-algebra consistency web** (`oracle_selfcert`): for symmetric `A`,
  `ν` is inner ⇒ `{}_1A_ν ≅ A` ⇒ `HH^•(A, {}_1A_ν) = HH^•(A, A)`; and
  `D(A) ≅ A` ⇒ `HH^•(A, D(A)) = HH^•(A, A)` both variances — combined with the
  duality identity this forces `dim HH^n = dim HH_n` (the classical symmetric-
  algebra fact), a simultaneous cross-check.
- **GF(p) bank oracle:** exact path ≡ ported `engine/bimodule.py`
  (`oracle_crossengine`, **degreewise** — dim vectors).
- **Relative:** `HH_•(A|kQ₀, M) ≡ absolute` degreewise (separability;
  `oracle_crossengine`), byte-identical single-vertex (`oracle_selfcert`).
- **Literature (both in-scope):** 1811.02211 gentle `HH¹`-with-coefficients
  (Task 8) **and** 2411.03080 radical-square-zero **relative-`HH¹`** — the latter
  is the **vertex-relative (`E = kQ₀`) setting P52 ships**, so it is an in-scope
  Task-8 oracle, not a general-`B` motivation (`oracle_literature`,
  **BLOCKED-until-transcribed** for the numeric values).
- **QPA:** QPA 1.37 has **no** Hochschild-with-coefficients / relative-HH surface
  (expected; probe live at impl time via the `NamesGVars()` sweep, exactly like
  Plan-35 products). The `tests/qpa/` battery **skips honestly and FAILS if that
  ever changes**; the covering oracles are the internal identities + the GF(p)
  bank + the literature pins.

---

## Mathematical preamble (the plan carries the mathematics)

Citations: Cartan–Eilenberg; Weibel *Homological Algebra* Ch. 9; Cibils,
*Cohomology of incidence algebras* / *rigid monomial algebras* (relative theory);
Chaparro–Schroll–Solotar (arXiv:1811.02211); Lindell–Rubio y Degrassi
(arXiv:2411.03080).

### P0. Setup

`A = kQ/I` over an exact `Domain`, `E = kQ₀ = ∏_v k e_v`, `A^e = A ⊗ A^op`. An
`A`-bimodule `M` = an `A^e`-module = a `k`-space with commuting left and right
`A`-actions: `(a·m)·a' = a·(m·a')`, `1·m = m = m·1`. Hochschild:
`HH^n(A, M) = Ext^n_{A^e}(A, M)`, `HH_n(A, M) = Tor_n^{A^e}(A, M)`. Take any
projective `A^e`-resolution `P_• → A`; the coefficient enters **only** through
`Hom_{A^e}(P_•, M)` (cohomology) and `M ⊗_{A^e} P_•` (homology).

### P1. The bar complex with coefficients (the reference; CS §2 normalized bar)

Normalized bar: `C_n(A, M) = M ⊗_k \bar A^{⊗n}`, `\bar A = A/k·1`, on the basis
`(s, J)` with `s` a basis index of `M` and `J ∈ {1..m-1}^n` (bar factors over the
reduced algebra basis). The coefficient occurs at exactly the outer two terms:

```
cohomology  (δ^n f)(a_1⊗…⊗a_{n+1}) = a_1 · f(a_2⊗…⊗a_{n+1})                 [LEFT action on M]
              + Σ_{i=1}^n (−1)^i f(a_1⊗…⊗a_i a_{i+1}⊗…)                      [algebra mult, interior]
              + (−1)^{n+1} f(a_1⊗…⊗a_n) · a_{n+1}                           [RIGHT action on M]
homology     b_n(m⊗a_1⊗…⊗a_n) = m·a_1 ⊗ a_2⊗…                              [RIGHT action on M]
              + Σ_{i=1}^{n-1} (−1)^i m⊗…⊗a_i a_{i+1}⊗…                       [algebra mult, interior]
              + (−1)^n a_n·m ⊗ a_1⊗…⊗a_{n-1}                                [LEFT action on M]
```

Concretely in `bar.py`: replace `A.T[K[0]][s]` → `M.Lact` applied (left, coh
term 0), `A.T[s][K[n]]` → `M.Ract` (right, coh last), `A.T[s][J[0]]` → `M.Ract`
(hom term 0), `A.T[J[n-1]][s]` → `M.Lact` (hom last); interior `A.T[·][·]` stays;
`cn = dim_M · (m-1)^n`. `M = A` (`Lact = Ract = A.T`, `dim_M = m`) recovers the
current code exactly.

### P2. The CS collapse with coefficients (`resolutions_cs`)

`P_n = ⊕_{σ∈S_n} A e_{o(σ)} ⊗ e_{t(σ)} A`; the differential `d_terms(n,σ) =
{(c_i, a_i, τ_i, b_i)}` (paths `a_i, b_i`, chain `τ_i ∈ S_{n-1}`) is
**coefficient-independent**. Collapse against `M`:

```
homology   C_n = ⊕_σ e_{t(σ)} M e_{o(σ)} :   (σ, w) ↦ Σ_i c_i · ( τ_i ,  b_i · w · a_i )
cohomology C^n = ⊕_σ e_{o(σ)} M e_{t(σ)} :   (δf)(τ) = Σ_i c_i · a_i · f(τ_i) · b_i
```

(`b·w·a` homology / `a·w·b` cohomology — the standard `A^e = A ⊗ A^op` op-twist).
For `M = A` this is `resolution.py::matrix` verbatim; the swap is exactly
`w ∈ e_t A e_o → e_t M e_o` and `A.multiply(…)` → `M.Lact`/`M.Ract`.

### P3. The minimal `A^e` collapse with coefficients (Plan-16 covariance, generalized)

The minimal resolution's boundary carries `A^e`-block coefficients `(uu, vv)`
(`uu` = left tensor factor, `vv` = right/op factor). For a coefficient basis
element `w`:

```
homology  (_contracted / _corner_contracted):  w ↦ Σ (e_vv · w) · e_uu   = b·w·a,  corner e_w A e_v  (tag (v,w))
cohomology(_cohomology / _corner_cohomology):   w ↦ Σ (e_uu · w) · e_vv   = a·w·b,  corner e_v A e_w  (SWAPPED tag)
```

**Generalization to `M` (the ultra-hard point):** the covariance and the
tag-swap describe how the *resolution's* `(uu, vv)` factors act — they are
**independent of the coefficient**. So for general `M`:

1. the coefficient block `cornerA[(i,j)] = e_j A e_i` becomes `e_j M e_i` (the
   `(j,i)` graded piece of `M`, image of `m ↦ e_j·m·e_i` over an `M`-basis; the
   local single-vertex path's "full copy of `A`, dim `m`" becomes "full copy of
   `M`, dim `dim_M`");
2. the two structure-constant multiplications become **`M`'s actions**: homology
   `T[vv,·]` → `M.Lact[vv]`, `T[·,uu]` → `M.Ract[uu]`; cohomology `T[uu,·]` →
   `M.Lact[uu]`, `T[·,vv]` → `M.Ract[vv]`;
3. `_solve_in_span` targets the corresponding `M`-corner; the image lands in its
   corner by `A^e`-linearity of the actions — **unchanged assertion**.

The `a·w·b` (coh) vs `b·w·a` (hom) covariance and the SWAPPED-tag block
accounting are **untouched**. This is why the record calls it a "shallow engine
extension": no resolution changes, only the collapse's coefficient slot.

### P4. The named coefficient bimodules

- **regular** `A`: `Lact = Ract = A.T`, `dim_M = m`. (`Bimodule.regular`.)
- **dual** `D(A) = Hom_k(A, k)`: `(a·f·b)(x) = f(b·x·a)`; on the dual basis the
  actions are transposes of `A`'s multiplications. `dim_M = m`. (`Bimodule.dual`;
  must match `engine/bimodule.py::dual_bimodule` over GF(p).)
- **twisted** `{}_φA_ψ`: underlying space `A`, `a·x·b = φ(a)·x·ψ(b)` for algebra
  automorphisms `φ, ψ` (`m×m` matrices, columns = images); `φ=ψ=id` = regular.
  `{}_1A_ν` = `twisted(A, psi=ν)` with `ν` the Nakayama automorphism; for
  self-injective `A`, `D(A) ≅ {}_1A_ν`. (`Bimodule.twisted` /
  `twisted_by_nakayama`; matches `engine/bimodule.py::twisted_bimodule`.)
- **quotient** `A/soc_{A^e}A`: `soc_{A^e}A = {x : (rad A)·x = 0 = x·(rad A)}`
  (basis via `families/trivial_extension.py::_bimodule_socle`); actions induced on
  the quotient; `dim_M = m − dim soc`. (`Bimodule.mod_socle`.)
- **generic** `from_actions(A, dim_M, left_maps, right_maps)` — one exact matrix
  per arrow per side, expanded (fold over paths) to `Lact`/`Ract`; validated.

### P5. The duality identity (proof, for the Task-6 pin)

`DA = Hom_k(A, k)`. Hom-tensor adjunction:
`Hom_{A^e}(P, Hom_k(A, k)) ≅ Hom_k(A ⊗_{A^e} P, k)` naturally in the `A^e`-module
`P`. Applied degreewise to `P_•`:
`Hom_{A^e}(P_•, DA) ≅ (A ⊗_{A^e} P_•)^*` (the `k`-linear dual complex). Over a
field, cohomology of the dual complex is the dual of homology:
`H^n(C_•^*) ≅ (H_n C_•)^*`. Hence `HH^n(A, DA) ≅ HH_n(A, A)^*`, so
**`dim HH^n(A, D(A)) = dim HH_n(A, A)`** for every `A`, every field, every `n`.
(Dually `dim HH_n(A, DA) = dim HH^n(A, A)`.)

### P6. `HH^0` / `HH_0` closed forms (Task-6 self-cert)

`HH^0(A, M) = M^A = {m ∈ M : a·m = m·a ∀a ∈ A}` — the joint kernel of
`(Lact[j] − Ract[j])` over the generators `j` (arrows + idempotents suffice).
`HH_0(A, M) = M/[A, M]`, `[A, M] = span{a·m − m·a}` = `M / (image of the same
maps)`. For `M = A`: `HH^0 = Z(A)`, `HH_0 = A/[A,A]`.

### P7. Relative `HH_•(A|E)` (Cibils reduced complex, separable `E`)

`Bar_n^E(A, M) = M ⊗_E \bar A_E^{⊗_E n}`, `\bar A_E = A/E`, tensor over
`E = kQ₀`. `E` separable ⇒ each term is `A^e`-projective ⇒ this computes the
absolute `HH_•(A, M)` (P0). The differential is the same alternating sum with the
`E`-relative tensor (only `E`-composable tuples survive). Single-vertex: `E = k`,
identical to P1.

---

## Global constraints

- **Venv is ALWAYS** `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python`.
  Repo root `/Users/marco/Desktop/HomologicalNetworks/quiverlab`. Branch
  `plan-52-hh-coefficients` off `dev`.
- **Every test command** exports `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2` and runs
  `-m pytest -q` from the repo root; ≤2 parallel processes; foreground.
- **Exact only.** The float-ban AST gate `tests/test_no_floats.py` stays green at
  every commit; no float/complex literals in `src/`. All coefficient actions are
  Domain matrices.
- **Loud typed refusals** at every scope boundary (`QuiverlabError` with `hint=`):
  `engine="fast"` + coefficients; general `relative_to`; CS coefficient on a
  presentation-less algebra; coefficient whose algebra ≠ `A`; inconsistent
  bimodule (`check()` fails); `_MAX_COEFF_DIM` allocation guard;
  `DepthLimitError` from the growth guard with the certified range.
- **The `coefficients=None` byte-identity gate is HARD.** Every existing HH
  golden, HH test, and `_runner_goldens.json` entry stays green with **zero**
  re-freezing. Task 1 pins the current `HHTable` attribute set and the hh block
  shape so drift is caught before any edit.
- **Left-to-right composition** everywhere (`a*b` = first `a` then `b`); the CS
  `b·w·a` / `a·w·b` op-twist is stated where used (P2).
- **Plan-32 oracle markers** on every new test; **canonical keys never re-keyed**
  (absent-block pop precedent); every merge updates `docs/verification.md` (new
  oracle rows + recounted class table, `tests/release/test_oracle_classes.py`
  green) + the README line.
- **Conventional commits**, green at every commit; trailer
  `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- **The bank is READ-ONLY** and `src/quiverlab/engine/bimodule.py` is treated as a
  read-only reference oracle (ported hanlab; do not modify it — copy its logic
  into the exact-Domain twin with attribution).

---

### Task 1: Interface freshness gate + branch

**Files:** branch `plan-52-hh-coefficients`; `tests/hochschild/test_coeff_interface_gate.py`.

**Interfaces:** consumes `core.Algebra.hochschild_{co}homology`,
`hochschild/bar.py`, `hochschild/table.py::HHTable`,
`engine/bimodule.py`, `Algebra.nakayama_automorphism`,
`families/trivial_extension.py::_bimodule_socle`. Produces the freshness gate.
STOP on any drift.

- [ ] **Step 1: Branch.** `git checkout -b plan-52-hh-coefficients` (off `dev`).

- [ ] **Step 2: Failing gate** `tests/hochschild/test_coeff_interface_gate.py`
  — pin the REAL surface this plan consumes, so a drift STOPs the plan:

```python
"""Plan-52 interface freshness gate. Pins the HH dispatch, HHTable, bar.py
reduction convention, the ported engine/bimodule reference, and the Nakayama /
bimodule-socle helpers. STOP on any drift."""
import inspect, pytest
import quiverlab as ql
from quiverlab.hochschild import bar, table


def test_hh_signatures_have_no_coefficients_yet():
    for name in ("hochschild_cohomology", "hochschild_homology"):
        params = inspect.signature(getattr(ql.Algebra, name)).parameters
        assert {"top", "max_cells", "engine", "auto_cs", "verbose", "trace"} <= set(params)
        assert "coefficients" not in params          # this task adds it (Task 3)
        assert "relative_to" not in params           # Task 7 adds it


def test_hhtable_attribute_set_frozen():
    t = ql.truncated_polynomial(2, field=ql.GF(7)).hochschild_cohomology(2)
    assert isinstance(t, table.HHTable)
    assert set(vars(t)) == {"dims", "kind", "top", "algebra_repr", "engine", "references"}
    #  ^ adding `coefficients` must default None and NOT appear here for the regular case


def test_bar_is_absolute_no_E_composability_filter():
    # SEMANTIC invariant (not a substring match): the reduced bar enumerates ALL
    # (m-1)^n multi-indices over the reduced basis {1..m-1} (unit index 0 dropped),
    # with NO E-composability filter. For a MULTI-VERTEX algebra an E-relative
    # complex (Task 7) would keep only E-composable tuples -- strictly fewer -- so
    # this exact count fails loudly if someone changes the reduction convention or
    # adds an E-filter to bar.py. (If _abar_tuples' signature has drifted, THIS
    # gate STOPs the plan -- which is the point.)
    A = ql.Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=ql.GF(7))
    m = A.dim                                          # kA2: dim 3 (e1, e2, a)
    for n in (1, 2, 3):
        tuples = list(bar._abar_tuples(m, n))
        assert len(tuples) == (m - 1) ** n             # 2, 4, 8 -- absolute, unfiltered
        assert all(len(t) == n and all(1 <= i < m for i in t) for t in tuples)


def test_ported_gfp_bimodule_reference_present():
    from quiverlab.engine import bimodule as eb
    for n in ("Bimodule", "regular_bimodule", "dual_bimodule", "twisted_bimodule",
              "check_bimodule", "hochschild_homology_with_coefficients",
              "hochschild_cohomology_with_coefficients"):
        assert hasattr(eb, n), f"engine.bimodule.{n} missing (reference oracle drift)"


def test_nakayama_generic_route_and_socle_helpers():
    # DD1b: twisted_by_nakayama uses the GENERIC trace-form route
    # (nakayama_automorphism_generic, A's OWN basis) -- NOT
    # Algebra.nakayama_automorphism(), which over GF(p) returns the ENGINE
    # unit-adapted basis (core/algebra.py:896-897). Pin the generic helper's
    # presence + A-basis shape, and the bimodule-socle helper for mod_socle.
    from quiverlab.invariants.frobenius import nakayama_automorphism_generic
    A = ql.truncated_polynomial(3, field=ql.CC)        # Frobenius => nu exists
    nu = nakayama_automorphism_generic(A)              # m x m, columns = images, A-basis
    assert len(nu) == A.dim and len(nu[0]) == A.dim
    from quiverlab.families.trivial_extension import _bimodule_socle
    assert callable(_bimodule_socle)
```

- [ ] **Step 3: Run.** Expected `5 passed`. Any FAIL ⇒ STOP; reconcile every
  Task-2..10 call site with the drifted surface in one commit, recording it.

- [ ] **Step 4: Commit** `test(hh): Plan-52 interface freshness gate`.

---

### Task 2: exact-`Domain` `Bimodule` container + constructors + validation

**Files:** create `src/quiverlab/hochschild/coefficients.py`;
`tests/hochschild/test_coefficients_bimodule.py`. Export `Bimodule` from
`quiverlab.__init__` (public).

**Interfaces:** consumes `core.Algebra` (`.T, .dim, .unit, .domain,
._basis_vec, .multiply, .quiver`, `.is_symmetric`, `.nakayama_automorphism` —
used only in the DD1b control assertion, never for construction),
`invariants/frobenius.py::nakayama_automorphism_generic` (the A-basis ν used to
build the twist, DD1b), `families/trivial_extension.py::_bimodule_socle`,
`modules/linalg_mod` (kernel / image / rank over Domain). Produces `Bimodule` +
constructors mirroring the GF(p) `engine/bimodule.py` over any Domain, **and the
chain-level `twisted_homology_classes`/`twisted_cohomology_classes` accessors**
(the P54 cross-plan contract). **Names the consumers** (P54
`twisted_by_nakayama` + `twisted_homology_classes`, P74 `twisted(phi=g)`, P72
`from_actions` on `B`).

- [ ] **Step 1: Failing tests** (`oracle_selfcert` — the container's
  compatibility certificate IS the assertion; plus the GF(p) parity with the
  ported reference is `oracle_crossengine`):

```python
"""Exact-Domain Bimodule. check() certifies the A^e-module axioms; the named
constructors reproduce engine/bimodule.py over GF(p) (the ported reference)."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.errors import QuiverlabError

pytestmark = pytest.mark.oracle_selfcert


def _A(field): return ql.truncated_polynomial(3, field=field)   # dual numbers-ish, self-injective


def test_regular_reproduces_algebra_multiplication():
    A = _A(ql.CC)
    M = Bimodule.regular(A)
    assert M.dim_M == A.dim
    M.check()                                    # loud QuiverlabError on any axiom failure
    # Lact/Ract are literally A.T
    assert M.Lact == [[list(row) for row in A.T[j]] for j in range(A.dim)]


def test_check_rejects_noncommuting_actions():
    A = _A(ql.CC)
    M = Bimodule.regular(A)
    M.Ract[1] = M.Lact[1]                         # break left/right commuting (generic tamper)
    with pytest.raises(QuiverlabError):
        M.check()


def test_dual_matches_ported_reference_over_gfp():
    # exact-Domain dual ≡ engine/bimodule.py::dual_bimodule (hanlab port) over
    # GF(p). SINGLE-VERTEX and already unit-adapted (k[x]/(x^3)) => A's basis ==
    # the engine basis, so Lact/Ract are directly comparable. (The twisted case is
    # NOT compared here -- its A-basis nu vs the engine-basis reference would need
    # a change-of-basis; it is validated by the multi-vertex discriminator below +
    # the Task-6 twisted≡dual/duality oracles.)
    from quiverlab.engine import bimodule as eb
    from quiverlab.engine.adapter import to_engine
    A = _A(ql.GF(32003)); eng = to_engine(A.unit_adapted())
    assert Bimodule.dual(A).same_actions_mod_p(eb.dual_bimodule(eng), 32003)


def test_twisted_nakayama_uses_A_basis_multivertex_gfp():
    # DD1b (M3): a 3-vertex NON-symmetric self-injective algebra over GF(p) where
    # unit_adapted() is a GENUINE change of basis and nu != id. twisted_by_nakayama
    # must build Lact/Ract with nu in A'S OWN basis: check() certifies the
    # A^e-module axioms and FAILS if the engine unit-adapted matrix were used.
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=ql.GF(32003))
    assert A.is_symmetric() is False               # nu != identity: nontrivial twist
    Bimodule.twisted_by_nakayama(A).check()        # loud if nu is in the wrong basis
    # CONTROL (adjust to reality): the raw ENGINE-basis matrix must NOT validate --
    # if it does, either this algebra is accidentally unit-adapted (pick another) or
    # check() is too lenient (a check() bug). nu_engine is A.nakayama_automorphism().
    nu_engine = A.nakayama_automorphism()          # engine unit-adapted basis (core/algebra.py:896)
    with pytest.raises(QuiverlabError):
        Bimodule.twisted(A, psi=nu_engine).check()


@pytest.mark.oracle_crossengine
def test_floats_refused_in_from_actions():
    A = ql.Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=ql.CC)
    with pytest.raises(QuiverlabError):
        Bimodule.from_actions(A, 2, left_maps={"a": [[0.5]]}, right_maps={"a": [[0]]})


def test_mod_socle_dim():
    A = _A(ql.CC)                                  # k[x]/(x^3): soc_{A^e}A = <x^2>, dim 1
    M = Bimodule.mod_socle(A)
    M.check()
    assert M.dim_M == A.dim - 1
```

- [ ] **Step 2: Run** — expect `ModuleNotFoundError: quiverlab.hochschild.coefficients`.

- [ ] **Step 3: Implement `coefficients.py`.** `Bimodule(A, dim_M, Lact, Ract,
  name)` (frozen-ish; `check()` verifies: each `Lact[j]`/`Ract[j]` is
  `dim_M×dim_M` over `A.domain`; `L` is an algebra map `Lact[i]·Lact[j] =
  Lact[A.T-product of i,j]`; `R` an anti-map `Ract[i]·Ract[j] = Ract[product]`;
  `Lact` and `Ract` commute; the unit acts as identity on both sides — checked on
  generators = arrows + idempotents, exact, loud). `corner(v,w)` = the
  `M`-basis indices fixed by `Lact[e_v]·(−)·Ract[e_w]`.
  - `regular(A)`: `dim_M=m`, `Lact[j][s] = A.T[j][s]`, `Ract[j][s] = A.T[s][j]`.
  - `dual(A)`: transpose actions per P4 (validate ≡ `engine/bimodule.py::dual_bimodule`
    mod p in the test).
  - `twisted(A, phi=None, psi=None)`: `Lact = A.T ∘ φ`, `Ract = A.T ∘ ψ` per the
    `a·x·b = φ(a)·x·ψ(b)` convention; `None` = identity. `phi`/`psi` are `m×m`
    Domain matrices **in `A`'s own basis** (DD1b).
  - `twisted_by_nakayama(A)`: `psi = invariants.frobenius.nakayama_automorphism_generic(A)`
    — the trace-form route **in `A`'s own basis** (DD1b), **NOT**
    `A.nakayama_automorphism()` (which over GF(p) returns the engine unit-adapted
    basis and is a silent bug against `A.T`-built actions). **(P54 entry point.)**
  - `mod_socle(A)`: quotient by `_bimodule_socle(A, rad, dom)`; induced actions.
  - `from_actions(A, dim_M, left_maps, right_maps)`: per-arrow matrices → fold
    over paths → `Lact`/`Ract`; `_valid_entry` (int / exact string, no float);
    `_MAX_COEFF_DIM` guard. **(P72 `H^*(B,M)`, P74 `{}_gA`, generic no-code.)**
  Add `same_actions_mod_p` and `describe()` (the HHTable provenance string).
  Export `Bimodule` in `__init__` (also re-export `QuiverlabError` there if the
  tests reference `ql.QuiverlabError`). **Do NOT add any
  `nakayama_automorphism_matrix_engine` helper — it does not and need not exist;
  `twisted_by_nakayama` takes the generic A-basis route (DD1b).**

- [ ] **Step 3b: The P54 chain-level accessor (cross-plan contract, W2).** Add
  `twisted_homology_classes(A, M, top)` (and its cohomology twin
  `twisted_cohomology_classes`) returning, per degree `n`: the twisted bar-chain
  **basis** (the `(s, J)` enumeration), the exact Domain **boundary matrices**
  `b_n`, and a **basis of class representatives** (`ker b_n / im b_{n+1}` lifted).
  Implement the CONVENTION BLOCK **verbatim from P54 §3** (Goal/consumer section
  above): twist on the LAST face map only, reduced bar (`\bar A = A/k·1`),
  `M`-index `s` outermost then bar multi-index `J` lexicographic. Add a
  self-cert test `test_twisted_classes_convention` pinning the basis ordering and
  `b_{n-1}∘b_n = 0` on `k[x]/(x^3)` with `{}_1A_ν`, and asserting the class count
  equals `A.hochschild_homology(top, coefficients=M).dims` (the dims surface and
  the chain surface agree). Marked `oracle_selfcert`.

- [ ] **Step 4: Run** — expect PASS. Run `tests/test_no_floats.py`.

- [ ] **Step 5: Commit** `feat(hh): exact-Domain Bimodule + constructors + twisted_homology_classes (P54 chain contract; P54/P72/P74 coefficients)`.

---

### Task 3: coefficient-aware bar collapse + public API kwarg + None byte-gate + M=A / GF(p)-bank oracles

**Files:** modify `src/quiverlab/hochschild/bar.py`, `src/quiverlab/hochschild/table.py`,
`src/quiverlab/core/algebra.py`; test `tests/hochschild/test_coeff_bar.py`.

**Interfaces:** consumes `Bimodule` (Task 2). Produces coefficient-aware
`hochschild_{co}homology_dims(A, top, coefficients=None, ...)` (bar) + the public
`coefficients=` kwarg on `Algebra.hochschild_{co}homology`, with the None
short-circuit.

- [ ] **Step 1: Failing tests:**

```python
"""HH(A,M) via the bar reference. None short-circuit is byte-identical; regular(A)
through the general path ≡ None (M=A oracle); exact path ≡ ported GF(p) reference."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.errors import QuiverlabError


@pytest.mark.oracle_crossengine
def test_regular_coefficient_equals_no_coefficient():
    A = ql.truncated_polynomial(3, field=ql.CC)
    base_h = A.hochschild_homology(5, engine="bar").dims
    base_c = A.hochschild_cohomology(5, engine="bar").dims
    Mreg = Bimodule.regular(A)
    assert A.hochschild_homology(5, engine="bar", coefficients=Mreg).dims == base_h
    assert A.hochschild_cohomology(5, engine="bar", coefficients=Mreg).dims == base_c


@pytest.mark.oracle_crossengine
def test_exact_bar_matches_ported_gfp_reference():
    from quiverlab.engine import bimodule as eb
    from quiverlab.engine.adapter import to_engine
    A = ql.truncated_polynomial(3, field=ql.GF(32003)); eng = to_engine(A.unit_adapted())
    Md = Bimodule.dual(A)
    ours = A.hochschild_cohomology(5, engine="bar", coefficients=Md).dims
    ref = eb.hochschild_cohomology_with_coefficients(eng, eb.dual_bimodule(eng), 5)[32003]
    assert ours == ref


@pytest.mark.oracle_selfcert
def test_none_is_byte_identical_short_circuit():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    t0 = A.hochschild_cohomology(3)
    assert set(vars(t0)) == {"dims", "kind", "top", "algebra_repr", "engine", "references"}
    # coefficients provenance is only surfaced for a non-None coefficient
    t1 = A.hochschild_cohomology(3, coefficients=Bimodule.dual(A))
    assert getattr(t1, "coefficients", None) == "D(A)"


def test_fast_refuses_coefficients():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    with pytest.raises(QuiverlabError):
        A.hochschild_homology(2, engine="fast", coefficients=Bimodule.dual(A))
```

- [ ] **Step 2: Run** — expect FAIL (`TypeError: unexpected keyword 'coefficients'`).

- [ ] **Step 3: Implement.**
  - `bar.py`: give `coboundary_matrix`/`boundary_matrix`/`hochschild_*_dims` a
    `coefficients=None` param. When None, the code path is **untouched**
    (byte-identical). When `M`, substitute per P1 at the four outer term sites,
    set `dim_M = M.dim_M`, `cn = dim_M · (m-1)^n`, and index the coefficient slot
    over `M`'s basis (fold `M.Lact[K[0]]`/`M.Ract[K[n]]` etc.). Interior stays
    `A.T`.
  - `table.py::HHTable`: add `coefficients=None` (last param); store; **only**
    include in `vars`/serialization when non-None (default None ⇒ current
    attribute set unchanged — Task-1 gate).
  - `core/algebra.py`: add `coefficients=None` to both public methods; thread to
    the bar route; set `table.coefficients = coefficients.describe()` when
    non-None; make `_use_fast_engine` return `False` when `coefficients is not
    None`; raise `QuiverlabError` on `engine="fast"` + coefficients; assert
    `coefficients.algebra is self` (or equal presentation).

- [ ] **Step 4: Run** — expect PASS. Then run the **hard byte gate**:
  `... -m pytest tests/hochschild tests/webapp/test_runner_delegation.py -q` —
  every existing golden green, zero re-freeze.

- [ ] **Step 5: Commit** `feat(hh): bar HH(A,M) with coefficients + public coefficients= kwarg (None byte-identical)`.

---

### Task 4: coefficient-aware CS collapse + `engine="cs"`/auto routing + bar ≡ CS

**Files:** modify `src/quiverlab/resolutions_cs/resolution.py` (`matrix`, `_basis`,
`dim_C`), `src/quiverlab/resolutions_cs/homology.py`
(`cs_{co}homology_dims(..., coefficients=None)`), `src/quiverlab/core/algebra.py`
(thread coefficients into `_route_to_cs`/`_cs_depth_fallback`); test
`tests/resolutions_cs/test_coeff_cs.py`.

**Interfaces:** consumes `Bimodule`, `AArith` (extend `corner`/`mul` to `M`).
Produces the CS coefficient route + the standing cross-engine oracle.

- [ ] **Step 1: Failing tests** (`tests/resolutions_cs/` → **deep**):

```python
"""CS HH(A,M) ≡ bar HH(A,M) degreewise over CC and GF(p), every named coefficient.
The resolution of A is coefficient-independent; only the collapse changes."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims

pytestmark = pytest.mark.oracle_crossengine


@pytest.mark.parametrize("field", [ql.CC, ql.GF(7)])
@pytest.mark.parametrize("mk", ["regular", "dual", "twisted_by_nakayama"])
def test_cs_matches_bar_with_coefficients(field, mk):
    A = ql.truncated_polynomial(3, field=field)     # presented, admissible (symmetric)
    M = getattr(Bimodule, mk)(A)
    top = 5
    assert cs_cohomology_dims(A, top, coefficients=M).dims == \
        A.hochschild_cohomology(top, engine="bar", coefficients=M).dims
    assert cs_homology_dims(A, top, coefficients=M).dims == \
        A.hochschild_homology(top, engine="bar", coefficients=M).dims


@pytest.mark.parametrize("mk", ["dual", "twisted_by_nakayama"])
def test_cs_matches_bar_nontrivial_nu(mk):
    # C2: bar≡CS on a NON-symmetric self-injective algebra (nu != id, so the
    # twist path is genuinely exercised, not ≅ regular). NakayamaAlgebra(3,3).
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=ql.GF(7))
    assert A.is_symmetric() is False
    M = getattr(Bimodule, mk)(A)
    assert cs_cohomology_dims(A, 5, coefficients=M).dims == \
        A.hochschild_cohomology(5, engine="bar", coefficients=M).dims


def test_public_engine_cs_routes_coefficients():
    A = ql.truncated_polynomial(3, field=ql.CC)
    M = Bimodule.dual(A)
    assert A.hochschild_cohomology(4, engine="cs", coefficients=M).dims == \
        A.hochschild_cohomology(4, engine="bar", coefficients=M).dims


def test_cs_needs_presentation():
    # structure-constant-only algebra (no _quiver): CS refuses (as today), bar
    # still works. Real ctor is Algebra(domain, T, unit, ...) -- a rebuild WITHOUT
    # a presentation, so resolutions_cs cannot form its reduction system.
    from quiverlab.errors import QuiverlabError
    A = ql.truncated_polynomial(3, field=ql.CC)
    Asc = ql.Algebra(A.domain, A.T, A.unit)          # presentation-less (_quiver=None)
    assert getattr(Asc, "quiver", None) is None      # confirm no presentation
    Md = Bimodule.dual(Asc)
    with pytest.raises(QuiverlabError):
        Asc.hochschild_cohomology(3, engine="cs", coefficients=Md)
    # bar still serves the same coefficient (presentation-free reference route)
    assert Asc.hochschild_cohomology(3, engine="bar", coefficients=Md).dims == \
        A.hochschild_cohomology(3, engine="bar", coefficients=Bimodule.dual(A)).dims
```

- [ ] **Step 2: Run** — expect FAIL (`coefficients` unknown to CS).

- [ ] **Step 3: Implement.**
  - `AArith`: add `corner_M(M, o, t, side)` returning `M`-basis indices of
    `e_t M e_o` / `e_o M e_t`, and route the collapse `mul`s through `M.Lact` /
    `M.Ract` (the paths `a_vec, c_vec` stay `A`-vectors acting on `M`).
  - `resolution.py::matrix`: parameterize by an optional coefficient `M`
    (default None = `A`, byte-identical). `_basis`/`dim_C` use `corner_M` when
    `M` is set. Collapse per P2: homology `M.Lact(b_i)·(M.Ract(a_i)·w)`,
    cohomology `M.Lact(a_i)·(M.Ract(b_i)·w)`.
  - `homology.py`: `cs_{co}homology_dims(A, top, coefficients=None, ...)`.
  - `core/algebra.py`: thread `coefficients` into the CS route and the Plan-34
    depth fallback (`_cs_depth_fallback` passes coefficients through).

- [ ] **Step 4: Run** — expect PASS. Confirm `d∘d=0` / order gates still hold
  (resolution unchanged) via the existing CS battery.

- [ ] **Step 5: Commit** `feat(cs): HH(A,M) collapse with coefficients (any Domain) + engine="cs"/auto routing`.

---

### Task 5: coefficient-aware minimal `A^e` collapse (GF(p), internal) + bar ≡ minimal

**Files:** modify `src/quiverlab/engine/resolutions_minimal.py`
(`_contracted_degree`, `_cohomology_degree`, `_corner_contracted_degree`,
`_corner_cohomology_degree`, `minimal_{co}homology_dims(..., coefficients=None)`,
`_CornerContext` corner blocks); test `tests/engine/test_coeff_minimal.py`.

**Interfaces:** internal (NOT public `engine=`). Consumes `Bimodule`; the
coefficient must be over GF(p) (engine int64). Produces the GF(p) cross-check +
P54's future GF(p) tool. Generalizes the Plan-16 covariance (P3).

- [ ] **Step 1: Failing tests** (`tests/engine/` → **deep**):

```python
"""Minimal A^e HH(A,M) over GF(p) ≡ bar HH(A,M), every named coefficient.
Plan-16 covariance generalized: hom b·w·a on e_w M e_v, coh a·w·b on the SWAPPED
e_v M e_w. Multi-vertex included (the corner path)."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.engine.resolutions_minimal import minimal_homology_dims, minimal_cohomology_dims

pytestmark = pytest.mark.oracle_crossengine


@pytest.mark.parametrize("mk", ["regular", "dual", "twisted_by_nakayama"])
def test_minimal_matches_bar_single_vertex(mk):
    A = ql.truncated_polynomial(3, field=ql.GF(32003))
    M = getattr(Bimodule, mk)(A)
    N = 6
    assert minimal_homology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_homology(N, engine="bar", coefficients=M).dims
    assert minimal_cohomology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_cohomology(N, engine="bar", coefficients=M).dims


def test_minimal_matches_bar_multivertex():
    # kA3 with a*b=0: exercises the corner path + SWAPPED-tag coh block
    A = ql.Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=ql.GF(32003))
    M = Bimodule.dual(A)
    N = 5
    assert minimal_cohomology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_cohomology(N, engine="bar", coefficients=M).dims


def test_minimal_matches_bar_nontrivial_nu_multivertex():
    # C2: multi-vertex NON-symmetric self-injective over GF(p) -- exercises the
    # twist path (nu != id) AND the corner/SWAPPED-tag block simultaneously, and
    # is the DD1b basis discriminator on the minimal engine.
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=ql.GF(32003))
    assert A.is_symmetric() is False
    M = Bimodule.twisted_by_nakayama(A)
    assert minimal_cohomology_dims(A, 5, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_cohomology(5, engine="bar", coefficients=M).dims
```

- [ ] **Step 2: Run** — expect FAIL.

- [ ] **Step 3: Implement per P3.** Coefficient block: `cornerA[(i,j)] = e_j A e_i`
  → `e_j M e_i` (a `Bimodule` corner); local "full copy of `A`" → "full copy of
  `M`". In each of the four collapse functions swap the two structure-constant
  multiplications for `M.Lact`/`M.Ract` at the sites named in P3, leaving the
  `uu ↔ vv` covariance and the SWAPPED-tag block accounting untouched.
  `minimal_{co}homology_dims` gain `coefficients=None` (None = A, byte-identical).
  Loud refusal if the coefficient's Domain is not GF(p) (engine int64).

- [ ] **Step 4: Run** — expect PASS. Run pure/numba parity spot-check on one case
  (`QUIVERLAB_NO_NUMBA=1`).

- [ ] **Step 5: Commit** `feat(engine): minimal A^e HH(A,M) over GF(p) (Plan-16 covariance generalized) + bar≡minimal oracle`.

---

### Task 6: `HH^0=M^A`, `HH_0=M/[A,M]`, the duality identity, the symmetric web

**Files:** test `tests/hochschild/test_coeff_identities.py` (no `src/` change if
the closed forms are computed in-test from the actions; add a small
`coefficients.py` helper `Bimodule.invariants()` / `.coinvariants()` if convenient).

**Interfaces:** the internal-identity oracles that carry certification while the
literature pins are BLOCKED.

- [ ] **Step 1: Failing/anchoring tests:**

```python
"""Coefficient identities: HH^0=M^A, HH_0=M/[A,M], the exact duality
dim HH^n(A,DA)=dim HH_n(A,A), and the symmetric-algebra consistency web."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule


@pytest.mark.oracle_selfcert
def test_degree0_closed_forms():
    A = ql.truncated_polynomial(3, field=ql.CC)
    M = Bimodule.dual(A)
    c0 = A.hochschild_cohomology(0, coefficients=M).dims[0]
    h0 = A.hochschild_homology(0, coefficients=M).dims[0]
    assert c0 == M.invariants_dim()          # dim M^A = joint ker(Lact[j]-Ract[j])
    assert h0 == M.coinvariants_dim()         # dim M/[A,M]


@pytest.mark.oracle_selfcert     # W4: a theory identity on one engine's own output
def test_duality_identity_smoke_symmetric():
    # SMOKE ONLY: truncated_polynomial(3) is symmetric, so D(A) ≅ A and the
    # identity is (nearly) tautological -- it exercises the plumbing, not content.
    A = ql.truncated_polynomial(3, field=ql.CC)
    assert A.is_symmetric() is True
    DA = Bimodule.dual(A)
    assert A.hochschild_cohomology(6, engine="bar", coefficients=DA).dims == \
        A.hochschild_homology(6, engine="bar").dims


@pytest.mark.oracle_selfcert     # W4
@pytest.mark.parametrize("field", [ql.CC, ql.GF(32003)])
def test_duality_identity_nonsymmetric_discriminator(field):
    # C3: the CONTENT case -- a NON-symmetric self-injective algebra where
    # D(A) ≇ A, so dim HH^n(A,DA)=dim HH_n(A,A) has real force (a wrong dual or a
    # wrong nu could not accidentally satisfy it). NakayamaAlgebra(3,3,cyclic).
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=field)
    assert A.is_symmetric() is False
    DA = Bimodule.dual(A)
    assert A.hochschild_cohomology(6, engine="bar", coefficients=DA).dims == \
        A.hochschild_homology(6, engine="bar").dims


@pytest.mark.oracle_selfcert
def test_symmetric_algebra_web():
    # A GENUINELY symmetric algebra (verified in-snippet): nu is inner => _1A_nu ≅ A,
    # D(A) ≅ A; all four HH agree, and dim HH^n = dim HH_n follows. NOTE: on a
    # symmetric algebra the twisted/dual tests are VACUOUS (twist ≅ regular); the
    # DISCRIMINATING nontrivial-nu coverage lives in
    # test_duality_identity_nonsymmetric_discriminator + Task-2's multi-vertex
    # twisted_by_nakayama gate.
    A = ql.truncated_polynomial(3, field=ql.CC)                # symmetric (verified below)
    assert A.is_symmetric() is True
    base = A.hochschild_cohomology(5).dims
    assert A.hochschild_cohomology(5, coefficients=Bimodule.dual(A)).dims == base
    assert A.hochschild_cohomology(5, coefficients=Bimodule.twisted_by_nakayama(A)).dims == base
    assert A.hochschild_homology(5).dims == base


@pytest.mark.oracle_selfcert
def test_twisted_equals_dual_nontrivial_nu():
    # C2/C3 discriminator: for self-injective A, D(A) ≅ _1A_nu, so twisted and dual
    # coefficients give IDENTICAL HH dims -- non-vacuously on NON-symmetric A
    # (nu != id). Fails if twisted_by_nakayama built nu in the wrong basis (DD1b).
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=ql.GF(32003))
    assert A.is_symmetric() is False
    assert A.hochschild_cohomology(5, coefficients=Bimodule.twisted_by_nakayama(A)).dims == \
        A.hochschild_cohomology(5, coefficients=Bimodule.dual(A)).dims
```

- [ ] **Step 2: Run** — expect FAIL (helpers absent) then implement
  `invariants_dim`/`coinvariants_dim` on `Bimodule` (exact kernel/quotient rank
  over Domain via `modules/linalg_mod`). **In-snippet symmetry is asserted with
  `A.is_symmetric()`**, so a mis-chosen family fails loudly rather than passing
  vacuously; `truncated_polynomial(3)` is symmetric and
  `NakayamaAlgebra(3,3,cyclic)` is NOT (both verified by the critic), giving one
  genuine symmetric web + one genuine nontrivial-ν discriminator. The Plan-31
  presented `T(A)` remains a symmetric fallback if a wider symmetric example is
  wanted.

- [ ] **Step 3–5:** Run → PASS → commit `test(hh): HH^0=M^A, HH_0=M/[A,M], duality (symmetric smoke + non-symmetric discriminator), twisted≡dual nontrivial-nu, symmetric web`.

---

### Task 7: relative `HH_•(A|kQ₀, M)` + separability oracle + general-`B` refusal

**Files:** modify `src/quiverlab/hochschild/bar.py` (or a sibling
`hochschild/relative.py`) for the `E`-relative reduced complex;
`src/quiverlab/core/algebra.py` (`relative_to=None` kwarg on both methods); test
`tests/hochschild/test_relative.py`.

**Interfaces:** consumes the Task-3 coefficient collapse (reused with the
`E`-relative tensor `\bar A_E = A/E`). Produces `HH_•(A|kQ₀, M)` and the free
oracle. Substrate note for P72/P73.

- [ ] **Step 1: Failing tests:**

```python
"""Relative HH over the vertex subalgebra E=kQ0. Separable E => relative ==
absolute degreewise (free strong oracle); byte-identical single-vertex; general B
refused loudly."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.errors import QuiverlabError

pytestmark = pytest.mark.oracle_crossengine


def test_relative_equals_absolute_multivertex():
    A = ql.Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=ql.CC)
    for k in range(6):
        assert A.hochschild_homology(k, relative_to="vertices").dims == \
            A.hochschild_homology(k).dims                         # separability


@pytest.mark.oracle_selfcert
def test_relative_byte_identical_single_vertex():
    A = ql.truncated_polynomial(3, field=ql.GF(7))                # E = k
    assert A.hochschild_cohomology(5, relative_to="vertices").dims == \
        A.hochschild_cohomology(5).dims


def test_relative_with_coefficients():
    A = ql.Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=ql.CC)
    M = Bimodule.dual(A)
    assert A.hochschild_homology(4, relative_to="vertices", coefficients=M).dims == \
        A.hochschild_homology(4, coefficients=M).dims


def test_general_B_refused():
    A = ql.truncated_polynomial(3, field=ql.CC)
    with pytest.raises(QuiverlabError):
        A.hochschild_homology(3, relative_to="some-subalgebra")
```

- [ ] **Step 2: Run** — expect FAIL.

- [ ] **Step 3: Implement.** Build `Bar_n^E` over `\bar A_E = A/E` (drop **all**
  vertex idempotents, keep length-≥1 paths), tensor over `E` (only
  `E`-composable tuples survive — the `t(a_i)=s(a_{i+1})` filter). Reuse the
  Task-3 coefficient collapse. `relative_to="vertices"` routes here;
  `relative_to=None` unchanged; anything else → `QuiverlabError` with the
  general-`B` follow-up hint. Cite `cibils_*` (relative theory).

- [ ] **Step 4: Run** — expect PASS.

- [ ] **Step 5: Commit** `feat(hh): relative HH(A|kQ0,M) via the E-relative reduced complex (separability oracle; general B refused)`.

---

### Task 8: literature oracle (BLOCKED-until-transcribed) + QPA honest probe

**Files:** test `tests/hochschild/test_coeff_literature.py`,
`tests/qpa/test_coeff_qpa.py`; modify `src/quiverlab/qpa/scripts.py`/`crosscheck.py`
only if the probe finds a verb (expected: none).

- [ ] **Step 1: Transcribe.** Read the PDF/TeX of **arXiv:1811.02211** and pick
  **one** concrete small gentle example with an explicit `dim HH¹(A, M)` for a
  named coefficient `M` (quiver + relations + `M` + integer + equation number);
  transcribe it verbatim into the test docstring. Likewise **arXiv:2411.03080**
  for a **radical-square-zero** relative-`HH¹` value (Task 7 substrate). **Until
  transcribed, the pin stays `xfail(strict=False, reason="BLOCKED: transcribe
  1811.02211 Ex. …")`** — never a fabricated number.

```python
"""Gentle HH^1-with-coefficients (arXiv:1811.02211) + relative HH^1
(arXiv:2411.03080). BLOCKED-until-transcribed: fill the value + eq. number from
the paper before flipping the xfail to a hard assert."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule

pytestmark = pytest.mark.oracle_literature


@pytest.mark.xfail(reason="BLOCKED: transcribe the 1811.02211 worked example", strict=False)
def test_gentle_hh1_with_coefficients():
    A = ...  # the paper's quiver/relations
    M = ...  # the paper's coefficient bimodule
    assert A.hochschild_cohomology(1, coefficients=M).dims[1] == ...  # paper value + eq #
```

- [ ] **Step 2: QPA probe** (`tests/qpa/`, `-m qpa`): `NamesGVars()` sweep for a
  Hochschild-with-coefficients / relative-HH verb (Plan-35 precedent). Expected:
  none. The test **skips honestly and FAILS if a matching verb ever appears**
  (then a real crosscheck is owed). Where QPA *can* compute ordinary HH of a
  presented `T(A)` (symmetric), cross-check `HH^n(A, D(A))` (= ordinary `HH_n`)
  against QPA's ordinary HH as an indirect anchor if convenient.

- [ ] **Step 3–5:** Run → PASS/xfail/skip as stated → commit
  `test(hh,qpa): coefficient literature pins (BLOCKED-until-transcribed) + QPA honest probe`.

---

### Task 9: GUI/webapp coefficient picker (builtin forms) end-to-end

**Files:** `src/quiverlab/hpc/spec.py` (`parse_request`, `ComputeRequest`,
`_dispatch` hh branch), `webapp/server/schema.py` (`CoefficientSpec`,
`ComputeRequest` field + `model_dump` pop-list + the module-gate `== 2`→`< 2` fix
+ `_schema_known` v1/v2/v3), `webapp/server/estimator.py` (`sizing_dim`
coefficient-dimension term, C4), `docs/gui/runner.py` (twin dispatch +
`ETA_MODEL`), `docs/gui/gui.js` == `webapp/static/gui/gui.js` (picker +
`renderBlock`), `webapp/static/app.js` (webapp renderer),
`webapp/server/i18n/{en,es,fr,zh}.json`, `src/quiverlab/trace/results_html.py`
(provenance line), `tests/webapp/_runner_goldens.json` +
`test_runner_delegation.py` (ONE documented fixture),
`docs/plans/DEEPER-ENGINES-BACKLOG.md` (the explicit-editor ledger entry). Tests:
`tests/webapp/test_coefficients_p52.py`.

**Interfaces:** `hh_cohomology`/`hh_homology` gain an optional `coefficients`
block; the builtin forms only in the GUI (explicit accepted server-side). Schema
`>= 3` gates the block; **conditional versioning** — a coefficients-less request
keeps sending schema 2 (or 1), so absent block canonicalizes byte-identically
(pop-list). The estimator sizes on the coefficient dimension (C4). Both runners
byte-identical.

- [ ] **Step 1: Failing cross-runner test** (unmarked, extras-gated dir; copy the
  `test_module_blocks_m0729.py` runner-pair fixture):

```python
def test_coefficient_block_shape(tmp_path):
    # request: compute ["hh_cohomology:0..4"], coefficients {"builtin":{"kind":"dual"}}, schema 3
    from quiverlab.hpc.spec import run_spec
    out = run_spec(_req_dual(), tmp_path)
    blk = _find_block(out, "HH^")
    assert blk["coefficients"] == "D(A)"
    assert blk["dims"] == _expected_dims_dual()

def test_twin_parity(tmp_path):
    # same request through docs/gui/runner.py; json.dumps(sort_keys=True) equal on the block

def test_absent_coefficients_keeps_canonical_key():
    from webapp.server.schema import ComputeRequest
    from webapp.server.cache import canonical_key
    r_old = ComputeRequest(**_req_plain_v2())      # no coefficients, schema_version 2
    # conditional versioning: a coefficients-less request stays schema 2, its
    # model_dump drops the absent block, and its key equals the pre-Plan-52 key.
    assert canonical_key(r_old.model_dump(by_alias=True), "x") == _frozen_plain_key()


def test_large_explicit_coefficient_does_not_route_instant():
    # C4: an explicit from_actions coefficient of large dim over a SMALL algebra
    # must size the job (bar cost is quadratic in dim_M) -> NOT the instant tier.
    from webapp.server.estimator import sizing_dim
    from webapp.server.schema import ComputeRequest
    small_algebra_dim = 3
    req = ComputeRequest(**_req_big_explicit_coeff(dim_M=60))   # schema 3
    assert sizing_dim(small_algebra_dim, req) >= 60             # coefficient dominates
```

- [ ] **Step 2: Run** — expect FAIL.

- [ ] **Step 3: Implement.**
  - `schema.py` (C5 — the current gate is unsatisfiable for a combined request):
    - `CoefficientSpec` (`{"builtin": {"kind":
      "regular"|"dual"|"twisted_nakayama"|"quotient_socle"}}` **or** the explicit
      `{dim, left_maps, right_maps}` accepted server-side but with **no GUI
      editor**); a `ComputeRequest.coefficients` field requiring
      `schema_version >= 3`; **add `"coefficients"` to the `model_dump`
      absent-block pop-list**; a validator gating it to hh kinds.
    - **change the module gate** (`schema.py:282-283`) from `self.schema_version
      != 2` to `self.schema_version < 2` (module/`ext_target`/`tor_target` blocks
      accept `schema_version >= 2`) — so a request may legally carry both a
      module and a coefficients block.
    - update `_schema_known` / the known-version set + its stale message to speak
      **v1 / v2 / v3**.
    - **CANONICAL-KEY RULE (verbatim, binding):** the GUI and every client send
      `schema_version = 3` **only when a `coefficients` block is present**; a
      request without coefficients keeps sending `schema_version = 2` (or 1), so
      its `model_dump` and canonical key are byte-identical to today. Clients must
      not gratuitously bump the version. Stdlib twin in `spec.py::parse_request`
      (accept schema 3) + `_parse_coefficients`.
  - `estimator.py` (C4): extend `sizing_dim` (`estimator.py:108-119`) with a
    `_coefficient_dim(req)` term inside its `max(...)` — an explicit
    `from_actions` coefficient of dim `D` over a small algebra drives a bar cost
    **quadratic in `D`** (`cn = D·(m-1)^n`), so it must size the job or it would
    mis-classify as instant. Mirror the `_module_dim` precedent (builtin
    coefficients size at `dim A`; explicit at their declared `dim`).
  - `spec.py::_dispatch` hh branch: build the `Bimodule` from the block
    (`Bimodule.regular/dual/twisted_by_nakayama/mod_socle` or `from_actions`),
    pass `coefficients=` to the method, add `block["coefficients"] =
    M.describe()`, append the citation keys.
  - `runner.py`: mirror byte-for-byte; add an `ETA_MODEL["scalars"]` note (hh
    with coefficients ≈ hh cost).
  - `gui.js` (both copies) + `app.js`: a **coefficient picklist** on the
    Hochschild kinds (regular default = no block emitted; dual / `_1A_ν` /
    `A/soc`); `renderBlock` prints the `coefficients:` provenance line.
  - i18n ×4: `inv.*`, `block.hh_*.coefficients`, `pick.kind.*` labels.
  - `results_html.py`: a "coefficients: D(A)" provenance line in the report.
  - ONE golden `hh_cohomology_dual_kA2` in `_runner_goldens.json`, documented in
    `test_runner_delegation.py`'s change-log; run the delegation test BEFORE
    adding it to confirm existing entries stay byte-identical.
  - `DEEPER-ENGINES-BACKLOG.md`: ledger entry "P52 explicit two-sided bimodule
    matrix editor deferred to P80 (Plan-26 editor is one-sided)."

- [ ] **Step 4: Run gates:** `... -m pytest tests/webapp/test_coefficients_p52.py
  tests/webapp/test_runner_delegation.py tests/hpc -q`; `node --check` both
  `gui.js`; `tests/webapp/test_js_parses.py`.

- [ ] **Step 5: Commit** `feat(gui,webapp,hpc): coefficient picker (builtin bimodules) on the Hochschild kinds; explicit editor ledgered`.

---

### Task 10: verification page, README, citations, suite gate

**Files:** `docs/verification.md`, `README.md`,
`src/quiverlab/citations/registry.py` + `references.bib`; existing release gates.

- [ ] **Step 1: Citations.** Add `_r(...)` + `references.bib` entries (verified
  BibTeX; a registry key with no `references.bib` entry makes `bibtex()` raise):
  - `chaparro_schroll_solotar` → `@article{ChaparroSchrollSolotar2020, …
    J. Algebra 558 (2020) 293–326, note={arXiv:1811.02211}}`.
  - `lindell_rubio_relative` → `@article{LindellRubio2024, … note={arXiv:2411.03080}}`.
  - reuse a `cibils_*` key for the relative reduced complex, or add
    `cibils_relative` if a primary source is BibTeX-verifiable.
  Append the keys to the hh-coefficient handler's `keys` list (both runners).
- [ ] **Step 2: Verification page.** Extend the hh subsystem row (bimodule
  coefficients + relative HH), and add oracle rows: `oracle_crossengine` (M=A ≡
  ordinary; bar ≡ CS ≡ minimal incl. nontrivial-ν; exact ≡ ported GF(p) bank;
  relative ≡ absolute); `oracle_selfcert` (`HH^0=M^A`, `HH_0=M/[A,M]`, the proven
  **duality** `dim HH^n(A,DA)`=dim `HH_n(A,A)` with its non-symmetric
  discriminator, twisted ≡ dual on nontrivial ν, symmetric web, single-vertex
  relative byte-identity, the P54 chain-contract self-check); `oracle_literature`
  (1811.02211 gentle-`HH¹` and the **in-scope vertex-relative** 2411.03080, both
  BLOCKED-until-transcribed).
  **Honest-scope entries:** (a) `engine="fast"` refuses coefficients (hard-wired
  regular bimodule); (b) **cyclic homology / Connes `B` with coefficients is out
  of v1 scope** — `connes_B_matrix` assumes the coefficient is `A` (rotates the
  unit into bar slots); the coefficient cyclic theory is P54's concern, refused
  here; (c) relative HH is `B = kQ₀` only, general `B` refused (the interesting
  non-separable case is the recorded follow-up); (d) CS coefficients need a
  presentation (structure-constant algebras → bar only); (e) QPA has no
  Hochschild-with-coefficients surface. Recount the class table
  (`tests/release/test_oracle_classes.py` drives the numbers).
- [ ] **Step 3: README.** One line: "Hochschild (co)homology with **arbitrary
  bimodule coefficients** (`D(A)`, twisted `{}_1A_ν`, `A/soc`, any no-code
  bimodule) and **relative HH over the vertices** — `coefficients=` on the
  Hochschild kinds."
- [ ] **Step 4: Full gate:** `... -m pytest tests/hochschild tests/resolutions_cs
  tests/engine -q` (deep touched dirs); `... -m pytest -q -m fast`;
  `... -m pytest tests/qpa -q -m qpa`; `... -m pytest tests/release -q`;
  `tests/citations/test_bib_structure.py`. All green.
- [ ] **Step 5: Commit** `docs(verification): Plan-52 coefficient + relative-HH oracle rows + honest scope + recounted classes`.

---

## Oracle table (consolidated)

| Oracle | Statement | Class | Where |
|---|---|---|---|
| M = A ≡ ordinary | `regular(A)` through the general path = `coefficients=None` dims, degreewise | `oracle_crossengine` | T3 |
| bar ≡ CS | equal dims, CC + GF(p), every named `M`; nontrivial-ν on `NakayamaAlgebra(3,3,cyclic)` | `oracle_crossengine` | T4 |
| bar ≡ minimal | equal dims over GF(p), incl. multi-vertex corner + nontrivial-ν | `oracle_crossengine` | T5 |
| exact ≡ ported GF(p) bank | exact-Domain bar = `engine/bimodule.py` over GF(p), **degreewise** (dim vectors) | `oracle_crossengine` | T2, T3 |
| `HH^0 = M^A` | dim = joint ker`(Lact−Ract)` | `oracle_selfcert` | T6 |
| `HH_0 = M/[A,M]` | dim = coinvariants | `oracle_selfcert` | T6 |
| duality | `dim HH^n(A,DA) = dim HH_n(A,A)` (proven, P5); non-symmetric discriminator (`NakayamaAlgebra(3,3,cyclic)`, `D(A)≇A`) | `oracle_selfcert` | T6 |
| twisted ≡ dual (nontrivial ν) | `HH^•(A,{}_1A_ν) = HH^•(A,DA)` on non-symmetric self-injective `A` | `oracle_selfcert` | T6 |
| symmetric web | symmetric `A` (verified `is_symmetric`): `HH^•(A,DA)=HH^•(A,{}_1A_ν)=HH^•(A,A)`, `dim HH^n=dim HH_n` | `oracle_selfcert` | T6 |
| relative ≡ absolute | `HH_•(A|kQ₀,M) = HH_•(A,M)` degreewise (separability, P7) | `oracle_crossengine` | T7 |
| relative single-vertex | byte-identical to absolute | `oracle_selfcert` | T7 |
| gentle `HH¹`-coeff | 1811.02211 worked value | `oracle_literature` (**BLOCKED**) | T8 |
| relative `HH¹` (vertex-relative, **in-scope**) | 2411.03080 rad²=0 value | `oracle_literature` (**BLOCKED**) | T8 |
| P54 chain contract | class count = dims + `b∘b=0`; convention block mirrors P54 §3 | `oracle_selfcert` | T2 |
| None byte-identity | every existing HH golden / runner golden unchanged | hard gate | T3, T9 |
| coefficient-less canonical key | absent block ⇒ key unchanged (conditional versioning) | hard gate | T9 |

## Honest scope (verification-page entries)

- **`engine="fast"` refuses coefficients** — the GF(p) bar-basis accelerator
  carries only `(m, T, unit)`; no coefficient object exists. Loud
  `QuiverlabError` pointing at `bar`/`cs`.
- **Cyclic homology / Connes `B` with coefficients is out of v1 scope** — Connes
  `B` (`hochschild/cyclic.py`) assumes the coefficient is `A` itself; general `M`
  is P54's BV concern.
- **Relative HH is `B = kQ₀` only** — the **vertex-relative (separable) case,
  which is exactly the Lindell–Rubio y Degrassi setting**, is IN scope; a general
  subalgebra `B` is refused, and the **non-separable / non-vertex** case
  (relative ≠ absolute) is the recorded follow-up (Cibils relative theory).
- **CS coefficients need a presentation** — presentation-less structure-constant
  algebras compute coefficients via `bar` only (CS refuses, as today).
- **QPA has no HH-with-coefficients / relative-HH surface** — the covering
  oracles are the internal identities + the ported GF(p) bank + the literature
  pins; the QPA probe FAILS if that ever changes.
- **Literature pins BLOCKED-until-transcribed** — the gentle-`HH¹` and
  relative-`HH¹` values from 1811.02211 / 2411.03080 are `xfail(strict=False)`
  until a concrete example is transcribed verbatim (value + equation number).

## Acceptance (Plan-52 definition of done)

1. `Bimodule` public, exact over any Domain, `check()`-certified; constructors
   `regular`/`dual`/`twisted`/`twisted_by_nakayama`/`mod_socle`/`from_actions`,
   matching `engine/bimodule.py` over GF(p). The P54/P72/P74 entry points exist
   and are named in the docstrings.
2. `Algebra.hochschild_{co}homology(top, coefficients=M, engine=…)` on **bar**
   (any Domain) and **cs** (presented); **minimal** internal over GF(p); **fast**
   refuses. `coefficients=None` is **byte-identical** — every existing golden and
   runner golden green with zero re-freezing.
3. The cross-engine web (M=A ≡ ordinary; bar ≡ CS ≡ minimal; exact ≡ ported GF(p)
   bank) green; `HH^0=M^A`, `HH_0=M/[A,M]`, the proven duality
   `dim HH^n(A,DA)=dim HH_n(A,A)`, and the symmetric web all pinned.
4. Relative `HH_•(A|kQ₀, M)` computes via the `E`-relative reduced complex;
   relative ≡ absolute degreewise (separability), byte-identical single-vertex;
   general `B` refused loudly.
5. Literature pins present as BLOCKED-until-transcribed `xfail`s (never
   fabricated); QPA probe skips honestly and FAILS on a future verb.
6. `coefficients` picker (builtin bimodules) clickable end-to-end (GUI canvas →
   block → report) in EN+ES+FR+ZH, both runners byte-identical, ONE documented
   golden, schema `>= 3`, coefficients-less requests byte-stable in the cache
   key; explicit-editor deferral ledgered.
7. `docs/verification.md` recounted (class table gated), README line added,
   citations resolve (`test_bib_structure.py` green); deep (touched dirs) + fast
   + qpa + release suites green.
