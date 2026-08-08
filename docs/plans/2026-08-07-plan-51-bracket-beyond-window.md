# Plan 51 — the Gerstenhaber bracket beyond the bar window (native homotopy liftings)

> **For agentic workers:** REQUIRED SUB-SKILL — use superpowers:subagent-driven-development
> (or superpowers:executing-plans) to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking. Do not touch `src/`/`tests/` until the
> matching task's failing test is written and run.

**Date:** 2026-08-07 · **Branch:** `plan-51-bracket-beyond-window` (off `dev`) ·
**Status:** spec+plan approved, awaiting implementation · **Record:** R1 (tier α,
Wave 1) of the computability-expansion program · **Size:** L ·
**Deps:** `resolutions_cs/diagonal.py` (Plan 20, shipped), `hochschild/products.py`
(Plan 35, shipped). Independent of the other Wave-1 plans.

**Change-log.** 2026-08-07 adversarial review: 3 majors + 7 minors applied
(coefficient matrix respecified off `d_terms` — a NEW per-corner `D_corner` matrix,
NOT `res.matrix`; consistency contract settled full-system / Lemma-2, phantom
inconsistent-ψ edge retracted; anchor ordering fixed — B.2 uses the existing
transport-only `bracket_of_cs_classes(u,v)`, engine= comparison moved to Task C;
QCI(2,2) named the odd-exponent sign discriminator; + cocycle guard, honest
presentation-less scope, version-scoped key claim, F.2 relabel, conservative build
margins).

**Goal.** Lift the Gerstenhaber bracket `[-,-]: HH^p ⊗ HH^q → HH^{p+q-1}` off the
bar comparison window. Today the bracket is served ONLY over GF(p) and ONLY inside
the bar-comparison window (`Comparison.bracket_of_cs_classes` transports to the bar,
caps at degree `2n+1`; the table route `Algebra.gerstenhaber_brackets` refuses `cs`
and off-GF(p) outright — `core/algebra.py:826-849`). This plan makes the bracket
**native on the Chouhy–Solotar resolution over any exact Domain, at any degree**,
via the Negron–Witherspoon / Volkov **homotopy-lifting** construction — the same
finite per-degree linear solve + `reduce_mod_nullspace` canonicalization that already
builds the Plan-20 diagonal Δ. In-window, over GF(p), the transported bracket is left
**byte-unchanged** and becomes the cross-engine sign anchor. This completes the
Tamarkin–Tsygan calculus surface (cup + cap already went native in Plans 20/21; the
bracket was the last window-bounded operation) and retires the now-outdated
"needs the CS brace/circle machinery" refusal in `comparison.py`.

**Architecture.** Two new modules under the top-level CS package, plus a thin table
builder and one class-level dispatch method — no new engine and no new mathematics;
ONE new assembled linear-algebra primitive (a single-complex corner differential
matrix built off `res.d_terms`, the sibling of the diagonal's `tensor_matrix`), and
the SAME degreewise `solve` + `reduce_mod_nullspace` lift-solve pattern the diagonal
uses, everything exact over the shipped Δ:

- `src/quiverlab/resolutions_cs/homotopy_lifting.py` — the homotopy-lifting tower
  `ψ_η: P → P[1-n]` for a CS cocycle η, built degreewise by the same lift-solve
  *pattern* as `diagonal.py::TensorComplex.diagonal` (coefficient matrix = a NEW
  per-corner **resolution-differential** matrix `d_P: P_{k-n+1} → P_{k-n}` assembled
  from `res.d_terms` — NOT `res.matrix`, which is the Hom-dual coboundary; RHS = the
  diagonal half-collapse `(η⊗1 − 1⊗η)Δ` plus the lower-degree ψ term; `solve` +
  `reduce_mod_nullspace` → canonical, byte-reproducible).
- `src/quiverlab/resolutions_cs/bracket.py` — `native_bracket(res, f_vec, p, g_vec, q)`
  assembling the bracket cochain directly from two homotopy liftings (Oke Thm 3.5),
  the sibling of `cup.py::native_cup` / `cap.py::native_cap`.
- `src/quiverlab/resolutions_cs/products.py` — `cs_bracket_tables(A, top, max_cells)`
  (the bracket sibling of `cs_product_tables`; reuses `_pairs`, `_class_coords`,
  `_capture_cs`, `HHProducts(kind="bracket")`).
- `src/quiverlab/resolutions_cs/comparison.py` — new
  `Comparison.bracket_of_cs_classes(u, v, engine="auto"/"native"/"transport")`,
  mirroring `cup_of_cs_classes`/`cap_of_cs_classes` verbatim (native past window,
  transported in-window byte-unchanged).
- `src/quiverlab/core/algebra.py` — `_product_dispatch` lifts the bracket refusals
  (a `cs` route + the DepthLimitError auto-fallback, previously excluded for
  `kind == "bracket"`).

**Tech stack.** Exact arithmetic through `res.dom` / `AArith` (the CS resolution's
Domain), sympy/int-mod-p linear algebra via `fields.linalg.solve` +
`reduce_mod_nullspace`. **No floats in `src/`** (the AST gate scans these files).
No numba, no engine imports on the CS route (Domain-generic).

---

## 1. Record (verbatim, R1)

Copied verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md §1`
(the durable adjudicated layer; any correction discovered during implementation is
folded back there with a dated note, per the metaplan §2):

> **R1 — Bracket beyond the bar window (homotopy liftings).** Object:
> `[-,-]: HH^p×HH^q → HH^{p+q-1}` natively on the CS/minimal/Bardzell resolution via
> homotopy-lifting maps ψ solving a finite linear system (the `_d_general`-style
> solve + `reduce_mod_nullspace` canonicalization). Refs: Volkov arXiv:1610.05741
> (Proc. Edinb. Math. Soc.); Oke arXiv:2103.12331; Witherspoon GSM 204;
> Negron–Witherspoon A∞-coderivations. Scope: needs the diagonal/comparison
> (shipped); new = past-window + any exact field. Oracles: in-window ≡ transported
> bracket (cross-engine); antisymmetry+Jacobi (self-cert); the 2103.12331
> Koszul-quiver tables; k[x]/(x^n) zero-entries. Size L. Deps:
> resolutions_cs/diagonal.py, hochschild/products.py.

Cross-cluster headline (research doc §0, H1): *"The Gerstenhaber bracket past the bar
window is the single most strategically valuable item. Two clusters independently
converged on it… quiverlab ships BOTH [a Plan-14 comparison and a Plan-20 diagonal Δ
on CS]. Scope: the new content is past-window and off-GF(p); in-window transported
bracket is the cross-engine anchor."*

Metaplan §5 card (P51): *"`[-,-]` natively on CS/minimal/Bardzell via homotopy
liftings … the finite linear solve + `reduce_mod_nullspace` canonicalization, same
pattern as the Plan-20 Δ. New content: past-window + any exact Domain … GUI:
`gerstenhaber_brackets` gains the past-window/off-GF(p) route + provenance gloss."*

---

## 2. Reference re-verification (mandatory, done at spec time 2026-08-07)

Per the metaplan standing rule (§1.5). Verified via arXiv + Cambridge Core + the
authors' pages; **zero discrepancies with the record's attributions**, and the two
primary refs are ALREADY in the citation registry (`registry.py:39-48`
`bracket_liftings`=NegronWitherspoon2016, `bracket_liftings_volkov`=Volkov2019).

| Ref | Status | What it gives | Used as |
|---|---|---|---|
| **Volkov, arXiv:1610.05741** — *Gerstenhaber bracket on the Hochschild cohomology via an arbitrary resolution* | CONFIRMED. Accepted **Proc. Edinburgh Math. Soc.** (submitted 2016-10-18, rev. 2018-11-15). | **Theorem 4**: `[f,g]_{φ,Δ} = (-1)^m f φ_g + (-1)^{m(n-1)} g φ_f`; **Definition 4** homotopy lifting `d φ_f = (f⊗1_P − 1_P⊗f)Δ_P`; **Lemma 2**: `(f⊗1 − 1⊗f)` is null-homotopic for a cocycle `f` (⇒ liftings exist). Uses a diagonal 2-approximation `Δ_P: P→P⊗_A P`. | the equation + existence proof; the SECOND, independent sign convention (arbitration input) |
| **Oke, arXiv:2103.12331** — *Bracket structure on Hochschild cohomology of Koszul quiver algebras using homotopy liftings* | CONFIRMED. Accepted **Comm. Algebra**. | **Definition 3.1**: `d(ψ_η) = (η⊗1_P − 1_P⊗η)Δ_P`, normalization `μ_P ψ_η ~ (-1)^{n-1} η ψ`; **Remark 3.4** (Koszul): the aux `ψ` can be `0`. **Theorem 3.5** (= Volkov): `[η,θ]_Δ = η ψ_θ − (-1)^{(m-1)(n-1)} θ ψ_η`. Worked family in **§7** (§7.1 deg-2 cocycles, §7.2 deg-1 cocycles, §7.3 derivation operators), some members Snashall–Solberg finite-generation counterexamples. | the IMPLEMENTED formulation + the literature-oracle tables (see §7 note below) |
| **Negron–Witherspoon, *An alternate approach to the Lie bracket on Hochschild cohomology*** | CONFIRMED. **Homology, Homotopy Appl. 18(1) (2016) 265–285.** | The homotopy-lifting definition and bracket theorem (the origin of Oke Def 3.1/Thm 3.5); "A∞-coderivations" is the companion (people.tamu.edu/~sjw/pub/infinity_arxiv.pdf). | primary attribution (`bracket_liftings`) |
| **Witherspoon, GSM 204** | CONFIRMED. *Hochschild Cohomology for Algebras*, **AMS Graduate Studies in Mathematics 204 (2019)**. | Textbook chapter deriving the homotopy-lifting bracket; running example HH*(k[x]/x²). | expository anchor; new citation key `witherspoon_gsm204` |

**Sign discrepancy is EXPECTED, not a defect.** Oke/NW give
`[η,θ] = η ψ_θ − (-1)^{(m-1)(n-1)} θ ψ_η` while Volkov gives
`(-1)^m f φ_g + (-1)^{m(n-1)} g φ_f` — the SAME bracket class under a different
homotopy-lifting normalization. This is precisely why the record and Design
Decision 4 make the sign **arbitrated, not assumed** (the Plan-21 cap precedent).

**Oke §7 tables — BLOCKED-until-transcribed.** The ar5iv HTML extraction confirmed
§7 exists and computes explicit homotopy liftings/brackets for degree-1 and degree-2
cocycles of a Koszul quiver-algebra family, but did NOT isolate a clean
quiver+relations+values table this session. Per the metaplan reference rule, **the
implementation MUST open arXiv:2103.12331 §7 (the PDF), transcribe the quiver, the
relations, and every explicit bracket value verbatim with equation/proposition
numbers into the test docstring BEFORE pinning them.** Fabricating pins is a
plan-review defect. Task F ships the k[x]/(x^n) and QuantumCI oracles unconditionally
and carries the Oke §7 pin as an `xfail(reason="Oke §7 not yet transcribed")` fence
that flips to a real assert once transcribed (the Plan-31 auto-flip precedent).

---

## 3. What exists today (audit) and what changes

`grep`-verified against the shipped tree:

| Layer | Today | After P51 |
|---|---|---|
| Class-level bracket | `Comparison.bracket_of_cs_classes(u,v)` — GF(p) only, transport only, `_check_window` raises `NotImplementedError` past the window (`comparison.py:561-577`) | gains `engine="auto"/"native"/"transport"`; native past-window (any Domain), transport in-window byte-unchanged |
| Table-level bracket | `Algebra.gerstenhaber_brackets(top, engine)` → `_product_dispatch("bracket",…)`: `cs`/auto-off-GF(p) → loud `QuiverlabError` (`algebra.py:827-836`); DepthLimitError fallback EXCLUDES bracket (`algebra.py:846`) | `cs` route + auto-off-GF(p) → `cs_bracket_tables`; DepthLimitError auto-fallback enabled for bracket |
| CS table builder | `cs_product_tables` hard-refuses `"bracket"` (`products.py:20-24`) | new `cs_bracket_tables` beside it (cup/cap builder untouched) |
| Native cochain | `cup.py::native_cup`, `cap.py::native_cap` | new `bracket.py::native_bracket` |
| Homotopy lifting | none | new `homotopy_lifting.py` |
| GF(p) in-window bar route | `gfp_product_tables` (`hochschild/products.py:124`), `gerstenhaber_bracket_matrix` (`tt_calculus.py:324`) | **UNCHANGED** — byte-identical (the anchor + the golden) |
| Doc string | `comparison.py:29-30,49-52` "the bracket needs the CS brace/circle machinery to go native (not built in v1)" | rewritten: the homotopy lifting REPLACES the brace; native route now shipped |

---

## 4. Design decisions (with rationale)

### DD1 — Formulation: Negron–Witherspoon / Oke homotopy liftings (not Volkov's formula directly)

The homotopy-lifting construction is the one that maps **one-to-one onto the shipped
Plan-20 machinery**, so we implement it and cite Volkov as the parallel/independent
formulation (and second sign convention).

**The equations to solve.** A CS cochain η representing a class in HH^n is a bimodule
cocycle `η: P_n → A`. Its **homotopy lifting** `ψ_η: P → P` is the degree `-(n-1)`
bimodule map (i.e. `ψ_η: P_k → P_{k-n+1}` on each generator) satisfying, as an
identity of Hom-complex differentials (Oke Def 3.1 / Volkov Def 4):

```
        d_P ∘ ψ_η  −  (−1)^{n-1} ψ_η ∘ d_P   =   (η ⊗ 1_P − 1_P ⊗ η) ∘ Δ_P
```

which, restricted to a free generator σ ∈ S_k (so `ψ_η(σ)` is the unknown PELT element
of the `(o(σ), t(σ))`-corner of `P_{k-n+1}`), is the **full corner linear system**

```
   D_corner(k-n+1, o, t) · ψ_η(σ)  =  [(η⊗1 − 1⊗η)Δ_k(σ)]  +  (−1)^{n-1} · ψ_η(d_k σ)   (★)
```

* **coefficient matrix `D_corner(k-n+1, o, t)`** = a **NEW assembled per-corner
  RESOLUTION-differential matrix** for `d_P: P_{k-n+1} → P_{k-n}` as a bimodule map,
  built by applying `res.d_terms` (glued via `pelt.apply_lower`, the `_d_general`
  idiom) to each free A^e-generator of the `(o,t)` corner of `P_{k-n+1}` and reading
  coordinates against the `(o,t)`-corner PELT basis of `P_{k-n}`. It is the
  **single-complex sibling of `diagonal.py::TensorComplex.tensor_matrix`** (which is
  the tensor-complex version, likewise built off `res.d_terms`). It is **NOT**
  `res.matrix(·, "coh")` — that is the Hom-dual coboundary δⁿ: Cⁿ→Cⁿ⁺¹ (global,
  `a·w·b` collapse, `resolution.py:194-205`), the wrong object and the wrong direction.
  No shipped primitive does this; Task A.0 builds it. (`_d_general` uses `apply_lower`
  on PELTs but exposes no reusable matrix; the diagonal's `tensor_matrix` is the
  closest template.)
* **RHS** = the diagonal half-collapse `(η⊗1 − 1⊗η)Δ_k(σ)` (computed from the shipped
  `diagonal(res, k)` double-PELT by evaluating η on the degree-n tensor factor and
  collapsing to A — the exact `cup._cochain_evaluator` + `ar.mul` arithmetic) PLUS the
  already-solved lower term `ψ_η(d_k σ)`;
* solved by `fields.linalg.solve` (the FULL corner system — like `diagonal.diagonal`,
  NOT a restricted lower-generator solve), then `reduce_mod_nullspace(x, M, dom)` pins
  the **free-variables-zero canonical representative** → ψ is byte-reproducible.

The RHS is a `d_P`-boundary (in `im D_corner`) exactly because η is a cocycle: Volkov
**Lemma 2** proves `(η⊗1 − 1⊗η)` is null-homotopic on `P⊗_A P`, so the lift exists at
every degree where Δ is built. The argument needs ONLY that Δ is a chain map lifting
the identity `A → A⊗_A A`, which quiverlab's degreewise-solved 2-approximation Δ
certifies exactly (`diagonal.py` raises `NotImplementedError` rather than return a
non-chain-map Δ; where Δ exists it IS a chain map). Hence the ψ-solve is **always
consistent at every built degree** — see DD1a (the consistency contract) and the
A.1/A.1a self-certs. This is structurally IDENTICAL to
`diagonal.py::TensorComplex._zeta` + `diagonal` — the only differences are that the
ambient is the single complex P (not P⊗P) and the RHS carries the diagonal
half-collapse. **No contracting homotopy, no brace machinery, no bar object.**

**The bracket cochain** (Oke Thm 3.5), for η ∈ C^p, θ ∈ C^q, is read off on
generators σ ∈ S_{p+q-1} with NO further solve:

```
   [η, θ](σ)  =  η( ψ_θ(σ) )  −  (−1)^{(p-1)(q-1)} · θ( ψ_η(σ) )
```

Degrees check: `ψ_θ: P_{p+q-1} → P_p`, then `η: P_p → A`; `ψ_η: P_{p+q-1} → P_q`, then
`θ: P_q → A`; both land a `(p+q-1)`-cochain. The result is a cocycle whose HH class
is choice-independent (the NW/Volkov theorem); the canonical ψ makes the returned
cochain representative byte-reproducible (the `native_cup` precedent — only the class
is invariant, the representative is pinned).

**Linear-system shapes / cost.** `ψ_η` is a per-generator tower: for each σ ∈ S_k a
solve whose matrix is `D_corner(k-n+1, o, t)` (rows = the `(o,t)`-corner PELT basis of
`P_{k-n}`, cols = the `(o,t)`-corner PELT basis of `P_{k-n+1}`). The dominant cost is Δ
itself (cached on `res` by
`diagonal.py`); a bracket table to `top` needs Δ up to degree `p+q-1 ≤ top` and two ψ
towers per class-pair (cached per class). This is comparable to `native_cup`'s Δ_{p+q}
cost and stays in the deep bucket for the tested fixtures (k[x]/xⁿ, straddle, QCI).

### DD1a — ψ-solve consistency contract: FULL corner system, always consistent by Lemma 2

**The ψ-solve is the FULL corner system** (like `diagonal.diagonal`), NOT a restricted
lower-generator solve (like `_d_general`'s `apply_lower` on a chosen generator subset):
at each degree we solve the whole `D_corner(k-n+1, o, t) · x = RHS` and canonicalize
with `reduce_mod_nullspace`. Rationale: uniformity with the shipped diagonal (same
`tensor_matrix`/`solve`/`reduce_mod_nullspace` shape), and the full system makes the
consistency argument crisp.

**The solve is ALWAYS consistent at every degree the diagonal was built.** Volkov
Lemma 2: for a cocycle η, `(η⊗1 − 1⊗η)` is null-homotopic on `P⊗_A P`, so
`(η⊗1 − 1⊗η)Δ` is a `d_P`-boundary — the RHS of (★) lies in `im D_corner` and `solve`
returns a solution. The argument needs only that Δ is a chain map lifting the identity
`A → A⊗_A A`; quiverlab's degreewise-solved 2-approximation Δ (`diagonal.py`) IS such a
chain map wherever it is built — it raises `NotImplementedError` rather than return a
non-chain-map Δ (the `diagonal.py:302-306` scope edge). So a bracket at degrees within
a built Δ inherits Δ's existence guarantee and NEVER hits an independent
inconsistent-ψ edge.

**Consequence for refusals (retract the phantom edge).** There is therefore **no
independent "inconsistent ψ-lift" scope edge** — the only refusal on the CS path is
Δ's own `NotImplementedError` (raised while BUILDING Δ, before any ψ-solve runs) plus
the cocycle guard (DD3/§6 Task A.4). The `solve is None` branch inside `homotopy_lift`
is a **defensive assertion** ("cannot occur once Δ is built — Lemma 2; if it fires it
is a bug, never an approximation"), NOT a user-facing scope boundary. §8 is corrected
accordingly. **A.1a** self-certifies ψ-consistency at every built degree on the
straddle + QCI fixtures.

### DD2 — Scope: CS in v1; minimal/Bardzell recorded as follow-up

The record names "CS/minimal/Bardzell". The homotopy-lifting method needs a
**diagonal Δ: P → P⊗_A P** on the resolution, and quiverlab ships Δ **only** on the
CS resolution (`resolutions_cs/diagonal.py`). `engine/resolutions_minimal.py` and
`engine/resolutions_bardzell.py` have NO diagonal — building one is Plan-20-scale
work per resolution.

**Decision: v1 delivers the CS-native bracket only.** Minimal/Bardzell are a recorded
follow-up, and the dependency is natural: Oke's own method is stated for **Koszul**
quiver algebras via the **GHMS comultiplicative minimal A^e resolution** — which is
exactly the minimal-resolution diagonal, and is exactly **Plan 75 (R10, GHMS
comultiplicative minimal resolution)**. So the minimal/Koszul native bracket rides on
Plan 75; the Bardzell (monomial) diagonal is its own smaller follow-up. Recorded in
§9 (honest scope) and appended to `docs/plans/DEEPER-ENGINES-BACKLOG.md`.

**Honest cost of the CS-only carrier.** The CS route needs a QUIVER PRESENTATION;
a presentation-less **structure-constants** algebra therefore gets NO native bracket
and the CS route refuses it loudly (the same presentation requirement `cup`/`cap`
already impose off GF(p), Plan 19/35). Such an algebra keeps only its EXISTING
in-window GF(p) transported bracket (the bar route `gfp_product_tables`, which needs
no presentation) — unchanged; past-window/off-GF(p) it is refused, not silently
served. This IS a genuine gap for presentation-less input: the minimal/Bardzell
engines accept structure-constants algebras, so the future minimal-diagonal bracket
(above) is what would close it. For every quiver-presented algebra — the vast
majority of user input, and everything the GUI can express — the native bracket is
fully reachable over any exact Domain via CS.

### DD3 — Circle/brace route NOT needed; the homotopy lifting IS the brace-free route

The bar bracket needs the circle (brace) product `f∘g` (insertion), which has no CS
analogue — that is why `comparison.py` said the CS bracket "would need the CS
brace/circle machinery". The homotopy lifting `ψ` **replaces** the insertion: the
bracket is assembled DIRECTLY from `η∘ψ_θ − sign·θ∘ψ_η` (DD1), no pre-antisymmetrized
circle product, no brace object. Task C rewrites the stale `comparison.py` doc/message
accordingly. (Degree-0 insertion action stays out of scope, same as Plan 35.)

### DD4 — Signs: pinned to Oke/NW, ARBITRATED by the in-window transported bracket

Two literature conventions exist (DD1 + §2). We pin the **Oke/NW** convention
(`[η,θ] = η ψ_θ − (−1)^{(p-1)(q-1)} θ ψ_η`, homotopy lifting
`d_P ψ − (−1)^{n-1} ψ d_P = (η⊗1−1⊗η)Δ`) as the implemented default, but the sign is
**decided by the anchor, not asserted a priori** — exactly the Plan-21 cap pattern
("the sign convention is ARBITRATED, not assumed"). The arbiter chain, all standing
oracles:

1. **THE ANCHOR (the ONLY sign-fixing oracle):** in-window, over GF(p), the native
   `native_bracket(res, …)` ≡ the transported `Comparison.bracket_of_cs_classes(u, v)`
   (the EXISTING transport-only signature, valid in-window — B.2 fixes the sign here,
   BEFORE the `engine=` kwarg exists; Task C re-verifies it via
   `engine="native"`/`"transport"`), **mod coboundary** on every nonzero HH^p × HH^q
   representative pair. The transported bracket is the classical Gerstenhaber sign (via
   `tt_calculus.gerstenhaber_bracket_cochain`), so this test FIXES the implemented ±.
   (`test_native_cup.py::_anchor` pattern.) **The anchor MUST exercise an ODD sign
   exponent** `(p-1)(q-1)` to actually constrain the `(−1)^{(p-1)(q-1)}` factor:
   k[x]/x² is sign-BLIND here (its only nonzero HH pairs have `p` or `q` even, so
   `(p-1)(q-1)` is even and both signs agree). **QuantumCI `q=1` over GF(5) at
   `(p,q)=(2,2)`** is the load-bearing discriminator — `(p-1)(q-1)=1` is odd, and
   `HH²(QCI/GF5)` is nonzero (dims `(2,2,1,0,2,…)`), so the `(2,2)` bracket is
   non-canceling and its sign is genuinely tested. B.2 asserts `(2,2)` is exercised
   and non-vanishing (never a vacuous 0≡0).
2. **graded Jacobi** (new self-cert) — a genuine sign constraint;
3. **Poisson/cup-Leibniz** `[η, θ∪φ] = [η,θ]∪φ + (−1)^{(p-1)q} θ∪[η,φ]` using the
   **native** cup (the GF(p) version is `test_bracket_cup_leibniz`; the native version
   is new and off-GF(p)-capable) — a genuine sign constraint.

**NOT a sign oracle:** graded antisymmetry `[η,θ] = −(−1)^{(p-1)(q-1)}[θ,η]` is a
**consistency check only** — it is TAUTOLOGICAL for the Oke/NW bracket formula
(`[η,θ] = ηψ_θ − sign·θψ_η` is manifestly antisymmetric under swap for ANY choice of
ψ and ANY ± convention, so it cannot distinguish the sign). It is kept (B.1/E.3,
extended to the native route) as a sanity/consistency check, not as a member of the
sign-fixing arbiter chain.

If Oke/NW literal signs disagree with the tt-facade convention (they may, by the
Volkov `(−1)^m` normalization), the anchor wins and the implemented sign is whatever
makes native ≡ transported hold; the discrepancy is folded back into the research doc.

### DD5 — Dispatch: table API stays `auto/bar/cs`; the native/transport switch lives at the class level

There are two engine vocabularies in the shipped tree, at two levels, and they must
NOT be conflated:

* **table level** (public `Algebra.cup_products/cap_products/gerstenhaber_brackets`):
  `engine ∈ {"auto","bar","cs"}` (`_product_dispatch`, `algebra.py:813`);
* **class level** (`Comparison.cup_of_cs_classes/cap_of_cs_classes`):
  `engine ∈ {"auto","native","transport"}` (the Plan-20/21 precedent).

**Decision:** `gerstenhaber_brackets` keeps the **`auto/bar/cs`** table vocabulary
(uniform with `cup_products`/`cap_products` — a whole table family is not uniformly
"in/out of window"; different `(p,q)` cross the window at different points, so
"native/transport" does not lift cleanly to a table family). The record's
`native/transport` precedent is honored where it belongs — the **class level** — by
the new `Comparison.bracket_of_cs_classes(u, v, engine="auto"/"native"/"transport")`,
byte-for-byte mirroring `cup_of_cs_classes`. Routing:

| `gerstenhaber_brackets(engine=)` | GF(p), in bar window | GF(p), past window / DepthLimit | any other exact Domain, presented | structure-constants only |
|---|---|---|---|---|
| `"auto"` | bar/tt (unchanged, `window=top`) | **CS native** (auto-fallback) | **CS native** | loud `QuiverlabError` (path basis needed) |
| `"bar"` | bar/tt (unchanged) | loud (honest wall) | loud (tt facade is GF(p)) | loud |
| `"cs"` | **CS native** | **CS native** | **CS native** | loud (CS needs a presentation) |

`"auto"` over GF(p) in-window is **byte-identical to today** (the golden + the anchor
guarantee it); `"auto"` only routes CS where the bar route used to raise
`DepthLimitError` (the Plan-34 dispatch-amendment precedent) or off GF(p) (where it
used to refuse). `engine="bar"` keeps its honest wall. This lifts the two bracket
refusals at `algebra.py:827-831` and `:846`, and adds `cs_bracket_tables` to the
`engine=="cs"`/auto-off-GF(p) branch and the fallback branch.

### DD6 — Result container + provenance (reuse Plan 35, no new types)

The CS route returns the existing frozen `HHProducts(kind="bracket")`
(`hochschild/products.py:36`) with:

* `engine = "Chouhy-Solotar native diagonal (homotopy lifting)"`,
  `basis = f"cs/{A.domain.name}"`, **`window = None`** (native — no bar bound; the
  block simply omits the window field, so the GUI's "served to degree window {n}
  (bar-transport bound)" line does NOT appear on native tables — a new provenance line
  does, DD7);
* `constants` exact strings (ints mod p on GF(p)-CS, Domain reprs otherwise — the AST
  float gate scans it);
* `references = ["bracket", "gerstenhaber", "bracket_liftings",
  "bracket_liftings_volkov", "chouhy_solotar", "oke_koszul"]` (the CS route appends the
  lifting refs, mirroring `cs_product_tables` appending `chouhy_solotar`);
* Plan-35 explicit representatives (`basis_classes`/`chain_basis`/`differentials`) via
  the existing `_capture_cs` (cohomology side only — the bracket lands in HH).

The GF(p) in-window bar route keeps `window=top`, `engine="hanlab engine (F_p fast
rank)"`, `basis="bar/GF(p)"` — byte-identical.

### DD7 — Canonical keys unchanged; one new golden

The webapp request encodes the kind string `"bracket:0..n"` and never the `engine`
(the runners always call `A.gerstenhaber_brackets(top)` with the default `auto` —
`spec.py:1548`, `runner.py:949`). So **at a fixed library `__version__`, every
existing bracket request keeps its byte-identical Plan-25 canonical key** — no
re-keying from this plan's changes, gated by `test_runner_delegation.py` (existing
goldens byte-identical). (`cache.py::canonical_key` hashes the library version into
the key BY DESIGN, so the normal v0.3.0→v1.0.0 bump rotates ALL keys anyway — this
plan simply adds no per-request re-keying of its own; the version bump is the P80
gate's concern, and the goldens are re-frozen there under the standing version-drift
gate.) ONE new golden is added: a native off-GF(p) or past-window bracket block (a
case that previously errored — so no golden moves, only one is born), documented in
the delegation test's change-log.

---

## 5. Global constraints

- Python is always `.venv/bin/python`; tests
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q …`.
- **No floats in `src/`** — the AST gate (`tests/test_no_floats.py`) scans the two new
  CS modules; all arithmetic through `res.dom`/`AArith`/`fields.linalg`.
- **Every scope boundary is a loud typed refusal** (metaplan §1.3): the ONE genuine CS
  scope edge is Δ's `NotImplementedError` (`diagonal.py:302-306`, the
  higher-CS-homotopy-correction edge), raised while BUILDING Δ before any ψ-solve runs;
  a non-cocycle input to `homotopy_lifting` → a DISTINCT typed error (a cocycle guard,
  never the diagonal's message, MINOR c / A.4); a presentation-less algebra →
  `QuiverlabError` (path basis needed); `engine="bar"` off GF(p) or past window → the
  honest wall. The ψ-solve's own `solve is None` branch is a **defensive assertion**
  ("cannot occur once Δ is built — Volkov Lemma 2; if it fires it is a bug"), NOT a
  user-facing scope boundary (DD1a). Never a silent fallback, never a fabricated value.
- **Buckets auto-assigned by directory** (`tests/conftest.py`): new tests in
  `tests/resolutions_cs/` and `tests/hochschild/` → **deep**; `tests/qpa/` → **qpa**.
  Run touched files by path in development; finish with a `-m deep` spot-run of the
  touched files + `-m fast` + `-m qpa` + `tests/release`.
- **Plan-32 oracle-class markers on every test** (§8): `oracle_selfcert` /
  `oracle_crossengine` / `oracle_literature` / `qpa`; overlap allowed (e.g. the anchor
  is both selfcert-adjacent and crossengine); the `qpa` bucket IS its class.
- **Every new oracle lands on `docs/verification.md`** with a recounted class table
  (`tests/release/test_oracle_classes.py` green) — Task H.
- **Byte-stability gates:** the GF(p) in-window bracket route is byte-unchanged
  (golden + anchor); no existing cache key, golden, or frozen pin moves except the ONE
  documented new golden.
- Conventional commits; green tests at every commit; branch off `dev`, merge back to
  `dev` (dev→main is the P80 gate). Commit/push only when asked.

---

## 6. TDD task list

Each task: write the failing test → run it (see it fail) → implement → run the task's
gates green → commit. Task granularity mirrors Plan 40.

### Task A — the homotopy-lifting tower `ψ_η`

**Files:** create `src/quiverlab/resolutions_cs/homotopy_lifting.py`;
test `tests/resolutions_cs/test_homotopy_lifting.py` (deep).

**Interface:**
```python
def homotopy_lifting(res, eta_vec, n):
    """The canonical homotopy lifting ψ_η of the CS cocycle η ∈ C^n (coordinates over
    res._basis(n,"coh")) as a HomotopyLifting bound to `res`: ψ_η(σ) for σ ∈ S_k is a
    PELT dict over P_{k-n+1}, built degreewise by the FULL corner lift-solve (★) —
    coefficient matrix D_corner(k-n+1, o, t) (the NEW per-corner RESOLUTION-differential
    matrix d_P built off res.d_terms; NOT res.matrix, which is the Hom-dual coboundary),
    RHS the diagonal half-collapse plus the lower ψ term, solve + reduce_mod_nullspace
    (canonical, byte-reproducible). Reuses diagonal(res,·) and cup._cochain_evaluator.
    Raises CocycleError if δη != 0 (a distinct typed guard); the solve is consistent at
    every degree Δ was built (Lemma 2) — a `solve is None` there is a defensive
    AssertionError ('bug, never an approximation'), not a scope refusal (DD1a)."""
```
A `HomotopyLifting` object caches its per-degree towers on `res` (like `TensorComplex`),
exposes `apply(k, generator) -> PELT` and `image_cochain(cocycle_deg, target_deg)`
helpers, and reuses `diagonal.TensorComplex` for Δ and the corner bookkeeping.

- [ ] **A.0 Build the per-corner resolution-differential matrix `D_corner` (the new
  primitive).** Failing test `test_d_corner.py` first: `D_corner(res, m, o, t)` is the
  matrix of `d_P: P_m → P_{m-1}` restricted to the `(o,t)` corner, columns = free
  A^e-generators of the corner of `P_m`, rows = the `(o,t)`-corner PELT basis of
  `P_{m-1}`, entries assembled by applying `res.d_terms(m, ·)` glued through
  `pelt.apply_lower` (the `_d_general` idiom, `resolution.py:114-126`). SELF-CERT:
  `D_corner(m-1) · D_corner(m) == 0` (d²=0) on every corner of the tested fixtures —
  this IS the correctness gate for the new primitive. Implement as the single-complex
  sibling of `diagonal.py::TensorComplex.tensor_matrix` (which builds off the same
  `res.d_terms`). `oracle_selfcert`. **Adjust to reality:** confirm the corner PELT
  basis ordering matches `pelt.py`/`tensor_matrix`'s deterministic order.
- [ ] **A.1 Failing test — the DEFINING equation is satisfied (self-cert).** For
  k[x]/x² over GF(5) and the straddle/QCI fixtures (the `test_native_cup.py` fixtures),
  build `ψ_η` for each basis cocycle η ∈ HH^n (n=1,2) and assert (★) EXACTLY on every
  generator σ ∈ S_k in the built range:
  `d_P(ψ_η(σ)) − (−1)^{n-1} ψ_η(d σ) == (η⊗1 − 1⊗η)Δ_k(σ)` as PELT dicts (with `d_P`
  applied via `res.d_terms`/`apply_lower`, the same primitive `D_corner` wraps). This is
  the homotopy-lifting analogue of `diagonal.py`'s chain-map identity
  (`test_diagonal.py`/`test_native_cup.py::test_deep_qq_diagonal_smoke`). `oracle_selfcert`.
- [ ] **A.1a Failing test — ψ CONSISTENCY at every built degree (DD1a).** On the
  straddle + QCI fixtures over GF(5), assert that for every basis cocycle η the ψ-solve
  is CONSISTENT (returns a solution, never hits the defensive `solve is None`) at every
  degree k where Δ_k is built — the concrete witness of the Lemma-2 guarantee. Also
  assert a NON-cocycle input raises the distinct `CocycleError`, never Δ's
  `NotImplementedError`. `oracle_selfcert`.
- [ ] **A.2 Failing test — Domain-generic smoke over QQ.** Build `ψ` for k[x]/x² over
  **QQ** (Comparison is GF(p)-gated and NOT used; drive the resolution directly like
  `test_deep_qq_diagonal_smoke`) and assert (★) exactly over QQ. `oracle_selfcert`.
- [ ] **A.3 Run → fail** (`ModuleNotFoundError: …homotopy_lifting`; then `D_corner`
  missing).
- [ ] **A.4 Implement** `D_corner` then `homotopy_lifting.py`. **Cocycle guard first:**
  `homotopy_lifting` checks `δη == 0` (apply `res.matrix(n,"coh")` to `eta_vec`) and
  raises a DISTINCT typed `CocycleError` ("homotopy_lifting needs a cocycle; δη ≠ 0")
  BEFORE any solve — never the diagonal's scope message (MINOR c). Degreewise recursion
  from the base degree (target `P_0`: the "no equations" branch of `diagonal.py:297-298`
  — choose the free-variables-zero lift, RHS must vanish; the normalization Oke Def 3.1
  clause 2 / Remark 3.4 folds into this canonical choice). The RHS half-collapse
  `(η⊗1 − 1⊗η)Δ`: from `diagonal(res,k)`, keep double-PELT terms whose η-side factor has
  degree n, evaluate η with `_cochain_evaluator`, `ar.mul`-collapse to A, land the
  surviving factor's PELT in `P_{k-n}`. The `solve is None` branch is a defensive
  `AssertionError` (DD1a) — NOT the diagonal's `NotImplementedError`. **Adjust to
  reality:** read `diagonal.py:225-314` + `resolution.py:114-126` and mirror the
  solve/`reduce_mod_nullspace` idiom.
- [ ] **A.5 Run A.0/A.1/A.1a/A.2 green; run `tests/resolutions_cs/test_diagonal.py`
  (unchanged).**
- [ ] **A.6 Commit** `feat(resolutions_cs): per-corner d_P matrix + homotopy-lifting tower ψ_η on the CS diagonal (Negron-Witherspoon/Volkov, canonical)`.

### Task B — `native_bracket` cochain

**Files:** create `src/quiverlab/resolutions_cs/bracket.py`;
test `tests/resolutions_cs/test_native_bracket.py` (deep). Model on `test_native_cup.py`.

**Interface:**
```python
def native_bracket(res, f_vec, p, g_vec, q):
    """The native CS Gerstenhaber bracket [f,g] of cochains f ∈ C^p, g ∈ C^q, as a
    coordinate vector over res._basis(p+q-1,"coh"): [f,g](σ) = f(ψ_g(σ)) − sign·g(ψ_f(σ))
    for σ ∈ S_{p+q-1}, sign = (−1)^{(p-1)(q-1)} (DD4, arbitrated). Builds two homotopy
    liftings (cached on res). No bar object — any degree, any Domain. Requires p,q ≥ 1."""
```

- [ ] **B.1 Failing tests — the arbiter chain at cochain level:**
  * **Descent:** `[f,g]` is a cocycle when f,g are cocycles (`δ[f,g] = 0` exactly) —
    over GF(5) on k[x]/x² and QCI, the full basis-cocycle grid for `p+q-1 ≤ 2`
    (Δ-degree budget of the deep bucket, the `test_native_cup.py` `_LIGHT` convention).
  * **Antisymmetry (cochain, mod coboundary):** `[f,g] ~ −(−1)^{(p-1)(q-1)}[g,f]`.
  * **`p<1` / `q<1` → `ValueError`** (degree-0 insertion out of scope, DD3).
  `oracle_selfcert`.
- [ ] **B.2 Failing test — IN-WINDOW ANCHOR (the sign-fixing crossengine oracle).**
  The `test_native_cup.py::_anchor` pattern for the bracket: on k[x]/x² over GF(32003)
  and QCI `q=1` over GF(5), for every nonzero HH^p × HH^q pair with `p+q-1 ≤ window`,
  assert `native_bracket(res, u.vec, p, v.vec, q) ≡ comp.bracket_of_cs_classes(u, v)`
  **mod coboundary** (`same_cohomology_class`). B.2 calls the **EXISTING**
  `bracket_of_cs_classes(u, v)` (the transport-only signature, valid in-window,
  `comparison.py:561`) — the `engine=`-routed native/transport comparison is a Task-C
  gate (C.1/C.2), AFTER C.4 adds the kwarg. This DECIDES the DD4 sign. **MINOR a
  (mandatory non-vacuity):** the QCI/GF(5) fixture MUST reach `(p,q)=(2,2)` — the
  odd-exponent `(p-1)(q-1)=1` discriminator — with a NON-canceling result (k[x]/x² is
  sign-blind); assert `(2,2)` is exercised and the bracket there is not identically the
  zero class (raise the window / pick reps so it lands). `oracle_crossengine` +
  `oracle_selfcert`.
- [ ] **B.3 Run → fail.**
- [ ] **B.4 Implement** `native_bracket` (sibling of `native_cup`/`native_cap`;
  `homotopy_lifting` + `_cochain_evaluator` + `ar.mul`). Pin the sign after B.2 decides
  (the QCI(2,2) odd-exponent case is what actually fixes it); if it disagrees with
  Oke/NW literal, record the discrepancy in the module docstring AND fold back into the
  research doc (dated note).
- [ ] **B.5 Run B.1/B.2 green.**
- [ ] **B.6 Commit** `feat(resolutions_cs): native_bracket — CS Gerstenhaber bracket past the bar window, any Domain`.

### Task C — `Comparison.bracket_of_cs_classes` native/transport dispatch

**Files:** modify `src/quiverlab/resolutions_cs/comparison.py`;
test `tests/resolutions_cs/test_comparison.py` (extend) / `test_native_bracket.py`.

- [ ] **C.1 Failing test — engine selector** (mirror
  `test_native_cup.py::test_cup_engine_selector_kx2`): on a TINY-window
  `Comparison(_kx2_gf5(), max_cells=8)` (window 0), `engine="native"` computes a
  past-window bracket, `engine="transport"` still raises `NotImplementedError`,
  `engine="bogus"` raises `ValueError` naming the three options; `engine="auto"` routes
  native past window and transport in-window (byte-identical to the old default there).
  `oracle_selfcert`.
- [ ] **C.2 Failing test — past-window delivery + bridge** (mirror
  `test_native_cup_past_window_kx2` / `test_native_cup_bridge_to_longer_transport_kx2`):
  a bracket past a tiny window computed natively equals — mod coboundary — the same
  bracket by transport on a wider-window instance (Plan-17 canonical CS bases make the
  coordinate vectors directly comparable; assert equal basis lengths first).
  `oracle_crossengine`.
- [ ] **C.3 Run → fail.**
- [ ] **C.4 Implement.** Add `bracket_of_cs_classes(self, u, v, engine="auto")`
  byte-for-byte mirroring `cup_of_cs_classes` (`comparison.py:522-559`): validate
  `engine ∈ {"auto","native","transport"}`; `native` or `auto`-past-window →
  `native_bracket(self._res, u.vec, p, v.vec, q)` after `self._ensure(p+q+1)` (the
  conservative margin matching the cup route's `_ensure(p+q+1)`, `comparison.py:551`;
  the bracket output degree is `p+q-1`, so this is two past the output — safe, and Δ is
  still built lazily only to `p+q-1`); else the
  existing transport body (keep `_check_window`, byte-unchanged). **Rewrite** the stale
  doc/message: the old `bracket_of_cs_classes` docstring (`:566-570`), the `_WINDOW_MSG`
  bracket clause (`:49-52`), and the module-header note (`:29-30`) now say the homotopy
  lifting delivers the bracket natively past the window (DD3) — the brace machinery is
  not needed.
- [ ] **C.5 Run C.1/C.2 + `test_comparison.py` green.**
- [ ] **C.6 Commit** `feat(resolutions_cs): Comparison.bracket_of_cs_classes engine=native/transport (native past window); retire the stale brace-machinery note`.

### Task D — table route: `cs_bracket_tables` + `_product_dispatch` lift

**Files:** modify `src/quiverlab/resolutions_cs/products.py`,
`src/quiverlab/core/algebra.py`;
test `tests/resolutions_cs/test_products_cs.py` (extend),
`tests/hochschild/test_products_api.py` (extend).

- [ ] **D.1 Failing tests:**
  * `Algebra.gerstenhaber_brackets(top, engine="cs")` on k[x]/x² over **QQ** returns an
    `HHProducts(kind="bracket")` with `engine`/`basis` naming the CS route, `window is
    None`, exact-string constants (previously a `QuiverlabError`). `oracle_selfcert`.
  * **Cross-engine table gate** (extend `test_products_identities.py`'s cup/cap gate to
    the bracket): `engine="bar"` vs `engine="cs"` in-window over GF(7)/GF(3) → same dims
    and RANK-equivalent flattened constants per bidegree (bases differ; compare only
    basis-independent data). `oracle_crossengine`.
  * `engine="bar"` off GF(p) → loud; presentation-less off GF(p) → loud.
    `oracle_selfcert`.
- [ ] **D.2 Run → fail** (`QuiverlabError: the Gerstenhaber bracket is served over GF(p)
  only`).
- [ ] **D.3 Implement.** In `products.py` add `cs_bracket_tables(A, top, max_cells)`
  (sibling of `cs_product_tables`: build the CS resolution to **`top+2`** — matching the
  shipped `cs_product_tables` margin conservatively, `products.py:39`; the resolution
  build is cheap and `native_bracket` still lazily builds Δ only to `p+q-1 ≤ top`, so no
  extra Δ degree is forced — `cs_hh_basis` coh per degree, iterate `_pairs("bracket",
  top)`, `native_bracket` per class-pair,
  `_class_coords` descent against `im δ^{p+q-2}`, `HHProducts(kind="bracket",
  engine="Chouhy-Solotar native diagonal (homotopy lifting)", basis=f"cs/{name}",
  window=None, references=…+lifting refs)`, `_capture_cs` coh side). In `algebra.py
  _product_dispatch`: (a) remove the bracket refusal at `:827-831` — route `cs`/auto-off-GF(p)
  bracket to `cs_bracket_tables`; (b) enable the DepthLimitError auto-fallback for
  bracket (drop `or kind == "bracket"` at `:846`, fall back to `cs_bracket_tables`).
- [ ] **D.4 Run D.1 green; run `test_products_api.py`, `test_products_containers.py`,
  `test_products_gfp.py` — the GF(p) route byte-unchanged.**
- [ ] **D.5 Commit** `feat(core,resolutions_cs): gerstenhaber_brackets CS route (any Domain, past window) + DepthLimit auto-fallback`.

### Task E — identity batteries: Jacobi + native Poisson/Leibniz (self-cert), extended

**Files:** modify `tests/hochschild/test_products_identities.py`;
new `tests/resolutions_cs/test_bracket_identities.py` (deep, native/off-GF(p)).

- [ ] **E.1 Graded Jacobi (new).** On k[x]/x^a over `PRIMES=(32003,2,3,5)` (table
  level, GF(p)) and on k[x]/x² over **QQ** (native, off-GF(p)) — the graded Jacobi
  identity `(-1)^{(p-1)(r-1)}[[f,g],h] + cyclic = 0` over the in-scope degree window,
  a content-bearing triple (choose degrees where all three tables are nonzero; document
  non-vacuity like `test_bracket_cup_leibniz`'s fixture note). `oracle_selfcert`.
- [ ] **E.2 Native Poisson/cup-Leibniz (new, off-GF(p)).** `[f, g∪h] = [f,g]∪h +
  (−1)^{(p-1)q} g∪[f,h]` using the **native** cup and native bracket, over QQ on k[x]/x²
  and over GF(5) on QCI past a tiny window — the class-level analogue of
  `test_bracket_cup_leibniz`. `oracle_selfcert`.
- [ ] **E.3 Antisymmetry extended.** The existing
  `test_bracket_antisymmetry` (table, GF(p)) plus a native/off-GF(p) row. `oracle_selfcert`.
- [ ] **E.4 Run green; run the whole `tests/hochschild/test_products_*` + touched
  `tests/resolutions_cs`.**
- [ ] **E.5 Commit** `test(products): graded Jacobi + native Poisson/Leibniz + off-GF(p) antisymmetry for the bracket`.

### Task F — literature oracles

**Files:** modify `tests/hochschild/test_products_literature.py` (extend);
new pins in `tests/resolutions_cs/test_bracket_literature_p51.py` (deep).

- [ ] **F.1 k[x]/x² bracket, off GF(p) (headline "any exact field").** Over **QQ**
  (char 0, `x` is a unit), the CS-native bracket reproduces the classical HH*(k[x]/x²)
  Gerstenhaber structure (Witherspoon GSM 204 running example): pin the KNOWN
  zero-entries (`k[x]/(x^n)` zero-brackets from the record) and the one nonzero
  bracket, all class-level (`same_cohomology_class`, never bytes). Extend the char-2
  survival to the native route (cross the tiny window). `oracle_literature`.
- [ ] **F.2 QuantumCI ties.** Two DISTINCT parts, honestly labelled:
  * **in-window** — CS-native bracket ≡ Plan-35 GF(p) transported bracket over
    GF(2)/GF(5), mod coboundary. This is `oracle_crossengine` (a real comparator).
  * **past-window** — the CS-native bracket COMPUTES past the window and the tables'
    DIMS line up with the BGMS/Bergh–Erdmann QuantumCI HH dims (via
    `test_qci_dims_line_up_with_bgms`, which pins **HH DIMENSIONS, not bracket
    VALUES**). This is `oracle_selfcert` ONLY: there is no transport comparator past the
    window and no published bracket-value table (the Oke §7 values are F.3, blocked). It
    certifies self-consistency + correct dims, NOT a literature bracket value.
- [ ] **F.3 Oke §7 tables — BLOCKED-until-transcribed (see §2).** Add the test skeleton
  with `@pytest.mark.xfail(reason="Oke arXiv:2103.12331 §7 quiver/relations/bracket
  values not yet transcribed from the PDF — DO NOT fabricate; open the PDF, transcribe
  verbatim with proposition/equation numbers, then flip to a real assert")`. The
  implementation session transcribes §7.1 (deg-2 cocycles) / §7.2 (deg-1 cocycles) and
  flips the fence. `oracle_literature`.
- [ ] **F.4 Run green (F.3 xfails as designed).**
- [ ] **F.5 Commit** `test(products): k[x]/x^n off-GF(p) bracket zeros + QuantumCI ties; Oke §7 fence (blocked-until-transcribed)`.

### Task G — GUI / webapp wiring (both runners) + QPA honest-skip

**Files:** `src/quiverlab/hpc/spec.py`, `docs/gui/runner.py` (twin),
`webapp/static/gui/gui.js` + `docs/gui/gui.js` (renderer), `webapp/static/app.js`,
`webapp/server/i18n/{en,es,fr,zh}.json`, `src/quiverlab/trace/results_html.py` +
`src/quiverlab/trace/products.py` (report gloss), `tests/webapp/_runner_goldens.json`
+ `tests/webapp/test_runner_delegation.py`, `tests/qpa/test_products_qpa.py`
(unchanged, re-affirmed);
tests `tests/webapp/test_bracket_native_p51.py`, `tests/hochschild/test_products_api.py`.

The runners already call `A.gerstenhaber_brackets(top)` (`spec.py:1548`,
`runner.py:949`); once Task D lifts the library refusal, the compute "just works" for
off-GF(p)/past-window — the wiring change is **provenance rendering + i18n + the report
gloss + the honest 4xx on a real refusal**.

- [ ] **G.1 Failing cross-runner test** (unmarked, extras-gated dir — copy the
  `test_module_blocks_m0729.py` runner-pair fixture): a `bracket:0..2` request on a
  presented algebra over **QQ** (a) is served by `hpc.spec` with the CS provenance and
  no `window` field, (b) is byte-identical through the Pyodide twin `runner.py`
  (`json.dumps(sort_keys=True)` equality on the block). A GF(p) in-window `bracket:0..n`
  request stays byte-identical to its existing golden. And a structure-constants-only
  off-GF(p) request surfaces the library refusal as a clean typed 4xx (never a 500).
- [ ] **G.2 Implement provenance rendering.** New i18n key
  `block.bracket.native` = "served natively on the Chouhy–Solotar resolution (homotopy
  liftings; any exact field, past the bar window)" in **en/es/fr/zh** (the four files
  already carry `block.bracket.window`); the block renderer in both `gui.js` copies +
  `app.js` shows `block.bracket.window` when `window` is present (GF(p) in-window) and
  `block.bracket.native` when it is absent (CS route). Report: `trace/products.py` +
  `results_html.py` gloss the two routes (the engine-gloss `_ENGINE_GLOSS` precedent —
  homotopy liftings on the CS diagonal). Both `gui.js` copies stay byte-identical
  (`test_js_parses.py`).
- [ ] **G.3 One new golden.** Add `bracket_native_qq` (the QQ past-window bracket block)
  to `_runner_goldens.json`; document it in `test_runner_delegation.py`'s change-log;
  run the delegation test BEFORE adding to confirm the existing bracket golden stays
  byte-identical (DD7).
- [ ] **G.4 QPA honest-skip re-affirmed.** `tests/qpa/test_products_qpa.py` is
  UNCHANGED — QPA 1.37 still has no Hochschild product/bracket surface (the live
  `NamesGVars()` sweep); the test FAILS if that ever changes. Add a one-line comment
  that P51 extends the bracket to CS-native and does not change the QPA scope. `qpa`.
- [ ] **G.5 Run** `tests/webapp/test_bracket_native_p51.py test_runner_delegation.py
  test_js_parses.py tests/hpc -q`; then `-m qpa`.
- [ ] **G.6 Commit** `feat(gui,webapp,hpc): gerstenhaber_brackets CS provenance gloss + i18n ×4; one native golden`.

### Task H — citations, verification page, README, suite gate

**Files:** `src/quiverlab/citations/registry.py` (+ `references.bib`),
`docs/verification.md`, `README.md`, `docs/plans/DEEPER-ENGINES-BACKLOG.md`.

- [ ] **H.1 Citations.** Add `oke_koszul` (Oke, *Bracket structure on Hochschild
  cohomology of Koszul quiver algebras using homotopy liftings*, arXiv:2103.12331,
  Comm. Algebra — BibTeX-verified) and `witherspoon_gsm204` (Witherspoon, *Hochschild
  Cohomology for Algebras*, AMS GSM 204, 2019), tagged `bracket`, via the `_r(...)`
  precedent (`registry.py:24`). REUSE the shipped `bracket_liftings`
  (NegronWitherspoon2016), `bracket_liftings_volkov` (Volkov2019), `bracket`,
  `gerstenhaber`, `chouhy_solotar`. Run `tests/citations/test_bib_structure.py`.
- [ ] **H.2 Verification page.** Add the P51 subsystem row (bracket | oracles:
  `oracle_crossengine` in-window native ≡ transported (incl. QCI(2,2)) + bar ≡ CS
  tables; `oracle_selfcert` `D_corner` d²=0 + homotopy-lifting equation (★) + ψ
  consistency + descent + graded Jacobi + native Poisson/cup-Leibniz (graded
  antisymmetry listed as a CONSISTENCY check, tautological, not a sign oracle);
  `oracle_literature` k[x]/xⁿ off-GF(p) zeros + Oke §7 [pending transcription] (QuantumCI
  DIMS-match is self-cert, not a literature bracket value); `qpa` honest-scope: QPA has
  no bracket surface). Honest-scope entries: (a) **native only on CS in v1** — minimal/Bardzell
  need a diagonal (follow-up on Plan 75 GHMS / a Bardzell diagonal), DD2; (b) the sign
  is **arbitrated** by the in-window anchor (the QCI(2,2) odd-exponent case), not
  assumed (DD4); (c) degree-0 insertion action out of scope (DD3); (d) the ONLY CS scope
  edge is Δ's `NotImplementedError` — the ψ-solve is always consistent by Lemma 2 (DD1a),
  no separate inconsistent-lift refusal; a non-cocycle input raises a distinct
  `CocycleError`; (e) Oke §7 pin **blocked-until-transcribed**. Recount the
  class table (`tests/release/test_oracle_classes.py` drives the numbers — collect,
  paste, re-run green).
- [ ] **H.3 README + backlog.** One feature line: "the Gerstenhaber bracket goes native
  on the Chouhy–Solotar resolution — past the bar window, over any exact field
  (homotopy liftings)". Append the minimal/Bardzell native-bracket follow-up to
  `DEEPER-ENGINES-BACKLOG.md` (deps: Plan 75 GHMS comultiplication / a Bardzell
  diagonal).
- [ ] **H.4 Full gate:** `-m pytest tests/resolutions_cs tests/hochschild -q` (deep),
  `-m fast`, `-m qpa`, `tests/release -q`, plus `QUIVERLAB_NO_NUMBA=1` parity on the
  touched deep files (CS is pure-Python but pin the exact-agreement contract). All
  green.
- [ ] **H.5 Commit** `docs(verification): P51 bracket-beyond-window oracle rows + honest scope (CS-only v1, arbitrated sign, Oke §7 pending) + recounted classes`.

---

## 7. Oracle table (Plan-32 classes)

| Oracle | Class | Where |
|---|---|---|
| `D_corner` per-corner d_P matrix: d²=0 (new primitive correctness) | `oracle_selfcert` | A.0 |
| Homotopy-lifting defining equation (★) holds exactly (GF(p) + QQ) | `oracle_selfcert` | A.1, A.2 |
| ψ-solve consistency at every built degree + `CocycleError` guard | `oracle_selfcert` | A.1a |
| Bracket cochain is a cocycle (descent) | `oracle_selfcert` | B.1 |
| Graded antisymmetry — CONSISTENCY CHECK only (tautological, NOT a sign oracle, DD4) | `oracle_selfcert` | B.1, E.3, existing `test_bracket_antisymmetry` |
| Graded Jacobi (genuine sign constraint) | `oracle_selfcert` | E.1 |
| Native Poisson/cup-Leibniz `[f,g∪h] = [f,g]∪h + (−1)^{(p-1)q} g∪[f,h]` (genuine sign constraint) | `oracle_selfcert` | E.2 |
| **In-window native ≡ transported bracket, incl. QCI(2,2) odd exponent** (THE SIGN ANCHOR) | `oracle_crossengine` | B.2 (transport sig), C.2 (engine=) |
| bar ≡ CS bracket tables in-window (dims + flattened rank) | `oracle_crossengine` | D.1 |
| k[x]/(x^n) known zero-brackets, off GF(p) | `oracle_literature` | F.1 |
| QuantumCI ties — IN-WINDOW native ≡ transported | `oracle_crossengine` | F.2 (in-window) |
| QuantumCI past-window computes + DIMS match BGMS (not values) | `oracle_selfcert` | F.2 (past-window) |
| Oke arXiv:2103.12331 §7 Koszul-quiver bracket tables | `oracle_literature` (BLOCKED-until-transcribed) | F.3 |
| QPA exposes no bracket surface (honest scope) | `qpa` | G.4 (unchanged) |

---

## 8. Honest-scope statement (contractual, on `docs/verification.md`)

1. **Native carrier = CS only in v1.** The homotopy-lifting bracket rides on the
   Plan-20 diagonal Δ, which lives only on the Chouhy–Solotar resolution. Minimal and
   Bardzell resolutions have no diagonal; their native bracket is a recorded follow-up
   (minimal/Koszul ← Plan 75 GHMS comultiplication; Bardzell ← a monomial diagonal).
   The CS route serves every **quiver-presented** algebra over every exact Domain, so
   the bracket is fully reachable for all presented input (and everything the GUI can
   express). A presentation-less **structure-constants** algebra is NOT served natively:
   the CS route refuses it loudly (same presentation requirement as `cup`/`cap` off
   GF(p)); it keeps only its existing in-window GF(p) transported bracket. Closing that
   gap needs the minimal/Bardzell diagonal (the follow-up above), since those engines
   accept structure-constants input. (DD2)
2. **The sign is arbitrated, not assumed** — fixed by the in-window native ≡ transported
   anchor over GF(p), then re-certified by antisymmetry/Jacobi/Poisson. Any deviation
   from the Oke/NW literal sign is recorded in the module docstring and folded back into
   the research doc. (DD4)
3. **Degree-0 insertion action out of scope** (`p,q ≥ 1`), same as Plan 35. (DD3)
4. **Scope edge = loud refusal (and the ψ-solve is NOT one).** The ONE genuine CS scope
   edge is Δ's `NotImplementedError` (the higher-CS-homotopy-correction edge, raised
   while BUILDING Δ, before any ψ-solve). The ψ-solve itself is ALWAYS consistent once Δ
   is built (Volkov Lemma 2, DD1a) — there is NO independent "inconsistent ψ-lift"
   refusal; its `solve is None` branch is a defensive `AssertionError` ("bug, never an
   approximation"), not a user-facing boundary. A non-cocycle input raises a distinct
   `CocycleError`; presentation-less algebras and `engine="bar"` off GF(p)/past window
   raise loudly. Never a silent fallback. (§5, DD1a)
5. **Oke §7 pin blocked-until-transcribed** — no fabricated literature values; the fence
   flips only after verbatim transcription from the PDF. (§2, F.3)
6. **QPA gap** — QPA 1.37 has no Hochschild product/bracket surface; the covering oracle
   is the identity battery + the literature pins + the in-window anchor. (G.4)

---

## 9. Citation keys

REUSED (shipped in `registry.py`): `bracket_liftings` (NegronWitherspoon2016),
`bracket_liftings_volkov` (Volkov2019), `bracket` (Gerstenhaber1963), `gerstenhaber`
(Gerstenhaber1963), `chouhy_solotar` (ChouhySolotar2015).
NEW: `oke_koszul` (Oke2021, arXiv:2103.12331, Comm. Algebra), `witherspoon_gsm204`
(Witherspoon, AMS GSM 204, 2019). All BibTeX-verifiable.

---

## 10. Risk register

- **The ψ base-degree ("no equations") branch.** At the smallest degree the target is
  `P_0` and `d_0 = 0`, so ψ(σ) is unconstrained by the solve; the RHS must vanish
  (degree reasons + Lemma 2) and the canonical (free-vars-zero) lift is chosen — the Oke
  Def 3.1 normalization folds in here. Mirror `diagonal.py:297-298` EXACTLY; a nonzero
  RHS there is IMPOSSIBLE for a cocycle with Δ built (DD1a), so it is a defensive
  `AssertionError` ("bug"), never a scope edge and never a silent zero. Self-certified by
  A.1 (the defining equation holds at every degree) and A.1a (consistency).
- **Sign arbitration must run BEFORE Task E/F pin any value, and MUST hit an odd
  exponent.** B.2 (the anchor) is the gate that fixes the sign, and it MUST reach
  QCI(2,2) (odd `(p-1)(q-1)`, non-canceling) — k[x]/x² alone is sign-blind (DD4). Do not
  pin Jacobi/Poisson/literature until B green.
- **Δ-degree cost.** The bracket to `top` needs Δ up to `p+q-1`; two ψ towers per
  class-pair. Keep the deep-bucket fixtures small-window (the `test_native_cup.py`
  Task-4/5 pattern: a tiny/zero comparison window makes "past-window" land at LOW
  absolute degree, so only the Δ build costs). Session-scope the heavy Δ fixture (the
  `qci_gf5_diag4` precedent).
- **Byte-stability of the GF(p) route.** The single largest regression risk is nudging
  `gfp_product_tables`/`tt_calculus`. Task D must NOT touch them; the golden + the
  anchor are the tripwires.
- **Oke §7 transcription.** If the PDF is unreadable at implementation time, F.3 stays
  xfail and the honest-scope entry names it — the plan does NOT block on it, and the
  k[x]/xⁿ + QuantumCI + anchor oracles carry the literature/crossengine burden.

---

## 11. Acceptance (definition of done)

1. The new per-corner resolution-differential primitive `D_corner` (off `res.d_terms`,
   d²=0-certified), `homotopy_lifting` (with the `CocycleError` guard), `native_bracket`,
   `cs_bracket_tables`, and `Comparison.bracket_of_cs_classes(engine="auto"/"native"/
   "transport")` public on the CS route; `Algebra.gerstenhaber_brackets(top,
   engine="cs"/"auto")` serves the bracket NATIVELY over any exact Domain and past the
   bar window. The ONLY CS scope edge is Δ's `NotImplementedError`; the ψ-solve is always
   consistent (Lemma 2, DD1a) — no phantom inconsistent-lift refusal.
2. `D_corner` d²=0 and the homotopy-lifting equation (★) are self-certified exactly
   (GF(p) + QQ); ψ consistency is self-certified at every built degree; the bracket
   descends (cocycle); graded Jacobi and native Poisson/cup-Leibniz hold (antisymmetry is
   a consistency check, not a sign oracle); the in-window native ≡ transported anchor —
   **exercising the QCI(2,2) odd exponent** — fixes the sign (crossengine); bar ≡ CS
   tables in-window.
3. Literature: k[x]/(x^n) off-GF(p) zero-brackets pinned; QuantumCI in-window tie pinned
   (crossengine) and past-window dims-match (self-cert, not a bracket value); the Oke
   §7 fence present (xfail until transcribed, never fabricated).
4. GF(p) in-window bracket route BYTE-UNCHANGED (golden + anchor); ONE new native golden
   documented; at a fixed library `__version__` every existing bracket request keeps its
   canonical key (the version bump rotates keys by design — a P80 concern).
5. `gerstenhaber_brackets` clickable end-to-end (GUI canvas → block → report) with the
   CS provenance gloss in EN/ES/FR/ZH, both runners byte-identical; the QPA honest-skip
   re-affirmed.
6. The stale `comparison.py` "needs the CS brace/circle machinery" note retired (DD3).
7. `docs/verification.md` recounted with the P51 rows + the six honest-scope entries;
   README + backlog updated; deep (touched dirs) + fast + qpa + release + no-numba parity
   green.
