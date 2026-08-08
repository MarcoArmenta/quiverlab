# Plan 78: L∞ structure / Maurer–Cartan formal deformations (P78 / R13)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking. Do the reference re-verification BEFORE
> writing any oracle pin — the `# PIN` fields must be resolved against the sources, not this
> plan's prose. **Three route/scope decisions are binding and benchmark-grounded (evidence in
> the Route decision and Live-verified sections):** (1) the computational backbone is
> **Chouhy–Solotar over ℚ** — the shipped native CS differential (ℓ₁) and native CS Gerstenhaber
> bracket (ℓ₂, P51) — **NOT** a char-0 Bardzell HH port; the benchmark shows CS-over-ℚ clears the
> Maurer–Cartan degrees (≤ 3) at oracle scale and the Bardzell speedup only helps at high
> degree / large dim that MC does not use. (2) The formal-deformation *functor* is computed on
> the Hochschild **DGLA** `C(A)` (Gerstenhaber's theorem — a genuine dg-Lie algebra, **not**
> "pretending dg-Lie"), where ℓ₂ + δ suffice and no higher bracket is ever needed; the
> `B(A)[1]` L∞ structure is the small-model *companion* whose only exact, certified instances in
> v1 are the **rad²=0 dg-Lie collapse** and a **feasibility-gated char-0 Bardzell ℓ₃** for the
> monomial classes RRB compute explicitly. (3) The classic `k[x]/(x²) ⤳ k[x]/(x²−t)`
> semisimplification is a **unit-direction** deformation that leaves the admissible `kQ/I` class
> — it is an MC/2-cocycle story on `C(A)`, **not** a presented-`A_α` feedback; the presented
> feedback covers the **radical** (admissible) sub-locus (RRRV, alg-closed).

**Goal.** A char-0 formal-deformation surface for a finite-dimensional `A = kQ/I`, built on the
first Hochschild-cohomology deformation theory and the L∞ machinery of Redondo–Rossi Bertone /
Müller–Redondo–Rossi Bertone–Suarez:

- **Infinitesimal deformations** = `HH²(A)` (the tangent space to the deformation functor), over
  any exact field, from the shipped CS engine.
- **The primary obstruction** `α ↦ [α,α] ∈ HH³(A)` — the ℓ₂ (Gerstenhaber) self-bracket, computed
  natively over ℚ by the P51 CS homotopy-lifting bracket. `α` lifts to second order iff
  `[α,α] = 0` in `HH³`. A deformation is **unobstructed** when the obstruction map vanishes as a
  class map on all of `HH²`; otherwise a **witness** obstructed direction is reported.
- **Maurer–Cartan elements = formal deformations.** In the **nilpotent regime** (MRRS Thm 5.4 —
  gentle under quiver hypotheses) the MC set **equals the 2-cocycles** `Z²` (every infinitesimal
  deformation integrates; unobstructed). Outside the regime, the honest MC picture is the DGLA
  equation on `C(A)` (`δμ + ½[μ,μ] = 0`) solved order-by-order to a **certified truncation
  order**, with the obstruction structure exposed.
- **The presented deformed algebra `A_α`** (RRRV) fed back into the engine: for a **radical**
  cocycle direction and a value `t`, build `A_α = kQ/I_α` via `Quiver.algebra`, **re-certify
  admissibility/flatness** (Gröbner; `dim A_α = dim A`), and recompute its invariants and its
  **Ext-algebra** (RRR 2202.01199, Plan 27) — the `A_α` handoff.
- **The `B(A)[1]` L∞ companion (honest L∞, char 0) — the worst case ships NO `ℓ_{≥3}` bracket.**
  The **guaranteed** deliverable is two things: (i) the Hochschild **DGLA** structure on `C(A)`
  (§1/§4 — the whole MC/obstruction backbone is `ℓ₂ + δ` on a genuine dg-Lie algebra, by
  Gerstenhaber's theorem), and (ii) the **rad²=0 ⇒ dg-Lie collapse certificate** (`ℓ_{≥3} ≡ 0`,
  RRB, decidable from the field-free Bardzell substrate). Beyond those, a **feasibility-gated**
  char-0 Bardzell **ℓ₃** for monomial `A` (validated against RRB's explicit truncated
  computations) ships **only if** the Task-4 spike lands; the route that would produce the higher
  `ℓ_n` is RRB's explicit homotopy-transfer / contracting-homotopy on Bardzell's complex
  (`rrb_linfty_bardzell`, arXiv:2008.08122). **If the spike does not land, no `ℓ_{≥3}` bracket
  ships at all** — a ledgered deferral, with the §4 DGLA-equivalence as the honest justification;
  that is a complete v1, not a gap. The induced-`ℓ₂` ≡ CS-bracket agreement is stated as the §4
  **model-independence theorem** (cited, not vacuously computed) and is promoted to a live
  crossengine oracle **only** when the spike builds the `B(A)` `ℓ₂` adapter (computing `B(A)`'s
  `ℓ₂` needs exactly that adapter — see Task 4 / the M1 note). **Never claim `B(A)` is dg-Lie
  outside rad²=0, and never claim `ℓ₄ = 0`** ("ℓ_n=0 for n≥5" is a *sufficient* collapse
  condition; ℓ₄ can be nonzero — honest L∞).

**GUI.** ONE new **algebra-only** budget-carrying compute kind `deformations`, served by all
three tiers (`hpc/spec.py` + the Pyodide twin `docs/gui/runner.py` byte-identical + the
`webapp/server/schema.py` grammar site), the deformation/recognizer panel row, all four locales
(en/es/fr/zh) with exact key parity, canonical keys stable, report rendering, and citations. v1
is **report-only**: the block reports `HH²`/the obstruction/the nilpotent verdict and **displays**
the presented `A_α` (quiver + deformed relations) for a canonical radical direction; the
GUI-**adopt** flow (load `A_α` back onto the canvas as a fresh input) is an explicit **GUI-deferral
ledger** entry for P80 (schema-heavy; needs Marco's sign-off).

**Architecture.** ONE new source module + one algebra-only kind; the deformation computations are
a thin exact layer over primitives already on `dev` (CS bracket/differential, `HH`, `Quiver.algebra`,
Gröbner, `ext_algebra`) plus — behind a feasibility gate — a char-0 adapter over the **field-free**
Bardzell combinatorics for ℓ₃:

- **`src/quiverlab/hochschild/deformations.py`** (new) — the formal-deformation surface, beside
  `products.py` (Plan 35, **already on `dev`**) and `lie_module.py` (Plan 71, **implemented-first**:
  only P71's *doc* is committed — `lie_module.py` is **not on `dev` yet**, so this "beside" is a
  landing-spot note, not a live dependency; P78's own compute imports neither — see the P70/P71
  dependency note in Global Constraints / Acceptance §7). Public:
  - `infinitesimal_deformations(A, *, engine="auto", max_cells=...) -> dict` — `HH²` tangent data.
  - `obstruction_map(A, *, engine="auto", max_cells=...) -> Obstruction` — the quadratic
    `α ↦ [α,α] ∈ HH³` (the (2,2) CS self-bracket) + unobstructed verdict + witness.
  - `maurer_cartan(A, *, order=..., engine="auto", max_cells=...) -> MCReport` — nilpotent-regime
    MC = `Z²`, else the order-by-order DGLA solution + certified truncation order.
  - `is_l_infinity_nilpotent(A) -> bool | None` — the MRRS quiver-hypothesis gate (char-0, gentle).
  - `dg_lie_certificate(A) -> bool | None` — rad²=0 ⇒ `ℓ_{≥3} ≡ 0` (monomial, RRB).
  - `l3_bracket(A, *, budget=...) -> L3Report` — feasibility-gated char-0 Bardzell ℓ₃ (monomial).
  - `deformed_algebra(A, direction, *, t="1") -> Algebra` — build the presented `A_α`, re-certify.
  - `deformation_structure(A, *, budget=..., engine="auto") -> Deformations` — the full report.
  - `deformations_block(A, *, budget=...) -> dict` — the runner block.
- `Algebra.deformation_structure(...)` / `Algebra.deformed_algebra(...)` /
  `Algebra.obstruction_map(...)` — thin lazy-import delegates (beside `gerstenhaber_brackets`).
- The `deformations` kind is wired into `hpc/spec.py::_dispatch` + `docs/gui/runner.py`
  (byte-identical twins) + `webapp/server/schema.py`, the GUI touchpoints, the i18n chains, and
  `trace/results_html.py`.

**Tech Stack.** Exact `Domain` arithmetic only; **no floats in `src/`** (AST-gated by
`tests/test_no_floats.py`): cocycle coordinates and obstruction constants are exact-`Domain`
strings; dims/orders `int`; verdicts `bool`/`str`. The CS bracket/differential and `HH` are the
substrate (`resolutions_cs`, over any exact Domain incl. ℚ); the ℓ₃ adapter reuses the field-free
`engine/resolutions_bardzell.MonomialPresentation` combinatorics with the differential's integer
±1 coefficients reinterpreted over ℚ (see Route decision). No new numeric kernels.

---

## Record (verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`, R13)

> **R13 — L∞ / Maurer–Cartan deformations on Bardzell's complex.** [A-scout P3;
> keep-with-corrections]
> Object: the ℓₙ brackets on B(A)[1] for monomial/gentle A; MC set = formal
> deformations; obstruction [α,α] ∈ HH³; the PRESENTED deformed algebra A_α + its
> Ext-algebra (feed back into the engine). Honest scope: char 0; "ℓₙ=0 for n≥5" is a
> SUFFICIENT condition and ℓ₄ can be nonzero (honest L∞, not DGLA); MC =
> 2-cocycles holds in the NILPOTENT regime (gentle, no parallel arrows, no
> oriented cycles — Thm 5.4). COST FLAG (critic): quiverlab's Bardzell engine is
> GF(p) int64 — this needs a char-0 Bardzell path (CS runs over any Domain and is
> the fallback). Refs: Redondo–Rossi Bertone arXiv:2008.08122 (JPAA 226(5) 2022);
> Müller–Redondo–Rossi Bertone–Suarez arXiv:2309.02582 (Comm. Alg. 2025);
> Redondo–Román–Rossi Bertone–Verdecchia arXiv:2003.10366 (Morita invariance);
> Redondo–Román–Rossi Bertone arXiv:2202.01199 (Ext-algebra of deformations —
> THREE authors, no Verdecchia); Chouhy arXiv:1708.02933 (degeneration link).
> Size L.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P78 (Wave 4, tier δ; **needs
P70/P71**; the char-0 Bardzell path or the CS fallback decided at spec time with a benchmark).
Size L.

---

## Reference re-verification (done at authoring; findings binding, `# PIN` = resolve at citation-add)

The metaplan standing rule (§1.5) requires the plan writer to re-verify its record's citations.
All five arXiv abstracts were fetched this authoring round; the two flagged items (RRR "three
authors"; the journal refs the card asserts) were resolved by search. Findings, all binding:

1. **Redondo, M. J. & Rossi Bertone, F. — "L∞-structure on Bardzell's complex for monomial
   algebras."** arXiv:2008.08122. **Two authors** (verified). **Abstract verified verbatim:** *"Let
   A be a monomial associative finite dimensional algebra over a field 𝕜 of characteristic zero …
   we describe an explicit L∞-structure on B(A) that induces a weak equivalence of L∞-algebras
   between B(A) and the Hochschild complex C(A) … This allows us to describe the Maurer–Cartan
   equation in terms of elements of degree 2 in B(A). Finally … we prove that Bardzell's complex
   for radical square zero algebras is in fact a dg-Lie algebra."* **Journal ref CONFIRMED (search):
   J. Pure Appl. Algebra 226 (2022), no. 5, Paper No. 106935.** MSC 16S80, 17B55, 18G35. (The arXiv
   title string carries a typo "Barzdell"; the journal title is "Bardzell's".) This grounds: the
   char-0 scope, the L∞-on-B(A), MC-in-degree-2, and the **rad²=0 ⇒ dg-Lie** certificate.
   `# PIN`: article number 106935 / DOI against the ScienceDirect landing page at citation-add.

2. **Müller, M., Redondo, M. J., Rossi Bertone, F. & Suarez, P. — "Maurer–Cartan equation for
   gentle algebras."** arXiv:2309.02582. **Four authors** (verified: Müller, Redondo, Rossi
   Bertone, Suarez — the record's "MRRS"). **Abstract verified verbatim:** *"Let A = 𝕜Q/I be a
   finite-dimensional gentle algebra. In this article, under some hypothesis on the quiver Q, we
   give conditions for nilpotency of the L∞-structure on the shifted Bardzell's complex B(A)[1].
   For nilpotent cases, we describe Maurer–Cartan elements."* **Journal ref CONFIRMED (search):
   Comm. Algebra 53 (2025), no. 9, 3984–4007, DOI 10.1080/00927872.2025.2473029.** Grounds: the
   nilpotent-regime gate and MC = 2-cocycles. `# PIN`: the **exact Thm-5.4 quiver hypotheses**
   (R13 states "gentle, no parallel arrows, no oriented cycles" — confirm the precise statement
   against §5 of the paper before writing the gate; do not guess).

3. **Redondo, M. J., Román, L., Rossi Bertone, F. & Verdecchia, M. — "Morita invariance for
   infinitesimal deformations."** arXiv:2003.10366. **Four authors** (verified: adds Verdecchia —
   the record's "RRRV"). **Abstract verified verbatim:** *"… we explicitly describe the transfer
   map connecting HH²(A) with HH²(B) … transfer Morita equivalence … to that between
   infinitesimal deformations … As an application, when 𝕜 is algebraically closed, we consider the
   quotient path algebra associated to A and describe the presentation by quiver and relations of
   the infinitesimal deformations of A."* Grounds: the **presented `A_α` by quiver-and-relations**
   (alg-closed) — the feedback build. No journal ref on the arXiv page → ship `@misc` with
   `note = {arXiv:2003.10366}` unless a journal ref is confirmed. `# PIN`: journal ref.

4. **Redondo, M. J., Román, L. & Rossi Bertone, F. — "The Ext-algebra for infinitesimal
   deformations."** arXiv:2202.01199. **THREE authors CONFIRMED** (Redondo, Román, Rossi Bertone —
   **no Verdecchia**; the card's "three authors" flag RESOLVED). **Abstract verified verbatim:**
   *"Let f be a Hochschild 2-cocycle and let A_f be an infinitesimal deformation of an associative
   finite dimensional algebra A over an algebraically closed field 𝕜. We investigate the algebra
   structure of the Ext-algebra of A_f and, under some conditions on f, we describe it in terms of
   the Ext-algebra of A … an explicit construction of minimal projective resolutions in mod A_f."*
   Grounds: the **Ext-algebra of `A_α`** handoff. No journal ref on arXiv → `@misc`. `# PIN`:
   journal ref (v2 2024-03-15).

5. **Chouhy, S. — "On geometric degenerations and Gerstenhaber formal deformations."**
   arXiv:1708.02933. **One author** (verified). **Abstract verified verbatim:** *"We study the
   degeneration relations on the varieties of associative and Lie algebra structures … and give a
   description of them in terms of Gerstenhaber formal deformations … For … finite dimensional
   associative algebras, we prove that the N-Koszul property is preserved under the degeneration
   relation for all N ≥ 2."* DOI `10.1112/blms.12277` (**Bull. Lond. Math. Soc.**) shown on the
   page. Grounds: the **degeneration ↔ formal-deformation** link (the geometric reading of `A ⤳ A_α`).
   `# PIN`: volume/year/pages against the Bull. LMS landing page (DOI confirmed).

**Reused (verify live in `registry.py` at implementation):** `chouhy_solotar` (CS), the Plan-51
bracket keys (`bracket_liftings`, `bracket_liftings_volkov`, `oke_koszul`), `gerstenhaber1963`
(added by P70), the Plan-27 `ext_algebra` citations (`gsz_book`/GSZ2001), `assem_book` (ASS2006).
**New keys to add:** `rrb_linfty_bardzell` (2008.08122), `mrrs_mc_gentle` (2309.02582),
`rrrv_morita_deform` (2003.10366), `rrr_ext_deform` (2202.01199), `chouhy_degeneration`
(1708.02933).

**Bonus references spotted (optional, NOT required — do not add unless an implementation task
uses them):** Artenstein–Redondo–… arXiv:2401.06429 "higher structures … toupie algebras" (a
follow-on L∞ note); Saorín–Solotar-school "semisimplicity of outer derivations of monomial
algebras" arXiv:1004.2820 (Strametz-adjacent).

---

## Mathematical foundation (the definitions, exactly — the plan's ground truth)

Throughout `A = kQ/I` is a finite-dimensional associative unital `k`-algebra, `k = A.domain` an
exact `Domain`, composition **left-to-right**. `C(A) = ⊕_n Hom_k(Ā^{⊗n}, A)` is the normalized
Hochschild cochain complex (`Ā = A/k·1`), `δ` its differential, `[·,·]` the Gerstenhaber bracket,
`HH^n(A) = H^n(C(A), δ)`. `B(A)` is Bardzell's (small) complex, quasi-isomorphic to `C(A)` and
computing the same `HH`.

### 1. The Hochschild DGLA and Gerstenhaber formal deformations (any field; the ℓ₂ theory)

`(C(A), δ, [·,·])` is a **differential graded Lie algebra** (Gerstenhaber 1963) after the shift
`C(A)[1]` (so `Hom(Ā^{⊗2}, A)` sits in degree 1). A **formal deformation** of `A` is a
Maurer–Cartan element `μ = Σ_{k≥1} t^k μ_k ∈ t·C²(A)[[t]]`:
```
δμ + ½[μ,μ] = 0     (the associativity of the deformed product a·b + Σ t^k μ_k(a,b))
```
solved order by order:
- **order 1:** `δμ₁ = 0` — the infinitesimal `μ₁ ∈ Z²(A)`; its class lives in `HH²(A)`;
- **order k ≥ 2:** `δμ_k = −½ Σ_{i+j=k, i,j≥1} [μ_i, μ_j]`.

The **primary obstruction** to lifting an infinitesimal `α ∈ HH²` to second order is the class
`[α,α] ∈ HH³(A)`: a lift `μ₂` exists iff `[α,α] = 0` in `HH³`. This is the **card's obstruction
`[α,α] ∈ HH³`**. Because `α` has (unshifted) degree 2, graded antisymmetry gives
`[α,α] = [α,α]` (the sign is `+1`), so the self-bracket is **not forced to vanish** — it is a
genuine quadratic obstruction map `HH² → HH³`.

**This ℓ₂ theory is field-general**; the interpretation as *formal* deformation needs the ½ (char
≠ 2) and, at higher order, the `1/n!` of the L∞ formulation (char 0). Deformation-theoretic
outputs are therefore **char-0 gated** (the papers' scope); the raw `HH²/HH³/[α,α]` compute over
any exact field.

### 2. The L∞ structure on `B(A)[1]` (RRB, monomial, char 0)

`B(A)` is much **smaller** than `C(A)`. To do MC theory on the small model, RRB transfer the DGLA
structure of `C(A)` along the quasi-isomorphism `B(A) ≃ C(A)`; the homotopy-transfer theorem yields
an **L∞ structure** `{ℓ_n}_{n≥1}` on `B(A)`:
- `ℓ₁` = Bardzell's differential;
- `ℓ₂` = the transferred (reduced) Gerstenhaber bracket; it **induces the Gerstenhaber bracket on
  `HH`** (so at the cohomology level `ℓ₂` = the shipped CS bracket — the anchor);
- `ℓ_{≥3}` are **generically nonzero** (honest L∞, **not** DGLA). *"ℓ_n = 0 for n ≥ 5"* is a
  **sufficient collapse condition, not automatic**, and **ℓ₄ can be nonzero** (R13, binding).
- **rad²=0 ⇒ dg-Lie:** for a monomial algebra with `rad² = 0`, RRB prove `B(A)` is a genuine
  **dg-Lie algebra** (`ℓ_{≥3} ≡ 0`). This is a *decidable structural* statement (is `rad² = 0`?)
  and the one exact ℓ_{≥3} claim v1 ships as a certificate.

The MC equation on `B(A)[1]` is `Σ_{n≥1} (1/n!) ℓ_n(α^{⊗n}) = 0` for `α` of degree 2 in `B(A)`
(degree 1 in `B(A)[1]`); RRB describe it **in terms of degree-2 elements**.

### 3. MC = 2-cocycles in the nilpotent regime (MRRS Thm 5.4, gentle, char 0)

MRRS prove: under quiver hypotheses on a gentle `A = kQ/I` (R13: *gentle, no parallel arrows, no
oriented cycles* — `# PIN` the exact Thm-5.4 statement), the L∞ structure on `B(A)[1]` is
**nilpotent**, and the **Maurer–Cartan set equals the 2-cocycles `Z²`**: every infinitesimal
deformation solves MC (the higher-bracket corrections vanish on `Z²`), hence **every first-order
deformation integrates and is unobstructed**. Outside this regime the honest picture is the full
MC equation (§1 on `C(A)`, or §2 with the `ℓ_n` corrections on `B(A)`), with the obstruction
`[α,α]` as the first non-vanishing term.

### 4. Why CS-over-ℚ suffices for the deformation *functor* (the route crux — binding)

`C(A)` is a **DGLA** (ℓ₂ only — Gerstenhaber's theorem, not conjecture); `B(A)` is only **L∞**
(the `ℓ_n`). They are **L∞-quasi-isomorphic**, so the deformation/MC *functors* agree (MC sets and
gauge-equivalence classes correspond under the quasi-iso). Therefore the plan's deformation
deliverables — infinitesimal `HH²`, the obstruction `[α,α]`, the DGLA order-by-order MC,
MC = 2-cocycles in the nilpotent regime, the presented `A_α`, the Ext-algebra — are computable
**on `C(A)` via CS-over-ℚ** using only `δ` (native CS differential) and `[·,·]` (native CS
Gerstenhaber bracket, P51), with **no higher bracket ever**. Working on `C(A)` is **not**
"pretending dg-Lie": `C(A)` *is* a DGLA by a theorem. The `B(A)[1]` L∞ structure (`ℓ_{≥3}`) is a
**small-model presentation** whose exact, certified v1 instances are the rad²=0 dg-Lie collapse
and the feasibility-gated ℓ₃ — it enriches the picture but is **not on the critical path** for any
MC/obstruction/feedback output.

### 5. The presented deformed algebra `A_α` (RRRV) and its Ext-algebra (RRR)

Given an infinitesimal 2-cocycle `f` (a deformation direction) and a value `t ∈ k`, the
first-order deformed product is `a ∗ b = ab + t·f(a,b)`. When `f` is supported on **paths of
length ≥ 2** (a **radical** direction), `A_f` stays in `rad²` and is presented by **the same
quiver `Q` with deformed relations `I_f`** — build it with `Quiver.algebra` and **re-certify**
admissibility + finiteness (Gröbner; the honest *flatness* check is `dim A_f = dim A`). Feed back:
recompute `HH•`, other invariants, and the **Ext-algebra `E(A_f)`** (RRR 2202.01199, via the
Plan-27 `ext_algebra`) — RRR describe `E(A_f)` in terms of `E(A)` under conditions on `f`.

**Honest scope (binding).** A **unit-direction** cocycle (e.g. `f(x,x) = 1` sending
`k[x]/(x²) ⤳ k[x]/(x²−t)`) leaves the admissible `kQ/I` class: the deformed relation acquires a
length-0 term and the radical *shrinks* (the fibre becomes semisimple, the Gabriel quiver
changes). Such deformations are expressed as an **MC/2-cocycle story on `C(A)`** (with the shipped
`HH²`), **not** as a presented-`A_α` feedback. RRRV's quiver-and-relations presentation holds
(over algebraically closed `k`) for the **admissible (radical) sub-locus** — the plan's
presented-feedback scope. Chouhy 1708.02933 supplies the geometric reading: `A ⤳ A_α` is a
degeneration; N-Koszulity is preserved.

---

## Route decision — the MANDATORY spec-time benchmark (char-0 Bardzell path vs CS fallback)

The card requires the char-0-Bardzell-vs-CS decision to be made **now, with a benchmark**. All
numbers below were measured live this authoring round on the **dev-tip engine** (the venv resolves
`import quiverlab` to the main-repo `dev` src at
`/Users/marco/Desktop/HomologicalNetworks/quiverlab/src`, reporting `0.3.0` — the bump waits for
P80; my worktree is fast-forwarded to dev tip `1fb1cda`), single machine, one run each (wall-clock,
±20%). Route A = the shipped CS resolution over ℚ; Route B (lower bound) = the GF(p) Bardzell HH
engine (`engine.hh_engine.hochschild_homology_dims` + `BardzellResolution`, single prime 32003).

**HH² / HH³ at the MC-relevant degree (`top = 3`).**

| algebra | dim | CS/ℚ (s) | CS/GF(p) (s) | Bardzell/GF(p) (s) |
|---|---|---|---|---|
| `k[x]/(x⁵)` | 5 | 0.035 | 0.020 | 0.0002 |
| `k[x]/(x¹²)` | 12 | 0.280 | 0.112 | 0.002 |
| `k[x]/(x²⁰)` | 20 | 2.203 | 0.731 | 0.006 |
| `kZ₆/J²` (gentle, rad²=0) | 12 | 0.104 | 0.054 | 0.004 |
| `kZ₁₀/J²` | 20 | 0.493 | 0.161 | 0.009 |
| `kZ₁₀/J³` | 30 | 1.536 | 0.435 | 0.011 |
| `kZ₁₆/J²` | 32 | 2.411 | 0.723 | 0.013 |
| `kZ₂₀/J²` | 40 | 5.626 | 1.538 | 0.012 |

**High-degree divergence (`k[x]/(x⁸)`, top 3→8):** CS/ℚ stays ≈ 0.07–0.12 s to degree 8 (its CS
resolution is minimal for a truncated algebra); Bardzell/GF(p) ≈ 0.0004–0.001 s. For *large-dim*
algebras (e.g. the dim-220 Nakayama zoo members) CS is infeasible while Bardzell is trivial — but
those are **out of MC scope** (MC lives in degree ≤ 3 at oracle-scale dim).

**The obstruction bracket (`ℓ₂`, the expensive piece — RE-BENCHMARKED this fix round, H2).**
`QuantumCI(0)/ℚ` `gerstenhaber_brackets(top=3, engine="cs")` (dim 4) took **18.1 s** on this machine
vs **2.37 s** over GF(32003) (single run, wall-clock ±20%) — the CS homotopy-lifting cost for ℚ
arithmetic (≈ 8×). The **HH dims themselves are ≈ 0.1 s**; the bracket dominates. **The earlier
"≈ 22 s / 2.8 s" figure is superseded** by this machine's measured **18.1 s / 2.37 s** (a prior
reviewer measured 5.27 s for the same ℚ computation — machine-dependent; per the fix-round rule
**this machine's measurement is the reference**). **Implementation note:** the obstruction needs
only the **(2,2)** self-bracket `HH² × HH² → HH³` (or a single-cocycle `native_bracket`), **not** the
full `(p,q)` table — but at `top = 3` the **degree-2 block already dominates**, so the targeted (2,2)
build (**17.9 s** measured) is essentially the whole cost (≈ the full-table 18.1 s here, **not**
"markedly cheaper" as the earlier draft claimed). The targeted build's real value is that it computes
**only** the obstruction and can **early-out** on the first nonzero class; it is what
`obstruction_map` must call. See the next subsection for the cost-driver benchmark.

**Bracket-cost re-benchmark — the driver is HH-richness, NOT dimension (H1, this fix round).** The
critic's finding (a dim-4 QuantumCI full bracket table outran a dim-20 `kZ₁₀/J²`) is confirmed and
sharpened: the **obstruction bracket** — the targeted (2,2) self-bracket `HH² × HH² → HH³` that
`obstruction_map` calls (CS-native, cached liftings) — was measured live over ℚ on this machine
against `(dim, dim HH², dim HH³)`. `HH^•` dims are the cheap `hochschild_cohomology(3, engine="cs")`
(≈ 0.1–0.2 s); the timing is the obstruction quadratic form:

| algebra | dim | dim HH² | dim HH³ | (2,2) pairs | targeted (2,2) obstruction | full `(p,q)` table |
|---|---|---|---|---|---|---|
| `kZ₆/J²` (rad²=0, `HH•=[1,1,0,0]`) | 12 | 0 | 0 | 0 | **0.006 s** | 0.227 s |
| `kZ₁₀/J²` (rad²=0, `HH•=[1,1,0,0]`) | 20 | 0 | 0 | 0 | **0.024 s** | 2.03 s |
| `QuantumCI(0)` = `k⟨x,y⟩/(x²,y²,xy)` (`HH•=[2,2,3,5]`) | 4 | 3 | 5 | 6 | **17.9 s** | 18.1 s |
| `QuantumCI(−1)` = `k⟨x,y⟩/(x²,y²,xy−yx)` (`HH•=[4,4,5,6]`) | 4 | 5 | 6 | 15 | **43.5 s** | 43.8 s |
| `QuantumCI(0,a=2,b=4)` = `k⟨x,y⟩/(x²,y⁴,xy)` (`HH•=[2,4,9,13]`, **dim-8 HH-rich**) | 8 | 9 | 13 | 45 | **> 600 s (did not finish in a 10-min synchronous run)** | not reached |

(Single run each, wall-clock ±20%, this machine. The last row's obstruction was killed at the 10-min
cap — that timeout **is** the measurement: dim 8 but HH²=9/HH³=13 → 45 lifting-solve pairs on a
larger CS resolution.)

**Reading (binding).** Dimension is **not** the predictor. Dim-20 `kZ₁₀/J²` (HH²=0) finishes the
obstruction in **0.024 s**; dim-4 `QuantumCI(−1)` (HH²=5) takes **43.5 s** — ~1800× slower at 1/5 the
dimension; and dim-8 `QuantumCI(0,2,4)` (HH²=9) does not finish in 10 minutes. The cost tracks the
**(2,2) work ≈ O((dim HH²)²) homotopy-lifting solves, each scaled by the CS-resolution size** — i.e.
**HH²/HH³ richness × resolution size**. rad²=0 gentle algebras (`kZ_n/J²`, HH²=0) get the obstruction
for free. At `top = 3` the targeted (2,2) build ≈ the full table because the degree-2 block dominates.

**Sizing the `deformations` kind — its OWN cap, not P70's (H1).** P70's `hh1_lie` uses
`DEFAULT_MAXDIM = 48`, calibrated for a **dim-driven** cost law; the deformation obstruction's cost
law is **HH-richness-driven** (above), so P78 does **not** inherit it. Instead: (1) a plan-specific
**HARD backstop `DEFORM_MAXDIM = 32`** on `A.dim` (a coarse DoS guard, below P70's 48 — it bounds the
cheap `HH²` pre-probe to a few seconds; `kZ₁₆/J²` HH ≈ 2.4 s at dim 32); and (2) the estimator's
**real routing** first computes the cheap `HH²`/`HH³` dims (≈ 0.1 s) and sizes the obstruction bracket
off **`(dim HH²)² × CS-resolution size`**, routing off the instant tier whenever that product is high
— **regardless of dim**. So an HH-poor dim-40 `kZ₂₀/J²` stays instant while an HH-rich dim-4/dim-8
QuantumCI is queued. The dim-220 Nakayama webapp examples carry **no** `deformations` (Plan-35
omission precedent; unchanged).

**What a char-0 Bardzell *port* would cost/save (survey-grounded).** The Bardzell **combinatorial
layer** (`MonomialPresentation`: `associated_paths`, `left/right_decomposition`, the Bardzell/Anick
chains) is **genuinely field-free** — pure tuple-of-arrow-id combinatorics, **no prime baked in**.
The **only** arithmetic is the `int64` **±1** structure signs in `BardzellResolution.
differential_matrix`, "kept UN-REDUCED mod p" — i.e. it already emits **integer** coefficients.
So a char-0 Bardzell HH path is *nearly free* to obtain (reinterpret the ±1 differential over
ℤ ⊂ ℚ). **But it is not needed for the backbone:** CS-over-ℚ already clears every MC-relevant
computation (degree ≤ 3, oracle dim ≤ ~40, ≤ 6 s for HH; the obstruction bracket sized off the
instant tier), and the Bardzell speedup (100–1000×) is entirely at **high degree / large dim**
that MC does not use.

### DECISION (binding)

1. **Backbone = CS-over-ℚ.** `HH²`, `HH³`, `ℓ₂` (the native CS Gerstenhaber bracket, P51), the
   obstruction `[α,α]`, the DGLA order-by-order MC on `C(A)`, MC = 2-cocycles in the nilpotent
   regime, the presented `A_α` feedback, and the Ext-algebra handoff — **all on the shipped CS
   engine over ℚ**. **No general char-0 Bardzell HH port** (its only advantage is off-scope). This
   resolves the card's binary in favour of the **CS fallback** for everything the deformation
   functor needs, benchmark-justified.
2. **`ℓ₃` is the plan's one genuinely-new engine question, and it lives on the char-0 Bardzell
   complex** (CS cannot produce `ℓ_{≥3}`). It is **feasibility-gated** (Task 4): ship (a) the
   **rad²=0 dg-Lie certificate** (`ℓ_{≥3} ≡ 0`, exact, cheap — the field-free Bardzell substrate
   trivially decides `rad² = 0`; **unconditional, no L∞ machinery**), (b) the **induced-ℓ₂ ≡
   CS-bracket model-independence statement** (§4 theorem, cited; promoted to a live crossengine
   self-cert **only if** the spike builds the `B(A)` ℓ₂ adapter — M1 ruling, see Task 4 Step 2b′),
   and (c) a **char-0 Bardzell `ℓ₃` for monomial `A`** validated against RRB's explicit truncated
   computations — (b-as-computed) and (c) both *only if* the bounded feasibility spike (Task 4
   Step 1) lands. If the spike shows RRB's `ℓ₃` transfer formula is beyond one plan, **freeze v1 ℓ₃
   scope to (a) + the (b) theorem statement** (no computed `B(A)` ℓ₂ without the adapter) and defer
   the general `ℓ₃`/`ℓ₄` to the GUI-deferral ledger with this benchmark + the §4
   DGLA-equivalence as the honest justification. Either outcome is a complete, honest v1.

---

## Live-verified facts (all recomputed in the venv at spec time)

Two independent computations back each pin: (A) the shipped CS engine over **ℚ**
(`hochschild_cohomology(engine="cs")`, `gerstenhaber_brackets(engine="cs")`, `ext_algebra`,
`Quiver.algebra`), and (B) hand cross-checks (dimension flatness, obstruction symmetry). Every
number below reproduced this authoring round.

**Pin 1 — the classic anchor `k[x]/(x²) ⤳ k[x]/(x²−t)` as an MC/2-cocycle story.**
`HH•(k[x]/(x²))/ℚ = [2, 1, 1, 1]` (CS). `HH² = 1` is the **1-dimensional deformation tangent**;
its generator is the 2-cocycle `f` with `f(x,x) = 1` (the unit direction), whose formal
integration is `k[x]/(x²−t)`, **semisimple at `t ≠ 0`** (`x²−t = (x−√t)(x+√t)`, `≅ k × k` over the
closure). **Honest scope, live-confirmed:** this deformation is *not* expressible as a presented
`kQ/I` — the deformed relation `x² = t·1` carries a **length-0 (unit) term**, which the relation
**parser refuses at parse level with `RelationError`** (live: `Quiver([1],{"x":(1,1)}).algebra(
["x*x - 1"], field=QQ)` → `RelationError("term '1' has no arrows  [hint: relations live in the
arrow ideal]")`) — the build never even reaches the rad² admissibility gate. (`AdmissibilityError`/
`NotFiniteDimensionalError` are the refusal path for a *parseable* but non-admissible / infinite
radical direction, a distinct case.) So the semisimplification is an `HH²`/MC story on `C(A)`,
**not** a presented-`A_α` feedback (§5 honest-scope, RRRV admissible sub-locus).

**Pin 2 — a REAL obstruction `[α,α] ≠ 0 ∈ HH³` (the card's "construct one" — DONE).**
`A = QuantumCI(0) = k⟨x,y⟩/(x², y², xy)` (**monomial**, in RRB's scope; `rad² = (yx) ≠ 0` so
**honest L∞**, not dg-Lie; **not gentle** — it has loops). `HH•/ℚ = [2, 2, 3, 5]` (CS). The
**(2,2) Gerstenhaber self-bracket quadratic form** (`native_bracket` on the `cs/ℚ` `HH²` basis,
`top = 3`) has some nonzero class — the *specific* entries are a **basis-dependent illustration**,
**NOT a value to pin**:
```
# ILLUSTRATION ONLY (basis cs/ℚ, representative-dependent -- do not hard-pin in code/tests).
# Two DIFFERENT cs/ℚ bases seen live, both a genuine obstruction:
#  * one basis (author's, dev@1fb1cda): diagonal ALL zero; nonzero [a0,a2]=1, [a1,a2]=2  -> witness a0+a2
#  * this fix round's basis (this worktree): [a1,a1]=2 != 0 (a BASIS cocycle self-obstructs),
#    plus [a0,a2]=1, [a1,a2]=2; [a0,a0]=[a0,a1]=[a2,a2]=0
```
Either way `[α, α] ≠ 0` in `HH³` for some `α ∈ HH²` — a **genuine second-order obstruction**: this
infinitesimal deformation does not lift. **What is pinned is basis-independent ONLY:** (i)
`unobstructed is False`; (ii) the obstruction class is **nonzero** and lives in **`HH³`** (degree 3),
of the recorded `HH³` dimension (3, 5). **CORRECTED this fix round (Minor-2, live):** the earlier
draft pinned "*all* diagonal self-brackets vanish, so the obstruction is a combination" — that is
**basis-DEPENDENT and FALSE on this engine** (`[α₁, α₁] = 2 ≠ 0` here; live-verified), so it is **NOT
pinned**. The implementation must scan the **whole quadratic form** `Q(c) = Σ c_i c_j [α_i, α_j]` for
a nonzero class — not because "the diagonal vanishes" (it need not) but for the basis-independent
reason that **unobstructedness requires `[α,α] = 0` for *all* `α`** (Open risk 2). The **specific
witness string** (`"α₀ + α₂"`, `"a_1"`, …) and the **specific constants** (`1`, `2`) are
basis/representative-dependent and **must never be hard-pinned** (Plan-35 basis rule; Minor-2).

**Pin 3 — gentle nilpotent examples are UNOBSTRUCTED (MRRS Thm 5.4 in action).**
Every small gentle/nilpotent example computed has vanishing obstruction: `kZ₃/J²`, `kZ₄/J²`,
`kZ₅/J²` (rad²=0, dg-Lie) have `HH³ = 0` (`HH• = [1,1,0,0]`); `kZ₃/J³` has `HH• = [1,1,1,1]` with
all self-brackets `0`. This is the **honest record the card asks for**: *the small gentle examples
are unobstructed* — exactly the MRRS nilpotent-regime prediction (MC = 2-cocycles). The genuine
obstruction (Pin 2) needed a **monomial-with-loops** (non-gentle, non-nilpotent) algebra.

**Pin 4 — the presented `A_α` feedback loop, live.** The `QuantumCI(q) = k⟨x,y⟩/(x², y², xy + q·yx)`
family is a **flat presented deformation** (dim 4 for **every** `q`), with `HH` moving under the
deformation:

| `q` | relations | dim | `HH⁰..³/ℚ` |
|---|---|---|---|
| `0` (monomial base) | `x², y², xy` | 4 | `[2, 2, 3, 5]` |
| `−1` (commutative) | `x², y², xy − yx` | 4 | `[4, 4, 5, 6]` |
| `1` (anticommutative) | `x², y², xy + yx` | 4 | `[2, 4, 6, 8]` |
| `2` | `x², y², xy + 2yx` | 4 | `[2, 2, 1, 0]` |
| `−1/2` | `x², y², xy − ½yx` | 4 | `[2, 2, 1, 0]` |

Deforming the monomial base by hand —
`Quiver([1],{"x":(1,1),"y":(1,1)}).algebra(["x*x","y*y","x*y - q*y*x"])` — rebuilds a **flat**
(`dim A_α = 4 = dim A`) presented algebra whose `HH` jumped `[2,2,3,5] → [4,4,5,6]` at **`q = −1`**
(the commutative point `k[x,y]/(x²,y²)`; under quiverlab's `xy + q·yx` convention the commutative
member is `q = −1`, the anticommutative/exterior member is `q = +1` → `[2,4,6,8]`). This is the
`A_α` feedback: build the deformed presented algebra, re-certify flatness, recompute invariants —
**all live over ℚ**.

**Pin 5 — the Ext-algebra handoff (RRR), live.** `QuantumCI(q).ext_algebra(3)` (Plan 27) returns a
`YonedaPresentation` for every `q ∈ {0, −1, 1}` (`2 generators, 1 relation, koszul=True`) — the
`E(A_α)` handoff computes on the deformed member. (RRR describe `E(A_f)` vs `E(A)`; the pin is that
the handoff is a live computation on the deformed presented algebra.)

**What I could NOT live-verify (labelled to-be-verified-at-implementation):**
- **The char-0 Bardzell `ℓ₃` values** — RRB give explicit `ℓ_n` for truncated algebras; I verified
  the *substrate* is field-free and reusable over ℚ (survey) and the *structural* claims (rad²=0 ⇒
  dg-Lie; MC in degree 2) from the abstracts, but I did **not** implement RRB's ℓ₃ formula. Task 4
  Step 1 is the bounded feasibility spike; the truncated ℓ₃ pin is `# PIN` against RRB §(concrete
  computations).
- **The exact MRRS Thm-5.4 quiver hypotheses** — R13 says "gentle, no parallel arrows, no oriented
  cycles"; the abstract says "under some hypothesis on Q". `# PIN` the precise conditions against
  §5 before writing `is_l_infinity_nilpotent`.
- **RRRV's presentation of `A_α` in the multi-relation gentle case** — I verified the QuantumCI
  radical family and the single-relation hand build; the general RRRV presentation (alg-closed) is
  the `# PIN` at implementation.

---

## Scope gates (contractual; every boundary is a loud typed refusal)

| surface | scope | refusal |
|---|---|---|
| `infinitesimal_deformations` (`HH²`), `obstruction_map` (`[α,α] ∈ HH³`) | **any exact field** (the raw `HH`/bracket compute); presented `A` with a CS-admissible presentation; `A.dim ≤ budget` | oversize (`A.dim > budget`) → `QuiverlabError` with the CS cost note; presentation-less structure-constant `A` → loud (CS needs a presentation) |
| the **deformation interpretation** (`maurer_cartan`, `is_l_infinity_nilpotent`, `dg_lie_certificate`, `l3_bracket`, the `A_α` narrative) | **`char == 0` only** (RRB/MRRS are char-0 theorems; the L∞ `1/n!` need char 0; the DGLA `½` needs char ≠ 2) | `char > 0` → loud `QuiverlabError("formal-deformation / L∞ outputs need characteristic 0 …")`; the field-general `HH²`/`[α,α]` block is still returned with the char-0 caveat |
| `is_l_infinity_nilpotent` / MC = 2-cocycles simplification | **gentle** `A` satisfying the MRRS Thm-5.4 quiver hypotheses (`# PIN`) | outside the hypotheses → verdict `None` + the honest note "nilpotent regime not certified; the full MC equation applies" (never a silent MC = Z² claim) |
| `dg_lie_certificate` / `l3_bracket` | **monomial** `A` (RRB's scope), char 0 | non-monomial → `l3_bracket` refuses loudly ("RRB's B(A) L∞ is monomial"); `dg_lie_certificate` returns `None` (only rad²=0 monomial is certified dg-Lie) |
| `l3_bracket` beyond the feasibility-spike scope | monomial `A` within the **shipped ℓ₃ window** (truncated + the classes the spike validated); else deferred | out of window → loud `status="deferred"` refusal citing the ledger entry + the §4 DGLA-equivalence (never a fabricated ℓ₃) |
| `deformed_algebra(A, direction, t)` | **radical** cocycle direction (paths length ≥ 2) that keeps `A_α` **admissible + finite** | **unit direction → parse-level `RelationError`** (the length-0 term: "term '1' has no arrows", Pin 1); parseable-but-non-admissible / infinite radical direction → `AdmissibilityError` / `NotFiniteDimensionalError` — all relayed as a loud `QuiverlabError`; the honest note names the unit-direction / radical-collapse scope (Pin 1) |
| the **MC solution beyond first order** (non-nilpotent regime) | solved on the `C(A)` DGLA to a **certified truncation order** `≤ order` (the `order` **parameter** is the truncation limit — the native CS bracket, P51, runs at *any* degree past the bar window, so there is no "bracket window" bounding the reachable order here) | past the certified order → the report states the order reached and that higher orders are not certified (never claims a complete formal solution it did not verify) |
| all | `A.dim ≤ DEFORM_MAXDIM` (**the plan's OWN cap = 32**, NOT P70's `DEFAULT_MAXDIM = 48` — the deformation obstruction has a **different cost law**: it tracks **HH²/HH³ richness × resolution size**, not `A.dim`, so a dim cap is a coarse backstop and the real routing is off HH-richness — H1 re-benchmark) | `A.dim > cap` → honest `status="budget"`, never a partial guess; the dim-220 Nakayama webapp examples carry **no** `deformations` (the Plan-35 products-omission precedent, same bracket cost) |

Bracket/obstruction **constants** always carry `basis="cs/<domain>"` provenance; cross-engine
comparisons use only **basis-independent** data (dims, verdicts, obstruction rank, the flattened
class-map rank) — the Plan-35 rule (never compare raw constants across engines/bases).

---

## API surface (public via `import quiverlab`; exact only)

```python
# src/quiverlab/hochschild/deformations.py
def infinitesimal_deformations(A, *, engine="auto", max_cells=4_000_000) -> dict
def obstruction_map(A, *, engine="auto", max_cells=4_000_000) -> Obstruction
def maurer_cartan(A, *, order=2, engine="auto", max_cells=4_000_000) -> MCReport
def is_l_infinity_nilpotent(A) -> bool | None      # MRRS Thm 5.4 gate (char-0, gentle)
def dg_lie_certificate(A) -> bool | None           # rad^2 = 0 => l_{>=3} = 0 (monomial)
def l3_bracket(A, *, budget=...) -> "L3Report"     # feasibility-gated char-0 Bardzell l3
def deformed_algebra(A, direction, *, t="1") -> "Algebra"   # build the presented A_alpha
def deformation_structure(A, *, budget=..., engine="auto") -> Deformations
def deformations_block(A, *, budget=...) -> dict   # the runner block

@dataclass(frozen=True)
class Obstruction:
    hh2_dim: int
    hh3_dim: int
    basis: str                 # "cs/<domain>"
    unobstructed: bool          # the quadratic map alpha |-> [alpha,alpha] is 0 on all of HH^2
    witness: str | None         # a basis-DEPENDENT display of an obstructed direction, e.g.
                                #   "a_0 + a_2" (Pin 2); NEVER hard-pinned by tests — tests assert
                                #   unobstructed is False + [witness,witness] != 0 in HH^3, not the string
    obstruction_constants: list # symmetric (2,2) self-bracket classes, EXACT strings (basis-dep;
                                #   provenance-tagged, never compared across bases — Plan-35 rule)
    characteristic: int
    references: list

@dataclass(frozen=True)
class Deformations:
    characteristic: int
    hh2_dim: int                # infinitesimal deformations (tangent)
    hh3_dim: int
    basis: str
    unobstructed: bool | None   # None over char p (field-general HH only)
    obstruction_witness: str | None
    nilpotent_regime: bool | None   # MRRS Thm 5.4 verdict (None if hypotheses not met/char p)
    mc_description: str          # "MC = Z^2 (all 2-cocycles, nilpotent regime)" | "obstructed ..." 
    mc_order_certified: int      # highest order the DGLA MC was solved+certified (>=1)
    dg_lie: bool | None          # rad^2=0 certificate (monomial); None otherwise
    l3_status: str               # "dg_lie" | "computed" | "deferred" | "n/a (char p / non-monomial)"
    a_alpha: dict | None         # presented A_alpha for a canonical radical direction:
                                 #   {"quiver": ..., "relations": [...], "t": "1", "flat": bool}
    ext_algebra_summary: str | None   # E(A_alpha) handoff (Plan 27), when a radical A_alpha exists
    base_change_note: str | None # RRRV presentation is alg-closed; honest wording (P62/P70 precedent)
    char0_note: str | None
    window_note: str | None      # the certified MC truncation order (the `order` parameter limit)
    status: str                  # "complete" | "budget" | "unsupported"
    note: str
    references: list

# Algebra delegates (core/algebra.py, thin lazy-import beside gerstenhaber_brackets):
Algebra.deformation_structure(budget=..., engine="auto") -> Deformations
Algebra.obstruction_map(engine="auto") -> Obstruction
Algebra.deformed_algebra(direction, t="1") -> Algebra
```

`Obstruction`/`Deformations`/`MCReport`/`L3Report` are frozen data reports (like Plan-35's
`HHProducts`); `__bool__` is NOT defined.

---

## Global Constraints

- Python always `.venv/bin/python`; tests `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2
  .venv/bin/python -m pytest -q ...`.
- **Reuse, do not reinvent (verified live in `dev` at spec time):**
  - **`A.hochschild_cohomology(top, engine="cs")`** / **`hochschild_homology`** — `HH²`/`HH³` over
    any exact Domain (ℚ included); self-certifies (`assert_dd_zero` + `assert_order_condition`).
  - **`A.gerstenhaber_brackets(top, engine="cs")`** (Plan 35 / **P51**) and the direct
    `resolutions_cs.bracket.native_bracket(res, f, p, g, q)` — the native CS Gerstenhaber bracket
    for any `p, q ≥ 1` over any Domain, **past the bar window** (`ℓ₂`; the obstruction `[α,α]` is
    the `(2,2)` case). Use a **targeted (2,2)** build for the obstruction — but note (re-benchmark,
    H1) the targeted (2,2) ≈ the full table at `top = 3` (both ≈ 18 s on `QuantumCI(0)/ℚ`, dim 4):
    the win is early-out + computing only the obstruction, not a cheaper block. **Cost tracks HH²/HH³
    richness × resolution size, not `A.dim`** (see the bracket-cost re-benchmark table).
  - **`resolutions_cs.homology.cs_hh_basis(A, n, side)`** — representative (co)cycles (for a chosen
    direction / the `A_α` build).
  - **`Quiver.algebra(relations=[...], field=...)`** (`field=`, **not** `domain=`; `QQ` from
    `quiverlab.fields`), **`groebner.build_reduction_system` / `reduction_system_of(A)`**, and the
    loud `RelationError` (parse-level: the unit's length-0 term) / `AdmissibilityError` /
    `NotFiniteDimensionalError` (`quiverlab.errors`) — the `A_α` build + re-certification (there is
    **no** `is_admissible` predicate; certify by building + catching all three).
  - **`A.ext_algebra(top)`** (Plan 27 `YonedaPresentation`) — the `E(A_α)` handoff.
  - **`engine.resolutions_bardzell.MonomialPresentation`** — the **field-free** Bardzell
    combinatorics (associated paths, left/right decompositions); the ℓ₃ adapter reuses it verbatim
    and reinterprets `BardzellResolution.differential_matrix`'s integer ±1 entries over ℚ.
  - **`A.is_gentle()`**, `A.quiver`, `A.relations` (parsed; `r.is_monomial`), `A.dim`, `A.domain`,
    `A.domain.characteristic`, `A.radical`-style probes for `rad² = 0`.
  - **`quiverlab.fields.linalg`**: `nullspace`, `rank`, `rref`, `solve` (exact over any Domain).
- **The char-0 gate is `A.domain.characteristic == 0`** for the deformation/L∞ interpretation (the
  papers' scope); the field-general `HH²`/`[α,α]` block stays live in every characteristic behind
  the caveat. Document at the gate; reference §1/§2 rationale, do not re-argue it.
- **No floats in `src/`.** Cocycle/obstruction constants are exact-`Domain` strings; dims/orders
  `int`; verdicts `bool`/`str`. Client-side numeric conversion only (`gui.js`, exempt).
- **Composition is left-to-right**; `δ` and `[·,·]` respect `A.multiply`. All refusals are
  `QuiverlabError` (relaying `RelationError`/`AdmissibilityError`/`NotFiniteDimensionalError` as
  typed 4xx in the webapp, never a 500).
- **Plan-32 markers** (orthogonal; `oracle_*` FORBIDDEN in `tests/{webapp,gui,hpc}`):
  - `oracle_literature`: `HH²(k[x]/x²) = 1` deformation tangent (Pin 1); the QuantumCI(q) flat
    family `HH` table (Pin 4); rad²=0 ⇒ dg-Lie (RRB) on `kZ_n/J²`; **the RRB truncated ℓ₃ pin**
    (`# PIN`, if Task 4 lands); the MRRS nilpotent-regime = unobstructed on gentle examples (Pin 3).
  - `oracle_crossengine`: `dim HH²/HH³ == A.hochschild_cohomology` (CS ≡ bar in-window degreewise);
    the obstruction `[α,α]` on **basis-independent** data (rank / class-nonvanishing) agrees between
    `engine="cs"` and the in-window transported bracket (GF(p)); `B(A)`-induced ℓ₂ ≡ the CS bracket
    on `HH` **only if the Task-4 spike lands** (the ℓ₃-adapter self-cert — without the adapter there
    is no `B(A)` ℓ₂ to compare, so this reduces to the §4 model-independence *theorem*, M1 ruling);
    `deformed_algebra` `HH` ≡ the direct family member
    (QuantumCI(q) built two ways).
  - `oracle_selfcert`: `[α,α]` is a cocycle (`δ[α,α] = 0`); the obstruction quadratic form is
    symmetric (degree-2 antisymmetry — a basis-independent property of `Q`, unlike the diagonal
    pattern); the DGLA MC recursion `δμ_k = −½ Σ[μ_i,μ_j]` is solvable iff the RHS class vanishes
    (the obstruction gate); `dim A_α = dim A` flatness; the char-p loud gate; the unit-direction
    `RelationError` / oversize loud refusals; `dg_lie_certificate is True` on rad²=0 monomial
    (`kZ_n/J²`) and `is False` on `QuantumCI(0)` (rad²≠0) — a decidable structural check (and on the
    rad²=0 members `HH²=0`, so the obstruction is trivially `unobstructed`).
  - `qpa`: **QPA/GAP have NO Hochschild deformation / L∞ / Maurer–Cartan surface** (a live
    `NamesGVars()` sweep confirms; an honest skip that **FAILS if it ever appears** — the Plan-35
    precedent). The nearest genuine QPA cross-check is **`dim HH²`/`dim HH³`** of `A` and of a built
    `A_α` via QPA's `HochschildCohomology` where QPA implements it, and the **Ext-algebra of `A_α`**
    via QPA's `ExtAlgebra*` (Plan-27 bridge) — i.e. the *inputs* to the deformation theory, not the
    deformation theory itself (honest-scope on the verification page).
- **Mid-merge-train counts.** The verification task recounts the oracle-class table at merge
  (`tests/release/test_oracle_classes.py`; paste live numbers, claim only this plan's deltas).
- Every merge updates `docs/verification.md` + adds citations to `references.bib` + `registry.py`
  (`bibtex()` hard-fails if the two disagree). Conventional commits; green at every commit.

**No-red-commits discipline.** Each Task ships a self-contained green slice: the infinitesimal +
obstruction surface before the MC report before the `A_α` feedback before the (gated) ℓ₃ before the
GUI. Char-0-only fields default to `None`/`"n/a"` so intermediate task boundaries are green.

---

## Task 0 — citations (do first; later tasks reference the keys)

**Files:** `src/quiverlab/citations/references.bib` + `registry.py`.

- [ ] **Step 1: Add BibTeX (verified fields only; `# PIN` the rest, never guess).**
```bibtex
@article{rrb_linfty_bardzell,
  author  = {Redondo, Mar\'ia Julia and Rossi Bertone, Fiorela},
  title   = {{$L_\infty$}-structure on {B}ardzell's complex for monomial algebras},
  journal = {J. Pure Appl. Algebra}, volume = {226}, number = {5}, year = {2022},
  pages   = {Paper No. 106935}, note = {arXiv:2008.08122}}          % article no./DOI # PIN
@article{mrrs_mc_gentle,
  author  = {M\"uller, Monique and Redondo, Mar\'ia Julia and Rossi Bertone, Fiorela and
             Suarez, Pamela},
  title   = {{M}aurer--{C}artan equation for gentle algebras},
  journal = {Comm. Algebra}, volume = {53}, number = {9}, pages = {3984--4007}, year = {2025},
  note    = {arXiv:2309.02582}}                                     % DOI 10.1080/00927872.2025.2473029 (confirmed)
@misc{rrrv_morita_deform,
  author = {Redondo, Mar\'ia Julia and Rom\'an, Lucrecia and Rossi Bertone, Fiorela and
            Verdecchia, Melina},
  title  = {{M}orita invariance for infinitesimal deformations},
  year   = {2020}, note = {arXiv:2003.10366}}                       % journal ref # PIN
@misc{rrr_ext_deform,
  author = {Redondo, Mar\'ia Julia and Rom\'an, Lucrecia and Rossi Bertone, Fiorela},
  title  = {The {E}xt-algebra for infinitesimal deformations},
  year   = {2022}, note = {arXiv:2202.01199}}                       % THREE authors; journal ref # PIN
@article{chouhy_degeneration,
  author  = {Chouhy, Sergio},
  title   = {On geometric degenerations and {G}erstenhaber formal deformations},
  journal = {Bull. Lond. Math. Soc.}, year = {2019}, note = {arXiv:1708.02933,
            DOI 10.1112/blms.12277}}                                % vol/pages/year # PIN
```
- [ ] **Step 2:** add the `registry.py` entries with annotations naming what each grounds
  (`rrb_linfty_bardzell` = L∞ on B(A), MC in degree 2, rad²=0 dg-Lie; `mrrs_mc_gentle` =
  nilpotent-regime MC = 2-cocycles; `rrrv_morita_deform` = presented `A_α` by quiver+relations;
  `rrr_ext_deform` = Ext-algebra of `A_α`; `chouhy_degeneration` = degeneration ↔ formal
  deformation). Verify `tests/citations/test_bib_structure.py` green.
- [ ] **Step 3:** commit `docs(citations): P78 deformation/L∞ references (RRB L∞-Bardzell, MRRS MC-gentle, RRRV Morita, RRR Ext-algebra [3 authors], Chouhy degeneration)`.

---

## Task 1 — infinitesimal deformations + the obstruction map (field-general HH / CS ℓ₂)

**Files:** create `src/quiverlab/hochschild/deformations.py`; test
`tests/hochschild/test_deformations.py` (deep bucket).

- [ ] **Step 1: Failing tests.**
  - `infinitesimal_deformations(k[x]/x²)/ℚ` reports `hh2_dim == 1` (Pin 1); `== A.hochschild_
    cohomology(2)` cross-check (crossengine).
  - `obstruction_map(QuantumCI(0))/ℚ` (monomial): `hh2_dim == 3`, `hh3_dim == 5`, `unobstructed is
    False`, and the returned `witness` yields `[witness, witness] ≠ 0` in `HH³` (Pin 2). **Pin only
    basis-independent facts** — assert `unobstructed is False`, and that the obstruction class the
    witness produces is nonzero and of degree 3. **Do NOT assert whether the diagonal vanishes** (it
    is basis-dependent — live: `[α₁,α₁] ≠ 0` on this engine, `= 0` on the author's — Minor-2); the
    reason `obstruction_map` scans the **whole quadratic form** `Q(c)=Σ c_i c_j [α_i,α_j]` is that
    unobstructedness needs `[α,α]=0` for *all* `α`, not that the diagonal vanishes (Open risk 2).
    **Do NOT assert the literal witness string or the raw constants** (basis/representative-dependent
    — Minor-2 ruling, Plan-35 basis rule). (literature + selfcert)
  - `obstruction_map` on `kZ₃/J³`, `kZ₄/J²`: `unobstructed is True` (Pin 3, nilpotent-regime).
  - `[α,α]` is a cocycle: `δ` of the returned obstruction chain is `0` (selfcert).
  - char-p: `obstruction_map(QuantumCI(0)/GF(5))` returns the field-general `HH²`/`[α,α]` block with
    a `char0_note` (the deformation *interpretation* is char-0 gated), not a crash.
- [ ] **Step 2: Implement** `infinitesimal_deformations` (thin over `hochschild_cohomology(2)`) and
  `obstruction_map` (a **targeted (2,2)** self-bracket: build `res = ChouhySolotarResolution`, get
  `HH²`/`HH³` reps via `cs_hh_basis`, call `native_bracket(res, r_i, 2, r_j, 2)` for the symmetric
  pairs, reduce to `HH³` classes; scan the quadratic form `Q(c) = Σ c_i c_j [α_i,α_j]` for a nonzero
  class → the witness). Provenance `basis="cs/<domain>"`. Char-0 caveat wired.
- [ ] **Step 3:** run `tests/hochschild/test_deformations.py -q`; **Step 4:** commit
  `feat(hochschild): P78 infinitesimal deformations HH^2 + obstruction map alpha|->[alpha,alpha] in HH^3 (CS native bracket, char-0 gated interpretation)`.

---

## Task 2 — the Maurer–Cartan report + the nilpotent-regime gate (char 0)

**Files:** `deformations.py`; `tests/hochschild/test_deformations_mc.py`.

- [ ] **Step 1: Failing tests.**
  - `is_l_infinity_nilpotent` on the pinned gentle examples returns `True` (within the `# PIN`
    MRRS hypotheses) and `maurer_cartan` reports `mc_description == "MC = Z^2 …"` with
    `nilpotent_regime is True` (literature).
  - `maurer_cartan(QuantumCI(0))/ℚ` (non-nilpotent): `nilpotent_regime is not True`,
    `mc_order_certified >= 1`, and the report exposes the obstruction (the order-2 lift **fails**
    for the witness direction — `δμ₂ = −½[α,α]` unsolvable because `[α,α] ≠ 0`) while the
    QuantumCI(q) direction **does** lift (selfcert: the DGLA recursion is solvable iff the RHS class
    vanishes).
  - char-p: `maurer_cartan(.../GF(5))` → loud `QuiverlabError` (deformation interpretation is char-0)
    while `infinitesimal_deformations`/`obstruction_map` still return their field-general block.
- [ ] **Step 2: Implement** `is_l_infinity_nilpotent` (char-0 + `is_gentle` + the `# PIN` MRRS
  quiver conditions — no parallel arrows, no oriented cycles, verified confluent-against §5) and
  `maurer_cartan`: in the nilpotent regime → `MC = Z²` (report `dim Z²`, "all 2-cocycles integrate,
  unobstructed"); else solve the **DGLA recursion on `C(A)`** order by order to `order` (default 2 —
  the first-order obstruction), each step a `solve` for `μ_k` against `−½ Σ[μ_i,μ_j]` (loud/record
  when the class is nonzero → obstructed), recording `mc_order_certified` and the `window_note` (the
  `order` **parameter** is the truncation limit — the native CS bracket has no degree window, so the
  reachable order is bounded by `order` and by `HH^{k+1}` staying computable, **not** by any bracket
  window). **Never claim a complete formal solution past the certified order.**
- [ ] **Step 3:** run; **Step 4:** commit `feat(hochschild): P78 Maurer-Cartan report -- nilpotent-regime MC=Z^2 gate (MRRS Thm 5.4) + order-by-order DGLA MC on C(A) with the [alpha,alpha] obstruction gate (char-0)`.

---

## Task 3 — the presented deformed algebra `A_α` + Ext-algebra handoff (feedback loop)

**Files:** `deformations.py`; `tests/hochschild/test_deformed_algebra.py`.

- [ ] **Step 1: Failing tests.**
  - `deformed_algebra(QuantumCI(0), direction=<yx>, t="1")` builds `k⟨x,y⟩/(x²,y²,xy − yx)` (the
    RRRV radical direction), `dim == 4` (**flat**), and its `HH` equals `QuantumCI(-1)` built
    directly (crossengine — the two-ways build, Pin 4).
  - a **unit** direction (`k[x]/x²`, `f(x,x)=1`) → loud `QuiverlabError` relaying **`RelationError`**
    (parse-level: the length-0 term "has no arrows", Pin 1), with the radical-collapse note; a
    parseable-but-non-admissible / infinite radical direction relays `AdmissibilityError` /
    `NotFiniteDimensionalError` instead (the distinct case).
  - `deformed_algebra(...).ext_algebra(3)` returns a `YonedaPresentation` (Pin 5, the RRR handoff);
    the `Deformations.ext_algebra_summary` is populated for a radical `A_α`.
  - a deformation making the algebra infinite-dimensional → loud `NotFiniteDimensionalError` relay.
- [ ] **Step 2: Implement** `deformed_algebra(A, direction, t)`: interpret `direction` as a radical
  2-cocycle (a target relation + a same/higher-degree radical perturbation, from `cs_hh_basis` or an
  explicit path spec), reconstruct the relation strings from `A.quiver`/`A.relations`, build
  `Quiver.algebra(deformed_relations, field=A.domain)`, and **re-certify**: catch **`RelationError`**
  (the parse-level unit-direction refusal — a length-0 term never reaches admissibility) **and**
  `AdmissibilityError`/`NotFiniteDimensionalError` (loud), assert `dim A_α == dim A` (flatness;
  `status` records
  a non-flat jump honestly). Populate `Deformations.a_alpha` (quiver + deformed relations + `t` +
  `flat`) and `ext_algebra_summary` for a **canonical** radical direction; `base_change_note` states
  RRRV's presentation is over algebraically closed `k` (honest wording, P62/P70 precedent).
- [ ] **Step 3:** run; **Step 4:** commit `feat(hochschild): P78 presented A_alpha feedback -- build+re-certify the deformed kQ/I (flatness dim A_alpha==dim A), loud on unit/non-admissible directions, Ext-algebra handoff (RRR)`.

---

## Task 4 — the `B(A)[1]` L∞ companion: rad²=0 dg-Lie + FEASIBILITY-GATED char-0 Bardzell ℓ₃

**Files:** `deformations.py` (+ possibly a small `hochschild/bardzell_linfty.py` adapter);
`tests/hochschild/test_linfty_l3.py`.

**Step 1 is a bounded FEASIBILITY SPIKE (the P79 precedent — scope frozen after the spike).**

- [ ] **Step 1 (spike, timeboxed):** implement RRB's ℓ₃ formula for the **truncated** algebra
  `k[x]/(xᵃ)` on a **char-0 Bardzell complex** — reuse the field-free `MonomialPresentation`
  combinatorics and reinterpret `BardzellResolution.differential_matrix`'s integer ±1 entries over
  ℚ (survey: the substrate is reusable verbatim). Validate `ℓ₃` against RRB's **explicit truncated
  computation** (`# PIN` the paper values). **Decision gate:** if it lands within the timebox →
  ship the monomial ℓ₃ (Step 2a) **and** promote the induced-ℓ₂ agreement to a live computed
  crossengine self-cert (Step 2b′, since the adapter now exists); if RRB's transfer formula is
  beyond one plan → **freeze v1 ℓ₃ scope** to the unconditional `dg_lie_certificate` (Step 2b) +
  the induced-ℓ₂ ≡ CS-bracket **model-independence theorem statement** (Step 2b′ — *not* a computed
  self-cert, since there is no `B(A)` ℓ₂ adapter to compare) and write the GUI/engine-ledger
  deferral for the general ℓ₃/ℓ₄ (Step 2c). **Record the outcome in this doc's Change log.**
- [ ] **Step 2a (if the spike lands): `l3_bracket`** for monomial `A` in the validated window —
  `ℓ₃` on `B(A)[1]` over ℚ, byte-reproducible (the CS-diagonal canonicalization precedent:
  `solve` + `reduce_mod_nullspace`). Oracles: RRB truncated pin (literature); `ℓ₃ ≡ 0` on rad²=0
  monomial (agrees with the certificate); honest loud refusal outside the window (never a fabricated
  value). **Do NOT claim `ℓ₄ = 0`** — report it as "not computed / may be nonzero" (honest L∞).
- [ ] **Step 2b (always): `dg_lie_certificate`** — monomial + `rad² = 0` ⇒ `ℓ_{≥3} ≡ 0` (RRB;
  decidable from the field-free Bardzell substrate, **no L∞ machinery needed**). Anchor: `kZ_n/J²`
  (rad²=0) ⇒ `dg_lie is True`; `QuantumCI(0)` (rad²=(yx)≠0) ⇒ `dg_lie is False`. **This certificate
  ships unconditionally.**
- [ ] **Step 2b′ (M1 ruling — the induced-ℓ₂ agreement):** the "`B(A)` ℓ₂ induces the CS
  Gerstenhaber bracket on `HH`" claim is, **unconditionally**, the §4 **model-independence theorem
  statement** (L∞-quasi-iso invariance — `rrb_linfty_bardzell` + §4), *cited, not computed*. It is
  promoted to a **live `oracle_crossengine` self-cert only if the Step-1 spike lands** — because
  computing `B(A)`'s ℓ₂ requires the char-0 Bardzell L∞ adapter, which is exactly the spike's
  deliverable. **Without the spike there is no `B(A)` ℓ₂ implementation to compare against**, so the
  crossengine test would be vacuous; the honest form in that case is the theorem statement in the
  report/verification page, never a trivial pass. (This is why M1's "vacuous-or-deferred" concern is
  resolved by gating the *computed* oracle behind the spike while keeping the *theorem* stated
  unconditionally.)
- [ ] **Step 2c (if the spike does NOT land): ledger deferral** — write the general-`ℓ₃`/`ℓ₄`
  deferral into `docs/plans/DEEPER-ENGINES-BACKLOG.md` (the GUI-deferral ledger P80 reconciles) with
  the benchmark + §4 DGLA-equivalence justification; `l3_bracket` returns `status="deferred"` with
  the loud honest note. **This is a complete, honest v1 outcome, not a failure.**
- [ ] **Step 3:** run `tests/hochschild/test_linfty_l3.py -q`; **Step 4:** commit
  `feat(hochschild): P78 B(A)[1] L-infinity companion -- rad^2=0 dg-Lie certificate (unconditional) + induced-l2 model-independence (theorem statement; computed self-cert iff the spike ships the B(A) l2 adapter) [+ char-0 Bardzell l3 for monomial A, per the feasibility spike]`.

---

## Task 5 — the `deformations` algebra-only GUI kind (all three tiers, i18n ×4)

An algebra-only budget-carrying kind (the `hh1_lie`/`tau_tilting` precedent — the budget caps
`A.dim` for the CS bracket). **v1 is report-only** (displays `A_α`; the GUI-adopt flow is a ledger
entry). GUI touchpoints follow `hh1_lie` (P70 Task 6).

> **Implemented-first dependency (Minor-5 ruling).** The `hh1_lie` GUI kind is **P70's** deliverable
> and is **NOT on `dev`** yet — only P70's *doc* is committed. This task **copies a pattern that must
> exist first**: the grammar-parse branch, `_dispatch` branch, `_snip`, the `gui.js` checkbox
> plumbing, and the estimator sizing are all "follow the `hh1_lie` precedent", which presupposes
> **P70's code is implemented and merged to `dev`**. So Task 5 is **sequenced after P70/P71 land**
> (Acceptance §7). P78's own *compute* (Tasks 1–4) consumes only P51 (native bracket) and P27
> (`ext_algebra`), both already on `dev`, and can be built in parallel; only this GUI task carries
> the P70 code dependency.

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (grammar parse — the `hh1_lie`/`tau_tilting` budget-suffix
  branch: `deformations` or `deformations:32` — + `_dispatch` branch + `_snip`),
  `docs/gui/runner.py` (the Pyodide twin — matching grammar + dispatch + `_snip` + ETA,
  shape-identical), **`webapp/server/schema.py`** (the third grammar-parse site), `docs/gui/gui.js`
  + `webapp/static/gui/gui.js` (checkbox `qlgui-deformations`, `S.ids`, push-list, `renderBlock`,
  `scheduleProbe`), `webapp/templates/index.html` (one checkbox), `webapp/server/i18n/{en,es,fr,zh}.
  json` (all four locales), `src/quiverlab/trace/results_html.py` (`_HEADINGS` + a render branch),
  `tests/webapp/_runner_goldens.json` + `test_runner_delegation.py` (ONE golden).
- Test: `tests/webapp/test_deformations_kind_p78.py`, `tests/gui/test_deformations_runner_twin_p78.py`.

**Block shape** (one shared `deformations_block(A, budget)`):
```python
{"kind": "deformations", "characteristic": int,
 "hh2_dim": int, "hh3_dim": int, "unobstructed": bool|None, "obstruction_witness": str|None,
 "nilpotent_regime": bool|None, "mc_description": str, "mc_order_certified": int,
 "dg_lie": bool|None, "l3_status": str,
 "a_alpha": {"quiver": ..., "relations": [...], "t": "1", "flat": bool}|None,
 "ext_algebra_summary": str|None,
 "base_change_note": str|None, "char0_note": str|None, "window_note": str|None,
 "status": str, "note": str|None,
 "references": ["rrb_linfty_bardzell","mrrs_mc_gentle","rrrv_morita_deform","rrr_ext_deform","chouhy_degeneration"]}
# refusal (char p interpretation / oversize budget / presentation-less) -> {"error": msg, "references":[...]}, never a 500.
```
The raw bracket **constants are NOT shipped in the block** (they explode and are basis-dependent);
the block reports `HH²`/`HH³` dims, the obstruction verdict + witness, the nilpotent verdict, the
dg-Lie certificate + ℓ₃ status, and the presented `A_α` (quiver + deformed relations) for a
canonical radical direction. The worked-steps report renders the deformation tangent, the
obstruction line (`[α,α]` unobstructed / witness), the MC description, and the `A_α` presentation
+ its Ext-algebra summary.

- [ ] **Step 1: Failing cross-runner tests** (unmarked — extras-gated dir): `QuantumCI(0)/ℚ` block
  has `hh2_dim==3`, `hh3_dim==5`, `unobstructed is False`, `"rrb_linfty_bardzell"` in the keys;
  `k[x]/(x²)/ℚ` block has `hh2_dim==1`; a `kZ₄/J²` block has `dg_lie is True`; a `.../GF(5)` block
  has a `char0_note`; twin parity (`json.dumps(sort_keys=True)` equality across `hpc/spec.py` and
  `docs/gui/runner.py`).
- [ ] **Step 2: Implement** the grammar parse in **all three** sites (budget suffix, default the
  **plan's OWN `DEFORM_MAXDIM = 32`**, NOT P70's `DEFAULT_MAXDIM = 48` — different cost law, H1;
  otherwise follow the `hh1_lie` branch **structurally** once P70 lands — see the implemented-first
  note), the `_dispatch` branch (`from quiverlab.hochschild.deformations import deformations_block`),
  the twin, the `gui.js` checkbox + `scheduleProbe` ETA (**estimator — measured, H1**: the
  obstruction bracket over ℚ is the driver — **18.1 s at dim 4 top 3** on this machine, and the cost
  tracks **HH²/HH³ richness × CS-resolution size, NOT `A.dim`** — dim-20 `kZ₁₀/J²` with HH²=0 is
  0.024 s while dim-4 `QuantumCI(−1)` with HH²=5 is 43.5 s. So the estimator **pre-computes the cheap
  `HH²`/`HH³` dims (≈ 0.1 s) and sizes the bracket off `(dim HH²)² × resolution size`**, routing off
  the instant tier when that product is high — regardless of dim; the dim-220 Nakayama examples get
  an honest `status="budget"` — no `deformations`, the Plan-35 omission precedent), the `_snip`
  recipe (`"deformations": "A.deformation_structure()"`), `_HEADINGS["deformations"] = "Formal
  deformations (L∞ / Maurer–Cartan)"` + the render branch.
- [ ] **Step 3: Add ONE golden** (`deformations_quantumci0`) to `_runner_goldens.json`; note it in
  `test_runner_delegation.py`'s change-log docstring; confirm existing goldens byte-identical first.
  **Canonical-key stability:** schema-v1 algebra-only, no `module` block → canonicalizes unchanged.
- [ ] **Step 4: GUI-deferral ledger entry** — the `A_α`-**adopt** flow (load the deformed algebra
  back onto the canvas as a fresh input; schema-heavy — a new quiver+relations spec the canvas
  ingests) is written into `docs/plans/DEEPER-ENGINES-BACKLOG.md` for P80 with Marco's sign-off note.
  v1 ships the presented `A_α` as **display-only**.
- [ ] **Step 5:** run `tests/webapp/test_deformations_kind_p78.py tests/webapp/test_runner_
  delegation.py tests/gui/test_deformations_runner_twin_p78.py tests/webapp/test_i18n.py tests/hpc -q`.
- [ ] **Step 6:** commit `feat(gui,webapp,hpc): deformations algebra-only kind (Plan 78) -- HH^2/obstruction/MC/dg-Lie + display-only A_alpha, both runners byte-identical, grammar in 3 sites, i18n x4, one golden; adopt-flow ledgered`.

**i18n keys — ALL FOUR LOCALES** (mirror `pick.kind.hh1_lie`):

| key | en | es | fr | zh |
|---|---|---|---|---|
| `pick.kind.deformations` | Formal deformations (L∞) | Deformaciones formales (L∞) | Déformations formelles (L∞) | 形变 (L∞) |
| `inv.deformations` | Deformations: HH², obstruction, Maurer–Cartan | Deformaciones: HH², obstrucción, Maurer–Cartan | Déformations : HH², obstruction, Maurer–Cartan | 形变：HH²、障碍、Maurer–Cartan |
| `block.deformations.title` | Formal deformations (L∞ / Maurer–Cartan) | Deformaciones formales | Déformations formelles | 形变 (L∞ / Maurer–Cartan) |
| `block.deformations.infinitesimal` | Infinitesimal deformations (HH²) | Deformaciones infinitesimales (HH²) | Déformations infinitésimales (HH²) | 无穷小形变 (HH²) |
| `block.deformations.obstruction` | Obstruction [α,α] ∈ HH³ | Obstrucción [α,α] ∈ HH³ | Obstruction [α,α] ∈ HH³ | 障碍 [α,α] ∈ HH³ |
| `block.deformations.unobstructed` | Unobstructed | Sin obstrucción | Sans obstruction | 无障碍 |
| `block.deformations.mc` | Maurer–Cartan description | Descripción de Maurer–Cartan | Description de Maurer–Cartan | Maurer–Cartan 描述 |
| `block.deformations.dglie` | rad²=0 ⇒ dg-Lie (ℓ≥3 = 0) | rad²=0 ⇒ dg-Lie | rad²=0 ⇒ dg-Lie | rad²=0 ⇒ dg-Lie |
| `block.deformations.aalpha` | Deformed algebra A_α | Álgebra deformada A_α | Algèbre déformée A_α | 形变代数 A_α |
| `block.deformations.charp` | Deformation/L∞ outputs need characteristic 0 | Requieren característica 0 | Exigent la caractéristique 0 | 需要特征 0 |

---

## Task 6 — verification page, README, suite gate

**Files:** `docs/verification.md`; `README.md`; `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick
P78); existing release gates.

- [ ] **Step 1: Verification page.** Add the P78 subsystem rows:
  - `hochschild/deformations.py` — `oracle_literature` (`HH²(k[x]/x²)=1`; the QuantumCI(q) flat
    family HH table; rad²=0 ⇒ dg-Lie on `kZ_n/J²`; the RRB truncated ℓ₃ pin **if Task 4 landed**;
    the MRRS nilpotent-regime = unobstructed on gentle); `oracle_crossengine`
    (`dim HH²/HH³ == hochschild_cohomology`; the obstruction on basis-independent data agrees
    CS ≡ transported in-window; `deformed_algebra` HH ≡ the direct family member; induced-ℓ₂ ≡ CS
    bracket **only if the Task-4 spike lands** — else it is the §4 model-independence theorem, not a
    computed test, M1); `oracle_selfcert` (`[α,α]` a cocycle; obstruction symmetric; DGLA MC recursion gate;
    flatness `dim A_α = dim A`; the char-p / unit-direction / oversize loud refusals); `qpa`
    (honest skip — no QPA deformation surface; the input HH²/HH³ + `E(A_α)` bridge only).
  - **Honest-scope entries (binding):**
    (a) **All deformation/L∞ outputs are char-0** (RRB/MRRS are char-0 theorems; the `1/n!`/`½`);
    the field-general `HH²`/`[α,α]` block stays live in every characteristic behind the caveat.
    (b) **The deformation FUNCTOR is computed on the Hochschild DGLA `C(A)` via CS-over-ℚ** (ℓ₂+δ);
    this is **not** "pretending dg-Lie" — `C(A)` *is* a DGLA (Gerstenhaber). The `B(A)[1]` L∞ is the
    **small-model companion**; **`ℓ₄` can be nonzero** and is **never claimed zero** ("ℓ_n=0 for n≥5"
    is a sufficient collapse condition, not automatic).
    (c) **`ℓ₃` scope** — the char-0 Bardzell `ℓ₃` (monomial) is **feasibility-gated**; v1 always
    ships the **rad²=0 dg-Lie certificate** (unconditional) + the induced-ℓ₂ ≡ CS-bracket **model-
    independence theorem statement** (§4, cited), and — **only if the spike lands** — the *computed*
    induced-ℓ₂ crossengine self-cert (needs the `B(A)` ℓ₂ adapter) plus the validated truncated `ℓ₃`;
    if the spike does not land, **no `ℓ_{≥3}` bracket and no computed `B(A)` ℓ₂ ship** — a **ledgered
    deferral** — recorded in the Change log. Neither is a fabricated general `ℓ₃` (M1/W1 rulings).
    (d) **MC = 2-cocycles is gated per instance** by the MRRS Thm-5.4 quiver hypotheses (`# PIN`);
    outside them, the honest order-by-order DGLA MC with the `[α,α]` obstruction, to a **certified
    truncation order** (bounded by the `order` **parameter**, not a bracket window — the native CS
    bracket runs at any degree; never a claimed complete formal solution).
    (e) **The presented `A_α` feedback covers the RADICAL (admissible) sub-locus.** Unit-direction
    deformations (e.g. `k[x]/(x²) ⤳ k[x]/(x²−t)`, the semisimplification) leave the admissible
    `kQ/I` class and are an `HH²`/MC story on `C(A)` — a **loud refusal** in `deformed_algebra`
    (the length-0 term is caught at parse level by `RelationError`, *before* the rad² admissibility
    gate; `AdmissibilityError`/`NotFiniteDimensionalError` cover the parseable non-admissible /
    infinite radical cases).
    RRRV's quiver-and-relations presentation is over **algebraically closed `k`** (honest base-change
    wording, P62/P70 precedent — no fabricated "algebraically closed" verdict).
    (f) **The bracket/obstruction constants are basis-dependent** (`basis="cs/<domain>"`);
    cross-engine comparison uses only basis-independent data (dims, verdicts, obstruction
    rank/class-nonvanishing) — the Plan-35 rule.
    (g) **A REAL obstruction exists and is pinned — on basis-independent facts only** (`QuantumCI(0)`:
    `unobstructed is False`, the obstruction class is **nonzero in `HH³`**). The `obstruction_map`
    scans the **whole quadratic form** `Q(c)=Σ c_i c_j [α_i,α_j]` because unobstructedness needs
    `[α,α]=0` for *all* `α`. **Whether the diagonal vanishes is basis-DEPENDENT** — live-verified this
    fix round: `[α₁,α₁] ≠ 0` on this engine (a basis cocycle self-obstructs), whereas the author's
    earlier basis had a vanishing diagonal with a combination witness — so it is **never pinned**
    (Minor-2 correction). The literal witness string and the raw constants are likewise
    basis/representative-dependent and are **never** hard-pinned (Plan-35 basis rule).
    Gentle nilpotent examples are unobstructed (Thm 5.4).
    (h) **Cost/budget (re-benchmarked, H1/H2)** — the CS obstruction bracket over ℚ is the driver
    (**18.1 s at dim 4, top 3** on this machine; HH dims ≈ 0.1 s), and its cost tracks **HH²/HH³
    richness × CS-resolution size, NOT `A.dim`** (dim-20 `kZ₁₀/J²` HH²=0 → 0.024 s vs dim-4
    `QuantumCI(−1)` HH²=5 → 43.5 s). The obstruction uses a **targeted (2,2)** build (early-out; at
    top 3 ≈ the full `(p,q)` table since the degree-2 block dominates). The plan's **OWN cap
    `DEFORM_MAXDIM = 32`** (NOT P70's 48 — different cost law) is a coarse backstop; the estimator's
    real routing pre-computes the cheap `HH²`/`HH³` dims and sizes off `(dim HH²)² × resolution size`,
    routing off the instant tier by HH-richness (not dim); the dim-220 Nakayama webapp examples carry
    **no** `deformations` (Plan-35 omission precedent).
    (i) **QPA/GAP have no Hochschild deformation / L∞ / MC surface** (live `NamesGVars()` sweep; an
    honest skip that FAILS if it appears); the genuine QPA cross-check is the **input** `HH²/HH³`
    and `E(A_α)` (Plan-27 bridge), not the deformation theory itself.
  - **Recount the class table** (`tests/release/test_oracle_classes.py` drives the numbers — run
    collection, paste LIVE counts, re-run to green; no guessed at-authoring number).
- [ ] **Step 2: README.** One features line (rep-theory-first per the standing positioning): "Formal
  deformations of `kQ/I` in characteristic 0 — the infinitesimal deformations `HH²`, the
  obstruction `[α,α] ∈ HH³`, Maurer–Cartan elements (= 2-cocycles in the nilpotent gentle regime,
  MRRS), the presented deformed algebra `A_α` fed back with its Ext-algebra, and the `B(A)[1]` L∞
  companion (rad²=0 dg-Lie certificate; honest L∞ with possibly-nonzero ℓ₄) — R13."
- [ ] **Step 3: Full gate:** `... tests/hochschild -q -m deep`, `... tests/webapp tests/gui
  tests/hpc -q -m fast`, `... tests/qpa -q -m qpa`, `... tests/release tests/citations -q`, plus a
  `QUIVERLAB_NO_NUMBA=1` parity spot-run of `tests/hochschild/test_deformations.py` — all green.
- [ ] **Step 4:** commit `docs(verification): P78 deformation/L∞ oracle rows + honest scope (char-0 gate, DGLA-on-C(A) vs L∞-on-B(A), l3 feasibility-gated, radical A_alpha sub-locus, real obstruction pin) + recounted classes; README line; tick P78`.

---

## Acceptance (Plan-78 definition of done)

1. **`quiverlab.hochschild.deformations` public:** `infinitesimal_deformations`,
   `obstruction_map`, `maurer_cartan`, `is_l_infinity_nilpotent`, `dg_lie_certificate`,
   `l3_bracket`, `deformed_algebra`, `deformation_structure`, `deformations_block`;
   `Algebra.deformation_structure`/`obstruction_map`/`deformed_algebra` delegates. The field-general
   `HH²`/`[α,α]` compute over any exact field; the deformation/L∞ interpretation is char-0 gated.
2. **Infinitesimal + obstruction pins green:** `HH²(k[x]/x²)/ℚ = 1`; `dim HH²/HH³ ==
   hochschild_cohomology` on the zoo; the **real obstruction** `QuantumCI(0)/ℚ` pinned on
   **basis-independent facts only** — `unobstructed is False` and a **nonzero obstruction class in
   `HH³`** (via the full quadratic-form scan). **The diagonal-vanishing is basis-DEPENDENT and is
   NOT pinned** (live: `[α₁,α₁] ≠ 0` on this engine; the witness string / raw constants are never
   hard-pinned — Minor-2); gentle nilpotent examples unobstructed (Thm 5.4); `[α,α]` a cocycle; the
   char-p field-general fallback + caveat.
3. **Maurer–Cartan green + gated:** the MRRS nilpotent-regime gate (char-0, gentle, `# PIN`
   hypotheses) ⇒ `MC = Z²`; outside it, the order-by-order DGLA MC on `C(A)` with the `[α,α]`
   obstruction gate and a certified truncation order (never a claimed complete formal solution);
   char-p loud refusal on the interpretation.
4. **Presented `A_α` feedback green:** `deformed_algebra` builds the radical `A_α` (flat
   `dim A_α = dim A`), its HH ≡ the direct family member (QuantumCI two-ways), and its Ext-algebra
   (Plan 27) is the RRR handoff; the unit direction **refuses loudly at parse level** via
   `RelationError` (the length-0 term), the parseable non-admissible / infinite cases via
   `AdmissibilityError`/`NotFiniteDimensionalError` (the radical-collapse honest scope); RRRV
   base-change wording honest.
5. **`B(A)[1]` L∞ companion green:** the **rad²=0 dg-Lie certificate** (`kZ_n/J²` True,
   `QuantumCI(0)` False) ships **unconditionally**; the **induced-ℓ₂ ≡ CS-bracket** agreement is the
   §4 **model-independence theorem statement** unconditionally, promoted to a *computed* crossengine
   self-cert **only if the spike builds the `B(A)` ℓ₂ adapter** (M1); and either the validated char-0
   Bardzell `ℓ₃` (monomial, RRB truncated pin — spike landed) or the **ledgered deferral** with **no
   `ℓ_{≥3}` shipped** (spike did not) — the Change log records which, and `ℓ₄` is never claimed zero.
6. **One algebra-only kind `deformations`** clickable end-to-end (canvas → block → report) in **all
   four locales** with exact key parity, schema v1, both runners byte-identical, grammar in **all
   three sites**, ONE golden with a documented change-log entry, canonical keys byte-stable, the
   estimator sizing off the bracket cost; the deformation/obstruction/MC/dg-Lie + display-only `A_α`
   render; the **adopt-flow ledgered** for P80.
7. `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the nine honest-scope
   entries (a)–(i); citations added and BibTeX-verified (three journal refs confirmed, two `@misc`
   PINs, RRR = three authors); README line; metaplan P78 ticked; deep + fast + qpa + release +
   citations green; a `QUIVERLAB_NO_NUMBA=1` parity spot-run green. **Prerequisite check (Minor-5
   ruling — "committed" means DOCS-only):** as of authoring **only the P70 and P71 *docs* are
   committed — their CODE is NOT on `dev`** (verified: `hochschild/lie_module.py` absent, no
   `hh1_lie` grammar branch in `hpc/spec.py`/`docs/gui/runner.py`/`webapp/server/schema.py`). So
   "needs P70/P71" is an **implemented-first** dependency: **P70/P71 must be implemented and merged
   to `dev`** before **Task 5 (the GUI kind, which copies P70's `hh1_lie` pattern)** runs. P78's own
   *compute* (Tasks 1–4) consumes only **P51's native bracket** and **P27's `ext_algebra`**, both
   already on `dev`, so it can be authored/implemented in parallel and merges after P70/P71 land.

---

## Methodology & assumptions

**Approach.** I (1) fast-forwarded my isolated worktree from the stale `b7fa566` (v0.3.0) to the
`dev` tip `1fb1cda` (P71's **doc** committed) after confirming `b7fa566` is an ancestor of `dev`;
the venv resolves `import quiverlab` to the **main-repo `dev` src** (`.../quiverlab/src`, version
`0.3.0` — bump waits for P80), so every benchmark and pin below ran on the **engine P78 actually
builds on** (P51 native CS bracket, P27 ext_algebra, the CS-over-any-Domain stack). **Note (Minor-5):
"P70/P71 committed" means their DOCS are committed — their CODE is NOT on `dev`** (`lie_module.py`
absent, no `hh1_lie` GUI kind); the fix round re-verified this in the worktree src. (2) Read the
metaplan P78 card + §§1–4 conventions and the R13 record verbatim; read the **committed P70 and P71
docs** (docs only — see the Minor-5 note) in full and honored their binding conventions (char-0 HARD
gate as a predicate; arithmetic-field
base-change wording, P62/P70 precedent; the entry-wise self-cert precedent — here the induced-ℓ₂ ≡
CS-bracket anchor and the `[α,α]`-is-a-cocycle self-cert; non-normative basis-dependent constants,
`basis="cs/<domain>"`, Plan-35 rule). (3) Verified all five citations by fetching arXiv abstracts
verbatim and resolved the two flagged items by search: **RRR 2202.01199 = three authors CONFIRMED**
(no Verdecchia); **RRRV 2003.10366 = four authors** (adds Verdecchia); the two journal refs the card
asserts are **both CONFIRMED** (JPAA 226(5) 2022 art. 106935; Comm. Algebra 53(9) 2025 3984–4007,
DOI confirmed). (4) Ran the **MANDATORY spec-time benchmark** live (CS/ℚ vs CS/GF(p) vs Bardzell/
GF(p) at the MC-relevant degree, the high-degree divergence, and the obstruction-bracket cost) and
took the route decision **CS-over-ℚ backbone; feasibility-gated char-0 Bardzell ℓ₃** with the numbers
recorded. (5) Ran a **survey subagent** over the Bardzell/gentle/build/Gröbner surfaces (the
field-free-combinatorics finding is load-bearing for the ℓ₃ route). (6) **Live-verified** every pin
(the `k[x]/x²` tangent, a **real obstruction** on `QuantumCI(0)`, the gentle-unobstructed record, the
flat `QuantumCI(q)` `A_α` feedback with moving HH, the `ext_algebra` handoff).

**Assumptions (stated, not silently relied on).** (i) The engine at `dev` tip is what P78 builds on
(true — venv resolves there; my worktree only holds the doc). (ii) RRB's ℓ₃ transfer formula is
implementable over a char-0 Bardzell complex — I verified the *substrate* is field-free and reusable
and the *structural* claims (rad²=0 dg-Lie; MC in degree 2) from the abstracts, but **NOT** the ℓ₃
formula itself; hence Task 4 is a **feasibility spike with a freeze**, not an unconditional ship. (iii)
The MRRS Thm-5.4 hypotheses are the R13 "gentle, no parallel arrows, no oriented cycles" — I have the
abstract, not §5, so the exact gate is `# PIN`. (iv) `C(A)`-DGLA and `B(A)`-L∞ compute the same
deformation functor (a theorem — L∞-quasi-iso invariance — so the CS-over-ℚ backbone is honest, not a
shortcut).

**What I deliberately did NOT check, and why it's safe.** (a) I did **not** implement RRB's ℓ₃
formula or fetch the full papers' §-level formulas — that is implementation work behind the Task-4
feasibility gate; the plan is honest that the truncated ℓ₃ value is `# PIN` and the general ℓ₃ may be
ledgered. (b) I did **not** run GAP/QPA (avoiding interference with concurrent `-m qpa` runs and
because QPA has no deformation surface — the honest-skip design is verified structurally, not by
running it). (c) I did **not** exhaustively search for a *gentle* real obstruction — the honest
record (gentle nilpotent ⇒ unobstructed, Thm 5.4) plus the *monomial* real obstruction (QuantumCI(0),
in RRB's scope) together satisfy the card's "construct one or honestly record unobstructed" without
needing a gentle counterexample (which the theorem says should not exist in the nilpotent regime).
(d) I did **not** micro-optimize the obstruction bracket (targeted (2,2) vs full table) beyond noting
the cost — the implementation note (targeted (2,2)) is a correctness-preserving speedup the worker
applies. (e) I did **not** re-run the whole suite — no `src/` changed; the plan writes only this doc.

---

## Open design risks (what a critic will likely attack)

1. **"You promised ℓ_n on B(A) and shipped mostly CS ℓ₂."** Mitigation: §4 (the DGLA-equivalence)
   shows the deformation *functor* needs only ℓ₂+δ on `C(A)` (a genuine DGLA — not "pretending");
   the `B(A)` ℓ_{≥3} is a small-model presentation whose exact v1 instances are the rad²=0 dg-Lie
   certificate and the feasibility-gated truncated ℓ₃; ℓ₄ is honestly never claimed zero. The card's
   "honest L∞, do NOT pretend dg-Lie" is honored (we only claim dg-Lie where RRB prove it: rad²=0).
2. **"Is scanning the diagonal enough for the obstruction?"** Mitigation: **no**, and the reason is
   basis-independent — unobstructedness needs `[α,α]=0` for *all* `α ∈ HH²`, so `obstruction_map`
   must scan the whole quadratic form `Q(c)=Σ c_i c_j [α_i,α_j]`. Whether a *basis* cocycle already
   self-obstructs is **basis-DEPENDENT** (live this fix round: `[α₁,α₁] ≠ 0` on this engine — a basis
   witness; the author's earlier basis had a vanishing diagonal and a *combination* witness `α₀+α₂`),
   so the machinery can rely on neither the diagonal alone nor a combination. An explicit test asserts
   (basis-independently) `unobstructed is False` and a nonzero obstruction class in `HH³`; it does
   **not** pin the witness string or the diagonal pattern (Minor-2).
3. **"The classic `k[x]/(x²−t)` example isn't a presented `A_α`."** Mitigation: it is pinned as an
   `HH²`/MC story (Pin 1) and `deformed_algebra` **refuses it loudly** — at *parse level* via
   `RelationError` on the unit's length-0 term (live-verified), never a silent build — with the
   radical-collapse honest scope; the presented feedback is validated on the **radical**
   QuantumCI(q) family instead. This is a feature (honest scope), not a gap.
4. **"Char-0 gate vs field-general HH."** Mitigation: the split is explicit — raw `HH²`/`[α,α]` any
   field (behind a caveat), deformation/L∞ interpretation char-0 (RBB/MRRS scope, the `1/n!`/`½`).
   The P70 precedent (char-0 HARD gate on the classification, field-general on the rest) is followed.
5. **"Feasibility-gated ℓ₃ is a hedge."** Mitigation: it is a **scope-freeze-after-spike**, the P79
   precedent the metaplan sanctions; both outcomes (ship truncated ℓ₃ / ledger the general ℓ₃ with
   **no `ℓ_{≥3}` shipped**) are complete honest v1 deliverables. The **rad²=0 dg-Lie certificate**
   ships unconditionally; the induced-ℓ₂ ≡ CS-bracket **model-independence theorem statement** is
   unconditional too, but its *computed* crossengine self-cert is gated behind the spike (M1 — no
   `B(A)` ℓ₂ adapter, nothing to compute, so it is never a vacuous pass).
6. **"Cost — and the earlier 22 s figure / dim-based sizing were wrong."** Mitigation (H1/H2
   re-benchmark): the ℚ bracket is **18.1 s at dim 4, top 3** on this machine (the "22 s" is
   superseded; a prior reviewer measured 5.27 s — machine-dependent, this machine is the reference),
   and the cost tracks **HH²/HH³ richness × resolution size, NOT `A.dim`** (dim-20 HH²=0 → 0.024 s vs
   dim-4 HH²=5 → 43.5 s). So the obstruction uses a targeted (2,2) early-out build, the plan's **OWN
   `DEFORM_MAXDIM = 32`** (not P70's 48) is a coarse backstop, the estimator pre-probes the cheap
   `HH²`/`HH³` and routes off HH-richness (not dim), and the dim-220 examples carry no `deformations`
   (Plan-35 precedent). The bracket-cost benchmark table is in the Route-decision section.
7. **"GUI adopt-flow deferral is silent."** Mitigation: it is an **explicit** GUI-deferral ledger
   entry (metaplan §1.2) that P80 reconciles with Marco's sign-off; v1 ships `A_α` display-only.
8. **"P70/P71 dependency — and 'committed' overstates it."** Mitigation: the metaplan marks P78
   *needs P70/P71*, and the fix round pins down that **only their DOCS are committed — their CODE is
   NOT on `dev`** (verified in the worktree src: no `lie_module.py`, no `hh1_lie` GUI kind). So the
   dependency is **implemented-first**: **Task 5 (the GUI kind, which copies P70's `hh1_lie` pattern)
   is sequenced after P70/P71 are implemented and merged to `dev`** (Acceptance §7). P78's own compute
   (Tasks 1–4) consumes only P51 (native bracket) and P27 (`ext_algebra`), both already on `dev` — so
   Tasks 1–4 can be authored/implemented in parallel, and the whole plan merges after P70/P71 land.

---

## Change log

- **2026-08-08 — plan authored** (P78 / R13). Worktree fast-forwarded to `dev` tip `1fb1cda`;
  benchmark + all pins run on the dev-tip engine (venv `0.3.0`). Route decided **CS-over-ℚ backbone
  + feasibility-gated char-0 Bardzell ℓ₃** with live numbers. Five citations verified (RRR = three
  authors; two journal refs confirmed). Real obstruction pinned (`QuantumCI(0)`). Task-4 spike
  outcome (truncated ℓ₃ shipped vs ledgered) **TO BE RECORDED HERE at implementation.**
- **2026-08-08 — critic NEEDS WORK → fix round applied** (adjudicated all-valid; every touched
  number RE-RUN LIVE this round on the worktree engine, single run wall-clock ±20%):
  - **H1 (MAJOR) — bracket cost re-benchmarked; driver is HH-richness, not dim.** New live table
    (Route-decision section) over `(dim, dim HH², dim HH³)` on 5 points incl. the **dim-8 HH-rich**
    `QuantumCI(0,a=2,b=4)` (`HH•=[2,4,9,13]`, HH²=9, 45 pairs) whose obstruction **did not finish in
    a 10-min synchronous run** — that timeout IS the datum. Confirmed dimension is not the predictor
    (dim-20 `kZ₁₀/J²` HH²=0 → 0.024 s vs dim-4 `QuantumCI(−1)` HH²=5 → 43.5 s). The `deformations`
    kind now gets its **OWN cap `DEFORM_MAXDIM = 32`** (NOT P70's `DEFAULT_MAXDIM = 48` — different
    cost law) + an estimator that pre-probes cheap `HH²`/`HH³` and sizes off `(dim HH²)² × resolution
    size`, routing off instant by HH-richness (scope-gate table, Task 5 Step 2, honest-scope (h),
    Open risk 6).
  - **H2 (MAJOR) — the "≈ 22 s at dim 4" figure fixed everywhere** to this machine's measured
    **18.1 s** (ℚ) / **2.37 s** (GF 32003); targeted (2,2) **17.9 s** ≈ full table (the "targeted is
    markedly cheaper" claim was false at top 3 — the degree-2 block dominates). Superseded at all 6
    sites (Route decision, Reuse list, scope-gate, Task 5, honest-scope (h), Open risk 6); a prior
    reviewer's 5.27 s noted as machine-dependent (this machine is the reference).
  - **M1 — the induced-ℓ₂ ≡ CS-bracket oracle** is now the §4 **model-independence THEOREM statement**
    unconditionally (cited), promoted to a *computed* crossengine self-cert **only if the Task-4
    spike builds the `B(A)` ℓ₂ adapter** (without it there is nothing to compute against — never a
    vacuous pass). New Task-4 Step 2b′; reconciled in the DECISION, marker lists, honest-scope (c),
    Acceptance §5, Open risk 5.
  - **W1 — L∞ deliverable re-captioned honestly:** worst case ships **NO `ℓ_{≥3}`**; the GUARANTEED
    content is the **DGLA on `C(A)` + the rad²=0 collapse certificate**; the higher-`ℓ_n` route
    (RRB's homotopy-transfer / contracting-homotopy on Bardzell's complex, `rrb_linfty_bardzell`
    2008.08122) is cited (Goal bullet).
  - **Minor — refusal class:** the `x²−t` unit direction is refused at **parse level by
    `RelationError`** ("term '1' has no arrows"), NOT the rad² `AdmissibilityError` gate
    (live-verified); fixed at Pin 1, scope-gate, Task 3, Reuse list, honest-scope (e), Acceptance §4,
    Open risk 3.
  - **Minor — witness basis-dependent + a LIVE CORRECTION:** the witness string / constants are never
    hard-pinned (pin `unobstructed is False`, nonzero class, degree/dim). **Re-running the (2,2) form
    disproved the earlier "diagonal all vanish" claim** — on this engine `[α₁,α₁] = 2 ≠ 0` (a basis
    cocycle self-obstructs); the diagonal pattern is basis-DEPENDENT and is now removed from every pin
    (Pin 2, Task 1, honest-scope (g), Acceptance §2, Open risk 2, oracle_selfcert). The scan-the-whole-
    quadratic-form guidance stands on the basis-independent reason (`[α,α]=0` for *all* `α`).
  - **Minor — `q = −1` not `q = 1`** for the `[4,4,5,6]` HH pin (re-verified: quiverlab's `xy+q·yx`
    convention puts the commutative `k[x,y]/(x²,y²)` at `q = −1`; `q = +1` → `[2,4,6,8]`); Pin 4 prose.
  - **Minor — MC order limit is the `order` PARAMETER**, not "the bracket window" (the native CS
    bracket runs at any degree past the bar window); fixed at the scope-gate, Task 2, honest-scope (d),
    field comment.
  - **Minor — "P70/P71 committed" = DOCS-only:** verified in the worktree src (no `lie_module.py`, no
    `hh1_lie` GUI kind) — their CODE is NOT on `dev`; the GUI **Task 5 is phrased implemented-first**
    (sequenced after P70/P71 land), and the module-placement note, methodology, Acceptance §7, and
    Open risk 8 say so. P78's own compute (Tasks 1–4) consumes only P51 + P27 (both on `dev`).
