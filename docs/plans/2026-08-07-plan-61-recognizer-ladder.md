# Plan 61: the quasi-tilted / shod / weakly-shod / laura / ada recognizer ladder (P61 / R18) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A single no-code button that classifies a representation-finite algebra `A`
against the **Assem-school recognizer ladder** — the five nested per-instance
certificates **quasi-tilted ⊂ shod ⊂ weakly shod ⊂ laura** and **ada** — and, for the
**ada** class over an algebraically closed field, turns quiverlab's shipped `HH¹` into a
**complete simple-connectedness oracle** (ACLV Theorem B: *A ada over algebraically
closed k ⟹ (A simply connected ⟺ HH¹(A) = 0)*). Every verdict is a definite
`True`/`False` with a **witness on `False`** and a **certificate on `True`**; the five
verdicts are computed off **ONE knit and ONE P55 `LeftRightAtlas` call** plus the shipped
`global_dimension`, and the **nesting theorems are asserted as a standing monotone
self-cert** (a `True` below a `False` above is a loud internal error). One no-code GUI
compute kind — `recognizer_ladder` — puts the whole ladder, the finite laura complement,
and the ada/HH¹ simple-connectedness block one click away. Exact only; loud typed refusals
at the rep-finite / algebraically-closed scope boundaries; no floats.

**Architecture:** ONE new module, a thin exact classifier layer **over the P55 substrate**
(no new AR machinery, no new homology engine — everything is read from `LeftRightAtlas` +
`global_dimension` + `hochschild_cohomology`):

- **`src/quiverlab/modules/recognizers_ladder.py`** (new) — the whole surface. It consumes
  the P55 `left_right_parts(A)` atlas (Plan 55, `modules/left_right.py`) and classifies; it
  is a `modules/` citizen exactly like `modules/left_right.py` (P55) and its battery lives
  in `tests/modules/` (the **deep** bucket, `tests/conftest.py`). Public surface:
  - `recognizer_ladder(A, *, budget=256) -> RecognizerLadder` — the primary compute: ONE
    `left_right_parts(A, budget=budget)` call, ONE `A.global_dimension()`, the five
    verdicts, the nesting + theorem-gate self-certs, and the ada/HH¹ block, in one pass.
  - `@dataclass(frozen=True) RecognizerLadder` — mirrors the P55 `LeftRightAtlas` /
    `GlobalDimension` honesty contract (`is_complete`/`status`, a certified value or a
    labelled refusal, never a fabricated verdict). Carries the five `LadderRung`s, the
    laura `complement`, the `SimpleConnectedness` block, and JSON-ready summaries.
  - `@dataclass(frozen=True) LadderRung` — `{name, verdict: bool, certificate: dict,
    witness: dict | None}` — one per class; `witness` is populated only on `verdict=False`.
  - `@dataclass(frozen=True) SimpleConnectedness` — the ACLV-Theorem-B block
    `{applicable, field_algebraically_closed, hh1_dim, verdict: bool | None, theorem, note}`.
  - `recognizer_ladder_block(A, *, budget=256) -> dict` — the JSON block the two runners
    share (strips `Module`s to names + dim-vectors).
  - Thin `Algebra` delegates in `core/algebra.py` beside `left_right_parts`
    (P55) / `ar_quiver` (P41): `Algebra.recognizer_ladder(budget=256)`,
    `Algebra.is_quasi_tilted(budget=256)`, `Algebra.is_shod(budget=256)`,
    `Algebra.is_weakly_shod(budget=256)`, `Algebra.is_laura(budget=256)`,
    `Algebra.is_ada(budget=256)` (each returns the corresponding `bool` off the same
    cached ladder; lazy-import to avoid the `modules → core` cycle).

No new math engines; no new AR machinery. Every rung is a bounded read off the finite
`LeftRightAtlas` (parts, complement, `projective_placement`/`injective_placement`) plus the
shipped `Algebra.global_dimension()` and — for the ada/HH¹ headline — the shipped
`Algebra.hochschild_cohomology(1)`. The one genuinely new sweep is the **weakly-shod
bounded-path test** (a finite SCC computation on the P41 AR-quiver arrows, §"weakly shod").

**Tech Stack:** exact `Domain` arithmetic throughout — P55 `left_right_parts` /
`LeftRightAtlas` (rep-finite scope, honest `status`), `Algebra.global_dimension() ->
GlobalDimension` (`.value`, `.exact`), `Algebra.hochschild_cohomology(1)` (any exact field;
`HHTable.dims[1]`), the P41 `ARQuiver.arrows` (irreducible-map digraph, trustworthy only
when `is_complete`), and the new **additive** `Domain.is_algebraically_closed` flag
(stamped `True` only on the CC working domain) for the Theorem-B gate — read off
`A.domain`, quiverlab's coefficient-domain slot (there is **no** `Algebra.field`
attribute). No floats in `src/` (AST-gated by
`tests/test_no_floats.py`); all dim-vectors `dict[vertex,int]`, all dims `int`, all verdicts
Boolean, all witnesses index/name sets.

---

## Record (R18, verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R18 — Quasi-tilted / shod / weakly-shod / laura / ada recognizer ladder.**
> [C-scout P3+P4+P11; keep, sequenced after R15] Object: the nested per-instance
> certificates — quasi-tilted (gl.dim ≤ 2 and every indec pd≤1 or id≤1), shod, weakly shod
> (bounded I⇝P paths), laura (L_A ∪ R_A cofinite; report the finite complement), ada (all
> projectives+ injectives in L_A ∪ R_A) — and for ada the THEOREM 'simply connected ⇔ HH¹ =
> 0' (verified in arXiv:1102.1188, algebraically closed k) turning quiverlab's HH¹ into a
> complete simple-connectedness oracle on that class. Rep-finite scope. Refs: Happel–Reiten–
> Smalø Mem. AMS 575; Coelho–Lanzilotta J. Algebra 265 (2003); Assem–Coelho J. Algebra 269
> (2003); Smith arXiv:math/0702562; Bordino–Fernández–Trepode arXiv:1404.5294 (Comm. Alg.
> 2017). Size M total.

Tier β, C-cluster spine; dependency chain (research doc §9): **R15 → R16 → {R17, R18, R21,
R37}**. P61 is the second consumer of the P55 substrate (metaplan §5: **P61 needs P55**),
sibling to P60 (tilted recognizer, R17). Branch `plan-61-recognizer-ladder` off `dev`
**after P55 merges** (its `modules/left_right.py::left_right_parts` is a hard prerequisite).

---

## Reference re-verification (mandatory standing rule; done at spec time)

Re-verified via web + the arXiv PDF of **1102.1188** (read pages 1–6, verbatim) and the
**"Organising the module category" survey PDF** (São Paulo J. Math. Sci. 16 (2022) 62–82,
§1 and §4, verbatim). **Findings — the ladder's ground truth (transcribe the `# PIN`s at
citation-add time):**

1. **Happel–Reiten–Smalø, Mem. Amer. Math. Soc. 120 (1996), no. 575** —
   *Tilting in Abelian Categories and Quasitilted Algebras* (ISBN 9780821804445). The
   **quasi-tilted characterisation** the record quotes is verbatim in the survey §4
   (attributing it to HRS): *"a quasitilted algebra A can be characterised by the
   properties: (QT1) gl.dim A ≤ 2; and (QT2) for each indecomposable module M, either
   pd_A M ≤ 1 or id_A M ≤ 1. They have also shown that (QT1) and (QT2) are independent and
   that (QT2) implies that gl.dim A ≤ 3."* `# PIN`: volume/number (120 / 575) confirmed via
   the AMS Memoirs catalogue at citation-add time.

2. **Coelho–Lanzilotta, J. Algebra 265 (2003), no. 2, 379–403** — *Weakly shod algebras*.
   The record's shod + weakly-shod reference. `shod` originates in Coelho–Lanzilotta,
   *Algebras with small homological dimensions*, Manuscripta Math. 100 (1999) 1–11 (the name
   = **s**mall **ho**mological **d**imension); the J. Algebra 265 (2003) paper recalls shod
   and introduces weakly shod. `# PIN`: pages 379–403 confirmed via the Semantic Scholar /
   Elsevier record; the worker either adds the 1999 Manuscripta entry as the shod-origin
   secondary or leaves the recorded J. Algebra 265 (2003) entry (house rule: no fabricated
   fields). **Weakly-shod condition (WSA), survey §4 verbatim:** *"there exists a positive
   integer n₀ such that any path of irreducible morphisms in ind A from an injective to a
   projective has length bounded by n₀"* — and *"for a shod algebra, [WSA] is true"* (shod
   ⟹ weakly shod). Survey §4 also: *"Reiten-Skowronski also studied weakly shod algebras
   under the name double tilted algebras."*

3. **Shod ⟺ complement empty — survey §4, Theorem 4.1 [40], verbatim.** *"The following are
   equivalent for an algebra A: (a) A is shod, that is, for each indecomposable module M,
   either pd M ≤ 1 or id M ≤ 1. (b) ind A = 𝓛_A ∪ 𝓡_A. …(f) any path from an indecomposable
   injective module to an indecomposable projective module can be refined to a path of
   irreducible morphisms and any such refinement has at most two hooks, and, in case there
   are two, they are consecutive."* ⟵ **(a)⟺(b) is the decisive P61 characterisation:
   `A is shod ⟺ atlas.complement == []`.**

4. **Assem–Coelho, J. Algebra 269 (2003), no. 2, 456–479** — *Two-sided gluings of tilted
   algebras* (the **laura** paper; laura was introduced here and independently by
   Reiten–Skowroński). **laura definition:** `A` is *laura* ⟺ `ind A ∖ (𝓛_A ∪ 𝓡_A)` is
   **finite** (i.e. `𝓛_A ∪ 𝓡_A` is cofinite). `# PIN`: title/pages confirmed via the
   Elsevier record (ScienceDirect `S0021869303004368`).

5. **Assem–Castonguay–Lanzilotta–Vargas, arXiv:1102.1188** — *Algebras determined by their
   supports* (the **ada** paper). Read verbatim from the PDF:
   - **§1.1 standing hypotheses (verbatim):** *"Throughout this paper, all our algebras are
     basic and connected artin algebras."* ⟵ **connectedness and basicness are hypotheses
     of Theorem B**; the ladder inherits them (single connected component; the P41 knit
     already refuses disconnected input downstream, and basicness is quiverlab's default).
   - **Definition 2.1 (verbatim):** *"An artin algebra A is called an ada algebra if
     A ⊕ DA ∈ add(𝓛_A ∪ 𝓡_A). Clearly, this is equivalent to requiring that, for every
     x ∈ A₀, we have both P_x and I_x lying in 𝓛_A ∪ 𝓡_A."* Also: *"Quasi-tilted algebras
     are clearly ada. We call **strict** an ada algebra which is not quasi-tilted."*
   - **Quasi-tilted ⟺ projectives-in-left (Introduction, verbatim):** *"an artin algebra A
     is quasi-tilted if and only if every indecomposable projective module lies in the …
     left part of the module category, or equivalently if and only if every indecomposable
     injective module lies in the right part."* ⟵ the **second** quasi-tilted route,
     computed from `atlas.projective_placement`.
   - **Corollary 2.6 (verbatim):** ada ⟹ (a) *"for any indecomposable module M, we have
     pd M ≤ 2 or id M ≤ 1"*, (b) *"gl.dim A ≤ 4"*. **Remark 2.7(a):** the bound is **sharp**
     (Example 2.2(b) has gl.dim 4).
   - **Example 2.2(a) (verbatim):** *"Let A be a shod algebra. Then ind A = 𝓛_A ∪ 𝓡_A.
     Therefore A is ada."* ⟵ **shod ⟹ ada.**
   - **Example 2.2(b) (verbatim):** the linear rad²=0 Nakayama `1 ← 2 ← 3 ← 4 ← 5` is a
     *"(representation-finite) ada algebra"* — **the P55 shared fixture** (gl.dim 4).
   - **Example 2.2(c) (verbatim):** `1 ⇉ 2 ⇉ 3 ⇉ 4` bound by rad²=0 is a
     *"(representation-infinite) ada algebra … an ada algebra may have infinitely many
     indecomposables which are not in 𝓛_A ∪ 𝓡_A"* — ⟵ **ada does NOT imply laura in
     general**; the honest rep-infinite refusal oracle (mathematically ada, refused because
     rep-infinite).
   - **Theorem B (VERBATIM — the headline):** *"Let A be an ada algebra over an algebraically
     closed field. Then A is simply connected if and only if HH¹(A) = 0. Moreover, if this
     is the case, then the Hochschild cohomology ring HH•(A) reduces to the base field."*
   `# PIN`: 1102.1188 ships as the **verified arXiv `@misc`** entry (published venue not
   verified here). **This entry already exists in the tree as `aclv_supports` (added by
   P55)** — P61 reuses it, does not re-add it.

6. **Smith, arXiv:math/0702562** — *Almost laura algebras* (published in J. Algebra;
   ScienceDirect `S0021869307004449`). Context for the laura landscape (almost-laura
   generalises laura). The **rep-finite laura triviality** the ladder reports is a *direct*
   consequence of the Assem–Coelho definition — `complement ⊆ ind A` and `ind A` is finite
   in our scope ⟹ the complement is finite ⟹ **every representation-finite algebra is
   laura**. We cite Smith for the class landscape, not for the (elementary) triviality.
   `# PIN`: authors/venue confirmed at citation-add time; ships as the verified arXiv entry.

7. **Bordino–Fernández–Trepode, arXiv:1404.5294** — *On the quiver with relations of a
   quasitilted algebra and applications*, Comm. Algebra 45 (2017), no. 9, 4050–4061.
   Context for the quiver-with-relations structure of quasi-tilted algebras; a secondary
   reference on the quasi-tilted rung. `# PIN`: volume/issue/pages confirmed at citation-add
   time (arXiv `@misc` + upgrade-to-`@article` if BibTeX-verifiable).

8. **Survey "Organising the module category"** (Alvares–Assem–Castonguay–Vargas, São Paulo
   J. Math. Sci. 16 (2022) 62–82) — the single source that states all five definitions and
   the hierarchy together (§4). **This entry already exists as `organising_module_category`
   (added by P55)** — P61 reuses it. Also verbatim §5.1: *"if A is a strict weakly shod
   algebra (that is, a weakly shod algebra which is not quasitilted), then Hⁱ(A) = 0 for
   i ≥ 2. The result cannot be extended to arbitrary weakly shod algebras since there are
   quasitilted algebras with the second Hochschild cohomology group nonzero."* ⟵ a bonus
   `HH^{≥2}` cross-oracle on the strict-weakly-shod rung (Task C, optional pin).

---

## The five classes — exact definitions and the characterisations P61 computes

All in **representation-finite scope** (the P55 atlas is complete iff `A` is rep-finite and
non-self-injective; §"Global Constraints"). `𝓛_A`, `𝓡_A` are the P55 left/right parts;
`complement = atlas.complement = ind A ∖ (𝓛_A ∪ 𝓡_A)`.

| rung | definition (verbatim source) | **what P61 computes** |
|------|------------------------------|-----------------------|
| **quasi-tilted** | HRS: (QT1) gl.dim ≤ 2 **and** (QT2) every indec has pd ≤ 1 or id ≤ 1 | `gl.dim ≤ 2 ∧ complement == []` (QT1 ∧ shod); **cross-checked** against ACLV `every P_x ∈ 𝓛_A` (all `atlas.projective_placement` ∈ {"L","both"}) |
| **shod** | Coelho–Lanzilotta: every indec has pd ≤ 1 **or** id ≤ 1 (⟹ gl.dim ≤ 3) | `atlas.complement == []` (survey Thm 4.1 (a)⟺(b)) |
| **weakly shod** | Coelho–Lanzilotta (WSA): ∃ `n₀` bounding every irreducible-morphism path from an indec injective to an indec projective | `no non-trivial AR-quiver SCC lies on a directed path from an injective to a projective` (finite-universe boundedness; §"weakly shod") |
| **laura** | Assem–Coelho: `ind A ∖ (𝓛_A ∪ 𝓡_A)` finite | **always `True`** in our rep-finite scope; report the finite `complement` (size + members) as the useful datum |
| **ada** | ACLV Def 2.1: every `P_x` **and** every `I_x` in `𝓛_A ∪ 𝓡_A` | `all(atlas.projective_placement[v] ≠ "neither") ∧ all(atlas.injective_placement[v] ≠ "neither")` |

**The nesting (each a theorem; the standing monotone self-cert):**

- `tilted ⟹ quasi-tilted` (HRS: End of a tilting module over a hereditary algebra is
  quasi-tilted). *P61 does not compute `tilted` — that is P60 (R17); its `# PIN` cross-check
  below.*
- `quasi-tilted ⟹ shod` (QT1 ∧ QT2 ⟹ QT2; survey §4 — QT2 **is** the shod condition).
- `shod ⟹ weakly shod` (survey §4: WSA holds for shod).
- `weakly shod ⟹ laura` (standard inclusion; the weakly-shod complement is finite).
- `quasi-tilted ⟹ ada` and `shod ⟹ ada` (ACLV: "quasi-tilted algebras are clearly ada";
  Example 2.2(a): "shod ⟹ ada").
- `ada ⟹ laura` **in rep-finite scope ONLY** (both are `True` here — the complement is
  finite). **NOT a general theorem** (ACLV Example 2.2(c): rep-infinite ada, infinite
  complement, not laura). We assert `ada ⟹ laura` as a self-cert *inside the rep-finite
  atlas* and annotate the honest general-scope boundary.

**Monotone consistency (self-cert, loud on violation).** A `True` at a **more special**
rung with a `False` at a **more general** rung is an internal error:
`quasi_tilted ⟹ shod ⟹ weakly_shod ⟹ laura`, and `quasi_tilted ⟹ ada`, `shod ⟹ ada`,
`ada ⟹ laura` (rep-finite). `recognizer_ladder` raises `QuiverlabError` (never emits) if
any implication fails — it means the P55 atlas or the sweep is wrong, not the mathematics.

**Theorem gates (self-cert, loud on violation).** Two cheap **necessary** consequences pin
the verdicts against the shipped `gl.dim`: `shod ⟹ gl.dim ≤ 3` (CL) and
`quasi_tilted ⟹ gl.dim ≤ 2` (QT1). A `True` verdict whose `gl.dim` exceeds the bound is a
loud internal error. (These are gates on **True** verdicts; they do not run pre-knit — the
knit refusal is inherited from P55.)

---

## The P55 shared fixture — the full five-verdict ladder, derived and PINNED in-plan

`A = kQ/rad²` for the linear Nakayama `1 ← 2 ← 3 ← 4 ← 5` (arrows `α_i : i+1 → i`), the
**ACLV Example 2.2(b)** algebra and the **P55 worked example** (P55 plan §"Mathematical
foundation": 9 indecomposables, `𝓛_A = {S₁,S₂,P₂,P₃}`, `𝓡_A = {S₄,S₅,P₄,P₅}`,
`𝓛_A ∩ 𝓡_A = ∅`, **`complement = {S₃}`** (pd 2, id 2), `e_λ = {1,2,3}`, `e_ρ = {3,4,5}`,
gl.dim 4). **Recomputed live in the venv with the shipped engine** (knit + `hom_dim`
predecessor closure + `betti(2)` pd/id sweep + placements + the AR-SCC weakly-shod route):

| rung | verdict | reason (all recomputed live) |
|------|---------|------------------------------|
| **quasi-tilted** | **False** | gl.dim = 4 > 2 (QT1 fails); **and** `complement = {S₃} ≠ ∅`; **and** ACLV route: `projective_placement = {1:L, 2:L, 3:L, 4:R, 5:R}` — `P₄,P₅ ∉ 𝓛_A` (both routes agree on `False`) |
| **shod** | **False** | `complement = {S₃} ≠ ∅`; `S₃` has pd 2 **and** id 2 (neither ≤ 1) |
| **weakly shod** | **True** | the AR quiver is **directed** (no non-trivial SCC — no cyclic module), so every I⇝P irreducible path is trivially bounded |
| **laura** | **True** | rep-finite ⟹ finite complement; `complement = {S₃}`, size 1 |
| **ada** | **True** | `projective_placement = {1:L,2:L,3:L,4:R,5:R}`, `injective_placement = {1:L,2:L,3:R,4:R,5:R}` — every `P_x` and every `I_x` in `𝓛_A ∪ 𝓡_A` (matches ACLV Example 2.2(b) verbatim) |

This is a **full five-way discrimination on ONE fixture**: `(qt, shod, wshod, laura, ada) =
(F, F, T, T, T)`. The nesting holds monotonically (`weakly_shod=T ⟹ laura=T` ✓;
`ada=T ⟹ laura=T` ✓; no `True`-below/`False`-above). The **ada/HH¹** block: the quiver is
the linear tree `A₅` ⟹ `π₁(Q,I) = 1` for every presentation ⟹ **simply connected**; ACLV
Theorem B then forces `HH¹(A) = 0` over an algebraically closed field — and indeed
**`dim HH¹(A) = 0` recomputed live over both QQ and CC**. This fixture is the **primary
Theorem-B agreement oracle** (ada `True`, alg-closed field, HH¹ = 0, simply connected — all
three consistent, and agreeing with P56's `is_simply_connected(A).verdict is True` on the
same tree). Pinned verbatim in `tests/modules/test_recognizer_ladder.py`.

---

## The ada / HH¹ simple-connectedness oracle (ACLV Theorem B) — the headline wiring

Theorem B makes `HH¹` a **complete decision procedure** for simple connectedness **on the
ada class over an algebraically closed field** — precisely where P56's `is_simply_connected`
must return the honest `None` (triviality of a finitely presented group is undecidable,
Adian–Rabin). P61 **resolves** that `None` on the ada+alg-closed class.

**The gate and the three regimes** (`SimpleConnectedness` block, populated only when the
ladder is complete and `ada.verdict is True`):

1. **`ada.verdict is False`** → `applicable = False`, `verdict = None`, `note` = "Theorem B
   applies only to ada algebras". (HH¹ is still reported if cheap, but no SC verdict.)
2. **`ada.verdict is True` and the field is NOT algebraically closed** (QQ / GF(p) /
   GF(pⁿ)) → `applicable = False`, `field_algebraically_closed = False`, `hh1_dim =
   dim HH¹(A)` (computed — it is field-defined and cheap), `verdict = None`, `note` = the
   **Theorem B statement** + "the verdict requires an algebraically closed field; recompute
   over CC for the simple-connectedness conclusion." **We do NOT emit an SC verdict off an
   algebraically closed field** — Theorem B's hypothesis is not met and the biconditional
   can genuinely fail; reporting `HH¹` + the theorem as documentation is the honest surface.
3. **`ada.verdict is True` and the field IS algebraically closed** (`A.domain
   .is_algebraically_closed` — the new flag, `True` only on the CC working domain; §"the
   algebraically-closed gate" below) → `applicable = True`, `field_algebraically_closed =
   True`, `hh1_dim = dim HH¹(A)`, **`verdict = (hh1_dim == 0)`** — the complete Theorem-B
   decision: `True` = simply connected, `False` = not simply connected. `note` records the
   theorem + (on `verdict is True`) the corollary *"HH•(A) reduces to the base field"*.

**The algebraically-closed gate (respecified — BLOCKING review fix H1).** `Algebra` stores
its coefficient field as **`self.domain`** — **there is no `Algebra.field` attribute**
(verified live: `hasattr(A, "field") is False`, `hasattr(A, "domain") is True`). And
`field=CC` is realised at construction as a `quiverlab.fields.complexfield.SympyExactDomain`
— **`ComplexField.make_domain` discards the `ComplexField` sentinel** (verified live:
`A.domain` for `field=CC` is `SympyExactDomain`, `repr = "CC (computing exactly in QQ)"`).
Crucially, **`field=QQi` (Gaussian rationals ℚ(i), NOT algebraically closed) ALSO yields a
`SympyExactDomain`** (verified live; distinguished only by `sdom = QQ_I` and its `name`), so
CC and QQi **cannot be told apart by class** — an `isinstance(..., SympyExactDomain)` gate
would wrongly mark ℚ(i) algebraically closed. The clean route is a new **additive** boolean:
- **`Domain.is_algebraically_closed`** — a property/attribute on the base `Domain`, default
  **`False`** (so `RationalField`, `PrimeField`, `FiniteField`, and the **QQi**
  `SympyExactDomain` all report `False`).
- **`SympyExactDomain.__init__(self, sdom, *, algebraically_closed=False)`** stores the flag.
- **`ComplexField.make_domain`** returns `SympyExactDomain(sdom, algebraically_closed=True)`
  (CC declares ℂ); **`GaussianRationalField.make_domain`** (QQi) keeps the default `False`.
This is a small additive property task in `src/quiverlab/fields/` (base `Domain` + the CC
path); no `src/` behaviour changes elsewhere. **The gate is `A.domain.is_algebraically_closed`
— never the `CC` sentinel.**

**Gate rationale (honest).** The working domain is a **finite extension of ℚ computing
exactly** (CC's `sdom` is `QQ` for rational entries, an algebraic extension otherwise) — it
is not literally ℂ. But `field=CC` is the user's **declaration of ℂ-intent**, and `dim HH¹`
is invariant under the base extension to ℂ (**flat base change** — the boundary matrices are
defined over the smaller field, and rank is preserved), so `HH¹ = 0 ⟺ simply connected`
computed on the exact working domain **IS** the ℂ statement of Theorem B. `dim HH¹` is
therefore char-0 field-independent (over CC = over QQ); *the theorem's validity, not the
number, is what the algebraically-closed hypothesis controls* — hence the gate is a
**declaration flag**, not a re-computation. Over QQ/GF(p)/GF(pⁿ)/QQi the biconditional is not
claimed. **Test mandate (H1):** the gate test constructs an algebra with `field=CC` and
asserts **`A.domain.is_algebraically_closed is True`**, and with `field=QQ` / `field=GF(p)` /
`field=QQi` asserts `False` — never asserting on the `CC` object itself.

**The P56 agreement battery (the cross-oracle interplay — Task B).** Where both `A` is
ada+alg-closed (so P61 emits a definite `verdict`) **and** P56's `is_simply_connected(A)`
returns a definite `True`/`False`, the two verdicts **must agree** — a disagreement is a
loud test failure (an engine bug in one of HH¹, the parts, or π₁). Where P56 returns `None`
(Adian–Rabin), P61's Theorem-B verdict is the **complete** answer that resolves it — the
headline. The battery runs on the ada members of the zoo over CC. **P56 alignment point
(flagged):** P56 (Wave 2, independent) may merge before or after P61; the agreement test
imports `quiverlab.invariants.coverings.is_simply_connected` behind an **auto-flipping
skip** (`pytest.importorskip` → real assert the moment P56 lands — the Plan-29/31 xfail-
fence precedent). The **core** Theorem-B verdict (regime 3) does **not** depend on P56 — it
is self-contained (HH¹ + the alg-closed gate); only the agreement battery does.

---

## weakly shod — the bounded-path sweep (the one genuinely new computation)

The WSA definition (survey §4): `A` is weakly shod ⟺ ∃ `n₀` bounding the length of every
path of **irreducible morphisms** in `ind A` from an indecomposable injective to an
indecomposable projective. On the **finite** AR quiver (rep-finite), the supremum of
path-lengths from the injective-set to the projective-set is infinite **iff** some vertex
lying on an oriented cycle is both reachable-from-an-injective and can-reach-a-projective.
The AR quiver is loopless (a translation quiver, survey §8), so a vertex lies on an oriented
cycle **iff its strongly-connected component (SCC) is non-trivial** (size ≥ 2). Hence:

> **`A` is weakly shod ⟺ no vertex `Z` in a non-trivial AR-quiver SCC satisfies
> `(∃ injective I : I ⇝ Z)` and `(∃ projective P : Z ⇝ P)`.**

*Computation (finite, off the P41 knit).* Build the AR-quiver arrow digraph
`G = ar.arrows` on the `N` indecomposables (trustworthy because `ar.is_complete` — the same
gate P55 uses for its cross-engine tie). Reflexive-transitive closure `reach` (Warshall,
`N ≤ budget`). SCCs by mutual reachability: `Z` cyclic ⟺ `∃ Y ≠ Z : reach[Z][Y] ∧
reach[Y][Z]`. Injective/projective indices located by the P55 `_index_in_U` helper
(dim-vector-prefiltered `is_isomorphic`, **QQ / char-0 decisive** — inherits P55's M1 scope,
§"Global Constraints"). Then `weakly_shod = not any(cyclic Z with (some inj ⇝ Z) and
(Z ⇝ some proj))`; the witness on `False` is `{Z, inj_source, proj_target}`.

*Directed ⟹ weakly shod (the common case).* A representation-**directed** algebra has no
non-trivial SCC, so `weakly_shod = True` immediately — verified live on the fixture (rad²=0
A₅ is directed). The knit refuses self-injective (the main source of cycles), so most
in-scope inputs are directed and this rung is `True`; the SCC route is the **definition-
faithful** computation that still gets a non-directed rep-finite input right.

*The `False` branch is exercised on a REAL algebra (M2 — verified live).*
**`NakayamaAlgebra(kupisch=[3,2,2])`** (non-self-injective cyclic Nakayama, `dim 7`) knits
**complete** with an oriented AR cycle: recomputed live, its entire module category is ONE
non-trivial SCC (all 7 indecomposables cyclic), and the cycle **does** sit on an injective→
projective route (`I₂ ⇝ cycle ⇝ P_x`), so **`weakly_shod = False`** with a real witness — the
inj/proj identification + reach plumbing runs on a genuinely non-empty-SCC atlas (it is fully
degenerate for the other rungs — everything in the complement, so `ada = shod = quasi_tilted =
False`, `laura = True`, gl.dim 3 — which is exactly the point: the seam is exercised
maximally). Pinned as `is_weakly_shod(kupisch=[3,2,2]) is False`. In ADDITION, a
**hand-constructed AR-digraph unit test** drives the pure `_weakly_shod_from_digraph(n, arrows,
inj, proj)` helper (a synthetic `I → 2-cycle → P` returns `False` + the `{cycle}` witness;
removing the `I→cycle` edge returns `True`) — the guaranteed presentation-free `False`-branch
coverage.

*Equivalence to the Hom-closure route (a REAL cross-engine oracle — M1(b)).* P55 pins
`Hom-closure predecessor ≡ AR-reachability` for rep-finite (rad^∞ = 0), so `reach` above
equals P55's `atlas._leq`. We implement **both** weakly-shod routes and assert they agree
(cheap, both exist): the **primary** `_weakly_shod(A, atlas, budget)` reads the AR-quiver
arrows (the WSA "irreducible morphisms" are literally the arrows) for the SCC + reach; the
**cross-check** `_weakly_shod_hom(atlas)` reads `atlas._leq` for reachability and
`end_dim(Z) ≥ 2` for a module's self-cyclicity (a non-brick indecomposable — the Hom-closure
analogue of a non-trivial SCC). `test_weakly_shod_two_routes_agree` asserts
`_weakly_shod(A, ...).verdict == _weakly_shod_hom(atlas).verdict` on the fixture (both `True`)
AND on `kupisch=[3,2,2]` (both `False`) — a non-vacuous tie that also re-validates P55's
arrow completeness. (If P55 does not expose `end_dim` conveniently, `_weakly_shod_hom` imports
`modules.hom.end_dim` directly — it is shipped.)

---

## Consumers / what P61 reads from the P55 atlas

P61 is a **pure consumer** of `LeftRightAtlas` (P55 plan §"Consumers" names exactly these):

- `atlas.complement` (`ind A ∖ (𝓛_A ∪ 𝓡_A)`) → **shod** (`== []`, survey Thm 4.1 a⟺b),
  **laura** (finite, reported), **quasi-tilted** (route 1: `[] ∧ gl.dim ≤ 2`).
- **`atlas.pd_le_1` / `atlas.id_le_1`** (per-indecomposable booleans, index-aligned with
  `_modules`; **P55 CONTRACT AMENDMENT — see below**) → the **shod QT2 route** (every indec
  has `pd_le_1[i] or id_le_1[i]`), an INDEPENDENT cross-check of the complement-empty route,
  and the honest **shod/complement witness** (`{module, pd_le_1, id_le_1}` for a complement
  module).
- `atlas.projective_placement` / `atlas.injective_placement` (`{v → "L"|"R"|"both"|
  "neither"}`) → **ada** (all `≠ "neither"`), **quasi-tilted** (route 2: all projective
  placements ∈ {"L","both"}).
- `atlas.is_complete` / `atlas.status` → P61 refuses in lockstep (never a partial ladder).
- `atlas._modules` / the P41 `ar.arrows` (via a fresh `ar_quiver` when the atlas does not
  expose the arrows) → the **weakly-shod** AR-SCC sweep. *Cross-plan contract note:* P55's
  `LeftRightAtlas` carries `_modules` and `_leq` but **not** the raw AR arrows; P61 gets the
  arrows by calling `A.ar_quiver(budget_modules=budget)` once (cheap, the knit is cached
  within the atlas build path) **or** — the minimal P55 contract extension — reads a new
  `atlas._ar_arrows` field if P55 chooses to expose it. **Decision:** P61 calls
  `ar_quiver` itself (no P55 change required); the plan names `atlas._ar_arrows` only as the
  optional optimisation P55 could add. The `_index_in_U` helper (QQ-scope) is imported from
  `modules.left_right` (do not re-implement).

`Algebra.global_dimension()` (P05/shipped) and `Algebra.hochschild_cohomology(1)` (shipped)
complete the inputs. **No new AR/homology machinery.**

**P55 contract amendment (review fix M1+W3+C1 — a dated addendum is appended to P55's plan
doc, committed together).** P55 already computes, per indecomposable, exactly `pd X ≤ 1` and
`id X ≤ 1` (its `_pd_le_1` / `_id_le_1` sweep, `betti(2) == 0`) — it just does not currently
expose them. P55's `LeftRightAtlas` gains two **additive** index-aligned fields
`pd_le_1: tuple[bool,...]` and `id_le_1: tuple[bool,...]` (assigning the vectors P55's sweep
already builds — **zero new mathematics**, no golden/key change if kept engine-internal).
This makes the shod QT2 route (§Task A) a genuinely independent computation (not a tautology
of complement-empty) and gives the shod/ada witnesses honest per-module pd/id flags. Until
the amendment lands, P61 falls back to re-running `_pd_le_1`/`_id_le_1` off `atlas._modules`
(the same helpers, imported from `modules.left_right`) — the QT2 route works either way; the
amendment is the clean, no-recompute path. **The `is_complete` placement contract (review fix
C2):** `atlas.projective_placement` / `injective_placement` are `dict`s (populated) **whenever
`atlas.is_complete`** and `None` only on the refusal path; `_ada` asserts they are populated
when `is_complete` and raises a typed `QuiverlabError` otherwise (a P55-contract violation,
never a `NoneType` crash).

---

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **P55 (left/right parts) is the hard prerequisite** — `modules/left_right.py::
  left_right_parts`, `LeftRightAtlas` (fields `left`, `right`, `intersection`, `complement`,
  `projective_placement`, `injective_placement`, `universe_size`, `is_complete`, `status`,
  `note`, `_modules`, `_leq`, **+ the P61 contract amendment `pd_le_1` / `id_le_1`**), and
  `_index_in_U` — merged on `dev` **before** this branch (the amendment addendum is appended to
  P55's doc and committed with this plan). Also shipped: `Algebra.global_dimension() ->
  GlobalDimension` (`.value: int | None`, `.exact: bool`), `Algebra.hochschild_cohomology(1,
  verbose=False, trace=False)` (`HHTable`, `[1]` = `dim HH¹`; the trace sink verified in the
  signature), `Algebra.ar_quiver(budget_modules=...) -> ARQuiver` (`.arrows`, `.is_complete`,
  `.status`), `modules.hom.end_dim`, `A.domain` (the coefficient domain — no `A.field`) +
  the new additive `Domain.is_algebraically_closed` (Task B Step 0). Branch
  `plan-61-recognizer-ladder` off
  `dev` **after P55 merges** (metaplan Wave 2; **needs P55**; sibling P60 also needs P55 —
  they may branch in parallel and merge sequentially, goldens/i18n/bib merge-order care).
- **Honest semi-decision contract (metaplan §1.3; inherited from P55/P41).**
  `recognizer_ladder` is complete **iff** `A` is representation-finite and non-self-injective
  (the P55 atlas is). A rep-infinite (`status="budget"`), self-injective
  (`status="unsupported"`), or knit-erroring (`status="error"`) `A` yields
  `is_complete=False` with that `status`/`note` and **no rungs** — never a partial ladder.
  `RecognizerLadder.is_complete`/`status` mirror `LeftRightAtlas` exactly. **ACLV Example
  2.2(c)** (rep-infinite ada `1⇉2⇉3⇉4`, rad²=0) is the honest boundary oracle: mathematically
  ada, refused because rep-infinite (`status="budget"`).
- **Field scope is load-bearing (M1 + Theorem-B gate).** The rungs (parts, complement,
  placements, gl.dim, AR-SCC) are exact over **every** `Domain`, but every step that
  *identifies* a module in the universe (`_index_in_U` for placements / injectives /
  projectives, imported from P55) calls `is_isomorphic`, which is **char-0 decisive** but
  **positive-only over large GF(p)/GF(pⁿ)** (raises when it cannot exhibit an isomorphism).
  So **all identification-touching batteries run over QQ** (decisive negatives), inheriting
  P55's M1 rule; a `GF(32003)` parity spot-check is kept only for the distinct-dim-vector
  `kA_n` family. The **ada/HH¹ Theorem-B verdict** additionally requires the field to be
  **algebraically closed** — emitted only when `A.domain.is_algebraically_closed` (the new
  additive flag, `True` only on the CC working domain; **never** the `CC` sentinel — there is
  no `Algebra.field`, and QQi shares CC's `SympyExactDomain` class, §"the algebraically-closed
  gate"); over QQ/GF(p)/GF(pⁿ)/QQi the block reports `dim HH¹` + the theorem note, no verdict.
- **`laura` is trivially `True` in scope — the honest framing.** Every representation-finite
  algebra is laura (finite `ind A` ⟹ finite complement). P61 reports `laura.verdict = True`
  **with the complement as the useful datum** (the interesting laura content is rep-infinite,
  out of knit scope; ACLV Example 2.2(c)). The `laura` rung's `note` states this honestly;
  the verification page carries it as an honest-scope entry. Cite `smith_almost_laura` for
  the laura landscape (not for the elementary triviality).
- **`tilted` is NOT a P61 rung.** `tilted ⟹ quasi-tilted` is the top of the nesting, but
  tiltedness is P60's `tilted_check` (R17). P61 ships a skipped auto-flipping
  `test_tilted_implies_quasi_tilted_PIN` (the P55 `test_support_components_are_tilted_PIN`
  precedent) that becomes a real assert when P60 lands.
- **No floats in `src/`.** Verdicts Boolean; complement/placements index/name sets;
  `hh1_dim`/`gl.dim` `int`. The only float conversion is client-side (`docs/gui/gui.js`,
  exempt).
- **Composition is left-to-right** (`a*b` = first `a` then `b`). All refusals are
  `QuiverlabError(message, hint=...)` (the two-space `[hint: ...]` convention). Monotone-
  consistency and theorem-gate violations raise `QuiverlabError` (loud internal error),
  never a silently-wrong verdict.
- **Plan-32 markers.** Monotone-consistency (nesting) + the **guard-RAISES** test + theorem-gate
  (gl.dim bounds) + the `_theorem_b_verdict` False-branch + the alg-closed gate predicate +
  ada/HH¹ "reduces to base field" corollary = `oracle_selfcert`; the **ACLV Example 2.2(b)
  five-verdict pin**, hereditary ⇒ all-five-`True` (`kA_n`), the **`kupisch=[3,2,2]` weakly-shod
  =False** pin, ACLV Example 2.2(c) rep-infinite-ada honest refusal, the ada/HH¹=0 tree fixture
  over CC = `oracle_literature`; the **shod two-routes** (complement-empty ≡ QT2 sweep), the
  **quasi-tilted two-routes agreement** (HRS `gl.dim≤2 ∧ complement=[]` ≡ ACLV
  `every P_x ∈ 𝓛_A`), the **weakly-shod AR-SCC ≡ Hom-closure** route tie, and the
  **P56 `is_simply_connected` agreement** on ada/CC members = `oracle_crossengine`; the QPA
  pointwise pd/id + gl.dim crosscheck = the `qpa` bucket. `oracle_*` markers FORBIDDEN in
  `tests/{webapp,gui,hpc}` (`test_oracle_classes.py`) — engine batteries in `tests/modules/`,
  cross-runner twins in `tests/webapp/` unmarked.
- **Mid-merge-train counts.** v1.0.0 lands ~29 subplans in waves. **Task F recounts the
  oracle-class table at merge time** by running `tests/release/test_oracle_classes.py` (paste
  the live numbers, never a guessed-at-authoring count) and claims only this plan's deltas.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table
  green) and adds its citations to `citations/references.bib` + `registry.py` (`bibtex()`
  hard-fails if the two disagree). Conventional commits; green tests at every commit.

**M-resolution (no red commits, no forward references).** The `RecognizerLadder`,
`LadderRung`, `SimpleConnectedness` dataclasses, `recognizer_ladder`, and the `Algebra`
delegates are assembled in **Task A**, with the `SimpleConnectedness` block defaulting to
`SimpleConnectedness(applicable=False, field_algebraically_closed=False, hh1_dim=None,
verdict=None, theorem=<the statement>, note="not computed")`. **Task B** populates it (the
ada/HH¹ verdict + P56 agreement). Because the field already exists (as a safe default) from
Task A, the Task-B tests reference a real attribute — every task boundary is green, no
forward reference, no `xfail` fence (the P55 M3 precedent).

---

### Task 0: citations (do first — later tasks reference the keys)

**Files:** Modify `src/quiverlab/citations/references.bib`,
`src/quiverlab/citations/registry.py`. Test: `tests/citations/test_bib_structure.py`
(existing gate — `bibtex(key)` resolvable, registry ↔ bib in sync).

**Interfaces:** the `_r(key, bibtex_key, kind, title, annotation, *tags)` helper
(`registry.py:24`). **Reuse (do NOT re-add):** `aclv_supports` (→ `ACLV2011`, ada +
Theorem B) and `organising_module_category` (→ `AACV2021`, survey) — both added by P55.

- [ ] **Step 1: add the new keys** (BibTeX-verify each before committing — Plan-29 rule):

```bibtex
@book{HRS1996,
  author    = {Happel, Dieter and Reiten, Idun and Smal{\o}, Sverre O.},
  title     = {Tilting in Abelian Categories and Quasitilted Algebras},
  series    = {Memoirs of the American Mathematical Society},
  volume    = {120},
  number    = {575},
  publisher = {American Mathematical Society},
  year      = {1996},
}
@article{CL2003weaklyshod,
  author  = {Coelho, Fl{\'a}vio U. and Lanzilotta, Marcelo A.},
  title   = {Weakly shod algebras},
  journal = {Journal of Algebra},
  volume  = {265},
  number  = {2},
  pages   = {379--403},
  year    = {2003},
}
@article{AC2003laura,
  author  = {Assem, Ibrahim and Coelho, Fl{\'a}vio U.},
  title   = {Two-sided gluings of tilted algebras},
  journal = {Journal of Algebra},
  volume  = {269},
  number  = {2},
  pages   = {456--479},
  year    = {2003},
}
@misc{Smith2007almostlaura,
  author        = {Smith, David},
  title         = {Almost laura algebras},
  year          = {2007},
  eprint        = {math/0702562},
  archivePrefix = {arXiv},
  primaryClass  = {math.RT},
}
@article{BFT2017quasitilted,
  author  = {Bordino, Natalia and Fern{\'a}ndez, Elsa and Trepode, Sonia},
  title   = {On the quiver with relations of a quasitilted algebra and applications},
  journal = {Communications in Algebra},
  volume  = {45},
  number  = {9},
  pages   = {4050--4061},
  year    = {2017},
}
```

  and in `registry.py`:

```python
_r("hrs_quasitilted", "HRS1996", "foundation",
   "Tilting in Abelian Categories and Quasitilted Algebras",
   "Happel-Reiten-Smalo: quasi-tilted = (QT1) gl.dim <= 2 AND (QT2) every indec pd<=1 or "
   "id<=1; QT2 alone => gl.dim <= 3. The definitional ground truth for Plan 61's "
   "quasi-tilted rung.", "recognizer"),
_r("coelho_lanzilotta_weakly_shod", "CL2003weaklyshod", "foundation",
   "Weakly shod algebras",
   "Coelho-Lanzilotta: shod = every indec pd<=1 or id<=1 (=> gl.dim <= 3); weakly shod = "
   "bounded irreducible-morphism paths from an injective to a projective (WSA). Plan 61's "
   "shod + weakly-shod rungs.", "recognizer"),
_r("assem_coelho_laura", "AC2003laura", "foundation",
   "Two-sided gluings of tilted algebras",
   "Assem-Coelho: laura = ind A minus (L_A u R_A) is finite. Plan 61's laura rung (trivially "
   "true in representation-finite scope; the finite complement is the reported datum).",
   "recognizer"),
_r("smith_almost_laura", "Smith2007almostlaura", "foundation",
   "Almost laura algebras",
   "Smith: the almost-laura generalisation of laura algebras -- context for the laura "
   "landscape; Plan 61 cites it for the class, not for the elementary rep-finite triviality.",
   "recognizer"),
_r("bft_quasitilted_quiver", "BFT2017quasitilted", "foundation",
   "On the quiver with relations of a quasitilted algebra and applications",
   "Bordino-Fernandez-Trepode: the quiver-with-relations structure of quasitilted algebras "
   "-- a secondary reference on Plan 61's quasi-tilted rung.", "recognizer"),
```

  **Spec-ambiguity resolution (recorded):** BibTeX keys `HRS1996` / `CL2003weaklyshod` /
  `AC2003laura` / `Smith2007almostlaura` / `BFT2017quasitilted`; snake-case registry keys as
  above (house convention). HRS ships as `@book` (the Memoirs volume). `Smith2007almostlaura`
  ships as the **verified arXiv `@misc`** (worker upgrades to `@article` iff BibTeX-
  verifiable at merge). `aclv_supports` (ada / Theorem B) + `organising_module_category`
  (survey) reused from P55 — verify they are present (`bibtex("aclv_supports")`) and do NOT
  duplicate. The shod-origin Manuscripta Math. 100 (1999) entry is a `# PIN` secondary the
  worker may add; the recorded `CL2003weaklyshod` covers the shod definition (it recalls it).
- [ ] **Step 2:** `... -m pytest tests/citations/test_bib_structure.py -q` green.
- [ ] **Step 3: Commit** — `docs(citations): recognizer-ladder references (HRS Mem AMS 575, Coelho-Lanzilotta weakly-shod, Assem-Coelho laura, Smith almost-laura, Bordino-Fernandez-Trepode)`

---

### Task A: `recognizers_ladder.py` — the five rungs + nesting/theorem self-certs + delegates

**Files:**
- Create: `src/quiverlab/modules/recognizers_ladder.py`
- Modify: `src/quiverlab/core/algebra.py` (thin delegates beside `left_right_parts`)
- Test: `tests/modules/test_recognizer_ladder.py`

**Interfaces:**
- Consumes: `modules.left_right.left_right_parts`, `modules.left_right._index_in_U`,
  `Algebra.global_dimension`, `Algebra.ar_quiver`.
- Produces:
  ```python
  @dataclass(frozen=True)
  class LadderRung:
      name: str            # "quasi_tilted" | "shod" | "weakly_shod" | "laura" | "ada"
      verdict: bool
      certificate: dict    # the positive datum (e.g. {"gldim":1,"complement_size":0})
      witness: dict | None # populated only on verdict=False (e.g. {"module":"S_3","pd":2,"id":2})

  @dataclass(frozen=True)
  class SimpleConnectedness:
      applicable: bool
      field_algebraically_closed: bool
      hh1_dim: int | None
      verdict: bool | None
      theorem: str
      note: str

  @dataclass(frozen=True)
  class RecognizerLadder:
      algebra: object
      rungs: dict                 # {name -> LadderRung}
      complement: tuple           # ({"index","name","dimvec"}, ...)  -- the laura datum
      simple_connectedness: SimpleConnectedness
      gldim: int | None
      universe_size: int
      is_complete: bool
      status: str
      note: str
      def verdict(self, name) -> bool: ...          # rungs[name].verdict

  def _shod(atlas) -> LadderRung
      # verdict = (complement == []); CROSS-CHECK the QT2 sweep every indec pd_le_1|id_le_1
      # (must equal complement==[] -- survey Thm 4.1 a<=>b, loud QuiverlabError if they differ).
      # witness on False = {"module","pd_le_1","id_le_1"} for a complement module (QT2 violator).
  def _quasi_tilted(atlas, gld) -> LadderRung         # gld<=2 AND complement==[]; xcheck _quasi_tilted_aclv_route
  def _quasi_tilted_aclv_route(A, atlas=None) -> bool # every P_x placement in {"L","both"} (ACLV)
  def _weakly_shod(A, atlas, budget) -> LadderRung    # PRIMARY: AR-arrow SCC bounded-path sweep
  def _weakly_shod_hom(atlas) -> LadderRung           # CROSS-CHECK: atlas._leq reach + end_dim>=2 cyclicity
  def _weakly_shod_from_digraph(n, arrows, inj, proj) -> tuple[bool, dict|None]  # pure SCC helper
  def _laura(atlas) -> LadderRung                     # always True (rep-finite); complement datum
  def _ada(atlas) -> LadderRung
      # every projective_placement/injective_placement != "neither"; ASSERTS the placements are
      # populated when is_complete (C2: typed QuiverlabError if None, never a NoneType crash).
  def _assert_nesting(rungs) -> None                  # monotone self-cert, loud QuiverlabError
  def _assert_theorem_gates(rungs, gld) -> None       # shod=>gld<=3, qt=>gld<=2, loud QuiverlabError
  def recognizer_ladder(A, *, budget=256) -> RecognizerLadder
  ```
  `Algebra.recognizer_ladder(budget=256)` + the five `is_*` delegates lazy-import
  `quiverlab.modules.recognizers_ladder`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_recognizer_ladder.py
"""The recognizer ladder quasi-tilted/shod/weakly-shod/laura/ada (Plan 61 / R18).
Literature: the ACLV Example 2.2(b) rad^2=0 A5 gives (qt,shod,wshod,laura,ada) =
(False,False,True,True,True) with complement={S3} -- a full five-way discrimination; a
hereditary kA_n gives all five True. Self-cert: monotone nesting (qt=>shod=>wshod=>laura,
qt/shod=>ada, ada=>laura), theorem gates (shod=>gldim<=3, qt=>gldim<=2). Cross-engine: the
two quasi-tilted routes (HRS gldim/complement vs ACLV projectives-in-left) agree. QQ-scope
(M1)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.recognizers_ladder import recognizer_ladder

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _radsq_nakayama_a5():
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


@lit
def test_aclv_22b_full_ladder():
    L = recognizer_ladder(_radsq_nakayama_a5())
    assert L.is_complete and L.universe_size == 9 and L.gldim == 4
    assert L.verdict("quasi_tilted") is False
    assert L.verdict("shod") is False
    assert L.verdict("weakly_shod") is True
    assert L.verdict("laura") is True
    assert L.verdict("ada") is True
    assert [r["name"] for r in L.complement] == ["S_3"]           # the laura datum
    assert L.rungs["shod"].witness["module"] == "S_3"             # pd 2, id 2 witness


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_hereditary_all_five_true(n):
    L = recognizer_ladder(linear_path_algebra(n, field=QQ))     # kA_n hereditary
    assert all(L.verdict(k) for k in
               ("quasi_tilted", "shod", "weakly_shod", "laura", "ada"))
    assert L.complement == ()                                    # empty complement


@xeng
def test_shod_two_routes_agree():
    # Thm 4.1 (a)<=>(b): the QT2 sweep (every indec pd<=1 or id<=1) == complement empty.
    # NON-tautological: QT2 reads atlas.pd_le_1/id_le_1 (independent of the L/R closure).
    # The monomial square (rep-finite, gldim 2) is NOT shod (S_3 lies in the complement? no --
    # its complement is non-empty; verified live: proj P_1 in neither, so complement != []);
    # kA3 (empty complement) IS shod.
    from quiverlab.modules.recognizers_ladder import _shod
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    for A, expect in ((Q.algebra(relations=["a*b"], field=QQ), False),        # monomial: not shod
                      (linear_path_algebra(3, field=QQ), True)):              # kA3: shod
        L = recognizer_ladder(A)
        assert L.verdict("shod") is expect
        assert (len(L.complement) == 0) == L.verdict("shod")        # complement-empty route
        # the QT2 route is a real independent computation inside _shod (asserted equal there);
        # a mismatch would have already raised QuiverlabError. Re-affirm the witness shape:
        if not L.verdict("shod"):
            w = L.rungs["shod"].witness
            assert set(w) == {"module", "pd_le_1", "id_le_1"} and not (w["pd_le_1"] or w["id_le_1"])


@xeng
def test_quasi_tilted_two_routes_agree():
    # HRS route (gldim<=2 AND complement==[]) == ACLV route (every P_x in L_A). PIN the
    # literature verdicts: comm-square-with-relation quasi-tilted True; rad^2=0 A5 False.
    from quiverlab.modules.recognizers_ladder import _quasi_tilted_aclv_route
    comm = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
                  ).algebra(relations=["a*b - c*d"], field=QQ)
    for A, expect in ((_radsq_nakayama_a5(), False), (comm, True)):
        L = recognizer_ladder(A)
        assert L.verdict("quasi_tilted") is expect                 # pinned literature verdict
        assert L.verdict("quasi_tilted") == _quasi_tilted_aclv_route(A)   # HRS route == ACLV route


@xeng
def test_weakly_shod_two_routes_agree():
    # PRIMARY (AR-arrow SCC) == CROSS-CHECK (Hom-closure atlas._leq + end_dim>=2). Non-vacuous:
    # the fixture is directed (both True); kupisch=[3,2,2] is one big SCC (both False).
    from quiverlab import NakayamaAlgebra
    from quiverlab.modules.recognizers_ladder import _weakly_shod, _weakly_shod_hom
    from quiverlab.modules.left_right import left_right_parts
    for A in (_radsq_nakayama_a5(), NakayamaAlgebra(kupisch=[3, 2, 2], field=QQ)):
        atlas = left_right_parts(A)
        assert _weakly_shod(A, atlas, 256).verdict == _weakly_shod_hom(atlas).verdict


@lit
def test_weakly_shod_false_on_kupisch_322():
    # M2: a REAL non-directed knittable algebra. kupisch=[3,2,2] knits complete, its whole
    # module category is one non-trivial SCC on an injective->projective route => weakly_shod
    # False (verified live; also fully degenerate -- ada/shod/quasi_tilted all False, laura True).
    from quiverlab import NakayamaAlgebra
    L = recognizer_ladder(NakayamaAlgebra(kupisch=[3, 2, 2], field=QQ))
    assert L.is_complete and L.universe_size == 7 and L.gldim == 3
    assert L.verdict("weakly_shod") is False
    assert L.rungs["weakly_shod"].witness is not None            # {Z, inj_source, proj_target}
    assert L.verdict("laura") is True and L.verdict("ada") is False and L.verdict("shod") is False


@selfcert
def test_monotone_nesting_holds_on_ladder():
    for A in (_radsq_nakayama_a5(), linear_path_algebra(4, field=QQ)):
        L = recognizer_ladder(A)
        v = L.verdict
        assert (not v("quasi_tilted")) or v("shod")
        assert (not v("shod")) or v("weakly_shod")
        assert (not v("weakly_shod")) or v("laura")
        assert (not v("quasi_tilted")) or v("ada")
        assert (not v("shod")) or v("ada")
        assert (not v("ada")) or v("laura")


@selfcert
def test_theorem_gates_gldim():
    for A in (_radsq_nakayama_a5(), linear_path_algebra(3, field=QQ)):
        L = recognizer_ladder(A)
        if L.verdict("shod"):
            assert L.gldim <= 3
        if L.verdict("quasi_tilted"):
            assert L.gldim <= 2


@selfcert
def test_selfinjective_refused_loudly():
    A = truncated_polynomial(3, field=QQ)                       # k[x]/(x^3), self-injective
    L = recognizer_ladder(A)
    assert L.is_complete is False and L.status == "unsupported"
    assert L.rungs == {}                                        # never a partial ladder


@selfcert
def test_rep_infinite_ada_22c_refused_loudly():
    # ACLV Example 2.2(c): 1=>2=>3=>4 bound by rad^2=0 -- mathematically ada, refused
    # because representation-infinite (the honest boundary; ada does NOT imply laura here).
    Q = Quiver([1, 2, 3, 4],
               {"a": (1, 2), "b": (1, 2), "c": (2, 3), "d": (2, 3), "e": (3, 4), "f": (3, 4)})
    A = RadicalSquareZero(Q, field=QQ)
    L = recognizer_ladder(A, budget=40)
    assert L.is_complete is False and L.status in ("budget", "error", "unsupported")


@selfcert
def test_delegates_match_free_function():
    A = _radsq_nakayama_a5()
    assert A.is_ada() is True and A.is_shod() is False
    assert A.is_quasi_tilted() is False and A.is_weakly_shod() is True and A.is_laura() is True
```

  Plus a **synthetic weakly-shod-False unit test** (no algebra — feed a hand-built AR digraph
  with an injective → 2-cycle → projective and assert the SCC sweep returns `False` + a
  witness; guarantees the `False` branch is covered even if no rep-finite non-directed algebra
  is pinned):

```python
@selfcert
def test_weakly_shod_sweep_detects_cycle_between_inj_and_proj():
    from quiverlab.modules.recognizers_ladder import _weakly_shod_from_digraph
    # 5 vertices: 0 = injective, 1,2 = a 2-cycle, 3 = projective, 4 = isolated.
    arrows = {(0, 1): 1, (1, 2): 1, (2, 1): 1, (2, 3): 1}
    ok, witness = _weakly_shod_from_digraph(n=5, arrows=arrows, inj={0}, proj={3})
    assert ok is False and set(witness["cycle"]) == {1, 2}
    # remove the injective->cycle edge: now bounded => weakly shod.
    ok2, _ = _weakly_shod_from_digraph(n=5, arrows={(1, 2): 1, (2, 1): 1, (2, 3): 1},
                                       inj={0}, proj={3})
    assert ok2 is True
```

  Plus a **guard-firing unit test** (M3 — the monotone self-cert RAISES on inconsistent rungs;
  mirrors the digraph unit-test pattern, no algebra needed):

```python
@selfcert
def test_monotone_guard_raises_on_inconsistent_rungs():
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.recognizers_ladder import _assert_nesting, LadderRung
    def rung(name, v): return LadderRung(name, v, {}, None if v else {"x": 1})
    good = {n: rung(n, v) for n, v in
            [("quasi_tilted", False), ("shod", False), ("weakly_shod", True),
             ("laura", True), ("ada", True)]}
    _assert_nesting(good)                                        # consistent: no raise
    bad = dict(good, shod=rung("shod", True))                    # shod True but weakly_shod... ok;
    bad["weakly_shod"] = rung("weakly_shod", False)              # shod True, weakly_shod False: BAD
    with pytest.raises(QuiverlabError, match="nesting|monoton|consist"):
        _assert_nesting(bad)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError:
  quiverlab.modules.recognizers_ladder`.
- [ ] **Step 3: Implement.** `recognizer_ladder`: one `left_right_parts(A, budget=budget)`;
  on `not atlas.is_complete` return `RecognizerLadder(A, {}, (), <default SC block>, None, 0,
  False, atlas.status, atlas.note)`. Else `gld = A.global_dimension().value`, build the five
  rungs (`_shod`, `_quasi_tilted`, `_weakly_shod`, `_laura`, `_ada`), `_assert_nesting`,
  `_assert_theorem_gates`, and assemble (the `SimpleConnectedness` field keeps its "not
  computed" default — Task B fills it; M-resolution). Expose `_quasi_tilted_aclv_route(A)`
  (the ACLV `every P_x ∈ 𝓛_A` route) and `_weakly_shod_from_digraph(n, arrows, inj, proj)`
  (the pure SCC helper the synthetic test drives). Refuse loudly (do not catch) if
  `is_isomorphic` raises off char-scope (M1) — the honest whole-compute refusal.
- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(modules): recognizer ladder (quasi-tilted/shod/weakly-shod/laura/ada) off the P55 atlas + monotone nesting & gl.dim theorem-gate self-certs; ACLV 2.2(b) (F,F,T,T,T)`

---

### Task B: the ada / HH¹ simple-connectedness block (ACLV Theorem B) + P56 agreement

**Files:**
- Modify: `src/quiverlab/fields/domain.py` (base `Domain.is_algebraically_closed = False`) +
  `src/quiverlab/fields/complexfield.py` (`SympyExactDomain.__init__` gains
  `algebraically_closed=False`; `ComplexField.make_domain` passes `True`; `GaussianRationalField`
  /QQi keeps the default `False`) — the additive alg-closed flag (§"the algebraically-closed
  gate")
- Modify: `src/quiverlab/modules/recognizers_ladder.py` (populate `simple_connectedness`;
  add `_theorem_b_verdict`)
- Test: `tests/modules/test_recognizer_ladder_hh1.py`, `tests/modules/test_ada_hh1_p56_agreement.py`
  (+ the gate predicate on constructed algebras, in `test_recognizer_ladder_hh1.py`)

**Interfaces:** `Algebra.hochschild_cohomology(1, verbose=False, trace=False)` (`HHTable`,
`[1] = dim HH¹`; the trace sink keeps the ladder from writing `Worked steps` files);
`A.domain.is_algebraically_closed` (the new additive alg-closed flag — the gate);
`quiverlab.invariants.coverings.is_simply_connected` (P56 — behind `importorskip`).

- [ ] **Step 0: the alg-closed flag.** Add `Domain.is_algebraically_closed` (default `False`)
  on the base `Domain`; thread `algebraically_closed` through `SympyExactDomain.__init__` and
  stamp `True` in `ComplexField.make_domain` only (QQi keeps `False`). Verify the mandate test
  (`test_alg_closed_gate_predicate_on_constructed_algebras`) is red before, green after; the
  `no_floats` AST gate is unaffected (a boolean, no float). Commit separately:
  `feat(fields): Domain.is_algebraically_closed (True only on the CC working domain; QQi stays False)`.
- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_recognizer_ladder_hh1.py
"""ACLV Theorem B: ada over an algebraically closed field => (simply connected <=> HH^1=0),
and HH.(A) reduces to k. Plan 61's headline oracle. Over CC the verdict is emitted; over
QQ/GF(p)/QQi only dim HH^1 + the theorem note (Theorem B's alg-closed hypothesis unmet).
The gate is A.domain.is_algebraically_closed -- there is NO Algebra.field, and QQi shares
CC's SympyExactDomain class, so the flag is the ONLY sound predicate."""
import pytest

from quiverlab import GF, Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ, CC, QQi
from quiverlab.modules.recognizers_ladder import recognizer_ladder

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _radsq_nakayama_a5(field):
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=field)


@selfcert
def test_alg_closed_gate_predicate_on_constructed_algebras():
    # H1 MANDATE: assert the gate on A.domain (the constructed algebra), NEVER the CC sentinel.
    # CC declares C => True; QQ/GF(p)/QQi are not algebraically closed => False (QQi shares CC's
    # SympyExactDomain class, so the flag -- not the type -- must distinguish them).
    assert linear_path_algebra(2, field=CC).domain.is_algebraically_closed is True
    for f in (QQ, GF(7), QQi):
        assert linear_path_algebra(2, field=f).domain.is_algebraically_closed is False


@lit
def test_fixture_simply_connected_over_CC():
    # rad^2=0 A5: ada True, tree quiver => simply connected => HH^1 = 0 (Theorem B).
    L = recognizer_ladder(_radsq_nakayama_a5(CC))
    sc = L.simple_connectedness
    assert sc.applicable is True and sc.field_algebraically_closed is True
    assert sc.hh1_dim == 0 and sc.verdict is True          # simply connected


@selfcert
def test_theorem_B_verdict_false_when_hh1_nonzero():
    # UNCONDITIONAL coverage of the verdict=False branch via the pure helper (a rep-finite
    # ada NOT-simply-connected instance is elusive -- ada strongly pushes rep-finite algebras
    # toward simple connectedness, see the in-plan note -- so the False branch is exercised
    # synthetically, mirroring the weakly-shod digraph unit test).
    from quiverlab.modules.recognizers_ladder import _theorem_b_verdict
    assert _theorem_b_verdict(ada=True, alg_closed=True, hh1_dim=1) is False
    assert _theorem_b_verdict(ada=True, alg_closed=True, hh1_dim=0) is True
    assert _theorem_b_verdict(ada=True, alg_closed=False, hh1_dim=0) is None   # gated off
    assert _theorem_b_verdict(ada=False, alg_closed=True, hh1_dim=0) is None   # not ada


@selfcert
def test_verdict_gated_off_non_algebraically_closed_field():
    # over QQ (not alg-closed): applicable False, hh1_dim reported, verdict None + theorem note.
    L = recognizer_ladder(_radsq_nakayama_a5(QQ))
    sc = L.simple_connectedness
    assert sc.field_algebraically_closed is False and sc.applicable is False
    assert sc.hh1_dim == 0 and sc.verdict is None and "HH" in sc.theorem


@lit
def test_not_ada_no_verdict():
    # A REAL in-scope non-ada witness (verified live): the triangle 1->2->3 + 1->3 with the
    # long route killed (I = <a*b>). Its projective P_1 lies in NEITHER part (P_1 has id >= 2
    # so P_1 not in R_A, and a bad-pd predecessor keeps it out of L_A), so ada = False; hence
    # Theorem B does not apply and the SC block is not applicable.
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3)}).algebra(
        relations=["a*b"], field=CC)
    L = recognizer_ladder(A)
    assert L.is_complete and L.verdict("ada") is False
    assert L.simple_connectedness.applicable is False and L.simple_connectedness.verdict is None
```

```python
# tests/modules/test_ada_hh1_p56_agreement.py
"""Plan 61 x Plan 56 interplay: on ada algebras over CC, the Theorem-B HH^1 verdict AGREES
with P56's is_simply_connected verdict wherever P56 is definite; where P56 returns None
(Adian-Rabin undecidable), Theorem B RESOLVES it. Auto-flips to a real assert when P56
lands."""
import pytest

pytest.importorskip("quiverlab.invariants.coverings")   # P56 alignment point
from quiverlab.invariants.coverings import is_simply_connected  # noqa: E402
from quiverlab import Quiver, RadicalSquareZero
from quiverlab.fields import CC
from quiverlab.modules.recognizers_ladder import recognizer_ladder

pytestmark = pytest.mark.oracle_crossengine


def _ada_cc_members():
    # two ada algebras over CC (both simply connected here): the tree fixture + the comm-square.
    Q5 = Quiver([1, 2, 3, 4, 5],
                {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    comm = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
                  ).algebra(relations=["a*b - c*d"], field=CC)
    return [RadicalSquareZero(Q5, field=CC), comm]


@pytest.mark.parametrize("A", _ada_cc_members())
def test_theorem_B_agrees_with_p56_on_ada_cc_members(A):
    # ada/CC members: Theorem-B HH^1 verdict AGREES with P56 wherever P56 is definite;
    # where P56 returns None (Adian-Rabin), Theorem B is the complete resolution (headline).
    L = recognizer_ladder(A)
    assert L.verdict("ada") is True and L.simple_connectedness.applicable is True
    p56 = is_simply_connected(A).verdict
    if p56 is not None:                                     # P56 definite => must match
        assert p56 == L.simple_connectedness.verdict
    else:                                                   # P56 undecided => P61 resolves it
        assert L.simple_connectedness.verdict is not None
```

- [ ] **Step 2: Run to verify failure** — `simple_connectedness` is the Task-A default
  ("not computed", `applicable=False`, `hh1_dim=None`), so `test_fixture_simply_connected_
  over_CC` fails (`hh1_dim is None != 0`) — the red that drives Step 3 (no forward reference:
  the attribute is real, just a default).
- [ ] **Step 3: Implement `_simple_connectedness(A, ada_verdict)`** + the pure decision helper
  `_theorem_b_verdict(ada, alg_closed, hh1_dim)` (the synthetic-test seam: returns `None` if
  `not ada`, `None` if `not alg_closed`, else `hh1_dim == 0`). In `_simple_connectedness`:
  `alg_closed = A.domain.is_algebraically_closed` (**the new flag off `A.domain` — never
  `isinstance(..., ComplexField)`, never `A.field`**; §"the algebraically-closed gate"). If
  `not ada_verdict` → default not-applicable block (`hh1_dim` may be left `None`). Else compute
  `hh1 = A.hochschild_cohomology(1, verbose=False, trace=False)[1]` — **pass the trace/verbose
  sink so the ladder does not spew `Worked steps: quiverlab_traces/*.html` files** (signature
  verified read-only: `hochschild_cohomology(top, max_cells=..., engine=..., auto_cs=...,
  verbose=None, trace=None)`; `verbose=False` suppressed the trace file in the venv — the
  worker confirms which of `verbose`/`trace` gates the file write and sets both quiet); honor
  the honest refusal if the HH engine raises (record `hh1_dim=None`, `verdict=None`,
  `note="HH^1 unavailable: <message>"`, never a crash — the Plan-30 per-entry precedent).
  `verdict = _theorem_b_verdict(ada_verdict, alg_closed, hh1)`; `applicable = (ada_verdict and
  alg_closed)`; `note` includes the "reduces to base field" corollary on `verdict is True`,
  and the Theorem B statement + "recompute over CC" on the gated-off path. The theorem string
  is a fixed constant citing `aclv_supports`.
- [ ] **Step 4: Run Tasks A+B, verify pass.**
- [ ] **Step 5: Commit** — `feat(modules): ada/HH^1 simple-connectedness block (ACLV Thm B, alg-closed CC gate) + P56 is_simply_connected agreement oracle`

---

### Task C: literature + honest-scope oracle battery

**Files:** Test only — `tests/modules/test_recognizer_ladder_oracles.py` (no `src/` change).

Covers (beyond Tasks A/B): the **strict-weakly-shod `HH^{≥2}=0`** bonus oracle (survey §5.1
[43]: strict weakly shod ⟹ `Hⁱ(A)=0` for `i ≥ 2` — assert on a weakly-shod-not-quasi-tilted
in-scope algebra such as the fixture: `dim HH²(A) == 0`; **but** the survey warns quasi-tilted
algebras CAN have `HH² ≠ 0`, so this is asserted only on the *strict* (non-quasi-tilted)
weakly-shod fixture — the fixture qualifies: weakly-shod `True`, quasi-tilted `False`); the
`ada ⟹ gl.dim ≤ 4` and `ada ⟹ (pd ≤ 2 or id ≤ 1)` **Corollary 2.6** self-certs (over the
in-scope zoo members that are ada); and the `# PIN` non-directed weakly-shod-`False` algebra
if the worker finds one. Mark `oracle_literature` (the survey/ACLV pins) and `oracle_selfcert`
(the corollary gates). Verify the fixture's `HH² == 0` live before pinning.

- [ ] **Step 1–3:** write, run (verify the survey `HH^{≥2}=0` value live), commit —
  `test(modules): recognizer-ladder literature + Corollary-2.6/strict-weakly-shod HH oracles`

---

### Task D: GUI/webapp/report — the `recognizer_ladder` compute kind

`recognizer_ladder` is an **ALGEBRA-level compute kind** (computed on the drawn `A`, like
P55 `left_right_parts` / `ar_quiver` — routed through `_dispatch`, **NOT** `_dispatch_module`;
**schema v1**, no module block; sizes on `A.dim`). Carries a budget:
`compute.push("recognizer_ladder:256")`.

**Files (the P55 `left_right_parts` touchpoints, mirrored):**
- Modify: `src/quiverlab/modules/recognizers_ladder.py` (`recognizer_ladder_block(A, budget)`)
- Modify: `src/quiverlab/hpc/spec.py` (`parse_compute_item` budget interceptor; `_dispatch`
  branch catching the M1 `is_isomorphic` / `DepthLimitError` refusal into
  `{"kind":"recognizer_ladder","error":"<loud message>"}`; `_snip` recipe)
- Modify: `webapp/server/schema.py` (parse twin — budget interceptor) +
  `webapp/server/estimator.py` (add `"recognizer_ladder"` to the `_max_degree` skip tuple
  beside `left_right_parts`/`ar_quiver` — its budget is not a homological degree)
- Modify: `docs/gui/runner.py` (Pyodide twin: `_parse_compute` interceptor, `compute_one`
  branch calling the **same** `recognizer_ladder_block`, `calls` snippet,
  `ETA_MODEL["scalars"]["recognizer_ladder"] ~2.5`)
- Modify: `docs/gui/gui.js` (+ re-copy byte-identical to `webapp/static/gui/gui.js`): the
  checkbox + budget picker, the `el` id list, the `buildRequest` push, a `renderBlock` branch
  `renderRecognizerLadder(div, b)` (the ✓/✗/— ladder table + witness-per-✗ + the laura
  complement chips + the ada/HH¹ verdict line), `KIND_CTRL {cb, budget:true}`, `THEMES`
  (`"structure"`), `SEARCH_INDEX`
- Modify: `webapp/server/i18n/{en,es,fr,zh}.json` (`inv.recognizer_ladder`,
  `pick.kind.recognizer_ladder` ×4 + `block.recognizer_ladder.*` row labels the renderer
  reads: the five class names, `verdict_yes/no`, `witness`, `complement`, `simply_connected`,
  `hh1`, `theorem_note`)
- Modify: `src/quiverlab/trace/results_html.py` (`_HEADINGS["recognizer_ladder"] =
  "Recognizer ladder"` + a `_block_html` branch `_recognizer_ladder_html(b)` — a rendered
  ✓/✗/— table with the witness column + the simple-connectedness paragraph)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE new fixture `recognizer_ladder_kA3`, existing entries byte-identical FIRST; docstring
  bullet), `tests/webapp/test_layout_picker.py` `ALL_KINDS`
- Test: `tests/webapp/test_recognizer_ladder_p61.py`, `tests/gui/test_recognizer_ladder_twin.py`

**Block shape (returned identically by both runners — the byte-parity contract):**
```python
{"kind": "recognizer_ladder",
 "n": int,                                      # |Q_0|
 "complete": bool, "status": "complete"|"budget"|"unsupported"|"error",
 "universe_size": int | None, "gldim": int | None,
 "ladder": {
   "quasi_tilted": {"verdict": bool, "certificate": {...}, "witness": {...}|null},
   "shod":         {...}, "weakly_shod": {...}, "ada": {...},
   "laura":        {"verdict": true, "complement_size": int}},
 "complement": [ {"name": "S_3"|..., "dimvec": {...}} , ... ],   # the laura datum (may be non-empty)
 "simple_connectedness": {"applicable": bool, "field_algebraically_closed": bool,
                          "hh1_dim": int|null, "verdict": true|false|null,
                          "theorem": str, "note": str},
 "note": str | null,
 "references": ["hrs_quasitilted","coelho_lanzilotta_weakly_shod","assem_coelho_laura",
                "aclv_supports","smith_almost_laura","bft_quasitilted_quiver",
                "organising_module_category"],
 "citations": [...]}
```

- [ ] **Step 1:** failing cross-runner test (unmarked — extras-gated dir; copy the P55
  `test_left_right_parts_p55.py` runner-pair fixture): `kA3` → all five `True`, empty
  complement, `simple_connectedness.applicable` False over the default field (or True over CC
  in a second case); the rad²=0 A5 → the (F,F,T,T,T) ladder + `complement == [S_3]`;
  self-injective → `status="unsupported"`, empty ladder, no crash; `test_twin_parity`
  (`json.dumps(sort_keys=True)` equality via the shared `recognizer_ladder_block`).
- [ ] **Step 2:** implement `recognizer_ladder_block(A, budget)` — call `recognizer_ladder`,
  strip `Module`s to names+dim-vectors, serialize the rungs + SC block, stamp `references`;
  on `is_complete=False` return the same shape with empty ladder + `status`/`note` (the
  honest refusal, like `left_right_parts_block`). The HH¹ call inside the ladder already
  passes `verbose=False, trace=False` (Task B Step 3) so the compute kind writes **no** stray
  `Worked steps` trace files — the webapp/HPC runners manage their own trace bundle. **§5
  disconnected-input note:** ACLV's standing hypothesis is a **connected** algebra (Theorem B
  and the whole ladder are stated for connected artin algebras); a disconnected input is
  handled the way the P41 knit / P55 atlas handle it — the worker verifies the knit's
  disconnected behavior read-only and either (a) the atlas returns `is_complete=False` /
  refuses (P61 mirrors it, no rungs) or (b) P61 raises a typed `QuiverlabError("recognizer
  ladder requires a connected algebra", hint=...)` up front — **never** a per-component
  verdict silently reported as the whole-algebra verdict.
- [ ] **Step 3:** wire the compute kind (the `left_right_parts`-style touchpoints above).
- [ ] **Step 4:** add ONE golden `recognizer_ladder_kA3`; verify existing entries
  byte-identical FIRST (schema-v1 budget kind adds no request field); add the change-log
  bullet.
- [ ] **Step 5:** run the gates —
  `... -m pytest tests/webapp/test_recognizer_ladder_p61.py tests/webapp/test_runner_delegation.py tests/gui/test_recognizer_ladder_twin.py tests/webapp/test_i18n.py tests/webapp/test_layout_picker.py tests/webapp/test_draw_page.py tests/webapp/test_js_parses.py tests/hpc -q`
- [ ] **Step 6:** commit — `feat(gui,webapp,hpc,trace): recognizer_ladder compute kind -- five-rung ladder + laura complement + ada/HH^1 simple-connectedness, both runners, i18n x4, one golden`

---

### Task E: QPA pointwise cross-oracle (probe-first) + honest fallback

**Files:** Test — `tests/qpa/test_recognizer_ladder_qpa.py`.

**Interfaces:** the QPA session (`session.should_skip_qpa()`, `session.libgap_handle()`).
**QPA has no quasi-tilted / shod / weakly-shod / laura / ada recognizer** — but its
`ProjDimensionOfModule` / `InjDimensionOfModule` / `GlobalDimension` corroborate the
**pointwise** ingredients (the pd/id data behind the parts, and the gl.dim gates).

- [ ] **Step 1:** probe live QPA (the Plan-35 fail-if-appears pattern): assert
  `not IsBoundGlobal("<name>")` for each of `("IsQuasiTiltedAlgebra","IsShodAlgebra",
  "IsWeaklyShodAlgebra","IsLauraAlgebra","IsAdaAlgebra")` — a trip-wire that this scope note
  is stale if QPA ever adds one. Corroborate the gl.dim gate: `A.global_dimension().value`
  matches QPA's `GlobalDimension` on the fixture (`shod ⟹ gl.dim ≤ 3`, `quasi-tilted ⟹
  gl.dim ≤ 2` are QPA-checkable bounds). Reuse the existing P05/P23 QPA pd/id crosscheck verb
  for the pointwise parts data (grep `qpa/crosscheck.py` for the verb name before wiring; if
  none clean, the probe + the gl.dim tie is the deliverable — honest, never a silent skip).
- [ ] **Step 2–3:** run live (`-m qpa`), commit — `test(qpa): recognizer-ladder probe (QPA has no quasi-tilted/shod/laura/ada verb, fail-if-appears) + gl.dim gate + pointwise pd/id`

---

### Task F: verification page, README, metaplan, suite gate

**Files:** Modify `docs/verification.md`, `README.md`,
`docs/plans/2026-08-07-metaplan-v1.0.0.md` (§6 ledger). Test: `tests/release/
test_oracle_classes.py`, `tests/citations/`.

- [ ] **Step 1: verification page.** Add the Plan-61 rows:
  - `modules/recognizers_ladder.py` (rungs) — `oracle_literature` (**ACLV Example 2.2(b):
    (qt,shod,wshod,laura,ada)=(F,F,T,T,T), complement={S₃}, gl.dim 4**; hereditary ⇒ all five
    `True`; shod ⟺ complement empty; **weakly_shod=False on `kupisch=[3,2,2]`**); `oracle_selfcert`
    (monotone nesting + the guard-RAISES test; gl.dim theorem gates; strict-weakly-shod
    `HH^{≥2}=0`; Corollary 2.6 `gl.dim ≤ 4` / `pd≤2 or id≤1`; the non-ada witness); `oracle_crossengine`
    (quasi-tilted two-routes agreement; **shod two-routes** complement-empty ≡ QT2 sweep;
    weakly-shod AR-SCC ≡ Hom-closure).
  - `modules/recognizers_ladder.py` (ada/HH¹) — `oracle_literature` (**ACLV Theorem B**: the
    fixture over CC — ada, HH¹=0, simply connected; the gate-predicate on constructed algebras);
    `oracle_selfcert` (the pure `_theorem_b_verdict` False-branch when HH¹≠0; the gated-off
    non-CC path); `oracle_crossengine` (**P56 `is_simply_connected` agreement** on the ada/CC
    members; the Adian-Rabin-`None` resolution is the headline note). **NOTE — no rep-finite
    ada + HH¹≠0 pin:** ada strongly constrains rep-finite algebras toward simple connectedness
    (the monomial square 1→2→4/1→3→4 with `a·b=0` is NOT ada — `P₁` lies in the complement,
    verified live — so it cannot serve; a genuine rep-finite ada-not-simply-connected instance
    is a `# PIN` the worker searches for, with the False branch covered unconditionally by the
    synthetic `_theorem_b_verdict` unit test).
  - **Honest-scope entries:** (a) **rep-finite + non-self-injective only** (inherited from P55
    /P41; ACLV 2.2(c) rep-infinite ada is refused despite being mathematically ada); (b)
    **laura is trivially `True` in scope** — the complement is the reported datum, the
    interesting laura content is rep-infinite (out of scope); (c) **identification is QQ/char-0
    decisive** — large GF(p)/GF(pⁿ) may raise a loud whole-compute refusal (M1); (d) the
    **ada/HH¹ verdict is emitted only over an algebraically closed field (CC)** — over QQ/GF(p)
    the block reports HH¹ + the theorem note, no verdict (Theorem B's hypothesis unmet); (e)
    **`tilted ⟹ quasi-tilted`** is `# PIN`'d for **P60** (skipped auto-flipping test); (f)
    **QPA cannot compare** the recognizer surface (fail-if-appears probe; the pd/id/gl.dim
    ingredients ARE QPA-checked). Recount the class table (run
    `tests/release/test_oracle_classes.py`; paste LIVE counts for
    `oracle_literature`/`oracle_crossengine`/`oracle_selfcert`/`qpa` + union; mid-merge-train
    honest).
- [ ] **Step 2: README.** One features line: "the Assem-school recognizer ladder
  (quasi-tilted / shod / weakly-shod / laura / ada) with witnesses, and — for ada algebras
  over an algebraically closed field — HH¹ as a complete simple-connectedness oracle (ACLV
  Theorem B) — no-code in the browser."
- [ ] **Step 3: metaplan.** Tick `P61` in the §6 progress ledger.
- [ ] **Step 4: full gate:** `... -m pytest tests/modules/test_recognizer_ladder*.py tests/modules/test_ada_hh1_p56_agreement.py -q` (deep); `... -m pytest -q -m fast`;
  `... -m pytest tests/webapp tests/gui -q`; `... -m pytest tests/qpa -q -m qpa`;
  `... -m pytest tests/release tests/citations -q`; `; echo EXIT=$?` — all green.
- [ ] **Step 5: commit** — `docs(verification): Plan-61 recognizer-ladder oracle rows + HRS/CL/AC/Smith/BFT citations + honest scope (rep-finite, laura-trivial, QQ id, CC-only Thm-B verdict, tilted PIN, no-QPA) + recounted classes`

---

## Acceptance (Plan-61 definition of done)

1. `recognizer_ladder`, `RecognizerLadder`, `LadderRung`, `SimpleConnectedness`,
   `recognizer_ladder_block` public in `src/quiverlab/modules/recognizers_ladder.py`; the
   `Algebra` delegates (`recognizer_ladder` + `is_quasi_tilted`/`is_shod`/`is_weakly_shod`/
   `is_laura`/`is_ada`) named; every verdict certified/witnessed or loudly refusing.
2. The five rungs computed off **ONE `left_right_parts` atlas** + `global_dimension` (+ the
   AR-SCC weakly-shod sweep): **shod ⟺ `complement == []`** (survey Thm 4.1 a⟺b);
   **quasi-tilted ⟺ `gl.dim ≤ 2 ∧ complement == []`** cross-checked against ACLV
   **`every P_x ∈ 𝓛_A`**; **weakly shod ⟺** no non-trivial AR-SCC on an injective→projective
   path (WSA); **laura ⟹ `True`** (rep-finite) with the finite `complement` reported;
   **ada ⟺** every `P_x`,`I_x` placement `≠ "neither"` (ACLV Def 2.1).
3. **ACLV Example 2.2(b)** pinned: `(quasi_tilted, shod, weakly_shod, laura, ada) =
   (False, False, True, True, True)`, `complement = {S₃}` (pd 2, id 2), `gl.dim = 4` — the
   full five-way discrimination on one fixture (verified live). Hereditary ⇒ all five `True`,
   empty complement (`kA_n`).
4. The **monotone nesting** (`qt ⟹ shod ⟹ weakly_shod ⟹ laura`, `qt/shod ⟹ ada`,
   `ada ⟹ laura` in scope) and the **gl.dim theorem gates** (`shod ⟹ gl.dim ≤ 3`,
   `quasi-tilted ⟹ gl.dim ≤ 2`) are standing self-certs — a violation is a loud
   `QuiverlabError`, never a silently-wrong verdict.
5. The **ada/HH¹ block (ACLV Theorem B)**: the gate is **`A.domain.is_algebraically_closed`**
   (the new additive flag — `True` only on the CC working domain; there is no `Algebra.field`,
   and QQi shares CC's `SympyExactDomain` class). On `ada.verdict is True` over an
   algebraically closed field, `verdict = (dim HH¹(A) == 0)` — the complete simple-connectedness
   decision, with the "HH•(A) reduces to base field" corollary noted; over a non-alg-closed
   field (QQ/GF(p)/GF(pⁿ)/QQi), `dim HH¹` + the theorem statement, **no verdict**; not applicable
   off ada. Verified live: the gate predicate on constructed algebras (`field=CC` ⇒ `True`;
   `QQ`/`GF(7)`/`QQi` ⇒ `False`), and the tree fixture over CC (ada, HH¹=0, simply connected);
   the `verdict=False` branch (HH¹≠0) is covered by the synthetic `_theorem_b_verdict` unit
   test (a rep-finite ada-not-simply-connected instance is a `# PIN` — ada strongly pushes
   rep-finite algebras toward simple connectedness).
6. The **P56 agreement oracle**: on ada/CC members, P61's Theorem-B verdict agrees with
   `is_simply_connected(A).verdict` wherever P56 is definite, and resolves P56's
   Adian-Rabin `None` — behind an auto-flipping `importorskip` (the P56 alignment point).
7. The **honest semi-decision contract**: complete iff rep-finite and non-self-injective;
   self-injective ⇒ `status="unsupported"`, rep-infinite (ACLV 2.2(c)) ⇒ `status="budget"` —
   no partial ladder; identical to the P55/P41 loud cap. Identification QQ/char-0-scoped
   (M1).
8. `recognizer_ladder` clickable end-to-end (GUI canvas → ✓/✗/— ladder table with
   witness-per-✗ + laura complement + ada/HH¹ line → report), algebra-level compute kind
   (schema v1, NO module block), both runners byte-identical via the shared
   `recognizer_ladder_block`, EN+ES+FR+ZH i18n, one golden with a documented change-log entry.
9. QPA probe green (`-m qpa`): no quasi-tilted/shod/weakly-shod/laura/ada verb
   (fail-if-appears); the gl.dim gate + pointwise pd/id corroborated.
10. `docs/verification.md` recounted (live numbers, mid-merge-train honest); the five new
    citations added and BibTeX-verified (`aclv_supports` + `organising_module_category`
    reused from P55, not re-added); README line + metaplan P61 tick; deep + fast + webapp/gui
    + qpa + release + citations suites green. Honest scope recorded: rep-finite +
    non-self-injective; laura-trivial-in-scope; QQ identification; CC-only Theorem-B verdict;
    `tilted ⟹ quasi-tilted` `# PIN`'d for P60; QPA cannot compare the recognizers.

---

## Change-log

- **2026-08-07** — plan authored (P61 / R18). Definitions re-verified verbatim from the ACLV
  1102.1188 PDF (Def 2.1, Theorem B, Cor 2.6, Examples 2.2(a–c)) and the "Organising the
  module category" survey PDF (§4: HRS QT1∧QT2 quasi-tilted; Thm 4.1 shod ⟺ ind A = 𝓛∪𝓡;
  WSA weakly shod; §5.1 strict-weakly-shod `HH^{≥2}=0`). **Verified live in the venv:** the
  rad²=0 A₅ fixture ladder `(F,F,T,T,T)` with `complement={S₃}`, gl.dim 4, `dim HH¹ = 0` over
  QQ and CC (simply connected, Theorem B); hereditary `kA_n` all-five-`True`; the comm-square
  `a·b−c·d` all-five-`True`. P55's `left_right_parts` is not yet on `dev` (unchecked in the
  metaplan ledger) — this plan designs against the committed P55 atlas contract and branches
  after P55 merges.
- **2026-08-07 adversarial review round 2 (blocking + 4 majors + 3 minors applied).**
  **Blocking H1 — the alg-closed gate was WRONG:** the authoring draft gated on
  `isinstance(A.field, ComplexField)`, but (verified live) `Algebra` has **no `.field`**
  (it stores `self.domain`), and `field=CC` yields a `SympyExactDomain` — the `ComplexField`
  sentinel is discarded — AND `field=QQi` (ℚ(i), not algebraically closed) yields the SAME
  `SympyExactDomain` class. The gate is respecified as a new additive
  **`Domain.is_algebraically_closed`** flag (`True` only on the CC working domain, verified
  live: `CC ⇒ True`, `QQ`/`GF(7)`/`QQi ⇒ False`), gated on `A.domain.is_algebraically_closed`;
  the earlier "verified live: `isinstance(CC, ComplexField)`" claim was a **mis-statement of
  what was tested** (the sentinel, not a constructed algebra) — corrected. **Major M1+W3+C1:**
  the shod two-routes oracle made REAL via the P55 atlas `pd_le_1`/`id_le_1` contract amendment
  (addendum appended to P55's doc, committed together — zero new math); the AR-SCC ≡ Hom-closure
  weakly-shod oracle implemented as a real two-route test; the placement-`None` contract stated
  (C2). **Major W1:** the monomial square is **NOT ada** (derived live — `P₁` in the complement),
  so it cannot be the ada-not-SC oracle; the vacuous `if L.verdict("ada")` guard is removed —
  the not-SC verdict-False branch is now covered by an unconditional synthetic
  `_theorem_b_verdict` unit test, `test_not_ada_no_verdict` uses a real non-ada witness (the
  triangle `1→2→3`+`1→3`, `a·b=0`, verified `ada=False`), and the quasi-tilted/shod two-routes
  tests pin the comm-square `True` verdict. **Major M2:** the weakly-shod seam exercised on the
  real non-directed `NakayamaAlgebra(kupisch=[3,2,2])` (knits complete, one big SCC on an I⇝P
  route, `weakly_shod=False`, verified live). **Minors:** the monotone-guard-RAISES test added
  (M3); the P56 agreement test parametrized over two ada/CC members (W2); the block builder
  passes the HH `verbose/trace` sink so the ladder writes no stray trace files, and the
  disconnected-input note added (§5). Open `# PIN`s: the shod-origin Manuscripta Math. 100 (1999)
  secondary citation; a genuine rep-finite ada-not-simply-connected algebra for the Theorem-B
  False branch (elusive — ada strongly pushes rep-finite algebras toward simple connectedness;
  synthetic unit test covers the branch meanwhile); the `tilted ⟹ quasi-tilted` auto-flipping
  test for P60.
