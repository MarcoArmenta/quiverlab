# Plan 69: Persistence / TDA bridge — barcodes over A_n, zigzag modules, and commutative ladders (R33)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking. Do the reference re-verification in
> the Record section BEFORE writing any oracle pin. **The CL(2)=11 / CL(3)=29 vertex counts
> are single-engine `oracle_selfcert` values, NOT literature pins — they rest on one engine
> (the AR knitter) until reconciled** (H3). Reconciling them against the Escolar–Hiraoka
> AR-quiver figures is a **BLOCKING task** (Task 6 Step 1a): count vertices from the
> RENDERED figure, never from `pdftotext` — the critic confirmed `pdftotext` linearizes the
> 2-D AR figures unreliably (reads ~30, ±1–2). **Two hard build gates that broke the
> author's own draft (verify first):** (1) build `CommutativeLadder` with **scalar** vertex
> names only — tuple vertices `(i,j)` crash `knit_ar_quiver` at `complex_reps.py:179`
> (H1, live-reproduced); (2) the `barcode` estimator tier is NOT "no change" — a CL barcode
> knits and must be pushed off the instant tier (H2).

**Goal.** Put topological data analysis (TDA) on quiverlab's representation-theory engine
— **representation theory FIRST, persistence as its dictionary** (Marco's positioning):

- **Barcodes = interval decompositions.** A *persistence module* over `A_n` (a type-`A`
  quiver of any orientation) is nothing but a finite-dimensional representation of `A_n`;
  by **Gabriel's theorem** its indecomposables are the **interval (thin) modules**, and its
  **barcode** is the multiset of the support intervals of `decompose(M)`'s summands
  (Krull–Schmidt, Plan 30). One orientation direction (`1→2→…→n`) is *ordinary*
  persistence; an *alternating* orientation is a **zigzag module** (Botnan–Crawley-Boevey).
  This path is **field-robust over every exact domain including GF(2)** — the TDA-native
  field — because interval modules are **bricks** (`dim End = 1`), so `decompose` certifies
  them in *any* characteristic (LIVE-VERIFIED, §Field honesty).

- **Commutative-ladder persistence.** `CL(n) = A_n □ A_2` (the box product, a bound quiver
  with **commuting squares**) is **representation-finite iff `n ≤ 4`** (Escolar–Hiraoka;
  "length < 5"). For `n ≤ 4` we deliver the **AR-quiver-indexed generalized persistence
  diagram** — `decompose(M)` + each summand matched to its position in the knitted AR
  quiver (Plan 41). For `n ≥ 5` the ladder is **representation-infinite**: a **loud typed
  refusal**, never a silent partial diagram.

- **`barcode` GUI kind.** A module-side no-code compute kind (it consumes the Plan-26
  module input — a persistence module IS such a module), rendered as an interval table +
  a simple HTML bar diagram (no new canvas), in all three tiers and four locales.

**Architecture.** One new module + one new family + the module-kind wiring; every number
is exact and either self-certifies or refuses loudly. No new math engine — a thin exact
layer over primitives already on `dev`:

- `src/quiverlab/families/commutative_ladder.py` — `CommutativeLadder(n, base_orientation
  ="forward", field=None) -> Algebra` (the presented `kQ/I` box product `A_n □ A_2` with
  commutativity relations, refusing `n ≥ 5` loudly), `is_commutative_ladder(A) -> (bool,
  n)`, and a `persistence_line(n, orientation="forward"|"zigzag"|dict) -> Quiver`
  convenience for the `A_n` base. Consumes `combinat/quiver.py::Quiver` +
  `Quiver.algebra` (the `IncidenceAlgebra` non-monomial commutativity idiom `"a*b - c*d"`).
- `src/quiverlab/modules/barcode.py` — `barcode(M, ...) -> Barcode` and `barcode_block(A,
  M)` (the compute kind). Consumes `modules/decompose.py::decompose` (the arbiter),
  `modules/hom.py::is_isomorphic`/`identify_standard`, `modules/ar.py::knit_ar_quiver`
  (CL AR-indexing only), `Module.dimension_vector()`.
- The `barcode` **module-side** compute kind wired into `hpc/spec.py` (`MODULE_KINDS`,
  `_MOD_REFS`, `_dispatch_module`), the Pyodide twin `docs/gui/runner.py` (`_MODULE_KINDS`,
  `_MOD_REFS`, `_module_block`), the pydantic schema `webapp/server/schema.py`
  (`MODULE_KINDS` — auto-guards via `_module_rules`), `trace/results_html.py` (`_HEADINGS`
  + `_block_html`), the two byte-identical `gui.js`, i18n ×4, one golden.

**Tech Stack.** Pure exact linear algebra over `Domain` (`modules/linalg_mod`,
`fields/linalg`) + exact combinatorics on the quiver; the AR knit and decompose are the
shipped engines. **No floats in `src/`** (AST-gated by `tests/test_no_floats.py`): barcode
endpoints are `int` (filtration = discrete vertex index); an optional per-index
`filtration_values` label is an exact `int`/`Fraction` (never a float — real-valued
filtrations are handled as exact rationals or refused, §Scope).

---

## Records (verbatim — the plan's charter)

> **R33 — Persistence modules over A_n and commutative ladders (TDA bridge).** [D-scout
> P9; keep] Object: barcodes = interval decompositions of A_n/zigzag modules (existing
> string/AR machinery); CL(n) = A_n □ A_2 for n ≤ 4 (rep-finite — "length < 5" verified
> exact): AR-quiver-indexed generalized persistence diagrams; CL(≥5) loud rep-infinite
> refusal. Refs: Escolar–Hiraoka arXiv:1404.7588 (DCG; Hiraoka-group provenance honestly
> flagged); Igusa–Rock–Todorov arXiv:1909.10499; Botnan–Crawley-Boevey. Oracles: A₅
> filtration barcode; Escolar–Hiraoka's CL(3) worked diagram. Size S–M.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P69 (Wave 3, tier γ,
independent). Tier placement in the research doc §9: "R33 (persistence/TDA — outreach
value)."

---

## Reference re-verification (mandatory — done at authoring, 2026-08-07)

Web-verified against the primary sources; findings recorded so the implementer does not
re-derive, and so the definitions are pinned against transcription error. **The
representation-finiteness THEOREM (`n ≤ 4`) is now literature-confirmed verbatim (see below);
the two figure-only Escolar–Hiraoka vertex COUNTS remain `# BLOCK` (single-engine
`oracle_selfcert` until reconciled); the durable oracles are the theorem + the engine's own
knit/decompose (self-cert).**

1. **Escolar–Hiraoka, arXiv:1404.7588** = *Persistence Modules on Commutative Ladders of
   Finite Type*, **Discrete & Computational Geometry 55(1) (2016) 100–157**, DOI
   `10.1007/s00454-015-9746-2` (abstract + search-verified journal metadata). Verified content:
   - **Theorem (headline) — LITERATURE-CONFIRMED VERBATIM (W4).** The critic extracted the
     PDF with `pdftotext`, which yields the theorem prose reliably (it is a text line, not a
     figure): *"representation-finite for `n ≤ 4` … representation-infinite for `n ≥ 5`"*,
     for **arbitrary orientation** (matching the abstract's *"commutative ladders of length
     less than 5 are representation-finite and their Auslander–Reiten quivers are explicitly
     shown"*). So `CL(n)` is representation-finite **iff `n ≤ 4`**. This is the **binding
     literature pin** and the whole scope boundary — no longer "could not verify".
   - **The commutative ladder** `CL_n(τ)`: the base `A_n` carries an orientation `τ`
     (a length-`n` string of `f`/`b`), the "ladder" direction `A_2` (the rungs) is fixed;
     the underlying **bound quiver is `A_n □ A_2` with every square commutative**. The
     representation TYPE (finite ⇔ `n ≤ 4`) does **not** depend on `τ`; the AR quiver and
     the **number of indecomposables DO** depend on `τ` (so a "CL(3) diagram" is
     orientation-specific — name the orientation in every pin).
   - **Not all indecomposables are intervals.** For `CL(n)` (unlike `A_n`) the AR quiver
     contains **non-interval** indecomposables, so the "generalized persistence diagram"
     is indexed by the **AR quiver**, not by interval endpoints (Escolar–Hiraoka's own
     framing). The interval summands are a distinguished sub-multiset.
   - **`# BLOCK` — the exact indecomposable counts / AR-quiver-figure vertex counts per
     orientation (single-engine `oracle_selfcert`, NOT `oracle_literature`).** The paper
     gives the AR quivers as **figures**; the vertex counts are not reliably machine-extractable.
     **The critic's finding (H3): `pdftotext` linearizes the 2-D AR figures unreliably — it
     read ~30 for the CL(3) figure, ±1–2 — so the worker MUST count vertices from the
     RENDERED figure (an image of Fig. 14), never from text extraction.** **Do NOT pin a
     paper-figure count from this plan's prose.** The counts `CL(2)=11` / `CL(3)=29` are
     asserted **against the same AR knitter that produces them** — they rest on ONE engine
     until reconciled, so they are `oracle_selfcert`, not `oracle_literature`. Reconciling
     against the paper figure is a **BLOCKING task** (Task 6 Step 1a); a cheaper second route
     is the **QPA cross-engine closer** (Task 4 — `DecomposeModuleWithMultiplicities` /
     rep-type on CL(3); see the QPA note there). **LIVE-VERIFIED on `dev` (2026-08-07 and
     re-verified 2026-08-08, `QUIVERLAB_NO_NUMBA=1`, QQ):** the fully-forward orientation
     `CL(n)`, built as the incidence algebra of the `[n]×[2]` grid poset **with SCALAR
     vertex names** (see the H1 tuple-vertex crash note in Task 1 — the natural tuple labels
     `(i,j)` crash the knit), knits **complete** (rep-finite) with `CL(2) → 11`,
     `CL(3) → 29` indecomposables (`dim CL(2) = 9`, `dim CL(3) = 18`). These are the shipped
     self-cert pins.
   - **`# BLOCK` — the CL(3) worked persistence-module example** (a specific block-matrix
     / persistence module and its indecomposable decomposition). The metaplan names "the
     CL(3) worked diagram"; the paper's specific example is figure/algorithm-level and was
     not extractable. Ship instead a **constructed** CL(3) module decomposed by the engine
     (self-cert: the summand dim-vectors sum to `M`), and mark the paper's own example as a
     `# PIN` to add at build if the worker can transcribe it from the PDF.

2. **Botnan–Crawley-Boevey**, *Decomposition of persistence modules*, **Proc. AMS 148(11)
   (2020) 4581–4596**, arXiv:1811.08946 (title/venue web-verified). The **structure
   theorem**: a pointwise-finite-dimensional persistence module (over a totally ordered or
   zigzag poset) decomposes into **interval modules**, uniquely (Krull–Remak–Schmidt–Azumaya).
   This is the theoretical justification that "barcode = interval decomposition" is
   well-defined for `A_n` and zigzag. (For a finite `A_n` it is Gabriel; Botnan–CB is the
   general PID-free statement, cited for the zigzag direction.)

3. **Igusa–Rock–Todorov**, arXiv:1909.10499 = *Continuous Quivers of Type A (I)* (title
   web-verified; a series continued in 1910.\*/2004.\*). The continuous-limit
   representation theory of type-`A` persistence — cited as the **conceptual bridge**
   (representation theory ⇄ persistence), not as a computed oracle (quiverlab is finite/exact).

4. **Gabriel's theorem** (the `A_n` indecomposables = positive roots = intervals, all
   **thin**, `End = k`). Classical (Gabriel, *Manuscripta Math.* 6 (1972) 71–103). Treated
   as ground truth; it is why the `A_n`/zigzag barcode is field-robust (bricks).

5. **`Asashiba–Escolar–Hiraoka–Takeuchi`**, *Matrix Method for Persistence Modules on
   Commutative Ladders of Finite Type*, **Japan J. Indust. Appl. Math. 36(1) (2019)
   97–130**, arXiv:1706.10027 (web-verified) — the block-matrix normal-form algorithm for
   `CL(n≤4)`. Optional supporting citation for the CL diagram; the AR-knit route is our
   primary engine, so this is documentation, not a computed oracle.

`# PIN` items (freeze BibTeX metadata at build, never guessed here): the full author lists
/ volume-page of **Escolar–Hiraoka 2016**, **Botnan–Crawley-Boevey 2020**,
**Igusa–Rock–Todorov 2019**, and (optional) **AEHT 2019** — each added to `references.bib`
only after the worker BibTeX-verifies it at merge (the P49/P59 precedent). The
paper-figure CL(3) counts and worked example are `# BLOCK` until reconciled with the PDF.

---

## TDA glossary — representation theory FIRST (Marco's positioning)

Every persistence object is a representation-theoretic object we already compute; the
report and GUI lead with the rep-theory name and gloss the TDA term.

| TDA term | Representation theory (what quiverlab computes) |
|---|---|
| filtration `X_1 ⊆ X_2 ⊆ … ⊆ X_n` | (context) the topological input; homology `H_k(X_•)` produces the module |
| **persistence module** | a finite-dimensional **representation `M` of `A_n`** (all arrows `1→2→…→n`) |
| **zigzag module** | a representation of `A_n` with an **alternating orientation** (e.g. `1→2←3→4`) |
| filtration **parameter / time** | the **vertex index** `i ∈ {1,…,n}` (discrete; the poset order) |
| **birth** `b` of a bar | the **smallest vertex** in the support interval of an interval summand |
| **death** `d` of a bar | the **largest vertex** in the support interval. **`essential` (alive at the top, `d = n`) is a FORWARD-persistence-only flag** (M2): for `kind="persistence"`, `essential = (d == n)`; for `kind="zigzag"` there is no single monotone "top", so *"alive at the top"* has no meaning — `essential` is **always `False`** on a zigzag bar (the support's `min`/`max` in the line order are still reported; a finite zigzag has no infinite/essential bars) |
| **bar `[b,d]`** | an **interval module `I[b,d]`** (thin: `dim = 1` on `b..d`, `0` elsewhere) |
| **barcode** | the **multiset of interval summands** of `decompose(M)` (Gabriel / Botnan–CB) |
| **persistence diagram** | the same data as `(birth, death)` points (barcode ⇄ diagram bijection) |
| **generalized persistence diagram** (CL) | `decompose(M)` **indexed by the AR quiver** of `CL(n≤4)` (Escolar–Hiraoka) |
| **elder rule** | Krull–Schmidt uniqueness: the older class survives a merge (the interval decomposition is canonical) |

The `barcode` report block states, in one line: *"A persistence module over `A_n` is a
representation of the quiver `A_n`; its barcode is the interval decomposition of `M`
(Gabriel/Botnan–Crawley-Boevey), i.e. the support intervals of the Krull–Schmidt
indecomposable summands."*

---

## Mathematical spec (the plan's ground truth)

### The persistence-line algebra and its indecomposables

Let `A_n` be a type-`A` quiver on the path graph `1 — 2 — … — n` with **any** orientation
of the edges (each edge `{i,i+1}` an arrow `(i,i+1)` or `(i+1,i)`). The **line order** is
the unique path order of the underlying graph (its two degree-1 endpoints; order from the
lower-labelled endpoint — deterministic and orientation-independent). `A = kA_n` (path
algebra, `relations=[]`) is hereditary and rep-finite; by **Gabriel** its indecomposables
are exactly the **interval modules** `I[b,d]` (`1 ≤ b ≤ d ≤ n`), each **thin** (`dim vector`
= the 0/1 indicator of `{b,…,d}` in the line order) and a **brick** (`End(I[b,d]) = k`).

**The `n = 1` edge case (M1).** A single-vertex quiver (`A_1 = k`, one vertex, no arrows) is
degenerate for the general line classifier — it has **zero** degree-1 endpoints, not two —
so it is handled as an **explicit special case** rather than refused: `A_1` is the line of
length 1, its `line_order` is that one vertex, and `barcode(M)` returns a single bar
`[1,1]` with `multiplicity = dim M` at that vertex (`kind = "persistence"`,
`essential = True` since `death == n_top == 1`). This is the honest choice over a loud
refusal (a length-1 filtration is a legitimate, if trivial, persistence module); a test
pins it (Task 2).

**Barcode.** For a representation `M` of `A_n` (any orientation), `decompose(M) =
[(I[b_i,d_i], m_i)]` (all summands are intervals, by Gabriel). The **barcode** is the
multiset `{ [b_i, d_i]^{m_i} }`, read off each summand's dimension vector (its support is a
contiguous interval in the line order — **asserted**; a non-contiguous support would mean
the input quiver is not `A_n`-shaped → loud error). Forward orientation ⇒ ordinary
persistence; alternating ⇒ a **zigzag barcode** (Botnan–Crawley-Boevey). For **forward**
`kind="persistence"` the `d_i = n` bars are flagged **essential** (alive at the top of the
filtration) — never printed as `∞` (no `∞` the finite computation did not produce; the
honest finite-filtration convention). **For `kind="zigzag"` the `essential` flag is not
applicable** (M2): a zigzag has no single monotone "top" direction, so *alive at the top*
carries no persistence meaning — every zigzag bar has `essential = False` (a finite zigzag
line has only finite bars; the endpoints `min`/`max` in the line order are still reported).
The convention is stated in the glossary, the `Bar` schema, and the report block.

**Self-cert identity.** `Σ_i m_i · dimvec(I[b_i,d_i]) == M.dimension_vector()` exactly.

### The rank-formula cross-check (a genuinely independent second route, forward only)

For **forward** `A_n` (`1→2→…→n`), the multiplicity of `[b,d]` has the classical closed
form via composite ranks `r_i^j = rank(M_i → M_j)` (`r_i^i = dim M_i`, out-of-range `= 0`):

    mult[b,d] = (r_b^d − r_{b−1}^d) − (r_b^{d+1} − r_{b−1}^{d+1}).

This uses **only `mat_rank` of composite arrow-action matrices** — no `decompose`, no
`Hom` — so its agreement with the `decompose` barcode is a real `oracle_crossengine`
(two independent implementations). Ship `_barcode_by_ranks(M)` (forward-only; refuses on a
non-forward/zigzag line, where the naive rank formula does not apply — zigzag uses the
`decompose` route + Botnan–CB) and assert equality on the forward oracles.

### Commutative ladders `CL(n) = A_n □ A_2`

The box product of the `A_n` base (orientation `τ`) and `A_2` (the rung `1→2`): vertices
`(i,j)` for `i ∈ {1,…,n}`, `j ∈ {1,2}`; arrows = the base arrows in each row (`(i,1)→(i+1,1)`
and `(i,2)→(i+1,2)` for forward `τ`) plus the rungs `(i,1)→(i,2)`; relations = **every
square commutes**, `a_i · c_{i+1} − c_i · b_i` (left-to-right composition, the
`IncidenceAlgebra` `"a*b - c*d"` idiom). For fully-forward `τ`, `CL(n)` is exactly the
**incidence algebra of the `[n]×[2]` grid poset** (`dim = 3·binom(n+1,2)`; `dim CL(3) = 18`,
LIVE-VERIFIED). **The grid poset must be built with SCALAR vertex names** (e.g. `"i_j"` or a
flat integer index), **never the natural tuple labels `(i,j)`** — tuple vertices crash the
downstream AR knit at `complex_reps.py:179` (H1; the full crash note is in Task 1). The
"LIVE-VERIFIED `dim 18`" run used the scalar-relabelled grid (W1).

**Escolar–Hiraoka:** `CL(n)` is **representation-finite iff `n ≤ 4`**. For `n ≤ 4` the AR
quiver is finite; the **generalized persistence diagram** of a module `M` is `decompose(M)`
with each indecomposable summand **matched to its AR-quiver vertex** (via a dimension-vector
prefilter + exact `is_isomorphic`, the `identify_standard` idiom generalized over *all*
knitted indecomposables). Interval summands (thin support) are flagged; non-interval
summands carry their AR index and dimension vector. For `n ≥ 5`: **loud rep-infinite
refusal** (the AR quiver is infinite; no complete diagram exists).

**Char scope for the CL route only.** The AR knit (`knit_ar_quiver`) and `decompose` of a
**non-brick** indecomposable both rest on the trace-form radical, rigorous over **char 0 or
char > dim** (`ar.py::_rad_end_basis`, `decompose.py::_certify_local`). CL indecomposables
need not be bricks, so **CL barcodes run over `QQ` (or `char > dim`)**; off scope,
`knit`/`decompose` refuse loudly and `barcode` surfaces that refusal. (Contrast the
`A_n`/zigzag route — all summands are bricks — which is field-robust; §Field honesty.)

**The char-scope refusal chain — verify the knit surfaces `"error"`, never a silently-drained
`"complete"` (M6).** Over `GF(2 ≤ dim)` the CL route refuses through a specific chain the
worker must confirm end-to-end: `knit_ar_quiver` internally calls `decompose` on a
**non-brick** indecomposable → `decompose.py::_certify_local` raises `QuiverlabError` (the
`r <= 1` brick short-circuit at `decompose.py:316` does NOT fire, so control reaches the char
gate at `decompose.py:321-322` and raises) → **`knit_ar_quiver` CATCHES the error and returns
`ARQuiver` with `status="error"`, `is_complete=False`** (a partially-knit quiver, NOT a raise)
→ `barcode`'s CL branch checks `not ar.is_complete` and **refuses loudly**. So the refusal is
`barcode`'s, keyed off `is_complete`, not an uncaught knit exception. **LIVE-VERIFIED
(2026-08-08, `QUIVERLAB_NO_NUMBA=1`):** the scalar-relabelled `CL(3)` over `GF(2)` knits to
`status="error"`, `is_complete=False`, 8 partial vertices (NOT `"complete"`/29). The worker
asserts `status == "error"` / `not is_complete` at the knit boundary in addition to the
`barcode` raise, so a future silent-drain regression (a knit that swallows the error and
reports `"complete"` with a truncated vertex set) is caught.

### The `Barcode` result object

```python
@dataclass(frozen=True)
class Barcode:
    kind: str                 # "persistence" (forward A_n) | "zigzag" | "commutative_ladder"
    line_order: tuple         # the A_n vertex order used to read intervals (or the CL base)
    bars: tuple               # tuple of Bar(birth, death, multiplicity, essential, dimvec)
                              #   -- for A_n/zigzag: interval summands
    diagram: tuple | None     # CL only: tuple of (ar_index|None, ar_name|None, dimvec,
                              #   multiplicity, is_interval); None for A_n/zigzag
    n: int                    # the ladder/line length
    field: str
    field_robust: bool        # True for A_n/zigzag (bricks, any char); False for CL (char-scoped)
    ar_status: str | None     # CL only: knit status ("complete" for n<=4); None otherwise
    references: tuple
```

`Bar` is `(birth:int, death:int, multiplicity:int, essential:bool, dimvec:dict)`. **The
`essential` field is FORWARD-persistence-only (M2):** it is `(death == n)` only for
`kind="persistence"`; for `kind="zigzag"` it is **always `False`** (a zigzag has no monotone
"top", so "alive at the top" has no meaning — see the glossary). For CL interval bars it is
likewise `(death == n)` on the base line only when the CL base orientation is forward,
`False` otherwise; the CL primary object is the `diagram`, not the bars.

### Public API surface (`import quiverlab`; exact only, no floats)

```python
# families/commutative_ladder.py
def CommutativeLadder(n, base_orientation="forward", field=None) -> Algebra
    # A_n [] A_2 presented kQ/I with commuting squares. base_orientation in
    # {"forward","backward", <dict edge->('f'|'b')>}. Refuses n >= 5 LOUDLY
    # (Escolar-Hiraoka: representation-infinite) and n < 1.
def is_commutative_ladder(A) -> tuple            # (bool, n) -- graph-shape recognizer
def persistence_line(n, orientation="forward") -> Quiver
    # the A_n base: "forward" (1->2->..->n), "zigzag" (1->2<-3->4<-..), or a dict.

# modules/barcode.py
def barcode(M, *, budget=512, budget_modules=256) -> Barcode
    # M a representation of an A_n line (any orientation) OR a CL(n<=4). A_n/zigzag:
    # interval decomposition (field-robust). CL(n<=4): AR-indexed generalized diagram
    # (char 0 / char > dim). n==1 (single-vertex A_1) is HANDLED as a special case, not
    # refused: one bar [1,1] with multiplicity = dim M (M1). Refuses loudly: non-A_n/non-CL
    # quiver, CL(n>=5), presentation-less algebra, char-undecidable decompose (propagated).
def barcode_block(A, M) -> dict                  # the `barcode` compute kind (Task 5)
```

`Algebra.barcode(M, ...)` thin delegate (lazy-import), beside `ar_quiver`.

---

## Field / exactness honesty (the load-bearing reconciliation)

**TDA practice is GF(2) (or a small GF(p)); our exact domains cover it — and the
`A_n`/zigzag barcode is certified over EVERY field, GF(2) included.** The reconciliation
with `decompose`'s char-`p` rigor window (`char 0` or `char > dim M`, from the trace-form
radical) is exact and was **LIVE-VERIFIED**:

- Every `A_n` (and zigzag `A_n`) indecomposable is a **thin interval module**, hence a
  **brick** (`dim End = 1`). `decompose`'s `_certify_local` (`decompose.py:311`) certifies
  a summand indecomposable via the **`dim End = 1 ⇒ End = k·id ⇒ local` branch, valid in
  ANY characteristic** — it never reaches the char-scoped trace-form fallback. The
  **splitting** step (`_try_split`) factors a candidate endomorphism's minimal polynomial,
  and factoring is **supported over `GF(p)`, `QQ`, and CC number fields**
  (`decompose.py::_factoring_supported`). So the whole interval decomposition is exact and
  certified over `GF(2)`/`GF(p)`/`QQ` alike.
- **LIVE EVIDENCE (2026-08-07, `dev`, `QUIVERLAB_NO_NUMBA=1`):** on `kA_5` (`1→2→3→4→5`)
  the module `P_1 ⊕ S_3 ⊕ I_5 ⊕ P_2` decomposes to the **identical** summand
  dimension-vector multiset over `QQ`, `GF(2)`, and `GF(3)`; the hand-computed filtration
  barcode (below) is byte-identical over `QQ` and `GF(2)`; a **zigzag** `1→2←3→4←5` module
  decomposes identically over `QQ` and `GF(2)`. So `Barcode.field_robust = True` for the
  `A_n`/zigzag route, and the batteries run a `GF(2)` cross-check that MUST match `QQ`.
- **The repeated-interval hedge was overstated (critic-verified).** The Task-2 "if
  `decompose` refuses on a pathological small-field repeated-interval budget trip" caveat is
  a **defensive contract, not an observed failure**: the critic verified that **repeated
  interval modules over `GF(2)` split fine up to multiplicity `m = 4`** (bricks split by the
  `dim End = 1` short-circuit regardless of multiplicity — the plan's own risk pointed
  against its own field-robustness claim and does not materialize). Keep the loud-on-refusal
  contract, but do not present it as a likely `A_n`/zigzag outcome.
- **The CL route is honestly char-scoped.** CL indecomposables need not be bricks, so a CL
  barcode certifies only over `char 0` or `char > dim`; `barcode(CL_module over GF(2 ≤ dim))`
  surfaces the `decompose`/`knit` `QuiverlabError` (never a guess). `field_robust = False`
  there; the report states the scope. **The refusal is via the knit's `is_complete=False` /
  `status="error"`, not an uncaught raise** (M6; live-verified).

**Real-number filtration values are OUT OF SCOPE as floats (exact-only house rule).** The
persistence PARAMETER is the **discrete vertex index** — the barcode combinatorics depend
only on the *order type* of the filtration, which the vertex order captures exactly. A user
who wants labelled births/deaths may attach an optional exact `filtration_values` list
(one `int`/`Fraction` per vertex, monotone; **display-only**, echoed on the bars); a float
label is **refused** (`QuiverlabError`, the `_valid_entry` precedent). No approximate
persistence, no float thresholds — by design.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **Prerequisites present on `dev` (verify before branching).** P30 (`modules/decompose.py
  ::decompose`/`is_indecomposable`, char guard at `decompose.py:319-329`), P41
  (`modules/ar.py::knit_ar_quiver`/`ARQuiver`, `_std_name`, char guard `_rad_end_basis`),
  P37 (`modules/hom.py::is_isomorphic`, `modules/morphism.py::direct_sum`), P23/24
  (`modules/hom.py::identify_standard`), P26 (the no-code module schema + both runners),
  the `IncidenceAlgebra` commutativity idiom (`families/incidence.py:32`). Verify:
  `python -c "import quiverlab.modules.decompose, quiverlab.modules.ar,
  quiverlab.modules.hom, quiverlab.families.incidence"`. **No dependency on other Wave-3
  plans (P63–P68)** — P69 merges independently.
- **Test buckets are auto-assigned by directory** (`tests/conftest.py`): `tests/modules/`
  and `tests/families/` → **deep**; `tests/webapp/`, `tests/gui/`, `tests/hpc/` → **fast**;
  `tests/qpa/` → **qpa**. Barcode batteries live in `tests/modules/test_barcode_p69.py`
  (deep — they decompose/knit); the ladder builder in
  `tests/families/test_commutative_ladder_p69.py` (deep); GUI/runner cross-tests in
  `tests/webapp/` + `tests/gui/` (fast); QPA in `tests/qpa/` (qpa). Run by path during
  development; finish each task with a `-m deep`/`-m fast`/`-m qpa` spot-run of touched
  files.
- **Loud typed refusals** (`GlobalDimension`/`enumerate_strings` precedent, all
  `QuiverlabError` with a `hint=`): `barcode` refuses on a non-`A_n`/non-CL quiver, a
  `CL(n ≥ 5)`, a presentation-less (structure-constant) algebra, or a char-undecidable
  `decompose`/`knit`; `CommutativeLadder(n ≥ 5)` refuses (rep-infinite). Never a silent
  partial barcode.
- **No floats in `src/`** (AST-gated). Bar endpoints/multiplicities are `int`; dim vectors
  `dict`; `filtration_values` `int`/`Fraction`; `essential` a `bool`. `∞` is NEVER used —
  a top-reaching FORWARD bar is `death = n` with `essential = True`; `essential` is a
  forward-only flag, always `False` for zigzag (M2).
- **Composition is left-to-right** (`a*b` = first `a` then `b`; the CL commutativity
  relations use this). Vertex/line order is deterministic (Quiver vertex order).
- **Convention/scope choices are ARBITRATED or GATED, not assumed** (house
  cup-sign/composition-order precedent): the barcode is arbitrated by TWO routes on the
  forward oracles (`decompose` ≡ rank-formula); the `field_robust` claim is gated by the
  `GF(2) ≡ QQ` battery; the CL char scope is hard-gated (loud off `char 0`/`char > dim`);
  the CL rep-finiteness boundary is the Escolar–Hiraoka theorem (`n ≤ 4`), a hard gate.
- **Plan-32 markers** (orthogonal): `oracle_selfcert` = the interval-sum identity (`Σ m_i·
  dimvec == M.dimvec`), each summand's support is a genuine interval, the CL AR-index
  match (`is_isomorphic` certificate), the essential-bar flag, loud refusals, **and the
  single-engine `CL(2)=11`/`CL(3)=29` vertex counts** (asserted against the same AR knitter
  that produces them — H3; they rest on ONE engine until reconciled against the EH figure /
  QPA);
  `oracle_crossengine` = `decompose` ≡ rank-formula (forward), `GF(2) ≡ QQ` field-robustness,
  barcode ≡ `decompose`, CL diagram ≡ `knit_ar_quiver` counts, **QPA
  `DecomposeModuleWithMultiplicities` parity on `A_n`** (and, if buildable at merge, the QPA
  CL(3) closer for the 29 count — Task 4);
  `oracle_literature` = the `A_5` filtration barcode `{[1,5],[2,2],[4,4]}`, and the
  Escolar–Hiraoka **rep-finiteness BOOLEAN** `CL(n)` rep-finite iff `n ≤ 4` (the `n ≤ 4`
  knit terminates complete AND `CL(5)` is refused — the theorem, literature-confirmed
  verbatim; **NOT the vertex counts**, which are self-cert above); QPA lives in
  `tests/qpa/` (bucket = the class, never double-marked). **Oracle-marked tests live in
  `tests/modules`/`tests/families`** — NEVER in `tests/webapp`/`gui` (the release gate
  forbids markers there).
- **Mid-merge-train counts.** v1.0.0 lands many subplans in overlapping waves; absolute
  suite counts drift. **Task 6 recounts the oracle-class table at merge time** by running
  `tests/release/test_oracle_classes.py` (paste the LIVE numbers, never a guessed count)
  and claims only the deltas this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table
  green, `tests/release/test_oracle_classes.py` green) and adds its citations to
  `citations/references.bib` + `registry.py` (`bibtex()` hard-fails if the two drift).
  Conventional commits; green at every commit; branch `plan-69-persistence-tda` off `dev`
  (do not commit/push until asked).

---

### Task 1: `CommutativeLadder` + `persistence_line` + recognizer (families)

**Files:**
- Create: `src/quiverlab/families/commutative_ladder.py`
- Modify (exports): `src/quiverlab/families/__init__.py` (export `CommutativeLadder`,
  `is_commutative_ladder`, `persistence_line`) **and** `src/quiverlab/__init__.py`
  (top-level export + `__all__`, beside the `BrauerGraph` line — the webapp discovery reads
  the top-level export).
- Modify (CATALOG): `src/quiverlab/families/discover.py` — add a `FamilyInfo` entry.
- Modify (BOTH discovery skip-sets — `CommutativeLadder` takes a scalar `n` and IS
  form-buildable, so it may NOT need skipping; but `persistence_line` returns a `Quiver`,
  not an `Algebra` — keep it OUT of the family catalog. Verify `discover.py` introspection
  accepts `CommutativeLadder(n, base_orientation="forward")` as a scalar form; if the
  `base_orientation` string arg trips the form-builder, add `"CommutativeLadder"` to the
  `hpc/spec.py:328` inline skip tuple and the `webapp/server/catalog.py:27`
  `_NON_FORM_FAMILIES` frozenset, and surface it as a GUI preset instead.)
- Test: `tests/families/test_commutative_ladder_p69.py`

**Interfaces:**
- Consumes: `combinat/quiver.py::Quiver` + `Quiver.algebra(relations, field)` (the
  `IncidenceAlgebra` `"a*b - c*d"` commutativity idiom), `families/nakayama.py` loud-hint
  validation template, `families/incidence.py` for the poset route (optional — the fully-
  forward `CL(n)` IS `IncidenceAlgebra(grid poset)`, but ONLY if the poset is relabelled to
  **scalar** names first; the natural tuple `(i,j)` labels crash the knit — see the H1
  "Adjust to reality" note), `invariants/recognizers.py` degree helpers for the recognizer.
- Produces `CommutativeLadder`, `is_commutative_ladder`, `persistence_line` (signatures in
  the spec above).

- [ ] **Step 1: Write the failing tests** (deep):

```python
# tests/families/test_commutative_ladder_p69.py
"""Commutative ladders CL(n) = A_n [] A_2 (Plan 69 / R33, Escolar-Hiraoka 2016).
Literature (oracle_literature): CL rep-finite iff n <= 4 -- the n<=4 knit terminates
complete AND CL(5) is refused (the theorem boundary, literature-confirmed verbatim).
Self-cert (oracle_selfcert): the fully-forward CL(n) = incidence algebra of [n]x[2],
dim = 3*C(n+1,2); the vertex COUNTS CL(2)=11 / CL(3)=29 are SINGLE-ENGINE values asserted
against the same AR knitter that produced them -- they rest on one engine until reconciled
against the Escolar-Hiraoka figure (BLOCKING, Task 6 Step 1a) / QPA (Task 4), so they are
oracle_selfcert NOT oracle_literature (H3); the recognizer accepts constructed ladders and
rejects non-ladders; VERTICES ARE SCALAR (tuple labels crash the knit, H1)."""
import pytest
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.families.commutative_ladder import (CommutativeLadder,
                                                   is_commutative_ladder, persistence_line)
from quiverlab.modules.ar import knit_ar_quiver

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
@pytest.mark.parametrize("n, dim", [(2, 9), (3, 18), (4, 30)])   # 3*C(n+1,2)
def test_forward_ladder_dimension(n, dim):
    A = CommutativeLadder(n, field=QQ)
    assert A.dim == dim
    ok, m = is_commutative_ladder(A)
    assert ok is True and m == n


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_cl_is_rep_finite_for_n_le_4(n):
    # Escolar-Hiraoka: rep-finite for n <= 4 -- the engine corroborates by terminating
    # the knit "complete". (The vertex COUNT is a separate SELF-CERT test below; the
    # literature content here is the finite/complete BOOLEAN, not the number.)
    A = CommutativeLadder(n, field=QQ)
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)
    assert ar.is_complete and ar.status == "complete"


@selfcert
@pytest.mark.parametrize("n, count", [(2, 11), (3, 29)])   # LIVE-VERIFIED engine counts
def test_cl_verified_counts_single_engine(n, count):
    # SINGLE-ENGINE self-cert: this count is asserted against the SAME AR knitter that
    # produces it -- it rests on one engine until reconciled against the Escolar-Hiraoka
    # AR figure (BLOCKING, Task 6 Step 1a; count from the RENDERED figure, NOT pdftotext --
    # which linearizes the 2-D figure unreliably, ~30 +-1..2) and/or QPA (Task 4). This is
    # oracle_selfcert, NOT oracle_literature (H3).
    A = CommutativeLadder(n, field=QQ)
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)
    assert ar.is_complete and len(ar.vertices) == count


@selfcert
def test_cl_vertices_are_scalar_and_knit_completes():
    # REGRESSION for H1: the natural [n]x[2] grid-poset labels (i,j) are TUPLES, which
    # crash knit_ar_quiver at complex_reps.py:179 (`"e_%s" % v`, TypeError). CommutativeLadder
    # MUST use scalar vertex names, so the knit that produced the 29 count above completes.
    A = CommutativeLadder(3, field=QQ)
    assert all(not isinstance(v, tuple) for v in A.quiver.vertices)   # no tuple labels
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)       # must NOT crash
    assert ar.is_complete and len(ar.vertices) == 29


@lit
def test_cl5_is_refused_rep_infinite():
    with pytest.raises(QuiverlabError):
        CommutativeLadder(5, field=QQ)   # Escolar-Hiraoka: rep-infinite for n >= 5


@selfcert
def test_persistence_line_orientations():
    fwd = persistence_line(5, "forward")
    assert fwd.is_acyclic() and len(fwd.arrows) == 4
    zz = persistence_line(5, "zigzag")     # 1->2<-3->4<-5
    assert zz.is_acyclic() and len(zz.arrows) == 4


@selfcert
def test_recognizer_rejects_non_ladder():
    from quiverlab import Quiver
    line = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=QQ)
    ok, _ = is_commutative_ladder(line)
    assert ok is False
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError:
  quiverlab.families.commutative_ladder`

- [ ] **Step 3: Implement.** Build the box-product quiver (rows `j=1,2`, base arrows per
  `base_orientation`, rungs `c_i:(i,1)→(i,2)`, one commutativity relation per square),
  present via `Quiver.algebra`. Certify `dim` for the fully-forward case
  (`3·binom(n+1,2)`, LIVE-VERIFIED `9/18/30`); with a non-forward `τ` accept the engine's
  `dim` and an optional `expected_dim=` (loud on mismatch, the `TrivialExtension`/`Brauer`
  precedent). `CommutativeLadder(n ≥ 5)` raises `QuiverlabError` **before** building
  ("Escolar–Hiraoka: CL(n) is representation-infinite for n ≥ 5; the AR-indexed diagram is
  undefined"). `is_commutative_ladder(A)`: recognize `A_n □ A_2` by the graph shape (two
  rows of `n`, `n` rungs, `n−1` commuting squares) — return `(True, n)`; refuse loudly on a
  presentation-less algebra. `persistence_line(n, orientation)` returns a `Quiver` (reuse
  `families/dynkin.py::dynkin_quiver("A"+n, orientation=dict)` for the zigzag dict, or build
  directly).

**Adjust to reality (Task 1):**
- **MANDATE SCALAR VERTEX LABELS (H1 — a hard build gate, live-reproduced).** The natural
  `[n]×[2]` grid-poset vertices are the **tuples `(i,j)`**, and an algebra with tuple vertex
  names **crashes the AR knit**: `knit_ar_quiver` computes almost-split sequences, which build
  Ext internally (`almost_split_sequence → end_action_on_ext1 → _ext1_data → ext_cocycle_data
  → _ext_complex → _terms_info`), and `_terms_info` does `labs.index("e_%s" % v)` at
  **`complex_reps.py:179`**; with `v = (i,j)` this is `"e_%s" % (i,j)` → **`TypeError: not all
  arguments converted during string formatting`** (the `%`-with-tuple trap). The sibling
  `_vertex_basis` at **`complex_reps.py:191`** (`N.action["e_%s" % v]`, the Ext/Tor path)
  crashes identically. **Subtlety that hid this from the author:** `A.projective(v)` builds
  FINE with tuple vertices (dim 6, live-checked) — only the KNIT crashes — so a build-only
  smoke test misses it; you must knit to catch it. **LIVE-REPRODUCED 2026-08-08**
  (`QUIVERLAB_NO_NUMBA=1`): tuple-labelled grid `IncidenceAlgebra` of `[3]×[2]` → `dim 18`,
  `projective((1,1))` OK, `knit_ar_quiver` → the `complex_reps.py:179` `TypeError`; the SAME
  grid with scalar labels `"i_j"` → knit **complete, 29** (and `[2]×[2]` scalar → 11). Fix:
  build `CommutativeLadder` with **scalar** names only — either hand-build the box-product
  quiver with names like `"1_1"`/`"1_2"`/… (a flat integer index is fine too), or **relabel
  the grid poset to scalars BEFORE delegating to `IncidenceAlgebra`**. The `@selfcert`
  `test_cl_vertices_are_scalar_and_knit_completes` above is the regression (asserts no tuple
  labels AND a complete 29-vertex knit).
- **Do NOT "delegate to `IncidenceAlgebra` of the grid poset" as-is (the removed advice).**
  The clean-looking route — `IncidenceAlgebra(grid_poset)` with tuple elements — is exactly
  what crashes (above). If you reuse `IncidenceAlgebra`, you MUST pass it a **scalar-relabelled**
  poset (map each `(i,j) ↦ "i_j"`); the recognizer must then accept both the delegated and the
  hand-built form.
- **Non-monomial relations route through Gröbner** (`Quiver.algebra` sends any non-monomial
  relation to `groebner_algebra`), so `CommutativeLadder` builds are heavier than monomial
  ones — fine at these sizes (`dim ≤ 30`), but keep the batteries at `n ≤ 4` and knit with
  `budget_modules ≥ 400`. **LIVE-VERIFIED (scalar-labelled):** `CL_2` dim 9 / knit complete
  11 indec; `CL_3` dim 18 / knit complete 29 indec (a SINGLE-ENGINE self-cert count — see H3
  and `test_cl_verified_counts_single_engine`).
- **`CL(4)` dim 30 is `# PIN`** — build it and confirm knit `status == "complete"` at
  implementation (rep-finite by the theorem; the count is `# PIN`, not asserted here since
  it was not live-run at authoring). Do NOT pin a `CL(4)` count you did not run.
- **CL(n≥5) is refused at CONSTRUCTION, not by a knit-budget probe** — the theorem is the
  gate (a budget-exhausted knit is unreliable evidence of rep-infiniteness). Never try to
  knit `CL(≥5)` to "discover" it is infinite.

- [ ] **Step 4: Run tests** — PASS
- [ ] **Step 5: Commit** (`feat(families): CommutativeLadder(n) = A_n [] A_2 presented
  kQ/I + persistence_line + is_commutative_ladder; rep-finite gate n<=4 (Escolar-Hiraoka),
  loud n>=5 refusal, dim certificate 3*C(n+1,2)`).

---

### Task 2: `barcode(M)` core — A_n / zigzag interval decomposition (field-robust)

**Files:**
- Create: `src/quiverlab/modules/barcode.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.barcode(M, ...)` thin delegate,
  beside `ar_quiver`, `algebra.py:493`)
- Test: `tests/modules/test_barcode_p69.py`

**Interfaces:**
- Consumes: `modules/decompose.py::decompose` (the arbiter), `Module.dimension_vector()`,
  `modules/hom.py::is_isomorphic`, `combinat/quiver.py` (the line-order / A_n-shape check),
  `families/commutative_ladder.py::is_commutative_ladder` (routing).
- Produces `Barcode`/`Bar` (spec above), `barcode(M, ...)`, `_barcode_by_ranks(M)`
  (forward cross-check), `barcode_block(A, M)` (Task 5).

- [ ] **Step 1: Write the failing tests** (deep):

```python
# tests/modules/test_barcode_p69.py
"""Barcodes = interval decompositions of A_n / zigzag modules (Plan 69 / R33).
Literature/self-cert: the A_5 filtration barcode {[1,5],[2,2],[4,4]} (H_0 of a filtration
with dims (1,2,1,2,1); HAND-COMPUTED + LIVE-VERIFIED, identical over QQ and GF(2)).
Cross-engine: decompose == rank-formula (forward); GF(2) == QQ (field-robust bricks);
barcode == decompose summands."""
import pytest
from quiverlab import Quiver, GF
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.barcode import barcode

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert

FWD = {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 5)}


def _H0(field):
    # H_0 of a filtration: dims (1,2,1,2,1); expected barcode {[1,5],[2,2],[4,4]}.
    # blocks V1={0},V2={1,2},V3={3},V4={4,5},V5={6}; full 7x7 arrow matrices.
    A = Quiver([1, 2, 3, 4, 5], FWD).algebra(relations=[], field=field)
    o, z = field.one(), field.zero()

    def M(pairs):
        m = [[z] * 7 for _ in range(7)]
        for (i, j) in pairs:
            m[i][j] = o
        return m
    return A.module({1: 1, 2: 2, 3: 1, 4: 2, 5: 1},
                    {"a": M([(1, 0)]), "b": M([(3, 1), (3, 2)]),
                     "c": M([(4, 3)]), "d": M([(6, 4), (6, 5)])}, name="H0")


@lit
@pytest.mark.parametrize("field", [QQ, GF(2)])
def test_a5_filtration_barcode(field):
    bc = barcode(_H0(field))
    got = sorted((b.birth, b.death, b.multiplicity) for b in bc.bars)
    assert got == [(1, 5, 1), (2, 2, 1), (4, 4, 1)]     # LIVE-VERIFIED
    assert bc.field_robust is True
    assert any(b.death == 5 and b.essential for b in bc.bars)   # [1,5] essential


@xeng
def test_gf2_equals_qq_field_robust():
    q = sorted((b.birth, b.death, b.multiplicity) for b in barcode(_H0(QQ)).bars)
    g = sorted((b.birth, b.death, b.multiplicity) for b in barcode(_H0(GF(2))).bars)
    assert q == g


@xeng
def test_barcode_equals_rank_formula_forward():
    from quiverlab.modules.barcode import _barcode_by_ranks
    M = _H0(QQ)
    a = sorted((b.birth, b.death, b.multiplicity) for b in barcode(M).bars)
    b = sorted(_barcode_by_ranks(M))                    # independent rank route
    assert a == b


@selfcert
def test_interval_sum_identity_and_support_contiguous():
    M = _H0(QQ)
    bc = barcode(M)
    acc = {v: 0 for v in [1, 2, 3, 4, 5]}
    for bar in bc.bars:
        assert list(range(bar.birth, bar.death + 1)) == sorted(bar.dimvec)  # contiguous
        for v, d in bar.dimvec.items():
            acc[v] += bar.multiplicity * d
    assert acc == M.dimension_vector()                  # Sigma m_i dimvec == M


@selfcert
def test_zigzag_barcode(field=QQ):
    # 1->2<-3->4<-5 : intervals are contiguous in the line order 1-2-3-4-5.
    Z = Quiver([1, 2, 3, 4, 5], {"a": (1, 2), "b": (3, 2), "c": (3, 4), "d": (5, 4)})
    A = Z.algebra(relations=[], field=field)
    from quiverlab.modules.morphism import direct_sum
    M = direct_sum(A.projective(1), A.simple(4), A.injective(2))[0]
    bc = barcode(M)
    assert bc.kind == "zigzag" and bc.field_robust is True
    # M2: the essential flag is forward-only -- every zigzag bar has essential=False,
    # even one whose death reaches the terminal line index (no monotone "top").
    assert all(b.essential is False for b in bc.bars)


@selfcert
def test_a1_single_vertex_handled():
    # M1: a single-vertex A_1 = k has ZERO degree-1 endpoints, so the general line
    # classifier can't see it -- barcode HANDLES it as a special case (one bar [1,1]
    # with multiplicity = dim M), it does NOT refuse.
    A = Quiver([1], {}).algebra(relations=[], field=QQ)
    M = A.module({1: 3}, {}, name="V3")          # a 3-dim vector space at the one vertex
    bc = barcode(M)
    assert bc.kind == "persistence" and bc.n == 1
    got = sorted((b.birth, b.death, b.multiplicity, b.essential) for b in bc.bars)
    assert got == [(1, 1, 3, True)]              # one bar [1,1], mult 3, essential (death==n)


@selfcert
def test_loud_refusals():
    # non-A_n quiver (a vertex of degree 3): not a persistence line.
    kD4 = Quiver([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)}).algebra(
        relations=[], field=QQ)
    with pytest.raises(QuiverlabError):
        barcode(kD4.simple(0))                          # not A_n / not CL
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.modules.barcode`

- [ ] **Step 3: Implement.** `barcode(M)`:
  1. `A = M.algebra`; refuse if `A.quiver is None` (presentation-less). Classify the quiver:
     (a0) **`n == 1` special case (M1):** a single vertex with **no arrows** is `A_1` — the
     general path test below finds ZERO degree-1 endpoints (not two), so handle it explicitly:
     `line_order = (that vertex,)`, `n_top = 1`, `kind = "persistence"`; do NOT fall through to
     the loud refusal. (a) underlying graph is a path `A_n` (`n ≥ 2`: every vertex degree ≤ 2,
     exactly two degree-1 endpoints, connected) → `kind = "persistence"` if all arrows forward
     in the line order else `"zigzag"`; (b) `is_commutative_ladder(A)` → CL route (Task 3);
     (c) else **loud refusal** (barcode is for persistence lines / commutative ladders).
  2. **A_n/zigzag:** `dec = decompose(M)`; for each `(I, m)` read `dv = I.dimension_vector()`;
     assert `support(dv)` is a **contiguous interval** in the line order (else loud — the
     input is not `A_n`-shaped); `birth = min`, `death = max`. **The `essential` flag is
     forward-only (M2):** for `kind == "persistence"` set `essential = (death == n_top)` where
     `n_top` is the line's terminal index; for `kind == "zigzag"` set `essential = False`
     unconditionally (no monotone "top"). Build `bars`, set `field_robust = True`,
     `diagram = None`. For the `n == 1` case, the sole bar is `[1,1]` with `multiplicity =
     dim M` (= the number of simple summands) and `essential = True`. Decompose refusals
     (char-undecidable) propagate loudly — but for genuine `A_n`/zigzag modules the brick
     certificate makes this unreachable.
  3. `_barcode_by_ranks(M)` (forward only): the rank-formula (spec) via `mat_rank` of
     composite arrow matrices; used only by the cross-engine test. Refuse on non-forward.

**Adjust to reality (Task 2):**
- **Line order is the path order, not the vertex labels.** Compute it from the underlying
  graph (two degree-1 endpoints; BFS from the lower-labelled endpoint). Read intervals in
  THAT order so zigzag/relabelled lines still give contiguous supports. **LIVE-VERIFIED:**
  the zigzag `1→2←3→4←5` module's summands have contiguous supports in `1-2-3-4-5` order.
- **Field-robustness is REAL, not assumed** — the interval summands are bricks (`dim End =
  1`), so `decompose` certifies them in any characteristic and the split step factors over
  `GF(p)`. **LIVE-VERIFIED** identical over `QQ`/`GF(2)`/`GF(3)`. Still, keep the honest
  contract: if `decompose` ever refuses (a pathological small-field repeated-interval
  budget trip), surface it loudly — never a partial barcode.
- **The rank-formula route is a genuine second implementation** (no `decompose`, no `Hom`)
  — it is what makes `test_barcode_equals_rank_formula_forward` a real `oracle_crossengine`,
  the analogue of P57's mesh-vs-linear-algebra route. Keep it forward-only (the naive rank
  formula does not apply to zigzag; zigzag is validated by the interval-sum identity +
  Botnan–CB).
- **No `∞`.** A **forward** bar with `death == n_top` is `essential = True`, `death = n`.
  Never emit `float('inf')` or a string `"inf"`. For zigzag, `essential` is always `False`
  (M2) — a finite zigzag line has no infinite/essential bars; the endpoints are still the
  support `min`/`max` in the line order.

- [ ] **Step 4: Run tests** — PASS
- [ ] **Step 5: Commit** (`feat(modules): barcode(M) -- interval decomposition of A_n /
  zigzag persistence modules (Gabriel/Botnan-Crawley-Boevey), field-robust over any exact
  domain incl GF(2); rank-formula cross-check; A_5 filtration barcode {[1,5],[2,2],[4,4]}`).

---

### Task 3: CL generalized persistence diagram (AR-indexed, char-scoped)

**Files:**
- Modify: `src/quiverlab/modules/barcode.py` (the CL branch of `barcode` + a
  `_generalized_diagram(A, M)` helper)
- Test: `tests/modules/test_barcode_cl_p69.py` (deep)

**Interfaces:**
- Consumes: `modules/ar.py::knit_ar_quiver` (+ `ARQuiver.vertices` `{"name","dimvec",
  "module"}`, `.is_complete`, `.status`), `modules/hom.py::is_isomorphic`,
  `modules/decompose.py::decompose`, `families/commutative_ladder.py::is_commutative_ladder`.

- [ ] **Step 1: Write the failing tests** (deep):

```python
# tests/modules/test_barcode_cl_p69.py
"""Commutative-ladder generalized persistence diagrams (Plan 69 / R33). Self-cert: the
diagram's summand dim-vectors sum to M; every summand matched to an AR-quiver vertex.
Cross-engine: the diagram's distinct-indecomposable count is bounded by the knit's
indecomposable count. Honest scope: CL(n<=4) only (char 0 / char > dim); CL(n>=5) refused
at construction; a specific CL(3) worked example is # PIN vs the Escolar-Hiraoka figure."""
import pytest
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ, GF
from quiverlab.families.commutative_ladder import CommutativeLadder
from quiverlab.modules.barcode import barcode

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@selfcert
def test_cl3_diagram_sums_to_module_and_indexes_ar():
    A = CommutativeLadder(3, field=QQ)
    M = A.projective(A.quiver.vertices[0])          # a concrete CL(3) module (P of a source)
    bc = barcode(M)
    assert bc.kind == "commutative_ladder" and bc.ar_status == "complete"
    acc = {}
    for entry in bc.diagram:                        # (ar_index, ar_name, dimvec, mult, is_interval)
        assert entry.ar_index is not None          # matched to an AR vertex
        for v, d in entry.dimvec.items():
            acc[v] = acc.get(v, 0) + entry.multiplicity * d
    assert acc == {v: d for v, d in M.dimension_vector().items() if d}


@xeng
def test_cl3_diagram_indecs_are_knit_vertices():
    from quiverlab.modules.ar import knit_ar_quiver
    A = CommutativeLadder(3, field=QQ)
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)
    M = A.injective(A.quiver.vertices[-1])
    bc = barcode(M)
    assert len({e.ar_index for e in bc.diagram}) <= len(ar.vertices)   # subset of the 29


def test_cl_char_scope_refused_off_char0():
    # CL indecomposables need not be bricks -> AR knit / decompose char-scoped.
    A = CommutativeLadder(3, field=GF(2))            # char 2 <= dim -> knit/decompose refuse
    with pytest.raises(QuiverlabError):
        barcode(A.projective(A.quiver.vertices[0]))


def test_cl_char_scope_knit_surfaces_error_not_silent_complete():
    # M6: verify the refusal CHAIN, not just that barcode raises. Over GF(2) the knit
    # CATCHES the internal decompose char-refusal and returns status="error",
    # is_complete=False -- it must NEVER silently drain to a truncated "complete".
    # (LIVE-VERIFIED 2026-08-08: CL(3)/GF(2) -> status="error", is_complete=False, 8 verts.)
    from quiverlab.modules.ar import knit_ar_quiver
    A = CommutativeLadder(3, field=GF(2))
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)
    assert ar.status == "error" and ar.is_complete is False
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** the CL branch: `ar = knit_ar_quiver(A)`; if not
  `ar.is_complete` → loud refusal. **This `is_complete` check IS the char-scope refusal
  (M6):** for `n ≤ 4` over `char 0`/`char > dim` the knit returns `status="complete"`, but
  over `char ≤ dim` (e.g. `GF(2)`) the knit's internal `decompose` of a non-brick raises,
  the knit CATCHES it and returns `status="error"`, `is_complete=False` (a partial quiver,
  NOT a raise) — so `not ar.is_complete` fires and `barcode` refuses. Do NOT assume the knit
  raises; check `is_complete` / `status`. Then `dec = decompose(M)`; for each summand match
  it to `ar.vertices` by a **dimension-vector
  prefilter then exact `is_isomorphic`** (the `identify_standard` idiom over all knitted
  indecomposables) → `ar_index`, `ar_name` (`vertices[i]["name"]` or `f"X{i}"`),
  `is_interval` (thin support = contiguous 0/1). Set `diagram`, `bars` = the interval
  sub-multiset (so a CL barcode still exposes its interval bars), `field_robust = False`,
  `ar_status = ar.status`. An unmatched summand (an `is_isomorphic` refusal on all
  candidates) is a **loud error** (the knit is complete, so every indecomposable summand
  MUST be a knit vertex — an unmatched one means a bug).

**Adjust to reality (Task 3):**
- **The CL route is char-scoped (`char 0`/`char > dim`)** — `knit_ar_quiver` and
  `decompose` of a non-brick indecomposable both lean on the trace-form radical. Batteries
  run over `QQ`; `GF(2 ≤ dim)` inherits the loud `QuiverlabError` (tested). Do NOT claim
  field-robustness for CL.
- **CL(≥5) never reaches `barcode`'s CL branch** — `CommutativeLadder` refuses at
  construction. If a user hand-builds a `CL(≥5)` quiver and calls `barcode`,
  `is_commutative_ladder` returns `(True, n≥5)` and `barcode` refuses loudly
  ("Escolar–Hiraoka: rep-infinite; no complete AR-indexed diagram").
- **The Escolar–Hiraoka CL(3) worked example is `# PIN`** — ship the constructed
  `P(source)`/`I(sink)` diagrams (self-cert: sums to `M`, indexes the knit) as the durable
  oracle; if the worker can transcribe the paper's specific block-matrix example + its
  decomposition from the PDF, add it as an `@lit` pin (never from this plan's prose). The
  live-verified `CL(3) = 29`-indecomposable knit is the standing rep-finiteness self-cert.

- [ ] **Step 4: Run tests** — PASS
- [ ] **Step 5: Commit** (`feat(modules): CL generalized persistence diagram -- decompose +
  AR-quiver-indexing for CL(n<=4) (Escolar-Hiraoka), char 0/char>dim scoped; loud n>=5
  refusal`).

---

### Task 4: batteries + QPA crosscheck

**Files:**
- Test: `tests/modules/test_barcode_battery_p69.py` (deep),
  `tests/qpa/test_barcode_qpa.py` (qpa)
- **No `src/` change** — QPA has no persistence surface; the real crosschecks are a
  **direct-session `DecomposeModuleWithMultiplicities` parity** on `A_n` (the barcode's
  interval summands == QPA's decomposition summands) **and the CL(3) cross-engine closer
  (H3)** — the same `DecomposeModuleWithMultiplicities` on a concrete CL(3) module, so the
  diagram's summands are confirmed by a SECOND engine (the 29 count is otherwise
  single-engine self-cert), plus a standing no-surface guard.

**Interfaces:**
- Consumes: `barcode`, `decompose`, a spread of `A_n`/zigzag/CL algebras;
  `qpa/session.py` (`should_skip_qpa`, `require_gap`, `run`, `libgap_handle`),
  `qpa/scripts.py` (the module-to-QPA bridge used by P30's decompose crosscheck —
  `qpa_module.py::graded_form`).

- [ ] **Step 1: Battery** — (a) barcode ≡ `decompose` summand dim-vectors on a zoo of
  `A_n` interval sums (n = 2..6) over `QQ` AND `GF(2)` (field-robust); (b) the rank-formula
  ≡ barcode on forward lines; (c) essential-bar flag correctness (only `death == n` bars);
  (d) CL(3)/CL(4) diagrams sum to `M` (`# PIN` CL(4) after live run).

- [ ] **Step 2: QPA scope.** QPA has **no** persistence/barcode/commutative-ladder verb
  (probe `NamesGVars()`; the guard FAILS if one ever appears). The real parity: barcode's
  `A_n` interval decomposition vs QPA `DecomposeModuleWithMultiplicities` — QPA CAN
  decompose a module (P30 already bridges this). Reuse the P30 decompose-QPA path
  (`qpa_module.py::graded_form` + a direct session call), asserting the multiset of summand
  dim-vectors matches barcode's bars.

```python
# tests/qpa/test_barcode_qpa.py
"""QPA scope for R33 (Plan 69). QPA has NO persistence/barcode surface; the crosscheck is
DecomposeModuleWithMultiplicities parity on A_n (barcode's interval bars == QPA's summands)
+ a standing guard that FAILS if a persistence surface ever appears."""
import pytest
from quiverlab.qpa import session
pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_has_no_persistence_surface():
    lg = session.libgap_handle()
    for name in ("PersistenceDiagram", "CommutativeLadder", "Barcode"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), \
            f"QPA now ships {name} -- add a real crosscheck (honest scope changed)"
# + a DecomposeModuleWithMultiplicities parity test on kA_5 (reuse the P30 bridge).
```

- [ ] **Step 2b: QPA CL(3) cross-engine closer (H3).** The `CL(3)=29` vertex count is
  single-engine self-cert; give the DIAGRAM's decompose a second engine by running QPA
  `DecomposeModuleWithMultiplicities` on a **concrete CL(3) module** (build the scalar-labelled
  CL(3) in QPA via the P30 `qpa_module.py::graded_form` bridge — it is a presented `kQ/I`, so
  QPA accepts it) and assert the summand dim-vector multiset == `barcode(M).diagram`'s. This
  closes the DIAGRAM (the AR-index route's `decompose`) across engines even though the FULL
  29-count enumeration is not a QPA verb (QPA has no whole-AR-quiver vertex count). **If** at
  build QPA is found to expose an AR-quiver / indecomposables enumeration for CL(3), also
  compare the count and record it as the count's cross-engine closer; otherwise the 29 count
  remains `oracle_selfcert` pending the BLOCKING figure reconciliation (Task 6 Step 1a). Do
  NOT claim a QPA count crosscheck that was not run.

- [ ] **Step 3: Run** the deep battery + `-m qpa` live (the `[qpa]` extra is installed).
- [ ] **Step 4: Commit** (`test(modules,qpa): R33 barcode battery (A_n/zigzag GF(2)==QQ
  field-robust, rank-formula cross-check, CL diagrams) + QPA DecomposeModuleWithMultiplicities
  parity on A_n AND the CL(3) diagram closer + no-persistence-surface guard`).

---

### Task 5: the `barcode` module-side GUI kind

A **module-side scalar kind** (consumes the Plan-26 no-code module block, like `decompose`
/`dimension_vector`). The touchpoints follow `decompose` verbatim (a scalar module kind).
Exact symbols verified on `dev` (see the wiring map below). **The estimator IS edited (H2 —
NOT "no change"):** a CL `barcode` runs a full `knit_ar_quiver`, which is knit-heavy — exactly
the cost class the estimator already special-cases for `ar_quiver`/`tau_tilting`/
`string_homological` (`estimator.py:64-82`, the ">7 min churn" rationale). Because `barcode`
carries no `hi` budget, `sizing_dim` would size it purely off the small algebra/module dim
and mislabel a small CL barcode **instant** — where the instant tier discards artifacts AND
runs `capture_reps=False`, and the knit can churn for minutes. The fix (spec below) pushes a
CL barcode off instant while leaving the cheap `A_n`/zigzag barcode (which only `decompose`s,
no knit) instant-eligible.

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` — `MODULE_KINDS` (line 72), `_MOD_REFS` (line 98; a
  KeyError if absent), `_dispatch_module` (line 2183, add a `barcode` branch modelled on
  `decompose` line 2219), `_snip` (line 2486; a reproduction snippet).
- Modify: `docs/gui/runner.py` — `_MODULE_KINDS` (line 37), `_MOD_REFS` (line ~449),
  `_module_block` (line 622, twin branch returning `"citations"`), cost dict (line ~1425).
  **Both dispatchers MUST emit byte-identical blocks.**
- Modify: `webapp/server/schema.py` — `MODULE_KINDS` (line 129; the `_module_rules`
  validator at line 351 then auto-requires a `module` block for `barcode` — no other schema
  change).
- Modify: `webapp/server/estimator.py` — **the H2 knit-heavy fix (Option B: only CL-shaped
  inputs are knit-heavy; `A_n`/zigzag stays module-sized).** In `classify` (`estimator.py:150`),
  INSIDE the would-be-instant branch (line 165, where `ops <= instant_ops_threshold and
  max_deg <= instant_max_degree` — so the algebra is already small, bounding the cost), add a
  check BEFORE returning `"instant"`: if `"barcode" in [parse_compute_item(r).kind for r in
  req.compute]` AND the request's algebra is CL-shaped, **upgrade the tier to `"queued"` with
  `reason="knit_heavy"`** — the SAME instant→queued upgrade pattern the `report_artifacts`
  gate already uses two lines below (line 178). Add a helper `_barcode_knit_heavy(req) ->
  bool` that (a) returns `False` unless `"barcode"` is requested, then (b) builds the algebra
  DEFENSIVELY — `try: build_algebra(req.algebra.model_dump()) ... except Exception: return
  False` (the exact `_algebra_b_dim` precedent at `estimator.py:110-123`, because
  `app.py`'s `classify(...)` is NOT wrapped) and returns
  `is_commutative_ladder(A)[0]` (lazy import from `families/commutative_ladder`). Cost is
  bounded: the build only runs when the request already classifies would-be-instant (small
  dim). A non-CL `A_n`/zigzag barcode returns `False` here and stays instant-eligible
  (module-sized, byte-identical classification to today). Also extend the `KNOWN LIMITATION`
  caveat comment (`estimator.py:68-81`) with one line noting `barcode`-on-CL is knit-heavy but,
  UNLIKE `ar_quiver`/`string_homological`, is caught by this `classify` upgrade (not
  mislabelled instant).
- Modify: `src/quiverlab/trace/results_html.py` — `_HEADINGS` (line 29, `"barcode":
  "Barcode / persistence diagram"`), `_block_html` (line 960, a `barcode` branch). Render
  the barcode as (i) an **interval table** (birth | death | multiplicity | essential |
  dim-vector) via `_dims_table`/`matrix_grid` idioms and the Coxeter `interval (%s,%s)`
  precedent (`results_html.py:418`), and (ii) a simple **HTML bar diagram** — one table row
  per bar with a light-grey cell run from `birth..death` (NO new canvas; a `<table>` of
  shaded cells over the `1..n` index axis). For CL, render the AR-indexed diagram table
  (ar_name | dim-vector | multiplicity | is_interval).
- Modify: `webapp/server/i18n/{en,es,fr,zh}.json` — add `pick.kind.barcode`, `inv.barcode`,
  and `mod.barcode_*` block-label keys to **ALL FOUR** locales (exact key parity is gated
  by `tests/webapp/test_i18n.py::test_key_parity_with_english`). Suggested strings:

  | key | en | es | fr | zh |
  |---|---|---|---|---|
  | `pick.kind.barcode` | Barcode (persistence) | Código de barras (persistencia) | Code-barres (persistance) | 条形码 (持续同调) |
  | `inv.barcode` | Barcode / persistence diagram | Diagrama de persistencia | Diagramme de persistance | 持续性条形码 |
  | `mod.barcode_heading` | Barcode | Código de barras | Code-barres | 条形码 |
  | `mod.barcode_birth` | Birth | Nacimiento | Naissance | 生成 |
  | `mod.barcode_death` | Death | Muerte | Mort | 消亡 |
  | `mod.barcode_interval` | Interval | Intervalo | Intervalle | 区间 |
  | `mod.barcode_diagram` | Generalized persistence diagram | Diagrama de persistencia generalizado | Diagramme de persistance généralisé | 广义持续性图 |

- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (**byte-identical** — edit both):
  the module-kinds checkbox panel (line ~152, add `qlgui-barcode`), `S.ids` (line ~260),
  `MOD_KIND_IDS` (line ~408), `buildRequest` scalar push-list (line ~913, add `"barcode"`),
  `renderBlock` (line ~2930, a `barcode` branch mirroring `_block_html`). Optionally a
  `persistence_line`/`CommutativeLadder` **preset** button (a forward `A_5` line seed on
  the canvas). The change-handler wiring must call `scheduleProbe()` (covered by the
  `MOD_KIND_IDS` pattern).
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE golden `barcode_a5`; note it in the change-log docstring; confirm existing goldens
  byte-identical BEFORE adding).
- Test: `tests/webapp/test_barcode_kind_p69.py`, `tests/gui/test_barcode_runner_twin_p69.py`,
  and the **H2 estimator test** in `tests/webapp/test_estimator.py` (the existing estimator
  suite): a small `CL(3)` `barcode` request is **NOT `"instant"`** (`decide_tier(...) !=
  "instant"`, `classify(...)["reason"] == "knit_heavy"`), while a small `A_5` `barcode`
  request (module-sized, no knit) **stays instant-eligible** (`decide_tier(...) == "instant"`
  absent a report request). This is unmarked (webapp dir).

**Block shape** (`barcode` → `barcode_block(A, M)`):
```python
{"kind": "persistence"|"zigzag"|"commutative_ladder", "n": int, "field": str,
 "field_robust": bool,
 "bars": [{"birth": int, "death": int, "multiplicity": int, "essential": bool,
           "dimvec": {...}}, ...],
 "diagram": [{"ar_name": str, "dimvec": {...}, "multiplicity": int,
              "is_interval": bool}, ...] | None,     # CL only
 "ar_status": str|None,
 "references": ["escolar_hiraoka", "botnan_crawley_boevey", "gabriel", "assem_book"]}
# a refusal (non-A_n/non-CL, CL n>=5, presentation-less, char-undecidable) ->
#   {"error": msg, "references": [...]}, never a 500 (the module-kind per-block precedent).
```

- [ ] **Step 1: Write failing cross-runner tests** (unmarked — extras-gated dir): a
  schema-v2 request (draw `A_5`, module dims/maps = the `_H0` example, `compute:
  ["barcode"]`) asserting `block["bars"]` == the 3 bars, `"escolar_hiraoka"` in the
  citation keys, and twin parity (`json.dumps(sort_keys=True)` equality between
  `hpc/spec.py` and `docs/gui/runner.py`). **Plus the H2 estimator test** (in
  `tests/webapp/test_estimator.py`): a small `CL(3)` `barcode` request is NOT `"instant"`
  (upgraded to `"queued"`, `reason="knit_heavy"`), while a small `A_5` `barcode` request
  stays instant-eligible.
- [ ] **Step 2: Implement** both dispatch branches (shared `barcode_block` builder so they
  can't drift), the `results_html.py` render, the i18n ×4, the two `gui.js` edits, the ETA
  cost, **and the H2 `estimator.py` knit-heavy upgrade** (`_barcode_knit_heavy` +
  the `classify` instant→queued `reason="knit_heavy"` branch).
- [ ] **Step 3: Add ONE golden** (`barcode_a5`); document it; confirm existing goldens
  byte-identical. Canonical-key stability: `barcode` is a module-side kind → it
  canonicalizes through the Plan-25 `canonical_key` with the `module` block already present
  (no new top-level field); confirm an existing module request's key is byte-identical.
- [ ] **Step 4: Run** `tests/webapp/test_barcode_kind_p69.py
  tests/webapp/test_runner_delegation.py tests/gui/test_barcode_runner_twin_p69.py
  tests/webapp/test_i18n.py tests/webapp/test_js_parses.py tests/webapp/test_estimator.py
  tests/hpc -q` — PASS.
- [ ] **Step 5: Commit** (`feat(gui,webapp,hpc): barcode module-side compute kind --
  interval table + HTML bar diagram, both runners byte-identical, i18n x4, one golden;
  estimator upgrades CL barcode instant->queued (knit-heavy)`).

---

### Task 6: verification page, citations, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick the P69 card)
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`)

- [ ] **Step 1: Citations** (BibTeX-VERIFIED only; the `_r(key, bibtex_key, kind, title,
  annotation, *tags)` registry format). New keys (verify volume/page/DOI at build; arXiv
  IDs given): `escolar_hiraoka` (EscolarHiraoka2016, DCG 55(1) 100–157, arXiv:1404.7588),
  `botnan_crawley_boevey` (BotnanCrawleyBoevey2020, Proc. AMS 148(11) 4581–4596,
  arXiv:1811.08946), `igusa_rock_todorov` (IgusaRockTodorov2019, arXiv:1909.10499),
  `gabriel` (Gabriel1972, *Manuscripta Math.* 6 71–103) if not already present (grep
  first). `assem_book` (ASS2006) already ships. Optional: `asashiba_matrix`
  (AEHT 2019, JJIAM 36(1) 97–130, arXiv:1706.10027). Registry annotations lead with the
  rep-theory statement, then the TDA gloss. Verify both files agree (`bibtex()` hard-fails
  on drift) and every key resolves.

- [ ] **Step 1a: BLOCKING — reconcile the `CL(2)=11`/`CL(3)=29` counts against the
  Escolar–Hiraoka AR figures (H3).** This is a **required** implementation task, NOT a
  deferrable `# PIN`: the counts are single-engine (asserted against the same AR knitter that
  produced them), so they MUST be checked against a second source before merge. The critic's
  finding: **`pdftotext` linearizes the paper's 2-D AR-quiver figures unreliably** (it read
  ~30 for the CL(3) figure, ±1–2), so **count vertices from the RENDERED figure image**
  (render Fig. 14 — and the CL(2) figure — from the PDF to an image and count by eye/pixel),
  never from text extraction. Record the figure count and any mismatch with 29/11 in
  `docs/verification.md` (honest-scope). The cheaper machine-checkable second route is the
  **QPA CL(3) diagram closer** (Task 4 Step 2b); if QPA also enumerates the count, note it.
  If the figure and engine disagree, DO NOT silently pin either — record the discrepancy as
  a documented honest-scope item (the CRS-2004 Example-2.20 precedent) and keep the engine
  value as `oracle_selfcert`.

- [ ] **Step 2: Verification page.** Add the P69 subsystem rows:
  - `modules/barcode.py` — `oracle_literature` (the `A_5` filtration barcode
    `{[1,5],[2,2],[4,4]}`, Gabriel/Botnan–CB; the Escolar–Hiraoka rep-finiteness **BOOLEAN**
    `CL(n)` rep-finite iff `n ≤ 4` — the theorem, literature-confirmed verbatim, NOT the
    counts); `oracle_crossengine` (`decompose` ≡ rank-formula on forward lines; `GF(2) ≡ QQ`
    field-robustness; CL diagram ⊆ knit vertices; QPA `DecomposeModuleWithMultiplicities`
    parity on `A_n` + the CL(3) diagram closer); `oracle_selfcert` (the interval-sum identity
    `Σ m_i·dimvec == M`; contiguous supports; CL AR-index match; essential-bar flag; loud
    refusals; **the single-engine `CL(2)=11`/`CL(3)=29` vertex counts** — H3, self-cert until
    reconciled, NOT literature).
  - `families/commutative_ladder.py` — `oracle_selfcert` (dimension `3·C(n+1,2)`; recognizer;
    **the `CL(2)=11`/`CL(3)=29` counts** — single-engine, H3); `oracle_literature` (the
    rep-finiteness BOOLEAN: rep-finite iff `n ≤ 4`; `CL(≥5)` refused).
  - `qpa` — `DecomposeModuleWithMultiplicities` parity on `A_n` + the standing
    `NamesGVars()` guard (QPA has no persistence surface).
  - **Honest-scope entries (binding):**
    (a) **Field scope split:** the `A_n`/zigzag barcode is **field-robust** (interval
    summands are bricks, `dim End = 1`, certified in any characteristic; LIVE-VERIFIED
    `GF(2) ≡ QQ ≡ GF(3)`); the **CL route is char-scoped** (`char 0` or `char > dim` — CL
    indecomposables need not be bricks; the AR knit + non-brick decompose use the trace-form
    radical), loud off scope.
    (b) **`CL(≥5)` is a loud refusal** (Escolar–Hiraoka: representation-infinite), enforced
    at `CommutativeLadder` construction and in `barcode`'s CL branch — never a partial diagram.
    (c) **Float filtration values are OUT OF SCOPE** (exact-only): the persistence parameter
    is the discrete vertex index; optional `filtration_values` are exact `int`/`Fraction`
    labels (display-only). No approximate/continuous persistence.
    (d) **Theorem confirmed; figure COUNTS single-engine (H3).** The Escolar–Hiraoka
    **rep-finiteness theorem** (`n ≤ 4` finite, `n ≥ 5` infinite, arbitrary orientation) is
    **literature-confirmed verbatim** (`pdftotext` reads the theorem prose reliably —
    *"representation-finite for `n ≤ 4` … representation-infinite for `n ≥ 5`"*), so it is a
    real `oracle_literature` pin. The **AR-quiver-figure indecomposable counts**
    (`CL(2)=11`/`CL(3)=29`) are **`oracle_selfcert`, single-engine** — asserted against the
    same AR knitter that produces them — because `pdftotext` linearizes the 2-D figures
    unreliably (~30, ±1–2); they are reconciled against the RENDERED figure as a **BLOCKING**
    task (Step 1a) and, where possible, against the QPA CL(3) closer (Task 4). Any
    figure-vs-engine mismatch is recorded as an honest-scope item, never silently pinned. The
    **paper's specific CL(3) worked example** stays `# PIN` (ship the constructed
    `P(source)`/`I(sink)` diagrams as the durable self-cert; add the paper's block-matrix
    example only if transcribed from the PDF at build). Igusa–Rock–Todorov is cited as the
    **conceptual bridge**, not a computed oracle (quiverlab is finite/exact, not continuous).
  - **Recount the class table** (`tests/release/test_oracle_classes.py` drives the numbers
    — run collection, paste the LIVE counts, re-run to green; do NOT guess an at-authoring
    number given the mid-merge-train drift).

- [ ] **Step 3: README.** One features line: "the persistence/TDA bridge — barcodes as
  interval decompositions of `A_n` and zigzag persistence modules (Gabriel/Botnan–Crawley-
  Boevey; field-robust over `GF(2)`), and AR-quiver-indexed generalized persistence
  diagrams for commutative ladders `CL(n) = A_n □ A_2` (`n ≤ 4`, rep-finite; Escolar–
  Hiraoka), representation theory first — the R33 TDA bridge."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/modules tests/families -q -m deep` (touched deep dirs),
  `... -m pytest tests/webapp tests/gui tests/hpc -q -m fast`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q`,
  and a citation-presence check (`escolar_hiraoka`/`botnan_crawley_boevey` resolve; the
  `barcode` block carries `escolar_hiraoka`) — all green.

- [ ] **Step 5: Commit** (`docs(verification): P69 R33 oracle rows + honest scope
  (field-robust A_n/zigzag vs char-scoped CL, CL(>=5) refusal, no-float-filtration, figure
  # PIN) + escolar_hiraoka/botnan_crawley_boevey/gabriel citations + recounted classes`).

---

## GUI / webapp wiring map (verified on `dev`, 2026-08-07)

A `barcode` **module-side** kind touches 9 layers; the invariant throughout: the wheel core
(`hpc/spec.py`) and the Pyodide twin (`docs/gui/runner.py`) emit **byte-identical** blocks,
and the two `gui.js` are **byte-identical** to each other.

| Layer | File | Symbol / line |
|---|---|---|
| wheel registry | `src/quiverlab/hpc/spec.py` | `MODULE_KINDS` (72), `_MOD_REFS` (98) — mandatory |
| wheel dispatch | `src/quiverlab/hpc/spec.py` | `_dispatch_module` (2183) — add branch (model on `decompose` @2219); `_snip` (2486) |
| wheel schema | `src/quiverlab/hpc/spec.py` | auto via `parse_request` (763) once in `MODULE_KINDS` |
| Pyodide twin | `docs/gui/runner.py` | `_MODULE_KINDS` (37), `_MOD_REFS` (449), `_module_block` (622); cost dict (1425) |
| pydantic schema | `webapp/server/schema.py` | `MODULE_KINDS` (129) — auto via `_module_rules` (351) |
| server runner | `webapp/server/runner.py` | none (pure delegation) |
| i18n (×4) | `webapp/server/i18n/{en,es,fr,zh}.json` | add `pick.kind.barcode` + `mod.barcode_*` to ALL 4 |
| report HTML | `src/quiverlab/trace/results_html.py` | `_HEADINGS` (29), `_block_html` (960) — add branch; `matrix_grid`/`_dims_table`; interval idiom @418 |
| goldens | `tests/webapp/_runner_goldens.json` + `test_runner_delegation.py` | append 1 entry + docstring bullet |
| GUI JS (×2 identical) | `docs/gui/gui.js` + `webapp/static/gui/gui.js` | panel (152), `S.ids` (260), `MOD_KIND_IDS` (408), `buildRequest` push (913), `renderBlock` (2930) |
| estimator | `webapp/server/estimator.py` | `classify` (150) — **H2 EDIT**: `_barcode_knit_heavy` + instant→queued `reason="knit_heavy"` for CL-shaped barcode (Option B); `A_n`/zigzag stays module-sized via `sizing_dim` (126) |

---

## Acceptance (Plan-69 definition of done)

1. `CommutativeLadder(n, base_orientation, field)` + `is_commutative_ladder(A)` +
   `persistence_line(n, orientation)` public in `quiverlab.families` (and top-level
   `quiverlab`): built with **SCALAR vertex names** (tuple `(i,j)` labels crash the knit at
   `complex_reps.py:179` — H1; regression-tested); the forward dimension certificate
   `3·C(n+1,2)` holds (`CL(2)=9`, `CL(3)=18`, `CL(4)=30`; `CL(3)` = incidence algebra of the
   scalar-relabelled `[3]×[2]`, LIVE-VERIFIED dim 18); `CL(n ≥ 5)` refuses loudly
   (Escolar–Hiraoka rep-infinite); the recognizer accepts constructed ladders and rejects
   non-ladders + presentation-less algebras.
2. `barcode(M)` public in `quiverlab.modules.barcode` (+ `Algebra.barcode`): for `A_n`
   (any orientation) it returns the **interval decomposition** (`kind ∈
   {"persistence","zigzag"}`, `field_robust = True`), and the **`A_5` filtration barcode is
   `{[1,5],[2,2],[4,4]}`** (LIVE-VERIFIED, identical over `QQ` and `GF(2)`), with the `[1,5]`
   bar flagged **essential** (forward only; zigzag bars are always `essential=False` — M2);
   the single-vertex `A_1` case is HANDLED (one bar `[1,1]`, multiplicity = dim M — M1), not
   refused; the interval-sum identity `Σ m_i·dimvec == M.dimvec` holds; supports are
   contiguous in the line order.
3. **Field-robustness is real and gated:** the `A_n`/zigzag barcode is byte-identical over
   `QQ`/`GF(2)`/`GF(3)` (bricks, `dim End = 1`); the CL route is honestly **char-scoped**
   (`char 0`/`char > dim`, loud off it).
4. **Two independent routes agree:** `decompose` barcode ≡ the `mat_rank`-only rank-formula
   on forward lines (`oracle_crossengine`).
5. **CL generalized persistence diagram:** for `CL(n ≤ 4)`, `barcode(M)` returns the
   AR-quiver-indexed diagram (`kind = "commutative_ladder"`, `ar_status = "complete"`),
   every summand matched to a knit vertex (an unmatched one is a loud bug), the diagram sums
   to `M`; `CL(n)` is rep-finite iff `n ≤ 4` (the theorem, literature-confirmed verbatim;
   engine-corroborated: `CL(2)` knit complete, `CL(3)` complete; `CL(5)` refused).
   `CommutativeLadder` uses **scalar** vertex names (tuple labels crash the knit — H1). The
   `CL(2)=11`/`CL(3)=29` vertex counts are **single-engine `oracle_selfcert`** (H3), and
   reconciling them against the Escolar–Hiraoka AR figure (from the RENDERED image, NOT
   `pdftotext`) is a **BLOCKING** task (Task 6 Step 1a), with the QPA CL(3) diagram closer as
   the cross-engine second route; the specific CL(3) worked example stays `# PIN`.
6. **No floats, no `∞`:** endpoints are `int`; a top-reaching **forward** bar is `death = n`,
   `essential = True`; `essential` is forward-only (always `False` for zigzag — M2); float
   filtration values refused; `filtration_values` are exact `int`/`Fraction` display labels.
7. `barcode` **module-side** compute kind clickable end-to-end (GUI canvas → no-code module
   → block → report) in **all four locales** with `pick.kind.barcode` + `mod.barcode_*`
   keys (exact parity), schema-guarded (auto via `MODULE_KINDS`), both runners
   byte-identical, ONE golden with a documented change-log entry; the interval table + HTML
   bar diagram (no new canvas) + CL diagram table render; canonical keys byte-stable.
8. Live QPA battery green (`-m qpa`): `DecomposeModuleWithMultiplicities` parity on `A_n`
   (barcode's interval bars == QPA's summand multiset) + the standing `NamesGVars()` guard
   (QPA has no persistence surface).
9. `docs/verification.md` recounted (live numbers) with the four honest-scope entries
   (field-robust `A_n`/zigzag vs char-scoped CL; `CL(≥5)` refusal; no-float filtration;
   theorem-confirmed-but-COUNTS-single-engine, with the BLOCKING figure reconciliation done —
   Task 6 Step 1a — and any figure-vs-engine mismatch recorded); `escolar_hiraoka` +
   `botnan_crawley_boevey` (+ `gabriel`, `igusa_rock_todorov`) citations added and
   BibTeX-verified; README line added; deep + fast + qpa + release + citations suites green.
   No dependency on other Wave-3 plans.

---

## Methodology & assumptions

**Approach.** I read the metaplan P69 card and the full R33 record verbatim; studied two
committed exemplar plans (P59, P57) for the house form; mapped the live machinery with two
parallel Explore agents (GUI/webapp module-kind wiring; AR/module/decompose/quiver APIs)
and by reading the load-bearing source directly (`decompose.py`, `module.py`, `hom.py`,
`algebra.py::module`, `quiver.py`); and **live-verified every numeric pin I could** with the
venv Python (pure-Python kernel, light probes).

**Live-verified (observed values):**
- `A_5` interval decomposition is **byte-identical over `QQ`, `GF(2)`, `GF(3)`** (summand
  dim-vectors of `P_1⊕S_3⊕I_5⊕P_2`) — confirming the brick / `dim End = 1` field-robustness.
- The **`A_5` filtration barcode** of the hand-built `H_0` module (dims `(1,2,1,2,1)`,
  maps id/merge/id/merge) is **`{[1,5],[2,2],[4,4]}`** over both `QQ` and `GF(2)` — matching
  my hand computation exactly (elder rule).
- **Zigzag** `1→2←3→4←5`: a standard-module direct sum decomposes identically over `QQ`
  and `GF(2)`; supports are contiguous in the line order.
- **`CommutativeLadder`** (fully-forward, built as `A_n □ A_2` presented `kQ/I` with
  commuting squares, **SCALAR vertex names** — the `[n]×[2]` grid poset relabelled from tuple
  `(i,j)` to `"i_j"`, W1/H1): `CL(2)` dim 9 / knit **complete, 11 indecomposables**; `CL(3)`
  dim 18 / knit **complete, 29 indecomposables**. (`CL(3)` = incidence algebra of the
  scalar-relabelled `[3]×[2]` grid, dim 18.) **The 11/29 counts are single-engine
  (self-cert), NOT literature.** **Live-reproduced the failure the original run masked
  (2026-08-08): the NATURAL tuple-labelled grid crashes `knit_ar_quiver` at
  `complex_reps.py:179` (`"e_%s" % (i,j)` → `TypeError`)** — the counts were obtained only
  after relabeling to scalars, which the first draft did not state.
- **Char-scope refusal chain (2026-08-08):** the scalar `CL(3)` over `GF(2)` knits to
  `status="error"`, `is_complete=False` (8 partial vertices), so `barcode`'s `not
  is_complete` refusal fires — the knit CATCHES the internal `decompose` char-refusal, it
  does not raise (M6).
- `decompose`'s char guard reads `trace_rigorous = (char == 0) or (char > n)` with a
  `dim End = 1 ⇒ local` branch for any characteristic (`decompose.py:311,319`); factoring is
  supported over `GF(p)`/`QQ`/CC (`_factoring_supported`). `Algebra.module(dimension_vector,
  arrow_action, side, name)` takes FULL `n×n` arrow matrices (the no-code schema expands
  small dim[t]×dim[s] blocks; verified by the failing-then-fixed probe).

**Could NOT verify (honestly flagged in-plan as `# PIN`/`# BLOCK`):**
- The **Escolar–Hiraoka paper-figure indecomposable counts** — the counts (11, 29) are
  shipped as **single-engine `oracle_selfcert`** (asserted against the same AR knitter that
  produces them), not literature. **The critic later extracted the PDF with `pdftotext`:** it
  reads the **theorem prose reliably** — *"representation-finite for `n ≤ 4` …
  representation-infinite for `n ≥ 5`"*, arbitrary orientation — but **linearizes the 2-D AR
  figures unreliably** (~30 for CL(3), ±1–2), so the counts must be reconciled from the
  RENDERED figure (BLOCKING, Task 6 Step 1a) and/or the QPA CL(3) closer. So the
  **rep-finiteness theorem `n ≤ 4` is now LITERATURE-CONFIRMED VERBATIM** (W4 — no longer
  "abstract + search only"); only the figure VERTEX COUNTS remain open, and only as the
  single-engine-vs-figure reconciliation, not as a theorem gap.
- The **specific CL(3) worked example** — figure/algorithm-level, not machine-extractable;
  ship the constructed `P(source)`/`I(sink)` diagrams as the durable self-cert, add the
  paper's block-matrix example only if transcribed at build (`# PIN`).
- `CL(4)` count (not live-run; only `CL(2)`/`CL(3)` were) — left as `# PIN`, with `CL(4)`
  rep-finite by the theorem and its dim `30` from the closed form.
- Exact BibTeX volume/page metadata — arXiv IDs and journal names verified; final BibTeX
  freeze is the worker's build-time step (the P49/P59 precedent).

**Assumptions.** (1) The Wave-3 prerequisites (P30/P41/P37/P23-24/P26) remain at their
verified `dev` signatures at branch time — the plan re-verifies with a one-line import
before Task 1. (2) The GUI wiring line numbers (from the Explore map) are stable enough to
locate the symbols; the worker greps the symbol name, not the line, if they have drifted.
(3) The fully-forward `CL(n)` = incidence algebra of the **scalar-relabelled** `[n]×[2]`
grid (tuple labels crash the knit — H1) is the canonical shipped orientation; other
orientations `τ` build via `base_orientation` and are rep-finite (n ≤ 4) by the same theorem,
with orientation-specific counts (`# PIN`).

**Why I believe the result is correct.** The two hardest claims a critic will probe — (i)
"decompose refuses over GF(2), so barcodes can't be exact in the TDA-native field" and (ii)
"CL(≥5) must be handled" — are settled by direct evidence: (i) is refuted by the
brick / `dim End = 1` certificate (path-independent of the trace-form window) and the
`QQ ≡ GF(2) ≡ GF(3)` live runs; (ii) is a hard construction-time refusal keyed to the
Escolar–Hiraoka theorem, not a knit-budget heuristic. The headline oracle (`A_5` barcode)
is both hand-derived AND engine-confirmed, and the second (rank-formula) route makes the
barcode cross-engine-checkable without `decompose`. The honest-scope boundaries (field
split, CL char scope, no-float filtration, figure `# PIN`) are all stated as loud refusals
or verification-page entries, per the metaplan's contractual honesty rule.

## Change log

- **2026-08-07 authoring.** Initial plan (R33 persistence/TDA bridge: barcodes over
  `A_n`/zigzag + commutative-ladder generalized persistence diagrams). References
  re-verified (abstract + search; PDF/HTML figures `# BLOCK`); `A_5` barcode hand-derived +
  live-verified `{[1,5],[2,2],[4,4]}` (`QQ ≡ GF(2)`); field-robustness of the interval route
  live-verified (`QQ ≡ GF(2) ≡ GF(3)`); `CL(2)=11`/`CL(3)=29` indecomposable counts
  live-verified; GUI/webapp module-kind wiring mapped on `dev`.
- **2026-08-08 fix round (critic NEEDS WORK → all findings adjudicated valid).** Applied 7
  adjudicated rulings, DOCUMENT-ONLY (no `src/` change). **H1:** mandate SCALAR vertex labels
  for `CommutativeLadder` — the natural `[n]×[2]` grid-poset tuple labels `(i,j)` crash
  `knit_ar_quiver` at `complex_reps.py:179`/`:191` (`"e_%s" % v` `TypeError`, reached via the
  knit's internal Ext-cocycle computation; live-reproduced 2026-08-08); deleted the "delegate
  to `IncidenceAlgebra` of the grid poset" advice and added a scalar-relabel regression test
  (`test_cl_vertices_are_scalar_and_knit_completes` → 29). **H2:** removed "no estimator
  change" — the CL `barcode` runs a full knit; spec'd the `classify` instant→queued upgrade
  (`reason="knit_heavy"`, Option B: CL-shaped only, `A_n`/zigzag stays module-sized) + a
  `tests/webapp/test_estimator.py` test. **H3:** re-marked `CL(2)=11`/`CL(3)=29` as
  single-engine `oracle_selfcert` (NOT literature) everywhere incl. the verification-page
  rows; made the EH-figure reconciliation a BLOCKING task (Step 1a; count from the RENDERED
  figure — `pdftotext` linearizes the 2-D figures unreliably, ~30 ±1–2); added the QPA CL(3)
  diagram cross-engine closer. **M1:** handle `n=1` (one bar `[1,1]·dim`, not a refusal) +
  test. **M2:** `essential` is forward-only — zigzag bars are always `essential=False`
  (glossary + `Bar` schema + Step 3). **M3:** documented + live-verified the CL char-scope
  refusal chain (knit catches internal `decompose` refusal → `status="error"`,
  `is_complete=False` → `barcode` refuses; never a silent "complete") + a knit-status test.
  **W1/W4:** corrected provenance — the `dim 18` / `29` runs used the scalar-relabelled grid;
  the rep-finiteness THEOREM (`n ≤ 4`) is now literature-confirmed verbatim (`pdftotext`
  theorem-prose extraction), only the figure counts remain open. All new pins live-verified
  with the venv Python.
