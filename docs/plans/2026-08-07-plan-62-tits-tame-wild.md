# Plan 62: Tits-form tame/wild certificates + Bongartz's criterion (R19) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The combinatorial **Tits quadratic form** `q_A` of a triangular algebra
`A = kQ/I`, and the EXACT decision of its **weak positivity** / **weak
nonnegativity** by bounded integer search — then the representation-type verdict:
for a **strongly simply connected** algebra over an **algebraically closed field**,

- `A` is representation-finite ⟺ `q_A` is weakly positive (Bongartz's criterion),
- `A` is tame ⟺ `q_A` is weakly nonnegative (Brüstle–de la Peña–Skowroński),

giving the clean trichotomy **rep-finite ⊂ tame ⊂ (weakly nonnegative); wild = not
weakly nonnegative**. Every "no" carries the **exact witness dimension vector** `d ≥ 0`
with `q_A(d) ≤ 0` (isotropic ⇒ tame boundary) or `q_A(d) < 0` (⇒ wild). One no-code
compute kind — `tame_wild` — exposes the form (as a matrix), the two boolean
form-verdicts, the certificate trail (the P56 strong-simple-connectivity certificate),
the representation-type verdict, and the witness, end-to-end in all four locales.

**Two layers, cleanly separated (the record's mandate).** The **form layer** —
computing `q_A` and deciding weak positivity/nonnegativity — is **field-free** integer
arithmetic and runs for any triangular presented algebra. The **verdict layer** —
turning those booleans into rep-finite/tame/wild — is **algebraically-closed-only
(CC)** and **gated on the P56 certificate**: the tame/wild axis needs
`is_strongly_simply_connected().verdict is True`; the rep-finite axis needs
`is_simply_connected().verdict is True` (Bongartz needs only simple connectivity).
Off scope the form is still computed and reported; only the *verdict* is refused, with
an honest note. P56's three-valued `None` (Adian–Rabin / budget / undecided_char)
**propagates** to a `None` verdict here — never a fabricated tame/wild claim.

**Architecture:** A new `src/quiverlab/invariants/tits.py` holds (1) the combinatorial
Tits form `tits_form_combinatorial(A, d)` / `tits_matrix_combinatorial(A)` and the
minimal-relation counts `r_{ij}`, (2) a small internal `UnitForm` value object (n +
integer symmetric Gram matrix, `.evaluate`, `.restrict`, `.underlying_graph`), (3) the
exact decisions `is_weakly_positive` (Ovsienko box, a cited complete decision) and
`is_weakly_nonnegative` (the **classified hypercritical list** as the primary certified
route + a sound branch-and-bound witness-finder), and (4) the verdict assembler
`tame_wild_certificate(A)`. The **classified critical / hypercritical lists** live in
`src/quiverlab/invariants/_tits_lists.py` — the data that *drives* the
weak-nonnegativity decision (not merely a cross-oracle). Thin lazy-import wrapper methods
on `Algebra` — `tits_form_combinatorial`, `is_weakly_positive`, `is_weakly_nonnegative`,
`tame_wild_certificate` — sit in the `# -- invariants` block. **Name clash RESOLVED
(W1):** `Algebra.tits_form` ALREADY EXISTS (P38, the homological Euler form via `C^{-1}`,
`core/algebra.py:588`); P62's public method is `Algebra.tits_form_combinatorial(d)` and
P38's method is left untouched. The no-code `tame_wild` kind is wired through the standard
seven-touchpoint pattern.

**The named cross-plan contract (R19's explicit reuse).** `r_{ij}` = the number of
**minimal relations** from `i` to `j` = `dim_k e_j·(I / (rad·I + I·rad))·e_i`. **P56's
`fundamental_group` already computes exactly this** block-by-block (the `I/(rad·I +
I·rad)` linear algebra, Plan-56 Task 2). P62 consumes it through a thin **public**
helper `minimal_relation_counts(A) -> dict[(src,tgt), int]` that P62 adds to P56's
`invariants/coverings.py` (a five-line extraction of the per-block
`dim W_{x,y} − dim R_{x,y}` count `fundamental_group` already forms). The **independent
homological count** `r_{ij} = dim Ext²_A(S_i, S_j)` (Bongartz; Butler–King) is the
**cross-engine oracle** the two must agree on — a genuine two-implementation pin, not
a tautology.

**Distinct from P38.** P38's `forms.py::tits_form` is the **homological Euler form**
`⟨d,d⟩` via `C^{-1}` (needs a *unimodular* Cartan; = `Σ(−1)^i dim Ext^i`). The
**combinatorial Tits form** here truncates at `Ext²` (arrows − relations only) and is
defined for **every** admissible presentation, including non-unimodular / infinite
global dimension. The two **coincide iff `gl.dim A ≤ 2`** (Bongartz: for `gl.dim ≤ 2`
the Euler form's symmetrization has no `Ext^{≥3}` contribution, so it equals the Tits
form). P38's `form_type` tests **definiteness over all of ℝⁿ** (positive
definite/semidefinite); **weak positivity/nonnegativity are strictly weaker** — `q > 0`
(resp. `≥ 0`) only on the **nonnegative cone** `d ≥ 0`. A form can be weakly positive
without being positive definite; the two agree for *hereditary* algebras (where Tits =
Euler and the cone/global tests align on Dynkin/Euclidean) but not in general. This
plan does **not** reuse `form_type` for the verdict — that would be a subtle bug.

**Two verdict surfaces, kept distinct in the GUI/report (W2).** P38's `form_type` is a
*definiteness heuristic* (positive definite/semidefinite of the symmetrized matrix; a
*proven* representation-type reading only for hereditary algebras) — surfaced everywhere
as **"form type (definiteness)"**. P62's `rep_type` is the *theorem-gated, certified*
representation type (rep-finite/tame/wild via Bongartz + BdlPS, under the P56 gate over
CC) — surfaced as **"representation type (certified)"**. They are computed independently
and must never share a label; the report prints P62's line via
`results_html._block_html`'s new `tame_wild` branch (Task 7), separate from P38's
`recognizers`/`form_type` line. Where BOTH are defined (a hereditary Dynkin/Euclidean
algebra over CC) they AGREE (finite↔rep-finite, tame↔tame, wild↔wild) — pinned as a
self-cert cross-check (Task 5) and as the `gl.dim ≤ 2` agreement test (Task 1).

**Tech Stack:** pure exact integer arithmetic (Python `int`; the Gram matrix is an
integer sympy matrix only for display). The minimal-relation counts come from P56's
exact field-linear-algebra (`fields.linalg`) via `coverings.minimal_relation_counts`;
the Ext² cross-check uses the existing `Algebra.ext`. No floats in `src/` — the form,
the box search, and the witness are all integer/`dict`.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **Bucket decision (stated, not guessed):** new tests live in **`tests/invariants/`**,
  which `tests/conftest.py::_bucket` (`_DEEP_DIRS`) does **not** list ⇒ they collect in
  the **fast** bucket (the P38/P49 `invariants/` precedent). Keep every battery fast:
  the weak-positivity box sweep runs only on the **small Dynkin/Euclidean oracle forms**
  (n ≤ 9: `E8`/`Ẽ8`); the weak-nonnegativity decision is a finite **classified-list
  match** (no sweep), so the wild 10-vertex `T_{2,3,7}` is decided by a list entry, not a
  sweep; the m-Kronecker ladder is n = 2; the strong-s.c.-gated end-to-end oracles use
  hand-sized trees so P56's convex sweep stays inside the fast budget. QPA
  honest-skip goes in `tests/qpa/` (bucket = the class); cross-runner webapp/gui tests
  go in `tests/webapp/` (unmarked — extras-gated dir, Plan-32 rule).
- New scalar kind acts on the **existing algebra block** ⇒ **schema stays v1**
  (Plan-26/38/56 precedent: request-schema bumps only for new request *blocks*).
- **Canonical keys are pinned to the frozen `_V` constant by design** (metaplan §3):
  adding `tame_wild` to the compute list must not re-key any existing request. Verify by
  running `test_runner_delegation.py` BEFORE and AFTER adding the golden; existing
  entries stay byte-identical.
- **Every citation key is BibTeX-verified before use** (Plan-29 rule); `_r(...)` registry
  entry + `references.bib` entry both land, gated by
  `tests/citations/test_bib_structure.py`.
- All refusals are loud (`QuiverlabError(message, hint=...)`, the two-space `[hint: ...]`
  convention, `errors.py`); presentation-less (structure-constant) algebras refuse the
  form and every verdict — the quiver is required — mirroring the `dynkin_type` guard.
- **P56 is a HARD prerequisite** (metaplan wave map: "P62 needs P56"). P56 merges to
  `dev` **before** P62 branches; P62 designs against the merged `invariants/coverings.py`
  — the `StrongSimpleConnectivity` dataclass (`verdict: bool|None`, `witness`, `reason`,
  `checked_convex`), `SimpleConnectivity` (`verdict`, `strongly`),
  `is_strongly_simply_connected`, `is_simply_connected`, and the new
  `minimal_relation_counts` helper P62 adds there. **If P56 has not merged when a worker
  picks this up, STOP and escalate** — do not fork P56's minimal-relation machinery.
- **P38 is also merged** (`invariants/forms.py`, `dynkin_type.py`, `roots.py`) and
  consumed read-only (naming/coexistence handled in Task 1).
- Conventional commits; green at every commit; branch `plan-62-tits-tame-wild` off `dev`.

---

## Record (verbatim — the mandate for this plan)

> **R19 — Tame/wild certificate via the Tits form (strongly simply connected).**
> [C-scout P7; keep-with-corrections] Object: weak positivity / weak nonnegativity of
> q_A decided by EXACT bounded integer search (primary route; the printed
> hypercritical/critical lists — de la Peña Banach Center 26 (1990), von Höhne, Unger,
> and the Barot–Jiménez-González–de la Peña 2019 book — as transcription-checked
> cross-oracles); verdicts: rep-finite ⇔ weakly positive, tame ⇔ weakly nonnegative.
> SCOPE: strongly simply connected (gate = R16) over ALGEBRAICALLY CLOSED fields — the
> verdict layer is CC-only, the form computations are field-free. Companion (critic's
> top find): Bongartz's criterion — Math. Ann. 269 (1984) 1–12 — rep-finite ⇔ weakly
> positive for simply connected, same machinery. Refs: Brüstle–de la Peña–Skowroński
> Adv. Math. 226 (2011) 887–951; arXiv:1905.06028. Oracles: m-Kronecker ladder (1/2/≥3),
> Dynkin/Euclidean trees, T_{2,3,7}. Size L.

---

## Reference re-verification (mandatory; findings recorded)

All statements below were re-read from primary/secondary sources during authoring
(WebSearch + Encyclopedia of Mathematics + the arXiv abstracts/pages of the cited
papers; the BJP 2019 book chapter and the compressed arXiv PDFs could not be
byte-extracted — those numbers carry a `# PIN`). Keep this section — the implementer
transcription-checks every `# PIN` against a primary source.

- **Tits form definition — VERIFIED (Encyclopedia of Mathematics, "Tits quadratic
  form"; corroborated by Kasjan–Skowroński arXiv:1905.06028 and the general search).**
  For a basic `A = kQ/I` with `I` admissible,
  `q_A(x) = Σ_{i∈Q_0} x_i² − Σ_{(i→j)∈Q_1} x_i x_j + Σ_{i,j∈Q_0} r_{ij} x_i x_j`,
  where `r_{ij}` = the number of **minimal relations** from `i` to `j` and **does not
  depend** on the chosen minimal generating set `R` of `I`. Homological identity:
  `r_{ij} = dim_k Ext²_A(S_i, S_j)` (the Encyclopedia writes `Ext²(S_j, S_i)` — that is
  the opposite-module index convention; in THIS package's convention arrows `i→j`
  correspond to `Ext¹(S_i, S_j)` (`forms.py` docstring, verified), so relations `i→j`
  correspond to `Ext²(S_i, S_j)` — we use the package convention throughout and pin the
  agreement by test). For a **triangular** (acyclic) `Q` there are no loops or oriented
  cycles, so `r_{ii} = 0` and `q_A` is a **unit form** (all diagonal coefficients `= 1`).
- **Bongartz's criterion — VERIFIED (Math. Ann. 269 (1984) 1–12, doi
  10.1007/BF01455993, "A criterion for finite representation type").** Bongartz
  introduced the Tits form of a basic algebra with acyclic quiver and proved: if `A` is
  **representation-finite** then `q_A` is **weakly positive**; and for a **simply
  connected** algebra the following are equivalent: (i) `A` is representation-finite;
  (ii) `q_A` is weakly positive; (iii) `A` has no convex subcategory that is critical.
  (Simple — not strongly-simple — connectivity suffices for the rep-finite axis.)
- **Brüstle–de la Peña–Skowroński — VERIFIED (Adv. Math. 226 (1) (2011) 887–951, "Tame
  algebras and Tits quadratic forms").** For a **strongly simply connected** algebra `A`
  (over an algebraically closed field): `A` is **tame** ⟺ `q_A` is **weakly
  nonnegative**. This is the theorem P62 consumes downstream of the P56
  strongly-simply-connected certificate. (Corroborated verbatim in
  arXiv:1905.06028 Kasjan–Skowroński §1, the P56-cited source.)
- **Ovsienko's theorem — VERIFIED (statement; sources: BJP 2019 book Ch. 5 "Weakly
  Positive Quadratic Forms"; von Höhne Comment. Math. Helv. 63 (1988) 312–336 "On weakly
  positive unit forms").** For a **weakly positive unit form** `q` in `n` variables,
  every **positive root** (`d ≥ 0`, `q(d) = 1`) has all coordinates `≤ 6`. **Decision
  (rested on the CITATION, not re-proved here):** a unit form `q` is weakly positive ⟺
  `q(d) > 0` for every `d ∈ {0,…,6}ⁿ, d ≠ 0` (BJP 2019 Ch. 5, the box test). **CORRECTION
  (adversarial review, 2026-08-07):** an earlier draft carried a closed "critical ⇒
  Euclidean isotropic radical, entries ≤ 6" proof of the box test — that step is **FALSE**.
  The critical unit forms (minimal not-weakly-positive) are NOT all Euclidean: the
  `m ≥ 3` Kronecker `q = x² + y² − m·xy` is critical for weak positivity (both
  1-variable restrictions are weakly positive) yet is **INDEFINITE with no isotropic
  radical** — its witness has `q(1,1) = 2 − m < 0`, not `0`. The box-6 decision is
  nonetheless the classical **cited** result (Ovsienko; BJP 2019 Ch. 5); we rest on the
  citation and do not reproduce a proof. Independent of any bound, a `False` is always a
  **sound exact witness** (`d ≥ 0`, `q(d) ≤ 0`).
- **von Höhne hypercritical classification — CORRECTED (adversarial review, 2026-08-07):
  the "≤ 9 variables" transcription was REFUTED.** A **hypercritical** unit form is not
  weakly nonnegative but every proper restriction is; they form a **finite
  explicitly-classified list** (von Höhne, Proc. London Math. Soc. (3) 73 (1996) 47–67).
  An earlier draft claimed "hypercritical unit forms have ≤ 9 variables" (mis-transcribed
  from a specific reduction-algorithm's subgraph bound). **This is false — refuted by
  this plan's own oracle `T_{2,3,7}`:** the star with arms `1,2,6` (10 vertices) is wild
  (`1/2 + 1/3 + 1/7 = 41/42 < 1`), yet every proper restriction is Dynkin or affine `Ẽ8`
  (weakly nonnegative) — so **`T_{2,3,7}` is itself a 10-variable hypercritical form**
  (verified in read-only Python this revision: `min q = −1`, at the **sincere** defect
  vector `(12,6,8,4,10,9,7,6,4,2)`, support all 10 vertices, entries up to 12). Hence
  **hypercritical forms reach at least 10 variables and their defect vectors need not be
  small**, so there is **NO safe universal support cap or entry-box** for weak
  nonnegativity. **Design consequence (the ruling):** the weak-nonnegativity decision is
  driven by the **classified finite hypercritical list**, NOT a guessed box; a certified
  `True` (tame) is emitted ONLY where the encoded list is complete for the algebra's
  vertex count; otherwise the verdict is honest `None`. Every support/entry bound used
  anywhere is **read off the encoded list**, never assumed. Reference: von Höhne 1996 /
  de la Peña Banach Center 26 (1990) / BJP 2019 — **# PIN** the full list (Task 3 encodes
  what is transcription-checked and records coverage).
- **Classified critical / hypercritical lists — the DECISION data (not merely a
  cross-oracle).** Critical (weak-positivity) list: the `m`-Kronecker family + the
  Euclidean forms `Ã_n, D̃_n, Ẽ_{6,7,8}` — but weak positivity is decided by Ovsienko's
  box-6 above, so this list is only a *cross-oracle* there. Hypercritical
  (weak-nonnegativity) list: von Höhne 1996 / de la Peña Banach-26 / BJP 2019 — **this is
  the PRIMARY decision data** for the tame axis. Encode every entry transcription-checked
  against a primary source (Gram matrix + defect vector, `q(defect) < 0` verified) and
  record `HYPERCRITICAL_COVERAGE` (the vertex counts / families for which the list is
  provably complete). The minimal wild TREES `T_{p,q,r}` with `1/p+1/q+1/r < 1` — e.g.
  `T_{2,3,7}` (10 vertices, defect above), `T_{2,4,5}` (9 vertices, `min q = −3`),
  `T_{3,3,4}` (8 vertices, `min q = −6`) — are transcribable now (defect vectors verified
  in read-only Python this revision). **# PIN** the non-tree hypercritical forms against a
  primary source; mark the rest deferred (coverage recorded, verdict honest-`None`
  outside coverage — never a guessed `True`).
- **The scope hypotheses — VERIFIED.** Both theorems are over an **algebraically closed
  field** (the tame/wild dichotomy, Drozd, is an algebraically-closed notion). In
  quiverlab the only algebraically closed field is **`CC`** (`fields/complexfield.py`);
  `GF(p)`, `GF(p^n)`, `QQ`, `QQi` are **not** — the verdict layer refuses over them
  (form still computed). The form coefficients `r_{ij}` are computed over the algebra's
  actual field via P56; the weak-positivity/nonnegativity **decision** of the resulting
  integer form is field-free.

---

### Task 0: citations (do first — later tasks reference the keys)

**Files:**
- Modify: `src/quiverlab/citations/references.bib`, `src/quiverlab/citations/registry.py`
- Test: `tests/citations/test_bib_structure.py` (existing gate)

**Interfaces:** the `_r(key, bibtex_key, kind, title, annotation, *tags)` registry helper
+ the `@article`/`@book`/`@incollection` bib patterns. `assem_book` (→ `ASS2006`) already
exists; do NOT re-add. `skowronski_ssc` and `kasjan_skowronski` (arXiv:1905.06028) **may
already have been added by P56** — check first and reuse, do not duplicate keys.
`_citation_pairs` (`hpc/spec.py`) resolves keys → `[key, formatted]` pairs for payloads.

- [ ] **Step 1: add the keys** (BibTeX-verify each before committing — Plan-29 rule; skip
  any P56 already added):
  - `bongartz_criterion` → Bongartz, *A criterion for finite representation type*, Math.
    Ann. **269** (1984) 1–12, doi 10.1007/BF01455993. (the rep-finite ⟺ weakly-positive
    companion; the Tits form's origin)
  - `bdps_tame_tits` → Brüstle, de la Peña, Skowroński, *Tame algebras and Tits quadratic
    forms*, Adv. Math. **226** (2011) 887–951, doi 10.1016/j.aim.2010.07.007. (tame ⟺
    weakly nonnegative, strongly simply connected)
  - `kasjan_skowronski` → Kasjan, Skowroński, *On tame strongly simply connected
    algebras*, arXiv:1905.06028. (the theorem statements P62 consumes; reuse if P56 added)
  - `ovsienko_forms` → Ovsienko, *Integral weakly positive forms*, in Schur Matrix
    Problems and Quadratic Forms, Inst. Mat. Akad. Nauk Ukrain. SSR, Preprint 78.25 (1978)
    3–17. (the ≤ 6 bound) — # verify the exact preprint pagination; if unresolved, cite
    via the BJP 2019 book Ch. 5 which restates it.
  - `vonhohne_wnn` → von Höhne, *On weakly non-negative unit forms and tame algebras*,
    Proc. London Math. Soc. (3) **73** (1996) 47–67. (the hypercritical classification —
    the PRIMARY weak-nonnegativity decision data; NOT bounded to ≤ 9 variables, see the
    `T_{2,3,7}` refutation in the reference section)
  - `delapena_banach26` → de la Peña, *Algebras with hypercritical Tits form*, in Topics
    in Algebra, Banach Center Publ. **26** Part 1, PWN (1990) 353–369. (printed
    hypercritical list) — # verify exact pages.
  - `bjp_quadratic_forms` → Barot, Jiménez-González, de la Peña, *Quadratic Forms:
    Combinatorics and Numerical Results*, Algebra and Applications **25**, Springer (2019),
    doi 10.1007/978-3-030-05627-8. (the search-box reference; critical/hypercritical
    exposition) — `@book`.
  - `unger_wild` → Unger, *The concealed algebras of the minimal wild hereditary
    algebras*, Bull. London Math. Soc. **22** (1990) 587–592. (Unger's wild-boundary
    tables) — # verify exact title/pages; OPTIONAL, only if a table is actually encoded
    in `_tits_lists.py`.
- [ ] **Step 2:** `... -m pytest tests/citations/test_bib_structure.py -q` green (registry
  ↔ bib in sync; every `_r` key resolvable via `bibtex(key)`).
- [ ] **Step 3: Commit** — `docs(citations): Tits-form tame/wild references (Bongartz 1984, BdlPS 2011, Ovsienko, von Höhne 1996, de la Peña Banach-26, BJP 2019)`

---

### Task 1: the combinatorial Tits form + `minimal_relation_counts` (the P56 contract)

**Files:**
- Create: `src/quiverlab/invariants/tits.py`
- Modify: `src/quiverlab/invariants/coverings.py` (add the public `minimal_relation_counts`
  extraction), `src/quiverlab/invariants/__init__.py`, `src/quiverlab/core/algebra.py`
- Test: `tests/invariants/test_tits_form.py`

**Interfaces:**
```python
# invariants/coverings.py  (P56's module — P62 adds this thin PUBLIC helper)
def minimal_relation_counts(A) -> dict:
    # {(src, tgt): count} for every ordered vertex pair with >= 1 minimal relation.
    # = dim_k e_tgt (I / (rad.I + I.rad)) e_src, block-by-block -- the SAME quantity
    # fundamental_group computes in Plan-56 Task 2. Presentation-less A: loud
    # QuiverlabError (needs the quiver). Empty dict for a hereditary (relation-free)
    # algebra. (W3: this is a small REFACTOR of P56, not a five-line extraction -- see
    # the honest note in Step 3 and the P56 addendum.)

# invariants/tits.py
def tits_form_combinatorial(A, d) -> int:      # W1: NOT tits_form (P38 owns that name)
    # q_A(d) = sum_i d_i^2 - sum_{arrows i->j} d_i d_j + sum_{(i,j)} r_ij d_i d_j
    # d a vertex-order list/tuple/dict; exact int. Field-free given r_ij.
def tits_matrix_combinatorial(A) -> "sympy.Matrix":
    # the integer symmetric Gram matrix G with q_A(d) = (1/2) d^T G d:
    #   G_ii = 2 ; G_ij = G_ji = r_ij + r_ji - a_ij - a_ji  (i != j)
    # (for triangular Q only one of the i->j / j->i directions is ever nonzero).
    # For DISPLAY / the payload; the decisions use UnitForm below.
def as_unit_form(A) -> "UnitForm":
    # RAISE loudly (QuiverlabError) unless A.quiver.is_acyclic() -- the Tits-form
    # tame/wild theory is defined for TRIANGULAR algebras only:
    #   "the combinatorial Tits form is a unit form only for a triangular (acyclic)
    #    quiver; this quiver has an oriented cycle/loop".
    # M2 (adversarial review): do NOT gate on G_ii == 2 -- that check WRONGLY PASSES on
    # k[x]/(x^3) (loop x: a_11 = 1 and the minimal relation x^3 gives r_11 = 1, so
    # G_11 = 2 - 2*1 + 2*1 = 2). is_acyclic() is the honest gate.

@dataclass(frozen=True)
class UnitForm:
    n: int
    gram: tuple          # tuple-of-tuples, symmetric integer, diagonal == 2 (acyclic-gated)
    labels: tuple        # vertex labels in index order (for witness readout)
    def evaluate(self, d) -> int          # (1/2) d^T gram d, exact int
    def restrict(self, support) -> "UnitForm"     # principal submatrix on `support`
    def underlying_graph(self) -> dict    # adjacency: i~j iff gram[i][j] != 0 (i!=j)
```

**Construction (field-free, exact):**
1. `minimal_relation_counts(A)` — a **small refactor** of P56 (W3): P56's
   `fundamental_group` forms the per-`(x,y)` count `dim W_{x,y} − dim R_{x,y}` but folds it
   straight into SNF rows without retaining a dict, so this is not a pure five-line
   extraction — P56 must RETAIN and expose the per-pair counts (committed in the P56
   addendum). Route `fundamental_group` through the retained dict (single source of truth,
   P56 byte-unchanged).
2. `tits_form_combinatorial` / `tits_matrix_combinatorial`: read `A.quiver.arrows` for
   `a_ij` and `minimal_relation_counts(A)` for `r_ij`; assemble the integer form directly.
   Guard `A.quiver is None` → loud refusal.
3. `as_unit_form` builds `UnitForm(n, gram, labels)` in `A.quiver.vertices` order; **raises
   unless `A.quiver.is_acyclic()`** (M2 — the honest triangular gate; the `G_ii == 2`
   heuristic is unsound, see above). For an acyclic quiver `G_ii = 2` automatically
   (no loops, no `i→i` relations).

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_tits_form.py
"""The combinatorial Tits form q_A = sum x_i^2 - sum_arrows + sum r_ij, r_ij the
minimal-relation count (P56's I/(rad.I+I.rad)) cross-checked against dim Ext^2(S_i,S_j).
Refs: Bongartz Math. Ann. 269 (1984); Encyclopedia of Mathematics 'Tits quadratic form'."""
import pytest
from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.coverings import minimal_relation_counts
from quiverlab.invariants.tits import as_unit_form, tits_form_combinatorial as titsc

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


@xeng
def test_r_ij_matches_ext2_of_simples():
    # the named cross-plan contract: combinatorial count (P56) == homological Ext^2.
    # 1->2->3 with the length-2 relation a*b=0 : ONE minimal relation 1->3.
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=GF(7))
    r = minimal_relation_counts(A)
    assert r.get((1, 3), 0) == 1 and sum(r.values()) == 1
    assert A.ext(A.simple(1), A.simple(3), 2) == 1          # independent engine
    assert A.ext(A.simple(1), A.simple(2), 2) == 0


@lit
def test_kronecker_tits_values():
    # m-Kronecker (hereditary): q(x,y) = x^2 + y^2 - m x y. r_ij all 0.
    for m, val11 in [(1, 1), (2, 0), (3, -1), (4, -2)]:
        A = Quiver([1, 2], {f"a{k}": (1, 2) for k in range(m)}).algebra(field=CC)
        assert minimal_relation_counts(A) == {}
        assert titsc(A, [1, 1]) == val11                     # 2 - m


@lit
def test_A2_and_nakayama_relation_form():
    A2 = Quiver([1, 2], {"a": (1, 2)}).algebra(field=CC)
    assert titsc(A2, [1, 1]) == 1                            # x^2+y^2-xy, a root
    # kA3 with a*b=0 : q = x1^2+x2^2+x3^2 - x1x2 - x2x3 + x1x3 (one relation 1->3)
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=CC)
    assert titsc(A, [1, 1, 1]) == 3 - 1 - 1 + 1              # == 2


@xeng
def test_gldim_le_2_agrees_with_P38_euler_symmetrization():
    # W2: where BOTH are defined (gl.dim <= 2), the combinatorial Tits form == P38's
    # homological Euler symmetrization (2 * C^-1 sym). The self-cert cross-check.
    from quiverlab.invariants.forms import tits_matrix as euler_sym
    A = linear_path_algebra(4, field=CC)                     # hereditary, gl.dim 1
    G = as_unit_form(A).gram
    import sympy as sp
    assert sp.Matrix(G) == euler_sym(A)                      # 2*(C^-1 sym) == our G


def test_non_triangular_refused():
    # M2: k[x]/(x^3) has a loop => G_11 == 2 by cancellation (a_11=1, r_11=1), so the
    # G_ii==2 heuristic would WRONGLY accept it. is_acyclic() is the honest gate.
    A = truncated_polynomial(3, field=CC)                    # k[x]/(x^3): a loop
    with pytest.raises(QuiverlabError, match="triangular|acyclic|cycle"):
        as_unit_form(A)


def test_presentationless_refused():
    from quiverlab.families import some_structure_constant_algebra  # or a hand SC algebra
    ...  # presentation-less Algebra; minimal_relation_counts + tits_form_combinatorial raise
```

- [ ] **Step 2: Run to verify failure** (`ModuleNotFoundError: tits` / `AttributeError:
  minimal_relation_counts`).
- [ ] **Step 3: Implement** the `coverings.minimal_relation_counts` refactor (W3 — have
  P56 RETAIN the per-`(x,y)` count dict it currently folds into SNF rows, and expose it;
  route `fundamental_group` through the retained dict — re-run the P56 suite
  `tests/coverings/` to confirm byte-unchanged; the P56 addendum records this contract),
  then `tits.py`. **Adjust to reality:** confirm `A.quiver.arrows` is `label → (s, t)`;
  confirm the Ext index convention on a hand example (arrows `i→j` ↔ `Ext¹(S_i,S_j)` per
  `forms.py`) and orient `r_ij` to match — the `test_r_ij_matches_ext2` pin is the arbiter,
  adjust the `(src,tgt)` orientation until it holds, never the pin. Fill the
  presentation-less fixture from an actual structure-constant `Algebra` in the tree.
- [ ] **Step 4: Run** — `... -m pytest tests/invariants/test_tits_form.py tests/coverings -q` green.
- [ ] **Step 5: Commit** — `feat(invariants): combinatorial Tits form + minimal_relation_counts (P56 I/(rad.I+I.rad) reuse) cross-checked vs Ext^2`

---

### Task 2: weak positivity (Ovsienko box) + weak nonnegativity (classified list)

**Files:**
- Modify: `src/quiverlab/invariants/tits.py`, `src/quiverlab/core/algebra.py`
- Test: `tests/invariants/test_weak_positivity.py`

**Interfaces:**
```python
@dataclass(frozen=True)
class FormVerdict:
    holds: bool | None       # True / False (witness found) / None (undecided --
                             #   budget_exceeded OR hypercritical list incomplete; never
                             #   a silent/guessed True)
    witness: tuple | None    # on False: the exact d >= 0 (vertex-order tuple) with
                             #   q(d) <= 0 (weak pos) / q(d) < 0 (weak nonneg)
    witness_value: int | None
    reason: str              # "weakly positive" | "witness d=.. q(d)=.." |
                             #   "budget_exceeded (n=..)" | "hypercritical list partially
                             #   transcribed -- no certified weakly-nonnegative verdict"
    checked: int             # lattice points evaluated (0 for a pure list decision)

def is_weakly_positive(form, *, budget=3_000_000) -> FormVerdict
    # Ovsienko (cited): q weakly positive iff q(d) > 0 for all 0 != d in {0..6}^n.
    # Realized as a pruned CONNECTED-SUPPORT BRANCH-AND-BOUND sweep (below). Complete
    # decision (True/False), or None on budget (never a guessed True).

def is_weakly_nonnegative(form, *, budget=3_000_000) -> FormVerdict
    # PRIMARY = the classified HYPERCRITICAL LIST (Task 3), NOT a guessed box:
    #   (1) run the bounded branch-and-bound WITNESS-FINDER; any d >= 0 found with
    #       q(d) < 0 is a SOUND False (exact witness, no bound assumption);
    #   (2) else, if some connected restriction q|_S matches a listed hypercritical
    #       form, return False with that entry's tabulated defect vector;
    #   (3) else, if _tits_lists.HYPERCRITICAL_COVERAGE covers every vertex count <= n
    #       (list provably complete here), return True (certified -- no hypercritical
    #       restriction exists);
    #   (4) else return None, reason "hypercritical list partially transcribed ...".
    # _WNN_BOX is DELETED as a decision mechanism: NO True from a guessed box, EVER.
```

**Branch-and-bound sweep (M1: real budget arithmetic).** The raw box `{0..6}ⁿ` is `7ⁿ`
(`≈ 2.8·10⁸` at `n = 10`) — infeasible naively. Two exact accelerations: (a) a
**minimal-support** witness has **connected** support (a disconnected support splits `q`
additively — one connected piece already witnesses), so enumerate connected vertex subsets
`S` (BFS grow-and-prune, dedup) and search only **sincere** `d ∈ {1..6}^S`; (b)
**branch-and-bound** — build `d` coordinate-by-coordinate over `S` and prune a partial
assignment whose best-case completion cannot reach `q ≤ 0` (a valid lower bound on `q` over
the remaining nonnegative coordinates: e.g. drop the already-fixed positive contributions
and bound the cross terms). **Costed on the named oracles:** `E8` (8 vertices, weakly
positive, *no* witness ⇒ worst case exhausts the sincere box, `6⁸ = 1 679 616`, but
positive-definiteness prunes hard — the completion lower bound rises to `> 0` early);
`Ẽ8` (9 vertices, only witness = the sincere isotropic radical `(6,3,4,2,5,4,3,2,1)`, max
entry 6, `6⁹ ≈ 1.0·10⁷` worst case but a center-first coordinate order finds it early).
So the **weak-positivity** budget must clear `E8`'s `1.68M` full product; default
`3_000_000` covers `E8` outright and `Ẽ8` with B&B pruning + early witness return. The
**weak-nonnegativity** path does NOT box-sweep the wild cases (they are decided by the
list) — its witness-finder is used only for SMALL witnesses (e.g. the `3`-Kronecker `(1,1)`),
so its budget is never the bottleneck.

**Soundness note (binding).** A `False` is always a *true* witness (an explicit `d ≥ 0`
with `q(d) ≤ 0` / `< 0`), independent of any bound — the bound only affects *completeness*.
Ovsienko's 6 makes weak **positivity** a complete decision (True/False, or None on budget).
Weak **nonnegativity** is complete **only where the encoded hypercritical list is complete
for the vertex count** (`HYPERCRITICAL_COVERAGE`); outside coverage the verdict is honest
`None` — **never a guessed True**. All support/entry bounds anywhere are read off the list,
never assumed (the `T_{2,3,7}` refutation: no universal cap exists).

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_weak_positivity.py
"""Weak positivity by the Ovsienko box (cited, branch-and-bound); weak nonnegativity by
the classified hypercritical list (primary) + a sound witness-finder. Oracles: Dynkin =>
weakly positive; Euclidean (extended Dynkin) => weakly nonnegative not weakly positive
(isotropic witness q=0); 3-Kronecker (small witness) and T_{2,3,7} (LIST entry, sincere
10-var witness) => not weakly nonnegative. Refs: Bongartz 1984; von Hohne 1996; BJP 2019."""
import pytest
from quiverlab import CC, Quiver, linear_path_algebra
from quiverlab.invariants.tits import (as_unit_form, is_weakly_nonnegative,
                                       is_weakly_positive)
from quiverlab.invariants.tits import tits_form_combinatorial as titsc

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _kron(m):
    return Quiver([1, 2], {f"a{k}": (1, 2) for k in range(m)}).algebra(field=CC)


def _star(arms):
    # star tree T_{p,q,r}: center 0 with `arms` = list of arm edge-lengths.
    verts, ars, nxt = [0], {}, 1
    for ai, length in enumerate(arms):
        prev = 0
        for _ in range(length):
            verts.append(nxt); ars[f"e{ai}_{nxt}"] = (prev, nxt); prev = nxt; nxt += 1
    return Quiver(verts, ars).algebra(field=CC)


@lit
def test_dynkin_weakly_positive():
    for A in (linear_path_algebra(5, field=CC), _star([1, 2, 4])):   # A5, E8 (n<=8)
        f = as_unit_form(A)
        assert is_weakly_positive(f).holds is True
        assert is_weakly_nonnegative(f).holds is True                # positive => nonneg


@lit
def test_euclidean_weakly_nonnegative_not_positive():
    # 2-Kronecker = ~A1 (isotropic (1,1)); affine E8 = ~E8 = star arms [1,2,5] (9 verts)
    for A in (_kron(2), _star([1, 2, 5])):
        f = as_unit_form(A)
        wp = is_weakly_positive(f)
        assert wp.holds is False and titsc(A, list(wp.witness)) <= 0   # isotropic root
        assert is_weakly_nonnegative(f).holds is True                  # certified (list complete)
    assert titsc(_kron(2), [1, 1]) == 0                                # the isotropic delta


@lit
def test_wild_not_weakly_nonnegative():
    # 3-Kronecker (n=2, small witness found by the finder) and T_{2,3,7}=star arms [1,2,6]
    # (n=10, decided by the LIST entry -- its witness is SINCERE on all 10 vertices,
    # (12,6,8,4,10,9,7,6,4,2), q=-1: unreachable by any support/box cap; verified in
    # read-only Python during the 2026-08-07 review).
    for A in (_kron(3), _star([1, 2, 6])):
        f = as_unit_form(A)
        v = is_weakly_nonnegative(f)
        assert v.holds is False
        assert titsc(A, list(v.witness)) < 0                          # strict negative


@selfcert
def test_witness_is_exact_and_nonnegative():
    f = as_unit_form(_kron(3))
    v = is_weakly_nonnegative(f)
    assert all(x >= 0 for x in v.witness) and any(x > 0 for x in v.witness)
    assert v.witness_value == titsc(_kron(3), list(v.witness))


@selfcert
def test_weak_positivity_budget_is_honest_None():
    # a tiny budget must return None (undecided), never a silent True.
    f = as_unit_form(linear_path_algebra(5, field=CC))
    v = is_weakly_positive(f, budget=1)
    assert v.holds in (True, None)
    if v.holds is None:
        assert v.reason.startswith("budget")


@selfcert
def test_weak_nonneg_outside_coverage_is_None_not_guessed_true():
    # a triangular algebra whose vertex count exceeds HYPERCRITICAL_COVERAGE, with NO
    # witness found: is_weakly_nonnegative must return None (honest), never True.
    # (Pick an input at implementation once coverage is fixed; assert holds is not True
    #  unless coverage includes n. The load-bearing assertion: no guessed True.)
    ...
```

- [ ] **Step 2: Run to verify failure**, implement, iterate.
- [ ] **Step 3: Implement** the `UnitForm` connected-support **branch-and-bound** sweep
  (grow connected subsets over `underlying_graph`, dedup; per support search sincere
  `{1..6}^S` coordinate-by-coordinate with the completion lower-bound prune; early-return on
  the first `q ≤ 0`/`< 0`; count evaluations, `None` at `budget`), used directly for
  `is_weakly_positive` and as the witness-finder step (1) of `is_weakly_nonnegative`; then
  wire steps (2)–(4) to `_tits_lists` (Task 3). **Adjust to reality:** `_star([1,2,4])` = E8
  (8 vertices, weakly positive), `_star([1,2,5])` = Ẽ8 (9 vertices, isotropic radical
  `(6,3,4,2,5,4,3,2,1)`), `_star([1,2,6])` = T_{2,3,7} (10 vertices, wild, LIST-decided) —
  recheck vertex counts via `1/p+1/q+1/r` and the isotropic-root readout (all verified in
  read-only Python this revision: min q = 1, 0, −1 respectively). Orient arms outward from
  the center (irrelevant to the symmetric form). Fill `test_weak_nonneg_outside_coverage`
  once `HYPERCRITICAL_COVERAGE` is fixed in Task 3.
- [ ] **Step 4: Run** — `... -m pytest tests/invariants/test_weak_positivity.py -q` green.
- [ ] **Step 5: Commit** — `feat(invariants): weak positivity (Ovsienko box, branch-and-bound) + weak nonnegativity (classified hypercritical list, honest None outside coverage), witness on failure`

---

### Task 3: the classified hypercritical list — the PRIMARY weak-nonnegativity decision

**Files:**
- Create: `src/quiverlab/invariants/_tits_lists.py`
- Modify: `src/quiverlab/invariants/tits.py` (the `is_weakly_nonnegative` steps (2)–(4)
  consume this)
- Test: `tests/invariants/test_tits_lists.py`

**Interfaces:**
```python
# _tits_lists.py -- the DECISION data (hypercritical drives the tame axis) + a matcher
CRITICAL_FORMS: tuple        # {"name": "~E8", "gram": (...), "radical": (...)} -- the
                             #   Euclidean list ~A_n/~D_n/~E_{6,7,8}; TRANSCRIBABLE now
                             #   (standard extended-Dynkin marks). q(radical)=0. This list
                             #   is a CROSS-ORACLE for weak positivity (decided by Ovsienko).
HYPERCRITICAL_FORMS: tuple   # {"name": "T237", "gram": (...), "defect": (...)} -- von Hohne
                             #   1996 / de la Pena Banach-26 / BJP 2019. The PRIMARY decision
                             #   data for weak nonnegativity. q(defect) < 0. Encode every
                             #   entry transcription-checked against a primary source; the
                             #   minimal wild TREES (T_{2,3,7}, T_{2,4,5}, T_{3,3,4}, ...) are
                             #   transcribable now (defects verified in read-only Python).
HYPERCRITICAL_COVERAGE: frozenset | dict
                             #   the vertex counts / families for which HYPERCRITICAL_FORMS
                             #   is PROVABLY complete. is_weakly_nonnegative returns a
                             #   certified True only when the algebra's n is covered; else
                             #   honest None. NEVER a guessed True.

def matches_hypercritical(unit_form) -> "entry | None"
    # does `unit_form` (a restriction q|_S) equal, up to a permutation of variables, a
    # listed hypercritical form? small-n permutation/relabel search. Returns the matched
    # entry (carrying the tabulated defect vector, zero-extended to the caller's index
    # space) or None. This is step (2) of is_weakly_nonnegative.
def dominates_hypercritical(unit_form) -> "entry | None"
    # optional: q|_S "contains" a hypercritical form as a sub-restriction (a listed form
    # on a subset of S's variables). Folded into the connected-support enumeration.
```

**The decision (verbatim, the ruling).** `is_weakly_nonnegative(form)`: (1) run the bounded
branch-and-bound **witness-finder** — any `d ≥ 0` found with `q(d) < 0` is a **sound
False** with an exact witness (no bound assumption). (2) Else consult the encoded
hypercritical **LIST**: if some connected restriction `q|_S` matches (permutation-equivalent
to) a listed hypercritical form, return **False** with that entry's tabulated defect vector
(zero-extended). (3) Else, if the encoded list is **complete for every vertex count ≤ n**
(`HYPERCRITICAL_COVERAGE ⊇ {1..n}`), return **True** (certified — no hypercritical
restriction exists). (4) Else return **None**, reason `"hypercritical list partially
transcribed — no certified weakly-nonnegative verdict for n vertices"`. **Never a True from
a guessed box, EVER**; support/entry bounds are read off the list, never assumed.

**Honest scope.** Weak positivity is decided by Ovsienko's box (Task 2, cited);
`CRITICAL_FORMS` is only a *cross-oracle* there (every listed Euclidean form is
`is_weakly_positive.holds is False` with `q(radical) = 0`, and `is_weakly_nonnegative.holds
is True`). Weak nonnegativity is decided by `HYPERCRITICAL_FORMS` (above). Encode what is
transcription-checkable, record `HYPERCRITICAL_COVERAGE`, and state on the verification page
exactly which vertex counts / families are covered (and therefore where a certified `tame`
verdict is emittable vs where the verdict is honest-`None`).

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_tits_lists.py
"""The classified critical (Euclidean) / hypercritical lists. Critical: q(radical)=0,
not weakly positive, weakly nonnegative (cross-oracle). Hypercritical: q(defect)<0, drives
the weak-nonnegativity decision -- including T_{2,3,7} (10 vars, defect (12,6,8,4,10,9,7,6,
4,2), q=-1, verified in read-only Python 2026-08-07). Refs: von Hohne 1996; de la Pena
Banach-26 (1990); BJP 2019."""
import pytest
from quiverlab.invariants._tits_lists import (CRITICAL_FORMS, HYPERCRITICAL_COVERAGE,
                                              HYPERCRITICAL_FORMS)
from quiverlab.invariants.tits import UnitForm, is_weakly_nonnegative, is_weakly_positive

pytestmark = pytest.mark.oracle_literature


@pytest.mark.parametrize("entry", CRITICAL_FORMS, ids=[c["name"] for c in CRITICAL_FORMS])
def test_critical_radical_isotropic_and_not_weakly_positive(entry):
    f = UnitForm(len(entry["gram"]), entry["gram"], tuple(range(len(entry["gram"]))))
    assert f.evaluate(entry["radical"]) == 0                 # radical: q(z) = 0
    assert is_weakly_positive(f).holds is False              # critical => not weakly pos
    assert is_weakly_nonnegative(f).holds is True            # Euclidean => weakly nonneg


@pytest.mark.parametrize("entry", HYPERCRITICAL_FORMS,
                         ids=[c["name"] for c in HYPERCRITICAL_FORMS] or ["<deferred>"])
def test_hypercritical_defect_negative_and_decided_false(entry):
    f = UnitForm(len(entry["gram"]), entry["gram"], tuple(range(len(entry["gram"]))))
    assert f.evaluate(entry["defect"]) < 0                   # q(defect) < 0
    assert is_weakly_nonnegative(f).holds is False           # list-decided False + witness


def test_T237_is_a_10_variable_hypercritical_entry():
    # THE refutation, encoded so no one reintroduces a <=9 cap: T_{2,3,7} is hypercritical
    # with 10 variables and a SINCERE defect vector.
    names = {c["name"] for c in HYPERCRITICAL_FORMS}
    assert "T237" in names
    e = next(c for c in HYPERCRITICAL_FORMS if c["name"] == "T237")
    assert len(e["gram"]) == 10 and all(x > 0 for x in e["defect"])   # sincere, 10 vars


def test_coverage_is_documented():
    # honest scope: coverage is an explicit datum the verification page cites.
    assert HYPERCRITICAL_COVERAGE is not None
```

- [ ] **Step 2: Run to verify failure**, encode `CRITICAL_FORMS` (extended-Dynkin marks;
  verify each `q(radical)=0`) and `HYPERCRITICAL_FORMS` (at minimum the minimal wild trees
  `T_{2,3,7}`/`T_{2,4,5}`/`T_{3,3,4}` with their read-only-verified defects, plus every
  non-tree entry transcription-checked against von Höhne 1996 / BJP 2019); set
  `HYPERCRITICAL_COVERAGE` to exactly the provably-complete range and `# PIN` the rest.
- [ ] **Step 3: Run** — `... -m pytest tests/invariants/test_tits_lists.py -q` green.
- [ ] **Step 4: Commit** — `feat(invariants): classified hypercritical list = primary weak-nonnegativity decision (T237 10-var refutation encoded); Euclidean critical list as weak-positivity cross-oracle`

---

### Task 4: the verdict layer — `tame_wild_certificate` (P56-gated, CC-only)

**Files:**
- Modify: `src/quiverlab/invariants/tits.py`, `src/quiverlab/core/algebra.py`
  (`Algebra.tame_wild_certificate`, and the thin `Algebra.tits_form_combinatorial` /
  `Algebra.is_weakly_positive` / `Algebra.is_weakly_nonnegative` delegates)
- Test: `tests/invariants/test_tame_wild.py`

**Interfaces:**
```python
@dataclass(frozen=True)
class TameWildCertificate:
    # -- form layer (field-free; always present for a triangular presented algebra) --
    gram: tuple                       # the integer Gram matrix (display)
    minimal_relation_counts: dict     # r_ij
    is_unit_form: bool
    weakly_positive: bool | None      # FormVerdict.holds (None = budget)
    weakly_nonnegative: bool | None
    witness: tuple | None             # the exact d >= 0 with q(d) <= 0 / < 0
    witness_value: int | None
    # -- P56 certificate trail --
    simply_connected: bool | None     # is_simply_connected(A).verdict
    strongly_simply_connected: bool | None   # is_strongly_simply_connected(A).verdict
    strong_certificate: object | None # the P56 StrongSimpleConnectivity (witness/reason/
                                       #   checked_convex); None if not computable
    # -- verdict layer (CC + gate) --
    field_alg_closed: bool
    rep_type: str | None              # "rep-finite" | "tame" | "wild" | None (refused/
                                       #   inconclusive); None NEVER means a guessed verdict
    reason: str                       # the full trail: which gate fired / why None
    scope_note: str                   # the honest scope sentence for the payload/report

def tame_wild_certificate(A, *, convex_budget=20000, search_budget=3_000_000) -> TameWildCertificate
```

**Decision logic (SOUND; a verdict NEVER outruns its hypotheses):**
1. **Guard:** `A.quiver is None` → loud refusal (needs the quiver).
2. **Form layer (always):** `f = as_unit_form(A)` (raises unless `A.quiver.is_acyclic()` —
   the Tits-form tame/wild theory needs a triangular quiver; M2 gate).
   `weakly_positive = is_weakly_positive(f, budget=search_budget)`,
   `weakly_nonnegative = is_weakly_nonnegative(f, budget=search_budget)`; carry the witness
   (prefer the weak-nonnegativity witness, the sharper one: `q < 0` ⇒ wild direction).
3. **P56 trail:** `sc = A.is_simply_connected()`; `ssc = A.is_strongly_simply_connected()`
   (P56 delegates). Record `sc.verdict`, `ssc.verdict`, and the `ssc.strongly` (or the
   `StrongSimpleConnectivity`) object.
4. **Field gate:** `field_alg_closed = _is_alg_closed(A.field)` (CC only; see below).
5. **Verdict assembly (the trichotomy, gated):**
   - If `not field_alg_closed`: `rep_type = None`, reason names the field
     ("verdict is algebraically-closed-only (CC); the form above is field-free and
     reported"). **Stop** (form still returned).
   - **tame/wild axis (BdlPS, needs STRONG s.c.):** if `ssc.verdict is True`:
     - `weakly_positive is True` → `rep_type = "rep-finite"` (⊂ tame).
     - `weakly_nonnegative is False` → `rep_type = "wild"` (a witness was FOUND — by the
       branch-and-bound finder or a list match — so this is ALWAYS sound, no
       list-completeness needed; the `q < 0` witness IS the wild certificate).
     - `weakly_nonnegative is True and weakly_positive is False` → `rep_type = "tame"`
       (properly tame; the isotropic witness records the tame direction). **Emittable ONLY
       because the hypercritical list is complete for `n` — see the honesty rule below.**
     - **else** (`weakly_nonnegative is None` — no witness found and the list is incomplete
       for `n`, OR `weakly_positive is None` on budget): `rep_type = None`, reason the
       list-incompleteness (`"hypercritical list partially transcribed — no certified tame
       verdict for n vertices"`) or the budget — the honest middle. **Never a guessed
       tame.**
   - **rep-finite-only fallback (Bongartz, needs SIMPLE s.c.):** elif `sc.verdict is True`:
     Bongartz still decides the rep-finite axis. `weakly_positive is True` →
     `rep_type = "rep-finite"`; `weakly_positive is False` → `rep_type` records
     `"not rep-finite"` (rep-infinite) **but not tame-vs-wild** (BdlPS needs *strong* s.c.,
     which is `ssc.verdict is not True` here) — reason states the tame/wild split is
     withheld pending the strong certificate. (Encode this honest middle state as
     `rep_type = None` with a precise reason, OR add an explicit `"rep-infinite"` value —
     pick `None` + reason to avoid a third label that the GUI must translate; **decide at
     implementation and pin it in the test**.)
   - **else** (`ssc.verdict` and `sc.verdict` both not `True` — `False` or `None`):
     `rep_type = None`. Reason distinguishes: `ssc.verdict is False` →
     "not strongly simply connected (P56 witness: …) — the Tits-form verdict is out of
     scope"; `ssc.verdict is None` → propagate P56's reason verbatim
     (Adian–Rabin / budget_exceeded / undecided_char). **The None of P56 is the None
     here** — never upgraded to a verdict.
6. `scope_note`: the standing honest sentence (form field-free; verdict CC + strongly
   simply connected; Bongartz/BdlPS attributions).

**Tame-emission honesty rule (VERBATIM — the ruling).** `rep_type = "tame"` is emitted
ONLY when `is_weakly_nonnegative` returned a **certified** `True` resting on a hypercritical
list that is **complete for the algebra's vertex count** (`HYPERCRITICAL_COVERAGE ⊇
{1..n}`), AND `weakly_positive is False`, AND `ssc.verdict is True`. While any relevant
hypercritical entry is untranscribed (`# PIN`), `is_weakly_nonnegative` returns `None` with
reason `"hypercritical list partially transcribed — no certified tame verdict"`, and
`tame_wild_certificate` sets `rep_type = None` (verdict withheld). By contrast `"wild"`
needs only a FOUND witness (branch-and-bound or list match — always sound) and `"rep-finite"`
needs only Ovsienko's box `True` (a cited complete decision). **No tame/wild verdict is ever
derived from a guessed box bound.**

**`_is_alg_closed(field)`** — **confirmed by grep this revision:** quiverlab's only
algebraically closed field is `CC` (`fields/complexfield.py::ComplexField`, `CC =
ComplexField()`); there is NO `is_algebraically_closed` attribute. Test
`isinstance(field, ComplexField)` (import lazily) — `GF`, `QQ`, `QQi`, `E` all return
`False`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_tame_wild.py
"""The Tits-form tame/wild verdict, gated on the P56 strong-simple-connectivity
certificate over CC. rep-finite <=> weakly positive (Bongartz, simply connected);
tame <=> weakly nonnegative (BdlPS, strongly simply connected). Off scope: verdict None,
form still computed. Refs: Bongartz Math. Ann. 269 (1984); BdlPS Adv. Math. 226 (2011);
Kasjan-Skowronski arXiv:1905.06028."""
import pytest
from quiverlab import CC, GF, Quiver, linear_path_algebra
from quiverlab.invariants.tits import tame_wild_certificate

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _star(arms):   # (same helper as test_weak_positivity)
    verts, ars, nxt = [0], {}, 1
    for ai, length in enumerate(arms):
        prev = 0
        for _ in range(length):
            verts.append(nxt); ars[f"e{ai}_{nxt}"] = (prev, nxt); prev = nxt; nxt += 1
    return Quiver(verts, ars).algebra(field=CC)


@lit
def test_dynkin_tree_rep_finite():
    c = tame_wild_certificate(linear_path_algebra(5, field=CC))   # A5 tree: ssc
    assert c.strongly_simply_connected is True
    assert c.weakly_positive is True and c.rep_type == "rep-finite"


@lit
def test_euclidean_tree_tame():
    c = tame_wild_certificate(_star([1, 2, 5]))                   # ~E8 tree: ssc, tame
    assert c.strongly_simply_connected is True
    assert c.weakly_nonnegative is True and c.weakly_positive is False
    assert c.rep_type == "tame" and c.witness is not None         # isotropic direction


@lit
def test_T237_tree_wild():
    c = tame_wild_certificate(_star([1, 2, 6]))                   # T_{2,3,7}: ssc, wild
    assert c.strongly_simply_connected is True
    assert c.weakly_nonnegative is False and c.rep_type == "wild"
    assert c.witness_value < 0                                    # the wild certificate


@selfcert
def test_kronecker_form_layer_only_not_simply_connected():
    # THE honest layer split: the m-Kronecker is NOT simply connected (parallel arrows,
    # pi1^ab = Z for m>=2), so the VERDICT is refused -- the m-Kronecker ladder exercises
    # the FORM layer, not the algebra verdict. The form booleans are still computed.
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=CC)   # 2-Kronecker
    c = tame_wild_certificate(K)
    assert c.weakly_nonnegative is True and c.weakly_positive is False # form: ~A1
    assert c.simply_connected is False                                 # not simply connected
    assert c.rep_type is None and "simply connected" in c.reason


@selfcert
def test_non_CC_refuses_verdict_keeps_form():
    A = linear_path_algebra(5, field=GF(7))                      # not algebraically closed
    c = tame_wild_certificate(A)
    assert c.weakly_positive is True                             # form computed (field-free)
    assert c.field_alg_closed is False and c.rep_type is None
    assert "algebraically closed" in c.reason or "CC" in c.reason


@selfcert
def test_P56_None_propagates():
    # if P56's strong verdict is None (budget / char), rep_type is None with P56's reason
    # -- never a fabricated tame/wild. (Construct/So mock a P56 None at implementation:
    # use an input where is_strongly_simply_connected returns None, e.g. force a tiny
    # convex_budget so P56 reports budget_exceeded.)
    A = linear_path_algebra(5, field=CC)
    c = tame_wild_certificate(A, convex_budget=0)
    assert c.rep_type in (None, "rep-finite")     # None if strong went budget; else Bongartz
    if c.rep_type is None:
        assert "budget" in c.reason.lower() or "simply connected" in c.reason
```

- [ ] **Step 2: Run to verify failure**, implement, iterate. **Decide-and-pin** the
  rep-infinite middle state (Bongartz says not-rep-finite but BdlPS withheld): the
  recommended choice is `rep_type = None` + an explicit reason; encode it and make
  `test_P56_None_propagates` / a dedicated `test_bongartz_rep_infinite_middle` assert the
  chosen contract.
- [ ] **Step 3: Implement** the assembler with the `_is_alg_closed(CC)` gate and the
  P56-`None` propagation. **Adjust to reality:** `linear_path_algebra(5, field=CC)` is a
  tree ⇒ P56's `is_strongly_simply_connected` returns `True` via the R1 tree route
  (cheap); confirm `convex_budget=0` actually forces a `None` (P56 maps budget to `None`)
  — if P56's tree route decides before consulting the budget, pick a NON-tree ssc input
  for the None test (a strongly-s.c. quiver whose convex sweep is nontrivial).
- [ ] **Step 4: Run** — `... -m pytest tests/invariants/test_tame_wild.py -q` green.
- [ ] **Step 5: Commit** — `feat(invariants): Tits-form tame/wild certificate -- rep-finite/tame/wild trichotomy gated on the P56 strong-s.c. certificate over CC, honest None propagation`

---

### Task 5: consistency oracles — the ladder, the trees, and the knit cross-check

**Files:**
- Test: `tests/invariants/test_tame_wild_oracles.py`

**Interfaces:** consumes only the public surface (Tasks 1–4) + P41 `A.ar_quiver()` (the
knit) + P56 verdicts. No new `src/`. This is the **ground-truth agreement** battery.

**The layer split (the record's key correction, made explicit in tests).**
- **Form layer (field-free), the m-Kronecker ladder:** `m = 1` weakly positive; `m = 2`
  weakly nonnegative not weakly positive; `m ≥ 3` not weakly nonnegative. **The Kronecker
  quiver with `m ≥ 2` is NOT simply connected** (Betti number `m − 1 ≥ 1`; parallel arrows
  ⇒ `π₁^{ab} = ℤ^{m−1} ≠ 0`), so it is **outside the verdict scope**: the ladder tests
  `is_weakly_positive`/`is_weakly_nonnegative` directly, and asserts `rep_type is None`
  (with `simply_connected is False`) for `m ≥ 2`. Only `m = 1` (= `A_2`, a tree) reaches
  the verdict layer and comes out rep-finite. This is the honest resolution of "which
  oracle exercises which layer".
- **Verdict layer (CC), the trees:** Dynkin trees (`A_n`, `D_n`, `E_{6,7,8}`) → rep-finite;
  Euclidean trees (`D̃_n`, `Ẽ_{6,7,8}` as star trees, simply connected ✓) → tame;
  `T_{2,3,7}` (star arms `1,2,6`, simply connected) → **wild, decided by the list entry**
  (its defect vector is sincere on all 10 vertices — no box/support cap reaches it).

**Expected cost + bucket per oracle (M1 residual, after branch-and-bound).** All in the
**fast** bucket: the m-Kronecker ladder is `n = 2` (trivial); Dynkin trees `A_5`/`E_8`
weak-positivity sweeps prune hard (positive-definite ⇒ the completion lower bound rises
early — far below the `6^8 = 1.68M` worst case); `Ẽ_8` (`n = 9`) returns `False` on the
center-first witness before nearing `6^9 ≈ 10M`; `T_{2,3,7}` and the `m ≥ 3` Kroneckers are
**list/finder-decided** (no full sweep); the tree verdicts add P56's convex sweep on ≤ 10
vertices (tree route = R1, `O(1)`). Each is well under a second; keep the star arms small
so P56 stays cheap.
- **Knit cross-check (cross-engine):** wherever P41's `ar_quiver()` **completes**
  (`is_complete` — rep-finite certified) on a strongly-s.c. CC instance, the verdict MUST
  be `"rep-finite"` and `weakly_positive is True`. A strong two-engine oracle: the AR knit
  (module-category enumeration) and the Tits form (quadratic-form arithmetic) are wholly
  independent implementations of "representation-finite". **Documented converse (not a
  test):** a `weakly_positive`-but-`ar_quiver` **budget-exhausted** instance is a
  *prediction* the knit can confirm with a raised budget — recorded in the plan/verification
  page, not asserted (the knit's incompleteness there is honest, not a disagreement).
- **P55/P56 shared fixture (end-to-end pin):** the `rad²=0` `A_5` fixture P55/P56 share —
  run `is_strongly_simply_connected` live; **if** it is strongly s.c. its rep-finite
  verdict must come out weakly positive (a genuine end-to-end pin joining P56 → P62).
  If P56 reports it NOT strongly s.c., assert `rep_type is None` with the ssc-False reason
  (equally a valid pin — the honesty is the point).

- [ ] **Step 1: Write the tests**

```python
# tests/invariants/test_tame_wild_oracles.py
"""Ground-truth agreement (Plan 62 / R19). m-Kronecker ladder = FORM layer (Kronecker is
NOT simply connected); Dynkin/Euclidean trees + T_{2,3,7} = VERDICT layer; the AR-knit
cross-check ties rep-finite to weakly positive on ssc CC instances. Refs: Bongartz 1984;
BdlPS 2011; Gabriel; Nazarova."""
import pytest
from quiverlab import CC, Quiver, linear_path_algebra
from quiverlab.invariants.tits import (is_weakly_nonnegative, is_weakly_positive,
                                       as_unit_form, tame_wild_certificate)

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _kron(m):
    return Quiver([1, 2], {f"a{k}": (1, 2) for k in range(m)}).algebra(field=CC)


def _star(arms):
    verts, ars, nxt = [0], {}, 1
    for ai, length in enumerate(arms):
        prev = 0
        for _ in range(length):
            verts.append(nxt); ars[f"e{ai}_{nxt}"] = (prev, nxt); prev = nxt; nxt += 1
    return Quiver(verts, ars).algebra(field=CC)


@lit
@pytest.mark.parametrize("m,wp,wnn", [(1, True, True), (2, False, True), (3, False, False)])
def test_kronecker_ladder_is_form_layer_only(m, wp, wnn):
    K = _kron(m)
    f = as_unit_form(K)
    assert is_weakly_positive(f).holds is wp
    assert is_weakly_nonnegative(f).holds is wnn
    c = tame_wild_certificate(K)
    if m == 1:
        assert c.rep_type == "rep-finite"                 # A2 tree: reaches the verdict
    else:
        assert c.simply_connected is False and c.rep_type is None   # not in verdict scope


@lit
@pytest.mark.parametrize("A,expect", [
    (linear_path_algebra(6, field=CC), "rep-finite"),     # A6
    (_star([1, 1, 3]), "rep-finite"),                     # D6 tree
    (_star([1, 2, 2]), "rep-finite"),                     # E6
    (_star([1, 2, 5]), "tame"),                           # ~E8
    (_star([1, 2, 6]), "wild"),                            # T_{2,3,7}
])
def test_tree_verdicts(A, expect):
    assert tame_wild_certificate(A).rep_type == expect


@xeng
def test_knit_agrees_with_weakly_positive():
    for A in (linear_path_algebra(5, field=CC), _star([1, 2, 2])):   # A5, E6: rep-finite
        c = tame_wild_certificate(A)
        ar = A.ar_quiver()
        if ar.is_complete:                                 # knit closed => rep-finite
            assert c.rep_type == "rep-finite" and c.weakly_positive is True
```

- [ ] **Step 2: Run**; recheck every tree encoding (`_star` arm lengths ↔ Dynkin/Euclidean
  type via `1/p+1/q+1/r`) and each verdict against the classical tables. A FAILURE means a
  wrong quiver encoding or a P56/P41 disagreement to escalate — fix the QUIVER, never the
  ground-truth verdict.
- [ ] **Step 3: Commit** — `test(invariants): tame/wild ground-truth oracles -- m-Kronecker form layer, Dynkin/Euclidean/T237 tree verdicts, AR-knit cross-check`

---

### Task 6: QPA honest-scope skip

**Files:**
- Create: `tests/qpa/test_tame_wild_qpa.py`

**Interfaces:** QPA 1.37 has **no** Tits-form tame/wild surface (no `TitsForm`,
`WeaklyPositive`, `WeaklyNonNegative`, `RepresentationType` verb). Probe live via
`NamesGVars()` (the `tests/qpa/test_products_qpa.py` / P56 precedent): the test SKIPS
honestly and **FAILS if QPA ever ships** any such verb, so the honest-scope claim on the
verification page cannot silently rot.

- [ ] **Step 1: Write the probe** (`pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
  reason=...)`; sweep `NamesGVars()` for the candidate names; `pytest.skip(...)` if absent,
  `assert False` (loud) if any appears).
- [ ] **Step 2: Run** — `... -m pytest tests/qpa/test_tame_wild_qpa.py -q -m qpa` (venv has
  `[qpa]`): expect a clean skip.
- [ ] **Step 3: Commit** — `test(qpa): honest no-Tits-tame/wild-surface probe (fails if QPA ships one)`

---

### Task 7: no-code exposure — the `tame_wild` compute kind

The organizing invariant is **byte-identity across the two runners**, achieved by routing
BOTH through a single **shared library builder** (the `recognizers_block` /
`coverings_block` precedent — each runner only adds `citations` locally). Anchors below are
current-tree; still "adjust to reality" if a merge shifts them.

**Files (the seven-touchpoint checklist):**
1. **Shared builder** — Create `src/quiverlab/invariants/tits_block.py`:
   `tame_wild_block(A) -> dict` returning the payload with a `"references"` list (both
   runners import it, so the block is byte-identical by construction).
2. **Server/HPC** — Modify `src/quiverlab/hpc/spec.py`: one `if kind == "tame_wild":`
   branch in `_dispatch` (the `strings`/`recognizers` pattern) calling the builder +
   `block["citations"] = _citation_pairs(block["references"])`; one `_snip` reproduce entry
   (`"tame_wild": lambda it: "A.tame_wild_certificate()"`). NOT in `MODULE_KINDS` (it is an
   algebra-level scalar routed through `_dispatch`).
3. **Pyodide twin** — Modify `docs/gui/runner.py`: one byte-identical `elif name ==
   "tame_wild":` branch in `compute_one`; one entry in the `calls` reproduce map
   (hard-indexed — KeyError if missing); one `ETA_MODEL["scalars"]` cost (`tame_wild`:
   ~3.0 — the bounded search + P56 convex sweep dominate; size it above the cheap scalars).
4. **GUI JS** — Modify `docs/gui/gui.js` AND `webapp/static/gui/gui.js` (diff-identical —
   apply every edit to both): checkbox `<input id="qlgui-tame_wild">` in the invariants
   block; the bare id in the id-registry array; the scalar push-list entry; `KIND_CTRL`
   entry `{cb: "tame_wild"}`; a `renderBlock` `else if (name === "tame_wild")` branch
   (render the trichotomy verdict, the two form booleans, the witness vector, the Gram
   matrix via the existing `matrixGrid`, the P56 certificate line, and the scope note)
   before the citations mount; the `scheduleProbe` change-listener entry; the `THEMES`
   `kinds` array (invariants theme) inside the `QLGUI-THEMES-BEGIN/END` sentinels.
   Optional: a structured renderer in `webapp/static/app.js`.
5. **i18n (FOUR locales) + picker template** — `webapp/server/i18n/` has
   `LANGS = ("en","es","fr","zh")`. Add to **all four** `{en,es,fr,zh}.json`:
   `pick.kind.tame_wild` (REQUIRED — gated by `test_layout_picker.py`), `inv.tame_wild`
   (checkbox label), `block.tame_wild.title` + per-row labels
   (`form`, `weakly_positive`, `weakly_nonnegative`, `witness`, `verdict_rep_finite`,
   `verdict_tame`, `verdict_wild`, `verdict_undecided`, `certificate`, `scope_note`). Add
   one `data-pick-kind-tame_wild="{{ t('pick.kind.tame_wild') }}"` line in
   `webapp/templates/draw.html`.
6. **Report renderer** — Modify `src/quiverlab/trace/results_html.py`: one `if kind ==
   "tame_wild":` branch in `_block_html` (render the form as a `matrix_grid`, the witness,
   the trichotomy, and the honest scope sentence + the P56 certificate trail) + one
   `_HEADINGS` entry. **Zero-safe / honest per the report passes:** the verdict `None`
   renders as "verdict withheld — <reason>", never a blank or a guessed type.
7. **Tests / gates** — add `tame_wild` to `ALL_KINDS` in `tests/webapp/test_layout_picker.py`
   (must equal the THEMES kind set, each once — `test_themes_cover_every_kind_exactly_once`
   + the four-locale `pick.kind.*` gate); add ONE golden to `tests/webapp/_runner_goldens.json`
   with a dated change-log bullet in `tests/webapp/test_runner_delegation.py` (verify
   existing entries byte-identical FIRST); add a twin-parity test (copy
   `test_quasi_hereditary_p47.py::test_twin_parity` / P56's `test_coverings_exposure`).
- Test files: `tests/webapp/test_tame_wild_exposure_p62.py`, `tests/hpc/test_p62_kinds.py`.

**Schema stays v1** (scalar kind on the existing algebra block). **Canonical keys are
request-derived** (`cache.py::canonical_key`, sha256 over the sorted-keys blob +
`library_version()`): a scalar-only request that does not use the new key keys
byte-identically to before, so existing goldens are untouched.

**The `tame_wild` block:**
```json
{"gram": [[...]], "minimal_relation_counts": {"1,3": 1},
 "is_unit_form": true, "weakly_positive": true|false|null,
 "weakly_nonnegative": true|false|null, "witness": [ints]|null, "witness_value": int|null,
 "simply_connected": true|false|null, "strongly_simply_connected": true|false|null,
 "rep_type": "rep-finite"|"tame"|"wild"|null, "reason": str, "scope_note": str,
 "latex": "q_A(x)=...", "references": ["bongartz_criterion","bdps_tame_tits",
 "kasjan_skowronski","ovsienko_forms","vonhohne_wnn","bjp_quadratic_forms"], "citations": ...}
```
`rep_type = null` (verdict withheld) renders as **"representation type undecided — <reason>"**
(the reason names: non-CC field / not (strongly) simply connected + P56 witness /
Adian–Rabin `None` / search budget) in GUI + report + all four locales. A refusal on this
input (presentation-less, or a non-unit form / loop) is reported as
`{"error": "<the loud message>"}`, never a silent default (the Plan-30 honest-per-entry
precedent).

- [ ] **Step 1: Write the failing cross-runner test** (copy the P56 / m0729 runner-pair
  fixture; unmarked file):

```python
# tests/webapp/test_tame_wild_exposure_p62.py
"""tame_wild kind: served by hpc.spec, mirrored by the Pyodide twin, verdict rendered
honestly, both runners byte-identical."""
import json


def test_tame_wild_block_shape(tmp_path):
    # request: A5 tree over CC, compute=["tame_wild"]:
    #   block["rep_type"] == "rep-finite"; block["weakly_positive"] is True
    #   "bongartz_criterion" in block["references"]; block["is_unit_form"] is True
    ...


def test_tame_wild_wild_witness(tmp_path):
    # T_{2,3,7} star over CC: rep_type "wild", witness present, witness_value < 0.
    ...


def test_tame_wild_off_scope_keeps_form(tmp_path):
    # 2-Kronecker over CC: weakly_nonnegative True, rep_type null, reason names simple conn.
    ...


def test_twin_parity(tmp_path):
    # run all requests through docs/gui/runner.py::compute_one and hpc.spec.run;
    # assert json.dumps(block, sort_keys=True) equality.
    ...
```

- [ ] **Step 2: Implement** the shared builder, then the server + twin dispatch, the GUI JS
  in both `gui.js` files, the four-locale i18n + `draw.html`, and the `results_html.py`
  branch.
- [ ] **Step 3: Gate the layout/i18n first** —
  `... -m pytest tests/webapp/test_layout_picker.py tests/webapp/test_js_parses.py -q`
  (THEMES == `ALL_KINDS`; `pick.kind.tame_wild` in all four locales).
- [ ] **Step 4: Add the golden** (`tame_wild_a5_cc`) to `_runner_goldens.json`; note it in
  the `test_runner_delegation.py` docstring change-log. Run the delegation test BEFORE
  adding to confirm existing entries stay byte-identical; on a rebase, merge the JSON
  **semantically** (dict-union), never textually.
- [ ] **Step 5: Run the gates** —
  `... -m pytest tests/webapp/test_tame_wild_exposure_p62.py tests/webapp/test_runner_delegation.py tests/hpc -q`.
  Expected PASS; both runners byte-identical.
- [ ] **Step 6: Commit** — `feat(gui,webapp,hpc): tame_wild compute kind -- Tits-form tame/wild verdict + witness clickable end-to-end (4 locales, v1 schema)`

---

### Task 8: verification page, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Test: existing release gates

- [ ] **Step 1: Verification page.** Add the Plan-62 subsystem row (`invariants/tits.py`
  | oracles: `oracle_literature` — the m-Kronecker ladder (form layer), the Dynkin trees
  rep-finite, the Euclidean trees tame, `T_{2,3,7}` wild, Bongartz's criterion, the BdlPS
  theorem; `oracle_crossengine` — `r_ij` (P56 count) == `dim Ext²(S_i,S_j)`, and the
  AR-knit == weakly-positive agreement; `oracle_selfcert` — the witness is exact and `≥ 0`
  with `q(witness) ≤ 0`, the budget→`None` honesty, the box-search vs printed-list
  agreement; `qpa` — the honest no-surface probe). Add the **honest-scope entries**:
  (a) the **verdict layer is CC-only** (algebraically closed; the tame/wild dichotomy is
  an algebraically-closed notion — Drozd); over other fields the form is computed, the
  verdict refused. (b) the **verdict is gated on the P56 strong-simple-connectivity
  certificate** (tame/wild axis, BdlPS) resp. **simple connectivity** (rep-finite axis,
  Bongartz); P56's three-valued `None` (Adian–Rabin / budget / undecided_char) **propagates**
  to a `None` verdict — never a fabricated type. (c) the **m-Kronecker ladder exercises the
  FORM layer only** — the Kronecker quiver is not simply connected (`m ≥ 2`), so it is out
  of the verdict scope; this is stated so no reader mistakes it for a verdict oracle.
  (d) **weak nonnegativity (the tame axis) is decided by the classified hypercritical
  list**, NOT a guessed box — a certified `tame` verdict is emitted only where
  `HYPERCRITICAL_COVERAGE` is complete for the vertex count `n`, otherwise the verdict is
  honest `None` ("hypercritical list partially transcribed"); the `≤ 9`-variable cap was
  REFUTED by `T_{2,3,7}` (a 10-variable hypercritical form) and no universal support/entry
  cap is claimed. State exactly which hypercritical vertex counts / families are
  transcription-checked (covered) and which are deferred (`# PIN`). (e) QPA has no
  Tits-form tame/wild surface. Recount the class table (`tests/release/test_oracle_classes.py`
  drives the numbers — collect, paste, re-run to green).
- [ ] **Step 2: README.** One features line: "representation-type certificates — the
  combinatorial Tits form `q_A`, exact weak positivity/nonnegativity, and the rep-finite /
  tame / wild verdict for strongly simply connected algebras (Bongartz; Brüstle–de la
  Peña–Skowroński) — clickable via `tame_wild`."
- [ ] **Step 3: Full gate** — `... -m pytest -q -m fast` green;
  `... -m pytest tests/invariants -q` green; `... -m pytest tests/qpa -q -m qpa` green;
  `... -m pytest tests/release -q` green.
- [ ] **Step 4: Commit** — `docs(verification): Plan-62 oracle rows + honest scope (CC-only verdict, P56 gate + None propagation, m-Kronecker form-layer-only, hypercritical-list coverage) + recounted classes`

---

## Acceptance (Plan-62 definition of done)

1. `A.tits_form_combinatorial(d)` (W1 — NOT `tits_form`, which is P38's untouched
   homological Euler form) and the combinatorial Gram matrix are public, exact, field-free,
   with `r_ij` = P56's minimal-relation count (`invariants.coverings.minimal_relation_counts`,
   a single retained-dict implementation shared with `fundamental_group`) **and**
   cross-checked against `dim Ext²(S_i, S_j)`. Non-triangular (loop/oriented-cycle,
   `is_acyclic` gate — M2) and presentation-less algebras refuse loudly. Where both are
   defined (`gl.dim ≤ 2`) the combinatorial form equals P38's Euler symmetrization.
2. `A.is_weakly_positive()` decides EXACTLY by Ovsienko's box `[0,6]` (cited complete
   decision, branch-and-bound sweep), `None` only on budget. `A.is_weakly_nonnegative()`
   decides by the **classified hypercritical list** (primary) + a sound branch-and-bound
   witness-finder: a `False` is a FOUND exact witness (`q < 0`); a `True` is emitted **only
   where the list is complete for `n`** (`HYPERCRITICAL_COVERAGE`), else honest `None`
   ("hypercritical list partially transcribed"). **No guessed box, no `_WNN_BOX` decision,
   NO True from a guess.** The `T_{2,3,7}` 10-variable hypercritical entry is encoded (the
   refutation of the old `≤ 9` cap).
3. `A.tame_wild_certificate()` returns the trichotomy **rep-finite / tame / wild** gated
   on the P56 certificate over CC: rep-finite ⟺ weakly positive (Bongartz, simple
   connectivity, cited-complete); **wild** = a FOUND non-weak-nonnegativity witness (always
   sound); **tame** = certified weakly nonnegative (list complete for `n`) and not weakly
   positive (BdlPS, strong simple connectivity). Off scope (non-CC / not (strongly) simply
   connected / P56 `None` / list incomplete / budget) the verdict is `None` with the honest
   reason and the form is still computed. **P56's `None` is P62's `None` — never upgraded;
   `tame` is never guessed.**
4. Ground-truth oracles green: the m-Kronecker ladder (`1/2/≥3`) on the **form layer**
   (asserting `rep_type is None`, not simply connected, for `m ≥ 2`); Dynkin trees
   rep-finite, Euclidean trees tame, `T_{2,3,7}` wild on the **verdict layer**; the AR-knit
   `is_complete` ⇒ `rep-finite` + weakly-positive cross-engine agreement; the P55/P56 shared
   `rad²=0 A_5` fixture end-to-end.
5. `tame_wild` clickable end-to-end (GUI canvas → block → report) in all four locales
   (en/es/fr/zh), both runners byte-identical via the shared builder, one golden added with
   a documented change-log entry, `test_layout_picker.py` `ALL_KINDS`/THEMES gate green,
   schema still v1, canonical keys unchanged.
6. QPA honest-scope probe skips (fails if QPA ever ships a Tits-form tame/wild surface);
   `docs/verification.md` recounted with the five honest-scope entries; README line added;
   fast + invariants + qpa + release suites green.

---

## Methodology & assumptions

**Approach.** I read the committed P56 plan end-to-end (its `StrongSimpleConnectivity`
dataclass — `verdict: bool|None`, `witness`, `reason`, `checked_convex` — and the
`is_simply_connected`/`is_strongly_simply_connected` surface are the gate this plan
consumes; the corrected transitive-closure separation convention is P56's, not mine), P38's
shipped `invariants/forms.py` (to avoid conflating the homological Euler form with the
combinatorial Tits form), and P49 (bounded-search + honest-budget patterns). I re-verified
every R19 reference on the web: the Tits-form definition and `r_ij = dim Ext²` (Encyclopedia
of Mathematics + Kasjan–Skowroński); Bongartz Math. Ann. 269 (1984) 1–12 (rep-finite ⟺
weakly positive, simply connected); BdlPS Adv. Math. 226 (2011) 887–951 (tame ⟺ weakly
nonnegative, strongly simply connected); Ovsienko's ≤ 6 bound; von Höhne's hypercritical
classification (Comment. Math. Helv. 63 (1988) + Proc. LMS 73 (1996)); and the BJP 2019
book (Algebra and Applications 25). I confirmed code-level facts against the tree:
`CC = ComplexField()` is the only algebraically closed field (no `is_algebraically_closed`
attribute); `Algebra.simple(v, side=)`, `Algebra.ext(M,N,n)`, `Quiver.is_acyclic()`,
`Algebra.loewy_length()` exist; P38's `forms.py` uses `C^{-1}` (unimodular Cartan) and
definiteness, both distinct from what P62 needs. During the 2026-08-07 adversarial-review
revision I ran read-only Python to classify the star trees: `E_8` `min q = 1`, `Ẽ_8`
`min q = 0` (radical max entry 6), `T_{2,3,7}` `min q = −1` at the sincere 10-vertex defect
`(12,6,8,4,10,9,7,6,4,2)`.

**Key correctness decisions.** (1) The plan builds a *new* combinatorial Tits form
(`tits_form_combinatorial`, W1), NOT a reuse of P38's `form_type`/`tits_form` — because weak
positivity/nonnegativity (positive cone) ≠ definiteness (all of ℝⁿ), and the Tits form
truncates at Ext² whereas the Euler form does not (they agree only at `gl.dim ≤ 2`). (2)
**Weak positivity** rests on Ovsienko's cited box-6 decision (branch-and-bound sweep) — the
earlier inline "∎" proof was WITHDRAWN in review (its "critical ⇒ Euclidean isotropic
radical" step is false: the `m ≥ 3` Kronecker is critical yet indefinite with no isotropic
radical). **Weak nonnegativity** is decided by the **classified hypercritical list**, not a
guessed box: the `≤ 9`-variable transcription was refuted by this plan's own `T_{2,3,7}`
(a 10-variable hypercritical form with a sincere large-entry defect), so there is no safe
universal cap; a certified `tame` is emitted only under recorded list-completeness, else
honest `None`. A `False` (wild) is always a sound found witness. (3) The unit-form gate is
`is_acyclic()` (M2) — the `G_ii == 2` heuristic wrongly passes `k[x]/(x³)`. (4) The verdict
is layered exactly as the two theorems require and P56's `None` propagates — the plan never
lets a search, a gate, or a partial list manufacture a verdict.

**Deliberately not checked / assumptions.** I did not run any P56 code (P56's
`coverings.py` is a committed *plan*, not yet implemented in `src/` — the stated
prerequisite; the plan flags "STOP and escalate if P56 has not merged", and the P56 addendum
commits P56 to retain + expose `minimal_relation_counts`). I did not byte-extract the BJP
2019 book chapter or the compressed arXiv PDFs (paywalled / FlateDecode); the full
hypercritical list is `# PIN`ned — the plan encodes what is transcription-checkable (the
minimal wild trees, defects verified in read-only Python) plus whatever the implementer
confirms, records `HYPERCRITICAL_COVERAGE`, and returns honest `None` outside coverage
(exactly as R19 instructs, "# PIN what you cannot transcribe now"). I did NOT retain any
guessed entry box for weak nonnegativity — `_WNN_BOX` as a *decision* mechanism is deleted;
the bounded sweep survives only as a sound witness-FINDER. I did not fix the `_star`-tree
orientations in code (irrelevant to the symmetric form; the implementer rechecks against
`1/p+1/q+1/r`). I recommend but did not force the Bongartz "rep-infinite-but-tame/wild-
withheld" middle state to `None` + reason (decide-and-pin at implementation).

**Why this is correct.** The two verdict equivalences are quoted theorems with verified
citations; their hypotheses (strongly/simply simply connected, algebraically closed) are
enforced as loud gates consuming P56's certificate. Weak positivity rests on Ovsienko's
cited box-6 (no home-grown proof after the review withdrew the flawed one); weak
nonnegativity rests on the classified hypercritical list with explicit coverage, emitting a
certified `tame` only inside coverage and honest `None` outside — the `T_{2,3,7}` refutation
proves no shortcut box is sound. Every "no" (wild) is a checkable exact witness. And the
independent oracles (Ext² vs the combinatorial count; the AR knit vs weak positivity; the
classified lists vs the sweep) cross-validate each layer from a genuinely different engine.

---

## Change log

- **2026-08-07 adversarial review: 3 blocking + 3 majors applied** — weak-nonnegativity
  decided by the classified hypercritical list (primary), no guessed boxes, tame only under
  recorded list-completeness, honest `None` otherwise; `T_{2,3,7}` 10-variable refutation
  documented (the `≤ 9`-variable transcription was false); `tits_form_combinatorial` rename
  (W1, P38's `tits_form` untouched); `is_acyclic` unit-form guard (M2, the `G_ii == 2`
  heuristic wrongly passes `k[x]/(x³)`); branch-and-bound + real budget arithmetic (M1,
  `E8` `6⁸ = 1.68M`, `Ẽ8` `6⁹ ≈ 10.1M`); the `form_type` (definiteness) vs `rep_type`
  (certified) GUI/report disambiguation (W2); Ovsienko box-6 rested on the citation and the
  flawed inline "∎" proof withdrawn; P56 addendum committing `minimal_relation_counts` (W3).
