# Plan 31 — TrivialExtension double-quiver presentation (Tier 1a)

- **Date:** 2026-07-26
- **Branch:** `plan-31-trivial-extension`
- **Backlog item:** Tier 1a "TrivialExtension double-quiver presentation" (found by
  Plan 29's `is_symmetric` fix).
- **Status:** in flight — executed multi-agent (3 research agents, then parallel
  implementation/test/docs agents; this doc is the contract they build against).

## Goal

`TrivialExtension(A)` of a quiver-presented `A` returns a **genuine
`kQ_T/I_T`-presented Algebra** (built via `Quiver.algebra`), so the path-basis
invariants (`is_symmetric`, `is_frobenius`, `is_weakly_symmetric`,
`is_selfinjective`, `cartan_matrix`, modules, `engine="cs"`, …) serve `T(A)`.
The Plan-29 trace-form certifier then returns `True` on every `T(A)` and the
`xfail(strict=False)` fences in `tests/invariants/test_symmetric_regression.py`
flip to real asserts.

## Research findings (2026-07-26, three parallel agents; all claims machine-verified)

1. **The construction is classical and was verified end-to-end.** Vertices of
   `Q_T` = `Q_0`. New arrows = duals of a corner-homogeneous basis of the
   **bimodule socle** `soc_{A^e}A = {a : (rad A)·a = 0 = a·(rad A)}`, with
   **reversed direction**: `w ∈ e_i A e_j` (a path-combination `i → j`, left-to-right
   convention) yields `β_w : j → i`. Proven consistent with the ⋉ coordinate
   formulas in the current `trivial_extension.py` (lines 25–34:
   `(a·f·b)(c) = f(b·c·a)`).
2. **QPA ships a native construction** — `TrivialExtensionOfQuiverAlgebra(A)`
   (QPA 1.37, `pathalgtensor.gd:67`), probed live: T(kA₂) → dim 6, quiver
   `{a: 1→2, te: 2→1}`, relations `{a·te·a, te·a·te}`, `IsSymmetricAlgebra = true`.
   This is a **construction oracle** (dim + arrow count + three self-injectivity
   predicates), stronger than a symmetry-verdict-only crosscheck. It retires the
   "T(A) is out of QPA scope" note.
3. **Worked anchors (frozen test pins, hand-derived AND machine-confirmed):**
   - `T(kA_n) ≅ kZ_n/J^{n+1}` (symmetric cyclic Nakayama / Brauer star,
     `n | (L−1)` with `L = n+1`); verified n = 2, 3, 4 (dims 6/12/20, Cartan,
     Loewy `n+1`, all four symmetry booleans).
   - `T(k[x]/(x^a)) = k⟨x,β⟩/(x^a, β², xβ−βx)` — the relation is the **plain
     commutator in every characteristic** (D(A) is an honest bimodule, no Koszul
     sign; "char 2" is moot). For a = 2 this is `k[x,y]/(x²,y²)`, reproducing the
     existing `HH_• = [4,4,5,6]` oracle on the verification page.
   - 2-Kronecker: 2 new arrows (bimodule socle `⟨a,b⟩`), dim 8, six quadratic
     relations (`a·q, b·p, p·b, q·a, a·p − b·q, p·a − q·b`).
   - Commutative square (`ac − bd`): one new arrow `4 → 1`, dim 18.
   - Cartan identity `C_T = C_A + C_Aᵀ` (in the repo convention
     `C[i][j] = dim e_i A e_j`); always symmetric.
4. **The decisive back-compat finding:** annotating the ⋉ structure constants
   with a quiver is NOT enough — `invariants/cartan.py`, `modules/builders.py`,
   `modules/module.py`, `modules/opposite.py`, `modules/tor.py` parse
   `basis_labels` as literal composable paths and **crash** on dual labels.
   The returned object must be a genuine `Quiver.algebra` output (real path
   labels). All consumed numbers are iso-invariant (confirmed: identical HH on
   both builds).
5. **`is_symmetric` needs no change:** a hand-built presented T(kA₂) already
   certifies `True` over QQ, CC, GF(2), GF(32003) through the existing
   Skowroński–Yamagata trace-form certifier. GF(2) is safe here — the positive
   branch is a field-agnostic full-rank Gram witness; the "loud on small fields"
   caveat only bites negatives.

## Design decisions

- **D1 — Algorithmic relations, certified per instance (not a transcribed
  theorem).** Compute `I_T` as a generating set of `ker(π: kQ_T → T_⋉)` by a
  length-lex mini-Gröbner kernel extraction (process words by increasing length;
  keep a reduced-echelon basis of images of normal words; extend only normal
  words; a word whose image reduces against prior normal forms emits a relation).
  Filtration-aware, so non-homogeneous base relations are handled. The
  Fernández–Platzeck closed form is **not** transcribed (primary source not in
  hand; repo rule: nothing unverified becomes an implementation or strict pin) —
  the classical special cases above serve as oracles instead.
- **D2 — The certificate is the dimension check.** π is surjective by
  construction and every emitted relation is in `ker π`, so
  `dim(kQ_T/I_T) ≥ 2·dim A` always, with equality iff iso. Build, then require
  `B.dim == 2*A.dim`; raise `QuiverlabError` loudly otherwise. `< 2·dim A` is
  impossible; over-generation of relations is harmless. Admissibility holds by
  construction (`I_T ⊆ J²` since arrows lift a basis of `rad T/rad²T`;
  `J^{L(T)+1} ⊆ I_T`). Kernel enumeration runs to length `L(A) + 2`.
- **D3 — Scope gate + fallback.** Presented route iff (a) `A` has a path-type
  basis (`invariants/pathbasis.py::path_type_basis` succeeds — probe, don't
  parse) and (b) the domain's relation coefficients are string-representable for
  `Quiver.algebra` (QQ → `Fraction` strings, GF(p) → ints; covers the whole
  battery). Otherwise return the **unchanged** ⋉ structure-constants build (kept
  verbatim as private `_trivial_extension_structure_constants`) — honest
  refusals downstream are preserved, and the old build doubles as the
  iso-invariance oracle. GF(pⁿ)/CC-algebraic with a non-monomial socle dual
  basis falls back (documented; the certificate never lets a wrong algebra
  through).
- **D4 — Dual-arrow naming.** Arrow names must match `^[A-Za-z_][A-Za-z0-9_]*$`
  (the current `x*` labels are illegal). Use `te0, te1, …` (corner-homogeneous
  socle-basis order), with a deterministic disjointness guard against existing
  arrow names (`Quiver`'s dict silently overwrites duplicates): while any `te{s}`
  collides, prefix underscores (`_te0, …`) until disjoint.
- **D5 — Citations.** No Fernández–Platzeck key (metadata not verifiable to
  BibTeX precision — the construction is per-instance-certified + QPA-oracled,
  so no theorem citation is load-bearing). Add `happel_trivial_extension`
  (Happel 1988, LMS Lecture Note Series 119, CUP — metadata already recorded in
  `docs/plans/2026-07-25-literature-oracles-deep-research.md`; verify there
  before adding). `_family_citations` of the presented build:
  `("assem_book", "skowronski_yamagata", "happel_trivial_extension",
  "cmrs_split")`; fallback keeps `("assem_book",)`.
- **D6 — QPA crosscheck, both variants.** (A) Feed our presented `T(A)` through
  the existing `crosscheck_symmetric` route (works today, zero new plumbing).
  (B) New `trivial_extension` crosscheck: QPA builds
  `TrivialExtensionOfQuiverAlgebra(A)` natively; compare dim, arrow count,
  `IsSymmetricAlgebra`, `IsWeaklySymmetricAlgebra`, `IsSelfinjectiveAlgebra`
  (`qpa/scripts.py::trivial_extension_script` +
  `qpa/crosscheck.py::crosscheck_trivial_extension`, registered in the
  dispatcher). QPA arrow labels differ (`te_a1_i_j`) — compare counts, not names.

## Implementation surface (file ownership per agent)

**Agent C1 — src/:**
- `src/quiverlab/families/trivial_extension.py`: the presented route (socle via
  `fields.linalg.nullspace`, dual covectors via `solve` — never a Gram inverse;
  kernel extraction; relation-string formatting with coefficients before arrows,
  signs folded into term separators; the D2 certificate), `_trivial_extension_structure_constants`
  fallback, D3 gate, D4 naming, D5 citations.
- `src/quiverlab/families/discover.py:39-41`: route `"structure-constant"` →
  `"general"`; citations tuple extended (must stay registry-valid —
  `test_catalog_citations_are_all_registered`).
- `src/quiverlab/citations/references.bib` + `registry.py`: the
  `happel_trivial_extension` key (D5).
- `src/quiverlab/qpa/scripts.py` + `crosscheck.py`: D6 variant B + dispatcher;
  fix the now-false "TrivialExtension cannot be fed to QPA" docstring clause.

**Agent C2 — tests/:**
- NEW `tests/families/test_trivial_extension_presented.py` (deep): structural
  (dim = 2·dim A on all six bases incl. `_validate()`; new-arrow count =
  bimodule-socle dim; `T(kA_n)` ≡ `NakayamaAlgebra(n, n+1, cyclic=True)`
  invariants n = 2,3,4; `T(k[x]/(x^a))` ≡ two-loop presentation a = 2,3 incl.
  the `[4,4,5,6]` cross-link; Cartan identity incl. zoo `line_abc_cde`);
  symmetry over GF(32003)/QQ/GF(2) (GF(2) only on the Nakayama/dual-number
  images — cut-list rank 1); iso-invariance (presented vs `_trivial_extension_structure_constants`
  bar-HH degreewise: dualnum `[4,4,5,6]` top 3, kA₂ `[3,1,1,1]` top 3,
  Kronecker `[3,4,6]` top 2; HHⁿ = HH_n); cross-engine CS ≡ bar on presented
  T(kA₂) (`[3,1,1]` top 2 — CS was refused pre-plan); negatives
  (presentation-less base → fallback, `quiver is None`, invariants raise
  `FieldError`, dim/bar-HH still serve; infinite-dim base refused upstream);
  Loewy/selfinjective (`loewy_length`: T(kA_n) = n+1, T(k[x]/x^a) = a+1,
  T(Kronecker) = 3); center = HH⁰ (n+1 for kA_n).
- `tests/invariants/test_symmetric_regression.py`: the fence flip — remove the
  `@pytest.mark.xfail`, rewrite the root-cause comment, strengthen the body
  (`quiver is not None`, + `is_weakly_symmetric`, + `is_selfinjective`).
- `tests/families/test_trivial_extension_hh1.py`: line 48 `quiver is None`
  inverts; the CS-refusal test (53–61) is rewritten as the CS-now-computes
  strengthening (or folded into the new file); stale docstrings refreshed.
- `tests/families/test_trivial_extension.py`: stale "no quiver" comments
  refreshed; optional `is_symmetric` tightening.
- NEW `tests/qpa/test_trivial_extension_qpa.py` (qpa): D6 variants A + B with
  the live-verified pins (kA₂ 6/2 arrows, kA₃ 12/3, Kronecker 8/4, dualnum 4/2,
  comm-square 18/5; all predicates true).
- `tests/qpa/test_symmetric_qpa.py`: retire the obsolete SCOPE NOTE (lines 12–16).
- Check `tests/families/test_freshness_gate.py` (route/interface snapshot may
  need re-pinning).

**Agent C3 — docs/ (+ bookkeeping):**
- `docs/verification.md`: retire the honest-scope refusal entry (~lines
  400–404) with the delivered-presentation text; add the Plan-31 oracle bullet
  (Class-1 literature list); drop the QPA-scope caveat and note the native
  `TrivialExtensionOfQuiverAlgebra` oracle; subsystem-row descriptions. Counts
  left to the final recount.
- `docs/internals/11-families-citations.md`: TrivialExtension leaves the
  structure-constant route (TensorProduct stays).
- `docs/plans/DEEPER-ENGINES-BACKLOG.md`: tick the Tier-1a item (Plan 31, date,
  branch).
- `CLAUDE.md`: status prose (line ~325 refusal note superseded; Plan-31 entry).
  Counts left to the final recount.

**Final (orchestrator):** run fast + families/invariants deep selections +
`-m qpa`; recount collect-only; sync counts across `README.md` (badge + prose),
`docs/verification.md`, `CLAUDE.md` (a release test pins badge ==
verification.md); commit.

## Acceptance

1. The four fence tests pass as **real asserts** (kA₂, kA₃, 2-Kronecker,
   commutative square over QQ): `is_frobenius`, `is_symmetric`,
   `is_weakly_symmetric`, `is_selfinjective` all `True`.
2. Every presented build certifies `dim == 2·dim A`; a failed certificate
   raises loudly; presentation-less bases keep the exact old behavior.
3. The battery above green (fast + deep + qpa), zero regressions elsewhere
   (`test_trivial_extension.py` dim/duality values, HH¹ ≠ 0 pins, discover
   catalog, freshness gate).
4. `docs/verification.md` carries the Plan-31 oracles + citations per the
   standing rule; counts synced badge == page == CLAUDE.md.
5. No floats in `src/`; engine internals untouched; merge/push only when Marco
   asks.
