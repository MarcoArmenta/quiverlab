# Plan 60: the tilted-algebra recognizer (P60 / R17) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A no-code **tilted-algebra recognizer**. Given a connected, basic, quiver-presented
`A = kQ/I`, decide whether `A` is **tilted** (i.e. `A ≅ End_H(T)` for a hereditary algebra `H`
and a tilting `H`-module `T`) by the **Liu–Skowroński criterion**: search the Auslander–Reiten
quiver `Γ_A` for a **faithful section `Σ`** (a slice) with `Hom_A(X, τY) = 0` for all
`X, Y ∈ Σ`. On **tilted**, return the **slice** (the section's modules), the **slice/tilting
module** `S = ⊕Σ`, the reconstructed **hereditary algebra** `H = End_A(S)` (presented as
`kQ/I`), and the **hereditary type** (the Dynkin diagram underlying `Σ`). On **not tilted**, a
genuine refutation within the representation-finite complete-knit + budget scope. Every verdict
is either a **theorem gate** (hereditary ⇒ tilted; non-semisimple self-injective ⇒ not tilted;
`gl.dim A > 2` ⇒ not tilted) or an **exhaustive faithful-section search** on the knitted
`Γ_A`, and it **self-certifies** (Ringel's slice theorem: `S` is a tilting `A`-module **and**
`End_A(S)` is hereditary) or **refuses loudly**. Exact only — no floats. This plan closes the
`# PIN` P55 left open: it flips P55's skipped `test_support_components_are_tilted_PIN` into a
real assert (each connected component of the left/right support algebras `A_λ`/`A_ρ` is
certified tilted by **this** recognizer).

**Architecture:** ONE new module, a thin exact layer over primitives that already exist
(P41 AR knitting, P44 tilting + Gabriel-quiver recovery, P37 hom/End, P55 the left/right atlas):

- **`src/quiverlab/modules/tilted.py`** (new) — the whole surface. Like `modules/left_right.py`
  (P55) and `modules/ar.py` (P41), this verb ENUMERATES and COMPARES whole iso-classes of
  modules on the knitted indecomposable universe, so it is a `modules/` citizen and its battery
  lives in `tests/modules/` (the **deep** bucket, `tests/conftest.py`). Public surface:
  - `tilted_check(A, *, budget_modules=256, budget_sections=4096) -> TiltedReport` — the primary
    compute: theorem gates first, then (rep-finite, non-self-injective, `gl.dim ≤ 2`) the
    faithful-section search on `Γ_A`, then the Ringel slice-theorem certificate. Two budgets: the
    **knit cap** (`budget_modules`, indecomposable count) and the **transversal cap**
    (`budget_sections`, a product of orbit sizes — a different, larger quantity; see W3).
  - `@dataclass TiltedReport` — mirrors P41 `ARQuiver` / P55 `LeftRightAtlas` honest
    `is_complete`/`status` contract; carries the verdict, the slice, the slice module, the
    reconstructed `H`, the hereditary type, and the certificate trail; JSON-ready summaries plus
    the actual `Module`s / `Algebra` for in-process consumers.
  - `tilted_check_block(A, *, budget_modules=256, budget_sections=4096) -> dict` — the JSON block
    the two runners share.
  - Internal seams the oracle battery calls directly: `_faithful_section_search(ar, A,
    budget_sections)` (returns the section indices or `None`), `_is_faithful(A, M)`,
    `_hom_tau_zero(A, section)`, `_certify_slice(A, S)` (Ringel 1.9(2): tilting + hereditary
    presented End), `_section_quiver(ar, indices)` + `_type_str` (over the shipped
    `invariants/dynkin_type.py::dynkin_type`).
  - Thin `Algebra` delegates in `core/algebra.py` beside `is_tilting_module` (P44,
    `core/algebra.py:469`): `Algebra.tilted_check(budget_modules=256, budget_sections=4096)` and
    the convenience `Algebra.is_tilted(budget_modules=256, budget_sections=4096) -> bool`
    (**exactly the method name P55's skipped test calls** — see Task C).

No new math engines; no new AR machinery. Everything is a bounded search over the finite
indecomposable universe the P41 knit already produces, plus the shipped exact `hom_dim`,
`is_tilting_module`, `end_algebra`, `presented_form`, `is_hereditary`, `global_dimension`,
`is_selfinjective`, and `Module.tau`.

**Tech Stack:** exact `Domain` linear algebra throughout — `modules/hom.py::hom_dim` (exact in
**every** characteristic — the `Hom(X, τY)=0` gate and faithfulness carry no char caveat),
`modules/duality.py::tau` (`Module.tau()`; `= 0` on projectives), `modules/tilting.py::
is_tilting_module` (P44 — the Bongartz count criterion; **char 0 / char > dim** via `decompose`),
`modules/endomorphism.py::end_algebra` (structure-constant `End_A(S)`),
`core/basic.py::presented_form` (P44 — Gabriel-quiver recovery of `End_A(S)` as `kQ/I`; **char 0
/ char > dim**), `invariants/recognizers.py::is_hereditary`, `modules/ext.py::global_dimension`,
`modules/ext.py::is_selfinjective`, `Algebra.ar_quiver` (P41), `left_right.py::left_right_parts`
(P55 — the shared knit + section seed). No floats in `src/` (AST-gated by
`tests/test_no_floats.py`); all dim-vectors are `dict[vertex,int]`, all dims `int`.

---

## Record (R17, verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R17 — Tilted-algebra recognizer (Liu–Skowroński faithful section).** [C-scout P5;
> keep-with-corrections] Object: find/refute a faithful section Σ with Hom(X, τY) = 0 (|Σ₀| = n
> is IMPLIED, state as consequence); returns slice, tilting module, hereditary type. Rep-finite
> via exhaustive section enumeration; the local criterion (arXiv:1409.2054) is the rep-infinite
> extension path. Refs: Liu, Arch. Math. 61 (1993) 12–19 (venue verified); Skowroński;
> arXiv:1409.2054. Oracles: hereditary ⇒ tilted; kZ₃/J² stable tube ⇒ not; cluster-tilted A₃ ⇒
> not. Size M.

Tier β, C-cluster recognizer ladder; dependency chain (research doc §9): **R15 → R16 → {R17,
R18, R21, R37}**. P60 (R17) `needs P55` (metaplan §4 Wave 2). **P55 (the left/right substrate) is
NOT yet merged** — only its plan is committed (`c08fdf7`); `src/quiverlab/modules/left_right.py`
does not exist, and the metaplan §6 ledger shows all of Wave 2 unimplemented. **P55 must be
implemented and merged FIRST**; `plan-60-tilted-recognizer` branches off `dev` **after** P55
lands (the binding **Task 0** STOP gate below enforces this).

---

## Reference re-verification (mandatory standing rule; done at spec time)

Re-verified via web + **the arXiv PDF of 1409.2054 read end to end (pages 1–9)**. **Findings:**

1. **Liu, S. — "Tilted algebras and generalized standard Auslander–Reiten components",
   Archiv der Mathematik (Basel) 61 (1993), no. 1, 12–19**, DOI `10.1007/BF01258050`. Venue,
   volume, year, and page range verified (Springer; Internet Archive TOC of Arch. Math. 61
   (1993)). This is the origin of the "faithful generalized-standard AR component / section"
   characterization the record names. `# PIN`: the exact issue number (no. 1) confirmed at
   citation-add time against the Springer landing page.

2. **Liu, Shiping — "Another characterization of tilted algebras", arXiv:1409.2054** (submitted
   6 Sep 2014; MSC 16G70, 16G20, 16E10). Title, sole author, and the mathematical content are
   verified **verbatim from the PDF** — the definitions and theorems below are transcribed from
   it. This is the record's "local criterion — the rep-infinite extension path": the paper's
   whole point is that a **cut** is a *finite/local* object (weakly convex rather than convex),
   so tiltedness is checkable **without** knowing the entire (possibly infinite) AR component.
   `# PIN`: the paper was later published (the arXiv v1 gives no journal ref); the worker either
   verifies and upgrades the entry to `@article` at merge, or ships the verified `@misc` arXiv
   entry — house rule, no fabricated bibliographic fields.

3. **Ringel's slice theorem** — Liu's Theorem 1.9 is quoted from his reference **[18]** (Ringel,
   *Tame algebras and integral quadratic forms*, LNM 1099, Springer 1984, (4.2)); the tilted
   class itself originates with **Happel–Ringel, "Tilted algebras", Trans. Amer. Math. Soc. 274
   (1982), no. 2, 399–443** (Liu's [7]). Both are the standard, BibTeX-verifiable references;
   the section/slice/tilted-algebra material is also in the shipped **ASS2006** (Assem–Simson–
   Skowroński, *Elements of the Representation Theory of Associative Algebras*, Vol. 1),
   Chapters VI (tilting) and VIII (tilted). `# PIN`: Happel–Ringel volume/pages and the Ringel
   LNM number confirmed at citation-add time.

4. **Skowroński** (record's third attribution) — the "Liu–Skowroński criterion" pairs Liu 1993
   with **Skowroński, "Generalized standard Auslander–Reiten components", J. Math. Soc. Japan 46
   (1994), 517–543** (found the characterization independently). `# PIN`: venue/pages confirmed
   at citation-add time; if not BibTeX-verifiable it is folded into the Liu1993 annotation
   rather than shipped as a fabricated entry.

### Definitions and theorems extracted verbatim from arXiv:1409.2054

Throughout, `A` is a (basic, connected) artin algebra, `mod A` its f.g. modules, `ind A` the
indecomposables, `Γ_A` the AR quiver, `τ = DTr`, `τ⁻ = TrD` (with the convention `τX = 0` if
`X` is projective, `τ⁻X = 0` if `X` is injective).

- **Tilting module (§1).** "A module `T` in mod `A` is called *tilting* if `pdim(T) ≤ 1`,
  `Ext¹_A(T, T) = 0`, and the number of non-isomorphic indecomposable direct summands of `T` is
  equal to the number of non-isomorphic simple `A`-modules." (⟵ this is P44's
  `is_tilting_module`, `n=1`, Bongartz count criterion.)
- **Tilted (§1.8 preamble).** "Recall that `A` is *tilted* if `A = End_H(T)`, where `H` is a
  hereditary artin algebra and `T` is a tilting module in mod `H`."
- **Faithful / sincere (§1).** "The *annihilator* of `Σ`, written `ann(Σ)`, is the intersection
  of all annihilators `ann(M)` with `M ∈ Σ`. One says that `Σ` is **faithful** if `ann(Σ) = 0`
  and **sincere** if every simple `A`-module is a composition factor of some module in `Σ`."
- **Section (§1).** "`Σ` is a *section* in a connected component `Γ` of `Γ_A` if `Σ` is a
  connected subquiver of `Γ`, which contains no oriented cycle, meets each `τ`-orbit in `Γ`
  exactly once, and is convex in `Γ`, that is, every path in `Γ` with end-points belonging to
  `Σ` lies entirely in `Σ`." (Liu's [13,(2.1)].)
- **Slice — Def 1.8 [18].** "A full subquiver `Δ` of `Γ_A` is called a *slice* if: (1) `Δ` is
  sincere and convex in `ind A`; (2) if `X ∈ Δ`, then `τX ∉ Δ`; (3) if `X → Y` is an arrow in
  `Γ_A` with `Y ∈ Δ`, then either `X` or `τ⁻X` belongs to `Δ`."
- **Theorem 1.9 [18] (Ringel — THE RECONSTRUCTION, transcribed verbatim).** Let `A` be an artin
  algebra, and let `Δ` be a full subquiver of `Γ_A`.
  1. If `A = End_H(T)` with `H` hereditary and `T` a tilting `H`-module, then `T` determines a
     slice in `Γ_A` generated by the direct summands of `Hom_A(T, D(H))`.
  2. **The subquiver `Δ` is a slice if and only if `S = ⊕_{X∈Δ} X` is a tilting module in
     `mod A` such that `H = End_A(S)` is hereditary. In this case, `D(S_H)` is a tilting
     `H`-module such that `A = End_H(D(S))` and `Δ` is the slice determined by `D(S)`.**
- **Cut — Def 2.1.** "A full subquiver `Δ` of `Γ_A` is called a *cut* if, for each arrow
  `X → Y`: (1) if `X ∈ Δ`, then either `Y` or `τY`, **but not both**, belongs to `Δ`; (2) if
  `Y ∈ Δ`, then either `X` or `τ⁻X`, **but not both**, belongs to `Δ`." (Remark: a section in a
  connected component is a cut — Liu's [15,(2.2)]; a cut is a *presection*.)
- **Proposition 2.5.** Let `Δ` be a cut of `Γ_A`. TFAE: (1) `Δ` is finite and weakly convex in
  `ind A`; (2) `Hom_A(X, τY) = 0` for all `X, Y ∈ Δ`; (3) `Hom_A(τ⁻X, Y) = 0` for all
  `X, Y ∈ Δ`. In this case `Δ` contains no oriented cycle.
- **Theorem 2.6 (THE MAIN CRITERION, transcribed verbatim).** "Let `A` be an artin algebra.
  Then `A` is tilted **if and only if** `Γ_A` contains a faithful cut `Δ` such that
  `Hom_A(X, τY) = 0` for all `X, Y ∈ Δ`; and in this case, `Δ` is a slice in `Γ_A`."
- **Theorem 2.7 (tilted quotients — noted, out of v1 scope).** If `Δ` is a cut with
  `Hom_A(X, τY) = 0` for all `X, Y ∈ Δ`, then `B = A/ann(Δ)` is a tilted algebra with `Δ` a
  slice of `Γ_B`. (Since our accepted `Δ` is **faithful**, `ann(Δ) = 0` and `B = A`; Theorem
  2.7's genuine extension — tilted *quotient* algebras of a non-tilted `A` — is deferred.)
- **The faithfulness-is-essential counterexample (Liu's Example after Thm 2.6, verbatim).** "Let
  `A` be an algebra with radical squared zero given by the quiver [3-vertex quiver `a, b, c`,
  `a` a source, `b—c` an edge]. … `Δ: P_b — S_b — P_a` is a **sincere** cut in `Γ_A` such that
  `Hom_A(X, τY) = 0` for all `X, Y ∈ Δ`. However, `A` is **not tilted**." ⟵ **the faithfulness
  of the cut cannot be weakened to sincereness**; this is the honest-scope oracle for `_is_
  faithful` (Task D). `# PIN`: the worker transcribes the exact arrow orientation of the 3-vertex
  `rad²=0` quiver from the PDF figure (p. 9) and verifies live that it is rep-finite,
  non-self-injective, carries the sincere-but-non-faithful cut, and `tilted_check` returns
  `not_tilted`.

**Two consequences the plan leans on (stated, with one-line proofs):**

- **`|Σ₀| = n` and the refutation-completeness chain (the record's "IMPLIED, state as
  consequence").** The chain is: **(i)** the Theorem-2.6 witness `Δ` is a **slice** (Thm 2.6:
  "*in this case, `Δ` is a slice in `Γ_A`*"). **(ii)** A **slice is a section** — Liu's Remark
  following Theorem 1.9 (transcribed from the PDF, p. 4): "*a slice `Δ` of `Γ_A` is necessarily
  connected … `Δ` is contained in a connected component `C` of `Γ_A`, which is actually a section
  in `C`; see [7](7.1) and [18](4.2)*" (Happel–Ringel, Trans. AMS (7.1); Ringel LNM 1099 (4.2)).
  `# PIN`: the parallel ASS VIII statement number ("every complete slice is a section", ASS
  Vol. 1, Ch. VIII) is not web-confirmed here — the worker cites the confirmable Liu Remark and
  optionally the ASS number at merge. **(iii)** A **section** "*meets each `τ`-orbit in `Γ`
  exactly once*" (Liu's §1 definition, quoted above). Chaining (i)–(iii): the witness is a
  **one-module-per-`τ`-orbit transversal** of a single connected component — which is exactly why
  **enumerating transversals is complete** (the witness, if it exists, is among them, so
  exhausting the transversals within budget is a genuine *refutation*). And since the slice module
  `S = ⊕Δ` is a tilting `A`-module with `n = |Q₀| = rk K₀(A)` non-isomorphic summands
  (tilting-module definition, §1: `#summands = #simples`), the hosting component has **exactly `n`
  `τ`-orbits**, so `|Δ₀| = n` and the **`#orbits == n` component prune is sound**. **A cut in the
  sense of Def 2.1 need NOT be a transversal** — the transversal property is a *slice/section*
  property, invoked only for the Theorem-2.6 witness; nothing in the search relies on cuts being
  transversals. This is a **search pruning + completeness** argument, **never** the definition.
- **The theorem gates (necessary conditions that make three oracles instant).**
  - **`tilted ⇒ gl.dim ≤ 2`.** `A = End_H(T)`, `H` hereditary (`gl.dim H ≤ 1`), `T` tilting, so
    `gl.dim A ≤ gl.dim H + 1 = 2` (the tilting theorem, Happel–Ringel / ASS VI.4). Contrapositive:
    **`gl.dim A ≥ 3` (finite or a certified lower bound) ⇒ not tilted.**
  - **Non-semisimple self-injective ⇒ not tilted.** A self-injective algebra of finite global
    dimension is semisimple (`gl.dim 0`): over a self-injective algebra a module has finite pd
    iff it is projective, so `gl.dim < ∞ ⇒` every module projective `⇒` semisimple. Combined with
    `tilted ⇒ gl.dim ≤ 2 < ∞`: **`A` self-injective and non-semisimple ⇒ `gl.dim A = ∞` ⇒ not
    tilted.** This is the clean gate for the `kZ₃/J²` oracle, which the P41 knit **refuses**
    (self-injective, `status="unsupported"`) — so the verdict comes from the gate, a **real
    `not_tilted`**, not a refusal.
  - **Hereditary ⇒ tilted.** `A = End_A(A)` with `A` a tilting module over itself and `H = A`
    hereditary; the slice is the postprojective section of the indecomposable projectives. This
    gate gives the correct `tilted` verdict even for **representation-infinite** hereditary
    algebras (e.g. the Kronecker quiver) that the knit would refuse.

---

## The P55 connection (the load-bearing dependency)

P55 (`src/quiverlab/modules/left_right.py`) — **which must be implemented and merged before P60
begins** (Wave-2 order; see Task 0) — will ship, per its committed plan (`c08fdf7`), an
auto-flipping fence: a **skipped** test in `tests/modules/test_left_right_support.py`:

```python
@pytest.mark.skip(reason="PIN: each support component is tilted -- VERIFIED in P60 "
                         "(tilted recognizer). Auto-flips to a real assert when P60 lands.")
def test_support_components_are_tilted_PIN():
    atlas = left_right_parts(linear_path_algebra(3, field=QQ))
    for comp in atlas.left_support.components:
        assert comp["algebra"].is_tilted()             # P60 API
```

This plan (Task C) makes `Algebra.is_tilted()` real and **flips this fence into a real assert**
— and, crucially, extends it to the **ACLV Example 2.2(b) `rad²=0 A₅`** shared fixture, where
P55 established `A_λ` on vertices `{1,2,3}` and `A_ρ` on `{3,4,5}` (each `= kA₃/rad²`). This
plan certifies **each connected component of `A_λ` and `A_ρ` is tilted** — the concrete content
of ACLV's "`A_λ` is a product of tilted algebras" theorem, which P55 `# PIN`'d for P60. P60 MAY
also consume `atlas.ext_injectives_left` / `atlas.ext_projectives_right` (the left/right
sections, ACLV Thm 3.1) as a **fast-path section seed** before falling back to the exhaustive
search (an optimization, guarded and documented — the exhaustive knit search stays authoritative).

**Live-verified facts (all recomputed in the venv at spec time, QQ):**

| algebra | `gl.dim` | self-inj | knit | verdict | via | slice / type |
|---|---|---|---|---|---|---|
| `kA₃` (`1→2→3`) hereditary | `1` (exact) | no | complete, 6 indec, 3 orbits | **tilted** | Gate H | projectives; type `A₃` |
| `kA₃/rad²` (`2→1, 3→2`) = P55 `A_λ` | `2` (exact) | no | complete, 5 indec, 3 orbits | **tilted** | search | slice `{S₂, P₂, P₃}` (shape `P₂→S₂→P₃`); `End = kA₃` hereditary; type `A₃` |
| `kZ₃/J²` = cluster-tilted `A₃` = `Jac(3-cycle, αβγ)` | `∞` (lower bd 32) | **yes** | refused (`unsupported`) | **not tilted** | Gate S | — |
| `rad²=0 A₅` (P55 fixture) | `4` (exact) | no | complete, 9 indec, 5 orbits | **not tilted** | Gate G / search | — |

**Cluster-tilted `A₃` = `kZ₃/J²` (a settled correction to the record).** The record lists
"`kZ₃/J²` ⇒ not" and "cluster-tilted `A₃` ⇒ not" as two oracles. Live computation shows the
type-`A₃` cluster-tilted algebra (the 3-cycle `1→2→3→1` with all length-2 paths zero — the
Jacobian algebra of the 3-cycle with potential `αβγ`, `∂W = {βγ, γα, αβ}`) **IS** `kZ₃/J²`, the
self-injective cyclic Nakayama algebra of Loewy length 2 (cluster-tilted algebras are Gorenstein
of dimension `≤ 1`, so being Gorenstein-dimension-0 self-injective is permitted). Both oracles
are therefore the **same** self-injective algebra, refuted by the **same** Gate S. The plan
records this honestly and adds a **genuinely rep-finite, non-self-injective search-refutation
oracle** so the search's `not_tilted` path is exercised by a real example (the `rad²=0 A₅` via
the direct `_faithful_section_search` seam, plus Liu's `# PIN` sincere-not-faithful example).

**Why the record's "cluster-tilted `A₃` ⇒ not tilted" must mean the non-hereditary
representative.** Type `A₃` has *two* cluster-tilted algebras up to isomorphism: the **hereditary
`kA₃`** itself (from the "line" cluster-tilting objects — a cluster-tilted algebra *can* be
hereditary), and the **non-hereditary 3-cycle** `Jac(3-cycle, αβγ) = kZ₃/J²`. The hereditary
`kA₃` **is tilted** (Gate H) — so if the record's oracle referred to *it*, the oracle would read
"tilted", contradicting the record. The oracle can therefore only be the **non-hereditary
representative**, `kZ₃/J²` — which is exactly why the record's two "not tilted" oracles
(`kZ₃/J²` and cluster-tilted `A₃`) collapse to one self-injective algebra under Gate S. The plan
states this so the merge is a *reasoned identification*, not an accident of a builder choice.

---

## Mathematical foundation (the verdict pipeline — the plan's ground truth)

`tilted_check(A, *, budget_modules=256, budget_sections=4096)` requires a **quiver presentation**
(`A.quiver is not None`; the GUI/webapp always provide one — a structure-constant-only `A`
refuses **loudly** with `QuiverlabError`, like P55's support build). It runs, in order:

1. **Gate H — hereditary ⇒ tilted (with a certified payload where computable).** If
   `is_hereditary(A)` (no relations and `Q` acyclic): verdict `tilted`, `reason="hereditary"`,
   `H = A`, hereditary type = `_type_str(dynkin_type(A.quiver))`. **For rep-finite hereditary
   `A`** (the knit completes within `budget_modules`): produce the projective section `Σ = {P_v}`
   AND run `_certify_slice(A, ⊕Σ)` on it — the slice module is `A_A`, `End_A(A_A) ≅ A^op` is
   hereditary, so the certificate holds; ship the certified slice + `hereditary_algebra`. **For
   rep-infinite hereditary `A`** (the knit refuses — e.g. the Kronecker quiver): the verdict
   `tilted` + type stand **without** the module-level slice; the slice is **omitted** with the
   honest note "slice omitted: representation-infinite hereditary; the postprojective section
   exists but is not knit-enumerable" (a record-level omission, stated in the honest-scope
   section). No knit is needed for the *verdict* — this is why rep-infinite hereditary algebras
   get the correct answer.
2. **Gate S — non-semisimple self-injective ⇒ not tilted.** Else if `A.is_selfinjective()`
   (exact, any field; note Gate H already ate the hereditary — hence semisimple — case, so here
   `A` is non-semisimple): verdict `not_tilted`, `reason="self_injective"`. Certificate:
   self-injective + non-hereditary ⇒ `gl.dim = ∞ ⇒ not tilted` (the theorem chain above). This
   is the only correct route for `kZ₃/J²`, which the knit refuses.
3. **Gate G — `gl.dim A ≥ 3` ⇒ not tilted (sound fast-path).** Else compute `gd =
   global_dimension(A)`. **The shipped contract (verified in `modules/ext.py`): `gd.value` is
   ALWAYS an `int` and `gd.exact` a `bool`; when `exact=False`, `gd.value` is a certified LOWER
   BOUND** (some simple has `pd ≥ gd.value`, unresolved within `bound=32`). So the gate is simply
   **`if gd.value >= 3: not_tilted`** — an exact `gl.dim ≥ 3` and a lower bound `≥ 3` both prove
   `gl.dim > 2`. (No `None`/`"infinite"`/symbolic-infinity case exists — the earlier draft's
   defensive `isinstance`/`None` guard was self-contradictory and is deleted.) `gd.value ≤ 2`
   never yields a verdict here — `gl.dim ≤ 2` is necessary, not sufficient (quasi-tilted algebras
   also have `gl.dim ≤ 2`; P61 separates those) — so control falls through to the search.
4. **The faithful-section search (rep-finite, non-self-injective, `gl.dim ≤ 2` or a `≤2` lower
   bound).** Knit `Γ_A` (P41 `A.ar_quiver(budget_modules=budget_modules)`). If `not
   ar.is_complete`, return the honest `unknown` verdict with `ar.status` (`"budget"` for
   rep-infinite, `"error"`) and the note that **the local criterion (arXiv:1409.2054, Thm 2.6 on
   a locally-computed cut) is the documented rep-infinite extension path** (ledger entry, Task G).
   Else run `_faithful_section_search(ar, A, budget_sections)` and certify.

### The search (`_faithful_section_search`)

Inputs from the complete knit: the indecomposable universe `U = [X₀ … X_{N−1}]` (each a
`Module`, with `name`/`dimvec`), the AR arrows `ar.arrows : {(i,j) → mult}`, the `τ`-orbits
`ar.tau_orbits` (index classes), and `n = |Q₀|`.

- **Components.** Split `U` into connected components under the **undirected** arrow graph.
- **Prune by the `#orbits == n` consequence (sound — see the completeness chain above).** Discard
  any component whose number of `τ`-orbits `≠ n`: the Theorem-2.6 witness is a slice, a slice is a
  section, a section meets each `τ`-orbit of its component exactly once and has exactly `n`
  modules; the only possible host is a "connecting component" with exactly `n` orbits. **(Cuts in
  general need not be transversals — we invoke the transversal property only for the
  slice/section witness.)**
- **Enumerate section candidates (budget-capped, exhaustive) — the `budget_sections` knob.** For
  a qualifying component, enumerate the one-module-per-`τ`-orbit **transversals** `Δ` (Cartesian
  product over the component's `n` orbits) with a hard cap **`budget_sections`** on the number of
  transversals examined. **The transversal count is the PRODUCT of the orbit sizes** — a different
  and typically much larger quantity than the module count `N` that `budget_modules` caps: e.g.
  (live-verified) `kA₃/rad²` → sizes `[3,1,1]` = **3** transversals, `rad²=0 A₅` → `[5,1,1,1,1]` =
  **5**, `kA₄`/`kA₅` (hereditary — for illustration only, Gate H short-circuits them) → `24`/`120`.
  For non-hereditary connecting components the product can grow combinatorially, hence a **separate
  cap**. The public default **`budget_sections=4096`** admits the entire tested zoo with wide
  headroom (all counts above are `≪ 4096`, reconciling the Task-A direct-search test which passes
  `budget_sections=4096`) while capping pathological products; `budget_modules=256` matches the
  P41/P55 knit default. If the cap trips before the component is exhausted, return
  `status="budget"` (an honest *unknown*, never a false `not_tilted`). **Optimization (P55 seed):**
  try `atlas.ext_projectives_right` / `ext_injectives_left` (the ACLV left/right sections) as the
  FIRST transversal — a certified hit short-circuits the enumeration. **Optimization (cut
  pre-filter):** a transversal that is not a **cut** (Def 2.1 — a cheap local arrow check using the
  `τ`/`τ⁻` links read off the orbit chains) cannot be a slice; skip it before the expensive
  filters.
- **Filter (Liu–Skowroński, Thm 2.6 — the record's operative condition).** For each transversal
  `Δ`: **faithful** `_is_faithful(A, ⊕Δ)` (`ann_A(⊕Δ) = 0`) **and** `_hom_tau_zero(A, Δ)`
  (`hom_dim(X, Y.tau()) == 0` for all `X, Y ∈ Δ` — `Y.tau()` is `0` when `Y` is projective, and
  `hom_dim(·, 0) = 0`). Both are **exact over every Domain** (no char caveat).
- **Certify (Ringel Thm 1.9(2) — the independent, airtight acceptance).** For the first `Δ`
  passing the filter, `S = ⊕Δ`; accept iff `_certify_slice(A, S)` holds:
  `is_tilting_module(S).is_tilting` **and** `is_hereditary(presented_form(end_algebra(S)))`.
  By Thm 1.9(2) this **proves** `Δ` is a slice, hence (Ringel) `A` is tilted. Return `Δ`.
- **Refute.** If every transversal of every qualifying component is examined within budget and
  none certifies, return `None` — a genuine `not_tilted` within the **rep-finite complete-knit +
  budget** scope.

Two theorems agree on every positive verdict: the **filter** is Liu–Skowroński (Thm 2.6, the
record's faithful-cut-with-`Hom(X,τY)=0`), and the **certificate** is Ringel (Thm 1.9(2),
tilting-module-with-hereditary-`End`). `is_tilting_module` already implies faithfulness, so the
explicit `_is_faithful` filter is the *cheap early prune that realizes the record's formulation*
(and the discriminator for Liu's sincere-not-faithful counterexample) — belt-and-suspenders,
never load-bearing alone.

### The reconstruction certificate + type (`_certify_slice`, `_section_quiver` + `dynkin_type`)

Given the accepted `S = ⊕Δ`:

- **`H = End_A(S)` presented.** `end_algebra(S)` is structure-constant (`quiver is None`);
  `global_dimension`/`is_hereditary` **refuse** structure-constant algebras (verified live:
  `_require_provenance` / `_require_quiver`). So we **present** `H` via P44
  `core/basic.py::presented_form(end_algebra(S))` (Gabriel-quiver recovery, `kQ/I`) and then
  `is_hereditary(H)` (`= not H.relations and H.quiver.is_acyclic()`) is the heredity certificate.
  This is Ringel's `H` **taken over `A`** (the algebra under test) — the reconstruction direction
  the record flags: `A = End_H(D(S))`, and `H = End_A(⊕Δ)`.
- **Hereditary type (via the shipped `dynkin_type`).** The **section `Σ = Δ`** — a connected,
  acyclic full subquiver of `Γ_A` — is a tree, so its type is read by the shipped
  `invariants/dynkin_type.py::dynkin_type(quiver)`, which **takes a `Quiver`** and returns a
  tuple `('A', 3)` / `('D', 4)` / `('E', 6)` / `('~A', n)` (a single cycle) or `None`. Two small
  seams (specified, self-contained):
  - **`_section_quiver(ar, section_indices) -> Quiver`** — the *section-graph-to-quiver adapter*:
    build a `Quiver` whose vertices are `section_indices` (relabeled `1..k`) and whose arrows are
    the AR arrows `(i, j) ∈ ar.arrows` with both endpoints in `section_indices` (multiplicity
    flattened to distinct arrow names). Since a section is connected and acyclic, the result is a
    tree; `dynkin_type` classifies it as `A/D/E` (rep-finite ⇒ always Dynkin).
  - **`_type_str(t) -> str`** — `('A', 3) → "A_3"`, `('~A', 2) → "~A_2"`, `None →
    "unclassified"` (on the search path a `None` would be a **loud internal error** — a valid
    slice's section is a Dynkin tree; assert it is not `None`).
  The section type is **cross-checked** against `_type_str(dynkin_type(H.quiver))` (the recovered
  Gabriel quiver of `H`): they must agree, or `_certify_slice` raises a loud internal
  `QuiverlabError` (a knit or recovery bug, never a silent mismatch). The section route is
  **char-free**; the `H`-Gabriel route inherits the P44 char caveat (see below). For **Gate H**,
  the type is `_type_str(dynkin_type(A.quiver))` directly (`A`'s own quiver).
- **Honest reconstruction wording.** `A = End_H(D(S))` is **guaranteed by Ringel's theorem**
  once `S` is a tilting module with hereditary `End_A(S)`. quiverlab has **no algebra-iso
  certifier**, so the report states this as *theorem-guaranteed* and additionally reports the
  checkable invariants `dim H`, `dim A`, Cartan-matrix agreement between `presented_form(H)`'s
  type and the section type, and (cross-engine, where cheap) a round-trip comparison of
  `dim`/Cartan of `A` against a rebuilt `End_H(D(S))`. It never claims a bare `≅` the engine
  did not check.

### `_is_faithful` (a redundant early prune) and the char scope

- **`_is_faithful(A, M)` — REDUNDANT for correctness, kept as a cheap prune.** `is_tilting_module`
  already enforces faithfulness (a **tilting module is faithful** — Bongartz), so the airtight
  `_certify_slice` acceptance never relies on this filter; Liu's sincere-not-faithful cut would
  fail `_certify_slice` regardless. `_is_faithful` exists (a) to realize the record's explicit
  Theorem-2.6 "**faithful** cut" formulation, and (b) as a cheap early prune that rejects the
  sincere-not-faithful transversal at the *filter* stage, before the expensive tilting/`End`
  checks. **Never load-bearing alone.** Spec: grep first for a shipped `annihilator`/`is_faithful`
  and reuse; **none exists at spec time** in `modules/`. Otherwise compute the **two-sided
  annihilator** `ann_A(M) = {a ∈ A : Ma = 0}` for the **right-module convention** (CLAUDE.md:
  right modules are the quiverlab default): `M` is faithful `⟺ ann_A(M) = 0 ⟺` the k-linear
  **action map** `A → ⊕_{v,w} Hom_k(M_v, M_w)` — sending each `A`-basis element (a path `p`) to
  its action block on `M` — is **injective** (`dim ann = dim A − rank`). **Do NOT hardcode the
  per-path action direction** (`M_{s(p)}→M_{t(p)}` vs the reverse): read each path's action block
  off the shipped `Module` representation (the direction is whatever `Module` uses for the right
  action). Exact over every Domain (no char caveat).
- **Char scope (load-bearing, inherited from P44/P55).** The pipeline's `Hom(X, τY)=0`,
  faithfulness, `is_selfinjective`, and `is_hereditary` checks are exact over **every** Domain.
  But the **certificate** uses `is_tilting_module` (summand count via `decompose`) and
  `presented_form` (Gabriel recovery) — both rigorous only over **char 0 or char > dim** and
  **loud** otherwise. Any module identification (the P55 seed lookup) is `is_isomorphic`,
  **QQ / char-0 decisive**, positive-only-and-raising over large `GF(p)`. Therefore **all
  batteries run over QQ**; a `GF(32003)` parity spot-check is kept **only** for `kA_n` (distinct
  dim-vectors, no `is_isomorphic` raise). Over `char ≤ dim` the recognizer inherits the loud
  `QuiverlabError` — never a silent wrong verdict.

---

## Consumers (design the API so the P55 flip + downstream P61 are served)

- **P55 (R15, must merge FIRST — see Task 0)** — Task C consumes `atlas.left_support.components`
  / `atlas.right_support.components` and certifies each factor tilted, flipping
  `test_support_components_are_tilted_PIN`. `Algebra.is_tilted()` is the exact method name that
  test calls.
- **P61 (R18, quasi-tilted/shod/laura/ada, `needs P55`)** reads `TiltedReport.verdict` /
  `.hereditary_type` when refining the ladder (tilted ⊂ quasi-tilted); P60 exposes both as
  first-class report fields and mirrors P55's `is_complete`/`status` contract so P61's
  rep-finite scope refuses in lockstep.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P41 (AR knitting), P44 (`is_tilting_module`, `bongartz_completion`, `core/basic.py::
  presented_form`/`gabriel_quiver`), P37 (`end_algebra`, `hom_dim`, `Module.tau`) are merged on
  `dev` — verified by live import at spec time.** **P55 (`left_right_parts`, `LeftRightAtlas`) is
  NOT merged** — its plan is committed (`c08fdf7`) but `src/quiverlab/modules/left_right.py` does
  not exist; **P55 must land before P60** (Wave-2 order; Task 0 STOP gate). The P55 signatures
  below are quoted from P55's **committed plan** (pending implementation), not from live code;
  the P41/P44/P37 signatures are as **verified in `dev`**. This plan consumes:
  - P41: `Algebra.ar_quiver(budget_modules=256) -> ARQuiver` (`core/algebra.py:421`); `ARQuiver`
    `.vertices` (`[{"name","dimvec","module"}]`), `.arrows` (`{(i,j): mult}`), `.tau_orbits`
    (list of index classes), `.is_complete`, `.status ∈ {"complete","budget","unsupported",
    "error"}`, `.note`. **Self-injective ⇒ `status="unsupported"`, empty vertices, NO raise**;
    rep-infinite ⇒ `"budget"`.
  - P44: `modules/tilting.py::is_tilting_module(T, n=1) -> TiltingReport` (`.is_tilting`, `.pd`,
    `.self_ext_vanishes`, `.num_summands`, `.num_vertices`, `.note`; **the count is char 0 /
    char > dim**); `core/basic.py::presented_form(A) -> Algebra` (Gabriel `kQ/I`, char 0 / char >
    dim, loud otherwise); `gabriel_quiver`, `basic_algebra`.
  - P37: `modules/endomorphism.py::end_algebra(M) -> Algebra` (structure-constant, `quiver is
    None`); `modules/hom.py::hom_dim(M, N)` / `is_isomorphic` (**QQ / char-scope-gated**);
    `modules/duality.py::tau(M)` / `Module.tau()` (`= 0` on projectives).
  - P55 (**quoted from P55's committed plan `c08fdf7`, pending implementation — Task 0 verifies
    it before P60 starts**): `modules/left_right.py::left_right_parts(A, budget=256) ->
    LeftRightAtlas` with `.left_support`/`.right_support` (`SupportAlgebra`: `.components =
    ({"vertices","algebra"}, …)`), `.ext_injectives_left`/`.ext_projectives_right` (the section
    seeds), `.is_complete`/`.status`, `._modules`. **If P55's final API names differ from its
    plan, Task 0 re-reads the merged code and this plan's P55 references are updated to match.**
  - Core/invariants (verified in `dev`): `invariants/recognizers.py::is_hereditary(A)` (**needs
    the quiver** — refuses structure-constant); `Algebra.is_selfinjective()`; `modules/ext.py::
    global_dimension(A, bound=32) -> GlobalDimension` — **`.value` is ALWAYS an `int`, `.exact` a
    `bool`; when `exact=False`, `.value` is a certified LOWER BOUND** (some simple has
    `pd ≥ value`); **needs the quiver**. The **shipped Dynkin classifier is
    `invariants/dynkin_type.py::dynkin_type(quiver) -> ('A', 3)`-style tuple | `None`** (takes a
    `Quiver`; classifies trees `E=V−1` as `A/D/E` and single cycles `E=V` as `('~A', v−1)`) — see
    the W2 adapter in the "hereditary type" foundation. `families/radical_square_zero.py::
    RadicalSquareZero`, `linear_path_algebra`.
  Branch `plan-60-tilted-recognizer` off `dev` **only after P55 has merged** (metaplan Wave 2 —
  the flip needs P55's skipped test present; Task 0 STOP gate enforces this).
- **Honest semi-decision contract (metaplan §1.3; P41/P55 loud-cap).** `tilted_check` is
  **complete** — a certified `tilted`/`not_tilted` — iff (a) a theorem gate fires (H/S/G), or
  (b) `A` is representation-finite (knit completes within `budget_modules`), non-self-injective,
  and the transversal enumeration finishes within `budget_sections`.
  A rep-infinite non-hereditary `A` (knit `"budget"`), a budget-tripped search, or a
  presentation-less `A` yields `verdict="unknown"` with `status ∈ {"budget","unsupported",
  "error"}` and the rep-infinite-extension-path note — **never a guessed verdict**.
  `TiltedReport.is_complete`/`status` mirror `ARQuiver` exactly.
- **Char scope is load-bearing (M1, inherited P44/P55).** Verdict batteries run over **QQ**; the
  `GF(32003)` parity check is used **only** for the distinct-dim-vector `kA_n` family. Over
  `char ≤ dim` the `is_tilting_module` count and `presented_form` recovery raise `QuiverlabError`
  (loud), never a silent wrong verdict.
- **No floats in `src/`.** Verdicts/reasons are strings; slices are index/name/dim-vector lists;
  memberships are Boolean/index sets. Client-side float conversion only (`docs/gui/gui.js`,
  exempt).
- **Composition is left-to-right** (`a*b` = first `a` then `b`). `_assert_comparable` guards
  every cross-module `hom_dim`. All refusals are `QuiverlabError`.
- **Plan-32 markers.** `oracle_selfcert`: Ringel Thm 1.9(2) two-way consistency (accepted `Δ` ⇒
  `⊕Δ` tilting **and** `End_A(⊕Δ)` presented hereditary; `|Δ| = n`; the section is acyclic +
  convex on the knit; the tilting-but-not-slice discriminator `A_A` of `kA₃/rad²`); the theorem
  gates' self-consistency (Gate S ⇒ `gl.dim` lower bound climbs; Gate G's lower-bound logic).
  `oracle_literature`: hereditary `kA_n`/`kD₄` ⇒ tilted with the projective section + Dynkin
  type; `kA₃/rad²` ⇒ tilted, slice `{S₂,P₂,P₃}`, `End = kA₃`, type `A₃`; `kZ₃/J²` (=
  cluster-tilted `A₃` = `Jac(3-cycle,αβγ)`) self-injective ⇒ not tilted; `rad²=0 A₅` ⇒ not
  tilted; Liu's sincere-not-faithful `rad²=0` counterexample ⇒ not tilted. `oracle_crossengine`:
  the reconstruction round-trip (`presented_form(End_A(S))` type ≡ section-graph type; and
  `dim`/Cartan of `A` vs a rebuilt `End_H(D(S))` where cheap); the search-refutation vs the
  gl.dim gate agreeing on `rad²=0 A₅`. `qpa`: QPA has **no** tilted-algebra recognizer
  (fail-if-appears probe), and the slice module is a **tilting module** crosschecked via QPA's
  tilting predicate (P44 precedent). `oracle_*` markers are FORBIDDEN in `tests/{webapp,gui,hpc}`
  (`test_oracle_classes.py`) — engine batteries in `tests/modules/`, cross-runner twins in
  `tests/webapp/` unmarked.
- **Mid-merge-train counts.** v1.0.0 lands ~29 subplans in waves; absolute suite counts drift.
  **Task G recounts the oracle-class table at merge time** by running
  `tests/release/test_oracle_classes.py` (paste live numbers, never a guessed-at-authoring
  count) and claims only this plan's deltas.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table
  green) and adds its citations to `citations/references.bib` + `registry.py` (`bibtex()`
  hard-fails if the two disagree). Conventional commits; green at every commit.

**M3 resolution (no red commits).** `TiltedReport`, `tilted_check`, the theorem gates, the
search, and the `Algebra` delegates are assembled in **Task A**, with the reconstruction fields
(`hereditary_algebra`, `reconstruction`) defaulting to `None` and `hereditary_type` populated
from the **section graph** (char-free). **Task B** populates the presented-`End` reconstruction
+ the Gabriel cross-check. Because the fields already exist (safe defaults) from Task A, every
task boundary is green — no forward reference, no `xfail` fence. Task C flips the P55 fence
(the only edit to an existing test file).

---

### Task 0: STOP gate — verify the P55 prerequisite is merged (binding, blocking)

**P55 is NOT yet merged** (verified at spec time: `src/quiverlab/modules/left_right.py` does not
exist; only P55's plan is committed at `c08fdf7`; the metaplan §6 ledger shows all of Wave 2
unimplemented). P60 consumes P55's `left_right_parts` atlas (the fast-path section seed) and, in
Task C, flips P55's skipped `test_support_components_are_tilted_PIN`. **Before any P60 work
begins, run the STOP gate:**

- [ ] **Step 1: Run the STOP gate — verbatim:**

  > Verify `from quiverlab.modules.left_right import left_right_parts` imports **AND**
  > `tests/modules/test_left_right_support.py::test_support_components_are_tilted_PIN` exists;
  > if either fails, **STOP — land P55 before P60.**

  ```bash
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -c \
    "from quiverlab.modules.left_right import left_right_parts; print('P55 import OK')"
  grep -q "def test_support_components_are_tilted_PIN" \
    tests/modules/test_left_right_support.py && echo "P55 fence present" || echo "MISSING"
  ```

  Both must succeed. (The module path `src/quiverlab/modules/left_right.py`, the function
  `left_right_parts`, and the test id `test_support_components_are_tilted_PIN` in
  `tests/modules/test_left_right_support.py` are quoted from P55's committed plan `c08fdf7`.)

- [ ] **Step 2: Reconcile the API.** If P55 merged with different final names (module path,
  `left_right_parts` signature, `LeftRightAtlas` fields `left_support`/`right_support`/
  `ext_injectives_left`/`ext_projectives_right`/`.components`, or the `is_tilted()` method name
  the fence calls), **re-read the merged P55 code** and update every P55 reference in this plan
  (the Global-Constraints consumes-list, the P55-connection section, Tasks A/C) to match the
  real API before writing code. There is **no forward-referencing tests** rule: Task C's flip
  targets the exact merged test id.

- [ ] **Step 3: Confirm P41/P44/P37 are still importable** (they were at spec time):
  ```bash
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -c \
    "from quiverlab.modules.tilting import is_tilting_module; \
     from quiverlab.core.basic import presented_form; \
     from quiverlab.modules.endomorphism import end_algebra; \
     from quiverlab.invariants.dynkin_type import dynkin_type; print('P41/P44/P37 OK')"
  ```

No commit for Task 0 (a gate, not a change). Proceed to Task A only when all three steps pass.

---

### Task A: `tilted.py` — theorem gates + faithful-section search + `TiltedReport` + delegates

**Files:**
- Create: `src/quiverlab/modules/tilted.py`
- Modify: `src/quiverlab/core/algebra.py` (thin `tilted_check`/`is_tilted` delegates beside
  `is_tilting_module`, `core/algebra.py:469`)
- Test: `tests/modules/test_tilted.py`

**Interfaces:**
- Consumes: `Algebra.ar_quiver`, `is_hereditary`, `Algebra.is_selfinjective`,
  `modules.ext.global_dimension`, `modules/hom.py::hom_dim`, `Module.tau`,
  `modules/morphism.py::direct_sum`, `modules/tilting.py::is_tilting_module`,
  `modules/endomorphism.py::end_algebra`, `core/basic.py::presented_form`,
  `left_right.py::left_right_parts` (seed, optional).
- Produces:
  ```python
  @dataclass
  class TiltedReport:
      algebra: object
      verdict: str            # "tilted" | "not_tilted" | "unknown"
      reason: str             # "hereditary"|"self_injective"|"gldim>2"|
                              #   "faithful_section_found"|"search_exhausted"|
                              #   "budget"|"unsupported"|"error"
      slice: list = ()        # [ {"index","name","dimvec"} , ... ]  (Sigma)
      slice_module: object = None       # the Module S = (+) Sigma
      hereditary_type: str = None       # "A_3" | "D_4" | ... | "~A_2" | "wild" | graph
      hereditary_algebra: object = None # H = End_A(S) presented (Task B)
      reconstruction: dict = None       # {"dim_A","dim_H","type","note"} (Task B)
      universe_size: int = 0
      is_complete: bool = False
      status: str = "error"             # mirrors ARQuiver / gate provenance
      note: str = ""
      def __bool__(self): return self.verdict == "tilted"   # so `assert A.is_tilted()` works
  def _section_quiver(ar, section_indices)            # -> Quiver (W2 adapter; relabel 1..k)
  def _type_str(t) -> str                             # ('A',3)->"A_3", ('~A',2)->"~A_2", None->...
  def _is_faithful(A, M) -> bool                      # redundant early prune (see W4)
  def _hom_tau_zero(A, section_modules) -> bool
  def _faithful_section_search(ar, A, budget_sections, seed=None) -> list[int] | None  # | "BUDGET"
  def _certify_slice(A, S) -> tuple[bool, object]     # (ok, presented H)  -- Task B fills H
  def tilted_check(A, *, budget_modules=256, budget_sections=4096) -> TiltedReport
  ```
  `Algebra.tilted_check(budget_modules=256, budget_sections=4096)` delegates (lazy-import to
  avoid the `modules → core` cycle); `Algebra.is_tilted(budget_modules=256, budget_sections=4096)
  -> bool` returns `bool(tilted_check(...))`. The type comes from the shipped
  `invariants/dynkin_type.py::dynkin_type` (a `Quiver → ('A',3)`-tuple | `None`) via
  `_section_quiver` + `_type_str` — there is **no** hand-rolled `_dynkin_type`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tilted.py
"""Tilted-algebra recognizer via the Liu-Skowronski faithful-section criterion (Plan 60 / R17,
Liu Arch. Math. 61 (1993); arXiv:1409.2054 Thm 2.6; Ringel LNM 1099 (4.2) = Liu Thm 1.9(2)).
Verdict pipeline: Gate H (hereditary => tilted), Gate S (non-semisimple self-injective => not),
Gate G (gl.dim > 2 => not), then the exhaustive faithful-section search on the knitted AR quiver.
QQ-scope (the tilting count + Gabriel recovery are char 0 / char > dim)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.tilted import tilted_check

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _radsq_a5():                                   # ACLV 2.2(b) / P55 shared fixture, gl.dim 4
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


def _kA3_radsq():                                  # P55 A_lambda shape: kA3/rad^2, TILTED, type A3
    return RadicalSquareZero(Quiver([1, 2, 3], {"a1": (2, 1), "a2": (3, 2)}), field=QQ)


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_hereditary_is_tilted(n):
    A = linear_path_algebra(n, field=QQ)           # kA_n hereditary, rep-finite
    rep = tilted_check(A)
    assert rep.verdict == "tilted" and rep.reason == "hereditary"
    assert rep.hereditary_type == f"A_{n}"
    assert bool(rep) is True and A.is_tilted() is True
    # M3: rep-finite hereditary carries a CERTIFIED projective slice (A_A), n summands.
    assert len(rep.slice) == n
    from quiverlab.modules.tilted import _certify_slice
    ok, _H = _certify_slice(A, rep.slice_module)   # End_A(A_A) hereditary => a genuine slice
    assert ok is True


@lit
def test_kA3_radsq_is_tilted_type_A3():
    A = _kA3_radsq()
    rep = tilted_check(A)
    assert rep.verdict == "tilted" and rep.reason == "faithful_section_found"
    assert len(rep.slice) == 3                      # |Sigma| = n (consequence)
    assert {r["name"] for r in rep.slice} == {"S_2", "P_2", "P_3"}
    assert rep.hereditary_type == "A_3"             # End_A(S) = kA3 (verified live)


@lit
def test_selfinjective_kZ3_J2_not_tilted():
    # kZ3/J^2 = cluster-tilted A3 (3-cycle) = Jac(3-cycle, a*b*c): self-injective, knit REFUSES;
    # the non-semisimple-self-injective theorem gate gives a REAL not_tilted, not a refusal.
    A = RadicalSquareZero(Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}), field=QQ)
    assert A.is_selfinjective() is True
    rep = tilted_check(A)
    assert rep.verdict == "not_tilted" and rep.reason == "self_injective"


@lit
def test_gldim_gt_2_not_tilted():
    A = _radsq_a5()                                 # gl.dim 4 (exact) -> Gate G
    assert A.global_dimension().value == 4
    rep = tilted_check(A)
    assert rep.verdict == "not_tilted" and rep.reason == "gldim>2"


@selfcert
def test_truncated_poly_selfinjective_not_tilted():
    A = truncated_polynomial(3, field=QQ)           # k[x]/(x^3) self-injective non-semisimple
    rep = tilted_check(A)
    assert rep.verdict == "not_tilted" and rep.reason == "self_injective"


@selfcert
def test_search_refutes_radsq_a5_directly():
    # exercise the SEARCH path (bypassing Gate G): the knittable rad^2=0 A5 (9 indec, 5 orbits,
    # NOT tilted) has NO faithful section -- _faithful_section_search returns None.
    from quiverlab.modules.tilted import _faithful_section_search
    A = _radsq_a5()
    ar = A.ar_quiver(budget_modules=128)
    assert ar.is_complete and len(ar.tau_orbits) == 5   # 5 orbits, sizes [5,1,1,1,1] => 5 transversals
    assert _faithful_section_search(ar, A, budget_sections=4096) is None


@selfcert
def test_slice_module_is_tilting_and_size_n():
    from quiverlab.modules.tilting import is_tilting_module
    rep = tilted_check(_kA3_radsq())
    assert is_tilting_module(rep.slice_module).is_tilting is True   # Ringel 1.9(2), tilting half
    assert len(rep.slice) == rep.universe_size - 2 == 3             # |Sigma| = n = 3 of 5 indec


@selfcert
def test_structure_constant_algebra_refused():
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.endomorphism import end_algebra
    A = linear_path_algebra(2, field=QQ)
    E = end_algebra(A.projective(1))                # quiver is None
    with pytest.raises(QuiverlabError):
        tilted_check(E)


@selfcert
def test_delegate_matches_free_function():
    A = _kA3_radsq()
    assert A.tilted_check().verdict == tilted_check(A).verdict
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.modules.tilted`.

  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_tilted.py -v`

- [ ] **Step 3: Implement the pipeline + search + `TiltedReport`** in `tilted.py`.

Structure (mirror P55's `left_right.py`):

```python
"""Tilted-algebra recognizer (Plan 60 / R17). Liu-Skowronski criterion (Liu Arch. Math. 61
(1993) 12-19; Liu arXiv:1409.2054 Thm 2.6): A is tilted iff Gamma_A has a faithful cut/section
Sigma with Hom_A(X, tau Y) = 0 for all X, Y in Sigma; certified by Ringel's slice theorem
(Thm 1.9(2), LNM 1099 (4.2)): Sigma is a slice iff S = (+) Sigma is a tilting A-module with
End_A(S) hereditary, and then A = End_H(D(S)). Representation-finite exhaustive scope (the knit
refuses self-injective + rep-infinite); three theorem gates (hereditary => tilted; non-semisimple
self-injective => not; gl.dim > 2 => not) extend and speed up the verdict. Exact; char 0 / char >
dim for the tilting count + Gabriel recovery (loud otherwise)."""
from __future__ import annotations
from dataclasses import dataclass, field

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.dynkin_type import dynkin_type
from quiverlab.invariants.recognizers import is_hereditary
from quiverlab.modules.hom import hom_dim
from quiverlab.modules.morphism import direct_sum


def tilted_check(A, *, budget_modules=256, budget_sections=4096):
    if getattr(A, "quiver", None) is None:
        raise QuiverlabError("tilted_check needs a quiver presentation",
                             hint="build A via Quiver.algebra(...); a tilted verdict recovers "
                                  "the hereditary type from the Gabriel quiver of End_A(S)")
    n = len(list(A.quiver.vertices))
    # Gate H
    if is_hereditary(A):
        typ = _type_str(dynkin_type(A.quiver))       # ('A',3) -> "A_3"
        slc, S = _projective_section(A, budget_modules)   # certified if rep-finite; else () / None
        note = ("hereditary: slice = the indecomposable projectives (A_A), End_A(A_A) hereditary"
                if S is not None else
                "hereditary but representation-infinite: slice omitted (the postprojective "
                "section exists but is not knit-enumerable); verdict + type stand")
        return TiltedReport(A, "tilted", "hereditary", slice=slc, slice_module=S,
                            hereditary_type=typ, is_complete=True, status="complete", note=note)
    # Gate S
    if A.is_selfinjective():
        return TiltedReport(A, "not_tilted", "self_injective", is_complete=True,
                            status="complete",
                            note="non-semisimple self-injective => gl.dim = infinity => not "
                                 "tilted (tilted => gl.dim <= 2)")
    # Gate G  --  gd.value is ALWAYS int; exact=False => it is a certified LOWER BOUND
    gd = A.global_dimension()
    if gd.value >= 3:                                 # exact >=3 OR lower bound >=3: gl.dim > 2
        return TiltedReport(A, "not_tilted", "gldim>2", is_complete=True, status="complete",
                            note=f"gl.dim {'=' if gd.exact else '>='} {gd.value} > 2; "
                                 f"tilted => gl.dim <= 2")
    # Search
    ar = A.ar_quiver(budget_modules=budget_modules)
    if not ar.is_complete:
        return TiltedReport(A, "unknown", ar.status, is_complete=False, status=ar.status,
                            note=(ar.note or "") + " | rep-infinite extension path: the local "
                                 "criterion arXiv:1409.2054 (Thm 2.6 on a locally-computed cut)")
    seed = _p55_seed(A, budget_modules)               # optional fast-path (guarded)
    found = _faithful_section_search(ar, A, budget_sections, seed=seed)
    if found == "BUDGET":
        return TiltedReport(A, "unknown", "budget", universe_size=len(ar.vertices),
                            is_complete=False, status="budget",
                            note="section (transversal) enumeration exceeded budget_sections")
    if found is None:
        return TiltedReport(A, "not_tilted", "search_exhausted", universe_size=len(ar.vertices),
                            is_complete=True, status="complete",
                            note="no faithful section with Hom(X, tau Y)=0 (rep-finite, complete "
                                 "knit, budget-exhaustive)")
    slc = [{"index": i, "name": ar.vertices[i]["name"], "dimvec": ar.vertices[i]["dimvec"]}
           for i in found]
    S = _direct_sum([ar.vertices[i]["module"] for i in found])
    typ = _type_str(dynkin_type(_section_quiver(ar, found)))
    return TiltedReport(A, "tilted", "faithful_section_found", slice=slc, slice_module=S,
                        hereditary_type=typ, universe_size=len(ar.vertices),
                        is_complete=True, status="complete")
    # hereditary_algebra / reconstruction stay None here; Task B fills them (M3).
```

`_faithful_section_search(ar, A, budget_sections, seed)`: components → prune to `#orbits == n` →
(seed first, then) `budget_sections`-capped transversal enumeration → cut pre-filter →
`_is_faithful` + `_hom_tau_zero` filter → `_certify_slice` acceptance (Task A ships the tilting
half of `_certify_slice`; Task B adds the presented-`End` half). Return the indices, `None`, or
the sentinel `"BUDGET"`.

**Adjust to reality (Task A):**
- **Gate ordering is load-bearing.** H before S (so semisimple, being hereditary, is `tilted`
  not mis-hit by S); S before G (self-injective's `gl.dim` is only a lower bound — S is the
  certain route); G before the knit (fast, and it gives verdicts for high-`gl.dim` rep-finite
  algebras the search would also refute). `gd.value >= 3` refutes on a **certified lower bound**
  too — a lower bound `≥ 3` already proves `gl.dim > 2`. Never refute on `gd.value <= 2`. **No
  `None`/`"infinite"` case exists** (`GlobalDimension.value` is always `int`) — do not add a
  defensive guard for one.
- **Gate H's certified payload (M3).** `_projective_section(A, budget_modules)`: knit
  `A.ar_quiver`; if complete (rep-finite), locate the indecomposable projectives `{P_v}` in `U`,
  form `S = ⊕P_v` (`= A_A`), run `_certify_slice(A, S)` (Task B) — it holds (`End_A(A_A)` is
  hereditary) — and return the certified `(slice_records, S)`. If the knit is **incomplete**
  (rep-infinite hereditary, e.g. Kronecker): return `([], None)` — the verdict `tilted` + type
  stand, the slice is omitted with the honest note (record-level omission, honest scope).
- **`Module.tau()` returns the zero module on projectives** (`.dim == 0`); `hom_dim(X, 0) == 0`,
  so `_hom_tau_zero` needs no projective special-case — but assert `hom_dim` accepts a
  zero-`dim` target (it does; grep `duality.tau`/`_zero_module`).
- **Type via the shipped classifier.** `dynkin_type(quiver)` returns `('A',3)`-style tuples or
  `None`; `_type_str` formats (`('~A',2) → "~A_2"`, `None → "unclassified"`). `_section_quiver`
  builds a `Quiver` from the section's AR arrows (relabel `1..k`). **No hand-rolled `_dynkin_type`
  — reuse `invariants/dynkin_type.py`.**
  For the rep-finite verdicts it is always Dynkin `A/D/E`.
- **`_projective_section`** (Gate H best-effort): if `A.ar_quiver` is complete, locate the
  indecomposable projectives in `U` and report them as the slice; if the knit is incomplete
  (rep-infinite hereditary), leave `slice=()` with the note — the `tilted` verdict + type stand
  without the module detail.
- **`gd.value >= 3` guard**: `GlobalDimension.value` may be `None`/`"infinite"` in some builds —
  normalize (`is not None and isinstance(value, int) and value >= 3`); a symbolic infinity is
  `>= 3`. Grep the `GlobalDimension` container for its infinity encoding and match it.

- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/tilted.py src/quiverlab/core/algebra.py tests/modules/test_tilted.py
git commit -m "feat(modules): tilted-algebra recognizer -- Liu-Skowronski faithful-section search + hereditary/self-injective/gl.dim theorem gates, rep-finite scope, honest unknown; kA3/rad^2 tilted type A3"
```

---

### Task B: the reconstruction certificate — presented `H = End_A(S)` + Gabriel type cross-check

**Files:**
- Modify: `src/quiverlab/modules/tilted.py` (finish `_certify_slice`; populate
  `hereditary_algebra` + `reconstruction` in `tilted_check`)
- Test: `tests/modules/test_tilted_reconstruct.py`

**Interfaces:**
- Consumes: `modules/endomorphism.py::end_algebra`, `core/basic.py::presented_form`,
  `invariants/recognizers.py::is_hereditary`, `modules/tilting.py::is_tilting_module`,
  `invariants/cartan.py::cartan_matrix` (round-trip invariants).
- Produces: `_certify_slice(A, S) -> (ok, H_presented)` — Ringel Thm 1.9(2): `ok =
  is_tilting_module(S).is_tilting and is_hereditary(presented_form(end_algebra(S)))`. On a
  `tilted` verdict, `tilted_check` sets `hereditary_algebra = H` and `reconstruction =
  {"dim_A","dim_H","type","note"}`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tilted_reconstruct.py
"""Ringel slice-theorem reconstruction (Plan 60 / R17, Liu Thm 1.9(2)): H = End_A(S) presented
via P44 is hereditary of the reported type; A = End_H(D(S)) is theorem-guaranteed and the
checkable invariants (dim, Cartan, section-graph type) are reported. The tilting-but-not-slice
discriminator: A_A of kA3/rad^2 is a tilting module whose End has a relation (NOT hereditary), so
it is NOT a slice -- Thm 1.9(2) rejects it. QQ-scope (presented_form is char 0 / char > dim)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.tilted import tilted_check, _certify_slice

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA3_radsq():
    return RadicalSquareZero(Quiver([1, 2, 3], {"a1": (2, 1), "a2": (3, 2)}), field=QQ)


@lit
def test_reconstructed_H_is_hereditary_kA3():
    rep = tilted_check(_kA3_radsq())
    H = rep.hereditary_algebra
    from quiverlab.invariants.recognizers import is_hereditary
    assert is_hereditary(H)                              # End_A(S) presented = hereditary kA3
    assert len(list(H.quiver.vertices)) == 3 and not H.relations
    assert rep.reconstruction["type"] == "A_3"
    assert rep.reconstruction["dim_H"] == H.dim


@selfcert
def test_tilting_but_not_slice_rejected():
    # A_A = S_1 (+) P_2 (+) P_3 of kA3/rad^2 IS a tilting module (pd 0) but End_A(A_A) = A has a
    # relation => NOT hereditary => NOT a slice (Ringel 1.9(2)). _certify_slice must return False.
    A = _kA3_radsq()
    from quiverlab.modules.morphism import direct_sum
    S = direct_sum(A.simple(1), A.projective(2), A.projective(3))[0]   # = A_A (S_1 = P_1)
    from quiverlab.modules.tilting import is_tilting_module
    assert is_tilting_module(S).is_tilting is True       # tilting ...
    ok, _H = _certify_slice(A, S)
    assert ok is False                                   # ... but not a slice


@xeng
def test_reconstruction_roundtrip_invariants():
    # A = End_H(D(S)) is theorem-guaranteed; check the checkable invariants agree.
    A = _kA3_radsq()
    rep = tilted_check(A)
    assert rep.reconstruction["dim_A"] == A.dim
    # section type (via _section_quiver + dynkin_type) == recovered-Gabriel type (routes agree)
    assert rep.reconstruction["type"] == rep.hereditary_type == "A_3"


@lit
def test_hereditary_reconstruction_is_self():
    A = linear_path_algebra(3, field=QQ)                 # H = A
    rep = tilted_check(A)
    assert rep.hereditary_type == "A_3"
```

- [ ] **Step 2: Run to verify failure** — `hereditary_algebra`/`reconstruction` are `None`
  (Task-A defaults), so `test_reconstructed_H_is_hereditary_kA3` fails (`None.quiver`). No
  forward reference — the attributes are real, just unpopulated (M3).

- [ ] **Step 3: Implement `_certify_slice` fully + wire the fields.** `_certify_slice(A, S)`:
  compute `is_tilting_module(S).is_tilting`; if `True`, `H = presented_form(end_algebra(S))` and
  return `(is_hereditary(H), H)`; else `(False, None)`. In `tilted_check`'s `tilted`/search
  branch, call `_certify_slice`, set `hereditary_algebra = H`, `reconstruction = {"dim_A": A.dim,
  "dim_H": H.dim, "type": rep.hereditary_type, "note": "A = End_H(D(S)) (Ringel Thm 1.9(2)); "
  "reported invariants are checkable, the isomorphism is theorem-guaranteed"}`. For Gate H,
  `H = A`, `reconstruction` self-referential.

**Adjust to reality (Task B):**
- **`_certify_slice` is the airtight acceptance** — a positive filter hit that FAILS
  `_certify_slice` is a search continuation (try the next transversal), NOT an accept. (In Task A
  the search calls the tilting-only half; wire the full `_certify_slice` here so the search's
  acceptance is Thm-1.9(2)-complete. The `kA₃/rad²` `A_A` case is exactly why: it is tilting but
  its `End` has a relation.)
- **`presented_form` char caveat** — over `char ≤ dim` it raises `QuiverlabError`; let it
  propagate as the loud whole-`tilted_check` refusal (batteries are QQ). Over QQ / GF(32003) it
  recovers the Gabriel `kQ/I`.
- **Type cross-check.** `rep.hereditary_type` (from `_type_str(dynkin_type(_section_quiver(...)))`,
  char-free) must equal `_type_str(dynkin_type(H.quiver))` (from the recovered `H`); assert
  equality inside `_certify_slice` and raise a loud internal `QuiverlabError` on divergence (never
  a silent mismatch — that would be a knit or recovery bug).

- [ ] **Step 4: Run Tasks A+B, verify pass.**
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/tilted.py tests/modules/test_tilted_reconstruct.py
git commit -m "feat(modules): tilted reconstruction certificate -- present End_A(S) via P44 (Ringel 1.9(2)), hereditary + Gabriel-type cross-check; tilting-but-not-slice discriminator"
```

---

### Task C: the P55 fence flip — support components are tilted (products-of-tilted theorem)

**Files:**
- Modify: `tests/modules/test_left_right_support.py` (flip `test_support_components_are_tilted_PIN`)
- Test: `tests/modules/test_tilted_p55_supports.py` (the ACLV-2.2(b)-scale assertion)

**Interfaces:** consumes `left_right.py::left_right_parts` (P55), `Algebra.is_tilted`.

- [ ] **Step 1: Flip P55's skipped fence** (the ONLY edit to an existing test). Replace the
  `@pytest.mark.skip(...)` decorator with `@pytest.mark.oracle_literature` and make the body a
  real assert (it already reads `comp["algebra"].is_tilted()`); add a one-line note that P60
  landed the API. Then write the fixture-scale battery:

```python
# tests/modules/test_tilted_p55_supports.py
"""P55/P60 tie: the left/right support algebras A_lambda, A_rho are PRODUCTS OF TILTED algebras
(ACLV Thm A for ada). On the shared ACLV 2.2(b) rad^2=0 A5 fixture, A_lambda (verts {1,2,3}) and
A_rho (verts {3,4,5}) are each kA3/rad^2; EACH connected component is certified tilted by the
Plan-60 recognizer (this is the content P55 PIN'd for P60). QQ-scope."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

pytestmark = pytest.mark.oracle_literature


def _radsq_a5():
    return RadicalSquareZero(Quiver([1, 2, 3, 4, 5],
        {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)}), field=QQ)


def test_hereditary_support_components_tilted():
    atlas = left_right_parts(linear_path_algebra(3, field=QQ))
    for comp in atlas.left_support.components:
        assert comp["algebra"].is_tilted() is True       # hereditary factor => tilted


def test_aclv_22b_support_components_tilted():
    atlas = left_right_parts(_radsq_a5())
    for support in (atlas.left_support, atlas.right_support):
        assert len(support.components) == 1               # each connected
        for comp in support.components:
            assert comp["algebra"].is_tilted() is True    # kA3/rad^2 factor => tilted (type A3)
```

- [ ] **Step 2: Run** — expected PASS (the recognizer certifies each `kA₃/rad²` factor tilted;
  verified live: slice `{S₂,P₂,P₃}`, `End = kA₃`). If a support component is presented with a
  vertex relabeling, `is_tilted()` is label-agnostic (it recovers the type from `End`).

- [ ] **Step 3: Commit**

```bash
git add tests/modules/test_left_right_support.py tests/modules/test_tilted_p55_supports.py
git commit -m "test(modules): flip P55 PIN -- support components A_lambda/A_rho certified tilted by the Plan-60 recognizer (products-of-tilted, ACLV Thm A)"
```

---

### Task D: literature + cross-engine oracle battery (incl. Liu's sincere-not-faithful example)

**Files:**
- Test: `tests/modules/test_tilted_oracles.py` (no `src/` change — pure oracle battery)

**Interfaces:** `tilted_check`, `_faithful_section_search`, `_is_faithful`, `Algebra.opposite`.

- [ ] **Step 1: Write the battery**

```python
# tests/modules/test_tilted_oracles.py
"""Literature + cross-engine oracles for the tilted recognizer (Plan 60 / R17). Literature:
hereditary Dynkin => tilted with the projective slice + correct type (kA_n, kD_4); kA3/rad^2 =>
tilted type A3; kZ3/J^2 (= cluster-tilted A3) => not tilted; Liu's rad^2=0 sincere-not-faithful
example => not tilted (faithfulness is essential -- arXiv:1409.2054 Ex. after Thm 2.6). Cross-
engine: the gl.dim gate and the direct section search agree on rad^2=0 A5; the reconstructed
Gabriel type equals the section-graph type. QQ-scope; GF(32003) parity for kA_n only."""
import pytest

from quiverlab import GF, Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.families.dynkin import ...            # kD_4 builder (grep: dynkin_path_algebra?)
from quiverlab.modules.tilted import tilted_check, _faithful_section_search, _is_faithful

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


@lit
def test_kD4_hereditary_tilted_type_D4():
    A = ...                                           # kD_4 (central vertex + 3 leaves), acyclic
    rep = tilted_check(A)
    assert rep.verdict == "tilted" and rep.hereditary_type == "D_4"


@lit
def test_kA3_parity_gf32003_tilted():
    # kA_n has distinct dim-vectors per indec -> is_isomorphic never enters the positive-only
    # branch -> GF(32003) is safe (M1). Verdict byte-identical to QQ.
    A = linear_path_algebra(3, field=GF(32003))
    assert tilted_check(A).verdict == "tilted"


@lit
def test_liu_sincere_not_faithful_counterexample():
    # arXiv:1409.2054 Ex. after Thm 2.6: a rad^2=0 3-vertex algebra with a SINCERE cut
    # P_b - S_b - P_a s.t. Hom(X,tau Y)=0, but NOT faithful => NOT tilted. Faithfulness cannot be
    # weakened to sincereness. PIN: transcribe the exact quiver from the PDF (p.9) and verify it
    # is rep-finite, non-self-injective, and that the sincere cut fails _is_faithful.
    A = RadicalSquareZero(Quiver([...], {...}), field=QQ)    # PIN exact orientation
    assert A.is_selfinjective() is False
    rep = tilted_check(A)
    assert rep.verdict == "not_tilted"


@xeng
def test_gate_and_search_agree_on_radsq_a5():
    A = RadicalSquareZero(Quiver([1, 2, 3, 4, 5],
        {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)}), field=QQ)
    # Gate G verdict:
    assert tilted_check(A).verdict == "not_tilted"
    # direct search verdict (bypassing gates):
    ar = A.ar_quiver(budget_modules=128)
    assert _faithful_section_search(ar, A, budget=4096) is None
```

- [ ] **Step 2: Pin Liu's counterexample.** Transcribe the exact `rad²=0` 3-vertex quiver from
  the arXiv:1409.2054 figure (p. 9) into the test; verify live: `is_selfinjective() is False`,
  the knit completes, the sincere cut `P_b—S_b—P_a` exists with `Hom(X, τY)=0` (so the Hom
  filter passes) but `_is_faithful(A, ⊕cut) is False` (so the section is rejected) and
  `tilted_check` returns `not_tilted`. If the exact orientation is ambiguous, try the few
  `rad²=0` orientations of the 3-vertex `a`-source `b—c` quiver and keep the one matching the
  paper's AR-quiver picture (`P_a, P_b, P_c` projectives, `S_a, S_b, S_c` simples, 6 indec).

- [ ] **Step 3: Commit**

```bash
git add tests/modules/test_tilted_oracles.py
git commit -m "test(modules): tilted recognizer oracles -- kD4 type D4, kA_n GF parity, Liu sincere-not-faithful counterexample, gate/search agreement on rad^2=0 A5"
```

---

### Task E: GUI/webapp/report — the `tilted_check` compute kind

`tilted_check` is an **ALGEBRA-level compute kind** (computed on the drawn `A`, like
`left_right_parts` / `ar_quiver` — routed through `_dispatch`, NOT `_dispatch_module`; **schema
v1**, no module block; sizes on `A.dim`). **Budget payload (W3):** the GUI/webapp compute string
carries the user-facing **knit** budget only — `compute.push("tilted_check:256")` maps to
`budget_modules`; **`budget_sections` is an internal knob with the documented default `4096`**,
exposed only via the Python API (a section/transversal cap is not a knit size and would confuse
the picker). The parse interceptors and `tilted_check_block` set
`budget_sections=4096` internally.

**Files:**
- Modify: `src/quiverlab/modules/tilted.py` (`tilted_check_block(A, budget_modules, budget_sections=4096)`)
- Modify: `src/quiverlab/hpc/spec.py` (`parse_compute_item` budget interceptor; `_dispatch`
  branch; `_snip` recipe)
- Modify: `webapp/server/schema.py` (parse twin) + `webapp/server/estimator.py` (`_max_degree`
  exclusion tuple: add `"tilted_check"`)
- Modify: `docs/gui/runner.py` (Pyodide twin: `_parse_compute` interceptor, `compute_one`
  branch calling the SAME `tilted_check_block`, `calls` snippet, `ETA_MODEL["scalars"]`)
- Modify: `docs/gui/gui.js` (+ re-copy byte-identical to `webapp/static/gui/gui.js`): checkbox +
  budget picker, `el` id list, `buildRequest` push, `renderBlock` branch `renderTiltedCheck`,
  `KIND_CTRL`, `THEMES` (`"structure"`), `SEARCH_INDEX`
- Modify: `webapp/server/i18n/{en,es,fr,zh}.json` (`inv.tilted_check` ×4 + any
  `block.tilted_check.*` strings the renderer reads — the `test_i18n.py` key-parity + no-
  placeholder gates require ALL FOUR)
- Modify: `src/quiverlab/trace/results_html.py` (`_HEADINGS["tilted_check"] = "Tilted-algebra
  check"` + a `_block_html` branch)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py` (ONE new
  fixture `tilted_check_kA3`; existing entries byte-identical; docstring bullet)
- Test: `tests/webapp/test_tilted_check_p60.py`

**Block shape (returned identically by both runners — the byte-parity contract):**
```python
{"kind": "tilted_check",
 "n": int,                                  # |Q_0|
 "complete": bool, "status": "complete"|"budget"|"unsupported"|"error",
 "verdict": "tilted"|"not_tilted"|"unknown",
 "reason": "hereditary"|"self_injective"|"gldim>2"|"faithful_section_found"|
           "search_exhausted"|"budget"|"unsupported"|"error",
 "slice": [ {"name": "S_2"|..., "dimvec": {...}} , ... ] | [],
 "slice_dimvec": {...} | None,              # dim-vector of S = (+) Sigma
 "hereditary_type": "A_3"|... | None,
 "hereditary_algebra": {"vertices":[...], "arrows":{...}, "relations":[...], "dim": int} | None,
 "reconstruction": {"dim_A": int, "dim_H": int, "type": str, "note": str} | None,
 "note": str | None,
 "references": ["liu_tilted_1993", "liu_another_2014", "happel_ringel_tilted", "assem_book"],
 "citations": [...]}
```

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir; copy the
  `tests/webapp/test_left_right_parts_p55.py` runner-pair fixture):

```python
# tests/webapp/test_tilted_check_p60.py
"""The tilted_check algebra-level compute kind: served by hpc.spec, mirrored byte-for-byte by
the Pyodide twin. kA3 -> tilted (hereditary, type A3); kA3/rad^2 -> tilted (faithful_section,
slice S_2/P_2/P_3, type A3); kZ3/J^2 -> not tilted (self_injective)."""


def test_tilted_check_hereditary_block(tmp_path):
    from quiverlab.hpc.spec import ComputeRequest, run
    out = run(ComputeRequest.model_validate(_tc_request(quiver=("kA3",), budget=256)), tmp_path)
    b = out["results"]["tilted_check"]
    assert b["complete"] and b["verdict"] == "tilted" and b["reason"] == "hereditary"
    assert b["hereditary_type"] == "A_3"


def test_tilted_check_faithful_section_block(tmp_path):
    # kA3/rad^2: verdict tilted, reason faithful_section_found, slice named, type A3, H recovered.
    ...


def test_tilted_check_not_tilted_selfinjective(tmp_path):
    # kZ3/J^2 (self-injective): verdict not_tilted, reason self_injective, complete True.
    ...


def test_twin_parity(tmp_path):
    # same request through docs/gui/runner.py; json.dumps(sort_keys=True) equality on the block.
    ...
```

- [ ] **Step 2: Implement `tilted_check_block(A, budget_modules, budget_sections=4096)`** — call
  `tilted_check(A, budget_modules=budget_modules, budget_sections=budget_sections)`, serialize the
  `TiltedReport` (strip `slice_module`/`hereditary_algebra`
  Module/Algebra objects to `slice`/`slice_dimvec`/`hereditary_algebra` dict summaries), stamp
  `references = ["liu_tilted_1993", "liu_another_2014", "happel_ringel_tilted", "assem_book"]`.
  On the `unknown`/refused path return the same shape with empty `slice` + `status`/`note` (the
  honest refusal, like `ar_quiver_block`). Catch the `presented_form`/`is_isomorphic` char-scope
  `QuiverlabError` into `{"kind": "tilted_check", "error": "<loud message>"}` (Plan-30 per-entry
  precedent — never a 500).

- [ ] **Step 3: Wire the compute kind** (the algebra-level, `left_right_parts`-style touchpoints,
  mirroring P55 Task E exactly):
  - `spec.py::parse_compute_item`: budget interceptor `if s == "tilted_check" or
    s.startswith("tilted_check:"): ... return ComputeItem("tilted_check", None, int(b) if b
    else None)`.
  - `spec.py::_dispatch` (before the catch-all): `if kind == "tilted_check": from
    quiverlab.modules.tilted import tilted_check_block; block = tilted_check_block(A,
    budget_modules=(item.hi if item.hi is not None else 256)); block["citations"] =
    _citation_pairs(block["references"]); return block, None`. (`budget_sections` keeps its
    internal default `4096`.)
  - `spec.py::_snip`: `"tilted_check": lambda it: (f"A.tilted_check(budget_modules="
    f"{it.hi if it.hi is not None else 256})")`.
  - `webapp/server/schema.py::parse_compute_item`: the same interceptor.
    `webapp/server/estimator.py::_max_degree`: add `"tilted_check"` to the skip tuple beside
    `"left_right_parts"`/`"ar_quiver"` (its budget is not a homological degree).
  - `docs/gui/runner.py`: `_parse_compute` interceptor, `compute_one` branch calling the SAME
    `tilted_check_block`, `calls` snippet, `ETA_MODEL["scalars"]["tilted_check"]`.
  - GUI (`docs/gui/gui.js`, then re-copy to `webapp/static/gui/gui.js` — byte-gate
    `test_draw_page.py`): checkbox + budget picker (mirror `left_right_parts`), `el` id entries
    `"tilted_check", "tilted_check-budget"`, the `buildRequest` push, a `renderBlock` branch
    `renderTiltedCheck(div, b)` (a verdict badge `tilted`/`not tilted`/`unknown` + the reason +
    the type + the slice chips `S_v/P_v/I_v` + the recovered `H` summary + `citesLine(b)`),
    `KIND_CTRL {cb, budget:true}`, `THEMES` add to `"structure"`, `SEARCH_INDEX` entry with
    multilingual keywords + an `example`.
  - i18n ×4 (`en/es/fr/zh.json`): `"inv.tilted_check"` in ALL FOUR + any `block.tilted_check.*`
    labels (verdict/reason glosses) the renderer reads.
  - `results_html.py`: `_HEADINGS["tilted_check"] = "Tilted-algebra check"` + a `_block_html`
    branch `_tilted_check_html(b)` (verdict + reason gloss + slice list + type + recovered `H`;
    reuse the `_derived_fingerprint_html` table idiom).

- [ ] **Step 4: Add ONE golden** `tilted_check_kA3` to `_runner_goldens.json` (full kA3 run:
  verdict tilted, reason hereditary, type A_3). Verify existing entries byte-identical BEFORE
  appending (a schema-v1 budget kind adds no request field — canonical keys untouched); add a
  docstring bullet to `test_runner_delegation.py`.

- [ ] **Step 5: Run the gates**

  `... -m pytest tests/webapp/test_tilted_check_p60.py tests/webapp/test_runner_delegation.py
  tests/webapp/test_i18n.py tests/webapp/test_draw_page.py tests/hpc -q`
  Expected PASS (both runners byte-identical; self-injective → verdict not_tilted; gui.js twins
  byte-identical; i18n key-parity green).

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc,trace): tilted_check compute kind -- verdict + slice + hereditary type + reconstruction, both runners, i18n x4, one golden"
```

---

### Task F: QPA cross-oracle (probe-first) + tilting-module crosscheck

**Files:**
- Test: `tests/qpa/test_tilted_qpa.py`

**Interfaces:** the QPA session (`session.should_skip_qpa()`, `session.libgap_handle()`) — mirror
the Plan-35/P55 probe. **QPA has no tilted-algebra recognizer** — but its tilting-module
predicate corroborates the slice module (P44 precedent, `tests/qpa/` tilting crosschecks).

- [ ] **Step 1: Probe live QPA** (the fail-if-appears pattern):

```python
# tests/qpa/test_tilted_qpa.py
"""QPA cross-check for the tilted recognizer (Plan 60). qpa-marked: skips locally, mandatory
under QUIVERLAB_REQUIRE_QPA=1. QPA has NO tilted-algebra recognizer -- this test documents that
(fail-if-appears) and corroborates the slice module: on kA3/rad^2, tilted_check's slice module is
a TILTING module, which QPA's tilting predicate confirms."""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_has_no_tilted_recognizer():
    lg = session.libgap_handle()
    for name in ("IsTiltedAlgebra", "IsTilted", "TiltedAlgebra", "SliceInModuleCategory"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), (
            f"QPA now exposes {name} -- wire a real crosscheck (this assert is the trip-wire "
            "that Plan 60's QPA scope note is stale)")


def test_slice_module_is_tilting_via_qpa():
    from quiverlab import Quiver, RadicalSquareZero
    from quiverlab.fields import QQ
    A = RadicalSquareZero(Quiver([1, 2, 3], {"a1": (2, 1), "a2": (3, 2)}), field=QQ)
    rep = A.tilted_check()
    assert rep.verdict == "tilted"
    A.crosscheck("is_tilting_module", rep.slice_module).assert_agree()   # P44 QPA verb (grep)
```

- [ ] **Step 2:** Confirm the tilting crosscheck verb name/signature in `qpa/crosscheck.py`
  (P44 shipped a QPA tilting comparison — grep before wiring; if it rides a different verb,
  adapt). If no clean two-argument tilting crosscheck exists, the probe
  (`test_qpa_has_no_tilted_recognizer`) is the deliverable on its own and the tilting
  corroboration falls back to the P44 QPA tilting battery — honest, never a silent skip.

- [ ] **Step 3: Run live** `... -m pytest tests/qpa/test_tilted_qpa.py -v` (the venv has
  `[qpa]`). Expected PASS live (probe absent; slice module confirmed tilting).

- [ ] **Step 4: Commit**

```bash
git add tests/qpa/test_tilted_qpa.py
git commit -m "test(qpa): tilted recognizer probe -- QPA has no tilted-algebra verb (fail-if-appears); slice module confirmed tilting via QPA"
```

---

### Task G: verification page, citations, README, metaplan, ledger, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick P60 in the §6 ledger)
- Modify: `docs/plans/DEEPER-ENGINES-BACKLOG.md` (the rep-infinite-extension ledger entry)
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`)

- [ ] **Step 1: Citations** (VERIFIED BibTeX only — `_r(key, bibtex_key, kind, title,
  annotation, *tags)`, `registry.py`; `bibtex()` hard-fails if `.bib` and registry disagree).
  Add to `references.bib`:

```bibtex
@article{Liu1993,
  author  = {Liu, Shiping},
  title   = {Tilted algebras and generalized standard {A}uslander--{R}eiten components},
  journal = {Archiv der Mathematik},
  volume  = {61},
  number  = {1},
  pages   = {12--19},
  year    = {1993},
  doi     = {10.1007/BF01258050},
}
@misc{Liu2014,
  author        = {Liu, Shiping},
  title         = {Another characterization of tilted algebras},
  year          = {2014},
  eprint        = {1409.2054},
  archivePrefix = {arXiv},
  primaryClass  = {math.RT},
}
@article{HappelRingel1982,
  author  = {Happel, Dieter and Ringel, Claus Michael},
  title   = {Tilted algebras},
  journal = {Transactions of the American Mathematical Society},
  volume  = {274},
  number  = {2},
  pages   = {399--443},
  year    = {1982},
}
```

  and in `registry.py`:

```python
_r("liu_tilted_1993", "Liu1993", "recognizer",
   "Tilted algebras and generalized standard Auslander-Reiten components",
   "Liu (independently with Skowronski): A is tilted iff Gamma_A has a faithful generalized-"
   "standard component with a section -- the criterion Plan 60 searches for. Venue verified.",
   "recognizer"),
_r("liu_another_2014", "Liu2014", "recognizer",
   "Another characterization of tilted algebras (arXiv:1409.2054)",
   "Liu: A is tilted iff Gamma_A contains a FAITHFUL CUT Delta with Hom(X, tau Y)=0 (Thm 2.6); "
   "the cut is a finite/local object (weakly convex), the documented rep-infinite extension "
   "path. Ringel's slice theorem (Thm 1.9(2)) is the reconstruction certificate Plan 60 uses.",
   "recognizer"),
_r("happel_ringel_tilted", "HappelRingel1982", "foundation",
   "Tilted algebras",
   "Happel-Ringel: the origin of tilted algebras A = End_H(T) (H hereditary, T tilting); "
   "tilted => gl.dim <= 2 -- the theorem gate Plan 60 uses to refute high-gl.dim algebras.",
   "recognizer"),
```

  **Spec-ambiguity resolution (recorded):** BibTeX keys `Liu1993`/`Liu2014`/`HappelRingel1982`;
  snake-case registry keys `liu_tilted_1993`/`liu_another_2014`/`happel_ringel_tilted` (house
  convention). `Liu2014` ships as the **verified arXiv `@misc`** (published-venue fields not
  verified — worker upgrades to `@article` iff BibTeX-verifiable at merge; no fabricated fields).
  Liu1993 issue-no. and Happel-Ringel vol/pages `# PIN`'d for merge-time confirmation. The Ringel
  slice theorem's LNM 1099 and Skowroński 1994 are folded into the annotations unless
  independently BibTeX-verifiable at merge (then add as `@book`/`@article`). `assem_book`
  (ASS2006, tilted-algebra chapters VI/VIII) is already shipped and cited by the block.

- [ ] **Step 2: Verification page.** Add the Plan-60 subsystem rows:
  - `modules/tilted.py` (verdict pipeline) — `oracle_literature`: hereditary `kA_n`/`kD₄` ⇒
    tilted (projective slice, Dynkin type); **`kA₃/rad²` ⇒ tilted, slice `{S₂,P₂,P₃}`,
    `End=kA₃`, type `A₃`**; **`kZ₃/J²` (= cluster-tilted `A₃` = `Jac(3-cycle,αβγ)`) self-injective
    ⇒ not tilted**; `rad²=0 A₅` (gl.dim 4) ⇒ not tilted; **Liu's `rad²=0` sincere-not-faithful
    example ⇒ not tilted (faithfulness essential, arXiv:1409.2054)**. `oracle_selfcert`: Ringel
    1.9(2) consistency (`⊕Σ` tilting **and** presented `End` hereditary; `|Σ|=n`; the tilting-
    but-not-slice `A_A` discriminator); the theorem gates. `oracle_crossengine`: reconstruction
    Gabriel-type ≡ section-graph type; gl.dim gate ≡ direct search on `rad²=0 A₅`.
  - Add the **honest-scope entries**: (a) **representation-finite + non-self-injective search
    scope** — the P41 knit refuses self-injective (`unsupported`) and rep-infinite (`budget`);
    the theorem gates (H/S/G) extend verdicts to hereditary (incl. rep-infinite) and self-
    injective and high-gl.dim, but a **rep-infinite non-hereditary `gl.dim-2` tilted algebra is
    an honest `unknown`** — the local criterion arXiv:1409.2054 (a finite/local cut) is the
    documented, **not-yet-implemented** extension path (ledger entry, Step 4); **rep-infinite
    hereditary is `tilted` (Gate H) but its slice is OMITTED at the record level** (the
    postprojective section is not knit-enumerable — stated, not faked); (b) the recognizer
    needs a **quiver presentation** (structure-constants-only ⇒ loud refusal); (c) **char scope**
    — the `Hom(X,τY)=0`/faithfulness/self-injective/hereditary checks are char-free, but the
    `is_tilting_module` count and `presented_form` Gabriel recovery are **char 0 / char > dim**
    (loud otherwise); batteries QQ, `GF(32003)` parity only for `kA_n`; (d) **`A = End_H(D(S))`
    is theorem-guaranteed** (Ringel), the report states the checkable invariants (dim, Cartan,
    type) — no bare algebra-`≅` the engine cannot certify; (e) **QPA cannot compare** the tilted-
    algebra recognizer (fail-if-appears probe); the slice module IS confirmed a tilting module by
    QPA. Recount the class table (`tests/release/test_oracle_classes.py` drives the numbers —
    run collection for each of `oracle_literature`/`oracle_crossengine`/`oracle_selfcert`/`qpa` +
    the union, paste the LIVE counts, re-run to green; mid-merge-train honest).

- [ ] **Step 3: README.** One features line: "tilted-algebra recognizer (Liu–Skowroński): decides
  whether `A ≅ End_H(T)` (`H` hereditary, `T` tilting) by a faithful-section search on the AR
  quiver, returning the slice, the reconstructed hereditary type, and a Ringel slice-theorem
  certificate — no-code in the browser."

- [ ] **Step 4: Ledger entry** in `DEEPER-ENGINES-BACKLOG.md`: the **rep-infinite tilted
  recognizer via the local cut criterion (arXiv:1409.2054 Thm 2.6 / 2.7)** — a locally-computed
  faithful cut in a bounded AR-neighborhood (no full knit), plus Thm 2.7's **tilted quotient
  algebras** `B = A/ann(Δ)` — deferred from P60 v1; P80 reconciles.

- [ ] **Step 5: Metaplan.** Tick `P60` in the §6 progress ledger.

- [ ] **Step 6: Full gate:**
  `... -m pytest tests/modules/test_tilted*.py tests/modules/test_left_right_support.py -q`
  (deep — the touched files), `... -m pytest -q -m fast`,
  `... -m pytest tests/webapp tests/gui -q`, `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q` — all green.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-60 tilted-recognizer oracle rows + Liu/Happel-Ringel citations + honest scope (rep-finite search, QQ certificate, theorem-guaranteed reconstruction, no-QPA) + rep-infinite ledger + recounted classes"
```

---

## Acceptance (Plan-60 definition of done)

0. **Task 0 STOP gate passed:** P55 is merged (`from quiverlab.modules.left_right import
   left_right_parts` imports and `tests/modules/test_left_right_support.py::
   test_support_components_are_tilted_PIN` exists), and any P55 API-name drift has been
   reconciled into this plan's references before code was written.
1. `tilted_check`, `TiltedReport`, `tilted_check_block` public in `src/quiverlab/modules/
   tilted.py`; `Algebra.tilted_check(budget_modules=256, budget_sections=4096)` and
   `Algebra.is_tilted(budget_modules=256, budget_sections=4096)` delegates named (two budgets:
   knit cap vs transversal cap); every verdict self-certified (a theorem gate or the
   Ringel-1.9(2) slice certificate) or loudly refusing.
2. The verdict pipeline is exactly: **Gate H** (hereditary ⇒ tilted, incl. rep-infinite, type =
   underlying graph of `Q`); **Gate S** (non-semisimple self-injective ⇒ not tilted); **Gate G**
   (`gl.dim ≥ 3`, exact or certified lower bound ⇒ not tilted); then the **faithful-section
   search** on the complete knit (Liu–Skowroński Thm 2.6 filter: faithful + `Hom(X,τY)=0`) with
   the **Ringel Thm 1.9(2) certificate** (`⊕Σ` a tilting `A`-module **and** `End_A(⊕Σ)`
   presented hereditary). `|Σ| = n` is a stated consequence used as the search pruning bound,
   never as the definition.
3. **Oracles pinned:** hereditary `kA_n`/`kD₄` ⇒ tilted with the correct Dynkin type;
   **`kA₃/rad²` ⇒ tilted, slice `{S₂,P₂,P₃}`, `End=kA₃` hereditary, type `A₃`** (verified live);
   **`kZ₃/J²` (= cluster-tilted `A₃` = the 3-cycle Jacobian) self-injective ⇒ not tilted** via
   Gate S; `rad²=0 A₅` ⇒ not tilted (Gate G and, directly, `_faithful_section_search` returns
   `None`); **Liu's `rad²=0` sincere-not-faithful example ⇒ not tilted** (faithfulness cannot be
   weakened to sincereness). The record's "cluster-tilted `A₃`" is documented as the **same
   self-injective algebra** as `kZ₃/J²` (a settled correction).
4. The **reconstruction certificate**: `H = End_A(⊕Σ)` presented via P44 `presented_form` is
   hereditary of the reported type; the section-graph type ≡ the recovered-Gabriel type (loud on
   divergence); `A = End_H(D(S))` reported as theorem-guaranteed with checkable invariants (dim,
   Cartan, type); the **tilting-but-not-slice** discriminator (`A_A` of `kA₃/rad²` is a tilting
   module whose `End` has a relation ⇒ rejected) is pinned.
5. **The P55 fence flipped:** `Algebra.is_tilted()` real; P55's skipped
   `test_support_components_are_tilted_PIN` is a live assert, extended so each connected component
   of the ACLV-2.2(b) `rad²=0 A₅` support algebras `A_λ`/`A_ρ` (each `kA₃/rad²`) is certified
   tilted (the "products of tilted algebras" content, ACLV Thm A).
6. The **honest semi-decision contract**: a certified verdict iff a gate fires or the rep-finite
   non-self-injective search finishes within budget; rep-infinite non-hereditary / budget-tripped
   / presentation-less ⇒ `verdict="unknown"` with `status ∈ {"budget","unsupported","error"}` and
   the arXiv:1409.2054 rep-infinite-extension note — never a guessed verdict. `TiltedReport.
   is_complete`/`status` mirror `ARQuiver`.
7. `tilted_check` clickable end-to-end (GUI canvas → verdict/slice/type/certificate block →
   report), algebra-level compute kind (schema v1, NO module block), both runners byte-identical
   via the shared `tilted_check_block`, EN+ES+FR+ZH i18n, one golden with a documented change-log
   entry.
8. QPA probe green (`-m qpa`): no tilted-algebra recognizer (fail-if-appears); the slice module
   confirmed a tilting module via QPA's tilting predicate.
9. **Char scope recorded:** verdict batteries QQ; the tilting count + Gabriel recovery are char 0
   / char > dim (loud otherwise); `GF(32003)` parity only for `kA_n` (distinct dim-vectors).
10. `docs/verification.md` recounted (live numbers, mid-merge-train honest); the three citations
    added and BibTeX-verified (Liu2014 as the verified arXiv entry; Liu1993 issue + Happel-Ringel
    vol/pages `# PIN`'d; Skowroński/Ringel-LNM folded into annotations unless verifiable);
    README line; metaplan P60 tick; the rep-infinite-extension ledger entry; deep
    (`test_tilted*`) + fast + webapp/gui + qpa + release + citations suites green.

---

## Change-log

- **2026-08-07** — plan authored (P60 / R17). Reference re-verification done at spec time:
  Liu Arch. Math. 61 (1993) 12–19 (venue verified); arXiv:1409.2054 read pages 1–9 (Def of
  section/faithful/slice/cut, Thm 1.9(2) reconstruction, Prop 2.5, Thm 2.6 main criterion,
  Thm 2.7 tilted quotients, and the sincere-not-faithful counterexample — all transcribed
  verbatim). Live-verified in the venv (QQ): `kA₃/rad²` is tilted with slice `{S₂,P₂,P₃}`,
  `End_A(S) = kA₃` (via `presented_form`), type `A₃`; `A_A` of `kA₃/rad²` is tilting but its
  `End` has a relation (not a slice — Thm 1.9(2) discriminator); `kZ₃/J²` (= the 3-cycle
  cluster-tilted / Jacobian algebra) is self-injective (Gate S); `rad²=0 A₅` has gl.dim 4
  (Gate G) and knits to 9 indec / 5 orbits (search-refutation seam); `global_dimension` and
  `is_hereditary` refuse structure-constant algebras, so End-heredity routes through P44
  `presented_form`. Settled correction: cluster-tilted `A₃` (3-cycle) IS `kZ₃/J²` — both oracles
  are one self-injective algebra, refuted by the same theorem gate. **What was live-verified:**
  P41/P44/P37 were confirmed by live import (`is_tilting_module`, `presented_form`, `end_algebra`,
  `ar_quiver`, `dynkin_type`, `global_dimension`); **P55 was read only as its committed plan
  doc (`c08fdf7`) — `modules/left_right.py` does not yet exist and was NOT live-verified.** The
  transversal counts were measured live (kA₃/rad²=3, rad²=0 A₅=5, kA₄=24, kA₅=120), and
  `dynkin_type(quiver)` was confirmed to return tuples like `('A',3)`.
- **2026-08-07 adversarial review: 2 blocking + 2 majors + 5 minors applied** (P55 dependency
  honesty + Task-0 STOP gate; the slice⇒section⇒transversal completeness chain cited correctly;
  budget split `budget_modules`/`budget_sections`; Gate-H certified payload; GlobalDimension
  contract corrected; `dynkin_type` adapter spec'd).
