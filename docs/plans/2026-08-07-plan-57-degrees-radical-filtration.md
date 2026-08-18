# Plan 57: Liu degree theory + AR-component invariants + the radical filtration of mod A (Liu–Chaio program) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Two records, one machine. **R37 (Marco's personal request — the
headline):** the radical filtration of `mod A` — the layer dimensions
`dim rad^n(X, Y)` for indecomposables `X, Y`, the **nilpotency index** of
`rad(mod A)` for representation-finite `A`, and the certificate
`rad^∞(mod A) = 0 ⇔ A representation-finite` (Auslander) as an operational
rep-finiteness gate dual to the P41 knitting bound. **R21 (the degree theory
that computes with those layers):** Liu's **left/right degrees** of irreducible
maps, sectional paths, the postprojective/preinjective/regular partition,
directing modules, the **representation-directed** recognizer (`Γ_A` acyclic),
and generalized-standard flags — all finite sweeps on the knitted AR quiver.
Everything is **exact linear algebra on the knitted category** (route (ii)
below); everything is either **fully certified on the rep-finite domain** or
**honestly refused / window-labeled** off it — never an overclaim. Two no-code
GUI compute kinds (`radical_filtration`, `ar_invariants`) put the whole surface
one click away.

**Architecture:** Two thin exact layers over the P41 AR machinery
(`modules/ar.py`: `knit_ar_quiver`/`ARQuiver`, `irreducible_maps`,
`ar_quiver_block`) and the P37 Hom/morphism glue (`hom_space`/`hom_dim`/
`is_isomorphic`/`hom_basis`, `ModuleHom.then`, `_rad_end_basis`, `_vec`/`_unvec`):

- **`src/quiverlab/modules/radical.py`** (new) — the R37 core. `radical_filtration(A)`
  → a `RadicalFiltration` object computed **exactly** by linear algebra on the
  knitted indecomposable universe: `rad(X_i,X_j)` as a coordinate subspace of
  `Hom(X_i,X_j)` (all of it when `i≠j`; `rad End(X_i)` when `i=j`), and
  `rad^{n+1}(X,Y) = Σ_Z rad^n(X,Z)∘rad(Z,X_j)` by exact matrix composition. It
  reports per-pair layer sequences `dim rad^n(X,Y)`, the **nilpotency index**
  `N = least n with rad^n ≡ 0`, and the `rad^∞ = 0` certificate (rep-finite ⇒ the
  finite `N` is the witness). Mirrors P41's `ARQuiver` semi-decision contract:
  `is_complete`/`status` — complete iff the knit closed (rep-finite); on a
  budget/self-injective knit it returns **window-restricted** layers on the
  discovered subcategory, loudly labeled as NOT `rad(mod A)`.
- **`src/quiverlab/modules/ar_invariants.py`** (new) — the R21 sweeps. Consumes
  the `RadicalFiltration` (degrees are defined by the radical powers):
  `left_degree`/`right_degree` (Liu, decidable finite sweep), `sectional_paths`,
  the `postprojective/preinjective/regular` partition, `directing_modules`,
  `is_representation_directed` (`Γ_A` acyclic), `is_generalized_standard`. Returns
  a frozen `ARInvariants` dataclass.
- **GUI (two new algebra-block kinds)** — `radical_filtration` and `ar_invariants`
  are **algebra kinds** (they consume the whole knitted category, exactly like the
  existing `ar_quiver` kind), dispatched in `hpc/spec.py::_dispatch` and its
  Pyodide twin `docs/gui/runner.py::compute_one`, both calling ONE shared library
  builder (`radical_filtration_block`/`ar_invariants_block`) so byte-parity holds
  by construction; rendered in `trace/results_html.py`; four-locale i18n; stable
  canonical keys.

No new math engines; every number is exact (`int`/`dict`, no floats). Every
result either **self-certifies** (`rad^{n+1} ⊆ rad^n`, `rad^N = 0` and
`rad^{N-1} ≠ 0`, composition closure, the degree-witness membership triple) or
**refuses loudly** (self-injective / budget / char-scope). The primary oracles
are the mesh combinatorics on `kA_n`/`D_4`, the Chaio-school nilpotency-index
formulas, and the internal degree-vs-layer consistency; QPA has **no**
module-category-radical or degree surface (probed fail-if-appears; the layer-1
dims `dim rad(X,Y)=dim Hom(X,Y)` are the one QPA-checkable slice).

**Tech Stack:** pure exact linear algebra over `Domain` (`modules/linalg_mod`
`matmul`/`mat_rank`/`kernel_columns`, `fields/linalg`), P41 `knit_ar_quiver`/
`ARQuiver`/`irreducible_maps`/`ar_quiver_block`/`_rad_end_basis`/`_vec`/`_unvec`,
P37 `hom_space`/`hom_dim`/`is_isomorphic`/`hom_basis`/`ModuleHom.then`, P30
`decompose`/`is_indecomposable`, `duality.tau`/`tau_minus`. No floats in `src/`
(AST-gated by `tests/test_no_floats.py`).

---

## The records (verbatim — the plan's charter)

**R21 — Liu degree theory + AR-component invariants; representation-directed
recognizer. [C-scout P10+P12; keep]** Object: left/right degrees of irreducible
maps, sectional paths, postprojective/preinjective/regular partition,
generalized-standard flags; directing modules and rep-directed (Γ_A acyclic) —
all finite sweeps on the knitted AR quiver (rep-finite). Refs: Liu JLMS 45 (1992)
32–54, JLMS 47 (1993) 405–416; Ringel LNM 1099; Bongartz CMH 57 (1982). Size M.
See R37 for the radical-filtration objects the degree theory computes with.

**R37 — The radical filtration of mod A and the infinite radical (Liu–Chaio
program). [added 2026-08-07 at Marco's request; anchors web-verified]** Object:
(a) the layer dimensions dim rad^n(X, Y) for indecomposables X, Y — on the
knitted rep-finite category rad(X,Y) is the non-isomorphism space and rad^n is
built by mesh-composites, with Liu's left/right degrees (R21) deciding exactly
when a composite of n irreducibles falls into rad^{n+1}; (b) the nilpotency index
of rad(mod A) for rep-finite A (the least n with rad^n = 0) — the Chaio-school
invariant relating it to AR-quiver structure ('Module Categories of Small Radical
Nilpotency', Alg. Rep. Theory 2023; 'On Sums of Compositions of Irreducible
Morphisms', ART 2018; 'Degrees of Irreducible Morphisms over Perfect Fields', ART
2019 / arXiv:1704.03933); (c) the certificate rad^∞(mod A) = 0 ⇔ A
representation-finite (Auslander) — an operational rep-finiteness gate dual to the
knitting bound; (d) for representation-INFINITE algebras, honest bounded-window
rad^n data plus the nilpotency taxonomy of rad^∞ as class oracles: (rad^∞)² = 0 ⇒
representation-finite (Coelho–Marcos–Merklen–Skowroński, 'Module categories with
infinite radical square zero are of finite type', Comm. Algebra 22(11) (1994));
the rep-infinite algebras with (rad^∞)³ = 0 include the tilted algebras of
Euclidean type (CMMS, 'Module categories with infinite radical cube zero', J.
Algebra (1996)); Kerner–Skowroński, 'On module categories with nilpotent infinite
radical', Compositio Math. (1991). Foundational bridge: Chaio–Liu, 'A note on the
radical of a module category' (Comm. Algebra 2013) — rep-finite representation
theory through radical nilpotency. Honest scope: fully certified on the rep-finite
domain (where quiverlab knits); on rep-infinite input the layers are
bounded-window computations and the (rad^∞)-nilpotency statements serve as class
pins, not per-instance deciders. Oracles: rad-layer tables on kA_n against the
mesh combinatorics; the nilpotency index on Nakayama zoo algebras; (rad^∞)²=0 ⇒
rep-finite as a discriminating battery; degree-vs-layer consistency with R21 (an
irreducible map of finite left degree d forces the predicted rad-layer drop).
Size M. Deps: AR knitting, Hom, R21 degrees.

---

## Reference re-verification (mandatory — done at authoring, 2026-08-07)

Web-verified against the primary sources; findings recorded so the implementer
does not re-derive and so the definitions are pinned against transcription error.

1. **Liu's left degree — the exact membership direction.** From *Degrees of
   Irreducible Morphisms over Perfect Fields* (arXiv:1704.03933, HTML mirror,
   §2), verbatim: the left degree `d_ℓ(f)` of `f ∈ rad^n \ rad^{n+1}` is "the
   least integer `m` such that there exists `(Z, g)` with `Z ∈ mod A`
   indecomposable and `g: Z → X` lying in `rad^m \ rad^{m+1}` and such that
   **`fg ∈ rad^{m+n+1}`**". The membership is **`∈` (IN), not `∉`.** For an
   **irreducible** `f`, `f ∈ rad \ rad²` so `n = 1`, and the condition is
   **`fg ∈ rad^{m+2}`**. **This corrects the P57 dispatch brief**, which wrote
   `fg ∉ rad^{n+2}` — that is a typo; the implementation uses `∈ rad^{m+2}`. The
   intuition confirms `∈`: for `g ∈ rad^m` and `f ∈ rad`, `fg ∈ rad^{m+1}`
   ALWAYS; the left degree detects the composite falling **deeper** than the
   generic `m+1`, i.e. into `rad^{m+2}` — a radical-layer DROP, which is exactly
   R37's "an irreducible map of finite left degree d forces the predicted
   rad-layer drop."
2. **Radical powers + infinite radical (same source, §2, verbatim):**
   `rad^0 = mod A`, `rad^{ℓ+1} = rad^ℓ·rad = rad·rad^ℓ`, and
   `rad^∞ = ⋂_{ℓ≥0} rad^ℓ`. Right degree `d_r(f)` is "defined dually" (mirror
   with `Z` on the target side: `g: Y → Z`, `gf ∈ rad^{m+n+1}`).
3. **Left-to-right convention.** quiverlab composes left-to-right (`a*b` = first
   `a` then `b`; ASS). Liu's `fg` is classical right-to-left (first `g`, then
   `f`). In quiverlab notation the left-degree condition reads: `g ∈ rad^m(Z,X)`,
   `f ∈ rad(X,Y)`, composite `g.then(f) ∈ rad^{m+2}(Z,Y)`; `ModuleHom.then`'s
   matrix is `f.matrix @ g.matrix` (a `Y.dim × Z.dim` matrix). The plan writes the
   left-to-right form throughout.
4. **Nilpotency index — definition and formulas.** From *On the nilpotency index
   of the radical of a module category* (Chaio–Guazzelli, arXiv:2003.04189, JPAA):
   the nilpotency index is "the minimal `m ≥ 1` such that `rad^m(mod A) = 0`".
   **Theorem 1.3 (general rep-finite `A ≅ kQ_A/I_A`):** the index is
   `max_{a ∈ Q_0} {r_a + 1}`, where `r_a` = the length of the nonzero path of
   irreducible morphisms between indecomposables from `P_a` to `I_a` going through
   `S_a`. **Theorem 1.5(a) (hereditary, underlying graph `A_n`):** the index is
   `n`. (Cross-checked by the hand derivation below: `index(kA_2) = 2`,
   `index(kA_3) = 3`.)
5. **Auslander's certificate** (`rad^∞ = 0 ⇔ rep-finite`). Not stated in
   1704.03933; it is the classical Auslander theorem, the operational content of
   R37(c). We treat it as the ground truth: on a **complete** knit (rep-finite)
   `rad^∞ = 0` and the finite nilpotency index is its witness; on an incomplete
   knit no `rad^∞` verdict is issued. **Foundational bridge (Chaio–Liu, Comm.
   Algebra 41(12) (2013) 4419–4424, web-confirmed):** characterizes rep-finiteness
   via the behaviour of projective covers / injective envelopes of simples under
   `rad^∞`, and shows the nilpotency of `rad(mod A)` in the rep-finite case is the
   maximal depth of the composites of those maps (independent of the maximal
   length of indecomposables) — the theorem behind Theorem 1.3.
6. **CMMS radical-square-zero (Comm. Algebra 22(11) (1994) 4511–4517,
   web-confirmed):** `(rad^∞)² = 0 ⇒ A of finite representation type`. **This is a
   CLASS oracle, NOT a per-instance decider here** — see the honest-scope section.
   CMMS cube-zero (J. Algebra 1996) and Kerner–Skowroński (Compositio 1991) are
   documentation-only class pins (rep-infinite `rad^∞ ≠ 0` — outside the certified
   domain).
7. **Liu 1993 title (BibTeX-verified):** "Semi-stable components of an
   Auslander–Reiten quiver", JLMS (2) 47(3) (1993) 405–416.
8. **Record citation-grouping correction (fold back to the research doc).** R37
   groups the nilpotency-index references as "(ART 2018/2019/2023,
   arXiv:1704.03933)". This conflates two distinct sources: **arXiv:1704.03933**
   (ART 2019) is *Degrees of Irreducible Morphisms over Perfect Fields* (the
   **degree** theory, R21), whereas the **nilpotency-index formula** (Thm 1.3
   `max_a{r_a+1}`, Thm 1.5(a) `index(kA_n)=n`) is **Chaio–Guazzelli,
   arXiv:2003.04189** (*On the nilpotency index of the radical of a module
   category*), which R37 does not cite. The plan uses 2003.04189 for the index
   and 1704.03933 only for the degree definition. **Implementer action (metaplan
   fold-back rule):** add a dated correction line to
   `docs/plans/2026-08-06-computability-expansion-deep-research.md` under R37:
   "2026-08-07 (P57): nilpotency-index formula source is arXiv:2003.04189
   (Chaio–Guazzelli), distinct from the degree paper 1704.03933; add 2003.04189 to
   R37's anchors."

`# PIN` items (leave the number/metadata to the implementer to freeze at build,
never guessed here): the full author lists / volume-page metadata of
**CMMS-1996 (cube zero)**, **Kerner–Skowroński 1991**, **Chaio–Guazzelli 2022
(JPAA, arXiv:2003.04189)**, **Chaio ART 2018 "On Sums of Compositions…"**, and
**Chaio ART 2023 "Module Categories of Small Radical Nilpotency"** — each added to
`references.bib` only after the worker BibTeX-verifies it at merge (the P49
precedent). The nilpotency-index *values* for the Nakayama zoo beyond the
hand-derived cases below are `# PIN`: computed by route (ii) and frozen at
implementation, cross-checked against Theorem 1.3.

---

## Mathematical foundation (derived — the plan's ground truth)

### The radical of the knitted category (route (ii), the arbiter)

Let `ind A = {X_1, …, X_r}` be the finite set of indecomposables produced by
`knit_ar_quiver` (rep-finite, `is_complete=True`). The radical of `mod A` is the
two-sided ideal `rad(X,Y) ⊆ Hom_A(X,Y)` of non-isomorphisms; on indecomposables:

- `rad(X_i, X_j) = Hom_A(X_i, X_j)` for `i ≠ j` (a map between **non-isomorphic**
  indecomposables is never an isomorphism, hence radical — this is exactly the
  P41 `_rad_basis` `i≠j` branch);
- `rad(X_i, X_i) = rad End_A(X_i)` (the non-units of the local ring `End(X_i)` —
  the P41 `_rad_end_basis` trace-form radical, char-scoped `char 0 or char > dim`).

Higher powers, in left-to-right composition (`p ∈ rad^n(X,Z)`, `q ∈ rad(Z,Y)`,
composite `p.then(q)`, matrix `q.matrix @ p.matrix`):

    rad^{n+1}(X, Y) = Σ_{Z ∈ ind A} { p.then(q) : p ∈ rad^n(X, Z), q ∈ rad(Z, Y) }.

**Why summing over INDECOMPOSABLE `Z` only loses nothing (Krull–Schmidt +
additivity, H2).** In general `rad^{n+1}(X,Y) = Σ_{W ∈ mod A} rad^n(X,W)∘rad(W,Y)`,
ranging over ALL modules `W`. By Krull–Schmidt `mod A = add(ind A)`: every `W`
splits `W ≅ ⊕_t Z_t` with `Z_t` indecomposable. `rad` and `Hom` are **additive
bifunctors**, so with the biproduct inclusions/projections `ι_t: Z_t → W`,
`π_t: W → Z_t` (`Σ_t ι_t π_t = 1_W`) any composite `p'∘q'` through `W` factors as
`Σ_t (p'∘ι_t)∘(π_t∘q')` with `p'∘ι_t ∈ rad^n(X,Z_t)` and `π_t∘q' ∈ rad(Z_t,Y)`
(radical is closed under composition with any morphism, and `ι_t`, `π_t` are the
structure maps of a biproduct of indecomposables). Hence every composite through
an arbitrary `W` is a sum of composites through the indecomposable summands `Z_t`,
and restricting the intermediate object to `ind A` reproduces `rad^{n+1}` exactly.
This is the same reduction that makes `rad^n(⊕X_i, ⊕Y_j) = ⊕_{i,j} rad^n(X_i,Y_j)`
(so the whole filtration is determined by its values on `ind A`).

**Bookkeeping (fixed once).** For each ordered pair `(i, j)` fix the basis
`H_{ij} = hom_basis(X_i, X_j)` of `Hom(X_i, X_j)` and identify a map with its
column-stacked `_vec`. Then each `rad^n(X_i, X_j)` is a **coordinate subspace** of
the `dim Hom(X_i, X_j)`-dimensional vec-space, stored as a spanning set of vecs;
its dimension is `mat_rank` of that set. `rad^1` is the coordinate subspace above.
`rad^{n+1}(X_i, X_j)` is the row-space (image) of the composites of `rad^n(X_i,X_k)`
against `rad^1(X_k, X_j)` over all `k`, computed by exact `matmul`. The composition
bilinear map is presentation-independent and exact — **route (ii) is the arbiter.**

**Nilpotency index.** `N := least n ≥ 1 with rad^n(X_i, X_j) = 0 for ALL i, j`
(equivalently the whole ideal `rad^n(mod A) = 0`). Because `ind A` is finite and
`rad^{n+1} ⊆ rad^n` strictly until it stabilizes at `0`, `N` is finite and
computable directly. `N` is Auslander's finite witness of `rad^∞ = 0`.

### Route (i) — mesh combinatorics (the crosscheck oracle, standardness-caveated)

On a **standard** component the mesh category `k(Γ_A)` is equivalent to `ind A`,
and `dim rad^n(X, Y)` equals the number of paths of length `n` from `X` to `Y` in
`Γ_A` modulo the mesh relations. For `kA_n` (and `D_4`) this gives closed-form
layer tables (derived below). **We do NOT ship a general functorial mesh engine
or an automatic standardness detector** (non-standard components exist —
Riedtmann's char-2 examples — and detecting them is a research-grade side
quest, DEFERRED as a named successor). Instead route (i) appears as: (a) the
`kA_n`/`D_4` closed-form layer/nilpotency values, pinned as **literature**
oracles; (b) a cheap **necessary-condition self-cert** that comes free from the
arrows — `dim rad^n(X,Y) ≠ 0 ⇒ there is a path of length ≥ n from X to Y in Γ_A`;
and (c) a `oracle_crossengine` assertion that the `kA_n`/`D_4` mesh closed form
**agrees with route (ii)** on the standard test algebras. Where a component were
non-standard the two could disagree — we state that boundary honestly on the
verification page; on the shipped test set they agree.

### Worked derivations (the hand-derived pins)

**`kA_2` (`1 → 2`, hereditary `A_2`).** `ind A = {S_1=[1], S_2=[2]=P_2,
P_1=[1,2]=I_2}`; AR quiver `S_2 → P_1 → S_1` (one AR sequence
`0 → S_2 → P_1 → S_1 → 0`, `τ S_1 = S_2`). Layer-1: `rad(S_2,P_1)`,
`rad(P_1,S_1)` each `1`; `rad(S_2,S_1) = Hom(S_2,S_1) = 0` (no path `2→1`). Layer-2:
the only candidate `rad²(S_2,S_1) = rad(S_2,P_1)∘rad(P_1,S_1)` — the composite of
the two mesh irreducibles is `0` (the AR-sequence composite). So `rad² ≡ 0`.
**Nilpotency index `N(kA_2) = 2`** (= `n`, Thm 1.5(a)). Thm 1.3 check: `r_1 = 1`
(`P_1 → S_1=I_1`), `r_2 = 1` (`P_2=S_2 → I_2=P_1`), `max{r_a+1} = 2`. ✓

**`kA_3` (`1 → 2 → 3`, hereditary `A_3`).** `ind A` = the 6 interval modules
`[1],[2],[3],[1,2],[2,3],[1,2,3]`. The `ZA_3` mesh has sectional composites of
length up to `2` that are nonzero; the longest nonzero radical composite gives
`rad² ≠ 0`, `rad³ = 0`. **Nilpotency index `N(kA_3) = 3`** (= `n`). Layer totals
(sum of `dim rad^n(X_i,X_j)` over ordered pairs) are pinned from route (ii) at
build (`# PIN`, cross-checked against the `ZA_3` path count).

**`kA_3/J²` (radical-square-zero truncated Nakayama, `1→2→3`, `ab=0`; NOT
self-injective, knittable).** `ind A = {S_1,S_2,S_3, P_1=[1,2], P_2=[2,3]}` (5
indecomposables — the P41 `test_ar_knit` Nakayama pin). Thm 1.3: `r_1 = 1`
(`P_1=[1,2] → S_1=I_1`), `r_2 = 2` (`P_2=[2,3] → S_2 → I_2=[1,2]`, through `S_2`),
`r_3 = 1` (`P_3=S_3 → I_3=[2,3]`). **Nilpotency index `= max{2,3,2} = 3`.**

**These are the `oracle_literature`/`oracle_crossengine` anchors.** `kA_4`
(`N = 4`) and further Nakayama-zoo indices are `# PIN` — computed by route (ii),
frozen at build, and cross-checked against Thm 1.3 by an in-plan `_thm13_index`
helper that walks the `P_a → S_a → I_a` path on the knit.

### Liu degrees on the knitted category (R21, the decidable form)

For an irreducible `f: X → Y` (`f ∈ rad(X,Y) \ rad²(X,Y)`, `X,Y ∈ ind A`):

    d_ℓ(f) = min over Z ∈ ind A, over d ≥ 1, of d such that
             { g ∈ rad^d(Z,X) : g.then(f) ∈ rad^{d+2}(Z,Y) }  ⊄  rad^{d+1}(Z,X);
             ∞ if no such (Z, d) with d < N exists.

i.e. there is a `g` of depth **exactly** `d` (`g ∈ rad^d \ rad^{d+1}`) whose
composite with `f` drops into `rad^{d+2}`. This is a finite decidable
linear-algebra sweep: for each `Z` and `d < N`, the set `{g ∈ rad^d(Z,X) :
g.then(f) ∈ rad^{d+2}(Z,Y)}` is the kernel of the (linear) map
`rad^d(Z,X) → rad^{d+1}(Z,Y)/rad^{d+2}(Z,Y)`, `g ↦ [g.then(f)]`; the witness
condition is that this kernel is **not contained** in `rad^{d+1}(Z,X)`.
Termination is guaranteed by `rad^N = 0`.

**Right degree — the formal dual (M2, stated once, no "flip and see").** `d_r(f)`
is the term-by-term dual with the target on the free side:

    d_r(f) = min over Z ∈ ind A, over d ≥ 1, of d such that
             { h ∈ rad^d(Y,Z) : f.then(h) ∈ rad^{d+2}(X,Z) }  ⊄  rad^{d+1}(Y,Z);
             ∞ if no such (Z, d) with d < N exists.

(the kernel of `rad^d(Y,Z) → rad^{d+1}(X,Z)/rad^{d+2}(X,Z)`, `h ↦ [f.then(h)]`,
not contained in `rad^{d+1}(Y,Z)`). This is exactly `d_ℓ` computed for `D f` over
`A^op` (`duality` is contravariant, swapping `X↔Y` and mono↔epi); the two agree by
construction and the `oracle_selfcert` `test_left_right_degree_op_symmetry` pins it
— there is no "flip once if it fails".

**Boundary witnesses make finite degrees the norm in rep-finite type (the review
insight).** At `d = N-1` the target layer `rad^{d+2} = rad^{N+1} = 0`, so the
witness condition `g.then(f) ∈ rad^{d+2}` becomes `g.then(f) = 0` — any nonzero
`g ∈ rad^{N-1}(Z,X)` (automatically `∉ rad^N = 0`) with `g.then(f) = 0` is a
legitimate witness. Because zero composites are ubiquitous on the knit (every
sink-directed pair `Hom(Z,Y) = 0` gives one), **most irreducible maps in a
rep-finite algebra have FINITE degree** — the naïve intuition "hereditary ⇒
infinite degree" is WRONG (see `kA₃` below, where three of six irreducibles have
finite left degree).

**Degree-vs-layer consistency (the R37 self-cert oracle).** When `d_ℓ(f) = d < ∞`
the routine returns the witness `(Z, g)`, and the test asserts the exact
membership TRIPLE **all three externally** (M3): `g ∈ rad^d(Z,X)`,
`g ∉ rad^{d+1}(Z,X)`, and `g.then(f) ∈ rad^{d+2}(Z,Y)` — the certificate that a
finite left degree forces the predicted layer drop.

### Worked degree derivations (kA₂, kA₃ — the in-plan pins, live-cross-checked)

The rad structure below was **verified live** at authoring against the engine
(`knit_ar_quiver` + `Algebra.hom`): all `kA₂`/`kA₃` AR arrows and all nonzero
radical Hom pairs match this hand derivation exactly (route (ii) will reproduce
them). Degrees are stated **by module name**, robust to the knit's BFS vertex
indexing. Composites use the verified mesh relations. `N(kA₂)=2`, `N(kA₃)=3`.

**`kA₂` (`1→2`), `ind A = {S₂=P₂, P₁=[1,2]=I₂, S₁=I₁}`, one AR sequence
`0→S₂→P₁→S₁→0`.** Two irreducibles:

| irreducible map        | type | `d_ℓ` | `d_r` |
|------------------------|------|-------|-------|
| `f₁: S₂ → P₁`          | mono | `∞`   | `1`   |
| `f₂: P₁ → S₁`          | epi  | `1`   | `∞`   |

*Derivation.* `rad = Hom` off-diagonal; the only nonzero radical Homs are
`f₁` (`Hom(S₂,P₁)=k`), `f₂` (`Hom(P₁,S₁)=k`); `Hom(S₂,S₁)=0`, so
`f₁.then(f₂)=0` (`rad²=0`). **`d_ℓ(f₂)=1`:** at `d=1`, `Z=S₂`, `g=f₁ ∈ rad(S₂,P₁)\rad²`,
`g.then(f₂)=f₁.then(f₂)=0 ∈ rad³=0`. ✓ **`d_ℓ(f₁)=∞`:** `S₂=P₂` is a projective
source, `rad(Z,S₂)=Hom(Z,S₂)=0` for all `Z` (no module maps into `S₂`), so no
witness exists. **`d_r`** is the mirror (`S₁=I₁` is an injective sink,
`rad(S₁,Z)=0` ⇒ `d_r(f₂)=∞`; `Z=S₁`, `h=f₂`, `f₁.then(f₂)=0` ⇒ `d_r(f₁)=1`). This
reproduces Liu's classical fact: the AR epi/mono of a sequence with **indecomposable
middle** has degree `1`.

**`kA₃` (`1→2→3`), `ind A =` the six interval modules; `N=3`.** rad¹ = the six
arrows; rad² = the three length-2 composites `S₃→P₁`, `P₂→I₂`, `P₁→S₁` (each dim
1); rad³ = 0. The three meshes give the relations `S₃→P₂→S₂ = 0`,
`P₂→P₁→I₂ + P₂→S₂→I₂ = 0`, `S₂→I₂→S₁ = 0`. The six irreducibles and their degrees:

| irreducible map          | type | `d_ℓ` | `d_r` |
|--------------------------|------|-------|-------|
| `S₃ → P₂`  (`[3]→[2,3]`)     | mono | `∞`   | `1`   |
| `P₂ → P₁`  (`[2,3]→[1,2,3]`) | mono | `∞`   | `2`   |
| `P₂ → S₂`  (`[2,3]→[2]`)     | epi  | `1`   | `∞`   |
| `P₁ → I₂`  (`[1,2,3]→[1,2]`) | epi  | `2`   | `∞`   |
| `S₂ → I₂`  (`[2]→[1,2]`)     | mono | `∞`   | `1`   |
| `I₂ → S₁`  (`[1,2]→[1]`)     | epi  | `1`   | `∞`   |

*Key derivations.* **`d_ℓ(P₂→S₂)=1`:** `Z=S₃`, `g=(S₃→P₂) ∈ rad\rad²`,
`g.then(f)=(S₃→P₂→S₂)`, and `Hom(S₃,S₂)=0` ⇒ `=0 ∈ rad³`. ✓ **`d_ℓ(I₂→S₁)=1`:**
`Z=S₂`, `g=(S₂→I₂)`, `g.then(f)=(S₂→I₂→S₁)=0` by the mesh relation at `S₁`. ✓
**`d_ℓ(P₁→I₂)=2` (NOT 1):** at `d=1` the only usable `Z` is `P₂` with `g=(P₂→P₁)`,
but `g.then(f)=(P₂→P₁→I₂) = −(P₂→S₂→I₂) ≠ 0` (mesh, nonzero); and `Z=S₃` fails
because `rad(S₃,P₁)=rad²(S₃,P₁)` (the map `S₃→P₁` is itself a rad² composite), so
`rad¹\rad²` is empty there — no `d=1` witness. At `d=2`, `Z=S₃`,
`g=(S₃→P₁) ∈ rad²\rad³`, `g.then(f)=(S₃→P₁→I₂)`, `Hom(S₃,I₂)=Hom([3],[1,2])=0`
⇒ `=0 ∈ rad⁴`. ✓ so `d_ℓ=2`. **`d_ℓ(monos)=∞`:** each mono `X→Y` has `X ∈
{S₃, P₂, S₂}`; `S₃=P₃` is a projective source (`rad(·,S₃)=0`), and for `P₂`, `S₂`
the only `d=1`/`d=2` candidates give nonzero composites (`P₂→P₁` gives the nonzero
`P₂→P₁ = S₃`-composite chain, etc.). The `d_r` column is the exact
opposite-algebra mirror (mono↔epi, `d_ℓ↔d_r`), verified directly for
`d_r(S₃→P₂)=1` and `d_r(P₂→P₁)=2`.

**Structural observation (a strong self-cert on `kA_n`).** For the hereditary
directed `kA_n`: **every irreducible mono has `d_ℓ=∞`, `d_r` finite; every
irreducible epi has `d_ℓ` finite, `d_r=∞`.** The plan pins this per-arrow
(`oracle_literature`) — it is the corrected replacement for the deleted, and
mathematically false, "all irreducibles have infinite left degree" stub.

### R21 component invariants (finite sweeps on the knit)

**What is genuinely non-degenerate in scope (W2, honest framing).** The
knittable rep-finite domain contains both **representation-directed** algebras
(`kA_n`, `D_4`: `Γ_A` acyclic) AND **non-directed** ones — the live-verified
witness `NakayamaAlgebra(kupisch=[3,2,2])` (cyclic `Z₃`, non-uniform Kupisch) is
**non-self-injective, knit-complete (7 indecomposables), and has an oriented cycle
in `Γ_A`**. So the recognizer's negative branch is **reachable and tested**, not
vacuous; the R21 objects that are non-degenerate in scope are the **degrees**
(finite and infinite both occur), the **sectional paths**, and the **directing /
representation-directed** recognizer (both verdicts realized). The τ-partition is
the clean classical trichotomy only on directed components (below).

- **τ-partition.** For `X ∈ ind A`, walk the `τ`-orbit within the knit
  (`duality.tau`/`tau_minus`, finite): `X` is **postprojective** if the orbit
  contains a projective, **preinjective** if it contains an injective, **regular**
  otherwise. On a **representation-directed** component (`kA_n`, `D_4`) this is the
  classical trichotomy: every indecomposable is postprojective = preinjective and
  there are **no regular** modules — the clean partition oracle. On a **non-directed**
  component (the `[3,2,2]` cyclic Nakayama) the τ-orbit classification still runs
  but postprojective/preinjective can overlap and the "regular" bucket flags the
  cyclic part; the block reports it **honestly labeled** ("τ-orbit classification;
  non-directed component — buckets may overlap"), never as the clean trichotomy.
- **Directing / representation-directed.** `X` is **directing** if it lies on no
  oriented cycle `X = X_0 → X_1 → … → X_r = X` (`r ≥ 1`) of the `Γ_A` arrow
  digraph (irreducible maps); `A` is **representation-directed** iff it is
  rep-finite (knit complete) and `Γ_A` is acyclic — every indecomposable
  directing. Computed by a Tarjan-SCC / digraph cycle test on `ARQuiver.arrows`
  (`kA_n`/`D_4` ⇒ `True`; `[3,2,2]` ⇒ `False`, the live-verified negative branch),
  refined by a route-(ii) check that no cycle carries a **nonzero** composite
  (Ringel: `X` directing ⇔ no cycle of nonzero non-isos through `X`). **Note the
  two flags are independent:** `[3,2,2]` is `representation_directed = False` yet
  `generalized_standard = True` (rep-finite ⇒ `rad^∞ = 0`) — the plan pins exactly
  this separation.
- **Sectional paths.** A path `X_0 →α_1 X_1 → … → X_r` in `Γ_A` is **sectional**
  if for no `i` is `X_{i-1} = τ X_{i+1}` (no mesh reversal). Sectional composites
  are nonzero (Bautista–Smalø): a sectional path of length `r` yields a nonzero
  element of `rad^r \ rad^{r+1}` — the "no-drop" certificate that pairs with the
  degree theory. We report a `is_sectional(path)` checker and the **maximal
  sectional-path length** per component.
- **Generalized standard.** A component `Γ` is generalized standard (Skowroński)
  if `rad^∞(X,Y) = 0` for all `X,Y ∈ Γ`. In our certified (rep-finite) scope
  `rad^∞ = 0` (Auslander), so **every** knitted component is generalized standard,
  witnessed by the finite nilpotency index — reported `True` with that witness and
  the honest note that the flag is only *informative* for rep-infinite algebras
  (outside the certified domain).

### Honest scope (metaplan §6 — stated on the verification page, never overclaimed)

- **Rep-finite only, certified.** The nilpotency index, `rad^∞ = 0`, the degrees,
  and every R21 invariant are certified **iff the knit closes** (`is_complete`).
  This inherits P41's contract exactly: `status ∈ {complete, budget, error,
  unsupported}`.
- **Self-injective refused.** `knit_ar_quiver` returns `status="unsupported"` for
  self-injective algebras (`modules/ar.py:621`, projective-seeded BFS would
  undercount). Therefore **cyclic Nakayama `kZ_n/J^ℓ` (self-injective) is out of
  scope** — `radical_filtration` refuses it loudly, and the Nakayama nilpotency
  zoo is the **non-self-injective** (truncated linear) `kA_n/J^ℓ`. The Chaio-school
  self-injective closed forms are cited as documentation, not computed.
- **Rep-infinite = window only, no verdict.** On a budget-exhausted knit,
  `radical_filtration` returns the layers of the **discovered subcategory**
  `add(𝒟)` — exact within the window but a **lower bound / partial view**, NOT
  `rad(mod A)` — with `is_complete=False`, a loud note, and **no** nilpotency-index
  or `rad^∞` claim. This is R37(d)'s "honest bounded-window rad^n data."
- **The `(rad^∞)²=0 ⇒ rep-finite` battery is a CLASS statement.** We CANNOT
  compute `rad^∞` on rep-infinite input (no complete knit), so we NEVER construct
  a rep-infinite instance and "verify CMMS." What we test: (1) every rep-finite
  zoo algebra has a finite nilpotency index (`rad^∞ = 0`, a fortiori
  `(rad^∞)² = 0`, and rep-finite — consistent with CMMS); (2) the **discrimination**
  — our gate returns a finite index (rep-finiteness certificate) on rep-finite
  input and an honest `status="budget"`/`"unsupported"` (no false finite index) on
  rep-infinite input (`m`-Kronecker, wild). CMMS/Kerner–Skowroński are the external
  class oracles that JUSTIFY the finite-index gate as a legitimate rep-finiteness
  certificate; they are documentation, not per-instance deciders.
- **Char scope.** `rad End(X_i)` (the `i=i` diagonal) uses the P41/P30 trace-form
  radical, rigorous over **char 0 or char > dim `X_i`** (`ar.py:354`,
  `decompose.py:321`). Batteries run over **QQ** (with a `GF(32003)` parity
  cross-check where cheap); `char ≤ dim` inherits the loud `QuiverlabError` refusal.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P41 is MERGED to `dev` and is the hard prerequisite** (verified at authoring:
  `modules/ar.py` present; `knit_ar_quiver`/`ARQuiver`/`irreducible_maps`/
  `ar_quiver_block`/`_rad_end_basis`/`_vec`/`_unvec` at the signatures in the P41
  plan). P37 (`hom_space`/`hom_dim`/`is_isomorphic`/`hom_basis`/`ModuleHom`) and
  P30 (`decompose`/`is_indecomposable`) are merged. Branch `plan-57-degrees-radical`
  off `dev`; **design against `dev`** (P51/P52/P53 run in separate worktrees — do
  not touch them).
- **`tests/modules/` auto-assigns to the deep bucket** (`tests/conftest.py`
  `_DEEP_DIRS`); `tests/qpa/` to the qpa bucket. Run new tests by path during
  development; finish each task with a `-m deep` (or `-m qpa`) spot-run of the
  touched files.
- **Consume P41 at its verified signatures:** `knit_ar_quiver(A, budget_modules=256,
  budget_dim=4096) -> ARQuiver` with `.vertices` (list of `{"name","dimvec",
  "module"}`, index = list position), `.arrows` (`dict[(i,j)->mult]`, populated
  only when complete), `.tau_orbits`, `.is_complete`, `.status ∈ {complete,budget,
  error,unsupported}`, `.note`; `irreducible_maps(M,N,within) -> int`;
  `_rad_end_basis(M) -> (H, rad_coords)` (char-scoped); `_vec`/`_unvec`;
  `ar_quiver_block(A, budget=512)` (the algebra-block serializer template).
- **Char scope is load-bearing** (above). Diagonal `rad End(X_i)` rigorous over
  char 0 / char > dim; batteries over QQ.
- **Honest semi-decision contract (metaplan §6):** `RadicalFiltration`/
  `ARInvariants` mirror P41's `ARQuiver.is_complete`/`.status` exactly — certified
  iff rep-finite; loud/window otherwise; never a silently truncated verdict.
- **No floats in `src/`.** Layer dims / indices / degrees are `int`; layer
  sequences `list[int]`; per-pair tables `dict`. `∞` is represented as `None` with
  a boolean `finite` flag (never `float('inf')`).
- **Composition is left-to-right** (`f.then(g)`). `_assert_comparable` guards every
  cross-module call. All refusals are `QuiverlabError`.
- Plan-32 markers: the filtration/degree certificates (`rad^{n+1}⊆rad^n`,
  `rad^N=0 & rad^{N-1}≠0`, composition closure, the degree-witness triple, poset/
  partition axioms) = `oracle_selfcert`; two-independent-route agreement (mesh
  route-(i) closed form ≡ route-(ii); Thm 1.3 index ≡ route-(ii) index;
  `dim rad(X,Y) ≡ dim Hom(X,Y)` for `X≇Y`) = `oracle_crossengine`; `kA_n`/`D_4`
  mesh layer tables + `N(kA_n)=n` + hand-derived Nakayama indices + "rep-directed,
  no regular, all directing on `kA_n`" = `oracle_literature`; QPA lives in
  `tests/qpa/` (bucket = the class, never double-marked). **Oracle-marked tests
  live in `tests/modules/`** (always collectible) — NEVER in `tests/webapp`/`gui`
  (the release gate forbids markers there).
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted
  class table green, `tests/release/test_oracle_classes.py` green) and adds its
  citations to `citations/references.bib` + `registry.py` (`bibtex()` hard-fails
  if the two drift). Conventional commits; green at every commit.

---

### Task 1: `modules/radical.py` — the radical-filtration engine (R37 core)

**Files:**
- Create: `src/quiverlab/modules/radical.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.radical_filtration(budget_modules=256, budget_dim=4096)` thin delegate, beside `ar_quiver`, `core/algebra.py:421`)
- Test: `tests/modules/test_radical_filtration.py`

**Interfaces:**
- Consumes: `knit_ar_quiver`/`ARQuiver` (P41), `_rad_end_basis`/`_vec`/`_unvec`
  (P41, `modules/ar.py`), `hom_basis` (P37 `morphism.py`), `is_isomorphic` (P37
  `hom.py`), `modules/linalg_mod` (`matmul`, `mat_rank`, `cols_to_matrix`).
- Produces:
  ```python
  @dataclass
  class RadicalFiltration:
      indecs        # list[Module], the ind A universe (== ARQuiver vertex modules)
      names         # list[str|None], identify_standard names
      layer         # dict[(i, j) -> list[int]]: layer[(i,j)][n-1] = dim rad^n(X_i, X_j)
      nilpotency_index  # int N (least n with rad^n ≡ 0) -- None if not is_complete
      rad_infinity_zero # True iff is_complete (rep-finite => Auslander); None otherwise
      is_complete   # bool: True iff the knit closed (rep-finite)
      status        # "complete" | "budget" | "unsupported" | "error"
      note          # str
      # private: _radspaces[(i,j)] = list of {n: spanning vecs of rad^n(X_i,X_j)}
      def pair_layers(self, X, Y) -> list[int]      # dim rad^n(X,Y) for the matched pair
      def layer_dim(self, i, j, n) -> int
      def is_generalized_standard(self) -> bool      # True in rep-finite scope (witness = N)
  def radical_filtration(A, *, budget_modules=256, budget_dim=4096) -> RadicalFiltration
  ```
  `Algebra.radical_filtration(...)` delegates (lazy-import).

**Algorithm (route (ii), exactly as the foundation section):**
1. `ar = knit_ar_quiver(A, budget_modules, budget_dim)`. If `not ar.is_complete`,
   compute the **window** layers on `ar.vertices` (the discovered submodule set) —
   exact within `add(𝒟)` — and return with `is_complete=False`,
   `nilpotency_index=None`, `rad_infinity_zero=None`, `status=ar.status`, and the
   loud window note. **No verdict.** (Self-injective `status="unsupported"` short-
   circuits before any layer work, pointing at the honest scope.)
2. Else `indecs = [v["module"] for v in ar.vertices]`; fix `H_{ij} = hom_basis`
   once; build `rad^1` coordinate subspaces (full Hom if `i≠j`; `rad End` via
   `_rad_end_basis` + `_combine`/`_vec` if `i=i` — inherit the P41 char guard).
3. Iterate `n = 1, 2, …`: `rad^{n+1}(X_i,X_j)` = the vecs of `q.matrix @ p.matrix`
   for `p ∈ rad^n(X_i,X_k)`, `q ∈ rad^1(X_k,X_j)`, all `k`; `dim` via `mat_rank`.
   Stop when all pairs give `0` → that `n` is `N`.
4. `layer[(i,j)] = [dim rad^n(X_i,X_j) for n in 1..N-1]`; `nilpotency_index = N`;
   `rad_infinity_zero = True`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_radical_filtration.py
"""The radical filtration of mod A (Plan 57 / R37). Self-cert: rad^{n+1} subset
rad^n, rad^N = 0 while rad^{N-1} != 0, composition closure. Literature: kA_n layer
tables + nilpotency index N(kA_n)=n; kA_3/J^2 index 3 (Thm 1.3). Cross-engine: the
mesh closed form == route (ii); Thm-1.3 index == route (ii). Honest scope:
self-injective refused (status 'unsupported'); rep-infinite -> window, no verdict."""
import pytest

from quiverlab import Quiver, linear_path_algebra, NakayamaAlgebra
from quiverlab.fields import QQ, GF
from quiverlab.modules.radical import radical_filtration

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_nilpotency_index_of_linear_a_n_is_n(n):
    rf = radical_filtration(linear_path_algebra(n, field=QQ))
    assert rf.is_complete and rf.status == "complete"
    assert rf.nilpotency_index == n            # Thm 1.5(a); hand-derived n=2,3
    assert rf.rad_infinity_zero is True         # Auslander witness = N


@lit
def test_ka3_layer_table_matches_the_mesh():
    A = linear_path_algebra(3, field=QQ)
    rf = radical_filtration(A)
    # S_2 -> [1,2] -> [1] is a length-2 mesh path; pin the specific pair layers
    # SHARPEN in Step 3 from the ZA_3 mesh (route (i)) and assert == route (ii).
    ...


@lit
def test_ka3_radsq_truncated_nakayama_index_is_3():
    # kA_3/J^2 (1->2->3, ab=0): NOT self-injective, 5 indecomposables, index 3.
    A = NakayamaAlgebra(n=3, l=2, cyclic=False, field=QQ)
    rf = radical_filtration(A)
    assert rf.is_complete and rf.nilpotency_index == 3   # Thm 1.3: max{2,3,2}


@selfcert
def test_layers_are_a_descending_filtration_and_terminate():
    rf = radical_filtration(linear_path_algebra(3, field=QQ))
    r = len(rf.indecs)
    N = rf.nilpotency_index
    for i in range(r):
        for j in range(r):
            seq = rf.layer[(i, j)]
            for n in range(1, len(seq)):
                assert seq[n] <= seq[n - 1]            # rad^{n+1} subset rad^n
    # rad^{N-1} != 0 somewhere, rad^N == 0 everywhere
    assert any(rf.layer_dim(i, j, N - 1) > 0 for i in range(r) for j in range(r))
    assert all(rf.layer_dim(i, j, N) == 0 for i in range(r) for j in range(r))


@selfcert    # NOT cross-engine: the Thm-1.3 walk and route (ii) both read the SAME
def test_thm13_index_equals_route_ii():
    # knit (one engine, two readings). This certifies the index computation is
    # internally consistent with the Chaio-Guazzelli formula on our own knit.
    from quiverlab.modules.radical import _thm13_index          # in-plan helper
    for A in (linear_path_algebra(3, field=QQ),
              NakayamaAlgebra(n=3, l=2, cyclic=False, field=QQ)):
        rf = radical_filtration(A)
        assert _thm13_index(A) == rf.nilpotency_index


@xeng
def test_mesh_layer_dims_match_route_ii_on_ka3():
    # GENUINELY INDEPENDENT route (i): count paths of length n in Gamma_A modulo
    # the mesh relations (a pure combinatorial walk on (arrows, tau) -- no Hom, no
    # matmul), and assert it reproduces route (ii)'s dim rad^n(X,Y) on kA_3
    # (standard directed => mesh count is exact). This is the mesh-vs-linear-algebra
    # agreement (a strong standardness self-cert on a standard component).
    from quiverlab.modules.radical import radical_filtration, _mesh_layer_dim
    from quiverlab.modules.ar import knit_ar_quiver
    A = linear_path_algebra(3, field=QQ)
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    r, N = len(rf.indecs), rf.nilpotency_index
    for i in range(r):
        for j in range(r):
            for n in range(1, N + 1):
                assert _mesh_layer_dim(ar, i, j, n) == rf.layer_dim(i, j, n)


@selfcert
def test_self_injective_refused_as_unsupported():
    # cyclic Nakayama kZ_3/J^2 is self-injective -> knit unsupported -> no verdict.
    A = NakayamaAlgebra(n=3, l=2, cyclic=True, field=QQ)
    rf = radical_filtration(A)
    assert rf.is_complete is False and rf.status == "unsupported"
    assert rf.nilpotency_index is None and rf.rad_infinity_zero is None


@selfcert
def test_representation_infinite_is_window_only_no_verdict():
    # 3-Kronecker (wild): knit budget-exhausts -> window layers, no index/rad^inf.
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2), "c": (1, 2)}).algebra(
        relations=[], field=QQ)
    rf = radical_filtration(A, budget_modules=40)
    assert rf.is_complete is False and rf.status in ("budget", "error")
    assert rf.nilpotency_index is None and rf.rad_infinity_zero is None
    assert "window" in rf.note.lower()          # loud: NOT rad(mod A)


@selfcert
def test_layer1_equals_hom_off_diagonal_field_parity():
    for field in (QQ, GF(32003)):
        rf = radical_filtration(linear_path_algebra(3, field=field))
        # rad(X_i, X_j) = Hom(X_i, X_j) for i != j (all maps radical)
        from quiverlab.modules.hom import hom_dim
        for i, Xi in enumerate(rf.indecs):
            for j, Xj in enumerate(rf.indecs):
                if i != j:
                    assert rf.layer_dim(i, j, 1) == hom_dim(Xi, Xj)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.modules.radical`

- [ ] **Step 3: Implement** per the algorithm. Key points / adjust-to-reality:
  - **The `i=i` diagonal is the only char-scoped piece.** Wrap the whole build in
    the P41 char guard: `_rad_end_basis` raises loudly over `char ≤ dim`; let it
    propagate (batteries run over QQ). Off-diagonal `rad = Hom` needs no char.
  - **`mat_rank`-based dims, not spanning-set sizes** (the composites over-span).
    Store each `rad^n(X_i,X_j)` as a reduced column basis (RREF via `mat_rank`'s
    pivots) so later composition is cheap and membership tests (`solve_columns`
    for the degree sweep in Task 2) are exact.
  - **Membership `v ∈ rad^{n+2}` is `solve_columns(basis_of_rad^{n+2}, v)`** — the
    same primitive Task 2's degree sweep reuses; expose an internal
    `_in_layer(rf, i, j, n, vec) -> bool`.
  - **`_thm13_index(A)`**: for each vertex `a`, find the projective `P_a`, simple
    `S_a`, injective `I_a` in the knit; walk the longest nonzero `Γ_A` path
    `P_a → … → S_a → … → I_a` (nonzero-composite verified by route (ii)); return
    `max_a{r_a + 1}`. **This reads the SAME knit/route-(ii) layers, so its
    agreement with `nilpotency_index` is `oracle_selfcert` (two readings of one
    engine), NOT cross-engine** — the internal-consistency check against the
    Chaio–Guazzelli formula.
  - **`_mesh_layer_dim(ar, i, j, n)` — the genuinely INDEPENDENT route (i)**: a
    pure combinatorial walk on `(ar.arrows, tau)` — count the paths of length `n`
    from `X_i` to `X_j` in `Γ_A` and quotient by the mesh relations (each mesh
    `Σ_{α: W→Z} (τα)·α = 0` identifies the two paths through a mesh; on the acyclic
    directed `kA_n` the surviving count = paths modulo those identifications). It
    touches **no `Hom`, no `matmul`** — hence its agreement with `layer_dim`
    (`test_mesh_layer_dims_match_route_ii_on_ka3`) is a real `oracle_crossengine`
    and the mesh-vs-linear-algebra standardness self-cert. Scoped to standard
    directed components (`kA_n`/`D_4`); on a non-standard component the two could
    diverge (the documented honest boundary — a general standardness detector is
    the named deferral).
  - **Window branch:** compute the layers on `ar.vertices` (the discovered set) by
    the SAME route-(ii) code but WITHOUT claiming an index; set the note to
    `"window-restricted to the {k} discovered indecomposables (add of the frontier);
    this is a lower bound on rad(mod A), NOT the category radical -- no nilpotency
    index or rad^inf verdict is issued (knit status: {status})"`.
  - **`is_generalized_standard`**: return `is_complete` (Auslander: `rad^∞=0` in
    rep-finite scope, witness `= nilpotency_index`); docstring states it is
    informative only for rep-infinite input (out of certified scope).

- [ ] **Step 4: Run tests** (sharpen the `kA_3` layer pin from the `ZA_3` mesh) — PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/radical.py src/quiverlab/core/algebra.py tests/modules/test_radical_filtration.py
git commit -m "feat(radical): radical filtration of mod A -- exact rad^n(X,Y) layers, nilpotency index, rad^inf=0 certificate (rep-finite); honest window on budget/self-injective"
```

---

### Task 2: `modules/ar_invariants.py` — Liu left/right degrees + sectional paths (R21)

**Files:**
- Create: `src/quiverlab/modules/ar_invariants.py`
- Test: `tests/modules/test_liu_degrees.py`

**Interfaces:**
- Consumes: `RadicalFiltration` (Task 1) + its `_in_layer`/`layer_dim`/basis
  accessors, `ARQuiver.arrows` (the irreducible-map arrows), `_vec`/`_unvec`
  (P41), `modules/linalg_mod` (`matmul`, `mat_rank`, `solve_columns`,
  `kernel_columns`), `duality.tau`/`tau_minus`.
- Produces:
  ```python
  def left_degree(rf, i, j, f_vec=None) -> tuple  # (d|None, finite, witness)
      # d_l of the irreducible-map class X_i -> X_j (f_vec: a rep of rad/rad^2;
      # default = the first basis vector). Returns (d, True, (z, g_vec, comp_vec))
      # -- comp_vec = coords of g.then(f) in Hom(X_z, X_j), so the caller can assert
      # ALL THREE memberships of the witness triple externally (M3) -- or
      # (None, False, None) for infinite degree (sweep exhausted, all d < N).
  def right_degree(rf, i, j, f_vec=None) -> tuple  # dual (h on the target side)
  def is_sectional(ar, path) -> bool               # no mesh reversal X_{k-1}=tau X_{k+1}
  def max_sectional_length(ar) -> int
  def degree_table(rf, ar) -> dict                 # {(i,j): {"d_l":..,"d_r":..}} over arrows
  ```

**Left-degree sweep (the foundation formula, left-to-right):** for `f` a rep of
`rad(X_i,X_j)/rad²`, for `d = 1 … N-1`, for each `Z = X_z`:
- `K = { g ∈ rad^d(Z,X_i) : g.then(f) ∈ rad^{d+2}(Z,X_j) }` = kernel of the linear
  map `rad^d(Z,X_i) → rad^{d+1}(Z,X_j)/rad^{d+2}(Z,X_j)`, `g ↦ [f.matrix @ g]`
  (build the map's matrix in the fixed bases; `K = kernel_columns`);
- witness iff `K ⊄ rad^{d+1}(Z,X_i)`, i.e. some `g ∈ K` has
  `_in_layer(rf, z, i, d+1, g) is False`;
- the least such `d` (over all `Z`) is `d_ℓ(f)`; return the witness `g`. Exhaust to
  `N-1` ⇒ infinite (`None, False, None`).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_liu_degrees.py
"""Liu left/right degrees on the knitted category (Plan 57 / R21+R37). Self-cert:
the degree-vs-layer witness triple (g in rad^d, g not in rad^{d+1}, g.then(f) in
rad^{d+2}) -- a finite left degree forces the layer drop. Literature: hand-derived
small degrees + sectional composites nonzero."""
import pytest

from quiverlab import linear_path_algebra, NakayamaAlgebra
from quiverlab.fields import QQ
from quiverlab.modules.ar import knit_ar_quiver
from quiverlab.modules.radical import radical_filtration
from quiverlab.modules.ar_invariants import (left_degree, right_degree,
                                             degree_table, max_sectional_length)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature

INF = None   # infinite degree is represented as None (no float inf; see constraints)


def _named_arrow(ar, src, tgt):
    """Look an arrow up by MODULE NAME (indexing-robust -- the knit orders vertices
    by BFS, so integer indices are not stable across refactors)."""
    idx = {v["name"]: i for i, v in enumerate(ar.vertices)}
    return idx[src], idx[tgt]


# The full HAND-DERIVED degree tables (verified live against the engine's rad
# structure at authoring), keyed (src_name, tgt_name) -> (d_l, d_r):
KA2 = {("S_2", "P_1"): (INF, 1),      # mono
       ("P_1", "S_1"): (1, INF)}      # epi
KA3 = {("S_3", "P_2"): (INF, 1),      # mono
       ("P_2", "P_1"): (INF, 2),      # mono
       ("P_2", "S_2"): (1, INF),      # epi
       ("P_1", "I_2"): (2, INF),      # epi (degree 2, NOT 1 -- decomposable mesh middle)
       ("S_2", "I_2"): (INF, 1),      # mono
       ("I_2", "S_1"): (1, INF)}      # epi


@selfcert
def test_finite_left_degree_witness_triple():
    # M3: assert ALL THREE memberships of the witness triple EXTERNALLY.
    A = NakayamaAlgebra(kupisch=[3, 2, 2], field=QQ)   # has finite-degree arrows
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    seen_finite = 0
    for (i, j) in ar.arrows:                            # every irreducible map
        d, finite, wit = left_degree(rf, i, j)
        if not finite:
            continue
        seen_finite += 1
        z, g, comp = wit                               # comp = coords of g.then(f)
        assert rf._in_layer(z, i, d, g) is True         # g in rad^d(Z, X_i)
        assert rf._in_layer(z, i, d + 1, g) is False    # g not in rad^{d+1}
        assert rf._in_layer(z, j, d + 2, comp) is True  # g.then(f) in rad^{d+2}(Z, X_j)
    assert seen_finite >= 1                             # not vacuous: a finite degree occurs


@lit
@pytest.mark.parametrize("n, table", [(2, "KA2"), (3, "KA3")])
def test_liu_degree_tables(n, table):
    # the corrected replacement for the (false) "all irreducibles have infinite
    # left degree" stub: the FULL hand-derived left+right degree tables, by name.
    A = linear_path_algebra(n, field=QQ)
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    expected = {2: KA2, 3: KA3}[n]
    assert len(expected) == len(ar.arrows)             # every irreducible covered
    for (src, tgt), (dl, dr) in expected.items():
        i, j = _named_arrow(ar, src, tgt)
        assert left_degree(rf, i, j)[0] == dl,  (src, tgt, "d_l")
        assert right_degree(rf, i, j)[0] == dr, (src, tgt, "d_r")


@lit
def test_mono_epi_degree_dichotomy_on_ka_n():
    # structural literature pin: on hereditary directed kA_n every irreducible MONO
    # has d_l = INF & d_r finite; every EPI has d_l finite & d_r = INF.
    for n in (2, 3, 4):
        A = linear_path_algebra(n, field=QQ)
        rf, ar = radical_filtration(A), knit_ar_quiver(A)
        for (i, j) in ar.arrows:
            dl = left_degree(rf, i, j)[0]
            dr = right_degree(rf, i, j)[0]
            # mono <=> d_l = INF (d_r finite); epi <=> d_r = INF (d_l finite):
            # exactly one side is infinite for every irreducible on kA_n.
            assert (dl is INF) ^ (dr is INF)


@lit
def test_sectional_composite_is_nonzero():
    # a sectional path of length r yields a nonzero element of rad^r \ rad^{r+1}.
    A = linear_path_algebra(3, field=QQ)
    ar = knit_ar_quiver(A)
    assert max_sectional_length(ar) >= 2          # ZA_3 has a length-2 sectional path


@selfcert
def test_left_right_degree_op_symmetry():
    # d_r(f) over A == d_l(D f) over A^op (the formal dual; NOT a "flip and see").
    A = linear_path_algebra(3, field=QQ)
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    Aop = A.opposite()
    rfop, arop = radical_filtration(Aop), knit_ar_quiver(Aop)
    for (i, j) in ar.arrows:
        dr = right_degree(rf, i, j)[0]
        # the dual arrow over A^op (D reverses; matched by dim-vector of D of each end)
        # SHARPEN the D-matching in Step 3; the identity d_r(f) == d_l(D f) is the pin.
        ...
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement.** Adjust-to-reality:
  - **The composite quotient map + the returned witness.** `rad^{d+1}(Z,X_j)/rad^{d+2}`
    — represent by the coordinates of `rad^{d+1}` reduced mod the `rad^{d+2}` basis
    (both are stored column bases in Task 1). `g ↦ [f.matrix @ g_matrix]` reduced
    mod `rad^{d+2}`. **On a hit, return the witness `(z, g_vec, comp_vec)` where
    `comp_vec` = the `_vec` of `f.matrix @ g_matrix` in `Hom(X_z, X_j)`** — so the
    caller asserts all three memberships externally (M3), nothing is hidden inside.
  - **`f_vec` default.** The arrow multiplicity is usually `1`, so `rad/rad²` is
    1-dimensional and `f` is the unique class; when `mult > 1`, iterate the
    `rad/rad²` basis and report per-class (degrees can differ). `degree_table`
    reports the min/list per arrow honestly.
  - **`is_sectional`.** `path` = list of vertex indices; sectional iff no
    consecutive triple `X_{k-1} = τ X_{k+1}` (match `duality.tau(X_{k+1})` to
    `X_{k-1}` by `is_isomorphic`, dim-vector prefilter). `max_sectional_length`
    DFS-enumerates sectional paths on `ar.arrows` (bounded by `N`).
  - **Right degree = the formal dual, implemented directly (M2 — no "flip once").**
    Per the foundation formula, `d_r(f)` sweeps `h ∈ rad^d(X_j, Z)` (the target `Y=X_j`
    on the free side) with the kernel condition `f.then(h) ∈ rad^{d+2}(X_i, Z)` not
    contained in `rad^{d+1}(X_j, Z)`. This is ONE code path parameterized by which
    argument is fixed (`_degree(rf, src, tgt, side)`), not a bookkeeping guess. The
    `oracle_selfcert` `test_left_right_degree_op_symmetry` INDEPENDENTLY confirms
    `d_r(f) == d_l(D f)` over `A^op` — a genuine cross-check of the formula, not the
    definition of the convention.

- [ ] **Step 4: Run tests** — PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ar_invariants.py tests/modules/test_liu_degrees.py
git commit -m "feat(ar-invariants): Liu left/right degrees (fg in rad^{m+2}, corrected) + sectional paths -- decidable sweep on the knit, degree-vs-layer witness certified"
```

---

### Task 3: `ar_invariants.py` — partition + directing + rep-directed + generalized-standard (R21)

**Files:**
- Modify: `src/quiverlab/modules/ar_invariants.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.ar_invariants(...)` delegate)
- Test: `tests/modules/test_ar_invariants.py`

**Interfaces:**
- Consumes: `ARQuiver` (`.vertices`/`.arrows`/`.tau_orbits`/`.is_complete`),
  `duality.tau`/`tau_minus`, `identify_standard`, Task 1 (`RadicalFiltration` for
  the generalized-standard witness + the nonzero-cycle directing refinement).
- Produces:
  ```python
  @dataclass
  class ARInvariants:
      partition     # dict[int -> "postprojective"|"preinjective"|"regular"] by vertex index
      directing     # set[int] of directing vertex indices
      is_representation_directed  # bool (rep-finite & Gamma_A acyclic)
      generalized_standard        # bool (True in rep-finite scope; witness = N)
      nilpotency_index            # int|None (carried from the filtration)
      max_sectional_length        # int
      degrees       # degree_table(rf, ar) (Task 2)
      is_complete   # bool
      status        # mirrors ARQuiver.status
      note          # str
  def ar_invariants(A, *, budget_modules=256, budget_dim=4096) -> ARInvariants
  ```
  `Algebra.ar_invariants(...)` delegates.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_ar_invariants.py
"""AR-component invariants (Plan 57 / R21). Literature: kA_n / D_4 are
representation-directed (Gamma_A acyclic, every indecomposable directing, no
regular modules); the NON-uniform cyclic Nakayama kupisch=[3,2,2] is knittable,
NOT self-injective, and NOT representation-directed (the reachable negative branch,
verified live) yet IS generalized standard (rep-finite => rad^inf = 0) -- the two
flags are independent. Self-cert: partition total; projective postprojective,
injective preinjective; rep-directed <=> acyclic; generalized_standard with the
nilpotency witness. Honest: self-injective (uniform cyclic) / rep-infinite refuse
loudly."""
import pytest

from quiverlab import Quiver, linear_path_algebra, NakayamaAlgebra
from quiverlab.fields import QQ
from quiverlab.families.dynkin import dynkin_quiver     # NOT a top-level export
from quiverlab.modules.ar_invariants import ar_invariants

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_linear_a_n_is_representation_directed(n):
    inv = ar_invariants(linear_path_algebra(n, field=QQ))
    assert inv.is_complete and inv.is_representation_directed
    assert all(p in ("postprojective", "preinjective") for p in inv.partition.values())
    assert "regular" not in inv.partition.values()          # Dynkin: no regular
    assert len(inv.directing) == n * (n + 1) // 2           # every indecomposable directing


@lit
def test_d4_is_representation_directed():
    inv = ar_invariants(dynkin_quiver("D4").algebra(relations=[], field=QQ))
    assert inv.is_complete and inv.is_representation_directed and len(inv.directing) == 12


@selfcert
def test_partition_places_projectives_and_injectives():
    A = linear_path_algebra(3, field=QQ)
    inv = ar_invariants(A)
    # SHARPEN in Step 3: assert P_v vertices are postprojective, I_v preinjective.
    ...


@selfcert
def test_generalized_standard_true_with_nilpotency_witness():
    inv = ar_invariants(linear_path_algebra(3, field=QQ))
    assert inv.generalized_standard is True and inv.nilpotency_index == 3


@lit
def test_cyclic_nakayama_is_knittable_but_NOT_representation_directed():
    # LIVE-VERIFIED WITNESS (the recognizer's reachable negative branch): the
    # non-uniform cyclic Nakayama kupisch=[3,2,2] on Z_3 is NOT self-injective, the
    # knit COMPLETES (7 indecomposables), and Gamma_A has an oriented cycle.
    A = NakayamaAlgebra(kupisch=[3, 2, 2], field=QQ)
    assert A.is_selfinjective() is False                     # verified live
    inv = ar_invariants(A)
    assert inv.is_complete and inv.status == "complete"      # knit accepts it
    assert len(inv.directing) < 7                            # some module lies on a cycle
    assert inv.is_representation_directed is False           # the negative branch
    # the two flags are INDEPENDENT: rep-finite => rad^inf = 0 => generalized standard
    assert inv.generalized_standard is True


@selfcert
def test_self_injective_and_wild_refuse_loudly():
    si = ar_invariants(NakayamaAlgebra(n=3, l=2, cyclic=True, field=QQ))  # uniform => self-inj
    assert si.is_complete is False and si.status == "unsupported"
    wild = ar_invariants(Quiver([1, 2], {"a": (1, 2), "b": (1, 2), "c": (1, 2)})
                         .algebra(relations=[], field=QQ), budget_modules=40)
    assert wild.is_complete is False and wild.status in ("budget", "error")
    assert wild.is_representation_directed is False          # never a false "directed"
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement.** Adjust-to-reality:
  - **Partition:** walk each `tau`-orbit within the knit; `postprojective` if it
    reaches a projective (`identify_standard(...)[0] == "projective"` OR
    `tau(X).dim == 0`), `preinjective` if it reaches an injective (`tau_minus(X).dim
    == 0`), else `regular`. `.tau_orbits` from P41 already groups the orbits — reuse.
  - **Directing / rep-directed:** Tarjan SCC on `ar.arrows` (index digraph); `X`
    directing iff its SCC is a singleton with no self-loop; rep-directed iff the
    whole digraph is acyclic (all SCCs singletons, no loops). **Refine with route
    (ii):** a cycle in `Γ_A` counts only if it carries a nonzero composite — for
    each nontrivial SCC, verify via the `RadicalFiltration` that the round-trip
    composite is nonzero before declaring `X` non-directing (Ringel). On the
    shipped `kA_n`/`D_4` test set `Γ_A` is acyclic so this is vacuous, but the
    refinement is correct and documented.
  - **Not-complete:** if `not ar.is_complete`, return `ARInvariants` with
    `is_representation_directed=False`, `generalized_standard` unknown (`False` with
    a note — never claim `True` off certified scope), `partition={}`,
    `directing=set()`, `status=ar.status`, and the loud note. **Never a false
    verdict.**

- [ ] **Step 4: Run tests** — PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ar_invariants.py src/quiverlab/core/algebra.py tests/modules/test_ar_invariants.py
git commit -m "feat(ar-invariants): postproj/preinj/regular partition + directing modules + rep-directed recognizer + generalized-standard flag -- finite sweeps, honest off-scope refusal"
```

---

### Task 4: GUI story — the `radical_filtration` + `ar_invariants` algebra kinds

**Files:**
- Modify: `src/quiverlab/modules/radical.py` (add `radical_filtration_block(A, budget=512)`)
- Modify: `src/quiverlab/modules/ar_invariants.py` (add `ar_invariants_block(A, budget=512)`)
- Modify: `src/quiverlab/hpc/spec.py` (the budget-carrying **parser** branch in `_parse` (beside `ar_quiver`, `spec.py:253`) AND the `_dispatch` branch (beside `ar_quiver`, `spec.py:1436`) + inline `references`/`citations`)
- Modify: `docs/gui/runner.py` (the twin: the budget-carrying **`_parse_compute`** branch (beside `ar_quiver`, `runner.py:167`) AND the `compute_one` branch (beside `ar_quiver`, `runner.py:835`) calling the SAME builders; the `_snip` recipe map, `runner.py:1184` twin `spec.py:2279`)
- Modify: `docs/gui/gui.js` **and** `webapp/static/gui/gui.js` (BYTE-IDENTICAL: pick-list checkboxes + `renderBlock` branches + `renderRadicalFiltration`/`renderARInvariants` + `SEARCH` entries)
- Modify: `webapp/static/app.js` (`renderResult` branches — server families page)
- Modify: `webapp/server/i18n/{en,es,fr,zh}.json` (new keys ×4) + `webapp/templates/draw.html` (`data-pick-kind-*` wiring)
- Modify: `src/quiverlab/trace/results_html.py` (`_HEADINGS` + `_block_html` branches + `_radical_filtration_html`/`_ar_invariants_html`)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py` (TWO new fixtures; existing byte-identical)
- Test: `tests/webapp/test_radical_filtration_p57.py`, `tests/gui/test_radical_runner_twin.py`

**These are ALGEBRA-block kinds** (they consume the whole knitted category),
wired EXACTLY like the existing `ar_quiver` kind — **NOT** in `MODULE_KINDS`. Each
requires TWO edits per runner (H3): (1) a **budget-carrying parser branch** so
`"radical_filtration"` / `"radical_filtration:512"` parses to a module budget (not
a homological degree, skips `MAX_DEGREE`) — `spec.py::_parse` at the `ar_quiver`
branch (`spec.py:253`, returning `ComputeItem(kind=..., lo=None, hi=budget)`) and
`docs/gui/runner.py::_parse_compute` at the `ar_quiver` branch (`runner.py:167`,
returning `(name, budget)`); (2) a **dispatch branch** — `spec.py::_dispatch`
(`spec.py:1436`) and `docs/gui/runner.py::compute_one` (`runner.py:835`), both
calling the shared builder with `budget=top if top is not None else 512`.
`ar_quiver_block` (`modules/ar.py:732`) is the serializer template.

**KIND-NAME disambiguation (ruling: keep `radical_filtration`).** Exact-match
dispatch is safe and `radical_filtration` is the natural search term, but it must
not be confused with the existing spectral-sequence **`radical_filtration_ss`**
(the Loewy/radical-series spectral sequence, a DIFFERENT object). Add a one-line
disambiguation to the i18n label (`pick.kind.radical_filtration` desc:
"the radical of the module category rad^n(X,Y) and its nilpotency index — NOT the
Loewy radical series (`radical_filtration_ss`)") and a `# NOTE: distinct from
radical_filtration_ss (Loewy series)` comment at BOTH dispatch sites.

- **Shared block builders** (byte-identical across both runners — both import the
  ONE library function; every value JSON-safe, `Module` objects dropped):
  ```python
  # modules/radical.py
  def radical_filtration_block(A, *, budget=512):
      rf = radical_filtration(A, budget_modules=budget)
      block = {"kind": "radical_filtration", "status": rf.status,
               "complete": rf.is_complete, "budget": budget,
               "num_indecomposables": len(rf.indecs),
               "nilpotency_index": rf.nilpotency_index,       # int | null
               "rad_infinity_zero": rf.rad_infinity_zero,     # true | null
               "generalized_standard": rf.is_generalized_standard() if rf.is_complete else None,
               "layer_profile": _layer_totals(rf),            # [sum_ij dim rad^n]_{n>=1} | null off-scope
               "vertices": [{"name": n, "dimvec": {str(v): int(c) for v, c in dv.items()}}
                            for n, dv in zip(rf.names, [X.dimension_vector() for X in rf.indecs])],
               "note": rf.note,
               "latex": r"\operatorname{rad}^{n}(X,Y)\ \text{and}\ N=\min\{n:\operatorname{rad}^n=0\}",
               "references": ["chaio_liu_radical", "liu_degrees", "cmms_radsq"],
               "citations": _citation_pairs(["chaio_liu_radical", "liu_degrees", "cmms_radsq"])}
      # optional pair layers when the request names a module M (matched to a knit vertex)
      return block
  ```
  `ar_invariants_block(A, budget=512)` mirrors it: `kind:"ar_invariants"`,
  `partition` (counts + per-vertex), `representation_directed`, `directing`
  count, `generalized_standard`, `max_sectional_length`, `degrees` (arrow →
  `{d_l, d_r}` with `∞` rendered as `null` + a `finite` flag), `status`/`complete`,
  refs `["liu_degrees", "liu_semistable", "ringel_tame"]`.
- **Off-scope (self-injective/budget)** blocks ship `complete:false`,
  `nilpotency_index:null`, `rad_infinity_zero:null`, and the loud `note` — the GUI
  renders the honest "window / unsupported" banner, never a fake index. A raised
  `QuiverlabError` (char scope) is caught into `{"kind": "...", "error": <message>}`
  (the per-block error precedent), never a 500.

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir;
  copy the `ar_quiver` / `derived_fingerprint` runner-pair fixture pattern):

```python
# tests/webapp/test_radical_filtration_p57.py
"""radical_filtration + ar_invariants algebra kinds: served by hpc.spec, mirrored
byte-identically by the Pyodide twin, honest off-scope banner."""


def test_radical_filtration_block_shape(tmp_path):
    # request: kA3 over QQ, compute ["radical_filtration"]. Assert:
    #   block["complete"] is True; block["nilpotency_index"] == 3;
    #   block["rad_infinity_zero"] is True; block["generalized_standard"] is True;
    #   "chaio_liu_radical" in [k for k, _ in block["citations"]].
    ...


def test_ar_invariants_block_shape(tmp_path):
    # kA3: block["representation_directed"] is True; partition has no "regular";
    #   block["degrees"] is a dict over arrows; "liu_degrees" cited.
    ...


def test_self_injective_off_scope_banner(tmp_path):
    # cyclic Nakayama kZ3/J2: block["complete"] is False, status "unsupported",
    #   nilpotency_index null -- NO fake verdict.
    ...


def test_twin_parity(tmp_path):
    # run the same two requests through docs/gui/runner.py; json.dumps(sort_keys=True)
    # equality on both blocks (shared library builders => byte-identical).
    ...
```

- [ ] **Step 2: Implement** the two shared builders; the **parser branches in BOTH
  parsers** (`spec.py::_parse` at `spec.py:253` and `docs/gui/runner.py::_parse_compute`
  at `runner.py:167`, mirroring the `ar_quiver` budget parse — a request
  `"radical_filtration"` / `"radical_filtration:512"` must reach dispatch with its
  module budget, or it will be rejected as an unknown/degree kind); the
  `spec.py::_dispatch` branches (`if kind == "radical_filtration": block =
  radical_filtration_block(A, budget=budget); return ...` with the
  `# NOTE: distinct from radical_filtration_ss` comment — algebra kinds hand-build
  `references`/`citations` inline, so put the key list in the builder as above);
  the `docs/gui/runner.py::compute_one` twin branches (same builders, inline
  comment `# byte-identical to spec._dispatch`); both `gui.js` copies (checkbox in
  `#qlgui-invariants`, a `-budget` input, `renderBlock` branches keyed on
  `name === "radical_filtration"` / `"ar_invariants"`, helper renderers using
  `matrixGrid` for any small layer matrix and a facts table for the index /
  partition / degrees, `SEARCH` catalog entries); `app.js` `renderResult`
  branches; the four i18n catalogs (`pick.kind.radical_filtration`,
  `pick.kind.ar_invariants`, `block.radical_filtration.title`,
  `block.radical_filtration.index`, `block.radical_filtration.radinf`,
  `block.radical_filtration.window`, `block.ar_invariants.title`,
  `block.ar_invariants.directed`, `block.ar_invariants.partition`,
  `block.ar_invariants.degrees` — **EN + ES + FR + ZH**, English fallback allowed)
  and the `data-pick-kind-radical_filtration` / `data-pick-kind-ar_invariants`
  lines in `draw.html`; and `results_html.py` (`_HEADINGS["radical_filtration"] =
  "Radical filtration"`, `_HEADINGS["ar_invariants"] = "AR-component invariants"`,
  the two `_block_html` branches + helpers — citations auto-render from
  `block["citations"]`).

- [ ] **Step 3: Add TWO golden fixtures** (`radical_filtration_kA3`,
  `ar_invariants_kA3`) to `_runner_goldens.json`; log them in the
  `test_runner_delegation.py` docstring change-log (the `almost_split_a3_s2` ADD
  entry is the model). **Verify existing goldens stay byte-identical BEFORE adding**
  (the `model_dump` override drops absent `module`/`ext_target` blocks, so every
  pre-existing key is unchanged — algebra-only requests carry no module block).

- [ ] **Step 4: Run the gates**

Run: `... -m pytest tests/webapp/test_radical_filtration_p57.py tests/webapp/test_runner_delegation.py tests/gui/test_radical_runner_twin.py tests/gui/test_tau_tilting_gui_wiring.py::test_gui_js_copies_byte_identical tests/hpc -q`
Expected: PASS (both runners byte-identical via the shared builders; the two
`gui.js` copies still byte-identical).

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): radical_filtration + ar_invariants algebra kinds -- nilpotency index / rad^inf / degrees / rep-directed one click, both runners byte-identical, i18n x4, two goldens"
```

---

### Task 5: QPA cross-oracle (probe-first) + honest fallback

**Files:**
- Modify: `src/quiverlab/qpa/scripts.py` / `crosscheck.py` (only if a live verb exists)
- Test: `tests/qpa/test_radical_qpa.py`

**Interfaces:**
- Consumes: the QPA session (`session.should_skip_qpa()`, `session.libgap_handle()`),
  the existing P37 module-Hom crosscheck (`crosscheck.py`, the `hom_glue` verb).
- **QPA has NO module-category-radical or degree surface.** QPA's
  `RadicalOfModule` is the Jacobson radical of a **single module** (`rad M`), NOT
  the category radical `rad(X,Y)`; there is no `LeftDegree`/`RightDegree`/
  `NilpotencyIndexOfRadical` verb. So the QPA leg is (a) a **fail-if-appears**
  probe and (b) the one checkable slice: `dim rad(X_i,X_j) = dim Hom_A(X_i,X_j)`
  for `X_i ≇ X_j` (layer-1, off-diagonal), crosschecked against
  `HomOverAlgebra`.

- [ ] **Step 1: Probe live QPA + Hom crosscheck** (the Plan-35/49 fail-if-appears
  pattern; primary oracles are the mesh/Chaio literature pins + the self-cert
  identities):

```python
# tests/qpa/test_radical_qpa.py
"""QPA probe for a module-category-radical / degree surface (Plan 57). qpa-marked:
skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1. QPA has NO category-radical
or degree verb (RadicalOfModule is the Jacobson radical of ONE module) -- documents
that honestly and FAILS if one appears; the layer-1 dims dim rad(X,Y)=dim Hom(X,Y)
(X !~ Y) are the one QPA-checkable slice via HomOverAlgebra."""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_has_no_category_radical_or_degree_surface():
    lg = session.libgap_handle()
    for name in ("LeftDegreeOfIrreducibleMorphism", "RightDegreeOfIrreducibleMorphism",
                 "NilpotencyIndexOfRadical", "RadicalOfModuleCategory"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), (
            f"QPA now exposes {name} -- wire a real crosscheck (this assert is the "
            "trip-wire that Plan 57's QPA scope note is stale)")


def test_layer1_dims_match_qpa_hom():
    # rad(X_i, X_j) = Hom(X_i, X_j) for X_i !~ X_j: dim-check the off-diagonal
    # layer-1 against QPA HomOverAlgebra on kA_3 indecomposables.
    from quiverlab import linear_path_algebra
    from quiverlab.fields import QQ
    from quiverlab.modules.radical import radical_filtration
    A = linear_path_algebra(3, field=QQ)
    rf = radical_filtration(A)
    for i, Xi in enumerate(rf.indecs):
        for j, Xj in enumerate(rf.indecs):
            if i != j:
                A.crosscheck("hom_glue", Xi, Xj).assert_agree()   # P37 Hom crosscheck
                assert rf.layer_dim(i, j, 1) == A.hom(Xi, Xj)     # rad = Hom off-diagonal
```

- [ ] **Step 2:** If (and only if) the probe finds a live verb, add a
  `crosscheck_radical`/`crosscheck_degree` mirroring `crosscheck_tau`. Otherwise
  the probe + the `dim rad(X,Y) = dim Hom(X,Y)` corroboration IS the QPA leg, and
  the verification page records honestly: "QPA has no module-category-radical or
  degree surface (`RadicalOfModule` is `rad M`, not `rad(X,Y)`); the layer-1
  off-diagonal dims are QPA-checked via `HomOverAlgebra`, the filtration/index/
  degrees are theory + self-cert."

- [ ] **Step 3: Run live** `... -m pytest tests/qpa/test_radical_qpa.py -v`
  (the venv has `[qpa]`). Expected: PASS live (probe absent; Hom agrees).

- [ ] **Step 4: Commit**

```bash
git add src/quiverlab/qpa/ tests/qpa/test_radical_qpa.py
git commit -m "test(qpa): radical/degree probe -- QPA has no category-radical/degree verb (fail-if-appears); layer-1 dims corroborated via HomOverAlgebra"
```

---

### Task 6: verification page, citations, README, metaplan, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick the P57 card)
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`)

- [ ] **Step 1: Citations** (VERIFIED BibTeX only — `_r(key, bibtex_key, kind,
  title, annotation, *tags)`, `registry.py:24`; `bibtex()` hard-fails if `.bib` and
  registry drift). Add to `references.bib` the **verified** entries and register
  them; the `# PIN` entries are added ONLY after the worker BibTeX-verifies the
  author/volume metadata at merge (the P49 precedent):

```bibtex
@article{Liu1992degrees,
  author  = {Liu, Shiping},
  title   = {Degrees of irreducible maps and the shapes of {A}uslander-{R}eiten quivers},
  journal = {Journal of the London Mathematical Society. Second Series},
  volume  = {45}, number = {1}, pages = {32--54}, year = {1992},
}
@article{Liu1993semistable,
  author  = {Liu, Shiping},
  title   = {Semi-stable components of an {A}uslander-{R}eiten quiver},
  journal = {Journal of the London Mathematical Society. Second Series},
  volume  = {47}, number = {3}, pages = {405--416}, year = {1993},
}
@article{ChaioLiu2013,
  author  = {Chaio, Claudia and Liu, Shiping},
  title   = {A note on the radical of a module category},
  journal = {Communications in Algebra},
  volume  = {41}, number = {12}, pages = {4419--4424}, year = {2013},
}
@article{CMMS1994radsq,
  author  = {Coelho, Fl{\'a}vio U. and Marcos, Eduardo N. and Merklen, H{\'e}ctor A. and Skowro{\'n}ski, Andrzej},
  title   = {Module categories with infinite radical square zero are of finite type},
  journal = {Communications in Algebra},
  volume  = {22}, number = {11}, pages = {4511--4517}, year = {1994},
}
@book{Ringel1984tame,
  author    = {Ringel, Claus Michael},
  title     = {Tame Algebras and Integral Quadratic Forms},
  series    = {Lecture Notes in Mathematics}, volume = {1099},
  publisher = {Springer}, year = {1984},
}
```

and in `registry.py` (mirroring `_r("assem_book", "ASS2006", "foundation", ...)`):

```python
_r("liu_degrees", "Liu1992degrees", "foundation",
   "Degrees of irreducible maps and the shapes of Auslander-Reiten quivers",
   "Liu's left/right degrees of irreducible morphisms and their control of the "
   "AR-quiver shape -- the R21 degree theory Plan 57 computes on the knit.", "ar"),
_r("liu_semistable", "Liu1993semistable", "foundation",
   "Semi-stable components of an Auslander-Reiten quiver",
   "Liu's component classification (sectional paths, semistable/directed "
   "components) underlying Plan 57's partition and rep-directed recognizer.", "ar"),
_r("chaio_liu_radical", "ChaioLiu2013", "foundation",
   "A note on the radical of a module category",
   "Chaio-Liu: rep-finiteness through the infinite radical; the nilpotency of "
   "rad(mod A) in the rep-finite case is the maximal depth of composites -- the "
   "theorem behind Plan 57's nilpotency index and rad^inf=0 gate.", "ar"),
_r("cmms_radsq", "CMMS1994radsq", "foundation",
   "Module categories with infinite radical square zero are of finite type",
   "CMMS: (rad^inf)^2 = 0 implies representation-finite -- the class oracle "
   "justifying Plan 57's finite-nilpotency-index rep-finiteness certificate "
   "(documented, not a per-instance decider off the rep-finite domain).", "ar"),
_r("ringel_tame", "Ringel1984tame", "foundation",
   "Tame Algebras and Integral Quadratic Forms",
   "Ringel LNM 1099: directing modules, the postprojective/regular/preinjective "
   "trichotomy -- the R21 component-invariant reference.", "book"),
```

  **`# PIN` (add at merge after BibTeX verification):** Chaio–Guazzelli 2022 (JPAA,
  arXiv:2003.04189 — the nilpotency-index Theorem 1.3/1.5(a) source, cite as
  `chaio_nilpotency_index` if the index formula is surfaced in the report),
  CMMS-1996 cube-zero (J. Algebra), Kerner–Skowroński 1991 (Compositio), Chaio ART
  2018/2023. These are **documentation-only class pins** (rep-infinite, out of the
  certified domain) — cite in the verification-page honest-scope prose, register
  only the ones actually referenced by a shipped block. `ARS1995` (`ars_book`) and
  `ASS2006` (`assem_book`) already exist and cover the AR/knitting background.

- [ ] **Step 2: Verification page.** Add the Plan-57 subsystem rows to the
  subsystem→oracle table (`docs/verification.md:441`):
  - `modules/radical.py` — `oracle_selfcert` (descending filtration, `rad^N=0 &
    rad^{N-1}≠0`, composition closure, layer-1 = Hom off-diagonal, field parity);
    `oracle_crossengine` (mesh route-(i) closed form ≡ route-(ii) on `kA_n`/`D_4`;
    Thm-1.3 index ≡ route-(ii); `dim rad = dim Hom` vs QPA `HomOverAlgebra`);
    `oracle_literature` (`N(kA_n)=n` Thm 1.5(a); `kA_3` layer table; `kA_3/J²`
    index 3 Thm 1.3).
  - `modules/ar_invariants.py` — `oracle_selfcert` (degree-witness triple —
    finite left degree forces the layer drop; partition totality; rep-directed ⇔
    acyclic; generalized-standard witness); `oracle_literature` (`kA_n`/`D_4`
    representation-directed, no regular, all directing; sectional composite
    nonzero; hand-derived small degrees).
  Add the **honest-scope entries** (metaplan §6): (a) certified **only** on the
  rep-finite (knit-complete) domain — `status ∈ {complete,budget,error,
  unsupported}`; (b) **self-injective refused** (`unsupported`) — cyclic Nakayama
  `kZ_n/J^ℓ` out of scope, the Chaio self-injective forms cited as documentation;
  (c) **rep-infinite = window-restricted layers, no nilpotency/`rad^∞` verdict**;
  (d) the **`(rad^∞)²=0 ⇒ rep-finite` battery is a CLASS statement** tested only on
  the rep-finite side + the honest discrimination (finite index vs honest refusal),
  CMMS/KS documentation not per-instance deciders; (e) **QPA cannot compare** the
  category radical or degrees (`RadicalOfModule` is `rad M`) — only the layer-1
  off-diagonal dims via `HomOverAlgebra`; (f) the **general functorial mesh engine
  + automatic standardness detection is DEFERRED** (a named successor) — route (i)
  ships as `kA_n`/`D_4` closed forms + a reachability necessary condition, and the
  route-(i)≡route-(ii) agreement is asserted only on the standard test set.
  Recount the class table (`tests/release/test_oracle_classes.py` drives the
  numbers — run `--collect-only -m <expr>` for each of the five marker classes +
  the union, **paste the LIVE counts**, re-run to green; mid-merge-train honest).

- [ ] **Step 3: README.** One features line: "the radical filtration of `mod A`
  (Liu–Chaio): exact `rad^n(X,Y)` layer dimensions, the nilpotency index of
  `rad(mod A)`, and the `rad^∞ = 0 ⇔ representation-finite` certificate; Liu's
  left/right degrees of irreducible maps, sectional paths, the
  postprojective/preinjective/regular partition, directing modules and the
  representation-directed recognizer — the R21+R37 axis (rep-finite; honest
  window off it)."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/modules/test_radical_filtration.py tests/modules/test_liu_degrees.py tests/modules/test_ar_invariants.py -q` (deep),
  `... -m pytest -q -m deep` (the touched bucket, spot),
  `... -m pytest -q -m fast`,
  `... -m pytest tests/webapp tests/gui tests/hpc -q`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q` — all green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-57 radical-filtration + Liu-degree oracle rows + Liu/Chaio-Liu/CMMS/Ringel citations + honest scope (rep-finite only, self-injective refused, window off-scope, no-QPA) + recounted classes"
```

---

## Acceptance (Plan-57 definition of done)

1. `radical_filtration` + `RadicalFiltration` in `modules/radical.py`;
   `left_degree`/`right_degree`/`is_sectional`/`max_sectional_length`/
   `degree_table` and `ar_invariants` + `ARInvariants` in
   `modules/ar_invariants.py`; the `Algebra.radical_filtration`/`ar_invariants`
   delegates named; every result self-certified or loudly refusing.
2. **R37 core:** exact `dim rad^n(X,Y)` layers by route (ii); the descending
   filtration + termination self-certified (`rad^{n+1}⊆rad^n`, `rad^N=0 &
   rad^{N-1}≠0`, composition closure, `layer-1 = Hom` off-diagonal); the
   **nilpotency index** with `N(kA_n)=n` (Thm 1.5(a)), `kA_3/J²=3` (Thm 1.3), and
   the Thm-1.3 walk ≡ route-(ii) cross-engine; the **`rad^∞=0`** certificate
   (rep-finite ⇒ finite `N` the witness).
3. **R21 degrees:** Liu left/right degrees as a decidable sweep, with the
   **corrected `fg ∈ rad^{m+2}`** definition (the brief's `∉` fixed and documented
   against 1704.03933); the **degree-vs-layer witness triple** certified (finite
   left degree forces the layer drop); sectional composites nonzero.
4. **R21 invariants:** the postprojective/preinjective/regular partition, directing
   modules, the **representation-directed** recognizer (`Γ_A` acyclic — `kA_n`/`D_4`
   rep-directed, no regular, all directing), and the generalized-standard flag
   (True in rep-finite scope, witness `= N`).
5. **Honest semi-decision contract** mirrored from P41: certified iff the knit
   closes; **self-injective refused** (`status="unsupported"` — cyclic Nakayama out
   of scope); **rep-infinite = window-restricted layers, no index/`rad^∞` verdict**;
   the **`(rad^∞)²=0 ⇒ rep-finite`** battery a CLASS statement tested only on the
   rep-finite side with honest discrimination — never a false verdict, never an
   overclaim.
6. `radical_filtration` + `ar_invariants` clickable end-to-end (GUI canvas → block
   → report) in **EN+ES+FR+ZH**, dispatched as algebra kinds in both runners
   byte-identically via the shared `radical_filtration_block`/`ar_invariants_block`,
   the two `gui.js` copies still byte-identical, TWO goldens added with documented
   change-log entries, canonical keys stable (algebra-only requests carry no module
   block — pre-existing keys byte-unchanged), the off-scope banner rendered.
7. QPA probe green (`-m qpa`): no category-radical/degree verb (fail-if-appears),
   the layer-1 off-diagonal dims corroborated via `HomOverAlgebra`.
8. `docs/verification.md` recounted (live numbers, mid-merge-train honest); the
   verified citations (Liu 1992/1993, Chaio–Liu 2013, CMMS 1994, Ringel LNM 1099)
   added and BibTeX-verified, the `# PIN` documentation-only pins added only after
   verification; README line + metaplan P57 card ticked; deep + fast + webapp/gui +
   qpa + release + citations suites green. Honest scope recorded in full.

---

## Change log

- **2026-08-07 implementation note:** `dim rad¹(kA₃) = 9` (ALL radical maps); the
  6 arrows are a basis of `rad/rad²` — the earlier "rad¹ = the six arrows" prose
  describes generators, not `rad¹` itself (the descending filtration puts each
  nonzero map into `rad^m` for every `m ≤ n₀`). Verified live and adjudicated by
  the adversarial critic (who independently re-derived `9` and `d_ℓ(P₁→I₂)=2`).
- **2026-08-07 implementation note (directing/rep-directed):** the plan's "no cycle
  carries a nonzero composite" refinement is WITHDRAWN — on `NakayamaAlgebra([3,2,2])`
  the round-trip composite is `0` in `rad^{≥N}` yet the module is genuinely
  non-directing. The implemented recognizer uses Ringel's definition directly (an
  oriented cycle in `Γ_A` is a cycle of nonzero non-isomorphisms because every
  `Γ_A` arrow is an irreducible map), via a Tarjan-SCC test. Critic-adjudicated in
  the implementation's favor.
- **2026-08-07 (authoring):** initial plan (R21 + R37; route (ii) arbiter, route
  (i) mesh crosscheck; six tasks).
- **2026-08-07 adversarial review: 3 majors + 7 minors applied** (kA₂/kA₃ degree
  tables hand-derived + wrong infinite-degree expectation replaced; non-directed
  witness resolved live or honestly descoped; import/parser/additivity/witness-
  triple/marker fixes). Specifically: **H1/W1** — full left+right Liu degree tables
  for kA₂ and kA₃ hand-derived in the foundation (live-cross-checked against the
  engine's rad structure) and pinned per-arrow by module name; the false
  `test_hereditary_a3_irreducibles_have_infinite_left_degree` stub deleted and
  replaced by `test_liu_degree_tables` + `test_mono_epi_degree_dichotomy_on_ka_n`;
  the boundary-witness insight (`d=N-1 ⇒ rad^{d+2}=0`) documented. **M1** — resolved
  in the affirmative: `NakayamaAlgebra(kupisch=[3,2,2])` is a LIVE-VERIFIED
  knittable (7 indecomposables, `status="complete"`), non-self-injective,
  non-representation-directed witness (`Γ_A` has an oriented cycle) that also
  separates `representation_directed=False` from `generalized_standard=True`;
  adopted as the recognizer's negative-branch test (no descope needed). **M2** —
  right degree stated as the explicit formal dual, "flip once if it fails" removed.
  **M3** — the witness-triple test asserts all three memberships externally
  (`left_degree` now returns `comp_vec`). **B1** — `dynkin_quiver` imported from
  `quiverlab.families.dynkin` (not top-level). **H2** — Krull–Schmidt/Hom-additivity
  justification for summing composites over indecomposables added to the
  foundation. **H3** — Task 4 wiring now names the budget-carrying parser branch in
  BOTH parsers (`spec.py:253`, `runner.py:167`) alongside the dispatch branches.
  **W2** — R21 section reframed to name what is non-degenerate in scope (degrees,
  sectional paths, directing/rep-directed with a realized negative branch); the
  partition's clean trichotomy scoped to directed components. **W3** —
  `test_thm13_index_equals_route_ii` relabeled `oracle_selfcert` (two readings of
  one engine) and a genuinely independent route-(i) `_mesh_layer_dim` crosscheck
  added as `oracle_crossengine`. **W4** — citation-grouping correction recorded
  (nilpotency index = arXiv:2003.04189, not 1704.03933) with a fold-back
  instruction to the research doc. **Kind-name ruling** — `radical_filtration`
  kept, with a one-line i18n/code-comment disambiguation vs `radical_filtration_ss`.
