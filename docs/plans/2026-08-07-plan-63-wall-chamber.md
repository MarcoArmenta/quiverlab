# Plan 63: Wall-and-chamber structure via bricks (R25) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote the **wall-and-chamber structure** of a finite-dimensional algebra
`A = kQ/I` to a first-class, no-code surface, built **via bricks**. For every brick `B`
the **wall** `D(B)` is computed as an **EXACT rational inequality system** over the
submodule dimension-vectors of `B` — `D(B) = {θ : θ·dim B = 0 and θ·dim N ≤ 0 for every
submodule N ⊆ B}` — a rational polyhedral cone of codimension ≥ 1 in the stability space
`K₀(proj A)_ℝ ≅ ℝ^{Q₀}`. The **chambers** are the g-vector cones of the support
τ-tilting pairs (the Plan-45 exchange graph). The **whole structure is a fan** (BST 2019):
chambers ↔ maximal g-cones, walls = their codim-1 boundaries, and `#chambers = #support
τ-tilting`. The construction is **certified complete iff `A` is brick-finite ⟺
τ-tilting-finite** (DIJ — decidable via the Plan-45 exchange-graph BFS); otherwise it
returns a **bounded region with honest truncation** (the P62 discipline: the discovered
sub-fan is exact as far as it goes, never claimed complete, no count asserted). One
algebra-level compute kind — `wall_chamber` — exposes the walls (inequality systems + the
brick per wall), the chambers (exact integer g-cone generators), the chamber↔wall
adjacency, the four counts, and a **2D drawing for rank ≤ 3** (a stereographic/affine slice
of the fan, exact rational coordinates server-side, pixels client-side) — a table
otherwise — end-to-end in all four locales.

**What is genuinely NEW over Plan 45 (the scope boundary — read this first).** Plan 45
already ships `wall_and_chamber_fan(A, budget)`: chambers (g-cones) and **walls-as-exchange
-edges** (a per-edge `normal` = the brick dim-vector, i.e. only the wall's HYPERPLANE
direction). For kA₂ that is **5 edge-facets**, and the fan is drawn for `n ∈ {2,3}` only,
returning **empty when the BFS is incomplete**. Plan 63 adds four things Plan 45 does not
have:

1. **The wall `D(B)` as a polyhedral CONE, not a hyperplane normal.** `D(B)` is the
   half-space-restricted stability locus — the inequality system cut by the submodule
   dim-vectors. For `P₁` over kA₂ this is a **RAY** `θ₁+θ₂=0 ∧ θ₂≤0` (direction `(1,-1)`),
   **not** the full line — Plan 45 only recorded the line's normal `(1,1)`. This is the
   headline object.
2. **Walls for ALL ranks.** The inequality-system payload is dimension-agnostic; rank ≤ 3
   also gets the drawing, rank ≥ 4 gets the inequality tables. (Plan 45's geometry stops at
   `n = 3`.)
3. **Walls grouped by brick — "5 chambers / 3 walls".** The three distinct bricks of kA₂
   give **3 walls** (`D(S₁)`, `D(S₂)` full lines; `D(P₁)` a ray), each a union of the
   exchange-edge facets that carry that brick label. (Plan 45 counts 5 edges.)
4. **The honest bounded region for brick-infinite `A`.** Plan 45's fan returns empty when
   incomplete; Plan 63 returns the discovered sub-fan with `complete=False`, a `truncation`
   note, and **no** count — the record's "bounded region with honest truncation".

**Architecture.** One new module `src/quiverlab/tautilting/wallchamber.py`, a thin
exact-rational layer over the merged Plan-45 machinery (no new math engines):

- `wallchamber.py::Wall` — a frozen value object for a single wall `D(B)`: the brick, its
  dim-vector + standard name, the equality normal `dim B`, the deduped submodule inequality
  set, `is_full_hyperplane`, `codim`, and (for `n ≤ 2`) the exact extreme rays.
- `wallchamber.py::wall_of_brick(B, *, budget)` — builds `Wall` from
  `stability._submodule_dimvecs(B)` (the shipped exact submodule enumerator). Reuses the
  shipped King convention (`θ·dim N ≤ 0`) so the wall is byte-consistent with
  `is_theta_semistable`.
- `wallchamber.py::wall_chamber_structure(A, *, budget)` — the full payload: chambers from
  `exchange_graph(A)` g-matrices, walls (one per distinct brick, via `bricks(A)` +
  `wall_of_brick`), chamber↔wall adjacency (exchange edges grouped by brick label), the four
  counts (complete case only), `render ∈ {"fan2d","fan3d","table"}`, the rank ≤ 3 drawing
  geometry, and the honest `complete`/`status`/`truncation` contract.
- `wallchamber.py::_is_green_path` (module-internal, oracle helper) — a maximal green
  sequence read as a monotone chamber path (source `(A,0)` → sink `(0,A)` across downward
  walls). Records "MGS = green paths".
- `Algebra.wall_chamber_structure(budget_pairs=512)` — thin lazy-import delegate beside
  `Algebra.exchange_graph` (`core/algebra.py`).

**Plan 45 is left byte-unchanged.** `wall_and_chamber_fan`, `is_theta_semistable`,
`_submodule_dimvecs`, `bricks`, `exchange_graph`, `tau_tilting_block` are consumed
read-only; Plan 63 calls `exchange_graph` directly (so it can render the bounded region
even when the BFS is incomplete — the fan cannot). No Plan-45 test moves.

**Tech Stack.** Pure exact rational geometry via `fractions.Fraction` (the Plan-45
`stability.py` precedent — the fan already ships `Fraction` strings; the JS renderer does
the only float conversion, `docs/gui/gui.js`, exempt). Submodule enumeration is the shipped
`stability._submodule_dimvecs` (exact over the `Domain`). No floats in `src/` (AST-gated by
`tests/test_no_floats.py`). Exact 2D orientation predicates (cross-product signs) for the
`n = 2` extreme-ray / tiling logic — never `atan2`.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **P45 (the C4 τ-tilting engine) is a HARD prerequisite, MERGED to `dev`.** This plan
  consumes, read-only: `tautilting.mutation.exchange_graph`/`ExchangeGraph`/`wall_normal`
  (the BFS + the honest `is_complete`/`status`/`n_regular` contract),
  `tautilting.stability._submodule_dimvecs`/`is_theta_semistable`/`is_theta_stable`/`_det`/
  `_l1_project`/`wall_and_chamber_fan`, `tautilting.torsion.bricks`/`hasse_orientation`,
  `tautilting.green.maximal_green_sequences`, `modules.hom.identify_standard`/
  `is_isomorphic`. **If P45 is not merged when a worker picks this up, STOP and escalate** —
  do not fork the exchange-graph or submodule machinery. (Wave-1 plans P51–P53+P55+P57+P58
  are also merged into `dev`; P63 is independent of them — no edge taken.)
- **Buckets are auto-assigned by directory (`tests/conftest.py`): `tests/modules/` →
  deep.** All Plan-63 engine tests live in `tests/modules/` as
  `tests/modules/test_wall_chamber*.py` (deep — they share the exchange-graph budget with
  the Plan-45 τ-tilting suite). `tests/qpa/` → qpa; `tests/webapp/`, `tests/gui/` → fast
  (extras-gated dirs, unmarked per the Plan-32 rule). Run new tests by path during
  development; finish each task with a `-m deep` spot-run of the touched files.
- **THE ENGINE RUNS OVER QQ BY DEFAULT** (the load-bearing Plan-45 char caveat). `bricks`,
  `is_isomorphic`, `identify_standard`, `decompose` and the exchange-graph BFS are rigorous
  only over **char 0 or char > dim** (Dickson/CIW). Over QQ the brick enumeration and the
  chamber count are decisive; the batteries pin over **QQ**. Over `char ≤ dim` the engine
  inherits the loud `QuiverlabError` refusal unchanged — never a silent wrong wall set.
- **Bricks decide over the implicit algebraically-closed / char-0 base.** A brick is
  `end_dim(B) == 1`; the wall count `= #bricks` reads correctly only when `End_A(B) = k`.
  Batteries pin over QQ; the GF(pⁿ) proper-division-ring caveat (`dim_k End(B) > 1`) is
  stated honestly on the verification page (inherited from Plan 45).
- **Self-injective input is IN SCOPE (M3).** The surface makes NO hereditary assumption — it
  reads only bricks + the exchange graph, which the Plan-45 engine computes for any presented
  f.d. algebra over QQ. The corrected KT Ex. 15 algebra (cyclic rad²-Nakayama N₃²) and
  kZ₂/rad² are BOTH self-injective, and their tests (14 chambers; 4 walls with P₁, P₂ opposite
  half-rays) double as the self-injective coverage. Bricks decide; no Auslander–Reiten /
  hereditary special-casing.
- **Honest semi-decision contract (metaplan §6; the P62 discipline).** The wall-and-chamber
  structure is **certified complete iff `A` is brick-finite ⟺ τ-tilting-finite** (DIJ).
  Completeness is decided by the Plan-45 exchange-graph BFS: `eg.is_complete` ⟺
  τ-tilting-finite within `budget_pairs`. When incomplete the payload carries
  `complete=False`, `status="budget"`, a `truncation` string, and **omits the counts**
  (`counts=None`, `green_count=None`) — a bounded region, never a claimed-complete fan and
  never a count that would be a partial lie. The discovered chambers and walls are each
  EXACT (a real chamber, a real wall); only totality/counts are withheld.
- **No floats in `src/`.** All geometry is exact `Fraction`. The inequality systems and the
  extreme rays ship as exact fraction/int data; the JS renderer does the only
  fraction→pixel conversion (`gui.js`, exempt). The `n = 2` ray-vs-line decision uses exact
  arithmetic (sign of a `Fraction`), never `atan2`.
- **Composition is left-to-right**; `_assert_comparable`-style guards inherited from Plan
  45's module calls. All refusals are `QuiverlabError(message, hint=...)`
  (the two-space `[hint: ...]` convention, `errors.py`). Presentation-less
  (structure-constant) algebras: the exchange graph / bricks already refuse loudly (the
  quiver is required) — Plan 63 surfaces that refusal as a clean typed 4xx, never a 500.
- **Every citation key is BibTeX-verified before use** (Plan-29 rule); `_r(...)` registry
  entry + `references.bib` entry both land, gated by `tests/citations/test_bib_structure.py`.
- New scalar/algebra kind acts on the **existing algebra block** ⇒ **schema stays v1**
  (the `tau_tilting`/`ar_quiver` precedent). **Canonical keys are request-derived**
  (`cache.py::canonical_key`, sha256 over the sorted-keys blob + `library_version()`): a
  request that does not use the new `wall_chamber` key keys byte-identically to before, so
  existing goldens are untouched. Verify by running `test_runner_delegation.py` BEFORE and
  AFTER adding the golden.
- Plan-32 markers: certificate/identity/dimension tests (`is_full_hyperplane`, codim,
  `θ·dim B = 0` on the wall, facet-shared-by-exactly-2-chambers, the `D(B)`-inequalities
  hold on the wall's rays) = `oracle_selfcert`; the kA₂ hand-derivation (5 chambers / 3
  walls, `D(P₁)` a ray), the DIJ brick-finite ⟺ τ-tilting-finite gate, the
  Kaipel–Treffinger Ex. 13/15 worked examples = `oracle_literature`; `#chambers = #sτ-tilt`, and
  the inequality-defined `D(B)` rays ≡ the grouped exchange-edge facets (two independent
  constructions) = `oracle_crossengine`; **QPA has NO wall-and-chamber surface** — the `qpa`
  bucket is an honest-scope probe that FAILS if QPA ever ships one.
- **Mid-merge-train counts:** the metaplan lands ~19 subplans in overlapping waves, so
  absolute suite counts drift. **Task 6 recounts the oracle-class table at merge time by
  running `tests/release/test_oracle_classes.py`** (paste-the-live-numbers, never a guessed
  count) and claims only the deltas this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table
  green) and the README line. Conventional commits; green at every commit; branch
  `plan-63-wall-chamber` off `dev`.

---

## Record (verbatim — the mandate for this plan)

> **R25 — Wall-and-chamber for all ranks via bricks.** [D-scout P7; keep-with-corrections]
> Object: walls D(B) from brick enumeration (exact linear inequalities over submodule
> dim-vectors), chambers = g-vector cones. GATE (critic): certified complete iff
> τ-tilting-finite/brick-finite — decidable; otherwise a bounded region with honest
> truncation. Refs: Brüstle–Smith–Treffinger arXiv:1805.01880 (Adv. Math.);
> Kaipel–Treffinger arXiv:2302.12699; Asai arXiv:1610.05860. Oracles: chamber count = #
> support τ-tilting (cross-engine); kA₂ = 5 chambers / 3 walls (verified — the P₁ wall is a
> ray); MGS = green paths. Size M.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P63.

---

## Reference re-verification (mandatory; findings recorded)

All statements below were re-read from the primary sources during authoring (WebFetch +
`pdftotext -layout` on the arXiv PDFs of BST 1805.01880v2, Kaipel–Treffinger 2302.12699v1,
Asai 1610.05860, DIJ 1503.00285; cross-checked against ar5iv HTML). Quotations are verbatim
from those extractions. The **empirical ground truth** (the kA₂ picture) was independently
confirmed by running the merged Plan-45 engine during authoring (recorded in Methodology).

- **The stability space, wall, chamber — VERIFIED VERBATIM (Brüstle–Smith–Treffinger,
  *Wall and Chamber Structure for finite-dimensional Algebras*, Adv. Math. **354** (2019)
  106746, arXiv:1805.01880).** BST draw `θ` and `[M]` in the **same** `ℝⁿ ≅ K₀(A)_ℝ`
  (`n` = #simples), pairing via the standard inner product `θ(M) := ⟨θ, dim M⟩` (NOT the
  dual). *Def. 3.1* (King [King 1994]): "a non-zero module `M` is **θ-stable** if … `θ(M) =
  0`, and `θ(L) < 0` for every proper submodule `L` of `M`. … `M` … is **θ-semistable** if
  `θ(L) ≤ 0` for every submodule `L` of `M`." *Def. 3.2*: "The stability space of `M` is
  `D(M) = {θ ∈ ℝⁿ : M is θ-semistable}`. … `D(M)` is a cone given by intersections of
  hyperplanes … We say `D(M)` … is a **wall** when `D(M)` has codimension one." BST intro:
  "`D(M)` … is contained in the hyperplane orthogonal to `θ` [i.e. `θ(dim M)=0`], but it
  could have smaller dimension." *Def. 3.3*: a **chamber** is a connected component of
  `ℝⁿ \ ⋃_{M} D(M)`.
- **Sign / sub-vs-quotient convention — VERIFIED, double-sourced.** BST *Def. 3.1* and KT
  *Def. 3* both use **King with submodules and `≤ 0`**: `θ·dim M = 0 ∧ θ·dim L ≤ 0 ∀ L ⊆ M`.
  KT state the equivalence explicitly (verbatim): "`⟨v, dim L⟩ ≤ 0` for any nonzero proper
  subobject `L` of `M` **or, equivalently**, `⟨v, dim N⟩ ≥ 0` for all nonzero proper
  quotients `N` of `M`." So the correct pairings are **(submodule, ≤ 0)** or **(quotient,
  ≥ 0)** — never "(quotient, ≤ 0)". The merged Plan-45 `stability.is_theta_semistable` uses
  **exactly BST/KT's (submodule, ≤ 0)**; Plan 63 reuses it verbatim, so the wall `D(B)` is
  byte-consistent with the shipped, tested predicate AND matches the papers character for
  character. (Confirmed in read-only Python: `_submodule_dimvecs(P₁) = {(0,0),(0,1),(1,1)}`
  ⇒ `D(P₁) = {θ₁+θ₂=0, θ₂≤0}` = the ray `(1,-1)`.)
- **Walls come from bricks — ATTRIBUTION CORRECTED (do not miscite BST).** BST prove one
  direction (*Example 3.4*, verbatim): "the dimension vector of a θ-stable module is a
  **brick** in mod A for every `θ`." BST do **NOT** have a single numbered theorem "the
  codim-1 walls are exactly the `D(B)` for bricks `B`". That packaging is
  **Demonet–Iyama–Jasso** *Thm 1.4* (τ-tilting finite ⟺ finitely many bricks) + *Thm 1.5*
  (bijection indecomposable **τ-rigid** modules ↔ bricks `X` with `T(X)` functorially
  finite) **together with** BST's "θ-stable ⇒ brick". Cite it as **DIJ Thm 1.4/1.5 + BST
  Ex. 3.4**, never as one BST theorem.
- **Chambers ↔ support τ-tilting pairs — VERIFIED VERBATIM (BST).** *Thm 1.2 / Cor. 3.29*:
  "there is an **injective** function `C` mapping the τ-tilting pair `(M, P)` onto a chamber
  `C_{(M,P)}` … if `A` is τ-tilting finite then `C` is also **surjective**." *Prop. 3.15*:
  the interior of the positive g-cone of `(M,P)` is a chamber. *Cor. 3.18*: a τ-tilting pair
  "induces a chamber … having **exactly `n` walls** `{D(N₁),…,D(Nₙ)}`". **The wall-grouping
  is BST *Remark 3.19*, verbatim:** "Not every wall is generated by the positive cone of
  some almost τ-tilting pair… In general **one wall `D(N)` can be made of more than one
  facet**, see … `D(S(2))` in Figure 2." This is exactly Plan 63's grouping of exchange-edge
  facets into one brick-wall `D(B)`. Hence **`#chambers = #support τ-tilting pairs`** (the
  cross-engine oracle: chamber count `= len(exchange_graph(A).vertices)`).
- **Asai *Semibricks* — SCOPE CORRECTED (indexing only, NOT the fan).** Asai
  (arXiv:1610.05860, IMRN 2020 no. 16, 4993–5054) *Def. 1.1* defines brick/semibrick; *Thm
  1.3(2)* gives the bijection `sτ-tilt A → f_L-sbrick A`; *Prop. 1.6* gives `sτ-tilt A →
  f-tors A`. **The words "g-vector", "fan", "cone", "wall", "chamber", "TF equivalence" do
  NOT appear anywhere in `Semibricks`.** Asai supplies the **brick/semibrick indexing** that
  LABELS walls and chambers — do **not** cite 1610.05860 for the geometric fan or for TF
  equivalence (those are in other Asai/Asai–Iyama work).
- **Brick-finite ⟺ τ-tilting-finite; the fan-completeness gate — VERIFIED (DIJ,
  arXiv:1503.00285, IMRN 2019 (3) 852–892).** *Thm 1.4*: "`A` is τ-tilting finite iff there
  are only finitely many … bricks." *Thm 1.2*: τ-tilting finite ⟺ every torsion class is
  functorially finite. *Cor. 2.9*: ⟺ finitely many indecomposable τ-rigid / basic support
  τ-tilting / f.f. torsion classes. *Thm 1.7*: when τ-tilting finite, the g-vector
  simplicial complex `Δ(A)` is homeomorphic to `S^{n-1}` (the fan is complete). **The
  decidable gate:** the Plan-45 BFS closes (`is_complete=True`) iff τ-tilting-finite within
  budget; else the structure is incomplete and Plan 63 returns the bounded region. **Note:**
  "the wall-and-chamber structure is a complete fan ⟺ τ-tilting-finite" is the COMBINATION
  **DIJ Thm 1.7** (`Δ(A) ≅ S^{n-1}`) **+ BST Cor. 3.29** (chambers ↔ pairs bijectively iff
  τ-tilting-finite) — not one verbatim theorem; cite the pair.
- **kA₂ = 5 chambers / 3 walls — VERIFIED VERBATIM (KT *Example 13*).** Kaipel–Treffinger
  compute, verbatim, for `A = kQ`, `Q: 1 → 2`: `D(1) = {v : ⟨v,(1,0)⟩=0} = {(0,y)}`;
  `D(2) = {v : ⟨v,(0,1)⟩=0} = {(x,0)}`; `D(1/2) = {v : ⟨v,(1,1)⟩=0 and ⟨v,(0,1)⟩≤0} =
  {(x,−x) : 0 ≤ x}`, "the last line contains two conditions, since … there also exists a
  non-trivial submodule `2 ↪ 1/2`." So (KT notation `1`=S₁, `2`=S₂, `1/2`=P₁):
  - `D(S₁) = {θ₁ = 0} = {(0,y)}` — **full line** (S₁ simple).
  - `D(S₂) = {θ₂ = 0} = {(x,0)}` — **full line** (S₂ simple).
  - `D(P₁) = {θ₁+θ₂ = 0, θ₂ ≤ 0} = {(x,−x) : x ≥ 0}` — the **RAY** (direction `(1,-1)`,
    θ₁≥0, θ₂≤0): P₁'s proper submodule S₂=(0,1) fires the extra inequality `θ₂≤0`, cutting
    the line in half. **This is the record's "the P₁ wall is a ray", now a verbatim KT pin.**
  The two axes' 4 quadrants + the ray bisecting the fourth ⇒ **5 chambers = |sτ-tilt kA₂| =
  Catalan(3) = 5** (BST Cor. 3.29; Plan-45-pinned). (KT print the equations + Figure 1.1;
  the integer "5" is read off the figure + the sτ-tilt count — `# PIN` the literal prose
  "5", the wall EQUATIONS are verbatim.) The engine confirms: 5 pairs, 5 edge-facets with
  brick labels `{S₁:2, S₂:2, P₁:1}`; grouping by brick gives **3 walls**.
- **KT *Example 15* = 14 chambers is the CYCLIC rad²-Nakayama N₃² (`1→2→3→1` mod rad²), NOT
  the linear one — engine-verified, a rank-3 oracle.** KT compute the wall-and-chamber
  structure of a radical-square-zero Nakayama algebra whose module list (Ex. 15/32) contains
  **`3/1`** (top S₃, socle S₁) — a module that requires an arrow `3→1`. So KT's quiver is the
  **3-CYCLE** `Q: 1→2→3→1` modulo the square of the arrow ideal (the self-injective Nakayama
  N₃², where P₁=`1/2`, P₂=`2/3`, P₃=`3/1`), and they conclude verbatim "the wall-and-chamber
  structure for `Q` consists of **14 chambers**." The algebra is
  `Quiver([1,2,3], {"a":(1,2),"b":(2,3),"c":(3,1)}).algebra(relations=["a*b","b*c","c*a"],
  field=QQ)`; the merged Plan-45 engine gives **14 support τ-tilting pairs (complete,
  n_regular), 6 bricks** — VERIFIED by running the engine on this exact algebra during the
  fix round (recorded in Methodology). **Provenance correction (H2):** an earlier draft
  mis-transcribed KT's algebra as the LINEAR `1→2→3` mod rad², which is a DIFFERENT algebra —
  the engine gives it **12** support τ-tilting pairs (5 bricks), NOT 14 (hand-count agrees).
  Per the metaplan's reference re-verification rule (run the engine on the pin's exact algebra
  at spec time) the pin is now against the cyclic algebra. The rank-3 wall equations (the
  length-2 brick walls as half-spaces, the three simple walls as full 2-planes) are **to be
  derived at implementation** from the cyclic algebra's submodule dim-vectors — the earlier
  draft's `D(1/2) = {(x,−x,z) : x ≥ 0}` was fabricated against the wrong algebra and is
  deleted. Do **not** call the cyclic count "= Catalan(4)": that is a numerical coincidence
  with hereditary kA₃'s Catalan(4)=14, which is a separate fact (see the linear kA₃/rad²
  entry). Plan 63's Task 3 pins the cyclic count as a rank-3 `lit` oracle
  (`num_chambers == 14`, `render=fan3d`).
- **Linear kA₃/rad² = 12 support τ-tilting pairs — engine-verified `oracle_crossengine`
  value, NOT attributed to KT.** `Quiver([1,2,3], {"a":(1,2),"b":(2,3)}).algebra(
  relations=["a*b"], field=QQ)` (the linear radical-square-zero Nakayama) has **12** chambers
  and **5** bricks, confirmed by the merged engine during the fix round. This is recorded as
  an honest consistency value (`num_chambers == len(exchange_graph(A).vertices) == 12`), NOT a
  KT number and NOT Catalan(4). (Hereditary kA₃ — `linear_path_algebra(3)`, no relations — is
  the algebra whose count IS Catalan(4)=14; that separate fact drives the `#chambers =
  #sτ-tilt` cross-engine parametrize.)
- **Maximal green sequences as green paths — VERIFIED (the Plan-45 `green.py` surface).**
  `maximal_green_sequences(A)` returns `{"sequences": [[pair_id,…]], "count", …}`; each
  sequence is a monotone directed path `(A,0) → (0,A)` in the Hasse-oriented exchange graph,
  i.e. a sequence of chambers each crossing one wall downward — **a green path** through the
  wall-and-chamber structure. kA₂ has exactly **2** (the two sides of the pentagon;
  Plan-45-pinned). Plan 63 records the path (not just the count) and a small verifier
  `_is_green_path` that each MGS is a source→sink monotone chamber path.
- **Kaipel–Treffinger, arXiv:2302.12699 — title + status VERIFIED.** Exact title:
  *Wall-and-chamber structures for finite-dimensional algebras and τ-tilting theory*
  (Kaipel, Treffinger — ICRA 2022 lecture notes, a proceedings chapter). Abstract
  (verbatim): "The wall-and-chamber structure is a geometric invariant that can be
  associated to any algebra. In this notes we give the definition of this object and we
  explain its relationship with torsion classes and τ-tilting theory." *Def. 11* (chamber =
  OPEN connected component of `ℝⁿ \ ⋃ D(M)`) matches BST *Def. 3.3*. It is the source of the
  Ex. 13 (kA₂) and Ex. 15 (14-chamber) oracles above. `# verify` the final journal
  reference if the proceedings volume has appeared by merge; the arXiv id is stable.

---

### Task 0: citations (do first — later tasks reference the keys)

**Files:**
- Modify: `src/quiverlab/citations/references.bib`, `src/quiverlab/citations/registry.py`
- Test: `tests/citations/test_bib_structure.py` (existing gate)

**Interfaces:** the `_r(key, bibtex_key, kind, title, annotation, *tags)` registry helper
(`registry.py:24`). `demonet_iyama_jasso` (→ `DIJ2019`) and `king_stability` (→ `King1994`)
**already exist** (Plan 45, `registry.py:488`,`493`) — reuse, do NOT re-add.
`_citation_pairs` (`hpc/spec.py`) resolves keys → `[key, formatted]` pairs for payloads.

- [ ] **Step 1: add the keys** (BibTeX-verify each before committing — Plan-29 rule):
  - `brustle_smith_treffinger` → Brüstle, Smith, Treffinger, *Wall and chamber structure
    for finite-dimensional algebras*, Adv. Math. **354** (2019) 106746, arXiv:1805.01880,
    doi 10.1016/j.aim.2019.106746. (the wall `D(M)`, the fan theorem, walls = `D(B)` for
    bricks)
  - `asai_semibricks` → Asai, *Semibricks*, Int. Math. Res. Not. IMRN (2020) no. 16,
    4993–5054, arXiv:1610.05860, doi 10.1093/imrn/rny150. (semibricks ↔ f.f. torsion
    classes ↔ support τ-tilting — the brick/semibrick INDEXING; NOT the geometric fan) —
    # verify the exact IMRN volume/pages.
  - `kaipel_treffinger` → Kaipel, Treffinger, *Wall-and-chamber structures for
    finite-dimensional algebras and τ-tilting theory*, arXiv:2302.12699 (ICRA 2022 lecture
    notes / proceedings chapter). (Examples 13 + 15 are Task 3's `lit` oracles) — # verify
    the final proceedings reference if it has appeared by merge.

  Registry entries (mirroring `_r("king_stability", "King1994", …)`). **The annotations
  reflect the corrected attributions (reference section): BST = the D(M)/chamber definitions
  + chambers↔τ-tilting pairs; DIJ = walls-from-bricks + the finiteness gate; Asai = the
  categorical INDEXING only, NOT the fan.**

```python
_r("brustle_smith_treffinger", "BST2019", "foundation",
   "Wall and Chamber Structure for finite-dimensional Algebras",
   "Brustle-Smith-Treffinger: the wall D(M) = {theta : M theta-semistable} (Def 3.1-3.3), "
   "chambers <-> support tau-tilting pairs (Thm 1.2 / Cor 3.29), one wall = many facets "
   "(Rem 3.19) -- the ground truth for Plan 63. NB walls=D(brick) is DIJ, not BST.",
   "tau-tilting", "stability"),
_r("asai_semibricks", "Asai2020", "foundation",
   "Semibricks",
   "Asai: semibricks <-> functorially finite torsion classes <-> support tau-tilting "
   "modules (Thm 1.3 / Prop 1.6) -- the brick/semibrick INDEXING that labels walls and "
   "chambers (Plan 63). Contains no g-vector/fan/wall/chamber content -- do not cite for "
   "the geometry.",
   "tau-tilting", "stability"),
_r("kaipel_treffinger", "KT2023", "context",
   "Wall-and-chamber structures for finite-dimensional algebras and tau-tilting theory",
   "Kaipel-Treffinger: the definition + torsion-class/tau-tilting relationship, with the "
   "worked kA2 (Ex 13: D(P1) a ray) and cyclic rad^2-Nakayama N3^2 (Ex 15: 14 chambers) "
   "examples -- Plan 63 "
   "literature oracles.",
   "tau-tilting", "stability"),
```

- [ ] **Step 2:** `... -m pytest tests/citations/test_bib_structure.py -q` green (registry ↔
  bib in sync; every `_r` key resolvable via `bibtex(key)`).
- [ ] **Step 3: Commit** — `docs(citations): wall-and-chamber references (Brustle-Smith-Treffinger 1805.01880, Asai Semibricks, Kaipel-Treffinger 2302.12699)`

---

### Task 1: `wallchamber.py` — the wall `D(B)` as an exact inequality system (THE headline)

**Files:**
- Create: `src/quiverlab/tautilting/wallchamber.py`
- Modify: `src/quiverlab/tautilting/__init__.py` (export `Wall`, `wall_of_brick`,
  `wall_chamber_structure` — the last added in Task 2)
- Test: `tests/modules/test_wall_chamber_walls.py`

**Interfaces:**
- Consumes: `stability._submodule_dimvecs(M, budget)` (the shipped exact submodule
  enumerator, returns a set of vertex-order dim-vector tuples),
  `stability.is_theta_semistable` (the King cross-check),
  `modules.hom.identify_standard` (the `S_v`/`P_v`/`I_v` naming).
- Produces:
  ```python
  @dataclass(frozen=True)
  class Wall:
      brick_dimvec: tuple          # dim B in vertex order (the equality NORMAL): theta . dim B = 0
      brick_name: str | None       # "S1"/"P2"/... via identify_standard, else None
      equality: tuple              # == brick_dimvec (the hyperplane {theta . dim B = 0})
      inequalities: tuple          # sorted tuple of dim-vector tuples d with theta . d <= 0,
                                   #   = submodule dim-vectors of B, MINUS the trivial 0 and
                                   #   the full dim B (both implied by the equality)
      is_full_hyperplane: bool     # True iff `inequalities` is empty (simple brick / no
                                   #   proper nontrivial submodule constraint) => D(B) is the
                                   #   whole hyperplane {theta . dim B = 0}
      codim: int                   # 1 (D(B) lies in a codim-1 hyperplane; the cone within
                                   #   may be lower-dimensional -- codim records the ambient)
      rays: tuple | None           # n <= 2 ONLY: the exact extreme rays of D(B) as
                                   #   vertex-order Fraction-string vectors ("p/q"); a full
                                   #   line = two opposite rays, a ray = one. None for n >= 3
                                   #   (the drawing groups exchange-edge facets -- Task 3).

  def wall_of_brick(B, *, budget=4096) -> Wall
      # D(B) = {theta : theta . dim B = 0 and theta . dim N <= 0 for every submodule N <= B}.
      # Reads _submodule_dimvecs(B); drops the 0 vector and the full dim B (implied by the
      # equality); is_full_hyperplane iff nothing remains; rays computed for n <= 2. Loud
      # QuiverlabError if _submodule_dimvecs exceeds its budget (large brick).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_wall_chamber_walls.py
"""The wall D(B) of a brick as an EXACT inequality system (Plan 63 / R25). D(B) =
{theta : theta.dim B = 0, theta.dim N <= 0 for every submodule N <= B}. Literature (BST
2019): D(B) is a rational polyhedral cone of codim >= 1; a simple brick gives the full
hyperplane, a non-simple brick a proper face. THE PIN: over kA2, D(S1)/D(S2) are full lines
but D(P1) is a RAY (the record's 'the P1 wall is a ray')."""
import pytest
from fractions import Fraction

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.tautilting.stability import is_theta_semistable
from quiverlab.tautilting.wallchamber import Wall, wall_of_brick

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@lit
def test_simple_brick_wall_is_full_hyperplane():
    A = _kA2()
    w1 = wall_of_brick(A.simple(1))
    assert w1.brick_dimvec == (1, 0) and w1.equality == (1, 0)
    assert w1.is_full_hyperplane is True and w1.inequalities == ()
    assert w1.codim == 1 and w1.brick_name == "S1"
    # full line {theta_1 = 0}: two opposite rays along the theta_2 axis
    assert set(map(tuple, [[Fraction(r[0]), Fraction(r[1])] for r in w1.rays])) == {
        (Fraction(0), Fraction(1)), (Fraction(0), Fraction(-1))}


@lit
def test_P1_wall_is_a_ray():
    # THE record pin. P1 = (1,1), unique proper submodule S2 = (0,1); D(P1) = {theta_1 +
    # theta_2 = 0, theta_2 <= 0} = the single ray (1,-1). NOT the full line.
    A = _kA2()
    w = wall_of_brick(A.projective(1))
    assert w.brick_dimvec == (1, 1) and w.equality == (1, 1)
    assert w.is_full_hyperplane is False           # a proper face -- a RAY
    assert (0, 1) in w.inequalities                # the S2 submodule constraint theta_2 <= 0
    rays = [(Fraction(r[0]), Fraction(r[1])) for r in w.rays]
    assert len(rays) == 1                          # a ray, not a line (one extreme ray)
    d = rays[0]
    assert d[0] > 0 and d[1] < 0 and d[0] + d[1] == 0   # direction (1,-1) up to scale


@selfcert
def test_wall_points_are_theta_semistable():
    # every theta in the relative interior of D(B) makes B theta-semistable (King) -- the
    # cross-check tying the inequality system to the shipped stability predicate.
    A = _kA2()
    B = A.projective(1)
    w = wall_of_brick(B)
    # a point on the ray (1,-1): theta = (1,-1)
    assert is_theta_semistable(B, [Fraction(1), Fraction(-1)]) is True
    # a point on the OPPOSITE ray (-1,1) is NOT in D(P1) (theta_2 = 1 > 0): not semistable
    assert is_theta_semistable(B, [Fraction(-1), Fraction(1)]) is False
    # the equality holds on every ray; every inequality holds on every ray
    for r in w.rays:
        theta = [Fraction(r[0]), Fraction(r[1])]
        assert sum(theta[i] * w.equality[i] for i in range(2)) == 0
        for d in w.inequalities:
            assert sum(theta[i] * d[i] for i in range(2)) <= 0


@selfcert
def test_codim_and_equality_normal():
    A = _kA2()
    for M in (A.simple(1), A.simple(2), A.projective(1)):
        w = wall_of_brick(M)
        dv = M.dimension_vector()
        assert w.equality == tuple(dv[v] for v in (1, 2))
        assert w.codim == 1
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.tautilting.wallchamber`
- [ ] **Step 3: Implement `wallchamber.py`.**
  - `wall_of_brick(B)`: `verts = list(B.algebra.quiver.vertices)`;
    `dv = tuple(B.dimension_vector()[v] for v in verts)`; `subs =
    stability._submodule_dimvecs(B, budget=budget)`; `ineqs = sorted(d for d in subs if
    any(d) and d != dv)` (drop `0` and the full `dim B`); `is_full_hyperplane = not ineqs`;
    `name = _std_name(B)` via `identify_standard`.
  - **The `n = 2` extreme rays (exact, no floats).** The hyperplane `θ·dv = 0` in `ℝ²` is a
    line through the origin with direction `t = (dv[1], -dv[0])` (rotate the normal 90°;
    exact ints). The line = `{s·t : s ∈ ℝ}`. Restrict by the inequalities: for each `d ∈
    ineqs`, `θ·d ≤ 0` becomes `s·(t·d) ≤ 0`, i.e. a sign constraint on `s`. Collect the
    feasible sign set: if no inequality constrains the sign (all `t·d == 0`), both signs are
    feasible ⇒ **full line** (rays `+t` and `-t`, primitivised). If some `t·d > 0` forces
    `s ≤ 0` and some `t·d < 0` forces `s ≥ 0` simultaneously ⇒ only `s = 0` (the origin, a
    degenerate wall — should not happen for a genuine brick; assert/raise). Otherwise a
    single feasible half ⇒ **one ray** (`+t` or `-t`). Primitivise each ray to the smallest
    integer direction (gcd), ship as `("p/q",…)` fraction strings (here integers). Use exact
    `Fraction`/`int` throughout.
  - **`n ≥ 3`: `rays = None`** (the drawing groups exchange-edge facets in Task 3; the
    inequality system is the all-rank payload). **`n = 1`:** the only stability space is
    `ℝ`; `dv = (d,)` with `d ≥ 1` ⇒ `θ·dv = 0` forces `θ = 0` ⇒ `D(B) = {0}`, `rays = ()`,
    `is_full_hyperplane = False`, `codim = 1`. (A degenerate single-vertex edge case; keep
    it honest — no crash.)
- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(tautilting): wall_of_brick + Wall -- D(B) as an exact submodule-inequality system, n<=2 extreme rays, the kA2 P1-ray pin`

**Adjust to reality (Task 1):**
- Confirm `_submodule_dimvecs` returns vertex-ORDER tuples matching `M.dimension_vector()`'s
  vertex order (it does — both iterate `M.algebra.quiver.vertices`); the `test_codim`
  equality-normal assertion is the arbiter.
- The `n=2` rotate-normal-90° direction `(dv[1], -dv[0])` must be checked against the
  submodule signs; the `test_P1_wall_is_a_ray` direction `(1,-1)` (with `t·(0,1) = -dv[0] =
  -1 < 0` forcing `s ≥ 0`) is the arbiter. Fix the sign bookkeeping until it holds, never
  the pin.
- Reuse of the private `stability._submodule_dimvecs` is intra-package and intentional; if a
  reviewer wants a public name, add a one-line public `submodule_dimvecs` re-export in
  `stability.py` and import that — behaviour byte-identical.

---

### Task 2: `wall_chamber_structure` — chambers + brick-walls + the honest gate

**Files:**
- Modify: `src/quiverlab/tautilting/wallchamber.py` (add `wall_chamber_structure`),
  `src/quiverlab/tautilting/__init__.py`, `src/quiverlab/core/algebra.py`
  (`Algebra.wall_chamber_structure(budget_pairs=512)` lazy delegate beside
  `Algebra.exchange_graph`)
- Test: `tests/modules/test_wall_chamber_structure.py`

**Interfaces:**
```python
def wall_chamber_structure(A, *, budget=512) -> dict:
    # {"n": n, "complete": bool, "status": "complete"|"budget",
    #  "render": "fan2d"|"fan3d"|"table",
    #  "num_chambers": int, "num_walls": int,
    #  "chambers": [ {"id": i, "g_matrix": [[..]], "rays": [["p/q",...],...],  # g-cone gens
    #                 "label": str, "support": [v,...], "is_initial": bool} , ... ],
    #  "walls": [ {"id": k, "brick_dimvec": {v:m}, "brick_name": str|None,
    #              "equality": [ints], "inequalities": [[ints],...],
    #              "is_full_hyperplane": bool, "codim": 1,
    #              "rays": [["p/q",...],...] | None,          # n<=2 (Task 3 fills n=3)
    #              "facets": [[i,j],...]} , ... ],            # exchange edges carrying this brick
    #  "counts": {"chambers": k, "walls": k, "bricks": k, "s_tau_tilt": k} | None,
    #  "green_count": int | None,
    #  "truncation": str | None,
    #  "references": [...]}
# Algebra.wall_chamber_structure(budget_pairs=512) delegates.
```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_wall_chamber_structure.py
"""The wall-and-chamber structure via bricks (Plan 63 / R25). Literature: kA2 = 5 chambers
/ 3 walls (the P1 wall a ray); DIJ brick-finite <=> tau-tilting-finite gate. Cross-engine:
#chambers = #support tau-tilting. Self-cert: each exchange-edge facet is shared by exactly
two chambers; every wall is coplanar (theta.dim B = 0 on its facets)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@lit
def test_kA2_five_chambers_three_walls():
    A = _kA2()
    wc = A.wall_chamber_structure()
    assert wc["complete"] and wc["status"] == "complete"
    assert wc["num_chambers"] == 5 and wc["num_walls"] == 3     # 3 bricks -> 3 walls
    assert wc["render"] == "fan2d"
    # exactly one wall is a proper ray (D(P1)); the other two are full lines (D(S1),D(S2))
    rays = [w for w in wc["walls"] if not w["is_full_hyperplane"]]
    lines = [w for w in wc["walls"] if w["is_full_hyperplane"]]
    assert len(rays) == 1 and len(lines) == 2
    assert sorted(rays[0]["brick_dimvec"].items()) == [(1, 1), (2, 1)]   # the P1 ray


@xeng
@pytest.mark.parametrize("n, chambers", [(2, 5), (3, 14)])   # HEREDITARY kA_n: Catalan(n+1)
def test_chamber_count_equals_support_tau_tilting(n, chambers):
    # linear_path_algebra(n) is HEREDITARY kA_n (no relations) -> Catalan(n+1) chambers. This
    # is a DIFFERENT algebra/fact from the rad^2-Nakayamas (linear = 12, cyclic = 14 chambers).
    from quiverlab.tautilting.mutation import exchange_graph
    A = linear_path_algebra(n, field=QQ)
    wc = A.wall_chamber_structure()
    eg = exchange_graph(A)
    assert wc["complete"]
    assert wc["num_chambers"] == len(eg.vertices) == chambers
    assert wc["counts"]["chambers"] == wc["counts"]["s_tau_tilt"] == chambers


@selfcert
def test_each_facet_bounds_exactly_two_chambers():
    # the fan property, FALSIFIABLY (W2): the grouped facets are in BIJECTION with the
    # exchange graph's edges -- every exchange edge appears as exactly one facet under exactly
    # one wall, every facet is a real edge between two DISTINCT chambers. A grouping that
    # dropped/duplicated an edge, or invented a non-edge facet, fails here (the old version
    # only checked i != j, which is tautological for an exchange edge).
    from quiverlab.tautilting.mutation import exchange_graph
    A = linear_path_algebra(3, field=QQ)
    wc = A.wall_chamber_structure()
    eg = exchange_graph(A)
    edges = {tuple(sorted(e)) for e in eg.arrows}           # eg.arrows keyed by (i<j)
    facets = [tuple(sorted(f)) for w in wc["walls"] for f in w["facets"]]
    for (i, j) in facets:
        assert i != j and 0 <= i < wc["num_chambers"] and 0 <= j < wc["num_chambers"]
    assert len(facets) == len(set(facets)) == len(edges)    # bijection: no drop, no dup
    assert set(facets) == edges                             # every facet is a real edge


@selfcert
def test_walls_are_coplanar_theta_dot_dimB_zero_on_facets():
    # every facet grouped under a wall lies in the hyperplane theta . dim B = 0: the facet
    # vector (a shared g-vector) is orthogonal to the brick dim-vector.
    A = _kA2()
    wc = A.wall_chamber_structure()
    verts = [1, 2]
    for w in wc["walls"]:
        nrm = [w["brick_dimvec"][v] for v in verts]
        # the wall's own extreme rays (n<=2) satisfy the equality + all inequalities
        for r in (w["rays"] or []):
            from fractions import Fraction
            theta = [Fraction(x) for x in r]
            assert sum(theta[k] * nrm[k] for k in range(2)) == 0


@lit
def test_brick_infinite_is_bounded_region_not_complete():
    # the 2-Kronecker is tau-tilting-infinite / brick-infinite (DIJ): the structure cannot
    # close -> honest truncation, a bounded region, NO counts.
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    wc = K.wall_chamber_structure(budget_pairs=40)
    assert wc["complete"] is False and wc["status"] == "budget"
    assert wc["counts"] is None and wc["green_count"] is None
    assert wc["truncation"]                          # a non-empty honest note
    assert wc["num_chambers"] > 0                    # but the bounded region IS shown


def _kZ2_radsq():
    # self-injective Nakayama N2^2: 1 <-> 2, BOTH length-2 paths killed (rad^2 = 0).
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


@lit
def test_kZ2_radsq_four_walls_two_share_a_hyperplane():
    # THE primary discriminating oracle (BST Remark 3.19). kZ2/rad^2 is self-injective with 6
    # support tau-tilting pairs and 4 bricks: S1=(1,0), S2=(0,1) simple, and TWO
    # non-isomorphic bricks P1, P2 BOTH of dim-vector (1,1) with submodule-dimvec sets
    # {(0,0),(1,0),(1,1)} and {(0,0),(0,1),(1,1)} -> D(P1), D(P2) are OPPOSITE half-rays of
    # the ONE line theta_1 + theta_2 = 0. A grouping keyed on dim-vector (up to scaling) would
    # MERGE them into 3 walls; the correct ISO-CLASS grouping keeps 4. (Engine-verified in the
    # fix round: 6 pairs, 4 bricks, the two submodule sets above.)
    from fractions import Fraction
    A = _kZ2_radsq()
    wc = A.wall_chamber_structure()
    assert wc["complete"] and wc["num_chambers"] == 6
    assert wc["num_walls"] == 4                       # NOT 3 -- P1, P2 stay separate
    # exactly two walls carry the dim-vector (1,1); both lie in the ONE line theta_1+theta_2=0
    dim11 = [w for w in wc["walls"] if sorted(w["brick_dimvec"].items()) == [(1, 1), (2, 1)]]
    assert len(dim11) == 2
    for w in dim11:
        assert list(w["equality"]) == [1, 1]         # the shared hyperplane normal
        assert w["is_full_hyperplane"] is False      # each is a proper RAY (half of the line)
        assert len(w["rays"]) == 1                    # one extreme ray each (a half-ray)
    # their extreme rays are the two OPPOSITE directions (1,-1) and (-1,1), NOT collapsed
    rays = {tuple(int(Fraction(x)) for x in w["rays"][0]) for w in dim11}
    assert rays == {(1, -1), (-1, 1)}                 # distinct opposite half-rays


@selfcert
@pytest.mark.parametrize("build", [
    _kA2,                                             # thin: kA2, 3 bricks (distinct dimvecs)
    _kZ2_radsq,                                       # non-thin: kZ2/rad^2, 4 bricks (P1,P2)
])
def test_counts_consistency_walls_equal_bricks(build):
    # M2 counts arbiter. ONE canonical brick set: torsion.bricks(A) -- the shipped
    # ISO-CLASS-deduped enumerator (chosen over an inline dim-vector dedup precisely because
    # it disambiguates same-dim-vector bricks by torsion-class membership, torsion.py:142, so
    # kZ2/rad^2's P1,P2 count as 2). In the complete case #walls == #bricks, and the walls'
    # brick dim-vector MULTISET matches torsion.bricks(A)'s (on kZ2 the (1,1) key has
    # multiplicity 2 -- a dim-vector-merged grouping would show 1 here and fail).
    from collections import Counter
    from quiverlab.tautilting.torsion import bricks as torsion_bricks
    A = build()
    wc = A.wall_chamber_structure()
    assert wc["complete"]
    ref = torsion_bricks(A)
    assert wc["counts"]["walls"] == wc["counts"]["bricks"] == len(ref) == wc["num_walls"]
    wall_dvs = Counter(tuple(sorted(w["brick_dimvec"].items())) for w in wc["walls"])
    ref_dvs = Counter(tuple(sorted(B.dimension_vector().items())) for B in ref)
    assert wall_dvs == ref_dvs
```

- [ ] **Step 2: Run to verify failure**, implement, iterate.
- [ ] **Step 3: Implement `wall_chamber_structure`.**
  - `eg = exchange_graph(A, budget_pairs=budget)` (called ONCE; thread it into the helpers to
    avoid the Plan-45 block's repeated BFS — see adjust-to-reality).
  - **Chambers** (always, even when incomplete): one per `eg.vertices[i]` — `g_matrix`,
    `label`, `support`, `is_initial`, and `rays` = the g-vector columns as `Fraction`
    strings (the Plan-45 fan idiom, `stability.wall_and_chamber_fan:223`).
  - **Walls — one per brick ISO-CLASS** (the canonical brick set is `torsion.bricks(A,
    budget)`, the shipped iso-class-deduped `Module` list; it is the arbiter for
    `counts["walls"] == counts["bricks"]`). For each brick `B`: `wall = wall_of_brick(B)`;
    attach `facets` = the exchange edges `(i,j)` whose `torsion._edge_brick(...)` is **iso to
    `B`** (`is_isomorphic`). **TRAP — group by ISO-CLASS, NEVER by dim-vector (or `wall_normal`
    up to scaling).** On a NON-thin algebra two non-isomorphic bricks can share a dim-vector:
    kZ₂/rad² has P₁, P₂ **both `(1,1)`**, on OPPOSITE half-rays of the ONE line `θ₁+θ₂=0` (BST
    Rem. 3.19). A dim-vector key (or a primitive-`wall_normal` key) MERGES them → 3 walls; the
    correct answer is **4**. The shipped `torsion._edge_brick` (`torsion.py:142`) already
    disambiguates same-dim-vector bricks by torsion-class membership (the DIRRT label), and
    `torsion.bricks` already iso-dedups (`torsion.py:223`) — reuse both; do NOT re-key on
    `brick_dimvec`. (An optional single-BFS optimization dedups the `eg.arrows` edge-bricks
    inline via `_edge_brick` + `is_isomorphic` — same iso-class grain — see adjust-to-reality;
    the DEFAULT is `torsion.bricks(A, budget)` for correctness + the M2 counts arbiter.)
  - **`render`** = `"fan2d"` (n==2), `"fan3d"` (n==3), `"table"` (n>=4 or n==1).
  - **Counts** (complete case only): `{"chambers": len(eg.vertices), "walls":
    len(walls), "bricks": len(bricks), "s_tau_tilt": len(eg.vertices)}`. When incomplete:
    `counts=None`, `green_count=None`, `truncation="A appears to be tau-tilting-infinite /
    brick-infinite (the exchange-graph BFS did not close within budget N); the region shown
    is the sub-fan explored from (A,0) -- each chamber and wall is exact, but the structure
    is NOT complete and no count is claimed (DIJ: brick-finite <=> tau-tilting-finite)."`.
  - **`green_count`** = `maximal_green_sequences(A, cap=budget)["count"]` (complete only).
  - `references = ["brustle_smith_treffinger", "demonet_iyama_jasso", "asai_semibricks",
    "king_stability", "kaipel_treffinger"]` (Task 3 pins the KT Ex. 13/15 oracles, so KT
    ships in the payload).
  - **The char caveat**: `bricks`/`is_isomorphic` refuse loudly off char-scope; let the
    `QuiverlabError` propagate (the compute kind catches it into a typed error block, Task 5).
  - `Algebra.wall_chamber_structure(budget_pairs=512)`: lazy-import delegate
    (`from quiverlab.tautilting.wallchamber import wall_chamber_structure; return
    wall_chamber_structure(self, budget=budget_pairs)`), beside `Algebra.exchange_graph`.
- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(tautilting): wall_chamber_structure + Algebra delegate -- chambers=g-cones, walls=D(B) per brick, DIJ brick-finite gate, honest bounded-region truncation`

**Adjust to reality (Task 2):**
- **Canonical brick set = `torsion.bricks(A, budget)` (the DEFAULT — chosen for the M2
  counts arbiter).** `torsion.bricks`, `semibricks`, `hasse_orientation` and
  `maximal_green_sequences` each call `exchange_graph` internally (the Plan-45 block accepts
  the recompute). Build the walls one-per-brick from `torsion.bricks(A, budget)`: it already
  iso-dedups (so kZ₂/rad²'s P₁, P₂ count as 2) and it is what `test_counts_consistency_walls_
  equal_bricks` compares against, so `counts["walls"] == counts["bricks"] == len(bricks(A))`
  holds by construction. **Optional single-BFS optimization:** compute `eg` once and derive
  the found bricks from `eg.arrows` inline (via `torsion._edge_brick(A, eg.vertices[i]["pair"],
  eg.vertices[j]["pair"], eg.arrows[(i,j)]["brick"], _torsion_universe(A, budget))` + an
  `is_isomorphic` dedup). This is allowed ONLY if it preserves the exact ISO-CLASS grain
  (`_edge_brick` + `is_isomorphic`, never a dim-vector key) so the counts arbiter still
  passes; note whichever path you took in the docstring.
- **The `facets` grouping arbiter** is `test_kA2_five_chambers_three_walls` (3 walls, facet
  multiset `{S1:2, S2:2, P1:1}`) AND `test_kZ2_radsq_four_walls_two_share_a_hyperplane` (the
  NON-THIN discriminator — 4 walls, P₁ and P₂ must NOT collapse) + `test_each_facet_bounds_
  exactly_two_chambers` (facets ↔ edges bijection). If two same-dim-vector bricks merge, the
  grouping is keying on the dim-vector / `wall_normal` instead of the ISO-CLASS — switch to
  `_edge_brick` + `is_isomorphic` (never fix by loosening the pin). If a brick's facets come
  out empty, the edge→brick iso match is wrong.
- **Incomplete-case walls.** When `eg.is_complete is False`, still build walls for the found
  edge-bricks (a bounded set) — this is the honest bounded region. Do NOT assert
  `#walls == #bricks(A)` in the incomplete case (there is no complete brick count).

---

### Task 3: rank ≤ 3 drawing geometry + MGS green paths + literature/consistency oracles

**Files:**
- Modify: `src/quiverlab/tautilting/wallchamber.py` (n=3 wall geometry via grouped facets +
  the L1 projection; `_is_green_path`)
- Test: `tests/modules/test_wall_chamber_geometry.py`,
  `tests/modules/test_wall_chamber_oracles.py`

**Interfaces:**
- Consumes: `stability._l1_project` (the shipped octahedron projection, `stability.py:185`),
  `stability.wall_and_chamber_fan` (the n≤3 chamber geometry cross-check),
  `green.maximal_green_sequences`, `torsion.hasse_orientation`.
- Produces (in `wall_chamber_structure`, for `n ∈ {2,3}`): each wall gains its drawing
  geometry — for `n=2` the `Wall.rays` from Task 1; for `n=3` the `rays` = the grouped
  exchange-edge facet vectors (each a shared g-vector on the wall), each L1-projected
  (`rays_l1`/`faces`/`net2d`), with a self-cert that every facet vector satisfies the
  inequality system. Chambers reuse the fan's L1 projection for `n=3`.
- `_is_green_path(eg, orient, seq) -> bool` — a maximal green sequence `seq` (list of
  pair-ids) is a monotone downward chamber path `(A,0) → (0,A)` crossing one wall per step.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_wall_chamber_geometry.py
"""Rank <= 3 drawing geometry + MGS green paths (Plan 63 / R25). Self-cert: the wall's
drawing rays satisfy its inequality system (theta.dim B = 0, theta.dim N <= 0); n=3 rays
carry an L1/octahedron projection; each MGS is a monotone source->sink chamber path."""
import pytest
from fractions import Fraction

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


@xeng
def test_n2_wall_rays_equal_grouped_exchange_facets():
    # TWO INDEPENDENT CONSTRUCTIONS (W1): the inequality-defined D(B) extreme rays (Task 1)
    # must equal, as SETS OF PRIMITIVE DIRECTIONS, the grouped exchange-edge facet vectors of
    # the Plan-45 fan (each facet a shared g-vector lying IN the wall). Not merely "both
    # nonempty" -- genuine set equality per brick. kA2's bricks have DISTINCT dim-vectors, so
    # the fan facets can be grouped by dim-vector here. Concretely: S1 -> {(0,1),(0,-1)},
    # S2 -> {(1,0),(-1,0)}, P1 -> {(1,-1)}.
    from math import gcd
    from quiverlab.tautilting.stability import wall_and_chamber_fan

    def prim(v):                                     # primitive integer direction (sign kept)
        ints = [int(Fraction(x)) for x in v]
        g = 0
        for x in ints:
            g = gcd(g, abs(x))
        return tuple(x // g for x in ints) if g else tuple(ints)

    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)
    wc = A.wall_chamber_structure()
    fan = wall_and_chamber_fan(A)
    fan_by_brick = {}
    for f in fan["walls"]:
        key = tuple(sorted(f["brick_dimvec"].items()))
        fan_by_brick.setdefault(key, set()).update(prim(v) for v in f["facet"])
    ineq_by_brick = {tuple(sorted(w["brick_dimvec"].items())): {prim(r) for r in w["rays"]}
                     for w in wc["walls"]}
    assert set(ineq_by_brick) == set(fan_by_brick)       # same brick set
    for key in ineq_by_brick:
        assert ineq_by_brick[key] == fan_by_brick[key]   # SET EQUALITY of directions
    assert ineq_by_brick[((1, 1), (2, 1))] == {(1, -1)}  # the P1 ray, single direction


@selfcert
def test_n3_walls_project_via_L1_and_lie_in_D_of_B():
    # n=3 self-cert (W2): each wall's drawing rays (grouped exchange-edge facet vectors)
    # actually LIE IN D(B) -- theta.dim B = 0 AND theta.dim N <= 0 for every submodule
    # inequality -- and carry the L1/octahedron projection. FALSIFIABLE: a facet violating an
    # inequality, or a wall carrying no rays at all, fails here (the old version's disjunction
    # was trivially satisfiable by w["rays"] being None).
    A = linear_path_algebra(3, field=QQ)             # hereditary kA3, n=3, tau-tilting-finite
    wc = A.wall_chamber_structure()
    assert wc["render"] == "fan3d"
    saw_wall_ray = False
    for w in wc["walls"]:
        rays = w["rays"] or []
        nrm = [w["brick_dimvec"][v] for v in (1, 2, 3)]
        assert len(w.get("rays_l1", rays)) == len(rays)          # one projection per ray
        for r in rays:
            saw_wall_ray = True
            theta = [Fraction(x) for x in r]
            assert sum(theta[k] * nrm[k] for k in range(3)) == 0          # in the hyperplane
            for d in w["inequalities"]:
                assert sum(theta[k] * d[k] for k in range(3)) <= 0        # in D(B)
    assert saw_wall_ray                                          # walls DID carry drawing rays
    for ch in wc["chambers"]:
        assert "net2d" in ch and ch["net2d"]          # chambers carry the L1 net position


@lit
def test_kA2_two_maximal_green_sequences_are_green_paths():
    from quiverlab.tautilting.wallchamber import _is_green_path
    from quiverlab.tautilting.mutation import exchange_graph
    from quiverlab.tautilting.torsion import hasse_orientation
    from quiverlab.tautilting.green import maximal_green_sequences
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)
    eg = exchange_graph(A)
    orient = hasse_orientation(eg)
    mgs = maximal_green_sequences(A)
    assert mgs["count"] == 2
    for seq in mgs["sequences"]:
        assert _is_green_path(eg, orient, seq)        # each MGS is a monotone chamber path
```

```python
# tests/modules/test_wall_chamber_oracles.py
"""Consistency + literature oracles (Plan 63 / R25). #chambers = #support tau-tilting on
the tau-tilting-finite zoo (cross-engine); the Kaipel-Treffinger Ex. 13 (kA2) + Ex. 15
(CYCLIC rad^2-Nakayama N3^2, 14 chambers) worked examples (lit); the linear kA3/rad^2 = 12
cross-engine value (NOT KT); the DIJ brick-finite gate."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ

xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


@xeng
@pytest.mark.parametrize("A", [
    linear_path_algebra(2, field=QQ),                                    # hereditary kA2 (5)
    linear_path_algebra(3, field=QQ),                                    # hereditary kA3 (14)
    Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(              # LINEAR rad^2-Nakayama
        relations=["a*b"], field=QQ),                                    # 12 pairs (NOT KT Ex.15)
])
def test_chamber_count_equals_exchange_graph(A):
    from quiverlab.tautilting.mutation import exchange_graph
    wc = A.wall_chamber_structure()
    if wc["complete"]:
        assert wc["num_chambers"] == len(exchange_graph(A).vertices)


@xeng
def test_linear_kA3_radsq_is_12_not_14():
    # H1 pin: the LINEAR 1->2->3 mod rad^2 (radical-square-zero Nakayama) has 12 support
    # tau-tilting pairs -- NOT 14. Engine-verified in the fix round (12 pairs, 5 bricks). This
    # is an honest cross-engine consistency value, NOT attributed to KT (KT Example 15 is the
    # CYCLIC N3^2, tested separately). Guards against the earlier mis-transcription that
    # pinned this algebra at 14.
    from quiverlab.tautilting.mutation import exchange_graph
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ)
    wc = A.wall_chamber_structure()
    assert wc["complete"]
    assert wc["num_chambers"] == len(exchange_graph(A).vertices) == 12   # NOT 14


@lit
def test_kaipel_treffinger_example_13_kA2_walls():
    # KT arXiv:2302.12699 Example 13 (VERBATIM): D(S1)={(0,y)} full line, D(S2)={(x,0)}
    # full line, D(P1)={(x,-x): x>=0} the RAY. 3 walls, 5 chambers.
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)
    wc = A.wall_chamber_structure()
    assert wc["num_chambers"] == 5 and wc["num_walls"] == 3
    rays = [w for w in wc["walls"] if not w["is_full_hyperplane"]]
    assert len(rays) == 1 and sorted(rays[0]["brick_dimvec"].items()) == [(1, 1), (2, 1)]


@lit
def test_kaipel_treffinger_example_15_cyclic_radsq_nakayama_14_chambers():
    # KT Example 15 is the CYCLIC radical-square-zero Nakayama N3^2 (Q: 1->2->3->1 mod rad^2,
    # self-injective): its module list contains 3/1 (top S3, socle S1), which needs the arrow
    # 3->1, so KT's algebra is the 3-CYCLE -- NOT the linear 1->2->3 (that one has 12 chambers,
    # see test_linear_kA3_radsq_is_12_not_14). Engine-verified in the fix round: 14 pairs
    # (complete, n_regular), 6 bricks. This test ALSO covers self-injective input (M3). A
    # rank-3 oracle (render fan3d). (NOT "= Catalan(4)": a numerical coincidence, not the
    # hereditary-kA3 Catalan fact.)
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}).algebra(
        relations=["a*b", "b*c", "c*a"], field=QQ)   # 3-cycle, all length-2 paths = 0 (rad^2)
    wc = A.wall_chamber_structure()
    assert wc["complete"] and wc["num_chambers"] == 14 and wc["render"] == "fan3d"
```

- [ ] **Step 2: Run to verify failure**, implement, iterate.
- [ ] **Step 3: Implement.**
  - **n=3 wall geometry.** For each wall, `rays` (n=3) = the DISTINCT facet vectors of its
    grouped exchange edges (each is a shared g-vector lying on the wall). L1-project each via
    `stability._l1_project` (exact fractions), attaching `rays_l1`/`faces`/`net2d`.
    Self-cert: every facet vector `θ` satisfies `θ·dim B = 0` and `θ·d ≤ 0` for each `d` in
    the wall's inequalities. (Directly computing the 2D extreme rays of a 3D wall cone from
    the inequality system is deferred — the grouped facets ARE the wall's rays in the
    complete fan, and the self-cert proves they lie in `D(B)`; recorded as honest scope.)
  - **Chambers (n=3):** reuse the fan's `_l1_project` on each g-vector column (mirror
    `wall_and_chamber_fan:227-231`).
  - **`_is_green_path(eg, orient, seq)`:** check `eg.vertices[seq[0]]["is_initial"]`; the
    last is the terminal `(0,A)` pair (support == all vertices); each consecutive
    `(seq[t], seq[t+1])` is an exchange edge oriented DOWN by `orient`; length `n` steps.
- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(tautilting): rank<=3 wall drawing geometry (n=3 L1 grouped facets) + MGS green-path verifier + #chambers=#stau-tilt oracles`

**Adjust to reality (Task 3):**
- The `n=3` direct extreme-ray computation from the inequality system is OPTIONAL (a 2D
  vertex enumeration inside the hyperplane). The shipped design uses the grouped facets +
  the `θ ∈ D(B)` self-cert, which is sound and matches the complete fan. If a reviewer wants
  the direct rays, add an exact 2D-cone solver inside the hyperplane (Fourier–Motzkin on 2
  free coordinates) — the grouped-facets result is the arbiter.
- **The Kaipel–Treffinger oracles (provenance, H2).** Ex. 13 (kA₂ wall equations, `D(P₁)` a
  ray) is transcription-checked verbatim AND engine-reproduced. Ex. 15 is the **CYCLIC**
  rad²-Nakayama N₃² (`Q: 1→2→3→1` mod rad², relations `["a*b","b*c","c*a"]`); its provenance
  is **engine-verified** (14 support τ-tilting pairs, run on this exact algebra during the fix
  round) PLUS the KT-PDF module-list identification (the module `3/1` forces the arrow `3→1`,
  so KT's quiver is the 3-cycle). **An earlier draft mis-transcribed the algebra as the linear
  `1→2→3` mod rad², which the engine gives 12 — not 14.** Per the metaplan's reference
  re-verification rule (run the engine on the pin's EXACT algebra at spec time) the pin is now
  against the cyclic algebra, and both counts (linear 12, cyclic 14) are engine-confirmed. On
  first run verify the relation strings kill every length-2 path; the count is the arbiter.
  The literal "5 chambers" for kA₂ is read off KT Figure 1.1 + the sτ-tilt count (the wall
  EQUATIONS are verbatim) — the pin asserts the count, not the figure.

---

### Task 4: QPA honest-scope skip

**Files:**
- Create: `tests/qpa/test_wall_chamber_qpa.py`

**Interfaces:** QPA 1.37 has **no** wall-and-chamber / brick-wall / stability-fan surface
(no `WallAndChamber`, `StabilitySpace`, `BrickWall`, `ChamberOfAlgebra` verb; and it has no
τ-tilting surface at all — the Plan-45 finding). Probe live via `NamesGVars()` (the
`tests/qpa/test_products_qpa.py` / Plan-45 precedent): the test SKIPS honestly and **FAILS
if QPA ever ships** any such verb, so the honest-scope claim on the verification page cannot
silently rot.

- [ ] **Step 1: Write the probe** (`pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
  reason=...)`; sweep `NamesGVars()` for the candidate names; `pytest.skip(...)` if absent,
  `assert False` (loud) if any appears). Also record that `#support τ-tilting` itself has no
  QPA cross-check (Plan-45 stated it), so `#chambers` inherits no QPA oracle — the external
  cross-check named is FD-Applet / DIJ tables, not a live QPA call.
- [ ] **Step 2: Run** — `... -m pytest tests/qpa/test_wall_chamber_qpa.py -q -m qpa` (venv
  has `[qpa]`): expect a clean skip.
- [ ] **Step 3: Commit** — `test(qpa): honest no-wall-and-chamber-surface probe (fails if QPA ships one)`

---

### Task 5: no-code exposure — the `wall_chamber` compute kind

The organizing invariant is **byte-identity across the two runners**, achieved by routing
BOTH through the single shared builder `wallchamber.wall_chamber_structure` (the
`tau_tilting_block` precedent — each runner only adds `citations` locally). Anchors below
were **re-grepped against the live tree in the fix round** (2026-08-08); still "adjust to
reality" if a later merge shifts them.

**`wall_chamber` is an ALGEBRA-level compute kind** carrying a PAIR BUDGET (not a degree
range), parsed exactly like `tau_tilting`: `wall_chamber` or `wall_chamber:512`. NOT a
`MODULE_KINDS` entry; sized on `A.dim` like `tau_tilting`.

**Files (the seven-touchpoint checklist):**
1. **Parse grammar** — `src/quiverlab/hpc/spec.py::parse_compute_item` (the `tau_tilting`
   branch, `spec.py:267`): add `if s == "wall_chamber" or s.startswith("wall_chamber:"):`
   returning `ComputeItem(kind="wall_chamber", lo=None, hi=budget)`. Mirror in
   `docs/gui/runner.py` parse (`runner.py:165`).
2. **Server/HPC dispatch** — `spec.py::_dispatch` (after the `tau_tilting` branch,
   `spec.py:1549`): `if kind == "wall_chamber": budget = item.hi or 512; from
   quiverlab.tautilting.wallchamber import wall_chamber_structure; block =
   wall_chamber_structure(A, budget=budget); block["kind"] = "wall_chamber";
   block["citations"] = _citation_pairs(block["references"]); return block, None`. Catch the
   char-caveat `QuiverlabError` into `{"error": "<loud message>"}` (the Plan-30
   honest-per-entry precedent — never a 500). Add the `_snippet` calls entry
   (`spec.py:2552`, the `_snippet` calls dict): `"wall_chamber": lambda it: ("A.wall_chamber_structure(budget_pairs="
   f"{it.hi if it.hi is not None else 512})")`.
3. **Pyodide twin** — `docs/gui/runner.py::compute_one` (after the `tau_tilting` branch,
   `runner.py:1078`): the byte-identical `elif name == "wall_chamber":` branch calling
   `wall_chamber_structure`; one entry in the `calls` reproduce map (`runner.py:1333`):
   `"wall_chamber": "A.wall_chamber_structure(budget_pairs=%d)"`; one `ETA_MODEL["scalars"]`
   cost (`runner.py:1458`) — `"wall_chamber": 2.5` (the exchange-graph BFS + the per-brick
   submodule enumeration; size it just above `tau_tilting`'s 2.0).
4. **GUI JS** — `docs/gui/gui.js` AND `webapp/static/gui/gui.js` (diff-identical — apply
   every edit to both): a checkbox `<input id="qlgui-wall_chamber">` + budget picker
   `qlgui-wall_chamber-budget` in the structure block (beside `tau_tilting`, `gui.js:133`);
   the two bare ids in the id-registry array (`gui.js:257`); the request push-list entry
   (`gui.js:892`) `if (el.wall_chamber.checked) compute.push("wall_chamber:" +
   el["wall_chamber-budget"].value)`; a `renderBlock` `else if (name === "wall_chamber")`
   branch (`gui.js:3370`) calling a new `renderWallChamber(div, b)`; the `THEMES` structure
   `kinds` array entry `"wall_chamber"` inside the `QLGUI-THEMES-BEGIN/END` sentinels
   (`gui.js:3826`); the `scheduleProbe`/`KIND_CTRL` entries (`KIND_CTRL` at `gui.js:4344`).
   **`renderWallChamber`**: the
   chambers table (id | pair | support | g-matrix, mirroring `renderTauTilting`), the WALLS
   table (brick | dim-vector | full-line-or-ray | inequality system), the counts row, and —
   for `render ∈ {fan2d,fan3d}` — a new **`renderWallChamberSVG(div, b)`** that draws the
   chamber g-cones (colored, like `renderWallAndChamber`, `gui.js:2887`) AND **overlays each
   brick-wall** `D(B)` as a distinct labeled ray/line (n=2: the `Wall.rays`; n=3: the
   L1/net2d facet rays), labelled by `brick_name`/dim-vector — the visual difference from
   Plan 45's fan. When incomplete, render the bounded region + the honest `truncation` note.
5. **i18n (FOUR locales) + picker template** — `webapp/server/i18n/{en,es,fr,zh}.json` (the
   shipped `LANGS = ("en","es","fr","zh")`): add `pick.kind.wall_chamber` (REQUIRED — gated
   by `test_layout_picker.py`), `struct.wall_chamber` (checkbox label), and
   `wc.*` block strings (`wc.title`, `wc.chambers`, `wc.walls`, `wc.wall_line`,
   `wc.wall_ray`, `wc.inequalities`, `wc.counts`, `wc.truncated`, `wc.budget`). Add one
   `data-pick-kind-wall_chamber="{{ t('pick.kind.wall_chamber') }}"` line in
   `webapp/templates/draw.html` (mirror `data-pick-kind-tau_tilting`).
6. **Report renderer** — `src/quiverlab/trace/results_html.py`: `_HEADINGS["wall_chamber"]`
   (`results_html.py:76`) + a `_wall_chamber_html(b)` branch in `_block_html`
   (`_block_html` at `results_html.py:960`, the `tau_tilting` branch at `results_html.py:1055`)
   rendering the chambers table, the WALLS table (brick |
   full-line/ray | equality | inequalities), the counts (complete) or the truncation note
   (incomplete), and the fan via a new `viz/tikz.py::tikz_wall_chamber(b)` TikZ twin (mirror
   `tikz_fan`, `tikz.py:72` — draw chamber rays + the labeled brick-walls for n≤3; an honest
   "table only (n≥4 or budget-capped)" node otherwise).
7. **Tests / gates** — add `wall_chamber` to `ALL_KINDS` in
   `tests/webapp/test_layout_picker.py:26` (must equal the THEMES kind set, each once —
   `test_themes_cover_every_kind_exactly_once` + the four-locale `pick.kind.*` gate); add
   ONE golden `wall_chamber_kA2` to `tests/webapp/_runner_goldens.json` with a dated
   change-log bullet in `tests/webapp/test_runner_delegation.py` (verify existing entries
   byte-identical FIRST; on rebase merge the JSON semantically, never textually); add a
   twin-parity test (copy `test_tau_tilting_p45.py::test_twin_parity`).
- Test files: `tests/webapp/test_wall_chamber_p63.py`, `tests/gui/test_wall_chamber_runner_twin.py`.

```python
# tests/webapp/test_wall_chamber_p63.py
"""The wall_chamber algebra-level compute kind: served by hpc.spec, mirrored by the Pyodide
twin, both runners byte-identical. kA2 -> 5 chambers / 3 walls, complete, render fan2d."""
import json


def test_wall_chamber_block_shape(tmp_path):
    from quiverlab.hpc.spec import ComputeRequest, run
    req = _wall_chamber_request(quiver=("kA2",), budget=512)
    out = run(ComputeRequest.model_validate(req), tmp_path)
    b = out["results"]["wall_chamber"]
    assert b["complete"] and b["num_chambers"] == 5 and b["num_walls"] == 3
    assert b["render"] == "fan2d"
    assert "brustle_smith_treffinger" in b["references"]


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py::compute_one and hpc.spec.run;
    # assert json.dumps(block, sort_keys=True) equality on the wall_chamber block.
    ...


def test_wild_budget_status(tmp_path):
    # 2-Kronecker with a small budget -> complete False, status "budget", truncation set,
    # counts null, no crash.
    ...
```

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir).
- [ ] **Step 2: Implement** the parse grammar (both runners), the server + twin dispatch, the
  GUI JS in both `gui.js` files, the four-locale i18n + `draw.html`, and the
  `results_html.py` + `viz/tikz.py` branches.
- [ ] **Step 3: Gate the layout/i18n first** —
  `... -m pytest tests/webapp/test_layout_picker.py tests/webapp/test_js_parses.py -q`
  (THEMES == `ALL_KINDS`; `pick.kind.wall_chamber` in all four locales).
- [ ] **Step 4: Add the golden** (`wall_chamber_kA2`) to `_runner_goldens.json`; note it in
  the `test_runner_delegation.py` docstring change-log. Run the delegation test BEFORE
  adding to confirm existing entries stay byte-identical.
- [ ] **Step 5: Run the gates** —
  `... -m pytest tests/webapp/test_wall_chamber_p63.py tests/webapp/test_runner_delegation.py tests/gui/test_wall_chamber_runner_twin.py tests/hpc -q`.
  Expected PASS; both runners byte-identical; the wild case returns `status="budget"` cleanly.
- [ ] **Step 6: Commit** — `feat(gui,webapp,hpc,trace): wall_chamber compute kind + LIVE brick-wall overlay SVG (n<=3) + tikz twin + one golden (4 locales, v1 schema)`

**Adjust to reality (Task 5):**
- The block already carries `"kind": "wall_chamber"` if `wall_chamber_structure` sets it;
  otherwise set it in `_dispatch` (the `tau_tilting_block` sets its own `"kind"`, mirror
  that — set it inside the builder so both runners agree).
- `test_layout_picker.py::ALL_KINDS` is a frozen set the THEMES array must equal exactly;
  add `wall_chamber` to BOTH or the gate fails. Confirm the THEMES sentinel block in BOTH
  `gui.js` files (vendored byte-identical) gets the new kind.
- If `estimator.sizing_dim` (`webapp/server/estimator.py`) needs an entry, size
  `wall_chamber` on `A.dim` exactly as `tau_tilting` is sized (algebra-level; not a module
  kind) — check the `tau_tilting` path and mirror it.

---

### Task 6: verification page, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick the P63 card)
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`)

- [ ] **Step 1: Verification page.** Add the Plan-63 subsystem row (`tautilting/wallchamber.py`):
  - `oracle_literature` — kA₂ = 5 chambers / 3 walls (the P₁ wall a ray); **kZ₂/rad² = 6
    chambers / 4 walls** (self-injective, non-thin — P₁, P₂ opposite half-rays of one
    hyperplane, BST Rem. 3.19); the DIJ brick-finite ⟺ τ-tilting-finite gate; the
    Kaipel–Treffinger worked examples (Ex. 13 kA₂; Ex. 15 the **CYCLIC** rad²-Nakayama N₃²,
    14 chambers — engine-verified against the exact algebra; provenance: an earlier draft
    mis-transcribed it as the LINEAR `1→2→3`/rad², which is 12, not 14).
  - `oracle_crossengine` — `#chambers = #support τ-tilting` (fan vs exchange graph;
    hereditary kA₂=5, kA₃=14); **linear kA₃/rad² = 12** (NOT KT, NOT Catalan); the
    inequality-defined `D(B)` rays ≡ the grouped exchange-edge facets **as sets of primitive
    directions** (two independent constructions — genuine set equality, not just non-empty).
  - `oracle_selfcert` — `is_full_hyperplane`/codim; `θ·dim B = 0` on every wall facet; the
    facets ↔ exchange-edges bijection (each facet a real edge between two chambers); every
    wall ray lies in `D(B)` (equality + all submodule inequalities); **counts consistency
    `#walls == #bricks == len(torsion.bricks(A))`** on a thin and a non-thin example; each
    MGS is a monotone chamber path.
  - the **honest semi-decision entry** (metaplan §6 ledger): certified complete iff
    brick-finite ⟺ τ-tilting-finite; on the 2-Kronecker (brick-infinite) the structure is a
    **bounded region** with `status="budget"`, no count claimed.
  - the **honest-scope entries** (inherited from Plan 45 + new): (a) rigorous only over
    char 0 / char > dim (QQ default; loud off scope); (b) bricks decide over the
    algebraically-closed / char-0 base (`end_dim=1`; the GF(pⁿ) proper-division-ring
    caveat); (c) the drawing is rank ≤ 3 only (n=3 via the L1/octahedron projection; n=3
    wall rays are the grouped exchange-edge facets, self-cert `θ ∈ D(B)`); rank ≥ 4 = tables;
    (d) **QPA cannot compare** — the honest no-surface probe (FD-Applet/DIJ tables are the
    named external cross-checks, not a live oracle); (e) **the Kaipel–Treffinger worked-example
    counts are engine-verified against the EXACT algebra** (kA₂ = 5; cyclic rad²-Nakayama N₃²
    Ex. 15 = 14) — no `# PIN` remains deferred; the record notes the corrected provenance
    (the earlier linear-vs-cyclic mis-transcription); (f) **self-injective input is in scope**
    (bricks decide, no hereditary assumption — the cyclic N₃² and kZ₂/rad² tests cover it).
  Recount the class table (`tests/release/test_oracle_classes.py` drives the numbers —
  collect, paste the LIVE counts, re-run to green; mid-merge-train honest).
- [ ] **Step 2: README.** One features line: "wall-and-chamber structure via bricks
  (Brüstle–Smith–Treffinger): the wall `D(B)` of every brick as an exact inequality system,
  chambers = g-vector cones, certified complete iff τ-tilting-finite — with a LIVE 2D fan
  drawing for rank ≤ 3, clickable via `wall_chamber`."
- [ ] **Step 3: Full gate:**
  `... -m pytest tests/modules/test_wall_chamber*.py -q` (deep, the touched files),
  `... -m pytest -q -m fast`, `... -m pytest tests/webapp tests/gui -q`,
  `... -m pytest tests/qpa -q -m qpa`, `... -m pytest tests/release tests/citations -q` —
  all green.
- [ ] **Step 4: Commit** — `docs(verification): Plan-63 wall-and-chamber oracle rows + BST/Asai/KT citations + honest scope (char caveat, bricks base, rank<=3 drawing, brick-finite gate, no-QPA) + recounted classes`

---

## Acceptance (Plan-63 definition of done)

1. `wall_of_brick(B)` returns `Wall` — the EXACT inequality system `D(B) = {θ·dim B = 0,
   θ·dim N ≤ 0 ∀ N ⊆ B}` over the submodule dim-vectors, with `is_full_hyperplane`, `codim`,
   and (n ≤ 2) the exact extreme rays; `wall_chamber_structure(A)` and
   `Algebra.wall_chamber_structure(budget_pairs)` return the full payload (chambers =
   g-cones, walls = one per brick, adjacency, counts, render, geometry) — all exact rational,
   float-free.
2. **kA₂ = 5 chambers / 3 walls** pinned, with `D(S₁)`/`D(S₂)` full lines and **`D(P₁)` a
   single ray** (the record's headline) — derived by hand in the Reference section and
   re-derived by the engine.
3. `#chambers = #support τ-tilting` (HEREDITARY kA_n Catalan(n+1): kA₂=5, kA₃=14) pinned
   cross-engine; linear kA₃/rad²=12 and cyclic rad²-Nakayama N₃²=14 engine-pinned separately;
   the inequality-defined `D(B)` rays ≡ the grouped exchange-edge facets **as primitive
   direction sets** (two independent constructions); the facets ↔ edges bijection (each facet
   shared by exactly two chambers); every wall coplanar (`θ·dim B = 0`); the **non-thin
   discriminator kZ₂/rad² = 4 walls** (P₁, P₂ not collapsed); `#walls == #bricks`.
4. The structure is **certified complete iff brick-finite ⟺ τ-tilting-finite** (DIJ); on the
   2-Kronecker it returns a **bounded region** with `complete=False`, `status="budget"`, a
   `truncation` note, and **no count** — the honest truncation, never a claimed-complete fan.
5. **MGS = green paths**: each maximal green sequence verified as a monotone source→sink
   chamber path (`_is_green_path`); kA₂ = 2 pinned.
6. `wall_chamber` clickable end-to-end (GUI canvas → block → **live brick-wall overlay SVG**
   for rank ≤ 3, table otherwise → report TikZ) in all four locales (en/es/fr/zh), both
   runners byte-identical via the shared builder, ONE golden added with a documented
   change-log entry, `test_layout_picker.py` `ALL_KINDS`/THEMES gate green, schema still v1,
   canonical keys unchanged.
7. QPA honest-scope probe skips (fails if QPA ever ships a wall-and-chamber surface);
   `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the honest
   scope; README line added; deep (wall-chamber files) + fast + webapp/gui + qpa + release +
   citations suites green. Plan 45 left byte-unchanged. Honest scope recorded: char 0 / char
   > dim (QQ default); bricks over the algebraically-closed base; rank ≤ 3 drawing (L1 for
   n=3); QPA cannot compare; `# PIN`ned reference numbers deferred.

---

## Methodology & assumptions

**Approach.** I read Plan 45 (the C4 τ-tilting engine) end-to-end — it is the substrate,
and P63 builds directly on its shipped surface. I then read the shipped
`src/quiverlab/tautilting/{stability,torsion,mutation,block}.py` and confirmed the exact
APIs P63 consumes: `exchange_graph`/`ExchangeGraph` (the honest `is_complete`/`status`
contract), `wall_normal`, `stability._submodule_dimvecs` (the shipped exact submodule
enumerator), `is_theta_semistable` (King, submodule ≤ 0), `wall_and_chamber_fan`,
`torsion.bricks` (edge-brick enumeration, iso-class deduped), `green.maximal_green_sequences`,
`modules.hom.identify_standard`. I read Plan 62 end-to-end for the house style and the
honest bounded-search / certified-completeness discipline P63 must match, and I mapped the
compute-kind wiring by grepping `hpc/spec.py` (`parse_compute_item`, `_dispatch`, `_snippet`
calls map), `docs/gui/runner.py` (`compute_one`, `calls`, `ETA_MODEL`), `docs/gui/gui.js`
(`renderTauTilting`/`renderWallAndChamber`, THEMES sentinels, id registry), the four-locale
i18n (`LANGS = ("en","es","fr","zh")`), `trace/results_html.py`, `viz/tikz.py`, and
`tests/webapp/test_layout_picker.py`.

**Empirical ground truth (run this authoring).** I ran the merged Plan-45 engine on kA₂ in
read-only Python: `exchange_graph` gives 5 pairs with g-matrices matching the 5 quadrant/
half-quadrant cones; the 5 edge-facets carry brick labels `{S₁:2, S₂:2, P₁:1}`;
`_submodule_dimvecs(P₁) = {(0,0),(0,1),(1,1)}`, `_submodule_dimvecs(S₁)={(0,0),(1,0)}`,
`_submodule_dimvecs(S₂)={(0,0),(0,1)}`. This gives, by hand, `D(P₁)={θ₁+θ₂=0, θ₂≤0}` (a
RAY), `D(S₁)={θ₁=0}`, `D(S₂)={θ₂=0}` (full lines) — 3 brick-walls, 5 chambers. The
literature agent then confirmed this **character-for-character** against Kaipel–Treffinger
Example 13 (verbatim wall equations) and BST Def. 3.1–3.3 / Remark 3.19.

**Fix-round re-verification (2026-08-08).** The rank-3 pins were re-run to completion on their
EXACT algebras: linear kA₃/rad² (`1→2→3`, `a*b=0`) = **12** pairs / 5 bricks; cyclic
rad²-Nakayama N₃² (`1→2→3→1`, all length-2 = 0) = **14** pairs / 6 bricks (n_regular);
hereditary kA₃ = **14**. The discriminating self-injective kZ₂/rad² (`1⇄2`, both length-2 = 0)
= **6** pairs / **4** bricks, its two dim-(1,1) bricks carrying `_submodule_dimvecs` =
`{(0,0),(1,0),(1,1)}` and `{(0,0),(0,1),(1,1)}` ⇒ D(P₁), D(P₂) the opposite half-rays
`(-1,1)` and `(1,-1)` of `θ₁+θ₂=0`. These correct the original draft's rank-3 pin, which had
mis-transcribed KT Ex. 15 as the linear algebra (12, not 14).

**Key correctness decisions.** (1) **The wall `D(B)` is a polyhedral CONE, not a hyperplane
normal** — the genuine new object over Plan 45, whose fan records only per-edge normals
(5 edges for kA₂) and stops the geometry at n=3. P63 computes `D(B)` from the submodule
inequalities (BST/KT Def), grouping exchange-edge facets into 3 brick-walls (BST Remark 3.19
verbatim: "one wall `D(N)` can be made of more than one facet"). (2) **The King convention
is reused, not re-chosen** — Plan 45's `is_theta_semistable` uses (submodule, ≤ 0), which is
EXACTLY BST Def. 3.1 / KT Def. 3; reusing it makes `D(B)` byte-consistent with the shipped,
tested predicate and matches the papers. (3) **The completeness gate is DIJ brick-finite ⟺
τ-tilting-finite**, decided by the Plan-45 BFS; the brick-infinite case returns a **bounded
region with honest truncation** (the record's mandate and the P62 discipline) — unlike Plan
45's fan, which returns empty. Each discovered chamber/wall is exact; only totality and the
counts are withheld. (4) **Plan 45 is left byte-unchanged** — P63 is a new `wallchamber.py`
+ one algebra-level compute kind, calling `exchange_graph` directly so it can render the
bounded region the fan cannot. (5) **Attribution is corrected per the agent:** walls=D(brick)
→ DIJ Thm 1.4/1.5 + BST Ex. 3.4 (not one BST theorem); Asai *Semibricks* is cited only for
the categorical INDEXING (it contains no fan/g-vector content); "complete fan ⟺
τ-tilting-finite" is DIJ Thm 1.7 + BST Cor. 3.29.

**Deliberately not checked / assumptions.** The deep rank-3 counts, which timed out under the
120s cap during the ORIGINAL authoring, were **run to completion in the fix round** (fresh
venv probes, recorded in the change log): linear kA₃/rad² = **12** pairs (5 bricks), cyclic
rad²-Nakayama N₃² = **14** pairs (6 bricks, n_regular), hereditary kA₃ = **14**, kZ₂/rad² =
**6** pairs (4 bricks with the two (1,1) submodule sets). So Ex. 15 no longer "rests on a
Plan-45 pin" — it is engine-verified against the exact (cyclic) algebra, and the earlier
linear mis-transcription (12, not 14) is corrected. I did not
compute the direct 2D extreme rays of an `n=3` wall cone from the inequality system (a
Fourier–Motzkin on 2 free coordinates) — the plan uses the grouped exchange-edge facets +
an `θ ∈ D(B)` self-cert for `n=3`, which is sound for the complete fan; the direct solver is
an optional refinement (recorded as honest scope). I did not verify `estimator.sizing_dim`
needs a `wall_chamber` entry (it may size algebra-level kinds generically like `tau_tilting`)
— flagged as an adjust-to-reality in Task 5. I did not byte-extract the BST literal-prose "5
chambers" (it is read off Figure 1 + the sτ-tilt count; the wall EQUATIONS are verbatim from
KT). The GF(pⁿ) proper-division-ring brick caveat and the char 0 / char > dim scope are
inherited from Plan 45 unchanged (batteries over QQ).

**Why this is correct.** The wall definition, the sign convention, the chamber↔τ-tilting
bijection, the brick-finite gate, and both worked examples (kA₂ 5/3; the CYCLIC rad²-Nakayama
N₃² Ex. 15 = 14) are verbatim-transcribed from BST, KT, and DIJ (double-sourced by the
reference agent) and independently reproduced by the merged Plan-45 engine on the EXACT
algebras in the fix round (kA₂ 5/3; cyclic N₃² 14; and the discriminating kZ₂/rad² 6/4). Every geometric object is
exact rational (reusing Plan 45's `Fraction` machinery); the completeness claim is gated on
the decidable DIJ criterion and refuses honestly (bounded region, no count) otherwise; and
the two independent constructions of each wall — the inequality-defined `D(B)` and the
grouped exchange-edge facets — cross-validate as an `oracle_crossengine` pin. `#chambers =
#support τ-tilting` ties the geometry back to the Plan-45 exchange graph from a different
route.

## Change log

- **2026-08-07 authoring** — plan drafted against merged `dev` (Plan 45 + Wave-1
  P51–P53+P55+P57+P58); scope boundary vs Plan 45 fixed (D(B) as a cone, all-rank
  inequalities, brick-grouped walls "5/3", honest bounded region); King convention reused
  from the shipped `is_theta_semistable`.
- **2026-08-07 reference re-verification folded in** — the literature agent's verbatim,
  double-sourced findings applied: BST Def. 3.1–3.3 / Cor. 3.18 / Rem. 3.19 transcribed; the
  attribution of "walls = D(brick)" corrected to DIJ Thm 1.4/1.5 + BST Ex. 3.4 (not a single
  BST theorem); Asai *Semibricks* scoped to the categorical INDEXING only (no fan content —
  citation annotation fixed); Kaipel–Treffinger title + Examples 13 (kA₂ wall equations,
  D(P₁) a ray) and 15 (rad²-Nakayama, 14 chambers) promoted from `# PIN` placeholders to
  firm `oracle_literature` pins; the "complete fan ⟺ τ-tilting-finite" packaging documented
  as DIJ Thm 1.7 + BST Cor. 3.29. **(SUPERSEDED: Ex. 15 was promoted against the LINEAR
  `1→2→3`/rad², which the fix round found is 12 chambers, not 14 — the correct KT Ex. 15
  algebra is the CYCLIC N₃²; see the 2026-08-08 entry.)**
- **2026-08-08 critic fix round (NEEDS WORK → adjudicated all-valid)** — eight rulings
  applied, DOCUMENT-ONLY:
  - **H1/H2 (blocking):** the rank-3 KT Ex. 15 pin was against the WRONG algebra — the draft
    built the LINEAR `1→2→3`/rad² (engine: **12** support τ-tilting pairs) and claimed 14.
    KT Ex. 15's module list contains `3/1`, forcing an arrow `3→1`, so KT's algebra is the
    **CYCLIC** N₃² (`1→2→3→1`/rad², self-injective) — engine-verified **14** pairs / 6 bricks.
    Swapped the test algebra + all prose to the cyclic one; deleted the false "linear = 14 =
    Catalan(4)" and the fabricated `D(1/2)={(x,−x,z):x≥0}` wall; added a separate
    engine-verified `oracle_crossengine` pin (linear kA₃/rad² = 12, NOT KT). Restored truthful
    provenance (engine-verified against the exact algebra + KT module-list ID; mis-transcription
    named; reference re-verification rule satisfied by running the engine at spec time).
  - **M1 (major):** added kZ₂/rad² as the PRIMARY discriminating oracle (self-injective, 6
    chambers / 4 walls; P₁, P₂ both (1,1) on opposite half-rays (−1,1)/(1,−1) of θ₁+θ₂=0, BST
    Rem. 3.19); fixed the Task-2 grouping pseudocode to group by ISO-CLASS (via
    `torsion._edge_brick` + `is_isomorphic`), never by dim-vector — named the trap explicitly
    (a dim-vector key merges P₁, P₂ → 3 walls; correct is 4). The shipped `torsion.bricks` /
    `_edge_brick` (`torsion.py:125–174`) already do this.
  - **W1 (major):** `test_n2_wall_rays_equal_grouped_exchange_facets` now primitivises both
    ray sets and asserts SET EQUALITY of directions per brick (was vacuously "both nonempty").
  - **M2 (minor):** added `test_counts_consistency_walls_equal_bricks` (canonical brick set =
    `torsion.bricks(A)`; asserts `#walls == #bricks == len(torsion.bricks(A))` on a thin and a
    non-thin example).
  - **M3 (minor):** stated self-injective input is in scope (global constraint); the cyclic
    N₃² + kZ₂/rad² tests double as that coverage.
  - **W2 (minor):** de-vacuumed `test_n3_walls_project_via_L1` (now asserts every wall ray
    lies in D(B): equality + all submodule inequalities) and
    `test_each_facet_bounds_exactly_two_chambers` (now the facets ↔ exchange-edges bijection).
  - **H3 (minor):** refreshed the Task-5 wiring anchors against the live tree
    (`runner.py` 1062→1078 / 1293→1333 / 1414→1458; `gui.js` 126→133 / 246→257 / 881→892 /
    3267→3370 / 3722→3826 / 2876→2887; `results_html.py` 72→76, `_block_html`@960 branch@1055;
    `spec.py` `_snippet`@2552; `gui.js` `KIND_CTRL`@4344); hedge kept.
  - Live-verified counts (fresh venv probes): linear kA₃/rad² = 12, cyclic N₃² = 14
    (n_regular), hereditary kA₃ = 14, hereditary kA₂ = 5, kZ₂/rad² = 6 (4 bricks, the two
    submodule sets confirmed). `gui.js` docs↔webapp confirmed byte-identical.
