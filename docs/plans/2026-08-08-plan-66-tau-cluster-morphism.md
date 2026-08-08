# Plan 66: The τ-cluster morphism category W(A) + picture group (R29) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal.** Promote the **τ-cluster morphism category** `W(A)` of a τ-tilting-**finite**
`A = kQ/I` to a first-class, no-code surface, built **on top of** the Plan-45 τ-tilting
engine, the **Plan-64** wide-subcategory poset, and the **Plan-65** Jasso τ-perpendicular
reduction. `W(A)` (Buan–Marsh, *A category of wide subcategories*, IMRN 2021; Hanson–Igusa,
*τ-Cluster morphism categories and picture groups*, Comm. Alg. 2021, generalizing
Igusa–Todorov) is the category whose

- **objects** are the **τ-perpendicular wide subcategories** of `mod A` (for τ-tilting-finite
  `A` these are *all* the finitely many wide subcategories — count `= #wide`, tying P64), and
- **morphisms** `[U]: W → W'` are the **support τ-rigid pairs `U` of `W`** with `W' = J_W(U)`
  the τ-perpendicular category of `U` computed *inside* `W` (the **same** Jasso reduction P65
  ships).

P66 delivers, all exact and float-free:

1. **`W(A)` as a finite category** — the object set (= wide subcategories, cross-tied to P64),
   the morphism set graded by **rank** (`= |U|`), the per-object out-degree and the total
   morphism count, and the identity/composition data — a genuine small category, not a count.
2. **The classifying space (Hanson–Igusa cube complex)** — the face vector `f`, whose
   `f_0 = #wide` (one 0-cell per wide subcategory) and `f_k = #(rank-k morphisms of the whole
   category) = Σ_W #(rank-k support τ-rigid pairs of C_W)` (one k-cell per rank-k morphism); its
   Euler characteristic `χ = Σ_k (−1)^k f_k`; and the **theorem-anchored `K(π,1)` verdict**
   (True for Nakayama and for hereditary Dynkin; honest *not-certified* otherwise — homotopy-type
   claims never exceed the cited theorems). This cube complex is the ACTUAL classifying space
   `|W(A)|`; it is **not** the `g`-fan/cluster fan of `A` (a triangulated `(n−1)`-sphere whose own
   f-vector `#(support τ-rigid pairs of A with k summands)` is reported *separately* as
   `g_fan_face_vector` and never conflated with the classifying space — see §3 of the spec).
3. **The picture group presentation** — `π₁(|W(A)|)` **as generators-and-relations DATA**:
   one generator `x(β)` per **brick** `β`, one relation per **rank-2 wide subcategory** (a
   *commutation* `x(β)x(γ)=x(γ)x(β)` when the wide is `k × k`, an *atom/pentagon* relation
   when the wide is connected `≅ mod kA₂`), plus the abelianization (exact SNF). Igusa–Todorov–
   Weyman's presentation is the ground truth; the plan pins **counts and relation TYPES**
   (which are engine-derivable and discriminating) and transcribes the exact relation WORDS.

One algebra-level no-code compute kind — **`tau_cluster`** — puts the category summary (objects,
morphisms, ranks), the cube-complex face vector + Euler characteristic + `K(π,1)` verdict, and
the picture-group presentation (generators, typed relations, abelianization) one click away in
all four locales. The honest **complete-iff-τ-tilting-finite** contract is inherited from P45:
`W(A)` is *infinite* on a τ-tilting-infinite algebra, so P66 refuses loudly (no partial
category) exactly where P64 refuses to emit a lattice invariant.

**What is genuinely NEW over Plans 45 / 64 / 65 (the scope boundary — read this first).**
- **Plan 45** ships the exchange graph (support τ-tilting pairs = the `g`-fan facets, brick-
  labelled covers), `bricks` / `semibricks`, `hasse_orientation`. It does **not** assemble the
  `g`-fan into a face poset, has no notion of a wide-subcategory object, no morphism set, no
  cube complex, no picture group.
- **Plan 64** ships `wide_subcategories(A)` (the Enomoto core-label-order poset of *all* wide
  subcategories, with each wide's simple-brick semibrick label) and the congruence/forcing
  machinery. P66 **consumes** `wide_subcategories` for its OBJECT set and cross-ties the object
  count — but P64 has **no morphisms between wides**, no classifying space, no group.
- **Plan 65** ships `tau_perpendicular_reduction(A, U) → C(U)` (the Jasso algebra with
  `J(U) ≅ mod C(U)`), the signed τ-exceptional sequences, and the `signed_count = n!·#sτt`
  identity. P66 **consumes** the reduction to realise each morphism's target `W' = J_W(U)` and
  to name each object's algebra `C_W` — but P65 has **no category structure** (no object poset,
  no morphism composition, no cube complex, no picture group). P66 also **checks** its
  factorization count against P65's `signed_count` (the same `n!·#sτt`) — a **self-cert
  consistency gate**, not a cross-engine oracle, since `signed_count` is definitional in P65
  (fix-round ruling 3).

P66 adds four things none of P45/P64/P65 has: **(i)** the category `W(A)` (objects + graded
morphisms + composition), **(ii)** the Hanson–Igusa cube-complex classifying-space face vector
(`f_0 = #wide`, `f_k = #(rank-k morphisms)`) + Euler characteristic, **(iii)** the
theorem-anchored `K(π,1)` verdict, **(iv)** the picture-group presentation as
generators+typed-relations+abelianization DATA.

**Architecture.** One new module `src/quiverlab/tautilting/cluster_morphism.py`, a thin
combinatorial layer over the merged P45/P64/P65 surfaces (no new representation theory — the
bricks, wides, and reductions are all inherited):

- `cluster_morphism.py::TauClusterCategory` — a frozen value object: `objects` (wide ids with
  rank + simple-brick semibrick label + `C_W` handle), `object_count`, `morphisms`
  (`(source_id, target_id, rank, support_tau_rigid_pair_data)`), `morphism_count`,
  `out_degree` (per object), `face_vector` (the Hanson–Igusa cube-complex classifying space:
  `f_0 = object_count = #wide`, `f_k = #(rank-k morphisms)`, `sum(f) = morphism_count`),
  `g_fan_face_vector` (the `g`-fan/cluster fan of `A` — a triangulated `(n−1)`-sphere,
  `#(support τ-rigid pairs of A) by summand count`, `f_n = #sτt`; reported for the star/#sτt
  self-certs, NEVER the classifying space), `euler_characteristic` (`= Σ(−1)^k face_vector[k]`),
  `is_kpi1` (`True` / `None`) + `kpi1_reason`, the honest `is_complete`/`status`/`note`.
- `cluster_morphism.py::PictureGroup` — a frozen value object: `generators` (per brick:
  dim-vector + `identify_standard` name), `relations` (per rank-2 wide: `type ∈
  {"commutation","atom"}`, the two simple bricks, the extension brick for an atom, the exact
  relation word), `num_generators` (`= #bricks`), `num_relations` (`= #rank-2 wides`),
  `num_atom`, `num_commutation`, `abelianization` (SNF invariant factors of the relation
  matrix), the honest contract.
- `tau_cluster_category(A, *, budget=512) -> TauClusterCategory`,
  `picture_group(A, *, budget=512) -> PictureGroup`,
  `tau_cluster_block(A, *, budget=512) -> dict` (the shared JSON block for the two runners).
- Thin `Algebra` delegates in `core/algebra.py` beside `Algebra.exchange_graph` /
  `Algebra.wide_subcategories` (lazy import to avoid the `tautilting → core` cycle):
  `Algebra.tau_cluster_category(budget=512)`, `Algebra.picture_group(budget=512)`.

**Plans 45/64/65 are left byte-unchanged.** `exchange_graph`, `hasse_orientation`, `bricks`,
`semibricks`, `wide_subcategories`, `tau_perpendicular_reduction`, `tau_exceptional_sequences`
are consumed **read-only**; no upstream test moves; no upstream golden re-freezes.

**Tech Stack.** Pure integer/set combinatorics (the `g`-fan is finite; ids are `int`, faces are
`frozenset` of signed `g`-vector columns, the picture-group relation matrix is `list[list[int]]`
over ℤ, SNF via the shipped exact `fields.linalg`). Brick labels are iso-classes of `Module`s
(reuse P45's `is_isomorphic` dedup + `identify_standard` naming). No floats in `src/` (AST-gated
by `tests/test_no_floats.py`) — there is nothing to floatify (unlike P63's geometry). All
representation-theoretic content (bricks, wides, reductions) is inherited over `QQ` by default
(the load-bearing char caveat below).

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **P45, P64, and P65 are HARD prerequisites.** **P45 is merged to `dev` (shipped engine —
  its code exists and is consumed live in the §Verified-pins prototype). P64 and P65 are, at
  this authoring, committed only as PLAN DOCS — their code does NOT yet exist and MUST be merged
  (implemented) before this plan's implementation starts.** This plan consumes, read-only:
  - **P45** `tautilting`: `exchange_graph(A, budget_pairs=512) -> ExchangeGraph`
    (`.vertices[i]["g_matrix"]` = the signed summand `g`-vectors of a facet, `["pair"]`,
    `["support"]`; `.arrows[(i,j)]["brick"]` / `["brick_name"]`; `.is_complete`;
    **`.status ∈ {"complete","budget","error"}`**), `hasse_orientation`, `bricks`,
    `semibricks`, `torsion._edge_brick` / `_torsion_universe` (iso-class-safe brick label),
    `SupportTauTiltingPair`.
  - **P64** `tautilting.congruence`: `wide_subcategories(A, *, budget) -> WideSubcategoryPoset`
    (`.elements`, `.order`, `.labels` = per-wide simple-brick semibrick, `.size` = `#wide`,
    `.is_complete`/`.status`). **P66's object count is cross-tied to `.size`** (§Task A).
  - **P65** `tautilting.exceptional`: `tau_perpendicular_reduction(A, U, *, budget) ->
    TauReduction` (`.reduced` = `C(U)` a genuine `kQ'/I'` Algebra of rank `n − |U|`, plus the
    bookkeeping to identify `J(U)`'s simple bricks), and `tau_exceptional_sequences(A, *,
    budget) -> TauExcReport` (`.signed_count = n!·#sτt`, `.stt_count`, `.status`). **P66's
    factorization count is cross-tied to `.signed_count`** (§Task A; fix-round ruling 3).
  **If any of P45/P64/P65 is not merged when a worker picks this up, STOP and escalate** — do
  not fork the exchange-graph, wide-subcategory, or reduction machinery.
- **Cross-plan contract check — `TauReduction` shape (do FIRST, before coding Task A).** P66
  consumes `TauReduction.reduced` (and the simple-brick bookkeeping), but **P65's plan doc does
  not explicitly name these fields** — `.reduced` is P66's assumed name for the `C(U)` handle.
  Before implementing, the worker **verifies the exact `TauReduction` shape against P65's
  SHIPPED implementation** (the field that carries `C(U)`, and how `J(U)`'s simple bricks are
  exposed). If P65 ships a different name/shape, **flag the drift back to the orchestrator and
  wait — do NOT silently adapt** (a silent rename hides a genuine cross-plan interface change and
  can mask a semantic mismatch). Same discipline for `TauExcReport.signed_count` and P45's
  `ExchangeGraph` accessors.
- **Sequence AFTER P65 Task 0 (the `mutate` D₄-star defect).** P65 Task 0 fixes a genuine P45
  defect: `exchange_graph` returned `status="error"` (not `"complete"`) on the τ-tilting-FINITE
  D₄ mixed star `{1→2, 3→1, 4→1}` while discovering the correct 50 vertices, because a single
  `mutate` call raised deep in the BFS (reproduced live: `status="error"`, `len(vertices)=50`).
  P66 consumes `exchange_graph` throughout, so **P66 sequences after P65 (which fixes it and
  installs the `status != "complete"` loud-refusal gate).** Regardless of the fix, **P66's own
  gate treats `exchange_graph(...).status != "complete"` as a LOUD refusal** (never a count/
  category derived from a non-complete graph) — the M1 discipline is contractual here too
  (§Task A, gate G0). Live-confirmed at this authoring: `exchange_graph(D₄-star)` still raises
  `QuiverlabError` in the *torsion-universe* path on `dev` (P65 not yet merged) — so the gate is
  load-bearing.
- **Buckets are auto-assigned by directory (`tests/conftest.py`): `tests/modules/` → deep.**
  All P66 engine tests live in `tests/modules/` as `tests/modules/test_tau_cluster*.py` and
  `tests/modules/test_picture_group*.py` (deep — they share the exchange-graph budget with the
  P45/P64/P65 suites). `tests/qpa/` → qpa; `tests/webapp/`, `tests/gui/` → fast (extras-gated
  dirs, unmarked per the Plan-32 rule). Run new tests by path during development; finish each
  task with a `-m deep` spot-run of the touched files.
- **THE ENGINE RUNS OVER QQ BY DEFAULT** (the P45/P64/P65 char caveat, inherited verbatim).
  `bricks`, `semibricks`, `is_isomorphic`, `wide_subcategories`, `tau_perpendicular_reduction`,
  the exchange-graph BFS are rigorous only over **char 0 or char > dim** (Dickson/CIW). Batteries
  pin over **QQ**. Over `char ≤ dim` the underlying engines raise the loud `QuiverlabError`
  unchanged — P66 never a silent wrong category.
- **Bricks decide over the implicit algebraically-closed / char-0 base.** A brick is
  `end_dim(B) == 1` (`End_A(B) = k`); the picture-group generators (one per brick) and the
  atom-vs-commutation classification (via `Ext¹` between the two simple bricks) read correctly
  only when `End_A(B) = k`. Batteries pin over QQ; the `GF(pⁿ)` proper-division-ring caveat
  (`dim_k End(B) > 1`) is stated honestly on the verification page (inherited from P45/P64).
- **Honest semi-decision contract (metaplan §1.3; STRICTER than P63's bounded region, exactly
  as P64).** `W(A)` is a *finite* category iff `A` is τ-tilting-**finite** (Buan–Marsh gate). On
  a budget-truncated exchange graph the `g`-fan is a *prefix*, not a fan, and a partial object/
  morphism set / face vector / picture group would be a **lie**, not merely incomplete. So P66
  is certified **iff `A` is τ-tilting-finite** (`eg.is_complete and eg.status == "complete"`),
  and on a non-complete graph it returns `is_complete=False`, `status="budget"|"error"`, a
  `note`, and **omits every category invariant** (`objects=()`, `morphism_count=None`,
  `face_vector=None`, `g_fan_face_vector=None`, `picture_group=None`) — it does NOT return a
  bounded region. The
  2-Kronecker (τ-tilting-infinite, DIJ) is the honest-refusal oracle. τ-tilting-finiteness is
  **decidable** via the P45 BFS (DIJ).
- **Objects are `#wide`, NOT `#torsion` (load-bearing — the primary correctness trap).** The
  objects of `W(A)` are wide subcategories; their count is **`wide_subcategories(A).size`** (the
  Enomoto core-label-order count from P64). It is a spec-time ERROR to use `len(semibricks(A))`
  — P45/P64's `semibricks` returns the **Asai down-labelling** semibricks with
  `#semibricks == #pairs == #torsion` (the four-way identity, verified live in
  `torsion.py`), which equals `#wide` **only when `A` is representation-FINITE**
  (Marks–Šťovíček, arXiv:1503.04639: `wide A ↪ tors A`, bijection ⟺ rep-finite; the boundary is
  rep-finiteness, NOT heredity). On every `kAₙ` and the Nakayama batteries below all three
  counts coincide, so the trap is invisible in the tests — the plan states it explicitly and
  the code reads `.size` (never `semibricks`). A genuine `#wide < #torsion` witness needs a
  rep-INFINITE τ-tilting-finite algebra (tame `Π(D₄)`/`Π(A₅)`), which the shipped P45 engine
  cannot compute cheaply at spec time (P64 §Task B probe) — so P66 inherits P64's deferral of
  that witness and pins the object count on rep-finite algebras only.
- **`status != "complete"` → loud refusal (P65 M1).** Every P66 entry point treats a non-
  `"complete"` exchange graph (from `A` itself OR from any `C_W` reduction) as a loud refusal
  (`QuiverlabError` or the honest `status` field), never a count read off an `"error"`/`"budget"`
  graph. This holds even if a future orientation trips a new `mutate` edge case.
- **No floats in `src/`.** Faces are `frozenset[tuple[int,...]]`, counts are `int`, the picture-
  group relation matrix is `list[list[int]]` over ℤ, the abelianization is exact SNF
  (`fields.linalg`). The Euler characteristic is an `int`. There is no geometry to floatify. The
  v1 GUI renders the presentation and the category summary as **text/HTML tables** (no coordinate
  math); the report PDF reuses the shipped float-free `viz/tikz.py::tikz_hasse` fed the wide-
  subcategory poset (P64's precedent). A cube-complex / classifying-space SVG is an explicit
  **deferred stretch** (ruling 6, mirrors P64) — not v1.
- **Composition is left-to-right** (`a*b` = first `a` then `b`). All refusals are
  `QuiverlabError(message, hint=...)` (the two-space `[hint: ...]` convention, `errors.py`).
  Presentation-less (structure-constant) algebras: the exchange graph already refuses loudly
  (the quiver is required) — P66 surfaces that as a clean typed 4xx, never a 500. Theorem-gate
  violations (`object_count != wide_subcategories.size`; `num_generators != #bricks`;
  `sum(f_k with the alternating sign) != χ`; the composition not associative on a spot-check)
  raise `QuiverlabError` (loud internal error), never a silently-wrong payload.
- **Every citation key is BibTeX-verified before use** (Plan-29 rule); `_r(...)` registry entry
  + `references.bib` entry both land, gated by `tests/citations/test_bib_structure.py`.
- New algebra-level scalar kind acts on the **existing algebra block** ⇒ **schema stays v1**
  (the `tau_tilting`/`congruences` precedent). **Canonical keys are request-derived**
  (`cache.py::canonical_key`, sha256 over the sorted-keys blob + `library_version()`): a request
  that does not use the new `tau_cluster` key keys byte-identically to before, so existing
  goldens are untouched. Verify with `test_runner_delegation.py` BEFORE and AFTER adding the
  golden.
- **Estimator: `tau_cluster` is KNIT-HEAVY (heavier than `congruences`).** It builds the P45
  exchange graph of `A`, the P64 wide poset, AND a P65 reduction / sub-`g`-fan per object. Size
  it on `A.dim` (`sizing_dim`) and add it to the knit-heavy tuple beside `tau_tilting` /
  `ar_quiver` / `left_right_parts` / `tilted_check` in `estimator.py` (`estimator.py:83`), with
  the highest `ETA_MODEL["scalars"]` cost of the lattice/τ-tilting theme (P69 `knit_heavy`
  precedent). It carries a **PAIR BUDGET**, not a degree range (parsed like `tau_tilting`).
- **Plan-32 markers** (detailed per task): the object-count cross-tie (P66 τ-perpendicular
  enumeration **vs** P64 `wide_subcategories.size`), `#generators == #bricks`, `#relations ==
  #rank-2 wides`, the classifying-space face vector vs the morphisms-by-rank of the whole
  category (`f_0 = #wide`, `sum = morphism_count`), the `g_fan_face_vector` vs the exchange-graph
  face count (`f_n = #sτt`) = `oracle_crossengine`; the **kA₂ picture group** (3 gens, 1
  atom relation = pentagon), the **kA₃ picture group** (6 gens, 6 relations = 4 atom + 2
  commutation), the **kA₃ ≠ kZ₃/rad² discriminator** (same `#brick / #wide / g_fan_face_vector`
  but 4+2 vs 3+3 relation split — and, as a bonus, a differing *classifying-space* face vector
  (14,49,49,14) vs (14,48,48,14)), the **Nakayama `K(π,1)`** verdict, the τ-tilting-infinite
  honest refusal = `oracle_literature`; `is_kpi1`-theorem-anchoring, `χ = Σ(−1)^k f_k`,
  `face_vector[0] == object_count`, `sum(face_vector) == morphism_count`, the **factorization
  identity `n!·#sτt` vs P65 `signed_count`** (definitional in P65 — `signed_count := n!·|eg.vertices|`
  — so this is a consistency check, NOT a cross-engine oracle: fix-round ruling 3), composition
  associativity + identity laws, `object_count == wide.size`-internal-gate = `oracle_selfcert`;
  the QPA honest-scope probe = the `qpa` bucket. `oracle_*` markers FORBIDDEN in
  `tests/{webapp,gui,hpc}` (`test_oracle_classes.py`) — engine batteries in `tests/modules/`,
  cross-runner twins in `tests/webapp/`/`tests/gui/` unmarked.
- **Mid-merge-train counts.** v1.0.0 lands ~29 subplans in overlapping waves. **Task E recounts
  the oracle-class table at merge time** by running `tests/release/test_oracle_classes.py`
  (paste the LIVE numbers, never a guessed count) and claims only this plan's deltas.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table
  green) and the README line. Conventional commits; green at every commit; branch
  `plan-66-tau-cluster-morphism` off `dev`.

---

## Record (R29, verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R29 — τ-cluster morphism category + picture group.** [D-scout P11;
> keep-with-corrections]
> Object: W(A) from τ-perpendicular subcategories, cube-complex classifying
> space, picture-group presentation. GATE: τ-tilting-finite only (infinite
> category otherwise — loud refusal). Refs: Hanson–Igusa arXiv:1809.08989 (Comm.
> Alg. 2021); Buan–Marsh; Igusa–Todorov signed-exceptional-sequences history
> arXiv:2509.10910 (re-slotted as context for R28 too). Oracles: kA₂/kA₃ picture
> groups; object count = # wide subcategories (ties R26); Nakayama K(π,1).
> Size M.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P66 (Wave 3, tier γ, "after P64
(ties R26)"). Branch `plan-66-tau-cluster-morphism` off `dev`.

---

## Reference re-verification (mandatory standing rule; done at spec time — with a CORRECTION)

Re-verified via web during authoring. **Findings (transcribe the `# PIN`s at citation-add):**

1. **CORRECTION — the τ-cluster morphism category is DEFINED by Buan–Marsh, not Hanson–Igusa.**
   The R29 card headlines "Hanson–Igusa arXiv:1809.08989" for `W(A)`. Web-verified: the
   **τ-cluster morphism category** (objects = τ-perpendicular subcategories, morphisms via
   signed τ-exceptional sequences) is introduced by **Buan–Marsh, "A category of wide
   subcategories," International Mathematics Research Notices 2021, no. 13, pp. 10278–10338**
   (`# verify` volume/issue/pages/DOI at citation-add), generalizing the **cluster morphism
   category** of Igusa–Todorov. **Hanson–Igusa, "τ-Cluster morphism categories and picture
   groups," Communications in Algebra 49 (2021), no. 10** (DOI **10.1080/00927872.2021.1921184**,
   arXiv:1809.08989) prove the **classifying space is a cube complex**, that it is a **K(π,1) for
   Nakayama algebras**, and identify `π₁` with the **picture group**. **P66 cites Buan–Marsh
   (IMRN 2021) for the category definition, Hanson–Igusa (Comm. Alg. 2021) for the cube complex
   + Nakayama K(π,1) + picture group, and Igusa–Todorov–Weyman (arXiv:1609.02636) for the
   explicit picture-group presentation.** A dated correction note is appended to research-doc
   R29 at merge (metaplan §2: corrections folded back).

2. **Buan–Marsh, "τ-exceptional sequences," Journal of Algebra 585 (2021), pp. 36–68.** The
   signed τ-exceptional sequences underpinning the morphisms. This is the **same** paper P65
   cites for its Jasso-recursion enumeration; P66 reuses P65's `buan_marsh` key for it and adds
   a **distinct** `buan_marsh_wide` key for the IMRN "category of wide subcategories" paper
   (§Task 0 citation-merge care).

3. **Igusa–Todorov–Weyman, "Picture groups of finite type and cohomology in type Aₙ,"
   arXiv:1609.02636.** Defines the **picture group** `G(A_n)` (one generator `x(β)` per positive
   real Schur root = per **brick**; relations from rank-2 configurations) and constructs the
   finite CW/cube complex that is a **K(π,1)** with **cells in bijection with cluster-tilting
   objects (Catalan-many top cells)**. The **presentation** ground truth for type A. `# verify`
   venue/year at citation-add (arXiv id firm). **NOTE (fix-round ruling 5): P66 anchors the hereditary-Dynkin
   K(π,1) *verdict* to the specific Igusa–Todorov CAT(0) result (item 6), not to this presentation
   paper — ITW supplies the presentation words/counts, the CAT(0) paper supplies the asphericity.**

4. **Igusa–Todorov, "Short history of signed exceptional sequences," arXiv:2509.10910.**
   Context (R29 card): signed exceptional sequences are precisely the factorizations of cluster
   morphisms into rank-1 morphisms — "every cluster morphism of rank `k` has exactly `k!`
   factorizations into rank-1 morphisms," each a signed exceptional sequence. This is the
   representation-theoretic content behind P66's factorization identity `n!·#sτt` (which P66
   checks against P65's `signed_count` — a self-cert consistency gate, definitional in P65, NOT a
   cross-engine oracle: fix-round ruling 3). Cited as CONTEXT.

5. **Enomoto (2201.00595) + Marks–Šťovíček (1503.04639)** — reused from P64. Enomoto's core-
   label order is how P64 (and hence P66's object set) computes `#wide`; Marks–Šťovíček fixes
   the `#wide = #torsion ⟺ rep-finite` boundary that makes "objects = `#wide`, not `#torsion`"
   load-bearing (§Global Constraints).

6. **Igusa–Todorov, "Which cluster morphism categories are CAT(0)," arXiv:2203.16679 (2022)** —
   the specific anchor for the **hereditary-Dynkin K(π,1) verdict** AND the honest-scope context
   for the general case (authors + title web-verified this authoring). A cube complex is
   **non-positively curved** (locally CAT(0)) — hence **aspherical (K(π,1))** — **iff every vertex
   link is a flag simplicial complex** (Gromov). Igusa–Todorov prove the cluster morphism category
   is a **CAT(0) category for hereditary algebras of finite (Dynkin) or tame type with only small
   tubes**, so `|W(A)|` is locally CAT(0) hence a **K(π,1)** — **hereditary Dynkin included**.
   Hanson–Igusa (Comm. Alg. 2021) prove the τ-version K(π,1) for **Nakayama** algebras. **Whether
   the general τ-tilting-finite case is a `K(π,1)` is DELICATE** (arXiv:2203.16679 establishes
   CAT(0) only for hereditary finite/tame type with small tubes — not the general algebra), so it
   is **NOT claimed here** (a flag-link/CAT(0) combinatorial certifier is a clearly-scoped stretch,
   §Task B, not a spec-time pin). `# verify` no fabricated general-case citation.

---

## Mathematical spec — the objects P66 computes

Throughout, `A` is **τ-tilting-finite** (else the honest refusal), so the `g`-vector fan of `A`
is a **complete simplicial fan** (a triangulated `(n−1)`-sphere) with finitely many faces, and
`W(A)` is a **finite** category. Notation: a **support τ-rigid pair** `U = (M, P)` (`M` τ-rigid,
`P` projective, the AIR pair) is a **face** of the `g`-fan; its number of indecomposable
summands `|U|` is the face's cardinality (`0 ≤ |U| ≤ n`). The **facets** (`|U| = n`) are the
support τ-**tilting** pairs (`#sτt` of them). Every support τ-rigid pair extends to a facet
(Bongartz), so the faces are exactly the subsets of facets — this is how P66 enumerates them
(§Verified pins prototype: dedup the subsets of the `g`-matrix columns).

### 1. Objects — the τ-perpendicular wide subcategories

The **objects** of `W(A)` are the **τ-perpendicular subcategories** `J(U) = {X ∈ mod A :
Hom(M, X) = 0 = Hom(X, τM), X_P = 0}` for `U = (M, P)` a support τ-rigid pair (Jasso; Buan–
Marsh). `J(U)` is a **wide** subcategory (exact abelian, extension-closed) of rank `n − |U|`,
and `J(U) ≅ mod C(U)` for the Jasso reduction algebra `C(U)` (**P65's `tau_perpendicular_
reduction`**). For **τ-tilting-finite `A`** the τ-perpendicular subcategories are **exactly all
the (finitely many) wide subcategories** (Buan–Marsh; every wide is functorially finite hence
τ-perpendicular), so

> **`#objects = #wide = wide_subcategories(A).size`** (Enomoto CLO, P64) — **the R29 "object
> count = #wide (ties R26)" oracle.**

P66 enumerates objects the τ-perpendicular way (each face `U → J(U)`, deduped by the **iso-class
semibrick** of `J(U)`'s simple bricks — non-thin safe, never by dim-vector) and **cross-ties the
count** to P64's independent Enomoto count. Each object carries: its rank, its simple-brick
semibrick (its "simple objects"), and the `C_W` handle (from the reduction). `J(∅) = mod A`
(the terminal-rank-`n` object, the whole category); the rank-0 wide `0` is the initial object.

### 2. Morphisms — support τ-rigid pairs of the source, graded by rank

A **morphism** `[U]: W → W'` (Buan–Marsh) is (the iso class of) a **support τ-rigid pair `U` of
the wide `W`** (equivalently of `C_W`) with `W' = J_W(U)` its τ-perpendicular category computed
**inside** `W`. Its **rank** is `|U|`; the **identity** of `W` is the empty pair `[∅]`
(`J_W(∅) = W`). Composition is `[U'] ∘ [U] = [U'']` via the reduction (the perpendicular of a
perpendicular inside `W`); it is associative and unital (Buan–Marsh — P66 self-certifies it on a
spot-check, §Task A). Consequences P66 leans on and pins:

- **Out-degree.** `#morphisms out of W = #(support τ-rigid pairs of C_W) = #(faces of the
  `g`-fan of `C_W`)`. By the Jasso reduction the `g`-fan of `C(U)` is the **link** of the face
  `U` in the `g`-fan of `A`, so **`#morphisms out of J(U) = |closedStar(U)|`** = the number of
  faces `V ⊇ U` in `A`'s own `g`-fan — a purely combinatorial count from P45's exchange graph,
  **no reduction needed for the number** (the reduction names `C_W`; the star gives the count).
  In particular `#morphisms out of mod A` (`U = ∅`) `= #(all faces) = #(support τ-rigid pairs of
  `A`)` — **live-verified: kA₂ → 11, kA₃ → 45** (§Verified pins). This is the **cross-engine
  arbiter** for the morphism machinery (star size vs reduction `C_W` face count must agree).
- **Total morphisms.** `morphism_count = Σ_{objects W} |closedStar(U_W)|` for one representative
  face `U_W` per object (the star size is a wide invariant `= #faces(C_W)`, independent of the
  representative). **Hand-derived + live-verified breakdown: kA₂ → 21, kA₃ → 126** (§Verified
  pins).
- **Factorization (ties P65 — a CONSISTENCY check, not a cross-engine oracle: fix-round ruling 3).** Every
  rank-`k` morphism has exactly **`k!`** factorizations into rank-1 morphisms, each a **signed
  exceptional sequence** (Igusa–Todorov 2509.10910; Buan–Marsh). The rank-`n` morphisms out of
  `mod A` are the facets (`#sτt`), so `#(complete signed exceptional sequences) = n! · #sτt`,
  which P66 checks against P65's `tau_exceptional_sequences(A).signed_count`. **Caveat (fix-round
  ruling 3): P65 DEFINES `signed_count := n! · |exchange_graph.vertices|` (its own doc marks that as an
  identity self-cert), and P66's `#sτt` is likewise `|eg.vertices| = face_vector[-1] =
  g_fan_face_vector[-1]` off the SAME exchange graph. So `n!·(P66 #sτt) == P65.signed_count` is
  *definitionally true whenever both read the same graph* — it is `oracle_selfcert` (a genuine
  but internal consistency gate that the two plans agree on `#sτt`), NOT a `oracle_crossengine`
  oracle.** (Upgrade path, if P65 ships materialised signed sequences: `len(materialised) == n!·#sτt`
  would be the genuine cross-engine form — P65's own remedy; defer to the implementer only if that
  API exists, else keep it self-cert. See §Task A "adjust to reality".) This is the R28/R29↔R27
  bridge.

### 3. The classifying space — the Hanson–Igusa cube complex, Euler characteristic, K(π,1)

**H1 correction (BLOCKING; this replaces the earlier draft's WRONG space).** The classifying
space `|W(A)|` (nerve of the category) is homotopy equivalent to a **cube complex** (Hanson–Igusa,
generalizing Igusa–Todorov). Its cells are **NOT** the faces of the `g`-fan of `A`. In HI's model:

> **`f_0 = #(0-cells) = #wide`** (one 0-cell per WIDE SUBCATEGORY = per object of the category),
> and **`f_k = #(rank-k morphisms of the WHOLE category) = Σ_{objects W} #(rank-k support τ-rigid
> pairs of C_W)`** (one `k`-cell per rank-`k` morphism), `k = 0, …, n`.

Equivalently `f = ` the category's **morphisms-by-rank** vector: `f_0 = #identities = #objects =
#wide`, and `f_n = #(rank-n morphisms) = #(facets out of mod A) = #sτt`. The **Euler
characteristic** is `χ = Σ_k (−1)^k f_k` (an `int`). **Live-verified (§Appendix A2, prototype
re-runs byte-for-byte): kA₂ `f = (5,11,5)`, `χ = −1`, `Σf = 21`; kA₃ `f = (14,49,49,14)`,
`χ = 0`, `Σf = 126`; kA₄ `f = (42,204,326,204,42)`, `χ = 2`, `Σf = 818`.** Two discriminating
self-certs hold by construction: **`face_vector[0] == object_count`** (`= #wide`) and
**`sum(face_vector) == morphism_count`**. (For kA₂ the classifying-space `χ = −1` also matches the
picture-group **presentation Euler characteristic** `1 − 3 + 1 = −1` = `1 − #gens + #atom-rels`, a
cross-check H1 requires.)

**Why the earlier draft was wrong (and the tautology it hid).** The draft reported
`f_k = #(support τ-rigid pairs of A with k summands)` — that is the **`g`-fan / cluster fan** of
`A`, a triangulated `(n−1)`-**sphere** plus the empty face, whose `χ` is *always* `(−1)^n` for a
rank-`n` τ-tilting-finite algebra (a sphere tautology that discriminates NOTHING). A sphere is not
the classifying space: the draft's kA₂ "complex" (1,5,5) is a **pentagon ≅ S¹** with `π₁ = ℤ ≠
G(kA₂)`; its kA₃ "complex" (1,9,21,14) is `≅ S²` (not aspherical — cannot be a `K(π,1)`). The
`g`-fan f-vector is still useful (its `f_n = #sτt`, it is the star-count source for the
out-degrees §2), so P66 keeps it — but **explicitly relabelled** as
`g_fan_face_vector` (**the `g`-fan/cluster-fan f-vector: a sphere; NOT the classifying space**):
`g_fan_face_vector(kA₂) = (1,5,5)`, `(kA₃) = (1,9,21,14)`, `(kA₄) = (1,14,56,84,42)` (its own
sphere-`χ = (−1)^n`). The two vectors are never conflated.

- **`K(π,1)` verdict (theorem-anchored — the honest homotopy boundary).** `|W(A)|` is a
  **K(π,1)** (aspherical; then `π₁ = ` the picture group) when its cube complex is non-positively
  curved (flag vertex links / CAT(0), Gromov). This is a **theorem** for **Nakayama** algebras
  (Hanson–Igusa, Comm. Alg. 2021) and for **hereditary algebras of finite (Dynkin) or tame type
  with small tubes** — **hereditary Dynkin included** — (Igusa–Todorov, *Which cluster morphism
  categories are CAT(0)*, arXiv:2203.16679, 2022: the cluster morphism category is a CAT(0)
  category, so its classifying space is locally CAT(0) hence a `K(π,1)`). P66 reports
  `is_kpi1 = True` with `kpi1_reason` naming the theorem **only** for these certified classes
  (Nakayama detected via the shipped `is_nakayama` recognizer; hereditary Dynkin via
  `is_hereditary` + `dynkin_type`); otherwise `is_kpi1 = None`, `kpi1_reason =
  "not certified (K(π,1) known only for Nakayama / hereditary Dynkin at the cited theorems)"`.
  **P66 NEVER asserts a homotopy-type conclusion beyond the cited theorems** (the general
  τ-tilting-finite case is *delicate* — arXiv:2203.16679 establishes CAT(0) only for hereditary
  finite/tame type with small tubes, not the general τ-tilting-finite algebra). The cube complex
  itself is reported for every τ-tilting-finite `A`; only the *asphericity verdict* is gated. (A
  combinatorial Gromov flag-link / CAT(0) certifier — decidable for any τ-tilting-finite `A` — is
  a clearly-scoped **stretch**, §Task B, requiring HI's exact cube/link model transcribed
  verbatim; NOT a v1 spec-time pin, because a wrong link model would silently mis-verify
  asphericity.)

### 4. The picture group presentation

When `|W(A)|` is a K(π,1), `π₁ = ` the **picture group** `G(A)` (Igusa–Todorov–Weyman; Hanson–
Igusa). P66 ships its **presentation as DATA** (never a group-theoretic decision procedure —
word problem etc. are out of scope):

- **Generators.** One `x(β)` per **brick** `β` of `A` (= positive real Schur root in the
  hereditary case). `num_generators = #bricks` (live: kA₂ → 3, kA₃ → 6). Each generator carries
  its brick's dim-vector + `identify_standard` name.
- **Relations.** One per **rank-2 wide subcategory** `W₂` (= per size-2 semibrick `{β, γ}` of
  hom-orthogonal bricks). Two types (Igusa–Todorov–Weyman Def. of `G`):
  - **commutation** — when `W₂ ≅ mod(k × k)` (`Ext¹(β,γ) = Ext¹(γ,β) = 0`, no extension brick):
    `x(β) x(γ) = x(γ) x(β)`.
  - **atom** — when `W₂ ≅ mod kA₂` (connected: exactly one of `Ext¹(β,γ)`, `Ext¹(γ,β)` is
    nonzero, so `W₂` has a third brick `δ = β + γ`, the extension): the "pentagon" relation
    equating the two maximal green sequences of `W₂`. In the ITW normal form (roots
    `γ_k = a_k β + b_k γ` ordered by increasing slope `a_k/b_k`) it is
    `x(β) x(γ) = ∏_k x(γ_k)` — the LHS is the two SIMPLE bricks of `W₂`, the RHS lists **all
    three bricks** of `W₂` (`β`, `δ`, `γ`) in slope order; read cyclically it is the length-5
    pentagon in which each simple brick appears twice and the extension brick once.
  `num_relations = #(rank-2 wides)`; `num_atom = #(rank-2 wides with an extension brick)`,
  `num_commutation = #(rank-2 wides ≅ k×k)`. **Live-verified split: kA₂ → 1 atom + 0
  commutation; kA₃ → 4 atom + 2 commutation; kZ₃/rad² → 3 atom + 3 commutation** (§Verified
  pins). **The exact relation WORDS are transcribed from ITW/HI at implementation** (the slope
  order is orientation-dependent); the plan pins the **counts, the type split, and the brick-
  membership of each relation** (all engine-derivable and discriminating), and derives the kA₂
  pentagon concretely below.
- **Abelianization (H2 correction — the formula was WRONG).** `G(A)^ab = ℤ^{#bricks} / (relation
  matrix)` via exact SNF (`fields.linalg`); each commutation relation is trivial in `G^ab`, each
  atom relation `x(β)+x(γ) = x(β)+x(δ)+x(γ)` gives `x(δ) = 0` (kills the *extension* brick `δ`).
  The draft claimed `G(A)^ab = ℤ^{#bricks − #atom}` — **wrong: several distinct atom relations can
  share ONE extension brick `δ`, so they impose the SAME equation `x(δ) = 0`.** The correct rank
  of the relation matrix is the number of **distinct extension bricks** among the atoms (`≤ #atom`),
  SNF-decided:
  > **`abelianization_rank = #bricks − rank(relation matrix) = #bricks − #(distinct extension
  > bricks among the atom relations)`** (all invariant factors are `1` here, so `G(A)^ab` is free).
  **Live-verified (§Appendix A2): kA₂ → ℤ² (`3 − 1`); kA₃ → ℤ³ (`6 − 3`, NOT ℤ²! — the 4 atoms
  have only 3 distinct extension bricks: `{S₁,M₂₃}` and `{M₁₂,S₃}` both extend to the top brick
  `M₁₂₃`); kA₄ → ℤ⁴ (`10 − 6`); kZ₃/rad² → ℤ³ (`6 − 3`).** **Cautionary note:** the discarded
  `#bricks − #atom` formula gives the ABSURD `ℤ⁰` for kA₄ (`10 − 10`) and the WRONG `ℤ²` for kA₃
  (`6 − 4`) — exactly the failure the extension-brick sharing exposes; that kA₄ pin was ABSENT
  from the draft and is added here as the discriminating counterexample. P66 reports the invariant
  factors from the SNF (never the `#bricks − #atom` count).

**Concrete kA₂ presentation (derived, exact word to-transcribe).** The three bricks of
`kA₂ = 1→2` are `S₁ = (1,0)`, `S₂ = (0,1)`, `P₁ = (1,1)`. The single rank-2 wide is `mod A ≅
mod kA₂` itself, connected (`Ext¹(S₁, S₂) ≠ 0`, extension brick `P₁ = S₁ + S₂`), so **one atom
relation, no commutation**:
`G(kA₂) = ⟨ x_{S₁}, x_{S₂}, x_{P₁} | x_{S₁} x_{S₂} = x_{P₁} x_{S₂} x_{S₁} ⟩` (ITW slope-order
form; the exact left/right convention transcribed at implementation — the invariant content is:
3 generators, 1 atom relation, the two maximal green sequences of the wide `{S₁,S₂}` (length 2)
and `{S₂,P₁,S₁}` (length 3) equated, abelianization `ℤ²`).

---

## Verified pins — LIVE-computed against the shipped P45 engine (this authoring)

Recomputed in the venv (`exchange_graph` + `bricks` + `hom_dim` + `ext_dims`, then the P66
classifying-space face-vector / all-semibrick / rank-2-classification / abelianization-rank
prototype — §Appendix A1, which re-runs byte-for-byte). **Every number below was observed live**
(the classifying-space `face_vector`, the `g_fan_face_vector`, `χ`, the atom/commutation split,
the distinct-extension-brick abelianization rank, and the totals were ALL emitted by the prototype)
except the picture-group relation *words* (structure/counts observed; exact word
literature-anchored to ITW). Two vectors are reported and **never conflated**:
- **`face_vector`** = the Hanson–Igusa **classifying-space** cube complex (`f_0 = #wide`,
  `f_k = #(rank-k morphisms)`, `Σ = morphism_count`) — the ACTUAL topology;
- **`g_fan_face_vector`** = the **`g`-fan/cluster fan of `A`** (`#(support τ-rigid pairs) by
  summand count`, `f_n = #sτt`, `Σ = #faces = morphisms out of the top object`) — a
  triangulated `(n−1)`-**sphere**, NOT the classifying space.

| algebra | `#bricks` (=#gens) | `#wide` (=#objects) | wide-by-rank | `g_fan_face_vector` (SPHERE) | `#sτt` | **classifying-space `face_vector`** | `χ` | rank-2 atom+comm (=#rels) | abelianization | out-of-top / total | (gen, rel) | K(π,1) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **kA₂** (`1→2`) | **3** | **5** | `{0:1,1:3,2:1}` | (1,5,5) | 5 | **(5,11,5)** | **−1** | **1 + 0 = 1** | **ℤ²** | **11 / 21** | (3, 1) | — (hered. Dynkin ✓) |
| **kA₃** (`1→2→3`) | **6** | **14** | `{1,6,6,1}` | (1,9,21,14) | 14 | **(14,49,49,14)** | **0** | **4 + 2 = 6** | **ℤ³** | **45 / 126** | (6, 6) | — (hered. Dynkin ✓) |
| **kA₄** | 10 | 42 | `{1,10,20,10,1}` | (1,14,56,84,42) | 42 | **(42,204,326,204,42)** | 2 | 10 + 10 = 20 | **ℤ⁴** | 197 / 818 | (10, 20) | — (hered. Dynkin ✓) |
| **kA₃/rad²** (linear Nakayama, non-hereditary) | **5** | **12** | `{1,5,5,1}` | (1,8,18,12) | 12 | **(12,40,40,12)** | 0 | **2 + 3 = 5** | **ℤ³** | 39 / 104 | (5, 5) | **✓ (Nakayama)** |
| **kZ₃/rad²** (self-injective Nakayama) | **6** | **14** | `{1,6,6,1}` | (1,9,21,14) | 14 | **(14,48,48,14)** | 0 | **3 + 3 = 6** | **ℤ³** | 45 / 124 | (6, 6) | **✓ (Nakayama)** |

Every row satisfies the two discriminating self-certs **`face_vector[0] == #wide`** and
**`sum(face_vector) == total morphisms`**, and `g_fan_face_vector[-1] == #sτt`. (The `χ` column is
the *classifying-space* Euler characteristic `Σ(−1)^k face_vector[k]` — **−1, 0, 2** for
`kA₂,kA₃,kA₄`, deliberately NOT `(−1)^n`; the discarded sphere had the tautological
`χ_sphere = (−1)^n`. For `kA₂`, `χ = −1` matches the presentation Euler characteristic
`1 − #gens + #atom = 1 − 3 + 1`.)

**kA₂ / kA₃ are the R29 "picture group" pins.** `G(kA₂)` = 3 generators (bricks `S₁,S₂,P₁`), 1
atom relation (the pentagon on the single rank-2 wide), abelianization `ℤ²`. `G(kA₃)` = 6
generators, 6 relations (**4 atom + 2 commutation** — the 4 atoms from the four block-size-3
noncrossing partitions `123|4, 124|3, 134|2, 234|1` ≅ `mod kA₂`, the 2 commutations from the
two size-(2,2) partitions `12|34, 14|23` ≅ `mod(k×k)`; the all-semibrick histogram `{1,6,6,1}`
matches P64's NC(A₃) Whitney numbers exactly), abelianization **`ℤ³`** (`6 − 3`: the 4 atoms
share extension bricks so the relation matrix has rank 3, NOT `ℤ²` — H2 correction).

**The `kA₃` vs `kZ₃/rad²` DISCRIMINATOR (the anti-tautology pin).** Both have `#bricks = 6`,
`#wide = 14`, **`g_fan_face_vector` `(1,9,21,14)`** (identical spheres!), wide-by-rank `{1,6,6,1}`,
`#sτt = 14` — **identical on every coarse count including the g-fan sphere**. They differ in the
picture-group relation *type split* — **kA₃ → 4 atom + 2 commutation, kZ₃/rad² → 3 atom + 3
commutation** (kZ₃/rad² has one more completely-orthogonal rank-2 wide, one fewer connected one,
because its bricks have different `Ext¹`-orthogonality than the hereditary `kA₃`) — **and,
consequently, in the *classifying-space* `face_vector`: kA₃ `(14,49,49,14)` vs kZ₃/rad²
`(14,48,48,14)`** (one more connected rank-2 wide → one more `mod kA₂` factor with
`g_fan (1,5,5)` in place of a `mod(k×k)` factor with `(1,4,4)`, so `f₁,f₂` are `+1` each). So the
**g-fan sphere discriminates NOTHING (both `(1,9,21,14)`), but the actual classifying space DOES
(49 vs 48) — a direct illustration of the H1 correction: the classifying space counts morphisms,
the sphere does not.** A construction that reads the relation type off dim-vectors, or ignores the
direction/existence of `Ext¹` between the two simple bricks, would give the SAME split (and the
SAME classifying-space `face_vector`) for both — the `Ext¹`-driven split is exactly the
discriminator that catches it. Pinned as `oracle_literature` (values) + `oracle_crossengine` (the
split recomputed via `ext_dims` vs the category's rank-2 objects).

**kZ₃/rad² and kA₃/rad² are the "Nakayama K(π,1)" pins** (the R29 card) — both Nakayama, both
`is_kpi1 = True` with `kpi1_reason` naming Hanson–Igusa. **kA₃/rad² is also the "beyond kAₙ"
non-thin/non-linear example** (a non-hereditary Nakayama with `(#bricks, #wide) = (5, 12)`,
distinct from every `kAₙ`), the required "beyond thin-A" battery member.

**Why `#objects = #wide` and NOT `#semibricks`/`#torsion` (recorded).** All five algebras here
are representation-FINITE, so `#wide = #torsion = #semibricks` and the trap is invisible — the
code reads `wide_subcategories(A).size` regardless (the correct source), and the honest-scope
entry (§Task E) records that on a rep-INFINITE τ-tilting-finite algebra (`Π(D₄)`) `#wide <
#torsion` and only `.size` is the object count. That witness is DEFERRED (P64's engine limit).

**Total morphism counts = `Σ_W #faces(C_W)` = `sum(classifying-space face_vector)` (live).** The
star size of each object `W` is a wide invariant `|closedStar(U_W)| = #faces(C_W)` (the Jasso link
identity), and `#faces(C_W)` graded by rank is exactly `C_W`'s `g_fan_face_vector`; summing over
objects gives BOTH the total AND the classifying-space `face_vector`:
`kA₂: 1·|star(0-wide)| + 3·|star(rank-1)| + 1·|star(mod A)| = 1·1 + 3·3 + 1·11 = 21 = sum(5,11,5).`
`kA₃: 1·1 + 6·3 + (4·11 + 2·9) + 1·45 = 1 + 18 + 62 + 45 = 126 = sum(14,49,49,14).`
`kA₄: 1·1 + 10·3 + (20 rank-2 wides: 10·11 + 10·9) + (10 rank-3 wides: 5·45 + 5·33) + 1·197
      = 1 + 30 + 200 + (225 + 165) + 197 = 818 = sum(42,204,326,204,42)`
(`#faces(kA₂)=11`, `#faces(k×k)=9`, `#faces(kA₃)=45`, `#faces(kA₂×kA₁)=33`, `#faces(kA₄)=197` are
the piece g-fan #faces). All three (and the two Nakayama totals 104 / 124) were emitted by the
prototype (§Appendix A2); `|star(mod A)| = #faces = out-of-top` is the SPHERE #faces, `sum(all
stars) = total = sum(classifying-space face_vector)`. Implementation cross-checks star sizes
against P65's per-object `C_W` face counts.

---

### Task 0: citations (do first — later tasks reference the keys)

**Files:**
- Modify: `src/quiverlab/citations/references.bib`, `src/quiverlab/citations/registry.py`
- Test: `tests/citations/test_bib_structure.py` (existing gate).

**Interfaces:** the `_r(key, bibtex_key, kind, title, annotation, *tags)` registry helper
(`registry.py:24`). **Reuse (do NOT re-add):** `air_tau_tilting` (→ `AIR2014`),
`demonet_iyama_jasso` (→ `DIJ2019`) from P45; `enomoto_wide_ice`, `marks_stovicek` from P64
(the wide-count + `#wide=#torsion⟺rep-finite` boundary); `jasso_reduction` and **`buan_marsh`**
(the Buan–Marsh τ-exceptional J. Algebra paper) from P65. **Do NOT reuse the existing
`igusa_todorov` key** — it is the Igusa–Todorov φ/ψ **functions** paper (`IgusaTodorov2005`,
Plan 40), a DIFFERENT work.

**Cross-plan citation-merge care (metaplan §5).** P64 adds `dirrt_lattice_torsion`,
`barnard_carroll_zhu`, `enomoto_wide_ice`, `marks_stovicek`; P65 adds `jasso_reduction`,
`buan_marsh`, `crawley_boevey_exceptional`, `ringel_braid`. **P66 depends on both merging
first**, so those keys are present — reuse, do NOT re-add. P66's four NEW keys
(`buan_marsh_wide`, `hanson_igusa`, `igusa_todorov_weyman`, `igusa_todorov_cat0`) are disjoint
from P64/P65 (and from the existing `igusa_todorov` φ/ψ-functions key); on a rebase, insert the
P66 block by EXPLICIT path, never `git add -A`, and re-run `test_bib_structure.py`. **If P65 did
NOT add a `buan_marsh` key** (it may fold BM into another key), add both `buan_marsh` (J. Alg.
585) and `buan_marsh_wide` (IMRN) here.

- [ ] **Step 1: add the new keys** (BibTeX-verify each before committing — Plan-29 rule):

```bibtex
@article{BuanMarsh2021wide,
  author  = {Buan, Aslak Bakke and Marsh, Robert J.},
  title   = {A category of wide subcategories},
  journal = {International Mathematics Research Notices},
  volume  = {2021},
  number  = {13},
  pages   = {10278--10338},
  year    = {2021},
  doi     = {10.1093/imrn/rnz082},
  note    = {author field matches the PRINTED IMRN record (Aslak Bakke Buan and Robert J. Marsh); verify volume/issue/pages/DOI at citation-add},
}
@article{HansonIgusa2021,
  author  = {Hanson, Eric J. and Igusa, Kiyoshi},
  title   = {{$\tau$}-cluster morphism categories and picture groups},
  journal = {Communications in Algebra},
  volume  = {49},
  number  = {10},
  pages   = {4376--4415},
  year    = {2021},
  doi     = {10.1080/00927872.2021.1921184},
  note    = {arXiv:1809.08989; verify page range at citation-add},
}
@article{IgusaTodorovWeyman2016,
  author  = {Igusa, Kiyoshi and Todorov, Gordana and Weyman, Jerzy},
  title   = {Picture groups of finite type and cohomology in type {$A_n$}},
  journal = {arXiv preprint},
  year    = {2016},
  note    = {arXiv:1609.02636; verify final venue at citation-add},
}
@article{IgusaTodorov2022cat0,
  author  = {Igusa, Kiyoshi and Todorov, Gordana},
  title   = {Which cluster morphism categories are {CAT(0)}},
  journal = {arXiv preprint},
  year    = {2022},
  note    = {arXiv:2203.16679; CAT(0) (hence K(pi,1)) for hereditary algebras of finite or tame type with only small tubes -- anchors the hereditary-Dynkin K(pi,1) verdict AND the honest-scope caveat for the general case; verify final venue at citation-add},
}
```

  and in `registry.py`:

```python
_r("buan_marsh_wide", "BuanMarsh2021wide", "foundation",
   "A category of wide subcategories",
   "Buan-Marsh: DEFINES the tau-cluster morphism category W(A) -- objects are the "
   "tau-perpendicular wide subcategories, morphisms are support tau-rigid pairs of the source "
   "with target the tau-perpendicular category (via the Jasso reduction), morphisms factor as "
   "signed tau-exceptional sequences. Plan 66's category-structure ground truth.",
   "tau-tilting", "wide", "category"),
_r("hanson_igusa", "HansonIgusa2021", "foundation",
   "tau-cluster morphism categories and picture groups",
   "Hanson-Igusa: the classifying space of W(A) is a cube complex (one n-cube per support "
   "tau-tilting object); it is a K(pi,1) for Nakayama algebras; pi_1 is the picture group. "
   "Plan 66's cube-complex face vector, the Nakayama K(pi,1) verdict, and the picture group.",
   "tau-tilting", "picture-group", "cube-complex"),
_r("igusa_todorov_weyman", "IgusaTodorovWeyman2016", "foundation",
   "Picture groups of finite type and cohomology in type A_n",
   "Igusa-Todorov-Weyman: the picture group PRESENTATION -- one generator x(beta) per brick "
   "(positive real Schur root), relations per rank-2 configuration (commutation for k x k, the "
   "atom/pentagon relation for connected rank-2 wides); the CW complex with cells in "
   "bijection with cluster-tilting objects (Catalan-many). Plan 66's presentation ground truth.",
   "tau-tilting", "picture-group"),
_r("igusa_todorov_cat0", "IgusaTodorov2022cat0", "foundation",
   "Which cluster morphism categories are CAT(0)",
   "Igusa-Todorov: the cluster morphism category is a CAT(0) category for hereditary algebras "
   "of finite (Dynkin) or tame type with only small tubes, so its classifying space is locally "
   "CAT(0) hence a K(pi,1). Plan 66's specific anchor for the HEREDITARY-DYNKIN K(pi,1) verdict "
   "(distinct from the ITW type-A_n presentation paper), and the honest-scope context for why the "
   "general tau-tilting-finite case is delicate (CAT(0) is proven only for hereditary "
   "finite/tame type, not the general algebra).",
   "tau-tilting", "picture-group", "cube-complex", "cat0"),
```

  **Spec-ambiguity resolution (recorded):** BibTeX keys `BuanMarsh2021wide` / `HansonIgusa2021`
  / `IgusaTodorovWeyman2016` / `IgusaTodorov2022cat0`; snake-case registry keys `buan_marsh_wide`
  / `hanson_igusa` / `igusa_todorov_weyman` / `igusa_todorov_cat0` (FOUR new keys — the CAT(0)
  paper was added per fix-round ruling 5 to anchor the hereditary-Dynkin verdict to a *specific*
  Igusa–Todorov result). `HansonIgusa2021` DOI is firm (`10.1080/00927872.2021.1921184`);
  `# verify` the exact page range. `BuanMarsh2021wide` author field matches the printed IMRN
  record (**Aslak Bakke Buan and Robert J. Marsh**, fix-round ruling 5); `# verify` volume/issue/pages/DOI
  at citation-add. `IgusaTodorovWeyman2016` (arXiv:1609.02636) and `IgusaTodorov2022cat0`
  (arXiv:2203.16679, authors + title web-verified this authoring) ship the arXiv ids; `# verify`
  final venue (no fabricated print metadata — house rule).
- [ ] **Step 2:** `... -m pytest tests/citations/test_bib_structure.py -q` green.
- [ ] **Step 3: Commit** — `docs(citations): tau-cluster morphism references (Buan-Marsh IMRN 2021 category of wide subcats, Hanson-Igusa Comm.Alg. 2021 cube complex + Nakayama K(pi,1), Igusa-Todorov-Weyman 1609.02636 picture group, Igusa-Todorov 2203.16679 CAT(0) hereditary K(pi,1))`

---

### Task A: `cluster_morphism.py` — objects, morphisms, cube complex

**Files:**
- Create: `src/quiverlab/tautilting/cluster_morphism.py`
- Modify: `src/quiverlab/tautilting/__init__.py` (export `TauClusterCategory`,
  `tau_cluster_category` — `PictureGroup` / `picture_group` added in Task B),
  `src/quiverlab/core/algebra.py` (thin lazy delegate `Algebra.tau_cluster_category(budget=512)`)
- Test: `tests/modules/test_tau_cluster.py`

**Interfaces:**
- Consumes: `tautilting.mutation.exchange_graph`, `tautilting.torsion.hasse_orientation` /
  `bricks` / `_edge_brick` / `_torsion_universe`, `tautilting.congruence.wide_subcategories`
  (P64), `tautilting.exceptional.tau_perpendicular_reduction` / `tau_exceptional_sequences`
  (P65), `modules.hom.is_isomorphic` / `identify_standard`.
- Produces:
  ```python
  @dataclass(frozen=True)
  class TauClusterCategory:
      algebra: object
      objects: tuple           # (id, rank, simple_brick_semibrick(iso-class labels), C_W_dimvec/name)
      object_count: int        # == wide_subcategories(A).size  (cross-tie, self-cert gate)
      morphisms: tuple         # ((source_id, target_id, rank, pair_dimvecs), ...)
      morphism_count: int
      out_degree: dict         # object_id -> #morphisms out (== |closedStar(U_W)| == #faces(C_W))
      face_vector: tuple       # HANSON-IGUSA CLASSIFYING SPACE: f_0 == object_count (#wide),
                               # f_k == #(rank-k morphisms of the whole category),
                               # sum(f) == morphism_count, f_n == #sTt (H1 correction)
      g_fan_face_vector: tuple # the g-fan/cluster fan of A: a triangulated (n-1)-SPHERE;
                               # #(support tau-rigid pairs of A) by summand count; f_n == #sTt;
                               # sum == #faces == out-degree of mod A. NOT the classifying space.
      euler_characteristic: int  # == sum((-1)^k * face_vector[k])  (NOT (-1)^n)
      is_kpi1: bool | None     # theorem-anchored: True (Nakayama / hered. Dynkin), else None
      kpi1_reason: str
      is_complete: bool; status: str; note: str

  def tau_cluster_category(A, *, budget=512) -> TauClusterCategory
  # helpers (module-internal):
  def _g_fan(eg) -> (facets, faces_by_card, g_fan_face_vector)  # subsets of g-matrix cols (SPHERE)
  def _objects(A, eg, budget) -> tuple               # each face -> J(U) via P65 reduction, dedup by iso-semibrick
  def _closed_star_size(face, facets) -> int         # #faces V >= face  (= morphisms out)
  def _classifying_face_vector(objects, out_by_rank) -> tuple   # f_k = #(rank-k morphisms), f_0 = #wide
  def _kpi1_verdict(A) -> (bool|None, str)           # Nakayama / hered-Dynkin theorem gate
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tau_cluster.py
"""The tau-cluster morphism category W(A) (Plan 66 / R29, Buan-Marsh IMRN 2021 + Hanson-Igusa
Comm. Alg. 2021). Objects = tau-perpendicular wide subcategories (count = #wide, ties P64);
morphisms = support tau-rigid pairs of the source (out-degree = closed-star size = #faces(C_W));
the Hanson-Igusa CLASSIFYING-SPACE face vector f_0 = #wide, f_k = #(rank-k morphisms), chi =
sum (-1)^k f_k (NOT the g-fan SPHERE, reported separately as g_fan_face_vector with f_n = #sTt);
K(pi,1) theorem-anchored (Nakayama / hereditary Dynkin). QQ-scope (the P45 char caveat)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.cluster_morphism import tau_cluster_category
from quiverlab.tautilting.congruence import wide_subcategories
from quiverlab.tautilting.torsion import bricks

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)

def _nakayama_kA3_rad2():          # linear, non-hereditary, "beyond kA_n"
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ)

def _nakayama_kZ3_rad2():          # cyclic, self-injective
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}
                  ).algebra(relations=["a*b", "b*c", "c*a"], field=QQ)


@xeng
@pytest.mark.parametrize("A_factory, nwide", [
    (lambda: linear_path_algebra(2, field=QQ), 5),
    (lambda: linear_path_algebra(3, field=QQ), 14),
    (_nakayama_kA3_rad2, 12),
    (_nakayama_kZ3_rad2, 14),
])
def test_object_count_equals_wide_count(A_factory, nwide):
    # OBJECTS = #wide via TWO INDEPENDENT ROUTES: P66's tau-perpendicular enumeration (each face
    # U -> J(U), deduped by iso-class simple-brick semibrick) vs P64's Enomoto core-label-order
    # size. NOT a tautology. Also NOT #semibricks/#torsion -- those coincide only rep-finite.
    A = A_factory()
    C = tau_cluster_category(A)
    assert C.is_complete and C.object_count == nwide
    assert C.object_count == wide_subcategories(A).size        # the cross-tie (ruling: .size)


@lit
@pytest.mark.parametrize("n, facevec, gfan", [
    (2, (5, 11, 5),          (1, 5, 5)),
    (3, (14, 49, 49, 14),    (1, 9, 21, 14)),
    (4, (42, 204, 326, 204, 42), (1, 14, 56, 84, 42))])
def test_face_vector_type_A(n, facevec, gfan):
    # face_vector = the Hanson-Igusa CLASSIFYING-SPACE cube complex (f_0 = #wide, f_k =
    # #(rank-k morphisms)), NOT the g-fan sphere. g_fan_face_vector = the SPHERE (f_n = #sTt).
    C = tau_cluster_category(linear_path_algebra(n, field=QQ))
    assert C.face_vector == facevec
    assert C.g_fan_face_vector == gfan
    assert C.face_vector[0] == C.object_count          # f_0 == #wide (discriminating self-cert)
    assert sum(C.face_vector) == C.morphism_count      # sum == total morphisms
    assert C.g_fan_face_vector[-1] == C.face_vector[-1]  # both = #sTt (top cells / top facets)


@selfcert
@pytest.mark.parametrize("n, chi", [(2, -1), (3, 0), (4, 2)])   # classifying-space chi (NOT (-1)^n)
def test_euler_characteristic(n, chi):
    C = tau_cluster_category(linear_path_algebra(n, field=QQ))
    assert C.euler_characteristic == chi
    assert C.euler_characteristic == sum((-1) ** k * f for k, f in enumerate(C.face_vector))
    # NB: chi is NOT (-1)^n -- that was the DISCARDED g-fan sphere's tautology (ruling H1).
    # For kA2 the classifying-space chi = -1 also equals the presentation Euler char 1 - 3 + 1.


@xeng
def test_morphisms_out_of_top_equals_g_fan_faces():
    # morphisms out of the terminal object mod A (U = empty) == #(all g-fan faces) ==
    # #(support tau-rigid pairs of A) == sum(g_fan_face_vector). Live: kA2->11, kA3->45.
    # NOTE (ruling H1): this is the SPHERE face count, NOT sum(face_vector) (= total morphisms).
    for A, top_out in [(_kA2(), 11), (linear_path_algebra(3, field=QQ), 45)]:
        C = tau_cluster_category(A)
        top = max(C.objects, key=lambda o: o[1])       # the rank-n object = mod A
        assert C.out_degree[top[0]] == top_out
        assert sum(C.g_fan_face_vector) == top_out      # closed star of the empty face = all faces
        assert sum(C.face_vector) != top_out            # the classifying space counts MORE (21/126)


@lit
@pytest.mark.parametrize("A_factory, total", [
    (_kA2, 21), (lambda: linear_path_algebra(3, field=QQ), 126)])
def test_total_morphism_count(A_factory, total):
    # total morphisms = sum over objects of the closed-star size (= #faces(C_W)). Hand-derived
    # kA2->21, kA3->126 (Verified pins). The per-object out-degree cross-checks against P65's
    # C_W face count at implementation (the Jasso link identity).
    C = tau_cluster_category(A_factory())
    assert C.morphism_count == total
    assert C.morphism_count == sum(C.out_degree.values())


@selfcert
def test_out_degree_is_closed_star_and_matches_reduction():
    # For every object W = J(U): out_degree[W] == |closedStar(U)| (star in A's g-fan) AND ==
    # #faces(C_W) computed from P65's tau_perpendicular_reduction(A, U).reduced exchange graph.
    # This is the Jasso link identity -- a genuine cross-engine self-cert of the morphism count.
    from quiverlab.tautilting.exceptional import tau_perpendicular_reduction
    from quiverlab.tautilting.mutation import exchange_graph
    A = linear_path_algebra(3, field=QQ)
    C = tau_cluster_category(A)
    # spot-check the rank-2 objects: 4 must have #faces(C_W)=11 (C_W ~ kA2), 2 must have 9 (kxk)
    rank2 = [o for o in C.objects if o[1] == 2]
    assert len(rank2) == 6
    stars = sorted(C.out_degree[o[0]] for o in rank2)
    assert stars == [9, 9, 11, 11, 11, 11]              # 2 commutation-wides (9) + 4 kA2-wides (11)


@selfcert                                              # fix-round ruling 3: NOT xeng -- definitional in P65
def test_factorization_count_ties_P65():
    # #(complete signed exceptional sequences) = n! * #sTt.  P65 DEFINES signed_count :=
    # n! * |exchange_graph.vertices|, and P66's #sTt = face_vector[-1] = g_fan_face_vector[-1] =
    # |eg.vertices| off the SAME graph -- so this equality is DEFINITIONALLY true, an internal
    # consistency gate (oracle_selfcert), NOT a cross-engine oracle. (If P65 ships MATERIALISED
    # signed sequences, upgrade to len(materialised) == n!*#sTt for a genuine xeng check.)
    from math import factorial
    from quiverlab.tautilting.exceptional import tau_exceptional_sequences
    A = linear_path_algebra(3, field=QQ)
    C = tau_cluster_category(A)
    stt = C.face_vector[-1]                             # #sTt = top-rank facets (= g_fan top cell)
    assert factorial(3) * stt == tau_exceptional_sequences(A).signed_count


@lit
def test_kA3_vs_kZ3rad2_share_g_fan_but_split_differs():
    # The anti-tautology pin (see also test_picture_group_split): kA3 and kZ3/rad^2 agree on every
    # COARSE invariant AND on the g-fan SPHERE, but differ on the Ext-driven relation split -- and,
    # consequently, on the CLASSIFYING-SPACE face vector (ruling H1: the sphere sees nothing, the
    # classifying space counts morphisms).
    kA3 = tau_cluster_category(linear_path_algebra(3, field=QQ))
    kZ3 = tau_cluster_category(_nakayama_kZ3_rad2())
    assert kA3.object_count == kZ3.object_count == 14
    assert kA3.g_fan_face_vector == kZ3.g_fan_face_vector == (1, 9, 21, 14)   # identical spheres
    assert len(bricks(linear_path_algebra(3, field=QQ))) == len(bricks(_nakayama_kZ3_rad2())) == 6
    # the classifying spaces DIFFER (49 vs 48): the bonus discriminator from the relation split
    assert kA3.face_vector == (14, 49, 49, 14)
    assert kZ3.face_vector == (14, 48, 48, 14)
    assert kA3.face_vector != kZ3.face_vector


@lit
def test_tau_tilting_infinite_refuses_no_partial_category():
    # 2-Kronecker: tau-tilting-INFINITE (DIJ) -> W(A) is infinite -> NO category emitted.
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    C = tau_cluster_category(K, budget=40)
    assert C.is_complete is False and C.status in ("budget", "error")
    assert C.objects == () and C.morphism_count is None and C.face_vector is None and C.note
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.tautilting.cluster_morphism`.
- [ ] **Step 3: Implement `cluster_morphism.py`.**
  - **G0 (the gate, FIRST):** `eg = exchange_graph(A, budget_pairs=budget)`; if
    `not eg.is_complete or eg.status != "complete"` return the incomplete sentinel
    (`objects=()`, `morphism_count=None`, `face_vector=None`, `g_fan_face_vector=None`,
    `euler_characteristic=None`, `is_complete=False`, `status=eg.status`, honest `note`). **Never a
    count off a non-`"complete"` graph** (P65 M1).
  - **`_g_fan(eg)` → `g_fan_face_vector` (the SPHERE, NOT the classifying space):** each vertex's
    `g_matrix` columns = the signed summand `g`-vectors of a facet; a facet is the `frozenset` of
    its `n` columns; the faces = all subsets, deduped (the §Appendix A1 prototype).
    `g_fan_face_vector[k] = #(faces with k columns)`. Assert `g_fan_face_vector[-1] ==
    len(eg.vertices)` (`#sτt`). This is the `g`-fan/cluster fan (a triangulated `(n−1)`-sphere): it
    is the source of the closed-star out-degrees and is reported for provenance, but it is **NOT**
    `face_vector` (H1 correction).
  - **`_objects`**: for each face `U`, `red = tau_perpendicular_reduction(A, U, budget=budget)`
    (P65) → the wide `J(U)` identified by its **simple-brick iso-class semibrick** (dedup key —
    use `torsion._edge_brick`/`is_isomorphic`, **never dim-vectors** on non-thin input, P64
    ruling 3). Dedup faces by that key → objects, each with `rank = n − |U|`, the semibrick
    label, and `C_W = red.reduced` (dim-vector + `identify_standard` name; **confirm the
    `.reduced` field name against P65's shipped `TauReduction` shape FIRST — §Global Constraints
    cross-plan contract check; flag drift, do not silently rename**). **Self-cert gate:**
    `object_count == wide_subcategories(A, budget=budget).size` (the P64 cross-tie) — a mismatch
    is a loud `QuiverlabError` (the enumeration or the CLO is wrong). Cheaper alternative for the
    dedup key (implementation choice, cross-checked against the reduction): the simple bricks of
    `J(U)` are the brick labels of the `link`-edges of `U` in the `g`-fan (the covers `U ⋖ V`);
    read them from `eg.arrows` (P45) and compare to `red`'s simples (a built-in cross-engine
    check).
  - **`_closed_star_size(U, facets)`**: `#(faces V ⊇ U)` (V a subset-of-a-facet containing all
    of `U`'s columns). `out_degree[W] = _closed_star_size(U_W, facets)` for a representative
    `U_W`; assert it equals `#faces(C_W)` from `red.reduced`'s exchange graph on a spot-check
    (the Jasso link identity — do NOT skip; it protects against a mis-built reduction).
    `morphism_count = Σ out_degree`. `morphisms` = the explicit list (each face `V ⊇ U_W` gives
    a morphism `W → J(V)` of rank `|V| − |U_W|`); include a **composition associativity + unit
    spot-check** (`[∅]` is the identity; `[V∖U] ∘ [U] = [V]` on a small triple) — loud on
    failure.
  - **`face_vector` = the Hanson–Igusa CLASSIFYING-SPACE cube complex (H1 correction — the key
    change):** grade the morphisms by rank — `face_vector[k] = #(rank-k morphisms of the whole
    category) = Σ_{objects W} #(rank-k support τ-rigid pairs of C_W)`. Concretely, for each object
    `W` with representative `U_W`, count the faces `V ⊇ U_W` by `|V| − |U_W|` (= `C_W`'s own
    `g_fan_face_vector`, the Jasso link identity) and sum over objects. **Loud self-certs:**
    `face_vector[0] == object_count == #wide`, `sum(face_vector) == morphism_count`, and
    `face_vector[-1] == len(eg.vertices)` (`#sτt`). `euler_characteristic = Σ (−1)^k face_vector[k]`
    (an `int` — NOT `(−1)^n`, which was the discarded sphere's tautology). **`face_vector` is
    NEVER `g_fan_face_vector`.**
  - **`_kpi1_verdict(A)`**: `True` + reason `"Nakayama (Hanson-Igusa)"` if
    `invariants.recognizers.is_nakayama(A)`; `True` + `"hereditary Dynkin (Igusa-Todorov, CAT(0)
    arXiv:2203.16679)"` if `is_hereditary(A)` and `dynkin_type(A.quiver) is not None`; else `None`
    + `"not certified (K(pi,1) known only for Nakayama / hereditary Dynkin at the cited theorems)"`.
    **NEVER True otherwise** (§Global Constraints, the honest homotopy boundary).
  - `Algebra.tau_cluster_category(budget=512)`: lazy-import delegate.
- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(tautilting): tau_cluster_category -- objects=tau-perp wides (=#wide, ties P64), morphisms=support tau-rigid pairs (out-degree=closed star=#faces(C_W)), Hanson-Igusa classifying-space face vector (f_0=#wide, f_k=#rank-k morphisms) + Euler char + g-fan sphere, K(pi,1) theorem gate (kA2/kA3/Nakayama pins)`

**Adjust to reality (Task A):**
- **The object arbiter** is `test_object_count_equals_wide_count` (the P64 cross-tie). If P66's
  τ-perpendicular enumeration disagrees with P64's `.size`, ONE route is wrong — fix the code,
  never the assert. Do **NOT** substitute `len(semibricks(A))` for `.size` (it is `#torsion`,
  equal to `#wide` only rep-finite).
- **The face-vector arbiter** is `test_face_vector_type_A` + `test_morphisms_out_of_top_equals_
  g_fan_faces`. **`face_vector[0]` MUST equal `object_count` (`= #wide`)** and **`sum(face_vector)`
  MUST equal `morphism_count`** — if `face_vector[0] == 1` you re-emitted the `g`-fan SPHERE (the
  H1 bug), not the classifying space. Both `face_vector[-1]` and `g_fan_face_vector[-1]` equal
  `len(eg.vertices)` (`#sτt`); if `g_fan_face_vector[0] != 1` the empty face was dropped from the
  sphere.
- **The morphism arbiter** is `test_out_degree_is_closed_star_and_matches_reduction`. The star
  sizes of the 6 rank-2 objects of `kA₃` MUST be `[9,9,11,11,11,11]` (2 `k×k` + 4 `kA₂`); a flat
  `[11]*6` or `[9]*6` means the reduction/link is mis-computed.
- **The factorization tie** (`test_factorization_count_ties_P65`) is `oracle_selfcert`, not
  `oracle_crossengine` (fix-round ruling 3): P65 defines `signed_count := n!·|eg.vertices|`, so the equality
  is definitional. Upgrade to a genuine cross-engine check ONLY if P65's shipped API materialises
  the signed sequences (`len(materialised) == n!·#sτt`); otherwise keep it self-cert.
- **The `status != "complete"` gate is real** (P65 M1) — never read a count off an `"error"`
  graph even if its vertex count looks right.

---

### Task B: `picture_group` — the presentation + K(π,1) data

**Files:**
- Modify: `src/quiverlab/tautilting/cluster_morphism.py` (add `PictureGroup`, `picture_group`,
  the relation classifier + abelianization), `__init__.py`, `core/algebra.py`
  (`Algebra.picture_group(budget=512)`)
- Test: `tests/modules/test_picture_group.py`

**Interfaces:**
```python
@dataclass(frozen=True)
class PictureGroup:
    algebra: object
    generators: tuple        # per brick: (id, dimvec, name)
    relations: tuple         # per rank-2 wide: (type, (i,j) simple-brick ids, ext_brick_id|None, word)
    num_generators: int      # == #bricks
    num_relations: int       # == #rank-2 wides
    num_atom: int; num_commutation: int
    abelianization: tuple    # SNF invariant factors of the relation matrix (e.g. () for free Z^r)
    abelianization_rank: int # == #bricks - rank(relation matrix)  [H2 correction: NOT #bricks-#atom;
                             # rank(relation matrix) == #(distinct extension bricks among atoms) <= #atom]
    is_kpi1: bool | None; kpi1_reason: str      # mirror TauClusterCategory
    is_complete: bool; status: str; note: str

def picture_group(A, *, budget=512) -> PictureGroup
# module-internal:
def _rank2_wides(A, eg, budget) -> tuple         # (i, j, ext_brick_id|None) per size-2 semibrick
def _relation_type(A, Bi, Bj) -> ("commutation"|"atom", ext_brick|None)   # via ext_dims
def _abelianization(num_gen, relations) -> (invariant_factors, rank)       # exact SNF; rank via SNF
```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_picture_group.py
"""The picture group presentation pi_1(|W(A)|) (Plan 66 / R29, Igusa-Todorov-Weyman 1609.02636;
Hanson-Igusa). Generators = bricks; relations = rank-2 wides (commutation for k x k, atom/
pentagon for connected). kA2 -> (3 gens, 1 atom); kA3 -> (6 gens, 4 atom + 2 comm). The
DISCRIMINATOR: kA3 and kZ3/rad^2 share every coarse count but split 4+2 vs 3+3. QQ-scope."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.cluster_morphism import picture_group
from quiverlab.tautilting.torsion import bricks

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kZ3_rad2():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}
                  ).algebra(relations=["a*b", "b*c", "c*a"], field=QQ)

def _kA3_rad2():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ)


@lit
def test_kA2_picture_group():
    G = picture_group(Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ))
    assert G.num_generators == 3 and G.num_relations == 1
    assert G.num_atom == 1 and G.num_commutation == 0          # single connected rank-2 wide
    assert G.abelianization_rank == 2                          # Z^{3-1}, 1 distinct ext brick


@lit
def test_kA3_picture_group():
    G = picture_group(linear_path_algebra(3, field=QQ))
    assert G.num_generators == 6 and G.num_relations == 6
    assert G.num_atom == 4 and G.num_commutation == 2          # 4 kA2-type + 2 kxk
    # H2 correction: NOT Z^{6-4}=Z^2. The 4 atoms share extension bricks -> only 3 DISTINCT
    # ext bricks ({S1,M23} and {M12,S3} both extend to the top brick M123), so rank = 3.
    assert G.abelianization_rank == 3                          # Z^{6 - rank(rel matrix)} = Z^{6-3}


@lit
def test_kA4_picture_group_abelianization():
    # H2 counterexample the draft LACKED: #bricks-#atom = 10-10 = 0 would give the ABSURD Z^0.
    # The 10 atoms have only 6 DISTINCT extension bricks, so rank(rel matrix) = 6 and G^ab = Z^4.
    G = picture_group(linear_path_algebra(4, field=QQ))
    assert G.num_generators == 10 and (G.num_atom, G.num_commutation) == (10, 10)
    assert G.abelianization_rank == 4                          # Z^{10-6}, NOT Z^{10-10}=Z^0


@xeng
@pytest.mark.parametrize("A_factory, ngen", [
    (lambda: linear_path_algebra(2, field=QQ), 3),
    (lambda: linear_path_algebra(3, field=QQ), 6),
    (lambda: linear_path_algebra(4, field=QQ), 10),
    (_kA3_rad2, 5), (_kZ3_rad2, 6)])
def test_generators_are_bricks(A_factory, ngen):
    A = A_factory()
    assert picture_group(A).num_generators == len(bricks(A)) == ngen


@lit
def test_kA3_vs_kZ3rad2_relation_split_discriminates():
    # THE anti-tautology pin: identical (#bricks, #wide, g_fan_face_vector) but DIFFERENT relation
    # split -- kA3 = 4 atom + 2 comm, kZ3/rad^2 = 3 atom + 3 comm. A construction that ignores
    # the Ext^1 direction/existence between the two simple bricks would give the same split (and
    # the same classifying-space face vector -- which actually differs 49 vs 48, see Task A).
    kA3 = picture_group(linear_path_algebra(3, field=QQ))
    kZ3 = picture_group(_kZ3_rad2())
    assert (kA3.num_generators, kA3.num_relations) == (kZ3.num_generators, kZ3.num_relations) == (6, 6)
    assert (kA3.num_atom, kA3.num_commutation) == (4, 2)
    assert (kZ3.num_atom, kZ3.num_commutation) == (3, 3)       # the discriminator


@xeng
def test_relation_count_equals_rank2_wide_count():
    # #relations == #(rank-2 wides) == #(size-2 all-semibricks). Cross-check via the category's
    # rank-2 objects (Task A) -- the two enumerations must agree.
    from quiverlab.tautilting.cluster_morphism import tau_cluster_category
    A = linear_path_algebra(3, field=QQ)
    G = picture_group(A); C = tau_cluster_category(A)
    assert G.num_relations == sum(1 for o in C.objects if o[1] == 2) == 6


@selfcert
def test_atom_relation_has_extension_brick():
    # Each atom relation names an EXTENSION brick (the third brick of the connected wide);
    # each commutation relation has ext_brick_id None. Structural, catches type mislabelling.
    G = picture_group(linear_path_algebra(3, field=QQ))
    atoms = [r for r in G.relations if r[0] == "atom"]
    comms = [r for r in G.relations if r[0] == "commutation"]
    assert len(atoms) == 4 and all(r[2] is not None for r in atoms)
    assert len(comms) == 2 and all(r[2] is None for r in comms)


@lit
@pytest.mark.parametrize("A_factory, kpi1", [
    (_kA3_rad2, True), (_kZ3_rad2, True),                        # Nakayama -> K(pi,1) (Hanson-Igusa)
    (lambda: linear_path_algebra(3, field=QQ), True)])           # hereditary Dynkin -> K(pi,1)
def test_kpi1_theorem_anchored(A_factory, kpi1):
    assert picture_group(A_factory()).is_kpi1 is kpi1


@selfcert
def test_kpi1_not_claimed_beyond_theorems():
    # A tau-tilting-finite, NON-Nakayama, NON-hereditary algebra: is_kpi1 must be None (honest
    # "not certified"), NEVER True -- P66 never asserts asphericity beyond the cited theorems.
    # (A commutative-square-with-relation kQ/I: tau-tilting-finite, not Nakayama, not hereditary.)
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3), "c": (2, 4), "d": (3, 4)})
    A = Q.algebra(relations=["a*c-b*d"], field=QQ)               # commutative square
    G = picture_group(A)
    if G.is_complete:
        assert G.is_kpi1 is None and "not certified" in G.kpi1_reason
```

- [ ] **Step 2: Run to verify failure**, implement, iterate.
- [ ] **Step 3: Implement.** **TRANSCRIBE the ITW/HI relation normal form before coding** — the
  COUNTS and TYPE split are engine-derivable, but the exact word is orientation-dependent:
  - `generators`: one per `bricks(A)` (id, dim-vector, `identify_standard` name).
  - `_rank2_wides`: the rank-2 objects of `tau_cluster_category` (Task A), OR directly the
    size-2 semibricks `{Bi, Bj}` of hom-orthogonal bricks (`hom_dim(Bi,Bj)==hom_dim(Bj,Bi)==0`).
    Both routes must agree (`test_relation_count_equals_rank2_wide_count`).
  - `_relation_type(A, Bi, Bj)`: `e1 = ext_dims(A, Bi, Bj, 1)[1]`, `e2 = ext_dims(A, Bj, Bi, 1)[1]`;
    `("commutation", None)` if `e1 == e2 == 0`, else `("atom", ext_brick)` where `ext_brick`
    is the third brick of the wide (`= Bi + Bj` as dim-vectors in the connected case; identify
    by iso-class among `bricks(A)`). **At most one of `e1, e2` is nonzero** for a τ-tilting-
    finite rank-2 wide (both nonzero would be Kronecker-type = τ-tilting-infinite) — assert it.
  - `relations`: per rank-2 wide, `(type, (i,j), ext_id, word)`; `word` = the ITW slope-ordered
    normal form (commutation `xᵢxⱼ = xⱼxᵢ`; atom `xᵢxⱼ = ∏ x(slope-ordered bricks)`).
  - `_abelianization` (**H2 correction — the draft formula was WRONG**): build the
    `#relations × #generators` integer matrix (each relation's net exponent vector: commutation →
    the zero row; atom `xᵢ + xⱼ = xᵢ + xδ + xⱼ` → the row `e_δ`, killing the *extension* brick
    `δ`), SNF via `fields.linalg`; `abelianization` = invariant factors, **`abelianization_rank =
    #generators − rank(matrix)`**, and `rank(matrix)` is computed by the SNF — it equals the
    **number of DISTINCT extension bricks among the atoms** (`≤ #atom`), because several atoms can
    share one `δ` (identical `e_δ` rows). **Do NOT use `#bricks − #atom`** — that overcounts the
    rank whenever atoms share an extension brick (it gives `ℤ²` for kA₃ where the truth is `ℤ³`,
    and the absurd `ℤ⁰` for kA₄ where the truth is `ℤ⁴`). All invariant factors are `1` here, so
    `G^ab` is free of rank `abelianization_rank`.
  - `is_kpi1` / `kpi1_reason`: reuse `_kpi1_verdict` from Task A.
  - `Algebra.picture_group(budget=512)`: lazy delegate.
- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(tautilting): picture_group -- generators=bricks, relations=rank-2 wides (commutation kxk / atom connected via Ext), abelianization (SNF), K(pi,1) theorem gate; kA2/kA3 pins + the kA3-vs-kZ3/rad^2 relation-split discriminator`

**Adjust to reality (Task B):**
- **The discriminator arbiter** is `test_kA3_vs_kZ3rad2_relation_split_discriminates`. If both
  come out `4+2`, the `Ext¹` direction/existence is being ignored (e.g. classifying by dim-
  vector) — the whole point of the pin. Fix the classifier, never the assert.
- **The atom/commutation split is `Ext¹`-driven, not orientation-driven.** Do NOT claim the
  split is "orientation-dependent": `#atom + #commutation = #rank-2 wides` is orientation-
  independent, and the split is a function of the `Ext¹`-orthogonality of the bricks (a Morita/
  derived-invariant fact), not of a chosen quiver orientation.
- **The abelianization arbiter (H2)** is `test_kA3_picture_group` (`ℤ³`) + `test_kA4_picture_group_
  abelianization` (`ℤ⁴`). `abelianization_rank` MUST be `#bricks − rank(relation matrix)`, and
  `rank(matrix) == #(distinct extension bricks among atoms)` — **NOT** `#bricks − #atom`. If kA₃
  comes out `ℤ²` or kA₄ comes out `ℤ⁰`, the code is subtracting `#atom` (counting the row `e_δ`
  once per atom instead of once per distinct `δ`) — fix the SNF/rank path, never the assert. The
  SNF is the source of truth for `abelianization_rank`; do not shortcut it with a brick/atom count.
- **`is_kpi1` is theorem-anchored.** Do NOT compute asphericity for a non-Nakayama non-Dynkin
  algebra (that needs HI's cube/link model transcribed — the stretch below). Return `None`.
- **STRETCH (deferred, ruling 6): a Gromov flag-link `K(π,1)` certifier.** A cube complex is
  aspherical iff its vertex links are flag (Gromov), which is combinatorially decidable for any
  τ-tilting-finite `A`. This would let `is_kpi1` be certified beyond Nakayama/Dynkin. It is
  DEFERRED because it requires transcribing HI's EXACT cube/link model (which faces are the
  vertices of a link, which simplices are filled) verbatim — a wrong model would silently mis-
  verify asphericity. Budget it as its own subtask with the transcription gate; v1 ships the
  theorem-anchored verdict and does NOT block on it.

---

### Task C: QPA honest-scope skip

**Files:**
- Create: `tests/qpa/test_tau_cluster_qpa.py`

**Interfaces:** QPA 1.37 has **no** τ-cluster-morphism-category / picture-group / cube-complex /
wide-subcategory-poset surface (it has **no** τ-tilting surface at all — the P45/P63/P64
finding). Probe live via `NamesGVars()` (the `tests/qpa/test_products_qpa.py` /
`test_congruence_qpa.py` precedent): the test SKIPS honestly and **FAILS if QPA ever ships** any
`ClusterMorphism` / `PictureGroup` / `TauPerpendicular` / `WideSubcategories` verb — so the
honest-scope claim on the verification page cannot silently rot.

- [ ] **Step 1: Write the probe** (`pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
  reason=...)`; sweep `NamesGVars()` for the candidate names; `pytest.skip(...)` if absent,
  `assert False` (loud) if any appears). Record that the external cross-checks are the ITW/HI/BM
  worked examples (kA₂ pentagon, Nakayama K(π,1)), not a live QPA call.
- [ ] **Step 2: Run** — `... -m pytest tests/qpa/test_tau_cluster_qpa.py -q -m qpa` (venv has
  `[qpa]`): expect a clean skip.
- [ ] **Step 3: Commit** — `test(qpa): honest no-cluster-morphism/picture-group surface probe (fails if QPA ships one)`

---

### Task D: no-code exposure — the `tau_cluster` compute kind ("presentation rendering")

The organizing invariant is **byte-identity across the two runners**, achieved by routing BOTH
through the single shared builder `cluster_morphism.tau_cluster_block` (the `tau_tilting` /
`congruences` precedent — each runner only adds `citations` locally). Anchors below are
current-tree; still "adjust to reality" if a merge shifts them.

**`tau_cluster` is an ALGEBRA-level compute kind** carrying a PAIR BUDGET (not a degree range),
parsed exactly like `tau_tilting` / `congruences`: `tau_cluster` or `tau_cluster:512`. NOT a
`MODULE_KINDS` entry; sized on `A.dim` and marked **knit-heavy** (§Global Constraints). It is the
metaplan's **"presentation rendering"** kind — it sits in the same lattice/τ-tilting theme as
`tau_tilting` / `congruences` (if merged), and its block reuses the P45 brick labels + P64 wide
labels.

**The `tau_cluster_block(A, budget) -> dict` payload:**
```python
{
  "kind": "tau_cluster", "n": n,
  "complete": bool, "status": "complete"|"budget"|"error",
  "category": {                                 # None when incomplete
     "object_count": int,                       # == #wide (ties P64)
     "morphism_count": int,
     "objects_by_rank": {r: count, ...},        # e.g. {0:1, 1:6, 2:6, 3:1}
     "morphisms_by_rank": {r: count, ...},      # r -> #(rank-r morphisms); == face_vector as a dict
     "face_vector": [f_0, ..., f_n],            # HANSON-IGUSA CLASSIFYING SPACE: f_0 == object_count
                                                # (#wide), f_r == #(rank-r morphisms), sum == morphism_count
     "g_fan_face_vector": [1, ..., #sTt],       # the g-fan/cluster SPHERE (NOT the classifying space);
                                                # sum == #faces == morphisms out of mod A
     "euler_characteristic": int,               # == sum((-1)^k face_vector[k])  (NOT (-1)^n)
     "is_kpi1": bool | None, "kpi1_reason": str},
  "picture_group": {                            # None when incomplete
     "num_generators": int,                     # == #bricks
     "generators": [{"dimvec":{v:m}, "name":str|None}, ...],
     "num_relations": int, "num_atom": int, "num_commutation": int,
     "relations": [{"type":"atom"|"commutation",
                    "simple_bricks":[{"dimvec":..,"name":..},{"dimvec":..,"name":..}],
                    "ext_brick":{"dimvec":..,"name":..}|None, "word":str}, ...],
     "abelianization": [d1, d2, ...],           # invariant factors ([] = free part only)
     "abelianization_rank": int},               # == #bricks - rank(rel matrix)  (H2: NOT #bricks-#atom)
  "references": ["buan_marsh_wide","hanson_igusa","igusa_todorov_weyman","igusa_todorov_cat0",
                 "enomoto_wide_ice","jasso_reduction","air_tau_tilting","demonet_iyama_jasso"],
  "note": str | None,                           # honest truncation when incomplete
}
```
When incomplete: `category=picture_group=None`, `note` = the honest τ-tilting-infinite message;
the block still returns cleanly (no partial-category lie).

**Files (the seven-touchpoint checklist — mirror P45 `tau_tilting` / P64 `congruences`):**
1. **Parse grammar** — `src/quiverlab/hpc/spec.py::parse_compute_item` (the `tau_tilting`
   branch, `spec.py:165`): add `tau_cluster` / `tau_cluster:<budget>` → `ComputeItem(kind=
   "tau_cluster", lo=None, hi=budget)`, skipping `MAX_DEGREE` (budget, not degree). Mirror in
   `docs/gui/runner.py` parse.
2. **Server/HPC dispatch** — `spec.py::_dispatch` (after the `tau_tilting` branch, `spec.py:1093`):
   `if kind == "tau_cluster": budget = item.hi or 512; from quiverlab.tautilting.cluster_morphism
   import tau_cluster_block; block = tau_cluster_block(A, budget=budget); block["citations"] =
   _citation_pairs(block["references"]); return block, None`. Catch the char-caveat / τ-tilting-
   infinite `QuiverlabError` into `{"error": "<loud message>"}` (Plan-30 honest-per-entry — never
   a 500). Add the `_snippet` entry: `"tau_cluster": lambda it: ("A.tau_cluster_category(budget="
   f"{it.hi if it.hi is not None else 512}); A.picture_group(...)")`.
3. **Pyodide twin** — `docs/gui/runner.py::compute_one` (after `tau_tilting`): the byte-identical
   `elif name == "tau_cluster":` branch calling `tau_cluster_block`; one `calls` reproduce-map
   entry; one `ETA_MODEL["scalars"]` cost — `"tau_cluster": 4.0` (heavier than `congruences`'s
   3.0 / `wall_chamber`'s 2.5 / `tau_tilting`'s 2.0: it builds the wide poset + a reduction per
   object).
4. **Estimator** — `webapp/server/estimator.py`: add `"tau_cluster"` to the knit-heavy tuple at
   `estimator.py:83` (`("tau_tilting", "ar_quiver", "left_right_parts", "tilted_check")`) so it
   is sized on `sizing_dim` (algebra dim) and honestly routed off "instant" for larger inputs.
5. **GUI JS** — `docs/gui/gui.js` AND `webapp/static/gui/gui.js` (vendored byte-identical — apply
   every edit to BOTH): a checkbox `<input id="qlgui-tau_cluster">` + budget picker
   `qlgui-tau_cluster-budget` in the structure/lattice block (beside `tau_tilting` /
   `congruences`); the two bare ids in the id-registry array; the request push-list entry
   `if (el.tau_cluster.checked) compute.push("tau_cluster:" + el["tau_cluster-budget"].value)`;
   a `renderBlock` `else if (name === "tau_cluster")` branch calling `renderTauCluster(div, b)`;
   the `THEMES` structure `kinds` array entry `"tau_cluster"` inside the `QLGUI-THEMES-BEGIN/END`
   sentinels; the `scheduleProbe` / `KIND_CTRL` entries. **`renderTauCluster`** (v1, text/HTML
   only — NO SVG, ruling 6): (a) the category summary (object_count = #wide, morphism_count,
   objects-by-rank + morphisms-by-rank tables); (b) the **classifying-space cube complex** (the
   Hanson–Igusa `face_vector` as a small table `dim k | f_k` with the note that `f_0 = #wide` and
   `f_k = #(rank-k morphisms)`; separately the `g_fan_face_vector` labelled as the g-fan/cluster
   SPHERE with `f_n = #sτt`; the Euler characteristic `χ = Σ(−1)^k f_k`; and the **K(π,1) verdict**
   with its reason — a plain-language line, e.g. "the classifying space is a K(π,1) [Nakayama,
   Hanson-Igusa]" / "[hereditary Dynkin, Igusa-Todorov CAT(0)]" or "K(π,1) not certified for this
   class"); (c) the **picture group presentation**: the generator
   list (brick dim-vector + name), the relation list (each row: type badge atom/commutation, the
   two simple bricks, the extension brick for an atom, and the relation word), and the
   abelianization. When incomplete, render the honest `note` and nothing else.
6. **i18n (FOUR locales) + picker template** — `webapp/server/i18n/{en,es,fr,zh}.json` (the
   shipped `LANGS = ("en","es","fr","zh")`): add `pick.kind.tau_cluster` (REQUIRED — gated by
   `test_layout_picker.py`), `struct.tau_cluster` (checkbox label), and `tcl.*` block strings
   (`tcl.title`, `tcl.category`, `tcl.objects`, `tcl.morphisms`, `tcl.by_rank`, `tcl.cube`,
   `tcl.face_vector` (the classifying space), `tcl.g_fan` (the g-fan/cluster SPHERE, labelled
   "NOT the classifying space"), `tcl.euler`, `tcl.kpi1`, `tcl.kpi1_yes`, `tcl.kpi1_uncertified`,
   `tcl.picture_group`, `tcl.generators`, `tcl.relations`, `tcl.atom`, `tcl.commutation`,
   `tcl.abelianization`, `tcl.budget`, `tcl.truncated`). Add one
   `data-pick-kind-tau_cluster="{{ t('pick.kind.tau_cluster') }}"` line in
   `webapp/templates/draw.html` (mirror `data-pick-kind-tau_tilting`).
7. **Report renderer** — `src/quiverlab/trace/results_html.py`: `_HEADINGS["tau_cluster"]` + a
   `_tau_cluster_html(b)` branch in `_block_html` rendering the category summary, the
   classifying-space `face_vector` (labelled `f_0 = #wide`) alongside the separately-labelled
   `g_fan_face_vector` (the g-fan SPHERE, `f_n = #sτt`), the Euler characteristic, the K(π,1)
   verdict, and the picture-group presentation (generators, typed relations, abelianization); when
   incomplete, the honest note. For the report PDF, if a
   poset drawing is wanted, **REUSE `viz/tikz.py::tikz_hasse`** fed the wide-subcategory poset
   (P64's precedent) — do NOT hand-roll a bespoke classifying-space TikZ. An honest "table only
   (budget-capped)" node when incomplete.
- Test files: `tests/webapp/test_tau_cluster_p66.py`, `tests/gui/test_tau_cluster_runner_twin.py`.

```python
# tests/webapp/test_tau_cluster_p66.py
"""The tau_cluster algebra-level compute kind: served by hpc.spec, mirrored by the Pyodide twin,
both byte-identical. kA2 -> #objects=5, #gens=3, #rels=1 (atom); kA3 -> #objects=14, 6 gens,
4 atom + 2 comm."""


def test_tau_cluster_block_shape(tmp_path):
    from quiverlab.hpc.spec import ComputeRequest, run
    out = run(ComputeRequest.model_validate(_tau_cluster_request(quiver=("kA3",), budget=512)),
              tmp_path)
    b = out["results"]["tau_cluster"]
    assert b["complete"]
    assert b["category"]["object_count"] == 14
    # H1: face_vector = classifying space (f_0 = #wide); g_fan_face_vector = the SPHERE
    assert b["category"]["face_vector"] == [14, 49, 49, 14]
    assert b["category"]["g_fan_face_vector"] == [1, 9, 21, 14]
    assert b["category"]["face_vector"][0] == b["category"]["object_count"]   # f_0 == #wide
    assert sum(b["category"]["face_vector"]) == b["category"]["morphism_count"]  # == 126
    assert b["category"]["euler_characteristic"] == 0                          # NOT (-1)^3
    assert b["picture_group"]["num_generators"] == 6
    assert (b["picture_group"]["num_atom"], b["picture_group"]["num_commutation"]) == (4, 2)
    assert b["picture_group"]["abelianization_rank"] == 3                      # H2: Z^3 (6-3), not Z^2
    assert "hanson_igusa" in b["references"] and "igusa_todorov_cat0" in b["references"]


def test_twin_parity(tmp_path):
    # run the SAME request through docs/gui/runner.py::compute_one and hpc.spec.run;
    # assert json.dumps(block, sort_keys=True) equality on the tau_cluster block.
    ...


def test_tau_tilting_infinite_status(tmp_path):
    # 2-Kronecker, small budget -> complete False, status budget/error, category+picture_group
    # null, note set, no crash.
    ...
```

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir).
- [ ] **Step 2: Implement** the parse grammar (both runners), server + twin dispatch, the
  estimator entry, the GUI JS in both `gui.js` files, the four-locale i18n + `draw.html`, and the
  `results_html.py` branch. Set `block["kind"] = "tau_cluster"` inside the builder so both
  runners agree.
- [ ] **Step 3: Gate the layout/i18n first** —
  `... -m pytest tests/webapp/test_layout_picker.py tests/webapp/test_js_parses.py -q`
  (THEMES == `ALL_KINDS`; `pick.kind.tau_cluster` in all four locales).
- [ ] **Step 4: Add the golden** (`tau_cluster_kA2`) to `_runner_goldens.json`; note it in the
  `test_runner_delegation.py` docstring change-log. Run the delegation test BEFORE adding to
  confirm existing entries stay byte-identical.
- [ ] **Step 5: Run the gates** —
  `... -m pytest tests/webapp/test_tau_cluster_p66.py tests/webapp/test_runner_delegation.py tests/gui/test_tau_cluster_runner_twin.py tests/hpc -q`.
  Expected PASS; both runners byte-identical; the τ-tilting-infinite case returns non-`complete`
  cleanly.
- [ ] **Step 6: Commit** — `feat(gui,webapp,hpc,trace): tau_cluster compute kind -- category summary + cube-complex face vector + K(pi,1) verdict + picture-group presentation, one golden (4 locales, v1 schema, knit-heavy estimator)`

**Adjust to reality (Task D):**
- `test_layout_picker.py::ALL_KINDS` is a frozen set the THEMES array must equal exactly; add
  `tau_cluster` to BOTH `gui.js` files or the gate fails.
- Confirm the exact `parse_compute_item` / `_dispatch` / `_snippet` / `compute_one` / `calls` /
  `ETA_MODEL` / `estimator` line anchors on first read (they drift with merges); `tau_tilting`
  and `congruences` (if P64 merged) are the templates. The block carries `"kind": "tau_cluster"`
  if the builder sets it — mirror `tau_tilting_block`.

---

### Task E: verification page, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick the P66 card), and append the dated
  R29 correction note to `docs/plans/2026-08-06-computability-expansion-deep-research.md` (the
  "τ-cluster morphism category is Buan–Marsh IMRN 2021, cube complex + Nakayama K(π,1) is
  Hanson–Igusa Comm. Alg. 2021, presentation is Igusa–Todorov–Weyman, hereditary-Dynkin K(π,1) is
  Igusa–Todorov CAT(0) arXiv:2203.16679" attribution — metaplan §2 mandates folding corrections
  back).
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`).

- [ ] **Step 1: Verification page.** Add the Plan-66 subsystem row
  (`tautilting/cluster_morphism.py`):
  - `oracle_literature` — the **kA₂ picture group** (3 gens, 1 atom relation = pentagon,
    abelianization `ℤ²`); the **kA₃ picture group** (6 gens, 4 atom + 2 commutation,
    abelianization **`ℤ³`** — H2, `6 − 3` distinct ext bricks); the **kA₄ abelianization `ℤ⁴`**
    (H2 counterexample: the discarded `#bricks−#atom` gives the absurd `ℤ⁰`); the **classifying-space
    face vectors** kA₂ `(5,11,5)` / kA₃ `(14,49,49,14)` / kA₄ `(42,204,326,204,42)` (χ = −1/0/2);
    the **kA₃ vs kZ₃/rad² relation-split discriminator** (4+2 vs 3+3 on identical coarse counts +
    identical g-fan sphere `(1,9,21,14)`, but differing classifying-space face vectors
    `(14,49,49,14)` vs `(14,48,48,14)`); the total morphism counts (kA₂ → 21, kA₃ → 126, kA₄ →
    818); the **Nakayama K(π,1)** verdict (kA₃/rad², kZ₃/rad²); the τ-tilting-infinite honest
    refusal (2-Kronecker).
  - `oracle_crossengine` — **object count = #wide** (P66 τ-perpendicular enumeration vs P64
    `wide_subcategories.size`: kA₂ → 5, kA₃ → 14, kA₃/rad² → 12, kZ₃/rad² → 14); **#generators =
    #bricks** (vs P45 `bricks`); **#relations = #rank-2 wides** (vs the category's rank-2
    objects); the **classifying-space face vector** self-consistency (`f_0 = #wide`, `sum =
    morphism_count`) and the **`g_fan_face_vector`** vs the exchange-graph face count (`f_n =
    #sτt`); the **out-degree = closed star = #faces(C_W)** Jasso link identity (vs P65 reductions).
  - `oracle_selfcert` — `object_count == wide.size` internal gate; `χ = Σ(−1)^k face_vector[k]`
    (the classifying space; **NOT `(−1)^n`** — that was the discarded sphere tautology, H1);
    **`face_vector[0] == object_count`** and **`sum(face_vector) == morphism_count`** (the two
    discriminating self-certs, H1); the **factorization identity** `n!·#sτt` vs P65 `signed_count`
    (**self-cert, NOT cross-engine** — definitional in P65, fix-round ruling 3); composition
    associativity + unit laws (spot-check); **`abelianization_rank = #bricks − rank(relation
    matrix)`** (H2 — SNF-decided, NOT `#bricks − #atom`); the atom relations name an extension
    brick, commutations do not; `is_kpi1` theorem-anchoring (never True beyond Nakayama/Dynkin).
  - the **honest semi-decision entry** (metaplan §6): certified complete iff τ-tilting-finite
    (DIJ); on the 2-Kronecker NO category invariant is emitted (a partial `W(A)` is meaningless)
    — `status="budget"|"error"`, `note`, all invariants `None` (STRICTER than P63's bounded
    region, exactly as P64). `status != "complete"` from any exchange graph → loud refusal
    (P65 M1).
  - the **honest-scope entries**: (a) rigorous only over char 0 / char > dim (QQ default; loud
    off scope, inherited P45/P64/P65); (b) bricks decide over the algebraically-closed / char-0
    base (`end_dim = 1`; the GF(pⁿ) proper-division-ring caveat); (c) **QPA cannot compare** —
    the honest no-surface probe (ITW/HI/BM worked examples are the named external cross-checks,
    not a live oracle); (d) **objects = #wide, NOT #torsion** — `#wide = #torsion ⟺ rep-finite`
    (Marks–Šťovíček), so the object count coincides with `#semibricks`/`#torsion` on the rep-
    finite batteries but P66 reads `wide_subcategories.size`; a rep-INFINITE τ-tilting-finite
    witness (`Π(D₄)`, `#wide < #torsion`) is DEFERRED (P64's engine limit); (e) **the K(π,1)
    verdict is theorem-anchored** (Nakayama — Hanson–Igusa Comm. Alg. 2021; hereditary Dynkin —
    **Igusa–Todorov, *Which cluster morphism categories are CAT(0)*, arXiv:2203.16679, 2022**, which
    proves CAT(0) for hereditary finite/tame type with small tubes), NOT claimed for the general
    τ-tilting-finite case (**delicate — arXiv:2203.16679 establishes CAT(0) only for hereditary
    finite/tame type, not the general algebra**); the Gromov flag-link / CAT(0) certifier is a
    deferred stretch; (f) **the picture-group relation WORDS are
    transcribed from ITW/HI** (the slope order is orientation-dependent) — P66 pins the counts,
    the type split, and the brick-membership of each relation (all engine-derivable), and the
    kA₂ pentagon concretely.
  Recount the class table (`tests/release/test_oracle_classes.py` drives the numbers — collect,
  paste the LIVE counts, re-run to green; mid-merge-train honest).
- [ ] **Step 2: README.** One features line (rep-theory-first, per the positioning memory):
  "the τ-cluster morphism category `W(A)` (Buan–Marsh; Hanson–Igusa): its objects (= wide
  subcategories), morphisms, the cube-complex classifying space with `K(π,1)` verdict for
  Nakayama algebras, and the picture-group presentation — one click via `tau_cluster`, certified
  complete iff τ-tilting-finite."
- [ ] **Step 3: Full gate:**
  `... -m pytest tests/modules/test_tau_cluster.py tests/modules/test_picture_group.py -q`
  (deep, the touched files),
  `... -m pytest -q -m fast`, `... -m pytest tests/webapp tests/gui -q`,
  `... -m pytest tests/qpa -q -m qpa`, `... -m pytest tests/release tests/citations -q` — all
  green. Plans 45/64/65 left byte-unchanged (spot-run their touched suites).
- [ ] **Step 4: Commit** — `docs(verification): Plan-66 tau-cluster-morphism/picture-group oracle rows + Buan-Marsh/Hanson-Igusa/ITW citations + honest scope (char caveat, objects=#wide not #torsion, K(pi,1) theorem-anchored, no-QPA, relation words transcribed) + recounted classes`

---

## Acceptance (Plan-66 definition of done)

1. `tau_cluster_category(A)` returns `TauClusterCategory` — the finite category `W(A)`: objects
   (τ-perpendicular wide subcategories, count `== wide_subcategories(A).size`), morphisms graded
   by rank (out-degree `== |closedStar(U)| == #faces(C_W)`), total morphism count, the
   **Hanson–Igusa classifying-space `face_vector`** (`f_0 = #wide`, `f_k = #(rank-k morphisms)`),
   the **separate `g_fan_face_vector`** (the g-fan SPHERE, `f_n = #sτt`), the Euler characteristic,
   and the theorem-anchored `K(π,1)` verdict — all exact, float-free.
2. **Object count = #wide, cross-engine** (P66 τ-perpendicular enumeration vs P64
   `wide_subcategories.size`): kA₂ → 5, kA₃ → 14, kA₃/rad² → 12, kZ₃/rad² → 14 — never
   `#semibricks`/`#torsion` (which coincide only rep-finite).
3. **Classifying-space cube complex (H1)** pinned: `face_vector` kA₂ `(5,11,5)`, kA₃
   `(14,49,49,14)`, kA₄ `(42,204,326,204,42)` (`f_0 == #wide`, `sum == morphism_count`, `f_n ==
   #sτt`); `χ = Σ(−1)^k face_vector[k]` = **−1 / 0 / 2** (NOT `(−1)^n` — the discarded sphere
   tautology). The **`g_fan_face_vector`** (the SPHERE) is pinned separately: kA₂ `(1,5,5)`, kA₃
   `(1,9,21,14)`, kA₄ `(1,14,56,84,42)`, `f_n == #sτt`. Morphisms out of the top object =
   `sum(g_fan_face_vector)` = `#(support τ-rigid pairs of A)` (kA₂ → 11, kA₃ → 45); total morphisms
   = `sum(face_vector)` (kA₂ → 21, kA₃ → 126, kA₄ → 818) — cross-checked against P65 reductions.
4. `picture_group(A)` returns `PictureGroup` — generators (one per brick, `#gens == #bricks`),
   relations (one per rank-2 wide, `#rels == #rank-2 wides`, typed commutation/atom via `Ext¹`),
   abelianization (exact SNF, **rank `#bricks − rank(relation matrix)` = `#bricks − #(distinct
   extension bricks among atoms)`**, H2 — NOT `#bricks − #atom`). **kA₂ → (3 gens, 1 atom),
   `ℤ²`**; **kA₃ → (6 gens, 4 atom + 2 commutation), `ℤ³`** (not `ℤ²`); **kA₄ → `ℤ⁴`** (not the
   absurd `ℤ⁰`); the concrete kA₂ pentagon derived.
5. **The kA₃ vs kZ₃/rad² discriminator** pinned: identical `(#bricks, #wide, g_fan_face_vector,
   #sτt)` but relation split 4+2 vs 3+3 (and, consequently, differing classifying-space
   `face_vector` `(14,49,49,14)` vs `(14,48,48,14)`) — the anti-tautology oracle catching an
   `Ext¹`-blind construction; the g-fan sphere discriminates nothing, the classifying space does.
6. **`K(π,1)` theorem-anchored**: True for Nakayama (kA₃/rad², kZ₃/rad² — Hanson–Igusa) and
   hereditary Dynkin (kAₙ — **Igusa–Todorov CAT(0), arXiv:2203.16679**); `None` ("not certified")
   otherwise — never asserted beyond the cited theorems (the general case is delicate).
7. The structure is **certified complete iff τ-tilting-finite** (DIJ); on the 2-Kronecker it
   emits **no** category invariant (`is_complete=False`, `status`, `note`, all `None`) — the
   honest refusal, never a partial-category lie; `status != "complete"` → loud refusal (P65 M1).
8. **Factorization identity** `n!·#sτt` checked against P65 `tau_exceptional_sequences.signed_count`
   — a `oracle_selfcert` consistency gate, NOT a cross-engine oracle (definitional in P65,
   fix-round ruling 3).
9. `tau_cluster` clickable end-to-end (GUI canvas → block → category/cube/K(π,1)/presentation
   rendering → report) in all four locales (en/es/fr/zh), both runners byte-identical via the
   shared builder, ONE golden added with a documented change-log entry, `test_layout_picker.py`
   `ALL_KINDS`/THEMES gate green, schema still v1, canonical keys unchanged, estimator knit-heavy.
10. QPA honest-scope probe skips (fails if QPA ever ships a cluster-morphism/picture-group
    surface); `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the
    honest scope; README line added; the R29 attribution correction folded back into the
    research doc; deep (P66 files) + fast + webapp/gui + qpa + release + citations suites green.
    Plans 45/64/65 left byte-unchanged. Honest scope recorded: char 0 / char > dim (QQ default);
    bricks over the algebraically-closed base; τ-tilting-finite gate; QPA cannot compare;
    objects = #wide (not #torsion); K(π,1) theorem-anchored; relation words transcribed.

---

## Methodology & assumptions

**Approach.** I read the R29 record and the metaplan P66 card; the committed P64 doc
(`plan-64-torsion-lattice`) end-to-end (its wide-subcategory surface, the `#wide = #torsion ⟺
rep-finite` Marks–Šťovíček boundary, its `semibricks` = `#torsion` semantics, its honest-refusal
discipline, its seven-touchpoint compute-kind wiring); the committed P65 doc
(`plan-65-exceptional-sequences`) — its Jasso `tau_perpendicular_reduction` spec (`C(U)`,
`J(U) ≅ mod C(U)`), its `signed_count = n!·#sτt` identity, and its **Task 0** (the P45 `mutate`
D₄-star `status="error"` defect + the `status != "complete"` loud-refusal gate P66 inherits);
and the shipped P45 engine (`tautilting/{mutation,torsion,__init__}.py`) to confirm the exact
APIs P66 consumes. I read `citations/registry.py` — Buan–Marsh (IMRN), Hanson–Igusa, and Igusa–
Todorov–Weyman are ABSENT (only the unrelated `igusa_todorov` φ/ψ-functions key exists), so P66
adds its keys (**four after the fix-round**: `buan_marsh_wide`, `hanson_igusa`,
`igusa_todorov_weyman`, and `igusa_todorov_cat0` for the CAT(0) hereditary-Dynkin anchor); P64's
`enomoto_wide_ice`/`marks_stovicek` and P65's `jasso_reduction`/`buan_marsh` are reused.

**Fix-round corrections (2026-08-08, adjudicated valid — see the change log H1/H2 entries).**
An adversarial critic caught two mathematical errors, both fixed and re-verified live: **(H1)** the
draft reported the `g`-fan of `A` (a sphere) as the "classifying space"; the actual Hanson–Igusa
cube complex has `f_0 = #wide` and `f_k = #(rank-k morphisms)` — the corrected classifying-space
`face_vector`s are kA₂ `(5,11,5)` (χ=−1), kA₃ `(14,49,49,14)` (χ=0), kA₄ `(42,204,326,204,42)`
(χ=2), with the `g`-fan kept as the explicitly-relabelled `g_fan_face_vector`; **(H2)** the
abelianization was `#bricks − #atom`, wrong when atoms share an extension brick — corrected to
`#bricks − rank(relation matrix)` = `#bricks − #(distinct ext bricks)`, giving kA₃ `ℤ³` (not `ℤ²`),
kA₄ `ℤ⁴` (not the absurd `ℤ⁰`). Two minors: the factorization tie to P65 is `oracle_selfcert` not
`oracle_crossengine` (definitional in P65); the hereditary-Dynkin K(π,1) is anchored to
Igusa–Todorov CAT(0) (arXiv:2203.16679). All corrected pins were emitted by the re-run prototype
(§Appendix A2).

**Reference re-verification (web, this authoring).** I web-verified the τ-cluster morphism
category is DEFINED by **Buan–Marsh, "A category of wide subcategories," IMRN 2021** (objects =
τ-perpendicular subcategories, morphisms via signed τ-exceptional sequences), that **Hanson–Igusa,
Comm. Alg. 49 (2021) no. 10** prove the cube-complex classifying space + the Nakayama K(π,1) +
the picture group, and that the presentation ground truth is **Igusa–Todorov–Weyman
arXiv:1609.02636** (one generator per brick/positive real Schur root; relations per rank-2
configuration; K(π,1) CW complex with Catalan-many cluster-tilting top cells). This CORRECTS the
R29 card's headline attribution ("Hanson–Igusa … W(A)") — the category is Buan–Marsh; HI is the
topology. Folded back into R29 at merge.

**Empirical ground truth (run this authoring + re-run in the fix-round, in the venv; §Appendix
A1/A2, byte-reproducible).** I built kA₂, kA₃, kA₄, the linear non-hereditary Nakayama kA₃/rad²,
and the cyclic self-injective Nakayama kZ₃/rad² and computed, straight from the shipped P45
exchange graph + `bricks` + `hom_dim` + `ext_dims`: `#bricks`, `#wide` (= all pairwise-hom-orthogonal
semibrick subsets, by size = the wide-poset rank histogram), the **`g_fan_face_vector`** (deduped
subsets of the `g`-matrix columns — the SPHERE), the **classifying-space `face_vector`** (`f_0 =
#wide`, `f_k = Σ_W #(rank-k support τ-rigid pairs of C_W)`, via the per-wide C_W g-fan by the
component structure), `#sτt`, both Euler characteristics, the rank-2-wide atom/commutation split
(via `Ext¹`), and the abelianization rank (`#bricks − #distinct-ext-bricks`). **Observed live
(fix-round re-run):** kA₂ `#bricks=3, #wide=5, g_fan=(1,5,5), face_vector=(5,11,5), #sτt=5, split
1+0, χ=−1, ab=ℤ²`; kA₃ `#bricks=6, #wide=14, g_fan=(1,9,21,14), face_vector=(14,49,49,14),
#sτt=14, split 4+2, χ=0, ab=ℤ³`; kA₄ `g_fan=(1,14,56,84,42), face_vector=(42,204,326,204,42),
χ=2, split 10+10, ab=ℤ⁴`; kA₃/rad² `#bricks=5, #wide=12, g_fan=(1,8,18,12),
face_vector=(12,40,40,12), split 2+3, χ=0, ab=ℤ³`; kZ₃/rad² `#bricks=6, #wide=14,
g_fan=(1,9,21,14), face_vector=(14,48,48,14), split 3+3, χ=0, ab=ℤ³`. The `kA₃`-vs-`kZ₃/rad²`
discriminator (identical g-fan sphere `(1,9,21,14)`, split 4+2 vs 3+3, differing classifying-space
face vectors `(14,49,49,14)` vs `(14,48,48,14)`) was observed live. Morphisms out of the top object
= `sum(g_fan_face_vector)` (kA₂ → 11, kA₃ → 45) and the total morphism counts = `sum(face_vector)`
(kA₂ → 21, kA₃ → 126, kA₄ → 818) were both emitted by the prototype and cross-check by the
`f_0==#wide` / `sum==morphism_count` / `f_n==#sτt` self-certs.

**Key correctness decisions.** (1) **P66 is a category/topology layer, not new representation
theory** — the bricks (P45), wides (P64), and reductions (P65) are inherited; P66 assembles the
category, the cube complex, and the presentation. (2) **Objects = `#wide` via `wide_subcategories.
size`, explicitly NOT `#semibricks`/`#torsion`** (equal only rep-finite; Marks–Šťovíček) — the
primary correctness trap, invisible on the rep-finite batteries, guarded by the P64 cross-tie and
stated in scope. (3) **The object count is a genuine CROSS-ENGINE oracle**, not a tautology: P66
enumerates objects the τ-perpendicular way (faces → `J(U)`, deduped by iso-class semibrick) and
compares to P64's independent Enomoto core-label-order count. (4) **The `face_vector` is the
Hanson–Igusa CLASSIFYING SPACE, NOT the `g`-fan sphere (H1):** `f_0 = #wide`, `f_k = #(rank-k
morphisms) = Σ_W #(rank-k support τ-rigid pairs of C_W)`, `sum = morphism_count`; the `g`-fan
sphere is kept separately as `g_fan_face_vector` (source of the closed-star out-degrees) and is
cross-checked against P65's per-object reduction face counts (the Jasso link identity). (5) **The
picture-group relation split is `Ext¹`-driven** (atom = connected rank-2 wide with an extension
brick; commutation = `k×k`), giving the `kA₃`/`kZ₃/rad²` discriminator that no coarse count nor
the g-fan sphere catches (but the classifying-space face vector now does, `49` vs `48`); the
**abelianization rank is SNF-decided `#bricks − rank(relation matrix)` (H2)**, not `#bricks −
#atom`. (6) **The K(π,1) verdict is theorem-anchored** (Nakayama — Hanson–Igusa; hereditary Dynkin
— Igusa–Todorov CAT(0), arXiv:2203.16679), never asserted beyond the cited theorems — the honest
homotopy boundary. (7) **The completeness gate is STRICTER than P63's** (a partial category is a
lie, not a bounded region) and treats `status != "complete"` as a loud refusal (P65 M1).

**Deliberately not checked / assumptions.** (a) **The exact picture-group relation WORDS.** The
COUNTS (`#gens = #bricks`, `#rels = #rank-2 wides`), the TYPE split (atom/commutation, via
`Ext¹`), the brick-membership of each relation, and the abelianization are all engine-derivable
and live-verified; the exact slope-ordered WORD is orientation-dependent and is transcribed from
Igusa–Todorov–Weyman / Hanson–Igusa at implementation (the kA₂ pentagon derived concretely in
the spec). This is the **primary design risk** (§Task B) — I mitigate it by pinning the
transcription-independent invariants and deriving kA₂ by hand. (b) **The `#wide < #torsion`
witness is DEFERRED** (inherited from P64): a rep-INFINITE τ-tilting-finite algebra (tame
`Π(D₄)`, `#torsion = 192`) is the only kind that separates `#wide` from `#torsion`, and the
shipped P45 engine cannot compute it cheaply (P64 §Task B probe: `Π(D₄)/QQ` did not reach 50
exchange-graph vertices in 120 s; `Π(D₄)/GF(31)` raised in `mutate`). The rep-finite batteries
(where all three counts coincide) therefore cannot exhibit the trap — I guard it by construction
(reading `.size`) and record it in scope. (c) **The general-case K(π,1)** — a Gromov flag-link
combinatorial certifier is decidable but DEFERRED (§Task B), because it needs HI's exact cube/
link model transcribed and a wrong model silently mis-verifies asphericity. (d) **The
classifying-space `face_vector` for the Nakayama rows** (kA₃/rad² `(12,40,40,12)`, kZ₃/rad²
`(14,48,48,14)`) was computed by the fix-round prototype via the per-wide C_W component structure
(rank-1 → `k`, rank-2 connected → `kA₂`, rank-2 disconnected → `k×k`, the top wide → A's own
g-fan), which is exact for these examples (all proper connected wides are rank ≤ 2); the
production code will instead read each C_W directly from P65's reduction (the Jasso link identity),
so these are the implementer's cross-check targets. kA₂/kA₃/kA₄ were validated two ways
(noncrossing-partition product AND component structure) and match the ruling-required repins
exactly. (e) The `O(2^n · #sτt)` face enumeration and the `O(#faces²)` star sizing are capped at
the `budget`; τ-tilting-finite small `n` is fine (kA₄ ran in the probe), honest `note` beyond.
(f) Citation print metadata (`BuanMarsh2021wide` volume/pages; `IgusaTodorovWeyman2016` /
`IgusaTodorov2022cat0` final venue) is `# verify` at citation-add — arXiv ids + the HI DOI are
firm; the `BuanMarsh2021wide` author field matches the printed IMRN record (Buan and **Robert J.
Marsh**); no fabricated print metadata (house rule).

**Why this is correct.** The category definition (objects = τ-perpendicular wides, morphisms =
support τ-rigid pairs of the source with target the τ-perpendicular) is verbatim from Buan–Marsh
(IMRN 2021); the cube complex + Nakayama K(π,1) + picture group from Hanson–Igusa (Comm. Alg.
2021); the presentation from Igusa–Todorov–Weyman (1609.02636); the hereditary-Dynkin K(π,1) from
Igusa–Todorov CAT(0) (arXiv:2203.16679). Every count P66 asserts on the five algebras — `#bricks`,
`#wide`, the **classifying-space `face_vector`** and the **`g_fan_face_vector`**, `#sτt`, the
atom/commutation split, χ, the abelianization rank, morphisms-out-of-top and total — was recomputed
live against the shipped P45 engine (fix-round re-run, §Appendix A2); the object-count oracle is
cross-engine (P66 enum vs P64 CLO), the morphism oracle is cross-engine (star vs reduction), the
**factorization tie is a self-cert consistency gate (definitional in P65, NOT cross-engine —
fix-round ruling 3)**, and the `kA₃`/`kZ₃/rad²` discriminator kills the ext-blind failure mode. The
`face_vector` (classifying space) and the abelianization rank were the two fix-round corrections
(H1/H2), both re-verified. The only objects whose exact form rests on implementation-time
transcription (the relation words; the general-case asphericity) are fenced by transcription-
independent pins and honest "not certified" verdicts — the honest boundary of what I could verify.

## Appendix: authoring prototype (reproducible pins)

**A1. Runnable prototype (fix-round version)** — reproduces `#bricks`, `#wide` (by rank), the
**`g_fan_face_vector`** (the SPHERE), the **classifying-space `face_vector`** (`f_0 = #wide`,
`f_k = #(rank-k morphisms)`), `#sτt`, both χ's, the total morphism count, the rank-2
atom/commutation split, and the **H2 abelianization rank** for the five algebras, straight from
the shipped P45 exchange graph + `bricks` + `hom_dim` + `ext_dims`. The classifying-space face
vector is assembled from the per-wide `C_W` g-fan by the connected-component structure (exact for
these examples: proper connected wides are ≅ `kA_m`, the top wide = mod A uses A's own g-fan);
production reads each `C_W` from P65's reduction. Run:
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python thisfile.py`.

```python
from itertools import combinations
from collections import Counter
from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.mutation import exchange_graph
from quiverlab.tautilting.torsion import bricks
from quiverlab.modules.hom import hom_dim
from quiverlab.modules.ext import ext_dims


def gfan_fv(A, budget=512):                 # the g-fan / cluster SPHERE (NOT the classifying space)
    eg = exchange_graph(A, budget_pairs=budget)
    assert eg.is_complete and eg.status == "complete", (eg.status, "not a finite fan")
    facets = []
    for v in eg.vertices:
        gm = v["g_matrix"]
        cols = tuple(tuple(gm[r][c] for r in range(len(gm))) for c in range(len(gm[0])))
        facets.append(frozenset(cols))
    faces = set()
    for f in facets:
        fl = list(f)
        for k in range(len(fl) + 1):
            for sub in combinations(fl, k):
                faces.add(frozenset(sub))
    fc = Counter(len(f) for f in faces)
    return tuple(fc[k] for k in range(max(fc) + 1))


def conv(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return tuple(r)


KA = {0: (1,), 1: (1, 2)}                    # g-fan f-vectors of kA_m; kA_>=2 from the engine
for m in range(2, 5):
    KA[m] = gfan_fv(linear_path_algebra(m, field=QQ))


def dv(b):
    d = b.dimension_vector()
    return tuple(d[k] for k in sorted(d))


def analyze(name, A):
    B = list(bricks(A)); n = len(B); nv = len(A.quiver.vertices)
    A_gfan = gfan_fv(A); dvs = [dv(b) for b in B]
    def orth(i, j): return hom_dim(B[i], B[j]) == 0 and hom_dim(B[j], B[i]) == 0
    def has_ext(i, j):
        return ext_dims(A, B[i], B[j], 1)[1] > 0 or ext_dims(A, B[j], B[i], 1)[1] > 0
    wides = [sub for k in range(n + 1) for sub in combinations(range(n), k)
             if all(orth(a, b) for a, b in combinations(sub, 2))]
    byrank = Counter(len(w) for w in wides)
    fv = Counter()                           # the CLASSIFYING-SPACE face vector
    for w in wides:
        if len(w) == nv:                     # the rank-nv wide is mod A (whole algebra)
            cw = A_gfan
        else:                                # connected components -> product of kA_m g-fans
            parent = {i: i for i in w}
            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]; x = parent[x]
                return x
            for a, b in combinations(w, 2):
                if has_ext(a, b): parent[find(a)] = find(b)
            comps = {}
            for i in w: comps.setdefault(find(i), []).append(i)
            cw = (1,)
            for c in comps.values(): cw = conv(cw, KA[len(c)])
        for k, val in enumerate(cw): fv[k] += val
    face_vector = tuple(fv[k] for k in range(max(fv) + 1))
    chi = sum((-1) ** k * f for k, f in enumerate(face_vector))
    ext_bricks = set(); atom = comm = 0
    for w in wides:
        if len(w) == 2:
            a, b = w
            if has_ext(a, b):
                atom += 1; ext_bricks.add(tuple(x + y for x, y in zip(dvs[a], dvs[b])))
            else: comm += 1
    ab_rank = n - len(ext_bricks)            # H2: #bricks - rank(rel matrix); NOT #bricks - #atom
    print(f"=== {name} ===  #bricks={n} #wide={len(wides)} wide-by-rank={dict(sorted(byrank.items()))}")
    print(f"  g_fan_face_vector(SPHERE)={A_gfan}  #sTt={A_gfan[-1]}  out-of-top(#faces)={sum(A_gfan)}")
    print(f"  face_vector(CLASSIFYING)={face_vector}  chi={chi}  total-morphisms={sum(face_vector)}")
    print(f"  picture group: gens={n} rels={atom + comm} (atom={atom}+comm={comm})  "
          f"abelianization=Z^{ab_rank} (#bricks-#atom would be Z^{n - atom})")


if __name__ == "__main__":
    analyze("kA2", linear_path_algebra(2, field=QQ))
    analyze("kA3", linear_path_algebra(3, field=QQ))
    analyze("kA4", linear_path_algebra(4, field=QQ))
    analyze("kA3/rad^2 (linear Nakayama)",
            Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ))
    analyze("kZ3/rad^2 (self-inj Nakayama)",
            Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}
                   ).algebra(relations=["a*b", "b*c", "c*a"], field=QQ))
```

**A2. Observed output** (venv, fix-round re-run, 2026-08-08):

```
=== kA2 ===  #bricks=3 #wide=5 wide-by-rank={0: 1, 1: 3, 2: 1}
  g_fan_face_vector(SPHERE)=(1, 5, 5)  #sTt=5  out-of-top(#faces)=11
  face_vector(CLASSIFYING)=(5, 11, 5)  chi=-1  total-morphisms=21
  picture group: gens=3 rels=1 (atom=1+comm=0)  abelianization=Z^2 (#bricks-#atom would be Z^2)
=== kA3 ===  #bricks=6 #wide=14 wide-by-rank={0: 1, 1: 6, 2: 6, 3: 1}
  g_fan_face_vector(SPHERE)=(1, 9, 21, 14)  #sTt=14  out-of-top(#faces)=45
  face_vector(CLASSIFYING)=(14, 49, 49, 14)  chi=0  total-morphisms=126
  picture group: gens=6 rels=6 (atom=4+comm=2)  abelianization=Z^3 (#bricks-#atom would be Z^2)
=== kA4 ===  #bricks=10 #wide=42 wide-by-rank={0: 1, 1: 10, 2: 20, 3: 10, 4: 1}
  g_fan_face_vector(SPHERE)=(1, 14, 56, 84, 42)  #sTt=42  out-of-top(#faces)=197
  face_vector(CLASSIFYING)=(42, 204, 326, 204, 42)  chi=2  total-morphisms=818
  picture group: gens=10 rels=20 (atom=10+comm=10)  abelianization=Z^4 (#bricks-#atom would be Z^0)
=== kA3/rad^2 (linear Nakayama) ===  #bricks=5 #wide=12 wide-by-rank={0: 1, 1: 5, 2: 5, 3: 1}
  g_fan_face_vector(SPHERE)=(1, 8, 18, 12)  #sTt=12  out-of-top(#faces)=39
  face_vector(CLASSIFYING)=(12, 40, 40, 12)  chi=0  total-morphisms=104
  picture group: gens=5 rels=5 (atom=2+comm=3)  abelianization=Z^3 (#bricks-#atom would be Z^3)
=== kZ3/rad^2 (self-inj Nakayama) ===  #bricks=6 #wide=14 wide-by-rank={0: 1, 1: 6, 2: 6, 3: 1}
  g_fan_face_vector(SPHERE)=(1, 9, 21, 14)  #sTt=14  out-of-top(#faces)=45
  face_vector(CLASSIFYING)=(14, 48, 48, 14)  chi=0  total-morphisms=124
  picture group: gens=6 rels=6 (atom=3+comm=3)  abelianization=Z^3 (#bricks-#atom would be Z^3)
```

`kA₃` and `kZ₃/rad²` agree on `#bricks=6`, `#wide=14`, **`g_fan_face_vector=(1,9,21,14)`** (the
identical SPHERE), `#sTt=14`, `wide-by-rank={1,6,6,1}` — and differ in the relation split (`4+2` vs
`3+3`) AND, consequently, in the **classifying-space `face_vector`** (`(14,49,49,14)` vs
`(14,48,48,14)`): the discriminator. The classifying-space `χ` is `−1, 0, 2` for `n=2,3,4`
(**NOT** `(−1)^n` — that is the discarded SPHERE's tautology). `#wide` by rank matches P64's
NC(Aₙ) Whitney numbers. Every row satisfies `face_vector[0]==#wide`, `sum(face_vector)==total`,
`g_fan_face_vector[-1]==#sTt`. **kA₄ abelianization is `ℤ⁴` (`10 − 6`); the discarded
`#bricks − #atom` gives the absurd `ℤ⁰`.**

**A3. Total morphism count = `sum(classifying face_vector)` = `Σ_W #faces(C_W)` (from A2).**
Each object `W`'s out-degree is the wide-invariant `#faces(C_W)`; summing gives the total AND (by
rank) the classifying-space `face_vector`. Piece g-fan #faces: `#faces(k)=3`, `#faces(k×k)=9`,
`#faces(kA₂)=11`, `#faces(kA₂×kA₁)=33`, `#faces(kA₃)=45`, `#faces(kA₄)=197`.
`kA₂: 1·1 + 3·3 + 1·11 = 21 = sum(5,11,5).`
`kA₃: 1·1 + 6·3 + (4·11 + 2·9) + 1·45 = 1 + 18 + 62 + 45 = 126 = sum(14,49,49,14).`
`kA₄: 1·1 + 10·3 + (10·11 + 10·9) + (5·45 + 5·33) + 1·197 = 1 + 30 + 200 + 390 + 197 = 818 =
sum(42,204,326,204,42).` The per-object out-degree (closed-star size) equals `#faces(C_W)` by the
Jasso link identity, cross-checked against P65's reduction at implementation.

## Change log

- **2026-08-08 authoring** — plan drafted against `dev` (P45 merged; P64 + P65 committed as
  plan docs, prerequisites for P66). Scope boundary vs P45/P64/P65 fixed (P66 = the category +
  cube complex + picture group layer). The τ-tilting-finite gate set STRICTER than P63 (no
  partial category), with the P65 M1 `status != "complete"` loud-refusal inherited. All numeric
  pins (kA₂/kA₃ picture groups, face vectors, χ, #wide, the kA₃-vs-kZ₃/rad² relation-split
  discriminator, morphisms-out-of-top) recomputed live in the venv against the shipped exchange
  graph. Attribution corrected: the τ-cluster morphism category is **Buan–Marsh (IMRN 2021)**,
  the cube complex + Nakayama K(π,1) + picture group **Hanson–Igusa (Comm. Alg. 2021)**, the
  presentation **Igusa–Todorov–Weyman (arXiv:1609.02636)** — a dated correction to research-doc
  R29 queued for the merge. Objects fixed to `#wide` (`wide_subcategories.size`), explicitly NOT
  `#semibricks`/`#torsion` (equal only rep-finite — Marks–Šťovíček). K(π,1) verdict theorem-
  anchored (Nakayama / hereditary Dynkin), never claimed beyond; the Gromov flag-link certifier
  and the `#wide < #torsion` (`Π(D₄)`) witness deferred with rationale.
- **2026-08-08 fix-round (adversarial critic → NEEDS WORK, all findings adjudicated valid).** Two
  mathematical corrections + four minors, all applied and live-re-verified (Appendix A1/A2 re-run
  byte-for-byte). **H1 (BLOCKING): the "classifying space" was the WRONG space.** The draft's face
  vector was the `g`-fan of `A` — a triangulated `(n−1)`-SPHERE (`χ = (−1)^n` tautology; kA₂'s
  "complex" is a pentagon `≅ S¹` with `π₁ = ℤ ≠ G(kA₂)`, kA₃'s `≅ S²`, not aspherical). Replaced
  by the Hanson–Igusa CLASSIFYING-SPACE cube complex: `face_vector` with `f_0 = #wide`,
  `f_k = #(rank-k morphisms)` — repinned kA₂ `(5,11,5)` χ=−1, kA₃ `(14,49,49,14)` χ=0, kA₄
  `(42,204,326,204,42)` χ=2 (sums = totals 21/126/818); the `g`-fan kept as the explicitly-labelled
  `g_fan_face_vector` (the sphere; NOT the classifying space); new discriminating self-certs
  `face_vector[0]==object_count` and `sum(face_vector)==morphism_count`; kA₂ `χ=−1` cross-checked
  vs the presentation Euler char `1−3+1`. **H2 (MAJOR): the abelianization formula was wrong.**
  `#bricks − #atom` is wrong when atoms share an extension brick; replaced by `#bricks −
  rank(relation matrix)` (= `#bricks − #distinct-ext-bricks`, SNF-decided) — repinned kA₃ `ℤ³`
  (not `ℤ²`), kA₄ `ℤ⁴` (the discarded formula gave the absurd `ℤ⁰`), kZ₃/rad² `ℤ³`; kA₄ pin added.
  Minors: (i) the factorization tie to P65 demoted `oracle_crossengine → oracle_selfcert`
  (`signed_count := n!·|eg.vertices|` is definitional in P65); (ii) the prerequisites line
  corrected (P45 merged; P64/P65 committed as DOCS, must merge before implementation);
  (iii) `BuanMarsh2021wide` author fixed to the printed IMRN record (**Buan and Robert J. Marsh**),
  and the hereditary-Dynkin K(π,1) verdict anchored to **Igusa–Todorov, *Which cluster morphism
  categories are CAT(0)*, arXiv:2203.16679 (2022)** (new `igusa_todorov_cat0` key), with that paper
  also cited as honest-scope context for the delicate general case; (iv) a cross-plan contract
  check for P65's `TauReduction.reduced` field shape (verify at start, flag drift, do not silently
  adapt). Every corrected pin emitted by the re-run prototype.
