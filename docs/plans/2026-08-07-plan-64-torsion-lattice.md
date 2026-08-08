# Plan 64: Torsion-lattice congruences, canonical joins, core label order = wide subcategories (R26) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote the **lattice theory of torsion classes** of a finite-dimensional algebra
`A = kQ/I` to a first-class, no-code surface, built **on top of the Plan-45 τ-tilting
engine** (which already ships the brick-labelled exchange graph = the Hasse quiver of
`tors A`, plus `bricks`, `semibricks`, `hasse_orientation`, `torsion_class_data`). P64 adds
the four objects the metaplan card names — **all lattice-theoretic, none re-deriving the
representation theory P45 already has**:

1. **The finite lattice `tors A` itself** — extracted from the P45 exchange graph as an
   abstract lattice (elements = functorially-finite torsion classes, order = containment,
   covers = the brick-labelled mutation edges): its join/meet, **join-irreducibles**,
   meet-irreducibles, and the **canonical join representation** of every torsion class
   (DIRRT / Barnard–Carroll–Zhu).
2. **`Con(tors A)`** — the **congruence lattice** of `tors A`, computed via the principal
   cover-congruences (the standard finite-lattice algorithm), distributive by
   Funayama–Nakayama, its join-irreducibles = the forcing-equivalence classes of covers.
3. **The forcing order on bricks** — DIRRT's representation-theoretic reading of the
   join-irreducible congruences: the forcing preorder on covers descends to a **partial
   order on `bricks A`**, and `Con(tors A) ≅` the lattice of order ideals of that brick
   poset.
4. **The core label order ≅ the poset of wide subcategories** — Enomoto's lattice-theoretic
   computation of `wide A` from `tors A` (via the core label order / the extended-κ order),
   turning the P45 brick labels into the **wide-subcategory poset** with no subcategory
   enumeration.

One algebra-level no-code compute kind — **`congruences`** — puts the torsion lattice
summary, the canonical join representations, `Con(tors A)` + the forcing order on bricks, and
the wide-subcategory poset one click away in all four locales. Exact only (the lattice is a
finite combinatorial object; every count is an `int`, no floats in `src/`); the honest
**complete-iff-τ-tilting-finite** contract inherited from P45 (a global lattice invariant
like `Con` or the forcing order is *undefined* on a truncated lattice, so P64 refuses to
emit them when the exchange graph did not close — never a partial-lattice lie).

**What is genuinely NEW over Plan 45 (the scope boundary — read this first).** Plan 45 ships
the *ingredients*: the exchange graph (`ExchangeGraph.vertices`/`.arrows`/`.adj`, the
brick per arrow), `hasse_orientation(eg)` (each edge oriented up/down by torsion-class
size), `bricks(A)` / `semibricks(A)` (iso-class-deduped), `torsion_class_data(pair)` (a
`Gen(M)` fingerprint, injective on pairs). Plan 45 does **not** assemble these into an
*abstract lattice* and computes **no** lattice-theoretic invariant. P64 adds four things P45
does not have:

- **The lattice order and its join/meet** (transitive closure of the oriented Hasse; the
  `is_lattice` self-cert), **join-irreducibles** (a cover of exactly one element),
  meet-irreducibles, and **canonical join representations** (each element's down-covers,
  labelled by its P45 semibrick).
- **`Con(tors A)`** — principal cover-congruences, the distributive congruence lattice, its
  size.
- **The forcing order on `bricks A`** (the poset whose down-set lattice is `Con`), and the
  DIRRT identification **join-irreducibles ↔ bricks**.
- **The wide-subcategory poset via the core label order** (Enomoto) — `wide A` from the
  brick-labelled lattice alone.

**Architecture.** One new module `src/quiverlab/tautilting/congruence.py`, a thin
combinatorial layer over the merged Plan-45 machinery (no new math engines, no new
homology, no subcategory enumeration):

- `congruence.py::TorsionLattice` — a frozen value object: `elements` (torsion-class ids,
  = exchange-graph vertex ids), `order` (the containment relation as reachable-sets),
  `covers` (brick-labelled), `top`/`bottom`, `join_irreducibles`, `meet_irreducibles`,
  `canonical_joins` (per element: its down-covers + their brick labels = a semibrick),
  `is_lattice`, `is_semidistributive`, `is_distributive`, `is_modular`, the honest
  `is_complete`/`status`/`note` contract.
- `congruence.py::CongruenceLattice` — `principal_congruences` (per cover), the distinct
  **join-irreducible congruences** (= forcing classes = bricks), `forcing_order` (the
  partial order on bricks), `size` (`|Con(tors A)|`), `is_distributive` (always True —
  self-cert), the honest contract.
- `congruence.py::WideSubcategoryPoset` — the core-label-order / κ-order model of `wide A`:
  `elements` (core-label sets), `order`, `size` (`#wide`), the two-constructions-agree
  self-cert, the honest contract.
- `torsion_lattice(A, *, budget=512) -> TorsionLattice`,
  `congruence_lattice(A, *, budget=512) -> CongruenceLattice`,
  `wide_subcategories(A, *, budget=512) -> WideSubcategoryPoset`,
  `congruences_block(A, *, budget=512) -> dict` (the shared JSON block for the two runners).
- Thin `Algebra` delegates in `core/algebra.py` beside `Algebra.exchange_graph` (lazy
  import to avoid the `tautilting → core` cycle): `Algebra.torsion_lattice(budget=512)`,
  `Algebra.congruence_lattice(budget=512)`, `Algebra.wide_subcategories(budget=512)`.

**Plan 45 is left byte-unchanged.** `exchange_graph`, `hasse_orientation`, `bricks`,
`semibricks`, `torsion_class_data`, `SupportTauTiltingPair` are consumed **read-only**; no
Plan-45 test moves; no Plan-45 golden re-freezes.

**Tech Stack.** Pure integer/set combinatorics (the lattice is finite; ids are `int`,
orders are `frozenset`/`dict[int, frozenset[int]]`, congruences are partitions =
`frozenset[frozenset[int]]`). Brick labels are iso-classes of `Module`s (reuse P45's
`is_isomorphic` dedup + `identify_standard` naming). No floats in `src/` (AST-gated by
`tests/test_no_floats.py`) — there is nothing to floatify (unlike P63's geometry). All
representation-theoretic content (bricks, semibricks, torsion-class order) is inherited from
Plan 45 over `QQ` by default (the load-bearing char caveat below).

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **P45 (the C4 τ-tilting engine) is a HARD prerequisite, MERGED to `dev`.** This plan
  consumes, read-only: `tautilting.mutation.exchange_graph` / `ExchangeGraph` (the BFS + the
  honest `is_complete`/`status`/`n_regular` contract; `.vertices[i]["pair"]`,
  `.arrows[(i,j)]["brick"]`, `.adj`), `tautilting.torsion.hasse_orientation` (the up/down
  cover orientation by `Gen`-size), `torsion.bricks` / `torsion.semibricks` (iso-class
  deduped), `torsion.torsion_class_data` (the `Gen(M)` fingerprint), `torsion._edge_brick` /
  `torsion._torsion_universe` (the disambiguated edge-brick label, non-thin-safe),
  `modules.hom.is_isomorphic` / `identify_standard`. **If P45 is not merged when a worker
  picks this up, STOP and escalate** — do not fork the exchange-graph or brick machinery.
  (Wave-3 sibling P63 `plan-63-wall-chamber` is independent of P64 — no edge taken; the two
  share only Plan-45 read-only surfaces and citation-merge care, §Task 0.)
- **Buckets are auto-assigned by directory (`tests/conftest.py`): `tests/modules/` → deep.**
  All P64 engine tests live in `tests/modules/` as `tests/modules/test_congruence*.py`
  (deep — they share the exchange-graph budget with the Plan-45 τ-tilting suite).
  `tests/qpa/` → qpa; `tests/webapp/`, `tests/gui/` → fast (extras-gated dirs, unmarked per
  the Plan-32 rule). Run new tests by path during development; finish each task with a
  `-m deep` spot-run of the touched files.
- **THE ENGINE RUNS OVER QQ BY DEFAULT** (the Plan-45 char caveat, inherited verbatim).
  `bricks`, `semibricks`, `is_isomorphic`, `identify_standard`, `torsion_class_data`, and the
  exchange-graph BFS are rigorous only over **char 0 or char > dim** (Dickson/CIW). Over QQ
  the brick enumeration, the torsion-class order, and every count are decisive; the batteries
  pin over **QQ**. Over `char ≤ dim` the underlying engine inherits the loud `QuiverlabError`
  refusal unchanged — P64 never a silent wrong lattice.
- **Bricks decide over the implicit algebraically-closed / char-0 base.** A brick is
  `end_dim(B) == 1` (`End_A(B) = k`); the join-irreducibles ↔ bricks bijection and the
  forcing order on bricks read correctly only when `End_A(B) = k`. Batteries pin over QQ; the
  `GF(pⁿ)` proper-division-ring caveat (`dim_k End(B) > 1`) is stated honestly on the
  verification page (inherited from Plan 45).
- **Honest semi-decision contract (metaplan §1.3; STRICTER than P63's bounded region).**
  Every P64 invariant — the join-irreducibles, `Con(tors A)`, the forcing order, the
  wide-subcategory poset — is a **global** function of the *whole* finite lattice. On a
  budget-truncated exchange graph the lattice is a *prefix*, not a lattice, and a partial
  `Con` / forcing order / `#wide` would be a **lie**, not merely incomplete. So P64 is
  certified **iff `A` is τ-tilting-finite** (`eg.is_complete`), and on an incomplete graph it
  returns `is_complete=False`, `status="budget"`, a `note`, and **omits every lattice
  invariant** (`join_irreducibles=None`, `congruence_size=None`, `forcing_order=None`,
  `wide_size=None`) — it does NOT return a bounded region (there is no meaningful bounded
  congruence lattice). The 2-Kronecker (τ-tilting-infinite, DIJ) is the honest-refusal
  oracle. τ-tilting-finiteness is **decidable** via the Plan-45 BFS (DIJ: τ-tilting-finite ⟺
  finitely many bricks ⟺ the exchange graph closes).
- **`tors A` is a lattice — asserted, never assumed (self-cert).** DIRRT prove `tors A` is a
  complete lattice for every f.d. algebra; when finite it is a finite lattice. P64 verifies
  `is_lattice` (every pair of elements has a unique join and a unique meet under the extracted
  order) and raises a loud `QuiverlabError` on violation — a violation means the exchange
  graph / orientation is wrong, not the mathematics. The oriented exchange graph **is** the
  Hasse quiver of `tors A` (AIR Thm 2.30 / DIRRT): mutation edges are exactly the covers.
- **No floats in `src/`.** The lattice is finite combinatorics; ids are `int`, orders/
  partitions are `frozenset`s, counts are `int`. There is no geometry to floatify. The v1 GUI
  Hasse rendering is a **text/HTML edge-list** (no coordinate math at all); the report PDF reuses
  the shipped float-free `viz/tikz.py::tikz_hasse` (`{p/q}` pgfmath coords). The only place any
  client-side coordinate math would appear is the **deferred** layered-DAG SVG stretch task
  (`docs/gui/gui.js`, exempt) — not in v1.
- **Composition is left-to-right** (`a*b` = first `a` then `b`). All refusals are
  `QuiverlabError(message, hint=...)` (the two-space `[hint: ...]` convention, `errors.py`).
  Presentation-less (structure-constant) algebras: the exchange graph already refuses loudly
  (the quiver is required) — P64 surfaces that refusal as a clean typed 4xx, never a 500.
  Monotone/theorem-gate violations (`Con` not distributive, `#J(Con) ≠ #bricks`, join-irred
  count ≠ #bricks) raise `QuiverlabError` (loud internal error), never a silently-wrong
  payload.
- **Every citation key is BibTeX-verified before use** (Plan-29 rule); `_r(...)` registry
  entry + `references.bib` entry both land, gated by `tests/citations/test_bib_structure.py`.
- New algebra-level scalar kind acts on the **existing algebra block** ⇒ **schema stays v1**
  (the `tau_tilting`/`wall_chamber` precedent). **Canonical keys are request-derived**
  (`cache.py::canonical_key`, sha256 over the sorted-keys blob + `library_version()`): a
  request that does not use the new `congruences` key keys byte-identically to before, so
  existing goldens are untouched. Verify with `test_runner_delegation.py` BEFORE and AFTER
  adding the golden.
- **Plan-32 markers** (detailed per task): certificate/identity tests (`is_lattice`,
  `Con` distributive, `#J(Con) == #bricks`, join-irred == #bricks, forcing-poset down-sets ==
  `|Con|`, two wide constructions agree) = `oracle_selfcert`; the **kA₂ pentagon** (`tors kA₂`
  = N5: 5 elements, semidistributive, NOT modular, NOT distributive), `|Con(tors kA₂)| = 5`,
  the **kA₃ congruence lattice** (`|tors kA₃| = 14`, `|Con(tors kA₃)| = 14`), the
  wide-subcategory counts (kA₂ → 5 ≅ M₃, kA₃ → 14 ≅ NC), the τ-tilting-infinite honest refusal
  = `oracle_literature`; **join-irreducibles ↔ bricks** (lattice vs P45 `bricks`), **canonical
  joins ↔ semibricks** (lattice down-covers vs P45 `semibricks`), **#J(Con) = #bricks** (the
  forcing order lives on bricks) = `oracle_crossengine`; the QPA honest-scope probe = the
  `qpa` bucket. `oracle_*` markers FORBIDDEN in `tests/{webapp,gui,hpc}`
  (`test_oracle_classes.py`) — engine batteries in `tests/modules/`, cross-runner twins in
  `tests/webapp/`/`tests/gui/` unmarked.
- **Mid-merge-train counts.** v1.0.0 lands ~29 subplans in overlapping waves. **Task E
  recounts the oracle-class table at merge time** by running
  `tests/release/test_oracle_classes.py` (paste the LIVE numbers, never a guessed count) and
  claims only this plan's deltas.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table
  green) and the README line. Conventional commits; green at every commit; branch
  `plan-64-torsion-lattice` off `dev`.

---

## Record (R26, verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R26 — Torsion-lattice congruences, canonical joins, core label order = wide
> subcategory poset.** [D-scout P8, reframed: brick LABELS already ship (P45)]
> New content only: Con(tors A), the forcing order on bricks, canonical join
> representations, core label order ≅ wide subcategories. Refs:
> Demonet–Iyama–Reading–Reiten–Thomas arXiv:1711.01785 (Trans. AMS B);
> Barnard–Carroll–Zhu. Oracles: join-irreducibles ↔ bricks (count cross-engine);
> kA₂ pentagon; kA₃ congruence lattice. Size S–M.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P64 (Wave 3, tier γ,
**independent**; "P66 needs P64" — the τ-cluster morphism category ties R26). Branch
`plan-64-torsion-lattice` off `dev`.

---

## Reference re-verification (mandatory standing rule; done at spec time — with a CORRECTION)

Re-verified via web + the arXiv/journal abstracts of the three primary sources during
authoring. **Findings (transcribe the `# PIN`s at citation-add time):**

1. **Demonet–Iyama–Reading–Reiten–Thomas, "Lattice theory of torsion classes: Beyond
   τ-tilting theory."** VERIFIED: published in **Transactions of the American Mathematical
   Society, Series B, vol. 10 (2023), pp. 542–612**, DOI **10.1090/btran/100**,
   arXiv:1711.01785. Abstract (verbatim, the mandate): *"…the set of torsion classes … is a
   complete lattice which enjoys very strong properties, as bialgebraicity and complete
   semidistributivity. … we introduce the brick labelling of its Hasse quiver and use it to
   study lattice congruences of `tors A`. In particular, we give a representation-theoretical
   interpretation of the so-called forcing order, and we prove that `tors A` is completely
   congruence uniform."* ⟵ **THE congruence + forcing + semidistributivity mandate.** The
   metaplan's "Trans. AMS B" is confirmed exactly.

2. **Barnard–Carroll–Zhu, "Minimal inclusions of torsion classes."** VERIFIED:
   **Algebraic Combinatorics, vol. 2, no. 5 (2019), pp. 879–901**, DOI **10.5802/alco.72**.
   Abstract (verbatim): *"…we characterize the cover relations in `torsΛ` by certain
   indecomposable modules. … First, we show that the completely join-irreducible torsion
   classes (torsion classes which cover precisely one element) are in bijection with bricks.
   Second, we characterize faces of the canonical join complex of `torsΛ` in terms of
   representation theory. … we study the torsion theory of a quotient of the preprojective
   algebra of type Aₙ. We show that its torsion class lattice is isomorphic to the weak order
   on Aₙ."* ⟵ **THE join-irreducibles ↔ bricks + canonical-join-complex mandate.**

3. **CORRECTION — "core label order ≅ wide subcategories" is NOT in Barnard–Carroll–Zhu.**
   The metaplan card cites "Barnard–Carroll–Zhu" for the core-label-order = wide-subcategory
   result, but BCZ (Algebraic Combinatorics 2 (2019), §2) proves **join-irreducibles ↔
   bricks** and the **canonical join complex** — it does **not** contain the "core label
   order" nor the wide-subcategory poset isomorphism. The precise result is
   **Enomoto, "From the lattice of torsion classes to the posets of wide subcategories and
   ICE-closed subcategories," arXiv:2201.00595** (published in **Algebras and
   Representation Theory**, DOI **10.1007/s10468-023-10214-0**; `# verify` the volume/year/
   pages at citation-add). Abstract (verbatim): *"…we compute the posets of wide
   subcategories and ICE-closed subcategories from the lattice of torsion classes … by using
   the kappa map in a completely semidistributive lattice. … for a completely
   semidistributive lattice, we give two poset structures on the set of elements with
   canonical join representations: the **kappa order** (defined using the extended kappa map
   of Barnard–Todorov–Zhu), and the **core label order** (generalizing the shard intersection
   order for congruence-uniform lattices). Then we show that these posets for the lattice of
   torsion classes coincide and are isomorphic to the poset of wide subcategories."* ⟵ **THE
   core-label-order ≅ wide-subcategories mandate.** **P64 cites Enomoto (2201.00595) for the
   wide/CLO result, BCZ (alco.72) for join-irreducibles↔bricks + canonical joins, and DIRRT
   (btran/100) for the lattice/forcing/congruence structure.** A dated correction note is
   appended to the research doc R26 at merge (metaplan §2: corrections folded back).

4. **Barnard–Todorov–Zhu (the extended κ map)** underpins Enomoto's κ-order. Enomoto cites
   it as "Barnard–Todorov–Zhu"; the paper is *"Dynamical combinatorics and torsion classes"*
   (J. Pure Appl. Algebra 225 (2021) — `# verify` title/venue/volume at citation-add). P64
   **may** add it as a secondary key on the κ-order construction; if it cannot be
   BibTeX-verified at merge, the κ-order is cited through **Enomoto** (which defines the
   extended κ map it uses) — no fabricated metadata (house rule).

5. **Con of a finite lattice — the standard theory P64 leans on (textbook, no new citation
   required, but recorded for the implementer).** *Funayama–Nakayama:* the congruence lattice
   `Con(L)` of any lattice is **distributive**. *Grätzer (finite-lattice theory):* the
   **join-irreducible congruences** of a finite lattice are exactly the principal congruences
   `con(a, b)` of prime intervals (covers `a ⋖ b`); by **Birkhoff**, `Con(L) ≅` the lattice
   of order ideals (down-sets) of the poset of join-irreducible congruences. DIRRT's
   contribution is the **representation-theoretic identification** of that poset with
   `(bricks A, forcing order)`. (These are stated in the DIRRT paper as the framework; cite
   `dirrt_lattice_torsion` for the whole package.)

6. **Marks–Šťovíček, "Torsion classes, wide subcategories and localisations."** VERIFIED
   (web, fix-round 2026-08-08): **Bulletin of the London Mathematical Society, vol. 49, no. 3
   (2017), pp. 405–416**, DOI **10.1112/blms.12033**, arXiv:1503.04639. The Ingalls–Thomas maps
   between torsion classes and wide subcategories give an **injection `wide A ↪ tors A`**
   (`#wide ≤ #torsion`), and the two are in **bijection iff `A` is representation-FINITE** — the
   boundary is **representation-finiteness, NOT heredity** (a prior draft's "coincides only for
   hereditary `A`" was wrong; every rep-finite algebra, hereditary or not, has `#wide = #torsion`).
   This is why the `#wide < #torsion` discriminator requires a rep-INFINITE τ-tilting-finite
   algebra (§Task B, ruling 1). Cite `marks_stovicek` for the bound and the honest-scope
   statement. **Preprojective rep-type** (also fix-round, for the deferred witness): `Π(Δ)` is
   representation-finite iff `Δ = Aₙ, n ≤ 4`, tame iff `Δ ∈ {A₅, D₄}`, wild otherwise
   (Dlab–Ringel; Jan Schröer, *FD-Atlas of finite-dimensional algebras*) — so `Π(D₄)` and
   `Π(A₅)` are tame ⇒ representation-infinite, the guaranteed strict-inequality witnesses.

---

## Mathematical spec — the four objects P64 computes

Throughout, `A` is **τ-tilting-finite** (else the honest refusal, §Global Constraints), so
`tors A` is a **finite** lattice and every torsion class is functorially finite (DIJ). Write
`L = tors A`. The **oriented exchange graph is the Hasse quiver of `L`** (AIR / DIRRT):
elements = support τ-tilting pairs `(M, P) ↔ Gen(M)`, covers = mutation edges, each cover
`T ⋗ T'` labelled by the **brick** `B(T, T')` (DIRRT brick labelling; the P45
`_edge_brick`, iso-class-disambiguated so non-thin algebras never collapse two bricks of the
same dim-vector).

### 1. The lattice `tors A` — order, join/meet, join-irreducibles, canonical joins

- **Order.** `T' ≤ T ⟺` there is a downward directed path `T ⇝ T'` in the oriented Hasse
  quiver (reflexive-transitive closure of `hasse_orientation`). Equivalently `Gen(M') ⊆
  Gen(M)`; the two agree because the exchange graph is the Hasse quiver. P64 computes the
  order by Warshall closure over the ≤ `budget` vertices.
- **Join/meet.** `T ∨ U =` the unique minimal common upper bound; `T ∧ U =` the unique
  maximal common lower bound. The **`is_lattice` self-cert** asserts uniqueness for every
  pair (raises `QuiverlabError` otherwise). Top = `(A, 0)` (= `mod A`), bottom = `(0, A)` (=
  `0`).
- **Join-irreducibles.** `T` is **join-irreducible** ⟺ `T` covers **exactly one** element
  (`T` has a unique lower cover `T_*`), i.e. `T ≠ bottom` and out-degree-in-the-downward-Hasse
  `= 1`. **Barnard–Carroll–Zhu:** the completely-join-irreducible torsion classes are **in
  bijection with `bricks A`** — the unique down-cover `T_* ⋖ T` is labelled by a brick, and
  `T ↦ B(T_*, T)` is the bijection `J(tors A) → bricks A`. (`# PIN`: kA₂ → 3, kA₃ → 6;
  live-verified, §Verified pins.) Meet-irreducibles dually (unique upper cover).
- **Canonical join representation.** In a (completely) semidistributive lattice every element
  `T` has a unique irredundant minimal join representation `T = ⋁ D(T)`, the **canonical join
  representation**; the joinands `D(T)` are join-irreducibles and their brick labels form a
  **semibrick** (pairwise Hom-orthogonal bricks). Combinatorially, `D(T) =` the set of
  **down-covers** of `T`, and the label set `= {B(T', T) : T ⋗ T'}` is **exactly the P45
  `semibricks` entry for `T`** (the Asai down-labelling). This gives the cross-engine oracle
  **canonical joins ↔ semibricks** (§Task B): the lattice down-cover count per element and the
  P45 semibrick agree, and `#semibricks = |L|` (each element one canonical join rep).

### 2. `Con(tors A)` — the congruence lattice

A **lattice congruence** `θ` is an equivalence relation on `L` compatible with `∨` and `∧`.
`Con(L)` (the set of all congruences, ordered by refinement) is a **distributive** lattice
(Funayama–Nakayama). P64 computes it by the standard finite-lattice algorithm:

- **Principal cover-congruence `con(T', T)`** for each cover `T' ⋖ T`: the smallest congruence
  collapsing `T' ≡ T`, computed by union–find closed under the compatibility rule *(if
  `x ≡ y` then `x∨z ≡ y∨z` and `x∧z ≡ y∧z` for all `z`)*, iterated to a fixed point.
- **Join-irreducible congruences = distinct `con(T', T)`** (Grätzer: exactly the
  prime-interval principal congruences; each is join-irreducible in `Con`). Deduped by
  partition equality → `J(Con(L))`.
- **`Con(L) ≅` down-sets of `J(Con(L))`** ordered by partition refinement (Birkhoff). P64
  reports `|Con(tors A)|` = the number of order ideals of that poset, and (small cases) the
  poset itself. (`# PIN`: `|Con(tors kA₂)| = 5`, `|Con(tors kA₃)| = 14`; live-verified.)
- **Self-cert:** `Con(L)` is distributive **by construction** (it is a down-set lattice), and
  `|J(Con(L))| == #bricks(A)` (the forcing order lives on bricks, next item) — both asserted.

### 3. The forcing order on bricks

The **forcing preorder** on covers: `(T'⋖T)` **forces** `(U'⋖U)` ⟺ every congruence
collapsing `T'≡T` also collapses `U'≡U` ⟺ `con(U',U) ⊆ con(T',T)`. DIRRT's
representation-theoretic reading: the forcing-equivalence classes of covers are **exactly the
bricks** (all covers with the same brick label force each other), so the forcing preorder
descends to a **partial order on `bricks A`**, and

> **`Con(tors A) ≅` the lattice of order ideals of `(bricks A, forcing order)`** (DIRRT).

P64 builds the forcing order as the partition-**refinement** order on the distinct
`con(T', T)` (= `J(Con)`): `θ_i ≤ θ_j` iff `θ_i ⊆ θ_j` as relations (every block of `θ_i` sits
inside a block of `θ_j`; `θ_i` finer), transported to bricks via the BCZ bijection
`J(tors A) → bricks A` (the join-irreducible torsion class of each cover's brick). **Direction
(definition-consistent — ruling 4).** With the forcing definition above, `(T'⋖T)` forces
`(U'⋖U)` ⟺ `con(U',U) ⊆ con(T',T)`, i.e. the *forced* cover has the *smaller* (contained,
finer) principal congruence — so in the forcing ORDER the forced brick sits **below** the
bricks that force it (`con(U',U) ⊆ con(T',T)` reads `θ(U'U) ≤ θ(T'T)` in `Con`). `# PIN`
(live-verified this authoring, appendix prototype): kA₂ → 3-brick poset with **2 relations,
`{θ_c ≤ θ_a, θ_c ≤ θ_b}`** — the single brick `θ_c` (whose `con` is contained in the other
two, so it is *forced by* them) sits **below** the other two: a "V" opening **upward** (one
minimum, two maxima); its 5 down-sets = `|Con| = 5`. kA₃ → 6-brick poset with **9 relations**
(14 down-sets = `|Con| = 14`).

### 4. The core label order ≅ the poset of wide subcategories

Enomoto computes `wide A` from `L` purely lattice-theoretically, via **two** poset structures
on the elements-with-canonical-join-representations (= all of `L` in the finite case) that
**coincide**:

- the **κ order** — the extended κ map (Barnard–Todorov–Zhu) sends each join-irreducible `j`
  to the meet-irreducible `κ(j)` (`κ(j) =` the largest `m ≥ j_*` with `m ⊉ j`), extended to
  canonical join representations;
- the **core label order (CLO)** — generalizing Reading's shard intersection order for
  congruence-uniform lattices, ordering the "core label sets" of intervals `[0̂, w]`.

Enomoto's theorem: **both posets, for `L = tors A`, are isomorphic to `(wide A, ⊆)`**, the
poset of wide subcategories of `mod A` ordered by inclusion. P64 computes **both** and
asserts they agree (`oracle_selfcert` — Enomoto's own theorem), reports `#wide` = the number
of core-label sets, and the poset. **Implementation MUST transcribe Enomoto's exact κ-map and
core-label definitions** (§Task B, and the design-risk flag in Methodology — my authoring
prototype reproduced the *counts* via a naive chain-label-intersection but got the *poset*
wrong, so the naive route is explicitly rejected). `# PIN`: kA₂ → `#wide = 5`, poset ≅ **M₃**
(bottom `0`, three incomparable atoms `add S₁`, `add S₂`, `add P₁`, top `mod A`); kA₃ →
`#wide = 14`, poset ≅ the noncrossing-partition lattice `NC(A₃)`. **Honest note (CORRECTED —
ruling 1).** The map `tors A → wide A` is a *bijection* **iff `A` is representation-FINITE** —
NOT merely iff `A` is hereditary (**Marks–Šťovíček**, arXiv:1503.04639, Bull. LMS 49 (2017):
"over representation-finite algebras torsion classes are in bijection with wide subcategories";
`wide A ↪ tors A` is always injective, so `#wide ≤ #torsion`, with **equality ⟺ `A` is
representation-finite**). Consequently EVERY `kAₙ` and *every* rep-finite algebra — including the
radical-square-zero Nakayama `1→2→3` and `kZ₂/rad²` that a prior draft named as candidates — has
`#wide = #torsion`, so **no rep-finite algebra can ever witness `#wide < #torsion`**. On the `kAₙ`
the count is additionally `Catalan(n+1)` (`5`, `14`), so the count alone does NOT discriminate;
the **poset structure** does (`M₃ ≇ N₅ = tors kA₂`, `NC(A₃) ≇ Tamari(A₃) = tors kA₃`) — this is
P64's **primary** wide-vs-torsion discriminator (§Task B, ruling 2). A genuine `#wide < #torsion`
witness needs a **representation-INFINITE but τ-tilting-finite** algebra (then Marks–Šťovíček
*forces* the strict inequality); the sourced candidates are the **tame** (hence rep-infinite)
preprojective algebras `Π(D₄)` (`#torsion = |W(D₄)| = 192`) and `Π(A₅)` (`|W(A₅)| = 720`) — but
the shipped P45 engine cannot compute either cheaply at spec time (probe evidence in §Task B), so
this strict-inequality oracle is **DEFERRED to implementation**, not shipped as a required
spec-time pin. Preprojective rep-type sourced: `Π(Δ)` is representation-finite iff `Δ = Aₙ, n ≤ 4`,
tame iff `Δ ∈ {A₅, D₄}`, wild otherwise (Dlab–Ringel; Jan Schröer, *FD-Atlas*).

---

## Verified pins — LIVE-computed against the shipped P45 engine (this authoring)

Recomputed in the venv (`exchange_graph` + `hasse_orientation` + `bricks` + `semibricks`,
then the P64 lattice extraction / principal-congruence / down-set-count prototype). **Every
number below was observed live** except the two labelled `# to-verify` (the wide/CLO poset
structure — count observed, structure literature-anchored):

| algebra | `|tors A|` | lattice properties | `#bricks` = `#J(tors A)` | `#semibricks` = canonical joins | `|Con(tors A)|` | `#J(Con)` = forcing classes | `#wide` (CLO) |
|---------|-----------|--------------------|--------------------------|------------------------------------|-----------------|-----------------------------|---------------|
| **kA₂** (`1→2`) | **5** (pentagon **N₅**) | semidistributive ✓, modular ✗, distributive ✗ | **3** ✓ | **5** ✓ | **5** ✓ | **3** (= #bricks) ✓ | **5** (≅ M₃ `# to-verify struct.`) |
| **kA₃** (`1→2→3`) | **14** | semidistributive ✓, modular ✗, distributive ✗ | **6** ✓ | **14** ✓ | **14** ✓ | **6** (= #bricks) ✓ | **14** (≅ NC(A₃) `# to-verify struct.`) |

**kA₂ is the headline "pentagon" pin (DIRRT / the R26 card).** `tors(kA₂)` is the 5-element
lattice **N₅** — the two maximal chains from `0` to `mod A` have lengths 3 and 2; it is
**semidistributive** (DIRRT: every `tors A` is) but **not modular** and **not distributive**
(N₅ is the smallest non-modular lattice, and it IS semidistributive — the two facts are
independent, which is exactly why `tors A` can be N₅). Its congruence lattice has **5**
elements; its **3** join-irreducibles are in bijection with the **3** bricks `S₁, S₂, P₁`; the
forcing order on those 3 bricks is a "V" opening **upward** — one brick sits **below** the other
two (that one's principal congruence is *contained in* the other two, so it is **forced by**
them, not the reverse; ruling 4) — whose 5 order ideals give `|Con| = 5`. This is a **full
lattice-theoretic discrimination on one fixture**. Pinned verbatim in
`tests/modules/test_congruence_lattice.py`.

**kA₃ is the "congruence lattice" pin (the R26 card).** `|tors kA₃| = 14`, `|Con(tors kA₃)| =
14`, `6` join-irreducibles = `6` bricks, forcing order a 6-brick poset with 9 relations. The
`14 = Catalan(4)` count of `tors` is Plan-45-pinned; the `|Con| = 14` and the `6`-brick
forcing poset are the NEW P64 pins.

**Why the `#wide` column equals `|tors A|` here (not a bug, and not a discriminator).** Both
`kA₂` and `kA₃` are **representation-finite**, so `#wide = #torsion` by Marks–Šťovíček (ruling
1) — the equality is *forced*, so these rows cannot witness `#wide < #torsion`. The wide-poset
*structure* (`M₃`, `NC(A₃)`) is what P64 pins here, and it IS discriminating (`M₃ ≇ N₅`,
`NC(A₃) ≇ Tamari(A₃)`; ruling 2). A strict-inequality witness would need a rep-INFINITE
τ-tilting-finite algebra (`Π(D₄)`/`Π(A₅)`); see §Task B for why that oracle is deferred. The two
`# to-verify struct.` wide-poset structures were NOT observed live at authoring (the `wide`
computation is exactly what P64 implements); they are literature-anchored (Ingalls–Thomas:
`wide(kAₙ) ≅ NC(Aₙ)`).

---

### Task 0: citations (do first — later tasks reference the keys)

**Files:**
- Modify: `src/quiverlab/citations/references.bib`, `src/quiverlab/citations/registry.py`
- Test: `tests/citations/test_bib_structure.py` (existing gate — `bibtex(key)` resolvable,
  registry ↔ bib in sync).

**Interfaces:** the `_r(key, bibtex_key, kind, title, annotation, *tags)` registry helper
(`registry.py:24`). **Reuse (do NOT re-add):** `air_tau_tilting` (→ `AIR2014`),
`demonet_iyama_jasso` (→ `DIJ2019`), `king_stability` (→ `King1994`) — all shipped by Plan 45
(`registry.py:493`–`508`). `_citation_pairs` (`hpc/spec.py`) resolves keys → `[key,
formatted]` pairs for payloads.

**Cross-plan citation-merge care (metaplan §5).** Wave-3 sibling P63 adds
`brustle_smith_treffinger`, `asai_semibricks`, `kaipel_treffinger`. If P63 merges FIRST, they
are present — do NOT re-add. P64's four new keys (`dirrt_lattice_torsion`,
`barnard_carroll_zhu`, `enomoto_wide_ice`, `marks_stovicek`) are disjoint from P63's; on a
rebase, merge `references.bib` / `registry.py` by inserting the P64 block, never `git add -A`,
and re-run `test_bib_structure.py`. (`marks_stovicek` was added by the fix-round: the
rep-finite `#wide = #torsion` boundary it establishes is now load-bearing for the honest-scope
statement — ruling 1.)

- [ ] **Step 1: add the new keys** (BibTeX-verify each before committing — Plan-29 rule):

```bibtex
@article{DIRRT2023,
  author  = {Demonet, Laurent and Iyama, Osamu and Reading, Nathan and Reiten, Idun and Thomas, Hugh},
  title   = {Lattice theory of torsion classes: Beyond {$\tau$}-tilting theory},
  journal = {Transactions of the American Mathematical Society, Series B},
  volume  = {10},
  pages   = {542--612},
  year    = {2023},
  doi     = {10.1090/btran/100},
}
@article{BCZ2019,
  author  = {Barnard, Emily and Carroll, Andrew and Zhu, Shijie},
  title   = {Minimal inclusions of torsion classes},
  journal = {Algebraic Combinatorics},
  volume  = {2},
  number  = {5},
  pages   = {879--901},
  year    = {2019},
  doi     = {10.5802/alco.72},
}
@article{Enomoto2023wide,
  author  = {Enomoto, Haruhisa},
  title   = {From the lattice of torsion classes to the posets of wide subcategories and {ICE}-closed subcategories},
  journal = {Algebras and Representation Theory},
  year    = {2023},
  doi     = {10.1007/s10468-023-10214-0},
  note    = {arXiv:2201.00595},
}
@article{MarksStovicek2017,
  author  = {Marks, Frederik and {\v{S}}{\v{t}}ov{\'\i}{\v{c}}ek, Jan},
  title   = {Torsion classes, wide subcategories and localisations},
  journal = {Bulletin of the London Mathematical Society},
  volume  = {49},
  number  = {3},
  pages   = {405--416},
  year    = {2017},
  doi     = {10.1112/blms.12033},
}
```

  and in `registry.py`:

```python
_r("dirrt_lattice_torsion", "DIRRT2023", "foundation",
   "Lattice theory of torsion classes: Beyond tau-tilting theory",
   "Demonet-Iyama-Reading-Reiten-Thomas: tors A is a complete, bialgebraic, completely "
   "semidistributive, completely congruence-uniform lattice; the brick labelling of its "
   "Hasse quiver; the representation-theoretic forcing order and the congruence lattice "
   "Con(tors A). The Plan-64 congruence + forcing ground truth.",
   "tau-tilting", "lattice"),
_r("barnard_carroll_zhu", "BCZ2019", "foundation",
   "Minimal inclusions of torsion classes",
   "Barnard-Carroll-Zhu: cover relations of tors A characterized by indecomposables; the "
   "completely join-irreducible torsion classes are in bijection with bricks; faces of the "
   "canonical join complex read representation-theoretically. Plan 64's join-irreducibles "
   "<-> bricks and canonical join representations.",
   "tau-tilting", "lattice"),
_r("enomoto_wide_ice", "Enomoto2023wide", "foundation",
   "From the lattice of torsion classes to the posets of wide subcategories and ICE-closed "
   "subcategories",
   "Enomoto: the kappa order (extended kappa map of Barnard-Todorov-Zhu) and the core label "
   "order on a completely semidistributive lattice coincide and are isomorphic to the poset "
   "of wide subcategories. Plan 64's core-label-order = wide-subcategory computation.",
   "tau-tilting", "lattice"),
_r("marks_stovicek", "MarksStovicek2017", "foundation",
   "Torsion classes, wide subcategories and localisations",
   "Marks-Stovicek: the Ingalls-Thomas maps between torsion classes and wide subcategories; "
   "wide A injects into tors A, and the two are in BIJECTION iff A is REPRESENTATION-FINITE "
   "(not merely hereditary). Plan 64's #wide <= #torsion bound and the honest-scope statement "
   "that the count discriminates only off the representation-finite case.",
   "tau-tilting", "lattice"),
```

  **Spec-ambiguity resolution (recorded):** BibTeX keys `DIRRT2023` / `BCZ2019` /
  `Enomoto2023wide` / `MarksStovicek2017`; snake-case registry keys `dirrt_lattice_torsion` /
  `barnard_carroll_zhu` / `enomoto_wide_ice` / `marks_stovicek` (house convention). `DIRRT2023`,
  `BCZ2019`, and `MarksStovicek2017` are fully BibTeX-verified (venue/volume/pages/DOI confirmed —
  `MarksStovicek2017` = Bull. LMS **49** (2017), no. 3, **405–416**, DOI **10.1112/blms.12033**,
  arXiv:1503.04639, web-verified at fix-round). `Enomoto2023wide` ships with the arXiv id in
  `note` and the AlgRepThy DOI; `# verify` the AlgRepThy volume/year/pages at merge (upgrade the
  entry if the print volume has appeared). The **Barnard–Todorov–Zhu** κ-map paper is an OPTIONAL secondary key
  (`# PIN` at citation-add; add `barnard_todorov_zhu` ONLY if BibTeX-verifiable, else the
  κ-order is cited through `enomoto_wide_ice`, which defines the extended κ map — no
  fabricated metadata).
- [ ] **Step 2:** `... -m pytest tests/citations/test_bib_structure.py -q` green.
- [ ] **Step 3: Commit** — `docs(citations): torsion-lattice references (DIRRT btran/100, Barnard-Carroll-Zhu alco.72, Enomoto wide/ICE 2201.00595, Marks-Stovicek blms.12033)`

---

### Task A: `congruence.py` — the lattice, join-irreducibles, canonical joins, `Con`, forcing

**Files:**
- Create: `src/quiverlab/tautilting/congruence.py`
- Modify: `src/quiverlab/tautilting/__init__.py` (export `TorsionLattice`,
  `CongruenceLattice`, `torsion_lattice`, `congruence_lattice` — `wide_subcategories` /
  `WideSubcategoryPoset` added in Task B), `src/quiverlab/core/algebra.py` (thin lazy
  delegates beside `Algebra.exchange_graph`)
- Test: `tests/modules/test_torsion_lattice.py`, `tests/modules/test_congruence_lattice.py`

**Interfaces:**
- Consumes: `tautilting.mutation.exchange_graph`, `tautilting.torsion.hasse_orientation` /
  `bricks` / `semibricks` / `_edge_brick` / `_torsion_universe`, `modules.hom.is_isomorphic`
  / `identify_standard`.
- Produces:
  ```python
  @dataclass(frozen=True)
  class TorsionLattice:
      algebra: object
      elements: tuple              # exchange-graph vertex ids 0..N-1 (= torsion classes)
      labels: dict                 # id -> human label (P45 _label; e.g. "(P1(+)S1,{2})")
      covers: tuple                # ((upper_id, lower_id, {"brick_dimvec":{v:m}, "brick_name":str|None}), ...)
      order: dict                  # id -> frozenset(ids <= it)   (the reflexive-transitive closure)
      top: int; bottom: int
      join_irreducibles: tuple     # ids with a unique lower cover
      meet_irreducibles: tuple     # ids with a unique upper cover
      canonical_joins: dict        # id -> tuple of (lower_id, brick_dimvec, brick_name)  (down-covers)
      is_lattice: bool
      is_semidistributive: bool; is_distributive: bool; is_modular: bool
      is_complete: bool; status: str; note: str

  @dataclass(frozen=True)
  class CongruenceLattice:
      algebra: object
      join_irreducible_congruences: tuple   # the distinct principal cover-congruences (partitions)
      forcing_order: tuple                  # ((brick_i_id, brick_j_id), ...) meaning brick_i <= brick_j
      forcing_bricks: tuple                 # brick dim-vector + name per J(Con) class (the BCZ transport)
      size: int                             # |Con(tors A)| = #order-ideals of the forcing poset
      is_distributive: bool                 # always True (self-cert)
      is_complete: bool; status: str; note: str

  def torsion_lattice(A, *, budget=512) -> TorsionLattice
  def congruence_lattice(A, *, budget=512) -> CongruenceLattice
  # helpers (module-internal):
  def _extract(A, budget) -> (eg, n, covers_down, covers_up, edge_brick, geq, leq) | incomplete
  def _lattice_ops(n, geq, leq) -> (is_lattice, join, meet)
  def _principal_congruence(cover, n, join, meet) -> frozenset[frozenset[int]]
  def _order_ideals_count(poset_pairs, m) -> int
  ```
- `Algebra.torsion_lattice(budget=512)` / `Algebra.congruence_lattice(budget=512)`:
  lazy-import delegates.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_torsion_lattice.py
"""The finite lattice tors A extracted from the P45 exchange graph (Plan 64 / R26).
Literature: tors(kA2) = the pentagon N5 (semidistributive, NOT modular, NOT distributive);
tors(kA3) has 14 elements. Cross-engine: #join-irreducibles == #bricks (BCZ); each element's
canonical-join down-covers carry exactly its P45 semibrick. Self-cert: is_lattice (unique
join/meet); (A,0) the unique top, (0,A) the unique bottom. QQ-scope (the P45 char caveat)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.congruence import torsion_lattice
from quiverlab.tautilting.torsion import bricks, semibricks

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@lit
def test_kA2_is_the_pentagon_N5():
    L = torsion_lattice(_kA2())
    assert L.is_complete and len(L.elements) == 5
    assert L.is_lattice and L.is_semidistributive          # DIRRT: every tors A is SD
    assert not L.is_modular and not L.is_distributive      # N5 is non-modular but SD
    assert len(L.join_irreducibles) == 3                   # = #bricks
    # the two maximal chains of N5 have lengths 3 and 2 (one long side, one short side)
    # (structural check via cover multiset -- 5 covers on 5 elements)
    assert len(L.covers) == 5


@lit
@pytest.mark.parametrize("n, size", [(2, 5), (3, 14)])     # Catalan(n+1)
def test_lattice_size(n, size):
    L = torsion_lattice(linear_path_algebra(n, field=QQ))
    assert L.is_complete and len(L.elements) == size and L.is_lattice


@xeng
@pytest.mark.parametrize("n, nbricks", [(2, 3), (3, 6)])
def test_join_irreducibles_biject_bricks(n, nbricks):
    # BCZ: completely join-irreducible torsion classes <-> bricks. Count cross-engine:
    # #join-irreducibles (lattice) == #bricks (P45 surface).
    A = linear_path_algebra(n, field=QQ)
    L = torsion_lattice(A)
    assert len(L.join_irreducibles) == len(bricks(A)) == nbricks


@xeng
def test_canonical_joins_are_semibricks():
    # STRUCTURAL (ruling 3): each element's canonical-join down-covers carry a set of pairwise
    # Hom-orthogonal bricks; the LABEL SETS THEMSELVES -- not merely their count -- must equal
    # the P45 semibricks. kA2 is thin, so a brick's dim-vector fingerprints it; compare the two
    # collections as MULTISETS of dim-vector label-sets (per-element, not just cardinality).
    A = _kA2()
    L = torsion_lattice(A)
    assert len(L.canonical_joins) == len(L.elements)       # one canonical join per element
    assert len(semibricks(A)) == len(L.elements)           # #semibricks == |L| (5 == 5)

    def _dv(d):
        return tuple(sorted((int(k), int(v)) for k, v in dict(d).items()))

    # canonical_joins[e] = tuple of (lower_id, brick_dimvec, brick_name)
    lattice_labelsets = sorted(
        tuple(sorted(_dv(dv) for (_lo, dv, _nm) in joinands))
        for joinands in L.canonical_joins.values())
    p45_labelsets = sorted(
        tuple(sorted(_dv(B.dimension_vector()) for B in sb))
        for sb in semibricks(A))
    assert lattice_labelsets == p45_labelsets              # per-element label sets AGREE
    # SCOPE: on a NON-thin algebra (kZ2/rad^2) the dim-vector does NOT fingerprint a brick;
    # the rigorous comparison is by MODULE iso-class (P45 torsion._same_iso_multiset -- the
    # same dedup semibricks() itself uses). The builder must therefore carry iso-class
    # identity, not just dim-vector+name, so this oracle stays exact off the thin case.


@selfcert
def test_unique_top_and_bottom():
    L = torsion_lattice(linear_path_algebra(3, field=QQ))
    # bottom = (0,A) covers nothing below; top = (A,0) initial pair, covered by nothing above
    assert L.order[L.bottom] == frozenset({L.bottom})
    assert all(L.bottom in L.order[e] for e in L.elements)  # bottom below everything
    assert all(L.top in L.order_above(e) for e in L.elements) if hasattr(L, "order_above") else True


@lit
def test_tau_tilting_infinite_refuses_no_partial_lattice():
    # the 2-Kronecker is tau-tilting-infinite (DIJ): NO lattice is emitted (a partial lattice
    # would make Con/forcing meaningless) -- honest is_complete False, status budget.
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    L = torsion_lattice(K, budget=40)
    assert L.is_complete is False and L.status == "budget"
    assert L.join_irreducibles == () and L.note
```

```python
# tests/modules/test_congruence_lattice.py
"""Con(tors A) + the forcing order on bricks (Plan 64 / R26). Literature: |Con(tors kA2)| = 5
(the pentagon), |Con(tors kA3)| = 14. Self-cert: Con is distributive (Funayama-Nakayama);
#join-irreducible congruences == #bricks (the forcing order lives on bricks, DIRRT).
Cross-engine: the forcing-poset order-ideal count == |Con|."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.congruence import congruence_lattice
from quiverlab.tautilting.torsion import bricks

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("n, con_size", [(2, 5), (3, 14)])
def test_congruence_lattice_size(n, con_size):
    C = congruence_lattice(linear_path_algebra(n, field=QQ))
    assert C.is_complete and C.size == con_size


@selfcert
def test_con_is_distributive():
    # Con(L) is ALWAYS distributive (Funayama-Nakayama); P64 builds it as a down-set lattice,
    # so distributivity is structural -- assert the flag and (small case) verify directly.
    C = congruence_lattice(linear_path_algebra(2, field=QQ))
    assert C.is_distributive is True


@xeng
@pytest.mark.parametrize("n, nbricks", [(2, 3), (3, 6)])
def test_join_irreducible_congruences_are_bricks(n, nbricks):
    # DIRRT: the forcing-equivalence classes of covers are the bricks; #J(Con) == #bricks.
    A = linear_path_algebra(n, field=QQ)
    C = congruence_lattice(A)
    assert len(C.join_irreducible_congruences) == len(bricks(A)) == nbricks
    assert len(C.forcing_bricks) == nbricks


@selfcert
def test_kA2_forcing_order_is_a_V():
    # the 3 bricks of kA2 with the forcing order: ONE brick BELOW the other two -- a "V"
    # opening UPWARD -- whose 5 order ideals give |Con| = 5. Ruling 4: pin the DIRECTION,
    # not just the relation count. forcing_order pairs are (i, j) meaning brick_i <= brick_j.
    C = congruence_lattice(Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ))
    assert len(C.forcing_bricks) == 3
    assert len(C.forcing_order) == 2                   # exactly two strict relations
    lowers = {i for (i, j) in C.forcing_order}         # the "below" (forced) endpoints
    uppers = {j for (i, j) in C.forcing_order}         # the "above" (forcing) endpoints
    assert len(lowers) == 1 and len(uppers) == 2       # ONE minimum below TWO maxima
    assert lowers.isdisjoint(uppers)                   # the minimum is not itself an upper
    # the minimum brick is the one whose principal congruence is CONTAINED in the other two
    # (con(min) subset con(a), con(b)) -- it is FORCED BY them, not the reverse. Observed live
    # this authoring: relations [(2, 0), (2, 1)] (brick 2 below bricks 0, 1; appendix prototype).
    assert C.size == 5
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.tautilting.congruence`.
- [ ] **Step 3: Implement `congruence.py`.**
  - **`_extract(A, budget)`**: `eg = exchange_graph(A, budget_pairs=budget)`; if
    `not eg.is_complete` return the incomplete sentinel (drives the honest refusal). Else
    `orient = hasse_orientation(eg)`; build `covers_down[a]`/`covers_up[b]` from each edge
    `(i,j)` (a covers b, a above b); record `edge_brick[(a,b)] =
    torsion._edge_brick(A, eg.vertices[a]["pair"], eg.vertices[b]["pair"],
    eg.arrows[(i,j)]["brick"], _torsion_universe(A, budget))` (the iso-class label, non-thin
    safe) + its dim-vector + `identify_standard` name. Warshall-close to `geq` (`id -> set of
    ids <= it`) and `leq`.
  - **`_lattice_ops`**: for every `(x, y)` the meet = the unique maximal common lower bound
    (`z` in `geq[x] & geq[y]` with every common lower bound `≤ z`), join dually; `is_lattice`
    = every pair has exactly one of each; **raise `QuiverlabError` if not** (a P45/orientation
    bug, not the mathematics — DIRRT guarantee `tors A` is a lattice).
  - **join/meet-irreducibles**: `join_irreducibles = [x for x in elements if
    len(covers_down[x]) == 1]`; meet dually. **`canonical_joins[x]` = the down-covers of `x`**
    with their brick labels (the Asai semibrick of `x`).
  - **Lattice properties** (small, `O(N³)`): `is_distributive` /`is_modular` /
    `is_semidistributive` via the standard identities over all triples (SD∨: `x∨y == x∨z ⟹
    x∨y == x∨(y∧z)`; SD∧ dually; distributive: `x∧(y∨z) == (x∧y)∨(x∧z)`; modular: `x ≤ z ⟹
    x∨(y∧z) == (x∨y)∧z`). Cap the triple sweep at the `budget` element count; note the
    `O(N³)` cost in the docstring (fine for τ-tilting-finite small `N`; the sweep is skipped
    for very large `N` with `is_semidistributive=None` + a note — but DIRRT guarantee SD, so
    the check is a self-cert, never load-bearing for correctness).
  - **`congruence_lattice`**: `_principal_congruence(cover)` for each cover via union–find
    closed under `∨`/`∧` compatibility to a fixed point (the standard algorithm; see the
    authoring prototype in Methodology). Dedup partitions → `J(Con)`. `forcing_order` = the
    partition-refinement inclusion among them (`θ_i ≤ θ_j` iff every block of `θ_i` is inside
    a block of `θ_j`). Transport each `J(Con)` class to its **brick** via the BCZ bijection
    (the cover's `_edge_brick`; all covers in a forcing class share the brick). `size =
    _order_ideals_count(forcing_order, |J(Con)|)` (brute-force down-set enumeration; capped —
    if `|J(Con)|` exceeds a safe bound, report `size=None` + a note, never hang). `is_distributive
    = True` (structural; the small-case test also verifies the distributive identity directly).
  - **Self-cert gates** (loud `QuiverlabError`): `len(J(Con)) == len(bricks(A))` (forcing
    lives on bricks) and `len(join_irreducibles) == len(bricks(A))` (BCZ) — a mismatch means
    the extraction is wrong.
  - `Algebra.torsion_lattice` / `Algebra.congruence_lattice`: lazy-import delegates beside
    `Algebra.exchange_graph`.
- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(tautilting): torsion_lattice + congruence_lattice -- lattice extraction, join-irreducibles=bricks, canonical joins=semibricks, Con via principal cover-congruences, forcing order on bricks (kA2 pentagon / kA3 |Con|=14 pins)`

**Adjust to reality (Task A):**
- **The order-extraction arbiter** is `test_kA2_is_the_pentagon_N5` (5 elements, SD, not
  modular/distributive, 3 join-irreducibles, 5 covers) and `test_join_irreducibles_biject_bricks`.
  If join-irreducibles come out `≠ #bricks`, the cover-orientation is wrong (down-covers must
  be counted in the *downward* Hasse; a join-irreducible has **one lower cover**, i.e.
  out-degree 1 in `covers_down`).
- **The principal-congruence fixed point** must iterate until no union changes (the authoring
  prototype confirms convergence on kA₂/kA₃); the arbiter is `test_congruence_lattice_size`
  (kA₂ → 5, kA₃ → 14) and `test_join_irreducible_congruences_are_bricks`. If `|Con|` comes
  out `6` on kA₂, a compatibility propagation step was dropped (the `x∧z`/`x∨z` closure must
  run over ALL currently-merged pairs and ALL `z`, to a fixed point — the authoring prototype
  is the reference; the `c∧b = 0` cross-forcing in N₅ is the subtle step).
- **The `is_lattice` guard** is real: if a pair has two maximal common lower bounds the
  orientation/closure is wrong. Never suppress it — it protects against a mis-extracted Hasse.
- The `order_above` helper referenced in `test_unique_top_and_bottom` is optional; if not
  added, drop that clause (the `bottom`/`order` clauses are the arbiter). Keep the test honest.

---

### Task B: `wide_subcategories` — the core label order ≅ the wide-subcategory poset

**Files:**
- Modify: `src/quiverlab/tautilting/congruence.py` (add `WideSubcategoryPoset`,
  `wide_subcategories`, the κ map + core-label constructions), `__init__.py`,
  `core/algebra.py` (`Algebra.wide_subcategories(budget=512)`)
- Test: `tests/modules/test_wide_subcategories.py`

**Interfaces:**
```python
@dataclass(frozen=True)
class WideSubcategoryPoset:
    algebra: object
    elements: tuple          # the distinct core-label sets (one per wide subcategory)
    order: dict              # inclusion order among them
    labels: tuple            # per element: the brick-label set (its "simple objects")
    size: int                # #wide
    constructions_agree: bool  # kappa order == core label order (Enomoto; self-cert)
    is_complete: bool; status: str; note: str

def wide_subcategories(A, *, budget=512) -> WideSubcategoryPoset
# module-internal:
def _kappa_map(L: TorsionLattice) -> dict         # extended kappa (Barnard-Todorov-Zhu)
def _kappa_order(L) -> (elements, order)
def _core_label_order(L) -> (elements, order)     # Enomoto / shard-intersection generalization
```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_wide_subcategories.py
"""The wide-subcategory poset via the core label order (Plan 64 / R26, Enomoto 2201.00595).
Literature: #wide(kA2) = 5 (poset M3: three incomparable atoms = the 3 bricks); #wide(kA3) =
14 (poset NC(A3)). Self-cert: the kappa order and the core label order coincide (Enomoto).
DISCRIMINATOR (ruling 1/2): #wide == #torsion iff A is representation-FINITE (Marks-Stovicek),
so on the rep-finite kA_n the COUNT never discriminates -- the poset STRUCTURE (M3, NC(A3))
does. The strict-inequality #wide < #torsion oracle needs a rep-INFINITE tau-tilting-finite
algebra and is DEFERRED to implementation (see the comment block at the end). QQ-scope."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.congruence import torsion_lattice, wide_subcategories

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("n, nwide", [(2, 5), (3, 14)])   # kA_n rep-finite => #wide == #torsion
def test_wide_count_type_A(n, nwide):
    # NOTE: on kA_n the count is Catalan(n+1) AND equals #torsion, because kA_n is
    # representation-FINITE (Marks-Stovicek: #wide == #torsion iff rep-finite). The count is
    # therefore NON-discriminating here; the poset STRUCTURE (below) is the arbiter.
    W = wide_subcategories(linear_path_algebra(n, field=QQ))
    assert W.is_complete and W.size == nwide


@selfcert
def test_kappa_order_equals_core_label_order():
    # Enomoto: the two poset constructions coincide (and both = wide subcategories).
    W = wide_subcategories(linear_path_algebra(2, field=QQ))
    assert W.constructions_agree is True


# --- poset helpers: W.order[e] = frozenset of elements <= e (down-set incl. e), matching the
#     TorsionLattice.order convention (id -> frozenset(ids <= it)). --------------------------
def _strict_below(W, e):
    return set(W.order[e]) - {e}

def _bottom(W):
    b = [e for e in W.elements if _strict_below(W, e) == set()]
    assert len(b) == 1, ("unique bottom", b)
    return b[0]

def _top(W):
    t = [e for e in W.elements if set(W.order[e]) == set(W.elements)]
    assert len(t) == 1, ("unique top", t)
    return t[0]

def _atoms(W, bot):
    return [e for e in W.elements if _strict_below(W, e) == {bot}]

def _lower_covers(W, x):
    below = _strict_below(W, x)
    return [c for c in below if not any(c in _strict_below(W, d) for d in below if d != c)]

def _incomparable(W, a, b):
    return a not in W.order[b] and b not in W.order[a]


@lit
def test_kA2_wide_poset_is_M3():
    # M3 (== NC(A2)): bottom 0, THREE pairwise-incomparable atoms (add S1, add S2, add P1),
    # top mod A; each atom covers bottom and is a lower cover of top. Ruling 2: assert the
    # SHAPE, not just the count -- and actually USE the derived bottom/top/atoms.
    W = wide_subcategories(Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ))
    assert W.size == 5 and len(W.elements) == 5 and len(W.labels) == 5
    bot, top = _bottom(W), _top(W)
    atoms = _atoms(W, bot)
    assert len(atoms) == 3                                   # exactly three atoms
    for i in range(3):                                       # pairwise incomparable
        for j in range(i + 1, 3):
            assert _incomparable(W, atoms[i], atoms[j])
    for a in atoms:                                          # each atom covers bottom
        assert _strict_below(W, a) == {bot}
    assert set(_lower_covers(W, top)) == set(atoms)          # top covers EXACTLY the 3 atoms
    # rank (= |label set| = #simple objects of the wide subcat) histogram of M3 is 1,3,1
    assert sorted(len(lab) for lab in W.labels) == [0, 1, 1, 1, 2]
    assert sum(1 for lab in W.labels if len(lab) == 1) == 3  # 3 singleton-brick atoms


@lit
def test_kA3_wide_poset_is_NC_A3():
    # NC(A3) == NC of [4]: 14 elements; rank (= |label set|) Whitney numbers 1,6,6,1 (Narayana
    # N(4,k)); 6 atoms each covering bottom, 6 coatoms each covered by top; the poset is
    # self-dual. Ruling 2: pin the STRUCTURE beyond the count 14.
    from collections import Counter
    W = wide_subcategories(linear_path_algebra(3, field=QQ))
    assert W.is_complete and W.size == 14 and len(W.elements) == 14
    bot, top = _bottom(W), _top(W)
    assert dict(Counter(len(lab) for lab in W.labels)) == {0: 1, 1: 6, 2: 6, 3: 1}
    assert len(_atoms(W, bot)) == 6                          # 6 atoms (single-brick wides)
    assert len(_lower_covers(W, top)) == 6                   # 6 coatoms (top's lower covers)
    # self-duality of NC(4): down-set-size and up-set-size histograms coincide
    downs = Counter(len(W.order[e]) for e in W.elements)
    ups = Counter(sum(1 for f in W.elements if e in W.order[f]) for e in W.elements)
    assert downs == ups
    # RIGOROUS OPTION: an order-isomorphism check of W.order against the hardcoded NC(4)
    # covering table DERIVED in the plan appendix (14 partitions of {1,2,3,4}, 28 covers) is
    # the definitive arbiter; add it at implementation if this Whitney/atom/coatom/self-dual
    # battery is ever insufficient to distinguish NC(4) from another 14-element graded poset.


# --- DEFERRED (ruling 1, path b): the strict-inequality #wide < #torsion discriminator. -----
# A rep-FINITE algebra can NEVER witness #wide < #torsion (Marks-Stovicek: #wide == #torsion
# iff rep-finite; wide ↪ tors is always injective). So the prior draft's candidates
# (kZ2/rad^2, the rad^2 Nakayama 1->2->3) are GUARANTEED #wide == #torsion and were REMOVED --
# the assert `W.size < len(L.elements)` on them would ALWAYS FAIL. A genuine witness needs a
# representation-INFINITE but tau-tilting-finite algebra: the TAME preprojective algebras
# Pi(D4) (#torsion = |W(D4)| = 192) or Pi(A5) (720). At authoring the shipped P45 engine could
# NOT compute either cheaply -- Pi(D4)/QQ: exchange_graph did not reach 50 vertices in 120s;
# Pi(D4)/GF(31): mutate() raised QuiverlabError deep in the BFS (status="error" after 6
# vertices). So this oracle is DEFERRED to implementation:
#     A = PreprojectiveAlgebra("D4", field=QQ)          # or field=GF(p), p > dim = 28, for speed
#     L, W = torsion_lattice(A, budget=256), wide_subcategories(A, budget=256)   # budget >= 192
#     if L.is_complete and W.is_complete:
#         assert W.size < len(L.elements)               # 192 = #torsion; #wide < 192 is a THEOREM
# marked oracle_literature (the strict inequality is Marks-Stovicek, not an empirical guess). If
# the engine still refuses Pi(D4)/Pi(A5), this stays a documented honest-scope gap and the
# M3/NC structural arbiters above carry the wide-vs-torsion discrimination.
```

- [ ] **Step 2: Run to verify failure**, implement, iterate.
- [ ] **Step 3: Implement.** **TRANSCRIBE Enomoto (2201.00595) §(κ map / core label order)
  VERBATIM before coding** — the naive "intersection of cover-labels over maximal chains" is
  WRONG (it reproduces the *count* on `kAₙ` but not the *poset*; see Methodology). Implement:
  - `_kappa_map(L)`: `κ(j)` for each join-irreducible `j` = the largest element `m` with
    `m ≥ j_*` and `j ⊄ m` (i.e. `m ⊉ j`); extend to canonical join representations per the
    extended-κ definition (Barnard–Todorov–Zhu, as Enomoto states it).
  - `_kappa_order(L)` and `_core_label_order(L)` per Enomoto; **assert they coincide**
    (`constructions_agree`; Enomoto's theorem — a real self-cert, loud on mismatch).
  - `WideSubcategoryPoset` from either (they agree): `elements` = the distinct core-label
    sets, `order` = inclusion, `labels` = the brick-label set per element, `size` = count.
  - `Algebra.wide_subcategories(budget=512)`: lazy delegate.
- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(tautilting): wide_subcategories via the core label order / kappa order (Enomoto) -- #wide + poset; kA2=M3, kA3=NC(A3) structural pins, kappa==CLO self-cert`

**Adjust to reality (Task B):**
- **The wide-COUNT pins on `kAₙ` are literature-anchored but NON-discriminating** (§Verified
  pins). `#wide = #torsion` **iff `A` is representation-FINITE** (Marks–Šťovíček, ruling 1),
  and every `kAₙ` is rep-finite, so on `kAₙ` the count `= #torsion = Catalan(n+1)` no matter
  what — it cannot witness `#wide < #torsion`. The **primary** wide-vs-torsion discriminator is
  therefore the poset STRUCTURE: `test_kA2_wide_poset_is_M3` and `test_kA3_wide_poset_is_NC_A3`
  (ruling 2). A genuine strict-inequality oracle needs a rep-INFINITE τ-tilting-finite algebra
  (`Π(D₄)`/`Π(A₅)`, sourced tame) and is **DEFERRED** (the end-of-file comment block, path b):
  the shipped engine could not compute either cheaply at spec time (QQ too slow; GF(31)
  `mutate()` errored deep in the BFS). Do **NOT** re-introduce `kZ₂/rad²` or a rad²-Nakayama —
  both are rep-finite, so `W.size < len(L.elements)` on them is GUARANTEED to FAIL.
- **The kA₂-M₃ and kA₃-NC(A₃) structural pins** (`test_kA2_wide_poset_is_M3`,
  `test_kA3_wide_poset_is_NC_A3`) are the arbiters that the poset — not just the count — is
  right; my authoring prototype's naive "chain-label-intersection" route reproduced the *counts*
  but produced a WRONG poset (a height-3 poset for kA₂, not M₃), which is exactly why the
  implementation must follow Enomoto's definition. Do NOT ship the naive chain-intersection.
  The `NC(A₃)` structural battery is Whitney numbers `1,6,6,1` + `6` atoms + `6` coatoms +
  self-duality; the definitive arbiter is an order-iso against the hardcoded NC(4) covering
  table DERIVED in the appendix (14 partitions, 28 covers) — add it if the battery ever proves
  insufficient.
- **`kappa == core label order` (Enomoto)** is a strong internal cross-check: if the two
  constructions disagree, ONE of them is mis-transcribed — fix the code, never the assert.
- If Enomoto's construction proves heavier than the S–M budget for `n = 3` (14 elements is
  fine; larger τ-tilting-finite inputs may be slow), cap the poset build at the `budget`
  element count and report `size` + a truncation note honestly; the `kAₙ` pins are small.

---

### Task C: QPA honest-scope skip

**Files:**
- Create: `tests/qpa/test_congruence_qpa.py`

**Interfaces:** QPA 1.37 has **no** torsion-lattice / congruence-lattice / brick-forcing /
core-label-order / wide-subcategory surface (it has **no** τ-tilting surface at all — the
Plan-45 finding, re-confirmed by the P63 probe). Probe live via `NamesGVars()` (the
`tests/qpa/test_products_qpa.py` / Plan-45 / P63 precedent): the test SKIPS honestly and
**FAILS if QPA ever ships** any `TorsionClasses` / `CongruenceLattice` / `WideSubcategories` /
`Bricks` / `ForcingOrder` verb — so the honest-scope claim on the verification page cannot
silently rot.

- [ ] **Step 1: Write the probe** (`pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
  reason=...)`; sweep `NamesGVars()` for the candidate names; `pytest.skip(...)` if absent,
  `assert False` (loud) if any appears). Record that `#torsion` / `#bricks` themselves have no
  QPA cross-check (Plan-45 stated it), so the lattice/congruence invariants inherit no QPA
  oracle — the external cross-checks named are the DIRRT/BCZ worked examples + Enomoto, not a
  live QPA call.
- [ ] **Step 2: Run** — `... -m pytest tests/qpa/test_congruence_qpa.py -q -m qpa` (venv has
  `[qpa]`): expect a clean skip.
- [ ] **Step 3: Commit** — `test(qpa): honest no-torsion-lattice/congruence surface probe (fails if QPA ships one)`

---

### Task D: no-code exposure — the `congruences` compute kind (the "lattice block extension")

The organizing invariant is **byte-identity across the two runners**, achieved by routing
BOTH through the single shared builder `congruence.congruences_block` (the `tau_tilting`
precedent — each runner only adds `citations` locally). Anchors below are current-tree; still
"adjust to reality" if a merge shifts them.

**`congruences` is an ALGEBRA-level compute kind** carrying a PAIR BUDGET (not a degree
range), parsed exactly like `tau_tilting`/`wall_chamber`: `congruences` or `congruences:512`.
NOT a `MODULE_KINDS` entry; sized on `A.dim` like `tau_tilting`. It is the metaplan's
**"lattice block extension"** — it sits in the same lattice/τ-tilting theme as `tau_tilting`
and (if merged) `wall_chamber`, and its block reuses the P45 brick/semibrick labels.

**The `congruences_block(A, budget) -> dict` payload:**
```python
{
  "kind": "congruences", "n": n,
  "complete": bool, "status": "complete"|"budget",
  "lattice": {                                  # None when incomplete
     "size": int, "num_covers": int,
     "is_semidistributive": bool, "is_modular": bool, "is_distributive": bool,
     "join_irreducibles": int, "meet_irreducibles": int,
     "hasse": [ {"from": up, "to": low, "brick_dimvec": {v:m}, "brick_name": str|None}, ...],
     "canonical_joins": [ {"element": id, "label": str,
                           "joinands": [{"brick_dimvec": {v:m}, "brick_name": str|None},...]},...]},
  "congruences": {                              # None when incomplete
     "size": int,                              # |Con(tors A)|
     "num_join_irreducibles": int,             # = #bricks
     "is_distributive": True,
     "forcing_order": {"bricks": [{"dimvec":{v:m},"name":str|None},...],
                       "relations": [[i, j], ...]}},   # brick_i <= brick_j
  "wide": {                                     # None when incomplete
     "size": int,                              # #wide
     "constructions_agree": True,
     "relations": [[a, b], ...], "labels": [[{"dimvec":..,"name":..},...],...]},
  "references": ["dirrt_lattice_torsion","barnard_carroll_zhu","enomoto_wide_ice",
                 "marks_stovicek","air_tau_tilting","demonet_iyama_jasso"],
  "note": str | None,                          # honest truncation when incomplete
}
```
When incomplete: `lattice=congruences=wide=None`, `note` = the honest τ-tilting-infinite
message; the block still returns cleanly (no partial-lattice lie).

**Files (the seven-touchpoint checklist — mirror P45 `tau_tilting` / P63 `wall_chamber`):**
1. **Parse grammar** — `src/quiverlab/hpc/spec.py::parse_compute_item` (the `tau_tilting`
   branch, ~`spec.py:267`): add `if s == "congruences" or s.startswith("congruences:"):`
   returning `ComputeItem(kind="congruences", lo=None, hi=budget)`. Mirror in
   `docs/gui/runner.py` parse (~`runner.py:165`).
2. **Server/HPC dispatch** — `spec.py::_dispatch` (after the `tau_tilting` branch): `if kind
   == "congruences": budget = item.hi or 512; from quiverlab.tautilting.congruence import
   congruences_block; block = congruences_block(A, budget=budget); block["citations"] =
   _citation_pairs(block["references"]); return block, None`. Catch the char-caveat
   `QuiverlabError` into `{"error": "<loud message>"}` (the Plan-30 honest-per-entry
   precedent — never a 500). Add the `_snippet` calls entry: `"congruences": lambda it:
   ("A.congruence_lattice(budget=" f"{it.hi if it.hi is not None else 512})")`.
3. **Pyodide twin** — `docs/gui/runner.py::compute_one` (after the `tau_tilting` branch): the
   byte-identical `elif name == "congruences":` branch calling `congruences_block`; one entry
   in the `calls` reproduce map (`"congruences": "A.congruence_lattice(budget=%d)"`); one
   `ETA_MODEL["scalars"]` cost — `"congruences": 3.0` (the exchange-graph BFS + the
   principal-congruence fixed points + the κ/CLO build; size it just above `wall_chamber`'s
   2.5 / `tau_tilting`'s 2.0).
4. **GUI JS** — `docs/gui/gui.js` AND `webapp/static/gui/gui.js` (vendored byte-identical —
   apply every edit to BOTH): a checkbox `<input id="qlgui-congruences">` + budget picker
   `qlgui-congruences-budget` in the structure/lattice block (beside `tau_tilting`); the two
   bare ids in the id-registry array; the request push-list entry `if (el.congruences.checked)
   compute.push("congruences:" + el["congruences-budget"].value)`; a `renderBlock` `else if
   (name === "congruences")` branch calling a new `renderCongruences(div, b)`; the `THEMES`
   structure `kinds` array entry `"congruences"` inside the `QLGUI-THEMES-BEGIN/END`
   sentinels; the `scheduleProbe`/`KIND_CTRL` entries. **`renderCongruences`**: (a) the
   lattice summary (size, covers, semidistributive/modular/distributive flags, #join-irred =
   #bricks); (b) the **Hasse quiver** of `tors A` rendered (v1) as a **text/HTML edge-list**
   (`from → to — brick label`), mirroring the SHIPPED `renderTauTilting` `<ul>` Hasse-edge
   precedent (`gui.js:~2871`) — **NO bespoke SVG layout in v1** (ruling 6; see the SVG stretch
   task below); the covers may also be grouped by `Gen`-size layer in the table for readability
   (a `<td>` layer column, no coordinate math); (c) the **canonical join representations** table
   (element | its semibrick joinands); (d) the **`Con(tors A)`** summary (`|Con|`, distributive,
   #J = #bricks) + the **forcing order on bricks** as a small edge-list table
   (`brick_i ≤ brick_j`); (e) the **wide-subcategory poset** (#wide + the poset as an edge-list
   / rank table). When incomplete, render the honest `note` and nothing else.
   - **STRETCH (deferred, ruling 6): a layered-DAG SVG Hasse diagram.** This is materially more
     work than the shipped radial wall-and-chamber fan (which just draws rays from the origin):
     it needs (1) layer assignment by `Gen`-size, (2) intra-layer x-positioning to reduce edge
     crossings, (3) edge routing between layers, (4) brick edge-label placement, (5) overlap/
     width management. Budget it as its own GUI subtask if pursued; it is the ONLY client-side
     coordinate math (float-exempt, like the fan). v1 ships the edge-list table above and does
     NOT block on the SVG.
5. **i18n (FOUR locales) + picker template** — `webapp/server/i18n/{en,es,fr,zh}.json` (the
   shipped `LANGS = ("en","es","fr","zh")`): add `pick.kind.congruences` (REQUIRED — gated by
   `test_layout_picker.py`), `struct.congruences` (checkbox label), and `cong.*` block
   strings (`cong.title`, `cong.lattice`, `cong.size`, `cong.covers`, `cong.semidistributive`,
   `cong.join_irred`, `cong.canonical_joins`, `cong.congruence_lattice`, `cong.forcing`,
   `cong.wide`, `cong.wide_poset`, `cong.truncated`, `cong.budget`). Add one
   `data-pick-kind-congruences="{{ t('pick.kind.congruences') }}"` line in
   `webapp/templates/draw.html` (mirror `data-pick-kind-tau_tilting`).
6. **Report renderer** — `src/quiverlab/trace/results_html.py`: `_HEADINGS["congruences"]` +
   a `_congruences_html(b)` branch in `_block_html` rendering the lattice summary, the
   canonical-join table, the `Con` summary + forcing-order poset (as a small table:
   brick_i ≤ brick_j), and the wide-subcategory poset; when incomplete, the honest note. For
   the report PDF Hasse diagram, **REUSE the existing `viz/tikz.py::tikz_hasse` + `poset_layout`**
   (a ranked, float-free Hasse tikzpicture already shipped and "reusable by P45", `tikz.py:53`)
   fed the torsion lattice as a poset — do NOT hand-roll a bespoke `tikz_torsion_lattice` layout.
   Brick edge labels are a small optional extension of `tikz_hasse` (an edge `node`); v1 may omit
   them from the TikZ and carry the labels in the HTML edge-list table. An honest "table only
   (budget-capped)" node when incomplete.
7. **Tests / gates** — add `congruences` to `ALL_KINDS` in
   `tests/webapp/test_layout_picker.py` (must equal the THEMES kind set, each once —
   `test_themes_cover_every_kind_exactly_once` + the four-locale `pick.kind.*` gate); add ONE
   golden `congruences_kA2` to `tests/webapp/_runner_goldens.json` with a dated change-log
   bullet in `tests/webapp/test_runner_delegation.py` (verify existing entries byte-identical
   FIRST; on rebase merge the JSON semantically, never textually); add a twin-parity test
   (copy `test_tau_tilting_p45.py::test_twin_parity`).
- Test files: `tests/webapp/test_congruences_p64.py`,
  `tests/gui/test_congruences_runner_twin.py`.

```python
# tests/webapp/test_congruences_p64.py
"""The congruences algebra-level compute kind: served by hpc.spec, mirrored by the Pyodide
twin, both runners byte-identical. kA2 -> |L|=5 (pentagon), |Con|=5, #wide=5, complete."""


def test_congruences_block_shape(tmp_path):
    from quiverlab.hpc.spec import ComputeRequest, run
    req = _congruences_request(quiver=("kA2",), budget=512)
    out = run(ComputeRequest.model_validate(req), tmp_path)
    b = out["results"]["congruences"]
    assert b["complete"]
    assert b["lattice"]["size"] == 5 and b["lattice"]["join_irreducibles"] == 3
    assert b["congruences"]["size"] == 5 and b["congruences"]["num_join_irreducibles"] == 3
    assert b["wide"]["size"] == 5
    assert "dirrt_lattice_torsion" in b["references"]


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py::compute_one and hpc.spec.run;
    # assert json.dumps(block, sort_keys=True) equality on the congruences block.
    ...


def test_tau_tilting_infinite_status(tmp_path):
    # 2-Kronecker with a small budget -> complete False, status "budget", note set,
    # lattice/congruences/wide all null, no crash.
    ...
```

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir).
- [ ] **Step 2: Implement** the parse grammar (both runners), the server + twin dispatch, the
  GUI JS in both `gui.js` files, the four-locale i18n + `draw.html`, and the
  `results_html.py` + `viz/tikz.py` branches. Set `block["kind"] = "congruences"` inside the
  builder so both runners agree.
- [ ] **Step 3: Gate the layout/i18n first** —
  `... -m pytest tests/webapp/test_layout_picker.py tests/webapp/test_js_parses.py -q`
  (THEMES == `ALL_KINDS`; `pick.kind.congruences` in all four locales).
- [ ] **Step 4: Add the golden** (`congruences_kA2`) to `_runner_goldens.json`; note it in the
  `test_runner_delegation.py` docstring change-log. Run the delegation test BEFORE adding to
  confirm existing entries stay byte-identical.
- [ ] **Step 5: Run the gates** —
  `... -m pytest tests/webapp/test_congruences_p64.py tests/webapp/test_runner_delegation.py tests/gui/test_congruences_runner_twin.py tests/hpc -q`.
  Expected PASS; both runners byte-identical; the τ-tilting-infinite case returns
  `status="budget"` cleanly.
- [ ] **Step 6: Commit** — `feat(gui,webapp,hpc,trace): congruences compute kind -- torsion lattice + Con + forcing order + wide subcategories, Hasse edge-list table + tikz_hasse reuse, one golden (4 locales, v1 schema)`

**Adjust to reality (Task D):**
- `test_layout_picker.py::ALL_KINDS` is a frozen set the THEMES array must equal exactly; add
  `congruences` to BOTH `gui.js` files (vendored byte-identical) or the gate fails.
- If `estimator.sizing_dim` (`webapp/server/estimator.py`) needs an entry, size `congruences`
  on `A.dim` exactly as `tau_tilting` is sized (algebra-level; not a module kind) — check the
  `tau_tilting` path and mirror it.
- Confirm the exact `parse_compute_item` / `_dispatch` / `_snippet` / `compute_one` / `calls`
  / `ETA_MODEL` line anchors on first read (they drift with merges); the `tau_tilting` and
  `wall_chamber` (if P63 merged) branches are the templates. The block already carries
  `"kind": "congruences"` if the builder sets it — mirror `tau_tilting_block`.

---

### Task E: verification page, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick the P64 card), and append the
  dated R26 correction note to `docs/plans/2026-08-06-computability-expansion-deep-research.md`
  (the "core label order = Enomoto, not BCZ" finding; metaplan §2 mandates folding corrections
  back).
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`).

- [ ] **Step 1: Verification page.** Add the Plan-64 subsystem row (`tautilting/congruence.py`):
  - `oracle_literature` — the **kA₂ pentagon** (`tors kA₂` = N₅: 5 elements, semidistributive,
    NOT modular, NOT distributive); `|Con(tors kA₂)| = 5`; the **kA₃ congruence lattice**
    (`|tors kA₃| = 14`, `|Con(tors kA₃)| = 14`); the wide-subcategory counts (kA₂ → 5 ≅ M₃,
    kA₃ → 14 ≅ NC(A₃)); the τ-tilting-infinite honest refusal (2-Kronecker).
  - `oracle_crossengine` — **join-irreducibles ↔ bricks** (lattice count vs P45 `bricks`:
    kA₂ → 3, kA₃ → 6); **canonical joins ↔ semibricks** (per-element label-set equality vs P45
    `semibricks`, not just cardinality); **#J(Con) = #bricks** (the forcing order lives on
    bricks). (The strict-inequality `#wide < #torsion` discriminator is NOT a cross-engine
    oracle here — it is DEFERRED, path b; see the honest-scope entry (d).)
  - `oracle_selfcert` — `is_lattice` (unique join/meet); `Con` distributive
    (Funayama–Nakayama, structural); the forcing-poset order-ideal count = `|Con|`; the κ
    order = the core label order (Enomoto); the `#J(Con) == #bricks` / `#join-irred ==
    #bricks` internal gates.
  - the **honest semi-decision entry** (metaplan §6): certified complete iff τ-tilting-finite
    (DIJ); on the 2-Kronecker (τ-tilting-infinite) NO lattice invariant is emitted (a partial
    congruence lattice is meaningless) — `status="budget"`, `note`, all invariants `None`
    (STRICTER than P63's bounded region, and stated as such).
  - the **honest-scope entries**: (a) rigorous only over char 0 / char > dim (QQ default;
    loud off scope, inherited from Plan 45); (b) bricks decide over the algebraically-closed /
    char-0 base (`end_dim = 1`; the GF(pⁿ) proper-division-ring caveat); (c) **QPA cannot
    compare** — the honest no-surface probe (DIRRT/BCZ/Enomoto worked examples are the named
    external cross-checks, not a live oracle); (d) the **wide/CLO poset structure** pins
    (kA₂ ≅ M₃, kA₃ ≅ NC(A₃)) are literature-anchored (Enomoto + Ingalls–Thomas'
    noncrossing-partition / wide-subcategory correspondence for type A), and are the PRIMARY
    wide-vs-torsion discriminators; (e) **`#wide = #torsion` iff `A` is representation-FINITE**
    (Marks–Šťovíček, arXiv:1503.04639 — NOT merely for hereditary `A`), so the count never
    discriminates on the rep-finite `kAₙ`, and the strict-inequality `#wide < #torsion` oracle is
    **DEFERRED**: it requires a rep-INFINITE τ-tilting-finite algebra (tame `Π(D₄)`/`Π(A₅)`),
    which the shipped engine could not compute at spec time (QQ too slow; GF(31) `mutate` error)
    — a documented honest-scope gap until the implementation-time run lands.
  Recount the class table (`tests/release/test_oracle_classes.py` drives the numbers — collect,
  paste the LIVE counts, re-run to green; mid-merge-train honest).
- [ ] **Step 2: README.** One features line (rep-theory-first, per the positioning memory):
  "the lattice theory of torsion classes (Demonet–Iyama–Reading–Reiten–Thomas): the congruence
  lattice `Con(tors A)`, the forcing order on bricks, canonical join representations, and the
  wide-subcategory poset (core label order, Enomoto) — one click via `congruences`, certified
  complete iff τ-tilting-finite."
- [ ] **Step 3: Full gate:**
  `... -m pytest tests/modules/test_torsion_lattice.py tests/modules/test_congruence_lattice.py tests/modules/test_wide_subcategories.py -q` (deep, the touched files),
  `... -m pytest -q -m fast`, `... -m pytest tests/webapp tests/gui -q`,
  `... -m pytest tests/qpa -q -m qpa`, `... -m pytest tests/release tests/citations -q` —
  all green. Plan 45 left byte-unchanged (spot-run `tests/modules/test_tau_tilting_*.py`).
- [ ] **Step 4: Commit** — `docs(verification): Plan-64 torsion-lattice/congruence oracle rows + DIRRT/BCZ/Enomoto citations + honest scope (char caveat, bricks base, tau-tilting-finite gate, no-QPA, wide-poset literature) + recounted classes`

---

## Acceptance (Plan-64 definition of done)

1. `torsion_lattice(A)` returns `TorsionLattice` — the finite lattice `tors A` extracted from
   the P45 oriented exchange graph: order, join/meet (`is_lattice` self-cert),
   join-irreducibles, meet-irreducibles, canonical join representations (down-covers +
   semibrick labels), and the SD/modular/distributive flags — all exact, float-free.
2. **kA₂ = the pentagon N₅** pinned (5 elements, semidistributive, NOT modular, NOT
   distributive, 3 join-irreducibles = 3 bricks, 5 covers) — derived in the Verified-pins
   section and re-derived by the engine.
3. `congruence_lattice(A)` returns `CongruenceLattice` — `Con(tors A)` via principal
   cover-congruences, distributive (Funayama–Nakayama, structural + verified), its
   join-irreducibles = the forcing classes = **bricks**, the **forcing order on bricks**, and
   `|Con(tors A)|` = the forcing poset's order-ideal count. **`|Con(tors kA₂)| = 5`,
   `|Con(tors kA₃)| = 14`** pinned; **#J(Con) = #bricks** (kA₂ → 3, kA₃ → 6) cross-engine.
4. **join-irreducibles ↔ bricks** (BCZ) pinned cross-engine (lattice vs P45 `bricks`);
   **canonical joins ↔ semibricks** pinned cross-engine (lattice down-covers vs P45
   `semibricks`; `#semibricks = |L|`).
5. `wide_subcategories(A)` returns `WideSubcategoryPoset` — the wide-subcategory poset via
   Enomoto's core label order AND κ order, **with the two constructions asserted to agree**
   (self-cert). `#wide(kA₂) = 5` (poset ≅ **M₃**, structurally pinned: 3 pairwise-incomparable
   atoms, each covering ⊥ and covered by ⊤) and `#wide(kA₃) = 14` (poset ≅ **NC(A₃)**,
   structurally pinned: Whitney `1,6,6,1`, 6 atoms, 6 coatoms, self-dual) — the poset STRUCTURE
   is the wide-vs-torsion discriminator (ruling 2). `#wide = #torsion` **iff `A` is
   representation-FINITE** (Marks–Šťovíček — NOT merely hereditary), so on the rep-finite `kAₙ`
   the count never discriminates; a strict-inequality `#wide < #torsion` witness needs a
   rep-INFINITE τ-tilting-finite algebra (tame `Π(D₄)`/`Π(A₅)`) and is **DEFERRED** to
   implementation (the shipped engine could not compute it at spec time — honest-scope gap).
6. The structure is **certified complete iff τ-tilting-finite** (DIJ, decided by the P45
   BFS); on the 2-Kronecker it emits **no** lattice invariant (`is_complete=False`,
   `status="budget"`, `note`, all `None`) — the honest refusal, never a partial-lattice lie.
7. `congruences` clickable end-to-end (GUI canvas → block → Hasse **edge-list table** +
   `Con`/forcing/wide rendering → report `tikz_hasse`) in all four locales (en/es/fr/zh), both
   runners byte-identical via the shared builder, ONE golden added with a documented change-log
   entry, `test_layout_picker.py` `ALL_KINDS`/THEMES gate green, schema still v1, canonical keys
   unchanged. (The layered-DAG SVG Hasse diagram is a deferred stretch task — ruling 6 — not a
   v1 acceptance requirement.)
8. QPA honest-scope probe skips (fails if QPA ever ships a torsion-lattice/congruence
   surface); `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the
   honest scope; README line added; the R26 "core label order = Enomoto, not BCZ" correction
   folded back into the research doc; deep (P64 files) + fast + webapp/gui + qpa + release +
   citations suites green. Plan 45 left byte-unchanged. Honest scope recorded: char 0 / char
   > dim (QQ default); bricks over the algebraically-closed base; τ-tilting-finite gate; QPA
   cannot compare; the wide/CLO poset-structure pins literature-anchored.

---

## Methodology & assumptions

**Approach.** I read the R26 record and the metaplan P64 card; then Plan 45 (the C4 τ-tilting
engine) end-to-end via its shipped code — `src/quiverlab/tautilting/{__init__,torsion,
mutation,pairs,block}.py` — and confirmed the exact APIs P64 consumes: `exchange_graph` /
`ExchangeGraph` (the `is_complete`/`status`/`n_regular` contract; `.vertices[i]["pair"]`,
`.arrows[(i,j)]["brick"]`, `.adj`), `hasse_orientation` (up/down by `Gen`-size),
`bricks` / `semibricks` (iso-class deduped), `torsion_class_data`, `_edge_brick` /
`_torsion_universe` (the non-thin-safe edge-brick label), `is_isomorphic` /
`identify_standard`. I read Plan 61 and Plan 63 (the closest siblings) end-to-end for the
house style, the honest-scope discipline, the citation/`_r` conventions, and the seven-touchpoint
compute-kind wiring, and I mapped that wiring by grepping `hpc/spec.py`, `docs/gui/runner.py`,
`docs/gui/gui.js`, the four-locale i18n, `trace/results_html.py`, `viz/tikz.py`, and
`tests/webapp/test_layout_picker.py`. I checked `citations/registry.py` — DIRRT / BCZ /
Enomoto are ABSENT (P45 ships only `air_tau_tilting` / `demonet_iyama_jasso` /
`king_stability`; P63 adds BST / Asai / KT, disjoint from P64), so P64 adds its three.

**Empirical ground truth (run this authoring, in the venv).** I built kA₂ and kA₃ via
`linear_path_algebra` and extracted the torsion lattice from the P45 exchange graph
(oriented Hasse → transitive closure → join/meet), then prototyped the principal-congruence /
forcing / order-ideal-count and a wide-count pass. **Observed live:** kA₂ → `|tors| = 5`,
`is_lattice`, **semidistributive True, modular False, distributive False** (the pentagon N₅),
3 join-irreducibles == 3 bricks, 5 semibricks, **`|Con| = 5`**, 3 join-irreducible congruences
(= 3 bricks) whose forcing poset is a "V" opening upward — **one brick below the other two**,
relations `[(2,0),(2,1)]` (the below-brick's `con` is contained in / forced by the other two;
ruling 4) — 2 relations, 5 order ideals; kA₃ → `|tors| = 14`, 6 join-irreducibles == 6 bricks,
14 semibricks, **`|Con| = 14`**, 6 join-irreducible congruences with a **9-relation** forcing
poset. Every table entry in "Verified pins" except the two `# to-verify` wide-poset *structures*
was observed live. **The prototype that reproduces these is shipped verbatim in the appendix
(ruling 5) and was re-run this fix-round — all numbers above confirmed byte-for-byte.**

**Key correctness decisions.** (1) **P64 is a lattice-theoretic layer, not new representation
theory** — the brick labels ship in P45; P64 assembles the abstract lattice and computes
`Con`, forcing, canonical joins, and the wide poset. (2) **The oriented exchange graph IS the
Hasse quiver of `tors A`** (AIR/DIRRT), so the lattice order = its transitive closure;
`is_lattice` is a real self-cert (DIRRT guarantee it holds). (3) **`Con` via principal
cover-congruences + Birkhoff down-sets** is the standard finite-lattice algorithm; DIRRT
supply the representation-theoretic identification `J(Con) ↔ (bricks, forcing order)`, giving
the cross-engine `#J(Con) = #bricks` pin. (4) **The completeness gate is STRICTER than P63's**
— a congruence lattice / forcing order / wide poset is a GLOBAL function of the whole finite
lattice, so on a truncated exchange graph P64 emits NO invariant (a partial `Con` is a lie,
not a bounded region). (5) **Attribution corrected (standing rule):** the metaplan card cites
"Barnard–Carroll–Zhu" for core-label-order = wide, but that result is **Enomoto (2201.00595)**
(BCZ 2019 gives join-irreducibles↔bricks + the canonical join complex; the κ map is
Barnard–Todorov–Zhu) — P64 cites all three correctly and folds the correction back into R26.

**Deliberately not checked / assumptions.** (a) **The wide/CLO POSET STRUCTURE.** My authoring
prototype's naive "intersection of cover-labels over maximal chains" reproduced the wide
*counts* (kA₂ → 5, kA₃ → 14) but produced the WRONG *poset* (a height-3 poset for kA₂, not
M₃) — so I explicitly REJECT the naive route and mandate transcribing Enomoto's κ-map / core-
label definitions verbatim at implementation, with the two-constructions-agree self-cert and
the **kA₂-M₃ and kA₃-NC(A₃) structural tests** as arbiters. This is the **primary design risk**
(flagged in Task B). (b) **The `#wide < #torsion` discriminator is DEFERRED (ruling 1, path b).**
The prior draft named `kZ₂/rad²` and a rad²-Nakayama as the discriminating algebras — both are
**representation-FINITE**, and `#wide = #torsion` *iff* `A` is rep-finite (Marks–Šťovíček,
arXiv:1503.04639 — the boundary is rep-finiteness, NOT heredity), so those candidates are
GUARANTEED `#wide = #torsion` and the strict inequality would ALWAYS FAIL. A genuine witness
needs a rep-INFINITE τ-tilting-finite algebra: the **tame** (hence rep-infinite; Dlab–Ringel /
Schröer FD-Atlas: `Π(Δ)` rep-finite iff `Δ = Aₙ, n ≤ 4`, tame iff `Δ ∈ {A₅, D₄}`, else wild)
preprojective algebras `Π(D₄)` (`#torsion = |W(D₄)| = 192`) and `Π(A₅)` (`720`), where
Marks–Šťovíček *forces* `#wide < #torsion`. I PROBED both at spec time: the shipped P45 engine
cannot compute either cheaply — `Π(D₄)/QQ` did not reach 50 exchange-graph vertices in 120 s;
`Π(D₄)/GF(31)` (char 31 > dim 28, so in P45 scope) raised `QuiverlabError` in `mutate()` deep in
the BFS (`status="error"` after 6 vertices) — and I found no published explicit `#wide` count for
a preprojective algebra. So I could not verify a witness or pin a literature value at spec time,
and per the ruling I take **path (b)**: DROP the strict-inequality *required* oracle, shift the
wide-vs-torsion discrimination onto the strengthened M₃/NC(A₃) structural arbiters (ruling 2),
and record `Π(D₄)`/`Π(A₅)` as the implementation-time candidates (to-be-verified). For every `kAₙ`
(rep-finite) `#wide = #torsion = Catalan(n+1)` — the COUNT is non-discriminating; the POSET is.
(c) I did not run the deep exchange
graph beyond kA₃ (the BFS is budget-capped; larger τ-tilting-finite inputs are the
implementer's spot-check). (d) The `Enomoto2023wide` AlgRepThy volume/year/pages and the
optional Barnard–Todorov–Zhu κ-map entry are `# verify`/`# PIN` at citation-add (arXiv id +
DOI are firm; no fabricated print metadata — house rule). (e) The `O(N³)` lattice-property
sweep and the order-ideal enumeration are capped at the `budget` element count (fine for
τ-tilting-finite small `N`; honest `None` + note beyond).

**Why this is correct.** The lattice / congruence / forcing / join-irreducible framework is
verbatim-transcribed from DIRRT (Trans. AMS B 10 (2023), abstract quoted) and BCZ (Alg.
Combin. 2 (2019), abstract quoted), the wide/CLO result from Enomoto (2201.00595, abstract
quoted); the finite-lattice `Con` algorithm is textbook (Funayama–Nakayama + Birkhoff +
Grätzer). Every count that P64 asserts on kA₂/kA₃ — `|tors|`, SD/non-modular/non-distributive,
`#join-irreducibles`, `#semibricks`, `|Con|`, `#J(Con)`, the forcing-poset shape — was
recomputed live against the shipped Plan-45 engine this authoring and agrees with the
literature (the pentagon; Catalan; DIRRT). The only object whose exact algorithm rests on
implementation-time transcription (Enomoto's core label order) is fenced by two independent
constructions that must agree and the structural **M₃ / NC(A₃)** arbiters — the honest boundary
of what I could verify at spec time.

## Appendix: authoring/fix-round prototype + NC(4) table (reproducible pins — ruling 5)

**A1. Runnable prototype** — reproduces `|tors|`, `is_lattice`, `#join-irreducibles = #bricks`,
`#semibricks = |tors|`, `|Con|`, `#J(Con) = #bricks`, and the forcing poset (relation count AND
direction) for kA₂/kA₃, straight from the shipped Plan-45 exchange graph. Run:
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python thisfile.py`.

```python
from itertools import product
from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.mutation import exchange_graph
from quiverlab.tautilting.torsion import hasse_orientation, bricks, semibricks


def extract(A, budget=512):
    eg = exchange_graph(A, budget_pairs=budget)
    assert eg.is_complete, "tau-tilting-infinite: no finite lattice"
    orient = hasse_orientation(eg)
    n = len(eg.vertices)
    covers = []                                   # (upper, lower, brick_dimvec, brick_name)
    for (i, j), d in orient.items():
        a, b = (i, j) if d == "down" else (j, i)  # a (upper) > b (lower)
        lab = eg.arrows[(i, j)]
        covers.append((a, b, tuple(sorted((lab["brick"] or {}).items())), lab["brick_name"]))
    down = {e: {e} for e in range(n)}             # down[e] = {f : f <= e}
    changed = True
    while changed:
        changed = False
        for (a, b, _, _) in covers:
            if not down[b] <= down[a]:
                down[a] |= down[b]; changed = True
    up = {e: {f for f in range(n) if e in down[f]} for e in range(n)}
    return n, covers, down, up


def meet(x, y, down):
    common = down[x] & down[y]
    cand = [z for z in common if all(w in down[z] for w in common)]
    assert len(cand) == 1, ("meet not unique", x, y, cand); return cand[0]


def join(x, y, up, down):
    common = up[x] & up[y]
    cand = [z for z in common if all(z in down[w] for w in common)]
    assert len(cand) == 1, ("join not unique", x, y, cand); return cand[0]


def is_lattice(n, down, up):
    try:
        for x in range(n):
            for y in range(n):
                meet(x, y, down); join(x, y, up, down)
        return True
    except AssertionError:
        return False


class UF:
    def __init__(self, n): self.p = list(range(n))
    def find(self, a):
        while self.p[a] != a: self.p[a] = self.p[self.p[a]]; a = self.p[a]
        return a
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb: self.p[ra] = rb; return True
        return False
    def partition(self, n):
        blocks = {}
        for e in range(n): blocks.setdefault(self.find(e), set()).add(e)
        return frozenset(frozenset(v) for v in blocks.values())


def principal_congruence(a0, b0, n, down, up):
    uf = UF(n); uf.union(a0, b0); changed = True
    while changed:
        changed = False
        reps = {}
        for e in range(n): reps.setdefault(uf.find(e), []).append(e)
        pairs = [(blk[i], blk[j]) for blk in reps.values()
                 for i in range(len(blk)) for j in range(i + 1, len(blk))]
        for (x, y) in pairs:
            for z in range(n):
                if uf.union(meet(x, z, down), meet(y, z, down)): changed = True
                if uf.union(join(x, z, up, down), join(y, z, up, down)): changed = True
    return uf.partition(n)


def refines(ti, tj):                               # ti <= tj iff ti subset tj (ti finer)
    return all(any(bi <= bj for bj in tj) for bi in ti)


def order_ideals_count(idx, leq):
    below = {e: {f for f in idx if leq(f, e) and f != e} for e in idx}
    return sum(1 for m in product([0, 1], repeat=len(idx))
               if all(below[e] <= {idx[k] for k in range(len(idx)) if m[k]}
                      for e in {idx[k] for k in range(len(idx)) if m[k]}))


def analyze(name, A):
    n, covers, down, up = extract(A)
    lc = {}
    for (a, b, dv, nm) in covers: lc.setdefault(a, []).append((b, dv, nm))
    join_irr = [a for a in range(n) if len(lc.get(a, [])) == 1]
    prin = [(principal_congruence(a, b, n, down, up), dv) for (a, b, dv, nm) in covers]
    JCon = []
    for (theta, dv) in prin:
        if not any(t == theta for (t, _) in JCon): JCon.append((theta, dv))
    idx = list(range(len(JCon)))
    leq = lambda i, j: refines(JCon[i][0], JCon[j][0])
    rels = [(i, j) for i in idx for j in idx if i != j and leq(i, j)]
    print(f"=== {name} ===")
    print(f"  |tors|={n}  is_lattice={is_lattice(n, down, up)}  "
          f"#join-irr={len(join_irr)} (#bricks={len(bricks(A))})  "
          f"#semibricks={len(semibricks(A))} (=|tors|? {len(semibricks(A)) == n})")
    print(f"  |Con|={order_ideals_count(idx, leq)}  #J(Con)={len(JCon)}  "
          f"#forcing_relations={len(rels)}  relations(i<=j)={rels}")


if __name__ == "__main__":
    analyze("kA2 (1->2)", Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ))
    analyze("kA3 (1->2->3)", linear_path_algebra(3, field=QQ))
```

**A2. Observed output** (re-run in the venv this fix-round, 2026-08-08):

```
=== kA2 (1->2) ===
  |tors|=5  is_lattice=True  #join-irr=3 (#bricks=3)  #semibricks=5 (=|tors|? True)
  |Con|=5  #J(Con)=3  #forcing_relations=2  relations(i<=j)=[(2, 0), (2, 1)]
=== kA3 (1->2->3) ===
  |tors|=14  is_lattice=True  #join-irr=6 (#bricks=6)  #semibricks=14 (=|tors|? True)
  |Con|=14  #J(Con)=6  #forcing_relations=9  relations(i<=j)=[(3, 0), (3, 1), (4, 1), (4, 2), (5, 0), (5, 1), (5, 2), (5, 3), (5, 4)]
```

The kA₂ relations `[(2,0),(2,1)]` (brick 2 the shared lower endpoint) pin the forcing DIRECTION
"one below two" (ruling 4); the kA₃ `9` relations pin the forcing poset size. `|Con|` = 5 / 14
and `#J(Con) = #bricks` = 3 / 6 reproduce exactly.

**A3. Derived NC(A₃) = NC(4) covering table** (the hardcoded arbiter for
`test_kA3_wide_poset_is_NC_A3`, ruling 2) — the 14 noncrossing partitions of `{1,2,3,4}` ordered
by refinement (bottom = all singletons; rank = 4 − #blocks = #simple objects of the wide subcat):

| id | partition | rank | lower covers |
|----|-----------|------|--------------|
| ⊥  | `1|2|3|4` | 0 | — |
| a1 | `12|3|4`  | 1 | ⊥ |
| a2 | `13|2|4`  | 1 | ⊥ |
| a3 | `14|2|3`  | 1 | ⊥ |
| a4 | `23|1|4`  | 1 | ⊥ |
| a5 | `24|1|3`  | 1 | ⊥ |
| a6 | `34|1|2`  | 1 | ⊥ |
| c1 | `123|4`   | 2 | a1, a2, a4 |
| c2 | `124|3`   | 2 | a1, a3, a5 |
| c3 | `134|2`   | 2 | a2, a3, a6 |
| c4 | `234|1`   | 2 | a4, a5, a6 |
| c5 | `12|34`   | 2 | a1, a6 |
| c6 | `14|23`   | 2 | a3, a4 |
| ⊤  | `1234`    | 3 | c1, c2, c3, c4, c5, c6 |

`13|24` and `24|13` are CROSSING and excluded (that is why a2/a5 have only two rank-2 covers and
c5/c6 only two lower covers). Whitney numbers by rank: `1, 6, 6, 1` (Narayana `N(4,k)`); 6 atoms,
6 coatoms; 28 covers total; self-dual. The wide-subcategory poset of `kA₃` is order-isomorphic to
this (Ingalls–Thomas: `wide(kAₙ) ≅ NC(Aₙ)`); an order-iso against this table is the definitive
kA₃ structural arbiter if the Whitney/atom/coatom/self-dual battery ever proves insufficient.

## Change log

- **2026-08-07 authoring** — plan drafted against merged `dev` (Plan 45 + Wave-1/2). Scope
  boundary vs Plan 45 fixed (P64 = the lattice-theoretic layer: `Con`, forcing order,
  canonical joins, wide poset; P45 ships the brick labels only). The τ-tilting-finite gate
  set STRICTER than P63 (no partial lattice). All numeric pins (kA₂ pentagon N₅, `|Con|` 5/14,
  join-irred=bricks 3/6) recomputed live in the venv against the shipped exchange graph.
- **2026-08-07 reference re-verification folded in** — DIRRT (Trans. AMS B 10 (2023) 542–612,
  btran/100) and BCZ (Alg. Combin. 2 (2019) 879–901, alco.72) confirmed verbatim; the
  metaplan's "core label order = Barnard–Carroll–Zhu" attribution CORRECTED to **Enomoto
  arXiv:2201.00595** (AlgRepThy, s10468-023-10214-0), with the κ map credited to
  Barnard–Todorov–Zhu — a dated correction to research-doc R26 is queued for the merge.
- **2026-08-08 critic fix-round (NEEDS WORK → adjudicated, 6 rulings applied; doc-only).**
  (1 BLOCKING) The `#wide < #torsion` "non-hereditary" discriminator was mathematically
  impossible as shipped: **Marks–Šťovíček** (arXiv:1503.04639, web-verified) makes `#wide =
  #torsion` hold **iff `A` is representation-FINITE** — the boundary is rep-finiteness, not
  heredity — so both named candidates (`kZ₂/rad²`, the rad²-Nakayama) are rep-finite and the
  `assert W.size < len(L.elements)` was GUARANTEED to fail. Fixed EVERY occurrence of the false
  framing (spec §4, Verified-pins note, Task-B tests/adjust-to-reality, verification rows,
  Acceptance #5, Methodology). Took **path (b)** (DROP the required strict-inequality oracle):
  probed the sourced witnesses `Π(D₄)`/`Π(A₅)` (tame ⇒ rep-infinite; Dlab–Ringel/Schröer
  classification web-verified) and found the shipped engine CANNOT compute them cheaply at spec
  time — `Π(D₄)/QQ` did not reach 50 exchange-graph vertices in 120 s; `Π(D₄)/GF(31)`
  (char 31 > dim 28) raised `QuiverlabError` in `mutate()` deep in the BFS (`status="error"`,
  6 vertices) — and no explicit literature `#wide` count surfaced, so the strict-inequality
  oracle is deferred to implementation with `Π(D₄)`/`Π(A₅)` named + rationale, and the
  discrimination burden shifted onto the structural arbiters. (2 MAJOR) `test_kA2_wide_poset_is_M3`
  now asserts the full M₃ shape (3 pairwise-incomparable atoms, each covering ⊥ and covered by ⊤,
  rank histogram `1,3,1`) and a new `test_kA3_wide_poset_is_NC_A3` pins NC(A₃) (Whitney `1,6,6,1`,
  6 atoms, 6 coatoms, self-dual) with a derived hardcoded NC(4) covering table (appendix A3) as
  the order-iso arbiter. (3 MINOR) `test_canonical_joins_are_semibricks` now checks per-element
  label-set equality vs P45 `semibricks`, not just cardinality. (4 MINOR) forcing DIRECTION
  reconciled to the definition-consistent "one brick below the other two" (the forced brick's
  `con` is contained in / forced BY the other two), pinned in `test_kA2_forcing_order_is_a_V`
  (one shared lower endpoint, two distinct uppers). (5 MINOR) the authoring prototype is shipped
  runnable in appendix A1 with re-run output A2 — `|Con|` 5/14, `#J(Con)=#bricks` 3/6, forcing
  relations 2 (`[(2,0),(2,1)]`) / 9 all reproduced this fix-round. (6 MINOR) the Hasse rendering
  is honestly costed: v1 is a **text edge-list** (mirroring the shipped `renderTauTilting` `<ul>`
  precedent) + report reuse of the existing `viz/tikz.py::tikz_hasse`; the layered-DAG SVG is a
  clearly-scoped deferred stretch task (layer assignment + positioning + edge labels enumerated).
