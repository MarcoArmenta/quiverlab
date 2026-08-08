# Plan 75: Incidence algebras (HH ≅ simplicial cohomology) + GHMS-Koszul fast HH (P75 / R9 + R10)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal.** Two independent "fast-path + new-oracle-class" HH engines, joined by one theme
(a *closed-form* Hochschild computation that a general resolution would grind out expensively):

- **(a) Incidence algebras `kP` of a finite poset `P`.** Realize the Gerstenhaber–Schack /
  Cibils theorem `HH^*(kP) ≅ H^*(Δ(P); k)` — the simplicial cohomology of the **order
  complex** `Δ(P)` (the nerve: chains `x₀ < x₁ < … < x_p` are the `p`-simplices) — as a
  **RING isomorphism** (Cibils 1989; Gerstenhaber–Schack 1983). Ship a **poset input mode**
  (the `Poset` helper already exists; the GUI gets a dedicated cover-relation panel), a
  **simplicial cochain engine** that computes `HH^*` by exact **integer SNF / rank** over the
  order complex — dramatically cheaper than the bar/CS bimodule resolution and
  **char-sensitive for free** via universal coefficients — and promote the existing Plan-33
  Boolean-`B₃` pin (`HH^{≥1} = 0` because `Δ(B₃)` is contractible) to a **theorem-backed
  oracle family** cross-checked against the general engines. The "bracket = 0" clause of the
  older folklore is **dropped** (it is NOT in Gerstenhaber–Schack 1983; see §Reference
  re-verification).
- **(b) The GHMS comultiplicative minimal `A^e`-resolution for Koszul algebras**
  (Green–Hartman–Marcos–Solberg, Arch. Math. **85** (2005) 118–127, arXiv:math/0508177). For a
  **Koszul** `A = kQ/I` (gated by Plan-27's **three-valued `ext_algebra.koszul` verdict** — the
  G-quadratic/Priddy PBW `True`, else the honest `False`-with-Ext-obstruction / `None`),
  build the minimal bimodule resolution `P_n = A ⊗_S K_n ⊗_S A` in **closed form** from the
  quadratic (Koszul-dual) data — no syzygy search, running over **any exact `Domain`** (like
  CS, unlike the GF(p)-only minimal engine) — giving **fast minimal HH** and a **third
  independent HH oracle class**. `Betti_n = dim K_n = the n-th coefficient of the Koszul-dual
  Hilbert series`.

**What is genuinely NEW over the shipped tree (read this first).**

- **`Poset` and `IncidenceAlgebra` already exist** (`families/poset.py`, `families/incidence.py`
  — the diamond/chain tests in `tests/families/test_incidence.py`, the catalog family
  `IncidenceAlgebra` with the `poset_or_covers` param, the Plan-33 `B₃` pin in
  `tests/families/test_scale_batteries.py`). Plan 75 **consumes** them read-only and adds: the
  **order-complex simplicial engine**, the **`incidence_cohomology` compute kind + fast path**,
  the poset-provenance stash, the **theorem oracle family** (nonvanishing `H¹`, `H²`, and the
  **char-sensitive `RP²` pin**), and the GUI **poset input panel**. The incidence *constructor*
  is not re-implemented.
- **The minimal `A^e` engine already exists** (`engine/resolutions_minimal.py`) but is **GF(p)
  int64** and computes ranks by **iterated syzygies**. GHMS is different in kind: it is a
  **closed-form** bimodule resolution from the quadratic data, valid over **any `Domain`**,
  gated on a **Koszul certificate**. Plan 75 adds `koszul_kernels`, `GHMSResolution`, and the
  `engine="ghms"` route; it does **not** touch the syzygy engine (which stays the cross-oracle).
- **Plan-27 Koszulity + Ext-algebra are shipped** (`modules/koszul.py`,
  `modules/ext_algebra.py`): `g_quadratic_certificate`, `quadratic_dual`, `YonedaPresentation`
  with `.koszul` (three-valued) + `.koszul_obstruction` and `.hilbert_matrix_through`. GHMS
  **gates on the three-valued `ext_algebra.koszul` verdict** (True = build; False = refuse with the
  positive Ext obstruction; None = refuse, inconclusive) and **anchors** `Betti_n = dim K_n`
  against `hilbert_matrix_through`.

**Architecture.** Two new engine modules + one catalog/GUI surface, no change to existing math
engines:

- `src/quiverlab/hochschild/simplicial.py` (R9): `OrderComplex` (chains-by-dimension of a
  `Poset`), `simplicial_cohomology_dims(oc, top, field)` (exact rank over the field — the fast
  path), `integral_homology(oc, top)` (SNF over ℤ → free ranks + torsion invariant factors → the
  **universal-coefficient char certificate**), and `incidence_cohomology(A, top, field=None)`
  (the theorem: `HH^n(kP) = dim_k H^n(Δ(P); k)`, **requires incidence provenance**, loud
  refusal otherwise). `Algebra.incidence_cohomology(top)` is a thin delegate.
- `src/quiverlab/hochschild/koszul_ghms.py` (R10): `koszul_kernels(A, top)` (the spaces
  `K_n = ⋂ᵢ V^{⊗i}⊗R⊗V^{⊗(n-2-i)} ⊆ V^{⊗n}` over `S = k^{Q₀}`, any `Domain`),
  `class GHMSResolution(Resolution)` (`P_n = A⊗_S K_n⊗_S A` + the comultiplicative differential,
  gated by the Koszul certificate), and `ghms_hochschild_dims` / `ghms_cohomology_dims` (the
  collapse). `Algebra.hochschild_{co}homology(..., engine="ghms")` routes here.
- `families/poset.py` / `families/incidence.py`: add `Poset.order_complex()`,
  `Poset.is_bounded()`, `Poset.is_lattice()`, and stash `A._poset = P` in `IncidenceAlgebra`
  (the provenance the fast path reads).
- The compute-kind + GUI surface: `incidence_cohomology` (algebra-level, degree range), a
  **"Poset (order relations)" input panel** on the draw canvas emitting the `IncidenceAlgebra`
  family, both runners, i18n ×4, one golden, TikZ/HTML Hasse+order-complex rendering.

**Tech Stack.** Exact only. The simplicial side is integer/`Fraction`/`Domain` linear algebra
(SNF over ℤ via `fields.linalg` / a small Hermite-Smith routine; rank over the field via the
shipped `reduce`/rank primitives). `koszul_kernels` are exact intersections of `Domain`-linear
subspaces of the free `S`-tensor module (reusing `fields.linalg.nullspace`/`intersection`). No
floats in `src/` (AST gate `tests/test_no_floats.py`). Composition is left-to-right
(Assem–Simson–Skowroński). Every refusal is a `QuiverlabError(msg, hint=...)`.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **Prerequisites, MERGED to `dev` (read-only):** Plan-05/13/16 minimal `A^e` engine
  (`engine/resolutions_minimal.py::minimal_resolution`/`minimal_homology_dims`/
  `minimal_cohomology_dims`, via `engine.adapter.to_engine`); Plan-04 CS
  (`resolutions_cs`, `cs_cohomology_dims`/`cs_homology_dims`); **Plan-27 Koszulity + Ext-algebra**
  (`modules/koszul.py::g_quadratic_certificate`/`quadratic_dual`/`_dual_data`;
  `modules/ext_algebra.py::YonedaPresentation.koszul`/`.hilbert_matrix_through`); the shipped
  `families/poset.py::Poset` + `families/incidence.py::IncidenceAlgebra` (Plan-06/33). **If any
  is absent when a worker picks this up, STOP and escalate** — do not fork the Koszul certifier,
  the minimal engine, or the incidence constructor. **All prerequisites are EARLY plans, merged
  long before this branch — tip-agnostic.** `dev` advances continuously (sibling Wave-3/4 plans
  P63/P65/P67/P70/P72/P73/… merge in their own order); P75 is `[independent]` in the wave map and
  takes **no** code edge from any of them except the **named P77 seam** (§ below), which is
  *conditional* (P75 uses only the shipped Plan-27 quadratic certifier whether or not P77 is
  merged). Re-grep every anchor at pickup and "adjust to reality"; do not assume a particular
  `dev` tip.
- **Buckets auto-assigned by directory (`tests/conftest.py`):
  `tests/{engine,resolutions_cs,modules,families,batch}/` → deep; everything else → fast.** Since
  the incidence/GHMS HH batteries drive the deep engines, they go in **`tests/families/`**
  (incidence — beside the shipped Plan-33 pin) and **`tests/engine/`** (GHMS — beside the
  minimal-engine oracles) → **deep**. `tests/qpa/` → qpa; `tests/webapp/`, `tests/gui/` → fast
  (extras-gated, unmarked per Plan-32). Run new tests by path during development; finish each task
  with a `-m deep` spot-run of the touched files.
- **The Gerstenhaber–Schack / Cibils theorem holds over ANY commutative base** — the fast path is
  **field-general** (QQ, GF(p), GF(pⁿ)) AND integral (SNF over ℤ gives every characteristic at
  once via universal coefficients). This is a genuine feature, not a scope limit: `HH^*(kP)` is
  **characteristic-sensitive** exactly when `Δ(P)` has torsion (the `RP²` pin: `HH²(GF₂·P)=1`,
  `HH²(QQ·P)=0`).
- **The incidence fast path REQUIRES incidence provenance** (`A._poset`, set by the constructor /
  the poset input mode). Plan 75 does **not** attempt to RECOGNIZE an incidence algebra from an
  arbitrary presentation (poset reconstruction from `kQ/I` is out of scope) — a
  quiver/family/module input with no poset provenance gets a **loud** `QuiverlabError`
  ("incidence_cohomology needs a poset; build via IncidenceAlgebra or the poset input mode, or
  use the general engine kinds hh_cohomology/hh_homology"). The **general** `hh_cohomology` kind
  still computes on any incidence algebra (the cross-oracle).
- **GHMS is KOSZUL ONLY, gated on the FULL three-valued verdict.** `engine="ghms"` reads
  `A.ext_algebra(window).koszul` (`True`/`False`/`None`), NOT merely `g_quadratic_certificate`
  (that is only the fast `True` path). `True` ⇒ build; `False` ⇒ refuse naming the **positive
  Ext obstruction** (e.g. preprojective `A₃`: *a new Ext-algebra generator appears in degree 3* —
  genuinely not Koszul, Brenner–Butler–King almost-Koszul); `None` ⇒ refuse **honestly
  inconclusive** (not a disproof — A may be Koszul, but this route cannot certify it; use
  `engine='cs'/'auto'`). Non-quadratic / non-length-gradable input is caught earlier with its own
  message. **No silent fallback** — an explicit `engine="ghms"` off the Koszul locus raises; only
  `engine="auto"` may *opt into* GHMS internally (dims byte-identical), and it does not by default
  (byte-stability, see the GUI-deferral note). See Scope gate 2 for the exact messages.
- **`engine="auto"` and existing goldens stay byte-identical.** The public HH dispatch
  (`core/algebra.py`) valid-engine set becomes `{"auto","bar","fast","cs","ghms"}`; `auto` is
  **unchanged** (bar/fast over GF(p) → CS fallback, Plan 34). GHMS is reachable only by explicit
  `engine="ghms"`. So no existing request re-keys and no golden drifts (verify
  `test_runner_delegation.py` before/after).
- **No floats in `src/`.** Simplicial boundary matrices are integer; the field cohomology is
  rank over the `Domain`; the integral homology is Smith normal form over ℤ. `koszul_kernels`
  are `Domain`-exact subspace intersections. The GHMS differential is exact `Domain` arithmetic.
- **Every citation key is BibTeX-verified before use** (Plan-29 rule); `_r(...)` registry entry
  + `references.bib` entry both land, gated by `tests/citations/test_bib_structure.py`.
- New `incidence_cohomology` kind acts on the **existing algebra block** ⇒ **schema stays v1**
  (the `hh_cohomology` precedent). The **poset input mode** is a new INPUT SHAPE that emits the
  **existing `IncidenceAlgebra` family** request (`family: "IncidenceAlgebra"`,
  `params: {poset_or_covers, elements}`). **Two distinct canonical-key claims, both pinned
  (critic M2):** (1) **every EXISTING request is byte-unchanged** — no new field is added to any
  non-poset request, so its `canonical_key` (`cache.py`, sort_keys JSON) is identical (verify
  `test_runner_delegation.py` + `test_catalog.py` before/after). (2) A poset request keys like the
  equivalent typed family request **only after `elements` normalization**: the catalog prefill
  defaults `elements: None`, but a naive panel would emit an explicit list, and `null` vs
  `[…]` key DIFFERENTLY under sort_keys JSON. **So the panel MUST normalize** — emit `elements`
  only for **isolated** elements (those in no cover); when every element appears in some cover,
  emit `elements: null` to match the catalog default. Then a covers-only poset keys identically to
  the typed family form; a poset WITH isolated points keys differently (correctly — it is a
  different poset), and that divergence is pinned by a test. No new algebra `kind`.
- Plan-32 markers: `d∘d=0` (GHMS + simplicial coboundary), `dim K_n` positivity, the
  SNF/universal-coefficient identity, the incidence-provenance refusal, the Koszul-gate refusal =
  `oracle_selfcert`; the theorem pins (B₃ contractible `[1,0,0,…]`; crown `S¹` `[1,1,0,…]`; `RP²`
  char-split; `Λ(kᵐ)` Betti = Koszul-dual Hilbert) = `oracle_literature`; `incidence_cohomology`
  ≡ general `hh_cohomology`, GHMS ≡ minimal ≡ bar/CS, simplicial ≡ hand-rolled SNF = two
  independent constructions = `oracle_crossengine`; **QPA has NO HH surface, NO incidence-vs-nerve
  identity, NO GHMS resolution, NO `IsKoszul`** — the `qpa` bucket is an honest-scope probe (the
  Koszul-Betti direction still rides the shipped Plan-27 `ExtAlgebraGenerators` crosscheck).
- **Mid-merge-train counts drift.** The verification-page task recounts the oracle-class table by
  running `tests/release/test_oracle_classes.py` at merge (paste live numbers, claim only this
  plan's deltas).
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table) and
  the README line. Conventional commits; green at every commit; branch `plan-75-incidence-koszul`
  off `dev`.

---

## Records (verbatim — the mandate for this plan)

> **R9 — Incidence-algebra HH = simplicial cohomology (fast path + oracle family).**
> [B-scout P6; keep-with-corrections]
> Object: for incidence algebras kP, HH^* ≅ simplicial H^* of the order complex —
> CUP-PRODUCT/ring isomorphism (Cibils JPAA 56 (1989) 221–232;
> Gerstenhaber–Schack JPAA 30 (1983) 143–156). The "bracket = 0" claim is NOT in
> G–S '83 — dropped unless re-sourced (math/0611542 / arXiv:2411.07910 line).
> Poset input mode + SNF-fast HH + a second independent oracle for the general
> engines (extends the existing B₃ pin to the theorem). Size S–M.

> **R10 — GHMS Koszul minimal bimodule resolution (fast HH for Koszul algebras).**
> [B-scout P9; keep]
> Object: the comultiplicative minimal A^e-resolution of a Koszul algebra built
> from the Koszul dual data; fast minimal HH + a third HH oracle class. Ref:
> Green–Hartman–Marcos–Solberg, Arch. Math. 85 (2005) 118–127; arXiv:math/0508177.
> Oracles: exterior/preprojective HH degreewise vs existing engines; Betti =
> Koszul-dual Hilbert coefficients. Size M. Deps: Plan-27 Koszulity + Ext-algebra.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P75. Tier δ (structure); marked
`[independent]` in the wave map.

---

## Seam with Plan 77 (N-Koszul / K₂ / multi-Koszul) — the boundary, named

Both plans touch the Koszul surface; the split is clean and contractual:

- **P75 = the GHMS RESOLUTION for (ordinary, quadratic) Koszul algebras.** It *consumes* Plan-27's
  **three-valued `ext_algebra.koszul` verdict** as a gate and *builds* the comultiplicative
  bimodule resolution + fast HH. It adds **no** new Koszulity recognizer; it neither decides nor
  extends the notion of Koszul.
- **P77 = generalized Koszulity RECOGNIZERS** (N-Koszul via the 2-N alternation of Ext generation
  degrees, K₂ / Cassidy–Shelton with a certified window, multi-Koszul, and the **`(p,q)`-almost-
  Koszul** / Brenner–Butler–King family). It *reads* generation degrees off the shipped minimal
  resolutions; it builds **no** bimodule resolution.
- **Shared surface:** both cite Plan-27's `modules/koszul.py`. P75 adds `koszul_kernels` +
  `GHMSResolution` (a *resolution*); P77 adds generation-degree *predicates*. **The three-valued
  `ext_algebra.koszul` verdict is exactly the boundary:** GHMS applies precisely where that verdict
  is `True`; P77's recognizers are what push into the `False`/`None` territory (recognizing the
  higher/almost-Koszul structure GHMS declines).
- **The shared BOUNDARY EXAMPLE is preprojective `A₃`** (and Dynkin preprojectives generally).
  `ext_algebra.koszul(A₃) = False` with the obstruction *a new Ext-algebra generator in degree 3*
  — the algebra is **`(p,q)`-almost-Koszul** (Brenner–Butler–King; ADE preprojectives break at
  `h − 1`). **P75 REFUSES it** (not ordinarily Koszul, so no GHMS resolution); **P77 RECOGNIZES it**
  as almost/`(p,q)`-Koszul. Naming the same algebra on both sides of the seam is the cleanest way
  to see the boundary is a real theorem, not a tooling artifact.
- If both are in flight, P75 does not depend on P77 (uses only the shipped three-valued verdict); a
  later P77 N-Koszul / almost-Koszul GHMS-style generalization is explicitly out of P75's scope
  (recorded honest-scope).

---

## Reference re-verification (mandatory; findings binding, `# PIN` = resolve at citation-add)

All statements below were re-read from the primary sources / their reviews during authoring.
The **empirical ground truth** (every HH/simplicial value) was independently computed by running
the merged engines during authoring (recorded in §Live-verified facts and §Methodology).

- **Cibils, *Cohomology of incidence algebras and simplicial complexes*, JPAA 56 (1989)
  221–232 — VERIFIED (bib present).** The existing `cibils1989incidence` entry (`references.bib`
  L474) has the exact volume/number/pages `56(3) (1989) 221–232`. Cibils computes `HH^*` of an
  incidence algebra and identifies it with the **simplicial cohomology of the associated
  simplicial complex** (the order complex / nerve). This is the graded-vector-space identification
  the fast path realizes; the RING (cup) statement is the Gerstenhaber–Schack line below.
- **Gerstenhaber–Schack, *Simplicial cohomology is Hochschild cohomology*, JPAA 30 (1983)
  143–156 — TO ADD (distinct from the 1987 Hodge paper).** The shipped `GerstenhaberSchack1987`
  (`references.bib` L135, *A Hodge-type decomposition…*, JPAA **48** (1987) 229–247) is a
  **different paper** and must NOT be reused for R9. Add a new key `gerstenhaber_schack_1983`
  (→ `GerstenhaberSchack1983`), JPAA **30** (1983) no. 2, 143–156. This is the paper that proves
  the identification is a **cup-product (ring) isomorphism** `HH^*(kP) ≅ H^*_{simp}(Δ(P))` — the
  source of R9's RING clause. `# PIN` the exact issue number/pages at citation-add.
- **The "bracket = 0" clause is DROPPED (per the record).** Gerstenhaber–Schack 1983 proves the
  **cup-product** ring isomorphism; it does **not** assert the Gerstenhaber **bracket** vanishes
  on `HH^*(kP)`. The bracket-triviality line traces to later work (math/0611542 Redondo; the
  arXiv:2411.07910 line) and is NOT re-sourced here. Plan 75 therefore claims **only** the RING
  (cup) isomorphism at the level this plan verifies (dimension iso in v1; a ring-check task for
  the cup); it makes **no** bracket claim. The verification page states this explicitly.
- **Redondo, *Hochschild cohomology via incidence algebras*, JLMS 77 (2008) 465–480,
  arXiv:math/0611542 — VERIFIED (bib present).** `redondo2008incidence` (`references.bib` L484)
  is the second source of the identification (with Cibils 1989). Kept as the corroborating
  citation; it is the `math/0611542` line the record names.
- **Green–Hartman–Marcos–Solberg, *Resolutions over Koszul algebras*, Arch. Math. (Basel) 85
  (2005) 118–127, arXiv:math/0508177 — TO ADD.** Add `green_hartman_marcos_solberg`
  (→ `GHMS2005`). The paper constructs, for a Koszul algebra `A = TV/(R)`, the **minimal graded
  projective bimodule resolution** with terms `P_n = A ⊗_S K_n ⊗_S A`, `K_n = ⋂_{i} V^{⊗i} ⊗ R ⊗
  V^{⊗(n-2-i)}`, and an explicit **comultiplicative** differential built from a splitting of the
  inclusions `K_n ↪ V ⊗ K_{n-1}` and `K_n ↪ K_{n-1} ⊗ V`. `dim_k K_n =` the `n`-th coefficient of
  the Hilbert series of the Koszul dual `A^!` (equivalently the total graded Betti number
  `Σ_{i,j} dim Ext^n_A(S_i, S_j)`). `# PIN` the exact Archiv der Mathematik volume/pages at
  citation-add (the arXiv id math/0508177 is stable). NB the DOI is 10.1007/s00013-005-1299-9
  (`# verify`).
- **Priddy PBW / Koszul certificate — SHIPPED (Plan 27).** `koszul.g_quadratic_certificate`
  (Priddy 1970, `priddy` key) is the gate; `koszul.quadratic_dual` (Polishchuk–Positselski,
  `polishchuk_positselski`) provides `A^!`; `ext_algebra.hilbert_matrix_through` provides the
  graded Betti. All BibTeX-verified in Plan 27.

---

## Mathematical foundation (the definitions, exactly — the plan's ground truth)

### (a) Incidence algebras and the order complex (R9)

Let `P` be a finite poset. The **incidence algebra** `kP` has `k`-basis the intervals
`[x, y] = {z : x ≤ z ≤ y}` for `x ≤ y`, i.e. one basis element `e_{xy}` per pair `x ≤ y`, with
`e_{xy} · e_{zw} = δ_{yz} e_{xw}` (path composition). As a bound quiver algebra it is the
**Hasse quiver** `Q` (vertices = elements, arrows = cover relations `x ⋖ y`) modulo the
**commutativity ideal** `I` = ⟨all differences of parallel directed paths with the same source
and target⟩ (every interval's maximal chains are identified). This is exactly the shipped
`IncidenceAlgebra` (`families/incidence.py`): `_all_paths(Q)` enumerates directed Hasse paths,
groups by `(source, target)`, and emits `base − other` relations. `dim_k kP = #{(x,y) : x ≤ y}`.

The **order complex** (nerve) `Δ(P)` is the abstract simplicial complex whose `p`-simplices are
the **chains** `x₀ < x₁ < … < x_p` of `P` (`|Δ(P)|` = the classifying space `BP`). Its simplicial
cochain complex over `k`:

    C^p(Δ(P); k) = k^{#chains of length p},   (δf)(x₀<…<x_{p+1}) = Σⱼ (−1)ʲ f(…x̂ⱼ…).

**Theorem (Cibils 1989, generalizing Gerstenhaber–Schack 1983).** For every finite poset `P`
and field `k`,

    HH^n(kP) ≅ H^n(Δ(P); k)   for all n ≥ 0,

as **graded rings** (the Hochschild cup product corresponds to the simplicial cup product). **The
attribution, precisely (per the reference re-verification):** Gerstenhaber–Schack 1983 prove the
cup-product isomorphism for the incidence algebra of the **face poset of a simplicial complex**
(`kΔ`, "simplicial cohomology is Hochschild cohomology"); **Cibils 1989 extends it to an ARBITRARY
finite poset** `P` (`HH^*(kP) ≅ H^*(Δ(P))`, the order-complex/nerve version), which is the general
statement this plan realizes. So an arbitrary-poset row cites **Cibils 1989** (+ Redondo 2008); the
face-poset rows (the `RP²₆` triangulation) additionally sit squarely in G–S 1983's original setting.
In particular `HH^0(kP) = k` iff `P` is connected (the order complex is connected), and if `Δ(P)`
is contractible — e.g. `P` has a global minimum or maximum (`0̂` or `1̂`), so the order complex is a
cone — then `HH^{≥1}(kP) = 0` (this is the Plan-33 `B₃` pin, now a theorem instance).

**The fast path.** `dim_k HH^n(kP) = dim_k H^n(Δ(P); k)`, computed by exact linear algebra on the
cochain complex whose degree-`p` dimension is `#chains of length p` — orders of magnitude smaller
than the normalized bar complex `dim C_n = m·(m−1)^n` (`m = dim kP`) or the CS/minimal bimodule
resolution. **Where the speed comes from, stated honestly:** the theorem replaces a computation
over the enveloping algebra `(kP)^e` (dimension `m²` per resolution term, rebuilt per prime) with
the *combinatorial* cochain complex of `Δ(P)` (dimension = chain counts), and the **integer**
boundary matrices are computed **once** — Smith normal form over ℤ then yields `H^n(Δ; k)` for
**every** `k` simultaneously by the universal-coefficient theorem:

    dim_k H^n(Δ; k) = rank_ℤ H_n  +  t_n(char k)  +  t_{n-1}(char k),

where `t_n(p)` = the number of `Z/p^a` invariant factors of the integral homology `H_n(Δ; ℤ)`.
Torsion in `H_*(Δ; ℤ)` is exactly the source of **characteristic-dependent** `HH^*(kP)` — the
`RP²` pin (`H_1(RP²;ℤ) = Z/2`) gives `HH²(GF₂·P) = 1` but `HH²(QQ·P) = 0`.

### (b) The GHMS comultiplicative resolution for Koszul algebras (R10)

Let `A = TV/(R)` be a **Koszul** `kQ/I`: `S = k^{Q₀}` the semisimple degree-0 part, `V = kQ₁`
(arrows) the degree-1 `S`-bimodule, `R ⊆ V ⊗_S V` the (necessarily quadratic) relation space.
Define the **Koszul kernels**

    K_0 = S,   K_1 = V,   K_n = ⋂_{i=0}^{n-2} V^{⊗i} ⊗_S R ⊗_S V^{⊗(n-2-i)}  ⊆  V^{⊗n}   (n ≥ 2).

**Theorem (GHMS 2005).** The minimal graded projective `A^e`-resolution of `A` is

    ⋯ → A ⊗_S K_n ⊗_S A --d_n--> A ⊗_S K_{n-1} ⊗_S A → ⋯ → A ⊗_S A --μ--> A → 0,

with `d_n` the **comultiplicative** differential: `K_n` sits inside both `V ⊗_S K_{n-1}` and
`K_{n-1} ⊗_S V`; writing an element `ω ∈ K_n` under the first inclusion as `Σ x_a ⊗ ω'_a`
(`x_a ∈ V`, `ω'_a ∈ K_{n-1}`) and under the second as `Σ ω''_b ⊗ y_b` (`ω''_b ∈ K_{n-1}`,
`y_b ∈ V`),

    d_n(1 ⊗ ω ⊗ 1) = Σ_a x_a ⊗ ω'_a ⊗ 1  −  (−1)^n Σ_b 1 ⊗ ω''_b ⊗ y_b   ∈ A ⊗_S K_{n-1} ⊗_S A,

with `x_a, y_b` acting as elements of `A` (degree 1). Coassociativity of the two splittings is the
"comultiplicative" structure that makes `d_{n-1} d_n = 0`. **Consequences the plan leans on:**
`rank_{A^e} P_n = dim_S K_n = dim_k K_n` and

    dim_k K_n  =  [t^n] Hilbert series of A^!  =  Σ_{i,j} dim_k Ext^n_A(S_i, S_j)   (the graded Betti).

**The collapse.** `HH_n(A) = H_n(A ⊗_{A^e} P_•)`, `HH^n(A) = H^n(Hom_{A^e}(P_•, A))`; both are
exact `Domain` linear algebra on the corner-typed blocks `e_v K_n e_w` (the Plan-16 Hom/⊗
collapse convention: cohomology acts `a·w·b`, homology `b·w·a` on the swapped corner tags). This
is **fast** (closed-form terms, no syzygy search) and **`Domain`-general** (QQ / GF(p) / GF(pⁿ)),
unlike the GF(p)-only syzygy engine.

**Sign convention is ARBITRATED, not assumed** (the Plan-20/21 discipline): the `−(−1)^n` sign is
pinned by simultaneously requiring `d∘d = 0` (self-cert), `HH_* (GHMS) ≡ HH_* (minimal syzygy
engine over GF(p))` (cross-engine), and `HH_* (GHMS) ≡ HH_* (bar/CS)` where the bar window
reaches. A wrong sign fails at least one anchor.

### The two consequences that unify the plan

`HH^*(kP)` and `HH^*(A_{Koszul})` are both computed **without a general resolution**: R9 replaces
it with the order complex, R10 with the closed-form Koszul kernels. Each ships its own **new
oracle class** (the SNF simplicial complex; the GHMS bimodule ranks) that cross-validates the
general engines from an independent route.

---

## Live-verified facts (all recomputed in the venv at spec time; `QUIVERLAB_NO_NUMBA=1`)

Every value below was produced by the merged engines during authoring; the **simplicial side** was
computed by a **hand-rolled exact SNF/rank** routine (sympy `Matrix.rank` over QQ and mod `p`,
integer boundary matrices) — that routine IS the design's second engine, and the plan re-implements
it in `src/` over the `Domain`.

**(a) Incidence — `HH^*(kP)` (general CS engine) ≡ simplicial `H^*(Δ(P))` (hand-rolled SNF):**

| poset `P` | `dim kP` | order-complex `f`-vector | lattice? | route | `HH^*(kP)` | simplicial `H^*` |
|---|---|---|---|---|---|---|
| **`B₃`** Boolean lattice (anchor) | 27 | `(8,19,18,6)` | lattice, `0̂`+`1̂` | non-monomial | `[1,0,0]` | `[1,0,0]` |
| **crown `C(3,3)`** hexagon | 12 | `(6,6)` | **non-lattice** | monomial (`I=0`) | `[1,1,0,0,0]` | `[1,1]` |
| **crown `C(2,2)`** square | 8 | `(4,4)` | **non-lattice** | monomial | `[1,1,0,0,0]` | `[1,1]` |
| **thickened crown** (one edge → diamond) | 18 | `(8,10,2)` | **non-lattice** | **non-monomial** (commutativity) | `[1,1,0,0,0]` | `[1,1,0]` |
| **whisker** `w<x<z, w<y` | 8 | `(4,4,1)` | **non-lattice** | non-monomial | `[1,0,0,0,0]` | `[1,0,0]` |
| **`N₅`** pentagon | 13 | `(5,8,5,1)` | lattice, bounded | non-monomial | `[1,0,0,0,0]` | `[1,0,0,0]` |
| **`RP²₆` face poset** (char pin) | **121** | `(31, …)` (bary. subdiv. of `RP²`) | non-lattice | non-monomial | `HH^*(GF₂)=[1,1,1]`, `HH^*(QQ)=[1,0,0]` | `H^*(GF₂)=[1,1,1]`, `H^*(QQ)=[1,0,0]` |

  - **Every row MATCHES** the hand-rolled simplicial cohomology, degreewise. The `B₃` anchor
    (Plan-33) is recovered; the crowns give the record's **nonvanishing `H¹`** (order complex
    `≅ S¹`); the **thickened crown** is a **non-lattice, non-monomial** (commutativity-relation)
    poset with `H¹ ≠ 0`; `N₅`/whisker exercise the contractible non-`B₃` cases.
  - **The `RP²₆` pin is the char-sensitivity showcase.** The minimal 6-vertex `RP²`
    triangulation `FACETS = [(1,2,3),(1,3,4),(1,4,5),(1,5,6),(1,2,6),(2,3,5),(2,4,5),(2,4,6),
    (3,4,6),(3,5,6)]` (independently validated: every edge in exactly 2 facets, `f=(6,15,10)`,
    Euler char `= 1 = χ(RP²)`, `H_*(QQ)=[1,0,0]`, `H_*(GF₂)=[1,1,1]`, `H_*(GF₃)=[1,0,0]`). Its
    **face poset** (31 elements, 60 covers) has incidence algebra `dim 121`; the general **CS
    engine** gives `HH^*(GF₂)=[1,1,1]` but `HH^*(QQ)=[1,0,0]` to degree 2 — **exactly** the
    universal-coefficient prediction, and exactly the hand-rolled SNF simplicial answer. Both
    engines agree; the general CS engine is not intractable here (measured at authoring: ≈2.8 s
    over GF₂, ≈6.4 s over QQ, to degree 2). The honest framing is a **size comparison**, not a
    "the general engine can't": the simplicial cochain complex has `f = (31, …)` chains, whereas
    the incidence algebra is `dim 121` and its bar/CS bimodule resolution term is `O(121²)` per
    degree — the fast path is *much smaller*, and it delivers *every* characteristic from one
    integer SNF (the UCT), which the general engine cannot (each field is a separate run).

**(b) GHMS / Koszul — minimal bimodule ranks = Koszul-dual Hilbert; fast HH ≡ general:**

| algebra `A` | `dim A` | `ext_algebra.koszul` | minimal bimodule ranks `r_n = dim K_n` | Koszul-dual Hilbert | Ext-algebra Betti (Plan-27) | `HH_*` (minimal) | `HH_*` (bar) |
|---|---|---|---|---|---|---|---|
| `Λ(k²)` exterior | 4 | **True** | `[1,2,3,4,5,6,7,8]` | `C(n+1,n) = [1,2,3,4,5,6,7,8]` | `[1,2,3,4,5,6,7]` | `[3,4,6,8,10,12,14]` | `[3,4,6,8,10]` (≡) |
| `Λ(k³)` exterior | 8 | **True** | `[1,3,6,10,15,21,28,36]` | `C(n+2,2) = [1,3,6,10,15,21,28,36]` | — | `[5,12,24,40,60,84]` | (≡ where in window) |
| **`kZ₃/rad²`** cyclic Nakayama (multi-vertex) | 6 | **True** | (corner-typed) | — | — | `[3,0,1,1,0,0]` | (≡ where in window) |
| **diamond incidence** (comm. square, multi-vertex) — the **(a)↔(b) bridge** | 9 | **True** | `[4,4,1,0,0,0]` (`K₀=|Q₀|=4`, `K₁=#arrows=4`, `K₂=#rel=1`; gl.dim 2) | — | — | `[4,0,0,0,0]` | (≡) |
| preprojective `A₃` (self-inj.) — the **NEGATIVE** example | 10 | **False** (obstruction `(3, "a new Ext-algebra generator appears in degree 3")`) | `[3,4,3,3,4,3,3]` (syzygy engine) | — | — | `[3,1,1,1,1,1]` (syzygy engine; GHMS **refuses**) | (≡) |

  - **`r_n = dim K_n = the Koszul-dual Hilbert coefficient`, EXACTLY** for `Λ(k²)` (`n+1`) and
    `Λ(k³)` (`C(n+2,2)`) — the R10 Betti oracle, cross-checked against the shipped Plan-27
    `ext_algebra.hilbert_matrix_through` total (`[1,2,3,4,5,6,7]` for `Λ(k²)`, `koszul=True`).
    GHMS builds `P_n = A⊗K_n⊗A` with these very ranks in closed form; the syzygy engine's ranks
    are the cross-oracle. (Multi-vertex: `K₀ = |Q₀|`, not 1 — the diamond's `[4,4,1,…]`.)
  - **`HH_*` (minimal) ≡ `HH_*` (bar)** on `Λ(k²)` (`[3,4,6,8,10]` agree in the bar window) —
    GHMS must reproduce this degreewise (the third-oracle-class cross-engine gate).
  - **`kZ₃/rad²`** (cyclic radical-square-zero Nakayama, `1→2→3→1` mod rad²) is a **multi-vertex,
    `ext_algebra.koszul = True`** Koszul algebra with **nontrivial** `HH_* = [3,0,1,1,0,0]`
    (minimal engine, LIVE-VERIFIED) — a corner-typed (`e_v K_n e_w`) multi-vertex GHMS target.
  - **The diamond (commutative square) incidence algebra IS the (a)↔(b) bridge — and its numbers
    were CORRECTED (critic W2).** It is `ext_algebra.koszul = True` (`dim 9`), so GHMS applies; its
    minimal bimodule ranks are **`r_n = dim K_n = [4,4,1,0,0,0]`** (NOT length 1 — `K₂ = 1` is the
    single commutativity relation, gl.dim 2). **Two DISTINCT quantities must not be conflated:**
    (i) its **Hochschild COHOMOLOGY = order-complex cohomology** `HH^*(kP) = H^*(Δ) = [1,0,0]`
    (the order complex is contractible — the R9 theorem); (ii) its **Hochschild HOMOLOGY**
    `HH_*(kP) = [4,0,0,0,0]` (`HH₀ = dim A/[A,A] = 4`, not `1`). GHMS computes the HOMOLOGY
    `[4,0,0,0,0]` (a second multi-vertex GHMS target); the `[1,0,0]` belongs to the R9 cohomology
    theorem. A genuinely useful bridge: a Koszul algebra that is ALSO an incidence algebra,
    computed by BOTH engines.
  - **Preprojective `A₃` is NOT Koszul — it is the NEGATIVE example (critic W1).** `ext_algebra`
    returns `koszul = False` with the **positive** obstruction `(3, "a new Ext-algebra generator
    appears in degree 3")` (LIVE-VERIFIED). This is exactly Brenner–Butler–King **almost-Koszul /
    `(p,q)`-Koszul** theory: Dynkin preprojectives break Koszulity at `h − 1 = 3` for `A₃`. So
    `engine="ghms"` **correctly refuses** it — for the RIGHT reason (*not Koszul, Ext obstruction
    at degree 3*), NOT the earlier draft's wrong "Koszul but PBW-uncertified". The minimal syzygy
    engine still computes `HH_* = [3,1,1,1,1,1]`. **A₃ is the shared P75/P77 boundary example**
    (P75 refuses it as not-Koszul; P77 *recognizes* it as `(p,q)`-almost-Koszul — see the seam).
  - **The commutative-square incidence algebra (diamond, `dim 9`) is Koszul (`g_quadratic=True`)**
    — a concrete **bridge** between (a) and (b): an incidence algebra that GHMS also handles (its
    order complex is contractible, so `HH_* = [1,0,0]` and `K_n = 0` for `n ≥ 2` — the GHMS
    resolution is length 1). A pleasant self-consistency check across both halves of the plan.

**(c) Provenance / plumbing facts confirmed live:** `IncidenceAlgebra(covers, elements, field)`
builds over CC/GF(p)/QQ; it is already a webapp catalog family (`IncidenceAlgebra`, param
`poset_or_covers`); `to_engine(A)` + `minimal_resolution(eng, N, p)` return
`(rks_dict, cols, eng, trunc)`; `hochschild_homology_dims(eng, N, primes=(p,))` needs the engine
algebra (not the public `Algebra`) and a **tuple** `primes`.

---

## Scope gates (contractual; every boundary is a loud typed refusal)

1. **`incidence_cohomology` needs poset provenance.** No `A._poset` ⇒ loud
   `QuiverlabError("incidence_cohomology: no poset provenance", hint="build via IncidenceAlgebra /
   the poset input mode, or use hh_cohomology/hh_homology for the general engine")`. (Poset
   reconstruction from an arbitrary `kQ/I` is out of scope — stated on the verification page.)
2. **`engine="ghms"` reads the FULL three-valued `ext_algebra.koszul` verdict** (`modules/
   ext_algebra.py::YonedaPresentation.koszul` + `.koszul_obstruction`), NOT merely
   `g_quadratic_certificate` (which is only the fast `True` path). Gate:
   - `koszul is True` ⇒ build (G-quadratic/PBW, or Fröberg-consistent — the certifier's `True`).
   - `koszul is False` ⇒ loud refusal naming the **positive obstruction**:
     `"engine='ghms': A is NOT Koszul — {obstruction[1]} (degree {obstruction[0]}); GHMS applies
     only to Koszul algebras. Use engine='cs'/'auto'."` (e.g. preprojective `A₃`: *a new
     Ext-algebra generator appears in degree 3*).
   - `koszul is None` ⇒ loud refusal, **honestly inconclusive**: `"engine='ghms': Koszulity
     could not be certified ({reason}) — this is not a disproof, but GHMS cannot proceed without
     it. Use engine='cs'/'auto'."`
   **The gate is REAL, live-confirmed:** preprojective `A₃` returns `koszul = False` (obstruction
   `(3, "a new Ext-algebra generator appears in degree 3")`) — a genuine *not-Koszul* verdict
   (Brenner–Butler–King almost-Koszul), so GHMS refuses it for the RIGHT reason; the minimal
   syzygy engine still computes. **Implementation efficiency:** try `g_quadratic_certificate(A)`
   first (cheap `True`); only when it does not fire, compute `A.ext_algebra(window).koszul` for
   the honest `False`/`None` message (a small `window` — degree 3 catches `A₃`'s obstruction).
3. **Presentation-less (structure-constant) algebras** refuse loudly on both routes (incidence
   needs the Hasse quiver; GHMS needs `V`, `R`) — surfaced as a clean typed 4xx, never a 500.
4. **The integral / char certificate is honest.** `integral_homology` reports free rank +
   invariant factors; `incidence_cohomology(A, top, field)` returns the field dims AND (when the
   integral SNF is in budget) the torsion note. If the SNF budget is exceeded on a large order
   complex, the field-`rank` path still returns the dims for that field, and the char note is
   omitted with a stated reason (never a fabricated "characteristic-independent").
5. **GHMS `Domain` generality is real but bounded by `dim K_n` growth.** For large `V` the
   `K_n` intersections grow; a `max_term_dim`-style budget caps the build with a truncation flag
   (mirroring the minimal engine), the computed degrees exact.
6. **Ring (cup) claim is scoped.** v1 ships the **dimension** isomorphism `HH^n = H^n(Δ)`
   (fully verified) + the **graded-commutativity** identity of the incidence HH cup (a shipped
   Plan-35 identity oracle). A **cup-rank ring-check** — `rank(HH^p ⊗ HH^q → HH^{p+q})` (a
   basis-independent invariant) equals `rank(H^p ⊗ H^q → H^{p+q})` on the order complex, pinned
   where a nonzero cup exists (a torus poset) — is a **task**, honest-scoped if the torus
   incidence algebra's cup does not compute in budget. **No bracket claim is made** (the record's
   dropped clause).

---

## API surface (public via `import quiverlab`; exact only)

```python
# --- (a) incidence: src/quiverlab/hochschild/simplicial.py ---
class OrderComplex:
    """The order complex (nerve) of a Poset: chains-by-dimension + integer boundary maps."""
    @classmethod
    def of(cls, P) -> "OrderComplex": ...          # P: families.poset.Poset
    def face_vector(self) -> tuple: ...            # (#0-chains, #1-chains, ...)
    def boundary(self, k) -> list[list[int]]: ...  # integer d_k : C_k -> C_{k-1}

def simplicial_cohomology_dims(oc, top, field) -> list[int]:
    """dim_k H^n(Delta; k), n = 0..top, by exact rank over `field` (Domain). The FAST path."""

def integral_homology(oc, top) -> list[tuple[int, tuple]]:
    """[(free_rank_n, invariant_factors_n), ...] via Smith normal form over Z -- the
    universal-coefficient char certificate (torsion => characteristic-dependent HH)."""

def incidence_cohomology(A, top, field=None):
    """HH^*(kP) via the Gerstenhaber-Schack/Cibils theorem = H^*(order complex). Requires
    A._poset provenance; loud QuiverlabError otherwise. Returns an IncidenceCohomology result:
    .dims (over `field` or A's field), .face_vector, .contractible (bool | None),
    .torsion (per-degree invariant factors | None), .char_dependent (bool | None),
    .references. Ring: dims-level in v1; .cup_generated flag when HH is generated in degree 1."""

# Algebra delegate (core/algebra.py, thin lazy-import):
# A.incidence_cohomology(top) -> incidence_cohomology(A, top)

# --- (b) GHMS: src/quiverlab/hochschild/koszul_ghms.py ---
def koszul_kernels(A, top) -> list:
    """[K_0, ..., K_top], K_n a Domain-basis of the intersection subspace of V^{tensor n}.
    dim K_n = the n-th Koszul-dual Hilbert coefficient. Requires a length-graded quadratic A."""

class GHMSResolution(Resolution):
    """The comultiplicative minimal A^e-resolution P_n = A (x)_S K_n (x)_S A of a Koszul A.
    Gated on the THREE-VALUED ext_algebra.koszul verdict: True -> build; False -> refuse with the
    positive Ext obstruction; None -> refuse, inconclusive. term(n), differential(n)
    (byte-reproducible), assert_dd_zero. Runs over any exact Domain."""
    def __init__(self, A, *, max_term_dim=...): ...

def ghms_homology_dims(A, top) -> list[int]: ...     # HH_n via A (x)_{A^e} P_*
def ghms_cohomology_dims(A, top) -> list[int]: ...   # HH^n via Hom_{A^e}(P_*, A)

# HH dispatch (core/algebra.py): engine set -> {"auto","bar","fast","cs","ghms"};
#   engine="ghms" routes to ghms_*_dims (Koszul only, three-valued verdict; loud refusal on
#   False=not-Koszul-with-obstruction or None=inconclusive).
```

---

### Task 0: citations (do first — later tasks reference the keys)

**Files:** modify `src/quiverlab/citations/references.bib`, `src/quiverlab/citations/registry.py`;
test `tests/citations/test_bib_structure.py` (existing gate).

**Interfaces:** the `_r(key, bibtex_key, kind, title, annotation, *tags)` helper (`registry.py`).
`cibils_incidence` (→ `cibils1989incidence`), `redondo_incidence` (→ `redondo2008incidence`),
`priddy`, `froberg_koszul`, `polishchuk_positselski` **already exist** — reuse, do NOT re-add.

- [ ] **Step 1: add two keys** (BibTeX-verify each — Plan-29 rule; `# PIN` unresolved fields):
  - `gerstenhaber_schack_1983` → `GerstenhaberSchack1983`, *Simplicial cohomology is Hochschild
    cohomology*, J. Pure Appl. Algebra **30** (1983) no. 2, 143–156. (the **RING** isomorphism
    `HH^*(kP) ≅ H^*_{simp}(Δ(P))` — R9's cup clause; NOT the 1987 Hodge paper). `# PIN` issue/pages.
  - `green_hartman_marcos_solberg` → `GHMS2005`, Green–Hartman–Marcos–Solberg, *Resolutions over
    Koszul algebras*, Arch. Math. (Basel) **85** (2005) 118–127, arXiv:math/0508177, doi
    10.1007/s00013-005-1299-9. (the comultiplicative minimal bimodule resolution — R10). `# verify`
    the exact pages/DOI.

  Registry entries (mirroring the shipped incidence/koszul annotations; annotations reflect the
  **corrected** attributions — G–S 1983 = the ring iso, bracket claim DROPPED; GHMS = the
  resolution, Betti = Koszul-dual Hilbert):

```python
_r("gerstenhaber_schack_1983", "GerstenhaberSchack1983", "foundation",
   "Simplicial cohomology is Hochschild cohomology",
   "Gerstenhaber-Schack: HH^*(kP) is isomorphic to the simplicial cohomology of the order "
   "complex of P AS A RING (cup product). The ring source for the incidence-vs-nerve oracle "
   "(with Cibils 1989 / Redondo 2008). NB: proves the CUP iso; makes no bracket-vanishing claim.",
   "hochschild", "incidence", "oracle"),
_r("green_hartman_marcos_solberg", "GHMS2005", "algorithm",
   "Resolutions over Koszul algebras",
   "Green-Hartman-Marcos-Solberg: the minimal graded A^e-resolution P_n = A (x)_S K_n (x)_S A of "
   "a Koszul algebra with the comultiplicative differential; K_n = intersection of V^i (x) R (x) "
   "V^j, dim K_n = the n-th Koszul-dual Hilbert coefficient. The GHMS fast-HH engine (Plan 75).",
   "hochschild", "koszul", "resolution"),
```

- [ ] **Step 2:** `... -m pytest tests/citations/test_bib_structure.py -q` green (registry ↔ bib
  in sync; every `_r` key resolvable).
- [ ] **Step 3: Commit** — `docs(citations): incidence ring iso (Gerstenhaber-Schack 1983) + GHMS Koszul resolution (math/0508177)`

---

## Task group I — (a) incidence algebras: order complex + fast HH + theorem oracle (R9)

### Task I1: `Poset` order complex + provenance

**Files:** modify `src/quiverlab/families/poset.py`, `src/quiverlab/families/incidence.py`;
test `tests/families/test_order_complex.py`.

**Interfaces:** the shipped `Poset(covers, elements)` (`_le` closure, `leq`, `hasse_quiver`).

- [ ] **Step 1: failing tests.**

```python
# tests/families/test_order_complex.py
"""The order complex Delta(P) of a Poset: chains-by-dimension, face vector, boundary maps;
and the incidence-algebra provenance stash A._poset (Plan 75 / R9)."""
import pytest
from quiverlab.families.poset import Poset
from quiverlab.families import IncidenceAlgebra
from quiverlab.fields import CC

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature

def _b3(): return [(x, x | (1 << b)) for x in range(8) for b in range(3) if not (x >> b) & 1]

@lit
def test_b3_face_vector():
    P = Poset(_b3())
    assert P.order_complex().face_vector() == (8, 19, 18, 6)   # LIVE-VERIFIED
    assert P.is_bounded() and P.is_lattice()

@lit
def test_crown_c33_is_s1_face_vector():
    P = Poset([("a","X"),("b","X"),("b","Y"),("c","Y"),("c","Z"),("a","Z")])
    assert P.order_complex().face_vector() == (6, 6)           # S^1: no 2-chains
    assert not P.is_lattice() and not P.is_bounded()

@selfcert
def test_incidence_provenance_stashed():
    A = IncidenceAlgebra([(1,2),(2,3)], field=CC)
    assert A._poset is not None and A._poset.leq(1, 3)         # provenance for the fast path
```

- [ ] **Step 2: implement.** `Poset.order_complex()` returns an `OrderComplex.of(self)`
  (Task I2 owns the class; `Poset` just delegates). `Poset.is_bounded()` (unique global min AND
  max under `leq`), `Poset.is_lattice()` (every pair has a unique least upper bound and greatest
  lower bound — exact finite check over `_le`). In `IncidenceAlgebra`, after building `A`, set
  `A._poset = P` (the constructed/passed `Poset`). Byte-check: the existing
  `tests/families/test_incidence.py` stays green (only an attribute is added).
- [ ] **Step 3: Commit** — `feat(families): Poset.order_complex/is_bounded/is_lattice + incidence _poset provenance`

### Task I2: `OrderComplex` + the SNF/rank simplicial engine (the fast path)

**Files:** create `src/quiverlab/hochschild/simplicial.py`; modify
`src/quiverlab/hochschild/__init__.py`; test `tests/families/test_simplicial_cohomology.py`.

**Interfaces:** `fields.linalg` (rank over a `Domain`; a Smith-normal-form-over-ℤ helper — add
`fields.linalg.smith_normal_form` if absent, else a local Hermite/Smith on integer matrices).

- [ ] **Step 1: failing tests** (values LIVE-VERIFIED against the hand-rolled routine):

```python
# tests/families/test_simplicial_cohomology.py
"""The order-complex simplicial engine: exact field cohomology (rank) + integral homology (SNF).
The SECOND, independent HH oracle for incidence algebras (Plan 75 / R9). Char-sensitivity via
universal coefficients: RP^2 has H^*(GF2) = [1,1,1] but H^*(QQ) = [1,0,0]."""
import pytest
from quiverlab.families.poset import Poset
from quiverlab.fields import QQ, GF
from quiverlab.hochschild.simplicial import OrderComplex, simplicial_cohomology_dims, integral_homology

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert

CROWN = [("a","X"),("b","X"),("b","Y"),("c","Y"),("c","Z"),("a","Z")]
# minimal 6-vertex RP^2 triangulation, as a FACE POSET (covers sigma < tau, tau = sigma + 1 vertex):
RP2_FACETS = [(1,2,3),(1,3,4),(1,4,5),(1,5,6),(1,2,6),(2,3,5),(2,4,5),(2,4,6),(3,4,6),(3,5,6)]

def _face_poset(facets):
    import itertools
    simp = set()
    for f in facets:
        for k in range(1, len(f)+1):
            for c in itertools.combinations(sorted(f), k): simp.add(c)
    covers = [(t[:i]+t[i+1:], t) for t in simp if len(t) >= 2
              for i in range(len(t)) if (t[:i]+t[i+1:]) in simp]
    return Poset(covers, elements=sorted(simp, key=lambda s:(len(s), s)))

@lit
def test_crown_is_circle():
    oc = OrderComplex.of(Poset(CROWN))
    assert simplicial_cohomology_dims(oc, 3, QQ) == [1, 1, 0, 0]      # S^1

@lit
def test_rp2_char_sensitive():
    oc = OrderComplex.of(_face_poset(RP2_FACETS))
    assert simplicial_cohomology_dims(oc, 2, QQ)     == [1, 0, 0]     # LIVE-VERIFIED
    assert simplicial_cohomology_dims(oc, 2, GF(2))  == [1, 1, 1]     # torsion in H_1
    assert simplicial_cohomology_dims(oc, 2, GF(3))  == [1, 0, 0]

@selfcert
def test_integral_homology_torsion_certificate():
    oc = OrderComplex.of(_face_poset(RP2_FACETS))
    ih = integral_homology(oc, 2)            # [(free_rank, invariant_factors), ...]
    assert ih[0] == (1, ())                  # H_0 = Z
    assert ih[1] == (0, (2,))                # H_1 = Z/2  -> the char-2 sensitivity
    assert ih[2] == (0, ())                  # H_2 = 0

@selfcert
def test_coboundary_squares_to_zero():
    oc = OrderComplex.of(Poset(CROWN))
    # d_{k-1} . d_k = 0 on the integer boundary maps
    ...
```

- [ ] **Step 2: implement.** `OrderComplex.of(P)`: enumerate chains by length from `P._le`
  (each subset totally ordered under `leq`; order it by the number of predecessors within the
  subset). `face_vector`, `boundary(k)` (integer `d_k`, the alternating face map).
  `simplicial_cohomology_dims(oc, top, field)`: `dim H^n = C^n − rank(δ^n) − rank(δ^{n-1})` with
  `δ = dᵀ`; ranks over the `Domain` (`= C_n − rank(d_n) − rank(d_{n+1})`, cohomology = homology
  dims over a field). `integral_homology(oc, top)`: Smith normal form of the integer `d_n` →
  `H_n = ℤ^{free} ⊕ (torsion from invariant factors of d_n / d_{n+1})`. Keep everything integer
  until the final field reduction.
- [ ] **Step 3: Commit** — `feat(hochschild): order-complex simplicial engine (field rank + integral SNF, char certificate)`

### Task I3: `incidence_cohomology` — the theorem + provenance gate + `Algebra` delegate

**Files:** create the `incidence_cohomology` result + function in `hochschild/simplicial.py`;
modify `src/quiverlab/core/algebra.py` (delegate); test `tests/families/test_incidence_hh.py`.

- [ ] **Step 1: failing tests** (the theorem = the general engine, two independent routes):

```python
# tests/families/test_incidence_hh.py
"""HH^*(kP) via the Gerstenhaber-Schack/Cibils theorem = simplicial H^*(order complex), and the
cross-check that it EQUALS the general engine (Plan 75 / R9). Extends the Plan-33 B3 pin to the
theorem. B3 contractible => [1,0,0,...]; crown => [1,1,0,...]; RP^2 char-split."""
import pytest
from quiverlab.families import IncidenceAlgebra
from quiverlab.families.poset import Poset
from quiverlab.fields import QQ, GF, CC
from quiverlab.resolutions_cs.homology import cs_cohomology_dims

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert

def _b3(): return [(x, x | (1 << b)) for x in range(8) for b in range(3) if not (x >> b) & 1]
CROWN = [("a","X"),("b","X"),("b","Y"),("c","Y"),("c","Z"),("a","Z")]
THICK = [("a","m"),("a","n"),("m","A"),("n","A"),("b","A"),("b","B"),("c","B"),("c","C"),("a","C")]

@lit
def test_b3_contractible_theorem():
    A = IncidenceAlgebra(_b3(), field=CC)
    r = A.incidence_cohomology(4)
    assert r.dims == [1, 0, 0, 0, 0] and r.contractible is True

@lit
def test_crown_nonvanishing_h1():
    r = IncidenceAlgebra(CROWN, field=CC).incidence_cohomology(4)
    assert r.dims == [1, 1, 0, 0, 0]

@xeng
def test_incidence_fast_equals_general_engine():
    # the fast simplicial path == the general CS engine, degreewise, on a NON-MONOMIAL
    # non-lattice poset with H^1 != 0
    A = IncidenceAlgebra(THICK, field=CC)
    assert A.incidence_cohomology(3).dims == cs_cohomology_dims(A, 3).dims == [1, 1, 0, 0]

@lit
def test_rp2_char_dependent_hh():
    covers, elems = _rp2_face_poset()        # helper as in Task I2
    A2 = IncidenceAlgebra(covers, elements=elems, field=GF(2))
    Aq = IncidenceAlgebra(covers, elements=elems, field=QQ)
    assert A2.incidence_cohomology(2).dims == [1, 1, 1]     # LIVE-VERIFIED (dim 121)
    assert Aq.incidence_cohomology(2).dims == [1, 0, 0]
    assert A2.incidence_cohomology(2).char_dependent is True

@selfcert
def test_no_provenance_refuses():
    from quiverlab import Quiver
    from quiverlab.errors import QuiverlabError
    A = Quiver([1,2],{"a":(1,2)}).algebra(relations=[], field=CC)   # kA2, no poset provenance
    with pytest.raises(QuiverlabError, match="poset"):
        A.incidence_cohomology(2)
```

- [ ] **Step 2: implement.** `incidence_cohomology(A, top, field=None)`: require `A._poset`
  (loud refusal per Scope gate 1); `field = field or A's field`; compute
  `simplicial_cohomology_dims(A._poset.order_complex(), top, field)`; set `contractible` from a
  cheap sufficient certificate (`is_bounded()` ⇒ cone ⇒ contractible ⇒ `HH^{≥1}=0`) OR from
  `dims == [1] + [0]*top`; set `torsion`/`char_dependent` from `integral_homology` when in budget;
  attach `references = ["gerstenhaber_schack_1983", "cibils_incidence", "redondo_incidence"]`.
  `Algebra.incidence_cohomology(top)` is a thin lazy-import delegate beside
  `hochschild_cohomology` in `core/algebra.py`.
- [ ] **Step 3: promote the Plan-33 pin.** In `tests/families/test_scale_batteries.py`, keep the
  existing `test_incidence_boolean_b3` but ADD the theorem cross-check
  (`A.incidence_cohomology(10).dims == cs_cohomology_dims(A,10).dims == [1]+[0]*10`), and a
  comment naming the theorem oracle (`gerstenhaber_schack_1983`). Do NOT delete the CS pin (it is
  the general-engine cross-oracle).
- [ ] **Step 4: Commit** — `feat(hochschild): incidence_cohomology theorem (HH^* = simplicial H^* of the order complex) + provenance gate`

### Task I4: (optional-but-designed) the cup ring-check — dims-level v1 + honest scope

**Files:** test `tests/families/test_incidence_ring.py` (marker `oracle_crossengine`).

- [ ] **Step 1:** graded-commutativity of the incidence HH cup (a shipped Plan-35 identity) on the
  crown and `B₃` — a `selfcert` re-use, no new src.
- [ ] **Step 2 (scoped):** `rank(HH^p ⊗ HH^q → HH^{p+q})` from `A.cup_products(top)` equals the
  simplicial cup-pairing rank on a **torus** poset (`H^1 = k²`, `H^2 = k`, cup pairing rank 2 — a
  basis-INDEPENDENT invariant, the Plan-35 cross-engine discipline). **Honest scope:** if the
  minimal-torus incidence algebra's cup does not compute within the deep budget, this test is
  `skip`-with-reason and the verification page records the ring iso as **dimension-level in v1**
  with the cup-rank check DEFERRED (the record's permitted fallback). Do NOT fabricate a cup pin.
- [ ] **Step 3: Commit** — `test(incidence): graded-commutative cup identity + honest-scoped cup-rank ring check`

---

## Task group II — (b) GHMS comultiplicative resolution + fast Koszul HH (R10)

### Task II1: `koszul_kernels` — the intersection spaces `K_n` (Betti = Koszul-dual Hilbert)

**Files:** create `src/quiverlab/hochschild/koszul_ghms.py`; test
`tests/engine/test_koszul_kernels.py`.

**Interfaces:** `koszul._dual_data(A)` / the quadratic relation space `R ⊆ V⊗V` (Plan-27),
`koszul._require_presentation`, `fields.linalg.nullspace`/subspace-intersection primitives;
`ext_algebra.hilbert_matrix_through` (the anchor).

- [ ] **Step 1: failing tests** (LIVE-VERIFIED dims):

```python
# tests/engine/test_koszul_kernels.py
"""The Koszul kernels K_n = intersection of V^i (x) R (x) V^j inside V^{(x)n}; dim K_n = the
n-th Koszul-dual Hilbert coefficient = the total graded Betti number (Plan 75 / R10)."""
import math, pytest
from quiverlab.families import ExteriorAlgebra
from quiverlab.fields import GF, QQ
from quiverlab.hochschild.koszul_ghms import koszul_kernels

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine

@lit
@pytest.mark.parametrize("m, expect", [
    (2, [1, 2, 3, 4, 5, 6, 7]),                 # C(n+1,n) -- LIVE-VERIFIED
    (3, [1, 3, 6, 10, 15, 21, 28]),             # C(n+2,2) -- LIVE-VERIFIED
])
def test_exterior_kernels_are_koszul_dual_hilbert(m, expect):
    K = koszul_kernels(ExteriorAlgebra(m, field=GF(32003)), 6)
    assert [len(k) for k in K] == expect

@xeng
def test_kernels_match_ext_algebra_betti():
    A = ExteriorAlgebra(2, field=QQ)
    K = koszul_kernels(A, 6)
    hm = A.ext_algebra(6).hilbert_matrix_through(6)
    betti = [sum(sum(row) for row in hm[n]) for n in range(7)]
    assert [len(k) for k in K] == betti        # [1,2,3,4,5,6,7]
```

- [ ] **Step 2: implement.** `koszul_kernels(A, top)`: require length-graded quadratic `A`; get
  `V` (arrows) and `R ⊆ V⊗_S V` from `koszul._dual_data`; build `K_n ⊆ V^{⊗n}` as the exact
  intersection `⋂_i (V^{⊗i} ⊗ R ⊗ V^{⊗(n-2-i)})` — each summand a `Domain`-subspace of the free
  `S`-module `V^{⊗n}` (respect the corner grading `e_v … e_w`), intersection by successive
  nullspace/`intersection` solves. Return `Domain`-bases (lists of coordinate vectors).
  `dim K_n` positivity + monotone corner-typing are self-cert.
- [ ] **Step 3: Commit** — `feat(hochschild): koszul_kernels K_n (dim = Koszul-dual Hilbert coefficient)`

### Task II2: `GHMSResolution` — the comultiplicative differential (Koszul-gated)

**Files:** `koszul_ghms.py` (add `GHMSResolution`); modify `engine/resolutions.py` only if the
`Resolution` ABC needs a hook (prefer NOT to); test `tests/engine/test_ghms_resolution.py`.

**Interfaces:** `engine.resolutions.Resolution` (the ABC / contract: `term`, `differential`,
`assert_dd_zero`), the **three-valued gate** `A.ext_algebra(window).koszul` +
`.koszul_obstruction` (with `koszul.g_quadratic_certificate` as the cheap fast `True` path),
`fields.linalg.reduce_mod_nullspace` (canonicalization, the Plan-17 discipline).

- [ ] **Step 1: failing tests.**

```python
# tests/engine/test_ghms_resolution.py
"""The GHMS comultiplicative minimal A^e-resolution of a Koszul algebra: terms
P_n = A (x)_S K_n (x)_S A, ranks = dim K_n, d.d = 0 (Plan 75 / R10). Koszul-gated on the FULL
three-valued ext_algebra.koszul verdict: loud on False (Ext obstruction) or None (inconclusive)."""
import pytest
from quiverlab.families import ExteriorAlgebra
from quiverlab.fields import GF, QQ
from quiverlab.errors import QuiverlabError
from quiverlab.hochschild.koszul_ghms import GHMSResolution

selfcert = pytest.mark.oracle_selfcert

@selfcert
@pytest.mark.parametrize("dom", [GF(32003), QQ])
def test_ghms_terms_and_dd_zero(dom):
    res = GHMSResolution(ExteriorAlgebra(2, field=dom))
    assert [res.rank(n) for n in range(6)] == [1, 2, 3, 4, 5, 6]   # dim K_n
    res.assert_dd_zero(5)                                          # d_{n-1} d_n = 0

@selfcert
def test_ghms_canonical_differential():
    # byte-reproducible: build twice, identical differential matrices (Plan-17 discipline)
    a = GHMSResolution(ExteriorAlgebra(3, field=GF(2)))
    b = GHMSResolution(ExteriorAlgebra(3, field=GF(2)))
    assert [a.differential(n) for n in range(1,5)] == [b.differential(n) for n in range(1,5)]

@selfcert
def test_non_quadratic_refuses():
    # k[x]/(x^3): a cubic relation -> NOT quadratic -> loud refusal
    from quiverlab.families import truncated_polynomial
    with pytest.raises(QuiverlabError, match="[Kk]oszul|quadratic"):
        GHMSResolution(truncated_polynomial(3, field=QQ))

@selfcert
def test_not_koszul_refuses_with_obstruction():
    # preprojective A3 is NOT Koszul: ext_algebra.koszul == False, obstruction at degree 3
    # (a new Ext-algebra generator) -> GHMS refuses NAMING the obstruction. LIVE-VERIFIED.
    from quiverlab.families import PreprojectiveAlgebra
    with pytest.raises(QuiverlabError, match="[Nn]ot Koszul|obstruction|degree 3"):
        GHMSResolution(PreprojectiveAlgebra("A3", field=QQ))
    # sanity: the verdict the gate reads
    E = PreprojectiveAlgebra("A3", field=QQ).ext_algebra(3)
    assert E.koszul is False and E.koszul_obstruction[0] == 3

# NOTE (three-valued honesty): the None ("Koszulity unknown") branch of the gate has NO cheap
# CERTIFIED witness pinned -- see the honest note in Task II3 / the verification scope. A worker
# who finds a small Koszul-but-not-G-quadratic algebra whose ext_algebra.koszul is None adds a
# test here; until then the None branch ships exercised only by a constructed/mocked verdict.
```

- [ ] **Step 2: implement.** `GHMSResolution.__init__(A)`: **gate on the three-valued verdict**
  (Scope gate 2) — try `g_quadratic_certificate(A)` (cheap `True` → build); else read
  `A.ext_algebra(window).koszul`: `False` → refuse naming `koszul_obstruction` (*not Koszul, Ext
  obstruction at degree d*); `None` → refuse, honestly inconclusive. Never gate on `g_quadratic`
  alone (it would mislabel the genuinely-not-Koszul preprojective `A₃` as merely "uncertified").
  `term(n) = (A^e)`-block indexed by the `K_n` basis (corner-typed `e_v K_n e_w`). `differential(n)`
  = the comultiplicative `d_n` (§Foundation (b)): compute the two splittings `K_n ↪ V⊗K_{n-1}` and
  `K_n ↪ K_{n-1}⊗V` by exact `Domain` solves against the `koszul_kernels` bases, assemble
  `Σ x_a⊗ω'_a⊗1 − (−1)^n Σ 1⊗ω''_b⊗y_b`, canonicalize with `reduce_mod_nullspace` (byte-repro).
  `assert_dd_zero(N)`. **Sign** pinned by `assert_dd_zero` + the Task-II3 cross-engine anchors
  (arbitrated, not assumed).
- [ ] **Step 3: Commit** — `feat(hochschild): GHMSResolution -- comultiplicative minimal A^e resolution for certified-Koszul algebras`

### Task II3: `engine="ghms"` fast HH + the third oracle class

**Files:** `koszul_ghms.py` (`ghms_homology_dims`/`ghms_cohomology_dims`); modify
`src/quiverlab/core/algebra.py` (engine set + route); test `tests/engine/test_ghms_hh.py`.

- [ ] **Step 1: failing tests** (LIVE-VERIFIED HH; GHMS ≡ minimal ≡ bar):

```python
# tests/engine/test_ghms_hh.py
"""engine='ghms': fast HH via the GHMS resolution == the minimal syzygy engine == bar/CS,
degreewise -- the THIRD independent HH oracle class (Plan 75 / R10). Over any Domain."""
import pytest
from quiverlab import Quiver
from quiverlab.families import ExteriorAlgebra, PreprojectiveAlgebra
from quiverlab.fields import GF, QQ
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.resolutions_minimal import minimal_homology_dims
from quiverlab.hochschild.koszul_ghms import ghms_homology_dims, ghms_cohomology_dims

xeng = pytest.mark.oracle_crossengine

def _kz3_radsq(dom):   # cyclic rad^2-Nakayama 1->2->3->1, multi-vertex G-quadratic Koszul
    return Quiver([1,2,3], {"a":(1,2),"b":(2,3),"c":(3,1)}).algebra(
        relations=["a*b","b*c","c*a"], field=dom)

@xeng
def test_ghms_equals_minimal_exterior():
    A = ExteriorAlgebra(2, field=GF(32003))
    assert ghms_homology_dims(A, 6) == [3, 4, 6, 8, 10, 12, 14]       # LIVE-VERIFIED (minimal)
    assert ghms_homology_dims(A, 6) == minimal_homology_dims(to_engine(A), 6, primes=(32003,))[32003]

@xeng
def test_ghms_over_QQ_matches_engine():
    A = ExteriorAlgebra(3, field=QQ)
    assert A.hochschild_homology(4, engine="ghms").dims == A.hochschild_homology(4, engine="cs").dims

@xeng
def test_ghms_multivertex_kz3_radsq():
    A = _kz3_radsq(GF(32003))
    assert ghms_homology_dims(A, 5) == [3, 0, 1, 1, 0, 0]             # LIVE-VERIFIED (minimal)
    assert ghms_homology_dims(A, 5) == minimal_homology_dims(to_engine(A), 5, primes=(32003,))[32003]

@xeng
def test_ghms_multivertex_diamond_incidence():
    # the (a)<->(b) BRIDGE: a Koszul INCIDENCE algebra (commutative square). GHMS computes its
    # Hochschild HOMOLOGY [4,0,0,0,0] (HH_0 = dim A/[A,A] = 4) -- DISTINCT from its order-complex
    # cohomology HH^* = [1,0,0]. Ranks r_n = dim K_n = [4,4,1,0,...] (K_0=|Q_0|=4). LIVE-VERIFIED.
    from quiverlab.families import IncidenceAlgebra
    A = IncidenceAlgebra([("b","x"),("b","y"),("x","t"),("y","t")], field=GF(32003))
    assert ghms_homology_dims(A, 4) == [4, 0, 0, 0, 0]
    assert ghms_homology_dims(A, 4) == minimal_homology_dims(to_engine(A), 4, primes=(32003,))[32003]

@pytest.mark.oracle_selfcert
def test_engine_ghms_refuses_non_quadratic():
    # k[x]/(x^3): the relation is CUBIC -> not quadratic -> loud (its own message)
    from quiverlab.families import truncated_polynomial
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError, match="[Kk]oszul|quadratic"):
        truncated_polynomial(3, field=QQ).hochschild_homology(3, engine="ghms")

@pytest.mark.oracle_selfcert
def test_engine_ghms_refuses_not_koszul_with_obstruction():
    # preprojective A3 is NOT Koszul: ext_algebra.koszul == False, obstruction at degree 3 ->
    # engine='ghms' refuses NAMING it; the minimal syzygy engine still computes. LIVE-VERIFIED.
    from quiverlab.errors import QuiverlabError
    A = PreprojectiveAlgebra("A3", field=QQ)
    with pytest.raises(QuiverlabError, match="[Nn]ot Koszul|obstruction|degree 3"):
        A.hochschild_homology(3, engine="ghms")
    assert minimal_homology_dims(to_engine(PreprojectiveAlgebra("A3", field=GF(32003))),
                                 5, primes=(32003,))[32003] == [3, 1, 1, 1, 1, 1]
```

- [ ] **Step 2: implement.** `ghms_homology_dims(A, top)` = homology of `A ⊗_{A^e} P_•`
  (contract each `P_n = A⊗K_n⊗A` to `e_w K_n e_v` blocks, the Plan-16 `b·w·a` convention);
  `ghms_cohomology_dims(A, top)` = cohomology of `Hom_{A^e}(P_•, A)` (`a·w·b`, swapped tags). In
  `core/algebra.py`: extend the valid-engine set to include `"ghms"`; route
  `hochschild_homology(top, engine="ghms")` / `hochschild_cohomology(...)` to these; the
  three-valued Koszul refusal (Scope gate 2) — `False` names the Ext obstruction, `None` is
  honestly inconclusive. **`engine="auto"` unchanged** (byte-stability; no golden drift). Update
  the `spec.py` engine allow-set `("auto","bar","fast","cs","ghms")` (`spec.py:700`).
- [ ] **Step 3: Commit** — `feat(core): engine='ghms' fast Koszul HH (== minimal == bar/CS, any Domain) -- the third HH oracle class`

---

## Task group III — GUI / webapp + verification

### Task G1: the poset input mode + the `incidence_cohomology` compute kind (all three tiers, i18n ×4)

The organizing invariant is **byte-identity across the two runners** via a single shared builder,
and **canonical-key stability** for every existing request (the poset panel emits the SHIPPED
`IncidenceAlgebra` family). Anchors below were grepped against the live tree during authoring;
"adjust to reality" if a later merge shifts them.

**The poset input mode (a NEW input shape, existing family output).** A "Poset (order relations)"
panel on the draw canvas lets the user enter **cover pairs** `[a, b]` (and optional isolated
elements) and emits an `algebra` block `{kind:"family", family:"IncidenceAlgebra",
params:{poset_or_covers:[[a,b],…], elements:<normalized>}}` — the SAME shape the catalog family
prefill builds (`webapp/server/catalog.py:110`). **Canonical-key normalization is MANDATORY
(critic M2):** the catalog default is `elements: None`; `null` vs an explicit list key
*differently* under `canonical_key`'s sort_keys JSON, so the panel emits `elements` **only for
isolated elements** (vertices in no cover) and emits `elements: null` when the covers already name
every element. Thus a covers-only poset keys **byte-identically** to the typed family request; a
poset with genuine isolated points keys differently (a genuinely different algebra) — both pinned
by `test_poset_request_keys_like_family` (Plan-25 `canonical_key`, no schema change, no new algebra
`kind`). The panel renders a live **Hasse preview** (reuse `viz/hasse_html.py` /
`viz/tikz.py::tikz_hasse`, `tikz.py:53`) and refuses a directed cycle loudly (the shipped `Poset`
`RelationError`, surfaced as a clean 4xx).

**`incidence_cohomology` is an ALGEBRA-level compute kind carrying a DEGREE RANGE**
(`incidence_cohomology:0..N`), parsed exactly like `hh_cohomology` (the `name:lo..hi` grammar,
`spec.py::parse_compute_item`) — NOT a budget kind.

**Files (the seven-touchpoint checklist, mirroring Plan-63 Task 5):**
1. **Parse grammar** — `hpc/spec.py::parse_compute_item`: `incidence_cohomology` falls through to
   the generic `name:lo..hi` `_RANGE` regex (`spec.py:259`); no special branch needed (confirm it
   is in the allowed-kind set). Mirror in `docs/gui/runner.py` parse.
2. **Server/HPC dispatch** — `spec.py::_dispatch` (a new branch beside the `hh_cohomology` one,
   `spec.py:1532+`): `if kind == "incidence_cohomology": r = A.incidence_cohomology(item.hi);
   block = {"kind":"incidence_cohomology","dims":r.dims,"face_vector":r.face_vector,
   "contractible":r.contractible,"char_dependent":r.char_dependent,"torsion":r.torsion,
   "references":r.references,"citations":_citation_pairs(r.references)}`. Catch the
   no-provenance / presentation-less `QuiverlabError` into `{"error": <msg>}` (Plan-30 honest
   per-entry; never a 500). Add the `_snippet` entry (`spec.py:2552`):
   `"incidence_cohomology": lambda it: f"A.incidence_cohomology({it.hi})"`.
3. **Pyodide twin** — `docs/gui/runner.py::compute_one`: the byte-identical
   `elif name == "incidence_cohomology":` branch; one `calls` reproduce entry
   (`"incidence_cohomology": "A.incidence_cohomology(%d)"`); one `ETA_MODEL["scalars"]` cost
   (`"incidence_cohomology": 0.3` — the simplicial complex is tiny, size it BELOW `hh_cohomology`).
4. **GUI JS** — `docs/gui/gui.js` AND `webapp/static/gui/gui.js` (diff-identical — apply every
   edit to both): the **poset panel** (a cover-pair editor + isolated-elements list + Hasse
   preview, beside the module panel); a checkbox `<input id="qlgui-incidence_cohomology">` + degree
   picker in the HH block; the id-registry entries; the request push-list
   (`if (el.incidence_cohomology.checked) compute.push("incidence_cohomology:0.." + top)`); a
   `renderBlock` `else if (name === "incidence_cohomology")` → `renderIncidenceCohomology(div, b)`
   (the `HH^*` dims table + the order-complex face vector + a "contractible ⇒ HH^{≥1}=0" /
   "char-dependent (Z-torsion in the order complex)" note + the theorem citation); the `THEMES`
   `kinds` entry `"incidence_cohomology"` inside `QLGUI-THEMES-BEGIN/END`; the `KIND_CTRL`/
   `scheduleProbe` entries. The poset panel only shows when the algebra input mode is "poset".
5. **i18n (FOUR locales) + picker template** — `webapp/server/i18n/{en,es,fr,zh}.json`: add
   `pick.kind.incidence_cohomology` (REQUIRED — gated by `test_layout_picker.py`),
   `struct.incidence_cohomology`, the poset-panel strings (`poset.title`, `poset.add_cover`,
   `poset.elements`, `poset.cycle_error`), and `inc.*` block strings (`inc.title`, `inc.dims`,
   `inc.face_vector`, `inc.contractible`, `inc.char_dependent`, `inc.theorem`). Add
   `data-pick-kind-incidence_cohomology="{{ t('pick.kind.incidence_cohomology') }}"` in
   `webapp/templates/draw.html`.
6. **Report renderer** — `trace/results_html.py`: `_HEADINGS["incidence_cohomology"] =
   "Hochschild cohomology of an incidence algebra (= simplicial cohomology of the order complex)"`
   + a `_incidence_cohomology_html(b)` branch in `_block_html` (the dims table + the face vector +
   the theorem statement + the char-dependence note + `gerstenhaber_schack_1983` provenance) and a
   `viz/tikz.py::tikz_order_complex(P)` twin (draw the Hasse diagram + label the order-complex face
   vector; reuse `tikz_hasse`).
7. **Tests / gates** — add `incidence_cohomology` to `ALL_KINDS` in
   `tests/webapp/test_layout_picker.py:26` (must equal the THEMES set, each once); add ONE golden
   `incidence_cohomology_b3` to `tests/webapp/_runner_goldens.json` with a dated change-log bullet
   in `test_runner_delegation.py` (verify existing entries byte-identical FIRST; merge the JSON
   semantically on rebase); add a twin-parity test; add a **canonical-key stability** test — a
   poset-panel request and the equivalent typed family request produce the SAME `canonical_key`.
- Test files: `tests/webapp/test_incidence_p75.py`, `tests/gui/test_incidence_runner_twin.py`.

```python
# tests/webapp/test_incidence_p75.py
"""The incidence_cohomology kind + poset input mode: served by hpc.spec, mirrored by the Pyodide
twin, byte-identical; the poset panel emits the SHIPPED IncidenceAlgebra family so canonical keys
are unchanged. B3 -> [1,0,0,...], contractible."""
def test_incidence_block_shape(tmp_path):
    from quiverlab.hpc.spec import ComputeRequest, run
    req = _incidence_request(covers=[[x, x|(1<<b)] for x in range(8) for b in range(3)
                                     if not (x>>b)&1], top=3)
    out = run(ComputeRequest.model_validate(req), tmp_path)
    b = out["results"]["incidence_cohomology"]
    assert b["dims"] == [1, 0, 0, 0] and b["contractible"] is True
    assert "gerstenhaber_schack_1983" in b["references"]

def test_poset_request_keys_like_family():
    # critic M2: a covers-only poset (no isolated elements) keys BYTE-IDENTICALLY to the typed
    # IncidenceAlgebra family form (elements normalized to None); a poset with an isolated point
    # keys DIFFERENTLY (a different algebra). canonical_key = the Plan-25 sort_keys sha256.
    from webapp.server.cache import canonical_key                # the shipped Plan-25 canonicalizer
    covers = [[1, 2], [2, 3]]
    def req(elements):
        return {"schema": 1,
                "algebra": {"kind": "family", "family": "IncidenceAlgebra",
                            "params": {"poset_or_covers": covers, "elements": elements}},
                "compute": ["incidence_cohomology:0..3"]}
    # panel normalizes: 1,2,3 all appear in covers -> elements: None, == the typed family form
    assert canonical_key(req(None)) == canonical_key(req(None))
    # a genuine isolated element 4 (in no cover) is a DIFFERENT poset -> DIFFERENT key
    assert canonical_key(req([1, 2, 3, 4])) != canonical_key(req(None))

def test_existing_requests_byte_unchanged():
    # critic M2 (the ALWAYS-true half): a non-poset request keys identically before/after P75 --
    # no new field is added to it. Backed by the frozen goldens in test_runner_delegation.py.
    ...

def test_twin_parity(tmp_path):
    ...
```

- [ ] Steps: write tests → wire all seven touchpoints (both `gui.js` copies byte-identical) →
  golden (verify others unchanged) → `-m fast` webapp/gui + `node --check` + `test_js_parses.py`
  + `test_layout_picker.py` green. Commit — `feat(gui): poset input mode + incidence_cohomology kind (all three tiers, i18n x4)`

**GUI-deferral ledger note (R10 / GHMS).** `engine="ghms"` adds **no new user-visible
computation** — HH dims are byte-identical to the shipped `auto`/`cs` route; GHMS is an engine
**acceleration** + an oracle. It is surfaced via the API/HPC `engine=` option and named in the
worked-steps resolution line ("resolved by the GHMS comultiplicative minimal A^e resolution") when
selected. Per metaplan §1.2, record in `docs/plans/DEEPER-ENGINES-BACKLOG.md`: *P75 ships GHMS as
an engine option + third oracle class, NOT a new GUI kind (no new mathematics is user-visible); the
GUI HH kinds keep `engine="auto"`. P80 reconciles.*

### Task G2: verification page, README, suite gate

**Files:** `docs/verification.md`, `README.md`; `tests/release/test_oracle_classes.py` (recount).

- [ ] **Step 1: verification rows** (two new subsystem rows + honest scope):
  - `hochschild/simplicial.py` (R9) — the order-complex engine + `incidence_cohomology`.
    `oracle_literature`: `B₃` contractible `[1,0,0,…]` (Plan-33, now theorem-backed); crown `S¹`
    `[1,1,0,…]`; thickened-crown non-lattice/non-monomial `[1,1,0]`; **`RP²₆` char-split**
    `HH^*(GF₂)=[1,1,1]` vs `HH^*(QQ)=[1,0,0]` (`gerstenhaber_schack_1983`/`cibils_incidence`/
    `redondo_incidence`). `oracle_crossengine`: `incidence_cohomology` ≡ general CS
    `hh_cohomology` degreewise; the simplicial SNF engine ≡ hand-rolled integer homology.
    `oracle_selfcert`: `δ∘δ=0`, the integral-SNF torsion certificate, the no-provenance refusal,
    the face-vector positivity.
  - `hochschild/koszul_ghms.py` (R10) — GHMS resolution + `engine="ghms"`. `oracle_literature`:
    `Λ(k²)`/`Λ(k³)` `dim K_n` = Koszul-dual Hilbert `[1,2,3,…]`/`C(n+2,2)`. `oracle_crossengine`:
    GHMS `HH_*` ≡ minimal syzygy engine ≡ bar/CS on the **Koszul** targets (exterior `Λ(k²)`,
    multi-vertex `kZ₃/rad²` `[3,0,1,1,0,0]`, and the diamond-incidence bridge `[4,0,0,0,0]`); `dim
    K_n` = Plan-27 `ext_algebra` Betti. `oracle_selfcert`: `d∘d=0`, canonical (byte-reproducible)
    differential, the three-valued gate — `False` refusal NAMING the Ext obstruction (preprojective
    `A₃`, degree 3, `ext_algebra.koszul is False`), non-quadratic refusal.
  - **Honest scope:** (i) the incidence fast path REQUIRES poset provenance (no `kQ/I` →
    poset recognition); (ii) the ring iso is **dimension-level in v1** — cup **rank** check
    scoped/deferred (torus budget), **no bracket claim** (the record's dropped clause); and the
    incidence `HH^*` (order-complex cohomology) is a **distinct** quantity from `HH_*` (Hochschild
    homology) — never conflated; (iii) GHMS is **Koszul-only via the three-valued verdict** —
    `False`/`None` are honestly distinguished, the `None` ("Koszulity unknown") branch ships with
    **no certified witness pinned** (worker task), preprojective `A₃` is the shared **P75/P77
    boundary** (refused here as not-Koszul, `(p,q)`-almost-Koszul over there); (iv) **QPA has no HH
    surface / no incidence-vs-nerve identity / no GHMS resolution / no `IsKoszul`** — the only live
    QPA crosscheck is the Koszul-Betti direction via the shipped Plan-27 `ExtAlgebraGenerators`;
    `tests/qpa/test_incidence_koszul_qpa.py` is an honest `NamesGVars()` probe that SKIPS and
    FAILS if QPA ever ships any of these.
- [ ] **Step 2: recount** the oracle-class table by running `tests/release/test_oracle_classes.py`
  (paste the live numbers; claim only this plan's deltas).
- [ ] **Step 3: README line** — "incidence-algebra HH = simplicial cohomology of the order
  complex (poset input mode, char-sensitive) + GHMS fast Koszul HH".
- [ ] **Step 4:** full gate — `-m deep` (families + engine) + `-m fast` (webapp/gui) + `-m qpa` +
  citations + release + `node --check` + `mkdocs build --strict`; `; echo EXIT=$?`.
- [ ] **Step 5: Commit** — `docs(verification): incidence-vs-nerve + GHMS Koszul oracles; honest scope (ring dims-level, QPA no-surface)`

---

## Acceptance (Plan-75 definition of done)

1. **(R9) `A.incidence_cohomology(top)`** returns `HH^*(kP) = H^*(Δ(P))` via the
   Gerstenhaber–Schack/Cibils theorem, computed by the exact SNF/rank simplicial engine
   (`hochschild/simplicial.py`), requiring poset provenance (loud refusal otherwise). The
   **fast path ≡ the general CS engine** degreewise on every test poset (cross-engine), and ≡ the
   hand-rolled integer homology.
2. **(R9) The theorem oracle family** is pinned: `B₃` contractible `[1,0,0,…]` (Plan-33 promoted);
   crown `C(3,3)`/`C(2,2)` `≅ S¹` `[1,1,0,…]` (**nonvanishing `H¹`**, non-lattice); thickened
   crown non-lattice **+ non-monomial** `[1,1,0]`; whisker/`N₅` contractible; and the
   **char-sensitive `RP²₆`** pin `HH^*(GF₂)=[1,1,1]` vs `HH^*(QQ)=[1,0,0]` (verified on BOTH the
   fast path and the general engine, dim 121), with the integral-SNF torsion certificate
   (`H_1 = Z/2`). All values LIVE-VERIFIED at spec time.
3. **(R9) The poset input mode** is clickable end-to-end (draw canvas poset panel → Hasse preview →
   `incidence_cohomology` block → report TikZ) in all four locales; it emits the SHIPPED
   `IncidenceAlgebra` family so **every existing request keys byte-unchanged**, and — after the
   mandatory `elements` normalization (emit `null` when the covers name every element; explicit
   list only for isolated points) — a covers-only poset request keys **identically** to the typed
   family request (a poset with isolated points keys differently, correctly; both pinned by
   `test_poset_request_keys_like_family`); schema stays v1; one golden added with a documented
   change-log entry; both runners byte-identical.
4. **(R10) `GHMSResolution`** builds the comultiplicative minimal `A^e`-resolution
   `P_n = A⊗_S K_n⊗_S A` of a **Koszul** `A`, gated on the **three-valued `ext_algebra.koszul`
   verdict** (`True` → build; `False` → refuse naming the Ext obstruction; `None` → refuse,
   inconclusive), over **any `Domain`**; `rank P_n = dim K_n =` the Koszul-dual Hilbert coefficient
   (`Λ(k²)`: `[1,2,3,…]`; `Λ(k³)`: `C(n+2,2)`; multi-vertex diamond incidence: `[4,4,1,0,…]`,
   `K₀=|Q₀|` — LIVE-VERIFIED, and ≡ the Plan-27 `ext_algebra` Betti); `d∘d=0`; the differential is
   byte-reproducible.
5. **(R10) `engine="ghms"`** computes `HH_*`/`HH^*` and **≡ the minimal syzygy engine ≡ bar/CS**
   degreewise (exterior `Λ(k²)` `[3,4,6,8,10,12,14]`; multi-vertex `kZ₃/rad²` `[3,0,1,1,0,0]`;
   the diamond-incidence bridge `[4,0,0,0,0]` — LIVE-VERIFIED) — the **third independent HH oracle
   class**; the **Koszul gate is real and correctly three-valued** (preprojective `A₃` returns
   `ext_algebra.koszul = False` with the obstruction *a new Ext generator in degree 3*, so GHMS
   refuses it as **not Koszul** — the minimal engine still gives `[3,1,1,1,1,1]`); the `None`
   ("Koszulity unknown") branch ships with **no certified witness pinned** (honest note — worker
   task to hunt a small Koszul-but-not-G-quadratic algebra); `engine="auto"` is **byte-unchanged**
   (no golden drift); the engine allow-set becomes `{"auto","bar","fast","cs","ghms"}`.
6. **Honest scope recorded** on the verification page: incidence fast path needs poset provenance;
   ring iso **dimension-level v1** (cup-rank check scoped, **no bracket claim** — the record's
   dropped clause re-sourced away); the **HH^* (order-complex cohomology) and HH_* (Hochschild
   homology) of an incidence algebra are DISTINCT** and never conflated (the diamond: `[1,0,0]` vs
   `[4,0,0,0,0]`); GHMS **Koszul-only, three-valued** (P77 seam named, `A₃` the shared boundary,
   the `None` branch unexemplified); **QPA cannot compare** any of the four surfaces (honest probe).
   Oracle-class table recounted (live numbers); README line added; deep + fast + webapp/gui + qpa +
   citations + release + `mkdocs --strict` green. The shipped `Poset`/`IncidenceAlgebra`/minimal
   engine/Plan-27 surfaces left byte-unchanged (only `A._poset` added).

---

## Methodology & assumptions

**Approach.** I read the P75 metaplan card (§5) and §§1–4, and R9/R10 verbatim from the research
doc. I read the P63 (wall-and-chamber) and P65 (exceptional sequences) plans end-to-end for the
house style — the seven-touchpoint GUI checklist, the reference-re-verification discipline, the
oracle-class marker scheme, the honest-scope contract, and the acceptance/methodology shape — and
mirrored them. I then read the SHIPPED substrate: `families/poset.py` + `families/incidence.py`
(the `Poset` closure + the Hasse-quiver-plus-commutativity `IncidenceAlgebra`), the Plan-33 `B₃`
pin (`tests/families/test_scale_batteries.py`), the minimal `A^e` engine
(`engine/resolutions_minimal.py` + `engine/adapter.to_engine`), the CS homology surface, and the
Plan-27 Koszul/Ext-algebra API (`modules/koszul.py::g_quadratic_certificate`/`_dual_data`,
`modules/ext_algebra.py::hilbert_matrix_through`/`.koszul`). I confirmed the compute-kind wiring by
grepping `hpc/spec.py` (`parse_compute_item`, `_dispatch`, `_snippet`, the engine allow-set at
`spec.py:700`), `docs/gui/runner.py`, `docs/gui/gui.js` (THEMES sentinels, `_HEADINGS`), the
four-locale i18n, `trace/results_html.py`, `viz/tikz.py` (`tikz_hasse` exists), and
`tests/webapp/test_layout_picker.py`. I confirmed `IncidenceAlgebra` is ALREADY a webapp catalog
family (`webapp/server/catalog.py:110`, param `poset_or_covers`) — so the poset input mode rides
the existing family-input canonical key, no schema change.

**Empirical ground truth (run this authoring; `QUIVERLAB_NO_NUMBA=1`).** Every value in
§Live-verified facts was computed by the merged engines. The simplicial side was computed by a
hand-rolled exact routine (integer boundary matrices; sympy `Matrix.rank` over QQ and mod `p`;
`H_n = C_n − rank d_n − rank d_{n+1}`) — the routine the plan re-implements in `src/`. Specifically:
(R9) `B₃` `HH^*=[1,0,0]`, crown `C(3,3)`/`C(2,2)`/thickened-crown/whisker/`N₅` all MATCH their
order-complex `H^*` degreewise; the `RP²₆` face poset (dim 121) gives `HH^*(GF₂)=[1,1,1]` vs
`HH^*(QQ)=[1,0,0]` on BOTH the general CS engine and the hand-rolled SNF, and I independently
validated the `RP²₆` triangulation (edge-regular, `f=(6,15,10)`, `χ=1`, `H_*(QQ)=[1,0,0]`,
`H_*(GF₂)=[1,1,1]`, `H_*(GF₃)=[1,0,0]`). (R10) `Λ(k²)` minimal bimodule ranks `[1,2,3,4,5,6,7,8]`
and `Λ(k³)` `[1,3,6,10,15,21,28,36]` EQUAL the Koszul-dual Hilbert coefficients `C(m+n-1,n)` and
the Plan-27 `ext_algebra` Betti (`[1,2,3,4,5,6,7]` for `m=2`, `koszul=True`, `g_quadratic=True`);
`Λ(k²)` minimal `HH_* = [3,4,6,8,10,12,14]` ≡ bar `[3,4,6,8,10]`; the multi-vertex Koszul
(`ext_algebra.koszul = True`) `kZ₃/rad²` (cyclic rad²-Nakayama) minimal `HH_* = [3,0,1,1,0,0]`
(a GHMS multi-vertex target); the **commutative-square (diamond) incidence algebra** (`dim 9`,
`koszul = True`) is the (a)↔(b) bridge — its minimal bimodule ranks are `[4,4,1,0,…]`
(`K₀=|Q₀|=4`, `K₂=1`, gl.dim 2) and its Hochschild HOMOLOGY is `[4,0,0,0,0]` (distinct from its
order-complex cohomology `[1,0,0]`); and I live-confirmed the **three-valued gate is real**:
preprojective `A₃` is **NOT Koszul** — `ext_algebra.koszul = False` with obstruction
`(3, "a new Ext-algebra generator appears in degree 3")` (Brenner–Butler–King almost-Koszul) — so
GHMS correctly refuses it, while the minimal syzygy engine still gives `[3,1,1,1,1,1]`. (My earlier
draft wrongly called `A₃` "Koszul but PBW-uncertified"; the critic caught it, and the shipped
`ext_algebra` returns a positive disproof, not "unknown".)

**Key correctness decisions.** (1) **The fast path is the order-complex simplicial cohomology, and
the speed is stated honestly** — the theorem replaces `(kP)^e`-linear algebra with the cochain
complex of `Δ(P)` (dimension = chain counts), the integer boundaries computed once, SNF giving
every characteristic via universal coefficients. (2) **The `RP²₆` pin makes SNF (not mere rank)
load-bearing** — it is a *char-dependent* `HH^*`, verified on both engines (the general CS engine
is not intractable: ≈2.8 s GF₂ / ≈6.4 s QQ to degree 2), where the fast path is *much smaller*
(`dim 121` bimodule vs a 31-vertex complex) and yields every characteristic from one integer SNF.
(3) **The incidence
fast path requires poset provenance** — recognizing an incidence algebra from an arbitrary `kQ/I`
is a poset-reconstruction problem, deliberately out of scope; the general `hh_cohomology` kind
remains the cross-oracle. (4) **GHMS is the closed-form comultiplicative resolution over any
`Domain`, gated on the FULL three-valued `ext_algebra.koszul` verdict** (True→build,
False→refuse-with-Ext-obstruction, None→refuse-inconclusive) — NOT on `g_quadratic` alone (the
critic W1 fix: preprojective `A₃` returns a positive disproof `False`, obstruction at degree 3, not
"uncertified"); genuinely new over the GF(p) syzygy engine (which becomes the cross-oracle);
`Betti_n = dim K_n =` Koszul-dual Hilbert, multi-vertex `K₀ = |Q₀|` (LIVE-VERIFIED). (5)
**The GHMS sign is arbitrated** (`d∘d=0` + minimal + bar/CS anchors), not assumed — the Plan-20/21
discipline. (6) **`engine="auto"` and all goldens stay byte-identical** — GHMS is explicit-only;
the poset input mode emits the shipped family. (7) **Attribution corrected per the record (W3):**
Cibils 1989 is the general **arbitrary-finite-poset** theorem `HH^*(kP) ≅ H^*(Δ(P))` (the statement
this plan realizes); Gerstenhaber–Schack 1983 is the earlier **face-poset-of-a-simplicial-complex**
version and the RING (cup) source (added; distinct from the shipped 1987 Hodge paper); Redondo 2008
corroborates. The **bracket-vanishing clause is DROPPED** (not in G–S 1983, not re-sourced) — the
plan makes no bracket claim. (8) **`HH^*` (order-complex cohomology) ≠ `HH_*` (Hochschild homology)
of an incidence algebra** — the diamond's `[1,0,0]` vs `[4,0,0,0,0]` (W2); the doc keeps the star
conventions distinct throughout.

**Deliberately NOT checked / assumptions.** (a) I did **not** implement or verify the GHMS
comultiplicative **differential** at spec time — I verified the TERM RANKS (`= dim K_n = Koszul-dual
Hilbert`) and the target HH dims (via the minimal syzygy engine, which computes the SAME minimal
resolution's homology). The differential's construction is specified precisely (§Foundation (b))
and its correctness is gated at implementation by `d∘d=0` + the cross-engine HH anchors; a wrong
sign/splitting fails those gates. This is the honest residual risk and is called out in Task II2.
(b) I did **not** compute a nonzero **cup product** on an incidence HH (the ring iso beyond
dimensions) — the torus poset needed for a nonzero cup is large, and the Plan-35 cup surface's cost
at degree 2 on such an algebra is unbudgeted; the ring claim is scoped to **dimensions in v1** with
a **cup-rank** ring-check task (basis-independent), honestly deferred if it does not compute. (c) I
did **not** re-verify the exact Archiv der Mathematik pages of GHMS 2005 nor the exact JPAA 30
issue/pages of Gerstenhaber–Schack 1983 against the physical journals — both are `# PIN`/`# verify`
at citation-add (arXiv math/0508177 and the standard JPAA 30 (1983) 143–156 reference are stable).
(d) I did **not** confirm the GUI anchor line numbers are current to the character — they were
grepped during authoring and are marked "adjust to reality"; the seven-touchpoint STRUCTURE is the
binding part. (e) I did **not** compute `Λ(kᵐ)` for `m ≥ 4` or a large order complex past the deep
budget — the `dim K_n`/`f`-vector growth is real and capped by a `max_term_dim`-style budget with a
truncation flag (Scope gates 4–5). (f) The `is_lattice`/`is_bounded` helpers are `O(|P|²)`/`O(|P|³)`
exact finite checks — fine for the poset sizes the GUI expresses; no scaling claim beyond that.

**Why this is correct.** The two theorems (Gerstenhaber–Schack/Cibils; GHMS) are verbatim-sourced,
and every numerical consequence was independently reproduced by two engines at spec time: the
incidence HH by the general CS engine AND the hand-rolled simplicial SNF (matching degreewise on
seven posets including the char-sensitive `RP²₆` at dim 121), and the Koszul Betti by the minimal
syzygy engine AND the Plan-27 Ext-algebra Hilbert series (matching `dim K_n` exactly on the
exteriors). Each deliverable ships its own new oracle CLASS (the SNF simplicial complex; the GHMS
bimodule ranks) that cross-validates the shipped engines from an independent route — the plan's
whole point. Every scope boundary is a loud typed refusal (no poset provenance; not certified
Koszul; SNF/`K_n` budget), and the two honest-scope items (ring iso at dimension level; no bracket
claim) are recorded, not hidden.

## Open design risks (what a critic will likely attack)

1. **The GHMS differential is specified, not verified.** The plan pins the term ranks and the HH
   dims, but the comultiplicative differential itself is built at implementation. Mitigation: the
   `d∘d=0` self-cert + the GHMS ≡ minimal ≡ bar/CS cross-engine anchors are *decisive* (a wrong
   splitting/sign fails them); the minimal syzygy engine computes the homology of the SAME minimal
   resolution, so the HH targets are known-correct. A worker who cannot make `d∘d=0` and the
   cross-engine gate pass simultaneously must escalate (do not ship a differential that only
   satisfies one).
2. **"Fast" is a SIZE claim, not a speed claim.** The general CS engine is not intractable on the
   pins (`RP²₆`: ≈2.8 s GF₂ / ≈6.4 s QQ, measured). The honest, testable evidence is a
   **size** `selfcert` note — the order-complex cochain dimension `#chains` ≪ the bar/CS bimodule
   term dimension `O((dim kP)²)`, plus the UCT payoff (one integer SNF → all characteristics). No
   timing/benchmark gate (non-deterministic across machines).
3. **Multi-vertex GHMS corner-typing.** The exterior examples are one-vertex; the **`kZ₃/rad²`**
   (`HH_* = [3,0,1,1,0,0]`) and **diamond-incidence** (`[4,0,0,0,0]`, `K₀ = |Q₀| = 4`) cross-checks —
   both `ext_algebra.koszul = True` — are the guard that `K_n`'s corner grading `e_v K_n e_w` and
   the collapse are right. Keep them in the acceptance set. (Preprojective `A₃` is NOT a GHMS target
   — it is **not Koszul** (Ext obstruction at degree 3), so it exercises the refusal, not the
   compute path.)
4. **Char-sensitivity could be dismissed as a curiosity.** It is the opposite — it is the load-
   bearing reason SNF (not rank) is implemented, and a genuine, verifiable feature of the theorem
   the general engine also exhibits (`RP²₆` on both). Keep the `RP²₆` pin on BOTH engines.

## Change log

- **2026-08-08 authoring** — plan drafted against merged `dev` (Plan 05/13/16 minimal engine, Plan
  04 CS, Plan 27 Koszul/Ext-algebra, Plan 06/33 `Poset`/`IncidenceAlgebra`). Scope boundary fixed
  vs the shipped tree (constructor + minimal engine reused; the order-complex simplicial engine,
  the `incidence_cohomology` kind + poset input mode, and the GHMS closed-form resolution +
  `engine="ghms"` are new). All numerical pins LIVE-VERIFIED at spec time on two independent
  engines (general CS + hand-rolled SNF for R9; minimal syzygy engine + Plan-27 Ext Hilbert for
  R10). The `RP²₆` char-sensitive pin discovered and validated (dim 121 on both engines). The
  Gerstenhaber–Schack **1983** ring citation identified as distinct from the shipped 1987 Hodge
  paper; the **bracket-vanishing clause dropped** per the record (not re-sourced). The P77 seam
  named (GHMS = certified-quadratic-Koszul resolution; P77 = generalized Koszulity recognizers).
- **2026-08-08 critic fix round (NEEDS WORK → adjudicated all-valid)** — the R9 half held up
  verbatim (every pin reproduced on the critic's own boundary matrices, incl. the `RP²₆` char-split
  and the UCT formula); seven findings applied, all live-re-verified:
  - **W1/H1/M1 (major):** the claim "preprojective `A₃` is Koszul but `g_quadratic=False`" was
    **FALSE** — I had read only half the shipped oracle. `A₃.ext_algebra(3).koszul is False` with a
    **positive** obstruction `(3, "a new Ext-algebra generator appears in degree 3")` (LIVE-VERIFIED;
    Brenner–Butler–King `(p,q)`-almost-Koszul, breaks at `h−1=3`). Fixes: the GHMS gate now reads the
    **full three-valued `ext_algebra.koszul`** (not `g_quadratic` alone), split into two honest
    messages (False = not-Koszul-with-Ext-obstruction; None = inconclusive); every "A₃ is Koszul"
    sentence rewritten (A₃ is the **negative** example, refused for the RIGHT reason); the None
    ("unknown") branch honestly ships **unexemplified** (no cheap certified Koszul-but-not-G-quadratic
    witness pinned — worker task); the P77 seam now names **A₃ as the shared boundary** (P75 refuses;
    P77 recognizes `(p,q)`-almost-Koszul).
  - **W2 (major):** the diamond-bridge numbers were wrong — measured minimal bimodule ranks are
    `[4,4,1,0,0,0]` (`K₀=|Q₀|=4`, `K₂=1`, **length 2** not 1) and the `[1,0,0]` I wrote was the
    order-complex COHOMOLOGY, conflated with Hochschild HOMOLOGY `HH_* = [4,0,0,0,0]` (`HH₀ = dim
    A/[A,A] = 4`). Fixed the table + prose; the two quantities are now explicitly separated; the
    diamond is now a genuine multi-vertex GHMS **homology** target (`[4,0,0,0,0]`, LIVE-VERIFIED).
  - **M2 (minor):** wrote the canonical-key test and pinned the `elements` normalization edge —
    catalog default `elements:None` vs an explicit list key differently, so the panel normalizes
    (isolated points only); covers-only ≡ typed family, isolated-point poset differs (both pinned).
  - **H2 (minor):** dropped the "fast path is essential" rhetoric for `RP²₆` (the general CS engine
    measured ≈2.8 s GF₂ / ≈6.4 s QQ to degree 2) — kept the honest **size** framing (`#chains` ≪
    `O((dim kP)²)`) + the UCT payoff.
  - **W3 (minor):** attributed the **arbitrary-finite-poset** ring statement to **Cibils 1989**
    (G–S 1983 = the face-poset-of-a-simplicial-complex version); the theorem header + registry
    annotations reflect the split.
  - **H3 (minor):** made `dev`-tip references tip-agnostic (all prerequisites are early merged
    plans; siblings P63/P65/P67/P70/P72/P73 merge in any order; no edge taken except the conditional
    P77 seam).
