# Plan 65: Exceptional sequences — τ-exceptional (Jasso recursion) + classical hereditary (braid action) (P65 / R27 + R28)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking. Do the reference re-verification in
> the Record/Reference sections BEFORE writing any oracle pin — the flagged
> attributions below (`# PIN`) must be resolved against the sources, not this plan's prose.

**Goal.** Two independent-but-cross-checked exceptional-sequence surfaces, each grounded on
machinery quiverlab already ships:

- **(a — R27) τ-exceptional sequences (Buan–Marsh).** Over a τ-tilting-finite algebra
  `A = kQ/I`, the **signed τ-exceptional sequences** via the **Jasso τ-perpendicular
  reduction** on the shipped P45 τ-tilting engine: a length-1 τ-exceptional object is an
  indecomposable **τ-rigid** module or a shifted indecomposable projective (the *sign*);
  a complete signed τ-exceptional sequence is built recursively by reducing to the
  **τ-perpendicular category** `J(U) ≅ mod C(U)`, an algebra of rank `n − |U|`. Enumeration
  goes through the **bijection with ordered support τ-tilting modules** (Buan–Marsh), NOT
  mutation-BFS — **mutation transitivity is proven only in rank 2** (Buan–Hanson–Marsh
  arXiv:2402.10301), an **honest note** the plan makes contractual. The count is
  cross-checked to `n! · #(support τ-tilting)` from P45.

- **(b — R28) Classical hereditary exceptional sequences (Crawley-Boevey / Ringel).** For a
  **hereditary** `A = kQ` (`Q` acyclic, `I = 0`): an **orthogonality recognizer**
  (`Hom(E_j,E_i)=0=Ext¹(E_j,E_i)` for `i<j`, each `E_i` exceptional), **braid mutation**
  `σ_i` by the canonical universal-extension / kernel / cokernel constructions, **braid-orbit
  BFS** (the braid group `B_n` acts **transitively** on complete exceptional sequences —
  Crawley-Boevey Ottawa 1992; Ringel), and **Dynkin closed-form counts** as literature
  oracles. Hereditary-only; enumeration is representation-finite (Dynkin) only, loud refusal
  otherwise (rep-infinite hereditary has an infinite braid orbit).

- **The cross-check (the two halves meet).** On a **hereditary representation-finite** `A`,
  the τ-exceptional sequences with **all objects genuine modules** (no negative shifts)
  coincide **as a set** with the classical exceptional sequences (Buan–Marsh: τ-exceptional
  generalizes exceptional, coinciding for hereditary). This is a discriminating oracle
  between (a) and (b) — a mismatch is a loud bug in one of them.

**GUI.** ONE new **algebra-only** budget-carrying compute kind `exceptional_sequences`
(the `ar_quiver`/`tau_tilting` precedent — expensive scalar kind, NOT the fast
`recognizers` block), served by **all three tiers** with a **budget** (`exceptional_sequences`
or `exceptional_sequences:512`), the recognizer-panel row, both runners byte-identical, **all
four locales (en/es/fr/zh)** with exact key parity, canonical keys stable, report rendering,
and citations.

**Architecture.** TWO new source modules + one algebra-only kind; everything else is a thin
exact layer over primitives already on `dev` (no new math engine):

- **`src/quiverlab/tautilting/exceptional.py`** (new — part a) — the τ-exceptional surface.
  Lives in `tautilting/` because it is built on the P45 τ-tilting engine (exchange graph,
  τ-rigid pairs, g-vectors, Bongartz completion in the τ sense). Public:
  - `tau_exceptional_objects(A, *, budget=512) -> list[TauExcObject]` — the length-1
    signed τ-exceptional objects (indecomposable τ-rigid modules + shifted projectives).
  - `tau_perpendicular_reduction(A, U, *, budget=512) -> TauReduction` — the **Jasso
    reduction**: `C(U)` (a genuine `kQ'/I'` Algebra of rank `n − |U|`) + the bookkeeping
    to recurse. The reusable primitive.
  - `tau_exceptional_sequences(A, *, budget=4096, want_sequences=True) -> TauExcReport` —
    complete signed τ-exceptional sequences via the ordered-sτ-tilt bijection / reduction
    recursion; τ-tilting-finite gate, loud refusal otherwise; honest budget cap.
  - `is_tau_exceptional_sequence(A, seq) -> bool` — the recognizer.
  - `tau_exceptional_block(A, *, budget=512) -> dict` — the JSON block for the runners.
- **`src/quiverlab/modules/exceptional.py`** (new — part b) — the classical hereditary
  surface (a `modules/` citizen: it knits AR + builds modules via kernels/cokernels/universal
  extensions). Public:
  - `is_exceptional_module(A, M) -> bool` — `End(M)` a division ring (`= k`, alg. closed /
    brick) and `Ext^{≥1}(M,M) = 0` (hereditary ⇒ only `Ext¹`).
  - `is_exceptional_sequence(A, seq) -> bool` — the orthogonality recognizer.
  - `braid_mutation(A, seq, i, *, direction="left") -> list` — `σ_i` / `σ_i^{-1}`.
  - `exceptional_sequences(A, *, budget=100000) -> ExcSeqReport` — Dynkin enumeration
    (direct orthogonality search) + the braid-orbit BFS transitivity certificate + the
    Dynkin closed-form count cross-check; hereditary-only, rep-finite-only, loud otherwise.
  - `c_matrix(A, seq) -> list[list[int]]` — the (signed) c-matrix of a complete sequence
    (hereditary-only; overlaps the P45 c-vectors / wall normals).
- The `exceptional_sequences` kind is wired into `hpc/spec.py::_dispatch` +
  `docs/gui/runner.py` (byte-identical twins) + `webapp/server/schema.py` (the grammar
  parse — a THIRD site, like `tau_tilting`/`ar_quiver`), the GUI touchpoints
  (checkbox / `S.ids` / push-list / `renderBlock` / `scheduleProbe`), the i18n chains, and
  `trace/results_html.py` — exactly following the `ar_quiver` / `tau_tilting` algebra-only
  kinds. One shared block builder dispatches: classical (if hereditary) + τ-exceptional (if
  τ-tilting-finite), each honest about applicability.

**Tech Stack.** Pure exact `Domain` linear algebra + exact combinatorics on the finite
indecomposable universe. **No floats in `src/`** (AST-gated by `tests/test_no_floats.py`):
c-matrices are integer matrices, verdicts/counts are `int`/`bool`/`str`, sequences are
module lists. All homology exact (`Algebra.ext`, `hom_dim`); τ via `Module.tau`; module
constructions via `ModuleHom.{kernel,image,cokernel}` + `modules.yoneda.baer_extension`.

---

## Records (verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R27 — τ-exceptional sequences + mutation.** [D-scout P5; keep-with-note] Object:
> (signed) τ-exceptional sequences via τ-perpendicular recursion (Jasso reduction on the
> shipped τ-tilting engine); bijection with ordered support τ-tilting. NOTE (critic):
> mutation transitivity for enumeration proven only in rank 2 (arXiv:2402.10301) —
> enumeration for rank ≥ 3 via the bijection, not via mutation-BFS. Refs: Buan–Marsh J.
> Algebra 585 (2021) 36–68; arXiv:2211.10428; Buan–Hanson–Marsh arXiv:2402.10301; Nakayama
> counting paper (s10468-021-10060-y) as oracle. Size M.

> **R28 — Classical exceptional sequences (hereditary): braid action, enumeration,
> c-matrices.** [D-scout P6; keep-with-corrections] Object: recognizer (Hom/Ext
> orthogonality), braid mutation σ_i (universal extension/kernel constructions), enumeration
> by braid-orbit BFS (transitive for hereditary — Crawley-Boevey; Ringel), Dynkin
> closed-form counts as oracles; c-matrix reading (overlaps existing τ-tilting c-vectors —
> hereditary-only scope). CORRECTED refs: cite Crawley-Boevey "Exceptional sequences of
> representations of quivers" (Ottawa 1992 proceedings) and Ringel (CMS Conf. Proc. 14)
> directly — arXiv:2102.04584 is Alvares–Marcos–Meltzer on weighted projective lines
> (mislabeled by the scout, dropped as anchor); Garver–Igusa–Matherne–Ostroff
> arXiv:1506.08927; the Dynkin count paper (eudml 284342). Size S–M.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P65 (Wave 3, tier γ,
**independent**). Sizes: (a) M, (b) S–M.

---

## Reference re-verification (done at authoring; findings binding, `# PIN` = resolve at citation-add)

The metaplan standing rule (§1.5) requires the plan writer to re-verify citations. Findings:

1. **Buan, A. B. & Marsh, R. J. — "τ-exceptional sequences", J. Algebra 585 (2021), 36–68.**
   Title, authors, venue verified. **arXiv:1802.01169** (submitted 2018-02-04, revised
   2021-06-02 "To appear in J. Algebra") — **verified this fix round** by fetching the arXiv
   abstract; this is THE load-bearing paper for the ordered-support-τ-tilting bijection.
   **Correction (this fix round):** the earlier draft's `note = {arXiv:2011.02068}` was
   WRONG — arXiv:2011.02068 is an unrelated Coptic-NLP paper (verified); ship the note as
   `arXiv:1802.01169` or drop it (never guess). `# PIN`: exact page range (36–68) and DOI
   (`10.1016/j.jalgebra.2021.05.001`) confirmed against the ScienceDirect landing page at
   citation-add; the **precise statement of the bijection theorem and the recursion
   ordering convention** (which end of the sequence is the reduced one) transcribed
   **verbatim from the PDF** before the recursion is coded (the enumeration correctness
   depends on it — see the Mathematical foundation).

2. **Buan, A. B., Hanson, E. J. & Marsh, R. J. — arXiv:2402.10301** ("Mutation of
   τ-exceptional pairs and sequences" — title verified this fix round via the arXiv listing).
   This is the record's transitivity anchor: **mutation transitivity is proven only in rank
   2** — the honest note that FORBIDS a mutation-BFS enumeration at rank ≥ 3. `# PIN`: verify
   sole content = the rank-2 result and that no rank ≥ 3 transitivity is claimed (read the
   abstract + statement); upgrade `@misc`→ `@article` only if a journal ref exists at merge
   (house rule, no fabricated fields).

3. **Jasso, G. — "Reduction of τ-tilting modules and torsion pairs", Int. Math. Res. Not.
   IMRN 2015, no. 16, 7190–7237.** DOI `10.1093/imrn/rnu163`. This is the τ-perpendicular
   reduction `J(U) ≅ mod C(U)`. `# PIN`: the **concrete construction of `C(U)`** (the
   plan uses the DIJ idempotent-quotient recipe — Bongartz-complete `U`, `B = End_A(T_U)`,
   `C = B/⟨e_U⟩`) is transcribed from Jasso Thm 1.4 **and** Demonet–Iyama–Jasso (already
   shipped as `demonet_iyama_jasso` = DIJ2019) §3, whichever states the algebra-level
   recipe most explicitly; the equivalence `F: mod C(U) → J(U)` is quoted before it is used
   to lift sequence terms.

4. **Crawley-Boevey, W. — "Exceptional sequences of representations of quivers", in
   Representations of algebras (Ottawa, ON, 1992), CMS Conf. Proc. 14, Amer. Math. Soc.,
   1993, pp. 117–124.** The transitivity of the braid action. `# PIN`: venue/pages
   confirmed at citation-add (CMS Conf. Proc. vol. 14).

5. **Ringel, C. M. — "The braid group action on the set of exceptional sequences of a
   hereditary Artin algebra", in Abelian group theory and related topics (Oberwolfach,
   1993), Contemp. Math. 171, Amer. Math. Soc., 1994, pp. 339–352.** **RESOLVED (this fix
   round, critic-adjudicated):** the metaplan card's "CMS Conf. Proc. 14" is **Crawley-
   Boevey's** venue (item 4), NOT Ringel's; Ringel's braid paper is **Contemp. Math. 171
   (1994), 339–352**. The former `# BLOCK` is lifted — ship the Ringel entry with this venue
   (see the `RingelBraid1994` BibTeX in Task G2). The braid transitivity is anchored by BOTH
   Crawley-Boevey (item 4) and Ringel; either alone suffices as the citable statement.

6. **Obaid, Nauman, Al-Shammakh, Fakieh, Ringel — "The number of complete exceptional
   sequences for a Dynkin algebra", Colloq. Math. 133 (2013), no. 2, 197–210** (eudml
   284342). The Dynkin closed-form counts. `# PIN`: author list, volume/number/pages, and
   the closed forms transcribed against the PDF. The **general formula** the plan pins is
   `#CES(Δ) = n! · h^n / |W(Δ)|` (`n` = rank, `h` = Coxeter number, `|W|` = Weyl group
   order), which specialises to `A_n → (n+1)^{n-1}` — **verify this is the paper's stated
   formula** (it is the standard one; the plan live-checked its VALUES for `A_{2..5}` and
   `D_4`, see below).

7. **Igusa–Todorov signed exceptional sequences** (context — arXiv:2509.10910 is re-slotted
   as context per the research doc; the operative signed-exceptional-sequence origin is
   Igusa–Todorov, "Signed exceptional sequences and the cluster morphism category"). `# PIN`:
   add only if load-bearing (it corroborates the `n!·#sτt` count via `n!·#clusters` for
   hereditary — a cross-reference, not a pin). Optional.

8. **Garver–Igusa–Matherne–Ostroff arXiv:1506.08927** and **Buan–Marsh arXiv:2211.10428**
   — supporting refs named in R27/R28. **Correction (this fix round):** arXiv:2211.10428 is a
   **DISTINCT, later** Buan–Marsh paper, "Mutating signed τ-exceptional sequences" (2022;
   verified via the arXiv listing) — it is NOT the arXiv version of the J. Algebra 585 (2021)
   bijection paper (that is arXiv:1802.01169, item 1), and the R27 record's inline
   "arXiv:2211.10428" beside "J. Algebra 585 (2021)" conflates the two. Add 2211.10428 as an
   `@misc` ONLY if cited in an annotation (it is the mutation companion); it is not
   load-bearing for any pin. Optional.

**Already shipped and reused (verified live in `registry.py`):** `air_tau_tilting`
(AIR2014), `demonet_iyama_jasso` (DIJ2019), `king_stability` (King1994), `bongartz_tilting`
(Bongartz1981), `ringel_tame` (Ringel1984tame), `assem_book` (ASS2006). **New keys to add:**
`buan_marsh_tau_exceptional` (arXiv:1802.01169 / J. Algebra 585), `buan_hanson_marsh`,
`jasso_reduction`, `crawley_boevey_exceptional`, `ringel_braid` (item 5 RESOLVED —
Contemp. Math. 171), `obaid_dynkin_count`.

---

## Mathematical foundation (the definitions, exactly — the plan's ground truth)

Throughout `A = kQ/I` is basic, connected, finite-dimensional over an exact `Domain`;
`n = |Q_0| = #simples = rk K_0(A)`; `τ = DTr`, `τ⁻ = TrD` (`Module.tau()`/`.tau_minus()`,
`τX = 0` on projectives); modules are **right** modules (the quiverlab default), composition
**left-to-right**.

### (b) Classical exceptional sequences (hereditary `A = kQ`, `I = 0`)

- **Exceptional module.** `M ∈ mod A` is **exceptional** if `End_A(M)` is a **division
  ring** and `Ext^i_A(M, M) = 0` for all `i ≥ 1`.
  - **Brick vs. exceptional — the recognizer's honest scope (M3).** `end_dim(M) = 1` is the
    **brick** criterion, i.e. `End_A(M) = k·id ≅ k`. Over ANY field a 1-dimensional
    endomorphism ring is `k`, a field, hence a division ring — so `end_dim(M) = 1` ⇒
    `End_A(M)` division ring **always** (no FALSE positives, even over QQ / GF(p)). The
    caveat is the **converse**: over a **non-algebraically-closed** field a module may be
    exceptional with `End_A(M)` a *larger* division ring (`end_dim > 1`, e.g. a field
    extension or a division algebra), and the brick recognizer returns **False** for it (a
    scope-limited FALSE negative). Over an **algebraically closed** field the only
    finite-dimensional division ring is `k`, so **brick ⟺ exceptional-module-with-no-self-ext**
    exactly. `is_exceptional_module` therefore implements the **brick** test
    (`end_dim(M) = 1 and A.ext(M, M, 1) = 0`) and is CONTRACTUALLY the brick criterion, which
    **coincides** with "exceptional" over algebraically closed `k` **and** on the plan's
    Dynkin/QQ battery (every Dynkin indecomposable is a brick). The scope table records this.
  For **hereditary** `A`, `gl.dim ≤ 1`, so only `Ext¹` matters:
  `M` brick-exceptional `⟺ end_dim(M) = 1` and `A.ext(M, M, 1) = 0`. **For rep-finite
  hereditary (Dynkin) `A`, every indecomposable is a brick and exceptional** (rigid brick —
  live-verified: `#indec = #exceptional` for `A_{2..5}`, `D_4`).
- **Exceptional sequence.** An ordered tuple `(E_1, …, E_r)` of exceptional modules with,
  **for all `i < j`**, `Hom_A(E_j, E_i) = 0` **and** `Ext¹_A(E_j, E_i) = 0` (no maps or
  extensions "backwards", from a later to an earlier term). **Convention fixed** (Crawley-
  Boevey/Ringel; the plan's convention — the opposite convention merely reverses each
  sequence, same count; the live checks below use THIS one). **Complete** `⟺ r = n`.
- **Braid action.** The braid group `B_n = ⟨σ_1,…,σ_{n-1}⟩` acts on complete exceptional
  sequences. `σ_i` acts on positions `(i, i+1)`: with `X = E_i`, `Y = E_{i+1}`,
  **left mutation** `σ_i : (…, X, Y, …) ↦ (…, L_X Y, X, …)` where `L_X Y` is defined by the
  canonical exact sequence, four cases on `(dim Hom(X,Y), dim Ext¹(X,Y))`:
  - both `0`: `L_X Y = Y` (orthogonal — a pure swap);
  - `Hom ≠ 0, Ext¹ = 0`: `0 → L_X Y → Hom(X,Y) ⊗_k X → Y → 0` (the **evaluation** map;
    `L_X Y = ker`);
  - `Hom = 0, Ext¹ ≠ 0`: `0 → Y → L_X Y → Ext¹(X,Y) ⊗_k X → 0` (the **universal
    extension**; `L_X Y = middle`);
  - **both `≠ 0` (case d — the concrete two-step, M2).** `L_X Y` is the fibre of the
    canonical evaluation `R Hom(X,Y) ⊗_k X → Y` (the derived/Bondal–Kapranov definition of
    left mutation); Ringel proves that over a **hereditary** algebra this fibre is again a
    module. Realise it as the **ordered** composite **universal extension FIRST, then
    evaluation-kernel** (this order — not the reverse):
    1. **Universal extension** of `Y` by copies of `X` along `Ext¹`:
       `0 → Y → Ỹ → Ext¹(X,Y) ⊗_k X → 0` (`baer_extension` on an `Ext¹(X,Y)` basis). The
       middle `Ỹ` satisfies `Ext¹(X, Ỹ) = 0` and `dim Hom(X, Ỹ) = dim Hom(X, Y)` (apply
       `Hom(X,-)`; `End(X)=k` and the connecting map is iso).
    2. **Evaluation-kernel** through the now-Ext-free `Ỹ`:
       `0 → L_X Y → Hom(X, Ỹ) ⊗_k X → Ỹ → 0`, `L_X Y = ker` of the evaluation.
    In `K₀` this gives `[L_X Y] = ⟨X,Y⟩·[X] − [Y]` with `⟨X,Y⟩ = dim Hom(X,Y) − dim
    Ext¹(X,Y)` the Euler form — the same class as the derived fibre (arithmetic-checked:
    step 1 adds `e·[X]`, step 2 subtracts it back and removes `[Y]`). **Refs:** Ringel,
    Contemp. Math. 171 (1994), 339–352, and Crawley-Boevey (Ottawa 1992) — `# PIN`: the
    exact module-representative statement (the fibre lands in `mod A`, and the sign/shift
    convention when the evaluation is not surjective) transcribed from Ringel before coding.
    The finite-universe `is_isomorphic` cross-identification (Task B2) and the orthogonality
    re-check catch a mis-ordered case-(d) construction loudly on the Dynkin battery.
  `σ_i^{-1}` is the **right mutation** `(…, X, Y, …) ↦ (…, Y, R_Y X, …)` (dual constructions).
  Each `L_X Y` / `R_Y X` is exceptional and `(…, L_X Y, X, …)` is again a complete
  exceptional sequence.
- **Worked braid mutation on `kA₂` (H4 — derived BY HAND, independent of the case table).**
  `A = kA₂ = 1→2` over QQ; indecomposables `S₁ = (1,0)`, `S₂ = (0,1) = P₂`, `P₁ = (1,1)`.
  Take the complete exceptional sequence `(E₁, E₂) = (P₁, S₁)` — valid because
  `Hom(S₁,P₁)=0` and `Ext¹(S₁,P₁)=0` (live-verified). Apply left mutation `σ₁`, so
  `X = P₁`, `Y = S₁`. Compute `Hom(P₁,S₁) = dim (S₁)₁ = 1 ≠ 0` (evaluation of `P₁` onto its
  top `S₁`) and `Ext¹(P₁,S₁) = 0` (`P₁` projective). This is **case (b)** (`Hom≠0, Ext¹=0`),
  so `L_X Y = ker` of the evaluation `Hom(P₁,S₁) ⊗ P₁ = P₁ → S₁`, i.e. the short exact
  sequence
  `0 → L_{P₁}S₁ → P₁ → S₁ → 0`, whose kernel is `rad P₁ = soc P₁ = S₂`, dim-vector `(0,1)`
  (`identify_standard = ('simple', 2)`). Hence
  `σ₁(P₁, S₁) = (L_{P₁}S₁, X) = (S₂, P₁)`,
  which is again a complete exceptional sequence (`Hom(P₁,S₂)=0`, `Ext¹(P₁,S₂)=0`
  live-verified). **All numbers live-verified this fix round** (QQ). Consistency check with
  transitivity: the three classical CES of `kA₂` — `(P₁,S₁)`, `(S₂,P₁)`, `(S₁,S₂)` — form a
  **single `B₂ = ⟨σ₁⟩` orbit of size 3 = `#CES(A₂)`** (verified via the derived-triangle
  mutation). NOTE — the OTHER steps of this orbit are not all the same clean case: e.g.
  mutating `(S₂,P₁)` has `Hom(S₂,P₁)=1` but the evaluation `S₂ → P₁` is the socle **inclusion**
  (INJECTIVE, not surjective), so its `L_XY` is the **cokernel** `S₁` (the module
  representative of the derived fibre), NOT `ker=0` — the case-(b) `L=ker` SES presentation
  presumes a surjective evaluation. The worked step above (`σ₁(P₁,S₁)`) is the clean surjective
  case; the injective sub-case is Ringel's dual construction (the implementation's Task-B2
  finite-universe `is_isomorphic` re-identification handles both).
- **Transitivity (Crawley-Boevey; Ringel).** `B_n` acts **transitively** on complete
  exceptional sequences of any hereditary `A`. For **Dynkin** `A` the set is **finite**;
  for rep-infinite hereditary it is **infinite** (infinitely many exceptional
  preprojectives) — so **enumeration is Dynkin-only** (the braid action is still defined,
  but the orbit is infinite — loud refusal for enumeration).
- **Dynkin closed-form counts (Obaid et al.).** `#CES(Δ) = n! · h^n / |W(Δ)|`. Specialises:
  `A_n → (n+1)^{n-1}`; `D_4 → 4!·6^4/192 = 162`. **Live-verified values** (this doc): see
  the table.
- **c-matrix (hereditary-only).** For a complete exceptional sequence `(E_1,…,E_n)` the
  **c-matrix** is the `n × n` integer matrix of the (signed) dimension vectors of the `E_i`
  in the basis of simples — overlapping the P45 **c-vectors** (wall normals / brick
  dim-vectors of the exchange graph). Honest-scope: **hereditary only** (the general
  τ-tilting c-vector already ships in P45; this is the exceptional-sequence reading of it).

### (a) τ-exceptional sequences (Buan–Marsh, any τ-tilting-finite `A`)

- **τ-rigid.** `M` is τ-rigid `⟺ Hom_A(M, τM) = 0` (P45 `is_tau_rigid`). A **support
  τ-rigid pair** `(M, P)`: `M` τ-rigid, `P` projective, `Hom(P, M) = 0`.
- **Signed τ-exceptional object** (length-1). Either `(M, 0)` with `M` an **indecomposable
  τ-rigid module** (positive sign), or `(0, P_v)` = the **shifted projective** `P_v[1]`
  (negative sign, one per vertex). These are the length-1 signed τ-exceptional sequences.
- **Jasso τ-perpendicular reduction (the recursion engine).** For a τ-rigid `U` (basic),
  Jasso's **τ-perpendicular category** `J(U)` is a wide subcategory equivalent to `mod C(U)`
  for an algebra `C(U)` of rank `n − |U|`. **Concrete construction (DIJ §3 idempotent
  quotient — the plan's recipe, `# PIN` transcribe before coding):**
  1. **Bongartz-complete** `U` in the τ sense to a support τ-tilting module `T_U = U ⊕ B_U`
     (the maximal support τ-tilting with `U ∈ add T_U`; in the rep-finite case this is a
     vertex of the P45 exchange graph containing `U`'s summands — computable).
  2. `B = end_algebra(T_U)` (structure-constant), presented `B̂ = presented_form(B)`
     (Gabriel `kQ'/I'`), with the vertices of `B̂` corresponding to `U`'s summands
     identified.
  3. `C(U) = B̂.quotient_by_idempotent(U-vertices)` (`= B̂ / B̂ e_U B̂`), rank `n − |U|`,
     a genuine `kQ''/I''` Algebra. **Self-certified:** `rk C(U) = n − |U|`, AND a **concrete
     iso pin** (H5): for `kA₃ = 1→2→3` and `U = S₁`, `C(S₁) ≅ kA₂` —
     **live-verified this fix round** (`S₁^⊥ = {X : Hom(S₁,X)=0, Ext¹(S₁,X)=0} = {X : X_a`
     iso`}` has exactly **3 indecomposables** `{P₁=(1,1,1), I₂=(1,1,0), S₃=(0,0,1)}`;
     3 indecomposables ⟹ the rank-2 reduction is **connected**, hence `kA₂` and not the
     disconnected `k×k` — which has only 2). The rank identity `rk C(U) = n − |U|` **alone**
     is insensitive to picking a **wrong** Bongartz completion or the wrong `U`-vertices
     (any single-vertex idempotent quotient of a rank-3 algebra has rank 2); the iso pin
     `C(S₁) ≅ kA₂` is what discriminates the correct completion + vertex identification from
     a plausible-but-wrong one, so the self-cert must assert BOTH.
  The equivalence `F: mod C(U) → J(U) ⊆ mod A` transports objects back to the ambient
  category (needed to lift deep sequence terms — see the honest-scope boundary).
- **Shifted-projective (negative-sign) reduction (H2 — the recipe applies ONLY to τ-rigid
  `U`; a shifted outer term reduces differently).** The `C(U) = End(T_U)/⟨e_U⟩` construction
  above requires `U` a **τ-rigid module**. When the outer signed object is a **shifted
  projective** `P_v[1]` (negative sign), its τ-perpendicular category is the module category
  of the **support quotient** `A/⟨e_v⟩`:
  ```
  C(P_v[1]) = A / ⟨e_v⟩ = A.quotient_by_idempotent([v])      # rank n-1, the shipped callable
  ```
  i.e. delete vertex `v` (and every arrow at it) — `J(P_v[1]) = {M ∈ mod A : M_v = 0} ≅
  mod(A/⟨e_v⟩)`, the wide subcategory of `A`-modules supported off `v`. **Live-verified this
  fix round:** for `kA₂`, `A.quotient_by_idempotent([1])` has one vertex `{2}` (dim 1) and
  `A.quotient_by_idempotent([2])` has one vertex `{1}` — rank `n−1 = 1` in each case. The
  reduction recursion, the recognizer `is_tau_exceptional_sequence`, and
  `tau_perpendicular_reduction(A, U)` must **dispatch on the sign of `U`**: τ-rigid module →
  the DIJ End-quotient; shifted projective `P_v[1]` → `A.quotient_by_idempotent([v])`. `F`
  for the shifted case is the tautological inclusion `mod(A/⟨e_v⟩) ↪ mod A` (cheap, no lift
  needed).
- **Complete signed τ-exceptional sequence.** `(M_1, …, M_n)` where (Buan–Marsh; the
  **exact ordering convention** is `# PIN`'d — transcribe from the PDF) one end is a signed
  τ-exceptional object `M` of `A`, and the remaining `n−1` terms form a complete signed
  τ-exceptional sequence of the reduction algebra `C(M)`. Recursion depth `n`.
- **Enumeration = the ordered-sτ-tilt bijection (Buan–Marsh).** There is a bijection
  between complete signed τ-exceptional sequences and **ordered support τ-tilting modules**
  (a support τ-tilting pair — always `n` signed summands — with a total order on its
  summands). Hence
  ```
  #(complete signed τ-exceptional sequences) = n! · #(support τ-tilting pairs).
  ```
  Enumeration walks the P45 exchange graph (τ-tilting-finite ⇒ complete; else **loud
  refusal**) and realises the bijection via the reduction recursion. **NO mutation-BFS**
  (transitivity only rank 2 — R27's honest note; the plan makes it a contractual comment).
- **Cross-check with (b).** On **hereditary rep-finite** `A`, the reduction coincides with
  the classical perpendicular category, and Buan–Marsh's τ-exceptional sequences with **all
  objects genuine modules (no negative shifts)** = the classical exceptional sequences. Set
  equality is a discriminating oracle.

- **The full `kA₂` signed list, BY HAND (H2 — the enumeration made concrete).**
  `A = kA₂ = 1→2`. The **5 length-1 signed τ-exceptional objects**: the 3 τ-rigid
  indecomposables `S₁=(1,0)`, `S₂=(0,1)`, `P₁=(1,1)` (sign `+`) and the 2 shifted
  projectives `P₁[1]`, `P₂[1]` (sign `−`, one per vertex). For EVERY such outer object `M`,
  the reduction `C(M)` has rank `n−|M| = 1` (a τ-rigid `M` gives `C(M) ≅ k` by the DIJ
  quotient; a shifted `P_v[1]` gives `A/⟨e_v⟩ ≅ k`), and a rank-1 algebra `k` has exactly
  **2** signed τ-exceptional objects (its simple `=` projective `G` with sign `+`, and `G[1]`
  with sign `−`). So the count is `5 outer × 2 inner = 10 = n!·#sτt = 2·5`. Writing each
  complete sequence as `(inner, outer)` with the **outer = the reduced-away end** (the exact
  end is `# PIN`'d to Buan–Marsh; reversing the convention just flips each pair — the SET
  below and every shift-count is convention-independent):
  - **3 all-module sequences (no shift)** — genuine-module outer, genuine-module inner.
    These are **exactly the 3 classical CES of `kA₂`** (the (a)↔(b) cross-check;
    `#classical = 3`, live-verified): `(S₁, S₂)`, `(S₂, P₁)`, `(P₁, S₁)`.
  - **7 sequences containing ≥ 1 shifted projective** `P_v[1]` (`= 10 − 3`), split as:
    - **4 with a shifted OUTER** (`outer ∈ {P₁[1], P₂[1]}`, 2 choices × 2 inner). Explicit:
      `C(P₂[1]) = A/⟨e_2⟩ = k` at vertex 1 lifts to `{S₁(+), P₁[1](−)}`, giving
      **`(S₁, P₂[1])`** [1 shift] and `(P₁[1], P₂[1])` [2 shifts]; symmetrically
      `C(P₁[1]) = A/⟨e_1⟩ = k` at vertex 2 lifts to `{S₂(+), P₂[1](−)}`, giving
      `(S₂, P₁[1])` [1 shift] and `(P₂[1], P₁[1])` [2 shifts].
    - **3 with a genuine-module outer but a shifted INNER** (`outer ∈ {S₁, S₂, P₁}`, inner
      `= G[1]`, the shift of `C(outer)`'s projective).
  **Pinned explicitly (H2):** `(S₁, P₂[1])` is a length-2 signed τ-exceptional sequence whose
  outer term is the shifted projective `P₂[1]` and whose inner term is the genuine module
  `S₁` — `is_tau_exceptional_sequence(A, (S₁, P₂[1]))` must return `True`, and the reduction
  of the outer `P₂[1]` is `A.quotient_by_idempotent([2])` (`= k` at vertex 1). Note
  `3 + 7 = 10` and `3 =` the classical count — a self-contained hand check of both the count
  `n!·#sτt` AND the (a)↔(b) unsigned-vs-classical set-equality on `kA₂`.

### The two consequences the plan leans on

1. **`length = n` and enumeration completeness.** A complete (signed) exceptional sequence
   has exactly `n` terms (rank drops by 1 per reduction / one exceptional module per
   `τ`-orbit position). For (b), the direct orthogonality search over the **finite**
   indecomposable universe is exhaustive (rep-finite), so its count is a genuine oracle. For
   (a), the bijection count `n!·#sτt` (P45, cross-engine) certifies completeness even where
   the reduction recursion's per-sequence lift is budget-bounded.
2. **The scope gates make refusals instant.** Non-hereditary ⇒ (b) refuses (`not
   hereditary`). τ-tilting-infinite ⇒ (a) refuses (exchange graph `status="budget"`).
   Rep-infinite ⇒ (b) enumeration refuses (infinite braid orbit; the knit/exchange graph
   caps loudly). `char ≤ dim` over `GF(p)` ⇒ the `decompose`/`is_isomorphic`/`presented_form`
   caveat raises loudly (never a silent wrong verdict).

---

## Live-verified facts (all recomputed in the venv at spec time, QQ)

| algebra | `n` | #indec | #exceptional | **#classical CES (b)** | `n!·h^n/|W|` | **#sτt (P45)** | **`n!·#sτt` = #signed τ-exc (a)** |
|---|---|---|---|---|---|---|---|
| `kA₂` (`1→2`) | 2 | 3 | 3 | **3** | `2·3²/6=3` | **5** | **10** |
| `kA₃` (`1→2→3`) | 3 | 6 | 6 | **16** | `6·4³/24=16` | **14** | **84** |
| `kA₄` | 4 | 10 | 10 | **125** | `(n+1)^{n-1}` | **42** | **1008** |
| `kA₅` | 5 | 15 | 15 | **1296** | `6⁴` | — | — |
| `kD₄` (subspace, `2,3,4→1`) | 4 | 12 | 12 | **162** | `24·6⁴/192=162` | **50** | **1200** |

- **Classical CES counts** (`#classical CES`) computed by **brute-force orthogonality
  enumeration** over the knitted indecomposable universe (the exact algorithm part (b)
  implements): `A_n = (n+1)^{n-1}` confirmed for `n = 2,3,4,5`; `D_4 = 162` confirmed. These
  are **HARD literature pins** (`oracle_literature`).
  - **`A₅ = 1296` runtime (W2 — re-verified this fix round).** The full count-enumeration
    (knit `ar_quiver` → cache the 15×15 Hom/Ext matrices → backtracking DFS over length-5
    tuples) ran in **~3.9 s total** over QQ (3.88 s to build the Hom/Ext cache, 0.01 s for
    the DFS); `A₄ = 125` in ~0.7 s. So `A₅` is CHEAP for the count leg and stays a HARD pin —
    added to the Task B3 parametrize. **Caveat:** the count is cheap; the **braid-orbit BFS
    transitivity** leg for `A₅` (BFS over all 1296 sequences applying `σ_i^{±1}`, each move
    building a module) is heavier and its runtime is to-be-verified-at-implementation — the
    Task B3 test asserts `count == closed_form` for `A₅` unconditionally but gates the
    `transitive is True` assertion for `A₅` behind an opt-in budget (see Task B3).
- **`#sτt`** computed from `tautilting.mutation.exchange_graph` (Catalan `C_{n+1}` for linear
  `A_n`: 5, 14, 42; `D_4` = 50). These are **HARD cross-engine pins** (`oracle_crossengine`:
  P45's `exchange_graph` vertex count vs. the literature Catalan/cluster number of the Dynkin
  type).
  - **`#sτt` is ORIENTATION-INDEPENDENT (H3 — corrected this fix round; the earlier
    "orientation-dependent" claim was FALSE).** `#sτt(Δ)` is the generalized Catalan /
    cluster number of the Dynkin type, the same derived invariant as the classical count —
    it does not depend on the orientation. **Live-verified this fix round:** `#sτt(A₃) = 14`
    for the linear `1→2→3`, the zigzag `1→2←3`, AND `1←2→3`; `#sτt(D₄) = 50` for the
    **subspace** star `2,3,4→1`, the **source** star `1→2,1→3,1→4`, AND the mixed star
    `1→2,3→1,4→1`. The `D_4` table row keeps the subspace orientation only as the concrete
    tested instance, not because the value depends on it. (The individual sequences and the
    modules in them DO differ across orientations — trivially, since the algebras differ —
    but the COUNT `#sτt`, and hence `n!·#sτt`, does not.)
- **`n!·#sτt`** = the Buan–Marsh signed τ-exceptional count.
  - **The identity `signed_count == n!·#sτt` is SELF-CERT, not cross-engine (H1).** The
    implementation DEFINES `signed_count = n! · len(exchange_graph.vertices)`, so a test that
    asserts `rep.signed_count == n!·len(eg.vertices)` checks the formula against itself
    (tautological — it passes even if zero sequences are ever materialised). It is pinned as
    `oracle_selfcert` (the formula is correctly wired + `stt_count` correctly read from P45).
  - **The REAL cross-engine oracle (H1) is MATERIALISATION.** With `want_sequences=True`,
    actually build the complete signed τ-exceptional sequences via the reduction recursion
    and assert `len(materialised) == signed_count`, all pairwise-**distinct**, and each passes
    `is_tau_exceptional_sequence`. Only that makes `n!·#sτt` a genuine cross-check between the
    formula and the enumerator (see Task A3). The `n!` multiplier itself is the **bijection
    theorem** (Buan–Marsh) — the theorem is airtight; the code realising it is what the
    materialisation test pins. Corroboration: `A_2 → 10` also equals Igusa–Todorov's
    `n!·#clusters` (`#clusters(A_2)=5`).

**What I could NOT live-verify (labelled to-be-verified-at-implementation):**
- The **actual τ-exceptional sequences** (the Jasso recursion output) — implementing the
  reduction is the plan's work, not a probe. The COUNT formula `n!·#sτt` is the oracle.
- **Braid-orbit transitivity via BFS** — I verified the total classical counts by direct
  orthogonality enumeration (not by braid moves). The orbit-BFS reaching all of them (the
  transitivity content) is verified at implementation against the pinned total.
- **E-type Dynkin counts** (`E_6 → 720·12^6/51840 = 41472`, etc.) — arithmetic-derived from
  the formula, not brute-forced (universe too large). Pinned as literature (Obaid et al.),
  live-check only `A_{≤5}`/`D_4`.

---

## Scope gates (contractual; every boundary is a loud typed refusal)

| surface | scope | refusal |
|---|---|---|
| (b) `is_exceptional_module` / `is_exceptional_sequence` / `braid_mutation` / `c_matrix` | **hereditary** `A = kQ`, `I = 0` | non-hereditary → `QuiverlabError("… requires a hereditary algebra …")` |
| (b) `is_exceptional_module` **certifies the BRICK criterion** (`end_dim=1 ∧ Ext¹=0`) (M3) | `= exceptional` over **algebraically closed `k`** and on the Dynkin/QQ battery (every Dynkin indec is a brick); over non-alg-closed `k` it is a **sufficient** (brick) test | never a false positive; a non-brick exceptional module (`End` a larger division ring) returns `False` — a scope-limited FALSE negative, honestly documented (does not occur on Dynkin) |
| (b) `exceptional_sequences` (enumeration + orbit BFS) | hereditary **+ representation-finite** (Dynkin) | rep-infinite hereditary → refusal (`status="budget"`; infinite braid orbit) |
| (a) `tau_exceptional_sequences` (enumeration) | **τ-tilting-finite** `A` | τ-tilting-infinite → refusal via `exchange_graph` `status="budget"`; **`exchange_graph` `status="error"` → LOUD refusal too (M1)** — the τ-side gate treats `status ∉ {"complete"}` as unusable and NEVER derives a count from a non-complete graph (an `"error"` graph can carry the right vertex count yet must not be trusted — see the M1 D₄-star defect + Task 0) |
| (a) `tau_perpendicular_reduction` / `is_tau_exceptional_sequence` | any `A` with a quiver presentation (reduction is per-instance) | presentation-less → loud; τ-tilting-infinite Bongartz completion → loud |
| both, over `GF(p)` | **char 0 or char > dim** (the `decompose`/`is_isomorphic`/`presented_form`/`is_tilting_module` caveat) | `char ≤ dim` → loud `QuiverlabError` from the shipped primitives (never a silent wrong verdict) |
| both | quiver-presented `A` (`A.quiver is not None`) | structure-constant-only → loud (like P45/P60) |

All batteries run over **QQ**; a `GF(32003)` parity spot-check is kept only for `kA_n`
(distinct dim-vectors, no `is_isomorphic` raise).

---

## API surface (public via `import quiverlab`; exact only)

```python
# --- part (b): classical hereditary (src/quiverlab/modules/exceptional.py) ---
def is_exceptional_module(A, M) -> bool
def is_exceptional_sequence(A, seq) -> bool                      # seq: list[Module]
def braid_mutation(A, seq, i, *, direction="left") -> list      # sigma_i / sigma_i^{-1}
def exceptional_sequences(A, *, budget=100_000) -> ExcSeqReport
def c_matrix(A, seq) -> list[list[int]]
@dataclass
class ExcSeqReport:
    algebra: object
    n: int
    dynkin_type: str | None          # "A_3" | "D_4" | ...
    sequences: list                  # [[Module, ...], ...] complete CES (budget-capped)
    count: int
    closed_form_count: int | None    # n!*h^n/|W| where the type is classified
    transitive: bool | None          # braid-orbit BFS from seq[0] reached all `count`
    is_complete: bool                # False iff budget-capped / rep-infinite
    status: str                      # "complete" | "budget" | "unsupported" | "error"
    note: str

# --- part (a): tau-exceptional (src/quiverlab/tautilting/exceptional.py) ---
def tau_exceptional_objects(A, *, budget=512) -> list           # [TauExcObject]
def tau_perpendicular_reduction(A, U, *, budget=512) -> TauReduction  # Jasso C(U)
def is_tau_exceptional_sequence(A, seq) -> bool
def tau_exceptional_sequences(A, *, budget=4096, want_sequences=True) -> TauExcReport
@dataclass
class TauExcReport:
    algebra: object
    n: int
    signed_count: int                # = n! * #sTt  (the bijection count)
    stt_count: int                   # #support tau-tilting pairs (P45)
    sequences: list | None           # reduction towers (None if want_sequences=False / budget)
    is_complete: bool                # False iff tau-tilting-infinite / budget
    status: str                      # "complete" | "budget" | "unsupported" | "error"
    note: str

# Algebra delegates (core/algebra.py, thin lazy-import beside is_tilting_module):
Algebra.exceptional_sequences(budget=100_000)                    -> ExcSeqReport
Algebra.tau_exceptional_sequences(budget=4096, want_sequences=True) -> TauExcReport
Algebra.tau_exceptional_objects(budget=512)                      -> list  # [TauExcObject]
Algebra.is_exceptional_sequence(seq)                             -> bool
Algebra.is_tau_exceptional_sequence(seq)                         -> bool
```

`TauExcObject` carries at least `sign ∈ {+1, −1}` (positive = τ-rigid module, negative =
shifted projective), the module (or `None` for a pure shift), and — for a shifted projective —
its `vertex` (so `A/⟨e_vertex⟩` is the reduction; used by `tau_perpendicular_reduction`'s
sign dispatch).

`ExcSeqReport.__bool__` / `TauExcReport.__bool__` are NOT defined (these are data reports,
not predicates) — mirroring `ARQuiver`/`ExchangeGraph`, whose truthiness is never relied on.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P45 (τ-tilting engine), P41 (AR knit), P44 (tilting + `presented_form`), P37
  (`end_algebra`, `hom_dim`, `Module.tau`), P47 (`quotient_by_idempotent`) are merged on
  `dev` — verified by live import at spec time.** This plan consumes, with the exact
  signatures verified in `dev`:
  - P45 `tautilting`: `exchange_graph(A, budget_pairs=512) -> ExchangeGraph` (`.vertices` =
    `[{"pair","g_matrix","label","support","summand_dimvecs","is_initial"}]`, `.arrows`,
    `.adj`, `.is_complete`, `.status ∈ {"complete","budget","error"}`); `SupportTauTiltingPair`
    (`.summands`, `.support`, `.g_key()`, `.g_matrix()`); `make_pair(A, summands, support,
    check=True)`; `initial_pair(A)`; `mutate(pair, k)`; `is_tau_rigid`; `g_columns` (module
    summands' g-vectors, then `-e_v` per support vertex in vertex order); `bricks`,
    `semibricks`, `torsion_class_data`.
  - P41: `Algebra.ar_quiver(budget_modules=256) -> ARQuiver` (`.vertices` =
    `[{"name","dimvec","module"}]`, `.arrows`, `.tau_orbits`, `.is_complete`, `.status`).
  - P37/P44: `modules/hom.py::hom_dim`/`hom_space`/`hom_basis`/`is_isomorphic`/`end_dim`/
    `identify_standard`; `modules/ext.py::ext_dims(A,M,N,top,with_reps=True)` (cocycle
    payload) / `ext`; `modules/morphism.py::ModuleHom.{kernel,image,cokernel,is_iso}` +
    `direct_sum`, `is_direct_summand`; `modules/yoneda.py::baer_extension` (universal
    extension middle from a cocycle); `modules/duality.py::tau`/`Module.tau()`;
    `modules/endomorphism.py::end_algebra`; `core/basic.py::presented_form`/`gabriel_quiver`.
  - P44/P47: `modules/tilting.py::bongartz_completion` (CLASSICAL, pd≤1 — **NOT** the τ
    completion; the τ Bongartz completion is built in Task A1 from the exchange graph);
    `Algebra.quotient_by_idempotent(vertices)`.
  - Invariants: `invariants/recognizers.py::is_hereditary(A)` (needs the quiver);
    `invariants/dynkin_type.py::dynkin_type(quiver) -> ('A',3)`-tuple | `None`.
  Branch `plan-65-exceptional-sequences` off `dev`.
- **Honest semi-decision contract (metaplan §1.3).** Every enumerator is **complete** iff
  its scope gate holds and its budget did not trip; otherwise `status ∈
  {"budget","unsupported","error"}` with a note — **never a guessed count/sequence**.
- **Char scope is load-bearing (inherited P44/P45).** Batteries run over **QQ**; the
  `GF(32003)` parity check is used only for the distinct-dim-vector `kA_n` family. Over
  `char ≤ dim` the shipped primitives raise loudly.
- **No floats in `src/`.** c-matrices are `list[list[int]]`; counts `int`; verdicts
  `bool`/`str`; sequences module lists. Client-side float conversion only (`gui.js`, exempt).
- **Composition is left-to-right** (`a*b` = first `a` then `b`). `_assert_comparable`
  guards every cross-module `hom_dim`/`ext`. All refusals are `QuiverlabError`.
- **Plan-32 markers** (orthogonal; `oracle_*` FORBIDDEN in `tests/{webapp,gui,hpc}` per
  `test_oracle_classes.py`):
  - `oracle_literature`: Dynkin CES counts (`A_{2..5}`, `D_4`); `#exceptional = #indec`
    on Dynkin; the closed-form `n!·h^n/|W|` values.
  - `oracle_crossengine`: (a) the **materialised** τ-sequence count `len(sequences) ==
    signed_count`, pairwise-distinct, each `is_tau_exceptional_sequence` on `A₂`/`A₃`/the
    non-hereditary `k(1↔2)/rad²` (the REAL check that `n!·#sτt` matches the enumerator — H1);
    (a)↔(b) set-equality on hereditary rep-finite; braid-orbit BFS count `==` direct-
    enumeration count `==` closed form; `c_matrix` rows `⊆` P45 c-vectors/wall normals.
  - `oracle_selfcert`: the τ formula identity `signed_count == n!·len(exchange_graph.vertices)`
    is SELF-CERT (tautological, NOT cross-engine — H1); each `braid_mutation` output is again
    a complete exceptional sequence (orthogonality re-checked) and `σ_i∘σ_i^{-1} = id`; each
    `tau_perpendicular_reduction` has `rk C(U) = n − |U|`, `C(U)` a valid Algebra, plus the
    iso pin `C(S₁) ≅ kA₂` for `kA₃` (H5) and the shifted `C(P_v[1]) = A/⟨e_v⟩` (H2); the
    `status ∉ {"complete"}` loud refusal (M1); every enumerated exceptional module
    self-certifies (`end_dim = 1` — the BRICK criterion, M3 — and `Ext¹(E,E) = 0`).
  - `qpa`: QPA 1.37 has **no** exceptional-sequence surface (live `NamesGVars()` guard that
    FAILS if that ever changes — the P35-products precedent); the exceptional MODULES are
    bricks with no self-ext, crosscheckable via QPA `Ext`/`EndOfModule` on the Dynkin zoo
    (an input-level anchor, not a verdict compare).
- **Mid-merge-train counts.** Absolute suite counts drift; the verification task **recounts
  the oracle-class table at merge time** (`tests/release/test_oracle_classes.py`; paste live
  numbers, claim only this plan's deltas).
- Every merge updates `docs/verification.md` (new oracle rows + recounted class table green)
  and adds citations to `references.bib` + `registry.py` (`bibtex()` hard-fails if the two
  disagree). Conventional commits; green at every commit.

**No-red-commits discipline.** Each Task ships a self-contained green slice: recognizers
before enumeration before GUI. Fields default to `None`/empty so intermediate task
boundaries are green (no forward reference, no `xfail`).

---

# Task 0 (BLOCKING for part a) — fix the P45 `mutate` spurious `status="error"` on the D₄ star (M1)

**Why FIRST.** Part (a)'s completeness certificate is `signed_count = n! · #sτt` with `#sτt =
len(exchange_graph(A).vertices)`. That is only trustworthy when the graph is `status="complete"`.
The critic found — and this fix round **reproduced live (QQ)** — a τ-tilting-**finite** input on
which P45's `exchange_graph` returns `status="error"` even though it discovers the CORRECT number
of vertices, so a naive gate that treats `complete`/`budget` as the only outcomes would either
crash or, worse, read a count off an `"error"` graph.

**Reproduced live (this fix round), pins to encode.** `A = kQ`, `Q = D₄` mixed star
`{a:1→2, b:3→1, c:4→1}` over QQ:
- `exchange_graph(A, budget_pairs=2000)` → **`status="error"`, `len(vertices)=50`,
  `is_complete=False`**.
- The 50 vertices are the CORRECT count: the same `Q` on `2,3,4→1` (subspace),
  `1→2,1→3,1→4` (source), and this mixed star all give **50** (`#sτt(D₄)=50`,
  orientation-independent — see H3). So the graph IS genuinely complete; the `"error"` is
  spurious.
- **Root cause (captured live).** Inside `tautilting/mutation.py::exchange_graph`, the BFS calls
  `mutate(pi, k)`; a `QuiverlabError` from any single `mutate` sets `status="error"` and
  `complete=False` (and `continue`s). On this input `mutate` raises at **summands 0 and 3** of
  ONE support-τ-tilting pair (summand dim-vectors `{(0,1,0,0),(1,1,1,0),(1,1,0,1),(1,1,1,1)}`,
  full support), message `"mutate: no valid exchange found at summand … — the 2-term silting
  cone/cocone produced no validating neighbour"`. It is **not** the `char ≤ dim` caveat (this is
  QQ). Those two neighbours EXIST (they are reached from other vertices), so the BFS still finds
  all 50 — `mutate` is simply **incomplete** for these two exchanges on this orientation. This is
  a genuine **P45 defect**, not a Plan-65 issue.

**Files:**
- Modify: `src/quiverlab/tautilting/mutation.py` (`mutate`'s 2-term silting cone/cocone search —
  find the valid neighbour these two exchanges are missing; if a summand genuinely has no
  exchange in a direction it must be reported as a **boundary/no-mutation**, distinct from an
  `"error"`, so `exchange_graph` does not conflate "no mutation here" with "engine failure").
- Test: `tests/modules/test_tau_exceptional.py` (regression: the D₄ mixed star
  `{1→2,3→1,4→1}` has `exchange_graph(...).status == "complete"` and `len(vertices) == 50`,
  matching the two clean orientations) — plus keep every existing P45 mutation test green.

- [ ] **Step 1: Reproduce** the probe above; confirm `status="error"` with 50 vertices, and
  capture the exact failing `mutate(pi, k)` (vertices/summands as pinned).
- [ ] **Step 2: Root-cause** the 2-term silting cone/cocone search in `mutate` for summands 0/3
  of the pinned pair; determine whether a valid exchange partner exists (it does — reached from
  a neighbouring vertex) and why the local search misses it.
- [ ] **Step 3: Fix** so the mutation returns the existing neighbour (turning `status` back to
  `"complete"`); if any exchange is genuinely a support boundary, classify it as no-mutation, not
  `"error"`. Add the regression test; run the full P45 mutation suite.
- [ ] **Step 4: DEFENSIVE gate regardless of the fix.** In part (a) (`tau_exceptional_sequences`
  and `tau_perpendicular_reduction`) treat `exchange_graph(...).status != "complete"` as a **loud
  refusal** (`status="error"`/`"budget"` → honest `TauExcReport` with that status, never a count
  derived from a non-complete graph). This must hold even if a future orientation trips a NEW
  `mutate` edge case (Task A3 depends on it).
- [ ] **Step 5:** commit
  `fix(tautilting): mutate finds the D_4-star exchange neighbours -- exchange_graph no longer spuriously status="error" on a tau-tilting-finite input (Plan 65 Task 0 / M1); regression pins #sTt(D4)=50 all orientations`.

**If Step 2 finds `mutate` is actually correct** (a real support boundary, not a miss): then the
fix is only the classification (`"error"` → a distinct no-mutation status) + Step 4's gate, and
this task documents the actual cause in the commit body instead of a `mutate` math change. (The
live evidence above — 50 vertices, matching two clean orientations — strongly indicates a genuine
miss, not a boundary.)

---

# Task group I — part (b): classical hereditary exceptional sequences (R28)

Ordered first: it needs **no** new math engine (only the shipped AR knit + hom/ext +
kernel/cokernel/universal-extension), the counts are live-verified, and it provides the
cross-check oracle that part (a) leans on.

### Task B1: `is_exceptional_module` + `is_exceptional_sequence` + `c_matrix` recognizers

**Files:** create `src/quiverlab/modules/exceptional.py`; test
`tests/modules/test_exceptional_classical.py`.

- [ ] **Step 1: Failing tests.**
```python
# tests/modules/test_exceptional_classical.py
"""Classical hereditary exceptional sequences (Plan 65 / R28; Crawley-Boevey Ottawa 1992;
Ringel). Recognizer (Hom/Ext orthogonality), braid mutation sigma_i, Dynkin enumeration,
closed-form counts. Hereditary + rep-finite scope; QQ (the knit/decompose char caveat)."""
import pytest
from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.basic import linear_path_algebra
from quiverlab.modules.exceptional import (is_exceptional_module, is_exceptional_sequence,
                                           c_matrix)
lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert

@selfcert
def test_all_dynkin_indecs_are_exceptional():
    A = linear_path_algebra(3, field=QQ)
    for v in A.ar_quiver().vertices:
        assert is_exceptional_module(A, v["module"]) is True   # every A_n indec is a rigid brick

@lit
def test_a2_exceptional_pair_orthogonality():
    A = linear_path_algebra(2, field=QQ)                       # 1->2 : S1, S2, P1
    S1, S2 = A.simple(1), A.simple(2); P1 = A.projective(1)
    # CONVENTION: for i<j need Hom(Ej,Ei)=0 AND Ext^1(Ej,Ei)=0 (no BACKWARD maps).
    # Live-verified (this fix round, QQ): Hom(S2,P1)=1 (the socle inclusion S2 = soc P1 -> P1),
    # Hom(P1,S2)=0, Ext^1(P1,S2)=0. So (S2, P1) IS exceptional (no backward Hom/Ext), while
    # (P1, S2) is NOT (backward Hom(S2,P1)=1). [Earlier draft had this pair reversed.]
    assert is_exceptional_sequence(A, [S2, P1]) is True
    assert is_exceptional_sequence(A, [P1, S2]) is False
    # a complete CES has length n=2; the three A2 CES are (S1,S2),(S2,P1),(P1,S1) -- Task B3.

@selfcert
def test_non_hereditary_refused():
    from quiverlab.errors import QuiverlabError
    from quiverlab.families.radical_square_zero import RadicalSquareZero
    A = RadicalSquareZero(Quiver([1,2,3], {"a":(2,1),"b":(3,2)}), field=QQ)  # kA3/rad^2
    with pytest.raises(QuiverlabError):
        is_exceptional_sequence(A, [A.simple(1)])
```
- [ ] **Step 2:** confirm failure (`ModuleNotFoundError`).
- [ ] **Step 3: Implement.** `_require_hereditary(A)` (via `is_hereditary`, loud else);
  `is_exceptional_module`: `end_dim(M) == 1 and A.ext(M, M, 1) == 0`;
  `is_exceptional_sequence`: each term exceptional AND for `i<j`
  `hom_dim(seq[j], seq[i]) == 0 and A.ext(seq[j], seq[i], 1) == 0` (the fixed convention);
  `c_matrix`: rows = the dimension vectors of the terms (integer, in vertex order — signed
  reading deferred to hereditary c-vectors, matching the P45 wall-normal sign convention;
  document that for module (unshifted) sequences all rows are the non-negative dim-vectors).
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(modules): classical exceptional module/sequence recognizer + c-matrix (Plan 65/R28, hereditary; Crawley-Boevey/Ringel)`.

### Task B2: `braid_mutation` — `σ_i` by the canonical constructions

**Files:** modify `src/quiverlab/modules/exceptional.py`; extend the test.

- [ ] **Step 1: Failing tests** — `σ_i` output is again a complete CES; `σ_i` then
  `σ_i^{-1}` recovers the input (as modules, up to iso via `is_isomorphic`); the mutated
  module is exceptional; dim-vector transforms as expected on a worked `A_2`/`A_3` pair.
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement `braid_mutation(A, seq, i, direction)`** by the four canonical
  cases (Mathematical foundation §b). Primitives:
  - `Hom ≠ 0, Ext¹ = 0`: build the evaluation map `Hom(X,Y) ⊗ X → Y` from `hom_basis(X,Y)`
    (a `ModuleHom` from `direct_sum` of `dim Hom` copies of `X`), `L_X Y = .kernel()`.
  - `Hom = 0, Ext¹ ≠ 0`: `L_X Y =` the universal extension middle via
    `ext_dims(A,X,Y,1,with_reps=True)` cocycles + `baer_extension` (iterate over an
    `Ext¹`-basis / the diagonal cocycle — Crawley-Boevey's universal extension).
  - both `0`: pure swap; both `≠ 0`: the two-step composite.
  - **Robustness cross-identification (rep-finite):** the mutated `L_X Y` is the UNIQUE
    indecomposable in the knitted universe making `(L_X Y, X)` a complete-extendable
    exceptional pair with the right dim-vector; assert the categorical construction's output
    `is_isomorphic` to that universe module (a self-cert that catches a construction bug).
  `direction="right"` = the dual (`R_Y X` via cokernel / the co-universal extension).
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(modules): braid mutation sigma_i via universal-extension/kernel/cokernel (Plan 65/R28); output re-certified as a complete CES + involution`.

### Task B3: `exceptional_sequences` — Dynkin enumeration + braid-orbit BFS transitivity

**Files:** modify `src/quiverlab/modules/exceptional.py`; the `Algebra` delegate in
`core/algebra.py`; extend the test.

- [ ] **Step 1: Failing tests (the HARD literature pins).**
```python
@lit
@pytest.mark.parametrize("n,expected", [(2,3),(3,16),(4,125),(5,1296)])  # (n+1)^(n-1)
def test_dynkin_An_ces_count(n, expected):
    # counts live-verified (this fix round, QQ): A2/3/4/5 = 3/16/125/1296; the A5 COUNT leg
    # (knit + 15x15 Hom/Ext cache + backtracking DFS) runs in ~3.9s -- cheap, stays a HARD pin.
    A = linear_path_algebra(n, field=QQ)
    rep = A.exceptional_sequences()
    assert rep.is_complete and rep.count == expected
    assert rep.closed_form_count == expected          # n! h^n / |W|
    assert rep.dynkin_type == f"A_{n}"
    if n <= 4:
        assert rep.transitive is True                 # braid-orbit BFS reached them all
    # A5 transitivity (BFS over 1296 sequences, each sigma move builds a module) is the heavy
    # leg -- runtime to-be-verified-at-implementation; assert it only if it lands within the
    # deep budget, else gate `transitive` for A5 behind an opt-in budget (keep count+closed_form).

@lit
@pytest.mark.slow   # opt-in: A5 braid-orbit transitivity certificate (may be minutes)
def test_dynkin_A5_transitive_optin():
    A = linear_path_algebra(5, field=QQ)
    rep = A.exceptional_sequences(budget=100_000)
    assert rep.count == 1296 and rep.closed_form_count == 1296
    assert rep.transitive is True                     # if this exceeds budget, mark xfail-slow, not a lie

@lit
def test_dynkin_D4_ces_count():
    A = Quiver([1,2,3,4], {"a":(2,1),"b":(3,1),"c":(4,1)}).algebra(field=QQ)  # subspace D4
    rep = A.exceptional_sequences()
    assert rep.is_complete and rep.count == 162 and rep.closed_form_count == 162
    assert rep.transitive is True

@selfcert
def test_rep_infinite_hereditary_refused():
    Kron = Quiver([1,2], {"a":(1,2),"b":(1,2)}).algebra(field=QQ)   # 2-Kronecker (tame)
    rep = Kron.exceptional_sequences(budget=200)
    assert rep.is_complete is False and rep.status in ("budget", "unsupported")
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement.** `_require_hereditary` + knit `ar_quiver` (rep-infinite ⇒
  `status="budget"`, honest `ExcSeqReport`). **Enumeration (primary):** direct orthogonality
  search — the exact algorithm the probe used: exceptional indecs `E` (`end_dim=1`,
  `ext(·,·,1)=0`), then all length-`n` tuples with the backward-orthogonality filter
  (budget-capped on tuples examined; `A_n`/`D_4` are well within any sane budget). **Braid-
  orbit BFS (transitivity certificate):** from one enumerated CES, BFS under all `σ_i^{±1}`
  (`braid_mutation`, Task B2), dedup by the tuple of dim-vectors + `is_isomorphic`;
  `transitive = (orbit size == count)`. **Closed form:** classify the type via
  `dynkin_type(A.quiver)`, compute `n!·h^n/|W|` from the type (`h`, `|W|` from a small
  per-type table — `A_n: h=n+1,|W|=(n+1)!`; `D_n: h=2n-2,|W|=2^{n-1}n!`; `E_6/7/8` pinned)
  and assert `count == closed_form_count`. `Algebra.exceptional_sequences` delegate added.
- [ ] **Step 4:** run (`tests/modules/test_exceptional_classical.py -v`); **Step 5:** commit
  `feat(modules): Dynkin exceptional-sequence enumeration + braid-orbit transitivity + closed-form n!h^n/|W| count (Plan 65/R28); pins A_{2..5}=3/16/125/1296, D_4=162`.

---

# Task group II — part (a): τ-exceptional sequences (R27)

### Task A1: `tau_perpendicular_reduction` — the Jasso reduction `C(U)` (the recursion engine)

**Files:** create `src/quiverlab/tautilting/exceptional.py`; test
`tests/modules/test_tau_exceptional.py` (deep bucket — mirrors where P45 engine tests live).

- [ ] **Step 1: Failing tests (self-cert).**
```python
# tests/modules/test_tau_exceptional.py
"""tau-exceptional sequences via the Jasso tau-perpendicular reduction (Plan 65 / R27;
Buan-Marsh J. Algebra 585 (2021); Jasso IMRN 2015; DIJ 2019). Enumeration via the ordered
support-tau-tilt bijection (count = n! * #sTt); NO mutation-BFS (transitivity only rank 2,
Buan-Hanson-Marsh 2402.10301). tau-tilting-finite + char-0/char>dim scope, QQ."""
import pytest
from quiverlab.combinat.quiver import Quiver
from quiverlab.fields import QQ
from quiverlab.families.basic import linear_path_algebra
from quiverlab.tautilting.exceptional import tau_perpendicular_reduction
selfcert = pytest.mark.oracle_selfcert

@selfcert
def test_reduction_rank_drops():
    A = linear_path_algebra(3, field=QQ)
    U = A.simple(1)                                   # a tau-rigid indecomposable
    red = tau_perpendicular_reduction(A, U)
    C = red.reduction_algebra
    assert len(list(C.quiver.vertices)) == 3 - 1      # rk C(U) = n - |U|

@selfcert
def test_reduction_is_kA2_not_just_rank2():
    """H5: rk = n-|U| ALONE is insensitive to a wrong Bongartz completion / wrong U-vertices
    (any single-vertex idempotent quotient of a rank-3 algebra has rank 2). Pin the ISO:
    C(S1) of kA3 must be kA2 (connected), not k x k. Live-verified (this fix round): the
    perpendicular category S1^perp = {P1=(1,1,1), I2=(1,1,0), S3=(0,0,1)} has 3 indecs, so the
    rank-2 reduction is connected -> kA2 (k x k has only 2 indecs)."""
    from quiverlab.modules.hom import is_isomorphic
    A = linear_path_algebra(3, field=QQ)
    C = tau_perpendicular_reduction(A, A.simple(1)).reduction_algebra
    kA2 = linear_path_algebra(2, field=QQ)
    assert len(list(C.quiver.arrows)) == 1                    # connected rank-2 = kA2, not k x k
    assert C.dim == kA2.dim                                   # kA2 has dim 3 (1+1+1 for 1->2)
    # (stronger, if a Morita/iso check is cheap: C is isomorphic to kA2 as an algebra)

@selfcert
def test_shifted_projective_reduction_is_support_quotient():
    """H2: the C(U)=End(T_U)/<e_U> recipe is ONLY for tau-rigid U; a shifted projective
    P_v[1] reduces to the support quotient A/<e_v>. Live-verified (this fix round): for kA2,
    A.quotient_by_idempotent([2]) is k at vertex 1 (rank n-1=1)."""
    A = linear_path_algebra(2, field=QQ)
    Pv_shift = next(o for o in A.tau_exceptional_objects() if o.sign < 0 and o.vertex == 2)
    red = tau_perpendicular_reduction(A, Pv_shift)            # dispatches on the sign
    assert list(red.reduction_algebra.quiver.vertices) == [1] # A/<e_2> = k at vertex 1
    assert red.reduction_algebra.dim == A.quotient_by_idempotent([2]).dim
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement the reduction — DISPATCH ON THE SIGN of `U` (H2).**
  - **τ-rigid module `U` (positive sign) — the DIJ recipe** (Mathematical foundation §a):
    τ-Bongartz-complete `U` (rep-finite: scan `exchange_graph` for a pair with `U`'s summands
    ∈ `add`; take the maximal — the initial-pair-side completion); `B = end_algebra(T_U)`;
    `B̂ = presented_form(B)`; identify the `U`-vertices; `C = B̂.quotient_by_idempotent(U-verts)`.
  - **Shifted projective `U = P_v[1]` (negative sign) — the support quotient (H2):**
    `C = A.quotient_by_idempotent([v])` (`= A/⟨e_v⟩`, rank `n−1`), and `F` is the
    tautological inclusion `mod(A/⟨e_v⟩) ↪ mod A` (`{M : M_v = 0}`) — no lift-solve needed.
  Return a `TauReduction(reduction_algebra=C, bongartz=T_U|None, u_vertices=..., equivalence=F)`
  where `F` (the `mod C → J(U)` object transport) is implemented for the **top layer** (the
  cheap direct case; tautological for the shifted case) and documented as the honest-scope
  boundary for deep lifts (see Acceptance). **Self-cert:** `rk C == n − |U|`; `C.quiver is not
  None`; **AND the iso pin `C(S₁) ≅ kA₂` for `kA₃`** (rank alone does not catch a wrong
  completion — H5).
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(tautilting): Jasso tau-perpendicular reduction -- C(U)=End(T_U)/<e_U> for tau-rigid U, A/<e_v> for shifted P_v[1] (Plan 65/R27; DIJ idempotent quotient); rk = n-|U|, C(S1 of kA3)=kA2 pinned`.

### Task A2: `tau_exceptional_objects` + `is_tau_exceptional_sequence`

**Files:** modify `src/quiverlab/tautilting/exceptional.py`; extend the test.

- [ ] **Step 1: Failing tests.** `tau_exceptional_objects(A)` returns the indecomposable
  τ-rigid modules + `n` shifted projectives (count self-cert); `is_tau_exceptional_sequence`
  accepts a hand-built length-1 τ-rigid object and rejects a non-τ-rigid module; on a
  hereditary `A_2` a classical exceptional sequence is accepted (the cross-check seed).
  **PLUS the shifted-projective recognizer pin (H2)** — a complete signed sequence containing
  a `P_v[1]` is recognized:
```python
@selfcert
def test_recognizer_accepts_shifted_projective_sequence():
    """H2 pin: on kA2 the signed sequence (S1, P2[1]) -- inner S1 (genuine module), outer the
    shifted projective P2[1] -- is a complete signed tau-exceptional sequence. Its outer term
    reduces via A/<e_2> = k at vertex 1, whose signed objects lift to {S1(+), P1[1](-)}, so S1
    is a valid inner. (Of kA2's 10 signed sequences, exactly 3 are all-module = the classical
    CES; the other 7 contain >= 1 shifted projective -- see the by-hand list in the Math
    foundation.)"""
    A = linear_path_algebra(2, field=QQ)
    objs = {(o.sign, getattr(o, "vertex", None)): o for o in A.tau_exceptional_objects()}
    P2shift = objs[(-1, 2)]                                  # the shifted projective P_2[1]
    S1 = A.simple(1)
    assert A.is_tau_exceptional_sequence([S1, P2shift]) is True   # (inner, outer) convention
    # and a non-tau-rigid / ill-ordered signed sequence is rejected (loud False), not crash.
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement.** `tau_exceptional_objects`: the τ-rigid indecomposables from the
  exchange-graph module universe (`is_tau_rigid` per indec) tagged sign `+`, plus each
  `A.projective(v)` tagged sign `−` (shift). `is_tau_exceptional_sequence`: the recursive
  check — the outer term is a signed τ-exceptional object of `A`, and the remaining terms
  (transported into `C(outer)` via Task A1's `F`) form a τ-exceptional sequence of `C(outer)`
  (recurse); base case length 0. **`# PIN`:** the exact ordering convention (which end is
  reduced) transcribed from Buan–Marsh before coding.
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(tautilting): tau-exceptional objects + recursive tau-exceptional-sequence recognizer (Plan 65/R27; Buan-Marsh)`.

### Task A3: `tau_exceptional_sequences` — the ordered-sτt bijection enumeration

**Files:** modify `src/quiverlab/tautilting/exceptional.py`; the `Algebra` delegate; extend
the test.

- [ ] **Step 1: Failing tests.**
```python
from math import factorial
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine

@selfcert   # H1: this is SELF-CERT, not cross-engine -- signed_count is DEFINED as n!*stt_count,
            # so asserting it == n!*len(eg.vertices) checks the formula against itself.
@pytest.mark.parametrize("n", [2, 3])            # A_2 -> 10, A_3 -> 84
def test_signed_count_formula_wired(n):
    from quiverlab.tautilting.mutation import exchange_graph
    A = linear_path_algebra(n, field=QQ)
    eg = exchange_graph(A, budget_pairs=2000)
    assert eg.status == "complete"                           # M1: never trust a non-complete graph
    rep = A.tau_exceptional_sequences(budget=4096, want_sequences=False)
    assert rep.is_complete
    assert rep.stt_count == len(eg.vertices)                 # #sTt correctly read from P45
    assert rep.signed_count == factorial(n) * len(eg.vertices)   # formula correctly wired (tautological)

@xeng   # H1: the REAL cross-engine oracle -- MATERIALISE the sequences and count them.
@pytest.mark.parametrize("factory,n", [
    ("A2", 2),                                               # signed_count = 2! * 5  = 10
    ("A3", 3),                                               # signed_count = 3! * 14 = 84 (deeper)
    ("nonhered", 2),                                         # k(1<->2)/rad^2 : 2! * 6 = 12, NON-hereditary tau-tilting-finite
])
def test_materialised_count_equals_signed_count(factory, n):
    """The formula n!*#sTt is only cross-checked when the enumerator actually BUILDS the
    sequences and their number matches. #sTt live-verified this fix round: A2=5, A3=14,
    k(1<->2)/rad^2 = 6 (non-hereditary, status=complete)."""
    if factory == "A2":  A = linear_path_algebra(2, field=QQ)
    elif factory == "A3": A = linear_path_algebra(3, field=QQ)
    else:
        from quiverlab.families.radical_square_zero import RadicalSquareZero
        A = RadicalSquareZero(Quiver([1,2], {"a":(1,2),"b":(2,1)}), field=QQ)  # non-hereditary ttf
    rep = A.tau_exceptional_sequences(budget=4096, want_sequences=True)
    assert rep.is_complete and rep.status == "complete"
    mats = rep.sequences
    assert len(mats) == rep.signed_count                     # materialised == n!*#sTt (the cross-check)
    # pairwise-DISTINCT (dedup by termwise (dimvec, sign) key + is_isomorphic on ties)
    assert len({_sig_key(s) for s in mats}) == len(mats)
    # and every materialised sequence is genuinely recognized
    assert all(A.is_tau_exceptional_sequence(s) for s in mats)

@selfcert
def test_tau_tilting_infinite_refused():
    Kron = Quiver([1,2], {"a":(1,2),"b":(1,2)}).algebra(field=QQ)
    rep = Kron.tau_exceptional_sequences(budget=200)
    assert rep.is_complete is False and rep.status in ("budget", "unsupported")

@selfcert
def test_exchange_graph_error_status_refused_loudly():
    """M1: a status='error' exchange graph must NOT yield a count (even when its vertex count
    is right). Until the Task-0 mutate fix lands, the D4 mixed star trips status='error'; the
    tau-side gate must refuse loudly rather than derive signed_count from it."""
    A = Quiver([1,2,3,4], {"a":(1,2),"b":(3,1),"c":(4,1)}).algebra(field=QQ)  # D4 mixed star
    rep = A.tau_exceptional_sequences(budget=4096)
    # AFTER Task 0: status='complete', signed_count == 24*50 = 1200.
    # BEFORE Task 0 (defensive gate): status='error'/'unsupported', is_complete False, no count.
    assert (rep.is_complete and rep.signed_count == 24 * 50) or (
        rep.is_complete is False and rep.status in ("error", "unsupported"))
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement.** `exchange_graph` — **GATE on `status` (M1):** only
  `status="complete"` yields a count; `status ∈ {"budget","error"}` ⇒ honest `TauExcReport`
  with that status and NO `signed_count`/`stt_count` derived from the graph (a `status="error"`
  graph can carry the right vertex count yet must never be trusted — the Task-0 D₄-star defect).
  On a complete graph: `stt_count = len(eg.vertices)`; `signed_count = n!·stt_count` (the
  bijection — the primary count oracle, computed WITHOUT enumerating, so it is trustworthy even
  when `want_sequences` is budget-capped). With `want_sequences=True`: realise the bijection by
  the reduction recursion (Task A1/A2) — each complete signed τ-exceptional sequence as a
  **reduction tower** `[(A, M_n, sign), (C_1, M_{n-1}, sign), …]`; a hard `budget` cap on towers
  materialised (honest `status="budget"` if tripped, but `signed_count` still exact via the
  formula). **The materialised list is the H1 cross-check** — when it completes within budget,
  `len(sequences) == signed_count`, all distinct, each `is_tau_exceptional_sequence` (Step-1
  `test_materialised_count_equals_signed_count`). `Algebra.tau_exceptional_sequences` delegate.
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(tautilting): tau-exceptional-sequence enumeration via the ordered-sTt bijection (Plan 65/R27); signed_count = n!*#sTt (self-cert formula) cross-checked by MATERIALISATION (A2/A3 + non-hereditary); status-gated (M1); tau-tilting-finite gate`.

### Task A4: the (a)↔(b) cross-check on hereditary rep-finite algebras

**Files:** extend `tests/modules/test_tau_exceptional.py` (an `oracle_crossengine` battery).

- [ ] **Step 1: The discriminating oracle.**
```python
@xeng
@pytest.mark.parametrize("n", [2, 3])
def test_tau_exc_unsigned_equals_classical(n):
    """On hereditary rep-finite A, tau-exceptional sequences with all objects genuine
    modules (no negative shifts) == classical exceptional sequences (Buan-Marsh: coincide
    for hereditary). Set equality (up to is_isomorphic termwise), a loud bug if they differ."""
    A = linear_path_algebra(n, field=QQ)
    classical = A.exceptional_sequences()                     # part (b)
    tau = A.tau_exceptional_sequences()                       # part (a)
    unsigned = _unsigned_module_sequences(tau)               # towers with all-positive signs, lifted
    assert _same_set_of_sequences(unsigned, classical.sequences)   # termwise is_isomorphic
```
- [ ] **Step 2–3:** implement the two helpers (`_unsigned_module_sequences` lifts the
  all-positive reduction towers to ambient A-modules via Task A1's `F` on hereditary — where
  `F` is the classical perpendicular equivalence and cheap; `_same_set_of_sequences` compares
  as sets, termwise `is_isomorphic`). **Honest-scope:** if `F`'s deep lift is not fully
  realised for a non-hereditary case, this cross-check is **hereditary-only** by design (the
  general case is covered by the `signed_count == n!·#sτt` count oracle) — state this on the
  verification page.
- [ ] **Step 4:** run; **Step 5:** commit
  `test(tautilting): (a)<->(b) cross-check -- unsigned tau-exceptional == classical CES on hereditary rep-finite (Plan 65)`.

---

# Task group III — GUI / webapp + verification

### Task G1: the `exceptional_sequences` algebra-only kind (all three tiers, i18n ×4)

An expensive **algebra-only** budget-carrying scalar kind (the `ar_quiver`/`tau_tilting`
precedent). The GUI touchpoints follow `tau_tilting` verbatim.

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (grammar parse — the `tau_tilting`/`ar_quiver`
  budget-suffix branch — + `_dispatch` branch + `_snip` recipe), `docs/gui/runner.py` (the
  Pyodide twin — matching grammar + dispatch + `_snip` + ETA, shape-identical),
  **`webapp/server/schema.py`** (the THIRD grammar-parse site — mirror the `tau_tilting`
  branch), `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox `qlgui-exceptional-
  sequences`, `S.ids`, push-list, `renderBlock`, `scheduleProbe`), `webapp/templates/
  index.html` (one checkbox), `webapp/server/i18n/{en,es,fr,zh}.json` (**all four locales**
  — the block key chain + `pick.kind.exceptional_sequences`; exact key parity gated by
  `tests/webapp/test_i18n.py`), `src/quiverlab/trace/results_html.py` (`_HEADINGS` + a
  render branch), `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_
  delegation.py` (ONE golden; existing byte-identical).
- Test: `tests/webapp/test_exceptional_kind_p65.py`, `tests/gui/test_exceptional_runner_twin_p65.py`.

**Block shape** (one shared `exceptional_sequences_block(A, budget)` dispatching both halves):
```python
{"kind": "exceptional_sequences", "n": int,
 "hereditary": bool, "dynkin_type": str|None,
 "classical": {"count": int, "closed_form_count": int|None, "transitive": bool|None,
               "complete": bool, "status": str} | None,          # None iff not hereditary
 "tau": {"signed_count": int, "stt_count": int, "complete": bool, "status": str} | None,
                                                                  # None iff tau-tilting-infinite
 "note": str|None,
 "references": ["buan_marsh_tau_exceptional", "crawley_boevey_exceptional",
                "jasso_reduction", "air_tau_tilting", "assem_book"]}
# refusal (presentation-less) -> {"error": msg, "references": [...]}, never a 500.
```
Sequences themselves are **NOT** shipped in the block by default (they explode — `D_4` has
162 classical / 1200 signed); the block reports **counts + status**, and the worked-steps
report (`results_html`) renders a small sample (first ≤ 8) with a stated elision.

- [ ] **Step 1: Failing cross-runner tests** (unmarked — extras-gated dir): `kA_3` block
  has `classical.count == 16`, `classical.closed_form_count == 16`, `tau.signed_count == 84`,
  `tau.stt_count == 14`, `"buan_marsh_tau_exceptional"` in the citation keys; twin parity
  (`json.dumps(sort_keys=True)` equality across `hpc/spec.py` and `docs/gui/runner.py`).
- [ ] **Step 2: Implement** the grammar parse in **all three** sites (budget suffix, skips
  the `name:0..N` degree grammar — copy the `tau_tilting` branch exactly), the `_dispatch`
  branch (`from quiverlab...exceptional_sequences_block`), the twin, the `gui.js` checkbox +
  `scheduleProbe` ETA (size like `tau_tilting`/`ar_quiver` — algebra-dim `sizing_dim`, **no
  estimator edit**), the `_snip` recipe
  (`"exceptional_sequences": "A.exceptional_sequences(); A.tau_exceptional_sequences()"`),
  `_HEADINGS["exceptional_sequences"] = "Exceptional sequences"` + the render branch (two
  count tables — classical / τ — each with its status, the closed-form column, the
  transitivity flag; the small sequence sample with elision).
- [ ] **Step 3: Add ONE golden** (`exceptional_sequences_kA3`) to `_runner_goldens.json`;
  note it in `test_runner_delegation.py`'s change-log docstring; confirm existing goldens
  byte-identical first. **Canonical-key stability:** schema-v1 algebra-only, no `module`
  block → canonicalizes unchanged (confirm an existing family request's key is byte-stable).
- [ ] **Step 4:** run `tests/webapp/test_exceptional_kind_p65.py
  tests/webapp/test_runner_delegation.py tests/gui/test_exceptional_runner_twin_p65.py
  tests/webapp/test_i18n.py tests/hpc -q`.
- [ ] **Step 5:** commit
  `feat(gui,webapp,hpc): exceptional_sequences algebra-only kind (Plan 65) -- classical + tau counts, both runners byte-identical, grammar in 3 sites, i18n x4, one golden`.

**i18n keys — ALL FOUR LOCALES** (mirror `pick.kind.tau_tilting`):

| key | en | es | fr | zh |
|---|---|---|---|---|
| `pick.kind.exceptional_sequences` | Exceptional sequences | Sucesiones excepcionales | Suites exceptionnelles | 例外序列 |
| `inv.exceptional_sequences` | Exceptional sequences (classical + τ) | Sucesiones excepcionales (clásicas + τ) | Suites exceptionnelles (classiques + τ) | 例外序列 (经典 + τ) |
| `block.exceptional_sequences.title` | Exceptional sequences | Sucesiones excepcionales | Suites exceptionnelles | 例外序列 |
| `block.exceptional_sequences.classical` | Classical (hereditary, braid orbit) | Clásicas (hereditaria, órbita de trenza) | Classiques (héréditaire, orbite de tresse) | 经典 (遗传, 辫群轨道) |
| `block.exceptional_sequences.tau` | τ-exceptional (signed, ordered sτ-tilt) | τ-excepcionales (con signo) | τ-exceptionnelles (signées) | τ-例外 (带符号) |
| `block.exceptional_sequences.count` | Count | Cantidad | Nombre | 数量 |
| `block.exceptional_sequences.transitive` | Braid-orbit transitive | Órbita de trenza transitiva | Orbite de tresse transitive | 辫群轨道可迁 |

### Task G2: verification page, citations, README, suite gate

**Files:** `src/quiverlab/citations/references.bib` + `registry.py`; `docs/verification.md`;
`README.md`; `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick P65); existing release gates.

- [ ] **Step 1: Citations (BibTeX-VERIFIED only).** Add to `references.bib`. **Resolved this
  fix round (W3/W4):** the `BuanMarsh2021` arXiv id is **arXiv:1802.01169** (NOT the earlier
  draft's `2011.02068`, which is an unrelated Coptic-NLP paper — verified via the arXiv
  abstract); **arXiv:2211.10428** is a DISTINCT later Buan–Marsh paper ("Mutating signed
  τ-exceptional sequences", 2022), not the 585 (2021) bijection paper; the Ringel-braid venue
  is **Contemp. Math. 171 (1994), 339–352** (item 5 RESOLVED — ship it). The remaining
  `doi`/`pages`/`number` fields stay `# PIN`'d against the source PDFs at citation-add
  (`Jasso2015`/`Obaid2013` volume data); the house rule forbids shipping an unverified field —
  delete any that cannot be confirmed rather than guess:
```bibtex
@article{BuanMarsh2021,
  author  = {Buan, Aslak Bakke and Marsh, Robert J.},
  title   = {{$\tau$}-exceptional sequences},
  journal = {J. Algebra}, volume = {585}, pages = {36--68}, year = {2021},
  doi = {10.1016/j.jalgebra.2021.05.001}, note = {arXiv:1802.01169}}   % arXiv id verified this fix round
@article{Jasso2015,
  author = {Jasso, Gustavo}, title = {Reduction of {$\tau$}-tilting modules and torsion pairs},
  journal = {Int. Math. Res. Not. IMRN}, year = {2015}, number = {16},
  pages = {7190--7237}, doi = {10.1093/imrn/rnu163}}
@incollection{CrawleyBoevey1993,
  author = {Crawley-Boevey, William},
  title = {Exceptional sequences of representations of quivers},
  booktitle = {Representations of algebras (Ottawa, ON, 1992)},
  series = {CMS Conf. Proc.}, volume = {14}, pages = {117--124},
  publisher = {Amer. Math. Soc.}, year = {1993}}
@incollection{RingelBraid1994,
  author = {Ringel, Claus Michael},
  title = {The braid group action on the set of exceptional sequences of a hereditary
           {A}rtin algebra},
  booktitle = {Abelian group theory and related topics (Oberwolfach, 1993)},
  series = {Contemp. Math.}, volume = {171}, pages = {339--352},
  publisher = {Amer. Math. Soc.}, year = {1994}}          % venue RESOLVED this fix round (item 5)
@misc{BuanHansonMarsh2024,
  author = {Buan, Aslak Bakke and Hanson, Eric J. and Marsh, Robert J.},
  title = {Mutation of {$\tau$}-exceptional pairs and sequences},
  year = {2024}, note = {arXiv:2402.10301}}               % title verified this fix round
@article{Obaid2013,
  author = {Obaid, Mustafa A. A. and Nauman, S. Khalid and Al-Shammakh, Wafaa S. M.
            and Fakieh, Wafaa M. and Ringel, Claus Michael},
  title = {The number of complete exceptional sequences for a {D}ynkin algebra},
  journal = {Colloq. Math.}, volume = {133}, number = {2}, pages = {197--210}, year = {2013}}
% Optional (annotation-only, NOT load-bearing): the DISTINCT later Buan-Marsh mutation paper.
@misc{BuanMarshMutating2022,
  author = {Buan, Aslak Bakke and Marsh, Robert J.},
  title = {Mutating signed {$\tau$}-exceptional sequences},
  year = {2022}, note = {arXiv:2211.10428}}
```
  and in `registry.py` (`_r(key, bibtex_key, kind, title, annotation, *tags)`):
```python
_r("buan_marsh_tau_exceptional", "BuanMarsh2021", "algorithm",
   "tau-exceptional sequences",
   "Buan-Marsh: signed tau-exceptional sequences via the Jasso tau-perpendicular "
   "reduction; bijection with ordered support tau-tilting modules (#signed = n!*#sTt).",
   "tau-tilting", "modules"),
_r("jasso_reduction", "Jasso2015", "foundation",
   "Reduction of tau-tilting modules and torsion pairs",
   "Jasso: the tau-perpendicular category J(U) ~ mod C(U), rank n-|U| (the DIJ idempotent "
   "quotient) -- the recursion engine for tau-exceptional sequences.", "tau-tilting"),
_r("crawley_boevey_exceptional", "CrawleyBoevey1993", "foundation",
   "Exceptional sequences of representations of quivers",
   "Crawley-Boevey: the braid group acts transitively on complete exceptional sequences "
   "of a hereditary algebra (Ottawa 1992).", "modules"),
_r("ringel_braid", "RingelBraid1994", "foundation",
   "The braid group action on the set of exceptional sequences of a hereditary Artin algebra",
   "Ringel: the braid B_n action on complete exceptional sequences (transitive), with the "
   "sigma_i mutation constructions -- the case-(d) two-step module realization. Contemp. "
   "Math. 171 (1994) (venue resolved).", "modules"),
_r("buan_hanson_marsh", "BuanHansonMarsh2024", "foundation",
   "Mutation of tau-exceptional pairs and sequences",
   "Buan-Hanson-Marsh: mutation transitivity proven only in rank 2 -- why enumeration "
   "at rank >= 3 goes through the ordered-sTt bijection, not mutation-BFS.", "tau-tilting"),
_r("obaid_dynkin_count", "Obaid2013", "foundation",
   "The number of complete exceptional sequences for a Dynkin algebra",
   "Obaid et al.: #CES(Delta) = n! h^n / |W|; A_n = (n+1)^{n-1}, D_4 = 162 -- the "
   "closed-form count oracle.", "modules"),
```
  Verify BibTeX presence (`tests/citations/`) and that every key resolves.
- [ ] **Step 2: Verification page.** Add the P65 subsystem rows:
  - `modules/exceptional.py` — `oracle_literature` (Dynkin CES counts `A_{2..5}` incl.
    `A_5=1296` (~4 s count leg), `D_4`; `#exceptional = #indec`; Obaid et al. closed form);
    `oracle_crossengine` (braid-orbit BFS count == direct-enumeration count == closed form;
    `c_matrix` rows ⊆ P45 c-vectors); `oracle_selfcert` (`σ_i` output re-certified CES +
    involution incl. the worked `σ₁(P₁,S₁)=(S₂,P₁)` on `kA₂` — H4; enumerated modules are
    rigid bricks — `is_exceptional_module` is the BRICK criterion, scope note (i)/M3).
  - `tautilting/exceptional.py` — `oracle_crossengine` (the **materialised** sequence count
    `len(sequences) == signed_count`, pairwise-distinct, each `is_tau_exceptional_sequence`, on
    `A₂` / `A₃` / the non-hereditary `k(1↔2)/rad²` — the REAL check that `n!·#sτt` matches the
    enumerator, H1; the (a)↔(b) unsigned-vs-classical set-equality on hereditary rep-finite);
    `oracle_selfcert` (**the formula identity `signed_count == n!·len(exchange_graph.vertices)`
    is SELF-CERT (tautological), NOT cross-engine — H1**; `rk C(U) == n − |U|` **plus** the iso
    pin `C(S₁) ≅ kA₂` for `kA₃` — H5; the shifted-projective reduction `C(P_v[1]) = A/⟨e_v⟩` —
    H2; the `status="error"`/`"budget"` loud refusal (never a count off a non-complete graph) —
    M1; τ-tilting-infinite loud refusal; recognizer recursion base case).
  - `qpa` — `NamesGVars()` guard that FAILS if QPA ever ships an exceptional-sequence
    surface; the exceptional MODULES crosschecked as bricks-with-no-self-ext via QPA
    `Ext`/`EndOfModule` on the Dynkin zoo (input-level anchor).
  - **Honest-scope entries (binding):**
    (a) **Mutation transitivity is proven only in rank 2** (Buan–Hanson–Marsh) — τ-exceptional
    enumeration goes through the **ordered-sτt bijection**, NOT mutation-BFS; the count
    `n!·#sτt` is the completeness certificate (theorem + P45 cross-engine).
    (b) The **`n!` multiplier of the τ count is the Buan–Marsh theorem**. The formula identity
    `signed_count == n!·#sτt` is SELF-CERT (tautological — it defines the count); the genuine
    cross-check is the **materialised** enumeration (`len(sequences) == signed_count`,
    distinct, each recognized) on `A₂` / `A₃` / the non-hereditary `k(1↔2)/rad²` (H1). Corrob:
    `A_2 = 10` also equals Igusa–Todorov's `n!·#clusters`.
    (c) The **(a)↔(b) set-equality cross-check is hereditary-only** (the general reduction's
    deep object-lift `F` is the scope boundary); the general case rests on the count oracle.
    (d) **Classical enumeration is Dynkin-only** (rep-infinite hereditary ⇒ infinite braid
    orbit, loud refusal); **τ enumeration is τ-tilting-finite-only** (loud otherwise).
    (e) **Char scope:** batteries over QQ; `char ≤ dim` refuses loudly (the shipped
    `decompose`/`is_isomorphic`/`presented_form` caveat).
    (f) **`#sτt` (hence the τ count) is ORIENTATION-INDEPENDENT** (corrected — the earlier
    "orientation-dependent" claim was FALSE): `#sτt(Δ)` is the generalized Catalan/cluster
    number of the Dynkin type, the same derived invariant as the classical count.
    Live-verified this fix round: `#sτt(A₃) = 14` for the linear, zigzag, and `1←2→3`
    orientations; `#sτt(D₄) = 50` for the subspace, source, and mixed stars. The `D_4` table
    row keeps the subspace orientation only as the concrete tested instance. (The individual
    sequences differ across orientations — trivially, different algebras — but not the count.)
    (g) **E-type Dynkin counts** are literature-pinned (Obaid et al.), not brute-forced.
    (h) **Ringel braid-paper venue** RESOLVED (item 5): **Contemp. Math. 171 (1994),
    339–352** (`RingelBraid1994`); the metaplan's "CMS Conf. Proc. 14" is Crawley-Boevey's
    venue. Transitivity is anchored on both.
    (i) **`is_exceptional_module` implements the BRICK criterion** (`end_dim=1 ∧ Ext¹=0`) — M3:
    `= exceptional` over algebraically closed `k` and on the Dynkin/QQ battery (every Dynkin
    indec is a brick), a sufficient test over non-alg-closed `k` (a non-brick exceptional
    module with a larger division-ring `End` returns False — a scope-limited false negative,
    never a false positive).
    (j) **P45 `mutate` D₄-star defect** (M1, Task 0): a genuine pre-existing defect (spurious
    `status="error"` with the correct 50 vertices on the mixed star) is fixed here with a
    regression pin; independently, the τ-side gate refuses any non-`"complete"` graph loudly.
  - **Recount the class table** (`tests/release/test_oracle_classes.py` drives the numbers —
    run collection, paste the LIVE counts, re-run to green; no guessed at-authoring number).
- [ ] **Step 3: README.** One features line: "exceptional sequences — the classical
  hereditary theory (orthogonality recognizer, braid mutation `σ_i`, braid-orbit
  transitivity, Dynkin closed-form counts `A_n = (n+1)^{n-1}`, `D_4 = 162`) and Buan–Marsh
  τ-exceptional sequences (Jasso τ-perpendicular reduction, the ordered-sτ-tilt bijection
  `#signed = n!·#sτt`) — R27+R28."
- [ ] **Step 4: Full gate:** `... -m pytest tests/modules -q -m deep`, `... tests/webapp
  tests/gui tests/hpc tests/invariants -q -m fast`, `... tests/qpa -q -m qpa`, `...
  tests/release tests/citations -q` — all green.
- [ ] **Step 5:** commit
  `docs(verification): P65 R27+R28 oracle rows + honest scope (rank-2 transitivity, hereditary-only cross-check, orientation-INDEPENDENT #sTt, brick criterion) + citations (Ringel Contemp. Math. 171) + recounted classes`.

---

## Acceptance (Plan-65 definition of done)

1. **Part (b) classical, `quiverlab.modules.exceptional`:** `is_exceptional_module`,
   `is_exceptional_sequence`, `braid_mutation`, `exceptional_sequences`, `c_matrix` public;
   hereditary-only (loud non-hereditary refusal), rep-finite enumeration (loud rep-infinite
   refusal). **HARD literature pins green:** `A_2=3`, `A_3=16`, `A_4=125`, `A_5=1296`,
   `D_4=162` (`count == closed_form_count == n!·h^n/|W|`); `#exceptional == #indec` on Dynkin.
2. **Braid action self-certifies:** every `braid_mutation` output is again a complete
   exceptional sequence (orthogonality re-checked) and `σ_i∘σ_i^{-1} = id` (termwise
   `is_isomorphic`); the braid-orbit BFS from one CES reaches **all** `count` sequences
   (`transitive is True`) — the Crawley-Boevey/Ringel transitivity certificate.
3. **Part (a) τ-exceptional, `quiverlab.tautilting.exceptional`:** `tau_exceptional_objects`,
   `tau_perpendicular_reduction`, `is_tau_exceptional_sequence`, `tau_exceptional_sequences`
   public; τ-tilting-finite gate (loud refusal otherwise; a non-`"complete"` exchange graph —
   `"budget"`/`"error"` — refuses loudly, never a count, M1). **The count `signed_count ==
   n!·#sτt`** with `#sτt` from the P45 `exchange_graph` (`A_2→10`, `A_3→84`) is the formula
   (SELF-CERT); the **cross-engine oracle is MATERIALISATION** (`len(sequences) ==
   signed_count`, distinct, each recognized) on `A₂`/`A₃`/the non-hereditary `k(1↔2)/rad²`
   (H1). `tau_perpendicular_reduction` self-certifies `rk C(U) == n − |U|` **and** the iso pin
   `C(S₁) ≅ kA₂` for `kA₃` (H5); a shifted `P_v[1]` reduces via `A/⟨e_v⟩` (H2).
4. **The (a)↔(b) cross-check holds:** on hereditary rep-finite `A`, the unsigned (all-module)
   τ-exceptional sequences equal the classical exceptional sequences as a set (termwise
   `is_isomorphic`) — a loud bug if they differ. Hereditary-only by design (honest-scope).
5. **Honest scope on the verification page (binding):** rank-2-only mutation transitivity
   (enumeration via the bijection, not BFS); the `signed_count == n!·#sτt` identity is
   self-cert, the materialised count is the cross-check; the `n!` multiplier is the
   Buan–Marsh theorem; the cross-check is hereditary-only; Dynkin-only / τ-tilting-finite-only
   enumeration; char-0/char>dim; **orientation-INDEPENDENT `#sτt`** (the earlier
   "orientation-dependent" claim was false); the BRICK-criterion scope of
   `is_exceptional_module` (M3); E-type counts literature-pinned.
5a. **Task 0 (M1) landed:** P45 `mutate`/`exchange_graph` no longer returns a spurious
   `status="error"` on the τ-tilting-finite D₄ mixed star (regression pins `#sτt(D₄)=50` all
   orientations), and the τ-side gate independently refuses any non-`"complete"` graph.
6. **One algebra-only kind `exceptional_sequences`** clickable end-to-end (GUI canvas →
   block → report) in **all four locales (en/es/fr/zh)** with `pick.kind.*` keys and exact
   key parity, schema v1, both runners byte-identical, the grammar parsed in **all three
   sites** (`hpc/spec.py`, `docs/gui/runner.py`, `webapp/server/schema.py`), ONE golden
   added with a documented change-log entry, canonical keys byte-stable (no `module` block);
   the count tables + transitivity flag + sequence sample render.
7. **QPA battery green (`-m qpa`):** the `NamesGVars()` guard (FAILS if QPA ships an
   exceptional-sequence surface) + the brick/no-self-ext input anchor on the Dynkin zoo.
8. `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the ten
   honest-scope entries (a)–(j); citations added and BibTeX-verified — `BuanMarsh2021` note
   `arXiv:1802.01169`, `RingelBraid1994` = Contemp. Math. 171 (1994) 339–352 (venue RESOLVED),
   the distinct `arXiv:2211.10428` noted separately; README line; deep + fast + qpa + release +
   citations green. No dependency on other Wave-3 plans (independent merge to `dev`).

---

## Methodology & assumptions

**Approach.** I (1) read the metaplan P65 card + the R27/R28 records verbatim; (2)
inventoried the live machinery by import/grep — confirming the full `tautilting` engine
(`exchange_graph`, `SupportTauTiltingPair`, `mutate`, `g_columns`, `bricks`, `is_tau_rigid`),
the module surfaces (`hom_dim`/`ext`/`end_dim`/`is_isomorphic`/`identify_standard`,
`ModuleHom.{kernel,image,cokernel}`, `baer_extension`, `end_algebra`, `presented_form`,
`quotient_by_idempotent`, `Module.tau`), and that **no** Jasso/τ-perpendicular/wide-
subcategory/exceptional machinery exists (so those are genuinely new); (3) **live-verified
the numeric pins** in the venv; (4) mirrored the house form from plan-59/plan-60 (records
verbatim, reference re-verification with `# PIN`, live-facts table, scope gates, Plan-32
markers, three-tier GUI, verification rows, TDD task list).

**Live-verified (venv, QQ) — observed values.** Classical complete exceptional sequence
counts by direct orthogonality enumeration over the knitted universe: **`kA_2 = 3`,
`kA_3 = 16`, `kA_4 = 125`, `kA_5 = 1296`, `kD_4` (subspace) `= 162`** — all matching
`n!·h^n/|W|` (and `A_n = (n+1)^{n-1}`); `#indec = #exceptional` on every Dynkin case tested.
Support-τ-tilting counts from `exchange_graph`: **`A_2 = 5`, `A_3 = 14`, `A_4 = 42`
(Catalan `C_{n+1}`), `D_4` (subspace) `= 50`**; hence the Buan–Marsh signed count
`n!·#sτt = 10 / 84 / 1008 / 1200`. I confirmed the i18n locale set is exactly
`en/es/fr/zh`, that P45 engine tests live in `tests/modules/` (deep bucket), and that the
compute grammar is parsed in **three** sites (`hpc/spec.py`, `docs/gui/runner.py`,
`webapp/server/schema.py`).

**Could NOT verify (labelled to-be-verified-at-implementation).** (1) The **actual
τ-exceptional sequences** from the Jasso recursion — implementing the reduction IS the
plan's work; the `n!·#sτt` count is the standing oracle, with the `#sτt` factor live-verified
and the `n!` multiplier resting on the Buan–Marsh bijection theorem (corroborated for `A_2`
by Igusa–Todorov's `n!·#clusters = 10`). (2) **Braid-orbit transitivity via BFS** — I
verified total counts by direct enumeration, not by braid moves; the orbit-BFS reaching all
of them is verified at implementation against the pinned totals. (3) **E-type Dynkin counts**
(arithmetic from the formula, universe too large to brute-force). (4) The **Kronecker
refusal** probe timed out at `budget_pairs=200` (too heavy for a concurrent-test-safe probe);
the loud budget-cap refusal is guaranteed by the `exchange_graph`/`ar_quiver` code I read
(`status="budget"` when the cap trips), so I rely on the code contract, not a re-run.

**Deliberately not checked, and why.** I did not run any test suite (other agents run tests
concurrently — I kept probes light and short). I did not read every line of the τ-tilting
engine (I read `pairs.py`, `mutation.py`, `torsion.py`, `__init__.py`, `rigid.py` headers —
enough to fix the g-vector column order and the exchange-graph contract). I did not resolve
the **Ringel braid-paper venue** bibliographically (flagged as the one hard `# BLOCK` — the
metaplan's "CMS Conf. Proc. 14" appears to belong to Crawley-Boevey, not Ringel; the worker
must resolve it or anchor transitivity on the verifiable Crawley-Boevey entry).

**Why I believe the result is correct.** The two HARD oracle families are independently
grounded and live-checked: the classical Dynkin counts match a closed-form literature
formula to five data points, and the τ count factors through P45's already-tested
`exchange_graph`. The design reuses only verified live APIs; every scope boundary the
research record names (rank-2 transitivity, hereditary-only, τ-tilting-finite, char>dim,
algebraically-closed for the brick reading) becomes a loud typed refusal. The one genuine
implementation risk — the Jasso reduction `C(U)` and its object-lift `F` — is isolated in
Task A1 with a self-certifying rank check and an honest-scope boundary (count oracle carries
completeness; the deep object-lift is hereditary-scoped for the set-equality cross-check).

**Fix-round addendum (2026-08-08).** The two authoring items above that were left open are now
closed: (1) the **Ringel braid-paper venue** IS resolved — Contemp. Math. 171 (1994), 339–352
(the metaplan's "CMS Conf. Proc. 14" is Crawley-Boevey's); and I additionally found the earlier
`BuanMarsh2021` arXiv note (`2011.02068`) was a WRONG id (a Coptic-NLP paper) — the correct id
is arXiv:1802.01169, with `arXiv:2211.10428` a distinct 2022 companion. New live-verifications
this round (venv, QQ): the shifted-projective reduction `A/⟨e_v⟩` (rank `n−1`), `C(S₁)≅kA₂` for
`kA₃` (perpendicular category has 3 indecs), the worked `σ₁(P₁,S₁)=(S₂,P₁)` mutation,
`#sτt` orientation-independence (`A₃=14`, `D₄=50` across orientations), the `A₅=1296` count
(~3.9 s), a non-hereditary τ-tilting-finite anchor (`k(1↔2)/rad²`, `#sτt=6`), and the M1 D₄-star
`status="error"`-with-50-vertices defect (root-caused to `mutate` missing two exchanges).
Still to-be-verified-at-implementation: the actual τ-exceptional sequences from the recursion
(the materialisation test is the standing cross-check), braid-orbit transitivity via BFS
(including `A₅`'s runtime), and E-type counts.

## Open design risks (what a critic will likely attack)

1. **The Jasso reduction `C(U)` construction and the equivalence `F`.** The DIJ idempotent-
   quotient recipe is stated but the object-transport `F: mod C(U) → J(U) ⊆ mod A` is the
   hardest, least-verified piece; deep-layer lifts to ambient A-modules are scope-bounded.
   *Mitigation:* the completeness oracle is the count `n!·#sτt` (needs no lift); the
   set-equality cross-check is hereditary-only where `F` is the cheap classical perpendicular.
2. **The Buan–Marsh recursion ordering convention** (which end is reduced) — `# PIN`'d to be
   transcribed from the PDF; a wrong convention silently mis-builds sequences (but not the
   count). *Mitigation:* the recognizer and the (a)↔(b) cross-check catch a wrong convention.
3. **`n!·#sτt` as the count** — rests on the bijection theorem; if the implementation's
   "signed"/"ordered" notion drifts from Buan–Marsh, the count mismatches. *Mitigation:* the
   `A_2 = 10` corroboration (`n!·#clusters`) + the hereditary unsigned == classical set-check.
4. **Braid mutation `L_X Y` construction** (four cases, universal extension via
   `baer_extension`) is fiddly; a wrong case gives a non-CES. *Mitigation:* the output is
   re-certified as a CES + the finite-universe `is_isomorphic` cross-identification + the
   involution self-cert.
5. ~~**Ringel braid-paper venue** — the one hard bibliographic block (item 5).~~ **RESOLVED
   (fix round):** Contemp. Math. 171 (1994), 339–352 (`RingelBraid1994`); "CMS Conf. Proc. 14"
   is Crawley-Boevey's.
6. **Enumeration cost / budget** — `D_4` signed = 1200, `E_6` classical = 41472; the block
   ships counts (cheap via formula / exchange graph) not sequences, and the enumeration is
   budget-capped with honest truncation. A critic may push on the estimator sizing for the
   worst case — sized like `ar_quiver`/`tau_tilting`, honest `status` on the cap.

## Change log

- **2026-08-07 authoring.** Initial plan (R27 τ-exceptional via Jasso reduction + the
  ordered-sτt bijection; R28 classical hereditary braid theory). References re-verified;
  classical Dynkin counts (`A_{2..5}`, `D_4`) and `#sτt` (`A_{2..4}`, `D_4`) live-verified in
  the venv; the `n!·#sτt` τ-count formula grounded (factor live-verified, multiplier =
  Buan–Marsh theorem). One hard bibliographic block flagged (Ringel braid-paper venue).
- **2026-08-08 fix round (critic NEEDS WORK → all findings adjudicated valid).** Applied,
  document-only, all live-verified in the venv (QQ) where cheap:
  - **H1** — the τ `signed_count == n!·#sτt` identity relabelled `oracle_selfcert`
    (tautological); added the REAL `oracle_crossengine` MATERIALISATION test (build the
    sequences, `len == signed_count`, distinct, each `is_tau_exceptional_sequence`) on
    `A₂`/`A₃` and the non-hereditary τ-tilting-finite `k(1↔2)/rad²` (`#sτt=6` → 12,
    live-verified). Verification row + honest-scope (b) updated.
  - **H2** — specified the shifted-projective reduction `C(P_v[1]) = A.quotient_by_idempotent
    ([v]) = A/⟨e_v⟩` (real callable, rank `n−1` live-verified); `tau_perpendicular_reduction`
    dispatches on sign; added the full by-hand `kA₂` signed list (10 = 3 all-module=classical
    + 7 with a shifted projective), pinned `(S₁, P₂[1])` + a recognizer test.
  - **H3** — corrected note (f): `#sτt` is **orientation-INDEPENDENT** (live-verified
    `#sτt(A₃)=14` for 3 orientations, `#sτt(D₄)=50` for 3 orientations); the earlier claim was
    false.
  - **H4** — added a fully hand-worked braid mutation `σ₁(P₁,S₁)=(S₂,P₁)` via
    `0→S₂→P₁→S₁→0` (case (b) kernel), identified `L=S₂=(0,1)=('simple',2)`, all live-verified.
  - **H5** — strengthened the reduction self-cert with `C(S₁) ≅ kA₂` for `kA₃` (live-verified:
    `S₁^⊥` has 3 indecs ⟹ connected rank-2), and stated rk alone is insensitive to a wrong
    completion.
  - **M1** — reproduced the D₄ mixed star `{1→2,3→1,4→1}` → `status="error"` WITH the correct
    50 vertices (real P45 `mutate` defect: fails at summands 0,3, "no valid exchange found",
    over QQ); added **Task 0** (P45 fix + regression) and a τ-side `status != "complete"` loud
    gate; scope table + honest-scope (j).
  - **M2** — replaced "the general canonical two-step" (case d) with the concrete ordered
    recipe (universal extension FIRST, then evaluation-kernel; `[L_XY]=⟨X,Y⟩[X]−[Y]`), cited
    Ringel 1994 / Crawley-Boevey 1993, `# PIN`'d the module-representative statement.
  - **M3** — scoped `is_exceptional_module` honestly as the BRICK criterion (`= exceptional`
    over alg-closed `k` and on Dynkin/QQ; sufficient otherwise); scope table + honest-scope
    (i).
  - **W2** — added `A₅=1296` to the Task B3 parametrize (count leg re-verified ~3.9 s);
    `transitive` for `A₅` gated to an opt-in slow test.
  - **W3/W4** — resolved the bibliography: `BuanMarsh2021` = arXiv:1802.01169 (the earlier
    `2011.02068` was a wrong id — a Coptic-NLP paper); `arXiv:2211.10428` = the DISTINCT 2022
    "Mutating signed τ-exceptional sequences" (noted separately); `RingelBraid1994` =
    Contemp. Math. 171 (1994), 339–352 (`# BLOCK` dropped); `BuanHansonMarsh2024` title fixed
    to "…pairs and sequences"; added the `ringel_braid` registry entry.
  - **Consistency fix (beyond the rulings, same file):** Task B1's `A₂` orthogonality test had
    the exceptional pair reversed — it is `(S₂, P₁)` that is exceptional (`Hom(P₁,S₂)=0`), not
    `(P₁, S₂)` (`Hom(S₂,P₁)=1` is a backward map); corrected to match the doc's own convention
    and the verified Hom/Ext (live-verified this fix round). Flagged so implementers know.
