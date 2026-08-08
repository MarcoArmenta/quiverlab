# Plan 73: Han's-conjecture bounded-extension recognizer + Jacobi–Zariski nearly-exact sequence (P73 / R7)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking. **Do the reference re-verification (below) BEFORE
> writing any oracle pin** — the flagged attributions (`# PIN`, and the corrected Jacobi–Zariski
> citation) must be resolved against the sources, not this plan's prose. The bounded-extension
> DECISION is an **honest capped semi-decision**: caps are stated, `"undecided"` is a first-class
> verdict, refusals are loud, no silent truncation — this is the plan's core value.

**Goal.** Decide, and certify, when an inclusion of finite-dimensional algebras `B ⊆ A` is a
**bounded extension** (Cibils–Lanzilotta–Marcos–Solotar), then **transport Han's conjecture**
across it — `B ⊨ Han ⇔ A ⊨ Han` — using the **Jacobi–Zariski long nearly exact sequence**
("exact twice in three") as the computational tool and the self-certification gate. Three
deliverables, each grounded on machinery quiverlab already ships plus one focused new substrate:

- **(1) The bounded-extension recognizer.** `B ⊂ A` is **left (resp. right) bounded** iff the
  `B`-bimodule `A/B` is (i) **`B`-tensor nilpotent**, (ii) of **finite projective dimension over
  `B^e`**, and (iii) **left (resp. right) `B`-projective** (CLMS 2101.02597, Introduction /
  Def. 2.3). The recognizer decides each leg **honestly**:
  - **(i) tensor-nilpotency** — the load-bearing leg — as a **capped semi-decision** on the direct
    tensor powers `(A/B)^{⊗_B m}`, PLUS a **certificate route** (relative cycles + `J`-interrupters,
    CLMS Thm 5.14/5.16 ⟹ certified `True`; a non-`J`-interrupted relative cycle whose tensor
    powers never die, CLMS Ex. 5.5 ⟹ certified `False`). Cap reached with no certificate ⟹
    `"undecided"` (never a fake verdict).
  - **(iii) one-sided projectivity** — `A/B` as a one-sided `B`-module, `pd_B = 0` on the shipped
    exact-Domain module resolution stack (decidable), with the CLMS Thm 5.20 combinatorial
    sufficient condition as a fast certificate.
  - **(ii) finite `pd_{B^e}(A/B)`** — **PRIMARY route: `gl.dim B < ∞`** (shipped
    `global_dimension`). By CLMS Ex. 6.1, `gl.dim B < ∞ ⟹ gl.dim B^e < ∞ ⟹ pd_{B^e}(A/B) < ∞` — the
    leg is **FREE**, no enveloping algebra needed (live-verified `gl.dim B = 2` for Ex. 5.3 AND
    Ex. 5.5). **FALLBACK (`gl.dim B = ∞`):** a capped semi-decision on the new enveloping algebra
    `B^e = B ⊗ B^op` (presented as a bound quiver algebra) with `A/B` as a right `B^e`-module, `pd`
    via the shipped `Module.projective_resolution(length).pd()` (finite within `length` ⟹ certified;
    else `"undecided"` / certified lower bound).
- **(2) Han transport — the injection/isomorphism ladder (CLMS Thm 3.1/4.6, read exactly).** The
  legs give a **graded** conclusion, NOT a blanket agreement — each row is honestly labelled with its
  exact hypotheses and coefficients:
  - **leg (i) alone (`A/B` tensor-nilpotent):** `H_*(B, A) ↪ H_*(A, A)` for `* ≫ 0` — an
    **injection** in homology **with coefficients in `A`** (mixed coefficients), NOT ordinary HH,
    NOT an isomorphism.
  - **legs (i)+(ii) (also finite `pd_{B^e}(A/B)`):** `HH_*(B) ↪ HH_*(A)` for `* ≫ 0` — an
    **injection**, ordinary coefficients (hence `HH_*(A)=0 ⟹ HH_*(B)=0`).
  - **all three legs = bounded (a fixed side):** `HH_*(B) ≅ HH_*(A)` for `* ≫ 0` — the
    **isomorphism**; only THIS delivers the full Han transport `B ⊨ Han ⇔ A ⊨ Han` (Thm 4.6).
  The `han_transport` certificate labels its claim by the exact row reached — `"bounded"` (iso),
  `"pd_injection"` (i+ii, ordinary), `"nilpotent_injection"` (i, mixed coeff), `"not_bounded"`
  (leg i `False`), `"undecided"` (a cap hit). **There is NO "≅ from leg (i) alone" verdict** — that
  is mathematically unavailable (the earlier "homology_agreement" verdict misread Thm 3.1; removed).
- **(3) The Jacobi–Zariski nearly-exact sequence as tool + self-cert gate.** Compute the relative
  Hochschild homology `HH_*(A|B, X)` via the **CLMS normalized relative bar complex** (Thm 2.2),
  **finite** once `A/B` is tensor-nilpotent (Cor 2.4). Assemble the JZ sequence
  `⋯ → HH_*(B) → HH_*(A) → HH_*(A|B) → HH_{*−1}(B) → ⋯` (CLMS 2009.05017 / the classical Kaygun
  origin, "exact twice in three") from `HH_*(A)` / `HH_*(B)` (shipped engines) and `HH_*(A|B)` (new).
  **Self-cert gate (the honest invariant):** `HH_*(A|B) = 0` for `* ≥ n` (Cor 2.4); exactness at the
  two guaranteed spots; and the **injection bound `dim HH_m(B) ≤ dim HH_m(A)`** for `* ≫ 0` — with
  **equality asserted ONLY under the full bounded certificate**. (The earlier "`HH_*(A) ≅ HH_*(B)`
  when tensor-nilpotent" gate was a false invariant — nilpotency gives an injection, not an iso —
  and is removed.)

**GUI.** ONE new **algebra-level certificate** compute kind **`han_transport`** carrying the
**new-arrow subset `F`** (the arrows whose removal from `A` defines `B`) as a request-level field,
served by **all three tiers** (schema guard + both runners byte-identical), the recognizer-panel
row + an **arrow-subset picker** on the drawn quiver, **all four locales (en/es/fr/zh)** with exact
key parity, canonical keys stable, report rendering, and citations.

**Architecture.** Three new source modules + one new enveloping constructor + one algebra-level
kind; everything else is a thin exact layer over shipped primitives:

- **`src/quiverlab/families/extension.py`** (new) — the **input model** ("extension by arrows and
  relations", CLMS Def. 5.2) and the relative-path combinatorics:
  - `arrow_removal_subalgebra(A, new_arrows) -> Extension` — build `B ⊆ A` by designating a subset
    `F ⊆ (Q_A)_1` of arrows as "new"; `B = kQ_B/I_B` with `Q_B = (Q_0, (Q_A)_1 ∖ F)` and
    `I_B = ker(π: kQ_B → A)` extracted + **dimension-certified** by the shipped
    `families/_present.py::present_from_pi` (the honest subalgebra — NOT "drop the `F`-relations",
    which can miss induced relations; see the Mathematical foundation). Validates `J ∩ B = 0`
    (the extension axiom) via the `dim B` certificate.
  - `Extension` dataclass: `A`, `B`, `new_arrows`, `relative_paths()` (the `A/B` basis = `A`-basis
    paths using `≥ 1` new arrow), `relative_cycles(cap)` (CLMS Def. 5.11), and the `B`-bimodule
    action data of `A/B`.
  - `enveloping_algebra(B) -> Algebra` — `B^e = B ⊗ B^op` presented as a **bound quiver algebra**
    (product quiver `Q_B ⊗ Q_B^op`, relations extracted by `present_from_pi` against the shipped
    `TensorProduct(B, B.opposite())` structure constants; certified `dim B^e = (dim B)^2`). The
    first first-class enveloping algebra in quiverlab (the Hochschild engines only handle `A^e`
    internally). `Algebra.enveloping()` delegate.
- **`src/quiverlab/invariants/han.py`** (new) — the recognizer, transport, and certificate:
  - `is_tensor_nilpotent(ext, *, cap=8) -> TensorNilpotency` — the capped semi-decision + the
    `J`-interrupter / no-relative-cycle certificate + the Ex-5.5 non-vanishing witness.
  - `one_sided_projective(ext, *, side="left") -> OneSidedProjectivity` — `pd_B(A/B) = 0` on the
    module stack + the Thm 5.20 combinatorial certificate.
  - `bounded_extension(A, new_arrows, *, side="auto", nilp_cap=8, pd_cap=16) -> BoundedCertificate`
    — the three-leg certificate.
  - `han_transport(A, new_arrows, *, ...) -> HanTransport` — the certificate kind + the self-cert
    gate; scales the transport claim to what is certified.
  - `han_transport_block(A, spec, *, budget=...) -> dict` — the JSON block for the runners.
- **`src/quiverlab/hochschild/jacobi_zariski.py`** (new) — the JZ computational tool:
  - `relative_homology(ext, top, *, coefficients=None) -> HHTable` — `HH_*(A|B, X)` via the CLMS
    normalized relative bar complex (Thm 2.2) over a **general** `B` (the shipped
    `hochschild/relative.py` does only `B = kQ_0`); finite when `A/B` is tensor-nilpotent.
  - `jacobi_zariski_sequence(ext, top) -> JZSequence` — the assembled nearly-exact sequence with
    per-degree exactness flags and the gap (`Ker/Im`) at the middle spots.
- The `han_transport` kind is wired into `hpc/spec.py::_dispatch` + `docs/gui/runner.py`
  (byte-identical twins) + `webapp/server/schema.py` (the request-field guard + grammar), the GUI
  touchpoints (checkbox / picker / `renderBlock` / `scheduleProbe`), the i18n chains, and
  `trace/results_html.py`, following the `recognizers` / `derived_compare` (second-input-carrying)
  precedents.

**Tech Stack.** Pure exact `Domain` linear algebra + exact combinatorics on the finite relative-path
universe. **No floats in `src/`** (AST-gated by `tests/test_no_floats.py`): dims/indices/verdicts
are `int`/`bool`/`str`, HH tables are integer dim lists, actions are exact-Domain matrices. HH via
the shipped `Algebra.hochschild_homology`/`hochschild_cohomology` (`engine="cs"` over any Domain,
`"auto"`/`"fast"` over `GF(p)`); module `pd` via `Module.projective_resolution`; subalgebra /
enveloping presentation via `families/_present.py::present_from_pi`.

---

## Records (verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`)

> **R7 — Han's-conjecture bounded-extension recognizer + Jacobi–Zariski sequence.** [B-scout P3;
> keep-with-corrections] Object: decide bounded extension B ⊆ A (tensor-nilpotency of A/B over B —
> an HONEST CAPPED SEMI-DECISION, powers can grow before dying, no a-priori index bound; finite pd
> of A/B over B^e; one-sided projectivity), then transport Han B ⊨ Han ⇔ A ⊨ Han and expose the
> Jacobi–Zariski nearly-exact sequence ("exact twice in three", per arXiv:1908.11130). Refs: CLMS
> arXiv:2101.02597 (J. Algebra 2022; J-interrupter Ex. 5.3/5.5 as oracles incl. the negative case);
> arXiv:1908.11130; arXiv:2409.00945 (recollement approach — AUTHOR ATTRIBUTION UNVERIFIED, re-check
> at plan time). Size M–L.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P73 (Wave 4, tier δ, **independent**).
Size M–L.

**Seam with P72 (R5+R6; sibling, `[needs P52]`).** P72 owns the **HH computation accelerators**:
(a) the split-extension long exact sequence `HH^*(B⋉M)` from `HH^*(B)+H^*(B,M)` (CMRS
math/0102194), and (b) certified **arrow removal/addition HH reductions** (CLMS 1812.07655). P73
owns the **recognizer + transport + JZ**: it decides whether an extension is *bounded* (possibly
non-split), transports Han's conjecture, and computes the JZ sequence. The **shared substrate** is
the "extension by arrows and relations" `B ⊆ A` model — the general **arrow-subset subalgebra
constructor** (`arrow_removal_subalgebra`) and the enveloping algebra `B^e`. **Shared-substrate
contract:** the arrow-removal constructor is owned by **whichever of P72/P73 lands on `dev` first**;
the later plan **consumes** it (does not re-implement). P73 specifies it in Task 0 with a stable
signature (`arrow_removal_subalgebra(A, new_arrows) -> Extension`) so P72 can depend on it, and vice
versa if P72 merges first — the second plan's Task 0 becomes a no-op that imports and re-uses.
(The P72 doc's critic is checking this collision; the contract is stated on both sides.) P73 does
NOT depend on P72 for merge (independent). Neither plan re-implements the other's HH formulas.

---

## Reference re-verification (done at authoring; findings binding, `# PIN` = resolve at citation-add)

The metaplan standing rule (§1.5) blocks P73 on the flagged `arXiv:2409.00945` attribution and
requires re-verifying every citation. **All four load-bearing arXiv abstracts were fetched at
authoring; two card/record citation errors were found and corrected.**

1. **`arXiv:2409.00945` — VERIFIED (the mandatory re-check).** Title **"A recollement approach to
   Han's conjecture"**; authors **Ren Wang, Xiaoxiao Xu, Jinbi Zhang, Guodong Zhou**; submitted
   **2 Sep 2024** (fetched from the arXiv abstract page). Abstract confirms it is on-topic: a
   recollement/derived-category **reduction** of Han's conjecture (Han holds for the middle ring iff
   for the two side rings), applied to Morita contexts, exact contexts, **skew-gentle algebras**,
   finite EI-category algebras, and GLS algebras. The record's "AUTHOR ATTRIBUTION UNVERIFIED" flag
   is **cleared** — the paper is real, correctly identified, and citable. It is a **corroborating /
   context** reference (an independent proof route to the same transport phenomenon; the skew-gentle
   result ties to P68). Ship it as `WangXuZhangZhou2024` (`@misc`, `arXiv:2409.00945`; upgrade to
   `@article` only if a journal ref exists at merge — no fabricated fields).

2. **The Jacobi–Zariski arXiv id is WRONG in the record — CORRECTED (with an integrity note on the
   card).** **What the metaplan card actually says** (§5 P73, verbatim): "the JZ nearly-exact
   sequence ('exact twice in three', 1908.11130)". The card does **NOT** write "Kaygun" — the name
   "Kaygun" appears nowhere in the card or in R7. (An earlier draft of this plan wrongly attributed
   "Kaygun 1908.11130" to the card; that was this plan-writer's error, corrected here. The Kaygun
   association is real but comes from **CLMS's own credit** — see below — not from the card.) The
   record's R7 lists `arXiv:1908.11130` as the JZ source. **The arXiv id is wrong, verified by
   fetching the abstracts:**
   - **`arXiv:1908.11130`** is **"Split bounded extension algebras and Han's conjecture"** by
     **Cibils, Lanzilotta, Marcos, Solotar**, **Pacific J. Math. 307 (2020) 63–77**, DOI
     `10.2140/pjm.2020.307.63`. It is the **split-case predecessor** of 2101.02597 — NOT the
     Jacobi–Zariski paper.
   - **The Jacobi–Zariski "long nearly exact / exact twice in three" sequence is
     `arXiv:2009.05017`**, **"Jacobi–Zariski long nearly exact sequences for associative algebras"**
     by **Cibils, Lanzilotta, Marcos, Solotar**, **Bull. Lond. Math. Soc. 54 (2022), no. 3**, DOI
     `10.1112/blms.12516` (year confirmed 2022, not 2021). Its abstract states verbatim: *"we
     obtain a Jacobi–Zariski long nearly exact sequence relating the Hochschild homologies of A and
     B, and the relative Hochschild homology … This long sequence is exact twice in three."* This is
     the operative source 2101.02597 relies on (cited there as "[12]").
   - **Kaygun IS credited — and IS cited (verified).** CLMS 2009.05017 §2 credit **A. Kaygun** as
     the first to obtain the noncommutative Jacobi–Zariski sequence in Hochschild/cyclic homology.
     That paper is **"Jacobi–Zariski Exact Sequence for Hochschild Homology and Cyclic (Co)Homology",
     Homology Homotopy Appl. 14 (2012), no. 1, 65–78, arXiv:1103.4377** (verified this fix round; a
     later Erratum exists). It is the **classical origin** (`B ⊆ A` with `A/B` flat); CLMS 2009.05017
     is the **"long nearly exact / exact twice in three"** refinement the plan actually computes with.
     **Ship both**, with the lineage in the annotation — this is the scholarly-correct resolution and
     it justifies the Kaygun name the card's source (CLMS) invokes.
   - **Correction to ship:** cite the JZ sequence as **`arXiv:2009.05017` (CLMS, Bull. LMS 54 (2022))**
     + **Kaygun `arXiv:1103.4377` (the classical origin)**; cite `arXiv:1908.11130` as the split-case
     predecessor (context). **Fold the arXiv-id correction back into research-doc R7 with a dated
     note** (metaplan §2; not done in this plan — flagged for the merge step).

3. **`arXiv:2101.02597` — VERIFIED (the primary oracle source).** Title **"Han's conjecture for
   bounded extensions"**; authors **Cibils, Lanzilotta, Marcos, Solotar**; v3 **4 Feb 2022**;
   **J. Algebra (2022)**, DOI **`10.1016/j.jalgebra.2022.01.022`** (exposed on the arXiv abstract
   page; verified this fix round). The bounded-extension definition (Introduction), the
   tensor-nilpotency machinery (§2, §5), the main transport theorem (Thm 4.6), and **Examples 5.3
   (bounded) / 5.5 (not bounded)** are transcribed in the Mathematical foundation from the fetched
   PDF. `# PIN`: exact **volume / page range** confirmed against the journal landing page at
   citation-add (the arXiv id + DOI are solid; do not guess the volume/pages — delete any field that
   cannot be confirmed).

4. **Han's conjecture — source VERIFIED.** **Y. Han, "Hochschild (co)homology dimension", J. London
   Math. Soc. (2) 73 (2006), no. 3, 657–668**, DOI `10.1112/S002461070602299X`. This is the paper
   stating the conjecture the whole plan transports (`HH_*(A,A) = 0` for `* ≫ 0` ⟹ `gl.dim A < ∞`).
   Ship as `Han2006`. (The scout's "Bull. LMS 38" was wrong — it is **J. LMS 73**; verified.)

5. **`arXiv:2301.07511` — VERIFIED (optional context).** "A survey on Han's conjecture" (CLMS,
   Latin American J. Math.) — a clean landscape reference; add as `@misc` context only if an
   annotation cites it (not load-bearing).

**Already shipped and reused (verify live in `registry.py` at implementation):** `assem_book`
(ASS2006) **and `han_conjecture` (Han2006) — already registered at `registry.py:210`, reuse it;
the P72 critic caught this doc listing it as new**. **New keys to add:** `clms_bounded_extensions`
(CLMS2101.02597), `clms_jacobi_zariski` (CLMS2009.05017 — the corrected JZ source),
`kaygun_jacobi_zariski` (Kaygun2012, arXiv:1103.4377 — the classical noncommutative-JZ origin CLMS
credit), `clms_split_bounded` (CLMS1908.11130 — split predecessor), `wang_recollement_han`
(2409.00945 — verified). Cross-check the P72 arrow-removal key `clms_arrow_removal` (1812.07655) if
P72 has merged; do not duplicate.

---

## Mathematical foundation (the definitions, exactly — the plan's ground truth)

Throughout `A = kQ_A/I_A` is basic, connected, finite-dimensional over an exact `Domain`; paths
compose **left-to-right** (quiverlab convention; the CLMS paper composes right-to-left, so every
transcribed path below is translated — see the live-verified examples). `HH_*` / `HH^*` are
Hochschild homology / cohomology; `B^e = B ⊗_k B^op`; a `B`-bimodule = a right `B^e`-module.

### The input model — "extension by arrows and relations" (CLMS Def. 5.1/5.2)

CLMS Thm 5.8 says: given `B = kQ/I` **first**, every finitely generated `B`-algebra `A` arises by
adding arrows + relations to `B`'s quiver. **The plan's tool goes the other direction** — the user
gives `A = kQ_A/I_A` and picks a subset `F ⊆ (Q_A)_1` of arrows to remove — so it can only express
extensions whose `B` sits on a **sub-quiver of `Q_A`** (`B` generated by kept arrows). **Honest
scope (M1):** the tool is general **among extensions by arrows and relations already presented on
`Q_A`**; a bounded extension `B ⊂ A` whose `B` is NOT generated by a subset of the user's arrows of
`A` is **inexpressible** in this input model (it would require re-presenting `A` on a quiver that
exhibits `B`'s arrows). This narrowing is disclosed as a loud scope-gate row (below), not silently
assumed to be "fully general". Concretely:

- Start from `A = kQ_A/I_A` (a shipped quiverlab `Algebra`). Designate a subset
  `F ⊆ (Q_A)_1` of arrows as **new** (`F` = the CLMS "new arrows"). Then `B` is the subalgebra of
  `A` generated by the **kept** arrows `(Q_A)_1 ∖ F` (and the vertices): `B = kQ_B / I_B`,
  `Q_B = (Q_0, (Q_A)_1 ∖ F)`, `I_B = ker(π: kQ_B → A)`. The extension axiom is `J ∩ B = 0`
  (CLMS Def. 5.2), automatic here because `B` is *defined* as `π`'s image — the certificate is
  `dim B` (the map is injective iff `dim(kQ_B/I_B) = dim(span of F-free A-basis paths)`).
- **`A/B` (the `B`-bimodule)** is spanned by the **relative paths** of positive `F`-length: the
  `A`-basis paths that use `≥ 1` new arrow (CLMS Def. 5.9–5.10; live-verified — the `A/B` basis is
  read directly off `A.basis_labels` by filtering on new-arrow usage).
- **Honest constructor (M-critical).** `I_B` is `ker(π)`, NOT `{relations of A avoiding F}`: an
  `F`-relation can induce a collapse among `F`-free paths. The constructor drives
  `present_from_pi(Q_B, img={kept arrow ↦ its coordinate in A}, T=A's structure constants, dim_
  expected=?)`. Because `dim B` is unknown a priori, use `present_from_pi` in its
  **dimension-discovering** mode (extract `I_B = ker` by length-lex, then set the algebra), and
  **certify** `dim B = #(F-free A-basis paths that are linearly independent in A)` (live-verified
  `= 20` for Ex. 5.3, `= 6` for Ex. 5.5). If the shipped `present_from_pi` requires a known
  `dim_expected`, Task 0 adds the two-pass (count `ker` dimension, then present).

### The bounded-extension definition (CLMS 2101.02597, Introduction; Def. 2.3)

> An extension `B ⊂ A` of `k`-algebras is called **left (resp. right) bounded** if the `B`-bimodule
> `A/B` is `B`-tensor nilpotent, its projective dimension is finite and it is left (resp. right)
> `B`-projective.

- **(i) Tensor nilpotency (Def. 2.3).** A `B`-bimodule `M` is **tensor nilpotent** if there exists
  `n` with `M^{⊗_B n} = 0`. The smallest such `n` is the nilpotency index. *No a-priori bound* —
  powers can grow before dying (the record's warning), hence the **capped semi-decision**.
- **(ii) Finite projective dimension:** `pd_{B^e}(A/B) < ∞`.
- **(iii) One-sided projectivity:** `A/B` is projective as a left (resp. right) `B`-module.

### Deciding tensor nilpotency (CLMS §5 — the certificate route)

- **Relative path** (Def. 5.9): a concatenation `β_n a_n ⋯ β_1 a_1` with `β_i` paths of `Q` not in
  `I` (nonzero paths of `B`) and `a_i ∈ F`. `F`-length = `#` new arrows. Relative paths of positive
  `F`-length span `A/B`.
- **Relative cycle** (Def. 5.11): a relative path `a_1 β_n a_n ⋯ β_1 a_1` that closes up
  (`t(β_n) = s(a_1)`).
- **`J`-interrupter** (Def. 5.13): an arrow `a_i` of the relative cycle such that the "product
  across it", `\overline{β_i a_i β_{i-1}}` (resp. `\overline{β_1 a_1 β_n}` for `i=1`), lies in `B`
  (i.e. the class in `A = T/J` reduces to an `F`-free element).
- **Certified nilpotent (Thm 5.14 + Thm 5.16).** If there are **no relative cycles**, or **every
  relative cycle has a `J`-interrupter**, then `A/B` is `B`-tensor nilpotent. The nilpotency index
  is bounded by the **length index** (Def. 5.19: the max `F`-length of a relative path in a
  `k`-basis of `A`, `+1` — finite since `A` is f.d.).
- **Certified NOT nilpotent (Ex. 5.5 pattern).** A relative cycle whose "self-product" survives in
  every tensor power is a non-nilpotency witness: e.g. a relative path `w` (a cycle at a vertex,
  `≥ 1` new arrow, no `J`-interrupter) with `w ⊗_B w ⊗_B ⋯ ⊗_B w ≠ 0` for all `m`. The recognizer
  produces `w` and certifies `False`.
- **The capped direct route (the semi-decision the card mandates).** Compute `(A/B)^{⊗_B m}` by
  relative-path concatenation for `m = 1, 2, …, cap`: the `m`-th power is spanned by concatenable
  `m`-tuples of relative paths, reduced modulo `J` and modulo the "tensorand lands in `B`" collapse
  (a summand with a factor of `F`-length 0 vanishes). If some `m ≤ cap` gives dimension `0` ⟹
  `True` (index `m`). If the cap is reached ⟹ **`"undecided"`** (unless the certificate route
  already decided). No guessed verdict.

### Deciding one-sided projectivity (CLMS Thm 5.20)

`A/B` is projective as a **left** `B`-module iff `pd_B(A/B \text{ as a left } B\text{-module}) = 0`
(decidable on the shipped module stack). CLMS Thm 5.20 gives a combinatorial **sufficient
certificate**: if the projection `J_{|0,n|}` of `J` on `T^{|0,n|}` is generated as a `B`-bimodule by
relations of the form `r = a·Σ λ_γ γ` (a single new arrow `a ∈ F` times relative paths, `λ_γ ∈ k`),
then `A/B` is left `B`-projective (verified in Ex. 5.3: `J = ⟨abcd − αβ⟩`, `abcd = a·(bcd)`).
Reversing the argument gives the right-projective certificate.

### The Jacobi–Zariski sequence (CLMS 2009.05017; recalled in 2101.02597 §2)

- **Normalized relative bar complex (Thm 2.2).** For an `A`-bimodule `X` and a `k`-section `σ` of
  `π: A → A/B`, `HH_*(A|B, X)` is the homology of
  `⋯ → X ⊗_{B^e} (A/B)^{⊗_B m} → ⋯ → X ⊗_{B^e} A/B → X_B → 0`
  with the differential `b` (Thm 2.2, `σ`-independent). **When `A/B` is tensor nilpotent with
  `(A/B)^{⊗_B n} = 0`, this complex has finitely many nonzero terms** (`m < n`) ⟹ `HH_*(A|B) = 0`
  for `* ≥ n` (Cor. 2.4). This makes `HH_*(A|B)` **exactly computable** by a finite complex on the
  relative-path bases (the same combinatorics as leg (i)).
- **The JZ long nearly exact sequence (Thm 3.1).** For `X` an `A`-bimodule there is a sequence
  `⋯ → HH_m(B, X) → HH_m(A, X) → HH_m(A|B, X) → HH_{m-1}(B, X) → ⋯`
  which is **"exact twice in three"** (Def. 2.5): exact at the `HH_*(B)` and `HH_*(A|B)` spots; its
  **gap** at the `HH_*(A)` spot is `Ker/Im`, a spectral sequence converging to it.
- **The transport consequence — an INJECTION, not an isomorphism, until bounded (H1, read exactly).**
  With `A/B` tensor nilpotent, `HH_*(A|B) = 0` for `* ≥ n`. Exactness AT the `HH_*(B)` spot then
  forces, for `* ≫ 0`:
  - **taking `X = A`** (the mixed-coefficient row Thm 3.1 gives from nilpotency alone):
    `H_*(B, A) ↪ H_*(A, A)` — an **injection**.
  - **adding finite `pd_{B^e}(A/B)`** (leg ii): `HH_*(B) ↪ HH_*(A)` — an injection in **ordinary**
    HH (so `HH_*(A)=0 ⟹ HH_*(B)=0`).
  - **adding one-sided projectivity too** (bounded): `HH_*(B) ≅ HH_*(A)` — the **isomorphism**.
  The gap at the `HH_*(A)` spot (Thm 3.1 is NOT exact there) is exactly why nilpotency alone gives
  only the injection, never surjectivity/iso: the ≅ needs the full bounded certificate. Do NOT claim
  `HH_*(A) ≅ HH_*(B)` from tensor-nilpotency alone.

### Han transport (CLMS Thm 4.6, via Thm 4.2/4.3)

> **Thm 4.6.** Let `B ⊂ A` be a **left or right bounded** extension of finite-dimensional
> `k`-algebras. Then `B` satisfies Han's conjecture **iff** `A` does. (No splitting assumed.)

The proof combines the JZ homology transport (needs only **tensor nilpotency**) with the
global-dimension transport: **Thm 4.3** (`A/B` tensor nilpotent + one-sided projective ⟹
`gl.dim A ≤ n + b − 1` when `gl.dim B = b`) and **Thm 4.2** (`pd_{B^e}(A/B) = r < ∞` + one-sided
projective ⟹ `gl.dim B ≤ r + a` when `gl.dim A = a`). So:
- **Forward** (`B ⊨ Han ⟹ A ⊨ Han`) uses tensor-nilpotency + one-sided projectivity.
- **Backward** (`A ⊨ Han ⟹ B ⊨ Han`) uses tensor-nilpotency + **finite `pd_{B^e}`** + one-sided
  projectivity.

**The recognizer's honest transport verdict — the injection/iso LADDER (H1; each row labelled with
its exact hypotheses and coefficients):**

| legs certified | `HanTransport.transport` | claim (for `* ≫ 0`) |
|---|---|---|
| (i)+(ii)+(iii), a fixed side (= **bounded**) | `"bounded"` | `HH_*(B) ≅ HH_*(A)` (iso) ⟹ full `iff` `B ⊨ Han ⇔ A ⊨ Han` (Thm 4.6) |
| (i)+(ii) | `"pd_injection"` | `HH_*(B) ↪ HH_*(A)` (**injection**, ordinary coeff) ⟹ `HH_*(A)=0 ⟹ HH_*(B)=0` |
| (i) only | `"nilpotent_injection"` | `H_*(B, A) ↪ H_*(A, A)` (**injection**, coefficients in `A`) — NOT ordinary HH, NOT an iso |
| (i) = `False` | `"not_bounded"` | not a bounded extension via this route (Ex. 5.5) |
| any leg `"undecided"` | `"undecided"` | honest cap reached; per-leg status reported |

**There is NO `≅`-from-leg-(i) verdict.** The self-cert gate checks `dim HH_m(B) ≤ dim HH_m(A)` (the
injection bound, computable on the shipped engines) and asserts **equality only under `"bounded"`**.

### The two consequences the plan leans on

1. **The relative-path universe is finite and directly readable.** `A/B` basis = `A`-basis paths
   using `≥ 1` new arrow (live-verified). So relative paths, relative cycles, `J`-interrupters, and
   tensor powers are exact finite combinatorics on `A.basis_labels` + `A.multiply` — no new HH
   engine, no guessing.
2. **The injection bound is the honest headline gate — an inequality, not an equality.** The JZ
   consequence from leg (i)+(ii) is `dim HH_m(B) ≤ dim HH_m(A)` for `* ≫ 0` (an injection), NOT
   `HH_*(A) ≅ HH_*(B)` (a false invariant from nilpotency alone). Equality holds only under the full
   bounded certificate. leg (ii) is **cheap when `gl.dim B < ∞`** (CLMS Ex. 6.1 — live-verified
   `gl.dim B = 2` for both Ex. 5.3 and Ex. 5.5), so the `pd_injection` row (and, with one-sided
   projectivity, `"bounded"`) is reachable without the enveloping algebra on the common
   finite-gl.dim case; the enveloping-algebra fallback is only for `gl.dim B = ∞`, honestly capped.

---

## Live-verified facts (all recomputed in the venv at spec time)

`quiverlab.__file__` confirmed `= …/quiverlab/src/quiverlab/__init__.py` (the dev-tip main-repo src,
which carries the merged P51–P62 including the P52 relative-HH / `Bimodule` coefficient surface;
this plan's worktree is pinned at v0.3.0 and is used only to hold this document — see Methodology).
Both example algebras were **built and computed** (`GF(32003)`, `engine="auto"`, deg 0..9):

| example | `dim A` | `dim B` | `dim A/B` | `gl.dim B` | `HH_*(A)` (0..9) | `HH_*(B)` (0..9) | verdict |
|---|---|---|---|---|---|---|---|
| **Ex. 5.3** (bounded) | **47** | **20** | **27** | **2** | `[5,0,0,0,0,0,0,0,0,0]` | `[5,0,0,0,0,0,0,0,0,0]` | `dim HH_m(B) = dim HH_m(A)` (iso — bounded) |
| **Ex. 5.5** (NOT bounded) | **10** | **6** | **4** | **2** | `[4,1,2,3,3,3,4,6,8,9]` | `[3,0,0,0,0,0,0,0,0,0]` | `HH_*(A) ≠ HH_*(B)` — extension NOT bounded |

- **Ex. 5.3 built** as `A = kQ_A/⟨alpha*beta, d*c*b*a − beta*alpha⟩` on
  `Q_A = {al:5→1, be:1→5, d:1→4, c:4→3, b:3→2, a:2→1}` (`F = {a}`); `B = kQ_B/⟨alpha*beta⟩` on
  `Q_B = Q_A ∖ {a}`. The **`a`-free `A`-basis has exactly 20 elements = `dim B`** (live-verified:
  includes `be*al`, the `1→5→1` cycle that equals `d*c*b*a` in `A`) ⟹ the inclusion `B ↪ A` is
  injective (`J ∩ B = 0`), the arrow-removal subalgebra is genuine. **`gl.dim B = 2` (live-verified,
  exact)** ⟹ by CLMS Ex. 6.1 `pd_{B^e}(A/B) < ∞` FREE (no enveloping algebra needed — B1). `HH_*(A)
  = HH_*(B) = [5,0,0,…]`: **both vanish for `* ≥ 1`**, so the injection `HH_*(B) ↪ HH_*(A)` is here
  an **equality** (consistent with `"bounded"` ⟹ iso; the self-cert gate asserts equality only under
  the full certificate). Both `A` and `B` satisfy Han's conjecture (finite gl.dim). **HARD literature
  pin** (`oracle_literature`): the recognizer must certify Ex. 5.3 **bounded** (tensor-nilpotent
  index 2, left `B`-projective, `gl.dim B = 2` ⟹ finite `pd_{B^e}`), verdict `"bounded"`.
- **Ex. 5.5 built** as `A = kQ_A/⟨a*b, d*a, b*d, d*c*d⟩` on
  `Q_A = {a:1→2, b:2→3, c:1→3, d:3→1}` (`F = {d}`); `B = kQ_B/⟨a*b⟩` on `Q_B = Q_A ∖ {d}`. The
  **`d`-free basis has exactly 6 elements = `dim B`**; the **`A/B` basis is `{d, c*d, d*c, c*d*c}`
  (dim 4)** (live-verified). `HH_*(A) = [4,1,2,3,3,3,4,6,8,9]` **grows** (`A` has infinite gl.dim)
  while `HH_*(B) = [3,0,0,…]` (`gl.dim B = 2`, live-verified) — they diverge because `A/B` is **not**
  tensor nilpotent (the relative path `d*c` (`3→3`) has `(d*c)^{⊗_B m} ≠ 0` for all `m`, CLMS
  Ex. 5.5's `cd ⊗ cd = cdc ⊗ d`). **This does NOT mean "Han transport fails" in the sense of Han's
  conjecture failing (M3):** both `A` and `B` DO satisfy Han's conjecture individually — `B` has
  finite gl.dim, and `A` satisfies it **vacuously** (its `HH_*` does not vanish, so the conjecture's
  premise is never met). What fails is that the **extension is not bounded**, so the
  bounded-extension transport machinery does not apply — and, consistently, `HH_*(A) ≠ HH_*(B)`.
  (Note the raw dims here happen to obey `dim HH_m(B) ≤ dim HH_m(A)`, but this is NOT the
  theorem-guaranteed injection — that injection is only claimed under tensor-nilpotency, which fails
  here; it does not exercise the strict-injection case, see Task III4.) **HARD literature pin**
  (`oracle_literature`): the recognizer must certify Ex. 5.5 **not bounded** (tensor-nilpotency
  `False`, witness `d*c`), and `HanTransport.transport = "not_bounded"`.
- **The `A/B` basis and every relative path are read directly off `A.basis_labels`** by filtering
  on new-arrow usage (live-verified for both examples) — the recognizer needs no new basis engine.
- **Bar HH blows up early** (Ex. 5.3 bar `b_2` is `2162 × 99452 > max_cells`); the batteries use
  `engine="auto"` over `GF(p)` (fast/minimal/Bardzell) and `engine="cs"` over `QQ`.

**Could NOT live-verify (to-be-verified-at-implementation).**
- The **tensor powers `(A/B)^{⊗_B m}` and the relative-cycle / `J`-interrupter machinery** — building
  them IS the plan's work; the dims + the `J`-interrupter combinatorics are hand-checkable on the
  two examples (Ex. 5.3 index 2; Ex. 5.5 witness `d*c`) and pinned there.
- **`pd_{B^e}(A/B)`** — for Ex. 5.3/5.5 this is **certified finite FREE** by `gl.dim B = 2`
  (live-verified; CLMS Ex. 6.1, B1) — the enveloping-algebra fallback is not exercised by the
  primary oracles. The `enveloping_algebra(B)` construction (does not exist today; see Design
  findings) is only needed for the `gl.dim B = ∞` case.
- **`relative_homology(ext, top)` (`HH_*(A|B)`)** — the CLMS finite complex is new; its output on
  Ex. 5.3 must satisfy `HH_*(A|B) = 0` for `* ≥ 2` (Cor. 2.4) and slot into the JZ sequence, whose
  injection `HH_*(B) ↪ HH_*(A)` is here an equality (`[5,0,…] = [5,0,…]`, consistent with bounded).
- **A tensor-nilpotent-but-UNBOUNDED example** (to exercise the strict injection `dim HH_m(B) < dim
  HH_m(A)`) — NOT constructed at spec time (see Task III4; honest "unexemplified" note there).

### Design findings (substrate gaps — the plan adds the constructors)

1. **No subalgebra representation exists.** No `remove_arrow`/`Subalgebra`/arrow-removal builder
   anywhere in `src/quiverlab` (`quotient_by_idempotent` deletes by *vertex* set and is a
   *quotient*, not a subalgebra). BUT `families/_present.py::present_from_pi(Q, img, T, dom,
   dim_expected, …)` is the shipped, reusable presenter (given a quiver + generator images in a
   target, it extracts `I = ker(π)` by length-lex and certifies the dimension). **Task 0 builds
   `arrow_removal_subalgebra` on top of it.**
2. **The finite-`pd_{B^e}` leg is FREE when `gl.dim B < ∞` (B1 — the primary route).** CLMS Ex. 6.1:
   `gl.dim B < ∞ ⟹ gl.dim B^e < ∞ ⟹ pd_{B^e}(A/B) < ∞` for every `B`-bimodule. The shipped
   `global_dimension` decides `gl.dim B` (live-verified `= 2` for Ex. 5.3 AND Ex. 5.5), so the leg
   needs **no enveloping algebra** on the common finite-gl.dim case. **Only the fallback
   (`gl.dim B = ∞`)** needs `enveloping_algebra(B)`: `TensorProduct(B, B.opposite())` yields a
   **quiver-less structure-constant** algebra the module `pd` stack refuses
   (`builders._require_provenance`); there is no product-quiver builder and no bimodule-as-`B^e`-module
   bridge — so the fallback is a new substrate task (Task Group II), demoted below the gl.dim route
   and honestly capped.
3. **Relative HH over a general `B` is refused** — the shipped `hochschild/relative.py` computes
   only `HH_*(A|kQ_0)` (the separable vertex subalgebra `E`, where relative `=` absolute); a general
   `B` raises `QuiverlabError("relative HH over B is implemented only for B = kQ_0 in v1")`. **P73
   delivers `HH_*(A|B)` for a general `B` via the CLMS finite complex (Task Group III).**
4. **Bonus oracle tie-in (CLMS Ex. 5.4).** For `E = kQ_0 ⊂ Λ` (the split vertex-subalgebra
   extension), `A/E = rad Λ` is `E`-tensor nilpotent **iff `Q` has no oriented cycles** — a purely
   combinatorial oracle that cross-checks leg (i) against the shipped `relative.py`'s separable
   route and against acyclicity.

---

## Scope gates (contractual; every boundary is a loud typed refusal or an honest `"undecided"`)

| surface | scope | refusal / verdict |
|---|---|---|
| **input model (M1)** — `B ⊆ A` by arrow subset | only extensions **whose `B` sits on a sub-quiver of `Q_A`** (`B` generated by kept arrows of the user's `A`) | a bounded extension whose `B` is NOT generated by a subset of `A`'s arrows is **inexpressible** — the tool discloses this loudly (it does not silently claim "fully general"); re-present `A` on a quiver exhibiting `B`'s arrows to bring it in scope |
| `arrow_removal_subalgebra(A, new_arrows)` | `A.quiver is not None`; `new_arrows ⊆ (Q_A)_1` | structure-constant `A`, or an unknown arrow name → loud `QuiverlabError` |
| `is_tensor_nilpotent` (capped route) | any f.d. `A` with a quiver | cap reached, no certificate → `status="undecided"` (never a guessed verdict) |
| `is_tensor_nilpotent` (certificate route) | any f.d. `A` with a quiver | no relative cycles / all `J`-interrupted → `"nilpotent"` (`True`); a non-`J`-interrupted surviving cycle → `"not_nilpotent"` (`False`, witness) |
| `one_sided_projective` | `A/B` as a one-sided `B`-module (module stack) | `char ≤ dim` over `GF(p)` → loud (the shipped `pd`/`is_isomorphic` caveat); presentation-less → loud |
| finite `pd_{B^e}` leg — **primary** (`gl.dim B < ∞`) | any `B` where `global_dimension` returns finite | `gl.dim B < ∞` → `pd_{B^e}(A/B) < ∞` certified FREE (CLMS Ex. 6.1); no enveloping algebra |
| finite `pd_{B^e}` leg — **fallback** `enveloping_algebra(B)` (`gl.dim B = ∞`) | `B` quiver-presented; product-quiver within `dim (dim B)^2 ≤` budget; `pd` within `pd_cap` | over budget / unresolved within `pd_cap` → `status="undecided"` + certified lower bound (never `∞` unproven) |
| `relative_homology(ext, top)` (`HH_*(A|B)`) | `A/B` tensor-nilpotent (finite complex) OR `top <` a stated bound | not-nilpotent + `top` past the truncation → loud (the complex is infinite; honest `status`) |
| `bounded_extension` / `han_transport` | quiver-presented `A`; a fixed side for legs (ii)/(iii) | any leg `"undecided"` → overall `"undecided"`; leg (i) `False` → `"not_bounded"` |
| both, over `GF(p)` | **char 0 or char > dim** where the module/`is_isomorphic` primitives are used | `char ≤ dim` → loud `QuiverlabError` from the shipped primitives (never a silent wrong verdict) |

All batteries run over **QQ** for the certificate legs (exactness) with a `GF(32003)` HH parity
spot-check on the examples (distinct dim-vectors, no `is_isomorphic` raise); HH over `QQ` routes
`engine="cs"`, over `GF(p)` routes `engine="auto"`.

---

## API surface (public via `import quiverlab`; exact only)

```python
# --- families/extension.py: the input model + relative-path combinatorics ---
def arrow_removal_subalgebra(A, new_arrows) -> Extension        # B <= A by removing new_arrows
def enveloping_algebra(B) -> Algebra                            # B^e = B (x) B^op, bound-quiver presented
@dataclass
class Extension:
    A: object
    B: object
    new_arrows: tuple                    # F
    dim_A: int; dim_B: int; dim_quotient: int
    def relative_paths(self) -> list     # A/B basis: A-paths using >= 1 new arrow (labels + coords)
    def relative_cycles(self, cap: int) -> list          # CLMS Def. 5.11
    def quotient_bimodule(self) -> object                # A/B as a hochschild.coefficients-style B-bimodule

# --- invariants/han.py: the recognizer + transport + certificate ---
def is_tensor_nilpotent(ext, *, cap=8) -> TensorNilpotency
def one_sided_projective(ext, *, side="left") -> OneSidedProjectivity
def bounded_extension(A, new_arrows, *, side="auto", nilp_cap=8, pd_cap=16) -> BoundedCertificate
def han_transport(A, new_arrows, *, side="auto", nilp_cap=8, pd_cap=16, hh_top=None) -> HanTransport
def han_transport_block(A, spec, *, budget=None) -> dict       # the runner JSON block
@dataclass
class TensorNilpotency:
    status: str                          # "nilpotent" | "not_nilpotent" | "undecided"
    index: int | None                    # nilpotency index n when "nilpotent"
    route: str                           # "no_relative_cycles" | "J_interrupter" | "direct_cap" | "witness"
    witness: object | None               # the surviving relative cycle when "not_nilpotent"
    cap: int; note: str
@dataclass
class BoundedCertificate:
    bounded: bool | None                 # True (a side) / False / None(undecided)
    side: str | None                     # "left" | "right" | None
    tensor_nilpotent: TensorNilpotency
    one_sided: OneSidedProjectivity
    pd_Be: object                        # {route:"gldim"|"enveloping", value:int|None, status:str} -- gldim primary (B1)
    note: str
@dataclass
class HanTransport:
    # the injection/iso LADDER (H1): NO "isomorphism from leg (i) alone" verdict exists.
    transport: str                       # "bounded"(iso) | "pd_injection" | "nilpotent_injection" | "not_bounded" | "undecided"
    certificate: BoundedCertificate
    han_B: bool | None; han_A: bool | None   # honest per-side Han status where derivable
    injection_from: int | None           # CERTIFIED LOWER BOUND (H3): a degree from which the injection/iso is
                                         # guaranteed; the THEOREM states only *>>0 (bounded-exactness degree,
                                         # NOT the nilpotency index) -- reported as a lower bound, not a precise threshold
    injection_bound_ok: bool | None      # self-cert gate: dim HH_m(B) <= dim HH_m(A) held over the checked range
    references: list; note: str

# --- hochschild/jacobi_zariski.py: the JZ computational tool ---
def relative_homology(ext, top, *, coefficients=None) -> object    # HHTable of HH_*(A|B, X)
def jacobi_zariski_sequence(ext, top) -> JZSequence
@dataclass
class JZSequence:
    hh_A: list; hh_B: list; hh_rel: list         # dim sequences
    exact_at_B: list[bool]; exact_at_rel: list[bool]     # the two guaranteed spots (self-cert)
    gap_at_A: list[int]                          # Ker/Im dim at the HH_*(A) spot
    note: str

# --- Algebra delegates (core/algebra.py, thin lazy-import) ---
Algebra.enveloping()                                           -> Algebra
Algebra.bounded_extension(new_arrows, *, side="auto", ...)     -> BoundedCertificate
Algebra.han_transport(new_arrows, *, ...)                      -> HanTransport
Algebra.relative_homology(new_arrows, top)                     -> HHTable      # HH_*(A|B)
```

`TensorNilpotency`/`BoundedCertificate`/`HanTransport`/`JZSequence` are **data reports** — no
`__bool__` (mirroring `ARQuiver`/`ExchangeGraph`); read the `status`/`transport` field.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **Consumes (verified live at spec time on dev-tip):** `Quiver.algebra(relations=…, field=…)`
  (Groebner accepts the inhomogeneous `abcd − αβ`); `Algebra.hochschild_homology/cohomology(top,
  engine=…, coefficients=…)` with `engine ∈ {"auto","bar","fast","cs"}`; `hochschild.coefficients.
  Bimodule` (regular/dual/twisted/`from_actions`); `Algebra.opposite()`;
  `families.tensor.TensorProduct`; `families/_present.py::present_from_pi`; `modules.module.Module`
  (`side="right"|"left"`) + `Module.projective_resolution(length).pd()` (int or `None`);
  `modules.tor`/`modules.ext`; `hochschild/relative.py` (the `E = kQ_0` separable route, as the
  Ex.-5.4 cross-oracle). Branch `plan-73-han-bounded-extension` off `dev`.
- **Honest semi-decision contract (metaplan §1.3).** Every leg is **decided** iff its scope holds
  and its cap did not trip; otherwise `status = "undecided"` with the cap + a note — **never a
  guessed index / verdict**. Refusals are `QuiverlabError`.
- **No floats in `src/`.** Dims/indices/verdicts `int`/`bool`/`str`; HH tables integer lists;
  actions exact-Domain matrices. Client-side float only (`gui.js`, exempt).
- **Composition is left-to-right** (`a*b` = first `a` then `b`) — the CLMS right-to-left paths are
  translated in the examples. `_assert_comparable`/`_require_provenance` guard the module + presenter
  calls.
- **Plan-32 markers** (orthogonal; `oracle_*` FORBIDDEN in `tests/{webapp,gui,hpc}`):
  - `oracle_literature`: Ex. 5.3 certified **bounded** (nilpotency index 2, left `B`-projective,
    finite `pd_{B^e}`, transport `iff`); Ex. 5.5 certified **not bounded** (nilpotency `False`,
    witness `d*c`); the dims `(47,20,27)` / `(10,6,4)`; Ex. 5.4 `rad` `E`-nilpotent ⇔ `Q` acyclic.
  - `oracle_crossengine`: the JZ transport witness `HH_*(A) ≅ HH_*(B)` for `* ≥ n` computed by the
    **shipped HH engines** vs the tensor-nilpotency prediction — agree on Ex. 5.3, correctly
    **disagree** on Ex. 5.5 (a discriminating oracle); `relative_homology(ext,top)` (`HH_*(A|B)`
    via the CLMS complex) `= 0` for `* ≥ n` and slots into the JZ sequence; `enveloping_algebra(B)`
    HH ≡ the shipped internal `A^e` route on small `B`.
  - `oracle_selfcert`: the `J`-interrupter certificate ⟹ `(A/B)^{⊗_B m} = 0` **directly recomputed**
    at `m = index`; the JZ sequence **"exact twice in three"** (exactness at the `HH_*(B)` and
    `HH_*(A|B)` spots, verified from the three dim sequences); `dim B` certified via `present_from_pi`;
    `dim B^e = (dim B)^2`; `arrow_removal_subalgebra` round-trip (`#(F-free A-basis) = dim B`); the
    cap/`"undecided"` honesty (a constructed input that trips the cap returns `"undecided"`, never a
    verdict).
  - `qpa`: QPA 1.37 has **no** bounded-extension / Han-conjecture / Jacobi–Zariski surface (live
    `NamesGVars()` guard that FAILS if that ever changes — the P35-products precedent); the
    **component** dims are QPA-anchored — `HH_*(A)`/`HH_*(B)` via QPA `HochschildHomology…` and the
    one-sided `pd_B(A/B)` via QPA `ProjectiveResolution` on the example algebras (input-level
    anchors, not a verdict compare).
- **Mid-merge-train counts.** Absolute suite counts drift; the verification task **recounts** the
  oracle-class table at merge time (`tests/release/test_oracle_classes.py`; paste live numbers,
  claim only this plan's deltas).
- Every merge updates `docs/verification.md` (new oracle rows + recounted class table green) and
  adds citations to `references.bib` + `registry.py` (`bibtex()` hard-fails if the two disagree).
  Conventional commits; green at every commit.

**No-red-commits discipline.** Each Task ships a self-contained green slice: the subalgebra
constructor before the recognizer before the JZ tool before the GUI. Fields default to
`None`/`"undecided"` so intermediate task boundaries are green (no forward reference, no `xfail`).

---

# Task 0 (BLOCKING) — the arrow-removal subalgebra constructor

**Why FIRST.** Every downstream leg reads `B`, `A/B`, and the relative paths off this object; the
dimension certificate is the extension axiom `J ∩ B = 0`.

**Files:** create `src/quiverlab/families/extension.py`; test `tests/families/test_extension.py`.

- [ ] **Step 1: Failing tests (self-cert).** Build Ex. 5.3 / Ex. 5.5 as `A`; call
  `arrow_removal_subalgebra(A, new_arrows=("a",))` / `("d",)`; assert `ext.dim_B == 20` / `6`,
  `ext.dim_quotient == 27` / `4`, `set(ext.new_arrows) == {"a"}` / `{"d"}`, and that
  `ext.relative_paths()` returns exactly the new-arrow-using basis paths (Ex. 5.5:
  `{"d","c*d","d*c","c*d*c"}`). A structure-constant `A` or an unknown arrow name → `QuiverlabError`.
- [ ] **Step 2:** confirm failure (`ModuleNotFoundError`).
- [ ] **Step 3: Implement.** `Q_B = Quiver(Q_0, {kept arrows})`; drive `present_from_pi(Q_B,
  img={kept arrow ↦ its coordinate in A}, T=A, dom, dim_expected=#(F-free lin-indep A-basis))` to
  get `B = kQ_B/ker(π)`, **certifying** `dim B`. If `present_from_pi` needs `dim_expected`
  up-front, compute it first by row-reducing the `F`-free `A`-basis coordinates (a shipped
  `fields.linalg` rank). `relative_paths()` filters `A.basis_labels` on new-arrow usage (tokenise
  on `*`; guard the `al` vs `a` prefix trap — match whole arrow tokens, live-verified). `Extension`
  caches `dim_A`, `dim_B`, `dim_quotient`, and the `B`-bimodule action of `A/B`
  (`quotient_bimodule()` via the shipped `hochschild.coefficients.Bimodule.from_actions`-style
  Lact/Ract restricted to `B`).
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(families): arrow-removal subalgebra B <= A (extension by arrows and relations, CLMS Def 5.2) via present_from_pi; certified dim B (J cap B = 0), relative-path A/B basis (Plan 73/R7)`.

---

# Task group I — tensor nilpotency (leg i): capped semi-decision + certificate (R7 core)

### Task I1: relative cycles + `J`-interrupters

**Files:** modify `src/quiverlab/families/extension.py`; create `src/quiverlab/invariants/han.py`;
test `tests/invariants/test_han_nilpotent.py`.

- [ ] **Step 1: Failing tests.** `ext.relative_cycles(cap=…)` on Ex. 5.3 returns the cycle
  `a·(bcd)·a` (one new arrow, closing at vertex 1); on Ex. 5.5 the cycle at vertex 3 through `d`.
  A `J`-interrupter predicate flags `a` in Ex. 5.3 (`(bcd)·a·(bcd) ∈ B`) and reports **no**
  `J`-interrupter for the `d*c` cycle in Ex. 5.5.
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement** `relative_cycles` (CLMS Def. 5.11 — enumerate concatenable
  relative-path necklaces, cap on `F`-length) and `_is_J_interrupter` (Def. 5.13 — the "product
  across the arrow" reduces into `B`, tested by `A.multiply` + the `relative_paths` membership).
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(invariants): relative cycles + J-interrupters (CLMS Def 5.11/5.13) for the Han bounded-extension recognizer (Plan 73/R7)`.

### Task I2: `is_tensor_nilpotent` — the capped semi-decision + certificate + witness

**Files:** modify `src/quiverlab/invariants/han.py`; extend the test.

- [ ] **Step 1: Failing tests (the HARD pins).**
```python
lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert

@lit
def test_ex53_tensor_nilpotent_index_2():
    ext = arrow_removal_subalgebra(build_ex53(QQ), ("a",))
    tn = is_tensor_nilpotent(ext, cap=8)
    assert tn.status == "nilpotent" and tn.index == 2       # length index 2 (CLMS Thm 5.20 disc.)
    assert tn.route in ("J_interrupter", "no_relative_cycles", "direct_cap")

@lit
def test_ex55_not_tensor_nilpotent_witness():
    ext = arrow_removal_subalgebra(build_ex55(QQ), ("d",))
    tn = is_tensor_nilpotent(ext, cap=8)
    assert tn.status == "not_nilpotent"                     # CLMS Ex 5.5
    assert tn.witness is not None                           # the surviving relative cycle (d*c)

@selfcert
def test_certificate_forces_vanishing():
    """J-interrupter route ⟹ (A/B)^{⊗_B index} recomputed directly == 0."""
    ext = arrow_removal_subalgebra(build_ex53(QQ), ("a",))
    tn = is_tensor_nilpotent(ext, cap=8)
    assert ext.tensor_power_dim(tn.index) == 0 and ext.tensor_power_dim(tn.index - 1) != 0

@selfcert
def test_cap_returns_exactly_undecided_not_a_lie():
    ext = arrow_removal_subalgebra(build_ex55(QQ), ("d",))  # (A/B)^{⊗ m} never vanishes
    tn = is_tensor_nilpotent(ext, cap=3, use_certificate=False)  # force the capped-only route
    # The cap-only route sees 3 nonzero powers and CANNOT conclude non-nilpotency (powers may
    # grow before dying) -- it must return EXACTLY "undecided", never "not_nilpotent", never True.
    assert tn.status == "undecided"
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement.** Certificate route first (no relative cycles / all `J`-interrupted ⟹
  `"nilpotent"`, index = length index; a surviving non-`J`-interrupted cycle ⟹ `"not_nilpotent"`
  with the witness). Direct route (`ext.tensor_power_dim(m)` by relative-path concatenation + `J`
  reduction + the `F`-length-0 collapse) as the cap-bounded cross-check / fallback; `"undecided"`
  at the cap. **The certificate and the direct route must AGREE where both decide** (a self-cert).
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(invariants): is_tensor_nilpotent -- capped semi-decision + J-interrupter certificate + Ex-5.5 non-vanishing witness (Plan 73/R7); pins Ex 5.3 index 2, Ex 5.5 not-nilpotent`.

---

# Task group II — one-sided projectivity (leg iii) + finite pd over B^e (leg ii)

### Task II1: finite `pd_{B^e}(A/B)` leg — the `gl.dim B` PRIMARY route (B1)

**Files:** modify `src/quiverlab/invariants/han.py`; test `tests/invariants/test_han_pd.py`.

- [ ] **Step 1: Failing tests (the PRIMARY leg-(ii) route — free when `gl.dim B < ∞`).**
```python
lit = pytest.mark.oracle_literature
@lit
def test_ex53_pd_Be_finite_via_gldim():
    ext = arrow_removal_subalgebra(build_ex53(QQ), ("a",))
    pd = finite_pd_Be(ext)                                   # CLMS Ex 6.1 route
    assert pd["route"] == "gldim" and pd["status"] == "finite"   # gl.dim B = 2 (live-verified) -> pd_{B^e} < inf
    # NO enveloping algebra is constructed on this path.
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement `finite_pd_Be(ext)`.** PRIMARY: `gd = B.global_dimension()`; if finite ⟹
  `{"route":"gldim","status":"finite","gldim_B":gd}` (CLMS Ex. 6.1: `gl.dim B < ∞ ⟹ gl.dim B^e < ∞
  ⟹ pd_{B^e}(A/B) < ∞`). Only when `gl.dim B = ∞` fall through to the Task-II2 enveloping fallback.
  Live-verified `gl.dim B = 2` (exact) for both Ex. 5.3 and Ex. 5.5.
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(invariants): finite pd_{B^e}(A/B) PRIMARY route via gl.dim B < inf (CLMS Ex 6.1) -- free, no enveloping algebra (Plan 73/R7 B1)`.

### Task II2: `enveloping_algebra(B)` — `B^e` as a bound quiver algebra (the `gl.dim B = ∞` FALLBACK)

**Files:** modify `src/quiverlab/families/extension.py` (or a new `families/enveloping.py`); the
`Algebra.enveloping` delegate; test `tests/families/test_enveloping.py`. **Demoted below Task II1**
— only reached when `gl.dim B = ∞`; the primary oracles (Ex. 5.3/5.5) do NOT exercise it.

- [ ] **Step 1: Failing tests (self-cert).** `Be = enveloping_algebra(B)` has `Be.quiver is not
  None`, `Be.dim == B.dim**2`, product-quiver vertices `= (Q_B)_0 × (Q_B)_0`; and the shipped
  internal `A^e` Hochschild route agrees with an `enveloping`-based `pd` on a tiny `B` (e.g.
  `k[x]/(x^2)`: `pd_{B^e}(B) =` Hochschild dim, cross-checked).
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement.** Product quiver `Q_B ⊗ Q_B^op` (vertices `(u,v)`; arrows
  `(α, e)` and `(e, β^op)`); present via `present_from_pi` against the shipped
  `TensorProduct(B, B.opposite())` structure constants, certifying `dim = (dim B)^2`. Loud refusal
  when `(dim B)^2` exceeds the presenter budget (honest — this is the heavy fallback leg).
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(families): enveloping_algebra(B) = B (x) B^op as a bound-quiver Algebra (product quiver, present_from_pi-certified dim (dim B)^2) -- first-class A^e, the gl.dim-infinite fallback (Plan 73/R7)`.

### Task II3: `one_sided_projective` + the combined finite-`pd_{B^e}` leg (primary + capped fallback)

**Files:** modify `src/quiverlab/invariants/han.py`; extend the tests.

- [ ] **Step 1: Failing tests.** `one_sided_projective(ext, side="left")` on Ex. 5.3 returns
  `projective=True` (CLMS Thm 5.20; live cross-check: `pd_B(A/B \text{ left}) == 0`). `finite_pd_Be`
  on Ex. 5.3 returns `route="gldim", status="finite"` (via Task II1); a constructed `gl.dim B = ∞`
  input routes to the enveloping fallback and returns a finite value within `pd_cap` (bounded) or an
  honest `"undecided"` + certified lower bound — **never `∞` unproven**.
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement.** `one_sided_projective`: build `A/B` as a one-sided `B`-module
  (`Module(B, dim_quotient, action, side=…)` from the restricted relative-path action),
  `pd_B = M.projective_resolution(bound).pd()`; `== 0` ⟹ projective. Also wire the Thm 5.20
  combinatorial certificate as a fast path. `finite_pd_Be` combines the two routes: gl.dim primary
  (Task II1), then the `enveloping_algebra(B)` fallback (`A/B` as a right `B^e`-module,
  `pd = M.projective_resolution(pd_cap).pd()`; `int` ⟹ certified finite; `None` ⟹ `"undecided"` +
  the reached-length lower bound).
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(invariants): one-sided B-projectivity (pd_B(A/B)=0 + CLMS Thm 5.20) + finite pd_{B^e}(A/B) combined (gl.dim primary, enveloping capped fallback) (Plan 73/R7)`.

---

# Task group III — Jacobi–Zariski tool + Han transport + self-cert gate

### Task III1: `relative_homology` — `HH_*(A|B)` via the CLMS finite complex

**Files:** create `src/quiverlab/hochschild/jacobi_zariski.py`; test
`tests/hochschild/test_jacobi_zariski.py`.

- [ ] **Step 1: Failing tests.** On Ex. 5.3 (`A/B` tensor-nilpotent, index 2),
  `relative_homology(ext, top=6)` returns `HH_*(A|B)` with `dim = 0` for `* ≥ 2` (Cor. 2.4). On a
  non-nilpotent input past the truncation → loud `QuiverlabError` (honest — the complex is
  infinite). `HH_0(A|B)` matches `dim (A/B)_B`-coinvariants.
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement** the CLMS normalized relative bar complex (Thm 2.2):
  terms `X ⊗_{B^e} (A/B)^{⊗_B m}`, differential `b` (`σ`-independent), homology by exact rank
  (shipped `fields.linalg`). `coefficients=None` ⟹ `X = A` regular. Finite when `A/B` nilpotent
  (`m < index`); else honest truncation to `top`.
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(hochschild): HH_*(A|B) over a GENERAL subalgebra B via the CLMS normalized relative bar complex (finite when A/B tensor-nilpotent) -- the Jacobi-Zariski tool (Plan 73/R7)`.

### Task III2: `jacobi_zariski_sequence` + `han_transport` + the self-cert gate

**Files:** modify `src/quiverlab/hochschild/jacobi_zariski.py` and `src/quiverlab/invariants/han.py`;
the `Algebra` delegates; test `tests/invariants/test_han_transport.py`.

- [ ] **Step 1: Failing tests (the transport pins).**
```python
lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert

@lit
def test_ex53_bounded_transport():
    ht = build_ex53(QQ).han_transport(("a",), hh_top=9)
    assert ht.transport == "bounded"                        # all three legs (a side) certified -> ISO
    assert ht.certificate.tensor_nilpotent.index == 2
    assert ht.certificate.pd_Be["route"] == "gldim"         # finite pd_{B^e} FREE via gl.dim B = 2 (B1)
    assert ht.injection_from is not None                    # certified LOWER BOUND (H3), not a precise threshold
    assert ht.injection_bound_ok is True                    # dim HH_m(B) <= dim HH_m(A); EQUALITY under "bounded"

@lit
def test_ex55_not_bounded():
    ht = build_ex55(QQ).han_transport(("d",), hh_top=9)
    assert ht.transport == "not_bounded"                    # Ex 5.5: A/B not tensor-nilpotent (leg i False)

@xeng
def test_injection_bound_self_cert_gate():
    """The HONEST JZ consequence on the SHIPPED HH engines is the INJECTION BOUND dim HH_m(B) <=
    dim HH_m(A) (NOT an isomorphism from nilpotency). Ex 5.3 is BOUNDED so it is an equality; on the
    NOT-bounded Ex 5.5 the sequences differ (and 5.5 is not a theorem-injection -- A/B not nilpotent)."""
    A53, B53 = build_ex53(GF(32003)), build_ex53_B(GF(32003))
    a53, b53 = A53.hochschild_homology(9).dims, B53.hochschild_homology(9).dims
    assert all(bi <= ai for ai, bi in zip(a53, b53))        # injection bound
    assert a53 == b53                                        # EQUALITY only because 5.3 is bounded (iso)
    A55, B55 = build_ex55(GF(32003)), build_ex55_B(GF(32003))
    assert A55.hochschild_homology(9).dims != B55.hochschild_homology(9).dims   # not bounded -> differ

@selfcert
def test_injection_ladder_logic():
    """W3(b): exercise the LADDER LOGIC itself (not just live HH). Given per-leg verdicts, the
    transport field is set per the H1 ladder -- independent of any concrete algebra's HH."""
    assert transport_verdict(nilp=True, pd_finite=True,  one_sided=True ) == "bounded"
    assert transport_verdict(nilp=True, pd_finite=True,  one_sided=False) == "pd_injection"
    assert transport_verdict(nilp=True, pd_finite=None,  one_sided=None ) == "nilpotent_injection"
    assert transport_verdict(nilp=False, pd_finite=None, one_sided=None ) == "not_bounded"
    assert transport_verdict(nilp=None,  pd_finite=None, one_sided=None ) == "undecided"

@selfcert
def test_jz_exact_twice_in_three():
    jz = jacobi_zariski_sequence(arrow_removal_subalgebra(build_ex53(QQ), ("a",)), top=6)
    assert all(jz.exact_at_B) and all(jz.exact_at_rel)      # the two guaranteed spots
    # NOTE: Ex 5.3 has trivial higher HH ([5,0,...]); a nontrivial-higher-HH nilpotent example is
    # exercised in Task III4 (Ex 5.5 CANNOT serve -- it is not tensor-nilpotent, so HH_*(A|B) is
    # an infinite complex; see Task III4's justification).
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement.** `transport_verdict(nilp, pd_finite, one_sided)` is the pure ladder
  function (H1): `nilp is False -> "not_bounded"`; `nilp is None -> "undecided"`; then
  `(pd_finite and one_sided) -> "bounded"`; `pd_finite and not one_sided -> "pd_injection"`; else
  `"nilpotent_injection"` (leg i only); any needed leg `"undecided"` -> `"undecided"`.
  `jacobi_zariski_sequence`: assemble `HH_*(A)`/`HH_*(B)` (shipped engines) + `HH_*(A|B)` (Task III1);
  verify exactness at the two guaranteed spots by rank bookkeeping; report the gap at the `HH_*(A)`
  spot. `han_transport`: run the three legs (`bounded_extension`), set `transport` via
  `transport_verdict`; the **self-cert gate** checks the **injection bound `dim HH_m(B) ≤ dim
  HH_m(A)`** over the computed range (`injection_bound_ok`) and asserts **equality only when
  `transport == "bounded"`**; record `injection_from` as a CERTIFIED LOWER BOUND (H3 — the theorem
  gives only `* ≫ 0`, governed by the bounded-exactness degree, not the nilpotency index; report a
  bound, never a precise threshold). Scale `han_A`/`han_B` honestly (claim the `iff` only when
  `"bounded"`). Delegates.
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(invariants,hochschild): Jacobi-Zariski sequence + han_transport injection/iso LADDER (CLMS Thm 3.1/4.6) with the injection-bound self-cert gate -- Ex 5.3 bounded (iso), Ex 5.5 not-bounded (Plan 73/R7 H1)`.

### Task III3: the CLMS Ex. 5.4 combinatorial oracle (bonus tie-in)

**Files:** extend `tests/invariants/test_han_nilpotent.py`.

- [ ] **Step 1:** `E = kQ_0 ⊂ Λ` (the vertex subalgebra; `new_arrows = all of (Q)_1`): assert
  `is_tensor_nilpotent` returns `"nilpotent"` **iff `Q` is acyclic** — pin an acyclic `Λ` (linear
  `A_3`) and a cyclic one (a loop / an oriented cycle). Cross-check against the shipped
  `hochschild/relative.py` separable route (`E`-relative `=` absolute for the acyclic case).
- [ ] **Step 2–4:** implement + run; **Step 5:** commit
  `test(invariants): CLMS Ex 5.4 -- rad is E-tensor-nilpotent iff Q acyclic (E=kQ_0 vertex subalgebra), cross-checked vs the shipped relative.py separable route (Plan 73)`.

### Task III4: the STRICT-injection / nontrivial-higher-HH example (H2 + W3c) — construct or document impossibility

**Why.** Ex. 5.3 (bounded) gives the *equality* case and Ex. 5.5 is *not tensor-nilpotent* (so it
exercises neither the theorem-guaranteed injection nor the JZ finite complex at higher degrees). The
ladder's **`"pd_injection"` / `"nilpotent_injection"` rows** and the **strict** injection
`dim HH_m(B) < dim HH_m(A)` are so far **unexemplified**, and the JZ exactness test runs only on a
trivial-higher-HH algebra. This task fills both gaps with ONE example if possible.

**Target.** A **tensor-nilpotent-but-UNBOUNDED** extension `B ⊂ A` (leg (i) holds; NOT bounded —
either `A/B` not one-sided projective, giving `"pd_injection"`/`"nilpotent_injection"`, so the
injection is **strict** at some degree), preferably with **`gl.dim B = ∞`** (nonzero higher
`HH_*(B)`) so the SAME example also exercises the JZ exactness at higher degrees and the
enveloping-algebra fallback (Task II2).

- [ ] **Step 1: Construct + live-verify (at implementation).** Candidate direction: a one-new-arrow
  monomial extension of a NON-finite-gl.dim `B` (e.g. `B` a symmetric Nakayama `kZ_n/J^L` or the
  dual numbers `k[x]/(x^2)`) whose single relative path is tensor-nilpotent (a `J`-interrupter kills
  it at some index) but is **not** left-`B`-projective. Compute `HH_*(A)`, `HH_*(B)`, `HH_*(A|B)`;
  assert `is_tensor_nilpotent == "nilpotent"`, `bounded_extension.bounded is False`, the ladder
  verdict `∈ {"pd_injection","nilpotent_injection"}`, and — the point — `dim HH_m(B) < dim HH_m(A)`
  **strictly** at some `m ≫ 0` (the injection is proper), with the JZ sequence exact at the two
  guaranteed spots in a degree where `HH_*(B)` is nonzero.
- [ ] **Step 2: If no such example is constructible within the arrow-extension model on small
  quivers, document the impossibility honestly.** Write the `oracle_selfcert` test as an explicit
  `pytest.mark.xfail(reason="unexemplified at spec time -- see plan Task III4")` **or** a passing
  test that asserts the ladder LOGIC (Task III2 `transport_verdict`) plus a docstring stating "no
  strict-injection instance found within scope; the ladder rows are logic-verified, the
  strict-injection case is theory-only" — NEVER a silent omission. Record the outcome on the
  verification page (honest-scope entry).
- [ ] **Step 3–4:** run; **Step 5:** commit
  `test(invariants): strict-injection / nontrivial-higher-HH JZ example -- tensor-nilpotent-but-unbounded (H2/W3c); construct + live-verify OR documented-impossibility xfail (Plan 73)`.

**Spec-time status (honest):** this example was **NOT constructed at authoring** — verifying
"tensor-nilpotent but not one-sided projective" needs the recognizer itself (unbuilt) and the
enveloping/relative machinery. The two live-built examples (5.3 bounded, 5.5 not-nilpotent) do not
cover the strict-injection case; the worker constructs it or documents impossibility per Step 2.

---

# Task group IV — GUI / webapp + verification

### Task G1: the `han_transport` algebra-level certificate kind (all three tiers, i18n ×4)

An **algebra-level certificate** kind carrying the **new-arrow subset `F`** as a request-level field
(the `derived_compare` second-input precedent — a kind that needs data beyond the algebra). GUI:
an **arrow-subset picker** on the drawn quiver (which arrows are "new").

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (dispatch branch for `han_transport` reading the
  `new_arrows`/`extension` field + `_snip` recipe), `docs/gui/runner.py` (the byte-identical twin),
  `webapp/server/schema.py` (a `new_arrows: list[str] | None` request field, guarded to the
  `han_transport` kind, canonicalizing through `canonical_key` — default `None` so every existing
  request's key is byte-stable), `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox
  `qlgui-han-transport`, the arrow-subset picker UI, `S.ids`, push-list, `renderBlock`,
  `scheduleProbe`, the `structure`/`invariants` THEME), `webapp/templates/draw.html`
  (`data-pick-kind-han_transport` + the arrow-picker), `webapp/server/i18n/{en,es,fr,zh}.json`
  (the block key chain + `pick.kind.han_transport`; parity gated by `tests/webapp/test_i18n.py`),
  `src/quiverlab/trace/results_html.py` (`_HEADINGS` + a render branch),
  `tests/webapp/_runner_goldens.json` + `test_runner_delegation.py` (ONE golden; existing
  byte-identical).
- Test: `tests/webapp/test_han_kind_p73.py`, `tests/gui/test_han_runner_twin_p73.py`.

**Block shape** (one shared `han_transport_block(A, spec, budget)`):
```python
{"kind": "han_transport", "new_arrows": ["a"],
 "dim_A": int, "dim_B": int, "dim_quotient": int,
 "tensor_nilpotent": {"status": str, "index": int|None, "route": str},
 "one_sided": {"side": str|None, "projective": bool|None},
 "pd_Be": {"route": "gldim"|"enveloping", "value": int|None, "status": str},   # gldim primary (B1)
 "transport": str,                  # "bounded" | "pd_injection" | "nilpotent_injection" | "not_bounded" | "undecided" (H1)
 "injection_from": int|None,        # certified lower bound (H3), NOT a precise threshold
 "note": str|None,
 "references": ["clms_bounded_extensions", "clms_jacobi_zariski", "kaygun_jacobi_zariski",
                "han_conjecture", "assem_book"]}
# refusal (structure-constant A / bad arrow name) -> {"error": msg, "references": [...]}, never a 500.
```
The JZ dim sequences + the per-degree exactness flags render in the worked-steps report
(`results_html`); the block itself carries counts + verdicts (cheap).

- [ ] **Step 1: Failing cross-runner tests** (unmarked — extras-gated dir): Ex. 5.3 block has
  `transport == "bounded"`, `tensor_nilpotent.index == 2`, `dim_B == 20`,
  `"clms_bounded_extensions"` in citation keys; Ex. 5.5 block `transport == "not_bounded"`; twin
  parity (`json.dumps(sort_keys=True)` across `hpc/spec.py` and `docs/gui/runner.py`); the
  `new_arrows` field canonical-key stability (an existing family request's key byte-unchanged).
- [ ] **Step 2: Implement** the dispatch + the `new_arrows` field in all three parse/schema sites,
  the twin, the `gui.js` checkbox + arrow-picker + `scheduleProbe` ETA (size like `recognizers` —
  algebra-dim `sizing_dim`; no estimator edit), the `_snip` recipe
  (`"han_transport": lambda it: "A.han_transport(%r)" % (it.new_arrows,)`),
  `_HEADINGS["han_transport"] = "Han transport (bounded extension)"` + the render branch (the
  three-leg table + the transport verdict + the JZ dim sequences).
- [ ] **Step 3: Add ONE golden** (`han_transport_ex53`) to `_runner_goldens.json`; note it in
  `test_runner_delegation.py`'s change-log docstring; confirm existing goldens byte-identical first.
- [ ] **Step 4:** run `tests/webapp/test_han_kind_p73.py tests/webapp/test_runner_delegation.py
  tests/gui/test_han_runner_twin_p73.py tests/webapp/test_i18n.py tests/hpc -q`.
- [ ] **Step 5:** commit
  `feat(gui,webapp,hpc): han_transport certificate kind (Plan 73) -- new-arrow-subset field, arrow picker, both runners byte-identical, i18n x4, one golden`.

**i18n keys — ALL FOUR LOCALES** (mirror `pick.kind.recognizers`):

| key | en | es | fr | zh |
|---|---|---|---|---|
| `pick.kind.han_transport` | Han transport (bounded extension) | Transporte de Han (extensión acotada) | Transport de Han (extension bornée) | Han 传递（有界扩张） |
| `block.han_transport.title` | Han transport | Transporte de Han | Transport de Han | Han 传递 |
| `block.han_transport.nilpotent` | A/B tensor-nilpotent | A/B tensor-nilpotente | A/B tenseur-nilpotent | A/B 张量幂零 |
| `block.han_transport.projective` | A/B one-sided B-projective | A/B B-proyectivo de un lado | A/B B-projectif d'un côté | A/B 单边 B-投射 |
| `block.han_transport.pd` | pd over B^e | pd sobre B^e | pd sur B^e | B^e 上的 pd |
| `block.han_transport.verdict` | Transport verdict | Veredicto del transporte | Verdict du transport | 传递判定 |
| `block.han_transport.newarrows` | New arrows (extension) | Flechas nuevas (extensión) | Flèches nouvelles (extension) | 新箭（扩张） |

### Task G2: verification page, citations, README, suite gate

**Files:** `src/quiverlab/citations/references.bib` + `registry.py`; `docs/verification.md`;
`README.md`; `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick P73); existing release gates.

- [ ] **Step 1: Citations (BibTeX-VERIFIED only).** Add the six keys resolved above. Ship the arXiv
  ids + DOIs (all verified this fix round); `# PIN` the volume/pages of `CLMSbounded2022` and the
  `Bull. LMS 54(3)` volume/number of `CLMSjacobiZariski2022` against the source landings at
  citation-add; delete any field that cannot be confirmed (house rule). **The JZ sequence is
  `2009.05017` (CLMS) + `1103.4377` (Kaygun, the classical origin) — NEVER `1908.11130`.**
```bibtex
@article{Han2006,
  author  = {Han, Yang},
  title   = {Hochschild (co)homology dimension},
  journal = {J. London Math. Soc. (2)}, volume = {73}, number = {3},
  pages   = {657--668}, year = {2006}, doi = {10.1112/S002461070602299X}}
@article{CLMSbounded2022,
  author  = {Cibils, Claude and Lanzilotta, Marcelo and Marcos, Eduardo N. and Solotar, Andrea},
  title   = {Han's conjecture for bounded extensions},
  journal = {J. Algebra}, year = {2022},
  doi = {10.1016/j.jalgebra.2022.01.022}, note = {arXiv:2101.02597}}   % # PIN volume/pages at citation-add
@article{CLMSjacobiZariski2022,
  author  = {Cibils, Claude and Lanzilotta, Marcelo and Marcos, Eduardo N. and Solotar, Andrea},
  title   = {Jacobi--Zariski long nearly exact sequences for associative algebras},
  journal = {Bull. Lond. Math. Soc.}, volume = {54}, number = {3}, year = {2022},
  doi = {10.1112/blms.12516}, note = {arXiv:2009.05017}}   % the JZ "exact twice in three" source (year 2022)
@article{Kaygun2012,
  author  = {Kaygun, Atabey},
  title   = {Jacobi--Zariski Exact Sequence for {H}ochschild Homology and Cyclic (Co)Homology},
  journal = {Homology Homotopy Appl.}, volume = {14}, number = {1}, pages = {65--78}, year = {2012},
  note    = {arXiv:1103.4377}}   % the classical noncommutative-JZ origin CLMS credit
@article{CLMSsplitBounded2020,
  author  = {Cibils, Claude and Lanzilotta, Marcelo and Marcos, Eduardo N. and Solotar, Andrea},
  title   = {Split bounded extension algebras and {H}an's conjecture},
  journal = {Pacific J. Math.}, volume = {307}, number = {1}, pages = {63--77},
  year = {2020}, doi = {10.2140/pjm.2020.307.63}, note = {arXiv:1908.11130}}
@misc{WangXuZhangZhou2024,
  author = {Wang, Ren and Xu, Xiaoxiao and Zhang, Jinbi and Zhou, Guodong},
  title  = {A recollement approach to {H}an's conjecture},
  year   = {2024}, note = {arXiv:2409.00945}}   % authorship VERIFIED this plan round
```
  and in `registry.py` (`_r(key, bibtex_key, kind, title, annotation, *tags)`):
```python
# NOTE: `han_conjecture` is ALREADY registered (registry.py:210) -- do NOT re-register; reuse the key.
_r("clms_bounded_extensions", "CLMSbounded2022", "algorithm",
   "Han's conjecture for bounded extensions",
   "CLMS: B subset A left/right bounded (A/B tensor-nilpotent, finite pd over B^e, one-sided "
   "B-projective) implies B satisfies Han iff A does (Thm 4.6). Examples 5.3 (bounded) / 5.5 "
   "(not bounded) are the P73 oracles.", "hochschild", "han"),
_r("clms_jacobi_zariski", "CLMSjacobiZariski2022", "foundation",
   "Jacobi-Zariski long nearly exact sequences for associative algebras",
   "CLMS: the Jacobi-Zariski long nearly exact sequence relating HH_*(A), HH_*(B), HH_*(A|B) "
   "('exact twice in three') -- the computational tool + self-cert gate for the Han transport.",
   "hochschild", "han"),
_r("kaygun_jacobi_zariski", "Kaygun2012", "foundation",
   "Jacobi-Zariski Exact Sequence for Hochschild Homology and Cyclic (Co)Homology",
   "Kaygun: the classical noncommutative Jacobi-Zariski sequence (B subset A with A/B flat) -- "
   "the origin CLMS credit for the noncommutative case; CLMS 2009.05017 is the 'long nearly "
   "exact / exact twice in three' refinement P73 computes with.", "hochschild"),
_r("clms_split_bounded", "CLMSsplitBounded2020", "foundation",
   "Split bounded extension algebras and Han's conjecture",
   "CLMS: the split-case predecessor of the bounded-extension theory (context).", "hochschild"),
_r("wang_recollement_han", "WangXuZhangZhou2024", "foundation",
   "A recollement approach to Han's conjecture",
   "Wang-Xu-Zhang-Zhou: an independent recollement/derived reduction of Han's conjecture "
   "(also proves Han for skew-gentle algebras -- ties to P68). Authorship verified.", "hochschild"),
```
  Verify BibTeX presence (`tests/citations/`) and that every key resolves.
- [ ] **Step 2: Verification page.** Add the P73 rows:
  - `families/extension.py` + `invariants/han.py` — `oracle_literature` (Ex. 5.3 bounded / index 2;
    Ex. 5.5 not-bounded / witness; dims `(47,20,27)` / `(10,6,4)`; Ex. 5.4 acyclicity);
    `oracle_selfcert` (`J`-interrupter ⟹ `(A/B)^{⊗ n}=0` recomputed; `dim B` / `dim B^e` certified;
    the cap `"undecided"` honesty).
  - `hochschild/jacobi_zariski.py` — `oracle_crossengine` (the **injection-bound** self-cert gate on
    the shipped engines: `dim HH_m(B) ≤ dim HH_m(A)`, an EQUALITY on the bounded Ex. 5.3, differing
    on the not-bounded Ex. 5.5; `HH_*(A|B)=0` for `* ≥ n`); `oracle_selfcert` (JZ "exact twice in
    three"; the `transport_verdict` ladder logic).
  - `qpa` — the `NamesGVars()` guard (FAILS if QPA ships a bounded-extension/Han/JZ surface) + the
    component HH / one-sided-`pd` input anchors on the example algebras.
  - **Honest-scope entries (binding):**
    (a) **Tensor-nilpotency is a capped semi-decision** with a certificate route; `"undecided"` is a
    first-class verdict (no a-priori index bound — powers can grow before dying).
    (b) **Finite `pd_{B^e}` — `gl.dim B < ∞` is the PRIMARY, free route** (CLMS Ex. 6.1;
    live-verified `gl.dim B = 2` on both examples); only `gl.dim B = ∞` falls to the capped
    enveloping-algebra fallback (`dim (dim B)^2`), `"undecided"` + certified lower bound over budget.
    (c) **The transport claim is an INJECTION/ISO LADDER, not a blanket agreement (H1):** leg (i)
    alone ⟹ `H_*(B,A) ↪ H_*(A,A)` (mixed coeff INJECTION); (i)+(ii) ⟹ `HH_*(B) ↪ HH_*(A)`
    (ordinary INJECTION); bounded ⟹ `HH_*(B) ≅ HH_*(A)` (ISO). **There is no "≅ from tensor-nilpotency
    alone"** — that would misread CLMS Thm 3.1. The self-cert gate is the injection bound
    `dim HH_m(B) ≤ dim HH_m(A)`, equality asserted only under `"bounded"`.
    (c′) **The `injection_from` threshold is a CERTIFIED LOWER BOUND (H3)**, not a precise degree:
    the theorem states `* ≫ 0` (bounded-exactness degree), not the nilpotency index.
    (c″) **The strict-injection case (`dim HH_m(B) < dim HH_m(A)` under a theorem-guaranteed
    injection) is UNEXEMPLIFIED at spec time (H2/W3c):** it requires a tensor-nilpotent-but-unbounded
    extension (Task III4); the worker constructs one and live-verifies, or documents impossibility
    within the arrow-extension model (the ladder LOGIC is `oracle_selfcert`-verified regardless).
    (d) **`HH_*(A|B)` over a general `B` is finite only when `A/B` is tensor-nilpotent** (Cor. 2.4);
    otherwise honestly truncated to `top` (the shipped `relative.py` `E=kQ_0` route is a separate,
    always-finite special case).
    (e) **Char scope:** module/`is_isomorphic` legs over `GF(p)` need char 0 or char > dim (loud
    otherwise); certificate legs run over `QQ` with a `GF(p)` HH parity spot-check.
    (f) **Citation correction (binding):** the Jacobi–Zariski sequence is **CLMS arXiv:2009.05017**
    (Bull. LMS 54 (2022), no. 3) + the classical origin **Kaygun arXiv:1103.4377** (HHA 14 (2012)) —
    NOT `arXiv:1908.11130` (the CLMS *split*-extension paper). The record's arXiv id is corrected
    (fold back into R7 with a dated note at merge). The card says only "1908.11130" (no "Kaygun");
    an earlier draft of THIS plan wrongly attributed "Kaygun" to the card — that plan-writer error is
    corrected in the reference section.
    (g) **`arXiv:2409.00945` authorship VERIFIED** (Wang–Xu–Zhang–Zhou, 2024) — cited as context.
    (h) **Input model narrowing (M1):** general only **among extensions by arrows and relations
    already presented on `Q_A`**; a `B` not generated by a subset of `A`'s arrows is inexpressible
    (loud scope-gate row), NOT "fully general".
  - **Recount the class table** (`tests/release/test_oracle_classes.py` drives the numbers — paste
    LIVE counts, re-run to green; no guessed at-authoring number).
- [ ] **Step 3: README.** One features line: "Han's conjecture across extensions — the bounded-
  extension recognizer (`A/B` tensor-nilpotent as an honest capped semi-decision + `J`-interrupter
  certificate, one-sided `B`-projectivity, finite `pd_{B^e}` free when `gl.dim B < ∞`), the Han
  transport as an **injection/isomorphism ladder** (`B ⊨ Han ⇔ A ⊨ Han` under a bounded extension),
  and the Jacobi–Zariski nearly-exact sequence — R7 (CLMS / Kaygun)."
- [ ] **Step 4: Full gate:** `... tests/invariants tests/families tests/hochschild -q -m deep`,
  `... tests/webapp tests/gui tests/hpc -q -m fast`, `... tests/qpa -q -m qpa`,
  `... tests/release tests/citations -q` — all green.
- [ ] **Step 5:** commit
  `docs(verification): P73 R7 oracle rows + honest scope (capped tensor-nilpotency, injection/iso ladder, gl.dim-primary pd leg, JZ citation correction 2009.05017 + Kaygun 1103.4377, input-model narrowing) + citations + recounted classes`.

---

## Acceptance (Plan-73 definition of done)

1. **Input model, `quiverlab.families.extension`:** `arrow_removal_subalgebra(A, new_arrows)`
   builds `B ⊆ A` with `dim B` certified via `present_from_pi` (the extension axiom `J ∩ B = 0`);
   `Extension.relative_paths()` = the `A/B` basis. **Pins:** Ex. 5.3 `(dim A, dim B, dim A/B) =
   (47, 20, 27)`; Ex. 5.5 `(10, 6, 4)` with `A/B` basis `{d, c*d, d*c, c*d*c}`.
2. **Tensor-nilpotency recognizer, `invariants.han.is_tensor_nilpotent`:** capped semi-decision +
   `J`-interrupter certificate + Ex-5.5 witness; `"undecided"` a first-class verdict. **Pins:**
   Ex. 5.3 `"nilpotent"` index 2; Ex. 5.5 `"not_nilpotent"` with witness; the certificate route's
   vanishing is directly recomputed (`(A/B)^{⊗ 2} = 0`).
3. **One-sided projectivity + finite `pd_{B^e}`:** `one_sided_projective` (`pd_B(A/B)=0` + Thm 5.20
   certificate); finite `pd_{B^e}` via the **`gl.dim B < ∞` PRIMARY route (CLMS Ex. 6.1, free —
   live-verified `gl.dim B = 2` on both examples)**, with `enveloping_algebra(B)` (`dim (dim B)^2`
   certified) as the capped `gl.dim B = ∞` fallback (`int` finite or honest `"undecided"` + lower
   bound, never `∞` unproven).
4. **Jacobi–Zariski tool, `hochschild.jacobi_zariski`:** `relative_homology` = `HH_*(A|B)` via the
   CLMS finite complex (`= 0` for `* ≥ n` on Ex. 5.3); `jacobi_zariski_sequence` "exact twice in
   three" (exact at the `HH_*(B)` and `HH_*(A|B)` spots).
5. **Han transport, `invariants.han.han_transport` — the injection/iso LADDER (H1):** the certificate
   labels its claim by the exact row — `"bounded"` (`HH_*(B) ≅ HH_*(A)`, full `iff`) on Ex. 5.3,
   `"pd_injection"` (i+ii, ordinary injection), `"nilpotent_injection"` (i, mixed-coeff injection),
   `"not_bounded"` on Ex. 5.5, `"undecided"` on a cap. **There is NO `≅`-from-leg-(i) verdict.**
   **Self-cert gate:** the **injection bound `dim HH_m(B) ≤ dim HH_m(A)`** on the shipped engines,
   with **equality asserted only under `"bounded"`** — an equality on Ex. 5.3 (`[5,0,…]=[5,0,…]`),
   differing on the not-bounded Ex. 5.5 (`[4,1,2,3,…] ≠ [3,0,…]`); the `transport_verdict` ladder
   logic is `oracle_selfcert`-verified independently of any algebra. `injection_from` is a certified
   lower bound (H3), not a precise threshold. **Ex. 5.5 wording (M3): both `A` and `B` satisfy Han's
   conjecture individually** (`B` finite gl.dim; `A` vacuously — `HH_*` does not vanish); what fails
   is that the extension is not bounded, so the transport does not apply.
6. **One certificate kind `han_transport`** clickable end-to-end (GUI canvas → arrow-subset picker →
   block → report) in **all four locales** with `pick.kind.*` parity, both runners byte-identical,
   the `new_arrows` request field guarded + canonical-key-stable, ONE golden added (others
   byte-identical); the three-leg table + transport verdict + JZ dims render.
7. **QPA battery green (`-m qpa`):** the `NamesGVars()` guard + the component HH / one-sided-`pd`
   input anchors on the example algebras.
8. `docs/verification.md` recounted (live numbers) with the honest-scope entries (a)–(h) — including
   the **injection/iso ladder** (no `≅`-from-leg-(i)), the **`injection_from` lower-bound** caveat,
   the **strict-injection unexemplified** note (Task III4), and the **input-model narrowing** (M1);
   citations added and BibTeX-verified — **JZ = arXiv:2009.05017 (Bull. LMS 54 (2022)) + Kaygun
   arXiv:1103.4377** (never 1908.11130), CLMS-bounded DOI `10.1016/j.jalgebra.2022.01.022`,
   2409.00945 authorship verified; README line; deep + fast + qpa + release + citations green. No
   dependency on other Wave-4 plans (independent merge to `dev`); the R7 JZ-arXiv-id correction
   folded back into the research doc with a dated note at merge.

---

## Methodology & assumptions

**Approach.** I (1) read the metaplan P73 card + §§1–4 conventions and the R7 record verbatim, and
the P72 card to fix the seam; (2) fetched and read the four load-bearing arXiv abstracts +
Han's-conjecture landscape, and read the CLMS 2101.02597 PDF (pages 1–14: Def. of bounded
extension, Thm 2.2 relative bar complex, Def. 2.3/2.5, Cor. 2.4, Thm 4.2/4.3/4.6, §5 relative
paths/cycles/`J`-interrupters Def. 5.9–5.20, and **Examples 5.3 and 5.5 verbatim**); (3) inventoried
the live substrate by three parallel Explore agents + direct `inspect` introspection — the P52
relative-HH / `Bimodule` surface, the module `pd`/`Tor`/`Ext` stack, the constructors, and the
kinds-wiring touchpoints; (4) **live-built both example algebras and computed the HH sequences +
arrow-removal dims in the venv**; (5) mirrored the house form from plan-65 (records verbatim,
reference re-verification with `# PIN`, live-facts table, scope gates, Plan-32 markers, three-tier
GUI, verification rows, TDD task list).

**Assumptions made.**
- **The venv resolves to the dev-tip src.** My worktree is pinned at `b7fa566` (v0.3.0), BEHIND
  `dev` (`d617108`, which carries P51–P62). The `.venv/bin/python` editable install imports
  `quiverlab` from the **main-repo** `src` (`__file__` confirmed) — i.e. the **dev tip**, the true
  substrate P73 builds on (with P52's relative-HH/`Bimodule` surface present). I deliberately did
  NOT set `PYTHONPATH` to my stale worktree src (it would test v0.3.0, missing P52). The plan file is
  written into my worktree (the deliverable); all live-verification used the dev-tip import.
- **Composition translation.** quiverlab composes left-to-right; CLMS right-to-left. I translated
  every transcribed path and **confirmed the builds are admissible and finite-dimensional** (the
  Groebner engine accepted the inhomogeneous `abcd − αβ`; dims 47/20 and 10/6 are stable), and that
  the `a`/`d`-free basis counts equal `dim B` (the arrow-removal subalgebra is genuine).
- **HH over `GF(32003)` stands in for the field-general statement.** Han's conjecture is a
  field-level statement; the illustration uses `GF(p)` (the fast engine) because the bar oracle
  blows up. The plan's batteries run the certificate legs over `QQ` (exact) with a `GF(p)` HH
  parity spot-check — the design does not depend on the illustration field.
- **`present_from_pi` can present the subalgebra and the enveloping algebra.** I verified it is the
  shipped, reused presenter (`_present.py`), but did NOT run it in subalgebra mode — Task 0 / Task
  II2 are its first such uses; if it strictly requires a known `dim_expected`, the tasks add the
  count-then-present two-pass (flagged).

**Deliberately NOT checked, and why.**
- **I ran no test suite and built no `src/` code** (this is a plan, and other agents run tests
  concurrently — I kept probes light: two algebra builds + HH to degree 9 + `gl.dim B` + a basis
  count). The tensor-power / relative-cycle / `J`-interrupter machinery, the `enveloping_algebra(B)`
  product quiver (fallback only), and `HH_*(A|B)` via the CLMS complex are the plan's implementation
  work — grounded by the hand-checkable Ex. 5.3 (index 2) / Ex. 5.5 (witness `d*c`) combinatorics and
  the live HH pins, but not themselves executed.
- **`pd_{B^e}(A/B)` was not computed directly**, but is **certified finite FREE** by the
  live-verified `gl.dim B = 2` on both examples (CLMS Ex. 6.1) — the enveloping-algebra fallback is
  not exercised by the primary oracles, only by a future `gl.dim B = ∞` input.
- **The exact J. Algebra volume/pages of CLMS 2101.02597** — the DOI is verified
  (`10.1016/j.jalgebra.2022.01.022`, exposed on the arXiv abstract page); only the volume/pages are
  `# PIN`'d (the journal ScienceDirect *article* page 403'd, but the DOI+arXiv suffice; the earlier
  "ScienceDirect 403 ⟹ no DOI" reason was superseded this fix round).
- **A. Kaygun's classical JZ paper is now VERIFIED and cited** — arXiv:1103.4377, Homology Homotopy
  Appl. 14 (2012), no. 1, 65–78 (the noncommutative-JZ origin CLMS credit); this fix round added it
  as `kaygun_jacobi_zariski`. The (distinct) later Erratum was not chased (not load-bearing).
- **The strict-injection example** (tensor-nilpotent-but-unbounded) was NOT constructed — deferred to
  Task III4 with an honest "construct or document impossibility" contract (the ladder LOGIC is
  self-cert-verified regardless).

**Why I believe the result is correct.** The two HARD oracles are independently grounded and
live-checked: Ex. 5.3 (bounded) shows `HH_*(A) = HH_*(B)` (`[5,0,0,…]`) — consistent with the
`"bounded"` ⟹ isomorphism row — and Ex. 5.5 (not bounded) shows them **diverging**
(`[4,1,2,3,3,3,4,6,8,9] ≠ [3,0,0,…]`), computed on the shipped engines. **The headline is the
injection/iso LADDER, honestly read from CLMS Thm 3.1** (leg (i) ⟹ mixed-coeff injection; (i)+(ii)
⟹ ordinary injection; bounded ⟹ iso — never "≅ from nilpotency alone"); the self-cert gate is the
injection bound `dim HH_m(B) ≤ dim HH_m(A)`, equality only under `"bounded"`. leg (ii) is **free when
`gl.dim B < ∞`** (CLMS Ex. 6.1; live-verified `gl.dim B = 2`), so `"bounded"` is reachable on the
common finite-gl.dim case without the enveloping algebra — the heaviest substrate is a fallback, not
the headline path. The recognizer's inputs (relative paths) are read directly off `A.basis_labels`
(live-verified `A/B` basis for Ex. 5.5, and the `a`/`d`-free counts `= dim B`), so leg (i)'s
combinatorics have no hidden engine. Every scope boundary is a loud verdict / honest `"undecided"`.
Finally, the mandatory 2409.00945 re-check is resolved (authorship verified, citable) and the JZ
arXiv-id error is corrected (2009.05017 + Kaygun 1103.4377, never 1908.11130) — both binding on the
citation task.

## Open design risks (what a critic will likely attack)

1. **`arrow_removal_subalgebra` correctness (`I_B = ker π`, not "drop `F`-relations").** An
   `F`-relation can induce a collapse among `F`-free paths; the honest constructor must extract
   `ker π` and certify `dim B`. *Mitigation:* the `dim B = #(F-free lin-indep A-basis)` certificate
   (live-verified 20 / 6) catches a wrong `I_B` loudly; Task 0 pins it.
2. **The finite-`pd_{B^e}` leg (enveloping algebra) — now a FALLBACK (B1).** The primary route is
   `gl.dim B < ∞` (CLMS Ex. 6.1, free — live-verified `gl.dim B = 2` on both examples); only
   `gl.dim B = ∞` needs `B^e` as a product quiver + `A/B` as a `B^e`-module, the heaviest substrate
   (`dim (dim B)^2` grows fast — Ex. 5.3: `400`). *Mitigation:* the primary oracles never touch it;
   the fallback is honestly capped (`"undecided"` + lower bound). A critic may still push to descope
   the fallback to the combinatorial Thm-5.20-style conditions — acceptable if the product quiver
   proves too heavy for v1.
3. **`HH_*(A|B)` via the CLMS complex.** The relative tensor `X ⊗_{B^e} (A/B)^{⊗_B m}` and the
   `σ`-independent differential `b` are fiddly. *Mitigation:* the finite-complex vanishing
   (`= 0` for `* ≥ n`, Cor. 2.4) + the JZ "exact twice in three" self-cert + the **injection-bound**
   transport witness `dim HH_m(B) ≤ dim HH_m(A)` (equality on the bounded Ex. 5.3, live-verified)
   triangulate a wrong differential loudly.
4. **The capped tensor-nilpotency semi-decision vs the certificate.** The two routes must agree
   where both decide; a wrong `J`-interrupter predicate would pass the certificate but fail the
   direct recompute. *Mitigation:* `test_certificate_forces_vanishing` recomputes `(A/B)^{⊗ n}`
   directly; the direct/certificate agreement is a self-cert.
5. **Citations.** The JZ arXiv-id error (record) and the 2409.00945 flag are both resolved here; a
   critic will check the bib ships **2009.05017 + Kaygun 1103.4377** for JZ (never 1908.11130) and
   does not fabricate the J. Algebra 2101.02597 volume (DOI shipped, volume/pages `# PIN`'d).
   *Mitigation:* the reference re-verification section + `# PIN`s + the `bibtex()`/registry gate.
7. **The injection/iso ladder (H1) must be read exactly.** A critic will check the plan never claims
   `HH_*(A) ≅ HH_*(B)` from tensor-nilpotency alone (only an injection), that the mixed-coefficient
   row is labelled as such, and that the self-cert gate is an inequality with equality only under
   `"bounded"`. *Mitigation:* the ladder table, the `transport_verdict` `oracle_selfcert` logic test,
   and the injection-bound gate all encode this; Ex. 5.5's HH divergence is explained as
   not-bounded (both algebras satisfy Han individually — M3), never "Han fails".
6. **GUI new-arrow field.** `han_transport` needs data beyond the algebra (the subset `F`); the
   canonical key must stay byte-stable for every non-`han_transport` request. *Mitigation:* the
   `new_arrows` field defaults `None` and is dropped from the canonical dump when absent (the
   Plan-26 `module`-block precedent); a golden pins an existing request's key unchanged.

## Change log

- **2026-08-08 authoring.** Initial plan (R7: bounded-extension recognizer + Han transport + the
  Jacobi–Zariski nearly-exact sequence). CLMS 2101.02597 read (Examples 5.3/5.5 verbatim, Def.
  bounded, Thm 2.2/3.1/4.6, §5 relative-path machinery). References re-verified: **arXiv:2409.00945
  authorship VERIFIED** (Wang–Xu–Zhang–Zhou — the mandatory check cleared, citable); **the
  Jacobi–Zariski sequence corrected to CLMS arXiv:2009.05017** (Bull. LMS, "exact twice in three")
  — the record/card's "1908.11130 / Kaygun" is a double error (1908.11130 is CLMS's *split*
  predecessor); Han's conjecture pinned to J. LMS 73 (2006). Both example algebras **live-built and
  HH-computed** (Ex. 5.3 bounded: `HH_*(A)=HH_*(B)=[5,0,…]`; Ex. 5.5 not bounded: `HH_*(A)=
  [4,1,2,3,3,3,4,6,8,9] ≠ HH_*(B)=[3,0,…]`); the arrow-removal subalgebra dims certified
  (`a`/`d`-free basis `= dim B`). Substrate design findings recorded: no subalgebra representation
  (Task 0 builds it on `present_from_pi`), no first-class enveloping algebra (Task II1 builds `B^e`
  as a product quiver), relative HH over a general `B` refused (Task III1 delivers it via the CLMS
  finite complex). The finite-`pd_{B^e}` leg flagged as the heaviest, honestly capped; the headline
  (`homology_agreement` via the JZ consequence) rests only on the tensor-nilpotency leg + the
  shipped HH engines.
- **2026-08-08 fix round (critic NEEDS WORK → all findings adjudicated VALID).** Applied,
  document-only; new numbers live-verified in the venv (dev-tip import) where cheap:
  - **H1 (BLOCKING)** — the headline was corrected from a false "`HH_*(A) ≅ HH_*(B)` from leg (i)
    alone" to the exact CLMS Thm 3.1 **injection/isomorphism ladder**: leg (i) ⟹ `H_*(B,A) ↪
    H_*(A,A)` (mixed-coeff injection); (i)+(ii) ⟹ `HH_*(B) ↪ HH_*(A)` (ordinary injection); bounded
    ⟹ `HH_*(B) ≅ HH_*(A)` (iso). The `"homology_agreement"` verdict is gone; the verdict enum is
    `"bounded" | "pd_injection" | "nilpotent_injection" | "not_bounded" | "undecided"`; the self-cert
    gate is the injection bound `dim HH_m(B) ≤ dim HH_m(A)`, equality only under `"bounded"`. Goal,
    Math foundation (transport consequence + verdict table), API dataclass, GUI block, tests,
    Acceptance, verification honest-scope all rewritten.
  - **H2** — the false "HH agree when tensor-nilpotent" gate removed; **Task III4 added** to
    construct a tensor-nilpotent-but-unbounded example exercising the strict injection (or document
    impossibility) — unexemplified at spec time, honest note.
  - **B1** — `pd_{B^e}(A/B)` leg reordered: **`gl.dim B < ∞` PRIMARY (CLMS Ex. 6.1, free)** — live-
    verified `gl.dim B = 2` (exact) for BOTH Ex. 5.3 and Ex. 5.5; `enveloping_algebra(B)` demoted to
    the `gl.dim B = ∞` fallback (Task II1 = gl.dim, Task II2 = enveloping, Task II3 = combined).
  - **W1** — integrity fix: the metaplan card says only "1908.11130" (NOT "Kaygun"); the earlier
    draft's "Kaygun 1908.11130 (card)" attribution was this plan-writer's error, corrected. **Kaygun
    IS now cited properly** — arXiv:1103.4377, HHA 14 (2012) 65–78 (the noncommutative-JZ origin CLMS
    credit), added as `kaygun_jacobi_zariski`.
  - **M1** — input-model overclaim removed: "fully general for f.d. inputs" → "general **among
    extensions by arrows and relations presented on `Q_A`**"; a loud scope-gate row added disclosing
    the narrowing (a `B` not on a sub-quiver of `Q_A` is inexpressible).
  - **W3** — (a) the cap-only test asserts EXACTLY `"undecided"`; (b) `test_injection_ladder_logic`
    + `transport_verdict` added (ladder logic, algebra-independent); (c) the JZ-exactness test notes
    it runs on trivial-higher-HH Ex. 5.3 and defers the nontrivial-higher-HH case to Task III4
    (Ex. 5.5 cannot serve — not tensor-nilpotent ⟹ infinite relative complex).
  - **H3** — `hh_agreement_from` renamed `injection_from`, respec'd as a CERTIFIED LOWER BOUND (the
    theorem gives `* ≫ 0` / bounded-exactness degree, not the nilpotency index).
  - **M3** — Ex. 5.5 wording fixed everywhere: both `A` and `B` satisfy Han individually; what fails
    is that the extension is not bounded (the transport does not apply), not "Han transport fails".
  - **W2** — `CLMSbounded2022` DOI filled (`10.1016/j.jalgebra.2022.01.022`, verified on the arXiv
    page); `CLMSjacobiZariski` year fixed to **2022** (Bull. LMS **54**, no. **3**); the false
    "ScienceDirect 403 ⟹ no DOI" methodology reason superseded.
  - **P72 seam** — the arrow-removal constructor restated as **shared substrate owned by whichever of
    P72/P73 lands first**, with an explicit consume-don't-reimplement contract.
