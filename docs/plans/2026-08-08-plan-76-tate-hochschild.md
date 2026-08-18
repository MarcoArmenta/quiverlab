# Plan 76: Tate–Hochschild (singular Hochschild) cohomology — positive AND negative degrees with cup product, via complete resolutions (P76 / R3)

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking. Do the reference re-verification in
> the Record/Reference sections BEFORE writing any oracle pin — the attributions below
> (`# PIN`) must be resolved against the sources, not this plan's prose. The four
> references were fetched and read at authoring (see Reference re-verification); the
> Bergh–Jorgensen theorems below are transcribed **verbatim from the PDF**.

**Goal.** Promote **Tate–Hochschild cohomology** `ĤH^*(A)` — a.k.a. **singular Hochschild
cohomology** `HH_sg^*(A)` — to a first-class, no-code surface: the graded ring
`ĤH^n(A) = Êxt^n_{A^e}(A, A)` defined for **all** `n ∈ ℤ` (**negative as well as
positive**) via a **complete resolution** `𝕋` of `A` over its enveloping algebra `A^e`,
built by **splicing the minimal `A^e`-projective resolution of `A` with its
(Nakayama-twisted) `A^e`-dual**. Above a stated degree the positive part **coincides with
ordinary `HH^*`** (Bergh–Jorgensen threshold); the negative part is new. The surface ships
with the **cup product** on the full `ℤ`-graded ring, an **eventual-periodicity certificate**
(the existence of an invertible homogeneous element — Usui — realised through the shipped
`Ω`-periodicity certificate of `modules/homdims.py`), and **honest scope**: the native ring
serves **self-injective** algebras `⇒` always (the dual-splice is well-defined there); not
(Iwanaga-)Gorenstein `⇒` LOUD refusal; **Gorenstein-non-self-injective (even eventually
periodic) `⇒` DEFERRED** (the `D(P_n)`-dual splice cannot build it — a periodicity-extension
is the named future task). One range compute kind `tate_hochschild` exposes the **positive AND
negative degree tables** end-to-end in all four locales.

**What is genuinely NEW over the shipped stack (the scope boundary — read this first).**
quiverlab already ships *ordinary* Hochschild cohomology (`Algebra.hochschild_cohomology`,
the bar/fast/CS engines) and the *nonnegative-degree* product surface (Plan 35 —
`Algebra.cup_products`). Plan 76 adds four things the stack does not have:

1. **Negative degrees.** `ĤH^{-1}, ĤH^{-2}, …` — the "stable" cohomology below degree 0,
   which ordinary `HH^*` (living only in degrees `≥ 0`) cannot express. This is the headline
   object and requires the **complete resolution**, new machinery.
2. **The complete-resolution splice.** The minimal `A^e`-resolution `P_•` (shipped
   `engine/resolutions_minimal.py`) is spliced with its Nakayama-twisted `A^e`-dual to give
   an acyclic two-sided complex `𝕋: … → T_1 → T_0 → T_{-1} → T_{-2} → …`. The dual half's
   corner-typed terms are the **swapped/`ν`-permuted corner tags** of the positive half (the
   Plan-16 note, made two-sided). `d∘d = 0` across the splice joint is a **mandatory
   self-cert**.
3. **`ĤH^0 ≠ HH^0` and the graded-`ℤ` ring.** The Tate `ĤH^0` is the *stable* centre (a
   quotient of `Z(A) = HH^0`), not `HH^0` itself; the cup product extends to the full
   `ℤ`-graded ring, where the eventual-periodicity generator is **invertible** (lives in a
   negative degree too).
4. **The eventual-periodicity certificate.** `Ω^{n+p}_{A^e}(A) ≅ Ω^n_{A^e}(A)` (the shipped
   `omega_periodicity` pattern, `is_isomorphic`-certified) ⟺ (for Gorenstein `A`, Usui) an
   invertible homogeneous element of `ĤH^*(A)` — the two are cross-checked.

**Architecture.** One new engine module + one new public surface module + one range compute
kind; everything else is a thin exact layer over primitives already on `dev`:

- **`src/quiverlab/engine/complete_resolution.py`** (new — the GF(p) splice engine, over
  int64 `GF(p)`, the minimal-`A^e` sibling). Public to the library (not to app code):
  - `complete_resolution(A, N, p) -> CompleteResolution` — build `T_n` for `-N-1 ≤ n ≤ N`:
    `T_n = P_n` (the minimal `A^e`-resolution, `n ≥ 0`); `T_{-n-1} = D_{A^e}(P_n)` (the
    Nakayama-twisted `A^e`-dual, corner tags `(v,w) ↦ (π w, π v)` where `π` is the VERTEX
    Nakayama permutation `soc P_v ≅ S_{π(v)}`, NOT the `d×d` algebra matrix — M3, `n ≥ 0`);
    the joint differential `T_0 → T_{-1}` is `(P_0 ↠ A ↪ D(P_0))`. **`assert_dd_zero`** across
    every joint (the mandatory self-cert), **`assert_acyclic`** in the served window. SELF-
    INJECTIVE input only in v1 (M4 — else the dual half is not projective; loud DEFERRED
    refusal).
  - `tate_cohomology_dims(A, N, p) -> dict[int,int]` — `dim ĤH^m = dim H^m(Hom_{A^e}(𝕋, A))`
    for `-N ≤ m ≤ N`; **positive-agreement self-cert** `ĤH^m ≡ HH^m` for `m ≥ agrees_from`.
  - `tate_homology_dims(A, N, p) -> dict[int,int]` — `dim ĤH_m = dim H_m(A ⊗_{A^e} 𝕋)`.
- **`src/quiverlab/hochschild/tate.py`** (new — the public surface; sibling of
  `hochschild/products.py`). Public via `import quiverlab`:
  - `tate_hochschild(A, top, *, engine="auto", max_cells=…) -> TateHochschild` — the frozen
    result object (`pos_dims`, `neg_dims`, `hat_hh0`, `agrees_from`, `period`,
    `periodicity_degree`, `scope`, `cup`, `references`, `note`); routes GF(p)-native /
    duality-dims / positive-threshold (see Mathematical foundation).
  - `tate_cup_products(A, top, *, engine="auto") -> HHProducts` — the structure-constant
    tables of the Tate cup `ĤH^p ⊗ ĤH^q → ĤH^{p+q}` on the recorded basis, reusing the
    Plan-35 `HHProducts` container; the positive part cross-checks `A.cup_products`.
  - `tate_periodicity(A, *, max_period=…, p=…) -> TatePeriodicity` — the eventual-periodicity
    certificate (Ω-periodicity of `A` over `A^e`, `is_isomorphic`-certified — the shipped
    `modules/homdims.py::omega_periodicity` pattern — cross-checked against Usui's invertible
    homogeneous element via the Tate cup).
  - `tate_hochschild_block(A, top, *, engine="auto") -> dict` — the JSON block for both
    runners.
- **`Algebra.tate_hochschild(top, engine="auto")`** — thin lazy-import delegate beside
  `Algebra.cup_products` (`core/algebra.py`).
- The `tate_hochschild` **range kind** (`tate_hochschild:0..N` — `N` = the positive top, the
  block reports the symmetric window `[-N, N]`) is wired into `hpc/spec.py::_dispatch` +
  `docs/gui/runner.py` (byte-identical twins) + `webapp/server/schema.py` (the THIRD
  grammar-parse site), the GUI touchpoints (`qlgui-tate` checkbox + `-top` input,
  `scheduleProbe`, `renderBlock`), the i18n chains, and `trace/results_html.py` — following
  the `hh_cohomology`/`cup` range-kind precedent.

**Tech Stack.** Exact `GF(p)` int64 linear algebra for the complete-resolution splice (the
minimal-`A^e` engine's `AeEngine`, `nullspace_mod_p`, `reduce_mod_nullspace`); exact `sympy`
/ `Domain` linear algebra for the duality-dims and positive-threshold routes over any exact
`Domain`. **No floats in `src/`** (AST-gated by `tests/test_no_floats.py`): dims are `int`,
cup constants are exact strings (Plan-35 `ProductTable`), verdicts/periods are `int`/`bool`/
`str`. Composition **left-to-right** (`a*b` = first `a` then `b`). All refusals are
`QuiverlabError`.

---

## Record (verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`, R3)

> **R3 — Tate–Hochschild (singular HH) for self-injective / eventually periodic Gorenstein.**
> [D-scout P3; keep] Object: ĤH^*(A) in negative and positive degrees with cup product, via
> complete resolutions (splice minimal A^e with its dual); eventual-periodicity certificate
> (invertible homogeneous element). Refs: Keller arXiv:1809.05121; Usui arXiv:2107.03326;
> Bergh–Jorgensen arXiv:1109.4019. Scope: eventually-periodic Gorenstein only (the class
> where one period makes it finite) — loud refusal otherwise. Oracles: positive degrees ≡
> ordinary HH; k[x]/(x^n) periodic Tate ring; symmetric Nakayama zoo periods; Tate ≡ ordinary
> above findim. Size M–L.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P76 (Wave 4, tier δ,
**independent**). Size M–L.

---

## Reference re-verification (done at authoring — the four papers were fetched and read; findings binding, `# PIN` = resolve at citation-add)

The metaplan standing rule (§1.5) requires the plan writer to re-verify its record's
citations. **All four references were fetched at authoring** (arXiv abstracts; the
Bergh–Jorgensen PDF read page-by-page). Findings:

1. **The DEFINITIONAL origin is Zhengfang Wang, NOT Keller (card correction, evidence).**
   The card names only "Keller arXiv:1809.05121". Verified: **arXiv:1508.00190 — Zhengfang
   Wang, "Singular Hochschild Cohomology and Gerstenhaber Algebra Structure"** is the paper
   that *introduces* singular Hochschild cohomology `HH_sg^i(A,A) = Hom_{D_sg(A⊗A^op)}(A, A[i])`
   for `i ∈ ℤ`, endows it with a Gerstenhaber (and, for symmetric `A`, BV) structure, and it
   is the origin of the theory this plan computes. **Add it as the primary definitional
   citation `wang_singular_hh`.** `# PIN`: journal ref (Adv. Math., verify at citation-add)
   or ship as `@misc{…, note={arXiv:1508.00190}}`.

2. **arXiv:1809.05121 — Bernhard Keller, "Singular Hochschild cohomology via the singularity
   category."** The card's "Keller-related … (VERIFY what this actually is)" is **RESOLVED**:
   it IS by Keller and it IS about Tate–Hochschild cohomology — it proves that for a
   noetherian algebra whose bounded dg derived category is smooth, the singular Hochschild
   cohomology `(= Tate–Hochschild cohomology)` is isomorphic as a graded algebra to the
   Hochschild cohomology of the dg singularity category. It is the **singularity-category
   identification** (building on Wang), not the computational recipe. C. R. Acad. Sci. Paris
   356 (2018), no. 11-12, 1106–1111. Key `keller_singular_hh`. `# PIN`: page range
   (1106–1111) at citation-add.

3. **arXiv:2107.03326 — Satoshi Usui, "Tate-Hochschild cohomology rings for eventually
   periodic Gorenstein algebras."** VERIFIED. Establishes: (i) Tate–Hochschild cohomology is
   defined on positive and negative degrees with a ring structure; (ii) **eventually periodic
   algebras are NOT necessarily Gorenstein**; (iii) for a **Gorenstein** algebra, eventual
   periodicity is characterised as **the existence of an invertible homogeneous element of
   the Tate–Hochschild cohomology ring** — *exactly* the card's "eventual-periodicity
   certificate (invertible homogeneous element)"; (iv) a construction of eventually periodic
   Gorenstein algebras via tensor algebras. SUT J. Math. Key `usui_tate_periodic`. `# PIN`:
   volume/pages of the SUT J. Math. published version at citation-add (ship
   `note={arXiv:2107.03326}` otherwise — never guess).

4. **arXiv:1109.4019 — Petter Andreas Bergh & David A. Jorgensen, "Tate-Hochschild homology
   and cohomology of Frobenius algebras," J. Noncommut. Geom. 7 (2013), no. 4, 907–937.**
   VERIFIED (PDF read). This is the **computational reference** — the definition via a
   complete resolution over `A^e`, the threshold, the duality, and an explicit computation.
   The load-bearing statements transcribed **verbatim** (see Mathematical foundation). Key
   `bergh_jorgensen_tate`. `# PIN`: DOI `10.4171/JNCG/139` and pages 907–937 at citation-add.

**Two card wordings corrected with evidence (both fold back into the research doc with a
dated note):**

- **"Tate ≡ ordinary above findim" is imprecise.** The exact Bergh–Jorgensen statement (PDF
  p. 4, verbatim): *"if the Gorenstein dimension of the enveloping algebra is `d`, then for
  every `n ≥ d+1` there are isomorphisms `ĤH^n(Λ,B) ≅ Ext^n_{Λ^e}(Λ,B)`"*, and
  `Ext^n_{Λ^e}(Λ,B) = HH^n(Λ,B)`. So the threshold is **`n ≥ d+1` with `d = the
  (Iwanaga-)Gorenstein dimension of the ENVELOPING algebra `A^e``**, not "findim of `A`".
  For self-injective `A` (the primary scope) `A^e` is self-injective, `d = 0`, so `ĤH^n ≡
  HH^n` for `n ≥ 1`. The plan pins the precise threshold and names `agrees_from = d+1`.

- **The scope "eventually-periodic Gorenstein" is refined with a live witness (attribution
  corrected, MINOR a).** The metaplan §5 R3 card says exactly *"scope = eventually-periodic
  Gorenstein only … loud refusal otherwise"* — it does **not** itself write "self-injective ⇒
  always". That reading is **my (plan-writer's) expansion**, matching the coordinator's scope
  guidance, and the live probes confirm it is the correct one, NOT "eventually periodic only":
  **`QuantumCI(2,2,2)` is self-injective but has `A^e`-Betti ranks `[1,2,3,4,…]` (complexity 2
  — NOT eventually periodic)**, yet Bergh–Jorgensen explicitly compute its Tate–Hochschild
  cohomology. So self-injectivity (`⇒ A^e` self-injective, Gorenstein dimension 0 `⇒` the
  complete resolution exists, and `D_{A^e}(P_n)` is projective so the dual-splice is
  well-defined) is what guarantees service in v1. **Eventual periodicity is NOT required for
  the self-injective case; it is the property that would let a FUTURE (deferred, M4)
  periodicity-extension construction serve the general Gorenstein NON-self-injective case**
  (where the dual-splice fails). Usui's "eventually periodic algebras need not be Gorenstein"
  is the dual caution: an eventually periodic non-Gorenstein algebra has no complete resolution
  over `A^e`, so `ĤH^*` via complete resolutions is **not defined** — LOUD refusal. **v1 native
  scope = self-injective; not-Gorenstein = loud refusal; Gorenstein-non-self-injective =
  DEFERRED.**

**New citation keys to add:** `wang_singular_hh` (arXiv:1508.00190), `keller_singular_hh`
(arXiv:1809.05121), `usui_tate_periodic` (arXiv:2107.03326), `bergh_jorgensen_tate`
(arXiv:1109.4019 / JNCG 7). **Already shipped and reused:** `quantum_ci` (BGMS2005),
`qci_hh_oracle` (BerghErdmann2008 — the ordinary-HH QCI oracle), `assem_book` (ASS2006),
`cyclic`, `bardzell`, `chouhy_solotar` (engine keys, surfaced when a route uses them).

---

## Mathematical foundation (the definitions, exactly — the plan's ground truth)

Throughout `A = kQ/I` is basic, connected, finite-dimensional over an exact `Domain`, with a
**path presentation** (`A.quiver is not None`); `A^e = A ⊗_k A^op` is the **enveloping
algebra**; `A` is the regular `A^e`-module (the diagonal bimodule); `ν` is the **Nakayama
automorphism** (`ν = id` iff `A` is symmetric).

### Bergh–Jorgensen: the definition (verbatim, PDF pp. 3–4)

For `Λ = A` two-sided Noetherian Gorenstein with `A^e` two-sided Noetherian Gorenstein of
Gorenstein dimension `d`, `A` admits a **complete resolution** `𝕋` over `A^e`:
```
𝕋 :  … → T_2 → T_1 → T_0 → T_{-1} → T_{-2} → …
```
an **acyclic** complex of finitely generated projective `A^e`-modules with (1) the dual
complex `𝕋^* = Hom_{A^e}(𝕋, A^e)` acyclic, and (2) a projective resolution `P_• → A` and a
chain map `𝕋 → P_•` bijective in degrees `n ≥ d` — so **`𝕋` is "eventually" the projective
resolution `P_•`**. Then, for `n ∈ ℤ`,
```
ĤH^n(A, A) = Êxt^n_{A^e}(A, A) = H^n( Hom_{A^e}(𝕋, A) ),
ĤH_n(A, A) = T̂or_n^{A^e}(A, A) = H_n( A ⊗_{A^e} 𝕋 ).
```

### Bergh–Jorgensen: the threshold (verbatim, PDF p. 4)

> Note that if the Gorenstein dimension of the enveloping algebra is `d`, then for every
> `n ≥ d+1` there are isomorphisms `ĤH^n(Λ, B) ≅ Ext^n_{Λ^e}(Λ, B)` and `ĤH_n(Λ, B) ≅
> Tor_n^{Λ^e}(B, Λ)`.

Since `Ext^n_{A^e}(A, A) = HH^n(A)` and `Tor_n^{A^e}(A, A) = HH_n(A)`:
```
ĤH^n(A) = HH^n(A)  and  ĤH_n(A) = HH_n(A)   for all n ≥ d+1,   d = Gor.dim(A^e).
```
`agrees_from := d+1`. **Self-injective ⇒ `A^e` self-injective ⇒ `d = 0` ⇒ agreement for
`n ≥ 1`.** (This is the "positive degrees ≡ ordinary HH" AND the "Tate ≡ ordinary above
findim" oracle — with the *precise* degree.)

### Bergh–Jorgensen: the Frobenius duality (verbatim, PDF p. 2)

> **Theorem.** Let `Λ` be a Frobenius algebra, with Nakayama automorphism `ν`. Then
> `dim_k ĤH^n(Λ, Λ) = dim_k ĤH^{-(n+1)}(Λ, _{ν²}Λ_1)` for all `n ∈ ℤ`, where `_{ν²}Λ_1`
> denotes the bimodule `Λ` twisted on the right by `ν²`.
> Thus Tate-Hochschild cohomology is symmetric when `ν` squares to the identity, and this is
> the case, for example, when `Λ` is a symmetric algebra or an exterior algebra.

For a **symmetric** algebra (`ν = id`, so `ν² = id` and `_{ν²}A_1 = A`):
```
dim ĤH^n(A) = dim ĤH^{-(n+1)}(A)   for all n ∈ ℤ,
      equivalently   dim ĤH^{-j}(A) = dim ĤH^{j-1}(A)   for all j ∈ ℤ.
```
**This confirms the card's `ĤH^{-n} ≅ ĤH^{n-1}` for symmetric algebras exactly** (set
`j = n`). The **homology** version is symmetric for *every* Frobenius algebra (PDF p. 2,
verbatim): `dim ĤH_n(A) = dim ĤH_{-(n+1)}(A)` for all `n ∈ ℤ`.

### Bergh–Jorgensen: the explicit quantum-complete-intersection computation (verbatim, PDF p. 2)

> In Section 4 we compute the Tate-Hochschild cohomology for the quantum complete
> intersection `A = k⟨X,Y⟩/(X^a, XY − qYX, Y^b)` with `a,b ≥ 2` and `q` not a root of unity
> in `k`, finding that
> `dim ĤH^n(A,A) = 1 if n=0, 2 if n=1, 1 if n=2, 0 if n ≠ 0,1,2`.

This is a Frobenius algebra whose Tate cohomology is **NOT symmetric** (`ν²≠id`, since `q` is
not a root of unity) — a discriminating oracle (see Live-verified facts). **Hypothesis note
(binding):** `q` not a root of unity ⇒ **char-0 only** (over `GF(p)` every nonzero `q` is a
root of unity, so this pin is `QQ`-only; quiverlab's `QuantumCI` HH oracle is already char-0,
Plan 33). **Convention note (minor):** quiverlab's `QuantumCI(q, a, b) = k⟨x,y⟩/(x^a, y^b,
xy + q·yx)` (verified in `families/quantum.py`; the established `+q·yx` convention, NOT the
textbook `XY − qYX` of Bergh–Jorgensen). For **generic `q`** (not a root of unity) both are
generic quantum complete intersections with the **same Hochschild (co)homology** (the
`families/quantum.py` docstring states this; they differ only at root-of-unity `q`), so
`QuantumCI(2,2,2)` over `QQ` is a valid instance of the Bergh–Jorgensen family for the HH/Tate
pins — with the sign folded honestly into the relation (`x*y + 2*y*x`, live-verified).

### The splice construction (the NEW machinery — specified precisely)

The plan builds `𝕋` in the window `-N-1 ≤ n ≤ N` from the shipped minimal-`A^e` engine:

1. **Positive half = the minimal `A^e`-resolution `P_•`** (`engine/resolutions_minimal.py`,
   `minimal_resolution(A, N, p)`): `T_n := P_n` for `0 ≤ n ≤ N`, corner-typed
   `P_n = ⊕_g A e_{v_g} ⊗ e_{w_g} A` (generator tags `(v_g, w_g)`), differentials `d_n`.
2. **Negative half = the Nakayama-twisted `A^e`-dual of `P_•`.** `T_{-n-1} := D_{A^e}(P_n)`
   for `0 ≤ n ≤ N`. Because `A^e` is self-injective on the served scope, the `A^e`-dual of a
   corner projective is again a corner projective with the **swapped, `π`-permuted tag**,
   where `π` is the **vertex-level Nakayama permutation** (`soc(P_v) ≅ S_{π(v)}`):
   ```
   D_{A^e}( A e_v ⊗ e_w A )  ≅  A e_{π(w)} ⊗ e_{π(v)} A.        # π = VERTEX permutation, not the algebra matrix
   ```
   This is the two-sided form of the **Plan-16 swapped-tag note** (`minimal_cohomology_dims`
   already acts `a·w·b` on the swapped corner block `e_v A e_w`).
   **`π` is the VERTEX permutation, NOT the `d×d` `nakayama_automorphism()` matrix (M3, live
   caveat).** The shipped `Algebra.nakayama_automorphism()` returns the algebra-level `d×d`
   matrix, which is **not identity even for a symmetric algebra** (live-verified: symmetric
   `kZ₂/J³` has a non-identity `6×6` `ν` matrix). What the tag rule needs is the permutation
   `π` on vertices; it is **identity iff `A` is *weakly symmetric*** (`soc P_v ≅ S_v`), which
   holds for every symmetric `A`. Task I1 adds a `nakayama_permutation(A) -> dict` helper
   extracting `π` from `soc(P_v)` (reusing the `is_selfinjective` socle logic —
   live-verified: `kZ₂/J³` → identity, `kZ₃/J²` → the 3-cycle `{1:2, 2:3, 3:1}`). For a
   **weakly symmetric / symmetric** `A` (`π = id`) the tag is exactly `(w, v)`; for a
   self-injective NON-symmetric `A` (e.g. `kZ₃/J²`) the `π`-permutation is load-bearing (the
   discriminating witness — see below). The negative differentials are the `A^e`-duals
   (corner-transposes) of the positive ones: `d_{-n-1} := D_{A^e}(d_{n+1})` with the
   `π`-permuted tag swap. **`# PIN`:** the exact side/direction of the `π`-twist (whether it
   is `π` or `π^{-1}`, and any Koszul sign) is derived from the projective-injective
   identification `I_v = D(A e_v)` at implementation and **arbitrated by `assert_dd_zero` /
   `assert_acyclic`** on `kZ₃/J²` (a wrong `π` fails them loudly — the P65-style
   convention-by-self-cert), not assumed from this prose.
3. **The splice joint `d_0 : T_0 → T_{-1}`** is the composite
   ```
   T_0 = P_0  --μ-->  A  --ι-->  D_{A^e}(P_0) = T_{-1},
   ```
   where `μ` is the augmentation `P_0 ↠ A` and `ι` is the coaugmentation `A ↪ I^0` into the
   injective envelope `I^0 = D_{A^e}(P_0)` (`= P_0`'s dual, injective since `A^e`
   self-injective). Concretely `ι` is `D_{A^e}` applied to `μ`.
4. **Mandatory self-certs.** `assert_dd_zero`: `d_{n} ∘ d_{n+1} = 0` for **every** joint,
   including the two joint-adjacent products `d_0 ∘ d_1 = 0` (`T_1 → T_0 → T_{-1}`) and
   `d_{-1} ∘ d_0 = 0` (`T_0 → T_{-1} → T_{-2}`) — a wrong splice sign/tag is caught here.
   `assert_acyclic`: `𝕋` is exact in the window `[-N, N]` (rank–nullity per degree). These
   are `oracle_selfcert` and are the primary internal oracle for the NEW machinery.

`ĤH^m` is then `H^m` of `Hom_{A^e}(𝕋, A)` (collapse `a·w·b`, the Plan-16 cohomology side) and
`ĤH_m` is `H_m` of `A ⊗_{A^e} 𝕋` (collapse `b·w·a`, the homology side, swapped tags).

### The three routes (honest dispatch)

- **`native` (GF(p), the primary NEW engine) — SELF-INJECTIVE ONLY in v1 (M4).** Build `𝕋`
  by the splice, compute `ĤH^m` for `-N ≤ m ≤ N` directly. Serves any **self-injective**
  path-presented algebra over `GF(p)` (`A^e` self-injective ⇒ `d = 0` ⇒ the dual `D_{A^e}(P_n)`
  IS projective, so the splice is well-defined) — **periodic or not** (the `QCI` witness is
  self-injective, complexity 2, still served as a bounded window). Ships `ĤH^0`, the negative
  degrees, reps, and the cup. All the `k[x]/(x^n)` (both char cases), symmetric-Nakayama, and
  the `kZ₃/J²` twist-witness pins pin here.
  - **DEFERRED (M4): the general eventually-periodic-Gorenstein (NON-self-injective) case.**
    The described `D_{A^e}(P_n)`-splice **cannot** build a complete resolution when `A` is
    Gorenstein but not self-injective: `D_{A^e}(P_n)` is projective **iff `A^e` is
    self-injective**, so off self-injective the dual half is not a complex of projectives.
    The correct construction there is a **periodicity-extension** complete resolution
    (splice `P_•` with a periodic repeat of a syzygy `Ω^{n+p} ≅ Ω^n`, per Usui), a **future
    task** (`P76-followup: periodicity-extension complete resolution for Gorenstein
    non-self-injective`). v1 native scope is **self-injective**; the eventual-periodicity
    certificate (`tate_periodicity`, Task A2) is still computed for every input (it is a
    property), but the native Tate RING for a non-self-injective Gorenstein algebra is
    deferred, not served.
- **`duality` (any exact `Domain`, symmetric `π = id` i.e. `ν² = id`) — degrees `|m| ≥ 1`,
  degree 0 native-only (M1).** The Bergh–Jorgensen symmetric duality reflects `n ↔ -(n+1)`
  (about `n = -½`), so it links `ĤH^{j-1} ↔ ĤH^{-j}`. Concretely:
  - **positive:** `ĤH^n = HH^n` for `n ≥ agrees_from = 1` (the threshold — the shipped
    ordinary `HH^*`). **Degree 0 is `HH^0`, which is NOT `ĤH^0`** — the positive route does
    not know `ĤH^0`.
  - **negative, `j ≥ 2` ONLY:** `dim ĤH^{-j} = dim ĤH^{j-1} = dim HH^{j-1}` (duality, then
    threshold, since `j-1 ≥ 1`). So the duality route supplies `ĤH^{-2}, ĤH^{-3}, …, ĤH^{-N}`.
  - **`ĤH^0` and `ĤH^{-1}` are a CLOSED 2-cycle** (`dim ĤH^0 = dim ĤH^{-1}`) **unlinked to any
    threshold-known positive** — the duality route CANNOT determine them. They are returned
    **`None` (native-only)**, **never `HH^0`** (the earlier draft's `ĤH^0 = HH^0` was wrong —
    the true value is `n-1 ≠ n = HH^0` for `k[x]/(x^n)`).
  So `duality` serves degrees `m ≥ 1` (`= HH`) and `m ≤ -2` (`= HH^{|m|-1}`), with `m ∈
  {0, -1}` honest-absent (`None`) unless the native GF(p) route runs. **Dims only**
  (theorem-backed, not an independent complete resolution) — `engine="duality"`. Refuses
  loudly when `π ≠ id` (`ν² ≠ id`, e.g. the QCI). Serves `k[x]/(x^n)`, symmetric Nakayama over
  `QQ` (in `|m| ≥ 1`, degree-0 pair absent).
- **`positive` (any exact `Domain`, always available) — degrees `n ≥ agrees_from` ONLY.** The
  Bergh–Jorgensen threshold alone: `ĤH^n = HH^n` for `n ≥ agrees_from` (`= 1` for
  self-injective) — a truthful *partial* answer that never claims `ĤH^0` (`= HH^0` is the
  ordinary value, NOT `ĤH^0` — M2) nor any negative degree. This makes "positive degrees ≡
  ordinary HH" a universal oracle and covers the QCI over `QQ` (`ĤH^{n≥1} = [2,1,0,…]`) even
  though its full ring is neither GF(p)-computable there (root of unity) nor
  duality-mirrorable (non-symmetric).

`engine="auto"` picks `native` over `GF(p)` when `A` is **self-injective** (path-presented),
else `duality` for symmetric (`π = id`) over any Domain, else `positive` with an honest note;
explicit `engine=` forces a route and refuses loudly when its hypothesis fails.

### The eventual-periodicity certificate (reuse the shipped Ω-periodicity — the precedent)

`A` is **eventually periodic (over `A^e`)** iff `Ω^{n+p}_{A^e}(A) ≅ Ω^n_{A^e}(A)` for some
`n ≥ 0`, `p ≥ 1` (`Ω` = syzygy). The certificate reuses `modules/homdims.py::omega_periodicity`
(least `k` with `Ω^k M ≅ M`, `is_isomorphic`-certified, honest `None`):

- **Cheap necessary signal:** `A.complexity(N) == 1` (bounded `A^e`-Betti — the shipped
  minimal-`A^e` growth invariant). Live: `k[x]/(x^n)`, symmetric Nakayama → complexity 1;
  `QCI` → 2; `kA₂` → 0.
- **Rigorous certificate:** `Ω^p_{A^e}(A) ≅ A` (periodic) / `Ω^{n+p} ≅ Ω^n` (eventually
  periodic) via `is_isomorphic` on the `A^e`-modules. `# PIN`: verify that
  `TensorProduct(A, A.opposite())` builds `A^e` as a first-class Algebra and that `A` can be
  presented as the regular `A^e`-module so the shipped `omega_periodicity` applies verbatim;
  if the path-basis guard trips on the tensor build, fall back to the engine-level syzygy
  certificate (compare the minimal-`A^e` differential blocks up to `is_isomorphic` of the
  syzygy corner modules — the same `is_isomorphic` certificate, at the engine layer).
- **Usui cross-check (the ring-theoretic equivalent):** for Gorenstein `A`, eventual
  periodicity ⟺ an invertible homogeneous `η ∈ ĤH^p` (`p ≥ 1`) with `η · η^{-1} = 1 ∈ ĤH^0`
  in the Tate cup ring. The `native` cup exhibits `η` (`periodicity_degree = p`); `η·η^{-1}=1`
  is an `oracle_selfcert` checked on `k[x]/(x^n)` (`p = 2`, see below).

**Multi-vertex gotcha (live-reproduced, binding).** The minimal-`A^e` engine's radical
detection **refuses `A.unit_adapted()` on multi-vertex algebras** (Plan-18 note:
unit-adaptation trips the Plan-13 radical guard) — reproduced live: `NakayamaAlgebra(n=2,l=3,
cyclic=True)` built via `Quiver.algebra` gives `A^e` ranks `[2,2,2,…]` **only when fed
`to_engine(A)` un-adapted**; `to_engine(A.unit_adapted())` raises `radical_basis: … not
nilpotent`. The complete-resolution engine MUST feed the minimal engine the **un-adapted**
algebra for multi-vertex input (single-loop `k[x]/(x^n)` is unaffected). Task I1 pins this.

---

## Live-verified facts (all recomputed in the venv at spec time; shipped engines)

**Ordinary `HH^*` / `HH_*` (0..8), the positive-degree targets and the threshold inputs:**

| algebra | field | `HH^*` = `HH_*` (0..8) | note |
|---|---|---|---|
| `k[x]/(x²)` | `QQ`, `GF(3)` (char∤2) | `[2,1,1,1,1,1,1,1,1]` | char ∤ n |
| `k[x]/(x²)` | `GF(2)` (char\|2) | `[2,2,2,2,2,2,2,2,2]` | char \| n |
| `k[x]/(x³)` | `QQ`, `GF(2)` (char∤3) | `[3,2,2,2,2,2,2,2,2]` | char ∤ n |
| `k[x]/(x³)` | `GF(3)` (char\|3) | `[3,3,3,3,3,3,3,3,3]` | char \| n |
| `k[x]/(x⁴)` | `QQ` (char∤4) | `[4,3,3,3,3,3,3,3,3]` | char ∤ n |
| `k[x]/(x⁴)` | `GF(2)` (char\|4) | `[4,4,4,4,4,4,4,4,4]` | char \| n |
| `kZ₂/J³` sym | `QQ` | `[3,1,1,1,1,1,1,1,1]` | self-inj, symmetric |
| `kZ₃/J⁴` sym | `QQ` | `[4,1,1,1,1,1,1,1,1]` | self-inj, symmetric |
| **`kZ₃/J²` (the ν-TWIST witness)** | `GF(32003)`, `QQ` | `HH^*=[1,1,0,0,0,0,1]` | self-inj, **NON-symmetric**, `π = (1 2 3)` |
| `QuantumCI(2,2,2)` | `QQ` | `HH^*=[2,2,1,0,0,0,0]`, `HH_*=[3,2,2,2,2,2,2]` | Frobenius, **NOT** symmetric |

**`A^e` minimal-resolution Betti ranks `r_0..r_10` (the periodicity signal):**

| algebra | `r_0..r_10` | eventually periodic? | `complexity` |
|---|---|---|---|
| `k[x]/(x²)`, `(x³)`, `(x⁴)` | `[1,1,1,…]` | **YES** (period-2 complex, single generator) | 1 |
| `kZ₂/J³`, `kZ₂/J²` | `[2,2,2,…]` | **YES** (periodic) | 1 |
| `kZ₃/J⁴` | `[3,3,3,…]` | **YES** (periodic) | 1 |
| `QuantumCI(2,2,2)` | `[1,2,3,4,5,6,7,8,9,10,11]` | **NO** (linear growth) | 2 |
| `kA₂` (hereditary) | (finite pd) | n/a | 0 |

**Certificates:** `k[x]/(x^n)` and `kZ_n/J^ℓ` (symmetric): `is_selfinjective=True`,
`is_symmetric=True`, vertex permutation `π = id`, `gorenstein_dimension = right 0 / left 0`.
`QuantumCI(2,2,2)`: `is_selfinjective=True`, `is_symmetric=False`, `π = id` (single vertex,
but `ν²≠id` at the algebra level). **`kZ₃/J²`: `is_selfinjective=True`, `is_symmetric=False`,
`is_weakly_symmetric=False`, vertex permutation `π = {1:2, 2:3, 3:1}` (the 3-cycle) —
live-verified over `GF(32003)` and `QQ`.** `kA₂`: `is_selfinjective=False`, `gorenstein =
right 1 / left 1` (Gorenstein of dim 1, NOT self-injective — the DEFERRED branch witness, M4).

**The ν-twist discriminating witness `kZ₃/J²` (M3).** Every symmetric example has `π = id`,
so the tag rule `(v,w) ↦ (π w, π v)` is indistinguishable from the plain swap `(w, v)` on
them — no test discriminates the `π`-permutation code. `kZ₃/J²` (self-injective,
NON-symmetric, `π =` the 3-cycle) is the **required discriminating witness**: the native
splice's `assert_dd_zero` / `assert_acyclic` hold **iff the vertex permutation `π` is used**
(a control with identity/plain-swap tags FAILS them) — a genuine discriminating self-cert (a
value: `dd=0` holds iff `π` is correct). Its positive threshold pins `ĤH^{n≥1} = HH^n =
[1,0,0,0,0,1]` (self-injective, `d=0`). The **exact untwisted negative values
`ĤH^{-j}(kZ₃/J², A)` are engine-computed** (non-symmetric ⇒ NOT the symmetric mirror
`HH^{j-1}` — a test asserts the engine does not naively mirror); the Bergh–Jorgensen `ν²`-twist
duality `dim ĤH^n(A,A) = dim ĤH^{-(n+1)}(A, A_{ν²})` is the LITERATURE statement they satisfy,
whose **twisted-coefficient recompute is deferred** (it needs the Plan-52 twisted-bimodule
coefficients — out of P76 scope; noted honestly).

**Derived Tate–Hochschild targets (threshold `n≥1` + duality + the explicit period-2 complete
resolution; each derivation shown):**

- **`k[x]/(x^n)` — the FULL closed form (a HARD literature pin; `ĤH^0` derived rigorously).**
  The minimal `A^e`-resolution is period-2 (differentials alternate `·(x⊗1−1⊗x)` and
  `·Σx^i⊗x^{n-1-i}`; `engine/resolutions.py::TruncatedPolynomialResolution`). Applying
  `Hom_{A^e}(−,A)=A`: the cochain differentials alternate `0` and `·(n·x^{n-1})`. Extending
  2-periodically to negative degrees, `ĤH^0 = ker(0)/im(n·x^{n-1})`. Hence:
  ```
  char k ∤ n :  dim ĤH^m(k[x]/(x^n)) = n − 1   for ALL m ∈ ℤ.
  char k | n :  dim ĤH^m(k[x]/(x^n)) = n       for ALL m ∈ ℤ.
  ```
  Instances (cross-checked against the ordinary-HH table + duality + `ĤH^0`):
  `k[x]/(x²)`: `1` in every degree (char∤2), `2` (char∣2, `GF(2)`);
  `k[x]/(x³)`: `2` every degree (char∤3), `3` (char∣3, `GF(3)`);
  `k[x]/(x⁴)`: `3` every degree (char∤4). Periodicity generator `η ∈ ĤH^2`, invertible
  (`periodicity_degree = 2`).
- **`QuantumCI(2,2,2)` over `QQ` — the Bergh–Jorgensen explicit non-symmetric pin (verbatim).**
  `dim ĤH^n = 1 (n=0), 2 (n=1), 1 (n=2), 0 (n∉{0,1,2})`. Cross-checks against the live HH:
  `ĤH^1=2=HH^1`, `ĤH^2=1=HH^2`, `ĤH^{≥3}=0=HH^{≥3}` (threshold `n≥1`, ✓); `ĤH^0=1 ≠ HH^0=2`
  (the Tate stable centre, the `ĤH^0≠HH^0` witness); `ĤH^{n<0}=0` (ALL negatives vanish —
  the **non-symmetric** signature, NOT the mirror of the positives). **Scope:** the *positive*
  part pins over `QQ` through the `positive` route (universal); the *negative vanishing* is a
  Bergh–Jorgensen literature statement the `native` engine reproduces only over `GF(p)` (where
  `q` is a root of unity and the algebra differs), so on `QQ` it is a `positive`-threshold pin
  plus a documented literature value (honest-scope entry (e)).
- **Symmetric Nakayama `kZ_n/J^ℓ` (self-inj, symmetric, periodic — the "symmetric Nakayama
  periods" oracle).** `ĤH^m = HH^m = 1` for `m ≥ 1` (threshold); `dim ĤH^{-j} = dim ĤH^{j-1}
  = HH^{j-1}` **for `j ≥ 2`** (symmetric duality + threshold — the negatives `m ≤ -2` are
  pinned); the `A^e`-resolution is periodic (constant ranks) — `period` engine-detected +
  cross-checked against the invertible element. `kZ₂/J³`: `HH^*=[3,1,1,…]`; `kZ₃/J⁴`:
  `HH^*=[4,1,1,…]`. **`ĤH^0` and `ĤH^{-1}` are the closed 2-cycle (M1): the duality reflects
  `n ↔ -(n+1)`, so `dim ĤH^0 = dim ĤH^{-1}` links them to EACH OTHER, NOT to any positive —
  the duality route CANNOT determine them.** They are `native`-only (engine-computed via the
  complete resolution); the `duality` route returns them `None`. Their FULL closed form is a
  `to-be-verified-at-implementation` value, honest-scope (f) — only `k[x]/(x^n)`'s `ĤH^0`
  (`= n-1`) is hand-derived, from its explicit period-2 resolution.

**What I could NOT live-verify (labelled to-be-verified-at-implementation):**
- The **actual complete resolution and its `ĤH` output** — building the splice IS the plan's
  work, not a probe. The derived targets above (threshold + duality + the period-2 `ĤH^0`
  derivation for `k[x]/(x^n)`) are the standing oracles.
- **`ĤH^0` of the symmetric Nakayama** (I derived it only for `k[x]/(x^n)` via the explicit
  period-2 resolution). The engine computes it; the duality self-consistency + threshold gate
  it.
- **The QCI negative-degree vanishing over `QQ`** (needs a char-0 complete resolution; the
  `native` engine is `GF(p)`). Pinned as literature (Bergh–Jorgensen) + a `QQ` positive-
  threshold pin; a stretch item (CS char-0 complete resolution) would upgrade it to a full
  recomputed pin.

---

## Scope gates (contractual; every boundary is a loud typed refusal)

| surface | scope | refusal / status |
|---|---|---|
| `tate_hochschild` (any route) | path-presented `A` (`A.quiver is not None`) | structure-constant-only → loud `QuiverlabError` (like the minimal/HH engines) |
| `tate_hochschild` (existence of `ĤH^*` at all) | **`A^e` (Iwanaga-)Gorenstein** — guaranteed by `is_selfinjective(A)`, else certified via `gorenstein_dimension` | **not Gorenstein → LOUD refusal** ("Tate–Hochschild cohomology needs the complete resolution, which requires `A^e` Gorenstein; this algebra is not") |
| `native` route (v1) | **self-injective** only (`⇒ A^e` self-inj, `d=0`, `D_{A^e}(P_n)` projective, splice well-defined) over **`GF(p)`**; periodic or not (`QCI` served as a bounded window) | non-self-injective (even if Gorenstein) → **DEFERRED (M4)**: `D_{A^e}(P_n)` is not projective off self-injective, so the dual-splice cannot build it — honest refusal pointing at the `P76-followup` periodicity-extension task; non-`GF(p)` → routes to `duality`/`positive` or refuses |
| **DEFERRED — Gorenstein NON-self-injective + eventually periodic** (M4) | the `native` ring is NOT served in v1 | needs a **periodicity-extension** complete resolution (splice `P_•` with a periodic syzygy repeat, per Usui), not the `D(P_n)`-dual splice — named `P76-followup`; the eventual-periodicity CERTIFICATE (`tate_periodicity`) is still computed for it |
| `duality` route | symmetric (**`π = id`, i.e. `ν² = id`**; `is_symmetric` / exterior) over any exact `Domain` — supplies degrees **`m ≥ 1`** (`= HH`) and **`m ≤ −2`** (`= HH^{|m|−1}`, M1) | `π ≠ id` (e.g. the QCI) → LOUD refusal; **degrees `0` and `−1` are a closed 2-cycle NOT supplied by duality → returned `None` (native-only), never `HH^0`** |
| `positive` route | any Gorenstein `A`, any `Domain` — returns `ĤH^{n≥agrees_from}=HH^n` ONLY (`agrees_from=1` self-injective) | never claims `ĤH^0` (`= HH^0` is ordinary, NOT `ĤH^0` — M2) nor any negative degree; `neg_dims=None`, degree-0 slot `None`, honest note |
| eventual-periodicity certificate | `complexity==1` signal + `is_isomorphic` syzygy cert | not periodic within `max_period` → `period=None` (honest, never guessed); `is_isomorphic` undecidable → the shipped loud refusal unchanged |
| over `GF(p)`, char `≤ dim` | the minimal-engine / `is_isomorphic` char caveat | `char ≤ dim` → loud `QuiverlabError` from the shipped primitives (never a silent wrong verdict) |
| the QCI literature pin | **char-0 (`QQ`)** — `q` not a root of unity | over `GF(p)` the pin does not hold (root of unity); the `native` GF(p) result is `oracle_selfcert`, not the literature value |

All `native`/self-cert batteries run over `GF(p)` primes `(32003, 2, 3)` (to exercise both
`char ∤ n` and `char | n` for `k[x]/(x^n)`); the `duality`/`positive` batteries run over
`QQ`; the QCI literature pin is `QQ`.

---

## API surface (public via `import quiverlab`; exact only)

```python
# --- src/quiverlab/hochschild/tate.py ---
def tate_hochschild(A, top, *, engine="auto", max_cells=4_000_000) -> TateHochschild
def tate_cup_products(A, top, *, engine="auto", max_cells=4_000_000) -> HHProducts   # Plan-35 container
def tate_periodicity(A, *, max_period=12, p=32003) -> TatePeriodicity
def tate_hochschild_block(A, top, *, engine="auto") -> dict                            # runner block

@dataclass(frozen=True)
class TateHochschild:
    top: int                       # N: the block reports ĤH^{-N..N}
    # pos_dims[m] = dim ĤH^m for degree m in 0..N, or None where the ROUTE does not know it.
    # HONESTY (M2): pos_dims[0] = dim ĤH^0 ONLY on the native route; it is None on the
    # positive/duality routes (ĤH^0 is native-only -- it is NOT HH^0, see ordinary_pos).
    pos_dims: tuple                # (ĤH^0|None, ĤH^1, ..., ĤH^N)
    # neg_dims[j-1] = dim ĤH^{-j} for j in 1..N, or None where unknown. HONESTY (M1): native
    # knows all; the duality route knows j>=2 only (ĤH^{-1} is None -- the {0,-1} closed
    # 2-cycle is native-only); the positive route sets neg_dims = None entirely.
    neg_dims: tuple | None         # (ĤH^{-1}|None, ĤH^{-2}, ..., ĤH^{-N}) or None
    hat_hh0: int | None            # dim ĤH^0 (the stable centre; != HH^0 in general). None off native (M1/M2).
    agrees_from: int               # d+1: ĤH^m == HH^m for m >= agrees_from  (d = Gor.dim A^e; 1 if self-inj)
    ordinary_pos: tuple            # (HH^0, ..., HH^N) -- ordinary HH, the agreement anchor. HH^0 is NOT ĤH^0.
    period: int | None             # the eventual period p (Ω^{n+p}≅Ω^n), or None
    periodicity_degree: int | None # the degree of the invertible element (== period when native cup ran)
    scope: str                     # "self_injective" | "positive_only"  (gorenstein_periodic ring DEFERRED, M4)
    engine: str                    # "native" | "duality" | "positive"
    cup: dict | None               # Plan-35 blocks() of the Tate cup, or None
    references: tuple
    note: str | None

@dataclass(frozen=True)
class TatePeriodicity:
    period: int | None             # least p with Ω^p_{A^e}(A) ≅ A (is_isomorphic-certified), or None
    eventually_periodic: bool      # a shift n>=0 with Ω^{n+p}≅Ω^n found within budget
    invertible_element_degree: int | None   # Usui: degree of an invertible homogeneous ĤH element
    complexity: int                # the shipped A^e-Betti complexity (1 signal)
    note: str

# --- core/algebra.py delegate (beside cup_products) ---
Algebra.tate_hochschild(top, engine="auto", max_cells=4_000_000) -> TateHochschild
```

`TateHochschild.__bool__` is NOT defined (a data report, not a predicate) — mirroring
`HHProducts`/`ARQuiver`.

---

## Global Constraints

- Python is always `.venv/bin/python`; tests
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **Prerequisites, all MERGED on `dev` (verified by live import at spec time):** the minimal
  `A^e` engine (`engine/resolutions_minimal.py`: `minimal_resolution`, `AeEngine`,
  `minimal_cohomology_dims`, `nullspace_mod_p`, the corner machinery); the engine adapter
  (`engine/adapter.py::to_engine`, `engine_cohomology_dims`); the Plan-35 product surface
  (`hochschild/products.py::HHProducts`/`ProductTable`, `Algebra.cup_products`,
  `resolutions_cs/products.py`) and the Plan-20/21 diagonal (`resolutions_cs/diagonal.py`);
  the duality/opposite machinery (`modules/duality.py`, `Algebra.opposite`,
  `families/tensor.py::TensorProduct`); the self-injective / Frobenius certificates
  (`modules/ext.py::is_selfinjective`, `invariants/frobenius.py::is_symmetric_generic`,
  `Algebra.nakayama_automorphism`); the homdim certificates
  (`modules/homdims.py::omega_periodicity`, `gorenstein_dimension`, `is_gorenstein`;
  `Algebra.complexity`). **If any is not merged when a worker picks this up, STOP and
  escalate** — do not fork the minimal engine or the product diagonal.
- **Buckets auto-assigned by directory** (`tests/conftest.py`, verified at spec time):
  `_DEEP_DIRS = (engine, resolutions_cs, modules, families, batch, specseq)` → deep;
  **`tests/hochschild/` is NOT in that list → fast**; `tests/webapp/`, `tests/gui/`,
  `tests/hpc/` → fast (extras-gated, unmarked per Plan-32). The heavy complete-resolution
  builds therefore live in **`tests/engine/test_complete_resolution_p76.py`** (deep — they
  share the minimal-`A^e` budget). The public-surface oracle tests in
  **`tests/hochschild/test_tate_hochschild_p76.py`** are **fast by default**, so any test that
  does a real GF(p) complete-resolution build MUST carry an explicit `@pytest.mark.deep`
  (explicit bucket wins in `conftest`) — never mislabel a heavy build as fast.
- **The `native` Tate engine runs over `GF(p)` int64** (the minimal-`A^e` sibling; primes
  `(32003, 2, 3)`). The `duality`/`positive` routes run over any exact `Domain`. `char ≤ dim`
  over `GF(p)` inherits the loud minimal-engine / `is_isomorphic` refusal unchanged.
- **Multi-vertex feeds the minimal engine UN-adapted** (the live-reproduced Plan-18 gotcha):
  `to_engine(A)` for multi-vertex, `to_engine(A.unit_adapted())` only for single-vertex —
  Task I1 encapsulates this in one helper with a regression pin.
- **`assert_dd_zero` across the splice joint is MANDATORY** — a Tate result is never returned
  without it (the NEW machinery's primary self-cert). `assert_acyclic` in-window likewise.
- **Honest semi-decision contract (metaplan §1.3).** Every route is **complete** iff its scope
  gate holds and its budget/window did not truncate; otherwise `status`/`note` says so and
  `neg_dims`/`period` is `None` — **never a guessed dim/period**. `ĤH^0 ≠ HH^0` is stated, not
  hidden.
- **No floats in `src/`.** Dims/periods `int`; verdicts/scope `bool`/`str`; cup constants exact
  strings (Plan-35 `ProductTable`). Client-side float only (`gui.js`, exempt).
- **Composition left-to-right**; all refusals `QuiverlabError`.
- **Plan-32 markers** (orthogonal; `oracle_*` FORBIDDEN in `tests/{webapp,gui,hpc}`):
  - `oracle_literature`: `k[x]/(x^n)` full closed form `ĤH^m = n−1` (char∤n) / `n` (char∣n)
    for `n∈{2,3,4}`, all `m∈[-N,N]`; the QCI `QQ` positive-threshold pin `ĤH^{n≥1}=[2,1,0,…]`
    (+ the Bergh–Jorgensen full-ring literature value as a documented statement); the
    threshold `agrees_from` values; the symmetric-Nakayama `ĤH^{m≥1}=1`.
  - `oracle_crossengine`: **positive-degree agreement** `ĤH^m == HH^m` for `m ≥ agrees_from`
    (the NEW `native` engine vs the shipped `hochschild_cohomology` — the headline
    cross-engine oracle); the **Tate cup positive part** vs the shipped `A.cup_products`
    (Plan-35) in-window; **`duality` dims ≡ `native` dims** for symmetric `k[x]/(x^n)` /
    Nakayama over shared `GF(p)`.
  - `oracle_selfcert`: `assert_dd_zero` across every splice joint + `assert_acyclic`
    in-window; the **symmetric duality self-consistency** `dim ĤH^{-j} == dim ĤH^{j-1}` on
    the native output (symmetric inputs); the **invertible-element identity** `η·η^{-1}=1 ∈
    ĤH^0` for `k[x]/(x^n)` (`periodicity_degree=2`); the `Ω`-periodicity `is_isomorphic`
    certificate; the not-Gorenstein / non-symmetric / structure-constant / char≤dim loud
    refusals.
  - `qpa`: **QPA 1.37 has NO Tate/singular-Hochschild surface** (live `NamesGVars()` sweep for
    `Tate`/`Singular`/`Complete` — the Plan-35 precedent; the test FAILS if that ever changes).
    The **ordinary** `HH^n` that `ĤH^n` agrees with for `n ≥ agrees_from` is QPA-crosschecked
    via the existing `HochschildCohomology` bridge on the Dynkin/Nakayama zoo (an input-level
    anchor for the agreement, not a Tate verdict compare).
- **Mid-merge-train counts.** Absolute suite counts drift; the verification task **recounts the
  oracle-class table at merge time** (`tests/release/test_oracle_classes.py`; paste live
  numbers, claim only this plan's deltas).
- Every merge updates `docs/verification.md` (new oracle rows + recounted class table) and adds
  citations to `references.bib` + `registry.py` (`bibtex()` hard-fails if they disagree).
  Conventional commits; green at every commit.

**No-red-commits discipline.** Each Task ships a self-contained green slice: the engine
splice + `assert_dd_zero` before dims before cup before GUI. Result fields default to
`None`/empty so intermediate task boundaries are green (no forward reference, no `xfail`).

---

# Task group I — the complete-resolution splice engine (GF(p))

### Task I1: `complete_resolution` — splice `P_•` with its Nakayama-twisted `A^e`-dual + `assert_dd_zero`

**Files:** create `src/quiverlab/engine/complete_resolution.py`; test
`tests/engine/test_complete_resolution_p76.py` (deep).

- [ ] **Step 1: Failing tests (self-cert, the NEW machinery's spine).**
```python
# tests/engine/test_complete_resolution_p76.py
"""Complete resolution of A over A^e by splicing the minimal A^e-resolution with its
Nakayama-twisted A^e-dual (Plan 76 / R3; Bergh-Jorgensen JNCG 7 (2013)). GF(p);
self-injective / eventually-periodic scope; d.d=0 across the splice joint is mandatory."""
import pytest
from quiverlab import truncated_polynomial, NakayamaAlgebra
from quiverlab.fields import GF
from quiverlab.engine.complete_resolution import complete_resolution
selfcert = pytest.mark.oracle_selfcert

@selfcert
@pytest.mark.parametrize("n,p", [(2,32003),(2,2),(3,32003),(3,3),(4,32003)])
def test_dd_zero_and_acyclic_kxn(n, p):
    A = truncated_polynomial(n, field=GF(p))
    T = complete_resolution(A, 5, p)
    T.assert_dd_zero()            # every joint, incl. d_0.d_1 and d_-1.d_0 (loud on failure)
    T.assert_acyclic(window=(-4, 4))

@selfcert
def test_dd_zero_multivertex_symmetric_nakayama():
    """The multi-vertex minimal engine must be fed UN-adapted (Plan-18 gotcha, live-repro):
    unit_adapted() trips radical_basis on multi-vertex; to_engine(A) does not."""
    A = NakayamaAlgebra(n=2, l=3, cyclic=True, field=GF(32003))   # dim 6, symmetric, periodic
    T = complete_resolution(A, 4, 32003)
    T.assert_dd_zero()
    assert T.term_rank(0) == T.term_rank(-1)     # T_{-1} = dual of T_0 (swapped tags)

@selfcert
def test_dd_zero_kz3j2_nu_twist_REQUIRES_pi(monkeypatch):
    """M3 discriminating witness: kZ3/J2 is self-injective, NON-symmetric, vertex permutation
    pi = (1 2 3). The correct pi-permuted dual passes dd=0/acyclic; forcing pi = identity (the
    plain swap) must FAIL them -- the value that discriminates the tag-permutation code."""
    A = NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(32003))   # dim 6, non-symmetric
    T = complete_resolution(A, 4, 32003)          # correct pi -> passes
    T.assert_dd_zero(); T.assert_acyclic(window=(-3, 3))
    # control: with identity tags (plain swap, ignoring pi) the splice is NOT a complex:
    import quiverlab.engine.complete_resolution as cr
    monkeypatch.setattr(cr, "nakayama_permutation", lambda A: {v: v for v in A.quiver.vertices})
    with pytest.raises(AssertionError):
        cr.complete_resolution(A, 4, 32003).assert_dd_zero()

@selfcert
def test_non_selfinjective_deferred():
    """M4: a Gorenstein-but-NOT-self-injective algebra (kA2, gldim 1) is DEFERRED -- the
    D(P_n)-dual splice cannot build it -> loud refusal (P76-followup periodicity-extension)."""
    from quiverlab.errors import QuiverlabError
    from quiverlab import linear_path_algebra
    A = linear_path_algebra(2, field=GF(32003))   # Gorenstein (gldim 1) but NOT self-injective
    with pytest.raises(QuiverlabError):
        complete_resolution(A, 4, 32003)

@selfcert
def test_not_gorenstein_refused():
    from quiverlab.errors import QuiverlabError
    from quiverlab.families.radical_square_zero import RadicalSquareZero
    from quiverlab import Quiver
    A = RadicalSquareZero(Quiver([1,2,3], {"a":(1,2),"b":(2,3)}), field=GF(32003))  # non-selfinj
    with pytest.raises(QuiverlabError):
        complete_resolution(A, 4, 32003)          # not self-injective -> refuse (M4)
```
- [ ] **Step 2:** confirm failure (`ModuleNotFoundError`).
- [ ] **Step 3: Implement.** `_engine_form(A)` = `to_engine(A if len(list(A.quiver.vertices))
  > 1 else A.unit_adapted())` (the multi-vertex gotcha). Build `P_•` via `minimal_resolution(_engine_
  form(A), N, p)` (corner tags on `eng.corner_tags`). Add a `nakayama_permutation(A) -> dict`
  helper (`π` from `soc(P_v)`, reusing the `is_selfinjective` socle logic — NOT the `d×d`
  `nakayama_automorphism()` matrix, M3). Build the dual half `T_{-n-1}` = the Nakayama-twisted
  `A^e`-dual of `P_n`: swap+`π`-permute the corner tags `(v,w)↦(π w, π v)` (`π` the VERTEX
  permutation; `π=id` iff weakly symmetric), transpose the `A^e` differential blocks (the
  corner-transpose already used by `minimal_cohomology_dims`). `# PIN` the exact `π`-direction/
  sign — arbitrated by `assert_dd_zero`/`assert_acyclic` on `kZ₃/J²`, not assumed. Splice joint
  `d_0 = D(μ)∘μ`. `CompleteResolution` carries `terms`, `diffs`, `corner_tags` (both signs),
  `assert_dd_zero` (rank/product mod p, loud), `assert_acyclic` (rank–nullity in window),
  `term_rank(n)`. **Scope gate (M4):** refuse loudly unless `is_selfinjective(A)` — the
  dual-splice is well-defined ONLY for self-injective `A` (`D_{A^e}(P_n)` projective iff `A^e`
  self-injective); a non-self-injective Gorenstein `A` gets the DEFERRED refusal pointing at
  `P76-followup` (the eventual-periodicity certificate is computed by Task A2 regardless).
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(engine): complete resolution of A over A^e by splicing the minimal A^e-resolution with its Nakayama-twisted dual (Plan 76/R3); d.d=0 across the joint + in-window acyclicity self-cert; GF(p), multi-vertex un-adapted`.

### Task I2: `tate_cohomology_dims` / `tate_homology_dims` + the positive-agreement self-cert

**Files:** modify `src/quiverlab/engine/complete_resolution.py`; extend the test.

- [ ] **Step 1: Failing tests.**
```python
@pytest.mark.oracle_literature
@pytest.mark.parametrize("n,p,expect", [
    (2,32003,1),(2,2,2),(3,32003,2),(3,3,3),(4,32003,3)])   # ĤH^m = n-1 (char∤n) / n (char|n) ALL m
def test_kxn_full_tate_ring_constant(n, p, expect):
    from quiverlab.engine.complete_resolution import tate_cohomology_dims
    A = truncated_polynomial(n, field=GF(p))
    dims = tate_cohomology_dims(A, 5, p)               # {m: dim} for m in [-5,5]
    for m in range(-5, 6):
        assert dims[m] == expect                        # constant across ALL degrees (period-2 sym.)

@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("n,p", [(3,32003),(2,2)])
def test_positive_agrees_with_ordinary_hh(n, p):
    """Bergh-Jorgensen threshold: ĤH^m == HH^m for m >= 1 (self-injective, d=0)."""
    from quiverlab.engine.complete_resolution import tate_cohomology_dims
    A = truncated_polynomial(n, field=GF(p))
    tate = tate_cohomology_dims(A, 6, p)
    hh = A.hochschild_cohomology(6).dims
    for m in range(1, 7):
        assert tate[m] == hh[m]                         # NEW engine vs shipped engine
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement.** `tate_cohomology_dims` = `H^m(Hom_{A^e}(𝕋, A))` (collapse
  `a·w·b`, the Plan-16 cohomology side, over both signs); `tate_homology_dims` =
  `H_m(A ⊗_{A^e} 𝕋)` (collapse `b·w·a`, swapped tags). Reuse the `minimal_cohomology_dims`
  contraction per degree, extended to negative degrees off the dual half. **Positive-agreement
  self-cert (inside the engine):** for `m ≥ agrees_from`, assert `tate[m] == HH^m` (drift gate
  — the complete resolution AGREES with the projective resolution in high degree by
  construction; a mismatch is a splice bug). `agrees_from = 1` for self-injective.
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(engine): tate_cohomology_dims / tate_homology_dims off the complete resolution (Plan 76/R3); k[x]/(x^n) full ring pins + positive-agreement drift gate vs shipped HH`.

---

# Task group II — the public Tate surface + the periodicity certificate + the cup

### Task A1: `tate_hochschild` — the three-route public surface + `Algebra` delegate

**Files:** create `src/quiverlab/hochschild/tate.py`; the `Algebra` delegate in
`core/algebra.py`; test `tests/hochschild/test_tate_hochschild_p76.py`.

- [ ] **Step 1: Failing tests.**
```python
# tests/hochschild/test_tate_hochschild_p76.py
"""Public Tate-Hochschild surface (Plan 76 / R3). Native GF(p) splice, duality dims
(symmetric, any Domain), positive threshold (any Gorenstein). Bergh-Jorgensen JNCG 7."""
import pytest
from quiverlab import truncated_polynomial, NakayamaAlgebra, QuantumCI, linear_path_algebra
from quiverlab.fields import QQ, GF
from quiverlab.errors import QuiverlabError
lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert

@lit
def test_kx3_native_full_ring():
    A = truncated_polynomial(3, field=GF(32003))
    r = A.tate_hochschild(5, engine="native")
    # NATIVE knows every degree, incl. ĤH^0: full ring = 2 in every degree (n-1=2).
    assert r.pos_dims == (2,2,2,2,2,2) and r.neg_dims == (2,2,2,2,2)   # ĤH^0..5 and ĤH^-1..-5
    assert r.hat_hh0 == 2 and r.agrees_from == 1 and r.scope == "self_injective"

@selfcert
def test_kz3j2_native_nu_twist_witness():
    """The ν-twist discriminating witness (M3): kZ3/J2 is self-injective, NON-symmetric,
    vertex permutation π = (1 2 3). The native splice's dd=0/acyclic hold ONLY with the
    π-permuted tags; positive threshold ĤH^{n>=1} = HH^n = [1,0,0,0,0,1] (self-inj, d=0);
    negatives are engine-computed and NOT the symmetric mirror (non-symmetric)."""
    from quiverlab.hochschild.tate import nakayama_permutation
    A = NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(32003))
    assert A.is_selfinjective() and not A.is_symmetric()
    assert nakayama_permutation(A) == {1: 2, 2: 3, 3: 1}          # live-verified 3-cycle
    r = A.tate_hochschild(6, engine="native")
    hh = A.hochschild_cohomology(6).dims                          # [1,1,0,0,0,0,1]
    assert tuple(r.pos_dims[1:]) == tuple(hh[1:])                 # ĤH^{n>=1} = HH^n (threshold)
    # non-symmetric => the negatives are NOT the naive symmetric mirror HH^{j-1}:
    assert not _is_symmetric_mirror(r.neg_dims, r.pos_dims)       # discriminating (see helper)

@xeng
def test_positive_threshold_universal_QQ():
    """positive route: ĤH^{n>=1} = HH^n for any Gorenstein A over any Domain.
    Degree 0 (pos_dims[0]) is None on the positive route -- ĤH^0 is native-only, NOT HH^0 (M2)."""
    A = truncated_polynomial(3, field=QQ)
    r = A.tate_hochschild(6, engine="positive")
    hh = A.hochschild_cohomology(6).dims
    assert r.pos_dims[0] is None                                  # ĤH^0 NOT claimed (native-only)
    assert r.pos_dims[1:] == tuple(hh[1:])                        # ĤH^{n>=1} = HH^n
    assert r.ordinary_pos == tuple(hh)                            # HH^0..6, the anchor (HH^0 != ĤH^0)
    assert r.neg_dims is None and r.agrees_from == 1 and r.hat_hh0 is None

@lit
def test_qci_positive_threshold_QQ():
    """QCI over QQ (q=2 not a root of unity; quiverlab xy+q·yx == BJ's XY-qYX for generic q):
    ĤH^{n>=1} = HH^n = [2,1,0,0,0]. Degree 0 is HH^0=2 in ordinary_pos but pos_dims[0]=None
    (BJ's Tate ĤH^0 = 1 != HH^0 = 2 -- honest-scope (g); the native GF(p) route cannot
    reproduce BJ's char-0 value, q being a root of unity there)."""
    A = QuantumCI(2, 2, 2, field=QQ)
    r = A.tate_hochschild(5, engine="positive")
    assert r.pos_dims[0] is None                                  # ĤH^0 NOT claimed (M2)
    assert r.pos_dims[1:] == (2, 1, 0, 0, 0)                      # ĤH^{n>=1} = HH^n
    assert r.ordinary_pos == (2, 2, 1, 0, 0, 0)                   # HH^0..5 (HH^0=2, the ordinary anchor)
    assert r.neg_dims is None

@xeng
def test_duality_route_degree0_absent_QQ():
    """Duality route (symmetric k[x]/x^3 over QQ): supplies m>=1 (=HH) and m<=-2 (=HH^{|m|-1});
    degrees 0 and -1 are a closed 2-cycle NOT supplied by duality -> None (M1). Never HH^0."""
    A = truncated_polynomial(3, field=QQ)
    r = A.tate_hochschild(5, engine="duality")
    hh = A.hochschild_cohomology(5).dims                          # [3,2,2,2,2,2]
    assert r.pos_dims[0] is None and r.hat_hh0 is None            # ĤH^0 native-only (NOT HH^0=3)
    assert r.pos_dims[1:] == tuple(hh[1:])                        # ĤH^{n>=1} = HH^n = 2..
    assert r.neg_dims[0] is None                                  # ĤH^{-1} native-only (closed cycle)
    assert r.neg_dims[1:] == tuple(hh[1:5])                       # ĤH^{-j}=HH^{j-1}, j=2..5 (all 2)

@selfcert
def test_qci_duality_refused():
    A = QuantumCI(2, 2, 2, field=QQ)                   # Frobenius but NOT symmetric (nu^2 != id, pi=id but algebra-nu^2!=id)
    with pytest.raises(QuiverlabError):
        A.tate_hochschild(4, engine="duality")

@selfcert
def test_native_non_selfinjective_deferred():
    """M4: kA2 is Gorenstein (gldim 1) but NOT self-injective -> the D(P_n)-dual splice cannot
    build it (D(P_n) projective iff A^e self-injective). v1 native scope is self-injective;
    the general Gorenstein-periodic ring is DEFERRED (P76-followup) -> loud refusal."""
    A = linear_path_algebra(2, field=GF(32003))
    assert not A.is_selfinjective() and A.is_gorenstein()
    with pytest.raises(QuiverlabError):
        A.tate_hochschild(4, engine="native")

@selfcert
def test_not_gorenstein_refused():
    from quiverlab.families.radical_square_zero import RadicalSquareZero
    from quiverlab import Quiver
    A = RadicalSquareZero(Quiver([1,2,3], {"a":(1,2),"b":(2,3)}), field=GF(32003))
    with pytest.raises(QuiverlabError):
        A.tate_hochschild(4)
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement `tate_hochschild(A, top, engine)`** dispatching the three routes
  (Mathematical foundation). `_require_gorenstein(A)` (self-injective OR `is_gorenstein` True;
  loud else). `native` (GF(p)): **require `is_selfinjective(A)`** (M4 — else loud refusal
  pointing at the DEFERRED `P76-followup` periodicity-extension task); `complete_resolution` +
  `tate_cohomology_dims`; `scope="self_injective"`; `pos_dims` (incl. `pos_dims[0]=ĤH^0`),
  `neg_dims` (all `-1..-N`), `hat_hh0` all filled. `duality` (`π=id` gate via `is_symmetric` /
  exterior; loud else): `pos_dims[0]=None`, `pos_dims[1:]=HH^{1..N}`, `neg_dims[0]=None`,
  `neg_dims[1:]=HH^{1..N-1}` (`ĤH^{-j}=HH^{j-1}`, `j≥2`), `hat_hh0=None`, `ordinary_pos=HH^{0..N}`.
  `positive`: `pos_dims=[None, HH^1..HH^N]`, `neg_dims=None`, `hat_hh0=None`, `ordinary_pos`
  filled, honest note. **Every route sets `ordinary_pos = HH^0..HH^N` (the anchor; `HH^0` is
  ordinary, never relabelled `ĤH^0`).** `agrees_from = 1` (self-injective; the only v1 native
  scope). Add a **`nakayama_permutation(A) -> dict`** helper (extract `π` from `soc(P_v)`,
  reuse the `is_selfinjective` socle logic) — used by the native tag rule and the M3 witness
  test. `Algebra.tate_hochschild` delegate added beside `cup_products`.
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(hochschild): public tate_hochschild -- native GF(p) splice + duality dims (symmetric) + positive threshold; Gorenstein/self-injective gate, loud refusals (Plan 76/R3)`.

### Task A2: `tate_periodicity` — the eventual-periodicity certificate (Ω-periodicity reuse + Usui invertible element)

**Files:** modify `src/quiverlab/hochschild/tate.py`; extend the test.

- [ ] **Step 1: Failing tests.**
```python
@selfcert
@pytest.mark.parametrize("A,p_expect", [
    (truncated_polynomial(3, field=GF(32003)), 2),           # k[x]/(x^3): period 2
    (NakayamaAlgebra(n=2,l=3,cyclic=True, field=GF(32003)), None),  # periodic; period engine-detected
])
def test_periodicity_certificate(A, p_expect):
    from quiverlab.hochschild.tate import tate_periodicity
    cert = tate_periodicity(A)
    assert cert.eventually_periodic is True and cert.complexity == 1
    if p_expect is not None:
        assert cert.period == p_expect

@selfcert
def test_qci_not_periodic():
    """QCI is self-injective (served) but NOT eventually periodic (complexity 2) -- the
    certificate is honest None, and the native route still serves it as a bounded window."""
    from quiverlab.hochschild.tate import tate_periodicity
    A = QuantumCI(2, 2, 2, field=GF(32003))
    cert = tate_periodicity(A)
    assert cert.period is None and cert.eventually_periodic is False and cert.complexity == 2
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement `tate_periodicity`.** Cheap signal `A.complexity(N)` (`==1` ⇒
  periodic candidate; `>=2` ⇒ not). Rigorous certificate: reuse the
  `modules/homdims.py::omega_periodicity` pattern for `A` over `A^e` — `# PIN`: if
  `TensorProduct(A, A.opposite())` builds `A^e` and `A` presents as its regular module,
  call `omega_periodicity` verbatim; else the engine-level syzygy `is_isomorphic` certificate
  (compare `Ω^{n+p}_{A^e}(A)` corner modules from the minimal-`A^e` resolution). Usui
  cross-check: when the native Tate cup ran, set `invertible_element_degree` = the least `p≥1`
  with an invertible `η∈ĤH^p` (`η·η^{-1}=1∈ĤH^0`); assert it equals `period`.
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(hochschild): tate_periodicity -- Omega-periodicity of A over A^e (is_isomorphic-certified, homdims reuse) + Usui invertible-element cross-check (Plan 76/R3)`.

### Task A3: `tate_cup_products` — the ℤ-graded Tate ring (native GF(p)) + positive cross-check

**Files:** modify `src/quiverlab/hochschild/tate.py`; the `Algebra` delegate; extend the test.

- [ ] **Step 1: Failing tests.**
```python
@xeng
def test_tate_cup_positive_matches_plan35():
    """The positive part of the Tate cup == the shipped ordinary cup (Plan 35) in-window."""
    A = truncated_polynomial(3, field=GF(32003))
    tate = A.tate_cup_products(2)               # HHProducts over the Tate basis
    ordinary = A.cup_products(2)
    # basis-independent compare (Plan-35 rule): dims + flattened rank mod p on p+q>=1 tables
    assert _basis_free_equal(tate, ordinary, only_nonneg=True)

@selfcert
def test_invertible_periodicity_element_kx2():
    """k[x]/(x^2): the degree-2 periodicity generator eta is invertible: eta . eta^{-1} = 1 in ĤH^0."""
    A = truncated_polynomial(2, field=GF(32003))
    r = A.tate_hochschild(3, engine="native")
    assert r.periodicity_degree == 2 and r.period == 2
```
- [ ] **Step 2:** confirm failure.
- [ ] **Step 3: Implement `tate_cup_products`.** Build the Tate cup on the complete resolution
  via a diagonal approximation `Δ: 𝕋 → 𝕋 ⊗_A 𝕋` extending the Plan-20 `resolutions_cs/
  diagonal.py` construction to the two-sided complex (the sign-free `a·w·b` collapse of
  Plan-20 `cup.py`), returning a Plan-35 `HHProducts` (kind `"cup"`, basis
  `"tate/GF(p)"`). **In-window anchor:** the `p+q ≥ agrees_from` part equals the shipped
  `A.cup_products` (basis-independent compare — the Plan-35 cross-engine rule). The invertible
  periodicity element is read off the cup (Task A2). **Honest scope:** the Tate cup is
  `native`/`GF(p)` only; over other Domains the surface ships dims (`duality`/`positive`)
  without the full ring — stated on the verification page. `Algebra.tate_cup_products` delegate.
- [ ] **Step 4:** run; **Step 5:** commit
  `feat(hochschild): tate_cup_products -- the Z-graded Tate ring via the two-sided diagonal (Plan 20/21 extended); positive part cross-checks Plan-35 cup; invertible periodicity element (Plan 76/R3)`.

---

# Task group III — GUI / webapp + verification

### Task G1: the `tate_hochschild` range kind (all three tiers, i18n ×4, negative-degree table)

A range kind `tate_hochschild:0..N` (the block reports `ĤH^{-N..N}`) following the
`hh_cohomology`/`cup` precedent (a `-top` number input, degree-range grammar, algebra-only).

**Grammar / canonical-key decision (MINOR c).** The form is pinned to exactly
`tate_hochschild:0..N`: the shared `parse_compute_item` **already rejects a nonzero `lo`
loudly** (`hpc/spec.py` lines 268–271 / `webapp/server/schema.py`, "compute range must start
at 0") and a bare `tate_hochschild` (no range) raises "needs a degree range" in dispatch, so
`0..5` is the ONLY accepted spelling — there is no `0..5` vs `1..5` duplicate (`1..5` never
parses). The canonical key is therefore a function of `(kind, hi)` only, byte-stable; no
normalization beyond the existing guard is needed. The plan adds a one-line regression pin
that `tate_hochschild:1..5` raises `SpecError` (the guard covers the new kind for free).

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (the `_RANGE` grammar already accepts
  `tate_hochschild:0..N`; add the `_dispatch` branch + `_snip` recipe), `docs/gui/runner.py`
  (the Pyodide twin — matching `compute_one` branch + `_snip` + cost, shape-identical),
  **`webapp/server/schema.py`** (the THIRD grammar-parse site — the `_RANGE` regex already
  covers it; add to any kind allow-list), `docs/gui/gui.js` + `webapp/static/gui/gui.js`
  (checkbox `qlgui-tate` + `qlgui-tate-top` number input, `S.ids`, push-list
  `"tate_hochschild:0.." + top`, `renderBlock`, `scheduleProbe`), `webapp/templates/
  index.html` (one checkbox), `webapp/server/i18n/{en,es,fr,zh}.json` (**all four locales** —
  the block key chain + `pick.kind.tate_hochschild`; exact key parity gated by
  `tests/webapp/test_i18n.py`), `src/quiverlab/trace/results_html.py` (`_HEADINGS` + a render
  branch — the **two-row degree table**: `ĤH^0..ĤH^N` and `ĤH^{-1}..ĤH^{-N}`),
  `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py` (ONE golden).
- Test: `tests/webapp/test_tate_kind_p76.py`, `tests/gui/test_tate_runner_twin_p76.py`.

**Block shape** (one shared `tate_hochschild_block(A, top, engine)`). Degrees the route does
not know are `null` (JSON), never a wrong number (M1/M2):
```python
{"kind": "tate_hochschild", "top": N, "engine": "native"|"duality"|"positive",
 "scope": "self_injective"|"positive_only",     # "gorenstein_periodic" ring DEFERRED (M4)
 # pos_dims[0] = dim ĤH^0 ONLY on native; null on positive/duality (ĤH^0 native-only, != HH^0):
 "pos_dims": [dim ĤH^0 | null, dim ĤH^1, ..., dim ĤH^N],
 # neg_dims: native all; duality [null, ĤH^-2..ĤH^-N] (ĤH^-1 the {0,-1} native-only cycle);
 # positive null:
 "neg_dims": [dim ĤH^{-1} | null, ..., dim ĤH^{-N}] | null,
 "hat_hh0": int | null,                          # null off native (M1/M2)
 "agrees_from": int, "ordinary_pos": [HH^0, ..., HH^N],   # HH^0 is ORDINARY, never labelled ĤH^0
 "period": int|null, "periodicity_degree": int|null,
 "note": str|null,
 "references": ["wang_singular_hh", "bergh_jorgensen_tate", "usui_tate_periodic",
                "keller_singular_hh", "assem_book"],
 "citations": [...]}
# refusal (not Gorenstein / non-self-injective native / structure-constant / char<=dim)
#   -> {"error": msg, "references": [...]}, never a 500.
```
The cup tables are **NOT** shipped in the block by default (they can be large); the block
reports the ring's dims + `periodicity_degree` (the invertible element's degree). The report
renders the two-row degree table (a `null` cell shows as "—", stated native-only) + the
agreement/threshold line + (native) a small cup sample. **The report/GUI label the degree-0
column `ĤH^0` only when non-null; the `ordinary_pos` HH row is labelled Hochschild cohomology
`HH^•`, distinct from `ĤH^•`.**

- [ ] **Step 1: Failing cross-runner tests** (unmarked — extras-gated dir): the `k[x]/(x^3)`
  `GF(32003)` NATIVE block (`tate_hochschild:0..2`) has `pos_dims == [2,2,2]`, `neg_dims ==
  [2,2]`, `hat_hh0 == 2`, `agrees_from == 1`, `scope == "self_injective"`,
  `"bergh_jorgensen_tate"` and `"wang_singular_hh"` in the citation keys; **a `positive`-route
  block (e.g. over `QQ`) has `pos_dims[0] is None`, `hat_hh0 is None`, `neg_dims is None`,
  `ordinary_pos[0] == HH^0` (M2 honesty in the block)**; twin parity
  (`json.dumps(sort_keys=True)` equality across `hpc/spec.py` and `docs/gui/runner.py`).
- [ ] **Step 2: Implement** the `_dispatch` branch (range kind — `top=item.hi`, loud if None,
  like `hh_cohomology`), the twin `compute_one` branch, the `gui.js` checkbox + `scheduleProbe`
  ETA (size like `hh_cohomology` with a ~3× cost factor — builds both halves + cup; a `_snip`
  `"tate_hochschild": "A.tate_hochschild(%d)"`), `_HEADINGS["tate_hochschild"] =
  "Tate–Hochschild cohomology"` + the two-row negative/positive degree-table render branch.
- [ ] **Step 3: Add ONE golden** (`tate_hochschild_kx3`) to `_runner_goldens.json`; note it in
  `test_runner_delegation.py`'s change-log docstring; confirm existing goldens byte-identical
  first. **Canonical-key stability:** schema-v1 algebra-only range kind, no `module` block →
  canonicalizes unchanged (confirm an existing family request's key is byte-stable).
- [ ] **Step 4:** run `tests/webapp/test_tate_kind_p76.py tests/webapp/test_runner_delegation.py
  tests/gui/test_tate_runner_twin_p76.py tests/webapp/test_i18n.py tests/hpc -q`.
- [ ] **Step 5:** commit
  `feat(gui,webapp,hpc): tate_hochschild range kind (Plan 76) -- positive + NEGATIVE degree tables, both runners byte-identical, grammar in 3 sites, i18n x4, one golden`.

**i18n keys — ALL FOUR LOCALES** (mirror `pick.kind.hh_cohomology`):

| key | en | es | fr | zh |
|---|---|---|---|---|
| `pick.kind.tate_hochschild` | Tate–Hochschild cohomology | Cohomología de Tate–Hochschild | Cohomologie de Tate–Hochschild | Tate–Hochschild 上同调 |
| `block.tate_hochschild.title` | Tate–Hochschild cohomology | Cohomología de Tate–Hochschild | Cohomologie de Tate–Hochschild | Tate–Hochschild 上同调 |
| `block.tate_hochschild.positive` | Positive degrees (ĤHⁿ, n ≥ 0) | Grados positivos | Degrés positifs | 正次数 |
| `block.tate_hochschild.negative` | Negative degrees (ĤH⁻ⁿ) | Grados negativos | Degrés négatifs | 负次数 |
| `block.tate_hochschild.agrees` | Agrees with ordinary HH from degree | Coincide con HH ordinaria desde el grado | Coïncide avec HH ordinaire à partir du degré | 从此次数起与普通 HH 一致 |
| `block.tate_hochschild.period` | Eventual period | Período eventual | Période éventuelle | 最终周期 |

### Task G2: verification page, citations, README, suite gate

**Files:** `src/quiverlab/citations/references.bib` + `registry.py`; `docs/verification.md`;
`README.md`; `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick P76); existing release gates.

- [ ] **Step 1: Citations (BibTeX-VERIFIED only; the four papers were read at authoring).**
  Add to `references.bib` (delete any `# PIN` field that cannot be confirmed at citation-add
  rather than guess):
```bibtex
@misc{Wang2015singular,
  author = {Wang, Zhengfang},
  title  = {Singular {H}ochschild cohomology and {G}erstenhaber algebra structure},
  year   = {2015}, note = {arXiv:1508.00190}}        % definitional origin; journal ref at citation-add
@article{Keller2018singular,
  author  = {Keller, Bernhard},
  title   = {Singular {H}ochschild cohomology via the singularity category},
  journal = {C. R. Math. Acad. Sci. Paris}, volume = {356}, number = {11-12},
  pages   = {1106--1111}, year = {2018}, note = {arXiv:1809.05121}}
@misc{Usui2021tate,
  author = {Usui, Satoshi},
  title  = {Tate--{H}ochschild cohomology rings for eventually periodic {G}orenstein algebras},
  year   = {2021}, note = {arXiv:2107.03326}}         % SUT J. Math. vol/pages at citation-add
@article{BerghJorgensen2013tate,
  author  = {Bergh, Petter Andreas and Jorgensen, David A.},
  title   = {Tate--{H}ochschild homology and cohomology of {F}robenius algebras},
  journal = {J. Noncommut. Geom.}, volume = {7}, number = {4}, pages = {907--937},
  year    = {2013}, doi = {10.4171/JNCG/139}}
```
  and in `registry.py` (`_r(key, bibtex_key, kind, title, annotation, *tags)`):
```python
_r("wang_singular_hh", "Wang2015singular", "foundation",
   "Singular Hochschild cohomology and Gerstenhaber algebra structure",
   "Wang: introduces singular Hochschild cohomology HH_sg^*(A) = Hom_{D_sg(A^e)}(A, A[*]) "
   "(= Tate-Hochschild cohomology) with its Gerstenhaber/BV structure -- the definitional origin.",
   "hochschild", "tate"),
_r("keller_singular_hh", "Keller2018singular", "foundation",
   "Singular Hochschild cohomology via the singularity category",
   "Keller: singular (= Tate) Hochschild cohomology is the Hochschild cohomology of the dg "
   "singularity category -- the singularity-category identification.", "hochschild", "tate"),
_r("usui_tate_periodic", "Usui2021tate", "foundation",
   "Tate-Hochschild cohomology rings for eventually periodic Gorenstein algebras",
   "Usui: for a Gorenstein algebra, eventual periodicity <=> an invertible homogeneous element "
   "of the Tate-Hochschild cohomology ring -- the periodicity certificate.", "hochschild", "tate"),
_r("bergh_jorgensen_tate", "BerghJorgensen2013tate", "algorithm",
   "Tate-Hochschild homology and cohomology of Frobenius algebras",
   "Bergh-Jorgensen: Tate-Hochschild (co)homology via complete resolutions over A^e; agrees "
   "with ordinary HH for n >= d+1 (d = Gor.dim A^e); the symmetric duality dim ĤH^{-j} = dim "
   "ĤH^{j-1}; the explicit quantum-complete-intersection computation.", "hochschild", "tate"),
```
  Verify BibTeX presence (`tests/citations/`) and that every key resolves.
- [ ] **Step 2: Verification page.** Add the P76 subsystem rows:
  - `engine/complete_resolution.py` — `oracle_selfcert` (`d.d=0` across every splice joint +
    in-window acyclicity; multi-vertex un-adapted regression; **the `kZ₃/J²` ν-twist
    discriminating self-cert — `d.d=0`/acyclic hold with the `π`-permuted tags and FAIL with
    identity tags**, M3); `oracle_crossengine` (positive-agreement `ĤH^m == HH^m` for `m ≥
    agrees_from`, the NEW engine vs shipped HH, incl. `kZ₃/J²` `[1,0,0,0,0,1]`).
  - `hochschild/tate.py` — `oracle_literature` (`k[x]/(x^n)` full ring `ĤH^m = n−1`(char∤n)/`n`
    (char∣n) for `n∈{2,3,4}` over `GF(p)`; the QCI `QQ` positive-threshold `[2,1,0,…]` + the
    Bergh-Jorgensen full-ring literature statement; symmetric-Nakayama `ĤH^{m≥1}=1`);
    `oracle_crossengine` (positive-threshold universal; Tate cup positive part ≡ Plan-35 cup;
    `duality` dims ≡ `native` dims on shared symmetric `GF(p)`); `oracle_selfcert` (symmetric
    duality self-consistency `dim ĤH^{-j}=dim ĤH^{j-1}`; the invertible-element identity
    `η·η^{-1}=1` on `k[x]/(x²)`, `periodicity_degree=2`; the Ω-periodicity `is_isomorphic`
    certificate; the not-Gorenstein / non-symmetric-duality / char≤dim / structure-constant
    loud refusals).
  - `qpa` — the `NamesGVars()` guard (FAILS if QPA ships a Tate/singular-HH surface); the
    ordinary `HH^n` agreement anchor crosschecked via the existing QPA `HochschildCohomology`
    bridge (input-level, not a Tate verdict compare).
  - **Honest-scope entries (binding):**
    (a) **The precise threshold is `n ≥ d+1`, `d = Gor.dim(A^e)`** (Bergh-Jorgensen), NOT
    "findim of `A`" (the card wording corrected); self-injective ⇒ `d=0` ⇒ agreement for
    `n ≥ 1`.
    (b) **The symmetric duality `dim ĤH^{-j} = dim ĤH^{j-1}`** is the Bergh-Jorgensen theorem
    (`ν²=id`); the `duality` route ships theorem-backed **dims only** (not an independent
    complete-resolution recomputation), labelled `engine="duality"`.
    (c) **Scope (M4):** self-injective ⇒ always served over `GF(p)` (`A^e` self-injective,
    `d=0`; `D_{A^e}(P_n)` projective so the dual-splice is well-defined; a bounded window even
    when NOT eventually periodic — the `QCI` witness). **Gorenstein NON-self-injective (even if
    eventually periodic) is DEFERRED in v1** — the `D_{A^e}(P_n)`-dual splice CANNOT build it
    (`D_{A^e}(P_n)` is projective iff `A^e` self-injective); it needs a **periodicity-extension**
    complete resolution (Usui), a named future task `P76-followup`. **Not Gorenstein ⇒ LOUD
    refusal** (no complete resolution; Usui: eventually periodic ⇏ Gorenstein). The
    eventual-periodicity CERTIFICATE (`tate_periodicity`) is still computed for every input.
    (d) **The `native` Tate ring (`ĤH^0`, negatives, cup) is `GF(p)` + self-injective only**
    (the minimal-`A^e` engine); over other Domains the surface ships `duality` dims (symmetric,
    degrees `|m| ≥ 1` with the `{0,-1}` pair native-only — M1) or the `positive` threshold
    (degrees `≥ 1` — M2); the full ℤ-graded ring off `GF(p)` and the Gorenstein-non-self-inj
    ring are out of scope (a CS char-0 complete resolution + the periodicity-extension are the
    future extensions).
    (e) **The QCI full-ring pin is char-0 (`QQ`) LITERATURE** (`q` not a root of unity —
    fails over `GF(p)`, where `q` IS a root of unity so **only self-cert (`d.d=0`/acyclic +
    threshold agreement) validates the GF(p) QCI in v1**, MINOR d); over `QQ` the plan pins the
    *positive* threshold part (`ĤH^{n≥1}=[2,1,0,…]`, `pos_dims[0]=None`) + records the
    Bergh-Jorgensen full ring (`ĤH^0=1`, negatives all 0) as a documented literature statement,
    NOT an engine-recomputed value. Convention: quiverlab `xy+q·yx` ≡ BJ `XY−qYX` for generic
    `q` (same HH; MINOR b).
    (f) **`ĤH^0` (and `ĤH^{-1}`) of the symmetric Nakayama are the closed 2-cycle** (M1): the
    duality links them to EACH OTHER (`dim ĤH^0 = dim ĤH^{-1}`), NOT to any positive, so they
    are `native`-only (the `duality` route returns them `None`, never `HH^0`). Only
    `k[x]/(x^n)`'s `ĤH^0` (`= n-1`) is hand-derived (its explicit period-2 resolution).
    (g) **`ĤH^0 ≠ HH^0`** in general (the Tate stable centre) — e.g. `k[x]/(x³)` has `ĤH^0=2`
    vs `HH^0=3`; QCI `ĤH^0=1` vs `HH^0=2`. The block/report carry `HH^0` only under
    `ordinary_pos` (labelled Hochschild), never in a `ĤH^•` slot; `pos_dims[0]` is `null` off
    the native route (M2). Stated, never conflated.
    (h) **The `ν`-twist tag rule uses the VERTEX permutation `π` (`soc P_v ≅ S_{π(v)}`), NOT
    the `d×d` `nakayama_automorphism()` matrix** (which is non-identity even for symmetric
    algebras — live-verified on `kZ₂/J³`). `π = id` iff weakly symmetric; the required
    discriminating witness is **`kZ₃/J²`** (self-injective, non-symmetric, `π = (1 2 3)`), on
    which `d.d=0`/acyclic hold ONLY with `π` (M3). The exact untwisted negative values there
    are engine-computed (non-symmetric ⇒ not the symmetric mirror); the BJ `ν²`-twist duality
    is the literature statement, its twisted-coefficient recompute DEFERRED (needs Plan-52).
  - **Recount the class table** (`tests/release/test_oracle_classes.py` drives the numbers —
    run collection, paste the LIVE counts, re-run to green; no guessed at-authoring number).
- [ ] **Step 3: README.** One features line: "Tate–Hochschild (singular Hochschild) cohomology
  `ĤH^*(A)` in positive AND negative degrees with cup product, via complete resolutions
  (splice the minimal `A^e`-resolution with its Nakayama-twisted dual) for **self-injective**
  algebras; an eventual-periodicity certificate (Usui's invertible homogeneous element);
  agreement with ordinary `HH^*` above the enveloping-algebra Gorenstein dimension — R3 (Wang;
  Keller; Usui; Bergh–Jorgensen)."
- [ ] **Step 4: Full gate:** `... tests/engine tests/hochschild -q -m deep`, `... tests/webapp
  tests/gui tests/hpc -q -m fast`, `... tests/qpa -q -m qpa`, `... tests/release tests/citations
  -q` — all green.
- [ ] **Step 5:** commit
  `docs(verification): P76 R3 Tate-Hochschild oracle rows + honest scope (threshold n>=d+1 corrected, native GF(p)-only ring, QCI char-0 literature, ĤH^0 != HH^0) + citations (Wang/Keller/Usui/Bergh-Jorgensen) + recounted classes`.

---

## Acceptance (Plan-76 definition of done)

1. **The complete-resolution engine, `quiverlab.engine.complete_resolution`:**
   `complete_resolution`, `tate_cohomology_dims`, `tate_homology_dims` build `𝕋` in
   `[-N-1, N]` by splicing the minimal `A^e`-resolution with its Nakayama-twisted `A^e`-dual;
   **`assert_dd_zero` across every splice joint** (incl. `d_0∘d_1`, `d_{-1}∘d_0`) and
   in-window `assert_acyclic` are MANDATORY and green. Multi-vertex is fed the minimal engine
   **un-adapted** (Plan-18 gotcha, regression-pinned). GF(p).
2. **The public surface, `quiverlab.hochschild.tate`:** `tate_hochschild`,
   `tate_cup_products`, `tate_periodicity`, `nakayama_permutation` public + `Algebra`
   delegates. Three routes (`native` GF(p) **self-injective** / `duality` symmetric any-Domain,
   degrees `|m|≥1` / `positive` any-Gorenstein, degrees `≥1`) with honest dispatch and loud
   refusals (not Gorenstein; **native on non-self-injective → DEFERRED refusal, M4**; `duality`
   on `π≠id`; structure-constant; `char ≤ dim`). **Degree-0 honesty (M1/M2):** `pos_dims[0]`/
   `hat_hh0` are `None` off native (`ĤH^0` is native-only, ≠ `HH^0`); the `duality` route
   returns `ĤH^{-1}` as `None` (the `{0,-1}` closed 2-cycle).
3. **HARD literature pins green:** `k[x]/(x^n)` FULL Tate ring `ĤH^m = n−1` (char∤n) / `n`
   (char∣n) for **all** `m∈[-N,N]`, `n∈{2,3,4}` (`GF(32003)` and `GF(2)`/`GF(3)` for the
   `char∣n` cases); the QCI `QQ` positive-threshold `ĤH^{n≥1}=[2,1,0,…]` (`pos_dims[0]=None`);
   the symmetric-Nakayama `ĤH^{m≥1}=1`; the `kZ₃/J²` ν-twist witness (self-inj, `π=(1 2 3)`,
   positive `ĤH^{n≥1}=[1,0,0,0,0,1]`, the `π`-required `d.d=0` discriminating self-cert).
4. **The cross-engine oracles hold:** `ĤH^m == HH^m` for `m ≥ agrees_from` (the NEW engine vs
   the shipped `hochschild_cohomology` — the headline cross-check); the Tate cup positive part
   ≡ the Plan-35 `cup_products` in-window; `duality` dims ≡ `native` dims on shared symmetric
   `GF(p)`.
5. **The self-certs hold:** `assert_dd_zero`/`assert_acyclic`; the symmetric duality
   self-consistency `dim ĤH^{-j} == dim ĤH^{j-1}`; the invertible-element identity
   `η·η^{-1}=1 ∈ ĤH^0` for `k[x]/(x²)` (`periodicity_degree = 2`); the `Ω`-periodicity
   `is_isomorphic` certificate (reuse of `homdims.omega_periodicity`), honest `None` when not
   periodic (the QCI: `period=None`, `complexity=2`, still served as a bounded self-injective
   window).
6. **Honest scope on the verification page (binding):** the corrected threshold `n ≥ d+1`
   (`d = Gor.dim A^e`); the `native` ring is `GF(p)` + **self-injective**-only; the `duality`
   route ships theorem-backed dims (symmetric, degrees `|m|≥1`, the `{0,-1}` pair native-only);
   the QCI full ring is char-0 LITERATURE (root-of-unity over `GF(p)` ⇒ only self-cert there);
   `ĤH^0 ≠ HH^0` (and `pos_dims[0]=None` off native); the `π` vertex-permutation tag rule with
   the `kZ₃/J²` witness; **self-injective-served vs Gorenstein-non-self-injective-DEFERRED vs
   not-Gorenstein-refused** (M4).
7. **One range kind `tate_hochschild`** clickable end-to-end (GUI canvas → block → report) in
   **all four locales (en/es/fr/zh)** with `pick.kind.*` keys and exact key parity, schema v1,
   both runners byte-identical, the grammar parsed in **all three sites**, ONE golden added
   with a documented change-log entry, canonical keys byte-stable; the **positive AND negative
   degree tables** + the period / `agrees_from` lines render.
8. **QPA battery green (`-m qpa`):** the `NamesGVars()` guard (FAILS if QPA ships a Tate
   surface) + the ordinary-`HH^n` agreement anchor on the Nakayama/Dynkin zoo.
9. `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the seven
   honest-scope entries (a)–(g); citations added and BibTeX-verified (Wang / Keller / Usui /
   Bergh–Jorgensen); README line; deep + fast + qpa + release + citations green. No dependency
   on other Wave-4 plans (independent merge to `dev`).

---

## Methodology & assumptions

**Approach.** I (1) read the metaplan P76 card + the R3 record verbatim; (2) fetched and read
**all four references** — the Bergh–Jorgensen PDF page-by-page (the definition, the threshold,
the Frobenius duality, and the QCI computation are transcribed verbatim), and the Wang /
Keller / Usui arXiv abstracts; (3) inventoried the live machinery by import/grep/probe —
confirming the minimal-`A^e` engine (`minimal_resolution`, `AeEngine`,
`minimal_cohomology_dims`, corner tags), the Plan-35 product surface + the Plan-20/21
diagonal, the duality/opposite/`TensorProduct` machinery, the Frobenius/self-injective
certificates, and the shipped `omega_periodicity`/`gorenstein_dimension`/`complexity` — and
that **no** Tate/singular/complete-resolution machinery exists (so it is genuinely new);
(4) **live-verified every numeric pin** in the venv (see below); (5) mirrored the house form
from Plan-65/Plan-63 (record verbatim, reference re-verification with `# PIN`, live-facts
tables, scope gates, Plan-32 markers, three-tier GUI, verification rows, TDD task list).

**Live-verified (venv; shipped engines):** ordinary `HH^*`/`HH_*` of `k[x]/(x^{2,3,4})` over
`QQ`/`GF(2)`/`GF(3)` (`[n, n−1, n−1, …]` for `char∤n`; `[n,n,n,…]` for `char∣n` — the hand
derivation confirmed exactly), of symmetric `kZ₂/J³` (`[3,1,1,…]`) and `kZ₃/J⁴` (`[4,1,1,…]`,
`is_selfinjective`/`is_symmetric` both `True`), and of `QuantumCI(2,2,2)` over `QQ`
(`HH^*=[2,2,1,0,…]`, `is_selfinjective=True`, `is_symmetric=False`); the `A^e`-Betti ranks
(`k[x]/(x^n)` → `[1,1,1,…]`; symmetric Nakayama → constant `[c,c,…]`, both **eventually
periodic**; `QCI` → `[1,2,3,…]`, **complexity 2, NOT periodic**); `complexity`
(`k[x]/(x^n)`/Nakayama `= 1`, `QCI = 2`, `kA₂ = 0`); `gorenstein_dimension` (`k[x]/(x^n)`,
Nakayama `= 0/0`; `kA₂ = 1/1`). From these + the verbatim Bergh–Jorgensen threshold
(`n ≥ d+1`) and symmetric duality (`dim ĤH^{-j} = dim ĤH^{j-1}`) I DERIVED the Tate targets:
`k[x]/(x^n)` FULL closed form `ĤH^m = n−1`/`n` for all `m` (with `ĤH^0` derived rigorously
from the explicit period-2 complete resolution), and the QCI Bergh–Jorgensen values
(`ĤH=[…0,1,2,1,0…]`). I confirmed the i18n locale set is exactly `en/es/fr/zh`, that the
compute grammar is parsed in **three** sites, and that the range-kind grammar (`name:0..N`)
already admits `tate_hochschild:0..N`.

**Could NOT live-verify (labelled to-be-verified-at-implementation).** (1) The **actual
complete resolution and its `ĤH` output** — building the splice IS the plan's work; the
derived targets (threshold + duality + the period-2 `ĤH^0`) are the standing oracles. (2)
**`ĤH^0` of the symmetric Nakayama** (I derived `ĤH^0` only for `k[x]/(x^n)` via the explicit
period-2 resolution; the engine computes the Nakayama `ĤH^0`, gated by the duality
self-consistency). (3) **The QCI negative-degree vanishing over `QQ`** (needs a char-0
complete resolution; the `native` engine is `GF(p)`, where `q` is a root of unity) — pinned as
Bergh–Jorgensen literature + a `QQ` positive-threshold pin. (4) Whether
`TensorProduct(A, A.opposite())` builds `A^e` cleanly enough for the shipped `omega_periodicity`
to apply verbatim (`# PIN` — the engine-level syzygy `is_isomorphic` certificate is the
fallback).

**Fix-round re-verification (2026-08-08, adversarial review NEEDS WORK → all findings VALID).**
Applied document-only, each live-verified in the venv where cheap: **M1** (duality degree-0):
re-derived the BJ reflection `n ↔ -(n+1)` — `ĤH^0`/`ĤH^{-1}` are a closed 2-cycle unlinked to
positives; the duality route now supplies `ĤH^{-j}=HH^{j-1}` for `j≥2` only and returns
`{0,-1}` as `None` (never `HH^0`). **M2** (`pos_dims[0]`): degree 0 off the native route is
`None` (`ĤH^0` native-only, `≠ HH^0`); `ordinary_pos` carries `HH^0` labelled ordinary; the
QCI/positive tests now assert `n≥1`. **M3** (ν-twist): **live-verified `kZ₃/J²` is
self-injective, NON-symmetric, vertex permutation `π = {1:2,2:3,3:1}` (the 3-cycle), `HH^* =
[1,1,0,0,0,0,1]`** over `GF(32003)` and `QQ`; added it as the required discriminating witness;
**live-confirmed `nakayama_automorphism()` returns a non-identity `6×6` matrix even for
symmetric `kZ₂/J³`** — so the tag rule uses the VERTEX permutation `π` (new
`nakayama_permutation` helper), not the algebra matrix. **M4** (Gorenstein-non-self-injective):
`D_{A^e}(P_n)` is projective iff `A^e` self-injective, so the dual-splice cannot build the
non-self-injective Gorenstein case — DEFERRED (`P76-followup` periodicity-extension); v1 native
scope narrowed to self-injective (`kA₂` refusal test added). **MINOR b:** `QuantumCI` is
`xy + q·yx` (live: `x*y + 2*y*x`), `≡` BJ's `XY−qYX` for generic `q`. **MINOR c:** the
`0..N`-only grammar already rejects nonzero `lo` (guard reused) → canonical key `(kind, hi)`,
no duplicate. **MINOR d:** over `GF(p)` the QCI is only self-cert-validated (root of unity).

**Deliberately not checked, and why.** I did not run any test suite (other agents run tests
concurrently — I kept probes light and short, GF(p) at moderate degrees). I did not build the
complete resolution myself (that is the plan's implementation, not a spec probe). I did not
resolve the exact minimal **period** of the symmetric Nakayama as a complex (constant `A^e`
ranks establish periodicity; the period value is engine-detected and cross-checked against the
invertible element — I only hand-pinned `k[x]/(x^n)`'s period 2). I did NOT hand-compute the
**untwisted negative `ĤH^{-j}(kZ₃/J²)` values** (they need the built complete resolution; the
discriminating oracle there is the `π`-required `d.d=0`/acyclic self-cert + non-mirror
assertion, and the exact BJ `ν²`-twist recompute is DEFERRED to Plan-52 twisted coefficients).
I did not verify the published journal volume/pages of Wang / Usui (shipped as
`note={arXiv:…}`, `# PIN` at citation-add) — the arXiv identities and titles ARE verified.

**Why I believe the result is correct.** The two HARD oracle families are independently
grounded and live-checked: the `k[x]/(x^n)` full Tate ring follows from a verbatim-verified
theorem (Bergh–Jorgensen threshold + symmetric duality) applied to live-computed ordinary HH,
with `ĤH^0` derived rigorously from the algebra's *known closed-form* period-2 `A^e`-resolution;
and the headline cross-engine oracle (`ĤH^m ≡ HH^m` for `m ≥ agrees_from`) checks the NEW
splice engine against the *already-tested* shipped HH engine, so a wrong splice cannot pass. The
design reuses only verified live APIs (the minimal-`A^e` engine, the Plan-20/21 diagonal, the
Plan-35 product container, `omega_periodicity`), and every scope boundary becomes a loud typed
refusal (**native = self-injective + `GF(p)`**; Gorenstein-non-self-injective DEFERRED;
not-Gorenstein refused; char-0-only QCI; `π=id`-only duality with the `{0,-1}` pair
native-only; `char ≤ dim`). The one genuine implementation risk — the `π`-permuted
`A^e`-dual and the splice-joint differential — is isolated in Task I1 behind the mandatory
`assert_dd_zero` + `assert_acyclic` self-certs, which catch a wrong tag/sign loudly on the
first `k[x]/(x²)` instance AND on the `kZ₃/J²` witness (where a wrong `π` — plain swap — fails
them, the discriminating test).

## Open design risks (what a critic will likely attack)

1. **The Nakayama-twisted `A^e`-dual and the splice joint.** Getting `D_{A^e}(A e_v ⊗ e_w A) ≅
   A e_{π(w)} ⊗ e_{π(v)} A` (`π` the VERTEX permutation, M3) and the joint `P_0 ↠ A ↪ D(P_0)`
   right (tags, signs, `π`-direction) is the hardest, least-verified piece. *Mitigation:*
   `assert_dd_zero` across every joint + in-window `assert_acyclic` are mandatory and catch a
   wrong construction on `k[x]/(x²)` (symmetric, `π=id`) AND the `kZ₃/J²` witness (`π=(1 2 3)`,
   where a plain swap ignoring `π` FAILS them — the discriminating test); the positive-agreement
   drift gate (`ĤH^m ≡ HH^m`, `m ≥ agrees_from`) catches a wrong positive half; the symmetric
   duality self-consistency catches a wrong negative half on symmetric inputs.
2. **The `native` engine is `GF(p)` while the QCI literature pin is char-0.** The most
   quotable Bergh–Jorgensen computation (QCI) cannot be reproduced by the `GF(p)` native
   engine (root of unity). *Mitigation:* pin the QCI *positive* part over `QQ` (universal
   threshold) + record the full ring as a documented literature value; the `k[x]/(x^n)`
   full-ring pin (which DOES work over `GF(p)`, both char cases) carries the negative-degree
   coverage. A CS char-0 complete resolution is named as the future upgrade, not a P76
   deliverable.
3. **`ĤH^0` and the stable centre.** `ĤH^0 ≠ HH^0` is a real, easily-miscomputed distinction.
   *Mitigation:* the `k[x]/(x^n)` `ĤH^0 = n−1` is hand-derived and pinned; the duality
   self-consistency (`dim ĤH^{-1} = dim ĤH^0`) gates the engine's `ĤH^0` on symmetric inputs.
4. **The periodicity certificate's `A^e`-module iso.** Reusing `omega_periodicity` needs `A`
   as an `A^e`-module; the `TensorProduct(A, A.opposite())` build may trip the path-basis
   guard. *Mitigation:* the engine-level syzygy `is_isomorphic` certificate is the fallback
   (same certificate, engine layer); `complexity==1` is the cheap necessary signal; Usui's
   invertible-element check (via the native cup) is the ring-theoretic cross-check.
5. **The multi-vertex unit-adaptation gotcha.** Feeding `A.unit_adapted()` to the minimal
   engine on multi-vertex input raises `radical_basis: not nilpotent` (live-reproduced).
   *Mitigation:* Task I1's `_engine_form` helper feeds un-adapted for multi-vertex, with a
   regression pin on `kZ₂/J³`.
6. **Cost / window sizing.** The complete resolution builds both halves; a big self-injective
   input with growing ranks (like a large QCI) has an expensive window. *Mitigation:* the kind
   is a bounded-window range kind (`:0..N`) sized like `hh_cohomology`×3, honest `status` on a
   budget cap; the block ships dims, not the (potentially large) cup tables.

## Change log

- **2026-08-08 authoring.** Initial plan (R3 Tate–Hochschild via the complete-resolution
  splice + the three routes + the periodicity certificate + the Tate cup). All four references
  fetched and read; the Bergh–Jorgensen definition/threshold/duality/QCI transcribed verbatim
  from the PDF. Card corrections: the definitional origin is **Wang** (1508.00190), not Keller
  (Keller 1809.05121 verified as the singularity-category identification); the threshold is
  **`n ≥ d+1`, `d = Gor.dim A^e`**, not "findim of `A`"; the scope is refined with the live
  `QCI` witness (self-injective ⇒ always, even when NOT eventually periodic). Live-verified in
  the venv (`GF(p)`/`QQ`): ordinary HH of `k[x]/(x^{2,3,4})`, symmetric `kZ₂/J³`/`kZ₃/J⁴`, and
  `QuantumCI(2,2,2)`; the `A^e`-Betti ranks (periodicity); `complexity`;
  `gorenstein_dimension`; the multi-vertex unit-adaptation gotcha. The `k[x]/(x^n)` full Tate
  ring (`ĤH^m = n−1`/`n` for all `m`, `ĤH^0` derived) and the QCI Bergh–Jorgensen values are
  the standing oracle targets. Still to-be-verified-at-implementation: the actual complete
  resolution output, the symmetric-Nakayama `ĤH^0`, the QCI negatives over `QQ`, and the
  `A^e`-module build for `omega_periodicity`.
- **2026-08-08 fix round (adversarial review NEEDS WORK → all four MAJORs + four MINORs
  adjudicated VALID; document-only, live-verified).**
  - **M1 (duality degree-0):** the BJ duality reflects `n ↔ -(n+1)` (about `n=-½`), so `ĤH^0`
    and `ĤH^{-1}` are a CLOSED 2-cycle unlinked to threshold-known positives. The `duality`
    route now supplies `ĤH^{-j}=HH^{j-1}` for `j≥2` ONLY, returns `{0,-1}` as `None` (never
    `HH^0`); the duality≡native cross-check is scoped to the supplied degrees; scope table +
    routes + honest-scope (f) + a `test_duality_route_degree0_absent_QQ` updated.
  - **M2 (`pos_dims[0]`):** degree 0 off the native route is `None` (`ĤH^0` native-only, `≠
    HH^0`); `ordinary_pos` carries `HH^0` labelled Hochschild; the QCI/positive tests assert
    `n≥1`; the `TateHochschild`/block fields + report labelling made honest.
  - **M3 (ν-twist):** the tag rule uses the VERTEX permutation `π` (`soc P_v≅S_{π(v)}`, new
    `nakayama_permutation` helper), NOT the `d×d` `nakayama_automorphism()` matrix (live: the
    matrix is non-identity even for symmetric `kZ₂/J³`). Added **`kZ₃/J²`** as the required
    discriminating witness (live: self-injective, non-symmetric, `π=(1 2 3)`,
    `HH^*=[1,1,0,0,0,0,1]`) with the `π`-required `d.d=0`/acyclic self-cert; exact untwisted
    negatives engine-computed, the BJ `ν²`-twist recompute DEFERRED (needs Plan-52).
  - **M4 (Gorenstein-non-self-injective):** the `D_{A^e}(P_n)`-dual splice cannot build it
    (`D(P_n)` projective iff `A^e` self-injective) — DEFERRED (`P76-followup`
    periodicity-extension); v1 native scope narrowed to **self-injective**; scope table + Goal
    + honest-scope (c) + a `test_native_non_selfinjective_deferred` (`kA₂`) aligned.
  - **MINORs:** (a) the "self-injective ⇒ always" reading is the plan-writer's expansion, NOT
    the R3 card's words (attribution corrected); (b) `QuantumCI` is `xy+q·yx` (live), `≡` BJ's
    `XY−qYX` for generic `q`; (c) the `0..N`-only grammar already rejects nonzero `lo` (guard
    reused) → canonical key `(kind, hi)`, stated; (d) over `GF(p)` the QCI is only
    self-cert-validated (root of unity). All numbers touched were re-verified in the venv.
- **2026-08-11 implementation (Tasks I1, I2, A1, A2, G1, G2 DELIVERED; A3 DEFERRED).**
  Branch `plan-76-tate-hochschild`; commits 47d540a (engine), a8ddef0 (public surface),
  23bac2a (three tiers + citations + verification). 54 tests.
  - **The construction was DERIVED, and two of its conventions differ from this plan's
    prose.** `T_{-n-1} = Psi(theta^*(P_n))` with `Psi` = the k-dual made a LEFT
    `A^e`-module through the swap anti-automorphism and `theta` = the twist with
    `theta^*(A) = D(A)`. (i) **The twist is `nu^{-1}`, not `nu`**: the identification
    `phi(a)(m) = lambda(am)` satisfies `phi(nu^{-1}(x) a y) = x phi(a) y`, so
    `D(A) = {}_{nu^{-1}}A_1`. (ii) **The tag rule is `(v,w) |-> (pi(w), v)`**, not the
    plan's sketched `(pi w, pi v)`; the two agree exactly when `pi^2 = id`, so every
    symmetric example is blind to the difference. Entries map
    `f_a (x) f_b |-> nu(f_b) (x) f_a` (transposed); the splice joint is the exact solve
    `Z^{(v)} = G^{-1} R^{(v)} G^{-1}`, `R^{(v)}[i][j] = lambda(nu^{-1}(f_j) f_i e_v)`,
    off the shipped `frobenius_form_generic`.
  - **The plan's Task-I1 control test was WRONG AS WRITTEN and is replaced.** It expected
    a wrong `pi` to fail `assert_dd_zero`. Measured: it does NOT. The negative half is
    carried as `A^e`-module maps in AMBIENT coordinates and the composite is taken there
    WITHOUT reference to the tags, so `d.d = 0` is tag-BLIND and a wrong `pi` sails
    through it while reporting wrong dimensions. The shipped battery pins that as a live
    negative result (`test_dd_zero_alone_does_not_discriminate_pi`) and arbitrates with
    **build-time corner typing** (tags from the socle-derived `pi` vs entries from the
    form-derived `nu` -- a real cross-check) and **`assert_acyclic`** (on `kZ3/J2` the
    `pi = identity` control is non-exact in every negative degree). The P75 lesson
    repeating: a self-cert that cannot see the thing it is meant to arbitrate.
  - Live values now pinned: `k[x]/(x^n)` `HHhat^m = n-1`/`n` for all `m in [-5,5]` over
    `32003/2/3`; `HHhat^0 = 2` vs `HH^0 = 3` on `k[x]/(x^3)`; `kZ2/J3` `HHhat = 1` in every
    degree; **`kZ3/J2` `HHhat[-3..3] = [0,0,0,1,1,0,0]`** (positives = `HH^{>=1}`, negatives
    all 0 -- NOT the symmetric mirror).
  - **Plan-doc correction:** `test_not_gorenstein_refused`'s witness (radical-square-zero on
    `1->2->3`) IS Gorenstein (`is_gorenstein() == True`), so it exercises the DEFERRED
    branch, not the not-Gorenstein one. Also: off self-injective input the threshold would
    need `gorenstein_dimension(A^e)`; rather than guess, `agrees_from = None` and the
    positive route claims nothing (so `kA2` on `auto` returns an honest empty report).
  - **Task A3 (the Tate cup) is DEFERRED to a follow-up, by decision.** The `Z`-graded ring
    needs a diagonal approximation on the TWO-SIDED complex -- the Plan-20/21 diagonal,
    which lives on the Chouhy-Solotar PELT resolution, rebuilt for the `GF(p)` corner
    complete resolution and extended past degree 0. That is plan-sized, not task-sized, and
    shipping only the `p, q >= 1` part would be vacuous (in degrees `>= 1` the Tate cochain
    complex IS the minimal resolution's, so "the positive cup matches Plan 35" has no
    content). `TateHochschild.cup` is therefore always `None`, `periodicity_degree` is the
    resolution's certified period cited as Usui's CRITERION rather than an inverse this
    library computed, and honest-scope entry (3) on `docs/verification.md` states it.
    Reconcile at P80.

