# Plan 68: Skew-Gentle Algebras (R32) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the **skew-gentle** world first-class and no-code, driven by the
user-facing **triple `(Q, I, Sp)`**. Deliver: a **recognizer** (`is_skew_gentle_triple`
+ a construction-tagged `recognizers_block` verdict); the **idempotent-split
constructor** `SkewGentleAlgebra(triple, field)` that presents the genuinely-admissible
**split algebra** `kQ̂/Î` (He–Zhou–Zhu / Chen), the internal object on which every
existing engine runs, dimension-certified against the **associated gentle algebra**;
the **special-string re-gluing** module classification (P46 strings/bands on the
associated gentle algebra, special strings re-glued into their `±` two forms — which in
the split model are literally the two split vertices); **support τ-tilting via the
engine** (P45 runs verbatim on the split algebra) with the **geometric orbifold model
as the cross-check oracle**; and the **brick-finite ⇔ representation-finite**
certificate (Demonet–Iyama–Jasso `brick-finite ⇔ τ-tilting-finite` composed with
Garcia–Lavoué's `rep-finite ⇔ brick-finite for skew-gentle, char ≠ 2`). Every
construction is dimension-certified (a NECESSARY check) and the mesh construction is
pinned by the AR / support-τ-tilting counts on a mesh-bearing example (the SUFFICIENCY
oracle); `Sp = ∅` **byte-reduces to the gentle engine** (the strongest self-cert **for the
no-special case that bypasses the split machinery entirely** — it says nothing about the
mesh/split construction, whose oracle is the counts); the char-2 boundary and the char ≠ 2
certificate scope are stated honestly.

**Architecture:** One new top-level package `src/quiverlab/skewgentle/` — a thin
combinatorial + presentation layer over primitives that already exist (no new math
engines), plus a families re-export for GUI/catalog discoverability:

- `skewgentle/triple.py` — `SkewGentleTriple(Q, I, Sp)` (validated dataclass: `Sp ⊆ Q₀`,
  every `i ∈ Sp` gets a fresh special loop, `I` a set of length-2 relation strings on
  `Q`), `is_skew_gentle_triple(...) -> bool` (the primary decidable recognizer), and
  `associated_gentle(triple) -> Algebra` (the genuinely-gentle `A^g = kQ^sp/⟨I ∪ {εᵢ²}⟩`
  with **nilpotent** loops — admissible, quiverlab-native, `dim A^g = dim A` by HZZ
  Lemma 1.5; the P46 string-classification substrate).
- `skewgentle/split.py` — `split_quiver(triple) -> (Quiver, relations)` (Chen §3: double
  each special vertex `i ↦ i⁺, i⁻`, split each ordinary arrow over its endpoints'
  copies, drop the special loop, emit the zero relations from `I` and the
  commutative/mesh relations at special middle vertices) and
  `SkewGentleAlgebra(triple, field=None) -> Algebra` (the split-presented **admissible**
  `kQ̂/Î`; **per-instance dimension-certified** `dim(split) == dim(associated_gentle)`;
  tagged `_family_citations` + a `_skew_gentle_triple` marker so `recognizers_block`
  reports it and `Sp = ∅` short-circuits to the plain gentle presentation).
- `skewgentle/modules.py` — `skew_gentle_indecomposables(triple, max_length, budget)`:
  P46 `enumerate_strings`/`find_bands` on the associated gentle algebra, each walk
  **typed** `(r, s) ∈ {u, p}²` by its special endpoints (Garcia–Lavoué §2.1), then
  **materialised as A-modules on the split algebra** via `Module.from_arrow_action`;
  a special string re-glues into its **two forms** (`X = i⁺` vs `i⁻`), the split-model
  incarnation of the classical char ≠ 2 `k[T]/(T²−1)` two 1-dim modules. Self-certified
  (`check_module` + `is_indecomposable`), count cross-checked against the AR quiver /
  exchange graph on rep-finite instances.
- `skewgentle/certificate.py` — `is_representation_finite(triple) -> bool | None` (the
  brick-finite ⇔ rep-finite certificate: DIJ `τ-tilting-finite ⇔ brick-finite` on the
  split algebra via P45 `exchange_graph(...).is_complete`, upgraded to rep-finite by
  Garcia–Lavoué Thm 3.1 under **char ≠ 2**; cross-checked against the associated
  gentle algebra's admissible-band census; honest `None` / narrowed verdict off scope).
- `skewgentle/block.py` — `skew_gentle_block(triple)` (the no-code compute payload:
  recognizer verdict + the split-quiver shape + dim law + the indecomposable/band
  classification counts + support τ-tilting via `tau_tilting_block(split)` + the
  rep-type certificate), mirroring `strings_block` / `tau_tilting_block`.
- `families/skew_gentle.py` — a thin re-export of `SkewGentleAlgebra` + a `discover.py`
  CATALOG `FamilyInfo` (non-scalar constructor, the `BrauerGraphAlgebra` precedent).

The whole thing composes on: the P38 recognizers (`is_gentle`, `is_special_biserial`,
`_len2_in_ideal`, `reduction_system_of`, `recognizers_block`); the presented backbone
`Quiver.algebra(relations, field, degree_bound)` (which accepts **binomial / mesh
relations** `"p - q"` — the `families/brauer.py` and `families/trivial_extension.py`
idiom); the P46 string subsystem (`enumerate_strings`, `find_bands`, `string_module`
via `Module.from_arrow_action`); the P45 τ-tilting engine (`exchange_graph`,
`tau_tilting_block`, `bricks`, `SupportTauTiltingPair`); P30 `decompose` /
`is_indecomposable`; P23/24 `Module` + `is_isomorphic`; and the block/runner/GUI wiring
of the `strings` and `tau_tilting` scalar kinds.

**Tech Stack:** pure exact combinatorics + exact presentation over `Domain`; **no
floats in `src/`** (AST-gated by `tests/test_no_floats.py`). The split algebra is a
genuine `kQ/I` — its matrices and relations are integer / exact-field data; the special
loops enter only as the vertex splitting, never as a float.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P38 (recognizers), P46 (gentle/strings) and P45 (τ-tilting) are MERGED to `dev`
  and are hard prerequisites.** This plan consumes `invariants/recognizers.py`
  (`is_gentle`/`is_special_biserial`/`_len2_in_ideal`/`recognizers_block`),
  `resolutions_cs/build.py::reduction_system_of`, `strings/walks.py`
  (`enumerate_strings`/`find_bands`/`is_valid_walk`/`letter_source`/`letter_target`),
  `strings/modules.py::string_module`, `tautilting/mutation.py::exchange_graph`,
  `tautilting/block.py::tau_tilting_block`, `tautilting/torsion.py::bricks`, and
  `modules/decompose.py`/`modules/hom.py` exactly as they read on `dev`. Branch
  `plan-68-skew-gentle` off `dev` **after P45+P46 have merged** (both present at
  authoring time — verify with
  `python -c "import quiverlab.strings.walks, quiverlab.tautilting.mutation, quiverlab.invariants.recognizers"`
  before starting). No dependency on any other v1.0.0 subplan — merges independently.
- **Test buckets are auto-assigned by directory** (`tests/conftest.py`):
  `tests/families/` and `tests/modules/` → **deep**; `tests/invariants/`,
  `tests/webapp/`, `tests/gui/`, `tests/hpc/` → **fast**; `tests/qpa/` → **qpa**. Per
  the metaplan surfaces decision the constructor/split/dim/certificate batteries live in
  **`tests/families/test_skew_gentle*.py`** (deep) and the re-gluing / τ-tilting /
  brick-finite batteries in **`tests/modules/test_skew_gentle*.py`** (deep); the
  webapp/GUI block tests in `tests/webapp/` + `tests/gui/` (fast/unmarked). Run new
  tests by path during development; finish each task with a `-m deep` / `-m fast` /
  `-m qpa` spot-run of the touched files.
- **The idempotent form `εᵢ² − εᵢ`, NOT the involution form `εᵢ² − eᵢ`.** HZZ (Def 1.3)
  and Chen use `εᵢ² = εᵢ` (a genuine idempotent); Garcia–Lavoué (Def 2.2, from
  Amiot–Brüstle) use `εᵢ² = eᵢ` (an involution, `k[T]/(T²−1)`). In char ≠ 2 the two are
  isomorphic; in **char 2** they differ (`k[ε]/(ε²−ε) ≅ k×k` is semisimple;
  `k[T]/(T²−1) = k[T]/((T−1)²)` is local). **quiverlab standardises on the idempotent
  form** because `εᵢ, eᵢ−εᵢ` are orthogonal idempotents in **every** characteristic
  (`εᵢ(eᵢ−εᵢ) = εᵢ² − εᵢ = 0`), so the split construction, the dim law, and the module
  classification are **characteristic-free** (Chen's whole point). The char-2 boundary
  bites only the brick-finite ⇔ rep-finite *certificate* (below), never the constructor.
- **The triple cannot be presented directly — it is NOT admissible (VERIFIED).**
  `Quiver.algebra(relations=["e*e - e"])` raises `AdmissibilityError` ("relation has a
  path of length 1: the ideal is not inside the square of the arrows"). `⟨εᵢ² − εᵢ⟩`
  contains `εᵢ ∈ rad \ rad²`, so it is not admissible. **The internal object is always
  the split algebra `kQ̂/Î` (admissible), never a direct presentation of the triple.**
  This is the load-bearing design decision (design question #1).
- **`decompose` / `is_isomorphic` / `is_indecomposable` char caveat is load-bearing.**
  They are rigorous only over char 0 or char > dim M (`modules/decompose.py`). Every
  battery that materialises + decomposes a skew-gentle module, spot-checks
  indecomposability, or runs the exchange-graph BFS defaults to **QQ** (with a
  large-prime `GF(32003)` parity cross-check where cheap). A **char-2** test exists to
  confirm the split construction + dim law are characteristic-free AND that the
  certificate honestly narrows its scope — never to decompose over char 2.
- **Loud refusals, never silent wrong answers** (the `BrauerGraphAlgebra` /
  `TrivialExtension` precedent): an invalid triple (`Sp ⊄ Q₀`; a special vertex already
  carrying a loop; `(Q^sp, I^sp)` not a gentle pair; a length-≠2 relation in `I`), a
  split algebra whose dimension certificate fails, a materialised module that fails
  `check_module`, a certificate asked over char 2 without honest narrowing — each raises
  `QuiverlabError` (or its `AdmissibilityError` subclass) with a `hint=`.
- **Certificate / mesh-sign convention is ARBITRATED, not assumed** (the house
  cup-sign / Brauer-direction precedent). The commutative/mesh relations at split
  vertices are emitted in the sign-free **commutative** form `p − q`; the sign is a
  gauge (rescale one arrow by `±1`) and is immaterial to the algebra's iso type (and
  trivially so in char 2). The **dimension certificate** `dim(split) ==
  dim(associated_gentle)` is a **NECESSARY** check that the relation set is right — an
  over/under-relation changes the dim (exactly the Brauer-graph pattern) — but it is
  **NOT sufficient**: two different relation sets can share the split dim (e.g. a mesh
  identification `p − q` and some wrong monomial both landing on the same count). The
  **sufficiency** oracle is the AR-quiver / support-τ-tilting counts on a
  **mesh-bearing** example (a special vertex in the middle of a length-2 relation forces
  a genuine binomial mesh relation, and the AR / τ-tilting counts pin its structure
  end-to-end — see the mesh oracle in Task 2/3/5). Do not hand-guess the relation set;
  start from Chen §3, let the dim check (necessary) + the AR/τ-tilting counts on the mesh
  example (sufficient) + the `Sp = ∅` byte-reduction + the Chen/QPA oracles decide.
- House conventions: path composition is **left-to-right** (`a*b` = first `a` then `b`,
  `target(a) = source(b)`); HZZ/Chen write `βα` for "first `α` then `β`" (right-to-left)
  — **flip when transcribing** their `G1–G4` / relation displays into quiverlab
  relation strings. Modules are **right** modules by default; `is_isomorphic`/`decompose`
  refuse loudly across sides/algebras. All refusals are `QuiverlabError`.
- Plan-32 markers: dimension/`check_module`/partition certificates + `Sp = ∅`
  byte-reduction + the split-shape self-identities = `oracle_selfcert`; split-model
  indecomposable count ≡ AR/exchange-graph count, geometric rank/count ≡ engine
  τ-tilting, `is_representation_finite` two-route agreement = `oracle_crossengine`; the
  HZZ/Chen/Garcia–Lavoué worked claims (dim law, `n_split = |Q₀|+|Sp|`, selfinjective
  iff, brick-finite ⇔ rep-finite, the `k×k` and small-example τ-tilting counts) =
  `oracle_literature`; QPA comparisons live in `tests/qpa/` (bucket = the class, never
  double-marked).
- **Mid-merge-train counts:** v1.0.0 lands many subplans in overlapping waves, so
  absolute suite counts drift between this plan's authoring and its merge. **Task 8
  recounts the oracle-class table at merge time by running
  `tests/release/test_oracle_classes.py`** (paste the live numbers, never a
  guessed-at-authoring count) and claims only the deltas this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class
  table green, `tests/release/test_oracle_classes.py` green) and adds its citations to
  `citations/references.bib` + `registry.py`. Conventional commits; green at every
  commit.

---

## The Record (R32, verbatim)

> **R32 — Skew-gentle algebras: recognizer, classification, τ-tilting.** [D-scout P10;
> keep] Object: skew-gentle triples (Q, I, Sp) recognizer; idempotent-split reduction to
> the associated gentle algebra, re-gluing special strings; support τ-tilting via the
> orbifold model; brick-finite ⇔ rep-finite certificate. Refs: He–Zhou–Zhu
> arXiv:2004.11136; Chen arXiv:2212.06467; Amiot arXiv:2107.02646; Garcia–Lavoué
> arXiv:2601.01744 (verified real, Jan 2026). Oracles: Sp = ∅ byte-reduces to the gentle
> engine; the 2212.06467 example counts; geometric-vs-engine τ-tilting cross-check.
> Size M.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P68.

---

## Reference re-verification (mandatory — done at authoring, 2026-08-07)

All four references fetched and read (arXiv PDFs, `pdftotext`); nothing taken on trust.

**He–Zhou–Zhu, arXiv:2004.11136** — *A geometric model for the module category of a
skew-gentle algebra* (Ping He, Yu Zhou, Bin Zhu).
- **Def 1.1 (gentle pair `(Q, I)`)** — (G1) ≤ 2 arrows in/out per vertex; (G2) every
  path in `I` has length 2; (G3) for each arrow `α`, ≤ 1 arrow `β` with `βα ∈ I` and
  ≤ 1 with `γα ∉ I`; (G4) dually for `αβ`. *(`βα` is right-to-left; = quiverlab `α*β`.)*
- **Def 1.3 (skew-gentle triple `(Q, Sp, I)`)** — `Sp ⊆ Q₀`; `Q^sp` = `Q` plus a loop
  `εᵢ` at each `i ∈ Sp`; `I^sp = I ∪ {εᵢ² | i ∈ Sp}`; the triple is *skew-gentle* iff
  `(Q^sp, I^sp)` is a gentle pair. The **skew-gentle algebra** is
  `A = kQ^sp/⟨I^sg⟩` where `I^sg = I ∪ {εᵢ² − εᵢ | i ∈ Sp}` — the loops are
  **specialised to idempotents** `εᵢ² = εᵢ`. A gentle algebra is the `Sp = ∅` case.
- **Remark 1.4** — `{eⱼ | j ∈ Q₀ \ Sp} ∪ {eᵢ − εᵢ, εᵢ | i ∈ Sp}` is a complete set of
  primitive orthogonal idempotents. **Each special vertex splits into two** → the split
  quiver has `|Q₀| + |Sp|` vertices.
- **Lemma 1.5** — `kQ^sp/⟨I^sg⟩ ≅ kQ^sp/⟨I^sp⟩` as **vector spaces**; hence
  `dim A = dim(associated gentle A^g)` where `A^g = kQ^sp/⟨I^sp⟩` (nilpotent loops).
  **This is the dimension certificate.**
- **§6 worked example** — quiver on vertices `1..7`, arrows `a,b,c,d,e,f,g,h`, loop `ε₁`
  (`Sp = {1}`), a *non-special* loop `ε₅` (`ε₅² ∈ I`, so vertex 5 is ordinary),
  `I = {ε₅², ab, ba, ed, hc}`. A basic support τ-tilting module has
  `|R| = 8 = |Q₀| + |Sp|` summands → the **rank oracle** `n_split = |Q₀| + |Sp|`.

**Chen, arXiv:2212.06467** — *A characteristic free approach to skew-gentle algebras*
(Yiping Chen, v2 2024).
- **§3 split construction (the binding source for `split.py`)** — split each special
  vertex `i ↦ i⁺, i⁻`; vertices `Q^A₀ = Q₀^{or} ⊔ {i⁺} ⊔ {i⁻}`. For each arrow
  `α : i → j` of `Q`: both ordinary → one arrow `α`; `i` special, `j` ordinary → two
  arrows `⁺α : i⁺→j`, `⁻α : i⁻→j`; `i` ordinary, `j` special → `α⁺ : i→j⁺`, `α⁻ : i→j⁻`;
  both special → four arrows (all copies). The special loop `εᵢ` yields **no arrow**
  (it becomes the vertex splitting). Relations: a length-2 path through an **ordinary**
  middle vertex → **zero relations** (doubled onto the relevant copies); through a
  **special** middle vertex → **commutative (anti-commutative) mesh relations** among
  the copies. Characteristic-free (`εᵢ, eᵢ−εᵢ` orthogonal idempotents in all char);
  char 2 is exactly the case this approach handles that Geiß–de la Peña / Amiot could
  not.
- **Cor 1.2(a)** — `Dsg(A(Q,I)) ≃ Dsg(A(Q,I,Sp))`; `gl.dim A(Q,I,Sp) < ∞ ⇔
  gl.dim A(Q,I) < ∞`.
- **Cor 1.2(c)** — an **indecomposable** skew-gentle `A(Q,I,Sp)` is **selfinjective iff
  `Sp = ∅` and the gentle `A(Q,I)` is selfinjective** (simple or selfinjective
  Nakayama). Clean boolean oracle. **The "indecomposable" (i.e. connected, not an
  algebra direct product) hypothesis is load-bearing (M2):** it is NOT true for
  decomposable algebras. Witness `k×k` (the `Sp = {1}` triple on a one-vertex quiver,
  split = two isolated vertices): it is **semisimple, hence selfinjective, with
  `Sp ≠ ∅`** — a direct counterexample to the *unqualified* "selfinjective ⇒ Sp = ∅".
  It is not a counterexample to the corollary as stated, because `k×k` is **decomposable**
  as an algebra (a product of two copies of `k`), so Cor 1.2(c)'s hypothesis excludes it.
  (VERIFIED live: the split of that triple has `dim == 2`, `is_selfinjective() == True`.)
  Every place this plan cites Cor 1.2(c) must keep the "indecomposable" qualifier, and
  **no test asserts the unqualified form** — the Task-2 selfinjective test uses only
  connected (indecomposable) algebras.
- **K-theory** — `Kᵢ(A(Q,I,Sp)) ≅ Kᵢ(A(Q,I)) ⊕ (|Sp| copies of Kᵢ(k))`; at `i = 0`,
  `K₀` rank `= |Q₀| + |Sp|` (agrees with the split-vertex count).
- **Thm 1.1(b)** — `A/AeA ≅ A(Q,I)` (the *plain* gentle, loops dropped), `e` = sum of
  ordinary-vertex idempotents. **Note:** this recollement quotient `A(Q,I)` differs from
  the associated gentle `A^g = kQ^sp/⟨I^sp⟩` (loops kept, nilpotent); we use `A^g` for
  the string classification and `A(Q,I)` only for the homological oracles.
- **`# PIN` on "the 2212.06467 example counts":** Chen's concrete, checkable content is
  **homological**, not enumerative (selfinjective-iff, `K₀` rank, gl.dim-finite-iff,
  Gorensteinness — Cor 1.2/1.3). We wire these as the "2212.06467 oracles"; the numeric
  **τ-tilting counts** are pinned from the **geometric model (HZZ §6 rank + small hand
  cases)**, attribution flagged on the verification page.

**Amiot, arXiv:2107.02646** — *Indecomposable objects in the derived category of a
skew-gentle algebra using orbifolds* (Claire Amiot; ICRA 2020 proceedings). Skew-gentle
algebras as `ℤ₂`-skew-group algebras of gentle algebras (char ≠ 2); a complete
description of indecomposable derived objects via curves on an **orbifold surface** and
its double cover. The citation backing the **orbifold/geometric model** that the
τ-tilting geometric cross-check rests on (module-category orbifold model is HZZ;
derived is Amiot).

**Garcia–Lavoué, arXiv:2601.01744** (v1, 5 Jan 2026) — *Brick-finite skew-gentle
algebras are representation-finite* (Monica Garcia, Léa Lavoué).
- **Main Thm 3.1** — *"Let `A ≅ kQ^sp/I^sp` be a skew-gentle algebra over a field `k`
  with `char(k) ≠ 2`. Then `A` is brick-finite if and only if it is
  representation-finite,"* generalising Plamondon's gentle result [Pla19]. **char ≠ 2**
  is the hypothesis → the certificate's hard scope.
- **Def 2.2 / §2.1 module classification (Table 1)** — indecomposables of `A` come from
  *admissible* strings/bands of `Q^sp`, each typed `(r, s) ∈ {u, p}²` (`p` = endpoint at
  a special vertex, not via the special loop). Associated algebras `A_x`: `(u,u) → k`;
  `(u,p),(p,u) → k[T]/(T²−1)`; `(p,p) → k⟨T,S⟩/(T²−1, S²−1)` with swap `ι`; band →
  `k[T,T⁻¹]`, `ι(T) = T⁻¹`. Surjection `Adm^(k)(Q^sp) → ind(A)` ([CB89, BTCB24]),
  `Mₓ(X) ≅ Mₓ′(X′) ⇔ (x′,X′) = (x⁻¹, X^ι)`. In the split model the `±` forms are the two
  split vertices — **characteristic-free**, sidestepping the char ≠ 2 `k[T]/(T²−1)` split.

**Composed certificate theorem (used by `certificate.py`):**
`brick-finite ⇔ τ-tilting-finite` (Demonet–Iyama–Jasso, any f.d. algebra) **∘**
`rep-finite ⇔ brick-finite` (Garcia–Lavoué Thm 3.1, skew-gentle, char ≠ 2) ⟹
for a skew-gentle `A` over char ≠ 2: `rep-finite ⇔ τ-tilting-finite ⇔
exchange_graph(split).is_complete`.

---

## Design decisions (the seven questions, settled)

1. **Non-admissibility → build the split.** The triple is the user input; the internal
   object is the admissible split algebra `kQ̂/Î` (Chen §3), on which every engine runs.
   VERIFIED: a direct `⟨εᵢ²−εᵢ⟩` presentation raises `AdmissibilityError`.
2. **Recognizer scope = the triple (primary) + construction-tag (secondary).**
   `is_skew_gentle_triple(Q, I, Sp)` is decidable: it reduces to `is_gentle(A^g)` on the
   associated gentle algebra plus the structural checks (`Sp ⊆ Q₀`, each special vertex
   gets a fresh loop, `I` length-2). Recognizing an **arbitrary presented algebra** as
   "isomorphic to some skew-gentle algebra" is the iso-problem — **not attempted**;
   `recognizers_block` reports `is_skew_gentle` only for algebras carrying the
   `SkewGentleAlgebra` construction marker (the `BrauerGraphAlgebra` honest-scope
   precedent).
3. **Associated gentle algebra + re-gluing.** `A^g = kQ^sp/⟨I ∪ {εᵢ²}⟩` (nilpotent
   loops) is genuinely gentle and admissible; P46 enumerates its strings/bands; special
   strings (typed `p`) re-glue into their two forms = the two split vertices `i⁺, i⁻`.
   Modules are materialised on the **split algebra** (that is where they are honest
   A-modules).
4. **τ-tilting = engine-primary, geometric-as-oracle.** The split algebra is a presented
   f.d. algebra → `exchange_graph`/`tau_tilting_block` (P45) run unchanged. The orbifold
   model (HZZ/Amiot) is the CROSS-CHECK oracle: the rank `n_split = |Q₀|+|Sp|` and small
   `|sτ-tilt|` counts (`# PIN`ned).
5. **brick-finite ⇔ rep-finite certificate.** DIJ (general) ∘ Garcia–Lavoué Thm 3.1
   (char ≠ 2). `is_representation_finite(triple)` = `exchange_graph(split).is_complete`,
   char ≠ 2-scoped for the rep-finite reading; cross-checked against the `A^g`
   admissible-band census.
6. **Oracles.** `Sp = ∅` byte-reduces (strongest self-cert — but only for the no-special
   case that bypasses the split machinery); dim law `dim(split) = dim(A^g)` (NECESSARY,
   not sufficient); the mesh-example AR / τ-tilting counts (the SUFFICIENCY oracle that the
   split relations are right); rank `n_split = |Q₀|+|Sp|`; Chen homological claims
   (selfinjective-iff, gl.dim-iff, `K₀` rank); geometric-vs-engine τ-tilting (rank +
   small counts); module `check_module`/`is_indecomposable` + count ≡ AR.
7. **GUI.** A `skew_gentle` compute kind + a skew-gentle **preset** (triple input:
   quiver + relations + special-loop picks); both runners byte-identical; i18n ×4; one
   golden; canonical keys; layout picker; a `recognizers_block` `is_skew_gentle` row.

---

### Task 1: `triple.py` — the validated triple, the recognizer, the associated gentle algebra

The user-facing entry. A `SkewGentleTriple` names `(Q, I, Sp)`; validation reduces
skew-gentle-ness to gentleness of the associated gentle algebra; `associated_gentle`
materialises that gentle algebra (the P46 substrate + the dim certificate reference).

**Files:**
- Create: `src/quiverlab/skewgentle/__init__.py`, `src/quiverlab/skewgentle/triple.py`
- Test: `tests/families/test_skew_gentle_triple.py`

**Interfaces:**
- Consumes: `combinat/quiver.py::Quiver` (`.vertices`, `.arrows` `{name: (src, tgt)}`,
  `.source`/`.target`), `Quiver.algebra(relations, field, degree_bound)`,
  `invariants/recognizers.py::is_gentle`, `resolutions_cs/build.py::reduction_system_of`,
  `errors.QuiverlabError`.
- Produces:
  ```python
  @dataclass(frozen=True)
  class SkewGentleTriple:
      quiver: Quiver          # the ORIGINAL Q (no special loops)
      relations: tuple        # length-2 relation strings on Q (the ideal I generators)
      special: frozenset      # Sp subset of Q_0 ; each gets a fresh special loop
      loop_names: dict        # {i: "eps_<i>"} the chosen fresh loop-name per special i
      def validate(self) -> None
          # Sp subset Q_0 ; no vertex in Sp already carries a loop in Q ; every string in
          # `relations` is a length-2 path in Q ; the associated gentle pair is gentle.
          # Loud QuiverlabError with a hint on each failure.
      def q_sp(self) -> Quiver          # Q plus one loop eps_i at each special i
      def i_sp(self) -> tuple           # relations + {"eps_i*eps_i"} (nilpotent form)

  def associated_gentle(triple, field=None) -> Algebra
      # A^g = kQ^sp / <I ∪ {eps_i^2}> -- a genuine gentle algebra (admissible), the P46
      # string-classification substrate. dim A^g == dim(skew-gentle A) (HZZ Lemma 1.5).

  def is_skew_gentle_triple(quiver, relations, special) -> bool
      # True iff (quiver, relations, special) is a valid skew-gentle triple: builds the
      # SkewGentleTriple, runs validate() inside a try (returns False on QuiverlabError,
      # never raises), i.e. iff the associated gentle pair is gentle + the structural
      # loop/Sp conditions hold. The PRIMARY decidable recognizer.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_skew_gentle_triple.py
"""The skew-gentle triple, the recognizer, and the associated gentle algebra (Plan 68).
Self-cert: a valid triple validates; the associated gentle algebra IS gentle; dim law
is set up (checked against the split in Task 2). Literature: HZZ Def 1.1/1.3 -- the
recognizer accepts the gentle-pair triples and refuses non-gentle ones."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import is_gentle
from quiverlab.skewgentle.triple import (SkewGentleTriple, associated_gentle,
                                         is_skew_gentle_triple)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _Q_arrow():                      # Q: 1 --a--> 2
    return Quiver([1, 2], {"a": (1, 2)})


def test_valid_triple_special_vertex_2():
    assert is_skew_gentle_triple(_Q_arrow(), relations=[], special={2})


@selfcert
def test_associated_gentle_is_gentle_and_admissible():
    t = SkewGentleTriple.make(_Q_arrow(), relations=[], special={2})
    Ag = associated_gentle(t, field=QQ)
    assert is_gentle(Ag)                              # A^g genuinely gentle (nilpotent loop)
    assert Ag.dim == 5                                # e1,e2,a,eps,a*eps  (HZZ Lemma 1.5)


@lit
def test_sp_empty_is_the_plain_gentle():
    # Sp = empty: the associated gentle algebra IS the plain gentle A(Q,I).
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    t = SkewGentleTriple.make(Q, relations=["a*b"], special=set())
    Ag = associated_gentle(t, field=QQ)
    plain = Q.algebra(relations=["a*b"], field=QQ)
    assert Ag.dim == plain.dim == 5


@selfcert
def test_special_vertex_outside_Q0_refused():
    with pytest.raises(QuiverlabError):
        SkewGentleTriple.make(_Q_arrow(), relations=[], special={99}).validate()


@selfcert
def test_non_gentle_associated_pair_refused():
    # a vertex with 3 arrows out cannot be a gentle pair -> not a skew-gentle triple.
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3), "c": (1, 4)})
    assert is_skew_gentle_triple(Q, relations=[], special=set()) is False


@selfcert
def test_length3_relation_refused():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)})
    with pytest.raises(QuiverlabError):
        SkewGentleTriple.make(Q, relations=["a*b*c"], special=set()).validate()
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.skewgentle.triple`

- [ ] **Step 3: Implement** `SkewGentleTriple` (+ a `make(...)` classmethod that picks
  fresh `loop_names` and freezes), `q_sp`/`i_sp`, `validate` (structural checks then
  `is_gentle(associated_gentle(self))`), `associated_gentle`, and `is_skew_gentle_triple`.

**Adjust to reality (Task 1):**
- **Loop-name freshness:** pick `eps_<i>` and guard against collision with an existing
  arrow name; a special vertex that already has a loop in `Q` is refused (HZZ Lemma 1.6:
  at most one loop per vertex). Store `loop_names` in the frozen dataclass so `split.py`
  and `modules.py` read the same names.
- **`i_sp` uses `"eps_i*eps_i"`** (the nilpotent square, admissible) — NOT `"... - eps_i"`
  (that is the non-admissible skew form and would raise). `associated_gentle` therefore
  presents cleanly; if it raises `AdmissibilityError`, the user's `I` was wrong, surface
  it with a `hint`.
- **`degree_bound`:** pass a generous bound (e.g. `2·(|Q_0| + |arrows| + |Sp|) + 2`) to
  `Quiver.algebra` so the Gröbner route terminates on longer gentle strings; mirror
  `families/preprojective.py`'s auto-bound if a tighter one is wanted.
- The recognizer is **convention-free and unconditional** — it must never raise; wrap
  `validate()` and return `False` on any `QuiverlabError` (the `is_string`/`is_gentle`
  boolean-recognizer contract).

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/skewgentle/__init__.py src/quiverlab/skewgentle/triple.py \
        tests/families/test_skew_gentle_triple.py
git commit -m "feat(skewgentle): SkewGentleTriple + is_skew_gentle_triple recognizer + associated_gentle (HZZ Def 1.3; dim law substrate)"
```

---

### Task 2: `split.py` — the idempotent-split algebra `kQ̂/Î` (THE CONSTRUCTOR)

The mathematical crux. Chen §3: double the special vertices, split the arrows, drop the
special loops, emit zero + mesh relations. Present via `Quiver.algebra`; **certify
`dim(split) == dim(associated_gentle)`** (HZZ Lemma 1.5); `Sp = ∅` short-circuits to the
plain gentle presentation (byte-identical).

**Files:**
- Create: `src/quiverlab/skewgentle/split.py`
- Modify: `src/quiverlab/skewgentle/__init__.py` (export `SkewGentleAlgebra`,
  `split_quiver`)
- Test: `tests/families/test_skew_gentle_split.py`

**Interfaces:**
- Consumes: `triple.py` (`SkewGentleTriple`, `associated_gentle`), `Quiver`,
  `Quiver.algebra(relations, field, degree_bound)` (accepts binomial `"p - q"` mesh
  relations — the `families/brauer.py` `_present` idiom), `Algebra.dim`,
  `reduction_system_of` / `_len2_in_ideal` (to read which length-2 paths of `Q` are in
  `I` and whether the middle vertex is special).
- Produces:
  ```python
  def split_quiver(triple) -> tuple[Quiver, list[str], dict]
      # (Q_hat, relations, meta). Q_hat vertices are ALL STRINGS: ordinary v -> str(v)
      # (kept), special i -> "{i}+", "{i}-" (so P45's sorted(pair.support) never mixes
      # int/str -- see H1). Arrows: each ordinary arrow a: s->t split over S(s) x S(t)
      # with S(x) = {""} (ordinary) or {"+","-"} (special) -- name a copy with an
      # IDENTIFIER-SAFE scheme `a` / `ap` / `pa` / `pam` etc. (never `a+`/`+a`: arrow
      # names must be identifiers -- W3); the special loop yields NO arrow. relations:
      # zero relations from I (doubled) + commutative mesh relations at special middle
      # vertices. meta records the vertex/arrow split maps for modules.py + block.py.
  def SkewGentleAlgebra(triple=None, *, quiver=None, relations=(), special=(),
                        field=None) -> Algebra
      # Build the split-presented admissible kQ_hat/I_hat. Accept either a
      # SkewGentleTriple or the (quiver, relations, special) triple pieces (the GUI path).
      # VALIDATES the triple first (loud). CERTIFIES dim == dim(associated_gentle) (loud
      # QuiverlabError otherwise). Tags A._skew_gentle_triple = triple and
      # A._family_citations. Sp = empty short-circuits: return quiver.algebra(relations,
      # field) UNCHANGED (byte-identical to the plain gentle algebra).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_skew_gentle_split.py
"""The idempotent-split skew-gentle algebra (Plan 68 / Chen 2212.06467 sec 3).
Self-cert: dim(split) == dim(associated gentle) (HZZ Lemma 1.5); n_split == |Q_0|+|Sp|
(HZZ Rmk 1.4 / K_0). Literature: Sp = empty BYTE-REDUCES to the plain gentle algebra;
Chen Cor 1.2(c) an INDECOMPOSABLE skew-gentle algebra is selfinjective iff Sp = empty and
the gentle A(Q,I) is selfinjective (the tests use connected algebras only). Loud: a broken
triple; a dimension-certificate failure."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import GF, QQ
from quiverlab.skewgentle.split import SkewGentleAlgebra, split_quiver
from quiverlab.skewgentle.triple import SkewGentleTriple, associated_gentle

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _t_arrow(special={2}):            # Q: 1 --a--> 2 , Sp = {2}
    return SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}),
                                 relations=[], special=special)


@selfcert
def test_dim_law_and_rank():
    t = _t_arrow()
    A = SkewGentleAlgebra(t, field=QQ)
    assert A.dim == associated_gentle(t, field=QQ).dim == 5      # HZZ Lemma 1.5
    Qhat, _rels, _meta = split_quiver(t)
    assert len(list(Qhat.vertices)) == 2 + 1                     # |Q_0| + |Sp| = 3
    assert A.dim == 5 and len(list(A.quiver.vertices)) == 3


@selfcert
@pytest.mark.parametrize("field", [QQ, GF(32003), GF(2)])       # GF(2): char-free!
def test_split_is_characteristic_free(field):
    A = SkewGentleAlgebra(_t_arrow(), field=field)
    assert A.dim == 5                                            # even in char 2


# --- the MESH oracle: a special vertex in the MIDDLE of a length-2 relation (H2a) ------
# Q: 1 --a--> 2 --b--> 3 , Sp = {2}. IMPORTANT: at a special MIDDLE vertex gentleness
# FORCES a*b in I (VERIFIED: a*b NOT in I gives is_gentle == False). So the mesh-bearing
# triple has a*b in I -- Chen sec 3 lifts that relation to a COMMUTATIVE mesh (identify
# the two composites ap*pb == am*mb) rather than a zero-split. This is where the dim law
# is only NECESSARY: a zero-split has dim 8, the mesh has dim 9 == dim A^g, and the
# sufficiency is the AR / tau-tilting counts below. (All three counts VERIFIED live on a
# hand-built dim-certified split with all-string vertex labels ['1','2+','2-','3'].)
def _t_mesh():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    return SkewGentleTriple.make(Q, relations=["a*b"], special={2})   # a*b IN I


@selfcert
def test_mesh_dim_law_is_necessary_not_sufficient():
    t = _t_mesh()
    A = SkewGentleAlgebra(t, field=QQ)
    assert A.dim == associated_gentle(t, field=QQ).dim == 9      # HZZ Lemma 1.5 (mesh)
    # the mesh split has 4 vertices (|Q_0| + |Sp| = 3 + 1) and identifies the two
    # composites through 2+ / 2- (a zero-split would drop the dim to 8).
    assert len(list(A.quiver.vertices)) == 4


@lit
def test_mesh_gentleness_forces_relation():
    # a*b NOT in I at a special middle vertex is NOT a gentle pair -> not a valid triple.
    from quiverlab.skewgentle.triple import is_skew_gentle_triple
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    assert is_skew_gentle_triple(Q, relations=["a*b"], special={2}) is True
    assert is_skew_gentle_triple(Q, relations=[], special={2}) is False    # VERIFIED


@xeng
def test_mesh_ar_and_tau_tilting_counts():
    # the SUFFICIENCY oracle for the mesh construction (dim alone is not enough).
    # VERIFIED live on the dim-certified hand-built split: 11 AR vertices, 46 s-tau-tilt
    # pairs. Confirm at implementation that SkewGentleAlgebra(_t_mesh()) reproduces these
    # (it is the same algebra up to iso -- dim-certified + the forced mesh relation).
    from quiverlab.modules.ar import knit_ar_quiver
    from quiverlab.tautilting.mutation import exchange_graph
    A = SkewGentleAlgebra(_t_mesh(), field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete and len(ar.vertices) == 11            # VERIFIED live
    eg = exchange_graph(A)
    assert eg.is_complete and len(eg.vertices) == 46            # VERIFIED live


@lit
def test_sp_empty_byte_reduces_to_plain_gentle():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    A = SkewGentleAlgebra(quiver=Q, relations=["a*b"], special=(), field=QQ)
    plain = Q.algebra(relations=["a*b"], field=QQ)
    # byte-identical presentation: same quiver, same relations, same Cartan + dim.
    assert A.dim == plain.dim
    assert A.cartan_matrix() == plain.cartan_matrix()
    assert sorted(A.quiver.arrows) == sorted(plain.quiver.arrows)


@lit
def test_chen_selfinjective_iff():
    # Chen Cor 1.2(c) (INDECOMPOSABLE hypothesis): an indecomposable skew-gentle algebra
    # is selfinjective iff Sp = empty AND the gentle A(Q,I) is selfinjective. BOTH algebras
    # below are connected (indecomposable), so the corollary applies as stated -- this is
    # NOT the unqualified form (k x k, a DECOMPOSABLE selfinjective with Sp != empty, is
    # excluded by the hypothesis; see test_kxk_* / the Cor 1.2(c) reference note).
    # the 2-cycle gentle kQ/(ab,ba) is selfinjective; Sp = empty keeps it so.
    Q2 = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)})
    self_inj = SkewGentleAlgebra(quiver=Q2, relations=["a*b", "b*a"], special=(),
                                 field=QQ)
    assert self_inj.is_selfinjective() is True
    # adding a special vertex breaks selfinjectivity; the split stays connected
    # (indecomposable) -- the two copies of vertex 1 are joined through vertex 2.
    not_self = SkewGentleAlgebra(quiver=Q2, relations=["a*b", "b*a"], special={1},
                                 field=QQ)
    assert not_self.is_selfinjective() is False


@selfcert
def test_broken_triple_refused():
    with pytest.raises(QuiverlabError):
        SkewGentleAlgebra(quiver=Quiver([1], {}), relations=[], special={99}, field=QQ)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.skewgentle.split`

- [ ] **Step 3: Implement** `split_quiver` (the Chen §3 construction) + `SkewGentleAlgebra`
  (validate → `Sp = ∅` short-circuit → present → dim-certify → tag). Mirror
  `families/brauer.py::_present` (quiver build + binomial relations + `Quiver.algebra` +
  dim certificate) and `families/trivial_extension.py` (the presented-vs-certificate
  idiom).

```python
# src/quiverlab/skewgentle/split.py  (skeleton -- the crux is _split_arrows + _relations)
"""The idempotent-split skew-gentle algebra kQ_hat/I_hat (Plan 68, Chen 2212.06467 sec 3;
HZZ 2004.11136 Def 1.3 / Rmk 1.4). The triple's ideal <eps^2 - eps> is NON-admissible,
so we present the split (basic, idempotents refined) algebra instead: each special
vertex splits into two, each ordinary arrow splits over its endpoints' copies, the
special loops become the vertex splitting, and I becomes zero + commutative mesh
relations. Per-instance dimension-certified dim == dim(associated gentle) (HZZ Lemma
1.5). Characteristic-free. Float-free / exact."""
from __future__ import annotations

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.skewgentle.triple import SkewGentleTriple, associated_gentle


def _copies(v, special):
    # ALL vertex labels are STRINGS. P45's exchange_graph does `sorted(pair.support)`
    # (tautilting/mutation.py::_label, also line 164), which raises
    # `TypeError: '<' not supported between instances of 'str' and 'int'` on a mix of
    # int and str labels (VERIFIED live -- the hand-built split with vertices
    # [1,'2+','2-',3] crashes exchange_graph at exactly this call). So an ordinary
    # vertex becomes `str(v)`, a special one `"{v}+"`/`"{v}-"`. Quiver accepts "2+" as
    # a VERTEX label (verified); it is only ARROW names that must be identifiers (W3).
    return (f"{v}+", f"{v}-") if v in special else (str(v),)   # S(v) -- all str


def split_quiver(triple):
    Q, Sp = triple.quiver, triple.special
    verts = []
    for v in Q.vertices:
        verts.extend(_copies(v, Sp))
    arrows, arrow_split = {}, {}          # arrow_split[a] = list of (name, s', t')
    for a, (s, t) in Q.arrows.items():
        pieces = []
        for s2 in _copies(s, Sp):
            for t2 in _copies(t, Sp):
                nm = _copy_name(a, s, t, s2, t2, Sp)      # a / ap / pa / pam ... (identifier-safe)
                arrows[nm] = (s2, t2)
                pieces.append((nm, s2, t2))
        arrow_split[a] = pieces
    Qhat = Quiver(verts, arrows)
    rels = _relations(triple, Qhat, arrow_split)          # zero + mesh (see note)
    meta = {"arrow_split": arrow_split, "special": sorted(Sp, key=repr)}
    return Qhat, rels, meta


def SkewGentleAlgebra(triple=None, *, quiver=None, relations=(), special=(), field=None):
    if triple is None:
        triple = SkewGentleTriple.make(quiver, relations, special)
    triple.validate()                                      # loud on a broken triple
    if not triple.special:                                # Sp = empty: plain gentle
        A = triple.quiver.algebra(relations=list(triple.relations), field=field)
        A._skew_gentle_triple = triple
        A._family_citations = ("he_zhou_zhu", "chen_skew_gentle", "assem_book")
        return A
    Qhat, rels, meta = split_quiver(triple)
    Ag = associated_gentle(triple, field=field)
    bound = getattr(Ag, "loewy_length", lambda: 4)() + 2
    A = Qhat.algebra(relations=rels, field=field, degree_bound=bound)
    if A.dim != Ag.dim:                                   # HZZ Lemma 1.5 -- NECESSARY check
        raise QuiverlabError(
            f"SkewGentleAlgebra: split dim {A.dim} != associated-gentle dim {Ag.dim} "
            "(HZZ Lemma 1.5) -- the split arrows/relations are wrong",
            hint="check the Chen sec-3 mesh relations at special middle vertices")
    A._skew_gentle_triple = triple
    A._skew_gentle_meta = meta
    A._family_citations = ("he_zhou_zhu", "chen_skew_gentle", "amiot_skew_gentle",
                           "assem_book")
    return A
```

**Adjust to reality (Task 2) — the split arrows/relations are the whole game:**
- **Vertex labels are ALL strings — the P45 sort constraint (H1, load-bearing).** P45's
  `exchange_graph` labels each pair with `sorted(pair.support)`
  (`tautilting/mutation.py::_label`, and again at `mutation.py:164`
  `"support": sorted(pair.support)`). A support set that mixes `int` and `str` raises
  `TypeError: '<' not supported between instances of 'str' and 'int'` — **VERIFIED live**
  on the hand-built split (vertices `[1, '2+', '2-', 3]` crash `exchange_graph` at exactly
  this call, while all-string `['1', '2+', '2-', '3']` runs clean and returns 46 pairs).
  So `_copies` **stringifies every vertex**: ordinary `v → str(v)`, special `i → "{i}+"`,
  `"{i}-"`. Quiver accepts `"2+"` as a **vertex** label (verified); the identifier
  restriction is on **arrow** names only (next bullet).
- **Arrow copy names + directions** are a CONVENTION checked (necessarily, not
  sufficiently) by the dim certificate + pinned by the AR/τ-tilting counts on the mesh
  example. Chen's `⁺α` / `α⁻` decorations: prefix = source copy, suffix = target copy.
  Emit **identifier-safe** deterministic names (`f"{a}"`, `f"{a}p"`, `f"p{a}"`,
  `f"p{a}m"` … — never `a+`/`+a`, which `Quiver` rejects with "arrow names must be
  identifiers") that never collide; keep `arrow_split` so `modules.py` and `block.py`
  read the same map. If the certificate fails, the copy set or a relation is wrong —
  **fix the construction, never weaken the certificate.**
- **The relation set (`_relations`) is the mathematical crux — start from Chen §3, let
  the dim check + the mesh-example counts arbitrate.** For each length-2 path `a*b` of `Q`
  (with `_len2_in_ideal` telling you whether it is in `I`) and each way it lifts to the
  split copies:
  - **middle vertex ordinary** and `a*b ∈ I` → **zero relation** on each lifted copy
    (monomial `"a_copy*b_copy"`);
  - **middle vertex special** (`i ∈ Sp`) → **the length-2 path `a*b` is necessarily
    `∈ I`.** (Gentleness of `(Q^sp, I^sp)` FORCES this: at a special vertex `i` the loop
    `εᵢ` occupies one in-slot and one out-slot, so `i` has at most one other in-arrow `a`
    and one other out-arrow `b`; among the out-arrows `{b, εᵢ}` of `a`, `a*εᵢ ∉ I`, so G3
    ["at most one arrow β with `a*β ∉ I`"] forces `a*b ∈ I`. **VERIFIED live:**
    `1→a→2→b→3, Sp={2}` with `a*b ∉ I` gives `is_gentle == False`; with `a*b ∈ I` it is
    gentle, `dim A^g == 9`.) Chen's construction lifts this relation NOT to a zero-split
    but to a **commutative mesh relation** (binomial `"p - q"`) that IDENTIFIES the two
    composites through `i⁺` and `i⁻`: for `a*b ∈ I`, `ap*pb - am*mb` (the only valid
    lifts; `ap*mb`/`am*pb` do not compose). Use the sign-free commutative form (the sign
    is a gauge, trivially so in char 2). **VERIFIED live:** a zero-split (`ap*pb = am*mb =
    0`) gives `dim == 8` (wrong, under by one); the mesh identification gives `dim == 9 ==
    dim A^g` (right). Chen's three special-middle cases (both/one/neither of the outer
    vertices special) fix which copies pair up — transcribe them, then let the dim check +
    the mesh-example AR/τ-tilting counts decide.
  There is genuine subtlety here; budget iteration. The dim certificate is a **necessary**
  check (it fires on any over/under-relation) but **not sufficient** (a wrong relation set
  can share the dim); the `Sp = ∅` byte-reduction pins the no-special path exactly, and
  Task 3's `is_indecomposable` materialisation + Task 5's τ-tilting counts on the mesh
  example are the **sufficiency** oracle that pins the mesh structure end-to-end.
- **`Sp = ∅` short-circuit is mandatory for the byte-reduction oracle** — do NOT route
  an empty-`Sp` triple through `split_quiver` (it would relabel nothing but might reorder
  arrows); return `quiver.algebra(relations, field)` verbatim so the presentation is
  byte-identical to the plain gentle algebra.
- **`degree_bound`:** derive from the associated gentle algebra's Loewy length (the split
  has the same dim, so the same-ish path lengths); widen once if `Quiver.algebra` raises
  `NotFiniteDimensionalError` (a genuine skew-gentle algebra IS finite-dimensional, so
  that signals a wiring bug, not an honest infinite case).
- `families/skew_gentle.py` re-export + `discover.py` CATALOG entry: append a
  `FamilyInfo(...)` tuple to the module-level `CATALOG` in `families/discover.py`,
  mirroring `BrauerGraphAlgebra` **exactly** — `route="general"`, its citation tuple, and
  a `summary` string that ANNOTATES the non-scalar nature in prose (BrauerGraphAlgebra's
  reads `"... Non-scalar constructor: not offered by the scalar form-builder."`). **There
  is no `_iter_families` and no skip-set in `discover.py`** (VERIFIED — the file has only
  the `FamilyInfo` dataclass, the `CATALOG` tuple, the `FamilyListing`
  `__iter__`/`names`/`by_name`/`to_dict` registry, and `families()`); the "not offered by
  the scalar form-builder" behaviour is carried by the prose annotation in the `summary`,
  not by any code-level exclusion list.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/skewgentle/split.py src/quiverlab/skewgentle/__init__.py \
        tests/families/test_skew_gentle_split.py
git commit -m "feat(skewgentle): SkewGentleAlgebra split constructor (Chen sec 3) -- admissible kQ_hat/I_hat, dim-certified vs associated gentle, Sp=empty byte-reduces, char-free"
```

---

### Task 3: `modules.py` — special-string re-gluing, indecomposable materialisation

The module classification: P46 strings/bands of the associated gentle algebra, typed by
special endpoints, re-glued into their `±` forms and materialised as honest A-modules on
the split algebra.

**Files:**
- Create: `src/quiverlab/skewgentle/modules.py`
- Test: `tests/modules/test_skew_gentle_modules.py`

**Interfaces:**
- Consumes: `triple.py`/`split.py`, `strings/walks.py`
  (`enumerate_strings`/`find_bands`/`letter_source`/`letter_target`/`invert`),
  `strings/modules.py::string_module`, `Module.from_arrow_action`,
  `modules/decompose.py::is_indecomposable`, `modules/hom.py::is_isomorphic`,
  `modules/ar.py::knit_ar_quiver` (P41, the rep-finite count oracle, if present).
- Produces:
  ```python
  @dataclass(frozen=True)
  class SkewGentleString:
      walk: tuple            # the associated-gentle walk (Task-1 A^g)
      type: tuple            # (r, s) in {"u","p"}^2  (Garcia-Lavoue 2.1)
      is_band: bool
      forms: tuple           # the split incarnations: 1 form for (u,u); 2 for a p-end; ...
  def classify(triple, max_length=8, budget=4096) -> list[SkewGentleString]
      # P46 census on A^g, each walk typed by whether an endpoint sits at a special
      # vertex (and is not the special loop).
      # CORRECTION (fix round, adjudicated): this is the LOOP-FREE A^g-walk census, a
      # documented STRICT SUBSET of the indecomposables at the MODULE level (headline 5 of
      # 6, mesh 8 of 11) -- NOT "complete iff rep-finite". The missing modules are the
      # loop-traversal / mixed-eigenvalue ones (e.g. the projective P_1); no loop-free
      # A^g-walk produces them, and the symmetric-string enumeration (Garcia-Lavoue Table
      # 1 / clan classification) that WOULD is a DEEPER-ENGINES-BACKLOG item.
      # skew_gentle_indecomposables (AR route) is the AUTHORITATIVE enumeration.
  def skew_gentle_module(triple, sgstring, form=0, field=None) -> Module
      # materialise the A-module on the SPLIT algebra for the chosen `form`. Uses the
      # split copies i_plus/i_minus as the two forms of a special string; the ordinary
      # string maps through the unique copies. Self-certifies via check_module. (A member
      # of the loop-free census -- the strict subset above.)
  def skew_gentle_indecomposables(triple, max_length=8, budget=4096, field=None) -> list[Module]
      # AUTHORITATIVE: ALL indecomposables via the P41 AR quiver, run over the char-FREE
      # split model (QQ) -- a presentation invariant, NOT the caller's field (M3, mirroring
      # the certificate; the char-2 regression). Falls back to the loop-free string census
      # (a sound partial sample) only on a loud non-complete AR status -- never a silent
      # read of "rep-infinite". Each is_indecomposable-checked over QQ.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_skew_gentle_modules.py
"""Skew-gentle module classification via special-string re-gluing (Plan 68; HZZ / GL
Table 1). Self-cert: every materialised module passes check_module (inside
from_arrow_action) and is indecomposable. Cross-engine: on a rep-finite skew-gentle
algebra the number of materialised indecomposables equals the AR-quiver vertex count.
Literature: a special (type-p) string has exactly TWO forms (the +/- split vertices);
an ordinary (u,u) string has one. Over QQ / GF(32003) (decompose char caveat)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.decompose import is_indecomposable
from quiverlab.skewgentle.modules import (classify, skew_gentle_indecomposables,
                                          skew_gentle_module)
from quiverlab.skewgentle.triple import SkewGentleTriple

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


def _t():                              # Q: 1 --a--> 2 , Sp = {2}
    return SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})


@selfcert
def test_materialised_modules_are_indecomposable():
    t = _t()
    for M in skew_gentle_indecomposables(t, field=QQ):
        assert is_indecomposable(M)                  # check_module already passed at build


@lit
def test_special_string_has_two_forms():
    t = _t()
    specials = [s for s in classify(t) if "p" in s.type and not s.is_band]
    assert specials, "expected at least one type-p string over Sp={2}"
    for s in specials:
        assert len(s.forms) == 2                     # the +/- forms = split vertices
    ordinary = [s for s in classify(t) if s.type == ("u", "u")]
    for s in ordinary:
        assert len(s.forms) == 1


@xeng
def test_indecomposable_count_equals_ar_quiver_count():
    # a rep-finite skew-gentle algebra: #materialised indecomposables == #AR vertices.
    from quiverlab.modules.ar import knit_ar_quiver
    from quiverlab.skewgentle.split import SkewGentleAlgebra
    t = _t()
    A = SkewGentleAlgebra(t, field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    mods = skew_gentle_indecomposables(t, field=QQ)
    assert len(mods) == len(ar.vertices)


@selfcert
def test_two_forms_are_non_isomorphic():
    from quiverlab.modules.hom import is_isomorphic
    t = _t()
    s = next(x for x in classify(t) if "p" in x.type and not x.is_band)
    M0 = skew_gentle_module(t, s, form=0, field=QQ)
    M1 = skew_gentle_module(t, s, form=1, field=QQ)
    assert not is_isomorphic(M0, M1)                  # the +/- forms are distinct A-modules
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.skewgentle.modules`

- [ ] **Step 3: Implement** `classify` (P46 census on `associated_gentle` + endpoint
  typing), `skew_gentle_module` (walk → split-copy arrow action → `from_arrow_action`),
  and `skew_gentle_indecomposables`.

**Adjust to reality (Task 3):**
- **The `(r, s)` typing (Garcia–Lavoué §2.1):** `r = p` iff `t(w₁) ∈ Sp` and `w₁` is not
  the special loop; `s = p` iff `s(wₙ) ∈ Sp` and `wₙ` is not the special loop (endpoints
  of the associated-gentle walk). A `p`-end contributes a two-way `±` choice = the two
  split vertices; count `forms = 2^{#p-ends}` for a string, and a band re-glues per the
  `k[T,T⁻¹]` / rotation rule. In the **split model** the `±` is literally which copy
  (`i⁺`/`i⁻`) the string's special end lands on — **characteristic-free**, sidestepping
  the classical char ≠ 2 `k[T]/(T²−1)` split (GL Table 1).
- **Materialisation on the split algebra:** map the associated-gentle walk to split-copy
  arrows via `split_quiver`'s `arrow_split`, then `Module.from_arrow_action` — its
  `check_module` is the self-certificate (a walk that crossed a relation is caught at
  build). For an ordinary string this reduces to `string_module` on the split algebra;
  reuse `strings.modules.string_module` where the walk lifts uniquely.
- **The count arbiter is `knit_ar_quiver` (P41) on a rep-finite instance** — if the
  materialised count over/under-shoots, a `p`-string is being given the wrong number of
  forms, or two forms collide (fix the `±` copy assignment). Do NOT weaken the count
  assertion. If P41 `knit_ar_quiver` is unavailable/incomplete on an instance, fall back
  to `decompose`-of-`⊕`-census parity or the P45 exchange-graph vertex count as the
  cross-oracle, and note it.
- **Char caveat:** `is_indecomposable`/`is_isomorphic` are rigorous over QQ / char > dim;
  batteries run over QQ (GF(32003) parity where cheap). The classification itself
  (`classify`) is pure combinatorics on `A^g` and is field-agnostic.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/skewgentle/modules.py tests/modules/test_skew_gentle_modules.py
git commit -m "feat(skewgentle): special-string re-gluing -- classify (u/p types) + skew_gentle_module on the split algebra; check_module + indecomposable, count = AR-quiver count"
```

---

### Task 4: `certificate.py` — brick-finite ⇔ representation-finite

The record's certificate. Compose Demonet–Iyama–Jasso (`brick-finite ⇔
τ-tilting-finite`, any algebra) with Garcia–Lavoué Thm 3.1 (`rep-finite ⇔ brick-finite`
for skew-gentle, char ≠ 2), cross-checked against the associated gentle algebra's
admissible-band census.

**Files:**
- Create: `src/quiverlab/skewgentle/certificate.py`
- Test: `tests/modules/test_skew_gentle_certificate.py`

**Interfaces:**
- Consumes: `split.py::SkewGentleAlgebra`, `tautilting/mutation.py::exchange_graph`
  (`.is_complete`/`.status` — τ-tilting-finite iff complete), `tautilting/torsion.py::
  bricks`, `strings/walks.py::find_bands`, `triple.py::associated_gentle`,
  `Algebra.field` (char).
- Produces:
  ```python
  def is_representation_finite(triple, budget=512, field=None) -> bool | None
      # char != 2: True iff exchange_graph(split).is_complete (tau-tilting-finite <=>
      # brick-finite <=> rep-finite, DIJ + Garcia-Lavoue Thm 3.1) -- route 1, PRIMARY and
      # rigorous. The associated-gentle band census (find_bands(A^g)) is a ONE-SIDED
      # cross-check ONLY: a band found => rep-infinite (sound), but it may miss SPECIAL
      # bands so "no-bands" does NOT certify rep-finite (W4); route 1 is authoritative.
      # char == 2: the rep-finite reading is WITHHELD (returns None) -- the rep-finite
      # equivalence is char != 2 only (Garcia-Lavoue); the tau-tilting-finiteness verdict
      # is computed on the char-free split model over QQ (M3, never decompose over char 2)
      # and lives in brick_finite_certificate["tau_tilting_finite"], scope-flagged. None
      # also iff the budget capped without a verdict (loud status carried, never a guess).
  def brick_finite_certificate(triple, budget=512, field=None) -> dict
      # {"tau_tilting_finite": bool|None, "num_bricks": int|None, "rep_finite": bool|None,
      #  "char": int, "scope": "char!=2"|"char==2 (narrowed)", "status": "complete"|"budget",
      #  "band_route": "no-bands"|"bands"|"unknown"}  -- the full two-route payload.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_skew_gentle_certificate.py
"""brick-finite <=> representation-finite for skew-gentle algebras (Plan 68).
Certificate = DIJ (brick-finite <=> tau-tilting-finite) o Garcia-Lavoue Thm 3.1
(rep-finite <=> brick-finite, char != 2). Cross-engine: the exchange-graph completeness
route (PRIMARY) and the associated-gentle band-census route (a ONE-SIDED no-false-negative
cross-check on rep-infiniteness -- find_bands(A^g) may miss special bands, W4). Literature:
a rep-finite skew-gentle algebra certifies True; a rep-infinite one (a band) certifies
False."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.skewgentle.certificate import (brick_finite_certificate,
                                              is_representation_finite)
from quiverlab.skewgentle.triple import SkewGentleTriple

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _rep_finite_triple():             # Q: 1 --a--> 2, Sp = {2}: split is hereditary A-type
    return SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})


def _rep_infinite_triple():           # a 2-cycle band + a special vertex: rep-infinite
    Q = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)})
    return SkewGentleTriple.make(Q, relations=["a*b", "b*a"], special={1})


@lit
def test_rep_finite_certifies_true():
    assert is_representation_finite(_rep_finite_triple(), field=QQ) is True


@lit
def test_rep_infinite_certifies_false():
    assert is_representation_finite(_rep_infinite_triple(), field=QQ) is False


@xeng
def test_band_route_is_a_sound_one_sided_check():
    # W4: find_bands(A^g) may be INCOMPLETE for SPECIAL bands (bands touching Sp), so the
    # band-census route is NOT a full cross-engine AGREEMENT oracle. It is a ONE-SIDED,
    # no-false-negative check on rep-INFINITENESS: a band FOUND certifies rep-infinite
    # (sound -- a real band is a real band), but "no-bands" does NOT on its own certify
    # rep-finite unless the census is extended to special bands. The exchange-graph route
    # (route 1) is the PRIMARY, rigorous (char != 2) decider; route 2 may only DOWNGRADE a
    # verdict to rep-infinite, never upgrade it to rep-finite.
    for t in (_rep_finite_triple(), _rep_infinite_triple()):
        cert = brick_finite_certificate(t, field=QQ)
        if cert["band_route"] == "bands":
            assert cert["rep_finite"] is False        # found band => genuinely rep-infinite
        # per-instance coincidence for these two hand-picked cases (documented, NOT a
        # general guarantee): the census happens to be complete here.
        assert cert["rep_finite"] == (cert["band_route"] == "no-bands")
        assert cert["tau_tilting_finite"] == cert["rep_finite"]   # char != 2 (route 1)


@lit
def test_char2_scope_is_narrowed_not_lying():
    from quiverlab.fields import GF
    cert = brick_finite_certificate(_rep_finite_triple(), field=GF(2))
    assert cert["char"] == 2 and cert["scope"] == "char==2 (narrowed)"
    # M3: the tau-tilting-finiteness verdict is computed on the char-FREE split model over
    # a decompose-rigorous field (QQ) -- NOT by running the exchange graph / decompose over
    # GF(2) (rigorous only char 0 or char > dim M). The rep-finite UPGRADE (char != 2,
    # Garcia-Lavoue) is WITHHELD -- never a guessed rep-finite over char 2.
    assert cert["tau_tilting_finite"] in (True, False)   # decided on the QQ model
    assert cert["rep_finite"] is None                    # withheld over char 2 (honest)
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement.** Build the split algebra, run `exchange_graph(A, budget)` →
  `is_complete` = τ-tilting-finite (route 1, PRIMARY); `len(bricks(A, budget))` where
  complete; `find_bands(associated_gentle(triple))` for the route-2 ONE-SIDED cross-check
  (band found ⟹ rep-infinite; incomplete for special bands, so never used to certify
  rep-finite — W4); over char 2 compute route 1 on the QQ split model and withhold the
  rep-finite upgrade (M3); compose per char.

**Adjust to reality (Task 4):**
- **The two routes (route 2 is a ONE-SIDED cross-check, not a full agreement oracle —
  W4).** (route 1, PRIMARY) `exchange_graph(split).is_complete` = τ-tilting-finite =
  brick-finite (DIJ) → rep-finite (Garcia–Lavoué, char ≠ 2) — rigorous over char ≠ 2.
  (route 2, CROSS-CHECK) a skew-gentle algebra is rep-finite iff its **admissible band**
  set is empty (HZZ / GL band classification), probed via
  `find_bands(associated_gentle(triple))`. **`find_bands(A^g)` may be INCOMPLETE for
  *special* bands** (bands touching `Sp` — in the skew-gentle algebra these re-glue
  through the `±` split copies, and a band of the non-monomial split has no direct
  `find_bands` enumerator). So route 2 is **only sound in one direction**: a band FOUND
  certifies rep-INFINITE (a real band is a real band), so it may only DOWNGRADE route 1's
  verdict — but "no-bands" does NOT by itself certify rep-finite. Therefore the
  `oracle_crossengine` test asserts the **sound implication** (`band found ⟹ rep-infinite`,
  never a false negative on rep-infiniteness), NOT full agreement; the two hand-picked
  cases happen to coincide (documented per-instance, not a general guarantee). **To
  upgrade route 2 to a complete agreement oracle** the census must be extended to special
  bands — enumerate the bands of `A^g` and, for each that touches `Sp`, add its `±`
  re-glued split forms (the special-band analogue of the string re-gluing in Task 3) — a
  backlog item, not required for v1; until then route 1 is authoritative.
- **char 2 honesty — reconciled with the `decompose` rigor window (M3).** DIJ
  (`brick-finite ⇔ τ-tilting-finite`) holds in all characteristic **as a theorem**, but
  the ENGINE that computes τ-tilting-finiteness — `exchange_graph` — leans on
  `is_isomorphic` / `is_indecomposable` / `decompose` (see `tautilting/mutation.py:12`
  "Runs over QQ by default (the decompose / is_isomorphic char caveat ...)" and the
  `mutation.py:61` hint "over char <= dim the decompose/is_isomorphic caveat can refuse;
  run over QQ", plus `tautilting/pairs.py`/`torsion.py` `is_isomorphic`/`is_indecomposable`
  usage). Those primitives are rigorous only over **char 0 or char > dim M**, so a
  τ-tilting-finiteness verdict computed **directly over GF(2)** is NOT certified once
  modules of `dim ≥ 2` appear (char 2 ≤ dim M). So the certificate must NOT claim a
  rigorous char-2 τ-tilting decision from a GF(2) BFS. Honest handling: when asked over
  char 2, compute the τ-tilting-finiteness verdict on the **char-free split model over a
  decompose-rigorous field (QQ)** — legitimate because the split *construction* is
  characteristic-free (Chen) — report that verdict, and withhold ONLY the rep-finite
  upgrade (char ≠ 2, Garcia–Lavoué) with `scope="char==2 (narrowed)"`. Never claim
  rep-finiteness over char 2, and never rely on `decompose` over char 2 (the project-wide
  "never decompose over char 2" rule). State this on the verification page.
- **Budget honesty:** if `exchange_graph` returns `status="budget"` (τ-tilting-infinite),
  the algebra is rep-infinite (over char ≠ 2) — `False`, not `None`. `None` is reserved
  for the genuinely-undecided case where a route refuses (e.g. the band route unavailable
  AND the exchange graph budget-capped without closing). Carry `status` in the payload.
- **char caveat:** the exchange-graph BFS + bricks run over QQ by default (P45 caveat);
  the batteries follow.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/skewgentle/certificate.py tests/modules/test_skew_gentle_certificate.py
git commit -m "feat(skewgentle): brick-finite <=> rep-finite certificate (DIJ o Garcia-Lavoue Thm 3.1) -- two-route (exchange-graph / band-census) cross-check, honest char-2 narrowing"
```

---

### Task 5: support τ-tilting via the engine + geometric-vs-engine cross-check

The split algebra is a presented f.d. algebra, so P45 runs on it verbatim. This task
wires the **engine-primary** τ-tilting path and pins the **geometric model** as the
cross-check oracle (rank + small counts).

**Files:**
- Modify: `src/quiverlab/skewgentle/__init__.py` (re-export helpers),
  `src/quiverlab/skewgentle/certificate.py` (a `support_tau_tilting(triple, budget)`
  thin wrapper over `tau_tilting_block(split)`)
- Test: `tests/modules/test_skew_gentle_tau_tilting.py`

**Interfaces:**
- Consumes: `split.py::SkewGentleAlgebra`, `tautilting/block.py::tau_tilting_block`,
  `tautilting/mutation.py::exchange_graph`.
- Produces:
  ```python
  def support_tau_tilting(triple, budget=512, field=None) -> dict
      # tau_tilting_block(SkewGentleAlgebra(triple)) -- the exchange graph, g-matrices,
      # brick-labelled Hasse edges, counts; honest complete-iff contract. Records
      # n = |Q_0| + |Sp| (the geometric rank; HZZ sec 6 |R| = |Q_0| + |Sp|).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_skew_gentle_tau_tilting.py
"""Support tau-tilting for skew-gentle algebras via the engine, geometric cross-check
(Plan 68 / HZZ sec 5-6). Engine-primary: P45 exchange_graph on the split algebra.
Literature/geometric: the pair rank equals |Q_0| + |Sp| (HZZ sec 6); k x k (one vertex,
Sp = {1}) has 4 support tau-tilting pairs (the Boolean square, geometrically the 4
generalized dissections)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.skewgentle.certificate import support_tau_tilting
from quiverlab.skewgentle.split import SkewGentleAlgebra
from quiverlab.skewgentle.triple import SkewGentleTriple
from quiverlab.tautilting.mutation import exchange_graph

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def test_pair_rank_is_Q0_plus_Sp():
    # HZZ sec 6: a basic support tau-tilting module has |Q_0| + |Sp| summands.
    t = SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})
    blk = support_tau_tilting(t, field=QQ)
    assert blk["n"] == 2 + 1                           # |Q_0| + |Sp|


@lit
def test_kxk_has_four_support_tau_tilting_pairs():
    # Q = one vertex, Sp = {1}: A = k[eps]/(eps^2-eps) = k x k; split = two isolated
    # vertices; stau-tilt(k x k) = 2 x 2 = 4 (the geometric count = 4 dissections).
    # NOTE (M2): k x k is semisimple -> selfinjective, WITH Sp != empty. This does NOT
    # contradict Chen Cor 1.2(c): that corollary is about INDECOMPOSABLE algebras, and
    # k x k is DECOMPOSABLE (a product of two copies of k). VERIFIED live: dim 2,
    # is_selfinjective() True, exchange graph complete with 4 pairs. This test pins the
    # tau-tilting count only; it deliberately does NOT assert the unqualified
    # selfinjective-iff (there is no such assertion anywhere in the battery).
    t = SkewGentleTriple.make(Quiver([1], {}), relations=[], special={1})
    A = SkewGentleAlgebra(t, field=QQ)
    eg = exchange_graph(A, budget_pairs=64)
    assert eg.is_complete and len(eg.vertices) == 4


def test_exchange_graph_on_split_does_not_raise():
    # H1 REGRESSION: split vertex labels must be a single orderable (string) type, else
    # P45's exchange_graph crashes at `sorted(pair.support)` (mutation.py::_label) with
    # `TypeError: '<' not supported between instances of 'str' and 'int'`. The headline
    # example Q = [1, 2], Sp = {2} splits vertex 2 into "2+"/"2-" while keeping 1 -- if 1
    # stays an int this raises. With ALL labels stringified it runs clean. (VERIFIED live:
    # mixed [1,'2+','2-',3] crashes here; all-string runs and returns a complete graph.)
    t = SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})
    A = SkewGentleAlgebra(t, field=QQ)
    eg = exchange_graph(A)                              # must NOT raise TypeError
    assert eg.is_complete


@xeng
def test_geometric_rank_matches_engine_pairs():
    t = SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})
    A = SkewGentleAlgebra(t, field=QQ)
    eg = exchange_graph(A)
    n = len(list(A.quiver.vertices))                  # = |Q_0| + |Sp|
    for rec in eg.vertices:                           # every pair has full rank n
        assert len(rec["summand_dimvecs"]) + len(rec["support"]) == n


@lit
@pytest.mark.deep
def test_hzz_section6_seven_vertex_rank_is_eight():
    # HZZ sec 6 worked example: vertices 1..7, arrows a..h, loop eps_1 (Sp = {1}), a
    # NON-special loop eps_5 (eps_5^2 in I so vertex 5 is ordinary), I = {eps_5^2, ab, ba,
    # ed, hc}. A basic support tau-tilting module has |R| = 8 = |Q_0| + |Sp| = 7 + 1
    # summands -- the rank oracle n_split = |Q_0| + |Sp| (HZZ sec 6). Pinned CHEAPLY via
    # the split vertex count (a complete support tau-tilting pair over a rank-n algebra has
    # exactly n summands, by definition), NOT by enumerating the 8-vertex exchange graph.
    Q = Quiver([1, 2, 3, 4, 5, 6, 7],
               {"a": (1, 2), "b": (2, 1), "c": (2, 3), "d": (3, 4), "e": (4, 3),
                "f": (4, 5), "g": (5, 6), "h": (6, 2), "eps5": (5, 5)})   # eps_5 ordinary
    t = SkewGentleTriple.make(Q, relations=["eps5*eps5", "a*b", "b*a", "e*d", "h*c"],
                              special={1})
    A = SkewGentleAlgebra(t, field=QQ)
    assert len(list(A.quiver.vertices)) == 7 + 1       # |R| = |Q_0| + |Sp| = 8 (HZZ sec 6)
    # The heavier confirmation -- one support tau-tilting pair actually has 8 summands via
    # the engine -- is the deep exchange-graph run; keep it OPTIONAL/behind budget if slow.
```

> **HZZ §6 caveat (to verify at implementation):** the arrow set / orientation above is a
> reconstruction consistent with the §6 prose (`Sp = {1}`, ordinary loop `ε₅`, relations
> `{ε₅², ab, ba, ed, hc}`); confirm it forms a valid skew-gentle triple against the paper
> when wiring the test (adjust the arrows so `is_skew_gentle_triple` holds), then the
> `|R| = |Q₀| + |Sp| = 8` rank pin is the literature oracle. The rank is
> presentation-independent (`|Q₀| + |Sp|`), so only the triple's *validity* needs the
> paper, not the arrow names.

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** `support_tau_tilting` (a thin `tau_tilting_block(split)`
  wrapper stamping `n = |Q₀| + |Sp|`).

**Adjust to reality (Task 5):**
- **Engine-primary, geometric-as-oracle** is the scope decision (design question #4):
  do NOT build the full orbifold/dissection enumerator in v1. The engine (P45) computes
  `sτ-tilt` on the split algebra; the geometric model contributes only *cross-check
  values* — the rank `n = |Q₀| + |Sp|` (always) and small `|sτ-tilt|` counts (pin the
  `k×k → 4` case; add one rep-finite `A_n`-type skew-gentle count if a transcribable
  geometric value exists, `# PIN` otherwise).
- **`# PIN` the τ-tilting counts, attribution honest:** the record cites "the 2212.06467
  [Chen] example counts", but Chen is homological — the numeric τ-tilting counts are
  geometric (HZZ). Wire Chen's homological claims as the "2212.06467 oracles" (Task 2's
  selfinjective/dim tests) and pin the τ-tilting counts from HZZ §6 (rank) + the `k×k`
  hand case; state the attribution split on the verification page.
- **τ-tilting-infinite honesty:** on a rep-infinite skew-gentle algebra the exchange
  graph is `status="budget"` (loud) — the block carries `complete=False`; never a silent
  partial graph (P45 contract).

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/skewgentle/certificate.py src/quiverlab/skewgentle/__init__.py \
        tests/modules/test_skew_gentle_tau_tilting.py
git commit -m "feat(skewgentle): support_tau_tilting via P45 on the split algebra -- geometric rank |Q_0|+|Sp| + k x k = 4 pins as the geometric-vs-engine cross-check"
```

---

### Task 6: `skew_gentle_block` + the `skew_gentle` GUI scalar kind + preset

The no-code surface: the triple input (quiver + relations + special-loop picks) → a
compute block (recognizer verdict, split shape, dim law, classification counts,
τ-tilting via the engine, rep-type certificate), clickable end-to-end.

**Files:**
- Create: `src/quiverlab/skewgentle/block.py`
- Modify: `src/quiverlab/hpc/spec.py` (a `skew_gentle` scalar-kind dispatch +
  `_snip`; the `family`/non-scalar-constructor route so the GUI can pass a triple),
  `docs/gui/runner.py` (the Pyodide twin: matching dispatch + `_snip` + ETA),
  `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox, `S.ids`, push-list,
  `renderBlock`, `scheduleProbe`, a **skew-gentle preset**: the special-loop picker on
  the canvas), `webapp/templates/index.html` (checkbox + preset),
  `webapp/server/i18n/en.json` + `es.json` (`inv.skew_gentle`, `block.skew_gentle.*`),
  `src/quiverlab/trace/results_html.py` (`_HEADINGS` + a `_skew_gentle_html` branch),
  `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE golden, existing byte-identical)
- Test: `tests/webapp/test_skew_gentle_block_p68.py`,
  `tests/gui/test_skew_gentle_runner_twin.py`

**Interfaces (mirror the `strings` / `tau_tilting` algebra-level kinds):**
- `skew_gentle_block(triple_or_algebra, budget=512)` returns a dict; each runner stamps
  `block["citations"] = _citation_pairs(block["references"])`. Shape:
  ```python
  {"kind": "skew_gentle",
   "is_skew_gentle": True,
   "triple": {"vertices": [...], "arrows": {...}, "relations": [...], "special": [...]},
   "split": {"num_vertices": n_split, "num_arrows": int, "dim": int},
   "dim_law": {"split_dim": int, "assoc_gentle_dim": int, "ok": bool},   # HZZ Lemma 1.5
   "rank": n_split,                                       # |Q_0| + |Sp|
   "classification": {"num_indecomposables": int|None, "num_special": int,
                      "has_bands": bool, "status": "complete"|"budget"},
   "tau_tilting": {"num_pairs": int, "complete": bool, "status": str},   # engine
   "rep_type": {"rep_finite": bool|None, "scope": str, "status": str},   # certificate
   "references": ["he_zhou_zhu", "chen_skew_gentle", "amiot_skew_gentle",
                  "garcia_lavoue", "assem_book"]}
  ```
  A non-triple input (a bare presented algebra without the construction marker) reports
  `is_skew_gentle` honestly (True only if `_skew_gentle_triple` is present) + a `note`.
- `skew_gentle` is an **algebra-level scalar kind** built from a **triple** request
  (the GUI passes the special-loop picks alongside the quiver+relations); route it like
  the non-scalar `family` constructors (`BrauerGraphAlgebra` precedent) that carry
  explicit quiver data, then dispatch `skew_gentle_block`. **Tier sizing (W2 — correct
  attribution):** `sizing_dim` lives in **`webapp/server/estimator.py`** (NOT
  `quiverlab.hpc`), with signature `sizing_dim(algebra_dim: int, req) -> int` — it takes a
  **pre-computed** `algebra_dim` and returns `max(algebra_dim, module dims, …)`. For a
  skew-gentle triple request the `algebra_dim` that flows in must be the **split
  algebra's** dim: `webapp/server/app.py` computes it in the classify preamble via
  `A = _build_or_error(req.algebra); dim = A.dim; classify(sizing_dim(dim, req), req, cfg)`
  (app.py ~line 187–195, run BEFORE dispatch and outside the wall net). So `build_algebra`
  (`hpc/spec.py` / `webapp/server/runner.py`) must route a skew-gentle triple request to
  `SkewGentleAlgebra(triple)`, whose `.dim` (the split dim, already dimension-certified) is
  the value passed to `sizing_dim` — the split is constructed once, up front, exactly as
  any other `family`/quiver request's primary algebra is, and no new estimator entry point
  is needed.

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir; copy
  the `strings`/`tau_tilting`-kind runner-pair fixture):

```python
# tests/webapp/test_skew_gentle_block_p68.py
"""The `skew_gentle` algebra-level scalar kind: served by hpc.spec from a triple,
mirrored by the Pyodide twin, byte-identical block."""


def test_skew_gentle_block_shape(tmp_path):
    # triple: Q = 1 --a--> 2, Sp = {2}, over QQ, compute ["skew_gentle"]. Assert:
    #   block["is_skew_gentle"] is True
    #   block["dim_law"]["ok"] is True and split_dim == assoc_gentle_dim == 5
    #   block["rank"] == 3                                 # |Q_0| + |Sp|
    #   block["rep_type"]["rep_finite"] is True
    #   "he_zhou_zhu" in [k for k, _ in block["citations"]]
    ...


def test_sp_empty_block_matches_gentle(tmp_path):
    # Sp = empty triple -> the block reports the plain gentle path (byte-reduction).
    ...


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py; json.dumps(sort_keys=True)
    # equality on the skew_gentle block (both runners byte-identical).
    ...
```

- [ ] **Step 2: Implement** `skewgentle/block.py` + the `spec.py` branch + the
  `docs/gui/runner.py` twin (shape-identical), the two `gui.js` touchpoints + the
  `index.html` checkbox **and the special-loop preset** (a canvas control to mark
  vertices as special), the ETA entry, the i18n keys
  (`inv.skew_gentle`, `block.skew_gentle.title`, `block.skew_gentle.dimlaw`,
  `block.skew_gentle.classification`, `block.skew_gentle.tau`, `block.skew_gentle.reptype`
  — EN and ES), the `_snip` recipe (`"skew_gentle": "SkewGentleAlgebra(triple)"` /
  `skew_gentle_block(triple)`), and the `results_html.py`
  `_HEADINGS["skew_gentle"] = "Skew-gentle"` + `_block_html` branch.

- [ ] **Step 3: Add ONE golden fixture** (`skew_gentle_arrow_sp2`) to
  `_runner_goldens.json`; note it in `test_runner_delegation.py`'s change-log. Run the
  delegation test BEFORE adding to confirm existing goldens stay byte-identical.
  **Canonical keys:** the triple request canonicalises through the Plan-25 `canonical_key`
  (special-set sorted, default explicit) — confirm a Sp-ordering permutation collides
  and a Sp change does not (mirror the module-block canonicalisation guard).

- [ ] **Step 4: Run the gates**

Run: `... -m pytest tests/webapp/test_skew_gentle_block_p68.py tests/webapp/test_runner_delegation.py tests/gui/test_skew_gentle_runner_twin.py tests/hpc -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): skew_gentle scalar kind + special-loop preset -- triple input, dim law + classification + tau-tilting + rep-type, both runners byte-identical, EN/ES, one golden"
```

---

### Task 7: QPA cross-oracle battery

**Files:**
- Modify: `src/quiverlab/qpa/scripts.py` / `crosscheck.py` (only if a new helper is
  needed; the recognizer + decompose + selfinjective crosschecks already exist)
- Test: `tests/qpa/test_skew_gentle_qpa.py`

**Interfaces:**
- Consumes: the QPA session (`session.should_skip_qpa()`, `session.run`,
  `session.libgap_handle()`), the existing `A.crosscheck(...)` bridge (P30/P38: the
  split algebra is a plain presented `kQ/I`, so `decompose` / `is_selfinjective` /
  `cartan_matrix` crosscheck through the existing machinery),
  `skewgentle/split.py::SkewGentleAlgebra`.
- Produces: `tests/qpa/test_skew_gentle_qpa.py`. **Honest scope: QPA has NO skew-gentle
  surface** (no `SkewGentleAlgebra`, no split constructor, no skew-gentle recognizer) —
  the crosschecks are therefore on the **split algebra as a plain `kQ/I`**: dimension /
  Cartan / `IsSelfinjectiveAlgebra` / `DecomposeModuleWithMultiplicities`, with a
  standing `NamesGVars()` guard that FAILS if QPA ever ships a skew-gentle surface (the
  Plan-35 skip-that-fails-if-appears precedent).

- [ ] **Step 1: Probe live QPA + write the battery**

```python
# tests/qpa/test_skew_gentle_qpa.py
"""QPA as the oracle for skew-gentle algebras (Plan 68). QPA has NO skew-gentle
constructor / recognizer -- so the crosschecks are on the SPLIT algebra as a plain
kQ/I: dim, Cartan, selfinjectivity (Chen Cor 1.2c), and module decompose. qpa-marked."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.qpa import session
from quiverlab.skewgentle.split import SkewGentleAlgebra
from quiverlab.skewgentle.triple import SkewGentleTriple

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _t():
    return SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})


def test_qpa_has_no_skew_gentle_surface():
    lg = session.libgap_handle()
    for name in ("SkewGentleAlgebra", "IsSkewGentleAlgebra", "SkewGentleTriple"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), \
            f"QPA now ships {name} -- add a real crosscheck (honest-scope changed)"


def test_split_algebra_dim_and_cartan_vs_qpa():
    A = SkewGentleAlgebra(_t(), field=QQ)
    A.crosscheck("cartan_matrix").assert_agree()         # split is a plain kQ/I
    A.crosscheck("dim").assert_agree()


def test_selfinjective_crosscheck():
    # Chen Cor 1.2(c): the 2-cycle gentle (Sp = empty) is selfinjective; QPA agrees.
    Q2 = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)})
    A = SkewGentleAlgebra(quiver=Q2, relations=["a*b", "b*a"], special=(), field=QQ)
    A.crosscheck("is_selfinjective").assert_agree()      # IsSelfinjectiveAlgebra
```

- [ ] **Step 2: Wire any missing crosscheck case.** `cartan_matrix` / `dim` /
  `is_selfinjective` / `decompose` crosschecks already exist (P30/P38/P29). If
  `A.crosscheck("dim")` is not a dispatch case, add it mirroring the existing scalar
  crosschecks; otherwise no `src/` change. The skew-gentle-surface probe is the standing
  guard.

- [ ] **Step 3: Run live** `... -m pytest tests/qpa/test_skew_gentle_qpa.py -v` (the venv
  has the `[qpa]` extra). Expected: PASS live; the surface probe confirms absence.

- [ ] **Step 4: Commit**

```bash
git add src/quiverlab/qpa tests/qpa/test_skew_gentle_qpa.py
git commit -m "test(qpa): skew-gentle split algebra crosschecks (dim/Cartan/selfinjective/decompose as a plain kQ/I) + honest no-skew-gentle-surface guard"
```

---

### Task 8: verification page, citations, README, family catalog, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `src/quiverlab/families/__init__.py` (export `SkewGentleAlgebra`),
  `src/quiverlab/families/discover.py` (append ONE `FamilyInfo(...)` to the `CATALOG`
  tuple, `route="general"`, a prose `summary` annotating the non-scalar constructor
  exactly like `BrauerGraphAlgebra`'s "Non-scalar constructor: not offered by the scalar
  form-builder." — there is no `_iter_families`/skip-set to touch),
  `src/quiverlab/skewgentle/__init__.py` (public surface),
  `src/quiverlab/__init__.py` (re-export `SkewGentleAlgebra` /
  `is_skew_gentle_triple` at top level, mirroring `BrauerGraphAlgebra`)
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick the P68 card)
- Test: existing release gates (`tests/release/test_oracle_classes.py`,
  `tests/citations/`)

- [ ] **Step 1: Citations** (BibTeX-VERIFIED at authoring; `_r(key, bibtex_key, kind,
  title, annotation, *tags)` registry precedent). Add:

```bibtex
@article{HeZhouZhu2020,
  author  = {He, Ping and Zhou, Yu and Zhu, Bin},
  title   = {A geometric model for the module category of a skew-gentle algebra},
  journal = {arXiv preprint},
  eprint  = {2004.11136}, archivePrefix = {arXiv}, primaryClass = {math.RT},
  year    = {2020},
}
@article{Chen2022skewgentle,
  author  = {Chen, Yiping},
  title   = {A characteristic free approach to skew-gentle algebras},
  journal = {arXiv preprint},
  eprint  = {2212.06467}, archivePrefix = {arXiv}, primaryClass = {math.RT},
  year    = {2024},
}
@article{Amiot2021skewgentle,
  author  = {Amiot, Claire},
  title   = {Indecomposable objects in the derived category of a skew-gentle algebra
             using orbifolds},
  journal = {arXiv preprint (ICRA 2020 proceedings)},
  eprint  = {2107.02646}, archivePrefix = {arXiv}, primaryClass = {math.RT},
  year    = {2021},
}
@article{GarciaLavoue2026,
  author  = {Garcia, Monica and Lavou{\'e}, L{\'e}a},
  title   = {Brick-finite skew-gentle algebras are representation-finite},
  journal = {arXiv preprint},
  eprint  = {2601.01744}, archivePrefix = {arXiv}, primaryClass = {math.RT},
  year    = {2026},
}
@article{GeissDeLaPena1999,
  author  = {Gei{\ss}, Christof and de la Pe{\~n}a, Jos{\'e} Antonio},
  title   = {Auslander-{R}eiten components for clans},
  journal = {Bolet\'in de la Sociedad Matem\'atica Mexicana. Tercera Serie},
  volume  = {5}, number = {2}, pages = {307--326}, year = {1999},
}
@article{CrawleyBoevey1989,
  author  = {Crawley-Boevey, William},
  title   = {Functorial filtrations {II}: clans and the {G}elfand problem},
  journal = {Journal of the London Mathematical Society},
  volume  = {40}, number = {1}, pages = {9--30}, year = {1989},
}
```

and in `registry.py` (mirror the `_r("assem_book", "ASS2006", "foundation", ...)` shape):

```python
_r("he_zhou_zhu", "HeZhouZhu2020", "family",
   "A geometric model for the module category of a skew-gentle algebra",
   "He-Zhou-Zhu: the skew-gentle triple (Q,Sp,I), the idempotent specialization "
   "eps^2=eps, the split/doubled quiver, and support tau-tilting via the orbifold "
   "model -- the primary skew-gentle source.", "families", "modules"),
_r("chen_skew_gentle", "Chen2022skewgentle", "family",
   "A characteristic free approach to skew-gentle algebras",
   "Chen: the char-free idempotent split construction (sec 3), the selfinjective "
   "classification (Sp=empty & gentle selfinjective), the K-theory/dim relation, and "
   "Gorensteinness -- the split constructor's binding source, valid in char 2.",
   "families"),
_r("amiot_skew_gentle", "Amiot2021skewgentle", "algorithm",
   "Indecomposable objects in the derived category of a skew-gentle algebra via orbifolds",
   "Amiot: skew-gentle as Z2-skew-group of a gentle algebra; the orbifold/double-cover "
   "geometric model backing the geometric-vs-engine tau-tilting cross-check.", "modules"),
_r("garcia_lavoue", "GarciaLavoue2026", "algorithm",
   "Brick-finite skew-gentle algebras are representation-finite",
   "Garcia-Lavoue Thm 3.1 (char != 2): brick-finite <=> rep-finite for skew-gentle -- "
   "composed with DIJ (brick-finite <=> tau-tilting-finite) gives the rep-type "
   "certificate.", "modules"),
_r("geiss_delapena", "GeissDeLaPena1999", "foundation",
   "Auslander-Reiten components for clans",
   "Geiss-de la Pena: the original skew-gentle / clan definition (char != 2).",
   "families"),
_r("crawley_boevey_clans", "CrawleyBoevey1989", "foundation",
   "Functorial filtrations II: clans and the Gelfand problem",
   "Crawley-Boevey: the classification of indecomposables for clans, underlying the "
   "special-string re-gluing (the +/- forms).", "modules"),
```

  Reuse the existing `demonet_iyama_jasso` (P45, brick-finite ⇔ τ-tilting-finite),
  `air_tau_tilting` (P45), `butler_ringel`/`avella_geiss` (P46), `assem_book`.
  **Spec-ambiguity resolution (recorded):** the record cites "the 2212.06467 [Chen]
  example counts"; Chen is homological, so the τ-tilting counts are pinned from the
  **geometric model (HZZ)** and the Chen citation backs the **homological** oracles
  (selfinjective-iff, dim, K₀ rank, gl.dim). State this attribution split on the
  verification page; do NOT invent an enumerative Chen example that the paper does not
  contain.

- [ ] **Step 2: Verification page.** Add the Plan-68 subsystem rows:
  - `skewgentle/triple.py` + `split.py` — `oracle_selfcert` (dim law
    `dim split = dim associated gentle`, HZZ Lemma 1.5, a **NECESSARY** check;
    `n_split = |Q₀|+|Sp|`; char-free split over QQ/GF(p)/GF(2); the H1 all-string
    vertex-label regression — `exchange_graph(split)` does not raise); `oracle_literature`
    (`Sp = ∅` byte-reduces to the gentle algebra; the mesh example `1→a→2→b→3, Sp={2},
    a*b ∈ I` — `dim A^g = 9` vs the zero-split's 8, gentleness forces `a*b ∈ I`; Chen Cor
    1.2(c) **indecomposable** selfinjective-iff).
  - `skewgentle/modules.py` — `oracle_selfcert` (`check_module` + indecomposable
    materialisation; a type-`p` string has two forms, distinct); `oracle_crossengine`
    (indecomposable count ≡ AR-quiver vertex count — the **SUFFICIENCY** oracle for the
    split relations; `11` on the mesh example).
  - `skewgentle/certificate.py` — `oracle_literature` (rep-finite/rep-infinite verdicts,
    Garcia–Lavoué Thm 3.1); `oracle_crossengine` (band-census route is a **ONE-SIDED**
    no-false-negative cross-check on rep-infiniteness — `find_bands(A^g)` may miss special
    bands, so the exchange-graph route is authoritative; NOT a full agreement oracle).
  - support τ-tilting — `oracle_literature` (rank `|Q₀|+|Sp|`, HZZ §6; `k×k → 4`; the HZZ
    §6 seven-vertex `|R| = 8`; the mesh example `46` pairs); `oracle_crossengine`
    (geometric rank ≡ engine pair rank).
  - `qpa` — split-algebra dim/Cartan/selfinjective/decompose parity as a plain `kQ/I`.
  Add the **honest-scope entries**: (a) the triple is NON-admissible — the internal
  object is always the split algebra; the dim law is NECESSARY, the AR/τ-tilting counts on
  the mesh example are the SUFFICIENCY oracle; (b) the brick-finite ⇔ rep-finite
  certificate is **char ≠ 2** (Garcia–Lavoué); over char 2 the τ-tilting-finiteness
  verdict is computed on the char-free split model over QQ (the exchange graph leans on
  `is_isomorphic`/`decompose`, rigorous only char 0 / char > dim, so it is NEVER run over
  GF(2) — M3) and the rep-finite upgrade is WITHHELD (`rep_finite is None`), while the
  split construction / dim law / module classification remain characteristic-free (Chen);
  (c) recognizing an arbitrary presented algebra as skew-gentle up to isomorphism is NOT
  attempted — the recognizer decides the TRIPLE (and construction-tagged algebras); (d)
  the τ-tilting numeric counts are geometric (HZZ), the 2212.06467 (Chen) oracles are
  homological (selfinjective/dim/K₀/gl.dim), attribution split; (e) the band-census
  route (route 2) is a one-sided check, not a full cross-engine agreement (W4 —
  `find_bands(A^g)` may be incomplete for special bands); (f) the `decompose`/
  exchange-graph char ≤ dim caveat (batteries over QQ / GF(32003)). **Recount the class table**
  (`tests/release/test_oracle_classes.py` drives the numbers — run collection, paste the
  LIVE counts, re-run to green; never a guessed-at-authoring number).

- [ ] **Step 3: README.** One features line: "skew-gentle algebras: the triple
  `(Q, I, Sp)` recognizer, the characteristic-free idempotent-split constructor,
  special-string module re-gluing, support τ-tilting via the engine (orbifold model as
  oracle), and the brick-finite ⇔ representation-finite certificate — R32."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/families tests/modules -q` (deep, touched dirs),
  `... -m pytest tests/invariants tests/webapp tests/gui tests/hpc -q -m fast`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q`,
  and a citation-presence check (`he_zhou_zhu`/`chen_skew_gentle`/`amiot_skew_gentle`/
  `garcia_lavoue` resolve; the `skew_gentle` block carries `he_zhou_zhu`) — all green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-68 skew-gentle oracle rows + honest scope (non-admissible triple, char!=2 certificate, geometric-vs-homological attribution) + citations + recounted classes"
```

---

## Acceptance (Plan-68 definition of done)

1. `SkewGentleTriple`/`is_skew_gentle_triple`/`associated_gentle` (triple),
   `split_quiver`/`SkewGentleAlgebra` (split), `classify`/`skew_gentle_module`/
   `skew_gentle_indecomposables` (modules), `is_representation_finite`/
   `brick_finite_certificate`/`support_tau_tilting` (certificate),
   `skew_gentle_block` (block) all public in `quiverlab.skewgentle`, with
   `SkewGentleAlgebra` re-exported through `quiverlab.families` and `quiverlab`; every
   constructor loudly-validated and dimension-certified.
2. **The triple is the user input; the split algebra is the internal object.** A direct
   `⟨εᵢ²−εᵢ⟩` presentation is impossible (`AdmissibilityError`, VERIFIED); every
   invariant is computed on the admissible split algebra `kQ̂/Î`, dimension-certified
   `dim(split) == dim(associated_gentle)` (HZZ Lemma 1.5) — a **NECESSARY** check; the
   mesh construction's **sufficiency** oracle is the AR / support-τ-tilting counts on a
   mesh-bearing example (`1→a→2→b→3, Sp={2}, a*b ∈ I`: `dim A^g = 9` vs the zero-split's
   8, `11` AR vertices, `46` τ-tilting pairs — all VERIFIED live). `n_split = |Q₀|+|Sp|`.
3. **`Sp = ∅` byte-reduces to the gentle engine:** `SkewGentleAlgebra(Sp=∅)` returns the
   plain gentle `A(Q,I)` presentation verbatim (byte-identical dim + Cartan + arrows +
   block) — the strongest self-cert **for the no-special case that bypasses the split
   machinery** (it does not exercise the mesh/split construction; the counts in #2 do).
4. **Special-string re-gluing:** the associated gentle algebra's strings/bands are typed
   `(r,s) ∈ {u,p}²`; a type-`p` string has exactly two forms (the split vertices
   `i⁺, i⁻`, distinct A-modules), an ordinary string one; every materialised module
   passes `check_module` and is `is_indecomposable`; the count equals the AR-quiver
   vertex count on rep-finite instances (over QQ / GF(32003)), pinned on the mesh example
   (`11` AR vertices).
5. **Support τ-tilting is engine-primary, geometric-as-oracle:** `tau_tilting_block` runs
   on the split algebra; the pair rank equals `|Q₀|+|Sp|` (HZZ §6), `k×k` has 4 support
   τ-tilting pairs, and the geometric rank matches the engine pair rank — the
   geometric-vs-engine cross-check.
6. **brick-finite ⇔ rep-finite certificate:** DIJ (`brick-finite ⇔ τ-tilting-finite`,
   general) composed with Garcia–Lavoué Thm 3.1 (`rep-finite ⇔ brick-finite`,
   char ≠ 2); the exchange-graph route is PRIMARY and the band-census route is a ONE-SIDED
   cross-check (a band found ⟹ rep-infinite; `find_bands(A^g)` may miss special bands so
   "no-bands" alone does not certify rep-finite — W4); a
   rep-finite triple certifies `True`, a rep-infinite one `False`; over char 2 the
   verdict is honestly narrowed to τ-tilting-finiteness (never a false rep-finite claim).
7. The `skew_gentle` algebra-level scalar kind is clickable end-to-end (GUI canvas +
   special-loop preset → block → report) in EN+ES, both runners byte-identical, ONE
   golden added with a documented change-log entry, canonicalising through the Plan-25
   key (Sp-order-invariant); the dim law, classification counts, τ-tilting, and rep-type
   render.
8. Live QPA battery green (`-m qpa`): split-algebra dim/Cartan/`IsSelfinjectiveAlgebra`/
   `DecomposeModuleWithMultiplicities` parity as a plain `kQ/I`, with the standing
   `NamesGVars()` guard that FAILS if QPA ever ships a skew-gentle surface (honest scope:
   QPA has none today).
9. `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the five
   honest-scope entries (non-admissible triple → split; char ≠ 2 certificate + char-free
   construction; recognizer decides the triple, not iso-type; τ-tilting counts geometric
   vs Chen homological attribution; decompose char caveat);
   `he_zhou_zhu`/`chen_skew_gentle`/`amiot_skew_gentle`/`garcia_lavoue`/`geiss_delapena`/
   `crawley_boevey_clans` citations added and BibTeX-verified; `demonet_iyama_jasso`/
   `air_tau_tilting`/`butler_ringel`/`avella_geiss` reused; README line added; deep
   (touched dirs) + fast + qpa + release + citations suites green. No dependency taken on
   any other v1.0.0 subplan (merges independently to `dev`).
