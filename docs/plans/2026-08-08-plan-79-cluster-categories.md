# Plan 79: Amiot–Keller cluster categories — the certified module-category slice (P79 / R31)

**Date:** 2026-08-08 · **Wave:** 5 (research-grade tail; the LAST record-plan before
the P80 release gate) · **Depends on:** P67 (silting/derived stack) — used only for the
*out-of-scope* boundary; P45 (τ-tilting exchange graph) + P63 (n-regularity recovery) —
the **computational backbone**; P44 (`families/jacobian.py`) — the cluster-tilted
End-algebra; P23/P24 (τ, AR translate) + P05 (module Ext) + the AR knitting engine
(`modules/ar.py`) — the module-category model. · **Trunk:** `dev`. · **Branch:**
`plan-79-cluster-categories`.

**One-line:** Ship the certified, finite, module-category slice of the Amiot–Keller
generalized cluster category `C_{(Q,W)}` — for the acyclic/Dynkin (equivalently:
representation-finite hereditary) case the whole cluster category is a **finite fundamental
domain** `ind(mod kQ) ⊔ {P_v[1]}` that shipped machinery already reaches, cluster-tilting
objects **are** support τ-tilting pairs (Adachi–Iyama–Reiten), the cluster-tilted
End-algebra **is** a Jacobian algebra (Buan–Marsh–Reiten / Amiot), and the 2-Calabi–Yau
duality reduces to the Auslander–Reiten formula. The dg/Ginzburg machinery, `D^b(Γ)`, and
the general non-acyclic `C_{(Q,W)}` are **out of scope, recorded** (verification page + a
P80 ledger entry).

---

## Status of the trees at spec time (read this first)

- **This worktree** was checked out behind `dev` (it did not contain `derived/silting.py`,
  the metaplan, or the research doc); the doc is a new file and is version-independent, so
  no fast-forward is needed. The implementer branches `plan-79-…` off `dev`.
- **The venv import resolves to the SHARED checkout** `/Users/marco/.../quiverlab/src/
  quiverlab/__init__.py`, which tracks **`dev`**. At spec time `dev` was **≈ Plan 74**;
  **`dev` is a MOVING TARGET** — P65/P67-adjacent work is landing concurrently — so this
  doc pins **no commit hash** and every scope claim is written as a **runtime property of
  the shipped machinery at compute time**, not a fact about a particular tree tip. Live
  numbers below were computed with `.venv/bin/python` against the dev tree at spec time;
  the acceptance tests re-pin every one against the branch tip.
- **P45 `mutate` root-cause dependency on P65.** The P45 exchange-graph `mutate` can raise a
  *spurious* `status="error"` on some support-τ-tilting pairs (the D₄ star). P65's "Task 0"
  fixes this root cause; **P65 is merging to `dev` concurrently with this plan**. This plan
  does NOT touch `mutation.py`; it reads the exchange graph through the **merged P63
  n-regularity recovery** (`tautilting/wallchamber.py::_closed_by_n_regularity`). The
  recovery certifies completeness only when the discovered graph **is n-regular** — which is
  a **runtime test, not a type guarantee**: live it holds for D₄ (recovers 50) but **NOT for
  D₅** (see the spike — the graph is not n-regular pre-P65, so the count is honestly refused,
  not fabricated). After P65's root fix lands, `mutate` is expected to stop raising and such
  types certify *natively* (`status="complete"`); the gate logic in this plan is unchanged
  either way. The P65 dependency is recorded in the ledger.

---

## Record (verbatim from `docs/plans/2026-08-06-computability-expansion-deep-research.md`, R31)

> **R31 — Amiot–Keller cluster categories from (Q,W) (stretch).** [D-critic omission;
> exploratory]
> Object: the generalized cluster category `C_{(Q,W)}` for the Jacobi-finite
> quivers-with-potential quiverlab already builds (surfaces → gentle Jacobians);
> cluster-tilting objects, 2-CY verification. Heavier dg machinery; certified scope to be
> established at spec time — candidate for a research-grade plan after R30. Size XL.

**Card (metaplan §5 P79):** `C_{(Q,W)}` for the Jacobi-finite `(Q,W)` quiverlab already
builds: the certified slice = cluster-tilting object verification + 2-CY checks on bounded
windows; the full dg machinery explicitly out of scope, recorded. Scope frozen in its plan
doc after a feasibility spike. GUI: only if the slice lands a stable kind; else ledger
entry.

**Card §1.4 (binding):** "R31 ships as a certified verifier-grade slice — scope fixed at
its plan spec, not open-ended."

---

## Reference re-verification (done live at authoring, 2026-08-08; all BibTeX-verifiable)

Every citation below was verified against its arXiv abstract / journal record this session.
The card flagged **Amiot 0805.1035** and **AIR 1210.1036** specifically — both confirmed.

| key (new unless noted) | ref | id | status |
|---|---|---|---|
| `amiot_cluster_category` | Amiot, *Cluster categories for algebras of global dimension 2 and quivers with potential*, **Ann. Inst. Fourier 59(6) (2009), 2525–2590** | arXiv:**0805.1035** | **VERIFIED**; abstract states the load-bearing theorem verbatim: *"when it is Jacobi-finite, [`C_{(Q,W)}`] is endowed with a cluster-tilting object whose endomorphism algebra is isomorphic to the Jacobian algebra."* |
| `bmrrt_cluster` | Buan–Marsh–Reineke–Reiten–Todorov, *Tilting theory and cluster combinatorics*, **Adv. Math. 204(2) (2006), 572–618**, doi 10.1016/j.aim.2005.06.003 | arXiv:math/0402054 | **VERIFIED**; the cluster category `D^b(kQ)/τ⁻¹[1]` + the fundamental-domain / module-category model. |
| `bmr_cluster_tilted` | Buan–Marsh–Reiten, *Cluster-tilted algebras*, **Trans. AMS 359 (2007), 323–332** | arXiv:math/0402075 | **VERIFIED**; `End_{C_Q}(T)` = the cluster-tilted algebra. |
| `keller_reiten_gorenstein` | Keller–Reiten, *Cluster-tilted algebras are Gorenstein and stably Calabi–Yau*, **Adv. Math. 211 (2007), 123–151** | arXiv:math/0512471 | **VERIFIED**; abstract: *"cluster-tilted algebras are Gorenstein of dimension at most one, and hereditary if they are of finite global dimension."* — the Gorenstein-≤1 self-cert. |
| `air_tau_tilting` (**EXISTS**, `AIR2014`) | Adachi–Iyama–Reiten, *τ-tilting theory*, **Compos. Math. 150(3) (2014), 415–452**, doi 10.1112/S0010437X13007422 | arXiv:**1210.1036** | **VERIFIED**; the support-τ-tilting ↔ cluster-tilting bijection backbone. Already in `references.bib` (DOI only; the card's arXiv id confirmed, may be added as a `note`). |
| `ginzburg_cy` | Ginzburg, *Calabi–Yau algebras* | arXiv:math/0612139 | **VERIFIED** — the Ginzburg dg algebra `Γ(Q,W)`; **out-of-scope machinery**, cited only in the boundary section. |
| `keller_yang_mutation` | Keller–Yang, *Derived equivalences from mutations of quivers with potential*, **Adv. Math. 226(3) (2011), 2118–2168** | arXiv:0906.0761 | **VERIFIED** — the "Keller" half of "Amiot–Keller"; **out-of-scope machinery**. |
| `derksen_weyman_zelevinsky` (**EXISTS**, `DWZ2008`) | Derksen–Weyman–Zelevinsky, *Quivers with potentials and their representations I*, Selecta Math. 2008 | — | already in `references.bib`/registry (P44). |
| `labardini` (**EXISTS**, `LabardiniFragoso2009`) | Labardini-Fragoso, surface QPs, Proc. LMS 2009 | — | already in `references.bib`/registry (P44/P48). |

**Non-load-bearing note:** the metaplan card names "Amiot arXiv:2107.02646", but that is the
**skew-gentle** Amiot paper (already `Amiot2021skewgentle`, used by P68), NOT the cluster
category paper. The correct cluster-category reference is **0805.1035** (this session's
correction — folded back to the research doc R31 with a dated note at merge).

---

## The feasibility spike (run live at spec time — the scope-freezing evidence)

Every value below was recomputed in the dev venv this session (`QQ` unless noted). The
scripts are reproduced by the acceptance tests.

### Spike 1 — the module-category model is a finite fundamental domain quiverlab already reaches

`C_Q = D^b(kQ)/τ⁻¹[1]` for `Q` acyclic (BMRRT): a fundamental domain for `F = τ⁻¹[1]` is
`ind(mod kQ) ⊔ {P_v[1] : v ∈ Q_0}`, so

```
#indec(C_Q) = #ind(mod kQ) + n = (#positive roots) + n = #almost-positive-roots = n(n+3)/2  (type A_n).
```

`modules/ar.py::knit_ar_quiver` gives `#ind(mod kQ)`; `tautilting.exchange_graph` gives the
cluster-tilting count. Live:

| `Q` | `dim kQ` | `#ind(mod kQ)` (knit) | `+n` | `#indec(C_Q)` | closed form | `#cluster-tilting` (exchange_graph) | expected |
|---|---|---|---|---|---|---|---|
| `A₂` | 3 | **3** (posroots) | +2 | **5** | `n(n+3)/2 = 5` | **5**, `status=complete`, `n_regular` | Catalan `C₃ = 5` |
| `A₃` | 6 | **6** | +3 | **9** | `9` | **14**, `complete` | Catalan `C₄ = 14` |
| `A₄` | 10 | **10** | +4 | **14** | `14` | **42**, `complete` | Catalan `C₅ = 42` |
| `A₅` | 15 | **15** | +5 | **20** | `20` | **132**, `complete` (≈84 s) | Catalan `C₆ = 132` |
| `D₄` | — | **12** (posroots, `complete`) | +4 | **16** | `n(n-1)+n = 16` | **50**, `status=error` → **recovered complete** (n-regular) | cluster number `D₄ = 50` (Buan–Marsh) |
| `D₅` | — | **20** (posroots, `complete`) | +5 | **25** | `n(n-1)+n = 25` | **182 discovered**, `status=error`, `_closed_by_n_regularity = False`, `complete = False` → **NOT certified → honest refusal** (≈284 s) | cluster number `D₅ = 182` (uncertifiable pre-P65) |

Both the `#sτt` = 5 / 14 values the card asked to cross-check against, and the D₄ = 50
Buan–Marsh number the P63 chamber count already pins, reproduce exactly. `#indec(C_Q)`
(a mod-`kQ` quantity) is certified for D₄ **and D₅** (the knit closes); the **cluster-tilting
count** is a *different* certificate that D₅ **fails** pre-P65 (see Spike 2). **⟹ the acyclic
cluster category's `#indec` and the certified cluster-tilting count are computable through
shipped finite machinery; the certification is a per-input runtime test, not a type list.**

### Spike 2 — cluster-tilting objects ARE support τ-tilting pairs (AIR); certification is a RUNTIME property

`exchange_graph(kQ)` **is** the cluster exchange graph: a support τ-tilting pair `(M,P)`
is the cluster-tilting object `T = M ⊕ P[1]` (AIR; every cluster-tilting object of `C_Q`
decomposes uniquely so). The counts match the cluster (associahedron) numbers exactly
(Spike 1). Whether the count is **certified** is a per-input runtime fact, live-characterized:

- **Type A₂–A₅** close cleanly: `status="complete"`, `n_regular=True`.
- **Type D₄**: the BFS discovers all **50** vertices but **2 of 200** `(vertex,direction)`
  mutations hit the P45 2-term-silting-cone boundary (`mutate: no valid exchange found …`),
  so the raw `status="error"`. The **P63 recovery** certifies completeness anyway BECAUSE
  the discovered graph **is n-regular**: `_closed_by_n_regularity(eg, 4) == True`, and
  `wall_chamber_structure(kD₄)` returns `complete=True` with `50` chambers and the honest
  provenance note. (Every edge is an involutive mutation rediscovered from its other
  endpoint, so an n-regular non-budget-capped graph is closed — a single spurious mutation
  failure is harmless.)
- **Type D₅ — the recovery FAILS** (critic's live run, adjudicated valid): the BFS discovers
  all **182** vertices but the graph is **NOT n-regular** after the P45 mutate boundary bites,
  so `_closed_by_n_regularity(eg, 5) == False`, `complete == False`. The cluster-tilting
  count is therefore **honestly refused** (G4), NOT reported as 182. This is the *natural*
  non-recoverable case — not a constructed hypothetical.

**Design consequence (frozen — a runtime property, NOT a type list):** the cluster-tilting
count is **certified** iff, *at compute time*, `status=="complete"` **OR** (`status=="error"`
**AND** `_closed_by_n_regularity`). Anything else — `status=="budget"`, or a non-n-regular
`"error"` — is **refused loudly** (G3/G4), never a fabricated finite number. With the
current pre-P65 `mutate`, this empirically certifies **A₂–A₆ + D₄** (the A_n that fit the
default budget; see Spike 2b) and refuses **D₅+** and the larger exceptional types; **after
P65's Task-0 root fix lands** (concurrently — the dependency is stated), `mutate` is expected
to stop raising so those types certify **natively** as `status="complete"` — and **the gate
logic here is unchanged either way**. The plan never claims a number the recovery could not
certify.

### Spike 2b — `status=="budget"` is NOT `rep-infinite` (critic-found; adjudicated valid)

`exchange_graph` returns `status="budget"` purely on `len(records) >= budget_pairs`
(default `512`). The cluster numbers **A₇=1430, D₆=672, E₆=833, E₇=4160, E₈=25080** all
exceed 512, so a naive "budget ⟹ infinite" reading would slander **representation-finite
Dynkin** algebras with a **false infiniteness claim**. But `budget` and `infinite` are
distinct facts, and the **Dynkin-type certificate ships** (P60's `dynkin_type` surface): a
`"budget"` stop on a **known-finite Dynkin** input means "rep-finite, but the cluster number
exceeds the pair budget — raise `budget_pairs`", with the exact cluster number stated where
the closed form is known (`C_{n+1}` for A_n; the Buan–Marsh/generalized-Catalan count for
D/E). A genuine **infiniteness** claim is made ONLY off the type certificate — an affine
(`Ã_n`) or wild quiver, where `mod kQ` truly has infinitely many indecomposables. **⟹ the
budget gate MUST split (G3a over-budget-finite vs G3b infinite) on the Dynkin-type
certificate; it never conflates the two.**

### Spike 3 — the cluster-tilted End-algebra is a Jacobian algebra (BMR/Amiot), verified instance-wise

Live: `Jac(1→2→3→1, αβγ)` (`families/jacobian.py`) `= kZ₃/J²` — `dim = 6`,
`cartan_matrix()` equal to the independent `NakayamaAlgebra(n=3, ℓ=2, cyclic=True)`,
`is_selfinjective == True`, `global_dimension ≥ 32` (certified infinite — consistent with
Keller–Reiten: cluster-tilted ⟹ *hereditary iff finite gl.dim*, and this one is not
hereditary). The Fomin–Zelevinsky **quiver** route is shipped and certified:
`matrix_mutation(exchange_matrix(kA₃), 1)` produces an oriented 3-cycle
(`[[0,-1,1],[1,0,-1],[-1,1,0]]`), i.e. the cluster-tilted quiver of that mutation; equipping
it with the type-A canonical potential `W = Σ(oriented 3-cycles)` and building `Jac`
reproduces `kZ₃/J²`. **⟹ the End-algebra construction (FZ quiver + Jacobian) and its
presentation verification are computable.** The *direct* orbit-category `Hom_{C}(T,T)` is
**not** (out of scope) — the equality `End_C(T) = Jac(Q_T,W_T)` is **cited** (BMR/Amiot) and
the target Jacobian's invariants are **verified**.

**Not-yet-spiked, spec'd in Task 3 (critic gap, valid):** the flagship is a *single*
oriented 3-cycle, so it does NOT exercise `_quiver_from_exchange_matrix`'s all-cycles
enumeration. Task 3 adds a **multi-3-cycle type-A instance** — an `A₄`/`A₅` mutation sequence
whose cluster-tilted quiver carries **two or more** oriented 3-cycles — verifying the full
sum-of-3-cycles potential path (dim + the FZ certificate) once the moving tree settles.

### Spike 4 — the 2-CY certificate reduces to the Auslander–Reiten formula (computable today)

For `X,Y ∈ mod kQ`, in `C_Q` one has (BMRRT) `Ext¹_C(X,Y) ≅ Ext¹_A(X,Y) ⊕ D Ext¹_A(Y,X)`,
so `dim Ext¹_C(X,Y) = dim Ext¹_A(X,Y) + dim Ext¹_A(Y,X)` — **manifestly symmetric**, which
is exactly the 2-Calabi–Yau duality `Ext¹_C(X,Y) ≅ D Ext¹_C(Y,X)`. The **genuine computable
content** is the **Auslander–Reiten formula** (Assem–Simson–Skowroński, *Elements…* Vol. 1,
**Thm IV.2.13**): `Ext¹_A(X,Y) ≅ D\overline{Hom}_A(Y, τX) ≅ D\underline{Hom}_A(τ⁻¹Y, X)`,
where `\overline{Hom}` = Hom modulo maps factoring through injectives and `\underline{Hom}` =
Hom modulo maps factoring through projectives. The certificate uses **ordinary** `Hom(Y,τX)`;
the theorem is about the **injectively-stable** `\overline{Hom}`. On a **representation-finite
hereditary** algebra the two coincide across the tested window — this is the **theorem-backed
proxy justification**, and the evidence is a live sweep with **0 mismatches**:

- `dim Ext¹_A(X,Y) == dim Hom_A(Y, τX)` held on **every ordered pair** across
  **A₂–A₅, D₄, and a zigzag-oriented A₄** (the critic's expanded run — MORE robust than the
  first draft's single `A₃` claim; 0 mismatches everywhere). On `A₃` the 5 nonzero-Ext pairs
  are `(P₂,I₂),(P₂,S₁),(S₃,I₂),(S₃,S₂),(S₂,S₁)`, each `Ext¹ = 1`.
- the symmetric assembly `dim Ext¹_C(X,Y) = dim Ext¹_C(Y,X)` held on all pairs.

The certificate payload records BOTH the pointwise AR result AND that the symmetry of
`dim Ext¹_C` is *by construction* (`a+b=b+a`) — the AR formula is the content, the symmetry
is bookkeeping; the payload says so. **⟹ the 2-CY Ext-symmetry certificate is computable on
the module window** (the shifted `P_v[1]` pairs are handled by citation — see the MODERATE
ruling in the API/Task 4 spec).

### What does NOT compute (the recorded out-of-scope, confirmed by the spike)

- The **Ginzburg dg algebra** `Γ(Q,W)`, `D^b(Γ)`, the derived orbit category's triangulated
  structure, and the **direct** `Hom_{C_{(Q,W)}}` as `⊕_i Hom_{D^b}(X, F^iY)` — no dg engine
  ships (`derived/` is `K^b(proj)` only, per P43/P67; the spike confirmed there is no
  `D^b(Γ)`). We instead use the **finite mod-`A` reductions** (BMRRT), valid **only** in the
  acyclic/hereditary case.
- **General Jacobi-finite non-acyclic `(Q,W)`**: when `(Q,W)` is *not* mutation-equivalent to
  an acyclic quiver, `C_{(Q,W)}` is not a hereditary cluster category and there is no
  finite mod-`A` model — the category-level invariants (`#indec`, exchange graph of
  `C_{(Q,W)}`) are **refused loudly**. Only the **Jacobian algebra** `Jac(Q,W)`'s own module
  theory (served by the rest of quiverlab) and the **algebra-level cluster-tilted
  certificate** (finite-dim + Gorenstein dim ≤ 1, Keller–Reiten) stand.
- **DWZ potential mutation / right-equivalence** (inferring `W_T` for a general mutated
  quiver): already deferred at P48.1 (`surfaces/flip.py`: "DWZ potential right-equivalence
  is deferred"). Out of scope here too; the **quiver** mutation is certified, the
  **potential** is constructed only where it is a sum of 3-cycles (type A), else honest
  partial.
- The **P45 `mutate` root-cause** (D₄ spurious `status="error"`): owned by **P65 Task 0**,
  not this plan. We use the P63 recovery and record the dependency.

---

## Frozen scope (the v1 certified slice — contractual)

**IN (certified, finite, shipped machinery):**

1. **The acyclic / representation-finite hereditary cluster category `C_Q`** (input: a
   Dynkin type string, an acyclic `Quiver`, or a hereditary `Algebra = kQ`, `I = 0`):
   - `#indec(C_Q) = #ind(mod kQ) + n` via `knit_ar_quiver` (certified when the knit closes —
     holds for D₄ **and** D₅); the indecomposables listed as `ind(mod kQ) ⊔ {P_v[1]}` (the
     almost-positive-roots model).
   - `#cluster-tilting objects` + the cluster exchange graph = `tautilting.exchange_graph`,
     with the certification a **per-input RUNTIME property** (Spike 2): certified iff, at
     compute time, `status=="complete"` **OR** (`status=="error"` **AND** the P63
     `_closed_by_n_regularity`); otherwise **refused loudly** (G3/G4 — D₅ is refused pre-P65,
     never reported as 182). With pre-P65 `mutate` this certifies **A₂–A₆ + D₄**; after
     P65's root fix, larger types certify natively — gate logic unchanged.
   - cluster-tilting objects **↔ support τ-tilting pairs** `(M,P)`, shown as `T = M ⊕ P[1]`
     (the **AIR bijection — the plan's backbone and its cross-engine oracle**).
2. **The cluster-tilted End-algebra** `End_{C_Q}(T) = Jac(Q_T, W_T)` (BMR/Amiot): the
   cluster-tilted **quiver** via shipped FZ matrix mutation (`surfaces/flip.py`,
   per-instance certified), the **algebra** as `JacobianAlgebra` for **type A** (canonical
   potential = sum of oriented 3-cycles — verified on BOTH a single-3-cycle flagship AND a
   multi-3-cycle A₄/A₅ instance, Task 3) with **presentation verification** (`dim`, Cartan,
   self-injective/Gorenstein); other types → certified quiver + honest potential-deferral.
3. **The 2-Calabi–Yau certificate** on the module window: the AR-duality anchor
   `dim Ext¹_A(X,Y) = dim Hom_A(Y, τX)` (ASS Thm IV.2.13; ordinary=stable proxy justified by
   the 0-mismatch A₂–A₅/D₄/zigzag-A₄ sweep) across the AR quiver + the symmetric `Ext¹_C`
   assembly (self-certifying). The certificate is **scoped to the module block**: the shifted
   `P_v[1]` pairs are asserted **by citation** (BMRRT 2-CY), NOT computed (see the MODERATE
   ruling) — the payload is honest that `is_2_calabi_yau` is a module-window verdict.
4. **Jacobi-finite `(Q,W)` input** (`from_potential`): build `Jac(Q,W)` (or the shipped
   loud `NotFiniteDimensionalError`); **certify the algebra-level Amiot statement** — it is
   finite-dimensional and Gorenstein of dimension ≤ 1 (Keller–Reiten). The **category
   combinatorics** are served **only** when `(Q,W)` is certified reducible to a hereditary
   model (the flagship `(3-cycle, αβγ) ≅ C_{A₃}`), else the category invariants are refused.

**OUT (recorded — verification page honest-scope + P80 ledger entry):** the dg/Ginzburg
machinery `Γ(Q,W)`/`D^b(Γ)`; direct orbit-category `Hom_C`; general non-acyclic
`C_{(Q,W)}` category invariants; DWZ potential mutation/right-equivalence; rep-infinite
(tame/wild) hereditary `#indec`/exchange graph (infinite — bounded-window refusal); the P45
`mutate` root-cause (P65).

**GUI decision (frozen from the spike):** the slice **lands a stable `cluster_category`
kind** — `#indec`, cluster-tilting count (with the recovery flag), a verified `End(T)`
presentation, the 2-CY certificate — so **GUI is IN** (per the card). No ledger deferral of
GUI; only the OUT items above go to the ledger.

---

## Mathematical foundation (the plan's ground truth)

Let `k` be a field, `Q` a finite acyclic quiver, `A = kQ` hereditary, `n = #Q_0`.

**(F1) The cluster category (BMRRT).** `C_Q := D^b(mod A)/F`, `F = τ_{D}⁻¹∘[1]` (the orbit
category), is a Hom-finite triangulated 2-Calabi–Yau category with suspension `[1]`. A
**fundamental domain** for `F` is `ind(mod A) ∪ (add A)[1]`; hence the indecomposables of
`C_Q` are `ind(mod A) ⊔ {P_v[1] : v∈Q_0}` and
`#indec(C_Q) = #ind(mod A) + n`. For `A` representation-finite this is finite and equals the
number of **almost-positive roots** of the underlying Dynkin diagram (`= n(n+3)/2` for
`A_n`, `= n(n-1)+n = n²` for `D_n`).

**(F2) Cluster-tilting objects ↔ support τ-tilting (AIR).** `T ∈ C_Q` is cluster-tilting iff
`Ext¹_C(T,T)=0` and `T` is maximal so. There is a bijection between basic cluster-tilting
objects of `C_Q` and support τ-tilting `A`-modules: `T ↦ (M,P)` where `T = M ⊕ P[1]`,
`M∈mod A`, `P∈proj A`, and `(M,P)` is a support τ-tilting pair. **Mutation of `T` = mutation
of `(M,P)`**, so `exchange_graph(A)` **is** the cluster exchange graph (associahedron for
`A_n`). Count = the cluster / generalized Catalan number.

**(F3) Cluster-tilted algebras (BMR / Amiot).** `B := End_{C_Q}(T)^{op}` is the
*cluster-tilted algebra* of `T`. Its quiver `Q_T` is the Fomin–Zelevinsky mutation of `Q`
along the mutation sequence reaching `T`. **Amiot:** for any Jacobi-finite quiver with
potential `(Q',W')`, the generalized cluster category `C_{(Q',W')}` has a cluster-tilting
object with `End = Jac(Q',W')`; specializing, `B = Jac(Q_T, W_T)`. For **type A** the
potential is the sum of the oriented 3-cycles of `Q_T` (a known instance). **Keller–Reiten:**
every cluster-tilted algebra is Gorenstein of dimension ≤ 1, and hereditary iff of finite
global dimension. Flagship: the once-mutated `A₃` gives `Q_T =` the 3-cycle and
`B = Jac(3-cycle, αβγ) = kZ₃/J²` (dim 6, self-injective).

**(F4) 2-Calabi–Yau ⟹ AR duality.** `C_Q` is 2-CY: `Ext¹_C(X,Y) ≅ D Ext¹_C(Y,X)` naturally.
On modules, `Ext¹_C(X,Y) ≅ Ext¹_A(X,Y) ⊕ D Ext¹_A(Y,X)`, and the classical **Auslander–
Reiten formula** (ASS *Elements…* Vol. 1, **Thm IV.2.13**)
`Ext¹_A(X,Y) ≅ D\overline{Hom}_A(Y, τX) ≅ D\underline{Hom}_A(τ⁻¹Y, X)` is the computable
anchor (`\overline{Hom}` = mod injectively-factoring, `\underline{Hom}` = mod
projectively-factoring). The symmetry of `dim Ext¹_C` is by construction (`a+b=b+a`); the
*non-trivial* verified content is the pointwise AR formula holding across the AR quiver
(ordinary Hom used as the theorem-backed proxy for `\overline{Hom}` — they coincide on the
tested rep-finite hereditary window, 0 mismatches). The two together certify the 2-CY
structure **on the module window**; the shifted-object pairs are asserted by citation.

---

## Live-verified facts (recomputed in the dev venv at spec time, `QQ`)

1. `#cluster-tilting(C_{A_n}) = Catalan(n+1)`: `A₂=5, A₃=14, A₄=42, A₅=132` — `= exchange_graph`
   vertex count, `status=complete`, `n_regular=True` (A₅ ≈84 s). **[xeng + lit]**
2. `#indec(C_{A_n}) = n(n+3)/2`: `A₂=5, A₃=9, A₄=14, A₅=20` — `knit_ar_quiver` count `+ n`.
   `#indec(C_{D₄}) = 16` (`12` posroots `+4`); `#indec(C_{D₅}) = 25` (`20` posroots `+5`),
   knit `complete`. **[lit]**
3. `#cluster-tilting(C_{D₄}) = 50` (Buan–Marsh) — `exchange_graph` discovers 50, raw
   `status=error`, **`_closed_by_n_regularity == True`** (n-regular) ⟹ `wall_chamber_structure`
   reports `complete=True, 50 chambers` with the honest recovery note. **[xeng + lit + selfcert]**
4. **`#cluster-tilting(C_{D₅})` REFUSED (not 182)** — `exchange_graph` discovers 182 vertices
   but `status=error`, **`_closed_by_n_regularity(eg,5) == False`**, `complete=False` (not
   n-regular pre-P65; ≈284 s). The gate refuses loudly; it never emits 182 as certified. The
   D₅ test asserts *certified-or-refused* in BOTH regimes (post-P65 it is expected to become
   `status=complete`). **[selfcert]**
5. `Jac(3-cycle, αβγ) = kZ₃/J²`: `dim=6`, `cartan == NakayamaAlgebra(3,2,cyclic)`,
   `is_selfinjective=True`, `global_dimension ≥ 32` (infinite; consistent with Keller–Reiten
   contrapositive). `matrix_mutation(exchange_matrix(kA₃),1)` contains an oriented 3-cycle.
   **[xeng + lit + selfcert]**
6. 2-CY: `dim Ext¹_A(X,Y) == dim Hom_A(Y, τX)` on **every ordered pair** across **A₂–A₅, D₄,
   and zigzag-A₄** (0 mismatches — the critic's expanded sweep, more robust than a single-A₃
   claim); symmetric `Ext¹_C` on all pairs. **[selfcert]**

---

## Scope gates (contractual — every boundary is a loud typed `QuiverlabError`, never a silent gap)

| # | Condition | Behaviour |
|---|---|---|
| G1 | input algebra not hereditary (`I ≠ 0`) and not a recognized reducible `(Q,W)` | category invariants **refused** (`cluster category model needs an acyclic / hereditary A or a certified-reducible (Q,W)`); the `Jac` algebra-level certificate is still offered for a `from_potential` call |
| G2 | hereditary `A = kQ` with `Q` NOT acyclic (has an oriented cycle) | **refused** (`kQ is not finite-dimensional / C_Q undefined for a quiver with cycles`) |
| G3a | `exchange_graph` `status="budget"` on a **known-finite Dynkin** input (P60 `dynkin_type` certifies A/D/E; cluster number > `budget_pairs`, e.g. A₇=1430, D₆=672, E₆=833) | cluster-tilting count **refused as OVER-BUDGET, NOT infinite** (`rep-finite (Dynkin <type>) but the cluster number <count-if-known> exceeds budget_pairs=<b> — raise the budget`); the exact cluster number is quoted where the closed form is known. **Never** an infiniteness claim. |
| G3b | acyclic `Q` genuinely **representation-infinite** (`dynkin_type` says affine `Ã`/wild) | `#indec` / full combinatorics **refused as INFINITE** (`C_Q has infinitely many indecomposables; bounded-window enumeration only`), claimed ONLY off the type certificate |
| G4 | `exchange_graph` `status="error"` and **NOT** `_closed_by_n_regularity` (the **D₅** case, live) | cluster-tilting count **refused** (`exchange graph did not close and n-regularity recovery failed — raise a P65-fixed mutate or wait for the root fix`); NEVER read as a fabricated finite count (D₅'s 182 is never emitted) |
| G5 | `cluster_tilted_algebra` for a mutated quiver whose potential is **not** a sum of 3-cycles (type D/E/general) | the **quiver** is returned (FZ-certified) with `algebra=None` + a note (`potential inference / DWZ right-equivalence is out of scope, P48.1`); NEVER a guessed potential |
| G6 | `from_potential(Q,W)` with `Jac(Q,W)` **infinite** | shipped `NotFiniteDimensionalError` propagates (Jacobi-infinite refusal) |
| G7 | general non-acyclic Jacobi-finite `(Q,W)`, category invariants requested | **refused** (`C_{(Q,W)} category invariants need D^b(Γ) — out of scope; only the Jacobian-algebra Amiot certificate is available`) |
| G8 | anything needing the Ginzburg dg algebra / `D^b(Γ)` / direct orbit-category `Hom_C` | **refused**, pointing at the recorded out-of-scope + the P80 ledger |

---

## API surface (public via `import quiverlab`; exact only)

```python
# src/quiverlab/cluster/category.py
class ClusterCategory:
    """The Amiot–Keller cluster category, certified module-category slice.
    Construct from a Dynkin type / acyclic Quiver / hereditary Algebra (the classical
    C_Q), or via `from_potential` (a Jacobi-finite (Q,W), algebra-level certificate)."""
    def __init__(self, source, field=None): ...          # normalizes to (A=kQ acyclic, n)
    @classmethod
    def from_potential(cls, Q, W, field=None): ...        # Jac(Q,W); category iff reducible

    n: int
    hereditary: bool
    dynkin_type: str | None                               # None if acyclic-but-not-Dynkin

    def num_indecomposables(self) -> int: ...             # #ind(mod A) + n  (G3b refusal)
    def indecomposables(self) -> list:                    # [("module", M) | ("shift", P_v[1])]
    def almost_positive_roots_count(self) -> int: ...     # = num_indecomposables (identity check)

    def exchange_graph(self):                             # the cluster exchange graph (P45)
    def num_cluster_tilting(self):                        # -> {"count"|None, "certified", "status", "note", "dynkin_type"}
    def cluster_tilting_objects(self) -> list:            # [{"pair": (M,P), "as_object": "M ⊕ P[1]"}]

    def cluster_tilted_algebra(self, mutation_seq):       # -> {"quiver", "algebra"|None, "verified", "note"}
    def two_cy_certificate(self):                         # -> {"ar_formula_holds", "ext1_C": {...}, "pairs_checked", "scope"}
    def is_2_calabi_yau(self):                            # -> {"verdict": bool, "scope": "module_window", "shifted_by_citation": True}

# core/algebra.py delegate (thin lazy-import, hereditary only):
#   Algebra.cluster_category(self) -> ClusterCategory(self)
```

- `num_cluster_tilting()` returns a **dict** so the recovery provenance and the
  budget-vs-infinite distinction are never lost. Three shapes:
  - **certified:** `{"count": 50, "certified": True, "status": "error"|"complete", "note":
    "…recovered by n-regularity (Plan 63)…" | "closed natively"}`.
  - **over-budget (G3a):** `{"count": None, "certified": False, "status": "budget",
    "dynkin_type": "D6", "note": "rep-finite (Dynkin D6); cluster number 672 exceeds
    budget_pairs=512 — raise the budget"}` — the count field is the KNOWN closed form when
    available, never a truncated discovered count, and **never** an infiniteness claim.
  - **not recoverable (G4, D₅) / infinite (G3b):** `{"count": None, "certified": False,
    "status": "error"|"budget", "note": "<the honest refusal>"}`.
- `is_2_calabi_yau()` returns a **scoped verdict**, NOT a bald categorical `True` (MODERATE
  ruling): `{"verdict": True, "scope": "module_window", "shifted_by_citation": True, "note":
  "2-CY verified on the module block via the AR formula (ASS IV.2.13) across the AR quiver;
  the shifted P_v[1] pairs hold by BMRRT (cited, not computed)"}`. The `P_v[1]`-pairs are
  never asserted as computed.
- `cluster_tilted_algebra` returns the **quiver** always (FZ-certified) and the **algebra**
  only where the potential is constructible (G5); `verified` records the presentation checks
  actually run (`dim`, `cartan`, `self_injective`).

---

## Global constraints (binding)

- **No floats in `src/`** — every payload exact (`int`/`Fraction`/domain elements). The AST
  gate `tests/test_no_floats.py` scans `src/`; the new `cluster/` package is in scope.
- **Engine internals stay internal.** `cluster/` composes the *public* surfaces
  (`tautilting`, `modules.ar`, `modules.ext`, `families.jacobian`, `surfaces.flip`); it does
  NOT reach into `engine.*`.
- **Reuse, never re-implement:** exchange graph = `tautilting.exchange_graph` +
  `wallchamber._closed_by_n_regularity`; τ/Ext/Hom = `modules` (P05/P23/P24); Jacobian =
  `families.jacobian`; FZ mutation = `surfaces.flip.matrix_mutation`/`exchange_matrix`. This
  plan writes **glue + certificates**, not new homological algebra.
- **`mutation.py` is NOT touched** (the P45 root-cause is P65 Task 0). The n-regularity
  recovery is read-only.
- **Oracle-class markers per Plan 32** on every test; new oracles land on
  `docs/verification.md` at acceptance (Plan-22 standing rule).
- Conventional commits; green tests at every commit; deep-bucket by directory
  (`tests/families/`, `tests/modules/` → deep; new `tests/cluster/` → **add to the deep
  bucket in `tests/conftest.py`**).

---

## Tasks (TDD — failing test first, then implement, then gate)

### Task 0 — dependency + reuse audit (no code)
- [ ] **Re-probe the moving tree at branch-cut** (P65 may have landed): record whether D₄/D₅
  `exchange_graph` return `status="error"` or (post-P65) `status="complete"`. Do NOT edit
  `tautilting/mutation.py` regardless; the plan reads `_closed_by_n_regularity`. Whatever the
  tree tip, the gate logic (Spike-2 runtime property) is unchanged — record the observed
  regime so the D₅ test (Task 2) picks the right expectation branch.

### Task 1 — `cluster/category.py`: the fundamental-domain model
`tests/cluster/test_category_indec.py` (deep):
- [ ] **Failing tests:** `ClusterCategory("A2").num_indecomposables() == 5`; `A3 == 9`;
  `A4 == 14`; `D4 == 16`; **`D5 == 25`** (`#indec` closes for D₅ even though its cluster-tilting
  count does not); the identity `num_indecomposables() == #ind(mod A) + n ==
  almost_positive_roots_count()`; `indecomposables()` returns exactly `#ind(mod A)` module
  entries + `n` `("shift", P_v[1])` entries; **G2** (cyclic `Q` refused), **G3b** (`Ã₂`
  / a wild input refused loudly as INFINITE off the type certificate — never a truncated
  count). `[oracle_literature`, `oracle_selfcert]`
- [ ] **Implement** input normalization (Dynkin string via `PathAlgebra`/`dynkin_quiver`;
  `Quiver` → `Q.algebra()`; `Algebra` accepted iff hereditary+acyclic — reuse the shipped
  hereditary/acyclic checks), `knit_ar_quiver` wrap (honor its `status`/budget → G3b),
  `dynkin_type` detection (**reuse P60's `dynkin_type` surface** — the same certificate that
  splits G3a/G3b).

### Task 2 — cluster-tilting = support τ-tilting (AIR); the recovery-certified count + the budget split
`tests/cluster/test_cluster_tilting_air.py` (deep):
- [ ] **Failing tests:** `num_cluster_tilting("A2")["count"] == 5` (`certified=True`,
  `status="complete"`); `A3 == 14`; `A4 == 42`; `A5 == 132`; **D4 `count == 50`,
  `certified=True`, `status=="error"`** (recovery) with the provenance note; the
  **cross-engine identity**: `count == len(exchange_graph(A).vertices)` AND (type A) `==
  Catalan(n+1)`; `cluster_tilting_objects()` yields pairs `(M,P)` with `|M|+|P| == n` for
  every object (basic, n-summand).
- [ ] **The D₅ BOTH-REGIMES test (MAJOR 1):** `num_cluster_tilting("D5")` must be
  **certified-or-refused-honestly** — assert `certified is False` with `count is None` and a
  loud note **OR** (post-P65) `certified is True, status=="complete", count==182`; the test
  **fails only** if a number is emitted as certified while `certified is False`, or if the
  refusal message claims infiniteness. A comment records the pre-/post-P65 expectation.
- [ ] **The budget-vs-infinite split test (MAJOR 2 / G3a):** with a small `budget_pairs`
  (e.g. 8) on `A5`, `num_cluster_tilting` returns `status=="budget"`, `certified is False`,
  `dynkin_type=="A5"`, and a note that says **"exceeds budget"** and quotes the cluster
  number (132) — and **never** contains "infinite"; contrast an affine `Ã` input whose note
  **does** claim infiniteness (G3b). `[oracle_crossengine`, `oracle_literature`,
  `oracle_selfcert]`
- [ ] **Implement** via `exchange_graph` + `_closed_by_n_regularity`, gated by the P60
  `dynkin_type` certificate for the budget/infinite split; expose the pairs as
  `T = M ⊕ P[1]` (read `pair.summands` + `pair.support`).

### Task 3 — the cluster-tilted End-algebra `Jac(Q_T,W_T)` (BMR/Amiot), verified
`tests/cluster/test_cluster_tilted_endalg.py` (deep):
- [ ] **Failing tests (single-3-cycle flagship):** `cluster_tilted_algebra` on the `A₃`
  mutation `[1]` returns a quiver equal (up to relabelling) to the 3-cycle **and** `algebra`
  with `dim == 6`, `cartan == NakayamaAlgebra(3,2,cyclic).cartan_matrix()`,
  `is_selfinjective == True`, and Gorenstein-dim ≤ 1 (Keller–Reiten self-cert); the FZ
  **quiver** step certified by `surfaces.flip` (`matrix_mutation` == the rebuilt quiver's
  exchange matrix).
- [ ] **Failing test — MULTI-3-cycle instance (minor b, exercises the all-cycles path):** an
  `A₄`/`A₅` mutation sequence whose cluster-tilted quiver carries **≥ 2** oriented 3-cycles;
  assert `_quiver_from_exchange_matrix` recovers **all** of them, the sum-of-3-cycles
  potential builds a finite `Jac`, `dim` matches the independent cluster-tilted-of-type-A
  count, and the FZ certificate holds. (The specific sequence is chosen at implementation
  once the tree settles — the acceptance is "≥ 2 three-cycles present and all enumerated".)
- [ ] **G5** (a D₄ mutation whose potential is not a sum of 3-cycles returns `algebra=None`
  + the deferral note, quiver present). `[oracle_crossengine`, `oracle_literature`,
  `oracle_selfcert]`
- [ ] **Implement** `_quiver_from_exchange_matrix(B, verts)` (positive entries → arrows),
  the type-A canonical potential `Σ oriented-3-cycles` (**enumerate ALL oriented 3-cycles**,
  not just one), `JacobianAlgebra` build + the verification battery; the End-identification
  `End_C(T)=Jac(Q_T,W_T)` is **cited, not computed** (the docstring says so).

### Task 4 — the 2-Calabi–Yau certificate (AR duality), scoped
`tests/cluster/test_two_cy.py` (deep):
- [ ] **Failing tests:** `two_cy_certificate("A3")["ar_formula_holds"] is True` with
  `pairs_checked == 36`; the same across **A₂–A₅, D₄, and zigzag-A₄** (0 mismatches — the
  expanded sweep is the recorded evidence for the ordinary=stable proxy); the returned
  `ext1_C` table symmetric (`ext1_C[X,Y]==ext1_C[Y,X]`) and equal to `e_XY + e_YX` on every
  pair; the 5 nonzero-`Ext¹_A` pairs of `A₃` reproduced; a **drift gate** (`ext1_C` dims
  agree with `modules.ext` recomputed independently).
- [ ] **MODERATE ruling — the verdict is SCOPED:** `is_2_calabi_yau()` returns
  `{"verdict": True, "scope": "module_window", "shifted_by_citation": True, …}`, NOT a bald
  `True`; the test asserts the `scope`/`shifted_by_citation` fields are present and that the
  payload never claims the `P_v[1]` pairs were computed. `[oracle_selfcert`,
  `oracle_crossengine]`
- [ ] **Implement** over `knit_ar_quiver` modules using `ext_dims`, `M.tau()`, `hom_space`;
  the certificate records BOTH the AR-formula pointwise result (ASS IV.2.13, ordinary Hom as
  the theorem-backed proxy for `\overline{Hom}`) and the by-construction symmetric assembly,
  and the honest module-window scope for the shifted pairs — all in the payload `note`.

### Task 5 — `from_potential`: the Jacobi-finite algebra-level certificate
`tests/cluster/test_from_potential.py` (deep):
- [ ] **Failing tests:** `ClusterCategory.from_potential(3-cycle, αβγ)` builds
  `Jac = kZ₃/J²`, certifies finite-dim + Gorenstein-dim ≤ 1 (Keller–Reiten) + self-injective;
  the flagship is **recognized reducible** to `C_{A₃}` so the category combinatorics
  (`num_indecomposables`, `num_cluster_tilting`) are served and agree with the `A₃` route;
  **G6** (a Jacobi-infinite `(Q,W)` — one vertex, two loops, `W=0` — raises
  `NotFiniteDimensionalError`); **G7** (a general non-acyclic non-reducible `(Q,W)` refuses
  category invariants while still exposing the algebra certificate). `[oracle_literature`,
  `oracle_selfcert]`
- [ ] **Implement** the `Jac` build + Gorenstein-dim ≤ 1 check (reuse P40 homdims / the
  self-injective certifier); the reducibility recognition kept **narrow and honest** (the
  documented 3-cycle↔A₃ instance; a general mutation-equivalence-to-acyclic test is
  explicitly *not* attempted — noted as future work).

### Task 6 — public export + `Algebra.cluster_category` delegate
- [ ] Export `ClusterCategory` from `src/quiverlab/cluster/__init__.py` and the top-level
  `quiverlab/__init__.py`; add the thin `Algebra.cluster_category()` (hereditary-only,
  lazy import). `tests/cluster/test_public_surface.py` (fast): `from quiverlab import
  ClusterCategory`; `PathAlgebra("A3").cluster_category().num_cluster_tilting()["count"]==14`.

### Task G1 — the `cluster_category` GUI/webapp/HPC kind (all three tiers, i18n ×4)
An **algebra-only, budget-carrying** kind (the `ar_quiver`/`wall_chamber` precedent —
`cluster_category` or `cluster_category:512`, the exchange-graph pair budget; bypasses the
`name:0..N` degree grammar).
- **Files (mirror the `wall_chamber` branch exactly):** `src/quiverlab/hpc/spec.py`
  (grammar parse + `_dispatch` + `_snip`), `docs/gui/runner.py` (the Pyodide twin —
  shape-identical), `webapp/server/schema.py` (the third grammar site),
  `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox `qlgui-cluster-category`,
  `S.ids`, push-list, `renderBlock`, `scheduleProbe` ETA), **`src/quiverlab/hpc/
  estimator.py`** (a `cluster_category` sizing branch — see the estimator note below; NOT a
  bare `sizing_dim`), `webapp/templates/index.html` (one checkbox),
  `webapp/server/i18n/{en,es,fr,zh}.json` (**all four** — block-key chain +
  `pick.kind.cluster_category`, parity gated by `tests/webapp/test_i18n.py`),
  `src/quiverlab/trace/results_html.py` (`_HEADINGS` + a render branch),
  `tests/webapp/_runner_goldens.json` + `test_runner_delegation.py` (ONE golden).
- **Estimator note (MAJOR 2 / risk 5 — the exchange-graph BFS is NOT algebra-dim-bound).**
  A bare `sizing_dim` on algebra dim under-sizes badly: `E₆` has `dim kQ = 36` but its
  exchange graph has 833 vertices and the BFS ran ~8 min in the critic's probe; `A₅` (dim 15)
  ran ≈84 s and `D₅` (—) ≈284 s. The `cluster_category` sizing MUST scale on the **expected
  cluster number** (the P60 `dynkin_type` closed form where known, else `budget_pairs`), not
  algebra dim, and route oversized inputs off the instant tier like a big family. State this
  in the estimator branch + a test that `E₆`/`A₇` size into the queued tier, not instant.
- **Block shape** (one shared `cluster_category_block(A, budget)`):
```python
{"kind": "cluster_category", "n": int, "hereditary": bool, "dynkin_type": str|None,
 "num_indec": int|None,                                   # None iff G3b (infinite)
 "num_cluster_tilting": {"count": int|None, "certified": bool, "status": str,
     "dynkin_type": str|None, "note": str|None} | None,   # count None on budget/error refusal
                                                          # note distinguishes over-budget vs infinite (G3a/G3b/G4)
 "cluster_tilted": {"quiver": [...], "dim": int|None, "self_injective": bool|None,
     "verified": [...], "note": str|None} | None,         # the flagship mutation, type A
 "two_cy": {"ar_formula_holds": bool, "pairs_checked": int,
     "scope": "module_window", "shifted_by_citation": True} | None,   # scoped verdict, not bald
 "note": str|None,
 "references": ["amiot_cluster_category", "bmrrt_cluster", "bmr_cluster_tilted",
                "keller_reiten_gorenstein", "air_tau_tilting", "assem_book"]}
# refusal (non-hereditary / rep-infinite / presentation-less) -> {"error": msg,
#   "references": [...]}, never a 500.  A budget/error stop is NOT an error block: it is a
#   populated num_cluster_tilting with count=None + certified=False + the honest note.
```
- [ ] **Step 1** failing cross-runner tests (unmarked, extras-gated dir): `kA₃` block has
  `num_indec == 9`, `num_cluster_tilting.count == 14` (`certified=True`),
  `cluster_tilted.dim == 6` (the `[1]` mutation), `two_cy.ar_formula_holds is True`,
  `"amiot_cluster_category"` in citation keys; **twin parity**
  (`json.dumps(sort_keys=True)` equality across `hpc/spec.py` and `docs/gui/runner.py`).
- [ ] **Step 2** implement grammar in **all three** sites + `_dispatch` + twin + `gui.js`
  checkbox/ETA + `_snip` (`"cluster_category": "A.cluster_category()"`) + `_HEADINGS
  ["cluster_category"] = "Cluster category"` + the render branch (indec count, the
  recovery-flagged cluster-tilting count table, the cluster-tilted `End(T)` presentation,
  the 2-CY certificate line, stated elisions for the big enumerations).
- [ ] **Step 3** add ONE golden `cluster_category_kA3`; confirm existing goldens
  byte-identical first; **canonical-key stability** (schema-v1 algebra-only, no `module`
  block → key byte-stable; confirm an existing family request's key unchanged).
- [ ] **Step 4** run `tests/webapp/test_cluster_category_kind_p79.py
  tests/webapp/test_runner_delegation.py tests/gui/test_cluster_runner_twin_p79.py
  tests/webapp/test_i18n.py tests/hpc -q`.
- [ ] **Step 5** commit `feat(gui,webapp,hpc): cluster_category algebra-only kind (Plan 79)
  — #indec + recovery-certified cluster-tilting count + cluster-tilted End(T) + 2-CY, both
  runners byte-identical, grammar in 3 sites, i18n x4, one golden`.

**i18n keys — ALL FOUR LOCALES** (mirror `pick.kind.wall_chamber`):

| key | en | es | fr | zh |
|---|---|---|---|---|
| `pick.kind.cluster_category` | Cluster category | Categoría de clúster | Catégorie amassée | 簇范畴 |
| `inv.cluster_category` | Cluster category (Amiot–Keller, acyclic slice) | Categoría de clúster (Amiot–Keller, caso acíclico) | Catégorie amassée (Amiot–Keller, cas acyclique) | 簇范畴 (Amiot–Keller, 无圈情形) |
| `block.cluster_category.title` | Cluster category | Categoría de clúster | Catégorie amassée | 簇范畴 |
| `block.cluster_category.indec` | Indecomposables (almost-positive roots) | Indescomponibles (raíces casi positivas) | Indécomposables (racines presque positives) | 不可分解对象 (几乎正根) |
| `block.cluster_category.tilting` | Cluster-tilting objects | Objetos clúster-basculantes | Objets amas-basculants | 簇倾斜对象 |
| `block.cluster_category.tilted` | Cluster-tilted algebra End(T) | Álgebra clúster-basculada End(T) | Algèbre amas-basculée End(T) | 簇倾斜代数 End(T) |
| `block.cluster_category.twocy` | 2-Calabi–Yau (AR duality) | 2-Calabi–Yau (dualidad AR) | 2-Calabi–Yau (dualité AR) | 2-Calabi–Yau (AR 对偶) |

### Task G2 — verification page, citations, README, metaplan tick, recount
**Files:** `src/quiverlab/citations/references.bib` + `registry.py`; `docs/verification.md`;
`README.md`; `docs/plans/2026-08-07-metaplan-v1.0.0.md` (tick P79 + the R31 arXiv correction
note); the P80 GUI-deferral ledger for the OUT items; existing release gates.
- [ ] **Citations (BibTeX-VERIFIED only, this session):** add `amiot_cluster_category`
  (0805.1035), `bmrrt_cluster` (math/0402054), `bmr_cluster_tilted` (math/0402075),
  `keller_reiten_gorenstein` (math/0512471), `ginzburg_cy` (math/0612139),
  `keller_yang_mutation` (0906.0761); add the `note = {arXiv:1210.1036}` to the existing
  `AIR2014`. Register each in `registry.py` with role tags. `derksen_weyman_zelevinsky` /
  `labardini` reused. Any field that cannot be confirmed is **deleted, not guessed**
  (house rule).
- [ ] **Verification page** — new subsystem block "Cluster categories (Amiot–Keller slice,
  Plan 79)" with the oracle→test map:
  - **oracle_literature:** `#indec = #almost-positive-roots` (A₂=5, A₃=9, A₄=14, A₅=20,
    D₄=16, **D₅=25**); `#cluster-tilting =` Catalan/cluster number (A₂=5, A₃=14, A₄=42,
    A₅=132, D₄=50); `Jac(3-cycle)=kZ₃/J²` dim 6.
  - **oracle_crossengine:** cluster-tilting count `== exchange_graph` (P45) `==` Catalan
    (the AIR bijection); `Jac(3-cycle).cartan == NakayamaAlgebra(3,2,cyclic).cartan`; the
    2-CY `Ext¹_C` table `==` independent `modules.ext` recompute.
  - **oracle_selfcert:** the almost-positive-roots identity; the AR-duality formula
    `Ext¹_A(X,Y)=dim Hom(Y,τX)` across the AR quiver **A₂–A₅/D₄/zigzag-A₄** (0-mismatch
    proxy sweep); the FZ quiver-mutation certificate (`surfaces.flip`); the multi-3-cycle
    potential path; Gorenstein-dim ≤ 1 for `Jac(3-cycle)`; the D₄ n-regularity recovery;
    **the D₅ certified-or-refused gate** and **the budget-vs-infinite split** (both loud).
  - **HONEST SCOPE (binding, verification page) — THREE distinct boundaries, each with its
    honest refusal:** (i) **budget ≠ infinite** — a `status="budget"` stop on a known-finite
    Dynkin input is refused as OVER-BUDGET with the cluster number quoted, NEVER as infinite;
    infiniteness is claimed only off the P60 `dynkin_type` certificate (`Ã`/wild). (ii) **the
    D₅ / non-recoverable `status="error"`** — pre-P65 the cluster-tilting count of D₅ (and
    larger error types) is REFUSED, not fabricated; after P65's Task-0 root fix it is expected
    to certify natively. (iii) **the shifted-pair 2-CY** is asserted by BMRRT citation, not
    computed (module-window verdict). Plus: **QPA 1.37 and Macaulay2 have no
    cluster-category / cluster-tilting-object construction** — no cross-system oracle for the
    category itself; add a QPA honest-**skip** gate (`tests/qpa/test_cluster_category_qpa.py`)
    that a **live `NamesGVars()` sweep** makes FAIL if QPA ever ships such a verb (P35
    precedent). Record the OUT items (dg/Ginzburg, `D^b(Γ)`, general non-acyclic `C_{(Q,W)}`,
    DWZ potential mutation, rep-infinite `#indec`) with the theorem/engine that would be needed.
- [ ] **README** feature line: "Amiot–Keller cluster categories — the certified acyclic
  (Dynkin) slice: #indec, cluster-tilting objects (= support τ-tilting, AIR), the
  cluster-tilted End-algebra as a Jacobian algebra, and the 2-Calabi–Yau certificate."
- [ ] **Metaplan** tick `[x] P79`; add the dated R31 correction note (0805.1035, not
  2107.02646).
- [ ] **P80 ledger entry** (`docs/plans/DEEPER-ENGINES-BACKLOG.md`): the OUT items as the
  research-grade continuation (Ginzburg dg / `D^b(Γ)` / general `C_{(Q,W)}` / DWZ potential
  mutation), plus the P45-mutate-root-cause dependency on P65 Task 0.
- [ ] **Recount** the oracle-class table live post-implementation and update
  `docs/verification.md`'s audited counts + the suite badge (`tests/release/
  test_oracle_classes.py`); run the full acceptance gate.

---

## Acceptance (Plan-79 definition of done)

- [ ] `ClusterCategory` public via `import quiverlab`; `Algebra.cluster_category()` delegate
  (hereditary-only) live.
- [ ] Every Live-verified fact (§) is a passing test with its Plan-32 marker; every Scope
  gate (G1, G2, **G3a, G3b**, G4–G8) is a passing loud-refusal test — including the D₅
  certified-or-refused gate and the budget-vs-infinite split.
- [ ] The `cluster_category` kind served by **all three tiers** (spec twin parity byte-equal;
  webapp schema guard; report render), i18n ×4 parity, one golden, canonical-key stability
  confirmed; the estimator sizes on the cluster number/budget (E₆/A₇ route to the queued
  tier, not instant).
- [ ] **Latency budget honesty:** the deep tests carrying long BFS (A₅ ≈84 s, D₅ ≈284 s) are
  in the deep bucket, marked, with the measured time in a comment; nothing long is `-m fast`.
- [ ] Citations BibTeX-verified + registered; `tests/citations/test_bib_structure.py` green.
- [ ] Verification page updated (oracle map + honest-scope + audited recount); README +
  metaplan + P80 ledger updated.
- [ ] Gate: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast`
  and `-m deep` green on **both** kernel paths (`QUIVERLAB_NO_NUMBA=1` too);
  `tests/webapp/test_js_parses.py`, `test_runner_delegation.py`, `test_i18n.py`,
  `tests/release/` green; `mkdocs build --strict` exit 0.
- [ ] `git grep -n '<<<<<<<\|>>>>>>>'` clean after any merge; `; echo EXIT=$?` on gate chains.

---

## Methodology & assumptions

**Approach.** I ran the feasibility spike **first and live** in the venv (which tracks a
moving `dev`, ≈ Plan 74 at spec time) before freezing any scope: (1) computed the
module-category model numbers (`knit_ar_quiver`, `exchange_graph`), the cluster-tilted
flagship (`JacobianAlgebra` vs `NakayamaAlgebra`), and the 2-CY AR-duality anchor
(`ext_dims`/`tau`/`hom_space`) directly; (2) characterized the D₄ `status="error"` boundary
and found the merged P63 `_closed_by_n_regularity` recovery certifies it; (3) verified all
eight citations against arXiv/journal records, correcting the card's R31 arXiv id
(2107.02646 is the skew-gentle Amiot paper; the cluster-category paper is 0805.1035). The
scope was then frozen to exactly what computed. **The adversarial fix round (2026-08-08)
then corrected two over-claims by live computation (see below), and I re-froze the scope as a
runtime property rather than a type list.**

**Assumptions.** (a) The BMRRT fundamental-domain description `ind(mod kQ) ⊔ {P_v[1]}` and
the AIR cluster-tilting ↔ support-τ-tilting bijection are the correct models for the acyclic
case — standard, and cross-checked numerically (Catalan/cluster numbers, D₄=50). (b) The
type-A cluster-tilted potential is the sum of oriented 3-cycles — a known instance, verified
by the `kZ₃/J²` dim-6 match; a multi-3-cycle instance is spec'd (Task 3) to exercise the
all-cycles path; I did NOT assume it for D/E (G5 defers there). (c) `dev` is a moving target;
every number is a **runtime property re-pinned by the acceptance tests at branch-cut**, and
the gate logic is written to be invariant across the pre-/post-P65 regimes.

**What I deliberately did NOT check, and why it is safe.** (1) **CORRECTED from the first
draft:** the first draft claimed D₅+ "recover … robustly covers all Dynkin A/D/E" WITHOUT
running the recovery on D₅ — that was an over-claim. The fix-round live run shows D₅'s
exchange graph is `status="error"`, `_closed_by_n_regularity == False`, `complete == False`
(182 discovered but NOT n-regular pre-P65), so the cluster-tilting count is **honestly
refused** (G4), not fabricated. The scope is now a runtime property (Spike 2), the D₅ test
asserts certified-or-refused in both regimes, and D₅ is the *natural* non-recoverable case in
risk 2. (2) I did not build the Ginzburg dg algebra or attempt `Hom_C` directly — recorded
out-of-scope; the spike confirmed no dg engine ships. (3) I did not attempt a general
"mutation-equivalent to acyclic" recognizer for `from_potential` — kept narrow/honest (the
documented 3-cycle↔A₃ instance) with G7 refusing the rest. (4) The AR-formula ordinary=stable
coincidence is a **theorem-backed proxy** (ASS IV.2.13) confirmed by the fix-round's expanded
0-mismatch sweep over **A₂–A₅, D₄, and zigzag-A₄** (more robust than the first draft's single
`A₃` claim); the payload states the exact `\overline{Hom}` theorem and the module-window
scope. (5) I did not re-run P45's internal correctness (its own suite owns that); I only read
`exchange_graph`/`_closed_by_n_regularity` as public, and this plan does not modify
`mutation.py`. (6) I did not re-run engine probes during the fix round — `dev` was
mid-merge; the two refutations (D₅ non-recovery; budget≠infinite) are the critic's live runs,
adjudicated valid, and are re-pinned by the acceptance tests.

**Why the result is correct.** Every headline number was recomputed live and matches an
independent literature value (Catalan/cluster numbers, almost-positive-root counts, the
Buan–Marsh D₄=50, the `kZ₃/J²` dim 6) AND an independent engine (`exchange_graph` vs Catalan;
`Jac` vs `NakayamaAlgebra`; `Ext¹_C` vs `modules.ext`). The honest-scope boundaries are not
asserted but **demonstrated** by the spike — the D₄ error+recovery, the **D₅
error+non-recovery refusal**, the **budget≠infinite** distinction, no dg engine, and the
Jacobian gl.dim ≥ 32 infinite.

---

## Open design risks (what a critic will likely attack)

1. **"The 2-CY certificate is a tautology (`a+b=b+a`)."** Addressed: the payload separates
   the *by-construction* symmetry from the *genuine* AR-duality content (`Ext¹_A(X,Y)=
   dim Hom(Y,τX)`, 0 mismatches over A₂–A₅/D₄/zigzag-A₄, ASS IV.2.13), and is a **scoped**
   verdict (module window; shifted pairs by citation), not a bald `True`.
2. **"D₅ / non-A cluster-tilting counts rest on a *recovered* error graph — and D₅ does NOT
   recover."** The fix round confirmed exactly this: D₅ is `status="error"`,
   `_closed_by_n_regularity == False` ⟹ **refused (G4)**, never fabricated. Certification is a
   per-input runtime property; `certified`/`status`/`count(None)` are surfaced. D₄ recovers
   (n-regular); D₅ does not. When P65's root fix lands, `mutate` stops raising and such types
   certify natively as `status="complete"` — the gate logic is unchanged, and the D₅ test
   accepts both regimes (certified-or-refused-honestly).
3. **"You claim `End_C(T)=Jac` but never compute `Hom_C`."** Correct and by design — the
   equality is a **cited theorem** (BMR/Amiot); we **verify the target Jacobian's
   presentation** (dim/Cartan/self-injective), which is the falsifiable content. The
   docstring and report say "cited, not computed".
4. **"Non-acyclic scope is thin."** Intentional (card §1.4 "scope fixed, not open-ended"):
   `from_potential` gives the algebra-level Amiot+Keller–Reiten certificate for any
   Jacobi-finite `(Q,W)`, and the category combinatorics only where a theorem reduces to the
   hereditary model. Everything else is a loud refusal + a P80 ledger entry.
5. **"GUI adds an expensive kind that the estimator under-sizes."** Fixed in the fix round:
   the exchange-graph BFS is NOT algebra-dim-bound (E₆ dim 36 → 833 vertices, ≈8 min), so the
   `cluster_category` estimator branch sizes on the **cluster number / budget** (P60
   `dynkin_type` closed form) and routes E₆/A₇ off the instant tier; rep-infinite inputs
   refuse before enumerating; over-budget finite inputs are refused as over-budget (not
   infinite).

---

## Change log
- 2026-08-08: initial plan. Feasibility spike run live (dev venv); scope frozen to the
  certified acyclic/Dynkin module-category slice + the Jacobi-finite algebra-level
  certificate; GUI IN (`cluster_category` kind). All 8 citations BibTeX-verified; R31 arXiv
  id corrected (0805.1035). Dg/Ginzburg/`D^b(Γ)`/general non-acyclic/DWZ-potential-mutation
  recorded OUT with a P80 ledger entry.
- 2026-08-08 (adversarial fix round, all rulings applied, read-only doc edits — venv
  mid-merge): **MAJOR 1** — "D/E recover; all Dynkin A/D/E" REFUTED (D₅: 182 discovered,
  `status=error`, `_closed_by_n_regularity=False`, `complete=False` ⟹ G4 refuses). Reframed
  the frozen scope as a **runtime property** (Spike 2), added the D₅ both-regimes test, fixed
  risk 2 to name D₅ as the natural non-recoverable case, and stated the P65-Task-0 concurrent
  merge + post-fix native-certification expectation. **MAJOR 2** — `status="budget"` ≠
  rep-infinite REFUTED (A₇/D₆/E₆/E₇/E₈ cluster numbers exceed the default budget). Split G3
  into G3a (over-budget on a known-finite Dynkin — quote the cluster number, never
  infinite) and G3b (genuine infiniteness off the P60 `dynkin_type` certificate); added
  Spike 2b, the budget-vs-infinite test, the estimator branch (size on cluster number/budget,
  not algebra dim), and the block-shape provenance. **MODERATE** — `is_2_calabi_yau` now a
  scoped `{verdict, scope:"module_window", shifted_by_citation}` payload, not a bald `True`.
  **MINORS** — AR formula relabelled with **ASS Thm IV.2.13** + the ordinary=stable
  theorem-backed proxy (0-mismatch A₂–A₅/D₄/zigzag-A₄ sweep recorded as evidence); a
  **multi-3-cycle A₄/A₅** instance spec'd (Task 3) to exercise `_quiver_from_exchange_matrix`;
  honest latencies recorded (A₅ ≈84 s, D₅ ≈284 s); the tree pointer made tip-agnostic (no
  commit hash; `dev` a moving target). Methodology updated with the D₅-contradiction
  correction and the fix-round provenance.
