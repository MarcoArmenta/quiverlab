# Metagoal-coverage deep research (2026-07-25)

**Provenance:** Marco's metagoals, 2026-07-25 — (1) **no code required** for
users; (2) **any computation done in representation theory**, so every
representation theorist can profit from the library and the GUI. Three parallel
deep-research agents mapped the gap between the current surface and metagoal
(2): a classical-textbook completeness audit (ASS vols 1–3, ARS, Ringel,
Benson, Schiffler), a software-systems inventory (QPA/QPA2, Meataxe/Magma,
SageMath, String Applet/SBStrips, historical tools), and a modern
research-practice survey (τ-tilting industry, torsion lattices, silting,
higher AR, geometric models, …). The distilled, prioritized coverage program
lives in `ROADMAP.md`; this file preserves the full evidence.

**Status:** all three LANDED (classical audit; modern practice; software systems).

---

## Report 1 — Classical-curriculum completeness audit (ASS 1–3, ARS, Ringel, Benson, Schiffler)

Every "has" verified against the code (files listed by the agent:
`core/algebra.py`, `modules/*`, `invariants/*`, `combinat/quiver.py`,
`families/*`). Public surface confirmed on `Algebra` and `Module` as of
2026-07-25 (incl. the Plan-27 `ext_algebra`).

### Quiver & algebra basics (ASS I.2–II.3)
| Computation | Status | Notes | Priority |
|---|---|---|---|
| Build kQ/I, admissibility, certified finite-dim, opposite, duality D | HAS | build-time gates + Plan 23/24 | — |
| **Gabriel quiver of an abstract algebra** (recover Q, I from structure constants; Ext-quiver; basic-ization) | **LACKS** | `from_structure_constants` carries no quiver; path-basis invariants refuse loudly | **P1** |
| Connectedness predicate / is_basic / is_sincere / is_faithful | LACKS | one-liners over existing data | P2 |

### Euler / Tits form & roots (ASS VII; Schiffler 2–3; Benson 4.5) — the largest CHEAP gap
| Computation | Status | Priority |
|---|---|---|
| Cartan/Coxeter matrix + polynomial + spectral | HAS | — |
| **Euler bilinear form ⟨d,e⟩ = d^T C^{−T} e; Tits form q(d)** | **LACKS** (2 lines over C^{-1}) | **P1** |
| **Positivity ⇒ Dynkin/Euclidean classification; type DETECTION from a quiver** | **LACKS** (constructor exists, no recognizer) | **P1** |
| Roots (real/imaginary), reflections/Weyl orbits, radical of q, defect | LACKS | P2–P3 |

### Grothendieck / composition / Loewy (ASS I.1, III.3)
| Computation | Status | Priority |
|---|---|---|
| dim vector, rad/top/soc (one step), algebra Loewy length | HAS | — |
| **Radical/socle SERIES of a module + Loewy-diagram display** | **LACKS** (iterate radical(); the picture users love) | **P1**/P2 |
| K0 object, composition-series listing | PARTIAL (dim vector = multiplicities for basic) | P2–P3 |

### Projectives/injectives/Nakayama functor (ASS I.5, III.2)
S/P/I both sides, projective cover, min proj+inj resolutions, pd/id — HAS.
Injective envelope (named), Nakayama functor ν as `nu(M)`, stable Hom — PARTIAL/LACKS (P2).

### AR theory (ARS V–VII; ASS IV) — highest-visibility PARTIAL
| Computation | Status | Priority |
|---|---|---|
| τ, τ⁻, Tr, D (QPA-crosschecked) | HAS | — |
| **Middle term E of 0→τM→E→M→0 (almost-split sequence construction)** | **LACKS** | **P1** |
| Irreducible maps / rad(M,N)/rad² | LACKS | P2 |
| **AR-quiver knitting** (rep-directed) | **LACKS** | **P1** |
| AR components; QPA-wired AR quiver as a feature | LACKS/PARTIAL | P2–P3 |

### Tilting (ASS VI; Ringel)
| Computation | Status | Priority |
|---|---|---|
| is_tilting(T) | PARTIAL (needs summand count ⇒ decomposition) | P2 |
| **End(M) AS a quiver algebra (tilted algebras)** | **LACKS** (only end_dim; ext_algebra is Yoneda-of-simples, not End(M)) | **P1** |
| Torsion pair from T, APR tilts/reflection functors, Brenner–Butler | LACKS | P3 |

### Constructors from module data
Trivial extension, tensor — HAS. **One-point (co)extension A[M]** — LACKS (P2). Repetitive Â — LACKS (P3).

### Representation type (ASS VII–VIII; Ringel)
Rep-finite/tame/wild decision, Bongartz criterion, rep-directed, knitting classification, components — ALL LACK (P2–P3; hereditary case unlocked by the Euler-form layer).

### Self-injective/Frobenius — HAS (is_weakly_symmetric one-liner missing, P3).

### Benson-side
Complexity HAS; periodicity PARTIAL; **kG → quiver bridge** LACKS (P3); support varieties LACKS (Tier-2 backlog); group cohomology ring — out-of-scope-ish.

### Simson-side
Species, coalgebras/comodules — **out of scope** (state explicitly). Koszulity/quadratic dual/Yoneda — HAS (Plan 27, a genuine strength).

### Decision procedures & the GLUE (limits "any computation" most)
| Computation | Status | Priority |
|---|---|---|
| **Is M indecomposable? (End(M) local)** | **LACKS** | **P1** |
| **Krull–Schmidt: decompose M into indecomposables** | **LACKS** (no inverse to _direct_sum; idempotent lifting) | **P1** |
| **Module maps as first-class objects; kernel/cokernel/image as Modules** | **LACKS public** (radtopsoc.submodule/quotient internals exist) | **P1** |
| Pushout/pullback, SES utilities, is-split | LACKS | P2 |
| Module iso test, dim Hom/End/Ext | HAS | — |

### Consolidated ranking (classical)
1. Krull–Schmidt decomposition + indecomposability (unblocks tilting, AR readout, multiplicities, knitting).
2. First-class module maps with kernel/cokernel/image + SES/is-split ("τ of the cokernel of this map" — largely a surfacing job).
3. Euler/Tits forms + definiteness + Dynkin/Euclidean type detection (cheapest high-payoff).
4. Almost-split sequence middle term; then knitting.
5. Gabriel-quiver recovery / basic-ization of abstract algebras.
6. End(M) as quiver algebra + is_tilting + one-point extension.
7. Loewy series + display; K0; ν named.
8. Representation-type decisions; tame-hereditary structure (defect/tubes).
9. kG bridge + support varieties.
10. Species/coalgebras: out of scope, stated.

## Report 2 — Modern research-practice survey (the last ~15 years)

**The single most important structural finding:** survey areas 1–4 (τ-tilting,
torsion classes/bricks/wide subcategories, stability/walls/c-vectors, 2-term
silting) are **not four features — they are one engine**. AIR + DIJ:
{support τ-tilting modules} ≅ {2-term silting complexes} ≅ {functorially finite
torsion classes} ≅ {chambers of the g-vector fan}, bricks = wall labels,
c-vectors = dual basis. quiverlab already has every primitive: minimal
projective presentations (g-vectors = [P₀]−[P₁] in K₀(proj)), `hom`, `τ` — the
τ-rigid test IS `Hom(M, τM) = 0`, computable today by composing two existing
calls. The only genuinely new code is the mutation / Bongartz-completion
enumerator. Honest decidability: the mutation BFS terminates iff the algebra is
τ-tilting finite (semi-algorithm; cap → "complete: N" or "cap hit,
τ-tilting-infinite witness (uncertified)").

| # | Area | Output | Algorithm (who) | Priority |
|---|------|--------|-----------------|----------|
| 1 | τ-tilting | support τ-tilting modules, g-vectors, is_τ_rigid, exchange graph, is_τ_tilting_finite | mutation BFS (Adachi–Iyama–Reiten 2014; Demonet–Iyama–Jasso 2019) | **P1** |
| 2 | Torsion/bricks/wide/semibricks | torsion-class lattice + brick-labelled Hasse, semibricks↔wide (Asai) | rides #1 (DIRRT 2023) | **P1** |
| 3 | Stability/walls/c-vectors | wall-and-chamber fan, θ-stable modules; **n=2,3 the drawable killer GUI demo** | King 1994; Bruestle–Smith–Treffinger 2019 | **P1** |
| 4 | 2-term silting (+exceptional/cluster) | =#1; hereditary exceptional sequences, clusters | AIR; Aihara–Iyama 2012; BMRRT 2006 | **P1**/P2 |
| 5 | AR quiver knitting; higher AR | AR quiver + sequences (rep-finite) | knitting (ASS); QPA covers today → native is P2 | P2/P3 |
| 6 | Representation type | finite/tame/wild; indec. dim vectors = roots | Gabriel 1972; **Kac 1980/82 root combinatorics = cheap exact P1**; Bongartz 1984; BGRS 1985; Drozd (no general decider — honest) | **P1**(hereditary/Kac)/P2/P3 |
| 7 | Gentle/string/special biserial; Brauer graph | ALL indecomposables as string/band modules, AR quiver, AAG derived invariant, surface model | Butler–Ringel 1987 (complete, combinatorial); AAG 2008; Opper–Plamondon–Schroll 2018 | **P1** (gentle)/P2 |
| 8 | findim/domdim/Igusa–Todorov/Gorenstein/delooping | homological-dimensions card; GProj/singularity (gentle combinatorial, Kalck 2015) | Igusa–Todorov 2005; Gélinas 2022 | **P1** (findim/IT)/P2 |
| 9 | Derived/stable equivalence | derived-fingerprint panel (HH/HC/Cartan/center ALREADY implemented); Brauer-graph + gentle classifications | Rickard 1989; Antipov–Zvonareva; Opper–Zvonareva | P2/P3 |
| 10 | Quiver geometry | Kac canonical decomposition (Schofield 1992 / Derksen–Weyman 2002), orbit dims + Voigt, degeneration order (Zwara 2000), Reineke Betti (P3) | exact integer LA + generic rank | P2/P3 |
| 11 | Species/valued quivers; group-algebra blocks; Hopf examples | B/C/F/G types (Dlab–Ringel 1976, needs a modulation layer over GF(p^n)); consume block basic algebras; curated Taft/quantum library | | P2/P3 |
| 12 | Hall numbers (exact GF(q) point counts); DT/scattering | | | P3 / out-of-scope (state) |

**One-line recommendation:** build the τ-tilting mutation engine first (four P1
areas at once, killer n=2/3 wall-picture GUI), then gentle/string combinatorics
(pairs with Nakayama), the findim/Igusa–Todorov family (pairs with
gldim/ext/complexity), and Kac/Dynkin root combinatorics as the cheap
high-visibility side-quest.

## Report 3 — Software-systems gap inventory (QPA/QPA2, MeatAxe/Magma/CREP, SageMath, SBStrips/String Applet/FD Applet)

Verified against quiverlab's public surface, the QPA `.gd` declaration files in
the local `[qpa]` install, and fetched docs. **Structural driver of all
priorities: quiverlab already owns the expensive primitives (τ/D/Tr/A^op,
Ext/Hom dims, minimal proj+inj resolutions, exact module iso, End dims) — the
missing top-tier features are COMPOSITION LAYERS over those primitives, not new
kernels.**

### QPA has / quiverlab lacks (from the .gd files)
- **Module category as a category (biggest cluster):** `HomOverAlgebra`/
  `EndOverAlgebra` return actual MAP BASES (we return dims);
  Kernel/Image/CoKernel of maps; split mono/epi tests; right-minimal versions;
  `DecomposeModule`(+multiplicities) via idempotent lifting;
  `IsIndecomposableModule`; is_projective/injective/simple; direct sums +
  `IsDirectSummand`/`CommonDirectSummand`; submodule lattices;
  `EndOfModuleAsQuiverAlgebra`; Trace/Reject; minimal generating set;
  annihilator; full Radical/Socle SERIES; ADR algebra.
- **AR theory:** `AlmostSplitSequence`, `PredecessorsOfModule` (knitting),
  `IrreducibleMorphisms*`.
- **Tilting/approximation:** `IsTiltingModule`/cotilting, complements +
  left/right mutation, minimal left/right approximations, `IsTauRigidModule`,
  `IsRigidModule`, ProjectiveCover/InjectiveEnvelope AS MAPS, Iyama generator,
  faithful dimension.
- **Dimensions:** dominant dimension, Gorenstein dimension/`IsGorenstein`,
  Ω/τ-periodicity tests.
- **Unit forms (CHEAP for us — Cartan exists):** Tits/Euler forms, positive
  roots, weak positivity/nonnegativity, reflections.
- **Complexes/derived:** user-facing complexes, chain maps, cones, homology,
  truncations, `YonedaProduct`, comparison liftings.
- **Functors:** Nakayama functor ν on modules, `StarOfModule`, module tensor
  products, restriction along algebra maps.
- **Constructors/recognizers:** blocks; `IsHereditary/Gorenstein/Semisimple/
  RadicalSquareZero/SpecialBiserial/String/Gentle/Nakayama/Canonical`;
  Canonical/Kronecker/BrauerConfiguration/Poset algebras; quiver predicates
  (connected/tree/Dynkin/components/subquiver); degeneration order + numeric AR
  for finite type.
- quiverlab EXCEEDS QPA on: all Hochschild/cyclic/TT-calculus, CS resolution,
  Koszulity/quadratic dual, deep engines, exact CC/GF(p^n), certification.

### String/gentle world (SBStrips, String Applet, FD Applet)
Complete SB/string layer absent from quiverlab: strings/bands, string-module
syzygies (SBStrips), SB AR quiver, delooping level; String Applet adds τ-tilting,
g-matrices, complete gentle derived classification (AAG + Nakayama permutation),
surface models. **FD Applet + String Applet are the UX template for the no-code
metagoal.** AAG invariant + surface/winding derived invariant have NO maintained
standalone package — genuine white space quiverlab could fill (ties to our HH).

### SageMath
Module theory essentially hereditary-only (no admissible kQ/I, no Ext/HH) — our
core is absent there. Sage has what we lack: **cluster mutation**
(ClusterQuiver/ClusterSeed: mutation type detection, g/c/d-vectors,
F-polynomials, g-vector fan) and morphism objects (QuiverRepHom with
kernel/cokernel — hereditary-only, silently wrong for kQ/I, an edge we already
handle right).

### MeatAxe/Magma/CREP — the decomposition flagship
GAP `MTX.Indecomposition` (finite fields, Las Vegas — randomness affects
runtime never correctness), Magma `IndecomposableSummands` (incl. number
fields, via maximal orders in End(M)), CREP (historical: Tits forms, AR
components). Math for our framework: decomposition = primitive orthogonal
idempotents of End(M) (Azumaya); GF(p)/GF(p^n) = idempotent lifting through the
radical; CC = clean once char-polys factor; Q/number fields decidable but
bottlenecked — matches our "exact, refuse loudly beyond budget" stance.

### τ-tilting
No maintained calculator anywhere general (FD Applet = rep-finite-ish UX;
feisele's `tautilting.g` = algorithmic blueprint; **QPA has no τ-tilting
mutation command**). Building it would make quiverlab MORE general than every
existing tool. Verification: QPA oracles almost every other new feature;
τ-tilting oracles against theory + FD Applet + the blueprint.

### Strategic summary
1. One lever unlocks the top tier: morphisms + End(M) + decomposition.
2. Two cheap immediate wins: unit forms; recognizer batch.
3. τ-tilting = highest-leverage modern feature (white space).
4. QPA extends the crosscheck harness to nearly all of it.
5. Genuine white space: AAG/surface derived invariant; wall-and-chamber +
   maximal green sequences; general-field no-code module decomposer.

Out of scope (with reasons): species/valued quivers (refoundation);
DWZ potential mutation/cluster categories (Jacobian-algebra CONSTRUCTOR is in
scope as a family); DT/BPS/moduli counting; group-MeatAxe at scale;
IsBasic/IsElementary (vacuous here).
