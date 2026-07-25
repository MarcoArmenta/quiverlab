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

**Status:** classical audit — LANDED; software systems — pending; modern
practice — pending.

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
