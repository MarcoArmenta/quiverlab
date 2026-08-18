# How quiverlab is verified

Every shipped feature of quiverlab is unit tested. This page says **how** — with
the highest rigour we can bring to it — and it is honest about the edges: where a
check is a cross-engine agreement, where it is a published number, where a live
external oracle can reach, and where it cannot.

The suite is **5338 tests** (live-recounted 2026-08-17 on the Plan-77 branch — R36 the
generalized Koszulity ladder (Berger N-Koszul, Cassidy–Shelton K₂, Brenner–Butler–King
`(p,q)`-almost-Koszul, the internal generation degrees of `Ext•(k,k)`) — cut from dev's
P76 tip; the previous recount was the Plan-75 branch (R9 incidence-vs-nerve `HH^*` + R10
the GHMS comultiplicative Koszul resolution) cut from dev's P73 tip, itself recounted at the Plan-73 top-up merge onto dev's P78 tip; collected
with the `[dev,fast,docs,web,qpa,hpc]` extras,
2026-08-05, after Plans 21–33, the Plan-48 marked-surface subsystem (marked surfaces →
ideal triangulations → gentle Jacobian algebras; +70 tests), the Wave-1 v0.2.0 trio — Plan-36 Macaulay2 fifth
oracle class, Plan-37 C1 categorical glue, Plan-38 C2 forms/type/positive-roots/
recognizer batteries + Koszulity (`ext_algebra`) exposure — the Plan-32 oracle-class markers + audit gate,
Marco's report-completeness pass, Marco's Cayley product-table render wave, the
Plan-35 Hochschild product surface, its
explicit-representatives capture, the Plan-35 UNIT-2/wave-3a rendering, the
Plan-35 wave-3b cyclic-homology explicit representatives, the Plan-35 wave-3c
Yoneda exact sequences + classical dictionary, and the Plan-35 wave-3d plain-HH
explicit representatives + element-wise dictionary read-offs — each product / Ext /
Tor / HC / HH class now ships its (co)cycle as a labeled term-sum + a coordinate
vector self-certified against the differential, every Ext class is CONSTRUCTED as
its explicit exact sequence, and the report/GUI lay it all out per degree with each
space's classical interpretation stated — HH⁰'s centre, HH¹'s derivations, HH²'s
deformation cochain, HH₀'s commutator residues read straight off the reps; and the
Plan-36 Macaulay2 oracle bridge — a fifth, external oracle class, `-m m2`). It
is not a pile of smoke tests: the mathematics is pinned by **two classes of
oracle**, and most numbers are checked by more than one. Every test is
**classifiable** into exactly this scheme — one of five oracle classes (literature,
cross-engine, self-certifying, live QPA, or live Macaulay2) or the
contract/infrastructure remainder — and since Plan 32 each oracle class is also a **standalone one-liner** a
reviewer can run (`pytest -m oracle_literature`, and so on), with the counts audited
against live collection (see [Oracle classes as runnable markers](#oracle-classes-as-runnable-markers)).

## The two oracle classes

**1. Theory and literature, on constructed examples.** We build many algebras the
literature (or a theorem we know) has already resolved, and assert quiverlab
reproduces the published value exactly. When no single published vector is at
hand, we cross-check against an *independent* computation path in the library and
say so inline. Nothing is a float; every equality is exact.

**2. Cross-engine and external (QPA/GAP) agreement.** quiverlab ships several
independent engines for the same invariants (the normalized bar complex; the
minimal `A^e`, Bardzell, and Chouhy–Solotar resolutions). Where two of them
overlap they must agree degreewise, over several characteristics. And wherever the
GAP package **QPA** implements a feature, we recompute with it and demand
equality. QPA does not implement everything quiverlab does — this page names
exactly where it is used and where a theory oracle stands in for it instead.

These are complementary. A regression that corrupted *two* of the library's own
engines identically would still be caught by a literature pin or by QPA; a
literature pin that was mis-transcribed is caught by the live cross-engine
agreement. The tests are written so that agreement is never a tautology.

---

## Class 1 — theory and literature oracles

### The normalized bar complex is the base oracle

Hochschild `HH^n` / `HH_n` from the exponential normalized bar complex
(`hochschild/bar.py`) is the ground truth every deeper engine is measured against
(Hochschild, *Ann. of Math.* 46 (1945), 58–67; registry key `bar`).
It is simple, slow, and structurally different from the resolutions, so degreewise
agreement with it is a genuine check. Its use is bounded only by the exponential
blow-up `dim C_n = a·(a-1)^n` — the whole point of the deeper engines is to reach
past that window, so on the overlap range the deeper engines are pinned to the bar
oracle, and past it they are pinned to each other and to closed forms.

### Multi-prime cross-checks

Engine-level batteries run over the prime set `{32003, 2, 3, 5}`: a large prime as
a characteristic-0 proxy, plus the small primes that probe characteristic
pathology (for `k[x]/(x^a)`, a prime dividing `a` collapses the norm map; for the
quantum complete intersection the small primes reshape the homology). A bug that
only shows up in one characteristic cannot hide.

### Named literature pins

`tests/resolutions_cs/test_battery_literature.py` and
`tests/engine/test_qpa_reference_validation.py` pin values that exist **outside**
the library, each with its provenance inline. Every source below has a full entry
in the [References](#references) section at the bottom of this page; where a work
carries a citations-registry key its bibliographic entry is the packaged
`references.bib`, cited by that key.

- **`k[x]/(x^n)`, classical.** In characteristic 0 (or whenever `n` is invertible)
  `dim HH_0 = n` and `dim HH_i = n-1` for `i ≥ 1`. Sources named in the tests:
  Loday, *Cyclic Homology*, and the "BACH" truncated-polynomial computations
  (CS battery); Happel (1989), registry key `happel_question`, for the same
  values (QPA-reference file) — and both are cross-confirmed live by the bar
  oracle. The char-5 pathology on `k[x]/(x^5)` (`p | n`) collapses both
  differentials, giving `HH_i = A` in every degree — pinned as `[5, 5, 5, 5, 5, 5, 5]`.
- **Quantum complete intersection** `k⟨x,y⟩/(x², y², yx − 2xy)`. The family (finite
  Hochschild cohomology with infinite global dimension) is
  Buchweitz–Green–Madsen–Solberg, *Math. Res. Lett.* 12 (2005), 805–816 (registry
  key `quantum_ci`); the explicit `HH_•` / `HH^•` values are Bergh–Erdmann,
  *Algebra & Number Theory* 2 (2008), 501–522 (key `qci_hh_oracle`). Homology
  persists, `HH_• = [3, 2, 2, …]`, while cohomology dies from degree 3,
  `HH^• = [2, 2, 1, 0, 0, …]` in characteristic 0 — the homology/cohomology
  asymmetry the deeper engines exist to detect.
- **Happel's theorem.** Happel, *Lecture Notes in Math.* 1404 (1989), 108–126
  (registry key `happel_question`): a hereditary algebra has `HH^i = 0` for
  `i ≥ 2`. Linear `A_3` is pinned to `HH^• = [1, 0, 0, 0, 0]`. (The tests state the
  theorem in words; no theorem number is asserted here, because none is recorded in
  the repository.)
- **Commutative complete intersection** `k[x,y]/(x², y²)`: `HH_• = 4, 4, 5, 6, 7,
  8, 9, …` — the value is attributed to Buchweitz–Green–Madsen–Solberg (2005) in
  the QPA-reference file, via the Künneth square of `k[x]/(x²)` (Künneth formula:
  Cartan–Eilenberg, *Homological Algebra*, 1956, registry key `tensor_product`).
  Being symmetric it also satisfies `HH^n = HH_n` — an internal cross-check QPA
  satisfies too.
- **Cyclic Nakayama** self-injective algebras `kZ_n/rad²`: low-degree dims frozen;
  the tests cite these as the standard QPA Nakayama-algebra examples (`references.bib`
  key `qpa`), the family itself per Assem–Simson–Skowroński (2006), key `nakayama`.
- **Gentle algebra** `kQ/(ab, ba)` on the 2-cycle: self-injective, so `HH^•` is
  nonzero in every degree (`[1, 1, 1, 1]`) — which makes the CS↔bar agreement a
  discriminating check rather than a run of zeros. Oracle: bar cross-check; the
  gentle-algebra notion is attributed in the test to Assem–Skowroński, *Algebra i
  Analiz* (1987).
- **Auslander–Reiten theory on worked examples** (Plans 23/24). For hereditary
  algebras the translate satisfies the Coxeter-transformation law
  `dim τM = Φ⁻ᵀ · dim M` (in quiverlab's `e_iAe_j` Cartan convention, calibrated
  against QPA's `DTr`) — an independent path-counting check of the op+D+Tr
  construction; plus the explicit kA₂/kA₃ AR tables, the Nakayama τ-orbit on
  `kZ_3/rad²`, `τ(projective) = 0` / `τ⁻(injective) = 0`, self-injective ⟺
  `inj.dim ∈ {0, ∞}` (`k[x]/(x²)`), `max_v inj.dim S_v = gl.dim` (commutative
  square = 2, hereditary kAₙ = 1), and the honest kA₂ side asymmetry (right
  `P(1) = e_1A` vs left `P(1) = Ae_1` have different dimension vectors). Source:
  Assem–Simson–Skowroński (2006), registry keys `assem_book`/`nakayama`, cited at
  chapter granularity (`tests/modules/test_duality_tau.py`, `test_injective.py`,
  `test_left_modules.py`).
- **The Plan-29 literature batteries** (2026-07-25; sources per test docstring +
  the registry keys added with them). Coxeter/spectral: the six exact Nakayama
  Coxeter polynomials + the `χ(r+7,r)` family and Coxeter numbers up to 126
  (Lenzing–Meltzer–Ruan, `lenzing_meltzer_ruan`), the Dynkin/affine/canonical
  tables and the star formula (Lenzing–de la Peña, `lenzing_delapena_spectral`),
  Lehmer's polynomial as the `[2,3,7]` Coxeter polynomial with the Mahler/spectral
  ordering (de la Peña, `delapena_mahler`), and the `χ(−1)`-perfect-square sweep.
  Identity oracles: Happel's trace formula `tr Φ = −Σ(−1)^i dim HH^i`
  (`happel_trace`, sign pinned on A₃), derived-invariance of HH^•/HH_•/HC_•
  across quiver orientations (`keller_cyclic_invariance`/`rickard_derived`),
  `HH¹(T(A)) ≠ 0` (`cmrs_split`), incidence ≅ nerve cohomology
  (`cibils_incidence`/`redondo_incidence`), acyclic HH_{≥1} = 0
  (`cibils_acyclic`), and the truncated finiteness boolean (`xhj_truncated`).
  Value batteries: Bergh–Erdmann quantum-CI cohomology `[2,2,1,0,…]` for
  general `(a,b)` over char 0 (`qci_hh_oracle`), the Redondo–Román
  triangular-string family with its degree-(2m+1) revival — anchored by the
  independent minimal-A^e engine (`redondo_roman_2014`), the radical-square-zero
  char-2 doubling (`cibils_radsq`), Taillefer's Taft-algebra cyclic homology
  (`taillefer_taft`), and canonical-algebra `HH² = t−3` (`schremmer_wpl`,
  attributing Happel) cross-linked to the trace identity. Symmetry: the
  Brauer-star battery (`skowronski_yamagata`, `n | (L−1)`) that exposed and now
  guards the `is_symmetric` fix. Module Tor: the duality anchor
  `dim Tor_n(M,N) = dim Ext^n(M, DN)` on every case, the resolve-either-side
  balance, and vanishing laws (`tensor_product`/`assem_book`).
- **The Plan-58 certified Coxeter spectral analysis** (2026-08-07; record R20).
  On top of the shipped exact Coxeter/spectral primitives, `A.coxeter_spectral()`
  is pinned to the record's oracle set: kA₂ = Φ₃ (finite Coxeter order 3), kA₄ = Φ₅
  (order 5), kA₅ = Φ₂Φ₃Φ₆ (three cyclotomic factors); the 3-Kronecker χ = t²−7t+1
  with ρ = M = (7+3√5)/2 as a **certified algebraic number** (minimal polynomial
  [1,−7,1], rational isolating interval (6,7), one root outside the unit circle);
  the m-Kronecker ladder x²−(m²−2)x+1; T_{2,3,7} = E₁₀ realizing **Lehmer's**
  degree-10 polynomial (interval (1,2), the smallest known Salem number); and the
  ADE/affine cyclotomic sweep (finite order = Coxeter number on Dynkin, `None`/Jordan
  on the affine defect — the periodic ⊊ cyclotomic-type separation). Sources:
  de la Peña, Adv. Math. 2014 (`dlPena2014mahler`) and 2013 arXiv:1310.1557
  (`dlPena2013cyclotomic`), and de la Peña–Takane, Arch. Math. 55 (1990) 120–134
  (`dlPenaTakane1990spectral`). The Lehmer dichotomy is reported as
  **class-conditional documentation only** — quiverlab emits no verdict on the open
  general Lehmer problem, and μ₀ is referenced by its minimal polynomial, never a
  float.
- **The Plan-31 trivial-extension presentation** (2026-07-26). The certified
  `kQ_T/I_T` build is pinned to classical special cases: `T(kA_n) ≅ kZ_n/J^{n+1}`
  (the symmetric Brauer star, verified n = 2, 3, 4 — dims 6/12/20, Loewy length
  n+1, all four symmetry booleans), and `T(k[x]/(x^a)) = k⟨x,y⟩/(x^a, y², xy−yx)`
  — the plain commutator in **every** characteristic (`D(A)` is an honest
  bimodule, no Koszul sign), whose a = 2 case is `k[x,y]/(x², y²)`, reproducing
  the existing `HH_• = [4,4,5,6]` pin above. The Cartan identity
  `C_T = C_A + C_Aᵀ` (repo convention: entry `dim e_i A e_j`) holds on every case
  including the zoo's `line_abc_cde`; the presented and ⋉ structure-constant
  builds agree degreewise on bar-HH (iso-invariance), and CS ≡ bar now serves the
  presented `T(kA₂)` (CS refused the old structure-constant build). Sources:
  Assem–Simson–Skowroński (2006, `assem_book`) and the symmetric-Nakayama
  criterion of Skowroński–Yamagata (`skowronski_yamagata`); the `HH¹(T(A)) ≠ 0`
  summand (`cmrs_split`) is the Plan-29 companion. The presentation itself cites
  no closed-form theorem — it is per-instance dimension-certified and QPA-oracled
  (see [Honest scope](#honest-scope)). `tests/families/test_trivial_extension_presented.py`.
- **The Plan-35 product surface** (2026-08-01). The public cup / cap / Gerstenhaber
  bracket tables and the induced Connes `B` are pinned on the dual numbers
  `k[x]/(x²)` and the commutative complete intersection. Over a char-0-shaped prime
  (`GF(32003)`, where `2` is a unit) the cup ring of `k[x]/(x²)` is the classical
  `HH^• = [2, 1, 1, 1, …]` — `HH^0 = Z(A) = A` of dimension `2`, then `k` in every
  positive degree — with the even generator composing to a nonzero even class and
  the **odd square vanishing** by graded commutativity; over `GF(2)` that odd square
  **survives**, the classical characteristic-2 phenomenon. *Correction pinned
  (the engine wins, per the CRS-2004 precedent):* the degree-0 dimension is `2`, not
  `1` — the implementation-plan brief stated `1`, but for a commutative algebra
  `HH^0 = Z(A) = A` and `HH_0 = A/[A,A] = A` both have dimension `dim A = 2`; the
  frozen value is the verified `2`, and the "`dim HH^n = 1`" statement holds only for
  `n ≥ 1`. The `QuantumCI(q=1)` cup-table dimensions reproduce the
  Buchweitz–Green–Madsen–Solberg commutative-CI vector `[4, 8, 12, …]` (`quantum_ci`;
  the Künneth square of `k[x]/(x²)`). Connes `B` on the dual numbers alternates
  iso/zero along the SBI pattern (`rank B_0 = 1`). Sources: `bar`, `cup`,
  `gerstenhaber`, `bracket`, `cyclic`, `quantum_ci`, `tensor_product`
  (`tests/hochschild/test_products_literature.py`, `test_connes_b.py`).

### The read-only bank as a byte-level oracle

`tests/resolutions_cs/test_battery_bank_oracle.py` pins the Plan-04 CS resolution
against the original hanlab bank's *hand-derived closed-form* CS differentials — a
wholly separate implementation of the Chouhy–Solotar formulae (Chouhy–Solotar,
*J. Algebra* 432 (2015), 22–61, arXiv:1406.2300; registry key `chouhy_solotar`)
for `k[x]/(x^a)` and the quantum CI. Two tiers: **HH-dimension equality** family
by family and prime by prime (rank-based, invariant under the correction's
nullspace freedom, yet swap-sensitive to any genuine differential bug), and, since
Plan 17, **entry-by-entry equality mod p** of the collapsed differentials once the
two generator orders are aligned — byte-reproducible by construction. The bank is
read-only project law: these tests import it by path and never modify it.

### Self-certifying internal identities

Some facts need no external oracle because the defining axioms are the gate:

- `d∘d = 0` and the CS **order condition** are asserted on every CS resolution
  before its homology is read (`resolutions_cs/homology.py`).
- The native deep-degree cup (Plan 20) uses **Leibniz as the sign arbiter** —
  exact over GF(p) — with the transported cup as the in-window anchor; plus
  `(d^{P⊗P})² = 0`, chain-map, graded-commutativity and associativity gates.
- The native deep-degree **cap** (Plan 21) reads the *same* lifted diagonal the
  homology way; its sign convention is **arbitrated, not assumed** — the exact
  unit cap `1 ∩ z = z`, exact cap-Leibniz, the module identity
  `(z ∩ f) ∩ g ~ z ∩ (f ∪ g)` (via the native cup), and in-window agreement with
  the transported cap all hold simultaneously, and the non-commutative quantum CI
  distinguishes the correct `b·w·a` collapse from its `a·w·b` mirror
  (`tests/resolutions_cs/test_native_cap.py`).
- The **Plan-35 product tables** (public cup/cap/bracket + Connes `B`) are gated by
  the axioms of the structure they realize, with no external oracle: the cup table
  is **graded-commutative and associative**, the bracket is **antisymmetric** and
  satisfies **cup-Leibniz**, the cap table obeys the **module law**
  `(z∩f)∩g = z∩(f∪g)`, and the induced Connes `B` satisfies **`B²=0`** at the
  induced level with rank consistent with the `(b,B)` cyclic dims (SBI). Each
  identity is checked entry-by-entry over the prime set `{32003, 2, 3, 5}`, with the
  fixtures chosen so the check is content-bearing rather than a vacuous `0 = 0` (the
  cup-Leibniz sign, for instance, is pinned on `GF(3)` where `±1` differ)
  (`tests/hochschild/test_products_identities.py`, `test_connes_b.py`). **Plan 51**
  extends this to the CS-**native** bracket (homotopy liftings, any exact Domain, past
  the bar window): the defining equation (★), ψ consistency, descent, graded Jacobi,
  the native Poisson/cup-Leibniz law, and off-GF(p) antisymmetry are all self-certified
  class-level over GF primes AND QQ (`tests/resolutions_cs/test_homotopy_lifting.py`,
  `test_native_bracket.py`, `test_bracket_identities.py`), and the sign is arbitrated
  by the in-window native ≡ transported anchor (see the honest-scope section).
- The module layer's functors self-certify: `(A^op)^op ≅ A`, `D∘D ≅ id`,
  `D(P_v` over `A^op) ≅ I_v`, and `τ⁻τM ≅ M` for non-projective indecomposables
  via an exact invertible-hom certificate (`tests/modules/test_opposite.py`,
  `test_duality_tau.py`, `test_module_iso.py`).
- The comparison maps (Plan 14) are gated by chain-map and roundtrip identities.
- CS **canonicalization** (Plan 17) is pinned by an *adversarial-solver* test:
  shifting the correction solve by a nullspace vector must not move a single byte.
- The Nakayama form and automorphism (Plan 19) **self-certify** across domains:
  the returned `λ` and `ν` satisfy `λ(ab) = λ(b·ν(a))`, `ν` multiplicative,
  `ν(1) = 1`, and the Gram matrix is nondegenerate — checked exactly, no oracle
  needed.
- The **Plan-35 wave-3c Yoneda exact sequences** are gated by exactness itself: each
  `Ext^n(M,N)` class is constructed as `0 → N → Q → … → M → 0` (the pushout middle
  module + the induced/spliced connecting maps) and `check_exact` verifies every map
  is an `A`-module map, the ends are injective/surjective, and `im = ker` by rank at
  each interior joint — no external oracle. The kA₂ Baer extension is additionally
  pinned to the projective cover `0 → S₂ → P₁ → S₁ → 0` by the library's **own**
  `is_isomorphic`/`identify_standard`, and a non-cocycle input is refused loudly
  (`tests/modules/test_yoneda_p35.py`).

### Second-model oracles

Where a computation has a genuinely different classical model, that model is coded
independently and required to agree:

- **Cyclic homology** rests on Connes' `B`-operator and mixed complex (Connes,
  *Publ. Math. IHÉS* 62 (1985), 41–144, registry key `cyclic`). Over a field
  containing `Q` it is cross-checked against the **λ-complex** second model
  (Loday, *Cyclic Homology*, **Theorem 2.1.5** — the theorem number is recorded in
  the test docstring) — unnormalized chains, a quotient model, disjoint from both
  `hochschild/bar.py` and `hochschild/cyclic.py`
  (`tests/invariants/test_cyclic_generic.py`), plus the mixed-complex identities
  `b² = 0`, `B² = 0`, `bB + Bb = 0` over QQ.
- **Complexity** off GF(p) is the relative-Tor (Cibils) Betti complex, with
  `H_n =` the minimal resolution's rank sequence over every field, gated against
  the GF(p) engine (`tests/invariants/test_betti_generic.py`).
- **Frobenius / symmetric** off GF(p) is the socle-permutation criterion
  (Skowroński–Yamagata) with a verified socle-dual form
  (`tests/invariants/test_frobenius_generic.py`).
- **The Plan-38 C2 batteries** (2026-08-05). Coxeter: the worked examples of
  Armenta, *The Coxeter transformation as an automorphism of the
  Tamarkin–Tsygan calculus* (arXiv:2606.15595, registry-less bib key
  `armenta_coxeter_calculus`) pinned on the EXISTING Coxeter surface — `D₄` vs
  `A₄` differ by exactly `t²`, `A₄` is Φ₅, and the two 8-vertex cospectral trees
  share `(t+1)⁴(t⁴−3t³+t²−3t+1)` (`tests/invariants/test_coxeter_paper_pins.py`).
  Forms: Gabriel's theorem — the Euler form is `dim Hom − dim Ext¹` on a
  hereditary example and the arrow formula, and the Tits form's definiteness gives
  finite/tame/wild on `A₂`/Kronecker/3-Kronecker (`test_forms.py`). Roots: the
  classical positive-root counts `A_n:n(n+1)/2, D_n:n(n−1), E6:36` (Bourbaki;
  `test_roots.py`). Recognizers: textbook gentle/string/special-biserial/Nakayama/
  hereditary examples (ASS; Butler–Ringel; Assem–Skowroński gentle papers,
  `assem_book`; `test_recognizers.py`).

### Closed-form and chain-count pins

`k[x]/(x^a)` has a closed-form period-2 resolution; cyclic Nakayama has closed
Bardzell chain counts. The minimal engine's corner Betti numbers were
independently re-derived as Bardzell chain counts (`6, 5, 2, 1, 0` on
`kQ/(abc, cde)`), catching the Plan-12 straddling-chain bug.

---

## Class 2 — cross-engine, multi-prime, and QPA agreement

### Cross-engine degreewise agreement

| Battery | Compares | Over | File |
|---|---|---|---|
| CS ≡ bar | Chouhy–Solotar vs normalized bar | CC, GF(2/3/5) | `test_battery_bar.py` |
| CS ≡ Bardzell | Chouhy–Solotar vs Bardzell | GF(32003) | `test_battery_bardzell.py` |
| minimal ≡ bar | minimal `A^e` vs bar | `{32003, 2, 3, 5}` | `tests/engine/` |
| minimal-coh ≡ CS-coh | Hom-collapse vs CS cohomology | to depth 8 | `test_minimal_cohomology.py` |
| bar-cup ≡ CS-cup | GF(p) bar/tt cup table vs the Domain-generic CS cup table (basis-independent: dims + flattened rank mod p) | GF(3/7), in-window | `test_products_identities.py` |
| bar-bracket ≡ CS-bracket (Plan 51) | GF(p) bar/tt Gerstenhaber bracket table vs the CS homotopy-lifting bracket table (dims + flattened rank mod p) | GF(3/7), in-window | `test_products_identities.py` |
| native ≡ transported bracket (Plan 51 — THE sign anchor) | the CS homotopy-lifting bracket vs the transported bar bracket, mod coboundary; the flipped sign breaks it; odd-exponent QCI (2,4)→HH⁵ non-vacuous | k[x]/x² GF(32003), QCI GF(5)/GF(3), GF(2) | `test_native_bracket.py`, `test_bracket_literature_p51.py` |

Both differentials in each pair are built by disjoint code, so equal ranks and HH
dimensions are a real cross-check, not a tautology. The engines are the bar complex
(`bar` → Hochschild 1945), Bardzell's minimal resolution (`bardzell` → Bardzell,
*J. Algebra* 188 (1997), 69–89), Chouhy–Solotar (`chouhy_solotar` → 2015), and the
minimal projective `A^e` resolution (`minimal_resolution` → Green–Solberg–Zacharia,
*Trans. Amer. Math. Soc.* 353 (2001), 2915–2939).

### numba vs pure-Python parity

Each numba kernel is unit-tested against its pure-Python twin
(`tests/engine/test_kernels.py`), and — decisively — the **entire deep suite runs
twice in CI**, once with numba and once with `QUIVERLAB_NO_NUMBA=1`. The two paths
must agree exactly; parity is gated, not assumed.

### Live QPA / GAP cross-check

`A.crosscheck(...)` (`src/quiverlab/qpa/`) drives the GAP package **QPA** to
recompute independently and refuses to silently disagree
(`CrosscheckReport.assert_agree()` raises on mismatch). It covers:

- **Hochschild cohomology dims** via the enveloping algebra, `HH^n =
  Ext^n_{A^e}(A, A)` — QPA ships no HH function, so we build the route from
  `EnvelopingAlgebra` + `ExtAlgebraGenerators`. Pinned on the commutative square
  (`HH^• = [1, 0, 0]`, Künneth) over GF(2), GF(3), and QQ.
- **Module self-Ext** `Ext^*(M, M)` via `ExtAlgebraGenerators`, pinned on the
  simple `S_1` of `kA_2` (`[1, 0, 0]`).
- **Module theory** (Plans 23/24): the AR translates **τ/τ⁻** via `DTr`/`TrD`
  with `IsomorphicModules` on a translated module (dimension vectors *and* iso
  class — `modules/qpa_module.py::graded_form` handles QPA's row convention),
  **projective/injective resolution terms** via `ProjectiveResolution` /
  `DualOfModule`, and **injective dimension** via `InjDimensionOfModule`
  (`false ↔ None`), across the zoo including the multi-vertex records. **Left**-side
  quantities are crosschecked by feeding QPA the opposite algebra — QPA is
  right-module native (`tests/qpa/test_module_ar_crosscheck.py`,
  `tests/qpa/test_left_modules_qpa.py`).
- **Trivial-extension construction** (Plan 31): QPA 1.37 builds `T(A)` natively
  via `TrivialExtensionOfQuiverAlgebra`, so the crosscheck compares its
  dimension, arrow count, and the `IsSymmetricAlgebra` / `IsWeaklySymmetricAlgebra`
  / `IsSelfinjectiveAlgebra` predicates against quiverlab's presented `T(A)`
  (QPA's arrow labels differ — `te_a1_i_j` vs `te0` — so counts, not names, are
  compared). Pinned on `T(kA₂)` (dim 6, 2 arrows, all predicates true), `T(kA₃)`
  (12/3), the 2-Kronecker (8/4), the dual numbers (4/2), and the commutative
  square (18/5) (`tests/qpa/test_trivial_extension_qpa.py`).

The live QPA suite is `-m qpa` (242 tests). GAP is heavy to install, so it runs in a
**weekly** CI job, not on every commit — but it is **never silently green**: under
`QUIVERLAB_REQUIRE_QPA=1` an absent or broken QPA is a hard failure of that job,
and locally the tests skip explicitly rather than pass vacuously.

Separately, `tests/engine/test_qpa_reference_validation.py` freezes the values QPA
*would* produce (as published in the literature) and requires the bar engine to
reproduce them — an independent check that runs in the normal matrix with no GAP
present.

### Where QPA cannot be compared — and what covers that ground

QPA's cross-check reaches Hochschild dims, module self-Ext, and — since Plans
23/24 — the AR translates, projective/injective resolutions, and injective
dimension, over **QQ or a prime field GF(p)** (number-field CC and `GF(p^n)` are
out of QPA scope, and raise loudly). Everything below is therefore covered by a
**theory oracle**, not QPA:

| Feature QPA does not cover | Theory oracle that covers it |
|---|---|
| Cup / cap / Gerstenhaber bracket + the induced Connes `B` (Plan 35 — QPA 1.37 has **no** Hochschild product surface: no `CupProduct`/`HochschildCohomologyRing*`, confirmed by a live `NamesGVars()` sweep with zero `Hochschild`/`Cup` name; its `ExtAlgebraGenerators`/`YonedaProduct` is *module* Ext, not `HH^*(A)`) | the Gerstenhaber identity batteries (graded commutativity, associativity, Jacobi/antisymmetry, cup-Leibniz, cap module law `(z∩f)∩g = z∩(f∪g)`, `B²=0`, SBI rank consistency) + the `k[x]/(x²)` and QuantumCI-BGMS literature pins + the bar↔CS in-window cross-engine gate |
| The BV operator `Δ` on `HH^*` (Plan 54 — QPA 1.37 exposes **no** BV / Δ surface: a fortiori, since it has no Hochschild product ring at all; confirmed by a live `NamesGVars()` sweep with zero `Batalin`/`Vilkovisky`/`BVOperator` name in `tests/qpa/test_bv_qpa.py`, an honest skip that FAILS if one ever appears) | `Δ²=0` + the seven-term relation `(7T)` + the perfect-pairing certificate + the twisted-Connes descent self-cert + **the decisive cross-engine arbiter** (the bracket recovered from Δ equals the independent Gerstenhaber bracket in-window, no Frobenius input) + the Tradler `k[x,y]/(x²,y²)` and BIKLZ `k[x]/(x^N)` literature pins |
| Cyclic homology | Connes λ-complex second model + mixed-complex identities |
| The Chouhy–Solotar resolution | CS ≡ bar, CS ≡ Bardzell, and the bank byte-level closed forms |
| Deep degrees past the bar window | bank closed forms + cross-engine + closed-form/chain-count pins |
| Frobenius / Nakayama / symmetry | self-certifying `λ`/`ν` identities + socle criterion |
| `HH` over CC and `GF(p^n)` | exact bar oracle + second-model oracles (field-generic) |
| Distinct-module `Ext(M, N)`, `M ≠ N` | flagged post-v1; self-Ext is the confirmed QPA idiom |
| **The Koszul verdict itself** (Plan 27 — QPA 1.37 has no `IsKoszul`/`KoszulDual`, confirmed by an exhaustive `NamesGVars()` sweep) | the G-quadratic certifier (Priddy PBW: confluent length-2-tip reduction system ⇒ Koszul) plus the generated-in-degree-1 falsifier and the Fröberg matrix identity `P(t)·C_A(−t)=I`; QPA validates every INPUT to the verdict — graded Ext dims and minimal-generator degrees (`ExtAlgebraGenerators`), quadraticity (`IsQuadraticIdeal`), and the quadratic perp (`QuadraticPerpOfPathAlgebraIdeal`) |
| Yoneda relations-by-degree (QPA exposes generator counts, not a presentation) | theory battery: `E(k[x]/x²)=k[y]`, `E(k[x]/xⁿ)=k[y,z]/(y²)` (char-independent, pinned over GF(2)/GF(3)/GF(32003)/char 0), hereditary `E=kQ/J²` with `as_algebra()` round-trip, rad²=0 `E=kQ` path counts, quantum CI `dim Eⁿ=n+1` (= the CS chain count), commutative square `E≅A` self-hosting; the monomial Anick chain-count gate; byte-reproducible lift products (Plan-17-style canonicalization) |
| **The Coxeter spectral surface** (Plan 58 — the Mahler measure, spectral radius, cyclotomic-type verdict, and Lehmer-class note: QPA 1.37 has **no** such surface — `Cyclotomic*` are GAP polynomial *constructors*, not algebra predicates; confirmed by a live `NamesGVars()` sweep in `tests/qpa/test_coxeter_spectral_qpa.py::test_qpa_has_no_spectral_surface`, which FAILS if one ever appears). QPA's `CoxeterPolynomial(A)` IS present and **is** crosschecked exactly (variable rename `x_1 → t`) — the foundation the spectral layer rests on. | literature pins (de la Peña 2014/2013 — `dlPena2014mahler`/`dlPena2013cyclotomic`; Lehmer/E₁₀) + cross-engine (the shipped exact `spectral_radius`/`mahler_measure` primitives) + self-cert (the minimal polynomial annihilates ρ/M, the Sturm rational interval brackets the root, ρ = M when exactly one root lies outside the circle) |

### Live Macaulay2 cross-check

`quiverlab.m2.crosscheck` (`src/quiverlab/m2/`) drives **Macaulay2** — a
genuinely different computer-algebra system — to recompute two things
independently and refuses to silently disagree (the same
`CrosscheckReport.assert_agree()` container the QPA bridge uses). It is the
**fifth oracle class** (`-m m2`), external and independent like `-m qpa`:

- **Single-vertex graded dimensions** of `kQ/I` over `GF(p)` via M2's
  `AssociativeAlgebras` package — an independent **noncommutative Gröbner (F4)**
  engine. The script builds `kk<|x,y,…|>/(rels)` and reads
  `dim_k B_n = numgens source ncBasis(n, B)` for `n = 0..top`; our side is the
  Hilbert data of the same algebra from its own reduction-system tips
  (`modules/koszul::_algebra_graded_matrices`). Because M2's Gröbner engine is
  written and maintained entirely separately from ours, an agreement is a
  cross-implementation check of the whole tip/normal-form stack.
- **Commutative Ext dimensions** — for a commutative example
  `k[x,y]/(relations)` (the presentation must carry the explicit commutator,
  checked textually and refused loudly otherwise), M2's `freeResolution` of the
  residue field over `ZZ/p[x,y]/(relations)` gives the graded Betti numbers
  `rank C_n`, compared against `A.ext_algebra(top).graded_dims_through(top)` —
  a fully independent homological route.

The transport is a **subprocess**, not an in-process library: each call writes
the script to a temp file and runs `M2 --script file`, parsing sentinel lines
`<<QL>> n v` back through exact-integer parsing (no floats cross the boundary).
Version policy: any Macaulay2 `≥ 1.24` with the bundled `AssociativeAlgebras` +
`Complexes` packages; the CI job pins the Ubuntu PPA build and fails (never
skips) when M2 is absent under `QUIVERLAB_REQUIRE_M2=1`. M2 sees **no**
multi-vertex algebra and **no** Hochschild anything — those requests are refused
loudly (see [Honest scope](#honest-scope)).

---

## Other structural gates

- **Float ban (AST gate).** `tests/test_no_floats.py` walks the AST of every file
  under `src/` and fails on any float/complex literal or `float()` call — and a
  planted-violation test proves the gate itself works. Non-algebraic glue
  (`webapp/`, `docs/gui/`) is exempt by design; it holds no exact algebra.
- **Golden-file worked-steps traces.** The verbose worked-steps documents
  (`trace/`) are golden-file tested: every dimension they print is *derived* from a
  recorded rank, so a printed claim can never drift from the computed value.
- **Admissibility certificates.** The Gröbner engine (`groebner/`) certifies
  finite-dimensionality (the `2L−1 ≤ D` bound plus the forbidden-word automaton);
  a non-admissible or infinite presentation fails loudly, never hangs
  (`tests/groebner/test_certificate.py`).
- **Refusal surfaces are tested.** Out-of-scope inputs (structure-constant
  algebras off GF(p) needing a path basis; an inconsistent CS correction solve;
  cross-mode checkpoint reuse) raise named errors — and those refusals are
  asserted, so they stay loud (`tests/invariants/test_refusal_surface.py`,
  `tests/engine/test_error_paths.py`).

---

## Subsystem → oracles → tests

Every `src/quiverlab/` subpackage, its test directory, the collected test count,
and the oracle class that guards it. Counts are `pytest --collect-only` with the
`[dev,fast,docs,web,qpa]` extras (2026-07-25 baseline, post-merge of Plans 21–26;
the `hochschild/` and `resolutions_cs/` rows are refreshed 2026-08-01 for the
Plan-35 product surface).

| Subsystem (`src/quiverlab/`) | Tests | Bucket | Primary oracle class |
|---|---:|---|---|
| `fields/` (QQ, GF(p), GF(p^n), exact CC = QQ_I) | 41 | fast | exact-arithmetic axioms; base-change invariance |
| `core/` + `combinat/` (Quiver, Algebra, relations, dispatch; the Plan-44 `basic.py` — `primitive_idempotents`/`basic_algebra`/`gabriel_quiver`/`presented_form`, the exact Wedderburn/trace-form recovery of a structure-constant algebra, batteried in `tests/families/test_gabriel_recovery.py`) | 43 | fast | structure-constant identities; left-to-right path law; the Gabriel-recovery certificates live in the deep `families/` bucket |
| `groebner/` (overlap completion, admissibility) | 50 | fast | admissibility certificate; finiteness; lowering |
| `hochschild/` (bar, cyclic; the Plan-34 auto→CS depth-fallback battery; the Plan-35 product surface — `products.py`: cup/cap/bracket tables + the induced Connes `B`, `basis_reps.py`: the explicit-representatives capture, `cyclic_reps.py`: the wave-3b cyclic-homology explicit representatives, and `hh_reps.py`: the wave-3d plain-HH explicit representatives) | 107 | fast | **the base bar oracle**; mixed-complex identities; dispatch-amendment pins; **the Gerstenhaber identity batteries** (graded-commutative + associative cup, antisymmetric bracket, cup-Leibniz, cap module law, `B²=0`, SBI rank) + the `k[x]/(x²)`/QuantumCI-BGMS product literature pins + the bar↔CS in-window cross-engine gate + **the explicit-reps self-certification** (every shipped product class satisfies `δ·v = 0` / `b·v = 0` from its shipped or note-rebuilt differential; hand-checked `k[x]/(x²)` labels; elision+rebuild path) + **the cyclic-homology explicit representatives** (every shipped HC class is a cycle of the (b,B) total complex — `D·v = 0` on both engines; GF(p)≡generic on prime 32003; hand-checked dual-numbers `HC_0 = A/[A,A]` + `Tot_2 = C_2 ⊕ C_0` column structure) + **the plain-HH explicit representatives** (the `hh_cohomology`/`hh_homology` dims blocks now carry per-degree reps over BOTH routes — the GF(p) bar and Chouhy-Solotar; every shipped vector annihilates its shipped differential; the `k[x]/(x²)` HH¹ `[x ↦ x]` = derivation `D(x)=x` hand-check ties the classical dictionary read-off to the captured representative) |
| `hochschild/bv/` (Plan-54 R2 — the Batalin–Vilkovisky operator `Δ: HHⁿ → HHⁿ⁻¹` on a Frobenius / symmetric algebra with semisimple Nakayama automorphism: `hypothesis.py` the decidable ν gate (minpoly-squarefree semisimplicity + the symmetric anchor + loud BIKLZ-blocked / non-Frobenius refusals), `transport.py` the Frobenius perfect pairing `(†)` + the adjoint transport `Δ_n = (P_{n-1}ᵀ)⁻¹ Bᵀ P_nᵀ` (symmetric **and** semisimple-ν routes), `twist.py` the P52 twisted-homology adapter (the cross-plan `twisted_homology_classes` contract), `twisted_connes.py` the twisted Connes `B_σ` that descends to twisted homology, `bracket.py` the bracket recovered from Δ via the BV relation + the arbiter) | 52 | fast | **`oracle_selfcert`**: `Δ²=0`, the seven-term relation `(7T)` (Δ is a differential operator of order ≤ 2), the perfect-pairing certificate (`P_n` square + invertible in-window), the twisted-Connes DESCENT self-cert (`B_σ` carries twisted cycles to cycles — a loud refusal otherwise), **GF(2) served** (Δ²=0 + BV-relation mod 2 + pairing-invertibility; the sign is inherited from the odd-prime anchor, not char-2-determined). **`oracle_crossengine`** (the decisive gate): the bracket recovered from Δ via `(BV)` equals the INDEPENDENT `gerstenhaber_brackets` table-for-table in-window over ODD primes (3/5/32003) — no Frobenius input in the independent bracket — which also SELECTS the twist direction (`QuantumCI(q=2/3)` semisimple diagonal ν, the non-trivial-twist flagship). **`oracle_literature`**: the Tradler `k[x,y]/(x²,y²)` ([4,4,5,6]) anchor, the `k[x]/(x^N)` char-sensitivity (dim HH⁰=N; HHⁿ=N−1 if char∤N else N) in BOTH regimes, and the BIKLZ §3.2 `kZ₁/J^N` value oracle (the exact char-sensitive dims + the Δ-rank profile `rank Δ_odd = N−1, rank Δ_even = 0` that the §3.2 formulas force — basis-independent, live). QPA has NO BV surface (`tests/qpa/test_bv_qpa.py`, honest skip that FAILS on change) |
| `engine/` (fast GF(p); minimal, Bardzell, periodic; TT-calculus; cyclic; Coxeter/Nakayama; Plan-29 literature/identity batteries) | 579 | deep | bar oracle; cross-engine; multi-prime; numba/pure parity; frozen QPA-literature values |
| `resolutions_cs/` (CS; comparison; diagonal; cup; cap; Plan-29 literature batteries; the Plan-35 Domain-generic CS product tables `products.py` — cup/cap on the CS basis over any exact Domain; the Plan-35 wave-3d `cs_hh_basis` 0-codomain fix; the Plan-51 native Gerstenhaber bracket — `homotopy_lifting.py` (the per-corner `D_corner` + the Negron–Witherspoon/Volkov tower ψ) + `bracket.py` (`native_bracket`) + `cs_bracket_tables`) | 275 | deep | CS ≡ bar, CS ≡ Bardzell; bank byte-level; literature pins; `d∘d=0` / order; Leibniz + cap identities (unit/module/transport anchors); canonicalization; the CS product unit-law + Domain-genericity self-cert; **the `cs_hh_basis` rep-count ≡ `cs_(co)homology_dims`** even when the top differential lands in a 0-dimensional space (the multi-vertex `kZ_3/J^2` witness — HH¹ = 1 with C² = 0 — that the old 0-row `nullspace` silently undercounted); **the Plan-51 bracket oracles** — `oracle_selfcert`: `D_corner` d²=0, the homotopy-lifting equation (★) exactly (GF(p)+QQ), ψ consistency at every built degree + the `CocycleError` guard, bracket descent, graded Jacobi, native Poisson/cup-Leibniz, off-GF(p) antisymmetry; `oracle_crossengine`: the in-window native ≡ transported bracket anchor (THE sign fix; QCI (2,4)→HH⁵ odd-exponent, flipped-sign breaks it) + bar ≡ CS bracket tables; `oracle_literature`: k[x]/xⁿ off-GF(p) zeros + `[α,β]≠0`, the QCI dims-line-up (self-cert, not a value), and the Oke §7 blocked-until-transcribed strict-xfail |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity; Plan-29 Tor; Plan-30 Krull–Schmidt decomposition; the retained injective-coresolution differentials certified exact; the Plan-35 wave-3a explicit Ext/Tor representatives — `complex_reps.py`; the Plan-35 wave-3c Yoneda exact sequences — `yoneda.py`; the Plan-37 C1 categorical glue — `morphism.py` first-class `ModuleHom` + kernel/image/cokernel, `ses.py` short exact sequences + split test + pushout/pullback, `endomorphism.py` `End(M)` as an Algebra, `direct_sum`/`is_direct_summand`, and covers/envelopes + radical/socle series + composition factors on `Module`; the Plan-39 C8 complex layer — `complexes.py` validated bounded `ChainComplex`/`ChainMap`, shift/truncation/homology, mapping cones + triangles + the derived-iso test, the Hom total complex `hyper_hom_dims`, and certified projective models `projective_model` + general `hyper_ext_dims`) | 317 | deep | AR/duality literature pins (ASS2006); **the Plan-37 categorical-glue self-certification** (`ModuleHom` validates the intertwining relations at construction; kernel/image/cokernel certified by rank-nullity + the epi–mono factorization `f = epi∘mono` + `f∘iota = 0` = `proj∘f`; SES exactness = the rank identity `im f = ker g`; split ⇔ a section solves; pushout/pullback squares certified by their universal-square identities; `End(M)` self-certified by `from_structure_constants(check=True)` with the regular-module `End(A_A) ≅ A` Loewy oracle; biproduct identities `proj_i∘incl_i = id`, `Σ incl_i∘proj_i = id`); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks; **the explicit Ext/Tor self-certification** (every shipped class satisfies `δ·v = 0` (Ext cocycle) / `d·v = 0` (Tor cycle) from its shipped differential; hand-checked kA₂ `Ext¹(S₁,S₂)` + loop `Tor₀ = M ⊗ N` cokernel labels; rep-count ≡ engine dims) + **the Yoneda exact-sequence self-certification** (every `Ext^n(M,N)` class is CONSTRUCTED as an `n`-fold exact sequence `0 → N → Q → … → M → 0` — the pushout middle module + connecting maps — and its exactness is self-certified at every joint: each map an `A`-module map, ends injective/surjective, `im = ker` by rank; the kA₂ Baer pin `0 → S₂ → P₁ → S₁ → 0` verified by the library's OWN `is_isomorphic`/`identify_standard`; a non-cocycle is refused loudly; the multi-vertex `_tor_boundary` collapse pinned on a rad²=0 Nakayama) + **the Plan-39 complex-layer oracles** — self-cert: `d∘d=0` refused at construction, the mapping-cone `d²=0` re-asserted under full validation, quasi-iso ⇔ cone acyclicity, shift/truncate degree+sign identities, and the Hom total-complex `δ∘δ=0` block-indexing certificate (the Weibel `−(−1)^n` sign is a documented convention — verified sign-independent for the dims: both `±` give isomorphic cochain complexes); cross-engine: stalk `hyper_hom_dims` ≡ module `Ext` degreewise (kA₃/(ab), all vertices), the resolution-augmentation quasi-iso, the Euler-characteristic triangle identity `χ(cone) = −χ(X)+χ(Y)`, the certified projective model (`projective_model` asserts `is_perfect` + `is_quasi_iso` before return — NEVER returns uncertified), and the sharpened two-term shift identity `hyper_ext([P₁→S₁])[n] = Ext^{n−1}(rad P₁, N)` (the derived-category source-shift arithmetic, replacing the plan's placeholder); live QPA Ch.10 (`qpa/`) |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity; Plan-29 Tor; Plan-30 Krull–Schmidt decomposition; the retained injective-coresolution differentials certified exact; the Plan-35 wave-3a explicit Ext/Tor representatives — `complex_reps.py`; the Plan-35 wave-3c Yoneda exact sequences — `yoneda.py`; the Plan-37 C1 categorical glue — `morphism.py` first-class `ModuleHom` + kernel/image/cokernel, `ses.py` short exact sequences + split test + pushout/pullback, `endomorphism.py` `End(M)` as an Algebra, `direct_sum`/`is_direct_summand`, and covers/envelopes + radical/socle series + composition factors on `Module`; the Plan-40 C6 homological-dimensions family — `homdims.py`: public `syzygy`/`cosyzygy` (byte-stable extraction from `minimal_resolution`), the Igusa–Todorov φ/ψ on the finite K₀, dominant + Gorenstein dimensions, Ω/τ-periodicity certificates, and finitistic-dimension bounds) | 322 | deep | AR/duality literature pins (ASS2006); **the Plan-37 categorical-glue self-certification** (`ModuleHom` validates the intertwining relations at construction; kernel/image/cokernel certified by rank-nullity + the epi–mono factorization `f = epi∘mono` + `f∘iota = 0` = `proj∘f`; SES exactness = the rank identity `im f = ker g`; split ⇔ a section solves; pushout/pullback squares certified by their universal-square identities; `End(M)` self-certified by `from_structure_constants(check=True)` with the regular-module `End(A_A) ≅ A` Loewy oracle; biproduct identities `proj_i∘incl_i = id`, `Σ incl_i∘proj_i = id`); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks; **the explicit Ext/Tor self-certification** (every shipped class satisfies `δ·v = 0` (Ext cocycle) / `d·v = 0` (Tor cycle) from its shipped differential; hand-checked kA₂ `Ext¹(S₁,S₂)` + loop `Tor₀ = M ⊗ N` cokernel labels; rep-count ≡ engine dims) + **the Yoneda exact-sequence self-certification** (every `Ext^n(M,N)` class is CONSTRUCTED as an `n`-fold exact sequence `0 → N → Q → … → M → 0` — the pushout middle module + connecting maps — and its exactness is self-certified at every joint: each map an `A`-module map, ends injective/surjective, `im = ker` by rank; the kA₂ Baer pin `0 → S₂ → P₁ → S₁ → 0` verified by the library's OWN `is_isomorphic`/`identify_standard`; a non-cocycle is refused loudly; the multi-vertex `_tor_boundary` collapse pinned on a rad²=0 Nakayama) + **the Plan-40 homological-dimensions oracles**: `oracle_selfcert` — the φ=ψ=pd identity for finite projective dimension, the Ω/τ-periodicity `is_isomorphic` certificates, and the decompose char-caveat propagation; `oracle_literature` — the Barrios–Mata truncated self-injective φ=ψ=0 closed form + projective additivity, the hereditary/self-injective dominant & Gorenstein values, and the cyclic-Nakayama period-from-Kupisch pins; `qpa` (`tests/qpa/test_homdims_qpa.py`) — live `GlobalDimensionOfAlgebra` / `DominantDimensionOfAlgebra` / `GorensteinDimensionOfAlgebra` agreement over kA₂ / kA₃(ab) / `line_abc_cde` / k[x]/(x³) (int or GAP `infinity`↔our infinite/unresolved marker) |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity; Plan-29 Tor; Plan-30 Krull–Schmidt decomposition; the retained injective-coresolution differentials certified exact; the Plan-35 wave-3a explicit Ext/Tor representatives — `complex_reps.py`; the Plan-35 wave-3c Yoneda exact sequences — `yoneda.py`; the Plan-37 C1 categorical glue — `morphism.py` first-class `ModuleHom` + kernel/image/cokernel, `ses.py` short exact sequences + split test + pushout/pullback, `endomorphism.py` `End(M)` as an Algebra, `direct_sum`/`is_direct_summand`, and covers/envelopes + radical/socle series + composition factors on `Module`; the Plan-41 C3 **Auslander–Reiten completion** — `ar.py`: the general chain-map lift, the Nakayama functor ν/ν⁻, stable Hom mod projectives, the End(M)-action on Ext¹, almost-split sequences `0 → τM → E → M → 0`, irreducible-map multiplicities `dim rad(M,N)/rad²`, and honest-semi-decision AR-quiver knitting) | 327 | deep | AR/duality literature pins (ASS2006, ARS1995); **the Plan-41 AR self-certification** (the chain-map lift asserts every square `d_n·φ_n = φ_{n-1}·d_n` and is byte-reproducible; ν tied to the trusted τ by `ker(νP₁→νP₀) ≅ τM` + `ν(P_v) ≅ I_v` / `ν⁻(I_v) ≅ P_v`; every almost-split sequence's class is identified by the **ARS socle theorem** — `soc_{End M} Ext¹(M,τM)` is a *simple* `End(M)`-module — pinned down computationally by the char-scoped trace-form `rad End(M)` (char 0 or char > dim M) plus the **socle-simplicity dimension certificate** `dim_k soc = dim_k End(M) − dim_k rad End(M)` (`= dim_k` the residue division algebra; refuses loudly if it fails); the exact (Yoneda + P37 SES) + non-split (P37 `is_split` False) + indecomposable-ends (Plan 30) checks are **necessary sanity checks, not sufficient** — they do NOT arbitrate the pick (over `k[x]/(x⁴)` with `M = k[x]/(x²)`, `dim Ext¹ = 2`, a non-socle class has an exact/non-split/indecomposable-ends extension whose middle is the *projective* `k[x]/(x⁴)`, not the true mesh middle `{1,3}` — the devil's-advocate non-brick regression, live-QPA-crosschecked)) + **the cross-engine AR formula** `dim Ext¹(M,N) = dim underline-Hom(τ⁻N, M)` + **literature Dynkin/Nakayama pins** (kA₂/kA₃/kA₄ = 3/6/10 indecomposables, D₄ = 12, the kA₃ mesh middle terms, Nakayama serial count) + **live QPA** `AlmostSplitSequence` middle-term (dim vector over QQ; summand multiset over GF(p)) and `PredecessorsOfModule`; **the Plan-37 categorical-glue self-certification** (`ModuleHom` validates the intertwining relations at construction; kernel/image/cokernel certified by rank-nullity + the epi–mono factorization `f = epi∘mono` + `f∘iota = 0` = `proj∘f`; SES exactness = the rank identity `im f = ker g`; split ⇔ a section solves; pushout/pullback squares certified by their universal-square identities; `End(M)` self-certified by `from_structure_constants(check=True)` with the regular-module `End(A_A) ≅ A` Loewy oracle; biproduct identities `proj_i∘incl_i = id`, `Σ incl_i∘proj_i = id`); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks; **the explicit Ext/Tor self-certification** (every shipped class satisfies `δ·v = 0` (Ext cocycle) / `d·v = 0` (Tor cycle) from its shipped differential; hand-checked kA₂ `Ext¹(S₁,S₂)` + loop `Tor₀ = M ⊗ N` cokernel labels; rep-count ≡ engine dims) + **the Yoneda exact-sequence self-certification** (every `Ext^n(M,N)` class is CONSTRUCTED as an `n`-fold exact sequence `0 → N → Q → … → M → 0` — the pushout middle module + connecting maps — and its exactness is self-certified at every joint: each map an `A`-module map, ends injective/surjective, `im = ker` by rank; the kA₂ Baer pin `0 → S₂ → P₁ → S₁ → 0` verified by the library's OWN `is_isomorphic`/`identify_standard`; a non-cocycle is refused loudly; the multi-vertex `_tor_boundary` collapse pinned on a rad²=0 Nakayama) |
| `derived/` (Plan-43 C8 **derived-category surface** — `homs.py`: reified hyper-Hom classes `hyper_hom_basis` (a basis of `H^n(Hom^•(X,Y))` as genuine chain maps `X → Y[n]`) + `ChainMap.then`; `tau.py`: the derived AR translate `τ_{D^b} = ν[−1]` / `τ⁻_{D^b}` on perfect complexes with the Happel finite-gl.dim gate; `_corner.py`: the shared corner-transpose `Hom_A(−,A)` (factored out of `duality._presentation_transpose`, imported by both); `tilting.py`: the tilting-complex verifier (rigidity DECIDED on the exact window; generation is honest **three-valued** — see the Plan-67 fix-round note below: `is_tilting` is a hard `True` only when generation is CERTIFIED (2-term/regular, IJY completion), `"unknown"` for a non-2-term rigid K₀-basis object since `rigid + det(g_proj)=±1 ⇒ tilting` is Rickard's OPEN rank question) + `End(T)` as a structure-constant algebra + `corner_cartan_of_complex` + `two_term_silting_from_presentation`; `fingerprint.py` + `block.py`: the necessary-condition derived fingerprint) | 24 | deep | **self-cert** — `hyper_hom_basis` reifies genuine chain maps (`ChainMap(check=True)` on every coset representative) with count ≡ `hyper_hom_dims`; `τ_{D^b}` output `d∘d=0` (ChainComplex check) and the `τ⁻∘τ` round-trip is a degreewise quasi-iso; `End(T)` self-certified by `from_structure_constants(check=True)`; the tilting rigidity window reported honestly; the **three-valued generation** — `generation == "certified"` for the regular object + APR/2-term tilts (hard `True`), `"k0_necessary_only"` ⇒ `is_tilting == "unknown"` for a genuine width-2 rigid K₀-basis complex over `k[1→2→3]/(ab)` (that same object verifies `"unknown"` under `is_silting_object` — the tilting rung and the K₀-only rung reconciled); **cross-engine** — `hyper_hom_basis` count ≡ module `Ext^n` on a projective-resolution source (kA₃/(ab), all vertices), and `τ_{D^b}` homology is concentrated in degree 0 and `≅` the trusted module `τ` over `kA_n`; **literature** — the K₀-action identity `χ(τ_{D^b}X) = c·χ(X)` with `c = −C·C⁻ᵀ` (the conjugate of P38's Coxeter matrix — same char poly, the dim-vector action), the kA₂ APR-tilt `P₁ ⊕ S₁` (tilting, `End(T)` = the reoriented A₂ = A^op, corner-Cartan `[[1,0],[1,1]]` — the theorem-anchored `End(A_A)=A` pin fixes the orientation), the D₄ vs A₄ Coxeter distinction and the 8-vertex cospectral-trees NON-distinction (the honest-scope demonstration); **live QPA** (`tests/qpa/test_derived_qpa.py`) — `τ_{D^b}` homology(0) vs `DTr(M)` (the documented module-level route: QPA's `TauOfComplex` on a `ProjectiveResolution` does not script through libgap — the P39 Ch.10 hazard, confirmed live) |
| `invariants/` (Cartan, Coxeter, spectral, Betti, cyclic, Frobenius incl. the Plan-29 trace-form symmetry certifier, scalar, sweep; Plan-29 Coxeter/identity literature batteries; the Plan-58 R20 **certified Coxeter spectral analysis** — `spectral.py::certify_real_algebraic`/`cyclotomic_factorization`/`off_circle_root_count` + `coxeter_spectral.py::coxeter_spectral`/`coxeter_spectral_block`) — `tests/invariants/test_coxeter_certified.py`, `test_coxeter_spectral.py`, `test_coxeter_spectral_literature.py` | 135 | fast | second models (λ-complex, relative-Tor Betti); self-certifying `λ`/`ν`; GF(p) engine parity; **the Plan-58 certified-spectral oracles** — `oracle_selfcert`: the certificate's minpoly round-trips `sympy.minimal_polynomial` (the independent anchoring is the literature pins), the rational Sturm interval brackets the root and the index is consistent, the Φ_n labelling is exact, the **product-Mahler M ≠ ρ branch** (two roots outside the circle: χ = (t²−7t+1)(t²−14t+1), degree-4 minpoly [1,−98,243,−98,1], interval (95,96)), the deterministic complex-dominant gate refuses ρ/M without ever calling `minimal_polynomial` (the measured 121 s anti-hang), and a singular/absent Cartan refuses per-field (never a crash); `oracle_literature`: kA₂ = Φ₃ (order 3), kA₄ = Φ₅ (order 5), kA₅ = Φ₂Φ₃Φ₆, the 3-Kronecker ρ = M = (7+3√5)/2 (minpoly [1,−7,1], interval (6,7)), the m-Kronecker ladder x²−(m²−2)x+1, T_{2,3,7} = E₁₀ = **Lehmer** (minpoly = the degree-10 polynomial, interval (1,2)), and the ADE/affine cyclotomic sweep (finite Coxeter order = Coxeter number on Dynkin, None/Jordan on affine) — pins from de la Peña 2014/2013 (`dlPena2014mahler`/`dlPena2013cyclotomic`) + de la Peña–Takane 1990 (`dlPenaTakane1990spectral`); `oracle_crossengine`: the report's ρ/M match the shipped `spectral_radius`/`mahler_measure` primitives |
| `invariants/geometry.py` (Plan-49 C8 — orbit dimension `dim O_M = Σ d_v² − dim End(M)`, Voigt rigidity `is_rigid`/`rigidity_codim`, the Kac `canonical_decomposition` over hereditary Dynkin, and the shared `orbit_geometry_block`) — `tests/invariants/test_geometry_orbit.py`, `test_geometry_canonical.py` | 21 | fast | **`oracle_selfcert`**: the orbit-dim identity `dim O_M = Σ d_v² − dim End(M)`, GF(p)↔QQ field parity, the canonical-decomposition sum-of-roots + per-instance rigidity certificate `Ext¹(G,G)=0`, the loud Euclidean-deferred / non-hereditary refusals. **`oracle_crossengine`**: the Voigt codim identity `dim Rep(Q,d) − dim O_M ≡ dim Ext¹(M,M)` on hereditary + the P38 `tits_form` tie; the canonical decomposition ≡ the Krull–Schmidt summands of the degeneration poset's maximum. **`oracle_literature`**: every Dynkin indecomposable is rigid (codim 0); `(2,1) = P₁ ⊕ S₁` over kA₂ (hand-derived Kac pin) |
| `modules/degeneration.py` (Plan-49 C8 — the Zwara–Bongartz degeneration = hom order poset for representation-finite algebras, `DegenerationPoset`) — `tests/modules/test_degeneration.py` | 6 | deep | **`oracle_literature`**: kA₂ (1,1) = the 2-chain `S₁⊕S₂ <_deg P₁`; kA₃ (1,1,1) = the diamond (orbit dims `[0,1,1,2]`, 4 covers, two incomparable middles) — both hand-derived. **`oracle_selfcert`**: the hom-order is a partial order (reflexive/antisymmetric), orbit dim strictly increases up every cover, the per-class orbit dim ≡ `geometry.orbit_dimension`, and the honest semi-decision cap (rep-infinite / self-injective ⇒ `is_complete=False` with a loud `status`, never a partial poset) |
| `modules/left_right.py` (Plan-55 R15 — the left/right parts `L_A`/`R_A` of the module category via the closed-under-predecessors pd/id ≤ 1 sweep on the knitted AR quiver, the finite complement `ind A ∖ (L_A ∪ R_A)`, the Ext-injectives of `add L_A` (+ dual Ext-projectives of `add R_A`), and the left/right support algebras `A_λ`/`A_ρ` as presented induced-convex-subquiver `Algebra`s with their connected-component factors; the P61-addendum per-indecomposable `pd_le_1`/`id_le_1` vectors) — `tests/modules/test_left_right_{parts,ext,support,oracles}.py`, `tests/qpa/test_left_right_qpa.py` | 24 + 2 qpa | deep + qpa | **`oracle_selfcert`** — `L_A` closed under predecessors / `R_A` under successors, intersection/complement consistency (`complement = U ∖ (L_A ∪ R_A)`), `gl.dim ≤ 1 ⇒ both parts total`; Ext-injectives ⊆ `L_A` and an injective in `L_A` is always Ext-injective; the support certificates — `e_λ`/`e_ρ` convex, hereditary ⇒ `A_λ = A_ρ = A` (connected), and the loud presentation-less refusal (structure-constants-only `A` ⇒ `QuiverlabError`, no fabricated quiver); the honest refusals — self-injective and **provably rep-infinite** (hereditary non-Dynkin, or rad²=0 with a non-Dynkin separated quiver — Gabriel, covering the 2-Kronecker and ACLV 2.2(c)) ⇒ `status="unsupported"` fast, any other rep-infinite input ⇒ the `budget_dim`-capped knit's `status="budget"`; the P61 `pd_le_1`/`id_le_1` index-alignment. **`oracle_literature`** — hereditary ⇒ `L_A = R_A = ind A`, empty complement (kA_n); **ACLV Example 2.2(b)** the rad²=0 linear-Nakayama A₅: `L_A = {S₁,S₂,P₂,P₃}`, `R_A = {S₄,S₅,P₄,P₅}`, `L_A ∩ R_A = ∅`, **complement = {S₃}** (pd 2, id 2 — the ada-with-non-empty-complement datum), `e_λ = {1,2,3}`, `e_ρ = {3,4,5}`, `gl.dim = 4`; kA_n Ext-injectives count = #injectives. **`oracle_crossengine`** — the Hom-nonzero transitive-closure predecessor relation ≡ the AR-quiver irreducible-map reachability closure (rad^∞ = 0, rep-finite); the ACLV duality `D R_A = L_{A^op}` / `D L_A = R_{A^op}` on the **proper-subset** A₅ (`\|R_A\| = 4` of 9, non-vacuous); the presented `A_λ` dim ≡ `end_algebra(⊕ P_x)` dim. **`qpa`** (`tests/qpa/test_left_right_qpa.py`) — the defining pd ≤ 1 / id ≤ 1 flags corroborated pointwise via QPA `ProjectiveResolution` / `InjDimensionOfModule` on the kA₃ indecomposables; a fail-if-appears `IsBoundGlobal` probe that trips if QPA ever ships a left/right-part or support-algebra verb |
| `modules/recognizers_ladder.py` (Plan-61 R18 — the quasi-tilted / shod / weakly-shod / laura / ada recognizer ladder off the P55 atlas + `global_dimension` + one AR-SCC weakly-shod sweep, each rung a definite verdict with a witness-on-`False` and a certificate-on-`True`; the monotone-nesting + gl.dim theorem-gate self-certs; and — for ada algebras over an algebraically closed field — the shipped `HH¹` as a COMPLETE simple-connectedness oracle via ACLV Theorem B, gated on the new additive `Domain.is_algebraically_closed` flag) — `tests/modules/test_recognizer_ladder{,_hh1,_oracles}.py`, `tests/modules/test_ada_hh1_p56_agreement.py`, `tests/qpa/test_recognizer_ladder_qpa.py` | 28 + 4 qpa | deep + qpa | **`oracle_literature`** — **ACLV Example 2.2(b)** the rad²=0 linear-Nakayama A₅: `(quasi_tilted, shod, weakly_shod, laura, ada) = (False, False, True, True, True)`, `complement = {S₃}` (pd 2, id 2), `gl.dim = 4` — a full five-way discrimination on ONE fixture; hereditary `kA_n` ⇒ all five `True`, empty complement; **`weakly_shod = False` on `NakayamaAlgebra(kupisch=[3,2,2])`** (one non-trivial AR-SCC on an injective→projective route, verified live — the seam maximally exercised: ada/shod/quasi_tilted all `False`, laura `True`, gl.dim 3); the ada/HH¹ tree fixture over CC (ada, HH¹ = 0, simply connected — Theorem B); the **REAL Theorem-B `False` pin** — the "square" `kQ/(a·b, c·d)` on `1→2→4`, `1→3→4` over CC is rep-finite and quasi-tilted (hence ada, gl.dim 2, empty complement) yet has `dim HH¹ = 1` and `π₁ = ℤ`, so it is NOT simply connected and Theorem B returns `verdict = False` (verified live — a genuine rep-finite ada + HH¹ ≠ 0 instance, not only the synthetic helper branch); the alg-closed gate predicate on constructed algebras (`field=CC` ⇒ `True`; `QQ`/`GF(7)`/`QQi` ⇒ `False`, off `A.domain` never the `CC` sentinel). **`oracle_selfcert`** — the monotone nesting (`quasi_tilted ⇒ shod ⇒ weakly_shod ⇒ laura`, `quasi_tilted/shod ⇒ ada`, `ada ⇒ laura` in scope) + the **guard-RAISES** unit test (an inconsistent rung set is a loud `QuiverlabError`); the gl.dim theorem gates (`shod ⇒ gl.dim ≤ 3`, `quasi_tilted ⇒ gl.dim ≤ 2`) — exercised AT the `shod ⇒ gl.dim ≤ 3` boundary on the strict-shod rad²=0 linear A₄ (`shod` `True` at gl.dim exactly 3, `quasi_tilted` `False`), not merely vacuously below it; the **inexact-gl.dim honesty** on `NakayamaAlgebra(kupisch=[3,3,2])` (rep-finite so the ladder completes, but `global_dimension` returns `value=32, exact=False` — a certified lower bound; the block carries `gldim_exact=False` and NEVER presents 32 as a definite gl.dim, and every gl.dim-thresholded verdict — `quasi_tilted`/`shod`/`ada` — resolves to `False` by its own logic since an inexact value is always ≥ the resolution bound `≥ 32 > 3 > 2`, never by trusting a spurious finite number); the strict-weakly-shod `HH^{≥2} = 0` (survey §5.1) + Corollary 2.6 (`ada ⇒ gl.dim ≤ 4`, `pd ≤ 2 or id ≤ 1`); the pure `_theorem_b_verdict` gate-logic branches across the four `(ada, alg_closed)` corners (the end-to-end `False`-branch is the REAL `oracle_literature` pin above); the self-injective / rep-infinite (ACLV 2.2(c)) loud refusals — no partial ladder; a real non-ada witness (the triangle `1→2→3`+`1→3`, `a·b=0`). **`oracle_crossengine`** — the quasi-tilted two routes agree (HRS `gl.dim ≤ 2 ∧ complement = []` ≡ ACLV `every P_x ∈ L_A`); the shod two routes agree (complement-empty ≡ the QT2 `pd_le_1`/`id_le_1` sweep); the weakly-shod AR-SCC ≡ Hom-closure (`atlas._leq`) route (both `True` on the fixture, both `False` on `kupisch=[3,2,2]`); the **P56 `is_simply_connected` agreement** on the ada/CC members — exercised on BOTH a `True` member (the tree fixture / commutative square, P56 = P61 = simply connected) AND a `False` member (the "square" `kQ/(a·b, c·d)`, P56 `is_simply_connected = False` == P61 Theorem-B `verdict = False`, `dim HH¹ = 1`), and where P56 returns `None` (Adian-Rabin) Theorem B resolves it — the headline; live since the P56 merge. **`qpa`** (`tests/qpa/test_recognizer_ladder_qpa.py`) — QPA 1.37 has NO quasi-tilted/shod/weakly-shod/laura/ada verb (a NamesGVars + `IsBoundGlobal` fail-if-appears probe); the gl.dim theorem gates corroborated via `GlobalDimensionOfAlgebra` (kA₃ = 1, rad²=0 A₅ = 4) and the pointwise id datum behind the shod refusal via `InjDimensionOfModule(S₃) = 2` |
| `families/` (catalog, zoo; Plan-29 trivial-extension/incidence batteries; Plan-31 certified trivial-extension presentation, `test_trivial_extension_presented.py`) | 166 | deep | closed-form family pins; zoo diversity gates; citations; Plan-31 special-case + Cartan + iso-invariance + CS≡bar pins |
| `strings/` (Plan-46 C5 gentle/string subsystem: reduced walks + σ/ε signs + string census + band detection; string/band module materialisation; string-τ by hooks/cohooks; the Avella-Alaminos–Geiss derived invariant; the `strings` block) — `tests/modules/test_strings_*.py`, `tests/invariants/test_ag_invariant.py` | 33 | deep + fast | Butler-Ringel `n(n+1)/2` interval count + Kronecker band existence (`oracle_literature`); string-τ ≡ engine τ + census count ≡ `knit_ar_quiver` vertex count (`oracle_crossengine`); `check_module` on every materialised string/band + permitted/forbidden thread partition of `Q_1` (`oracle_selfcert`); AAG-2008 pins reproduced verbatim (Nakaoka `arXiv:1811.00775` Example 2.15 = `{(3,2),(2,2),(0,3)}`) |
| `families/brauer.py` (Plan-46 Brauer graph algebra constructor from a ribbon graph + multiplicities) — `tests/families/test_brauer.py` | 10 | deep | `dim = Σ_v m_v·val(v)²` per-instance certificate + `is_symmetric` (`oracle_literature`); Brauer-star ≡ symmetric Nakayama `NakayamaAlgebra(n, mn+1, cyclic=True)` byte-equal Cartan (`oracle_crossengine`) |
| `modules/string_homological.py` (Plan-59 R34 — the homological string-algebra test `homological_string_test`, a three-valued semi-decision over the knitted indecomposables) — `tests/modules/test_string_homological_p59.py`, `test_string_homological_battery_p59.py`, `tests/qpa/test_string_homological_qpa.py` | 18 | deep + qpa | **`oracle_literature`**: `kD₄` (subspace) is not a string algebra with the live-verified 3-summand AR middle `E = P₁⊕P₂⊕P₃` (dim vector `{0:3,1:1,2:1,3:1}`), `kAₙ` never refutes (Suárez-Álvarez 2023). **`oracle_crossengine`** (the discriminating battery): a returned `not_string` (a ≥3-summand middle, k̄-sound) always has `is_string` False, and the refute side RAISES iff a ≥3 witness meets `is_string` True (one engine is a bug); the string side (kAₙ / gentle) never yields a false `not_string`, the non-string Dynkin side (kD₄/kD₅) always finds a witness. **`oracle_selfcert`**: every realized middle self-certifies (`baer_extension` cocycle + `ShortExactSequence` exactness) and the witness is a genuine ≥3 decomposition; rep-infinite / self-injective / presentation-less input refuses loudly. **`qpa`** (`tests/qpa/test_string_homological_qpa.py`): direct-session `IsSpecialBiserialAlgebra` parity on the R34 zoo + a standing `IsBoundGlobal` guard that FAILS if QPA ever ships a homological-string surface (it has none today) |
| `families/toupie.py` (Plan-59 R35 — `ToupieAlgebra` presented `kQ/I` + `is_toupie` recognizer + branch/direct-arrow counts + HH oracle + char-0 `sl_a` lower bound) — `tests/families/test_toupie_p59.py`, `test_toupie_hh_p59.py` | 23 | deep | **`oracle_literature`**: the `a`-Kronecker `HH^• = [1, a²−1, 0, …]` recomputed over QQ (ALS 2020, char-independent) for `a ∈ {2,3,4}`; `dim sl_a = a²−1` with `a` = # direct source→sink arrows (char 0, ALS Thm 6.5). **`oracle_crossengine`**: toupie HH bar ≡ CS degreewise where bar survives (incl. the non-hereditary commutative `[2,2]` at `top=2`, full CS `[1,0,0,0,0]`); `is_toupie` ≡ hand quiver shape. **`oracle_selfcert`**: the relation-free dimension certificate `dim = 2 + Σ(lᵢ−1) + Σ lᵢ(lᵢ+1)/2`, the connected + acyclic + degree-pattern recognizer (rejects `kD₄` and the path ⊔ oriented-cycle trap), `branch_count = out_deg(source) == in_deg(sink)`, and loud refusals (bad branches, length-1 cross relation, presentation-less, off-char-0 `sl_a`) |
| `modules/tilted.py` (Plan-60 R17 — the tilted-algebra recognizer `tilted_check` / `Algebra.is_tilted`: the Liu–Skowroński faithful-section criterion — a faithful section Σ with `Hom_A(X, τY) = 0` — certified by Ringel's slice theorem (Thm 1.9(2): `⊕Σ` a tilting A-module with `End_A(⊕Σ)` presented hereditary, then `A = End_H(D(S))`), with three theorem gates (hereditary ⇒ tilted; non-semisimple self-injective ⇒ not; `gl.dim > 2` ⇒ not) and the rep-finite exhaustive transversal search on the knitted `Γ_A`; the hereditary type from the section graph cross-checked with the recovered Gabriel quiver of `H`) — `tests/modules/test_tilted.py`, `test_tilted_reconstruct.py`, `test_tilted_oracles.py`, `test_tilted_p55_supports.py`, `tests/qpa/test_tilted_qpa.py` | 25 + 2 qpa | deep + qpa | **`oracle_literature`** — hereditary `kAₙ`/`kD₄` ⇒ tilted with the projective slice + Dynkin type; **rep-INFINITE hereditary (the Kronecker quiver `~A₁`, the acyclic non-commutative square `~A₃`) ⇒ tilted with the Euclidean type and the slice OMITTED, decided INSTANTLY from the Dynkin type of `A`'s own quiver — the AR knit is never built (bounded-runtime regression fence against the pre-fix hang)**; **`kA₃/rad²` ⇒ tilted, slice `{S₂,P₂,P₃}`, `End_A(S) = kA₃`, type `A₃`** (verified live); **`kZ₃/J²` (= cluster-tilted `A₃` = `Jac(3-cycle, αβγ)`) self-injective ⇒ not tilted** via Gate S (the knit REFUSES self-injective — the verdict comes from the theorem gate, a REAL not_tilted); `rad²=0 A₅` (gl.dim 4) ⇒ not tilted; **Liu's exact cut `P_b—S_b—P_a` of `kZ₃/J²` is sincere with `Hom(X, τY)=0` yet NOT faithful** — a STANDALONE `_is_faithful` demonstration (returns False; faithfulness cannot be weakened to sincereness — arXiv:1409.2054 Ex. after Thm 2.6), whereas the algebra `kZ₃/J²` itself is refuted by Gate S (self-injective), not by this faithfulness prune; the P55 fence flip — each `A_λ`/`A_ρ` support component of the ACLV-2.2(b) `rad²=0 A₅` is a tilted `kA₃/rad²` factor (products-of-tilted, ACLV Thm A). **`oracle_selfcert`** — Ringel Thm 1.9(2) two-way consistency (`⊕Σ` tilting AND presented `End_A(⊕Σ)` hereditary; `|Σ| = n`; the tilting-but-not-slice discriminator `A_A` of `kA₃/rad²`, a tilting module whose `End` has a relation, so NOT a slice); the theorem gates (self-injective `k[x]/(x³)` ⇒ not; the `gl.dim ≥ 3` lower-bound gate); the direct `_faithful_section_search` returns `None` on the rep-finite `rad²=0 A₅`; a **disconnected** `A` (`kA₂ ⊕ kA₂`) ⇒ not tilted (`reason="disconnected"`, up front) — tilted algebras are connected by definition (ASS2006), and the P55 support surface feeds the recognizer its connected components individually; the loud presentation-less refusal. **`oracle_crossengine`** — the reconstruction round-trip (`presented_form(End_A(S))` Gabriel type ≡ the section-graph type, loud on divergence; the checkable dim A / dim H invariants) and the gl.dim gate ≡ the direct search agreeing on `rad²=0 A₅`. **`qpa`** (`tests/qpa/test_tilted_qpa.py`) — QPA has NO tilted-algebra recognizer (fail-if-appears `IsBoundGlobal` probe); the slice module IS confirmed a tilting module by QPA's computational `TiltingModule(S, 1)` |
| `skewgentle/` (Plan-68 R32 — the triple `(Q, I, Sp)` + `is_skew_gentle_triple` recognizer + `associated_gentle`; the char-free idempotent-split constructor `SkewGentleAlgebra` (Chen §3); special-string re-gluing `classify`/`skew_gentle_module`/`skew_gentle_indecomposables`; the `is_representation_finite`/`brick_finite_certificate`/`support_tau_tilting` certificate; the `skew_gentle` block) — `tests/families/test_skew_gentle_{triple,split}.py`, `tests/modules/test_skew_gentle_{modules,certificate,tau_tilting}.py`, `tests/qpa/test_skew_gentle_qpa.py` | 38 | deep + qpa | **`oracle_literature`**: `Sp = ∅` byte-reduces to the plain gentle algebra (dim + Cartan + arrows); the mesh example `1→a→2→b→3, Sp={2}` where gentleness FORCES `a·b ∈ I` (`dim A^g = 9` vs the zero-split's 8); Chen Cor 1.2(c) **indecomposable** selfinjective-iff; the geometric rank `|R| = |Q₀|+|Sp|` (HZZ §6, `k×k → 4`, the 7-vertex `= 8`); Garcia–Lavoué rep-finite/rep-infinite verdicts. **`oracle_crossengine`**: the loop-free string census is a documented STRICT SUBSET of the AR-quiver indecomposables (headline `1→2, Sp={2}`: `5` of `6`; mesh `1→a→2→b→3, Sp={2}`: `8` of `11`) with EVERY census module embedding in the AR list up to iso — two independent enumerations (string materialisation vs P41 AR knit) agree on the overlap, and the SUFFICIENCY oracle for the mesh relations is the full AR vertex count (`11`) + `46` τ-tilting pairs (both engine-self-derived, VERIFIED live — no external number exists); the geometric rank ≡ engine pair rank; the band-census route is a ONE-SIDED no-false-negative cross-check on rep-INFINITENESS (`find_bands(A^g)` may miss special bands — W4 — so the exchange-graph route is authoritative, NOT a full agreement oracle). **`oracle_selfcert`**: the dim law `dim(split) == dim(associated gentle)` (HZZ Lemma 1.5, a **NECESSARY** check); `n_split = |Q₀|+|Sp|`; the char-free split over QQ/GF(p)/GF(2); the H1 all-string vertex-label regression (`exchange_graph(split)` does not raise on `Q=[1,2], Sp={2}`); every materialised module passes `check_module` + is indecomposable, a named string materialises to its pinned dimension vector (`M(a⁻¹)` → `{1:1, 2⁺:1, 2⁻:0}`), a type-`p` string has two distinct forms; loud refusals (invalid triple, special vertex already looped, length-≠2 relation, dim-certificate failure). **`qpa`** (`tests/qpa/test_skew_gentle_qpa.py`): the split algebra's `Dimension`/`IsSelfinjectiveAlgebra` parity as a plain `kQ/I` (QPA has NO skew-gentle surface — the standing `IsBoundGlobal` guard FAILS if that changes) |
| `surfaces/` (Plan-48 marked surfaces → ideal triangulations → gentle Jacobian algebras — `marked.py` `MarkedSurface`, `triangulation.py` `Triangulation` + fan/annulus/hexagon/once-punctured-torus, `qp.py` `quiver_of`/`potential_of`/`jacobian_of`, `flip.py` `flip`/`certify_flip_mutation`, `block.py` `surface_block`) — `tests/families/test_surfaces_*.py`, `tests/invariants/test_surfaces_arccount.py`, `tests/qpa/test_surfaces_qpa.py` | 70 | deep + fast + qpa | **`oracle_literature`**: the derived FST arc count `n = 6g−6+3(b+p)+Σkᵢ` on disc(n+3)→Aₙ / annulus(n,m)→n+m / once-punctured-torus→3, the FST admissibility exclusion list (monogon/digon/triangle, spheres with ≤3 punctures, once-punctured monogon), the hexagon-internal-triangle Jacobian `dim 6` (P44's pin), the disc-fan `kAₙ`, and the small annulus `C(2,1)` acyclic affine-`Ã₂`; **`oracle_crossengine`** (the P44+P46+P48 cross-subsystem tie): the disc-fan-`Aₙ` **orientation arbiter** (`quiver_of(fan((n+3))) = 1→2→…→n` exactly), `is_gentle(jacobian_of(T))` True across the disc/annulus/hexagon zoo (ABCP/LFS), flip ≡ Fomin–Zelevinsky matrix mutation on **every** interior arc, and `surface_block` AG invariant ≡ `strings.ag.ag_invariant`; **`oracle_selfcert`**: the two side-counting identities `3t=2n+c` / `p−n+t=χ`, the arc-adjacency + arc-count self-cert on every constructor (interior arcs in 2 triangles, boundary segments in 1), the self-folded refusal, and flip involution + `μₖ∘μₖ=id`; **`qpa`** (`tests/qpa/test_surfaces_qpa.py`): `IsGentleAlgebra`/`IsSpecialBiserialAlgebra` parity on the surface Jacobians + the standing `IsBoundGlobal` guard that FAILS if QPA ever ships a surface/triangulation constructor |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity; Plan-29 Tor; Plan-30 Krull–Schmidt decomposition; the retained injective-coresolution differentials certified exact; the Plan-35 wave-3a explicit Ext/Tor representatives — `complex_reps.py`; the Plan-35 wave-3c Yoneda exact sequences — `yoneda.py`; the Plan-37 C1 categorical glue — `morphism.py` first-class `ModuleHom` + kernel/image/cokernel, `ses.py` short exact sequences + split test + pushout/pullback, `endomorphism.py` `End(M)` as an Algebra, `direct_sum`/`is_direct_summand`, and covers/envelopes + radical/socle series + composition factors on `Module`; the Plan-40 C6 homological-dimensions family — `homdims.py`: public `syzygy`/`cosyzygy` (byte-stable extraction from `minimal_resolution`), the Igusa–Todorov φ/ψ on the finite K₀, dominant + Gorenstein dimensions, Ω/τ-periodicity certificates, and finitistic-dimension bounds; the Plan-44 C7 slice — `approximations.py` minimal left/right add(M)-approximations and `tilting.py` `is_tilting_module`/`is_cotilting_module` + self-certified `bongartz_completion`) | 341 | deep | AR/duality literature pins (ASS2006); **the Plan-37 categorical-glue self-certification** (`ModuleHom` validates the intertwining relations at construction; kernel/image/cokernel certified by rank-nullity + the epi–mono factorization `f = epi∘mono` + `f∘iota = 0` = `proj∘f`; SES exactness = the rank identity `im f = ker g`; split ⇔ a section solves; pushout/pullback squares certified by their universal-square identities; `End(M)` self-certified by `from_structure_constants(check=True)` with the regular-module `End(A_A) ≅ A` Loewy oracle; biproduct identities `proj_i∘incl_i = id`, `Σ incl_i∘proj_i = id`); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks; **the explicit Ext/Tor self-certification** (every shipped class satisfies `δ·v = 0` (Ext cocycle) / `d·v = 0` (Tor cycle) from its shipped differential; hand-checked kA₂ `Ext¹(S₁,S₂)` + loop `Tor₀ = M ⊗ N` cokernel labels; rep-count ≡ engine dims) + **the Yoneda exact-sequence self-certification** (every `Ext^n(M,N)` class is CONSTRUCTED as an `n`-fold exact sequence `0 → N → Q → … → M → 0` — the pushout middle module + connecting maps — and its exactness is self-certified at every joint: each map an `A`-module map, ends injective/surjective, `im = ker` by rank; the kA₂ Baer pin `0 → S₂ → P₁ → S₁ → 0` verified by the library's OWN `is_isomorphic`/`identify_standard`; a non-cocycle is refused loudly; the multi-vertex `_tor_boundary` collapse pinned on a rad²=0 Nakayama) + **the Plan-40 homological-dimensions oracles**: `oracle_selfcert` — the φ=ψ=pd identity for finite projective dimension, the Ω/τ-periodicity `is_isomorphic` certificates, and the decompose char-caveat propagation; `oracle_literature` — the Barrios–Mata truncated self-injective φ=ψ=0 closed form + projective additivity, the hereditary/self-injective dominant & Gorenstein values, and the cyclic-Nakayama period-from-Kupisch pins; `qpa` (`tests/qpa/test_homdims_qpa.py`) — live `GlobalDimensionOfAlgebra` / `DominantDimensionOfAlgebra` / `GorensteinDimensionOfAlgebra` agreement over kA₂ / kA₃(ab) / `line_abc_cde` / k[x]/(x³) (int or GAP `infinity`↔our infinite/unresolved marker) |
| `families/` (catalog, zoo; Plan-29 trivial-extension/incidence batteries; Plan-31 certified trivial-extension presentation, `test_trivial_extension_presented.py`; the Plan-44 C7 constructions — `one_point.py` `OnePointExtension`, `repetitive.py` `repetitive_slice`, `jacobian.py` `Potential`/`JacobianAlgebra`/`cyclic_derivative`, and the Task-C Gabriel-recovery battery `test_gabriel_recovery.py` over `core/basic.py`) | 187 | deep | closed-form family pins; zoo diversity gates; citations; Plan-31 special-case + Cartan + iso-invariance + CS≡bar pins; **the Plan-44 construction oracles** — `oracle_literature`: the one-point Cartan block `[[1, dim-vector M],[0,C_A]]` + `pd(S_ω)=pd_A(M)+1`, `repetitive_slice` `copies=1==A` + the `dim==(2·copies−1)·dim A` slice certificate, the hand-derived Jacobian triangle `dim=6`, and the `M₂(k)→k` / `kA₂` Gabriel round-trip; `oracle_crossengine`: `Jac(3-cycle, abc) ≅ cyclic Nakayama kZ₃/J²` (dim + Cartan) and `presented_form(End(⊕P_v)) ` recovers `kA₃` (tied to P37 `regular_corner_dims`); `oracle_selfcert`: complete-orthogonal primitive idempotents + the per-instance dimension/multiplicativity recovery certificate, the cyclic-derivative identities, and the loud char/split + `NotFiniteDimensionalError` refusals |
| `modules/quasihereditary.py` + `modules/recollement.py` (Plan-47 C-series: standard/costandard modules Δ(i)/∇(i), the quasi-heredity test `QHReport`, good-filtration multiplicities + BGG reciprocity, the characteristic tilting module + Ringel dual; `Recollement(A,S)` — the corner algebra `eAe`, the quotient `A/AeA`, and the six functors) — `tests/modules/test_quasihereditary_*.py`, `tests/modules/test_recollement.py` | 42 | deep | **`oracle_literature`** — Dlab–Ringel: `kA_n` natural-order Δ(i)=S_i / ∇(i)=[1..i], opposite-order Δ(i)=P(i); `kA_n` quasi-hereditary for BOTH orders; `k[x]/(x²)` NOT quasi-hereditary (loud note); Ringel: T = D(A) (natural) / T = A (opposite); the double-Ringel-dual Cartan Smith-form identity; **the discriminating NON-hereditary oracle** — the commutative square `1→2→4, 1→3→4, ab=cd` (gl.dim 2, every Δ simple) has T = D(A) (dim 9, dim-vector {1:4,2:2,3:2,4:1}), which the classical pd≤1 tilting certificate rejects and the single-pass `T(j)`-extension missed; the **`eAe`-vs-subquiver trap** (`kA₃`, S={1,3} → `eAe = kA₂` dim 3, NOT the subquiver `k×k` dim 2) + the worked `kA₃`, S={2} recollement. **`oracle_crossengine`** — BGG reciprocity `(P(i):Δ(j)) = [∇(j):S(i)]` across both orders. **`oracle_selfcert`** — top Δ(i)=S_i + [Δ(i):S(i)]=1 + socle ∇(i)=S_i; the greedy Δ-peel certificate (loud `certified=False` when no filtration); the characteristic-tilting arbiter — the genuine Ringel iteration lands in `F(Δ)∩F(∇)` (`Ext¹(Δ(j),T)=0` AND `Ext¹(T,∇(j))=0` for all j) with `is_tilting_module(T, n=gl.dim A)` (the exact global dimension, not the classical 1; gl.dim non-exact refuses loudly); the six-functor **adjunction dim identities** `dim Hom_A(j_!X,M)=dim Hom_{eAe}(X,j^*M)` etc., the counit isos `j^*j_!≅id`/`j^*j_*≅id`, and — through the ACTUAL functor outputs (each natural map a `ModuleHom` re-certified as an A-map) — BOTH BBD exact sequences: the counit `j_!j^*M → M → i_*i^*M → 0` (`im(counit)=ker(unit)`, unit epi) and the unit `0 → i_*i^!M → M → j_*j^*M` (mono, `im=ker`), over `kA₃`/the commutative square/`kA₅`/`kA₄`; the degenerate-`S` refusals (`S`=all → e=1, `S`=∅ → e=0); a `GF(2)` cell in each battery (Δ/∇/qh/recollement are char-clean). **QPA white space** — QPA has NO quasi-hereditary / recollement surface, so there is no `qpa` row here (stated in Honest scope) |
| `tautilting/` (Plan-45 C4 τ-tilting engine — `rigid.py` g-vectors + τ-rigidity, `pairs.py` certified support τ-tilting pairs, `mutation.py` the 2-term silting exchange + `exchange_graph` BFS (`_twoterm.py` = the K^b(proj A) cone/cocone + minimal-complex reduction engine), `torsion.py` the torsion lattice + Hasse orientation + bricks/semibricks, `stability.py` King θ-stability + the wall-and-chamber fan, `green.py` maximal green sequences, `silting.py` the 2-term silting bridge, `block.py` the algebra-level payload) — `tests/modules/test_tau_tilting_*.py` | 34 | deep | **`oracle_literature`** — `#sτ-tilt(kA_n) = Catalan(n+1)` (2/5/14), exchange-graph n-regularity, the AIR **four-way count identity** `#sτ-tilt = #f.f. torsion = #2-term silting = #semibricks` on kA₂/kA₃, hereditary `τ-rigid ⇔ rigid`, kA₂ = 2 maximal green sequences, and the **non-thin gate** kZ₂/rad² (symmetric Nakayama, dim 4): FOUR pairwise-non-isomorphic bricks — S₁, S₂ and the two projective-injectives P₁, P₂ *both* of dim-vector (1,1) — six semibricks, and the four-way identity 6 = 6 = 6 = 6 (the identity a dim-vector-keyed brick count silently broke: bricks → 3, semibricks → 5); **`oracle_selfcert`** — `g^{P_v}=e_v` + additivity, the four-axiom pair certification, mutation is an **involution** swapping exactly one g-column, every chamber g-matrix unimodular (det ±1), the n=2 fan **tiles R²** by an exact angular sweep (no atan2), the **n=3 L1/octahedron unfolding sanity** (nondegenerate net triangles, one per chamber — a rendering check, NOT a 3D tiling certificate), the **iso-class wall labelling** (kZ₂/rad²'s two same-dim-vector (1,1) walls carry distinct P₁ vs P₂ labels — four brick iso-classes across the walls, never a first-dim-vector collapse), King θ-stability on the worked kA₂ example, unique Hasse source/sink, brick `end_dim=1`, semibrick Hom-orthogonality, and the **honest semi-decision contract** (the 2-Kronecker is τ-tilting-infinite → loud `status="budget"`); **`oracle_crossengine`** — pair ↔ `Gen(M)` torsion-class injectivity, the fan's wall brick-normals ⊥ the shared g-facet (King), and the **no-proactive-char-guard** pin (GF(2) kZ₂/rad² agrees with the certified QQ counts — char ≤ dim computes where every module is a brick/splits). The webapp/GUI cross-runner τ-tilting tests are UNMARKED (extras-gated dirs, Plan-32 ruling). **QPA cannot compare — there is no `qpa` battery for τ-tilting** (see the honest-scope section). |
| `tautilting/congruence.py` (Plan-64 R26 lattice theory of torsion classes — the finite lattice `tors A` from the P45 oriented exchange graph, `Con(tors A)` via principal cover-congruences, the forcing order on bricks, canonical join representations, and the wide-subcategory poset via Enomoto's κ order = core label order; the algebra-level `congruences` no-code block) — `tests/modules/test_torsion_lattice.py`, `test_congruence_lattice.py`, `test_wide_subcategories.py` | 19 | deep | **`oracle_literature`** — the **kA₂ pentagon** (`tors kA₂` = N₅: 5 elements, semidistributive ✓, NOT modular, NOT distributive), `|tors kAₙ| = Catalan(n+1)` (5/14), `|Con(tors kA₂)| = 5` and the **kA₃ congruence lattice** `|Con(tors kA₃)| = 14`, the wide-subcategory counts (kA₂ → 5 ≅ **M₃**: 3 pairwise-incomparable atoms each covering ⊥ and covered by ⊤, rank histogram 1,3,1; kA₃ → 14 ≅ **NC(A₃)**: Whitney 1,6,6,1, 6 atoms, 6 coatoms, self-dual — Ingalls–Thomas `wide(kAₙ) ≅ NC(Aₙ)`), and the **τ-tilting-infinite honest refusal** (the 2-Kronecker emits NO lattice invariant); **`oracle_crossengine`** — **join-irreducibles ↔ bricks** (BCZ: lattice #join-irreducibles == P45 `bricks`, kA₂ → 3, kA₃ → 6), **canonical joins ↔ semibricks** (per-element down-cover label-sets == P45 `semibricks` as multisets, `#semibricks = |L|`), **#J(Con) = #bricks** (the forcing order lives on bricks, DIRRT); **`oracle_selfcert`** — `is_lattice` (unique join/meet per pair, loud `QuiverlabError` on violation — a DIRRT guarantee), `Con` distributive (Funayama–Nakayama, structural), the forcing-poset order-ideal count == `|Con|`, the **κ order == core label order** (Enomoto's theorem, asserted on every element — loud on mismatch), the `#J(Con) == #bricks` / `#join-irred == #bricks` internal gates, and the **kA₂ forcing DIRECTION** (the "V" `[(2,0),(2,1)]`: one brick forced-by two, ruling 4). Honest semi-decision STRICTER than a bounded region: certified complete iff `A` is τ-tilting-finite (DIJ); on a truncated exchange graph NO lattice invariant is emitted (`status="budget"`, all `None`). The webapp/GUI cross-runner `congruences` tests are UNMARKED (extras-gated dirs, Plan-32 ruling). **QPA cannot compare — no torsion-lattice/congruence/wide surface** (`tests/qpa/test_congruence_qpa.py` is an honest skip that FAILS if QPA ever ships one; the external cross-checks are the DIRRT/BCZ/Enomoto worked examples, not a live call). |
| `tautilting/cluster_morphism.py` (Plan-66 R29 — the τ-cluster morphism category W(A) (Buan–Marsh IMRN 2021) + its Hanson–Igusa cube-complex classifying space + the Igusa–Todorov–Weyman picture group, a thin combinatorial layer over the P45 exchange graph / P64 wide poset / P65 Jasso reduction: objects = τ-perpendicular wides (grouped from the g-fan faces by the bricks in `J(U)`), morphisms = support τ-rigid pairs of the source graded by rank, `face_vector` via the closed-star = C_W g-fan Jasso link identity; the `tau_cluster` no-code block) — `tests/modules/test_tau_cluster.py`, `test_picture_group.py`, `tests/qpa/test_tau_cluster_qpa.py` | 32 + 1 qpa | deep + qpa | **`oracle_literature`** — the **kA₂ picture group** (3 generators, 1 atom/pentagon relation, abelianization ℤ²), the **kA₃ picture group** (6 generators, **4 atom + 2 commutation**, abelianization **ℤ³** — H2, `6 − 3` distinct extension bricks, NOT ℤ²), the **kA₄ abelianization ℤ⁴** (H2 counterexample — the discarded `#bricks−#atom` gives the absurd ℤ⁰); the **classifying-space face vectors** kA₂ `(5,11,5)` / kA₃ `(14,49,49,14)` / kA₄ `(42,204,326,204,42)` with χ = **−1/0/2** (NOT `(−1)^n` — the discarded g-fan SPHERE's tautology); the **kA₃ vs kZ₃/rad² relation-split discriminator** (identical `#bricks=6 / #wide=14 / g_fan_face_vector=(1,9,21,14)` but **4+2 vs 3+3** relation split, hence differing classifying-space face vectors `(14,49,49,14)` vs `(14,48,48,14)`); the total morphism counts (kA₂ → 21, kA₃ → 126); the **Nakayama K(π,1)** verdict (kA₃/rad², kZ₃/rad² — Hanson–Igusa); the τ-tilting-infinite honest refusal (2-Kronecker emits NO category). **`oracle_crossengine`** — **object count = #wide** via TWO independent routes (P66 τ-perpendicular enumeration deduped by the `J(U)` brick set **vs** P64 `wide_subcategories(A).size`: kA₂→5, kA₃→14, kA₃/rad²→12, kZ₃/rad²→14); **#generators = #bricks** (vs P45 `bricks`, kA₂→3 … kA₄→10); **#relations = #rank-2 wides** (vs the category's rank-2 objects); the **morphisms out of `mod A` = Σ(g_fan_face_vector) = #support τ-rigid pairs** (kA₂→11, kA₃→45). **`oracle_selfcert`** — `object_count == wide_subcategories.size` internal gate (loud `QuiverlabError` on mismatch); **χ = Σ(−1)^k face_vector[k]** (the classifying space, NOT `(−1)^n`); **`face_vector[0] == object_count`** and **`sum(face_vector) == morphism_count`** and **`face_vector[-1] == #sτt`** (the discriminating H1 self-certs); the **factorization identity** `n!·#sτt == P65 signed_count` (a SELF-CERT consistency gate — `signed_count := n!·|eg.vertices|` is definitional in P65, fix-round ruling 3, NOT cross-engine); the composition rank-drop + unit spot-check; the atom relations name an extension brick (commutations do not); the **abelianization rank = #bricks − rank(relation matrix)** (SNF-decided, H2); `is_kpi1` theorem-anchoring (True only for Nakayama / hereditary Dynkin, else `None` — never asserted beyond the cited theorems). Honest semi-decision STRICTER than a bounded region: certified complete iff `A` is τ-tilting-finite (DIJ); on a non-complete exchange graph NO category invariant is emitted (`status`/`note`, all `None`) — `status != "complete"` from any graph → loud refusal (P65 M1). The webapp/GUI cross-runner `tau_cluster` tests are UNMARKED (extras-gated dirs, Plan-32 ruling). **QPA cannot compare — no cluster-morphism/picture-group/cube-complex/τ-perpendicular surface** (`tests/qpa/test_tau_cluster_qpa.py` is an honest skip that FAILS if QPA ever ships one; the external cross-checks are the ITW/HI/BM worked examples — the kA₂ pentagon, the Nakayama K(π,1) — not a live call). |
| `tautilting/wallchamber.py` (Plan-63 R25 — the wall `D(B)` of every brick as an EXACT rational inequality system over the submodule dim-vectors `D(B) = {θ : θ·dim B = 0 and θ·dim N ≤ 0 ∀ N ⊆ B}` (BST 2019 / King 1994), chambers = the g-vector cones of the support τ-tilting pairs, walls grouped one-per-brick ISO-CLASS via `torsion.bricks` + `_edge_brick`, the rank ≤ 3 drawing (n=2 exact extreme rays, n=3 L1/octahedron-projected grouped facets), `_is_green_path`; the `wall_chamber` compute kind) — `tests/modules/test_wall_chamber_*.py`, `tests/qpa/test_wall_chamber_qpa.py` | 27 + 1 qpa | deep + qpa | **`oracle_literature`** — kA₂ = 5 chambers / 3 walls with `D(S₁)`/`D(S₂)` full lines and **`D(P₁)` a single RAY** (direction (1,−1); the record's headline; Kaipel–Treffinger Ex. 13 verbatim wall equations); **kZ₂/rad² = 6 chambers / 4 walls** — the two dim-(1,1) bricks P₁, P₂ on OPPOSITE half-rays (1,−1)/(−1,1) of the ONE line θ₁+θ₂=0 (BST Rem. 3.19; self-injective coverage, M3); the **cyclic rad²-Nakayama N₃²** (1→2→3→1 mod rad²) = **14 chambers** (Kaipel–Treffinger Ex. 15, engine-verified against the EXACT algebra — an earlier draft mis-transcribed it as the LINEAR 1→2→3/rad², which is 12); the DIJ brick-finite ⟺ τ-tilting-finite gate (the 2-Kronecker is brick-infinite → a bounded region, `status="budget"`, NO count); kA₂ = 2 maximal green sequences, each a monotone source→sink chamber path; the **D₄ star** (hereditary underlying type D₄, {1→2, 3→1, 4→1}) = **50 chambers / 12 walls** (Cat(D₄)=50 support τ-tilting modules, 12 positive roots = 12 bricks) — RECOVERED as complete despite a spurious exchange-engine `mutate` error, NEVER the false "τ-tilting-infinite" truncation (the fix-round pin; asserts the OUTCOME so it outlives Plan 65's mutate root fix). **`oracle_crossengine`** — `#chambers = #support τ-tilting` (hereditary kA₂=5, kA₃=14; fan vs exchange graph); **linear kA₃/rad² = 12** (NOT KT, NOT Catalan — an honest consistency value); the inequality-defined `D(B)` extreme rays ≡ the grouped exchange-edge facet vectors **as sets of primitive directions** (two independent constructions — genuine set equality per brick, not merely non-empty). **`oracle_selfcert`** — `is_full_hyperplane`/`codim`, `θ·dim B = 0` on every wall ray + every submodule inequality `θ·dim N ≤ 0` holds on every ray (the D(B) membership, tying the system to the shipped `is_theta_semistable`); the facets ↔ exchange-edges **bijection** (each facet a real edge between two distinct chambers, no drop/dup/invention); every n=3 wall ray lies in D(B) and carries an L1 projection; **counts consistency** `#walls == #bricks == len(torsion.bricks(A))` on a thin (kA₂) and a NON-thin (kZ₂/rad²) example — the (1,1) dim-vector key has multiplicity 2, which a dim-vector-merged grouping would collapse; the **status="error" vs "budget" distinction** (fix-round) — the n-regularity completeness certificate (every discovered chamber has exactly n neighbours ⇒ the graph is closed under mutation), a non-n-regular error graph refused LOUDLY as unreliable (never a bounded sub-fan), and a genuine budget cap reporting `num_walls=None` with the discovered groups under `partial_wall_groups`; the **n=1 degenerate wall** (`D(B)={0}` is a point ⇒ `is_full_hyperplane=False`, no rays, `num_walls=0`). The webapp/GUI cross-runner `wall_chamber` tests are UNMARKED (extras-gated dirs, Plan-32 ruling). **QPA cannot compare — there is no `qpa` battery for wall-and-chamber** (`tests/qpa/test_wall_chamber_qpa.py` is an honest `NamesGVars()` probe that SKIPS and FAILS if QPA ever ships one; #chambers = #support τ-tilting inherits no live QPA cross-check either — FD-Applet / DIJ tables are the named external checks). |
| `derived/silting.py` (Plan-67 R30 **silting theory** on the P43 derived stack — the silting-object verifier `is_silting_object` (presilting = `Hom_{D^b}(T,T[n])=0` on the EXACT positive window `[1, n_max]`, weaker than tilting's `n≠0`; three-valued generation via the sufficiency ladder tilting/2-term-IJY/local/K0-basis-only/False), the co-t-structure dictionary `co_t_structure_of` (coheart `add(T)`, AI Prop 2.23(b)), single silting mutation `silting_mutate`/`silting_neighbors` (AI Def 2.34, one approximation triangle = the P45 `_twoterm` K^b(proj) minimal add-approx + cone/cocone + minimal-complex reduction, reused; the mutant is re-verified silting), the bounded-radius exploration `bounded_silting_exploration` (LOUD status, complete only for local), and `derived/block.py::silting_block` + the `silting` compute kind; **PLUS the Task-0 fix to `derived/tilting.py`** — `g_proj` computes the K₀ g-matrix in the PROJECTIVE basis `(+)_v Z[P_v]`, not the composition-factor basis, so `is_tilting_complex`'s generation leg is `det(g_proj)` not `det(Cartan·g_proj)` — a systematic P43 false-negative on non-unimodular Cartan; **PLUS the Plan-67 fix round (H1)** — `is_tilting_complex`'s generation is now honest three-valued (`generation` sibling field), reconciled with the silting K₀-only rung — see the honest-scope note) — `tests/modules/test_silting_verifier.py`, `test_silting_mutation.py`, `test_silting_exploration.py`, `test_silting_oracles.py`, `tests/qpa/test_silting_qpa.py` (+ the Task-0 regression AND the H1 three-valued-generation pins in `test_derived_tilting.py`) | 26 + 1 qpa | deep + qpa | **`oracle_selfcert`** — the Task-0 g_proj regression (`det g_proj = 1` for the regular object `A = ⊕P_v` over EVERY algebra, incl. the non-unimodular-Cartan AI Example 2.47 `det C = 0` and self-injective `kZ₃/J² det C = 2`); the presilting positive-window scan; the **silting-vs-tilting separation** (`P₂ ⊕ P₁[1]` over kA₂ is presilting + silting via the IJY 2-term rung but NOT tilting — nonzero `Hom_{D^b}(T,T[-1])`) and the **presilting-vs-silting separation** (`P₁` alone is presilting but `#summands ≠ #simples` ⇒ not silting); the **H1 three-valued-generation reconciliation** — a genuine width-2 rigid K₀-basis tilting complex over `k[1→2→3]/(ab)` returns `is_tilting == "unknown"` (generation `"k0_necessary_only"`, Rickard-open) AND `is_silting_object` returns `"unknown"` on the same input (the tilting rung and the K₀-only rung agree); every mutant re-verified **PRESILTING** (the decidable positive-window half, not full generation) + shares exactly `n-1` summand PROFILES (the fingerprint-ambiguity vs degeneration branches are distinct, loud); the `μ⁻∘μ⁺ = id` involution on kA₂ AND non-hereditary kZ₃/J² + a multi-vertex mutation chain; the `_cx_to_pc`/`_pc_to_cx` bridge round-trip on a non-2-term complex (degrees/dims/differentials byte-identical); `co_t_structure_of` refuses on non-silting input; the non-perfect-summand refusal. **`oracle_crossengine`** — the 2-term slice ≡ P45 τ-tilting is a REAL BIDIRECTIONAL check (H2): the SILTING engine INDEPENDENTLY enumerates the 2-term silting objects (BFS via `silting_mutate`+`is_silting_object`+a `g_proj` fingerprint in the canonical {0,1} window) and its SET of g-fingerprints equals P45's `exchange_graph` vertex set BOTH ways, the `silting_neighbors` edges equal the `exchange_graph` edges, pinned on kA₂ (5) AND the NON-hereditary self-injective kZ₃/J² (14); `End(μT)` underlying quiver ≡ Oppermann Thm 1.1 (Example 2.47 left mutation at vertex 1 → the single-arrow linear A₂, `dim End = 3`, corner-Cartan `[[1,0],[1,1]]` — the graded dg degrees are `# PIN`'d, the underlying (source,target) count is verified). **`oracle_literature`** — `k[x]/(x²)` silting = shifts (Thm 2.26, local-complete at radius 0); the kA₂ first mutation ring (Example 2.45); Example 2.47 cones (`X = cone(P₁→P₂)`); the local-complete vs kA₂-infinite exploration status. **QPA cannot compare — no `qpa` battery** (QPA 1.37 has no silting surface; `tests/qpa/test_silting_qpa.py` SKIPS and FAILS if one ever appears — see the honest-scope section). |
| `modules/exceptional.py` (Plan-65 R28 — the CLASSICAL hereditary exceptional-sequence surface: `is_exceptional_module` (the BRICK criterion `end_dim=1 ∧ Ext¹=0`), `is_exceptional_sequence` (backward Hom/Ext orthogonality), `braid_mutation` σ_i by the canonical universal-extension / kernel / cokernel constructions (four cases incl. Ringel's case-(d) two-step), `exceptional_sequences` (backward-orthogonality DFS + braid-orbit BFS transitivity + the Obaid closed form), `c_matrix`) — `tests/modules/test_exceptional_classical.py` | (see class table) | deep | **`oracle_literature`** — the Dynkin CES counts `A_n = (n+1)^{n-1}` for `n=2,3,4,5` (3/16/125/1296) and `D_4 = 162`, `count == closed_form_count == n!·h^n/|W|` (Obaid et al.), and `#exceptional == #indec` on every Dynkin case (rigid bricks); **`oracle_crossengine`** — the braid-orbit BFS orbit size `==` the direct-enumeration count `==` the closed form (two independent enumerations agree, the transitivity certificate); **the c-matrix rows `⊆` the P45 wall normals** (every dim-vector of a CES term is a brick labelling an exchange-graph edge — the hereditary c-vector reading overlaps the P45 wall-and-chamber surface, `A₂`/`A₃`); **`oracle_selfcert`** — every `braid_mutation` output is again a complete CES (orthogonality re-checked) and `σ_i∘σ_i⁻¹ = id` (incl. the worked `σ₁(P₁,S₁)=(S₂,P₁)` on kA₂, H4); the finite-universe `is_isomorphic` cross-identification of `L_XY`; every enumerated module is a rigid brick (`is_exceptional_module`); the rep-infinite hereditary refusal (2-Kronecker → `status="budget"`, instant Gabriel gate — no infinite knit). Scope: `is_exceptional_module` certifies the BRICK criterion = exceptional over algebraically closed `k` and on the Dynkin/QQ battery (M3, honest-scope (i)). **QPA** (`tests/qpa/test_exceptional_qpa.py`) — the exceptional-SEQUENCE surface has NO QPA verb (a fail-if-appears `NamesGVars()` guard), so the sequence enumerators are theory/self-cert-oracled; but the exceptional-MODULE criterion IS QPA-anchored per indecomposable on the kA₃ zoo — `dim End_A(M) = 1` via QPA `HomOverAlgebra` (the brick half) AND `dim Ext¹_A(M,M) = 0` via QPA `ExtOverAlgebra` (the rigid half), both agreeing with our engine. |
| `tautilting/exceptional.py` (Plan-65 R27 — the TAU-exceptional surface on the P45 engine: `tau_exceptional_objects` (indec τ-rigids + n shifted projectives), `tau_perpendicular_reduction` (the Jasso `C(U)`: sign-dispatched — τ-rigid `U` → the DIJ idempotent quotient `End(T_U)/⟨e_U⟩` built at the structure-constant level then presented; shifted `P_v[1]` → `A/⟨e_v⟩` with the tautological restriction), `is_tau_exceptional_sequence` (the recursive recognizer — tower or ambient), `tau_exceptional_sequences` (the ordered-sτ-tilt bijection `signed_count = n!·#sTt` + materialisation), `exceptional_sequences_block`) — `tests/modules/test_tau_exceptional.py` | (see class table) | deep | **`oracle_crossengine`** — the **materialised** τ-sequence count `len(sequences) == signed_count`, pairwise-distinct (per-level index disambiguates non-iso objects with equal dim-vector — `P₁`/`P₂` of `k(1↔2)/rad²`), each `is_tau_exceptional_sequence`, on `A₂`(10)/`A₃`(84)/the non-hereditary `k(1↔2)/rad²`(12) — the REAL check that `n!·#sτt` matches the enumerator (H1); the (a)↔(b) all-module-tower count `==` classical CES count on hereditary rep-finite (A₂=3, A₃=16, two independent enumerations); **the M-1 DEEP materialisation pins** — `kD₄` (subspace) signed `1200 = 4!·50` MATERIALISED with the all-module-tower count `162 = ` the classical CES count, and REORDER-INVARIANT (reversing the object enumeration at every level yields the same 1200 / 162), and `kA₄` signed `1008 = 4!·42` materialised pairwise-distinct where completion-multiplicity is real; **`oracle_selfcert`** — the formula identity `signed_count == n!·len(exchange_graph.vertices)` is SELF-CERT (tautological, NOT cross-engine — H1); `rk C(U) == n−|U|` **plus** the iso pin `C(S₁) ≅ kA₂` for `kA₃` (dim 3, connected — discriminates the correct completion from `k×k`, H5); **the H-1 completion-independence certificate** — `C(U)` as an ALGEBRA is completion-DEPENDENT (different tau-tilting completions give different dim/quiver, witnessed live on kD₄), but the enumeration invariant `#sTt(C(U))` is completion-INDEPENDENT (DIJ), self-certified across completions inside the reduction (loud raise on mismatch); **the H-2 tower chain check** — `_verify_tower` recomputes the reduction chain and rejects a forged tower (a hand-built `k×k` middle where the true reduction is `kA₂`: `#sTt` 4 ≠ 5), no longer trusting the rung labels; the shifted-projective reduction `C(P_v[1]) = A/⟨e_v⟩` (H2, incl. the recognized `(S₁,P₂[1])` sequence); the `status ∉ {"complete"}` loud refusal (M1 — the D₄ mixed star post-Task-0 completes to `signed=1200`); the τ-tilting-infinite loud refusal (Kronecker, instant). The webapp/GUI cross-runner `exceptional_sequences` tests are UNMARKED (extras-gated dirs). **QPA has NO τ-exceptional surface.** |
| `modules/homdims.py` + `modules/fractional_cy.py` (Plan-53 R23/R24 — φdim/ψdim as ALGEBRA invariants + the φ-spectrum & gaps + Lat-Igusa-Todorov finitistic certificates; the stable-category fractional Calabi–Yau dimension of self-injective algebras) — `tests/modules/test_phidim.py`, `test_phi_spectrum.py`, `test_lit_certificate.py`, `test_fractional_cy.py`, `tests/qpa/test_phidim_fcy_qpa.py` | 18 + 9 qpa | deep + qpa | **`oracle_selfcert`** — the standing chain `findim ≤ φdim ≤ ψdim ≤ gldim` on computed terms + the `gldim`-finite collapse (all four equal); the Barrios–Mata–Rama endpoints `{0, 1, m−1, m} ⊆ Spec_φ` when `0 < φdim = m < ∞`; the LIT bound sits in `[findim_lower, gl.dim]`; the LIT `None → 0` self-injective flip of `finitistic_dimension_bounds`; the fractional-CY `Ω^{m+ℓ} ν^ℓ ≅ id` self-identity + the loud not-self-injective refusal + the decompose char-caveat propagation; the family-4 `ψ_𝒟(V)+n+1` machinery arbitration gate (never below a known findim). **`oracle_literature`** — kA₃/J² φdim=ψdim=2, hereditary kA₃ φdim=1, self-injective k[x]/(x⁴) φdim=0; `k[x]/(x^a)` stable CY (1,1) + dual-numbers shift-trivial (0,1) + `Π(A₃)` 2-CY (2,1) (Ivanov–Volkov / Geiß–Leclerc–Schröer). **`oracle_crossengine`** — φdim via ⊕-of-all ≡ `max φ` over the knit indecomposables; the ℓ=1 fractional-CY certificate `S ≅ Σⁿ` ≡ the Ivanov–Volkov integer form `Ω^{n+1} ≅ ν⁻¹` on the simples. **`qpa`** (`tests/qpa/test_phidim_fcy_qpa.py`) — the prerequisites the values rest on: `IsSelfinjectiveAlgebra` ↔ `is_selfinjective` (the fractional-CY scope gate), `IsSymmetricAlgebra` ↔ `is_symmetric` (ν=id ⇔ symmetric, the k[x]/(x^a) derivation), and `NthSyzygy` of the simple ↔ `Module.syzygy()` (the `Ω²≅id` fact behind (1,1)); the φdim/CY VALUE probe SKIPS honestly (QPA has no such surface) and FAILS if QPA ever ships one. |
| `invariants/coverings.py` (Plan-56 R14/R16 — the presentation fundamental group π₁(Q,I) + π₁^ab by exact ℤ SNF, the Hom(π₁,k⁺)↪HH¹ Hurewicz bound, the three-valued `is_simply_connected`, the separation condition + `is_strongly_simply_connected` R16 recognizer that GATES P62, and `minimal_relation_counts` = the Tits-form rᵢⱼ; the `combinat/quiver.py` graph primitives) — `tests/coverings/*.py`, `tests/qpa/test_pi1_qpa.py` | 34 + 1 qpa | fast + qpa | **`oracle_literature`** — the CORRECTED commutative square (WITH the commutativity relation ⇒ π₁=1 / WITHOUT ⇒ ℤ), Le Meur Example 1 (monomial `⟨da⟩`⇒ℤ vs binomial `⟨da−dcb⟩`⇒0), trees/single-loop/multi-loop free ranks, the **two-independent-pairs discriminator** (a non-minimal generating set is NOT over-glued: π₁^ab=ℤ not 0 — the minimal relations are exact `I/(rad·I+I·rad)` linear algebra, never a greedy drop), the Zito example (simply connected but NOT strongly, separation fails at vertex 2), the star/deeper-tree/branching-source-in-a-non-tree separation discriminators (the corrected Qₐ = delete a + its transitive predecessor closure), the hereditary square genuinely NOT separated, and the loud intrinsic-π₁ refusal. **`oracle_crossengine`** — Hom(π₁,k⁺) ≤ dim HH¹ over the triangular zoo, equality on the Schurian square (the Hurewicz iso). **`oracle_selfcert`** — the graph primitives (spanning forest Betti number, components, induced subquiver), the annulus π₁^ab=ℤ (Betti of a 4-cycle), the `minimal_relation_counts` total == Σ block counts, and the char-caveat undecided path (decompose refusal ⇒ `undecided_char` ⇒ verdict None, never a mid-sweep raise). **`qpa`** (`tests/qpa/test_pi1_qpa.py`) — QPA 1.37 has NO fundamental-group / simple-connectivity surface (live `NamesGVars()` sweep); the probe SKIPS honestly and FAILS if QPA ever ships `FundamentalGroup`/`IsSimplyConnected`/`SeparationCondition`. |
| `invariants/tits.py` + `invariants/_tits_lists.py` (Plan-62 R19 — the COMBINATORIAL Tits form `q_A = Σx_i² − Σ_{arrows} + Σ r_ij x_i x_j` with `r_ij` = the P56 `minimal_relation_counts` (the `I/(rad·I+I·rad)` count, cross-checked vs `dim Ext²(S_i,S_j)`); `is_weakly_positive` (Ovsienko box-6, positive-definite / isotropic-radical / Euclidean-restriction fast certificates + connected-support branch-and-bound); `is_weakly_nonnegative` (positive-SEMIdefinite/Euclidean True certificate — with `HYPERCRITICAL_COVERAGE` empty by design, PSD is the ONLY route to a `True`/tame verdict — + the classified hypercritical list — K3/T334/T245/**T237 the 10-variable refutation** — as the primary False route + a sound witness-finder, honest `None` outside coverage, NEVER a guessed True); `tame_wild_certificate` — the rep-finite/tame/wild trichotomy gated on the P56 strong-simple-connectivity certificate over an algebraically closed base field (the CC working domain, or — absent P61's `is_algebraically_closed` flag — any characteristic-0 field read by base change to the algebraic closure; `GF(p)`/`GF(p^n)` refused), P56's `None` propagated) — `tests/invariants/test_tits_form.py`, `test_weak_positivity.py`, `test_tits_lists.py`, `test_tame_wild.py`, `test_tame_wild_oracles.py`, `tests/qpa/test_tame_wild_qpa.py` | 44 + 9 qpa | fast + qpa | **`oracle_literature`** — the m-Kronecker ladder `q(1,1)=2−m` (1/2/≥3 = weakly positive / Euclidean isotropic / wild), the Dynkin trees rep-finite, the Euclidean trees (`~A_n`/`~D_n`/`~E_{6,7,8}`) weakly nonnegative not weakly positive with `q(radical)=0`, `T_{2,3,7}` wild with the SINCERE 10-vertex defect `(12,6,8,4,10,9,7,6,4,2)`, `q=−1` (the classified-list refutation of the `≤9`-variable cap), the Bongartz + Brüstle–de la Peña–Skowroński theorem verdicts (Dynkin→rep-finite, Euclidean→tame, `T_{2,3,7}`→wild). **`oracle_crossengine`** — `r_ij` (P56 count) == `dim Ext²(S_i,S_j)`; the combinatorial Tits Gram == P38's Euler symmetrization `E+Eᵀ` where both are defined (`gl.dim ≤ 2`); the AR-knit `is_complete` ⇒ `rep-finite` + `weakly_positive` (module-category enumeration vs quadratic-form arithmetic). **`oracle_selfcert`** — every `False` witness is exact and `≥ 0` with `q(witness) ≤ 0`/`< 0`, the budget→`None` honesty (no guessed True), the m-Kronecker form-layer split (not simply connected ⇒ verdict None), the non-CC form-computed/verdict-refused split, P56's `None` propagation, the relabelling-invariant hypercritical matcher, the **two-level verdict field gate** (P61's `is_algebraically_closed` flag when present, else `characteristic == 0` — pinned to hold on both the pre-P61 base and post-merge dev), the **honest base-change wording** (over QQ/QQ(i) the reason never fabricates "over an algebraically closed field"), and the **Bongartz-certain rep-infinite** surfacing (`certified == "rep_infinite"` when the form is not weakly positive on a simply connected algebra, `rep_type` staying `None` because tame-vs-wild is genuinely undecided). **`qpa`** (`tests/qpa/test_tame_wild_qpa.py`) — QPA 1.37 SHIPS a Tits-form surface: `IsWeaklyPositiveUnitForm`/`IsWeaklyNonnegativeUnitForm` of `TitsUnitFormOfAlgebra(A)` agree with our decisions on A5/E8/`~E8`/`T_{2,3,7}`/K2/K3/kA₃-rel (a real two-implementation oracle, incl. the `T_{2,3,7}` list-decided `wnn=False` confirmed independently); QPA has NO representation-type (tame/wild) VERDICT verb, so the verdict layer stays theory-oracled (a fail-if-appears probe). |
| `invariants/hh1_lie.py` (Plan-70 R11 — HH¹ = Der(A)/Inn(A) as a Lie algebra: the field-general Leibniz-null-space `derivations` + `inner_derivations`, the matrix-commutator bracket structure constants in the `der_inn` basis, the derived / lower-central series + solvable / nilpotent / abelian / perfect verdicts over ANY exact Domain; the char-0 classification behind a hard `characteristic == 0` gate — the Killing form `κ`, the solvable radical `rad = [L,L]^⊥` (de Graaf), the Levi factor `S = L/rad`, the sl₂-count via the invariant-symmetric-forms simple-ideal count, toral rank + type; the RSS Ext-quiver solvability certificate; the Gerstenhaber-bracket cross-engine helper) — `tests/invariants/test_hh1_lie.py`, `tests/qpa/test_hh1_lie_qpa.py` | 39 + 3 qpa | fast + qpa | **`oracle_literature`** — the `k[x]/(x^n)` dichotomy: solvable with derived-series dims `[2,1,0]/[3,2,0]/[4,3,1,0]` for `char ∤ n`, and `W₁` (Jacobson–Witt) at `n = char = p` (dim `p`, perfect + simple for `p ≥ 3`, 2-dim solvable at `p = 2`), plus the `char | n, n ≠ p` refinement (`k[x]/(x⁴)`, `k[x]/(x⁶)` over `GF(2)` SOLVABLE, `k[x]/(x⁶)` over `GF(3)` NOT); `kK₂` (Kronecker) `≅ sl₂` (dim 3, perfect, simple, `levi_type "A1"`, sl₂-count 1, toral rank 1); `T(kK₂) = k ⋉ sl₂` (dim 4, radical 1, sl₂-count 1, NOT simple — the card refinement); the **default-field `CC` pin** (formally algebraically closed, exact QQ arithmetic → computes, equals the explicit-QQ run, `base_change_note` present); the **`GF(4) = GF(2²)` non-prime-field row**; the RSS no-loops/no-parallel ⇒ solvable pins (`k(1↔2)/rad²`, 3-cycle`/rad²` over QQ/GF(3)/GF(5)). **`oracle_crossengine`** — `dim(Der/Inn) == A.hochschild_cohomology(top=2).dims[1]` on the zoo; the Plan-35/P51 degree-(1,1) `gerstenhaber_brackets` route agrees with the Der/Inn commutator on BASIS-INDEPENDENT data (dim, solvable, nilpotent, derived-series dims) over GF(p) in-window; `rss_solvable_certificate ⇒ is_solvable_hh1`. **`oracle_selfcert`** — `Inn ⊆ Der` **checked computationally** (`[D, ad] ∈ span(Inn)` for every Der×Inn pair), the **Jacobi identity** on the HH¹ bracket structure constants across the pin zoo (incl. the char-p Witt case and sl₂ ⊕ sl₂), the char-p loud gate on radical/Levi/sl₂/toral (`require_char0=True` raises), the base-change note gated on the ARITHMETIC field not the formal `is_algebraically_closed` flag (present over CC AND QQ), the oversize-budget loud refusal, abelian/perfect flags, the presentation-less RSS refusal. **`qpa`** (`tests/qpa/test_hh1_lie_qpa.py`) — QPA has NO HH¹-Lie surface, but GAP's core Lie library (`LieAlgebraByStructureConstants` + `IsLieSolvable` / `LieDerivedSeries` / `SemiSimpleType`) recomputes our verdicts from the shipped structure constants (`sl₂ → "A1"`, `IsLieSolvable == false`; `k[x]/x³/QQ` solvable, derived dims `[2,1,0]`; the MULTI-FACTOR `sl₂ ⊕ sl₂` where GAP's space-separated `"A1 A1"` is compared to our `"A1+A1"` as a normalized sorted factor multiset — the path that feeds P71); an honest skip that FAILS if the Lie functions disappear. |
| `hochschild/lie_module.py` + `modules/decompose.py::decompose_representation` (Plan-71 R12 — HH•(A) as a graded Lie module over HH¹(A): the field-general action `ρ_n(D)` = the Gerstenhaber degree-1 Lie derivative `L_D f = D∘f − Σ f(…,Da_i,…)` on the normalized bar cochain complex over ANY exact Domain (no resolution engine), the module axiom `ρ_n([D,E]) = [ρ_n(D),ρ_n(E)]` + inner-acts-zero self-certs; over `char == 0` the from-scratch maximal-torus weight decomposition (Cartan subalgebra by de Graaf Engel/Fitting → ad-semisimple part → radical-toral extension → rational simultaneous diagonalization); the indecomposable Lie-module summands via the matrix-level Fitting core `decompose_representation` driven by the directly-computed commutant `End_{HH¹}(HH^n)`) — `tests/hochschild/test_lie_module.py`, `tests/qpa/test_hh_lie_module_qpa.py` | 32 + 2 qpa | fast + qpa | **`oracle_literature`** — `kK₂` (an ALS toupie, a=2) `HH¹ = sl₂` acting on `HH^1` as the IRREDUCIBLE adjoint `L(2)` (one summand dim 3, weights a symmetric sl₂-string `{−c,0,c}`, `c ≠ 0` — the PATTERN, not the literal integers; `HH^0 = L(0)`, one weight 0); the 3-Kronecker (ALS toupie a=3) `HH• = [1,8]`, `HH¹ = sl₃` irreducible adjoint (dim 8); `k[x]/(x^n)` the truncated-Witt anchor (`HH¹ = ⟨x^{i+1}∂⟩`, a Virasoro-subquotient *analogue*) with each `HH^m` a single indecomposable and an equal-gap grading-weight progression; the `k[x]/(x²)` weight ladder; `k[x,y]/(x,y)²` `HH^n` multi-summand `gl₂`-modules. **`oracle_crossengine`** — `hh_dims == A.hochschild_cohomology` on the zoo; the **STRONG entry-wise** check `L_D` on the engine degree-n cochain basis ≡ the shipped `engine.tt_calculus.gerstenhaber_bracket_cochain(·,1,n,·)` cochain-for-cochain over GF(p) (the `p=1` Koszul-sign collapse is exact, not just up to invariants); the weaker in-window `gerstenhaber_brackets` degree-`(1,n)` induced-map arbiter on BASIS-INDEPENDENT data (dim HH^n, dim associative envelope `B_n`, dim commutant `End`). **`oracle_selfcert`** — the module axiom `ρ_n([D,E]) = [ρ_n(D),ρ_n(E)]` and `ρ_n(ad_x) = 0` (the action factors through `HH¹ = Der/Inn`) on `k[x]/(x^n)`/`k[x,y]/(x,y)²`/`kK₂`; the torus generators commute; `Σ_λ dim HH^n_λ = dim HH^n`; each certified summand reassembles `HH^n`; the char-p LOUD gate on **weights/torus** (`require_char0=True` raises) while the SUMMAND decomposition is governed INDEPENDENTLY by decompose's own guard (`char 0 or char > d_n`) — so summands appear at `char p > d_n` while `weights` is `None`; `decompose_representation` reproduces `decompose(M)` on a standard module (byte-stable). **`qpa`** (`tests/qpa/test_hh_lie_module_qpa.py`) — QPA has NO Hochschild-Lie-module surface (an honest `NamesGVars` skip that FAILS if one appears); the substantive oracle is GAP's own MeatAxe (`GModuleByMats` + `MTX.CompositionFactors`) on the fed `ρ_n` matrices — the MULTI-summand `k[x,y]/(x,y)²` `gl₂` case (per-degree constituent dims `[1,3]`/`[2,4]` == ours), with the sl₂ = `kK₂` case a trivial irreducible sentinel. |
| `hochschild/deformations.py` (Plan-78 R13 — formal deformation theory: the infinitesimal space `HH²(A,A)`, the Gerstenhaber obstruction `α ↦ [α,α] ∈ HH³` on the native CS bracket, the order-by-order Maurer–Cartan report with the nilpotent-regime gate (MRRS Thm 5.4 `MC = Z²`), the presented deformed algebra `A_α` re-certified flat (`dim A_α == dim A`) with an Ext-algebra handoff, and the `B(A)[1]` L∞ companion with an UNCONDITIONAL `rad² = 0` dg-Lie certificate; the `deformations` compute kind across all three tiers) — `tests/hochschild/test_deformations{,_mc}.py`, `test_deformed_algebra.py`, `test_linfty_l3.py`, `tests/gui/test_deformations_runner_twin_p78.py`, `tests/webapp/test_deformations_kind_p78.py` | (see class table) | fast + deep | **`oracle_literature`** — the classical dictionary `HH² = infinitesimal deformations`, `HH³ = obstructions` (Gerstenhaber 1964); the commutative `k[x,y]/(x²,y²)` (= `QuantumCI(-1)` under this codebase's `xy + q·yx` convention) `HH^• = [4,4,5,6]`, so `HH² = 5` with a genuine OBSTRUCTED direction. **cross-runner contract (UNMARKED by the Plan-32 extras-gated ruling — `tests/gui` / `tests/webapp` never carry oracle-class markers, so the audited counts stay environment-independent)** — the wheel's `hpc.spec` core and the Pyodide `docs/gui/runner.py` emit the `deformations` block byte-for-byte identical (obstructed `QuantumCI(-1)` on the `deep` leg, the cheap `rad²=0` 3-cycle on `fast`). **`oracle_selfcert`** — `[α,α] = 0` scanned over EVERY `HH²` basis direction (never a diagonal-vanishing shortcut — a basis cocycle can self-obstruct, `[α₁,α₁] = 2 ≠ 0` live); the `rad² = 0` ⇒ dg-Lie certificate; flatness `dim A_α == dim A` on the rebuilt presentation; loud refusals (char-p interpretation gate, `A.dim > DEFORM_MAXDIM = 32`, unit / non-admissible direction, presentation-less input) as clean typed blocks, never a 500. |
| `modules/radical.py` + `modules/ar_invariants.py` (Plan-57 R37/R21 — the radical filtration of `mod A`: exact `dim rad^n(X,Y)` layers on the knitted indecomposables (route (ii), left-to-right composition), the nilpotency index of `rad(mod A)`, the `rad^∞ = 0 ⇔ representation-finite` (Auslander) gate; Liu's left/right degrees of irreducible maps (`fg ∈ rad^{m+2}`, corrected), sectional paths, the postprojective/preinjective/regular partition, directing modules, the representation-directed recognizer (`Γ_A` acyclic), the generalized-standard flag) — `tests/modules/test_radical_filtration.py`, `test_liu_degrees.py`, `test_ar_invariants.py`, `tests/qpa/test_radical_qpa.py` | 29 + 2 qpa | deep + qpa | **`oracle_selfcert`** — the descending filtration `rad^{n+1} ⊆ rad^n`, `rad^N = 0` while `rad^{N-1} ≠ 0`, composition closure, `rad^1 = Hom` off the diagonal (QQ↔GF(32003) parity); the Thm-1.3 `max_a{r_a+1}` index walk ≡ route-(ii) `nilpotency_index` (two readings of ONE knit); the **degree-vs-layer witness triple** externally asserted — a finite left degree returns `(Z, g)` with `g ∈ rad^d`, `g ∉ rad^{d+1}`, `g.then(f) ∈ rad^{d+2}` (the certified layer drop); the `d_r(f) ≡ d_l(Df)` opposite-algebra multiset symmetry; the partition totality; the honest semi-decision contract (self-injective ⇒ `status="unsupported"`, no verdict; budget ⇒ window-restricted layers, `nilpotency_index=None`, no `rad^∞` claim; never a false `representation_directed`). **`oracle_literature`** — `N(kA_n) = n` (Thm 1.5(a)), the `kA_3` layer table (total `dim rad^1 = 9`, `rad^2 = 3`, `N = 3`), `kA_3/J^2` index 3 (Thm 1.3 `max{2,3,2}`), the hand-derived `kA_2`/`kA_3` left+right degree tables by module name (mono ⇒ `d_l = ∞`, epi ⇒ `d_r = ∞` on directed `kA_n`), `kA_n`/`D_4` representation-directed with no regular modules and every indecomposable directing, sectional composites nonzero; the live-verified non-directed witness `NakayamaAlgebra(kupisch=[3,2,2])` (non-self-injective, 7 indecomposables, oriented `Γ_A` cycle) separating `representation_directed=False` from `generalized_standard=True`. **`oracle_crossengine`** — the genuinely independent route-(i) `ZA_n` mesh closed form (pure interval combinatorics on dimension vectors, no Hom/matmul) ≡ route-(ii) `dim rad^n(X,Y)` on `kA_3`/`kA_4`. **`qpa`** (`tests/qpa/test_radical_qpa.py`) — QPA has NO module-category-radical or degree surface (`RadicalOfModule` is the Jacobson radical of ONE module, not `rad(X,Y)`; no `LeftDegree`/`RightDegree`/`NilpotencyIndexOfRadical`): a fail-if-appears `IsBoundGlobal` probe, plus the one checkable slice `dim rad(X,Y) = dim Hom(X,Y)` (`X ≇ Y`) corroborated against `HomOverAlgebra`. |
| `invariants/han.py` + `families/arrow_removal_subalgebra` / `enveloping_algebra` (Plan-73 R7 — Han's conjecture transported across an arrow-removal subalgebra `B ⊆ A` (CLMS `2009.05017`): the subalgebra by arrow+relation extension (Def 5.2, certified `dim B` via `J ∩ B = 0`, P72-seam-gated `remove_arrows == subalgebra` on inert `F`), `is_tensor_nilpotent` (the length-index theorem Def 5.19 + J-interrupter route + honest capped semi-decision), the finite `pd_{B^e}(A/B)` PRIMARY route via `gl.dim B < ∞` (Ex 6.1, free — no enveloping algebra) with the `B ⊗ B^op` fallback, the Jacobi–Zariski relative bar complex `HH_*(A|B)` (Thm 2.2, finite via Cor 2.4), and `han_transport` — the FIVE-ROW injection/iso ladder; the `han_transport` compute kind) — `tests/invariants/test_han_*.py`, `tests/gui/test_han_runner_twin_p73.py` | (see class table) | fast + deep | **`oracle_literature`** — CLMS Ex 5.3 tensor-nilpotency index 2; Ex 5.5 NOT nilpotent (the `not_bounded` witness); Ex 5.4 `rad` is `E`-tensor-nilpotent iff `Q` acyclic (`E = kQ_0`); Ex 6.1 the `gl.dim B < ∞ ⇒ pd_{B^e}(A/B) < ∞` route. **`oracle_crossengine`** — a 3-WAY cross-check of the relative bar complex: this plan's `HH_*(A|B)` == `relative.py` == the fast engine on `kA₃` (`HH_0 = k^3`); the live `HH_*` agreement/divergence pair (Ex 5.3/5.5: agrees in every degree `[5,0,…]` vs diverges `[4,1,2,3,3,3,4,6,8,9] ≠ [3,0,…]`); the cross-runner `han_transport` twin. **`oracle_selfcert`** — the INJECTION BOUND `dim HH_m(B) ≤ dim HH_m(A)` asserted on every run, with EQUALITY asserted ONLY under the `bounded` row (a `bounded` verdict whose dims disagree raises loudly); the certified `dim B`; the enveloping fallback exercised on a `gl.dim B = ∞` loop extension. |
| `hochschild/split_extension.py` + `hochschild/arrow_removal.py` (Plan-72 R5/R6 — (a) the CMRS split-extension / trivial-extension Hochschild LES `0→M→L→B→0` in `HH^•(L,−)`: `Bimodule.inflate` (the `L`-bimodule by inflation along `π: L↠B`), `split_extension_cohomology`/`_homology` assembling `HH^•(T(B))` from the flanks + the SNAKE connecting map + the direct cross-check, `hh1_grading_witness` (CMRS Thm 5.5); (b) certified arrow removal/addition (CLMS `1812.07655`): `inert_arrows` (Def 3.1), `remove_arrows`/`add_arrows`, `arrow_removal` (the clean `HH_{≥2}` iso Thm 3.2 + the cohomology Ext-correction Thm 4.2)) — `tests/families/test_split_extension_les_p72.py` (26), `tests/families/test_arrow_removal_p72.py` (19) | 45 | deep | **`oracle_literature`** — the frozen `HH^•(T(kA₂))=[3,1,1,1,1]` / `T(kA₃)=[4,1,1,1,1]` / `T(kD₄)=[5,1,0,0,2]`, the flank pins (`T(kA₂)`: `flank_M=[2,1,0,1,2,1,0]`, `flank_B=[1,0,1,2,1,0]`; the `T(kD₄)` **boundary-rule** pin `flank_M=[4,1,0,0,1,4]` with `HH^5(L,M)=4≠0`), the `HH^1(T(B))≠0` sweep + the directed formula `HH^1(T(B))=1+HH^1(B)` (kA₂/kA₃/kD₄→1, 2-Kronecker→4, N/A on `k[x]/(x^a)`); the CLMS `HH_{≥2}(A)=HH_{≥2}(B)` iso on the P1–P4 removed-arrow pairs and the **P2 `n=2` Ext-correction `=1`** (the homology/cohomology asymmetry). **`oracle_crossengine`** — **assembled `HH^•(L)` == direct `HH^•(L)`** (the LES-over-the-CS-Hom-complex route vs. the standalone bar/CS engine; the `T(kD₄)` case EXERCISES the `top+1` M-flank boundary rule with `δ^4=2`); `arrow_removal` `HH_n(A)==HH_n(B)` for `n≥2` == direct on P1–P4; the leading-piece **identity `flank_B[0]==HH^0(B)`** + the **summand inequality `flank_B[n]≥HH^n(B)`**, tight/non-vacuous on `T(2-Kronecker)` (`flank_B[1]=3=HH^1(B)`, NOT an equality tautology — `D(B)` is bimodule-projective only iff `B` self-injective); the `Bimodule.inflate(D(B))` flank == the report's `flank_M` (inflate ≅ the ideal-`M` sub-bimodule). **`oracle_selfcert`** — the LES exactness / degreewise SES dim identity `dim C^n(L,M)+dim C^n(L,B)=dim C^n(L,L)`; the snake `δ` well-defined (the adversarial-lift test — shifting the lift by any `M`-cochain leaves `δ[φ]`'s class fixed); `Bimodule.inflate(...).check()` (incl. the multi-vertex `kD₄`); the `HH^1` grading cocycle is a cocycle with nonzero (outer) class; `dim T(B)=2·dim B`; `remove_arrows` dim drop + inert witnesses (incl. the P4 **binomial** relation `a*c−b*d`); the addition round-trip `add_arrows(remove_arrows(A,D),D)≅A`; `HH_0` provably invariant under inert removal on P1–P4; the loud refusals (non-inert arrow with witness, relative-cycle addition, unpresentable `L`). **QPA has NO split-extension-LES or arrow-removal HH surface** (see the honest-scope entry). |
| `modules/barcode.py` + `families/commutative_ladder.py` (Plan-69 R33 — the persistence/TDA bridge, representation theory FIRST: `barcode(M)` = the interval decomposition of an `A_n`/zigzag persistence module (Gabriel / Botnan–Crawley-Boevey — the support intervals of the Krull–Schmidt summands), field-robust over EVERY exact domain incl. GF(2) (interval modules are bricks, `dim End = 1`); the forward-only `mat_rank` rank-formula second route; `CommutativeLadder(n) = A_n □ A_2` (scalar `"i_j"` vertex names — tuple labels crash the AR knit) with the AR-quiver-indexed generalized persistence diagram for `n ≤ 4` (char-scoped: char 0 / char > dim), rep-finite iff `n ≤ 4` (Escolar–Hiraoka); the `barcode` module-side compute kind) — `tests/modules/test_barcode_p69.py`, `test_barcode_cl_p69.py`, `test_barcode_battery_p69.py`, `tests/families/test_commutative_ladder_p69.py`, `tests/qpa/test_barcode_qpa.py` | (see class table) | deep + qpa | **`oracle_literature`** — the `A_5` filtration barcode `{[1,5],[2,2],[4,4]}` (hand-derived + engine-confirmed, identical over QQ and GF(2)); the Escolar–Hiraoka **rep-finiteness BOOLEAN** `CL(n)` rep-finite iff `n ≤ 4` (theorem, literature-confirmed VERBATIM — the `n=2,3` knits terminate `"complete"` AND `CL(5)` is refused at construction; `CL(4)` is theorem-certified + dim-30-built, its knit excluded from CI as impractically slow — see honest scope); **the `CL(2)=11` / `CL(3)=29` indecomposable counts** — figure-confirmed against the Escolar–Hiraoka AR-quiver figures (Fig. 13 = `3+5+3=11` all thin; Fig. 14 = `1+6+11+6+5=29`, exactly 2 non-thin, matching the engine's 27+2 thin/non-thin split; the equioriented ladder = incidence algebra of the `[n]×[2]` grid poset, iso-invariant), the plan's BLOCKING reconciliation (Task 6 Step 1a) DONE — kept ALSO `oracle_selfcert` (the same knit reproduces them, Plan-32 overlap). **`oracle_crossengine`** — `decompose` barcode ≡ the `mat_rank`-only rank formula on forward lines; `GF(2) ≡ QQ` field-robustness (bricks); the CL diagram's distinct-indecomposable count ⊆ the knit's vertex count; **QPA** `DecomposeModuleWithMultiplicities` parity on `A_5` (barcode's interval bars == QPA's summand dim-vector multiset) + the CL(3) diagram DECOMPOSE closer (our decompose of a CL(3) module ≡ QPA over GF(23)). **`oracle_selfcert`** — the interval-sum identity `Σ m_i·dimvec == M.dimension_vector()`; every summand's support is a contiguous interval in the line order; the CL AR-index match (`is_isomorphic` certificate, an unmatched summand a loud bug); the `essential` forward-only flag (always `False` for zigzag — M2; and always `False` for a commutative-ladder interval bar as well — a ladder base is not a single monotone persistence line, and the CL primary object is the AR-indexed diagram, not the bars — a deliberate simplification of the plan's Bar-schema prose); the `n=1` single-vertex handling (one bar `[1,1]·dim M` — M1); the char-scope refusal CHAIN over GF(2 ≤ dim) (the knit CATCHES the internal `decompose` char-refusal → `status="error"`, `is_complete=False` → `barcode` refuses, never a silent truncated `"complete"` — M6); the loud refusals (non-`A_n`/non-CL, presentation-less); the forward dimension certificate `3·C(n+1,2)` (9/18/30) and the scalar-vertex + graph-shape recognizer; **the `CL(2)=11` / `CL(3)=29` indecomposable counts** (the SAME AR knitter that produces them — kept `oracle_selfcert` in addition to the figure-confirmed `oracle_literature` pin, Plan-32 overlap; the thin/non-thin 27+2 split reproduces Fig. 14). **QPA has NO persistence/barcode/commutative-ladder surface** (`tests/qpa/test_barcode_qpa.py` — a fail-if-appears `IsBoundGlobal` guard on `PersistenceDiagram`/`CommutativeLadder`/`Barcode`); the whole-AR-quiver 29-vertex COUNT is not a QPA verb (QPA has no AR-quiver enumeration), so it stays `oracle_selfcert`. |
| `families/skew_group.py` + `hochschild/skew_group.py` (Plan-74 R8 — the skew group algebra `A⋊G` (smash `A#kG`): `QuiverAutomorphism` (a per-instance-certified algebra automorphism from a quiver `(π,ρ,s)` triple) + `GroupAction` (BFS closure, Cayley table, conjugacy classes, centralizers), `skew_group_algebra(A, action)` (structure constants `(a·g)(b·h)=a·g(b)·gh`, dim `|G|·dim A`, char-agnostic), `Algebra.skew_group`, `is_free_action`; the Ştefan conjugacy-class decomposition `stefan_decomposition` (`HH^n(A⋊G) ≅ ⊕_{[g]} HH^n(A,{}_gA)^{Z(g)}` via P52 `Bimodule.twisted` + the NEW `Z(g)`-action transport `_autom_action_on_classes` + Reynolds invariants), the `SkewGroupAlgebra` construction family + the `skew_group_hh` compute kind) — `tests/families/test_skew_group_action.py` (9), `test_skew_group_algebra.py` (5), `test_skew_group_free_action.py` (3), `tests/hochschild/test_skew_group_summands.py` (3), `test_skew_group_decomposition.py` (6), `tests/hochschild/test_skew_group_transport_basis.py` (5), `tests/qpa/test_skew_group_qpa.py` (3) | (see class table) | fast + deep + qpa | **`oracle_literature`** — `HH^•(k[x]/(x²))=[2,1,1,1]`; the **dim law** `dim(A⋊G)=|G|·dim A`; the orbit-Nakayama `HH^•=[1,1,1,1]`/`HH_•=[2,1,1,1]` (self-injective Nakayama); the **non-abelian `S₃`-on-`k³`** `HH^0=2` (semisimple, Morita `k×k`, dim `18`). **`oracle_crossengine`** — **the R8 oracle** `stefan_decomposition.dims == HH^•(A⋊G) direct` (`Z/2` on the dual numbers, both sides `[1,1,1,1]`/`[2,1,1,1]`); `HH^•(A⋊G) direct == HH^•(orbit Nakayama)`; the **genuinely FREE `Z/2`-swap** on `k(1⇄2)/rad²` → `presented_form == k[x]/(x²)`, `HH^•==[2,1,1]` (QQ); `Bimodule.regular ≡ ordinary` (P52 `M=A`); the twisted summand `HH^•(A,{}_σA)`; the **HOM `n=1` per-summand split** `[{e:0},{σ:1}]` (a NONZERO twisted invariant at `n≥1`, findings 1a/2b — a dead/sign-wrong bridge fails it even though the total is unchanged); the `Z/3` on `k[x]/(x³)`/GF(7) non-involution direct (`[1,1,1]`/`[3,2,2]`, `HH_0` split `1+1+1`). **`oracle_selfcert`** — the dim-law certificate; `GroupAction.check` (each generator an algebra automorphism preserving `I`; closure a group of the claimed order over the field); trivial-`G` `T`-byte-identity (`A⋊{1}.T == A.T`); `presented_form` dim + multiplicativity; the Reynolds idempotent `e²=e`; the identity-summand monomorphism `dim HH^n(A)^G ≤ dim HH^n(A⋊G)` (Marcos–Martínez-Villa); the degree-0 hand-pin `HH^0(A⋊G)=1`; the **twist-symmetry** `HH^•(A,phi=g)==HH^•(A,psi=g⁻¹)` (order-2 AND order-3, so the twist convention is NOT summand-arbitrable); the modular `char \| |G|` loud refusal; the canonical-key generator-order normalization; **the transport's TWO-BASES arbiter** — on a `hbar` that is neither symmetric nor self-inverse (the `Z/3` 3-cycle rotation), EXACTLY ONE of the eight (coefficient-slot, bar-block) combinations is a chain map on each side, plus the fixture guard that fails if `hbar` ever goes symmetric/self-inverse (scope fact (k)); the vertex-permuting `Z/2`-swap decomposition ≡ DIRECT on both sides (the regression: the unit-adapted matrix on the coefficient slot made the transport leave the cycle span and the whole decomposition refuse). **QPA (`-m qpa`)** — a `NamesGVars()` fail-if-appears guard (QPA 1.37 has NO skew-group / smash / conjugacy-class-decomposition surface) + the input-level DIRECT `HH^•(A⋊G)` crosscheck: build `A⋊G`'s presented `kQ'/I'` (`presented_form`) and compare its QPA `HochschildCohomology` dims to our direct engine (`Z/2` on the dual numbers over QQ → `[1,1,1,1]`; `Z/3` on `k[x]/(x³)` over GF(13) → `[1,1,1]` — GF(13) so `char > dim 9` for `presented_form`). The DECOMPOSITION itself has no QPA counterpart (honest scope). |
| `hochschild/simplicial.py` + `families/poset.py` (Plan-75 R9 — the ORDER-COMPLEX route to `HH^*` of an incidence algebra: `OrderComplex.of(P)` (chains by dimension, face vector, integer boundary `d(x_0..x_k) = Σ(−1)^j (x_0..x̂_j..x_k)`), `simplicial_cohomology_dims` (Domain ranks), `integral_homology` (exact-integer Smith normal form) and `cohomology_dims_from_integral` (universal coefficients — ONE integer SNF answers EVERY characteristic), plus `Poset.order_complex` / `is_lattice` / `has_global_bound`; `Algebra.incidence_cohomology(top)` and the `incidence_cohomology` compute kind across all three tiers with the poset input mode) — `tests/families/test_order_complex.py` (7), `test_simplicial_cohomology.py` (10), `test_incidence_hh.py` (9), `tests/trace/test_incidence_report_p75.py` (7), `tests/webapp/test_incidence_p75.py`, `tests/gui/test_incidence_runner_twin.py`, `tests/qpa/test_incidence_koszul_qpa.py` | (see class table) | deep + fast + qpa | **`oracle_literature`** — the theorem `HH^n(kP) = H^n(Δ(P); k)` (Gerstenhaber–Schack 1983 for the face poset of a simplicial complex; Cibils 1989 for an ARBITRARY finite poset) reproduced on: the Boolean lattice `B₃` `[1,0,0,…]` (the Plan-33 pin PROMOTED from a coincidence to a theorem instance — `B₃` has a global bound, so `Δ(B₃)` is a cone), the crowns `C(3,3)`/`C(2,2)` `≅ S¹` `[1,1,0,…]` (**nonvanishing `H¹`** on NON-lattice posets), and the **char-sensitive `RP²₆`** split `HH^*(GF₂) = [1,1,1]` vs `HH^*(QQ) = [1,0,0]` with the integral certificate `H₁(RP²;ℤ) = ℤ/2`; the validity of the 6-vertex `RP²` triangulation itself (every edge in exactly 2 facets, `f = (6,15,10)`, `χ = 1`). **`oracle_crossengine`** — the fast path ≡ the general CS engine (`cs_cohomology_dims`) degreewise on every test poset incl. the `RP²₆` char split at dim 121; universal coefficients from the ONE integer SNF ≡ the direct field rank in char 0 / 2 / 3. **`oracle_selfcert`** — `∂∘∂ = 0` on the integer boundaries; the Euler characteristic `Σ(−1)^p f_p` == the alternating Betti sum; chains strictly increasing and enumerated once; the `has_global_bound` ⇒ cone certificate proved WEAKER than `is_bounded` (either bound alone cones the complex — the plan's Task-I3 sketch said both, and the weaker hypothesis is the correct, strictly more general one); `contractible` **three-valued** (`True` only when PROVED by a global bound, `None` for "vanished in the computed range but not proved", `False` otherwise — never guessed from an observation); `torsion`/`char_dependent` `None` when over the integral budget, never estimated; and the **provenance HARD GATE** — an algebra without `_poset` refuses loudly (a plain `kA₂` raises), because recognizing an arbitrary `kQ/I` as an incidence algebra is a different and harder problem and applying the theorem without checking its hypothesis would be a silent lie. **cross-runner contract (UNMARKED, Plan-32 extras-gated ruling)** — both runners emit the `incidence_cohomology` block byte-for-byte identical from the shared `incidence_cohomology_block`, and the poset panel's `elements` normalization keys a covers-only poset BYTE-IDENTICALLY to the typed `IncidenceAlgebra` family request (an isolated point keys differently — and is different mathematics, `H⁰ = k²`). **QPA** — a fail-if-appears probe: QPA 1.37 has no Hochschild surface and no order-complex/nerve surface at all (verified live via `NamesGVars()`). |
| `hochschild/koszul_ghms.py` (Plan-75 R10 — the GHMS comultiplicative minimal `A^e`-resolution of a KOSZUL algebra: `koszul_kernels(A, top)` / `koszul_betti` (the intersections `K_n = ⋂_i V^i ⊗ R ⊗ V^j` inside `V^{⊗n}`, built by the recursion `K_n = (V⊗K_{n−1}) ∩ (K_{n−1}⊗V)` — exact, ONE subspace meet per degree, coordinates on COMPOSABLE PATHS so the corner grading is automatic and `dim K_0 = |Q_0|`), `GHMSResolution` with the comultiplicative differential `d_n(1⊗ω⊗1) = Σ x_a⊗ω'_a⊗1 − (−1)^n Σ 1⊗ω''_b⊗y_b`, and `engine="ghms"` on both public HH methods — the third independent HH oracle class) — `tests/engine/test_koszul_kernels.py` (10), `test_ghms_resolution.py` (11), `test_ghms_hh.py` (13) | (see class table) | deep | **`oracle_literature`** — `dim K_n` is the Koszul-dual Hilbert coefficient: `Λ(k²)` `[1,2,3,4,5,6,7]`, `Λ(k³)` `[1,3,6,10,15,21,28] = C(n+2,2)`. **`oracle_crossengine`** — `engine="ghms"` `HH_*` ≡ the minimal syzygy engine (`Λ(k²)` `[3,4,6,8,10,12,14]`, `Λ(k³)` `[5,12,24,40,60]`, multi-vertex `kZ₃/rad²` `[3,0,1,1,0,0]`, the diamond incidence `[4,0,0,0,0]`) and `HH^*` ≡ CS over QQ on all four; `dim K_n` == the Plan-27 `ext_algebra` graded Betti numbers; the diamond incidence algebra is the explicit BRIDGE between the two halves of this plan (`K_n = [4,4,1,0,0,0]`, `HH_* = [4,0,0,0,0]` — NOT to be conflated with its order-complex `HH^* = [1,0,0]`). **`oracle_selfcert`** — `d∘d = 0` over several Domains; the differential is canonical (byte-reproducible); the **sign is ARBITRATED, not assumed** — `d∘d = 0` pins the ALTERNATION (a constant sign IS caught) but is BLIND to a global flip, and the AUGMENTATION `μ∘d₁ = 0` (`d₁(1⊗a⊗1) = a⊗1 − 1⊗a`) is what pins the global sign, so both probes are tests: shipping the plan's written `−(−1)^n` passed every `d∘d` check and still gave `HH₀(Λ(k²)) = 2` instead of 3; the THREE-VALUED Koszul gate — `engine="ghms"` refuses on `False` **naming the Ext obstruction** and refuses on `None` ("Koszulity unknown"), never on `g_quadratic` alone; the non-quadratic refusal; and `engine="auto"` deliberately NOT routed to GHMS (byte stability of every shipped result, pinned by a test). |
| `cluster/category.py` (Plan-79 R31 — the **Amiot–Keller cluster category**, certified module-category slice: the finite fundamental domain `ind(C_Q) = ind(mod kQ) ⊔ {P_v[1]}` (BMRRT) giving `#indec = #ind(mod kQ) + n`; cluster-tilting objects **as** support τ-tilting pairs (AIR) counted through the shipped P45 exchange graph with a RUNTIME certification gate; the cluster-tilted End-algebra `End_C(T) = Jac(Q_T, W_T)` (BMR/Amiot) via shipped FZ matrix mutation; the module-window 2-CY certificate via the AR formula; `from_potential`'s algebra-level Amiot/Keller–Reiten certificate; `Algebra.cluster_category()` and the `cluster_category` kind across all three tiers) — `tests/cluster/test_category_indec.py` (22), `test_cluster_tilting_air.py` (14), `test_cluster_tilted_endalg.py` (10), `test_two_cy.py` (12), `test_from_potential.py` (8), `test_public_surface.py` (3), `tests/webapp/test_cluster_category_kind_p79.py`, `tests/gui/test_cluster_runner_twin_p79.py`, `tests/qpa/test_cluster_category_qpa.py` | (see class table) | deep + fast + qpa | **`oracle_literature`** — `#indec(C_Q)` = the almost-positive-root count: `A₂`=5, `A₃`=9, `A₄`=14, `A₅`=20 (`= n(n+3)/2`), `D₄`=16, `D₅`=25 (`= n²`); the **cluster numbers** `A₂`=5, `A₃`=14, `A₄`=42, `A₅`=132 (`= Catalan(n+1)`) and the **Buan–Marsh** `D₄`=50, `D₅`=182; `Jac(3-cycle, αβγ) = kZ₃/J²` with dim 6 and self-injective; **Keller–Reiten** Gorenstein-dim ≤ 1 on every cluster-tilted algebra built — exercised at BOTH values (0 on the self-injective flagship, **1** on the multi-3-cycle `A₅[1,3]` instance, so the bound is not just the self-injective case restated); the flagship's infinite global dimension (Keller–Reiten's hereditary-iff-finite contrapositive); `A₃`'s five nonzero `Ext¹` pairs. **`oracle_crossengine`** — the **AIR bijection is the backbone**: the P45 exchange graph is built by τ-tilting mutation and knows nothing about cluster combinatorics, yet its vertex count IS the cluster number, so `count == len(exchange_graph.vertices) ==` Catalan/Buan–Marsh is a genuine two-route agreement; the flagship's Cartan matrix equals the independent `NakayamaAlgebra(3,2,cyclic)`'s **up to a vertex permutation** (see honest scope — raw equality is FALSE and the test asserts the raw matrices differ, so the weakening cannot pass unnoticed); the 2-CY `Ext¹` table equals an INDEPENDENT `modules.ext` recompute outside the certificate's own loop; the `from_potential` flagship's category invariants equal the direct `A₃` route (9 and 14, computed by different entry points). **`oracle_selfcert`** — the almost-positive-roots identity `#indec = #ind(mod kQ) + n`; the **Auslander–Reiten formula** `dim Ext¹_A(X,Y) = dim Hom_A(Y, τX)` with **0 mismatches** across `A₂`–`A₅`, `D₄` and a **zigzag-oriented `A₄`** (the evidence for the ordinary-Hom proxy, REPORTED via a `mismatches` list rather than assumed); the symmetric `Ext¹_C = e_XY + e_YX` and the payload's statement that this half is BY CONSTRUCTION; the **FZ quiver round-trip** (rebuild from the mutated matrix, re-read its exchange matrix, require the original back); the **two readings of BBK-style `q`** — the `D₄`/`D₅` n-regularity RECOVERY is load-bearing, not decoration; the **budget-vs-infinite split** (G3a quotes the cluster number and never says infinite; G3b claims infiniteness ONLY off the Gabriel/Dynkin certificate); the **D₅ certified-or-refused gate**; `count is None` on EVERY refusal with the literature value quarantined in `expected_count`; `cluster_tilting_objects()` refusing rather than handing out a partial enumeration; the G5 never-guess-a-potential refusal; the G7 `D^b(Γ)` refusal that still exposes the algebra certificate. **cross-runner contract (UNMARKED, Plan-32 extras-gated ruling)** — both runners emit the `cluster_category` block byte-for-byte identical from the shared `cluster_category_block` (proved by driving the Pyodide twin's OWN `run_build` + `compute_one`); the budget grammar parses identically in ALL THREE grammar sites; the estimator sizes on the CLUSTER NUMBER, not `dim kQ` (`A₃`/`A₄` instant, `A₅`/`A₇`/`D₅` queued, non-cluster requests unmoved); a refusal is a typed error block or a populated field-level refusal, never a 500. **QPA and Macaulay2 have NO cluster-category surface** — no cluster category, no orbit category, no cluster-tilting enumeration, no cluster-tilted-algebra builder, no quiver-with-potential Jacobian constructor (probed live via `NamesGVars()`, failing loudly if that changes). QPA offers classical tilting one level below (`TiltingModule`, `DTr`). |
| `engine/complete_resolution.py` + `hochschild/tate.py` (Plan-76 R3 — **Tate–Hochschild (singular Hochschild) cohomology in EVERY integer degree**: the complete resolution of `A` over `A^e` built by splicing the minimal `A^e`-resolution `P_•` with its Nakayama-twisted `A^e`-dual `T_{-n-1} = Ψ(θ^*(P_n))` (`Ψ` = the k-dual made a LEFT `A^e`-module through the swap anti-automorphism, `θ` = the twist with `θ^*(A) = D(A)`), tags `(v,w) ↦ (π(w), v)`, entries `f_a⊗f_b ↦ ν(f_b)⊗f_a` transposed, and the splice joint `P_0 ↠ A ↪ T_{-1}` solved exactly as `Z = G^{-1}RG^{-1}` off the Frobenius form; `tate_cohomology_dims` / `tate_homology_dims`; the public three-route `Algebra.tate_hochschild(top, engine=…)` (`native` / `duality` / `positive`), `tate_periodicity`, and the `tate_hochschild` compute kind on all three tiers) — `tests/engine/test_complete_resolution_p76.py` (26), `tests/hochschild/test_tate_hochschild_p76.py` (20), `tests/webapp/test_tate_kind_p76.py` (8) | (see class table) | deep + fast | **`oracle_literature`** — the `k[x]/(x^n)` FULL RING `dim ĤH^m = n−1` (char ∤ n) / `n` (char ∣ n) for **every** `m ∈ [−5,5]`, over primes `32003 / 2 / 3` (`n = 2,3,4`), derived from the explicit period-2 `A^e`-resolution; the **`ĤH^0 ≠ HH^0` witness** (`k[x]/(x³)`: `ĤH^0 = 2` against `HH^0 = dim Z(A) = 3` — the Tate degree 0 is the STABLE centre, a proper quotient); the Bergh–Jorgensen threshold `agrees_from = d+1 = 1` for self-injective input; the QuantumCI `QQ` positive-threshold pin `ĤH^{n≥1} = HH^n = [2,1,0,0,0]` (Bergh–Jorgensen's full-ring value `1,2,1` in degrees `0,1,2` and `0` elsewhere is a char-0 literature statement — see honest scope). **`oracle_crossengine`** — **positive-degree agreement** `ĤH^m == HH^m` for `m ≥ 1` between the NEW splice engine and the shipped bar/fast HH engine, on single-vertex, symmetric multi-vertex (`kZ₂/J³`) and NON-symmetric multi-vertex (`kZ₃/J²`) input; the `duality` route ≡ the `native` route in every degree the duality claims. **`oracle_selfcert`** — `d∘d = 0` across every joint including both splice joints; **in-window acyclicity** (rank–nullity per degree — the DEFINING property of a complete resolution); build-time **corner typing** (each entry lies in `ε_source·A^e·ε_target`, a genuine cross-check because the tags come from the socle-derived `π` and the entries from the form-derived `ν`); the symmetric duality self-consistency `dim ĤH^{-j} == dim ĤH^{j-1}`; the Frobenius HOMOLOGY duality `dim ĤH_n == dim ĤH_{-(n+1)}` (valid for every Frobenius algebra, checked on the non-symmetric witness too); the non-symmetric negatives are NOT the symmetric mirror; the periodicity certificate read off the DIFFERENTIALS (equal Betti numbers are necessary and not sufficient) with an honest `None` outside budget; and the loud refusals (DEFERRED Gorenstein-non-self-injective, non-Gorenstein, non-symmetric on `duality`, non-`GF(p)` on `native`, presentation-less, unknown engine, negative top). **THE ARBITER IS ACYCLICITY, NOT `d∘d = 0`** — measured and pinned as a live negative result (`test_dd_zero_alone_does_not_discriminate_pi`): the negative half is carried as `A^e`-module maps in AMBIENT coordinates and the composite is taken there WITHOUT reference to the tags, so a wrong Nakayama permutation passes `assert_dd_zero` and still reports wrong dimensions; on `kZ₃/J²` the `π = identity` control is non-exact in EVERY negative degree. **cross-runner contract (UNMARKED, Plan-32 extras-gated ruling)** — both runners emit the `tate_hochschild` block byte-for-byte identical from the shared `tate_hochschild_block`, a route's unknown degrees are `null` (rendered as an em dash), and a refusal is a typed `error` block, never a 500. **QPA has NO Tate / singular-Hochschild / complete-resolution surface** (honest scope). |
| `modules/nkoszul.py` (Plan-77 R36 — the generalized Koszulity ladder BEYOND the quadratic case, a thin recognizer layer over the shipped Plan-05 minimal resolutions and the Plan-27 Yoneda engine: `generation_degrees` (the INTERNAL path-length degrees `ℓ_i(n)` of the minimal resolution of each simple — the primitive Plan 27 does not have), `n_homogeneous_degree`, `n_koszul_certificate` (Berger's 2-N alternation), `k2_certificate` (Cassidy–Shelton, through an explicit certified window), `almost_koszul_certificate` (Brenner–Butler–King `(p,q)`), `multi_koszul_certificate` (Herscovich, scoped), `koszul_profile[_block]` and `Algebra.koszul_profile(top)`, plus the `koszul` compute kind across all three tiers) — `tests/modules/test_nkoszul_generation_degrees.py` (10), `test_nkoszul_ncertificate.py` (12), `test_nkoszul_k2.py` (10), `test_nkoszul_almost.py` (24), `test_nkoszul_multi.py` (6), `test_nkoszul_profile.py` (8), `tests/webapp/test_koszul_p77.py`, `tests/gui/test_koszul_runner_twin.py`, `tests/qpa/test_nkoszul_qpa.py` | (see class table) | deep + fast + qpa | **`oracle_literature`** — Berger's closed form `δ(n) = (N/2)n` (`n` even) / `(N/2)(n−1)+1` (`n` odd) reproduced EXACTLY as the internal generation degrees of `k[x]/(x^N)` for `N = 2,3,4,5` (`[0,1,2,3,…]`, `[0,1,3,4,6,7,9,10,12]`, `[0,1,4,5,8,9,…]`, `[0,1,5,6,10,11,…]`), and on the multi-vertex finite-gl.dim `kA₄/J³` (`S₁: [0,1,3]`, `N = 3`); the Cassidy–Shelton K₂ criterion (`E(A)` generated in cohomological degrees 1,2) on `k[x]/x^N` (generators in `[1,2]`) and `kA₄/J³` (complete `True`); the **BBK preprojective classification** `(h−2, 2)` reproduced via `p = ` top degree and `q = e − p` for `Π(A₃) = (2,2)`, `Π(A₄) = (3,2)`, `Π(A₅) = (4,2)`, `Π(D₄) = (4,2)`, with `p = h−2` verified independently for `Π(A₂..A₆)` and `Π(D₄)`; `Π(A_{n≥3})` is K₂-**False** (a genuine Yoneda generator in cohomological degree 3) and `Π(A₂)` is the Koszul boundary. **`oracle_crossengine`** — Berger's `N ≥ 3` characterization makes the TWO N-Koszul certificates independent readings of one property, and they are required to agree live (`k2_agrees`): the purity side reads INTERNAL degrees off the Plan-05 minimal resolution, the K₂ side reads HOMOLOGICAL generator degrees off the Plan-27 Yoneda engine; a DISAGREEMENT is reported as an inconsistency rather than resolved to one side. `k2_certificate`'s `generator_degrees` / `window` are Plan-27's `generators_by_degree` / `certified_through_degree` VERBATIM, and the whole quadratic half of the profile (`koszul`, `koszul_reason`, `koszul_obstruction`) is Plan 27's verdict verbatim — the NAMED overlap, pinned equal, with Plan 27 left byte-unchanged. **`oracle_selfcert`** — `ℓ_i(n) ≥ n` always (minimality ⇒ every differential has path length ≥ 1), with equality exactly on a linear strand; the extraction's per-column homogeneity (`_NonPure` stays dormant on graded-pure input — kept as a DEFENSIVE assert, since any re-derivation shares its single method, so Berger's closed form is the real anchor); Berger's alternation `1, N−1, 1, N−1, …` from the closed form; the **two readings of BBK's `q` must AGREE** (`q = e − p` and `q = n★ − 1` — BBK puts the error of internal degree `p+q` at homological step `q+1`), and `verdict=True` is GATED on that agreement; the **P75 seam** (`ext_algebra.koszul`'s obstruction degree == the almost-Koszul break degree, both 3 on `Π(A₃)`/`Π(A₄)` — a structural-consistency invariant, ONE first-failure-of-linearity degree read two ways, not two independent computations); K₂ vacuity on hereditary input; the length-graded refusal (`x³ − x²` is admissible but inhomogeneous and raises); the multi-Koszul verdict is `None` on EVERY shipped shape (no branch may invent one); and the `latex` summary must SAY "through degree W" whenever the claim is window-bounded. **cross-runner contract (UNMARKED, Plan-32 extras-gated ruling)** — both runners emit the `koszul` block byte-for-byte identical from the shared `koszul_profile_block` (proved by driving the Pyodide twin's OWN `run_build` + `compute_one`, not by re-calling the library), `koszul` adds no request field so pre-P77 requests key byte-unchanged, and a presentation-less algebra becomes a typed `error` block, never a 500. **QPA has NO generalized Koszulity surface** — no `IsKoszul` (the Plan-27 finding, re-probed so it cannot rot), no `IsNKoszul`/`IsDKoszul`, no K₂, no almost-Koszul, no multi-Koszul, and nothing reporting internal generation degrees (verified live via `NamesGVars()`, bound names only). It DOES ship `ExtAlgebraGenerators`, and since `k2_certificate` reports Plan-27's `generators_by_degree` verbatim, the existing Plan-27 `ExtAlgebraGenerators` crosscheck IS a live QPA leg under K₂'s input data. |
| `batch/` (labdb port, open-zone scans) | 11 | deep | labdb port equality; scan-surface checks |
| `citations/` (registry, bibliography) | 12 | fast | packaged-bib resolution; result references |
| `trace/` (worked-steps incl. the Plan-30 module events, the kA₂ replay golden, the 2026-07-29 report-completeness battery, the Plan-35 UNIT-2 HH explicit-reps rendering, the Plan-35 wave-3a Ext/Tor explicit-reps rendering, the Plan-35 wave-3b cyclic-homology explicit-reps rendering — the total-complex `Tot_n = C_n ⊕ C_{n-2} ⊕ …` column heading, per-degree classes + verification; and the Plan-35 wave-3c Yoneda-sequence + classical-dictionary rendering — `interpretations.py`; and the Plan-35 wave-3d plain-HH explicit-reps + element-wise dictionary rendering — `hh_element_interpretation`/`hh_reps_sections`) | 250 | fast | golden-file equality (dims derived from ranks); **the per-degree explicit-reps layout** (each product/Connes class rendered as term-sum + coordinate vector under a stable anchor, with the annihilating differential + a one-line verification sentence; the bar AND Chouhy-Solotar HH worked-steps carry each (co)chain term's ordered basis, length-guarded against the recorded term dim; module resolution `term_basis` lengths match the differential row/col dims, injective order pinned against the transposed proj-resolution-of-DM; the degree anchors are linked from every product table) + **the module Ext/Tor per-degree sections** (ordered Hom/tensor basis → classes → differential + verification, `cr-`/`ws-` anchors, the `ExtReps` worked-steps event, Tor₀ = M ⊗ N cokernel note) + **the Yoneda-sequence + dictionary rendering** (each Ext class' constructed exact sequence — sequence line, middle module, exactness verified — under `cr-ext-yoneda-deg-n`; the shared classical-dictionary framing on the ext/tor/HH/cyclic blocks; the HH¹ derivation read-off; matrix-grid double zebra striping is structure-safe) + **the plain-HH element-wise dictionary + per-degree reps** (HH⁰'s central elements, HH¹'s `D(arrow)=value` derivations + the inner-derivation subspace dimension `rank δ⁰`, HH²'s deformation 2-cocycle, HH₀'s commutator residues — read straight off the captured term-sums; the per-degree explicit-reps sections under `cr-hh_cohomology`/`cr-hh_homology` anchors; both gui.js copies mirror it) + the missing-fields tolerance + the two-runner `term_basis`/reps/interpretation equality + **the renderer-return gate** (2026-08-06: every per-kind `*_html` helper in `results_html.py` must end in an explicit `return` — an AST scan plus a live `derived_fingerprint` render through `results_section`; regression for the fall-through that returned `None` and killed the whole worked-steps bundle for any request containing `derived_fingerprint`) |
| `specseq/` (Plan-42 spectral sequences — `filtered.py`/`double.py` filtered & double complexes, `pages.py` the Weibel-5.4 page engine, `convergence.py` the standing self-certificate, `presets.py` the four presets, `block.py` the `ss_hochschild` no-code block) | 33 | deep | **self-cert** (`d_r∘d_r=0`, `E_{r+1}=H(E_r,d_r)`, `E_∞` totals == total homology on every construction, canonical-rep reproducibility, the radical-SS converges to `H(X)`, the subcomplex-filtration + double-complex anticommutation gates); **cross-engine** (the Hochschild `(b,B)` `E_∞` total == `A.cyclic_homology`, the Cartan–Eilenberg/Grothendieck `E_∞` total == module `A.ext` on several instances incl. a multi-vertex one + NONZERO pins); **literature** (ground-field `HC=[1,0,1,0,…]`, `k[x]/(x²)` HC, the arbitrated Koszul `E_2` degeneration); **`m2`** (the commutative Koszul total-complex `E_∞` totals vs Macaulay2 `Complexes` homology) |
| `viz/` (draw, tikz; the Plan-49 generic Hasse twin — `layout.poset_layout` + `tikz_hasse` + `hasse_html.hasse_svg`, tested in `tests/modules/test_degeneration_render.py`) | 18 | fast | exact `int`/`Fraction` layout; TikZ; **the Hasse twin** (`oracle_selfcert`): the layout ranks the poset (minimum at rank 0, every cover spans exactly one rank), and `tikz_hasse`/`hasse_svg` emit non-empty markup naming every class + drawing every cover, float-free |
| `qpa/` (GAP/QPA crosscheck) | 158 | 137 qpa + 21 fast | **live GAP/QPA** (HH dims, self-Ext, τ/τ⁻, proj/inj resolutions, inj dim, Plan-31 native trivial-extension construction — left side via `A^op`; Plan-38 `IsSpecialBiserialAlgebra`/`IsGentleAlgebra`; the Plan-37 **hom-glue battery** — `Length(HomOverAlgebra)` vs our `hom_basis` dim, and for a canonical dim-1 hom the kernel/image/cokernel dimension vectors vs QPA `KernelInclusion`/`ImageInclusion`/`CoKernelProjection` over kA₂/kA₃(ab)/`line_abc_cde`; the Plan-39 **complexes battery** — QPA 1.37's Ch.10 `StalkComplex`/`FiniteComplex`/`HomologyOfComplex`/`Shift` vs our `ChainComplex`: stalk homology, our mapping cone `[M→N]` via the equivalent `FiniteComplex` — QPA's `MappingCone` object is not homology-scriptable through libgap, the documented fallback — and `Shift` bookkeeping under QPA's OPPOSITE `−k` convention); script builders + guards run without GAP |
| `qpa/` (GAP/QPA crosscheck) | 154 | 133 qpa + 21 fast | **live GAP/QPA** (HH dims, self-Ext, τ/τ⁻, proj/inj resolutions, inj dim, Plan-31 native trivial-extension construction — left side via `A^op`; the Plan-37 **hom-glue battery** — `Length(HomOverAlgebra)` vs our `hom_basis` dim, and for a canonical dim-1 hom the kernel/image/cokernel dimension vectors vs QPA `KernelInclusion`/`ImageInclusion`/`CoKernelProjection` over kA₂/kA₃(ab)/`line_abc_cde`; the Plan-41 **AR battery** — `AlmostSplitSequence` middle-term dimension vector (over QQ) + summand multiset (over GF(p), `DecomposeModule`) and `PredecessorsOfModule` immediate predecessors, on kA₃ and the linear Nakayama algebra, **plus the devil's-advocate non-brick case** `k[x]/(x⁴)` with `M = k[x]/(x²)` (`dim Ext¹(M,τM) = 2`), whose middle `{1,3}` is crosschecked against QPA `AlmostSplitSequence`/`DecomposeModule` over GF(32003)); script builders + guards run without GAP |
| `webapp/` (server tier + result cache + offline GUI — non-algebraic glue) | 426 | fast | API / schema / cache canonicalizer (replay-safety rests on exactness) / isolation / artifacts; all math delegated to the library; Plan-28 runner delegation pinned **byte-identical** (frozen goldens + unchanged `canonical_key`) |
| `qpa/` (GAP/QPA crosscheck) | 165 | 144 qpa + 21 fast | **live GAP/QPA** (HH dims, self-Ext, τ/τ⁻, proj/inj resolutions, inj dim, Plan-31 native trivial-extension construction — left side via `A^op`; Plan-38 `IsSpecialBiserialAlgebra`/`IsGentleAlgebra`; the Plan-37 **hom-glue battery** — `Length(HomOverAlgebra)` vs our `hom_basis` dim, and for a canonical dim-1 hom the kernel/image/cokernel dimension vectors vs QPA `KernelInclusion`/`ImageInclusion`/`CoKernelProjection` over kA₂/kA₃(ab)/`line_abc_cde`; the Plan-39 **complexes battery** — QPA 1.37's Ch.10 `StalkComplex`/`FiniteComplex`/`HomologyOfComplex`/`Shift` vs our `ChainComplex`: stalk homology, our mapping cone `[M→N]` via the equivalent `FiniteComplex` — QPA's `MappingCone` object is not homology-scriptable through libgap, the documented fallback — and `Shift` bookkeeping under QPA's OPPOSITE `−k` convention); script builders + guards run without GAP; the Plan-44 tilting/approximation battery (`test_tilting_qpa.py`) — `is_tilting_module` vs the computational `TiltingModule(T, n) <> false`, and `right_add_approximation`/`left_add_approximation` vs `MinimalRightAddMApproximation`/`MinimalLeftAddMApproximation` (source/range dimension vectors) |
| `hpc/` (headless CLI + spec core + container assets — non-algebraic glue, Plan 28) | 64 | fast (checkpoint-resume: deep) | **CLI ≡ public-API parity** on fixture configs; renderer golden tokens (LaTeX/HTML/text ladder); checkpoint-resume end-to-end equals the uninterrupted run; import-boundary + exit-code contract; sbatch/Dockerfile/workflow asset gates |
| `docs/gui/` (Pyodide GUI + no-code module panel — non-algebraic glue) | 81 | fast | runner artifacts / invariants; build hook; freshness; the two-runner Ext/Tor reps equality + the two-runner cyclic-homology reps equality + the Plan-44 two-runner `tilting_check` math-subkey parity |
| release + top-level (`test_no_floats`, `test_errors`, `test_quickstart`; the Plan-32 `test_oracle_classes` audit gate) | 57 | fast (audit gates: deep) | **float-ban AST gate**; error taxonomy; packaging; docs-nav coverage; **oracle-class count audit** (page == live collection) |

Non-algebraic glue (`webapp/`, `docs/gui/`) carries no oracle *because it holds no
mathematics of its own* — it calls `import quiverlab` and is tested for correct
plumbing, not for algebra.

## Buckets and the CI matrix

Test buckets are auto-assigned by directory in `tests/conftest.py` (an explicit
marker wins); the partition is disjoint and exhaustive, enforced by a partition
test. Markers (`pyproject.toml`): `fast`, `deep`, `slow` (implies `deep`), `qpa`,
`m2`; plus the **orthogonal** oracle-class markers below (which never change a bucket).

| Bucket | Tests | Runs where |
|---|---:|---|
| `fast` | 1630 | every CI cell: `{ubuntu, macos, windows} × py{3.10, 3.11, 3.12, 3.13}` |
| `deep` | 1676 | one Linux · py3.12 cell, **twice**: numba and pure (`QUIVERLAB_NO_NUMBA=1`) |
| `qpa` | 189 | weekly Linux · py3.12 job with GAP + QPA (`QUIVERLAB_REQUIRE_QPA=1`) |
| `m2` | 11 | Linux · py3.12 job with Macaulay2 (`QUIVERLAB_REQUIRE_M2=1`) |
| `slow` | 0 | opt-in (`-m slow`); rides the deep leg |

The `lint` CI job runs the float-gate and release-metadata tests standalone. The
docs site is built `--strict` in its own workflow, so any internals chapter or
page missing from the nav fails the build.

## Oracle classes as runnable markers

The two-oracle narrative above is also carried by **orthogonal pytest markers**
(Plan 32), so a reviewer can run each oracle class as a one-liner. These markers
classify *how* a test verifies, not *where it runs* — they are independent of the
`fast`/`deep`/`qpa` runtime buckets (adding them changed no bucket: the sweep was
byte-identical), and a test may carry **more than one** (a battery pins a literature
value *and* asserts cross-engine agreement in the same test, so it carries both).
The assignment lives at module level in each battery file; see
`docs/plans/2026-07-26-plan-32-oracle-markers.md` for the class boundary and every
edge-case ruling.

- **`oracle_literature`** — the pass criterion is a value or identity from the
  literature or classical theory, frozen as a constant the engine must reproduce:
  paper-pinned dims and closed forms, theorem identities (Coxeter/spectral tables,
  Happel trace, Theorem B/C, Cartan identities, symmetric ⇒ `HH^n = HH_n`,
  `dim Tor_n = dim Ext^n`, the `kZ_n/J^L` symmetry classification, Künneth, the
  `k[x]/(x²)` cup ring and the QuantumCI-BGMS cup dims), and the
  read-only bank's closed-form differentials. The marker face of **Class 1**.
- **`oracle_crossengine`** — two *independent* implementations are run and required
  to agree live: CS ≡ bar ≡ Bardzell ≡ minimal degreewise, numba ≡ pure and
  sparse ≡ dense parity, presented ≡ ⋉ iso-invariance, native ≡ transported
  cup/cap, generic-Domain ≡ GF(p) engine, the bar ≡ CS product tables (Plan 35), and
  the Connes λ-complex second model. The library-internal face of **Class 2**.
- **`oracle_selfcert`** — an internal mathematical certificate *is* the assertion:
  `d∘d = 0`, the CS order condition, canonicalization / adversarial-solver
  byte-reproducibility, dimension and iso certificates, the self-certifying
  Nakayama `λ`/`ν` identities, and the unit/Leibniz/module identities that arbitrate
  a sign convention (including the Plan-35 Gerstenhaber-algebra product batteries —
  graded commutativity, associativity, antisymmetry, cup-Leibniz, the cap module law,
  and `B²=0`; plus the 2026-08-03 report-presentation contract: a zero differential is STATED (`d = 0`), never drawn or echoed, Ext/Tor name the resolved module + resolution before any number, engine provenance lines gloss themselves, and the worked-resolution-steps chapter names the A^e-resolution it walks -- `tests/trace/test_report_fixes_m0803.py`; pass 2: tensor separators are ⊗-only with the bar-tensor vs CS-generator semantics stated in the typing paragraphs, product sections declare their basis and warn when it differs from the HH sections' route, Ext/Tor show the resolution of M before the data, every max_cells mention glosses itself, and the A^e-resolution chapters precede the computed results -- `tests/trace/test_report_notation_m0803b.py`). (These are the "self-certifying internal identities" of Class 1,
  surfaced as their own runnable class.)
- **`qpa`** — the existing bucket marker *is* the fourth oracle class: our value ≡
  live GAP/QPA. It needs no new marker; the live-QPA face of **Class 2**.
- **`m2`** — the Plan-36 bucket marker *is* the fifth oracle class: our value ≡ live
  **Macaulay2** (single-vertex nc graded dims via `AssociativeAlgebras`, commutative
  Ext via `freeResolution`), driven as a subprocess. Like `qpa` it is an external
  system, never double-marked with an `oracle_*` mark; the second live-external face
  of **Class 2**.

Everything else is **contract & infrastructure** (unmarked): refusal/error
surfaces, API and protocol contracts, the float-ban AST gate, freshness/interface
gates, the Gröbner admissibility certificate, worked-steps golden plumbing, the
foundational field/algebra/linear-algebra datatype contracts, and the
GUI/webapp/HPC/release/docs tiers.

The counts below are **audited against live collection** by
`tests/release/test_oracle_classes.py` (the badge==page doctrine, cf. the buckets):
if a future plan adds a battery and forgets to bump a number here, that test fails.
They overlap by design, so the union is smaller than their sum.

| Oracle class | Run | Tests | What agreement means |
|---|---|---:|---|
| Literature / theory pins | `-m oracle_literature` | 1271 | the engine reproduces a value/identity that exists outside the library |
| Cross-engine agreement | `-m oracle_crossengine` | 751 | two independent implementations compute the same thing and match live |
| Self-certifying certificates | `-m oracle_selfcert` | 1655 | an internal axiom (d∘d=0, canonicality, an arbitration identity) holds by construction |
| Live QPA / GAP | `-m qpa` | 249 | an independent external system (QPA) recomputes and agrees |
| Live Macaulay2 | `-m m2` | 11 | an independent external system (Macaulay2) recomputes and agrees |
| Any oracle class (union) | `-m "oracle_literature or oracle_crossengine or oracle_selfcert or qpa or m2"` | 3315 | the test is pinned by at least one oracle (the remaining tests are contract/infrastructure) |
Recounted 2026-08-17 on the Plan-79 branch (R31 — the Amiot–Keller cluster category,
certified acyclic slice), cut from dev's P77 tip. The preceding recount was the Plan-77 branch (R36 — the generalized Koszulity ladder:
Berger N-Koszul, Cassidy–Shelton K₂ through an explicit certified window,
Brenner–Butler–King `(p,q)`-almost-Koszul, Herscovich's scoped multi-Koszul, and the
internal generation degrees of `Ext•(k,k)`), cut from dev's P76 tip; the Plan-77
additions (lit +25 / xeng +8 / selfcert +43 / qpa +1 / m2 +0, union +70 — suite
5201 → **5292**, fast 2539 / deep 2493) fold into these live-collected totals. The
**P80 release gate** then added the three-tier compute-grammar parity gate (46 unmarked
contract tests, so no oracle class moved): suite **5338**, fast **2585**. The
previous recount was the Plan-77 branch (R36 generalized Koszulity), which added
lit +37 / xeng +8 / selfcert +44 / qpa +1, union +90 (suite 5087 → 5201). The
previous recount was the Plan-75 branch (R9 incidence-vs-nerve `HH^*` + R10 the GHMS
comultiplicative Koszul resolution and `engine="ghms"`), which added
lit +17 / xeng +16 / selfcert +44 / qpa +1 / m2 +0, union +78 (suite 4925 → 5087). The
previous recount was the Plan-73 top-up merge (R7 Han's conjecture transported across an
arrow-removal subalgebra B ⊆ A — the bounded-extension certificate, the Jacobi–Zariski
relative complex, and the five-row injection/iso ladder), which added
lit +17 / xeng +4 / selfcert +16 / qpa +0 / m2 +0, union +37.

Counts as of the P71 merge (HH• as a graded Lie module over HH¹ R12, after P69 persistence/TDA R33 / P72 split-extension+arrow-removal R5+R6 / P70 HH¹-Lie / P64 congruences / P65 exceptional sequences / P67 silting / P54 BV / P68 skew-gentle / P63 wall-and-chamber);
sibling plans in the v0.2.0

Collected 2026-08-05 (through the Wave-1 v0.2.0 merges: Plans 36, 37, 38). The oracle markers live only on the
Collected 2026-08-05 (Plan 40, C6 homological-dimensions family; recounted from a live
collection on the `plan-40-homdims` branch -- the pre-Plan-40 numbers were mid-merge-train
Plan-39 complex layer + the Plan-42 spectral-sequence engine). **Mid-merge-train
pure-library `engine` / `resolutions_cs` / `hochschild` / `modules` / `invariants` /
`families` / `batch` / `trace` / `specseq` suites (the `trace` renderer tests import the
pure-library serializers only, and their `hpc.spec` uses are function-local), so these
counts do **not** depend on the `[web]`/`[hpc]` extras.

## The standing rule

**Every future plan adds its new oracles to this page as part of its acceptance** —
exactly as every plan already updates the "Under the hood" internals chapters. When
a plan ships a new engine, invariant, or operation, its acceptance task extends the
tables above with the oracle that guards it and the test file that runs it — and
updates the audited counts here and in the README tests badge (a release test pins
the two numbers equal, so a stale badge fails the suite). This
page is the single living record of how each shipped feature is verified, and it is
kept honest: if a subsystem lacks an oracle, this page says so rather than implying
one.

**Every literature oracle carries its citation.** We cite the literature we test
against, at the precision the repository can actually verify — author, year, venue,
and a theorem / example / proposition number *only when it is actually recorded in
a test, docstring, plan doc, or the read-only bank's attribution* — never a guessed
number. Where a source already has a citations-registry key
(`src/quiverlab/citations/`), the entry is the packaged `references.bib`, cited by
that key so it stays consistent with `quiverlab.bibliography(...)` and the
[References page](bibliography.md); where it does not, the source is named at the
verified precision and listed below as such.

## Honest scope

- **Cluster categories (Plan 79) — six honest boundaries.** (1) **The dg machinery is
  out of scope**: the Ginzburg dg algebra `Γ(Q,W)`, `D^b(Γ)`, the derived orbit
  category's triangulated structure and the DIRECT `Hom_{C}` are never formed — no dg
  engine ships. Everything computed here rides the FINITE mod-`kQ` model (BMRRT), valid
  only in the acyclic/hereditary case; a general non-acyclic `C_{(Q,W)}` refuses loudly
  (G7) while still exposing the Jacobian algebra's Amiot/Keller–Reiten certificate.
  (2) **`End_C(T) = Jac(Q_T, W_T)` is CITED, not computed** (BMR/Amiot). What is verified
  per instance is the FZ quiver mutation (a round-trip certificate) and the Jacobian
  algebra's own presentation — never the isomorphism itself. (3) **The potential is never
  guessed**: the canonical sum-of-oriented-3-cycles is the TYPE-A potential, so outside
  type A (or where no 3-cycle exists) the FZ-certified quiver is returned with
  `algebra=None` and a note — inferring a general potential needs DWZ mutation /
  right-equivalence, deferred since P48.1. (4) **`budget` ≠ `infinite`**: a budget stop on
  a known-finite Dynkin input is refused as OVER-BUDGET with the cluster number quoted in
  `expected_count`, NEVER as an infiniteness claim; infiniteness is asserted only off the
  Gabriel/Dynkin type certificate. And on ANY refusal `count` is `None` — the discovered
  vertex count is never emitted as if certified (the `D₅` gate). (5) **The 2-CY verdict is
  a MODULE-WINDOW verdict**, not the categorical statement: the AR formula is checked on
  every ordered pair of `ind(mod kQ)`, but the pairs involving the shifted `P_v[1]` hold
  by BMRRT — **cited, not computed** — and the payload always says so. The AR check itself
  uses ORDINARY `Hom(Y, τX)` as a theorem-backed proxy for the injectively-stable
  `\overline{Hom}` of ASS Thm IV.2.13; on a representation-finite hereditary algebra the
  two coincide across the tested window, and the evidence is the reported 0-mismatch
  sweep, not an assumption. (6) **The `from_potential` reduction is narrow and certified**:
  only the documented 3-cycle instance (gated on the BUILT algebra's dimension, not the
  quiver's shape) is recognized as reducible to a hereditary model; a general
  mutation-equivalence-to-acyclic test is deliberately not attempted. Plus: **the P45
  `mutate` root cause is NOT owned here** — `D₄` still reports `status="error"` on the
  current tree, so the Plan-63 n-regularity recovery remains load-bearing.

- **Generalized Koszulity (Plan 77) — six honest boundaries.** (1) The internal-degree
  recognizers (`generation_degrees`, N-Koszul, almost-Koszul) require a **length-graded**
  `A`: if `I` is inhomogeneous the minimal resolution is not graded and an internal degree
  simply does not exist, so they refuse loudly rather than report a pattern. **K₂ needs
  only homological data**, so it is defined without that hypothesis. (2) **The certified
  window is contractual.** When `gl.dim A` is finite and exact, `E(A)` is
  finite-dimensional and the verdicts are decisive; otherwise every affirmative verdict is
  three-valued `None` meaning *"holds through degree `W`"* and is NEVER printed as an
  unconditional `True` — this is why `k[x]/(x^N)` (self-injective) reports N-Koszul as
  `None`, not `True`, even though its internal degrees match Berger's `δ(n)` exactly. A
  definitive `False` (a generator or pattern break **inside** `W`) stays `False`. (3)
  **`almost_koszul_certificate` recognizes the SIGNATURE, not BBK's full definition.** It
  checks a linear prefix followed by a single break to internal degree `e = p + q` with
  `q ≥ 2` — what BBK's definition implies pointwise — and does **not** verify the
  existence of the finite linear complex of projectives with one error term, nor the
  `2(h−1)` periodicity. The reported `jump_spacing` is **window-observed** and is
  deliberately NOT called a period: for `Π(A₃)` the observed spacing is 3 while BBK's
  cited periodicity is `2(h−1) = 6`. (4) **`q = 1` is refused as almost-Koszul** — that is
  BBK's degenerate Koszul-type boundary (the whole `k[x]/x^N` family), caught by the
  N-Koszul recognizer instead; and an algebra with **no break at all** (`Π(A₂)`,
  `k[x]/x²` — both radical-square-zero, hence genuinely Koszul) reports `None` with `q`
  `None`, never a fabricated label. (5) **Multi-Koszul is connected-graded only.**
  Herscovich's notion assumes `A₀ = k`, i.e. a local single-vertex `kQ/I`; multi-vertex
  `kQ/I` has `A₀ = k^{Q₀}` semisimple and refuses loudly, pointing at K₂ — the transfer
  Herscovich Prop. 3.30 pins (multi-Koszul ⟹ K₂). The multi-Koszul **decision** itself
  requires Tor/Ext-vanishing of the minimal graded BIMODULE resolution (Herscovich §3.2),
  a homological engine of its own: it is **scoped out** with `verdict` always `None`, a
  named backlog entry, and no guess — and Herscovich's canonical examples (Yang–Mills,
  super-Yang–Mills) are **infinite-dimensional**, outside the f.d. engine, so f.d. local
  multi-Koszul examples are genuinely sparse. (6) **QPA cannot compare**: QPA 1.37 ships
  no Koszul surface of any kind (no `IsKoszul`, `IsNKoszul`/`IsDKoszul`, K₂,
  almost-Koszul, multi-Koszul, and nothing reporting internal degrees), probed live and
  failing loudly if that ever changes. The one live QPA leg is one level down —
  `ExtAlgebraGenerators`, which is exactly K₂'s input data.

- `complexity` is a lower-bound estimate that can under-report; it is **exact only
  on local / single-vertex inputs**. The Betti-complex identity it rests on
  (`H_n =` the minimal resolution's ranks) is gated exactly over every field.
- `is_symmetric` off GF(p) decides the definitional "`ν` is inner" by a
  Schwartz–Zippel sweep — **loud when inconclusive**, never a silent wrong answer.
- **Tate–Hochschild (Plan 76) — three honest boundaries.** (1) The `native` complete-resolution
  ring serves **self-injective** algebras over `GF(p)` only: `D_{A^e}(P_n)` is projective
  exactly when `A^e` is self-injective, so the dual-splice is not even well defined off that
  scope. A Gorenstein but NOT self-injective algebra is a **loud DEFERRED refusal** naming the
  periodicity-extension construction (Usui) as the follow-up — never a silent partial answer;
  its eventual-periodicity certificate is still computed, since that is a property of the
  algebra rather than of a route. (2) **`ĤH^0` and `ĤH^{-1}` are native-only.** They form a
  closed 2-cycle under the Bergh–Jorgensen reflection `n ↔ −(n+1)`, which links them to each
  other and to no threshold-known positive degree, so the `duality` and `positive` routes
  return them absent. Reporting `HH^0` there would be **wrong, not merely imprecise**: for
  `k[x]/(x^n)` the true `ĤH^0` is `n−1` while `HH^0` is `n`. Off self-injective input the
  threshold degree itself would need `gorenstein_dimension(A^e)`; rather than guess it, the
  surface reports `agrees_from = None` and claims no Tate value at all. (3) **No Tate cup
  product ships in this slice** — the `ℤ`-graded ring structure needs a diagonal approximation
  on the two-sided complex (the Plan-20/21 diagonal extended past degree 0), which is deferred;
  `TateHochschild.cup` is always `None`, and `periodicity_degree` is the degree of the
  resolution's certified period, cited as Usui's *criterion* for an invertible homogeneous
  element, NOT as an inverse this library computed. The Bergh–Jorgensen QuantumCI full-ring
  value is likewise a char-0 **literature statement**: `q` is a root of unity over every
  `GF(p)`, so the native engine cannot reproduce it there, and over `QQ` only the
  positive-threshold part is pinned.
- **QPA 1.37 has no Tate / singular-Hochschild / complete-resolution surface at all**, so every
  oracle above is a literature pin, a cross-engine agreement against the shipped HH engines, or
  an internal certificate.
- Live QPA cross-checks run **weekly**, not per commit (GAP is heavy). They are not
  silently skipped: the dedicated job makes an absent QPA a hard failure, and the
  frozen-value validation (`test_qpa_reference_validation.py`) runs in every matrix
  cell as the always-on stand-in.
- **Macaulay2 cannot see multi-vertex algebras or Hochschild anything** (Plan 36) —
  its `AssociativeAlgebras` package has no quiver / vertex-idempotent type, so the M2
  bridge is single-vertex `kQ/I` graded dimensions plus commutative-example Ext only;
  multi-vertex and every Hochschild quantity stay with QPA + the theory oracles. The
  bridge **refuses those inputs loudly** (multi-vertex and any non-`{graded_dims,
  commutative_ext}` subject raise `QuiverlabError`), never silently narrows scope. The
  live M2 bucket (`-m m2`) skips cleanly without a local Macaulay2 and is a hard
  failure in the dedicated CI job under `QUIVERLAB_REQUIRE_M2=1`.
- The `webapp/` and `docs/gui/` tiers are verified as software (plumbing,
  isolation, artifacts), not as mathematics — they compute nothing themselves.
- **HH with bimodule coefficients + relative HH (Plan 52):**
  - **`engine="fast"` refuses a coefficient loudly** — the GF(p) bar-basis
    accelerator carries only `(m, T, unit)`; the regular bimodule is hard-wired, so
    no coefficient object exists there. The refusal points at `bar` (any field) /
    `cs` (presented). `engine="auto"` never routes fast with a coefficient.
  - **Cyclic homology / Connes `B` with a coefficient is out of v1 scope** — Connes
    `B` assumes the coefficient is `A` itself (it rotates the unit into bar slots);
    a general `M` is P54's BV concern, not shipped here.
  - **Relative HH is `B = kQ₀` only** — the vertex-relative (separable) case, exactly
    the Lindell–Rubio y Degrassi setting, is IN scope (`relative_to="vertices"`); a
    general / non-separable subalgebra `B` is **refused loudly** and is the recorded
    Cibils-relative follow-up.
  - **CS coefficients need a presentation** — a presentation-less structure-constant
    algebra computes coefficients via `bar` only; CS refuses (as today) because it
    cannot form a reduction system.
  - **The minimal `A^e` coefficient path is GF(p) int64 only** (engine-internal, the
    cross-check); off GF(p) the coefficient routes through `bar`/`cs`.
  - **QPA has no HH-with-coefficients / relative-HH surface** (probed live via a
    `NamesGVars()` sweep; `tests/qpa/test_coeff_qpa.py` skips honestly and FAILS if a
    matching verb ever appears). The covering oracles are the internal identities +
    the ported GF(p) bank + the CS/minimal cross-engine web + the literature pins.
  - **The literature pins are BLOCKED-until-transcribed** — the gentle-`HH¹`
    (1811.02211) and radical-square-zero relative-`HH¹` (2411.03080) values are
    `xfail(strict=False)` until a concrete example is transcribed VERBATIM (value +
    equation number); a fetch of the 2411.03080 ar5iv HTML errored out at
    implementation time, so no number was fabricated.
- **Gentle / string subsystem (Plan 46):**
  - For a **string** algebra the Butler–Ringel classification is complete **iff there
    are no bands** (rep-finite). When bands exist the algebra is rep-infinite and
    `enumerate_strings` returns a **length-capped sample** (`status="budget"`), never a
    `status="complete"` list. (The plan's original 2-cycle `kQ/(ab,ba)` "band" example
    is a self-injective Nakayama algebra — rep-**finite**, no band; any algebra on the
    2-cycle quiver is Nakayama. The genuine minimal gentle band algebra is the
    **Kronecker quiver**, band `a·b⁻¹`; both are pinned.)
  - The **AG invariant is a DERIVED invariant, provably NOT complete** — completeness
    needs the graded Opper–Plamondon–Schroll geometric data (out of scope). Never
    claim it separates all derived-equivalence classes. Implemented via Asashiba's
    blossoming form (Nakaoka `arXiv:1811.00775`), reproducing that paper's Example
    2.15 verbatim; the convention-free self-cert is that permitted & forbidden threads
    each partition `Q_1`, with `Σn = 2|Q_0|−|Q_1|` and `Σm = |Q_1|`.
  - **String-τ is engine-arbitrated.** `string_tau`/`string_tau_minus` compute the
    Butler–Ringel hook/cohook combinatorics as the primary method but VERIFY the
    result against the trusted Plan-23/41 engine translate (`Module.tau`) on every
    call; at a few AR-quiver boundary strings (the translate of a non-injective
    projective at a source/sink — a structural jump, not a local arm move) the result
    is completed from the engine translate. It never returns a walk whose module is
    not `is_isomorphic` to the engine translate, and never a guessed convention.
  - **QPA has NO string/band enumeration and NO AG surface** — the QPA crosschecks are
    recognizer-level (`IsGentleAlgebra`/`IsSpecialBiserialAlgebra`) + module-level
    (`decompose` of a sum of string modules); a standing `IsBoundGlobal` probe FAILS
    if QPA ever ships one. **SBStrips / String-Applet are not installed and are not
    oracles here** — the honest oracles are AAG 2008 (literature), our own
    bar/CS/AR engines (cross-engine), and QPA recognizers.
  - **Band modules need the eigenvalue in the field** (loud otherwise). The
    `decompose`-based indecomposability spot-checks carry the `char ≤ dim` caveat, so
    the string/band batteries run over **QQ / GF(32003)** (`char > dim`).
- **Marked-surface subsystem (Plan 48):**
  - **v1 = UNPUNCTURED surfaces with non-empty boundary only.** This is the ABCP/LFS
    regime where every arc-adjacency is clean, there are no self-folded triangles, and
    `Jac(Q(T),W(T))` is gentle (hence finite, certified three ways: the FST arc count,
    P44's finiteness certificate, P46's `is_gentle`). Punctured surfaces, closed
    surfaces, and self-folded configurations **refuse loudly** (`quiver_of`/`jacobian_of`
    name the successor P48.1: puncture potentials + self-folded triangles + the
    once-punctured-torus / Markov quiver). The **once-punctured torus is the pinned
    loud-refusal oracle** — it constructs as a valid `Triangulation` but `quiver_of`
    refuses it.
  - **Flip ↔ mutation is certified at the QUIVER level** — `certify_flip_mutation`
    compares `quiver_of(flip(T,a))` against the exact Fomin–Zelevinsky skew-symmetric
    matrix mutation `μₐ` on every interior arc. Full DWZ **potential** right-equivalence
    under mutation is a named successor, **not attempted**; for the gentle v1 scope the
    quiver-level certificate plus `is_gentle` on both sides is the shipped guarantee.
  - **The angle→arrow orientation is ARBITRATED, not assumed.** The disc oracle fixes it:
    the fan of the `(n+3)`-gon must give the linear `Aₙ` quiver `1→2→…→n`. The naive
    anticlockwise reading `sᵢ→sᵢ₊₁` gave the reversed chain `n→…→1`, so v1 ships the
    flipped convention `sᵢ₊₁→sᵢ` (documented in `qp.py`). The annulus orientation is
    doubly-guarded: the P44 finiteness certificate would refuse a fully-oriented cycle
    (the affine `Ã` quiver is acyclic).
  - **The AG invariant is a DERIVED invariant, NOT complete** (inherited from P46, above)
    — `surface_block` carries it only when `is_gentle` is True; never claim completeness.
  - **QPA has NO surface / triangulation / marked-surface constructor** (`IsBoundGlobal`
    sweep), so the QPA crosschecks are at the resulting **gentle-algebra** level
    (`IsGentleAlgebra`/`IsSpecialBiserialAlgebra`), mirroring P46, with a standing guard
    that FAILS if QPA ever ships one — a no-code surface *input* method is white space
    even in QPA.
  - **The free-form draw-a-surface canvas is deferred** (named successor). Surfaces are an
    **input method**, not a new compute kind: v1 ships three build-time presets (disc fan
    `A₃`, annulus `C(2,2)`, hexagon-with-internal-triangle) and catalogs the surface
    constructors (skipped in the webapp scalar form, the `zoo`/non-scalar precedent); the
    produced gentle algebra flows through **every** existing compute kind (hh, resolutions,
    modules, products, …).
- **CRS-2004 Example 2.20 does not reproduce** (Plan 29): the paper states
  `HH¹ = 0` for its Z₅-cycle monomial example, but the validated bar oracle
  robustly gives `dim HH¹ = 1` (an explicit surviving oriented 5-cycle; both
  orientations, all 2-relation variants, CC and GF(32003)). The test pins the
  VERIFIED value and documents the discrepancy — no literature number is frozen
  against a live engine disagreement.
- **Recognizer batteries (Plan 59):**
  - **R34 is a SEMI-DECISION, not a decision.** `homological_string_test` returns a
    definitive `not_string` (a ≥3-summand middle witnessed — k̄-sound and
    field-independent, since a direct-sum decomposition base-changes up to k̄) OR
    `string_over_ground_field` (a **finite ground field** exhaustive confirm that is
    **INCONCLUSIVE w.r.t. the Suárez-Álvarez k̄ theorem** — over k̄ there may be more
    indecomposables / non-rational classes a finite field cannot see) OR `inconclusive`.
    There is **no bare `string` verdict**; the definitive k̄ "string" answer is the
    decidable syntactic `is_string` (P38), which the test carries alongside. The **only**
    loud raise is refute-side: a ≥3 witness meeting `is_string(A) is True` (one engine is
    a bug). A `string_over_ground_field`-vs-`is_string` mismatch is a recorded
    **k̄-gap** (`kbar_gap_note`), never raised. The alternative "middle-term-candidate"
    confirm route is a **scoped non-goal** — realizability of a fixed `M'` as an extension
    is an existence question over an infinite space, undecidable without algebraic-geometry
    machinery the library does not have.
  - **R34 char-scope.** The refute PIPELINE (`knit_ar_quiver` / `almost_split_sequence` /
    `decompose`) leans on char 0 or char > dim (the trace-form radical + Fitting locality);
    only the terminal split step is char-robust. So refute runs over **QQ** in practice
    (small primes only on tiny algebras where char > dim holds throughout), and the
    finite-field exhaustive confirm runs over **small primes with char > dim of the
    middles**. A char-undecidable middle makes that instance `inconclusive`, never a guess.
    Battery scope (perf): the AR knit is ~3 s/module and grows with dim, so the discriminating
    battery pins `kD₄` + `kD₅` (both live-verified 3-summand meshes) and defers `kE₆`
    (minutes to knit) — the mesh mechanism is identical.
  - **R35 pins only the a-Kronecker HH from the literature.** `HH^• = [1, a²−1, 0, …]`
    (ALS 2020) is the sole transcribable closed-form; ALS defers the general toupie HH
    *dimensions* to prior work and gives bases/Gerstenhaber structure, so all other toupie
    HH values are **cross-engine** (bar ≡ CS where bar survives). **Example 7.4.1 is
    DEFERRED** (the quiver is a figure only and the stated vertex/arrow counts did not
    reconcile in the P29 reconstruction — not pinned).
  - **R35 `sl_a` is char 0 only, at the dimension level.** ALS Thm 6.5 gives
    `HH¹(A) ⊇ sl_a` for `k = ℂ`, with `a` = the number of **direct source→sink arrows**
    (NOT the branch count — verified: the commutative `[2,2]` toupie has 2 branches but
    `HH¹ = 0`, i.e. `a_direct = 0`). `toupie_sl_a_lower_bound` returns `a²−1` and is
    **hard-gated to char 0** (loud off it); the claim is `dim HH¹ ≥ a²−1` (equality on the
    a-Kronecker), NOT a Lie-algebra summand decomposition of `HH¹`.
  - **QPA has NO homological-string / toupie surface** (live `IsBoundGlobal` sweep): the
    R34 crosscheck is a direct-session `IsSpecialBiserialAlgebra` parity, and a standing
    guard FAILS if QPA ever ships one.
- **Skew-gentle algebras (Plan 68 / R32):**
  - **The triple `(Q, I, Sp)` is NOT admissible — the internal object is always the
    split algebra.** `⟨εᵢ² − εᵢ⟩` contains `εᵢ ∈ rad \ rad²`, so a direct presentation
    raises `AdmissibilityError` (VERIFIED). Every invariant is computed on the admissible
    idempotent-split algebra `kQ̂/Î` (Chen §3), which is ISOMORPHIC to the skew-gentle
    algebra. The **dim law** `dim(split) == dim(associated gentle)` (HZZ Lemma 1.5) is a
    **NECESSARY** check that the split relations are right — NOT sufficient (a wrong
    relation set can share the dim). The **SUFFICIENCY** oracle is the AR-quiver /
    support-τ-tilting counts on the mesh example (`1→a→2→b→3, Sp={2}, a·b ∈ I`: dim 9 vs
    the zero-split's 8, `11` AR vertices, `46` τ-tilting pairs — all VERIFIED live, and
    all **engine-self-derived**: no external table pins the mesh `11`/`46`, so they are an
    internal two-engine agreement (AR knit vs τ-tilting exchange graph), stated as such).
  - **The brick-finite ⇔ rep-finite certificate is char ≠ 2 (Garcia–Lavoué Thm 3.1).**
    Over char 2 the rep-finite UPGRADE is WITHHELD (`rep_finite is None`,
    `scope = "char==2 (narrowed)"`); the τ-tilting-finiteness verdict is still computed on
    the char-FREE split model over QQ (the exchange-graph BFS leans on
    `is_isomorphic`/`decompose`, rigorous only char 0 / char > dim M, so it is NEVER run
    over GF(2) — M3), while the split construction / dim law remain characteristic-free
    (Chen). Batteries run over QQ (GF(32003) parity where cheap).
  - **The classification / τ-tilting COUNTS also route through the char-free QQ model
    (M3, the char-2 regression).** `skew_gentle_indecomposables` and `support_tau_tilting`
    (hence the `skew_gentle` block's `num_indecomposables` / `tau_tilting.num_pairs`) are
    presentation invariants and are computed over QQ regardless of the caller's field —
    over GF(2) the P41 AR knit returns `is_complete=False` (`status='error'`, the char-2
    `decompose` caveat), which must NOT be read as "rep-infinite". A regression pins the
    GF(2) block reporting the SAME `6` indecomposables as QQ (not the incomplete loop-free
    `5`-sample), with an explicit `classification.note` naming the char-free QQ model.
  - **The recognizer decides the TRIPLE, not the iso-type.** `is_skew_gentle_triple`
    reduces to gentleness of the associated gentle pair + the structural loop/Sp checks;
    recognizing an arbitrary presented algebra as skew-gentle up to isomorphism is the
    iso-problem and is NOT attempted (the `BrauerGraphAlgebra` honest-scope precedent).
  - **The τ-tilting numeric counts are geometric (HZZ), not Chen.** The record cites "the
    2212.06467 [Chen] example counts", but Chen's checkable content is HOMOLOGICAL
    (selfinjective-iff, `K₀` rank, gl.dim-finite-iff); the numeric τ-tilting counts are
    pinned from the geometric model (HZZ §6 rank + the `k×k → 4` hand case). Attribution
    split, stated here by design — no enumerative Chen example is invented. The two
    τ-tilting `oracle_literature` pins are **arithmetic-shaped** — the Boolean square
    `k×k → 4 = 2×2` (the four generalized dissections) and the presentation-independent
    rank `|R| = |Q₀|+|Sp|` (pinned `= 3` on `1→2, Sp={2}` and `= 8` on the 7-vertex line)
    — not published enumerations; the mesh `46` τ-tilting pairs is **engine-self-derived**
    (internal two-engine agreement with the AR count, no external number exists), labelled
    `oracle_crossengine`, not `oracle_literature`.
  - **The module string layer is a LOOP-FREE census, a documented STRICT SUBSET — NOT a
    complete classification.** `classify` / `skew_gentle_module` enumerate the loop-free
    `A^g`-walks (headline `5` of `6`, mesh `8` of `11`); the missing modules are the
    symmetric / loop-traversal (mixed-eigenvalue) indecomposables (e.g. the projective
    `P_1`), which the theory DOES describe (Garcia–Lavoué Table 1 / clan classification) —
    so the incompleteness is in the implementation's loop-free enumeration, NOT the theory.
    `skew_gentle_indecomposables` (the P41 AR quiver) is the AUTHORITATIVE enumeration and
    the counts above are AR-authoritative; completing the symmetric-string enumeration is a
    DEEPER-ENGINES-BACKLOG item (it also upgrades route 2 of the certificate to a full
    special-band agreement oracle).
  - **The band-census route (route 2) is a ONE-SIDED check, not a full cross-engine
    agreement.** `find_bands(A^g)` may be INCOMPLETE for SPECIAL bands (bands touching
    `Sp`), so a band FOUND certifies rep-INFINITE (sound, may only DOWNGRADE the verdict)
    but "no-bands" does NOT by itself certify rep-finite — the exchange-graph route is
    authoritative (W4). Extending route 2 to special bands is a backlog item.
  - **QPA has NO skew-gentle surface** (live `IsBoundGlobal` sweep): the crosschecks are
    on the split algebra as a plain `kQ/I` (`Dimension` / `IsSelfinjectiveAlgebra`), and a
    standing guard FAILS if QPA ever ships a skew-gentle constructor/recognizer.
  - **Plan-doc deviations (documented in the tests):** the plan's `_rep_infinite_triple`
    (2-cycle with `{ab, ba}` + `Sp={1}`) is actually REP-FINITE (its split is AR-complete
    with 10 indecomposables, verified live) and additionally trips a P45 `exchange_graph`
    `status == 'error'`; the certificate battery uses the Kronecker `1⇉2` (the canonical
    rep-infinite gentle algebra) instead. The plan's exact HZZ §6 seven-vertex arrow
    reconstruction is NOT a valid triple (`is_skew_gentle_triple` False); the rank pin uses
    a valid 7-vertex `Sp={1}` triple (linear `A₇`), which pins the same
    presentation-independent rank `|R| = |Q₀|+|Sp| = 8`.
- **QPA has no native Tor** (probed live): the module-Tor crosscheck computes
  `Ext^n(M, DN)` inside QPA by dimension-shifting through `NthSyzygy` and uses
  the duality identity as the bridge — plus quiverlab's own self-certifying
  duality/balance anchors.
- **The Plan-35 Hochschild product surface has no external oracle** (2026-08-01;
  `A.cup_products`, `A.cap_products`, `A.gerstenhaber_brackets`,
  `A.connes_differentials`). QPA 1.37 exposes **no** Hochschild product surface at
  all — no `CupProduct`, no `HochschildCohomologyRing*` (a live `NamesGVars()`
  sweep finds zero `Hochschild`/`Cup` name; its `ExtAlgebraGenerators`/`YonedaProduct`
  is the *module* Ext algebra `Ext^*_A(M,M)`, a different object from
  `HH^*(A) = Ext^*_{A^e}(A,A)`), so `tests/qpa/test_products_qpa.py` is an honest
  **skip** that FAILS loudly should a future QPA ever grow the surface. The covering
  oracles are therefore internal: the **identity batteries** (graded commutativity,
  associativity, Jacobi/antisymmetry, cup-Leibniz, the cap module law `(z∩f)∩g =
  z∩(f∪g)`, `B²=0`, and SBI rank consistency) and the **literature pins** on
  `k[x]/(x²)` (with the corrected dimension-2 degree-0) and the QuantumCI-BGMS cup
  dims. Two scope facts are binding: (i) **the Gerstenhaber bracket now goes native
  on the Chouhy–Solotar resolution (Plan 51, 2026-08-07)** — the Negron–Witherspoon /
  Volkov homotopy liftings deliver `[-,-]` over **any exact Domain, at any degree
  (past the bar window)** for every quiver-presented algebra, exactly like the
  CS-native cup and cap; over GF(p) in-window the bar/tt route is unchanged (the result
  object records the window) and is the transported anchor. The only honest gap is a
  **presentation-less structure-constants algebra off GF(p)** (the CS route needs a
  presentation, same as cup/cap): it keeps only its in-window GF(p) transported
  bracket, and closing that needs a minimal/Bardzell diagonal (recorded follow-up,
  Plan 75 GHMS / a monomial diagonal). The degree-0 insertion action stays out of
  scope. (ii) The **structure constants are basis-dependent** — they are read on the
  recorded HH basis (bar/GF(p) or the CS class basis), and each product object
  records which basis (`HHProducts.basis`); the cross-engine gate therefore compares
  only basis-independent data (dims and flattened rank), never the raw constants.
- **The native bracket sign is ARBITRATED, not assumed (Plan 51).** The
  `(−1)^{(p-1)(q-1)}` factor (Oke Thm 3.5 / Volkov / classical Gerstenhaber, matching
  the transported route by construction) is FIXED by the in-window native ≡ transported
  anchor over GF(p), which is a genuine cross-engine comparator (homotopy liftings vs
  the bar circle product). It genuinely discriminates the sign: a flipped sign breaks
  the anchor on a nonzero class (`test_flip_breaks_anchor_kx2`). CORRECTION to the
  R1/plan brief: the plan named QuantumCI (2,2) as the odd-exponent discriminator, but
  `HH³(QCI/GF5)=0` (dims `[2,2,1,0,2,4]`) makes the (2,2) bracket the zero class —
  VACUOUS. The genuine non-vacuous odd-exponent anchor is **QuantumCI (2,4) → HH⁵**
  (dim 4, `(p-1)(q-1)=3` odd): the class is nonzero, native ≡ transported, and the
  flipped sign disagrees (`test_odd_exponent_sign_anchor_qci_2_4`, `slow`; needs Δ₅
  ≈ 2 min). Graded antisymmetry is a consistency check (it forces the relative minus
  between the two terms but not the overall exponent). The **Oke arXiv:2103.12331 §7
  Koszul-quiver bracket tables are BLOCKED-until-transcribed** — a strict-`xfail` fence
  (`test_oke_section7_koszul_bracket_tables`) that flips to a real assert only after
  the PDF's quiver/relations/values are transcribed verbatim; the plan never fabricates
  them, and the k[x]/xⁿ + QuantumCI + anchor oracles carry the burden meanwhile.
- **The BV operator `Δ` on `HH^*` (Plan 54, R2):** (a) **GF(p), in-window in v1.**
  The bracket arbiter that certifies correctness is itself GF(p)-window-bounded
  (Plan 35), so v1 serves `Δ` over GF(p) in the bar window and refuses loudly off
  GF(p) / past window (`engine="cs"` raises `"not available until P51"` — the P51
  native past-window enhancer is an optional follow-up, never a silent fallback).
  (b) **The general BIKLZ non-semisimple construction is refused** (`e ≥ 2`
  self-injective Nakayama with non-semisimple ν, `char | ord ν`): the loud
  BIKLZ-blocked refusal. The `e = 1` case `k[x]/(x^N)` is *symmetric* (ν = id), so its
  char-sensitivity + the BIKLZ §3.2 Δ-value oracle are LIVE via the symmetric route
  in every characteristic — nothing about the `e = 1` values is deferred. The BIKLZ
  §3.2 value oracle asserts the char-sensitive dims + the basis-independent Δ-RANK
  profile the formulas force (`rank Δ_odd = N−1`, `rank Δ_even = 0`); the exact
  per-representative Δ entries are representative-dependent (a Plan non-goal) and are
  pinned indirectly-but-rigorously by the cross-engine bracket arbiter. (c) **QPA has
  NO BV surface** — the covering oracle is the theory battery (`Δ²=0` + the seven-term
  relation + the perfect-pairing certificate + the derived-bracket == independent-
  Gerstenhaber arbiter + the Tradler/BIKLZ literature pins). (d) **ν is defined only
  up to inner automorphism**; the gate tests the concrete `nakayama_automorphism()`
  representative, and the bracket arbiter is the per-instance CORRECTNESS GATE — but
  note precisely what it certifies: it pins Δ **modulo cup-derivations** (the data the
  BV relation `[a,b] = ε(Δ(a∪b) − Δa∪b − (−1)^p a∪Δb)` constrains) together with
  `Δ²=0`, **not** every last coordinate (a Δ and a cup-derivation-shifted Δ are
  indistinguishable to it). A wrong ν-representative fails the arbiter and refuses,
  never returns a wrong Δ. The twist DIRECTION (ν vs ν⁻¹) is determined by the arbiter,
  not assumed; on the flagship `QuantumCI(q=2)` the two directions even give different
  twisted-homology dims, so the perfect-pairing dim-match already rejects ν⁻¹ and the
  arbiter confirms ν. (e) **char 2 is SERVED** but its certification is `Δ²=0` + the BV
  relation mod 2 + pairing invertibility, NOT sign-determination — the sign/transpose
  convention is inherited from the odd-prime (3/5/32003) arbiter, applied uniformly.
  (f) **The weakly-symmetric semisimple-ν "exterior class" is REFUSED loudly (engine
  limitation, not a theorem gap).** The semisimple-ν route carries the ORDINARY twisted
  Connes operator `B_σ = s∘N` (the naive norm-then-degeneracy), whose descent to twisted
  homology is a PER-INSTANCE certified property, **NOT** a consequence of semisimplicity:
  it holds only when the paracyclic defect `(1 − T)` (with `T = t_σ^{n+1}` the diagonal
  σ-action) vanishes EXACTLY on the twisted cycle representatives. For a diagonal ν of
  order 3/4 (the `QuantumCI(q)` instances) the reps are fixed by `T` on the nose and the
  descent certifies; but for a weakly-symmetric ν of **order 2 acting by −1 on odd-degree
  factors** — the exterior algebra `ExteriorAlgebra(2) = Λ(k²)` over GF(5)/GF(7), where
  `ν = diag(1,−1,−1,1) = −id` — `B_σ` of a degree-2 cycle is a cycle only MODULO
  boundaries, so the chain-level descent self-cert fires and `bv_operator` refuses
  loudly rather than return a wrong Δ. LZZ (arXiv:1405.5325) guarantees `HH^*(Λ(k²))`
  IS a BV algebra; the fix is the **general LZZ twisted Connes operator** (a per-class
  correction solve `b·w = (1 − T)·z`, `B̃z = Bz − w`), which is **backlogged**
  (DEEPER-ENGINES-BACKLOG). A fix-round bounded probe confirmed that no cheap
  strengthening of `s∘N` (the `(1 − t)`-corrected forms, the sign-normalization variants,
  or the norm over the paracyclic orbit of the correct order `r(n+1)`) recovers an exact
  chain-level descent that ALSO reproduces the independent Gerstenhaber bracket: the only
  variant passing the exterior descent (`s(1−t)N = s(1−T)`) gives the WRONG bracket there
  and zeroes the QuantumCI Δ. **Presentation-dependence (pre-existing `is_symmetric`
  behavior), stated plainly:** the SAME algebra Λ(k²) presented as `QuantumCI(q=−1)` is
  detected `is_symmetric = True` and is SERVED via the symmetric (Tradler) route, whereas
  presented as `ExteriorAlgebra(2)` it is detected non-symmetric and routes to the
  refused semisimple-ν branch — same mathematics, different route/outcome, because
  `is_symmetric` reads the concrete presentation's trace form. Pinned by
  `tests/hochschild/test_bv_twisted.py` (the exterior refusal as a known boundary, and
  the QuantumCI(−1)-symmetric-served vs ExteriorAlgebra(2)-refused contrast).
- **Silting theory (Plan 67) is honestly scoped on five axes.** (i) The **silting
  quiver can be infinite** — kA₂ already is (AI Example 2.45) — and
  mutation-transitivity is proven only for **local / hereditary / canonical** algebras
  (AI Thm 1.2; it FAILS for a symmetric algebra [AGI]). So `bounded_silting_exploration`
  makes **no general BFS / enumeration claim**: it is a bounded-radius walk with a loud
  `status` (`complete`/`radius`/`budget`) and `finite_class ∈ {"local", None}` —
  certified **complete only for local**. There is deliberately **no `"two_term"`
  finite_class**: a general mutation walk leaves the 2-term slice, so 2-term finiteness
  is P45's and is cross-checked *directly* AND **bidirectionally** (Plan-67 fix round, H2):
  the silting engine INDEPENDENTLY enumerates the 2-term slice (a `silting_mutate` BFS in
  the canonical {0,1} window, dedup by `g_proj` fingerprint) and its SET of objects **and
  its edges** equal P45's `exchange_graph` both ways — pinned on kA₂ (5) and the
  non-hereditary self-injective kZ₃/J² (14) — never P45 re-verified against itself, and
  never claimed by the walk. (ii)
  **Generation is three-valued** — `is_silting` is `True` only on the tilting / 2-term
  (IJY) / local classes, `"unknown"` on the K₀-basis-only case, `False` otherwise;
  `det(g_proj) = ±1` in the **projective** K₀ basis `(+)_v Z[P_v]` is NECESSARY (AI Thm
  2.27) but **not sufficient** in general (thick subcategories are not K₀-classified —
  Krah phantom, `arXiv:2302.12502`), so a silent `True` is never emitted off the
  decidable classes. The projective basis is load-bearing: the P43 `is_tilting_complex`
  used the composition-factor χ (`det(Cartan·g_proj)`) and was a systematic
  false-negative on non-unimodular Cartan (self-injective/symmetric; fixed in Task 0,
  `g_proj`). **The same three-valued honesty now governs `is_tilting_complex` itself
  (Plan-67 fix round, H1):** the P43 surface previously returned a hard `is_tilting = True`
  whenever the object was two-sided rigid with `det(g_proj)=±1`, but
  `rigid + (#summands = rk K₀) ⇒ tilting (generation)` is exactly **Rickard's rank
  QUESTION — an OPEN problem** (only partial answers, e.g. Zhang for self-orthogonal
  τ-tilting modules of finite pd). Generation is now CERTIFIED (`generation == "certified"`,
  hard `True`) only where a completion theorem reaches — a **2-term** self-orthogonal
  K₀-basis object is 2-term silting (IJY) and, two-sided rigid, tilting (covering the
  regular object and every APR / 2-term tilt, all backward-compatible) — and a **wide**
  (non-2-term) rigid K₀-basis object is Rickard-open (`"k0_necessary_only"` ⇒
  `is_tilting == "unknown"`). This is the SAME bare-K₀ inference the silting K₀-only rung
  refuses, so the tilting rung and the silting K₀-only rung now agree exactly in the
  Rickard-open regime (a genuine width-2 rigid tilting complex over `k[1→2→3]/(ab)`
  returns `"unknown"` under BOTH verifiers). (iii) The **co-t-structure is a documentation record** (coheart `add(T)` +
  aisle/coaisle descriptor strings + references), NOT a computed subcategory — the aisles
  are infinite; it refuses on non-silting input. (iv) `End(μT)` **vs Oppermann is verified
  at the underlying-quiver level** (the degree-0 (source,target) arrow multiset — Example
  2.47 left-at-vertex-1 → the single-arrow linear A₂); the full graded dg degree
  bookkeeping (Oppermann Thm 1.1's rotation / composition / anti-composition arrows) is
  **`# PIN`'d, out of the engine's scope**. (v) **QPA 1.37 has no silting surface** (a
  live `NamesGVars()` sweep finds no `SiltingObjects`/`SiltingMutation`/`SiltingQuiver`);
  `tests/qpa/test_silting_qpa.py` is an honest **skip** that FAILS loudly should a future
  QPA ever grow one. Derived-discrete recognition is out of scope (no shipped
  recognizer). The covering oracles are therefore the self-cert batteries (positive-window
  presilting, mutant re-verification of PRESILTING + shares-`n−1`, the `μ⁻∘μ⁺` involution
  on kA₂ AND non-hereditary kZ₃/J², the ChainComplex↔PComplex bridge round-trip, the H1
  tilting↔silting three-valued reconciliation) + the AI/Oppermann literature pins + the
  **bidirectional** P45 2-term cross-check (both set directions + edges, incl. the
  non-hereditary kZ₃/J²).
- **The delooping level (Plan 40, Task F) is DEFERRED, not shipped.** Gélinas's
  `dell(M) = inf{ n : Ωⁿ M is a direct summand of Ω^{n+1} N for SOME f.d. `N` }`
  (Gélinas, *Adv. Math.* 394, 2022 — cited as `gelinas_delooping`) has an
  existential quantifier over `N` with no crisp bounded decision procedure from
  `syzygy` + `is_direct_summand` alone (it needs the injective-side cosyzygy tower
  or an a-priori representation-dimension bound, neither of which this plan builds).
  Shipping a heuristic that silently fixes a finite candidate set for `N` could
  return a wrong finite `dell` when the true witness lies outside the probe — a
  house-honesty violation. The one implementable special case (`findim(A) < ∞ ⇒
  dell(A) ≤ findim(A)+1`) is a bound, not the value. The named successor is a
  future C6-extension plan that builds the injective-side cosyzygy tower and the
  summand-membership check, at which point `delooping_level_bound(A, probe_depth)`
  becomes crisply implementable and this deferral flips.
- **`is_gorenstein` is three-valued True/None — never a bare `False`** (Plan 40).
  It is `True` when both the right and left injective dimensions of the regular
  module resolve finite within the bounded engine, and `None` when either is only
  a certified lower bound. A `False` verdict would require a *proof of infinite
  injective dimension* that the syzygy/injective engines never furnish (they only
  certify "resolved / not resolved within depth N"; only a periodicity certificate
  proves infinity, and that is not wired into `is_gorenstein`).
- **The finitistic-dimension UPPER bound degrades honestly to `None` when the
  global dimension is infinite** (Plan 40). The lower bound is always rigorous (a
  finite pd actually found); the upper bound is `gl.dim` when that is exact-finite
  (`findim = gl.dim`), otherwise `None`. The Igusa–Todorov per-module theorem
  `pd M ≤ ψ(ΩM)+1` (Igusa–Todorov 2005) is genuine, but the aggregate
  `ψ(⊕_v ΩS_v)+1` is **not** a certifiable general `findim` bound — a finite
  general upper bound computed from the presentation would resolve the OPEN
  finitistic dimension conjecture — so no folklore number is emitted (a numeric
  upper is always ≥ lower, gated).
- **QPA has no Igusa–Todorov surface** (Plan 40, probed live 2026-08-05): a
  `NamesGVars()` sweep finds no `Igusa`/`Todorov`/`phiDimension`/`psiDimension`
  name, so the φ/ψ functions have no external QPA oracle — their coverage is the
  Task-B literature battery (φ=pd for finite pd, the Barrios–Mata self-injective
  closed form, projective additivity). `tests/qpa/test_homdims_qpa.py`'s IT probe
  **skips that comparison honestly and FAILS should QPA ever ship an IT surface**;
  the same file's `GlobalDimensionOfAlgebra`/`DominantDimensionOfAlgebra`/
  `GorensteinDimensionOfAlgebra` crosschecks ARE live external oracles.
- **The two deep curated examples carry no products** (Plan 35 §5, Task-12
  feasibility probe). The seeded webapp examples `nakayama-kz20-deep` and
  `nakayama-kz24-deep` (dim ≥ 220) omit the entire product surface: the products
  route through the bar/tt calculus (`to_engine` + cochain bases), whose setup alone
  is ~290 s on kZ₂₀ and whose degree-2 cochain basis is 10.5M cells (over
  `max_cells`, forcing the CS route or OOM), so no product finishes the ~120 s
  probe box at any degree — confirmed directly by a **1500 s (25-minute)** in-process
  `cup:0..2` probe on kZ₂₀ that timed out. The bracket now has a CS-native route too
  (Plan 51), but at dim ≥ 220 its diagonal builds are just as far past the probe box,
  and Connes `B` has no CS-diagonal route at all. Every trim and omission is recorded
  per-example in
  `webapp/precomputed/manifest.yaml`; the four tractable examples carry the full
  surface (`tests/webapp/test_curated_reachability.py`).
- **AR-quiver knitting is an honest SEMI-decision** (Plan 41, C3). `knit_ar_quiver`
  (`Algebra.ar_quiver`) is COMPLETE — it closes with `status="complete"` and every
  indecomposable — **iff the algebra is representation-finite** on the knitted
  component; on a wild or large algebra it hits the module/dimension budget and
  refuses LOUDLY with `status="budget"` and `is_complete=False`, never a silently
  truncated "AR quiver". (Rep-finiteness is undecidable in general, so a complete
  closure is the certificate of finite type on that component, and the budget cap is
  the honest non-answer otherwise.) **Self-injective input is refused up front**
  (`status="unsupported"`, `is_complete=False`): the projective-seeded BFS is only
  valid on an algebra with a postprojective slice, and a self-injective algebra has
  none — every indecomposable projective is injective, so `τ⁻` of each seed is `0` and
  the BFS drains immediately (it used to return `status="complete"` while grossly
  undercounting: `k[x]/(x³)` gave 1 vertex vs the true 3, cyclic `kZ₃/rad²` gave 3 vs
  6; the stable AR component of a self-injective algebra is a periodic tube reachable
  only by stable-component knitting, deferred). The individual `almost_split_sequence`
  is certified whenever it returns by the **ARS socle theorem + char-scoped trace-form
  `rad End(M)` + the socle-simplicity dimension certificate** (`dim_k soc = dim_k
  End(M) − dim_k rad End(M)`) — its exact / non-split / indecomposable-ends checks are
  necessary sanity checks, not the arbiter (see the `modules/` row) — and it refuses
  loudly for a projective, decomposable, or (char ≤ dim, over GF(p)) undecidable input;
  the AR batteries run over QQ or GF(32003)/GF(7 with dim < 7) so both the Fitting
  split search and the trace-form locality certificate decide.
- **Representation geometry — five binding scope facts (Plan 49, C8).**
  (a) `rigidity_codim(M) = dim Ext¹(M,M)` is the codimension of the orbit closure in
  `Rep(Q,d)` **only on a hereditary algebra** (Voigt, `Rep` smooth); on a general
  `kQ/I` it is an **upper bound** (the module variety is cut by the relations) — the
  block states which, and never claims equality off hereditary. (b) `is_rigid` (⇒ open
  orbit) and `orbit_dimension` hold over **every exact Domain** (`dim End` + `Ext` are
  exact everywhere). (c) `canonical_decomposition` is **Dynkin only**: hereditary
  Euclidean/wild is DEFERRED to a named successor (the general Schofield /
  Derksen–Weyman recursion — imaginary Schur roots, isotropic multiplicities), and
  non-hereditary `kQ/I` is refused — both loud; the Kronecker `δ=(1,1)` is a refusal
  oracle. It is rigorous only over char 0 or char > dim (it leans on
  `decompose`/`identify_standard`/the AR knit). (d) `degeneration_order` is
  **representation-finite only** — it inherits the Plan-41 AR-knit semi-decision, so a
  representation-infinite (a hereditary non-Dynkin input is caught up front by the P38
  `form_type` type-check, before the — for tame algebras pathologically slow — knit) or
  self-injective algebra returns `is_complete=False` with a loud `status`, never a
  partial poset. (e) **Hall numbers are out of scope** (they are the P3 axis), and
  **QPA has no orbit / canonical-decomposition / degeneration surface** — the
  `tests/qpa/test_geometry_qpa.py` probe is a fail-if-appears trip-wire; what QPA CAN
  corroborate is `dim End(M)` (the orbit-dim factor) via `HomOverAlgebra`.
- **The Plan-42 spectral-sequence engine — five binding scope facts.**
  second line. (d) **QPA has NO spectral-sequence surface** — there is no `qpa`
  spectral-sequence functions); the covering oracles are cross-engine (`HC` / module
  `SpectralSequences` package is NOT scriptable under M2 1.26** — it rides the
  convergence) against M2's `Complexes` homology of the same total complex; the `E_2`
  only `ss_hochschild`** (the `(b,B)` sequence, algebra-only, schema v1); the
  pattern as P39's GUI deferral. The Koszul degeneration statement is arbitrated, not
  forced: for a Koszul algebra the radical-filtration SS of the minimal
  simple-resolution **degenerates at E_2** (the observed provable page — the folklore
  E_2 collapse, pinned on kA₃ and kA₄).
- **`TrivialExtension(A)` is now a certified quiver presentation** (Plan 31; was
  a silent wrong `False`, then a loud refusal). For a presented `A` over QQ or
  GF(p), `T(A)` is returned as a genuine `kQ_T/I_T` — the quiver of `A` plus one
  arrow dual to each corner-homogeneous basis element of the bimodule socle
  `soc_{A^e}A` (direction reversed), with relations extracted algorithmically
  from the ⋉ structure by a length-lex kernel enumeration. It carries **no
  closed-form theorem citation**: the Fernández–Platzeck presentation was not
  obtained to BibTeX precision, so nothing is transcribed. Each instance is
  instead **self-certified** by the dimension identity `dim kQ_T/I_T = 2·dim A`
  (a `QuiverlabError` otherwise) and **QPA-oracled** against the native
  `TrivialExtensionOfQuiverAlgebra` — consistent with this page's
  no-unverified-pins doctrine. `is_symmetric`, `is_weakly_symmetric`,
  `is_frobenius`, and `is_selfinjective` now return `True` on every `T(A)`
  through the unchanged Plan-29 trace-form certifier, and the four former
  `xfail` fences in `tests/invariants/test_symmetric_regression.py` are real
  asserts. A base with no usable path presentation falls back to the unchanged
  ⋉ structure-constants build (honest refusals preserved, doubling as the
  iso-invariance oracle); the per-instance certificate never lets a wrong
  algebra through.
- **The Plan-33 scale batteries** (2026-07-26). Plan 29 pinned the small
  directions; Plan 33 pushes the *same* oracles to scale, each value cited or
  explicitly cross-engine (the honest-scope labels below are binding). Quantum
  CI: the generalized `k⟨x,y⟩/(x^a, y^b, yx − q·xy)` for
  `(a,b) ∈ {(2,4),(3,4),(4,4),(2,5),(5,5)}` over CC, Bergh–Erdmann cohomology
  `[2,2,1,0,…]` verified independent of `(a,b)` and pushed past degree 8, homology
  `[a+b−1, a+b−2, …]` (`qci_hh_oracle`; char-0 branch only — the small-prime
  reshapes need infinite fields we lack, documented). Preprojective algebras
  Π(A₄/A₅/D₄/D₅) (dims 20/35/28/60): the structural pins `dim`,
  `is_selfinjective`, and Loewy length = h−1 (Coxeter numbers 5/6/6/8)
  (`preprojective`, `assem_book`; Erdmann–Snashall for self-injectivity and
  Loewy length). Depth on monomial self-injective algebras: the Bardzell
  resolution of the cyclic Nakayama algebra kZ₂₀/J¹¹ (dim 220) reaches Hochschild
  degree 300 — the depth showcase, guarded by a context-managed recursion-limit
  raise and a degree-300 regression test (`bardzell`) — and the symmetric Brauer
  stars kZ₄/J⁹, kZ₅/J¹¹ (dims 36/55), whose symmetry booleans the Plan-29
  trace-form fix now certifies (`skowronski_yamagata`, `n | (L−1)`). Taft
  algebras Λ₅/Λ₆ = kZ_n/J^n: `HH_• = [n, n−1, n−1, …]`, with the cyclic-homology
  alternation `HC_{2c}=n, HC_{2c+1}=n−1` (`taillefer_taft`; the HC alternation is
  pinned on the small Λ₂/Λ₃ where the mixed complex is feasible over CC / a char-0
  GF(p) proxy). Canonical algebra C(2,2,2,2,2) (dim 19, 7 vertices):
  `HH² = t−3 = 2`, the first ≥2 case (`schremmer_wpl`, attributing Happel),
  cross-linked to the Happel trace identity. Boolean-lattice B₃ incidence algebra
  (dim 27): `HH^{≥1} = 0` at every depth because the order complex is contractible
  (`0̂` and `1̂` present) — nerve vanishing at scale, `HH_0 = 8` = #elements
  (`cibils_incidence`/`redondo_incidence`). Presented trivial extensions
  T(kD₄)/T(kA₅)/T(kA₆) (dims 18/30/42): the four symmetry booleans,
  `C_T = C_A + C_Aᵀ`, and `HH¹ ≠ 0` (`cmrs_split`), per-instance certified. The
  wild m-Kronecker (m = 3, 4): `HH^• = [1, m²−1, 0, 0]` (Happel, `happel_question`)
  with Coxeter polynomial `t² − (m²−2)t + 1` (`lenzing_delapena_spectral`).
  Exterior algebras Λ(k³)/Λ(k⁴) (dims 8/16): Koszul via
  `g_quadratic_certificate`, self-injective, Loewy = n+1 (`priddy`,
  `froberg_koszul`).
- **Redondo–Román 2018 cup-nonvanishing is deferred**: the paper presents
  `HH^n` as combinatorial sets, not integer vectors; without a clean bar
  anchor the exact nonzero products are convention-risky, so the predicate is
  a documented skip, not a pin.
- **Krull–Schmidt decomposition is certificate-bounded** (Plan 30): Fitting
  splits are exact everywhere, but the LOCALITY certificate (End/rad via the
  natural trace form) is rigorous exactly when char 0 or char > dim M — in the
  small-char regime with no split found, `decompose`/`is_indecomposable` raise
  loudly (naming QPA / a larger characteristic) rather than guess. The τ/τ⁻
  result blocks therefore attach the indecomposability certificate only when
  it is certifiable, and omit it honestly otherwise. Oracles: live QPA
  `DecomposeModuleWithMultiplicities`/`IsIndecomposableModule` (dim-vector
  multisets + multiplicities, GF(7)), constructed direct-sum round-trips,
  Krull–Schmidt uniqueness, and τ-additivity.
- **The worked-steps bundle is replayable by construction** (Plan 30): the
  kA₂ golden asserts every differential of the S₁ resolution appears verbatim
  in the `.html`; larger objects render as stated shape+rank elisions, never
  silent omissions; the `trace_steps.html` source itself is a served artifact.
- **The report is the session's complete record, and its presentation never
  hides or fabricates** (2026-07-29, Marco's desktop-app pass —
  `tests/trace/test_report_completeness_m0729.py`,
  `tests/webapp/test_module_blocks_m0729.py`). Self-certifying gates: the page
  contains no `overflow` rule at all (nothing is clipped behind a scrollbar —
  an over-wide matrix is typeset a size down by a shrink-only, integer-valued
  rule); an arrow acting as the exact zero map is *named* rather than printed,
  but never silently dropped; a differential identical to one already shown is
  *referenced*, and an **elided** differential is never matched as a repeat
  (its body was not recorded, so claiming equality would be a fabrication);
  every result block the two runners produce reaches the report, a failed
  computation included. Honesty pins: a homological dimension whose resolution
  did not terminate by the probed depth renders as the certified lower bound
  `pd M > 32`, never a bare `∞`; the Chouhy–Solotar resolution's terms are
  named as projective bimodules `C_n = ⊕_{s∈S_n} A e_{o(s)} ⊗ e_{t(s)} A` from
  the recorded generator corners, and a term whose corners were not recorded
  (the bar resolution over a structure-constants algebra, which is not
  vertex-graded) claims no decomposition at all. Cross-runner: the block shapes
  are asserted identical for `quiverlab.hpc.spec` and its Pyodide twin
  `docs/gui/runner.py`, so the served page and the desktop app cannot disagree
  about the same computation.
- **The report describes the modules, not just their dimension vectors**
  (2026-07-29 second pass). "The modules" section gives each module the
  computation was about — `M`, and `N` when a second module was named — as its
  Loewy series with top and socle plus the exact matrix of every arrow, through
  the same `module_blocks` serializer the no-code panel consumes (so a printed
  module can be typed straight back in). A Krull–Schmidt summand isomorphic to a
  standard indecomposable is NAMED `S_v` / `P_v` / `I_v` and its matrices omitted
  (`modules/hom.py::identify_standard`: dimension-vector prefilter, then the
  exact `is_isomorphic` certificate; an undecidable case leaves the summand
  unnamed and shown in full, never guessed); every other summand carries its full
  action. Section headings name what they hold — Hochschild homology /
  cohomology, Ext, Tor — and the (co)homology table is printed once, not both in
  the computed results and again under a heading that says "Result".
- **Matrices are indexed grids, and artifacts are written as UTF-8** (2026-07-29
  third pass). Every displayed matrix carries an extra header row of column
  indices and header column of row indices over a light-grey rule, so an entry
  is readable by position; entries are copied verbatim and HTML-escaped, and a
  zero-dimensional matrix renders as the symbol `0` rather than an empty box
  (`tests/trace/_matrix_grid.py` reads matrices back out of the rendered page, so
  the renderer tests assert ENTRIES, not a presentation). Separately, a real
  **Windows** defect is now gated: `Path.write_text` defaults to the *locale*
  codec, so the report's em dashes were written as cp1252 bytes and every utf-8
  reader raised — the entire Windows CI matrix failed on it while macOS/Linux
  (whose locale codec is utf-8) stayed green. `tests/trace/test_artifact_encoding.py`
  is a source-level AST scan (no text I/O in `src/quiverlab`, `webapp`, `docs/gui`
  may omit `encoding=`) plus a live round-trip under a forced cp1252 locale.
- **The offline desktop app has no time limit** (2026-07-30, Marco: a user may
  start a real computation and leave the machine overnight). The deployed
  server's 15-minute wall cap and its "too big, use the email tier" refusal are
  DoS protection and cost-gating for a SHARED public service; on the user's own
  laptop neither applies, so the offline config sets `job_wall_seconds = 0` (both
  the parent deadline kill and the child's `RLIMIT_CPU` are disarmed) and lifts
  the queued-tier thresholds so every request the GUI can express is queued and
  run. The MEMORY ceiling stays. `tests/webapp/test_offline_no_time_limit.py`
  pins all of it, including that the DEPLOYED defaults are unchanged, that an
  explicit `QLWEB_*` override still wins, and that quitting the app now FAILS the
  interrupted job rather than requeueing it (with no wall cap, requeue-on-launch
  would restart it forever). It also pins the `or`-vs-`is None` fix: a job row
  carrying an explicit `wall_seconds = 0` used to have that swapped for the config
  cap, because 0 is falsy.
- The Plan-28 container tier: what pytest verifies is the **wheel-side story**
  (CLI ≡ public-API parity, renderer goldens, checkpoint-resume, byte-stable
  runner delegation, asset-file gates) plus the CI image smoke (build → run a
  tiny config → render → text-extract, on every tagged release). **Real
  Apptainer on a real cluster is a manual release-checklist step**; the local
  drac-local emulator exercises only the no-container (venv-fallback)
  orchestration path, and `--mem`/OOM behaviour is validated by the host
  `deepen` memory-guard tests, not by the emulator (which records but does not
  enforce memory).
- **Preprojective and exterior-algebra Hochschild values are cross-engine-only**
  (Plan 33): no published Hochschild table was consulted for the preprojective
  algebras Π(Aₙ)/Π(Dₙ) or the exterior algebras Λ(kⁿ), so their `HH`/`HH^•`
  dimensions are labeled **xeng** — supported by CS ≡ bar agreement in the low
  degrees the bar complex reaches (e.g. Π(A₄): `HH_• = [4,2]`, bar ≡ CS at degree
  1) and by QPA where it computes, but by no literature pin. The **structural**
  pins on the same algebras (dimension, self-injectivity, Loewy length = h−1,
  Koszulity) are theory-pinned, as are the Taft homology and the canonical-algebra
  `HH²`.
- **The genuinely deep Plan-33 computations are cluster-scale** and deferred to the
  `SUBMISSION.md` step-4 list, not run in CI: preprojective HH at scale (Π(D₅),
  dim 60, past the shallow degrees the laptop CS reaches — the D/E-type CS
  reduction system needs a larger Gröbner bound than the default, and the
  constructor's `degree_bound` does not propagate into the CS engine), Λ(kⁿ≥4)
  Hochschild depth, `ext_algebra`/Koszul certification at scale, `decompose`
  ≳ dim 50, and dim-30+ non-monomial HH past ~degree 10.
- **Recorded but not built** (Plan 33, build risk or cost): the (D,A)-stacked
  Example 1.2 and the Cassidy non-Koszul witness (Plan-27 feeders); the
  Π(E₆)/Π(D₆) preprojective builds (`AdmissibilityError` at the tested bounds, the
  certification cost growing past them); the incidence algebra ≅ S²; the Toupie
  figure-only example (never a pin); and Redondo–Román 2018 (also deferred above).
- **Deciding derived equivalence is NOT algorithmic** (Plan 43, C8). The derived
  surface ships **verifiers** and **necessary-condition** invariants, never a
  decider. `is_tilting_complex` DECIDES rigidity (on the exact, honestly-reported
  window outside which hyper-Hom is provably zero) but reports **generation** honest
  three-valued for a *given* candidate — K₀-unimodularity (`det g_proj = ±1`) is
  NECESSARY, certified sufficient only on the 2-term/regular class (IJY), else
  `"unknown"` (Rickard's rank question is open); it does not search for tilting
  complexes; `End(T)` is the Rickard derived-equivalent algebra of that candidate. `derived_fingerprint` /
  `compare_fingerprints` speak in **"distinguished / not distinguished by these
  invariants"** — never "(in)equivalent": equal fingerprints do **not** imply a
  derived equivalence, and the 8-vertex cospectral trees are the standing pinned
  counterexample (equal Coxeter polynomial / Cartan / HH / centre, yet not derived
  equivalent). The Cartan Smith factors are the `GL_n(ℤ)`-equivalence class (a
  necessary condition, coarser than ℤ-congruence — the docstring does not claim
  congruence). Classifying `D^b` indecomposables for wild algebras is out of scope.
- **`τ_{D^b}` refuses loudly at infinite global dimension** (Plan 43, Happel). The
  Serre functor / AR triangles of `D^b(mod A)` exist iff `gl.dim A < ∞`; `tau_Db`
  raises a `QuiverlabError` otherwise (`k[x]/(x²)` is the pinned negative case),
  never returning a bogus complex. The `cyclic_dims` field of the fingerprint is an
  **honest per-field `{error}`** when the generic `(b,B)` mixed complex blows up
  (no CS route for cyclic homology off GF(p), so a ≥4-vertex algebra over CC/QQ hits
  the `max_cells` guard) — captured, never crashed, and skipped in the comparison.
- **The two-algebra derived-compare panel is DEFERRED to a post-v0.2.0 successor**
  (Plan 43). The single-algebra `derived_fingerprint` scalar kind ships now on all
  three tiers (schema v1, both runners byte-identical), and `compare_fingerprints`
  exists at the library level; the side-by-side compare panel needs a second-algebra
  request field (a schema change) and is not built in the v0.2.0 release gate (P50 =
  integration + docs) — the same GUI-deferral discipline as the P39 complex layer.
  See the [v0.2.0 GUI-deferral ledger](#v020-gui-deferral-ledger).
- **Basic-ization / Gabriel recovery is char-scoped and split-only** (Plan 44 C7).
  `primitive_idempotents`/`basic_algebra`/`gabriel_quiver`/`presented_form` rest on the
- **`repetitive_slice` ships certified FINITE slices only** (Plan 44 C7). The full
  repetitive algebra `hat(A)` is infinite-dimensional; only `copies`-block slices are
- **`JacobianAlgebra` refuses Jacobian-infinite inputs** (Plan 44 C7) with
  `NotFiniteDimensionalError` (an under-constrained potential — e.g. two loops with
  of Labardini-Fragoso (whose Jacobian-finiteness is a per-surface theorem) are the
- **The preprojective algebra of a Dynkin quiver is NOT a Jacobian algebra of its
  commutator potential `Σ_a (a a* − a* a)` is identically 0 (`∂_a(a a* − a* a) = a* − a*
  cubic potential exists to produce the quadratic mesh relations. The Jacobian
- **Tilting complement MUTATION is deferred** (Plan 44 C7): `bongartz_completion` ships
  the one-complement (Bongartz) case; iterated mutation of an almost-complete tilting
  module is a named successor. QPA DOES expose a computational tilting/approximation
  surface — `TiltingModule(T, n)` (a `false`/`[true, …]` verdict, not the stored
  `IsTiltingModule` PROPERTY) and `MinimalRight/LeftAddMApproximation` — so, contrary to
  test_tilting_qpa.py`).
- **Quasi-hereditary structure + recollements are QPA white space** (Plan 47). QPA has
  no quasi-hereditary / recollement surface at all, so the oracle class here is theory
  pins (Dlab–Ringel, Ringel, CPS) + internal self-certificates — there is no `qpa`
  agreement row, by construction (a `tests/qpa` battery would be an honest skip; none
  ships).
- **Quasi-heredity is order-dependent** (Plan 47). `is_quasi_hereditary(order)` and the
  standard/costandard modules depend on the chosen vertex order; the GUI/report
  `quasi_hereditary` block reports the NATURAL vertex order only, stated in-block.
- **The characteristic tilting summand count and PRESENTED Ringel duals inherit the
  char 0 / char > dim caveat** (Plan 47). `characteristic_tilting`'s `is_tilting_module`
  self-cert uses the P30 `decompose` summand count, and `ringel_dual`'s `presented_form`
  uses the P44 trace-form basic-ization — both rigorous only over char 0 or char > dim,
  refusing loudly off scope. Everything else — Δ/∇ construction, the quasi-heredity test,
  the Δ-filtration multiplicities, and the entire `Recollement` (corner structure
  constants, `A/AeA`, six functors, adjunction dims) — is char-clean pure linear algebra,
  proved by a `GF(2)` cell in each battery.
- **φdim / ψdim are EXACT only for representation-finite input** (Plan 53 R23a;
  `docs/plans/2026-08-07-plan-53-invariants-phidim-fcy.md`). The ⊕-of-all-indecomposables
  theorem (add-monotonicity) needs the AR knit to CLOSE (`knit_ar_quiver.is_complete`).
  When the knit caps at budget, `phi_dim`/`psi_dim` return a **certified LOWER bound**
  over the discovered prefix ∪ the simples (`exact=False`, `status="budget"`) — never a
  claimed sup. A genuine knit **error** RAISES loudly (a lower bound on a broken knit is
  not certified); only `status="budget"` degrades softly. Self-injective short-circuits to
  the exact `0` (Plan-40 φ≡0) before any knit. **Rep-infinite termination is governed by
  the AR knit's budget semantics** (P53 critic): the knit's inner almost-split loop is not
  hard-step-capped, so a LARGE `budget_modules` on a representation-infinite algebra (e.g.
  the 2-Kronecker) can be slow before the budget trips — pass an explicit SMALL
  `budget_modules` for rep-infinite input (the lower bound is honest regardless). A hard
  inner-loop step cap in `knit_ar_quiver` is a named DEEPER-ENGINES-BACKLOG successor.
- **LIT finitistic certificates: no known decision procedure in general** (Plan 53 R23c) —
  explicitly NOT "proven undecidable". Families 1–3 (self-injective / Iwanaga-Gorenstein /
  finite-φdim) emit certified EXECUTING findim upper bounds; **family 4 (finite one-sided
  `id(A_A)`) ships the `ψ_𝒟(V)+n+1` machinery + `# PIN` and is unit-covered on constructed
  `(𝒟, n)` data, but has NO end-to-end executing test** — its "one-sided-finite, not
  Gorenstein" precondition is not certifiable by the bounded engine, which never PROVES an
  injective dimension infinite. The `finitistic_dimension_bounds` `None → certified` flip
  happens exactly in the applicable families; elsewhere the honest `None` is preserved.
- **Fractional Calabi–Yau is certified at the WEAK-ON-GENERATORS tier** (Plan 53 R24) —
  tier 3 of three (strong = a functor isomorphism `S ≅ Σⁿ`; weak = object-wise on ALL
  objects; weak-on-generators = object-wise on the non-projective simples + the `Ω^i ν^j`
  orbit reps). This is a NECESSARY condition for the weak (hence strong) CY property, never
  a functor isomorphism, and it can UNDER-report BOTH the numerator `m` (defined mod the
  Σ-period, reported as `sigma_period`) AND the denominator `ℓ` (a spurious small `ℓ` may
  pass on the generators while failing on a non-generator object). The `tier` field labels
  every payload; "as functors" appears nowhere in what is certified. The tier-1 functorial
  certificate `Ω^{n+1}_{Aᵉ}(A) ≅ (A^∨)_φ` (Ivanov–Volkov Thm 1.8, via the minimal `A^e`
  engine) is the named successor. Complete only for finite-representation-type
  self-injective algebras (periodicity); a rep-infinite self-injective input caps at
  `status="budget"` (never "not fractionally CY").
- **Ivanov–Volkov Table 1 is BLOCKED-until-transcribed** (Plan 53 R24). The full stable CY
  dimensions of standard/nonstandard self-injective algebras of finite type are keyed by
  Asashiba type `(Δ, f, t)`; mapping a type triple to a quiverlab presentation is
  non-trivial and deferred. The fractional-CY acceptance rests on the two independently
  solid pins — `k[x]/(x^a)` = (1,1) (self-derived from the Serre-functor formula) and
  `Π(Δ)` = (2,1) (Geiß–Leclerc–Schröer) — plus the internal-consistency self-cert (one
  `(m, ℓ)` for all generators) and the `ℓ=1 ≡` Ivanov–Volkov integer-form cross-check.
- **QPA has NO φdim / ψdim / Igusa-Todorov-dimension / stable-Calabi-Yau surface** (Plan
  53). `tests/qpa/test_phidim_fcy_qpa.py` scans `NamesGVars()` and SKIPS the value probe
  honestly (FAILS if QPA ever ships one); the covering oracle for the VALUES is the
  literature battery. What QPA DOES crosscheck live are the prerequisites the values rest
  on — `IsSelfinjectiveAlgebra`, `IsSymmetricAlgebra` (ν=id ⇔ symmetric), and `NthSyzygy`
  of the simples (the `Ω²≅id` fact behind (1,1)).
- **The INTRINSIC fundamental group is refused loudly — not bounded-computable** (Plan 56
  R14). `A.fundamental_group()` returns the PRESENTATION group π₁(Q,I) (walks mod the
  relation-homotopy along a spanning tree) and its abelianization by exact ℤ Smith normal
  form — always computable. `A.intrinsic_fundamental_group()` ALWAYS raises: the intrinsic
  group (the inverse limit over connected gradings, Cibils–Redondo–Solotar arXiv:0906.3069)
  is not bounded-computable in general. The char-p oracle π₁(k[x]/(xᵖ)) = ℤ × C_p is the
  INTRINSIC group; the presentation group of the monomial loop is ℤ (which is what quiverlab
  computes and pins) — the C_p torsion is invisible to any presentation, which is exactly
  why the intrinsic group is refused rather than faked.
- **`is_simply_connected` is three-valued and `None` is honest per Adian–Rabin** (Plan 56
  R14). `True` is emitted ONLY via a decidable sufficient criterion (the underlying graph is
  a tree; or triangular with no bypasses so π₁ is presentation-independent, Le Meur, and a
  Tietze-lite substitution trivialises the presentation; or the R16 separation certificate
  holds). `False` is emitted ONLY via a decidable witness (disconnected / oriented cycle /
  nontrivial π₁^ab of the stored presentation). Otherwise the verdict is `None` — triviality
  of a finitely presented group is undecidable, so a `True` is NEVER manufactured from a
  failed search. The `Hom(π₁,k⁺) ↪ HH¹` bound is asserted for TRIANGULAR algebras only
  (Assem–de la Peña; iso for Schurian, CRS) — the sweep excludes non-triangular zoo members,
  whose π₁^ab is still emitted.
- **The strongly-simply-connected / separation recognizer requires triangular input and is
  char/budget-honest** (Plan 56 R16). `separation_condition` / `is_strongly_simply_connected`
  RAISE loudly only on non-triangular top input. The Krull–Schmidt `decompose` of `rad P_a`
  is rigorous only over char 0 or char > dim (Plan 30); a per-vertex refusal is CAUGHT as
  `undecided_char` and the convex sweep maps it (and any budget overrun) to `verdict = None`
  — NEVER a mid-sweep raise and NEVER a silent `True`. The oracles run over QQ / GF(32003)
  where decompose decides. The separation subquiver convention is the CORRECTED one
  (`Q_a` = delete `a` together with its transitive predecessor closure), so trees come out
  separated — the star/deeper-tree/branching-source discriminators pin it.
- **QPA has NO fundamental-group / simple-connectivity surface** (Plan 56).
  `tests/qpa/test_pi1_qpa.py` sweeps `NamesGVars()` and SKIPS honestly (FAILS if QPA ever
  ships `FundamentalGroup`/`IsSimplyConnected`/`SeparationCondition`); the covering oracles
  are the Plan-56 literature/theory pins + the Hom(π₁,k⁺) ≤ dim HH¹ cross-engine bound.
- **The Tits-form tame/wild verdict — five binding scope facts** (Plan 62 R19;
  `docs/plans/2026-08-07-plan-62-tits-tame-wild.md`). (a) **The verdict layer requires an
  algebraically closed base field.** The tame/wild dichotomy (Drozd) is an
  algebraically-closed notion; the two theorems (Bongartz 1984; Brüstle–de la
  Peña–Skowroński 2011) are over an algebraically closed field. The gate `_is_alg_closed`
  reads the domain's `is_algebraically_closed` flag when present (P61 — `True` only on the
  `CC` working domain; `QQ`, `QQ(i)` and every `GF(p)`/`GF(p^n)` are `False`) and otherwise
  falls back to `characteristic == 0`. The combinatorial Tits form and the representation
  type are FIELD-INDEPENDENT in characteristic 0, so a char-0 verdict is an honest base
  change to the algebraic closure (`ℂ̄`), which fixes the quiver + relations and hence `q_A`;
  a char-0 field that is not flagged closed (`QQ`, `QQ(i)`) is a refusal once the flag is
  present, and the reason string NEVER calls such a field "algebraically closed". Over
  `GF(p)`/`GF(p^n)` the FORM is still computed (field-free) and the verdict is refused. (The
  earlier "`CC` and `QQ` collapse to one indistinguishable char-0 field" wording was wrong:
  they are distinct domains — `CC` → `SympyExactDomain` flagged closed, `QQ` →
  `RationalField` flagged not-closed — and P61's flag distinguishes even the entry-level
  look-alikes `CC` ↔ `QQ(i)` by construction route.) (b) **The verdict is gated on the
  P56 certificate:** the tame/wild axis needs `is_strongly_simply_connected().verdict is
  True` (BdlPS); the rep-finite axis needs `is_simply_connected().verdict is True`
  (Bongartz). P56's three-valued `None` (Adian–Rabin / budget / undecided_char)
  **propagates** to a `None` verdict — never a fabricated type. (c) **The m-Kronecker ladder
  exercises the FORM layer only** — the Kronecker quiver (`m ≥ 2`) is NOT simply connected
  (`π₁^ab = ℤ^{m−1}`), so it is out of the verdict scope (`rep_type = None`,
  `simply_connected = False`); it pins `is_weakly_positive`/`is_weakly_nonnegative`, not the
  algebra verdict. (d) **Weak nonnegativity (the tame axis) is decided by an exact
  positive-SEMIdefinite certificate + the classified hypercritical list**, NOT a guessed box.
  Because `HYPERCRITICAL_COVERAGE` is EMPTY by design (the full non-tree hypercritical
  transcription is not verified here), the ONLY route to a `True`/tame verdict shipped today
  is the exact positive-SEMIdefinite (Euclidean) certificate (a PSD form provably has no
  hypercritical restriction — a genuine completeness certificate); the list-completeness
  route is present but dormant until the coverage set is populated. Outside PSD the verdict is
  honest `None` ("hypercritical list partially transcribed"), never a guessed `True`. A
  `False`/wild is always a FOUND exact witness (`q < 0`). The `≤ 9`-variable cap was REFUTED
  by `T_{2,3,7}` (a 10-variable hypercritical form with the sincere large-entry defect
  `(12,6,8,4,10,9,7,6,4,2)`), so NO universal support/entry cap is claimed; every bound is
  read off the encoded list. Transcription-checked (covered) entries: the minimal wild trees
  `T_{2,3,7}`/`T_{2,4,5}`/`T_{3,3,4}` and the 3-Kronecker (defects verified in exact integer
  arithmetic); the full non-tree hypercritical list is `# PIN`ned/deferred. (e) **QPA has a
  Tits-form FORM surface but NO representation-type VERDICT verb** (verified live, QPA 1.37):
  `TitsUnitFormOfAlgebra` / `IsWeaklyPositiveUnitForm` / `IsWeaklyNonnegativeUnitForm` are
  crosschecked exactly (`tests/qpa/test_tame_wild_qpa.py`, a genuine two-implementation
  oracle for weak positivity/nonnegativity, incl. `T_{2,3,7}`); there is no
  `IsTameAlgebra`/`IsWildAlgebra`/`RepresentationType` verb, so the rep-finite/tame/wild
  VERDICT is theory-oracled (Bongartz; BdlPS) — a fail-if-appears probe guards the claim.
  Distinct from P38's `form_type` (definiteness heuristic) — the two surfaces never share a
  label (`form type (definiteness)` vs `representation type (certified)`).
- **Exceptional sequences (classical + τ-exceptional) — ten binding scope facts** (Plan 65
  R27+R28; `docs/plans/2026-08-07-plan-65-exceptional-sequences.md`). (a) **τ-exceptional
  mutation transitivity is proven only in rank 2** (Buan–Hanson–Marsh arXiv:2402.10301) — so
  the τ-exceptional enumeration goes through the **ordered-sτ-tilt bijection**, NOT a
  mutation-BFS; the count `n!·#sτt` is the completeness certificate (the Buan–Marsh theorem +
  P45's already-tested `exchange_graph`). (b) **The `n!` multiplier of the τ count is the
  Buan–Marsh theorem.** The formula identity `signed_count == n!·#sτt` is SELF-CERT
  (tautological — the code DEFINES `signed_count = n!·len(exchange_graph.vertices)`); the
  genuine cross-check is the **materialised** enumeration (`len(sequences) == signed_count`,
  pairwise-distinct, each recognized) on `A₂`(10)/`A₃`(84)/the non-hereditary `k(1↔2)/rad²`(12),
  and — the M-1 deep pins where completion-multiplicity is real — `kD₄`(1200, reorder-invariant,
  all-module count 162 = classical) and `kA₄`(1008). Corroboration: `A₂ = 10` also equals
  Igusa–Todorov's `n!·#clusters`. (b′) **The reduction algebra `C(U)` is completion-DEPENDENT**
  (H-1): different tau-tilting completions of the same τ-rigid `U` give `End(T_U)/⟨e_U⟩` of
  different dim/quiver (live-proved on kD₄); only the category `J(U)`, hence `#sTt(C(U))`, is
  completion-INDEPENDENT (DIJ) — and the enumeration relies ONLY on `#sTt`. The reduction picks
  the first completion in BFS order and self-certifies `#sTt`-invariance across completions (a
  loud raise otherwise); the returned algebra is one valid presentation, its `#sTt` the
  trustworthy datum. The tower recognizer's H-2 chain check compares successors by this same
  `#sTt` invariant. (c) **The (a)↔(b)
  cross-check is hereditary-only, realized as a COUNT equality** — the number of all-module
  (all-positive) τ-exceptional towers equals the classical CES count (A₂=3, A₃=16), two
  independent enumerations agreeing; the termwise `is_isomorphic` lift of a deep reduction
  tower to ambient A-modules needs the general τ-perpendicular object-lift `F` (the DIJ
  equivalence `mod C(U) → J(U)`), which is the plan's scope boundary — the general case rests
  on the count oracle `n!·#sτt`. (d) **Classical enumeration is Dynkin-only** (rep-infinite
  hereditary ⇒ an infinite braid orbit — a loud `status="budget"` refusal via the instant
  Gabriel gate, never an infinite knit); **τ enumeration is τ-tilting-finite-only** (loud
  otherwise). (e) **Char scope:** batteries over QQ (a GF(p) parity spot-check only where
  dim-vectors are distinct); `char ≤ dim` refuses loudly (the shipped
  `decompose`/`is_isomorphic`/`presented_form` caveat). (f) **`#sτt` (hence the τ count) is
  ORIENTATION-INDEPENDENT** — the generalized Catalan/cluster number of the Dynkin type
  (live-verified `#sτt(A₃)=14` for three orientations, `#sτt(D₄)=50` for the subspace, source
  and mixed stars); the individual sequences differ across orientations (different algebras),
  the count does not. (g) **E-type Dynkin counts are literature-pinned** (Obaid et al.), not
  brute-forced (`E₆ → 41472`, etc.); live-check only `A_{≤5}`/`D₄`. (h) **Citations resolved
  (BibTeX-verified):** `BuanMarsh2021` note `arXiv:1802.01169` (the earlier `2011.02068` was a
  WRONG id — a Coptic-NLP paper), `RingelBraid1994` = Contemp. Math. 171 (1994) 339–352 (venue
  RESOLVED; the metaplan's "CMS Conf. Proc. 14" is Crawley-Boevey's), the distinct
  `arXiv:2211.10428` noted separately (not shipped, not load-bearing). (i)
  **`is_exceptional_module` implements the BRICK criterion** (`end_dim=1 ∧ Ext¹=0`) — M3: this
  is `= exceptional` over algebraically closed `k` and on the Dynkin/QQ battery (every Dynkin
  indec is a rigid brick), a SUFFICIENT test over non-alg-closed `k` (a non-brick exceptional
  module with a larger division-ring `End` returns False — a scope-limited false negative,
  never a false positive). (j) **The P45 `mutate` D₄-star defect** (M1, Task 0): a genuine
  pre-existing defect (a spurious `status="error"` with the correct 50 vertices on the mixed
  star `{1→2,3→1,4→1}`) is fixed here with a regression pin (`#sτt(D₄)=50`, `status="complete"`
  on all three orientations); independently, the τ-side gate refuses ANY non-`"complete"` graph
  loudly (never a count off it). **QPA has no exceptional-SEQUENCE surface** (a fail-if-appears
  `NamesGVars()` guard, `tests/qpa/test_exceptional_qpa.py`), so the sequence enumerators are
  theory/self-cert-oracled; the exceptional-MODULE criterion IS QPA-anchored per indecomposable
  on the kA₃ zoo — `dim End_A(M)=1` via QPA `HomOverAlgebra` (brick) and `dim Ext¹_A(M,M)=0`
  via QPA `ExtOverAlgebra` (rigid), both agreeing with our engine.
- **Split-extension LES + certified arrow removal — five binding scope facts** (Plan 72
  R5+R6; `docs/plans/2026-08-08-plan-72-split-extension-arrow-removal.md`). (a) **The snake
  connecting map `δ` is computed on the Chouhy–Solotar `L^e`-projective resolution, cited as
  the cup with the extension class (CMRS Thm 4.1).** The bar route the plan's Task A2 names is
  intractable at the plan's own pinned depths (`T(kD₄)`'s `δ^4` lands in `HH^5(L,M)` whose bar
  chain group has `>2·10^6` rows); CS is another `L^e`-projective resolution, so the SES snake
  and its ranks are IDENTICAL (validated: every pin reproduced, incl. the `T(kD₄)` boundary
  case). `δ` is realized as the coupling block of the coefficient-split-block-triangular CS
  coboundary and validated by the exactness self-cert + assembled≡direct — NOT by an
  independent cup recomputation. (b) **A literal coefficient-cup recomputation of `δ` via the
  Plan-35 surface is out of scope (DD-A2):** Plan-35's cup is coefficient-`A` only
  (`HH^•(A)⊗HH^•(A)→HH^•(A)`), whereas Thm 4.1's cup pairs *different* coefficient bimodules;
  building that coefficient-cup is a ledger deferral (P80). (c) **The `⊗_B` tensor-power
  decomposition (CMRS Cor 3.2) is deferred and does NOT fully split here (DD-A3):** the card's
  "from `HH^•(B) + H^•(B,M)`" is realized as the `p=0` **leading piece** surfaced from the
  `L`-flanks, pinned as the `n=0` identity + the summand INEQUALITY — a degreewise EQUALITY
  `flank_B==HH^•(B)` would be FALSE (and, as `leading_B := HH^•(B)`, a tautology), because Cor
  3.2's full splitting needs `M=D(B)` one-sided `B`-projective, i.e. `B` self-injective (false
  for `kA_n`/`kD₄`). (d) **Cohomology arrow removal is NOT a clean isomorphism.** The homology
  `HH_{≥2}(A)≅HH_{≥2}(B)` (Thm 3.2) is the only iso claim; cohomology carries the Thm-4.2 Ext
  correction (`coh_correction`, nonzero `HH^2` on P2) — reported separately from the `n=0,1`
  center/disconnection deltas (`coh_low_delta`; `HH_0` is provably invariant, an inert arrow
  never lies on a cycle). (e) **QPA has NO split-extension-LES or arrow-removal HH surface**
  (the P35 precedent — a fail-if-appears `NamesGVars()` guard); the INPUTS are QPA-crosscheckable
  (`TrivialExtensionOfQuiverAlgebra`, `HochschildCohomologyHomology` of `A`/`B`) but the LES/
  reduction verdicts are theory/self-cert/cross-engine oracled. General `B⋉M` (arbitrary
  bimodule `M`) is refused loudly — only the trivial-extension flagship `M=D(B)` (whose `L` is
  presentable via `TrivialExtension`) is in scope (DD-A1).
- **Skew group algebras `A⋊G` + the Ştefan HH decomposition — ten binding scope facts**
  (Plan 74 R8; `docs/plans/2026-08-08-plan-74-skew-group.md`). (a) **The DIRECT
  `HH^•(A⋊G)` (the shipped engine on the constructed algebra) is the trustworthy primary;
  the decomposition is the value-add, cross-checked degreewise against it** — a mismatch is
  a loud bug in the bridge, never a silent number. (b) **Char scope is on the DECOMPOSITION
  + covering ONLY, not the constructor**: the smash structure constants + the direct HH are
  char-agnostic (the build never divides by `|G|`); only Reynolds/Maschke need `char k ∤ |G|`.
  (c) **The presented `kQ'/I'` view is the BASIC algebra of `A⋊G`** (via the shipped
  `presented_form`, Morita-equivalent ⇒ HH-faithful, NOT dimension-faithful — the dim law is
  a statement about the full structure-constant `A⋊G`, which the constructor always returns).
  (d) **The covering-reduction HH transport is v1-deferred** (`is_free_action` detection + the
  FREE-orbit oracle `Z/2`-swap → `k[x]/(x²)` shipped; the general Galois-`G`-covering HH
  reduction is ledgered — the Ştefan route already gives HH within scope); the FLAGSHIP
  (`Z/2` on the dual numbers) is `kG`-idempotent-splitting `k[Z/2]≅k×k`, NOT a free orbit
  (`is_free_action False` — the vertex is fixed). (e) **QPA has NO skew-group-HH-decomposition
  surface** (a `NamesGVars()` fail-if-appears guard, the P35 precedent); the input-level DIRECT
  `HH^•(A⋊G)` is QPA-crosscheckable by building the smash by hand, but the isotypic/`Z(g)`-
  decomposition is theory/self-cert/cross-engine oracled. (f) **The left-vs-right twist convention
  is SOURCE-PINNED from Shepler–Witherspoon, NOT summand-arbitrable** — the summand dims are
  twist-symmetric (`{}_φA_ψ ≅ {}_1A_{φ⁻¹ψ}`, live-verified order-2 AND order-3), so the only
  discriminator is the assembled total + the per-summand split against DIRECT. (g) **Deeper
  degrees / larger `|G|` are budget-capped, direct-oracle-pinned** (the dim-9 `Z/3` caps at
  degree 2, the dim-18 `S₃` at degree 0 — honest `status="budget"`); **non-abelian `G` is IN
  scope** (the `S₃`-on-`k³` pin, dim 18, `HH^0=2`, classes `{e}/3 transp/2 3-cyc` with proper
  centralizers `Z(g)=S₃/Z/2/Z/3`). (h) **`arXiv:1804.02223` is WITHDRAWN and mis-attributed in
  the source card** (Cibils–Marcos, two authors — NOT "CMRS"); it is NOT anchored — the
  decomposition rests on Ştefan 1995 (JPAA 103, the spectral-sequence origin) + Shepler–
  Witherspoon 2012 (J. Algebra 351, the additive decomposition). (i) **Modular `char \| |G|`:
  `stefan_decomposition` refuses loudly** (`"modular"`); the direct `A⋊G` is still computable,
  but its `presented_form`/`is_selfinjective` then describe a GENUINELY DIFFERENT (modular)
  algebra — distinct from the `char ≤ dim` `presented_form` refusal seen live on the dim-8 free
  example over GF(7). (j) **The graded Lie / Gerstenhaber-bracket structure of `HH^•(A⋊G)` is
  OUT of scope** — this plan ships only the decomposition of DIMENSIONS; the bracket/Lie content
  is downstream of P70/P71/P51 running on `A⋊G` as an ordinary input (Shepler–Witherspoon 0911.0938
  is fundamentally about the bracket — the P70 seam). (k) **The `Z(g)`-transport reads TWO
  DIFFERENT BASES, and the arbiter that pins them is non-vacuous by construction**
  (`tests/hochschild/test_skew_group_transport_basis.py`). The bar slots are unit-adapted, but
  the COEFFICIENT slot is not: P52's `bar._coeff_in_unit_basis` transports a bimodule with
  `Bimodule.change_of_basis`, which by contract re-indexes only the A-element slot and leaves
  `M`'s own basis alone — so `M` stays in `A`'s ORIGINAL basis and the coefficient slot takes
  `action.matrix_of(h, A)` UNTRANSPORTED. The two matrices coincide on every local algebra
  (`A` already unit-adapted) and under any vertex-FIXING automorphism, so the distinction is
  invisible to the dual-numbers / `k[x]/(x³)` / arrow-swap batteries; it bites exactly when a
  VERTEX-PERMUTING automorphism acts on a multi-vertex algebra. The convention is therefore
  arbitrated against the chain-map property `δ∘H_n = H_{n+1}∘δ` on the `Z/3` 3-cycle rotation,
  whose `hbar` is neither symmetric nor self-inverse — so `hbar`, `hbarᵀ`, `hbar⁻¹`,
  `(hbar⁻¹)ᵀ` are four DISTINCT matrices and the test demands that EXACTLY ONE of the eight
  (coefficient-slot, bar-block) combinations is a chain map on each side: `(hbar⁻¹)ᵀ` for
  cohomology, `hbar` for homology, both on the untransported slot. A test asserting only that
  the shipped combination works would pass under a symmetric involution; a fixture guard fails
  if `hbar` ever goes symmetric or self-inverse. Observed failure mode of the wrong slot is a
  LOUD refusal (the transported class leaves the cycle span), not a wrong number — but that is
  the observed behaviour on these inputs, NOT a theorem, which is why the end-to-end
  vertex-permuting agreement vs DIRECT is pinned too.
- **The radical filtration of `mod A` + Liu degrees — six binding scope facts** (Plan 57
  R37/R21; `docs/plans/2026-08-07-plan-57-degrees-radical-filtration.md`). (a) **Certified
  only on the representation-finite (knit-complete) domain.** `radical_filtration` /
  `ar_invariants` inherit P41's semi-decision contract exactly (`status ∈
  {complete, budget, error, unsupported}`); the nilpotency index, `rad^∞ = 0`, the degrees
  and every R21 invariant are certified iff the knit closes. (b) **Self-injective input is
  refused** (`status="unsupported"`): `knit_ar_quiver` cannot knit a self-injective algebra
  from its projectives, so cyclic Nakayama `kZ_n/J^ℓ` is out of scope — the Chaio-school
  self-injective closed forms are cited as documentation, not computed. (c) **Rep-infinite
  input yields window-restricted layers, NO verdict**: a budget-exhausted knit returns the
  layers of the discovered subcategory (a lower bound on `rad(mod A)`, loudly labelled),
  with `nilpotency_index=None` and no `rad^∞` claim. (d) **The `(rad^∞)²=0 ⇒ rep-finite`
  battery (CMMS 1994) is a CLASS statement, not a per-instance decider**: `rad^∞` is not
  computable on rep-infinite input (no complete knit), so we NEVER construct a rep-infinite
  instance and "verify CMMS" — what is tested is that every rep-finite zoo algebra has a
  finite nilpotency index (`rad^∞=0`, consistent with CMMS) and the honest discrimination
  (a finite index on rep-finite input, `status="budget"/"unsupported"` on rep-infinite —
  never a false finite index); CMMS-1996 cube-zero and Kerner–Skowroński 1991 are
  documentation-only class pins. (e) **QPA cannot compare** the module-category radical or
  the degrees: `RadicalOfModule` is the Jacobson radical of ONE module (`rad M`), not the
  category radical `rad(X,Y)`, and there is no `LeftDegree`/`RightDegree`/
  `NilpotencyIndexOfRadical` verb — `tests/qpa/test_radical_qpa.py` is a fail-if-appears
  `IsBoundGlobal` probe; the one QPA-checkable slice is the layer-1 off-diagonal
  `dim rad(X,Y) = dim Hom(X,Y)` via `HomOverAlgebra`. (f) **A general functorial mesh
  engine + automatic standardness detection is DEFERRED (a named successor)**: route (i)
  ships only as the `kA_n`/`D_4` `ZA_n` closed form (`_mesh_layer_dim`, pure interval
  combinatorics) — the route-(i) ≡ route-(ii) agreement is asserted only on the standard
  test set, and on a non-standard component (Riedtmann's char-2 examples) the two could
  diverge; route (ii) exact linear algebra is the arbiter throughout.
- **Stratifications beyond the quasi-hereditary case are a named successor** (Plan 47).
  This plan ships the quasi-hereditary highest-weight toolkit; general standardly
  stratified / properly stratified algebras and good-filtration *dimensions* past the
  directed oracles are not shipped.
- **The Plan-45 C4 τ-tilting engine — four binding scope facts.** (a) **The mutation
  BFS is an honest SEMI-decision.** `exchange_graph` (`Algebra.exchange_graph`) is
  COMPLETE — closes with `status="complete"`, n-regular, every support τ-tilting pair —
  **iff the algebra is τ-tilting-finite** (AIR Cor 2.38: the exchange graph is connected);
  on a τ-tilting-infinite algebra (e.g. the 2-Kronecker) it hits the pair budget and
  refuses LOUDLY with `status="budget"`, `is_complete=False`, never a silently truncated
  graph — identical to the AR-quiver loud-cap contract. Every downstream enumeration
  (`bricks`, `semibricks`, `maximal_green_sequences`, `wall_and_chamber_fan`, the four-way
  counts) inherits the same complete-iff honesty and omits its value (never a partial one)
  when capped. (b) **Rigorous over char 0 / char > dim; over small characteristic it
  INHERITS `decompose` / `is_isomorphic`'s refusal at the point a certificate is needed — it
  does NOT refuse proactively.** The BFS, every g-matrix dedup, and every brick/semibrick
  enumeration lean on `decompose` / `is_isomorphic` / the trace-form radical, which are
  rigorous over char 0 or char > dim (Dickson/CIW). There is **no proactive char guard**:
  where every module involved is a brick or splits, the engine computes correctly even over
  char ≤ dim — GF(2) kZ₂/rad² (char 2 ≤ dim 4) agrees with the certified QQ result
  (`tests/modules/test_tau_tilting_torsion.py`). It refuses **loudly** — the same
  `QuiverlabError` `decompose` raises — only when a decomposition/iso certificate is actually
  required and the trace-form radical is unreliable at char ≤ dim; it never returns a silent
  wrong pair set or count. **The batteries run over QQ** (with a GF(32003) cross-check where
  cheap, plus the GF(2) no-proactive-refusal pin). (c) **Bricks decide over the algebraically-closed / char-0
  base:** a brick is `end_dim(B) == 1`, which reads "`End_A(B) = k`" only over an
  algebraically closed base (or char 0 with no proper division-ring endomorphisms); the
  GF(p^n) division-ring caveat (`End(B)` a proper division ring, `dim_k > 1`) is stated
  honestly and the brick/semibrick batteries pin over QQ. (d) **The wall-and-chamber fan is
  drawn for n = 2, 3 only** — n = 2 is the exact angular sweep (no floats; the JS does the
  only fraction→pixel conversion), n = 3 uses the **L1/octahedron projection** (the antipodal
  `(0,A)` cone has coordinate-sum < 0 and projects off any single affine plane, so the naive
  `x+y+z=1` plane is the *positive-sum chart only*); for n > 3 the fan payload is `None`.
  **The TILING certificate is n = 2 only** — the exact angular sweep proves the 2D cones
  partition ℝ² (cover the circle once, no gaps or overlaps). **The n = 3 L1/octahedron
  unfolding is a RENDERING, not a certified tiling**: it is certified only PER CHAMBER (every
  chamber g-matrix is unimodular, det ±1) plus a cheap sanity check that each chamber's
  projected net is a nondegenerate 2D triangle and the net count equals the chamber count
  (`test_n3_l1_unfolding_is_a_sane_rendering`). There is **no** proof that the projected faces
  tile the octahedron net without gaps or overlaps — a full 3D fan-tiling certificate is out
  of scope.
  **QPA CANNOT COMPARE — there is no `qpa` battery for τ-tilting**: QPA 1.37 exposes no
  support-τ-tilting / mutation / g-vector surface, so the covering oracles are the AIR
  four-way count identity + Catalan/n-regularity (literature), the involution / unimodularity
  / n=2 fan-tiling + n=3 per-chamber-unimodular unfolding-sanity / King certificates
  (self-cert), and the pair↔`Gen(M)` + fan-normal cross-checks (cross-engine). The external cross-checks NAMED (not run live) are the
  Demonet–Iyama–Jasso tables and Iyama's `fd-applet`; neither is wired as a live oracle.
- **The Plan-63 wall-and-chamber structure via bricks — six binding scope facts.** (a) **The
  wall `D(B)` is a polyhedral CONE, not a hyperplane normal** — the genuinely new object over
  Plan 45, whose fan records only per-edge normals. `D(B) = {θ : θ·dim B = 0 and θ·dim N ≤ 0
  for every submodule N ⊆ B}` is cut from the shipped exact `stability._submodule_dimvecs`
  reusing King's (submodule, ≤ 0) verbatim, so it is byte-consistent with the tested
  `is_theta_semistable`. A simple brick gives the full hyperplane; a non-simple brick a proper
  face (`D(P₁)` over kA₂ is a RAY). (b) **Certified COMPLETE iff brick-finite ⟺
  τ-tilting-finite** (DIJ), decided by the Plan-45 exchange-graph BFS closing; on a
  brick-infinite algebra (the 2-Kronecker) it returns a **bounded region** with
  `complete=False`, `status="budget"`, a `truncation` note and NO count — each discovered
  chamber/wall is exact, only totality is withheld (the P62 discipline; the discovered
  wall-normal groups are counted under `partial_wall_groups`, `num_walls=None` — never
  presented as a definitive wall count). **A spurious exchange-engine `mutate` error
  (`status="error"`) is NEVER read as τ-tilting-infinite** (fix-round, adjudicated
  2026-08-07): the P45 `mutate` can fail on one edge on some finite algebras (the D₄ star; the
  ROOT is fixed by Plan 65 Task 0 in `mutation.py`, untouched here). Since every support
  τ-tilting pair has exactly n mutations and each edge is rediscovered from its other endpoint,
  an n-REGULAR non-budget-capped graph is provably CLOSED — so wall_chamber recomputes the
  n-regularity certificate (it cannot read `eg.n_regular`, which `mutation.py` gates on
  `is_complete`) and, if it holds, RECOVERS the full complete structure (`complete=True`, exact
  counts, an honest `note` recording the provenance; the maximal-green-sequence count is omitted
  rather than a false 0). If the certificate FAILS (a genuinely non-n-regular error graph) the
  computation is refused LOUDLY as unreliable — never the "τ-tilting-infinite" text, never a
  bounded sub-fan. This is a self-cert oracle (D₄ star recovered to 50/12; the loud path and the
  certificate falsifiably tested). (c) **Walls are
  grouped by brick ISO-CLASS, never by dim-vector** — kZ₂/rad²'s two non-isomorphic (1,1)
  bricks P₁, P₂ stay separate (4 walls, not 3), via `torsion.bricks` + `_edge_brick` +
  `is_isomorphic`. (d) **Rigorous over char 0 / char > dim** (QQ default; the brick /
  `is_isomorphic` char caveat is inherited from Plan 45; off scope the loud `QuiverlabError`
  propagates, surfaced as a clean typed error block, never a wrong wall set); bricks decide
  over the algebraically-closed / char-0 base (`end_dim = 1`; the GF(pⁿ) proper-division-ring
  caveat is honest-scope, inherited). (e) **The drawing is rank ≤ 3 only** — n = 2 the exact
  extreme rays from the inequality system, n = 3 the grouped exchange-edge facet vectors
  (each self-certified to lie in `D(B)`) L1/octahedron-projected; the DIRECT 2D extreme-ray
  computation of an n = 3 wall cone (Fourier–Motzkin) is NOT done — the grouped facets are the
  wall's rays in the complete fan and the `θ ∈ D(B)` self-cert proves membership; rank ≥ 4
  gets the inequality tables (`render = "table"`). (f) **QPA CANNOT COMPARE** — there is no
  wall-and-chamber / stability-fan surface in QPA 1.37 (`tests/qpa/test_wall_chamber_qpa.py`
  is an honest `NamesGVars()` probe that skips and FAILS if one ever ships), and since QPA has
  no support-τ-tilting surface either, `#chambers = #support τ-tilting` inherits no live QPA
  cross-check — the named external checks are FD-Applet / Demonet–Iyama–Jasso tables.
- **HH¹ as a Lie algebra (Plan 70 / R11) — ten binding scope facts**
  (`docs/plans/2026-08-08-plan-70-hh1-lie.md`). (a) **Radical / Levi / sl₂-count / toral rank
  are characteristic-0 ONLY** — a loud `QuiverlabError` over `GF(p)` (or `None` + a `char0_note`
  without `require_char0`), on the honest grounds that Cartan's / Levi's / Weyl's theorems have
  **NO positive-characteristic analogue**, so the classification is unjustified for any
  `char > 0` regardless of `dim` (NOT because `char > dim` — that bound does not even apply to
  `W₁`, where `char = dim = p`). Solvability and nilpotency stay live in every characteristic.
  §4 of the plan records the rich char-p phenomenon this walls off: `W₁ = HH¹(k[x]/(x^p))`, on
  which `κ` is degenerate. (b) **sl₂-count / toral rank carry base-change provenance gated on the
  ARITHMETIC field, not the formal `is_algebraically_closed` flag.** Every quiverlab char-0
  domain computes in exact **QQ** — including the DEFAULT `CC` (flag `True`, arithmetic QQ) — so
  the count may undercount an anisotropic sl₂-form that splits only over `k̄`; the
  `base_change_note` is attached for CC / QQ / `GF(p)` / `GF(pⁿ)` alike (the CC case worded
  "formally closed, exact QQ arithmetic"). **MAJOR fix**: gating on the formal flag would suppress
  the note on the DEFAULT field while the arithmetic is QQ — a silent undercount under a "closed"
  banner; no fabricated "algebraically closed" claim (the P62 precedent). (c) **Toral rank = rank
  of the semisimple Levi factor** — the broader maximal torus of `HH¹` including ad-semisimple
  radical elements (e.g. the `k[x]/(x^n)` grading operator) is out of scope. (d) **The bracket
  structure constants are basis-dependent** (`basis="der_inn"`); cross-engine comparison uses only
  basis-independent data (dim, verdicts, series-dim profiles, sl₂-count) — the Plan-35 rule.
  (e) **Card refinement — the sl₂ example is `kK₂`**, `HH¹ ≅ sl₂` exactly (dim 3); `T(kK₂)` gives
  `k ⋉ sl₂` (dim 4, sl₂-count 1, NOT simple) — the card's "T(Kronecker) ≅ sl₂" is the Levi factor
  of `T(kK₂)`, exactly sl₂ only for `kK₂`. (f) **Card refinement — the Witt dichotomy is sharp at
  `n = char = p`**, not the looser "char | n" (`k[x]/(x⁴)`, `k[x]/(x⁶)` over `GF(2)` are solvable;
  `k[x]/(x⁶)` over `GF(3)` is not — but none is the simple `W₁`). (g) **The Gerstenhaber-bracket
  cross-check is GF(p) in-window** (Plan 35 is GF(p); P51 extends to any Domain in-window; degree
  (1,1) is in-window). The **Der/Inn commutator is the primary, field-general route** and needs no
  bracket engine — the whole point of R11. (h) **QPA has no HH¹-Lie surface**; the oracle is GAP's
  core Lie library via the `[qpa]` libgap backend (an honest skip that FAILS if it disappears).
  (i) **Cost/budget** — the Der solve is `d²` unknowns / `d³` equations, empirically **≈ d^5.4 over
  QQ** (measured: dim 16 = 4.3 s, dim 24 = 36 s; `GF(p)` faster); `DEFAULT_MAXDIM = 48`, `A.dim > 48`
  refuses loudly (`status="budget"` at the block layer), and the estimator sizes `hh1_lie` on
  `A.dim` (its budget is NOT a homological degree) so anything large routes off the instant tier.
  The dim-220 Nakayama webapp examples (`nakayama-kz20/24-deep`) carry **no** `hh1_lie` (`≈d⁵`
  infeasible — the Plan-35 products-omission precedent); every oracle (max dim 8) and the CC/GF(4)
  rows fit trivially. (j) **`sl2_count` / `levi_type` beyond one factor are NOT exercised by the
  pins** (all pinned cases have `dim S ∈ {0, 3}`); the general many-factor simple-ideal
  decomposition is the honest-scope fallback — the sl₂-count reads `dim S / 3` ONLY after certifying
  `S` is a sum of sl₂'s (via the invariant-symmetric-forms simple-ideal count), else `status`
  records the incomplete decomposition (`"levi_incomplete"`, `sl2_count`/`levi_type` left `None`),
  never a guessed type.
- **Formal deformations / L∞ / Maurer–Cartan — five binding scope facts** (Plan 78 R13;
  `docs/plans/2026-08-08-plan-78-linfinity-deformations.md`). (a) **No `ℓ_{≥3}` ships.**
  The guaranteed L∞ deliverable is the DGLA on `C(A)` plus the UNCONDITIONAL `rad² = 0`
  dg-Lie certificate on `B(A)[1]`; the `ℓ₃` feasibility spike is FROZEN and ledgered, and
  `l3_status` says so rather than implying a higher bracket was computed. (b) **The
  obstruction is scanned over every basis direction, never inferred from a diagonal.** A
  basis cocycle CAN self-obstruct (`[α₁,α₁] = 2 ≠ 0`, found live during the plan's own fix
  round, correcting an earlier "diagonal vanishing" claim), so `unobstructed` means
  `[α,α] = 0` for ALL `α`, and `obstruction_witness` is a basis-DEPENDENT exhibit whose
  existence is not. (c) **Cost tracks HH-RICHNESS × resolution size, NOT `A.dim`** —
  live-measured: dim-20 `kZ₁₀/J²` with `HH² = 0` is 0.024 s, dim-4 `QuantumCI(-1)` with
  `HH² = 5` is 43.5 s. So `DEFORM_MAXDIM = 32` (the plan's OWN cap — NOT P70's 48, a
  different cost law) is only a coarse backstop, and the webapp estimator routes off a
  cheap `HH²` pre-probe (`estimator._deformations_dim`, effective size
  `dim × (1 + dim HH²)²`); the pre-probe is worth its seconds precisely because the bracket
  it gates costs minutes. The same inversion is why the two cross-runner twins are split
  across buckets (85 s obstructed pin on `deep`, 0.03 s `rad²=0` twin on `fast`). (d) **The
  `q` convention is `xy + q·yx`** (`families.quantum.QuantumCI`), so `QuantumCI(-1)` is the
  COMMUTATOR — the commutative `k[x,y]/(x²,y²)`, `HH^• = [4,4,5,6]` — while the common
  literature convention `yx = q·xy` would read `q = -1` as ANTI-commuting, which here is
  `QuantumCI(1)` with `HH^• = [2,4,6,8]`. The Plan-78 fixtures spell the relation out
  rather than naming the family, so the pinned algebra cannot be misread. (e) **The
  char-p interpretation is gated, not silently reused** — the deformation-theoretic reading
  of `HH²`/`HH³` is stated over char 0; in char p the block still computes the dimensions
  but the interpretation notes say what is and is not claimed.
- **HH• as a graded Lie module over HH¹ (Plan 71 / R12) — eight binding scope facts**
  (`docs/plans/2026-08-08-plan-71-hh-lie-module.md`). (a) **The weight/torus decomposition is
  characteristic-0 ONLY** (P70's gate — no positive-char Cartan/Weyl; reference P70's rationale, do
  not re-argue): a loud `QuiverlabError` (`require_char0=True`) or `weights = None` + a `char0_note`
  over `char > 0`. The **indecomposable-summand decomposition is NOT tied to that gate** — it follows
  `decompose`'s OWN char guard independently (entry (b)); the field-general `ρ_n` / module-axiom /
  graded-module surface stays live in every characteristic. (b) **The summand decomposition inherits
  `decompose.py`'s char guard** (`char 0 or char > d_n`; a per-degree loud refusal surfaced as a note
  at `char p ≤ d_n` with no `dim End = 1` certificate) — so summands ARE returned at `char p > d_n`
  even while `weights` is `None`; never a guessed decomposition. (c) **The weight torus is a
  DIFFERENT object from P70's `toral_rank` — NET-NEW code P70 scoped OUT** (not a reuse of P70
  scaffolding): a maximal ad-diagonalizable abelian subalgebra of ALL of `HH¹`, including radical
  toral elements like the `k[x]/x^n` grading (MNPRS's "Virasoro degree"). P70 exposes only ints/bools
  (no Cartan/torus matrices) and explicitly scopes this torus out, so P71 computes it from scratch
  (the Cartan by de Graaf Engel/Fitting → ad-semisimple part → radical-toral extension → rational
  simultaneous diagonalization). The exact weight LABELS are **basis-dependent / NON-NORMATIVE**
  (scalar-doubling a torus generator doubles every weight; the torus-normalization provenance is
  recorded — the P70 sl₂-triple precedent); NORMATIVE = distinct-weight count, per-weight dims, and
  the weight PATTERN (sl₂-string `{−c,0,c}` / equal-gap progression). (d) **The action structure
  constants are basis-dependent** (`basis="der_inn/bar"`); cross-engine / cross-run comparison uses
  only basis-independent data (dims, the induced-map dim `B_n` / commutant `End`, weight multisets) —
  the Plan-35 rule. The action constants are NOT shipped in the block (they explode). (e) **The sign
  arbiter is GF(p) in-window** (Plan 35 is GF(p); P51 extends to any Domain in-window; degree `(1,n)`
  is in-window). **The field-general `L_D` = the shipped `(1,n)` bracket by a source-level sign
  identity** (`p = 1` collapses all Koszul signs) — the primary route needs NO bracket engine (the
  whole point of R12), verified cochain-for-cochain over GF(p). (f) **QPA has no Hochschild-Lie-module
  surface** (the Plan-35 products precedent); the substantive oracle is GAP's own MeatAxe on the fed
  `ρ_n` matrices (the multi-summand `gl₂` case; composition factors == indecomposable summands for
  these SEMISIMPLE reductive-Lie modules), plus the internal identity batteries. (g) **Cost/budget** —
  the action is the P70 Der solve (`≈ d^5.4` over QQ, bounded by `DEFAULT_MAXDIM = 48` on `A.dim`)
  PLUS the exponential bar cochain complex in `top`; **for `top ≥ 2` the binding limiter is
  `max_cells`, not `DEFAULT_MAXDIM`** (the `d·(d−1)^{top+1}` bar-cell wall `A.hochschild_cohomology`
  already enforces, reused). Every oracle (max `dim` 5, `top ≤ 4`) fits trivially; the dim-220
  Nakayama webapp examples carry NO `hh_lie_module` (the Plan-35 omission precedent). (h) **The MNPRS
  / larger-ALS decomposition tables at scale are DEFERRED, not invented** — the pinned small cases
  (`kK₂`/3-Kronecker ALS toupies, the `k[x]/x^n` truncated-Witt, `k[x,y]/(x,y)²` `gl₂`) reproduce the
  *structure* (irreducible/uniserial/multi-summand); transcribing a specific published MNPRS
  special-biserial decomposition table verbatim needs the paper in hand and is a `# PIN` follow-up (no
  fabricated numbers — the house rule).
- **The Plan-55 left/right parts — five binding scope facts.** (a) **Representation-finite and
  non-self-injective only.** `left_right_parts` is complete **iff** `A` is rep-finite and not
  self-injective. A self-injective algebra (`k[x]/(xⁿ)`) is refused by the P41 knit with
  `status="unsupported"`. A **provably** rep-infinite algebra is refused up front by a
  shipped-theorem certificate, also `status="unsupported"`: the **2-Kronecker** by Gabriel's
  Dynkin criterion on its hereditary quiver, and **ACLV Example 2.2(c)** (`1⇉2⇉3⇉4` bound by
  rad²=0 — *mathematically ada* yet rep-infinite, and **non-hereditary** so the hereditary
  route misses it) by Gabriel's **separated-quiver** criterion (rad²=0 whose separated quiver
  is not a disjoint union of Dynkin diagrams). Both fire **instantly** — the knit's per-module
  almost-split cost makes even a small `budget_dim` cap take minutes on 2.2(c) (measured ~120s
  at `budget_dim=16`), so a fast sufficient certificate is required, not the budget route. Any
  **other** rep-infinite input — not covered by the two certificates — falls through to the
  `budget_dim`-capped knit (per-module dimension cap, default 64) and returns a loud
  `status="budget"` in bounded time, never a partial atlas. `is_complete`/`status` mirror the
  `ARQuiver` loud cap. (b) **The support-algebra build needs a quiver presentation.** A
  structure-constants-only `A` (`A.quiver is None`) cannot present the induced subquiver, so
  `_support_algebra` raises `QuiverlabError` up front — never a fabricated quiver (the
  GUI/webapp always feed a quiver-presented `A`). (c) **Identification is QQ / char-0
  decisive.** Every step that *locates* a module in the universe (`_index_in_U` for `τ⁻¹X` /
  `P_x` / `I_x`) calls `is_isomorphic`, which is decisive over char 0 but **positive-only and
  RAISES** over large GF(p)/GF(p^n) when it cannot exhibit an isomorphism; on an in-scope
  algebra with two non-isomorphic indecomposables sharing a dimension vector this propagates a
  **loud whole-compute refusal** — so the identification batteries run over **QQ**, and the
  `GF(p)` route is used only where every indecomposable has a distinct dimension vector
  (kA_n, the rad²=0 linear Nakayama), where the dim-vector prefilter never enters the
  positive-only branch. Never a silent wrong part or support. (d) **The "product of tilted
  algebras" property is REPORTED, not certified here.** `A_λ`/`A_ρ` are a product of
  quasi-tilted algebras in general (tilted for ada — ACT [5, §2.3] / ACLV Thm A); P55 ships the
  connected-component factors and pins the per-factor *tiltedness* with a
  `skipif(not hasattr(Algebra, "is_tilted"))`-guarded test (`test_support_components_are_tilted_PIN`)
  that GENUINELY auto-activates into a real assert the moment **P60** ships `Algebra.is_tilted`
  (no manual un-skip). (e) **QPA CANNOT COMPARE the left/right-part surface** — QPA
  1.37 exposes no `L_A`/`R_A`/support-algebra verb (a fail-if-appears `IsBoundGlobal` probe is
  the trip-wire), so the covering oracles are the ACLV literature pins + the predecessor-closure
  / support / duality certificates; the *defining* pd ≤ 1 / id ≤ 1 flags ARE corroborated
  pointwise by QPA `ProjectiveResolution` / `InjDimensionOfModule`.
- **The Plan-61 recognizer ladder — six binding scope facts.** (a) **Representation-finite and
  non-self-injective only** (inherited from P55/P41): `recognizer_ladder` is complete iff the
  P55 atlas is, else a loud `status` and NO rungs — never a partial ladder. **ACLV Example
  2.2(c)** (`1⇉2⇉3⇉4`, rad²=0) is the honest boundary oracle: *mathematically ada* yet refused
  because rep-infinite (ada does NOT imply laura in general). (b) **`laura` is trivially `True`
  in scope** — every rep-finite algebra is laura (finite `ind A` ⇒ finite complement), so the
  useful datum is the reported finite complement; the interesting laura content is rep-infinite,
  out of knit scope. (c) **Identification is QQ / char-0 decisive** — the weakly-shod injective/
  projective location and the ada/quasi-tilted placements call `is_isomorphic`, positive-only
  over large GF(p)/GF(pⁿ); a loud whole-compute refusal propagates, so the batteries run over
  **QQ**. (d) **The ada/HH¹ Theorem-B verdict is emitted only over an algebraically closed
  field** — gated on the new additive `A.domain.is_algebraically_closed` flag (`True` only on
  the CC working domain; there is **no** `Algebra.field`, and QQi shares CC's `SympyExactDomain`
  class, so the flag — never the class — is the sound predicate). Over QQ/GF(p)/GF(pⁿ)/QQi the
  block reports `dim HH¹` + the theorem statement, **no verdict** (Theorem B's hypothesis
  unmet). `dim HH¹` is char-0 field-independent (flat base change), so the number computed on
  the exact working domain IS the ℂ statement; the alg-closed hypothesis controls the theorem's
  *validity*, hence the gate is a declaration flag, not a recomputation. The Theorem-B
  `False` branch has a **REAL pin**: the "square" `kQ/(a·b, c·d)` on `1→2→4`, `1→3→4` (both
  length-2 routes killed) over CC is rep-finite and quasi-tilted (hence ada, gl.dim 2, empty
  complement) yet has `dim HH¹ = 1` and `π₁ = ℤ`, so it is NOT simply connected and Theorem B
  returns `verdict = False` — cross-agreeing with P56's `is_simply_connected = False` (both
  `False`, verified live); the pure `_theorem_b_verdict` unit test keeps the gate logic covered
  across all four `(ada, alg_closed)` corners. (e)
  **`tilted ⇒ quasi-tilted`** is `# PIN`'d for **P60** — a skipped auto-flipping test that
  becomes a real assert the moment P60 ships `Algebra.is_tilted`. (f) **QPA CANNOT COMPARE the
  recognizer surface** — QPA 1.37 exposes no quasi-tilted/shod/weakly-shod/laura/ada verb (a
  fail-if-appears `IsBoundGlobal` probe is the trip-wire); the gl.dim gates + the pointwise
  pd/id ingredients ARE QPA-checked (`GlobalDimensionOfAlgebra` / `InjDimensionOfModule`).

- **The Plan-60 tilted-algebra recognizer — five binding scope facts.** (a)
  **Representation-finite + non-self-injective SEARCH scope; theorem gates extend it.**
  `tilted_check` is a **complete** verdict iff a theorem gate fires (**H** hereditary ⇒ tilted,
  including rep-infinite hereditary like the Kronecker quiver; **S** non-semisimple self-injective
  ⇒ not tilted; **G** `gl.dim ≥ 3`, exact or a certified lower bound, ⇒ not tilted) **or** `A` is
  rep-finite, non-self-injective, and the transversal enumeration finishes within
  `budget_sections`. A **rep-infinite non-hereditary `gl.dim ≤ 2`** algebra is an honest
  `verdict="unknown"` (`status="budget"`), and the **local criterion arXiv:1409.2054 (Thm 2.6 on
  a locally-computed finite cut) is the documented, not-yet-implemented rep-infinite extension
  path** (ledger: `docs/plans/DEEPER-ENGINES-BACKLOG.md`). **Rep-infinite hereditary is `tilted`
  (Gate H) but its module-level slice is OMITTED** — Gate H reads the Dynkin type of `A`'s **own**
  quiver first (finite ADE ⇒ rep-finite ⇒ the projective slice is knit-enumerable; Euclidean
  `~…`/wild ⇒ rep-infinite), so for rep-infinite input the AR knit is **never built** (it would
  not terminate) — the Kronecker quiver returns `tilted`/`~A₁` **instantly**; the postprojective
  section exists but is not knit-enumerable, so the slice is omitted, not faked.
  `is_complete`/`status` mirror the `ARQuiver` loud cap. **A disconnected `A` is refuted up front**
  (`verdict="not_tilted"`, `reason="disconnected"`) — a tilted algebra is a connected End-algebra
  over a connected hereditary algebra (ASS2006); the P55 support surface feeds the recognizer its
  connected components individually, so the guard never fires there.
  (b) **Needs a quiver presentation.** A structure-constants-only `A` (`A.quiver is None`) refuses
  loudly with `QuiverlabError` (the GUI/webapp always feed a quiver-presented `A`). (c) **Char
  scope is load-bearing.** The `Hom(X, τY)=0` / faithfulness / self-injective / hereditary checks
  are **exact over every Domain**, but the Ringel certificate uses `is_tilting_module` (summand
  count via `decompose`) and `presented_form` (Gabriel recovery) — **rigorous only over char 0 or
  char > dim, LOUD otherwise** — and any module identification (`is_isomorphic`) is QQ/char-0
  decisive, positive-only-and-raising over large GF(p). So the verdict batteries run over **QQ**;
  a **`GF(32003)` parity check is kept only for `kAₙ`** (distinct dim-vectors, no `is_isomorphic`
  raise). Never a silent wrong verdict over `char ≤ dim`. (d) **`A = End_H(D(S))` is
  theorem-guaranteed, not certified by a bare iso.** quiverlab has no algebra-iso certifier, so on
  a `tilted` verdict the report states `A = End_H(D(S))` as **Ringel-theorem-guaranteed** and
  reports the checkable invariants (`dim A`, `dim H`, the section-graph type ≡ the recovered
  Gabriel type of `H`, loud on divergence) — never a `≅` the engine did not check. (e) **QPA
  CANNOT COMPARE the tilted-algebra recognizer** — QPA 1.37 exposes no `IsTiltedAlgebra` / slice
  search (a fail-if-appears `IsBoundGlobal` probe is the trip-wire), so the covering oracles are
  the Liu/Happel–Ringel/ASS literature pins + the Ringel-1.9(2) self-certificate; the slice module
  IS confirmed a **tilting module** by QPA's computational `TiltingModule(S, 1)`. **Settled
  correction:** the record's "cluster-tilted `A₃` ⇒ not tilted" and "`kZ₃/J²` ⇒ not tilted" are
  the **same** self-injective algebra (the type-`A₃` non-hereditary cluster-tilted algebra
  `Jac(3-cycle, αβγ)` IS `kZ₃/J²`), refuted by the same Gate S.
- **The Plan-69 persistence/TDA bridge — four binding scope facts** (R33;
  `docs/plans/2026-08-07-plan-69-persistence-tda.md`). (a) **Field scope SPLIT.** The
  `A_n`/zigzag barcode is **field-robust** — interval summands are **bricks** (`dim End = 1`),
  which `decompose` certifies in ANY characteristic (the `dim End = 1 ⇒ local` short-circuit,
  never the char-scoped trace-form fallback), so the barcode is byte-identical over
  `QQ`/`GF(2)`/`GF(3)` (LIVE-VERIFIED). The **commutative-ladder route is char-scoped**
  (`char 0` or `char > dim`): CL indecomposables need not be bricks, so the AR knit + non-brick
  `decompose` lean on the trace-form radical; over `GF(2 ≤ dim)` the knit CATCHES its internal
  `decompose` char-refusal and returns `status="error"` / `is_complete=False`, and `barcode`
  refuses LOUDLY off that flag — never a silently truncated `"complete"` (M6, live-verified
  CL(3)/GF(2) → 8 partial vertices). (b) **`CL(n ≥ 5)` is a loud refusal** (Escolar–Hiraoka:
  representation-infinite), enforced at `CommutativeLadder` construction AND in `barcode`'s CL
  branch — never a partial diagram, never a knit-budget "discovery" of infiniteness. (c) **Float
  filtration values are OUT OF SCOPE** (exact-only house rule): the persistence parameter is the
  discrete vertex index; optional `filtration_values` are exact `int`/`Fraction` display labels;
  a float is refused. `∞` is NEVER emitted — a top-reaching FORWARD bar is `death = n` with
  `essential = True` (`essential` is forward-only, always `False` for zigzag — M2). (d) **Theorem
  AND figure counts literature-confirmed (H3 — the BLOCKING reconciliation is DONE).** The
  Escolar–Hiraoka **rep-finiteness theorem** (`n ≤ 4` finite / `n ≥ 5` infinite, arbitrary
  orientation) is a real `oracle_literature` pin. The **AR-quiver-figure indecomposable counts**
  `CL(2)=11` / `CL(3)=29` are **now `oracle_literature` too** — the plan's BLOCKING figure
  reconciliation (Task 6 Step 1a) was **COMPLETED**: the Escolar–Hiraoka paper (arXiv:1404.7588,
  DCG 55(1) 2016) AR-quiver figures were extracted (`pdftotext -layout` + the rendered pages 46–47
  at 150 dpi) and counted from the RENDERED figures — **Figure 13** = `CL(f)` = the plan's
  equioriented `CL(2)`: `3 + 5 + 3 = 11` indecomposables, all thin; **Figure 14** = `CL(ff)` = the
  plan's equioriented `CL(3)`: `1 + 6 + 11 + 6 + 5 = 29`, with **exactly 2 non-thin** entries (the
  `121/011` and `110/121` dim-vectors) — matching the engine's live **27 thin + 2 non-thin** split.
  The orientation identification is airtight: the equioriented ladder IS the incidence algebra of
  the `[n]×[2]` grid poset (iso-invariant), and the figure dimension vectors `2×2` / `2×3` confirm
  the `4` / `6` vertices. So both counts are **figure-confirmed literature pins**, kept
  ALSO as `oracle_selfcert` (Plan-32 overlap is allowed) since the engine's own knit reproduces
  them. The QPA closer independently confirms the DIAGRAM's `decompose` engine across engines (QPA
  has no whole-AR-quiver vertex-count verb). The **`CL(4)` count is `# PIN`** and its **AR knit is
  EXCLUDED from CI** — the dim-30 CL(4) knit did not complete within 30+ min in-session
  (impractical for the deep bucket), so CL(4) rep-finiteness rests on the theorem + the dim-30
  build certificate, and the knit corroboration runs only at `n = 2, 3`. Igusa–Rock–Todorov is
  cited as the **conceptual bridge**, not a computed oracle (quiverlab is finite/exact, not
  continuous). **QPA CANNOT COMPARE** — QPA 1.37 has no persistence/barcode/commutative-ladder
  surface (a fail-if-appears `IsBoundGlobal` guard); the covering cross-engine oracle is
  `DecomposeModuleWithMultiplicities` parity on `A_n` + the CL(3) diagram decompose closer.

- **Han transport across `B ⊆ A` — four binding scope facts** (Plan 73 R7;
  `docs/plans/2026-08-08-plan-73-han-bounded-extension.md`). (a) **The verdict is a FIVE-ROW
  LADDER because CLMS Thm 3.1 makes distinctions a single boolean would erase.** The first
  design read "tensor-nilpotency ⇒ `HH_*(A) ≅ HH_*(B)`", which MISREADS the theorem, and the
  shipped ladder keeps the legs apart: leg (i) alone gives only
  `H_*(B,A) ↪ H_*(A,A)` — an injection with coefficients in `A`, NOT ordinary Hochschild
  homology (`nilpotent_injection`); legs (i)+(ii) give the ordinary-coefficient injection
  `HH_*(B) ↪ HH_*(A)` (`pd_injection`); only FULL boundedness gives the isomorphism and
  hence the two-way `B ⊨ Han ⟺ A ⊨ Han` (`bounded`, Thm 4.6); and `not_bounded` /
  `undecided` report per-leg status honestly. The Han verdict for `A` is claimed ONLY on
  the `bounded` row. (b) **The self-cert gate matches that asymmetry**: the injection bound
  `dim HH_m(B) ≤ dim HH_m(A)` is checked on every run, but EQUALITY is asserted only under
  `bounded` — where a mismatch raises loudly rather than being reported. (c) **`pd_finite`
  is three-valued and conservative**: `True` only when `pd_{B^e}(A/B)` is CERTIFIED finite,
  otherwise `None` = unknown — never `False` — so an uncertified or infinite-looking `pd`
  falls through to `undecided` instead of letting the ladder claim a row. (d) **The
  tensor-nilpotency decision is a capped SEMI-decision**: the length-index theorem +
  J-interrupter routes decide it outright where they apply, and the cap is reported as
  `undecided` with the reached index, never silently treated as "not nilpotent"; the
  reported `injection_from` index is a certified LOWER BOUND, not a precise degree.
- **Incidence-vs-nerve `HH^*` and the GHMS Koszul engine — four binding scope facts**
  (Plan 75 R9/R10; `docs/plans/2026-08-08-plan-75-incidence-koszul-hh.md`).
  (a) **The fast path REQUIRES poset provenance — it is never inferred.** `HH^n(kP) =
  H^n(Δ(P); k)` is a theorem *about incidence algebras*, so applying it needs to KNOW that
  `A` is one. quiverlab does not attempt to recognize an arbitrary `kQ/I` as an incidence
  algebra (a genuinely harder problem, and not obviously decidable here up to isomorphism):
  `incidence_cohomology` reads the `_poset` stashed by `families.IncidenceAlgebra` and
  **refuses loudly** without it — a plain `kA₂` raises, and the webapp/GUI surface it as a
  clean per-entry error, never a number. There is no `kQ/I → poset` recognizer in v1.
  (b) **The ring isomorphism is DIMENSION-LEVEL in v1.** Gerstenhaber–Schack prove an
  isomorphism of graded RINGS; what is pinned here is the graded-vector-space half (the
  order-complex dims ≡ the CS cup engine's table dims, degreewise) plus the
  Gerstenhaber-algebra identities holding on incidence inputs. The cup **pairing-rank**
  check against the simplicial cup product is **DEFERRED and budget-MEASURED**: it needs a
  poset whose order complex is a closed surface, and both smallest candidates were run live
  and produced nothing — the 7-vertex (Császár) torus face poset over QQ (`dim kP = 168`,
  the CS cup route) ran past a 25-minute box, and `RP²₆` over `GF(2)` (`dim kP = 121`) routes
  through the BAR cup instead, which already costs `~10s` at `dim 12` (the crown `C(3,3)`,
  measured over `GF(2)`/`GF(3)`/`GF(32003)`) and returned no result at `dim 121` — while the
  order-complex route answers BOTH in `0.02s`. (That gap is the point of the fast path, and
  also why the ring check is the expensive half.) **No cup value is fabricated**; the check
  ships as an opt-in test (`QUIVERLAB_INCIDENCE_CUP_RANK=1`,
  `tests/families/test_incidence_ring.py`). There is likewise **no bracket claim** on this
  route. And the two quantities are never conflated: this route computes `HH^*`
  (cohomology of `Δ(P)`), which is a DIFFERENT graded space from `HH_*` (Hochschild
  homology) — the diamond incidence algebra has `HH^* = [1,0,0]` and `HH_* = [4,0,0,0,0]`.
  (c) **GHMS is Koszul-only, gated on the FULL three-valued verdict.** `engine="ghms"`
  refuses on `koszul is False` **naming the `Ext`-algebra obstruction**, and refuses
  separately on `koszul is None` ("Koszulity unknown") — the two are honestly
  distinguished, and the gate is never taken on `g_quadratic` alone (which would mislabel a
  genuinely non-Koszul algebra as merely uncertified). The `None` branch ships with **no
  certified witness pinned**. The sharp negative example is the preprojective `Π(A₃)`,
  whose Ext algebra grows a new generator in degree 3: its Koszul kernels VANISH from
  degree 3 (`[3,4,3,0,0,0]`) while the true minimal ranks are `[3,4,3,3,4,3,3]`, so a GHMS
  resolution built from them would claim `gl.dim 2` for a SELF-INJECTIVE algebra of infinite
  global dimension. `koszul_kernels` itself asserts nothing about Koszulity (the kernels are
  defined for any quadratic algebra) — the gate belongs on the RESOLUTION claim, not on the
  kernels. `Π(A₃)` is the shared **P75/P77 boundary**: refused here as not-Koszul, treated
  there as `(p,q)`-almost-Koszul.
  (d) **QPA cannot compare any of this.** QPA 1.37 has no Hochschild surface, no
  order-complex/nerve or simplicial-(co)homology surface, no comultiplicative bimodule
  resolution, and no `IsKoszul` (only `IsQuadraticIdeal`) — all four verified live by a
  `NamesGVars()` sweep in `tests/qpa/test_incidence_koszul_qpa.py`, which SKIPS with that
  finding and **FAILS if QPA ever ships one of them**. The single live QPA leg is the
  Koszul-**Betti** direction, already wired through the Plan-27 `ExtAlgebraGenerators`
  crosschecks that `dim K_n` must match. Everything else is covered by the theory pins
  (`B₃` / crowns / the `RP²₆` characteristic split) and the cross-engine agreements (order
  complex ≡ CS; GHMS ≡ minimal syzygy ≡ bar/CS).
- **The Plan-64 R26 torsion-lattice congruence surface — five binding scope facts.**
  (a) **Rigorous over char 0 / char > dim (QQ default).** The bricks, semibricks, the
  torsion-class order and every count inherit the Plan-45 char caveat verbatim: decisive over
  char 0 or char > dim (Dickson/CIW); off scope the underlying engine refuses **loudly**
  (never a silent wrong lattice). The batteries pin over **QQ**. (b) **Bricks decide over the
  algebraically-closed / char-0 base** (`end_dim(B) = 1`, `End_A(B) = k`); the
  join-irreducibles ↔ bricks bijection and the forcing order on bricks read correctly only
  then — the `GF(pⁿ)` proper-division-ring caveat (`dim_k End(B) > 1`) is inherited from
  Plan 45. (c) **Certified complete iff τ-tilting-finite, STRICTER than a bounded region.**
  Every P64 invariant (`Con`, the forcing order, `#wide`, the join-irreducibles) is a GLOBAL
  function of the *whole* finite lattice; on a budget-truncated exchange graph the lattice is
  a prefix, not a lattice, so P64 emits **no** lattice invariant — `is_complete=False`,
  `status="budget"`, a `note`, everything `None` — rather than a partial-lattice lie. The
  2-Kronecker (τ-tilting-infinite, DIJ) is the honest-refusal oracle. (d) **QPA CANNOT
  COMPARE — no torsion-lattice / congruence / wide-subcategory surface** (QPA has no
  τ-tilting surface at all); `tests/qpa/test_congruence_qpa.py` is an honest skip that FAILS
  if QPA ever ships one. The external cross-checks NAMED (not run live) are the DIRRT
  (Trans. AMS B 10 (2023)), Barnard–Carroll–Zhu (Alg. Combin. 2 (2019)) and Enomoto
  (arXiv:2201.00595) worked examples; none is wired as a live oracle. (e) **`#wide =
  #torsion` iff `A` is representation-FINITE** (Marks–Šťovíček, Bull. LMS 49 (2017), NOT
  merely for hereditary `A`), so on the rep-finite `kAₙ` the count `= Catalan(n+1)` never
  discriminates — the **poset structure** does (kA₂ ≅ M₃ ≇ N₅ = tors kA₂; kA₃ ≅ NC(A₃) ≇
  Tamari = tors kA₃; the M₃/NC-structure pins are literature-anchored via Ingalls–Thomas
  `wide(kAₙ) ≅ NC(Aₙ)`). The strict-inequality `#wide < #torsion` oracle needs a
  rep-INFINITE τ-tilting-finite algebra (the tame preprojective `Π(D₄)`/`Π(A₅)`), which the
  shipped P45 engine could not compute cheaply at spec time (QQ too slow; GF(31) `mutate`
  errored deep in the BFS) — a **documented honest-scope gap**, DEFERRED to a later
  implementation run, with the M₃/NC(A₃) structural arbiters carrying the wide-vs-torsion
  discrimination in the interim.

- **The Plan-66 R29 τ-cluster morphism category + picture group — six binding scope facts.**
  (a) **Rigorous over char 0 / char > dim (QQ default).** The bricks, wides, the Jasso
  reductions and every count inherit the Plan-45/64/65 char caveat verbatim: decisive over
  char 0 or char > dim (Dickson/CIW); off scope the underlying engine refuses **loudly**
  (never a silent wrong category). The batteries pin over **QQ**. (b) **Bricks decide over the
  algebraically-closed / char-0 base** (`end_dim(B) = 1`); the picture-group generators (one
  per brick) and the atom-vs-commutation split (via `Ext¹` between the two simple bricks) read
  correctly only then — the `GF(pⁿ)` proper-division-ring caveat is inherited. (c) **Certified
  complete iff τ-tilting-finite, STRICTER than a bounded region.** `W(A)` is a *finite* category
  iff `A` is τ-tilting-finite (Buan–Marsh); on a budget-truncated (or `status != "complete"`)
  exchange graph the g-fan is a prefix, not a fan, so P66 emits **no** category invariant
  (`is_complete=False`, `status`, `note`, all `None`) — never a partial-category lie. The
  2-Kronecker is the honest-refusal oracle; `status != "complete"` from `A` itself OR any `C_W`
  reduction is a loud refusal (P65 M1). (d) **The K(π,1) verdict is THEOREM-ANCHORED, not
  computed.** `is_kpi1 = True` only for **Nakayama** (Hanson–Igusa, Comm. Alg. 49 (2021) no. 10)
  and **hereditary Dynkin** (Igusa–Todorov, *Which cluster morphism categories are CAT(0)*,
  arXiv:2203.16679, 2022 — CAT(0) for hereditary finite/tame type with small tubes); otherwise
  `is_kpi1 = None` ("not certified"). The general τ-tilting-finite case is **delicate**
  (arXiv:2203.16679 establishes CAT(0) only for hereditary finite/tame type, NOT the general
  algebra), so it is NEVER claimed; a Gromov flag-link / CAT(0) combinatorial certifier
  (decidable but requiring HI's exact cube/link model transcribed) is a **deferred stretch**.
  (e) **Objects = #wide, NOT #torsion/#semibricks.** `#wide = #torsion ⟺ rep-finite`
  (Marks–Šťovíček), so on the rep-finite batteries all three coincide and the trap is
  invisible — the code reads `wide_subcategories(A).size` (never `semibricks(A)`, which is the
  Asai down-labelling `= #torsion`). A rep-INFINITE τ-tilting-finite `#wide < #torsion` witness
  (`Π(D₄)`) is DEFERRED (inherited from the P64 engine limit). The factorization tie
  `n!·#sτt == P65 signed_count` is a **self-cert** consistency gate (definitional in P65), NOT a
  cross-engine oracle (fix-round ruling 3). (f) **QPA CANNOT COMPARE — no cluster-morphism /
  picture-group / cube-complex / τ-perpendicular surface** (QPA has no τ-tilting surface at
  all); `tests/qpa/test_tau_cluster_qpa.py` is an honest skip that FAILS if QPA ever ships one.
  The external cross-checks NAMED (not run live) are the Igusa–Todorov–Weyman (arXiv:1609.02636),
  Hanson–Igusa and Buan–Marsh (IMRN 2021) worked examples — the kA₂ pentagon, the Nakayama
  K(π,1); the picture-group relation **WORDS** are transcribed from ITW/HI (the slope order is
  orientation-dependent), and P66 pins only the transcription-independent invariants (the
  counts, the type split, the brick-membership of each relation, the abelianization).

### v1.0.0 deferral ledger (P51–P79, reconciled at the P80 release gate)

The v1.0.0 program (P51–P79, the 33 implementation plans over records R1–R37) shipped
every record's **mathematics**. What follows is the complete list of what it did **not**
ship, reconciled here so nothing is a silent gap. Each entry names what is missing, what
IS available instead, and the theorem or engine that would be needed to close it. All of
them are also live entries in `docs/plans/DEEPER-ENGINES-BACKLOG.md`.

**Carried — engine-sized, each needs machinery beyond one plan:**

- **The Tate cup product** (Plan 76 / R3). `TateHochschild.cup` is always `None`. The
  Z-graded ring needs a diagonal approximation on the TWO-SIDED complete complex — the
  Plan-20/21 diagonal rebuilt for the GF(p) corner complete resolution and extended past
  degree 0, which is plan-sized. Shipping only the positive part would have been
  **vacuous**: in degrees ≥ 1 the Tate cochain complex IS the minimal resolution's, so
  "the positive cup matches Plan 35" has no content. `periodicity_degree` is the
  resolution's certified period cited as **Usui's criterion**, never an inverse computed.
- **The multi-Koszul DECISION** (Plan 77 / R36). `multi_koszul_certificate(...)["verdict"]`
  is always `None`. Herscovich's multi-Koszul is not a single definition to evaluate: §3.2
  builds it from Tor/Ext-vanishing on the minimal graded BIMODULE resolution — a
  homological engine, not a recognizer over shipped surfaces. What v1.0.0 reports instead
  is the **generation-degree table** and the **K₂ verdict that multi-Koszul implies**
  (Prop. 3.30, pinned). Its canonical examples (Yang–Mills) are infinite-dimensional and
  outside the f.d. engine regardless.
- **The Ginzburg dg engine** (Plan 79 / R31): `Γ(Q,W)`, `D^b(Γ)`, the direct
  orbit-category `Hom_C`, general non-acyclic `C_{(Q,W)}`, and DWZ potential mutation /
  right-equivalence. v1.0.0 ships the certified ACYCLIC slice on BMRRT's finite
  fundamental domain; `End_C(T) = Jac(Q_T,W_T)` is cited, with the FZ quiver step and the
  Jacobian presentation verified per instance.
- **The char-0 Bardzell `ℓ₃` / higher L∞ brackets** (Plan 78). Needs the full
  homotopy-transfer machinery (an explicit contracting homotopy of Bardzell's complex plus
  the tree-summed transferred bracket, RRB arXiv:2008.08122). v1.0.0 ships the
  unconditional **rad²=0 dg-Lie certificate** (`ℓ_{≥3} ≡ 0`) and the induced-`ℓ₂` ≡
  CS-bracket model-independence **theorem statement** (cited, explicitly NOT a computed
  self-certificate).
- **The P72 split-extension deferrals DD-A2 / DD-A3.** A literal coefficient-cup
  recomputation of the CMRS connecting map (Plan 35's cup is coefficient-`A` only, so
  pairing different coefficient bimodules needs a new coefficient-cup surface; `δ` today
  is the CS-Hom-complex snake, cited as the cup and validated by exactness plus
  assembled ≡ direct), and a standalone tensor-power + bimodule-`Ext` engine for the CMRS
  graded decomposition.
- **The general LZZ twisted Connes operator `B_σ`** (Plan 54): the exterior class
  (`ν` of order 2) still refuses loudly rather than returning a BV structure it cannot
  certify.

**Carried — GUI surfaces whose mathematics is already reachable:**

- **The `deformations` `A_α`-adopt flow** (Plan 78) — the deformed algebra ships
  **display-only**; loading it back onto the canvas as a fresh input is a new schema
  surface. The metaplan flagged this one as needing Marco's sign-off to carry rather than
  build: **signed off 2026-08-17 — carried, documented.**
- **The P52 explicit two-sided bimodule matrix editor.** The library and the server both
  ACCEPT the explicit `{dim, left_maps, right_maps}` coefficient form; only the canvas
  editor is missing, because the Plan-26 module editor is one-sided (a right action per
  arrow) while a bimodule needs both actions per generator. The builtin named bimodules
  (regular / dual / twisted / `A/soc`) are a pick-list today. The backlog had provisionally
  assigned this editor to P80; **Marco's release-gate decision (2026-08-17) is to carry it
  documented**, since the mathematics is already reachable through the pick-list builtins
  and the explicit API form.
- **`engine="ghms"` is an ENGINE OPTION, not a GUI kind** (Plan 75) — deliberate, not an
  omission: it adds no new user-visible computation, its HH dims being byte-identical to
  the shipped routes. Its value is as a third independent oracle class.

**Open defects carried with their blast radius measured:**

- **The bibliography parser drops the LAST field of a `.bib` entry** (found at Plan 77).
  27 of 187 entries measured; Plan 77 and Plan 79 fixed their own by adding a trailing
  comma, leaving **20 pre-existing** (e.g. `Usui2021tate` loses its arXiv id,
  `SuarezAlvarez2023` its doi). Not fixed in-band because repairing the parser re-formats
  those citation strings, and result blocks EMBED citation payloads, so it churns every
  frozen runner golden that cites one — a standalone change with its own gate.
- **The P45 `mutate` root cause** (owned by Plan 65 Task 0) is still open: Plan 79's
  Task-0 probe MEASURED that `D₄` reports `status="error"` even post-P65, so the Plan-63
  n-regularity recovery remains load-bearing rather than transitional. Both `D₄` (50) and
  `D₅` (182) are certified only through that recovery.

### v0.2.0 GUI-deferral ledger

v0.2.0 ships the whole C1–C8 mathematics, but a handful of no-code GUI *surfaces*
(and one whole research axis) are deliberately held back to named successors — the
mathematics is reachable now via the library / HPC-config tiers, only the point-and-click
front is deferred. This is the same discipline as the P39 complex layer: every plan ships
a no-code story, and the rest is a named successor, never a silent gap. Each entry names
its plan-doc pointer, and the last entry mirrors the metaplan §8 backlog ledger exactly.

- **Spectral-sequence GUI presets beyond `ss_hochschild`** (Plan 42;
  `docs/plans/2026-08-05-plan-42-spectral-sequences.md`). The engine ships all four
  presets as a library surface — `hochschild_bB_ss`, `radical_filtration_ss`,
  `cartan_eilenberg_ss`, and the `grothendieck_double_complex` builder it wraps — and the
  only no-code compute kind is `ss_hochschild` (the `(b, B)` sequence, algebra-only,
  schema v1). The Cartan–Eilenberg / Grothendieck / radical presets are **API +
  HPC-config accessible** this release; their no-code GUI needs new request fields (a
  second module/algebra plus a preset selector — a schema change) and is deferred to a
  post-release successor. Independently, the *general* Grothendieck sequence (an arbitrary
  `(B, A)`-bimodule via Eilenberg–Watts) is deferred even at the library level: only the
  `U = B` change-of-rings / Cartan–Eilenberg specialization is implemented, and
  `grothendieck_double_complex` refuses `U is not B` loudly.
- **The two-algebra derived-compare panel** (Plan 43;
  `docs/plans/2026-08-05-plan-43-derived-category.md`). The single-algebra
  `derived_fingerprint` scalar kind ships now on all three tiers, and
  `compare_fingerprints` exists at the library level; the side-by-side compare panel needs
  a second-algebra request field (a schema change) and is not built in the v0.2.0 release
  gate (P50 = integration + docs) — deferred to a post-release successor.
- **The free-form draw-a-surface canvas, plus punctures / self-folded triangles and DWZ
  potential right-equivalence** (Plan 48 → successor **P48.1**;
  `docs/plans/2026-08-05-plan-48-surfaces.md`). Surfaces are a no-code *input* method: v1
  ships three build-time presets (disc fan `A₃`, annulus `C(2,2)`, hexagon-with-internal-
  triangle) and the produced gentle algebra flows through every existing compute kind. The
  free-form "draw a surface, triangulate on the canvas" flagship is deferred. On the
  mathematics side, punctured surfaces, closed surfaces, and self-folded triangles refuse
  loudly (P48.1 = puncture potentials + self-folded triangles + the once-punctured-torus /
  Markov quiver), and `certify_flip_mutation` certifies flip ↔ mutation at the quiver
  (Fomin–Zelevinsky skew-symmetric matrix) level only — full DWZ potential
  right-equivalence under mutation is also P48.1.
- **σ_A / τ-Hochschild machinery** (deferred by Marco's explicit choice, not by
  feasibility, to `docs/plans/DEEPER-ENGINES-BACKLOG.md` Tier 2; see the metaplan §8
  ledger). The classical Coxeter matrix / polynomial itself ships (Plan 38, exact
  Cartan-derived), but the Tamarkin–Tsygan-calculus automorphism σ_A of
  arXiv:2606.15595 (its per-degree matrix on `HH_•`, the Thm B/C verification), the
  per-HH-degree "higher Coxeter polynomials" (literature white space — nobody has defined
  them), and τ-Hochschild (co)homology of arXiv:2607.10913 are all held for a future
  release; Chen–Ruan–Yang arXiv:2509.12984 is recorded there as a candidate external
  oracle pending a human read.

---

## References

The literature these oracles test against. Entries with a **registry key** are
rendered from the single packaged `src/quiverlab/citations/references.bib` and also
appear, grouped and annotated, on the [References page](bibliography.md) (cite them
in code via `quiverlab.bibliography(...)`). Entries **without** a registry key are
cited only from test comments or the read-only bank's attributions, at the
precision the repository verifies — no bibliographic detail is invented, and no
theorem number is asserted unless it is actually recorded.

**In the citations registry** (key → work):

- `bar` — Hochschild, G. (1945). On the cohomology groups of an associative
  algebra. *Annals of Mathematics* 46, 58–67.
- `bardzell` — Bardzell, M. J. (1997). The alternating syzygy behavior of monomial
  algebras. *Journal of Algebra* 188, 69–89.
- `chouhy_solotar` — Chouhy, S.; Solotar, A. (2015). Projective resolutions of
  associative algebras and ambiguities. *Journal of Algebra* 432, 22–61.
  arXiv:1406.2300.
- `cup`, `bracket`, `gerstenhaber` — Gerstenhaber, M. (1963). The cohomology
  structure of an associative ring. *Annals of Mathematics* (2) 78, 267–288. (The
  associative cup product and the graded Lie bracket that together make `HH^•` a
  Gerstenhaber algebra — the definitional source for the Plan-35 product surface.)
- `happel_question` — Happel, D. (1989). Hochschild cohomology of
  finite-dimensional algebras. *Lecture Notes in Mathematics* 1404, 108–126.
- `happel_trivial_extension` — Happel, D. (1988). *Triangulated Categories in
  the Representation Theory of Finite Dimensional Algebras.* London
  Mathematical Society Lecture Note Series 119, Cambridge University Press.
  (The trivial extension `T(A) = A ⋉ DA` is symmetric for every
  finite-dimensional `A`; the repetitive-algebra framework.)
- `happel_triangulated` — Happel, D. (1988). *Triangulated Categories in the
  Representation Theory of Finite Dimensional Algebras.* London Mathematical
  Society Lecture Note Series 119, Cambridge University Press. (The Serre functor
  / AR triangles of `D^b(mod A)` exist iff `gl.dim < ∞`, and `τ_{D^b} = ν[−1]` —
  the ground truth for the Plan-43 derived surface; a clean derived alias of the
  same book, distinct from the trivial-extension use above.)
- `rickard_derived` — Rickard, J. (1989). Morita theory for derived categories.
  *Journal of the London Mathematical Society* (2) 39, 436–456. (Derived-equivalent
  algebras share Hochschild/cyclic homology and the centre; `End(T)` of a tilting
  complex `T` is the derived-equivalent algebra — the Plan-43 tilting/fingerprint
  ground truth.)
- `quantum_ci` — Buchweitz, R.-O.; Green, E. L.; Madsen, D.; Solberg, Ø. (2005).
  Finite Hochschild cohomology without finite global dimension. *Mathematical
  Research Letters* 12, 805–816. arXiv:math/0407108.
- `qci_hh_oracle` — Bergh, P. A.; Erdmann, K. (2008). Homology and cohomology of
  quantum complete intersections. *Algebra & Number Theory* 2, 501–522.
- `tensor_product` — Cartan, H.; Eilenberg, S. (1956). *Homological Algebra.*
  Princeton University Press. (The Künneth formula for Hochschild (co)homology.)
- `cyclic` — Connes, A. (1985). Non-commutative differential geometry.
  *Publications Mathématiques de l'IHÉS* 62, 41–144.
- `minimal_resolution`, `module_ext` — Green, E. L.; Solberg, Ø.; Zacharia, D.
  (2001). Minimal projective resolutions. *Transactions of the American
  Mathematical Society* 353, 2915–2939.
- `assem_book`, `nakayama`, `path_algebra` — Assem, I.; Simson, D.; Skowroński, A.
  (2006). *Elements of the Representation Theory of Associative Algebras, Vol. 1.*
  Cambridge University Press.
- `han_conjecture` — Han, Y. (2006). Hochschild (co)homology dimension. *Journal of
  the London Mathematical Society* 73, 657–668. arXiv:math/0408402.
- `qpa` (software; in `references.bib`, no registry key) — Green, E. L.; Solberg,
  Ø. *QPA — Quivers, path algebras and representations*, a GAP package.

**Cited in tests, no registry key** (verified only at the precision shown):

- Loday, J.-L. *Cyclic Homology.* — the classical `k[x]/(x^n)` Hochschild homology
  values and, at **Theorem 2.1.5**, the Connes λ-complex model
  (`tests/invariants/test_cyclic_generic.py`,
  `tests/resolutions_cs/test_battery_literature.py`). The tests name author, title,
  and (for the λ-complex) theorem number; no publication year is asserted here.
- "BACH" — named in `test_battery_literature.py` as the source of the
  truncated-polynomial computations; the repository gives only this token, so
  nothing further is claimed. Those `k[x]/(x^n)` values are additionally attributed
  to Happel (1989) in `test_qpa_reference_validation.py` and cross-confirmed live by
  the bar oracle.
- Assem, I.; Skowroński, A. (1987). *Algebra i Analiz.* — the gentle-algebra notion
  used to build the self-injective 2-cycle pin (`test_battery_literature.py`); cited
  at author, venue, and year, exactly as the test states.
