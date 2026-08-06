# How quiverlab is verified

Every shipped feature of quiverlab is unit tested. This page says **how** — with
the highest rigour we can bring to it — and it is honest about the edges: where a
check is a cross-engine agreement, where it is a published number, where a live
external oracle can reach, and where it cannot.

The suite is **3432 tests** (collected with the `[dev,fast,docs,web,qpa,hpc]` extras,
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
  (`tests/hochschild/test_products_identities.py`, `test_connes_b.py`).
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

The live QPA suite is `-m qpa` (189 tests). GAP is heavy to install, so it runs in a
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
| Cyclic homology | Connes λ-complex second model + mixed-complex identities |
| The Chouhy–Solotar resolution | CS ≡ bar, CS ≡ Bardzell, and the bank byte-level closed forms |
| Deep degrees past the bar window | bank closed forms + cross-engine + closed-form/chain-count pins |
| Frobenius / Nakayama / symmetry | self-certifying `λ`/`ν` identities + socle criterion |
| `HH` over CC and `GF(p^n)` | exact bar oracle + second-model oracles (field-generic) |
| Distinct-module `Ext(M, N)`, `M ≠ N` | flagged post-v1; self-Ext is the confirmed QPA idiom |
| **The Koszul verdict itself** (Plan 27 — QPA 1.37 has no `IsKoszul`/`KoszulDual`, confirmed by an exhaustive `NamesGVars()` sweep) | the G-quadratic certifier (Priddy PBW: confluent length-2-tip reduction system ⇒ Koszul) plus the generated-in-degree-1 falsifier and the Fröberg matrix identity `P(t)·C_A(−t)=I`; QPA validates every INPUT to the verdict — graded Ext dims and minimal-generator degrees (`ExtAlgebraGenerators`), quadraticity (`IsQuadraticIdeal`), and the quadratic perp (`QuadraticPerpOfPathAlgebraIdeal`) |
| Yoneda relations-by-degree (QPA exposes generator counts, not a presentation) | theory battery: `E(k[x]/x²)=k[y]`, `E(k[x]/xⁿ)=k[y,z]/(y²)` (char-independent, pinned over GF(2)/GF(3)/GF(32003)/char 0), hereditary `E=kQ/J²` with `as_algebra()` round-trip, rad²=0 `E=kQ` path counts, quantum CI `dim Eⁿ=n+1` (= the CS chain count), commutative square `E≅A` self-hosting; the monomial Anick chain-count gate; byte-reproducible lift products (Plan-17-style canonicalization) |

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
| `hochschild/` (bar, cyclic; the Plan-34 auto→CS depth-fallback battery; the Plan-35 product surface — `products.py`: cup/cap/bracket tables + the induced Connes `B`, `basis_reps.py`: the explicit-representatives capture, `cyclic_reps.py`: the wave-3b cyclic-homology explicit representatives, and `hh_reps.py`: the wave-3d plain-HH explicit representatives) | 103 | fast | **the base bar oracle**; mixed-complex identities; dispatch-amendment pins; **the Gerstenhaber identity batteries** (graded-commutative + associative cup, antisymmetric bracket, cup-Leibniz, cap module law, `B²=0`, SBI rank) + the `k[x]/(x²)`/QuantumCI-BGMS product literature pins + the bar↔CS in-window cross-engine gate + **the explicit-reps self-certification** (every shipped product class satisfies `δ·v = 0` / `b·v = 0` from its shipped or note-rebuilt differential; hand-checked `k[x]/(x²)` labels; elision+rebuild path) + **the cyclic-homology explicit representatives** (every shipped HC class is a cycle of the (b,B) total complex — `D·v = 0` on both engines; GF(p)≡generic on prime 32003; hand-checked dual-numbers `HC_0 = A/[A,A]` + `Tot_2 = C_2 ⊕ C_0` column structure) + **the plain-HH explicit representatives** (the `hh_cohomology`/`hh_homology` dims blocks now carry per-degree reps over BOTH routes — the GF(p) bar and Chouhy-Solotar; every shipped vector annihilates its shipped differential; the `k[x]/(x²)` HH¹ `[x ↦ x]` = derivation `D(x)=x` hand-check ties the classical dictionary read-off to the captured representative) |
| `engine/` (fast GF(p); minimal, Bardzell, periodic; TT-calculus; cyclic; Coxeter/Nakayama; Plan-29 literature/identity batteries) | 579 | deep | bar oracle; cross-engine; multi-prime; numba/pure parity; frozen QPA-literature values |
| `resolutions_cs/` (CS; comparison; diagonal; cup; cap; Plan-29 literature batteries; the Plan-35 Domain-generic CS product tables `products.py` — cup/cap on the CS basis over any exact Domain; the Plan-35 wave-3d `cs_hh_basis` 0-codomain fix) | 228 | deep | CS ≡ bar, CS ≡ Bardzell; bank byte-level; literature pins; `d∘d=0` / order; Leibniz + cap identities (unit/module/transport anchors); canonicalization; the CS product unit-law + Domain-genericity self-cert; **the `cs_hh_basis` rep-count ≡ `cs_(co)homology_dims`** even when the top differential lands in a 0-dimensional space (the multi-vertex `kZ_3/J^2` witness — HH¹ = 1 with C² = 0 — that the old 0-row `nullspace` silently undercounted) |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity; Plan-29 Tor; Plan-30 Krull–Schmidt decomposition; the retained injective-coresolution differentials certified exact; the Plan-35 wave-3a explicit Ext/Tor representatives — `complex_reps.py`; the Plan-35 wave-3c Yoneda exact sequences — `yoneda.py`; the Plan-37 C1 categorical glue — `morphism.py` first-class `ModuleHom` + kernel/image/cokernel, `ses.py` short exact sequences + split test + pushout/pullback, `endomorphism.py` `End(M)` as an Algebra, `direct_sum`/`is_direct_summand`, and covers/envelopes + radical/socle series + composition factors on `Module`; the Plan-39 C8 complex layer — `complexes.py` validated bounded `ChainComplex`/`ChainMap`, shift/truncation/homology, mapping cones + triangles + the derived-iso test, the Hom total complex `hyper_hom_dims`, and certified projective models `projective_model` + general `hyper_ext_dims`) | 317 | deep | AR/duality literature pins (ASS2006); **the Plan-37 categorical-glue self-certification** (`ModuleHom` validates the intertwining relations at construction; kernel/image/cokernel certified by rank-nullity + the epi–mono factorization `f = epi∘mono` + `f∘iota = 0` = `proj∘f`; SES exactness = the rank identity `im f = ker g`; split ⇔ a section solves; pushout/pullback squares certified by their universal-square identities; `End(M)` self-certified by `from_structure_constants(check=True)` with the regular-module `End(A_A) ≅ A` Loewy oracle; biproduct identities `proj_i∘incl_i = id`, `Σ incl_i∘proj_i = id`); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks; **the explicit Ext/Tor self-certification** (every shipped class satisfies `δ·v = 0` (Ext cocycle) / `d·v = 0` (Tor cycle) from its shipped differential; hand-checked kA₂ `Ext¹(S₁,S₂)` + loop `Tor₀ = M ⊗ N` cokernel labels; rep-count ≡ engine dims) + **the Yoneda exact-sequence self-certification** (every `Ext^n(M,N)` class is CONSTRUCTED as an `n`-fold exact sequence `0 → N → Q → … → M → 0` — the pushout middle module + connecting maps — and its exactness is self-certified at every joint: each map an `A`-module map, ends injective/surjective, `im = ker` by rank; the kA₂ Baer pin `0 → S₂ → P₁ → S₁ → 0` verified by the library's OWN `is_isomorphic`/`identify_standard`; a non-cocycle is refused loudly; the multi-vertex `_tor_boundary` collapse pinned on a rad²=0 Nakayama) + **the Plan-39 complex-layer oracles** — self-cert: `d∘d=0` refused at construction, the mapping-cone `d²=0` re-asserted under full validation, quasi-iso ⇔ cone acyclicity, shift/truncate degree+sign identities, and the Hom total-complex `δ∘δ=0` block-indexing certificate (the Weibel `−(−1)^n` sign is a documented convention — verified sign-independent for the dims: both `±` give isomorphic cochain complexes); cross-engine: stalk `hyper_hom_dims` ≡ module `Ext` degreewise (kA₃/(ab), all vertices), the resolution-augmentation quasi-iso, the Euler-characteristic triangle identity `χ(cone) = −χ(X)+χ(Y)`, the certified projective model (`projective_model` asserts `is_perfect` + `is_quasi_iso` before return — NEVER returns uncertified), and the sharpened two-term shift identity `hyper_ext([P₁→S₁])[n] = Ext^{n−1}(rad P₁, N)` (the derived-category source-shift arithmetic, replacing the plan's placeholder); live QPA Ch.10 (`qpa/`) |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity; Plan-29 Tor; Plan-30 Krull–Schmidt decomposition; the retained injective-coresolution differentials certified exact; the Plan-35 wave-3a explicit Ext/Tor representatives — `complex_reps.py`; the Plan-35 wave-3c Yoneda exact sequences — `yoneda.py`; the Plan-37 C1 categorical glue — `morphism.py` first-class `ModuleHom` + kernel/image/cokernel, `ses.py` short exact sequences + split test + pushout/pullback, `endomorphism.py` `End(M)` as an Algebra, `direct_sum`/`is_direct_summand`, and covers/envelopes + radical/socle series + composition factors on `Module`; the Plan-40 C6 homological-dimensions family — `homdims.py`: public `syzygy`/`cosyzygy` (byte-stable extraction from `minimal_resolution`), the Igusa–Todorov φ/ψ on the finite K₀, dominant + Gorenstein dimensions, Ω/τ-periodicity certificates, and finitistic-dimension bounds) | 322 | deep | AR/duality literature pins (ASS2006); **the Plan-37 categorical-glue self-certification** (`ModuleHom` validates the intertwining relations at construction; kernel/image/cokernel certified by rank-nullity + the epi–mono factorization `f = epi∘mono` + `f∘iota = 0` = `proj∘f`; SES exactness = the rank identity `im f = ker g`; split ⇔ a section solves; pushout/pullback squares certified by their universal-square identities; `End(M)` self-certified by `from_structure_constants(check=True)` with the regular-module `End(A_A) ≅ A` Loewy oracle; biproduct identities `proj_i∘incl_i = id`, `Σ incl_i∘proj_i = id`); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks; **the explicit Ext/Tor self-certification** (every shipped class satisfies `δ·v = 0` (Ext cocycle) / `d·v = 0` (Tor cycle) from its shipped differential; hand-checked kA₂ `Ext¹(S₁,S₂)` + loop `Tor₀ = M ⊗ N` cokernel labels; rep-count ≡ engine dims) + **the Yoneda exact-sequence self-certification** (every `Ext^n(M,N)` class is CONSTRUCTED as an `n`-fold exact sequence `0 → N → Q → … → M → 0` — the pushout middle module + connecting maps — and its exactness is self-certified at every joint: each map an `A`-module map, ends injective/surjective, `im = ker` by rank; the kA₂ Baer pin `0 → S₂ → P₁ → S₁ → 0` verified by the library's OWN `is_isomorphic`/`identify_standard`; a non-cocycle is refused loudly; the multi-vertex `_tor_boundary` collapse pinned on a rad²=0 Nakayama) + **the Plan-40 homological-dimensions oracles**: `oracle_selfcert` — the φ=ψ=pd identity for finite projective dimension, the Ω/τ-periodicity `is_isomorphic` certificates, and the decompose char-caveat propagation; `oracle_literature` — the Barrios–Mata truncated self-injective φ=ψ=0 closed form + projective additivity, the hereditary/self-injective dominant & Gorenstein values, and the cyclic-Nakayama period-from-Kupisch pins; `qpa` (`tests/qpa/test_homdims_qpa.py`) — live `GlobalDimensionOfAlgebra` / `DominantDimensionOfAlgebra` / `GorensteinDimensionOfAlgebra` agreement over kA₂ / kA₃(ab) / `line_abc_cde` / k[x]/(x³) (int or GAP `infinity`↔our infinite/unresolved marker) |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity; Plan-29 Tor; Plan-30 Krull–Schmidt decomposition; the retained injective-coresolution differentials certified exact; the Plan-35 wave-3a explicit Ext/Tor representatives — `complex_reps.py`; the Plan-35 wave-3c Yoneda exact sequences — `yoneda.py`; the Plan-37 C1 categorical glue — `morphism.py` first-class `ModuleHom` + kernel/image/cokernel, `ses.py` short exact sequences + split test + pushout/pullback, `endomorphism.py` `End(M)` as an Algebra, `direct_sum`/`is_direct_summand`, and covers/envelopes + radical/socle series + composition factors on `Module`; the Plan-41 C3 **Auslander–Reiten completion** — `ar.py`: the general chain-map lift, the Nakayama functor ν/ν⁻, stable Hom mod projectives, the End(M)-action on Ext¹, almost-split sequences `0 → τM → E → M → 0`, irreducible-map multiplicities `dim rad(M,N)/rad²`, and honest-semi-decision AR-quiver knitting) | 327 | deep | AR/duality literature pins (ASS2006, ARS1995); **the Plan-41 AR self-certification** (the chain-map lift asserts every square `d_n·φ_n = φ_{n-1}·d_n` and is byte-reproducible; ν tied to the trusted τ by `ker(νP₁→νP₀) ≅ τM` + `ν(P_v) ≅ I_v` / `ν⁻(I_v) ≅ P_v`; every almost-split sequence's class is identified by the **ARS socle theorem** — `soc_{End M} Ext¹(M,τM)` is a *simple* `End(M)`-module — pinned down computationally by the char-scoped trace-form `rad End(M)` (char 0 or char > dim M) plus the **socle-simplicity dimension certificate** `dim_k soc = dim_k End(M) − dim_k rad End(M)` (`= dim_k` the residue division algebra; refuses loudly if it fails); the exact (Yoneda + P37 SES) + non-split (P37 `is_split` False) + indecomposable-ends (Plan 30) checks are **necessary sanity checks, not sufficient** — they do NOT arbitrate the pick (over `k[x]/(x⁴)` with `M = k[x]/(x²)`, `dim Ext¹ = 2`, a non-socle class has an exact/non-split/indecomposable-ends extension whose middle is the *projective* `k[x]/(x⁴)`, not the true mesh middle `{1,3}` — the devil's-advocate non-brick regression, live-QPA-crosschecked)) + **the cross-engine AR formula** `dim Ext¹(M,N) = dim underline-Hom(τ⁻N, M)` + **literature Dynkin/Nakayama pins** (kA₂/kA₃/kA₄ = 3/6/10 indecomposables, D₄ = 12, the kA₃ mesh middle terms, Nakayama serial count) + **live QPA** `AlmostSplitSequence` middle-term (dim vector over QQ; summand multiset over GF(p)) and `PredecessorsOfModule`; **the Plan-37 categorical-glue self-certification** (`ModuleHom` validates the intertwining relations at construction; kernel/image/cokernel certified by rank-nullity + the epi–mono factorization `f = epi∘mono` + `f∘iota = 0` = `proj∘f`; SES exactness = the rank identity `im f = ker g`; split ⇔ a section solves; pushout/pullback squares certified by their universal-square identities; `End(M)` self-certified by `from_structure_constants(check=True)` with the regular-module `End(A_A) ≅ A` Loewy oracle; biproduct identities `proj_i∘incl_i = id`, `Σ incl_i∘proj_i = id`); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks; **the explicit Ext/Tor self-certification** (every shipped class satisfies `δ·v = 0` (Ext cocycle) / `d·v = 0` (Tor cycle) from its shipped differential; hand-checked kA₂ `Ext¹(S₁,S₂)` + loop `Tor₀ = M ⊗ N` cokernel labels; rep-count ≡ engine dims) + **the Yoneda exact-sequence self-certification** (every `Ext^n(M,N)` class is CONSTRUCTED as an `n`-fold exact sequence `0 → N → Q → … → M → 0` — the pushout middle module + connecting maps — and its exactness is self-certified at every joint: each map an `A`-module map, ends injective/surjective, `im = ker` by rank; the kA₂ Baer pin `0 → S₂ → P₁ → S₁ → 0` verified by the library's OWN `is_isomorphic`/`identify_standard`; a non-cocycle is refused loudly; the multi-vertex `_tor_boundary` collapse pinned on a rad²=0 Nakayama) |
| `derived/` (Plan-43 C8 **derived-category surface** — `homs.py`: reified hyper-Hom classes `hyper_hom_basis` (a basis of `H^n(Hom^•(X,Y))` as genuine chain maps `X → Y[n]`) + `ChainMap.then`; `tau.py`: the derived AR translate `τ_{D^b} = ν[−1]` / `τ⁻_{D^b}` on perfect complexes with the Happel finite-gl.dim gate; `_corner.py`: the shared corner-transpose `Hom_A(−,A)` (factored out of `duality._presentation_transpose`, imported by both); `tilting.py`: the tilting-complex verifier (rigidity on the exact window + K₀-unimodular generation) + `End(T)` as a structure-constant algebra + `corner_cartan_of_complex` + `two_term_silting_from_presentation`; `fingerprint.py` + `block.py`: the necessary-condition derived fingerprint) | 22 | deep | **self-cert** — `hyper_hom_basis` reifies genuine chain maps (`ChainMap(check=True)` on every coset representative) with count ≡ `hyper_hom_dims`; `τ_{D^b}` output `d∘d=0` (ChainComplex check) and the `τ⁻∘τ` round-trip is a degreewise quasi-iso; `End(T)` self-certified by `from_structure_constants(check=True)`; the tilting rigidity window reported honestly; **cross-engine** — `hyper_hom_basis` count ≡ module `Ext^n` on a projective-resolution source (kA₃/(ab), all vertices), and `τ_{D^b}` homology is concentrated in degree 0 and `≅` the trusted module `τ` over `kA_n`; **literature** — the K₀-action identity `χ(τ_{D^b}X) = c·χ(X)` with `c = −C·C⁻ᵀ` (the conjugate of P38's Coxeter matrix — same char poly, the dim-vector action), the kA₂ APR-tilt `P₁ ⊕ S₁` (tilting, `End(T)` = the reoriented A₂ = A^op, corner-Cartan `[[1,0],[1,1]]` — the theorem-anchored `End(A_A)=A` pin fixes the orientation), the D₄ vs A₄ Coxeter distinction and the 8-vertex cospectral-trees NON-distinction (the honest-scope demonstration); **live QPA** (`tests/qpa/test_derived_qpa.py`) — `τ_{D^b}` homology(0) vs `DTr(M)` (the documented module-level route: QPA's `TauOfComplex` on a `ProjectiveResolution` does not script through libgap — the P39 Ch.10 hazard, confirmed live) |
| `invariants/` (Cartan, Coxeter, spectral, Betti, cyclic, Frobenius incl. the Plan-29 trace-form symmetry certifier, scalar, sweep; Plan-29 Coxeter/identity literature batteries) | 115 | fast | second models (λ-complex, relative-Tor Betti); self-certifying `λ`/`ν`; GF(p) engine parity |
| `invariants/geometry.py` (Plan-49 C8 — orbit dimension `dim O_M = Σ d_v² − dim End(M)`, Voigt rigidity `is_rigid`/`rigidity_codim`, the Kac `canonical_decomposition` over hereditary Dynkin, and the shared `orbit_geometry_block`) — `tests/invariants/test_geometry_orbit.py`, `test_geometry_canonical.py` | 21 | fast | **`oracle_selfcert`**: the orbit-dim identity `dim O_M = Σ d_v² − dim End(M)`, GF(p)↔QQ field parity, the canonical-decomposition sum-of-roots + per-instance rigidity certificate `Ext¹(G,G)=0`, the loud Euclidean-deferred / non-hereditary refusals. **`oracle_crossengine`**: the Voigt codim identity `dim Rep(Q,d) − dim O_M ≡ dim Ext¹(M,M)` on hereditary + the P38 `tits_form` tie; the canonical decomposition ≡ the Krull–Schmidt summands of the degeneration poset's maximum. **`oracle_literature`**: every Dynkin indecomposable is rigid (codim 0); `(2,1) = P₁ ⊕ S₁` over kA₂ (hand-derived Kac pin) |
| `modules/degeneration.py` (Plan-49 C8 — the Zwara–Bongartz degeneration = hom order poset for representation-finite algebras, `DegenerationPoset`) — `tests/modules/test_degeneration.py` | 6 | deep | **`oracle_literature`**: kA₂ (1,1) = the 2-chain `S₁⊕S₂ <_deg P₁`; kA₃ (1,1,1) = the diamond (orbit dims `[0,1,1,2]`, 4 covers, two incomparable middles) — both hand-derived. **`oracle_selfcert`**: the hom-order is a partial order (reflexive/antisymmetric), orbit dim strictly increases up every cover, the per-class orbit dim ≡ `geometry.orbit_dimension`, and the honest semi-decision cap (rep-infinite / self-injective ⇒ `is_complete=False` with a loud `status`, never a partial poset) |
| `families/` (catalog, zoo; Plan-29 trivial-extension/incidence batteries; Plan-31 certified trivial-extension presentation, `test_trivial_extension_presented.py`) | 166 | deep | closed-form family pins; zoo diversity gates; citations; Plan-31 special-case + Cartan + iso-invariance + CS≡bar pins |
| `strings/` (Plan-46 C5 gentle/string subsystem: reduced walks + σ/ε signs + string census + band detection; string/band module materialisation; string-τ by hooks/cohooks; the Avella-Alaminos–Geiss derived invariant; the `strings` block) — `tests/modules/test_strings_*.py`, `tests/invariants/test_ag_invariant.py` | 33 | deep + fast | Butler-Ringel `n(n+1)/2` interval count + Kronecker band existence (`oracle_literature`); string-τ ≡ engine τ + census count ≡ `knit_ar_quiver` vertex count (`oracle_crossengine`); `check_module` on every materialised string/band + permitted/forbidden thread partition of `Q_1` (`oracle_selfcert`); AAG-2008 pins reproduced verbatim (Nakaoka `arXiv:1811.00775` Example 2.15 = `{(3,2),(2,2),(0,3)}`) |
| `families/brauer.py` (Plan-46 Brauer graph algebra constructor from a ribbon graph + multiplicities) — `tests/families/test_brauer.py` | 10 | deep | `dim = Σ_v m_v·val(v)²` per-instance certificate + `is_symmetric` (`oracle_literature`); Brauer-star ≡ symmetric Nakayama `NakayamaAlgebra(n, mn+1, cyclic=True)` byte-equal Cartan (`oracle_crossengine`) |
| `surfaces/` (Plan-48 marked surfaces → ideal triangulations → gentle Jacobian algebras — `marked.py` `MarkedSurface`, `triangulation.py` `Triangulation` + fan/annulus/hexagon/once-punctured-torus, `qp.py` `quiver_of`/`potential_of`/`jacobian_of`, `flip.py` `flip`/`certify_flip_mutation`, `block.py` `surface_block`) — `tests/families/test_surfaces_*.py`, `tests/invariants/test_surfaces_arccount.py`, `tests/qpa/test_surfaces_qpa.py` | 70 | deep + fast + qpa | **`oracle_literature`**: the derived FST arc count `n = 6g−6+3(b+p)+Σkᵢ` on disc(n+3)→Aₙ / annulus(n,m)→n+m / once-punctured-torus→3, the FST admissibility exclusion list (monogon/digon/triangle, spheres with ≤3 punctures, once-punctured monogon), the hexagon-internal-triangle Jacobian `dim 6` (P44's pin), the disc-fan `kAₙ`, and the small annulus `C(2,1)` acyclic affine-`Ã₂`; **`oracle_crossengine`** (the P44+P46+P48 cross-subsystem tie): the disc-fan-`Aₙ` **orientation arbiter** (`quiver_of(fan((n+3))) = 1→2→…→n` exactly), `is_gentle(jacobian_of(T))` True across the disc/annulus/hexagon zoo (ABCP/LFS), flip ≡ Fomin–Zelevinsky matrix mutation on **every** interior arc, and `surface_block` AG invariant ≡ `strings.ag.ag_invariant`; **`oracle_selfcert`**: the two side-counting identities `3t=2n+c` / `p−n+t=χ`, the arc-adjacency + arc-count self-cert on every constructor (interior arcs in 2 triangles, boundary segments in 1), the self-folded refusal, and flip involution + `μₖ∘μₖ=id`; **`qpa`** (`tests/qpa/test_surfaces_qpa.py`): `IsGentleAlgebra`/`IsSpecialBiserialAlgebra` parity on the surface Jacobians + the standing `IsBoundGlobal` guard that FAILS if QPA ever ships a surface/triangulation constructor |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity; Plan-29 Tor; Plan-30 Krull–Schmidt decomposition; the retained injective-coresolution differentials certified exact; the Plan-35 wave-3a explicit Ext/Tor representatives — `complex_reps.py`; the Plan-35 wave-3c Yoneda exact sequences — `yoneda.py`; the Plan-37 C1 categorical glue — `morphism.py` first-class `ModuleHom` + kernel/image/cokernel, `ses.py` short exact sequences + split test + pushout/pullback, `endomorphism.py` `End(M)` as an Algebra, `direct_sum`/`is_direct_summand`, and covers/envelopes + radical/socle series + composition factors on `Module`; the Plan-40 C6 homological-dimensions family — `homdims.py`: public `syzygy`/`cosyzygy` (byte-stable extraction from `minimal_resolution`), the Igusa–Todorov φ/ψ on the finite K₀, dominant + Gorenstein dimensions, Ω/τ-periodicity certificates, and finitistic-dimension bounds; the Plan-44 C7 slice — `approximations.py` minimal left/right add(M)-approximations and `tilting.py` `is_tilting_module`/`is_cotilting_module` + self-certified `bongartz_completion`) | 341 | deep | AR/duality literature pins (ASS2006); **the Plan-37 categorical-glue self-certification** (`ModuleHom` validates the intertwining relations at construction; kernel/image/cokernel certified by rank-nullity + the epi–mono factorization `f = epi∘mono` + `f∘iota = 0` = `proj∘f`; SES exactness = the rank identity `im f = ker g`; split ⇔ a section solves; pushout/pullback squares certified by their universal-square identities; `End(M)` self-certified by `from_structure_constants(check=True)` with the regular-module `End(A_A) ≅ A` Loewy oracle; biproduct identities `proj_i∘incl_i = id`, `Σ incl_i∘proj_i = id`); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks; **the explicit Ext/Tor self-certification** (every shipped class satisfies `δ·v = 0` (Ext cocycle) / `d·v = 0` (Tor cycle) from its shipped differential; hand-checked kA₂ `Ext¹(S₁,S₂)` + loop `Tor₀ = M ⊗ N` cokernel labels; rep-count ≡ engine dims) + **the Yoneda exact-sequence self-certification** (every `Ext^n(M,N)` class is CONSTRUCTED as an `n`-fold exact sequence `0 → N → Q → … → M → 0` — the pushout middle module + connecting maps — and its exactness is self-certified at every joint: each map an `A`-module map, ends injective/surjective, `im = ker` by rank; the kA₂ Baer pin `0 → S₂ → P₁ → S₁ → 0` verified by the library's OWN `is_isomorphic`/`identify_standard`; a non-cocycle is refused loudly; the multi-vertex `_tor_boundary` collapse pinned on a rad²=0 Nakayama) + **the Plan-40 homological-dimensions oracles**: `oracle_selfcert` — the φ=ψ=pd identity for finite projective dimension, the Ω/τ-periodicity `is_isomorphic` certificates, and the decompose char-caveat propagation; `oracle_literature` — the Barrios–Mata truncated self-injective φ=ψ=0 closed form + projective additivity, the hereditary/self-injective dominant & Gorenstein values, and the cyclic-Nakayama period-from-Kupisch pins; `qpa` (`tests/qpa/test_homdims_qpa.py`) — live `GlobalDimensionOfAlgebra` / `DominantDimensionOfAlgebra` / `GorensteinDimensionOfAlgebra` agreement over kA₂ / kA₃(ab) / `line_abc_cde` / k[x]/(x³) (int or GAP `infinity`↔our infinite/unresolved marker) |
| `families/` (catalog, zoo; Plan-29 trivial-extension/incidence batteries; Plan-31 certified trivial-extension presentation, `test_trivial_extension_presented.py`; the Plan-44 C7 constructions — `one_point.py` `OnePointExtension`, `repetitive.py` `repetitive_slice`, `jacobian.py` `Potential`/`JacobianAlgebra`/`cyclic_derivative`, and the Task-C Gabriel-recovery battery `test_gabriel_recovery.py` over `core/basic.py`) | 187 | deep | closed-form family pins; zoo diversity gates; citations; Plan-31 special-case + Cartan + iso-invariance + CS≡bar pins; **the Plan-44 construction oracles** — `oracle_literature`: the one-point Cartan block `[[1, dim-vector M],[0,C_A]]` + `pd(S_ω)=pd_A(M)+1`, `repetitive_slice` `copies=1==A` + the `dim==(2·copies−1)·dim A` slice certificate, the hand-derived Jacobian triangle `dim=6`, and the `M₂(k)→k` / `kA₂` Gabriel round-trip; `oracle_crossengine`: `Jac(3-cycle, abc) ≅ cyclic Nakayama kZ₃/J²` (dim + Cartan) and `presented_form(End(⊕P_v)) ` recovers `kA₃` (tied to P37 `regular_corner_dims`); `oracle_selfcert`: complete-orthogonal primitive idempotents + the per-instance dimension/multiplicativity recovery certificate, the cyclic-derivative identities, and the loud char/split + `NotFiniteDimensionalError` refusals |
| `modules/quasihereditary.py` + `modules/recollement.py` (Plan-47 C-series: standard/costandard modules Δ(i)/∇(i), the quasi-heredity test `QHReport`, good-filtration multiplicities + BGG reciprocity, the characteristic tilting module + Ringel dual; `Recollement(A,S)` — the corner algebra `eAe`, the quotient `A/AeA`, and the six functors) — `tests/modules/test_quasihereditary_*.py`, `tests/modules/test_recollement.py` | 42 | deep | **`oracle_literature`** — Dlab–Ringel: `kA_n` natural-order Δ(i)=S_i / ∇(i)=[1..i], opposite-order Δ(i)=P(i); `kA_n` quasi-hereditary for BOTH orders; `k[x]/(x²)` NOT quasi-hereditary (loud note); Ringel: T = D(A) (natural) / T = A (opposite); the double-Ringel-dual Cartan Smith-form identity; **the discriminating NON-hereditary oracle** — the commutative square `1→2→4, 1→3→4, ab=cd` (gl.dim 2, every Δ simple) has T = D(A) (dim 9, dim-vector {1:4,2:2,3:2,4:1}), which the classical pd≤1 tilting certificate rejects and the single-pass `T(j)`-extension missed; the **`eAe`-vs-subquiver trap** (`kA₃`, S={1,3} → `eAe = kA₂` dim 3, NOT the subquiver `k×k` dim 2) + the worked `kA₃`, S={2} recollement. **`oracle_crossengine`** — BGG reciprocity `(P(i):Δ(j)) = [∇(j):S(i)]` across both orders. **`oracle_selfcert`** — top Δ(i)=S_i + [Δ(i):S(i)]=1 + socle ∇(i)=S_i; the greedy Δ-peel certificate (loud `certified=False` when no filtration); the characteristic-tilting arbiter — the genuine Ringel iteration lands in `F(Δ)∩F(∇)` (`Ext¹(Δ(j),T)=0` AND `Ext¹(T,∇(j))=0` for all j) with `is_tilting_module(T, n=gl.dim A)` (the exact global dimension, not the classical 1; gl.dim non-exact refuses loudly); the six-functor **adjunction dim identities** `dim Hom_A(j_!X,M)=dim Hom_{eAe}(X,j^*M)` etc., the counit isos `j^*j_!≅id`/`j^*j_*≅id`, and — through the ACTUAL functor outputs (each natural map a `ModuleHom` re-certified as an A-map) — BOTH BBD exact sequences: the counit `j_!j^*M → M → i_*i^*M → 0` (`im(counit)=ker(unit)`, unit epi) and the unit `0 → i_*i^!M → M → j_*j^*M` (mono, `im=ker`), over `kA₃`/the commutative square/`kA₅`/`kA₄`; the degenerate-`S` refusals (`S`=all → e=1, `S`=∅ → e=0); a `GF(2)` cell in each battery (Δ/∇/qh/recollement are char-clean). **QPA white space** — QPA has NO quasi-hereditary / recollement surface, so there is no `qpa` row here (stated in Honest scope) |
| `tautilting/` (Plan-45 C4 τ-tilting engine — `rigid.py` g-vectors + τ-rigidity, `pairs.py` certified support τ-tilting pairs, `mutation.py` the 2-term silting exchange + `exchange_graph` BFS (`_twoterm.py` = the K^b(proj A) cone/cocone + minimal-complex reduction engine), `torsion.py` the torsion lattice + Hasse orientation + bricks/semibricks, `stability.py` King θ-stability + the wall-and-chamber fan, `green.py` maximal green sequences, `silting.py` the 2-term silting bridge, `block.py` the algebra-level payload) — `tests/modules/test_tau_tilting_*.py` | 34 | deep | **`oracle_literature`** — `#sτ-tilt(kA_n) = Catalan(n+1)` (2/5/14), exchange-graph n-regularity, the AIR **four-way count identity** `#sτ-tilt = #f.f. torsion = #2-term silting = #semibricks` on kA₂/kA₃, hereditary `τ-rigid ⇔ rigid`, kA₂ = 2 maximal green sequences, and the **non-thin gate** kZ₂/rad² (symmetric Nakayama, dim 4): FOUR pairwise-non-isomorphic bricks — S₁, S₂ and the two projective-injectives P₁, P₂ *both* of dim-vector (1,1) — six semibricks, and the four-way identity 6 = 6 = 6 = 6 (the identity a dim-vector-keyed brick count silently broke: bricks → 3, semibricks → 5); **`oracle_selfcert`** — `g^{P_v}=e_v` + additivity, the four-axiom pair certification, mutation is an **involution** swapping exactly one g-column, every chamber g-matrix unimodular (det ±1), the n=2 fan **tiles R²** by an exact angular sweep (no atan2), the **n=3 L1/octahedron unfolding sanity** (nondegenerate net triangles, one per chamber — a rendering check, NOT a 3D tiling certificate), the **iso-class wall labelling** (kZ₂/rad²'s two same-dim-vector (1,1) walls carry distinct P₁ vs P₂ labels — four brick iso-classes across the walls, never a first-dim-vector collapse), King θ-stability on the worked kA₂ example, unique Hasse source/sink, brick `end_dim=1`, semibrick Hom-orthogonality, and the **honest semi-decision contract** (the 2-Kronecker is τ-tilting-infinite → loud `status="budget"`); **`oracle_crossengine`** — pair ↔ `Gen(M)` torsion-class injectivity, the fan's wall brick-normals ⊥ the shared g-facet (King), and the **no-proactive-char-guard** pin (GF(2) kZ₂/rad² agrees with the certified QQ counts — char ≤ dim computes where every module is a brick/splits). The webapp/GUI cross-runner τ-tilting tests are UNMARKED (extras-gated dirs, Plan-32 ruling). **QPA cannot compare — there is no `qpa` battery for τ-tilting** (see the honest-scope section). |
| `batch/` (labdb port, open-zone scans) | 11 | deep | labdb port equality; scan-surface checks |
| `citations/` (registry, bibliography) | 12 | fast | packaged-bib resolution; result references |
| `trace/` (worked-steps incl. the Plan-30 module events, the kA₂ replay golden, the 2026-07-29 report-completeness battery, the Plan-35 UNIT-2 HH explicit-reps rendering, the Plan-35 wave-3a Ext/Tor explicit-reps rendering, the Plan-35 wave-3b cyclic-homology explicit-reps rendering — the total-complex `Tot_n = C_n ⊕ C_{n-2} ⊕ …` column heading, per-degree classes + verification; and the Plan-35 wave-3c Yoneda-sequence + classical-dictionary rendering — `interpretations.py`; and the Plan-35 wave-3d plain-HH explicit-reps + element-wise dictionary rendering — `hh_element_interpretation`/`hh_reps_sections`) | 250 | fast | golden-file equality (dims derived from ranks); **the per-degree explicit-reps layout** (each product/Connes class rendered as term-sum + coordinate vector under a stable anchor, with the annihilating differential + a one-line verification sentence; the bar AND Chouhy-Solotar HH worked-steps carry each (co)chain term's ordered basis, length-guarded against the recorded term dim; module resolution `term_basis` lengths match the differential row/col dims, injective order pinned against the transposed proj-resolution-of-DM; the degree anchors are linked from every product table) + **the module Ext/Tor per-degree sections** (ordered Hom/tensor basis → classes → differential + verification, `cr-`/`ws-` anchors, the `ExtReps` worked-steps event, Tor₀ = M ⊗ N cokernel note) + **the Yoneda-sequence + dictionary rendering** (each Ext class' constructed exact sequence — sequence line, middle module, exactness verified — under `cr-ext-yoneda-deg-n`; the shared classical-dictionary framing on the ext/tor/HH/cyclic blocks; the HH¹ derivation read-off; matrix-grid double zebra striping is structure-safe) + **the plain-HH element-wise dictionary + per-degree reps** (HH⁰'s central elements, HH¹'s `D(arrow)=value` derivations + the inner-derivation subspace dimension `rank δ⁰`, HH²'s deformation 2-cocycle, HH₀'s commutator residues — read straight off the captured term-sums; the per-degree explicit-reps sections under `cr-hh_cohomology`/`cr-hh_homology` anchors; both gui.js copies mirror it) + the missing-fields tolerance + the two-runner `term_basis`/reps/interpretation equality |
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
| `fast` | 1556 | every CI cell: `{ubuntu, macos, windows} × py{3.10, 3.11, 3.12, 3.13}` |
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
| Literature / theory pins | `-m oracle_literature` | 882 | the engine reproduces a value/identity that exists outside the library |
| Cross-engine agreement | `-m oracle_crossengine` | 528 | two independent implementations compute the same thing and match live |
| Self-certifying certificates | `-m oracle_selfcert` | 1080 | an internal axiom (d∘d=0, canonicality, an arbitration identity) holds by construction |
| Live QPA / GAP | `-m qpa` | 189 | an independent external system (QPA) recomputes and agrees |
| Live Macaulay2 | `-m m2` | 11 | an independent external system (Macaulay2) recomputes and agrees |
| Any oracle class (union) | `-m "oracle_literature or oracle_crossengine or oracle_selfcert or qpa or m2"` | 2118 | the test is pinned by at least one oracle (the remaining tests are contract/infrastructure) |
Counts as of the P43 merge (the derived-category surface); sibling plans in the v0.2.0

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

- `complexity` is a lower-bound estimate that can under-report; it is **exact only
  on local / single-vertex inputs**. The Betti-complex identity it rests on
  (`H_n =` the minimal resolution's ranks) is gated exactly over every field.
- `is_symmetric` off GF(p) decides the definitional "`ν` is inner" by a
  Schwartz–Zippel sweep — **loud when inconclusive**, never a silent wrong answer.
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
  dims. Two scope facts are binding: (i) the **Gerstenhaber bracket is GF(p)-only and
  window-bounded** — it is served on the bar/tt route (the result object records the
  window), there is no CS route for it, and the degree-0 insertion action is out of
  scope; the CS-native cup and cap tables compute over **any exact Domain** (the
  bracket refuses on the CS route, `tests/resolutions_cs/test_products_cs.py`).
  (ii) The **structure constants are basis-dependent** — they are read on the
  recorded HH basis (bar/GF(p) or the CS class basis), and each product object
  records which basis (`HHProducts.basis`); the cross-engine gate therefore compares
  only basis-independent data (dims and flattened rank), never the raw constants.
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
  `cup:0..2` probe on kZ₂₀ that timed out; and bracket/Connes `B` have no CS route at
  all. Every trim and omission is recorded per-example in
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
  decider. `is_tilting_complex` decides rigidity (on the exact, honestly-reported
  window outside which hyper-Hom is provably zero) and K₀-unimodular generation for
  a *given* candidate — it does not search for tilting complexes; `End(T)` is the
  Rickard derived-equivalent algebra of that candidate. `derived_fingerprint` /
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
