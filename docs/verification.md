# How quiverlab is verified

Every shipped feature of quiverlab is unit tested. This page says **how** — with
the highest rigour we can bring to it — and it is honest about the edges: where a
check is a cross-engine agreement, where it is a published number, where a live
external oracle can reach, and where it cannot.

The suite is **2099 tests** (collected with the `[dev,fast,docs,web,qpa,hpc]` extras,
2026-07-26 post-merge of Plans 21–31; Plan 32 added the oracle-class markers + the
audit gate). It
is not a pile of smoke tests: the mathematics is pinned by **two classes of
oracle**, and most numbers are checked by more than one. Every test is
**classifiable** into exactly this scheme — one of four oracle classes (literature,
cross-engine, self-certifying, or live QPA) or the contract/infrastructure
remainder — and since Plan 32 each oracle class is also a **standalone one-liner** a
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

The live QPA suite is `-m qpa` (112 tests). GAP is heavy to install, so it runs in a
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
| Cup / cap / Gerstenhaber bracket | Leibniz sign arbiter + transported-anchor + associativity/commutativity gates |
| Cyclic homology | Connes λ-complex second model + mixed-complex identities |
| The Chouhy–Solotar resolution | CS ≡ bar, CS ≡ Bardzell, and the bank byte-level closed forms |
| Deep degrees past the bar window | bank closed forms + cross-engine + closed-form/chain-count pins |
| Frobenius / Nakayama / symmetry | self-certifying `λ`/`ν` identities + socle criterion |
| `HH` over CC and `GF(p^n)` | exact bar oracle + second-model oracles (field-generic) |
| Distinct-module `Ext(M, N)`, `M ≠ N` | flagged post-v1; self-Ext is the confirmed QPA idiom |
| **The Koszul verdict itself** (Plan 27 — QPA 1.37 has no `IsKoszul`/`KoszulDual`, confirmed by an exhaustive `NamesGVars()` sweep) | the G-quadratic certifier (Priddy PBW: confluent length-2-tip reduction system ⇒ Koszul) plus the generated-in-degree-1 falsifier and the Fröberg matrix identity `P(t)·C_A(−t)=I`; QPA validates every INPUT to the verdict — graded Ext dims and minimal-generator degrees (`ExtAlgebraGenerators`), quadraticity (`IsQuadraticIdeal`), and the quadratic perp (`QuadraticPerpOfPathAlgebraIdeal`) |
| Yoneda relations-by-degree (QPA exposes generator counts, not a presentation) | theory battery: `E(k[x]/x²)=k[y]`, `E(k[x]/xⁿ)=k[y,z]/(y²)` (char-independent, pinned over GF(2)/GF(3)/GF(32003)/char 0), hereditary `E=kQ/J²` with `as_algebra()` round-trip, rad²=0 `E=kQ` path counts, quantum CI `dim Eⁿ=n+1` (= the CS chain count), commutative square `E≅A` self-hosting; the monomial Anick chain-count gate; byte-reproducible lift products (Plan-17-style canonicalization) |

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
`[dev,fast,docs,web,qpa]` extras (2026-07-25, post-merge of Plans 21–26).

| Subsystem (`src/quiverlab/`) | Tests | Bucket | Primary oracle class |
|---|---:|---|---|
| `fields/` (QQ, GF(p), GF(p^n), exact CC = QQ_I) | 37 | fast | exact-arithmetic axioms; base-change invariance |
| `core/` + `combinat/` (Quiver, Algebra, relations, dispatch) | 43 | fast | structure-constant identities; left-to-right path law |
| `groebner/` (overlap completion, admissibility) | 50 | fast | admissibility certificate; finiteness; lowering |
| `hochschild/` (bar, cyclic) | 11 | fast | **the base bar oracle**; mixed-complex identities |
| `engine/` (fast GF(p); minimal, Bardzell, periodic; TT-calculus; cyclic; Coxeter/Nakayama; Plan-29 literature/identity batteries) | 559 | deep | bar oracle; cross-engine; multi-prime; numba/pure parity; frozen QPA-literature values |
| `resolutions_cs/` (CS; comparison; diagonal; cup; cap; Plan-29 literature batteries) | 213 | deep | CS ≡ bar, CS ≡ Bardzell; bank byte-level; literature pins; `d∘d=0` / order; Leibniz + cap identities (unit/module/transport anchors); canonicalization |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity; Plan-29 Tor; Plan-30 Krull–Schmidt decomposition) | 235 | deep | AR/duality literature pins (ASS2006); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks |
| `invariants/` (Cartan, Coxeter, spectral, Betti, cyclic, Frobenius incl. the Plan-29 trace-form symmetry certifier, scalar, sweep; Plan-29 Coxeter/identity literature batteries) | 112 | fast | second models (λ-complex, relative-Tor Betti); self-certifying `λ`/`ν`; GF(p) engine parity |
| `families/` (catalog, zoo; Plan-29 trivial-extension/incidence batteries; Plan-31 certified trivial-extension presentation, `test_trivial_extension_presented.py`) | 139 | deep | closed-form family pins; zoo diversity gates; citations; Plan-31 special-case + Cartan + iso-invariance + CS≡bar pins |
| `batch/` (labdb port, open-zone scans) | 11 | deep | labdb port equality; scan-surface checks |
| `citations/` (registry, bibliography) | 12 | fast | packaged-bib resolution; result references |
| `trace/` (worked-steps incl. the Plan-30 module events + kA₂ replay golden) | 60 | fast | golden-file equality (dims derived from ranks) |
| `viz/` (draw, tikz) | 18 | fast | exact `int`/`Fraction` layout; TikZ |
| `qpa/` (GAP/QPA crosscheck) | 119 | 112 qpa + 7 fast | **live GAP/QPA** (HH dims, self-Ext, τ/τ⁻, proj/inj resolutions, inj dim, Plan-31 native trivial-extension construction — left side via `A^op`); script builders + guards run without GAP |
| `webapp/` (server tier + result cache + offline GUI — non-algebraic glue) | 313 | fast | API / schema / cache canonicalizer (replay-safety rests on exactness) / isolation / artifacts; all math delegated to the library; Plan-28 runner delegation pinned **byte-identical** (frozen goldens + unchanged `canonical_key`) |
| `hpc/` (headless CLI + spec core + container assets — non-algebraic glue, Plan 28) | 53 | fast (checkpoint-resume: deep) | **CLI ≡ public-API parity** on fixture configs; renderer golden tokens (LaTeX/HTML/text ladder); checkpoint-resume end-to-end equals the uninterrupted run; import-boundary + exit-code contract; sbatch/Dockerfile/workflow asset gates |
| `docs/gui/` (Pyodide GUI + no-code module panel — non-algebraic glue) | 61 | fast | runner artifacts / invariants; build hook; freshness |
| release + top-level (`test_no_floats`, `test_errors`, `test_quickstart`; the Plan-32 `test_oracle_classes` audit gate) | 53 | fast (audit gates: deep) | **float-ban AST gate**; error taxonomy; packaging; docs-nav coverage; **oracle-class count audit** (page == live collection) |

Non-algebraic glue (`webapp/`, `docs/gui/`) carries no oracle *because it holds no
mathematics of its own* — it calls `import quiverlab` and is tested for correct
plumbing, not for algebra.

## Buckets and the CI matrix

Test buckets are auto-assigned by directory in `tests/conftest.py` (an explicit
marker wins); the partition is disjoint and exhaustive, enforced by a partition
test. Markers (`pyproject.toml`): `fast`, `deep`, `slow` (implies `deep`), `qpa`;
plus the **orthogonal** oracle-class markers below (which never change a bucket).

| Bucket | Tests | Runs where |
|---|---:|---|
| `fast` | 804 | every CI cell: `{ubuntu, macos, windows} × py{3.10, 3.11, 3.12, 3.13}` |
| `deep` | 1183 | one Linux · py3.12 cell, **twice**: numba and pure (`QUIVERLAB_NO_NUMBA=1`) |
| `qpa` | 112 | weekly Linux · py3.12 job with GAP + QPA (`QUIVERLAB_REQUIRE_QPA=1`) |
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
  `dim Tor_n = dim Ext^n`, the `kZ_n/J^L` symmetry classification, Künneth), and the
  read-only bank's closed-form differentials. The marker face of **Class 1**.
- **`oracle_crossengine`** — two *independent* implementations are run and required
  to agree live: CS ≡ bar ≡ Bardzell ≡ minimal degreewise, numba ≡ pure and
  sparse ≡ dense parity, presented ≡ ⋉ iso-invariance, native ≡ transported
  cup/cap, generic-Domain ≡ GF(p) engine, and the Connes λ-complex second model.
  The library-internal face of **Class 2**.
- **`oracle_selfcert`** — an internal mathematical certificate *is* the assertion:
  `d∘d = 0`, the CS order condition, canonicalization / adversarial-solver
  byte-reproducibility, dimension and iso certificates, the self-certifying
  Nakayama `λ`/`ν` identities, and the unit/Leibniz/module identities that arbitrate
  a sign convention. (These are the "self-certifying internal identities" of Class 1,
  surfaced as their own runnable class.)
- **`qpa`** — the existing bucket marker *is* the fourth oracle class: our value ≡
  live GAP/QPA. It needs no new marker; the live-QPA face of **Class 2**.

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
| Literature / theory pins | `-m oracle_literature` | 670 | the engine reproduces a value/identity that exists outside the library |
| Cross-engine agreement | `-m oracle_crossengine` | 396 | two independent implementations compute the same thing and match live |
| Self-certifying certificates | `-m oracle_selfcert` | 604 | an internal axiom (d∘d=0, canonicality, an arbitration identity) holds by construction |
| Live QPA / GAP | `-m qpa` | 112 | an independent external system (QPA) recomputes and agrees |
| Any oracle class (union) | `-m "oracle_literature or oracle_crossengine or oracle_selfcert or qpa"` | 1227 | the test is pinned by at least one oracle (the remaining tests are contract/infrastructure) |

Collected 2026-07-26 (Plan 32). The oracle markers live only on the pure-library
`engine` / `resolutions_cs` / `hochschild` / `modules` / `invariants` / `families` /
`batch` suites, so these counts do **not** depend on the `[web]`/`[hpc]` extras.

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
- The `webapp/` and `docs/gui/` tiers are verified as software (plumbing,
  isolation, artifacts), not as mathematics — they compute nothing themselves.
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
  in the `.tex`; larger objects render as stated shape+rank elisions, never
  silent omissions; the `.tex` source itself is a served artifact.
- The Plan-28 container tier: what pytest verifies is the **wheel-side story**
  (CLI ≡ public-API parity, renderer goldens, checkpoint-resume, byte-stable
  runner delegation, asset-file gates) plus the CI image smoke (build → run a
  tiny config → render → text-extract, on every tagged release). **Real
  Apptainer on a real cluster is a manual release-checklist step**; the local
  drac-local emulator exercises only the no-container (venv-fallback)
  orchestration path, and `--mem`/OOM behaviour is validated by the host
  `deepen` memory-guard tests, not by the emulator (which records but does not
  enforce memory).

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
- `happel_question` — Happel, D. (1989). Hochschild cohomology of
  finite-dimensional algebras. *Lecture Notes in Mathematics* 1404, 108–126.
- `happel_trivial_extension` — Happel, D. (1988). *Triangulated Categories in
  the Representation Theory of Finite Dimensional Algebras.* London
  Mathematical Society Lecture Note Series 119, Cambridge University Press.
  (The trivial extension `T(A) = A ⋉ DA` is symmetric for every
  finite-dimensional `A`; the repetitive-algebra framework.)
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
