# How quiverlab is verified

Every shipped feature of quiverlab is unit tested. This page says **how** — with
the highest rigour we can bring to it — and it is honest about the edges: where a
check is a cross-engine agreement, where it is a published number, where a live
external oracle can reach, and where it cannot.

The suite is **1583 tests** (collected with the `[dev,fast,docs,web,qpa]` extras,
2026-07-25 post-merge of Plans 21–26). It
is not a pile of smoke tests: the mathematics is pinned by **two classes of
oracle**, and most numbers are checked by more than one.

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

The live QPA suite is `-m qpa` (31 tests). GAP is heavy to install, so it runs in a
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
| `engine/` (fast GF(p); minimal, Bardzell, periodic; TT-calculus; cyclic; Coxeter/Nakayama) | 543 | deep | bar oracle; cross-engine; multi-prime; numba/pure parity; frozen QPA-literature values |
| `resolutions_cs/` (CS; comparison; diagonal; cup; cap) | 197 | deep | CS ≡ bar, CS ≡ Bardzell; bank byte-level; literature pins; `d∘d=0` / order; Leibniz + cap identities (unit/module/transport anchors); canonicalization |
| `modules/` (Ext, Hom, resolutions; `A^op`, `D`, τ/τ⁻, injectives, left/right sides; Plan-27 Yoneda Ext-algebra + Koszulity) | 188 | deep | AR/duality literature pins (ASS2006); functorial self-certification (`D∘D`, `(A^op)^op`, `τ⁻τ`); live QPA τ/resolutions/inj-dim crosschecks; Yoneda 7-oracle battery (Priddy/Fröberg/Polishchuk–Positselski-cited) + monomial Anick gate + live `ExtAlgebraGenerators`/`IsQuadraticIdeal` crosschecks |
| `invariants/` (Cartan, Coxeter, spectral, Betti, cyclic, Frobenius, scalar, sweep) | 53 | fast | second models (λ-complex, relative-Tor Betti); self-certifying `λ`/`ν`; GF(p) engine parity |
| `families/` (catalog, zoo) | 72 | deep | closed-form family pins; zoo diversity gates; citations |
| `batch/` (labdb port, open-zone scans) | 11 | deep | labdb port equality; scan-surface checks |
| `citations/` (registry, bibliography) | 12 | fast | packaged-bib resolution; result references |
| `trace/` (worked-steps) | 46 | fast | golden-file equality (dims derived from ranks) |
| `viz/` (draw, tikz) | 18 | fast | exact `int`/`Fraction` layout; TikZ |
| `qpa/` (GAP/QPA crosscheck) | 38 | 31 qpa + 7 fast | **live GAP/QPA** (HH dims, self-Ext, τ/τ⁻, proj/inj resolutions, inj dim — left side via `A^op`); script builders + guards run without GAP |
| `webapp/` (server tier + result cache — non-algebraic glue) | 227 | fast | API / schema / cache canonicalizer (replay-safety rests on exactness) / isolation / artifacts; all math delegated to the library |
| `docs/gui/` (Pyodide GUI + no-code module panel — non-algebraic glue) | 51 | fast | runner artifacts / invariants; build hook; freshness |
| release + top-level (`test_no_floats`, `test_errors`, `test_quickstart`) | 50 | fast | **float-ban AST gate**; error taxonomy; packaging; docs-nav coverage |

Non-algebraic glue (`webapp/`, `docs/gui/`) carries no oracle *because it holds no
mathematics of its own* — it calls `import quiverlab` and is tested for correct
plumbing, not for algebra.

## Buckets and the CI matrix

Test buckets are auto-assigned by directory in `tests/conftest.py` (an explicit
marker wins); the partition is disjoint and exhaustive, enforced by a partition
test. Markers (`pyproject.toml`): `fast`, `deep`, `slow` (implies `deep`), `qpa`.

| Bucket | Tests | Runs where |
|---|---:|---|
| `fast` | 584 | every CI cell: `{ubuntu, macos, windows} × py{3.10, 3.11, 3.12, 3.13}` |
| `deep` | 968 | one Linux · py3.12 cell, **twice**: numba and pure (`QUIVERLAB_NO_NUMBA=1`) |
| `qpa` | 31 | weekly Linux · py3.12 job with GAP + QPA (`QUIVERLAB_REQUIRE_QPA=1`) |
| `slow` | 0 | opt-in (`-m slow`); rides the deep leg |

The `lint` CI job runs the float-gate and release-metadata tests standalone. The
docs site is built `--strict` in its own workflow, so any internals chapter or
page missing from the nav fails the build.

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
