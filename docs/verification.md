# How quiverlab is verified

Every shipped feature of quiverlab is unit tested. This page says **how** — with
the highest rigour we can bring to it — and it is honest about the edges: where a
check is a cross-engine agreement, where it is a published number, where a live
external oracle can reach, and where it cannot.

The suite is **1377 tests** (collected with the `[dev,fast,docs,web]` extras). It
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
(`hochschild/bar.py`) is the ground truth every deeper engine is measured against.
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
the library, each with its provenance inline:

- **`k[x]/(x^n)`, classical.** In characteristic 0 (or whenever `n` is invertible)
  `dim HH_0 = n` and `dim HH_i = n-1` for `i ≥ 1` (Loday, *Cyclic Homology*; the
  BACH truncated-polynomial computations). The char-5 pathology on `k[x]/(x^5)`
  (`p | n`) collapses both differentials, giving `HH_i = A` in every degree —
  pinned as the vector `[5, 5, 5, 5, 5, 5, 5]`.
- **Quantum complete intersection** `k⟨x,y⟩/(x², y², yx − 2xy)`
  (Buchweitz–Green–Madsen–Solberg, *Math. Res. Lett.* 2005; Bergh–Erdmann,
  *Algebra & Number Theory* 2008). Homology persists, `HH_• = [3, 2, 2, …]`, while
  cohomology dies from degree 3, `HH^• = [2, 2, 1, 0, 0, …]` in characteristic 0 —
  the homology/cohomology asymmetry the deeper engines exist to detect.
- **Happel's theorem** (Happel, LNM 1404, 1989): a hereditary algebra has
  `HH^i = 0` for `i ≥ 2`. Linear `A_3` is pinned to `HH^• = [1, 0, 0, 0, 0]`.
- **Commutative complete intersection** `k[x,y]/(x², y²)`: `HH_• = 4, 4, 5, 6, 7,
  8, 9, …` (the Künneth square of `k[x]/(x²)`), and being symmetric it satisfies
  `HH^n = HH_n` — an internal cross-check QPA also satisfies.
- **Cyclic Nakayama** self-injective algebras `kZ_n/rad²`: low-degree dims frozen.
- **Gentle algebra** `kQ/(ab, ba)` on the 2-cycle: self-injective, so `HH^•` is
  nonzero in every degree (`[1, 1, 1, 1]`) — which makes the CS↔bar agreement a
  discriminating check rather than a run of zeros (oracle: bar cross-check).

### The read-only bank as a byte-level oracle

`tests/resolutions_cs/test_battery_bank_oracle.py` pins the Plan-04 CS resolution
against the original hanlab bank's *hand-derived closed-form* CS differentials — a
wholly separate implementation of the Chouhy–Solotar formulae (arXiv:1406.2300)
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

- **Cyclic homology** over a field containing `Q` is cross-checked against the
  Connes **λ-complex** second model (Loday, Thm 2.1.5) — unnormalized chains, a
  quotient model, disjoint from both `hochschild/bar.py` and
  `hochschild/cyclic.py` (`tests/invariants/test_cyclic_generic.py`), plus the
  mixed-complex identities `b² = 0`, `B² = 0`, `bB + Bb = 0` over QQ.
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
dimensions are a real cross-check, not a tautology.

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

The live QPA suite is `-m qpa` (5 tests). GAP is heavy to install, so it runs in a
**weekly** CI job, not on every commit — but it is **never silently green**: under
`QUIVERLAB_REQUIRE_QPA=1` an absent or broken QPA is a hard failure of that job,
and locally the tests skip explicitly rather than pass vacuously.

Separately, `tests/engine/test_qpa_reference_validation.py` freezes the values QPA
*would* produce (as published in the literature) and requires the bar engine to
reproduce them — an independent check that runs in the normal matrix with no GAP
present.

### Where QPA cannot be compared — and what covers that ground

QPA's cross-check reaches only Hochschild dims and module self-Ext, over **QQ or a
prime field GF(p)** (number-field CC and `GF(p^n)` are out of QPA scope, and raise
loudly). Everything below is therefore covered by a **theory oracle**, not QPA:

| Feature QPA does not cover | Theory oracle that covers it |
|---|---|
| Cup / cap / Gerstenhaber bracket | Leibniz sign arbiter + transported-anchor + associativity/commutativity gates |
| Cyclic homology | Connes λ-complex second model + mixed-complex identities |
| The Chouhy–Solotar resolution | CS ≡ bar, CS ≡ Bardzell, and the bank byte-level closed forms |
| Deep degrees past the bar window | bank closed forms + cross-engine + closed-form/chain-count pins |
| Frobenius / Nakayama / symmetry | self-certifying `λ`/`ν` identities + socle criterion |
| `HH` over CC and `GF(p^n)` | exact bar oracle + second-model oracles (field-generic) |
| Distinct-module `Ext(M, N)`, `M ≠ N` | flagged post-v1; self-Ext is the confirmed QPA idiom |

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
and the oracle class that guards it. Counts are `pytest --collect-only` in the
worktree venv (`[dev,fast,docs,web]`).

| Subsystem (`src/quiverlab/`) | Tests | Bucket | Primary oracle class |
|---|---:|---|---|
| `fields/` (QQ, GF(p), GF(p^n), exact CC = QQ_I) | 37 | fast | exact-arithmetic axioms; base-change invariance |
| `core/` + `combinat/` (Quiver, Algebra, relations, dispatch) | 43 | fast | structure-constant identities; left-to-right path law |
| `groebner/` (overlap completion, admissibility) | 50 | fast | admissibility certificate; finiteness; lowering |
| `hochschild/` (bar, cyclic) | 11 | fast | **the base bar oracle**; mixed-complex identities |
| `engine/` (fast GF(p); minimal, Bardzell, periodic; TT-calculus; cyclic; Coxeter/Nakayama) | 543 | deep | bar oracle; cross-engine; multi-prime; numba/pure parity; frozen QPA-literature values |
| `resolutions_cs/` (CS; comparison; diagonal; cup) | 164 | deep | CS ≡ bar, CS ≡ Bardzell; bank byte-level; literature pins; `d∘d=0` / order; Leibniz; canonicalization |
| `modules/` (Ext, Hom, resolutions) | 57 | deep | exact Ext/Hom; live QPA self-Ext crosscheck |
| `invariants/` (Cartan, Coxeter, spectral, Betti, cyclic, Frobenius, scalar, sweep) | 53 | fast | second models (λ-complex, relative-Tor Betti); self-certifying `λ`/`ν`; GF(p) engine parity |
| `families/` (catalog, zoo) | 72 | deep | closed-form family pins; zoo diversity gates; citations |
| `batch/` (labdb port, open-zone scans) | 11 | deep | labdb port equality; scan-surface checks |
| `citations/` (registry, bibliography) | 12 | fast | packaged-bib resolution; result references |
| `trace/` (worked-steps) | 46 | fast | golden-file equality (dims derived from ranks) |
| `viz/` (draw, tikz) | 18 | fast | exact `int`/`Fraction` layout; TikZ |
| `qpa/` (GAP/QPA crosscheck) | 11 | 5 qpa + 6 fast | **live GAP/QPA** (HH dims, self-Ext); script builders + guards run without GAP |
| `webapp/` (server tier — non-algebraic glue) | 161 | fast | API / schema / isolation / artifacts; all math delegated to the library |
| `docs/gui/` (Pyodide GUI — non-algebraic glue) | 39 | fast | runner artifacts / invariants; build hook; freshness |
| release + top-level (`test_no_floats`, `test_errors`, `test_quickstart`) | 49 | fast | **float-ban AST gate**; error taxonomy; packaging; docs-nav coverage |

Non-algebraic glue (`webapp/`, `docs/gui/`) carries no oracle *because it holds no
mathematics of its own* — it calls `import quiverlab` and is tested for correct
plumbing, not for algebra.

## Buckets and the CI matrix

Test buckets are auto-assigned by directory in `tests/conftest.py` (an explicit
marker wins); the partition is disjoint and exhaustive, enforced by a partition
test. Markers (`pyproject.toml`): `fast`, `deep`, `slow` (implies `deep`), `qpa`.

| Bucket | Tests | Runs where |
|---|---:|---|
| `fast` | 504 | every CI cell: `{ubuntu, macos, windows} × py{3.10, 3.11, 3.12, 3.13}` |
| `deep` | 868 | one Linux · py3.12 cell, **twice**: numba and pure (`QUIVERLAB_NO_NUMBA=1`) |
| `qpa` | 5 | weekly Linux · py3.12 job with GAP + QPA (`QUIVERLAB_REQUIRE_QPA=1`) |
| `slow` | 0 | opt-in (`-m slow`); rides the deep leg |

The `lint` CI job runs the float-gate and release-metadata tests standalone. The
docs site is built `--strict` in its own workflow, so any internals chapter or
page missing from the nav fails the build.

## The standing rule

**Every future plan adds its new oracles to this page as part of its acceptance** —
exactly as every plan already updates the "Under the hood" internals chapters. When
a plan ships a new engine, invariant, or operation, its acceptance task extends the
tables above with the oracle that guards it and the test file that runs it. This
page is the single living record of how each shipped feature is verified, and it is
kept honest: if a subsystem lacks an oracle, this page says so rather than implying
one.

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
