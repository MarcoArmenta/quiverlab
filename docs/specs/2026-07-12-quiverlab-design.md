# quiverlab — Design Specification (v1)

**Date:** 2026-07-12
**Status:** Approved design, pending Marco's review of this written spec
**Owner:** Marco Armenta
**License:** MIT

---

## 0. Decision log (settled during brainstorming, 2026-07-12)

| # | Decision | Marco's choice |
|---|----------|----------------|
| D1 | Scope vs NN machinery | **Pure algebra only, no NN code at all.** The bridge programme becomes an external consumer. |
| D2 | Ground fields | **ℂ and all finite fields.** ℂ realized by exact algebraic arithmetic (see §4.2); all GF(p^n) exact. |
| D3 | Float policy | **Loud failure on any floating-point contamination.** No silent approximation anywhere, ever. |
| D4 | Viewing the quiver | **`draw()` (PNG/SVG) and TikZ export are v1 must-haves.** Plain-text `repr` included as a freebie. |
| D5 | v1 depth | **v1 must include the full general Chouhy–Solotar resolution** (not just monomial Bardzell). |
| D6 | Deliverable | **PyPI + documentation site + tutorial notebooks + citable software paper (JOSS).** |
| D7 | Architecture | **Hybrid:** pure-Python MIT core; optional GPL QPA backend via `pip install quiverlab[qpa]`. |
| D8 | Name | **quiverlab** (PyPI availability to be verified at implementation start; fallback `quiver-lab`). |
| D9 | Worked-steps traces | **`verbose` defaults to True**; computations write a step-by-step worked-example document (PDF when a LaTeX toolchain exists, HTML+MathJax fallback) so beginners learn and experts audit. |

---

## 1. Purpose and audience

quiverlab is a pure-Python library for computing with finite-dimensional associative
algebras presented as quivers with relations, over exact fields, aimed at **research
algebraists who barely code**. Three lines from install to a Hochschild cohomology
table; every number exact and certified; every computation optionally accompanied by
a human-readable worked-steps document.

The library unifies and generalizes Marco's existing research software (hanlab and the
bridge-lab obstruction/module engines) into a public tool that is independent of the
neural-network research that produced it.

**Success criteria for v1**

1. `pip install quiverlab` works on macOS (incl. M-series), Linux, Windows with zero
   external systems.
2. A non-coder can build a custom `kQ/I`, draw it, and get HH^0..HH^n with cup products
   and Gerstenhaber brackets from the README example alone.
3. Every capability in §3 is implemented, exact, and covered by the four-ring test
   strategy of §10.
4. JOSS paper submitted.

---

## 2. Prior art and the gap (research findings, 2026-07-12)

Full ecosystem research was performed (web, verified against current docs). Summary:

- **QPA v1.37 (GAP, GPL-2.0, active, 2026-05-13):** the strongest existing system —
  `kQ/I` by admissible ideals (Gröbner bases delegated to GBNP), Green–Solberg–Zacharia
  minimal module resolutions, Ext, AR theory, homological-dimension bounds, rich family
  library and property recognizers. **No packaged Hochschild cohomology** (only
  assemblable by hand via `EnvelopingAlgebra` + module Ext), **no cup product as a
  product, no Gerstenhaber bracket at all.** Install requires GAP (non-trivial for
  non-coders; native Windows historically Cygwin/WSL).
- **QPA2:** dormant work-in-progress, no releases. Never build on it.
- **SageMath:** native `PathAlgebra` is the **free** path algebra only — no
  quotient-by-relations bound-quiver object. `HochschildComplex(A, M)` exists for any
  algebra-with-basis but is unreduced-bar-complex only (exponential; toy sizes), with
  no cup/bracket. GB-scale install.
- **Magma Basic Algebras (Carlson):** strong Ext-algebra machinery; proprietary,
  non-redistributable, no bracket.
- **Macaulay2 / Singular (Letterplace, Plural):** graded/Ore-type noncommutative
  Gröbner engines; not f.d. path-algebra representation theory; no Hochschild product.
- **QuiverTools (Belmans–Franzen–Petrella v1.0):** moduli of quiver representations
  (stability, HN strata); no relations, no higher Hochschild.
- **Python native:** **nothing exists** — no quivers-with-relations, no Hochschild, no
  bracket, anywhere on PyPI/GitHub. SymPy offers only non-commutative symbols.

**The gap quiverlab fills:** no system on earth ships Hochschild cohomology **with its
operations** (cup product, Gerstenhaber bracket) for finite-dimensional algebras; no
system ships the Chouhy–Solotar resolution; nothing at all is pip-installable for
non-coders. quiverlab's differentiators, in order: (1) Hochschild theory *with the
operations*, deep degrees; (2) first full Chouhy–Solotar implementation; (3) trivial
installation + non-coder API; (4) worked-steps trace documents; (5) exact-only
semantics with loud failure.

Key algorithm sources:
- Bardzell, *The alternating syzygy behavior of monomial algebras* (1997) — minimal
  bimodule resolution, monomial case. Already implemented in hanlab.
- Chouhy–Solotar, *Projective resolutions of associative algebras and ambiguities*,
  J. Algebra 432 (2015), arXiv:1406.2300 — general `kQ/I` via reduction systems.
- Negron–Witherspoon, arXiv:1406.0036; Volkov, arXiv:1610.05741 — Gerstenhaber
  bracket on arbitrary resolutions via homotopy liftings / comparison morphisms.
- Green–Solberg–Zacharia (Trans. AMS 2001) — module-level minimal resolutions (QPA's
  algorithm; ours for module Ext).

---

## 3. User-facing API

Design lessons adopted from SnapPy (census strings, one-line objects), GAP SmallGroups
(indexed families), networkx (generators + draw), SymPy/Sage (field as explicit first
argument, exact defaults).

### 3.1 One import

```python
from quiverlab import *
```

Everything a working algebraist needs is exported flat: `Quiver`, field constructors
(`CC`, `GF`), family constructors (§3.4), and helper functions. Power users can import
submodules explicitly.

### 3.2 Fields

Two user-facing fields (D2):

```python
CC          # "the complex numbers", exact — see §4.2
GF(q)       # any prime power: GF(2), GF(7), GF(4), GF(27), ...
```

`field=CC` is the default everywhere. Exact ℂ entries accepted: Python ints,
`Fraction`, exact strings (`"1/3"`, `"i"`, `"sqrt(2)"`), and `E(n)` = primitive n-th
root of unity (GAP convention). Floats raise `ExactnessError` (D3).

### 3.3 Custom algebras: quivers with relations

```python
Q = Quiver(vertices=[1, 2, 3],
           arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
A = Q.algebra(relations=["a*b - c"], field=CC)     # kQ/I
```

- Loops and multiple parallel arrows allowed.
- Relations are strings: linear combinations of parallel paths with exact
  coefficients (`"a*b - 2*c"`, `"x^2"`, `"a*b*a"`). Non-parallel summands are
  rejected loudly with the source/target mismatch shown.
- On construction quiverlab runs noncommutative Gröbner completion (§5, `groebner`)
  and certifies finite-dimensionality; otherwise `NotFiniteDimensionalError` /
  `AdmissibilityError` (never a hang, never a guess).
- Monomial presentations are detected automatically (engine dispatch, §5).
- Escape hatch for experts: `Algebra.from_structure_constants(T, unit, field=...)`.

### 3.4 Families (v1 catalog)

One-liners, field always an explicit keyword with default `CC`:

```python
NakayamaAlgebra([3, 2, 2])              # by Kupisch series
NakayamaAlgebra(n=4, l=3, cyclic=True)  # cyclic quiver, rad^l = 0
PathAlgebra("D5")                       # hereditary Dynkin/Euclidean, any orientation option
TruncatedPathAlgebra("A5", 2)           # kQ / rad^2
RadicalSquareZero(Q)                    # from any quiver
IncidenceAlgebra(poset)                 # from Hasse-diagram data
PreprojectiveAlgebra("A3")
QuantumCI(q="i")                        # quantum complete intersection k<x,y>/(x^2, y^2, xy+q yx)
ExteriorAlgebra(3)
TrivialExtension(A)
TensorProduct(A, B)
zoo(dim_max=12)                         # iterator over a curated exact zoo (from hanlab open_zoo)
families()                              # discoverability: list all of the above with signatures
```

### 3.5 Computations on an algebra

```python
A.dimension();  A.basis()               # irreducible-path basis
A.hochschild_cohomology(6)              # HHTable object: dims of HH^0..HH^6
A.hochschild_homology(6)
A.cyclic_homology(6)                    # Connes mixed complex
A.hh_basis(n)                           # explicit cocycle representatives
A.cup(u, v)                             # cup product of classes
A.bracket(u, v)                         # Gerstenhaber bracket
A.cap(z, u)                             # cap action
A.hochschild_cohomology(6, coefficients=M)   # bimodule coefficients, incl. twisted
A.global_dimension()                    # exact value or certified bound (labeled)
A.cartan_matrix(); A.coxeter_matrix(); A.coxeter_polynomial()
A.nakayama_automorphism(); A.is_frobenius(); A.is_symmetric(); A.is_selfinjective()
A.loewy_length(); A.complexity(n); A.center()
```

`HHTable` prints as a readable table, exposes `.dims`, `.latex()`, and states
loudly which degrees were computed and which engine certified them.

### 3.6 Modules

```python
S = A.simple(1);  P = A.projective(2);  I = A.injective(1)
M.radical(); M.top(); M.socle(); M.dimension_vector()
A.hom(M, N); A.ext(M, N, n)             # module-level Ext via minimal resolutions
M.projective_resolution(5)              # inspectable resolution object
```

### 3.7 Looking at the algebra (D4)

```python
A.draw()                                # matplotlib PNG/SVG; loops, multiplicities, relation list
A.draw(file="A.svg")
A.tikz()                                # TikZ code, same layout coordinates, paste into paper
print(A)                                # plain-text: vertices, arrows, relations
```

### 3.8 Worked-steps traces (D9)

```python
quiverlab.verbose        # True by default
A.hochschild_cohomology(4)
# -> computes, then writes ./quiverlab_traces/HH_A_<hash>.pdf (or .html) and prints:
#    "Worked steps: quiverlab_traces/HH_A_3f2a.pdf (12 pp)"
A.hochschild_cohomology(4, verbose=False)   # per-call off
quiverlab.verbose = False                   # global off
```

Trace contents: the chosen resolution and why (dispatch decision), each term with its
generators, each differential as a matrix over the stated field (elided above a size
threshold with an explicit "elided" note), each rank computation, the resulting
dimensions; for cup/bracket: the representatives, comparison/lifting steps, the
product cocycle and its reduction. Renderers: LaTeX→PDF when `pdflatex`/`tectonic`
is on PATH; otherwise HTML+MathJax with a loud one-line explanation of both facts.
Golden-file tests assert trace claims equal computed values (§10).

### 3.9 Moving variables (D2 "all of the moving variables")

```python
sweep(NakayamaAlgebra, [3, 2, 2], fields=[CC, GF(2), GF(3), GF(5)])
# -> table: invariant × field; the characteristic-dependence view
```

`sweep` rebuilds the algebra over each field and tabulates any requested invariants.
Batch/zoo scans with SQLite persistence (lifted from hanlab `labdb`) ship as
`quiverlab.batch` for experts; not part of the beginner surface.

---

## 4. Semantics and guarantees

### 4.1 Exactness contract (D3)

- No engine has a floating-point code path. CI enforces structurally (no float dtypes
  in core modules; lint gate) plus debug-mode runtime asserts.
- Floats in any user input raise `ExactnessError` with a fix-it message.
- Every returned number is exact over the stated field. Partial computations are
  labeled with their certified range.

### 4.2 What "ℂ" means (D2)

An algebra "over ℂ" must have structure constants in an exact subfield of ℂ:
entries generate a number field ℚ(α) (constructed via sympy's algebraic-number
arithmetic from ints, fractions, `i`, `sqrt`, `E(n)`). All computed invariants
(dims of HH/Ext, gl.dim, Cartan/Coxeter data) are invariant under the flat base
change ℚ(α) → ℂ, so the answers **are** the ℂ answers, provably, with no floats.
This covers every algebra with algebraic structure constants — essentially all of
the literature. Genuinely transcendental entries are out of scope for v1 and fail
loudly with an explanation.

### 4.3 Finite fields

GF(p) uses the hanlab exact mod-p kernel stack. GF(p^n) uses polynomial arithmetic
over GF(p) with a bundled Conway-polynomial table (Lübeck's tables, covering all q
of practical size); a user-supplied irreducible polynomial is accepted beyond the
table; anything else fails loudly.

---

## 5. Architecture

Twelve components, four layers. Internal currency: hanlab's structure-constant
`Algebra` `(m, T, unit)` — unit-adapted basis, generalized from int64/F_p to a
domain-generic exact element type.

**Foundations**
1. `fields` — domain protocol (add/mul/inv/rank/nullspace/solve). Backends: ℚ fast
   path (`Fraction`, Bareiss fraction-free elimination), number fields ℚ(α) (sympy),
   GF(p) (hanlab pure/sparse/numba triple, equality-gated, `QUIVERLAB_NO_NUMBA`
   guard), GF(p^n) (§4.3).
2. `combinat` — `Quiver`, paths/words, relation parser (§3.3), reduction systems.
3. `groebner` — noncommutative overlap completion for path-algebra ideals with degree
   bound + **admissibility certificate** (finite irreducible-word basis ⟹ certified
   f.d.); emits the reduction system consumed by Chouhy–Solotar.

**Core**
4. `core` — `Algebra` (structure constants over a domain); `Quiver.algebra()` lowers
   `kQ/I` via the Gröbner basis (basis = irreducible paths). Unifies hanlab's three
   presentation front-ends (`MonomialPresentation`, `ReductionSystem`, ad-hoc
   builders) into one construction path that always emits an `Algebra`.

**Engines**
5. `resolutions` — one interface, four backends, auto-dispatch with loud reporting:
   `bar` (oracle, small), `minimal` (hanlab A^e syzygy stepper, memory/term guards,
   checkpointed deepen), `bardzell` (monomial, deep; proven to degree 1702 in hanlab),
   `chouhy_solotar` (general relations, deep — §6, **new**). Dispatch: monomial →
   Bardzell; general admissible → CS; bar/minimal always available as cross-checks.
6. `hochschild` — HH^•/HH_• dims and bases from any backend; bimodule/twisted
   coefficients; cyclic homology (Connes B, mixed complex); Tamarkin–Tsygan calculus
   (cup, cap, Gerstenhaber bracket) lifted from hanlab `tt_calculus`, computed on bar
   cochains; comparison maps (homotopy liftings, Negron–Witherspoon/Volkov style)
   connect CS-resolution classes to bar representatives. Operations report the degree
   window in which they are certified.
7. `modules` — right modules, simples/projectives/injectives, radical/top/socle,
   Hom/End, module Ext^n via minimal resolutions. Engine = bridge
   `obstruction/module_ext.py` **generalized** (arbitrary vertex sets, any domain —
   removes the hardcoded 4-vertex diamond) and cross-checked against hanlab.
8. `invariants` — gl.dim (exact when certifiable, else labeled verified bounds),
   Cartan/Coxeter (matrix, polynomial, spectrum, cyclotomicity/Mahler tests from
   hanlab `coxeter_spectrum`), Nakayama automorphism, Frobenius/symmetric/
   selfinjective, Loewy length, complexity, center; `sweep` (§3.9).
9. `families` — §3.4 catalog + `zoo()` (hanlab `open_zoo` lift) + `families()`
   discoverability.

**Surface**
10. `viz` — `draw()`: pure-matplotlib layered layout by default (graphviz optional
    nicety, never required); `tikz()` shares the same computed coordinates; relation
    list rendered under the diagram.
11. `trace` — engines emit typed step events (`ResolutionTerm`, `Differential`,
    `RankStep`, `LiftStep`, `Dispatch`); renderers: LaTeX→PDF, HTML+MathJax, text.
    Default on (D9); size-eliding rules with explicit elision notes.
12. `qpa` (extra) — `pip install quiverlab[qpa]` → passagemath-gap wheels (macOS/
    Linux; loud, graceful message on Windows). Provides: AR quiver of an algebra,
    QPA property recognizers (gentle/string/special-biserial/…), and
    `A.crosscheck(...)` — independent recomputation of module Ext and HH dims (via
    the enveloping-algebra route scripted in GAP) for validation workflows.

Engineering notes carried from the inventory: proper relative imports (hanlab's
flat-namespace hack is not carried over); numba strictly optional (`[fast]` extra,
pure path default, both paths equality-tested); no multiprocessing in core; no
hardcoded absolute paths; thread counts default-throttled.

---

## 6. The Chouhy–Solotar engine (v1 flagship, D5)

Goal: a small projective bimodule resolution for arbitrary admissible `kQ/I`
(specializing to Bardzell's minimal one in the monomial case), enabling deep exact
Hochschild theory beyond monomial algebras — never fully implemented in any system.

Plan (arXiv:1406.2300):
1. From the certified Gröbner/reduction system (§5, component 3): the S-sequence of ambiguity
   chains (generalizing Bardzell's AP_n of overlaps).
2. Resolution terms A ⊗ kS_n ⊗ A; differentials constructed recursively with an
   explicit contracting homotopy (the CS recursion); all over the domain-generic
   exact arithmetic.
3. Comparison morphisms to the bar resolution (homotopy liftings) to transport
   Tamarkin–Tsygan operations (§5, component 6).

Validation battery (all in CI):
- **CS ≡ Bardzell** exactly (terms, ranks, HH dims) on every monomial algebra in the
  test zoo — CS specializes to Bardzell there.
- **CS ≡ bar ≡ minimal** degreewise on a mixed zoo of non-monomial algebras within
  the overlap window those engines can reach.
- **Literature oracle:** published HH^• of quantum complete intersections
  (Buchweitz–Green–Madsen–Solberg), gentle cases, k[x]/(x^n) across characteristics,
  hereditary Dynkin (Happel).

Risk register:
- Ambiguity-chain combinatorics for non-graded algebras is research-grade → v1 gates
  on the admissibility certificate; anything outside raises `NotImplementedError`
  at the exact boundary rather than risking a wrong answer.
- Term growth for wild examples → hanlab-style memory/term guards + checkpointed
  resume; guards fail loudly with the certified range (D3).

---

## 7. Errors

`QuiverlabError` hierarchy, every message mathematician-readable with a fix-it hint:

- `ExactnessError` — float or non-exact input; "write '1/3' or Fraction(1, 3)".
- `AdmissibilityError` / `NotFiniteDimensionalError` — Gröbner completion cannot
  certify a finite basis; shows the offending cycles/irreducible-word growth.
- `RelationError` — non-parallel paths in a relation; shows source/target mismatch.
- `DepthLimitError` — memory/term guard hit; states certified degrees and the
  checkpoint-resume command to push deeper.
- `FieldError` — unsupported field/entry (e.g. transcendental, GF beyond table
  without user polynomial).

Never hang, never guess, never silently truncate (D3).

---

## 8. Testing (four rings + two batteries)

1. **Lifted suite:** hanlab's 602 tests travel with the code; pure/numba dual paths
   equality-gated across primes {2, 3, 5, 32003}.
2. **Literature oracle battery:** published values recomputed and matched, citations
   in the test files (Happel; k[x]/(x^n) small-p pathologies; BGMS quantum CIs;
   gentle/zigzag; preprojective).
3. **QPA cross-check (CI-only, Linux job, dev dependency):** random zoo samples —
   module Ext^n vs QPA; HH dims vs the GAP-scripted enveloping-algebra route.
4. **Property-based laws:** Gerstenhaber graded Jacobi + Leibniz, cup graded
   commutativity, basis-order invariance, field-extension invariance (ℚ vs ℚ(i)),
   HH^n/HH_n dualities for symmetric algebras — randomized small algebras.

Plus: the CS battery (§6), and golden-file trace tests (trace claims ≡ computed
values). Tutorials execute in CI (§11).

---

## 9. Packaging and distribution (D6, D7)

- Pure-Python wheel. Hard deps: `numpy`, `sympy`, `matplotlib`. Extras: `[fast]`
  (numba), `[qpa]` (passagemath-gap; macOS/Linux). Python ≥ 3.10.
- Repo: `/Users/marco/Desktop/quiverlab/`, fresh git, this spec as first commit;
  public GitHub under Marco's account when he says so. hanlab and the banks stay
  frozen and untouched; code is lifted with attribution (both MIT, both Marco's;
  the one lifted bridge engine, `module_ext.py`, gains an MIT header — bridge/ is
  currently unlicensed).
- CI: GitHub Actions — macOS/Linux/Windows × Python 3.10–3.13; docs build; QPA
  cross-check job (Linux); float-free lint gate.
- Versioning: semver; 0.x during battle-testing; 1.0 at JOSS acceptance.

## 10. Documentation

Docs site (mkdocs-material) with auto-generated API reference and
CI-executed tutorial notebooks, non-coder first:

1. Your first algebra in 5 minutes (install → build → draw → HH table).
2. The families tour.
3. Hochschild operations (cup, bracket, cyclic).
4. Characteristic sweeps ("all the moving variables").
5. Reading the worked-steps PDF.
6. Experts: depth, checkpoints, zoo scans, structure-constant escape hatch.

## 11. Paper (D6)

- **v1 milestone: JOSS submission** — refereed, fast, citable; scope = the library.
- Later: a mathematics-audience paper (first full Chouhy–Solotar implementation +
  a headline computation unavailable in any other system, e.g. Gerstenhaber
  structure across a family under a characteristic sweep).

## 12. Non-goals (v1)

Native AR-quiver computation (available via `[qpa]`; native later) · Morita/derived
equivalence · A∞-structures · Ext-algebra presentations · infinite-dimensional
algebras · non-admissible ideals · group cohomology per se · anything neural-network.

## 13. Lift map (inventory → quiverlab)

From hanlab (`…/HansConjecture/hanlab/`), lifted nearly as-is into the named
components: `hh_engine` (bar HH, `Algebra`, `rank_mod_p`) → core/5/6; `scan3`
(cohomology, complexity) → 6/8; `bimodule` (coefficients, twisted, dual) → 6;
`tt_calculus` (cup/cap/bracket) → 6; `cyclic` (Connes B, HC) → 6; `coxeter`,
`coxeter2`, `coxeter_spectrum` (Cartan/Coxeter/Nakayama + mod-p linear algebra)
→ 8/1; `resolutions_bardzell` → 5; `resolutions_minimal` + `deepen` (guards,
checkpoints) → 5; `resolutions_cs` (2 closed-form families) → subsumed by §6;
`reduction_algebra`, `open_zoo` → 4/9; `linalg_fast`, `_kernels` (+ no-numba guard)
→ 1; `labdb` → `quiverlab.batch`. From bridge lab: `obstruction/module_ext.py`
→ 7 (generalized per §5, component 7). Left behind: all NN-coupled code (quiver_obs,
theory_obs, barcode, bimodule/reps, differential_*, mnist1d_falsify, rho_battery,
knowledgematrix — D1).

## 14. Build-order sketch (detail in the implementation plan)

fields + core → port hanlab engines (tests green) → groebner + CS (§6) →
operations bridge (liftings) → modules + invariants → families → viz + trace →
qpa extra → docs + paper.
