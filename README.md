# quiverlab

**Quivers with relations and Hochschild theory, exactly, for algebraists.**

quiverlab computes with finite-dimensional algebras presented by quivers with
relations, over the complex numbers (exactly — no floating point, ever) and
over all finite fields. Floats fail loudly by design.

## Install (development preview)

```bash
pip install -e '.[dev]'
```

## Three lines to a Hochschild table

```python
from quiverlab import Quiver, CC, GF

Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
A = Q.algebra(relations=["a*b"], field=CC)
print(A.hochschild_cohomology(3))
```

Change one word to move the characteristic:

```python
A2 = Q.algebra(relations=["a*b"], field=GF(2))
print(A2.hochschild_cohomology(3))
```

## The classic characteristic pathology, in one loop

```python
from quiverlab import truncated_polynomial, CC, GF

for field in (CC, GF(2), GF(3)):
    print(field, truncated_polynomial(2, field=field).hochschild_cohomology(4).dims)
# CC     [2, 1, 1, 1, 1]
# GF(2)  [2, 2, 2, 2, 2]
# GF(3)  [2, 1, 1, 1, 1]
```

## General quivers with relations (kQ/I)

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3, 4],
           arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
A = Q.algebra(relations=["a*b - c*d"], field=CC)   # commutative square, exact
print(A.dim)                                        # 9
print(A.hochschild_cohomology(1))                   # HH^0 = 1  HH^1 = 0
```

Non-monomial relations are completed with an exact noncommutative Gröbner
(Buchberger–Mora overlap) engine and certified finite-dimensional; a
non-admissible or infinite presentation fails loudly with `AdmissibilityError`
or `NotFiniteDimensionalError`, never a hang.

## Modules and invariants

```python
from quiverlab import Quiver, CC

A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
           ).algebra(relations=["a*b - c*d"], field=CC)   # commutative square

S1, S4 = A.simple(1), A.simple(4)
A.projective(1).dimension_vector()      # {1: 1, 2: 1, 3: 1, 4: 1}
A.ext(S1, S4, 2)                        # 1     (Ext^2 of simples)
int(A.global_dimension())              # 2
A.loewy_length()                       # 3
A.simple(1).projective_resolution(4)   # P_1 <- P_2(+)P_3 <- P_4 <- 0
```

Every module is a right A-module over the stated exact field; Ext, Hom, and the
projective resolution are exact. Exact `spectral_radius`/`mahler_measure`, `center()`,
`complexity()` (a lower-bound estimate — can under-report, exact only on local /
single-vertex inputs), and `sweep()` (invariant × field) round out the invariant surface.

## Families and citations

```python
from quiverlab import NakayamaAlgebra, QuantumCI, families, bibliography

A = NakayamaAlgebra([3, 2, 2])          # cyclic Nakayama, dim 7
print(A.hochschild_cohomology(0))       # HH^0 = 1
print(A.citations())                    # ('nakayama', 'assem_book', 'bar')

print(families())                       # the whole v1 catalog with signatures
print(bibliography(A.citations()))      # grouped, annotated references
```

## Status

Engine, module, and families phase (Plans 01–06 delivered, together with the
Plan-04 Chouhy–Solotar resolution). On top of the foundations — monomial presentations,
exact fields, bar-complex Hochschild (co)homology — the hanlab deep engine is now
ported and wired in:

- **A fast GF(p) engine** behind the field interface: `hochschild_cohomology`
  and `hochschild_homology` take `engine="auto" | "bar" | "fast"`. `auto` picks
  the numpy mod-p rank engine over prime fields and the exact bar path
  everywhere else; both agree exactly where both can run. The fast engine still
  builds the exponential bar basis, so it guards its depth loudly (raise
  `max_cells` deliberately) — the depth *unlock* lives in the resolutions below.
- **Deep monomial resolutions.** The minimal (Bardzell) and periodic bimodule
  resolutions reach degrees the bar complex never could — k[x]/(x^a) and cyclic
  Nakayama to depth 40 instantly — and certify structural facts (a finite global
  dimension shows up as vanishing generators), cross-checked exactly against the
  bar oracle over primes {32003, 2, 3, 5} on the overlap range.
- **The Chouhy–Solotar resolution** (`resolutions_cs`, `engine="cs"`). The
  domain-generic CS projective bimodule resolution for admissible kQ/I — its
  HH•/HH^• dimensions and representative (co)cycles reach Hochschild degrees the
  bar oracle cannot, with CS↔bar comparison maps; it specializes to Bardzell's
  minimal resolution on monomial algebras (operation transport is certified
  inside the bar-buildable window).
- **Tamarkin–Tsygan calculus** at the engine level: cup product, cap product,
  and the Gerstenhaber bracket; plus **cyclic homology** (Connes' mixed complex).
- **Invariants:** the integer **Cartan** matrix, the **Coxeter** matrix and its
  characteristic polynomial (all fields, exact via sympy); and, over GF(p), the
  **Nakayama** automorphism with the **Frobenius** and **symmetric** tests
  (loud `FieldError` off a prime field).
- **Modules, scalar invariants, and the exact spectral layer.** Right A-modules
  with exact **Ext**, **Hom**, and minimal **projective resolutions**; the scalar
  invariants **Loewy length**, **center**, and **complexity** (GF(p); the last a
  lower-bound estimate that can under-report, exact only on local / single-vertex
  inputs); and the
  exact **spectral radius** / **Mahler measure** of the Coxeter polynomial as
  sympy algebraic numbers — no floats, ever.
- **Algebra families and citations.** A curated catalog of named families
  (`NakayamaAlgebra`, `QuantumCI`, `ExteriorAlgebra`, `IncidenceAlgebra`,
  `PreprojectiveAlgebra`, `TrivialExtension`, `TensorProduct`, …) with `families()`
  discovery and the `zoo` iterator, each stamped with the literature it comes from;
  `A.citations()` and `bibliography(...)` resolve those keys to grouped, annotated
  references, plus a batch scan surface for family sweeps.

Everything is exact — no floating point, ever — and the full test suite runs
green on both the numba kernel path and the pure-Python path
(`QUIVERLAB_NO_NUMBA=1`).

Honest scope note: the calculus lives at the *engine* level today. A classy
`A.cup(u, v)` on named cohomology classes awaits the cohomology-classes
machinery of a later phase (see `docs/plans/ROADMAP.md`).

Coming next (see `docs/plans/ROADMAP.md`): full operation transport, drawing and
TikZ export, worked-steps PDFs, and an optional QPA backend.

MIT © 2026 Marco Armenta
