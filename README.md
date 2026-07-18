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

## Status

Engine phase (Plans 01–02 delivered). On top of the foundations — monomial
presentations, exact fields, bar-complex Hochschild (co)homology — the hanlab
deep engine is now ported and wired in:

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
- **Tamarkin–Tsygan calculus** at the engine level: cup product, cap product,
  and the Gerstenhaber bracket; plus **cyclic homology** (Connes' mixed complex).
- **Invariants:** the integer **Cartan** matrix, the **Coxeter** matrix and its
  characteristic polynomial (all fields, exact via sympy); and, over GF(p), the
  **Nakayama** automorphism with the **Frobenius** and **symmetric** tests
  (loud `FieldError` off a prime field).

Everything is exact — no floating point, ever — and the full test suite runs
green on both the numba kernel path and the pure-Python path
(`QUIVERLAB_NO_NUMBA=1`).

Honest scope note: the calculus lives at the *engine* level today. A classy
`A.cup(u, v)` on named cohomology classes awaits the cohomology-classes
machinery of a later phase (see `docs/plans/ROADMAP.md`).

Coming next (see `docs/plans/ROADMAP.md`): general relations via noncommutative
Groebner bases, the first full Chouhy–Solotar resolution and operation
transport, module Ext and the remaining invariants, family catalogs and batch,
drawing and TikZ export, worked-steps PDFs, and an optional QPA backend.

MIT © 2026 Marco Armenta
