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

Foundations phase (Plan 01). Monomial presentations, exact fields, bar-complex
Hochschild (co)homology. Coming next (see `docs/plans/ROADMAP.md`): the hanlab
deep engines, general relations via noncommutative Groebner bases, the first
full Chouhy–Solotar resolution, cup products and Gerstenhaber brackets, module
Ext, drawing and TikZ export, worked-steps PDFs, and an optional QPA backend.

MIT © 2026 Marco Armenta
