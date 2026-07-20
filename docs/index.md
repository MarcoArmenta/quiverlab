# quiverlab

**Quivers with relations and Hochschild theory, exactly, for algebraists.**

Finite-dimensional algebras `kQ/I` over the complex numbers (exactly — no floating
point, ever) and all finite fields: certified finite-dimensionality, Hochschild
(co)homology with cup products and Gerstenhaber brackets, the first full
Chouhy–Solotar resolution, module Ext, Cartan/Coxeter invariants, drawings, and
worked-steps traces.

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
A = Q.algebra(relations=["a*b"], field=CC)
print(A.hochschild_cohomology(3))
```

- **Tutorials** — start here (executable notebooks).
- **Under the hood** — how each object is represented and each number produced.
- **API Reference** — every public function and class.
- **Web GUI** — compute without installing anything (Plan 09).
- **Cite** — see the JOSS paper and `CITATION.cff`.
