# quiverlab

**Quivers with relations and Hochschild theory, exactly, for algebraists.**

Draw a quiver below, type relations, pick a field, and compute — right here,
with nothing to install. Everything runs **exactly** (never floating point) in
your browser, on the same engine the Python library ships.

<div id="qlgui">
  <noscript><p><strong>The interactive GUI needs JavaScript.</strong>
  The Python library works without it: <code>pip install quiverlab</code>.</p></noscript>
  <p>Loading the GUI…</p>
</div>

## Prefer code?

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

- **Tutorials** — start here (executable notebooks):
  [1 — Exact fields](tutorials/01-exact-fields/) ·
  [2 — Quivers and algebras](tutorials/02-quivers-and-algebras/) ·
  [3 — Hochschild theory](tutorials/03-hochschild/).
- **[Under the hood](internals/)** — how each object is represented and each
  number produced.
- **[API Reference](reference/)** — every public function and class.
- **Web GUI** — the form at the top of this page runs in your browser; a
  server-backed tier for big jobs is planned (Plan 09).
- **Cite** — see the JOSS paper and
  [`CITATION.cff`](https://github.com/MarcoArmenta/quiverlab/blob/main/CITATION.cff).
