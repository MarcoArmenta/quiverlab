# quiverlab

**Quivers with relations and Hochschild theory, exactly, for algebraists.**

Compute with the **containerized application** — a zero-code GUI plus a batch
CLI, fully offline, using your machine's real cores and RAM — or with the
self-hostable **web tier** for shared and queued jobs, or with the Python
library. Everything runs **exactly** (never floating point).

- **[Offline app](offline-app/)** — pull one image, run `gui`, open localhost.
- **[Run on your HPC cluster](hpc/)** — the same image runs batch configs,
  with checkpointed resume and rendered HTML reports.

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
- **No-code interfaces** — the [offline app](offline-app/) and the
  self-hostable web tier (`webapp/` in the repository).
- **Cite** — see the JOSS paper and
  [`CITATION.cff`](https://github.com/MarcoArmenta/quiverlab/blob/main/CITATION.cff).
