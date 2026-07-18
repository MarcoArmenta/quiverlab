# quiverlab tutorials

Short, executable notebooks for algebraists new to the library. They use only
the public `import quiverlab` surface and run end to end in a few seconds.
Read them in order:

1. [`01-exact-fields.ipynb`](01-exact-fields.ipynb) — the field zoo (`CC`,
   `GF(p)`, `GF(p^n)`, `E(n)`), exact arithmetic, and why floats raise
   `ExactnessError` by design.
2. [`02-quivers-and-algebras.ipynb`](02-quivers-and-algebras.ipynb) — quivers,
   left-to-right path composition, the relation parser, certified finiteness of
   `kQ/I`, and the starter families.
3. [`03-hochschild.ipynb`](03-hochschild.ipynb) — Hochschild (co)homology via
   the normalized bar complex, the `char | n` pathology, hereditary vanishing,
   symmetric-algebra duality, and reading an `HHTable`.

The notebooks are committed with their outputs so they render on GitHub. To
re-execute: `jupyter nbconvert --execute --inplace docs/tutorials/*.ipynb`
(with `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2`).
