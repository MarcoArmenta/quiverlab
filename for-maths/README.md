# quiverlab v0.1.0 — complete worked examples for mathematicians

## Example 2 (`witt-*`): non-trivial Gerstenhaber brackets at high degree

**The algebra.** `A = GF(3)[x]/(x³)` — the modular truncated polynomial
algebra, characteristic dividing the exponent. This is the classical
"restricted Lie theory enters Hochschild cohomology" example: `HH¹(A)` is
the 3-dimensional Jacobson–Witt algebra `W(1;1)`, and — unlike generic
quantum complete intersections, whose `HH` vanishes above degree 2 —
**every** `HH^n(A)` is non-zero, so the Gerstenhaber bracket has room to
act in every degree.

**Headline results** (all exact over `GF(3)`, `quiverlab 0.1.0`):

- `HH^•(A) = HH_•(A) = [3, 3, 3, 3, …]` out to degree 12.
- The Gerstenhaber bracket `[·,·]: HH^p × HH^q → HH^{p+q−1}` on the full
  certified window (total degree ≤ 8, recorded in the output): **30 of the
  36 tables are non-zero**, up to `[HH⁸, HH¹]` and `[HH¹, HH⁸]` — with the
  graded-Lie parity pattern visible (even–even pairs vanish, odd arguments
  act everywhere).
- Cyclic homology climbs the characteristic-3 staircase
  `HC_• = [3, 1, 4, 2, 5, 3, 6, 4, 7]`, Connes `B` ranks `2,0,2,0,…`.
- 55 cup and 55 cap tables (all pairs of total degree ≤ 9).

Files mirror Example 1: `witt-report-v0.1.0.html` (full report),
`witt-worked-steps-v0.1.0.html` (every differential, rank and
justification), `witt-gui-v0.1.0.png` (the GUI set up for this run),
`witt-result/trace-v0.1.0.json`, `witt-quiver.tikz.tex`,
`witt-config.yaml` (the complete input).

## Example 1 (`quantumci-*`): the quantum complete intersection

**The algebra.** The quantum complete intersection
`A = k⟨x,y⟩ / (x², y², xy + 2·yx)` over `k = GF(32003)` — one vertex, two
loops, dim A = 4. Everything below was computed **exactly** (no floating
point anywhere) by the released `quiverlab 0.1.0` (`pip install quiverlab`),
reproducible from `quantumci-config.yaml` via
`quiverlab-hpc run quantumci-config.yaml`.

**Headline results** (all in the report):

- `HH^•(A)` to degree 12: `[2, 2, 1, 0, 0, …]` — the Bergh–Erdmann pattern
  for a quantum CI whose parameter is not a root of unity.
- `HH_•(A)` to degree 12: `[3, 2, 2, 2, …]` (BGMS).
- Cyclic homology `HC_•` to degree 5: `[3, 0, 3, 0, 3, 0]`, with the induced
  Connes differential `B: HH_n → HH_{n+1}` of ranks `2, 0, 2, 0, 2` —
  consistent with the SBI sequence.
- The full Tamarkin–Tsygan product surface: **28 cup-product tables**
  `HH^p ⊗ HH^q → HH^{p+q}` and **28 cap-product tables**
  `HH^p ⊗ HH_n → HH_{n−p}` for all pairs up to total degree 6 — computed
  past the bar resolution's reach on the Chouhy–Solotar resolution — plus
  the Gerstenhaber bracket on its certified window (recorded honestly in
  the result: `window: 4`).
- Cartan matrix, center, Coxeter polynomial, and the projective-dimension
  probe (self-injective ⇒ the resolution never terminates; the report
  states the certified bound rather than guessing `∞`).

**The files.**

| file | what it is |
|---|---|
| `quantumci-gui-v0.1.0.png` | the no-code GUI, set up for exactly this computation |
| `quantumci-report-v0.1.0.html` | the full report: every result block, rendered matrices, citations |
| `quantumci-worked-steps-v0.1.0.html` | the "homework-grade" worked steps: the resolution, every differential, pivot/rank lines, justifications |
| `quantumci-result-v0.1.0.json` | the exact machine-readable results (versioned schema) |
| `quantumci-trace-v0.1.0.json` | the exact event stream behind the worked steps |
| `quantumci-quiver.tikz.tex` | the quiver as TikZ, ready for a paper |
| `quantumci-config.yaml` | the complete input — this is *all* a user writes |

Every number is exact arithmetic over `GF(32003)`; every scope boundary
(e.g. the bar complex's exponential growth, the bracket's certified window)
is stated in the output rather than silently truncated. The same
computation runs in the browser GUI (screenshot), on a laptop CLI, or on an
HPC cluster, from the same config.
