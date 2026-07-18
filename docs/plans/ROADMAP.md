# quiverlab v1 — Implementation Roadmap

Spec: `docs/specs/2026-07-12-quiverlab-design.md` (approved 2026-07-12).
Each plan below produces working, testable software on its own and ends with all
tests green. Plans are written when their phase starts (not speculatively); this
roadmap fixes scope and interfaces between them.

| # | Plan | Delivers | Spec sections |
|---|------|----------|---------------|
| 01 | **Foundations** (`2026-07-12-plan-01-foundations.md`, written) | Package skeleton; exact fields CC/GF(p)/GF(p^n) + generic exact linear algebra; Quiver + relation parser; structure-constant `Algebra` (unit-adapted); monomial `kQ/I` with certified finiteness; normalized bar complex → HH^n/HH_n dims; two starter builders; float-ban AST gate | §3.2, §3.3 (monomial), §4, §5 (components 1, 2, 4; bar of 5), §7 |
| 02 | hanlab engine port | Fast GF(p) kernel stack (pure/sparse/numba, equality-gated) behind the Domain interface; `minimal` A^e resolution + guards + checkpointed deepen; `bardzell`; Tamarkin–Tsygan calculus (cup/cap/bracket) + cyclic homology + bimodule coefficients; Cartan/Coxeter/Nakayama; hanlab test suite travels along. Source is read from the read-only bank `HomologicalAlgebra/HansConjecture/` (copy with attribution; never modify the bank) | §5 (1, 5, 6, 8 partial), §8 ring 1 |
| 03 | Gröbner + general kQ/I | Noncommutative overlap completion, degree bound, admissibility certificate; general relations lower to `Algebra`; reduction systems emitted for CS | §5 (component 3), §3.3 (general) |
| 04 | Chouhy–Solotar + operation transport | The general CS resolution; comparison morphisms/homotopy liftings; CS validation battery (CS≡Bardzell on monomial; CS≡bar/minimal degreewise; literature oracle) | §6, §8 battery |
| 05 | Modules + invariants | `module_ext` generalized (any vertex set, any domain); simples/projectives/injectives, Hom/End, Ext^n; gl.dim, remaining invariants; `sweep` | §5 (7, 8), §3.5–§3.6, §3.9 |
| 06 | Families + batch | Full §3.4 catalog, `zoo()`, `families()`; `quiverlab.batch` (labdb lift) | §5 (9), §3.4 |
| 07 | Viz + trace | `draw()` (matplotlib hard dep added here) + `tikz()`; the worked-steps trace subsystem (PDF/HTML/text renderers, `verbose=True` default, eliding rules, golden-file tests) | §5 (10, 11), §3.7, §3.8, D9 |
| 08 | QPA extra + release | `[qpa]` extra (passagemath-gap), crosscheck oracle + CI job; GitHub Actions matrix; docs site + CI-run tutorials; PyPI packaging; JOSS paper draft | §5 (12), §8 ring 3, §9–§11 |

Standing constraints for every plan: exact arithmetic only, floats fail loudly
(AST gate enforces); read-only banks are never modified; ≤2 parallel agents during
execution; path composition is left-to-right (`a*b` = first `a`, then `b`,
requiring `target(a) = source(b)` — Assem–Simson–Skowroński convention); internal
currency is the unit-adapted structure-constant `Algebra`.

Interface freeze between plans: later plans consume `Domain`
(`coerce/add/sub/neg/mul/inv/is_zero/eq/characteristic`), `linalg.rank/nullspace/rref/solve`,
`Quiver`, `Relation`, `Algebra` (`.domain, .dim, .T, .unit, .multiply, .unit_adapted,
.basis_labels`), and `HHTable` exactly as defined in Plan 01.

**Standing constraint added 2026-07-18 (Marco):** `docs/internals/` — "under the
hood for algebraists" — explains at low coding level how each object is
represented and how each computation runs (chapters for Plans 01–02 exist; the
format is fixed there). Every subsequent plan's acceptance task adds/updates the
chapters for what that plan introduces (Plan 03: groebner + ReductionSystem;
Plan 04: CS resolution + comparison maps; Plan 05: modules — module
representation and module resolutions; Plan 06: families/batch; Plan 07:
viz/trace; Plan 09: webapp). Plans 03/04/09 were committed before this rule:
their executors receive the obligation at dispatch. Chapters ship on the Plan 08
docs site as an "Under the hood" section.
