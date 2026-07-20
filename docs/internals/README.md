# quiverlab under the hood

These chapters explain *how* quiverlab represents quivers, algebras, and Hochschild
theory inside the computer, and *how* each number it prints is actually produced. They
are written for algebraists who know the mathematics and want to trust (or audit, or
extend) the code, but who do not spend their days programming.

## How to read these chapters

You do not need to know Python. Wherever a Python construct appears we name it in one
clause the first time it shows up. The four you will meet constantly:

- a **list** — an ordered sequence written with square brackets, `[a, b, c]`; you index
  it as `L[0]`, `L[1]`, ... (counting starts at 0, not 1);
- a **tuple** — the same idea but *immutable* (it cannot be edited after creation) and
  written with round brackets, `(a, b, c)`; quiverlab uses tuples for things that must be
  usable as dictionary keys, such as paths;
- a **dict** (dictionary) — a lookup table from keys to values, written `{key: value}`;
  `arrows["a"]` retrieves the value stored under the key `"a"`;
- a **class** — a blueprint bundling data together with the operations on it; an
  *instance* of the class `Quiver` is one particular quiver.

Every chapter has the same five sections:

1. **The mathematics** — one paragraph, the object as you already know it.
2. **How it is represented** — the actual in-memory data, written out for a small example.
3. **How the computation runs** — a numbered walk keyed to the real function names.
4. **A worked micro-example** — one mathematical object, its data, one computation step.
5. **Where to look in the code** — a table of concept -> file -> function/class.

We cite files and symbol names but never line numbers (they drift as code changes).
Every representational claim here was checked by reading the source and, where marked,
by running the code.

## The chapters

- **01 — Exact fields.** The `Domain` protocol; how QQ, GF(p), GF(p^n), and exact CC
  represent field elements; why floats are structurally forbidden; exact linear algebra.
- **02 — Quivers and relations.** Vertices, arrows, paths read left-to-right, the
  relation parser, and the finiteness certificate (the forbidden-word automaton).
- **03 — The algebra.** The structure-constant `Algebra`: basis labels, the
  multiplication table `T`, the unit vector, unit-adaptation, and how kQ/I becomes a table.
- **04 — Hochschild via the bar complex.** What a cochain is in memory, how the
  coboundary matrix is filled entry by entry, and how dimensions fall out of ranks.
- **05 — Resolutions.** The `Resolution` protocol; the minimal A^e resolution by
  syzygies; the Bardzell resolution for monomial algebras; checkpointed deepening; the
  fast GF(p) engine and its pure-Python fallback.
- **06 — Invariants.** Cartan and Coxeter matrices, the exact charpoly, and the
  engine-backed GF(p) extras (Nakayama automorphism, Frobenius/symmetric tests).
- **07 — Dispatch.** How `engine="auto"|"bar"|"fast"` chooses a path, and the
  cross-check philosophy that two independent implementations gives correctness.
- **08 — Gröbner.** General (non-monomial) relations: how a relation becomes a
  reduction rule, Buchberger–Mora completion, the finiteness certificate
  (2L−1 ≤ D plus the forbidden-word automaton), and how a general kQ/I becomes
  structure constants.
- **09 — The Chouhy–Solotar resolution.** The general CS projective bimodule resolution
  for admissible kQ/I: the ambiguity S-sequence (Bardzell associated paths of the tip
  monomial algebra), the CS differential (Bardzell leading map + order-condition-pinned
  correction), the two collapse maps, and why it reaches Hochschild degrees the bar
  complex cannot. Landed with Plan 04.
- **10 — Modules.** How a finite-dimensional right A-module is stored (the `action`
  dict, column-vector / anti-homomorphism convention), simples / projectives /
  injectives read off the multiplication table, radical / top / socle, Hom and Ext,
  and the minimal projective resolution by iterated projective covers.
- **11 — Families and citations.** The v1 family catalogue (Nakayama through
  trivial extension), the three construction routes (monomial / general /
  structure-constant), the curated zoo, and the citation registry that stamps
  every algorithm and family with the paper behind it.
- **12 — Visualization and worked-steps traces (viz + trace).** The exact
  `int`/`Fraction` layered layout shared by `A.draw()` and `A.tikz()`, the typed
  trace-event taxonomy, the recorder and its elision guards, and how a verbose
  computation renders a LaTeX→PDF (or self-contained no-JS HTML) worked-steps
  document whose dimensions are *derived* from the recorded ranks. Landed with Plan 07.

## Honest coverage statement

This tree carries **Plans 01–06** together with the **Plan-04 Chouhy–Solotar
resolution**. What is documented here is what is on disk *now*:

- The **fast GF(p) engine** (`engine/`) is ported and live: bar homology/cohomology over
  a prime field, the minimal and Bardzell resolutions, cyclic homology, and the
  Coxeter/Nakayama layer.
- **General (non-monomial) relations** are built (Chapter 08): `Quiver.algebra` routes a
  monomial presentation through the Plan-01 path and a non-monomial one through the
  noncommutative Gröbner engine (`groebner/`), which completes the relations, certifies
  finite-dimensionality, and lowers kQ/I to a structure-constant `Algebra`.
- The **Chouhy–Solotar resolution** (`resolutions_cs`) has **landed with Plan 04**
  (Chapter 09 — The Chouhy–Solotar resolution): the domain-generic CS bimodule
  resolution, its HH•/HH^• dimensions and representative (co)cycles, the CS↔bar
  comparison maps, and the `engine="cs"` dispatch path (which reaches Hochschild
  degrees the bar oracle cannot).
- The **module + invariants surface** has **landed with Plan 05**: right A-modules with
  simples / projectives / injectives, radical / top / socle, Hom / Ext, minimal
  projective resolutions and global dimension (Chapter 10), plus the exact
  **spectral-radius / Mahler-measure** layer, `loewy_length`, `center`, `complexity`, and
  the `sweep` (invariant × field) table (documented in the "exact spectral layer" section
  of Chapter 06).
- The full algebra **family catalogue** (`families/`) and the **citation registry**
  (`citations/`) have **landed with Plan 06** (Chapter 11): Nakayama through trivial
  extension, the three construction routes (monomial / general / structure-constant),
  the curated zoo, and `bibliography()`.
- The **visualization and worked-steps trace** surface (`viz/`, `trace/`) has
  **landed with Plan 07** (Chapter 12): the exact `int`/`Fraction` layout shared by
  `A.draw()` and `A.tikz()`, the typed trace-event taxonomy and bounded recorder, and
  the LaTeX→PDF (else self-contained no-JS HTML) worked-steps documents whose
  dimensions are derived from the recorded ranks and cite through `bibliography()`.

Where a chapter describes something whose full form is still to come, it says so inline.

## Release infrastructure and the QPA cross-check

Plan 08 (release) adds **no internals chapter**: CI, packaging, the docs site, and
the JOSS paper are infrastructure, not core mathematics. The one mathematical addition
— the optional QPA cross-check that recomputes Hochschild dimensions via the
enveloping algebra `HH^n(A) = Ext^n_{A^e}(A,A)` — is documented on the docs
[Development page](../development/release.md), not here.
