# 07 — Dispatch: choosing a computation path

## The mathematics

The Hochschild dimensions are mathematical facts that do not depend on how you compute
them. quiverlab has two genuinely independent implementations of them — the pure,
exact-over-any-field *bar oracle* (Chapter 04) and the fast numpy-int64 *engine* over F_p
(Chapter 05) — plus the specialised resolution backends. The dispatch layer decides which
path runs for a given request, and the design principle is that the *agreement* of two
independent implementations is the correctness argument.

## How it is represented — the `engine` switch

`hochschild_cohomology` and `hochschild_homology` on `Algebra` take a keyword
`engine="auto" | "bar" | "fast"` (`core/algebra.py`). An unknown value raises
`QuiverlabError` immediately. The decision is a single predicate, `_use_fast_engine`:

```
_use_fast_engine(engine) := (engine == "fast")
                            or (engine == "auto" and domain is a PrimeField)
```

Reading off the three cases:

- **`"bar"`** — always the pure bar oracle, over whatever field the algebra uses (QQ,
  GF(p), GF(p^n), or exact CC). This is the universal path.
- **`"fast"`** — always the engine. If the field is not a prime field the engine's adapter
  (`to_engine`, Chapter 05) raises `FieldError` loudly; `"fast"` never silently degrades.
- **`"auto"`** — the engine **exactly when the domain is a `PrimeField`**, otherwise the
  bar oracle.

## Why "auto" picks the fast engine exactly for prime fields

The fast engine is a GF(p) accelerator and nothing else: its whole speed comes from numpy
int64 matrices and modular-inverse rank mod p, which only make sense over F_p. So `"auto"`
routes to it precisely when — and only when — the algebra lives over a `PrimeField`. The
subtlety worth internalising: **GF(p^n) for n ≥ 2 is a `FiniteField`, not a `PrimeField`**
(Chapter 01 — its elements are coefficient tuples, not single ints). It is *not* an instance
of `PrimeField`, so `isinstance(domain, PrimeField)` is false and `"auto"` correctly sends
GF(4), GF(27), ... to the bar oracle rather than to an engine that cannot handle them.
Every non-prime field — QQ, GF(p^n), CC — takes the bar path under `"auto"`.

The result object records which path ran: an `HHTable` from the bar path is tagged
`normalized bar complex`, one from the engine is tagged `hanlab engine (F_p fast rank)`, so
you can always see after the fact which implementation produced a number.

## The cross-check philosophy

Two independently written implementations that always agree is a far stronger correctness
claim than one implementation that is merely tested against expected values. quiverlab
leans on this at several layers:

1. **engine ≡ bar equality batteries.** `HHTable` defines equality as "same kind and same
   dimension list" (`hochschild/table.py`), so `A.hochschild_cohomology(N, engine="fast")
   == A.hochschild_cohomology(N, engine="bar")` is a single, meaningful assertion. The
   engine test files run exactly these comparisons across primes and algebras.
2. **Backend ≡ oracle.** Every faster/smaller resolution — minimal A^e, Bardzell,
   Chouhy–Solotar — is required by contract to reproduce the bar oracle's numbers on the
   overlap range before it is trusted past the bar wall (Chapter 05). `BarResolution` is
   the designated ground truth.
3. **Kernels ≡ pure Python.** Within the engine, the numba-compiled kernels have
   pure-Python twins that are the permanent oracle; a parity test toggles
   `QUIVERLAB_NO_NUMBA` and checks both produce identical output (Chapter 05).

No single path is privileged as "the truth"; the truth is where the independent paths meet.
When they would disagree, the machinery is built to fail loudly (a `FieldError`, a
`DepthLimitError`, a raised `RuntimeError` inside `quotient_induced`) rather than quietly pick a
winner.

## A worked micro-example — dispatching k[x]/(x^3)

Build k[x]/(x^3) three ways and ask for HH^0..HH^3:

- **over GF(5)** (a `PrimeField`): `engine="auto"` produces an `HHTable` tagged
  `hanlab engine (F_p fast rank)`, while `engine="bar"` is tagged
  `normalized bar complex`; both give `[3, 2, 2, 2]` and compare equal — the cross-check in
  miniature.
- **over GF(4)** (a `FiniteField`): `engine="auto"` is tagged `normalized bar complex` —
  it fell back to the oracle because GF(4) is not a prime field. Forcing `engine="fast"`
  here raises `FieldError` ("the fast engine computes over prime fields only").

(All of these tags, dimensions, and the raised error were produced by running the code.)

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| the engine switch + dispatch predicate | `core/algebra.py` | `hochschild_cohomology`, `hochschild_homology`, `_use_fast_engine` |
| loud refusal of non-prime fields by the engine | `engine/adapter.py` | `to_engine` |
| pure bar path (any field) | `hochschild/bar.py` | `hochschild_cohomology_dims`, `hochschild_homology_dims` |
| engine path (F_p) | `engine/adapter.py`, `engine/scan3.py`, `engine/hh_engine.py` | `engine_cohomology_dims`, `engine_homology_dims` |
| equality used by the cross-checks | `hochschild/table.py` | `HHTable.__eq__` |
| numba-vs-pure fallback flag | `engine/_kernels.py` | `USE_KERNELS` (`QUIVERLAB_NO_NUMBA`) |
