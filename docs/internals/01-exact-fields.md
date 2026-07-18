# 01 — Exact fields

## The mathematics

Every computation in quiverlab happens over a field k. The library supports the
rationals QQ, the finite fields GF(p) and GF(p^n), and "the complex numbers" CC handled
*exactly* (computing in the algebraic subfield of C that your entries generate). A single
uniform interface — add, multiply, invert, test-for-zero — lets the rest of the code be
written once and run over any of them. Crucially, k is never approximated: there are no
floating-point numbers anywhere, because a rounding error would silently corrupt a rank
and hence a homology dimension.

## How it is represented

A field is a **class** (a blueprint bundling data with operations) called `Domain`. A
particular field is an *instance* of it — an object you can hand to the arithmetic. The
protocol is the set of methods every field must provide (`fields/domain.py`):

```
coerce(x)      # turn a raw user entry into an element of this field
zero(), one()  # the additive and multiplicative identities
add, neg, sub, mul, inv   # field arithmetic
is_zero(a), eq(a, b)      # exact equality tests
```

Each concrete field chooses how an *element* is stored — and the choices differ:

| field | constructor | an element is a... | example: how `x` prints |
|---|---|---|---|
| QQ (rationals) | `QQ` | Python `Fraction` (exact ratio of integers) | `Fraction(1, 3)` |
| GF(p), p prime | `GF(7)` | plain `int` in `0..p-1` | `5` |
| GF(p^n) | `GF(4)` | `tuple` of ints, little-endian coeffs mod p | `(0, 1)` = the generator x |
| CC (exact) | `CC` | a sympy algebraic-field element | `sqrt(2)` |

- **QQ** (`fields/rationals.py`) stores elements as Python's built-in `Fraction`. `add`
  is literally `a + b`, `inv` is `Fraction(1) / a`. That is all — Fractions are already
  exact.
- **GF(p)** (`fields/primefield.py`) stores elements as ordinary integers reduced modulo
  p. `mul(a, b)` is `(a * b) % p`; `inv(a)` uses `pow(a, -1, p)` (the modular inverse).
  The constructor rejects a non-prime p loudly with a `FieldError`.
- **GF(p^n)** (`fields/finitefield.py`) stores an element as a **tuple** of n integers —
  the coefficients of a polynomial over GF(p), little-endian (constant term first).
  Arithmetic is polynomial arithmetic *modulo a fixed monic irreducible polynomial*. That
  modulus comes from a bundled table (`fields/conway.py`, Conway-style polynomials from
  Lübeck's tables). The table is not trusted blindly: at construction the code runs
  `poly_is_irreducible` (trial division by every monic polynomial up to half the degree)
  and refuses a reducible modulus with a `FieldError`.
- **CC** (`fields/complexfield.py`) is special: `CC` itself is *not* a `Domain`. It first
  inspects *all* your entries, then hands back the concrete working field they generate.
  `make_domain(entries)` calls sympy's `construct_domain(..., extension=True)` to build
  the smallest exact algebraic extension of QQ containing your entries (e.g. QQ(sqrt(2)),
  or QQ(i)); elements then live in that sympy field. Genuinely transcendental input is
  refused loudly. The helper `E(n)` gives the exact primitive n-th root of unity
  exp(2*pi*i/n), following GAP's convention.

## How the computation runs — banning floats

Exactness is enforced at two independent layers.

1. **The runtime gate.** `reject_inexact(x)` (`fields/domain.py`) is called on the way in
   to every field. It raises `ExactnessError` if `x` is a `float` or `complex`, if a
   string carries a decimal point (`"0.5"`) or scientific notation (`"15e-1"` — which
   Python's `Fraction` would otherwise silently swallow), or even if `x` is a Python
   `bool` (because `True` is secretly the integer 1 and that is a footgun). Every message
   states the problem and a fix ("write '1/3' or Fraction(1, 3), never 0.333").
2. **The static gate (the "AST gate").** A test, `tests/test_no_floats.py`, parses every
   `.py` file under `src/quiverlab/` into its abstract syntax tree (the parsed form of the
   source) and walks it looking for any float or complex *literal* (like `0.5` or `1j`) or
   any call to `float(...)`. If it finds one, the whole test suite fails. So a float
   cannot even be *written* in the library, let alone reach a computation. The test also
   plants a known-bad file to confirm the detector still fires.

## Exact linear algebra over a Domain

Ranks, kernels, and solves are the workhorses of homology, and they are done exactly over
whatever `Domain` you are using (`fields/linalg.py`). `rref(rows, dom)` is textbook
Gauss–Jordan elimination, but every arithmetic step goes *through the domain*:
`dom.is_zero(...)` chooses pivots, `dom.inv(...)` normalises a pivot row, and
`dom.sub(x, dom.mul(f, y))` eliminates. Because the domain is exact, a pivot is either
exactly zero or exactly invertible — there is no tolerance, no "nearly singular". From
`rref` the module derives `rank` (count of pivot columns), `nullspace` (a basis of the
solution space of the homogeneous system, built by back-substituting each free column),
and `solve` (which returns `None`, i.e. "no solution", when the augmented column becomes a
pivot — the signal of inconsistency). This is the correctness-first path; over GF(p) a
much faster numeric path exists (Chapter 05), cross-checked against this one.

## A worked micro-example — GF(4)

`GF(4)` factors as 2^2, so quiverlab builds `FiniteField(2, 2)` with the bundled modulus
`[1, 1, 1]` = 1 + x + x^2 (irreducible over GF(2)). An element is a length-2 tuple
`(c0, c1)` meaning c0 + c1*x. Running the code:

- the generator `x` is `(0, 1)`; `one()` is `(1, 0)`; `zero()` is `(0, 0)`;
- `x * x` = `(1, 1)` — because x^2 = x + 1 after reducing modulo 1 + x + x^2;
- `x * x * x` = `(1, 0)` = 1 — confirming x has multiplicative order 3, as it must in the
  cyclic group of order 3 = 4 - 1;
- `coerce(3)` = `(1, 0)` — the integer 3 is reduced mod 2 to 1, then embedded as 1 + 0*x.

(These four outputs were produced by running the code.)

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| the field protocol | `fields/domain.py` | `Domain`, `reject_inexact`, `parse_rational` |
| rationals QQ | `fields/rationals.py` | `RationalField`, `QQ` |
| prime field GF(p) | `fields/primefield.py` | `PrimeField` |
| GF(p^n) polynomial arithmetic | `fields/finitefield.py` | `FiniteField`, `GF`, `poly_is_irreducible` |
| bundled irreducible moduli | `fields/conway.py` | `CONWAY` |
| exact CC | `fields/complexfield.py` | `ComplexField`, `CC`, `SympyExactDomain`, `E` |
| the exactness exception | `errors.py` | `ExactnessError` |
| the static float ban (AST gate) | `tests/test_no_floats.py` | `_violations`, `test_no_float_literals_or_calls_in_src` |
| exact rref / rank / nullspace / solve | `fields/linalg.py` | `rref`, `rank`, `nullspace`, `solve` |
