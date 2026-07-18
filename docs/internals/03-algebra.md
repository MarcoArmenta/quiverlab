# 03 — The algebra

## The mathematics

A finite-dimensional associative unital algebra A over k is fixed, up to isomorphism, by a
basis b_0, ..., b_{m-1} together with its structure constants: the coefficients expressing
each product b_i * b_j back in the basis, plus the coordinates of the unit 1_A. Every
downstream computation — Hochschild homology, Cartan matrices, the bar complex — needs
only these numbers. So quiverlab's central object is not a quiver and not a presentation;
it is the *structure-constant algebra*, and everything else is a way of producing one.

## How it is represented

The **class** `Algebra` (`core/algebra.py`) carries:

- `domain` — the field k (a `Domain` from Chapter 01);
- `dim` — the dimension m;
- `T` — the multiplication table, described below;
- `unit` — a **list** of m field elements, the coordinates of 1_A;
- `basis_labels` — an optional list of human-readable strings, one per basis vector;
- `quiver`, `relations` — kept when the algebra came from a presentation (used by the
  Cartan machinery), otherwise `None`.

**What `T[i][j]` literally is.** `T` is a list of lists of lists. `T[i][j]` is *the
coordinate vector of the product b_i * b_j* — a list of m field elements, the coefficients
c such that b_i * b_j = sum_t c[t] * b_t. It is a **list of coefficients, not a dict**, and
it always has length m (mostly zeros). So `T[i][j][t]` is the coefficient of b_t in
b_i * b_j. (Verified by reading `build_monomial_algebra` and by executing the example
below.)

`unit` is *unit-adapted* when b_0 = 1_A exactly — that is, `unit[0]` equals `one()` and
every other coordinate is zero. The flag `is_unit_adapted` records this; the constructor
computes it if you do not supply it. The bar complex (Chapter 04) *requires* this form,
because it identifies the reduced space A/k.1 with "drop coordinate 0". An algebra not in
this form is silently fixed by `unit_adapted()` before any Hochschild computation.

## How the computation runs

### Multiplying two elements

`multiply(u, v)` takes two coordinate vectors and returns their product's coordinate
vector, purely from `T`: for every nonzero u[i] and nonzero v[j] it forms the scalar
`c = u[i]*v[j]` and adds `c * T[i][j][t]` into output coordinate t. It skips zero
coordinates for speed but is otherwise the literal bilinear extension of the table.

### Building the table from kQ/I

`build_monomial_algebra(quiver, relations, field)` (`core/monomial.py`) turns a certified
monomial presentation into an `Algebra`:

1. **Admissibility.** Each relation must have all paths of length >= 2 (the ideal sits
   inside the square of the arrow ideal); a length-1 path raises `AdmissibilityError`.
2. **Field placement.** The field is materialised on `{0, 1}` via `make_domain` (Chapter
   01), and every relation coefficient is checked nonzero in that field — a relation whose
   coefficient vanishes mod p is refused as "0 = 0".
3. **Basis.** `irreducible_paths` (Chapter 02) returns the finite list of irreducible
   words. The basis is the trivial paths `("e", v)` for each vertex, followed by the path
   words `("p", w)`; `index` is a dict mapping each basis element (those `("e", v)` / `("p", w)` tuples) to its position.
4. **The product rule.** A local function `prod(x, y)` implements multiplication of two
   basis elements: `None` (i.e. zero) if their endpoints do not meet; the other factor if
   one is a trivial path; otherwise the concatenated word — unless that word now contains a
   forbidden subword, in which case it is again `None`.
5. **Filling `T`.** For every ordered pair `(i, j)` the code computes `prod(b_i, b_j)`; if
   it is a basis element, `T[i][j]` is the indicator vector with `one()` in that slot,
   else the zero vector.
6. **The unit.** `unit` is the sum of the trivial idempotents e_v — so for a one-vertex
   quiver, `unit` already has `one()` in slot 0 and the algebra is born unit-adapted.

### Why the code insists on unit-adaptation

`from_structure_constants` (the hand-built entry point) runs `_validate`: it checks the
supplied unit really is a two-sided identity on every basis vector, and that the table is
associative — `(b_i b_j) b_k = b_i (b_j b_k)` for all i, j, k — raising `QuiverlabError`
otherwise. quiverlab never guesses a structure constant. To *reach* unit-adaptation when
1_A is spread across several idempotents (a multi-vertex algebra), `unit_adapted()`
performs a change of basis: it builds an invertible matrix P whose column 0 is the unit
vector, then `change_of_basis(P)` recomputes `T` and `unit` in the new basis by solving
the linear systems `P x = (b_i b_j)` exactly (via `solve` over the domain). The result is
an isomorphic algebra whose basis vector 0 is 1_A.

## A worked micro-example — the literal `T` for k[x]/(x^3)

Built over QQ, `Q.algebra(relations=["x^3"], field=QQ)` gives dim 3, basis labels
`["e_1", "x", "x*x"]`, so b_0 = 1, b_1 = x, b_2 = x^2, and `unit = [1, 0, 0]`
(already unit-adapted). The actual table (Fractions written as plain integers):

```
T[0] = [ [1,0,0], [0,1,0], [0,0,1] ]      # 1*1=1,   1*x=x,   1*x^2=x^2
T[1] = [ [0,1,0], [0,0,1], [0,0,0] ]      # x*1=x,   x*x=x^2, x*x^2=0
T[2] = [ [0,0,1], [0,0,0], [0,0,0] ]      # x^2*1=x^2, x^2*x=0, x^2*x^2=0
```

Read one cell: `T[1][1] = [0, 0, 1]` says x * x = 0*1 + 0*x + 1*x^2 = x^2. And
`T[1][2] = [0, 0, 0]` says x * x^2 = 0 — because the concatenated word `xxx` contains the
forbidden relation and collapses to zero. (This table was produced by running the code.)

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| the structure-constant algebra | `core/algebra.py` | `Algebra` |
| the meaning of `T[i][j]` | `core/algebra.py` | `Algebra.__init__`, `multiply` |
| hand-built table + validation | `core/algebra.py` | `from_structure_constants`, `_validate` |
| unit-adaptation and base change | `core/algebra.py` | `unit_adapted`, `change_of_basis` |
| kQ/I -> table | `core/monomial.py` | `build_monomial_algebra` |
| the associativity/unit error | `errors.py` | `QuiverlabError`, `AdmissibilityError` |
