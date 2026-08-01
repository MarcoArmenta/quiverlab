# 04 — Hochschild (co)homology via the bar complex

## The mathematics

For A unital over k, the normalized bar complex computes Hochschild (co)homology. Using
the reduced space Abar = A/k.1, the chain space is C_n = A ⊗ Abar^{⊗n} with the Hochschild
boundary b, and the cochain space is C^n = Hom_k(Abar^{⊗n}, A) with the coboundary d. Then
dim HH_n = dim C_n − rank(b_n) − rank(b_{n+1}) and dim HH^n = dim C^n − rank(d^n) −
rank(d^{n−1}). This is the exact but *exponential* "oracle" backend: dim C_n = m·(m−1)^n
grows fast, so it is used on small algebras and over any field, and as the ground truth
against which faster engines are checked.

## How it is represented

Work in a unit-adapted algebra (Chapter 03), so 1_A = b_0 and Abar = span of basis vectors
1, ..., m−1. A **cochain basis element** is a `tuple` `(s, J)` (`hochschild/bar.py`,
`_cochain_basis`):

- `s` is an integer in `0..m−1` — the *output* slot, the basis vector b_s that this
  cochain returns;
- `J` is a tuple of n integers, each in `1..m−1` — the *input* reduced tensor
  b_{J[0]} ⊗ ... ⊗ b_{J[n−1]}.

So `(s, J)` is the elementary cochain that sends the reduced tensor J to b_s and every
other reduced tensor to 0. For a homology chain the same tuple shape `(s, J)` is reused,
now meaning b_s ⊗ b_{J[0]} ⊗ ... (`boundary_matrix`). Verified by execution: for k[x]/(x^2)
(m = 2, reduced basis just {x} = index 1),

```
C^0 basis = [(0, ()), (1, ())]
C^1 basis = [(0, (1,)), (1, (1,))]
C^2 basis = [(0, (1,1)), (1, (1,1))]
```

The coboundary is stored as a plain matrix `D` — a list of lists of field elements — with
one **row per C^{n+1} basis element** and one **column per C^n basis element**. A dict
`row_index` maps each output cochain to its row.

## How the computation runs

`coboundary_matrix(A, n, max_cells)` fills `D` one contribution at a time. Reading the
normalized coboundary
(d f)(a_1,...,a_{n+1}) = a_1·f(a_2,...,a_{n+1}) + Σ_{i=1}^{n} (−1)^i f(...,a_i a_{i+1},...)
+ (−1)^{n+1} f(a_1,...,a_n)·a_{n+1}, the code loops over each input cochain column `(s, J)`
and each candidate output tensor `K` (an (n+1)-fold reduced tuple), and adds the three
kinds of term into `D` via a helper `bump(t, K, ci, val)` that writes `val` into row
`(t, K)`, column `ci`:

1. **Front term** `b_{K0}·f(K[1:])`: fires when `K[1:] == J`. It multiplies b_{K[0]} into
   the output using the table row `A.T[K[0]][s]`, and each nonzero coordinate `t` bumps
   row `(t, K)`.
2. **Interior terms** i = 1..n: fire when K matches J everywhere except that positions
   i−1, i of K contract to J's position i−1. The coefficient is the *Abar component* of
   b_{K[i-1]}·b_{K[i]} — read straight from the table `A.T[K[i-1]][K[i]][J[i-1]]` — with
   sign (−1)^i (the unit coordinate is excluded, i.e. the projection to Abar).
3. **Back term** `(−1)^{n+1} f(K[:n])·b_{K[n]}`: fires when `K[:n] == J`, multiplying by
   b_{K[n]} on the right via `A.T[s][K[n]]`, with the alternating sign.

`hochschild_cohomology_dims(A, top, max_cells)` first calls `unit_adapted()`, then for each
n takes `rank(D, dom)` over the domain (the exact rref of Chapter 01), and finally does the
Euler-characteristic bookkeeping: `dim HH^n = m·(m−1)^n − rank(d^n) − rank(d^{n−1})`. The
homology twin (`boundary_matrix`, `hochschild_homology_dims`) is identical in spirit with b
in place of d. Results are returned as an `HHTable` (`hochschild/table.py`) — a small
object that prints as `HH^0 = ...  HH^1 = ...` and compares equal to another table with the
same kind and dimensions (this equality is what the cross-check batteries of Chapter 07
use).

### The `max_cells` guard

Because dim C_n is exponential, a differential matrix can be astronomically large.
`_check_cells(rows, cols, max_cells, what)` computes `rows * cols` *before* allocating and
raises `DepthLimitError` if it exceeds `max_cells` (default 4,000,000). The message states
exactly which differential and how big it would have been, and the hint explains that the
bar oracle is exponential and deeper engines (Bardzell, minimal, Chouhy–Solotar) are the
right tool past this wall. So `DepthLimitError` is never a crash — it is the code telling
you the certified range ran out and pointing at the faster path.

## A worked micro-example — one coboundary entry for k[x]/(x^3)

m = 3, unit-adapted basis {1, x, x^2}, reduced indices {1, 2} = {x, x^2}. Consider the
1-cochain f = `(0, (1,))`: it sends the reduced input x to b_0 = 1, and x^2 to 0. We compute
the column of d^1 for this f, at the output tensor K = (x, x) i.e. `(1, 1)`.

By the formula, (d^1 f)(x ⊗ x) = x·f(x) − f(x·x) + f(x)·x. Now f(x) = 1 and f(x^2) = 0, so
this is x·1 − f(x^2) + 1·x = x − 0 + x = 2x = 2·b_1. In the code:

- the **front term** fires (`K[1:] = (1,) = J`): b_1·f(x) uses `A.T[1][0] = [0,1,0]`,
  bumping +1 into row `(1, (1,1))`;
- the **interior term** i = 1 needs the Abar component of x·x = x^2 at slot J[0] = x
  (index 1): `A.T[1][1][1] = 0`, so it contributes nothing (indeed f(x^2) = 0);
- the **back term** fires (`K[:1] = (1,) = J`): f(x)·b_1 uses `A.T[0][1] = [0,1,0]` with
  sign (−1)^{2} = +1, bumping +1 into row `(1, (1,1))`.

The two +1's add to the entry `D[row (1,(1,1))][col (0,(1,))] = 2` — exactly the coefficient
of x in 2x. (This matrix was printed by running the code. Restricted to the output tensor
(x,x), the column has this single nonzero 2; the full column for `(0,(1,))` also carries
`D[(2,(1,2))] = 1` and `D[(2,(2,1))] = 1` from the other output tensors.)

## The product surface (Plan 35) — cup, cap, bracket, Connes B

On top of the (co)homology *dimensions*, the Tamarkin–Tsygan calculus is exposed as
**structure-constant tables**: for each pair of degrees, the product of a basis class
on the left with a basis class on the right, written back in the basis of the target
degree. Four public methods on `Algebra` return frozen result objects:

| method | operation | returns |
|---|---|---|
| `cup_products(top, engine="auto", max_cells=4_000_000)` | `HH^p ⊗ HH^q → HH^{p+q}`, all `p+q ≤ top` | `HHProducts` |
| `cap_products(top, engine="auto", max_cells=4_000_000)` | `HH^p ⊗ HH_n → HH_{n−p}`, `p ≤ n ≤ top` | `HHProducts` |
| `gerstenhaber_brackets(top, engine="auto", max_cells=4_000_000)` | `HH^p ⊗ HH^q → HH^{p+q−1}`, `p,q ≥ 1` | `HHProducts` |
| `connes_differentials(top, max_cells=4_000_000)` | induced `B : HH_n → HH_{n+1}`, `0 ≤ n < top` | `ConnesB` |

`HHProducts.tables` is a `dict` keyed by the degree pair `(p, q)` (for cap, `(p, n)`)
whose values are `ProductTable`s; a `ProductTable` carries `dims = (dim_left,
dim_right, dim_out)` and `constants[k][i][j]` — the coefficient of the `k`-th output
class in `left_i · right_j`, as an **exact string** (int mod p on the GF(p) route,
a Domain repr on the CS route; the float ban forbids anything else). `ConnesB` carries
`hh_dims`, the per-degree matrices (rows indexed by `HH_{n+1}`, columns by `HH_n`),
and their `ranks`.

**Engine routing** (`core/algebra.py::_product_dispatch`). `engine="auto"` uses the
GF(p) bar/tt-calculus route over a prime field (`hochschild/products.py::
gfp_product_tables` → `engine/tt_calculus.py`), and for a *presented* algebra over any
other exact Domain the **Chouhy–Solotar native diagonal** route
(`resolutions_cs/products.py`); over GF(p) it also falls back to CS when the bar route
hits `DepthLimitError` (except for the bracket, which has no CS route). `engine="cs"`
forces the CS route; `engine="bar"` forces the GF(p) tt facade (loud off GF(p)). The
**Gerstenhaber bracket is GF(p)-only and window-bounded** — the result records the
served `window`, and the degree-0 insertion action is out of scope. Connes `B` needs
no engine choice: GF(p) via the fast `(b,B)` engine, any other exact Domain via the
generic bar `(b,B)` mixed complex.

**Basis-provenance caveat.** The structure constants are **basis-dependent** — they
are read on whatever HH basis the route produced, and that basis is *not* canonical
across routes. Every product object therefore **records which basis it used**:
`HHProducts.basis` is a provenance string (`"bar/GF(p)"` or `"cs/…"`) and `ConnesB`
carries its `engine` string. Do not compare raw `constants` across engines; the
cross-engine gate compares only basis-independent data (dims and the flattened rank of
the constants tensor mod p). Within one object the constants are exact and the
Gerstenhaber-algebra identities (graded-commutative associative cup, antisymmetric
cup-Leibniz bracket, cap module law, `B²=0`) hold on the nose — see the identity
batteries in `tests/hochschild/test_products_identities.py`.

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| cochain / chain basis element `(s, J)` | `hochschild/bar.py` | `_cochain_basis`, `_abar_tuples` |
| public product methods | `core/algebra.py` | `cup_products`, `cap_products`, `gerstenhaber_brackets`, `connes_differentials` |
| product result objects | `hochschild/products.py` | `HHProducts`, `ProductTable`, `ConnesB` |
| GF(p) bar/tt product tables | `hochschild/products.py` | `gfp_product_tables` |
| CS-native (any Domain) product tables | `resolutions_cs/products.py` | `cs_product_tables` |
| induced Connes `B` | `hochschild/products.py` | `connes_b_tables` |
| coboundary matrix, entry by entry | `hochschild/bar.py` | `coboundary_matrix` (the `bump` helper) |
| boundary matrix (homology) | `hochschild/bar.py` | `boundary_matrix` |
| dims from ranks | `hochschild/bar.py` | `hochschild_cohomology_dims`, `hochschild_homology_dims` |
| the size guard | `hochschild/bar.py` | `_check_cells` |
| the result table | `hochschild/table.py` | `HHTable` |
| depth-limit exception | `errors.py` | `DepthLimitError` |
