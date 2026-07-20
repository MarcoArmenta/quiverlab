# 06 — Invariants: Cartan, Coxeter, and the GF(p) extras

## The mathematics

For a bound quiver algebra A = kQ/I the Cartan matrix C records
C_{ij} = dim e_i A e_j — the number of basis paths from vertex i to vertex j — and is an
integer matrix independent of the ground field. When C is invertible over Q the Coxeter
transformation is Phi = −C^{−T} C, and its characteristic polynomial (the Coxeter
polynomial) is a classical invariant tying the algebra to Dynkin/Euclidean type. A second
family of invariants lives over the enveloping algebra and is only computed via the fast
GF(p) engine: whether A is Frobenius or symmetric, and, when it is, its Nakayama
automorphism — the twist by which the dualizing bimodule differs from A itself.

## How it is represented

The Cartan matrix is a plain **list of lists of ints** (an integer matrix); the Coxeter
matrix is either a list of lists of ints or, when it must, a list of lists of exact sympy
rationals; the Coxeter polynomial is a sympy `Poly` in the symbol `t` (an exact polynomial
object, not a float approximation). The Nakayama automorphism is returned as a list-of-lists
**integer matrix** whose columns are the images of the basis vectors, in the unit-adapted
basis, with entries taken mod p.

## How the computation runs

### Cartan from quiver provenance

`cartan_matrix(A)` (`invariants/cartan.py`) does not re-derive anything homological — it
literally counts the path basis that Chapter 03 already produced. It requires the algebra
to remember its presentation: if `A.quiver` or `A.basis_labels` is `None` (a hand-built
structure-constant algebra carries no path basis) it raises `QuiverlabError` with a hint to
build via `Quiver.algebra`. Otherwise it walks `A.basis_labels`: a label `"e_v"` is a
trivial path and bumps the diagonal `C[v][v]`; any other label is a path word, and bumps
`C[source][target]` by reading the word's endpoints off the quiver. The result is exact and
field-independent.

### Coxeter matrix and polynomial, and when they fail loudly

`coxeter_matrix(A)` lifts C into a sympy `Matrix` and first checks `C.det()`. If the
determinant is **zero**, Phi = −C^{−T} C is undefined and the code raises `QuiverlabError`
("Cartan matrix is singular"), with the hint that this happens e.g. at infinite global
dimension when |det C| ≠ 1. Otherwise it forms `Phi = -C.inv().T * C` exactly over Q.
Here is the **det ≠ ±1 caveat**: classically the Coxeter matrix is an *integer* matrix, and
that holds exactly when det C = ±1 (so C is unimodular). When det C is not ±1, Phi can have
genuine rational entries; rather than round or lie, the code detects non-integer entries
(`x.q != 1`) and returns the *exact rationals* (via `sympy.nsimplify`), falling back to ints
only when every entry really is integral. `coxeter_polynomial(A)` likewise refuses a
singular C, then returns `Phi.charpoly(t)` as an exact `Poly` — no numerical root-finding
(the exact spectral-radius and Mahler-measure layer lands with Plan 05; see "The exact
spectral layer" below). As of Plan 05 this is documented on `coxeter_polynomial` itself:
the domain follows the coefficients (ZZ when integral, QQ when genuinely rational), so a
non-unimodular Cartan that yields a rational Coxeter transformation is called out, not
silent (see below).

### The engine-backed GF(p) extras

`nakayama_automorphism`, `is_frobenius`, `is_symmetric`, and `cyclic_homology` are methods
on `Algebra` that route through the fast engine. Each first calls `_require_prime_field`,
which raises `FieldError` unless the domain is a `PrimeField` — because these invariants are
wired only through the numpy-int64/mod-p engine (`engine/coxeter.py`), which does its linear
algebra over F_p with modular inverses. Over QQ or GF(p^n) they are not yet available and
say so.

- `is_frobenius` searches for a non-degenerate Frobenius form: `frobenius_form` tries a
  deterministic sequence of covectors λ (coordinate functionals first, then all-ones, then
  seeded pseudo-random ones), forms the Gram matrix G_{ij} = λ(e_i e_j), and accepts the
  first λ whose G is full-rank mod p. Frobenius ⇔ self-injective for a finite-dimensional
  algebra, so a `None` result means "not Frobenius".
- `nakayama_automorphism` takes that Frobenius form and returns N = G^{−1} G^T as an integer
  matrix (columns = images); it raises `ValueError` if the algebra is not Frobenius. N is
  the identity exactly when G is symmetric — i.e. when A is *symmetric*.
- `is_symmetric` is "Frobenius **and** the Nakayama automorphism is the identity mod p".

`complexity(A, n)` (`invariants/scalar.py`) is engine-backed the same way — it reads the
minimal A^e resolution's term-dimension growth through `complexity_of`, so it too calls
`_require_prime_field` and is GF(p)-only — but the number it returns can silently
**under-report**, exact only on local / single-vertex inputs. Two unguarded caveats travel
with it: the resolution's radical logic is **local-only**, so a multi-vertex algebra of
genuinely *infinite* complexity can still come back *finite* (read it as a lower bound), and
a memory-truncated build adds a silent prefix because the truncation marker is not consulted.
No vertex-count guard is added on purpose — it would break the correct
`linear_path_algebra(2) → 0`.

## A worked micro-example — A_2 and k[x]/(x^2)

For `linear_path_algebra(2)` (vertices 1, 2, arrow a: 1 → 2), the basis labels are
`["e_1", "e_2", "a1"]`. Counting: `e_1` bumps `C[0][0]`, `e_2` bumps `C[1][1]`, and the
path `a1` (source 1, target 2) bumps `C[0][1]`. So
`C = [[1, 1], [0, 1]]`. Then Phi = −C^{−T} C = `[[-1, -1], [1, 0]]` (integer, since
det C = 1), and the Coxeter polynomial is `t**2 + t + 1` — the A_2 Coxeter polynomial.

Over GF(5): `k[x]/(x^2)` reports `is_frobenius = True`, `is_symmetric = True`, and its
Nakayama matrix is the 2×2 identity `[[1, 0], [0, 1]]` — as it must be for a symmetric
algebra. Asking `A_2` (not self-injective) for its Nakayama automorphism raises
`ValueError` ("not Frobenius"). (All five outputs above were produced by running the code.)

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| Cartan from the path basis | `invariants/cartan.py` | `cartan_matrix` |
| Coxeter matrix −C^{−T}C, singular guard | `invariants/cartan.py` | `coxeter_matrix` |
| exact Coxeter polynomial | `invariants/cartan.py` | `coxeter_polynomial` |
| public method wrappers + GF(p) gate | `core/algebra.py` | `cartan_matrix`, `coxeter_matrix`, `_require_prime_field`, `nakayama_automorphism`, `is_frobenius`, `is_symmetric` |
| Frobenius form / Nakayama over F_p | `engine/coxeter.py` | `frobenius_form`, `is_frobenius`, `nakayama_automorphism` |
| cyclic homology (GF(p)) | `core/algebra.py`, `engine/cyclic.py` | `cyclic_homology`, `cyclic_homology_dims` |
| the field-gate exception | `errors.py` | `FieldError`, `QuiverlabError` |

## The exact spectral layer (Plan 05)

`invariants/spectral.py` computes the **spectral radius** and **Mahler measure** of the
Coxeter polynomial exactly, reimplementing the hanlab float layer (which used mpmath
`nroots`) with exact sympy algebraic numbers. `spectral_radius(p)` is `max_i |α_i|` and
`mahler_measure(p)` is `|lc|·∏_{|α|>1}|α|`; both short-circuit to the exact integer 1 when
`is_cyclotomic_product(p)` (all roots on the unit circle). No floats are used — magnitudes
are `sympy.Abs` of `CRootOf` roots and comparisons are `.is_positive`. The subtle part is
**soundness for complex roots**: `real_roots` alone is unsound (real-roots-suffice is a
theorem only for hereditary quivers), so the code forms the non-cyclotomic part `q` and,
via the self-inversive `y = z + 1/z` substitution, uses a Sturm real-root count of `Q(y)` —
with no complex-root isolation — to decide between fast `real_roots(q)` (Branch A:
hereditary/Salem/Lehmer) and correct `all_roots(q)` (Branch B: non-hereditary,
complex-dominant). The Lehmer star T(2,3,7) = `star_quiver([1,2,6])` carries Lehmer's
polynomial (Branch A), whose spectral radius is Lehmer's number 1.17628…, the smallest known
Mahler measure > 1; the trivial extension T(A) = A ⋉ DA collapses the Coxeter polynomial to
(t+1)^v regardless of representation type (Cartan C ↦ C + Cᵀ, Φ = −I).

The **non-unimodular caveat** (§ "Coxeter matrix and polynomial") is now **documented**:
`coxeter_polynomial`'s docstring names the `det C ∉ {0, ±1}` case, where Φ may be rational
so the Coxeter polynomial is over **QQ** (e.g. t² + 3t/2 + 1 for C = [[2,1],[0,1]]) — exact,
but not the classical integral Coxeter transformation. The domain follows the actual
**coefficients** (sympy's inference), not det C: a non-unimodular Cartan can still be integral
(k[x]/(x²) → t+1 over ZZ; diag(1,2) → (t+1)² over ZZ). The `coxeter_matrix` sibling surfaces
the same fact via its rational-entry branch.
