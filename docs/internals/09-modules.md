# 09 — Modules: representation and resolutions

## The mathematics

For a bound quiver algebra A = kQ/I, a finite-dimensional **right** A-module M is a
vector space carrying a right action of A. Every such module is glued from **simples**
S_v (one per vertex), and each simple sits atop an **indecomposable projective**
P_v = e_v A (the paths starting at v) and beneath an **indecomposable injective**
I_v = D(A e_v) (the dual of the paths ending at v). The radical filtration
M ⊇ rad M ⊇ rad² M ⊇ ... records the Loewy structure; the **top** M/rad M lists the
generators and the **socle** the essential simple submodules. Homological questions —
Hom_A(M, N), Ext^n_A(M, N), the projective dimension of M — are answered from a
**minimal projective resolution** ... → Q_1 → Q_0 → M → 0.

## How it is represented

A module is the class `Module` (`modules/module.py`). It stores its k-dimension `dim`
and a **dict** `action` mapping each algebra basis label (a **string**, e.g. `"e_1"`,
`"a"`, `"a*b"`) to a `dim × dim` **matrix** (a list of lists of exact Domain elements)
of **right** multiplication. The convention is stated in every docstring: an element m is
a **column vector** in a fixed k-basis of M, and m·b = `action[b] @ m`; because paths
compose left-to-right, the action is an *anti-homomorphism*,
`action["x*y"] = action["y"] @ action["x"]`. The vertex subspace M·e_v is the image of
`action["e_v"]`, so the **dimension vector** — a dict vertex → int — is just
`rank(action["e_v"])` per vertex, and these ranks sum to `dim`. A simple S_v has `dim 1`,
`action["e_v"] = [[1]]` and every other label acting as `[[0]]`. A projective P_v is read
directly off the algebra's own multiplication table: its basis is the sub-list of
`A.basis_labels` whose path starts at v, and `action[b]` is right multiplication by b
restricted to that sub-basis (column k = coordinates of pathₖ·b). An injective I_v is the
**transpose** of the left-multiplication action on the paths ending at v (the k-dual).

All matrix work goes through `modules/linalg_mod.py`, which is a thin exact layer over
`fields.linalg` (`rank`, `nullspace`, `rref`, `solve`): matrix product, Kronecker
product, kernel, column-space pivots, "express these columns in that basis", and the
greedy "independent-modulo-a-subspace" selector that picks minimal generators. Nothing
here is field-specific, so a module lives over ℚ, ℚ(α), GF(p), or GF(p^n) identically —
the same code, different Domain.

## How the computation runs

### Radical, top, socle

`rad M` (`modules/radtopsoc.py`) is the sum of the images of the arrow-action matrices
(`M·rad A = Σ_α image(action[α])`), reduced to an independent column basis; `top M` is
the quotient `M/rad M`; `soc M` is the intersection of the arrow-action kernels
(`∩_α ker(action[α])`). A submodule (radical, socle) is rebuilt as its own `Module` by
**restricting** the action — for each generator b, apply `action[b]` to each submodule
basis column and solve for the coordinates back in that basis (`solve_columns`); a
quotient (top) is rebuilt by choosing coset representatives independent modulo the
submodule and reading the action modulo it. Iterating `rad` gives the radical series, and
its length is the Loewy length.

### Hom and Ext

`Hom_A(M, N)` (`modules/hom.py`) is the space of matrices φ (dim N × dim M) commuting
with the action: `N.action[b] @ φ = φ @ M.action[b]` for every generator b (arrows and
idempotents suffice). Column-stacking turns this into one homogeneous linear system
`(I ⊗ N.action[b] − M.action[b]ᵀ ⊗ I) vec(φ) = 0`; the nullspace is Hom, its dimension is
`A.hom(M, N)`. `Ext^n` (`modules/ext.py`) applies Hom_A(−, N) to a minimal resolution of
M and takes cohomology: with the Hom-space bases H_n, the coboundary δ^n is precomposition
with d_{n+1}, and dim Ext^n = dim H_n − rank δ^n − rank δ^{n-1}.

### The minimal projective resolution

`modules/resolution.py` (generalized from the bridge `obstruction/module_ext.py`, which
was hardcoded to a 4-vertex diamond over ℚ; here MIT-headered, over any vertex set and any
Domain) builds ... → Q_1 → Q_0 → M → 0 by **iterated projective covers**
(Green–Solberg–Zacharia). One step:

1. **Top generators.** Compute `top M = M/rad M` and choose a k-basis, each vector lifted
   into a single vertex block M·e_v (homogeneous generators). If vertex v supplies t_v
   generators, the cover is Q_0 = ⊕_v P_v^{t_v}.
2. **Cover map.** d_0 : Q_0 → M sends the canonical generator of each P_v summand (its
   e_v basis vector) to the lifted top-generator g; a general P_v basis vector for path p
   maps to g·p = `action_M[p] @ g`.
3. **Syzygy.** Ω_1 = ker(d_0), a submodule of Q_0, rebuilt by restriction; cover it to get
   Q_1 and d_1 = (cover of Ω_1) followed by (inclusion Ω_1 ↪ Q_0). Repeat.

Minimality is guaranteed because generators are chosen independent **modulo the radical**
(the `independent_modulo` selector), so d_n(Q_n) ⊆ rad Q_{n-1} and the summand counts are
the true Betti numbers. The resolution **terminates** at length n exactly when Ω_n = 0
(projective dimension n); a term dim blow-up past `max_term_dim` raises `DepthLimitError`
with the certified length. `M.projective_resolution(k)` returns a `ProjectiveResolution`:
`.term(n)` (the summand vertices), `.betti(n)`, `.differential(n)`, `.pd()`, and a readable
`P_1 <- P_2 <- 0` repr. `global_dimension(A) = sup_v pd(S_v)` — exact when every simple
resolves within the depth budget, else a labeled certified lower bound.

## A worked micro-example — S_1 over the A₂ path algebra

`linear_path_algebra(2)` is Q: 1 → 2 with basis `["e_1", "e_2", "a1"]` (the arrow is
auto-named `a1`). The simple S_1 has `action["e_1"] = [[1]]`,
`action["e_2"] = action["a1"] = [[0]]`. Its cover: top S_1 = S_1, one generator at vertex
1, so Q_0 = P_1 = e_1 A (basis `[e_1, a1]`, dimvec {1:1, 2:1}) and d_0 : P_1 → S_1 kills a1.
The syzygy Ω_1 = ker(d_0) = span{a1} ≅ S_2 = P_2 (basis `[e_2]`,
dimvec {1:0, 2:1}), which is projective, so Q_1 = P_2 and Ω_2 = 0. The resolution is
0 → P_2 → P_1 → S_1 → 0, `pd(S_1) = 1`, and since pd(S_2) = 0 the algebra is hereditary,
`global_dimension = 1`. Then Hom(S_1, S_2) = 0 and Ext^1(S_1, S_2) = 1 (the single arrow
1 → 2). (These values were produced by running the code.)

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| the module object, dimension vector | `modules/module.py` | `Module`, `dimension_vector`, `from_arrow_action` |
| exact matrix layer over a Domain | `modules/linalg_mod.py` | `matmul`, `kron`, `kernel_columns`, `independent_modulo`, `solve_columns` |
| simples / projectives / injectives | `modules/builders.py` | `simple`, `projective`, `injective` |
| radical / top / socle, submodule/quotient | `modules/radtopsoc.py` | `radical`, `top`, `socle`, `submodule`, `quotient` |
| Hom / End | `modules/hom.py` | `hom_space`, `hom_dim`, `end_dim` |
| minimal projective resolution | `modules/resolution.py` | `minimal_resolution`, `projective_cover`, `ProjectiveResolution` |
| module Ext^n, global dimension | `modules/ext.py` | `ext`, `ext_dims`, `global_dimension` |
| public methods on the algebra | `core/algebra.py` | `simple`, `projective`, `injective`, `hom`, `ext`, `global_dimension` |
