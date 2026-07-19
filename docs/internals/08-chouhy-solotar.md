# 08 — The Chouhy–Solotar resolution

## What this computes
Given `A = kQ/I` (any admissible presentation, any exact field), the CS engine builds a
small projective resolution of `A` as an `A`-bimodule and reads Hochschild (co)homology
off it — reaching degrees the bar complex never can.

## The objects (how they are represented in code)
- **Reduction system `rs`** (`quiverlab.groebner.build_reduction_system`): the confluent
  rewriting rules `s → f_s`. `rs.leading_words()` are the tips `S`; `rs.irreducibles` a
  basis of `A`. Words are tuples of arrow NAMES, left to right.
- **Ambiguity chain `Chain`** (`resolutions_cs.terms`): an element of `S_n`. `word` is the
  path; `blocks = (u_0,…,u_{n-1})` its unique CS decomposition; `o,t` its endpoints.
  `S_0`=vertices, `S_1`=arrows, `S_2`=tips, `S_n (n≥3)` = `(n-1)`-fold overlaps of tips.
  Computed by reusing Bardzell's `associated_paths` on the tip monomial algebra `A_S`.
- **Resolution term `P_n = ⨁_{σ∈S_n} A e_{o(σ)} ⊗ e_{t(σ)} A`.** Tensoring/homming down,
  each `σ` contributes the corner `e_tAe_o` (homology) or `e_oAe_t` (cohomology).
- **A term** `(coeff, a, τ, c)` means `coeff · (a ⊗ τ ⊗ c)`; `a,c` are paths, `τ ∈ S_{n-1}`.

## The computation, step by step (on k[x]/(x²))
1. Tips `S = {xx}`. `S_n = {x^n}`; `P_n ≅ A^e`; collapsed `C_n = A = ⟨1,x⟩`.
2. Leading map: `d_n` odd = `x⊗1 − 1⊗x`; even = `1⊗x + x⊗1`. No correction (monomial).
3. Collapse to `C_n`: odd `↦ 0`; even `↦` multiply-by-`2x` = `[[0,0],[2,0]]`, rank 1.
4. `HH_n = dim C_n − rank d_n − rank d_{n+1} = [2,1,1,1,…]` (char 0); `[2,2,…]` (char 2).

## The differential in general (the one subtle step)
The leading map `δ_n` is Bardzell's (depends only on the tips). The tails add a correction
strictly below `σ` in the reduction order; its coefficients are the solution of the linear
system `d_{n-1}∘d_n = 0` (CS Theorems 4.1/4.2 — the same trick CS use in §6). Two gates
certify every run: `assert_dd_zero` and `assert_order_condition`.

## Worked non-monomial example — the commutative square (HH^• = (1,0,0))
`A = kQ/(ab − cd)`, `Q`: `1→2→4`, `1→3→4` (arrows `a,b,c,d`), tip `cd`, `dim A = 9`.
`S_0 = {e_1,e_2,e_3,e_4}`, `S_1 = {a,b,c,d}`, `S_2 = {cd}`, `S_n = ∅ (n≥3)`.
Cohomology terms: `C^0 = ⟨ê_1..ê_4⟩` (dim 4), `C^1 = ⟨â,b̂,ĉ,d̂⟩` (dim 4), `C^2 = ⟨âb⟩` (dim 1).

`δ^0` sends `(λ_1,λ_2,λ_3,λ_4) ↦ ((λ_2−λ_1)a, (λ_4−λ_2)b, (λ_3−λ_1)c, (λ_4−λ_3)d)`:

```
        e1  e2  e3  e4
   a  [ -1   1   0   0 ]
   b  [  0  -1   0   1 ]        rank 3,  ker = <(1,1,1,1)> = the centre  ⟹  HH^0 = 1
   c  [ -1   0   1   0 ]
   d  [  0   0  -1   1 ]
```

`δ^1` sends a 1-cochain `g(a)=αa, g(b)=βb, g(c)=γc, g(d)=δd` to `(δ+γ−β−α)·ab`:

```
        a   b   c   d
  cd  [ -1  -1   1   1 ]        rank 1  ⟹  HH^2 = 1 − 1 = 0
```

`HH^0 = 4 − rank δ^0 = 1`,  `HH^1 = (4 − rank δ^1) − rank δ^0 = 3 − 3 = 0`,  `HH^2 = 1 − rank δ^1 = 0`.
Three cross-checks agree: Euler `4 − 4 + 1 = 1`; Gerstenhaber–Schack (order complex of the diamond
poset = two triangles glued on `{1,4}` = contractible ⟹ `HH^{>0}=0`); Künneth (`A = kA_2 ⊗ kA_2`,
`kA_2` a hereditary tree ⟹ `HH^0=k`, `HH^{≥1}=0`). Homology instead gives `HH_• = (4,0,0)` (the
quiver is acyclic, so `C_1 = C_2 = 0` and `HH_0 = A/[A,A] = k^{|Q_0|} = k^4`); `A` is not symmetric,
so `HH_0 = 4 ≠ 1 = HH^0`.

## What is certified, and what is not
- Deep dims: certified for `k[x]/(x^a)` and the quantum CI (byte-oracle) to any depth.
- General `kQ/I`: computed, certified per instance by `d²=0` + order gate + bar window,
  **restricted to quadratic tips or monomial presentations**; a non-quadratic non-monomial
  presentation raises `NotImplementedError` (the `right_decomposition` stretch item lifts this).
- Operations (cup/bracket): transported to bar; certified only in the bar window.
  (Cap products exist on the bar engine `engine/tt_calculus.py` only — no CS transport wrapper yet.)
