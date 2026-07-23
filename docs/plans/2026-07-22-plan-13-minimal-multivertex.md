# Plan 13: Minimal A^e engine — multi-vertex (corner-typed) support

**Goal:** `engine/resolutions_minimal.py` silently returned a zero resolution on any
multi-vertex algebra. Fix it by building the minimal **projective** bimodule resolution
with corner terms `A^e·(e_v ⊗ e_w)`; local algebras remain the one-corner special case on
the unchanged kernel-accelerated path.

## Root cause (confirmed empirically, 2026-07-22)

`radical_basis` assumed A local: rad(A) = span{f_i : i ≠ t}. On multi-vertex input that
span contains the other vertex idempotents (verified: on kA_2 it returns `[e_2, a]` with
`e_2·e_2 = e_2`), so the fake rad(A^e) has codimension 1 and `rad(A^e)·ker ⊇` essentially
all of ker ⇒ the Nakayama generator selection returns 0 generators ⇒
`status="terminated"`, `rks = {0:1, 1:0, …}` — silently wrong (ker μ ≠ 0).

**Why "fix the radical" alone is NOT enough:** over a multi-vertex algebra, minimal *free*
A^e-covers do not terminate — the kernel of `A^e ↠ A` contains whole off-diagonal corner
projectives `Ae_v ⊗ e_w A` (v ≠ w), and free covers of projectives spawn projective junk
forever (even for `k × k`). Correct Tor but exponentially growing terms = useless engine.
The real fix is corner-typed **projective** covers.

## Design

- **Vertex detection (engine `Algebra`)**: `vertex_indices = {i : unit[i] == 1}` in the
  original basis, validated as orthogonal idempotents via T (`e_i e_j = δ_ij e_i`); the
  basis change to the f-basis leaves non-t basis vectors untouched, so path-type bases
  keep single-coordinate radical vectors (the `_rad_ab_pairs` invariant).
- **`radical_basis` (corrected)**: standard basis vectors at all non-vertex indices.
  Last-line-of-defense guard (any basis, any provenance): the multiplicative closure of
  the candidate span must be nilpotent within `m` steps, else raise `QuiverlabError`
  loudly (path-type basis required — rebuild via `Quiver.algebra`).
- **Corner machinery (pure Python, ≥ 2 vertices only)**:
  - corner basis `C_{vw}` = independent columns of right-mult by `ε_v ⊗ ε_w` on A^e
    (engine convention `(a⊗b)(p⊗q) = ap ⊗ qb` ⇒ right-mult by `(x⊗y)` is `(px)⊗(yq)`);
  - `P_0 = ⊕_v A^e·(ε_v⊗ε_v)` (ambient `(A^e)^{#vertices}`, generator `g_v = ε_v⊗ε_v`);
  - state `cur` = matrix from P_n **corner coordinates** to P_{n-1} flat ambient; kernel
    over corner coordinates (this is what excludes the off-corner junk);
  - kernel vectors flattened through `C` blocks → radK via the existing (kernel-capable)
    `_build_radK`; candidates = **corner components** `(ε_i⊗ε_j)·k` of kernel vectors
    (K = ⊕ ε K as left modules), greedy independent-modulo-radK with tags;
  - contraction: `A ⊗_{A^e} A^e·(ε_v⊗ε_w) ≅ e_w A e_v`; Dbar blocks are corner bases of
    A, entries by the existing `e_vv · α · e_uu` collapse generalized to corner vectors,
    coordinates recovered through per-corner expression matrices (reconstruction
    asserted — loud on failure).
- **Local path untouched**: `len(vertex_indices) == 1` routes through the existing code
  (kernels, goldens, deepen, batch identical).
- **`deepen`**: corner mode refuses loudly (`NotImplementedError`) — it is the
  checkpointing wrapper for the local open-zone scans; its manual state rebuild has no
  corner data. Unchanged for local.
- **`complexity` (invariants/scalar.py)**: caveat (a) rewritten — multi-vertex is now
  exact; the `linear_path_algebra(2) → 0` pin must now pass for the *right* reason
  (finite corner resolution `6,5,2,1,0`-style termination, not instant bogus termination).

## Validation matrix (tests/engine/test_minimal_multivertex.py)

1. `kA_2`: HH ≡ bar (4 primes), termination, `complexity → 0`.
2. Commutative square `kQ/(ab − cd)` (dim 9): HH_• = `[4,0,0]` ≡ bar live.
3. Line quiver `kQ/(abc, cde)` (dim 16): corner Betti = Bardzell chain counts
   `{0:6, 1:5, 2:2, 3:1, 4:0}` — independent re-derivation of Plan 12's straddle chain
   `abcde` by pure syzygy linear algebra.
4. Cyclic Nakayama `kZ_3/rad²` (dim 6, HH nonzero in all degrees): HH ≡ bar (4 primes).
5. Loud-guard test: a structure-constant algebra whose non-unit basis vector is
   idempotent (`k×k` with basis `{1, (1,0)}`) raises `QuiverlabError`.
6. Local regression: full existing engine suite unchanged.

## Status

Executed 2026-07-22 in-session (branch `plan-13-minimal-multivertex`); amendments
recorded inline in commits.
