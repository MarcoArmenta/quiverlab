# 05 — Resolutions

## The mathematics

Hochschild homology is Tor over the enveloping algebra: HH_n(A) = Tor_n^{A^e}(A, A), where
A^e = A ⊗ A^op. To compute it you need a projective resolution of A as an A^e-bimodule,
apply A ⊗_{A^e} (−), and take homology. Any resolution gives the same answer, so the game
is to use the *smallest* one that a given algebra admits. quiverlab ships several: the bar
complex (Chapter 04, huge but universal), the **minimal** A^e-resolution built by
syzygies (works for every algebra, smallest possible terms), and **Bardzell's**
resolution (closed-form and minimal for monomial algebras, reaching enormous depth).

## How it is represented — the `Resolution` protocol

A resolution backend is a **class** implementing a fixed interface `Resolution`
(`engine/resolutions.py`). The two primitives every homology backend must expose:

- `term_basis(alg, n)` — an ordered **list** of hashable generators of the n-th term; its
  length is dim P_n;
- `differential_matrix(alg, n, basis_n, index_nm1)` — the matrix of d_n as a numpy
  **int64 array** (a fixed-width-integer matrix) of shape `(dim P_{n-1}, dim P_n)`.

Two firm contract rules: matrices are **integer, never pre-reduced mod p** (a small prime
p dividing an entry carries genuine signal — torsion — that reducing early would erase);
and the differential is expressed in exactly the ordering that `term_basis` returned. Two
more primitives (`cochain_basis`, `coboundary_matrix`) serve cohomology, and two
(`sigma_chain_matrix`, `sigma_cochain_matrix`) serve the Nakayama-action layer of Chapter
06. `BarResolution` is the reference backend that delegates to the bar-complex code and is
the oracle every other backend is checked against.

## The minimal A^e-resolution (the syzygy stepper)

`engine/resolutions_minimal.py` builds the minimal projective bimodule resolution one
degree at a time, valid for *any* finite-dimensional path-basis A given by its structure
constants. The engine `AeEngine` precomputes, over F_p, the m^2 × m^2 matrices of
left-multiplication by each A^e basis element e_a ⊗ e_b (using
(e_a ⊗ e_b)·(e_p ⊗ e_q) = (e_a e_p) ⊗ (e_q e_b)).

**Local vs multi-vertex terms (Plan 13).** For a local algebra projective = free, so the
terms are free: `... -> (A^e)^{r_n} -> ... -> A^e -> A -> 0` (the original,
kernel-accelerated path). Over a multi-vertex algebra a minimal *free* resolution does
not exist — ker(A^e ↠ A) contains whole off-diagonal corner projectives Ae_v ⊗ e_wA, and
free covers of projectives spawn projective junk forever — so the engine builds
**corner-typed** terms `P_n = ⊕_j A^e·(ε_{v_j} ⊗ ε_{w_j})` (`_CornerContext`, pure
Python): kernels are computed over corner coordinates, generator candidates are the
corner components (ε_i ⊗ ε_j)·k of kernel vectors, and each generator carries its corner
tag. The vertex idempotents are read off the unit's 1-coordinates (validated against T);
a nilpotent-closure guard in `radical_basis` refuses non-path-type bases loudly instead
of returning a silently wrong resolution (the pre-Plan-13 multi-vertex failure mode).
Validation: HH ≡ bar on kA_2 / the commutative square / kZ_3/rad², and on the monomial
line quiver kQ/(abc, cde) the corner Betti numbers equal Bardzell's chain counts
6, 5, 2, 1, 0 — with the 1 the straddling overlap `abcde` (an independent syzygy-side
re-derivation of the Plan-12 chain).

`_advance_resolution` computes the next degree n from the current differential `cur`:

1. **Syzygies.** `nullspace_mod_p(cur, p)` finds `ker` = a basis of the kernel of the
   current map over F_p (these are the relations among the previous generators). If the
   kernel is empty the resolution *terminates* — A has finite Hochschild dimension.
2. **Minimal generators.** The kernel must be cut down to generators *modulo* the radical
   of A^e (minimality). `_build_radK` forms rad(A^e)·ker (the non-minimal part), and
   `_independent_modulo` greedily keeps exactly the kernel vectors that are independent of
   it — via one incremental row-reduction, not repeated rank tests. Those chosen vectors
   `gens` are the columns of d_n and their count is r_n. rad(A) is the span of the
   non-idempotent basis vectors (vertex idempotents read off the unit; for a local
   algebra this degenerates to "all non-unit basis vectors") — computed
   characteristic-independently by `radical_basis`, with a nilpotent-closure guard.
3. **Assemble d_n.** `_build_Dn` writes each generator into a column block, giving the new
   `cur`; the loop repeats.

Then `_contracted_complex` applies A ⊗_{A^e} (−) to every degree, and
`minimal_homology_dims` reads off dim HH_n = dim(A ⊗ P_n) − rank(dbar_n) − rank(dbar_{n+1}).

**Guards.** The host is small, so two budgets stop the build gracefully instead of being
killed by the operating system: `max_term_dim` caps a term's k-dimension `m^2·r_n`, and
`max_transient_bytes` *predicts* the size of the large transient `radK` array *before*
allocating it and stops if it would blow the budget. Either way `truncated_at` records the
last fully-known degree, and the returned dimensions are exact up to `truncated_at − 1`.
Because a minimal resolution depends on the characteristic, the whole thing is rebuilt per
prime; the default large prime 32003 is the faithful char-0 proxy and small primes expose
torsion.

## Checkpointed deepening

`engine/deepen.py` drives the same stepper for very deep runs on a cluster, one degree per
step, writing an **atomic checkpoint after every degree** so a timeout or node failure
costs at most one degree of recompute. A checkpoint is a `pickle` file (Python's binary
save format) whose payload is a dict containing literally: `n` (the last completed degree),
`cur` and `cur_r` (the current differential and its generator count), `rks` (the r_n so
far), `last_gens` (the previous degree's generator columns, needed to finalize the next
HH), `HH` (the dimensions computed so far), and `per_degree` (per-degree timing and memory
records). Resuming reads the latest checkpoint and rebuilds only the cheap `AeEngine`. A
`latest.txt` pointer only ever advances (a stale second writer cannot regress it), and
`deepen` predicts the next degree's walltime as a multiple of the last one's, stopping
*before* a degree that would overrun the time budget rather than being SIGKILL'd mid-degree
with a stale result. `finalize_only=True` re-emits the summary from the newest checkpoint —
the recovery path when a job died before writing its JSON.

## The Bardzell resolution (monomial algebras)

`engine/resolutions_bardzell.py` implements Bardzell's minimal projective bimodule
resolution for a monomial algebra kQ/I, exposed as a `Resolution` backend over a
`MonomialPresentation` (vertices, arrows as `(id, source, target)` tuples, relations as
path tuples). Its combinatorics are the graded **associated paths** AP^n (equivalently
Anick chains / n-ambiguities):

- AP^0 = vertices, AP^1 = arrows, AP^2 = the minimal relations;
- AP^n (n ≥ 2) = paths admitting a unique *left decomposition* p = u_0 u_1 ... u_{n−1}
  (u_0 a single arrow): each next block is the **shortest extension making the consecutive
  pair reducible** — the witness relation ends at the pair's end and *straddles* the block
  boundary; with mixed-length relations it can be a proper suffix of the pair rather than
  the pair itself (CS §3 / Anick minimality, corrected in Plan 12). Built recursively by
  `associated_paths`. Example, relations `{xx, yy, xyx}`: `xyxx = (x)(yx)(x)` — the pair
  `yx·x` is not a relation but contains `xx`; the exact-pair shortcut missed this chain
  (and `xyxyx`), giving AP^3 = 3 instead of the true 5 = the minimal-A^e Betti number.
- every associated path also has a unique *right decomposition* (`right_decomposition`,
  last block a single arrow, mirror-greedy): same set (CS §3 Prop.), different blocks —
  for `xyxx` it is `(xy)(x)(x)`.

The n-th term is P_n = ⊕_{p ∈ AP^n} A e_{o(p)} ⊗ e_{t(p)} A; after A ⊗_{A^e} (−) each
summand collapses to the loop space e_{o(p)} A e_{t(p)} (paths closing p), so a basis
element of the contracted term is a pair `(p, w)` — an associated path and a nonzero path
closing it (`term_basis`, `loops`). The differential `differential_matrix` alternates
between Bardzell's "big" map (even n: sum over (n−1)-associated subpaths) and "small" map
(odd n ≥ 3: `(+1)` drop the leftmost **right**-block `v_top`, `(−1)` drop the rightmost
**left**-block — CS §4 `f_n` even; using `u_0` in the first slot is correct only for
quadratic/palindromic relation sets), with n = 1 the commutator/augmentation map.
Because the resolution is *minimal* and its terms are the small associated-path sets rather
than the bar complex's m·(m−1)^n, it runs far past the bar wall — the bank record is exact
homology to **degree 1702** (a hanlab-bank record; not reproduced or verified in this
port). It is cross-checked entry-exact against the bar oracle on
k[x]/(x^a), cyclic Nakayama algebras, the radical-square-zero families, and (Plan 12)
the straddling-overlap presentation k⟨x,y⟩/(xx, yy, xyx) against bar **and** the
minimal-A^e engine simultaneously.

## The periodicity detector, and eventually-periodic families

For self-injective and other eventually-periodic algebras the HH sequence eventually
repeats. Two homology-only *family* backends in `engine/resolutions_periodic.py` give a
caller the small resolution by name: `CyclicNakayamaResolution(n, ell)` (a thin verified
wrapper over Bardzell) and `QuantumCIResolution(c)` (a wrapper over the Chouhy–Solotar
closed form — **dormant in this port; it lands with Plan 04**). Separately, once a
sequence of dimensions is in hand, `complexity_diagnostic` (`engine/scan2.py`) is the
*detector*: it takes finite differences to spot polynomial growth, checks for trailing
zeros (eventual vanishing, the Han-counterexample signature, complexity 0), and tests for
eventual periodicity `seq[i] == seq[i−p]` returning a label like `eventually_periodic_p2`.
It is honest — it reports only what the *computed* degrees show, never a proof of the tail.

## The fast GF(p) engine, and its pure-Python fallback

Over a prime field the engine replaces exact-domain rref with numpy int64 matrices and
rank mod p. `rank_mod_p` (`engine/hh_engine.py`) is dense Gaussian elimination over F_p
(modular inverse via `pow(a, p−2, p)`); `sparse_rank_mod_p` (`engine/linalg_fast.py`) does
the same on a dictionary-of-columns for the very sparse deep differentials, and
`rank_mod_p_auto` dispatches to whichever fits. The genuinely hot loops (nullspace, rank,
the radK and Dn matvecs — matrix-times-vector products) also have **numba-compiled kernels** (`engine/_kernels.py`) —
just-in-time-compiled machine code. The flag `USE_KERNELS` is true only when numba is
importable and not disabled; setting the environment variable `QUIVERLAB_NO_NUMBA=1`
forces the pure-Python twins. Those twins are not just a fallback — they are the permanent
oracle, and a parity test runs both. (Overflow is controlled: the matvec kernels
accumulate up to m^2 products before one reduction, safe while m^2·(p−1)^2 < 2^63, which
holds comfortably for p = 32003 at the target dimensions.)

## The adapter — public Algebra into engine form

The engine speaks numpy int64, the public `Algebra` speaks field elements. `to_engine`
(`engine/adapter.py`) bridges them: it **refuses any non-prime field loudly** (a
`FieldError` — the fast engine is a GF(p) accelerator, the bar path serves all fields),
then extracts integer coefficients by looping `T[i, j, t] = int(vec[t])` over the public
table and `int(c)` over the unit, producing an `engine/hh_engine.py` `Algebra` (which
additionally changes to its own unit-adapted basis so that A/k.1 is "drop coordinate t").
`engine_cohomology_dims` / `engine_homology_dims` guard the exponential bar-basis size
against `max_cells` exactly as the pure oracle does, then return plain `list[int]`.

## A worked micro-example — a Bardzell term for k[x]/(x^3)

`MonomialPresentation.truncated_polynomial(3)` has one vertex `"v"`, one loop arrow `0`,
and the single relation `(0, 0, 0)` = x^3. Running it: AP^1 = `[(0,)]` (the arrow x),
AP^2 = `[(0, 0, 0)]` (the relation x^3), AP^3 = `[(0,0,0,0)]`. The contracted degree-0
term basis is `[(('v','v'), ()), (('v','v'), (0,)), (('v','v'), (0,0))]` — the vertex loop
paired with the loops 1, x, x^2 closing it (dim 3 = dim A). Each higher term keeps this
constant size, which is exactly why Bardzell reached degree 1702 in the hanlab bank where the bar complex —
dim C_n = 3·2^n — chokes. (These bases were produced by running the code.)

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| the resolution interface | `engine/resolutions.py` | `Resolution`, `BarResolution`, `_default` |
| minimal A^e engine, syzygy step | `engine/resolutions_minimal.py` | `AeEngine`, `_advance_resolution`, `minimal_homology_dims` |
| minimal-resolution guards | `engine/resolutions_minimal.py` | `minimal_resolution` (`max_term_dim`, `max_transient_bytes`) |
| checkpointed deepening | `engine/deepen.py` | `deepen`, `_save_ckpt`, `_load_ckpt` |
| Bardzell resolution, associated paths | `engine/resolutions_bardzell.py` | `BardzellResolution`, `MonomialPresentation.associated_paths` |
| eventually-periodic family backends | `engine/resolutions_periodic.py` | `CyclicNakayamaResolution`, `QuantumCIResolution` |
| periodicity / complexity detector | `engine/scan2.py`, `engine/scan3.py` | `complexity_diagnostic`, `complexity_of` |
| fast rank over F_p | `engine/hh_engine.py`, `engine/linalg_fast.py` | `rank_mod_p`, `sparse_rank_mod_p`, `rank_mod_p_auto` |
| numba kernels + fallback flag | `engine/_kernels.py` | `USE_KERNELS`, `rank_mod_p_kernel`, `nullspace_kernel` |
| public Algebra -> engine | `engine/adapter.py` | `to_engine`, `engine_cohomology_dims`, `engine_homology_dims` |
