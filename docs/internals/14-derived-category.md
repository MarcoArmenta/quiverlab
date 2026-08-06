# 14 — The derived-category surface

## What this computes

`derived/` (Plan 43) is the bounded derived category `D^b(mod A)` made computable, for a
finite-dimensional `A = kQ/I`. It reifies morphisms `Hom_{D^b}(X, Y[n])` as honest chain
maps, computes the derived Auslander–Reiten translate `τ_{D^b} = ν[−1]` on perfect
complexes, **verifies** tilting complexes and builds their endomorphism algebra `End(T)`
(the Rickard derived-equivalent algebra), and assembles a *necessary-condition* derived
fingerprint. It is a thin exact-linear-algebra layer over P37 (`End`/`ModuleHom`), P38
(Cartan/Coxeter), P39 (complexes / hyper-Hom) and P41 (Nakayama / corner-transpose) — no
new math engine, public surface only.

The whole surface is honest about a hard fact: **deciding derived equivalence is not
algorithmic.** Everything here is a verifier or a necessary-condition invariant; nothing
claims to decide equivalence.

## Complex conventions (`modules/complexes.py`, the P39 substrate)

`ChainComplex(terms, dmats, check=True)` is **homological**: `d_n: C_n → C_{n-1}`, with
each differential matrix written **rows = target `C_{n-1}`, columns = source `C_n`** —
byte-identical to `modules.resolution`. A missing degree is the zero module; cohomological
indexing is presentation-only (`C^n := C_{-n}`). Construction is self-certifying: `check`
validates each `d_n` as a `ModuleHom` and asserts every composite `d_n ∘ d_{n+1}` is zero
(loud otherwise). `shift(k)` sends degree `n` to `n − k` and multiplies each differential
by `(-1)^k`; `truncate`, `homology_dims`, and `homology(n) = Z_n/B_n` follow. A
`ChainMap(src, tgt, components, check=True)` validates each component as a module map **and**
that every square commutes (`d^tgt f_n = f_{n-1} d^src`); `.then(g)` is left-to-right
composition. The Weibel sign for the Hom total complex is
`Hom^n(X, Y) = ⊕_p Hom_A(X_p, Y_{p-n})` with `(δf)_p = d^Y f_p − (-1)^n f_{p-1} d^X`; the
homology dims are sign-independent (Weibel `ε = -1` and `ε = +1` give isomorphic cochain
complexes), and the code takes `ε = -1`.

## Reified hyper-Hom (`homs.py`)

`hyper_hom_basis(X, Y, n)` returns a **list of `ChainMap`s** that is a basis of
`H^n(Hom^•(X, Y))`; for `X` perfect this is `Hom_{D^b(mod A)}(X, Y[n])`. The reification is
exact: a class in `Hom^n(X, Y)` has components `f_p: X_p → Y_{p-n}` landing in
`Y.shift(n).term(p)`, and its cocycle condition is *precisely* the chain-map square for
`X → Y[n]` — so `ChainMap(X, Y.shift(n), comps, check=True)` passes on a cocycle and fails
on a non-cocycle. Each returned map is built from a **canonical coset representative** of
`ker(δ^n) / im(δ^{n-1})`: `_coset_reps` picks classes with `independent_modulo`, then
`_reduce_mod_span` subtracts only coboundary-span vectors (an RREF reduction that keeps a
cocycle a cocycle) — the free-variables-zero canonicalization adapted to a *column* span,
deliberately **not** `reduce_mod_nullspace(r, transpose(B))`, which would move `r` out of
the kernel. A self-cert closes the loop: the reified class count must equal
`hyper_hom_dims(X, Y, n)`, else it raises.

## The derived AR translate (`tau.py`)

`τ_{D^b} = ν[−1]` where `ν = D Hom_A(−, A)` is the Nakayama functor, applied **termwise** on
a perfect (projective) complex: each projective term `⊕P_v` maps to
`ν(⊕P_v) = D Hom_A(⊕P_v, A) = D(⊕Ae_v) = ⊕I_v`, each differential to its `ν`-image via the
shared **corner-transpose** (`_corner.py::corner_transpose`, factored out of
`duality._presentation_transpose`), and the whole complex is then shifted by `−1`.
`tau_Db_minus` is the inverse, `ν^{-1}[+1]`, on a perfect complex of injectives (the output
shape of `tau_Db`); the round-trip `tau_Db_minus(tau_Db(X))` is certified a quasi-iso.

Two design facts:

- **Happel's finite-gl.dim guard.** The Serre functor / AR triangles of `D^b(mod A)` exist
  iff `gl.dim A < ∞`. `_require_finite_gldim(A)` reads `global_dimension(A)` and refuses
  loudly unless it is certified finite (`g.exact`); `k[x]/(x²)` is the pinned negative
  case. `tau_Db` also refuses a non-perfect input (resolve with `projective_model` first).
- **Basis discipline (the P39/P41 mismatch rule).** The injective terms are built as
  `dualize(Hom_A(term, A))` — exactly what `nakayama_functor` produces — so their k-basis
  is byte-consistent with the corner-transpose differentials; a `builders.injective` basis
  is never mixed with a corner-transpose basis. Kind-tagged provenance
  (`_term_provenance = ("injective"|"projective", …)`) stops a projective complex from
  masquerading as injective for `tau_Db_minus`, and vice versa. Correctness never rests on
  the bookkeeping: the output self-certifies (`ChainComplex(check=True)` gives `d∘d = 0`)
  and is pinned against the trusted module `τ`.

The K₀ arbiter is the identity `χ(τ_{D^b} X) = c · χ(X)` with `c = −C · C^{-ᵀ}` (the
K₀-*action* Coxeter matrix). This is deliberately **not** P38's `coxeter_matrix`
(`= −C^{-ᵀ} C`, its conjugate — same characteristic polynomial, different action), and it
pins the whole termwise `ν` bookkeeping against the concrete module `τ` on kA₂/kA₃.

## Tilting verifier and `End(T)` (`tilting.py`)

`is_tilting_complex(summands)` is a **verifier** for a given list of perfect summands — it
decides, it never searches. It returns a `TiltingReport` (`is_tilting`, `rigid`,
`generates`, `window`, `g_matrix`, `det`) built from two checks:

1. **Rigidity** — `Hom_{D^b}(T, T[n]) = 0` for all `n ≠ 0`, scanned over the *exact* window
   outside which hyper-Hom is provably the zero cochain group; the window is reported
   honestly in `.window`.
2. **Generation** of `K^b(proj A)` — the K₀ g-matrix (rows = summand Euler characteristics
   `χ`) must be square and unimodular (`det = ±1`).

`is_tilting = rigid and generates` (Rickard). `end_algebra_of_complex(summands)` builds
`End_{D^b}(T) = ⊕_{i,j} Hom_{D^b}(T_i, T_j)` as a **structure-constant `Algebra`** — the
degree-0 hyper-Hom classes composed with `ChainMap.then`, reduced to canonical homotopy
representatives, handed to `Algebra.from_structure_constants(..., check=True)` (associativity
and unit validated). `corner_cartan_of_complex` is the corner-Cartan of `End(T)` (the entry
at row i, column j is `dim_k Hom_{D^b}(T_j, T_i)`), whose orientation is pinned by the
`T = A` oracle (`corner_cartan_of_complex(A) = cartan_matrix(A)`,
`End(A_A) ≅ A`). `two_term_silting_from_presentation(M)` returns the 2-term complex
`P_1 → P_0` of `M` with its rigidity report — the bridge the P45 τ-tilting engine consumes
(a single object has `generates = False` by design; the consumer reads `.rigid` plus the
g-vector, not `.is_tilting`).

## The fingerprint — what it is, and is not (`fingerprint.py`, `block.py`)

`derived_fingerprint(A, top=4)` returns a dict of **necessary-condition** derived
invariants: `coxeter_polynomial`, `cartan_det`, `cartan_smith` (the `GL_n(ℤ)`-equivalence
invariant factors — coarser than ℤ-congruence, and the docstring says so), the Hochschild
`hh_cohomology_dims` / `hh_homology_dims`, `cyclic_dims`, `center_dim`, and `gl_dim`. Every
field is wrapped by `_field`, which captures a `QuiverlabError` (including a `DepthLimitError`
blow-up of the generic cyclic-homology mixed complex off GF(p)) as `{"error": …}` — an honest
per-field non-answer, never a crash.

`compare_fingerprints(fa, fb)` speaks only in `"distinguished"` /
`"not distinguished by these invariants"` — **never** "(in)equivalent". Equal fingerprints
are a necessary, not sufficient, condition: the **8-vertex cospectral trees** are the pinned
counterexample (equal Coxeter polynomial, Cartan, HH and centre, yet *not* derived
equivalent). A field that errored on either side is surfaced in `incomparable_fields`, never
silently dropped. The fingerprint is **never a decider** — that is the whole honest-scope
point.

`block.py::derived_fingerprint_block(A, top=4)` is the single-algebra `derived_fingerprint`
scalar compute kind (schema v1), built once so all three runners return byte-identical
blocks; its `scope` string states "equal values are a necessary condition for derived
equivalence, not a proof". The **two-algebra compare panel is deferred** to a post-release
successor (it needs a second-algebra request field — a schema change); see the
[v0.2.0 GUI-deferral ledger](../verification.md#v020-gui-deferral-ledger).

## The oracles

- **Self-cert** — `ChainMap(check=True)` on every reified class, `ChainComplex(check=True)`
  for every `τ_{D^b}` output (`d∘d = 0`), the class-count `= hyper_hom_dims` identity, the
  `tau_Db_minus ∘ tau_Db` quasi-iso round-trip, the honest per-field fingerprint errors,
  and `from_structure_constants(check=True)` for `End(T)`.
- **Cross-engine** — `hyper_hom_basis` count `= A.ext(M, N, n)` on a projective-resolution
  source; `τ_{D^b}` of a non-projective module concentrated in degree 0 and isomorphic to
  the trusted module `τ` over kA_n.
- **Literature** — the K₀ identity above; the kA₂ APR tilt `P_1 ⊕ S_1` is tilting with
  `End(T)` the reoriented A₂ (`= A^op`); D₄ vs A₄ distinguished by the Coxeter polynomial;
  the cospectral-trees non-distinction.
- **Live QPA** — a probe-first battery: `TauOfComplex` on a `ProjectiveResolution` does not
  script cleanly through libgap (the P39 Ch. 10 complex-scripting hazard), so the crosscheck
  falls back to the documented module-level route (`τ_{D^b}` homology in degree 0 vs QPA
  `DTr(M)`), recording exactly what was compared — never a silent skip.

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| bounded complexes, chain maps, cones | `modules/complexes.py` | `ChainComplex`, `ChainMap`, `hyper_hom_dims`, `projective_model` |
| reified hyper-Hom classes | `derived/homs.py` | `hyper_hom_basis` |
| the corner-transpose `Hom_A(−, A)` | `derived/_corner.py` | `corner_transpose` |
| derived AR translate `τ_{D^b}` / inverse | `derived/tau.py` | `tau_Db`, `tau_Db_minus`, `_require_finite_gldim` |
| tilting verifier, `End(T)`, corner-Cartan | `derived/tilting.py` | `is_tilting_complex`, `TiltingReport`, `end_algebra_of_complex`, `corner_cartan_of_complex`, `two_term_silting_from_presentation` |
| the necessary-condition fingerprint | `derived/fingerprint.py` | `derived_fingerprint`, `compare_fingerprints` |
| the `derived_fingerprint` block | `derived/block.py` | `derived_fingerprint_block` |
