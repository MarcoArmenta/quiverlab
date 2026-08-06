# 13 — Spectral sequences

## What this computes

Given a **bounded filtered complex** or a **bounded double complex** of modules over an
exact `Domain`, the `specseq/` engine (Plan 42) builds the associated spectral sequence:
the pages `E^r`, their induced differentials `d^r`, and a convergence report that says
which page stabilizes and whether it degenerates. Everything is exact linear algebra over
`fields.linalg` — no floats (the `src/` AST gate) — and every page representative is
byte-reproducible. Four presets wrap the engine around named homological questions: the
Hochschild `(b, B)` bicomplex, the radical / associated-graded filtration, and the
Cartan–Eilenberg / Grothendieck change-of-rings sequence.

## Conventions (the one place they are fixed)

The engine is **homological** throughout: an increasing filtration `F_p` (Weibel, *An
Introduction to Homological Algebra*, 5.4.6), total degree `n = p + q`, and every
differential written `d_n: V_n → V_{n-1}` with matrices **rows = target, columns =
source** — byte-identical to the `modules.complexes.ChainComplex` and
`modules.resolution` layout. A **cohomological** source (the Grothendieck Hom double
complex) is stored with *negated* total degree, position `(-p, -q)` — the P39
`C^n := C_{-n}` discipline — so it lands in the same homological machinery and never
forces a reindex on the caller.

## The objects (how they are represented)

Every subspace in the engine is a list of **coordinate columns over the total-degree
basis** `Tot_n`; the shared exact-subspace helpers live in `_subspace.py` (`colmat`,
`image`, `intersect`, `preimage_selecting`, `reduce_to_independent`, `span_dim`).

- **`FilteredComplex`** (`filtered.py`) — a bounded homological complex plus an
  increasing, exhaustive (`F_top = whole`), Hausdorff (`F_{-large} = 0`) subcomplex
  filtration; `from_chain_complex(X, filt, lo)` wraps a P39 `ChainComplex`. Accessors
  `piece(n, p) = F_p C_n`, `dmat(n)`, `levels()`, and `total_homology_dims()` (the rank
  formula `dim H_n = dim V_n − rank d_n − rank d_{n+1}`).
- **`DoubleComplex`** (`double.py`) — a bounded homological double complex with horizontal
  `d_h` and vertical `d_v`; `__init__(..., check=True)` gates anticommutativity
  (`d_h d_v + d_v d_h = 0`) as an exact matrix identity. `total()` assembles
  `Tot_n = ⊕_{p+q=n} D_{p,q}`; `column_filtration()` / `row_filtration()` return the two
  associated `FilteredComplex`es.
- **`SpectralSequence`** (`pages.py`) — the page engine of a `FilteredComplex`. Holds the
  filtration, `width` (number of levels), `height`, and a memoized page cache.
- **`Page`** / **`Subquotient`** — one page `E^r` and one position `E^r_{p,q}` (its `dim`
  and canonical `reps`). `Page.grid()` prints an M2-`netPage`-style ASCII grid (`p`
  across, `q` up).
- **`ConvergenceReport`** (`convergence.py`) — the standing certificate's result.

## The pages, step by step (`pages.py`)

The subquotient formulas are pinned verbatim in the module header (Weibel 5.4.6):

    Z^r_{p,q}    = { x in F_p C_{p+q} : d(x) in F_{p-r} C_{p+q-1} }
    Bdry^r_{p,q} = F_p C_{p+q} ∩ d( F_{p+r-1} C_{p+q+1} )   (= d Z^{r-1}_{p+r-1,q-r+2})
    E^r_{p,q}    = Z^r_{p,q} / ( Z^{r-1}_{p-1,q+1} + Bdry^r_{p,q} )
    d^r : E^r_{p,q} → E^r_{p-r, q+r-1},  induced by d.

- `_Zr(p, q, r)` builds `Z^r` by restricting `d_n` to `F_p C_n`, then selecting the
  preimage that lands inside `F_{p-r} C_{n-1}` (`preimage_selecting`); a vanishing `d_n`
  makes the condition vacuous and returns the whole piece.
- `_bdry(p, q, r)` intersects `F_p C_n` with the image of `d_{n+1}` off `F_{p+r-1} C_{n+1}`.
- `_cell(p, q, r)` returns `(reps, denom)`: the denominator is
  `Z^{r-1}_{p-1,q+1} + Bdry^r_{p,q}`, and the representatives are the columns of `Z^r`
  that grow the rank over the denominator, picked in deterministic `rref` column order by
  `modules.linalg_mod.independent_modulo` — the CS/Plan-17 canonicalization mandate applied
  to pages, so the reps are **byte-reproducible** run to run.

**The one arbitrated index.** The boundary sub-object uses `F_{p+r-1}`, not the
plan-brief's `F_{p+r-2}` (which is one filtration step short and makes `E^1_{0,0}` of the
trivial one-step filtration come out `dim V_0` instead of `H_0`). The choice is
*arbitrated*, not assumed: both the trivial-filtration base-case test and the standing
`E_∞ == H(Tot)` self-certificate fail under the short form and pass under this one.

**The differential.** `_dr_matrix(p, q, r)` (exposed as `Page.differential(p, q)`) is the
lift-apply-reduce map `E^r_{p,q} → E^r_{p-r,q+r-1}`, **rows = target reps, columns = source
reps**: for each source representative it applies the honest complex differential, expresses
the image in the target's rep-plus-denominator basis by an exact `solve`, and canonicalizes
the coefficients with `reduce_mod_nullspace` (the unique free-variables-zero coset
representative) before keeping the target-rep coordinates.

## The standing self-certificate — `E_∞` totals equal `H(Tot)`

`SpectralSequence.__init__` calls `certify_convergence(self)` (`convergence.py`), which
checks the rank identity

    Σ_{p+q=n} dim E_∞^{p,q}  ==  dim H_n(Tot)   for every total degree n

and raises loudly on any mismatch (a page/filtration bookkeeping bug). The `E_∞` page is
reached by `e_infinity_page = max(width, height) + 1` (generous for a bounded filtration).
`ConvergenceReport` carries `e_infinity_page`, `degenerates_at` (the least `r ≥ 1` whose
per-cell dims already equal `E_∞`'s — degeneration is decidable by rank), and `abutment`
(`{n: dim H_n(Tot)}`), with `collapse()` (`degenerates_at ∈ {1, 2}`) and `prose()`.

Distinct from the stabilization page is `certified_window` (a per-instance attribute set
only by the presets that certify a finite abutment window against an *external* oracle):
the `(b, B)` and Cartan–Eilenberg presets build one degree deeper than they report, and
`certified_abutment(n)` **refuses loudly** to read an abutment degree outside the certified
window — a truncated double complex is silently wrong out there, so the reads are gated,
not guessed.

## The four presets (`presets.py`)

1. **`hochschild_bB_ss(A, top)`** — the first-quadrant `(b, B)` bicomplex
   `D_{p,q} = C_{q-p}` on the unit-adapted bar basis, vertical `b`, horizontal Connes `B`.
   The two differentials already anticommute (the mixed-complex identity `bB + Bb = 0`), so
   no sign adjustment is needed; the abutment is cyclic homology `HC_*(A)`. The exponential
   bar basis (`dim C_n = m(m-1)^n`) is guarded up front by length arithmetic
   (`_guard_bB_cells`), raising `DepthLimitError` before any matrix is built.
2. **`radical_filtration_ss(X)`** — the associated-graded filtration
   `F_p X_n = X_n · rad^{max(0,-p)}` of a P39 `ChainComplex`, radical powers iterated in
   place (`_rad_powers`). A semisimple complex collapses at `E_1`; for a Koszul algebra the
   minimal resolution is linear and the sequence degenerates early — the exact page is
   **arbitrated per instance** (`degenerates_at` by rank), pinned at `E_2` on kA₃/kA₄,
   never forced.
3. **`cartan_eilenberg_ss(A, B, M, N)`** — the change-of-rings sequence
   `E_2^{p,q} = Ext_B^p(M, Ext_A^q(B, N)) ⇒ Ext_A^{p+q}(M|_A, N)` for an **admissible
   quotient** `B = A/I'` (same quiver, `rel(A) ⊆ rel(B)`; `_assert_change_of_rings` gates
   it).
4. **`grothendieck_double_complex(M, U, N, p_len, q_len)`** — the double complex the
   Cartan–Eilenberg sequence is built on.

**CE is the `U = B` Grothendieck case.** `cartan_eilenberg_ss` builds its double complex by
calling `grothendieck_double_complex(M, B, N, …)`. The general `(B, A)`-bimodule / Eilenberg–
Watts Grothendieck sequence is **not implemented this release** — `grothendieck_double_complex`
refuses `U is not B` loudly. In the `U = B` case the term
`Hom_B(Q_p, Hom_A(B, J^q))` collapses by the change-of-rings adjunction to
`Hom_A(res_A Q_p, J^q)`, and *that* is exactly what the builder assembles:
`D^{p,q} = Hom_A(res_A Q_p, J^q(N))` for `Q_•` the minimal `B`-projective resolution of `M`
and `J^•` the minimal `A`-injective coresolution of `N`. It is stored cohomologically at
`(-p, -q)`, with horizontal = precompose with `d^Q` and vertical = `(-1)^p ·` postcompose
with `d^J` (the Koszul sign the anticommutativity gate confirms).

**Grothendieck acyclicity is a per-instance hypothesis check.** After building the sequence,
`cartan_eilenberg_ss` certifies the abutment against the module-`Ext` oracle degree by degree
over the window `[0, min(p_len, q_len) − 1]`; if the `E_∞` total differs from
`Ext_A^n(M|_A, N)` it **refuses loudly** — either the change-of-rings acyclicity hypothesis
`Ext_B^{>0}(M, Hom_A(B, J^q)) = 0` fails for this instance, or the truncation is too shallow —
never a wrong abutment. On success it sets `certified_window`.

## The no-code surface

`block.py::specseq_block(A, top)` is the shared `ss_hochschild` compute kind (algebra-only,
schema v1) driven by both runners (`hpc/spec.py` and the Pyodide twin `docs/gui/runner.py`);
it returns the `E_∞` page, the trimmed grid, the abutment, the degeneration page, and the
convergence prose, catching the bar-basis `DepthLimitError` as a clean error block. The
Cartan–Eilenberg / Grothendieck / radical presets are library + HPC-config accessible this
release; their no-code GUI is a named post-release successor (see the verification page's
[v0.2.0 GUI-deferral ledger](../verification.md#v020-gui-deferral-ledger)).

## The oracles

- **Self-cert** — `E_∞ == H(Tot)` at every construction, and the per-cell rank identity
  `dim E_{r+1} = dim E_r − rank(d_r out) − rank(d_r in)`; these arbitrate the boundary
  index and the double-complex sign.
- **Cross-engine** — the `(b, B)` `E_∞` total against `HC_*`, and the Cartan–Eilenberg
  abutment against `modules.ext.ext_dims` (the acyclicity certificate above).
- **Literature** — the closed-form `k[x]/(x^a)` (which forces a nonzero higher `d_r` for
  `a ≥ 3`), the ground-field `HC = 1, 0, 1, 0, …`, and the arbitrated Koszul `E_2`
  collapse on kA₃/kA₄.
- **Macaulay2 (`m2`)** — the commutative Koszul double complex's `E_∞` totals against M2's
  total-complex homology (M2 1.26's `SpectralSequences` package is unscriptable — it rides
  the removed `ChainComplex` type — so only the convergence target is compared, not the
  page grid).
- **QPA** — none: QPA has no spectral-sequence surface (stated in honest scope).

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| exact subspace layer (column spans) | `specseq/_subspace.py` | `colmat`, `image`, `intersect`, `preimage_selecting`, `reduce_to_independent` |
| filtered complex + total homology | `specseq/filtered.py` | `FilteredComplex`, `piece`, `from_chain_complex`, `total_homology_dims` |
| double complex + total / filtrations | `specseq/double.py` | `DoubleComplex`, `total`, `column_filtration`, `row_filtration` |
| pages, subquotients, `d^r` | `specseq/pages.py` | `SpectralSequence`, `Page`, `Subquotient`, `_Zr`, `_bdry`, `_cell`, `_dr_matrix` |
| the convergence certificate | `specseq/convergence.py` | `certify_convergence`, `ConvergenceReport` |
| the four presets | `specseq/presets.py` | `hochschild_bB_ss`, `radical_filtration_ss`, `cartan_eilenberg_ss`, `grothendieck_double_complex` |
| the `ss_hochschild` block | `specseq/block.py` | `specseq_block` |
