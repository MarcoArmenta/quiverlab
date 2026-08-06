# Plan 42: Spectral-Sequence Engine + Four Presets — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A first-class **spectral-sequence engine** over `fields.linalg`: bounded
filtered complexes and double complexes, exact pages `E_r = Z_r/(Z_{r-1}+B_{r-1})`
(Weibel §5.4) with `reduce_mod_nullspace`-canonical (byte-reproducible)
representatives, `d_r` differentials by lift-apply-reduce, an M2-`netPage`-style
page grid, and a convergence report — plus **the four presets Marco chose**
(metaplan §5 P42): Cartan–Eilenberg change-of-rings, Grothendieck (general, via a
`(B,A)`-bimodule), associated-graded / radical filtration, and the Hochschild
`(b,B)` bicomplex. Custom filtered / double complexes are the public substrate.
Every construction carries **the standing self-certificate**: `E_∞` totals per
total degree equal the homology of the total complex — a rank identity, loud on
mismatch. P43 (derived category) reuses this engine's page surface.

**Architecture:** New package `src/quiverlab/specseq/` (there is NO existing
spectral-sequence or double-complex type — 2026-08-05 exploration). The engine's
**internal currency is a Domain vector-space complex**: `(support degrees,
dim[n], dmat[n])` with `dmat[n]` the **homological** `d_n: V_n -> V_{n-1}`,
rows=target/cols=source — byte-for-byte P39's `ChainComplex._dmats` layout, and
engine/cyclic's total-differential layout. A P39 `ChainComplex` of `A`-modules is
consumed by **forgetting** its module structure to this vector-space complex
(`term(n).dim` + `differential(n).matrix`); a raw double complex produces one
directly. **Homological throughout** (matching P39 and `engine/cyclic.py`);
cohomological presets (Grothendieck's Hom double complex) negate their indices at
construction — presentation-only, exactly P39's `C^n := C_{-n}` discipline. Pages,
`Z_r`/`B_r`, intersections, preimages, and `d_r` are all exact subspace linear
algebra over `fields.linalg` (`rref`/`rank`/`nullspace`/`solve`/
`reduce_mod_nullspace`); canonical page reps use the greedy independent-modulo
filter of `hochschild/cyclic.py::_capture_generic_reps` (the tested precedent).

**Tech Stack:** `fields.linalg` (exact Gaussian elimination + the canonical
coset representative), `modules/linalg_mod` (`matmul`/`mat_rank`/`kernel_columns`/
`cols_to_matrix`/`solve_columns`), P39 `modules/complexes.py`
(`ChainComplex`/`ChainMap`), the module engines `modules/{ext,hom,injective,
resolution,tor}.py`, `hochschild/cyclic.py` (`connes_B_matrix`/`boundary_matrix`
over any Domain), `families/{trivial_extension,tensor}.py`. No floats in `src/`.

## Global Constraints

- **P39 must be merged to `dev` before this plan runs** — `FilteredComplex`
  consumes `ChainComplex`/`ChainMap` and the radical-filtration preset takes a
  `ChainComplex` input. Branch `plan-42-spectral-sequences` off `dev` AFTER P39
  lands. (P39 is the critical-path predecessor: metaplan §4 `P37 -> P39 -> P42`.)
- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **Test bucket = deep, via one conftest edit (Task 10, Step 1, done FIRST).**
  A new top-level `tests/specseq/` lands in the *fast* bucket by
  `tests/conftest.py` default; the SS engine runs real resolutions/injectives, so
  add `"specseq"` to `_DEEP_DIRS` (`tests/conftest.py:50`). The partition audit
  `tests/release/test_markers.py` is **self-adjusting** — it derives every bucket
  from a live `pytest --collect-only -m <expr>` sweep (no hardcoded dir list), so
  the conftest edit alone keeps it green; the task step is to make the edit and
  RUN `tests/release` to confirm. Do this edit before writing any test so the
  batteries collect deep from the first commit.
- **Homological convention is pinned** (verbatim, do not silently reindex):
  every vector-space complex, filtered complex, and double-complex total is
  `d_n: V_n -> V_{n-1}`, rows=target. A cohomological source (the Grothendieck Hom
  double complex) is stored with negated total degree; its oracle compares the
  abutment at the negated degree, stated at the call site.
- **The standing self-certificate runs on EVERY `SpectralSequence` construction**
  (`convergence.py`): `sum_{p+q=n} dim E_∞^{p,q} == dim H_n(Tot)` for every total
  degree `n` — a rank identity that must hold by construction; a mismatch raises
  `QuiverlabError` (a bug in the page/filtration bookkeeping), never a silent
  wrong page. This is the strong-convergence guarantee (bounded filtration of a
  bounded complex — metaplan §2 SS brief) made into a live gate.
- All validation loud (`QuiverlabError`); `check=False` fast paths only for
  internally-constructed data whose certificate is asserted separately.
- **Grothendieck-SS acyclicity is a per-instance hypothesis check** with a loud
  refusal path (metaplan §6): the preset probes `Ext_B^{>0}(M, Hom_A(U, J^q)) = 0`
  and refuses with a `NotImplementedError`-style message naming the failed
  hypothesis when it does not hold — never returns a wrong abutment.
- Plan-32 oracle-class markers, per battery file at module level:
  `oracle_selfcert` for the rank identities (`d_r∘d_r=0`, `E_{r+1}=H(E_r,d_r)`,
  `E_∞` totals == total homology, radical-SS converges to `H(X)`);
  `oracle_crossengine` for the CE/Grothendieck `E_∞` total == module `Ext` and the
  Hochschild `E_∞` total == `HC` batteries; `oracle_literature` for the closed-form
  `k[x]/(x^a)`, the ground-field `HC = [1,0,1,0,...]`, and the arbitrated Koszul
  degeneration; `m2` (the fifth class) for the one Macaulay2 `SpectralSequences`
  crosscheck; `qpa` — **none** (QPA has no spectral-sequence surface; stated on the
  verification page, no test claims it).
- **GUI story (metagoal 1, honest scope):** ONE new scalar/range compute kind
  `ss_hochschild` (algebra-only, schema v1) — the `(b,B)` SS pages/dims +
  convergence prose, via `specseq_block`, through all seven touchpoints, ONE
  golden. The CE / Grothendieck / radical presets are **API + HPC-config
  accessible** this release; their no-code GUI needs new request fields (a
  second/third module, an admissible quotient, a filtration choice) and is
  **DEFERRED to P50 integration**, the reason recorded on the verification page —
  the same pattern as P39's GUI deferral (the standing "every plan ships a no-code
  story" rule is satisfied by `ss_hochschild` now + the named successor for the
  rest, never a silent skip).
- Every merge updates `docs/verification.md` (new oracle rows + recounted class
  table, `tests/release/test_oracle_classes.py` green) and the README line.
  **Mid-merge-train count note:** the suite counts below are illustrative — the
  orchestrator recounts from a live `pytest --collect-only` at merge time and P50
  does the final recount; do not hardcode a number this plan cannot verify.
- Conventional commits; green at every commit.

---

### Task 1: `FilteredComplex` — the filtered vector-space complex

**Files:**
- Create: `src/quiverlab/specseq/__init__.py`, `src/quiverlab/specseq/filtered.py`
- Test: `tests/specseq/test_filtered.py`

**Interfaces:**
- Consumes: `fields.linalg` (`rank`, `nullspace`, `rref`), `modules.linalg_mod`
  (`matmul`, `mat_rank`, `cols_to_matrix`, `identity`, `col`), P39
  `modules.complexes.ChainComplex` (`degrees`, `term(n).dim`,
  `differential(n).matrix`).
- Produces:
  ```python
  class FilteredComplex:
      """A bounded homological Domain vector-space complex with an increasing,
      exhaustive, Hausdorff filtration given per degree by column-span bases.

      terms  : {n: int}           dim V_n (missing degree = 0)
      dmats  : {n: matrix}        d_n : V_n -> V_{n-1}, rows=target (homological)
      filt   : {n: list[cols]}    filt[n][j] = a column-basis (list of coordinate
                                  columns over V_n's basis) for F_{lo+j} V_n, an
                                  INCREASING chain: F_lo <= F_{lo+1} <= ... <= V_n.
                                  Level index j runs 0..len-1; the filtration
                                  degree p = lo + j (lo may be negative).
      lo     : int                the least filtration degree present.

      With check=True: (1) every d_n is a genuine chain differential
      (d_{n-1} . d_n == 0, loud); (2) each filt[n] is a nested increasing chain
      whose top level spans all of V_n (exhaustive) -- QuiverlabError otherwise;
      (3) d_n(F_p V_n) <= F_p V_{n-1} for every p (the filtration is a subcomplex
      filtration) -- QuiverlabError naming (n, p) otherwise.
      """
      def __init__(self, terms, dmats, filt, lo, dom, check=True)
      @classmethod
      def from_chain_complex(cls, X, filt, check=True) -> "FilteredComplex"
          # forget X (a P39 ChainComplex) to (dims, dmats) via term(n).dim +
          # differential(n).matrix; attach the caller's per-degree filtration.
      def degrees(self) -> list[int]
      def dim(self, n) -> int
      def dmat(self, n) -> list            # zeros(dim V_{n-1}, dim V_n) if absent
      def levels(self) -> list[int]        # filtration degrees p present, ascending
      def piece(self, n, p) -> list        # column-basis of F_p V_n (clamped: below
                                           # lo -> [], at/above top -> full basis)
      def total_homology_dims(self) -> dict # {n: dim H_n(V_.)} via the rank formula
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/specseq/test_filtered.py
"""FilteredComplex: validation (closed under d, exhaustive) + total homology.
Self-certifying (d.d=0 gate, filtration-subcomplex gate, rank identity)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import ChainComplex
from quiverlab.specseq.filtered import FilteredComplex

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def _trivial_filt(terms):
    # the one-step filtration F_0 V_n = V_n (a valid, exhaustive, degenerate filt)
    return {n: [[[1 if i == j else 0 for i in range(d)] for j in range(d)]]
            for n, d in terms.items()}, 0


def test_forget_chaincomplex_total_homology_matches():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    dims = {n: X.term(n).dim for n in X.degrees()}
    filt, lo = _trivial_filt(dims)
    F = FilteredComplex.from_chain_complex(X, {n: filt[n] for n in dims})
    # the forgotten complex has the SAME homology as X (module structure ignored)
    assert F.total_homology_dims() == X.homology_dims()


def test_non_closed_filtration_is_refused():
    # a filtration piece F_0 V_1 whose image under d escapes F_0 V_0 must raise.
    A = _a3()
    P1, S1 = A.projective(1), A.simple(1)
    from quiverlab.modules.resolution import projective_cover
    Q0, d0, _ = projective_cover(S1)              # d0: Q0 -> S1, surjective
    terms = {1: Q0.dim, 0: S1.dim}
    dmats = {1: d0}
    # F_0 V_1 = a 1-dim piece mapping to a nonzero vector, but F_0 V_0 = {0}:
    bad = {1: [[[1] + [0] * (Q0.dim - 1)]], 0: [[]]}   # F_0 V_0 empty -> not closed
    with pytest.raises(QuiverlabError, match="closed|subcomplex|F_"):
        FilteredComplex(terms, dmats, bad, lo=0, dom=A.domain)


def test_non_exhaustive_filtration_is_refused():
    A = _a3()
    terms = {0: 2}
    dmats = {}
    filt = {0: [[[1, 0]]]}                          # top level spans only a line
    with pytest.raises(QuiverlabError, match="exhaustive|span"):
        FilteredComplex(terms, dmats, filt, lo=0, dom=A.domain)


def test_dd_nonzero_is_refused():
    A = _a3()
    # d_1 = d_2 = identity on a 1-dim space => d.d != 0
    terms = {2: 1, 1: 1, 0: 1}
    dmats = {2: [[1]], 1: [[1]]}
    filt, lo = _trivial_filt(terms)
    with pytest.raises(QuiverlabError, match="d.*d|differential"):
        FilteredComplex(terms, dmats, filt, lo, A.domain)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.specseq.filtered`

- [ ] **Step 3: Implement** `filtered.py` per the Interfaces block. Notes:
  - `from_chain_complex`: `dims = {n: X.term(n).dim for n in X.degrees()}`;
    `dmats = {n: X.differential(n).matrix for n in X.degrees() if X.term(n-1).dim and X.term(n).dim}`
    (skip zero-target/source maps; `differential` already returns a zero matrix
    when absent).
  - `_validate`: (a) `d.d=0` via `lm.matmul(dmats[n-1], dmats[n], dom)` all-zero;
    (b) exhaustive = `mat_rank(cols_to_matrix(piece(n, top)), dom) == dim(n)`;
    (c) nested = each level's span contains the previous
    (`rank([prev | cur]) == rank(cur)`); (d) closed = for each `n, p`, every
    column `c` of `piece(n, p)` has `d_n · c ∈ span(piece(n-1, p))`, i.e.
    `solve_columns(cols_to_matrix(piece(n-1, p)), cols_to_matrix([d_n·c]), dom)`
    is not `None` — loud naming `(n, p)`.
  - `piece(n, p)`: clamp `p` below `lo` to `[]`, at/above the top level to the
    full identity basis; else return `filt[n][p - lo]`.
  - `total_homology_dims`: `dim H_n = dim V_n - rank d_n - rank d_{n+1}` (the
    P39 / `ext_dims` rank formula).

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/specseq/__init__.py src/quiverlab/specseq/filtered.py tests/specseq/test_filtered.py
git commit -m "feat(specseq): FilteredComplex -- bounded filtered vector-space complex, subcomplex-filtration + exhaustiveness gates"
```

---

### Task 2: `DoubleComplex` — total complex + row/column filtrations

**Files:**
- Create: `src/quiverlab/specseq/double.py`
- Test: `tests/specseq/test_double.py`

**Interfaces:**
- Consumes: `FilteredComplex` (Task 1), `modules.linalg_mod`, the block-assembly
  pattern of `engine/cyclic.py::_total_differential` (L91-127: a block lands at
  its target's row-offset in `Tot_{n-1}`) and its generic mirror
  `hochschild/cyclic.py::cyclic_homology_dims` (L82-113).
- Produces:
  ```python
  class DoubleComplex:
      """A bounded HOMOLOGICAL first-quadrant-style double complex over a Domain.

      terms : {(p, q): int}           dim D_{p,q} (missing = 0)
      d_h   : {(p, q): matrix}        D_{p,q} -> D_{p-1,q}, rows=target
      d_v   : {(p, q): matrix}        D_{p,q} -> D_{p,q-1}, rows=target

      With check=True: d_h.d_h == 0, d_v.d_v == 0, and the ANTICOMMUTATION
      d_h . d_v + d_v . d_h == 0 at every (p, q) -- QuiverlabError naming (p, q)
      otherwise. (This is the sign convention that makes Tot's differential square
      to zero; a strictly commuting bicomplex must be sign-adjusted by the caller
      -- the standard s_{p,q} = (-1)^p on one differential -- BEFORE construction.
      The Hochschild preset's (b, B) already anticommute (mixed-complex identity);
      the Hom double complex is built with the Koszul sign in Task 7.)
      """
      def __init__(self, terms, d_h, d_v, dom, check=True)
      def total(self) -> tuple          # (terms:{n:int}, dmats:{n:matrix}, dom):
          # Tot_n = (+)_{p+q=n} D_{p,q}; D = d_h + d_v : Tot_n -> Tot_{n-1};
          # each block placed at its target's row-offset (the engine/cyclic layout).
      def column_filtration(self) -> FilteredComplex   # F_j Tot_n = (+)_{p<=j} D_{p,n-p}
      def row_filtration(self) -> FilteredComplex       # F_j Tot_n = (+)_{q<=j} D_{n-q,q}
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/specseq/test_double.py
"""DoubleComplex: anticommutation gate, total-complex assembly, and both
filtrations. Self-certifying (Tot d.d=0 by construction; filtration is a
subcomplex filtration -- FilteredComplex re-checks it)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.specseq.double import DoubleComplex
from quiverlab.specseq.filtered import FilteredComplex

pytestmark = pytest.mark.oracle_selfcert


def _dom():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(7)).domain


def test_anticommute_gate_refuses_commuting_square():
    dom = _dom()
    # D_{1,1}=D_{0,1}=D_{1,0}=D_{0,0}=1-dim; d_h=d_v=1 everywhere => they COMMUTE
    # (d_h d_v = d_v d_h = 1), so d_h d_v + d_v d_h = 2 != 0 over GF(7): refused.
    terms = {(1, 1): 1, (0, 1): 1, (1, 0): 1, (0, 0): 1}
    d_h = {(1, 1): [[1]], (1, 0): [[1]]}
    d_v = {(1, 1): [[1]], (0, 1): [[1]]}
    with pytest.raises(QuiverlabError, match="anticommut|d_h.*d_v"):
        DoubleComplex(terms, d_h, d_v, dom)


def test_signed_square_totals_correctly():
    dom = _dom()
    # put the (-1)^p sign on d_v so the square anticommutes; Tot of the 2x2 square
    # [P (+) P -> P] has the homology of an exact square (acyclic in the middle).
    neg = dom.neg(dom.one())
    terms = {(1, 1): 1, (0, 1): 1, (1, 0): 1, (0, 0): 1}
    d_h = {(1, 1): [[1]], (1, 0): [[1]]}
    d_v = {(1, 1): [[neg]], (0, 1): [[1]]}     # sign on the p=1 column
    dc = DoubleComplex(terms, d_h, d_v, dom)
    terms_tot, dmats_tot, _ = dc.total()
    assert terms_tot == {2: 1, 1: 2, 0: 1}     # Tot_2=D11, Tot_1=D01(+)D10, Tot_0=D00
    # the total of a commuting (sign-fixed) 2x2 square with all maps iso is acyclic
    F = dc.column_filtration()
    assert F.total_homology_dims() == {2: 0, 1: 0, 0: 0}


def test_both_filtrations_are_valid_subcomplex_filtrations():
    # column_filtration / row_filtration must both pass FilteredComplex's closed +
    # exhaustive gates (they are built to; this pins that they do).
    dom = _dom()
    neg = dom.neg(dom.one())
    terms = {(1, 1): 1, (0, 1): 1, (1, 0): 1, (0, 0): 1}
    d_h = {(1, 1): [[1]], (1, 0): [[1]]}
    d_v = {(1, 1): [[neg]], (0, 1): [[1]]}
    dc = DoubleComplex(terms, d_h, d_v, dom)
    assert isinstance(dc.column_filtration(), FilteredComplex)   # constructs (checks)
    assert isinstance(dc.row_filtration(), FilteredComplex)
```

- [ ] **Step 2: Run to verify failure** — `ImportError: DoubleComplex`

- [ ] **Step 3: Implement** `double.py`. Notes:
  - `total()`: order the columns of `Tot_n` by descending `p` (mirror
    `engine/cyclic._tot_degrees`); `col_off[p]`, `row_off` over `Tot_{n-1}`;
    `d_h[(p,q)]` block into row-block `(p-1)` of `Tot_{n-1}`, `d_v[(p,q)]` into
    row-block `p` (its `q` drops). Reuse the exact offset-and-place loop of
    `hochschild/cyclic.py:99-112`.
  - `column_filtration()`: `F_j Tot_n` = the span of the columns whose source
    block has `p <= j`; expressed as identity columns at those offsets. `row`
    symmetric on `q`. Build `filt[n]` as the nested identity-column chain and hand
    it to `FilteredComplex(check=True)` — which re-validates closed+exhaustive
    (the double-complex differential respects both filtrations by construction, so
    a failure here is a real bug and should surface loudly).
  - `_validate`: `d_h∘d_h`, `d_v∘d_v` zero, and per `(p,q)`
    `matmul(d_h[(p-1,q)], d_v[(p,q)]) + matmul(d_v[(p-1,q)], d_h[(p,q)])` all-zero
    (guard absent blocks as zero).

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/specseq/double.py tests/specseq/test_double.py
git commit -m "feat(specseq): DoubleComplex -- anticommutation gate, total complex (engine/cyclic layout), row/column filtrations"
```

---

### Task 3: `SpectralSequence` + `Page` — the Weibel §5.4 page engine

**Files:**
- Create: `src/quiverlab/specseq/pages.py`
- Test: `tests/specseq/test_pages.py`

**Interfaces:**
- Consumes: `FilteredComplex` (Task 1), `fields.linalg`
  (`rref`/`rank`/`nullspace`/`solve`/`reduce_mod_nullspace`), `modules.linalg_mod`
  (`matmul`/`cols_to_matrix`/`col`), the canonical-rep greedy filter of
  `hochschild/cyclic.py::_capture_generic_reps` (L127-165 — the tested
  independent-modulo pick over `fields.linalg`).
- Produces the SS-page API mandated by the P36 design note
  (`docs/plans/2026-08-05-plan-36-m2-design-notes.md` §3, quoted verbatim):
  > `E^r` is a page object; `E^r_{p,q}` returns the actual module at position
  > `(p,q)` ... `page(r)[p, q]` returns the subquotient with its chosen basis,
  > `page(r).dim(p,q)` the number ... `spots` enumerates the non-zero positions ...
  > `netPage` prints the page as a 2-D grid (p across, q up) with the module
  > dimensions in the cells.
  ```python
  class Subquotient:
      p: int; q: int; dim: int
      reps: list          # canonical coordinate columns over Tot_{p+q}'s basis
                          # (reduce_mod_nullspace / greedy-independent-modulo pins
                          # them byte-reproducibly -- the CS canonicalization rule)

  class Page:
      r: int
      def __getitem__(self, pq) -> Subquotient      # page(r)[p, q]
      def dim(self, p, q) -> int
      @property
      def spots(self) -> list                        # sorted (p, q) with dim > 0
      def differential(self, p, q) -> list           # the d_r matrix
                                                     # E^r_{p,q} -> E^r_{p-r,q+r-1}
      def grid(self) -> str                          # netPage-style: p->cols, q->rows,
                                                     # dims in cells, "." for empty

  class SpectralSequence:
      def __init__(self, F: FilteredComplex)
      def page(self, r) -> Page                       # r >= 0 (E_0), memoized
      @property
      def total_homology_dims(self) -> dict           # {n: dim H_n(Tot)} (from F)
      @property
      def width(self) -> int                          # #distinct filtration degrees p
      @property
      def height(self) -> int                         # #distinct complementary q
      # convergence.ConvergenceReport is attached by Task 4 and RUN at __init__.
  ```

**Weibel §5.4.6 formulas (homological, increasing filtration `F_p`), pinned
verbatim in the module docstring** — total degree `n = p + q`, all subspaces held
as column-span coordinate matrices over `Tot_n`'s basis:

```
Z^r_{p,q} = { x in F_p Tot_{p+q} : d(x) in F_{p-r} Tot_{p+q-1} }
B^r_{p,q} = F_p Tot_{p+q}  ∩  d( F_{p+r-1} Tot_{p+q+1} )
E^r_{p,q} = Z^r_{p,q} / ( Z^{r-1}_{p-1,q+1} + B^{r-1}_{p,q} )
d^r : E^r_{p,q} -> E^r_{p-r, q+r-1},  induced by d.
E^∞_{p,q} = F_p H_{p+q}(Tot) / F_{p-1} H_{p+q}(Tot)   (associated graded of homology)
```

- [ ] **Step 1: Write the failing tests**

```python
# tests/specseq/test_pages.py
"""Spectral-sequence pages (Weibel 5.4). Self-certifying: d_r . d_r = 0 at every
page, E_{r+1} = H(E_r, d_r) (rank identity), and canonical reps are stable across
runs. A trivial (one-step) filtration collapses at E_1 = E_inf."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.complexes import ChainComplex
from quiverlab.specseq.filtered import FilteredComplex
from quiverlab.specseq.pages import SpectralSequence

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def _trivial(X):
    dims = {n: X.term(n).dim for n in X.degrees()}
    filt = {n: [[[1 if i == j else 0 for i in range(d)] for j in range(d)]]
            for n, d in dims.items()}
    return FilteredComplex.from_chain_complex(X, filt)


def test_trivial_filtration_collapses_at_E1():
    # one-step filtration => E^0 = Tot, d^0 = d, E^1 = H(Tot) = E^inf immediately.
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    ss = SpectralSequence(_trivial(X))
    E1 = ss.page(1)
    # E^1 total per degree == H_n(Tot); on a single filtration level q=0, p=n
    for n, h in X.homology_dims().items():
        assert E1.dim(n, 0) == h
    # already E_inf: d^1 is zero, page 2 == page 1
    assert ss.page(2).dim(0, 0) == E1.dim(0, 0)


def test_dr_squares_to_zero_and_E_next_is_homology():
    # a genuinely 2-column filtration (built from a double complex's total in
    # test_double; here reuse the stupid filtration of a length-2 complex split by
    # degree) -- pin d_r.d_r = 0 and dim E_{r+1} = dim ker d_r - rank d_r ... rank.
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(2), length=3)
    # split the filtration by homological degree: F_p Tot_n = Tot_n if n <= p else 0
    dims = {n: X.term(n).dim for n in X.degrees()}
    lo = min(dims)
    filt = {}
    for n, d in dims.items():
        full = [[1 if i == j else 0 for i in range(d)] for j in range(d)]
        # level j (p = lo+j): full basis once p >= n, else empty
        filt[n] = [([] if (lo + j) < n else full) for j in range(len(dims))]
    F = FilteredComplex.from_chain_complex(X, filt)
    ss = SpectralSequence(F)
    for r in (0, 1, 2, 3):
        pg = ss.page(r)
        for (p, q) in pg.spots:
            dr = pg.differential(p, q)
            # d_r . d_r = 0 : compose with the differential out of the target cell
            tgt = (p - r, q + r - 1)
            dr2 = ss.page(r).differential(*tgt)
            from quiverlab.modules import linalg_mod as lm
            comp = lm.matmul(dr2, dr, A.domain) if (dr and dr2 and dr[0] and dr2[0]) else []
            assert not any(not A.domain.is_zero(x) for row in comp for x in row)


def test_canonical_reps_are_reproducible():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=3)
    r1 = SpectralSequence(_trivial(X)).page(1)
    r2 = SpectralSequence(_trivial(X)).page(1)
    for (p, q) in r1.spots:
        assert r1[p, q].reps == r2[p, q].reps       # byte-identical reps
```

- [ ] **Step 2: Run to verify failure** — `ImportError: SpectralSequence`

- [ ] **Step 3: Implement** `pages.py`. Notes (real linear algebra over
  `fields.linalg`; write small subspace helpers and cite Weibel 5.4.6 verbatim):
  - Subspace helpers on column-span matrices `U`, `W` over `Tot_n`:
    `_dim(U) = rank(cols_to_matrix(U))`; `_sum(U, W) = U + W` (list concat) with
    `_dim_sum = rank([U|W])`; `_dim_inter = _dim(U) + _dim(W) - _dim_sum`;
    `_intersect(U, W)` = `{ Σ c_i u_i : (c;-c') ∈ null([U|W]) }` (combine the
    `U`-part of each nullspace vector).
  - `Z^r_{p,q}`: `{x ∈ span F_p Tot_n : d_n·x ∈ span F_{p-r} Tot_{n-1}}` — stack
    `x = F_p·y` and solve `d_n·F_p·y ∈ span(F_{p-r})` via
    `solve_columns` / a nullspace of the combined system; return the `x`-basis.
  - `B^r_{p,q}`: `_intersect(piece(n, p), image(d_{n+1}, piece(n+1, p+r-1)))`
    where `image(dmat, cols) = matmul(dmat, cols_to_matrix(cols))` columns.
  - `E^r` cell **dim**: `_dim(Zr) - _dim_sum(Z^{r-1}_{p-1,q+1}, B^{r-1}_{p,q})`
    (both subspaces sit inside `Z^r`, so the numerator/denominator dims give the
    quotient dim). Guard `r=0`: `Z^0_{p,q} = F_p Tot_n`, `B^0 = 0`,
    `Z^{-1} := Z^0` shifted / the whole `F_{p-1}` piece — pin the `r=0`/`r=1` base
    cases against `test_trivial_filtration_collapses_at_E1` (the arbiter).
  - **Canonical reps** (`Subquotient.reps`): mirror
    `hochschild/cyclic._capture_generic_reps` exactly — greedily pick columns of
    `Zr` that are independent *modulo* the denominator span (rank strictly
    increases when appended to the denominator basis). `rref` is deterministic so
    the pick is byte-reproducible; this is the CS canonicalization mandate applied
    to pages. `test_canonical_reps_are_reproducible` is the gate.
  - `d^r` matrix (`Page.differential`): for each rep `x` of `E^r_{p,q}`, `d_n·x`
    lands in `F_{p-r} Tot_{n-1}` (that is what `Z^r` guarantees); express it in the
    target cell `E^r_{p-r,q+r-1}`'s rep basis by `solve` modulo the target's
    denominator, canonicalized with `reduce_mod_nullspace`. The self-cert
    `d_r∘d_r=0` and `E_{r+1}=H(E_r,d_r)` are the arbiters: if a cell dim or the
    induced-map coordinates are mis-indexed, one of them fails — fix the indexing
    once, never special-case a degree (the P39 sign-arbiter discipline).
  - `grid()`: format like M2 `netPage` — `p` across (columns), `q` up (rows),
    `dim` in each cell, `.` for empty; wrap in a fenced monospace block. Copy the
    layout from the P36 note §3.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/specseq/pages.py tests/specseq/test_pages.py
git commit -m "feat(specseq): SpectralSequence + Page engine -- Weibel 5.4 subquotient pages, d_r, canonical reps, netPage grid"
```

---

### Task 4: `ConvergenceReport` — the standing self-certificate

**Files:**
- Create: `src/quiverlab/specseq/convergence.py`
- Modify: `src/quiverlab/specseq/pages.py` (run the certificate at
  `SpectralSequence.__init__`; expose `.convergence`)
- Test: `tests/specseq/test_convergence.py`

**Interfaces:**
- Produces:
  ```python
  @dataclass
  class ConvergenceReport:
      e_infinity_page: int      # max(width, height) + 1 -- E_r stabilizes by here
                                # for a bounded filtration of a bounded complex
      degenerates_at: int | None  # least r with E_r == E_{r+1} (all d_r zero),
                                  # decided by rank/dimension comparison; None if it
                                  # only stabilizes at e_infinity_page
      abutment: dict            # {n: dim H_n(Tot)} (the target of convergence)
      def collapse(self) -> bool          # degenerates_at in (1, 2)
      def prose(self) -> str              # human sentence for the worked-steps report
  def certify_convergence(ss) -> ConvergenceReport
      # THE self-certificate: assert  sum_{p+q=n} dim E_inf^{p,q} == abutment[n]
      # for every n (E_inf = the page at e_infinity_page). Raises QuiverlabError
      # "spectral sequence does not converge to its abutment" on any mismatch --
      # this can only fire on a page/filtration bookkeeping bug.
  ```
  `certify_convergence` is called at the end of `SpectralSequence.__init__` and
  the report stored as `ss.convergence`; every downstream construction (the four
  presets) therefore ships pre-certified.

- [ ] **Step 1: Write the failing tests**

```python
# tests/specseq/test_convergence.py
"""Convergence: E_inf totals == total homology (the standing self-certificate),
E_inf page bound, degeneration decided by rank. A one-step filtration degenerates
at E_1."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import ChainComplex
from quiverlab.specseq.filtered import FilteredComplex
from quiverlab.specseq.pages import SpectralSequence

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def _trivial(X):
    dims = {n: X.term(n).dim for n in X.degrees()}
    filt = {n: [[[1 if i == j else 0 for i in range(d)] for j in range(d)]]
            for n, d in dims.items()}
    return FilteredComplex.from_chain_complex(X, filt)


def test_einf_totals_equal_total_homology():
    A = _a3()
    for v in (1, 2, 3):
        X = ChainComplex.from_projective_resolution(A.simple(v), length=4)
        ss = SpectralSequence(_trivial(X))            # __init__ certifies
        rep = ss.convergence
        totals = {}
        Einf = ss.page(rep.e_infinity_page)
        for (p, q) in Einf.spots:
            totals[p + q] = totals.get(p + q, 0) + Einf.dim(p, q)
        for n, h in ss.total_homology_dims.items():
            assert totals.get(n, 0) == h
        assert rep.abutment == ss.total_homology_dims


def test_trivial_filtration_degenerates_at_E1():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=3)
    ss = SpectralSequence(_trivial(X))
    assert ss.convergence.degenerates_at == 1
    assert ss.convergence.collapse() is True


def test_bookkeeping_bug_is_caught():
    # feed a FilteredComplex whose declared filtration is a valid subcomplex filt
    # but hand-mangle the abutment expectation: the certificate compares against
    # the complex's OWN homology, so a genuine construction can never mismatch.
    # This test instead asserts the certificate is WIRED (runs at __init__) by
    # confirming ss.convergence exists and is a ConvergenceReport.
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(2), length=2)
    ss = SpectralSequence(_trivial(X))
    from quiverlab.specseq.convergence import ConvergenceReport
    assert isinstance(ss.convergence, ConvergenceReport)
```

- [ ] **Step 2: Run to verify failure** — `ImportError: certify_convergence`

- [ ] **Step 3: Implement.** `e_infinity_page = max(ss.width, ss.height) + 1`
  (metaplan §2 SS brief: bounded filtration of a bounded complex reaches `E_∞` by
  `max(width, height)+1`). `degenerates_at`: the least `r >= 1` with every `d_r`
  zero (dimension-stable page: `ss.page(r)` cell dims == `ss.page(r+1)` cell dims
  across `spots`) — the metaplan's "degeneration decidable by rank". `collapse` =
  `degenerates_at in (1, 2)`. `certify_convergence`: build the `E_∞` page totals
  per `n`, assert equality with `ss.total_homology_dims`, raise loud on any `n`.
  Call it at the tail of `SpectralSequence.__init__`.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/specseq/convergence.py src/quiverlab/specseq/pages.py tests/specseq/test_convergence.py
git commit -m "feat(specseq): ConvergenceReport + standing self-certificate (E_inf totals == total homology) run at every construction"
```

---

### Task 5: Preset 4 — Hochschild `(b,B)` bicomplex

**Files:**
- Create: `src/quiverlab/specseq/presets.py`
- Modify: `src/quiverlab/core/algebra.py` (thin `Algebra.hochschild_bB_ss` wrapper)
- Test: `tests/specseq/test_preset_hochschild.py`

**Interfaces:**
- Consumes: `hochschild/cyclic.py::connes_B_matrix`/`boundary_matrix` over any
  Domain, `A.unit_adapted()`, `DoubleComplex`/`SpectralSequence`,
  `A.cyclic_homology(top)` (the cross-engine oracle).
- Produces:
  ```python
  def hochschild_bB_ss(A, top, max_cells=4_000_000) -> SpectralSequence
      # the first-quadrant (b, B) bicomplex: D_{p,q} = C_{q-p} (the cyclic
      # bicomplex, q >= p >= 0), d_v = b (: C_k -> C_{k-1}, lowers q), d_h = B
      # (: C_k -> C_{k+1}, ... reindexed to lower p). Both filtrations available;
      # its total complex is engine/cyclic's Tot, so E_inf total == HC.
      # max_cells guards the exponential bar basis (reuse the DepthLimitError guard
      # of hochschild/cyclic.py -- do NOT re-materialize a basis to count it).
  ```
  `Algebra.hochschild_bB_ss(self, top, ...)` delegates.

**Index map (pin verbatim; arbitrate against the HC oracle).** The `(b,B)`
bicomplex `BC` has `BC_{p,q} = C_{q-p}` for `q >= p >= 0`, vertical `b` and
horizontal `B`; `Tot_n = (+)_{p} C_{n-2p}` computes `HC_n` — exactly
`engine/cyclic._tot_degrees(n) = range(n, -1, -2)` and `_total_differential`. Set
the `DoubleComplex` entries and both differentials so that `dc.total()` reproduces
`cyclic_homology_dims`'s total complex; the oracle
`E_∞ total == A.cyclic_homology(top)` is the arbiter of the exact `(p,q)` layout
and the `b`/`B` sign that makes them anticommute (the mixed-complex identity
`bB + Bb = 0` already holds on the matrices — `hochschild/cyclic.py` docstring).

- [ ] **Step 1: Write the failing tests**

```python
# tests/specseq/test_preset_hochschild.py
"""Hochschild (b,B) spectral sequence. Cross-engine: E_inf total == HC dims
(the SBI/free-oracle pins), incl. the ground-field HC = [1,0,1,0,...] closed
form. Self-cert convergence rides along (run at construction)."""
import pytest

from quiverlab import GF, Quiver, truncated_polynomial
from quiverlab.specseq.presets import hochschild_bB_ss

pytestmark = pytest.mark.oracle_crossengine


def _ground(p):
    # k as a quiver algebra: one vertex, no arrows => A = k.
    return Quiver([1], {}).algebra(relations=[], field=GF(p))


def test_ground_field_hc_closed_form():
    A = _ground(7)
    top = 6
    ss = hochschild_bB_ss(A, top)
    Einf = ss.page(ss.convergence.e_infinity_page)
    totals = {}
    for (p, q) in Einf.spots:
        totals[p + q] = totals.get(p + q, 0) + Einf.dim(p, q)
    hc = A.cyclic_homology(top)
    assert [totals.get(n, 0) for n in range(top + 1)] == list(hc.dims)
    assert list(hc.dims) == [1, 0, 1, 0, 1, 0, 1]        # HC_*(k) closed form


@pytest.mark.oracle_literature
def test_dual_numbers_hc():
    # k[x]/(x^2): the (b,B) SS abuts to HC_*(k[x]/(x^2)); pin E_inf total == HC.
    A = truncated_polynomial(2, field=GF(5))
    top = 5
    ss = hochschild_bB_ss(A, top)
    Einf = ss.page(ss.convergence.e_infinity_page)
    totals = {}
    for (p, q) in Einf.spots:
        totals[p + q] = totals.get(p + q, 0) + Einf.dim(p, q)
    hc = A.cyclic_homology(top)
    assert [totals.get(n, 0) for n in range(top + 1)] == list(hc.dims)


def test_max_cells_guard_is_loud():
    from quiverlab.errors import DepthLimitError
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}).algebra(
        relations=[], field=GF(5))            # multi-vertex: bar basis blows up
    with pytest.raises(DepthLimitError):
        hochschild_bB_ss(A, 6, max_cells=10_000)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.specseq.presets`

- [ ] **Step 3: Implement** `hochschild_bB_ss` in `presets.py`, reusing
  `hochschild/cyclic.py`'s `boundary_matrix`/`connes_B_matrix` and its cell-count
  guard (raise `DepthLimitError` BEFORE materializing any basis, via the same
  `m*(m-1)^k` length arithmetic — do not enumerate). Build the `DoubleComplex`
  with `D_{p,q} = C_{q-p}`, place `b` and `B`; certify by construction (the
  `SpectralSequence.__init__` self-cert) and let the HC cross-engine test arbitrate
  the layout. Thin wrapper in `core/algebra.py`.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/specseq/presets.py src/quiverlab/core/algebra.py tests/specseq/test_preset_hochschild.py
git commit -m "feat(specseq): preset -- Hochschild (b,B) bicomplex SS, E_inf total == HC (ground-field [1,0,1,0,...] pinned)"
```

---

### Task 6: Preset 3 — associated-graded / radical filtration

**Files:**
- Modify: `src/quiverlab/specseq/presets.py`
- Test: `tests/specseq/test_preset_radical.py`

**Interfaces:**
- Consumes: P39 `ChainComplex` (input `X`), `Module.action` (arrow matrices) for
  the radical subspaces, `FilteredComplex.from_chain_complex`, `SpectralSequence`.
- Produces:
  ```python
  def radical_filtration_ss(X) -> SpectralSequence
      # X : a P39 ChainComplex of A-modules. F_p X_n = X_n . rad^{max(0,-p)}:
      # F_0 = X_n (whole), F_{-1} = rad X_n, F_{-2} = rad^2 X_n, ... -- an
      # INCREASING (in p), exhaustive (F_0 = X_n), Hausdorff (rad nilpotent =>
      # F_{-large} = 0) filtration. Converges to H(X) by the standing self-cert.
  ```

**Radical subspaces (exact, no module re-instantiation).** For a basic algebra,
`rad M = sum over arrows a of image(M.action[a])`; iterate for `rad^i`. Compute
the descending chain of column-spans INSIDE each `X_n`'s basis directly:
`_rad_step(cols) = column-span of { M.action[a] . c : a in arrows, c in cols }`,
starting from the full identity basis, until it stabilizes at `{0}`. Reverse into
the increasing `filt[n]`. (Cite `modules/radtopsoc.py::radical` for the
single-step definition; the iterate stays in `X_n`'s coordinates, which
`radical()` alone would not — it returns a fresh submodule.)

- [ ] **Step 1: Write the failing tests**

```python
# tests/specseq/test_preset_radical.py
"""Radical-filtration spectral sequence. Self-cert: converges to H(X) always
(the standing certificate). Koszul tie-in is an ARBITRATED oracle hypothesis:
for a Koszul algebra the SS on the Hom complex of the minimal resolution of k
degenerates early -- pin the provable page, correct the statement if it differs."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.modules.complexes import ChainComplex
from quiverlab.specseq.presets import radical_filtration_ss

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_converges_to_homology():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    ss = radical_filtration_ss(X)                 # __init__ self-certifies
    assert ss.convergence.abutment == X.homology_dims()


def test_semisimple_input_degenerates_immediately():
    # a complex of SEMISIMPLE modules (rad = 0) has a trivial radical filtration:
    # F_0 = whole, F_{-1} = 0, so E_1 = E_inf (collapse).
    A = Quiver([1, 2], {}).algebra(relations=[], field=GF(7))   # semisimple k x k
    X = ChainComplex.stalk(A.simple(1), 0)
    ss = radical_filtration_ss(X)
    assert ss.convergence.collapse() is True


@pytest.mark.oracle_literature
def test_koszul_degeneration_arbitrated():
    # kA_n is Koszul (hereditary => the minimal resolution of a simple is linear).
    # Build X = the minimal projective resolution of S_1 as a ChainComplex, run the
    # radical-filtration SS, and PIN the degeneration page the computation proves.
    # ARBITRATE: if the precise page differs from "E_2 collapse", correct THIS
    # assertion to the provable statement and record it in the plan/verification
    # page (Design-decision 3: never force the oracle). Placeholder to sharpen:
    A = linear_path_algebra(3, field=GF(5))       # hereditary kA3, Koszul
    X = ChainComplex.from_projective_resolution(A.simple(1), length=3)
    ss = radical_filtration_ss(X)
    # the resolution is already a complex of projectives; the radical filtration's
    # degeneration page is computed and pinned here -- SHARPEN during implementation:
    assert ss.convergence.degenerates_at is not None
    assert ss.convergence.abutment == X.homology_dims()   # always true (self-cert)
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** `radical_filtration_ss`. Build the per-degree
  descending radical chain via `_rad_step`, reverse into the increasing
  `filt[n]` (level `p = -i` for `rad^i`, with `lo = -(max nilpotency)`), and hand
  to `FilteredComplex.from_chain_complex`. **Sharpen the Koszul arbiter:** during
  implementation, compute `degenerates_at` for the kA_n case and replace the
  placeholder assertion with the exact provable degeneration statement; if it is
  not the folklore "E_2 collapse", CORRECT the assertion and add a one-line
  honest-scope note (the P40 Igusa–Todorov-constant arbitration discipline). The
  self-cert (`abutment == H(X)`) always holds and is the load-bearing oracle.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/specseq/presets.py tests/specseq/test_preset_radical.py
git commit -m "feat(specseq): preset -- radical-filtration SS, converges to H(X); Koszul degeneration arbitrated"
```

---

### Task 7: Presets 1 & 2 — Grothendieck + Cartan–Eilenberg change-of-rings

**Files:**
- Modify: `src/quiverlab/specseq/presets.py`
- Test: `tests/specseq/test_preset_grothendieck.py`

**Interfaces:**
- Consumes: `modules/resolution.py::minimal_resolution` (the `B`-projective
  resolution `Q_•` of `M`), `modules/injective.py::injective_resolution` (the
  `A`-injective coresolution `J^•` of `N`), `modules/hom.py::hom_space`,
  `modules/ext.py::ext_dims` (the acyclicity probe + the cross-engine oracle),
  `DoubleComplex`/`SpectralSequence`, the `Module` re-instantiation idiom
  (`Module(A, dim, action, side=...)`) for restriction along an admissible
  quotient.
- Produces:
  ```python
  def grothendieck_double_complex(M, U, N, p_len, q_len, max_term_dim=200000) -> DoubleComplex
      # D_{p,q} (cohomological, stored with NEGATED total degree) built as
      # Hom_B(Q_p(M), Hom_A(U, J^q(N))) = Hom_A(res_A Q_p, J^q(N)) for the U=B
      # (change-of-rings) case; general U is a (B,A)-bimodule given as module data.
      # PER-INSTANCE ACYCLICITY CHECK: probe Ext_B^{>0}(M, Hom_A(U, J^q)) == 0 for
      # each injective term J^q; refuse loudly (NotImplementedError-style, naming
      # the failed hypothesis) when it does not hold -- never a wrong abutment.
  def cartan_eilenberg_ss(A, B, M, N, p_len=6, q_len=6) -> SpectralSequence
      # the U = B special case for an ADMISSIBLE QUOTIENT B = A/I' of A: A and B
      # share the quiver, rel(A) subset of rel(B); a B-module restricts to A by
      # REINTERPRETING its action over A (identical matrices; valid because I' >= I).
      # E_2^{p,q} = Ext_B^p(M, Ext_A^q(B, N))  =>  Ext_A^{p+q}(M, N)  (the abutment).
      # SCOPE (honest, Design-decision 2): B is restricted to quotient algebras
      # kQ/I -> kQ/I' with I subset I'; full-generality change-of-rings for an
      # arbitrary algebra map is out of scope this release (stated on the
      # verification page). M is a B-module; N is an A-module.
  ```

**The construction, pinned (arbitrated by the cross-engine oracle).** For the
change-of-rings case `U = B`, the adjunction `Hom_B(Q_p, Hom_A(B, J^q)) ≅
Hom_A(Q_p ⊗_B B, J^q) = Hom_A(res_A Q_p, J^q)` collapses the double complex to
`D^{p,q} = Hom_A(res_A Q_p, J^q(N))`, horizontal differential from `d^Q_p`
(pre-composition, raises `p`), vertical from the `J^•` coresolution
(post-composition, raises `q`), Koszul sign on one so they anticommute. **Two
oracles arbitrate every index/side/sign choice:**
1. **Self-cert** (always, from `SpectralSequence.__init__`): `E_∞` totals ==
   `H^{p+q}(Tot D)`.
2. **Cross-engine** (the abutment): `H^n(Tot D) == Ext_A^n(M, N) = A.ext(M, N, n)`
   — holds under the per-instance acyclicity hypothesis; the probe guards it, and a
   failing probe is a loud refusal, not a wrong number. Compare at the negated
   homological degree (cohomological source stored negated).

The **degenerate pin** `B = A` (trivial quotient `I' = I`): `Ext_A^q(A, N) = N` for
`q = 0` and `0` for `q > 0`, so `E_2^{p,0} = Ext_A^p(M, N)` and the SS collapses at
`E_2 = E_∞` — a clean self- and cross-engine anchor.

- [ ] **Step 1: Write the failing tests**

```python
# tests/specseq/test_preset_grothendieck.py
"""Cartan-Eilenberg / Grothendieck change-of-rings SS. Cross-engine: E_inf total
== module Ext on 3+ instances incl. a multi-vertex one; degenerate B=A collapse;
acyclicity-failure loud refusal."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.specseq.presets import cartan_eilenberg_ss

pytestmark = pytest.mark.oracle_crossengine


def _a3rel():                                       # A = kA3 / (a*b)
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def _kA3():                                         # hereditary kA3 (I = 0)
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=GF(5))


def _einf_totals(ss, top):
    Einf = ss.page(ss.convergence.e_infinity_page)
    totals = {}
    for (p, q) in Einf.spots:
        totals[p + q] = totals.get(p + q, 0) + Einf.dim(p, q)
    return [totals.get(n, 0) for n in range(top + 1)]


def test_degenerate_B_equals_A_collapses_to_ext():
    # B = A: E_2 = E_inf, abutment == Ext_A(M, N).
    A = _kA3()
    B = _kA3()                                      # same presentation => I' = I
    M, N = B.simple(1), A.simple(3)
    ss = cartan_eilenberg_ss(A, B, M, N, p_len=4, q_len=4)
    got = _einf_totals(ss, 3)
    assert got == [A.ext(A.simple(1), A.simple(3), n) for n in range(4)]
    assert ss.convergence.collapse() is True


def test_change_of_rings_abuts_to_A_ext():
    # A = kA3 (hereditary), B = kA3/(a*b) (an admissible quotient): a B-module M
    # and an A-module N; E_inf total == Ext_A(M|_A, N).
    A = _kA3()
    B = _a3rel()                                    # rel(A)=∅ subset of {a*b}
    M, N = B.simple(2), A.simple(1)
    ss = cartan_eilenberg_ss(A, B, M, N, p_len=5, q_len=5)
    got = _einf_totals(ss, 4)
    from quiverlab.modules.module import Module
    M_over_A = Module(A, M.dim, M.action, side=M.side)   # restriction
    assert got == [A.ext(M_over_A, N, n) for n in range(5)]


def test_multivertex_instance():
    # the required multi-vertex change-of-rings pin (Design-decision 5a).
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=[], field=GF(7))
    B = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=["a*b", "b*c"], field=GF(7))
    M, N = B.simple(1), A.simple(4)
    ss = cartan_eilenberg_ss(A, B, M, N, p_len=5, q_len=5)
    got = _einf_totals(ss, 4)
    from quiverlab.modules.module import Module
    M_over_A = Module(A, M.dim, M.action, side=M.side)
    assert got == [A.ext(M_over_A, N, n) for n in range(5)]


def test_acyclicity_failure_refuses_loudly():
    # a constructed instance whose per-instance acyclicity probe fails must raise
    # a loud NotImplementedError-style refusal naming the hypothesis, NOT return a
    # wrong abutment. (Sharpen the witness during implementation: pick M, N where
    # Ext_B^{>0}(M, Hom_A(B, J^q)) != 0 for some q.)
    from quiverlab.errors import QuiverlabError
    A = _kA3()
    B = _a3rel()
    M, N = B.simple(1), A.simple(1)
    try:
        ss = cartan_eilenberg_ss(A, B, M, N, p_len=3, q_len=3)
    except QuiverlabError as exc:
        assert "acyclic" in str(exc).lower() or "hypothesis" in str(exc).lower()
    else:
        # if THIS instance happens to satisfy the hypothesis, the abutment must
        # still be correct (cross-engine) -- never silently wrong.
        got = _einf_totals(ss, 2)
        from quiverlab.modules.module import Module
        M_over_A = Module(A, M.dim, M.action, side=M.side)
        assert got == [A.ext(M_over_A, N, n) for n in range(3)]
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement.** Build `Q_• = minimal_resolution(M, p_len)` over `B`,
  `J^• = injective_resolution(N, q_len)` over `A`; restrict each `Q_p` to `A` by
  `Module(A, Q_p.dim, Q_p.action, side=Q_p.side)` (assert `A`, `B` share the
  quiver and `rel(A) ⊆ rel(B)` first — reinterpretation is valid exactly then, and
  the reconstructed module is re-validated by the `Module` constructor). `D^{p,q}
  = hom_space(res_A Q_p, J^q_module)` as vector spaces; horizontal/vertical
  differentials by pre-/post-composition (the `ext._delta_matrix` coordinate
  pattern), Koszul sign on one. **Per-instance acyclicity probe:** for each `q`,
  compute `ext_dims(B, M, Hom_A(B, J^q), k)` for `k >= 1` and refuse loudly if any
  is nonzero. Store the cohomological complex with negated total degree; the
  cross-engine test arbitrates the exact side/sign. `cartan_eilenberg_ss` calls
  `grothendieck_double_complex` with `U = B` then wraps the chosen filtration in a
  `SpectralSequence`.
  - **Sharpen `test_acyclicity_failure_refuses_loudly`:** find a concrete `(A, B,
    M, N)` whose probe genuinely fails and pin the loud refusal (remove the
    try/else hedge into a definite `pytest.raises` once the witness is known); if
    no small witness exists in the quotient-scope, KEEP the hedge and record why
    on the verification page (honest — the refusal path is still implemented and
    unit-tested with a synthetic nonzero probe).

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/specseq/presets.py tests/specseq/test_preset_grothendieck.py
git commit -m "feat(specseq): presets -- Grothendieck/Cartan-Eilenberg change-of-rings SS (quotient-scoped), E_inf total == A.ext, per-instance acyclicity refusal"
```

**Happel-LES oracle — DEFERRED (recorded).** Design-decision 5's NOTE and the
metaplan P42 card name a Happel one-point-extension LES oracle. The CE preset is
scoped to **admissible quotient maps** `kQ/I -> kQ/I'` (Design-decision 2), which
do NOT reach one-point extensions (those change the quiver — a new vertex). The
Happel-LES oracle is therefore **DEFERRED to P44** (C7 constructions: one-point
`(co)extensions A[M]` "feeds P42's Happel-LES oracle", metaplan §5 P44), with the
pointer recorded on the verification page. `families.TrivialExtension` is NOT a
substitute (it doubles the algebra via a socle-dual quiver, not a change-of-rings
along an algebra map), so this plan does not misuse it as the CE oracle instance.

---

### Task 8: Macaulay2 crosscheck — one commutative page table

**Files:**
- Modify: `src/quiverlab/m2/scripts.py` (`commutative_ss_script`),
  `src/quiverlab/m2/crosscheck.py` (`crosscheck_commutative_ss`, registered)
- Test: `tests/m2/test_spectral_sequence_m2.py`

**Interfaces:**
- Consumes: `m2/session.py::run_script`/`should_skip_m2`, `scripts.SENTINEL`/
  `parse_sentinels`, the `CrosscheckReport(what, ours, qpa, agree).assert_agree()`
  container (`m2/crosscheck.py`), M2's `SpectralSequences` package (probed live —
  the M2 design note §3 drove it at plan time under M2 1.26). Our side: a small
  `DoubleComplex` over `GF(p)` matching the M2 filtered complex.
- Produces a live `m2`-gated crosscheck: build ONE small commutative filtered
  complex both systems can express (e.g. a two-step filtration of a short Koszul
  complex over `ZZ/p[x,y]`), print the sentinel page-cell dimensions from M2's
  `SpectralSequences`, and compare to our `SpectralSequence` page dims on the SAME
  data. **Convention arbitration (honest):** M2's page indexing / filtration
  orientation may differ from ours; the compared quantity is the **`E_∞` totals
  per total degree** (== the total-complex homology — convention-robust, both
  systems must agree) PLUS, where the page indexing is confirmed to match live,
  the `E_2` cell grid. Record exactly what was compared in the test docstring
  (the P39 "QPA `Shift` uses the opposite `-k` convention" precedent — arbitrate
  live, never a silent skip). Honest scope: M2 `SpectralSequences` is
  commutative-only, like every M2 oracle (Plan-36 scope note).

- [ ] **Step 1: Probe M2 + write the script/crosscheck + the battery**

```python
# tests/m2/test_spectral_sequence_m2.py
"""Live M2 SpectralSequences crosscheck (commutative-only, m2-gated). Compares the
convention-robust E_inf totals (== total-complex homology) of a small filtered
complex both systems build; the E_2 grid is compared only where the page indexing
is confirmed to match live (recorded here)."""
import pytest

from quiverlab.m2 import crosscheck as cc
from quiverlab.m2 import session

pytestmark = pytest.mark.skipif(session.should_skip_m2(),
                                reason="Macaulay2 not installed")


def test_commutative_ss_einf_totals_match_m2():
    # a small filtered complex over ZZ/7[x,y]; compare E_inf totals (robust).
    cc.crosscheck_commutative_ss(p=7, top=4).assert_agree()
```

- [ ] **Step 2: Run live** (`-m m2`); Expected: PASS (or an honest skip when M2 is
  absent). If M2's page objects turn out not to be scriptable through
  `--script`/sentinel stdout the way `Complexes` is, fall back to comparing the
  `E_∞` totals (== M2's `homology` of the filtered complex, which DOES script) and
  record that in the docstring — honest coverage, never a silent skip. (`m2` bucket
  = the fifth oracle class; the test FAILS loudly under `QUIVERLAB_REQUIRE_M2=1`.)

- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/m2/ tests/m2/test_spectral_sequence_m2.py
git commit -m "test(m2): SpectralSequences crosscheck -- commutative filtered-complex E_inf totals vs M2 (live, arbitrated)"
```

---

### Task 9: the `ss_hochschild` no-code compute kind (GUI + webapp + HPC)

**Files:**
- Create: `src/quiverlab/specseq/block.py` (`specseq_block`, the shared builder)
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch`: kind `ss_hochschild` after the
  `cyclic_homology` branch ~L1155-1169; `_snippet` lambda map ~L1737)
- Modify: `docs/gui/runner.py` (`compute_one` twin branch after `cyclic_homology`
  L618-634; `ETA_MODEL["scalars"]` L886-897 add `"ss_hochschild"`; the
  `python_snippet` `calls` dict L799-826 add an entry)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (byte-identical twins):
  checkbox `qlgui-ss_hochschild`(+`-top`) near the `cyclic_homology` row L70;
  `S.ids` L141-157; range push in `buildRequest` L655-656; `renderBlock` branch
  after the `cyclic_homology` renderer L2271-2288
- Modify: `webapp/server/i18n/en.json` + `es.json` (add `inv.ss_hochschild` in
  BOTH; a checkbox `<label>` in `webapp/templates/index.html` near the
  `cyclic_homology` compute row ~L49)
- Modify: `src/quiverlab/trace/results_html.py` (`_HEADINGS` L29-56 add
  `ss_hochschild`; `_block_html` branch after `cyclic_homology` L301-320)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE new fixture `ss_hochschild_dualnumbers`, changelog bullet)
- Test: `tests/webapp/test_ss_hochschild_p42.py` (cross-runner contract),
  `tests/gui/test_ss_hochschild_gui_wiring.py` (checkbox/push/i18n/byte-identity)

**Interfaces:**
- `webapp/server/schema.py` needs **no change** — `parse_compute_item`'s
  `_RANGE` regex already accepts `ss_hochschild:0..4`; `MODULE_KINDS` is untouched
  (this is an **algebra-only** kind, schema v1). `webapp/static/app.js` needs **no
  change** — algebra scalar/range blocks surface via its JSON dump
  (`renderResult` L484), exactly as `cyclic_homology`/`ext_algebra` do today.
- `specseq_block(A, top, max_cells=...)` — the SHARED builder both runners call
  (byte-identical blocks, the `recognizers_block`/`ext_algebra_block` precedent):
  ```python
  def specseq_block(A, top, max_cells=4_000_000) -> dict:
      """The no-code Hochschild (b,B) spectral-sequence block: E_inf page dims per
      (p,q), the convergence prose, and the abutment (== HC dims). Loud DepthLimit
      guard on the exponential bar basis is caught and reported as {"error": ...}
      (the Plan-30 per-entry honesty precedent), never a 500."""
      # returns {"kind": "ss_hochschild", "top": top,
      #          "einf": [[p, q, dim], ...],          # E_inf page support
      #          "grid": "<netPage text>",
      #          "abutment": [dim H_0, ...],           # == A.cyclic_homology dims
      #          "degenerates_at": int | None, "collapse": bool,
      #          "prose": "<convergence sentence>",
      #          "references": ["cyclic", "weibel_homological"], "citations": [...]}
      # OR {"kind": "ss_hochschild", "error": "<the loud DepthLimitError message>"}
  ```

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated
  dir per the Plan-32 ruling; copy the `_server`/`_pyodide` pair from
  `tests/webapp/test_koszul_exposure_p38.py`):

```python
# tests/webapp/test_ss_hochschild_p42.py
"""ss_hochschild scalar/range kind: served by hpc.spec, mirrored byte-identically
by the Pyodide twin (docs/gui/runner.py)."""
import importlib.util
import json
from pathlib import Path

from quiverlab.hpc import spec


def _req():
    # k[x]/(x^2) over GF(5), compute = ["ss_hochschild:0..4"] -- copy the exact
    # request-dict shape from test_koszul_exposure_p38's builder (schema v1).
    ...


def _server(req, tmp_path):
    return spec.run(req, tmp_path)["results"]


def _pyodide(req):
    p = Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    m = importlib.util.spec_from_file_location("gui_runner", p)
    mod = importlib.util.module_from_spec(m)
    m.loader.exec_module(mod)
    mod.run_build(json.dumps(req))
    out = {}
    for s in req["compute"]:
        out[s.split(":")[0]] = mod.compute_one(s)["block"]
    return out


def test_block_shape(tmp_path):
    req = _req()
    block = next(b for b in _server(req, tmp_path)
                 if b.get("kind") == "ss_hochschild")
    assert block["abutment"] == list(
        __import__("quiverlab").truncated_polynomial(2, field=__import__("quiverlab").GF(5))
        .cyclic_homology(4).dims)
    assert "grid" in block and "prose" in block


def test_twin_parity(tmp_path):
    req = _req()
    srv = next(b for b in _server(req, tmp_path) if b.get("kind") == "ss_hochschild")
    twin = _pyodide(req)["ss_hochschild"]
    assert json.dumps(srv, sort_keys=True) == json.dumps(twin, sort_keys=True)
```

- [ ] **Step 2: Run to verify failure** — the block is unknown to both runners.

- [ ] **Step 3: Implement the handler + wiring.** `specseq_block` in
  `specseq/block.py` (catch `DepthLimitError` into `{"error": ...}`). In
  `hpc/spec.py::_dispatch`, after the `cyclic_homology` branch, add
  `if kind == "ss_hochschild": ...` (parse `item.hi`, call `specseq_block`, attach
  `_citation_pairs`); mirror byte-for-byte in `docs/gui/runner.py::compute_one`
  after its `cyclic_homology` elif. Then the GUI wiring: checkbox +`S.ids` + push
  in both `gui.js` files, `renderBlock` branch (dims/grid table + engine note +
  the convergence prose), `ETA_MODEL["scalars"]["ss_hochschild"] = 2.0`, the
  `python_snippet` `calls` entry `"ss_hochschild": "A.hochschild_bB_ss(%d)"`, the
  `_snippet` lambda `"ss_hochschild": lambda it: f"A.hochschild_bB_ss({it.hi})"`,
  the `results_html.py` `_HEADINGS` + `_block_html` branch (mirror the
  `cyclic_homology` block: a `_dims_table("dim E_inf total", abutment)` + the grid
  in a monospace block + the prose), and `inv.ss_hochschild` in en/es + the
  `index.html` checkbox label. Copy `gui.js` → the webapp vendored twin, verified
  by `tests/webapp/test_draw_page.py::test_gui_js_and_css_are_byte_identical_to_docs_gui`.

- [ ] **Step 4: Add ONE golden fixture** `ss_hochschild_dualnumbers` to
  `_runner_goldens.json` (schema v1, `compute:["ss_hochschild:0..4"]`, k[x]/(x^2)
  over GF(5)); add the changelog bullet to `test_runner_delegation.py`'s docstring.
  Run the delegation test BEFORE adding it to confirm the existing entries stay
  byte-identical (the Plan-38 `ext_algebra_exterior_gf7` precedent).

- [ ] **Step 5: Write + run the GUI-wiring test** (mirror
  `tests/gui/test_products_gui_wiring.py`: checkbox ids present, push order,
  `inv.ss_hochschild` in BOTH locales, gui.js byte-identity), then the gates:

Run: `... -m pytest tests/webapp/test_ss_hochschild_p42.py tests/webapp/test_runner_delegation.py tests/gui/test_ss_hochschild_gui_wiring.py tests/hpc -q`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): ss_hochschild compute kind -- (b,B) spectral-sequence block, both runners byte-identical, one golden"
```

---

### Task 10: buckets, exports, citations, verification page, README, gate

**Files:**
- Modify: `tests/conftest.py` (`_DEEP_DIRS` += `"specseq"`)
- Modify: `src/quiverlab/__init__.py` (export `SpectralSequence`, `DoubleComplex`,
  `FilteredComplex`), `src/quiverlab/specseq/__init__.py` (re-export the public
  surface + the four presets)
- Modify: `src/quiverlab/citations/registry.py` (+ `references.bib`)
- Modify: `docs/verification.md`, `docs/plans/2026-08-05-metaplan-v0.2.0.md`
  (tick the P42 card's delivery note)
- Modify: `README.md`
- Test: existing release gates + one top-level-import test in `tests/specseq/`

- [ ] **Step 1: Bucket edit FIRST** (do this before Task 1 in practice — repeated
  here as the audit checkpoint). Add `"specseq"` to `_DEEP_DIRS` at
  `tests/conftest.py:50`. The partition audit `tests/release/test_markers.py` is
  self-adjusting (derives buckets from a live collect-only sweep — no hardcoded
  dir list), so no edit there; run `... -m pytest tests/release -q` and confirm
  `test_buckets_partition_the_suite` stays green with `tests/specseq/*` now under
  `-m deep`.

- [ ] **Step 2: Exports + one import test.**
  `from quiverlab import SpectralSequence, DoubleComplex, FilteredComplex` works;
  `from quiverlab.specseq import hochschild_bB_ss, radical_filtration_ss,
  cartan_eilenberg_ss, grothendieck_double_complex`; `A.hochschild_bB_ss(4)`
  round-trips. Add a `tests/specseq/test_public_surface.py` (`oracle_selfcert`).

- [ ] **Step 3: Citations.** Add (verified BibTeX, the `_r(...)` registry
  precedent `registry.py:24`, only keys that resolve): `weibel_homological`
  (Weibel, *An Introduction to Homological Algebra*, CUP 1994 — the §5.4 page
  formulas + strong-convergence theorem), `barakat_homalg` (Barakat,
  *Spectral Filtrations via Generalized Morphisms* / homalg, arXiv:0904.0240 — the
  Grothendieck-over-module-categories framing, metaplan §2). REUSE existing keys:
  `cartan_eilenberg` → add a registry alias to the present `CartanEilenberg1956`
  bibtex (the change-of-rings SS) if not already aliased; `cyclic` (the `(b,B)`
  bicomplex); `priddy`/`froberg_koszul` (the Koszul tie-in); `happel_question`
  (Happel1989 — cited ONLY in the deferred Happel-LES honest-scope note);
  `assem_book` (generic). Do NOT cite the "ISH spectral sequence" (metaplan §2:
  unverified — not standard).

- [ ] **Step 4: Verification page.** Add the Plan-42 subsystem row
  (`specseq/` | oracles: `oracle_selfcert` `d_r∘d_r=0` + `E_{r+1}=H(E_r,d_r)` +
  `E_∞` totals == total homology + radical-SS converges to `H(X)`;
  `oracle_crossengine` CE/Grothendieck `E_∞` total == module `Ext` (incl.
  multi-vertex) + Hochschild `E_∞` total == `HC`; `oracle_literature` ground-field
  `HC=[1,0,1,0,...]` + `k[x]/(x^a)` + arbitrated Koszul degeneration; `m2` the
  commutative `SpectralSequences` crosscheck). Add the **honest-scope entries**:
  (a) the CE preset is scoped to **admissible quotient maps** `kQ/I -> kQ/I'`;
  full-generality change-of-rings for an arbitrary algebra map is out of scope
  this release; (b) the **Happel-LES oracle is DEFERRED to P44** (one-point
  extensions change the quiver, unreachable by quotient maps; `TrivialExtension`
  is not a substitute — recorded pointer); (c) **Grothendieck-SS acyclicity is a
  per-instance hypothesis with a loud refusal path**, not a proven-for-all
  theorem; (d) **QPA has NO spectral-sequence surface** — no `qpa` oracle for this
  subsystem (the covering oracles are cross-engine + M2); (e) the GUI exposes only
  `ss_hochschild`; the **CE/Grothendieck/radical no-code surfaces are DEFERRED to
  P50** (they need new request fields). Recount the class table
  (`tests/release/test_oracle_classes.py` drives the numbers — run collection,
  paste, re-run to green). Tick the P42 metaplan card.

- [ ] **Step 5: README.** One features line: "spectral sequences — filtered &
  double complexes, exact `E_r` pages with canonical representatives + a
  convergence certificate, and four presets (Cartan–Eilenberg change-of-rings,
  Grothendieck, radical filtration, Hochschild `(b,B)`); the `(b,B)` SS is
  clickable via `ss_hochschild`."

- [ ] **Step 6: Full gate:**
  `... -m pytest tests/specseq -q` (deep, the new dir),
  `... -m pytest tests/modules tests/hochschild -q` (the P39/cyclic neighbours —
  byte-unchanged), `... -m pytest -q -m fast`,
  `... -m pytest tests/m2 -q -m m2` (skips without M2),
  `... -m pytest tests/release -q` — all green.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-42 oracle rows + honest scope (CE quotient-scoped, Happel-LES deferred to P44, no QPA/GUI-preset surface) + specseq bucket + exports + citations"
```

---

## Acceptance (Plan-42 definition of done)

1. `FilteredComplex`, `DoubleComplex`, `SpectralSequence`/`Page`,
   `ConvergenceReport` public and validated: filtered complexes are subcomplex
   filtrations (closed under `d`, exhaustive, loud), double complexes anticommute
   (loud), pages are the Weibel §5.4 subquotients with byte-reproducible
   (`reduce_mod_nullspace` / greedy-independent-modulo) representatives, `netPage`
   grids, and `d_r` differentials.
2. **The standing self-certificate runs at every `SpectralSequence` construction**
   (`E_∞` totals == total homology) — a live rank-identity gate, never bypassed.
   `d_r∘d_r=0` and `E_{r+1}=H(E_r,d_r)` are pinned (`oracle_selfcert`).
3. All four presets ship: Hochschild `(b,B)` (`E_∞` total == `HC`, ground-field
   `[1,0,1,0,...]` + `k[x]/(x^a)` pinned); radical filtration (converges to `H(X)`
   always; Koszul degeneration arbitrated, not forced); Cartan–Eilenberg /
   Grothendieck change-of-rings (quotient-scoped, `E_∞` total == module `Ext` on
   3+ instances incl. a multi-vertex one, per-instance acyclicity loud refusal,
   degenerate `B=A` collapse pin).
4. One `m2`-gated `SpectralSequences` crosscheck green live (or an honest skip
   without M2; convention arbitrated and recorded in the docstring). No `qpa`
   oracle claimed (QPA has no SS surface — stated).
5. `ss_hochschild` clickable end-to-end (GUI canvas → block → report) in EN+ES,
   both runners byte-identical, ONE golden with a documented change-log entry,
   schema v1; `app.js`/`schema.py` correctly untouched (algebra-only kind).
6. Honest scope on the verification page: CE quotient-scoping, Happel-LES deferred
   to P44 (`TrivialExtension` not misused as its oracle), Grothendieck acyclicity
   as a per-instance hypothesis, and the CE/Grothendieck/radical no-code GUI
   deferred to P50 — every deferral satisfied by a named successor, not a silent
   skip.
7. `tests/specseq/` in the deep bucket (`_DEEP_DIRS` edited, partition audit
   green); `docs/verification.md` recounted; README line added; P39/cyclic
   neighbours byte-unchanged; deep (new dir) + fast + m2 + release suites green.
   P43 can consume `SpectralSequence`/`Page` exactly as declared in the Interfaces
   blocks.
