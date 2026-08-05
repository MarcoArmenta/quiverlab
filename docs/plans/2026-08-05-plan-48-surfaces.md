# Plan 48: Surfaces — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the **marked-surface → triangulation → quiver-with-potential → algebra**
pipeline first-class and no-code. A representation theorist picks (or draws) a surface,
triangulates it, and reads off the associated **gentle** algebra on the canvas — a new
*input* method absent from QPA. The v1 scope is the clean, fully-certifiable case:
**unpunctured surfaces with non-empty boundary.** For those, every arc-adjacency is
clean, `Jac(Q(T), W(T))` is a *gentle* algebra (Assem–Brüstle–Charbonneau-Jodoin–
Plamondon 2010, "gentle algebras arising from surface triangulations", henceforth ABCP;
Labardini-Fragoso 2009), and so every construction we ship self-certifies THREE ways at
once: the P44 Jacobian **finiteness certificate**, the P46 **`is_gentle` recognizer**, and
the FST **arc-count** identity. **Flip ↔ Fomin–Zelevinsky matrix mutation** is certified
per instance (exact skew-symmetric combinatorics — implementable; full DWZ
right-equivalence of potentials is NOT attempted, recorded as the deferral). Punctures,
closed surfaces, and self-folded triangles are **deferred loudly** with the successor
named — never a silent wrong answer.

**Architecture:** One new package `src/quiverlab/surfaces/`, five thin exact-combinatorics
modules over primitives that already exist (P44's Jacobian machinery, P46's `is_gentle` /
`ag_invariant`, the `Quiver`/`Quiver.algebra` backbone) — **no new math engines**:

- `src/quiverlab/surfaces/marked.py` — `MarkedSurface(genus, boundary_marked=(k1,...,kb),
  punctures=p)`: the combinatorial surface datum. Loud validation (each `k_i >= 1`; the
  FST admissibility list — refuse the excluded degenerate surfaces). Euler-characteristic
  invariants and the derived **ideal-arc count** `n = 6g - 6 + 3(b + p) + Σ k_i` (interior
  arcs, NOT boundary segments — derivation in Task 1).
- `src/quiverlab/surfaces/triangulation.py` — `Triangulation(surface, triangles)`: an
  ideal triangulation as ordered triangles of *sides* (interior arcs + boundary segments).
  Arc-adjacency validation (every interior arc in exactly 2 triangles; every boundary
  segment in exactly 1; the arc-count self-cert). Canonical constructors
  `fan_triangulation(...)`, `annulus_triangulation(n, m)`, `once_punctured_torus()`.
- `src/quiverlab/surfaces/qp.py` — `quiver_of(T)` (one vertex per interior arc, one arrow
  per angle between consecutive arcs in a triangle, followed by the standard 2-cycle
  reduction; orientation ARBITRATED by the disc oracle), `potential_of(T)` (one 3-cycle
  term per **internal** triangle — a triangle all of whose sides are arcs), and
  `jacobian_of(T, field)` = P44's `JacobianAlgebra(quiver_of(T), potential_of(T), field)`.
  Self-folded / puncture handling scoped OUT of v1 with a loud refusal.
- `src/quiverlab/surfaces/flip.py` — `flip(T, arc)` (the combinatorial quadrilateral flip;
  loud on a boundary segment or a non-flippable arc) and `certify_flip_mutation(T, arc)`
  (compare `quiver_of(flip(T, arc))` against the exact Fomin–Zelevinsky matrix mutation of
  `quiver_of(T)` at that vertex — the per-instance certification of the flip↔mutation
  theorem).
- `src/quiverlab/surfaces/block.py` — `surface_block(T)`, the algebra-adjacent descriptor
  (surface invariants + the resulting quiver/relations + the P46 `is_gentle` verdict + the
  AG invariant when gentle), for the report and the preset provenance. Surfaces are an
  **input method**, so there is no new webapp *compute* kind — the produced algebra flows
  through EVERY existing compute kind (hh, resolutions, modules, …).

The GUI/no-code story is **input-side** (metagoal 1, honest): three canvas presets
generated at build time (fan-disc-A3, annulus(2,2), hexagon-with-internal-triangle) via the
`scripts/gui_build_hook.py::_preset_algebras()` mechanism, plus catalog registration
(`families/discover.py` `FamilyInfo` + `webapp/server/catalog.py` skip-set). The free-form
"draw a surface, triangulate on the canvas" flagship is **DEFERRED to a post-release
successor** — recorded on the verification page.

**Tech Stack:** pure exact combinatorics (int / tuple / frozenset) + the P44 Jacobian
route (Gröbner/monomial via `Quiver.algebra`). **No floats in `src/`** (AST-gated by
`tests/test_no_floats.py`, which scans `src/`): every surface datum is exact — arc counts,
Euler characteristics, angle-orderings, and the FZ exchange matrix are integers; the
Jacobian coefficients are exact scalars, the potential coefficients `+1` (`int`).

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P44 (constructions) and P46 (gentle/string) are MERGED to `dev` and are hard
  prerequisites.** This plan consumes:
  - `families/jacobian.py::Potential(quiver, terms)` / `cyclic_derivative(W, arrow)` /
    `JacobianAlgebra(Q, W, field=None, degree_bound=None)` and the honest
    `NotFiniteDimensionalError` propagation (P44) — exactly as their signatures read on
    `dev`;
  - `invariants/recognizers.py::is_gentle(A)` (already on `dev` from P38, reaffirmed by
    P44/P46) and `strings/ag.py::ag_invariant(A)` / `permitted_threads` /
    `forbidden_threads` (P46);
  - the presented backbone `combinat/quiver.py::Quiver` (`.vertices`, `.arrows` a
    `{name: (src, tgt)}` dict, `.source`/`.target`, `.compose_ok`) + `Quiver.algebra(
    relations, field, degree_bound)`.

  Branch `plan-48-surfaces` off `dev` **after P44 and P46 have merged** — verify with
  `python -c "import quiverlab.families.jacobian, quiverlab.strings.ag; from quiverlab.invariants.recognizers import is_gentle"`
  before starting. **No dependency on P42/P43/P45/P47/P49** — this branch merges
  independently.
- **Test buckets are auto-assigned by directory** (`tests/conftest.py`): `tests/families/`
  → **deep**; `tests/invariants/`, `tests/gui/`, `tests/webapp/` → **fast**; `tests/qpa/`
  → **qpa**. Per the metaplan surfaces decision, the main batteries live in
  **`tests/families/test_surfaces_*.py`** (deep — `quiver_of`/`jacobian_of`/`flip` build
  algebras and mutate matrices), the pure combinatorial arc-count / Euler self-cert lives
  in **`tests/invariants/test_surfaces_arccount.py`** (fast — no algebra), and the
  `IsGentleAlgebra` crosscheck in **`tests/qpa/test_surfaces_qpa.py`** (qpa). Run new tests
  by path during development; finish each task with a `-m deep` / `-m fast` / `-m qpa`
  spot-run of the touched files.
- **v1 SCOPE (binding): unpunctured surfaces with non-empty boundary only** —
  `punctures == 0` and `b == len(boundary_marked) >= 1`. This is the ABCP/LFS regime where
  (i) there are **no self-folded triangles** (a self-folded triangle is a loop enclosing a
  puncture — unpunctured ⇒ none), (ii) every arc-adjacency is clean, and (iii)
  `Jac(Q(T), W(T))` is **gentle**, hence finite-dimensional and certifiable via P46's
  `is_gentle` + P44's finiteness certificate. `MarkedSurface` itself represents ANY
  FST-admissible surface (so it can carry punctured/closed data and refuse it downstream),
  but `quiver_of` / `potential_of` / `jacobian_of` **refuse loudly** on `punctures > 0`,
  `b == 0` (closed), or a detected self-folded configuration, each `QuiverlabError` naming
  the successor (P48.1 — punctured/closed surfaces + puncture potentials + self-folded
  triangles).
- **No floats, all combinatorics exact** (AST-gated `tests/test_no_floats.py`, `src/` only).
- **Loud refusals, never silent wrong answers** (the `GlobalDimension` / `TrivialExtension`
  / `JacobianAlgebra` precedent): an inadmissible `MarkedSurface`, a `Triangulation` whose
  arc-adjacency count is wrong, a flip of a boundary segment, an out-of-v1-scope `quiver_of`
  — each raises `QuiverlabError` with a `hint=`, never returns a plausible-looking wrong
  object. A Jacobian the presented backbone cannot certify finite propagates
  `NotFiniteDimensionalError` unchanged.
- **Convention flips are ARBITRATED, not assumed** (the house cup-sign / composition-order
  / Gabriel-direction precedent): the angle→arrow **orientation** (clockwise vs
  anticlockwise around each triangle) is pinned by the **disc-A_n arbiter** (the fan
  triangulation of the `(n+3)`-gon must give the linear `A_n` quiver `1→2→…→n` exactly).
  The FZ-mutation vertex-identification (which arc-index the flipped arc inherits) is pinned
  by the flip↔mutation battery. Write both candidates, let the oracle decide, record the
  outcome in the docstring.
- **House conventions:** path composition is **left-to-right** (`a*b` = first `a` then `b`,
  `target(a)=source(b)`); a potential term is an oriented **cycle** (P44 `Potential`
  validates this). All refusals are `QuiverlabError` (or its `NotFiniteDimensionalError`
  subclass).
- Plan-32 markers: arc-count / Euler / adjacency / flip-involution certificates and the
  surface self-identities = `oracle_selfcert`; the disc-fan-`A_n` identity, the annulus-Ã
  small-case pin, the hexagon-internal-triangle `dim 6` Jacobian, the FST admissibility
  refusals = `oracle_literature`; `is_gentle(jacobian_of(T))` True across the disc/annulus/
  2-holed zoo (cross-subsystem P44+P46+P48), and flip ≡ FZ-mutation = `oracle_crossengine`;
  QPA `IsGentleAlgebra`/`IsSpecialBiserialAlgebra` on the surface algebras live in
  `tests/qpa/` (bucket = the class, never double-marked).
- **Mid-merge-train counts:** v0.2.0 lands ~15 subplans in overlapping waves, so the
  absolute suite counts drift between this plan's authoring and its merge. **Task 7 recounts
  the oracle-class table at merge time by running `tests/release/test_oracle_classes.py`**
  (paste the live numbers, never a guessed-at-authoring count) and claims only the deltas
  this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class table
  green) and adds its citations to `citations/references.bib` + `registry.py`. Conventional
  commits; green tests at every commit.

---

### Task 1: `MarkedSurface` — the surface datum, arc-count, FST admissibility

The combinatorial surface. The one genuinely mathematical decision here is the **ideal-arc
count**, derived below and pinned as a self-certificate on every downstream construction.

**Mathematical grounding (derived in this plan; FST 2008 is the binding oracle).**
Let `(S, M)` be a compact connected oriented surface of genus `g` with `b` boundary
components; `M` its marked points, split into `p` **punctures** (interior) and, on boundary
component `i`, `k_i >= 1` boundary marked points. Write `c = Σ k_i` (total boundary marked
points). An **ideal triangulation** is a maximal collection of pairwise non-crossing
**arcs** — curves between marked points, not isotopic to a boundary segment or a point —
cutting `S` into ideal triangles. **Only interior arcs count** (`n`); boundary segments are
not arcs (not mutable, no vertex). Give `S` the CW structure of the triangulation:
- `V = p + c` (all marked points are 0-cells);
- `E = n + c` (the `n` interior arcs plus the `c` boundary segments — component `i` with
  `k_i` marked points has exactly `k_i` boundary segments, summing to `c`);
- `F = t` (the ideal triangles).

Two identities pin `n`:
1. **Euler characteristic.** `V - E + F = χ(S) = 2 - 2g - b`, i.e.
   `(p + c) - (n + c) + t = 2 - 2g - b`, so `p - n + t = 2 - 2g - b`.  **(I)**
2. **Side counting.** Each triangle has 3 sides; each interior arc borders 2 triangles,
   each boundary segment borders 1: `3t = 2n + c`, so `t = (2n + c)/3`.  **(II)**

Substitute (II) into (I): `p - n + (2n + c)/3 = 2 - 2g - b`. Multiply by 3:
`3p - 3n + 2n + c = 6 - 6g - 3b`, i.e. `3p - n + c = 6 - 6g - 3b`. Solving,

> **`n = 6g - 6 + 3b + 3p + c = 6g - 6 + 3(b + p) + Σ k_i`**, and `t = (2n + c)/3`.

**Sanity checks** (pinned as literature oracles in Task 3): a **disc with `(n+3)` marked
points** (`g=0, b=1, p=0, c=n+3`) gives `6·0 - 6 + 3·1 + 0 + (n+3) = n` arcs and
`t = (2n + n + 3)/3 = n + 1` triangles — exactly the diagonals and triangles of a
triangulated `(n+3)`-gon → type `A_n`. An **annulus `C(n, m)`** (`g=0, b=2, p=0, c=n+m`)
gives `-6 + 6 + (n+m) = n + m` arcs → affine `Ã_{n+m-1}`. A **once-punctured torus**
(`g=1, b=0, p=1`) gives `6 - 6 + 3 = 3` arcs (FST-admissible, but out of v1 scope).

**FST admissibility (the exclusion list; Fomin–Shapiro–Thurston 2008, §2).** `(S, M)`
admits an ideal triangulation iff each boundary component has `>= 1` marked point AND
`(S, M)` is none of: a **sphere with 1, 2, or 3 punctures**; an **unpunctured or
once-punctured monogon** (`b=1, k=(1,)`); an **unpunctured digon** (`b=1, k=(2,)`); an
**unpunctured triangle** (`b=1, k=(3,)`). Equivalently these are exactly the cases with
`n <= 0` (monogon `n=-2`, digon `n=-1`, triangle `n=0`) plus the two `n>0` degeneracies the
list names explicitly (sphere-with-3-punctures `n=3`; once-punctured monogon `n=1`).
`MarkedSurface` refuses all of them loudly.

**Files:**
- Create: `src/quiverlab/surfaces/__init__.py`, `src/quiverlab/surfaces/marked.py`
- Test: `tests/families/test_surfaces_marked.py`,
  `tests/invariants/test_surfaces_arccount.py`

**Interfaces:**
- Consumes: `quiverlab.errors::QuiverlabError`.
- Produces:
  ```python
  @dataclass(frozen=True)
  class MarkedSurface:
      genus: int
      boundary_marked: tuple      # (k1, ..., kb), each k_i >= 1; () for a closed surface
      punctures: int = 0
      # __post_init__ validates: genus >= 0, punctures >= 0, each k_i >= 1, and the FST
      # admissibility list (loud QuiverlabError on an excluded / n<1 surface).
      @property
      def num_boundary_components(self) -> int          # b = len(boundary_marked)
      @property
      def total_marked(self) -> int                     # p + Σ k_i
      def euler_characteristic(self) -> int             # 2 - 2g - b (compact, with boundary)
      def arc_count(self) -> int                        # 6g - 6 + 3(b + p) + Σ k_i  (>= 1)
      def triangle_count(self) -> int                   # (2n + c)/3  (exact integer)
      def in_v1_scope(self) -> bool                     # punctures == 0 and b >= 1
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_surfaces_arccount.py
"""FST ideal-arc count + Euler invariants (Plan 48). Pure combinatorics -- fast bucket.
Self-cert: n = 6g-6+3(b+p)+Σk_i and t=(2n+c)/3 are integers matching the derivation;
literature: disc(n+3)->A_n arc count, annulus(n,m)->n+m, once-punctured torus->3."""
import pytest

from quiverlab.surfaces.marked import MarkedSurface

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


@lit
@pytest.mark.parametrize("marked, n", [(4, 1), (5, 2), (6, 3), (7, 4), (10, 7)])
def test_disc_arc_count_is_A_n(marked, n):
    # disc with `marked` boundary points -> A_{marked-3}: n = marked - 3 arcs.
    S = MarkedSurface(genus=0, boundary_marked=(marked,), punctures=0)
    assert S.arc_count() == n == marked - 3
    assert S.triangle_count() == n + 1          # (marked-2) triangles
    assert S.euler_characteristic() == 1        # 2 - 0 - 1


@lit
@pytest.mark.parametrize("n, m", [(1, 1), (2, 1), (2, 2), (3, 1)])
def test_annulus_arc_count_is_n_plus_m(n, m):
    S = MarkedSurface(genus=0, boundary_marked=(n, m), punctures=0)
    assert S.arc_count() == n + m
    assert S.triangle_count() == n + m          # (2(n+m)+(n+m))/3
    assert S.euler_characteristic() == 0        # 2 - 0 - 2


@lit
def test_once_punctured_torus_has_three_arcs():
    S = MarkedSurface(genus=1, boundary_marked=(), punctures=1)
    assert S.arc_count() == 3 and S.triangle_count() == 2
    assert not S.in_v1_scope()                  # closed + punctured -> out of v1 scope


@selfcert
@pytest.mark.parametrize("g, bm, p", [(0, (5,), 0), (0, (2, 2), 0), (1, (1,), 0),
                                      (0, (4, 4), 0), (2, (1,), 0)])
def test_side_counting_identity_holds(g, bm, p):
    S = MarkedSurface(genus=g, boundary_marked=bm, punctures=p)
    n, t, c = S.arc_count(), S.triangle_count(), sum(bm)
    assert 3 * t == 2 * n + c                    # identity (II)
    assert p - n + t == S.euler_characteristic() # identity (I)
```

```python
# tests/families/test_surfaces_marked.py
"""MarkedSurface validation + the FST admissibility exclusion list (Plan 48, deep).
Loud refusals: the monogon/digon/triangle (n<=0), the small punctured spheres, the
once-punctured monogon; bad marked-point counts."""
import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.marked import MarkedSurface

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("g, bm, p", [
    (0, (1,), 0),        # unpunctured monogon (n = -2)
    (0, (1,), 1),        # once-punctured monogon (n = 1 but FST-excluded)
    (0, (2,), 0),        # unpunctured digon (n = -1)
    (0, (3,), 0),        # unpunctured triangle (n = 0)
    (0, (), 1),          # sphere, 1 puncture
    (0, (), 2),          # sphere, 2 punctures
    (0, (), 3),          # sphere, 3 punctures (n = 3 but FST-excluded)
])
def test_fst_inadmissible_surfaces_refused(g, bm, p):
    with pytest.raises(QuiverlabError):
        MarkedSurface(genus=g, boundary_marked=bm, punctures=p)


@selfcert
@pytest.mark.parametrize("g, bm, p", [(-1, (4,), 0), (0, (0,), 0), (0, (2, 0), 0),
                                      (0, (4,), -1)])
def test_bad_counts_refused(g, bm, p):
    with pytest.raises(QuiverlabError):
        MarkedSurface(genus=g, boundary_marked=bm, punctures=p)


@lit
def test_admissible_surfaces_construct():
    # square (n=1), pentagon (n=2), annulus (n+m), 4-punctured sphere, once-punctured torus
    assert MarkedSurface(0, (4,), 0).arc_count() == 1
    assert MarkedSurface(0, (2, 2), 0).arc_count() == 4
    assert MarkedSurface(0, (), 4).arc_count() == 6        # 4-punctured sphere is admissible
    assert MarkedSurface(1, (), 1).arc_count() == 3
```

- [ ] **Step 2: Run to verify failure**

Run: `... -m pytest tests/invariants/test_surfaces_arccount.py tests/families/test_surfaces_marked.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.surfaces.marked`

- [ ] **Step 3: Implement**

```python
# src/quiverlab/surfaces/marked.py
"""Marked bordered surfaces (Plan 48). A MarkedSurface is the combinatorial datum
(genus, boundary marked-point counts, punctures) of a compact oriented surface with
marked points. The ideal-arc count n = 6g-6+3(b+p)+Σk_i (Fomin-Shapiro-Thurston 2008)
and the Euler invariants are derived exactly; construction refuses the FST-inadmissible
surfaces loudly. Float-free (all counts are ints)."""
from __future__ import annotations

from dataclasses import dataclass, field as _dc_field

from quiverlab.errors import QuiverlabError

# FST 2008 explicit exclusions with n > 0 (n <= 0 is refused by the arc-count check):
#   (genus, sorted-boundary-tuple, punctures)
_FST_EXCLUDED = {
    (0, (), 1), (0, (), 2), (0, (), 3),      # sphere with 1/2/3 punctures
    (0, (1,), 1),                            # once-punctured monogon
}


@dataclass(frozen=True)
class MarkedSurface:
    genus: int
    boundary_marked: tuple
    punctures: int = 0

    def __post_init__(self):
        object.__setattr__(self, "boundary_marked", tuple(self.boundary_marked))
        if self.genus < 0:
            raise QuiverlabError("MarkedSurface: genus must be >= 0",
                                 hint=f"got genus={self.genus}")
        if self.punctures < 0:
            raise QuiverlabError("MarkedSurface: punctures must be >= 0",
                                 hint=f"got punctures={self.punctures}")
        for k in self.boundary_marked:
            if not isinstance(k, int) or isinstance(k, bool) or k < 1:
                raise QuiverlabError(
                    "MarkedSurface: each boundary component needs >= 1 marked point",
                    hint=f"got boundary_marked={self.boundary_marked}")
        b, p, c = self.num_boundary_components, self.punctures, sum(self.boundary_marked)
        if b == 0 and p == 0:
            raise QuiverlabError("MarkedSurface: no marked points -- cannot triangulate",
                                 hint="a closed surface needs >= 1 puncture")
        key = (self.genus, tuple(sorted(self.boundary_marked)), self.punctures)
        n = 6 * self.genus - 6 + 3 * (b + p) + c
        if n < 1 or key in _FST_EXCLUDED:
            raise QuiverlabError(
                "MarkedSurface: not FST-admissible (no nontrivial ideal triangulation)",
                hint="excluded (FST 2008): the monogon/digon/triangle, spheres with "
                     "<=3 punctures, and the once-punctured monogon")

    @property
    def num_boundary_components(self) -> int:
        return len(self.boundary_marked)

    @property
    def total_marked(self) -> int:
        return self.punctures + sum(self.boundary_marked)

    def euler_characteristic(self) -> int:
        return 2 - 2 * self.genus - self.num_boundary_components

    def arc_count(self) -> int:
        b, p, c = self.num_boundary_components, self.punctures, sum(self.boundary_marked)
        return 6 * self.genus - 6 + 3 * (b + p) + c

    def triangle_count(self) -> int:
        n, c = self.arc_count(), sum(self.boundary_marked)
        assert (2 * n + c) % 3 == 0                       # exact by the derivation
        return (2 * n + c) // 3

    def in_v1_scope(self) -> bool:
        return self.punctures == 0 and self.num_boundary_components >= 1
```

**Adjust to reality (Task 1):**
- `frozen=True` forbids `self.x = ...`; validate/normalise in `__post_init__` via
  `object.__setattr__` (the idiom above). Keep `boundary_marked` a `tuple` so the dataclass
  is hashable (it is used as a preset key later).
- The `_FST_EXCLUDED` set covers only the `n > 0` degeneracies FST names explicitly; the
  `n < 1` check absorbs the monogon/digon/triangle. Cite `fomin_shapiro_thurston` in the
  docstring; the exclusion parametrisation IS the literature oracle.
- `bool` is an `int` subclass — the `isinstance(k, bool)` guard rejects `True`/`False`
  smuggled in as a "1"/"0" count (the `catalog._params_of` bool-before-int precedent).

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/surfaces/__init__.py src/quiverlab/surfaces/marked.py \
        tests/families/test_surfaces_marked.py tests/invariants/test_surfaces_arccount.py
git commit -m "feat(surfaces): MarkedSurface -- derived FST arc count 6g-6+3(b+p)+Σk_i + Euler invariants, FST admissibility refusals"
```

---

### Task 2: `Triangulation` — ideal triangulation, adjacency validation, constructors

**Files:**
- Create: `src/quiverlab/surfaces/triangulation.py`
- Test: `tests/families/test_surfaces_triangulation.py`

**Interfaces:**
- Consumes: `marked.py::MarkedSurface`, `quiverlab.errors::QuiverlabError`.
- Produces:
  ```python
  # A SIDE is either an interior arc id (a positive int) or a boundary segment id (a
  # string like "b0_2" -- component 0, segment 2). A TRIANGLE is an ordered 3-tuple of
  # sides in ANTICLOCKWISE cyclic order (consistent with the surface orientation) -- the
  # order is what fixes the angle/arrow convention downstream.
  Side = int | str
  Triangle = tuple           # (s0, s1, s2)

  @dataclass(frozen=True)
  class Triangulation:
      surface: MarkedSurface
      triangles: tuple        # tuple of Triangle
      # __post_init__ validates: every interior arc (int side) appears in exactly 2
      # triangle-sides; every boundary segment (str side) in exactly 1; the number of
      # distinct arcs == surface.arc_count() (the arc-count self-cert); the number of
      # triangles == surface.triangle_count(); no side appears TWICE in one triangle
      # (a self-folded configuration -- refused loudly, out of v1 scope).
      def arcs(self) -> tuple                 # sorted distinct interior-arc ids
      def boundary_segments(self) -> tuple    # sorted distinct boundary-segment ids
      def triangles_containing(self, arc) -> tuple   # the (<=2) triangles with `arc` a side

  def fan_triangulation(marked) -> Triangulation
      # the fan of the (marked)-gon from vertex 0: arcs 0-2, 0-3, ..., 0-(marked-2).
  def annulus_triangulation(n, m) -> Triangulation
      # a standard triangulation of the annulus C(n, m) (n outer, m inner marked points).
  def once_punctured_torus() -> Triangulation
      # g=1, p=1: 3 loop-arcs at the puncture. FST-admissible; OUT of v1 scope (Task 3).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_surfaces_triangulation.py
"""Ideal triangulations (Plan 48, deep). Self-cert: arc-adjacency counts + the
arc-count identity on every constructor; loud on a self-folded configuration and on a
triangulation whose arc count contradicts the surface."""
import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.marked import MarkedSurface
from quiverlab.surfaces.triangulation import (Triangulation, annulus_triangulation,
                                             fan_triangulation, once_punctured_torus)

selfcert = pytest.mark.oracle_selfcert


@selfcert
@pytest.mark.parametrize("marked", [4, 5, 6, 7, 8])
def test_fan_adjacency_and_arc_count(marked):
    T = fan_triangulation(marked)
    assert len(T.arcs()) == T.surface.arc_count() == marked - 3
    assert len(T.triangles) == T.surface.triangle_count() == marked - 2
    for a in T.arcs():                                   # every interior arc in 2 triangles
        assert len(T.triangles_containing(a)) == 2
    for seg in T.boundary_segments():                    # every boundary segment in 1
        assert sum(seg in tri for tri in T.triangles) == 1


@selfcert
@pytest.mark.parametrize("n, m", [(2, 1), (2, 2), (3, 1)])
def test_annulus_adjacency_and_arc_count(n, m):
    T = annulus_triangulation(n, m)
    assert len(T.arcs()) == T.surface.arc_count() == n + m
    for a in T.arcs():
        assert len(T.triangles_containing(a)) == 2


@selfcert
def test_wrong_arc_count_refused():
    S = MarkedSurface(0, (5,), 0)                        # pentagon: needs exactly 2 arcs
    with pytest.raises(QuiverlabError):
        Triangulation(S, triangles=((1, "b0_0", "b0_1"),))   # one triangle, one arc: wrong


@selfcert
def test_self_folded_refused():
    S = MarkedSurface(0, (4,), 0)
    with pytest.raises(QuiverlabError):
        Triangulation(S, triangles=((1, 1, "b0_0"),))        # arc 1 twice in one triangle
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.surfaces.triangulation`

- [ ] **Step 3: Implement**

```python
# src/quiverlab/surfaces/triangulation.py
"""Ideal triangulations of marked surfaces (Plan 48). A triangulation is a tuple of
ordered triangles, each an anticlockwise 3-tuple of SIDES (interior-arc ints + boundary-
segment strings). Construction validates the FST arc-adjacency: interior arcs in exactly
2 triangles, boundary segments in exactly 1, arc count == surface.arc_count(). Canonical
constructors for the disc fan, the annulus band, and the once-punctured torus. Float-free."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.marked import MarkedSurface


def _is_arc(side) -> bool:
    return isinstance(side, int) and not isinstance(side, bool)


@dataclass(frozen=True)
class Triangulation:
    surface: MarkedSurface
    triangles: tuple

    def __post_init__(self):
        object.__setattr__(self, "triangles",
                           tuple(tuple(t) for t in self.triangles))
        for tri in self.triangles:
            if len(tri) != 3:
                raise QuiverlabError("Triangulation: every triangle needs 3 sides",
                                     hint=f"got {tri!r}")
            if len(set(tri)) != 3:                        # a side repeated => self-folded
                raise QuiverlabError(
                    "Triangulation: a side appears twice in one triangle (self-folded)",
                    hint="self-folded triangles need a puncture and are out of v1 scope "
                         "(successor P48.1)")
        occ = Counter(s for tri in self.triangles for s in tri)
        for side, k in occ.items():
            want = 2 if _is_arc(side) else 1
            if k != want:
                raise QuiverlabError(
                    f"Triangulation: side {side!r} borders {k} triangle(s), expected {want}",
                    hint="interior arcs must border exactly 2 triangles, boundary "
                         "segments exactly 1")
        n_arcs = len(self.arcs())
        if n_arcs != self.surface.arc_count():
            raise QuiverlabError(
                f"Triangulation: {n_arcs} arcs but surface has {self.surface.arc_count()}",
                hint="the arc-count self-certificate failed (bad triangulation)")
        if len(self.triangles) != self.surface.triangle_count():
            raise QuiverlabError(
                f"Triangulation: {len(self.triangles)} triangles but surface has "
                f"{self.surface.triangle_count()}", hint="arc-count self-cert failed")

    def arcs(self) -> tuple:
        return tuple(sorted({s for tri in self.triangles for s in tri if _is_arc(s)}))

    def boundary_segments(self) -> tuple:
        return tuple(sorted({s for tri in self.triangles for s in tri if not _is_arc(s)},
                            key=str))

    def triangles_containing(self, arc) -> tuple:
        return tuple(tri for tri in self.triangles if arc in tri)


def fan_triangulation(marked) -> Triangulation:
    """The fan of the `marked`-gon from vertex 0 (vertices 0..marked-1 anticlockwise).
    Arc `i` (i = 1..marked-3) is the diagonal 0-(i+1). Triangle i (i=0..marked-3):
    (arc into vertex i+1, boundary seg (i+1)-(i+2), arc out to vertex i+2), with the two
    extreme triangles using boundary segments in place of a missing diagonal."""
    S = MarkedSurface(genus=0, boundary_marked=(marked,), punctures=0)
    # diagonals 0-(j) for j = 2..marked-2 -> arc id (j-1) in 1..marked-3.
    def diag(j):                                          # the arc id of diagonal 0-j
        return j - 1
    def bseg(i):                                          # boundary segment i-(i+1)
        return f"b0_{i}"
    tris = []
    for i in range(1, marked - 1):                        # triangle (0, i, i+1)
        left = bseg(0) if i == 1 else diag(i)             # side 0-i  (bd 0-1 for i==1)
        base = bseg(i)                                    # side i-(i+1)
        right = bseg(marked - 1) if i == marked - 2 else diag(i + 1)  # side (i+1)-0
        tris.append((left, base, right))                  # anticlockwise 0->i->(i+1)->0
    return Triangulation(S, tuple(tris))


def annulus_triangulation(n, m) -> Triangulation:
    """A standard band triangulation of the annulus C(n, m). Arcs cross between the outer
    boundary (n marked points) and the inner boundary (m marked points); the resulting
    quiver is affine Ã_{n+m-1} (Task 3). See the module-level derivation note."""
    S = MarkedSurface(genus=0, boundary_marked=(n, m), punctures=0)
    ...   # explicit band of n+m triangles / n+m arcs (see "Adjust to reality")
    return Triangulation(S, tuple(tris))


def once_punctured_torus() -> Triangulation:
    """The once-punctured torus (g=1, p=1): 3 loop-arcs at the single puncture, 2
    triangles. FST-admissible; OUT of v1 scope -- quiver_of/jacobian_of refuse it (Task 3)."""
    S = MarkedSurface(genus=1, boundary_marked=(), punctures=1)
    # the classical triangulation: arcs {1,2,3}, two triangles each using all three arcs.
    tris = ((1, 2, 3), (1, 3, 2))
    return Triangulation(S, tris)
```

**Adjust to reality (Task 2):**
- **The fan constructor is the arbiter's fixture** — its arc→vertex labelling (arc `j-1` =
  diagonal `0-j`) must make `quiver_of` produce the *linear* `A_n` quiver `1→2→…→n` in
  Task 3. If the disc arbiter test there fails, the fix is here (the diagonal-to-arc-id map
  or the anticlockwise side order), NOT a weakened assertion. Verify the two extreme
  triangles' boundary/diagonal sides by hand for the pentagon and hexagon before moving on.
- **`annulus_triangulation` is the one constructor needing care.** Write the band as
  `n + m` triangles wrapping the annulus: label the `n` outer segments and `m` inner
  segments, run a "zig-zag" of arcs between them, and let `Triangulation.__post_init__`'s
  arc-adjacency + arc-count self-cert be the loud guard that you wired it correctly (a
  mis-wired band fails the `== n+m` arc count or the "exactly 2 triangles per arc" check).
  For the **degenerate `C(1,1)`** (two arcs sharing both endpoints — a bigon of triangles),
  the raw adjacency is still 2-per-arc, but the angle bookkeeping in Task 3 is degenerate;
  ship `C(1,1)` only if the Task-3 2-cycle reduction certifies it against the Kronecker
  quiver, else restrict the annulus presets to `n + m >= 3` and note `C(1,1)` as a named
  edge. **Pin the smallest non-degenerate annulus (`C(2,1)` or `C(2,2)`) by hand.**
- `once_punctured_torus` is deliberately in-scope for *construction* (valid `Triangulation`)
  and out of scope for `quiver_of` — it is the LOUD-REFUSAL oracle (Task 3). Its arcs are
  loops (both endpoints the puncture); the arc-adjacency check still passes (each arc in 2
  triangles). Do NOT special-case it here.
- The boundary-segment id scheme `f"b{component}_{k}"` must be consistent between the
  surface's `Σ k_i` boundary-segment count and the triangulation; the arc-count self-cert
  only checks arcs, so add an internal assertion that the distinct boundary-segment count
  equals `c` if it helps catch labelling bugs early.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/surfaces/triangulation.py tests/families/test_surfaces_triangulation.py
git commit -m "feat(surfaces): Triangulation + arc-adjacency/arc-count self-cert + fan/annulus/once-punctured-torus constructors (self-folded refused)"
```

---

### Task 3: `quiver_of` / `potential_of` / `jacobian_of` — THE CRUX (convention + gentle tie)

The heart of the plan: turn a triangulation into a quiver-with-potential and hand it to
P44's Jacobian constructor. Two conventions are pinned by oracles — the angle→arrow
**orientation** (disc-A_n arbiter) and the potential's **internal-triangle** rule
(hexagon dim-6 pin) — and the whole thing is tied to P44+P46 by the **gentle**
cross-subsystem battery.

**Convention (arbitrated, not assumed).** `quiver_of(T)` has one vertex per interior arc.
For each triangle `Δ = (s0, s1, s2)` in **anticlockwise** cyclic order, and each cyclically
consecutive ordered pair `(s_i, s_{i+1})` of sides that are BOTH arcs, add an arrow
`s_i → s_{i+1}`. Then remove maximal sets of oriented **2-cycles** (the standard QP
reduction; for unpunctured surfaces this only fires on the minimal double-arc annulus).
The clockwise-vs-anticlockwise choice is exactly what the disc oracle fixes: *whichever
direction makes the fan of the `(n+3)`-gon give the linear `A_n` quiver `1→2→…→n` is the one
v1 ships.* Worked check (pentagon, `A_2`): fan arcs `α = 0-2` (id 1), `β = 0-3` (id 2); the
only triangle with two arcs is `(α, ·, β)`, giving `α → β`, i.e. `1 → 2` — linear `A_2`. ✔

**Potential.** `potential_of(T)` emits one 3-cycle term (coefficient `+1`) per **internal**
triangle — a triangle all three of whose sides are arcs — as the cycle `s0·s1·s2` in the
anticlockwise order. Its P44 cyclic derivatives are the three length-2 paths, the gentle
relations. For the **fan** (no internal triangle) `W = 0` and `Jac = kA_n` (hereditary). For
the **hexagon-with-internal-triangle** (arcs `{13, 35, 51}` forming a central triangle),
`W` is one 3-cycle and `Jac` is the triangle algebra of **dim 6** — exactly P44's
`test_triangle_jacobian_dimension`.

**Gentle (ABCP/LFS, the cross-subsystem oracle).** For every unpunctured-with-boundary
`T`, `Jac(quiver_of(T), potential_of(T))` is a **gentle** algebra — so
`is_gentle(jacobian_of(T))` is `True` across the disc / annulus / 2-holed zoo. This single
battery ties P44 (Jacobian) + P46 (`is_gentle`) + P48 (surface constructor).

**Files:**
- Create: `src/quiverlab/surfaces/qp.py`
- Modify: `src/quiverlab/surfaces/__init__.py` (export the public surface API)
- Test: `tests/families/test_surfaces_qp.py`

**Interfaces:**
- Consumes: `triangulation.py::Triangulation`, `combinat/quiver.py::Quiver`,
  `families/jacobian.py::Potential`/`JacobianAlgebra` (P44),
  `invariants/recognizers.py::is_gentle` (P46/P38), `quiverlab.errors::QuiverlabError`.
- Produces:
  ```python
  def quiver_of(T) -> Quiver
      # one vertex per interior arc (vertices = T.arcs()), one arrow per anticlockwise
      # arc-to-arc angle, then 2-cycle reduced. Loud on out-of-v1-scope T (punctures /
      # closed / self-folded).
  def potential_of(T) -> Potential
      # a P44 Potential over quiver_of(T): one +1 3-cycle per internal triangle.
  def jacobian_of(T, field=None, degree_bound=None) -> Algebra
      # JacobianAlgebra(quiver_of(T), potential_of(T), field). NotFiniteDimensionalError
      # propagates (should not fire in v1 scope -- gentle => finite -- but stays honest).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_surfaces_qp.py
"""Q(T)/W(T)/Jac(T) (Plan 48, deep). ARBITER: fan of the (n+3)-gon -> linear A_n quiver
EXACTLY (fixes the angle orientation). Literature: hexagon-internal-triangle Jac dim 6.
Cross-subsystem (P44+P46+P48): is_gentle(jacobian_of(T)) True across the surface zoo.
Loud: punctures / closed surfaces are out of v1 scope."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import is_gentle
from quiverlab.surfaces.qp import jacobian_of, potential_of, quiver_of
from quiverlab.surfaces.triangulation import (annulus_triangulation, fan_triangulation,
                                             once_punctured_torus, Triangulation)

arb = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _linear_A(n):
    arrows = {f"a{i}": (i, i + 1) for i in range(1, n)}
    return Quiver(list(range(1, n + 1)), arrows)


@arb
@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_disc_fan_is_linear_A_n(n):
    # THE ARBITER: fan of the (n+3)-gon -> linear A_n quiver 1->2->...->n.
    Q = quiver_of(fan_triangulation(n + 3))
    verts = list(Q.vertices)
    assert len(verts) == n
    # underlying arrows == a linear chain on the arc ids (up to the arc-id bijection):
    ends = sorted((Q.source(a), Q.target(a)) for a in Q.arrows)
    chain = sorted((v, v + 1) for v in verts[:-1]) if n >= 2 else []
    assert len(ends) == max(n - 1, 0)
    # the arc ids 1..n are consecutive and each consecutive pair carries exactly one arrow
    assert {tuple(sorted(e)) for e in ends} == {tuple(sorted(c)) for c in chain}


@lit
def test_hexagon_internal_triangle_jacobian_dim_6():
    # hexagon (6 marked pts) triangulated with the central triangle {13,35,51}: one
    # internal triangle -> one 3-cycle -> Jac = triangle algebra dim 6 (P44's pin).
    T = _hexagon_internal_triangle()                     # helper below
    J = jacobian_of(T, field=QQ)
    assert J.dim == 6
    assert is_gentle(J)


@lit
def test_disc_fan_jacobian_is_kA_n():
    # fan => no internal triangle => W = 0 => Jac = kA_n (hereditary). D_5 (pentagon) -> kA_2.
    J = jacobian_of(fan_triangulation(5), field=QQ)      # pentagon
    assert J.dim == 3                                    # kA_2: 2 vertices + 1 arrow
    assert is_gentle(J)


@arb
@pytest.mark.parametrize("factory", [
    lambda: fan_triangulation(5), lambda: fan_triangulation(7),
    lambda: annulus_triangulation(2, 1), lambda: annulus_triangulation(2, 2),
    lambda: _hexagon_internal_triangle(),
])
def test_surface_jacobians_are_gentle(factory):
    # ABCP/LFS: unpunctured-with-boundary surface Jacobians are gentle. P44+P46+P48 tie.
    J = jacobian_of(factory(), field=QQ)
    assert is_gentle(J)


@lit
def test_annulus_C21_is_affine_A2_shape():
    # smallest non-degenerate annulus: C(2,1) -> 3 arcs, affine Ã_2 (acyclic), hereditary,
    # gentle. Pin the exact arrow multiset by hand in the helper (arbitrated by the
    # disc-fixed convention + finiteness).
    Q = quiver_of(annulus_triangulation(2, 1))
    assert len(list(Q.vertices)) == 3
    assert len(Q.arrows) == 3                            # affine Ã_2 cycle, acyclically oriented
    J = jacobian_of(annulus_triangulation(2, 1), field=QQ)
    assert is_gentle(J)


@selfcert
def test_punctured_and_closed_refused():
    with pytest.raises(QuiverlabError):
        quiver_of(once_punctured_torus())                # p>0 out of v1 scope
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.surfaces.qp`

- [ ] **Step 3: Implement**

```python
# src/quiverlab/surfaces/qp.py
"""Quiver with potential of an ideal triangulation (Plan 48). For an UNPUNCTURED surface
with boundary: quiver_of(T) has one vertex per arc and one arrow per anticlockwise angle
between consecutive arcs in a triangle (2-cycle reduced); potential_of(T) is the sum of
3-cycles over internal triangles; jacobian_of(T) is P44's JacobianAlgebra. By ABCP 2010
/ Labardini 2009 the Jacobian is a GENTLE algebra (self-certified via P46's is_gentle).
Punctures / closed surfaces / self-folded triangles are OUT of v1 scope (loud). Float-free."""
from __future__ import annotations

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.families.jacobian import JacobianAlgebra, Potential


def _require_v1_scope(T):
    S = T.surface
    if S.punctures > 0:
        raise QuiverlabError(
            "surfaces: punctured surfaces are out of v1 scope",
            hint="v1 handles unpunctured surfaces with boundary; puncture potentials "
                 "and self-folded triangles are the successor (P48.1)")
    if S.num_boundary_components == 0:
        raise QuiverlabError(
            "surfaces: closed surfaces are out of v1 scope",
            hint="v1 requires non-empty boundary (successor P48.1)")


def _is_arc(side):
    return isinstance(side, int) and not isinstance(side, bool)


def _raw_arrows(T):
    """The angle arrows before 2-cycle reduction: for each anticlockwise triangle, each
    cyclically-consecutive ordered pair of ARC sides (a, b) yields an arrow a -> b."""
    edges = []                                            # list of (src_arc, tgt_arc)
    for tri in T.triangles:
        for i in range(3):
            a, b = tri[i], tri[(i + 1) % 3]
            if _is_arc(a) and _is_arc(b):
                edges.append((a, b))
    return edges


def _reduce_two_cycles(edges):
    """Cancel oriented 2-cycles pairwise: while both (a,b) and (b,a) occur, drop one of
    each. Returns the surviving multiset of directed edges (FST reduced quiver)."""
    from collections import Counter
    cnt = Counter(edges)
    for (a, b) in list(cnt):
        if a < b:                                         # handle each unordered pair once
            k = min(cnt[(a, b)], cnt[(b, a)])
            cnt[(a, b)] -= k
            cnt[(b, a)] -= k
    return [e for e, c in cnt.items() for _ in range(c)]


def quiver_of(T) -> Quiver:
    _require_v1_scope(T)
    arcs = list(T.arcs())
    edges = _reduce_two_cycles(_raw_arrows(T))
    arrows, k = {}, 0
    for (a, b) in sorted(edges):
        k += 1
        arrows[f"m{k}"] = (a, b)
    return Quiver(arcs, arrows)


def _internal_triangles(T):
    return [tri for tri in T.triangles if all(_is_arc(s) for s in tri)]


def potential_of(T) -> Potential:
    Q = quiver_of(T)
    # map each internal triangle's anticlockwise arc-cycle to the arrow word of Q. Because
    # Q is 2-cycle reduced, an internal triangle's three angle-arrows survive (no boundary
    # side => no cancellation with a boundary), so the 3-cycle is a genuine cyclic word.
    terms = []
    by_ends = {(Q.source(a), Q.target(a)): a for a in Q.arrows}
    for tri in _internal_triangles(T):
        word = tuple(by_ends[(tri[i], tri[(i + 1) % 3])] for i in range(3))
        terms.append((1, word))
    return Potential(Q, terms)


def jacobian_of(T, field=None, degree_bound=None):
    Q = quiver_of(T)
    W = potential_of(T)
    return JacobianAlgebra(Q, W, field=field, degree_bound=degree_bound)
```

**Adjust to reality (Task 3):**
- **The disc arbiter is `test_disc_fan_is_linear_A_n`.** If the fan gives the *reversed*
  chain (`n→…→1`) or a non-chain, flip the anticlockwise pair to `(s_{i+1}, s_i)` in
  `_raw_arrows` ONCE and document the outcome in the docstring (the house "oracle decides
  the orientation" pattern, P37/P40/P44). Do NOT relax the assertion — the fan MUST be
  linear `A_n`.
- `potential_of`'s `by_ends` lookup assumes each internal triangle's three angle-arrows
  survive reduction with distinct endpoints — true because an internal triangle shares no
  side with a boundary and (in v1 scope, no double-arc) its three arcs are pairwise
  distinct with distinct endpoints. If a `KeyError` fires, a 2-cycle reduction wrongly
  cancelled an internal-triangle arrow (only possible in a degenerate annulus) — that is
  the signal to restrict the annulus presets to the non-degenerate range, not to guess.
- **`JacobianAlgebra` is the finiteness certificate** (P44): in v1 scope it must NOT raise
  `NotFiniteDimensionalError` (gentle ⇒ finite). If it does on the annulus, the quiver has
  a *fully-oriented* cycle with `W=0` — meaning the orientation convention put all annulus
  arrows the same way round (wrong; the affine `Ã` quiver is acyclic, one source + one
  sink). The finiteness certificate thus doubly-guards the annulus orientation.
- `_hexagon_internal_triangle()` (test helper): `MarkedSurface(0, (6,), 0)`, arcs
  `{1: 1-3, 2: 3-5, 3: 5-1}` forming the central triangle `(1, 2, 3)` plus three "ear"
  triangles each with one arc + two boundary segments. Assemble the six boundary segments
  `b0_0..b0_5` and the four triangles; the `Triangulation` arc-count self-cert (3 arcs)
  guards the assembly. This helper is the hexagon preset's ground truth (Task 6).
- `is_gentle(J)` is a P46/P38 recognizer over the presentation — field-agnostic, so the
  gentle battery can run over `QQ` (exact). Keep the presets (Task 6) over `CC`/`GF` per the
  existing preset conventions; the gentle *recognizer* verdict is identical.

- [ ] **Step 4: Run tests** — Expected: PASS (orientation flipped once if needed)
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/surfaces/qp.py src/quiverlab/surfaces/__init__.py \
        tests/families/test_surfaces_qp.py
git commit -m "feat(surfaces): quiver_of (disc-A_n-arbitrated angles + 2-cycle reduction) / potential_of (internal-triangle 3-cycles) / jacobian_of -- gentle (ABCP) tie to P44+P46"
```

---

### Task 4: `flip` + `certify_flip_mutation` (flip ↔ Fomin–Zelevinsky mutation)

**Files:**
- Create: `src/quiverlab/surfaces/flip.py`
- Modify: `src/quiverlab/surfaces/__init__.py` (export `flip`, `certify_flip_mutation`)
- Test: `tests/families/test_surfaces_flip.py`

**Interfaces:**
- Consumes: `triangulation.py::Triangulation`, `qp.py::quiver_of`,
  `combinat/quiver.py::Quiver`.
- Produces:
  ```python
  def flip(T, arc) -> Triangulation
      # the combinatorial quadrilateral flip: replace interior arc `arc` (a diagonal of the
      # quad formed by its 2 adjacent triangles) with the other diagonal. Loud on a boundary
      # segment (not flippable) or a non-flippable arc (self-folded -- out of scope).
  def exchange_matrix(Q) -> list          # skew-symmetric B: B[i][j] = #(i->j) - #(j->i)
  def matrix_mutation(B, k) -> list       # Fomin-Zelevinsky mutation μ_k (exact integers)
  def certify_flip_mutation(T, arc) -> bool
      # True iff exchange_matrix(quiver_of(flip(T, arc))) == matrix_mutation(
      #   exchange_matrix(quiver_of(T)), arc), under the arc-index identification (the
      # flipped arc inherits `arc`'s vertex slot). The per-instance flip<->mutation
      # certificate (FST 2008). Full DWZ potential right-equivalence is NOT attempted.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_surfaces_flip.py
"""Flip <-> Fomin-Zelevinsky mutation (Plan 48, deep). Cross-engine: quiver_of(flip(T,a))
matches the exact matrix mutation of quiver_of(T) at vertex a, on every interior arc of
the disc/annulus zoo. Self-cert: flip is an involution."""
import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.flip import (certify_flip_mutation, exchange_matrix, flip,
                                     matrix_mutation)
from quiverlab.surfaces.qp import quiver_of
from quiverlab.surfaces.triangulation import annulus_triangulation, fan_triangulation

xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


@selfcert
@pytest.mark.parametrize("factory", [lambda: fan_triangulation(6),
                                     lambda: fan_triangulation(7),
                                     lambda: annulus_triangulation(2, 2)])
def test_flip_is_an_involution(factory):
    T = factory()
    for a in T.arcs():
        T2 = flip(flip(T, a), a)
        assert set(T2.arcs()) == set(T.arcs())
        # same quiver back (up to arrow renaming): compare exchange matrices
        assert exchange_matrix(quiver_of(T2)) == exchange_matrix(quiver_of(T))


@xeng
@pytest.mark.parametrize("factory", [lambda: fan_triangulation(6),
                                     lambda: fan_triangulation(7),
                                     lambda: annulus_triangulation(2, 1)])
def test_every_flip_certifies_against_matrix_mutation(factory):
    T = factory()
    for a in T.arcs():
        assert certify_flip_mutation(T, a)


@xeng
def test_matrix_mutation_is_an_involution_on_B():
    Q = quiver_of(fan_triangulation(6))
    B = exchange_matrix(Q)
    k = 1
    assert matrix_mutation(matrix_mutation(B, k), k) == B


@selfcert
def test_flip_of_boundary_segment_refused():
    T = fan_triangulation(6)
    with pytest.raises(QuiverlabError):
        flip(T, "b0_0")                                  # boundary segment: not flippable
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.surfaces.flip`

- [ ] **Step 3: Implement**

```python
# src/quiverlab/surfaces/flip.py
"""Quadrilateral flips + the flip<->mutation certificate (Plan 48). Flipping an interior
arc replaces it with the other diagonal of the quadrilateral formed by its two adjacent
triangles. The Fomin-Shapiro-Thurston theorem: quiver_of(flip(T, a)) equals the Fomin-
Zelevinsky matrix mutation of quiver_of(T) at vertex a. We certify that identity per
instance (exact skew-symmetric integer combinatorics); full DWZ right-equivalence of the
potentials is a named successor, not attempted here. Float-free."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.qp import quiver_of
from quiverlab.surfaces.triangulation import Triangulation, _is_arc


def flip(T, arc) -> Triangulation:
    if not _is_arc(arc):
        raise QuiverlabError(f"flip: {arc!r} is a boundary segment, not flippable",
                             hint="only interior arcs flip")
    tris = T.triangles_containing(arc)
    if len(tris) != 2:
        raise QuiverlabError(f"flip: arc {arc!r} does not border 2 distinct triangles",
                             hint="non-flippable (self-folded configuration, out of scope)")
    # the two triangles share `arc`; their other four sides form a quadrilateral. The new
    # arc joins the two "opposite corners". Rebuild the two triangles across the OTHER
    # diagonal, keeping `arc`'s id for the new diagonal (mutation keeps the vertex slot).
    ...   # combinatorial quad flip; see "Adjust to reality"
    return Triangulation(T.surface, new_tris)


def exchange_matrix(Q):
    verts = list(Q.vertices)
    idx = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    B = [[0] * n for _ in range(n)]
    for a in Q.arrows:
        i, j = idx[Q.source(a)], idx[Q.target(a)]
        B[i][j] += 1
        B[j][i] -= 1
    return B


def matrix_mutation(B, k):
    """Fomin-Zelevinsky μ_k on a skew-symmetric integer matrix (indices are VERTEX
    POSITIONS -- caller aligns k to the arc's position)."""
    n = len(B)
    Bp = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == k or j == k:
                Bp[i][j] = -B[i][j]
            else:
                Bp[i][j] = B[i][j] + (abs(B[i][k]) * B[k][j]
                                      + B[i][k] * abs(B[k][j])) // 2
    return Bp


def certify_flip_mutation(T, arc) -> bool:
    Q = quiver_of(T)
    verts = list(Q.vertices)
    k = verts.index(arc)
    B_mut = matrix_mutation(exchange_matrix(Q), k)
    # quiver_of(flip) keeps the same vertex SET (the flipped arc reuses `arc`'s id), so the
    # exchange matrices are index-aligned by the shared vertex ordering.
    Qf = quiver_of(flip(T, arc))
    if list(Qf.vertices) != verts:
        # re-index Qf onto verts' order (same set, possibly re-sorted) before comparing
        B_flip = _exchange_on_order(Qf, verts)
    else:
        B_flip = exchange_matrix(Qf)
    return B_flip == B_mut
```

**Adjust to reality (Task 4):**
- **The quad-flip is the one structural routine.** Given the two triangles sharing `arc`,
  identify the four other sides and the two "opposite corner" endpoints; the new arc runs
  between them, and the two new triangles are `(corner-A-sides..., new_arc)` and
  `(corner-B-sides..., new_arc)` in anticlockwise order. The **arbiter is
  `certify_flip_mutation`** — if a flip's quiver does not match the FZ mutation, the
  anticlockwise side-ordering of the two rebuilt triangles is wrong (swap once), exactly the
  house oracle-decides pattern. Keep `arc`'s id for the new diagonal so the exchange
  matrices stay index-aligned (this is why the FST theorem reads as `μ_k`).
- `_exchange_on_order(Q, order)`: the `exchange_matrix` built with a supplied vertex order
  (a 5-line variant) — needed because `flip` may return arcs in a different sort order.
- The FZ mutation is the exact `Cluster algebras I` (Fomin–Zelevinsky 2002) formula — cite
  `fomin_zelevinsky_ca1`. Its involutivity (`μ_k∘μ_k = id`) is the self-cert;
  flip↔mutation is the cross-engine oracle; the FST flip theorem is the literature anchor
  (`fomin_shapiro_thurston`).
- **Honest scope (record in the docstring + verification page):** this certifies the QUIVER
  under flip↔mutation. The DWZ **potential** right-equivalence under mutation (that
  `μ_k(Q, W)` is right-equivalent to `(quiver_of(flip), potential_of(flip))`) is NOT
  attempted — the named successor. For the gentle v1 scope the quiver-level certificate
  plus `is_gentle` on both sides is the shipped guarantee.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/surfaces/flip.py src/quiverlab/surfaces/__init__.py \
        tests/families/test_surfaces_flip.py
git commit -m "feat(surfaces): flip (quadrilateral) + certify_flip_mutation vs Fomin-Zelevinsky matrix mutation (per-instance flip<->mutation; DWZ potential deferred) + involution"
```

---

### Task 5: `surface_block` descriptor + the AG (gentle-dictionary) tie

Surfaces are an **input method**, not a new computation — so there is no new webapp compute
kind (the produced algebra flows through every existing kind). `surface_block` is a
lightweight descriptor for the report and the preset provenance: the surface invariants,
the resulting quiver/relations, the `is_gentle` verdict, and — because v1 Jacobians are
gentle — the P46 **AG invariant** (Avella-Alaminos–Geiss), tying the gentle-dictionary
loop (metaplan §2: "gentle ↔ dissection dictionary").

**Files:**
- Create: `src/quiverlab/surfaces/block.py`
- Modify: `src/quiverlab/surfaces/__init__.py` (export `surface_block`)
- Test: `tests/families/test_surfaces_block.py`

**Interfaces:**
- Consumes: `qp.py::quiver_of`/`jacobian_of`, `invariants/recognizers.py::is_gentle`
  (P46/P38), `strings/ag.py::ag_invariant` (P46).
- Produces:
  ```python
  def surface_block(T, field=None) -> dict
      # {"kind": "surface",
      #  "genus", "boundary_marked", "punctures", "arc_count", "triangle_count",
      #  "euler_characteristic",
      #  "vertices": [...], "arrows": {name: [s, t]}, "relations": [...],
      #  "dim": int, "is_gentle": bool,
      #  "ag_invariant": [[n, m], ...]  # present iff gentle (a multiset, sorted),
      #  "references": ["fomin_shapiro_thurston", "labardini", "abcp"],
      #  "citations": [...]}
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_surfaces_block.py
"""surface_block descriptor + the AG (gentle-dictionary) tie (Plan 48, deep).
Self-cert: the block's invariants match the surface; gentle Jacobians carry an AG
invariant equal to strings.ag.ag_invariant of the SAME algebra."""
import pytest

from quiverlab.fields import QQ
from quiverlab.strings.ag import ag_invariant
from quiverlab.surfaces.block import surface_block
from quiverlab.surfaces.qp import jacobian_of
from quiverlab.surfaces.triangulation import annulus_triangulation, fan_triangulation

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@selfcert
def test_block_invariants_match_surface():
    T = fan_triangulation(6)                              # A_3
    b = surface_block(T, field=QQ)
    assert b["kind"] == "surface" and b["arc_count"] == 3
    assert b["euler_characteristic"] == 1 and b["is_gentle"] is True
    assert b["dim"] == jacobian_of(T, field=QQ).dim


@xeng
@pytest.mark.parametrize("factory", [lambda: fan_triangulation(7),
                                     lambda: annulus_triangulation(2, 2)])
def test_block_ag_invariant_matches_engine(factory):
    T = factory()
    b = surface_block(T, field=QQ)
    J = jacobian_of(T, field=QQ)
    engine_ag = [list(pair) for pair in ag_invariant(J).pairs]
    assert b["ag_invariant"] == engine_ag
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.surfaces.block`

- [ ] **Step 3: Implement**

```python
# src/quiverlab/surfaces/block.py
"""surface_block: the algebra-adjacent descriptor of a triangulated surface (Plan 48).
Packages the surface invariants + the resulting gentle algebra's presentation, gentle
verdict, and AG invariant (P46) for the report and preset provenance. Not a webapp compute
kind -- surfaces are an INPUT method; the produced algebra runs through every existing
compute kind. Float-free."""
from __future__ import annotations

from quiverlab.invariants.recognizers import is_gentle
from quiverlab.surfaces.qp import jacobian_of, quiver_of

_REFS = ("fomin_shapiro_thurston", "labardini", "abcp")


def surface_block(T, field=None) -> dict:
    S = T.surface
    Q = quiver_of(T)
    J = jacobian_of(T, field=field)
    gentle = is_gentle(J)
    block = {
        "kind": "surface",
        "genus": S.genus, "boundary_marked": list(S.boundary_marked),
        "punctures": S.punctures, "arc_count": S.arc_count(),
        "triangle_count": S.triangle_count(),
        "euler_characteristic": S.euler_characteristic(),
        "vertices": list(Q.vertices),
        "arrows": {a: list(Q.arrows[a]) for a in Q.arrows},
        "relations": [str(r) for r in (J.relations or ())],
        "dim": J.dim, "is_gentle": gentle,
        "references": list(_REFS), "citations": list(_REFS),
    }
    if gentle:
        from quiverlab.strings.ag import ag_invariant
        block["ag_invariant"] = [list(pair) for pair in ag_invariant(J).pairs]
    return block
```

**Adjust to reality (Task 5):**
- Read `strings/ag.py::AGInvariant` (P46) FIRST — confirm the attribute is `.pairs` (a
  sorted tuple of `(n, m)`); if P46 named it differently, mirror the actual field. The AG
  invariant is asked only when `is_gentle` is True (P46 refuses non-gentle input loudly), so
  the `if gentle:` guard is load-bearing.
- `surface_block` is deliberately NOT wired as a webapp/GUI compute kind (no `MODULE_KINDS`
  / `_dispatch` edit). It is consumed by the report's provenance line for a surface-derived
  preset and is available programmatically. State on the verification page that the
  interactive surface *compute* surface is folded into the existing algebra kinds (surfaces
  are input), and the free-form draw-a-surface canvas is the named successor.
- No `viz` change: `viz/draw.py::draw_quiver` and `tikz.py::tikz_quiver` already render the
  Jacobian's quiver+relations (`A.draw()`/`A.tikz()`); a triangulation/potential geometric
  renderer is out of scope (there is no such function in `viz` today) and is recorded as a
  successor.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/surfaces/block.py src/quiverlab/surfaces/__init__.py \
        tests/families/test_surfaces_block.py
git commit -m "feat(surfaces): surface_block descriptor + AG-invariant gentle-dictionary tie (P46)"
```

---

### Task 6: GUI/no-code story — catalog registration, three presets, QPA gentle parity

Surfaces are **input-side**: the three presets let a user pick a triangulated surface's
algebra on the canvas with zero code, and the catalog lists the surface constructors for
`families()` discoverability. **No new compute kind** (the resulting algebra flows through
all existing kinds) — so no GUI `MODULE_KINDS`/i18n/`gui.js` edits beyond the data-driven
preset dropdown (which needs no JS change).

**Files:**
- Modify: `src/quiverlab/__init__.py` (top-level export the surface API:
  `MarkedSurface`, `Triangulation`, `fan_triangulation`, `annulus_triangulation`,
  `once_punctured_torus`, `quiver_of`, `potential_of`, `jacobian_of`, `flip`,
  `certify_flip_mutation`, `surface_block`)
- Modify: `src/quiverlab/families/discover.py` (three `FamilyInfo` CATALOG entries)
- Modify: `webapp/server/catalog.py` (`_iter_families` skip-set — the surface
  constructors take non-scalar args, like the P44 `JacobianAlgebra`/`OnePointExtension`
  precedent and the `"zoo"` skip)
- Modify: `scripts/gui_build_hook.py::_preset_algebras()` (three drawable presets) +
  regenerate `webapp/static/gui/presets.json`
- Create: `tests/qpa/test_surfaces_qpa.py`
- Test: `tests/families/test_surfaces_catalog.py`, the existing
  `tests/gui/test_build_hook.py::test_presets_round_trip_through_runner`

- [ ] **Step 1: Catalog + exports.** Add to `families/discover.py::CATALOG`:

```python
FamilyInfo("fan_triangulation", "fan_triangulation(marked)",
           "surface", ("fomin_shapiro_thurston", "abcp", "labardini"),
           "Fan triangulation of a disc with `marked` boundary points (type A_{marked-3})."),
FamilyInfo("annulus_triangulation", "annulus_triangulation(n, m)",
           "surface", ("fomin_shapiro_thurston", "abcp", "labardini"),
           "Band triangulation of the annulus C(n, m) (affine A~ shape)."),
FamilyInfo("jacobian_of", "jacobian_of(Triangulation, field=...)",
           "surface", ("fomin_shapiro_thurston", "abcp", "derksen_weyman_zelevinsky"),
           "Gentle Jacobian algebra of an ideal triangulation (unpunctured, boundary)."),
```

The `"surface"` route is a NEW route string — extend the `FamilyInfo.route` doc-comment
(`"monomial" | "general" | "structure-constant" | "iterator" | "surface"`). Each `name` must
be a top-level `quiverlab` export (the `getattr(quiverlab, name)` contract in
`webapp/server/catalog.py::_iter_families`), so export all three from
`src/quiverlab/__init__.py`.

- [ ] **Step 2: Skip-set.** In `webapp/server/catalog.py::_iter_families`, extend the
  single-name skip (currently `if name == "zoo": continue`) to also skip the surface
  constructors, which take `Triangulation`/`MarkedSurface`/int args and do NOT fit the
  scalar auto-form (`_params_of` only classifies `bool`/`int`/`str`):

```python
_SURFACE_INPUT = {"fan_triangulation", "annulus_triangulation", "jacobian_of"}
...
    if name == "zoo" or name in _SURFACE_INPUT:
        continue
```

- [ ] **Step 3: Presets.** Add three entries to
  `scripts/gui_build_hook.py::_preset_algebras()` — built via the library (never
  hand-authored JSON; `generate_presets()` reads `A.quiver`/`A.relations` back out, and the
  round-trip test rebuilds each preset through the runner):

```python
import quiverlab as ql
from quiverlab.surfaces.qp import jacobian_of
from quiverlab.surfaces.triangulation import fan_triangulation, annulus_triangulation
...
entries.append(("Surface: disc fan A3 (CC)",
                jacobian_of(fan_triangulation(6), field=ql.CC), {"kind": "CC"}))
entries.append(("Surface: annulus C(2,2) affine A~3 (CC)",
                jacobian_of(annulus_triangulation(2, 2), field=ql.CC), {"kind": "CC"}))
entries.append(("Surface: hexagon w/ internal triangle -> 3-cycle gentle (CC)",
                jacobian_of(_hexagon_internal_triangle(), field=ql.CC), {"kind": "CC"}))
```

`_hexagon_internal_triangle()` is promoted from the Task-3 test helper into a small module
function (e.g. `surfaces/triangulation.py::hexagon_with_internal_triangle()`) so both the
test and the build hook build the same triangulation. Then **regenerate the committed
`webapp/static/gui/presets.json`** to match `generate_presets()` output (by hand — no sync
script; the round-trip test `tests/gui/test_build_hook.py` asserts `len >= 4` and that every
preset rebuilds as a valid quiver through the runner). The canvas lays vertices on a circle
regardless of geometry (no per-preset coordinate hint in the schema) — acceptable for v1;
the geometric surface-drawing canvas is the deferred successor.

- [ ] **Step 4: Catalog + preset tests.**

```python
# tests/families/test_surfaces_catalog.py
"""Surface constructors are catalog-listed + top-level exported, skipped in the webapp
family form (non-scalar args, the zoo precedent) (Plan 48)."""
import quiverlab as ql
from quiverlab.families.discover import families


def test_surface_constructors_catalogued_and_exported():
    names = families().names()
    for nm in ("fan_triangulation", "annulus_triangulation", "jacobian_of"):
        assert nm in names
        assert getattr(ql, nm, None) is not None          # top-level export contract


def test_webapp_form_skips_surface_inputs():
    from webapp.server.catalog import _iter_families
    listed = {name for name, _ in _iter_families()}
    assert "fan_triangulation" not in listed and "jacobian_of" not in listed
```

- [ ] **Step 5: QPA gentle parity** (probe-first, the P46 precedent — QPA has
  `IsGentleAlgebra`/`IsSpecialBiserialAlgebra`; it has NO surface/triangulation surface, so
  we crosscheck the RESULTING gentle algebra, not the surface):

```python
# tests/qpa/test_surfaces_qpa.py
"""Surface Jacobians are gentle -- crosscheck our is_gentle against QPA's IsGentleAlgebra
on the surface algebras (Plan 48). QPA has no triangulation surface (NamesGVars sweep),
so the comparison is at the resulting-algebra level, mirroring P46's recognizer parity."""
import pytest

from quiverlab import qpa
from quiverlab.fields import GF
from quiverlab.invariants.recognizers import is_gentle
from quiverlab.surfaces.qp import jacobian_of
from quiverlab.surfaces.triangulation import (annulus_triangulation, fan_triangulation)

pytestmark = pytest.mark.skipif(qpa.session.should_skip_qpa(), reason="QPA unavailable")


@pytest.mark.parametrize("factory", [lambda: fan_triangulation(6),
                                     lambda: annulus_triangulation(2, 1)])
def test_surface_jacobian_gentle_matches_qpa(factory):
    J = jacobian_of(factory(), field=GF(3))
    ours = is_gentle(J)
    theirs = qpa.crosscheck.is_gentle_via_qpa(J)          # existing P46 QPA bridge
    assert ours == theirs is True
    # standing guard: QPA has no surface-triangulation surface today
    assert "TriangulationOfSurface" not in qpa.session.names_gvars()
```

- [ ] **Step 6: Run the gates**

Run: `... -m pytest tests/families/test_surfaces_catalog.py tests/gui/test_build_hook.py -q`
then `... -m pytest tests/qpa/test_surfaces_qpa.py -q -m qpa` (venv has `[qpa]`).
Expected: PASS; the QPA surface-surface guard holds (honest scope: QPA has none today).

- [ ] **Step 7: Commit**

```bash
git add src/quiverlab/__init__.py src/quiverlab/families/discover.py \
        webapp/server/catalog.py scripts/gui_build_hook.py \
        webapp/static/gui/presets.json src/quiverlab/surfaces/triangulation.py \
        tests/families/test_surfaces_catalog.py tests/qpa/test_surfaces_qpa.py
git commit -m "feat(surfaces,gui,webapp,qpa): catalog + top-level exports + 3 drawable presets (disc/annulus/hexagon) + webapp form skip + IsGentleAlgebra QPA parity"
```

---

### Task 7: verification page, citations, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-05-metaplan-v0.2.0.md` (tick the P48 card delivery note)
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`)

- [ ] **Step 1: Citations** (BibTeX-VERIFIED only; `_r(key, bibtex_key, kind, title,
  annotation, *tags)` registry precedent). Reuse P44's `labardini`
  (LabardiniFragoso2009), `derksen_weyman_zelevinsky` (DWZ 2008), and `assem_book`. Add:

```bibtex
@article{FominShapiroThurston2008,
  author  = {Fomin, Sergey and Shapiro, Michael and Thurston, Dylan},
  title   = {Cluster algebras and triangulated surfaces. {P}art {I}:
             {C}luster complexes},
  journal = {Acta Mathematica},
  volume  = {201}, number = {1}, pages = {83--146}, year = {2008},
}
@article{FominZelevinsky2002,
  author  = {Fomin, Sergey and Zelevinsky, Andrei},
  title   = {Cluster algebras {I}: {F}oundations},
  journal = {Journal of the American Mathematical Society},
  volume  = {15}, number = {2}, pages = {497--529}, year = {2002},
}
@article{ABCP2010,
  author  = {Assem, Ibrahim and Br{\"u}stle, Thomas and
             Charbonneau-Jodoin, Gabrielle and Plamondon, Pierre-Guy},
  title   = {Gentle algebras arising from surface triangulations},
  journal = {Algebra \& Number Theory},
  volume  = {4}, number = {2}, pages = {201--229}, year = {2010},
}
```

and in `registry.py` (mirror the `_r("assem_book", "ASS2006", "foundation", ...)` shape):

```python
_r("fomin_shapiro_thurston", "FominShapiroThurston2008", "foundation",
   "Cluster algebras and triangulated surfaces I",
   "The arc/triangulation combinatorics: the ideal-arc count n=6g-6+3(b+p)+Σk_i, the "
   "admissibility exclusion list, and flip<->mutation -- the ground truth for the surface "
   "subsystem.", "families", "surfaces"),
_r("fomin_zelevinsky_ca1", "FominZelevinsky2002", "foundation",
   "Cluster algebras I: Foundations",
   "The skew-symmetric matrix mutation mu_k that flip is certified against.", "families"),
_r("abcp", "ABCP2010", "family",
   "Gentle algebras arising from surface triangulations",
   "For an unpunctured surface with boundary the Jacobian Jac(Q(T),W(T)) is a GENTLE "
   "algebra -- the v1 certifiability theorem (with Labardini 2009).", "families", "surfaces"),
```

- [ ] **Step 2: Verification page.** Add the Plan-48 subsystem rows:
  - `surfaces/marked.py` — `oracle_literature` (disc/annulus/torus arc counts, the FST
    admissibility refusals); `oracle_selfcert` (the two side-counting identities).
  - `surfaces/triangulation.py` — `oracle_selfcert` (arc-adjacency + arc-count self-cert on
    every constructor; self-folded refusal).
  - `surfaces/qp.py` — `oracle_crossengine` (the disc-fan-`A_n` ARBITER; `is_gentle(
    jacobian_of(T))` True across the disc/annulus/2-holed zoo — the P44+P46+P48 tie);
    `oracle_literature` (hexagon-internal-triangle `dim 6`; disc-fan `kA_n`; the small
    annulus-`Ã` pin).
  - `surfaces/flip.py` — `oracle_crossengine` (flip ≡ FZ matrix mutation on every interior
    arc); `oracle_selfcert` (flip involution; `μ_k∘μ_k = id`).
  - `surfaces/block.py` — `oracle_crossengine` (block AG invariant ≡ `strings.ag.ag_invariant`).
  - `qpa` — `IsGentleAlgebra`/`IsSpecialBiserialAlgebra` parity on the surface algebras.

  Add the **honest-scope entries**: (a) **v1 = unpunctured surfaces with boundary only** —
  punctures, closed surfaces, and self-folded triangles are refused loudly (successor
  P48.1: puncture potentials + self-folded triangles + the once-punctured-torus / Markov
  quiver); the once-punctured torus is the pinned loud-refusal oracle. (b) **flip↔mutation
  is certified at the QUIVER level** (exact FZ matrix mutation) — full DWZ *potential*
  right-equivalence is a named successor, not attempted. (c) The **AG invariant is a derived
  invariant, NOT complete** (inherited from P46; never claim completeness). (d) **QPA has no
  surface/triangulation surface** (`NamesGVars()` sweep) — crosschecks are at the resulting
  gentle-algebra level (`IsGentleAlgebra`), with the standing guard that FAILS if QPA ever
  ships one. (e) The **free-form draw-a-surface canvas is deferred** — v1 ships three
  presets (input-side); surfaces are an input method, not a new compute kind. **Recount the
  class table** (`tests/release/test_oracle_classes.py` drives the numbers — run collection,
  paste the LIVE counts, re-run to green; do NOT guess an at-authoring number given the
  mid-merge-train drift).

- [ ] **Step 3: README.** One features line: "marked surfaces → ideal triangulations →
  gentle Jacobian algebras (Fomin–Shapiro–Thurston / Labardini / ABCP), with flip ↔
  cluster mutation certified per instance — draw or pick a surface and get the algebra, a
  no-code input method absent from QPA (unpunctured-with-boundary v1)."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/families tests/surfaces -q` (deep, touched dirs) — note the surface
  batteries live under `tests/families/`; `... -m pytest tests/invariants tests/gui
  tests/webapp -q -m fast`; `... -m pytest tests/qpa -q -m qpa`;
  `... -m pytest tests/release tests/citations -q`, plus a citation-presence check
  (`fomin_shapiro_thurston`/`abcp`/`fomin_zelevinsky_ca1` resolve; the surface block carries
  `fomin_shapiro_thurston`) — all green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-48 surfaces oracle rows + honest scope (v1 unpunctured-with-boundary, flip<->mutation quiver-level, AG not complete, no QPA surface, free-form canvas deferred) + citations + recounted classes"
```

---

## Acceptance (Plan-48 definition of done)

1. `MarkedSurface`, `Triangulation` (+ `fan_triangulation`/`annulus_triangulation`/
   `once_punctured_torus`), `quiver_of`/`potential_of`/`jacobian_of`, `flip`/
   `certify_flip_mutation`, and `surface_block` all public in `quiverlab.surfaces` (and
   top-level exported), all loudly-validated, all self-certified (arc-count identity +
   arc-adjacency on every constructor; finiteness + gentle on every Jacobian).
2. The **ideal-arc count `n = 6g - 6 + 3(b + p) + Σ k_i`** is derived (Task 1), pinned on
   every constructor (disc→`A_n`, annulus→`n+m`, once-punctured-torus→3), and the FST
   admissibility exclusion list (monogon/digon/triangle, small punctured spheres,
   once-punctured monogon) refuses loudly.
3. The **angle→arrow orientation is arbitrated** by the disc oracle: the fan of the
   `(n+3)`-gon yields the linear `A_n` quiver `1→2→…→n` EXACTLY (flipped once and documented
   if the anticlockwise convention came out reversed); the 2-cycle reduction is part of
   `quiver_of`.
4. **Cross-subsystem tie (P44+P46+P48):** `is_gentle(jacobian_of(T))` is `True` across the
   disc/annulus/2-holed surface zoo (ABCP/LFS); the hexagon-with-internal-triangle Jacobian
   is the triangle algebra of **dim 6** (P44's pin); the disc fan gives `kA_n`; a small
   annulus is pinned as an acyclic affine-`Ã` gentle algebra.
5. **Flip ↔ mutation certified per instance:** on every interior arc of the disc/annulus
   zoo, `quiver_of(flip(T, a))` matches the exact Fomin–Zelevinsky matrix mutation at
   vertex `a`; flip is an involution; DWZ potential right-equivalence is deferred to a named
   successor (recorded).
6. **Honest scope enforced:** punctured surfaces, closed surfaces, and self-folded
   triangles refuse loudly (`quiver_of` naming P48.1); the once-punctured torus is the
   pinned loud-refusal oracle; the AG invariant is labelled derived-but-not-complete.
7. **GUI/no-code (input-side):** three drawable presets (disc fan `A_3`, annulus `C(2,2)`,
   hexagon-with-internal-triangle) built via the library and served through the canvas
   dropdown; the surface constructors are catalog-listed and skipped in the webapp family
   form (the `zoo`/non-scalar precedent); the free-form draw-a-surface canvas is the named
   successor. No new compute kind — the produced algebra flows through every existing kind.
8. QPA `IsGentleAlgebra`/`IsSpecialBiserialAlgebra` parity green live (`-m qpa`) on the
   surface algebras, with the standing `NamesGVars()` guard that FAILS if QPA ever ships a
   surface/triangulation surface.
9. `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the five
   honest-scope entries; `fomin_shapiro_thurston`/`abcp`/`fomin_zelevinsky_ca1` citations
   added and BibTeX-verified; README line added; deep (touched dirs) + fast + qpa + release
   + citations suites green. No dependency taken on P42/P43/P45/P47/P49 (merged
   independently to `dev`).
