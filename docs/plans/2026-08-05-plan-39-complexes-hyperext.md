# Plan 39: Complexes, Chain Maps, Cones, Hyper-Ext — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A first-class bounded-chain-complex layer over `mod A`: validated
complexes and chain maps, shift/truncation/homology, mapping cones and
triangles, the derived-iso test, and `Hom_{D^b}(X, Y[n])` (hyper-Ext) via a
certified projective model — the substrate P42 (spectral sequences) and P43
(derived category) build on. C8's "user-facing complex/chain-map API" item.

**Architecture:** New `src/quiverlab/modules/complexes.py` (there is NO
existing generic complex type — 2026-08-05 exploration). Convention pinned
**homological**: `d_n: C_n -> C_{n-1}`, matrices rows=target/cols=source
(exactly `modules/resolution.py::minimal_resolution`'s `dmats` layout;
cohomological indexing is presentation-only: `C^n := C_{-n}`). Terms are
`Module`s over one algebra and one side; differentials validated as
`ModuleHom`-style intertwiners with `d∘d = 0` loud at construction.
Homology dims via the rank formula (`dim H_n = dim ker d_n − rank d_{n+1}`),
homology *modules* via `radtopsoc.submodule`+`quotient` on cycle/boundary
columns. The projective model of a complex is a **twisted totalization**:
termwise `minimal_resolution`s glued by a degreewise lift-solve
(`solve_columns` + `fields.linalg.reduce_mod_nullspace`, canonical — the
Chouhy–Solotar correction-solve precedent), **self-certified**: `d²=0`
asserted exactly AND quasi-isomorphism certified by cone acyclicity, so
correctness never depends on trusting the construction.

**Tech Stack:** `modules/linalg_mod.py` primitives, `ModuleHom` (Plan 37 —
this plan REQUIRES P37 merged), `fields.linalg.reduce_mod_nullspace`.
No floats in `src/`.

## Global Constraints

- Python is always `.venv/bin/python`; tests via
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- Tests live in `tests/modules/test_complexes_*.py` → **deep** bucket
  automatically (`_DEEP_DIRS` includes "modules"); do NOT create a new
  top-level tests dir.
- Convention hazard (from the exploration): module projective resolutions
  are homological (`d_n: Q_n→Q_{n-1}`), injective/Ext machinery is
  cohomological. `ChainComplex` is homological-only; every cohomological
  consumer re-indexes at ITS boundary. State this in the module docstring.
- All validation loud (`QuiverlabError`); `check=False` fast paths only for
  internally-constructed data (cones, models) whose certificates are
  asserted separately.
- Plan-32 markers: identity/certificate tests `oracle_selfcert`; QPA
  comparisons in `tests/qpa/`; the stalk⇒module-Ext battery is
  `oracle_crossengine` (two independent routes to the same numbers).
- GUI story: **explicitly deferred to P42/P43** (this plan is their
  substrate; the metaplan assigns the user-facing SS/derived surfaces
  there). Record this deferral in `docs/verification.md`'s Plan-39 note and
  in the plan's final commit message — the standing rule is honored by the
  named successor, not silently skipped.
- Conventional commits; green at every commit; branch
  `plan-39-complexes` off `dev` (AFTER P37 has merged).

---

### Task 1: `ChainComplex` — validated bounded complexes

**Files:**
- Create: `src/quiverlab/modules/complexes.py`
- Test: `tests/modules/test_complexes_core.py`

**Interfaces:**
- Consumes: `Module`, `ModuleHom` (P37 `modules/morphism.py`),
  `linalg_mod` (`mat_rank`, `kernel_columns`, `matmul`),
  `radtopsoc.submodule/quotient`,
  `modules/resolution.py::minimal_resolution(M, length, ...) -> (terms, dmats)`
  (`terms[n]._term_data.module`, dmats rows=target).
- Produces:
  ```python
  class ChainComplex:
      def __init__(self, terms: dict[int, Module], dmats: dict[int, list], check=True)
          # terms: degree -> Module (bounded, missing degrees = zero module);
          # dmats[n] is d_n: terms[n] -> terms[n-1], rows=target.
          # check: every d_n is a module map (ModuleHom validation) and
          # every composite d_{n} . d_{n+1} == 0 -- QuiverlabError otherwise.
      algebra, side, domain                       # from the terms; must agree
      def degrees(self) -> list[int]              # sorted support
      def term(self, n) -> Module                 # zero module off-support
      def differential(self, n) -> ModuleHom
      def shift(self, k) -> "ChainComplex"        # X[k]_n = X_{n-k}, d -> (-1)^k d
      def truncate(self, lo, hi) -> "ChainComplex"  # brutal truncation
      def homology_dims(self, lo=None, hi=None) -> dict[int, int]
      def homology(self, n) -> Module             # ker/im as a Module
      def is_acyclic(self) -> bool
      def total_dim(self) -> int
      @classmethod
      def stalk(cls, M, degree=0) -> "ChainComplex"
      @classmethod
      def from_projective_resolution(cls, M, length) -> "ChainComplex"
          # P_length -> ... -> P_0 (degrees length..0); homology = M at 0
      def is_perfect(self) -> bool                # every term projective:
          # decompose-free check: term ≅ (+) P_v read off from construction
          # provenance when available, else identify via P37 is_direct_summand
          # against projectives -- if undecidable, QuiverlabError (loud).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_complexes_core.py
"""ChainComplex core: validation, shift/truncate identities, homology.
Self-certifying (d.d=0 gates, rank identities) + cross-engine (resolution
round-trip reproduces the module)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import ChainComplex

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_stalk_homology_is_the_module():
    A = _a3()
    S1 = A.simple(1)
    X = ChainComplex.stalk(S1, degree=0)
    assert X.homology_dims() == {0: 1}
    assert X.homology(0).dimension_vector() == S1.dimension_vector()


def test_resolution_roundtrip_homology_concentrated_in_zero():
    A = _a3()
    S1 = A.simple(1)
    X = ChainComplex.from_projective_resolution(S1, length=4)
    h = X.homology_dims()
    assert h.get(0) == S1.dim
    assert all(v == 0 for n, v in h.items() if n != 0)


def test_dd_nonzero_is_refused():
    A = _a3()
    P1 = A.projective(1)
    idm = P1.identity_hom().matrix
    with pytest.raises(QuiverlabError, match="d.*d|square"):
        ChainComplex({1: P1, 0: P1, -1: P1}, {1: idm, 0: idm})


def test_non_module_map_differential_refused():
    A = _a3()
    P1, S1 = A.projective(1), A.simple(1)
    bad = [[1, 1]]                        # not an intertwiner (P37 pin)
    with pytest.raises(QuiverlabError, match="module map"):
        ChainComplex({1: P1, 0: S1}, {1: bad})


def test_shift_moves_degrees_and_signs():
    A = _a3()
    S1 = A.simple(1)
    X = ChainComplex.from_projective_resolution(S1, length=2)
    Y = X.shift(3)
    assert Y.homology_dims().get(3) == S1.dim
    Z = Y.shift(-3)
    assert Z.homology_dims() == X.homology_dims()
    # d.d = 0 still validated after odd shift (sign flip is consistent)
    ChainComplex({n: Y.term(n) for n in Y.degrees()},
                 {n: Y.differential(n).matrix for n in Y.degrees()
                  if Y.term(n - 1).dim and Y.term(n).dim})


def test_truncate_brutal():
    A = _a3()
    S1 = A.simple(1)
    X = ChainComplex.from_projective_resolution(S1, length=4)
    T = X.truncate(0, 1)
    assert set(T.degrees()) <= {0, 1}


def test_mixed_sides_refused():
    A = _a3()
    R = A.simple(1)
    L = A.simple(1, side="left")
    with pytest.raises(QuiverlabError):
        ChainComplex({0: R, 1: L}, {1: [[0]]})
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_complexes_core.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.modules.complexes`

- [ ] **Step 3: Implement** — `complexes.py` per the Interfaces block.
  Implementation notes (bind to reality, read the named files first):
  - Zero modules off-support: build via `Module(A, 0, {}, ...)` — read how
    the zero module is represented elsewhere (`decompose` handles dim 0;
    mirror it; if no precedent, a `Module` with `dim=0` and empty action —
    verify `check_module` accepts it, else add the trivial guard there).
  - `homology_dims`: `dim ker d_n − rank d_{n+1}` (absent differentials =
    zero maps). `homology(n)`: cycle columns via `kernel_columns(d_n)`
    (whole term if `d_n` absent), boundary columns = columns of `d_{n+1}`;
    `submodule(term, cycles)` then `quotient` by boundaries expressed in
    the submodule — reuse the exact pattern of
    `complex_reps._reps_from_complex` (read it; it is the tested homology-
    of-two-maps primitive) rather than reinventing coordinates.
  - `shift(k)`: degree reindex + multiply every differential by `(-1)^k`
    (domain-exact scalar).
  - `from_projective_resolution`: consume `minimal_resolution(M, length)`;
    `terms[n].module`; note `dmats[0]` is the augmentation `Q_0 -> M` —
    EXCLUDE it (the complex is `Q_*`, homology at 0 is M by exactness);
    read the exact indexing in `modules/resolution.py:100-130` first.
  - `is_perfect`: constructor provenance flag (`from_projective_resolution`
    sets it) OR the P37 `identify_standard`/`is_direct_summand` route;
    loud when undecidable.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/complexes.py tests/modules/test_complexes_core.py
git commit -m "feat(complexes): validated bounded ChainComplex -- homological convention, shift/truncate/homology"
```

---

### Task 2: `ChainMap`, mapping cone, triangles, derived-iso test

**Files:**
- Modify: `src/quiverlab/modules/complexes.py`
- Test: `tests/modules/test_complexes_cone.py`

**Interfaces:**
- Produces:
  ```python
  class ChainMap:
      def __init__(self, src: ChainComplex, tgt: ChainComplex,
                   components: dict[int, list], check=True)
          # components[n]: matrix src.term(n) -> tgt.term(n); every square
          # d^Y_n f_n == f_{n-1} d^X_n validated (loud).
      def cone(self) -> ChainComplex
          # cone(f)_n = X_{n-1} (+) Y_n, d = [[-d_X, 0], [f, d_Y]]
          # (d^2 = 0 automatic -- but STILL constructed with check=True in
          #  tests once per shape as the self-certificate)
      def is_quasi_iso(self) -> bool              # cone acyclic
      def triangle(self) -> tuple                  # (X, Y, cone, maps):
          # inclusion Y -> cone, projection cone -> X[1]
  def identity_chain_map(X) -> ChainMap
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_complexes_cone.py
"""Cones and triangles. Self-certifying: cone d.d=0 under full validation;
LES rank identities; quasi-iso <=> cone acyclic on known cases."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import (ChainComplex, ChainMap,
                                         identity_chain_map)

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def test_cone_of_identity_is_acyclic():
    A = _a2()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=3)
    f = identity_chain_map(X)
    C = f.cone()
    assert C.is_acyclic()
    assert f.is_quasi_iso()


def test_cone_of_zero_map_is_direct_sum_shift():
    A = _a2()
    S1 = A.simple(1)
    X = ChainComplex.stalk(S1, 0)
    Y = ChainComplex.stalk(S1, 0)
    z = ChainMap(X, Y, {0: [[0]]})
    C = z.cone()
    # cone(0) = X[1] (+) Y: homology in degrees 1 and 0
    assert C.homology_dims().get(1) == 1 and C.homology_dims().get(0) == 1


def test_resolution_augmentation_is_quasi_iso():
    # the augmentation Q_* -> M (stalk) is the canonical quasi-iso
    A = _a2()
    S1 = A.simple(1)
    Q = ChainComplex.from_projective_resolution(S1, length=3)
    M = ChainComplex.stalk(S1, 0)
    from quiverlab.modules.resolution import projective_cover
    Q0, d0, _ = projective_cover(S1)
    aug = ChainMap(Q, M, {0: d0})
    assert aug.is_quasi_iso()


def test_bad_square_refused():
    A = _a2()
    P1 = A.projective(1)
    X = ChainComplex.from_projective_resolution(A.simple(1), length=2)
    with pytest.raises(QuiverlabError, match="square|chain map"):
        ChainMap(X, X, {n: [[1] * X.term(n).dim for _ in range(X.term(n).dim)]
                        for n in X.degrees() if X.term(n).dim})


def test_triangle_euler_characteristic():
    # LES consequence: chi(cone) = chi(X[1]) + chi(Y) = -chi(X) + chi(Y)
    A = _a2()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=2)
    Y = ChainComplex.stalk(A.simple(2), 0)
    z = ChainMap(X, Y, {0: [[0] * X.term(0).dim]})
    C = z.cone()
    chi = lambda Z: sum((-1) ** n * d for n, d in Z.homology_dims().items())
    assert chi(C) == -chi(X) + chi(Y)
```

- [ ] **Step 2: Run to verify failure** — `ImportError: ChainMap`

- [ ] **Step 3: Implement** per the Interfaces block. The cone assembles
  block matrices with the sign on `-d_X`; the triangle maps are the block
  inclusion/projection. `test_resolution_augmentation_is_quasi_iso` fixes
  the augmentation convention — if `projective_cover`'s `d0` shape
  disagrees (read it), adapt the test comment, not the convention.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/complexes.py tests/modules/test_complexes_cone.py
git commit -m "feat(complexes): ChainMap + mapping cone + triangles + quasi-iso via cone acyclicity"
```

---

### Task 3: Hom total complex — `hyper_hom_dims` for perfect sources

**Files:**
- Modify: `src/quiverlab/modules/complexes.py`
- Test: `tests/modules/test_complexes_hyperhom.py`

**Interfaces:**
- Consumes: `hom.py::hom_space` (the single-Hom primitive),
  `ext.py::_delta_matrix` (the φ ↦ φ∘d induced-map pattern — read it; the
  Hom-total differential generalizes it with the Y-side post-composition
  term and the Koszul sign).
- Produces:
  ```python
  def hyper_hom_dims(X: ChainComplex, Y: ChainComplex, lo: int, hi: int) -> dict[int, int]
      # dim Hom_{K(A)}(X, Y[-n])? NO -- pin the meaning precisely:
      # returns {n: dim H^n(Hom^•(X, Y))} for n in lo..hi, where
      # Hom^n = (+)_{q-p=n... } -- SIGN/INDEX CONVENTION (document verbatim):
      # Hom^n(X, Y) = prod_p Hom_A(X_p, Y_{p-n}),
      # (delta f)_p = d^Y . f_p - (-1)^n f_{p-1}? -- use the standard
      # convention (Weibel 2.7.4): (delta f) = d_Y f - (-1)^n f d_X.
      # For X PERFECT (complex of projectives) H^n(Hom^•) computes
      # Hom_{D^b}(X, Y[n]).
      # QuiverlabError if X is not certified perfect (Task-4 lifts this).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_complexes_hyperhom.py
"""Hyper-Hom of a perfect complex: pinned degreewise against module
Ext/Hom -- two independent engines computing the same numbers."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.complexes import ChainComplex, hyper_hom_dims

pytestmark = pytest.mark.oracle_crossengine


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_projective_stalk_hyper_hom_is_hom():
    A = _a3()
    P1 = A.projective(1)
    for v in (1, 2, 3):
        N = ChainComplex.stalk(A.simple(v), 0)
        X = ChainComplex.stalk(P1, 0)
        hh = hyper_hom_dims(X, N, 0, 3)
        assert hh[0] == A.hom(P1, A.simple(v))
        assert all(hh[n] == 0 for n in (1, 2, 3))


def test_resolution_hyper_hom_computes_ext():
    # THE pin: X = proj resolution of M as a perfect complex, Y = stalk N
    # => H^n(Hom(X, Y)) = Ext^n(M, N) for n within the resolution window.
    A = _a3()
    for v in (1, 2, 3):
        M, N = A.simple(v), A.simple(1)
        X = ChainComplex.from_projective_resolution(M, length=5)
        got = hyper_hom_dims(X, ChainComplex.stalk(N, 0), 0, 4)
        for n in range(0, 5):
            assert got[n] == A.ext(M, N, n), (v, n)


def test_shift_bookkeeping():
    A = _a3()
    M, N = A.simple(2), A.simple(1)
    X = ChainComplex.from_projective_resolution(M, length=4)
    Y = ChainComplex.stalk(N, 0)
    base = hyper_hom_dims(X, Y, 0, 3)
    shifted = hyper_hom_dims(X, Y.shift(1), -1, 2)
    assert all(shifted[n - 1] == base[n] for n in range(0, 3))


def test_nonperfect_source_refused():
    A = _a3()
    X = ChainComplex.stalk(A.simple(2), 0)     # simple: not projective
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError, match="perfect"):
        hyper_hom_dims(X, X, 0, 1)
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement.** Per total degree n build the product space
  from the `hom_space(X_p, Y_{p-n})` bases; the differential's two blocks
  are post-composition with `d_Y` and pre-composition with `d_X` (the
  `_delta_matrix` coordinate pattern, plus the `(-1)^n` sign); dims by the
  rank formula. Write the chosen sign convention VERBATIM in the docstring
  and cite Weibel 2.7.4 — `test_resolution_hyper_hom_computes_ext` is the
  arbiter: if it fails on odd degrees, the sign is wrong, flip it once,
  document, never special-case per degree.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/complexes.py tests/modules/test_complexes_hyperhom.py
git commit -m "feat(complexes): Hom total complex -- hyper-Hom dims for perfect sources, Ext-pinned signs"
```

---

### Task 4: certified projective model + general hyper-Ext

**Files:**
- Modify: `src/quiverlab/modules/complexes.py`
- Test: `tests/modules/test_complexes_model.py`

**Interfaces:**
- Consumes: `minimal_resolution`, `solve_columns`,
  `fields.linalg.reduce_mod_nullspace` (read its exact signature in
  `src/quiverlab/fields/linalg.py` — the CS correction-solve precedent),
  `ChainMap.is_quasi_iso` (Task 2 — the certificate).
- Produces:
  ```python
  def projective_model(X: ChainComplex, length: int) -> tuple[ChainComplex, ChainMap]
      # (P, eps) with P perfect (degrees up to max(X)+length), eps: P -> X
      # a chain map, CERTIFIED: eps.is_quasi_iso() asserted before return
      # (QuiverlabError "model failed certification" otherwise -- never
      # return an uncertified model).
      # Construction: twisted totalization -- termwise minimal resolutions
      # P_{p,q} of X_p; vertical d_v = resolution differentials; horizontal
      # d_h lifted degreewise (solve_columns, canonicalized via
      # reduce_mod_nullspace); higher corrections d_k solved so that the
      # total differential squares to zero, degree by degree (each equation
      # is solvable because the columns are projective and the vertical
      # complexes are exact in positive degrees; the solve is the same
      # shape as resolutions_cs's _d_general correction solve).
  def hyper_ext_dims(X, Y, lo, hi, length=8) -> dict[int, int]
      # projective_model(X, length) then hyper_hom_dims -- the general case.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_complexes_model.py
"""Projective models: certified quasi-iso by construction; hyper-Ext on
stalks reproduces module Ext (cross-engine); two-term complex LES pin."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.complexes import (ChainComplex, ChainMap,
                                         hyper_ext_dims, projective_model)

pytestmark = pytest.mark.oracle_crossengine


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_model_of_stalk_is_certified_and_computes_ext():
    A = _a3()
    for v in (1, 2, 3):
        M = A.simple(v)
        X = ChainComplex.stalk(M, 0)
        P, eps = projective_model(X, length=5)
        assert P.is_perfect() and eps.is_quasi_iso()
        got = hyper_ext_dims(X, ChainComplex.stalk(A.simple(1), 0), 0, 4,
                             length=6)
        for n in range(0, 5):
            assert got[n] == A.ext(M, A.simple(1), n), (v, n)


def test_model_of_two_term_complex():
    # X = [P1 --f--> S1] (f the cover, degrees 1,0): quasi-iso to the stalk
    # of ker f shifted -- hyper-Ext must match Ext of rad P1 with a shift.
    A = _a3()
    P1, S1 = A.projective(1), A.simple(1)
    from quiverlab.modules.resolution import projective_cover
    Q0, d0, _ = projective_cover(S1)          # Q0 = P1 here
    X = ChainComplex({1: Q0, 0: S1}, {1: d0})  # NOTE: check orientation --
    # this complex has homology ker d0 = rad P1 in degree 1 and 0 in degree 0
    P, eps = projective_model(X, length=4)
    assert eps.is_quasi_iso()
    R = Q0.radical()
    got = hyper_ext_dims(X, ChainComplex.stalk(A.simple(2), 0), -1, 3)
    for n in range(0, 3):
        # Hom_{D}(R[1], N[n]) = Ext^{n+1}(R, N)... verify the exact shift
        # ARITHMETIC against first principles in the implementation session
        # and pin the passing identity with a comment deriving it.
        assert got[n] == A.ext(R, A.simple(2), n + 1) or True  # sharpen: see Step 3


def test_model_length_window_is_honest():
    # dims past the model window must not be silently wrong: the function
    # raises when [lo, hi] exceeds what `length` certifies.
    A = _a3()
    X = ChainComplex.stalk(A.simple(1), 0)
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError, match="length|window"):
        hyper_ext_dims(X, X, 0, 30, length=3)
```

**Sharpening rule for `test_model_of_two_term_complex`:** the `or True`
placeholder MUST be removed during implementation — derive the correct
shift identity on paper (the complex `[P1 -> S1]` with the cover map is
quasi-isomorphic to `rad P1` placed in degree 1), write the derivation as
a comment, assert the exact equality. A test that cannot be sharpened
means the model is wrong — stop and debug, do not weaken.

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** per the Interfaces block; sharpen the
  two-term test. The window guard: a model of length L certifies
  hyper-Ext degrees `hi <= L - 1` (one degree of slack for the rank
  formula's `d_{n+1}`) — derive and enforce, loud message naming `length`.

- [ ] **Step 4: Run tests** — Expected: PASS with the placeholder removed

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/complexes.py tests/modules/test_complexes_model.py
git commit -m "feat(complexes): certified projective models (twisted totalization) + general hyper-Ext"
```

---

### Task 5: QPA crosscheck battery

**Files:**
- Modify: `src/quiverlab/qpa/scripts.py`, `src/quiverlab/qpa/crosscheck.py`
- Test: `tests/qpa/test_complexes_qpa.py`

**Interfaces:**
- Consumes: QPA Ch. 10 (`Complex`/`FiniteComplex` constructors,
  `MappingCone`, `HomologyOfComplex`, `Shift`) — probe live for the exact
  constructor signature QPA 1.37 exposes (the `NamesGVars()` sweep
  precedent in `tests/qpa/test_products_qpa.py`); `scripts.module_decl`.
- Produces: live crosschecks — (a) homology dims of a small explicit
  complex (both systems build the SAME complex over kA3-with-relation);
  (b) cone homology dims of a chain map both systems build;
  (c) `Shift` bookkeeping. Standard skipif header. If QPA's complex
  constructors turn out not to be scriptable through libgap (streams/
  categories oddities), fall back to comparing OUR cone/homology against
  QPA's `ProjectiveResolutionOfComplex`-free primitives that DO script,
  and record exactly what was compared in the test docstring — honest
  coverage, never a silent skip.

- [ ] **Step 1: Probe + write the battery** (concrete tests in the
  `test_tor_qpa.py` style)
- [ ] **Step 2: Run live** (`-m qpa`); Expected: PASS
- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/qpa/ tests/qpa/test_complexes_qpa.py
git commit -m "test(qpa): complexes battery -- homology/cone/shift vs QPA Ch.10 (live)"
```

---

### Task 6: public surface, verification page, gate

**Files:**
- Modify: `src/quiverlab/core/algebra.py` or `src/quiverlab/__init__.py`
  (export: `from quiverlab import ChainComplex, ChainMap` — read
  `src/quiverlab/__init__.py` for the export idiom and add the two names;
  plus `Algebra.chain_complex(terms, dmats)` convenience wrapper)
- Modify: `docs/verification.md`, `docs/plans/2026-08-05-metaplan-v0.2.0.md`
  (tick the P39 card's delivery note + the GUI-deferral pointer at P42/P43)
- Test: extend `tests/modules/test_complexes_core.py` with one top-level-import test

- [ ] **Step 1:** exports + the one test
  (`from quiverlab import ChainComplex` works; `A.chain_complex(...)`
  round-trips).
- [ ] **Step 2:** verification page: Plan-39 subsystem row
  (`modules/complexes.py` | oracles: self-cert d²=0/cone + cross-engine
  Ext pins + live QPA), the honest note "GUI surface arrives with P42/P43",
  recounted class counts (`tests/release/test_oracle_classes.py` green).
- [ ] **Step 3:** Full gate:
  `... -m pytest tests/modules -q` (deep, the touched files),
  `... -m pytest -q -m fast`, `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release -q` — all green.
- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "feat(complexes): public exports + verification page; GUI story deferred to P42/P43 (recorded)"
```

---

## Acceptance (Plan-39 definition of done)

1. `ChainComplex`/`ChainMap`/cone/triangle/quasi-iso/hyper-Hom/
   `projective_model`/`hyper_ext_dims` public, validated, certified
   (models NEVER returned uncertified).
2. Stalk hyper-Ext ≡ module Ext degreewise (cross-engine battery);
   augmentation quasi-iso pin; Euler-characteristic triangle pin.
3. Live QPA complex battery green; honest docstring if constructor
   scripting forced the fallback comparison.
4. Sign conventions documented verbatim in docstrings (homological d,
   Weibel Hom-complex sign), each pinned by a named arbiter test.
5. Verification page updated with the GUI-deferral note; all gates green.
   P42/P43 can consume the API exactly as specified in Interfaces blocks.
