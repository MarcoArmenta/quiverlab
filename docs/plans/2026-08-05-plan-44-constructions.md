# Plan 44: C7 Tilting & Constructions — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The construction toolkit a representation theorist reaches for once modules
are first-class — **minimal left/right add(M)-approximations**; **`is_tilting`/`is_cotilting`
+ Bongartz completion**; the crux **basic-ization / Gabriel-quiver recovery of a
structure-constant algebra** (which finally makes `End(M)` / `End(T)` readable as a
`kQ/I` — tilted algebras, derived-equivalent algebras); **one-point extensions
`A[M]`** (feeding P42's Happel-LES oracle); **finite repetitive-algebra slices**; and a
**Jacobian-algebra constructor from a quiver-with-potential `(Q, W)`** with a
per-instance finite-dimensionality certificate (feeding P48). Every construction is
CERTIFIED per instance (dimension identity + a checkable structural oracle), and every
honesty edge — the basic-ization char caveat, an infinite Jacobian algebra, the
finite-slice-only scope of the repetitive algebra — refuses loudly, never a silent
wrong answer.

**Architecture:** Four new modules, each a thin exact-linear-algebra / mini-Gröbner
layer over EXISTING primitives (no new math engines):

- `src/quiverlab/modules/approximations.py` — minimal add(M)-approximations over
  `hom_basis` (P37) + `linalg_mod.solve_columns`/`independent_modulo` + the
  `end_algebra` (P37) radical.
- `src/quiverlab/modules/tilting.py` — `is_tilting_module`/`is_cotilting_module` over
  `ProjectiveResolution.pd()` + `ext_dims` + `decompose` (P30) + `is_direct_summand`
  (P37); `bongartz_completion` over `yoneda.baer_extension` (the existing Ext¹
  realizer) + `direct_sum`.
- `src/quiverlab/core/basic.py` — the crux. `primitive_idempotents` /
  `basic_algebra` / `gabriel_quiver` / `presented_form` via the trace-form radical
  (the `decompose.py::_trace_form_rank` / `endomorphism.py::_with_radical_labels`
  idiom, the ONLY exact radical route that exists), a semisimple quotient, Newton-style
  idempotent lifting, and the `TrivialExtension.extract_relations` length-lex
  mini-Gröbner kernel enumeration for the relations. Thin `Algebra` wrappers.
- `src/quiverlab/families/one_point.py`, `.../repetitive.py`, `.../jacobian.py` — the
  three named constructions, each presented via `Quiver.algebra` with a per-instance
  dimension certificate, mirroring `families/trivial_extension.py::_presented_trivial_extension`.

The whole thing composes on the P37 hub (`ModuleHom`/`direct_sum`/`is_direct_summand`/
`end_algebra`/`ses`/`hom_basis`), the P30 `decompose`, the P38 forms (`euler_form`,
`cartan_matrix`), and the presented-algebra backbone
(`Quiver.algebra`, `combinat/relations.py::parse_relations`,
`groebner/certificate.py`'s `NotFiniteDimensionalError`).

**Tech Stack:** pure exact linear algebra over `Domain` (`modules/linalg_mod`,
`fields/linalg`); sympy only where `decompose` already uses it (univariate factoring in
`decompose._factor_min_poly`). No floats in `src/` (AST-gated by
`tests/test_no_floats.py`).

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P37 (categorical glue) and P38 (forms) are MERGED to `dev` and are hard
  prerequisites** — this plan consumes `modules/morphism.py`
  (`ModuleHom`/`hom_basis`/`direct_sum`/`is_direct_summand`), `modules/ses.py`,
  `modules/endomorphism.py::end_algebra`, `modules/yoneda.py::baer_extension`,
  `modules/decompose.py`, and `invariants/forms.py::euler_form`. **This plan has NO
  dependency on P39, P40, or P41** — it can execute in parallel with the rest of Wave-2
  / Wave-3. Branch `plan-44-constructions` off `dev` (do NOT wait on P39/P40/P41).
- Buckets are auto-assigned by directory (`tests/conftest.py`): `tests/modules/` and
  `tests/families/` → **deep**; `tests/core/` and `tests/invariants/` → **fast**;
  `tests/qpa/` → **qpa**; `tests/webapp/` and `tests/gui/` → **fast**. Run new deep
  tests by path during development, finish with a `-m deep` spot-run of the touched
  files. **The Task-C basic-recovery tests live in `tests/families/` (deep)** so the
  heavier construct-and-factor batteries share the deep budget with the constructions
  that consume them.
- **The `decompose` char-caveat is load-bearing for Task C and the tilting count.**
  `decompose`/`is_isomorphic` and the trace-form radical (`_trace_form_rank`,
  `_with_radical_labels`) are only rigorous over **char 0 or char > dim**. So
  `primitive_idempotents`/`basic_algebra`/`gabriel_quiver`/`presented_form`, and the
  `is_tilting_module` summand COUNT, run their batteries over **QQ** or a **large prime
  GF(32003)**. Over `char ≤ dim` these functions inherit the loud
  `QuiverlabError` refusal unchanged — **never a silent wrong idempotent set or count**.
- **The honesty pattern is mandatory** (the `GlobalDimension` precedent,
  `modules/ext.py:63-82`, and the `TrivialExtension` per-instance certificate,
  `families/trivial_extension.py:298-308`): every construction is CERTIFIED per instance
  (a dimension identity plus a structural oracle) or REFUSES loudly. A Jacobian algebra
  the presented backbone cannot certify finite raises `NotFiniteDimensionalError`; a
  basic-ization off char-scope raises `QuiverlabError`; a repetitive *slice* is finite
  and certified, the full repetitive algebra is out of scope and stated so.
- **Composition is left-to-right like paths** (`f.then(g)`; `a*b` = first `a` then `b`).
  Never overload `*` on morphisms. `_assert_comparable` guards every cross-module/
  cross-side/cross-algebra call.
- Plan-32 markers: dimension/round-trip/idempotent certificates and the
  approximation/tilting self-identities = `oracle_selfcert`; ASS-VI tilting worked
  examples, the `M₂(k)→k` / `kA_n` Gabriel-recovery pins, the one-point Cartan-block +
  `pd(S_ω)=pd_A(M)+1` identity, and the hand-derived Jacobian dimensions =
  `oracle_literature`; the preprojective-as-Jacobian and `End(⊕P_v)`-recovers-`kA_n`
  agreements = `oracle_crossengine`; QPA comparisons live in `tests/qpa/` (bucket = the
  class, never double-marked).
- **Mid-merge-train counts:** v0.2.0 lands ~15 subplans in overlapping waves, so the
  absolute suite counts drift between this plan's authoring and its merge. **Task H
  recounts the oracle-class table at merge time by running
  `tests/release/test_oracle_classes.py`** (paste-the-live-numbers, never a
  guessed-at-authoring count) and only claims the deltas this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class
  table green) and the README line. Conventional commits; green at every commit.

---

### Task A: minimal left/right add(M)-approximations

**Files:**
- Create: `src/quiverlab/modules/approximations.py`
- Test: `tests/modules/test_approximations.py`

**Interfaces:**
- Consumes: `morphism.py::hom_basis(M, C) -> list[ModuleHom]`,
  `morphism.py::direct_sum(*mods) -> (D, incls, projs)`, `morphism.py::ModuleHom`,
  `endomorphism.py::end_algebra(M)` (P37 — its structure-constant `End_A(M)` and the
  trace-form radical labels), `linalg_mod::solve_columns`/`independent_modulo`/
  `cols_to_matrix`, `hom.py::_assert_comparable`, `duality.py::dualize`.
- Produces:
  ```python
  def right_add_approximation(M, C) -> ModuleHom
      # f: M^k ->> "onto Hom(M,-)" C, a minimal right add(M)-approximation:
      #   (approx)  Hom_A(M, M^k) -> Hom_A(M, C) is surjective  (every M->C factors
      #             through f)
      #   (minimal) k = # minimal generators of Hom_A(M,C) as a right End_A(M)-module,
      #             so dropping any summand of M^k breaks surjectivity.
      # Returns the ModuleHom f (src = the k-fold direct sum of M, tgt = C).
  def left_add_approximation(M, C) -> ModuleHom
      # dual: g: C >-> M^k, a minimal left add(M)-approximation, built over the
      # opposite algebra via dualize (D of a right add(DM)-approximation of DC).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_approximations.py
"""Minimal add(M)-approximations (Plan 44 / C7). Self-certifying: the approximation
property is checked by solve_columns surjectivity of Hom(M, M^k) -> Hom(M, C); minimality
by no-proper-summand-restriction. The projective-cover tie is the literature anchor."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.approximations import (left_add_approximation,
                                              right_add_approximation)
from quiverlab.modules.morphism import hom_basis

pytestmark = pytest.mark.oracle_selfcert


def _kA2():
    # QQ so the End(M) trace-form radical (the minimality selector) is rigorous.
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def _factors_every_hom(f, M, C):
    """Every h in Hom(M, C) equals f-after-(some M -> M^k): the approximation property.
    Reduce to a solve over the columns { (h_j composed into f's coordinates) }."""
    from quiverlab.modules import linalg_mod as lm
    dom = C.domain
    homs = hom_basis(M, C)
    # columns spanning im(Hom(M, M^k) -> Hom(M, C)) via f: each basis phi: M -> M^k
    # gives f.then?  no -- phi.then(f): M -> C.  build the induced-map image.
    from quiverlab.modules.morphism import hom_basis as hb
    induced = [phi.then(f) for phi in hb(M, f.src)]        # M -> M^k -> C
    def vec(g):  # g: M -> C, flatten
        return [g.matrix[i][j] for j in range(M.dim) for i in range(C.dim)]
    B = lm.cols_to_matrix([vec(g) for g in induced]) if induced else []
    for h in homs:
        if lm.solve_columns(B, [vec(h)], dom) is None:
            return False
    return True


def test_projective_cover_is_the_min_add_P_approximation():
    # add(P1) contains P(S1) = P1, so the min right add(P1)-approximation of S1 IS the
    # projective cover P1 ->> S1 (ASS I.5). Anchors approximation = cover.
    A = _kA2()
    P1, S1 = A.projective(1), A.simple(1)
    f = right_add_approximation(P1, S1)
    assert f.is_epi() and f.tgt is S1
    K, _ = f.kernel()
    assert K.dimension_vector() == P1.radical().dimension_vector()   # cover: ker in rad


def test_radical_inclusion_is_the_min_add_S2_approximation():
    # Hom(S2, P1) = k (S2 = rad P1 >-> P1); the min right add(S2)-approx of P1 is that
    # mono, NOT surjective (P1/S2 = S1 != 0). A genuine non-epi approximation.
    A = _kA2()
    S2, P1 = A.simple(2), A.projective(1)
    f = right_add_approximation(S2, P1)
    assert f.is_mono() and not f.is_epi() and f.src.dim == 1


@pytest.mark.parametrize("mv, cv", [(1, 1), (1, 2), (2, 1)])
def test_approximation_property_and_minimality_battery(mv, cv):
    A = _kA2()
    pool = {1: A.projective(1), 2: A.projective(2)}
    M, C = pool[mv], A.simple(cv)
    f = right_add_approximation(M, C)
    assert _factors_every_hom(f, M, C)                     # approximation property
    # minimality: dropping any single copy of M breaks surjectivity.
    from quiverlab.modules.morphism import direct_sum
    k = f.src.dim // M.dim if M.dim else 0
    if k >= 2:
        # rebuild f restricted to k-1 copies via the same constructor on a proper
        # submodule of Hom(M,C) generators would fail surjectivity -- the constructor
        # already returns the minimal k, so assert k is minimal by the End-module
        # generator count (checked inside; here we assert k is not inflated).
        assert k == right_add_approximation(M, C).src.dim // M.dim


def test_left_is_dual_of_right():
    A = _kA2()
    M, C = A.projective(1), A.simple(2)
    g = left_add_approximation(M, C)
    assert g.src is C                                      # C >-> M^k
    # every C -> M factors through g (dual approximation property)
    from quiverlab.modules.morphism import hom_basis as hb
    from quiverlab.modules import linalg_mod as lm
    dom = C.domain
    induced = [g.then(psi) for psi in hb(g.tgt, M)]        # C -> M^k -> M
    def vec(h):
        return [h.matrix[i][j] for j in range(C.dim) for i in range(M.dim)]
    B = lm.cols_to_matrix([vec(h) for h in induced]) if induced else []
    for h in hb(C, M):
        assert lm.solve_columns(B, [vec(h)], dom) is not None


def test_cross_algebra_refused():
    A, B = _kA2(), _kA2()
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError):
        right_add_approximation(A.simple(1), B.simple(1))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_approximations.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.modules.approximations`

- [ ] **Step 3: Implement `src/quiverlab/modules/approximations.py`**

```python
"""Minimal left/right add(M)-approximations (Plan 44 / C7).

A right add(M)-approximation f: X -> C (X in add M) means every M -> C factors through
f (Hom(M, -) surjectivity). It is right MINIMAL iff X has no summand mapping to 0 -- the
minimal one is unique, with X = M^k where k = the number of generators of Hom_A(M, C) as
a right End_A(M)-module (= dim Hom(M,C) / Hom(M,C)*rad End(M)). Dual: left add(M)-approx
via D over A^op. Float-free, exact linear algebra over the Domain (ASS I.2, I.5)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import _assert_comparable
from quiverlab.modules.morphism import ModuleHom, direct_sum, hom_basis


def _min_add_generators(M, C):
    """Return the sublist of hom_basis(M, C) that minimally generates Hom_A(M, C) as a
    right End_A(M)-module: a basis modulo Hom(M,C)*rad End(M). Uses the End(M) radical
    (trace form; char 0 / char > dim M) -- inherits its loud refusal off scope."""
    homs = hom_basis(M, C)
    if not homs:
        return homs
    dom = C.domain
    def vec(g):
        return [g.matrix[i][j] for j in range(M.dim) for i in range(C.dim)]
    # right End(M)-action generators: g . rho = rho.then(g) for rho a radical endo of M.
    from quiverlab.modules.endomorphism import _radical_endos   # NEW small helper, below
    rad = _radical_endos(M)                                     # list[ModuleHom], may be []
    sub = []
    for g in homs:
        for rho in rad:
            sub.append(vec(rho.then(g)))                        # in Hom(M,C)*rad End(M)
    sub_cols = sub
    cand = [vec(g) for g in homs]
    keep = lm.independent_modulo(cand, sub_cols, dom)           # minimal generators
    return [homs[i] for i in keep]


def right_add_approximation(M, C):
    _assert_comparable(M, C, "add-approximation")
    gens = _min_add_generators(M, C)
    k = len(gens)
    if k == 0:                                                 # Hom(M, C) = 0: X = 0
        dom = C.domain
        Z = M.with_side(M.side)                                # a zero source; see note
        return ModuleHom(_zero_like(M, 0), C,
                         lm.zeros(C.dim, 0, dom), check=False)
    D, incls, _projs = direct_sum(*([M] * k))
    dom = C.domain
    # f: M^k -> C is [g_0 | ... | g_{k-1}] on the k blocks: f-after-incl_i == g_i.
    cols = [[dom.zero()] * D.dim for _ in range(C.dim)]        # C.dim x D.dim
    for i, g in enumerate(gens):
        # place g's columns into the block incl_i occupies
        for j in range(M.dim):
            src_col = incls[i].matrix                          # M -> D, D.dim x M.dim
            # column of D corresponding to basis j of the i-th M block:
            dcol = [src_col[r][j] for r in range(D.dim)]
            pivot = next(r for r in range(D.dim) if not dom.is_zero(dcol[r]))
            for r in range(C.dim):
                cols[r][pivot] = g.matrix[r][j]
    f = ModuleHom(D, C, cols, check=True)
    return f


def left_add_approximation(M, C):
    _assert_comparable(M, C, "add-approximation")
    from quiverlab.modules.duality import dualize
    DM, DC = dualize(M), dualize(C)
    Df = right_add_approximation(DM, DC)                        # DC's approx: (DM)^k -> DC
    g = dualize_hom(Df)                                         # C -> M^k, see note
    return g
```

**Adjust to reality (Task A):**
- `endomorphism.py` currently exposes `_with_radical_labels`; add a tiny public-ish
  helper `_radical_endos(M) -> list[ModuleHom]` there (or inline in
  `approximations.py`): the radical of `End_A(M)` is the trace-form nullspace
  (`decompose._trace_form_rank`'s Gram matrix, then `lm.kernel_columns`), each nullspace
  vector re-expressed as an endomorphism via the `hom_basis(M, M)` change-of-basis `B`.
  Reuse `endomorphism._structure_constants(M)`'s `B`. **This inherits the `char ≤ dim M`
  loud refusal** (return the same `QuiverlabError` decompose raises) — the batteries run
  over QQ.
- The zero-source construction (`Hom(M,C)=0`) needs a genuine zero `Module`; mirror the
  zero-module idiom in `resolution.py`/`duality.py::_zero_module` (a helper `_zero_like`
  or `duality._zero_module(M.algebra, side=M.side)`); DO NOT hand-roll it.
- `dualize_hom` (dualizing a `ModuleHom`): `D` is contravariant, so
  `dualize(Df): D(DC) -> D((DM)^k)`, i.e. `C -> M^k` after the canonical `D²≅id`. If
  `duality.py` has no map-level `dualize`, add one (transpose the matrix through the
  `dualize`-on-objects basis) — the P41 AR-completion plan may already have added it;
  check `duality.py` first and reuse.
- Minimality certificate (the `oracle_selfcert` anchor): the returned `k` equals
  `len(_min_add_generators(M, C))` BY CONSTRUCTION; the test's "dropping a summand breaks
  surjectivity" is the checkable witness. Write the drop-a-column-then-`solve_columns`
  check as a real sub-assertion if `k >= 2` appears in the tested zoo.

- [ ] **Step 4: Run tests, verify pass**

Run: `... -m pytest tests/modules/test_approximations.py -v` — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/approximations.py src/quiverlab/modules/endomorphism.py \
        tests/modules/test_approximations.py
git commit -m "feat(modules): minimal left/right add(M)-approximations (End-module minimal generators, projective-cover-anchored)"
```

---

### Task B: `is_tilting` / `is_cotilting` + Bongartz completion

**Files:**
- Create: `src/quiverlab/modules/tilting.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.is_tilting_module`, `Algebra.bongartz_completion` thin wrappers beside `global_dimension`, `core/algebra.py:448`)
- Test: `tests/modules/test_tilting.py`

**Interfaces:**
- Consumes: `Module.projective_resolution(bound).pd()` (P05/P30 honesty pattern),
  `ext.py::ext_dims(A, T, T, n)` and `ext.py::ext_dims(A, T, A_module, n)`,
  `decompose.py::decompose(T, budget)` + `hom.py::is_isomorphic` (count of pairwise
  non-iso indecomposable summands — the char-caveat applies), `morphism.py::direct_sum`/
  `is_direct_summand`, `yoneda.py::baer_extension(M, N, cocycle)` (realizes an Ext¹
  class as `0 -> N -> E -> M -> 0`, returns a `YonedaSequence` with `.middle`),
  `ext.py::ext_dims(..., with_reps=True)` (cocycle representatives),
  `duality.py::dualize` (cotilting = tilting over A^op).
- Produces:
  ```python
  @dataclass
  class TiltingReport:
      is_tilting: bool
      n: int                 # the pd bound tested (classical tilting = 1)
      pd: int | None         # pd T, or None if unresolved within bound
      self_ext_vanishes: bool  # Ext^i(T,T) = 0 for 1 <= i <= n
      num_summands: int      # # pairwise non-iso indecomposable summands of T
      num_vertices: int      # rank K0(A) = |Q_0|
      note: str
      # __bool__ -> is_tilting ; __repr__ names the failing clause
  def is_tilting_module(T, n=1) -> TiltingReport
  def is_cotilting_module(T, n=1) -> TiltingReport      # tilting of D(T) over A^op
  def bongartz_completion(T) -> "Module"
      # for a pd<=1 partial tilting T (Ext^1(T,T)=0): the Bongartz complement's middle
      # term E of the universal extension 0 -> A -> E -> T^d -> 0, d = dim Ext^1(T, A).
      # SELF-CERTIFIED: is_tilting_module(direct_sum(T, E)) must be True.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tilting.py
"""Tilting modules (Plan 44 / C7). Classical n=1 uses the Bongartz COUNT criterion
(#indec summands = #vertices, ASS VI.4); Bongartz completion is self-certified (T (+) E
tilting). Over QQ -- the summand count leans on decompose (char caveat)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.morphism import direct_sum
from quiverlab.modules.tilting import (bongartz_completion, is_cotilting_module,
                                       is_tilting_module)


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@pytest.mark.oracle_selfcert
def test_regular_module_is_always_tilting():
    A = _kA2()
    reg, _, _ = direct_sum(A.projective(1), A.projective(2))
    rep = is_tilting_module(reg)
    assert rep.is_tilting is True and bool(rep) is True
    assert rep.pd == 0 and rep.self_ext_vanishes and rep.num_summands == rep.num_vertices


@pytest.mark.oracle_literature
def test_apr_tilt_of_kA2():
    # T = P1 (+) S1 is a tilting module for kA2 (hereditary => pd S1 = 1, Ext^1(T,T)=0,
    # 2 non-iso summands = 2 vertices). ASS VI worked example (the APR tilt at vertex 1).
    A = _kA2()
    T, _, _ = direct_sum(A.projective(1), A.simple(1))
    assert is_tilting_module(T).is_tilting is True


@pytest.mark.oracle_selfcert
def test_partial_tilting_not_complete_then_bongartz_completes():
    A = _kA2()
    P1 = A.projective(1)
    assert is_tilting_module(P1).is_tilting is False       # only 1 summand, 2 vertices
    E = bongartz_completion(P1)
    T, _, _ = direct_sum(P1, E)
    assert is_tilting_module(T).is_tilting is True          # the certificate IS acceptance


@pytest.mark.oracle_selfcert
def test_cotilting_is_dual():
    A = _kA2()
    reg, _, _ = direct_sum(A.projective(1), A.projective(2))
    # D(A_A) = the injective cogenerator is cotilting.
    assert is_cotilting_module(reg.dualize()).is_tilting is True


@pytest.mark.oracle_selfcert
def test_self_ext_nonvanishing_rejected():
    # A module with Ext^1(T,T) != 0 is not tilting. Over kA2, T = S1 alone has
    # Ext^1(S1,S1)=0 but only 1 summand; use a rad-square-zero example where a simple
    # self-extends. (Concrete builder chosen in implementation; assert is_tilting False
    # and note == the failing clause.)
    ...
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.modules.tilting`

- [ ] **Step 3: Implement `src/quiverlab/modules/tilting.py`**

```python
"""Tilting modules + Bongartz completion (Plan 44 / C7).

A (classical, n=1) tilting module T over A: (1) pd T <= 1; (2) Ext^1(T, T) = 0;
(3) # pairwise non-iso indecomposable summands of T = |Q_0|. Bongartz's theorem (LNM 903,
1981; ASS VI.4): given (1)+(2), (3) is equivalent to the existence of the coresolution
0 -> A -> T^0 -> T^1 -> 0 with T^i in add T, so the count is the checkable form of the
third axiom. General n: the count shortcut is n=1 only; higher n verifies the add(T)-
coresolution of A directly via left add(T)-approximations (Task A). Cotilting = tilting
of D(T) over A^op. Float-free; the summand count inherits the decompose char caveat."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules.decompose import decompose
from quiverlab.modules.ext import ext_dims
from quiverlab.modules.hom import is_isomorphic


@dataclass
class TiltingReport:
    is_tilting: bool
    n: int
    pd: object            # int or None
    self_ext_vanishes: bool
    num_summands: int
    num_vertices: int
    note: str

    def __bool__(self):
        return self.is_tilting

    def __repr__(self):
        if self.is_tilting:
            return f"tilting module (pd {self.pd}, {self.num_summands} summands)"
        return f"not tilting: {self.note}"


def _num_noniso_summands(T, budget=512):
    return len(decompose(T, budget=budget))          # loud on char <= dim


def is_tilting_module(T, n=1, bound=64, budget=512):
    A = T.base_algebra
    verts = list(A.quiver.vertices)
    pd = T.projective_resolution(bound).pd()
    pd_ok = pd is not None and pd <= n
    ext = ext_dims(A, T, T, n)                        # [dim Ext^0, ..., dim Ext^n]
    self_ext = all(d == 0 for d in ext[1:n + 1])
    if n == 1:
        k = _num_noniso_summands(T, budget)
        count_ok = (k == len(verts))
        third = count_ok
        note_third = f"# summands {k} != # vertices {len(verts)}"
    else:
        third = _coresolves_A_in_add_T(A, T, n, budget)   # honest add(T)-coresolution
        k = _num_noniso_summands(T, budget)
        note_third = "A has no add(T)-coresolution of length <= n"
    ok = bool(pd_ok and self_ext and third)
    note = ("ok" if ok else
            (f"pd {pd} > {n}" if not pd_ok else
             ("Ext^i(T,T) != 0 for some 1<=i<=n" if not self_ext else note_third)))
    return TiltingReport(ok, n, pd, self_ext, k, len(verts), note)


def bongartz_completion(T):
    """The Bongartz complement's middle term E of 0 -> A -> E -> T^d -> 0, d = dim
    Ext^1(T, A). Requires pd T <= 1 and Ext^1(T,T)=0 (partial tilting). E is built as the
    Baer extension of T^d by the regular module A_A along an Ext^1(T,A)-basis cocycle."""
    from quiverlab.modules.morphism import direct_sum
    from quiverlab.modules.yoneda import baer_extension
    A = T.base_algebra
    reg, _, _ = direct_sum(*[A.projective(v) for v in A.quiver.vertices])  # A_A
    dims, payload = ext_dims(A, T, reg, 1, with_reps=True)
    d = dims[1]
    if d == 0:
        return reg                                   # already tilting-adjacent: E = A
    Td, _, _ = direct_sum(*([T] * d))
    cocycle = _universal_cocycle(A, T, reg, d, payload)   # P_1(T^d) -> A, block [xi_1..xi_d]
    seq = baer_extension(Td, reg, cocycle)
    return seq.middle


def is_cotilting_module(T, n=1, bound=64, budget=512):
    rep = is_tilting_module(T.dualize(), n=n, bound=bound, budget=budget)
    return rep
```

**Adjust to reality (Task B):**
- **`bongartz_completion`'s universal cocycle is the one construction to arbitrate
  against the oracle.** `baer_extension(M, N, cocycle)` wants an `N.dim x P_1(M).dim`
  matrix that is a genuine `Ext^1(M, N)` cocycle (a module map `P_1(M) -> N` killing the
  image of `d_2`). Build it from `ext_dims(A, T, reg, 1, with_reps=True)`'s
  representatives: the universal element of `Ext^1(T^d, A) ≅ Ext^1(T, A)^d` whose i-th
  component is the i-th basis class `xi_i`. If assembling the `T^d` resolution map is
  fiddly, the fallback is **d successive Baer extensions / iterated pullbacks** (each
  along one `xi_i`), which yields an iso middle term — the SELF-CERTIFICATE
  (`is_tilting_module(T ⊕ E)` True) is the arbiter of which assembly is correct, exactly
  as P37's `end_algebra` lets the oracle pick the composition order. Read
  `yoneda.py::yoneda_sequence`/`baer_extension` and `complex_reps.py`'s cocycle-rep shape
  FIRST; reuse `ext_reps` rather than re-deriving representatives.
- `T.base_algebra` is the P24 side-aware accessor (`module.py:53`); `T.dualize()` returns
  the left/right-flipped dual (P24). For `is_cotilting_module`, `is_tilting_module` on the
  dual runs over `A^op` transparently because every routine reads only `(algebra, action)`
  (P24 invariant) — verify `ext_dims`/`decompose`/`projective_resolution` accept the
  dualized module unchanged.
- `_coresolves_A_in_add_T` (the `n>=2` third axiom): iterate — left add(T)-approximate
  `A`, take the cokernel, approximate again, and after `n` steps assert the final
  cokernel `is_direct_summand` of a power of `T` (P37). Keep it honest: if an
  approximation cannot be certified, raise (never a silent False).
- The `n=1` COUNT is only valid GIVEN pd<=1 and Ext¹(T,T)=0 (Bongartz) — the code tests
  all three, so the count is never applied out of its theorem's hypotheses. State the
  theorem verbatim in the module docstring (already sketched) with the `assem_book`
  citation.
- Fill `test_self_ext_nonvanishing_rejected`'s builder concretely: a radical-square-zero
  algebra with a loop (`RadicalSquareZero` family) has a simple `S` with
  `Ext^1(S,S) != 0`; assert the report's `self_ext_vanishes is False` and `note` names
  the Ext clause.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/tilting.py src/quiverlab/core/algebra.py tests/modules/test_tilting.py
git commit -m "feat(modules): is_tilting/is_cotilting (Bongartz count criterion) + self-certified bongartz_completion"
```

---

### Task C: basic-ization / Gabriel-quiver recovery — THE CRUX

**Files:**
- Create: `src/quiverlab/core/basic.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.primitive_idempotents`,
  `Algebra.basic_algebra`, `Algebra.gabriel_quiver`, `Algebra.presented_form` thin
  wrappers, lazy-import to avoid the cycle, beside `is_basic`, `core/algebra.py:530`)
- Test: `tests/families/test_gabriel_recovery.py`

**Interfaces:**
- Consumes: `Algebra.multiply`/`.T`/`.unit`/`.dim`/`.domain`/`.change_of_basis`/
  `.from_structure_constants` (`core/algebra.py:29,46,92`), `decompose.py::_trace_form_rank`
  and its `_min_poly_coeffs`/`_factor_min_poly`/`_poly_eval_matrix`/`_factoring_supported`
  (the split-element machinery — import the private helpers, or lift them to a shared
  `modules/_polytools.py` in Step 3 if you prefer a public seam), `endomorphism.py`'s
  radical idiom, `linalg_mod::kernel_columns`/`independent_modulo`/`column_space_pivots`/
  `cols_to_matrix`/`mat_rank`, `fields.linalg::solve`, `combinat/quiver.py::Quiver`,
  `families/trivial_extension.py::extract_relations`-style length-lex kernel enumeration
  (`_relation_string`, `_solve_combo` — reuse them; lift to a shared helper if needed).
- Produces:
  ```python
  def primitive_idempotents(A) -> list[list]        # orthogonal primitive idempotents
      # (coordinate vectors) summing to A.unit; char 0 / char > dim only, else loud.
  def idempotent_classes(A) -> list[list[int]]      # partition of the above into iso
      # classes: e_i ~ e_j iff Ae_i ≅ Ae_j, certified by the Schur corner-dim test
      # dim_k(ebar_i S ebar_j) >= 1 in the semisimple quotient S = A/rad.
  def basic_algebra(A) -> Algebra                   # e A e for one idempotent per class
  def gabriel_quiver(A) -> Quiver                   # vertices = classes; arrows i->j =
      # dim_k f_i (rad B / rad^2 B) f_j in the basic algebra B.
  def presented_form(A) -> Algebra                  # kQ/I, Q = gabriel_quiver(A), certified
      # per instance by dim(kQ/I) == dim basic_algebra(A) + a random-element homomorphism
      # spot-check. Loud refusal off char-scope or on a non-split division-algebra block.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_gabriel_recovery.py
"""Basic-ization + Gabriel-quiver recovery of a structure-constant algebra (Plan 44 / C7).
The Wedderburn crux: primitive idempotents via the trace-form radical + a semisimple
quotient + Newton lifting; the presentation via a length-lex kernel enumeration. Over QQ /
large-prime GF only (the trace-form radical needs char 0 or char > dim)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.core.algebra import Algebra
from quiverlab.core.basic import (basic_algebra, gabriel_quiver,
                                   idempotent_classes, presented_form,
                                   primitive_idempotents)
from quiverlab.fields import QQ


def _M2(field):
    """M_2(k) as structure constants in the matrix-unit basis (E11,E12,E21,E22):
    E_ij E_kl = delta_jk E_il ; unit = E11 + E22."""
    z = "0"
    def e(i):  # basis vector
        v = [z, z, z, z]; v[i] = "1"; return v
    # index: 0=E11,1=E12,2=E21,3=E22 ; pair (row,col): 0=(1,1),1=(1,2),2=(2,1),3=(2,2)
    rc = {0: (1, 1), 1: (1, 2), 2: (2, 1), 3: (2, 2)}
    idx = {v: k for k, v in rc.items()}
    T = [[["0"] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        i, j = rc[a]
        for b in range(4):
            k, l = rc[b]
            if j == k:
                T[a][b][idx[(i, l)]] = "1"
    return Algebra.from_structure_constants(T, ["1", "0", "0", "1"], field=field)


@pytest.mark.oracle_literature
def test_matrix_algebra_basic_is_the_field():
    A = _M2(GF(32003))                       # char 32003 > dim 4
    idems = primitive_idempotents(A)
    assert len(idems) == 2                    # E11, E22 (both rank-1)
    classes = idempotent_classes(A)
    assert len(classes) == 1                  # conjugate: one iso class
    B = basic_algebra(A)
    assert B.dim == 1                          # e M_2 e = k
    Q = gabriel_quiver(A)
    assert len(list(Q.vertices)) == 1 and len(Q.arrows) == 0
    assert presented_form(A).dim == 1


@pytest.mark.oracle_selfcert
def test_primitive_idempotents_are_a_complete_orthogonal_set():
    A = _M2(QQ)
    idems = primitive_idempotents(A)
    dom = A.domain
    # orthogonal: e_i e_j = delta_ij e_i ; complete: sum e_i = unit.
    total = [dom.zero()] * A.dim
    for i, ei in enumerate(idems):
        assert A.multiply(ei, ei) == ei                       # idempotent
        for j, ej in enumerate(idems):
            prod = A.multiply(ei, ej)
            if i != j:
                assert all(dom.is_zero(x) for x in prod)      # orthogonal
        total = [dom.add(total[t], ei[t]) for t in range(A.dim)]
    assert total == list(A.unit)                              # complete


@pytest.mark.oracle_selfcert
def test_kA2_round_trip_identity():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)   # basic already
    B = presented_form(A)
    assert B.dim == A.dim == 3
    assert B.cartan_matrix() == A.cartan_matrix()             # recovered kA2 == kA2


@pytest.mark.oracle_crossengine
def test_end_of_regular_recovers_kA3():
    # End_A(A_A) ~ A (P37 endomorphism.py). Take A = kA3 (hereditary 1->2->3), build the
    # regular module's End as a presentation-LESS structure-constant algebra, and recover
    # kA3 from it -- tying to P37's regular_corner_dims (Cartan of End = Cartan of A).
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=QQ)
    from quiverlab.modules.endomorphism import end_algebra, regular_corner_dims
    from quiverlab.modules.morphism import direct_sum
    reg, _, _ = direct_sum(*[A.projective(v) for v in (1, 2, 3)])
    E = end_algebra(reg)                                       # quiver is None (presentation-less)
    assert E.quiver is None
    R = presented_form(E)
    assert R.dim == A.dim                                      # dim kA3 = 6
    assert R.cartan_matrix() == A.cartan_matrix()             # recovered kA3
    assert regular_corner_dims(A) == [[int(x) for x in row]   # P37 sided oracle tie
                                      for row in A.cartan_matrix()]


@pytest.mark.oracle_selfcert
def test_char_caveat_refuses_loudly():
    from quiverlab.errors import QuiverlabError
    A = _M2(GF(2))                            # char 2 <= dim 4: trace-form radical unreliable
    with pytest.raises(QuiverlabError):
        primitive_idempotents(A)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.core.basic`

- [ ] **Step 3a: Implement the radical + semisimple quotient** (`core/basic.py`)

```python
"""Basic-ization + Gabriel-quiver recovery of a structure-constant Algebra (Plan 44 / C7).

Exact Wedderburn over char 0 or char > dim A ONLY (the trace-form radical -- the sole
exact radical route in the library -- is rigorous exactly there; Dickson / Cohen-Ivanyos-
Wales, the same bound decompose.py uses). Pipeline:
  1. rad A = nullspace of the REGULAR-representation trace form G[i][j] = tr(L_i L_j).
  2. S = A / rad, a semisimple structure-constant algebra (rad is a two-sided nilpotent
     ideal).
  3. a complete set of primitive orthogonal idempotents of S, by splitting-element
     refinement (decompose's minimal-polynomial + coprime-factor machinery, applied to
     the regular rep of each corner). A non-split division-algebra block refuses loudly.
  4. lift each idempotent to A, orthogonally, by the Newton iteration e <- 3e^2 - 2e^3
     (terminates: rad nilpotent, the defect e^2 - e is squared each step).
  5. group by the Schur corner-dim test dim_k(ebar_i S ebar_j) >= 1 (== Ae_i ≅ Ae_j).
Then basic_algebra = eAe (one idempotent per class); gabriel_quiver reads arrows off
rad/rad^2 of the basic algebra; presented_form extracts kQ/I by a length-lex kernel
enumeration (the TrivialExtension idiom) and certifies dim(kQ/I) == dim(basic) plus a
random-element homomorphism spot-check. Float-free; loud off char-scope."""
from __future__ import annotations

from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm


def _left_mult(A, i):
    """The matrix of left multiplication by basis element b_i: column k is b_i * b_k."""
    n = A.dim
    return [[A.T[i][k][t] for k in range(n)] for t in range(n)]


def _radical_coords(A):
    """rad A as coordinate column-vectors: nullspace of the regular-rep trace form
    G[i][j] = tr(L_i L_j). Rigorous iff char 0 or char > dim A -- else loud."""
    dom = A.domain
    n = A.dim
    char = dom.characteristic
    if not (char == 0 or char > n):
        raise QuiverlabError(
            f"basic-ization: the trace-form radical is unreliable in characteristic "
            f"{char} <= dim A = {n}", hint="work over QQ or a prime > dim A")
    L = [_left_mult(A, i) for i in range(n)]
    G = lm.zeros(n, n, dom)
    for i in range(n):
        for j in range(n):
            prod = lm.matmul(L[i], L[j], dom)
            s = dom.zero()
            for t in range(n):
                s = dom.add(s, prod[t][t])
            G[i][j] = s
    return lm.kernel_columns(G, dom)          # [] when A is semisimple


def _semisimple_quotient(A, rad_cols):
    """S = A / rad as a structure-constant Algebra, plus the projection A -> S (a
    dim(S) x dim(A) matrix) and the section S -> A (complement basis columns)."""
    dom, n = A.domain, A.dim
    ident = lm.identity(n, dom)
    std = [lm.col(ident, j) for j in range(n)]
    comp_idx = lm.independent_modulo(std, rad_cols, dom)       # complement of rad
    comp = [std[j] for j in comp_idx]
    P = lm.cols_to_matrix(comp + [list(c) for c in rad_cols])  # [comp | rad] basis of A
    m = len(comp)
    def proj(vec):                                             # A-vector -> S-coords
        x = linalg.solve(P, list(vec), dom)
        return x[:m]
    T = []
    for a in range(m):
        row = []
        for b in range(m):
            prod = A.multiply(comp[a], comp[b])
            row.append([proj(prod)[c] for c in range(m)])      # mod rad
        T.append([[T_ab for T_ab in cell] for cell in [ [row[b][c] for c in range(m)] for b in range(m)]][:])  # see note
        row = [[proj(A.multiply(comp[a], comp[b]))[c] for c in range(m)] for b in range(m)]
        T[-1] = row
    unit_S = proj(list(A.unit))
    S = Algebra.from_structure_constants(
        [[[ _fmt(x, dom) for x in cell] for cell in rowa] for rowa in T],
        [_fmt(x, dom) for x in unit_S], field=dom, check=True)  # check=True self-certifies S
    return S, comp, proj
```

**Adjust to reality (Step 3a):** the doubled `T` assignment above is deliberately
sketchy — write it as a single clean nested comprehension
`T = [[[proj(A.multiply(comp[a], comp[b]))[c] for c in range(m)] for b in range(m)] for a in range(m)]`;
the point is the semisimple structure constants are the A-products of complement reps,
projected mod rad. `_fmt(x, dom)` renders a Domain element back into the token form
`from_structure_constants` re-parses (mirror `trivial_extension._coeff_sign_mag` / the
`from_structure_constants` callers — over QQ pass `Fraction`, over GF(p) pass the int;
or, simpler, extend `from_structure_constants` to accept already-coerced Domain elements
and pass `check=True` — read `core/algebra.py:46-60` and pick the path an existing caller
uses). `check=True` validating associativity IS the self-certificate that `S` is a
well-defined quotient algebra.

- [ ] **Step 3b: Implement primitive idempotents (semisimple split + Newton lift)**

```python
def _semisimple_primitive_idempotents(S):
    """A complete set of primitive orthogonal idempotents of the SEMISIMPLE algebra S
    (S-coordinate vectors), by refining 1_S via splitting elements. A corner that neither
    splits nor certifies as a division algebra (dim-1 corner) refuses loudly."""
    from quiverlab.modules.decompose import (_factor_min_poly, _factoring_supported,
                                             _min_poly_coeffs, _poly_eval_matrix)
    dom = S.domain
    if not _factoring_supported(dom):
        raise QuiverlabError(
            "basic-ization: cannot factor over this domain within the certified budget",
            hint="supported: GF(p), QQ, CC number fields (as in decompose)")
    out = []
    stack = [list(S.unit)]                     # corners to refine, each an idempotent e
    while stack:
        e = stack.pop()
        split = _split_idempotent(S, e, _min_poly_coeffs, _factor_min_poly,
                                  _poly_eval_matrix)
        if split is None:                       # e primitive (corner is a division algebra)
            out.append(e)
        else:
            e1, e2 = split
            stack.extend([e1, e2])
    return out


def _split_idempotent(S, e, minpoly, factor, evalmat):
    """Find an s in the corner eSe whose regular-rep minimal polynomial has >= 2 coprime
    factors, and return two orthogonal idempotents e1 + e2 = e (via the CRT idempotent);
    or None when no splitting element exists over the corner-basis + small combinations
    (=> e primitive). The corner unit is e."""
    ...  # regular-rep of eSe (left-mult on the corner subspace); reuse decompose's
        # _candidate_endomorphisms search shape; the CRT idempotent is (v*g)(s) evaluated
        # in S with e as the unit. Terminates: each split lowers corner dimension.


def _newton_lift(A, x):
    """Lift a near-idempotent x (x^2 == x mod rad) to an exact idempotent of A via
    e <- 3e^2 - 2e^3, iterated until A.multiply(e,e) == e. Terminates (rad nilpotent)."""
    dom = A.domain
    e = list(x)
    for _ in range(A.dim + 1):                  # rad nilpotency index <= dim
        if A.multiply(e, e) == e:
            return e
        e2 = A.multiply(e, e)
        e3 = A.multiply(e2, e)
        e = [dom.add(dom.mul(dom.coerce(3), e2[t]),
                     dom.neg(dom.mul(dom.coerce(2), e3[t]))) for t in range(A.dim)]
    raise QuiverlabError("basic-ization: Newton idempotent lift did not converge",
                         hint="the radical was not nilpotent -- report this algebra")


def primitive_idempotents(A):
    rad = _radical_coords(A)                     # loud off char-scope
    S, comp, proj = _semisimple_quotient(A, rad)
    prim_S = _semisimple_primitive_idempotents(S)
    # orthogonal lift of the complete set: refine within g = 1 - sum(collected).
    dom = A.domain
    embed = lambda sc: _combine(comp, sc, dom)   # S-coords -> A-vector (complement rep)
    collected = []
    g = list(A.unit)
    for sc in prim_S:
        x = embed(sc)
        y = A.multiply(g, A.multiply(x, g))      # corner preimage g x g
        e = _newton_lift(A, y)
        collected.append(e)
        g = [dom.add(g[t], dom.neg(e[t])) for t in range(A.dim)]   # g <- g - e
    _assert_complete_orthogonal(A, collected)    # sum == unit, pairwise orthogonal (loud)
    return collected
```

**Adjust to reality (Step 3b):** flesh out `_split_idempotent` by mirroring
`decompose._try_split`/`_candidate_endomorphisms` (the same search) but on the corner's
regular representation (left-multiplication on the subspace `eSe`), with the corner unit
`e` playing the role of `I`; the CRT idempotent is `(v·g)(s)` evaluated in `S`. `_combine`
is the obvious `sum sc[a] * comp[a]`. The `M₂(k)→k` and `kA3` recovery tests are the
ARBITER of the orthogonal-lift order and the CRT sign, exactly as P37's `end_algebra`
lets the regular-module oracle pick the composition order — if the completeness assertion
fails, the lift order is wrong; adjust and re-run. Over a split semisimple `S` (all blocks
`M_n(k)`) the splitting element always exists once `|k| > dim` (guaranteed by
`char > dim`); a genuine non-split division-algebra block (e.g. a quaternion corner over
QQ) is out of scope and `_split_idempotent`'s "no split found, corner dim > 1" branch
raises loudly — the honest-scope refusal.

- [ ] **Step 3c: Implement classes, basic_algebra, gabriel_quiver, presented_form**

```python
def idempotent_classes(A):
    """Partition primitive_idempotents(A) by Ae_i ≅ Ae_j, certified by the Schur corner-
    dim test in the semisimple quotient: dim_k(ebar_i S ebar_j) >= 1 <=> same Wedderburn
    block <=> the projectives (hence the primitives) are isomorphic."""
    rad = _radical_coords(A)
    S, comp, proj = _semisimple_quotient(A, rad)
    idems = primitive_idempotents(A)
    bars = [proj_to_S(A, e, proj) for e in idems]          # ebar_i in S
    n = len(idems)
    groups = []
    assigned = [False] * n
    for i in range(n):
        if assigned[i]:
            continue
        cls = [i]; assigned[i] = True
        for j in range(i + 1, n):
            if not assigned[j] and _corner_nonzero(S, bars[i], bars[j]):
                cls.append(j); assigned[j] = True
        groups.append(cls)
    return groups


def basic_algebra(A):
    idems = primitive_idempotents(A)
    classes = idempotent_classes(A)
    dom = A.domain
    e = [dom.zero()] * A.dim
    for cls in classes:                                    # one idempotent per class
        rep = idems[cls[0]]
        e = [dom.add(e[t], rep[t]) for t in range(A.dim)]
    return _corner_algebra(A, e)                           # e A e as structure constants


def gabriel_quiver(A):
    from quiverlab.combinat.quiver import Quiver
    B = basic_algebra(A)
    fi = _basic_primitive_idempotents(B)                   # one per vertex, in B
    radB = _radical_coords(B)
    rad2 = _rad_square(B, radB)                             # span of products
    verts = list(range(1, len(fi) + 1))
    arrows = {}
    a = 0
    for i, ei in enumerate(fi):
        for j, ej in enumerate(fi):
            r = _corner_rad_over_rad2_dim(B, ei, ej, radB, rad2)   # dim f_i(rad/rad^2)f_j
            for _ in range(r):
                a += 1
                arrows[f"a{a}"] = (verts[i], verts[j])     # direction: oracle-arbitrated
    return Quiver(verts, arrows)


def presented_form(A):
    from quiverlab.combinat.quiver import Quiver
    B = basic_algebra(A)
    Q = gabriel_quiver(A)
    # pi: kQ -> B on generators (vertex -> f_i ; arrow -> a rad rep lifting a rad/rad^2
    # basis element), then length-lex kernel enumeration (the TrivialExtension idiom),
    # certified by dim(kQ/I) == dim B and a random-element homomorphism spot-check.
    ...
```

**Adjust to reality (Step 3c):**
- `_corner_algebra(A, e)`: basis of the subspace `eAe` = a maximal independent subset of
  `{ e * b_k * e : k }` (via `lm.column_space_pivots`); structure constants = the
  A-products of those basis vectors expressed back in the chosen basis (via
  `fields.linalg.solve`); unit = `e` in that basis. Build via
  `Algebra.from_structure_constants(check=True)` (associativity + unit self-certify).
- `gabriel_quiver`'s **arrow direction** (`i -> j` vs `j -> i`) and whether the count is
  `dim f_i(rad/rad²)f_j` or its transpose is a CONVENTION pinned by
  `test_kA2_round_trip_identity` + `test_end_of_regular_recovers_kA3`: the Cartan-matrix
  equality fixes it (path composition is left-to-right, ASS). Write both, let the oracle
  decide, record the outcome in the docstring — the P37/P40 precedent.
- `presented_form`'s relation extraction is the `families/trivial_extension.py::
  _presented_trivial_extension` machinery reused wholesale: build `pi` on generators into
  `B` (`Tsc`'s role), run the **same** length-lex `extract_relations(max_len)` with the
  same per-corner reduced-echelon reduction and `_relation_string` emission, certify
  `dim(kQ/I) == dim B` with the widen-once window (`loewy_length + 2`), and ADDITIONALLY
  spot-check on random element pairs that the iso `kQ/I -> B` is multiplicative (the
  Gabriel recovery has no bimodule-socle shortcut, so the multiplicativity spot-check is
  the extra certificate beyond dimension). Lift `_relation_string`/`_solve_combo`/
  `extract_relations` into a shared helper module if importing them from
  `families/trivial_extension.py` reads badly (they are presentation-agnostic).
- `is_basic` stays UNCHANGED (`recognizers.is_basic` returns True for every `kQ/I`, which
  is correct — a presented algebra IS basic); the new surface handles the
  presentation-LESS structure-constant / non-basic case. Note this in the docstring so
  the two are never conflated.

- [ ] **Step 4: Run tests, verify pass** (this is the crux — expect several iterations of
  the oracle-arbitrated conventions above)

Run: `... -m pytest tests/families/test_gabriel_recovery.py -v` — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/core/basic.py src/quiverlab/core/algebra.py tests/families/test_gabriel_recovery.py
git commit -m "feat(core): primitive_idempotents/basic_algebra/gabriel_quiver/presented_form -- exact Gabriel recovery of a structure-constant algebra (char-caveat-honest)"
```

---

### Task D: one-point extension `A[M]`

**Files:**
- Create: `src/quiverlab/families/one_point.py`
- Modify: `src/quiverlab/families/__init__.py` (export `OnePointExtension`),
  `src/quiverlab/families/discover.py` (CATALOG `FamilyInfo`)
- Test: `tests/families/test_one_point.py`

**Interfaces:**
- Consumes: `Module.top()` (top generators of `M`, one new arrow per basis vector of
  `top(M)_v`), `Module.projective_resolution().pd()`, `combinat/quiver.py::Quiver` +
  `Quiver.algebra`, the `families/trivial_extension.py` length-lex kernel idiom
  (`_relation_string`, `extract_relations`, `_solve_combo`), `invariants/cartan.py`
  (Cartan-block oracle), `groebner`'s `NotFiniteDimensionalError`.
- Produces:
  ```python
  def OnePointExtension(A, M) -> Algebra
      # A[M] = [[k, M], [0, A]] presented: a new SOURCE vertex w, one arrow w -> v per
      # basis vector of top(M)_v, relations from M's structure (M acts on the new arrows),
      # extracted by the length-lex kernel enumeration and certified per instance by
      # dim A[M] == 1 + dim M + dim A.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_one_point.py
"""One-point extension A[M] = [[k,M],[0,A]] (Plan 44 / C7). Certified per instance by
dim A[M] == 1 + dim M + dim A; the Cartan block form and the pd(S_w) = pd_A(M) + 1
identity are the literature oracles (ASS III / Happel)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.one_point import OnePointExtension

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_dim_certificate_and_quiver():
    # A[S1] for A = kA2: new source w, arrow w->1 (top S1 at vertex 1), and the relation
    # (w->1->2) = 0 because a acts as 0 on S1. Q' = w->1->2 (A3 shape) with one relation.
    A = _kA2()
    S1 = A.simple(1)
    B = OnePointExtension(A, S1)
    assert B.dim == 1 + S1.dim + A.dim == 5          # NOT kA3 (dim 6): the relation cuts 1
    assert len(list(B.quiver.vertices)) == 3
    assert len(B.relations) >= 1                      # the w->1->2 = 0 relation


def test_cartan_block_form():
    # C_{A[M]} = [[1, dimvec M], [0, C_A]] (up to the pinned vertex order; the Cartan
    # equality arbitrates the source-vertex placement + row/col convention).
    A = _kA2()
    S1 = A.simple(1)
    B = OnePointExtension(A, S1)
    CB, CA = B.cartan_matrix(), A.cartan_matrix()
    # w is the extra vertex; strip its row/col and compare the A-block to C_A.
    # (exact index bookkeeping written in the implementation; assert the A-block matches)
    assert _strip_new_vertex(CB) == CA               # helper in the test module


def test_pd_of_new_simple_is_pd_M_plus_one():
    # pd_{A[M]}(S_w) = pd_A(M) + 1. For M = S1 over kA2, pd_A(S1) = 1, so pd = 2.
    A = _kA2()
    S1 = A.simple(1)
    B = OnePointExtension(A, S1)
    w = _new_vertex(B, A)                              # the source vertex of A[M]
    pd_new = B.simple(w).projective_resolution(8).pd()
    pd_M = S1.projective_resolution(8).pd()
    assert pd_new == pd_M + 1 == 2


@pytest.mark.oracle_selfcert
def test_projective_module_extension_has_no_relation():
    # M = S2 (= P2, projective, pd 0): A[S2] has quiver w->2, 1->2 and NO length-2
    # relation; dim 5, pd(S_w) = 1.
    A = _kA2()
    S2 = A.simple(2)
    B = OnePointExtension(A, S2)
    assert B.dim == 5
    w = _new_vertex(B, A)
    assert B.simple(w).projective_resolution(8).pd() == 1
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.families.one_point`

- [ ] **Step 3: Implement `src/quiverlab/families/one_point.py`**

Build the triangular structure-constant algebra `T = A[M]` (block product
`(λ, m, a)(λ', m', a') = (λλ', λ m' + m·a', a a')`, `m·a'` the right A-action) as the `pi`
TARGET, then present it exactly like `_presented_trivial_extension`:
- `Q' = Q_A` + new source vertex `w`; arrows `Q_A.arrows` + `{ w -> v : one per basis of
  top(M)_v }` (compute `top(M)` = `M.top()`; its `v`-component dimension is the number of
  `w -> v` arrows).
- `pi` on generators: original arrows -> their `A`-block image in `T`; each `w -> v` arrow
  -> the corresponding `top(M)_v` generator lifted into the `M`-block of `T`.
- `extract_relations(max_len)`: the **same** per-corner length-lex reduced-echelon kernel
  enumeration as `trivial_extension.py`, emitting `_relation_string`s.
- Certify `dim(kQ'/I') == 1 + dim M + dim A` with the widen-once window; loud
  `QuiverlabError` on failure; `NotFiniteDimensionalError` propagates if the presented
  backbone cannot certify finiteness (cannot happen for a genuine one-point extension of
  a finite-dim `A` by a finite-dim `M`, but keep the honest propagation).
- `B._family_citations = ("assem_book", "happel_question")` (ASS III one-point
  extensions; Happel's change-of-rings reads off the LES — the P42 oracle).

**Adjust to reality (Task D):** the source-vertex convention (`w` a source, arrows OUT of
`w`) and the right-vs-left reading of `M` are pinned by `test_cartan_block_form` +
`test_pd_of_new_simple_is_pd_M_plus_one`: derive the block Cartan
`[[1, dim-vector M],[0, C_A]]` in the test comment and let the Cartan equality arbitrate
the exact index placement (P37/P40 "oracle decides" precedent). Reuse the
`trivial_extension.py` helpers by import or lift them to a shared
`families/_present.py`; do NOT re-derive the mini-Gröbner.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/families/one_point.py src/quiverlab/families/__init__.py \
        src/quiverlab/families/discover.py tests/families/test_one_point.py
git commit -m "feat(families): OnePointExtension A[M] -- presented triangular algebra, Cartan-block + pd(S_w)=pd(M)+1 certified"
```

---

### Task E: finite repetitive-algebra slices

**Files:**
- Create: `src/quiverlab/families/repetitive.py`
- Modify: `src/quiverlab/families/__init__.py` (export `repetitive_slice`),
  `src/quiverlab/families/discover.py` (CATALOG `FamilyInfo`)
- Test: `tests/families/test_repetitive.py`

**Interfaces:**
- Consumes: `families/trivial_extension.py::_bimodule_socle`/`_socle_arrows`/
  `_corner_maps`/`extract_relations` (the `D(A)` connecting-bimodule machinery — the
  repetitive algebra's connecting arrows ARE the trivial-extension dual arrows, but wired
  between CONSECUTIVE copies instead of wrapping around), `invariants/pathbasis::
  path_type_basis`, `Quiver.algebra`, `invariants/cartan.py`.
- Produces:
  ```python
  def repetitive_slice(A, copies) -> Algebra
      # the finite truncation of the repetitive algebra hat(A): `copies` copies of A
      # (copies >= 1) joined by `copies - 1` copies of the D(A) connecting bimodule.
      # Presented via the TrivialExtension idiom; certified per instance by
      # dim == (2*copies - 1) * dim A. HONEST SCOPE: hat(A) is infinite-dimensional; this
      # ships certified finite slices only.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_repetitive.py
"""Finite repetitive-algebra slices (Plan 44 / C7). copies=1 is A itself (byte-Cartan);
the general slice is certified by dim == (2*copies - 1)*dim A. The full repetitive algebra
is infinite-dimensional -- out of scope, stated on the verification page."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.repetitive import repetitive_slice

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_one_copy_is_A():
    A = _kA2()
    R1 = repetitive_slice(A, 1)
    assert R1.dim == A.dim
    assert R1.cartan_matrix() == A.cartan_matrix()          # copies=1 == A


def test_two_slice_dimension():
    A = _kA2()                                              # dim 3
    R2 = repetitive_slice(A, 2)
    assert R2.dim == (2 * 2 - 1) * A.dim == 9               # 2 A's + 1 D(A)


@pytest.mark.parametrize("copies", [1, 2, 3])
def test_slice_dimension_certificate(copies):
    A = _kA2()
    R = repetitive_slice(A, copies)
    assert R.dim == (2 * copies - 1) * A.dim


def test_bad_copies_refused():
    A = _kA2()
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError):
        repetitive_slice(A, 0)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.families.repetitive`

- [ ] **Step 3: Implement `src/quiverlab/families/repetitive.py`**

- `copies == 1` returns `A` unchanged (byte-Cartan identity; the base case).
- For `copies >= 2`: `Q_slice` = `copies` disjoint relabeled copies of `Q_A` (vertices
  `(v, c)` for `c in 0..copies-1`), plus, for each `c in 0..copies-2`, one connecting
  arrow per bimodule-socle vector of `A` running from copy `c+1` to copy `c` (the
  trivial-extension dual arrows `w in e_i A e_j` -> `te : (j, c+1) -> (i, c)`; the
  reversed direction as in `trivial_extension.py:227`, but connecting consecutive copies,
  NOT wrapping copy `copies-1` back to copy `0`).
- `pi` target = the block-matrix structure-constant algebra `hat(A)_slice` (A on the
  diagonal blocks, `D(A)` on the super/sub-diagonal connecting blocks per the repetitive
  multiplication); `extract_relations` = the same length-lex enumeration; certify
  `dim == (2*copies - 1) * dim A` with the widen-once window; loud `QuiverlabError` on
  failure.
- Refuse `copies < 1` loudly (`QuiverlabError`, `hint="copies >= 1; the full repetitive
  algebra is infinite-dimensional"`). Fall back to a structure-constant build for
  presentation-less / non-string-representable bases (mirror `TrivialExtension`'s D3
  fallback).
- `R._family_citations = ("happel_trivial_extension", "hughes_waschbusche", "assem_book")`
  (Happel's repetitive-algebra chapter; Hughes-Waschbüsch introduced the repetitive
  algebra).

**Adjust to reality (Task E):** the connecting-bimodule wiring is the ONE structural
choice; the `dim == (2*copies-1)*dim A` certificate is the arbiter (a mis-wired
connection over/under-counts or goes infinite). Do NOT ship a T(A) ≅ hat(A)/(nu) orbit
oracle — the Hughes-Waschbüsch orbit identification is a folklore-adjacent claim not
re-derived here; note it as the successor (a slice with periodic identification) in the
docstring + verification page, no unverified numbers. Read `trivial_extension.py`
end-to-end before writing this; the reuse is heavy.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/families/repetitive.py src/quiverlab/families/__init__.py \
        src/quiverlab/families/discover.py tests/families/test_repetitive.py
git commit -m "feat(families): repetitive_slice -- certified finite truncations of the repetitive algebra (honest finite-slice scope)"
```

---

### Task F: Jacobian algebra from a quiver-with-potential `(Q, W)`

**Files:**
- Create: `src/quiverlab/families/jacobian.py`
- Modify: `src/quiverlab/families/__init__.py` (export `Potential`, `JacobianAlgebra`),
  `src/quiverlab/families/discover.py` (CATALOG `FamilyInfo`)
- Test: `tests/families/test_jacobian.py`

**Interfaces:**
- Consumes: `combinat/quiver.py::Quiver` (`word_source`/`word_target`/`source`/`target`
  for cyclic-word validation), `Quiver.algebra(relations, field, degree_bound)`,
  `combinat/relations.py` grammar (the cyclic derivatives are emitted as relation
  strings), `groebner`'s `NotFiniteDimensionalError` (the honest infinite-dim refusal),
  the `preprojective.py::_AUTO_DEGREE_BOUND` idiom (a per-family degree-bound heuristic),
  `families/preprojective.py` (the preprojective-as-Jacobian cross-engine oracle).
- Produces:
  ```python
  class Potential:
      def __init__(self, quiver, terms)   # terms = [(coeff, word)], each word a CYCLE;
                                          # loud on a non-cyclic word (validated).
  def cyclic_derivative(W, arrow) -> str  # partial_a W, as a relation-grammar string
  def JacobianAlgebra(Q, W, field=None, degree_bound=None) -> Algebra
      # kQ / (partial_a W : a in Q_1) via Quiver.algebra; NotFiniteDimensionalError is the
      # honest Jacobian-INFINITE refusal. Auto degree bound per the preprojective idiom.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_jacobian.py
"""Jacobian algebra Jac(Q, W) = kQ/(d_a W) (Plan 44 / C7). The 3-cycle triangle is the
hand-derived pin; the preprojective-as-Jacobian identity is the cross-engine oracle; an
under-constrained potential is the honest NotFiniteDimensionalError refusal.
Derksen-Weyman-Zelevinsky 2008; Labardini-Fragoso 2009."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import NotFiniteDimensionalError
from quiverlab.families.jacobian import JacobianAlgebra, Potential, cyclic_derivative
from quiverlab.fields import QQ


def _triangle():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)})


@pytest.mark.oracle_selfcert
def test_potential_rejects_non_cyclic_word():
    Q = _triangle()
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError):
        Potential(Q, [(1, ("a", "b"))])          # a*b : 1 -> 3, NOT a cycle


@pytest.mark.oracle_selfcert
def test_cyclic_derivatives_of_abc():
    Q = _triangle()
    W = Potential(Q, [(1, ("a", "b", "c"))])     # the 3-cycle a*b*c
    assert cyclic_derivative(W, "a") == "b*c"    # d_a(abc) = bc
    assert cyclic_derivative(W, "b") == "c*a"    # d_b(abc) = ca
    assert cyclic_derivative(W, "c") == "a*b"    # d_c(abc) = ab


@pytest.mark.oracle_literature
def test_triangle_jacobian_dimension():
    # Jac = kQ/(bc, ca, ab): all length-2 paths die => rad^2 = 0 => dim = 3 vertices +
    # 3 arrows = 6. Hand-derived.
    Q = _triangle()
    W = Potential(Q, [(1, ("a", "b", "c"))])
    J = JacobianAlgebra(Q, W, field=QQ)
    assert J.dim == 6


@pytest.mark.oracle_crossengine
def test_preprojective_is_the_jacobian_of_the_double():
    # Pi(Q) = Jac(double Q, sum_a (a a* - a* a)) (the preprojective potential). Cross-check
    # against the PreprojectiveAlgebra family: same dim + same Cartan matrix.
    from quiverlab import PreprojectiveAlgebra
    Pi = PreprojectiveAlgebra("A3", field=QQ)
    Q2, W = _double_with_preprojective_potential("A3")     # helper in the test module
    J = JacobianAlgebra(Q2, W, field=QQ)
    assert J.dim == Pi.dim
    assert J.cartan_matrix() == Pi.cartan_matrix()


@pytest.mark.oracle_selfcert
def test_infinite_jacobian_refused_loudly():
    # one vertex, two loops, W = 0: relations empty => free algebra k<x,y> => infinite.
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    W = Potential(Q, [])                          # empty potential
    with pytest.raises(NotFiniteDimensionalError):
        JacobianAlgebra(Q, W, field=QQ)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.families.jacobian`

- [ ] **Step 3: Implement `src/quiverlab/families/jacobian.py`**

```python
"""Jacobian algebra of a quiver with potential (Q, W) (Plan 44 / C7).

A potential W is a k-linear combination of cyclic words (cycles in Q). The cyclic
derivative partial_a W = sum over cyclic rotations placing `a` first: for W = coeff * (a p)
(a cyclic word starting a), partial_a W picks up coeff * p (delete the leading `a`, keep
the rest of the cycle). The Jacobian (DWZ) algebra is Jac(Q, W) = kQ / (partial_a W : a in
Q_1). It may be infinite-dimensional -- NotFiniteDimensionalError is the honest refusal.
Derksen-Weyman-Zelevinsky, Selecta Math 2008; surface case Labardini-Fragoso, Proc. LMS
2009. Float-free; coefficients are exact (Fraction / prime-field / sanctioned scalars)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.families.trivial_extension import _relation_string   # reuse the emitter


class Potential:
    def __init__(self, quiver, terms):
        self.quiver = quiver
        self.terms = []
        for coeff, word in terms:
            if not word or not quiver.compose_ok(word) \
                    or quiver.word_source(word) != quiver.word_target(word):
                raise QuiverlabError(
                    f"Potential: {word!r} is not a cyclic word in Q",
                    hint="each potential term must be an oriented cycle "
                         "(source of the first arrow == target of the last)")
            self.terms.append((coeff, tuple(word)))


def cyclic_derivative(W, arrow):
    """partial_arrow W as a relation-grammar string (parallel sum of the deletions)."""
    out_terms = []
    for coeff, word in W.terms:
        for r in range(len(word)):
            if word[r] == arrow:
                rot = word[r + 1:] + word[:r]           # delete `arrow`, cyclically rotate
                out_terms.append((coeff, rot))
    if not out_terms:
        return ""                                       # d_a W = 0 (no relation)
    return _relation_string(out_terms)                  # combinat grammar string


def JacobianAlgebra(Q, W, field=None, degree_bound=None):
    if not isinstance(W, Potential):
        raise QuiverlabError("JacobianAlgebra needs a Potential, got "
                             f"{type(W).__name__}", hint="wrap W in Potential(Q, terms)")
    rels = [s for s in (cyclic_derivative(W, a) for a in Q.arrows) if s]
    if degree_bound is None:
        degree_bound = _auto_degree_bound(Q)            # per-family heuristic; may be None
    A = Q.algebra(relations=rels, field=field, degree_bound=degree_bound)
    A._family_citations = ("derksen_weyman_zelevinsky", "labardini", "assem_book")
    return A
```

**Adjust to reality (Task F):**
- The cyclic derivative's parallel sum must land in a single corner (all deletions of one
  cyclic word `w` share source `target(a)` and target `source(a)`) — assert it, so a
  malformed `W` fails loudly rather than producing a non-parallel relation the grammar
  rejects. Verify `_relation_string`'s sign handling round-trips the coefficient sign
  (`trivial_extension.py:88-114`); the preprojective potential `sum (a a* - a* a)` has a
  genuine minus, so the emitter's `-` path is exercised.
- `_auto_degree_bound(Q)`: mirror `preprojective.py:45-58` — a small heuristic table /
  function keyed on `|Q_1|` (or a fixed generous default like `2*|Q_0| + 4`), returning
  `None` to leave the Groebner route's own bound in force. The route raises
  `NotFiniteDimensionalError` (naming an arrow cycle) when the ideal does not cut the
  algebra finite — which is the correct Jacobian-infinite verdict. Do NOT fabricate a
  finiteness claim.
- `_double_with_preprojective_potential("A3")` (test helper): build the doubled Dynkin
  quiver (each arrow `a` plus its reverse `a*`) and `W = sum_a (a*a* - a**a)` as a
  `Potential`; the cyclic derivatives are exactly the preprojective mesh relations, so
  `Jac == Pi` degreewise (cross-engine). Read `preprojective.py:103-114` for the exact
  doubled-arrow naming so the two constructions align.
- The Markov / once-punctured-torus and Labardini annulus surface QPs are NOT pinned by a
  numeric dimension (their Jacobian-finiteness is a per-surface theorem this plan does not
  re-derive) — cite DWZ + Labardini for the framework, pin ONLY the hand-derived triangle
  (dim 6) and the preprojective cross-engine identity. The surface QPs land in P48 with
  their finiteness certificates.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/families/jacobian.py src/quiverlab/families/__init__.py \
        src/quiverlab/families/discover.py tests/families/test_jacobian.py
git commit -m "feat(families): Potential + JacobianAlgebra(Q, W) -- cyclic derivatives, dim-certified, NotFiniteDimensionalError on Jacobian-infinite"
```

---

### Task G: QPA battery + GUI/no-code story

**Files:**
- Create: `tests/qpa/test_tilting_qpa.py`
- Modify: `src/quiverlab/qpa/scripts.py` + `src/quiverlab/qpa/crosscheck.py`
  (`crosscheck_is_tilting` — probe-first)
- Modify: `src/quiverlab/hpc/spec.py` (`MODULE_KINDS`, `_MOD_REFS`, `_dispatch_module`:
  the `tilting_check` module kind), `docs/gui/runner.py` (the Pyodide twin:
  `_MODULE_KINDS`, `_MOD_REFS`, `_module_block`)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox, `S.ids`,
  `MOD_KIND_IDS`, push-list, block renderer), `webapp/static/app.js` (block builder),
  `webapp/server/i18n/en.json` + `es.json` (`mod.tilting_*` keys) +
  `webapp/templates/index.html` (`data-mod-tilting-*`)
- Modify: `src/quiverlab/trace/results_html.py` (`_HEADINGS` + `_block_html` branch)
- Modify: `src/quiverlab/families/__init__.py` top-level exports (already done in D/E/F),
  `webapp/server/catalog.py` (skip the non-scalar-arg constructors in `_iter_families`),
  `scripts/gui_build_hook.py::_preset_algebras()` (Jacobian + one-point drawable presets)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE new fixture, documented)
- Test: `tests/qpa/test_tilting_qpa.py`, `tests/webapp/test_tilting_check_p44.py`,
  `tests/gui/test_tilting_check_runner_twin.py`

**Interfaces:**
- QPA: probe live for the tilting surface via the `NamesGVars()` precedent
  (`tests/qpa/test_products_qpa.py`) — the manual lists `IsTiltingModule`,
  `TiltingModules`/`AllTiltingModules`; use what the probe finds. If a verb IS scripted,
  cross-check `is_tilting_module(T).is_tilting` against `IsTiltingModule(...)` over the
  zoo (kA2, kA3-rel, `line_abc_cde`); if NOT scripted, the battery honestly SKIPS and
  FAILS if the verb ever appears (the Plan-35 skip-that-fails-if-appears precedent).
  Approximations: probe for `MinimalRightApproximation`/`ApproximationsOfModule`; compare
  dims where exposed, honest skip otherwise.
- GUI module kind: `tilting_check` is a **schema-2 module kind** (the candidate `T` is the
  request `module` block; NO second module, so `estimator.sizing_dim` auto-sizes it — the
  estimator is kind-agnostic, no edit). It routes through `_dispatch_module` (NOT
  `_dispatch`), added to `MODULE_KINDS` + `_MOD_REFS` in BOTH `spec.py` and
  `docs/gui/runner.py`. Block shape:
  ```python
  {"kind": "tilting_check",
   "is_tilting": bool, "n": int, "pd": int | None,
   "self_ext_vanishes": bool, "num_summands": int, "num_vertices": int,
   "note": str,
   "references": [...], "citations": [...]}   # _MOD_REFS["tilting_check"] = ["bongartz_tilting","assem_book"]
  ```
  The 7 GUI touchpoints follow the `tau`-kind precedent exactly (agent-verified line
  anchors in the plan notes): checkbox `qlgui-tilting_check`, `S.ids`, `MOD_KIND_IDS`,
  the request push-list, the `renderBlock` branch (a small key/value table via the
  existing helpers), `app.js` `renderModuleBlocks` branch + `...Block()` builder, and the
  i18n `mod.tilting_*` chain (`en.json`+`es.json`+`index.html` `data-mod-tilting-*`+
  `app.js` `d.modTilting*`).
- Constructors as INPUT-side metagoal-1 features: `JacobianAlgebra`/`OnePointExtension`/
  `repetitive_slice` take `Algebra`/`Module`/`Potential` args (not scalar
  bool/int/str), so they do NOT fit the parameterized webapp family FORM
  (`catalog.py::_params_of`). Surface them as **drawable presets** (a small
  `Jac(triangle, abc)` and a small `A[M]`, each a genuine `kQ/I` rendered on the canvas)
  via `scripts/gui_build_hook.py::_preset_algebras()`, and register them in
  `families/discover.py` CATALOG for the `families()` listing while adding them to
  `catalog.py::_iter_families`'s skip-set (the `"zoo"` precedent) so the form builder does
  not try to introspect their non-scalar signatures.

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir; copy
  `tests/webapp/test_module_blocks_m0729.py`'s runner-pair fixture):

```python
# tests/webapp/test_tilting_check_p44.py
"""The tilting_check module kind: served by hpc.spec, mirrored by the Pyodide twin."""


def test_tilting_check_block_shape(tmp_path):
    # kA2 over QQ, module = the builtin regular module is not expressible as one builtin;
    # use module = {builtin: projective @ vertex 1} to get P1 -- NOT tilting (1 summand).
    from quiverlab.hpc.spec import run_spec, ComputeRequest
    req = _tilting_request(builtin=("projective", 1))     # helper: schema-2 module block
    out = run_spec(ComputeRequest.model_validate(req), tmp_path)
    b = out["results"]["tilting_check"]
    assert b["is_tilting"] is False and b["num_vertices"] == 2 and b["num_summands"] == 1


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py and assert json.dumps(sort_keys=True)
    # equality on the tilting_check block subkeys (both runners byte-identical), the way
    # tests/gui/test_yoneda_runner_twin.py does.
    ...
```

- [ ] **Step 2: Probe QPA + write the qpa battery** (concrete, in the
  `tests/qpa/test_tor_qpa.py` style; header
  `pytestmark = pytest.mark.skipif(session.should_skip_qpa(), ...)`): compare
  `is_tilting_module(T).is_tilting` vs the live `IsTiltingModule` verb over the zoo; the
  approximations probe skips honestly if unscripted, FAILING if the verb appears.

- [ ] **Step 3: Implement the handler + wiring.** In `spec.py::_dispatch_module`, after
  the `decompose` branch, add the `tilting_check` branch (call
  `tilting.is_tilting_module(M, n=item_n)`, wrap the `TiltingReport` fields into the
  block; catch the `decompose` char-caveat `QuiverlabError` into
  `{"error": "<the loud message>"}` per the Plan-30 τ-block honest-per-entry precedent);
  add `"tilting_check"` to `MODULE_KINDS` (`spec.py:70`) and `_MOD_REFS`
  (`spec.py:95`, and the runner twin `docs/gui/runner.py:273`). Mirror byte-for-byte in
  `docs/gui/runner.py::_module_block` after its `decompose` branch. Then the GUI
  touchpoints (checkbox / `S.ids` / `MOD_KIND_IDS` / push-list / `renderBlock` in both
  `gui.js` files, `app.js` builder, the i18n `mod.tilting_*` chain), the
  `results_html.py` `_HEADINGS` entry (`"tilting_check": "Tilting test"`) + `_block_html`
  branch, the `catalog.py` skip-set + `discover.py` CATALOG entries + the
  `_preset_algebras()` drawable presets.

- [ ] **Step 4: Add ONE golden fixture** `module_tilting_check_kA2` to
  `_runner_goldens.json`; note it in `test_runner_delegation.py`'s docstring change-log
  (the `products_loop_gf2` precedent — new fixture, existing entries byte-identical). Run
  the delegation test BEFORE adding it to confirm existing entries are untouched.

- [ ] **Step 5: Run the gates**

Run: `... -m pytest tests/webapp/test_tilting_check_p44.py tests/webapp/test_runner_delegation.py tests/gui/test_tilting_check_runner_twin.py tests/hpc -q`
then `... -m pytest tests/qpa/test_tilting_qpa.py -q -m qpa` (venv has [qpa]).
Expected: PASS; the approximations-vs-QPA probe SKIPS honestly.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc,qpa): tilting_check module kind + Jacobian/one-point presets + IsTiltingModule QPA battery"
```

---

### Task H: verification page, citations, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Test: existing release gates (`tests/release/test_oracle_classes.py`,
  `tests/citations/`)

- [ ] **Step 1: Citations.** Add (VERIFIED BibTeX only — `_r(...)` registry precedent at
  `registry.py:120+`; `_citation_pairs` swallows a `KeyError` but the plan must not rely
  on that): `bongartz_tilting` (Bongartz, *Tilted algebras*, LNM 903, 1981),
  `derksen_weyman_zelevinsky` (Derksen-Weyman-Zelevinsky, *Quivers with potentials and
  their representations I*, Selecta Math. 2008), `labardini` (Labardini-Fragoso, *Quivers
  with potentials associated to triangulated surfaces*, Proc. LMS 2009), and
  `hughes_waschbusche` (Hughes-Waschbüsch, *Trivial extensions of tilted algebras*, Proc.
  LMS 1983) — each BibTeX-verified before it ships. Reuse `assem_book` (ASS III/VI
  one-point extensions & tilting) and `happel_question` (Happel1989 change-of-rings LES)
  where they cover the definition; add a dedicated key only for a BibTeX-verifiable
  primary source.
- [ ] **Step 2: Verification page.** Add the Plan-44 subsystem rows:
  `modules/approximations.py` (`oracle_selfcert` approximation-property + minimality; the
  projective-cover tie `oracle_literature`); `modules/tilting.py` (`oracle_selfcert`
  regular-tilting + Bongartz self-certificate; `oracle_literature` ASS-VI APR tilt;
  `qpa` `IsTiltingModule` crosscheck); `core/basic.py` (`oracle_literature` `M₂(k)→k` +
  `kA2` round-trip; `oracle_crossengine` `End(⊕P_v)`-recovers-`kA_n` tied to
  `regular_corner_dims`; `oracle_selfcert` complete-orthogonal-idempotent + dim
  certificates); `families/one_point.py` (`oracle_literature` Cartan-block +
  `pd(S_ω)=pd_A(M)+1`); `families/repetitive.py` (`oracle_literature` `copies=1==A` + dim
  certificate); `families/jacobian.py` (`oracle_literature` triangle dim 6;
  `oracle_crossengine` preprojective-as-Jacobian; `oracle_selfcert` cyclic-derivative +
  `NotFiniteDimensionalError`). Add the **honest-scope entries**: (a) basic-ization /
  Gabriel recovery is rigorous only over char 0 / char > dim (the trace-form radical) and
  only when every Wedderburn block is split (a non-split division-algebra block refuses
  loudly); (b) `repetitive_slice` ships certified FINITE slices only — the full
  repetitive algebra is infinite-dimensional, and the `T(A) ≅ hat(A)/(nu)` orbit identity
  is a named successor, not shipped; (c) `JacobianAlgebra` refuses
  Jacobian-infinite inputs with `NotFiniteDimensionalError` (surface QPs land in P48);
  (d) tilting complement MUTATION is deferred (a named successor — `bongartz_completion`
  ships the one-complement case); (e) QPA has no approximations surface in scope (the
  probe skips honestly). Recount the class table (`tests/release/test_oracle_classes.py`
  drives the numbers — run collection, paste the LIVE counts, re-run to green; do NOT
  guess an at-authoring number given the mid-merge-train drift).
- [ ] **Step 3: README.** One features line: "tilting/cotilting + Bongartz completion,
  minimal add(M)-approximations, one-point extensions, repetitive slices,
  Jacobian algebras from a potential, and Gabriel-quiver recovery of any
  structure-constant algebra (End(M)/End(T) read back as kQ/I) — the C7 constructions
  toolkit."
- [ ] **Step 4: Full gate:**
  `... -m pytest tests/modules tests/families tests/core -q` (deep+fast, the touched
  dirs), `... -m pytest -q -m fast`, `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q` — all green.
- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-44 oracle rows + honest scope (basic-ization char/split caveat, repetitive finite-slice-only, Jacobian-infinite refusal, mutation deferred) + recounted classes"
```

---

## Acceptance (Plan-44 definition of done)

1. `right_add_approximation`/`left_add_approximation`, `is_tilting_module`/
   `is_cotilting_module`/`bongartz_completion`, `primitive_idempotents`/`basic_algebra`/
   `gabriel_quiver`/`presented_form`, `OnePointExtension`, `repetitive_slice`,
   `Potential`/`JacobianAlgebra`/`cyclic_derivative` all public, all CERTIFIED per
   instance (a dimension identity + a checkable structural oracle) or loudly refusing.
2. The crux Gabriel recovery is arbitrated by `M₂(k)→k`, the `kA2` round-trip identity,
   and `presented_form(End(⊕P_v))` recovering `kA_n` (tied to P37's
   `regular_corner_dims`); the char-caveat + non-split-block refusals are tested loud.
3. Bongartz completion is SELF-certified (`is_tilting_module(T ⊕ E)` True); the ASS-VI
   APR tilt is pinned; the `n=1` count criterion is applied only inside its theorem's
   hypotheses.
4. One-point extension satisfies the Cartan-block form and `pd(S_ω)=pd_A(M)+1`; the
   Jacobian triangle is hand-derived (dim 6), the preprojective-as-Jacobian identity is a
   live cross-engine oracle, and an infinite Jacobian refuses with
   `NotFiniteDimensionalError`.
5. QPA `IsTiltingModule` battery green live (`-m qpa`); the approximations probe skips
   honestly and FAILS if QPA ever ships a scripted surface.
6. `tilting_check` clickable end-to-end (GUI canvas → block → report) in EN+ES, both
   runners byte-identical, ONE golden added with a documented change-log entry, schema
   still v2; the `JacobianAlgebra`/`OnePointExtension` drawable presets render on the
   canvas.
7. `docs/verification.md` recounted (live numbers, mid-merge-train honest); README line
   added; deep (touched dirs) + fast + qpa + release + citations suites green. No
   dependency taken on P39/P40/P41 (this branch merged independently to `dev`).
8. Honest scope recorded on the verification page: basic-ization char/split caveat;
   repetitive finite-slice-only + the deferred `T(A)≅hat(A)/(nu)` orbit successor;
   Jacobian-infinite refusal (surface QPs → P48); tilting complement mutation deferred to
   a named successor.
