# Plan 47: Quasi-hereditary algebras & recollements — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The highest-weight-category toolkit a representation theorist reaches for once
modules, morphisms (P37) and constructions (P44) exist — **standard/costandard modules
Δ(i)/∇(i)** for a chosen vertex order; a **quasi-heredity test** per Dlab–Ringel with
three-valued honesty; **good (Δ-)filtration multiplicities** + the **BGG reciprocity**
oracle; the **characteristic tilting module** and its **Ringel dual**; and **recollements
from an idempotent** — the corner algebra `eAe`, the quotient `A/AeA`, and the **six
functors** as bimodule operations (these also hand P42's Grothendieck preset its worked
examples). Everything is CERTIFIED per instance (a dimension identity plus a checkable
structural oracle — Δ(i) has top S_i and End Δ(i)=k, P(i) is Δ-filtered, T is tilting with
`Ext¹(T,∇)=0`, the six functors satisfy their adjunction dim identities and the canonical
exact sequences) or REFUSES loudly (an algebra not quasi-hereditary in the given order, an
uncertifiable Δ-filtration greedy peel). **QPA has no quasi-hereditary / recollement
surface at all** — this is white space (stated on the verification page), so the oracle
class here is theory pins + internal self-certificates, not QPA agreement.

**Architecture:** Two new modules, each a thin exact-linear-algebra layer over EXISTING
primitives (no new math engines):

- `src/quiverlab/modules/quasihereditary.py` — `trace_module` (the sum-of-hom-images
  submodule, over `ModuleHom.image` + `column_space_pivots` + `radtopsoc.submodule`),
  `standard_modules`/`costandard_modules` (Δ = `P(i)/trace{P(j):j▷i}` via
  `radtopsoc.quotient`; ∇ = D of the opposite-algebra Δ via `duality.dualize` +
  `Algebra.opposite`), `is_quasi_hereditary` (`QHReport`, per-index brick + Δ-filtration
  certificates + the necessary finite-gl.dim check), `delta_multiplicities` (the greedy
  top-down Δ-peel with a loud uncertified verdict), `characteristic_tilting` (the
  Dlab–Ringel/Ringel universal-extension construction over `yoneda.baer_extension`,
  self-certified by `is_tilting_module` (P44) + `Ext¹(T,∇)=0`), and `ringel_dual`
  (`end_algebra(T).opposite()` presented via P44 `presented_form`, loud degrade otherwise).
- `src/quiverlab/modules/recollement.py` — `class Recollement(A, S)` building `eAe` on the
  corner path-type basis via structure constants (certified `dim = Σ` corner dims) and
  `A/AeA` on the complement subquiver (certified `dim = dim A − dim A e_S A`), plus the six
  functors `i_*`, `i^*`, `i^!`, `j^*`, `j_!`, `j_*` as methods. The self-certs — the four
  adjunction dim identities, the two canonical exact sequences at each joint, and
  `j^* j_! ≅ id`/`j^* j_* ≅ id` — ARE the test battery.

Thin `Algebra` wrappers (`standard_modules`, `is_quasi_hereditary`, `characteristic_tilting`,
`ringel_dual`, `corner_algebra`, `quotient_by_idempotent`, `recollement`) and ONE
algebra-scalar GUI compute kind `quasi_hereditary` (the `recognizers`/`ext_algebra`
scalar-kind precedent) round it out. The whole thing composes on the P37 hub
(`ModuleHom`/`hom_basis`/`image`/`kernel`/`cokernel`/`direct_sum`/`is_isomorphic`/
`end_algebra`), P30 `decompose`/`identify_standard`, the P44 constructions
(`is_tilting_module`, `presented_form`), `modules/tor.py`'s induction pattern
(`_induced`/`_vertex_basis`/`_left_action`), and the presented-algebra backbone
(`Quiver.algebra`, `invariants/pathbasis.py::path_type_basis`,
`families/trivial_extension.py`'s length-lex relation extraction).

**Tech Stack:** pure exact linear algebra over `Domain` (`modules/linalg_mod`,
`fields/linalg`); sympy only for the Ringel-dual Cartan Smith-normal-form oracle
(`sympy.Matrix.invariant_factors`, in the TEST, not `src/`). No floats in `src/`
(AST-gated by `tests/test_no_floats.py`).

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P37 (categorical glue) and P44 (C7 constructions) are the hard prerequisites and are
  MERGED to `dev` before this plan executes.** This plan consumes `modules/morphism.py`
  (`ModuleHom` with `.image()`/`.kernel()`/`.cokernel()`/`.then()`/`.is_epi`/`.is_mono`,
  `hom_basis`, `direct_sum`, `is_direct_summand`), `modules/hom.py`
  (`hom_space`/`is_isomorphic`/`identify_standard`/`_assert_comparable`),
  `modules/endomorphism.py::end_algebra`, `modules/yoneda.py::baer_extension`,
  `modules/duality.py::dualize`, `modules/radtopsoc.py` (`submodule`/`quotient`/
  `top`/`radical`/`socle`), `modules/ext.py` (`ext_dims`, `global_dimension`,
  `GlobalDimension`), and from **P44**: `modules/tilting.py::is_tilting_module` and
  `core/basic.py::presented_form`. **The Interfaces of P44 are frozen in
  `docs/plans/2026-08-05-plan-44-constructions.md`** (`is_tilting_module` lines ~367–381,
  `primitive_idempotents`/`basic_algebra`/`gabriel_quiver`/`presented_form` lines
  ~605–616) — read them there, do not re-derive. Branch `plan-47-quasihereditary` off
  `dev` **after P44 has merged**.
- Buckets are auto-assigned by directory (`tests/conftest.py`): `tests/modules/` → **deep**,
  `tests/core/`/`tests/invariants/` → **fast**, `tests/webapp/`/`tests/gui/` → **fast**,
  `tests/qpa/` → **qpa**. **All Δ/∇/qh/tilting/recollement tests live in `tests/modules/`
  (deep).** Run new tests by path during development; finish each task with a `-m deep`
  spot-run of the touched files.
- **The char scope is NARROW here and must be stated precisely — most of this plan is
  char-clean.** Δ/∇ construction (`trace_module` + quotient), `is_quasi_hereditary`
  (`hom_space` brick test + dim-vector greedy peel + `GlobalDimension`), `delta_multiplicities`,
  and the entire `Recollement` (corner structure constants, `A/AeA`, six functors, adjunction
  dims) are **pure exact linear algebra over any exact `Domain`** — no `decompose`, no
  trace-form radical, NO char caveat. The caveat bites in exactly TWO places, both inherited:
  (1) the `is_tilting_module(T)` self-certificate inside `characteristic_tilting` uses the
  P30 `decompose` summand COUNT (rigorous over char 0 / char > dim), and (2) a **presented**
  Ringel dual goes through P44 `presented_form` (rigorous over char 0 / char > dim, split
  Wedderburn blocks). So the tilting/Ringel-dual batteries run over **QQ** (or a large prime
  `GF(32003)`), and off-scope those two functions inherit P30/P44's loud `QuiverlabError`
  refusal — never a silent wrong T or dual. Δ/∇/qh/recollement have NO such restriction and
  their batteries include a `GF(2)` cell to prove it.
- **Order convention (pinned once, verbatim in every docstring).** `order` is a sequence of
  the vertices listed **lowest → highest**; `rank(v)` = its index in that list; `j ▷ i`
  ("j above i") means `rank(j) > rank(i)`. `order=None` is the **natural order**
  `sorted(A.quiver.vertices)` (ascending vertex labels). `Δ(i) = P(i) / trace{P(j) : j ▷ i}`.
  **Quasi-heredity is order-dependent** (a fact the GUI note states); the Dlab–Ringel theorem
  that a hereditary algebra is quasi-hereditary for EVERY order is a two-orders oracle (Task C).
- **The `eAe`-vs-subquiver trap is load-bearing and must be stated in the module docstring +
  the verification page.** For a vertex subset `S`, the corner algebra `eAe` (`e = Σ_{v∈S} e_v`)
  is **NOT** the full-subquiver algebra on `S`: a path `p` may leave `S` in its interior yet
  satisfy `e p e = p ≠ 0` (worked example: `kA₃` with `S={1,3}`, the path `ab: 1→3` travels
  through vertex `2 ∉ S` but `e(ab)e = ab`, so `eAe = kA₂` of dim 3, not the subquiver `k×k`
  of dim 2). Build `eAe` on the corner path-type basis via **structure constants**, certified
  `dim = Σ_{v,w∈S} dim e_v A e_w`. Conversely `A/AeA` IS presentable on the complement full
  subquiver, certified by `dim A/AeA = dim A − dim(A e_S A)` computed from the path-basis span.
- **The honesty pattern is mandatory** (the `GlobalDimension`/`TiltingReport` three-valued
  precedent): every construction is CERTIFIED per instance (a dim identity + a structural
  oracle) or REFUSES loudly. `is_quasi_hereditary` returns a `QHReport` naming the failing
  clause (never a bare `False`); `delta_multiplicities` returns `certified=False` (loud) when
  the greedy peel cannot realize a Δ-filtration; `characteristic_tilting` raises when the
  self-cert (`is_tilting_module(T)` True AND `Ext¹(T,∇(j))=0`) fails; `ringel_dual` degrades
  from `kQ/I` to the structure-constant form loudly (a `note`), never silently.
- **Composition is left-to-right** (`f.then(g)`; `a*b` = first `a` then `b`, ASS). Never
  overload `*` on morphisms. `_assert_comparable` guards every cross-module/cross-side/
  cross-algebra Hom.
- Plan-32 markers: the Δ(i)-top / End Δ(i)=k / P(i)-Δ-filtered / adjunction-dim /
  exact-sequence-at-each-joint / `j^*j_!≅id` certificates and the tilting/Ringel-dual
  self-certs = `oracle_selfcert`; the Dlab–Ringel theory pins (`kA_n` natural-order Δ=S_i,
  the two-orders quasi-heredity, the `k[x]/(x²)` negative, the hand-derived `T`), the BGG
  reciprocity battery, and the double-Ringel-dual Cartan Smith-form identity =
  `oracle_literature`; any two-independent-route dim agreement (e.g. `dim ∇(j)_i =
  (P(i):Δ(j))`) = `oracle_crossengine`. **No `tests/qpa/` here** — QPA white space.
- **Mid-merge-train counts:** v0.2.0 lands ~15 subplans in overlapping waves, so absolute
  suite counts drift between authoring and merge. **Task G recounts the oracle-class table
  at merge time** by running `tests/release/test_oracle_classes.py` (paste the live numbers,
  never a guessed-at-authoring count) and claims only the deltas this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + honest-scope entries +
  recounted class table green) and the README line. Conventional commits; green at every commit.

---

### Task A: `trace_module` — the sum-of-hom-images submodule

The one net-new primitive Δ needs: the trace `tr_{add S}(M) = Σ_{N∈S, f∈Hom(N,M)} im f`, an
A-submodule of M. Every other Δ/∇/qh/recollement piece is composed from it and P37.

**Files:**
- Create: `src/quiverlab/modules/quasihereditary.py`
- Test: `tests/modules/test_quasihereditary_trace.py`

**Interfaces:**
- Consumes: `morphism.py::hom_basis(N, M) -> list[ModuleHom]`,
  `ModuleHom.image() -> (I, epi, mono)` (`morphism.py:115` — `mono: I >-> M`, its columns
  are the image basis expressed in `M`'s coordinates), `radtopsoc.py::submodule(M, cols,
  name, side)` (`radtopsoc.py:27`), `linalg_mod::column_space_pivots`/`cols_to_matrix`/`col`,
  `hom.py::_assert_comparable`.
- Produces:
  ```python
  def trace_module(sources, M, name="trace") -> tuple["Module", "ModuleHom"]:
      # tr(M) = sum over N in `sources`, f in Hom_A(N, M) of im(f): the smallest A-submodule
      # of M containing every hom image (A-stable, being a sum of module-map images).
      # Returns (T, iota) with iota: T >-> M the (mono) inclusion; T.dim == 0 gives the
      # zero submodule + a 0-column mono. Every N must be comparable to M (same algebra/side).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_quasihereditary_trace.py
"""The trace submodule tr_{add S}(M) (Plan 47). Self-certifying: iota is mono, its image is
A-stable, every Hom(N,M) factors through iota, and tr is idempotent-monotone. Over QQ AND
GF(2) -- trace_module is pure linear algebra, no char caveat."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.morphism import hom_basis
from quiverlab.modules.quasihereditary import trace_module

pytestmark = pytest.mark.oracle_selfcert


def _a3(field=QQ):
    return linear_path_algebra(3, field=field)


def test_trace_of_P2_in_P1_is_the_radical():
    # kA3 (1->2->3): P1 = [1..3], P2 = [2..3]. Hom(P2, P1) = P1 e_2 = e_1 A e_2 = {a} (1-dim),
    # its image is the submodule generated by a = [2..3] = rad P1. So tr_{P2}(P1) = rad P1.
    A = _a3()
    P1, P2 = A.projective(1), A.projective(2)
    T, iota = trace_module([P2], P1)
    assert iota.is_mono()
    assert T.dimension_vector() == P1.radical().dimension_vector()   # {2:1, 3:1}


def test_trace_factors_every_hom():
    A = _a3()
    P1 = A.projective(1)
    sources = [A.projective(2), A.projective(3)]
    T, iota = trace_module(sources, P1)
    # every h: N -> P1 (N in sources) lands inside im(iota): solve iota . x = h columnwise.
    dom = A.domain
    imgcols = [lm.col(iota.matrix, j) for j in range(T.dim)]
    B = lm.cols_to_matrix(imgcols) if imgcols else []
    for N in sources:
        for h in hom_basis(N, P1):
            for j in range(N.dim):
                col = lm.col(h.matrix, j)
                assert lm.solve_columns(B, lm.cols_to_matrix([col]), dom) is not None


def test_trace_is_A_stable_over_GF2():
    # pure linear algebra: no char restriction. tr is a genuine submodule (action-closed).
    A = _a3(GF(2))
    P1 = A.projective(1)
    T, iota = trace_module([A.projective(2)], P1)
    # A-stability: for each action label b, b maps im(iota) into itself -> the composite
    # iota.then(action-as-hom) stays in im(iota). Checked structurally by submodule(): if it
    # returned, the span was closed. Assert the recorded dim matches rad P1.
    assert T.dim == 2


def test_empty_sources_is_zero_submodule():
    A = _a3()
    P1 = A.projective(1)
    T, iota = trace_module([], P1)
    assert T.dim == 0 and iota.is_mono()


def test_cross_algebra_refused():
    from quiverlab.errors import QuiverlabError
    A, B = _a3(), _a3()
    with pytest.raises(QuiverlabError):
        trace_module([A.projective(2)], B.projective(1))
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_quasihereditary_trace.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.modules.quasihereditary`

- [ ] **Step 3: Implement**

```python
"""Quasi-hereditary structure: standard/costandard modules, the quasi-heredity test,
good-filtration multiplicities, the characteristic tilting module and its Ringel dual
(Plan 47). Right modules; exact over any Domain. Δ/∇/qh/filtration are char-clean pure
linear algebra; only the tilting summand COUNT (P30 decompose) and a PRESENTED Ringel dual
(P44 presented_form) inherit the char 0 / char > dim caveat. Order convention: `order` is a
sequence of vertices lowest->highest, rank = index, j ▷ i iff rank(j) > rank(i),
Δ(i) = P(i)/trace{P(j): j ▷ i} (Dlab-Ringel, Illinois J. Math. 1989)."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules import radtopsoc
from quiverlab.modules.builders import _require_provenance
from quiverlab.modules.hom import _assert_comparable, hom_space, is_isomorphic
from quiverlab.modules.morphism import direct_sum, hom_basis


def trace_module(sources, M, name="trace"):
    dom = M.domain
    cols = []
    for N in sources:
        _assert_comparable(N, M, "trace")
        for f in hom_basis(N, M):
            _I, _epi, mono = f.image()               # mono: im(f) >-> M, cols in M-coords
            for j in range(_I.dim):
                cols.append(lm.col(mono.matrix, j))
    from quiverlab.modules.morphism import ModuleHom
    if not cols:
        T = radtopsoc.submodule(M, [], name=name, side=M.side)
        return T, ModuleHom(T, M, lm.zeros(M.dim, 0, dom), check=False)
    piv = lm.column_space_pivots(lm.cols_to_matrix(cols), dom)
    basis = [cols[j] for j in piv]
    T = radtopsoc.submodule(M, basis, name=name, side=M.side)     # A-stable span
    iota = ModuleHom(T, M, lm.cols_to_matrix(basis), check=False)
    return T, iota
```

**Adjust to reality (Task A):**
- `ModuleHom.image()` returns `(I, epi, mono)` with `mono.matrix` = the image basis columns
  in `M`'s coordinates (`morphism.py:124`). The union of those over all `(N, f)` spans the
  trace; `column_space_pivots` deduplicates to an independent basis.
- `radtopsoc.submodule(M, basis_cols, name, side)` expects an **action-closed** spanning set.
  The trace IS action-closed (a sum of module-map images, each a submodule), so the
  column-space basis is a valid submodule generating set. If `submodule` re-verifies closure
  (read `radtopsoc.py:27`), it passes; if it CLOSES a raw set it is idempotent on ours. Pass
  `side=M.side` so the returned `T` is a module on the same side as `M`.
- `_assert_comparable(N, M, "trace")` fires the cross-algebra/cross-side refusal BEFORE any
  `hom_basis` work — the `test_cross_algebra_refused` arbiter.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/quasihereditary.py tests/modules/test_quasihereditary_trace.py
git commit -m "feat(modules): trace_module -- sum-of-hom-images A-submodule (Δ's building block, char-clean)"
```

---

### Task B: standard & costandard modules Δ(i)/∇(i)

**Files:**
- Modify: `src/quiverlab/modules/quasihereditary.py`
- Test: `tests/modules/test_quasihereditary_standard.py`

**Interfaces:**
- Consumes: Task A (`trace_module`), `builders.projective(A, v)`,
  `radtopsoc.quotient(M, sub_cols, name, side) -> (Q, proj)` (`radtopsoc.py:47`),
  `Algebra.opposite()` (`core/algebra.py:398`), `duality.dualize(M)` (side-flipping D,
  `duality.py:31`), `builders.injective`, `Module.top()`/`.dimension_vector()`/
  `.composition_factors()`.
- Produces:
  ```python
  def _order_ranks(A, order) -> tuple[dict, list]      # (rank map, order list); validated
  def standard_modules(A, order=None) -> dict          # vertex -> Δ(i) Module (right A-mod)
  def costandard_modules(A, order=None) -> dict         # vertex -> ∇(i) = D(Δ_{A^op}(i))
  ```

**Math pinned (the hand derivation — write it verbatim in the docstring).**
`kA_n = linear_path_algebra(n)` has arrows `a_i: i→i+1` (forward), RIGHT modules, so
`P(v) = e_v A` = paths STARTING at `v` = the interval module `[v..n]`, and
`Hom(P(j),P(i)) = P(i)e_j = e_i A e_j` = paths `i→j`, nonzero iff `i ≤ j`, with image the
interval submodule `[j..n] ⊆ [i..n]`.

- **Natural order** `1◁2◁…◁n` (so `rank(v)=v`): `trace{P(j):j▷i} = Σ_{j>i}[j..n] = [i+1..n]
  = rad P(i)`, hence **`Δ(i) = P(i)/rad P(i) = S_i`** (the simple = degenerate interval
  `[i..i]`). Dually `∇(i) = D(Δ_{A^op}(i)) = D(A e_i) = I(i) = [1..i]` (the interval from 1
  to `i`, dim `i`, socle `S_i`, top `S_1`). BGG: `(P(i):Δ(j)) = [i≤j] = [∇(j):S(i)]`. This is
  the order the GUI reports.
- **Opposite order** `n◁…◁2◁1` (so `rank(v) = n−v`): now `j▷i ⟺ j<i` and `Hom(P(j),P(i))≠0
  ⟺ i≤j` fails for `j<i`, so `trace = 0` and **`Δ(i) = P(i) = [i..n]`** (genuine interval
  modules, non-simple), `∇(i) = S_i`. Both orders are quasi-hereditary (Dlab–Ringel) — the
  Task-C two-orders oracle. And `Δ(i)=P(i)` ⇒ `T = A` (Task D).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_quasihereditary_standard.py
"""Standard/costandard modules Δ(i)/∇(i) (Plan 47). Literature (Dlab-Ringel): kA_n natural
order Δ(i)=S_i, ∇(i)=[1..i]; opposite order Δ(i)=P(i). Self-cert: top Δ(i)=S_i, [Δ(i):S(i)]=1,
socle ∇(i)=S_i. Char-clean (a GF(2) cell)."""
import pytest

from quiverlab import GF, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.hom import is_isomorphic
from quiverlab.modules.quasihereditary import costandard_modules, standard_modules

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _a3(field=QQ):
    return linear_path_algebra(3, field=field)


@lit
def test_natural_order_standards_are_simple():
    A = _a3()
    D = standard_modules(A)                       # order=None -> natural
    for v in (1, 2, 3):
        assert D[v].dim == 1                       # Δ(v) = S_v
        assert is_isomorphic(D[v], A.simple(v))


@lit
def test_natural_order_costandards_are_intervals():
    A = _a3()
    N = costandard_modules(A)
    for v in (1, 2, 3):
        assert N[v].dimension_vector() == {w: (1 if w <= v else 0) for w in (1, 2, 3)}
        assert is_isomorphic(N[v], A.injective(v))  # ∇(v) = I(v) = [1..v]


@lit
def test_opposite_order_standards_are_projectives():
    A = _a3()
    order = [3, 2, 1]                              # n ◁ ... ◁ 1
    D = standard_modules(A, order)
    for v in (1, 2, 3):
        assert is_isomorphic(D[v], A.projective(v))  # Δ(v) = P(v) = [v..3]


@selfcert
def test_delta_top_is_simple_and_mult_one():
    A = _a3()
    for order in (None, [3, 2, 1], [2, 1, 3]):
        D = standard_modules(A, order)
        for v in (1, 2, 3):
            assert D[v].top().dimension_vector() == {w: (1 if w == v else 0)
                                                     for w in (1, 2, 3)}   # top Δ(i)=S_i
            cf = D[v].composition_factors()
            assert cf.count(v) == 1                # [Δ(i):S(i)] = 1


@selfcert
def test_standard_char_clean_gf2():
    A = _a3(GF(2))
    D = standard_modules(A)
    assert all(D[v].dim == 1 for v in (1, 2, 3))   # no char caveat on Δ


@selfcert
def test_bad_order_refused():
    from quiverlab.errors import QuiverlabError
    A = _a3()
    with pytest.raises(QuiverlabError):
        standard_modules(A, [1, 2])                # not a permutation of the vertices
```

- [ ] **Step 2: Run to verify failure** — `ImportError: standard_modules`

- [ ] **Step 3: Implement**

```python
def _order_ranks(A, order):
    _require_provenance(A, "standard_modules")
    verts = list(A.quiver.vertices)
    if order is None:
        order = sorted(verts, key=lambda v: (isinstance(v, str), v))
    order = list(order)
    if sorted(order, key=lambda v: (isinstance(v, str), str(v))) != \
            sorted(verts, key=lambda v: (isinstance(v, str), str(v))):
        raise QuiverlabError(
            f"quasi-hereditary order {order!r} is not a permutation of the vertices {verts!r}",
            hint="pass each vertex exactly once, lowest->highest; None = natural order")
    return {v: r for r, v in enumerate(order)}, order


def standard_modules(A, order=None):
    ranks, order = _order_ranks(A, order)
    P = {v: A.projective(v) for v in order}
    out = {}
    for i in order:
        higher = [P[j] for j in order if ranks[j] > ranks[i]]
        T, iota = trace_module(higher, P[i], name=f"tr>{i}")
        subcols = [lm.col(iota.matrix, j) for j in range(T.dim)]
        Q, _proj = radtopsoc.quotient(P[i], subcols, name=f"Delta_{i}", side=P[i].side)
        out[i] = Q
    return out


def costandard_modules(A, order=None):
    from quiverlab.modules.duality import dualize
    Aop = A.opposite()
    dstd = standard_modules(Aop, order)             # Δ over the opposite algebra
    return {i: dualize(dstd[i]) for i in dstd}       # ∇_A(i) = D(Δ_{A^op}(i))
```

**Adjust to reality (Task B):**
- `radtopsoc.quotient(M, sub_cols, name, side)` returns `(Q, proj)` with `proj: M ->> Q`
  (`radtopsoc.py:47`); pass the trace's `iota.matrix` columns as `sub_cols`.
- `A.opposite()` reverses the quiver and transposes the structure constants (P23/P24), so its
  `projective(v)` and the same vertex `order` are well-defined; `dualize` (P24, side-aware D)
  turns the right `A^op`-module `Δ_{A^op}(i)` into a right `A`-module `∇(i)`. The `∇(v) = I(v)`
  and `socle ∇(v) = S_v` self-certs (Task B tests) ARBITRATE the side bookkeeping — if `dualize`
  landed on the wrong side, `is_isomorphic(N[v], A.injective(v))` fails; fix the side tag once
  and pin it (the P24 "side is a presentation tag" precedent), never weaken the test.
- The costandard route reuses `standard_modules` on `A.opposite()` — do NOT re-implement a
  dual trace; the D-of-opposite identity is exactly the Dlab–Ringel definition of ∇.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/quasihereditary.py tests/modules/test_quasihereditary_standard.py
git commit -m "feat(modules): standard/costandard modules Δ(i)/∇(i) -- trace-quotient + D-of-opposite, order-parametrized"
```

---

### Task C: `is_quasi_hereditary` (`QHReport`) + `delta_multiplicities` + BGG reciprocity

**Files:**
- Modify: `src/quiverlab/modules/quasihereditary.py`
- Test: `tests/modules/test_quasihereditary_test.py`

**Interfaces:**
- Consumes: Task B (`standard_modules`, `costandard_modules`, `_order_ranks`),
  `hom.py::hom_space` (dim End Δ(i) = brick test), `ext.py::global_dimension(A) ->
  GlobalDimension` (`.exact` = finite gl.dim, `ext.py:85`), `builders.projective`,
  `Module.top()`/`.dimension_vector()`/`.composition_factors()`,
  `morphism.py::hom_basis` (surjection search), `ModuleHom.kernel()`.
- Produces:
  ```python
  @dataclass
  class QHReport:
      is_quasi_hereditary: bool
      order: list
      gl_dim: object          # GlobalDimension (qh => finite, Dlab-Ringel)
      per_index: dict         # vertex -> {"brick": bool, "delta_filters_P": bool, "note": str}
      note: str               # names the failing clause, "quasi-hereditary" when ok
      # __bool__ -> is_quasi_hereditary ; __repr__ names the failing clause
  def is_quasi_hereditary(A, order=None) -> QHReport
  def delta_multiplicities(M, deltas, order=None) -> tuple[dict, bool]
      # greedy top-down Δ-peel: (mult: vertex -> (M:Δ(i)), certified). certified=False (LOUD)
      # when the peel cannot realize a Δ-filtration -- three-valued honesty like TiltingReport.
  ```

**Certification (Dlab–Ringel, Illinois J. Math. 33 (1989)):** `A` with order `▷` is
quasi-hereditary iff for every `i`: (1) `End_A(Δ(i)) = k` (Δ(i) is a brick — `dim hom_space
= 1`), and (2) `P(i)` has a Δ-filtration whose top layer is `Δ(i)` and the rest are `Δ(j)`,
`j ▷ i`. A **necessary** classical consequence is `gl.dim A < ∞` (checked via
`GlobalDimension.exact`, cited). The report is three-valued: `True`, or `False` with `note`
naming the first failing clause, or (inside `delta_multiplicities`) an uncertified peel.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_quasihereditary_test.py
"""Quasi-heredity test + Δ-multiplicities + BGG reciprocity (Plan 47). Literature: kA_n is
quasi-hereditary for BOTH the natural and opposite orders (Dlab-Ringel: hereditary => qh for
every order); k[x]/(x^2) is NOT (infinite gl.dim / End Δ != k). BGG: (P(i):Δ(j))=[∇(j):S(i)].
Char-clean (a GF(2) cell)."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.quasihereditary import (costandard_modules, delta_multiplicities,
                                               is_quasi_hereditary, standard_modules)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _a3(field=QQ):
    return linear_path_algebra(3, field=field)


def _dual_numbers(field=QQ):
    # k[x]/(x^2): one vertex, one loop x, x^2 = 0. Local, self-injective, gl.dim = infinity.
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=field)


@lit
@pytest.mark.parametrize("order", [None, [3, 2, 1]])
def test_ka3_is_quasi_hereditary_for_two_orders(order):
    A = _a3()
    rep = is_quasi_hereditary(A, order)
    assert rep.is_quasi_hereditary is True and bool(rep) is True
    assert rep.gl_dim.exact and int(rep.gl_dim) == 1          # hereditary
    assert all(rep.per_index[v]["brick"] for v in (1, 2, 3))
    assert all(rep.per_index[v]["delta_filters_P"] for v in (1, 2, 3))


@lit
def test_dual_numbers_not_quasi_hereditary():
    A = _dual_numbers()
    rep = is_quasi_hereditary(A)
    assert rep.is_quasi_hereditary is False and bool(rep) is False
    # Δ(1) = P(1) = k[x]/(x^2) has End = k[x]/(x^2) (dim 2 != 1) AND gl.dim is infinite;
    # the note names a concrete failing clause (not a bare False).
    assert rep.note and ("gl.dim" in rep.note or "End" in rep.note or "1" in rep.note)


@selfcert
def test_qh_char_clean_gf2():
    A = _a3(GF(2))
    assert is_quasi_hereditary(A).is_quasi_hereditary is True   # no char caveat on the test


@selfcert
def test_delta_filtration_of_P_is_certified():
    A = _a3()
    D = standard_modules(A)                        # natural order: Δ(i)=S_i
    for v in (1, 2, 3):
        mult, ok = delta_multiplicities(A.projective(v), D)
        assert ok is True
        # P(v) = [v..3] is uniserial with factors S_v..S_3, each Δ=S once:
        assert mult == {w: (1 if w >= v else 0) for w in (1, 2, 3)}


@selfcert
def test_non_filtered_module_is_uncertified_loudly():
    # a module whose top is not a Δ-top cannot be Δ-peeled: certified=False, not a wrong count.
    A = _a3()
    D = {1: A.simple(1)}                            # deliberately incomplete Δ set
    mult, ok = delta_multiplicities(A.projective(1), D)   # P(1) needs Δ(2),Δ(3) too
    assert ok is False


@xeng
def test_bgg_reciprocity():
    A = _a3()
    for order in (None, [3, 2, 1]):
        D = standard_modules(A, order)
        N = costandard_modules(A, order)
        for i in (1, 2, 3):
            mult, ok = delta_multiplicities(A.projective(i), D, order)
            assert ok
            for j in (1, 2, 3):
                # [∇(j):S(i)] = multiplicity of vertex i in ∇(j)'s composition factors
                comp = N[j].composition_factors().count(i)
                assert mult[j] == comp, (order, i, j, mult[j], comp)
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```python
def _hom_dim(M, N):
    return len(hom_space(M, N))


def delta_multiplicities(M, deltas, order=None):
    """Greedy top-down Δ-peel. At each step take the HIGHEST-in-order Δ(i) whose top S_i
    appears in top(M) and which admits an epi M ->> Δ(i); quotient by its kernel; recurse.
    certified=False (LOUD) if the peel stalls with M != 0 (no Δ-filtration)."""
    A = deltas[next(iter(deltas))].base_algebra if deltas else M.base_algebra
    ranks, order = _order_ranks(A, order)
    mult = {i: 0 for i in deltas}
    cur = M
    guard = 0
    while cur.dim > 0:
        guard += 1
        if guard > M.dim + 1:                        # cannot exceed dim M genuine layers
            return mult, False
        tops = cur.top().dimension_vector()
        cands = sorted((i for i in deltas if tops.get(i, 0) > 0),
                       key=lambda i: ranks.get(i, -1), reverse=True)
        peeled = False
        for i in cands:
            epi = _find_surjection(cur, deltas[i])
            if epi is None:
                continue
            K, _iota = epi.kernel()
            mult[i] += 1
            cur = K
            peeled = True
            break
        if not peeled:
            return mult, False                        # uncertified: no Δ peels off the top
    return mult, True


def _find_surjection(M, D):
    """An epi M ->> D, or None. Since top(D) = S_i and S_i occurs in top(M), a hom realizing
    D as a quotient exists iff D is a top quotient of M; scan the Hom basis and small
    combinations for a surjective one (rank(matrix) == D.dim)."""
    dom = M.domain
    homs = hom_basis(M, D)
    for f in homs:
        if f.is_epi():
            return f
    # try combinations coordinate-searched cheaply (the Δ's are small); combine basis homs.
    ...  # see Adjust to reality: a bounded search over {0,1,-1} (and field units) combos;
        # for the directed/uniserial oracles a single basis hom already surjects.
    return None


def is_quasi_hereditary(A, order=None):
    from quiverlab.modules.ext import global_dimension
    ranks, order = _order_ranks(A, order)
    deltas = standard_modules(A, order)
    gld = global_dimension(A)                         # qh => finite (Dlab-Ringel)
    per, ok = {}, True
    for i in order:
        brick = (_hom_dim(deltas[i], deltas[i]) == 1)     # End Δ(i) = k
        _mult, filt = delta_multiplicities(A.projective(i), deltas, order)
        note = ("ok" if (brick and filt) else
                ("End Δ != k" if not brick else "P has no Δ-filtration"))
        per[i] = {"brick": brick, "delta_filters_P": filt, "note": note}
        ok = ok and brick and filt
    ok = ok and gld.exact
    if ok:
        note = "quasi-hereditary"
    elif not gld.exact:
        note = f"gl.dim not finite ({gld!r})"
    else:
        note = "; ".join(f"{i}: {per[i]['note']}" for i in order if per[i]["note"] != "ok")
    return QHReport(ok, order, gld, per, note)
```

**Adjust to reality (Task C):**
- `_find_surjection`: for the directed/uniserial oracles (`kA_n`, both orders) a single basis
  hom already surjects onto Δ(i); implement the general case as a bounded search over small
  field-coefficient combinations of `hom_basis(M, D)` (units first) testing `.is_epi()`, and
  if none is found return `None` (the honest "cannot peel Δ(i) here"). Do NOT fabricate a
  surjection — the `certified=False` verdict is the loud contract (`test_non_filtered_module_
  is_uncertified_loudly`). The greedy top-down order (highest Δ first) is the ARBITER: if it
  produced a wrong count, `test_delta_filtration_of_P_is_certified` fails — re-check the order
  direction, do not weaken.
- `deltas[i].base_algebra` recovers the algebra without an extra argument (P24 accessor,
  `module.py:50`). `delta_multiplicities` takes `deltas` (not the algebra) so it also serves
  the BGG oracle and Task D's filtration checks with a possibly-restricted Δ set.
- The `gl.dim finite` clause is NECESSARY, not sufficient; it is ANDed after the per-index
  brick + Δ-filtration certificates so a finite-gl.dim non-qh order still fails on a clause the
  `note` names. Cite `dlab_ringel` for `qh ⇒ gl.dim < ∞` in the docstring.
- `QHReport` mirrors `TiltingReport`'s three-valued `__bool__`/`__repr__` (P44, `tilting.py`);
  `__repr__` returns `"quasi-hereditary (order …)"` or `"not quasi-hereditary: {note}"`.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/quasihereditary.py tests/modules/test_quasihereditary_test.py
git commit -m "feat(modules): is_quasi_hereditary (QHReport) + delta_multiplicities greedy peel + BGG reciprocity oracle"
```

---

### Task D: characteristic tilting module + Ringel dual

**Files:**
- Modify: `src/quiverlab/modules/quasihereditary.py`
- Test: `tests/modules/test_quasihereditary_ringel.py`

**Interfaces:**
- Consumes: Task B/C (`standard_modules`, `costandard_modules`, `_order_ranks`,
  `is_quasi_hereditary`), `ext.py::ext_dims(A, M, N, 1, with_reps=True)` (Ext¹ + cocycle
  reps), `yoneda.py::baer_extension(M, N, cocycle) -> YonedaSequence` (`.middle`),
  `morphism.py::direct_sum`, `endomorphism.py::end_algebra(M)` (P37 — presentation-less
  `End_A(M)`), `Algebra.opposite()`, **P44** `tilting.py::is_tilting_module(T) ->
  TiltingReport` and `core/basic.py::presented_form(A) -> Algebra`.
- Produces:
  ```python
  def characteristic_tilting(A, order=None) -> "Module"
      # T = ⊕_i T(i), the Dlab-Ringel characteristic tilting module: T(i) built inductively
      # low->high by universal extensions of Δ(i) through the lower T(j) via baer_extension.
      # SELF-CERTIFIED: Ext^1(T, ∇(j)) = 0 for all j AND is_tilting_module(T).is_tilting.
  def ringel_dual(A, order=None) -> "Algebra"
      # R(A) = End_A(T)^op, presented as kQ/I via P44 presented_form when char permits;
      # loud degrade to the structure-constant form (a `note`) otherwise.
  ```

**Math pinned (the hand derivation — verbatim in the docstring).** `T(i)` is the unique
indecomposable in `F(Δ) ∩ F(∇)` with `Δ(i)` as its bottom Δ-layer; Ringel's construction
builds it as the universal extension `0 → (⊕_{j▷... } T(j)^{d_j}) → T(i) → Δ(i) → 0` with
`d_j = dim Ext¹(Δ(i), T(j))`, iterated over the T(j) already built for `j ◁ i`. For
`kA_n` natural order (`Δ(i)=S_i`, `∇(i)=I(i)`): every module is trivially Δ-filtered
(Δ = simples), so `T(i) = ∇(i) = I(i)` and **`T = ⊕_i I(i) = D(A)`** (the injective
cogenerator). For `kA_n` opposite order (`Δ(i)=P(i)`): every module is Δ-filtered by
projectives, so `T(i) = P(i)` and **`T = A`** (the regular module). `R(A) = End_A(T)^op`;
by Ringel's theorem the double dual `R(R(A))` is Morita-equivalent to `A`, so their Cartan
matrices have equal invariant factors (Smith normal form) — the safe cross-check.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_quasihereditary_ringel.py
"""Characteristic tilting module + Ringel dual (Plan 47). Literature (Ringel, Math. Z. 1991):
kA_n natural order T = D(A) (injective cogenerator); opposite order T = A (regular). Self-cert:
Ext^1(T,∇)=0 AND is_tilting_module(T). Oracle: double Ringel dual has the same Cartan Smith
normal form as A. Over QQ (decompose/presented_form char-clean at char 0)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.ext import ext_dims
from quiverlab.modules.hom import is_isomorphic
from quiverlab.modules.morphism import direct_sum
from quiverlab.modules.quasihereditary import (characteristic_tilting, costandard_modules,
                                               ringel_dual)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _a3():
    return linear_path_algebra(3, field=QQ)


@lit
def test_natural_order_T_is_D_of_A():
    A = _a3()
    T = characteristic_tilting(A)                       # natural order
    DA, _, _ = direct_sum(*[A.injective(v) for v in (1, 2, 3)])   # D(A) = ⊕ I(v)
    assert is_isomorphic(T, DA)


@lit
def test_opposite_order_T_is_regular():
    A = _a3()
    T = characteristic_tilting(A, [3, 2, 1])
    reg, _, _ = direct_sum(*[A.projective(v) for v in (1, 2, 3)])  # A_A
    assert is_isomorphic(T, reg)


@selfcert
def test_T_is_tilting_and_ext_perp_nabla():
    from quiverlab.modules.tilting import is_tilting_module     # P44
    A = _a3()
    T = characteristic_tilting(A)
    assert is_tilting_module(T).is_tilting is True
    N = costandard_modules(A)
    for j in (1, 2, 3):
        assert ext_dims(A, T, N[j], 1)[1] == 0                  # Ext^1(T, ∇(j)) = 0


@lit
def test_double_ringel_dual_cartan_smith_form():
    import sympy
    A = _a3()
    R2 = ringel_dual(ringel_dual(A))
    def invf(M):
        S = sympy.Matrix([[int(x) for x in row] for row in M.cartan_matrix()])
        return list(sympy.Matrix(S).invariant_factors())
    assert invf(R2) == invf(A)                                  # Ringel duality involution


@selfcert
def test_ringel_dual_of_ka3_is_a_path_algebra():
    A = _a3()
    R = ringel_dual(A)
    assert R.quiver is not None                                 # presented (char 0 => kQ/I)
    assert R.dim == A.dim                                       # Morita to a kA3-shaped algebra
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```python
def characteristic_tilting(A, order=None):
    ranks, order = _order_ranks(A, order)
    deltas = standard_modules(A, order)
    nablas = costandard_modules(A, order)
    from quiverlab.modules.ext import ext_dims
    from quiverlab.modules.yoneda import baer_extension
    T = {}                                          # vertex -> T(i)
    for i in order:                                 # low -> high
        cur = deltas[i]
        for j in order:
            if ranks[j] >= ranks[i]:
                break                               # only lower T(j) absorb into T(i)
            d = ext_dims(A, deltas[i], T[j], 1)[1]  # dim Ext^1(Δ(i), T(j))
            for _ in range(d):                      # universal extension, one class at a time
                _dims, payload = ext_dims(A, cur, T[j], 1, with_reps=True)
                if _dims[1] == 0:
                    break
                cocycle = _one_ext_cocycle(A, cur, T[j], payload)   # P_1(cur) -> T[j]
                cur = baer_extension(cur, T[j], cocycle).middle
        T[i] = cur
    Tmod, _, _ = direct_sum(*[T[i] for i in order])
    # SELF-CERTIFY (the arbiter of the assembly, P37/P44 precedent):
    from quiverlab.modules.tilting import is_tilting_module
    for j in order:
        if ext_dims(A, Tmod, nablas[j], 1)[1] != 0:
            raise QuiverlabError(
                f"characteristic_tilting: Ext^1(T, ∇({j})) != 0 -- the universal extension "
                f"assembly is wrong for order {order!r}", hint="report this presentation")
    if not is_tilting_module(Tmod).is_tilting:
        raise QuiverlabError(
            "characteristic_tilting: T is not tilting (self-certificate failed)",
            hint="the Δ set or the extension order is inconsistent")
    return Tmod


def ringel_dual(A, order=None):
    T = characteristic_tilting(A, order)
    from quiverlab.modules.endomorphism import end_algebra
    R = end_algebra(T).opposite()                   # End_A(T)^op
    try:
        from quiverlab.core.basic import presented_form
        return presented_form(R)                    # kQ/I (char 0 / char > dim, split blocks)
    except QuiverlabError:
        R._ringel_note = ("Ringel dual returned in structure-constant form: presented_form "
                          "refused (char <= dim or a non-split block)")
        return R
```

**Adjust to reality (Task D):**
- `_one_ext_cocycle(A, M, N, payload)` reconstructs ONE `Ext¹(M, N)` class as a `N.dim ×
  P_1(M).dim` cocycle that `baer_extension` accepts — reuse the tested reconstruction from
  P41's `ar.py::_ext1_data` / `complex_reps._reconstruct_cocycle` (or P44's Bongartz
  `_universal_cocycle`), do NOT re-derive representatives. The iterated single-class extension
  (rather than a one-shot block cocycle) is the honest fallback; the SELF-CERTIFICATE
  (`Ext¹(T,∇)=0` AND `is_tilting_module(T)`) is the arbiter of whether the assembly is right —
  a wrong order/assembly fails the cert and raises loudly, exactly as P44's Bongartz completion
  lets its `is_tilting_module(T⊕E)` cert arbitrate.
- `characteristic_tilting` runs `is_tilting_module` (P30 `decompose` count) and `ringel_dual`
  runs `presented_form` (P44) — both inherit the char 0 / char > dim scope; the batteries use
  QQ. Off scope they raise the SAME loud `QuiverlabError` P30/P44 raise (never a silent wrong
  T/dual). Δ/∇ inside are char-clean, so the refusal is localized and honestly worded.
- The double-dual Cartan Smith-form oracle needs `ringel_dual(A)` to be PRESENTED (so
  `ringel_dual(ringel_dual(A))` can build projectives) — over QQ `presented_form` succeeds and
  yields `kQ/I`. If `presented_form` degraded, `ringel_dual(ringel_dual(...))` would refuse at
  the second `characteristic_tilting`; the QQ battery avoids that. State this dependency in the
  test docstring.

- [ ] **Step 4: Run tests** — Expected: PASS (expect a couple of iterations on the extension
  order — the self-cert is the arbiter)
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/quasihereditary.py tests/modules/test_quasihereditary_ringel.py
git commit -m "feat(modules): characteristic_tilting (Ringel universal extension, self-certified) + ringel_dual (End(T)^op, presented)"
```

---

### Task E: recollements from an idempotent — `Recollement(A, S)`

**Files:**
- Create: `src/quiverlab/modules/recollement.py`
- Test: `tests/modules/test_recollement.py`

**Interfaces:**
- Consumes: `invariants/pathbasis.py::path_type_basis(A) -> (idem, rad, src, tgt)`
  (`pathbasis.py:17` — idempotent/radical basis-index lists + per radical index its source and
  target vertex), `Algebra.multiply`/`.T`/`.unit`/`.dim`/`.domain`/`.basis_labels`/`._basis_vec`/
  `.from_structure_constants` (`core/algebra.py`), `Quiver.algebra` +
  `families/trivial_extension.py` length-lex relation extraction (`_relation_string`,
  `extract_relations`, `_solve_combo` — reuse for the `A/AeA` presentation), `linalg_mod`
  (`column_space_pivots`/`cols_to_matrix`/`kernel_columns`/`mat_rank`/`matvec`),
  `fields.linalg.solve`, `modules/tor.py::_induced`/`_vertex_basis`/`_left_action` (the
  induction/tensor pattern for `j_!`), `morphism.py::hom_basis`/`ModuleHom`,
  `hom.py::hom_space`/`is_isomorphic`, `Module`/`Module.with_side`/`.vertex_projection`.
- Produces:
  ```python
  class Recollement:
      def __init__(self, A, S)          # S: a vertex subset; builds eAe (certified) + A/AeA
      eAe        # Algebra: corner e_S A e_S on the corner path-type basis (structure consts)
      quotient   # Algebra: A/A e_S A, presented on the complement subquiver
      # the six functors (ASS/CPS/BBD recollement of (mod A/AeA, mod A, mod eAe)):
      def i_star(self, N)          # inflation: an A/AeA-module as an A-module
      def i_upper_star(self, M)    # M ↦ M/M·(AeA)                 (left adjoint of i_*)
      def i_upper_shriek(self, M)  # M ↦ the AeA-annihilated submodule (right adjoint of i_*)
      def j_upper_star(self, M)    # M ↦ Me with the eAe-action    (the "corner restriction")
      def j_shriek(self, X)        # X ↦ X ⊗_{eAe} eA              (left adjoint of j^*)
      def j_star(self, X)          # X ↦ Hom_{eAe}(Ae, X)          (right adjoint of j^*)
  ```

**Math pinned (the `eAe`-vs-subquiver trap + the worked `kA₃`, `S={2}` example — verbatim in
the docstring).** `kA₃` (1→2→3, basis `e₁,e₂,e₃,a,b,ab`, dim 6):
- **`S={2}`, the canonical recollement.** `A e₂ A` = span of every basis path visiting vertex
  2 = `{e₂, a, b, ab}` (dim 4), so `A/Ae₂A = k×k` on vertices `{1,3}` with no arrows (dim 2).
  `eAe = e₂Ae₂ = k` (dim 1, since `e₂ a e₂ = e₂ b e₂ = e₂ (ab) e₂ = 0`). Worked functors:
  `j^*(P(1)) = P(1)e₂ = [1..3]e₂` = the 1-dim vertex-2 part (a `k`-space);
  `i^*(P(1)) = [1..3]/[1..3]·Ae₂A = [1..3]/[2..3] = S₁` as an `A/Ae₂A`-module.
- **The trap, `S={1,3}`.** The full subquiver on `{1,3}` has NO arrows (both `a,b` touch 2), so
  the subquiver algebra is `k×k` (dim 2). But `e(ab)e = ab ≠ 0` because the path `ab:1→3`
  leaves `S` only in its INTERIOR (vertex 2), so **`eAe = span{e₁, e₃, ab} = kA₂` (dim 3)** —
  NOT the subquiver algebra. This is why `eAe` is built on the corner path-type basis via
  structure constants, certified `dim eAe = Σ_{v,w∈S} dim e_v A e_w = 1+1+1 = 3`, while `A/AeA`
  IS the complement subquiver, certified `dim = dim A − dim A e_S A`.

**Self-certs (the oracle class — no QPA):** the four adjunction dim identities, the two
canonical exact sequences at each joint (rank identities), and `j^*j_!≅id`/`j^*j_*≅id`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_recollement.py
"""Recollement from an idempotent (Plan 47). The eAe-vs-subquiver trap is the pinned math;
the six functors' adjunction dim identities + the canonical exact sequences + j^*j_! ≅ id
are the self-cert oracle (QPA white space). Char-clean (a GF(2) cell). Over kA3, the
commutative square, and a multi-vertex zoo algebra."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.hom import hom_space, is_isomorphic
from quiverlab.modules.recollement import Recollement

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _a3(field=QQ):
    return linear_path_algebra(3, field=field)


def _square(field=QQ):
    # commutative square 1->2->4, 1->3->4 with ab = cd (one relation).
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b-c*d"], field=field)


@lit
def test_eae_is_not_the_subquiver_algebra():
    # S = {1,3}: the path ab:1->3 travels through 2 ∉ S but survives in eAe.
    A = _a3()
    R = Recollement(A, [1, 3])
    assert R.eAe.dim == 3                          # kA2, NOT the subquiver k×k (dim 2)
    assert len(list(R.eAe.quiver.vertices)) == 2 if R.eAe.quiver else True


@lit
def test_worked_ka3_S2():
    A = _a3()
    R = Recollement(A, [2])
    assert R.eAe.dim == 1                          # e2 A e2 = k
    assert R.quotient.dim == 2                     # A/Ae2A = k×k on {1,3}
    P1 = A.projective(1)
    assert R.j_upper_star(P1).dim == 1             # P(1) e2 = vertex-2 part (1-dim)
    assert is_isomorphic(R.i_upper_star(P1), R.i_star(R.quotient.simple(1)))  # = S1 inflated


@selfcert
def test_dim_certificates():
    A = _a3()
    for S in ([1], [2], [3], [1, 3], [1, 2]):
        R = Recollement(A, S)
        # eAe corner-dim certificate + A/AeA complement-dim certificate are asserted inside
        # __init__; re-check the public identity here.
        assert R.eAe.dim == R._corner_dim_sum()
        assert R.quotient.dim == A.dim - R._aeA_dim()


@selfcert
@pytest.mark.parametrize("build", [_a3, _square])
def test_j_adjunctions_and_counit_iso(build):
    A = build()
    R = Recollement(A, [2] if A.dim == 6 else [2, 3])
    dom = A.domain
    # (j_!, j^*):  dim Hom_A(j_! X, M) == dim Hom_eAe(X, j^* M)
    # (j^*, j_*):  dim Hom_A(M, j_* X) == dim Hom_eAe(j^* M, X)
    X = R.eAe.simple(next(iter(R.eAe.quiver.vertices))) if R.eAe.quiver else None
    for M in [A.projective(v) for v in A.quiver.vertices]:
        if X is not None:
            assert len(hom_space(R.j_shriek(X), M)) == len(hom_space(X, R.j_upper_star(M)))
            assert len(hom_space(M, R.j_star(X))) == len(hom_space(R.j_upper_star(M), X))
    # counit iso: j^* j_! X ≅ X and j^* j_* X ≅ X
    if X is not None:
        assert is_isomorphic(R.j_upper_star(R.j_shriek(X)), X)
        assert is_isomorphic(R.j_upper_star(R.j_star(X)), X)


@selfcert
def test_i_adjunctions():
    A = _a3()
    R = Recollement(A, [2])
    B = R.quotient
    for N in [B.simple(v) for v in B.quiver.vertices]:
        for M in [A.projective(v) for v in A.quiver.vertices]:
            # (i^*, i_*):  dim Hom_A(M, i_* N) == dim Hom_B(i^* N?, ...) -- pin the pair
            assert len(hom_space(M, R.i_star(N))) == len(hom_space(R.i_upper_star(M), N))
            assert len(hom_space(R.i_star(N), M)) == len(hom_space(N, R.i_upper_shriek(M)))


@selfcert
def test_canonical_sequence_exact_at_each_joint():
    # j_! j^* M --> M --> i_* i^* M --> 0 exact: second is epi, im(first)=ker(second).
    A = _a3()
    R = Recollement(A, [2])
    for M in [A.projective(v) for v in (1, 2, 3)]:
        assert R._counit_sequence_exact(M)         # rank identities inside; returns bool


@selfcert
def test_recollement_char_clean_gf2():
    A = _a3(GF(2))
    R = Recollement(A, [1, 3])
    assert R.eAe.dim == 3                           # no char caveat on the corner build
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.modules.recollement`

- [ ] **Step 3a: Implement the two algebras `eAe` and `A/AeA`**

```python
"""Recollements from an idempotent (Plan 47). For a vertex subset S with e = Σ_{v∈S} e_v:
the corner algebra eAe (NOT the subquiver algebra -- a path may leave S in its interior yet
survive e·p·e; built on the corner path-type basis via structure constants, certified
dim = Σ_{v,w∈S} dim e_v A e_w), the quotient A/AeA (presentable on the complement full
subquiver, certified dim = dim A − dim A e_S A), and the six functors of the recollement
(mod A/AeA, mod A, mod eAe) as bimodule operations. Char-clean exact linear algebra.
Cline-Parshall-Scott (J. Reine Angew. Math. 1988); recollement: Beilinson-Bernstein-Deligne
(Asterisque 100, 1982)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm


class Recollement:
    def __init__(self, A, S):
        from quiverlab.modules.builders import _require_provenance
        _require_provenance(A, "Recollement")
        self.A = A
        self.S = list(S)
        self.Sset = set(self.S)
        self._e = self._idempotent_vec()            # e_S = Σ_{v∈S} e_v, in A-coords
        self._aeA_cols = self._aeA_span()            # basis of A e_S A (column span)
        self.eAe = self._build_corner()              # certified inside
        self.quotient = self._build_quotient()       # certified inside

    def _idempotent_vec(self):
        A, dom = self.A, self.A.domain
        e = [dom.zero()] * A.dim
        for i, lab in enumerate(A.basis_labels):
            if lab.startswith("e_") and any(lab == f"e_{v}" for v in self.S):
                e[i] = dom.one()
        return e

    def _aeA_span(self):
        """A e_S A = span{ b_k · e_S · b_l : all k, l }: every element visiting some v∈S."""
        A, dom = self.A, self.A.domain
        cols = []
        for k in range(A.dim):
            bk_e = A.multiply(A._basis_vec(k), self._e)
            for l in range(A.dim):
                cols.append(A.multiply(bk_e, A._basis_vec(l)))
        piv = lm.column_space_pivots(lm.cols_to_matrix(cols), dom)
        return [cols[j] for j in piv]

    def _aeA_dim(self):
        return len(self._aeA_cols)

    def _corner_dim_sum(self):
        """Σ_{v,w∈S} dim e_v A e_w from the path-type basis (the corner certificate target)."""
        from quiverlab.invariants.pathbasis import path_type_basis
        idem, rad, src, tgt = path_type_basis(self.A, "Recollement corner")
        idlabel = {i: v for i, v in [(i, self._vertex_of(i)) for i in idem]}
        n = 0
        for i in idem:                               # idempotents e_v with v∈S
            if self._vertex_of(i) in self.Sset:
                n += 1
        for i in rad:                                # radical paths f_i with src,tgt ∈ S
            if self._vertex_of(src[i]) in self.Sset and self._vertex_of(tgt[i]) in self.Sset:
                n += 1
        return n

    def _build_corner(self):
        """eAe on the basis = column_space_pivots of { e·b_k·e : k } (the corner elements),
        structure constants = A-products re-expressed in that basis. Certified dim ==
        _corner_dim_sum(). Carries e_v labels for v∈S so it stays path-type when possible."""
        A, dom = self.A, self.A.domain
        raw = [A.multiply(A.multiply(self._e, A._basis_vec(k)), self._e) for k in range(A.dim)]
        piv = lm.column_space_pivots(lm.cols_to_matrix(raw), dom)
        basis = [raw[j] for j in piv]
        P = lm.cols_to_matrix(basis)
        m = len(basis)
        def coords(vec):
            return linalg.solve(P, list(vec), dom)
        T = [[[coords(A.multiply(basis[a], basis[b]))[c] for c in range(m)]
              for b in range(m)] for a in range(m)]
        unit = coords(self._e)
        eAe = A.from_structure_constants(
            [[[_fmt(x, dom) for x in cell] for cell in rowa] for rowa in T],
            [_fmt(x, dom) for x in unit], field=dom, check=True)   # check=True self-certifies
        if eAe.dim != self._corner_dim_sum():
            raise QuiverlabError(
                f"Recollement: eAe dim {eAe.dim} != corner sum {self._corner_dim_sum()}",
                hint="the corner path-type basis is inconsistent -- report this presentation")
        return eAe

    def _build_quotient(self):
        """A/A e_S A: presented on the complement full subquiver (vertices V∖S, arrows between
        them), certified dim == dim A − dim A e_S A. Reuses the trivial-extension length-lex
        relation extraction; falls back to a structure-constant quotient build when the basis
        is not string-representable (the TrivialExtension D3 fallback)."""
        A = self.A
        target = A.dim - self._aeA_dim()
        ...  # complement subquiver Q' = (V∖S, arrows with both ends in V∖S); pi: kQ' -> A/AeA
             # via the quotient projection; extract_relations(max_len) length-lex; certify
             # dim(kQ'/I') == target with the widen-once window; loud QuiverlabError otherwise.
             # Structure-constant fallback: quotient A by the two-sided ideal _aeA_cols directly.
```

- [ ] **Step 3b: Implement the six functors**

```python
    # --- j-side (corner eAe) -------------------------------------------------
    def j_upper_star(self, M):
        """M ↦ Me_S with the eAe-action: the underlying space is im(Σ_{v∈S} M.action[e_v]);
        an eAe basis element (= e·b·e) acts by M.action[b] restricted to Me_S."""
        ...  # build a Module over self.eAe on the corner subspace of M; e_v-projections give
             # the vertex grading, the corner path elements give the arrow action.

    def j_shriek(self, X):
        """X ↦ X ⊗_{eAe} eA, a right A-module (induction). eA is an (eAe, A)-bimodule; reuse
        the modules.tor induction pattern (_vertex_basis / _left_action / _induced) applied to
        the bimodule eA over eAe. Self-cert: j^* j_! X ≅ X (counit iso)."""
        ...

    def j_star(self, X):
        """X ↦ Hom_{eAe}(Ae, X), a right A-module (coinduction, right adjoint of j^*). Ae is
        an (A, eAe)-bimodule; build the Hom-space over eAe with the residual A-action.
        Self-cert: j^* j_* X ≅ X."""
        ...

    # --- i-side (quotient A/AeA) ---------------------------------------------
    def i_star(self, N):
        """An A/AeA-module N as an A-module via A ->> A/AeA (inflation, exact fully faithful):
        A-action label b acts as its image in A/AeA; any path visiting S acts as 0."""
        ...

    def i_upper_star(self, M):
        """M ↦ M/M·(AeA) = M ⊗_A A/AeA (left adjoint of i_*): quotient M by the submodule
        M·(A e_S A) = span{ M.action[g] applied : g in _aeA_cols }; view as an A/AeA-module."""
        ...

    def i_upper_shriek(self, M):
        """M ↦ Hom_A(A/AeA, M) = the largest submodule of M annihilated by AeA (right adjoint
        of i_*): the joint kernel ∩_{g ∈ AeA} ker(M.action[g]); view as an A/AeA-module."""
        ...

    # --- the certificate helpers (used by the tests) -------------------------
    def _counit_sequence_exact(self, M):
        """j_! j^* M --> M --> i_* i^* M --> 0 exact at each joint (rank identities)."""
        ...
```

**Adjust to reality (Task E):**
- `_vertex_of(i)` maps a `path_type_basis` idempotent-index to its vertex label by matching
  `A.basis_labels[i] == f"e_{v}"`; factor it once. The corner certificate `_corner_dim_sum`
  reads the src/tgt of each radical basis element — a path with `src∈S` and `tgt∈S` lands in
  `eAe` REGARDLESS of its interior (this IS the trap), so counting by endpoints is correct for
  the corner dim while the `A/AeA` dim uses the visiting-span `_aeA_dim` (interior-aware).
- `from_structure_constants(check=True)` validating associativity + unit IS the self-cert that
  `eAe` is a well-defined algebra; `_fmt(x, dom)` renders a Domain element into the token form
  `from_structure_constants` re-parses (mirror `families/trivial_extension.py` /
  `core/basic.py::_semisimple_quotient`'s `_fmt` — over QQ a `Fraction`, over GF(p) the int; or
  extend `from_structure_constants` to accept coerced Domain elements — pick the path an
  existing caller uses, read `core/algebra.py:46`).
- `j_shriek` (the tensor `− ⊗_{eAe} eA`) is the one heavy functor. The `modules/tor.py`
  `_induced`/`_vertex_basis`/`_left_action` idiom (`tor.py:67,79,95`) is the template applied
  to the bimodule `eA` over `eAe`; the SELF-CERT `j^* j_! X ≅ X` (`test_j_adjunctions_and_
  counit_iso`) is the arbiter — if the tensor grading/action is wrong the counit iso fails,
  fix and re-run (never weaken to a dim-only check). Likewise `j_star`'s `Hom_{eAe}(Ae, −)` is
  arbitrated by `j^* j_* X ≅ X`.
- The four adjunction identities and the two canonical exact sequences are the self-cert
  BATTERY (`oracle_selfcert`); they hold simultaneously only when all six functors are correct,
  exactly the "multiple identities pin one convention" discipline of P35's product surface.
  Include the multi-vertex zoo algebra (`_square`, and one `line_abc_cde`-style record) so the
  identities decide on a non-directed example, not only `kA₃`.
- `i_star` on an `A/AeA`-module needs the label-to-label action map `A → A/AeA`; the quotient
  projection built in `_build_quotient` supplies it. Keep `i_*` exact and fully faithful (its
  image = the AeA-torsion-free... = modules supported off `S`); the `test_i_adjunctions`
  identities arbitrate the direction of each adjunction.

- [ ] **Step 4: Run tests** — Expected: PASS (expect iteration on the `j`-side functor
  gradings — the counit isos arbitrate)
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/recollement.py tests/modules/test_recollement.py
git commit -m "feat(modules): Recollement(A,S) -- eAe (corner, NOT subquiver) + A/AeA + six functors, adjunction/exactness self-certified"
```

---

### Task F: `Algebra` wrappers + the `quasi_hereditary` GUI scalar kind

**Files:**
- Modify: `src/quiverlab/core/algebra.py` (thin wrappers, lazy-import to dodge the cycle,
  beside `global_dimension`/`is_basic`)
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch` scalar-kind branch + `_snip`) and
  `docs/gui/runner.py` (the byte-identical twin)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox + `S.ids` + push-list +
  `renderBlock` branch), `webapp/static/app.js` (block builder)
- Modify: `webapp/server/i18n/en.json` + `es.json` (`inv.quasi_hereditary`,
  `block.quasi_hereditary.*`) + `webapp/templates/index.html` (`data-*`)
- Modify: `src/quiverlab/trace/results_html.py` (`_HEADINGS` + `_quasi_hereditary_html`)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE new fixture, documented)
- Test: `tests/modules/test_quasihereditary_wrappers.py`,
  `tests/webapp/test_quasi_hereditary_p47.py`,
  `tests/gui/test_quasi_hereditary_runner_twin.py`

**Interfaces:**
- `Algebra` wrappers (each a two-line lazy-import delegate, char scope inherited):
  ```python
  def standard_modules(self, order=None)         # -> dict vertex -> Δ(i)
  def is_quasi_hereditary(self, order=None)      # -> QHReport
  def characteristic_tilting(self, order=None)   # -> Module
  def ringel_dual(self, order=None)              # -> Algebra
  def corner_algebra(self, vertices)             # -> Recollement(self, vertices).eAe
  def quotient_by_idempotent(self, vertices)     # -> Recollement(self, vertices).quotient
  def recollement(self, vertices)                # -> Recollement(self, vertices)
  ```
- `quasi_hereditary` is an **algebra-scalar** compute kind (schema v1, NO module block — the
  `recognizers`/`ext_algebra` precedent, `spec.py:1212,1220`), NOT a module kind. The block
  builder lives once in `quasihereditary.py::quasi_hereditary_block(A, order=None)` and BOTH
  runners call it (byte-identical). It reports the NATURAL order with the honest order-dependence
  note. Block shape:
  ```python
  {"kind": "quasi_hereditary",
   "is_quasi_hereditary": bool, "order": [...], "order_note": "quasi-heredity is order-dependent;"
       " this report uses the natural vertex order",
   "gl_dim": {"value": int, "exact": bool},
   "per_index": {v: {"brick": bool, "delta_filters_P": bool, "note": str}},
   "standard_dims": {v: {"dim": int, "dimvec": {...}}},   # the Δ dim-vectors table
   "note": str, "references": [...], "citations": [...]}   # refs = ["dlab_ringel","assem_book"]
  ```

- [ ] **Step 1: Write the failing wrapper test**

```python
# tests/modules/test_quasihereditary_wrappers.py
"""Algebra-level wrappers for the quasi-hereditary surface (Plan 47)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.hom import is_isomorphic

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return linear_path_algebra(3, field=QQ)


def test_wrappers_delegate():
    A = _a3()
    assert A.is_quasi_hereditary().is_quasi_hereditary is True
    assert all(A.standard_modules()[v].dim == 1 for v in (1, 2, 3))
    R = A.recollement([2])
    assert A.corner_algebra([1, 3]).dim == 3               # the eAe trap via the wrapper
    assert A.quotient_by_idempotent([2]).dim == 2
    assert R.eAe.dim == 1
    T = A.characteristic_tilting()
    from quiverlab.modules.morphism import direct_sum
    DA, _, _ = direct_sum(*[A.injective(v) for v in (1, 2, 3)])
    assert is_isomorphic(T, DA)
```

- [ ] **Step 2: Write the failing cross-runner test** (unmarked — extras-gated dir; copy the
  `tests/webapp/test_module_blocks_m0729.py` runner-pair fixture):

```python
# tests/webapp/test_quasi_hereditary_p47.py
"""The quasi_hereditary algebra-scalar kind: schema-1, served by hpc.spec, mirrored by the
Pyodide twin, both byte-identical on the block."""


def test_quasi_hereditary_block_shape(tmp_path):
    from quiverlab.hpc.spec import ComputeRequest, run_spec
    req = _qh_request()                              # kA3 over QQ, compute ["quasi_hereditary"]
    out = run_spec(ComputeRequest.model_validate(req), tmp_path)
    b = out["results"]["quasi_hereditary"]
    assert b["is_quasi_hereditary"] is True
    assert b["gl_dim"] == {"value": 1, "exact": True}
    assert set(b["standard_dims"]) == {1, 2, 3} or set(b["standard_dims"]) == {"1", "2", "3"}
    assert "dlab_ringel" in [k for k, _ in b["citations"]]
    assert "order-dependent" in b["order_note"]


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py; assert json.dumps(sort_keys=True)
    # equality on the quasi_hereditary block (both runners byte-identical).
    ...
```

- [ ] **Step 3: Implement**
- Add the seven `Algebra` wrappers (lazy `from quiverlab.modules.quasihereditary import ...`
  / `from quiverlab.modules.recollement import Recollement` inside each, to avoid the
  `core → modules → core` import cycle — the `global_dimension`/`is_tilting_module` wrapper
  precedent, `core/algebra.py:448`).
- Add `quasi_hereditary_block(A, order=None)` to `quasihereditary.py` (build the block dict
  above from `is_quasi_hereditary` + `standard_modules`; `references = ["dlab_ringel",
  "assem_book"]`).
- `spec.py::_dispatch`: after the `recognizers` branch (`spec.py:1220`), add
  ```python
  if kind == "quasi_hereditary":
      from quiverlab.modules.quasihereditary import quasi_hereditary_block
      block = quasi_hereditary_block(A)
      block["citations"] = _citation_pairs(block["references"])
      return block, None
  ```
  and the `_snip` recipe (`"quasi_hereditary": "A.is_quasi_hereditary()"`, `spec.py:1755`).
  Mirror byte-for-byte in `docs/gui/runner.py` after ITS `recognizers` branch
  (`runner.py:695`).
- `results_html.py`: `_HEADINGS["quasi_hereditary"] = "Quasi-hereditary structure"`
  (`results_html.py:29`) + a `_quasi_hereditary_html(b)` branch (`results_html.py:337` ladder)
  rendering the flags, per-index certificates, the Δ dim-vectors table, and the
  order-dependence note.
- GUI (7 touchpoints, the `recognizers` scalar-kind precedent, `gui.js:78,150,663`): checkbox
  `qlgui-quasi_hereditary`, the `S.ids` list, the plain-kind push-list (`recognizers`,
  `homological_profile` sit there), the `renderBlock` branch (a flags + Δ-dims table via the
  existing helpers) in BOTH `gui.js` files, the `app.js` builder, and the i18n
  `inv.quasi_hereditary` / `block.quasi_hereditary.*` chain (EN + ES + `index.html` `data-*`).

- [ ] **Step 4: Add ONE golden fixture** `algebra_quasi_hereditary_kA3` to
  `_runner_goldens.json`; note it in `test_runner_delegation.py`'s docstring change-log (the
  `products_loop_gf2` precedent — new fixture, existing entries byte-identical). Run the
  delegation test BEFORE adding to confirm existing entries are untouched.

- [ ] **Step 5: Run the gates**

Run: `... -m pytest tests/modules/test_quasihereditary_wrappers.py tests/webapp/test_quasi_hereditary_p47.py tests/webapp/test_runner_delegation.py tests/gui/test_quasi_hereditary_runner_twin.py tests/hpc -q`
Expected: PASS; both runners byte-identical.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(core,gui,webapp,hpc): Algebra qh/recollement wrappers + quasi_hereditary scalar compute kind (natural order, order-dependence note)"
```

---

### Task G: verification page, citations, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-05-metaplan-v0.2.0.md` (tick the P47 card delivery note)
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`)

- [ ] **Step 1: Citations** (VERIFIED BibTeX only — the `_r(...)` registry precedent,
  `registry.py:24,128`). Add:
  - `dlab_ringel` → `DlabRingel1989` (Dlab, Ringel, *Quasi-hereditary algebras*, Illinois
    J. Math. **33** (1989), no. 2, 280–291).
  - `cps` → `CPS1988` (Cline, Parshall, Scott, *Finite-dimensional algebras and highest weight
    categories*, J. Reine Angew. Math. **391** (1988), 85–99).
  - `ringel_dual` → `Ringel1991` (Ringel, *The category of modules with good filtrations over
    a quasi-hereditary algebra has almost split sequences*, Math. Z. **208** (1991), 209–223).
  - `bbd` → `BBD1982` (Beĭlinson, Bernstein, Deligne, *Faisceaux pervers*, Astérisque **100**,
    Soc. Math. France, 1982) — the recollement origin.

  ```bibtex
  @article{DlabRingel1989,
    author  = {Dlab, Vlastimil and Ringel, Claus Michael},
    title   = {Quasi-hereditary algebras},
    journal = {Illinois Journal of Mathematics},
    volume  = {33}, number = {2}, pages = {280--291}, year = {1989},
  }
  ```
  and (mirroring `_r("assem_book", "ASS2006", "foundation", ...)`):
  ```python
  _r("dlab_ringel", "DlabRingel1989", "foundation",
     "Quasi-hereditary algebras",
     "The definition of quasi-hereditary algebras, standard modules, and the "
     "quasi-heredity test used in Plan 47; qh => finite gl.dim.", "article"),
  _r("ringel_dual", "Ringel1991", "foundation",
     "Good filtrations and the characteristic tilting module",
     "The characteristic tilting module and the Ringel dual (Plan 47).", "article"),
  _r("cps", "CPS1988", "foundation",
     "Finite-dimensional algebras and highest weight categories",
     "Highest weight categories and the idempotent recollement (eAe, A/AeA) of Plan 47.",
     "article"),
  _r("bbd", "BBD1982", "foundation",
     "Faisceaux pervers",
     "The origin of recollement and the six-functor formalism (Plan 47).", "book"),
  ```
  Reuse `assem_book` (ASS2006) where it covers the general module theory. Each BibTeX entry is
  VERIFIED before it ships (`tests/citations/` resolves the keys).

- [ ] **Step 2: Verification page.** Add the Plan-47 subsystem rows:
  `modules/quasihereditary.py` — `oracle_literature`: `kA_n` natural-order Δ(i)=S_i /
  ∇(i)=[1..i], the two-orders quasi-heredity (Dlab–Ringel), the `k[x]/(x²)` negative, the
  hand-derived `T = D(A)` (natural) / `T = A` (opposite), the double-Ringel-dual Cartan
  Smith-form identity; `oracle_crossengine`: BGG reciprocity `(P(i):Δ(j)) = [∇(j):S(i)]`;
  `oracle_selfcert`: top Δ(i)=S_i / [Δ(i):S(i)]=1, `Ext¹(T,∇)=0` + `is_tilting_module(T)`,
  the delta-peel certificate. `modules/recollement.py` — `oracle_literature`: the
  `eAe`-vs-subquiver trap dims (`S={1,3}` → `kA₂`, not `k×k`) + the worked `kA₃` `S={2}`;
  `oracle_selfcert`: the four adjunction dim identities, the two canonical exact sequences at
  each joint, `j^*j_!≅id`/`j^*j_*≅id`. Add the **honest-scope entries**: (a) **QPA white
  space** — QPA has no quasi-hereditary / recollement surface, so the oracle class is theory
  pins + internal self-certificates, no QPA agreement; (b) quasi-heredity is
  **order-dependent** and the GUI reports the natural order only (stated in the block);
  (c) the tilting summand-count self-cert (P30) and PRESENTED Ringel duals (P44) inherit the
  char 0 / char > dim caveat — Δ/∇/qh/recollement are char-clean (a `GF(2)` cell proves it);
  (d) `stratifications` beyond the quasi-hereditary case are a named successor (not shipped).
  Recount the class table (`tests/release/test_oracle_classes.py` drives the numbers — run
  collection, paste the LIVE counts, re-run to green; NEVER a guessed-at-authoring number
  given the mid-merge-train drift).

- [ ] **Step 3: README.** One features line: "standard/costandard modules, a quasi-heredity
  test (Dlab–Ringel), good-filtration multiplicities + BGG reciprocity, the characteristic
  tilting module and its Ringel dual, and recollements from an idempotent (eAe, A/AeA, the six
  functors) — quasi-hereditary algebras, white space in QPA."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/modules -q` (deep, the touched files),
  `... -m pytest -q -m fast`,
  `... -m pytest tests/release tests/citations -q` — all green
  (no `tests/qpa/` this plan — QPA white space).

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-47 quasi-hereditary + recollement oracle rows + honest scope (QPA white space, order-dependence, eAe trap, char-clean core) + recounted classes + Dlab-Ringel/Ringel/CPS/BBD citations"
```

---

## Acceptance (Plan-47 definition of done)

1. `trace_module`, `standard_modules`/`costandard_modules`, `is_quasi_hereditary`
   (`QHReport`) + `delta_multiplicities`, `characteristic_tilting` + `ringel_dual`, and
   `Recollement(A, S)` (with `eAe`, `A/AeA`, and the six functors) all public, plus the seven
   `Algebra` wrappers (`standard_modules`, `is_quasi_hereditary`, `characteristic_tilting`,
   `ringel_dual`, `corner_algebra`, `quotient_by_idempotent`, `recollement`). Every
   construction CERTIFIED per instance (a dim identity + a structural oracle) or refusing
   loudly (three-valued `QHReport`/`delta_multiplicities`; loud `characteristic_tilting`/
   `ringel_dual` refusals).
2. The `eAe`-vs-subquiver trap is pinned: `Recollement(kA₃, {1,3}).eAe` is `kA₂` (dim 3), NOT
   the subquiver `k×k` (dim 2), certified `dim = Σ_{v,w∈S} dim e_v A e_w`; `A/AeA` is the
   complement subquiver, certified `dim = dim A − dim A e_S A`; the worked `kA₃`, `S={2}`
   recollement is derived by hand in the plan and pinned in a test.
3. `kA_n` is quasi-hereditary for TWO orders (Dlab–Ringel), with Δ(i)=S_i (natural) /
   Δ(i)=P(i) (opposite) hand-derived; `k[x]/(x²)` is NOT quasi-hereditary with a loud note;
   BGG reciprocity `(P(i):Δ(j)) = [∇(j):S(i)]` holds across both orders.
4. The characteristic tilting module is `T = D(A)` (natural) / `T = A` (opposite) as
   hand-derived, SELF-certified by `Ext¹(T,∇(j))=0` AND `is_tilting_module(T)` (P44); the
   double Ringel dual has the same Cartan invariant factors (Smith normal form) as `A`.
5. The six recollement functors satisfy their four adjunction dim identities, the two
   canonical exact sequences at each joint, and `j^*j_!≅id`/`j^*j_*≅id` — simultaneously,
   across `kA₃`, the commutative square, and a multi-vertex zoo algebra; these self-certs ARE
   the oracle class (QPA white space, stated).
6. Δ/∇/qh/recollement are char-clean (a `GF(2)` cell in each battery); only the tilting
   summand count (P30) and PRESENTED Ringel duals (P44) inherit the char 0 / char > dim
   caveat and refuse loudly off scope — never a silent wrong answer.
7. `quasi_hereditary` clickable end-to-end (GUI canvas → block → report) in EN+ES, schema v1
   (algebra-scalar), both runners byte-identical, ONE golden added with a documented
   change-log entry, reporting the natural order with the order-dependence note; the seven
   `Algebra` wrappers all reachable.
8. `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the honest-scope
   entries (QPA white space; order-dependence; the char-clean core vs the P30/P44-inherited
   tilting/Ringel caveat; stratifications deferred); README line added; Dlab–Ringel / Ringel /
   CPS / BBD citations BibTeX-verified; deep (touched dirs) + fast + release + citations suites
   green. No dependency taken on P39/P40/P41/P42/P43 (this branch merges independently to `dev`
   after P37 + P44).

---

## Methodology & assumptions

- **Structure copied verbatim from the two templates.** This plan mirrors
  `2026-08-05-plan-44-constructions.md` and `2026-08-05-plan-41-ar-completion.md` exactly:
  Goal / Architecture / Tech Stack / Global Constraints, then per-task **Files → Interfaces
  (Consumes/Produces with real signatures) → Step 1 failing tests (real pytest) → Step 2 run
  to fail → Step 3 implement (real code) → "Adjust to reality" → Step 4 run → Step 5 commit**,
  closing with an Acceptance checklist. The "oracle decides the convention" arbiter discipline
  (P37/P44/P35 precedent) is used for every ambiguous assembly: `dualize` side tag (∇), the
  greedy Δ-peel order, the Ringel universal-extension order, and the six-functor gradings are
  each pinned by a self-certificate the test asserts, never by an unchecked assumption.
- **Codebase grounding is verified, not assumed.** Every interface cited was read in the live
  tree: `ModuleHom.image() -> (I, epi, mono)` (`morphism.py:115`), `radtopsoc.submodule`/
  `quotient` (`:27/:47`), `duality.dualize` + `Algebra.opposite` (`duality.py:31`,
  `algebra.py:398`), `ext.global_dimension` + `GlobalDimension.exact` (`ext.py:85`),
  `path_type_basis -> (idem, rad, src, tgt)` (`pathbasis.py:17`), `builders.projective` =
  `e_v A` = paths STARTING at `v` (`builders.py:45`), `linear_path_algebra` = `i→i+1`
  (`families/basic.py:15`), the `recognizers`/`ext_algebra` algebra-scalar GUI kind
  (`spec.py:1212,1220`, `runner.py:695`, `results_html.py:29,337`, `gui.js:78,150,663`). P44's
  `is_tilting_module`/`presented_form` are consumed from the frozen Interfaces in the P44 plan
  doc (that plan merges first).
- **The two flagged ambiguities are resolved in the plan text, not deferred.**
  (1) *The hand-derived Δ/T for `kA_n`* — with the codebase's forward orientation (`i→i+1`),
  RIGHT modules, `P(v)=[v..n]`, `Hom(P(j),P(i))≠0 ⟺ i≤j` with image `[j..n]`: the **natural
  order** gives `Δ(i)=S_i`, `∇(i)=I(i)=[1..i]`, and `T=⊕I(v)=D(A)`; the **opposite order**
  gives `Δ(i)=P(i)`, `∇(i)=S_i`, and `T=A`. Both are quasi-hereditary (Dlab–Ringel two-orders
  oracle). The GUI reports the natural order. (2) *The `eAe`-vs-subquiver trap* — `eAe` is
  built on the corner path-type basis via structure constants because a path (e.g. `ab:1→3`
  for `S={1,3}`) can leave `S` in its interior yet satisfy `e·p·e=p≠0`; the subquiver algebra
  would miss it (`kA₂` dim 3 vs subquiver `k×k` dim 2). `A/AeA` IS the complement subquiver,
  certified by `dim A/AeA = dim A − dim A e_S A` from the visiting-span (interior-aware) column
  space. Both facts are pinned as `oracle_literature` tests and stated in the module docstring
  + the verification page.
- **Open risks.** (a) The Ringel universal-extension assembly (`characteristic_tilting`) is the
  riskiest construction — the plan uses iterated single-class Baer extensions and leans on the
  `Ext¹(T,∇)=0` + `is_tilting_module(T)` self-cert as the arbiter; if the cert cannot be met
  for some order the function refuses loudly (honest), and the shipped oracles are the two
  hand-derived `kA_n` cases where `T` is known independently. (b) `j_!` (`− ⊗_{eAe} eA`) reuses
  the `modules/tor.py` induction pattern on a bimodule — the grading/action is arbitrated by
  the counit iso `j^*j_!≅id`; a wrong grading fails the iso, not a silent dim match.
  (c) `presented_form` (P44) may degrade off char scope, blocking the double-Ringel-dual oracle
  — the battery pins that oracle over QQ only and states the dependency. (d) `_find_surjection`
  in the greedy Δ-peel is a bounded combination search; for non-directed algebras it may
  return `None` where a surjection exists via a larger combination — the honest verdict is
  `certified=False` (loud), and the shipped filtration oracles are directed/uniserial where a
  single basis hom surjects. (e) Mid-merge-train count drift is handled by Task G recounting
  the class table live at merge, never a guessed authoring number.
- **What is deliberately NOT in scope** (named successors, stated on the verification page):
  stratifications beyond the quasi-hereditary case; good-filtration *dimensions* past the
  directed oracles; a QPA comparison (there is no QPA surface to compare against); tilting
  MUTATION for the characteristic tilting module (P44 ships the single-complement Bongartz
  case). These are recorded honestly, never half-implemented.
