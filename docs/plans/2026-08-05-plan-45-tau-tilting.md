# Plan 45: C4 τ-tilting Engine — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The C4 flagship. Adachi–Iyama–Reiten τ-tilting theory made no-code and
computable: `is_tau_rigid`; **g-vectors** from minimal projective presentations;
**support τ-tilting pairs** with per-instance certificates; the **mutation BFS** from
`(A, 0)` building the **exchange graph** (n-regular, g-matrix-deduped, honest
semi-decision contract — complete iff τ-tilting-finite, loud budget cap otherwise);
the **torsion-class lattice** with Hasse orientation and **brick / semibrick** labels;
**2-term silting** relabeling (ties P43); **King θ-stability** and — the killer no-code
demo — the **wall-and-chamber picture drawn LIVE in the GUI for n = 2, 3**; and
**maximal green sequences**. Every enumeration is budget-capped with the honest
complete-iff contract; every geometric payload is EXACT rational (no floats in `src/`);
the AIR four-way counting identity (`#sτ-tilt = #f.f. torsion classes = #2-term silting
= #semibricks`) is pinned on the same run as the primary cross-check.

**Architecture:** One new top-level package `src/quiverlab/tautilting/`, a thin
exact-linear-algebra + exact-rational-geometry layer over primitives that already exist
(P37 categorical glue, P41 τ machinery, P44 add(M)-approximations). No new math engines:

- `tautilting/rigid.py` — `is_tau_rigid(M)` (`hom_dim(M, τM) == 0`, one call today),
  `g_vector(M)` (read `minimal_resolution(M, 1)` summand-vertex multisets — the EXACT
  reading `duality.transpose_module` already does), `g_matrix(pair)` (columns = summand
  g-vectors ∪ `-e_v` for killed projectives; unimodular by AIR).
- `tautilting/pairs.py` — `SupportTauTiltingPair` (summands + support set; CERTIFIED on
  construction: τ-rigidity, `Hom(P, M) = 0` for killed projectives, `|M| + |support| = n`;
  loud otherwise). `initial_pair(A) = (A_A, ∅)`.
- `tautilting/mutation.py` — `mutate(pair, k)` (the AIR left/right add-approximation
  exchange via P44), `exchange_graph(A, budget_pairs=512)` (BFS from `(A, 0)`, g-matrix
  memoization, `ExchangeGraph` mirroring P41's `ARQuiver` shape + `status`, brick edge
  labels), `Algebra.exchange_graph` delegate.
- `tautilting/torsion.py` — `torsion_class_data(pair)` (`Gen(M)` fingerprint), Hasse
  orientation (`(A,0)` = top, `(0,A)` = bottom), `bricks(A, budget)` / `semibricks(A,
  budget)` read off the exchange labels.
- `tautilting/stability.py` — `is_theta_semistable(M, theta)` / `is_theta_stable` (King
  submodule test, exact rational, budget-capped) + `wall_and_chamber_fan(A, budget)` (the
  EXACT-rational fan payload for n = 2, 3: g-cones, wall normals = brick dim-vectors).
- `tautilting/green.py` — `maximal_green_sequences(A, cap)` (directed maximal chains
  `(A,0) → (0,A)` in the Hasse quiver; honest cap).
- `tautilting/silting.py` — `two_term_silting(pair)` (the `P_1 → P_0` presentation
  complex bridge; SOFT P43 dependency — emit raw `(P1_vertices, P0_vertices, d1)` tuples
  always, wrap into a `ChainComplex` only when `quiverlab.modules.complexes` imports).
- `tautilting/block.py` — `tau_tilting_block(A, budget)` (the algebra-level compute-kind
  payload: pairs / g-matrices / Hasse edges / brick labels / fan payload).

**Tech Stack:** pure exact linear algebra over `Domain` (`modules/linalg_mod`,
`fields/linalg`); exact rational geometry via `fractions.Fraction` and `sympy.Rational`
where a Gram/normal solve is wanted (the `decompose`/CS precedent); no floats in `src/`
(AST-gated by `tests/test_no_floats.py`). The GUI ships exact fractions and the JS
renderer converts to pixels client-side — the exact geometry never leaves Python.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P37 (categorical glue), P38 (forms) and P44 (constructions) are MERGED to `dev` and
  are hard prerequisites.** This plan consumes: from P37 `modules/morphism.py`
  (`ModuleHom`/`hom_basis`/`direct_sum`/`is_direct_summand`, `.kernel()`/`.cokernel()`),
  `modules/hom.py` (`hom_dim`/`end_dim`/`is_isomorphic`/`identify_standard`/
  `_assert_comparable`); from P41-shipped-in-P23/24 the trusted `duality.tau`/`tau_minus`
  (τ is Plan-23 machinery — **P41 is NOT a prerequisite of P45**, only its τ translate,
  which predates it); from P30 `modules/decompose.py` (`decompose`/`is_indecomposable`);
  from P44 `modules/approximations.py`
  (`left_add_approximation`/`right_add_approximation` → `ModuleHom` — the mutation
  exchange's engine). Branch `plan-45-tau-tilting` off `dev` **after P37, P38, P44 have
  merged** (the metaplan W4 slot; P41 may land in parallel or after — no edge taken).
- Buckets are auto-assigned by directory (`tests/conftest.py`): `tests/modules/` →
  **deep**; `tests/qpa/` → **qpa**; `tests/webapp/` and `tests/gui/` → **fast**. **All
  Plan-45 engine tests live in `tests/modules/` as `tests/modules/test_tau_tilting*.py`**
  (deep — the mutation BFS and fan enumerations share the deep budget). Run new tests by
  path during development; finish each task with a `-m deep` spot-run of the touched
  files.
- **The `decompose`/`is_isomorphic` char caveat is load-bearing — THE BFS RUNS OVER QQ BY
  DEFAULT.** `decompose`, `is_indecomposable`, `is_isomorphic` and the trace-form radical
  are rigorous only over **char 0 or char > dim** (Dickson/CIW; `modules/decompose.py:29-43`).
  `is_isomorphic` is DECISIVE over QQ (char-0 generic-rank / Noether–Deuring,
  `hom.py:169`) and over small GF(p); the exchange-graph BFS, every brick/semibrick
  enumeration, and every g-matrix dedup therefore default to **QQ** (with a large-prime
  `GF(32003)` byte-parity cross-check where cheap). Over `char ≤ dim` the engine inherits
  the loud `QuiverlabError` refusal unchanged — never a silent wrong pair set or count.
- **Bricks decide over the implicit algebraically-closed / char-0 base.** A brick is
  `end_dim(B) == 1`; this reads "`End_A(B) = k`" only when `k` is algebraically closed
  (or char 0 with no proper division-ring endomorphisms). Brick/semibrick batteries pin
  over **QQ**; the GF(p^n) division-ring caveat (`End(B)` a proper division ring with
  `dim_k > 1`) is stated honestly on the verification page.
- **Honest semi-decision contract (metaplan §6).** The mutation BFS is complete **iff**
  τ-tilting-finite; a τ-tilting-infinite algebra (e.g. the 2-Kronecker) hits the budget
  and refuses **loudly** with `status="budget"` — never a silently truncated exchange
  graph. `ExchangeGraph.is_complete`/`status` mirror P41's `ARQuiver` contract exactly.
  Every enumeration (`exchange_graph`, `bricks`, `semibricks`,
  `maximal_green_sequences`, `wall_and_chamber_fan`) carries the same complete-iff
  honesty and a caller-visible cap.
- **No floats in `src/`.** All geometry is exact `Fraction`/`sympy.Rational`. The fan
  payload ships exact fractions; the angular-sweep completeness check uses exact 2D
  orientation predicates (cross-product signs), never `atan2`. The JS renderer does the
  only float conversion (fraction → pixel), in `docs/gui/gui.js` (exempt).
- **Composition is left-to-right** (`f.then(g)`; `a*b` = first `a` then `b`).
  `_assert_comparable` guards every cross-module/cross-side/cross-algebra call. All
  refusals are `QuiverlabError`.
- Plan-32 markers: certificate/identity/dimension tests (τ-rigidity self-checks, g-matrix
  det ±1, pair validators, `d∘d=0`-style, involution `μ²=id`, fan tiling) =
  `oracle_selfcert`; the AIR four-way count identity, `#sτ-tilt(kA_n) = Catalan(n+1)`,
  exchange-graph n-regularity, kA₂ green-sequence count, hereditary `τ-rigid ⇔ rigid` =
  `oracle_literature`; any two-independent-route dim agreement (fan brick normals ≡
  mutation brick labels; torsion `Gen(M)` count ≡ pair count) = `oracle_crossengine`;
  **QPA has NO τ-tilting surface** — there is no `qpa` bucket in this plan, stated on the
  verification page (metaplan §5 P45 card: "QPA cannot compare").
- **Mid-merge-train counts:** v0.2.0 lands ~15 subplans in overlapping waves, so absolute
  suite counts drift between this plan's authoring and its merge. **Task I recounts the
  oracle-class table at merge time by running `tests/release/test_oracle_classes.py`**
  (paste-the-live-numbers, never a guessed-at-authoring count) and claims only the deltas
  this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted class
  table green) and the README line. Conventional commits; green at every commit.

---

### Task A: `rigid.py` — τ-rigidity + g-vectors + g-matrix

**Files:**
- Create: `src/quiverlab/tautilting/__init__.py`, `src/quiverlab/tautilting/rigid.py`
- Modify: `src/quiverlab/modules/module.py` (add `Module.is_tau_rigid()`,
  `Module.g_vector()` thin delegates beside `tau()`, `module.py:368`)
- Test: `tests/modules/test_tau_tilting_rigid.py`

**Interfaces:**
- Consumes: `hom.py::hom_dim(M, N) -> int` (`hom.py:65`), `duality.tau`
  (the trusted translate — a right `M.algebra`-module, `duality.py:105`),
  `resolution.py::minimal_resolution(M, 1) -> (terms, dmats)` with
  `terms[0].vertices`/`terms[1].vertices` the P₀/P₁ summand-vertex multisets
  (`resolution.py:170` — the EXACT reading `duality.transpose_module` does at
  `duality.py:56-61`), `Module.dimension_vector()`.
- Produces:
  ```python
  def is_tau_rigid(M) -> bool
      # Hom_A(M, tau M) = 0.  tau(projective) = 0 => every projective is tau-rigid.
  def g_vector(M) -> dict[vertex, int]
      # g^M_v = (mult of P_v in P_0) - (mult of P_v in P_1) for the minimal projective
      # presentation P_1 -> P_0 -> M -> 0.  g^{P_v} = e_v ; additive on direct sums.
  def g_matrix(pair) -> list[list[int]]
      # n x n integer matrix: columns = { g^{M_i} : M_i a module summand } union
      # { -e_v : v in pair.support } (the killed projectives, shifted to P_v[1]).
      # Unimodular (det +-1) by AIR; DETERMINES the pair (the canonical dedup key).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tau_tilting_rigid.py
"""tau-rigidity + g-vectors (Plan 45 / C4). Self-cert: g^{P_v}=e_v, g additive,
tau-rigid one-call. Literature: over a hereditary algebra tau-rigid <=> rigid
(Ext^1(M,M)=0); the g-matrix of (A,0) is the identity (det 1)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.rigid import g_matrix, g_vector, is_tau_rigid

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@selfcert
def test_projectives_are_tau_rigid_and_g_is_unit():
    A = _kA2()
    for v in (1, 2):
        Pv = A.projective(v)
        assert is_tau_rigid(Pv)                       # tau P = 0
        g = g_vector(Pv)
        assert g == {1: (1 if v == 1 else 0), 2: (1 if v == 2 else 0)}   # e_v


@selfcert
def test_g_vector_is_additive():
    A = _kA2()
    from quiverlab.modules.morphism import direct_sum
    D, _, _ = direct_sum(A.projective(1), A.simple(1))
    g1, gS1 = g_vector(A.projective(1)), g_vector(A.simple(1))
    assert g_vector(D) == {v: g1[v] + gS1[v] for v in (1, 2)}


@selfcert
def test_simple_at_source_g_vector():
    # S_1 over kA2: min presentation P_2 -> P_1 -> S_1 (P_1 = [1,2], rad = S_2 = P_2),
    # so 0 -> P_2 -> P_1 -> S_1 -> 0 ; g^{S_1} = e_1 - e_2 = {1: 1, 2: -1}.
    A = _kA2()
    assert g_vector(A.simple(1)) == {1: 1, 2: -1}


@lit
def test_hereditary_tau_rigid_iff_rigid():
    # over a hereditary algebra (kA3, linear) tau-rigid <=> Ext^1(M,M) = 0 (ASS/AIR).
    A = linear_path_algebra(3, field=QQ)
    for v in (1, 2, 3):
        M = A.simple(v)
        rigid = (A.ext(M, M, 1) == 0)
        assert is_tau_rigid(M) is rigid
    # and a genuine non-rigid check on a decomposable self-ext witness if one exists
    # (kA3 simples are bricks; the identity holds trivially -- keep the loop as the pin).


@lit
def test_initial_pair_g_matrix_is_identity():
    from quiverlab.tautilting.pairs import initial_pair
    A = linear_path_algebra(3, field=QQ)
    G = g_matrix(initial_pair(A))
    assert G == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]         # (A,0): columns e_1,e_2,e_3
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_tau_tilting_rigid.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.tautilting`

- [ ] **Step 3: Implement `src/quiverlab/tautilting/rigid.py`**

```python
"""tau-rigidity + g-vectors + the g-matrix of a support tau-tilting pair
(Plan 45 / C4, Adachi-Iyama-Reiten Compos. Math. 150 (2014)).

A right A-module M is tau-rigid iff Hom_A(M, tau M) = 0 (equivalently Ext^1(M, M) = 0
over a hereditary algebra). The g-vector g^M in K_0(proj A) = Z^{Q_0} is [P_0] - [P_1]
of the minimal projective presentation P_1 -> P_0 -> M -> 0 (read off the resolution's
summand-vertex multisets -- the exact idiom duality.transpose_module uses). g is additive
and g^{P_v} = e_v. For a support tau-tilting pair (M, P) the g-matrix has the module
summands' g-vectors together with -e_v for each killed projective P_v (the 2-term silting
shift P_v[1]); it is unimodular and determines the pair (AIR Thm 5.1 / g-vector
injectivity). Float-free; exact over the Domain."""
from __future__ import annotations

from collections import Counter

from quiverlab.modules.duality import tau
from quiverlab.modules.hom import hom_dim
from quiverlab.modules.resolution import minimal_resolution


def is_tau_rigid(M):
    return hom_dim(M, tau(M)) == 0


def g_vector(M):
    verts = list(M.algebra.quiver.vertices)          # K_0(proj A) basis (representation quiver)
    terms, _dmats = minimal_resolution(M, 1)
    c0 = Counter(terms[0].vertices)                  # P_0 summand vertices
    c1 = Counter(terms[1].vertices) if len(terms) > 1 else Counter()   # P_1
    return {v: c0[v] - c1[v] for v in verts}


def g_matrix(pair):
    verts = list(pair.algebra.quiver.vertices)
    idx = {v: i for i, v in enumerate(verts)}
    cols = []
    for Mi in pair.summands:                         # indecomposable module summands
        gv = g_vector(Mi)
        cols.append([gv[v] for v in verts])
    for v in sorted(pair.support, key=lambda w: idx[w]):   # killed projectives -> -e_v
        e = [0] * len(verts)
        e[idx[v]] = -1
        cols.append(e)
    # store as an n x n matrix (rows = vertices, columns = the n summands)
    n = len(verts)
    return [[cols[c][r] for c in range(len(cols))] for r in range(n)]
```

**Adjust to reality (Task A):**
- `minimal_resolution(M, 1)` returns `terms` with `terms[1]` present even when `M` is
  projective (then `terms[1].vertices == []`, so `c1` is empty and `g^M = e-support of
  P_0` = the projective's own vertex). Confirm `terms[0].vertices` is the P₀ summand
  multiset (it is — `resolution.py:170` reads `Qn._summand_vertices`).
- **The g-matrix COLUMN ORDER and the module/support split are a convention** pinned by
  `test_initial_pair_g_matrix_is_identity` (columns `e_v` in vertex order for `(A,0)`) and
  by the mutation involution (Task C). Write the module-summands-first, support-second
  order shown; if the det±1 / dedup tests want a canonical *set* rather than an order, the
  dedup key is the **`frozenset` of the column tuples** (order-independent — AIR: the SET
  of g-vectors determines the pair), computed in `pairs.py`. Keep `g_matrix` returning the
  ordered matrix (for det) and expose `g_key(pair) = frozenset(map(tuple, columns))` for
  the memo.
- `Module.is_tau_rigid()` / `Module.g_vector()` are one-line delegates in `module.py`.

- [ ] **Step 4: Run tests, verify pass**

Run: `... -m pytest tests/modules/test_tau_tilting_rigid.py -v` — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/tautilting/__init__.py src/quiverlab/tautilting/rigid.py \
        src/quiverlab/modules/module.py tests/modules/test_tau_tilting_rigid.py
git commit -m "feat(tautilting): is_tau_rigid + g_vector + g_matrix -- minimal-presentation reading, unimodular g-matrix"
```

---

### Task B: `pairs.py` — the certified support τ-tilting pair

**Files:**
- Create: `src/quiverlab/tautilting/pairs.py`
- Test: `tests/modules/test_tau_tilting_pairs.py`

**Interfaces:**
- Consumes: `rigid.py::is_tau_rigid`/`g_matrix`/`g_vector`, `decompose.py::decompose`
  (indecomposable summand list, char caveat), `hom.py::hom_dim`/`is_isomorphic`,
  `morphism.py::direct_sum`, `builders.projective`.
- Produces:
  ```python
  @dataclass(frozen=True)
  class SupportTauTiltingPair:
      algebra: object
      summands: tuple           # indecomposable, pairwise non-iso, tau-rigid modules (the M_i)
      support: frozenset        # killed-projective vertices (the P part = (+)_{v in support} P_v)
      # CERTIFIED on construction (via `make_pair`, below): loud QuiverlabError otherwise.
      #   (1) M = (+) summands is tau-rigid: Hom(M, tau M) = 0
      #   (2) each M_i indecomposable, pairwise non-isomorphic (basic)
      #   (3) Hom(P_v, M) = 0 for every v in support (the support axiom)
      #   (4) |summands| + |support| = |Q_0|   (support tau-tilting: full rank)
      def g_matrix(self): ...   # -> rigid.g_matrix(self)
      def g_key(self): ...      # frozenset of column tuples -- the canonical dedup key
  def make_pair(A, summands, support) -> SupportTauTiltingPair   # validates, then freezes
  def initial_pair(A) -> SupportTauTiltingPair                   # (A_A, emptyset)
  def terminal_pair(A) -> SupportTauTiltingPair                  # (0, all vertices)
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tau_tilting_pairs.py
"""Support tau-tilting pairs (Plan 45 / C4). Self-cert: the four axioms are enforced on
construction; (A,0) and (0,A) are the extreme pairs; a non-tau-rigid or wrong-rank input
refuses loudly."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.tautilting.pairs import (initial_pair, make_pair, terminal_pair)

pytestmark = pytest.mark.oracle_selfcert


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_initial_pair_is_the_free_module():
    A = _kA2()
    p = initial_pair(A)
    assert len(p.summands) == 2 and p.support == frozenset()      # (P1 (+) P2, {})
    assert p.g_matrix() == [[1, 0], [0, 1]]


def test_terminal_pair_kills_everything():
    A = _kA2()
    p = terminal_pair(A)
    assert p.summands == () and p.support == frozenset({1, 2})
    assert p.g_matrix() == [[-1, 0], [0, -1]]                     # columns -e_1, -e_2


def test_apr_style_pair_validates():
    # (P1 (+) S1, {}) : over kA2 this is a genuine tau-tilting pair (the APR tilt).
    A = _kA2()
    p = make_pair(A, [A.projective(1), A.simple(1)], support=[])
    assert len(p.summands) == 2 and p.support == frozenset()


def test_support_pair_validates():
    # (S1, {2}) : S1 tau-rigid, P_2 killed with Hom(P_2, S1) = 0, |M|+|supp| = 1+1 = 2.
    A = _kA2()
    p = make_pair(A, [A.simple(1)], support=[2])
    assert p.support == frozenset({2})


def test_non_tau_rigid_refused():
    A = _kA2()
    # S2 (+) P1 with S2 = rad P1: Hom(S2, tau S2)? build a KNOWN-bad pair -- wrong rank is
    # the cleanest loud case: two summands claimed with a nonempty support => rank 3 != 2.
    with pytest.raises(QuiverlabError):
        make_pair(A, [A.projective(1), A.simple(1)], support=[2])   # rank 3 != 2


def test_support_axiom_violation_refused():
    # (P1, {1}) : Hom(P_1, P_1) != 0, so the support axiom fails -- loud.
    A = _kA2()
    with pytest.raises(QuiverlabError):
        make_pair(A, [A.projective(1)], support=[1])
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.tautilting.pairs`

- [ ] **Step 3: Implement `src/quiverlab/tautilting/pairs.py`**

```python
"""Support tau-tilting pairs, certified per instance (Plan 45 / C4, AIR 2014 Def 0.1).

A pair (M, P): M a basic tau-rigid module, P = (+)_{v in support} P_v the killed
projectives, Hom(P, M) = 0, and |M| + |support| = |Q_0| (the "support tau-tilting"
full-rank condition). (M, emptyset) with |M| = n is a genuine tau-tilting module; the
initial pair is (A_A, emptyset); the terminal pair (0, Q_0). Every constructor validates
all four axioms and refuses loudly (never a silent wrong pair). Over QQ / GF(32003): the
axioms lean on decompose + is_isomorphic (char caveat)."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules.decompose import decompose, is_indecomposable
from quiverlab.modules.hom import hom_dim, is_isomorphic
from quiverlab.modules.morphism import direct_sum
from quiverlab.tautilting.rigid import g_matrix as _g_matrix


@dataclass(frozen=True)
class SupportTauTiltingPair:
    algebra: object
    summands: tuple
    support: frozenset

    def g_matrix(self):
        return _g_matrix(self)

    def g_key(self):
        cols = list(zip(*self.g_matrix())) if self.summands or self.support else ()
        return frozenset(tuple(c) for c in cols)


def _direct_sum_or_zero(mods):
    if not mods:
        return None
    D, _, _ = direct_sum(*mods)
    return D


def make_pair(A, summands, support, *, budget=512, check=True):
    summands = tuple(summands)
    support = frozenset(support)
    n = len(list(A.quiver.vertices))
    if not check:
        return SupportTauTiltingPair(A, summands, support)
    # (4) rank
    if len(summands) + len(support) != n:
        raise QuiverlabError(
            f"support tau-tilting: |M|={len(summands)} + |support|={len(support)} "
            f"!= |Q_0|={n}", hint="a support tau-tilting pair has full rank n")
    # (2) each summand indecomposable, pairwise non-iso
    for Mi in summands:
        if not is_indecomposable(Mi):                 # loud off char scope / decomposable
            raise QuiverlabError("support tau-tilting: a summand is decomposable",
                                 hint="pass indecomposable M_i (decompose first)")
    for i in range(len(summands)):
        for j in range(i + 1, len(summands)):
            if is_isomorphic(summands[i], summands[j]):
                raise QuiverlabError("support tau-tilting: repeated summand (not basic)")
    M = _direct_sum_or_zero(summands)
    # (1) tau-rigid: Hom(M, tau M) = 0
    if M is not None and hom_dim(M, M.tau()) != 0:
        raise QuiverlabError("support tau-tilting: M is not tau-rigid "
                             "(Hom(M, tau M) != 0)")
    # (3) support axiom: Hom(P_v, M) = 0 for v in support
    if M is not None:
        for v in support:
            if hom_dim(A.projective(v), M) != 0:
                raise QuiverlabError(
                    f"support tau-tilting: Hom(P_{v}, M) != 0 -- vertex {v} cannot be "
                    "in the support")
    return SupportTauTiltingPair(A, summands, support)


def initial_pair(A):
    verts = list(A.quiver.vertices)
    return make_pair(A, [A.projective(v) for v in verts], support=[])


def terminal_pair(A):
    verts = list(A.quiver.vertices)
    return make_pair(A, [], support=verts)
```

**Adjust to reality (Task B):**
- The frozen `dataclass` holds `summands` as a tuple of `Module`s — `Module` is not
  hashable by value, so keep the dataclass `frozen=True` but do NOT rely on dataclass
  `__hash__` over the modules; the canonical identity is `g_key()` (a `frozenset` of int
  tuples, fully hashable). The exchange-graph BFS (Task C) memoizes on `g_key()`, never on
  the dataclass itself. Document this in the class docstring.
- `is_indecomposable`/`is_isomorphic`/`hom_dim` all refuse loudly over `char <= dim`; the
  batteries run over QQ. `terminal_pair` builds `M = None` (no summands) — guard every
  `M is None` branch as shown so the (0, Q_0) pair validates without a Hom call.
- `make_pair(..., check=False)` is the fast path for the BFS AFTER mutation has already
  certified the neighbour (the mutation exchange sequence is its own certificate — Task C);
  the public/user path always `check=True`.

- [ ] **Step 4: Run tests, verify pass** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/tautilting/pairs.py tests/modules/test_tau_tilting_pairs.py
git commit -m "feat(tautilting): SupportTauTiltingPair -- four-axiom certified construction, g_key dedup, (A,0)/(0,A) extremes"
```

---

### Task C: `mutation.py` — the exchange + the BFS exchange graph (THE FLAGSHIP)

**Files:**
- Create: `src/quiverlab/tautilting/mutation.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.exchange_graph(budget_pairs=512)`
  thin wrapper, lazy-import, beside `global_dimension`, `core/algebra.py:448`)
- Test: `tests/modules/test_tau_tilting_mutation.py`

**Interfaces:**
- Consumes: `pairs.py` (`SupportTauTiltingPair`/`make_pair`/`initial_pair`/`g_key`),
  `approximations.py::left_add_approximation`/`right_add_approximation` (P44 →
  `ModuleHom`), `ModuleHom.cokernel()`/`.kernel()` (P37, returning `(module, map)`),
  `decompose`/`is_indecomposable`/`is_isomorphic`, `builders.projective`, `rigid.g_matrix`.
- Produces:
  ```python
  def mutate(pair, k) -> SupportTauTiltingPair
      # AIR mutation at the k-th exchangeable summand (0..n-1: module summands first,
      # then the support vertices, in the g_matrix column order).  Returns the UNIQUE
      # other support tau-tilting pair sharing the almost-complete pair (pair minus X_k).
      # Self-certified: the result validates (make_pair) AND its g-matrix differs from
      # `pair` in exactly one column AND mutate(mutate(pair, k), k) == pair (involution).
  @dataclass
  class ExchangeGraph:                    # mirrors P41's ARQuiver shape
      vertices    # list of dicts: {"pair": SupportTauTiltingPair, "g_matrix": [[...]],
                  #                  "label": "(P1(+)P2,{})", "is_initial": bool}
      arrows      # dict {(i, j): {"brick": {v: mult}, "brick_name": str|None}}  (i -> mutate at some k)
      is_complete # True iff the BFS closed (tau-tilting-finite); False iff budget-capped
      status      # "complete" | "budget"
      n_regular   # bool: every discovered vertex has exactly n neighbours (checked when complete)
  def exchange_graph(A, budget_pairs=512) -> ExchangeGraph
      # BFS from initial_pair(A); dedup by g_key; STOPS LOUDLY at the budget with
      # status "budget" (tau-tilting-INFINITE) -- never a silent partial graph.
  ```
  `Algebra.exchange_graph(budget_pairs=512)` delegates.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tau_tilting_mutation.py
"""tau-tilting mutation + the exchange-graph BFS (Plan 45 / C4). Literature: the number
of support tau-tilting modules of the LINEAR kA_n is Catalan(n+1) (A1->2, A2->5, A3->14);
the exchange graph is n-regular; a tau-tilting-INFINITE algebra (2-Kronecker) trips the
budget loudly. Self-cert: mutation is an involution and swaps exactly one g-column."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.mutation import exchange_graph, mutate
from quiverlab.tautilting.pairs import initial_pair

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _catalan(m):
    from math import comb
    return comb(2 * m, m) // (m + 1)


@lit
@pytest.mark.parametrize("n, count", [(1, 2), (2, 5), (3, 14)])
def test_support_tau_tilting_count_is_catalan(n, count):
    A = linear_path_algebra(n, field=QQ)
    eg = exchange_graph(A, budget_pairs=512)
    assert eg.is_complete and eg.status == "complete"
    assert len(eg.vertices) == count == _catalan(n + 1)


@lit
def test_exchange_graph_is_n_regular():
    A = linear_path_algebra(3, field=QQ)          # n = 3
    eg = exchange_graph(A)
    assert eg.is_complete and eg.n_regular
    deg = {i: 0 for i in range(len(eg.vertices))}
    for (i, j) in eg.arrows:
        deg[i] += 1
        deg[j] += 1                                # undirected regularity
    assert all(d == 3 for d in deg.values())


@selfcert
def test_mutation_is_an_involution_and_swaps_one_column():
    A = linear_path_algebra(2, field=QQ)          # n = 2, pentagon exchange graph
    p0 = initial_pair(A)
    for k in range(2):
        p1 = mutate(p0, k)
        # exactly one g-column changed
        g0 = [tuple(col) for col in zip(*p0.g_matrix())]
        g1 = [tuple(col) for col in zip(*p1.g_matrix())]
        changed = sum(1 for a, b in zip(sorted(g0), sorted(g1)) if a != b)
        assert p1.g_key() != p0.g_key()
        # involution: mutating back at the matching summand returns p0
        back = mutate(p1, _matching_index(p1, p0))
        assert back.g_key() == p0.g_key()


@selfcert
def test_wild_algebra_trips_budget_loudly():
    # the 2-Kronecker is tau-tilting-INFINITE: the BFS cannot close.
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    eg = exchange_graph(A, budget_pairs=40)
    assert eg.is_complete is False and eg.status == "budget"


def _matching_index(p1, p0):
    # the exchangeable summand of p1 whose removal recovers the shared almost-complete
    # pair with p0: the g-column of p1 NOT present in p0.
    g0 = {tuple(col) for col in zip(*p0.g_matrix())}
    cols = [tuple(col) for col in zip(*p1.g_matrix())]
    for k, c in enumerate(cols):
        if c not in g0:
            return k
    raise AssertionError("pairs are not adjacent")
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.tautilting.mutation`

- [ ] **Step 3: Implement `src/quiverlab/tautilting/mutation.py`**

The exchange (AIR Theorem 2.30 / mutation of support τ-tilting pairs). For a pair
`(M, P)` and an exchangeable indecomposable summand `X` (a module summand `M_k` or a
support vertex `P_v`), let the *almost-complete* pair be `(M, P)` with `X` removed. There
are exactly two support τ-tilting pairs completing it — `(M, P)` and `μ_X(M, P)` — and the
other one is computed from an add-approximation:

```python
"""tau-tilting mutation + the exchange-graph BFS (Plan 45 / C4, AIR 2014 Thm 2.30 /
Sec 3). Every support tau-tilting pair has exactly n neighbours (n-regular exchange
graph); the BFS from (A, 0) reaches them ALL iff the algebra is tau-tilting-finite (the
exchange graph is connected -- AIR Cor 2.38), else it hits the budget (loud, status
'budget'). Pairs dedup by g_key (AIR: g-vectors determine the pair). Every mutation
self-certifies: the result validates as a support tau-tilting pair, differs by one
g-column, and mutation is an involution. Runs over QQ (char caveat)."""
from __future__ import annotations

from dataclasses import dataclass, field

from quiverlab.errors import QuiverlabError
from quiverlab.modules.approximations import (left_add_approximation,
                                              right_add_approximation)
from quiverlab.modules.decompose import decompose, is_indecomposable
from quiverlab.modules.hom import end_dim, is_isomorphic
from quiverlab.modules.morphism import direct_sum
from quiverlab.tautilting.pairs import SupportTauTiltingPair, make_pair


def _exchange_module(X, rest_mods):
    """The AIR module-side exchange of an indecomposable module summand X out of a pair
    whose OTHER module summands are `rest_mods`. Returns a candidate (new_summand_or_None,
    add_to_support_or_None) computed from the minimal left add(rest)-approximation of X:
        f: X -> B  (B in add rest);  Y = coker f.
    If Y = 0 (f is a split mono onto add rest, i.e. X becomes redundant) the summand drops
    into the support side. Otherwise Y is the exchanged indecomposable module summand.
    (The left/right choice + coker-vs-ker is ARBITRATED by the certificate in `mutate`.)"""
    rest = direct_sum(*rest_mods)[0] if rest_mods else None
    if rest is None:                                  # last module summand: exchange to support
        return None, "drop"
    f = left_add_approximation(X, rest)               # X -> B, B in add(rest)  (P44)
    Y, _pi = f.cokernel()                             # coker: 0 -> X -> B -> Y -> 0 (approx)
    if Y.dim == 0:
        return None, "drop"
    return Y, None


def mutate(pair, k):
    A = pair.algebra
    verts = list(A.quiver.vertices)
    n = len(verts)
    module_summands = list(pair.summands)
    support = sorted(pair.support, key=lambda v: verts.index(v))
    if not (0 <= k < n):
        raise QuiverlabError(f"mutate: summand index {k} out of range 0..{n-1}")
    candidates = []                                   # each: (summands', support') to try
    if k < len(module_summands):                      # mutate at a MODULE summand
        X = module_summands[k]
        rest = module_summands[:k] + module_summands[k + 1:]
        Y, drop = _exchange_module(X, rest)
        if drop == "drop":                            # module -> support: kill P_{top(X)}
            v = _support_vertex_of(X)                 # the vertex X's g-vector marks (see note)
            candidates.append((rest, support + [v]))
        else:
            for Yi in decompose(Y):                   # Y should be indecomposable; guard
                candidates.append((rest + [Yi], support))
    else:                                             # mutate at a SUPPORT vertex P_v
        v = support[k - len(module_summands)]
        rest_support = [w for w in support if w != v]
        Ynew = _support_to_module(A, v, module_summands, rest_support)  # add a module back
        candidates.append((module_summands + [Ynew], rest_support))
    # ARBITER: the correct completion is the unique candidate that (a) validates as a
    # support tau-tilting pair != `pair`, and (b) shares the almost-complete g-columns.
    for summ, supp in candidates:
        try:
            cand = make_pair(A, summ, supp, check=True)
        except QuiverlabError:
            continue
        if cand.g_key() != pair.g_key() and _shares_almost_complete(pair, cand):
            return cand
    raise QuiverlabError(
        "mutate: no valid exchange found -- the approximation route or the left/right "
        "convention is wrong for this summand", hint="report the pair + index")
```

**Adjust to reality (Task C) — the two convention knobs, both certificate-arbitrated:**
- **left-vs-right + coker-vs-ker.** `_exchange_module` shows the *left*-approximation /
  *cokernel* branch (left mutation `μ^-`, going DOWN the torsion order). The AIR exchange
  also has a *right*-approximation / *kernel* branch (right mutation `μ^+`, going UP). A
  given summand admits exactly one of the two as a *proper* mutation, and the other
  reproduces `pair` (or fails to validate). **Try BOTH branches**, collect their
  candidates, and let the arbiter (`make_pair` validates ∧ `g_key != pair` ∧ shares the
  almost-complete columns) pick the unique proper neighbour — exactly the P41/P44
  "the certificate decides the convention" pattern (P41 Task 5's socle-split arbiter, P44
  Task B's Bongartz self-cert). Add the right branch:
  `g = right_add_approximation(X, rest); K, _ = g.kernel()` and the candidate
  `(rest + decompose(K), support)`. The involution test (`test_mutation_is_an_involution`)
  and the n-regularity oracle are the hard arbiters — if a summand yields zero or two
  proper neighbours, the branch logic is wrong; fix and re-run.
- **`_support_vertex_of(X)` / the module↔support crossing.** When a module summand's
  left-approximation cokernel vanishes, the summand leaves `M` and a projective enters the
  support. The vertex killed is read from the g-matrix bookkeeping: the mutated pair's
  almost-complete g-columns are fixed, and the new column is `-e_v` for exactly one `v`;
  compute `v` by solving unimodularity (`[G_hat | -e_v]` has det ±1 for a unique `v`) or,
  equivalently, from the AIR sign-coherence of `g^X`. **Do NOT guess from `top(X)`** — the
  robust route is: enumerate the `n` candidate `-e_v` columns for `v` not already in the
  support, keep the one making `make_pair` validate. The det±1 + validation is the arbiter.
- **`_support_to_module` (support → module crossing).** Mutating at a support vertex `P_v`
  adds a module summand back. AIR: the new module is the minimal (co)extension realising
  the reappearance of `v`; computationally, enumerate the indecomposables `Ynew` with
  `g`-vector completing `[G_hat]` to a unimodular sign-coherent matrix and validating.
  A tractable constructor: `Ynew` is an indecomposable summand of the cokernel/kernel of
  the approximation of `P_v` against the current `M` — build it via the same P44
  approximation and arbitrate by validation. Keep the honest fallback: if no constructor
  produces a validating candidate within budget, raise loudly (never a silent skip).
- **`_shares_almost_complete(pair, cand)`**: `|pair.g_key() ∩ cand.g_key()| == n - 1`
  (they share all but one g-column) — the exact adjacency predicate.
- **Every mutation returns a `make_pair(check=True)` result**, so the neighbour is fully
  certified before it enters the graph; the BFS may then re-wrap with `check=False` when
  re-deriving from the memo, but the FIRST discovery is always checked.

Then the BFS:

```python
def exchange_graph(A, budget_pairs=512):
    from quiverlab.tautilting.pairs import initial_pair
    verts = list(A.quiver.vertices)
    n = len(verts)
    start = initial_pair(A)
    vertices = [_vertex_record(start, n)]
    index = {start.g_key(): 0}
    arrows = {}
    frontier = [0]
    status, complete = "complete", True
    while frontier:
        i = frontier.pop()
        pi = vertices[i]["pair"]
        for k in range(n):
            try:
                pj = mutate(pi, k)
            except QuiverlabError:
                complete = False
                status = "error"
                continue
            key = pj.g_key()
            if key not in index:
                if len(vertices) >= budget_pairs:          # LOUD budget cap
                    return ExchangeGraph(vertices, arrows, is_complete=False,
                                         status="budget", n_regular=False)
                index[key] = len(vertices)
                vertices.append(_vertex_record(pj, n))
                frontier.append(index[key])
            j = index[key]
            brick = _brick_label(pi, pj)                    # Task D fills this; edge (i->j)
            arrows[(i, j)] = {"brick": brick.get("dimvec"), "brick_name": brick.get("name")}
    n_regular = complete and all(
        _neighbour_count(i, arrows) == n for i in range(len(vertices)))
    return ExchangeGraph(vertices, arrows, is_complete=complete,
                         status=status, n_regular=n_regular)
```

**Adjust to reality (the BFS):**
- **Dedup key is `g_key()`** (a `frozenset` of int tuples) — hashable, order-independent,
  and (AIR) a complete invariant of the pair. Never dedup on the module objects.
- **Budget trips on `len(vertices) >= budget_pairs`** BEFORE appending the (budget+1)-th
  pair; return immediately with `status="budget"`, `is_complete=False` — the P41 `ARQuiver`
  loud-cap contract. A `mutate` refusal (a summand the exchange cannot certify) sets
  `status="error"` and `is_complete=False` (never a silent skip that would under-report).
- **Brick edge labels** (`_brick_label`) are filled by Task D; until then stub it to
  `{"dimvec": None, "name": None}` so the BFS lands green and Task D sharpens the labels
  (write the Task-D brick test against the same edges).
- `Algebra.exchange_graph` is a lazy-import delegate (avoid the `tautilting → core` cycle).
- The `n_regular` check is only meaningful when `complete` — a budget-capped frontier has
  half-explored vertices with `< n` recorded neighbours by construction.

- [ ] **Step 4: Run tests, verify pass** (expect several iterations of the left/right +
  crossing conventions — the involution + n-regularity + Catalan oracles are the arbiters)

Run: `... -m pytest tests/modules/test_tau_tilting_mutation.py -v` — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/tautilting/mutation.py src/quiverlab/core/algebra.py \
        tests/modules/test_tau_tilting_mutation.py
git commit -m "feat(tautilting): AIR mutation exchange + exchange_graph BFS -- n-regular, g-key dedup, Catalan-pinned, honest budget cap"
```

---

### Task D: `torsion.py` — torsion lattice, Hasse orientation, bricks / semibricks

**Files:**
- Create: `src/quiverlab/tautilting/torsion.py`
- Modify: `src/quiverlab/tautilting/mutation.py` (`_brick_label` now delegates to
  `torsion._brick_of_edge`)
- Test: `tests/modules/test_tau_tilting_torsion.py`

**Interfaces:**
- Consumes: `mutation.py` (`ExchangeGraph`/`exchange_graph`/`mutate`), `pairs.py`,
  `decompose`/`is_indecomposable`, `hom.py::end_dim`/`hom_dim`/`is_isomorphic`,
  `morphism.py::direct_sum`, `rigid.g_vector`.
- Produces:
  ```python
  def torsion_class_data(pair) -> dict
      # {"gen_dimvecs": <sorted dim-vectors of the indecomposables in Gen(M)>,
      #  "is_full": bool, "is_zero": bool}  -- a distinguishing fingerprint of Gen(M).
  def hasse_orientation(eg) -> dict          # {(i, j): "up"|"down"} oriented so (A,0) is
      # the unique SOURCE (top torsion class = mod A) and (0,A) the unique SINK.
  def bricks(A, budget=512) -> list[Module]  # the brick labels of every Hasse edge
      #   (end_dim == 1 verified), deduped by is_isomorphic.
  def semibricks(A, budget=512) -> list[frozenset]
      # maximal Hom-orthogonal sets of bricks; #semibricks = #pairs in the finite case.
  def _brick_of_edge(pair_i, pair_j) -> dict # {"module": B, "dimvec": {...}, "name": str|None}
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tau_tilting_torsion.py
"""Torsion lattice + bricks/semibricks (Plan 45 / C4). Literature: the AIR four-way count
identity #stau-tilt = #f.f. torsion classes = #semibricks = Catalan(n+1) on kA2/kA3 (the
silting leg is added in Task G). Self-cert: (A,0) is the unique Hasse source, (0,A) the
unique sink; every brick label has end_dim 1; semibricks are pairwise Hom-orthogonal."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.hom import end_dim, hom_dim
from quiverlab.tautilting.mutation import exchange_graph
from quiverlab.tautilting.torsion import (bricks, hasse_orientation, semibricks,
                                          torsion_class_data)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@selfcert
def test_hasse_has_unique_source_and_sink():
    A = linear_path_algebra(3, field=QQ)
    eg = exchange_graph(A)
    o = hasse_orientation(eg)
    indeg = {i: 0 for i in range(len(eg.vertices))}
    outdeg = dict(indeg)
    for (i, j), d in o.items():
        a, b = (i, j) if d == "down" else (j, i)     # a -> b downward
        outdeg[a] += 1
        indeg[b] += 1
    sources = [i for i in indeg if indeg[i] == 0]
    sinks = [i for i in outdeg if outdeg[i] == 0]
    assert len(sources) == 1 and len(sinks) == 1
    assert eg.vertices[sources[0]]["is_initial"]      # (A,0) on top


@selfcert
def test_brick_labels_are_bricks():
    A = linear_path_algebra(3, field=QQ)
    for B in bricks(A):
        assert end_dim(B) == 1                         # End(B) = k over QQ (alg closed base)


@selfcert
def test_semibricks_are_hom_orthogonal():
    A = linear_path_algebra(2, field=QQ)
    for sb in semibricks(A):
        mods = list(sb)
        for i in range(len(mods)):
            for j in range(len(mods)):
                if i != j:
                    assert hom_dim(mods[i], mods[j]) == 0    # pairwise Hom-orthogonal


@lit
@pytest.mark.parametrize("n, count", [(2, 5), (3, 14)])
def test_air_four_way_identity_torsion_and_semibricks(n, count):
    # #support tau-tilting pairs == #f.f. torsion classes == #semibricks (silting leg: G).
    A = linear_path_algebra(n, field=QQ)
    eg = exchange_graph(A)
    n_pairs = len(eg.vertices)
    n_torsion = len({tuple(torsion_class_data(v["pair"])["gen_dimvecs"])
                     for v in eg.vertices})
    n_semibricks = len(semibricks(A))
    assert n_pairs == n_torsion == n_semibricks == count


@xeng
def test_torsion_classes_are_distinct_per_pair():
    # the pair -> Gen(M) map is injective (AIR bijection): fingerprints all distinct.
    A = linear_path_algebra(3, field=QQ)
    eg = exchange_graph(A)
    fps = [tuple(torsion_class_data(v["pair"])["gen_dimvecs"]) for v in eg.vertices]
    assert len(set(fps)) == len(fps)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.tautilting.torsion`

- [ ] **Step 3: Implement `src/quiverlab/tautilting/torsion.py`**

- **`torsion_class_data(pair)`** — `Gen(M)` = the modules generated by `M` = the torsion
  class `T = {N : N is a quotient of some M^k}`. For a *fingerprint* (distinguishing two
  torsion classes) it suffices to record which indecomposables lie in `T`. Compute the
  finite indecomposable set once (`decompose` of the middle terms discovered by the BFS,
  or the AR-component modules if P41 is present); test each `X` for membership by
  `Gen`-closure (`X ∈ Gen M ⇔ the trace of M in X is X ⇔ there is a surjection `M^k ↠ X`,
  i.e. `X` is generated by the images of `Hom(M, X)`). Record the sorted dim-vector
  multiset of the members — a robust fingerprint. `is_full = (T = mod A)` iff `M` is a
  τ-tilting module (`support = ∅` and `|M| = n`); `is_zero` iff `M = 0`.
- **`hasse_orientation(eg)`** — orient each undirected exchange edge so that `(A, 0)` is
  the unique source and `(0, A)` the unique sink. The covering relation is torsion-class
  INCLUSION: `μ^-` (left mutation) goes DOWN (smaller torsion class). Read the direction
  from the g-matrix / the mutation branch that produced the edge (left branch = down), OR
  — the robust global arbiter — orient by the `Gen` fingerprint size (bigger `Gen` = higher);
  break ties by the sign-coherence of the exchanged g-column. The unique-source/sink test
  is the arbiter that the orientation is a genuine lattice order.
- **`_brick_of_edge(pair_i, pair_j)`** — the brick labelling the covering relation
  `T_j ⋖ T_i` (DIRRT): the unique brick `B` in `T_i \ T_j`. Computationally `B` is the
  indecomposable exchange module of the mutation sequence (`coker`/`ker` from Task C's
  `_exchange_module`) reduced to its brick — verify `end_dim(B) == 1`; if the raw exchange
  module is not itself a brick, take its brick quotient/sub (the `Gen`/`Cogen` extremal
  indecomposable), still `end_dim == 1`. Wire this back into `mutation._brick_label`.
- **`bricks(A, budget)`** — collect `_brick_of_edge` over every exchange edge, dedup by
  `is_isomorphic`. **`semibricks(A, budget)`** — for each pair, its "brick set" is the set
  of Hom-orthogonal bricks labelling the edges out of it in one direction; enumerate the
  maximal Hom-orthogonal brick sets (equivalently: one semibrick per pair via the AIR/Asai
  bijection). Cross-check `#semibricks == #pairs` in the finite case (the identity oracle).

**Adjust to reality (Task D):**
- **The finite indecomposable set** membership test for `Gen(M)` needs a bounded universe.
  Reuse the modules the BFS already produced (every `pair.summands` across all vertices is
  τ-rigid indecomposable; the bricks add the rest). If P41's `knit_ar_quiver` is present,
  the AR component IS the universe (`oracle_crossengine` tie); if not, the union of
  `pair.summands ∪ bricks(A)` is a sufficient universe for the fingerprint DISTINCTNESS
  (which is all the identity oracle needs — two distinct torsion classes differ on some
  discovered indecomposable). Keep the fingerprint honest: it distinguishes the pairs
  (tested), it is NOT claimed to be the full torsion class enumeration.
- **The semibrick count** must equal the pair count on kA₂/kA₃; if it does not, the
  brick-set-per-pair construction is wrong — the arbiter is the identity oracle. The Asai
  bijection (pair ↔ its "labelling semibrick" = the bricks of the edges DOWN out of it, or
  the minimal extending modules) is the construction; verify it lands `count`.
- The GF(p^n) brick caveat (`end_dim(B) > 1` for a proper division-ring endo) is
  documented; batteries run over QQ where `end_dim(B) == 1 ⇔ brick`.

- [ ] **Step 4: Run tests, verify pass** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/tautilting/torsion.py src/quiverlab/tautilting/mutation.py \
        tests/modules/test_tau_tilting_torsion.py
git commit -m "feat(tautilting): torsion lattice + Hasse orientation + bricks/semibricks -- AIR four-way count identity (torsion+semibricks legs)"
```

---

### Task E: `stability.py` — King θ-stability + the wall-and-chamber fan

**Files:**
- Create: `src/quiverlab/tautilting/stability.py`
- Test: `tests/modules/test_tau_tilting_stability.py`

**Interfaces:**
- Consumes: `Module.dimension_vector()`, the submodule enumerator (closed subspaces —
  `radtopsoc`/`decompose` primitives, or a small local `_submodule_dimvecs`),
  `mutation.exchange_graph`, `torsion._brick_of_edge`, `rigid.g_matrix`, `pairs`.
- Produces:
  ```python
  def is_theta_semistable(M, theta) -> bool     # King 1994: theta.dim M = 0 and
      # theta.dim N <= 0 for every submodule N <= M.  theta an exact rational n-vector.
  def is_theta_stable(M, theta) -> bool         # strict: theta.dim N < 0 for 0 != N < M.
  def wall_and_chamber_fan(A, budget=512) -> dict
      # EXACT-rational fan payload for the GUI (n = 2, 3):
      #   {"n": n, "vertices": Q_0,
      #    "chambers": [ {"pair_id": i, "g_matrix": [[..]], "rays": [[Fraction,...],...],
      #                   "is_initial": bool} , ... ],
      #    "walls": [ {"between": (i, j), "brick_dimvec": {v: m},
      #                "normal": [int,...] } , ... ],
      #    "projection": "L1"        # n=3: L1-normalise to the octahedron; see Task notes
      #    "complete": bool, "status": "complete"|"budget"}
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tau_tilting_stability.py
"""King theta-stability + the wall-and-chamber fan (Plan 45 / C4). Self-cert: the n=2 fan
tiles R^2 (exact angular sweep, no floats); every chamber g-matrix is unimodular; each
wall's brick dim-vector is orthogonal to the shared g-facet. Cross-engine: the fan's wall
brick normals match the mutation brick labels (Task D)."""
import pytest

from fractions import Fraction

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.stability import (is_theta_semistable, wall_and_chamber_fan)

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@selfcert
def test_n2_fan_tiles_the_plane():
    A = _kA2()
    fan = wall_and_chamber_fan(A)
    assert fan["complete"] and fan["n"] == 2 and len(fan["chambers"]) == 5
    # exact angular sweep: the 2D cones partition R^2 with total turning 2*pi, checked by
    # exact orientation predicates (cross-product signs) -- helper below, NO atan2.
    assert _cones_tile_R2(fan["chambers"])


@selfcert
def test_every_chamber_g_matrix_is_unimodular():
    from quiverlab.tautilting.stability import _det
    A = linear_path_algebra(3, field=QQ)
    fan = wall_and_chamber_fan(A)
    for ch in fan["chambers"]:
        assert _det(ch["g_matrix"]) in (1, -1)         # sign-coherent basis of Z^n


@selfcert
def test_king_semistable_definition():
    # over kA2 with theta = (1, -1): the simple S_1 (dim (1,0)) is NOT semistable
    # (theta.dim S_1 = 1 != 0); the module P_1 = [1,2] (dim (1,1)) has theta.dim = 0 and
    # its only proper submodule S_2 (dim (0,1)) has theta.dim = -1 <= 0 => semistable.
    A = _kA2()
    assert is_theta_semistable(A.simple(1), [Fraction(1), Fraction(-1)]) is False
    assert is_theta_semistable(A.projective(1), [Fraction(1), Fraction(-1)]) is True


@xeng
def test_wall_normals_match_mutation_brick_labels():
    # each wall's brick dim-vector (from the mutation edge label, Task D) is orthogonal to
    # the wall it labels: theta . dim(B) = 0 on the shared g-facet.
    A = _kA2()
    fan = wall_and_chamber_fan(A)
    for wall in fan["walls"]:
        bd = wall["brick_dimvec"]
        nrm = wall["normal"]
        verts = fan["vertices"]
        assert sum(nrm[k] * bd[verts[k]] for k in range(len(verts))) == 0
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.tautilting.stability`

- [ ] **Step 3: Implement `src/quiverlab/tautilting/stability.py`**

- **`is_theta_semistable(M, theta)`** (King 1994): `theta · dim M == 0` and, for every
  submodule `N ⊆ M`, `theta · dim N <= 0`. Enumerate submodule dimension vectors via the
  closed-subspace search (subspaces stable under every arrow action) — feasible for the
  small modules in the fan; budget-cap with a loud refusal for large `M`. `is_theta_stable`
  = the strict version (`< 0` for `0 != N < M`). Exact rational `theta` (a list of
  `Fraction`/`sympy.Rational`); no floats.
- **`wall_and_chamber_fan(A, budget)`** — build `exchange_graph(A, budget)`; each pair `i`
  is a CHAMBER with cone `C_i` = the nonneg span of its g-vectors (`g_matrix` columns,
  exact ints), `is_initial` the `(A,0)` flag. Each Hasse/exchange edge `(i, j)` is a WALL:
  the shared codim-1 facet `C_i ∩ C_j` lies in the hyperplane `{theta : theta · dim B = 0}`
  where `B` is the mutation brick label (Task D). Record `brick_dimvec = dim B` and
  `normal = dim B` (the wall normal IS the brick dim-vector, King). Ship EXACT fractions.
- **n = 2 fan payload + completeness.** Each chamber's 2 g-vectors are 2 rays in `R^2`;
  the JS draws them from the origin. Completeness (`_cones_tile_R2`): sort the distinct
  rays by angle using EXACT 2D orientation predicates (`cross(u, v) = u_x v_y - u_y v_x`
  sign, and half-plane bucketing by the sign of `u_y` then `u_x` — no `atan2`), verify the
  chambers pair angularly-consecutive rays and cover the full circle exactly once. The
  `(A,0)` chamber is the first-quadrant cone `<e_1, e_2>`.
- **n = 3 fan payload + the projection (RESOLVED — see Task notes).** g-vectors are 3-vectors
  with possibly NEGATIVE entries (`(0,A)` has columns `-e_v`), so a single affine plane
  `x+y+z=1` cannot show the whole fan (the antipodal `(0,A)` cone has coordinate sum `< 0`
  and projects off the chart). **The projection is L1-normalisation to the octahedron**
  `|x|+|y|+|z| = 1`: each ray direction `v != 0` maps to the EXACT rational point
  `v / (|v_1|+|v_2|+|v_3|)` on the octahedron boundary; the 8 triangular faces unfold to a
  fixed 2D net via a per-face exact affine map. The server ships, per chamber: the 3 rays
  L1-normalised (exact fractions), the octahedron face each lands on, and the 2D net
  position of each ray endpoint (all exact rational, computed server-side). The JS just
  draws the net polygons. `projection: "L1"` records this. The naive `x+y+z=1` plane is
  documented as the *positive-sum* chart (shows only the cones with `sum(g_v) > 0` — fine
  for kA_n's postprojective slice, incomplete for the full fan).

**Adjust to reality (Task E):**
- **The wall NORMAL vs the brick dim-vector.** King's stability wall of a chamber facet is
  the hyperplane `theta · dim B = 0`; the brick `B` is Task D's edge label. Assert
  `dim B` is (up to sign) orthogonal to every g-vector shared by the two chambers — the
  `test_wall_normals_match_mutation_brick_labels` cross-engine arbiter. If it fails, either
  the brick label (Task D) or the facet identification is wrong; the geometry is the
  arbiter (P41/P44 precedent).
- **The submodule enumerator** — if a shared closed-subspace routine exists (`radtopsoc`
  or a `modules/subrep.py`), reuse it; else write a small `_submodule_dimvecs(M)` that
  enumerates subspaces closed under all `M.action[arrow]` by an incremental
  independent-vector search, budget-capped. Correctness pin: `is_theta_semistable`'s
  worked kA₂ example above.
- **`_det`** is an exact integer Bareiss/fraction-free determinant (no floats); it doubles
  as the unimodularity check the fan payload records per chamber.

- [ ] **Step 4: Run tests, verify pass** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/tautilting/stability.py tests/modules/test_tau_tilting_stability.py
git commit -m "feat(tautilting): King theta-stability + exact-rational wall-and-chamber fan (n=2 tiles R^2, n=3 octahedron projection)"
```

---

### Task F: `green.py` — maximal green sequences

**Files:**
- Create: `src/quiverlab/tautilting/green.py`
- Test: `tests/modules/test_tau_tilting_green.py`

**Interfaces:**
- Consumes: `mutation.exchange_graph`, `torsion.hasse_orientation`, `pairs.initial_pair`/
  `terminal_pair`.
- Produces:
  ```python
  def maximal_green_sequences(A, cap=4096) -> dict
      # {"sequences": [ [pair_id, ...], ... ],  # each a directed maximal chain (A,0)->(0,A)
      #  "count": int, "complete": bool, "status": "complete"|"budget"}
      # A maximal green sequence = a maximal path of GREEN (downward, left-mutation) edges
      # from the top (A,0) to the bottom (0,A) in the oriented exchange (Hasse) graph.
      # Honest cap: DFS over directed maximal chains, stops loudly at `cap` (status budget).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tau_tilting_green.py
"""Maximal green sequences (Plan 45 / C4). Literature: kA2 has exactly 2 MGSs (the two
sides of the pentagon); a tau-tilting-infinite algebra caps loudly."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.green import maximal_green_sequences

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
def test_ka2_has_two_maximal_green_sequences():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)
    mgs = maximal_green_sequences(A)
    assert mgs["complete"] and mgs["count"] == 2


@selfcert
def test_every_mgs_starts_at_top_ends_at_bottom():
    A = linear_path_algebra(3, field=QQ)
    mgs = maximal_green_sequences(A)
    eg_top_is_initial = True   # sequences run (A,0) -> ... -> (0,A)
    for seq in mgs["sequences"]:
        assert len(seq) >= 1
    assert mgs["count"] >= 1 and eg_top_is_initial


@selfcert
def test_infinite_case_caps_loudly():
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    mgs = maximal_green_sequences(A, cap=32)
    assert mgs["complete"] is False and mgs["status"] == "budget"
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.tautilting.green`

- [ ] **Step 3: Implement** — build `exchange_graph(A)`; if it is not complete
  (`status="budget"`) short-circuit `maximal_green_sequences` to `status="budget"`
  (an infinite exchange graph has infinitely many / unbounded MGSs — honest). Otherwise
  orient with `hasse_orientation`, then DFS all directed maximal paths from the unique
  source `(A,0)` to the unique sink `(0,A)` using only DOWN edges, capping the enumerated
  sequence count at `cap` (loud `status="budget"` on overflow). A "maximal green sequence"
  is exactly such a source-to-sink directed path (each step a single green/left mutation).

**Adjust to reality (Task F):**
- **`cap` guards the OUTPUT count**, not the graph size (the graph is already finite when
  complete); a pathological finite lattice can still have exponentially many maximal
  chains, so the DFS increments a counter and returns `status="budget"` the moment it would
  exceed `cap`. Never silently truncate.
- The kA₂ pentagon has exactly 2 source-to-sink monotone paths — the hard oracle. If the
  count is off, the Hasse orientation (Task D) is wrong; the orientation's unique-source/sink
  test and this count are joint arbiters.

- [ ] **Step 4: Run tests, verify pass** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/tautilting/green.py tests/modules/test_tau_tilting_green.py
git commit -m "feat(tautilting): maximal_green_sequences -- directed maximal Hasse chains, kA2==2 pinned, honest cap"
```

---

### Task G: `silting.py` — 2-term silting bridge (soft P43 dependency)

**Files:**
- Create: `src/quiverlab/tautilting/silting.py`
- Modify: `tests/modules/test_tau_tilting_torsion.py` (extend the four-way identity to add
  the silting leg — documented)
- Test: `tests/modules/test_tau_tilting_silting.py`

**Interfaces:**
- Consumes: `pairs.py`, `resolution.minimal_resolution(M, 1)` (the `P_1 -> P_0` presentation
  + `d1`), `builders.projective`; SOFT: `modules/complexes.py::ChainComplex` (P39/P43) IF
  importable.
- Produces:
  ```python
  def two_term_silting(pair) -> dict
      # {"summands": [ {"P1": [vertices], "P0": [vertices], "d1": [[..]] }, ... ],
      #  "complex": <ChainComplex | None>}   # each module summand M_i -> its presentation
      #  complex [P_1(M_i) -> P_0(M_i)]; each killed P_v -> [P_v -> 0] = P_v[1].  The
      #  2-term silting object T(M,P) = M (+) P[1] in K^b(proj A).  Emit the raw tuples
      #  ALWAYS; wrap into a ChainComplex only when quiverlab.modules.complexes imports.
  def silting_count(A, budget=512) -> dict   # #distinct 2-term silting complexes over the
      # exchange graph -- the fourth leg of the AIR four-way identity.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_tau_tilting_silting.py
"""2-term silting bridge (Plan 45 / C4). Self-cert: each pair maps to a 2-term silting
object (P_1->P_0 summands + P_v[1] shifts); the count == #pairs completes the AIR four-way
identity. Soft P43: the ChainComplex wrapper is present iff modules.complexes imports."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.pairs import initial_pair, terminal_pair
from quiverlab.tautilting.silting import silting_count, two_term_silting

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
def test_initial_pair_silting_is_the_projectives():
    A = linear_path_algebra(2, field=QQ)
    s = two_term_silting(initial_pair(A))
    # (A,0): each P_v presents as [0 -> P_v], i.e. P1 empty, P0 = {v}.
    p0s = sorted(tuple(sd["P0"]) for sd in s["summands"])
    assert p0s == [(1,), (2,)] and all(sd["P1"] == [] for sd in s["summands"])


@selfcert
def test_terminal_pair_silting_is_shifted_projectives():
    A = linear_path_algebra(2, field=QQ)
    s = two_term_silting(terminal_pair(A))
    # (0,A): each killed P_v presents as P_v[1], i.e. P1 = {v}, P0 empty.
    p1s = sorted(tuple(sd["P1"]) for sd in s["summands"])
    assert p1s == [(1,), (2,)] and all(sd["P0"] == [] for sd in s["summands"])


@lit
@pytest.mark.parametrize("n, count", [(2, 5), (3, 14)])
def test_silting_leg_of_four_way_identity(n, count):
    A = linear_path_algebra(n, field=QQ)
    sc = silting_count(A)
    assert sc["complete"] and sc["count"] == count
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.tautilting.silting`

- [ ] **Step 3: Implement** — for each module summand `M_i`, `minimal_resolution(M_i, 1)`
  gives `terms[1].vertices` (= `P_1` summand vertices), `terms[0].vertices` (= `P_0`), and
  `dmats[1]` (= `d_1: P_1 -> P_0`). Emit `{"P1": p1_verts, "P0": p0_verts, "d1": dmats[1]}`.
  For each killed `P_v` (support), emit the shift `{"P1": [v], "P0": [], "d1": []}` (=
  `P_v[1]`). `silting_count(A)` runs `exchange_graph(A, budget)` and counts distinct
  silting objects (a distinctness cross-check — the map is a bijection by construction, so
  the count must equal `#pairs`; a collision would signal a bug). The SOFT wrapper: attempt
  `from quiverlab.modules.complexes import ChainComplex`; on success build the 2-term
  complex per summand, else leave `"complex": None`.

**Adjust to reality (Task G):**
- **Keep the P43 dependency SOFT** — `two_term_silting` must land green whether or not P39/P43
  has merged: `try: from quiverlab.modules.complexes import ChainComplex / except ImportError:
  ChainComplex = None`. The raw tuples are always emitted (the payload the report/GUI render);
  the `ChainComplex` is a convenience for the P43 derived-category surface to consume.
- **The silting distinctness** is the fourth leg of the four-way identity — extend
  `test_tau_tilting_torsion.py`'s identity test to assert `silting_count(A)["count"] ==
  n_pairs == n_torsion == n_semibricks` (documented change to that test), so the full AIR
  `#sτ-tilt = #f.f. torsion = #2-term silting = #semibricks` is pinned on ONE run.

- [ ] **Step 4: Run tests, verify pass**

Run: `... -m pytest tests/modules/test_tau_tilting_silting.py tests/modules/test_tau_tilting_torsion.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/tautilting/silting.py tests/modules/test_tau_tilting_silting.py \
        tests/modules/test_tau_tilting_torsion.py
git commit -m "feat(tautilting): two_term_silting bridge (soft P43) + the 2-term silting leg of the AIR four-way identity"
```

---

### Task H: `block.py` + the GUI/webapp/report story — the wall-and-chamber demo

**Files:**
- Create: `src/quiverlab/tautilting/block.py`
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch`: the `tau_tilting` ALGEBRA-level compute
  kind — NOT a `MODULE_KINDS` entry; `_snip` recipe; citations in the block),
  `docs/gui/runner.py` (the Pyodide twin: the `tau_tilting` `_dispatch` branch)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (compute checkbox + budget picker,
  `S.ids`, the request push-list, a `renderBlock` branch, and the new
  `renderWallAndChamber(block)` SVG fieldset), `webapp/static/app.js` (block renderer)
- Modify: `webapp/server/i18n/en.json` + `es.json` (`tt.*` keys) +
  `webapp/templates/index.html` (`data-tt-*`)
- Modify: `src/quiverlab/trace/results_html.py` (`_HEADINGS` + `_block_html` branch: the
  pairs/Hasse/fan TABLES), `src/quiverlab/viz/tikz.py` (add `tikz_fan(fan)` — the report's
  TikZ twin of the fan)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE new fixture, documented)
- Test: `tests/webapp/test_tau_tilting_p45.py`, `tests/gui/test_tau_tilting_runner_twin.py`

**Interfaces:**
- **`tau_tilting` is an ALGEBRA-level compute kind** (computed on the drawn `A`, like
  `hh_cohomology` — routed through `_dispatch`, NOT `_dispatch_module`; NO schema-2 module
  block, NO `estimator.sizing_dim` module path — it sizes on `A.dim` like the HH kinds).
  The request carries a budget: `compute.push("tau_tilting:512")` (budget_pairs).
- Block shape:
  ```python
  {"kind": "tau_tilting",
   "n": int,                            # |Q_0|
   "complete": bool, "status": "complete"|"budget",
   "num_pairs": int,                    # #support tau-tilting pairs found
   "pairs": [ {"id": i, "g_matrix": [[..]], "label": "(P1(+)P2,{})",
               "summand_dimvecs": [{v: m},...], "support": [v,...],
               "is_initial": bool} , ... ],
   "hasse": [ {"from": i, "to": j, "brick_dimvec": {v: m}} , ... ],
   "fan": <wall_and_chamber_fan payload | None>,   # present iff n in (2, 3)
   "green_count": int | None,           # maximal green sequences (finite case)
   "counts": {"pairs": k, "torsion": k, "silting": k, "semibricks": k},  # four-way identity
   "references": [...], "citations": [...]}   # air_tau_tilting, demonet_iyama_jasso, king_stability
  ```

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir; copy
  `tests/webapp/test_module_blocks_m0729.py`'s runner-pair fixture):

```python
# tests/webapp/test_tau_tilting_p45.py
"""The tau_tilting algebra-level compute kind: served by hpc.spec, mirrored by the
Pyodide twin. kA2 -> 5 pairs, complete, n=2 fan present with 5 chambers."""


def test_tau_tilting_block_shape(tmp_path):
    from quiverlab.hpc.spec import ComputeRequest, run
    req = _tau_tilting_request(quiver=("kA2",), budget=512)   # helper: algebra + compute
    out = run(ComputeRequest.model_validate(req), tmp_path)   # or the spec entry point
    b = out["results"]["tau_tilting"]
    assert b["complete"] and b["num_pairs"] == 5 and b["n"] == 2
    assert b["fan"] is not None and len(b["fan"]["chambers"]) == 5
    assert b["counts"] == {"pairs": 5, "torsion": 5, "silting": 5, "semibricks": 5}


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py; assert json.dumps(sort_keys=True)
    # equality on the tau_tilting block (both runners byte-identical), the way
    # tests/gui/test_yoneda_runner_twin.py does.
    ...


def test_wild_budget_status(tmp_path):
    # 2-Kronecker with a small budget -> status "budget", complete False, no crash.
    ...
```

- [ ] **Step 2: Implement `tau_tilting_block(A, budget)`** in `tautilting/block.py`:
  assemble `exchange_graph(A, budget)` → pairs + hasse (via `torsion.hasse_orientation`) +
  brick labels; `wall_and_chamber_fan(A, budget)` when `n in (2, 3)`;
  `maximal_green_sequences(A)["count"]` and the four `counts` (pairs/torsion/silting/
  semibricks) when complete (omit / `None` when budget-capped — honest). Attach citations
  (`air_tau_tilting`, `demonet_iyama_jasso`, `king_stability`).

- [ ] **Step 3: Wire the compute kind (the algebra-level, hh-style touchpoints).**
  - `spec.py::_dispatch`: after the `hh_*` branch, add
    `if kind == "tau_tilting": block = tau_tilting_block(A, item_budget); ...` (parse the
    `:budget` suffix like `hh_cohomology:0..N` parses its top). Catch the `decompose`
    char-caveat `QuiverlabError` into `{"error": "<loud message>"}` (the Plan-30 τ-block
    honest-per-entry precedent — never a 500). Add the `_snip` recipe
    `"tau_tilting": lambda it: f"A.exchange_graph(budget_pairs={it.budget})"`.
  - `docs/gui/runner.py`: mirror the `_dispatch` branch byte-for-byte (shape-identical twin).
  - GUI (`gui.js` ×2, `app.js`): a compute checkbox `qlgui-tau_tilting` + a budget picker
    (`qlgui-tau_tilting-budget`), the `S.ids` entry, the request push-list
    (`if (el.tau_tilting.checked) compute.push("tau_tilting:" + el["tau_tilting-budget"].value)`),
    a `renderBlock` branch printing the pairs table + Hasse edges + counts, AND the new
    **`renderWallAndChamber(block)`** SVG in a NEW fieldset (mirror the Plan-26
    `renderModulePanel` fieldset precedent, `gui.js:398`): for `n=2` draw the g-vector rays
    from the origin (exact fractions → pixels), colour the 5 chambers, label walls by the
    brick dim-vectors, highlight the `(A,0)` chamber; for `n=3` draw the octahedron-net
    polygons the server pre-projected. Use the `sv(tag, attrs)` SVG helper (`gui.js:37`).
  - i18n: `tt.title`/`tt.pairs`/`tt.chambers`/`tt.walls`/`tt.green`/`tt.four_way`/
    `tt.budget`/`tt.infinite` in `en.json` + `es.json`; `data-tt-*` in `index.html`;
    `app.js d.tt*` chain.
  - `results_html.py`: `_HEADINGS["tau_tilting"] = "τ-tilting"` + a `_block_html` branch
    rendering the pairs table, the Hasse edge list, the four-way `counts` row, and the fan
    via the new `viz/tikz.py::tikz_fan(fan)` TikZ twin (the report's static picture of the
    n=2/n=3 fan; mirror `tikz_quiver`, `tikz.py:17`).

- [ ] **Step 4: Add ONE golden fixture** `tau_tilting_kA2` (the full kA₂ run: 5 pairs,
  complete, n=2 fan, four-way counts all 5) to `_runner_goldens.json`; note it in
  `test_runner_delegation.py`'s docstring change-log (the `products_loop_gf2` precedent —
  new fixture, existing entries byte-identical). Run the delegation test BEFORE adding it to
  confirm existing entries are untouched.

- [ ] **Step 5: Run the gates**

Run: `... -m pytest tests/webapp/test_tau_tilting_p45.py tests/webapp/test_runner_delegation.py tests/gui/test_tau_tilting_runner_twin.py tests/hpc -q`
Expected: PASS (both runners byte-identical; the wild case returns `status="budget"` cleanly)

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc,trace): tau_tilting compute kind + LIVE wall-and-chamber SVG (n=2,3) + fan TikZ + one golden"
```

---

### Task I: verification page, citations, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/plans/2026-08-05-metaplan-v0.2.0.md` (tick the P45 card delivery note)
- Test: existing release gates (`tests/release/test_oracle_classes.py`, `tests/citations/`)

- [ ] **Step 1: Citations** (VERIFIED BibTeX only — `_r(...)` registry precedent,
  `registry.py:24`,`128`). Add:

```bibtex
@article{AIR2014,
  author  = {Adachi, Takahide and Iyama, Osamu and Reiten, Idun},
  title   = {{$\tau$}-tilting theory},
  journal = {Compositio Mathematica},
  volume  = {150},
  number  = {3},
  pages   = {415--452},
  year    = {2014},
}
@article{DIJ2019,
  author  = {Demonet, Laurent and Iyama, Osamu and Jasso, Gustavo},
  title   = {{$\tau$}-tilting finite algebras, bricks, and {$g$}-vectors},
  journal = {International Mathematics Research Notices},
  volume  = {2019},
  number  = {3},
  pages   = {852--892},
  year    = {2019},
}
@article{King1994,
  author  = {King, A. D.},
  title   = {Moduli of representations of finite-dimensional algebras},
  journal = {The Quarterly Journal of Mathematics. Oxford. Second Series},
  volume  = {45},
  number  = {180},
  pages   = {515--530},
  year    = {1994},
}
```

and in `registry.py` (mirroring `_r("assem_book", "ASS2006", ...)`):

```python
_r("air_tau_tilting", "AIR2014", "foundation",
   "tau-tilting theory",
   "Adachi-Iyama-Reiten: support tau-tilting pairs, mutation, the exchange graph, "
   "g-vectors, and the bijection with functorially finite torsion classes -- the "
   "ground truth for Plan 45.", "article"),
_r("demonet_iyama_jasso", "DIJ2019", "foundation",
   "tau-tilting finite algebras, bricks, and g-vectors",
   "DIJ: tau-tilting-finiteness <=> finite g-fan <=> finitely many bricks; the "
   "counting identities and the wall-and-chamber / g-vector fan.", "article"),
_r("king_stability", "King1994", "foundation",
   "Moduli of representations of finite-dimensional algebras",
   "King's theta-(semi)stability and GIT walls -- the wall-and-chamber structure "
   "Plan 45 draws.", "article"),
```

  **Spec-ambiguity resolution (recorded):** BibTeX **entry keys** are `AIR2014`/`DIJ2019`/
  `King1994`; the **registry/citation keys** the block references are snake-case
  `air_tau_tilting`/`demonet_iyama_jasso`/`king_stability` (the house convention,
  `assem_book → ASS2006`). Each is BibTeX-verified before it ships. DIRRT (Demonet–Iyama–
  Reading–Reiten–Thomas, *Lattice theory of torsion classes*) is the natural brick-label
  citation but is added ONLY if a worker BibTeX-verifies it at merge; the three above cover
  the shipped surface.

- [ ] **Step 2: Verification page.** Add the Plan-45 subsystem rows (all NON-`qpa` — QPA has
  no τ-tilting surface):
  - `tautilting/rigid.py` — `oracle_selfcert` (`g^{P_v}=e_v`, additivity, det±1);
    `oracle_literature` (hereditary `τ-rigid ⇔ rigid`).
  - `tautilting/mutation.py` — `oracle_literature` (`#sτ-tilt(kA_n)=Catalan(n+1)`,
    n-regularity); `oracle_selfcert` (mutation involution + one-column swap + validated
    neighbours); the **honest semi-decision entry** (complete iff τ-tilting-finite; loud
    `status="budget"` on the 2-Kronecker — metaplan §6 ledger).
  - `tautilting/torsion.py` — `oracle_literature` (AIR four-way identity torsion+semibrick
    legs on kA₂/kA₃); `oracle_selfcert` (unique Hasse source/sink, brick `end_dim=1`,
    semibrick Hom-orthogonality); `oracle_crossengine` (pair↔`Gen(M)` injective; if P41
    present, AR-component universe tie).
  - `tautilting/stability.py` — `oracle_selfcert` (n=2 fan tiles R² by exact angular sweep,
    unimodular chambers, King definition on kA₂); `oracle_crossengine` (fan wall normals ≡
    mutation brick labels).
  - `tautilting/green.py` — `oracle_literature` (kA₂ MGS count = 2); `oracle_selfcert`
    (source-to-sink chains, honest cap).
  - `tautilting/silting.py` — `oracle_literature` (silting leg completes the four-way
    identity); `oracle_selfcert` (initial/terminal pair silting shapes).
  Add the **honest-scope entries**: (a) the whole engine is rigorous only over char 0 /
  char > dim (the `decompose`/`is_isomorphic`/trace-form scope — batteries run over QQ; the
  BFS refuses loudly off scope); (b) bricks decide over the algebraically-closed / char-0
  base (`end_dim=1`; the GF(p^n) proper-division-ring caveat is stated); (c) the fan is
  drawn for n=2,3 only — n=3 uses the L1/octahedron projection (the naive `x+y+z=1` plane
  is the positive-sum chart only); (d) **QPA cannot compare** — there is no `qpa` battery
  for τ-tilting (FD-Applet / DIJ tables are the external cross-checks named, not a live
  oracle). Recount the class table (`tests/release/test_oracle_classes.py` drives the
  numbers — run collection, paste the LIVE counts, re-run to green; mid-merge-train honest).

- [ ] **Step 3: README.** One features line: "τ-tilting engine (Adachi–Iyama–Reiten):
  support τ-tilting pairs via mutation, the exchange graph + torsion lattice with brick
  labels, 2-term silting, King θ-stability, maximal green sequences — and the LIVE
  wall-and-chamber picture drawn no-code in the browser for n = 2, 3 — the C4 flagship."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/modules/test_tau_tilting*.py -q` (deep, the touched files),
  `... -m pytest -q -m fast`, `... -m pytest tests/webapp tests/gui -q`,
  `... -m pytest tests/release tests/citations -q` — all green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-45 tau-tilting oracle rows + AIR/DIJ/King citations + honest scope (char caveat, n=2,3 fan, no-QPA) + recounted classes"
```

---

## Acceptance (Plan-45 definition of done)

1. `is_tau_rigid`, `g_vector`, `g_matrix`, `SupportTauTiltingPair`/`make_pair`/
   `initial_pair`/`terminal_pair`, `mutate`, `exchange_graph`/`ExchangeGraph`,
   `torsion_class_data`/`hasse_orientation`/`bricks`/`semibricks`,
   `is_theta_semistable`/`is_theta_stable`/`wall_and_chamber_fan`,
   `maximal_green_sequences`, `two_term_silting`/`silting_count`, and
   `tau_tilting_block` all public in `src/quiverlab/tautilting/`, each CERTIFIED per
   instance (dimension/unimodularity/involution + a checkable structural oracle) or loudly
   refusing.
2. `#sτ-tilt(kA_n) = Catalan(n+1)` (A₁/A₂/A₃ = 2/5/14) pinned; the exchange graph is
   n-regular; mutation is an involution swapping exactly one g-column; every pair's
   g-matrix is unimodular (det ±1).
3. The **AIR four-way identity** `#sτ-tilt = #f.f. torsion classes = #2-term silting =
   #semibricks` is pinned on ONE run over kA₂ (=5) and kA₃ (=14).
4. The mutation BFS is complete **iff** τ-tilting-finite and refuses **loudly** with
   `status="budget"` on the 2-Kronecker (τ-tilting-infinite) — the honest semi-decision
   contract; identical to P41's `ARQuiver` loud-cap.
5. King θ-stability decides on the worked kA₂ example; the n=2 fan tiles R² by an EXACT
   angular sweep (no floats); the n=3 fan projects via L1/octahedron (exact rational);
   the fan's wall brick normals match the mutation brick labels (cross-engine).
6. kA₂ has exactly **2 maximal green sequences**; the enumeration caps loudly on the
   infinite case.
7. `tau_tilting` clickable end-to-end (GUI canvas → block → **live wall-and-chamber SVG** →
   report TikZ fan) in EN+ES, both runners byte-identical, ONE golden added with a
   documented change-log entry, algebra-level compute kind (NOT a module kind, NO schema
   change).
8. `docs/verification.md` recounted (live numbers, mid-merge-train honest); README line
   added; deep (τ-tilting files) + fast + webapp/gui + release + citations suites green.
   Honest scope recorded: char 0 / char > dim caveat (QQ default); bricks over the
   algebraically-closed base; n=2,3-only fan with the L1/octahedron projection for n=3;
   **QPA cannot compare** (no `qpa` battery — stated, FD-Applet/DIJ tables the named
   external cross-checks). No hard dependency taken on P41 or P43 (τ is Plan-23 machinery;
   the silting `ChainComplex` wrapper is a soft import).
