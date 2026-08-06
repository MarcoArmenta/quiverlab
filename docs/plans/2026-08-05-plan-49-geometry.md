# Plan 49: C8 Geometry of Representations — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The C8 axis: the geometry of the representation variety, no-code and
exact. Orbit dimensions in `Rep(Q, d)` (`dim O_M = dim GL(d) − dim End(M)`);
**Voigt rigidity** (`is_rigid(M) ⇔ Ext¹(M,M) = 0`, and `rigidity_codim` with
HONEST semantics — the Voigt codim equality on hereditary, an upper bound on
`kQ/I`); the **Kac canonical decomposition** of a dimension vector over a
hereditary DYNKIN algebra (Schofield / Derksen–Weyman root theory, certified per
instance by rigidity of the generic module); and the **degeneration / hom order**
for representation-finite algebras (Zwara–Bongartz, a Hasse poset built over the
P41 AR indecomposables, mirroring `ARQuiver` with viz twin renderers). One no-code
GUI compute kind — `orbit_geometry` — puts orbit dim, the rigidity verdict, the
honest codimension, and (on hereditary Dynkin) the canonical decomposition one
click away. **Hall numbers stay out** (they are the P3 axis). Every enumeration
is budget-capped with the honest complete-iff contract; every number is exact
(int/dict — no floats).

**Architecture:** Two thin exact layers over primitives that already exist (P30
`decompose`, P37 `hom`/`morphism`, P38 forms/roots, P41 AR knitting):

- **`src/quiverlab/invariants/geometry.py`** (new) — the single-module /
  dimension-vector INVARIANTS: `orbit_dimension(M)`, `group_dim`,
  `representation_variety_dim(A, d)`, `is_rigid(M)`, `rigidity_codim(M)`,
  `canonical_decomposition(A, d)`, and the shared GUI block builder
  `orbit_geometry_block(M)`. These size on a single module / a single dimension
  vector and are cheap on the oracle inputs, so their batteries live in
  `tests/invariants/` (the **fast** bucket). The orbit/rigidity core computes
  over EVERY exact `Domain` (`end_dim` + `Algebra.ext` are exact everywhere,
  Plan 19/29); `canonical_decomposition` is a hereditary-Dynkin notion and leans
  on `positive_roots` + `ar_quiver` + `is_isomorphic` (QQ / char-scope).
- **`src/quiverlab/modules/degeneration.py`** (new) — `degeneration_order(A, d)`
  + the `DegenerationPoset` object. **Justification for the `modules/` home
  (not `invariants/`):** this verb ENUMERATES and COMPARES whole iso-classes of
  modules — it consumes the P41 `knit_ar_quiver` indecomposable universe and the
  `hom_dim`/`direct_sum` module machinery, and is naturally heavy (the AR BFS +
  a multiset enumeration + a pairwise Hom matrix). `tests/modules/` is the
  **deep** bucket (`tests/conftest.py`), which is exactly where this belongs,
  and it sits alongside `modules/ar.py` whose `ARQuiver` shape it mirrors. The
  single-module invariants (orbit dim, rigidity) are genuinely `invariants/`
  citizens; the poset over the whole module category is a `modules/` citizen.
- **`src/quiverlab/viz/`** (extend) — a generic Hasse/poset layout
  (`viz/layout.py::poset_layout`) + its two renderers `viz/tikz.py::tikz_hasse`
  (report/PDF) and `viz/hasse_html.py::hasse_svg` (HTML/GUI) — the "viz twin"
  the `DegenerationPoset` (and, later, P45's exchange lattice) render through.
  All coordinates exact `int`/`Fraction` (the `viz/layout.py` precedent); the
  only float conversion is client-side pixels.

No new math engines; every result self-certifies (orbit-dim identity, Voigt codim
identity on hereditary, canonical rigidity, poset axioms) and is cross-checked
against P38 roots, our own `Ext`/`Hom` engines, hand-derived posets, and a
probe-first QPA battery.

**Tech Stack:** pure exact linear algebra over `Domain` (`modules/hom.py`
`hom_dim`/`end_dim`, `modules/morphism.py` `direct_sum`, `Algebra.ext`),
P38 `positive_roots`/`form_type`/`tits_form`/`is_hereditary`, P41
`knit_ar_quiver`/`ARQuiver`, P30 `decompose`/`is_indecomposable`; exact geometry
via `fractions.Fraction` in the viz layout (the `viz/layout.py` precedent). No
floats in `src/` (AST-gated by `tests/test_no_floats.py`).

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P30, P37, P38, P41 are MERGED to `dev` and are the hard prerequisites**
  (verified at authoring: `f261878 Merge plan-41-ar-completion`,
  `87f3f0a Merge plan-38-forms-recognition`, `80db0b4 Merge plan-37-categorical-glue`,
  P30 in `modules/decompose.py`). This plan consumes, at the signatures verified
  in `dev`:
  - P38 (via `Algebra` delegates, `core/algebra.py`): `Algebra.positive_roots()`
    → sorted list of **vertex-order tuples** (raises `QuiverlabError` unless
    hereditary of Dynkin A/D/E type); `Algebra.form_type()` →
    `"finite"|"tame"|"wild"`; `Algebra.dynkin_type()` → tuple `("A",n)` / `("~A",n)`
    / `None`; `Algebra.tits_form(d)` = `⟨d,d⟩` (d a vertex-order list);
    `Algebra.euler_form(d, e)`; `Algebra.is_hereditary()`.
  - P41: `Algebra.ar_quiver(budget_modules=256, budget_dim=4096)` →
    `ARQuiver` with `.vertices` (list of dicts `{"name","dimvec","module"}` —
    **each vertex carries the actual `Module`**), `.arrows` `(i,j)→mult`,
    `.is_complete`, `.status ∈ {"complete","budget","error","unsupported"}`.
    **`knit_ar_quiver` REFUSES self-injective input** (`status="unsupported"`) —
    `degeneration_order` inherits that scope (loud).
  - P37: `modules/hom.py` `hom_dim(M,N)` / `end_dim(M)` / `is_isomorphic` /
    `identify_standard` / `_assert_comparable`; `modules/morphism.py`
    `direct_sum(*mods) → (D, incls, projs)`, `hom_basis`, `ModuleHom`.
  - P30: `modules/decompose.py` `decompose(M, budget) → [(M_i, m_i), …]`
    (**tuples**, not bare modules), `is_indecomposable(M, budget)`.
  - Core: `Algebra.ext(M, N, n) → int`, `Algebra.simple/projective/injective(v, side=)`,
    `Module.dimension_vector() → dict`, `A.quiver.vertices`, `A.quiver.arrows`
    (dict `label → (s, t)`).
  Branch `plan-49-geometry` off `dev` (the metaplan W5 slot; P49 needs P38, and
  the finite-type constructions here additionally consume the P41 AR universe —
  both merged, no wave edge broken).
- **Char scope is load-bearing.** `orbit_dimension`, `is_rigid`,
  `rigidity_codim` are exact over EVERY `Domain` (`end_dim`, `Algebra.ext`).
  `canonical_decomposition` and `degeneration_order` lean on
  `is_isomorphic`/`identify_standard`/`decompose` and the AR knit, rigorous only
  over **char 0 or char > dim** (Dickson/CIW, `decompose.py:321`); their
  batteries run over **QQ** (with a `GF(32003)` byte-parity cross-check on the
  orbit/rigidity core where cheap). Over `char ≤ dim` the engine inherits the
  loud `QuiverlabError` refusal — never a silent wrong decomposition or poset.
- **Honest semi-decision contract (metaplan §6).** `degeneration_order` is
  complete **iff** the algebra is representation-finite (`ar_quiver` closes); a
  representation-infinite or self-injective input returns
  `is_complete=False` with `status ∈ {"budget","unsupported","error"}` — never a
  silently truncated poset. `DegenerationPoset.is_complete`/`status` mirror
  P41's `ARQuiver` contract exactly. `canonical_decomposition` is complete for
  hereditary DYNKIN and refuses **loudly** for hereditary Euclidean/wild (a named
  successor) and for non-hereditary `kQ/I`. Every enumeration carries a
  caller-visible budget.
- **No floats in `src/`.** Orbit dims / codims / multiplicities are `int`;
  dimension vectors `dict[vertex,int]`; the Hasse layout uses
  `fractions.Fraction`. The JS renderer does the only fraction→pixel conversion
  (`docs/gui/*.js`, exempt).
- **Composition is left-to-right** (`f.then(g)`). `_assert_comparable` guards
  every cross-module/cross-side/cross-algebra call. All refusals are
  `QuiverlabError`.
- Plan-32 markers: the orbit-dim identity, Voigt codim identity, canonical
  rigidity certificate, and the poset axioms (reflexive/antisymmetric/transitive,
  cover monotonicity) = `oracle_selfcert`; two-independent-route agreement
  (poset maximum ≡ `canonical_decomposition`; codim via `dim Rep − dim O_M` ≡
  `dim Ext¹(M,M)` on hereditary; canonical components ≡ `positive_roots`) =
  `oracle_crossengine`; the hand-derived posets (kA₂/kA₃), the `(2,1)=P₁⊕S₁`
  Kac pin, and "every Dynkin indecomposable is rigid" = `oracle_literature`;
  the QPA probe lives in `tests/qpa/` (bucket = the class, never double-marked).
- **Mid-merge-train counts:** v0.2.0 lands ~15 subplans in overlapping waves, so
  absolute suite counts drift between authoring and merge. **Task 7 recounts the
  oracle-class table at merge time by running
  `tests/release/test_oracle_classes.py`** (paste the live numbers, never a
  guessed-at-authoring count) and claims only the deltas this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted
  class table green, `tests/release/test_oracle_classes.py` green) and adds its
  citations to `citations/references.bib` + `registry.py` (the `bibtex()` helper
  hard-fails if the two are out of sync). Conventional commits; green at every
  commit.

### Mathematical foundation (derived — the plan's ground truth)

*Orbit dimension (quiver reps).* For a quiver `Q`, dimension vector `d`, the
representation variety is `Rep(Q, d) = ⊕_{a: i→j} Mat_{d_j × d_i}(k)`, so
`dim Rep(Q, d) = Σ_{a: i→j} d_i d_j`. The base-change group
`GL(d) = ∏_v GL_{d_v}` has `dim GL(d) = Σ_v d_v²` and acts with stabiliser
`Aut(M)`, a Zariski-open subset of the affine space `End(M)`. Hence
`dim O_M = dim GL(d) − dim End(M) = Σ_v d_v² − dim_k End_A(M)`. This formula holds
verbatim for `kQ/I`: an `A`-module endomorphism IS a quiver-rep endomorphism
respecting the relations, so `End_A(M)` is the same object.

*Codimension & Voigt.* `codim_{Rep(Q,d)} O_M = dim Rep − dim O_M
= dim End(M) − ⟨d,d⟩`, where `⟨d,d⟩ = Σ_v d_v² − Σ_{a:i→j} d_i d_j` is the Euler
form. Over a **hereditary** algebra `⟨d,d⟩ = dim End(M) − dim Ext¹(M,M)`, so
`codim O_M = dim Ext¹(M,M)` (Voigt). For a general `kQ/I` the module variety
`mod_A(d) ⊆ Rep(Q,d)` is cut by the relations; the scheme tangent-space bound
gives `codim_{mod_A(d)} O_M ≤ dim Ext¹_A(M,M)` (equality iff smooth at `M`).
Thus **`rigidity_codim(M) := dim Ext¹(M,M)` is the codim on hereditary and an
UPPER BOUND on `kQ/I`** — never claim equality off hereditary. Voigt's lemma
(`Ext¹(M,M)=0 ⇒ O_M open`) holds in general, so `is_rigid ⇒ open orbit` always;
the converse needs smoothness (holds on hereditary).

*Kac canonical decomposition (Dynkin scope).* For hereditary `A` of Dynkin type,
`Rep(Q,d)` is a smooth irreducible affine space and a finite union of orbits
(finite type), so there is a UNIQUE dense (open) orbit — the generic module,
which is rigid. Its decomposition `d = Σ m_i β_i` into positive (real Schur)
roots with `ext(β_i, β_j)=0` for `i≠j` (Kac) IS the canonical decomposition. For
Dynkin every positive root is a real Schur root with a unique indecomposable
`M_β` (Gabriel), so the generic ext `ext(β_i, β_j)` equals the actual
`dim Ext¹_A(M_{β_i}, M_{β_j})` (one module per root), and rigidity of
`G = ⊕ M_{β_i}^{m_i}` is exactly the Kac condition — a computable, self-certifying
arbiter. **Scope: Dynkin only; Euclidean/wild hereditary is deferred (named
successor); non-hereditary is refused.**

*Degeneration = hom order (rep-finite).* `M ≤_deg N :⇔ O_M ⊆ \overline{O_N}`
(`M` is a degeneration of `N`; `M` more special). Zwara 2000 (`≤_deg = ≤_ext`)
and Bongartz 1996 (`≤_deg = ≤_hom` for representation-finite) give
`M ≤_deg N ⇔ dim Hom(M, X) ≥ dim Hom(N, X)` for all indecomposable `X`. The
generic module (open orbit) is the unique MAXIMUM; the semisimple module the
MINIMUM. Convention fixed for this plan: **`Hom` OUT of the class** (`[M,X]`),
degenerate = lower = bigger hom-dims (matches the metaplan card's
`S₁⊕S₂ <_deg P₁`).

---

### Task 1: `invariants/geometry.py` — orbit dimension + Voigt rigidity

**Files:**
- Create: `src/quiverlab/invariants/geometry.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.orbit_dimension(M)`,
  `Algebra.is_rigid(M)`, `Algebra.rigidity_codim(M)` thin delegates, beside
  `positive_roots`, `core/algebra.py:544`)
- Test: `tests/invariants/test_geometry_orbit.py`

**Interfaces:**
- Consumes: `modules/hom.py::end_dim(M)`, `Algebra.ext(M, N, n) -> int`,
  `Module.dimension_vector() -> dict`, `A.quiver.arrows` (dict `label → (s,t)`),
  `A.quiver.vertices`, `Algebra.tits_form(d)` (P38, `d` a vertex-order list),
  `Algebra.is_hereditary()`.
- Produces:
  ```python
  def group_dim(dimvec) -> int              # dim GL(d) = sum_v d_v^2
  def representation_variety_dim(A, dimvec) -> int   # sum_{a:i->j} d_i d_j (ambient Rep(Q,d))
  def orbit_dimension(M) -> int             # dim O_M = group_dim(d) - end_dim(M)
  def is_rigid(M) -> bool                    # Ext^1_A(M, M) == 0 (Voigt: => open orbit)
  def rigidity_codim(M) -> int               # dim Ext^1_A(M, M); = codim O_M iff hereditary, else upper bound
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_geometry_orbit.py
"""Orbit dimension + Voigt rigidity (Plan 49 / C8). Self-cert: the orbit-dim
identity dim O_M = sum d_v^2 - dim End(M); the Voigt codim identity
dim Rep - dim O_M == dim Ext^1(M,M) on hereditary; is_rigid one call.
Literature: every indecomposable over a Dynkin path algebra is rigid (codim 0)."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.invariants.geometry import (group_dim, is_rigid,
                                           orbit_dimension, rigidity_codim,
                                           representation_variety_dim)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@selfcert
def test_orbit_dim_identity():
    A = _kA2()
    for M in (A.projective(1), A.simple(1), A.simple(2)):
        dv = M.dimension_vector()
        assert orbit_dimension(M) == group_dim(dv) - A.hom(M, M)   # end_dim == hom(M,M)


@xeng
def test_voigt_codim_identity_on_hereditary():
    # hereditary => codim O_M in Rep(Q,d) == dim Ext^1(M, M) (Voigt).
    A = linear_path_algebra(3, field=QQ)
    for v in (1, 2, 3):
        for M in (A.simple(v), A.projective(v), A.injective(v)):
            dv = M.dimension_vector()
            codim = representation_variety_dim(A, dv) - orbit_dimension(M)
            assert codim == A.ext(M, M, 1) == rigidity_codim(M)
            # and the P38 tie: dim Rep - dim GL(d) == -<d,d> = -tits_form
            dlist = [dv[w] for w in A.quiver.vertices]
            assert representation_variety_dim(A, dv) - group_dim(dv) == -A.tits_form(dlist)


@lit
def test_every_dynkin_indecomposable_is_rigid():
    A = linear_path_algebra(4, field=QQ)          # kA4, hereditary Dynkin
    ar = A.ar_quiver()
    assert ar.is_complete
    for vtx in ar.vertices:                        # all 10 indecomposables
        M = vtx["module"]
        assert is_rigid(M)                         # real-root modules: Ext^1(M,M)=0
        dv = M.dimension_vector()
        assert representation_variety_dim(A, dv) - orbit_dimension(M) == 0   # codim 0 (open)


@selfcert
def test_non_rigid_has_positive_codim():
    # the 2-Kronecker generic (1,1)-module is NOT rigid: Ext^1(M,M) = 1 (isotropic delta).
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    M = A.module({1: 1, 2: 1}, {"a": [[1]], "b": [[1]]})   # a generic pencil point
    assert is_rigid(M) is False
    assert rigidity_codim(M) == 1


@selfcert
def test_orbit_dim_field_parity():
    # orbit dim / rigidity are field-independent for these examples (QQ vs GF(32003)).
    for field in (QQ, GF(32003)):
        A = linear_path_algebra(3, field=field)
        M = A.projective(1)
        assert orbit_dimension(M) == group_dim(M.dimension_vector()) - A.hom(M, M)
        assert is_rigid(M) is True
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/invariants/test_geometry_orbit.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.invariants.geometry`

- [ ] **Step 3: Implement `src/quiverlab/invariants/geometry.py`**

```python
"""Geometry of representations (Plan 49 / C8): orbit dimensions in the
representation variety, Voigt rigidity, and the Kac canonical decomposition.
Right modules; the orbit/rigidity layer is exact over EVERY Domain
(dim End + Ext are exact); the canonical decomposition is a hereditary-Dynkin
notion with a loud refusal off scope. Float-free (int/dict)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules.hom import end_dim


def group_dim(dimvec):
    """dim GL(d) = sum_v d_v^2 -- the base-change group acting on Rep(Q, d)."""
    return sum(int(n) * int(n) for n in dimvec.values())


def representation_variety_dim(algebra, dimvec):
    """dim Rep(Q, d) = sum over arrows a: i -> j of d_i * d_j (the AMBIENT quiver
    representation variety). For kQ/I the module variety mod_A(d) is the closed
    subvariety cut by the relations; this is its ambient space -- stated honestly
    by every codim consumer."""
    total = 0
    for (s, t) in algebra.quiver.arrows.values():
        total += int(dimvec[s]) * int(dimvec[t])
    return total


def orbit_dimension(M):
    """dim of the GL(d)-orbit of M in Rep(Q, d):
        dim O_M = dim GL(d) - dim_k End_A(M) = sum_v d_v^2 - dim End_A(M).
    Aut(M) is Zariski-open in the affine space End_A(M), so
    dim Stab(M) = dim End_A(M). Holds verbatim for kQ/I."""
    return group_dim(M.dimension_vector()) - end_dim(M)


def is_rigid(M):
    """Voigt: M is rigid iff Ext^1_A(M, M) = 0, and then O_M is OPEN in the
    module variety (rigid => open orbit, in general). One Ext call."""
    return M.algebra.ext(M, M, 1) == 0


def rigidity_codim(M):
    """dim Ext^1_A(M, M). HONEST semantics (stated by every caller):
      * HEREDITARY A: EQUALS codim of the orbit closure in Rep(Q, d) -- Voigt,
        since Rep(Q, d) is smooth: codim O_M = dim End(M) - <d,d> = dim Ext^1(M,M).
      * general kQ/I: only an UPPER BOUND -- codim_{mod_A(d)} O_M <= dim Ext^1(M,M)
        (equality iff mod_A(d) is smooth at M). NEVER claim equality off hereditary.
    """
    return M.algebra.ext(M, M, 1)
```

**Adjust to reality:** `Algebra.ext(M, M, 1)` is the canonical Ext-dim call
(`core/algebra.py:416`); `end_dim(M) == A.hom(M, M)` (both route through
`hom_space`). `representation_variety_dim` reads `A.quiver.arrows.values()`
(each `(s, t)`); confirm the vertex keys of `dimvec` match `A.quiver.vertices`
(they do — `dimension_vector` keys on `A.quiver.vertices`). The `Algebra`
delegates are one-liners (`return orbit_dimension(M)` etc., lazy-importing
`quiverlab.invariants.geometry`).

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/geometry.py src/quiverlab/core/algebra.py tests/invariants/test_geometry_orbit.py
git commit -m "feat(geometry): orbit dimension + Voigt rigidity -- dim O_M = dim GL(d) - dim End(M); codim honest (=Ext^1 hereditary, upper bound on kQ/I)"
```

---

### Task 2: `canonical_decomposition` — the Kac decomposition (Dynkin)

**Files:**
- Modify: `src/quiverlab/invariants/geometry.py`
- Modify: `src/quiverlab/core/algebra.py` (add
  `Algebra.canonical_decomposition(d, budget=4096)` delegate)
- Test: `tests/invariants/test_geometry_canonical.py`

**Interfaces:**
- Consumes: `Algebra.positive_roots()` (P38 — vertex-order tuples, raises off
  hereditary Dynkin), `Algebra.form_type()` (`"finite"|"tame"|"wild"`),
  `Algebra.dynkin_type()`, `Algebra.is_hereditary()`, `Algebra.ar_quiver()`
  (P41 — `.vertices[i]["module"]` / `["dimvec"]`), `Algebra.ext(M, N, 1)`,
  `modules/hom.py::identify_standard`.
- Produces:
  ```python
  def canonical_decomposition(A, d, *, budget=4096) -> list[dict]
      # Kac canonical decomposition of d over a HEREDITARY DYNKIN algebra:
      #   [{"root": (v-order tuple), "multiplicity": m, "name": "P_1"|"S_2"|None}, ...]
      # d = sum m_i * root_i ; ext(root_i, root_j) = 0 for all i,j in the support
      # (i.e. the generic module (+) M_{root_i}^{m_i} is RIGID). Certified per
      # instance (Ext^1(G, G) == 0). Loud refusal off scope (see below).
  ```

**Scope (loud `QuiverlabError` otherwise), the resolved ambiguity:**
1. **non-hereditary** `A` (has relations / gl.dim > 1): the canonical
   decomposition is a hereditary (Kac) notion — refused, pointing at
   `orbit_dimension`/`is_rigid` which DO apply to `kQ/I`.
2. **hereditary, but not Dynkin** (Euclidean/tame or wild): the general
   Schofield generic-ext + Derksen–Weyman recursion (imaginary Schur roots,
   isotropic-root multiplicities) is **DEFERRED to a named successor plan** —
   refused loudly naming it. (`positive_roots()` already refuses off Dynkin, so
   the finite-root machinery is unavailable anyway; the explicit `form_type()`
   check gives the honest message BEFORE `positive_roots` fires.)
3. Dynkin (finite type): full support, rigidity-certified.

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_geometry_canonical.py
"""Kac canonical decomposition over hereditary Dynkin (Plan 49 / C8). Literature:
(2,1) over kA2 = P1 (+) S1 (hand-derived). Self-cert: sum m_i root_i == d, each
root a positive root, pairwise ext 0 (=> generic module rigid). Refusals: the
Euclidean 2-Kronecker (deferred) and the non-hereditary k[x]/(x^2)."""
import pytest

from quiverlab import Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.geometry import canonical_decomposition

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@lit
def test_kA2_21_is_P1_plus_S1():
    A = _kA2()
    cd = canonical_decomposition(A, {1: 2, 2: 1})
    # (2,1) = (1,1) + (1,0) = P1 (+) S1 (the rank-1 generic map k^2 -> k)
    got = sorted((c["root"], c["multiplicity"]) for c in cd)
    assert got == [((1, 0), 1), ((1, 1), 1)]
    names = sorted(c["name"] for c in cd)
    assert names == ["P_1", "S_1"]


@selfcert
@pytest.mark.parametrize("d", [{1: 1, 2: 1}, {1: 2, 2: 1}, {1: 2, 2: 2}, {1: 3, 2: 2}])
def test_kA2_decomposition_is_certified(d):
    A = _kA2()
    cd = canonical_decomposition(A, d)
    verts = list(A.quiver.vertices)
    # sum of components == d
    total = {v: 0 for v in verts}
    for c in cd:
        for k, v in enumerate(verts):
            total[v] += c["multiplicity"] * c["root"][k]
    assert total == d
    # each component is a positive root
    roots = {tuple(r) for r in A.positive_roots()}
    assert all(tuple(c["root"]) in roots for c in cd)


@xeng
def test_canonical_top_equals_generic_over_kA3():
    # cross-engine tie: the canonical decomposition of d is the decomposition of
    # the MAXIMUM of degeneration_order(A, d) (the generic/open-orbit module).
    from quiverlab.modules.degeneration import degeneration_order
    A = linear_path_algebra(3, field=QQ)
    d = {1: 1, 2: 1, 3: 1}
    cd = canonical_decomposition(A, d)
    poset = degeneration_order(A, d)
    top = max(poset.vertices, key=lambda x: x["orbit_dim"])
    cd_multiset = sorted((c["root"], c["multiplicity"]) for c in cd)
    top_multiset = sorted((s_root, m) for (s_root, m) in _summand_roots(top, A))
    assert cd_multiset == top_multiset          # SHARPEN _summand_roots in Step 3


@selfcert
def test_euclidean_is_deferred():
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    with pytest.raises(QuiverlabError, match="Dynkin|finite|deferred|Euclidean"):
        canonical_decomposition(A, {1: 1, 2: 1})     # delta = (1,1) isotropic, tame


@selfcert
def test_non_hereditary_refused():
    A = truncated_polynomial(2, field=QQ)            # k[x]/(x^2): not hereditary
    with pytest.raises(QuiverlabError, match="hereditary"):
        canonical_decomposition(A, {1: 1})
```

(`test_canonical_top_equals_generic_over_kA3` uses Task 3; write it in Task 3's
commit if Task 3 lands after — the tie is the arbiter either way. `_summand_roots`
is a one-line test helper reading a poset vertex's `summands` back to root tuples.)

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** (append to `geometry.py`)

```python
def canonical_decomposition(algebra, d, *, budget=4096):
    from quiverlab.modules.hom import identify_standard
    verts = list(algebra.quiver.vertices)
    dvec = tuple(int(d[v]) for v in verts) if isinstance(d, dict) else tuple(int(x) for x in d)
    if not algebra.is_hereditary():
        raise QuiverlabError(
            "canonical_decomposition: A is not hereditary -- the Kac canonical "
            "decomposition is a hereditary (path-algebra) notion",
            hint="orbit_dimension and is_rigid still apply to kQ/I; the canonical "
                 "decomposition does not")
    if algebra.form_type() != "finite":
        raise QuiverlabError(
            f"canonical_decomposition: A is hereditary of {algebra.dynkin_type()} "
            "type (Euclidean/wild, not Dynkin/finite) -- the general Schofield / "
            "Derksen-Weyman recursion (imaginary Schur roots, isotropic multiplicities) "
            "is DEFERRED to a successor plan",
            hint="Dynkin (finite type) is fully supported here")
    roots = [tuple(int(x) for x in r) for r in algebra.positive_roots()]  # Dynkin-gated
    # one indecomposable per positive root (Gabriel): read them off the AR quiver
    ar = algebra.ar_quiver()
    if not ar.is_complete:
        raise QuiverlabError(f"canonical_decomposition: AR knitting did not close "
                             f"(status {ar.status})")
    mod_of = {}
    for vtx in ar.vertices:
        key = tuple(vtx["dimvec"][v] for v in verts)
        mod_of[key] = vtx["module"]
    missing = [r for r in roots if r not in mod_of]
    if missing:
        raise QuiverlabError(f"canonical_decomposition: no indecomposable found for "
                             f"root(s) {missing} (AR/root mismatch)")
    # generic ext(beta, gamma) = actual Ext^1 between the unique indecomposables
    ext = {(a, b): algebra.ext(mod_of[a], mod_of[b], 1) for a in roots for b in roots}
    best = _search_canonical(roots, dvec, ext, budget)
    if best is None:
        raise QuiverlabError("canonical_decomposition: no rigid decomposition within "
                             "budget", hint="raise `budget`")
    # certificate: the direct sum of the chosen indecomposables is rigid (Kac)
    out = []
    for root, mult in sorted(best):
        std = identify_standard(mod_of[root])
        name = f"{std[0][0].upper()}_{std[1]}" if std else None
        out.append({"root": root, "multiplicity": mult, "name": name})
    return out


def _search_canonical(roots, target, ext, budget):
    """DFS for the (unique, over Dynkin) multiset of roots summing to `target`
    whose direct sum is RIGID: for every pair beta, gamma in the chosen support
    ext(beta, gamma) = ext(gamma, beta) = 0 (real roots have ext(beta,beta)=0).
    Since dim Ext^1(G, G) = sum_{b,g} m_b m_g ext(b, g) >= 0, this is exactly
    Ext^1(G, G) = 0 (the Kac criterion)."""
    n = len(target)
    order = sorted(roots, key=lambda r: -sum(r))       # bigger roots first
    steps = {"n": 0}

    def compatible(a, b):
        return ext[(a, b)] == 0 and ext[(b, a)] == 0

    def dfs(idx, remaining, chosen):
        steps["n"] += 1
        if steps["n"] > budget:
            raise QuiverlabError("canonical_decomposition: search budget exceeded")
        if all(x == 0 for x in remaining):
            return list(chosen.items())
        if idx >= len(order):
            return None
        beta = order[idx]
        if any(beta[k] > remaining[k] for k in range(n)):
            return dfs(idx + 1, remaining, chosen)
        maxm = min(remaining[k] // beta[k] for k in range(n) if beta[k] > 0)
        for m in range(maxm, -1, -1):
            if m > 0:
                if ext[(beta, beta)] != 0:              # imaginary root: skip (Dynkin: never)
                    continue
                if any(not compatible(beta, g) for g in chosen):
                    continue
            newrem = tuple(remaining[k] - m * beta[k] for k in range(n))
            newchosen = dict(chosen)
            if m > 0:
                newchosen[beta] = m
            res = dfs(idx + 1, newrem, newchosen)
            if res is not None:
                return res
        return None

    return dfs(0, target, {})
```

**Adjust to reality (Task 2):**
- **The rigidity certificate is the arbiter, not the search order.** After
  `_search_canonical` returns `best`, ASSERT the certificate by building
  `G = direct_sum` of the chosen indecomposables (with multiplicity) and checking
  `algebra.ext(G, G, 1) == 0` — a wrong pick fails loudly. Equivalently (cheaper,
  used inside the search) `Σ m_b m_g ext[(b,g)] == 0`; keep the direct-sum
  double-check in the DYNKIN test as the independent oracle.
- **`ar_quiver` gives one module per root** for Dynkin (rep-finite, not
  self-injective — closes). If `dynkin_type()` is E₆/₇/₈ the root set is larger
  but still finite; the `budget` guards the multiset search, not the roots.
- **Uniqueness** (Dynkin): the rigid decomposition is unique. The DFS returns the
  first; a paranoid variant can continue and raise if a second rigid solution
  appears (there is none over Dynkin) — optional, gated behind a debug flag.
- **`positive_roots()` returns real roots for Dynkin** and `ext[(β,β)] == 0` for
  all of them (real Schur roots); the `ext[(β,β)] != 0` guard is dead code on
  Dynkin but documents the imaginary-root boundary the deferred Euclidean case
  crosses (kept for the successor).
- The `Algebra.canonical_decomposition(d, budget=4096)` delegate lazy-imports
  `quiverlab.invariants.geometry`.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/geometry.py src/quiverlab/core/algebra.py tests/invariants/test_geometry_canonical.py
git commit -m "feat(geometry): Kac canonical decomposition (Dynkin) -- root enumeration + rigidity certificate; Euclidean/wild deferred, non-hereditary refused"
```

---

### Task 3: `modules/degeneration.py` — the degeneration / hom order poset

**Files:**
- Create: `src/quiverlab/modules/degeneration.py`
- Modify: `src/quiverlab/core/algebra.py` (add
  `Algebra.degeneration_order(d, budget=256)` delegate)
- Test: `tests/modules/test_degeneration.py`

**Interfaces:**
- Consumes: `Algebra.ar_quiver(budget_modules=...)` (P41 — the finite
  indecomposable universe; `.is_complete`/`.status`; each vertex's `["module"]`
  and `["dimvec"]`), `modules/hom.py::hom_dim`/`identify_standard`,
  `modules/morphism.py::direct_sum`, `invariants/geometry.py::group_dim` (for the
  per-class orbit dim).
- Produces:
  ```python
  @dataclass
  class DegenerationPoset:                 # mirrors P41's ARQuiver shape + contract
      vertices    # list of dicts: {"index", "dimvec", "summands": [(name, mult)],
                  #                  "orbit_dim": int, "module": Module, "is_generic": bool}
      covers      # list of (lower_i, upper_j) Hasse cover pairs (i degenerates from j)
      is_complete # True iff rep-finite (AR knit closed); False iff budget/unsupported
      status      # "complete" | "budget" | "unsupported" | "error"
      note        # optional string
  def degeneration_order(A, d, *, budget=256) -> DegenerationPoset
      # For representation-FINITE A: all iso-classes of dimension vector d, ordered
      # by the hom (= degeneration) order. M <=_deg N iff dim Hom(M,X) >= dim Hom(N,X)
      # for all indecomposable X (Bongartz 1996, Zwara 2000). STOPS LOUDLY at the
      # budget / on self-injective (knit unsupported) -- never a silent partial poset.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_degeneration.py
"""Degeneration = hom order for representation-finite algebras (Plan 49 / C8).
Literature (hand-derived): kA2 d=(1,1) is the 2-chain S1(+)S2 <_deg P1; kA3
d=(1,1,1) is the diamond (semisimple bottom, two incomparable middles, uniserial
top). Self-cert: the relation is a partial order; orbit dim strictly increases up
each cover; the unique maximum is the generic (rigid) module."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.geometry import is_rigid, orbit_dimension
from quiverlab.modules.degeneration import degeneration_order

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@lit
def test_kA2_11_is_a_two_chain():
    A = _kA2()
    P = degeneration_order(A, {1: 1, 2: 1})
    assert P.is_complete and len(P.vertices) == 2      # {P1}, {S1 (+) S2}
    assert len(P.covers) == 1                          # a single cover
    lo, hi = P.covers[0]
    # the lower (degenerate) class is S1 (+) S2 (two summands); the upper is P1
    assert len(P.vertices[lo]["summands"]) == 2
    assert len(P.vertices[hi]["summands"]) == 1
    assert P.vertices[hi]["orbit_dim"] > P.vertices[lo]["orbit_dim"]
    assert P.vertices[hi]["is_generic"] and is_rigid(P.vertices[hi]["module"])


@lit
def test_kA3_111_is_the_diamond():
    A = linear_path_algebra(3, field=QQ)
    P = degeneration_order(A, {1: 1, 2: 1, 3: 1})
    assert P.is_complete and len(P.vertices) == 4       # [1,3]; [1,2](+)S3; S1(+)[2,3]; S1(+)S2(+)S3
    # diamond: 4 covers (bottom->2 middles, 2 middles->top), 2 incomparable middles
    assert len(P.covers) == 4
    orbit = sorted(v["orbit_dim"] for v in P.vertices)
    assert orbit == [0, 1, 1, 2]                        # Z=0, A=B=1, G=2 (open)
    tops = [v for v in P.vertices if v["is_generic"]]
    assert len(tops) == 1 and tops[0]["orbit_dim"] == 2


@selfcert
def test_hom_order_is_a_partial_order():
    A = linear_path_algebra(3, field=QQ)
    P = degeneration_order(A, {1: 1, 2: 1, 3: 1})
    leq = _reachability(P)                              # test-local: reflexive-transitive closure
    n = len(P.vertices)
    for a in range(n):
        assert leq[a][a]                               # reflexive
    for a in range(n):
        for b in range(n):
            if a != b and leq[a][b]:
                assert not leq[b][a]                   # antisymmetric


@selfcert
def test_orbit_dim_monotone_up_covers():
    A = linear_path_algebra(3, field=QQ)
    P = degeneration_order(A, {1: 1, 2: 1, 3: 1})
    for lo, hi in P.covers:
        assert P.vertices[lo]["orbit_dim"] < P.vertices[hi]["orbit_dim"]


@selfcert
def test_representation_infinite_refused():
    # 2-Kronecker is representation-infinite: knit cannot close -> loud, not a poset.
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    P = degeneration_order(A, {1: 2, 2: 2}, budget=40)
    assert P.is_complete is False and P.status in ("budget", "error", "unsupported")
```

(`_reachability` is a one-line test helper: BFS the reflexive-transitive closure
of `covers`. The `orbit == [0,1,1,2]` and `covers == 4` pins are the HAND-DERIVED
diamond — see the plan's foundation section.)

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement `src/quiverlab/modules/degeneration.py`**

```python
"""Degeneration (closure) order for representation-finite algebras (Plan 49 / C8),
computed as the equivalent HOM order (Bongartz 1996 Adv. Math. 121; Zwara 2000
Compos. Math. 121):
  M <=_deg N  <=>  O_M subset closure(O_N)  <=>  dim Hom(M, X) >= dim Hom(N, X)
                   for every indecomposable X.
The generic module (open orbit) is the MAXIMUM; the semisimple module the MINIMUM.
Convention: Hom OUT of the class ([M,X]); degenerate = lower = bigger hom-dims.
Rep-finite only -- the finite indecomposable universe is P41 knit_ar_quiver; a
rep-infinite or self-injective input refuses loudly (never a partial poset)."""
from __future__ import annotations

from dataclasses import dataclass, field

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.geometry import group_dim
from quiverlab.modules.hom import hom_dim, identify_standard


@dataclass
class DegenerationPoset:
    vertices: list
    covers: list
    is_complete: bool
    status: str
    note: str = ""


def _std_name(std):
    return f"{std[0][0].upper()}_{std[1]}" if std else None


def _enumerate_classes(idimvecs, target, budget):
    """All multisets (multiplicity vectors over the indecomposables) whose
    dim vectors sum to `target`. Budget-capped; returns None on overflow."""
    n = len(target)
    m = len(idimvecs)
    out = []

    def dfs(idx, remaining, chosen):
        if len(out) > budget:
            return False
        if all(x == 0 for x in remaining):
            out.append(tuple(chosen))
            return True
        if idx >= m:
            return True
        dv = idimvecs[idx]
        if any(dv[k] > remaining[k] for k in range(n)):
            return dfs(idx + 1, remaining, chosen)
        maxc = min(remaining[k] // dv[k] for k in range(n) if dv[k] > 0)
        for c in range(maxc, -1, -1):
            newrem = tuple(remaining[k] - c * dv[k] for k in range(n))
            chosen[idx] = c
            if not dfs(idx + 1, newrem, chosen):
                chosen[idx] = 0
                return False
            chosen[idx] = 0
        return True

    ok = dfs(0, target, [0] * m)
    return out if (ok and len(out) <= budget) else None


def _hasse_covers(leq):
    """Cover pairs (i, j) with i < j (i.e. leq[i][j], i != j) and NO k strictly
    between (i < k < j)."""
    n = len(leq)
    covers = []
    for i in range(n):
        for j in range(n):
            if i == j or not leq[i][j]:
                continue
            if any(k not in (i, j) and leq[i][k] and leq[k][j] for k in range(n)):
                continue
            covers.append((i, j))
    return covers


def degeneration_order(algebra, d, *, budget=256):
    from quiverlab.modules.morphism import direct_sum
    verts = list(algebra.quiver.vertices)
    dvec = {v: int(d[v]) for v in verts} if isinstance(d, dict) else \
        {v: int(d[i]) for i, v in enumerate(verts)}
    target = tuple(dvec[v] for v in verts)
    ar = algebra.ar_quiver(budget_modules=budget)
    if not ar.is_complete:
        return DegenerationPoset([], [], is_complete=False, status=ar.status,
                                 note="not representation-finite / knit did not close")
    indecs = [vtx["module"] for vtx in ar.vertices]
    idimvecs = [tuple(vtx["dimvec"][v] for v in verts) for vtx in ar.vertices]
    classes = _enumerate_classes(idimvecs, target, budget)
    if classes is None:
        return DegenerationPoset([], [], is_complete=False, status="budget",
                                 note="too many iso-classes of dimension d")
    # pairwise Hom matrix over the indecomposables (computed once)
    H = [[hom_dim(indecs[i], indecs[j]) for j in range(len(indecs))]
         for i in range(len(indecs))]

    def homvec(a):
        return tuple(sum(a[i] * H[i][j] for i in range(len(a)))
                     for j in range(len(indecs)))

    homvecs = [homvec(a) for a in classes]
    if len(set(homvecs)) != len(homvecs):
        raise QuiverlabError("degeneration_order: two iso-classes share a hom-vector "
                             "(rep-finite => impossible; report the input)")
    n = len(classes)
    leq = [[all(homvecs[a][j] >= homvecs[b][j] for j in range(len(indecs)))
            for b in range(n)] for a in range(n)]       # a <=_deg b: a has bigger homs
    covers = _hasse_covers(leq)
    gd = group_dim(dvec)
    max_orbit = None
    vertices = []
    for a in range(n):
        end_a = sum(classes[a][i] * classes[a][j] * H[i][j]
                    for i in range(len(indecs)) for j in range(len(indecs)))
        orbit = gd - end_a
        max_orbit = orbit if max_orbit is None else max(max_orbit, orbit)
        summ = [(indecs[i], classes[a][i]) for i in range(len(indecs)) if classes[a][i]]
        mods = []
        for (m, mult) in summ:
            mods.extend([m] * mult)
        module = direct_sum(*mods)[0] if len(mods) > 1 else mods[0]
        vertices.append({
            "index": a,
            "dimvec": dict(dvec),
            "summands": [(_std_name(identify_standard(m)), mult) for m, mult in summ],
            "orbit_dim": orbit,
            "module": module,
        })
    for v in vertices:
        v["is_generic"] = (v["orbit_dim"] == max_orbit)
    return DegenerationPoset(vertices, covers, is_complete=True, status="complete")
```

**Adjust to reality (Task 3):**
- **`end_a` via the Hom matrix** is `dim End(⊕ M_i^{a_i}) = Σ a_i a_j [M_i, M_j]`,
  so `orbit = group_dim(d) − end_a` matches `invariants.geometry.orbit_dimension`
  on the built `module` — add a per-vertex `assert orbit == orbit_dimension(module)`
  in the DYNKIN self-cert test as the independent oracle (do NOT assert inside the
  loop for every input — it doubles the End solve).
- **Antisymmetry rests on distinct hom-vectors.** For rep-finite algebras the
  hom-vector `([M, X])_X` determines `M` up to iso (Auslander), so distinct
  classes have distinct hom-vectors — the loud `QuiverlabError` fires only on a
  genuine bug or an input the AR knit mis-reported.
- **`covers` orientation** is `(lower, upper)` = `(i, j)` with `i <_deg j`
  (`i` more degenerate). The Hasse renderer (Task 4) draws `lower` at the bottom.
- **Self-injective** (e.g. a cyclic Nakayama) makes `ar_quiver` return
  `status="unsupported"` (P41) — surfaced as `is_complete=False` with that status,
  never a partial poset. State this on the verification page as an honest-scope
  boundary alongside "rep-finite only".
- **The unique maximum ties Task 2:** for hereditary Dynkin the `is_generic`
  vertex's `summands` are exactly `canonical_decomposition(A, d)`
  (`test_canonical_top_equals_generic_over_kA3`, `oracle_crossengine`).
- The `Algebra.degeneration_order(d, budget=256)` delegate lazy-imports
  `quiverlab.modules.degeneration`.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/degeneration.py src/quiverlab/core/algebra.py tests/modules/test_degeneration.py
git commit -m "feat(degeneration): hom-order (=degeneration) poset for rep-finite algebras -- kA2 2-chain, kA3 diamond, honest rep-finite cap"
```

---

### Task 4: viz twin renderers — the Hasse poset (TikZ + HTML/SVG)

**Files:**
- Modify: `src/quiverlab/viz/layout.py` (add `poset_layout(nodes, covers)`)
- Modify: `src/quiverlab/viz/tikz.py` (add `tikz_hasse(poset)`)
- Create: `src/quiverlab/viz/hasse_html.py` (`hasse_svg(poset)`)
- Test: `tests/modules/test_degeneration_render.py`

**Interfaces:**
- Consumes: `DegenerationPoset` (Task 3; `.vertices`, `.covers`), the exact
  `fractions.Fraction` layout idiom already in `viz/layout.py`.
- Produces:
  ```python
  def poset_layout(nodes, covers) -> dict   # {index: (x: Fraction, y: int)}
      # y = rank (longest upward chain from a minimum); x spreads nodes per rank,
      # exact Fraction centred. Reusable for P45's exchange lattice.
  def tikz_hasse(poset, label=lambda v: ...) -> str    # a standalone tikzpicture
  def hasse_svg(poset) -> str                          # an inline <svg> string (GUI/report)
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_degeneration_render.py
"""The Hasse twin renderers (Plan 49 / C8). Self-cert: the layout ranks the poset
(minimum at rank 0, cover endpoints one rank apart); tikz_hasse / hasse_svg emit
non-empty markup naming every class and drawing every cover. Float-free."""
import pytest

from fractions import Fraction

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.degeneration import degeneration_order
from quiverlab.viz.hasse_html import hasse_svg
from quiverlab.viz.layout import poset_layout
from quiverlab.viz.tikz import tikz_hasse

pytestmark = pytest.mark.oracle_selfcert


def _diamond():
    return degeneration_order(linear_path_algebra(3, field=QQ), {1: 1, 2: 1, 3: 1})


def test_layout_ranks_the_poset():
    P = _diamond()
    nodes = [v["index"] for v in P.vertices]
    pos = poset_layout(nodes, P.covers)
    assert all(isinstance(pos[i][0], (int, Fraction)) for i in nodes)   # exact x
    for lo, hi in P.covers:
        assert pos[hi][1] == pos[lo][1] + 1        # a cover spans exactly one rank
    ys = {pos[i][1] for i in nodes}
    assert min(ys) == 0 and max(ys) == 2           # diamond has 3 ranks


def test_tikz_and_svg_are_nonempty_and_name_classes():
    P = _diamond()
    tk = tikz_hasse(P)
    svg = hasse_svg(P)
    assert tk.strip().startswith(r"\begin{tikzpicture}") and r"\end{tikzpicture}" in tk
    assert svg.strip().startswith("<svg") and "</svg>" in svg
    # every cover is drawn (4 edges in the diamond)
    assert tk.count("--") >= len(P.covers)
    assert svg.count("<line") >= len(P.covers)
```

- [ ] **Step 2: Run to verify failure** — `ImportError: poset_layout` / `hasse_svg`

- [ ] **Step 3: Implement** — `poset_layout`: rank each node by the longest chain
  of covers from a minimum (a topological longest-path on the cover DAG, exact
  ints); within a rank, spread nodes at exact `Fraction` x-positions centred on
  0. `tikz_hasse`: emit `\begin{tikzpicture}` … `\node`/`\draw` per vertex/cover
  (mirror `tikz_quiver`'s `_coord`, `tikz.py:10`). `hasse_svg`: emit an
  `<svg>` with `<line>` covers + `<text>`/`<circle>` nodes (fraction→int pixel
  via a fixed exact scale, done in Python since the string is static markup, or —
  preferred — leave fractions as the coordinate and multiply by an integer scale
  so no float appears). The label defaults to the class's `summands`
  (`"P_1"`, `"S_1 ⊕ S_2"`).

**Adjust to reality (Task 4):**
- **No float ever** — `poset_layout` x-positions are `Fraction`; `hasse_svg`
  multiplies by an integer pixel scale (`Fraction * int` stays exact, then
  `int(...)` for the attribute is an int cast of an exact rational, allowed —
  `math.floor`/`int` on an exact value is not a float literal, but confirm the
  no-floats gate: cast via `int(round(...))` is fine only if `round` gets an int;
  keep the multiply-then-`.numerator//.denominator` idiom to stay provably int).
- **Reuse, not fork** — `poset_layout` is deliberately generic (nodes + cover
  pairs) so P45's exchange-graph Hasse can render through the same twin; note this
  in the docstring.
- **This is a library renderer, not a GUI compute kind.** `degeneration_order`
  is not wired to a live `orbit_geometry` panel (that kind reports the canonical
  decomposition, not the whole poset); the Hasse twin is available for reports and
  for any future degeneration compute kind. State this scope honestly (Task 7).

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/viz/layout.py src/quiverlab/viz/tikz.py src/quiverlab/viz/hasse_html.py tests/modules/test_degeneration_render.py
git commit -m "feat(viz): generic Hasse poset layout + tikz_hasse + hasse_svg twin renderers (exact Fraction, reusable by P45)"
```

---

### Task 5: GUI story — the `orbit_geometry` module compute kind

**Files:**
- Modify: `src/quiverlab/invariants/geometry.py` (add the shared builder
  `orbit_geometry_block(M)`)
- Modify: `src/quiverlab/hpc/spec.py` (`MODULE_KINDS`, `_MOD_REFS`,
  `_dispatch_module`)
- Modify: `docs/gui/runner.py` (the twin: `_MODULE_KINDS`, `_MOD_REFS`,
  `_module_block`)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (pick-list entry +
  `renderBlock` branch + `renderOrbitGeometry`)
- Modify: `webapp/static/app.js` (block renderer)
- Modify: `webapp/server/i18n/en.json`, `es.json` (labels)
- Modify: `src/quiverlab/trace/results_html.py` (`_HEADINGS` + `_block_html` branch)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE new fixture, existing byte-identical)
- Test: `tests/webapp/test_orbit_geometry_p49.py`, `tests/gui/test_orbit_geometry_runner_twin.py`

**Interfaces (the `tau` module kind is the template — follow its schema-2 gating
exactly, `spec.py:1602`; `almost_split` is the closest recent add):**
- `orbit_geometry` is a **MODULE kind**: add it to `MODULE_KINDS` (`spec.py:70`),
  so a request naming it **requires schema 2** with a `module` block (the same
  guard as `tau`, `spec.py:630-632`). It is NOT a range kind, NOT a trace kind.
  `_MOD_REFS["orbit_geometry"] = ["voigt_rigidity", "kac_canonical",
  "schofield_general_reps", "derksen_weyman_canonical"]`.
- **Shared block builder** (byte-identical across both runners — both import the
  ONE library function):
  ```python
  def orbit_geometry_block(M):
      A = M.algebra
      dv = M.dimension_vector()
      ext1 = A.ext(M, M, 1)
      hereditary = A.is_hereditary()
      block = {
          "kind": "orbit_geometry", "side": M.side,
          "dim_vector": {str(v): int(n) for v, n in dv.items()},
          "group_dim": group_dim(dv),                       # dim GL(d)
          "rep_variety_dim": representation_variety_dim(A, dv),  # dim Rep(Q,d) (ambient)
          "end_dim": end_dim(M),                            # dim End_A(M)
          "orbit_dim": orbit_dimension(M),                  # dim O_M
          "rigid": ext1 == 0,
          "ext1_self": ext1,                                # dim Ext^1(M,M) = rigidity_codim
          "hereditary": hereditary,
          "codim_semantics": "hereditary" if hereditary else "general",
          "latex": r"\dim \mathcal{O}_M = \sum_v d_v^2 - \dim_k \operatorname{End}_A(M)",
      }
      try:                                                  # Dynkin-only extra
          block["canonical_decomposition"] = canonical_decomposition(A, dv)
      except QuiverlabError as e:
          block["canonical_decomposition"] = None
          block["canonical_note"] = str(e)
      return block
  ```
- `_dispatch_module` branch: `block = _with_refs(orbit_geometry_block(M), "orbit_geometry")`.
  **No hard refusal** — orbit dim / rigidity / codim compute for ANY module over
  ANY `kQ/I`; the canonical decomposition is the conditional extra (present on
  hereditary Dynkin, else `null` + `canonical_note`). If `end_dim`/`ext` raise
  (a char-scope edge), catch into the honest `{"kind": "orbit_geometry",
  "error": <loud message>}` per the Plan-30 per-entry precedent — never a 500.
- GUI: an `orbit_geometry` pick-list entry (alongside `tau`/`almost_split`); the
  renderer prints the orbit dim / `dim GL(d)` / `dim Rep(Q,d)` / `dim End` line,
  the **rigidity verdict** with the honest codim gloss (`hereditary`: "= codim of
  the orbit closure (Voigt)"; `general`: "≤ codim (upper bound; the module
  variety of kQ/I is cut by relations)"), and — when present — the canonical
  decomposition `d = ⊕ β_i^{m_i}` with names; MathJax + citations like the
  existing module blocks.

- [ ] **Step 1: Write the failing cross-runner test** (unmarked, extras-gated dir
  — copy `tests/webapp/test_almost_split_p41.py` / the m0729 runner-pair fixture):

```python
# tests/webapp/test_orbit_geometry_p49.py
"""The orbit_geometry module kind: schema-2 gated, served by hpc.spec, mirrored by
the Pyodide twin. Reports orbit dim + rigidity + honest codim; canonical
decomposition present on hereditary Dynkin."""


def test_orbit_geometry_block_shape(tmp_path):
    # schema-2 request: kA3 over QQ, module = builtin simple S_2, compute ["orbit_geometry"].
    # Assert: block["orbit_dim"] is an int; block["rigid"] is True (S_2 a real-root brick);
    #   block["ext1_self"] == 0; block["hereditary"] is True;
    #   block["canonical_decomposition"] is a non-empty list (dim (0,1,0) = S_2 itself);
    #   "voigt_rigidity" in [k for k, _ in block["citations"]].
    ...


def test_orbit_geometry_codim_gloss_general(tmp_path):
    # a NON-hereditary algebra (k[x]/(x^2), module S) => block["hereditary"] False,
    #   block["codim_semantics"] == "general", block["canonical_decomposition"] is None,
    #   block["canonical_note"] mentions "hereditary".
    ...


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py; assert json.dumps(sort_keys=True)
    # equality on the orbit_geometry block (both runners byte-identical via the shared
    # orbit_geometry_block builder).
    ...
```

(Fill the `...` by copying the almost_split fixture; the assertions listed are the
contract. `orbit_geometry` needs `schema: 2` + a `module` block — follow the tau
kind's request shape verbatim.)

- [ ] **Step 2: Implement** the shared builder + the `spec.py` branch + the
  `docs/gui/runner.py` twin (both call `orbit_geometry_block`, keep the wrapping
  shape-identical), the two `gui.js` renderers (`renderOrbitGeometry(div, b)` +
  the `renderBlock` branch, `gui.js:2284`) + `app.js`
  (`orbitGeometryBlock(block, d)` + the `renderResult` dispatch, `app.js:567`),
  i18n keys (`inv.orbit_geometry`, `block.orbit_geometry.title`,
  `block.orbit_geometry.orbit`, `block.orbit_geometry.rigid`,
  `block.orbit_geometry.codim_hered`, `block.orbit_geometry.codim_general`,
  `block.orbit_geometry.canonical` — EN and ES), and the `results_html.py`
  `_HEADINGS["orbit_geometry"] = "Orbit geometry"` + `_block_html` branch
  (a small facts table + the canonical-decomposition row; reuse `_math`).

- [ ] **Step 3: Add ONE golden fixture** (`orbit_geometry_kA3_s2`) to
  `_runner_goldens.json`; note it in the `test_runner_delegation.py` docstring
  change-log (the `almost_split_a3_s2` ADD entry is the model). Verify existing
  goldens stay byte-identical BEFORE adding.

- [ ] **Step 4: Run the gates**

Run: `... -m pytest tests/webapp/test_orbit_geometry_p49.py tests/webapp/test_runner_delegation.py tests/gui/test_orbit_geometry_runner_twin.py tests/hpc -q`
Expected: PASS (both runners byte-identical via the shared builder)

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): orbit_geometry module kind -- orbit dim + Voigt rigidity + honest codim + Dynkin canonical decomposition, one golden"
```

---

### Task 6: QPA cross-oracle (probe-first) + honest fallback

**Files:**
- Modify: `src/quiverlab/qpa/scripts.py` / `crosscheck.py` (only if a live verb exists)
- Test: `tests/qpa/test_geometry_qpa.py`

**Interfaces:**
- Consumes: the QPA session (`session.should_skip_qpa()`, `session.libgap_handle()`)
  — mirror the Plan-35 products probe (`tests/qpa/test_products_qpa.py`).
- **The metaplan P49 card names "QPA degeneration" as an oracle, but QPA's
  degeneration/orbit surface is UNVERIFIED here — probe live first.** QPA has
  `DTr`/`AlmostSplitSequence`/module homs (used by P41), but no confirmed
  `orbit dimension` / `canonical decomposition` / `degeneration order` verb.

- [ ] **Step 1: Probe live QPA** for any geometry surface (the Plan-35
  fail-if-appears pattern — the primary oracles are the hand-derived posets + the
  self-cert identities + the Voigt/Kac literature pins; QPA is opportunistic):

```python
# tests/qpa/test_geometry_qpa.py
"""QPA probe for a representation-geometry surface (Plan 49). qpa-marked: skips
locally, mandatory under QUIVERLAB_REQUIRE_QPA=1. QPA has NO confirmed orbit /
canonical-decomposition / degeneration verb -- this test documents that honestly
and FAILS if one ever appears (forcing a real crosscheck)."""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_has_no_degeneration_surface():
    lg = session.libgap_handle()
    # candidate names a QPA geometry verb might use; none are expected present.
    for name in ("OrbitDimension", "CanonicalDecomposition",
                 "DegenerationOrder", "DegenerateModule"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), (
            f"QPA now exposes {name} -- wire a real crosscheck (this assert is the "
            "trip-wire that Plan 49's QPA scope note is stale)")


def test_end_dim_matches_qpa_hom():
    # what QPA CAN corroborate: dim End(M) (=> orbit dim via dim GL(d) - dim End).
    # cross-check dim Hom_A(M, M) against QPA's HomOverAlgebra on kA3 simples.
    from quiverlab import linear_path_algebra
    from quiverlab.fields import QQ
    A = linear_path_algebra(3, field=QQ)
    for v in (1, 2, 3):
        M = A.simple(v)
        A.crosscheck("hom_glue", M, M).assert_agree()  # existing P37 hom crosscheck verb
```

**Adjust to reality (Task 6):** the crosscheck dispatch (`qpa/crosscheck.py`) has
NO `"hom"` verb — the module-Hom crosscheck is `"hom_glue"` (`crosscheck.py:542`).
Confirm its argument signature before wiring; if `HomOverAlgebra` cannot be driven
as a two-module crosscheck cleanly, the probe test (`test_qpa_has_no_degeneration_surface`)
is the deliverable on its own and the `dim End` corroboration falls back to the
`tests/qpa/test_ar_qpa.py` (P41) module-Hom coverage — honest, never a silent skip.

- [ ] **Step 2:** If (and only if) the probe finds a live geometry verb, add a
  `crosscheck_orbit`/`crosscheck_degeneration` in `crosscheck.py` mirroring
  `crosscheck_tau`. Otherwise the probe + the `end_dim`-via-`HomOverAlgebra`
  corroboration IS the QPA leg (orbit dimension is `dim GL(d) − dim End`, and
  `dim End` is QPA-verifiable), and the verification page records "QPA has no
  degeneration/orbit surface — the orbit-dim's `dim End` factor is QPA-checked,
  the rest is theory/self-cert" honestly.

- [ ] **Step 3: Run live** `... -m pytest tests/qpa/test_geometry_qpa.py -v`
  (the venv has `[qpa]`). Expected: PASS live (probe absent; `HomOverAlgebra` agrees).

- [ ] **Step 4: Commit**

```bash
git add src/quiverlab/qpa/ tests/qpa/test_geometry_qpa.py
git commit -m "test(qpa): geometry probe -- QPA has no orbit/canonical/degeneration verb (fail-if-appears); dim End corroborated via HomOverAlgebra"
```

---

### Task 7: verification page, citations, README, metaplan, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-05-metaplan-v0.2.0.md` (tick the P49 card)
- Test: existing release gates (`tests/release/test_oracle_classes.py`,
  `tests/citations/`)

- [ ] **Step 1: Citations** (VERIFIED BibTeX only — the `_r(...)` registry
  precedent is `_r(key, bibtex_key, kind, title, annotation, *tags)`,
  `registry.py:24`; `bibtex()` hard-fails if `.bib` and registry disagree). Add
  to `references.bib`:

```bibtex
@article{Kac1980,
  author  = {Kac, Victor G.},
  title   = {Infinite root systems, representations of graphs and invariant theory},
  journal = {Inventiones Mathematicae},
  volume  = {56}, number = {1}, pages = {57--92}, year = {1980},
}
@article{Schofield1992,
  author  = {Schofield, Aidan},
  title   = {General representations of quivers},
  journal = {Proceedings of the London Mathematical Society. Third Series},
  volume  = {65}, number = {1}, pages = {46--64}, year = {1992},
}
@article{DerksenWeyman2002,
  author  = {Derksen, Harm and Weyman, Jerzy},
  title   = {On the canonical decomposition of quiver representations},
  journal = {Compositio Mathematica},
  volume  = {133}, number = {3}, pages = {245--265}, year = {2002},
}
@article{Zwara2000,
  author  = {Zwara, Grzegorz},
  title   = {Degenerations of finite-dimensional modules are given by extensions},
  journal = {Compositio Mathematica},
  volume  = {121}, number = {2}, pages = {205--218}, year = {2000},
}
@article{Bongartz1996,
  author  = {Bongartz, Klaus},
  title   = {On degenerations and extensions of finite dimensional modules},
  journal = {Advances in Mathematics},
  volume  = {121}, number = {2}, pages = {245--287}, year = {1996},
}
@book{Voigt1977,
  author    = {Voigt, Detlef},
  title     = {Induzierte Darstellungen in der Theorie der endlichen, algebraischen Gruppen},
  series    = {Lecture Notes in Mathematics},
  volume    = {592}, publisher = {Springer}, year = {1977},
}
```

and in `registry.py` (mirroring `_r("assem_book", "ASS2006", "foundation", ...)`):

```python
_r("kac_canonical", "Kac1980", "foundation",
   "Infinite root systems, representations of graphs and invariant theory",
   "Kac's roots, Schur roots, and the canonical decomposition of a dimension "
   "vector -- the ground truth for Plan 49's canonical_decomposition.", "geometry"),
_r("schofield_general_reps", "Schofield1992", "foundation",
   "General representations of quivers",
   "Schofield's generic hom/ext and the general-representation identities "
   "hom - ext = <a,b> underlying the canonical decomposition.", "geometry"),
_r("derksen_weyman_canonical", "DerksenWeyman2002", "foundation",
   "On the canonical decomposition of quiver representations",
   "The Derksen-Weyman recursive algorithm for the canonical decomposition "
   "(Plan 49 ships the Dynkin case; Euclidean/wild is the named deferral).", "geometry"),
_r("zwara_degenerations", "Zwara2000", "foundation",
   "Degenerations of finite-dimensional modules are given by extensions",
   "Zwara: the degeneration order equals the extension order for Artin algebras "
   "-- half of Plan 49's degeneration_order theorem.", "geometry"),
_r("bongartz_degenerations", "Bongartz1996", "foundation",
   "On degenerations and extensions of finite dimensional modules",
   "Bongartz: degeneration = hom order for representation-finite algebras -- the "
   "computable form Plan 49's degeneration_order uses.", "geometry"),
_r("voigt_rigidity", "Voigt1977", "foundation",
   "Induzierte Darstellungen ... (Voigt's lemma)",
   "Voigt's lemma: Ext^1(M,M) = 0 => the orbit of M is open (rigid => open orbit); "
   "the codim = dim Ext^1(M,M) equality on hereditary algebras.", "geometry"),
```

  **Spec-ambiguity resolution (recorded):** BibTeX **entry keys** are
  `Kac1980`/`Schofield1992`/`DerksenWeyman2002`/`Zwara2000`/`Bongartz1996`/
  `Voigt1977`; the **registry/citation keys** the block references are snake-case
  `kac_canonical`/`schofield_general_reps`/`derksen_weyman_canonical`/
  `zwara_degenerations`/`bongartz_degenerations`/`voigt_rigidity` (the house
  convention, `assem_book → ASS2006`). All six are BibTeX-verified before they
  ship. The `kind` field is `"foundation"` (matching AIR/ARS); the free tag is
  `"geometry"`. Crawley-Boevey's *Lectures on representations of quivers* is the
  natural accessible secondary but is unpublished (no stable BibTeX) — added ONLY
  if a worker BibTeX-verifies a citable version at merge; the six above cover the
  shipped surface.

- [ ] **Step 2: Verification page.** Add the Plan-49 subsystem rows:
  - `invariants/geometry.py` (orbit/rigidity) — `oracle_selfcert` (orbit-dim
    identity, field parity); `oracle_crossengine` (Voigt codim identity
    `dim Rep − dim O_M ≡ dim Ext¹(M,M)` on hereditary + the P38 `tits_form` tie);
    `oracle_literature` (every Dynkin indecomposable rigid, codim 0).
  - `invariants/geometry.py` (canonical decomposition) — `oracle_literature`
    (`(2,1) = P₁⊕S₁` over kA₂); `oracle_selfcert` (sum-of-roots + rigidity
    certificate); `oracle_crossengine` (components ≡ `positive_roots`; canonical
    ≡ degeneration maximum).
  - `modules/degeneration.py` — `oracle_literature` (kA₂ 2-chain, kA₃ diamond,
    hand-derived); `oracle_selfcert` (partial-order axioms, orbit-dim
    monotonicity up covers); the **honest semi-decision entry** (complete iff
    rep-finite; loud `status` on the 2-Kronecker / self-injective — metaplan §6).
  - `viz` Hasse twin — `oracle_selfcert` (layout ranks, renderers emit every
    class/cover).
  Add the **honest-scope entries**: (a) `rigidity_codim` is the orbit-closure
  codim only on hereditary; on `kQ/I` it is an UPPER BOUND (the module variety is
  cut by relations) — stated, never claimed as equality; (b) `canonical_decomposition`
  is **Dynkin only** — Euclidean/wild hereditary is DEFERRED (Schofield /
  Derksen–Weyman general recursion, the named successor), non-hereditary refused;
  (c) `degeneration_order` is **representation-finite only** (self-injective
  refused via the P41 knit `unsupported` status); (d) **Hall numbers are out of
  scope** (the P3 axis); (e) **QPA cannot compare** the geometry verbs — QPA has
  no orbit/canonical/degeneration surface (the probe is a fail-if-appears
  trip-wire; the orbit-dim's `dim End` factor IS QPA-checked via
  `HomOverAlgebra`). Recount the class table (`tests/release/test_oracle_classes.py`
  drives the numbers — run collection, paste the LIVE counts, re-run to green;
  mid-merge-train honest).

- [ ] **Step 3: README.** One features line: "geometry of representations
  (Kac/Voigt): orbit dimensions in the representation variety, Voigt rigidity
  with honest codimension, the Kac canonical decomposition of a dimension vector
  (Dynkin), and the Zwara–Bongartz degeneration/hom-order poset for
  representation-finite algebras — the C8 axis."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/invariants/test_geometry*.py -q` (fast),
  `... -m pytest tests/modules/test_degeneration*.py -q` (deep),
  `... -m pytest -q -m fast`, `... -m pytest tests/webapp tests/gui -q`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q` — all green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-49 geometry oracle rows + Kac/Schofield/DW/Zwara/Bongartz/Voigt citations + honest scope (Dynkin canonical, rep-finite degeneration, no-QPA) + recounted classes"
```

---

## Acceptance (Plan-49 definition of done)

1. `orbit_dimension`, `group_dim`, `representation_variety_dim`, `is_rigid`,
   `rigidity_codim`, `canonical_decomposition`, `orbit_geometry_block` public in
   `invariants/geometry.py`; `degeneration_order` + `DegenerationPoset` in
   `modules/degeneration.py`; the `Algebra` delegates
   (`orbit_dimension`/`is_rigid`/`rigidity_codim`/`canonical_decomposition`/
   `degeneration_order`) named; every result self-certified or loudly refusing.
2. `dim O_M = Σ d_v² − dim End_A(M)` verified; the Voigt codim identity
   `dim Rep(Q,d) − dim O_M ≡ dim Ext¹(M,M)` pinned cross-engine on hereditary,
   with the P38 `tits_form` tie; `rigidity_codim` documented as the orbit-closure
   codim on hereditary and an UPPER BOUND on `kQ/I` (never equality off
   hereditary); every Dynkin indecomposable rigid (codim 0).
3. The Kac canonical decomposition is Dynkin-scoped and rigidity-certified:
   `(2,1) = P₁ ⊕ S₁` over kA₂ pinned; components are positive roots and pairwise
   ext-orthogonal; Euclidean/wild hereditary DEFERRED (named successor),
   non-hereditary refused — all loud.
4. `degeneration_order` reproduces the hand-derived posets (kA₂ (1,1) = 2-chain
   `S₁⊕S₂ <_deg P₁`; kA₃ (1,1,1) = diamond, orbit dims `[0,1,1,2]`, 4 covers),
   is a certified partial order with orbit-dim strictly monotone up covers, its
   unique maximum ≡ the canonical decomposition (cross-engine), and refuses
   **loudly** with `status ∈ {"budget","unsupported","error"}` on the 2-Kronecker
   / self-injective — the honest semi-decision contract.
5. The Hasse twin renderers (`poset_layout`/`tikz_hasse`/`hasse_svg`) rank and
   draw the poset, exact `Fraction`, float-free, reusable by P45.
6. `orbit_geometry` clickable end-to-end (GUI canvas → block → report) in EN+ES,
   schema-2 gated like `tau`, both runners byte-identical via the shared
   `orbit_geometry_block`, ONE golden added with a documented change-log entry;
   the honest codim gloss (hereditary vs general) and the Dynkin canonical
   decomposition rendered.
7. QPA probe green (`-m qpa`): no orbit/canonical/degeneration verb (fail-if-
   appears), the orbit-dim's `dim End` factor corroborated via `HomOverAlgebra`.
8. `docs/verification.md` recounted (live numbers, mid-merge-train honest); the
   six citations added and BibTeX-verified; README line + metaplan P49 card
   ticked; fast (`tests/invariants/test_geometry*.py`) + deep
   (`tests/modules/test_degeneration*.py`) + webapp/gui + qpa + release +
   citations suites green. Honest scope recorded: Voigt codim equality on
   hereditary only; Dynkin-only canonical decomposition (Euclidean/wild deferred);
   rep-finite-only degeneration (self-injective refused); Hall numbers out; QPA
   cannot compare the geometry verbs.
