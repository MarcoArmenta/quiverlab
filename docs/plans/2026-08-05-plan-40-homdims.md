# Plan 40: C6 Homological-Dimensions Family — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The homological-dimension surface every representation theorist reaches
for, exact and honest: public **syzygy/cosyzygy** operators, the **Igusa–Todorov
functions φ/ψ** (finite K₀ + syzygy bookkeeping), **finitistic-dimension bounds**,
**dominant dimension**, **Gorenstein dimension** + `is_gorenstein`, and **Ω/τ-
periodicity certificates** — all composed on the existing minimal-resolution /
Ext / Plan-30-decompose engines, each result carrying the `GlobalDimension`-style
certified-value-or-lower-bound honesty (never a bare number when the engine did
not resolve). Delooping level (Gélinas) is scoped with an explicit decision gate
and deferred with a recorded reason (Task F). One no-code `homological_profile`
compute kind exposes the whole family end-to-end.

**Architecture:** A new `src/quiverlab/modules/homdims.py` holds the Igusa–Todorov
functions, dominant/Gorenstein dimensions, the periodicity certificates, and the
finitistic bounds, mirroring `modules/ext.py::global_dimension`'s home (library
function + `core/algebra.py` thin wrapper); the honesty dataclasses live beside
`GlobalDimension`'s pattern. `syzygy`/`cosyzygy` are extracted as **public**
functions in `modules/resolution.py` (the natural home — the syzygy is already
computed inline in `minimal_resolution`'s loop) and are byte-stable refactors
gated by the existing engine/resolution suites. Everything is thin exact linear
algebra over the EXISTING primitives: `resolution.py::projective_cover`/
`minimal_resolution`, `duality.py::dualize`/`tau`/`tau_minus`,
`injective.py::injective_resolution`/`injective_dimension`,
`decompose.py::decompose`, `hom.py::is_isomorphic`, `ext.py::GlobalDimension`.
The K₀ bookkeeping (φ/ψ) is over **ℤ**, independent of the algebra's field —
integer ranks via sympy, exact. No new math engines.

**Tech Stack:** pure exact linear algebra over `Domain` (`modules/linalg_mod`)
for the module-level operations; **sympy integer/rational matrix rank** for the
Igusa–Todorov K₀ (a ℤ-module computation, not a `Domain` computation); Plan-37
`direct_sum` for the regular module. No floats in `src/`.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P37 must be merged before this plan runs** — the regular module is built with
  `modules/morphism.py::direct_sum`, and the deferred delooping task (F) would
  consume `is_direct_summand`. Branch `plan-40-homdims` off `dev` AFTER P37 lands.
- Buckets are auto-assigned by directory (`tests/conftest.py`): `tests/modules/`
  and `tests/families/` → **deep**; `tests/invariants/` → **fast**; `tests/qpa/`
  → **qpa**. Run new deep tests by path during development, finish with a
  `-m deep` spot-run of the touched files.
- **The decompose char-caveat is load-bearing here.** `decompose`/`is_isomorphic`
  RAISE loudly over GF(p) when `char ≤ dim M` (the trace-form radical is
  unreliable). The Igusa–Todorov batteries — which lean on `decompose` +
  `is_isomorphic` for the K₀ classes — therefore run over **QQ** (char 0) or a
  **large prime GF(32003)** so both routines decide. Any φ/ψ test over a small
  GF(p) with `char ≤ dim` is a bug in the test, not the engine. The library
  functions inherit `decompose`'s loud refusal — never a silent wrong φ/ψ.
- **The honesty pattern is mandatory** (the `GlobalDimension` precedent,
  `modules/ext.py:63-98`): a homological dimension that the bounded engine did not
  resolve is returned as a *labeled certified lower bound* (or an explicit
  `infinite`/`undecided` marker), NEVER a bare int and NEVER a fabricated `∞`.
  The syzygy/injective engines only ever certify "resolved within depth N" or
  "not resolved within depth N"; they never prove infinity by termination, so
  `is_gorenstein` is three-valued **True / None** (a `False` verdict would require
  a proof of infinite injective dimension the bounded engine does not furnish —
  Task E's periodicity certificate is the only route to such a proof, and is not
  wired into `is_gorenstein` here).
- Plan-32 markers: the φ=ψ=pd identity, periodicity `is_isomorphic` certificates,
  and dimension-additivity identities are `oracle_selfcert`; the Barrios–Mata
  truncated-path φ/ψ closed forms, hereditary/self-injective domdim & Gorenstein
  values, and cyclic-Nakayama period-from-Kupisch pins are `oracle_literature`;
  QPA comparisons live in `tests/qpa/` (bucket = the class, never double-marked).
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted
  class table, `tests/release/test_oracle_classes.py` green) and the README line.
- Conventional commits; green at every commit.

---

### Task A: public `syzygy` / `cosyzygy`

**Files:**
- Modify: `src/quiverlab/modules/resolution.py` (extract `syzygy`, add `cosyzygy`)
- Modify: `src/quiverlab/modules/module.py` (add `Module.syzygy()`, `Module.cosyzygy()`)
- Test: `tests/modules/test_syzygy.py`

**Interfaces:**
- Consumes: `resolution.py::projective_cover(M) -> (Q0, d0, gens)`,
  `linalg_mod::kernel_columns`, `radtopsoc::submodule`, `duality.py::dualize`.
- Produces:
  ```python
  def syzygy(M) -> "Module"        # Omega M = ker(projective_cover(M) -> M),
                                   # a submodule of Q0; minimality => no projective
                                   # summands. Omega of a projective is the zero module.
  def cosyzygy(M) -> "Module"      # Omega^{-1} M via the dual route: D(syzygy(D M)).
  ```
  `Module.syzygy(self)` / `Module.cosyzygy(self)` delegate.
  `minimal_resolution`'s inline syzygy step is refactored to call `syzygy`
  (byte-stable: the existing engine/resolution suites gate identical
  `(terms, dmats)`).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_syzygy.py
"""Public syzygy / cosyzygy. Self-certifying: the syzygy is the kernel of the
projective cover (the exact submodule the minimal resolution already builds)."""
import pytest

from quiverlab import GF, Quiver, truncated_polynomial
from quiverlab.modules.hom import is_isomorphic
from quiverlab.modules.resolution import cosyzygy, syzygy

pytestmark = pytest.mark.oracle_selfcert


def _a3rel():
    # kA3 with a*b = 0 (Quiver 1->2->3, radical-square-zero on the 1->3 path).
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_syzygy_of_simple_is_radical_of_cover():
    A = _a3rel()
    S1 = A.simple(1)
    Om = syzygy(S1)                          # Omega S1 = rad P1 = S2 here
    assert Om.dimension_vector() == A.simple(2).dimension_vector()
    assert is_isomorphic(Om, A.simple(2))


def test_syzygy_of_projective_is_zero():
    A = _a3rel()
    assert syzygy(A.projective(2)).dim == 0


def test_truncated_polynomial_syzygy_shift():
    # k[x]/(x^3): Omega(k[x]/(x^i)) = k[x]/(x^{3-i}); Omega S = M_2 (dim 2).
    A = truncated_polynomial(3, field=GF(7))
    S = A.simple(1)
    Om = syzygy(S)
    assert Om.dim == 2
    assert syzygy(Om).dim == 1               # Omega^2 S = M_1 = S again (dim 1)


def test_cosyzygy_matches_dual_route():
    # k[x]/(x^3) is self-injective: cosyzygy(S) = Omega^{-1} S = M_2 (dim 2).
    A = truncated_polynomial(3, field=GF(7))
    S = A.simple(1)
    co = cosyzygy(S)
    assert co.dim == 2
    from quiverlab.modules.duality import dualize
    assert co.dimension_vector() == dualize(syzygy(dualize(S))).dimension_vector()


def test_module_method_delegates():
    A = _a3rel()
    assert A.simple(1).syzygy().dim == syzygy(A.simple(1)).dim
    assert A.simple(1).cosyzygy().dim == cosyzygy(A.simple(1)).dim
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_syzygy.py -v`
Expected: FAIL — `ImportError: cannot import name 'syzygy'`

- [ ] **Step 3: Implement** — extract the syzygy step in `resolution.py`:

```python
def syzygy(M):
    """Omega M = ker(projective_cover(M) -> M) as a submodule of Q0. Minimality of
    the cover => Omega M has no projective summands. Omega of a projective (or the
    zero module) is the zero module."""
    Q0, d0, _ = projective_cover(M)
    if not (d0 and d0[0]):                   # cover is an iso / M projective / zero
        from quiverlab.modules.module import Module
        dom = M.domain
        return Module(M.algebra, 0, {lab: lm.zeros(0, 0, dom) for lab in M.action},
                      name=f"Omega({M.name})", side=M.side)
    ker_cols = lm.kernel_columns(d0, M.domain)
    if not ker_cols:
        from quiverlab.modules.module import Module
        dom = M.domain
        return Module(M.algebra, 0, {lab: lm.zeros(0, 0, dom) for lab in M.action},
                      name=f"Omega({M.name})", side=M.side)
    return submodule(Q0, ker_cols, name=f"Omega({M.name})")


def cosyzygy(M):
    """Omega^{-1} M via the dual route: D(Omega(D M)). Cosyzygy of an injective
    (or the zero module) is the zero module (dualizes Omega of a projective)."""
    from quiverlab.modules.duality import dualize
    out = dualize(syzygy(dualize(M)))
    out.name = f"Omega^-1({M.name})"
    return out
```

Refactor `minimal_resolution`'s `n == 1` syzygy so the first `Omega_1` is
literally `syzygy(M)` (identical columns — `projective_cover` + `kernel_columns`
+ `submodule` is exactly what the loop already does at `resolution.py:107-115`).
The deeper `Omega_n` steps stay as-is (they syzygy the *previous syzygy* through
the already-computed cover); DO NOT rewrite them to re-call `projective_cover` if
that changes the byte output — the acceptance gate is the untouched
`tests/engine/` + `tests/modules/test_*resolution*` suites staying green. Add
`Module.syzygy` / `Module.cosyzygy` in `module.py` delegating to these.

**Adjust to reality:** read `resolution.py:100-130` and mirror its EXACT zero-
module construction (it builds `Module(M.algebra, 0, {lab: zeros...}, ...)`); the
zero-syzygy branch above copies that idiom — verify the action-label set matches.

- [ ] **Step 4: Run tests, verify pass; run the byte-stability neighbors**

Run: `... -m pytest tests/modules/test_syzygy.py tests/modules/test_resolution.py tests/engine -q`
Expected: PASS — resolutions byte-unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/resolution.py src/quiverlab/modules/module.py tests/modules/test_syzygy.py
git commit -m "feat(modules): public syzygy/cosyzygy extracted from minimal_resolution (byte-stable refactor)"
```

---

### Task B: Igusa–Todorov functions φ / ψ

**Files:**
- Create: `src/quiverlab/modules/homdims.py`
- Modify: `src/quiverlab/modules/module.py` (add `Module.igusa_todorov_phi()` / `_psi()`)
- Test: `tests/modules/test_igusa_todorov.py`

**Interfaces:**
- Consumes: `syzygy` (Task A), `decompose.py::decompose(M, budget)`,
  `hom.py::is_isomorphic`, `Module.projective_resolution(bound).pd()`, sympy
  (`Matrix.rank` over ℤ — the K₀ is an abelian-group rank, NOT a `Domain` rank).
- Produces:
  ```python
  def igusa_todorov_phi(M, budget=512, bound=64) -> int
  def igusa_todorov_psi(M, budget=512, bound=64) -> int
  ```
  Both terminate (Fitting's lemma); both RAISE the `decompose`/`is_isomorphic`
  loud refusal unchanged when the K₀ classes cannot be certified (char caveat).
  `Module.igusa_todorov_phi(self, budget=512)` / `_psi` delegate.

  **Definitions pinned verbatim (Igusa–Todorov, Fields Inst. Commun. 45 (2005),
  201–204).** Let `K₀` be the free abelian group on iso-classes of indecomposable
  NON-projective f.d. modules ([projective] = 0). `⟨M⟩` ⊆ K₀ is the subgroup
  generated by the classes of the non-projective indecomposable summands of `M`;
  its rank is the number of distinct such summands. `L: K₀ → K₀` is the syzygy
  operator, `L[X] = [ΩX]` (ΩX decomposed, projective summands dropped). The rank
  sequence `rank(Lⁿ⟨M⟩)` is monotone non-increasing (L maps each subgroup onto
  the next), so it stabilizes.
  - **φ(M) = the least n ≥ 0 with `rank(Lⁿ⟨M⟩) = rank(Lⁿ⁺¹⟨M⟩)`** (Fitting: once
    two consecutive ranks agree, L is injective on the free group `Lⁿ⟨M⟩` and the
    rank is fixed forever after).
  - **ψ(M) = φ(M) + fpd(Ω^{φ(M)} M)**, where `fpd(N) = max{ pd(Y) : Y a summand of
    N with pd(Y) < ∞ }` (0 if there are none).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_igusa_todorov.py
"""Igusa-Todorov phi/psi. Literature closed forms (Barrios-Mata 2019 for
truncated path algebras) + the phi=psi=pd identity for finite proj. dimension
(Igusa-Todorov 2005) + additivity. Over QQ/GF(32003) -- the decompose char
caveat forbids char <= dim M."""
import pytest

from quiverlab import GF, TruncatedPathAlgebra, linear_path_algebra, \
    truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import igusa_todorov_phi, igusa_todorov_psi

pytestmark = pytest.mark.oracle_literature


def test_phi_equals_pd_when_finite():
    # kA3/rad^2 (a*b=0): Nakayama, gl.dim 2. pd S1 = 2, pd S2 = 1, pd S3 = 0.
    A = TruncatedPathAlgebra("A3", 2, field=QQ)
    for v, expect in ((1, 2), (2, 1), (3, 0)):
        S = A.simple(v)
        pd = S.projective_resolution(8).pd()
        assert pd == expect
        assert igusa_todorov_phi(S) == pd              # phi = pd for finite pd
        assert igusa_todorov_psi(S) == pd              # psi = pd for finite pd


def test_self_injective_truncated_all_zero():
    # k[x]/(x^a) is self-injective: every module has phi = psi = 0
    # (Omega permutes {M_1..M_{a-1}}, so ranks are stable at n=0; no finite-pd
    # non-projective summand, so fpd = 0). Barrios-Mata 2019.
    A = truncated_polynomial(4, field=QQ)
    for i in (1, 2, 3):
        M = A.simple(1)
        for _ in range(i - 1):
            M = M.syzygy()
        assert igusa_todorov_phi(M) == 0
        assert igusa_todorov_psi(M) == 0


def test_psi_geq_phi_and_projective_additivity():
    A = linear_path_algebra(3, field=QQ)               # hereditary kA3
    from quiverlab.modules.morphism import direct_sum
    S1, P3 = A.simple(1), A.projective(3)
    assert igusa_todorov_psi(S1) >= igusa_todorov_phi(S1)
    D, _, _ = direct_sum(S1, P3)                        # add a projective summand
    assert igusa_todorov_phi(D) == igusa_todorov_phi(S1)   # phi(M (+) P) = phi(M)


@pytest.mark.oracle_selfcert
def test_char_caveat_refuses_loudly():
    # over GF(2) with dim M >= 2 the decompose trace form is unreliable ->
    # the K0 bookkeeping must inherit the loud refusal, never a silent phi.
    A = truncated_polynomial(3, field=GF(2))
    from quiverlab.errors import QuiverlabError
    M = A.simple(1).syzygy()                            # dim 2 over GF(2): char <= dim
    with pytest.raises(QuiverlabError):
        igusa_todorov_phi(M)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.modules.homdims`

- [ ] **Step 3: Implement `src/quiverlab/modules/homdims.py`** (the K₀ core):

```python
"""Homological-dimension family (Plan 40 / C6): Igusa-Todorov phi/psi, finitistic
bounds, dominant + Gorenstein dimension, Omega/tau-periodicity certificates.

Composed on the existing minimal-resolution / Ext / decompose engines. Every
dimension follows the GlobalDimension honesty rule: a certified value, or a
labeled lower bound / infinite / undecided marker -- never a bare or fabricated
number. The Igusa-Todorov K0 is an ABELIAN-GROUP (Z) computation: integer ranks
via sympy, independent of the algebra's field."""
from __future__ import annotations

from quiverlab.modules.decompose import decompose
from quiverlab.modules.hom import is_isomorphic
from quiverlab.modules.resolution import syzygy


def _is_projective(X, bound=1):
    return X.dim == 0 or X.projective_resolution(bound).pd() == 0


class _K0:
    """Registry of non-projective indecomposable iso-classes discovered on the fly;
    maps a module to its multiplicity vector over the registry (a Z-vector)."""

    def __init__(self, budget):
        self.reps = []
        self.budget = budget

    def _index(self, X):
        for i, R in enumerate(self.reps):
            if R.dim == X.dim and is_isomorphic(R, X):    # loud if undecidable
                return i
        self.reps.append(X)
        return len(self.reps) - 1

    def classvec(self, M):
        vec = {}
        for S, mult in decompose(M, budget=self.budget):  # loud on char caveat
            if _is_projective(S):
                continue
            k = self._index(S)
            vec[k] = vec.get(k, 0) + mult
        return vec

    def rank(self, modules):
        import sympy
        vecs = [self.classvec(X) for X in modules if X.dim]
        width = len(self.reps)
        if not vecs or width == 0:
            return 0
        rows = [[v.get(j, 0) for j in range(width)] for v in vecs]
        return sympy.Matrix(rows).rank()


def igusa_todorov_phi(M, budget=512, bound=64):
    K = _K0(budget)
    seq = [S for S, _ in decompose(M, budget=budget) if not _is_projective(S)]
    prev = K.rank(seq)
    for n in range(1, bound + 1):
        seq = [syzygy(X) for X in seq]
        seq = [X for X in seq if X.dim]
        cur = K.rank(seq)
        if cur == prev:
            return n - 1                        # least n with rank L^n = rank L^{n+1}
        prev = cur
    raise QuiverlabError(                        # Fitting guarantees stabilization
        f"igusa_todorov_phi: rank did not stabilize within depth {bound}",
        hint="raise `bound`; a genuine non-stabilization indicates a K0 bug")


def igusa_todorov_psi(M, budget=512, bound=64):
    phi = igusa_todorov_phi(M, budget=budget, bound=bound)
    W = M
    for _ in range(phi):
        W = syzygy(W)
    fpd = 0
    for Y, _ in decompose(W, budget=budget):
        pdY = Y.projective_resolution(bound).pd()
        if pdY is not None:
            fpd = max(fpd, pdY)
    return phi + fpd
```

(Add `from quiverlab.errors import QuiverlabError` at the top.) Wire
`Module.igusa_todorov_phi` / `_psi` in `module.py`.

**Adjust to reality:** confirm `sympy.Matrix([]).rank()` / an all-zero matrix
returns `0` (the `not vecs or width == 0` guard covers the empty case); confirm
`decompose` groups by iso-class so `classvec`'s multiplicities are correct.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/homdims.py src/quiverlab/modules/module.py tests/modules/test_igusa_todorov.py
git commit -m "feat(modules): Igusa-Todorov phi/psi on the finite K0 + syzygy bookkeeping (char-caveat-honest)"
```

---

### Task C: dominant dimension

**Files:**
- Modify: `src/quiverlab/modules/homdims.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.dominant_dimension`)
- Test: `tests/modules/test_dominant_dimension.py`

**Interfaces:**
- Consumes: `modules/morphism.py::direct_sum` (regular module = ⊕_v P_v),
  `injective.py::injective_resolution(regular, bound)` (its `.terms[n]` are the
  E^n modules), `Module.projective_resolution(1).pd()` (term projective ⇔ pd 0),
  `ext.py::is_selfinjective` (the ∞ short-circuit).
- Produces:
  ```python
  @dataclass
  class DominantDimension:
      value: int | None      # count of leading projective injective terms
      exact: bool            # value is the true domdim (a non-proj term was reached)
      infinite: bool         # self-injective => domdim = infinity (certified)
      # __int__, __eq__(int) mirroring GlobalDimension; __repr__:
      #   infinite -> "dom.dim = ∞ (self-injective; every injective coresolvent
      #                of A is projective)"
      #   exact    -> "dom.dim = {value}"
      #   else     -> ">= {value} (certified lower bound; not resolved within depth N)"
  def dominant_dimension(A, bound=32) -> DominantDimension
  ```
  domdim A = the index of the first term in a minimal injective coresolution of
  the regular right module `A_A = ⊕_v P_v` that is NOT projective (∞ iff every
  such term is projective ⇔ A self-injective).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_dominant_dimension.py
"""Dominant dimension: leading projective terms of the injective coresolution of
the regular module. Self-injective => infinity; hereditary kA2 => 1."""
import pytest

from quiverlab import Quiver, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import dominant_dimension

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_self_injective_is_infinite():
    A = truncated_polynomial(3, field=QQ)              # k[x]/(x^3): self-injective
    dd = dominant_dimension(A)
    assert dd.infinite is True
    assert dd.value is None


def test_hereditary_kA2_domdim_one():
    # A = P1 (+) P2; E(A) = P1 (+) P1 (projective, counts), E^1 = S1 (not
    # projective, stops) => domdim = 1.
    A = _kA2()
    dd = dominant_dimension(A)
    assert dd.value == 1 and dd.exact is True and dd.infinite is False
    assert int(dd) == 1 and dd == 1


def test_wrapper_matches():
    A = _kA2()
    assert A.dominant_dimension() == dominant_dimension(A)
```

- [ ] **Step 2: Run to verify failure** — `ImportError: dominant_dimension`

- [ ] **Step 3: Implement** (append to `homdims.py`; wrapper in `core/algebra.py`):

```python
def _regular_module(A, side="right"):
    from quiverlab.modules.morphism import direct_sum
    Ps = [A.projective(v, side=side) for v in A.quiver.vertices]
    D, _, _ = direct_sum(*Ps)
    return D


def dominant_dimension(A, bound=32):
    from quiverlab.modules.ext import is_selfinjective
    from quiverlab.modules.injective import injective_resolution
    if is_selfinjective(A):
        return DominantDimension(None, exact=False, infinite=True)
    reg = _regular_module(A)
    res = injective_resolution(reg, bound)
    count = 0
    for n in range(bound + 1):
        term = res.terms[n] if n < len(res.terms) else None
        if term is None or term.dim == 0:                 # coresolution ended proj.
            return DominantDimension(count, exact=True, infinite=False)
        if term.projective_resolution(1).pd() == 0:        # E^n projective: counts
            count += 1
        else:
            return DominantDimension(count, exact=True, infinite=False)
    return DominantDimension(count, exact=False, infinite=False)   # lower bound
```

Define `DominantDimension` beside `GlobalDimension`'s honesty pattern (store
`_bound` for the repr message if you want the depth named; keep the three-way
`__repr__` above). `Algebra.dominant_dimension(self, bound=32)` imports and
delegates (mirror `Algebra.global_dimension` at `core/algebra.py:432`), citing
`assem_book` + the dominant-dimension reference (see the citation note in the
final task).

**Adjust to reality:** confirm `InjectiveResolution.terms[n]` is the E^n *module*
(it is — `injective.py:64-71`), not the vertex list (`.term(n)` is the vertex
list). The regular module must accept `injective_resolution` (it dualizes then
resolves over A^op — a plain `Module`, works).

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/homdims.py src/quiverlab/core/algebra.py tests/modules/test_dominant_dimension.py
git commit -m "feat(modules): dominant_dimension -- leading projective injective coresolvents, self-injective => infinity"
```

---

### Task D: Gorenstein dimension + `is_gorenstein`

**Files:**
- Modify: `src/quiverlab/modules/homdims.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.gorenstein_dimension`, `Algebra.is_gorenstein`)
- Test: `tests/modules/test_gorenstein.py`

**Interfaces:**
- Consumes: `injective.py::injective_dimension(regular, bound)`, `_regular_module`
  (Task C), `Algebra.opposite()` (left side = right over A^op).
- Produces:
  ```python
  @dataclass
  class GorensteinDimension:
      right_id: int | None      # inj.dim of A_A within bound (None = unresolved)
      left_id: int | None       # inj.dim of _A A = inj.dim of the regular A^op-module
      is_gorenstein: bool | None
      # __repr__: both finite -> "Gorenstein: right inj.dim {r}, left inj.dim {l}";
      #           either None  -> "undecided within depth N (an injective dimension
      #                            did not resolve; infinity is not proven)"
  def gorenstein_dimension(A, bound=32) -> GorensteinDimension
  def is_gorenstein(A, bound=32) -> bool | None
  ```
  A is (Iwanaga-)Gorenstein iff both the right and left injective dimensions of
  the regular module are finite. `is_gorenstein` is **three-valued True/None**:
  True when both resolve finite; None (with a loud repr) when either is a mere
  lower bound — the bounded engine never proves infinity, so a bare `False` is
  never emitted (see the Global-Constraints honesty note).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_gorenstein.py
"""Gorenstein dimension via the injective dimension of the regular module on both
sides. Self-injective => 0; hereditary => gl.dim; both cases Gorenstein."""
import pytest

from quiverlab import Quiver, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import gorenstein_dimension, is_gorenstein

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_self_injective_gorenstein_zero():
    A = truncated_polynomial(3, field=QQ)              # symmetric => id(A) = 0
    gd = gorenstein_dimension(A)
    assert gd.right_id == 0 and gd.left_id == 0
    assert gd.is_gorenstein is True
    assert is_gorenstein(A) is True


def test_hereditary_gorenstein_equals_gldim():
    A = _kA2()                                          # gl.dim 1
    gd = gorenstein_dimension(A)
    assert gd.right_id == 1 and gd.left_id == 1
    assert gd.is_gorenstein is True
    assert int(A.global_dimension()) == 1


def test_wrapper_matches():
    A = _kA2()
    assert A.is_gorenstein() is True
```

- [ ] **Step 2: Run to verify failure** — `ImportError: gorenstein_dimension`

- [ ] **Step 3: Implement** (append to `homdims.py`; two wrappers in `core/algebra.py`):

```python
def _regular_injective_dimension(A, bound):
    from quiverlab.modules.injective import injective_dimension
    return injective_dimension(_regular_module(A), bound=bound)


def gorenstein_dimension(A, bound=32):
    right = _regular_injective_dimension(A, bound)
    left = _regular_injective_dimension(A.opposite(), bound)
    both_finite = right is not None and left is not None
    verdict = True if both_finite else None            # never a bare False
    return GorensteinDimension(right, left, verdict)


def is_gorenstein(A, bound=32):
    return gorenstein_dimension(A, bound).is_gorenstein
```

`Algebra.gorenstein_dimension` / `Algebra.is_gorenstein` delegate (mirror
`global_dimension`). Citations: `assem_book` + the Gorenstein reference.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/homdims.py src/quiverlab/core/algebra.py tests/modules/test_gorenstein.py
git commit -m "feat(modules): gorenstein_dimension + is_gorenstein (three-valued True/None, honesty-gated)"
```

---

### Task E: Ω- and τ-periodicity certificates

**Files:**
- Modify: `src/quiverlab/modules/homdims.py`
- Modify: `src/quiverlab/modules/module.py` (add `Module.omega_periodicity()` / `Module.tau_periodicity()`)
- Test: `tests/modules/test_periodicity.py`

**Interfaces:**
- Consumes: `syzygy` (Task A), `Module.tau()` (`duality.py`),
  `hom.py::is_isomorphic` (the certificate — loud when undecidable).
- Produces:
  ```python
  def omega_periodicity(M, max_period=12, bound=64) -> int | None
      # least k in 1..max_period with Omega^k M ~ M (is_isomorphic certified);
      # None when no such k within max_period (incl. when some Omega^i M = 0,
      # i.e. M has finite projective dimension => not Omega-periodic).
  def tau_periodicity(M, max_period=12) -> int | None
      # least k in 1..max_period with tau^k M ~ M; None otherwise (a tau^i M = 0
      # => tau-orbit terminates => not periodic).
  ```
  Both RAISE the `is_isomorphic` loud refusal unchanged when an iso comparison is
  undecidable (never a silent None where the answer is unknown).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_periodicity.py
"""Omega/tau-periodicity certificates (is_isomorphic-certified). Cyclic Nakayama
kZ_n/J^2 has Omega-periodic simples of period n (from the Kupisch series);
hereditary => not Omega-periodic (finite pd)."""
import pytest

from quiverlab import NakayamaAlgebra, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.homdims import omega_periodicity, tau_periodicity

pytestmark = pytest.mark.oracle_literature


def _cyclic_rad2(n):
    # kZ_n / J^2: Kupisch [2]*n, self-injective; Omega(S_i) = S_{i-1}, period n.
    return NakayamaAlgebra(n=n, l=2, cyclic=True, field=QQ)


def test_cyclic_nakayama_omega_period_is_n():
    A = _cyclic_rad2(3)
    S = A.simple(list(A.quiver.vertices)[0])
    p = omega_periodicity(S)
    assert p == 3
    # certificate self-check: Omega^p S ~ S
    Om = S
    for _ in range(p):
        Om = Om.syzygy()
    from quiverlab.modules.hom import is_isomorphic
    assert is_isomorphic(Om, S)


def test_hereditary_not_omega_periodic():
    A = linear_path_algebra(3, field=QQ)               # hereditary: pd finite
    assert omega_periodicity(A.simple(1)) is None


def test_self_injective_tau_periodic():
    A = _cyclic_rad2(3)
    S = A.simple(list(A.quiver.vertices)[0])
    p = tau_periodicity(S)
    assert p is not None and 1 <= p <= 3
    T = S
    for _ in range(p):
        T = T.tau()
    from quiverlab.modules.hom import is_isomorphic
    assert is_isomorphic(T, S)


@pytest.mark.oracle_selfcert
def test_projective_not_omega_periodic():
    A = linear_path_algebra(3, field=QQ)
    assert omega_periodicity(A.projective(1)) is None  # Omega P = 0
```

- [ ] **Step 2: Run to verify failure** — `ImportError: omega_periodicity`

- [ ] **Step 3: Implement** (append to `homdims.py`):

```python
def omega_periodicity(M, max_period=12, bound=64):
    if M.dim == 0 or _is_projective(M):
        return None
    cur = M
    for k in range(1, max_period + 1):
        cur = syzygy(cur)
        if cur.dim == 0:                     # finite projective dimension
            return None
        if cur.dim == M.dim and is_isomorphic(cur, M):   # loud if undecidable
            return k
    return None


def tau_periodicity(M, max_period=12):
    if M.dim == 0 or _is_projective(M):
        return None
    cur = M
    for k in range(1, max_period + 1):
        cur = cur.tau()
        if cur.dim == 0:                     # tau-orbit terminated (hit a projective)
            return None
        if cur.dim == M.dim and is_isomorphic(cur, M):
            return k
    return None
```

Wire `Module.omega_periodicity` / `Module.tau_periodicity` in `module.py`.

**Adjust to reality:** `tau(projective) = 0` (`duality.py:106`), so the `cur.dim
== 0` guard catches the terminating orbit; confirm `is_isomorphic`'s dim-vector
prefilter makes the `cur.dim == M.dim` guard redundant-but-cheap (it short-
circuits before building Hom).

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/homdims.py src/quiverlab/modules/module.py tests/modules/test_periodicity.py
git commit -m "feat(modules): Omega/tau-periodicity certificates (is_isomorphic-certified, loud when undecidable)"
```

---

### Task F: delooping level — DECISION GATE (deferred with recorded reason)

**Files:**
- Modify: `docs/plans/2026-08-05-plan-40-homdims.md` (this doc — the deferral note below)
- Modify (final task): `docs/verification.md` (honest-scope entry)

**Decision (made at authoring time; re-evaluate only if the criterion below is
met during implementation):** The delooping level of Gélinas
(J. Pure Appl. Algebra 2022) is `dell(M) = inf{ n ≥ 0 : Ωⁿ M is a direct summand
of Ω^{n+1} N for SOME f.d. module N }`. The existential quantifier over `N` has
no crisp bounded decision procedure from `syzygy` + `is_direct_summand` alone: it
requires either the cosyzygy-of-the-injectives characterization (a construction
this plan does not build) or an a-priori bound from the representation dimension.
Shipping a heuristic that silently fixes a finite candidate set for `N` would
violate the house honesty rule (it could return a wrong finite `dell` when the
true witness lies outside the probed set). Per Design-decision-7's decision gate,
**Task F is DEFERRED**, recorded here and on the verification page, with:
- the one implementable special case NOTED for a successor plan: when
  `findim(A) < ∞`, Gélinas gives `dell(A) ≤ findim(A) + 1` and the two agree in
  the Iwanaga-Gorenstein case — a bound, not the value;
- the successor entry point: a future C6-extension plan builds the injective-side
  cosyzygy tower and the summand-membership check against it, at which point
  `delooping_level_bound(A, probe_depth)` becomes crisply implementable via
  `is_direct_summand` (P37) and this gate flips.

**No code, no test, no commit for Task F.** The deferral is delivered by the
verification-page note in the final task (the standing "every plan records its
honest scope" rule is satisfied by the named successor, not a silent skip).

---

### Task G: finitistic-dimension bounds

**Files:**
- Modify: `src/quiverlab/modules/homdims.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.finitistic_dimension_bounds`)
- Test: `tests/modules/test_finitistic.py`

**Interfaces:**
- Consumes: `ext.py::global_dimension(A, bound)` (the exact-or-lower-bound
  primitive), `syzygy` (Task A), `igusa_todorov_psi` (Task B),
  `Module.projective_resolution(bound).pd()`, `radtopsoc` radical layers.
- Produces:
  ```python
  @dataclass
  class FinitisticBounds:
      lower: int          # max finite pd found over the probe set (a genuine
                          # lower bound: findim(A) = sup{pd N : pd N < inf})
      upper: int | None   # gl.dim.value when gl.dim exact-finite (rigorous:
                          # findim <= gl.dim, equality when finite); else the
                          # Igusa-Todorov psi-bound OR None (honest degrade)
      exact: bool         # True iff lower == upper both finite (gl.dim resolved)
      note: str           # provenance of `upper`
      # __repr__: exact -> "findim = {lower}"; else "in [{lower}, {upper}]" /
      #           ">= {lower} (upper bound undetermined)" when upper is None
  def finitistic_dimension_bounds(A, bound=32) -> FinitisticBounds
  ```

  **Lower bound (rigorous, always):** the maximum finite projective dimension over
  the probe set `{ S_v } ∪ { Ω S_v } ∪ { radical layers of each P_v }`. By
  definition `findim(A) = sup{ pd N : pd N < ∞ }`, so any finite pd found is a
  valid lower bound.

  **Upper bound:**
  - If `gl.dim A` is **exact-finite**: `upper = gl.dim.value`, `exact = True`
    (findim ≤ gl.dim with equality when gl.dim < ∞ — a theorem, no folklore).
  - Else (gl.dim only a lower bound): the intended Igusa–Todorov bound is
    `upper = ψ(⊕_v Ω S_v) + 1` (Igusa–Todorov 2005, the ψ-finitistic bound —
    Design-decision-8's spec). **This is the one theorem this plan does not
    author from first principles; the implementer MUST, in the implementation
    session:** (1) pin the exact statement (module = the syzygies of the simples;
    the additive constant) verbatim from Igusa–Todorov 2005 in the docstring,
    (2) ARBITRATE the constant against the Barrios–Mata truncated-path closed
    forms where `findim` is known (the oracle decides `+1` vs `+0`, exactly as
    P39 arbitrates its Hom-complex sign), and (3) HONEST-DEGRADE: if the exact
    statement cannot be certified as a genuine theorem within the session, set
    `upper = None`, `exact = False`, and `note` = "upper bound undetermined
    (Igusa–Todorov ψ-bound not certified for this presentation)" — a `None`
    upper, never a folklore number. A safety gate the implementation MUST enforce:
    a numeric `upper` is always `≥ lower` (a violated inequality means the bound
    is mis-implemented — raise, do not clamp).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_finitistic.py
"""Finitistic-dimension bounds. Exact when gl.dim is finite (findim = gl.dim);
honest [lower, upper] otherwise. Barrios-Mata truncated closed forms."""
import pytest

from quiverlab import Quiver, TruncatedPathAlgebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import finitistic_dimension_bounds

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_hereditary_exact():
    A = _kA2()                                          # gl.dim 1, findim 1
    fb = finitistic_dimension_bounds(A)
    assert fb.lower == 1 and fb.upper == 1 and fb.exact is True


def test_nakayama_exact_equals_gldim():
    A = TruncatedPathAlgebra("A3", 2, field=QQ)         # gl.dim 2, findim 2
    fb = finitistic_dimension_bounds(A)
    assert fb.lower == 2 and fb.upper == 2 and fb.exact is True


def test_self_injective_lower_zero_upper_valid():
    # k[x]/(x^3): findim = 0 (only projectives have finite pd), gl.dim infinite.
    A = truncated_polynomial(3, field=QQ)
    fb = finitistic_dimension_bounds(A)
    assert fb.lower == 0
    assert fb.exact is False
    # upper is a valid bound (>= lower) or an honest None
    assert fb.upper is None or fb.upper >= fb.lower


def test_wrapper_matches():
    A = _kA2()
    fb = A.finitistic_dimension_bounds()
    assert fb.exact and fb.lower == 1
```

- [ ] **Step 2: Run to verify failure** — `ImportError: finitistic_dimension_bounds`

- [ ] **Step 3: Implement** (append to `homdims.py`; wrapper in `core/algebra.py`):

```python
def _finite_pd_probe_lower(A, bound):
    """Max finite pd over { S_v } U { Omega S_v } U { radical layers of P_v }."""
    lower = 0
    probes = []
    for v in A.quiver.vertices:
        Sv = A.simple(v)
        probes.append(Sv)
        probes.append(syzygy(Sv))
        layer = A.projective(v)                # walk rad^k P_v until zero
        while layer.dim:
            probes.append(layer)
            layer = layer.radical()
    for N in probes:
        if N.dim == 0:
            continue
        pd = N.projective_resolution(bound).pd()
        if pd is not None:
            lower = max(lower, pd)
    return lower


def finitistic_dimension_bounds(A, bound=32):
    from quiverlab.modules.ext import global_dimension
    lower = _finite_pd_probe_lower(A, bound)
    g = global_dimension(A, bound=bound)
    if g.exact:
        upper = g.value
        # findim <= gl.dim, equality when gl.dim finite; lower already >= any
        # finite pd, so lower == upper == gl.dim here (assert the sanity gate).
        assert upper >= lower, "finitistic lower bound exceeded gl.dim (bug)"
        return FinitisticBounds(lower, upper, exact=True,
                                note="gl.dim exact and finite => findim = gl.dim")
    # gl.dim only a lower bound: the Igusa-Todorov psi-bound (see the Interfaces
    # note -- pin the theorem verbatim, arbitrate the constant, degrade honestly).
    upper = _igusa_todorov_finitistic_upper(A, bound)   # int or None
    if upper is not None and upper < lower:
        raise QuiverlabError(                           # the mandated sanity gate
            "finitistic upper bound < lower bound: the Igusa-Todorov bound is "
            "mis-implemented for this presentation")
    return FinitisticBounds(
        lower, upper, exact=False,
        note=("Igusa-Todorov psi-bound" if upper is not None
              else "upper bound undetermined (psi-bound not certified)"))
```

`_igusa_todorov_finitistic_upper(A, bound)` computes `ψ(⊕_v Ω S_v) + const` per
the pinned theorem, or returns `None` on the honest degrade. Define
`FinitisticBounds` beside the other dataclasses. Wrapper in `core/algebra.py`.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/homdims.py src/quiverlab/core/algebra.py tests/modules/test_finitistic.py
git commit -m "feat(modules): finitistic_dimension_bounds -- rigorous lower + gl.dim/Igusa-Todorov upper, honesty-degraded"
```

---

### Task H: QPA battery + the `homological_profile` no-code compute kind

**Files:**
- Create: `tests/qpa/test_homdims_qpa.py`
- Modify: `src/quiverlab/qpa/scripts.py`, `src/quiverlab/qpa/crosscheck.py`
  (`crosscheck_dominant_dimension`, `crosscheck_gorenstein` — probe-first)
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch`: kind `homological_profile`, `_snip` ~line 1713)
- Modify: `docs/gui/runner.py` (twin handler + ETA `"scalars"` entry ~line 865)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox, `S.ids`, push list, block renderer)
- Modify: `webapp/static/app.js` (webapp block renderer)
- Modify: `webapp/server/i18n/en.json` + `es.json`
- Modify: `src/quiverlab/trace/results_html.py` (report renderer branch)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py` (ONE new fixture, documented)
- Test: `tests/qpa/test_homdims_qpa.py`, `tests/webapp/test_homological_profile_p40.py`

**Interfaces:**
- Consumes (QPA): probe live for the exact verb names via the `NamesGVars()`
  precedent (`tests/qpa/test_products_qpa.py`) — expected
  `GlobalDimensionOfAlgebra`, `DominantDimensionOfAlgebra`,
  `IsGorensteinAlgebra` (QPA manual Ch. 3/6); use what the probe finds. The
  `crosscheck_inj_dimension` template (`qpa/crosscheck.py:186`) is the shape:
  build the algebra script, run the verb, `int(val)` or `None` on GAP `false`.
  QPA has NO Igusa–Todorov surface — the probe SKIPS honestly for φ/ψ (the
  Task-B literature battery is their oracle), FAILING only if that ever changes.
- Consumes (GUI): the `_dispatch` scalar-handler pattern
  (`spec.py:1185` `global_dimension` — payload + `references` +
  `citations`); the runner twin (`docs/gui/runner.py:641`); the seven-touchpoint
  checklist (P38 Task 6). `homological_profile` is an **algebra-level scalar
  kind** — it routes through `_dispatch`, NOT `_dispatch_module`, and is NOT
  added to `MODULE_KINDS`. **Schema stays v1** (no request block).
- Produces: one `homological_profile` block aggregating the family, each entry
  honest-noted:
  ```python
  {"kind": "homological_profile",
   "global_dimension": {"text": str, "exact": bool, "value": int},
   "finitistic": {"lower": int, "upper": int | None, "exact": bool, "note": str},
   "dominant": {"value": int | None, "exact": bool, "infinite": bool, "text": str},
   "gorenstein": {"right_id": int | None, "left_id": int | None,
                  "is_gorenstein": bool | None, "text": str},
   "igusa_todorov": {"module": "(+)_v S_v", "phi": int, "psi": int}
                    | {"error": "<the loud decompose/iso message>"},
   "references": ["igusa_todorov", "assem_book", ...], "citations": [...]}
  ```
  The `igusa_todorov` entry of `⊕_v S_v` may RAISE over `char ≤ dim` — caught and
  reported as `{"error": ...}` (the Plan-30 τ-block honest-per-entry precedent),
  never a silent omission. Both runners byte-identical on the block.

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated
  dir; copy `tests/webapp/test_module_blocks_m0729.py`'s runner-pair fixture):

```python
# tests/webapp/test_homological_profile_p40.py
"""homological_profile scalar kind: served by hpc.spec, mirrored by the Pyodide
twin, honest per-entry notes."""
import json


def _req_kA2():
    # the standard request-dict fixture for kA2 over GF(7) with
    # compute = ["homological_profile"] (copy the shape from
    # test_module_blocks_m0729's request builder)
    ...


def test_profile_block_shape(tmp_path):
    from quiverlab.hpc.spec import run_spec           # or the module's run entry
    out = run_spec(_req_kA2(), tmp_path)
    block = _find_block(out, "homological_profile")
    assert block["global_dimension"]["value"] == 1
    assert block["finitistic"]["exact"] is True
    assert block["dominant"]["value"] == 1
    assert block["gorenstein"]["is_gorenstein"] is True
    assert "igusa_todorov" in block


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py's dispatch the way
    # test_module_blocks_m0729 does; assert json.dumps(sort_keys=True) equality
    # on the homological_profile block (both runners byte-identical).
    ...
```

(Write the `...` bodies concretely by copying the m0729 runner-pair fixture —
that file is THE pattern for cross-runner block contracts. kA2 over GF(7): φ/ψ of
⊕S_v are over char 7 > dim of each simple, so the IT entry is a value, not an
error.)

- [ ] **Step 2: Probe QPA + write the qpa battery** (concrete tests in the
  `tests/qpa/test_tor_qpa.py` style; header
  `pytestmark = pytest.mark.skipif(session.should_skip_qpa(), ...)`):
  compare `A.global_dimension()`, `A.dominant_dimension()`, `A.is_gorenstein()`
  against the live verbs the probe found, over the zoo (kA2/kA3-rel/
  `line_abc_cde`); the IT probe skips honestly (no QPA surface).

- [ ] **Step 3: Implement the handler + wiring.** In `spec.py::_dispatch`, after
  the `global_dimension` branch (`spec.py:1185`), add the `homological_profile`
  handler (call the five library functions, wrap the honest markers into the
  block, catch the IT loud refusal into `{"error": ...}`); mirror byte-for-byte
  in `docs/gui/runner.py` after its `global_dimension` elif
  (`docs/gui/runner.py:641`). Then the GUI wiring: checkbox id
  `qlgui-homological_profile` + `S.ids` + push-list in both `gui.js` files,
  block renderer in both `gui.js` + `app.js`, ETA scalar
  (`"homological_profile": 2.0` in `runner.py`'s `ETA_MODEL["scalars"]`), i18n
  keys (`inv.homological_profile`, `block.homological_profile.title` and per-row
  labels — EN + ES), `_snip` recipe
  (`"homological_profile": lambda it: "A.global_dimension(), A.dominant_dimension(), A.gorenstein_dimension(), A.finitistic_dimension_bounds()"`),
  and the `results_html.py` report branch.

- [ ] **Step 4: Add ONE golden fixture** `homological_profile_kA2` to
  `_runner_goldens.json`; note it in `test_runner_delegation.py`'s docstring
  change-log. Run the delegation test BEFORE adding it to confirm existing
  entries stay byte-identical.

- [ ] **Step 5: Run the gates**

Run: `... -m pytest tests/webapp/test_homological_profile_p40.py tests/webapp/test_runner_delegation.py tests/hpc -q`
then `... -m pytest tests/qpa/test_homdims_qpa.py -q -m qpa` (venv has [qpa]).
Expected: PASS; the IT-vs-QPA probe SKIPS honestly.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc,qpa): homological_profile compute kind + domdim/Gorenstein QPA battery"
```

---

### Task I: verification page, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Modify: `src/quiverlab/citations/registry.py` (+ `references.bib`) — add the
  keys this plan cites
- Test: existing release gates

- [ ] **Step 1: Citations.** Add (verified BibTeX, `_r(...)` registry precedent
  at `registry.py:120-164`, only keys that resolve — `_citation_pairs` swallows a
  `KeyError` but the plan must not rely on that): `igusa_todorov`
  (Igusa–Todorov, Fields Inst. Commun. 45, 2005), `barrios_mata`
  (Barrios–Mata, on the φ function for truncated path algebras, 2019),
  `gelinas_delooping` (Gélinas, JPAA 2022 — cited only in the Task-F honest-scope
  note), and a dominant-dimension / Iwanaga–Gorenstein reference (reuse
  `assem_book` where it covers the definition; add a dedicated key only if a
  primary source is BibTeX-verifiable).
- [ ] **Step 2: Verification page.** Add the Plan-40 subsystem row
  (`modules/homdims.py` | oracles: `oracle_selfcert` φ=ψ=pd identity +
  periodicity `is_isomorphic` certificates + additivity; `oracle_literature`
  Barrios–Mata truncated φ/ψ, hereditary/self-injective domdim & Gorenstein,
  cyclic-Nakayama period-from-Kupisch; `qpa` domdim/Gorenstein/gl.dim
  crosschecks). Add the **honest-scope entries**: (a) Task F delooping level is
  DEFERRED with the recorded reason (existential-`N` quantifier has no crisp
  bounded procedure here; successor plan builds the injective-side tower);
  (b) `is_gorenstein` is three-valued True/None — a `False` verdict needs a proof
  of infinite injective dimension the bounded engine does not furnish; (c) QPA
  has NO Igusa–Todorov surface (the literature battery is its oracle). Recount
  the class table (`tests/release/test_oracle_classes.py` drives the numbers —
  run collection, paste, re-run to green).
- [ ] **Step 3: README.** One features line: "finitistic / dominant / Gorenstein
  dimensions, Igusa–Todorov φ/ψ, Ω/τ-periodicity certificates — the C6
  homological-dimensions family, clickable via `homological_profile`."
- [ ] **Step 4: Full gate:**
  `... -m pytest tests/modules -q` (deep, the touched files),
  `... -m pytest -q -m fast`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release -q` — all green.
- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-40 oracle rows + honest-scope (delooping deferred, is_gorenstein three-valued) + recounted classes"
```

---

## Acceptance (Plan-40 definition of done)

1. `syzygy`/`cosyzygy` public (byte-stable refactor — engine/resolution suites
   green); `igusa_todorov_phi`/`_psi`, `dominant_dimension`,
   `gorenstein_dimension`/`is_gorenstein`, `omega_periodicity`/`tau_periodicity`,
   `finitistic_dimension_bounds` all public on `Algebra`/`Module`, all loudly-
   validated, all following the `GlobalDimension` certified-value-or-honest-marker
   pattern (never a bare number when unresolved; never a fabricated ∞).
2. φ=ψ=pd (finite pd), Barrios–Mata truncated φ/ψ closed forms, projective
   additivity, hereditary/self-injective domdim & Gorenstein values, and
   cyclic-Nakayama period-from-Kupisch all pinned (`oracle_literature` /
   `oracle_selfcert`); the decompose char-caveat inherited and tested (loud
   refusal over `char ≤ dim`).
3. QPA domdim/Gorenstein/gl.dim battery green live (`-m qpa`); the IT-vs-QPA
   probe skips honestly and FAILS if QPA ever ships an IT surface.
4. `homological_profile` clickable end-to-end (GUI canvas → block → report) in
   EN+ES, both runners byte-identical, ONE golden added with a documented
   change-log entry, schema still v1.
5. The Igusa–Todorov finitistic UPPER bound is either the pinned theorem
   (arbitrated against Barrios–Mata) or an honest `None` — never folklore; the
   `upper ≥ lower` sanity gate is live.
6. Task F (delooping level) DEFERRED with the recorded reason on the verification
   page (honest scope satisfied by the named successor, not a silent skip);
   `is_gorenstein`'s three-valued honesty documented there too.
7. `docs/verification.md` recounted; README line added; deep (touched dirs) +
   fast + qpa + release suites green.
</content>
</invoke>
