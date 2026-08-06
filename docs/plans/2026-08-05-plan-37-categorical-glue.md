# Plan 37: C1 Categorical Glue — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Module homomorphisms become first-class objects — `ModuleHom` with
kernel/image/cokernel as Modules-with-maps, SES + split tests +
pushout/pullback, `End(M)` as an `Algebra`, direct sums with
inclusion/projection maps, covers/envelopes as maps, radical/socle series —
the composition layer every later plan (P39, P41, P44, P45, P47) builds on.

**Architecture:** A new `src/quiverlab/modules/morphism.py` (the `ModuleHom`
class + hom-basis constructors) and `src/quiverlab/modules/ses.py` (SES,
split test, pushout/pullback), both thin exact-linear-algebra layers over
the EXISTING primitives: `hom.py::hom_space` (dense `tgt.dim × src.dim`
matrices), `radtopsoc.py::submodule/quotient`, `linalg_mod.py`
(`kernel_columns`, `column_space_pivots`, `solve_columns`),
`yoneda.py::_quotient_with_maps` (cokernel with proj+lift — the template),
and `Algebra.from_structure_constants` for `End(M)`. No new math engines.

**Tech Stack:** pure exact linear algebra over `Domain` (`linalg_mod`),
sympy only where already used. No floats in `src/`.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- `tests/modules/` auto-assigns to the **deep** bucket — run new tests by
  path during development, and finish with a `-m deep` spot-run of the
  touched files.
- House conventions: a module-map matrix is dense `tgt.dim × src.dim` over
  `M.domain`, columns = source ambient basis (exactly `hom_space`'s layout;
  `tgt.action[b] @ f == f @ src.action[b]` for every generator).
  Hom/comparisons refuse loudly across sides/algebras via
  `hom.py::_assert_comparable`.
- Composition is written **left-to-right like paths**: `f.then(g)` for
  `M --f--> N --g--> P` (matrix `g.matrix @ f.matrix`). Never overload `*`.
- Plan-32 markers: exactness/certificate tests = `oracle_selfcert`; QPA
  comparisons live in `tests/qpa/` (bucket = the class, never double-marked);
  theory pins = `oracle_literature`.
- Every plan merge updates `docs/verification.md` (new oracle rows + counts).
- Conventional commits; green at every commit; branch `plan-37-categorical-glue` off `dev`.

---

### Task 1: `ModuleHom` — the first-class morphism

**Files:**
- Create: `src/quiverlab/modules/morphism.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.hom_basis(M, N)` beside `hom` at ~line 405)
- Modify: `src/quiverlab/modules/module.py` (add `Module.identity_hom()`)
- Test: `tests/modules/test_morphism.py`

**Interfaces:**
- Consumes: `hom.py::hom_space(M, N)` (list of dense matrices),
  `hom.py::_assert_comparable(M, N, op)`,
  `yoneda.py::_is_module_map(f, src, tgt, dom)` (the validation predicate —
  import it, do not reimplement), `linalg_mod` (`mat_rank`, `matmul`).
- Produces:
  ```python
  class ModuleHom:
      def __init__(self, src, tgt, matrix, check=True)  # QuiverlabError if not a module map
      src, tgt, matrix, domain                          # attributes
      def then(self, g) -> "ModuleHom"                  # self then g (left-to-right)
      def rank(self) -> int
      def is_zero(self) -> bool
      def is_mono(self) -> bool                         # rank == src.dim
      def is_epi(self) -> bool                          # rank == tgt.dim
      def is_iso(self) -> bool
      def __repr__(self)                                # "Hom(M -> N, rank r)"
  def hom_basis(M, N) -> list[ModuleHom]                # wraps hom_space
  def zero_hom(M, N) -> ModuleHom
  ```
  `Algebra.hom_basis(M, N)` delegates to `morphism.hom_basis`.
  `Module.identity_hom()` returns the identity `ModuleHom`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_morphism.py
"""First-class module homomorphisms (Plan 37 / C1). Self-certifying:
the constructor validates the intertwining relations exactly."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.morphism import ModuleHom, hom_basis, zero_hom

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def test_hom_basis_wraps_hom_space_with_validation():
    A = _a2()
    P1, S2 = A.projective(1), A.simple(2)
    basis = hom_basis(P1, S2)
    assert len(basis) == A.hom(P1, S2)          # dims agree with the old API
    for f in basis:
        assert f.src is P1 and f.tgt is S2


def test_constructor_rejects_non_module_map():
    A = _a2()
    P1, S1 = A.projective(1), A.simple(1)       # dims 2 and 1
    # [[1, 1]]: P1 -> S1 hitting the vertex-2 coordinate -- breaks the e_1
    # idempotent intertwining relation, so the constructor must refuse.
    with pytest.raises(QuiverlabError, match="module map"):
        ModuleHom(P1, S1, [[1, 1]], check=True)


def test_then_composes_left_to_right():
    A = _a2()
    P1 = A.projective(1)
    top = P1.top()                               # = S1
    # P1 ->> top(P1): from the projective-cover direction we at least have
    # SOME epi in the hom basis; compose with an endo of top.
    fs = [f for f in hom_basis(P1, top) if f.is_epi()]
    assert fs, "an epi P1 ->> top(P1) must exist"
    f = fs[0]
    idt = top.identity_hom()
    g = f.then(idt)
    assert g.matrix == f.matrix and g.src is P1 and g.tgt is top


def test_then_refuses_mismatched_middle():
    A = _a2()
    P1, S2 = A.projective(1), A.simple(2)
    f = zero_hom(P1, S2)
    with pytest.raises(QuiverlabError, match="compose"):
        f.then(zero_hom(P1, S2))                 # tgt S2 != src P1


def test_mono_epi_iso_flags():
    A = _a2()
    S2 = A.simple(2)
    idm = S2.identity_hom()
    assert idm.is_mono() and idm.is_epi() and idm.is_iso()
    z = zero_hom(S2, S2)
    assert z.is_zero() and not z.is_mono()


def test_cross_algebra_refused():
    A, B = _a2(), _a2()
    with pytest.raises(QuiverlabError):
        hom_basis(A.simple(1), B.simple(1))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_morphism.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.modules.morphism`

- [ ] **Step 3: Implement `src/quiverlab/modules/morphism.py`**

```python
"""First-class module homomorphisms (Plan 37 / C1).

A ModuleHom is (src, tgt, matrix): matrix is dense tgt.dim x src.dim over
the shared Domain, satisfying tgt.action[b] @ matrix == matrix @ src.action[b]
for every generator label -- validated at construction (check=True).
Composition is written left-to-right like paths: f.then(g)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import _assert_comparable, hom_space
from quiverlab.modules.yoneda import _is_module_map


class ModuleHom:
    def __init__(self, src, tgt, matrix, check=True):
        _assert_comparable(src, tgt, "Hom")
        self.src, self.tgt = src, tgt
        self.domain = src.domain
        self.matrix = [[self.domain.coerce(x) for x in row] for row in matrix]
        if len(self.matrix) != tgt.dim or (tgt.dim and any(
                len(r) != src.dim for r in self.matrix)):
            raise QuiverlabError(
                f"not a module map: expected {tgt.dim}x{src.dim} matrix")
        if check and not _is_module_map(self.matrix, src, tgt, self.domain):
            raise QuiverlabError(
                "not a module map: the matrix does not intertwine the actions")

    def then(self, g: "ModuleHom") -> "ModuleHom":
        """self then g (left-to-right): src --self--> tgt --g--> g.tgt."""
        if g.src is not self.tgt:
            raise QuiverlabError("cannot compose: middle modules differ "
                                 f"({self.tgt.name} vs {g.src.name})")
        return ModuleHom(self.src, g.tgt,
                         lm.matmul(g.matrix, self.matrix, self.domain),
                         check=False)

    def rank(self) -> int:
        return lm.mat_rank(self.matrix, self.domain)

    def is_zero(self) -> bool:
        z = self.domain.zero
        return all(x == z for row in self.matrix for x in row)

    def is_mono(self) -> bool:
        return self.rank() == self.src.dim

    def is_epi(self) -> bool:
        return self.rank() == self.tgt.dim

    def is_iso(self) -> bool:
        return self.src.dim == self.tgt.dim and self.rank() == self.src.dim

    def __repr__(self):
        return (f"Hom({self.src.name} -> {self.tgt.name}, "
                f"rank {self.rank()})")


def hom_basis(M, N):
    """Basis of Hom_A(M, N) as validated ModuleHom objects."""
    _assert_comparable(M, N, "Hom")
    return [ModuleHom(M, N, mat, check=False) for mat in hom_space(M, N)]


def zero_hom(M, N):
    _assert_comparable(M, N, "Hom")
    z = N.domain.zero
    return ModuleHom(M, N, [[z] * M.dim for _ in range(N.dim)], check=False)
```

Wire the two delegations: `Algebra.hom_basis(self, M, N)` (import inside the
method, mirroring `Algebra.hom` at `core/algebra.py:405`), and
`Module.identity_hom(self)` returning
`ModuleHom(self, self, identity-matrix, check=False)`.
**Adjust to reality:** if `Domain` has no `.coerce`/`.zero` with these exact
names, use whatever `linalg_mod`/`module.py::_coerce_matrix` actually use —
read those first; `_coerce_matrix` is the precedent for entry coercion.

- [ ] **Step 4: Run tests to verify they pass**

Run: `... -m pytest tests/modules/test_morphism.py -v` — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/morphism.py src/quiverlab/core/algebra.py \
        src/quiverlab/modules/module.py tests/modules/test_morphism.py
git commit -m "feat(modules): ModuleHom -- validated first-class morphisms with left-to-right composition"
```

---

### Task 2: kernel / image / cokernel as Modules-with-maps

**Files:**
- Modify: `src/quiverlab/modules/morphism.py`
- Test: `tests/modules/test_kernel_image_cokernel.py`

**Interfaces:**
- Consumes: `linalg_mod::kernel_columns`, `column_space_pivots`,
  `solve_columns`, `cols_to_matrix`; `radtopsoc::submodule, quotient`;
  `yoneda::_quotient_with_maps(ambient, sub_cols, dom, name)` → `(Q, proj, lift)`.
- Produces (methods on `ModuleHom`):
  ```python
  def kernel(self) -> tuple["Module", "ModuleHom"]      # (K, iota: K -> src), iota mono
  def image(self) -> tuple["Module", "ModuleHom", "ModuleHom"]
      # (I, epi: src -> I, iota: I -> tgt) with epi.then(iota) == self
  def cokernel(self) -> tuple["Module", "ModuleHom"]    # (C, proj: tgt -> C), proj epi
  ```
  Names: `K = f"ker({src.name}->{tgt.name})"`, `im(...)`, `coker(...)`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_kernel_image_cokernel.py
"""Kernel/image/cokernel of a ModuleHom, self-certified by rank-nullity,
mono/epi flags, and the epi-mono factorization identity."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.morphism import hom_basis

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def _all_homs(A):
    mods = [A.projective(v) for v in (1, 2, 3)] + [A.simple(v) for v in (1, 2, 3)]
    for M in mods:
        for N in mods:
            yield from hom_basis(M, N)


def test_rank_nullity_and_factorization_battery():
    A = _a3()
    seen = 0
    for f in _all_homs(A):
        K, iota = f.kernel()
        I, epi, mono = f.image()
        C, proj = f.cokernel()
        assert K.dim + f.rank() == f.src.dim              # rank-nullity
        assert I.dim == f.rank()
        assert C.dim == f.tgt.dim - f.rank()
        assert iota.is_mono() and epi.is_epi() and mono.is_mono() and proj.is_epi()
        assert epi.then(mono).matrix == f.matrix          # f = mono . epi
        assert iota.then(f).is_zero()                     # f . iota = 0
        assert f.then(proj).is_zero()                     # proj . f = 0
        seen += 1
    assert seen >= 20                                      # battery is nonempty


def test_kernel_of_projective_cover_is_radical_syzygy():
    A = _a3()
    S1 = A.simple(1)
    from quiverlab.modules.resolution import projective_cover
    Q0, d0, _ = projective_cover(S1)
    from quiverlab.modules.morphism import ModuleHom
    f = ModuleHom(Q0, S1, d0)
    K, _ = f.kernel()
    # ker(P(S1) ->> S1) = rad P(S1)
    assert K.dimension_vector() == Q0.radical().dimension_vector()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `... -m pytest tests/modules/test_kernel_image_cokernel.py -v`
Expected: FAIL — `AttributeError: 'ModuleHom' object has no attribute 'kernel'`

- [ ] **Step 3: Implement** (append to `morphism.py`)

```python
    def kernel(self):
        """(K, iota) with K = ker(self) a submodule of src, iota mono."""
        from quiverlab.modules import radtopsoc
        cols = lm.kernel_columns(self.matrix, self.domain)
        K = radtopsoc.submodule(self.src, cols,
                                name=f"ker({self.src.name}->{self.tgt.name})")
        iota = ModuleHom(K, self.src, lm.cols_to_matrix(cols, self.src.dim,
                                                        self.domain),
                         check=False)
        return K, iota

    def image(self):
        """(I, epi, mono): I = im(self) as a submodule of tgt,
        self == epi.then(mono)."""
        from quiverlab.modules import radtopsoc
        # image columns: the columns of the matrix, reduced to a basis
        all_cols = [[row[j] for row in self.matrix] for j in range(self.src.dim)]
        pivots = lm.column_space_pivots(self.matrix, self.domain)
        cols = [all_cols[j] for j in pivots]
        I = radtopsoc.submodule(self.tgt, cols,
                                name=f"im({self.src.name}->{self.tgt.name})")
        mono = ModuleHom(I, self.tgt, lm.cols_to_matrix(cols, self.tgt.dim,
                                                        self.domain),
                         check=False)
        # epi: express each f(basis_j of src) in the chosen image basis
        coeffs = lm.solve_columns(mono.matrix, all_cols, self.domain)
        epi = ModuleHom(self.src, I,
                        [[coeffs[j][i] for j in range(self.src.dim)]
                         for i in range(I.dim)], check=False)
        return I, epi, mono

    def cokernel(self):
        """(C, proj) with C = tgt/im(self), proj epi."""
        from quiverlab.modules.yoneda import _quotient_with_maps
        all_cols = [[row[j] for row in self.matrix] for j in range(self.src.dim)]
        C, proj_mat, _lift = _quotient_with_maps(
            self.tgt, all_cols, self.domain,
            name=f"coker({self.src.name}->{self.tgt.name})")
        return C, ModuleHom(self.tgt, C, proj_mat, check=False)
```

**Adjust to reality:** read `yoneda._quotient_with_maps`'s actual return
types (`proj` may already be a raw matrix, and it may not take zero-column
input) and `radtopsoc.submodule`'s handling of an empty spanning set (the
zero module) — mirror whatever `radtopsoc.socle` does for the empty case.

- [ ] **Step 4: Run tests, verify pass; run the neighbors**

Run: `... -m pytest tests/modules/test_morphism.py tests/modules/test_kernel_image_cokernel.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/morphism.py tests/modules/test_kernel_image_cokernel.py
git commit -m "feat(modules): kernel/image/cokernel as Modules with certified connecting maps"
```

---

### Task 3: short exact sequences, split test, pushout/pullback

**Files:**
- Create: `src/quiverlab/modules/ses.py`
- Test: `tests/modules/test_ses.py`

**Interfaces:**
- Consumes: `ModuleHom` (Tasks 1–2), `hom.py::hom_space`,
  `yoneda.py::YonedaSequence.check_exact` (the rank-identity template —
  reuse its logic shape, not the class), `linalg_mod::solve_columns`.
- Produces:
  ```python
  class ShortExactSequence:                 # 0 -> L --f--> M --g--> N -> 0
      def __init__(self, f: ModuleHom, g: ModuleHom, check=True)
          # QuiverlabError unless f mono, g epi, im f = ker g (rank identity
          # dim M = dim L + dim N plus g.then? NO -- f.then(g).is_zero()
          # and rank f + rank g == M.dim)
      L, M, N, f, g                          # attributes
      def is_split(self) -> bool             # a section s: N->M with s.then(g)=id exists
      def __repr__(self)
  def pushout(f: ModuleHom, g: ModuleHom) -> tuple  # f: A->B, g: A->C shared src
      # returns (P, inB: B->P, inC: C->P) with f.then(inB) == g.then(inC)
  def pullback(f: ModuleHom, g: ModuleHom) -> tuple # f: B->D, g: C->D shared tgt
      # returns (P, prB: P->B, prC: P->C) with prB.then(f) == prC.then(g)
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_ses.py
"""SES objects: exactness certified at construction; split test; pushout/
pullback squares certified by their universal-square identities."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.morphism import ModuleHom, hom_basis, zero_hom
from quiverlab.modules.ses import ShortExactSequence, pullback, pushout

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def _rad_ses(A):
    """0 -> rad P1 -> P1 -> S1 -> 0 for kA2 (rad P1 = S2, dim 1)."""
    P1, S1 = A.projective(1), A.simple(1)
    R = P1.radical()
    iota = next(f for f in hom_basis(R, P1) if f.is_mono())
    epi = next(f for f in hom_basis(P1, S1) if f.is_epi())
    return ShortExactSequence(iota, epi)


def test_ses_certifies_exactness():
    ses = _rad_ses(_a2())
    assert ses.M.dim == ses.L.dim + ses.N.dim


def test_ses_rejects_non_exact():
    A = _a2()
    S1, S2 = A.simple(1), A.simple(2)
    with pytest.raises(QuiverlabError, match="exact"):
        ShortExactSequence(zero_hom(S1, S2), zero_hom(S2, S1))


def test_split_iff_ext_vanishes():
    A = _a2()
    ses = _rad_ses(A)
    # Ext^1(S1, S2) = k for kA2 (the arrow) => the rad sequence is NOT split
    assert A.ext(ses.N, ses.L, 1) == 1
    assert ses.is_split() is False


def test_direct_sum_ses_splits():
    A = _a2()
    S1, S2 = A.simple(1), A.simple(2)
    from quiverlab.modules.morphism import direct_sum          # Task 4 name
    D, (i1, i2), (p1, p2) = direct_sum(S2, S1)
    ses = ShortExactSequence(i1, p2)
    assert ses.is_split() is True


def test_pushout_square_commutes():
    A = _a2()
    P1 = A.projective(1)
    R = P1.radical()
    f = next(h for h in hom_basis(R, P1) if h.is_mono())
    g = R.identity_hom()
    P, inB, inC = pushout(f, g)
    assert f.then(inB).matrix == g.then(inC).matrix
    assert P.dim == P1.dim + R.dim - R.dim                  # pushout along id: P ~ P1


def test_pullback_square_commutes():
    A = _a2()
    P1, S1 = A.projective(1), A.simple(1)
    f = next(h for h in hom_basis(P1, S1) if h.is_epi())
    g = S1.identity_hom()
    P, prB, prC = pullback(f, g)
    assert prB.then(f).matrix == prC.then(g).matrix
```

(NOTE: `test_direct_sum_ses_splits` consumes Task 4's `direct_sum` — mark it
`@pytest.mark.skip(reason="Task 4")` when running Task 3 standalone and
unskip in Task 4, OR reorder: implement Task 4 before this one test. The
plan keeps the test here so the SES/split surface is specified in one
place; executors: leave it skipped until Task 4 lands, then unskip.)

- [ ] **Step 2: Run to verify failure** — Expected: `ModuleNotFoundError: ses`

- [ ] **Step 3: Implement `src/quiverlab/modules/ses.py`**

```python
"""Short exact sequences, split test, pushout/pullback (Plan 37 / C1)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.morphism import ModuleHom, hom_basis


class ShortExactSequence:
    """0 -> L --f--> M --g--> N -> 0, exactness certified at construction:
    f mono, g epi, g.f = 0, rank f + rank g = dim M (=> im f = ker g)."""

    def __init__(self, f: ModuleHom, g: ModuleHom, check=True):
        if f.tgt is not g.src:
            raise QuiverlabError("not composable: f.tgt is not g.src")
        self.f, self.g = f, g
        self.L, self.M, self.N = f.src, f.tgt, g.tgt
        if check:
            ok = (f.is_mono() and g.is_epi() and f.then(g).is_zero()
                  and f.rank() + g.rank() == self.M.dim)
            if not ok:
                raise QuiverlabError("sequence is not exact "
                                     "(mono/epi/g.f=0/rank identity failed)")

    def is_split(self) -> bool:
        """True iff a section s: N -> M with s.then(g) = id_N exists.
        Linear in s: expand s over the Hom(N, M) basis and solve."""
        basis = hom_basis(self.N, self.M)
        if not basis:
            return self.N.dim == 0
        dom = self.M.domain
        n = self.N.dim
        # columns: vec(s_k.then(g).matrix) for each basis element; target: vec(id)
        cols = []
        for s in basis:
            comp = s.then(self.g).matrix
            cols.append([comp[i][j] for j in range(n) for i in range(n)])
        target = [dom.one if i == j else dom.zero
                  for j in range(n) for i in range(n)]
        sol = lm.solve_columns(lm.cols_to_matrix(cols, n * n, dom),
                               [target], dom)
        return sol is not None


def pushout(f: ModuleHom, g: ModuleHom):
    """Pushout of B <--f-- A --g--> C: P = (B (+) C) / <(f(a), -g(a))>."""
    if f.src is not g.src:
        raise QuiverlabError("pushout needs a shared source")
    from quiverlab.modules.morphism import direct_sum
    from quiverlab.modules.yoneda import _quotient_with_maps
    D, (iB, iC), _ = direct_sum(f.tgt, g.tgt)
    dom = D.domain
    A = f.src
    diag_cols = []
    for j in range(A.dim):
        fb = [f.matrix[i][j] for i in range(f.tgt.dim)]
        gc = [dom.negate(g.matrix[i][j]) for i in range(g.tgt.dim)]
        diag_cols.append(fb + gc)
    P, proj_mat, _ = _quotient_with_maps(D, diag_cols, dom, name="pushout")
    proj = ModuleHom(D, P, proj_mat, check=False)
    return P, iB.then(proj), iC.then(proj)


def pullback(f: ModuleHom, g: ModuleHom):
    """Pullback of B --f--> D <--g-- C: P = ker(B (+) C -> D)."""
    if f.tgt is not g.tgt:
        raise QuiverlabError("pullback needs a shared target")
    from quiverlab.modules.morphism import direct_sum
    D, (iB, iC), (pB, pC) = direct_sum(f.src, g.src)
    diff = ModuleHom(D, f.tgt,
                     [[*f.matrix[i], *[D.domain.negate(x) for x in g.matrix[i]]]
                      for i in range(f.tgt.dim)], check=False)
    P, iota = diff.kernel()
    return P, iota.then(pB), iota.then(pC)
```

(`direct_sum` and `Domain.negate` naming: Task 4 defines `direct_sum`;
check the Domain API for negation — if entries negate via
`dom.sub(dom.zero, x)` or sympy `-x`, use that. Executors implement Task 4
FIRST if running sequentially — the task order 3↔4 may be swapped freely;
both orders keep every commit green if the skip-note in Step 1 is honored.)

- [ ] **Step 4: Run tests** — Expected: PASS (minus the one Task-4 skip)

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ses.py tests/modules/test_ses.py
git commit -m "feat(modules): ShortExactSequence with certified exactness, split test, pushout/pullback"
```

---

### Task 4: direct sums with maps; `is_direct_summand`

**Files:**
- Modify: `src/quiverlab/modules/morphism.py`
- Test: `tests/modules/test_direct_sum.py`

**Interfaces:**
- Consumes: `yoneda.py::_direct_sum2(P, N)` (existing 2-ary direct sum —
  read it; reuse or generalize), `decompose.py::decompose`,
  `hom.py::is_isomorphic`.
- Produces:
  ```python
  def direct_sum(*modules) -> tuple
      # (D, inclusions, projections): D = M1 (+) ... (+) Mk,
      # inclusions[i]: Mi -> D, projections[i]: D -> Mi,
      # certified: proj_i . incl_i = id, sum incl_i . proj_i = id_D
  def is_direct_summand(N, M, budget=512) -> bool
      # Krull-Schmidt route: decompose both (Plan 30, certified), compare
      # multiplicities via is_isomorphic. Loud QuiverlabError when
      # decompose cannot certify within budget.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_direct_sum.py
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.morphism import direct_sum, is_direct_summand

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def test_direct_sum_identities():
    A = _a2()
    P1, S2, S1 = A.projective(1), A.simple(2), A.simple(1)
    D, incls, projs = direct_sum(P1, S2, S1)
    assert D.dim == P1.dim + S2.dim + S1.dim
    for i, (inc, prj) in enumerate(zip(incls, projs)):
        assert inc.then(prj).matrix == inc.src.identity_hom().matrix
    # sum of the idempotents is the identity on D
    dom = D.domain
    total = [[dom.zero] * D.dim for _ in range(D.dim)]
    for inc, prj in zip(incls, projs):
        e = prj.then(inc)          # D -> Mi -> D
        for r in range(D.dim):
            for c in range(D.dim):
                total[r][c] = dom.add(total[r][c], e.matrix[r][c])
    assert total == D.identity_hom().matrix


def test_is_direct_summand_true_and_false():
    A = _a2()
    P1, S1, S2 = A.projective(1), A.simple(1), A.simple(2)
    D, _, _ = direct_sum(P1, S1)
    assert is_direct_summand(S1, D) is True
    assert is_direct_summand(P1, D) is True
    # S2 = rad P1 is a SUBmodule of P1 but P1 is indecomposable:
    assert is_direct_summand(S2, P1) is False
```

- [ ] **Step 2: Run to verify failure** — `ImportError: direct_sum`

- [ ] **Step 3: Implement** (append to `morphism.py`; block-diagonal
  assembly — read `yoneda._direct_sum2` first and reuse its action-assembly
  idiom, generalized to k summands, then build the inclusion/projection
  matrices as the obvious block columns/rows; `is_direct_summand`):

```python
def is_direct_summand(N, M, budget=512):
    """Krull-Schmidt test: N | M iff every indecomposable summand of N
    appears in M with >= multiplicity. Certified by Plan-30 decompose;
    loud when decompose cannot certify."""
    from quiverlab.modules.decompose import decompose
    from quiverlab.modules.hom import is_isomorphic
    if N.dim == 0:
        return True
    dn = decompose(N, budget=budget)
    dm = list(decompose(M, budget=budget))
    for (ni, nmult) in dn:
        hit = None
        for k, (mi, mmult) in enumerate(dm):
            if ni.dim == mi.dim and is_isomorphic(ni, mi):
                hit = k
                break
        if hit is None or dm[hit][1] < nmult:
            return False
        dm[hit] = (dm[hit][0], dm[hit][1] - nmult)
    return True
```

Also add `Domain`-correct `add`/`negate` shims if the Domain API differs —
single chokepoint helpers `_dadd(dom, x, y)` / `_dneg(dom, x)` at module
top, used by Task 3 too.

- [ ] **Step 4: Run tests** (and UNSKIP `test_direct_sum_ses_splits` in
  `test_ses.py` if Task 3 landed first)

Run: `... -m pytest tests/modules/test_direct_sum.py tests/modules/test_ses.py -v`
Expected: PASS, no skips remaining

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/morphism.py tests/modules/
git commit -m "feat(modules): k-ary direct sums with certified inclusion/projection maps; Krull-Schmidt is_direct_summand"
```

---

### Task 5: `End(M)` as an `Algebra`

**Files:**
- Create: `src/quiverlab/modules/endomorphism.py`
- Modify: `src/quiverlab/modules/module.py` (add `Module.end_algebra()`)
- Test: `tests/modules/test_end_algebra.py`

**Interfaces:**
- Consumes: `hom_basis(M, M)` (Task 1), `Algebra.from_structure_constants(T, unit, field=..., check=True, basis_labels=...)`
  (`core/algebra.py:45` — presentation-less is legal; `check=True` validates
  associativity + unit, which is the self-certificate),
  `linalg_mod::solve_columns` (to express a product in the basis).
- Produces: `end_algebra(M) -> Algebra` — structure constants of
  `End_A(M)` in the `hom_basis(M, M)` basis, **composition order matching
  the house left-to-right rule**: the product `f*g` of basis elements is
  `f.then(g)` (document this in the docstring — it makes `End(A_A) ≅ A`,
  not `A^op`, under right-module conventions; VERIFY this in the test and
  flip to `g.then(f)` if the oracle says otherwise — the test is the
  arbiter, the docstring records the outcome).
  `Module.end_algebra()` delegates.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_end_algebra.py
"""End(M) as a structure-constant Algebra. Self-certified by
from_structure_constants(check=True) (associativity + unit); the
regular-module oracle pins the composition-order convention."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.endomorphism import end_algebra

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def test_end_of_simple_is_the_field():
    A = _a2()
    E = end_algebra(A.simple(1))
    assert E.dim == 1


def test_end_dim_matches_hom_dim():
    A = _a2()
    P1 = A.projective(1)
    E = end_algebra(P1)
    assert E.dim == A.hom(P1, P1)
    E._validate()                      # associativity + unit, loud on failure


def test_end_of_regular_module_has_algebra_dimension():
    # End_A(A_A) ~ A as a k-algebra (dim check + Loewy length agree)
    A = _a2()
    from quiverlab.modules.morphism import direct_sum
    regular, _, _ = direct_sum(A.projective(1), A.projective(2))
    E = end_algebra(regular)
    assert E.dim == 3                  # dim kA2 = 3
    assert E.loewy_length() == A.loewy_length() == 2
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: endomorphism`

- [ ] **Step 3: Implement `src/quiverlab/modules/endomorphism.py`**

```python
"""End_A(M) as a structure-constant Algebra (Plan 37 / C1).

Basis = hom_basis(M, M); the product of basis elements b_i * b_j is
b_i.then(b_j) (left-to-right, the path convention). The unit is id_M.
The returned Algebra is presentation-less: arithmetic, center,
loewy_length, decompose-style analysis work; .simple/.projective need a
quiver presentation and refuse loudly (builders._require_provenance)."""
from __future__ import annotations

from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.morphism import hom_basis


def end_algebra(M):
    if M.dim == 0:
        raise QuiverlabError("End of the zero module is the zero ring -- "
                             "not an Algebra with 1 here; refused.")
    basis = hom_basis(M, M)
    dom = M.domain
    n = len(basis)
    d = M.dim
    vec = lambda f: [f.matrix[i][j] for j in range(d) for i in range(d)]
    B = lm.cols_to_matrix([vec(f) for f in basis], d * d, dom)
    T = []
    for bi in basis:
        row = []
        for bj in basis:
            prod = bi.then(bj)
            coeffs = lm.solve_columns(B, [vec(prod)], dom)
            if coeffs is None:
                raise QuiverlabError("End(M) product left the Hom basis -- "
                                     "hom_space is inconsistent (bug)")
            row.append(list(coeffs[0]))
        T.append(row)
    unit = lm.solve_columns(B, [vec(M.identity_hom())], dom)[0]
    return Algebra.from_structure_constants(
        T, list(unit), field=dom, check=True,
        basis_labels=None)
```

**Adjust to reality:** `from_structure_constants(field=...)` — check
whether it takes a `Domain` instance or a field spec (read
`core/algebra.py:45-60`); pass whatever an existing caller passes.

- [ ] **Step 4: Run tests** — Expected: PASS. If
  `test_end_of_regular_module...` fails on `loewy_length`, the composition
  order is flipped: change `bi.then(bj)` to `bj.then(bi)`, update the
  module docstring's convention note, re-run — the oracle decides.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/endomorphism.py src/quiverlab/modules/module.py tests/modules/test_end_algebra.py
git commit -m "feat(modules): End(M) as a validated structure-constant Algebra"
```

---

### Task 6: covers/envelopes as maps; radical & socle series; composition factors

**Files:**
- Modify: `src/quiverlab/modules/module.py`
- Modify: `src/quiverlab/trace/modules.py` (delegate `_radical_layers` to the new public API)
- Test: `tests/modules/test_series_and_covers.py`

**Interfaces:**
- Consumes: `resolution.py::projective_cover(M) -> (Q0, d0, gens)`,
  `injective.py::injective_resolution(M, 1)` (term 0 + `iota` = differential 0),
  `radtopsoc::radical, socle`, `trace/modules.py::_radical_layers` (the
  existing layer computation — MOVE its logic, keep a delegating shim).
- Produces (methods on `Module`):
  ```python
  def projective_cover_hom(self) -> ModuleHom      # P(M) ->> M
  def injective_envelope_hom(self) -> ModuleHom    # M >-> E(M)
  def radical_series(self) -> list["Module"]       # [M, rad M, rad^2 M, ...] until 0
  def socle_series(self) -> list["Module"]         # [0, soc M, soc^2 M, ...] until M
  def loewy_layers(self) -> list[dict]             # top-to-bottom composition-factor
                                                   # multiplicity dicts (the moved logic)
  def composition_factors(self) -> dict            # total multiplicities {S_v label: m}
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_series_and_covers.py
import pytest

from quiverlab import GF, Quiver

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_projective_cover_hom_is_epi_with_superfluous_kernel():
    A = _a3()
    S1 = A.simple(1)
    f = S1.projective_cover_hom()
    assert f.is_epi() and f.tgt.dim == S1.dim
    K, _ = f.kernel()
    # cover kernel lies in rad(P): K + rad P = rad P
    assert K.dimension_vector() == f.src.radical().dimension_vector()


def test_injective_envelope_hom_is_mono():
    A = _a3()
    S2 = A.simple(2)
    f = S2.injective_envelope_hom()
    assert f.is_mono() and f.src.dim == S2.dim
    assert f.tgt.socle().dimension_vector() == S2.dimension_vector()


def test_radical_series_strictly_decreases_to_zero():
    A = _a3()
    P1 = A.projective(1)
    series = P1.radical_series()
    dims = [X.dim for X in series]
    assert dims[0] == P1.dim and dims[-1] == 0
    assert all(a > b for a, b in zip(dims, dims[1:]))
    assert len(series) - 1 <= A.loewy_length()


def test_socle_series_reaches_the_module():
    A = _a3()
    P1 = A.projective(1)
    socs = P1.socle_series()
    assert socs[0].dim == 0 and socs[-1].dim == P1.dim


def test_loewy_layers_sum_to_composition_factors():
    A = _a3()
    P1 = A.projective(1)
    layers = P1.loewy_layers()
    total = {}
    for layer in layers:
        for k, v in layer.items():
            total[k] = total.get(k, 0) + v
    assert total == P1.composition_factors()
    # all simples here are 1-dimensional, so factor count == module dim
    assert sum(total.values()) == P1.dim
```

(Layer keys are the str composition-factor labels `_radical_layers`
already emits — keep that exact format when the logic moves; the trace
renderers depend on it.)

- [ ] **Step 2: Run to verify failures** — `AttributeError`s

- [ ] **Step 3: Implement** — `projective_cover_hom` wraps
  `projective_cover(M)`'s `(Q0, d0, ...)` as `ModuleHom(Q0, M, d0)`;
  `injective_envelope_hom` wraps `injective_resolution(M, 1)`'s term-0 +
  `differential(0)` (`iota: M -> E^0` per `injective.py:32-35`) as
  `ModuleHom(M, E0, iota)`; `radical_series` iterates `radical()` until
  `dim == 0`; `socle_series` iterates preimages of socles (`soc^{n+1}/soc^n
  = soc(M/soc^n)` — build via `radtopsoc.submodule` on accumulated columns);
  `loewy_layers` = the MOVED `trace/modules.py::_radical_layers` body
  (public now; `trace` keeps `_radical_layers = Module.loewy_layers`-style
  delegation so its renderers are untouched); `composition_factors` sums
  the layers. Fix the last test's final assertion to
  `sum(m * <dim of S from its label> ...) == P1.dim` with the real label
  format.

- [ ] **Step 4: Run tests + the trace neighbors** (the moved logic must not
  change report bytes):

Run: `... -m pytest tests/modules/test_series_and_covers.py tests/trace/test_report_completeness_m0729.py -q`
Expected: PASS both.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/module.py src/quiverlab/trace/modules.py tests/modules/test_series_and_covers.py
git commit -m "feat(modules): covers/envelopes as maps, radical/socle series, composition factors"
```

---

### Task 7: QPA cross-oracle battery

**Files:**
- Modify: `src/quiverlab/qpa/scripts.py` (hom/kernel dim script builders)
- Modify: `src/quiverlab/qpa/crosscheck.py` (`crosscheck_hom_glue`)
- Test: `tests/qpa/test_hom_glue_qpa.py`

**Interfaces:**
- Consumes: the QPA session (`session.run`), `scripts.module_decl` (existing
  module serializer — read its signature), QPA GAP functions
  `HomOverAlgebra(M, N)` (list of homs; `Length` = dim),
  `KernelOfWhat...` — **first probe live QPA for the exact kernel/image
  function names** (`NamesGVars()` grep, the Plan-35 precedent in
  `tests/qpa/test_products_qpa.py`): expected `Kernel`/`Image` on module
  homomorphisms per QPA Ch. 7 — use what the probe finds, and if kernels
  of homs are NOT exposed, the battery honestly covers dims of Hom only
  and says so in the test docstring.
- Produces: `tests/qpa/test_hom_glue_qpa.py` — for the zoo pairs (kA2, kA3
  with relation, `line_abc_cde` from Plan 18): `len(hom_basis(M, N))` vs
  `Length(HomOverAlgebra(...))`, and (if exposed) kernel/image dims of a
  nonzero hom vs QPA. File header: the standard
  `pytestmark = pytest.mark.skipif(session.should_skip_qpa(), ...)`.

- [ ] **Step 1: Probe live QPA** for the hom-kernel surface; write the
  battery accordingly (test code follows the existing
  `tests/qpa/test_tor_qpa.py` shape — build both sides, compare ints).
- [ ] **Step 2: Run** `... -m pytest tests/qpa/test_hom_glue_qpa.py -v`
  (the venv has the [qpa] extra installed per project memory). Expected: PASS live.
- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/qpa/ tests/qpa/test_hom_glue_qpa.py
git commit -m "test(qpa): hom-glue battery -- Hom dims (and kernel/image dims where QPA exposes them)"
```

---

### Task 8: GUI story + verification page + suite gate

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (`rad_top_soc` block gains a `series` field)
- Modify: `docs/gui/runner.py` (the byte-identical twin)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (render the Loewy stack)
- Modify: `webapp/static/app.js` (same rendering for the webapp)
- Modify: `webapp/server/i18n/en.json`, `es.json` (`block.rad_top_soc.series` label)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py` (documented re-freeze)
- Modify: `docs/verification.md` (Plan-37 oracle rows + counts)
- Test: `tests/webapp/test_loewy_series_p37.py`, extend `tests/trace/` renderer checks

**Interfaces:**
- Consumes: `Module.loewy_layers()` (Task 6), the `rad_top_soc` handler in
  `spec.py::_dispatch_module` and its runner twin (read both, keep
  shape-identical), the `module_blocks` serializer.
- Produces: the `rad_top_soc` block carries
  `"series": [ {factor-label: multiplicity}, ... ]` (top-to-bottom), both
  GUI renderers + the report print it as a stacked Loewy diagram (one row
  per layer, factors as `S_v^m`), i18n'd. Golden `module_left_a2`-family
  entries re-frozen with the additive key — documented in the
  `test_runner_delegation.py` docstring per house precedent.

- [ ] **Step 1: Write the failing cross-runner test** (unmarked, per the
  Plan-32 extras-gated ruling — follow `tests/webapp/test_module_blocks_m0729.py`'s
  header and fixtures):

```python
# tests/webapp/test_loewy_series_p37.py
"""The rad_top_soc block carries the Loewy series; both runners agree."""


def test_rad_top_soc_block_has_series(tmp_path):
    # build the Plan-26-style request for a small module (copy the
    # module_left_a2 golden request), run quiverlab.hpc.spec.run_spec,
    # and assert block["series"] is a list of dicts whose totals match
    # block-level top/socle data; then assert the docs/gui/runner.py twin
    # returns the same series (import it the way test_module_blocks_m0729
    # does).
```

(Write it concretely by copying `test_module_blocks_m0729.py`'s runner-pair
fixture — that file is the pattern for cross-runner block contracts.)

- [ ] **Step 2: Implement** — in `spec.py`'s `rad_top_soc` handler add
  `block["series"] = [dict(layer) for layer in M.loewy_layers()]`; mirror
  in `docs/gui/runner.py`; renderers: a `loewySeries(block)` helper in both
  `gui.js` files + `app.js` and a `<table class="ql-loewy">` in
  `results_html.py`'s rad_top_soc branch; i18n keys.
- [ ] **Step 3: Re-freeze goldens** — run the delegation test, see the
  byte diff is EXACTLY the additive `series` key, re-freeze, document in
  the docstring change-log.
- [ ] **Step 4: Verification page** — add Plan-37 rows (morphism/SES/End
  self-cert batteries, the QPA hom-glue battery) to the subsystem table +
  recount the audited class counts (`tests/release/test_oracle_classes.py`
  must pass).
- [ ] **Step 5: Full gate + commit**

Run: `... -m pytest tests/modules tests/webapp tests/trace tests/release -q`
(deep bucket by path is fine); Expected: green.

```bash
git add -A && git commit -m "feat(gui,report): Loewy-series display for rad_top_soc; verification page Plan-37 rows"
```

---

## Acceptance (Plan-37 definition of done)

1. `ModuleHom` + kernel/image/cokernel + SES/split/pushout/pullback +
   `direct_sum`/`is_direct_summand` + `end_algebra` + series/covers all
   public, all loudly-validated, all `oracle_selfcert`-batteried.
2. QPA hom-glue battery green live (`-m qpa`).
3. The rad_top_soc GUI/report Loewy story shipped, goldens re-frozen with
   a documented reason, EN/ES keys present.
4. `docs/verification.md` updated; `tests/release/` green; deep bucket
   green on the touched directories.
5. P39/P41/P44/P45/P47 can consume `ModuleHom` (signatures in this doc are
   the contract).
