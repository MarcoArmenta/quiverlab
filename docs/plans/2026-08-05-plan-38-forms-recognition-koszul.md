# Plan 38: C2 Forms, Roots, Recognition + Coxeter Pins + Koszulity Exposure — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The C2 pillar: Euler/Tits forms with definiteness classification,
Dynkin/Euclidean type detection, positive-root combinatorics, the recognizer
batch (`is_hereditary` … `is_gentle`), oracle pins from Marco's Coxeter
paper on the EXISTING Coxeter surface, and the metagoal-1 debt: Koszulity
(`ext_algebra`) exposed in every no-code tier.

**Architecture:** CRITICAL scoping fact from the 2026-08-05 audit:
`cartan_matrix` / `coxeter_matrix` / `coxeter_polynomial` **already exist**
(`src/quiverlab/invariants/cartan.py`, surfaced at `core/algebra.py:444-457`,
compute kinds `cartan`/`coxeter_polynomial` already dispatched in
`hpc/spec.py:1174-1184` with GUI checkboxes), and `invariants/spectral.py`
already has exact `spectral_radius`/`mahler_measure`. **This plan adds NO
new Coxeter machinery** — only literature pins on it. New code: forms +
type detection + roots + recognizers in `src/quiverlab/invariants/`
(`forms.py`, `dynkin_type.py`, `recognizers.py`), and two new scalar
compute kinds (`ext_algebra`, `recognizers`) wired through the standard
seven-touchpoint pattern (spec dispatch, runner twin, gui.js ×2, i18n,
ETA, `_snip`, golden).

**Tech Stack:** sympy (exact Rational linear algebra — already a dep),
the Gröbner reduction system (`resolutions_cs.build.reduction_system_of`)
for ideal-membership of length-2 paths. No floats in `src/`.

## Global Constraints

- Python is always `.venv/bin/python`; tests via
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- `tests/invariants/` collects in the **fast** bucket; QPA comparisons go
  in `tests/qpa/` (bucket = class); new webapp/gui tests are unmarked
  (extras-gated dirs carry no oracle markers — Plan-32 rule).
- New scalar kinds act on the existing algebra block ⇒ **schema stays v1**
  (Plan-26 precedent: version bumps only for new request BLOCKS).
- Every new citation key is BibTeX-verified before use (Plan-29 rule).
- All refusals loud (`QuiverlabError`), never silent defaults;
  presentation-less algebras refuse recognizers that need the quiver.
- Conventional commits; green at every commit; branch
  `plan-38-forms-recognition` off `dev`.

---

### Task 1: Coxeter oracle pins from arXiv:2606.15595 (no new machinery)

**Files:**
- Modify: `src/quiverlab/citations/references.bib` (add key `armenta_coxeter_calculus`, arXiv:2606.15595, 2026 — verify the BibTeX before committing)
- Test: `tests/invariants/test_coxeter_paper_pins.py`

**Interfaces:**
- Consumes: `A.coxeter_polynomial()` (existing, sympy `Poly` in `t`),
  `families/dynkin.py::dynkin_quiver(type_str, orientation=...)`.
- Produces: frozen literature pins; the citation key later tasks and
  reports reference.

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_coxeter_paper_pins.py
"""Worked examples of Armenta, 'The Coxeter transformation as an
automorphism of the Tamarkin-Tsygan calculus' (arXiv:2606.15595), pinned
on the existing Coxeter surface. D4 vs A4: isomorphic TT calculi,
distinguished by the Coxeter polynomial; the 8-vertex cospectral trees:
equal Coxeter polynomials, not derived equivalent (honest-scope demo:
these invariants are necessary conditions only)."""
import sympy as sp
import pytest

from quiverlab import Quiver

pytestmark = pytest.mark.oracle_literature

t = sp.Symbol("t")


def test_d4_coxeter_polynomial():
    # D4: three arrows into the central vertex 1 (paper Ex 5.7)
    A = Quiver([1, 2, 3, 4], {"a": (2, 1), "b": (3, 1), "c": (4, 1)}).algebra()
    p = A.coxeter_polynomial().as_expr()
    assert sp.expand(p - (t**4 + t**3 + t + 1)) == 0
    assert sp.factor(p) == sp.factor((t + 1) ** 2 * (t**2 - t + 1))


def test_a4_coxeter_polynomial():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra()
    p = A.coxeter_polynomial().as_expr()
    assert sp.expand(p - (t**4 + t**3 + t**2 + t + 1)) == 0   # Phi_5


def test_d4_a4_differ_by_t_squared():
    D4 = Quiver([1, 2, 3, 4], {"a": (2, 1), "b": (3, 1), "c": (4, 1)}).algebra()
    A4 = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra()
    diff = sp.expand(A4.coxeter_polynomial().as_expr()
                     - D4.coxeter_polynomial().as_expr())
    assert diff == t**2                                        # paper Ex 5.7


def _spider():
    # 8 vertices: center 1 with arms of lengths 3,1,1,1,1 (degree seq 5,2,2,1x5)
    return Quiver(list(range(1, 9)),
                  {"a1": (1, 2), "a2": (2, 3), "a3": (3, 4),
                   "b": (1, 5), "c": (1, 6), "d": (1, 7), "e": (1, 8)}).algebra()


def _double_star():
    # 8 vertices: centers 1,2 joined; three leaves on each (degree seq 4,4,1x6)
    return Quiver(list(range(1, 9)),
                  {"m": (1, 2), "a": (1, 3), "b": (1, 4), "c": (1, 5),
                   "d": (2, 6), "e": (2, 7), "f": (2, 8)}).algebra()


def test_cospectral_trees_share_coxeter_polynomial():
    ps = _spider().coxeter_polynomial().as_expr()
    pd = _double_star().coxeter_polynomial().as_expr()
    target = sp.expand((t + 1) ** 4 * (t**4 - 3 * t**3 + t**2 - 3 * t + 1))
    assert sp.expand(ps - target) == 0 and sp.expand(pd - target) == 0
```

- [ ] **Step 2: Run to verify** — these should PASS immediately if the
  existing surface is correct and the tree encodings match the paper; a
  FAILURE here means the encoding is wrong (recheck the degree sequences
  5,2,2,1,1,1,1,1 and 4,4,1,1,1,1,1,1 — orientation is irrelevant for the
  char poly, tree shape is not). Fix the QUIVERS, never the pins.

Run: `... -m pytest tests/invariants/test_coxeter_paper_pins.py -v`

- [ ] **Step 3: Add the citation** to `references.bib` (author Armenta,
  title as above, eprint 2606.15595, year 2026); run
  `... -m pytest tests/trace/test_references.py -q` to confirm the bib
  still parses.

- [ ] **Step 4: Commit**

```bash
git add tests/invariants/test_coxeter_paper_pins.py src/quiverlab/citations/references.bib
git commit -m "test(invariants): pin arXiv:2606.15595 worked examples on the existing Coxeter surface"
```

---

### Task 2: Euler + Tits forms with definiteness classification

**Files:**
- Create: `src/quiverlab/invariants/forms.py`
- Modify: `src/quiverlab/invariants/__init__.py`, `src/quiverlab/core/algebra.py` (methods `euler_form`, `tits_form`, `form_type`)
- Test: `tests/invariants/test_forms.py`

**Interfaces:**
- Consumes: `invariants/cartan.py::cartan_matrix(A)` (int rows), sympy.
- Produces:
  ```python
  def euler_form_matrix(A) -> sympy.Matrix     # (C^{-1})^T, exact Rational;
                                               # QuiverlabError if det C == 0
  def euler_form(A, d, e) -> sympy.Rational    # d * E * e^T (d, e int vectors
                                               # in vertex order)
  def tits_matrix(A) -> sympy.Matrix           # E + E^T (the symmetrization)
  def tits_form(A, d) -> sympy.Rational        # q(d) = euler_form(A, d, d)
  def form_type(A) -> str                      # "finite" | "tame" | "wild":
      # positive definite / positive semidefinite (not definite) / neither,
      # decided EXACTLY on tits_matrix via sympy over Rationals.
      # Docstring states the rep-type meaning is proven for HEREDITARY
      # algebras (Gabriel / Donovan-Freislich / Nazarova); for others it is
      # just the form's signature.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_forms.py
"""Euler/Tits forms off the Cartan matrix. Literature: Gabriel's theorem
(finite type <=> Dynkin <=> positive definite Tits form), ASS Ch. VII."""
import sympy as sp
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.forms import (euler_form, euler_form_matrix,
                                        form_type, tits_form)

pytestmark = pytest.mark.oracle_literature


def _kronecker(m):
    arrows = {f"a{i}": (1, 2) for i in range(m)}
    return Quiver([1, 2], arrows).algebra()


def test_hereditary_euler_form_is_the_arrow_formula():
    # <d,e> = sum d_v e_v - sum_{a: s->t} d_s e_t   (hereditary)
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra()
    d, e = [1, 2, 1], [3, 1, 2]
    expected = (1*3 + 2*1 + 1*2) - (1*1 + 2*2)
    assert euler_form(A, d, e) == expected


def test_euler_form_computes_hom_minus_ext_dimension():
    # homological meaning, pinned on a hereditary example:
    # <dim M, dim N> = dim Hom(M,N) - dim Ext^1(M,N)
    A = Quiver([1, 2], {"a": (1, 2)}).algebra()
    S1, S2 = A.simple(1), A.simple(2)
    d1 = [1, 0]; d2 = [0, 1]
    assert euler_form(A, d1, d2) == A.hom(S1, S2) - A.ext(S1, S2, 1)
    assert euler_form(A, d2, d1) == A.hom(S2, S1) - A.ext(S2, S1, 1)


def test_form_type_dynkin_euclidean_wild():
    assert form_type(Quiver([1, 2], {"a": (1, 2)}).algebra()) == "finite"   # A2
    assert form_type(_kronecker(2)) == "tame"                               # ~A1
    assert form_type(_kronecker(3)) == "wild"                               # 3-Kronecker


def test_tits_form_values():
    # A2: q(d) = d1^2 + d2^2 - d1 d2; q(1,1) = 1 (a root!)
    A = Quiver([1, 2], {"a": (1, 2)}).algebra()
    assert tits_form(A, [1, 1]) == 1
    # Kronecker: q(1,1) = 0 (the isotropic imaginary root)
    assert tits_form(_kronecker(2), [1, 1]) == 0


def test_singular_cartan_refused():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"])   # k[x]/(x^2)
    with pytest.raises(QuiverlabError, match="[Cc]artan"):
        euler_form_matrix(A)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: forms`

- [ ] **Step 3: Implement `src/quiverlab/invariants/forms.py`**

```python
"""Euler bilinear form, Tits quadratic form, definiteness classification
(Plan 38 / C2). All exact: sympy Rational off the integer Cartan matrix.

The Euler form matrix is E = (C^{-1})^T so <d, e> = d E e^T; for finite
global dimension <dim M, dim N> = sum (-1)^i dim Ext^i(M, N). The Tits
form is q(d) = <d, d>; its symmetrized matrix is E + E^T. The
finite/tame/wild reading of definiteness is a THEOREM only for hereditary
algebras (Gabriel; Donovan-Freislich/Nazarova) -- form_type computes the
signature for any invertible Cartan and the docstring says exactly this."""
from __future__ import annotations

import sympy as sp

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.cartan import cartan_matrix


def euler_form_matrix(A) -> sp.Matrix:
    C = sp.Matrix(cartan_matrix(A))
    if C.det() == 0:
        raise QuiverlabError(
            "Euler form needs an invertible Cartan matrix (finite global "
            "dimension); det C = 0 here.")
    return C.inv().T


def euler_form(A, d, e):
    E = euler_form_matrix(A)
    dv, ev = sp.Matrix(1, E.rows, list(d)), sp.Matrix(E.rows, 1, list(e))
    return sp.nsimplify((dv * E * ev)[0, 0])


def tits_matrix(A) -> sp.Matrix:
    E = euler_form_matrix(A)
    return E + E.T


def tits_form(A, d):
    return euler_form(A, d, d)


def form_type(A) -> str:
    """'finite' / 'tame' / 'wild' by exact definiteness of the Tits form
    (positive definite / semidefinite-not-definite / neither). The
    representation-type meaning is proven for hereditary algebras only."""
    Q = tits_matrix(A)
    if Q.is_positive_definite:
        return "finite"
    if Q.is_positive_semidefinite:
        return "tame"
    return "wild"
```

Surface as `Algebra.euler_form(d, e)`, `Algebra.tits_form(d)`,
`Algebra.form_type()` beside `coxeter_polynomial` at
`core/algebra.py:444-457` (same lazy-import style), and re-export the four
functions from `invariants/__init__.py`.
**Adjust to reality:** confirm `sympy.Matrix.is_positive_semidefinite`
exists in the pinned sympy and is exact over Rationals (it is, via
charpoly/minors); if the attribute is a method in this version, call it.

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/ src/quiverlab/core/algebra.py tests/invariants/test_forms.py
git commit -m "feat(invariants): Euler/Tits forms + exact definiteness classification"
```

---

### Task 3: Dynkin / Euclidean type detection

**Files:**
- Create: `src/quiverlab/invariants/dynkin_type.py`
- Modify: `src/quiverlab/invariants/__init__.py`, `core/algebra.py` (`Algebra.dynkin_type()`), `src/quiverlab/combinat/quiver.py` (`Quiver.is_connected()`)
- Test: `tests/invariants/test_dynkin_type.py`

**Interfaces:**
- Consumes: `Quiver.vertices`, `Quiver.arrows`;
  `families/dynkin.py::dynkin_quiver(type_str, orientation=...)` (the
  generator — round-trip oracle).
- Produces:
  ```python
  def dynkin_type(quiver) -> tuple | None
      # ("A", n) | ("D", n) | ("E", 6|7|8) | ("~A", n) | ("~D", n) | ("~E", 6|7|8)
      # | None (anything else). Orientation-blind (underlying multigraph).
      # ("~A", 1) is the double edge (Kronecker).
  def is_connected(quiver) -> bool
  ```
  Algorithm: undirected multigraph; count E vs V: tree (E = V-1, no
  multi-edge) → ADE by degree signature (path → A; exactly one degree-3
  vertex → arm lengths sorted: (1,1,n-2)→D, (1,2,2)→E6, (1,2,3)→E7,
  (1,2,4)→E8); E = V with single cycle covering all vertices and no
  multi-edge → ~A_{V-1}; V=2 with a double edge → ~A1; trees with two
  degree-3 vertices or one degree-4 → ~D_n / ~E_n by the affine degree
  signatures (~D: two deg-3 vertices each with two leaf-arms, or D4-star
  with 4 leaves = ~D4 has a degree-4 center; ~E6: center degree 3, arms
  (2,2,2); ~E7: arms (1,3,3); ~E8: arms (1,2,5)). Everything else → None.

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_dynkin_type.py
"""Type detection, round-tripped against the families/dynkin generators
and cross-checked against the Tits-form signature (Gabriel)."""
import pytest

from quiverlab import Quiver
from quiverlab.families.dynkin import dynkin_quiver
from quiverlab.invariants.dynkin_type import dynkin_type, is_connected
from quiverlab.invariants.forms import form_type

pytestmark = pytest.mark.oracle_crossengine


FINITE = [("A", n) for n in range(1, 9)] + [("D", n) for n in range(4, 9)] \
         + [("E", 6), ("E", 7), ("E", 8)]
AFFINE = [("~A", n) for n in range(2, 7)] + [("~D", n) for n in range(4, 8)] \
         + [("~E", 6), ("~E", 7), ("~E", 8)]


@pytest.mark.parametrize("typ,n", FINITE, ids=[f"{t}{n}" for t, n in FINITE])
def test_roundtrip_finite(typ, n):
    Q = dynkin_quiver(f"{typ}{n}")
    assert dynkin_type(Q) == (typ, n)


@pytest.mark.parametrize("typ,n", AFFINE, ids=[f"{t}{n}" for t, n in AFFINE])
def test_roundtrip_affine(typ, n):
    Q = dynkin_quiver(f"{typ.replace('~', '~')}{n}")     # generator syntax "~A3"
    assert dynkin_type(Q) == (typ, n)


def test_kronecker_is_affine_a1():
    Q = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)})
    assert dynkin_type(Q) == ("~A", 1)


def test_agreement_with_tits_signature():
    # Dynkin <=> finite form; affine <=> tame form (hereditary, Gabriel)
    for typ, n in [("A", 4), ("D", 5), ("E", 6)]:
        A = dynkin_quiver(f"{typ}{n}").algebra()
        assert dynkin_type(A.quiver) is not None
        assert dynkin_type(A.quiver)[0] in ("A", "D", "E")
        assert form_type(A) == "finite"


def test_wild_star_is_none():
    Q = Quiver([1, 2, 3, 4, 5, 6],
               {"a": (2, 1), "b": (3, 1), "c": (4, 1), "d": (5, 1), "e": (6, 1)})
    assert dynkin_type(Q) is None          # 5-star: not ADE, not affine


def test_disconnected_is_none_and_flagged():
    Q = Quiver([1, 2, 3], {"a": (1, 2)})
    assert is_connected(Q) is False
    assert dynkin_type(Q) is None
```

**Adjust to reality:** check `families/dynkin.py::_TYPE` for the exact
affine syntax it generates (`"~A3"` per the regex `^(~|t)?([ADE])(\d+)$`)
and what index convention the generator uses for affine types (vertex
count = n+1); make the round-trip use the generator's own convention.

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** `dynkin_type.py` per the algorithm block in
  Interfaces (pure combinatorics, ~120 lines: build adjacency lists with
  multi-edge counts; helpers `_degrees`, `_arm_lengths(center)`;
  the classification table exactly as specified). `is_connected` = BFS on
  the undirected graph; also add the thin method
  `Quiver.is_connected()` delegating to it.

- [ ] **Step 4: Run tests** — Expected: PASS (all ~30 parametrized ids)

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/dynkin_type.py src/quiverlab/invariants/__init__.py \
        src/quiverlab/core/algebra.py src/quiverlab/combinat/quiver.py tests/invariants/test_dynkin_type.py
git commit -m "feat(invariants): orientation-blind Dynkin/Euclidean type detection + is_connected"
```

---

### Task 4: positive roots (hereditary Dynkin)

**Files:**
- Create: `src/quiverlab/invariants/roots.py`
- Modify: `invariants/__init__.py`, `core/algebra.py` (`Algebra.positive_roots()`)
- Test: `tests/invariants/test_roots.py`

**Interfaces:**
- Consumes: `tits_matrix`/`tits_form` (Task 2), `dynkin_type` (Task 3).
- Produces:
  ```python
  def positive_roots(A) -> list[tuple]
      # all positive roots of the Tits form as dimension vectors (vertex
      # order), for hereditary A with dynkin_type != None (finite type).
      # QuiverlabError otherwise ("positive roots enumerated only for
      # Dynkin hereditary type; affine/wild have infinitely many").
      # Algorithm: reflection closure of the simple roots under the Weyl
      # reflections s_i(d) = d - (d, e_i) e_i  with (,-) the symmetrized
      # form; equivalently BFS on {d : q(d) = 1, d >= 0}. Both finite for
      # Dynkin -- implement the reflection closure (exact, no search cap).
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_roots.py
"""Root counts pinned to the classical tables (Bourbaki; ASS VII):
A_n: n(n+1)/2, D_n: n(n-1), E6: 36, E7: 63, E8: 120. Gabriel: for a
Dynkin quiver these are exactly the dimension vectors of the
indecomposables."""
import pytest

from quiverlab.families.dynkin import dynkin_quiver
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.roots import positive_roots

pytestmark = pytest.mark.oracle_literature


@pytest.mark.parametrize("typ,count", [
    ("A2", 3), ("A3", 6), ("A5", 15), ("D4", 12), ("D5", 20), ("E6", 36),
])
def test_root_counts(typ, count):
    A = dynkin_quiver(typ).algebra()
    roots = positive_roots(A)
    assert len(roots) == count
    assert len(set(roots)) == count            # no duplicates
    from quiverlab.invariants.forms import tits_form
    assert all(tits_form(A, list(r)) == 1 for r in roots)


def test_a2_roots_explicit():
    A = dynkin_quiver("A2").algebra()
    assert sorted(positive_roots(A)) == [(0, 1), (1, 0), (1, 1)]


def test_highest_root_d4():
    A = dynkin_quiver("D4").algebra()
    assert max(positive_roots(A), key=sum) == (2, 1, 1, 1) or \
           sorted(max(positive_roots(A), key=sum)) == [1, 1, 1, 2]
    # (vertex order depends on the generator -- assert the multiset)


def test_non_dynkin_refused():
    from quiverlab import Quiver
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra()
    with pytest.raises(QuiverlabError, match="Dynkin"):
        positive_roots(K)
```

(Sharpen `test_highest_root_d4` at implementation time to the generator's
actual vertex order — assert the sorted tuple only, drop the `or`.)

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** — reflection closure: start from the unit
  vectors; repeatedly apply all simple reflections
  `s_i(d)_j = d_j` for `j != i`, `s_i(d)_i = -d_i + sum_{j~i} m_ij d_j`
  (`m_ij` = number of edges between i and j in the underlying graph);
  keep vectors with all entries `>= 0`, stop at closure (finite for
  Dynkin). Guard: `is_hereditary` not yet available until Task 5 —
  require `A.relations` empty + `dynkin_type(A.quiver)` in the finite
  list, refusing loudly otherwise (then Task 5's `is_hereditary` replaces
  the inline check — leave a one-line TODO-free direct call swap in
  Task 5 Step 3).

- [ ] **Step 4: Run tests** — Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/roots.py src/quiverlab/invariants/__init__.py \
        src/quiverlab/core/algebra.py tests/invariants/test_roots.py
git commit -m "feat(invariants): positive-root enumeration for Dynkin hereditary type"
```

---

### Task 5: the recognizer batch

**Files:**
- Create: `src/quiverlab/invariants/recognizers.py`
- Modify: `invariants/__init__.py`, `core/algebra.py` (one thin method per recognizer)
- Test: `tests/invariants/test_recognizers.py`, `tests/qpa/test_recognizers_qpa.py`

**Interfaces:**
- Consumes: `A.quiver`, `A.relations`, `A.loewy_length()` (exact,
  presentation-less-safe), `combinat/relations.py::is_monomial`,
  `combinat/quiver.py::is_acyclic`, `is_connected` (Task 3),
  `resolutions_cs.build.reduction_system_of(A)` for ideal membership of
  length-2 paths (read its API first: the reduction system's normal form /
  `leading_words()`; a length-2 path `a*b` is in I iff its normal form is 0
  — find the exact call in `resolutions_cs/build.py` and use it; if no
  public normal-form exists, reduce via the Gröbner rules directly).
- Produces (each takes `A`, returns `bool`, raises `QuiverlabError` with
  "needs the quiver presentation" on presentation-less input where the
  quiver is required — mirror `builders._require_provenance`):
  ```python
  is_semisimple(A)            # loewy_length == 1        (presentation-less OK)
  is_radical_square_zero(A)   # loewy_length <= 2        (presentation-less OK)
  is_hereditary(A)            # presented: relations == [] and quiver acyclic
                              # (admissible I != 0 => never hereditary)
  is_basic(A)                 # presented kQ/I: always True (docstring: the
                              # e_v are a complete set of primitive orthogonal
                              # idempotents with pairwise non-isomorphic
                              # projectives); presentation-less: loud refusal
  is_nakayama(A)              # underlying quiver is a linear A_n or a single
                              # oriented cycle (every vertex in-degree <= 1
                              # and out-degree <= 1)
  is_special_biserial(A)      # (1) every vertex has <= 2 in-arrows and <= 2
                              # out-arrows; (2) for every arrow b: at most one
                              # arrow c with b*c not in I, at most one arrow a
                              # with a*b not in I
  is_string(A)                # special biserial + I monomial
  is_gentle(A)                # string + I generated by length-2 paths + the
                              # dual conditions: for every arrow b, at most
                              # one c with b*c in I, at most one a with a*b in I
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_recognizers.py
"""Recognizer batch pinned on textbook examples (ASS; Butler-Ringel;
Assem-Skowronski gentle papers)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.invariants import recognizers as rec

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra()


def _gentle_a3():
    # 1 --a--> 2 --b--> 3 with a*b = 0: the standard gentle example
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"])


def _sb_not_string():
    # special biserial with a non-monomial relation: commutative square
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"])


def test_semisimple_and_radsq():
    S = Quiver([1, 2], {}).algebra()
    assert rec.is_semisimple(S) is True
    A = _kA2()
    assert rec.is_semisimple(A) is False
    assert rec.is_radical_square_zero(A) is True
    B = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x*x"], field=GF(5))
    assert rec.is_radical_square_zero(B) is False


def test_hereditary():
    assert rec.is_hereditary(_kA2()) is True
    assert rec.is_hereditary(_gentle_a3()) is False        # relations present
    loop = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(5))
    assert rec.is_hereditary(loop) is False                # cycle


def test_nakayama():
    assert rec.is_nakayama(_kA2()) is True
    cyc = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b*a", "b*a*b"], field=GF(5))
    assert rec.is_nakayama(cyc) is True
    assert rec.is_nakayama(_sb_not_string()) is False      # vertex 1 out-deg 2


def test_gentle_string_special_biserial_hierarchy():
    G = _gentle_a3()
    assert rec.is_gentle(G) and rec.is_string(G) and rec.is_special_biserial(G)
    SB = _sb_not_string()
    assert rec.is_special_biserial(SB) is True
    assert rec.is_string(SB) is False                      # binomial relation
    assert rec.is_gentle(SB) is False
    # NOT special biserial: three arrows out of one vertex
    W = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3), "c": (1, 4)}).algebra()
    assert rec.is_special_biserial(W) is False


def test_hereditary_path_algebra_is_gentle_iff_a_n():
    # kA_n with no relations is gentle; the 3-star is not (condition 1)
    assert rec.is_gentle(_kA2()) is True
```

```python
# tests/qpa/test_recognizers_qpa.py
"""Live QPA crosscheck: IsGentleAlgebra / IsSpecialBiserialAlgebra."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.invariants import recognizers as rec
from quiverlab.qpa import scripts, session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


CASES = [
    ("gentle_a3", lambda: Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))),
    ("square", lambda: Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4),
                                             "c": (1, 3), "d": (3, 4)}).algebra(
        relations=["a*b - c*d"], field=GF(5))),
    ("three_star", lambda: Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3),
                                                 "c": (1, 4)}).algebra(field=GF(5))),
]


@pytest.mark.parametrize("name,build", CASES, ids=[c[0] for c in CASES])
def test_gentle_and_sb_match_qpa(name, build):
    A = build()
    decl = scripts.quiver_and_algebra_script(A)      # read the real name/signature
    ours_sb = rec.is_special_biserial(A)
    ours_g = rec.is_gentle(A)
    qpa_sb = bool(session.run(decl + "\nIsSpecialBiserialAlgebra(A);"))
    qpa_g = bool(session.run(decl + "\nIsGentleAlgebra(A);"))
    assert (ours_sb, ours_g) == (qpa_sb, qpa_g)
```

**Adjust to reality:** the GAP variable name the script builder binds
(`A` assumed — read `scripts.quiver_and_algebra_script` and use its actual
output variable), and how `session.run`'s GAP booleans coerce.

- [ ] **Step 2: Run to verify failures** (fast file only)

- [ ] **Step 3: Implement `recognizers.py`** per the Interfaces contract
  (~150 lines; the only subtle part is length-2 ideal membership for
  special-biserial/gentle — one helper `_len2_in_ideal(A, a, b) -> bool`
  built on the reduction system, used by both conditions). Swap Task 4's
  inline hereditary check to `is_hereditary`. Add the `Algebra` methods.

- [ ] **Step 4: Run tests** — fast file + the QPA file live
  (`... -m pytest tests/invariants/test_recognizers.py tests/qpa/test_recognizers_qpa.py -v`).
  Expected: PASS both (venv has [qpa]).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/invariants/recognizers.py src/quiverlab/invariants/__init__.py \
        src/quiverlab/core/algebra.py tests/invariants/ tests/qpa/test_recognizers_qpa.py
git commit -m "feat(invariants): recognizer batch (semisimple..gentle) + live QPA crosscheck"
```

---

### Task 6: Koszulity + recognizer no-code exposure (the GUI story)

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch`: kinds `ext_algebra`, `recognizers`; `_snip` recipe map ~line 1713)
- Modify: `docs/gui/runner.py` (twin handlers + ETA `"scalars"` entries ~line 865)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkboxes ~line 71-73, `S.ids` ~142, push list ~651, block renderers)
- Modify: `webapp/static/app.js` (webapp block renderer)
- Modify: `webapp/server/i18n/en.json` + `es.json`
- Modify: `src/quiverlab/trace/results_html.py` (report renderer branch)
- Modify: `tests/webapp/_runner_goldens.json` + `test_runner_delegation.py` (NEW fixtures, existing entries byte-identical)
- Test: `tests/webapp/test_koszul_exposure_p38.py`, `tests/hpc/test_p38_kinds.py`

**Interfaces:**
- Consumes: `A.ext_algebra(top)` → `YonedaPresentation` with `.koszul`
  (True/False/None), `.koszul_obstruction`, `._koszul_reason`,
  `.certified_through_degree`, `.graded_dims_through(d)`,
  `.generators_by_degree`, `.relations_by_degree`; the recognizers (Task 5),
  `dynkin_type`/`form_type` (Tasks 2–3); the handler pattern of
  `spec.py:1174-1184` (`cartan`/`coxeter_polynomial` — copy its shape:
  payload + `"references"` + `"citations": _citation_pairs(keys)`).
- Produces two new scalar kinds (schema v1 — no new request block):
  - `ext_algebra` block: `{"koszul": true|false|null, "koszul_reason": str,
    "obstruction": [deg, reason]|null, "certified_through_degree": int,
    "graded_dims": [ints], "generators_by_degree": {deg: count},
    "relations_by_degree": {deg: count}, "latex": ..., "references":
    ["priddy", "froberg_koszul", "polishchuk_positselski"], "citations": ...}`
    with a `top` default of 6 taken from the request's existing top field
    the way other kinds read it (read how `_dispatch` receives `top`).
  - `recognizers` block: `{"flags": {"is_semisimple": ..., ...all eight...,
    "is_selfinjective": ..., "is_symmetric": ...}, "dynkin_type": "A5"|null,
    "form_type": "finite"|"tame"|"wild"|null (null when Cartan singular),
    "references": ["assem_book"], "citations": ...}` — flags that REFUSE on
    this input (e.g. presentation-less) are reported as
    `{"error": "<the loud message>"}` per-flag, never a silent False
    (the Plan-30 τ-block precedent for honest per-entry errors).
- Both runners byte-identical on these blocks; three-valued `koszul`
  rendered as "Koszul / not Koszul (obstruction at degree d) / undecided
  through degree N" in the GUI + report + EN/ES.

- [ ] **Step 1: Write the failing cross-runner test** (concrete pattern:
  copy `tests/webapp/test_module_blocks_m0729.py`'s runner-pair fixture)

```python
# tests/webapp/test_koszul_exposure_p38.py
"""ext_algebra + recognizers kinds: served by hpc.spec, mirrored by the
Pyodide twin, three-valued verdict rendered honestly."""
# (unmarked file -- extras-gated dir)


def test_ext_algebra_block_shape(tmp_path):
    # request: the koszul poster child k<x,y>/(x^2, y^2, x*y+y*x) over GF(7)
    # via the standard request-dict fixture; assert:
    #   block["koszul"] is True
    #   block["graded_dims"][0:3] == [1, 2, 3]     # quadratic dual of exterior = symmetric
    #   "priddy" in block["references"]
    ...


def test_recognizers_block_flags(tmp_path):
    # request: gentle A3 (a*b = 0). assert flags is_gentle/is_string/
    # is_special_biserial all True, is_hereditary False,
    # dynkin_type == "A3", form_type == "finite"
    ...


def test_twin_parity(tmp_path):
    # run the same requests through docs/gui/runner.py's dispatch the way
    # test_module_blocks_m0729 does; assert json.dumps(sort_keys=True)
    # equality on both new blocks.
    ...
```

(Write the `...` bodies concretely by copying the m0729 fixture; the three
assertions listed are the contract. `graded_dims` for the exterior-algebra
case: E(Λ) = symmetric algebra, dims 1, 2, 3, 4, ... — a literature-true
pin worth a comment.)

- [ ] **Step 2: Implement the two handlers in `spec.py`**, then the twin,
  then the GUI wiring (checkbox ids `qlgui-ext_algebra`,
  `qlgui-recognizers`; push-list entries; `S.ids`), ETA scalars
  (`"ext_algebra": 2.0, "recognizers": 0.1`), i18n keys
  (`inv.ext_algebra`, `inv.recognizers`, `block.ext_algebra.title`,
  `block.ext_algebra.koszul_yes/no/undecided`, `block.recognizers.title` —
  EN and ES), `_snip` recipes, `results_html.py` branch.

- [ ] **Step 3: Add the two golden fixtures** to `_runner_goldens.json`
  (e.g. `ext_algebra_exterior_gf7`, `recognizers_gentle_a3`) and note the
  addition in the `test_runner_delegation.py` docstring change-log.
  Existing golden entries must remain byte-identical — verify by running
  the delegation test BEFORE adding the fixtures.

- [ ] **Step 4: Run the gates**

Run: `... -m pytest tests/webapp/test_koszul_exposure_p38.py tests/webapp/test_runner_delegation.py tests/hpc -q`
Expected: PASS; then `... -m pytest tests/gui -q` (runner twin's own suite).

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): ext_algebra + recognizers compute kinds -- Koszulity clickable end-to-end (metagoal-1 debt)"
```

---

### Task 7: verification page, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Test: existing release gates

- [ ] **Step 1:** Verification page: add Plan-38 rows to the subsystem
  table (`invariants/` row gains forms/type/roots/recognizers with their
  oracle classes; the `qpa/` row's count grows), add the arXiv:2606.15595
  pins to the Class-1 literature list, recount the class table
  (`tests/release/test_oracle_classes.py` drives the numbers — run
  collection, paste, re-run to green). README: one line in the features
  list ("type detection, Euler/Tits forms, positive roots, gentle/string
  recognizers; Koszulity in the GUI").
- [ ] **Step 2:** Full gate:
  `... -m pytest -q -m fast` green;
  `... -m pytest tests/qpa -q -m qpa` green;
  `... -m pytest tests/release -q` green.
- [ ] **Step 3: Commit**

```bash
git add docs/verification.md README.md
git commit -m "docs(verification): Plan-38 oracle rows + recounted classes"
```

---

## Acceptance (Plan-38 definition of done)

1. Euler/Tits forms, `form_type`, `dynkin_type`, `positive_roots`, and the
   eight recognizers public on `Algebra`, exact, loud on refusal.
2. arXiv:2606.15595 pins green on the existing Coxeter surface; citation
   key verified.
3. QPA recognizer crosscheck green live; round-trip type detection green
   against every `families/dynkin.py` generator.
4. `ext_algebra` + `recognizers` clickable end-to-end (GUI canvas → block →
   report) in EN+ES, both runners byte-identical, goldens added with a
   documented change-log entry, schema still v1.
5. Verification page recounted; fast + qpa + release suites green.
