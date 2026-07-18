# quiverlab Plan 01 — Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A pip-installable package where a user builds a monomial quiver algebra over exact ℂ or GF(q) in three lines and gets certified Hochschild (co)homology dimensions from the normalized bar complex.

**Architecture:** Exact-arithmetic Domain layer (ℚ via `Fraction`, GF(p), GF(p^n), exact-ℂ via sympy algebraic fields) under a structure-constant `Algebra` (hanlab's `(m, T, unit)` currency, domain-generic). A monomial `kQ/I` front-end certifies finite-dimensionality with an automaton and emits `Algebra`. The normalized bar complex computes HH^n/HH_n over any Domain. Everything later (hanlab fast kernels, Gröbner, Chouhy–Solotar) slots behind these interfaces.

**Tech Stack:** Python ≥ 3.10, sympy (only hard dep this phase), pytest, setuptools src-layout.

## Global Constraints

- Repo root: `/Users/marco/Desktop/quiverlab`. All paths below are relative to it.
- Python ≥ 3.10. Hard dependency this phase: `sympy>=1.12` only. numpy arrives in Plan 02, matplotlib in Plan 07.
- **No floats anywhere in `src/`** — enforced by the AST gate test (Task 1). No float or complex literals, no `float()` calls. Decimal literals in user strings raise `ExactnessError`.
- Every quiverlab exception message states the problem AND a fix-it hint (spec §7).
- Path composition is **left-to-right**: `a*b` means first `a` then `b`, requires `target(a) == source(b)` (Assem–Simson–Skowroński convention). Document in every relevant docstring.
- The read-only banks under `/Users/marco/Desktop/HomologicalNetworks/` are never touched in this plan (not even read — Plan 02 reads hanlab).
- Run tests with `python -m pytest -q` from the repo root; must be green at every commit.
- Commit messages: conventional prefixes (`feat:`, `test:`, `chore:`); append your harness's standard co-author trailer.
- License: MIT, copyright Marco Armenta 2026.

---

### Task 1: Package skeleton + float-ban AST gate

**Files:**
- Create: `pyproject.toml`, `LICENSE`, `README.md`, `src/quiverlab/__init__.py`, `tests/test_no_floats.py`, `tests/__init__.py` (empty)

**Interfaces:**
- Consumes: nothing.
- Produces: importable package `quiverlab` with `__version__: str`; the float gate every later task must pass.

- [ ] **Step 1: Write the failing test**

`tests/test_no_floats.py`:

```python
"""Structural enforcement of the loud-exactness contract (spec D3, §4.1):
no float/complex literals and no float() calls anywhere under src/."""
import ast
import pathlib

SRC = pathlib.Path(__file__).resolve().parent.parent / "src" / "quiverlab"


def _violations(path: pathlib.Path) -> list[str]:
    tree = ast.parse(path.read_text(), filename=str(path))
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, (float, complex)):
            out.append(f"{path.name}:{node.lineno}: literal {node.value!r}")
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "float"
        ):
            out.append(f"{path.name}:{node.lineno}: float() call")
    return out


def test_package_importable():
    import quiverlab
    assert isinstance(quiverlab.__version__, str)


def test_no_float_literals_or_calls_in_src():
    assert SRC.is_dir(), "src/quiverlab missing"
    bad = [v for f in SRC.rglob("*.py") for v in _violations(f)]
    assert bad == [], "floats are banned in quiverlab core:\n" + "\n".join(bad)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_no_floats.py -q`
Expected: FAIL / error — `ModuleNotFoundError: No module named 'quiverlab'`.

- [ ] **Step 3: Write minimal implementation**

`pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "quiverlab"
version = "0.1.0.dev0"
description = "Quivers with relations and Hochschild theory, exactly, for algebraists"
readme = "README.md"
requires-python = ">=3.10"
license = { text = "MIT" }
authors = [{ name = "Marco Armenta", email = "drmarcoarmenta@gmail.com" }]
dependencies = ["sympy>=1.12"]

[project.optional-dependencies]
dev = ["pytest>=8"]

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

`LICENSE`: the standard MIT license text with the line `Copyright (c) 2026 Marco Armenta`.

`README.md`:

```markdown
# quiverlab

Quivers with relations and Hochschild theory, exactly, for algebraists.
Work in progress — see `docs/specs/` and `docs/plans/`.
```

`src/quiverlab/__init__.py`:

```python
"""quiverlab: quivers with relations and Hochschild theory, exactly."""

__version__ = "0.1.0.dev0"
```

- [ ] **Step 4: Install editable and verify tests pass**

Run: `python -m pip install -e '.[dev]' && python -m pytest -q`
Expected: `2 passed`.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml LICENSE README.md src tests
git commit -m "chore: package skeleton with float-ban AST gate"
```

---

### Task 2: Error hierarchy

**Files:**
- Create: `src/quiverlab/errors.py`, `tests/test_errors.py`

**Interfaces:**
- Produces: `QuiverlabError(Exception)`; subclasses `ExactnessError`, `FieldError`, `RelationError`, `AdmissibilityError`, `NotFiniteDimensionalError`, `DepthLimitError` — all importable from `quiverlab.errors` and re-exported from `quiverlab`.

- [ ] **Step 1: Write the failing test**

`tests/test_errors.py`:

```python
import pytest
import quiverlab
from quiverlab.errors import (
    QuiverlabError, ExactnessError, FieldError, RelationError,
    AdmissibilityError, NotFiniteDimensionalError, DepthLimitError,
)


def test_hierarchy():
    for cls in (ExactnessError, FieldError, RelationError, AdmissibilityError,
                NotFiniteDimensionalError, DepthLimitError):
        assert issubclass(cls, QuiverlabError)


def test_reexported_from_package():
    assert quiverlab.ExactnessError is ExactnessError


def test_messages_carry_hint():
    err = ExactnessError("0.5 is a float", hint="write '1/2' instead")
    assert "0.5 is a float" in str(err) and "write '1/2' instead" in str(err)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_errors.py -q`
Expected: FAIL — `ImportError` (no `quiverlab.errors`).

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/errors.py`:

```python
"""quiverlab exceptions (spec §7). Every message states the problem and a fix-it hint."""


class QuiverlabError(Exception):
    def __init__(self, message: str, hint: str | None = None):
        self.hint = hint
        super().__init__(message if hint is None else f"{message}  [hint: {hint}]")


class ExactnessError(QuiverlabError):
    """A float (or other non-exact input) tried to enter quiverlab."""


class FieldError(QuiverlabError):
    """Unsupported field, or an entry that does not live in the stated field."""


class RelationError(QuiverlabError):
    """A relation string is malformed, non-composable, or not parallel."""


class AdmissibilityError(QuiverlabError):
    """The ideal cannot be certified admissible."""


class NotFiniteDimensionalError(QuiverlabError):
    """The presented algebra is (or cannot be certified) finite-dimensional."""


class DepthLimitError(QuiverlabError):
    """A guard stopped a computation; the certified range is stated in the message."""
```

Append to `src/quiverlab/__init__.py`:

```python
from quiverlab.errors import (  # noqa: E402,F401
    QuiverlabError, ExactnessError, FieldError, RelationError,
    AdmissibilityError, NotFiniteDimensionalError, DepthLimitError,
)
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/errors.py src/quiverlab/__init__.py tests/test_errors.py
git commit -m "feat: QuiverlabError hierarchy with fix-it hints"
```

---

### Task 3: Domain protocol + exact rationals

**Files:**
- Create: `src/quiverlab/fields/__init__.py`, `src/quiverlab/fields/domain.py`, `src/quiverlab/fields/rationals.py`, `tests/fields/__init__.py` (empty), `tests/fields/test_rationals.py`

**Interfaces:**
- Consumes: `quiverlab.errors`.
- Produces: abstract `Domain` with methods `coerce(x)`, `zero()`, `one()`, `add(a,b)`, `sub(a,b)`, `neg(a)`, `mul(a,b)`, `inv(a)`, `is_zero(a)`, `eq(a,b)`, `to_str(a)`, attributes `name: str`, `characteristic: int`, and hooks `parse_entry(x)` (default: return x) and `make_domain(entries)` (default: return self). Singleton `QQ` (elements are `fractions.Fraction`). **Every Domain.coerce accepts:** domain elements, Python `int`, `fractions.Fraction`, and strings like `"2"`, `"-1/3"`; floats raise `ExactnessError`; decimal strings raise `ExactnessError`.

- [ ] **Step 1: Write the failing test**

`tests/fields/test_rationals.py`:

```python
from fractions import Fraction

import pytest
from quiverlab.errors import ExactnessError
from quiverlab.fields import QQ


def test_coerce_and_arithmetic():
    a = QQ.coerce("1/3")
    b = QQ.coerce(2)
    assert QQ.eq(QQ.add(a, b), Fraction(7, 3))
    assert QQ.eq(QQ.mul(a, QQ.inv(a)), QQ.one())
    assert QQ.is_zero(QQ.sub(b, QQ.coerce(Fraction(2))))
    assert QQ.characteristic == 0
    assert QQ.to_str(a) == "1/3"


def test_floats_fail_loudly():
    with pytest.raises(ExactnessError):
        QQ.coerce(0.5)
    with pytest.raises(ExactnessError):
        QQ.coerce("0.5")


def test_bool_rejected():
    with pytest.raises(ExactnessError):
        QQ.coerce(True)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/fields/test_rationals.py -q`
Expected: FAIL — `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/fields/domain.py`:

```python
"""The Domain protocol: exact field arithmetic behind one interface (spec §5, component 1)."""
from fractions import Fraction

from quiverlab.errors import ExactnessError

_FLOAT_HINT = "quiverlab is exact-only: write '1/3' or Fraction(1, 3), never 0.333"


def reject_inexact(x):
    """Loud gate for obviously non-exact Python inputs (spec D3)."""
    if isinstance(x, bool):
        raise ExactnessError(f"{x!r} is a bool, not a scalar", hint="use 0 or 1")
    if isinstance(x, (float, complex)):
        raise ExactnessError(f"{x!r} is a float", hint=_FLOAT_HINT)
    if isinstance(x, str) and any(
        ch == "." and (i + 1 < len(x) and x[i + 1].isdigit() or i > 0 and x[i - 1].isdigit())
        for i, ch in enumerate(x)
    ):
        raise ExactnessError(f"decimal literal in {x!r}", hint=_FLOAT_HINT)
    return x


class Domain:
    """Abstract exact field. Elements are plain Python objects; ops go through the domain."""

    name: str = "?"
    characteristic: int = 0

    # -- construction hooks used by Algebra builders ------------------------
    def parse_entry(self, x):
        """Pre-parse a raw user entry (overridden by CC). Default: pass through."""
        return reject_inexact(x)

    def make_domain(self, entries):
        """Return the concrete Domain for these entries (overridden by CC)."""
        return self

    # -- required field interface ------------------------------------------
    def coerce(self, x):
        raise NotImplementedError

    def zero(self):
        raise NotImplementedError

    def one(self):
        raise NotImplementedError

    def add(self, a, b):
        raise NotImplementedError

    def neg(self, a):
        raise NotImplementedError

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        raise NotImplementedError

    def inv(self, a):
        raise NotImplementedError

    def is_zero(self, a) -> bool:
        raise NotImplementedError

    def eq(self, a, b) -> bool:
        return self.is_zero(self.sub(a, b))

    def to_str(self, a) -> str:
        return str(a)

    def __repr__(self):
        return self.name


def parse_rational(x) -> Fraction:
    """int | Fraction | 'a/b' string -> Fraction, loudly exact."""
    reject_inexact(x)
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, Fraction):
        return x
    if isinstance(x, str):
        try:
            return Fraction(x.strip())
        except (ValueError, ZeroDivisionError) as exc:
            raise ExactnessError(
                f"cannot read {x!r} as an exact rational", hint="use forms like '2' or '-1/3'"
            ) from exc
    raise ExactnessError(f"cannot read {type(x).__name__} {x!r} as an exact scalar",
                         hint=_FLOAT_HINT)
```

`src/quiverlab/fields/rationals.py`:

```python
"""Exact rationals: the char-0 workhorse domain (elements: fractions.Fraction)."""
from fractions import Fraction

from quiverlab.fields.domain import Domain, parse_rational


class RationalField(Domain):
    name = "QQ"
    characteristic = 0

    def coerce(self, x):
        if isinstance(x, Fraction) and not isinstance(x, bool):
            return x
        return parse_rational(x)

    def zero(self):
        return Fraction(0)

    def one(self):
        return Fraction(1)

    def add(self, a, b):
        return a + b

    def neg(self, a):
        return -a

    def mul(self, a, b):
        return a * b

    def inv(self, a):
        return Fraction(1) / a

    def is_zero(self, a):
        return a == 0


QQ = RationalField()
```

`src/quiverlab/fields/__init__.py`:

```python
from quiverlab.fields.domain import Domain, reject_inexact, parse_rational  # noqa: F401
from quiverlab.fields.rationals import QQ  # noqa: F401
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass (float gate still green — no float literals were added).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/fields tests/fields
git commit -m "feat: Domain protocol and exact rational field QQ"
```

---

### Task 4: Generic exact linear algebra over a Domain

**Files:**
- Create: `src/quiverlab/fields/linalg.py`, `tests/fields/test_linalg.py`

**Interfaces:**
- Consumes: `Domain`.
- Produces: `rref(rows, dom) -> tuple[list[list], list[int]]` (reduced row echelon + pivot columns), `rank(rows, dom) -> int`, `nullspace(rows, dom) -> list[list]` (basis of right kernel), `solve(A, b, dom) -> list | None` (one solution of A x = b, None if inconsistent). Matrices are `list[list[element]]`, vectors `list[element]`. These four names are frozen for all later plans.

- [ ] **Step 1: Write the failing test**

`tests/fields/test_linalg.py`:

```python
from fractions import Fraction as F

from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace, rank, rref, solve


def _m(rows):
    return [[QQ.coerce(x) for x in r] for r in rows]


def test_rank_and_rref():
    A = _m([[1, 2, 3], [2, 4, 6], [1, 0, 1]])
    assert rank(A, QQ) == 2
    R, piv = rref(A, QQ)
    assert piv == [0, 1]
    assert R[0][0] == F(1) and R[1][1] == F(1)


def test_nullspace_is_kernel():
    A = _m([[1, 2, 3], [2, 4, 6], [1, 0, 1]])
    N = nullspace(A, QQ)
    assert len(N) == 1
    v = N[0]
    for row in A:
        s = QQ.zero()
        for a, x in zip(row, v):
            s = QQ.add(s, QQ.mul(a, x))
        assert QQ.is_zero(s)


def test_solve():
    A = _m([[1, 1], [1, -1]])
    b = [QQ.coerce(3), QQ.coerce(1)]
    x = solve(A, b, QQ)
    assert x == [F(2), F(1)]
    A2 = _m([[1, 1], [1, 1]])
    assert solve(A2, [QQ.coerce(0), QQ.coerce(1)], QQ) is None


def test_empty_matrix():
    assert rank([], QQ) == 0
    assert nullspace([], QQ) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/fields/test_linalg.py -q`
Expected: FAIL — `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/fields/linalg.py`:

```python
"""Exact Gaussian elimination over any Domain. Correctness first; fast GF(p)
kernels replace this behind the same names in Plan 02."""


def rref(rows, dom):
    A = [list(r) for r in rows]
    nr = len(A)
    nc = len(A[0]) if nr else 0
    pivots = []
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if not dom.is_zero(A[i][c])), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = dom.inv(A[r][c])
        A[r] = [dom.mul(inv, x) for x in A[r]]
        for i in range(nr):
            if i != r and not dom.is_zero(A[i][c]):
                f = A[i][c]
                A[i] = [dom.sub(x, dom.mul(f, y)) for x, y in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
        if r == nr:
            break
    return A, pivots


def rank(rows, dom):
    return len(rref(rows, dom)[1])


def nullspace(rows, dom):
    if not rows:
        return []
    R, pivots = rref(rows, dom)
    nc = len(rows[0])
    free = [c for c in range(nc) if c not in pivots]
    basis = []
    for fc in free:
        v = [dom.zero()] * nc
        v[fc] = dom.one()
        for r, pc in enumerate(pivots):
            v[pc] = dom.neg(R[r][fc])
        basis.append(v)
    return basis


def solve(A, b, dom):
    if not A:
        return [] if all(dom.is_zero(x) for x in b) else None
    nc = len(A[0])
    aug = [list(row) + [x] for row, x in zip(A, b)]
    R, pivots = rref(aug, dom)
    if nc in pivots:  # pivot in the augmented column: inconsistent
        return None
    x = [dom.zero()] * nc
    for r, pc in enumerate(pivots):
        x[pc] = R[r][nc]
    return x
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/fields/linalg.py tests/fields/test_linalg.py
git commit -m "feat: exact rref/rank/nullspace/solve over any Domain"
```

---

### Task 5: Prime fields GF(p)

**Files:**
- Create: `src/quiverlab/fields/primefield.py`, `tests/fields/test_primefield.py`
- Modify: `src/quiverlab/fields/__init__.py`

**Interfaces:**
- Consumes: `Domain`, `parse_rational`, `FieldError`.
- Produces: `PrimeField(p)` (elements: `int` in `0..p-1`); `PrimeField.coerce` additionally maps rationals `a/b` to `a * b^{-1} mod p` (raising `FieldError` when `p | b`). Not exported to users directly — the `GF(q)` dispatcher (Task 6) is the public door.

- [ ] **Step 1: Write the failing test**

`tests/fields/test_primefield.py`:

```python
import pytest
from quiverlab.errors import ExactnessError, FieldError
from quiverlab.fields.primefield import PrimeField


def test_arithmetic_mod_5():
    F5 = PrimeField(5)
    assert F5.characteristic == 5
    assert F5.add(3, 4) == 2
    assert F5.mul(3, 4) == 2
    assert F5.inv(3) == 2
    assert F5.is_zero(F5.sub(2, 2))
    assert F5.coerce(-1) == 4
    assert F5.coerce("1/3") == F5.mul(1, F5.inv(3))


def test_denominator_divisible_by_p():
    F5 = PrimeField(5)
    with pytest.raises(FieldError):
        F5.coerce("1/5")


def test_nonprime_rejected():
    with pytest.raises(FieldError):
        PrimeField(6)


def test_float_rejected():
    with pytest.raises(ExactnessError):
        PrimeField(5).coerce(0.2)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/fields/test_primefield.py -q`
Expected: FAIL — `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/fields/primefield.py`:

```python
"""GF(p): exact prime-field arithmetic on plain ints (spec §4.3)."""
import sympy

from quiverlab.errors import FieldError
from quiverlab.fields.domain import Domain, parse_rational


class PrimeField(Domain):
    characteristic: int

    def __init__(self, p: int):
        if not (isinstance(p, int) and not isinstance(p, bool) and sympy.isprime(p)):
            raise FieldError(f"GF({p!r}): {p!r} is not a prime",
                             hint="use GF(p) with p prime, or GF(p**n) for prime powers")
        self.p = p
        self.characteristic = p
        self.name = f"GF({p})"

    def coerce(self, x):
        if isinstance(x, int) and not isinstance(x, bool):
            return x % self.p
        q = parse_rational(x)
        if q.denominator % self.p == 0:
            raise FieldError(
                f"{q} has denominator divisible by {self.p}, so it is not in GF({self.p})",
                hint="clear denominators in your presentation or choose another characteristic",
            )
        return (q.numerator * pow(q.denominator, -1, self.p)) % self.p

    def zero(self):
        return 0

    def one(self):
        return 1 % self.p

    def add(self, a, b):
        return (a + b) % self.p

    def neg(self, a):
        return (-a) % self.p

    def mul(self, a, b):
        return (a * b) % self.p

    def inv(self, a):
        if a % self.p == 0:
            raise ZeroDivisionError(f"0 has no inverse in GF({self.p})")
        return pow(a, -1, self.p)

    def is_zero(self, a):
        return a % self.p == 0
```

Append to `src/quiverlab/fields/__init__.py`:

```python
from quiverlab.fields.primefield import PrimeField  # noqa: F401
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/fields/primefield.py src/quiverlab/fields/__init__.py tests/fields/test_primefield.py
git commit -m "feat: exact prime fields GF(p)"
```

---

### Task 6: Finite fields GF(p^n) and the GF dispatcher

**Files:**
- Create: `src/quiverlab/fields/finitefield.py`, `src/quiverlab/fields/conway.py`, `tests/fields/test_finitefield.py`
- Modify: `src/quiverlab/fields/__init__.py`, `src/quiverlab/__init__.py`

**Interfaces:**
- Consumes: `PrimeField`, `Domain`, `FieldError`.
- Produces: public **`GF(q, modulus=None)`** — for prime `q` returns `PrimeField`; for `q = p^n` returns `FiniteField(p, n, modulus)` (elements: `tuple[int, ...]` of length `n`, little-endian coefficients; the generator `x` is `(0, 1, 0, ...)`). `FiniteField.gen()` returns the generator. Bundled modulus table `CONWAY[(p, n)]` (validated irreducible at construction; a wrong bundled entry is caught by the tests below and must be replaced from Lübeck's Conway-polynomial tables). `GF` is exported from `quiverlab`.

- [ ] **Step 1: Write the failing test**

`tests/fields/test_finitefield.py`:

```python
import pytest
import sympy
from quiverlab.errors import FieldError
from quiverlab.fields import GF
from quiverlab.fields.conway import CONWAY
from quiverlab.fields.finitefield import FiniteField, poly_is_irreducible


def test_dispatcher():
    assert GF(7).characteristic == 7          # PrimeField
    F4 = GF(4)
    assert isinstance(F4, FiniteField)
    assert F4.characteristic == 2 and F4.q == 4
    with pytest.raises(FieldError):
        GF(6)
    with pytest.raises(FieldError):
        GF(1)


def test_gf4_arithmetic():
    F4 = GF(4)
    x = F4.gen()
    x2 = F4.mul(x, x)
    # modulus x^2 + x + 1: x^2 = x + 1
    assert x2 == F4.add(x, F4.one())
    assert F4.mul(x, F4.inv(x)) == F4.one()
    # Frobenius sanity: a^4 = a for all a in GF(4)
    for a in F4.elements():
        assert F4.mul(F4.mul(a, a), F4.mul(a, a)) == a


def test_every_bundled_conway_entry_is_irreducible_and_primitive():
    for (p, n), coeffs in CONWAY.items():
        assert poly_is_irreducible(coeffs, p), f"CONWAY[{(p, n)}] not irreducible: replace from Lübeck's table"
        F = GF(p**n)
        x = F.gen()
        q = p**n
        # primitivity: order of x is exactly q - 1
        for ell in sympy.factorint(q - 1):
            assert F.pow(x, (q - 1) // ell) != F.one(), \
                f"CONWAY[{(p, n)}] generator not primitive: replace from Lübeck's table"
        assert F.pow(x, q - 1) == F.one()


def test_user_modulus_and_bad_modulus():
    # x^2 + 1 is irreducible over GF(3)
    F9 = GF(9, modulus=[1, 0, 1])
    assert F9.mul(F9.gen(), F9.gen()) == F9.neg(F9.one())
    with pytest.raises(FieldError):
        GF(9, modulus=[2, 0, 1])  # x^2 + 2 = (x+1)(x+2) mod 3: reducible


def test_beyond_table_needs_modulus():
    with pytest.raises(FieldError):
        GF(101**3)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/fields/test_finitefield.py -q`
Expected: FAIL — `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/fields/conway.py`:

```python
"""Bundled modulus polynomials for GF(p^n), little-endian coefficient lists
[c0, ..., cn] with cn = 1 (monic). Sourced from Lübeck's Conway-polynomial
tables; each entry is machine-validated irreducible at field construction and
irreducible+primitive in the test suite — a failing entry must be replaced
from the table, never trusted. Exact Conway normalization (cross-system
embedding compatibility) is a non-goal in this phase."""

CONWAY = {
    (2, 2): [1, 1, 1],
    (2, 3): [1, 1, 0, 1],
    (2, 4): [1, 1, 0, 0, 1],
    (2, 5): [1, 0, 1, 0, 0, 1],
    (2, 6): [1, 1, 0, 1, 1, 0, 1],
    (2, 7): [1, 1, 0, 0, 0, 0, 0, 1],
    (2, 8): [1, 0, 1, 1, 1, 0, 0, 0, 1],
    (2, 9): [1, 0, 0, 0, 1, 0, 0, 0, 0, 1],
    (2, 10): [1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1],
    (3, 2): [2, 2, 1],
    (3, 3): [1, 2, 0, 1],
    (3, 4): [2, 0, 0, 2, 1],
    (3, 5): [1, 2, 0, 0, 0, 1],
    (5, 2): [2, 4, 1],
    (5, 3): [3, 3, 0, 1],
    (5, 4): [2, 4, 4, 0, 1],
    (7, 2): [3, 6, 1],
    (7, 3): [4, 0, 6, 1],
    (11, 2): [2, 7, 1],
    (13, 2): [2, 12, 1],
    (17, 2): [3, 16, 1],
    (19, 2): [2, 18, 1],
    (23, 2): [5, 21, 1],
}
```

`src/quiverlab/fields/finitefield.py`:

```python
"""GF(p^n): exact arithmetic on coefficient tuples modulo a monic irreducible
polynomial over GF(p) (spec §4.3). Elements: tuple[int, ...] of length n,
little-endian."""
import itertools

import sympy

from quiverlab.errors import FieldError
from quiverlab.fields.conway import CONWAY
from quiverlab.fields.domain import Domain, parse_rational
from quiverlab.fields.primefield import PrimeField


def _poly_trim(c):
    while c and c[-1] == 0:
        c = c[:-1]
    return c


def _poly_mod(a, mod, p):
    """a mod (monic) mod, coefficients mod p."""
    a = [x % p for x in a]
    dm = len(mod) - 1
    while len(_poly_trim(a)) - 1 >= dm:
        a = _poly_trim(a)
        d = len(a) - 1
        lead = a[-1]
        for i, c in enumerate(mod):
            a[d - dm + i] = (a[d - dm + i] - lead * c) % p
    return _poly_trim(a)


def _poly_mul(a, b, p):
    out = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
    return _poly_trim(out)


def poly_is_irreducible(coeffs, p):
    """Trial division by every monic polynomial of degree <= deg/2 over GF(p).
    Fine for the bundled table sizes; loudly the wrong tool for huge q."""
    f = _poly_trim([c % p for c in coeffs])
    n = len(f) - 1
    if n < 1 or f[-1] != 1:
        return False
    for d in range(1, n // 2 + 1):
        for tail in itertools.product(range(p), repeat=d):
            g = list(tail) + [1]  # monic degree d
            if not _poly_mod(f, g, p):
                return False
    return True


class FiniteField(Domain):
    def __init__(self, p: int, n: int, modulus=None):
        if modulus is None:
            if (p, n) not in CONWAY:
                raise FieldError(
                    f"GF({p}^{n}) is beyond the bundled modulus table",
                    hint=f"pass an irreducible monic modulus: GF({p}**{n}, modulus=[c0, ..., 1])",
                )
            modulus = CONWAY[(p, n)]
        modulus = [c % p for c in modulus]
        if len(modulus) != n + 1 or modulus[-1] != 1:
            raise FieldError(f"modulus must be monic of degree {n} over GF({p})",
                             hint="little-endian coefficients [c0, ..., 1]")
        if not poly_is_irreducible(modulus, p):
            raise FieldError(f"modulus {modulus} is reducible over GF({p})",
                             hint="pick an irreducible polynomial (see Lübeck's Conway tables)")
        self.p, self.n, self.q = p, n, p**n
        self.modulus = list(modulus)
        self.characteristic = p
        self.name = f"GF({p}^{n})"
        self._prime = PrimeField(p)

    def _tup(self, coeffs):
        c = list(coeffs)[: self.n] + [0] * max(0, self.n - len(coeffs))
        return tuple(x % self.p for x in c)

    def gen(self):
        return self._tup([0, 1])

    def elements(self):
        return (self._tup(t) for t in itertools.product(range(self.p), repeat=self.n))

    def coerce(self, x):
        if isinstance(x, tuple) and len(x) == self.n:
            return self._tup(x)
        return self._tup([self._prime.coerce(x)])

    def zero(self):
        return self._tup([])

    def one(self):
        return self._tup([1])

    def add(self, a, b):
        return tuple((x + y) % self.p for x, y in zip(a, b))

    def neg(self, a):
        return tuple((-x) % self.p for x in a)

    def mul(self, a, b):
        return self._tup(_poly_mod(_poly_mul(list(a), list(b), self.p), self.modulus, self.p))

    def pow(self, a, k: int):
        out = self.one()
        base = a
        while k:
            if k & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            k >>= 1
        return out

    def inv(self, a):
        if self.is_zero(a):
            raise ZeroDivisionError(f"0 has no inverse in {self.name}")
        return self.pow(a, self.q - 2)

    def is_zero(self, a):
        return all(x % self.p == 0 for x in a)

    def to_str(self, a):
        terms = []
        for i, c in enumerate(a):
            if c:
                terms.append(f"{c}" if i == 0 else (f"x^{i}" if c == 1 else f"{c}*x^{i}"))
        return " + ".join(terms) if terms else "0"


def GF(q, modulus=None):
    """The public finite-field constructor: GF(p) or GF(p^n) (spec §3.2)."""
    if not (isinstance(q, int) and not isinstance(q, bool) and q >= 2):
        raise FieldError(f"GF({q!r}): argument must be a prime power >= 2",
                         hint="examples: GF(2), GF(7), GF(4), GF(27)")
    fac = sympy.factorint(q)
    if len(fac) != 1:
        raise FieldError(f"GF({q}): {q} is not a prime power",
                         hint="finite fields exist only for prime powers")
    (p, n), = fac.items()
    if n == 1:
        if modulus is not None:
            raise FieldError(f"GF({p}) takes no modulus", hint="modulus is for GF(p^n), n >= 2")
        return PrimeField(p)
    return FiniteField(p, n, modulus)
```

Append to `src/quiverlab/fields/__init__.py`:

```python
from quiverlab.fields.finitefield import GF, FiniteField  # noqa: F401
```

Append to `src/quiverlab/__init__.py`:

```python
from quiverlab.fields import GF  # noqa: E402,F401
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass. If a `CONWAY` entry fails the irreducibility/primitivity test, replace that entry with the value from Lübeck's Conway-polynomial table (https://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/) — do not weaken the test.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/fields/finitefield.py src/quiverlab/fields/conway.py src/quiverlab/fields/__init__.py src/quiverlab/__init__.py tests/fields/test_finitefield.py
git commit -m "feat: finite fields GF(p^n) with validated bundled moduli and GF dispatcher"
```

---

### Task 7: Exact complex numbers CC

**Files:**
- Create: `src/quiverlab/fields/complexfield.py`, `tests/fields/test_complexfield.py`
- Modify: `src/quiverlab/fields/__init__.py`, `src/quiverlab/__init__.py`

**Interfaces:**
- Consumes: `Domain`, `ExactnessError`, `FieldError`, sympy.
- Produces: singleton **`CC`** (a field *spec*: `CC.parse_entry(x) -> sympy.Expr` exact; `CC.make_domain(entries) -> Domain` wrapping the sympy field generated by the entries — spec §4.2); helper **`E(n)`** = primitive n-th root of unity; internal `SympyExactDomain`. `CC` and `E` exported from `quiverlab`. Entries accepted: int, `Fraction`, exact strings (`"1/3"`, `"i"`, `"sqrt(2)"`, `"E(3)"`), sympy exact expressions. Floats, decimal strings, and sympy `Float` atoms raise `ExactnessError`. Entries sympy cannot place in an algebraic extension raise `FieldError` (loud, spec §4.2).

- [ ] **Step 1: Write the failing test**

`tests/fields/test_complexfield.py`:

```python
from fractions import Fraction

import pytest
import sympy
from quiverlab.errors import ExactnessError, FieldError
from quiverlab.fields import CC, E


def test_rational_entries_fast_path():
    dom = CC.make_domain([1, Fraction(1, 2), "2/3"])
    a = dom.coerce("2/3")
    assert dom.characteristic == 0
    assert dom.eq(dom.add(a, a), dom.coerce("4/3"))
    assert dom.is_zero(dom.sub(a, a))


def test_i_and_sqrt2():
    dom = CC.make_domain(["i", "sqrt(2)"])
    i = dom.coerce("i")
    s = dom.coerce("sqrt(2)")
    assert dom.eq(dom.mul(i, i), dom.coerce(-1))
    assert dom.eq(dom.mul(s, s), dom.coerce(2))
    assert dom.eq(dom.mul(i, dom.inv(i)), dom.one())


def test_root_of_unity():
    dom = CC.make_domain([E(3)])
    w = dom.coerce(E(3))
    w3 = dom.mul(dom.mul(w, w), w)
    assert dom.eq(w3, dom.one())
    # 1 + w + w^2 = 0
    assert dom.is_zero(dom.add(dom.add(dom.one(), w), dom.mul(w, w)))


def test_floats_fail_loudly_everywhere():
    with pytest.raises(ExactnessError):
        CC.parse_entry(0.5)
    with pytest.raises(ExactnessError):
        CC.parse_entry("0.5")
    with pytest.raises(ExactnessError):
        CC.parse_entry(sympy.Float("0.5"))
    with pytest.raises(ExactnessError):
        CC.parse_entry(1 + 2j)


def test_non_number_rejected():
    with pytest.raises(FieldError):
        CC.parse_entry("x + 1")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/fields/test_complexfield.py -q`
Expected: FAIL — `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/fields/complexfield.py`:

```python
"""CC: 'the complex numbers', exactly (spec §4.2). Entries must be exact
algebraic expressions; the actual computation domain is the exact subfield of C
they generate (sympy algebraic field). Dimensions computed there provably equal
the C answers by flat base change."""
from fractions import Fraction

import sympy
from sympy.polys.constructor import construct_domain

from quiverlab.errors import ExactnessError, FieldError
from quiverlab.fields.domain import Domain, reject_inexact

_FLOAT_HINT = "quiverlab is exact-only: write '1/2', 'sqrt(2)', 'i', E(3), never 0.5"


def E(n: int):
    """Primitive n-th root of unity exp(2*pi*i/n), exact (GAP's E(n) convention)."""
    if not (isinstance(n, int) and not isinstance(n, bool) and n >= 1):
        raise FieldError(f"E({n!r}): n must be a positive integer", hint="e.g. E(3), E(8)")
    return sympy.exp(2 * sympy.pi * sympy.I / n)


class SympyExactDomain(Domain):
    """Wraps a sympy field domain (QQ or QQ<alpha>) behind the quiverlab Domain protocol."""

    characteristic = 0

    def __init__(self, sdom):
        self.sdom = sdom
        self.name = f"CC (computing exactly in {sdom})"

    def coerce(self, x):
        expr = CC.parse_entry(x)
        try:
            return self.sdom.from_sympy(expr)
        except Exception as exc:
            raise FieldError(
                f"{expr} does not lie in the exact field {self.sdom} generated by this "
                f"algebra's entries",
                hint="include this number among the algebra's entries at construction",
            ) from exc

    def zero(self):
        return self.sdom.zero

    def one(self):
        return self.sdom.one

    def add(self, a, b):
        return a + b

    def neg(self, a):
        return -a

    def mul(self, a, b):
        return a * b

    def inv(self, a):
        if self.is_zero(a):
            raise ZeroDivisionError("0 has no inverse")
        return self.sdom.exquo(self.sdom.one, a)

    def is_zero(self, a):
        return a == self.sdom.zero

    def to_str(self, a):
        return str(self.sdom.to_sympy(a))


class ComplexField:
    """Field spec for 'the complex numbers'. Not itself a Domain: it inspects all
    entries first, then hands back the exact working Domain (make_domain)."""

    name = "CC"
    characteristic = 0

    def parse_entry(self, x):
        if isinstance(x, sympy.Expr):
            expr = x
        else:
            reject_inexact(x)
            if isinstance(x, int):
                expr = sympy.Integer(x)
            elif isinstance(x, Fraction):
                expr = sympy.Rational(x.numerator, x.denominator)
            elif isinstance(x, str):
                try:
                    expr = sympy.sympify(
                        x, locals={"i": sympy.I, "E": E}, rational=False, evaluate=True
                    )
                except (sympy.SympifyError, TypeError, SyntaxError) as exc:
                    raise FieldError(f"cannot read {x!r} as an exact complex number",
                                     hint="examples: '1/3', 'i', 'sqrt(2)', 'E(3)'") from exc
            else:
                raise FieldError(f"cannot read {type(x).__name__} {x!r} as an exact scalar",
                                 hint=_FLOAT_HINT)
        if expr.atoms(sympy.Float):
            raise ExactnessError(f"{expr} contains a floating-point number", hint=_FLOAT_HINT)
        if not expr.is_number:
            raise FieldError(f"{expr} is not a number", hint="entries must be exact scalars")
        return expr

    def make_domain(self, entries):
        exprs = [self.parse_entry(e) for e in entries]
        try:
            sdom, _ = construct_domain(exprs or [sympy.Integer(0)], field=True, extension=True)
        except Exception as exc:
            raise FieldError(
                "sympy cannot place these entries in an exact algebraic extension of QQ "
                "(spec §4.2: transcendental or unsupported entries fail loudly)",
                hint="stick to rationals, i, radicals and roots of unity E(n)",
            ) from exc
        return SympyExactDomain(sdom)


CC = ComplexField()
```

Append to `src/quiverlab/fields/__init__.py`:

```python
from quiverlab.fields.complexfield import CC, E, SympyExactDomain  # noqa: F401
```

Append to `src/quiverlab/__init__.py`:

```python
from quiverlab.fields import CC, E  # noqa: E402,F401
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass. Note: `test_root_of_unity` exercises sympy's minimal-polynomial machinery on `exp(2*pi*I/3)`; if a sympy version cannot handle an `E(n)`, the code path already raises the loud `FieldError` — the test uses `E(3)`, which current sympy handles.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/fields/complexfield.py src/quiverlab/fields/__init__.py src/quiverlab/__init__.py tests/fields/test_complexfield.py
git commit -m "feat: exact complex field CC via sympy algebraic extensions"
```

---

### Task 8: Quiver and paths

**Files:**
- Create: `src/quiverlab/combinat/__init__.py`, `src/quiverlab/combinat/quiver.py`, `tests/combinat/__init__.py` (empty), `tests/combinat/test_quiver.py`

**Interfaces:**
- Consumes: `RelationError`.
- Produces: **`Quiver(vertices, arrows)`** with `.vertices: list`, `.arrows: dict[str, tuple]`, `.source(name)`, `.target(name)`, `.is_acyclic() -> bool`, `.word_source(word) / .word_target(word)` for arrow-name tuples, `.compose_ok(word) -> bool`, `__repr__`. Arrow names must match `[A-Za-z_][A-Za-z0-9_]*`. `Quiver.algebra(...)` is attached in Task 11. Exported from `quiverlab`.

- [ ] **Step 1: Write the failing test**

`tests/combinat/test_quiver.py`:

```python
import pytest
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import RelationError


def _q():
    return Quiver(vertices=[1, 2, 3],
                  arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})


def test_basic_accessors():
    Q = _q()
    assert Q.source("a") == 1 and Q.target("a") == 2
    assert Q.word_source(("a", "b")) == 1 and Q.word_target(("a", "b")) == 3
    assert Q.compose_ok(("a", "b")) and not Q.compose_ok(("b", "a"))
    assert Q.is_acyclic()


def test_loops_and_cycles():
    L = Quiver(vertices=[1], arrows={"x": (1, 1)})
    assert not L.is_acyclic()
    assert L.compose_ok(("x", "x", "x"))


def test_validation():
    with pytest.raises(RelationError):
        Quiver(vertices=[1], arrows={"a": (1, 2)})       # endpoint not a vertex
    with pytest.raises(RelationError):
        Quiver(vertices=[1], arrows={"a*b": (1, 1)})     # bad arrow name
    with pytest.raises(RelationError):
        Quiver(vertices=[1, 1], arrows={})               # duplicate vertex


def test_repr_shows_arrows():
    s = repr(_q())
    assert "1 --a--> 2" in s
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/combinat/test_quiver.py -q`
Expected: FAIL — `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/combinat/quiver.py`:

```python
"""Quiver = finite directed multigraph with named arrows. Paths are tuples of
arrow names read LEFT TO RIGHT: ('a', 'b') means first a, then b, and requires
target(a) == source(b) (Assem-Simson-Skowronski convention)."""
import re

from quiverlab.errors import RelationError

_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class Quiver:
    def __init__(self, vertices, arrows):
        vertices = list(vertices)
        if len(set(vertices)) != len(vertices):
            raise RelationError("duplicate vertices", hint="each vertex must appear once")
        self.vertices = vertices
        vset = set(vertices)
        self.arrows = {}
        for name, ends in dict(arrows).items():
            if not (isinstance(name, str) and _NAME.match(name)):
                raise RelationError(
                    f"bad arrow name {name!r}",
                    hint="arrow names must be identifiers like a, b2, alpha (they appear in relation strings)",
                )
            try:
                s, t = ends
            except (TypeError, ValueError):
                raise RelationError(f"arrow {name!r}: endpoints must be a (source, target) pair",
                                    hint="e.g. arrows={'a': (1, 2)}") from None
            if s not in vset or t not in vset:
                raise RelationError(f"arrow {name!r}: endpoint not a vertex",
                                    hint=f"vertices are {vertices}")
            self.arrows[name] = (s, t)

    # -- accessors -----------------------------------------------------------
    def source(self, name):
        return self.arrows[name][0]

    def target(self, name):
        return self.arrows[name][1]

    def word_source(self, word):
        return self.arrows[word[0]][0] if word else None

    def word_target(self, word):
        return self.arrows[word[-1]][1] if word else None

    def compose_ok(self, word) -> bool:
        return all(self.target(a) == self.source(b) for a, b in zip(word, word[1:]))

    def is_acyclic(self) -> bool:
        adj = {v: [] for v in self.vertices}
        for s, t in self.arrows.values():
            adj[s].append(t)
        WHITE, GRAY, BLACK = 0, 1, 2
        color = {v: WHITE for v in self.vertices}
        for start in self.vertices:
            if color[start] != WHITE:
                continue
            stack = [(start, iter(adj[start]))]
            color[start] = GRAY
            while stack:
                v, it = stack[-1]
                nxt = next(it, None)
                if nxt is None:
                    color[v] = BLACK
                    stack.pop()
                elif color[nxt] == GRAY:
                    return False
                elif color[nxt] == WHITE:
                    color[nxt] = GRAY
                    stack.append((nxt, iter(adj[nxt])))
        return True

    def __repr__(self):
        lines = [f"Quiver with vertices {self.vertices} and arrows:"]
        lines += [f"  {s} --{a}--> {t}" for a, (s, t) in self.arrows.items()]
        return "\n".join(lines)
```

`src/quiverlab/combinat/__init__.py`:

```python
from quiverlab.combinat.quiver import Quiver  # noqa: F401
```

Append to `src/quiverlab/__init__.py`:

```python
from quiverlab.combinat import Quiver  # noqa: E402,F401
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/combinat tests/combinat src/quiverlab/__init__.py
git commit -m "feat: Quiver with left-to-right path words"
```

---

### Task 9: Relation parser

**Files:**
- Create: `src/quiverlab/combinat/relations.py`, `tests/combinat/test_relations.py`
- Modify: `src/quiverlab/combinat/__init__.py`

**Interfaces:**
- Consumes: `Quiver`, `RelationError`, `parse_rational`.
- Produces: **`Relation`** dataclass: `.terms: list[tuple[Fraction, tuple[str, ...]]]` (coefficient, arrow word), `.source`, `.target`, `.is_monomial: bool` (single term), `.max_length`, `.min_length`; **`parse_relation(s: str, quiver) -> Relation`** and **`parse_relations(strings, quiver) -> list[Relation]`**. Grammar this phase: terms joined by `+`/`-`; each term = optional rational coefficient (`2*`, `-1/3*`) followed by `*`-joined arrow names; `name^k` = k-fold repetition. Coefficients beyond rationals raise `RelationError` (general exact coefficients arrive with the Gröbner phase, Plan 03).

- [ ] **Step 1: Write the failing test**

`tests/combinat/test_relations.py`:

```python
from fractions import Fraction as F

import pytest
from quiverlab.combinat import Quiver
from quiverlab.combinat.relations import parse_relation
from quiverlab.errors import RelationError


def _q():
    return Quiver(vertices=[1, 2, 3],
                  arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3), "x": (1, 1)})


def test_monomial():
    r = parse_relation("a*b", _q())
    assert r.is_monomial and r.terms == [(F(1), ("a", "b"))]
    assert (r.source, r.target) == (1, 3)


def test_power():
    r = parse_relation("x^3", _q())
    assert r.terms == [(F(1), ("x", "x", "x"))]


def test_linear_combination_parallel():
    r = parse_relation("a*b - 2*c", _q())
    assert not r.is_monomial
    assert r.terms == [(F(1), ("a", "b")), (F(-2), ("c",))]
    assert (r.source, r.target) == (1, 3)


def test_fraction_coefficient():
    r = parse_relation("1/2*a*b + c", _q())
    assert r.terms[0][0] == F(1, 2)


def test_not_composable():
    with pytest.raises(RelationError) as ei:
        parse_relation("b*a", _q())
    assert "target" in str(ei.value)


def test_not_parallel():
    with pytest.raises(RelationError) as ei:
        parse_relation("a - c", _q())
    assert "parallel" in str(ei.value)


def test_unknown_arrow():
    with pytest.raises(RelationError):
        parse_relation("a*z", _q())


def test_decimal_coefficient_fails():
    with pytest.raises(Exception):  # ExactnessError via parse_rational
        parse_relation("0.5*a*b", _q())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/combinat/test_relations.py -q`
Expected: FAIL — `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/combinat/relations.py`:

```python
"""Relation strings -> exact linear combinations of parallel paths (spec §3.3).
Grammar (this phase): terms joined by + and -; term = [rational*] arrows;
'p^k' repeats an arrow k times. Paths read LEFT TO RIGHT."""
import re
from dataclasses import dataclass
from fractions import Fraction

from quiverlab.errors import RelationError
from quiverlab.fields.domain import parse_rational

_COEFF = re.compile(r"^[+-]?\d+(/\d+)?$")
_POW = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\^(\d+)$")


@dataclass(frozen=True)
class Relation:
    terms: tuple  # tuple[tuple[Fraction, tuple[str, ...]], ...]
    source: object
    target: object

    @property
    def is_monomial(self) -> bool:
        return len(self.terms) == 1

    @property
    def max_length(self) -> int:
        return max(len(w) for _, w in self.terms)

    @property
    def min_length(self) -> int:
        return min(len(w) for _, w in self.terms)

    def __repr__(self):
        bits = []
        for c, w in self.terms:
            path = "*".join(w)
            if c == 1:
                bits.append(path)
            elif c == -1:
                bits.append(f"-{path}")
            else:
                bits.append(f"{c}*{path}")
        return " + ".join(bits).replace("+ -", "- ")


def _split_terms(s: str):
    """'a*b - 2*c + d' -> ['a*b', '-2*c', '+d'] (sign attached)."""
    s = s.replace(" ", "")
    if not s:
        raise RelationError("empty relation string", hint="e.g. 'a*b - c' or 'x^2'")
    out, cur = [], ""
    for ch in s:
        if ch in "+-" and cur:
            out.append(cur)
            cur = ch
        else:
            cur += ch
    out.append(cur)
    return out


def _parse_term(term: str, quiver):
    sign = Fraction(1)
    if term and term[0] in "+-":
        sign = Fraction(-1) if term[0] == "-" else Fraction(1)
        term = term[1:]
    factors = term.split("*")
    coeff = sign
    word: list[str] = []
    for f in factors:
        if not f:
            raise RelationError(f"malformed term {term!r}", hint="factors are joined by single *")
        if _COEFF.match(f):
            if word:
                raise RelationError(f"coefficient {f!r} appears after arrows in {term!r}",
                                    hint="write coefficients first: '2*a*b'")
            coeff = coeff * parse_rational(f)
            continue
        m = _POW.match(f)
        if m:
            name, k = m.group(1), int(m.group(2))
            reps = [name] * k
        else:
            name, reps = f, [f]
        if name not in quiver.arrows:
            raise RelationError(f"unknown arrow {name!r} in relation term {term!r}",
                                hint=f"arrows are {sorted(quiver.arrows)}")
        word.extend(reps)
    if not word:
        raise RelationError(f"term {term!r} has no arrows",
                            hint="relations live in the arrow ideal")
    w = tuple(word)
    if not quiver.compose_ok(w):
        bad = next((i for i, (a, b) in enumerate(zip(w, w[1:]))
                    if quiver.target(a) != quiver.source(b)))
        raise RelationError(
            f"path {'*'.join(w)} is not composable: target({w[bad]}) = "
            f"{quiver.target(w[bad])} but source({w[bad + 1]}) = {quiver.source(w[bad + 1])}",
            hint="paths compose left to right: a*b needs target(a) == source(b)",
        )
    return coeff, w


def parse_relation(s: str, quiver) -> Relation:
    if not isinstance(s, str):
        raise RelationError(f"relations are strings, got {type(s).__name__}",
                            hint="e.g. 'a*b - c'")
    terms = [_parse_term(t, quiver) for t in _split_terms(s)]
    terms = [(c, w) for c, w in terms if c != 0]
    if not terms:
        raise RelationError(f"relation {s!r} is identically zero", hint="remove it")
    srcs = {quiver.word_source(w) for _, w in terms}
    tgts = {quiver.word_target(w) for _, w in terms}
    if len(srcs) > 1 or len(tgts) > 1:
        raise RelationError(
            f"terms of {s!r} are not parallel (sources {sorted(srcs, key=str)}, "
            f"targets {sorted(tgts, key=str)})",
            hint="every summand of a relation must share source and target",
        )
    return Relation(tuple(terms), srcs.pop(), tgts.pop())


def parse_relations(strings, quiver) -> list:
    return [parse_relation(s, quiver) for s in strings]
```

Append to `src/quiverlab/combinat/__init__.py`:

```python
from quiverlab.combinat.relations import Relation, parse_relation, parse_relations  # noqa: F401
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/combinat tests/combinat/test_relations.py
git commit -m "feat: relation parser (exact coefficients, parallel-path validation)"
```

---

### Task 10: The Algebra core

**Files:**
- Create: `src/quiverlab/core/__init__.py`, `src/quiverlab/core/algebra.py`, `tests/core/__init__.py` (empty), `tests/core/test_algebra.py`

**Interfaces:**
- Consumes: `Domain`, `linalg.solve`, errors.
- Produces: **`Algebra`** with attributes `.domain`, `.dim: int`, `.T: list[list[list]]` (`T[i][j]` = vector of `b_i * b_j` in the basis), `.unit: list`, `.basis_labels: list[str] | None`, `.is_unit_adapted: bool`; methods `.multiply(u, v) -> list`, `.unit_adapted() -> Algebra` (basis 0 becomes 1_A), `.change_of_basis(P) -> Algebra` (columns of P = new basis in old coordinates); classmethod **`Algebra.from_structure_constants(T, unit, field=CC, check=True, basis_labels=None)`** (coerces entries via `field.make_domain`, validates unit and associativity when `check=True`). This exact shape is the frozen internal currency for Plans 02–08.

- [ ] **Step 1: Write the failing test**

`tests/core/test_algebra.py`:

```python
from fractions import Fraction as F

import pytest
from quiverlab.core.algebra import Algebra
from quiverlab.errors import ExactnessError, FieldError, QuiverlabError
from quiverlab.fields import CC, GF


def _dual_numbers(field=CC):
    # basis (1, x), x^2 = 0
    T = [
        [[1, 0], [0, 1]],
        [[0, 1], [0, 0]],
    ]
    return Algebra.from_structure_constants(T, unit=[1, 0], field=field)


def test_dual_numbers_multiply():
    A = _dual_numbers()
    dom = A.domain
    one = [dom.coerce(1), dom.coerce(0)]
    x = [dom.coerce(0), dom.coerce(1)]
    assert A.multiply(x, x) == [dom.zero(), dom.zero()]
    assert A.multiply(one, x) == x
    assert A.is_unit_adapted


def test_associativity_check_catches_garbage():
    T = [
        [[1, 0], [0, 1]],
        [[0, 1], [1, 0]],  # x*x = 1 ... with x*1 = x this IS associative (k[x]/(x^2-1)); tweak:
    ]
    # break unit instead: claim unit = [0, 1]
    with pytest.raises(QuiverlabError):
        Algebra.from_structure_constants(T, unit=[0, 1], field=CC)


def test_nonassociative_rejected():
    # b1*b1 = b1 with b1 the unit, b2*b2 = b1, b2*b1 = b2, but b1*b2 = 0: not unital/associative
    T = [
        [[1, 0], [0, 0]],
        [[0, 1], [1, 0]],
    ]
    with pytest.raises(QuiverlabError):
        Algebra.from_structure_constants(T, unit=[1, 0], field=CC)


def test_floats_rejected_in_T():
    T = [
        [[1, 0], [0, 1]],
        [[0, 1], [0.0, 0]],
    ]
    with pytest.raises(ExactnessError):
        Algebra.from_structure_constants(T, unit=[1, 0], field=CC)


def test_gf_algebra():
    A = _dual_numbers(field=GF(2))
    assert A.domain.characteristic == 2
    assert A.dim == 2


def test_unit_adapted_transform():
    # k x k with basis (e1, e2): unit = e1 + e2 is NOT a basis vector
    T = [
        [[1, 0], [0, 0]],
        [[0, 0], [0, 1]],
    ]
    A = Algebra.from_structure_constants(T, unit=[1, 1], field=CC)
    assert not A.is_unit_adapted
    B = A.unit_adapted()
    assert B.is_unit_adapted
    dom = B.domain
    # in the new basis, b0 is the unit
    e0 = [dom.one()] + [dom.zero()] * (B.dim - 1)
    v = [dom.coerce(3), dom.coerce(-2)]
    assert B.multiply(e0, v) == v and B.multiply(v, e0) == v
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/core/test_algebra.py -q`
Expected: FAIL — `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/core/algebra.py`:

```python
"""The structure-constant Algebra: quiverlab's internal currency (spec §5).
T[i][j] is the coordinate vector of b_i * b_j. 'Unit-adapted' means b_0 = 1_A
(hanlab's convention), which the bar complex requires."""
from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import solve


class Algebra:
    def __init__(self, domain, T, unit, basis_labels=None, is_unit_adapted=None, _quiver=None,
                 _relations=None):
        self.domain = domain
        self.T = T
        self.unit = unit
        self.dim = len(T)
        self.basis_labels = basis_labels
        self.quiver = _quiver
        self.relations = _relations
        if is_unit_adapted is None:
            one = domain.one()
            is_unit_adapted = (
                not domain.is_zero(unit[0])
                and domain.eq(unit[0], one)
                and all(domain.is_zero(c) for c in unit[1:])
            )
        self.is_unit_adapted = is_unit_adapted

    # -- arithmetic -----------------------------------------------------------
    def multiply(self, u, v):
        dom = self.domain
        out = [dom.zero()] * self.dim
        for i, ui in enumerate(u):
            if dom.is_zero(ui):
                continue
            for j, vj in enumerate(v):
                if dom.is_zero(vj):
                    continue
                c = dom.mul(ui, vj)
                for t, w in enumerate(self.T[i][j]):
                    if not dom.is_zero(w):
                        out[t] = dom.add(out[t], dom.mul(c, w))
        return out

    # -- construction ---------------------------------------------------------
    @classmethod
    def from_structure_constants(cls, T, unit, field=None, check=True, basis_labels=None):
        if field is None:
            from quiverlab.fields import CC
            field = CC
        m = len(T)
        raw = [x for row in T for vec in row for x in vec] + list(unit)
        parsed = [field.parse_entry(x) for x in raw]
        dom = field.make_domain(parsed)
        Tc = [[[dom.coerce(field.parse_entry(x)) for x in T[i][j]] for j in range(m)]
              for i in range(m)]
        unit_c = [dom.coerce(field.parse_entry(x)) for x in unit]
        A = cls(dom, Tc, unit_c, basis_labels=basis_labels)
        if check:
            A._validate()
        return A

    def _basis_vec(self, i):
        dom = self.domain
        v = [dom.zero()] * self.dim
        v[i] = dom.one()
        return v

    def _validate(self):
        dom = self.domain
        for i in range(self.dim):
            bi = self._basis_vec(i)
            left = self.multiply(self.unit, bi)
            right = self.multiply(bi, self.unit)
            if left != bi or right != bi:
                raise QuiverlabError(
                    f"the given unit vector is not a two-sided unit (fails on basis {i})",
                    hint="check the structure constants and the unit coordinates",
                )
        for i in range(self.dim):
            for j in range(self.dim):
                ij = self.T[i][j]
                for k in range(self.dim):
                    lhs = self.multiply(ij, self._basis_vec(k))
                    rhs = self.multiply(self._basis_vec(i), self.T[j][k])
                    if lhs != rhs:
                        raise QuiverlabError(
                            f"structure constants are not associative: (b{i}·b{j})·b{k} != b{i}·(b{j}·b{k})",
                            hint="re-derive the multiplication table; quiverlab never guesses",
                        )

    # -- base change ----------------------------------------------------------
    def change_of_basis(self, P):
        """New algebra in the basis whose j-th vector has old coordinates column j of P."""
        dom = self.domain
        m = self.dim
        cols = [[P[i][j] for i in range(m)] for j in range(m)]
        newT = []
        for i in range(m):
            row = []
            for j in range(m):
                prod_old = self.multiply(cols[i], cols[j])
                x = solve(P, prod_old, dom)
                if x is None:
                    raise QuiverlabError("change of basis matrix is singular",
                                         hint="columns must form a basis")
                row.append(x)
            newT.append(row)
        new_unit = solve(P, list(self.unit), dom)
        if new_unit is None:
            raise QuiverlabError("change of basis matrix is singular",
                                 hint="columns must form a basis")
        return Algebra(dom, newT, new_unit, basis_labels=None,
                       _quiver=self.quiver, _relations=self.relations)

    def unit_adapted(self):
        """Return an isomorphic copy whose basis vector 0 is 1_A (spec §5, component 4)."""
        if self.is_unit_adapted:
            return self
        dom = self.domain
        m = self.dim
        j = next(i for i, c in enumerate(self.unit) if not dom.is_zero(c))
        P = [[dom.one() if r == c else dom.zero() for c in range(m)] for r in range(m)]
        for r in range(m):
            P[r][j] = self.unit[r]
        if j != 0:
            for r in range(m):
                P[r][0], P[r][j] = P[r][j], P[r][0]
        out = self.change_of_basis(P)
        labels = None
        if self.basis_labels is not None:
            labels = list(self.basis_labels)
            old0 = labels[j]
            labels[j] = old0 if j == 0 else labels[0]
            labels[0] = "1"
            if j == 0:
                labels[0] = "1"
        out.basis_labels = labels
        out.is_unit_adapted = True
        return out

    def __repr__(self):
        base = f"Algebra of dimension {self.dim} over {self.domain.name}"
        if self.basis_labels:
            base += "\nbasis: " + ", ".join(self.basis_labels)
        return base
```

`src/quiverlab/core/__init__.py`:

```python
from quiverlab.core.algebra import Algebra  # noqa: F401
```

Append to `src/quiverlab/__init__.py`:

```python
from quiverlab.core import Algebra  # noqa: E402,F401
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/core tests/core src/quiverlab/__init__.py
git commit -m "feat: structure-constant Algebra with validation and unit adaptation"
```

---

### Task 11: Monomial kQ/I with certified finiteness

**Files:**
- Create: `src/quiverlab/core/monomial.py`, `tests/core/test_monomial.py`
- Modify: `src/quiverlab/combinat/quiver.py` (add `Quiver.algebra`), `src/quiverlab/core/__init__.py`

**Interfaces:**
- Consumes: `Quiver`, `Relation`, `parse_relations`, `Algebra`, errors.
- Produces: **`Quiver.algebra(relations=(), field=CC) -> Algebra`** — parses relations; monomial presentations are lowered to `Algebra` (basis = trivial paths + irreducible paths, labels like `"e_1"`, `"a*b"`); non-monomial relations raise `NotImplementedError` naming Plan 03; relations of length < 2 raise `AdmissibilityError`; infinite irreducible growth raises `NotFiniteDimensionalError` **showing an arrow cycle**. Also `irreducible_paths(quiver, forbidden) -> list[tuple[str, ...]]` (the certifying automaton, reused by Bardzell in Plan 02).

- [ ] **Step 1: Write the failing test**

`tests/core/test_monomial.py`:

```python
import pytest
from quiverlab import CC, GF, Quiver
from quiverlab.errors import AdmissibilityError, NotFiniteDimensionalError


def test_dual_numbers_from_quiver():
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x^2"], field=GF(2))
    assert A.dim == 2                      # e_1, x
    assert "x" in (A.basis_labels or [])


def test_truncated_loop():
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x^4"])
    assert A.dim == 4                      # e, x, x^2, x^3


def test_a3_rad_square_zero():
    Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
    A = Q.algebra(relations=["a*b"])
    assert A.dim == 5                      # e1, e2, e3, a, b


def test_hereditary_no_relations():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    A = Q.algebra()
    assert A.dim == 3                      # e1, e2, a


def test_loop_without_relations_is_infinite():
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    with pytest.raises(NotFiniteDimensionalError) as ei:
        Q.algebra()
    assert "x" in str(ei.value)            # the offending cycle is named


def test_two_loops_one_relation_still_infinite():
    # k<x, y>/(xy): the words y^j x^i are all irreducible -> infinite-dimensional,
    # and the automaton must name a cycle (a loop on x or on y).
    Q = Quiver(vertices=[1], arrows={"x": (1, 1), "y": (1, 1)})
    with pytest.raises(NotFiniteDimensionalError) as ei:
        Q.algebra(relations=["x*y"])
    msg = str(ei.value)
    assert ("x" in msg) or ("y" in msg)


def test_short_relation_not_admissible():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    with pytest.raises(AdmissibilityError):
        Q.algebra(relations=["a"])


def test_nonmonomial_waits_for_plan03():
    Q = Quiver(vertices=[1, 2, 3],
               arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
    with pytest.raises(NotImplementedError):
        Q.algebra(relations=["a*b - c"])


def test_multiplication_table_is_path_concatenation():
    Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
    A = Q.algebra(relations=[], field=CC)  # kA_3: e1,e2,e3,a,b,a*b -> dim 6
    assert A.dim == 6
    labels = A.basis_labels
    ia, ib, iab = labels.index("a"), labels.index("b"), labels.index("a*b")
    dom = A.domain
    va = [dom.one() if i == ia else dom.zero() for i in range(6)]
    vb = [dom.one() if i == ib else dom.zero() for i in range(6)]
    assert A.multiply(va, vb)[iab] == dom.one()
    assert all(dom.is_zero(c) for i, c in enumerate(A.multiply(vb, va)))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/core/test_monomial.py -q`
Expected: FAIL — `Quiver` has no `algebra` / `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/core/monomial.py`:

```python
"""Monomial kQ/I -> Algebra with a certified finite basis (spec §3.3, §5 component 4).
Basis = trivial paths e_v + irreducible paths (words avoiding every forbidden
word as a contiguous subword). Finiteness is decided by a suffix-window
automaton; an infinite family is reported with an explicit arrow cycle."""
from collections import deque

from quiverlab.core.algebra import Algebra
from quiverlab.errors import AdmissibilityError, NotFiniteDimensionalError


def _contains_forbidden(word, forbidden):
    return any(
        word[i: i + len(f)] == f
        for f in forbidden
        for i in range(len(word) - len(f) + 1)
    )


def _automaton(quiver, forbidden):
    """States (vertex, window) reachable by irreducible words; window = last r-1 arrows."""
    r = max((len(f) for f in forbidden), default=1)
    starts = [(v, ()) for v in quiver.vertices]
    graph = {}
    seen = set(starts)
    dq = deque(starts)
    while dq:
        state = dq.popleft()
        v, w = state
        outs = []
        for a, (s, t) in quiver.arrows.items():
            if s != v:
                continue
            word = w + (a,)
            if any(len(f) <= len(word) and word[-len(f):] == f for f in forbidden):
                continue
            ns = (t, word[-(r - 1):] if r > 1 else ())
            outs.append((a, ns))
            if ns not in seen:
                seen.add(ns)
                dq.append(ns)
        graph[state] = outs
    return starts, graph


def _find_cycle(graph):
    """Return the arrow labels of a cycle in the state graph, or None."""
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {s: WHITE for s in graph}
    for root in graph:
        if color[root] != WHITE:
            continue
        stack = [(root, iter(graph[root]), None)]
        path = []  # arrows along the current DFS path
        color[root] = GRAY
        while stack:
            state, it, _ = stack[-1]
            step = next(it, None)
            if step is None:
                color[state] = BLACK
                stack.pop()
                if path:
                    path.pop()
                continue
            arrow, ns = step
            if color[ns] == GRAY:
                # unwind to where ns sits on the stack
                idx = next(i for i, fr in enumerate(stack) if fr[0] == ns)
                return path[idx:] + [arrow]
            if color[ns] == WHITE:
                color[ns] = GRAY
                stack.append((ns, iter(graph[ns]), arrow))
                path.append(arrow)
    return None


def irreducible_paths(quiver, forbidden):
    """All irreducible words, sorted by (length, word). Raises
    NotFiniteDimensionalError (with a cycle) when infinitely many exist."""
    forbidden = [tuple(f) for f in forbidden]
    starts, graph = _automaton(quiver, forbidden)
    cycle = _find_cycle(graph)
    if cycle is not None:
        raise NotFiniteDimensionalError(
            "kQ/I is infinite-dimensional: irreducible paths grow forever along the "
            "cycle " + " -> ".join(cycle),
            hint="add relations killing a power of this cycle (monomial), or check your quiver",
        )
    words = []
    for st in starts:
        stack = [(st, ())]
        while stack:
            state, word = stack.pop()
            for arrow, ns in graph[state]:
                nw = word + (arrow,)
                words.append(nw)
                stack.append((ns, nw))
    return sorted(set(words), key=lambda w: (len(w), w))


def build_monomial_algebra(quiver, relations, field):
    """relations: parsed Relation objects, all monomial, lengths >= 2."""
    for rel in relations:
        if rel.min_length < 2:
            raise AdmissibilityError(
                f"relation {rel!r} has a path of length {rel.min_length}: the ideal is "
                "not inside the square of the arrow ideal",
                hint="admissible relations use paths of length >= 2",
            )
    forbidden = [w for rel in relations for _, w in rel.terms]  # monomial: one word each
    words = irreducible_paths(quiver, forbidden)
    trivial = [("e", v) for v in quiver.vertices]
    basis = trivial + [("p", w) for w in words]
    index = {b: i for i, b in enumerate(basis)}
    m = len(basis)

    def src(b):
        return b[1] if b[0] == "e" else quiver.word_source(b[1])

    def tgt(b):
        return b[1] if b[0] == "e" else quiver.word_target(b[1])

    def prod(x, y):
        if tgt(x) != src(y):
            return None
        if x[0] == "e":
            return y
        if y[0] == "e":
            return x
        w = x[1] + y[1]
        if _contains_forbidden(w, [tuple(f) for f in forbidden]):
            return None
        return ("p", w)

    parsed_pool = [field.parse_entry(0), field.parse_entry(1)]
    dom = field.make_domain(parsed_pool)
    zero, one = dom.zero(), dom.one()
    T = [[[zero] * m for _ in range(m)] for _ in range(m)]
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            p = prod(bi, bj)
            vec = [zero] * m
            if p is not None:
                vec[index[p]] = one
            T[i][j] = vec
    unit = [zero] * m
    for v in quiver.vertices:
        unit[index[("e", v)]] = one
    labels = [f"e_{b[1]}" if b[0] == "e" else "*".join(b[1]) for b in basis]
    return Algebra(dom, T, unit, basis_labels=labels, _quiver=quiver,
                   _relations=list(relations))
```

Add to `src/quiverlab/combinat/quiver.py` (method on `Quiver`, at the end of the class):

```python
    def algebra(self, relations=(), field=None):
        """Build kQ/I over the field (default CC). Monomial presentations are
        certified and lowered; general relations arrive with the Groebner engine."""
        from quiverlab.combinat.relations import parse_relations
        from quiverlab.core.monomial import build_monomial_algebra

        if field is None:
            from quiverlab.fields import CC
            field = CC
        rels = parse_relations(list(relations), self)
        if all(r.is_monomial for r in rels):
            return build_monomial_algebra(self, rels, field)
        raise NotImplementedError(
            "general (non-monomial) relations arrive with the Groebner engine "
            "(Plan 03); this build certifies monomial presentations only"
        )
```

Append to `src/quiverlab/core/__init__.py`:

```python
from quiverlab.core.monomial import build_monomial_algebra, irreducible_paths  # noqa: F401
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass, including the named-cycle assertions.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/core src/quiverlab/combinat/quiver.py tests/core/test_monomial.py
git commit -m "feat: monomial kQ/I with automaton-certified finite dimension"
```

---

### Task 12: Normalized bar complex — Hochschild cohomology dimensions

**Files:**
- Create: `src/quiverlab/hochschild/__init__.py`, `src/quiverlab/hochschild/bar.py`, `src/quiverlab/hochschild/table.py`, `tests/hochschild/__init__.py` (empty), `tests/hochschild/test_bar_cohomology.py`
- Modify: `src/quiverlab/core/algebra.py` (add `hochschild_cohomology` method)

**Interfaces:**
- Consumes: `Algebra` (unit-adapted), `linalg.rank`, `DepthLimitError`.
- Produces: **`HHTable`** (`.dims: list[int]`, `.kind: str` in `{"HH^", "HH_", "HC_"}`, `.top: int`, pretty `__repr__`, `__getitem__`); **`hochschild_cohomology_dims(A, top, max_cells=4_000_000) -> HHTable`**; method **`A.hochschild_cohomology(top)`**. Cochain spaces: normalized `C^n = Hom(Abar^{⊗n}, A)`, `dim = m(m-1)^n`; differential per the standard formula with products projected to `Abar` in middle terms.

- [ ] **Step 1: Write the failing test**

`tests/hochschild/test_bar_cohomology.py`:

```python
import pytest
from quiverlab import CC, GF, Quiver
from quiverlab.errors import DepthLimitError


def _dual(field):
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=["x^2"], field=field)


def test_dual_numbers_char0():
    # classical: HH^0 = 2 (commutative), HH^n = 1 for n >= 1 in char 0
    A = _dual(CC)
    t = A.hochschild_cohomology(4)
    assert t.dims == [2, 1, 1, 1, 1]
    assert t.kind == "HH^"


def test_dual_numbers_char2_pathology():
    # char 2: every differential vanishes -> HH^n = 2 for all n
    A = _dual(GF(2))
    assert A.hochschild_cohomology(4).dims == [2, 2, 2, 2, 2]


def test_dual_numbers_gf4_matches_gf2():
    A = _dual(GF(4))
    assert A.hochschild_cohomology(3).dims == [2, 2, 2, 2]


def test_semisimple_k_times_k():
    Q = Quiver(vertices=[1, 2], arrows={})
    A = Q.algebra()
    assert A.hochschild_cohomology(3).dims == [2, 0, 0, 0]


def test_hereditary_kA2():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    A = Q.algebra()
    assert A.hochschild_cohomology(3).dims == [1, 0, 0, 0]


def test_guard_fails_loudly():
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x^4"], field=GF(3))
    with pytest.raises(DepthLimitError):
        A.hochschild_cohomology(30, max_cells=1000)


def test_repr_is_a_table():
    A = _dual(CC)
    s = repr(A.hochschild_cohomology(2))
    assert "HH^0" in s and "HH^2" in s
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/hochschild/test_bar_cohomology.py -q`
Expected: FAIL — `ImportError` / missing method.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/hochschild/table.py`:

```python
"""Human-readable dimension tables (spec §3.5)."""


class HHTable:
    def __init__(self, dims, kind, algebra_repr, engine="normalized bar complex"):
        self.dims = list(dims)
        self.kind = kind
        self.top = len(dims) - 1
        self.algebra_repr = algebra_repr
        self.engine = engine

    def __getitem__(self, n):
        return self.dims[n]

    def __iter__(self):
        return iter(self.dims)

    def __eq__(self, other):
        if isinstance(other, HHTable):
            return self.kind == other.kind and self.dims == other.dims
        return NotImplemented

    def __repr__(self):
        head = f"{self.kind}n dimensions for {self.algebra_repr} (engine: {self.engine})"
        cells = "  ".join(f"{self.kind}{n} = {d}" for n, d in enumerate(self.dims))
        return head + "\n" + cells
```

`src/quiverlab/hochschild/bar.py`:

```python
"""Normalized bar (co)chain complexes over any Domain (spec §5, component 5:
the 'bar' oracle backend — exponential, small algebras only, always exact)."""
import itertools

from quiverlab.errors import DepthLimitError
from quiverlab.fields.linalg import rank
from quiverlab.hochschild.table import HHTable

_GUARD_HINT = ("the bar oracle is exponential; deeper engines (Bardzell, minimal, "
               "Chouhy-Solotar) arrive in later phases — raise max_cells only if you "
               "know what you are doing")


def _abar_tuples(m, n):
    return list(itertools.product(range(1, m), repeat=n))


def _cochain_basis(m, n):
    return [(s, J) for s in range(m) for J in _abar_tuples(m, n)]


def _check_cells(rows, cols, max_cells, what):
    if rows * cols > max_cells:
        raise DepthLimitError(
            f"{what}: differential matrix would have {rows} x {cols} entries "
            f"(> max_cells = {max_cells})",
            hint=_GUARD_HINT,
        )


def coboundary_matrix(A, n, max_cells):
    """Matrix of d: C^n -> C^{n+1} for unit-adapted A. Rows: C^{n+1} basis, cols: C^n."""
    dom = A.domain
    m = A.dim
    cols = _cochain_basis(m, n)
    rows = _cochain_basis(m, n + 1)
    _check_cells(len(rows), len(cols), max_cells, f"bar coboundary d^{n}")
    row_index = {b: i for i, b in enumerate(rows)}
    D = [[dom.zero()] * len(cols) for _ in range(len(rows))]

    def bump(t, K, ci, val):
        if not dom.is_zero(val):
            r = row_index[(t, K)]
            D[r][ci] = dom.add(D[r][ci], val)

    for ci, (s, J) in enumerate(cols):
        for K in _abar_tuples(m, n + 1):
            # term 0: b_{K0} * f(K[1:])
            if K[1:] == J:
                vec = A.T[K[0]][s]
                for t in range(m):
                    bump(t, K, ci, vec[t])
            # middle terms i = 1..n: (-1)^i f(..., proj(K[i-1] K[i]), ...)
            for i in range(1, n + 1):
                pre, post = K[: i - 1], K[i + 1:]
                if pre == J[: i - 1] and post == J[i - 1 + 1:]:
                    x = J[i - 1]
                    coef = A.T[K[i - 1]][K[i]][x]  # Abar component only (x >= 1)
                    if not dom.is_zero(coef):
                        val = coef if i % 2 == 0 else dom.neg(coef)
                        bump(s, K, ci, val)
            # last term: (-1)^{n+1} f(K[:n]) * b_{K[n]}
            if K[:n] == J:
                vec = A.T[s][K[n]]
                for t in range(m):
                    val = vec[t] if (n + 1) % 2 == 0 else dom.neg(vec[t])
                    bump(t, K, ci, val)
    return D, len(cols), len(rows)


def hochschild_cohomology_dims(A, top, max_cells=4_000_000):
    B = A.unit_adapted()
    dom = B.domain
    m = B.dim
    ranks = []  # rank of d^n for n = 0..top
    for n in range(top + 1):
        D, ncols, nrows = coboundary_matrix(B, n, max_cells)
        ranks.append(rank(D, dom) if nrows and ncols else 0)
    dims = []
    for n in range(top + 1):
        cn = m * (m - 1) ** n
        prev = ranks[n - 1] if n >= 1 else 0
        dims.append(cn - ranks[n] - prev)
    return HHTable(dims, "HH^", repr(A).splitlines()[0])
```

`src/quiverlab/hochschild/__init__.py`:

```python
from quiverlab.hochschild.bar import hochschild_cohomology_dims  # noqa: F401
from quiverlab.hochschild.table import HHTable  # noqa: F401
```

Add to `src/quiverlab/core/algebra.py` (method on `Algebra`, at the end of the class):

```python
    def hochschild_cohomology(self, top, max_cells=4_000_000):
        """Dimensions of HH^0..HH^top via the normalized bar complex (exact)."""
        from quiverlab.hochschild.bar import hochschild_cohomology_dims
        return hochschild_cohomology_dims(self, top, max_cells=max_cells)
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass. The dual-numbers values (2,1,1,1,1) / (2,2,2,2,2) are hand-derived anchors — if they fail, the bug is in the differential, not the test.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/hochschild tests/hochschild src/quiverlab/core/algebra.py
git commit -m "feat: Hochschild cohomology dims via normalized bar complex"
```

---

### Task 13: Bar chains — Hochschild homology dimensions

**Files:**
- Create: `tests/hochschild/test_bar_homology.py`
- Modify: `src/quiverlab/hochschild/bar.py`, `src/quiverlab/hochschild/__init__.py`, `src/quiverlab/core/algebra.py`

**Interfaces:**
- Consumes: same as Task 12.
- Produces: **`hochschild_homology_dims(A, top, max_cells=4_000_000) -> HHTable`** (`kind="HH_"`); method **`A.hochschild_homology(top)`**. Chains: `C_n = A ⊗ Abar^{⊗n}`, boundary with the cyclic last term `(-1)^n (a_n a_0) ⊗ ...`.

- [ ] **Step 1: Write the failing test**

`tests/hochschild/test_bar_homology.py`:

```python
from quiverlab import CC, GF, Quiver


def _dual(field):
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=["x^2"], field=field)


def test_dual_numbers_char0_homology():
    # HH_0 = 2, HH_n = 1 for n >= 1 in char 0
    assert _dual(CC).hochschild_homology(4).dims == [2, 1, 1, 1, 1]


def test_dual_numbers_char2_homology():
    assert _dual(GF(2)).hochschild_homology(4).dims == [2, 2, 2, 2, 2]


def test_kA2_homology_vanishes_positively():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    A = Q.algebra()
    # HH_0 = k^{#vertices} for an acyclic monomial algebra; higher vanish (hereditary, acyclic)
    assert A.hochschild_homology(3).dims == [2, 0, 0, 0]


def test_symmetric_algebra_duality_smoke():
    # k[x]/(x^2) is symmetric: HH^n and HH_n dimensions agree in every characteristic
    for field in (CC, GF(2), GF(3)):
        A = _dual(field)
        assert A.hochschild_cohomology(3).dims == A.hochschild_homology(3).dims
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/hochschild/test_bar_homology.py -q`
Expected: FAIL — no `hochschild_homology`.

- [ ] **Step 3: Write minimal implementation**

Append to `src/quiverlab/hochschild/bar.py`:

```python
def boundary_matrix(A, n, max_cells):
    """Matrix of b: C_n -> C_{n-1}, n >= 1, for unit-adapted A.
    C_n basis: (s, J) with s in 0..m-1 (the A slot), J in Abar^{⊗n}."""
    dom = A.domain
    m = A.dim
    cols = _cochain_basis(m, n)        # same index shape: (s, J)
    rows = _cochain_basis(m, n - 1)
    _check_cells(len(rows), len(cols), max_cells, f"bar boundary b_{n}")
    row_index = {b: i for i, b in enumerate(rows)}
    D = [[dom.zero()] * len(cols) for _ in range(len(rows))]

    def bump(t, K, ci, val):
        if not dom.is_zero(val):
            r = row_index[(t, K)]
            D[r][ci] = dom.add(D[r][ci], val)

    for ci, (s, J) in enumerate(cols):
        # term 0 (+): (b_s b_{J0}) ⊗ J[1:]
        vec = A.T[s][J[0]]
        for t in range(m):
            bump(t, J[1:], ci, vec[t])
        # middle terms i = 1..n-1, sign (-1)^i: merge J[i-1], J[i], project to Abar
        for i in range(1, n):
            merged = A.T[J[i - 1]][J[i]]
            for x in range(1, m):
                coef = merged[x]
                if not dom.is_zero(coef):
                    val = coef if i % 2 == 0 else dom.neg(coef)
                    bump(s, J[: i - 1] + (x,) + J[i + 1:], ci, val)
        # last term, sign (-1)^n: (b_{J[n-1]} b_s) ⊗ J[:n-1]
        vec = A.T[J[n - 1]][s]
        for t in range(m):
            val = vec[t] if n % 2 == 0 else dom.neg(vec[t])
            bump(t, J[: n - 1], ci, val)
    return D, len(cols), len(rows)


def hochschild_homology_dims(A, top, max_cells=4_000_000):
    B = A.unit_adapted()
    dom = B.domain
    m = B.dim
    ranks = [0]  # rank of b_n, with b_0 = 0
    for n in range(1, top + 2):
        D, ncols, nrows = boundary_matrix(B, n, max_cells)
        ranks.append(rank(D, dom) if nrows and ncols else 0)
    dims = []
    for n in range(top + 1):
        cn = m * (m - 1) ** n
        dims.append(cn - ranks[n] - ranks[n + 1])
    return HHTable(dims, "HH_", repr(A).splitlines()[0])
```

Append to `src/quiverlab/hochschild/__init__.py`:

```python
from quiverlab.hochschild.bar import hochschild_homology_dims  # noqa: F401
```

Add to `src/quiverlab/core/algebra.py` (method on `Algebra`):

```python
    def hochschild_homology(self, top, max_cells=4_000_000):
        """Dimensions of HH_0..HH_top via the normalized bar complex (exact)."""
        from quiverlab.hochschild.bar import hochschild_homology_dims
        return hochschild_homology_dims(self, top, max_cells=max_cells)
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/hochschild tests/hochschild src/quiverlab/core/algebra.py
git commit -m "feat: Hochschild homology dims via normalized bar chains"
```

---

### Task 14: Starter families + flat exports

**Files:**
- Create: `src/quiverlab/families/__init__.py`, `src/quiverlab/families/basic.py`, `tests/families/__init__.py` (empty), `tests/families/test_basic.py`
- Modify: `src/quiverlab/__init__.py`

**Interfaces:**
- Consumes: `Quiver.algebra`, fields.
- Produces: **`truncated_polynomial(n, field=CC) -> Algebra`** (k[x]/(x^n), loop quiver, label `x`) and **`linear_path_algebra(n, field=CC) -> Algebra`** (path algebra of the linear A_n quiver `1 -> 2 -> ... -> n`, no relations). Both exported flat. These are building blocks; the full spec §3.4 catalog is Plan 06.

- [ ] **Step 1: Write the failing test**

`tests/families/test_basic.py`:

```python
import pytest
from quiverlab import GF, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError


def test_truncated_polynomial_dims():
    A = truncated_polynomial(4)
    assert A.dim == 4
    A2 = truncated_polynomial(2, field=GF(5))
    assert A2.hochschild_cohomology(2).dims == [2, 1, 1]  # char 5 behaves like char 0 here


def test_linear_path_algebra_dim():
    # kA_n has dimension n(n+1)/2
    assert linear_path_algebra(2).dim == 3
    assert linear_path_algebra(4).dim == 10


def test_bad_arguments_fail_loudly():
    with pytest.raises(QuiverlabError):
        truncated_polynomial(1)
    with pytest.raises(QuiverlabError):
        linear_path_algebra(0)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/families/test_basic.py -q`
Expected: FAIL — `ImportError`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/families/basic.py`:

```python
"""Starter algebra families (building blocks; the full catalog is Plan 06)."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError


def truncated_polynomial(n: int, field=None):
    """k[x]/(x^n) as the one-loop quiver with the monomial relation x^n, n >= 2."""
    if not (isinstance(n, int) and not isinstance(n, bool) and n >= 2):
        raise QuiverlabError(f"truncated_polynomial({n!r}): need an integer n >= 2",
                             hint="n = 1 is the ground field; use n >= 2")
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=[f"x^{n}"], field=field)


def linear_path_algebra(n: int, field=None):
    """The hereditary path algebra of 1 -> 2 -> ... -> n (dimension n(n+1)/2)."""
    if not (isinstance(n, int) and not isinstance(n, bool) and n >= 1):
        raise QuiverlabError(f"linear_path_algebra({n!r}): need an integer n >= 1",
                             hint="e.g. linear_path_algebra(3)")
    arrows = {f"a{i}": (i, i + 1) for i in range(1, n)}
    Q = Quiver(vertices=list(range(1, n + 1)), arrows=arrows)
    return Q.algebra(relations=[], field=field)
```

`src/quiverlab/families/__init__.py`:

```python
from quiverlab.families.basic import linear_path_algebra, truncated_polynomial  # noqa: F401
```

Append to `src/quiverlab/__init__.py`:

```python
from quiverlab.families import linear_path_algebra, truncated_polynomial  # noqa: E402,F401
```

- [ ] **Step 4: Run tests to verify pass**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/families tests/families src/quiverlab/__init__.py
git commit -m "feat: starter families truncated_polynomial and linear_path_algebra"
```

---

### Task 15: Quickstart acceptance test + README

**Files:**
- Create: `tests/test_quickstart.py`
- Modify: `README.md`

**Interfaces:**
- Consumes: everything above.
- Produces: the README's example, executed verbatim in CI; the Plan-01 exit criterion.

- [ ] **Step 1: Write the failing test**

`tests/test_quickstart.py`:

```python
"""The README example, verbatim: this is the Plan-01 acceptance test.
A non-coder builds a custom monomial algebra and reads off exact Hochschild
dimensions with a characteristic sweep."""
from quiverlab import CC, GF, Quiver


def test_readme_example():
    Q = Quiver(vertices=[1, 2, 3],
               arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
    A = Q.algebra(relations=["a*b"], field=CC)

    assert A.dim == 6            # e1, e2, e3, a, b, c
    assert A.basis_labels == ["e_1", "e_2", "e_3", "a", "b", "c"]

    hh = A.hochschild_cohomology(3)
    hh_gf2 = Q.algebra(relations=["a*b"], field=GF(2)).hochschild_cohomology(3)

    # both are exact, certified answers; printing them shows readable tables
    assert hh.dims[0] >= 1 and len(hh.dims) == 4
    assert len(hh_gf2.dims) == 4
    assert "HH^0" in repr(hh)


def test_characteristic_sweep_dual_numbers():
    from quiverlab import truncated_polynomial

    table = {}
    for name, field in [("CC", CC), ("GF(2)", GF(2)), ("GF(3)", GF(3)), ("GF(4)", GF(4))]:
        table[name] = truncated_polynomial(2, field=field).hochschild_cohomology(4).dims
    assert table["CC"] == [2, 1, 1, 1, 1]
    assert table["GF(2)"] == [2, 2, 2, 2, 2]      # the char-2 pathology, exactly
    assert table["GF(3)"] == [2, 1, 1, 1, 1]
    assert table["GF(4)"] == [2, 2, 2, 2, 2]
```

- [ ] **Step 2: Run test to verify it fails or passes honestly**

Run: `python -m pytest tests/test_quickstart.py -q`
Expected: PASS if Tasks 1–14 are correct (this task adds no code). If it fails, a previous task has a bug — fix there, not here.

- [ ] **Step 3: Write the README**

Replace `README.md` with:

````markdown
# quiverlab

**Quivers with relations and Hochschild theory, exactly, for algebraists.**

quiverlab computes with finite-dimensional algebras presented by quivers with
relations, over the complex numbers (exactly — no floating point, ever) and
over all finite fields. Floats fail loudly by design.

## Install (development preview)

```bash
pip install -e '.[dev]'
```

## Three lines to a Hochschild table

```python
from quiverlab import Quiver, CC, GF

Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
A = Q.algebra(relations=["a*b"], field=CC)
print(A.hochschild_cohomology(3))
```

Change one word to move the characteristic:

```python
A2 = Q.algebra(relations=["a*b"], field=GF(2))
print(A2.hochschild_cohomology(3))
```

## The classic characteristic pathology, in one loop

```python
from quiverlab import truncated_polynomial, CC, GF

for field in (CC, GF(2), GF(3)):
    print(field, truncated_polynomial(2, field=field).hochschild_cohomology(4).dims)
# CC     [2, 1, 1, 1, 1]
# GF(2)  [2, 2, 2, 2, 2]
# GF(3)  [2, 1, 1, 1, 1]
```

## Status

Foundations phase (Plan 01). Monomial presentations, exact fields, bar-complex
Hochschild (co)homology. Coming next (see `docs/plans/ROADMAP.md`): the hanlab
deep engines, general relations via noncommutative Groebner bases, the first
full Chouhy–Solotar resolution, cup products and Gerstenhaber brackets, module
Ext, drawing and TikZ export, worked-steps PDFs, and an optional QPA backend.

MIT © 2026 Marco Armenta
````

- [ ] **Step 4: Run the full suite**

Run: `python -m pytest -q`
Expected: all pass (float gate included).

- [ ] **Step 5: Commit**

```bash
git add tests/test_quickstart.py README.md
git commit -m "feat: quickstart acceptance test and README"
```

---

## Plan self-review notes (done at writing time)

- Spec coverage of this phase: §3.2 (CC, GF), §3.3 (monomial subset + loud
  NotImplementedError boundary), §4.1–§4.3, §5 components 1, 2, 4 and the `bar`
  backend of 5, §7 (all six error classes exercised in tests). Deliberately
  deferred (roadmap): components 3, 6–12, §3.4–§3.9 beyond the starters, §6, §8
  rings 2–4, §9–§11.
- The dual-numbers HH values used as test anchors are hand-derived in the plan's
  design notes and double-checked against the classical literature values.
- Type consistency: `Domain`, `linalg` names, `Relation.terms`, `Algebra`
  attributes, and `HHTable` fields are used identically across Tasks 3–15 and
  frozen in ROADMAP.md for later plans.
