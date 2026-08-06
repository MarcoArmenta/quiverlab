# Plan 36: Macaulay2 Oracle Bridge — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Macaulay2 becomes quiverlab's fifth oracle class (`m2`): an
independent cross-check for single-vertex graded `kQ/I` dimensions/Hilbert
series (via `AssociativeAlgebras`, a genuinely different nc-Gröbner engine)
and for commutative examples' Ext/Betti data (via `freeResolution`), wired
into the Plan-32 class audit and the verification page.

**Architecture:** `src/quiverlab/m2/` mirrors `src/quiverlab/qpa/`
(`session.py` availability + run, `scripts.py` pure string builders,
`crosscheck.py` typed comparisons reusing `CrosscheckReport`) — except the
transport is a **subprocess** (`M2 --script` per call, sentinel-line stdout
parsing) instead of in-process libgap. Tests live in a new extras-gated
bucket `tests/m2/`, exactly parallel to `tests/qpa/`.

**Tech Stack:** Python stdlib `subprocess`/`shutil.which` (no new Python
deps — M2 is a system binary, not pip-installable, so there is **no `[m2]`
extra**), Macaulay2 ≥ 1.24 with bundled `AssociativeAlgebras` + `Complexes`.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- No floats in `src/` (AST gate `tests/test_no_floats.py`): all parsed M2
  output must go through exact-int parsing (mirror `qpa/crosscheck.py::_read_int_list`).
- Every text-mode `subprocess.run` passes `encoding="utf-8"`
  (`tests/trace/test_artifact_encoding.py` scans the shipping tree).
- The `m2` marker is a bucket AND the oracle class — tests in `tests/m2/`
  are **never** double-marked with `oracle_*` markers (Plan-32 rule, same as
  `qpa`).
- Honest scope everywhere: M2 sees no multi-vertex algebras, no Hochschild —
  builders raise `QuiverlabError` on multi-vertex input, docs say so.
- Conventional commits; suite green at every commit; work on branch
  `plan-36-m2-oracle` off `dev`.

**Prerequisite (manual, before Task 4):** a local M2:
`brew tap Macaulay2/tap && brew install M2` (macOS). If unavailable, Tasks
4–6's live steps run in CI only — the plan still completes (skip predicate).

---

### Task 1: `m2/session.py` — availability probe, runner, skip predicate

**Files:**
- Create: `src/quiverlab/m2/__init__.py`
- Create: `src/quiverlab/m2/session.py`
- Modify: `src/quiverlab/errors.py` (add `M2UnavailableError` beside `QpaUnavailableError`)
- Test: `tests/m2/test_session_guard.py`, `tests/m2/__init__.py` (empty)

**Interfaces:**
- Produces: `m2_available() -> bool` (cached), `require_m2() -> None`
  (raises `M2UnavailableError` with a brew/apt hint),
  `run_script(script: str, timeout: int = 120) -> str` (writes the script to
  a temp file, runs `[M2, "--script", path]`, returns stdout; raises
  `M2UnavailableError` if absent, `RuntimeError` on nonzero exit with stderr
  in the message), `should_skip_m2() -> bool`
  (`not m2_available() and os.environ.get("QUIVERLAB_REQUIRE_M2") != "1"`),
  `m2_version() -> str`.
- Consumes: nothing from other tasks.

- [ ] **Step 1: Write the failing tests** (no M2 needed — these are the
  always-run guards, `pytestmark = pytest.mark.fast` like
  `tests/qpa/test_session_guard.py`)

```python
# tests/m2/test_session_guard.py
"""Session guards for the Macaulay2 bridge -- run everywhere, no M2 needed."""
import os
import pytest

from quiverlab.errors import M2UnavailableError
from quiverlab.m2 import session

pytestmark = pytest.mark.fast


def test_unavailable_raises_with_hint(monkeypatch):
    monkeypatch.setattr(session, "_which_m2", lambda: None)
    session.m2_available.cache_clear()
    with pytest.raises(M2UnavailableError) as e:
        session.require_m2()
    assert "brew" in str(e.value) or "apt" in str(e.value)
    session.m2_available.cache_clear()


def test_skip_predicate_env_override(monkeypatch):
    monkeypatch.setattr(session, "_which_m2", lambda: None)
    session.m2_available.cache_clear()
    monkeypatch.delenv("QUIVERLAB_REQUIRE_M2", raising=False)
    assert session.should_skip_m2() is True
    monkeypatch.setenv("QUIVERLAB_REQUIRE_M2", "1")
    assert session.should_skip_m2() is False   # CI: absence must FAIL, not skip
    session.m2_available.cache_clear()


def test_import_does_not_probe(monkeypatch):
    # importing quiverlab.m2 must not shell out (mirror qpa laziness)
    import importlib
    import quiverlab.m2
    importlib.reload(quiverlab.m2)   # no error even with no M2 on PATH
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/m2/test_session_guard.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.m2` / `ImportError: M2UnavailableError`

- [ ] **Step 3: Implement**

In `src/quiverlab/errors.py`, next to `QpaUnavailableError` (same base
class as it uses):

```python
class M2UnavailableError(QuiverlabError):
    """The Macaulay2 binary is not on PATH (or refused to run)."""
```

`src/quiverlab/m2/session.py`:

```python
"""Macaulay2 subprocess session: availability probe + script runner.

Unlike the QPA bridge (in-process libgap), M2 is driven as a subprocess:
each call writes the script to a temp file and runs ``M2 --script file``.
Importing this module never probes the binary (lazy, like qpa.session).
"""
from __future__ import annotations

import functools
import os
import shutil
import subprocess
import tempfile

from quiverlab.errors import M2UnavailableError

_HINT = ("Macaulay2 not found. Install it: macOS `brew tap Macaulay2/tap && "
         "brew install M2`; Ubuntu `sudo add-apt-repository ppa:macaulay2/macaulay2 "
         "&& sudo apt install macaulay2`. The m2 test bucket skips without it "
         "unless QUIVERLAB_REQUIRE_M2=1.")


def _which_m2():
    return shutil.which("M2")


@functools.lru_cache(maxsize=1)
def m2_available() -> bool:
    """True iff the M2 binary is on PATH and answers --version. Cached."""
    exe = _which_m2()
    if exe is None:
        return False
    try:
        out = subprocess.run([exe, "--version"], capture_output=True,
                             encoding="utf-8", timeout=30)
    except (OSError, subprocess.SubprocessError):
        return False
    return out.returncode == 0


def require_m2() -> None:
    """Raise M2UnavailableError (with a fix-it hint) unless M2 is live."""
    if not m2_available():
        raise M2UnavailableError(_HINT)


def m2_version() -> str:
    """The M2 version string (requires M2)."""
    require_m2()
    out = subprocess.run([_which_m2(), "--version"], capture_output=True,
                         encoding="utf-8", timeout=30)
    return out.stdout.strip() or out.stderr.strip()


def run_script(script: str, timeout: int = 120) -> str:
    """Run an M2 script headless; return raw stdout. Loud on failure."""
    require_m2()
    with tempfile.NamedTemporaryFile("w", suffix=".m2", delete=False,
                                     encoding="utf-8") as f:
        f.write(script)
        path = f.name
    try:
        out = subprocess.run([_which_m2(), "--script", path],
                             capture_output=True, encoding="utf-8",
                             timeout=timeout)
    finally:
        os.unlink(path)
    if out.returncode != 0:
        raise RuntimeError(f"M2 script failed (exit {out.returncode}):\n"
                           f"{out.stderr}")
    return out.stdout


def should_skip_m2() -> bool:
    """Skip predicate for tests/m2 (collection-time, like should_skip_qpa):
    skip when M2 is absent -- unless QUIVERLAB_REQUIRE_M2=1 (CI), where
    absence must fail loudly, not skip."""
    return not m2_available() and os.environ.get("QUIVERLAB_REQUIRE_M2") != "1"
```

`src/quiverlab/m2/__init__.py`:

```python
"""Macaulay2 oracle bridge (subprocess). Importing this package does NOT
probe the M2 binary -- everything is lazy, mirroring quiverlab.qpa."""
from quiverlab.m2.session import m2_available, require_m2, should_skip_m2

__all__ = ["m2_available", "require_m2", "should_skip_m2", "crosscheck"]


def crosscheck(algebra, what, *args, **kwargs):   # lazy import, like qpa
    from quiverlab.m2.crosscheck import crosscheck as _cc
    return _cc(algebra, what, *args, **kwargs)
```

(`crosscheck` lands in Task 3; the lazy wrapper keeps this import-safe now —
add a placeholder module raise if needed, or defer the `__init__` export to
Task 3. Prefer: export only session names now, add `crosscheck` in Task 3.)

- [ ] **Step 4: Run tests to verify they pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/m2/test_session_guard.py -v`
Expected: PASS (3 tests). Note: until Task 4 wires the bucket, these
collect as `fast` via the explicit `pytestmark` — fine.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/m2/ src/quiverlab/errors.py tests/m2/
git commit -m "feat(m2): session layer -- availability probe, --script runner, skip predicate"
```

---

### Task 2: `m2/scripts.py` — pure string builders (no M2 needed)

**Files:**
- Create: `src/quiverlab/m2/scripts.py`
- Test: `tests/m2/test_scripts.py`

**Interfaces:**
- Consumes: an `Algebra` built as
  `Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(relations=[...], field=GF(p))`;
  reads exactly the public accessors the QPA builders read
  (`algebra.quiver` — vertices/arrows, `algebra.relations` — relation
  strings like `"x*x"`, `"x*y-y*x"`, `"x^4"`, `algebra.domain.characteristic`).
  Read `src/quiverlab/qpa/scripts.py::quiver_and_algebra_script` first and
  reuse its accessor idioms verbatim.
- Produces:
  - `graded_dims_script(A, top: int) -> str` — M2 code that builds the
    single-vertex free-algebra quotient over GF(p) and prints one sentinel
    line `<<QL>> n d_n` per degree `0..top`.
  - `commutative_ext_script(p: int, variables: list[str], relations: list[str], top: int) -> str`
    — M2 code for `ZZ/p[vars]/(relations)`, printing `<<QL>> n b_n` where
    `b_n = rank` of the n-th term of the minimal free resolution of the
    residue field.
  - `parse_sentinels(stdout: str) -> list[int]` — extracts `<<QL>> n v`
    lines, checks n = 0,1,2,… contiguity, exact `int()` on v (ValueError on
    anything non-integer — the no-floats discipline at the bridge boundary).
  - `SENTINEL = "<<QL>>"`.

- [ ] **Step 1: Write the failing tests** (`pytestmark = pytest.mark.fast`,
  mirroring `tests/qpa/test_scripts.py` — string assertions only)

```python
# tests/m2/test_scripts.py
"""M2 script builders are pure string generation -- exercised without M2."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.m2 import scripts

pytestmark = pytest.mark.fast


def _dual_numbers():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(7))


def test_graded_dims_script_shape():
    s = scripts.graded_dims_script(_dual_numbers(), top=4)
    assert "ZZ/7" in s
    assert "x*x" in s
    assert scripts.SENTINEL in s
    assert "ncBasis" in s


def test_multi_vertex_refused():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))
    with pytest.raises(QuiverlabError, match="single-vertex"):
        scripts.graded_dims_script(A, top=3)


def test_caret_relations_translate():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x^3"], field=GF(5))
    s = scripts.graded_dims_script(A, top=3)
    assert "x^3" in s        # M2 shares the caret power syntax


def test_parse_sentinels_roundtrip():
    out = "junk\n<<QL>> 0 1\n<<QL>> 1 2\n<<QL>> 2 1\nnoise"
    assert scripts.parse_sentinels(out) == [1, 2, 1]


def test_parse_sentinels_rejects_noninteger():
    with pytest.raises(ValueError):
        scripts.parse_sentinels("<<QL>> 0 1.5")


def test_parse_sentinels_rejects_gap_in_degrees():
    with pytest.raises(ValueError):
        scripts.parse_sentinels("<<QL>> 0 1\n<<QL>> 2 1")


def test_commutative_ext_script_shape():
    s = scripts.commutative_ext_script(7, ["x", "y"], ["x^2", "y^2"], top=6)
    assert "ZZ/7[x,y]" in s.replace(" ", "")
    assert "freeResolution" in s and "LengthLimit" in s
    assert scripts.SENTINEL in s
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `... -m pytest tests/m2/test_scripts.py -v`
Expected: FAIL — `ImportError: cannot import name 'scripts'`

- [ ] **Step 3: Implement `src/quiverlab/m2/scripts.py`**

```python
"""Pure M2 script builders. String generation only -- no subprocess here.

Translation notes (single-vertex kQ/I over GF(p) -> AssociativeAlgebras):
loops become free-algebra generators; quiverlab relation strings already use
`*` and `^`, which M2 shares, so relations pass through verbatim after a
whitespace strip. Multi-vertex input is refused loudly: M2 has no
vertex-idempotent notion (Plan-36 honest scope).
"""
from __future__ import annotations

from quiverlab.errors import QuiverlabError

SENTINEL = "<<QL>>"


def _single_vertex_data(A):
    quiver = A.quiver
    if len(quiver.vertices) != 1:
        raise QuiverlabError(
            "M2 oracle scope is single-vertex algebras only: "
            "AssociativeAlgebras has no quiver/idempotent type "
            f"(got {len(quiver.vertices)} vertices).")
    p = A.domain.characteristic
    if not p:
        raise QuiverlabError("M2 oracle batteries run over GF(p); "
                             "characteristic 0 input not supported here.")
    gens = sorted(quiver.arrows)          # loop names, deterministic order
    rels = [r.replace(" ", "") for r in A.relations]
    return p, gens, rels


def graded_dims_script(A, top: int) -> str:
    """M2 script printing `<<QL>> n dim_k(B_n)` for n = 0..top, where B is
    the degree-`top`-truncated nc Groebner quotient of the free algebra."""
    p, gens, rels = _single_vertex_data(A)
    gen_list = ",".join(gens)
    ideal = ",".join(rels) if rels else ""
    lines = [
        f"kk = ZZ/{p}",
        f"F = kk<|{gen_list}|>",
    ]
    if ideal:
        lines += [
            f"I = ideal {{{ideal}}}",
            f"Igb = NCGB(I, {2 * top + 2})",
            "B = F/I",
        ]
    else:
        lines += ["B = F"]
    lines += [
        f'for n from 0 to {top} do print("{SENTINEL} " | toString n | " " '
        "| toString numgens source ncBasis(n, B))",
    ]
    return "\n".join(lines) + "\n"


def commutative_ext_script(p: int, variables, relations, top: int) -> str:
    """M2 script printing `<<QL>> n rank(P_n)` for the minimal graded free
    resolution of the residue field of ZZ/p[variables]/(relations)."""
    vs = ",".join(variables)
    rels = ",".join(r.replace(" ", "") for r in relations)
    return "\n".join([
        "needsPackage \"Complexes\"",
        f"R = ZZ/{p}[{vs}]/({rels})",
        f"C = freeResolution(coker vars R, LengthLimit => {top})",
        f'for n from 0 to {top} do print("{SENTINEL} " | toString n | " " '
        "| toString rank C_n)",
    ]) + "\n"


def parse_sentinels(stdout: str) -> list:
    """Extract the `<<QL>> n v` values in degree order; exact ints only."""
    got = {}
    for line in stdout.splitlines():
        line = line.strip()
        if not line.startswith(SENTINEL):
            continue
        _, n_s, v_s = line.split(None, 2)
        got[int(n_s)] = int(v_s)          # int() raises on "1.5" -- wanted
    if sorted(got) != list(range(len(got))):
        raise ValueError(f"non-contiguous degrees in M2 output: {sorted(got)}")
    return [got[n] for n in range(len(got))]
```

**Implementation note:** confirm the exact accessor names against
`src/quiverlab/qpa/scripts.py` (`algebra.quiver`, `algebra.relations`,
`algebra.domain.characteristic` — adjust to whatever that file actually
uses; it is the source of truth for reading a presentation).

- [ ] **Step 4: Run tests to verify they pass**

Run: `... -m pytest tests/m2/test_scripts.py -v` — Expected: PASS (8 tests)

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/m2/scripts.py tests/m2/test_scripts.py
git commit -m "feat(m2): script builders -- graded dims + commutative ext, sentinel parsing"
```

---

### Task 3: `m2/crosscheck.py` — typed comparisons + dispatcher

**Files:**
- Create: `src/quiverlab/m2/crosscheck.py`
- Modify: `src/quiverlab/m2/__init__.py` (export `crosscheck`)
- Test: `tests/m2/test_crosscheck_offline.py`

**Interfaces:**
- Consumes: `session.run_script`, `scripts.*` (Tasks 1–2);
  `quiverlab.qpa.crosscheck.CrosscheckReport` (reused dataclass — fields
  `what, ours, qpa, agree`, method `assert_agree()`; the `qpa` field name is
  historical, treat it as "theirs");
  `quiverlab.modules.koszul._algebra_graded_matrices(A, through)` → our
  graded dims (`sum(sum(row) for row in C_b)` per degree);
  `A.ext_algebra(top).graded_dims_through(top)` → our Ext(k,k) dims.
- Produces:
  - `crosscheck_graded_dims(A, top=6) -> CrosscheckReport` — ours = totals
    of `_algebra_graded_matrices(A, top)`, theirs =
    `parse_sentinels(run_script(graded_dims_script(A, top)))`.
  - `crosscheck_commutative_ext(A, variables, relations, top=6) -> CrosscheckReport`
    — ours = `A.ext_algebra(top=top).graded_dims_through(top)`, theirs =
    the M2 Betti ranks. The caller asserts commutativity by supplying the
    polynomial presentation; the function first verifies every commutator
    `x*y-y*x` is in (the ideal closure of) `A.relations` textually — if a
    commutator string is not among the normalized relations, raise
    `QuiverlabError("not presented as commutative")`.
  - `crosscheck(algebra, what: str, *args, **kwargs)` dispatcher with
    `what in {"graded_dims", "commutative_ext"}` (loud `KeyError`-style
    `QuiverlabError` otherwise), mirroring `qpa/crosscheck.py:332`.

- [ ] **Step 1: Write the failing offline tests** (fast; monkeypatch
  `session.run_script` so no M2 is needed)

```python
# tests/m2/test_crosscheck_offline.py
"""Crosscheck plumbing with a canned M2 transcript -- no M2 binary."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.m2 import crosscheck as cc
from quiverlab.m2 import session

pytestmark = pytest.mark.fast


def _dual_numbers():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(7))


def test_graded_dims_agree_path(monkeypatch):
    # k[x]/(x^2): graded dims 1,1,0,0,0
    monkeypatch.setattr(session, "run_script",
                        lambda s, timeout=120: "<<QL>> 0 1\n<<QL>> 1 1\n"
                                               "<<QL>> 2 0\n<<QL>> 3 0\n<<QL>> 4 0\n")
    rep = cc.crosscheck_graded_dims(_dual_numbers(), top=4)
    assert rep.agree and rep.ours == rep.qpa == [1, 1, 0, 0, 0]
    rep.assert_agree()


def test_disagreement_is_loud(monkeypatch):
    monkeypatch.setattr(session, "run_script",
                        lambda s, timeout=120: "<<QL>> 0 1\n<<QL>> 1 2\n"
                                               "<<QL>> 2 0\n<<QL>> 3 0\n<<QL>> 4 0\n")
    rep = cc.crosscheck_graded_dims(_dual_numbers(), top=4)
    assert not rep.agree
    with pytest.raises(AssertionError):
        rep.assert_agree()


def test_dispatcher_unknown_subject():
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError, match="graded_dims"):
        cc.crosscheck(_dual_numbers(), "hochschild")   # honest: no M2 HH


def test_commutative_guard():
    from quiverlab.errors import QuiverlabError
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y"], field=GF(7))       # free product, NOT comm.
    with pytest.raises(QuiverlabError, match="commutative"):
        cc.crosscheck_commutative_ext(A, ["x", "y"], ["x^2", "y^2"], top=4)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `... -m pytest tests/m2/test_crosscheck_offline.py -v`
Expected: FAIL — no `crosscheck` module

- [ ] **Step 3: Implement `src/quiverlab/m2/crosscheck.py`**

```python
"""Typed M2 comparisons. Reuses the QPA CrosscheckReport container (pure
dataclass; its `qpa` field is read as "the external system's value")."""
from __future__ import annotations

from itertools import combinations

from quiverlab.errors import QuiverlabError
from quiverlab.m2 import scripts, session
from quiverlab.qpa.crosscheck import CrosscheckReport


def _our_graded_dims(A, top):
    from quiverlab.modules.koszul import _algebra_graded_matrices
    mats = _algebra_graded_matrices(A, top)
    return [sum(sum(row) for row in mn) for mn in mats]


def crosscheck_graded_dims(A, top=6):
    ours = _our_graded_dims(A, top)
    theirs = scripts.parse_sentinels(
        session.run_script(scripts.graded_dims_script(A, top)))
    return CrosscheckReport(what=f"graded_dims:0..{top}", ours=ours,
                            qpa=theirs, agree=ours == theirs)


def _assert_commutative_presentation(A):
    rels = {r.replace(" ", "") for r in A.relations}
    gens = sorted(A.quiver.arrows)
    for a, b in combinations(gens, 2):
        if f"{a}*{b}-{b}*{a}" not in rels and f"{b}*{a}-{a}*{b}" not in rels:
            raise QuiverlabError(
                f"not presented as commutative: missing commutator for "
                f"({a},{b}) -- the M2 side would be a different algebra.")


def crosscheck_commutative_ext(A, variables, relations, top=6):
    _assert_commutative_presentation(A)
    ours = A.ext_algebra(top=top).graded_dims_through(top)
    theirs = scripts.parse_sentinels(session.run_script(
        scripts.commutative_ext_script(A.domain.characteristic,
                                       variables, relations, top)))
    return CrosscheckReport(what=f"commutative_ext:0..{top}", ours=ours,
                            qpa=theirs, agree=ours == theirs)


_SUBJECTS = {
    "graded_dims": crosscheck_graded_dims,
    "commutative_ext": crosscheck_commutative_ext,
}


def crosscheck(algebra, what, *args, **kwargs):
    """Dispatch an M2 crosscheck. Valid subjects: graded_dims,
    commutative_ext. Anything else (e.g. hochschild) is refused loudly --
    M2 has no Hochschild theory (honest scope)."""
    if what not in _SUBJECTS:
        raise QuiverlabError(
            f"no M2 oracle for {what!r}; available: "
            f"{sorted(_SUBJECTS)} (M2 has no quiver or Hochschild support).")
    return _SUBJECTS[what](algebra, *args, **kwargs)
```

Add to `src/quiverlab/m2/__init__.py` `__all__` and define the lazy
`crosscheck` wrapper exactly as sketched in Task 1 Step 3.

- [ ] **Step 4: Run tests to verify they pass**

Run: `... -m pytest tests/m2/ -v` — Expected: PASS (all offline tests)

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/m2/ tests/m2/test_crosscheck_offline.py
git commit -m "feat(m2): crosscheck layer -- graded dims + commutative ext vs M2, loud scope refusals"
```

---

### Task 4: bucket + marker + oracle-class wiring

**Files:**
- Modify: `pyproject.toml` (markers list, ~line 95-110)
- Modify: `tests/conftest.py` (`_bucket`, `_BUCKETS`)
- Modify: `tests/release/test_markers.py` (partition check)
- Modify: `tests/release/test_oracle_classes.py` (`CANON`, union slices, extras-guard dir list)
- Modify: `docs/verification.md` (class table + subsystem row + Class-2 subsection + honest scope)

**Interfaces:**
- Consumes: the existing audit mechanics — `CANON` tuple
  (`tests/release/test_oracle_classes.py:36-42`), `_page_counts` parsing of
  ``| ... | `-m <expr>` | <count> |`` rows, `_bucket` in `tests/conftest.py`.
- Produces: `tests/m2/` is its own collection bucket `m2` (like `qpa`);
  `-m m2` is the fifth oracle class; the verification page documents it and
  the audit gates it.

- [ ] **Step 1: Write the wiring** (the failing "test" here IS the release
  gate — run it before touching anything to see the baseline pass, then
  after each edit):

1. `pyproject.toml` markers — add beside the qpa line:
   `"m2: requires a local Macaulay2 binary; the CI M2 job only. Doubles as the fifth oracle class (external M2 recomputes and agrees).",`
2. `tests/conftest.py`:
   - `_BUCKETS = {"fast", "deep", "qpa", "m2"}`
   - in `_bucket`, after the qpa branch: `if top == "m2": return "m2"`
3. `tests/release/test_markers.py` — fold `m2` into the
   disjoint/exhaustive partition exactly as `qpa` is folded (read the file;
   add `_ids("m2")` beside `_ids("qpa")` in both the disjointness pairs and
   the union-equals-everything assertion).
4. `tests/release/test_oracle_classes.py`:
   - `CANON` becomes the 6-tuple: the four existing single classes, then
     `"m2"`, then the union string extended to
     `"oracle_literature or oracle_crossengine or oracle_selfcert or qpa or m2"`.
   - Update the parts/union slicing: parts = `CANON[:5]`, union = `CANON[5]`.
   - `test_no_oracle_markers_in_extras_gated_dirs`: add `"m2"` to the
     scanned dirs tuple (tests/m2 must never carry `oracle_*` marks).
5. `docs/verification.md`:
   - Class table (~line 481): insert before the union row:
     `| Live Macaulay2 | \`-m m2\` | <count> | an independent external system (Macaulay2) recomputes and agrees |`
     and extend the union row's expression + count. Fill `<count>` from the
     live collection in Step 2.
   - Subsystem table (~line 407): add
     ``| `m2/` (Macaulay2 crosscheck) | <n> | <k> m2 + <j> fast | **live Macaulay2** -- nc graded dims (single-vertex); commutative Ext/Betti |``
   - New `### Live Macaulay2 cross-check` subsection under `## Class 2`
     (after the QPA subsection, ~line 339): 2-3 paragraphs — what M2
     recomputes (single-vertex graded dims via an independent nc-Gröbner
     F4 engine; commutative-example Ext dims via `freeResolution`), the
     subprocess transport, and the version pin policy.
   - `## Honest scope` (~line 517): add the entry "Macaulay2 cannot see
     multi-vertex algebras or Hochschild anything — `AssociativeAlgebras`
     has no quiver/idempotent type; multi-vertex and HH stay with QPA +
     theory oracles. The bridge refuses those inputs loudly."

- [ ] **Step 2: Get the live counts and finish the page**

Run: `... -m pytest -q --collect-only -m m2 2>/dev/null | tail -2` and the
union expression likewise; paste the exact numbers into the two tables.
(At this point the m2 count = the live-battery tests that exist so far; the
audit test recomputes at every future merge, so Task-5's additions will
require updating the number again — do the final count AFTER Task 5 and
say so in that task.)

- [ ] **Step 3: Run the release gates**

Run: `... -m pytest tests/release/test_markers.py tests/release/test_oracle_classes.py -q`
Expected: PASS. Also run `... -m pytest tests/m2/ -q` (still green) and a
10-file smoke of fast tests to confirm no collection breakage:
`... -m pytest tests/fields tests/core -q`.

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml tests/conftest.py tests/release/ docs/verification.md
git commit -m "feat(m2): m2 bucket + fifth oracle class -- conftest, markers, Plan-32 audit, verification page"
```

---

### Task 5: live batteries (skip-gated; the actual oracle tests)

**Files:**
- Create: `tests/m2/test_graded_dims_m2.py`
- Create: `tests/m2/test_commutative_ext_m2.py`
- Modify: `docs/verification.md` (final counts from Task 4 Step 2)

**Interfaces:**
- Consumes: `quiverlab.m2.crosscheck` (Task 3), the zoo builders
  (`truncated_polynomial`, `QuantumCI`, `Quiver(...).algebra(...)`).
- Produces: the live m2 oracle battery. Every file header:
  `pytestmark = pytest.mark.skipif(session.should_skip_m2(), reason="Macaulay2 not installed")`.

- [ ] **Step 1: Write the batteries**

```python
# tests/m2/test_graded_dims_m2.py
"""Live M2 battery: graded dims of single-vertex kQ/I over GF(p) --
an independent nc-Groebner engine recomputes our Hilbert data."""
import pytest

from quiverlab import GF, Quiver, truncated_polynomial
from quiverlab.m2 import crosscheck as cc
from quiverlab.m2 import session

pytestmark = pytest.mark.skipif(session.should_skip_m2(),
                                reason="Macaulay2 not installed")


CASES = [
    ("truncated_x5", lambda: truncated_polynomial(5, field=GF(32003)), 8),
    ("dual_numbers", lambda: Quiver([1], {"x": (1, 1)}).algebra(
        relations=["x*x"], field=GF(7)), 6),
    ("two_loops_radsq", lambda: Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "x*y", "y*x", "y*y"], field=GF(5)), 6),
    ("exterior_2", lambda: Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y+y*x"], field=GF(32003)), 6),
    ("quantum_ci_q2", lambda: Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x-2*x*y"], field=GF(32003)), 6),
    ("straddle_monomial", lambda: Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y*x"], field=GF(3)), 8),
]


@pytest.mark.parametrize("name,build,top", CASES, ids=[c[0] for c in CASES])
def test_graded_dims_match_m2(name, build, top):
    cc.crosscheck_graded_dims(build(), top=top).assert_agree()
```

```python
# tests/m2/test_commutative_ext_m2.py
"""Live M2 battery: Ext_A(k,k) dims of commutative examples vs the minimal
free resolution M2 computes -- a fully independent homological route."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.m2 import crosscheck as cc
from quiverlab.m2 import session

pytestmark = pytest.mark.skipif(session.should_skip_m2(),
                                reason="Macaulay2 not installed")


def _comm(relations, extra, p):
    """k[x,y]/(relations) as a quiver algebra: loops + commutator."""
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=relations + ["x*y-y*x"] + extra, field=GF(p))


CASES = [
    ("kxy_x2y2_gf7",   ["x*x", "y*y"],        [], 7, ["x^2", "y^2"], 6),
    ("kxy_x2y2_gf3",   ["x*x", "y*y"],        [], 3, ["x^2", "y^2"], 6),
    ("kxy_x3y2",       ["x*x*x", "y*y"],      [], 5, ["x^3", "y^2"], 6),
    ("kxy_x2y3",       ["x*x", "y*y*y"],      [], 32003, ["x^2", "y^3"], 6),
]


@pytest.mark.parametrize("name,rels,extra,p,m2rels,top", CASES,
                         ids=[c[0] for c in CASES])
def test_ext_k_matches_m2(name, rels, extra, p, m2rels, top):
    A = _comm(rels, extra, p)
    cc.crosscheck_commutative_ext(A, ["x", "y"], m2rels, top=top).assert_agree()
```

- [ ] **Step 2: Run the live battery** (requires local M2 — install per
  the prerequisite; otherwise push and let the CI job from Task 6 be the
  first live run, and say so in the PR/commit message)

Run: `... -m pytest tests/m2/ -m m2 -v`
Expected: PASS. **If any case disagrees, STOP — that is either a real bug
in our Gröbner/Ext stack or an M2 translation error: debug with the
systematic-debugging skill before weakening anything.** Timebox slow M2
GB cases by lowering `top`, never by dropping the case silently.

- [ ] **Step 3: Refresh the audited counts**

Rerun the Task-4 Step-2 collection counts, update the two
`docs/verification.md` tables, run
`... -m pytest tests/release/test_oracle_classes.py -q` → PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/m2/ docs/verification.md
git commit -m "test(m2): live batteries -- 6 nc graded-dims cases + 4 commutative Ext cases"
```

---

### Task 6: CI job + docs + design-mining note

**Files:**
- Create: `.github/workflows/m2.yml` (model: `.github/workflows/qpa.yml` — read it first, copy its trigger/paths/venv steps)
- Create: `docs/plans/2026-08-05-plan-36-m2-design-notes.md`
- Modify: `README.md` (verification section: mention the fifth oracle class)

**Interfaces:**
- Consumes: everything above.
- Produces: a CI job that FAILS (not skips) when M2 is absent, via
  `QUIVERLAB_REQUIRE_M2=1`; the design-mining note that P39/P42 read.

- [ ] **Step 1: Write the workflow** — copy `qpa.yml`'s structure; the
  M2-specific job steps:

```yaml
# .github/workflows/m2.yml (core steps; mirror qpa.yml's triggers/caching)
      - name: Install Macaulay2
        run: |
          sudo add-apt-repository -y ppa:macaulay2/macaulay2
          sudo apt-get update
          sudo apt-get install -y macaulay2
          M2 --version
      - name: Run the m2 oracle bucket
        env:
          QUIVERLAB_REQUIRE_M2: "1"
          NUMBA_NUM_THREADS: "2"
          OMP_NUM_THREADS: "2"
        run: .venv/bin/python -m pytest -q -m m2
```

Trigger on `push` to `dev`/`main` and on the plan branch, plus
`workflow_dispatch` (again: mirror whatever qpa.yml does — consistency
beats invention).

- [ ] **Step 2: Write the design-mining note** —
  `docs/plans/2026-08-05-plan-36-m2-design-notes.md`, three short sections
  addressed to P39/P42's authors: (1) `Complexes` ergonomics to mirror
  (`C.dd` single indexable differential object; `C_i` terms;
  `LengthLimit=>n` as first-class truncation — maps to our `top`);
  (2) the ASCII `betti` grid format (`total:` row, dots for zeros, columns
  = homological degree) as the display target for resolution tables;
  (3) the `SpectralSequences` page API (`E^r`, `E^r_{p,q}` returning the
  actual module, `spots`, `netPage` grid printing) as the SS-surface shape.
  Cite the M2 doc URLs from the 2026-08-05 research brief (they are in the
  metaplan §2 sources).

- [ ] **Step 3: README** — in the verification/oracles paragraph, extend
  "four oracle classes" wording to five, adding one clause: "and a live
  Macaulay2 bridge recomputes nc graded dimensions and commutative Ext data
  (single-vertex scope; `-m m2`)". Keep the existing badge/counts flow —
  run `... -m pytest tests/release/test_readme.py -q` after editing.

- [ ] **Step 4: Full local gate**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast`
Expected: green (the new fast guards ride along; nothing else changed
buckets). Then `... -m pytest tests/m2/ tests/release/ -q` → green.

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/m2.yml docs/plans/2026-08-05-plan-36-m2-design-notes.md README.md
git commit -m "ci(m2): required-M2 CI job; docs: fifth oracle class + P39/P42 design-mining note"
```

---

## Acceptance (Plan-36 definition of done)

1. `tests/m2/` green locally with M2 installed; skips cleanly without it;
   CI `m2.yml` green with `QUIVERLAB_REQUIRE_M2=1`.
2. `tests/release/test_markers.py` + `test_oracle_classes.py` green — the
   five-class table on the verification page matches live collection.
3. Multi-vertex and Hochschild requests refuse loudly with the honest-scope
   message; the verification page's Honest scope section says why.
4. Fast suite green; no `oracle_*` markers inside `tests/m2/`.
5. Design-mining note exists and is referenced by the metaplan's P39/P42
   cards (add the pointer to the metaplan doc in the final commit).
