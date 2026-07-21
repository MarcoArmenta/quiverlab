# quiverlab GUI on the docs landing page (Plan 10) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The docs landing page (`docs/index.md`) leads with a zero-server GUI: draw a quiver on an SVG canvas, type relations, pick a field, tick invariants, compute in the visitor's browser (Pyodide running this repo's own wheel), read MathJax-rendered results, print the worked-steps report, download TikZ/JSON, and copy a reproducing Python snippet.

**Architecture:** All execution logic lives in one Python file, `docs/gui/runner.py`, shipped verbatim into the site and imported by a Pyodide Web Worker (`docs/gui/worker.js`); it is unit-tested natively under CPython. The editor (`docs/gui/gui.js`) is dependency-free vanilla JS + SVG that emits the Plan-09 request schema (`"schema": 1`, `kind: "quiver"`). A mkdocs hook (`scripts/gui_build_hook.py`) packages the engine payload into the built site: the checkout's wheel, a manifest, and library-extracted presets.

**Tech Stack:** Python ≥ 3.10 (`quiverlab` public surface + three sanctioned trace helpers), Pyodide (pinned CDN) with its prebuilt numpy/sympy/matplotlib, vanilla JS/SVG/CSS, mkdocs-material (existing site, existing MathJax), pytest.

**Approved spec:** `docs/specs/2026-07-21-quiverlab-gui-landing-design.md`. Work happens on the existing `plan-10-gui-landing` branch.

## Global Constraints

- **Venv (always):** `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python`. System python is 3.8 — never use it.
- **Every test command** runs from the repo root as `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q ...`. No parallel processes in tests.
- **No node toolchain, no npm, no bundler.** JS is hand-written static files under `docs/gui/`. No new third-party JS beyond the pinned Pyodide CDN (same posture as the site's pinned MathJax CDN in `mkdocs.yml`).
- **Pyodide pin:** `v0.28.3` full distribution from `https://cdn.jsdelivr.net/pyodide/v0.28.3/full/`. Exact version, never a floating tag. Task 7 verifies the URL resolves at execution time; if a newer stable exists, bump the pin **in worker.js and in this line** in the same commit.
- **`runner.py` import policy:** the public `quiverlab` surface plus exactly three sanctioned internals — `quiverlab.trace.render_html.render_html`, `quiverlab.trace.provenance.references_for`, `quiverlab.trace.provenance.resolve_references` (there is no public HTML-report-as-string API). `quiverlab.engine.*` stays forbidden. The Task-1 freshness gate pins all of it.
- **No trace files, ever:** `runner.run_build` sets `quiverlab.verbose = False`; worked steps are captured via the `trace=` event-sink parameter and rendered to an HTML **string**.
- **Float-literal AST gate:** `tests/test_no_floats.py` scans `src/` only. Nothing in this plan touches `src/`; `docs/gui/runner.py` computes no algebraic values itself.
- **Degree cap:** the GUI and `runner.py` both cap HH degree ranges at `MAX_DEGREE = 10` (`0..N`, N ≤ 10).
- **CI buckets:** `tests/gui/` auto-assigns to `fast` via `tests/conftest.py::_bucket` (top dir not in `_DEEP_DIRS`). The wheel-build test carries an explicit `@pytest.mark.deep` (explicit bucket wins). Do not name any test file `test_acceptance.py`, `test_complete.py`, `test_deepen.py`, `test_properties.py`, or with prefixes `test_cs_`/`test_bardzell`/`test_minimal` — those filenames force the deep bucket.
- **Golden values in Tasks 2–4 were computed with the actual library on 2026-07-21** (loop `x`, relation `x*x*x`, GF(2): dim 3, HH⁰..² = [3, 2, 2]; Kronecker over CC: dim 4, HH⁰..³ = [1, 3, 0, 0], HH₀..₂ = [2, 0, 0], Cartan [[1, 2], [0, 1]], Coxeter LaTeX `t^{2} - 2 t + 1`, gl.dim exactly 1, center dim 1). If any golden fails, STOP and report drift — do not "fix" the expected value.
- **Commits:** conventional commits; green tests at every commit; end every commit message with:
  ```
  Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
  ```

---

### Task 1: Interface-freshness gate (+ spec import-policy amendment)

The runner consumes an exact library surface, verified live on 2026-07-21. Pin it so drift fails loudly before any GUI code runs. **If any assertion fails, STOP and report — amend the plan, do not patch around it.**

**Files:**
- Create: `tests/gui/__init__.py` (empty), `tests/gui/test_interface_freshness.py`
- Modify: `docs/specs/2026-07-21-quiverlab-gui-landing-design.md` (§7 first paragraph)

**Interfaces:**
- Consumes: the released library in the venv (`pip install -e .` already done).
- Produces: a pinned contract every later task may rely on without re-checking: `Quiver(vertices, arrows)`; `Q.algebra(relations=(), field=None, ...)`; `A.dim` (int), `A.quiver.vertices` (list), `A.quiver.arrows` (dict name → (s, t)), `A.relations` (list, `str()`-able); `A.hochschild_cohomology(top, ..., verbose=None, trace=None)` → table with `.dims`/`.kind`/`.engine`/`.references`; same for `hochschild_homology`; `A.cartan_matrix()` → list-of-lists of int; `A.coxeter_polynomial()` → sympy Poly (`.as_expr()`); `A.global_dimension()` → object with `.exact`/`.value` and informative `str()`; `A.center()` → `(int, list)`; `A.tikz()` → str; `A.citations()` → tuple of registry keys; `GF(q, modulus=None)`, `CC`; `zoo(dim_max=12, field=None)` generator of Algebras; `bibliography()` iterable of entries with `.key`/`.bibtex_key`/`.formatted`; `render_html(events, title, references)`, `references_for(events)`, `resolve_references(keys)`; error types `QuiverlabError, FieldError, RelationError, NotFiniteDimensionalError, AdmissibilityError, DepthLimitError, ExactnessError`.

- [ ] **Step 1: Write the gate test**

`tests/gui/__init__.py`: empty file.

`tests/gui/test_interface_freshness.py`:

```python
"""Plan-10 freshness gate: pin the exact quiverlab surface docs/gui/runner.py
consumes (public API + the three sanctioned trace helpers). A failure here means
library drift: STOP and amend the plan; never patch the runner around it."""
import inspect

import quiverlab as ql


def test_quiver_and_algebra_build_surface():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    assert list(inspect.signature(Q.algebra).parameters) == [
        "relations", "field", "degree_bound", "trace"]
    A = Q.algebra(relations=["x*x*x"], field=ql.GF(2))
    assert isinstance(A.dim, int) and A.dim == 3
    assert list(A.quiver.vertices) == [1]
    assert dict(A.quiver.arrows) == {"x": (1, 1)}
    assert [str(r) for r in A.relations] == ["x*x*x"]


def test_invariant_surface():
    A = ql.Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
        relations=["x*x*x"], field=ql.GF(2))
    for name in ("hochschild_cohomology", "hochschild_homology"):
        params = inspect.signature(getattr(A, name)).parameters
        assert "verbose" in params and "trace" in params, name
    ev = []
    t = A.hochschild_cohomology(2, verbose=False, trace=ev)
    assert list(t.dims) == [3, 2, 2]
    assert t.kind == "HH^" and isinstance(t.engine, str)
    assert tuple(t.references) == ("bar",)
    assert len(ev) >= 1  # engines fill the explicit event sink
    m = A.cartan_matrix()
    assert m == [[3]]
    p = A.coxeter_polynomial()
    assert hasattr(p, "as_expr")
    g = A.global_dimension()
    assert hasattr(g, "exact") and hasattr(g, "value") and str(g)
    dim_z, basis = A.center()
    assert isinstance(dim_z, int) and isinstance(basis, list)
    assert isinstance(A.tikz(), str) and A.tikz().startswith(r"\begin{tikzpicture}")
    assert isinstance(A.citations(), tuple)


def test_fields_zoo_bibliography():
    assert list(inspect.signature(ql.GF).parameters) == ["q", "modulus"]
    assert str(ql.GF(8)) == "GF(2^3)"  # q = p**n spelling the runner uses
    zoo = list(ql.zoo(dim_max=12))
    assert zoo and all(hasattr(a, "dim") and hasattr(a, "quiver") for a in zoo)
    entries = list(ql.bibliography())
    assert entries and all(
        hasattr(e, "key") and hasattr(e, "bibtex_key") and hasattr(e, "formatted")
        for e in entries)


def test_sanctioned_trace_helpers():
    from quiverlab.trace.provenance import references_for, resolve_references
    from quiverlab.trace.render_html import render_html
    A = ql.Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
        relations=["x*x"], field=ql.CC)
    ev = []
    A.hochschild_cohomology(1, verbose=False, trace=ev)
    keys = references_for(ev)
    pairs = resolve_references(keys)
    assert all(len(p) == 2 for p in pairs)
    html = render_html(list(ev), title="t", references=pairs)
    assert isinstance(html, str) and "<html" in html.lower()
    assert isinstance(render_html([], title="t", references=()), str)


def test_error_types_exist():
    for name in ("QuiverlabError", "FieldError", "RelationError",
                 "NotFiniteDimensionalError", "AdmissibilityError",
                 "DepthLimitError", "ExactnessError"):
        assert isinstance(getattr(ql, name), type), name
```

- [ ] **Step 2: Run the gate**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/test_interface_freshness.py`
Expected: **5 passed**. (This is a gate, not TDD: it must pass against the current library. Any failure = drift = STOP.)

- [ ] **Step 3: Amend spec §7 import policy**

In `docs/specs/2026-07-21-quiverlab-gui-landing-design.md`, replace the sentence:

`runner.py` uses only the **public** `quiverlab` surface (`import quiverlab`) and
the stdlib.

with:

`runner.py` uses the **public** `quiverlab` surface (`import quiverlab`), the
stdlib, and exactly three sanctioned trace helpers for the worked-steps report
string (`quiverlab.trace.render_html.render_html`,
`quiverlab.trace.provenance.references_for`,
`quiverlab.trace.provenance.resolve_references`) — there is no public
HTML-report-as-string API; the freshness gate pins these three and
`quiverlab.engine.*` stays forbidden.

Also in §7, replace the three sketched signatures (`build(algebra_spec: dict) ->
dict`, `compute(invariant: str, params: dict) -> dict`, `render_trace() -> str`)
with the final postMessage-friendly API (all JSON strings in/out):

- `run_build(request_json: str) -> str` — validates the schema, builds the
  algebra into module state, returns basic info (dimension, counts) or a typed
  error.
- `compute_one(spec: str) -> str` — runs one Plan-09 compute string
  (`"hh_cohomology:0..4"`, `"cartan"`, …) against the built algebra; worked-steps
  events accumulate in module state.
- `trace_html() -> str`, `tikz() -> str`, `python_snippet() -> str`,
  `result_bundle() -> str` — the artifacts.

- [ ] **Step 4: Commit**

```bash
git add tests/gui/ docs/specs/2026-07-21-quiverlab-gui-landing-design.md
git commit -m "test(gui): Plan-10 interface-freshness gate + spec import-policy amendment"
```

---

### Task 2: `runner.py` — schema validation and algebra build

**Files:**
- Create: `docs/gui/runner.py`, `tests/gui/conftest.py`
- Test: `tests/gui/test_runner_build.py`

**Interfaces:**
- Consumes: the Task-1 pinned surface.
- Produces (for Tasks 3–4 and worker.js): module-level `_state` dict (`algebra`, `request`, `events`, `results`); `RequestError(Exception)`; `run_build(request_json: str) -> str` returning JSON `{"ok": true, "dim": int, "n_vertices": int, "n_arrows": int, "algebra": str}` or `{"ok": false, "error": {"type": str, "message": str}[, "detail": str]}`; helper `_fail(exc) -> dict`; `_field_from_spec(spec) -> field`; constants `SCHEMA_VERSION = 1`, `MAX_DEGREE = 10`. **All public runner functions take and return JSON strings** (postMessage-friendly; no PyProxy leaks).

- [ ] **Step 1: Write the failing tests**

`tests/gui/conftest.py`:

```python
"""Load docs/gui/runner.py as a FRESH module per test (clean _state), exactly
the file the browser gets — no copies, no sys.path games."""
import importlib.util
import pathlib

import pytest

RUNNER_PATH = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"


@pytest.fixture()
def runner():
    spec = importlib.util.spec_from_file_location("gui_runner", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
```

`tests/gui/test_runner_build.py`:

```python
import json

LOOP_GF2 = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": [], "artifacts": {"pdf": False, "tikz": True},
}


def _build(runner, req):
    return json.loads(runner.run_build(json.dumps(req)))


def test_build_loop_gf2(runner):
    out = _build(runner, LOOP_GF2)
    assert out["ok"] is True
    assert out["dim"] == 3 and out["n_vertices"] == 1 and out["n_arrows"] == 1
    assert out["algebra"].startswith("Algebra of dimension 3")


def test_build_gf_prime_power(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["algebra"]["field"] = {"kind": "GF", "p": 2, "n": 2}   # GF(4)
    req["algebra"]["relations"] = ["x*x"]
    out = _build(runner, req)
    assert out["ok"] is True and out["dim"] == 2


def test_family_kind_points_at_plan_09(runner):
    req = {"schema": 1, "algebra": {"kind": "family", "family": "truncated_polynomial",
                                    "params": {"n": 2}, "field": {"kind": "CC"}},
           "compute": []}
    out = _build(runner, req)
    assert out["ok"] is False
    assert out["error"]["type"] == "RequestError"
    assert "Plan 09" in out["error"]["message"]


def test_bad_schema_version_rejected(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["schema"] = 2
    out = _build(runner, req)
    assert out["ok"] is False and out["error"]["type"] == "RequestError"


def test_field_error_surfaces_verbatim(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["algebra"]["field"] = {"kind": "GF", "p": 6, "n": 1}   # 6 is not prime
    out = _build(runner, req)
    assert out["ok"] is False
    assert out["error"]["type"] == "FieldError"
    assert out["error"]["message"]                       # library text, verbatim


def test_relation_error_surfaces_verbatim(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["algebra"]["relations"] = ["x*nosucharrow"]
    out = _build(runner, req)
    assert out["ok"] is False
    assert out["error"]["type"].endswith("Error")
    assert out["error"]["type"] != "InternalError"       # a NAMED library error


def test_malformed_arrows_rejected(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["algebra"]["arrows"] = {"x": [1]}                # not a [source, target] pair
    out = _build(runner, req)
    assert out["ok"] is False and out["error"]["type"] == "RequestError"


def test_verbose_forced_off(runner):
    import quiverlab
    quiverlab.verbose = True
    _build(runner, LOOP_GF2)
    assert quiverlab.verbose is False    # the GUI must never write trace files
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/test_runner_build.py`
Expected: FAIL — `FileNotFoundError` (no `docs/gui/runner.py`) or `AttributeError: run_build`.

- [ ] **Step 3: Write the implementation**

`docs/gui/runner.py`:

```python
"""quiverlab GUI runner (Plan 10): executes landing-page GUI requests.

Runs IDENTICALLY under CPython (pytest, tests/gui/) and Pyodide (the docs-site
Web Worker) — this exact file is shipped into the built site and imported in
the browser. Import policy: the public quiverlab surface + three sanctioned
trace helpers (render_html / references_for / resolve_references); pinned by
tests/gui/test_interface_freshness.py. quiverlab.engine.* is forbidden.

All public functions take and return JSON STRINGS (postMessage-friendly)."""
import json
import os
import traceback

os.environ.setdefault("MPLBACKEND", "Agg")   # never let matplotlib probe for a display

import quiverlab

SCHEMA_VERSION = 1
MAX_DEGREE = 10

_state = {"algebra": None, "request": None, "events": None, "results": None}


class RequestError(Exception):
    """Invalid GUI request (schema violation) — reported like a library error."""


def _fail(exc):
    if isinstance(exc, (quiverlab.QuiverlabError, RequestError)):
        return {"ok": False,
                "error": {"type": type(exc).__name__, "message": str(exc)}}
    # Unexpected: generic message to the page, full traceback for the console.
    return {"ok": False,
            "error": {"type": "InternalError",
                      "message": "unexpected engine error — details in the browser console"},
            "detail": traceback.format_exc()}


def _field_from_spec(spec):
    kind = spec.get("kind") if isinstance(spec, dict) else None
    if kind == "CC":
        return quiverlab.CC
    if kind == "GF":
        p, n = spec.get("p"), spec.get("n", 1)
        if not (isinstance(p, int) and isinstance(n, int) and p >= 2 and n >= 1):
            raise RequestError("field GF needs integers p >= 2 and n >= 1")
        return quiverlab.GF(p ** n)   # FieldError (p not prime, ...) surfaces verbatim
    raise RequestError("unknown field kind %r (expected 'CC' or 'GF')" % (kind,))


def run_build(request_json):
    """Parse + validate a schema-1 request, build the algebra, reset all state."""
    _state.update(algebra=None, request=None, events=[], results=[])
    quiverlab.verbose = False   # the GUI renders its own report; never write trace files
    try:
        req = json.loads(request_json)
        if req.get("schema") != SCHEMA_VERSION:
            raise RequestError("unsupported schema %r (this GUI speaks schema 1)"
                               % (req.get("schema"),))
        alg = req.get("algebra") or {}
        kind = alg.get("kind")
        if kind == "family":
            raise RequestError("algebra kind 'family' is the server tier (Plan 09); "
                               "this GUI submits kind 'quiver' only")
        if kind != "quiver":
            raise RequestError("unknown algebra kind %r (expected 'quiver')" % (kind,))
        vertices = alg.get("vertices")
        if not (isinstance(vertices, list) and vertices
                and all(isinstance(v, int) for v in vertices)):
            raise RequestError("algebra.vertices must be a non-empty list of integers")
        arrows = alg.get("arrows")
        if not (isinstance(arrows, dict) and all(
                isinstance(st, list) and len(st) == 2
                and all(isinstance(x, int) for x in st)
                for st in arrows.values())):
            raise RequestError("algebra.arrows must map names to [source, target] pairs")
        relations = alg.get("relations", [])
        if not (isinstance(relations, list)
                and all(isinstance(r, str) for r in relations)):
            raise RequestError("algebra.relations must be a list of strings")
        field = _field_from_spec(alg.get("field"))
        Q = quiverlab.Quiver(vertices=vertices,
                             arrows={k: (s, t) for k, (s, t) in arrows.items()})
        A = Q.algebra(relations=relations, field=field)
        _state.update(algebra=A, request=req)
        out = {"ok": True, "dim": A.dim, "n_vertices": len(vertices),
               "n_arrows": len(arrows), "algebra": repr(A).splitlines()[0]}
    except Exception as exc:
        out = _fail(exc)
    return json.dumps(out)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/test_runner_build.py`
Expected: **8 passed**.

- [ ] **Step 5: Commit**

```bash
git add docs/gui/runner.py tests/gui/conftest.py tests/gui/test_runner_build.py
git commit -m "feat(gui): runner.py — Plan-09 schema validation + algebra build (kind: quiver)"
```

---

### Task 3: `runner.py` — invariants with golden tests

**Files:**
- Modify: `docs/gui/runner.py` (append)
- Test: `tests/gui/test_runner_invariants.py`

**Interfaces:**
- Consumes: Task 2's `_state`, `run_build`, `_fail`, `RequestError`, `MAX_DEGREE`.
- Produces (for Task 4 and worker.js): `compute_one(spec: str) -> str` where `spec` is one Plan-09 compute string (`"hh_cohomology:0..4"`, `"hh_homology:0..3"`, `"cartan"`, `"coxeter_polynomial"`, `"global_dimension"`, `"center"`). Returns JSON `{"ok": true, "invariant": spec, "block": {...}}` or `{"ok": false, "invariant": spec, "error": {...}}`. Blocks: HH → `{kind, top, dims, engine, citations}`; cartan → `{matrix, latex, citations}`; coxeter → `{latex, text, citations}`; gl.dim → `{text, exact, value, citations}`; center → `{dim, basis, citations}` (basis entries are exact-value *strings* — over CC the library returns sympy MPQ rationals, which JSON cannot carry natively; amended 2026-07-21 during execution). `citations` is a list of `[bibtex_key, formatted]` pairs. Successful blocks are appended to `_state["results"]` as `{"invariant": spec, **block}`; HH events accumulate in `_state["events"]`. Also `_parse_compute(spec) -> (name, top_or_None)`.

- [ ] **Step 1: Write the failing tests**

`tests/gui/test_runner_invariants.py`:

```python
import json

LOOP_GF2 = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": [], "artifacts": {"pdf": True, "tikz": True},
}
KRONECKER_CC = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1, 2],
                "arrows": {"a": [1, 2], "b": [1, 2]},
                "relations": [], "field": {"kind": "CC"}},
    "compute": [], "artifacts": {"pdf": True, "tikz": True},
}


def _ready(runner, req):
    out = json.loads(runner.run_build(json.dumps(req)))
    assert out["ok"], out
    return out


def _one(runner, spec):
    return json.loads(runner.compute_one(spec))


def test_goldens_loop_gf2(runner):
    _ready(runner, LOOP_GF2)
    hh = _one(runner, "hh_cohomology:0..2")
    assert hh["ok"] and hh["block"]["dims"] == [3, 2, 2]
    assert hh["block"]["kind"] == "HH^" and hh["block"]["top"] == 2
    assert hh["block"]["citations"] and len(hh["block"]["citations"][0]) == 2
    assert _one(runner, "cartan")["block"]["matrix"] == [[3]]
    cox = _one(runner, "coxeter_polynomial")["block"]
    assert cox["latex"] == "t + 1"
    g = _one(runner, "global_dimension")["block"]
    assert g["text"].startswith(">=") and g["exact"] is False
    assert _one(runner, "center")["block"]["dim"] == 3


def test_goldens_kronecker_cc(runner):
    _ready(runner, KRONECKER_CC)
    assert _one(runner, "hh_cohomology:0..3")["block"]["dims"] == [1, 3, 0, 0]
    assert _one(runner, "hh_homology:0..2")["block"]["dims"] == [2, 0, 0]
    assert _one(runner, "cartan")["block"]["matrix"] == [[1, 2], [0, 1]]
    assert _one(runner, "coxeter_polynomial")["block"]["latex"] == "t^{2} - 2 t + 1"
    g = _one(runner, "global_dimension")["block"]
    assert g["exact"] is True and g["value"] == 1
    assert _one(runner, "center")["block"]["dim"] == 1


def test_cartan_latex_is_pmatrix(runner):
    _ready(runner, KRONECKER_CC)
    latex = _one(runner, "cartan")["block"]["latex"]
    assert latex == r"\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}"


def test_results_and_events_accumulate(runner):
    _ready(runner, KRONECKER_CC)
    _one(runner, "hh_cohomology:0..1")
    _one(runner, "cartan")
    assert [r["invariant"] for r in runner._state["results"]] == [
        "hh_cohomology:0..1", "cartan"]
    assert runner._state["events"]          # HH filled the shared event sink


def test_unknown_invariant_and_missing_range(runner):
    _ready(runner, LOOP_GF2)
    out = _one(runner, "determinant")
    assert not out["ok"] and out["error"]["type"] == "RequestError"
    out = _one(runner, "hh_cohomology")     # range is mandatory for HH
    assert not out["ok"] and out["error"]["type"] == "RequestError"


def test_degree_cap_enforced(runner):
    _ready(runner, LOOP_GF2)
    out = _one(runner, "hh_cohomology:0..11")
    assert not out["ok"] and "cap" in out["error"]["message"]


def test_compute_before_build_is_loud(runner):
    out = _one(runner, "cartan")
    assert not out["ok"] and out["error"]["type"] == "RequestError"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/test_runner_invariants.py`
Expected: FAIL with `AttributeError: ... 'compute_one'`.

- [ ] **Step 3: Append the implementation**

Append to `docs/gui/runner.py`:

```python
def _parse_compute(spec):
    name, _, rng = spec.partition(":")
    top = None
    if rng:
        lo, _, hi = rng.partition("..")
        if lo != "0" or not hi.isdigit():
            raise RequestError("bad compute range %r (expected 'name:0..N')" % (spec,))
        top = int(hi)
        if top > MAX_DEGREE:
            raise RequestError("degree cap is %d (got %d)" % (MAX_DEGREE, top))
    return name, top


def _citation_pairs(keys):
    from quiverlab.trace.provenance import resolve_references
    return [list(p) for p in resolve_references(tuple(keys))]


def _latex_matrix(rows):
    body = r" \\ ".join(" & ".join(str(x) for x in row) for row in rows)
    return r"\begin{pmatrix} %s \end{pmatrix}" % body


def compute_one(spec):
    """Run ONE Plan-09 compute string against the built algebra."""
    A = _state["algebra"]
    try:
        if A is None:
            raise RequestError("no algebra built (run_build first)")
        name, top = _parse_compute(spec)
        if name in ("hh_cohomology", "hh_homology"):
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            method = (A.hochschild_cohomology if name == "hh_cohomology"
                      else A.hochschild_homology)
            table = method(top, verbose=False, trace=_state["events"])
            block = {"kind": table.kind, "top": top, "dims": list(table.dims),
                     "engine": table.engine,
                     "citations": _citation_pairs(table.references)}
        elif name == "cartan":
            m = A.cartan_matrix()
            block = {"matrix": m, "latex": _latex_matrix(m),
                     "citations": _citation_pairs(A.citations())}
        elif name == "coxeter_polynomial":
            import sympy
            p = A.coxeter_polynomial()
            block = {"latex": sympy.latex(p.as_expr()), "text": str(p.as_expr()),
                     "citations": _citation_pairs(A.citations())}
        elif name == "global_dimension":
            g = A.global_dimension()
            block = {"text": str(g), "exact": g.exact, "value": g.value,
                     "citations": _citation_pairs(A.citations())}
        elif name == "center":
            dim_z, basis = A.center()
            # Basis entries are exact ints/rationals (sympy MPQ over CC) — not
            # JSON-serializable; ship them as exact strings.
            block = {"dim": dim_z,
                     "basis": [[str(x) for x in row] for row in basis],
                     "citations": _citation_pairs(A.citations())}
        else:
            raise RequestError("unknown invariant %r" % (name,))
        _state["results"].append(dict(block, invariant=spec))
        out = {"ok": True, "invariant": spec, "block": block}
    except Exception as exc:
        out = _fail(exc)
        out["invariant"] = spec
    return json.dumps(out)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/test_runner_invariants.py tests/gui/test_runner_build.py`
Expected: **15 passed** (goldens are pre-verified; a golden failure = drift = STOP).

- [ ] **Step 5: Commit**

```bash
git add docs/gui/runner.py tests/gui/test_runner_invariants.py
git commit -m "feat(gui): runner.py — invariant computation (HH ranges, Cartan, Coxeter, gl.dim, center) with goldens"
```

---

### Task 4: `runner.py` — artifacts (worked-steps HTML, TikZ, snippet, bundle)

**Files:**
- Modify: `docs/gui/runner.py` (append)
- Test: `tests/gui/test_runner_artifacts.py`

**Interfaces:**
- Consumes: Tasks 2–3 (`_state`, `_parse_compute`).
- Produces (for worker.js): `trace_html() -> str` (empty string when no events — JS disables the Print button); `tikz() -> str`; `python_snippet() -> str` (runnable reproduction); `result_bundle() -> str` (JSON: schema, request, quiverlab_version, results).

- [ ] **Step 1: Write the failing tests**

`tests/gui/test_runner_artifacts.py`:

```python
import io
import json
import contextlib

KRONECKER_CC = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1, 2],
                "arrows": {"a": [1, 2], "b": [1, 2]},
                "relations": [], "field": {"kind": "CC"}},
    "compute": ["hh_cohomology:0..2", "cartan"],
    "artifacts": {"pdf": True, "tikz": True},
}


def _run_all(runner, req):
    out = json.loads(runner.run_build(json.dumps(req)))
    assert out["ok"], out
    for spec in req["compute"]:
        assert json.loads(runner.compute_one(spec))["ok"]


def test_trace_html_after_hh(runner):
    _run_all(runner, KRONECKER_CC)
    html = runner.trace_html()
    assert "<html" in html.lower() and "References" in html


def test_trace_html_empty_without_hh(runner):
    req = dict(KRONECKER_CC, compute=["cartan"])
    _run_all(runner, req)
    assert runner.trace_html() == ""


def test_tikz(runner):
    _run_all(runner, KRONECKER_CC)
    assert runner.tikz().startswith(r"\begin{tikzpicture}")


def test_python_snippet_reproduces(runner):
    _run_all(runner, KRONECKER_CC)
    snippet = runner.python_snippet()
    assert 'Quiver(vertices=[1, 2], arrows={"a": (1, 2), "b": (1, 2)})' in snippet
    assert "A.hochschild_cohomology(2)" in snippet and "A.cartan_matrix()" in snippet
    # The GUI-to-library bridge must actually RUN.
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        exec(snippet, {"__name__": "__snippet__"})
    assert "4" in buf.getvalue()     # print(A.dim) for the Kronecker algebra


def test_result_bundle(runner):
    _run_all(runner, KRONECKER_CC)
    bundle = json.loads(runner.result_bundle())
    assert bundle["schema"] == 1
    assert bundle["request"]["algebra"]["kind"] == "quiver"
    assert bundle["quiverlab_version"]
    assert [r["invariant"] for r in bundle["results"]] == KRONECKER_CC["compute"]


def test_artifacts_before_build_are_empty(runner):
    assert runner.trace_html() == "" and runner.tikz() == ""
    assert runner.python_snippet() == ""
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/test_runner_artifacts.py`
Expected: FAIL with `AttributeError: ... 'trace_html'`.

- [ ] **Step 3: Append the implementation**

Append to `docs/gui/runner.py`:

```python
def trace_html():
    """The worked-steps report as an HTML string ('' when nothing was traced)."""
    events = _state["events"] or []
    if not events:
        return ""
    from quiverlab.trace.provenance import references_for, resolve_references
    from quiverlab.trace.render_html import render_html
    title = "Worked steps — %s" % repr(_state["algebra"]).splitlines()[0]
    return render_html(list(events), title=title,
                       references=resolve_references(references_for(events)))


def tikz():
    return "" if _state["algebra"] is None else _state["algebra"].tikz()


def python_snippet():
    """Copy-paste reproduction of the GUI computation (the GUI-to-library bridge)."""
    req = _state["request"]
    if req is None:
        return ""
    alg = req["algebra"]
    f = alg["field"]
    if f["kind"] == "CC":
        field_name, field_expr = "CC", "CC"
    else:
        q = f["p"] ** f.get("n", 1)
        field_name, field_expr = "GF", "GF(%d)" % q
    arrows = ", ".join('"%s": (%d, %d)' % (k, s, t)
                       for k, (s, t) in alg["arrows"].items())
    lines = [
        "from quiverlab import Quiver, %s" % field_name,
        "",
        "Q = Quiver(vertices=%r, arrows={%s})" % (alg["vertices"], arrows),
        "A = Q.algebra(relations=%r, field=%s)" % (list(alg.get("relations", [])),
                                                   field_expr),
        "print(A.dim)",
    ]
    calls = {"hh_cohomology": "A.hochschild_cohomology(%d)",
             "hh_homology": "A.hochschild_homology(%d)",
             "cartan": "A.cartan_matrix()", "coxeter_polynomial": "A.coxeter_polynomial()",
             "global_dimension": "A.global_dimension()", "center": "A.center()"}
    for spec in req.get("compute", []):
        name, top = _parse_compute(spec)
        call = calls[name] % top if "%d" in calls[name] else calls[name]
        lines.append("print(%s)" % call)
    return "\n".join(lines) + "\n"


def result_bundle():
    return json.dumps({"schema": SCHEMA_VERSION, "request": _state["request"],
                       "quiverlab_version": quiverlab.__version__,
                       "results": _state["results"] or []}, indent=1)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/`
Expected: **26 passed** (freshness 5 + build 8 + invariants 7 + artifacts 6).

- [ ] **Step 5: Commit**

```bash
git add docs/gui/runner.py tests/gui/test_runner_artifacts.py
git commit -m "feat(gui): runner.py — worked-steps HTML, TikZ, Python snippet, result bundle"
```

---

### Task 5: mkdocs build hook — wheel, manifest, presets

**Files:**
- Create: `scripts/gui_build_hook.py`
- Modify: `mkdocs.yml` (add `hooks:`, `extra_css:`, append to `extra_javascript:`)
- Test: `tests/gui/test_build_hook.py`

**Interfaces:**
- Consumes: the library (installed in the docs env — mkdocs-jupyter already imports it) and `pip` in the running interpreter.
- Produces (for gui.js/worker.js): in the built site, `gui/quiverlab-<version>-py3-none-any.whl`, `gui/manifest.json` = `{"schema": 1, "wheel": <filename or null>, "quiverlab_version": str}`, `gui/presets.json` = list of `{"label", "vertices", "arrows", "relations", "field"}`. Env `QLGUI_SKIP_WHEEL=1` skips only the wheel (fast `mkdocs serve` iterations → `"wheel": null`). Functions: `generate_presets() -> list[dict]`, `build_wheel(gui_dir) -> str`, `on_post_build(config)`.

- [ ] **Step 1: Write the failing tests**

`tests/gui/test_build_hook.py`:

```python
import importlib.util
import json
import pathlib

import pytest

HOOK_PATH = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "gui_build_hook.py"


def _hook():
    spec = importlib.util.spec_from_file_location("gui_build_hook", HOOK_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_presets_round_trip_through_runner(runner):
    presets = _hook().generate_presets()
    assert len(presets) >= 4
    labels = [p["label"] for p in presets]
    assert any("Kronecker" in l for l in labels)
    for p in presets:
        req = {"schema": 1,
               "algebra": {"kind": "quiver", "vertices": p["vertices"],
                           "arrows": p["arrows"], "relations": p["relations"],
                           "field": p["field"]},
               "compute": []}
        out = json.loads(runner.run_build(json.dumps(req)))
        assert out["ok"], (p["label"], out)


@pytest.mark.deep     # shells out to a full PEP-517 wheel build (~10s+)
def test_build_wheel_and_manifest(tmp_path):
    hook = _hook()
    gui = tmp_path / "site" / "gui"
    gui.mkdir(parents=True)
    name = hook.build_wheel(gui)
    assert name.startswith("quiverlab-") and name.endswith("-py3-none-any.whl")
    assert (gui / name).exists()
    hook.on_post_build({"site_dir": str(tmp_path / "site")})
    manifest = json.loads((gui / "manifest.json").read_text())
    assert manifest["schema"] == 1 and manifest["wheel"] == name
    assert manifest["quiverlab_version"]
    presets = json.loads((gui / "presets.json").read_text())
    assert isinstance(presets, list) and presets


def test_skip_wheel_env(tmp_path, monkeypatch):
    hook = _hook()
    monkeypatch.setenv("QLGUI_SKIP_WHEEL", "1")
    hook.on_post_build({"site_dir": str(tmp_path / "site")})
    manifest = json.loads((tmp_path / "site" / "gui" / "manifest.json").read_text())
    assert manifest["wheel"] is None
    assert (tmp_path / "site" / "gui" / "presets.json").exists()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/test_build_hook.py`
Expected: FAIL — `FileNotFoundError` for `scripts/gui_build_hook.py`.

- [ ] **Step 3: Write the hook**

`scripts/gui_build_hook.py`:

```python
"""mkdocs build hook (Plan 10): package the GUI engine payload into the site.

on_post_build writes, under <site>/gui/:
  quiverlab-<version>-py3-none-any.whl   pip wheel of THIS checkout (the GUI
                                         always runs the exact code being documented)
  manifest.json                          {"schema": 1, "wheel": ..., "quiverlab_version": ...}
  presets.json                           curated examples EXTRACTED via the library
                                         (A.quiver / A.relations), never hand-written

QLGUI_SKIP_WHEEL=1 skips only the (slow) wheel build for fast `mkdocs serve`
iterations; manifest.json then carries {"wheel": null} and the GUI shows an
"engine payload not built" chip instead of loading Pyodide."""
import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

REPO = pathlib.Path(__file__).resolve().parents[1]


def _preset_algebras():
    import quiverlab as ql
    ql.verbose = False
    entries = [
        ("Kronecker quiver (CC)",
         ql.Quiver(vertices=[1, 2], arrows={"a": (1, 2), "b": (1, 2)})
           .algebra(field=ql.CC),
         {"kind": "CC"}),
        ("A3 path, zero relation a*b (CC)",
         ql.Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
           .algebra(relations=["a*b"], field=ql.CC),
         {"kind": "CC"}),
        ("Truncated polynomial GF(2)[x]/(x^3)",
         ql.truncated_polynomial(n=3, field=ql.GF(2)),
         {"kind": "GF", "p": 2, "n": 1}),
    ]
    for i, A in enumerate(ql.zoo(dim_max=12)):
        if i >= 2:
            break
        entries.append(("Exact zoo #%d — dim %d (CC)" % (i + 1, A.dim), A,
                        {"kind": "CC"}))
    return entries


def generate_presets():
    return [{"label": label,
             "vertices": list(A.quiver.vertices),
             "arrows": {k: list(v) for k, v in A.quiver.arrows.items()},
             "relations": [str(r) for r in A.relations],
             "field": field_spec}
            for label, A, field_spec in _preset_algebras()]


def build_wheel(gui_dir):
    """pip-wheel this checkout into gui_dir; return the wheel filename."""
    gui_dir = pathlib.Path(gui_dir)
    # pip builds in-tree: a stale build/ from a previous run collides
    # ([Errno 17] on dist-info), so every rebuild must start clean.
    shutil.rmtree(REPO / "build", ignore_errors=True)
    with tempfile.TemporaryDirectory() as td:
        proc = subprocess.run(
            [sys.executable, "-m", "pip", "wheel", "--no-deps", "-w", td, str(REPO)],
            capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError("GUI wheel build failed:\n%s" % proc.stderr)
        wheels = list(pathlib.Path(td).glob("quiverlab-*.whl"))
        if len(wheels) != 1:
            raise RuntimeError("expected exactly one quiverlab wheel, got %r"
                               % ([w.name for w in wheels],))
        shutil.copy2(wheels[0], gui_dir / wheels[0].name)
        return wheels[0].name


def on_post_build(config):
    import quiverlab
    gui = pathlib.Path(config["site_dir"]) / "gui"
    gui.mkdir(parents=True, exist_ok=True)
    wheel = None
    if not os.environ.get("QLGUI_SKIP_WHEEL"):
        wheel = build_wheel(gui)
    (gui / "manifest.json").write_text(json.dumps(
        {"schema": 1, "wheel": wheel, "quiverlab_version": quiverlab.__version__}))
    (gui / "presets.json").write_text(json.dumps(generate_presets(), indent=1))
```

- [ ] **Step 4: Wire it into `mkdocs.yml`**

In `mkdocs.yml`, add a top-level `hooks:` block directly above `markdown_extensions:`, add `extra_css:`, and append `gui/gui.js` to `extra_javascript:`. **Amendment (2026-07-21, found live):** also add `ignore: ["gui/*.py"]` to the `mkdocs-jupyter` plugin config — its default `include` covers `*.py`, so without the ignore it renders `docs/gui/runner.py` as a notebook page and the raw file 404s for the worker:

```yaml
hooks:
  - scripts/gui_build_hook.py       # Plan 10: wheel + manifest + presets into site/gui/

extra_css:
  - gui/gui.css                     # Plan 10: landing-page GUI (no-ops elsewhere)
```

and change the `extra_javascript:` list to:

```yaml
extra_javascript:
  - javascripts/mathjax.js
  # PINNED third-party CDN (never a floating @3 tag) from a reputable source, same
  # supply-chain hygiene as stripping polyfill.io from the trace HTML. Fully-offline
  # option (Marco's call): vendor this asset under docs/javascripts/ and reference the
  # local path instead, so the published docs load no third-party JS at all.
  - https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-chtml.js
  - gui/gui.js                      # Plan 10: landing-page GUI (exits unless #qlgui exists)
```

(`docs/gui/gui.css` and `docs/gui/gui.js` are created in Task 6; mkdocs only
resolves these paths at build time, and Task 6 lands before any docs build in
this plan's acceptance. Do not run `mkdocs build` at this task.)

- [ ] **Step 5: Run tests to verify they pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/test_build_hook.py`
Expected: **3 passed** (the deep-marked wheel test runs too — it is only *bucketed* deep for CI).

- [ ] **Step 6: Commit**

```bash
git add scripts/gui_build_hook.py tests/gui/test_build_hook.py mkdocs.yml
git commit -m "build(docs): GUI engine-payload hook — checkout wheel, manifest, library-extracted presets"
```

---

### Task 6: GUI shell — SVG editor, controls, landing page embed

No pytest here (JS + markdown); the deliverable is verified with the manual checks in Step 5. Write `gui.js` exactly as shown — it contains one clearly marked `/* TASK 7 */` stub (`startWorker`) that Task 7 replaces wholesale with the real worker plumbing; everything else in the file is final.

**Files:**
- Create: `docs/gui/gui.css`, `docs/gui/gui.js`
- Modify: `docs/index.md` (full replacement shown below)

**Interfaces:**
- Consumes: `gui/presets.json`, `gui/manifest.json` (Task 5, at runtime).
- Produces (for Task 7): DOM with ids `qlgui-*`; state object `S` (`vertices: [{id,x,y}]`, `arrows: [{name,s,t}]`, `selected`, `busy`, `worker`, `engineReady`, `artifacts`); `buildRequest() -> object` (the Plan-09 schema-1 request); `setStatus(text, cls)`; `render()`; a `startWorker()` **stub** that Task 7 replaces (marked `/* TASK 7 */`). Exposed as `window.QLGUI = {S, buildRequest}` for console/manual testing.

- [ ] **Step 1: Write `docs/gui/gui.css`**

```css
/* Plan 10: landing-page GUI. Everything namespaced under #qlgui.
   Light palette by default; slate overrides at the bottom. */
#qlgui { border: 1px solid var(--md-default-fg-color--lightest); border-radius: 8px;
         padding: 12px 14px; margin: 1em 0 2em; }
#qlgui .qlgui-row { display: flex; flex-wrap: wrap; gap: 10px; align-items: center;
                    margin: 6px 0; }
#qlgui label { font-size: .72rem; }
#qlgui select, #qlgui input[type="number"], #qlgui input[type="text"] {
  font: inherit; font-size: .72rem; padding: 2px 6px;
  background: var(--md-default-bg-color); color: var(--md-default-fg-color);
  border: 1px solid var(--md-default-fg-color--lighter); border-radius: 4px; }
#qlgui input[type="number"] { width: 4.2em; }
#qlgui-relations { flex: 1 1 260px; }
#qlgui-canvas-wrap { position: relative; }
#qlgui-canvas { width: 100%; height: 340px; display: block; cursor: crosshair;
  background: var(--md-code-bg-color); border: 1px solid var(--md-default-fg-color--lightest);
  border-radius: 6px; }
#qlgui-canvas .qlv circle { fill: var(--md-default-bg-color);
  stroke: var(--md-primary-fg-color); stroke-width: 2; cursor: pointer; }
#qlgui-canvas .qlv text { font-size: 14px; fill: var(--md-default-fg-color);
  pointer-events: none; user-select: none; }
#qlgui-canvas .qla path { fill: none; stroke: var(--md-default-fg-color--light);
  stroke-width: 1.8; cursor: pointer; }
#qlgui-canvas .qla text { font-size: 13px; font-style: italic;
  fill: var(--md-default-fg-color); cursor: pointer; user-select: none; }
#qlgui-canvas .sel circle, #qlgui-canvas .sel path { stroke: #e53935; }
#qlgui-canvas .qlgui-rubber { stroke: var(--md-primary-fg-color);
  stroke-dasharray: 4 3; stroke-width: 1.5; pointer-events: none; }
#qlgui-rename { position: absolute; display: none; width: 5em; z-index: 5; }
#qlgui .qlgui-hint { font-size: .64rem; color: var(--md-default-fg-color--light);
  margin: 2px 0 8px; }
#qlgui button { font: inherit; font-size: .72rem; padding: 5px 14px; border-radius: 5px;
  border: 1px solid var(--md-primary-fg-color); cursor: pointer;
  background: var(--md-primary-fg-color); color: var(--md-primary-bg-color); }
#qlgui button.qlgui-secondary { background: transparent;
  color: var(--md-primary-fg-color); }
#qlgui button:disabled { opacity: .45; cursor: not-allowed; }
#qlgui-status { font-size: .68rem; padding: 3px 10px; border-radius: 999px;
  background: var(--md-code-bg-color); color: var(--md-default-fg-color--light); }
#qlgui-status.ok  { color: #2e7d32; }
#qlgui-status.err { color: #c62828; }
#qlgui-results .qlgui-block { border-left: 3px solid var(--md-primary-fg-color);
  padding: 4px 12px; margin: 10px 0; }
#qlgui-results .qlgui-error { border-left-color: #c62828; color: #c62828; }
#qlgui-results table { font-size: .72rem; border-collapse: collapse; }
#qlgui-results td, #qlgui-results th { border: 1px solid
  var(--md-default-fg-color--lightest); padding: 2px 10px; text-align: center; }
#qlgui-results .qlgui-cites { font-size: .62rem;
  color: var(--md-default-fg-color--light); margin-top: 4px; }
[data-md-color-scheme="slate"] #qlgui-canvas .sel circle,
[data-md-color-scheme="slate"] #qlgui-canvas .sel path { stroke: #ef5350; }
[data-md-color-scheme="slate"] #qlgui-status.ok  { color: #81c784; }
[data-md-color-scheme="slate"] #qlgui-status.err,
[data-md-color-scheme="slate"] #qlgui-results .qlgui-error { color: #ef9a9a; }
```

- [ ] **Step 2: Write `docs/gui/gui.js`**

```javascript
/* quiverlab landing-page GUI (Plan 10). Vanilla JS + SVG, zero dependencies.
   Loaded site-wide via extra_javascript: exits immediately unless #qlgui exists
   (the landing page only). Compute lives in gui/worker.js (Pyodide); this file
   is the editor, the request builder, and the results renderer. */
(function () {
  "use strict";
  var root = document.getElementById("qlgui");
  if (!root) return;

  var SVGNS = "http://www.w3.org/2000/svg";
  var R = 16;                       // vertex radius (px in SVG user units)
  var S = { vertices: [], arrows: [], nextId: 1, selected: null, dragFrom: null,
            dragMoved: false, dragOrigin: null,
            worker: null, engineReady: false, manifest: null, busy: false,
            artifacts: { tikz: "", snippet: "", bundle: "", traceHtml: "" } };

  // ---------- tiny DOM helpers ----------
  function h(tag, attrs) {
    var el = document.createElement(tag);
    for (var k in (attrs || {})) {
      if (k === "text") el.textContent = attrs[k];
      else el.setAttribute(k, attrs[k]);
    }
    for (var i = 2; i < arguments.length; i++) el.appendChild(arguments[i]);
    return el;
  }
  function sv(tag, attrs) {
    var el = document.createElementNS(SVGNS, tag);
    for (var k in (attrs || {})) el.setAttribute(k, attrs[k]);
    return el;
  }

  // ---------- static shell ----------
  root.innerHTML =
    '<div class="qlgui-row">' +
    '  <label>Preset <select id="qlgui-preset"><option value="">— build your own —</option></select></label>' +
    '  <label>Field <select id="qlgui-field"><option value="CC">CC</option><option value="GF">GF(p^n)</option></select></label>' +
    '  <label id="qlgui-p-wrap" style="display:none">p <input type="number" id="qlgui-p" value="2" min="2"></label>' +
    '  <label id="qlgui-n-wrap" style="display:none">n <input type="number" id="qlgui-n" value="1" min="1"></label>' +
    '  <button id="qlgui-clear" class="qlgui-secondary" type="button">Clear</button>' +
    '  <span id="qlgui-status">engine loads on first use</span>' +
    '</div>' +
    '<div id="qlgui-canvas-wrap">' +
    '  <svg id="qlgui-canvas" viewBox="0 0 800 340" preserveAspectRatio="xMidYMid meet"></svg>' +
    '  <input type="text" id="qlgui-rename">' +
    '</div>' +
    '<p class="qlgui-hint">Click empty space: add a vertex. Drag vertex → vertex: add an arrow ' +
    '(onto itself: a loop). Click: select. Double-click an arrow label: rename. Delete key: remove.</p>' +
    '<div class="qlgui-row"><label style="flex:1 1 260px">Relations ' +
    '<input type="text" id="qlgui-relations" placeholder="e.g. a*b - c, x*x*x"></label></div>' +
    '<div class="qlgui-row" id="qlgui-invariants">' +
    '  <label><input type="checkbox" id="qlgui-hhc" checked> HH^0..<select id="qlgui-hhc-top"></select></label>' +
    '  <label><input type="checkbox" id="qlgui-hhh"> HH_0..<select id="qlgui-hhh-top"></select></label>' +
    '  <label><input type="checkbox" id="qlgui-cartan" checked> Cartan matrix</label>' +
    '  <label><input type="checkbox" id="qlgui-coxeter_polynomial"> Coxeter polynomial</label>' +
    '  <label><input type="checkbox" id="qlgui-global_dimension"> gl.dim</label>' +
    '  <label><input type="checkbox" id="qlgui-center"> center</label>' +
    '  <label><input type="checkbox" id="qlgui-trace" checked> worked-steps report</label>' +
    '</div>' +
    '<div class="qlgui-row">' +
    '  <button id="qlgui-compute" type="button" disabled>Compute</button>' +
    '  <button id="qlgui-cancel" class="qlgui-secondary" type="button" disabled>Cancel</button>' +
    '  <button id="qlgui-print" class="qlgui-secondary" type="button" disabled>Print report (PDF)</button>' +
    '  <button id="qlgui-tikz" class="qlgui-secondary" type="button" disabled>TikZ</button>' +
    '  <button id="qlgui-json" class="qlgui-secondary" type="button" disabled>JSON</button>' +
    '  <button id="qlgui-snippet" class="qlgui-secondary" type="button" disabled>Copy Python</button>' +
    '</div>' +
    '<div id="qlgui-results"></div>';

  var el = {};
  ["preset", "field", "p-wrap", "n-wrap", "p", "n", "clear", "status", "canvas",
   "rename", "relations", "hhc", "hhc-top", "hhh", "hhh-top", "cartan",
   "coxeter_polynomial", "global_dimension", "center", "trace", "compute",
   "cancel", "print", "tikz", "json", "snippet", "results"]
    .forEach(function (id) { el[id] = document.getElementById("qlgui-" + id); });
  [el["hhc-top"], el["hhh-top"]].forEach(function (sel) {
    for (var i = 0; i <= 10; i++) sel.appendChild(h("option", { text: String(i) }));
    sel.value = "4";
  });

  function setStatus(text, cls) {
    el.status.textContent = text;
    el.status.className = cls || "";
  }

  // ---------- editor state ops ----------
  function nextArrowName() {
    var used = {};
    S.arrows.forEach(function (a) { used[a.name] = 1; });
    for (var suffix = 0; ; suffix++) {
      for (var c = 97; c <= 122; c++) {
        var n = String.fromCharCode(c) + (suffix ? String(suffix) : "");
        if (!used[n]) return n;
      }
    }
  }
  function vertexAt(id) {
    return S.vertices.filter(function (v) { return v.id === id; })[0];
  }
  function removeSelected() {
    if (!S.selected) return;
    if (S.selected.type === "vertex") {
      S.vertices = S.vertices.filter(function (v) { return v.id !== S.selected.key; });
      S.arrows = S.arrows.filter(function (a) {
        return a.s !== S.selected.key && a.t !== S.selected.key; });
    } else {
      S.arrows = S.arrows.filter(function (a) { return a.name !== S.selected.key; });
    }
    S.selected = null;
    render();
  }

  // ---------- geometry + render ----------
  function siblings(a) {   // arrows sharing the same UNORDERED vertex pair
    return S.arrows.filter(function (b) {
      return (b.s === a.s && b.t === a.t) || (b.s === a.t && b.t === a.s);
    });
  }
  function arrowPath(a) {
    var p = vertexAt(a.s), q = vertexAt(a.t);
    if (a.s === a.t) {                       // loop(s), stacked above the vertex
      var loops = S.arrows.filter(function (b) { return b.s === a.s && b.t === a.s; });
      var i = loops.indexOf(a), off = 54 + 34 * i;
      return { d: "M " + (p.x - 11) + " " + (p.y - 12) +
                  " C " + (p.x - 42) + " " + (p.y - off) + ", " +
                          (p.x + 42) + " " + (p.y - off) + ", " +
                          (p.x + 11) + " " + (p.y - 12),
               lx: p.x, ly: p.y - off + 6 };
    }
    var sib = siblings(a), i = sib.indexOf(a), n = sib.length;
    // Perpendicular from the CANONICAL (min,max) order so opposite-direction
    // arrows never collapse onto the same curve.
    var c0 = vertexAt(Math.min(a.s, a.t)), c1 = vertexAt(Math.max(a.s, a.t));
    var dx = c1.x - c0.x, dy = c1.y - c0.y, len = Math.sqrt(dx * dx + dy * dy) || 1;
    var px = -dy / len, py = dx / len;
    var k = 36 * (i - (n - 1) / 2);
    var mx = (p.x + q.x) / 2 + k * px, my = (p.y + q.y) / 2 + k * py;
    var ux = (q.x - p.x), uy = (q.y - p.y), ul = Math.sqrt(ux * ux + uy * uy) || 1;
    ux /= ul; uy /= ul;
    var sx = p.x + ux * (R + 2), sy = p.y + uy * (R + 2);
    var ex = q.x - ux * (R + 6), ey = q.y - uy * (R + 6);
    return { d: "M " + sx + " " + sy + " Q " + mx + " " + my + " " + ex + " " + ey,
             lx: 0.25 * sx + 0.5 * mx + 0.25 * ex,
             ly: 0.25 * sy + 0.5 * my + 0.25 * ey - 5 };
  }
  function render() {
    var svg = el.canvas;
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    var defs = sv("defs");
    var marker = sv("marker", { id: "qlgui-arrowhead", viewBox: "0 0 10 10",
      refX: "9", refY: "5", markerWidth: "7", markerHeight: "7", orient: "auto-start-reverse" });
    var head = sv("path", { d: "M 0 0 L 10 5 L 0 10 z" });
    head.style.fill = "currentColor";
    marker.appendChild(head); defs.appendChild(marker); svg.appendChild(defs);
    S.arrows.forEach(function (a) {
      var g = sv("g", { "class": "qla" +
        (S.selected && S.selected.type === "arrow" && S.selected.key === a.name ? " sel" : "") });
      var geo = arrowPath(a);
      var path = sv("path", { d: geo.d, "marker-end": "url(#qlgui-arrowhead)" });
      var label = sv("text", { x: geo.lx, y: geo.ly, "text-anchor": "middle" });
      label.textContent = a.name;
      [path, label].forEach(function (node) {
        node.addEventListener("click", function (e) {
          e.stopPropagation(); S.selected = { type: "arrow", key: a.name }; render();
        });
        node.addEventListener("dblclick", function (e) {
          e.stopPropagation(); startRename(a, geo);
        });
      });
      g.appendChild(path); g.appendChild(label); svg.appendChild(g);
    });
    S.vertices.forEach(function (v) {
      var g = sv("g", { "class": "qlv" +
        (S.selected && S.selected.type === "vertex" && S.selected.key === v.id ? " sel" : "") });
      var c = sv("circle", { cx: v.x, cy: v.y, r: R });
      c.addEventListener("mousedown", function (e) {
        e.preventDefault(); e.stopPropagation();
        S.dragFrom = v.id; S.dragMoved = false;
        S.dragOrigin = [e.clientX, e.clientY];
      });
      c.addEventListener("mouseup", function (e) {
        if (S.dragFrom === null) return;
        e.stopPropagation();
        // A press-release with no movement on the SAME vertex is a click
        // (select), not a self-loop; loops need a real drag gesture.
        if (S.dragFrom !== v.id || S.dragMoved) {
          S.arrows.push({ name: nextArrowName(), s: S.dragFrom, t: v.id });
          S.dragFrom = null; S.dragMoved = false; render();
        } else {
          S.dragFrom = null; S.dragMoved = false;
        }
      });
      c.addEventListener("click", function (e) {
        e.stopPropagation(); S.selected = { type: "vertex", key: v.id }; render();
      });
      var t = sv("text", { x: v.x, y: v.y + 5, "text-anchor": "middle" });
      t.textContent = String(v.id);
      g.appendChild(c); g.appendChild(t); svg.appendChild(g);
    });
    el.compute.disabled = !(S.engineReady && S.vertices.length && !S.busy);
  }

  function canvasPoint(e) {
    var pt = el.canvas.createSVGPoint();
    pt.x = e.clientX; pt.y = e.clientY;
    return pt.matrixTransform(el.canvas.getScreenCTM().inverse());
  }
  el.canvas.addEventListener("click", function (e) {
    if (e.target !== el.canvas) return;
    var p = canvasPoint(e);
    S.vertices.push({ id: S.nextId++, x: p.x, y: p.y });
    S.selected = null; render();
  });
  el.canvas.addEventListener("mousemove", function (e) {
    var old = el.canvas.querySelector(".qlgui-rubber");
    if (old) old.remove();
    if (S.dragFrom === null) return;
    if (S.dragOrigin && Math.hypot(e.clientX - S.dragOrigin[0],
                                   e.clientY - S.dragOrigin[1]) > 8) {
      S.dragMoved = true;
    }
    var v = vertexAt(S.dragFrom), p = canvasPoint(e);
    el.canvas.appendChild(sv("line", { "class": "qlgui-rubber",
      x1: v.x, y1: v.y, x2: p.x, y2: p.y }));
  });
  document.addEventListener("mouseup", function () {
    S.dragFrom = null; S.dragMoved = false; S.dragOrigin = null;
    var old = el.canvas.querySelector(".qlgui-rubber");
    if (old) old.remove();
  });
  document.addEventListener("keydown", function (e) {
    if ((e.key === "Delete" || e.key === "Backspace") &&
        document.activeElement.tagName !== "INPUT") { removeSelected(); }
  });

  // ---------- arrow rename ----------
  function startRename(arrow, geo) {
    var box = el.canvas.getBoundingClientRect();
    var scaleX = box.width / 800, scaleY = box.height / 340;
    el.rename.style.display = "block";
    el.rename.style.left = (geo.lx * scaleX - 24) + "px";
    el.rename.style.top = (geo.ly * scaleY - 12) + "px";
    el.rename.value = arrow.name;
    el.rename.focus(); el.rename.select();
    el.rename.onkeydown = function (e) {
      if (e.key === "Enter") commit();
      if (e.key === "Escape") {
        // Explicit blur: a display:none input silently keeps keyboard focus,
        // which makes the document-level Delete guard eat later deletions.
        el.rename.onblur = null;
        el.rename.style.display = "none";
        el.rename.blur();
      }
    };
    el.rename.onblur = commit;
    function commit() {
      var name = el.rename.value.trim();
      var taken = S.arrows.some(function (b) { return b !== arrow && b.name === name; });
      if (/^[A-Za-z][A-Za-z0-9_]*$/.test(name) && !taken) {
        arrow.name = name;
        el.rename.onblur = null;      // avoid blur() re-entering commit
        el.rename.style.display = "none";
        el.rename.blur();
        render();
      } else { el.rename.style.borderColor = "#c62828"; }
    }
  }

  el.clear.addEventListener("click", function () {
    S.vertices = []; S.arrows = []; S.nextId = 1; S.selected = null;
    el.relations.value = ""; el.results.innerHTML = ""; render();
  });
  el.field.addEventListener("change", function () {
    var gf = el.field.value === "GF";
    el["p-wrap"].style.display = gf ? "" : "none";
    el["n-wrap"].style.display = gf ? "" : "none";
  });

  // ---------- presets ----------
  fetch("gui/presets.json").then(function (r) { return r.ok ? r.json() : []; })
    .then(function (presets) {
      presets.forEach(function (p, i) {
        el.preset.appendChild(h("option", { value: String(i), text: p.label }));
      });
      el.preset.addEventListener("change", function () {
        if (el.preset.value === "") return;
        var p = presets[parseInt(el.preset.value, 10)];
        S.vertices = p.vertices.map(function (id, i) {
          var angle = 2 * Math.PI * i / p.vertices.length - Math.PI / 2;
          var rad = p.vertices.length === 1 ? 0 : 110;
          return { id: id, x: 400 + rad * Math.cos(angle), y: 185 + rad * Math.sin(angle) };
        });
        S.nextId = Math.max.apply(null, p.vertices.concat([0])) + 1;
        S.arrows = Object.keys(p.arrows).map(function (name) {
          return { name: name, s: p.arrows[name][0], t: p.arrows[name][1] };
        });
        el.relations.value = p.relations.join(", ");
        el.field.value = p.field.kind === "CC" ? "CC" : "GF";
        el.field.dispatchEvent(new Event("change"));
        if (p.field.kind === "GF") { el.p.value = p.field.p; el.n.value = p.field.n || 1; }
        S.selected = null; el.results.innerHTML = ""; render();
      });
    }).catch(function () { /* presets are a convenience; the editor still works */ });

  // ---------- request ----------
  function buildRequest() {
    var arrows = {};
    S.arrows.forEach(function (a) { arrows[a.name] = [a.s, a.t]; });
    var relations = el.relations.value.split(",")
      .map(function (s) { return s.trim(); }).filter(Boolean);
    var field = el.field.value === "CC" ? { kind: "CC" }
      : { kind: "GF", p: parseInt(el.p.value, 10) || 0, n: parseInt(el.n.value, 10) || 1 };
    var compute = [];
    if (el.hhc.checked) compute.push("hh_cohomology:0.." + el["hhc-top"].value);
    if (el.hhh.checked) compute.push("hh_homology:0.." + el["hhh-top"].value);
    ["cartan", "coxeter_polynomial", "global_dimension", "center"].forEach(function (k) {
      if (el[k].checked) compute.push(k);
    });
    return { schema: 1,
             algebra: { kind: "quiver",
                        vertices: S.vertices.map(function (v) { return v.id; }),
                        arrows: arrows, relations: relations, field: field },
             compute: compute,
             artifacts: { pdf: el.trace.checked, tikz: true } };
  }

  /* TASK 7 replaces this stub: worker startup, compute wiring, results
     rendering, artifact buttons. Until then the GUI is editor-only. */
  function startWorker() {
    setStatus("engine wiring lands in Task 7", "err");
  }

  window.QLGUI = { S: S, buildRequest: buildRequest };
  render();
  // Engine loads on FIRST INTENT (whole-branch review decision): pure readers
  // never pay the ~60 MB download; the first GUI touch starts it.
  var engineStarted = false;
  function ensureEngine() {
    if (engineStarted) return;
    engineStarted = true;
    startWorker();
  }
  el.canvas.addEventListener("mousedown", ensureEngine, true); // capture: circle handlers stopPropagation
  el.preset.addEventListener("change", ensureEngine);
  el.relations.addEventListener("focus", ensureEngine);
})();
```

- [ ] **Step 3: Replace `docs/index.md`**

Full new content of `docs/index.md`:

```markdown
# quiverlab

**Quivers with relations and Hochschild theory, exactly, for algebraists.**

Draw a quiver below, type relations, pick a field, and compute — right here,
with nothing to install. Everything runs **exactly** (never floating point) in
your browser, on the same engine the Python library ships.

<div id="qlgui">
  <noscript><p><strong>The interactive GUI needs JavaScript.</strong>
  The Python library works without it: <code>pip install quiverlab</code>.</p></noscript>
  <p>Loading the GUI…</p>
</div>

## Prefer code?

Finite-dimensional algebras `kQ/I` over the complex numbers (exactly — no floating
point, ever) and all finite fields: certified finite-dimensionality, Hochschild
(co)homology with cup products and Gerstenhaber brackets, the first full
Chouhy–Solotar resolution, module Ext, Cartan/Coxeter invariants, drawings, and
worked-steps traces.

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
A = Q.algebra(relations=["a*b"], field=CC)
print(A.hochschild_cohomology(3))
```

- **Tutorials** — start here (executable notebooks).
- **Under the hood** — how each object is represented and each number produced.
- **API Reference** — every public function and class.
- **Web GUI** — the form at the top of this page runs in your browser; a
  server-backed tier for big jobs is planned (Plan 09).
- **Cite** — see the JOSS paper and `CITATION.cff`.
```

- [ ] **Step 4: Sanity-run the existing gui suite (nothing broken)**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/`
Expected: **29 passed** (unchanged — this task adds no Python).

- [ ] **Step 5: Manual editor checks (fast serve, wheel skipped)**

Run: `QLGUI_SKIP_WHEEL=1 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/mkdocs serve`
(Notebook execution makes the first build take a few minutes — existing site behavior.)

In the browser at `http://127.0.0.1:8000/quiverlab/` verify:
1. GUI shell renders at the top of the landing page; status chip reads "engine wiring lands in Task 7".
2. Click empty canvas 3× → vertices 1, 2, 3 appear.
3. Drag 1→2, 1→2 again, 2→2 → two curved parallel arrows `a`, `b` and one loop `c`.
4. Double-click a label → rename to `alpha` commits; renaming another arrow to `alpha` refuses (red border).
5. Select vertex 2, press Delete → vertex and its incident arrows disappear.
6. Preset dropdown lists ≥ 4 entries; choosing "Kronecker quiver (CC)" lays out 2 vertices + 2 curved arrows and clears/loads relations & field.
7. `QLGUI.buildRequest()` in the devtools console returns `{schema: 1, algebra: {kind: "quiver", ...}}` matching the canvas.
8. Toggle the site's dark mode → canvas, chips, inputs all legible.
9. Any other docs page (e.g. Tutorials) shows no GUI and no console errors.

- [ ] **Step 6: Commit**

```bash
git add docs/gui/gui.css docs/gui/gui.js docs/index.md
git commit -m "feat(gui): SVG quiver editor + landing-page embed (schema-1 request builder)"
```

---

### Task 7: Pyodide worker — in-browser compute wired to the editor

**Files:**
- Create: `docs/gui/worker.js`
- Modify: `docs/gui/gui.js` (replace the `/* TASK 7 */` stub section)

**Interfaces:**
- Consumes: `runner.py` functions (Task 2–4 signatures), `gui/manifest.json`, the Task-6 DOM/state/`buildRequest`.
- Produces: worker protocol — main→worker `{cmd: "init", manifest}` then `{cmd: "run", request}`; worker→main `{type: "ready", version}`, `{type: "built", data}`, `{type: "result", data}`, `{type: "trace", html}`, `{type: "artifacts", tikz, snippet, bundle}`, `{type: "done"}`, `{type: "fatal", message}`. Cancel = `worker.terminate()` + fresh `startWorker()`.

- [ ] **Step 1: Verify the Pyodide pin resolves**

Run: `curl -sI https://cdn.jsdelivr.net/pyodide/v0.28.3/full/pyodide.js | head -1`
Expected: `HTTP/2 200`. If 404: pick the latest stable `v0.28.x` from https://github.com/pyodide/pyodide/releases, and use it consistently in Step 2 **and** update the Global Constraints pin line of this plan in the same commit.

- [ ] **Step 2: Write `docs/gui/worker.js`**

```javascript
/* quiverlab GUI Web Worker (Plan 10): Pyodide + the site's own wheel.
   PINNED Pyodide version — never a floating tag (supply-chain posture of
   mkdocs.yml). Protocol: see gui.js. All runner I/O is JSON strings. */
"use strict";
var PYODIDE_BASE = "https://cdn.jsdelivr.net/pyodide/v0.28.3/full/";
var pyodide = null, runner = null;

self.onmessage = function (e) {
  var m = e.data;
  var job = m.cmd === "init" ? init(m.manifest)
          : m.cmd === "run" ? run(m.request)
          : Promise.reject(new Error("unknown cmd " + m.cmd));
  job.catch(function (err) {
    self.postMessage({ type: "fatal", message: String(err && err.message || err) });
  });
};

async function init(manifest) {
  importScripts(PYODIDE_BASE + "pyodide.js");
  pyodide = await loadPyodide({ indexURL: PYODIDE_BASE });
  await pyodide.loadPackage(["micropip", "numpy", "sympy", "matplotlib"]);
  var micropip = pyodide.pyimport("micropip");
  await micropip.install(new URL(manifest.wheel, self.location.href).href);
  var resp = await fetch(new URL("runner.py", self.location.href));
  pyodide.FS.writeFile("runner.py", await resp.text());
  runner = pyodide.pyimport("runner");
  self.postMessage({ type: "ready", version: manifest.quiverlab_version });
}

async function run(request) {
  // Pass the FULL request: python_snippet/result_bundle read request.compute
  // from the stored request (run_build itself only validates schema+algebra).
  var built = JSON.parse(runner.run_build(JSON.stringify(request)));
  self.postMessage({ type: "built", data: built });
  if (built.ok) {
    for (var i = 0; i < request.compute.length; i++) {
      var res = JSON.parse(runner.compute_one(request.compute[i]));
      self.postMessage({ type: "result", data: res });
    }
    self.postMessage({ type: "trace",
                       html: request.artifacts && request.artifacts.pdf
                             ? runner.trace_html() : "" });
    self.postMessage({ type: "artifacts", tikz: runner.tikz(),
                       snippet: runner.python_snippet(),
                       bundle: runner.result_bundle() });
  }
  self.postMessage({ type: "done" });
}
```

- [ ] **Step 3: Replace the Task-7 stub in `docs/gui/gui.js`**

Replace the block from `/* TASK 7 replaces this stub ... */` through the end of the `startWorker` stub (keep everything after it) with:

```javascript
  // ---------- engine (Pyodide worker) ----------
  function startWorker() {
    if (S.worker) S.worker.terminate();
    S.engineReady = false;
    S.worker = new Worker("gui/worker.js");
    S.worker.onmessage = onWorkerMessage;
    S.worker.onerror = function (e) {
      setStatus("engine error — see console", "err");
      console.error("qlgui worker:", e);
    };
    if (S.manifest) {                    // restart path (after Cancel)
      setStatus("engine reloading…");
      S.worker.postMessage({ cmd: "init", manifest: S.manifest });
      return;
    }
    fetch("gui/manifest.json")
      .then(function (r) { if (!r.ok) throw new Error("no manifest"); return r.json(); })
      .then(function (m) {
        if (!m.wheel) { setStatus("engine payload not built (QLGUI_SKIP_WHEEL)", "err"); return; }
        S.manifest = m;
        setStatus("engine loading… (~60 MB once, then cached)");
        S.worker.postMessage({ cmd: "init", manifest: m });
      })
      .catch(function () { setStatus("engine manifest missing — editor-only preview", "err"); });
  }

  function setBusy(b) {
    S.busy = b;
    el.cancel.disabled = !b;
    render();
  }

  function onWorkerMessage(e) {
    var m = e.data;
    if (m.type === "ready") {
      S.engineReady = true;
      setStatus("engine ready — quiverlab " + m.version, "ok");
      render();
    } else if (m.type === "built") {
      if (m.data.ok) {
        el.results.appendChild(h("div", { "class": "qlgui-block",
          text: m.data.algebra + "  (dim = " + m.data.dim + ")" }));
      } else { renderError(m.data); }
    } else if (m.type === "result") {
      if (m.data.ok) renderBlock(m.data); else renderError(m.data);
    } else if (m.type === "trace") {
      S.artifacts.traceHtml = m.html;
      el.print.disabled = !m.html;
    } else if (m.type === "artifacts") {
      S.artifacts.tikz = m.tikz; S.artifacts.snippet = m.snippet;
      S.artifacts.bundle = m.bundle;
      el.tikz.disabled = el.json.disabled = el.snippet.disabled = false;
    } else if (m.type === "done") {
      setBusy(false);
      setStatus("engine ready — quiverlab " + (S.manifest ? S.manifest.quiverlab_version : ""), "ok");
    } else if (m.type === "fatal") {
      console.error("qlgui engine:", m.message);
      setStatus("engine failed — see console", "err");
      setBusy(false);
    }
  }

  function renderError(res) {
    if (res.detail) console.error(res.detail);
    var div = h("div", { "class": "qlgui-block qlgui-error",
      text: res.error.type + ": " + res.error.message });
    el.results.appendChild(div);
  }

  function citesLine(block) {
    return h("div", { "class": "qlgui-cites",
      text: (block.citations || []).map(function (c) { return c[1]; }).join(" · ") });
  }

  function renderBlock(res) {
    var b = res.block, name = res.invariant.split(":")[0];
    var div = h("div", { "class": "qlgui-block" });
    if (name === "hh_cohomology" || name === "hh_homology") {
      var sup = name === "hh_cohomology";
      var head = h("tr"), row = h("tr");
      head.appendChild(h("th", { text: "n" }));
      row.appendChild(h("th", { text: sup ? "dim HH^n" : "dim HH_n" }));
      b.dims.forEach(function (d, n) {
        head.appendChild(h("td", { text: String(n) }));
        row.appendChild(h("td", { text: String(d) }));
      });
      div.appendChild(h("p", { text: sup ? "Hochschild cohomology" : "Hochschild homology" }));
      div.appendChild(h("table", {}, head, row));
      div.appendChild(h("div", { "class": "qlgui-cites", text: b.engine }));
    } else if (name === "cartan") {
      div.appendChild(h("p", { text: "Cartan matrix:" }));
      // className EXACTLY "arithmatex": the site's MathJax config matches
      // class patterns against the full className string, so a combined
      // "qlgui-block arithmatex" is silently skipped (found live).
      div.appendChild(h("p", { "class": "arithmatex", text: "\\[ C = " + b.latex + " \\]" }));
    } else if (name === "coxeter_polynomial") {
      div.appendChild(h("p", { "class": "arithmatex", text: "\\[ \\chi(t) = " + b.latex + " \\]" }));
    } else if (name === "global_dimension") {
      div.appendChild(h("p", { text: b.text }));
    } else if (name === "center") {
      div.appendChild(h("p", { "class": "arithmatex", text: "\\( \\dim Z(A) = " + b.dim + " \\)" }));
    }
    div.appendChild(citesLine(b));
    el.results.appendChild(div);
    if (window.MathJax && window.MathJax.typesetPromise) {
      // Full-page sweep, NOT typesetPromise([div]): with explicit roots the
      // walker never consults the root's own class, so the site's
      // ignoreHtmlClass ".*|" config silently skips the block (found live).
      window.MathJax.typesetPromise();
    }
  }

  // ---------- buttons ----------
  el.compute.addEventListener("click", function () {
    if (S.busy || !S.engineReady) return;
    el.results.innerHTML = "";
    S.artifacts = { tikz: "", snippet: "", bundle: "", traceHtml: "" };
    el.print.disabled = el.tikz.disabled = el.json.disabled = el.snippet.disabled = true;
    setBusy(true);
    setStatus("computing…");
    S.worker.postMessage({ cmd: "run", request: buildRequest() });
  });
  el.cancel.addEventListener("click", function () {
    if (!S.busy) return;
    setBusy(false);
    startWorker();                             // sets its own transient status…
    setStatus("cancelled — engine reloading…"); // …so the acknowledgment must land last
  });
  function download(name, text, type) {
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], { type: type }));
    a.download = name;
    a.click();
    URL.revokeObjectURL(a.href);
  }
  el.print.addEventListener("click", function () {
    var url = URL.createObjectURL(new Blob([S.artifacts.traceHtml], { type: "text/html" }));
    var w = window.open(url, "_blank");
    if (w) w.addEventListener("load", function () { w.print(); });
  });
  el.tikz.addEventListener("click", function () {
    download("quiver.tex", S.artifacts.tikz, "text/x-tex");
  });
  el.json.addEventListener("click", function () {
    download("quiverlab-result.json", S.artifacts.bundle, "application/json");
  });
  el.snippet.addEventListener("click", function () {
    navigator.clipboard.writeText(S.artifacts.snippet).then(function () {
      setStatus("Python snippet copied", "ok");
    });
  });
```

(The trailing `window.QLGUI = ...`, `render()`, and `requestIdleCallback` lines from Task 6 stay unchanged below this block.)

- [ ] **Step 4: Manual end-to-end checks (full build, wheel included)**

Run: `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/mkdocs serve` (no skip env — the hook builds the wheel; needs network for the Pyodide CDN).

At `http://127.0.0.1:8000/quiverlab/` verify:
1. Status chip starts "engine loads on first use"; after the FIRST GUI interaction (canvas touch / preset pick / relations focus) it goes "engine loading…" → "engine ready — quiverlab 0.1.0.dev0" (first load downloads Pyodide; ~1–2 min on ordinary broadband, instant when cached).
2. Preset "Kronecker quiver (CC)", Compute → algebra header "(dim = 4)"; with the default HH^ range 0..4 the table reads `1, 3, 0, 0, 0`; Cartan `pmatrix` typeset by MathJax.
3. Draw a loop `x` on one vertex, relation `x*x*x`, field GF p=2 n=1, all boxes ticked, Compute → results appear one at a time; gl.dim shows the `>= 32` certified-bound text.
4. Relation typo (`x*noarrow`) → red error block naming the library error type, verbatim message; no crash; next Compute works.
5. GF with p=6 → verbatim `FieldError` block.
6. Long run (zoo preset, HH^0..10) → Cancel mid-run → "cancelled — engine reloading…" → "engine ready" again; Compute works after.
7. "Print report (PDF)" opens the worked-steps HTML in a new tab and the print dialog appears (Save as PDF works).
8. TikZ downloads `quiver.tex` starting `\begin{tikzpicture}`; JSON downloads the bundle; "Copy Python" then pasting into a REPL reproduces the computation.
9. Devtools console: no uncaught errors through all of the above.

- [ ] **Step 5: Commit**

```bash
git add docs/gui/worker.js docs/gui/gui.js
git commit -m "feat(gui): Pyodide worker — in-browser compute, cancel/restart, MathJax results, artifacts"
```

---

### Task 8: Acceptance — full suites, strict docs build, changelog, Plan-09 amendment

**Files:**
- Modify: `CHANGELOG.md`, `README.md`, `docs/plans/2026-07-18-plan-09-web.md` (amendment note at top)

**Interfaces:**
- Consumes: everything above.
- Produces: a merge-ready branch.

- [ ] **Step 1: Full test suites**

Run (buckets must stay a partition; the markers test shells out, so run it explicitly too):
```bash
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q -m fast
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q tests/gui/ tests/release/test_markers.py
```
Expected: all green; `test_buckets_partition_the_suite` passes with the new `tests/gui/` files.

- [ ] **Step 2: Strict docs build with the GUI payload**

Run: `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/mkdocs build --strict`
Expected: exit 0. Then verify the payload:
```bash
ls site/gui/
```
Expected: `gui.css  gui.js  manifest.json  presets.json  quiverlab-0.1.0.dev0-py3-none-any.whl  runner.py  worker.js` (plus nothing else).

- [ ] **Step 3: CHANGELOG entry**

In `CHANGELOG.md`, under the unreleased/top section (match the file's existing heading style), add:

```markdown
- **In-browser GUI on the docs landing page (Plan 10).** Draw a quiver, type
  relations, pick CC/GF(p^n), and compute HH ranges, Cartan/Coxeter, gl.dim and
  the center — no installation; the page runs the repo's own wheel under
  Pyodide in a Web Worker. Worked-steps report printable to PDF; TikZ/JSON
  downloads; copy-paste Python reproduction. GUI requests use the Plan-09
  schema (`kind: "quiver"`), and `docs/gui/runner.py` is the execution
  semantics a future server tier reuses.
```

- [ ] **Step 4: README pointer**

In `README.md`, in the feature-tour/bullets section (match surrounding style), add one bullet:

```markdown
- **In-browser GUI** — the [docs landing page](https://marcoarmenta.github.io/quiverlab/)
  computes examples with nothing installed (Pyodide running the same exact engine).
```

- [ ] **Step 5: Plan-09 amendment note**

At the top of `docs/plans/2026-07-18-plan-09-web.md`, directly under the H1, insert:

```markdown
> **Amendment (2026-07-21, Plan 10):** the docs landing page now ships an
> in-browser GUI (`docs/gui/`) that emits this plan's request schema with
> `algebra.kind == "quiver"` — the v2 canvas editor exists. When this server
> tier is built: (a) the worker should adopt `docs/gui/runner.py`'s execution
> semantics (schema validation, per-invariant blocks, citation pairs, trace
> HTML) rather than reimplementing them; (b) v1 must accept `kind: "quiver"`
> in addition to `kind: "family"`.
```

- [ ] **Step 6: Full manual browser checklist (final pass)**

Repeat Task 7 Step 4's checks 1–9 once on the `mkdocs build` output served statically:
```bash
cd site && python3 -m http.server 8000
```
Browse `http://127.0.0.1:8000/` (site root here — no `/quiverlab/` prefix needed for this server) and run the checklist.

- [ ] **Step 7: Commit**

```bash
git add CHANGELOG.md README.md docs/plans/2026-07-18-plan-09-web.md
git commit -m "docs(gui): Plan-10 acceptance — changelog, README pointer, Plan-09 amendment"
```
