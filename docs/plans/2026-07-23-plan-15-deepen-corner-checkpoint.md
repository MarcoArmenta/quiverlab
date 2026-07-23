# Plan 15 — Corner-mode checkpoint format for `deepen` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `engine/deepen.py` (the checkpointed, resumable, incremental-N driver) works
for multi-vertex algebras — lifting the Plan-13 `NotImplementedError` boundary and
unlocking cluster-scale multi-vertex scans.

**Architecture:** Corner data (`_CornerContext`, `gens0`, `rad_ab_pairs`, the engine
itself) is deterministic from `(A, prime)`, so the checkpoint persists only what the
stepper accumulates: `cur`/`cur_r`/`rks` (as today) **plus the per-degree corner
`tags`** (tiny — one vertex pair per generator). Resume rebuilds the full state via
`_init_resolution(A, prime)` and overlays the persisted fields — this also removes the
current resume path's duplicate `AeEngine` build. The deepen loop's per-degree HH
finalization grows a corner branch mirroring `minimal_homology_dims`:
`_corner_contracted_degree` (needs tags at degrees k−2, k−1, k) and
`dimn = Σ corner_dim_A(tags[k−1])` instead of `_contracted_degree` and `m·r_{k−1}`.
A mode-mismatch guard (`QuiverlabError`) refuses a corner algebra resumed against a
local checkpoint dir and vice versa (detected via presence of `"tags"` in the
payload — old local checkpoints stay loadable).

**Tech Stack:** pure Python + numpy int64 mod p (the corner path is pure Python by
design; the local free path keeps its kernel acceleration bit-for-bit). Tests live in
`tests/engine/` → auto-bucketed **deep**.

## Global Constraints

- No float literals anywhere in `src/` (AST gate `tests/test_no_floats.py`).
- Python is always `.venv/bin/python`; tests run as
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q …`.
- Engines must agree exactly on numba and pure paths (`QUIVERLAB_NO_NUMBA=1`).
- Oracles live, never hardcoded dims: `minimal_homology_dims` (Plan-13-certified).
- Local (one-vertex) deepen behavior must be bit-for-bit unchanged (existing
  `tests/engine/test_deepen.py` is the regression gate).
- Conventional commits; green tests at every commit; merge/push only when Marco asks.

---

### Task 1: Corner-mode deepen — fresh runs, HH finalization, checkpoint payload

**Files:**
- Create: `tests/engine/test_deepen_corner.py`
- Modify: `src/quiverlab/engine/deepen.py`

**Interfaces:**
- Consumes: `_init_resolution` / `_advance_resolution` /
  `_corner_contracted_degree(eng, ctx, gens_n, tags_n, tags_nm1, p)` /
  `_CornerContext.corner_dim_A(tag)` from `engine/resolutions_minimal.py`;
  `QuiverlabError` from `quiverlab.errors`.
- Produces: `deepen(A, ckpt_dir, …)` accepting multi-vertex `A`; corner checkpoint
  payload = local payload + `"tags"` key (dict degree → list of `(i, j)` vertex
  pairs). Task 2 relies on that payload shape for resume.

- [ ] **Step 1: Write the failing tests** — `tests/engine/test_deepen_corner.py`:

```python
"""Plan 15: corner-mode checkpoints for deepen (multi-vertex resumable driver).

Before this plan `deepen` refused multi-vertex algebras (NotImplementedError, the
Plan-13 boundary): its checkpoint payload and manual resume-state rebuild had no
corner data.  Corner data (_CornerContext, gens0, rad_ab_pairs, the engine) is
deterministic from (A, prime), so the checkpoint persists only what the stepper
accumulates: cur / cur_r / rks / tags (+ the rolling last_gens and the HH /
per_degree records).

Oracle: minimal_homology_dims (the Plan-13-certified corner engine) — live,
never hardcoded."""
import pytest

import quiverlab as ql
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.deepen import deepen, _load_ckpt
from quiverlab.engine.hh_engine import truncated_polynomial
from quiverlab.engine.resolutions_minimal import minimal_homology_dims
from quiverlab.errors import QuiverlabError

PRIMES = (32003, 2, 3, 5)


def _eng(vertices, arrows, relations, p=32003):
    Q = ql.Quiver(vertices, arrows)
    return to_engine(Q.algebra(relations=relations, field=ql.GF(p)))


def _cn32(p=32003):
    """kZ_3/rad^2: periodic, HH nonzero in every degree (the resume fixture)."""
    return _eng([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)},
                ["a*b", "b*c", "c*a"], p=p)


def test_deepen_corner_matches_minimal_homology_dims(tmp_path):
    """CN(3,2) over four primes: deepen's HH_* == the batch corner engine."""
    for p in PRIMES:
        A = _cn32(p=p)
        out = deepen(A, str(tmp_path / ("ck%d" % p)), prime=p, max_degree=4)
        ref = minimal_homology_dims(A, 4, primes=(p,))[p]
        assert out["HH"] == ref, "p=%d: %s != %s" % (p, out["HH"], ref)


def test_deepen_corner_termination_ka2(tmp_path):
    """kA_2 (hereditary): the corner resolution terminates; deepen must report it
    (stop_reason='terminated', hochschild_dim=1) with HH == the batch engine."""
    A = _eng([1, 2], {"a": (1, 2)}, [])
    out = deepen(A, str(tmp_path / "ck"), prime=32003)
    assert out["stop_reason"] == "terminated" and out["terminated"]
    assert out["hochschild_dim"] == 1
    ref = minimal_homology_dims(A, 1, primes=(32003,))[32003]
    assert out["HH"] == ref


def test_deepen_corner_nonmonomial_square(tmp_path):
    """kQ/(ab - cd) (dim 9, non-monomial multi-vertex): HH_0 = 4 (Plan-13 pin)."""
    A = _eng([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)},
             ["a*b - c*d"])
    out = deepen(A, str(tmp_path / "ck"), prime=32003, max_degree=2)
    ref = minimal_homology_dims(A, 2, primes=(32003,))[32003]
    assert out["HH"] == ref[:len(out["HH"])]
    assert out["HH"][0] == 4


def test_deepen_corner_checkpoint_has_tags(tmp_path):
    """The corner checkpoint persists per-degree tags, consistent with rks."""
    A = _cn32()
    ck = str(tmp_path / "ck")
    deepen(A, ck, prime=32003, max_degree=3)
    payload = _load_ckpt(ck)
    assert payload is not None and "tags" in payload
    for n_, tgs in payload["tags"].items():
        assert len(tgs) == payload["rks"][n_]
        assert all(len(tg) == 2 for tg in tgs)


def test_deepen_corner_memory_wall(tmp_path):
    """A 1-byte transient budget stops with stop_reason='memory' at degree 1 and
    the same (empty) exact HH prefix as the batch engine under the same budget."""
    A = _cn32()
    out = deepen(A, str(tmp_path / "ck"), prime=32003, max_transient_bytes=1)
    assert out["stop_reason"] == "memory"
    assert out["wall_radK_bytes"] > 1 and out["wall_degree"] == 1
    ref = minimal_homology_dims(A, 6, primes=(32003,), max_transient_bytes=1)[32003]
    assert out["HH"] == ref


def test_deepen_corner_finalize_only(tmp_path):
    """finalize_only re-emits the summary from the corner checkpoint unchanged."""
    A = _cn32()
    ck = str(tmp_path / "ck")
    full = deepen(A, ck, prime=32003, max_degree=3)
    fin = deepen(A, ck, prime=32003, finalize_only=True)
    assert fin["stop_reason"] == "checkpoint"
    assert fin["HH"] == full["HH"]
    assert fin["max_degree_reached"] == full["max_degree_reached"]
```

(The resume and mode-mismatch tests are Task 2 — they need the resume path.)

- [ ] **Step 2: Run the new tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_deepen_corner.py -p no:cacheprovider`
Expected: FAIL — every test raises `NotImplementedError: deepen supports local algebras only …`

- [ ] **Step 3: Implement corner mode in `deepen.py`**

(a) Imports — replace the `resolutions_minimal` import block and drop the now-unused
`AeEngine`; add `QuiverlabError`:

```python
from quiverlab.engine.resolutions_minimal import (
    _init_resolution, _advance_resolution, _contracted_degree,
    _corner_contracted_degree)
from quiverlab.engine.hh_engine import rank_mod_p
from quiverlab.errors import QuiverlabError
```

(b) Replace the fresh-start/resume block (both `NotImplementedError` guards and the
manual resume-state rebuild) with a single `_init_resolution` + overlay:

```python
    st = _init_resolution(A, prime)      # engine, rad pairs (+ corner ctx, gens0):
    corner = st.get("corner") is not None   # deterministic from (A, prime) -- rebuilt, never pickled
    if ck is None:
        last_gens = None                                    # cols[0] is None
        HH = []
        per_degree = []
        log("deepen: fresh start (dim=%d, prime=%d%s)"
            % (A.m, prime, ", corner mode" if corner else ""))
    else:
        if corner != ("tags" in ck):
            raise QuiverlabError(
                "deepen: ckpt_dir holds a %s-mode checkpoint but this algebra needs "
                "%s mode" % ("corner" if "tags" in ck else "local",
                             "corner" if corner else "local"),
                hint="each ckpt_dir belongs to one (algebra, prime) run -- point "
                     "this run at its own directory")
        st.update({"cur": ck["cur"], "cur_r": ck["cur_r"], "rks": ck["rks"],
                   "n": ck["n"], "cols": {}})
        if corner:
            st["tags"] = ck["tags"]
        last_gens = ck["last_gens"]
        HH = list(ck["HH"])
        per_degree = list(ck["per_degree"])
        log("deepen: resumed at degree n=%d (HH so far: %s)" % (ck["n"], HH))
```

(c) HH finalization — corner branch mirroring `minimal_homology_dims`'s corner path
(the contracted complex has corner blocks `e_w A e_v`, so both the differentials and
the term dimension change):

```python
        if k >= 1:
            r_km1 = st["rks"].get(k - 1, 0)
            if corner:
                ctx = st["corner"]
                tags = st["tags"]
                dbar_km1 = (_corner_contracted_degree(
                                st["eng"], ctx, last_gens or [], tags.get(k - 1, []),
                                tags.get(k - 2, []), prime)
                            if (k - 1) >= 1 and r_km1 > 0 else None)
                dbar_k = (_corner_contracted_degree(
                              st["eng"], ctx, gens or [], tags.get(k, []),
                              tags.get(k - 1, []), prime)
                          if st["rks"].get(k, 0) > 0 else None)
                dimn = sum(ctx.corner_dim_A(tg) for tg in tags.get(k - 1, []))
            else:
                r_km2 = st["rks"].get(k - 2, 0)
                dbar_km1 = (_contracted_degree(st["eng"], last_gens or [], r_km2, k - 1)
                            if (k - 1) >= 1 and r_km1 > 0 else None)
                dbar_k = (_contracted_degree(st["eng"], gens or [], r_km1, k)
                          if st["rks"].get(k, 0) > 0 else None)
                dimn = m * r_km1
            rn = rank_mod_p(dbar_km1, prime) if dbar_km1 is not None else 0
            rnp1 = rank_mod_p(dbar_k, prime) if dbar_k is not None else 0
            HH.append(int(dimn - rn - rnp1))                # HH_{k-1}
```

(d) Checkpoint payload — persist `tags` in corner mode (everything else unchanged):

```python
        last_gens = gens                                    # roll forward
        payload = {"n": st["n"], "cur": st["cur"], "cur_r": st["cur_r"],
                   "rks": st["rks"], "last_gens": last_gens,
                   "HH": HH, "per_degree": per_degree}
        if corner:
            payload["tags"] = st["tags"]        # tiny: one vertex pair per generator
        _save_ckpt(ckpt_dir, payload)
```

(e) Docstrings — module docstring and `deepen()` docstring: state that both local
(free path) and multi-vertex (Plan-13 corner path) algebras are supported, and that
corner checkpoints persist only the extra per-degree `tags` (corner data is
deterministic from `(A, prime)`).

- [ ] **Step 4: Run the new tests to verify they pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_deepen_corner.py -p no:cacheprovider`
Expected: PASS (7 tests; the two Task-2 tests are not written yet)

- [ ] **Step 5: Run the local regression gate (both kernel paths)**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_deepen.py -p no:cacheprovider`
Run: `QUIVERLAB_NO_NUMBA=1 .venv/bin/python -m pytest -q tests/engine/test_deepen.py tests/engine/test_deepen_corner.py -p no:cacheprovider`
Expected: PASS — local deepen bit-for-bit unchanged; corner path identical pure/numba.

### Task 2: Resume + mode-mismatch guard

**Files:**
- Modify: `tests/engine/test_deepen_corner.py` (append two tests)

**Interfaces:**
- Consumes: Task 1's payload shape (`"tags"` key present iff corner mode) and the
  `QuiverlabError` mismatch guard already implemented in Task 1 step 3(b).
- Produces: certified resume semantics for corner checkpoints.

- [ ] **Step 1: Write the failing/blank tests** — append to `tests/engine/test_deepen_corner.py`:

```python
def test_deepen_corner_resume(tmp_path):
    """A stop-early run then a continue run == one full fresh run == the oracle."""
    A = _cn32()
    ck = str(tmp_path / "ck")
    deepen(A, ck, prime=32003, max_degree=2)                # stop early
    out = deepen(A, ck, prime=32003, max_degree=5)          # resume -> continue
    ref = minimal_homology_dims(A, 5, primes=(32003,))[32003]
    assert out["HH"] == ref
    assert out["max_degree_reached"] >= 3                   # actually went past the cap
    fresh = deepen(A, str(tmp_path / "ck2"), prime=32003, max_degree=5)
    assert out["HH"] == fresh["HH"]


def test_deepen_mode_mismatch_refuses(tmp_path):
    """A corner algebra resumed against a local checkpoint dir (and vice versa)
    must refuse loudly -- a silent overlay would corrupt the resolution state."""
    local, corner = truncated_polynomial(3), _cn32()
    ck_local = str(tmp_path / "ck_local")
    deepen(local, ck_local, prime=32003, max_degree=2)
    with pytest.raises(QuiverlabError):
        deepen(corner, ck_local, prime=32003, max_degree=3)
    ck_corner = str(tmp_path / "ck_corner")
    deepen(corner, ck_corner, prime=32003, max_degree=2)
    with pytest.raises(QuiverlabError):
        deepen(local, ck_corner, prime=32003, max_degree=3)
```

- [ ] **Step 2: Run them**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_deepen_corner.py -p no:cacheprovider`
Expected: PASS if Task 1's implementation is complete (the resume overlay and guard
are Task 1 step 3(b)); any FAIL here is a real resume bug — fix before proceeding.

- [ ] **Step 3: Commit (code + tests)**

```bash
git add src/quiverlab/engine/deepen.py tests/engine/test_deepen_corner.py
git commit -m "feat(engine): deepen supports multi-vertex algebras -- corner-mode checkpoints (Plan 15)"
```

### Task 3: Docs + status updates

**Files:**
- Modify: `CLAUDE.md` (status paragraph: drop the "deepen stays local-only" boundary,
  add Plan 15 to the delivered list)
- Modify: `docs/plans/ROADMAP.md` (add DELIVERED row 15)
- Modify: `docs/internals/05-resolutions.md` (update the deepen local-only mention)
- Modify: `docs/plans/DEEPER-ENGINES-BACKLOG.md` (already ticked at branch start —
  verify wording: checkbox + plan number + date)

- [ ] **Step 1: Make the doc edits** (each is a one-paragraph factual update; the
  deepen boundary sentence in CLAUDE.md becomes: Plan 15 delivered — corner-mode
  checkpoints persist per-degree `tags` only, everything else rebuilt from
  `(A, prime)`; mode-mismatch guard refuses cross-mode ckpt dirs)

- [ ] **Step 2: Full deep + fast suites**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m deep -p no:cacheprovider`
Expected: PASS (~19 min)
Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast -p no:cacheprovider`
Expected: PASS (includes the `test_no_floats.py` AST gate over the edited `src/` file)

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md docs/plans/ROADMAP.md docs/internals/05-resolutions.md docs/plans/DEEPER-ENGINES-BACKLOG.md
git commit -m "docs: Plan-15 status -- deepen corner-mode checkpoints delivered"
```

## Validation matrix (tests/engine/test_deepen_corner.py)

1. CN(3,2) `kZ_3/rad²` over 4 primes: deepen HH ≡ `minimal_homology_dims` (periodic,
   HH nonzero every degree — the strongest dim cross-check).
2. `kA_2`: finite corner resolution → `stop_reason='terminated'`, `hochschild_dim=1`.
3. Commutative square `kQ/(ab−cd)` (non-monomial multi-vertex): HH_0 = 4 pin + oracle.
4. Checkpoint payload: `tags` present, `len(tags[n]) == rks[n]`, pairs.
5. Memory wall: `stop_reason='memory'`, wall at degree 1, HH ≡ oracle under same budget.
6. `finalize_only` re-emits the checkpoint summary unchanged.
7. Resume: stop-early + continue ≡ fresh full run ≡ oracle.
8. Mode-mismatch (both directions) → `QuiverlabError`.
9. Local regression: `tests/engine/test_deepen.py` unchanged-green on numba AND
   `QUIVERLAB_NO_NUMBA=1` paths.

## Status

- [ ] Executed (fill in on completion)
