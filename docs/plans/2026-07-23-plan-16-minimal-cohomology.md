# Plan 16 — HH cohomology from the minimal/corner A^e resolution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** the minimal A^e engine — today homology-only — also computes deep
`dim HH^n(A; F_p)` for **any** f.d. algebra (local and multi-vertex), giving a second
deep-degree oracle against CS.

**Architecture:** apply `Hom_{A^e}(-, A)` to the SAME minimal projective resolution
`minimal_resolution` already builds. For a free term `(A^e)^{r_n}`,
`Hom_{A^e}((A^e)^{r_n}, A) ≅ A^{r_n}`; the coboundary `δ^{n-1} : A^{r_{n-1}} → A^{r_n}`
is precomposition with `d_n`, whose block entries act **two-sidedly the other way
round** than the homology collapse: coefficient `cf` at `e_u ⊗ e_v` sends
`α ↦ cf · e_u · α · e_v` (homology's `_contracted_degree` computes `e_v · α · e_u` —
the CS code calls these `a·w·b` vs `b·w·a`). On the corner path (Plan 13), the
cochain block of a generator tagged `(v, w)` is `Hom_{A^e}(A^e·(ε_v⊗ε_w), A) ≅
e_v A e_w` — the **opposite corner** of the homology target `e_w A e_v`, i.e. the
existing `ctx.cornerA[(i, j)] = e_j A e_i` dict looked up with the tag **swapped**.
`dim HH^n = dim C^n − rank δ^n − rank δ^{n-1}`.

**Tech Stack:** numpy int64 mod p, pure-Python loops (the corner path is pure by
design; `rank_mod_p` carries the numba/pure dispatch). Tests in `tests/engine/` →
deep bucket.

## Global Constraints

- No float literals in `src/` (AST gate `tests/test_no_floats.py`).
- Python is always `.venv/bin/python`; tests via
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q …`.
- Exact numba/pure agreement: run new tests under `QUIVERLAB_NO_NUMBA=1` too.
- Oracles live, never hardcoded: `scan3.hochschild_cohomology_dims` (dual normalized
  bar, returns `{p: [dim HH^0..HH^N]}`), `cs_cohomology_dims(A, top).dims` (CS coh
  side, core-Algebra input), plus the Happel/Künneth `[1,0,0]` literature pins.
- Homology side stays bit-for-bit untouched (`tests/engine/test_minimal_*.py` green).
- Truncation semantics identical to `minimal_homology_dims`: exact through
  `truncated_at − 1` (δ^t needs the unknown `d_{t+1}`).
- Conventional commits; green at every commit; merge/push only when Marco asks.

---

### Task 1: Local (free) cochain path — `_cohomology_degree` + `minimal_cohomology_dims`

**Files:**
- Create: `tests/engine/test_minimal_cohomology.py`
- Modify: `src/quiverlab/engine/resolutions_minimal.py` (append after
  `minimal_homology_dims`)

**Interfaces:**
- Consumes: `minimal_resolution(A, N, p, max_term_dim, max_transient_bytes)` →
  `(rks, cols, eng, truncated_at)`; `rank_mod_p(M, p)`; `eng.corner_ctx` /
  `eng.corner_tags` (None on the local path).
- Produces: `_cohomology_degree(eng, gens_n, r_nm1, n) -> np.ndarray` of shape
  `(m*r_n, m*r_nm1)` (the map `C^{n-1} → C^n` built from `d_n`);
  `minimal_cohomology_dims(A, N, primes=(32003,), max_term_dim=20000,
  max_transient_bytes=None) -> {p: [dim HH^0, ...]}`. Task 2 fills the corner
  branch with `_corner_cohomology_degree(eng, ctx, gens_n, tags_n, tags_nm1, p)`.

- [ ] **Step 1: Write the failing local-path tests** — create
  `tests/engine/test_minimal_cohomology.py`:

```python
"""Plan 16: HH cohomology from the minimal/corner A^e resolution (Hom-collapse).

The minimal engine was homology-only; Hom_{A^e}(-, A) on the SAME resolution
gives deep HH^. for any f.d. algebra. The coh-side collapse acts a.w.b (the
homology side is b.w.a), and on the corner path the block of a generator tagged
(v, w) is e_v A e_w -- the OPPOSITE corner of the homology target e_w A e_v.

Oracles: the dual normalized bar complex (scan3.hochschild_cohomology_dims,
live), the CS coh side (second deep engine), and the Happel/Kunneth [1,0,0]
pins."""
import pytest

import quiverlab as ql
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.hh_engine import truncated_polynomial
from quiverlab.engine.scan3 import quantum_ci, hochschild_cohomology_dims
from quiverlab.engine.resolutions_minimal import minimal_cohomology_dims

PRIMES = (32003, 2, 3, 5)


def _eng(vertices, arrows, relations, p=32003):
    Q = ql.Quiver(vertices, arrows)
    return to_engine(Q.algebra(relations=relations, field=ql.GF(p)))


def test_local_zoo_matches_bar():
    """k[x]/x^3 and quantum_ci(2) over four primes: minimal HH^. == bar HH^."""
    for make in (lambda: truncated_polynomial(3), lambda: quantum_ci(2)):
        A = make()
        for p in PRIMES:
            mc = minimal_cohomology_dims(A, 4, primes=(p,))[p]
            bc = hochschild_cohomology_dims(A, 4, primes=(p,))[p]
            assert mc == bc[:len(mc)], f"{A.name} p={p}: {mc} != {bc}"


def test_deep_cross_oracle_vs_cs():
    """Past the bar window: minimal coh == CS coh degreewise to depth 8 on the
    quantum CI over GF(3) -- two INDEPENDENT deep engines agreeing."""
    from quiverlab.resolutions_cs.homology import cs_cohomology_dims
    Q = ql.Quiver([1], {"x": (1, 1), "y": (1, 1)})
    A = Q.algebra(relations=["x*x", "y*y", "x*y + y*x"], field=ql.GF(3))
    cs = cs_cohomology_dims(A, 8).dims
    mc = minimal_cohomology_dims(to_engine(A), 8, primes=(3,))[3]
    assert mc == cs


def test_truncation_semantics():
    """A tiny term cap truncates: the dims list is SHORTER and an exact prefix,
    never padded or wrong (delta^t needs the unknown d_{t+1})."""
    A = quantum_ci(2)
    full = minimal_cohomology_dims(A, 6, primes=(32003,))[32003]
    trunc = minimal_cohomology_dims(A, 6, primes=(32003,), max_term_dim=40)[32003]
    assert len(trunc) < len(full)
    assert trunc == full[:len(trunc)]
```

- [ ] **Step 2: Run to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_minimal_cohomology.py -p no:cacheprovider`
Expected: FAIL at import — `cannot import name 'minimal_cohomology_dims'`

- [ ] **Step 3: Implement the local path** — append to
  `src/quiverlab/engine/resolutions_minimal.py` (after `minimal_homology_dims`,
  before `hochschild_dimension`):

```python
def _cohomology_degree(eng, gens_n, r_nm1, n):
    """Hom_{A^e}(-, A) applied to d_n: one degree of the cochain complex.
    `gens_n` = cols[n] (r_n generator blocks over (A^e)^{r_nm1}), `r_nm1` = rks[n-1].
    Returns delta^{n-1} : C^{n-1} = A^{r_nm1} -> C^n = A^{r_n}, alpha |-> alpha o d_n,
    of shape (m*r_n, m*r_nm1).  A coefficient cf at e_u (x) e_v in block blk of
    generator j sends e_a |-> cf * (e_u e_a) e_v -- the coh-side two-sided action;
    the homology collapse (_contracted_degree) is (e_v e_a) e_u."""
    m = eng.m
    m2 = eng.m2
    p = eng.p
    T = eng.T
    r_n = len(gens_n)
    M = np.zeros((m * r_n, m * r_nm1), dtype=np.int64)
    for j, g in enumerate(gens_n):
        for blk in range(r_nm1):
            w = g[blk * m2:(blk + 1) * m2]
            for a in range(m):
                col = blk * m + a
                acc = np.zeros(m, dtype=np.int64)
                for uu in range(m):
                    for vv in range(m):
                        cf = w[uu * m + vv]
                        if cf % p == 0:
                            continue
                        mid = T[uu, a, :]              # e_u * e_a
                        for s in np.nonzero(mid % p)[0]:
                            acc = (acc + cf * mid[s] * T[s, vv, :]) % p   # (e_u e_a) e_v
                M[j * m:(j + 1) * m, col] = (M[j * m:(j + 1) * m, col] + acc) % p
    return M


def minimal_cohomology_dims(A, N, primes=(32003,), max_term_dim=20000,
                            max_transient_bytes=None):
    """dim HH^n(A; F_p) for n=0..N via Hom_{A^e}(-, A) on the SAME minimal A^e
    resolution the homology side uses (rebuilt per prime).  Returns
    {p: [dim HH^0, ..., dim HH^M]} with M = N unless a budget truncated the build
    at degree t (then the list stops at t-1 and is exact: delta^t needs the
    unknown d_{t+1}).  Corner path (multi-vertex): the cochain block of a
    generator tagged (v, w) is e_v A e_w -- the homology dict cornerA[(i, j)] =
    e_j A e_i read with the tag SWAPPED."""
    out = {}
    for p in primes:
        rks, cols, eng, trunc = minimal_resolution(
            A, N, p, max_term_dim=max_term_dim, max_transient_bytes=max_transient_bytes)
        m = eng.m
        last = (trunc - 1) if trunc is not None else N
        dims = []
        ctx = getattr(eng, "corner_ctx", None)
        if ctx is not None:
            tags = eng.corner_tags
            D = {n: _corner_cohomology_degree(eng, ctx, cols.get(n, []) or [],
                                              tags.get(n, []), tags.get(n - 1, []), p)
                 for n in range(1, N + 2)}
            for n in range(0, last + 1):
                dimn = sum(ctx.corner_dim_A((tg[1], tg[0])) for tg in tags.get(n, []))
                rn = (rank_mod_p(D[n + 1], p)                       # rank delta^n
                      if rks.get(n + 1, 0) > 0 and rks.get(n, 0) > 0 else 0)
                rnm1 = (rank_mod_p(D[n], p)                         # rank delta^{n-1}
                        if n >= 1 and rks.get(n, 0) > 0 and rks.get(n - 1, 0) > 0 else 0)
                dims.append(int(dimn - rn - rnm1))
            out[p] = dims
            continue
        D = {n: _cohomology_degree(eng, cols.get(n, []) or [], rks.get(n - 1, 0), n)
             for n in range(1, N + 2)}
        for n in range(0, last + 1):
            dimn = m * rks.get(n, 0)
            rn = (rank_mod_p(D[n + 1], p)
                  if rks.get(n + 1, 0) > 0 and rks.get(n, 0) > 0 else 0)
            rnm1 = (rank_mod_p(D[n], p)
                    if n >= 1 and rks.get(n, 0) > 0 and rks.get(n - 1, 0) > 0 else 0)
            dims.append(int(dimn - rn - rnm1))
        out[p] = dims
    return out
```

(`_corner_cohomology_degree` lands in Task 2; the local tests never reach the
corner branch. Do not stub it — a corner call in this state must NameError loudly.)

- [ ] **Step 4: Run to verify the local tests pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_minimal_cohomology.py -p no:cacheprovider`
Expected: PASS (3 tests)

- [ ] **Step 5: Homology-side regression + pure path**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_minimal_resolution.py tests/engine/test_minimal_multivertex.py tests/engine/test_minimal_memory_guard.py -p no:cacheprovider`
Run: `QUIVERLAB_NO_NUMBA=1 .venv/bin/python -m pytest -q tests/engine/test_minimal_cohomology.py -p no:cacheprovider`
Expected: PASS both.

- [ ] **Step 6: Commit**

```bash
git add src/quiverlab/engine/resolutions_minimal.py tests/engine/test_minimal_cohomology.py
git commit -m "feat(engine): minimal A^e cohomology -- Hom-collapse on the minimal resolution (Plan 16, local path)"
```

### Task 2: Corner (multi-vertex) cochain path — `_corner_cohomology_degree`

**Files:**
- Modify: `src/quiverlab/engine/resolutions_minimal.py` (insert before
  `minimal_cohomology_dims`)
- Modify: `tests/engine/test_minimal_cohomology.py` (append)

**Interfaces:**
- Consumes: `_CornerContext` (`ctx.cornerA[(i, j)]` = `(m, d)` basis of `e_j A e_i`;
  `ctx.corner_dim_A(tag)`; `ctx.m`, `ctx.m2`), `_solve_in_span(B, v, p)`,
  `eng.T`, and Task 1's corner branch in `minimal_cohomology_dims` (which calls
  this with tags of degrees n and n−1).
- Produces: `_corner_cohomology_degree(eng, ctx, gens_n, tags_n, tags_nm1, p) ->
  np.ndarray` — the corner-typed `δ^{n-1} : C^{n-1} → C^n` with **coh corners**
  (`cornerA` read with swapped tags on both sides).

- [ ] **Step 1: Write the failing corner tests** — append to
  `tests/engine/test_minimal_cohomology.py`:

```python
def test_ka2_happel_pin():
    """kA_2 hereditary: HH^. = [1, 0, 0] (Happel). The coh corner of the P_1
    tag (1,2) is e_1 A e_2 (dim 1), NOT the homology corner e_2 A e_1 (dim 0):
    the tag swap is load-bearing here."""
    for p in PRIMES:
        E = _eng([1, 2], {"a": (1, 2)}, [], p=p)
        assert minimal_cohomology_dims(E, 2, primes=(p,))[p] == [1, 0, 0]


def test_commutative_square_kunneth_pin():
    """kQ/(ab - cd) = kA_2 (x) kA_2: HH^. = [1, 0, 0] (Kunneth; the qpa
    crosscheck fixture) -- non-monomial multi-vertex."""
    E = _eng([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)},
             ["a*b - c*d"])
    assert minimal_cohomology_dims(E, 2, primes=(32003,))[32003] == [1, 0, 0]


def test_cyclic_nakayama_matches_bar():
    """kZ_3/rad^2 over four primes: corner coh == bar coh degreewise (the
    strongest multi-vertex cross-check; nonzero in high degrees)."""
    for p in PRIMES:
        E = _eng([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)},
                 ["a*b", "b*c", "c*a"], p=p)
        mc = minimal_cohomology_dims(E, 3, primes=(p,))
        bc = hochschild_cohomology_dims(E, 3, primes=(p,))
        assert mc[p] == bc[p][:len(mc[p])], f"CN(3,2) p={p}: {mc[p]} != {bc[p]}"
```

- [ ] **Step 2: Run to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_minimal_cohomology.py -p no:cacheprovider`
Expected: the three new tests FAIL with `NameError: _corner_cohomology_degree`;
the Task-1 tests still PASS.

- [ ] **Step 3: Implement the corner path** — insert into
  `src/quiverlab/engine/resolutions_minimal.py` directly before
  `minimal_cohomology_dims`:

```python
def _corner_cohomology_degree(eng, ctx, gens_n, tags_n, tags_nm1, p):
    """One degree of the corner-typed cochain complex: Hom_{A^e}(A^e.(eps_v (x)
    eps_w), A) ~ e_v A e_w, the OPPOSITE corner of the homology collapse -- both
    row and column blocks read cornerA with the tag swapped.  Columns run over
    the C^{n-1} blocks (tags_nm1), rows over the C^n blocks (tags_n = one per
    generator); entries by the coh-side two-sided action e_u . alpha . e_v,
    coordinates recovered in the row block's coh-corner basis (reconstruction
    asserted -- loud on failure)."""
    m, m2 = ctx.m, ctx.m2
    T = eng.T % p
    row_offs, off = [], 0
    for tg in tags_n:
        row_offs.append(off)
        off += ctx.corner_dim_A((tg[1], tg[0]))
    nrows = off
    cols_out = []
    for blk, tgp in enumerate(tags_nm1):
        Bcol = ctx.cornerA[(tgp[1], tgp[0])]          # e_v A e_w for source tag (v, w)
        for c in range(Bcol.shape[1]):
            alpha = Bcol[:, c]
            col = np.zeros(nrows, dtype=np.int64)
            for j, (g, tg) in enumerate(zip(gens_n, tags_n)):
                w = g[blk * m2:(blk + 1) * m2]
                acc = np.zeros(m, dtype=np.int64)
                for uu in range(m):
                    for vv in range(m):
                        cf = w[uu * m + vv]
                        if cf % p == 0:
                            continue
                        ua = np.zeros(m, dtype=np.int64)
                        for a in np.nonzero(alpha)[0]:      # e_uu . alpha
                            ua = (ua + alpha[a] * T[uu, a, :]) % p
                        outv = np.zeros(m, dtype=np.int64)
                        for s in np.nonzero(ua)[0]:         # ... . e_vv
                            outv = (outv + ua[s] * T[s, vv, :]) % p
                        acc = (acc + cf * outv) % p
                x = _solve_in_span(ctx.cornerA[(tg[1], tg[0])], acc, p)
                assert x is not None, "corner cochain image left its corner (bug)"
                col[row_offs[j]:row_offs[j] + x.shape[0]] = x
            cols_out.append(col)
    if not cols_out:
        return np.zeros((nrows, 0), dtype=np.int64)
    return np.stack(cols_out, axis=1) % p
```

- [ ] **Step 4: Run the full new battery on both kernel paths**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_minimal_cohomology.py -p no:cacheprovider`
Run: `QUIVERLAB_NO_NUMBA=1 .venv/bin/python -m pytest -q tests/engine/test_minimal_cohomology.py -p no:cacheprovider`
Expected: PASS (6 tests) on both.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/engine/resolutions_minimal.py tests/engine/test_minimal_cohomology.py
git commit -m "feat(engine): corner-typed cochain path -- multi-vertex minimal HH^. (Plan 16)"
```

### Task 3: Docs + suites + backlog

**Files:**
- Modify: `CLAUDE.md` (minimal-engine bullet gains `minimal_cohomology_dims`; status
  paragraph gains Plan 16)
- Modify: `docs/plans/ROADMAP.md` (DELIVERED row 16)
- Modify: `docs/internals/05-resolutions.md` (cohomology paragraph in the minimal
  section: Hom-collapse, `a·w·b` vs `b·w·a`, swapped corner)
- Modify: `docs/plans/DEEPER-ENGINES-BACKLOG.md` (tick item 3 — done at branch start)

- [ ] **Step 1: Make the doc edits** (facts as in this plan's Architecture block)

- [ ] **Step 2: Full suites**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m deep -p no:cacheprovider`
Expected: PASS (~28 min)
Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast -p no:cacheprovider`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md docs/plans/ROADMAP.md docs/internals/05-resolutions.md docs/plans/DEEPER-ENGINES-BACKLOG.md
git commit -m "docs: Plan-16 status -- minimal/corner HH cohomology delivered"
```

## Validation matrix (tests/engine/test_minimal_cohomology.py)

1. Local zoo (`k[x]/x^3`, `quantum_ci(2)`) ≡ bar coh over 4 primes.
2. Deep cross-oracle: minimal ≡ CS coh to depth 8 on the GF(3) quantum CI —
   two independent deep engines, degreewise.
3. Truncation: capped run yields a strictly shorter exact prefix.
4. `kA_2` Happel pin `[1,0,0]` over 4 primes — the corner-swap is load-bearing
   (coh corner `e_1 A e_2` is 1-dim where the homology corner is 0).
5. Commutative square Künneth pin `[1,0,0]` (non-monomial multi-vertex).
6. `kZ_3/rad²` ≡ bar over 4 primes (corner path, nonzero deep degrees).
7. Homology-side regression: existing `test_minimal_*` suites untouched-green on
   both kernel paths.

## Status

- [ ] Executed (fill in on completion)
