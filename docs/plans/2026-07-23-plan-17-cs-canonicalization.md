# Plan 17 — CS canonicalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** the CS differential's correction solve is reduced to an explicit normal form
modulo its nullspace, making CS differentials **byte-reproducible by construction** —
which justifies flipping the 7 `xfail(strict=False)` coefficient pins strict
(Plan-04 stretch item E2, backlog Tier-1 item 4).

**Architecture:** `_d_general` solves `M·γ = rhs` for the order-condition correction
γ, unique only modulo `Null(M)` (on the quantum-CI family the nullity grows with
degree). Today `fields.linalg.solve` happens to return the free-variables-zero
particular solution and that *accidentally* matches the bank's hand-derived closed
forms — the `test_battery_bank_oracle.py` docstring explicitly forbids strict pins
until the representative is pinned. The fix: an explicit
`reduce_mod_nullspace(x, A, dom)` in `fields/linalg.py` — the unique coset
representative with zero coordinates at every free (non-pivot) column of `A`'s RREF —
applied to γ in `_d_general` right after the solve. It is a **no-op on today's
solver output** (so all byte pins keep passing) and a **guarantee against any future
solver change** (pivot order, free-variable convention). Coset moves preserve
`M·γ = rhs` (d²=0) and the order condition (every generator in `_lower_generators`
is `≺ σ` by construction, so *any* coefficient vector over them satisfies it).

**Tech Stack:** exact Domain arithmetic (`fields/linalg.py` `rref`/`nullspace`/
`solve` conventions). Tests in `tests/resolutions_cs/` → deep bucket.

## Global Constraints

- No float literals in `src/` (AST gate).
- Python is always `.venv/bin/python`; tests via
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q …`.
- The canonicalization must be a **no-op on current outputs**: every existing byte
  pin (bank oracle, paper d₂ pin) must pass *unchanged* — a diff means the
  canonical form was defined wrong, not that the pins need editing.
- d²=0 + order-condition gates stay per-instance (`assert_dd_zero`,
  `assert_order_condition`) and must stay green.
- After this plan the deep suite reports **0 xpassed** from the CS pins (they
  become plain passes; the suite-wide count drops from 7 to 0).
- Conventional commits; green at every commit; merge/push only when Marco asks.

---

### Task 1: `reduce_mod_nullspace` + canonicalization unit battery

**Files:**
- Create: `tests/resolutions_cs/test_canonicalization.py`
- Modify: `src/quiverlab/fields/linalg.py` (append after `solve`)

**Interfaces:**
- Consumes: `rref(rows, dom) -> (R, pivots)`, `nullspace(rows, dom)`,
  `solve(A, b, dom)` — all existing in `fields/linalg.py`.
- Produces: `reduce_mod_nullspace(x, A, dom) -> list` — the canonical coset
  representative; Task 2 wires it into `_d_general`.

- [ ] **Step 1: Write the failing tests** — create
  `tests/resolutions_cs/test_canonicalization.py`:

```python
"""Plan 17: the CS correction solve is canonicalized modulo its nullspace.

reduce_mod_nullspace(x, A, dom) is the unique element of x + Null(A) with zero
coordinates at every free (non-pivot) column of A's RREF.  solve() already
returns that representative (free variables set to 0), so canonicalization is a
NO-OP today -- the point is that it is now an explicit, tested guarantee instead
of a solver-convention accident (see the WARNING block that used to live in
test_battery_bank_oracle.py).  The adversarial test in Task 2 proves the CS
differential bytes no longer depend on WHICH solution the solver returns."""
import pytest

from quiverlab import Quiver, GF
from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace, reduce_mod_nullspace, solve


def _dom_cases():
    return (GF(5), QQ)


def test_coset_invariance_and_idempotence():
    """Every solution of A y = b canonicalizes to the SAME vector; applying the
    reduction twice equals applying it once."""
    for dom in _dom_cases():
        i = dom.coerce
        # rank-2 system with a 2-dim nullspace (4 unknowns)
        A = [[i(1), i(2), i(0), i(1)],
             [i(0), i(0), i(1), i(3)],
             [i(1), i(2), i(1), i(4)]]          # row3 = row1 + row2 (dependent)
        b = [i(1), i(2), i(3)]
        x0 = solve(A, b, dom)
        assert x0 is not None
        canon = reduce_mod_nullspace(x0, A, dom)
        assert canon == reduce_mod_nullspace(canon, A, dom)      # idempotent
        for v in nullspace(A, dom):
            shifted = [dom.add(a_, b_) for a_, b_ in zip(x0, v)]
            assert reduce_mod_nullspace(shifted, A, dom) == canon  # coset-invariant


def test_solver_output_is_already_canonical():
    """solve()'s free-variables-zero particular solution IS the canonical
    representative -- the no-op property that keeps every byte pin passing."""
    for dom in _dom_cases():
        i = dom.coerce
        A = [[i(1), i(2), i(0), i(1)],
             [i(0), i(0), i(1), i(3)]]
        b = [i(4), i(1)]
        x0 = solve(A, b, dom)
        assert reduce_mod_nullspace(x0, A, dom) == x0


def test_full_rank_is_untouched():
    """No nullspace -> the reduction returns the input unchanged."""
    dom = GF(7)
    i = dom.coerce
    A = [[i(1), i(1)], [i(0), i(1)]]
    b = [i(3), i(2)]
    x0 = solve(A, b, dom)
    assert reduce_mod_nullspace(x0, A, dom) == x0
```

- [ ] **Step 2: Run to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/resolutions_cs/test_canonicalization.py -p no:cacheprovider`
Expected: FAIL at import — `cannot import name 'reduce_mod_nullspace'`

- [ ] **Step 3: Implement** — append to `src/quiverlab/fields/linalg.py` after
  `solve`:

```python
def reduce_mod_nullspace(x, A, dom):
    """The canonical representative of the coset x + Null(A): the unique element
    with ZERO coordinates at every free (non-pivot) column of A's RREF.

    Applying this to ANY solution of A·y = b yields the same vector, so a linear
    solve's answer becomes independent of the solver's particular-solution
    convention (pivot order, free-variable assignment).  solve() already returns
    this representative (free variables set to 0); this function is the explicit,
    order-pinned guarantee the CS correction solve canonicalizes through
    (Plan 17 -- byte-reproducible CS differentials)."""
    if not A or not x:
        return list(x)
    R, pivots = rref(A, dom)
    nc = len(A[0])
    y = list(x)
    for fc in (c for c in range(nc) if c not in pivots):
        c = y[fc]
        if dom.is_zero(c):
            continue
        # subtract c * v_fc where v_fc is nullspace()'s basis vector for the free
        # column fc: v_fc[fc] = 1, v_fc[pc] = -R[r][fc] at each pivot column pc
        y[fc] = dom.zero()
        for r, pc in enumerate(pivots):
            y[pc] = dom.add(y[pc], dom.mul(c, R[r][fc]))
    return y
```

- [ ] **Step 4: Run to verify they pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/resolutions_cs/test_canonicalization.py -p no:cacheprovider`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/fields/linalg.py tests/resolutions_cs/test_canonicalization.py
git commit -m "feat(fields): reduce_mod_nullspace -- canonical coset representative for underdetermined solves (Plan 17)"
```

### Task 2: Wire into `_d_general` + adversarial-solver battery

**Files:**
- Modify: `src/quiverlab/resolutions_cs/resolution.py` (`_d_general`, after the
  `gamma is None` check)
- Modify: `tests/resolutions_cs/test_canonicalization.py` (append)

**Interfaces:**
- Consumes: Task 1's `reduce_mod_nullspace(x, A, dom)`;
  `ChouhySolotarResolution._solve(M, rhs, ncols)` / `._d_cache` / `.d_terms(n, chain)`
  / `.matrix(n, side)`; `fields.linalg.nullspace`.
- Produces: canonicalized `gamma` inside `_d_general` — CS differential bytes
  independent of the solver's choice. Task 3 relies on this to flip the pins.

- [ ] **Step 1: Write the failing adversarial tests** — append to
  `tests/resolutions_cs/test_canonicalization.py`:

```python
def _qci(field=None):
    from quiverlab.groebner import build_reduction_system
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    f = field or GF(5)
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", "y*x - 2*x*y"]
    return ChouhySolotarResolution(Q.algebra(relations=rels, field=f),
                                   build_reduction_system(Q, rels, f), max_degree=7)


def test_correction_nullity_is_nonzero_somewhere():
    """The adversarial battery is NOT vacuous: on the quantum CI some degree has a
    correction system with genuine nullspace freedom (the bank docstring measured
    nullity growing with degree)."""
    res = _qci()
    assert max(_correction_nullities(res, upto=6)) > 0


def _correction_nullities(res, upto):
    """Nullities of the correction systems actually solved while building d_1..d_upto."""
    from quiverlab.resolutions_cs.pelt import terms_to_pelt, apply_lower
    dom = res.dom
    out = [0]
    for n in range(1, upto + 1):
        for chain in res.ss.S(n):
            gens = res._lower_generators(n, chain)
            if not gens:
                continue
            delta = res.delta_terms(n, chain)
            rhs_pe = apply_lower(res, n, terms_to_pelt(res, delta))
            cols = [apply_lower(res, n, terms_to_pelt(res, [g])) for g in gens]
            keys = sorted(set(rhs_pe) | {k for col in cols for k in col})
            M = [[col.get(k, dom.zero()) for col in cols] for k in keys]
            out.append(len(nullspace(M, dom)) if M else len(gens))
    return out


def test_cs_bytes_survive_an_adversarial_solver(monkeypatch):
    """THE Plan-17 theorem: shift the solver's answer by a nullspace vector (a
    DIFFERENT valid solution) -- the built differentials must be byte-identical
    anyway, because _d_general canonicalizes.  Before this plan the bytes moved
    (the bank-oracle WARNING); now they cannot."""
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    baseline = _qci()
    ref = {n: [sorted((baseline.to_int(c), a, t, cc)
                      for (c, a, t, cc) in baseline.d_terms(n, q))
               for q in baseline.ss.S(n)] for n in range(1, 7)}

    real_solve = ChouhySolotarResolution._solve
    shifted = {"count": 0}

    def adversarial_solve(self, M, rhs, ncols):
        sol = real_solve(self, M, rhs, ncols)
        if sol is None or not M:
            return sol
        basis = nullspace(M, self.dom)
        if not basis:
            return sol
        shifted["count"] += 1
        return [self.dom.add(s, v) for s, v in zip(sol, basis[0])]

    monkeypatch.setattr(ChouhySolotarResolution, "_solve", adversarial_solve)
    perturbed = _qci()
    got = {n: [sorted((perturbed.to_int(c), a, t, cc)
                      for (c, a, t, cc) in perturbed.d_terms(n, q))
               for q in perturbed.ss.S(n)] for n in range(1, 7)}
    assert shifted["count"] > 0            # the adversary actually fired
    assert got == ref                      # ...and the bytes did not move
    perturbed.assert_dd_zero(upto=6, side="hom")
    perturbed.assert_order_condition(upto=6)
```

- [ ] **Step 2: Run to verify the adversarial test fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/resolutions_cs/test_canonicalization.py -p no:cacheprovider`
Expected: `test_cs_bytes_survive_an_adversarial_solver` FAILS on `got == ref`
(the shifted solution changes bytes — canonicalization not wired yet);
the nullity probe and Task-1 tests PASS.

- [ ] **Step 3: Wire the canonicalization** — in
  `src/quiverlab/resolutions_cs/resolution.py::_d_general`, directly after the
  `gamma is None` raise:

```python
        from quiverlab.fields.linalg import reduce_mod_nullspace
        # canonical representative: γ is unique only mod Null(M); pin the
        # free-variables-zero coset element so d_n is byte-reproducible no matter
        # which particular solution _solve returned (Plan 17; flips the bank-oracle
        # byte pins strict)
        gamma = reduce_mod_nullspace(gamma, M, dom)
```

- [ ] **Step 4: Run to verify everything passes**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/resolutions_cs/test_canonicalization.py -p no:cacheprovider`
Expected: PASS (5 tests)

- [ ] **Step 5: CS regression battery (bytes must be a no-op)**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/resolutions_cs/ -p no:cacheprovider`
Expected: PASS with the 7 pins still XPASS (no byte moved). Any xfail here means
the canonical form disagrees with the solver's representative — fix the
implementation, do NOT touch the pins.

- [ ] **Step 6: Commit**

```bash
git add src/quiverlab/resolutions_cs/resolution.py tests/resolutions_cs/test_canonicalization.py
git commit -m "feat(cs): canonicalize the correction solve modulo its nullspace -- byte-reproducible differentials (Plan 17)"
```

### Task 3: Flip the 7 pins strict + docs

**Files:**
- Modify: `tests/resolutions_cs/test_battery_bank_oracle.py` (drop the 2 xfail
  markers; replace the WARNING block)
- Modify: `tests/resolutions_cs/test_differential.py` (drop the xfail marker;
  update the docstring)
- Modify: `CLAUDE.md`, `docs/plans/ROADMAP.md`, `docs/internals/09-chouhy-solotar.md`,
  `docs/plans/DEEPER-ENGINES-BACKLOG.md` (item 4 tick — done at branch start)

- [ ] **Step 1: Remove the two `@pytest.mark.xfail(strict=False, reason=
  "canonicalization pending (correction unique mod nullspace)")` decorators in
  `test_battery_bank_oracle.py`; replace the docstring's EMPIRICAL FINDING /
  WARNING block with a short "canonicalization landed (Plan 17)" note: the pins
  are now strict because `_d_general` reduces γ to the free-variables-zero coset
  representative (`reduce_mod_nullspace`), so byte identity is a guarantee, not a
  tie-breaking coincidence; the adversarial-solver test in
  `test_canonicalization.py` enforces it.

- [ ] **Step 2: Remove the xfail decorator on
  `test_qci_d3_correction_matches_paper` in `test_differential.py`; reword its
  docstring: the pin is now binding — the canonicalized correction reproduces the
  paper's d₂ verbatim.**

- [ ] **Step 3: Verify the flips**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/resolutions_cs/ -p no:cacheprovider`
Expected: PASS with **0 xpassed** (the 7 former pins are plain passes now).

- [ ] **Step 4: Doc edits** (facts per the Architecture block): CLAUDE.md CS bullet
  gains "differentials canonical/byte-reproducible (Plan 17)" + status paragraph;
  ROADMAP DELIVERED row 17; internals 09 gains a "canonical form" paragraph
  (coset representative, no-op-by-construction, adversarial gate).

- [ ] **Step 5: Full suites**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m deep -p no:cacheprovider`
Expected: PASS, **0 xpassed** (was 7).
Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast -p no:cacheprovider`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add tests/resolutions_cs/test_battery_bank_oracle.py tests/resolutions_cs/test_differential.py CLAUDE.md docs/plans/ROADMAP.md docs/internals/09-chouhy-solotar.md docs/plans/DEEPER-ENGINES-BACKLOG.md
git commit -m "test(cs): flip the 7 byte pins strict -- canonicalization guarantees the representative (Plan 17)"
```

## Validation matrix (tests/resolutions_cs/test_canonicalization.py + flipped pins)

1. Coset invariance + idempotence of `reduce_mod_nullspace` over GF(5) and QQ.
2. `solve()`'s output is already canonical (the no-op property).
3. Full-rank systems pass through untouched.
4. Non-vacuity: the qci correction systems really have nullity > 0.
5. **Adversarial solver**: shifting the solve by a nullspace vector leaves every
   d_n byte-identical, and d²=0/order still hold.
6. The 7 former xfail pins (2×3 bank byte pins + the paper d₂ pin) pass as
   **strict** tests; deep suite xpassed count drops 7 → 0.
7. Full CS battery + deep/fast suites green.

## Status

- [x] Executed 2026-07-23 in-session (branch `plan-17-cs-canonicalization`).
  Canonicalization is a verified no-op on prior outputs (full CS battery byte-stable);
  the adversarial-solver gate holds; the 7 pins are strict (CS battery: 114 passed,
  0 xpassed); full `-m deep` and `-m fast` suites green.
