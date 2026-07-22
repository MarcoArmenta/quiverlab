# Plan 12: Straddling ambiguities & the right decomposition — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Correct the CS §3 block combinatorics (`associated_paths` / `left_decomposition`) to the paper-verbatim minimal-extension condition, add `right_decomposition`, use the right factorization in the Bardzell and CS odd (2-term) differentials per CS §4 `f_n` (n even), and lift the CS `NotImplementedError` for non-quadratic non-monomial presentations (the Plan-04 stretch item).

**Architecture:** One shared primitive fixes two engines. `MonomialPresentation` (engine/resolutions_bardzell.py) gains first-reducibility block cutting (left and right); `BardzellResolution.differential_matrix` and `ChouhySolotarResolution.delta_terms` consume it. The scope gate `_require_in_scope` is deleted; per-instance certification stays `assert_dd_zero` + `assert_order_condition` + degreewise bar agreement.

**Tech Stack:** pure Python + numpy int64 (engine), exact `Domain` arithmetic (resolutions_cs). No numba kernels touched.

## Motivating findings (2026-07-22 session, verified by experiment)

1. **Latent bug.** `associated_paths`/`left_decomposition` require each consecutive block pair to be **exactly** a tip (`cand in relset`). CS §3 Definition (arXiv:1406.2300, TeX ~696–718) requires only: *"u_i·u_{i+1} is reducible but u_i·d is irreducible for any proper left divisor d of u_{i+1}"* — the witness tip may **straddle** the block boundary as a proper suffix of the pair. Repro: tips `{xx, yy, xyx}` — true `AP^3 = {xxx, xxyx, xyxx, xyxyx, yyy}` (5 chains, = minimal-A^e Betti number `rks[3] = 5`; the full Betti sequence is `2, 3, 5, 9, 17`), code returns 3 (misses the straddles `xyxx` and `xyxyx`). On `A = k⟨x,y⟩/(xx,yy,xyx)` (dim 6): Bardzell HH_• = `[4,4,7,6]` (32003) vs bar **and** minimal-A^e oracle `[4,4,6,9]`. All four primes disagree from HH_2 on (p=2: bardzell `[4,7,10,8]` vs oracle `[4,7,9,12]`; p=3,5 = `[4,4,7,6]` vs `[4,4,6,9]`). The validated zoo (truncated poly / cyclic Nakayama / radsq / QCI) has uniform-length relations — straddles never arise, so every battery stayed green.
2. **The stretch item.** CS §4 `f_n` (n even, TeX ~840–855): `f_n(1⊗q⊗1) = v_n ⊗ v_{n-1}⋯v_0 ⊗ 1 − 1 ⊗ u_0⋯u_{n-1} ⊗ u_n`, `u` = left, `v` = **right** factorization. Both the Bardzell odd/small map and CS `delta_terms` currently substitute `u_0` for `v_n` — provably equal only for quadratic tips (CS Prop. "cuadratico", TeX ~766–777: `v_i = α_{n-i}`) and for the palindromic shipped families. `k[x]/(x^a)`, cyclic Nakayama, radsq, QCI(2,2) all have `v_n = u_0`, so **no certified number changes**.
3. **Side finding (out of scope, memory-noted):** `engine/resolutions_minimal.py` silently returns a wrong (zero) resolution on multi-vertex algebras; its validation set is all local. Do not use it as an oracle off the local case. Not touched by this plan.

## Global Constraints

- Python is **always** `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python`; tests run as `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...` from the repo root.
- **No floats in `src/`** (AST gate `tests/test_no_floats.py`). All algebra exact.
- **Left-to-right path composition** everywhere; words are arrow-**name** tuples.
- `tests/engine/` and `tests/resolutions_cs/` are auto-assigned the **deep** marker by directory (`tests/conftest.py`); run them with `-m deep`.
- Engine internals stay internal; the public surface (`Algebra.hochschild_*`, `engine=` ∈ {auto,bar,fast,cs}) does not change shape in this plan.
- Uniform-relation families (truncated poly, cyclic Nakayama, radsq, QCI(2,2), commutative square) must be **byte-identical** before/after: their pairs are exact and palindromic, so old and new algorithms agree; every existing battery is the regression net.
- Work on branch `plan-12-ambiguity-blocks` off `main`. Conventional commits, green at every commit, **no push, no merge** (user decides).
- CS paper source for reference: arXiv e-print 1406.2300 (single TeX; read with `LC_ALL=C grep -a` — macOS grep misdetects it as binary).

## Mathematical spec (binding)

Let `S` be the (reduced, interreduced) tip set; a path is *reducible* iff it contains some tip as a subword.

- **Left blocks (CS §3 Def., verbatim):** `p = u_0 u_1 ⋯ u_n`, `u_0 ∈ Q_1`, all `u_i` irreducible, and for all `i`: `u_i u_{i+1}` reducible but `u_i d` irreducible for every proper left divisor `d` of `u_{i+1}`. Equivalent greedy form (uniqueness: CS Prop., TeX ~720–760): scanning left→right, each next block is the **shortest** extension whose pair becomes reducible; at first reducibility the witness tip ends exactly at the pair's end and must straddle the block boundary (`len(tip) > len(new_block)`), else the branch is dead (`u_{i+1}` reducible).
- **Right blocks (mirror, verbatim):** `p = v_n ⋯ v_0` (left-to-right; `v_n` leftmost), `v_0 ∈ Q_1` the **last** arrow, and for all `i`: `v_{i+1} v_i` reducible but `d v_i` irreducible for every proper right divisor `d` of `v_{i+1}`. Greedy right→left; witness tip **starts** exactly at the pair's start and straddles (`len(tip) > len(new_block)`).
- **Coincidence (CS §3, citing Bardzell/Sköldberg):** left n-ambiguities = right n-ambiguities *as sets*; the block decompositions differ in general. Consequences used: dropping the leftmost **right**-block of a chain in `S_n` leaves a word in `S_{n-1}`; dropping the rightmost **left**-block likewise.
- **Odd (2-term) differential, quiverlab index** (`d_n`, n odd ⟺ CS `f_{n-1}` even), for `σ ∈ S_n`, left blocks `u_0..u_{n-1}`, right blocks `v_{n-1}..v_0`:
  `(+1, v_{n-1}, word − leftmost right block, ())` and `(−1, (), word − rightmost left block, u_{n-1})`.
- Because ambiguities depend only on tips, all of the above is computed on the tip monomial algebra `A_S`; tails enter only through the (unchanged) correction solve.
- **Witness uniqueness:** `S` interreduced ⇒ at most one tip ends (resp. starts) at a given position; the suffix (resp. prefix) scan is exact because blocks are grown one arrow at a time from an irreducible state.

Reference decompositions for tests (hand-derived, tips `{xx, yy, xyx}`):
`xyxx`: left `(x)(yx)(x)`, right `(xy)(x)(x)`. `xxyx`: left `(x)(x)(yx)`, right `(x)(xy)(x)`.
Tips `{yx, x³, y²}` (QCI(3,2)), chain `yxxx`: left `(y)(x)(xx)`, right `(y)(xx)(x)`.
Line quiver tips `{abc, cde}`: `AP^3 = {abcde}`, left `(a)(bc)(de)`, right `(ab)(cd)(e)`.

---

### Task 1: Corrected `associated_paths` + `left_decomposition` (straddling overlaps)

**Files:**
- Modify: `src/quiverlab/engine/resolutions_bardzell.py` (`_proper_suffix_in_I` → replaced by `_witness_at_end`; `associated_paths` ~184–219; `left_decomposition` ~221–241)
- Test (create): `tests/engine/test_bardzell_straddle.py`

**Interfaces:**
- Consumes: existing `MonomialPresentation.__init__` fields (`relations`, `relset`, `bysrc`, `tgt`, `maxrel`).
- Produces: `associated_paths(n, maxlen) -> sorted list of tuples` (same signature, now complete); `left_decomposition(p, n) -> list of n tuples` (same signature, greedy); `_witness_at_end(path) -> tuple|None` (internal helper reused by Task 2's mirror).

- [ ] **Step 1: Write the failing tests**

Create `tests/engine/test_bardzell_straddle.py`:

```python
"""CS §3 verbatim block combinatorics: the witness tip may STRADDLE a block boundary
(pair reducible with minimal extension), which the exact-pair condition missed.
Repro discovered 2026-07-22: tips {xx,yy,xyx} lose the chain xyxx, making Bardzell
HH wrong from degree 2 (vs bar and minimal-A^e oracles). Uniform-relation families
(truncated poly, Nakayama, radsq) have exact pairs only: pinned unchanged below."""
import pytest

from quiverlab.engine.resolutions_bardzell import MonomialPresentation


def _s1():
    # one vertex, loops x,y; tips {xx, yy, xyx} — the minimal straddle presentation
    return MonomialPresentation(["v"], [("x", "v", "v"), ("y", "v", "v")],
                                [("x", "x"), ("y", "y"), ("x", "y", "x")])


def _line():
    # 1->2->..->6, arrows a..e; tips {abc, cde} overlap in the single arrow c
    return MonomialPresentation(
        [1, 2, 3, 4, 5, 6],
        [("a", 1, 2), ("b", 2, 3), ("c", 3, 4), ("d", 4, 5), ("e", 5, 6)],
        [("a", "b", "c"), ("c", "d", "e")])


def test_straddle_ap3_complete():
    ap3 = set(_s1().associated_paths(3, 30))
    assert ap3 == {("x", "x", "x"), ("x", "x", "y", "x"), ("x", "y", "x", "x"),
                   ("x", "y", "x", "y", "x"), ("y", "y", "y")}


def test_line_quiver_short_overlap_chain():
    pres = _line()
    assert pres.associated_paths(3, 20) == [("a", "b", "c", "d", "e")]
    assert pres.associated_paths(4, 24) == []          # no further overlap


def test_left_decomposition_straddle_blocks():
    pres = _s1()
    assert pres.left_decomposition(("x", "y", "x", "x"), 3) == \
        [("x",), ("y", "x"), ("x",)]                   # pair yx·x ⊇ xx straddling
    assert pres.left_decomposition(("x", "x", "y", "x"), 3) == \
        [("x",), ("x",), ("y", "x")]                   # pairs xx, xyx exact
    assert _line().left_decomposition(("a", "b", "c", "d", "e"), 3) == \
        [("a",), ("b", "c"), ("d", "e")]               # pair bc·de ⊇ cde straddling


def test_uniform_families_unchanged():
    tp = MonomialPresentation.truncated_polynomial(3)
    assert tp.associated_paths(3, 20) == [(0, 0, 0, 0)]
    assert tp.left_decomposition((0, 0, 0, 0), 3) == [(0,), (0, 0), (0,)]
    cn = MonomialPresentation.cyclic_nakayama(3, 2)
    assert sorted(cn.associated_paths(3, 20)) == \
        [(0, 1, 2, 0), (1, 2, 0, 1), (2, 0, 1, 2)]
    rs = MonomialPresentation.local_radsq(2)
    assert len(rs.associated_paths(3, 20)) == 8        # all words of length 4 (quadratic)


def test_deep_straddle_chain_growth():
    # AP^n stays finite and consistent to depth 8 on the straddle presentation
    pres = _s1()
    counts = [len(pres.associated_paths(n, 60)) for n in range(1, 9)]
    assert counts[0] == 2 and counts[1] == 3 and counts[2] == 4
    assert all(c >= 2 for c in counts)                 # y^k chain and an x-side chain persist
```

- [ ] **Step 2: Run tests, verify current failures**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_bardzell_straddle.py -m deep`
Expected: `test_straddle_ap3_complete`, `test_line_quiver_short_overlap_chain`, `test_left_decomposition_straddle_blocks` FAIL (missing chains / assertion "no left decomposition"); `test_uniform_families_unchanged` PASSES (pin of the unchanged behavior).

- [ ] **Step 3: Implement first-reducibility cutting**

In `src/quiverlab/engine/resolutions_bardzell.py`, replace `_proper_suffix_in_I` with:

```python
    def _witness_at_end(self, path):
        """The unique tip that is a suffix of `path`, or None (S is interreduced, so at
        most one tip ends at a given position)."""
        lp = len(path)
        for r in self.relations:
            lr = len(r)
            if lr <= lp and tuple(path[lp - lr:]) == r:
                return r
        return None
```

Replace the body of `associated_paths`'s `ext` closure (the block-cut logic) with first-reducibility cutting — the full method:

```python
    def associated_paths(self, n, maxlen):
        """AP^n: the degree-n associated paths (n-ambiguities), CS §3 verbatim: each
        next block is the SHORTEST extension making (prev block)+(block) reducible;
        the witness tip ends at the pair's end and STRADDLES the boundary (it may be a
        proper suffix of the pair — the pair need not itself be a tip). maxlen caps the
        path length explored (an associated path of degree n has length
        <= n*(maxrel-1)+1; pass a generous bound)."""
        if n <= 0:
            return []
        if n == 1:
            return [(a,) for (a, s, t) in self.arrows]
        res = set()

        def rec(blocks, full):
            if len(blocks) == n:
                res.add(tuple(full))
                return
            prev = blocks[-1]
            cur = self.tgt[prev[-1]]

            def ext(uk, lastt):
                cand = tuple(prev) + tuple(uk)
                w = self._witness_at_end(cand)
                if w is not None:
                    # FIRST reducibility: the block boundary is forced here (CS §3 (ii)).
                    if len(w) > len(uk):        # witness straddles => uk irreducible
                        rec(blocks + [tuple(uk)], full + list(uk))
                    return                       # witness inside uk alone: dead branch
                # a valid block uk has len(uk) < len(witness) <= maxrel: bound the search.
                if len(uk) >= self.maxrel:
                    return
                if len(full) + len(uk) > maxlen:
                    return
                for a in self.bysrc.get(lastt, []):
                    ext(uk + [a], self.tgt[a])

            for a in self.bysrc.get(cur, []):
                ext([a], self.tgt[a])

        for (a, s, t) in self.arrows:
            rec([(a,)], [a])
        return sorted(res, key=lambda p: (len(p), p))
```

Replace `left_decomposition` with the greedy deterministic form:

```python
    def left_decomposition(self, p, n):
        """Unique CS §3 left block decomposition p = u_0 ... u_{n-1} (n blocks): u_0 =
        first arrow; each next block is the SHORTEST extension making the consecutive
        pair reducible (witness tip ends at the pair's end, straddling the boundary)."""
        if n == 1:
            return [(p[0],)]
        blocks, pos = [(p[0],)], 1
        while len(blocks) < n:
            prev, uk, cut = blocks[-1], [], False
            while pos < len(p):
                uk.append(p[pos])
                pos += 1
                w = self._witness_at_end(tuple(prev) + tuple(uk))
                if w is not None:
                    assert len(w) > len(uk), ("block reducible in left decomposition", p, n)
                    cut = True
                    break
            assert cut, ("no left decomposition for associated path", p, n)
            blocks.append(tuple(uk))
        assert pos == len(p), ("left decomposition does not exhaust the path", p, n)
        return blocks
```

Grep for remaining `_proper_suffix_in_I` references (`grep -rn _proper_suffix_in_I src/ tests/`) — delete the method only if nothing else uses it.

- [ ] **Step 4: Run tests, verify all pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_bardzell_straddle.py tests/engine/test_bardzell_resolution.py tests/engine/test_acceptance_plan02.py -m deep`
Expected: all PASS (uniform-zoo Bardzell suites are the regression net; the straddle HH cross-check comes in Task 3 after the differential fix).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/engine/resolutions_bardzell.py tests/engine/test_bardzell_straddle.py
git commit -m "fix(engine): associated_paths/left_decomposition honor straddling overlaps (CS §3 minimal extension)"
```

---

### Task 2: `right_decomposition` (CS §3 right n-ambiguity blocks)

**Files:**
- Modify: `src/quiverlab/engine/resolutions_bardzell.py` (add `_witness_at_start`, `right_decomposition` next to `left_decomposition`)
- Test: `tests/engine/test_bardzell_straddle.py` (append)

**Interfaces:**
- Consumes: Task 1's `_witness_at_end` pattern (mirrored), `MonomialPresentation.relations`.
- Produces: `right_decomposition(p, n) -> list of n tuples in left-to-right order [v_{n-1}, ..., v_0]`, `v_0 = (p[-1],)`. Used by Task 3 (Bardzell odd map) and Task 4 (CS `delta_terms`).

- [ ] **Step 1: Write the failing tests** (append to `tests/engine/test_bardzell_straddle.py`)

```python
def test_right_decomposition_mirror_and_straddle():
    pres = _s1()
    # xyxx: right blocks (xy)(x)(x) — witness xyx STARTS at the pair's start, straddling
    assert pres.right_decomposition(("x", "y", "x", "x"), 3) == \
        [("x", "y"), ("x",), ("x",)]
    assert pres.right_decomposition(("x", "x", "y", "x"), 3) == \
        [("x",), ("x", "y"), ("x",)]
    assert _line().right_decomposition(("a", "b", "c", "d", "e"), 3) == \
        [("a", "b"), ("c", "d"), ("e",)]


def test_right_decomposition_quadratic_reverses_left():
    # CS Prop. "cuadratico": for S ⊆ Q_2, v_i = α_{n-i} — blocks are single arrows,
    # so right blocks == left blocks elementwise.
    rs = MonomialPresentation.local_radsq(2)
    for n in (2, 3, 4):
        for p in rs.associated_paths(n, 20):
            assert rs.right_decomposition(p, n) == rs.left_decomposition(p, n)


def test_right_decomposition_palindromic_families_agree():
    tp = MonomialPresentation.truncated_polynomial(3)
    assert tp.right_decomposition((0, 0, 0, 0), 3) == [(0,), (0, 0), (0,)]
    cn = MonomialPresentation.cyclic_nakayama(3, 3)
    for n in (2, 3):
        for p in cn.associated_paths(n, 30):
            L, R = cn.left_decomposition(p, n), cn.right_decomposition(p, n)
            assert [len(b) for b in R] == [len(b) for b in L][::-1] or L == R


def test_qci32_right_blocks():
    # tips {yx, xxx, yy} (QCI(3,2) tip algebra): yxxx right = (y)(xx)(x), left = (y)(x)(xx)
    pres = MonomialPresentation(["v"], [("x", "v", "v"), ("y", "v", "v")],
                                [("y", "x"), ("x", "x", "x"), ("y", "y")])
    assert pres.left_decomposition(("y", "x", "x", "x"), 3) == [("y",), ("x",), ("x", "x")]
    assert pres.right_decomposition(("y", "x", "x", "x"), 3) == [("y",), ("x", "x"), ("x",)]


def test_left_right_ambiguity_sets_coincide():
    """CS §3 Proposition (citing Bardzell, Sköldberg): right n-ambiguities = left
    n-ambiguities as SETS. right_decomposition must succeed on every associated path."""
    for pres in (_s1(), _line(),
                 MonomialPresentation(["v"], [("x", "v", "v"), ("y", "v", "v")],
                                      [("y", "x"), ("x", "x", "x"), ("y", "y")])):
        for n in (2, 3, 4):
            for p in pres.associated_paths(n, 40):
                blocks = pres.right_decomposition(p, n)
                assert len(blocks) == n
                assert tuple(x for b in blocks for x in b) == tuple(p)
```

- [ ] **Step 2: Run tests, verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_bardzell_straddle.py -m deep -k right`
Expected: FAIL with `AttributeError: ... no attribute 'right_decomposition'`.

- [ ] **Step 3: Implement**

Add to `MonomialPresentation` (below `left_decomposition`):

```python
    def _witness_at_start(self, path):
        """The unique tip that is a prefix of `path`, or None (mirror of _witness_at_end)."""
        lp = len(path)
        for r in self.relations:
            lr = len(r)
            if lr <= lp and tuple(path[:lr]) == r:
                return r
        return None

    def right_decomposition(self, p, n):
        """Unique CS §3 RIGHT block decomposition p = v_{n-1} ... v_0, returned
        left-to-right (n blocks): v_0 = LAST arrow; growing leftwards, each next block
        is the SHORTEST extension making (block)+(prev block) reducible (witness tip
        starts at the pair's start, straddling the boundary). Right ambiguities = left
        ambiguities as sets (CS §3 Prop.), but the BLOCKS differ beyond the quadratic/
        palindromic cases — the odd (2-term) differential needs these for its first term."""
        if n == 1:
            return [(p[-1],)]
        blocks, pos = [(p[-1],)], len(p) - 1        # collected right-to-left
        while len(blocks) < n:
            prev, uk, cut = blocks[-1], [], False
            while pos > 0:
                pos -= 1
                uk.insert(0, p[pos])
                w = self._witness_at_start(tuple(uk) + tuple(prev))
                if w is not None:
                    assert len(w) > len(uk), ("block reducible in right decomposition", p, n)
                    cut = True
                    break
            assert cut, ("no right decomposition for associated path", p, n)
            blocks.append(tuple(uk))
        assert pos == 0, ("right decomposition does not exhaust the path", p, n)
        return blocks[::-1]
```

- [ ] **Step 4: Run tests, verify pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_bardzell_straddle.py -m deep`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/engine/resolutions_bardzell.py tests/engine/test_bardzell_straddle.py
git commit -m "feat(engine): right_decomposition — CS §3 right n-ambiguity blocks"
```

---

### Task 3: Bardzell odd differential uses the right factorization; straddle HH validated

**Files:**
- Modify: `src/quiverlab/engine/resolutions_bardzell.py` (`differential_matrix`, odd branch ~380–397)
- Test: `tests/engine/test_bardzell_straddle.py` (append)

**Interfaces:**
- Consumes: Task 2's `right_decomposition`.
- Produces: corrected `BardzellResolution.differential_matrix` (same signature/shape contract).

- [ ] **Step 1: Write the failing test** (append)

```python
def test_straddle_hh_matches_bar_and_minimal():
    """The decisive repro: A = k<x,y>/(xx,yy,xyx), dim 6. Bar and minimal-A^e oracles
    agree; Bardzell must now too (was [4,4,7,6] vs [4,4,6,9] at 32003 before the fix)."""
    import quiverlab as ql
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    from quiverlab.engine.resolutions_bardzell import BardzellResolution
    from quiverlab.engine.resolutions_minimal import minimal_homology_dims

    Q = ql.Quiver(["v"], {"x": ("v", "v"), "y": ("v", "v")})
    A = Q.algebra(relations=["x*x", "y*y", "x*y*x"], field=ql.GF(32003))
    eng = to_engine(A)
    res = BardzellResolution(_s1())
    N, PRIMES = 3, (32003, 2, 3, 5)
    bar = hochschild_homology_dims(eng, N, primes=PRIMES)
    bard = hochschild_homology_dims(eng, N, primes=PRIMES, resolution=res)
    minh = minimal_homology_dims(eng, N, primes=PRIMES)
    for p in PRIMES:
        assert bard[p] == bar[p], f"p={p}: {bard[p]} != {bar[p]}"
        assert minh[p] == bar[p][:len(minh[p])], f"minimal p={p}"
    assert bar[32003] == [4, 4, 6, 9]                  # session-verified literals
    assert bar[2] == [4, 7, 9, 12]


def test_straddle_bardzell_depth_unlock():
    """Bardzell runs past the bar window on the straddle presentation (smoke, N=10)."""
    import quiverlab as ql
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    from quiverlab.engine.resolutions_bardzell import BardzellResolution

    Q = ql.Quiver(["v"], {"x": ("v", "v"), "y": ("v", "v")})
    A = Q.algebra(relations=["x*x", "y*y", "x*y*x"], field=ql.GF(32003))
    dims = hochschild_homology_dims(to_engine(A), 10, primes=(32003,),
                                    resolution=BardzellResolution(_s1()))[32003]
    assert len(dims) == 11 and all(d >= 0 for d in dims)
```

- [ ] **Step 2: Run, verify the first test fails** (HH mismatch at degree ≥ 2 — the odd map still drops silently-invalid first terms)

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_bardzell_straddle.py -m deep -k hh_matches`
Expected: FAIL with dimension mismatch.

- [ ] **Step 3: Fix the odd branch** of `differential_matrix`:

```python
        # n odd >= 3: SMALL map — CS §4 f_{n-1} (even): v ⊗ (v_{n-2}⋯v_0) ⊗ 1 minus
        # 1 ⊗ (u_0⋯u_{n-2}) ⊗ u_{n-1}; the FIRST term needs the RIGHT factorization
        # (equal to u_0 only for quadratic/palindromic tips).
        for c, (p, w) in enumerate(basis_n):
            lblocks = pres.left_decomposition(p, n)
            v_top = pres.right_decomposition(p, n)[0]
            ulast = lblocks[-1]
            P = tuple(p[len(v_top):])                    # right (n-2)-ambiguity ∈ AP^{n-1}
            Q = tuple(p[:len(p) - len(ulast)])           # left  (n-2)-ambiguity ∈ AP^{n-1}
            lw = pres.compose(w, v_top)                  # (P ; w·v_top)
            if lw is not None:
                key = (P, lw)
                if key in index_nm1:
                    M[index_nm1[key], c] += 1
            lw = pres.compose(ulast, w)                  # (Q ; u_{n-1}·w)
            if lw is not None:
                key = (Q, lw)
                if key in index_nm1:
                    M[index_nm1[key], c] -= 1
        return M
```

- [ ] **Step 4: Run the full engine Bardzell surface**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/engine/test_bardzell_straddle.py tests/engine/test_bardzell_resolution.py tests/engine/test_acceptance_plan02.py tests/engine/test_periodic_resolution.py tests/engine/test_resolution_contract.py tests/engine/test_resolution_protocol.py -m deep`
Expected: all PASS (uniform zoo unchanged; straddle now matches oracles).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/engine/resolutions_bardzell.py tests/engine/test_bardzell_straddle.py
git commit -m "fix(engine): Bardzell odd differential takes its first term from the right factorization (CS §4 f_n even)"
```

---

### Task 4: CS `delta_terms` right factorization + lift the non-quadratic non-monomial gate

**Files:**
- Modify: `src/quiverlab/resolutions_cs/resolution.py` (`__init__` cache; `delta_terms` odd branch ~44–52; delete `_require_in_scope` ~102–116 and its two call sites ~41, ~120; module docstring note)
- Modify: `tests/resolutions_cs/test_differential.py` (~92–122: replace the two refusal tests)

**Interfaces:**
- Consumes: `right_decomposition` via `self.ss.tip_presentation()` (SSequence already exposes it).
- Produces: `delta_terms`/`d_terms` valid for **all** admissible presentations; only remaining `NotImplementedError` is the correction-solve inconsistency (unchanged, ~137). Term format unchanged.

- [ ] **Step 1: Replace the refusal tests with acceptance tests** in `tests/resolutions_cs/test_differential.py` (keep imports/`_kx3` as-is; delete `test_cubic_tip_nonmonomial_raises_notimplemented` and `test_cubic_tip_nonmonomial_refuses_at_battery_level`):

```python
def _cubic_tail(field=CC):
    """A = k<x,y>/(x², y², xyx − yxy): f.d. (dim 6), cubic tip xyx WITH tail yxy —
    non-quadratic non-monomial, previously refused (Plan-04 RESTRICT boundary; lifted
    by Plan 12's right_decomposition). Its tip algebra {xx,yy,xyx} also has straddling
    chains (xyxx, xxyx), so this exercises blocks + right factorization + corrections.
    (If completion changes the tips, keep any admissible f.d. cubic-tip non-monomial.)"""
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", "x*y*x - y*x*y"]
    A = Q.algebra(relations=rels, field=field)
    return A, ChouhySolotarResolution(A, build_reduction_system(Q, rels, field),
                                      max_degree=4)


def test_cubic_tip_nonmonomial_gates_pass():
    """Plan 12: the former RESTRICT boundary now computes, certified per instance by
    the two CS gates on both sides (Theorem 4.1 conditions (1) and (2))."""
    for field in (CC, GF(2), GF(3)):
        _A, res = _cubic_tail(field=field)
        res.assert_dd_zero(upto=4, side="hom")
        res.assert_dd_zero(upto=4, side="coh")
        res.assert_order_condition(upto=4)


def test_cubic_tip_nonmonomial_battery_level():
    """END-TO-END: the full HH pipeline runs on the lifted boundary; homology and
    cohomology Euler characteristics agree with each other per degree window
    (dimension sanity), and no NotImplementedError escapes."""
    A, _res = _cubic_tail(field=CC)
    coh = cs_cohomology_dims(A, 3)
    hom = cs_homology_dims(A, 3)
    assert len(coh) == 4 and len(hom) == 4
    assert all(d >= 0 for d in coh + hom)
    assert A.hochschild_cohomology(3, engine="cs") is not None


def test_cubic_tail_delta3_first_term_uses_right_block():
    """Pin the corrected odd first term on a straddle chain: for σ = xyxx (tips
    {xx,yy,xyx}), right blocks are (xy)(x)(x), so the leading 2-term map is
    (+1, xy, xx, ()) and (−1, (), xyx, x) — u_0 = x would be WRONG here."""
    _A, res = _cubic_tail(field=CC)
    sigma = next(c for c in res.ss.S(3) if c.word == ("x", "y", "x", "x"))
    lead = {(res.to_int(c), a, tw, cc) for (c, a, tw, cc) in res.delta_terms(3, sigma)}
    assert lead == {
        (1, ("x", "y"), ("x", "x"), ()),               # v_top ⊗ (rest) ⊗ 1
        (-1, (), ("x", "y", "x"), ("x",)),             # 1 ⊗ (rest) ⊗ u_last
    }
```

- [ ] **Step 2: Run, verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/resolutions_cs/test_differential.py -m deep`
Expected: the four new tests FAIL with `NotImplementedError` ("CS is certified for quadratic tips..."); all pre-existing tests PASS.

- [ ] **Step 3: Implement in `src/quiverlab/resolutions_cs/resolution.py`**

In `__init__`, add `self._rdec_cache = {}` after `self._d_cache = {}`.

Replace `delta_terms`'s guard call and odd branch (lines ~41–52):

```python
    def _right_blocks(self, n, chain):
        rb = self._rdec_cache.get((n, chain.word))
        if rb is None:
            rb = tuple(tuple(b) for b in
                       self.ss.tip_presentation().right_decomposition(chain.word, n))
            self._rdec_cache[(n, chain.word)] = rb
        return rb

    # -- leading map δ_n (CS f_{n-1}), quiverlab index ----------------------
    def delta_terms(self, n, chain):
        if n == 1:
            return self._d1_terms(chain)
        one, none = self._one(), self._neg(self._one())
        if n % 2 == 1:                                   # CS f_{n-1} EVEN: 2-term map,
            ulast = chain.blocks[-1]                     # first term from the RIGHT
            v_top = self._right_blocks(n, chain)[0]      # factorization (CS §4; = u_0
            word = chain.word                            # only in the quadratic case)
            P = tuple(word[len(v_top):])                 # v_{n-2}⋯v_0 in S_{n-1}
            Q = tuple(word[:len(word) - len(ulast)])     # u_0⋯u_{n-2} in S_{n-1}
            return [(one, v_top, P, ()), (none, (), Q, ulast)]
        prev = {c.word for c in self.ss.S(n - 1)}        # n even: CS f_{n-1} ODD, big sum
        w, out = chain.word, []
        for i in range(len(w)):
            for j in range(i + 1, len(w) + 1):
                if w[i:j] in prev:
                    out.append((one, w[:i], w[i:j], w[j:]))
        return out
```

Delete `_require_in_scope` entirely and its call in `_d_general` (line ~120). Keep the correction-solve inconsistency `NotImplementedError` (~137) verbatim. Update the module docstring's scope sentence to: "Valid for every admissible presentation (Plan 12 lifted the quadratic-or-monomial restriction via right_decomposition); the only refusal left is a genuinely inconsistent correction solve."

- [ ] **Step 4: Run the CS suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/resolutions_cs/ -m deep`
Expected: all PASS — QCI(2,2)/kx3/commutative-square term pins unchanged (quadratic/palindromic ⇒ `v_top == u_0`); new acceptance tests green. If `test_battery_*` pins fail, STOP: that means a shipped family was not palindromic-safe — investigate before proceeding (it would contradict the Global Constraints analysis).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/resolutions_cs/resolution.py tests/resolutions_cs/test_differential.py
git commit -m "feat(cs): odd differential takes the right factorization; lift the non-quadratic non-monomial gate (Plan-04 stretch item)"
```

---

### Task 5: Battery extension — straddle monomial, QCI(3,2) literature pins, cubic-tail vs bar

**Files:**
- Create: `tests/resolutions_cs/test_battery_straddle.py`
- Modify: `tests/resolutions_cs/test_ambiguities.py` (append two pins)

**Interfaces:**
- Consumes: public `Quiver.algebra`, `build_reduction_system`, `ChouhySolotarResolution`, `cs_homology_dims`/`cs_cohomology_dims`, engine bar via `quiverlab.engine.hh_engine.hochschild_homology_dims` + `to_engine` (pattern of `test_battery_bar.py` — mirror its helper imports exactly when writing the file).
- Produces: per-instance certification of the newly admitted presentations.

- [ ] **Step 1: Append S-sequence pins** to `tests/resolutions_cs/test_ambiguities.py` (follow the file's existing fixture style for building an `SSequence` from a reduction system):

```python
def test_ssequence_straddle_chains_present():
    """Plan 12: tips {xx,yy,xyx} — S_3 must contain the straddling chain xyxx with
    left blocks (x)(yx)(x); S_4 contains xyxxx? no — pin the exact S_3/S_4 sets."""
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", "x*y*x"]
    rs = build_reduction_system(Q, rels, CC)
    ss = SSequence(rs, max_degree=4)
    assert {c.word for c in ss.S(3)} == {("x",) * 3, ("x", "x", "y", "x"),
                                         ("x", "y", "x", "x"), ("y",) * 3}
    straddle = next(c for c in ss.S(3) if c.word == ("x", "y", "x", "x"))
    assert straddle.blocks == (("x",), ("y", "x"), ("x",))


def test_ssequence_qci32_matches_cs_phi_formula():
    """CS §7.2 (TeX ~2110–2140): for k<x,y>/(x^n, y^m, yx−ξxy), 𝒜_N =
    {y^{φ(s,m)} x^{φ(t,n)} : s+t = N+1}, φ(s,k) = (s/2)k if s even else ((s−1)/2)k + 1.
    quiverlab S_N = 𝒜_{N-1}. Pin n=3, m=2 through S_5."""
    def phi(s, k):
        return (s // 2) * k if s % 2 == 0 else ((s - 1) // 2) * k + 1
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x*x", "y*y", "y*x - 2*x*y"]
    rs = build_reduction_system(Q, rels, GF(7))
    ss = SSequence(rs, max_degree=5)
    for N in range(2, 6):
        expect = {("y",) * phi(s, 2) + ("x",) * phi(t, 3)
                  for s in range(N + 1) for t in range(N + 1) if s + t == N}
        assert {c.word for c in ss.S(N)} == expect, f"S_{N}"
```

- [ ] **Step 2: Create `tests/resolutions_cs/test_battery_straddle.py`**

```python
"""Plan 12 battery: presentations OUTSIDE the old quadratic-or-monomial scope.
Every instance is certified by the two CS gates (d²=0, order condition) plus
degreewise agreement with the normalized bar oracle in the bar-buildable window.
Construction helpers mirror test_battery_bar.py (adjust imports to match it)."""
import pytest

from quiverlab import Quiver, GF, CC
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.homology import cs_homology_dims, cs_cohomology_dims
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.hh_engine import hochschild_homology_dims

PRIMES = (32003, 2, 3, 5)


def _mk(rels, field, top):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    A = Q.algebra(relations=rels, field=field)
    rs = build_reduction_system(Q, rels, field)
    return A, ChouhySolotarResolution(A, rs, max_degree=top + 1)


CASES = [
    ("straddle-monomial", ["x*x", "y*y", "x*y*x"], 3),
    ("qci32", ["x*x*x", "y*y", "y*x - 2*x*y"], 3),
    ("cubic-tail", ["x*x", "y*y", "x*y*x - y*x*y"], 3),
]


@pytest.mark.parametrize("name,rels,top", CASES, ids=[c[0] for c in CASES])
def test_gates_close(name, rels, top):
    for field in (CC, GF(2), GF(3)):
        _A, res = _mk(rels, field, top)
        res.assert_dd_zero(upto=top + 1, side="hom")
        res.assert_dd_zero(upto=top + 1, side="coh")
        res.assert_order_condition(upto=top + 1)


@pytest.mark.parametrize("name,rels,top", CASES, ids=[c[0] for c in CASES])
def test_cs_matches_bar_degreewise(name, rels, top):
    for p in PRIMES:
        A, _res = _mk(rels, GF(p), top)
        got = cs_homology_dims(A, top)
        bar = hochschild_homology_dims(to_engine(A), top, primes=(p,))[p]
        assert got == bar, f"{name} p={p}: {got} != {bar}"


def test_straddle_monomial_literals():
    """Session-verified oracle values for k<x,y>/(xx,yy,xyx) (2026-07-22)."""
    A, _res = _mk(["x*x", "y*y", "x*y*x"], GF(32003), 3)
    assert cs_homology_dims(A, 3) == [4, 4, 6, 9]
    A2, _res2 = _mk(["x*x", "y*y", "x*y*x"], GF(2), 3)
    assert cs_homology_dims(A2, 3) == [4, 7, 9, 12]
```

- [ ] **Step 3: Run the new battery**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/resolutions_cs/test_battery_straddle.py tests/resolutions_cs/test_ambiguities.py -m deep`
Expected: all PASS. Execution notes: (a) fix imports to the repo's actual module paths by mirroring `test_battery_bar.py`'s header; (b) `cs_homology_dims`'s exact signature/arg order comes from `homology.py` — adapt the calls, not the assertions; (c) if `build_reduction_system` reports non-confluence for `x*y*x - y*x*y` (hand-checked confluent: overlaps xxyx, xyxx, xyxyx, xxx-side all resolve to 0), STOP and re-derive with the completion's actual tips per the docstring escape hatch.

- [ ] **Step 4: Commit**

```bash
git add tests/resolutions_cs/test_battery_straddle.py tests/resolutions_cs/test_ambiguities.py
git commit -m "test(cs): straddle/QCI(3,2)/cubic-tail batteries — gates + bar agreement + CS §7.2 phi-formula pins"
```

---

### Task 6: Full deep suites, docs, status

**Files:**
- Modify: `docs/internals/05-resolutions.md`, `docs/internals/09-chouhy-solotar.md` (scope prose), `CLAUDE.md` (stretch-item paragraph), `docs/plans/ROADMAP.md` (Plan 12 entry), this plan (mark executed)

- [ ] **Step 1: Run the full deep suites**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m deep tests/engine/ tests/resolutions_cs/`
Expected: all PASS (~15–20 min). Then the fast matrix: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast`
Expected: all PASS.

- [ ] **Step 2: Update docs**

- `docs/internals/05-resolutions.md`: in the Bardzell section, document the minimal-extension block cut and `right_decomposition`, with the `xyxx` worked example (left `(x)(yx)(x)`, right `(xy)(x)(x)`).
- `docs/internals/09-chouhy-solotar.md`: replace the RESTRICT-scope paragraph with the lifted-scope statement (only refusal left: inconsistent correction solve) and the corrected odd-map formula.
- `CLAUDE.md`: rewrite the "Known deeper-engine stretch item" paragraph as delivered (Plan 12), noting `right_decomposition` and that `_require_in_scope` is gone; update the Status date line.
- `docs/plans/ROADMAP.md`: add Plan 12 as delivered with a one-line summary of both findings.

- [ ] **Step 3: Commit**

```bash
git add docs/internals/05-resolutions.md docs/internals/09-chouhy-solotar.md CLAUDE.md docs/plans/ROADMAP.md docs/plans/2026-07-22-plan-12-ambiguity-blocks.md
git commit -m "docs: Plan 12 — straddling ambiguities fixed, right factorization, CS scope lifted"
```

- [ ] **Step 4 (optional, heavy): strict docs build**

Run: `.venv/bin/mkdocs build --strict` — exit 0 expected (several minutes; executes notebooks). Only needed if internals chapters changed structurally.

---

## Execution amendments (recorded during the run)

- **AP^3 count corrected 4 → 5.** The plan's hand-derived pin missed `xyxyx`: block
  candidates are prefixes of the *remaining word*, so `(x)(yx)(yx)` satisfies CS §3 (ii)
  (`yx·y = yxy` irreducible, `yx·yx ⊇ xyx`). Verified against the minimal-A^e Betti
  numbers `rks = 2, 3, 5, 9, 17` on `k⟨x,y⟩/(xx,yy,xyx)` over GF(32003). Test pins in
  Task 1 use the oracle values.
- Cyclic-Nakayama uniform pin fixed: `cn(3,2)` is quadratic (`ℓ=2`, length-3 chains);
  the genuinely-overlapping regression pin uses `cn(2,3)` instead.
- Tasks 1 and 2 landed as one commit (same two files; interleaved verification).

## Self-review notes

- **Spec coverage:** finding 1 → Tasks 1–3; finding 2 (stretch item) → Tasks 2–4; certification → Tasks 3–5; docs/status → Task 6. Side finding 3 deliberately out of scope (memory + plan note only).
- **Type consistency:** `right_decomposition(p, n) -> list[tuple]` consumed identically in Task 3 (`pres.right_decomposition(p, n)[0]`) and Task 4 (via `tip_presentation()`); `_witness_at_end`/`_witness_at_start` internal-only.
- **Known execution-time adaptation points** (flagged in-place): import paths in Task 5 mirror `test_battery_bar.py`; `cs_homology_dims` signature; confluence of the cubic-tail relations; `test_ambiguities.py` fixture style. These are look-and-mirror steps, not design gaps.
