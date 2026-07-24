# Plan 19 — Field generality of engine-backed invariants Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** the five engine-backed invariants (`cyclic_homology`, `complexity`,
`is_frobenius`, `nakayama_automorphism`, `is_symmetric`) compute over EVERY exact
Domain (QQ, CC, QQi, GF(p^n)) — not just GF(p) — and the one residual refusal left
(structure-constants algebras without a quiver, off GF(p), for the invariants that
need a path-type basis) is reworded honestly: no more "later phase that
generalizes this invariant" promise (backlog Tier-1 item 6).

**Architecture:** GF(p) keeps routing to the engine byte-unchanged (int64 fast
rank). Off GF(p), three generic exact paths are added, each with a rigorous
correctness gate rather than trust:
(1) **Cyclic homology** — Connes' B on the normalized bar basis over any Domain
(`hochschild/cyclic.py`, mirroring `engine/cyclic.py`'s conventions exactly), HC
via the (b,B) total complex; gated by GF(p) parity with the engine, the
mixed-complex identities over QQ, and a *second chain model* (Connes' λ-complex,
valid over char 0 — Loday Thm 2.1.5) implemented as a test oracle.
(2) **Complexity** — Betti numbers of the minimal A^e resolution over any Domain
via the E-relative (Cibils) complex `r^{⊗_E •}` with the middle-face differential
(`invariants/betti.py`): for a path-type basis A = E ⊕ r, dim H_n of this small
complex = the engine's `rks[n]` generator counts, over every field (minimality
kills the induced differential on E⊗E-coefficients); gated by exact GF(p) parity
with `minimal_resolution` incl. multi-vertex + straddling zoo records, and char-0
closed forms.
(3) **Frobenius / Nakayama / symmetry** — the exact socle criterion for basic
split algebras (Nakayama; Skowroński–Yamagata *Frobenius Algebras I*): A is
Frobenius ⟺ every soc(e_vA) is simple and v ↦ vertex(soc(e_vA)) is a permutation
(pure Domain linear algebra, conclusive both ways); the form is the socle-dual
functional, *verified* nondegenerate (Gram rank); ν = G⁻¹Gᵀ, self-certified in
tests by the defining identity λ(ab) = λ(b·ν(a)) + multiplicativity; symmetry =
Frobenius ∧ ν inner, decided by (a) nontrivial Nakayama vertex permutation → NOT
symmetric (inner automorphisms fix primitive-idempotent classes), else (b) an
exact Schwartz–Zippel grid sweep for an invertible element of the twisted
centralizer U = {u : ν(a)u = ua} — conclusive whenever the Domain supplies
> dim A distinct samples (always, off small finite fields), LOUD
inconclusive-error otherwise (never a silent wrong answer).
The shared path-type-basis split (idempotent labels `e_v`, radical source/target
via multiplication — no label parsing) lives in `invariants/pathbasis.py` and
raises the reworded FieldError when the basis is not path-type.

**Tech Stack:** `fields.linalg` (rank/nullspace/solve over any Domain),
`hochschild/bar.py` machinery (normalized-bar basis + `boundary_matrix` already
generic), `engine/*` untouched (parity oracles only), `families` + Plan-18 zoo
records as the test zoo. Tests in `tests/invariants/` → fast bucket (sizes kept
tiny; heaviest single rank ≈ 364² over GF(5)).

## Global Constraints

- No float literals in `src/` (AST gate `tests/test_no_floats.py`).
- Python is always `.venv/bin/python`; tests via
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q …`.
- Exact only: all new arithmetic through `Domain` ops / `fields.linalg`.
- GF(p) behavior byte-unchanged: `isinstance(domain, PrimeField)` keeps the exact
  current engine code paths for all five invariants.
- Engine internals stay internal: generic paths import NOTHING from
  `quiverlab.engine.*` except in the GF(p) dispatch branches that already do, and
  `engine.scan3.complexity_of` (pure growth-labeling helper, no field content).
- Oracles live, never hardcoded from memory: engine parity over GF(p); the
  λ-complex second model over QQ; closed forms only where classical and
  char-independent (semisimple HC, monomial Betti, hereditary termination).
- Char-2 trap: `QuantumCI(q)` degenerates when q ≡ 0 (e.g. q=2 over GF(4) gives
  the monomial k⟨x,y⟩/(x²,y²,xy), socle 2-dim, NOT Frobenius) — per-field
  expectations in tests are explicit, never "same value for every field".
- Conventional commits; green at every commit; merge/push only when Marco asks.
- Branch `plan-19-field-generality` stacks on `plan-18-zoo-diversity` (unmerged;
  the Betti/Frobenius zoo reuses Plan-18 catalog records `comm_square`, `cn_3_2`,
  `straddle_xx_yy_xyx`).

---

### Task 1: Generic cyclic homology (`hochschild/cyclic.py`)

**Files:**
- Create: `src/quiverlab/hochschild/cyclic.py`
- Modify: `src/quiverlab/core/algebra.py` (`cyclic_homology` dispatch; keep the
  GF(p) body verbatim)
- Test: `tests/invariants/test_cyclic_generic.py` (new)
- Modify: `tests/invariants/test_engine_backed.py` (flip `test_cyclic_homology_cc_loud`)

**Interfaces:**
- Consumes: `hochschild.bar._cochain_basis(m, n)`, `hochschild.bar.boundary_matrix(A, n, max_cells)`
  (returns `(rows_list_of_lists, ncols, nrows)`, rows = C_{n−1} basis),
  `fields.linalg.rank(rows, dom)`, `HHTable(dims, kind, algebra_repr, engine=…)`.
- Produces: `connes_B_matrix(A, n, max_cells) -> (M, ncols, nrows)` (M over
  A.domain, rows = C_{n+1} basis, unit-adapted A) and
  `cyclic_homology_dims(A, top, max_cells=4_000_000) -> HHTable` (kind `"HC_"`;
  unit-adapts internally). Task 4's structure-constants test relies on this NOT
  requiring a quiver.

- [ ] **Step 1: Write the failing tests** — create
  `tests/invariants/test_cyclic_generic.py`:

```python
"""Generic-Domain cyclic homology (Plan 19): GF(p) parity with the engine,
mixed-complex identities over QQ, the char-0 Connes lambda-complex second
model (Loday Thm 2.1.5), and closed-form pins."""
import itertools

from quiverlab import CC, GF, Quiver, truncated_polynomial
from quiverlab.families import QuantumCI
from quiverlab.fields import QQ, linalg


def _k(field):
    return Quiver([1], {}).algebra(relations=[], field=field)


def _kxk(field):
    return Quiver([1, 2], {}).algebra(relations=[], field=field)


def test_hc_generic_matches_engine_over_gfp():
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.cyclic import cyclic_homology_dims as engine_hc
    from quiverlab.hochschild.cyclic import cyclic_homology_dims as generic_hc
    cases = ((truncated_polynomial(2, field=GF(5)), 3),
             (truncated_polynomial(3, field=GF(3)), 3),
             (QuantumCI(2, field=GF(5)), 2))
    for A, top in cases:
        p = A.domain.p
        want = engine_hc(to_engine(A.unit_adapted()), top, primes=(p,))[p]
        got = generic_hc(A, top)
        assert got.kind == "HC_"
        assert got.dims == [int(d) for d in want]


def _matmul(X, Y, dom):
    """(a x b) @ (b x c) over dom, list-of-rows convention."""
    if not X or not Y:
        return []
    c = len(Y[0])
    out = [[dom.zero()] * c for _ in range(len(X))]
    for i, Xi in enumerate(X):
        for j, xv in enumerate(Xi):
            if dom.is_zero(xv):
                continue
            Yj = Y[j]
            for k in range(c):
                if not dom.is_zero(Yj[k]):
                    out[i][k] = dom.add(out[i][k], dom.mul(xv, Yj[k]))
    return out


def test_mixed_complex_identities_over_qq():
    from quiverlab.hochschild.bar import boundary_matrix
    from quiverlab.hochschild.cyclic import connes_B_matrix
    A = truncated_polynomial(2, field=QQ).unit_adapted()
    dom = A.domain
    b = {k: boundary_matrix(A, k, 10 ** 6)[0] for k in range(1, 5)}
    B = {k: connes_B_matrix(A, k, 10 ** 6)[0] for k in range(0, 4)}

    def zero(M):
        return all(dom.is_zero(x) for row in M for x in row)

    for n in range(0, 3):
        if n >= 1:
            assert zero(_matmul(b[n], b[n + 1], dom)), f"b^2 != 0 at {n}"
        assert zero(_matmul(B[n + 1], B[n], dom)), f"B^2 != 0 at {n}"
        bB = _matmul(b[n + 1], B[n], dom)
        if n >= 1:
            Bb = _matmul(B[n - 1], b[n], dom)
            S = [[dom.add(bB[i][j], Bb[i][j]) for j in range(len(bB[0]))]
                 for i in range(len(bB))]
        else:
            S = bB
        assert zero(S), f"bB + Bb != 0 at {n}"


def _lambda_hc(A, top):
    """Char-0 SECOND MODEL: HC_n = H_n(C^lambda) where C^lambda_n =
    A^{(x)(n+1)} / im(1 - t), t = (-1)^n cyclic rotation, b = the FULL cyclic
    Hochschild boundary (n+1 faces, wrap-around last).  Valid over fields
    containing Q (Loday, Cyclic Homology, Thm 2.1.5).  Independent of
    hochschild/bar.py AND hochschild/cyclic.py: unnormalized chains, no
    unit-adaptation, quotient model instead of bicomplex."""
    dom = A.domain
    m = A.dim

    def b_cols(n):
        """Columns of b: C_n -> C_{n-1} as vectors (rank is transpose-invariant)."""
        colws = list(itertools.product(range(m), repeat=n + 1))
        rowws = list(itertools.product(range(m), repeat=n))
        ridx = {w: i for i, w in enumerate(rowws)}
        cols = []
        for w in colws:
            v = [dom.zero()] * len(rowws)
            for i in range(n + 1):
                if i < n:
                    prod, key_of = A.T[w[i]][w[i + 1]], (lambda t, i=i: w[:i] + (t,) + w[i + 2:])
                else:
                    prod, key_of = A.T[w[n]][w[0]], (lambda t: (t,) + w[1:n])
                for t, cf in enumerate(prod):
                    if dom.is_zero(cf):
                        continue
                    val = cf if i % 2 == 0 else dom.neg(cf)
                    r = ridx[key_of(t)]
                    v[r] = dom.add(v[r], val)
            cols.append(v)
        return cols

    def v_gens(n):
        """Spanning set of im(1 - t) on C_n."""
        ws = list(itertools.product(range(m), repeat=n + 1))
        idx = {w: i for i, w in enumerate(ws)}
        out = []
        for w in ws:
            rot = (w[-1],) + w[:-1]
            v = [dom.zero()] * len(ws)
            v[idx[w]] = dom.add(v[idx[w]], dom.one())
            s = dom.neg(dom.one()) if n % 2 == 0 else dom.one()   # -(-1)^n
            v[idx[rot]] = dom.add(v[idx[rot]], s)
            out.append(v)
        return out

    qdim, rkbar = {}, {0: 0}
    for n in range(top + 2):
        vg = v_gens(n)
        qdim[n] = m ** (n + 1) - (linalg.rank(vg, dom) if vg else 0)
        if n >= 1:
            vprev = v_gens(n - 1)
            rprev = linalg.rank(vprev, dom) if vprev else 0
            rkbar[n] = linalg.rank(b_cols(n) + vprev, dom) - rprev
    return [qdim[n] - rkbar[n] - rkbar[n + 1] for n in range(top + 1)]


def test_hc_lambda_complex_second_model_qq():
    for A, top in ((_k(QQ), 4), (_kxk(QQ), 3),
                   (truncated_polynomial(2, field=QQ), 3)):
        assert A.cyclic_homology(top).dims == _lambda_hc(A, top)


def test_hc_closed_forms_char0_and_gf4():
    for field in (QQ, GF(4)):
        assert _k(field).cyclic_homology(4).dims == [1, 0, 1, 0, 1]
        assert _kxk(field).cyclic_homology(4).dims == [2, 0, 2, 0, 2]


def test_hc_dual_numbers_cc_matches_qq():
    # dims are invariant under field extension QQ -> CC (flat base change)
    got = truncated_polynomial(2, field=CC).cyclic_homology(2)
    assert got.kind == "HC_"
    assert got.dims == truncated_polynomial(2, field=QQ).cyclic_homology(2).dims
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/invariants/test_cyclic_generic.py -p no:cacheprovider`
Expected: FAIL / ERROR — `ModuleNotFoundError: quiverlab.hochschild.cyclic`, and
the dispatch tests raise FieldError.

- [ ] **Step 3: Implement `src/quiverlab/hochschild/cyclic.py`**

```python
"""Cyclic homology over any exact Domain (Plan 19): Connes' B on the
normalized bar complex + the (b, B) bicomplex — the generic mirror of
engine/cyclic.py (which stays authoritative over GF(p) via the fast rank).

Conventions (identical to the engine; parity is gated in
tests/invariants/test_cyclic_generic.py):

  * C_n = A (x) Abar^{(x)n} on the unit-adapted basis: basis (s, J) with
    s in 0..m-1, J in {1..m-1}^n (hochschild/bar.py's shapes; basis 0 = 1_A).
  * b = the normalized bar boundary (hochschild/bar.py::boundary_matrix).
  * B(a_0 (x) a_1 (x) ... (x) a_n) =
        sum_{i=0}^{n} (-1)^{n i} 1 (x) a_i (x) ... (x) a_n (x) a_0 (x) ... (x) a_{i-1}
    — insert 1_A in the A-slot, rotate all n+1 entries; any rotation that
    puts the unit (basis index 0) in a bar slot dies in Abar = A/k.1.
  * HC_n = dim Tot_n - rank D_n - rank D_{n+1},
    Tot_n = C_n (+) C_{n-2} (+) ..., D = b + B.
"""
from quiverlab.errors import DepthLimitError
from quiverlab.fields.linalg import rank
from quiverlab.hochschild.bar import _cochain_basis, boundary_matrix
from quiverlab.hochschild.table import HHTable

_GUARD_HINT = ("the (b, B) bicomplex on the bar basis is exponential; over "
               "GF(p) use the fast engine — raise max_cells only if you know "
               "what you are doing")


def connes_B_matrix(A, n, max_cells):
    """Matrix of B : C_n -> C_{n+1} over A.domain for unit-adapted A.
    Returns (list-of-rows, ncols, nrows); rows indexed by the C_{n+1} basis."""
    dom = A.domain
    m = A.dim
    cols = _cochain_basis(m, n)
    rows = _cochain_basis(m, n + 1)
    if len(rows) * len(cols) > max_cells:
        raise DepthLimitError(
            f"Connes B_{n}: matrix would have {len(rows)} x {len(cols)} "
            f"entries (> max_cells = {max_cells})", hint=_GUARD_HINT)
    row_index = {b: i for i, b in enumerate(rows)}
    M = [[dom.zero()] * len(cols) for _ in range(len(rows))]
    one, negone = dom.one(), dom.neg(dom.one())
    for ci, (s, J) in enumerate(cols):
        entries = (s,) + J                       # (a_0, r_1..r_n), length n+1
        for i in range(n + 1):
            rotated = entries[i:] + entries[:i]
            if 0 in rotated:                     # the unit died in a bar slot
                continue
            r = row_index[(0, rotated)]          # new A-slot = the unit
            sign = one if (n * i) % 2 == 0 else negone
            M[r][ci] = dom.add(M[r][ci], sign)
    return M, len(cols), len(rows)


def _tot_degrees(n):
    """Chain degrees in Tot_n = C_n (+) C_{n-2} (+) ... (descending)."""
    return list(range(n, -1, -2))


def cyclic_homology_dims(A, top, max_cells=4_000_000):
    """HHTable of dim HC_0..HC_top over A.domain via the (b, B) bicomplex.

    Works for ANY unital algebra (no quiver needed): only the unit-adapted
    basis is required. Exponential in top like the bar oracle; max_cells
    guards every assembled matrix."""
    B0 = A.unit_adapted()
    dom = B0.domain
    m = B0.dim
    maxdeg = top + 1
    dims = {k: m * (m - 1) ** k for k in range(maxdeg + 2)}
    bmats = {k: boundary_matrix(B0, k, max_cells)[0] for k in range(1, maxdeg + 1)}
    Bmats = {k: connes_B_matrix(B0, k, max_cells)[0] for k in range(0, maxdeg)}
    ranks = {}
    for n in range(top + 2):
        src, tgt = _tot_degrees(n), _tot_degrees(n - 1)
        if not tgt:
            ranks[n] = 0
            continue
        row_off, off = {}, 0
        for d in tgt:
            row_off[d] = off
            off += dims[d]
        nrows = off
        ncols = sum(dims[d] for d in src)
        if nrows * ncols > max_cells:
            raise DepthLimitError(
                f"cyclic total differential D_{n}: {nrows} x {ncols} entries "
                f"(> max_cells = {max_cells})", hint=_GUARD_HINT)
        D = [[dom.zero()] * ncols for _ in range(nrows)]
        c0 = 0
        for d in src:
            for r0_key, blk in ((d - 1, bmats.get(d) if d >= 1 else None),
                                (d + 1, Bmats.get(d))):
                if blk is None or r0_key not in row_off:
                    continue
                r0 = row_off[r0_key]
                for r, rowvec in enumerate(blk):
                    Dr = D[r0 + r]
                    for c, val in enumerate(rowvec):
                        if not dom.is_zero(val):
                            Dr[c0 + c] = dom.add(Dr[c0 + c], val)
            c0 += dims[d]
        ranks[n] = rank(D, dom) if nrows and ncols else 0
    out = []
    for n in range(top + 1):
        tot = sum(dims[d] for d in _tot_degrees(n))
        out.append(tot - ranks[n] - ranks[n + 1])
    return HHTable(out, "HC_", repr(A).splitlines()[0],
                   engine=f"bar (b,B) mixed complex over {dom.name}")
```

- [ ] **Step 4: Dispatch** — in `src/quiverlab/core/algebra.py`, replace the body
  of `cyclic_homology` (keep the GF(p) branch verbatim; drop the
  `_require_prime_field` call — the method loses ALL refusals):

```python
    def cyclic_homology(self, top, max_cells=4_000_000):
        """Dimensions of HC_0..HC_top (Connes (b, B) mixed complex).

        GF(p): the fast engine (int64 rank). Any other exact Domain: the
        generic mixed complex on the normalized bar basis (exponential —
        max_cells guards the blow-up). Works for any unital algebra."""
        from quiverlab.fields.primefield import PrimeField
        if isinstance(self.domain, PrimeField):
            from quiverlab.engine.adapter import to_engine
            from quiverlab.engine.cyclic import cyclic_homology_dims
            from quiverlab.hochschild.table import HHTable
            p = self.domain.p
            out = cyclic_homology_dims(to_engine(self.unit_adapted()), top, primes=(p,))
            dims = [int(d) for d in out[p]]
            return HHTable(dims, "HC_", repr(self).splitlines()[0],
                           engine="hanlab engine (F_p fast rank)")
        from quiverlab.hochschild.cyclic import cyclic_homology_dims
        return cyclic_homology_dims(self, top, max_cells=max_cells)
```

- [ ] **Step 5: Flip the loud pin** — in `tests/invariants/test_engine_backed.py`
  replace `test_cyclic_homology_cc_loud` with:

```python
def test_cyclic_homology_cc_computes():
    # Plan 19: off GF(p) the generic (b, B) mixed complex serves any Domain
    t = truncated_polynomial(2, field=CC).cyclic_homology(2)
    assert t.kind == "HC_"
    assert len(t.dims) == 3 and all(isinstance(d, int) for d in t.dims)
```

- [ ] **Step 6: Run to verify pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/invariants/test_cyclic_generic.py tests/invariants/test_engine_backed.py -p no:cacheprovider`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/quiverlab/hochschild/cyclic.py src/quiverlab/core/algebra.py tests/invariants/test_cyclic_generic.py tests/invariants/test_engine_backed.py
git commit -m "feat(hochschild): cyclic homology over any exact Domain -- generic (b,B) mixed complex (Plan 19)"
```

### Task 2: Path-type basis split + generic Betti numbers / complexity

**Files:**
- Create: `src/quiverlab/invariants/pathbasis.py`
- Create: `src/quiverlab/invariants/betti.py`
- Modify: `src/quiverlab/invariants/scalar.py` (`complexity` routes)
- Test: `tests/invariants/test_betti_generic.py` (new)
- Modify: `tests/invariants/test_scalar.py` (flip `test_complexity_cc_loud`)

**Interfaces:**
- Consumes: `Algebra.T` / `.basis_labels` / `.quiver` / `.unit` / `.domain`,
  `fields.linalg.rank/nullspace`, `engine.scan3.complexity_of(seq)` (plain int
  list in, label out), Plan-18 `families.zoo.load_catalog/build_from_record`.
- Produces: `path_type_basis(A, what) -> (idem, rad, src, tgt)` (index lists +
  per-radical-index source/target idempotent index; raises the REWORDED
  FieldError — Task 3 reuses this), and
  `relative_betti_numbers(A, top, max_cells=4_000_000) -> list[int]`
  (= engine `rks[0..top]`).

- [ ] **Step 1: Write the failing tests** — create
  `tests/invariants/test_betti_generic.py`:

```python
"""Generic-Domain Betti numbers of the minimal A^e resolution (Plan 19):
exact GF(p) parity with engine/resolutions_minimal (incl. multi-vertex and
straddling Plan-18 records), char-0/GF(4) closed forms, complexity routing."""
import pytest

from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.families import QuantumCI
from quiverlab.families.zoo import build_from_record, load_catalog
from quiverlab.fields import QQ


def _rec(name):
    return next(r for r in load_catalog() if r.get("name") == name)


def test_betti_matches_engine_rks_over_gfp():
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import minimal_resolution
    from quiverlab.invariants.betti import relative_betti_numbers
    zoo = [truncated_polynomial(3, field=GF(5)),
           QuantumCI(2, field=GF(5)),
           linear_path_algebra(2, field=GF(5)),
           build_from_record(_rec("straddle_xx_yy_xyx"), field=GF(5)),
           build_from_record(_rec("comm_square"), field=GF(5)),
           build_from_record(_rec("cn_3_2"), field=GF(3))]
    for A in zoo:
        p = A.domain.p
        rks, _cols, _e, _trunc = minimal_resolution(to_engine(A), 5, p)
        want = [max(0, rks[k]) for k in sorted(rks)]
        got = relative_betti_numbers(A, 5)
        assert got == want, f"{getattr(A, 'zoo_name', repr(A))}: {got} != {want}"


def test_betti_char0_closed_forms():
    from quiverlab.invariants.betti import relative_betti_numbers
    # monomial: char-independent closed forms
    assert relative_betti_numbers(truncated_polynomial(3, field=QQ), 5) == [1] * 6
    assert relative_betti_numbers(linear_path_algebra(2, field=QQ), 4) == [2, 1, 0, 0, 0]
    # commutative complete intersection k[x,y]/(x^2, y^2): Betti_n = n + 1
    assert relative_betti_numbers(QuantumCI(-1, field=QQ), 5) == [1, 2, 3, 4, 5, 6]
    # gldim-2 multi-vertex: (vertices, arrows, relations, 0, 0, ...)
    assert relative_betti_numbers(
        build_from_record(_rec("comm_square"), field=QQ), 5) == [4, 4, 1, 0, 0, 0]


def test_betti_dd_zero_over_qq():
    # the middle-face differential is a complex: d o d = 0, asserted exactly
    from quiverlab.invariants.betti import _relative_complex
    A = QuantumCI(-1, field=QQ)
    dom = A.domain
    chains, mats = _relative_complex(A, 4, 4_000_000)
    for n in range(3, 5):
        M1, M2 = mats[n], mats[n - 1]          # d_n, d_{n-1}
        for ci in range(len(chains[n])):
            col = [M1[r][ci] for r in range(len(chains[n - 1]))]
            acc = [dom.zero()] * len(chains[n - 2])
            for j, cv in enumerate(col):
                if dom.is_zero(cv):
                    continue
                for r in range(len(chains[n - 2])):
                    acc[r] = dom.add(acc[r], dom.mul(cv, M2[r][j]))
            assert all(dom.is_zero(x) for x in acc), f"d o d != 0 at {n}"


def test_complexity_field_generality():
    assert truncated_polynomial(2, field=CC).complexity(4) == 1
    assert truncated_polynomial(2, field=GF(4)).complexity(4) == 1
    assert QuantumCI(-1, field=QQ).complexity(5) == 2
    assert linear_path_algebra(2, field=QQ).complexity(4) == 0
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/invariants/test_betti_generic.py -p no:cacheprovider`
Expected: FAIL — `ModuleNotFoundError: quiverlab.invariants.betti`; the
complexity calls raise FieldError.

- [ ] **Step 3: Implement `src/quiverlab/invariants/pathbasis.py`**

```python
"""Path-type basis extraction: the split A = E (+) r read off basis_labels.

A *path-type basis* (every Quiver(...).algebra(...) instance) has labels
'e_v' for a complete set of orthogonal idempotents spanning a Wedderburn
complement E, and radical labels spanning r = rad A, each with a unique
source and target idempotent. This is the structural hypothesis behind the
generic any-Domain engine-backed invariants (Plan 19); it is VERIFIED here
— structurally, via multiplication, never by parsing path labels — and the
refusal for non-path-type bases is the honest one (no 'later phase')."""
from quiverlab.errors import FieldError

_HINT = ("present the algebra via Quiver(...).algebra(...) (path-type basis), "
         "or construct it over a prime field GF(p), where the fast engine "
         "serves any basis")


def path_type_basis(A, what="this invariant"):
    """(idem, rad, src, tgt): idempotent and radical basis-index lists plus,
    per radical index i, the unique v, w in idem with e_v f_i = f_i = f_i e_w.
    Raises FieldError if the basis is not path-type."""
    dom = A.domain
    labels = A.basis_labels
    if A.quiver is None or not labels:
        raise FieldError(
            f"{what} over {dom.name} needs a quiver presentation "
            f"(path-type basis); this algebra has none", hint=_HINT)

    def is_basis_vec(vec, j):
        return all(dom.eq(c, dom.one() if t == j else dom.zero())
                   for t, c in enumerate(vec))

    idem = [i for i, lab in enumerate(labels) if lab.startswith("e_")]
    rad = [i for i, lab in enumerate(labels) if not lab.startswith("e_")]
    bad = FieldError(
        f"{what} over {dom.name}: basis is not path-type "
        f"(idempotent/vertex structure not visible in the labels)", hint=_HINT)
    # complete orthogonal idempotents summing to 1
    for v in idem:
        for w in idem:
            prod = A.T[v][w]
            if v == w:
                if not is_basis_vec(prod, v):
                    raise bad
            elif any(not dom.is_zero(c) for c in prod):
                raise bad
    for t, c in enumerate(A.unit):
        want = dom.one() if t in idem else dom.zero()
        if not dom.eq(c, want):
            raise bad
    src, tgt = {}, {}
    for i in rad:
        sv = [v for v in idem if is_basis_vec(A.T[v][i], i)]
        tw = [w for w in idem if is_basis_vec(A.T[i][w], i)]
        if len(sv) != 1 or len(tw) != 1:
            raise bad
        src[i], tgt[i] = sv[0], tw[0]
    return idem, rad, src, tgt
```

- [ ] **Step 4: Implement `src/quiverlab/invariants/betti.py`**

```python
"""Betti numbers of the minimal A^e resolution over ANY exact Domain
(Plan 19), via the E-relative (Cibils) complex.

For A = E (+) r split basic (path-type basis, E separable), the relative
bar resolution A (x)_E r^{(x)_E n} (x)_E A is an A^e-projective resolution
of A; applying (x)_{A^e} (E (x) E) kills the two outer faces (their factor
acts through A/r on E), leaving the SMALL complex

    T_n = r^{(x)_E n},   d(r_1 (x)...(x) r_n)
        = sum_{i=1}^{n-1} (-1)^i r_1 (x)..(x) r_i r_{i+1} (x)..(x) r_n

whose H_n = Tor_n^{A^e}(A, E (x) E). Minimality of the true minimal
resolution kills ITS induced differential on E (x) E coefficients, so
dim H_n = the number of A^e-generators of P_n — the engine's rks[n]
(engine/resolutions_minimal.py) — over every field. GF(p) parity is gated
in tests/invariants/test_betti_generic.py."""
from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.fields.linalg import rank
from quiverlab.invariants.pathbasis import path_type_basis

_GUARD_HINT = ("the relative-Tor chain count grows with the algebra's own "
               "complexity; raise max_cells only if you know what you are doing")


def _relative_complex(A, top, max_cells):
    """(chains, mats): chains[n] = composable radical index tuples (n >= 1;
    chains[0] is the list of idempotent indices standing for E), mats[n] =
    the matrix of d_n : T_n -> T_{n-1} for n >= 2 (d_1 = 0)."""
    dom = A.domain
    idem, rad, src, tgt = path_type_basis(A, "complexity")
    idemset = set(idem)
    chains = {0: [(v,) for v in idem], 1: [(i,) for i in rad]}
    for n in range(2, top + 2):
        chains[n] = [ch + (j,) for ch in chains[n - 1]
                     for j in rad if src[j] == tgt[ch[-1]]]
    mats = {}
    for n in range(2, top + 2):
        rows, cols = len(chains[n - 1]), len(chains[n])
        if rows * cols > max_cells:
            raise DepthLimitError(
                f"relative Tor d_{n}: {rows} x {cols} entries "
                f"(> max_cells = {max_cells})", hint=_GUARD_HINT)
        ridx = {ch: k for k, ch in enumerate(chains[n - 1])}
        M = [[dom.zero()] * cols for _ in range(rows)]
        for ci, ch in enumerate(chains[n]):
            for i in range(1, n):
                prod = A.T[ch[i - 1]][ch[i]]
                for t, cf in enumerate(prod):
                    if dom.is_zero(cf):
                        continue
                    if t in idemset:
                        raise QuiverlabError(
                            "product of radical basis elements has an "
                            "idempotent coordinate — basis is not path-type")
                    r = ridx[ch[:i - 1] + (t,) + ch[i + 1:]]
                    val = cf if i % 2 == 0 else dom.neg(cf)
                    M[r][ci] = dom.add(M[r][ci], val)
        mats[n] = M
    return chains, mats


def relative_betti_numbers(A, top, max_cells=4_000_000):
    """[rks_0, ..., rks_top]: A^e-generator counts of the minimal resolution
    of A as a bimodule, over A.domain (any exact field)."""
    dom = A.domain
    chains, mats = _relative_complex(A, top, max_cells)
    ranks = {0: 0, 1: 0}
    for n in range(2, top + 2):
        M = mats[n]
        ranks[n] = rank(M, dom) if M and M[0] else 0
    return [len(chains[n]) - ranks[n] - ranks[n + 1] for n in range(top + 1)]
```

- [ ] **Step 5: Route `complexity`** — in `src/quiverlab/invariants/scalar.py`
  replace the `complexity` body's guard line (keep the engine branch verbatim)
  and extend the docstring:

```python
def complexity(A, n):
    """Apparent complexity of A from the minimal A^e (bimodule) resolution's
    term-dimension growth up to degree n. GF(p): the fast engine. Any other
    exact Domain: the E-relative (Cibils) Betti complex on a path-type basis
    (Plan 19) — same generator counts, gated by GF(p) parity tests. Returns
    complexity_of's honest label (int / None / '>=2').

    (…keep the existing Plan-13 and SILENT TRUNCATION PREFIX paragraphs…)
    """
    from quiverlab.engine.scan3 import complexity_of
    from quiverlab.fields.primefield import PrimeField
    if not isinstance(A.domain, PrimeField):
        from quiverlab.invariants.betti import relative_betti_numbers
        return complexity_of(relative_betti_numbers(A, n))
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import minimal_resolution
    eng = to_engine(A)
    p = A.domain.p
    rks, cols, _e, _trunc = minimal_resolution(eng, n, p)
    seq = [max(0, rks[k]) for k in sorted(rks)]
    return complexity_of(seq)
```

- [ ] **Step 6: Flip the loud pin** — in `tests/invariants/test_scalar.py`
  replace `test_complexity_cc_loud` with:

```python
def test_complexity_cc_computes():
    # Plan 19: off GF(p) complexity runs on the relative-Tor Betti complex
    assert truncated_polynomial(2, field=CC).complexity(4) == 1
```

- [ ] **Step 7: Run to verify pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/invariants/test_betti_generic.py tests/invariants/test_scalar.py -p no:cacheprovider`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add src/quiverlab/invariants/pathbasis.py src/quiverlab/invariants/betti.py src/quiverlab/invariants/scalar.py tests/invariants/test_betti_generic.py tests/invariants/test_scalar.py
git commit -m "feat(invariants): complexity over any exact Domain -- relative-Tor Betti complex on path-type bases (Plan 19)"
```

### Task 3: Generic Frobenius / Nakayama / symmetry

**Files:**
- Create: `src/quiverlab/invariants/frobenius.py`
- Modify: `src/quiverlab/core/algebra.py` (`is_frobenius`,
  `nakayama_automorphism`, `is_symmetric` dispatch; GF(p) branches verbatim)
- Test: `tests/invariants/test_frobenius_generic.py` (new)
- Modify: `tests/invariants/test_engine_backed.py` (flip `test_frobenius_cc_loud`)

**Interfaces:**
- Consumes: `path_type_basis(A, what)` (Task 2), `fields.linalg.rank/nullspace/solve`.
- Produces: `is_frobenius_generic(A) -> bool`,
  `frobenius_form_generic(A) -> (lam, G)` (covector as list of Domain scalars +
  Gram matrix, VERIFIED nondegenerate),
  `nakayama_automorphism_generic(A) -> list[list[domain element]]` (columns =
  images), `is_symmetric_generic(A, budget=4096) -> bool`.

- [ ] **Step 1: Write the failing tests** — create
  `tests/invariants/test_frobenius_generic.py`:

```python
"""Generic-Domain Frobenius / Nakayama / symmetry (Plan 19): socle-criterion
knowns over QQ / CC / GF(4), self-certification of the form and nu, the
inner-automorphism symmetry semantics, GF(p) parity with the engine.

Char-2 trap (global constraint): QuantumCI(q) degenerates when q == 0 in the
field — expectations here are per-field, never blanket."""
import pytest

from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.families import QuantumCI
from quiverlab.families.zoo import build_from_record, load_catalog
from quiverlab.fields import QQ


def _rec(name):
    return next(r for r in load_catalog() if r.get("name") == name)


def _commxy(field):
    """k[x,y]/(x^2, y^2, xy): socle span(x, y) is 2-dim => NOT Frobenius."""
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    return Q.algebra(relations=["x^2", "y^2", "x*y", "y*x"], field=field)


def test_frobenius_knowns_any_domain():
    for f in (QQ, GF(4)):
        assert truncated_polynomial(4, field=f).is_frobenius() is True
        assert linear_path_algebra(2, field=f).is_frobenius() is False
        assert _commxy(f).is_frobenius() is False
    assert truncated_polynomial(2, field=CC).is_symmetric() is True


def test_exterior_frobenius_not_symmetric_qq():
    # Lambda(x, y): nu = diag(1, -1, -1, 1) is NOT inner (trivial vertex
    # permutation but no invertible u with nu(a)u = ua) — the case that
    # separates the inner test from the permutation shortcut.
    E2 = QuantumCI(1, field=QQ)
    assert E2.is_frobenius() is True
    assert E2.is_symmetric() is False


def test_commutative_ci_symmetric_qq():
    assert QuantumCI(-1, field=QQ).is_symmetric() is True


def test_quantum_ci_nakayama_qq():
    A = QuantumCI(2, field=QQ)
    dom = A.domain
    assert A.is_frobenius() is True
    assert A.is_symmetric() is False
    N = A.nakayama_automorphism()
    m = A.dim
    ix = A.basis_labels.index("x")
    iy = A.basis_labels.index("y")
    for i in range(m):
        for j in range(m):
            if i != j:
                assert dom.is_zero(N[i][j]), "nu of the quantum CI is diagonal"
    two = dom.coerce(2)
    half = dom.inv(two)
    assert ((dom.eq(N[ix][ix], half) and dom.eq(N[iy][iy], two)) or
            (dom.eq(N[ix][ix], two) and dom.eq(N[iy][iy], half)))


def test_cn_3_2_frobenius_not_symmetric():
    # multi-vertex: kZ_3/rad^2 — Nakayama permutation is the 3-cycle, so
    # Frobenius but NOT symmetric (inner autos fix vertex classes)
    for f in (QQ, GF(4)):
        A = build_from_record(_rec("cn_3_2"), field=f)
        assert A.is_frobenius() is True
        assert A.is_symmetric() is False


def test_nakayama_self_certifies_across_domains():
    """The returned form and nu satisfy their defining equations EXACTLY:
    lambda(f_i f_j) = G[i][j] nondegenerate, lambda(ab) = lambda(b nu(a)),
    nu multiplicative, nu(1) = 1. No oracle needed — the axioms are the gate."""
    from quiverlab.fields import linalg
    from quiverlab.invariants.frobenius import (frobenius_form_generic,
                                                nakayama_automorphism_generic)
    zoo = [truncated_polynomial(3, field=QQ),
           truncated_polynomial(2, field=GF(4)),
           QuantumCI(2, field=QQ),
           QuantumCI(1, field=QQ),
           build_from_record(_rec("cn_3_2"), field=QQ)]
    for A in zoo:
        dom = A.domain
        m = A.dim
        lam, G = frobenius_form_generic(A)
        assert linalg.rank(G, dom) == m
        N = nakayama_automorphism_generic(A)

        def nu(vec):
            out = [dom.zero()] * m
            for j, c in enumerate(vec):
                if dom.is_zero(c):
                    continue
                for i in range(m):
                    out[i] = dom.add(out[i], dom.mul(c, N[i][j]))
            return out

        def lam_of(vec):
            acc = dom.zero()
            for t, c in enumerate(vec):
                acc = dom.add(acc, dom.mul(lam[t], c))
            return acc

        for a in range(m):
            fa = A._basis_vec(a)
            nua = nu(fa)
            for b in range(m):
                fb = A._basis_vec(b)
                # lambda(ab) = lambda(b nu(a))
                assert dom.eq(lam_of(A.multiply(fa, fb)),
                              lam_of(A.multiply(fb, nua)))
                # nu(ab) = nu(a) nu(b)
                got = nu(A.multiply(fa, fb))
                want = A.multiply(nua, nu(fb))
                assert all(dom.eq(x, y) for x, y in zip(got, want))
        one = nu(A.unit)
        assert all(dom.eq(x, y) for x, y in zip(one, A.unit))


def test_gfp_parity_with_engine():
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.coxeter import is_frobenius as eng_frob
    from quiverlab.invariants.frobenius import (is_frobenius_generic,
                                                is_symmetric_generic)
    zoo = [truncated_polynomial(3, field=GF(5)),
           QuantumCI(2, field=GF(5)),
           QuantumCI(1, field=GF(7)),
           _commxy(GF(5)),
           linear_path_algebra(2, field=GF(5))]
    for A in zoo:
        p = A.domain.p
        want = bool(eng_frob(to_engine(A.unit_adapted()), p))
        assert is_frobenius_generic(A) is want, repr(A)
    # symmetry parity where the engine's nu = id semantics and the inner
    # semantics provably coincide (nu diagonal / commutative cases)
    assert is_symmetric_generic(truncated_polynomial(3, field=GF(5))) is True
    assert is_symmetric_generic(QuantumCI(2, field=GF(5))) is False
    A5 = truncated_polynomial(3, field=GF(5))
    assert A5.is_symmetric() is True          # dispatch: engine branch
    assert QuantumCI(2, field=GF(5)).is_symmetric() is False
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/invariants/test_frobenius_generic.py -p no:cacheprovider`
Expected: FAIL — `ModuleNotFoundError: quiverlab.invariants.frobenius`; the
dispatched calls off GF(p) raise FieldError.

- [ ] **Step 3: Implement `src/quiverlab/invariants/frobenius.py`**

```python
"""Frobenius / Nakayama / symmetry over ANY exact Domain (Plan 19).

DECISION (exact, any field): a basic split algebra — path-type basis — is
self-injective, equivalently Frobenius, iff every soc(e_v A) is simple and
v |-> vertex(soc(e_v A)) is a permutation (Nakayama; Skowronski–Yamagata,
Frobenius Algebras I). Both directions are Domain linear algebra:
soc(e_v A) = {x in e_v A : x r = 0}, and sufficiency is the dimension count
P_v -> I(soc P_v) with sum dim P_v = dim A = sum dim I(S_v).

FORM: the socle-dual covector is tried first, then deterministic fallbacks;
nondegeneracy is VERIFIED (Gram rank = dim A), so a returned form is never
wrong. nu = G^{-1} G^T (columns = images); the defining identity
lambda(ab) = lambda(b nu(a)) and multiplicativity are asserted by the test
battery.

SYMMETRY: A symmetric iff Frobenius and nu is INNER. (a) nontrivial
Nakayama vertex permutation -> False (inner automorphisms fix
primitive-idempotent classes); (b) else search span(U),
U = {u : nu(a) u = u a}, for an invertible element on a grid of
min(|field|, dim A + 1) sample values per coordinate: det of the
left-multiplication is a polynomial of degree <= dim A on U, so
(Schwartz–Zippel) vanishing on the whole grid with > dim A distinct values
per variable forces vanishing on U — the sweep is CONCLUSIVE. If the Domain
cannot supply enough distinct samples within `budget` combinations, raise
loudly (never a silent wrong answer)."""
import itertools

from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import nullspace, rank, solve
from quiverlab.invariants.pathbasis import path_type_basis


def _right_socles(A, idem, rad, src):
    """{v: basis of soc(e_v A) as full-coordinate vectors}."""
    dom = A.domain
    out = {}
    for v in idem:
        mine = [v] + [i for i in rad if src[i] == v]
        rows = []
        for b in rad:
            for t in range(A.dim):
                rows.append([A.T[i][b][t] for i in mine])
        if rows:
            ker = nullspace(rows, dom)
        else:                                   # no radical at all: soc = e_v A
            ker = [[dom.one() if k == pos else dom.zero()
                    for k in range(len(mine))] for pos in range(len(mine))]
        socs = []
        for vec in ker:
            full = [dom.zero()] * A.dim
            for pos, i in enumerate(mine):
                full[i] = vec[pos]
            socs.append(full)
        out[v] = socs
    return out


def nakayama_data(A):
    """(ok, perm, gens): ok = the socle criterion holds; perm maps each
    vertex-idempotent index v to the vertex of soc(e_v A); gens[v] = the
    1-dim socle generator (full coordinates)."""
    dom = A.domain
    idem, rad, src, tgt = path_type_basis(A, "the Frobenius test")
    socs = _right_socles(A, idem, rad, src)
    perm, gens = {}, {}
    for v in idem:
        if len(socs[v]) != 1:
            return False, None, None
        s = socs[v][0]
        nz = [i for i, c in enumerate(s) if not dom.is_zero(c)]
        targets = {tgt[i] if i in tgt else i for i in nz}
        if len(targets) != 1:
            raise QuiverlabError(
                "1-dimensional socle generator is not corner-homogeneous — "
                "impossible for a path-type basis (bug)")
        perm[v] = targets.pop()
        gens[v] = s
    if len(set(perm.values())) != len(idem):
        return False, None, None
    return True, perm, gens


def is_frobenius_generic(A):
    """Exact, conclusive, any Domain (socle criterion on a path-type basis)."""
    ok, _perm, _gens = nakayama_data(A)
    return ok


def _lam_of(lam, vec, dom):
    acc = dom.zero()
    for t, c in enumerate(vec):
        if not dom.is_zero(c):
            acc = dom.add(acc, dom.mul(lam[t], c))
    return acc


def frobenius_form_generic(A):
    """(lam, G): a VERIFIED-nondegenerate Frobenius covector (lam(x) =
    sum_t lam[t] x_t) and its Gram matrix G[i][j] = lam(f_i f_j).
    Raises QuiverlabError if A is not Frobenius, and loudly if no candidate
    passes the exact rank check (never returns an unverified form)."""
    dom = A.domain
    m = A.dim
    ok, _perm, gens = nakayama_data(A)
    if not ok:
        raise QuiverlabError(
            "algebra is not Frobenius (socle criterion), so it has no "
            "Frobenius form / Nakayama automorphism")

    def gram(lam):
        return [[_lam_of(lam, A.T[i][j], dom) for j in range(m)] for i in range(m)]

    socdual = [dom.zero()] * m
    for s in gens.values():
        for t, c in enumerate(s):
            socdual[t] = dom.add(socdual[t], c)
    cands = [socdual]
    cands.extend([[dom.one() if t == c else dom.zero() for t in range(m)]
                  for c in range(m)])
    cands.append([dom.one()] * m)
    for lam in cands:
        G = gram(lam)
        if rank(G, dom) == m:
            return lam, G
    raise QuiverlabError(
        "algebra IS Frobenius (socle criterion) but no deterministic "
        "candidate covector was nondegenerate — please report this algebra",
        hint="the socle-dual functional is expected to work for every "
             "path-type basis")


def nakayama_automorphism_generic(A):
    """nu as a matrix of Domain elements (columns = images): N = G^{-1} G^T,
    i.e. the solution of G N = G^T, from the defining identity
    lam(ab) = lam(b nu(a))."""
    lam, G = frobenius_form_generic(A)
    dom = A.domain
    m = A.dim
    N = [[dom.zero()] * m for _ in range(m)]
    for j in range(m):
        col = solve(G, [G[j][i] for i in range(m)], dom)   # (G^T) column j
        for i in range(m):
            N[i][j] = col[i]
    return N


def _sample_values(dom, need):
    """Up to `need` pairwise-distinct Domain elements, deterministically:
    integer coercions first (0, 1, 2, ...). Returns fewer than `need` only
    when the Domain's prime subfield is exhausted (small finite fields)."""
    vals = []
    k = 0
    while len(vals) < need and k < 4 * need + 8:
        cand = dom.coerce(k)
        if all(not dom.eq(cand, v) for v in vals):
            vals.append(cand)
        k += 1
    return vals


def is_symmetric_generic(A, budget=4096):
    """A symmetric iff Frobenius and nu is inner — exact; loud when the
    sample sweep cannot be made conclusive within `budget`."""
    ok, perm, _gens = nakayama_data(A)
    if not ok:
        return False
    if any(v != w for v, w in perm.items()):
        return False                    # inner autos fix vertex classes
    dom = A.domain
    m = A.dim
    N = nakayama_automorphism_generic(A)
    # U = {u : nu(f_a) u = u f_a for all a}
    rows = []
    for a in range(m):
        nua = [N[i][a] for i in range(m)]
        for t in range(m):
            row = []
            for k in range(m):
                left = dom.zero()       # coeff of u_k in (nu(f_a) u)_t
                for i in range(m):
                    if not dom.is_zero(nua[i]):
                        left = dom.add(left, dom.mul(nua[i], A.T[i][k][t]))
                row.append(dom.sub(left, A.T[k][a][t]))
            rows.append(row)
    U = nullspace(rows, dom)
    if not U:
        return False

    def invertible(u):
        L = [[dom.zero()] * m for _ in range(m)]   # left multiplication by u
        for i, c in enumerate(u):
            if dom.is_zero(c):
                continue
            for j in range(m):
                for t, x in enumerate(A.T[i][j]):
                    if not dom.is_zero(x):
                        L[t][j] = dom.add(L[t][j], dom.mul(c, x))
        return rank(L, dom) == m

    from quiverlab.fields.primefield import PrimeField
    samples = _sample_values(dom, m + 1)
    # Conclusive iff > dim A distinct samples (Schwartz–Zippel), OR the sweep
    # enumerates the WHOLE coefficient space — which integer coercions give
    # only over the prime field itself. Over GF(p^n), n > 1, with p <= dim A,
    # integer samples span just the prime subfield of each coefficient: that
    # is NOT the whole space, so it must refuse loudly, never guess.
    whole_field = isinstance(dom, PrimeField) and len(samples) == dom.characteristic
    if len(samples) < m + 1 and not whole_field:
        raise QuiverlabError(
            "symmetry test: the Domain supplies too few distinct sample "
            "values for a conclusive Schwartz–Zippel sweep",
            hint="decide symmetry over GF(p) via the engine, or over a "
                 "larger field")
    if len(samples) ** len(U) > budget:
        raise QuiverlabError(
            f"symmetry test: the conclusive sweep needs "
            f"{len(samples)}^{len(U)} > budget={budget} combinations",
            hint="raise the budget")
    for coeffs in itertools.product(samples, repeat=len(U)):
        u = [dom.zero()] * m
        for cf, base in zip(coeffs, U):
            if dom.is_zero(cf):
                continue
            for t, x in enumerate(base):
                u[t] = dom.add(u[t], dom.mul(cf, x))
        if invertible(u):
            return True
    return False
```

- [ ] **Step 4: Dispatch** — in `src/quiverlab/core/algebra.py` replace the three
  method bodies (GF(p) branches verbatim; `_require_prime_field` calls dropped):

```python
    def nakayama_automorphism(self):
        """Nakayama automorphism nu as a matrix (columns = images) in the
        algebra's basis. GF(p): integer matrix via the engine (unit-adapted
        basis). Other exact Domains: Domain-element matrix on the path-type
        basis (Plan 19). Loud if not Frobenius."""
        from quiverlab.fields.primefield import PrimeField
        if isinstance(self.domain, PrimeField):
            from quiverlab.engine.adapter import to_engine
            from quiverlab.engine.coxeter import nakayama_automorphism
            S, _ = nakayama_automorphism(to_engine(self.unit_adapted()), self.domain.p)
            return [[int(S[i, j]) for j in range(S.shape[1])] for i in range(S.shape[0])]
        from quiverlab.invariants.frobenius import nakayama_automorphism_generic
        return nakayama_automorphism_generic(self)

    def is_frobenius(self):
        """Is the algebra Frobenius? GF(p): engine form search. Other exact
        Domains: the exact socle criterion on a path-type basis (Plan 19)."""
        from quiverlab.fields.primefield import PrimeField
        if isinstance(self.domain, PrimeField):
            from quiverlab.engine.adapter import to_engine
            from quiverlab.engine.coxeter import is_frobenius
            return bool(is_frobenius(to_engine(self.unit_adapted()), self.domain.p))
        from quiverlab.invariants.frobenius import is_frobenius_generic
        return is_frobenius_generic(self)

    def is_symmetric(self):
        """Is the algebra symmetric? GF(p): Frobenius with identity Nakayama
        automorphism (engine). Other exact Domains: Frobenius with INNER
        Nakayama automorphism (Plan 19, the definitional test)."""
        from quiverlab.fields.primefield import PrimeField
        if isinstance(self.domain, PrimeField):
            if not self.is_frobenius():
                return False
            from quiverlab.engine.adapter import to_engine
            from quiverlab.engine.coxeter import is_identity, nakayama_automorphism
            E = to_engine(self.unit_adapted())
            S, _ = nakayama_automorphism(E, self.domain.p)
            return bool(is_identity(S, self.domain.p))
        from quiverlab.invariants.frobenius import is_symmetric_generic
        return is_symmetric_generic(self)
```

  Then DELETE the now-unused `_require_prime_field` method entirely.

- [ ] **Step 5: Flip the loud pin** — in `tests/invariants/test_engine_backed.py`
  replace `test_frobenius_cc_loud` with:

```python
def test_frobenius_cc_computes():
    # Plan 19: off GF(p) the exact socle criterion serves any Domain
    A = truncated_polynomial(2, field=CC)
    assert A.is_frobenius() is True
    assert A.is_symmetric() is True
```

- [ ] **Step 6: Run to verify pass**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/invariants/test_frobenius_generic.py tests/invariants/test_engine_backed.py -p no:cacheprovider`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/quiverlab/invariants/frobenius.py src/quiverlab/core/algebra.py tests/invariants/test_frobenius_generic.py tests/invariants/test_engine_backed.py
git commit -m "feat(invariants): Frobenius/Nakayama/symmetry over any exact Domain -- socle criterion + inner-nu test (Plan 19)"
```

### Task 4: Residual refusal surface, honestly worded

**Files:**
- Test: `tests/invariants/test_refusal_surface.py` (new)

**Interfaces:**
- Consumes: `Algebra.from_structure_constants(T, unit, field=…)`,
  the Tasks 1–3 dispatches, `path_type_basis`'s FieldError.
- Produces: the pinned refusal contract Task 5's docs describe.

- [ ] **Step 1: Write the tests** — create
  `tests/invariants/test_refusal_surface.py`:

```python
"""Plan 19 refusal contract: off GF(p), only path-basis-needing invariants
on quiver-less algebras refuse — with an honest message (no 'later phase'
promise anywhere in src/). cyclic_homology refuses NOWHERE."""
import pathlib

import pytest

import quiverlab
from quiverlab import CC, truncated_polynomial
from quiverlab.core.algebra import Algebra
from quiverlab.errors import FieldError


def _dual_numbers_sc(field=CC):
    """k[x]/(x^2) via raw structure constants: NO quiver, NO e_ labels."""
    T = [[[1, 0], [0, 1]],
         [[0, 1], [0, 0]]]
    return Algebra.from_structure_constants(T, [1, 0], field=field)


def test_structure_constants_refusals_are_honest():
    A = _dual_numbers_sc(CC)
    for call in (A.is_frobenius, A.is_symmetric, A.nakayama_automorphism,
                 lambda: A.complexity(3)):
        with pytest.raises(FieldError) as ei:
            call()
        msg = str(ei.value)
        assert "later phase" not in msg
        assert "quiver" in msg or "path-type" in msg


def test_cyclic_homology_needs_no_quiver():
    # the bar mixed complex serves ANY unital algebra: same dims from the
    # structure-constants presentation and the quiver presentation
    got = _dual_numbers_sc(CC).cyclic_homology(2)
    want = truncated_polynomial(2, field=CC).cyclic_homology(2)
    assert got.dims == want.dims


def test_no_later_phase_promise_left_in_src():
    root = pathlib.Path(quiverlab.__file__).parent
    hits = [str(p) for p in root.rglob("*.py")
            if "later phase that generalizes" in p.read_text()]
    assert hits == []
```

- [ ] **Step 2: Run to verify pass** (Tasks 1–3 already implement the contract;
  this task pins it — expect immediate PASS; if any assertion fails, fix the
  message in `invariants/pathbasis.py`, never weaken the test)

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/invariants/test_refusal_surface.py -p no:cacheprovider`
Expected: PASS

- [ ] **Step 3: Run the full fast bucket**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast -p no:cacheprovider`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add tests/invariants/test_refusal_surface.py
git commit -m "test(invariants): pin the Plan-19 refusal contract -- honest messages, no 'later phase' promise"
```

### Task 5: Docs + suites + backlog

**Files:**
- Modify: `docs/internals/06-invariants.md` (Plan-19 sections: generic HC mixed
  complex; the relative-Tor Betti complex + why H_n = rks_n; socle criterion +
  the dimension-count proof sketch; form/ν construction; the symmetry decision
  incl. the Schwartz–Zippel argument and the loud budgets; the dispatch table
  GF(p)-engine vs generic vs honest refusal)
- Modify: `CLAUDE.md` (Status paragraph gains Plan 19; the "engine-backed
  invariants are GF(p)-only" claim in the `invariants/scalar.py` docstring quote
  is superseded — reword the Status line accordingly)
- Modify: `docs/plans/ROADMAP.md` (DELIVERED row 19)
- Modify: `docs/plans/DEEPER-ENGINES-BACKLOG.md` (tick item 6 with plan number +
  date)

- [ ] **Step 1: Doc edits** (facts per the Architecture block; internals chapter
  follows the fixed format of `docs/internals/06-invariants.md`'s existing
  sections; state the is_symmetric semantic note: engine = "ν literally id for
  its form", generic = "ν inner" — coinciding on the validated zoo, gated by
  `test_gfp_parity_with_engine`)

- [ ] **Step 2: Full suites**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast -p no:cacheprovider`
Expected: PASS
Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m deep -p no:cacheprovider`
Expected: PASS (~30 min; run detached per the long-suite protocol)
Run: `.venv/bin/mkdocs build --strict`
Expected: exit 0

- [ ] **Step 3: Commit**

```bash
git add docs/internals/06-invariants.md CLAUDE.md docs/plans/ROADMAP.md docs/plans/DEEPER-ENGINES-BACKLOG.md
git commit -m "docs: Plan-19 status -- engine-backed invariants over any exact Domain"
```

## Validation matrix

1. HC: generic ≡ engine over GF(3)/GF(5) (three algebras); (b,B) identities
   exact over QQ; λ-complex second model agrees over QQ; semisimple closed forms
   over QQ + GF(4); CC ≡ QQ dims.
2. Betti: generic ≡ engine `rks` over GF(3)/GF(5) on six algebras incl.
   multi-vertex (`comm_square`, `cn_3_2`) and straddling monomial; d∘d = 0
   exact; char-0 closed forms (monomial, CI, gldim-2); complexity labels over
   CC / GF(4) / QQ.
3. Frobenius: socle-criterion knowns over QQ/CC/GF(4) incl. the 2-dim-socle
   refuter and the multi-vertex Nakayama 3-cycle; form + ν self-certify
   (defining identity, multiplicativity, ν(1)=1) over QQ and GF(4); exterior
   algebra separates inner-test from permutation shortcut; GF(p) parity with
   the engine.
4. Refusal contract: quiver-less + off-GF(p) refuses ONLY where a path-type
   basis is genuinely needed, message honest; `cyclic_homology` never refuses;
   no "later phase that generalizes" string survives in `src/`.
5. Full fast + deep suites green; `mkdocs build --strict` exit 0.

## Status

- [ ] Not yet executed.
