# Plan 35 Implementation Plan — HH products: cup, cap, Gerstenhaber bracket, Connes B

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Surface the existing cup/cap/bracket/Connes-B implementations as public
`Algebra` methods, four compute kinds across all tiers + Pyodide GUI, worked-steps
report chapter, and product blocks in the six curated cache examples — gating v0.1.0.

**Architecture:** New assembly module `src/quiverlab/hochschild/products.py` builds
structure-constant tables from the GF(p) bar engine (`engine/tt_calculus.py`) or the
CS route (`resolutions_cs/products.py`, new, using `cs_hh_basis` + `native_cup`/
`native_cap`); `core/algebra.py` gains four dispatching methods mirroring
`hochschild_cohomology`; `hpc/spec.py` + `docs/gui/runner.py` serve four new kinds;
GUI/report layers render the blocks.

**Tech Stack:** pure Python + numpy int64 (GF(p) engine path), `fields.linalg`
(Domain-generic path), pytest with the Plan-32 oracle markers.

## Global Constraints

- **No floats in `src/`** — AST gate `tests/test_no_floats.py`; all constants exact
  (int mod p, or Domain elements stringified via `str()` at the block boundary).
- **Engine internals stay internal** — `webapp/`, `docs/gui/` never import
  `quiverlab.engine.*`; only `src/quiverlab/` modules may.
- Python is **always** `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q`.
- Test buckets auto-assign by directory (`tests/conftest.py`): `tests/engine/`,
  `tests/resolutions_cs/` → deep; `tests/hochschild/`, `tests/webapp/`,
  `tests/trace/`, `tests/hpc/` → fast. Place accordingly and mark oracle classes
  (`oracle_selfcert` / `oracle_crossengine` / `oracle_literature` / `qpa`).
- Both vendored GUI copies stay byte-identical:
  `docs/gui/gui.js` ≡ `webapp/static/gui/gui.js` (gated by an existing test).
- Every user-facing webapp string gets EN + ES i18n keys.
- Conventional commits; suite green at every commit.
- **Spec deviation to encode (discovered during plan research, amend the spec in
  Task 1):** the Gerstenhaber bracket is **GF(p)-only** — the transport route
  (`Comparison`) itself requires `PrimeField` (`comparison.py:66`), so over QQ/CC
  the bracket has NO route and refuses loudly (the spec's §2 row-2 bracket cell
  overpromised "CS transport" for QQ/CC). Also bracket tables cover pairs with
  `p, q >= 1` only (`circle_cochain` requires `p >= 1`; the degree-0 insertion
  action is out of scope, stated in the block and docstring).

---

### Task 1: Result containers in `hochschild/products.py` + spec amendment

**Files:**
- Create: `src/quiverlab/hochschild/products.py`
- Create: `tests/hochschild/test_products_containers.py`
- Modify: `docs/plans/2026-08-01-plan-35-hh-products-surface.md` (spec §2 bracket
  cell + §1 bracket note, per Global Constraints)

**Interfaces:**
- Produces: `ProductTable(kind, degrees, out_degree, dims, constants)` (frozen
  dataclass; `constants` is a nested tuple `[k][i][j]` of `str` exact entries);
  `HHProducts(kind, top, tables, engine, basis, window, references)` with
  `.blocks() -> dict`; `ConnesB(top, hh_dims, matrices, ranks, engine, references)`
  with `.blocks() -> dict`. Later tasks (5, 8, 9, 11) consume exactly these names.

- [ ] **Step 1: Write the failing test**

```python
# tests/hochschild/test_products_containers.py
"""Container contract: exact-string constants, degenerate degrees present,
canonical .blocks() shape shared by both runners (Plan 35)."""
import pytest

from quiverlab.hochschild.products import ProductTable, HHProducts, ConnesB


def test_product_table_is_frozen_and_stringly_exact():
    t = ProductTable(kind="cup", degrees=(1, 1), out_degree=2,
                     dims=(2, 2, 1), constants=((("3", "0"), ("0", "1")),))
    assert t.constants[0][0][0] == "3"
    with pytest.raises(Exception):
        t.kind = "cap"


def test_hhproducts_blocks_shape():
    t = ProductTable(kind="cup", degrees=(0, 0), out_degree=0,
                     dims=(1, 1, 1), constants=((("1",),),))
    hp = HHProducts(kind="cup", top=0, tables={(0, 0): t},
                    engine="hanlab engine (F_p fast rank)",
                    basis="bar/GF(7)", window=None, references=["cup"])
    b = hp.blocks()
    assert b["kind"] == "cup" and b["top"] == 0
    assert b["basis"] == "bar/GF(7)"
    assert b["tables"][0]["degrees"] == [0, 0]
    assert b["tables"][0]["constants"] == [[["1"]]]
    assert "window" not in b            # cup carries no window key


def test_bracket_blocks_carry_window():
    hp = HHProducts(kind="bracket", top=3, tables={},
                    engine="hanlab engine (F_p fast rank)",
                    basis="bar/GF(7)", window=2, references=["bracket"])
    assert hp.blocks()["window"] == 2


def test_connesb_blocks_shape():
    cb = ConnesB(top=2, hh_dims=[1, 1, 1],
                 matrices={0: [["1"]], 1: [["0"]]}, ranks={0: 1, 1: 0},
                 engine="engine (b,B) GF(7)", references=["cyclic"])
    b = cb.blocks()
    assert b["kind"] == "connes_b" and b["ranks"] == {"0": 1, "1": 0}
    assert b["matrices"]["1"] == [["0"]]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_products_containers.py -q`
Expected: FAIL — `ModuleNotFoundError: quiverlab.hochschild.products`

- [ ] **Step 3: Write the containers**

```python
# src/quiverlab/hochschild/products.py
"""Plan 35 -- the HH product surface: structure-constant tables for the cup
product, the cap module action, the Gerstenhaber bracket, and the induced
Connes differential, as frozen result objects with one canonical block
serialization (consumed identically by hpc/spec.py and docs/gui/runner.py).

Constants are ALWAYS exact strings at the boundary (`str(entry)`): ints mod p
on the GF(p) routes, Domain reprs on the CS route. No floats can appear (the
AST gate scans this file)."""
from dataclasses import dataclass


@dataclass(frozen=True)
class ProductTable:
    kind: str            # "cup" | "cap" | "bracket"
    degrees: tuple       # (p, q) -- for cap, (p, n): HH^p (x) HH_n -> HH_{n-p}
    out_degree: int
    dims: tuple          # (dim_left, dim_right, dim_out)
    constants: tuple     # [k][i][j] -> str: left_i * right_j = sum_k c^k_ij out_k

    def as_dict(self):
        return {"degrees": list(self.degrees), "out_degree": self.out_degree,
                "dims": list(self.dims),
                "constants": [[[c for c in row] for row in mat]
                              for mat in self.constants]}


class HHProducts:
    """A family of product tables up to `top`. kind in {"cup","cap","bracket"}."""

    def __init__(self, kind, top, tables, engine, basis, window, references):
        self.kind = kind
        self.top = top
        self.tables = dict(tables)     # {(p, q): ProductTable}
        self.engine = engine
        self.basis = basis
        self.window = window           # int for bracket (served window), else None
        self.references = list(references)

    def blocks(self):
        out = {"kind": self.kind, "top": self.top, "engine": self.engine,
               "basis": self.basis,
               "tables": [self.tables[k].as_dict()
                          for k in sorted(self.tables)],
               "references": list(self.references)}
        if self.window is not None:
            out["window"] = self.window
        return out

    def __repr__(self):
        return (f"<HHProducts {self.kind} top={self.top} "
                f"tables={len(self.tables)} basis={self.basis!r}>")


class ConnesB:
    """Induced Connes differentials B: HH_n -> HH_{n+1} for 0 <= n < top."""

    def __init__(self, top, hh_dims, matrices, ranks, engine, references):
        self.top = top
        self.hh_dims = list(hh_dims)   # dim HH_0..HH_top
        self.matrices = dict(matrices) # {n: rows of str, shape hh_dims[n+1] x hh_dims[n]}
        self.ranks = dict(ranks)       # {n: int}
        self.engine = engine
        self.references = list(references)

    def blocks(self):
        return {"kind": "connes_b", "top": self.top,
                "hh_dims": list(self.hh_dims),
                "matrices": {str(n): self.matrices[n] for n in sorted(self.matrices)},
                "ranks": {str(n): self.ranks[n] for n in sorted(self.ranks)},
                "engine": self.engine, "references": list(self.references)}

    def __repr__(self):
        return f"<ConnesB top={self.top} ranks={self.ranks}>"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_products_containers.py -q`
Expected: 4 passed

- [ ] **Step 5: Amend the spec** — in
`docs/plans/2026-08-01-plan-35-hh-products-surface.md` §2 table, change the
row-2 bracket cell to: *"no route — loud refusal off GF(p) (the transport
itself is the GF(p) tt facade, `comparison.py:66`)"*; add to §1 under
`gerstenhaber_brackets`: *"tables cover pairs `p, q >= 1` (the degree-0
insertion action is out of scope, stated in the block)"*.

- [ ] **Step 6: Commit**

```bash
git add src/quiverlab/hochschild/products.py tests/hochschild/test_products_containers.py docs/plans/2026-08-01-plan-35-hh-products-surface.md
git commit -m "feat(hochschild): Plan-35 product result containers + spec bracket-scope amendment"
```

---

### Task 2: GF(p) bar-route table builders

**Files:**
- Modify: `src/quiverlab/hochschild/products.py` (append)
- Create: `tests/hochschild/test_products_gfp.py`

**Interfaces:**
- Consumes: `engine.tt_calculus.cup_product_matrix / cap_product_matrix /
  gerstenhaber_bracket_matrix(alg, p, q, prime) -> (C, dl, dr, dout)` where `C`
  is int64 of shape `(dout, dl, dr)`; `engine.adapter.to_engine(A.unit_adapted())`.
- Produces: `gfp_product_tables(A, kind, top, max_cells) -> HHProducts` — the
  GF(p) route Task 5 dispatches to. Raises `DepthLimitError` past the bar wall
  (propagating the underlying guard), `QuiverlabError` never (routing errors live
  in Task 5).

- [ ] **Step 1: Write the failing test**

```python
# tests/hochschild/test_products_gfp.py
"""GF(p) bar-route tables vs direct tt_calculus calls (Plan 35).

k[x]/(x^2) over GF(7): dim HH^n = 1 for all n (the classic pin); the cup ring
in char != 2 has HH^odd squaring to 0 and HH^even polynomial -- we assert the
table entries EQUAL the raw tt_calculus tensor, and the unit law on (0, q)."""
import pytest

import quiverlab as ql

pytestmark = [pytest.mark.oracle_crossengine]


@pytest.fixture(scope="module")
def A():
    return ql.TruncatedPolynomial(2, field=ql.GF(7))


def test_cup_tables_match_tt(A):
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.tt_calculus import cup_product_matrix
    from quiverlab.hochschild.products import gfp_product_tables
    hp = gfp_product_tables(A, "cup", 3, max_cells=4_000_000)
    assert hp.kind == "cup" and hp.basis == "bar/GF(7)"
    E = to_engine(A.unit_adapted())
    for (p, q), table in hp.tables.items():
        C, dl, dr, dout = cup_product_matrix(E, p, q, 7)
        assert table.dims == (dl, dr, dout)
        for k in range(dout):
            for i in range(dl):
                for j in range(dr):
                    assert table.constants[k][i][j] == str(int(C[k, i, j]))


def test_cup_pairs_cover_exactly_p_plus_q_le_top(A):
    from quiverlab.hochschild.products import gfp_product_tables
    hp = gfp_product_tables(A, "cup", 2, max_cells=4_000_000)
    assert sorted(hp.tables) == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)]


def test_cap_pairs_cover_p_le_n_le_top(A):
    from quiverlab.hochschild.products import gfp_product_tables
    hp = gfp_product_tables(A, "cap", 2, max_cells=4_000_000)
    assert sorted(hp.tables) == [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    t = hp.tables[(1, 2)]
    assert t.out_degree == 1            # HH^1 (x) HH_2 -> HH_1


def test_bracket_pairs_need_p_and_q_ge_1(A):
    from quiverlab.hochschild.products import gfp_product_tables
    hp = gfp_product_tables(A, "bracket", 3, max_cells=4_000_000)
    assert all(p >= 1 and q >= 1 for (p, q) in hp.tables)
    assert (1, 1) in hp.tables and hp.tables[(1, 1)].out_degree == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_products_gfp.py -q`
Expected: FAIL — `ImportError: cannot import name 'gfp_product_tables'`

- [ ] **Step 3: Implement the builder** (append to `products.py`)

```python
def _pairs(kind, top):
    """The degree pairs a `kind` table family covers up to `top`."""
    if kind == "cup":
        return [(p, q) for p in range(top + 1) for q in range(top + 1 - p)]
    if kind == "cap":
        return [(p, n) for n in range(top + 1) for p in range(n + 1)]
    if kind == "bracket":
        return [(p, q) for p in range(1, top + 2) for q in range(1, top + 2 - p + 1)
                if p + q - 1 <= top]
    raise ValueError(f"unknown product kind {kind!r}")


def gfp_product_tables(A, kind, top, max_cells):
    """The GF(p) bar-route table family: tt_calculus structure constants on the
    bar HH basis. A must be over a prime field (the caller routes)."""
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine import tt_calculus as TT
    prime = A.domain.p
    E = to_engine(A.unit_adapted())
    fn = {"cup": TT.cup_product_matrix, "cap": TT.cap_product_matrix,
          "bracket": TT.gerstenhaber_bracket_matrix}[kind]
    out_deg = {"cup": lambda p, q: p + q, "cap": lambda p, n: n - p,
               "bracket": lambda p, q: p + q - 1}[kind]
    tables = {}
    for (p, q) in _pairs(kind, top):
        C, dl, dr, dout = fn(E, p, q, prime)
        tables[(p, q)] = ProductTable(
            kind=kind, degrees=(p, q), out_degree=out_deg(p, q),
            dims=(dl, dr, dout),
            constants=tuple(tuple(tuple(str(int(C[k, i, j])) for j in range(dr))
                                  for i in range(dl)) for k in range(dout)))
    window = top if kind == "bracket" else None
    return HHProducts(kind=kind, top=top, tables=tables,
                      engine="hanlab engine (F_p fast rank)",
                      basis=f"bar/GF({prime})", window=window,
                      references=_REFERENCES[kind])


_REFERENCES = {"cup": ["cup", "gerstenhaber"],
               "cap": ["cup", "gerstenhaber"],
               "bracket": ["bracket", "gerstenhaber"],
               "connes_b": ["cyclic"]}
```

Note: the bar route has no separate wall check here — `cochain_basis` growth is
bounded by the same `max_cells` discipline as HH itself; thread `max_cells` into
a pre-check identical to `hochschild/bar.py::_check_cells` semantics:
before each `fn(E, p, q, prime)` call compute
`len(cochain_basis(E, p)) * len(cochain_basis(E, q))` and raise
`DepthLimitError` when it exceeds `max_cells` (import from `quiverlab.errors`,
hint text: `"raise max_cells or lower top"`). Include this in the implementation
(it is what Task 5's CS fallback catches).

- [ ] **Step 4: Run test to verify it passes**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_products_gfp.py -q`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/hochschild/products.py tests/hochschild/test_products_gfp.py
git commit -m "feat(hochschild): GF(p) bar-route product tables (cup/cap/bracket)"
```

---

### Task 3: CS-route table builders (`resolutions_cs/products.py`)

**Files:**
- Create: `src/quiverlab/resolutions_cs/products.py`
- Create: `tests/resolutions_cs/test_products_cs.py`

**Interfaces:**
- Consumes: `resolutions_cs.homology.cs_hh_basis(A, n, side, max_cells)` →
  list of Domain coordinate vectors; `resolutions_cs.cup.native_cup(res, f_vec,
  p, g_vec, q)`; `resolutions_cs.cap.native_cap(res, f_vec, p, z_vec, n)`;
  `resolutions_cs.build.reduction_system_of(A)`;
  `resolutions_cs.resolution.ChouhySolotarResolution(A, rs, max_degree, max_cells)`;
  `fields.linalg.solve`, `fields.linalg.nullspace`.
- Produces: `cs_product_tables(A, kind, top, max_cells) -> HHProducts` for
  `kind in ("cup", "cap")` over ANY exact Domain (bracket is NOT served here —
  Global Constraints). Basis string: `"cs/" + A.domain.name`.

- [ ] **Step 1: Write the failing test**

```python
# tests/resolutions_cs/test_products_cs.py
"""CS-route product tables (Plan 35): Domain-generic cup/cap on the CS basis.

Deep-bucket (tests/resolutions_cs -> deep). Oracles: the QQ smoke is
self-certifying (unit row + dims equal cs dims); the GF(p) cross-engine
equality against the bar tables lives in tests/hochschild/test_products_identities.py
(Task 6) -- HERE we pin shape, unit law, and Domain-genericity."""
import pytest

import quiverlab as ql

pytestmark = [pytest.mark.oracle_selfcert]


@pytest.fixture(scope="module")
def A_qq():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1), "y": (1, 1)})
    return Q.algebra(relations=["x*x", "y*y", "x*y + y*x"], field=ql.QQ)


def test_cs_cup_over_qq_unit_law(A_qq):
    from quiverlab.resolutions_cs.homology import cs_cohomology_dims
    from quiverlab.resolutions_cs.products import cs_product_tables
    hp = cs_product_tables(A_qq, "cup", 2, max_cells=4_000_000)
    assert hp.basis.startswith("cs/")
    dims = list(cs_cohomology_dims(A_qq, 2).dims)
    # unit law: the (0, q) table with the identity class of HH^0 = Z(A) acting
    # as identity -- for the class index of 1_A, the table column is the identity.
    t = hp.tables[(0, 1)]
    assert t.dims == (dims[0], dims[1], dims[1])


def test_cs_cap_over_qq_shapes(A_qq):
    from quiverlab.resolutions_cs.homology import cs_homology_dims, cs_cohomology_dims
    from quiverlab.resolutions_cs.products import cs_product_tables
    hp = cs_product_tables(A_qq, "cap", 2, max_cells=4_000_000)
    hdims = list(cs_homology_dims(A_qq, 2).dims)
    cdims = list(cs_cohomology_dims(A_qq, 2).dims)
    t = hp.tables[(1, 2)]
    assert t.dims == (cdims[1], hdims[2], hdims[1])


def test_cs_bracket_refuses():
    from quiverlab.resolutions_cs.products import cs_product_tables
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=ql.QQ)
    with pytest.raises(ql.QuiverlabError):
        cs_product_tables(A, "bracket", 2, max_cells=4_000_000)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/resolutions_cs/test_products_cs.py -q`
Expected: FAIL — `ModuleNotFoundError: quiverlab.resolutions_cs.products`

- [ ] **Step 3: Implement**

```python
# src/quiverlab/resolutions_cs/products.py
"""Plan 35 -- Domain-generic CS product tables: cup and cap on the CS HH basis
via the Plan-20/21 native collapses of the lifted diagonal. Any exact Domain,
any degree (no bar window). The bracket is NOT served here: its only route is
the GF(p) tt facade (see the Plan-35 spec amendment)."""
from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import solve
from quiverlab.hochschild.products import HHProducts, ProductTable, _pairs, _REFERENCES


def _class_coords(vec, reps, image_cols, dom):
    """Coordinates of `vec` in the class basis `reps`, modulo `image_cols`:
    solve [image | reps] x = vec and read off the reps segment. Loud when the
    vector is not a (co)cycle representative (descent failure)."""
    cols = list(image_cols) + list(reps)
    if not cols:
        if any(not dom.is_zero(dom.coerce(v)) for v in vec):
            raise QuiverlabError("product failed to land in the zero space")
        return []
    Mat = [[dom.coerce(cols[c][r]) for c in range(len(cols))]
           for r in range(len(vec))]
    x = solve(Mat, [dom.coerce(v) for v in vec], dom)
    if x is None:
        raise QuiverlabError("product failed to descend to (co)homology "
                             "(not in the cycle span) -- this is a bug, report it")
    return [str(c) for c in x[len(image_cols):]]


def cs_product_tables(A, kind, top, max_cells):
    if kind == "bracket":
        raise QuiverlabError(
            "the Gerstenhaber bracket has no CS-native route",
            hint="the bracket is served over GF(p) within the bar window only; "
                 "construct the algebra over GF(p)")
    if kind not in ("cup", "cap"):
        raise QuiverlabError(f"unknown product kind {kind!r}")
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.homology import cs_hh_basis, _require_admissible
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    from quiverlab.resolutions_cs.cup import native_cup
    from quiverlab.resolutions_cs.cap import native_cap

    rs = reduction_system_of(A)
    _require_admissible(rs)
    dom = A.domain
    res = ChouhySolotarResolution(A, rs, max_degree=top + 2, max_cells=max_cells)
    coh = {n: cs_hh_basis(A, n, "coh", max_cells=max_cells) for n in range(top + 1)}
    hom = ({n: cs_hh_basis(A, n, "hom", max_cells=max_cells) for n in range(top + 1)}
           if kind == "cap" else {})

    def _image(n, side):
        M = res.matrix(n - 1, "coh") if side == "coh" else res.matrix(n + 1, "hom")
        if side == "coh" and n == 0:
            return []
        if not M or not M[0]:
            return []
        return [[M[r][c] for r in range(len(M))] for c in range(len(M[0]))]

    tables = {}
    for (p, q) in _pairs(kind, top):
        if kind == "cup":
            left, right, out_n, side = coh[p], coh[q], p + q, "coh"
            prod = lambda f, z: native_cup(res, f, p, z, q)
            out_reps = coh.get(out_n)
        else:
            left, right, out_n, side = coh[p], hom[q], q - p, "hom"
            prod = lambda f, z: native_cap(res, f, p, z, q)
            out_reps = hom.get(out_n)
        dl, dr, dout = len(left), len(right), len(out_reps)
        img = _image(out_n, side)
        consts = [[[None] * dr for _ in range(dl)] for _ in range(dout)]
        for i in range(dl):
            for j in range(dr):
                coords = _class_coords(prod(left[i], right[j]), out_reps, img, dom)
                for k in range(dout):
                    consts[k][i][j] = coords[k]
        tables[(p, q)] = ProductTable(
            kind=kind, degrees=(p, q), out_degree=out_n, dims=(dl, dr, dout),
            constants=tuple(tuple(tuple(row) for row in mat) for mat in consts))
    return HHProducts(kind=kind, top=top, tables=tables,
                      engine="Chouhy-Solotar native diagonal",
                      basis=f"cs/{A.domain.name}", window=None,
                      references=_REFERENCES[kind] + ["chouhy_solotar"])
```

Implementation notes for the engineer (verify while wiring, all are existing
facts): `res.matrix(n, side)` returns list-of-rows over the Domain;
`native_cup` needs the resolution built to `p+q+1` (`comparison.py:548-551` —
`max_degree=top+2` covers every pair); `cs_hh_basis` builds its own resolution
internally — if profiling shows that doubling matters, pass-through reuse is a
follow-up, NOT this plan.

- [ ] **Step 4: Run test to verify it passes**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/resolutions_cs/test_products_cs.py -q`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/resolutions_cs/products.py tests/resolutions_cs/test_products_cs.py
git commit -m "feat(resolutions_cs): Domain-generic CS product tables via the native diagonal"
```

---

### Task 4: Induced Connes differentials

**Files:**
- Modify: `src/quiverlab/hochschild/products.py` (append)
- Create: `tests/hochschild/test_connes_b.py`

**Interfaces:**
- Consumes: GF(p): `engine.cyclic.connes_B_matrix(alg, n, basis_n, index_np1)`
  (numpy, shape `(dim C_{n+1}, dim C_n)`), `engine.tt_calculus.homology_classes
  (alg, n, p) -> _Quotient` (`.reps` int64 `(dimC, dim)`, `.coords(v)`,
  `.dim`); generic: `hochschild.cyclic.connes_B_matrix(A, n, max_cells) ->
  (rows, ncols, nrows)` over `A.domain`, `hochschild.bar.boundary_matrix(A, n,
  max_cells)`, `fields.linalg.nullspace/solve`.
- Produces: `connes_b_tables(A, top, max_cells) -> ConnesB`.

- [ ] **Step 1: Write the failing test**

```python
# tests/hochschild/test_connes_b.py
"""Induced Connes B on HH (Plan 35): B^2 = 0 at the induced level, rank
consistency with the (b,B) cyclic dims, and GF(p)/generic agreement."""
import pytest

import quiverlab as ql

pytestmark = [pytest.mark.oracle_selfcert]


@pytest.fixture(scope="module", params=[ql.GF(7), ql.QQ], ids=["GF7", "QQ"])
def A(request):
    return ql.TruncatedPolynomial(2, field=request.param)


def test_induced_B_squares_to_zero(A):
    from quiverlab.hochschild.products import connes_b_tables
    cb = connes_b_tables(A, 3, max_cells=4_000_000)
    # composite HH_n -> HH_{n+1} -> HH_{n+2} must vanish entry-wise
    dom = A.domain
    for n in range(2):
        M1, M2 = cb.matrices[n], cb.matrices[n + 1]
        rows, mid, cols = len(M2), len(M1), (len(M1[0]) if M1 else 0)
        for r in range(rows):
            for c in range(cols):
                acc = dom.zero()
                for k in range(mid):
                    acc = dom.add(acc, dom.mul(dom.coerce(M2[r][k]),
                                               dom.coerce(M1[k][c])))
                assert dom.is_zero(acc), f"(B∘B)[{r}][{c}] != 0 at n={n}"


def test_ranks_and_dims_recorded(A):
    from quiverlab.hochschild.products import connes_b_tables
    cb = connes_b_tables(A, 2, max_cells=4_000_000)
    hh = list(A.hochschild_homology(2, verbose=False).dims)
    assert cb.hh_dims == hh
    assert set(cb.matrices) == {0, 1} and set(cb.ranks) == {0, 1}
    for n in (0, 1):
        assert len(cb.matrices[n]) == hh[n + 1]
        assert 0 <= cb.ranks[n] <= min(hh[n], hh[n + 1])


def test_gfp_and_generic_ranks_agree():
    from quiverlab.hochschild.products import connes_b_tables
    A7 = ql.TruncatedPolynomial(3, field=ql.GF(32003))
    Aq = ql.TruncatedPolynomial(3, field=ql.QQ)
    r7 = connes_b_tables(A7, 3, max_cells=4_000_000).ranks
    rq = connes_b_tables(Aq, 3, max_cells=4_000_000).ranks
    assert r7 == rq          # char-0-shaped p: ranks agree (32003 is the big prime)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_connes_b.py -q`
Expected: FAIL — `ImportError: cannot import name 'connes_b_tables'`

- [ ] **Step 3: Implement** (append to `products.py`)

```python
def _generic_homology_quotient(A, n, max_cells):
    """(reps, image_cols) for HH_n over A.domain from the bar boundary matrices
    (the Domain-generic sibling of engine.tt_calculus.homology_classes)."""
    from quiverlab.fields.linalg import nullspace, rank as _rank
    from quiverlab.hochschild.bar import boundary_matrix
    dom = A.domain
    if n == 0:
        bnames = A.dim
        cycles = [[dom.one() if i == j else dom.zero() for j in range(bnames)]
                  for i in range(bnames)]
    else:
        M, _, _ = _as_matrix(boundary_matrix(A, n, max_cells))
        cycles = nullspace(M, dom)
    Mn1, _, _ = _as_matrix(boundary_matrix(A, n + 1, max_cells))
    image = [[Mn1[r][c] for r in range(len(Mn1))] for c in range(len(Mn1[0]))] \
        if Mn1 and Mn1[0] else []
    # greedy rank filter: keep cycles independent modulo the image
    reps, base = [], list(image)
    r0 = _rank(_cols_to_matrix(base), dom) if base else 0
    for v in cycles:
        rr = _rank(_cols_to_matrix(base + [v]), dom)
        if rr > r0:
            reps.append(v); base.append(v); r0 = rr
    return reps, image
```

(`_as_matrix` normalizes `boundary_matrix`'s `(rows, ncols, nrows)`-style return
— check the actual return of `hochschild.bar.boundary_matrix` at implementation
time: `bar.py:96`; `_cols_to_matrix` transposes a column list into row-major for
`fields.linalg.rank`. Both are 6-line private helpers in `products.py`; write
them with the real return shape in front of you and unit-test them via the
public tests above.)

```python
def connes_b_tables(A, top, max_cells=4_000_000):
    """Induced Connes differentials B: HH_n -> HH_{n+1}, 0 <= n < top."""
    from quiverlab.fields.primefield import PrimeField
    dom = A.domain
    if isinstance(dom, PrimeField):
        import numpy as np
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.cyclic import connes_B_matrix
        from quiverlab.engine.hh_engine import cn_basis
        from quiverlab.engine.tt_calculus import homology_classes
        p = dom.p
        E = to_engine(A.unit_adapted())
        H = {n: homology_classes(E, n, p) for n in range(top + 1)}
        matrices, ranks = {}, {}
        for n in range(top):
            idx = {g: i for i, g in enumerate(cn_basis(E, n + 1))}
            B = connes_B_matrix(E, n, cn_basis(E, n), idx)
            rows = []
            for i in range(H[n].dim):
                img = (B @ H[n].reps[:, i]) % p
                rows.append([int(x) for x in H[n + 1].coords(img)])
            # rows[i] = coords of B(e_i): store as matrix hh_{n+1} x hh_n
            matrices[n] = [[str(rows[i][k]) for i in range(H[n].dim)]
                           for k in range(H[n + 1].dim)]
            ranks[n] = _int_rank_mod_p(rows, H[n + 1].dim, p)
        hh = [H[n].dim for n in range(top + 1)]
        engine = f"engine (b,B) GF({p})"
    else:
        from quiverlab.fields.linalg import rank as _rank
        from quiverlab.hochschild.cyclic import connes_B_matrix
        AU = A.unit_adapted()
        quots = {n: _generic_homology_quotient(AU, n, max_cells)
                 for n in range(top + 1)}
        matrices, ranks = {}, {}
        for n in range(top):
            M, ncols, nrows = connes_B_matrix(AU, n, max_cells)
            reps_n, _ = quots[n]
            reps_n1, img_n1 = quots[n + 1]
            cols = []
            for v in reps_n:
                w = [dom.zero()] * nrows
                for r in range(nrows):
                    acc = dom.zero()
                    for c in range(ncols):
                        acc = dom.add(acc, dom.mul(M[r][c], dom.coerce(v[c])))
                    w[r] = acc
                cols.append(_class_coords_generic(w, reps_n1, img_n1, dom))
            matrices[n] = [[str(cols[i][k]) for i in range(len(reps_n))]
                           for k in range(len(reps_n1))]
            ranks[n] = _rank([[dom.coerce(cols[i][k]) for i in range(len(reps_n))]
                              for k in range(len(reps_n1))], dom) if reps_n else 0
        hh = [len(quots[n][0]) for n in range(top + 1)]
        engine = f"generic (b,B) mixed complex / {dom.name}"
    return ConnesB(top=top, hh_dims=hh, matrices=matrices, ranks=ranks,
                   engine=engine, references=_REFERENCES["connes_b"])
```

(`_class_coords_generic` is the same solve-in-`[image|reps]` as
`resolutions_cs/products.py::_class_coords` — put ONE shared copy in
`hochschild/products.py`, import it from `resolutions_cs/products.py`, returning
Domain elements here and strings there via a `stringify=` flag. `_int_rank_mod_p`
row-reduces a small int list mod p — reuse
`engine.coxeter.rref_mod_p` on the numpy array and count pivots.)

- [ ] **Step 4: Run test to verify it passes**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_connes_b.py -q`
Expected: 5 passed (3 tests, one parametrized x2)

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/hochschild/products.py tests/hochschild/test_connes_b.py
git commit -m "feat(hochschild): induced Connes differentials on HH (GF(p) + generic Domain)"
```

---

### Task 5: Public `Algebra` methods + routing + refusals

**Files:**
- Modify: `src/quiverlab/core/algebra.py` (append after `cyclic_homology`, ~line 520)
- Create: `tests/hochschild/test_products_api.py`

**Interfaces:**
- Consumes: `gfp_product_tables`, `cs_product_tables`, `connes_b_tables` (Tasks 2–4).
- Produces: `Algebra.cup_products(top, engine="auto", max_cells=4_000_000)`,
  `Algebra.cap_products(...)`, `Algebra.gerstenhaber_brackets(...)`,
  `Algebra.connes_differentials(top, max_cells=4_000_000)` — Tasks 8/9/11 call
  exactly these.

- [ ] **Step 1: Write the failing test**

```python
# tests/hochschild/test_products_api.py
"""Public product surface: routing, refusal matrix, provenance (Plan 35 §2)."""
import pytest

import quiverlab as ql
from quiverlab.errors import QuiverlabError

pytestmark = [pytest.mark.oracle_selfcert]


def test_gfp_routes_to_bar():
    A = ql.TruncatedPolynomial(2, field=ql.GF(7))
    hp = A.cup_products(2)
    assert hp.basis == "bar/GF(7)"


def test_quiver_presented_qq_routes_to_cs():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=ql.QQ)
    hp = A.cup_products(2)
    assert hp.basis.startswith("cs/")


def test_structure_constants_off_gfp_refuse():
    A = ql.TruncatedPolynomial(2, field=ql.QQ)   # presentation-less over QQ?
    # TruncatedPolynomial IS quiver-presented; build a genuinely presentation-less
    # algebra from structure constants:
    B = ql.Algebra.from_structure_constants(
        [[[1, 0], [0, 1]], [[0, 1], [0, 0]]], field=ql.QQ)  # k[x]/(x^2) as SC
    with pytest.raises(QuiverlabError):
        B.cup_products(2)


def test_bracket_refuses_off_gfp():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=ql.QQ)
    with pytest.raises(QuiverlabError):
        A.gerstenhaber_brackets(2)


def test_unknown_engine_refuses():
    A = ql.TruncatedPolynomial(2, field=ql.GF(7))
    with pytest.raises(QuiverlabError):
        A.cup_products(2, engine="fast")     # products know auto/bar/cs only


def test_explicit_cs_on_gfp_serves_cs_basis():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=ql.GF(7))
    hp = A.cup_products(2, engine="cs")
    assert hp.basis == "cs/GF(7)"


def test_connes_serves_both_fields():
    for F in (ql.GF(7), ql.QQ):
        A = ql.TruncatedPolynomial(2, field=F)
        cb = A.connes_differentials(2)
        assert set(cb.matrices) == {0, 1}
```

Check `Algebra.from_structure_constants`'s real name before running (grep
`from_structure_constants\|structure_constants` in `core/algebra.py`; the
Plan-19 refusal tests already construct such an algebra — copy their
construction verbatim if the name differs).

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_products_api.py -q`
Expected: FAIL — `AttributeError: 'Algebra' object has no attribute 'cup_products'`

- [ ] **Step 3: Implement the four methods** (in `core/algebra.py`, after
`cyclic_homology`; follow the `hochschild_cohomology` docstring style)

```python
    def _product_dispatch(self, kind, top, engine, max_cells):
        """Shared Plan-35 routing: GF(p) -> bar/tt tables; quiver-presented ->
        CS native (any Domain; also the DepthLimitError fallback target);
        presentation-less off GF(p) -> loud refusal. One engine end-to-end."""
        from quiverlab.errors import DepthLimitError
        from quiverlab.fields.primefield import PrimeField
        from quiverlab.hochschild.products import gfp_product_tables
        if engine not in ("auto", "bar", "cs"):
            raise QuiverlabError(
                f"unknown engine {engine!r} for {kind} tables",
                hint="choose 'auto', 'bar', or 'cs'")
        is_gfp = isinstance(self.domain, PrimeField)
        presented = self.quiver is not None and self.relations is not None
        if engine == "cs" or (engine == "auto" and not is_gfp):
            if kind == "bracket":
                raise QuiverlabError(
                    "the Gerstenhaber bracket is served over GF(p) only "
                    "(bar window; no CS-native brace machinery in v1)",
                    hint="construct the algebra over GF(p)")
            if not presented:
                raise QuiverlabError(
                    f"{kind} tables off GF(p) need a quiver presentation "
                    "(the CS route); this algebra has structure constants only",
                    hint="build the algebra via Quiver.algebra, or use GF(p)")
            from quiverlab.resolutions_cs.products import cs_product_tables
            return cs_product_tables(self, kind, top, max_cells)
        if not is_gfp:            # engine == "bar" explicitly, off GF(p)
            raise QuiverlabError(
                f"engine='bar' {kind} tables need GF(p) (the tt facade)",
                hint="use engine='auto' (routes CS for presented algebras)")
        try:
            return gfp_product_tables(self, kind, top, max_cells)
        except DepthLimitError:
            if engine != "auto" or not presented or kind == "bracket":
                raise
            from quiverlab.resolutions_cs.products import cs_product_tables
            return cs_product_tables(self, kind, top, max_cells)

    def cup_products(self, top, engine="auto", max_cells=4_000_000):
        """Structure-constant tables of the cup product HH^p (x) HH^q ->
        HH^{p+q} for every p+q <= top, on the recorded basis. Exact. engine:
        'auto' (GF(p) -> bar/tt, else CS for presented algebras, with the CS
        depth fallback), 'bar' (GF(p) tt facade, loud otherwise), 'cs'
        (Chouhy-Solotar native diagonal, presented algebras, any Domain)."""
        return self._product_dispatch("cup", top, engine, max_cells)

    def cap_products(self, top, engine="auto", max_cells=4_000_000):
        """Structure-constant tables of the cap action HH^p (x) HH_n ->
        HH_{n-p} for p <= n <= top. Same engine semantics as cup_products."""
        return self._product_dispatch("cap", top, engine, max_cells)

    def gerstenhaber_brackets(self, top, engine="auto", max_cells=4_000_000):
        """Structure-constant tables of the Gerstenhaber bracket HH^p (x)
        HH^q -> HH^{p+q-1} for pairs p, q >= 1 with p+q-1 <= top. GF(p) only
        and window-bounded (the result records the served window); the
        degree-0 insertion action is out of scope."""
        return self._product_dispatch("bracket", top, engine, max_cells)

    def connes_differentials(self, top, max_cells=4_000_000):
        """Induced Connes differentials B : HH_n -> HH_{n+1} (matrices +
        ranks) for 0 <= n < top. GF(p) via the engine (b,B); any other exact
        Domain via the generic mixed complex — no engine choice to make."""
        from quiverlab.hochschild.products import connes_b_tables
        return connes_b_tables(self, top, max_cells=max_cells)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_products_api.py -q`
Expected: 7 passed

- [ ] **Step 5: Run the float gate and fast bucket sanity**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/test_no_floats.py tests/hochschild/ -q`
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add src/quiverlab/core/algebra.py tests/hochschild/test_products_api.py
git commit -m "feat(core): public cup/cap/bracket/Connes product surface with engine routing"
```

---

### Task 6: Identity + cross-engine oracle batteries

**Files:**
- Create: `tests/hochschild/test_products_identities.py`

**Interfaces:** consumes only the public API (Task 5).

- [ ] **Step 1: Write the battery** (these MUST pass immediately — they are
oracles for already-committed code; any failure is a Task-2/3/4 bug: STOP and
fix there)

```python
# tests/hochschild/test_products_identities.py
"""Plan 35 identity batteries. Table-level: the identities of a Gerstenhaber
algebra + the cap module structure + B^2=0, on small algebras over several
primes, plus bar<->CS table equality in-window (the cross-engine gate)."""
import itertools
import pytest

import quiverlab as ql

PRIMES = (32003, 2, 3, 5)


def _mat(table):
    """constants as int tensor [k][i][j]."""
    return [[[int(c) for c in row] for row in mat] for mat in table.constants]


def _compose(hp, p, q, r, prime):
    """(f_i ∪ g_j) ∪ h_l and f_i ∪ (g_j ∪ h_l) as coordinate tensors — equality
    is associativity at the table level."""
    Tpq, Tqr = hp.tables[(p, q)], hp.tables[(q, r)]
    Tpq_r, Tp_qr = hp.tables[(p + q, r)], hp.tables[(p, q + r)]
    dl, dm, dr_ = Tpq.dims[0], Tpq.dims[1], Tqr.dims[1]
    dout = Tpq_r.dims[2]
    left = [[[[0] * dout for _ in range(dr_)] for _ in range(dm)] for _ in range(dl)]
    right = [[[[0] * dout for _ in range(dr_)] for _ in range(dm)] for _ in range(dl)]
    A_, B_, C_, D_ = _mat(Tpq), _mat(Tpq_r), _mat(Tqr), _mat(Tp_qr)
    for i, j, l in itertools.product(range(dl), range(dm), range(dr_)):
        for k in range(dout):
            left[i][j][l][k] = sum(A_[m][i][j] * B_[k][m][l]
                                   for m in range(Tpq.dims[2])) % prime
            right[i][j][l][k] = sum(C_[m][j][l] * D_[k][i][m]
                                    for m in range(Tqr.dims[2])) % prime
    return left, right


@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("prime", PRIMES)
def test_cup_graded_commutative_and_associative(prime):
    A = ql.TruncatedPolynomial(3, field=ql.GF(prime))
    hp = A.cup_products(3)
    for (p, q), t in hp.tables.items():
        if (q, p) not in hp.tables:
            continue
        s = hp.tables[(q, p)]
        sign = 1 if (p * q) % 2 == 0 else prime - 1
        M, N = _mat(t), _mat(s)
        for k in range(t.dims[2]):
            for i in range(t.dims[0]):
                for j in range(t.dims[1]):
                    assert M[k][i][j] % prime == (sign * N[k][j][i]) % prime
    for p, q, r in [(0, 0, 1), (0, 1, 1), (1, 1, 1)]:
        left, right = _compose(hp, p, q, r, prime)
        assert left == right, f"associativity fails at {(p, q, r)}"


@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("prime", PRIMES)
def test_bracket_antisymmetry(prime):
    A = ql.TruncatedPolynomial(3, field=ql.GF(prime))
    hb = A.gerstenhaber_brackets(3)
    for (p, q), t in hb.tables.items():
        if (q, p) not in hb.tables:
            continue
        s, M, N = hb.tables[(q, p)], _mat(t), _mat(hb.tables[(q, p)])
        sign = 1 if ((p - 1) * (q - 1)) % 2 == 0 else prime - 1
        for k in range(t.dims[2]):
            for i in range(t.dims[0]):
                for j in range(t.dims[1]):
                    assert M[k][i][j] % prime == (-sign * N[k][j][i]) % prime, \
                        f"[f,g] != -(-1)^((p-1)(q-1))[g,f] at {(p, q)}"


@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("prime", (7, 3))
def test_bar_vs_cs_cup_tables_in_window(prime):
    """The cross-engine gate: same DIMS and RANK-equivalent tables. Bases
    differ, so we compare basis-independent data: dims and, for each (p,q),
    the RANK of the flattened constants matrix (dl*dr x dout) mod p."""
    import numpy as np
    from quiverlab.engine.coxeter import rref_mod_p
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=ql.GF(prime))
    bar = A.cup_products(2, engine="bar")
    cs = A.cup_products(2, engine="cs")
    assert sorted(bar.tables) == sorted(cs.tables)
    for key in bar.tables:
        tb, tc = bar.tables[key], cs.tables[key]
        assert tb.dims == tc.dims
        def flat_rank(t):
            dl, dr, dout = t.dims
            if 0 in t.dims:
                return 0
            M = np.array([[int(t.constants[k][i][j]) for k in range(dout)]
                          for i in range(dl) for j in range(dr)], dtype=np.int64)
            _, piv = rref_mod_p(M % prime, prime)
            return len(piv)
        assert flat_rank(tb) == flat_rank(tc), f"table rank differs at {key}"
```

Also add (same file) the cap-module identity `(z∩f)∩g = z∩(f∪g)` at table
level over GF(7) on `TruncatedPolynomial(2)` — the composition helper mirrors
`_compose` with one cap and one cup table; and cup-Leibniz for the bracket
(`[f, g∪h] = [f,g]∪h + (-1)^{(p-1)q} g∪[f,h]`) on the (1,1,1) triple. Write
them as two additional `oracle_selfcert` tests with the same `_mat` helper —
the formulas are exactly the two displayed equations, composed mod p.

- [ ] **Step 2: Run**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_products_identities.py -q`
Expected: all pass on first run (oracles over committed code). A failure is a
REAL BUG in Tasks 2–4: do not adjust the test; fix the builder.

- [ ] **Step 3: Commit**

```bash
git add tests/hochschild/test_products_identities.py
git commit -m "test(hochschild): Gerstenhaber/cap/B identity + cross-engine table batteries"
```

---

### Task 7: Literature battery + QPA probe

**Files:**
- Create: `tests/hochschild/test_products_literature.py`
- Create: `tests/qpa/test_products_qpa.py`

- [ ] **Step 1: Literature pins**

```python
# tests/hochschild/test_products_literature.py
"""Plan 35 literature pins. k[x]/(x^2): HH^* in char 0 is k[u]/(u^2) x k[e]
shaped -- concretely over GF(32003) (char-0-shaped): dim HH^n = 1 for all n,
even-degree classes multiply freely (cup of generators in degrees 2a, 2b is a
unit multiple of the degree 2a+2b generator), odd*odd lands in even with the
char-0 vanishing u^2 = 0. QuantumCI q=1 over GF(2) is the (2,2) commutative CI
of the BGMS pin [4,8,12,16]."""
import pytest

import quiverlab as ql

pytestmark = [pytest.mark.oracle_literature]


def _entry(t, k, i, j):
    return int(t.constants[k][i][j])


def test_dual_numbers_cup_ring_char0_shape():
    A = ql.TruncatedPolynomial(2, field=ql.GF(32003))
    hp = A.cup_products(4)
    dims = {n: hp.tables[(0, n)].dims[1] for n in range(5)}
    assert dims == {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
    # even generators compose to a nonzero multiple of the even basis class
    assert _entry(hp.tables[(2, 2)], 0, 0, 0) != 0
    # odd squares vanish in char != 2 (graded commutativity forces 2x = 0)
    assert _entry(hp.tables[(1, 1)], 0, 0, 0) == 0
    assert _entry(hp.tables[(3, 1)], 0, 0, 0) == 0


def test_dual_numbers_char2_odd_squares_survive():
    A = ql.TruncatedPolynomial(2, field=ql.GF(2))
    hp = A.cup_products(2)
    # char 2: HH^* of k[x]/(x^2) is the polynomial-like ring where the odd
    # generator squares NONZERO (the classical char-2 phenomenon)
    assert _entry(hp.tables[(1, 1)], 0, 0, 0) != 0


def test_qci_dims_line_up_with_bgms():
    A = ql.QuantumCI(field=ql.GF(2), q=1)
    hp = A.cup_products(2)
    assert [hp.tables[(0, n)].dims[1] for n in range(3)] == [4, 8, 12]


def test_connes_b_rank_dual_numbers():
    A = ql.TruncatedPolynomial(2, field=ql.GF(32003))
    cb = A.connes_differentials(4)
    # HH_n(k[x]/(x^2)) has dim 1 in every degree (char-0 shape); Connes B
    # alternates iso/zero along the SBI pattern: rank B_0 = 1 here.
    assert cb.hh_dims == [1, 1, 1, 1, 1]
    assert cb.ranks[0] == 1
```

**Verify the char-2 and B-rank pins against the engine BEFORE freezing** — run
each computation interactively first; if a stated value disagrees with the
engine, the ENGINE (already oracle-gated by `tests/test_tt_calculus.py` /
`test_gerstenhaber.py` pins) wins: correct the test's expected value and cite
the actual computed number in a comment (the CRS-2004 precedent — pin the
verified value, note the source discrepancy).

- [ ] **Step 2: QPA probe** — `tests/qpa/test_products_qpa.py` with the `qpa`
marker: session-probe for any of `CupProduct`, `HochschildCohomologyRing`,
GAP-level multiplication of `ExtAlgebra` generators; if none exists (expected),
`pytest.skip("QPA 1.37 exposes no Hochschild product surface")` and the
verification page (Task 13) records the honest-scope entry naming the Task-6
identity batteries as the covering oracle. Follow the probe pattern of
`tests/qpa/test_left_modules_qpa.py` (session fixture + `libgap.eval`
single-statements).

- [ ] **Step 3: Run + commit**

```bash
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/test_products_literature.py -q
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/qpa/test_products_qpa.py -q -m qpa
git add tests/hochschild/test_products_literature.py tests/qpa/test_products_qpa.py
git commit -m "test(hochschild): literature pins for the product surface + QPA probe"
```

---

### Task 8: `hpc/spec.py` dispatch + schema + estimator

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch`, `_snip` at :1522, the
  ALGEBRA-kind acceptance in `parse_compute_item` consumers)
- Modify: `webapp/server/schema.py` (no new validation needed — `parse_compute_item`
  accepts any `[a-z_]+` kind; the runner refuses unknowns — but ADD the four
  kinds to any explicit whitelist if one exists; grep `unsupported computation`)
- Modify: `webapp/server/estimator.py::sizing_dim` (products size like
  `hh_cohomology` at the same top)
- Create: `tests/hpc/test_products_kinds.py`
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE new golden entry; document the addition in the docstring re-freeze list)

**Interfaces:**
- Consumes: `A.cup_products(top)` etc. (Task 5), `.blocks()` (Task 1).
- Produces: result blocks under keys `"cup"`, `"cap"`, `"bracket"`,
  `"connes_b"` in `results`, each `= method(top).blocks() + {"citations": _citation_pairs(refs)}`.
  Tasks 9/10/11 rely on these exact keys and the block fields
  `kind/top/engine/basis/tables/window/matrices/ranks/hh_dims/references/citations`.

- [ ] **Step 1: Write the failing test**

```python
# tests/hpc/test_products_kinds.py
"""The four product kinds through the spec runner (Plan 35)."""
import pytest

from quiverlab.hpc.spec import run as spec_run, parse_config


def _req(compute):
    return parse_config({
        "schema": 2,
        "algebra": {"kind": "quiver", "vertices": [1],
                    "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                    "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": compute,
        "artifacts": {"pdf": False, "tikz": False}})


def test_cup_block(tmp_path):
    res = spec_run(_req(["cup:0..2"]), tmp_path)
    b = res["results"]["cup"]
    assert b["kind"] == "cup" and b["top"] == 2
    assert b["basis"].startswith(("bar/", "cs/"))
    assert b["citations"] and b["references"]
    degs = [t["degrees"] for t in b["tables"]]
    assert [0, 0] in degs and [1, 1] in degs


def test_all_four_kinds_serve(tmp_path):
    res = spec_run(_req(["cup:0..2", "cap:0..2", "bracket:0..2",
                         "connes_b:0..2"]), tmp_path)
    assert set(res["results"]) == {"cup", "cap", "bracket", "connes_b"}
    assert res["results"]["bracket"]["window"] == 2
    assert res["results"]["connes_b"]["ranks"].keys() == {"0", "1"}


def test_range_required(tmp_path):
    from quiverlab.hpc.spec import ComputeError
    with pytest.raises(Exception):
        spec_run(_req(["cup"]), tmp_path)     # ComputeError -> typed 4xx upstream
```

(Adjust `parse_config`/`spec_run` names to the module's real entry points —
`grep -n "^def run\|def parse" src/quiverlab/hpc/spec.py` — the webapp tests
`tests/webapp/test_runner_delegation.py` show the canonical calling shape via
`webapp.server.runner.run_spec`; using THAT wrapper here is equally fine and
tests the same dispatch.)

- [ ] **Step 2: Run to verify failure** — expected `ComputeError: unsupported
computation 'cup'`.

- [ ] **Step 3: Implement in `_dispatch`** (before the final `raise ComputeError`):

```python
    if kind in ("cup", "cap", "bracket", "connes_b"):
        top = item.hi
        if top is None:
            raise ComputeError("SchemaError",
                               f"{kind} needs a degree range, e.g. '{kind}:0..4'")
        method = {"cup": A.cup_products, "cap": A.cap_products,
                  "bracket": A.gerstenhaber_brackets,
                  "connes_b": A.connes_differentials}[kind]
        obj = method(top)
        block = obj.blocks()
        keys = list(block["references"])
        block["citations"] = _citation_pairs(keys)
        return block, None
```

Add snippet lines to `_snip` (spec.py:1522 area):

```python
        "cup": lambda it: f"A.cup_products({it.hi})",
        "cap": lambda it: f"A.cap_products({it.hi})",
        "bracket": lambda it: f"A.gerstenhaber_brackets({it.hi})",
        "connes_b": lambda it: f"A.connes_differentials({it.hi})",
```

Estimator (`webapp/server/estimator.py`): in `sizing_dim`/wherever
`hh_cohomology` tops feed the ops estimate, treat the four kinds identically to
`hh_cohomology` with the same `hi` (grep `hh_cohomology` there; extend the kind
set in place).

- [ ] **Step 4: Run tests + freeze ONE new runner golden**

Run the new test; then extend `tests/webapp/_runner_goldens.json` with one
entry `products_loop_gf2` — body: the `_req(["cup:0..2", "connes_b:0..2"])`
request above; freeze `result_json` + `canonical_key` by running
`webapp.server.runner.run_spec` once and dumping
`json.dumps(result, sort_keys=True, default=str)` (the exact recipe in
`test_runner_delegation.py`); add the dated entry to the docstring re-freeze
list ("2026-08-01 (products_loop_gf2): ADDED for Plan 35 — new kinds, existing
six entries untouched").

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hpc/test_products_kinds.py tests/webapp/test_runner_delegation.py -q`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/hpc/spec.py webapp/server/estimator.py tests/hpc/test_products_kinds.py tests/webapp/_runner_goldens.json tests/webapp/test_runner_delegation.py
git commit -m "feat(hpc,webapp): serve cup/cap/bracket/connes_b compute kinds"
```

---

### Task 9: Pyodide runner twin

**Files:**
- Modify: `docs/gui/runner.py` (`compute_one` at :532, the snippet map at :686,
  and the block-postprocess region at :804)
- Create: `tests/gui/test_products_runner_twin.py`

**Interfaces:** consumes Task 5 methods; produces the SAME block dicts as Task 8
(the cross-runner contract test asserts key-for-key equality).

- [ ] **Step 1: Write the failing test**

```python
# tests/gui/test_products_runner_twin.py
"""docs/gui/runner.py serves the product kinds SHAPE-IDENTICALLY to
quiverlab.hpc.spec (the Plan-26 twin contract, extended to Plan 35)."""
import json
import sys
import types


def _load_gui_runner():
    # the runner twin is import-loaded the way tests/gui does it today --
    # copy the loader helper from the existing tests/gui/test_runner.py.
    import importlib.util, pathlib
    p = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    spec = importlib.util.spec_from_file_location("gui_runner", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_cup_block_shape_matches_spec_runner(tmp_path):
    gui = _load_gui_runner()
    body = {"schema": 2,
            "algebra": {"kind": "quiver", "vertices": [1],
                        "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                        "field": {"kind": "GF", "p": 2, "n": 1}},
            "compute": ["cup:0..2", "connes_b:0..2"],
            "artifacts": {"pdf": False, "tikz": False}}
    out = json.loads(gui.run_request(json.dumps(body)))
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    assert out["results"]["cup"] == ref["results"]["cup"]
    assert out["results"]["connes_b"] == ref["results"]["connes_b"]
```

(Adapt `run_request`/loader names from the existing `tests/gui/` runner tests —
copy their import scaffold verbatim; the assertion is the contract.)

- [ ] **Step 2–4: Fail → implement → pass.** Implementation in `compute_one`:
the same 12-line branch as Task 8 Step 3 (method map + `.blocks()` +
citations), and the snippet map gains the same four lambda lines. Where the
runner postprocess at :804 special-cases `hh_*`, the product kinds need no
special case (their blocks are complete).

- [ ] **Step 5: Commit**

```bash
git add docs/gui/runner.py tests/gui/test_products_runner_twin.py
git commit -m "feat(gui): product kinds in the Pyodide runner twin"
```

---

### Task 10: GUI pick-list + rendering + i18n

**Files:**
- Modify: `docs/gui/gui.js` AND `webapp/static/gui/gui.js` (BYTE-IDENTICAL —
  edit one, copy over the other, the existing vendoring test gates it)
- Modify: `webapp/static/app.js` + `webapp/server/i18n/en.json` + `es.json`
  (families-page invariants list)
- Create: `tests/gui/test_products_gui_wiring.py` (static assertions on the JS
  sources — the pattern of existing `tests/gui` string checks)

- [ ] **Step 1: gui.js** — four new checkboxes with degree pickers, inserted
directly after the `hhh` block (gui.js:65 area, ids `qlgui-cup`, `qlgui-cap`,
`qlgui-bracket`, `qlgui-connes_b`, each with a `-top` numeric input defaulting
to 2); register ids in the `el` list (:132-134); compose the request **in this
exact order, immediately after `hh_homology`** (gui.js:631 — this order is the
Task-12 curated-request order rule):

```javascript
    if (el.cup.checked) compute.push("cup:0.." + el["cup-top"].value);
    if (el.cap.checked) compute.push("cap:0.." + el["cap-top"].value);
    if (el.bracket.checked) compute.push("bracket:0.." + el["bracket-top"].value);
    if (el.connes_b.checked) compute.push("connes_b:0.." + el["connes_b-top"].value);
```

Rendering (the block renderer switch at gui.js:1157 area): for
`cup/cap/bracket` render one sub-section per table — heading
`HH^p ∪ HH^q → HH^{p+q}` (cap: `HH^p ∩ HH_n → HH_{n-p}`; bracket:
`[HH^p, HH^q] → HH^{p+q-1}`), then the equation list: for each `(i, j)` the
line `e_i ∙ e_j = Σ_k c·f_k` built from `constants[k][i][j]`, skipping zero
terms, `matIsZero`-style single line "all products vanish" when the whole
table is zero. `bracket.window` renders the one-line note "served to degree
window N (bar-transport bound)". For `connes_b` render per-`n` `matrixGrid`
of `matrices[n]` + "rank B_n = r" line.

- [ ] **Step 2: app.js families page** — add the four to the invariants
checkbox list with `inv.cup` etc.; i18n keys:

```json
  "inv.cup": "Cup products (ring structure) 0..N",
  "inv.cap": "Cap products (module structure) 0..N",
  "inv.bracket": "Gerstenhaber brackets 0..N",
  "inv.connes_b": "Connes differentials 0..N",
  "block.cup.title": "Cup product tables",
  "block.cap.title": "Cap product tables",
  "block.bracket.title": "Gerstenhaber bracket tables",
  "block.connes_b.title": "Connes differentials",
  "block.bracket.window": "served to degree window {n} (bar-transport bound)",
  "block.products.zero": "all products vanish in this degree pair"
```

with the ES twins (`"inv.cup": "Productos cup (estructura de anillo) 0..N"`,
`"inv.cap": "Productos cap (estructura de módulo) 0..N"`, `"inv.bracket":
"Corchetes de Gerstenhaber 0..N"`, `"inv.connes_b": "Diferenciales de Connes
0..N"`, `"block.bracket.window": "calculado hasta la ventana de grado {n}
(cota del transporte bar)"`, `"block.products.zero": "todos los productos se
anulan en este par de grados"`, titles: `"Tablas de productos cup"`, `"Tablas
de productos cap"`, `"Tablas de corchetes de Gerstenhaber"`, `"Diferenciales
de Connes"`).

- [ ] **Step 3: the static test** asserts: both gui.js copies byte-equal
(existing gate covers it — run it), the four `compute.push` lines exist in
order after the `hh_homology` push, all four i18n keys present in BOTH locale
files (load the JSONs, assert key sets).

- [ ] **Step 4: Run** `tests/gui/ tests/webapp/ -q` — all pass. **Manually
smoke** the offline GUI: `.venv/bin/quiverlab-hpc gui --data-dir /tmp/ql-p35`
→ draw the `x^3` loop, tick Cup + Connes, Compute, verify the tables render.

- [ ] **Step 5: Commit**

```bash
git add docs/gui/gui.js webapp/static/gui/gui.js webapp/static/app.js webapp/server/i18n/en.json webapp/server/i18n/es.json tests/gui/test_products_gui_wiring.py
git commit -m "feat(gui): product pick-list, table rendering, EN/ES labels"
```

---

### Task 11: Worked-steps products chapter

**Files:**
- Modify: `src/quiverlab/trace/events.py` (one new event), `src/quiverlab/trace/render_html.py`
  + `src/quiverlab/trace/results_html.py` (chapter renderer), `src/quiverlab/hpc/spec.py`
  (trace-selection extension at :750-845)
- Create: `src/quiverlab/trace/products.py` (narration builder)
- Create: `tests/trace/test_products_chapter.py`

**Interfaces:**
- Produces: `trace/products.py::products_chapter(A, kind, obj) -> list[event]`
  where `obj` is the Task-1 result object; a new frozen dataclass
  `trace/events.py::ProductStep(kind, degrees, heading, lines, matrix, note)`
  (matrix = list-of-rows of str or None); renderer branch in `render_html` that
  turns `ProductStep` into the definitional preamble + `matrix_grid` output.
  The ALL_EVENTS gate in the renderer must learn the new event type (it REFUSES
  foreign stream objects — extending the union IS the task).

- [ ] **Step 1: Write the failing test**

```python
# tests/trace/test_products_chapter.py
"""The products chapter: events render, dims drift-gate, full record (Plan 35)."""
import pytest

import quiverlab as ql

pytestmark = [pytest.mark.oracle_selfcert]


def test_chapter_events_render_and_carry_tables(tmp_path):
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    A = ql.TruncatedPolynomial(2, field=ql.GF(7))
    hp = A.cup_products(2)
    events = products_chapter(A, "cup", hp)
    assert events, "chapter must not be empty"
    html = render_html(events, table=None, algebra=A, kind="cup", top=2)
    assert "cup product" in html.lower()
    assert "HH" in html
    # every nonzero constant appears in the page
    for t in hp.tables.values():
        for mat in t.constants:
            for row in mat:
                for c in row:
                    if c != "0":
                        assert c in html


def test_drift_gate_fires_on_dim_mismatch():
    from quiverlab.trace.products import products_chapter
    A = ql.TruncatedPolynomial(2, field=ql.GF(7))
    hp = A.cup_products(1)
    # sabotage a dim -- the builder must refuse to narrate inconsistent data
    hp.tables[(0, 0)] = hp.tables[(0, 0)].__class__(
        kind="cup", degrees=(0, 0), out_degree=0, dims=(99, 1, 1),
        constants=hp.tables[(0, 0)].constants)
    with pytest.raises(Exception):
        products_chapter(A, "cup", hp)
```

(`render_html`'s real signature: mirror how `tests/trace/` call it today —
copy the call shape from `tests/trace/test_report_completeness_m0729.py`.)

- [ ] **Step 2–4: Fail → implement → pass.** `products_chapter` emits, per
kind: one `StepNote` with the definition text (cup: the Gerstenhaber-1963
formula `(f ∪ g)(a_1..a_{p+q}) = f(a_1..a_p)·g(a_{p+1}..a_{p+q})` in the TeX
the renderer already typesets; bracket: circle + graded commutator; cap: the
`a_0·f(a_1..a_p) ⊗ ...` collapse; B: the cyclic-rotation sum), one `ResultDims`
carrying the per-degree HH dims (the drift gate: assert
`table.dims` match the recorded `HH` dims — raise `QuiverlabError` otherwise),
then one `ProductStep` per table with the equation lines and (for `connes_b`)
the matrix. `spec.run`: extend the trace-selection so when NO hh/module trace
was chosen and a product kind is present with `req.artifacts.pdf`, the first
product kind backs the worked-steps bundle: mirror the `module_trace` pattern
at spec.py:750-758 (`product_trace = (kind, item.hi, obj)`) and add a
`_write_product_worked_steps` sibling of `_write_module_worked_steps` calling
`products_chapter` and `write_trace` — copy that function's structure line for
line, swapping the events source.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/trace/ src/quiverlab/hpc/spec.py tests/trace/test_products_chapter.py
git commit -m "feat(trace): worked-steps products chapter with drift gates"
```

---

### Task 12: Curated cache examples gain the products

**Files:**
- Modify: `webapp/precomputed/examples/*/request.json` (6 files) + `manifest.yaml`
  comment + `container/examples/*.yaml` (the 6 CLI twins — keep them in sync)
- Create: `tests/webapp/test_curated_reachability.py`
- Regenerate: `webapp/precomputed/examples/*/{result.json,tikz.tex}`

- [ ] **Step 1: The reachability test FIRST** (it pins the order rule):

```python
# tests/webapp/test_curated_reachability.py
"""Every curated request must be reachable from the GUI compose path: kinds in
the GUI's push order (hh_cohomology, hh_homology, cup, cap, bracket, connes_b,
cartan, coxeter_polynomial, global_dimension, center, then module kinds), and
its canonical key equal to the seed a fresh seeding run would write (Plan 35 §5)."""
import json
import pathlib

import pytest

from webapp.server.cache import canonical_key, library_version
from webapp.server.schema import ComputeRequest, parse_compute_item

_EX = pathlib.Path(__file__).resolve().parents[2] / "webapp" / "precomputed" / "examples"
_GUI_ORDER = ["hh_cohomology", "hh_homology", "cup", "cap", "bracket",
              "connes_b", "cartan", "coxeter_polynomial", "global_dimension",
              "center", "dimension", "dimension_vector", "rad_top_soc", "tau",
              "tau_minus", "projective_resolution", "injective_resolution",
              "projective_dimension", "injective_dimension", "decompose",
              "ext", "tor"]


@pytest.mark.parametrize("bundle", sorted(p.name for p in _EX.iterdir() if p.is_dir()))
def test_compute_list_is_in_gui_order(bundle):
    body = json.loads((_EX / bundle / "request.json").read_text(encoding="utf-8"))
    kinds = [parse_compute_item(s).kind for s in body["compute"]]
    order = {k: i for i, k in enumerate(_GUI_ORDER)}
    assert kinds == sorted(kinds, key=order.__getitem__), \
        f"{bundle}: compute list must follow the GUI compose order"
    assert {"cup", "cap", "connes_b"} <= set(kinds), \
        f"{bundle}: Plan 35 -- the curated examples carry the products"


@pytest.mark.parametrize("bundle", sorted(p.name for p in _EX.iterdir() if p.is_dir()))
def test_canonical_key_matches_model_roundtrip(bundle):
    body = json.loads((_EX / bundle / "request.json").read_text(encoding="utf-8"))
    req = ComputeRequest.model_validate(body)
    assert canonical_key(req.model_dump(by_alias=True), library_version()) == \
        canonical_key(body, library_version()), \
        f"{bundle}: raw request and model round-trip must produce ONE key"
```

Cross-check `_GUI_ORDER` against the REAL gui.js push order after Task 10
(read gui.js:629-651 in the final tree; the module kinds' relative order comes
from the module-panel push block — adjust the list to the source of truth,
which is gui.js, and leave a comment pointing at the line range).

- [ ] **Step 2: Edit the six `request.json`** — insert, right after the
`hh_homology` entry (or after `hh_cohomology` where no homology is requested),
the four kinds with tops equal to that example's HH top T (`"cup:0..T"`,
`"cap:0..T"`, `"bracket:0..T"`, `"connes_b:0..T"`). `bracket` only where the
field is GF(p) (all six are GF(p) — verify each `field` block). Keep the six
`container/examples/*.yaml` twins in sync (same insertion in their `compute`
lists). Run the reachability test — green.

- [ ] **Step 3: Recompute the bundles** — reuse the 2026-07-31 release-prep
scripts (scratchpad `regen_remaining.py` / `replace_bundles_v3.py` pattern,
committed nowhere — rewrite them inline in `/tmp` scratchpad): compute each
new `request.json` through `webapp.server.runner.run_spec` (deep ones under
`nohup … & disown` + Monitor — kz24-deep exceeded 90 min BEFORE products;
expect longer; if a deep example's product degrees are infeasible, TRIM that
example's product tops (e.g. `cup:0..4`) and note it in `manifest.yaml` beside
the entry — the spec's stated fallback). Gate: every PRE-EXISTING result block
byte-identical; the ONLY additions are the four product blocks; `tikz.tex`
byte-identical. Replace `result.json`; update the `manifest.yaml` comment to
name the product blocks.

- [ ] **Step 4: Run the webapp + seed suites**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/webapp/ -q`
Expected: all pass (the six existing runner goldens are UNTOUCHED — they pin
runner behavior, not bundle bytes).

- [ ] **Step 5: Commit**

```bash
git add webapp/precomputed/ container/examples/ tests/webapp/test_curated_reachability.py
git commit -m "feat(webapp): curated examples carry the product blocks (keys move, reachability gated)"
```

---

### Task 13: Verification page, docs, marker audit, final gates

**Files:**
- Modify: `docs/verification.md` (products row: subsystem → oracles → tests;
  audited class counts re-run; honest-scope entry: "QPA exposes no Hochschild
  product surface — covered by the identity batteries (graded commutativity,
  associativity, Jacobi, cup-Leibniz, cap module law, B²=0) and the
  k[x]/(x²)/QuantumCI literature pins"; bracket GF(p)/window scope named)
- Modify: `docs/internals/` — add a products section to the natural chapter
  (grep for the TT-calculus mention; likely `docs/internals/` HH chapter) with
  the four public signatures + the basis-provenance caveat
- Modify: `README.md` — one line in the capabilities list ("cup/cap products,
  Gerstenhaber bracket, Connes differentials — exact, with worked-steps
  reports")
- Modify: CLAUDE.md — Plan 35 status paragraph (delivered summary, suite counts)

- [ ] **Step 1: Re-audit the oracle-class counts** — run the Plan-32 counting
recipe (`tests/release/test_oracle_classes.py` documents it: collect with each
marker expression, compare to the verification-page table); update the page
numbers AND the gate's expected counts in one commit (they are cross-gated).

- [ ] **Step 2: Full local gates**

```bash
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/hochschild/ tests/resolutions_cs/test_products_cs.py tests/release -q
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/qpa -q -m qpa
.venv/bin/python -m build && .venv/bin/python -m twine check dist/*
```

Deep bucket in full runs on CI (the two legs); locally run at least
`tests/engine/test_tt_calculus.py tests/engine/test_gerstenhaber.py
tests/resolutions_cs/ -q` under nohup + Monitor.

- [ ] **Step 3: Commit docs, merge choreography**

```bash
git add docs/verification.md docs/internals/ README.md CLAUDE.md
git commit -m "docs: Plan 35 verification rows, internals section, README line"
# merge (after review):
git checkout main && git merge --ff-only plan-35-hh-products && git push origin main
# CI green on the merge commit -> Marco tags v0.1.0 (NOT the agent).
```

---

## Self-Review (run after writing, fixed inline)

1. **Spec coverage:** §1→Tasks 1–5; §2→Tasks 2–5 (+ the bracket amendment,
   Task 1); §3→Tasks 8–10; §4→Task 11; §5→Task 12; §6→Tasks 6–8, 12, 13;
   §7→Task 13; §8 non-goals honored (no representatives, no native bracket,
   no new curated example). Covered.
2. **Placeholders:** the plan directs verification of four call-shapes against
   the live tree (`from_structure_constants` name, `spec.run` entry names,
   `render_html` signature, gui.js module-kind order) — each names the exact
   file/line to read; no TBDs.
3. **Type consistency:** `HHProducts.tables: dict[(p,q)] -> ProductTable`,
   `constants[k][i][j] -> str` used identically in Tasks 1, 2, 3, 6, 8, 10,
   11; `ConnesB.matrices[n]` rows = HH_{n+1}-indexed, columns = HH_n-indexed —
   stated in Tasks 1, 4, and rendered as such in Task 10.
