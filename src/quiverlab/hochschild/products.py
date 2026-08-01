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
    from quiverlab.engine.scan3 import cochain_basis
    from quiverlab.errors import DepthLimitError
    prime = A.domain.p
    E = to_engine(A.unit_adapted())
    fn = {"cup": TT.cup_product_matrix, "cap": TT.cap_product_matrix,
          "bracket": TT.gerstenhaber_bracket_matrix}[kind]
    out_deg = {"cup": lambda p, q: p + q, "cap": lambda p, n: n - p,
               "bracket": lambda p, q: p + q - 1}[kind]
    tables = {}
    for (p, q) in _pairs(kind, top):
        cells = len(cochain_basis(E, p)) * len(cochain_basis(E, q))
        if cells > max_cells:
            raise DepthLimitError(
                f"{kind} product table ({p}, {q}): the two cochain bases pair "
                f"{cells} cells (> max_cells = {max_cells})",
                hint="raise max_cells or lower top")
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


# ---------------------------------------------------------------------------
# The one shared "coordinates in the class basis" solve. resolutions_cs/products
# imports this (stringify=True -> exact-string coords, its product-table form);
# the induced-B assembly below calls it with stringify=False (Domain elements).
# ---------------------------------------------------------------------------
def _class_coords(vec, reps, image_cols, dom, stringify=True):
    """Coordinates of `vec` in the class basis `reps`, modulo `image_cols`:
    solve [image | reps] x = vec and read off the reps segment. Loud when the
    vector is not a (co)cycle representative (descent failure).

    The reps segment is well-defined for any particular solution: the reps are
    linearly independent modulo span(image_cols) (that is how the class basis is
    picked), so two solutions differ only in the image block. `stringify=True`
    returns exact-string coordinates; False keeps Domain elements."""
    from quiverlab.errors import QuiverlabError
    from quiverlab.fields.linalg import solve
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
    tail = x[len(image_cols):]
    return [str(c) for c in tail] if stringify else list(tail)


# ---------------------------------------------------------------------------
# Induced Connes differential B : HH_n -> HH_{n+1}
# ---------------------------------------------------------------------------
def _as_matrix(triple):
    """Normalize a differential return `(rows, ncols, nrows)` to
    `(row-major matrix, ncols, nrows)`. boundary_matrix already yields the
    matrix first; this pins the unpacking convention in one place."""
    M, ncols, nrows = triple
    return M, ncols, nrows


def _cols_to_matrix(cols):
    """A list of column vectors -> the row-major matrix whose columns they are
    (rows indexed by the ambient coordinate), ready for fields.linalg.rank."""
    if not cols:
        return []
    L = len(cols[0])
    return [[cols[c][r] for c in range(len(cols))] for r in range(L)]


def _int_rank_mod_p(rows, ncols, p):
    """Rank over F_p of the small integer matrix `rows` (each row length
    `ncols`), via engine.coxeter.rref_mod_p pivot count."""
    if not rows or ncols == 0:
        return 0
    import numpy as np
    from quiverlab.engine.coxeter import rref_mod_p
    _, piv = rref_mod_p(np.array(rows, dtype=np.int64) % p, p)
    return len(piv)


def _generic_homology_quotient(A, n, max_cells):
    """(reps, image_cols) for HH_n over A.domain from the bar boundary matrices
    (the Domain-generic sibling of engine.tt_calculus.homology_classes). reps
    are cycles in C_n kept independent modulo image_cols = the columns of
    b_{n+1} (the boundaries), so reps and image share one orientation."""
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
            reps.append(v)
            base.append(v)
            r0 = rr
    return reps, image


def connes_b_tables(A, top, max_cells=4_000_000):
    """Induced Connes differentials B: HH_n -> HH_{n+1}, 0 <= n < top.

    Returns a ConnesB whose matrices[n] is the hh_{n+1} x hh_n matrix of the
    induced B in the class bases (rows indexed by HH_{n+1}); B^2 = 0 holds at
    the induced level because bB + Bb = 0 makes B descend to homology and B^2
    vanishes on chains. GF(p): the fast (b,B) engine; any other exact Domain:
    the generic bar (b,B) mixed complex (max_cells guards the blow-up)."""
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
            # rows[i] = coords of B(e_i); store as matrix hh_{n+1} x hh_n
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
                cols.append(_class_coords(w, reps_n1, img_n1, dom, stringify=False))
            matrices[n] = [[str(cols[i][k]) for i in range(len(reps_n))]
                           for k in range(len(reps_n1))]
            ranks[n] = _rank([[dom.coerce(cols[i][k]) for i in range(len(reps_n))]
                              for k in range(len(reps_n1))], dom) if reps_n else 0
        hh = [len(quots[n][0]) for n in range(top + 1)]
        engine = f"generic (b,B) mixed complex / {dom.name}"
    return ConnesB(top=top, hh_dims=hh, matrices=matrices, ranks=ranks,
                   engine=engine, references=_REFERENCES["connes_b"])
