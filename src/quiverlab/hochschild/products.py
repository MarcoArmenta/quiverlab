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


def _by_side(mapping):
    """Serialize a ``{(side, n): value}`` dict to the JSON-safe nested
    ``{side: {str(n): value}}`` shape (Plan 35 explicit-reps block fields)."""
    out = {}
    for (side, n), value in sorted(mapping.items()):
        out.setdefault(side, {})[str(n)] = value
    return out


class HHProducts:
    """A family of product tables up to `top`. kind in {"cup","cap","bracket"}."""

    def __init__(self, kind, top, tables, engine, basis, window, references,
                 basis_classes=None, chain_basis=None, differentials=None):
        self.kind = kind
        self.top = top
        self.tables = dict(tables)     # {(p, q): ProductTable}
        self.engine = engine
        self.basis = basis
        self.window = window           # int for bracket (served window), else None
        self.references = list(references)
        # Plan 35 explicit representatives: {(side, n): [class dict, ...]} the ACTUAL
        # (co)cycles that produced the constants; the ordered enumeration they index
        # into; and the annihilating differential (self-certification). None when a
        # legacy caller omits them (blocks() then omits the fields; renderers fall back).
        self.basis_classes = dict(basis_classes) if basis_classes else None
        self.chain_basis = dict(chain_basis) if chain_basis else None
        self.differentials = dict(differentials) if differentials else None

    def blocks(self):
        out = {"kind": self.kind, "top": self.top, "engine": self.engine,
               "basis": self.basis,
               "tables": [self.tables[k].as_dict()
                          for k in sorted(self.tables)],
               "references": list(self.references)}
        if self.window is not None:
            out["window"] = self.window
        if self.basis_classes is not None:
            out["basis_classes"] = _by_side(self.basis_classes)
        if self.chain_basis is not None:
            out["chain_basis"] = _by_side(self.chain_basis)
        if self.differentials is not None:
            out["differentials"] = _by_side(self.differentials)
        return out

    def __repr__(self):
        return (f"<HHProducts {self.kind} top={self.top} "
                f"tables={len(self.tables)} basis={self.basis!r}>")


class ConnesB:
    """Induced Connes differentials B: HH_n -> HH_{n+1} for 0 <= n < top."""

    def __init__(self, top, hh_dims, matrices, ranks, engine, references,
                 basis_classes=None, chain_basis=None, differentials=None):
        self.top = top
        self.hh_dims = list(hh_dims)   # dim HH_0..HH_top
        self.matrices = dict(matrices) # {n: rows of str, shape hh_dims[n+1] x hh_dims[n]}
        self.ranks = dict(ranks)       # {n: int}
        self.engine = engine
        self.references = list(references)
        # Plan 35 explicit representatives: the homology cycle bases z^n_j (0..top),
        # their ordered chain enumeration, and the boundary b_n that annihilates them.
        self.basis_classes = dict(basis_classes) if basis_classes else None
        self.chain_basis = dict(chain_basis) if chain_basis else None
        self.differentials = dict(differentials) if differentials else None

    def blocks(self):
        out = {"kind": "connes_b", "top": self.top,
               "hh_dims": list(self.hh_dims),
               "matrices": {str(n): self.matrices[n] for n in sorted(self.matrices)},
               "ranks": {str(n): self.ranks[n] for n in sorted(self.ranks)},
               "engine": self.engine, "references": list(self.references)}
        if self.basis_classes is not None:
            out["basis_classes"] = _by_side(self.basis_classes)
        if self.chain_basis is not None:
            out["chain_basis"] = _by_side(self.chain_basis)
        if self.differentials is not None:
            out["differentials"] = _by_side(self.differentials)
        return out

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
    from quiverlab.hochschild import basis_reps as BR
    prime = A.domain.p
    AU = A.unit_adapted()
    E = to_engine(AU)
    labels = BR.labels_of(AU)
    out_deg = {"cup": lambda p, q: p + q, "cap": lambda p, n: n - p,
               "bracket": lambda p, q: p + q - 1}[kind]

    # Class caches: cohomology/homology reps computed ONCE and shared -- the SAME
    # objects both produce the constants (passed into the tt matrix builders) and
    # are captured as the explicit representatives (Plan 35).
    _coh, _hom = {}, {}

    def coh(n):
        if n not in _coh:
            _coh[n] = TT.cohomology_classes(E, n, prime)
        return _coh[n]

    def hom(n):
        if n not in _hom:
            _hom[n] = TT.homology_classes(E, n, prime)
        return _hom[n]

    tables = {}
    for (p, q) in _pairs(kind, top):
        cells = len(cochain_basis(E, p)) * len(cochain_basis(E, q))
        if cells > max_cells:
            raise DepthLimitError(
                f"{kind} product table ({p}, {q}): the two cochain bases pair "
                f"{cells} cells (> max_cells = {max_cells})",
                hint="raise max_cells or lower top")
        if kind == "cup":
            C, dl, dr, dout = TT.cup_product_matrix(E, p, q, prime,
                                                    coh(p), coh(q), coh(p + q))
        elif kind == "bracket":
            C, dl, dr, dout = TT.gerstenhaber_bracket_matrix(
                E, p, q, prime, coh(p), coh(q), coh(p + q - 1))
        else:  # cap: (p, n) -> HH^p (x) HH_n -> HH_{n-p}
            C, dl, dr, dout = TT.cap_product_matrix(E, p, q, prime,
                                                    coh(p), hom(q), hom(q - p))
        tables[(p, q)] = ProductTable(
            kind=kind, degrees=(p, q), out_degree=out_deg(p, q),
            dims=(dl, dr, dout),
            constants=tuple(tuple(tuple(str(int(C[k, i, j])) for j in range(dr))
                                  for i in range(dl)) for k in range(dout)))

    bc, cb, diffs = _capture_gfp(E, labels, kind, top, coh, hom)
    window = top if kind == "bracket" else None
    return HHProducts(kind=kind, top=top, tables=tables,
                      engine="hanlab engine (F_p fast rank)",
                      basis=f"bar/GF({prime})", window=window,
                      references=_REFERENCES[kind],
                      basis_classes=bc, chain_basis=cb, differentials=diffs)


def _gfp_coh_diff(E, n):
    """(shape, build, note) for the coboundary delta^n: C^n -> C^{n+1} (engine),
    the differential annihilating a degree-n cohomology class."""
    m = E.m
    shape = (m * (m - 1) ** (n + 1), m * (m - 1) ** n)

    def build():
        from quiverlab.engine.scan3 import cochain_basis, coboundary_matrix
        idx = {g: i for i, g in enumerate(cochain_basis(E, n + 1))}
        return coboundary_matrix(E, n, cochain_basis(E, n), idx)

    note = ("quiverlab.engine.scan3.coboundary_matrix(E, %d, cochain_basis(E, %d), "
            "{g: i for i, g in enumerate(cochain_basis(E, %d))}), "
            "E = to_engine(A.unit_adapted())" % (n, n, n + 1))
    return shape, build, note


def _gfp_hom_diff(E, n):
    """(shape, build, note) for the boundary b_n: C_n -> C_{n-1} (engine), which
    annihilates a degree-n homology class. b_0 = 0 (no map out of C_0)."""
    m = E.m
    if n == 0:
        return (0, m), (lambda: []), "b_0 = 0 (every 0-chain is a cycle)"
    shape = (m * (m - 1) ** (n - 1), m * (m - 1) ** n)

    def build():
        from quiverlab.engine.hh_engine import cn_basis, differential_matrix
        idx = {g: i for i, g in enumerate(cn_basis(E, n - 1))}
        return differential_matrix(E, n, cn_basis(E, n), idx)

    note = ("quiverlab.engine.hh_engine.differential_matrix(E, %d, cn_basis(E, %d), "
            "{g: i for i, g in enumerate(cn_basis(E, %d))}), "
            "E = to_engine(A.unit_adapted())" % (n, n, n - 1))
    return shape, build, note


def _capture_gfp(E, labels, kind, top, coh, hom):
    """Serialize the GF(p) explicit representatives: per degree the class list, the
    ordered enumeration labels, and the annihilating differential. Cohomology side
    for every kind; homology side too for cap."""
    from quiverlab.hochschild import basis_reps as BR
    bc, cb, diffs = {}, {}, {}
    for n in range(top + 1):
        coh_elems = BR.engine_coh_elements(E, n, labels)
        cols = [coh(n).reps[:, i] for i in range(coh(n).reps.shape[1])]
        bc[("coh", n)] = BR.classes_from_columns(cols, coh_elems, n, "cochain", None)
        cb[("coh", n)] = BR.enumeration_labels(coh_elems, "cochain")
        shape, build, note = _gfp_coh_diff(E, n)
        diffs[("coh", n)] = BR.serialize_differential(shape, build, note, None)
        if kind == "cap":
            hom_elems = BR.engine_hom_elements(E, n, labels)
            hcols = [hom(n).reps[:, i] for i in range(hom(n).reps.shape[1])]
            bc[("hom", n)] = BR.classes_from_columns(hcols, hom_elems, n, "chain", None)
            cb[("hom", n)] = BR.enumeration_labels(hom_elems, "chain")
            shape, build, note = _gfp_hom_diff(E, n)
            diffs[("hom", n)] = BR.serialize_differential(shape, build, note, None)
    return bc, cb, diffs


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
    from quiverlab.hochschild import basis_reps as BR
    dom = A.domain
    AU = A.unit_adapted()
    labels = BR.labels_of(AU)
    if isinstance(dom, PrimeField):
        import numpy as np
        from quiverlab.errors import DepthLimitError
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.cyclic import connes_B_matrix
        from quiverlab.engine.hh_engine import cn_basis
        from quiverlab.engine.tt_calculus import homology_classes
        p = dom.p
        E = to_engine(AU)
        # Guard the bar (b,B) blow-up BEFORE building any matrix. homology_classes(n)
        # materializes the DENSE boundary matrices b_n (cn[n-1] x cn[n]) and b_{n+1}
        # (cn[n] x cn[n+1]); the bar chain basis grows exponentially with degree, so at
        # high top these are gigabytes (e.g. QuantumCI connes_b:0..7 needs b_8, an
        # 8748 x 26244 int64 matrix ~= 1.8 GB, whose rref/nullspace copies peaked ~6 GB
        # and SIGKILLed the memory-capped worker). This mirrors the cup/cap/bracket
        # cochain-pair guard: over GF(p) connes has NO Chouhy-Solotar route, so -- like
        # an explicit bar engine -- it refuses LOUDLY rather than silently OOMing. The
        # max_cells parameter was already honoured on the generic-Domain branch below;
        # the GF(p) branch simply dropped it (the bug). Dims/other invariants stand.
        cdims = [len(cn_basis(E, n)) for n in range(top + 2)]
        for n in range(top + 1):
            cells = max(cdims[n] * cdims[n + 1],
                        cdims[n - 1] * cdims[n] if n else 0)
            if cells > max_cells:
                raise DepthLimitError(
                    f"connes_b: the bar (b,B) boundary at degree {n} pairs {cells} "
                    f"cells (> max_cells = {max_cells})",
                    hint="raise max_cells or lower top")
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
        bc, cb, diffs = {}, {}, {}
        for n in range(top + 1):
            elems = BR.engine_hom_elements(E, n, labels)
            cols = [H[n].reps[:, i] for i in range(H[n].reps.shape[1])]
            bc[("hom", n)] = BR.classes_from_columns(cols, elems, n, "chain", None)
            cb[("hom", n)] = BR.enumeration_labels(elems, "chain")
            shape, build, note = _gfp_hom_diff(E, n)
            diffs[("hom", n)] = BR.serialize_differential(shape, build, note, None)
    else:
        from quiverlab.fields.linalg import rank as _rank
        from quiverlab.hochschild.cyclic import connes_B_matrix
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
        bc, cb, diffs = {}, {}, {}
        m = AU.dim
        for n in range(top + 1):
            elems = BR.bar_chain_elements(m, n, labels)
            bc[("hom", n)] = BR.classes_from_columns(quots[n][0], elems, n, "chain", dom)
            cb[("hom", n)] = BR.enumeration_labels(elems, "chain")
            shape = ((m * (m - 1) ** (n - 1)) if n else 0, m * (m - 1) ** n)
            note = ("quiverlab.hochschild.bar.boundary_matrix(A.unit_adapted(), %d, "
                    "max_cells)[0]" % n)
            build = (lambda nn=n: _generic_boundary_rows(AU, nn, max_cells))
            diffs[("hom", n)] = BR.serialize_differential(shape, build, note, dom)
    return ConnesB(top=top, hh_dims=hh, matrices=matrices, ranks=ranks,
                   engine=engine, references=_REFERENCES["connes_b"],
                   basis_classes=bc, chain_basis=cb, differentials=diffs)


def _generic_boundary_rows(AU, n, max_cells):
    """Row-major boundary matrix b_n over the Domain (empty for n=0: b_0 = 0)."""
    if n == 0:
        return []
    from quiverlab.hochschild.bar import boundary_matrix
    M, _, _ = boundary_matrix(AU, n, max_cells)
    return M
