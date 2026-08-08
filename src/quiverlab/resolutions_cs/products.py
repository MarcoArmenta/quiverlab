"""Plan 35 -- Domain-generic CS product tables: cup and cap on the CS HH basis
via the Plan-20/21 native collapses of the lifted diagonal. Any exact Domain,
any degree (no bar window). Plan 51 adds the Gerstenhaber BRACKET as a SEPARATE
sibling builder ``cs_bracket_tables`` (the Negron-Witherspoon / Volkov homotopy
lifting on the SAME diagonal); ``cs_product_tables`` keeps serving only cup/cap
and points a stray "bracket" call at the bracket builder."""
from quiverlab.errors import QuiverlabError
# _class_coords: the shared solve-in-[image|reps] lives in hochschild.products
# (stringify=True is exactly this module's former string-returning behavior).
from quiverlab.hochschild.products import (
    HHProducts, ProductTable, _class_coords, _pairs, _REFERENCES)


def _columns(M):
    """The columns of a list-of-rows matrix, each as a coordinate vector -- byte
    identical to ``resolutions_cs.homology._columns`` (the image basis the reps
    were reduced against, so reps and image share one orientation)."""
    return [[row[c] for row in M] for c in range(len(M[0]))] if M and M[0] else []


def cs_product_tables(A, kind, top, max_cells):
    if kind == "bracket":
        raise QuiverlabError(
            "cs_product_tables serves only cup/cap; the Gerstenhaber bracket has its "
            "own CS-native builder",
            hint="call cs_bracket_tables(A, top, max_cells) or "
                 "Algebra.gerstenhaber_brackets(top, engine='cs')")
    if kind not in ("cup", "cap"):
        raise QuiverlabError(f"unknown product kind {kind!r}")
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.homology import cs_hh_basis, _require_admissible
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    from quiverlab.resolutions_cs.cup import native_cup
    from quiverlab.resolutions_cs.cap import native_cap
    from quiverlab.hochschild import basis_reps as BR

    rs = reduction_system_of(A)
    _require_admissible(rs)
    dom = A.domain
    # native_cup/native_cap need S built one past the pair (p+q for cup, n for cap);
    # top+2 covers every pair (out degree <= top, diagonal reads one degree past).
    res = ChouhySolotarResolution(A, rs, max_degree=top + 2, max_cells=max_cells)
    coh = {n: cs_hh_basis(A, n, "coh", max_cells=max_cells) for n in range(top + 1)}
    hom = ({n: cs_hh_basis(A, n, "hom", max_cells=max_cells) for n in range(top + 1)}
           if kind == "cap" else {})

    def _image(out_n, side):
        """Columns of the differential whose image is the (co)boundary space at
        `out_n`: B^{out_n} = im delta^{out_n-1} (coh) / B_{out_n} = im b_{out_n+1}
        (hom). Both live in C^{out_n} and match the image used by cs_hh_basis."""
        if side == "coh":
            if out_n == 0:                       # B^0 = 0 (delta^{-1} = 0); no matrix(-1)
                return []
            M = res.matrix(out_n - 1, "coh")
        else:
            M = res.matrix(out_n + 1, "hom")
        return _columns(M)

    tables = {}
    for (p, q) in _pairs(kind, top):
        if kind == "cup":
            left, right, out_n, side = coh[p], coh[q], p + q, "coh"
            prod = lambda f, z, p=p, q=q: native_cup(res, f, p, z, q)
            out_reps = coh.get(out_n, [])
        else:
            left, right, out_n, side = coh[p], hom[q], q - p, "hom"
            prod = lambda f, z, p=p, q=q: native_cap(res, f, p, z, q)
            out_reps = hom.get(out_n, [])
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

    bc, cb, diffs = _capture_cs(A, res, coh, hom, kind, top)
    return HHProducts(kind=kind, top=top, tables=tables,
                      engine="Chouhy-Solotar native diagonal",
                      basis=f"cs/{A.domain.name}", window=None,
                      references=_REFERENCES[kind] + ["chouhy_solotar"],
                      basis_classes=bc, chain_basis=cb, differentials=diffs)


def cs_bracket_tables(A, top, max_cells):
    """Plan 51: the Gerstenhaber bracket table family HH^p (x) HH^q -> HH^{p+q-1}
    on the CS HH basis, over any exact Domain, PAST the bar window -- the sibling of
    cs_product_tables built from resolutions_cs.bracket.native_bracket (two
    Negron-Witherspoon / Volkov homotopy liftings on the shipped diagonal).  Same
    descent-in-[image|reps] class coordinates and explicit-reps capture; window=None
    (native, no bar bound)."""
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.homology import cs_hh_basis, _require_admissible
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    from quiverlab.resolutions_cs.bracket import native_bracket
    from quiverlab.resolutions_cs.homotopy_lifting import homotopy_lifting

    rs = reduction_system_of(A)
    _require_admissible(rs)
    dom = A.domain
    # native_bracket needs S/Delta up to p+q-1 <= top; top+2 covers the pair AND the
    # cocycle guard's one-past coboundary read (mirrors cs_product_tables' margin).
    res = ChouhySolotarResolution(A, rs, max_degree=top + 2, max_cells=max_cells)
    coh = {n: cs_hh_basis(A, n, "coh", max_cells=max_cells) for n in range(top + 1)}

    def _image(out_n):
        if out_n == 0:                       # B^0 = 0 (delta^{-1} = 0)
            return []
        return _columns(res.matrix(out_n - 1, "coh"))

    # DD1: one homotopy lifting per (degree, class index), reused across the whole
    # row/column of every table instead of rebuilt per (i, j) pair. Canonical solve ->
    # byte-identical to the uncached path (test_native_bracket_psi_cache_byte_identical).
    _lift = {}

    def lift(deg, idx):
        obj = _lift.get((deg, idx))
        if obj is None:
            obj = homotopy_lifting(res, coh[deg][idx], deg)
            _lift[(deg, idx)] = obj
        return obj

    tables = {}
    for (p, q) in _pairs("bracket", top):
        out_n = p + q - 1
        left, right = coh[p], coh[q]
        out_reps = coh.get(out_n, [])
        dl, dr, dout = len(left), len(right), len(out_reps)
        img = _image(out_n)
        consts = [[[None] * dr for _ in range(dl)] for _ in range(dout)]
        for i in range(dl):
            for j in range(dr):
                coords = _class_coords(
                    native_bracket(res, left[i], p, right[j], q,
                                   psi_f=lift(p, i), psi_g=lift(q, j)),
                    out_reps, img, dom)
                for k in range(dout):
                    consts[k][i][j] = coords[k]
        tables[(p, q)] = ProductTable(
            kind="bracket", degrees=(p, q), out_degree=out_n, dims=(dl, dr, dout),
            constants=tuple(tuple(tuple(row) for row in mat) for mat in consts))

    bc, cb, diffs = _capture_cs(A, res, coh, {}, "bracket", top)
    return HHProducts(
        kind="bracket", top=top, tables=tables,
        engine="Chouhy-Solotar native diagonal (homotopy lifting)",
        basis=f"cs/{A.domain.name}", window=None,
        references=_REFERENCES["bracket"] + [
            "bracket_liftings", "bracket_liftings_volkov", "chouhy_solotar", "oke_koszul"],
        basis_classes=bc, chain_basis=cb, differentials=diffs)


def _capture_cs(A, res, coh, hom, kind, top):
    """Serialize the CS explicit representatives (Plan 35): per degree the class
    list (coordinate vectors over res._basis(n, side)), the ordered CS chain
    enumeration, and the annihilating differential res.matrix(n, side). Cohomology
    side always; homology side too for cap."""
    from quiverlab.hochschild import basis_reps as BR
    dom = A.domain
    labels = BR.labels_of(res.ar.A)
    bc, cb, diffs = {}, {}, {}
    for n in range(top + 1):
        coh_elems = BR.cs_elements(res, n, "coh", labels)
        bc[("coh", n)] = BR.classes_from_columns(coh[n], coh_elems, n, "cochain", dom)
        cb[("coh", n)] = BR.enumeration_labels(coh_elems, "cochain")
        shape = (res.dim_C(n + 1, "coh"), res.dim_C(n, "coh"))
        note = ("ChouhySolotarResolution(A, reduction_system_of(A), "
                "max_degree>=%d).matrix(%d, 'coh')" % (n + 1, n))
        diffs[("coh", n)] = BR.serialize_differential(
            shape, (lambda nn=n: res.matrix(nn, "coh")), note, dom)
        if kind == "cap":
            hom_elems = BR.cs_elements(res, n, "hom", labels)
            bc[("hom", n)] = BR.classes_from_columns(hom[n], hom_elems, n, "chain", dom)
            cb[("hom", n)] = BR.enumeration_labels(hom_elems, "chain")
            hshape = (res.dim_C(n - 1, "hom") if n else 0, res.dim_C(n, "hom"))
            hnote = ("ChouhySolotarResolution(A, reduction_system_of(A), "
                     "max_degree>=%d).matrix(%d, 'hom')" % (n + 1, n))
            diffs[("hom", n)] = BR.serialize_differential(
                hshape, (lambda nn=n: res.matrix(nn, "hom") if nn else []), hnote, dom)
    return bc, cb, diffs
