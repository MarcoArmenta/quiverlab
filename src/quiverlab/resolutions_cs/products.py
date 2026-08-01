"""Plan 35 -- Domain-generic CS product tables: cup and cap on the CS HH basis
via the Plan-20/21 native collapses of the lifted diagonal. Any exact Domain,
any degree (no bar window). The bracket is NOT served here: its only route is
the GF(p) tt facade (see the Plan-35 spec amendment)."""
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
