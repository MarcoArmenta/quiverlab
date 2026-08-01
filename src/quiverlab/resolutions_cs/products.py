"""Plan 35 -- Domain-generic CS product tables: cup and cap on the CS HH basis
via the Plan-20/21 native collapses of the lifted diagonal. Any exact Domain,
any degree (no bar window). The bracket is NOT served here: its only route is
the GF(p) tt facade (see the Plan-35 spec amendment)."""
from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import solve
from quiverlab.hochschild.products import HHProducts, ProductTable, _pairs, _REFERENCES


def _columns(M):
    """The columns of a list-of-rows matrix, each as a coordinate vector -- byte
    identical to ``resolutions_cs.homology._columns`` (the image basis the reps
    were reduced against, so reps and image share one orientation)."""
    return [[row[c] for row in M] for c in range(len(M[0]))] if M and M[0] else []


def _class_coords(vec, reps, image_cols, dom):
    """Coordinates of `vec` in the class basis `reps`, modulo `image_cols`:
    solve [image | reps] x = vec and read off the reps segment. Loud when the
    vector is not a (co)cycle representative (descent failure).

    The reps segment is well-defined for any particular solution: the reps are
    linearly independent modulo span(image_cols) (that is how ``cs_hh_basis``
    picks them), so two solutions differ only in the image block."""
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
    return HHProducts(kind=kind, top=top, tables=tables,
                      engine="Chouhy-Solotar native diagonal",
                      basis=f"cs/{A.domain.name}", window=None,
                      references=_REFERENCES[kind] + ["chouhy_solotar"])
