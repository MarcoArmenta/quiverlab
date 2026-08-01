"""Plan 35 wave 3d -- explicit representatives of the PLAIN ``hh_cohomology`` /
``hh_homology`` dims blocks.

The HH-DIMS sibling of ``hochschild.basis_reps`` (HH products), ``modules.complex_reps``
(Ext / Tor) and ``hochschild.cyclic_reps`` (HC): here we LABEL and SERIALIZE the
explicit HH^n / HH_n class representatives that the SAME dims path exposes, so a reader
can read off each class -- the central elements of ``HH^0``, the derivations of ``HH^1``,
the deformation cochains of ``HH^2``, the commutator residues of ``HH_0`` -- AND verify
each is a genuine (co)cycle of the bar / Chouhy-Solotar complex.

Two routes, matched to the engine that produced the dims (the SAME-arrays /
SAME-resolution discipline, never a second independent computation that could pick a
different basis and silently mismatch the numbers the report prints):

  * the GF(p) bar route (``engine.tt_calculus.cohomology_classes`` /
    ``homology_classes``) when the dims came from the hanlab F_p fast engine (or, as a
    fallback, any prime field);
  * the Chouhy-Solotar route (``resolutions_cs.homology.cs_hh_basis``) when the dims
    came from CS (or, as a fallback, any admissible quiver-presented algebra over a
    Domain that has no F_p route).

Presentation / domain combos with NEITHER route (a structure-constants algebra off a
prime field, an inadmissible reduction system, or a degree beyond the tractable bar
cell cap) keep the dims ONLY -- the reps fields are simply absent (the renderers
tolerate that). A CAPTURED degree whose rep count disagrees with the table dimension is
a genuine capture bug and is refused loudly (never shipped).

Each captured class carries the three coherent views ``basis_reps.serialize_class``
builds (a labelled term-sum, a sparse coordinate vector, the degree) over the ordered
(co)chain enumeration the ``chain_basis`` lists. Alongside every degree we ship the
dimension of the (co)boundary subspace (``inner_dims[str(n)]``) -- for HH^1 this is the
dimension of the inner derivations, ``rank(delta^0)``.

Float-free: every coefficient / matrix entry is an int or an exact Domain-element
string (the AST gate scans this file).
"""
from quiverlab.errors import QuiverlabError
from quiverlab.hochschild import basis_reps as BR


def _is_prime_field(dom):
    from quiverlab.fields.primefield import PrimeField
    return isinstance(dom, PrimeField)


def _chain_kind(side):
    return "cochain" if side == "coh" else "chain"


def _check_degree(count, dims, n):
    """A captured degree MUST have exactly the table's dimension (basis independence).
    A drift is a capture bug -- refuse loudly rather than ship reps that disagree with
    the number the report prints."""
    if n < len(dims or []) and count != int(dims[n]):
        raise QuiverlabError(
            "hh reps: captured %d classes in degree %d but the table dim is %d -- the "
            "explicit representatives drifted from the computed dimensions"
            % (count, n, int(dims[n])))


# --------------------------------------------------------------------------- #
# Public entry: the additive reps block, or None (dims-only, honestly).
# --------------------------------------------------------------------------- #
def hh_reps_blocks(A, kind, top, dims, engine, max_cells=4_000_000):
    """The additive per-degree explicit-representatives block for a plain
    ``hh_cohomology`` / ``hh_homology`` dims table -- a dict with single-side
    ``{str(degree): ...}`` ``basis_classes`` / ``chain_basis`` / ``differentials`` plus
    ``inner_dims`` -- or ``None`` when no representative route applies (the block then
    stays dims-only)."""
    side = "coh" if kind == "hh_cohomology" else "hom"
    eng = str(engine or "").lower()
    from quiverlab.errors import AdmissibilityError, DepthLimitError
    try:
        if "chouhy" in eng or "solotar" in eng:
            return _capture_cs(A, side, top, dims, max_cells)
        if "hanlab" in eng or "f_p" in eng or "fast" in eng:
            return _capture_gfp(A, side, top, dims, max_cells)
        # Bardzell / minimal / generic bar dims: prefer a matching representative route.
        if _is_prime_field(A.domain):
            return _capture_gfp(A, side, top, dims, max_cells)
        return _capture_cs(A, side, top, dims, max_cells)
    except (AdmissibilityError, DepthLimitError, NotImplementedError):
        return None                                  # no route -- dims-only, honestly


def _payload(bc, cb, diffs, inner):
    if not bc:
        return None
    return {"basis_classes": bc, "chain_basis": cb, "differentials": diffs,
            "inner_dims": inner}


# --------------------------------------------------------------------------- #
# GF(p) bar route (engine.tt_calculus): the SAME arrays the F_p fast dims path uses.
# --------------------------------------------------------------------------- #
def _capture_gfp(A, side, top, dims, max_cells):
    if not _is_prime_field(A.domain):
        return None
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine import tt_calculus as TT
    from quiverlab.engine.scan3 import cochain_basis
    from quiverlab.engine.hh_engine import cn_basis
    from quiverlab.hochschild.products import _gfp_coh_diff, _gfp_hom_diff
    prime = A.domain.p
    AU = A.unit_adapted()
    E = to_engine(AU)
    labels = BR.labels_of(AU)
    ckind = _chain_kind(side)
    bc, cb, diffs, inner = {}, {}, {}, {}
    for n in range(top + 1):
        # Capturing degree n builds the annihilating differential; refuse an intractable
        # bar blow-up (the fast engine reaches these degrees WITHOUT the full bar
        # complex, so a partial reps capture is honest -- the dims still stand).
        if side == "coh":
            cells = len(cochain_basis(E, n + 1)) * len(cochain_basis(E, n))
        else:
            cells = len(cn_basis(E, n)) * (len(cn_basis(E, n - 1)) if n else 1)
        if cells > max_cells:
            break
        if side == "coh":
            Q = TT.cohomology_classes(E, n, prime)
            elems = BR.engine_coh_elements(E, n, labels)
            shape, build, note = _gfp_coh_diff(E, n)
        else:
            Q = TT.homology_classes(E, n, prime)
            elems = BR.engine_hom_elements(E, n, labels)
            shape, build, note = _gfp_hom_diff(E, n)
        cols = [Q.reps[:, i] for i in range(Q.reps.shape[1])]
        _check_degree(len(cols), dims, n)
        bc[str(n)] = BR.classes_from_columns(cols, elems, n, ckind, None)
        cb[str(n)] = BR.enumeration_labels(elems, ckind)
        diffs[str(n)] = BR.serialize_differential(shape, build, note, None)
        inner[str(n)] = int(Q.nim)                   # dim of the (co)boundary subspace
    return _payload(bc, cb, diffs, inner)


# --------------------------------------------------------------------------- #
# Chouhy-Solotar route (resolutions_cs): representative (co)cycles over the CS basis,
# the SAME resolution the CS dims path uses (canonical since Plan 17 -> reproducible).
# --------------------------------------------------------------------------- #
def _capture_cs(A, side, top, dims, max_cells):
    from quiverlab.fields.linalg import rank
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.homology import _require_admissible, cs_hh_basis
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    if A.quiver is None or A.relations is None:
        return None                                   # presentation-less -> no CS route
    rs = reduction_system_of(A)
    _require_admissible(rs)                           # AdmissibilityError if not
    dom = A.domain
    AU = A.unit_adapted()
    labels = BR.labels_of(AU)
    res = ChouhySolotarResolution(A, rs, max_degree=top + 1, max_cells=max_cells)
    key = side                                        # "coh" / "hom"
    ckind = _chain_kind(side)
    bc, cb, diffs, inner = {}, {}, {}, {}
    for n in range(top + 1):
        elems = BR.cs_elements(res, n, key, labels)
        reps = cs_hh_basis(A, n, key, max_cells=max_cells)
        _check_degree(len(reps), dims, n)
        bc[str(n)] = BR.classes_from_columns(reps, elems, n, ckind, dom)
        cb[str(n)] = BR.enumeration_labels(elems, ckind)
        diffs[str(n)] = _cs_differential(res, n, side, dom)
        inner[str(n)] = _cs_inner_dim(res, n, side, dom, rank)
    return _payload(bc, cb, diffs, inner)


def _cs_differential(res, n, side, dom):
    """The annihilating differential of a degree-n CS class, serialized (capped 250k,
    elided + rebuild-note over the cap): delta^n : C^n -> C^{n+1} (coh) / b_n : C_n ->
    C_{n-1} (hom; b_0 = 0)."""
    if side == "coh":
        rows, cols = res.dim_C(n + 1, "coh"), res.dim_C(n, "coh")
        M = res.matrix(n, "coh")
        note = ("res.matrix(%d, 'coh'), res = ChouhySolotarResolution(A, "
                "reduction_system_of(A))" % n)
    else:
        if n == 0:
            return {"shape": [0, int(res.dim_C(0, "hom"))], "rows": [],
                    "note": "b_0 = 0 (every 0-chain is a cycle; HH_0 = A/[A,A])"}
        rows, cols = res.dim_C(n - 1, "hom"), res.dim_C(n, "hom")
        M = res.matrix(n, "hom")
        note = ("res.matrix(%d, 'hom'), res = ChouhySolotarResolution(A, "
                "reduction_system_of(A))" % n)
    return BR.serialize_differential((rows, cols), (lambda: M), note, dom)


def _cs_inner_dim(res, n, side, dom, rank):
    """The dimension of the (co)boundary subspace landing in degree n: rank(delta^{n-1})
    (coh; the inner derivations for n = 1) / rank(b_{n+1}) (hom)."""
    if side == "coh":
        return int(rank(res.matrix(n - 1, "coh"), dom)) if n >= 1 else 0
    return int(rank(res.matrix(n + 1, "hom"), dom))
