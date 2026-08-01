"""Plan 35 wave 3b -- explicit cyclic-homology (HC) representatives + self-cert data.

The HC sibling of ``hochschild.basis_reps`` (HH products) and ``modules.complex_reps``
(Ext / Tor): here we LABEL and SERIALIZE the explicit HC class representatives that the
two ``cyclic_homology_dims`` engines expose (GF(p) fast rank / generic (b, B) mixed
complex), so a reader can read off each class AND verify it is a genuine cycle of the
total complex.

Cyclic homology is the homology of the total complex of Connes' (b, B) bicomplex:

    Tot_n = C_n (+) C_{n-2} (+) C_{n-4} (+) ... ,   D = b + B : Tot_n -> Tot_{n-1}

(``C_k = A (x) Abar^{(x)k}`` the normalized bar chains). A class is a coset of
ker D_n modulo im D_{n+1}. The ordered basis of ``Tot_n`` is the CONCATENATION of the
per-column bar chain bases, so the coordinate vector of a class has an explicit column
structure -- which slice lives in which ``C_k`` -- carried alongside as
``column_structure`` and stamped into every enumeration label as ``col C_k: <chain>``.

Two coherent views of every class (mirroring the siblings):

  * ``terms``  -- the labelled term-sum ``Sigma c * (col C_d: v (x) w1 (x) ... )``;
  * ``vector`` -- SPARSE coordinates ``[[index, coeff_str], ...]`` over the ordered
    ``Tot_n`` enumeration the ``chain_basis`` lists.

The reps are the engine's own quotient columns (never recomputed against a different
basis); the captured rep count is asserted equal to the HHTable HC dims, so the reps
and the numbers the report prints are the same computation.

Float-free: every coefficient / matrix entry is an int or an exact Domain-element
string (the AST gate scans this file).
"""
from quiverlab.errors import QuiverlabError
from quiverlab.hochschild import basis_reps as BR
from quiverlab.hochschild.basis_reps import (
    MATRIX_CELL_CAP,
    _coeff_str,
    _is_zero,
    serialize_differential,
)


def _tot_degrees(n):
    """Chain degrees in Tot_n = C_n (+) C_{n-2} (+) ... (descending)."""
    return list(range(n, -1, -2))


# --------------------------------------------------------------------------- #
# The ordered Tot_n basis: the concatenated per-column bar chain bases, each element
# stamped with its column degree so the reader knows which C_k slice it lives in.
# --------------------------------------------------------------------------- #
def _tot_elements(n, chain_elements_of):
    """Ordered ``(col_degree, word_labels, value_label)`` enumeration of Tot_n."""
    out = []
    for d in _tot_degrees(n):
        for (w, v) in chain_elements_of(d):
            out.append((d, tuple(w), v))
    return out


def column_structure(n, col_dims):
    """The explicit column layout of Tot_n: ``{'total': dim, 'columns': [{'degree',
    'offset', 'dim'}, ...]}`` -- which coordinate slice lives in which ``C_k``."""
    cols, off = [], 0
    for d in _tot_degrees(n):
        dd = int(col_dims[d])
        cols.append({"degree": d, "offset": off, "dim": dd})
        off += dd
    return {"total": off, "columns": cols}


def enumeration_labels(elements):
    """The ordered Tot_n enumeration as display labels ``col C_d: v (x) w1 (x) ...``
    (what the coordinate vectors index into). Capped: an over-long enumeration ships an
    elided marker with the length and a rebuild pointer -- never silently truncated."""
    if len(elements) > MATRIX_CELL_CAP:
        return {"elided": True, "length": len(elements),
                "note": "enumeration exceeds the 250000-cell display cap; "
                        "rebuild from the ordered total-complex basis"}
    return ["col C_%d: %s" % (d, BR.element_label(list(w), v, "chain"))
            for (d, w, v) in elements]


def serialize_class(column, elements, n, dom):
    """One HC class = one coordinate column over the ordered Tot_n enumeration; the
    term-sum carries each term's column degree (``[coeff, col_degree, word, value]``)."""
    terms, vector = [], []
    for idx, coeff in enumerate(column):
        if _is_zero(coeff, dom):
            continue
        cs = _coeff_str(coeff, dom)
        d, w, v = elements[idx]
        terms.append([cs, d, list(w), v])
        vector.append([idx, cs])
    return {"kind": "cyclic", "degree": n, "terms": terms, "vector": vector}


# --------------------------------------------------------------------------- #
# The total differential D_n = b + B : Tot_n -> Tot_{n-1}, serialized (capped 250k,
# elided + rebuild-note over the cap; the note names the cyclic engine to rebuild).
# --------------------------------------------------------------------------- #
def _serialize_total_diff(D, n, dom, ncols, nrows, note_fmt):
    if nrows == 0:                               # D_0: every 0-chain is a cycle
        return {"shape": [0, int(ncols)], "rows": [],
                "note": "D_0 = 0 on Tot_0 = C_0; every 0-chain is a cycle "
                        "(HC_0 = A/[A,A])"}
    note = note_fmt % (n, n)
    return serialize_differential((nrows, ncols), (lambda: (D or [])), note, dom)


# --------------------------------------------------------------------------- #
# The shared payload builder. `chain_elements_of(d)` returns the ordered bar chain
# basis of C_d as (word_labels, value_label) pairs (engine int-index labels or the
# generic bar labels); `dom` is None on the GF(p) route (int coeffs), a Domain
# otherwise. `raw` is the engine's exposed quotient (reps columns + total matrices).
# --------------------------------------------------------------------------- #
def _build_payload(dims, top, chain_elements_of, col_dims, raw, dom, note_fmt):
    bc, cb, diffs, colstruct, out_dims = {}, {}, {}, {}, []
    for n in range(top + 1):
        elements = _tot_elements(n, chain_elements_of)
        cs_n = column_structure(n, col_dims)
        ncols = cs_n["total"]
        nrows = column_structure(n - 1, col_dims)["total"] if n >= 1 else 0
        cols = raw["reps"].get(n) or []
        out_dims.append(len(cols))
        bc[str(n)] = [serialize_class(col, elements, n, dom) for col in cols]
        cb[str(n)] = enumeration_labels(elements)
        colstruct[str(n)] = cs_n
        diffs[str(n)] = _serialize_total_diff(raw["totmats"].get(n), n, dom,
                                              ncols, nrows, note_fmt)
    _cross_check(out_dims, dims)
    return {"basis_classes": bc, "chain_basis": cb, "differentials": diffs,
            "column_structure": colstruct}


def _cross_check(captured, dims):
    """The rep pass must return EXACTLY the HC table dims (basis independence). A drift
    is a capture bug -- refuse loudly rather than ship reps that disagree with the number
    the report prints."""
    if list(captured) != list(dims):
        raise QuiverlabError(
            "cyclic reps: captured dims %s disagree with HC dims %s -- the explicit "
            "representatives drifted from the computed dimensions" % (captured, dims))


_GFP_NOTE = ("D_%d = b + B on Tot_%d; rebuild via "
             "quiverlab.engine.cyclic.cyclic_homology_dims(to_engine(A.unit_adapted()), "
             "top, primes=(p,), with_reps=True)")
_GEN_NOTE = ("D_%d = b + B on Tot_%d; rebuild via "
             "quiverlab.hochschild.cyclic.cyclic_homology_dims(A, top, max_cells, "
             "with_reps=True)")


def gfp_payload(E, AU, top, raw, dims):
    """Label + serialize the GF(p) engine's exposed HC quotient. ``E`` is the engine
    algebra, ``AU`` the unit-adapted source (its basis labels name the values), ``raw``
    the ``engine.cyclic.cyclic_homology_dims(..., with_reps=True)`` payload, ``dims`` the
    HC dimension list (for the drift cross-check)."""
    labels = BR.labels_of(AU)
    return _build_payload(
        dims, top,
        lambda d: BR.engine_hom_elements(E, d, labels),
        raw["col_dims"], raw, None, _GFP_NOTE)


def generic_payload(A, top, raw, dims):
    """Label + serialize the generic (b, B) engine's exposed HC quotient over
    ``A.domain``; the ordered bar chain basis names the value/word labels."""
    AU = A.unit_adapted()
    labels = BR.labels_of(AU)
    m = AU.dim
    return _build_payload(
        dims, top,
        lambda d: BR.bar_chain_elements(m, d, labels),
        raw["col_dims"], raw, A.domain, _GEN_NOTE)
