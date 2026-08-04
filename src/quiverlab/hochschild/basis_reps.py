"""Plan 35 explicit-representatives -- capture the ACTUAL (co)cycle representatives
that produced the HH product structure constants, label them with the algebra's
basis, and serialize them so a reader can verify each class is a genuine (co)cycle.

Each captured class carries three coherent views of the same representative:

  * ``terms``  -- the labeled term-sum  ``Sigma c * [w1|...|wp -> v]`` (cohomology,
    bar-bracket form) / ``Sigma c * (v (x) w1 (x) ... (x) wn)`` (homology);
  * ``vector`` -- SPARSE coordinates ``[[index, coeff_str], ...]`` over the ordered
    (co)chain basis of the degree (the SAME enumeration ``chain_basis`` lists);
  * ``kind`` / ``degree`` -- so the reader knows which differential annihilates it.

Reps are CAPTURED at table-build time from the same class objects that produce the
constants, never recomputed independently (recomputation could pick a different
basis and silently mismatch the constants). The three routes -- the GF(p) engine
(int64 cochain/chain bases), the generic bar mixed complex (Domain), and
Chouhy-Solotar (Domain, CS chain words) -- feed the SAME serializer here through a
uniform "labeled basis element" abstraction.

Float-free: every coefficient and matrix entry is an int or an exact Domain-element
string (the AST gate scans this file).
"""

# Plan-34 recorder backstop: a (co)chain enumeration list, or a differential
# matrix, larger than this many cells is ELIDED with a rebuild pointer rather than
# shipped inline. The self-cert path rebuilds the elided matrix from the pointer.
MATRIX_CELL_CAP = 250_000


def labels_of(A):
    """Basis labels of ``A`` (already unit-adapted at the call sites that need it),
    with an honest index fallback when the algebra carries no labels."""
    labs = getattr(A, "basis_labels", None)
    if labs is None:
        return ["e_%d" % i for i in range(A.dim)]
    return list(labs)


def element_label(word_labels, value_label, kind):
    """One ordered-basis element as a display string.

    cochain: ``[w1 (x) ... (x) wp -> v]`` (degree 0: just ``v``, an element of
    A = C^0); chain: ``v (x) w1 (x) ... (x) wn`` (degree 0: just ``v``). The
    separator is ALWAYS the k-tensor ``(x)`` (displayed ⊗) -- Marco 2026-08-03:
    ``|`` read as something else, and a tensor w (x) w is nonzero even when the
    product w*w vanishes in A."""
    if kind == "cochain":
        if not word_labels:
            return value_label
        return "[%s -> %s]" % (" (x) ".join(word_labels), value_label)
    if not word_labels:
        return value_label
    return " (x) ".join((value_label,) + tuple(word_labels))


# --------------------------------------------------------------------------- #
# Coefficient stringification. dom is None on the GF(p) engine (int64 coeffs);
# a Domain otherwise (generic bar / CS). Exact strings only.
# --------------------------------------------------------------------------- #
def _is_zero(coeff, dom):
    if dom is None:
        return int(coeff) == 0
    return dom.is_zero(dom.coerce(coeff))


def _coeff_str(coeff, dom):
    if dom is None:
        return str(int(coeff))
    return str(dom.coerce(coeff))


# --------------------------------------------------------------------------- #
# One class = one column vector over an ordered, labeled basis enumeration.
# `elements` is the list of (word_labels_tuple, value_label); `kind` in
# {"cochain","chain"}; `dom` selects the coefficient arithmetic.
# --------------------------------------------------------------------------- #
def serialize_class(column, elements, n, kind, dom):
    terms, vector = [], []
    for idx, coeff in enumerate(column):
        if _is_zero(coeff, dom):
            continue
        cs = _coeff_str(coeff, dom)
        word_labels, value_label = elements[idx]
        terms.append([cs, list(word_labels), value_label])
        vector.append([idx, cs])
    return {"kind": kind, "degree": n, "terms": terms, "vector": vector}


def classes_from_columns(columns, elements, n, kind, dom):
    """A list of column vectors -> the serialized class list of a degree."""
    return [serialize_class(col, elements, n, kind, dom) for col in columns]


def enumeration_labels(elements, kind):
    """The ordered basis enumeration as display labels (what the vector indices
    point into). Capped: an over-long enumeration ships an elided marker with the
    length and a rebuild pointer -- never silently truncated."""
    if len(elements) > MATRIX_CELL_CAP:
        return {"elided": True, "length": len(elements),
                "note": "enumeration exceeds the 250000-cell display cap; "
                        "rebuild from the ordered (co)chain basis"}
    return [element_label(w, v, kind) for (w, v) in elements]


# --------------------------------------------------------------------------- #
# Differential serialization. `shape` is (rows, cols) known WITHOUT building the
# matrix (so the cap can refuse an expensive build); `build` is a thunk returning
# the row-major matrix (int64 rows on the engine, Domain rows otherwise); `note`
# is the exact quiverlab API call to rebuild an elided matrix.
# --------------------------------------------------------------------------- #
def serialize_differential(shape, build, note, dom):
    rows, cols = int(shape[0]), int(shape[1])
    if rows * cols > MATRIX_CELL_CAP:
        return {"elided": True, "shape": [rows, cols], "note": note}
    M = build()
    return {"shape": [rows, cols],
            "rows": [[_coeff_str(e, dom) for e in row] for row in M]}


# --------------------------------------------------------------------------- #
# Route element builders: each returns the ordered list of
# (word_labels_tuple, value_label) for one degree, matching the SAME enumeration
# the reps index into.
# --------------------------------------------------------------------------- #
def engine_coh_elements(E, n, labels):
    """Engine cochain basis C^n: (w, j), w a tuple of engine indices, value e_j."""
    from quiverlab.engine.scan3 import cochain_basis
    return [(tuple(labels[x] for x in w), labels[j])
            for (w, j) in cochain_basis(E, n)]


def engine_hom_elements(E, n, labels):
    """Engine chain basis C_n: (a0, r1..rn) = b_{a0} (x) b_{r1} (x) ..."""
    from quiverlab.engine.hh_engine import cn_basis
    out = []
    for gen in cn_basis(E, n):
        a0, tail = gen[0], gen[1:]
        out.append((tuple(labels[x] for x in tail), labels[a0]))
    return out


def bar_chain_elements(m, n, labels):
    """Generic bar chain basis (hochschild.bar._cochain_basis): (s, J), value b_s,
    tensor word J -- the shape the generic-Domain Connes reps live over."""
    from quiverlab.hochschild.bar import _cochain_basis
    return [(tuple(labels[x] for x in J), labels[s])
            for (s, J) in _cochain_basis(m, n)]


def _word_str(word):
    """A CS chain word (a path = sequence of arrow labels, or a string) rendered as
    ONE composition string; empty path -> '' (a degree-0 generator, no word)."""
    if isinstance(word, str):
        return word
    return "·".join(word)               # arrows joined left-to-right (a·b·c)


def cs_elements(res, n, side, labels):
    """CS basis res._basis(n, side): (chain, j); Marco -- the chain's whole path
    word is ONE label (a 1-tuple), value e_j. A degree-0 generator has no word."""
    out = []
    for (ch, j) in res._basis(n, side):
        ws = _word_str(ch.word)
        out.append(((ws,) if ws else (), labels[j]))
    return out
