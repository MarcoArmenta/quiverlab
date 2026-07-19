"""Simples, projectives, injectives from the quiver presentation (spec §3.6, §5 c.7).

RIGHT modules. P_v = e_v A (basis paths STARTING at v); S_v = top P_v (1-dim at v,
arrows act 0); I_v = D(A e_v) = k-dual of the LEFT projective A e_v (basis paths ENDING
at v), right action = transpose of left multiplication. All exact over any Domain."""
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module


def _require_provenance(A, what):
    if A.quiver is None or A.basis_labels is None:
        raise QuiverlabError(
            f"{what} needs the quiver presentation",
            hint="construct the algebra via Quiver.algebra(...); structure-constant "
                 "algebras carry no path basis",
        )


def _label_vertex_source(A, label):
    if label.startswith("e_"):
        return next(v for v in A.quiver.vertices if f"e_{v}" == label)
    return A.quiver.word_source(tuple(label.split("*")))


def _label_vertex_target(A, label):
    if label.startswith("e_"):
        return next(v for v in A.quiver.vertices if f"e_{v}" == label)
    return A.quiver.word_target(tuple(label.split("*")))


def simple(A, v):
    _require_provenance(A, "simple(v)")
    dom = A.domain
    o, z = dom.one(), dom.zero()
    action = {}
    for w in A.quiver.vertices:
        action[f"e_{w}"] = [[o if w == v else z]]
    for label in A.basis_labels:
        if not label.startswith("e_"):
            action[label] = [[z]]        # rad acts as 0 on a simple
    return Module(A, 1, action, name=f"S_{v}")


def projective(A, v):
    """P_v = e_v A. Basis = the algebra basis labels whose path STARTS at v, in the
    algebra's basis order. Right action of a basis element b: right multiplication,
    read from the algebra's structure constants restricted to this sub-basis."""
    _require_provenance(A, "projective(v)")
    dom = A.domain
    # sub-basis: indices of basis labels starting at v (idempotent e_v included)
    sub = [i for i, lab in enumerate(A.basis_labels) if _label_vertex_source(A, lab) == v]
    pos = {gi: k for k, gi in enumerate(sub)}
    n = len(sub)
    action = {}
    for blab in A.basis_labels:
        bi = A.basis_labels.index(blab)
        Mb = lm.zeros(n, n, dom)               # column k = coords of (p_k * b)
        for k, gi in enumerate(sub):
            prod = A.multiply(A._basis_vec(gi), A._basis_vec(bi))  # p_k * b in A-coords
            for gj in sub:
                Mb[pos[gj]][k] = prod[gj]
        action[blab] = Mb
    mod = Module(A, n, action, name=f"P_{v}")
    # Representation hooks for the resolution engine (private surface, leading
    # underscore): the ordered basis labels (paths starting at v) and this summand's
    # vertex. The canonical generator of P_v is its e_v basis vector = the trivial
    # path, the first entry of _pv_basis_labels.
    mod._pv_basis_labels = [A.basis_labels[i] for i in sub]
    mod._pv_vertex = v
    return mod


def injective(A, v):
    """I_v = D(A e_v). Left projective A e_v has basis the labels whose path ENDS at v.
    Left multiplication L_b (b * p) gives the left action; the right action on the dual
    is its transpose. dim e_w I_v = # basis paths w -> v = C[w][v]."""
    _require_provenance(A, "injective(v)")
    dom = A.domain
    sub = [i for i, lab in enumerate(A.basis_labels) if _label_vertex_target(A, lab) == v]
    pos = {gi: k for k, gi in enumerate(sub)}
    n = len(sub)
    action = {}
    for blab in A.basis_labels:
        bi = A.basis_labels.index(blab)
        Lb = lm.zeros(n, n, dom)               # column k = coords of (b * p_k)
        for k, gi in enumerate(sub):
            prod = A.multiply(A._basis_vec(bi), A._basis_vec(gi))  # b * p_k
            for gj in sub:
                Lb[pos[gj]][k] = prod[gj]
        action[blab] = lm.transpose(Lb)        # right action = transpose of left mult
    return Module(A, n, action, name=f"I_{v}")
