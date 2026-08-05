"""One-point extension A[M] = [[k, M], [0, A]] (Plan 44 / C7).

The triangular matrix algebra of a right A-module M read as a (k, A)-bimodule: a new
SOURCE vertex w, one arrow ``w -> v`` per basis vector of ``top(M)_v`` (the projective-cover
generators of M), and the relations forced by M's structure -- extracted by the length-lex
kernel enumeration (``families/_present.py``) and CERTIFIED per instance by
``dim A[M] == 1 + dim M + dim A``. The Cartan block form ``[[1, dim-vector M], [0, C_A]]``
(w first) and ``pd_{A[M]}(S_w) = pd_A(M) + 1`` (Happel's change-of-rings LES) are the
literature oracles. Float-free (ASS III one-point extensions)."""
from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError


def _fresh_vertex(vertices):
    """A label not among ``vertices`` for the new source vertex (an int for int-labelled
    quivers, else a ``w``-string)."""
    verts = list(vertices)
    if verts and all(isinstance(v, int) and not isinstance(v, bool) for v in verts):
        return max(verts) + 1
    existing = set(verts) | set(map(str, verts))
    cand = "w"
    while cand in existing:
        cand += "_"
    return cand


def _one_point_structure_constants(A, M):
    """A[M] as raw structure constants (dim 1 + dim M + dim A). Basis: index 0 = the k-part
    idempotent ``e_w``; ``1..dim M`` = M; ``1+dim M ..`` = A. Block multiplication
    ``(l, m, a)(l', m', a') = (l l', l m' + m.a', a a')`` (``m.a'`` the right A-action)."""
    dom = A.domain
    da, dm = A.dim, M.dim
    m_ext = 1 + dm + da
    zero, one = dom.zero(), dom.one()

    def midx(i):
        return 1 + i

    def aidx(j):
        return 1 + dm + j

    T = [[[zero] * m_ext for _ in range(m_ext)] for _ in range(m_ext)]
    T[0][0][0] = one                                    # e_w * e_w = e_w
    for i in range(dm):                                 # e_w * m_i = m_i
        T[0][midx(i)][midx(i)] = one
    labels = A.basis_labels
    for j in range(da):                                 # m_i * a_j = (col i of the right action)
        Maj = M.action[labels[j]]
        for i in range(dm):
            for t in range(dm):
                T[midx(i)][aidx(j)][midx(t)] = Maj[t][i]
    for j in range(da):                                 # a_j * a_j' in the A-part
        for jp in range(da):
            prod = A.multiply(A._basis_vec(j), A._basis_vec(jp))
            for t in range(da):
                T[aidx(j)][aidx(jp)][aidx(t)] = prod[t]
    unit = [zero] * m_ext
    unit[0] = one
    for t in range(da):
        unit[aidx(t)] = A.unit[t]
    # check=True validates associativity/unit -- the self-certificate that A[M] is a genuine
    # algebra (i.e. that M is a genuine (k, A)-bimodule).
    return Algebra.from_structure_constants(T, unit, field=dom, check=True), midx, aidx


def OnePointExtension(A, M):
    """The one-point extension ``A[M]`` as a genuine ``kQ'/I'``-presented Algebra.

    ``M`` is read as a RIGHT A-module (the (k, A)-bimodule). Certified per instance by the
    dimension identity ``dim A[M] == 1 + dim M + dim A``; loud ``QuiverlabError`` otherwise.
    """
    from quiverlab.combinat.quiver import Quiver
    from quiverlab.families._present import present_from_pi
    from quiverlab.modules.resolution import _homogeneous_top_generators
    if A.quiver is None or A.basis_labels is None:
        raise QuiverlabError(
            "OnePointExtension needs a presented base algebra A",
            hint="build A via Quiver.algebra(...); structure-constant algebras carry no "
                 "path basis to attach the new vertex to")
    if M.algebra is not A or M.side != "right":
        raise QuiverlabError(
            "OnePointExtension(A, M): M must be a RIGHT A-module over the SAME algebra A",
            hint="A[M] reads M as the (k, A)-bimodule; pass a right module over A")
    dom = A.domain
    da, dm = A.dim, M.dim
    T, midx, aidx = _one_point_structure_constants(A, M)

    w = _fresh_vertex(A.quiver.vertices)                # new source vertex, placed FIRST
    verts = [w] + list(A.quiver.vertices)
    arrows = dict(A.quiver.arrows)
    img = {}
    label_index = {lab: t for t, lab in enumerate(A.basis_labels)}
    for a in A.quiver.arrows:                           # original arrow -> A-part image
        vec = [dom.zero()] * T.dim
        vec[aidx(label_index[a])] = dom.one()
        img[a] = vec

    gens = _homogeneous_top_generators(M)               # [(vertex, column-in-M)] -- the w -> v arrows
    prefix = "n"
    while any((prefix + str(s)) in arrows for s in range(len(gens))):
        prefix += "_"
    for s, (v, gcol) in enumerate(gens):
        name = prefix + str(s)
        arrows[name] = (w, v)
        vec = [dom.zero()] * T.dim
        for i in range(dm):
            vec[midx(i)] = gcol[i]                       # pi(w->v) = the top generator in the M-part
        img[name] = vec

    Q = Quiver(verts, arrows)
    base_bound = A.loewy_length() + 2
    return present_from_pi(Q, img, T, dom, 1 + dm + da, base_bound,
                           citations=("assem_book", "happel_question"))
