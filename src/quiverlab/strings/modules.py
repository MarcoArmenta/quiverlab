"""String and band module materialisation (Plan 46 / C5).

Right modules built in the WALK-order basis (one basis block per visited vertex, in
the order the walk visits them) and validated by ``Module.check_module`` -- the
self-certificate. We build the idempotent projections ourselves rather than route
through ``Module.from_arrow_action``, because that classmethod orders the ambient
basis by VERTEX BLOCK, whereas a string/band basis is ordered by WALK POSITION (the
two orders differ whenever a walk visits vertices out of vertex order or revisits a
vertex). Float-free: matrices are 0/1 + exact field scalars.

Orientation convention (ARBITRATED, not assumed -- the house "oracle decides the
orientation" pattern): for a RIGHT module the vertex idempotents force
``action[a] == P_{t(a)} @ action[a] @ P_{s(a)}`` (``Module.check_module`` (ii)).
For a DIRECT letter ``a`` at position ``i`` (``a: verts[i] -> verts[i+1]``) that
sends the source basis vector ``z_i`` to the target vector ``z_{i+1}``:
``action[a][i+1][i] = 1``. An INVERSE letter transposes the slot:
``action[a][i][i+1] = 1``. Pinned by
``test_length_one_direct_string_is_the_projective_arm`` (the walk ``(a)`` over gentle
kA3/(ab) is ``P_1 = [1;2]``, whose engine-built ``action[a]`` is ``[[0,0],[1,0]]``)
-- the transpose of the plan skeleton, flipped once per the grading oracle."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module
from quiverlab.strings.walks import (_is_trivial, is_valid_walk, letter_source,
                                     letter_target)


def _materialise(A, vertex_of_index, arrow_action, name):
    """Build + self-certify a right A-module whose ambient basis is indexed by
    ``vertex_of_index`` (basis vector ``i`` sits at vertex ``vertex_of_index[i]``),
    with the given per-arrow action matrices. Raises loudly if ``check_module`` fails."""
    dom = A.domain
    n = len(vertex_of_index)
    action = {}
    for v in A.quiver.vertices:
        P = lm.zeros(n, n, dom)
        for i, vi in enumerate(vertex_of_index):
            if vi == v:
                P[i][i] = dom.one()
        action[f"e_{v}"] = P
    for aname, mat in arrow_action.items():
        action[aname] = mat
    M = Module(A, n, action, name=name)
    M._extend_to_basis_labels()
    ok, why = M.check_module()
    if not ok:
        raise QuiverlabError(
            f"strings: {name} is not a valid module: {why}",
            hint="the walk/band does not present a module over A (a relation or the "
                 "vertex grading is violated)")
    return M


def _walk_vertices(A, walk):
    """The vertices z_0..z_n visited by ``walk`` (z_i between letters). A trivial walk
    ``((None, v),)`` yields the single vertex ``[v]``."""
    Q = A.quiver
    if _is_trivial(walk):
        return [walk[0][1]]
    verts = [letter_source(Q, walk[0])]
    for ell in walk:
        verts.append(letter_target(Q, ell))
    return verts


def _walk_name(walk):
    if _is_trivial(walk):
        return f"M(e_{walk[0][1]})"
    parts = [name if d > 0 else f"{name}^-1" for name, d in walk]
    return "M(" + " ".join(parts) + ")"


def string_module(A, walk, name=None):
    """The string module ``M(walk)``: dim = #vertices on the walk = len(walk)+1, basis
    ``z_0..z_n`` one per visited vertex, each direct/inverse letter a 0/1 partial map.
    Self-certifies via ``check_module``."""
    if A.quiver is None or A.relations is None:
        raise QuiverlabError("strings: string_module needs a quiver-presented algebra",
                             hint="build the algebra via Quiver.algebra(...)")
    if not _is_trivial(walk) and not is_valid_walk(A, walk):
        raise QuiverlabError(f"strings: {walk!r} is not a valid string over A",
                             hint="use enumerate_strings(A) to get valid walks")
    Q, dom = A.quiver, A.domain
    verts = _walk_vertices(A, walk)
    n = len(verts)
    action = {a: lm.zeros(n, n, dom) for a in Q.arrows}
    if not _is_trivial(walk):
        for i, (nm, d) in enumerate(walk):
            if nm is None:
                continue
            if d > 0:                                  # direct: z_i |-> z_{i+1}
                action[nm][i + 1][i] = dom.one()
            else:                                      # inverse: transposed slot
                action[nm][i][i + 1] = dom.one()
    return _materialise(A, verts, action, name or _walk_name(walk))


def _jordan(mult, lam, dom):
    """The Jordan block ``J_mult(lam)``: ``lam`` on the diagonal, ``1`` on the
    superdiagonal (a single block => the band monodromy is indecomposable)."""
    J = lm.zeros(mult, mult, dom)
    for i in range(mult):
        J[i][i] = lam
    for i in range(mult - 1):
        J[i][i + 1] = dom.one()
    return J


def band_module(A, walk, eigenvalue, mult=1, name=None):
    """The band module ``M(walk, lambda, mult)``: a cyclic string of length ``L``;
    dim = ``mult * L``. Identity ``mult x mult`` blocks for all letters but the closing
    one, which carries ``J_mult(lambda)`` -- the loop monodromy is a single Jordan
    block, so the module is indecomposable for any nonzero ``lambda`` in the field.
    ``eigenvalue`` MUST be a nonzero field element (coerced via ``A.domain``) -- loud
    otherwise. Self-certifies via ``check_module``."""
    if A.quiver is None or A.relations is None:
        raise QuiverlabError("strings: band_module needs a quiver-presented algebra",
                             hint="build the algebra via Quiver.algebra(...)")
    if mult < 1:
        raise QuiverlabError("strings: band multiplicity must be >= 1")
    dom = A.domain
    try:
        lam = dom.coerce(eigenvalue)
    except Exception as exc:
        raise QuiverlabError(
            f"strings: band eigenvalue {eigenvalue!r} is not in the field",
            hint="the eigenvalue must be a nonzero element of A.domain") from exc
    if dom.is_zero(lam):
        raise QuiverlabError("strings: band eigenvalue must be nonzero")
    # primitive-band guard (devil's-advocate fix, 2026-08-05): a proper power
    # b^k materialises as a DECOMPOSABLE module -- silently wrong under the
    # docstring's indecomposability promise. Refuse loudly instead.
    from quiverlab.strings.walks import _is_band_walk, _is_proper_power
    if _is_proper_power(tuple(walk)):
        raise QuiverlabError(
            "strings: the walk is a proper power of a shorter band -- a band "
            "module needs a PRIMITIVE band",
            hint="pass the primitive band once and use mult= for multiplicity")
    if not _is_band_walk(A, tuple(walk)):
        raise QuiverlabError(
            "strings: the walk is not a band (cyclic, reduced, composable, "
            "relation-avoiding at the closure)",
            hint="find_bands(A) enumerates the primitive bands")
    Q = A.quiver
    L = len(walk)
    if L < 1:
        raise QuiverlabError("strings: a band walk must have length >= 1")
    verts = [letter_source(Q, walk[k]) for k in range(L)]
    dim = L * mult
    vertex_of_index = [verts[k] for k in range(L) for _ in range(mult)]
    action = {a: lm.zeros(dim, dim, dom) for a in Q.arrows}
    for k, (nm, d) in enumerate(walk):
        if d > 0:                                      # arrow: block k -> block (k+1)
            src_blk, tgt_blk = k, (k + 1) % L
        else:                                          # inverse letter
            src_blk, tgt_blk = (k + 1) % L, k
        B = _jordan(mult, lam, dom) if k == L - 1 else lm.identity(mult, dom)
        for r in range(mult):
            for c in range(mult):
                action[nm][tgt_blk * mult + r][src_blk * mult + c] = B[r][c]
    return _materialise(A, vertex_of_index, action, name or "band")
