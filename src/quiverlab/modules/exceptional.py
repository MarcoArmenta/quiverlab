"""Classical exceptional sequences of a hereditary algebra (Plan 65 / R28).

For a **hereditary** ``A = kQ`` (``Q`` acyclic, ``I = 0``, so ``gl.dim <= 1``):

- an **exceptional module** ``M`` has ``End_A(M)`` a division ring and
  ``Ext^{>=1}_A(M, M) = 0``. Over a hereditary algebra only ``Ext^1`` can be nonzero.
  We implement the **brick** criterion ``end_dim(M) = 1 and Ext^1(M, M) = 0`` (M3): a
  1-dimensional endomorphism ring is ``k`` -- a field, hence a division ring -- so this is
  never a FALSE positive; over a non-algebraically-closed field a module can be exceptional
  with a *larger* division-ring ``End`` (``end_dim > 1``), which the brick test returns
  ``False`` for (a scope-limited false negative that does NOT occur on the Dynkin/QQ
  battery -- every Dynkin indecomposable is a rigid brick). Over algebraically closed ``k``,
  brick <=> exceptional exactly.
- an **exceptional sequence** ``(E_1, ..., E_r)`` is a tuple of exceptional modules with, for
  every ``i < j``, ``Hom_A(E_j, E_i) = 0`` AND ``Ext^1_A(E_j, E_i) = 0`` (no maps or
  extensions BACKWARD, from a later term to an earlier one -- the Crawley-Boevey / Ringel
  convention; the opposite convention merely reverses each sequence, same count). Complete
  iff ``r = n = #simples``.
- the **braid group** ``B_n`` acts on complete exceptional sequences by mutation ``sigma_i``
  (Task B2), transitively (Crawley-Boevey Ottawa 1992; Ringel Contemp. Math. 171 (1994)).
- for **Dynkin** ``A`` the set of complete exceptional sequences is finite, of size
  ``#CES = n! * h^n / |W|`` (Obaid et al.); ``A_n = (n+1)^{n-1}``, ``D_4 = 162``.

Scope (contractual, loud typed refusals): **hereditary only** (non-hereditary refuses);
enumeration is **representation-finite (Dynkin) only** (rep-infinite hereditary has an
infinite braid orbit -> honest ``status`` cap). Over ``GF(p)`` the shipped
``decompose``/``is_isomorphic`` char caveat (char 0 or char > dim) applies. Exact only --
c-matrices are integer matrices, counts ``int``, verdicts ``bool``/``str``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from quiverlab.errors import QuiverlabError
from quiverlab.modules.ext import ext
from quiverlab.modules.hom import end_dim, hom_dim


def _require_hereditary(A):
    """Refuse loudly unless ``A`` is a hereditary path algebra ``kQ`` (``I = 0``). Delegates
    to :func:`invariants.recognizers.is_hereditary`, which also refuses a structure-constant-
    only algebra (no quiver)."""
    from quiverlab.invariants.recognizers import is_hereditary
    if not is_hereditary(A):
        raise QuiverlabError(
            "classical exceptional sequences require a hereditary algebra kQ (I = 0); "
            "this algebra has relations (gl.dim >= 2) or no quiver presentation",
            hint="the tau-exceptional surface (tau_exceptional_sequences) covers "
                 "non-hereditary tau-tilting-finite algebras")


# --------------------------------------------------------------------------- #
# recognizers
# --------------------------------------------------------------------------- #
def is_exceptional_module(A, M) -> bool:
    """True iff ``M`` is a brick with no self-extension (the BRICK criterion, M3):
    ``end_dim(M) == 1 and Ext^1_A(M, M) == 0``. Over a hereditary algebra ``gl.dim <= 1``, so
    higher self-Ext vanish automatically and ``Ext^1`` is the only obstruction. Coincides
    with "exceptional" over algebraically closed ``k`` and on the Dynkin/QQ battery."""
    _require_hereditary(A)
    return end_dim(M) == 1 and ext(A, M, M, 1) == 0


def is_exceptional_sequence(A, seq) -> bool:
    """True iff ``seq = (E_1, ..., E_r)`` is an exceptional sequence: each ``E_i`` is an
    exceptional module and, for every ``i < j``, ``Hom_A(E_j, E_i) = 0`` and
    ``Ext^1_A(E_j, E_i) = 0`` (no backward maps/extensions -- the fixed convention)."""
    _require_hereditary(A)
    seq = list(seq)
    for E in seq:
        if not (end_dim(E) == 1 and ext(A, E, E, 1) == 0):
            return False
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if hom_dim(seq[j], seq[i]) != 0:
                return False
            if ext(A, seq[j], seq[i], 1) != 0:
                return False
    return True


def c_matrix(A, seq):
    """The c-matrix of an exceptional sequence: rows = the dimension vectors of the terms in
    vertex order (integers). For a MODULE (unshifted) sequence every row is the non-negative
    dim-vector of ``E_i``; this is the hereditary exceptional-sequence reading of the P45
    c-vectors (wall normals / brick dim-vectors). Hereditary-only."""
    _require_hereditary(A)
    verts = list(A.quiver.vertices)
    rows = []
    for E in seq:
        dv = E.dimension_vector()
        rows.append([int(dv[v]) for v in verts])
    return rows


# --------------------------------------------------------------------------- #
# Task B2: braid mutation sigma_i (and sigma_i^{-1}) by the canonical constructions
# --------------------------------------------------------------------------- #
_UNIVERSE_CACHE = {}


def _universe(A):
    """The knitted indecomposable universe of a rep-finite hereditary ``A`` (cached by
    ``id(A)``). Refuses loudly if ``A`` is representation-INFINITE (infinite braid orbit --
    mutation identification is Dynkin-only)."""
    key = id(A)
    u = _UNIVERSE_CACHE.get(key)
    if u is None:
        arq = A.ar_quiver()
        if not arq.is_complete:
            raise QuiverlabError(
                "braid_mutation: the AR quiver did not close (representation-infinite "
                "hereditary algebra -- the braid orbit is infinite); mutation "
                "identification is Dynkin (rep-finite) only")
        u = [v["module"] for v in arq.vertices]
        _UNIVERSE_CACHE[key] = u
    return u


def _canonicalize(A, L, universe):
    """Return the canonical universe representative isomorphic to the constructed module
    ``L`` (rep-finite Dynkin cross-identification, Task B2). Asserts a match exists --
    a self-cert that catches a construction bug (a wrong ``L`` matches NO universe module)."""
    from quiverlab.modules.hom import is_isomorphic
    if universe is None:
        universe = _universe(A)
    for u in universe:
        if u.dim == L.dim and u.dimension_vector() == L.dimension_vector() \
                and is_isomorphic(u, L):
            return u
    raise QuiverlabError(
        "braid_mutation: the constructed mutation module matches no indecomposable in the "
        "knitted universe -- the categorical construction is inconsistent",
        hint="report the algebra + sequence + position; this is a construction-side bug, "
             "never returned as a silent wrong module")


def _evaluation_mutation(A, X, Y, e):
    """The evaluation fibre ``L`` of ``ev: Hom(X,Y) (x) X -> Y`` (``e = dim Hom(X,Y) > 0``):
    ``L = ker(ev)`` when ``ev`` is surjective, ``L = coker(ev)`` when injective (Ringel's
    dual sub-case -- H4's socle-inclusion example)."""
    from quiverlab.modules.morphism import ModuleHom, direct_sum, hom_basis
    dom = A.domain
    hb = hom_basis(X, Y)
    Xe, _incls, _projs = direct_sum(*([X] * e))
    ev_mat = [[dom.zero()] * (e * X.dim) for _ in range(Y.dim)]
    for j, f in enumerate(hb):
        for a in range(Y.dim):
            fa = f.matrix[a]
            for b in range(X.dim):
                ev_mat[a][j * X.dim + b] = fa[b]
    ev = ModuleHom(Xe, Y, ev_mat, check=True)
    if ev.is_epi():
        L, _iota = ev.kernel()
        return L
    if ev.is_mono():
        C, _proj = ev.cokernel()
        return C
    raise QuiverlabError(
        "braid_mutation: the evaluation map is neither surjective nor injective "
        "(outside the hereditary rep-finite scope)")


def _coevaluation_mutation(A, X, Y, e):
    """The co-evaluation fibre ``R`` of ``coev: X -> Hom(X,Y)^* (x) Y`` (the dual of
    :func:`_evaluation_mutation`; ``e = dim Hom(X,Y) > 0``): ``R = coker(coev)`` when injective,
    ``R = ker(coev)`` when surjective."""
    from quiverlab.modules.morphism import ModuleHom, direct_sum, hom_basis
    dom = A.domain
    hb = hom_basis(X, Y)
    Ye, _incls, _projs = direct_sum(*([Y] * e))
    coev_mat = [[dom.zero()] * X.dim for _ in range(e * Y.dim)]
    for j, f in enumerate(hb):
        for a in range(Y.dim):
            fa = f.matrix[a]
            for b in range(X.dim):
                coev_mat[j * Y.dim + a][b] = fa[b]
    coev = ModuleHom(X, Ye, coev_mat, check=True)
    if coev.is_mono():
        C, _proj = coev.cokernel()
        return C
    if coev.is_epi():
        K, _iota = coev.kernel()
        return K
    raise QuiverlabError(
        "braid_mutation: the co-evaluation map is neither injective nor surjective "
        "(outside the hereditary rep-finite scope)")


def _ext1_cocycles(A, X, Y):
    """A basis of ``Ext^1_A(X, Y)`` as explicit cocycle MATRICES ``f: P_1(X) -> Y`` (each
    ``Y.dim x P_1(X).dim``), for the universal-extension constructions. Cocycles
    ``Z^1 = {phi in Hom(P_1, Y) : phi . d_2 = 0}`` modulo coboundaries
    ``B^1 = {psi . d_1 : psi in Hom(P_0, Y)}`` (the standard ``H^1`` of ``Hom(P_*, Y)``)."""
    from quiverlab.modules import linalg_mod as lm
    from quiverlab.modules.hom import hom_space
    from quiverlab.modules.resolution import minimal_resolution
    dom = A.domain
    terms, dmats = minimal_resolution(X, 2)
    P1 = terms[1].module if len(terms) > 1 else None
    if P1 is None or P1.dim == 0:
        return []
    H1 = hom_space(P1, Y)
    if not H1:
        return []

    def flat(m):
        return [x for row in m for x in row]

    B1 = lm.cols_to_matrix([flat(h) for h in H1])          # coord frame for Hom(P_1, Y)
    # cocycles: kernel of phi -> phi . d_2 (in H1-coords)
    d2 = dmats[2] if len(dmats) > 2 else None
    P2 = terms[2].module if len(terms) > 2 else None
    if d2 and d2[0] and P2 is not None and P2.dim:
        cols = [flat(lm.matmul(h, d2, dom)) for h in H1]
        cocyc = lm.kernel_columns(lm.cols_to_matrix(cols), dom)
    else:
        cocyc = [lm.col(lm.identity(len(H1), dom), j) for j in range(len(H1))]
    # coboundaries: image of psi -> psi . d_1 (expressed in H1-coords)
    P0 = terms[0].module
    d1 = dmats[1] if len(dmats) > 1 else None
    cobound = []
    if d1 and d1[0] and P0 is not None and P0.dim:
        for psi in hom_space(P0, Y):
            comp = lm.matmul(psi, d1, dom)
            x = lm.solve_columns(B1, lm.cols_to_matrix([flat(comp)]), dom)
            if x is not None:
                cobound.append(x[0])
    chosen = lm.independent_modulo(cocyc, cobound, dom) if cocyc else []
    out = []
    for j in chosen:
        c = cocyc[j]
        M = lm.zeros(Y.dim, P1.dim, dom)
        for k, ck in enumerate(c):
            if dom.is_zero(ck):
                continue
            Hk = H1[k]
            for a in range(Y.dim):
                Ma = M[a]
                Hka = Hk[a]
                for b in range(P1.dim):
                    Ma[b] = dom.add(Ma[b], dom.mul(ck, Hka[b]))
        out.append(M)
    return out


def _universal_extension_middle(A, X, Y, e):
    """The universal-extension middle ``L`` of ``0 -> Y -> L -> X^e -> 0`` (``e = dim
    Ext^1(X,Y) > 0``): the diagonal class in ``Ext^1(X^e, Y)``. Self-certifies exactness."""
    from quiverlab.modules.morphism import direct_sum
    from quiverlab.modules.yoneda import baer_extension
    dom = A.domain
    basis = _ext1_cocycles(A, X, Y)
    if len(basis) != e:
        raise QuiverlabError("braid_mutation: Ext^1 cocycle basis size mismatch (universal)")
    if e == 1:
        M, cocycle = X, basis[0]
    else:
        M, _i, _p = direct_sum(*([X] * e))
        w = len(basis[0][0])                               # P_1(X).dim
        cocycle = [[dom.zero()] * (e * w) for _ in range(Y.dim)]
        for j, bj in enumerate(basis):
            for a in range(Y.dim):
                for b in range(w):
                    cocycle[a][j * w + b] = bj[a][b]
    seq = baer_extension(M, Y, cocycle)
    seq.assert_exact()
    return seq.middle


def _couniversal_extension_middle(A, X, Y, e):
    """The co-universal-extension middle ``R`` of ``0 -> Y^e -> R -> X -> 0`` (``e = dim
    Ext^1(X,Y) > 0``): the diagonal class in ``Ext^1(X, Y^e)``. Self-certifies exactness."""
    from quiverlab.modules.morphism import direct_sum
    from quiverlab.modules.yoneda import baer_extension
    dom = A.domain
    basis = _ext1_cocycles(A, X, Y)
    if len(basis) != e:
        raise QuiverlabError("braid_mutation: Ext^1 cocycle basis size mismatch (co-universal)")
    if e == 1:
        N, cocycle = Y, basis[0]
    else:
        N, _i, _p = direct_sum(*([Y] * e))
        w = len(basis[0][0])                               # P_1(X).dim
        cocycle = [[dom.zero()] * w for _ in range(e * Y.dim)]
        for j, bj in enumerate(basis):
            for a in range(Y.dim):
                for b in range(w):
                    cocycle[j * Y.dim + a][b] = bj[a][b]
    seq = baer_extension(X, N, cocycle)
    seq.assert_exact()
    return seq.middle


def _left_mutation(A, X, Y, universe):
    """``L_X Y`` (left mutation): the four canonical cases on ``(dim Hom(X,Y), dim
    Ext^1(X,Y))`` -- pure swap / evaluation fibre / universal extension / (case d) the
    ordered universal-extension-then-evaluation composite (Ringel Contemp. Math. 171)."""
    from quiverlab.modules.hom import hom_dim
    eh, ee = hom_dim(X, Y), ext(A, X, Y, 1)
    if eh == 0 and ee == 0:
        return Y                                   # pure swap (orthogonal)
    if ee == 0:
        L = _evaluation_mutation(A, X, Y, eh)
    elif eh == 0:
        L = _universal_extension_middle(A, X, Y, ee)
    else:                                          # case d: univ ext FIRST, then eval kernel
        Yt = _universal_extension_middle(A, X, Y, ee)
        eh2 = hom_dim(X, Yt)
        L = _evaluation_mutation(A, X, Yt, eh2)
    return _canonicalize(A, L, universe)


def _right_mutation(A, X, Y, universe):
    """``R_Y X`` (right mutation, ``sigma_i^{-1}``): the dual four cases -- pure swap /
    co-evaluation fibre / co-universal extension / (case d) the dual composite."""
    from quiverlab.modules.hom import hom_dim
    eh, ee = hom_dim(X, Y), ext(A, X, Y, 1)
    if eh == 0 and ee == 0:
        return X                                   # pure swap
    if ee == 0:
        R = _coevaluation_mutation(A, X, Y, eh)
    elif eh == 0:
        R = _couniversal_extension_middle(A, X, Y, ee)
    else:                                          # case d dual: co-univ ext then co-eval
        Xt = _couniversal_extension_middle(A, X, Y, ee)
        eh2 = hom_dim(Xt, Y)
        R = _coevaluation_mutation(A, Xt, Y, eh2)
    return _canonicalize(A, R, universe)


def braid_mutation(A, seq, i, *, direction="left", _universe=None):
    """The braid mutation ``sigma_i`` (``direction="left"``) or ``sigma_i^{-1}``
    (``direction="right"``) acting on positions ``(i, i+1)`` of the exceptional sequence
    ``seq`` (Plan 65 / R28; Crawley-Boevey / Ringel). With ``X = seq[i]``, ``Y = seq[i+1]``:

    - **left** ``sigma_i``: ``(..., X, Y, ...) -> (..., L_X Y, X, ...)``;
    - **right** ``sigma_i^{-1}``: ``(..., X, Y, ...) -> (..., Y, R_Y X, ...)``.

    Each new module is built by the canonical universal-extension / kernel / cokernel
    constructions, then cross-identified against the knitted rep-finite universe (a self-cert
    that catches a construction bug). Hereditary + rep-finite scope (loud otherwise)."""
    _require_hereditary(A)
    seq = list(seq)
    if not (0 <= i < len(seq) - 1):
        raise QuiverlabError(
            f"braid_mutation: position {i} out of range 0..{len(seq) - 2} for a length-"
            f"{len(seq)} sequence (sigma_i acts on (i, i+1))")
    X, Y = seq[i], seq[i + 1]
    if direction == "left":
        L = _left_mutation(A, X, Y, _universe)
        return seq[:i] + [L, X] + seq[i + 2:]
    if direction == "right":
        R = _right_mutation(A, X, Y, _universe)
        return seq[:i] + [Y, R] + seq[i + 2:]
    raise QuiverlabError(
        f'braid_mutation: direction must be "left" or "right", got {direction!r}')
