"""Tilting modules + Bongartz completion (Plan 44 / C7).

A (classical, n=1) tilting module ``T`` over ``A``: (1) ``pd T <= 1``; (2)
``Ext^1(T, T) = 0``; (3) the number of pairwise non-iso indecomposable summands of ``T``
equals ``|Q_0|``. Bongartz's theorem (LNM 903, 1981; ASS VI.4): given (1)+(2), (3) is
equivalent to the existence of the add(T)-coresolution ``0 -> A -> T^0 -> T^1 -> 0``, so
the COUNT is the checkable form of the third axiom -- applied ONLY inside its theorem's
hypotheses (pd<=1 and Ext^1(T,T)=0 are tested first). For general ``n`` the count shortcut
is n=1 only; higher ``n`` verifies the add(T)-coresolution of ``A`` directly via left
add(T)-approximations (Task A). Cotilting = tilting of ``D(T)`` over ``A^op`` (Plan 24).

Bongartz completion of a partial tilting ``T`` (pd<=1, Ext^1(T,T)=0): the middle term
``E`` of the universal extension ``0 -> A -> E -> T^d -> 0``, ``d = dim Ext^1(T, A_A)``,
built as ``d`` SUCCESSIVE Baer extensions along an ``Ext^1(T, A_A)`` basis (each pushed
forward along the accumulated ``A -> E_i`` inclusion; ``Ext^1(T,T)=0`` splits every
intermediate quotient so ``E_d / A = T^d``, and the connecting-map argument gives
``Ext^1(T, E) = 0``, the universal property). SELF-CERTIFIED: ``is_tilting_module(T (+) E)``
is the arbiter. Float-free; the summand count inherits the ``decompose`` char caveat
(rigorous over char 0 or char > dim), refusing loudly off scope.
"""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.decompose import decompose
from quiverlab.modules.ext import ext_dims
from quiverlab.modules.resolution import minimal_resolution


@dataclass
class TiltingReport:
    is_tilting: bool
    n: int
    pd: object            # int or None (unresolved within bound)
    self_ext_vanishes: bool
    num_summands: int
    num_vertices: int
    note: str

    def __bool__(self):
        return self.is_tilting

    def __repr__(self):
        if self.is_tilting:
            return f"tilting module (pd {self.pd}, {self.num_summands} summands)"
        return f"not tilting: {self.note}"


def _num_noniso_summands(T, budget):
    return len(decompose(T, budget=budget))          # loud on char <= dim (decompose caveat)


def _coresolves_A_in_add_T(A, T, n, budget):
    """Third tilting axiom for ``n >= 2``: does ``A_A`` have an ``add(T)``-coresolution
    ``0 -> A -> T^0 -> ... -> T^{<=n} -> 0``? Iterate the minimal left ``add(T)``-
    approximation and its cokernel; the coresolution exists (length <= n) iff the
    cokernel vanishes within ``n`` steps. Each approximation is certified or raises
    (never a silent False); a non-injective approximation means no coresolution."""
    from quiverlab.modules.approximations import left_add_approximation
    from quiverlab.modules.morphism import direct_sum
    reg, _, _ = direct_sum(*[A.projective(v) for v in A.quiver.vertices])
    cur = reg
    for _step in range(n + 1):
        if cur.dim == 0:
            return True
        g = left_add_approximation(T, cur, budget=budget)   # cur >-> T^k (minimal)
        if not g.is_mono():
            return False                             # A does not embed in add T here
        cur, _proj = g.cokernel()
    return cur.dim == 0


def is_tilting_module(T, n=1, bound=64, budget=512):
    """A :class:`TiltingReport` for whether ``T`` is an ``n``-tilting module. ``__bool__``
    is ``is_tilting``; ``note`` names the first failing clause. n=1 (classical) applies the
    Bongartz count criterion only after pd<=1 and Ext^1(T,T)=0 hold."""
    A = T.base_algebra
    verts = list(A.quiver.vertices)
    pd = T.projective_resolution(bound).pd()
    pd_ok = pd is not None and pd <= n
    ext = ext_dims(A, T, T, n)                        # [dim Ext^0, ..., dim Ext^n]
    self_ext = all(d == 0 for d in ext[1:n + 1])
    k = _num_noniso_summands(T, budget)
    if n == 1:
        third = (k == len(verts))
        note_third = f"# summands {k} != # vertices {len(verts)}"
    else:
        third = _coresolves_A_in_add_T(A, T, n, budget)
        note_third = f"A has no add(T)-coresolution of length <= {n}"
    ok = bool(pd_ok and self_ext and third)
    note = ("ok" if ok else
            (f"pd {pd} > {n}" if not pd_ok else
             ("Ext^i(T,T) != 0 for some 1<=i<=n" if not self_ext else note_third)))
    return TiltingReport(ok, n, pd, self_ext, k, len(verts), note)


def is_cotilting_module(T, n=1, bound=64, budget=512):
    """``T`` is ``n``-cotilting over ``A`` iff ``D(T)`` is ``n``-tilting over ``A^op``
    (Plan 24: ``D`` flips the side over the same base algebra)."""
    return is_tilting_module(T.dualize(), n=n, bound=bound, budget=budget)


def _ext1_cocycles(T, N, terms, dmats):
    """Cocycle representatives (matrices ``f: P_1(T) -> N``, ``N.dim x P_1.dim``) of a
    basis of ``Ext^1_A(T, N)``: each ``f`` satisfies ``f @ d_2 = 0`` (a genuine cocycle)
    and the classes are independent modulo the coboundaries ``im(delta^0)``. Uses the
    SAME minimal resolution ``(terms, dmats)`` the caller passes to ``baer_extension``, so
    the ``P_1`` basis lines up exactly."""
    from quiverlab.modules.ext import _delta_matrix
    from quiverlab.modules.hom import hom_space
    dom = T.domain
    P0 = terms[0].module
    P1 = terms[1].module
    if P1 is None or P1.dim == 0:
        return []                                    # T projective: Ext^1(T, -) = 0
    P2 = terms[2].module if len(terms) > 2 else None
    H0 = hom_space(P0, N) if (P0 is not None and P0.dim) else []
    H1 = hom_space(P1, N)
    H2 = hom_space(P2, N) if (P2 is not None and P2.dim) else []
    d1 = dmats[1]                                    # P1 -> P0
    d2 = dmats[2] if len(dmats) > 2 else None        # P2 -> P1
    delta0 = (_delta_matrix(H0, H1, d1, dom) if (d1 and d1[0])
              else lm.zeros(len(H1), len(H0), dom))
    delta1 = (_delta_matrix(H1, H2, d2, dom) if (d2 and d2[0])
              else lm.zeros(len(H2), len(H1), dom))
    if delta1 and delta1[0]:
        cocyc = lm.kernel_columns(delta1, dom)
    else:                                            # no P2 constraint: all of H1 are cocycles
        cocyc = [lm.col(lm.identity(len(H1), dom), j) for j in range(len(H1))]
    bdys = ([lm.col(delta0, j) for j in range(len(H0))]
            if (delta0 and delta0[0]) else [])
    chosen = lm.independent_modulo(cocyc, bdys, dom)
    out = []
    for coords in (cocyc[i] for i in chosen):
        f = lm.zeros(N.dim, P1.dim, dom)
        for idx, c in enumerate(coords):
            if dom.is_zero(c):
                continue
            h = H1[idx]
            for r in range(N.dim):
                fr, hr = f[r], h[r]
                for j in range(P1.dim):
                    fr[j] = dom.add(fr[j], dom.mul(c, hr[j]))
        out.append(f)
    return out


def _assert_partial_tilting(T, bound):
    """Bongartz precondition: ``T`` is a PARTIAL tilting module -- ``pd T <= 1`` (exact)
    and ``Ext^1(T, T) = 0``. These are exactly the ``is_tilting_module`` clauses MINUS the
    summand count. Raises :class:`QuiverlabError` loudly otherwise (the Kronecker-regular
    module with ``Ext^1(T, T) != 0``, and any ``pd > 1`` or unresolved ``pd``, are
    refused rather than fed into the universal-extension machinery)."""
    A = T.base_algebra
    pd = T.projective_resolution(bound).pd()
    if pd is None:
        raise QuiverlabError(
            "bongartz_completion needs a partial tilting module: pd T is unresolved "
            f"within bound {bound} (not certified <= 1)")
    if pd > 1:
        raise QuiverlabError(
            f"bongartz_completion needs a partial tilting module: pd T = {pd} is not <= 1")
    ext1 = ext_dims(A, T, T, 1)[1]                    # dim Ext^1(T, T)
    if ext1 != 0:
        raise QuiverlabError(
            "bongartz_completion needs a partial tilting module: "
            f"Ext^1(T, T) has dimension {ext1} != 0 (T has a self-extension)")


def bongartz_completion(T, bound=64):
    """The Bongartz complement middle term ``E`` (see the module docstring). Returns a
    :class:`~quiverlab.modules.module.Module`. PRECONDITION: ``T`` is partial tilting
    (``pd T <= 1`` and ``Ext^1(T, T) = 0``) -- asserted up front, refused loudly
    otherwise. SELF-CERTIFICATE: ``is_tilting_module(direct_sum(T, E))`` is asserted True
    before returning (never returns an ``E`` that fails to complete ``T``).
    ``d = dim Ext^1(T, A_A) = 0`` (e.g. ``T`` projective) gives ``E = A_A``."""
    from quiverlab.modules.morphism import direct_sum
    from quiverlab.modules.yoneda import baer_extension
    A = T.base_algebra
    dom = T.domain
    _assert_partial_tilting(T, bound)                # precondition (Medium 2)
    reg, _, _ = direct_sum(*[A.projective(v) for v in A.quiver.vertices])
    terms, dmats = minimal_resolution(T, 2)
    xis = _ext1_cocycles(T, reg, terms, dmats)       # list of matrices P_1(T) -> reg
    if not xis:
        E = reg                                      # already tilting-adjacent: E = A_A
    else:
        cur = reg
        a_to_cur = lm.identity(reg.dim, dom)         # module map A_A -> cur (cur.dim x reg.dim)
        for xi in xis:
            pushed = lm.matmul(a_to_cur, xi, dom)    # cur.dim x P_1.dim : P_1(T) -> cur
            seq = baer_extension(T, cur, pushed, terms=terms, dmats=dmats)
            iota = seq.maps[0]                       # cur -> E (E.dim x cur.dim)
            a_to_cur = lm.matmul(iota, a_to_cur, dom)  # A_A -> cur -> E
            cur = seq.middle
        E = cur
    TE, _, _ = direct_sum(T, E)                       # self-certificate (Medium 2)
    rep = is_tilting_module(TE)
    if not rep.is_tilting:
        raise QuiverlabError(
            "bongartz_completion self-certificate failed: T (+) E is not a tilting "
            f"module ({rep.note})")
    return E
