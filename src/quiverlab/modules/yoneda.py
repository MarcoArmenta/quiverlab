"""Yoneda exact sequences realizing ``Ext^n_A(M, N)`` classes (Plan 35 wave 3c).

Marco: "we know that for ext the bases ARE exact sequences, and we need to SHOW
these exact sequences." A class of ``Ext^n_A(M, N)`` is not merely a number: it is
(the class of) an ``n``-fold exact sequence

    0 -> N -> Q -> P_{n-2} -> ... -> P_0 -> M -> 0

with ``M`` on the right, ``N`` on the left, and the middle a genuine module. This file
CONSTRUCTS that sequence from the captured cocycle and SELF-CERTIFIES its exactness.


Derivation (the one genuinely new algorithm here)
=================================================
Fix the minimal projective resolution of the right ``A``-module ``M``:

    ... --d_2--> P_1 --d_1--> P_0 --eps--> M -> 0,   Omega^k M := im(d_k) subset P_{k-1}.

A degree-``n`` Ext class is a cocycle ``f: P_n -> N`` in ``Hom_A(P_*, N)``, i.e. an
``A``-module map with ``delta^n f = f . d_{n+1} = 0`` -- so ``f`` VANISHES on
``ker(d_n) = im(d_{n+1})``. Because ``d_n: P_n ->> Omega^n M`` is surjective with kernel
exactly ``ker(d_n)``, ``f`` FACTORS through it: there is a unique module map

    phi : Omega^n M -> N   with   phi . (d_n : P_n ->> Omega^n M) = f.                (*)

Concretely: pick a basis of ``Omega^n M = im(d_n)`` as the PIVOT COLUMNS of ``d_n`` at
column indices ``p_1, ..., p_r`` -- then the basis vector ``w_j = d_n(e_{p_j})`` has the
standard-basis preimage ``e_{p_j} in P_n``, so ``phi(w_j) = f(e_{p_j}) = col_{p_j}(f)``.
That ``phi`` agrees on a basis with the canonical map of (*), hence equals it (both
linear), and is a module map iff ``f`` is a genuine cocycle -- which we re-certify.

Now form the PUSHOUT of ``P_{n-1} <--i-- Omega^n M --phi--> N`` (``i`` the inclusion):

              i
    Omega^n M ---> P_{n-1}
        |             |
     phi|             | j
        v             v
        N  --------->  Q            Q := (P_{n-1} (+) N) / im(psi),
              iota

where ``psi(x) = (i(x), -phi(x))``. Since ``i`` and ``phi`` are module maps, ``psi`` is a
module map, so ``im(psi)`` is a SUBMODULE and ``Q`` is a genuine module. Set

    iota : N -> Q,        n |-> [(0, n)]            (the pushout leg),
    del  : Q -> P_{n-2},  [(p, n)] |-> d_{n-1}(p)   (n >= 2; well-defined: d_{n-1} kills
                                                     i(Omega^n M) = im d_n = ker d_{n-1}),
    eps  : Q -> M,        [(p, n)] |-> eps(p)       (n == 1 -- the Baer extension end).

The spliced sequence is ``iota``, ``del``, then the ORIGINAL resolution differentials
``d_{n-2}, ..., d_1`` and ``eps = d_0`` -- exact at every joint:

  * at ``N``:      ``iota`` is injective (``i`` injective => ``psi`` injective);
  * at ``Q``:      ``ker(del) = im(iota)`` -- the standard pushout computation
                   ``[(i(x), m)] = [(0, m + phi(x))] = iota(...)``;
  * at ``P_{n-2}``: ``im(del) = im(d_{n-1}) = ker(d_{n-2})``;
  * at ``P_{n-3}..P_0``: the resolution's own exactness ``im d_{k+1} = ker d_k``;
  * at ``M``:      ``eps`` is surjective.

For ``n == 1`` there are no ``P``'s in the middle: ``0 -> N -> Q -> M -> 0`` (the Baer
extension), ``Q`` the pushout, exact at ``N`` / ``Q`` / ``M`` by the same three facts.

Nothing above assumes ``M``/``N`` monomial, quadratic, or over a prime field: it is
exact linear algebra over the module's Domain. Float-free (the AST gate scans this file).
"""
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module
from quiverlab.modules.resolution import minimal_resolution


# --------------------------------------------------------------------------- #
# Direct sum of two arbitrary modules over the same algebra (block-diagonal on
# every shared action label). `resolutions._direct_sum` wants each summand to carry
# a `_pv_vertex` (projective summands only); N here is arbitrary, so we roll our own.
# --------------------------------------------------------------------------- #
def _direct_sum2(P, N):
    A = P.algebra
    dom = A.domain
    if N.algebra is not A:
        raise QuiverlabError("yoneda: the two summands are over different algebras")
    n = P.dim + N.dim
    action = {}
    for label, Pb in P.action.items():
        Nb = N.action.get(label)
        if Nb is None:
            raise QuiverlabError(
                "yoneda: module N is missing the action of %r that P carries; the two "
                "summands must share the algebra's basis labels" % label)
        M = lm.zeros(n, n, dom)
        for i in range(P.dim):
            for j in range(P.dim):
                M[i][j] = Pb[i][j]
        for i in range(N.dim):
            for j in range(N.dim):
                M[P.dim + i][P.dim + j] = Nb[i][j]
        action[label] = M
    return Module(A, n, action, name="P(+)N", side=P.side)


# --------------------------------------------------------------------------- #
# Quotient of an ambient module by an A-stable subspace, WITH the two structural maps
# the pushout needs: the projection ambient ->> Q and a section (representative lift)
# Q -> ambient. Self-contained (mirrors radtopsoc.quotient's rep selection) so the
# maps are guaranteed consistent with the returned module.
# --------------------------------------------------------------------------- #
def _quotient_with_maps(ambient, sub_cols, dom, name="Q"):
    n_amb = ambient.dim
    ident = lm.identity(n_amb, dom)
    std = [lm.col(ident, j) for j in range(n_amb)]
    rep_idx = lm.independent_modulo(std, sub_cols, dom)
    reps = [std[i] for i in rep_idx]
    n = len(reps)
    whole = [list(c) for c in sub_cols] + reps          # a basis of the ambient space
    W = lm.cols_to_matrix(whole) if whole else lm.zeros(n_amb, 0, dom)
    s = len(sub_cols)

    def coords(vec):
        sol = lm.solve_columns(W, lm.cols_to_matrix([vec]), dom)
        if sol is None:                                  # unreachable: W spans the ambient
            raise QuiverlabError("yoneda: internal basis W failed to span the ambient")
        return sol[0]

    action = {}
    for label, Ab in ambient.action.items():
        cols = []
        for r in reps:
            img = lm.matvec(Ab, r, dom)
            cols.append(coords(img)[s:])                 # drop the submodule part
        action[label] = lm.cols_to_matrix(cols) if cols else lm.zeros(n, n, dom)
    Q = Module(ambient.algebra, n, action, name=name, side=ambient.side)

    proj = lm.cols_to_matrix([coords(std[j])[s:] for j in range(n_amb)]) \
        if n else lm.zeros(0, n_amb, dom)                # ambient -> Q
    lift = lm.cols_to_matrix(reps) if reps else lm.zeros(n_amb, 0, dom)   # Q -> ambient
    return Q, proj, lift


# --------------------------------------------------------------------------- #
# Module-map checks and the exactness certifier.
# --------------------------------------------------------------------------- #
def _generators(A):
    return [f"e_{v}" for v in A.quiver.vertices] + list(A.quiver.arrows)


def _is_module_map(f, src, tgt, dom):
    """f (tgt.dim x src.dim) is an A-module map src -> tgt iff for every generator b
    tgt.action[b] @ f == f @ src.action[b]."""
    for b in _generators(src.algebra):
        lhs = lm.matmul(tgt.action[b], f, dom)
        rhs = lm.matmul(f, src.action[b], dom)
        if not _eq(lhs, rhs, dom):
            return False
    return True


def _eq(A, B, dom):
    if len(A) != len(B):
        return False
    for ra, rb in zip(A, B):
        if len(ra) != len(rb) or any(not dom.is_zero(dom.sub(x, y)) for x, y in zip(ra, rb)):
            return False
    return True


def _is_zero(M, dom):
    return all(dom.is_zero(x) for row in M for x in row)


def _rank(M, dom):
    return lm.mat_rank(M, dom) if (M and M[0]) else 0


class YonedaSequence:
    """An ``n``-fold exact sequence ``0 -> modules[0] -> ... -> modules[L] -> 0``
    realizing an ``Ext^n_A(M, N)`` class. ``modules[0] = N``, ``modules[-1] = M``,
    ``maps[i] : modules[i] -> modules[i+1]``; ``middle`` is the constructed pushout
    module ``Q`` (index 1). ``roles[i]`` labels each module (``"sub"``, ``"middle"``,
    ``"resolution_term"``, ``"quotient"``)."""

    def __init__(self, degree, modules, maps, roles, term_index, dom):
        self.degree = degree
        self.modules = modules
        self.maps = maps
        self.roles = roles
        self.term_index = term_index          # resolution-term index of each module (or None)
        self.dom = dom
        self.middle = modules[1]

    # -- self-certification --------------------------------------------------- #
    def check_exact(self):
        """Return ``(True, facts)`` when the sequence is exact and every map is an
        ``A``-module map, else ``(False, reason)``. ``facts`` is a per-node list of the
        verified rank identities (what the report states)."""
        dom = self.dom
        mods, maps = self.modules, self.maps
        for i, f in enumerate(maps):
            if not _is_module_map(f, mods[i], mods[i + 1], dom):
                return False, "map %d (%s -> %s) is not an A-module map" % (
                    i, self.roles[i], self.roles[i + 1])
        # left end: mono
        r0 = _rank(maps[0], dom)
        if r0 != mods[0].dim:
            return False, "not left-exact: %s -> %s is not injective" % (
                self.roles[0], self.roles[1])
        # right end: epi
        rl = _rank(maps[-1], dom)
        if rl != mods[-1].dim:
            return False, "not right-exact: %s -> %s is not surjective" % (
                self.roles[-2], self.roles[-1])
        facts = [{"node": self.roles[0], "fact": "injective", "rank": r0,
                  "dim": mods[0].dim}]
        # interior nodes: im(prev) = ker(next)
        for i in range(1, len(mods) - 1):
            comp = lm.matmul(maps[i], maps[i - 1], dom)
            if comp and comp[0] and not _is_zero(comp, dom):
                return False, "d.d != 0 at node %d (%s)" % (i, self.roles[i])
            r_prev = _rank(maps[i - 1], dom)
            r_next = _rank(maps[i], dom)
            if r_prev + r_next != mods[i].dim:
                return False, ("not exact at node %d (%s): rank in %d + rank out %d != "
                               "dim %d" % (i, self.roles[i], r_prev, r_next, mods[i].dim))
            facts.append({"node": self.roles[i], "fact": "im=ker",
                          "rank_in": r_prev, "rank_out": r_next, "dim": mods[i].dim})
        facts.append({"node": self.roles[-1], "fact": "surjective", "rank": rl,
                      "dim": mods[-1].dim})
        return True, facts

    def assert_exact(self):
        ok, info = self.check_exact()
        if not ok:
            raise QuiverlabError("yoneda: constructed sequence failed exactness: " + info)
        return info


# --------------------------------------------------------------------------- #
# The construction.
# --------------------------------------------------------------------------- #
def yoneda_sequence(M, N, cocycle, n, terms=None, dmats=None):
    """Construct the ``n``-fold exact sequence realizing the ``Ext^n_A(M, N)`` class
    whose cocycle is ``cocycle`` (an ``N.dim x P_n.dim`` matrix, the module map
    ``f: P_n -> N``). ``n >= 1``. Pass a precomputed minimal resolution
    ``(terms, dmats)`` (as from :func:`minimal_resolution`) to avoid recomputing it.

    Returns a :class:`YonedaSequence`. The construction is EXACT by derivation; call
    ``.assert_exact()`` (or ``.check_exact()``) to self-certify a given instance."""
    if n < 1:
        raise QuiverlabError("yoneda_sequence: degree must be >= 1 (Ext^0 = Hom(M,N) is "
                             "not realized by an extension)")
    A = M.algebra
    dom = A.domain
    if terms is None:
        terms, dmats = minimal_resolution(M, n + 1)
    if n >= len(terms) or terms[n].module is None or terms[n].dim == 0:
        raise QuiverlabError(
            "yoneda_sequence: the minimal resolution of M terminates before degree %d, "
            "so Ext^%d(M, N) = 0 (no class to realize)" % (n, n))
    P_prev = terms[n - 1].module           # P_{n-1}
    d_n = dmats[n]                          # P_n -> P_{n-1}
    if not (d_n and d_n[0]):
        raise QuiverlabError("yoneda_sequence: d_%d is empty, Ext^%d(M, N) = 0" % (n, n))
    if len(cocycle) != N.dim or len(cocycle[0]) != len(d_n[0]):
        raise QuiverlabError(
            "yoneda_sequence: cocycle shape %dx%d does not match f: P_%d -> N "
            "(expected %dx%d)" % (len(cocycle), len(cocycle[0]) if cocycle else 0,
                                  n, N.dim, len(d_n[0])))

    # Omega^n M = im(d_n) subset P_{n-1}: pivot columns of d_n give a basis, and their
    # column indices are standard-basis preimages in P_n (see the module docstring).
    piv = lm.column_space_pivots(d_n, dom)
    B = lm.cols_to_matrix([lm.col(d_n, p) for p in piv])            # inclusion i (P_prev.dim x r)
    phi = lm.cols_to_matrix([lm.col(cocycle, p) for p in piv])      # phi: Omega -> N (N.dim x r)
    r = len(piv)

    PN = _direct_sum2(P_prev, N)
    # psi(x) = (i(x), -phi(x)); its columns (independent, since i is injective) span im(psi).
    psi_cols = []
    for j in range(r):
        top = lm.col(B, j)
        bot = [dom.neg(x) for x in lm.col(phi, j)]
        psi_cols.append(list(top) + bot)
    # im(psi) is a SUBMODULE iff psi is a module map iff `cocycle` is a genuine cocycle.
    # Verify it (the quotient action would silently be ill-defined otherwise, since the
    # ambient basis spans everything and the coordinate solve never fails on its own).
    _assert_submodule(PN, psi_cols, dom, n)
    Q, proj, lift = _quotient_with_maps(PN, psi_cols, dom,
                                        name=("E" if n == 1 else "Q"))

    # iota: N -> Q  = proj . (inclusion of N as the second summand of PN)
    iota = lm.cols_to_matrix([lm.col(proj, P_prev.dim + k) for k in range(N.dim)]) \
        if Q.dim else lm.zeros(0, N.dim, dom)

    modules = [N, Q]
    maps = [iota]
    roles = ["sub", "middle"]
    term_index = [None, None]

    if n == 1:
        # end map eps = d_0 : P_0 -> M, extended by 0 on N, descended through the lift.
        eps = dmats[0]
        outmap = _augment_zero(eps, N.dim, dom)            # [eps | 0]
        pi = lm.matmul(outmap, lift, dom) if lift and lift[0] else lm.zeros(M.dim, 0, dom)
        modules.append(M)
        maps.append(pi)
        roles.append("quotient")
        term_index.append(None)
    else:
        # del: Q -> P_{n-2} = [d_{n-1} | 0] . lift
        P_below = terms[n - 2].module
        d_below = dmats[n - 1]                              # P_{n-1} -> P_{n-2}
        outmap = _augment_zero(d_below, N.dim, dom)
        delta = lm.matmul(outmap, lift, dom) if lift and lift[0] else \
            lm.zeros(P_below.dim, 0, dom)
        modules.append(P_below)
        maps.append(delta)
        roles.append("resolution_term")
        term_index.append(n - 2)
        # the original resolution tail: d_{n-2}, ..., d_1, then eps = d_0.
        for k in range(n - 2, 0, -1):
            modules.append(terms[k - 1].module)
            maps.append(dmats[k])
            # every P_{k-1} in the resolution tail (down to and including P_0) is a
            # resolution term; the only non-term in the tail is the end object M, which
            # is appended below with role "quotient".
            roles.append("resolution_term")
            term_index.append(k - 1)
        modules.append(M)
        maps.append(dmats[0])
        roles.append("quotient")
        term_index.append(None)

    return YonedaSequence(n, modules, maps, roles, term_index, dom)


def _assert_submodule(ambient, sub_cols, dom, n):
    """Refuse loudly unless ``span(sub_cols)`` is ``A``-stable in ``ambient`` -- i.e.
    unless the pushout map ``psi`` is a module map, i.e. unless ``cocycle`` is a genuine
    cocycle. (The quotient's own coordinate solve cannot catch this: its basis spans the
    whole ambient space.)"""
    if not sub_cols:
        return
    S = lm.cols_to_matrix(sub_cols)
    for b in _generators(ambient.algebra):
        img = lm.matmul(ambient.action[b], S, dom)
        if lm.solve_columns(S, img, dom) is None:
            raise QuiverlabError(
                "yoneda_sequence: the induced map Omega^%d M -> N is not an A-module "
                "map, so the given cochain is not a cocycle -- no exact sequence "
                "realizes it" % n)


def _augment_zero(mat, extra_cols, dom):
    """[mat | 0] -- append `extra_cols` zero columns (the map is 0 on the N summand)."""
    rows = len(mat)
    out = []
    for i in range(rows):
        out.append(list(mat[i]) + [dom.zero()] * extra_cols)
    return out


# -- friendly aliases the plan names ---------------------------------------- #
def baer_extension(M, N, cocycle, terms=None, dmats=None):
    """The Baer extension ``0 -> N -> E -> M -> 0`` of an ``Ext^1_A(M, N)`` class
    (``n = 1`` case of :func:`yoneda_sequence`). ``E`` is the pushout module."""
    return yoneda_sequence(M, N, cocycle, 1, terms=terms, dmats=dmats)


def spliced_sequence(M, N, cocycle, n, terms=None, dmats=None):
    """The spliced ``n``-fold exact sequence ``0 -> N -> Q -> P_{n-2} -> ... -> P_0 ->
    M -> 0`` of an ``Ext^n_A(M, N)`` class (``n >= 2``)."""
    if n < 2:
        raise QuiverlabError("spliced_sequence is for n >= 2; use baer_extension for n = 1")
    return yoneda_sequence(M, N, cocycle, n, terms=terms, dmats=dmats)
