"""Bounded chain complexes over ``mod A`` and their maps (Plan 39 / C8).

A first-class layer over the module surface: validated bounded complexes
(:class:`ChainComplex`), chain maps (:class:`ChainMap`), mapping cones and
triangles, the derived-iso test, the Hom total complex (hyper-Hom) and its
general form via a certified projective model (hyper-Ext). P42 (spectral
sequences) and P43 (derived category) build on exactly the surface declared here.

Convention -- pinned **homological** (verbatim, do not silently reindex):

    d_n : C_n -> C_{n-1},   matrices  rows = target C_{n-1},  cols = source C_n

which is byte-for-byte the ``dmats`` layout of
``modules.resolution.minimal_resolution`` (rows=target). Cohomological indexing
is presentation-only: ``C^n := C_{-n}``. Every cohomological consumer
(injectives, Ext, the Hom total complex) re-indexes at ITS OWN boundary -- this
file is homological-only and never does that reindex for the caller.

Homology dims come from the rank formula ``dim H_n = dim ker d_n - rank d_{n+1}``;
homology *modules* from ``radtopsoc.submodule`` (the cycles) then ``quotient``
(by the boundaries expressed in the cycle basis), the tested
homology-of-two-maps pattern of ``complex_reps._reps_from_complex``.

All validation is loud (:class:`~quiverlab.errors.QuiverlabError`); ``check=False``
fast paths exist only for internally-constructed data (cones, models) whose
certificates (cone acyclicity, ``d^2 = 0``) are asserted separately. Float-free:
every scalar goes through the Domain (the AST gate scans this file).
"""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module, _coerce_matrix


class ChainComplex:
    """A validated bounded chain complex of ``A``-modules on one fixed side.

    ``terms``: ``degree -> Module`` (bounded; a missing degree is the zero module).
    ``dmats[n]`` is ``d_n : terms[n] -> terms[n-1]`` with rows=target (the
    ``minimal_resolution`` layout). With ``check=True`` every ``d_n`` is validated
    as an ``A``-module map (the P37 :class:`ModuleHom` predicate) and every composite
    ``d_n . d_{n+1}`` is asserted zero -- both refuse loudly."""

    def __init__(self, terms, dmats, check=True):
        raw_terms = {int(n): M for n, M in terms.items() if M is not None}
        reals = list(raw_terms.values())
        if not reals:
            raise QuiverlabError(
                "ChainComplex: at least one nonzero term is required to fix the "
                "algebra and side (the empty complex has no home)")
        A = reals[0].algebra
        side = reals[0].side
        for M in reals:
            if M.algebra is not A:
                raise QuiverlabError(
                    "ChainComplex: all terms must be modules over the same algebra")
            if M.side != side:
                raise QuiverlabError(
                    f"ChainComplex: mixing a {side} term and a {M.side} term "
                    "(a complex lives on one fixed side)")
        self.algebra = A
        self.side = side
        self.domain = A.domain
        self._terms = raw_terms
        self._ref_labels = list(reals[0].action.keys())
        dom = self.domain
        self._dmats = {int(n): _coerce_matrix(mat, dom) for n, mat in dmats.items()}
        self.__zero = None
        if check:
            self._validate()

    # -- validation ---------------------------------------------------------- #
    def _validate(self):
        from quiverlab.modules.morphism import ModuleHom
        dom = self.domain
        for n, mat in self._dmats.items():
            src, tgt = self.term(n), self.term(n - 1)
            # ModuleHom validates the shape AND the intertwining (loud "module map")
            ModuleHom(src, tgt, mat, check=True)
        for n in sorted(self._dmats):
            d_n = self._dmats.get(n)
            d_n1 = self._dmats.get(n + 1)
            if d_n is None or d_n1 is None:
                continue
            comp = lm.matmul(d_n, d_n1, dom)
            if comp and any(not dom.is_zero(x) for row in comp for x in row):
                raise QuiverlabError(
                    f"d.d != 0: the composite d_{n} ∘ d_{n + 1} is nonzero, so this "
                    "is not a chain complex")

    # -- structure ----------------------------------------------------------- #
    def _zero_module(self, name="0"):
        dom = self.domain
        return Module(self.algebra, 0,
                      {lab: lm.zeros(0, 0, dom) for lab in self._ref_labels},
                      name=name, side=self.side)

    def degrees(self):
        return sorted(self._terms)

    def term(self, n):
        M = self._terms.get(int(n))
        if M is not None:
            return M
        if self.__zero is None:
            self.__zero = self._zero_module()
        return self.__zero

    def differential(self, n):
        from quiverlab.modules.morphism import ModuleHom
        src, tgt = self.term(n), self.term(n - 1)
        mat = self._dmats.get(int(n))
        if mat is None:
            mat = lm.zeros(tgt.dim, src.dim, self.domain)
        return ModuleHom(src, tgt, mat, check=False)

    def _dmat(self, n):
        return self._dmats.get(int(n))

    def total_dim(self):
        return sum(self.term(n).dim for n in self.degrees())

    # -- shift / truncate ---------------------------------------------------- #
    def shift(self, k):
        """``X[k]_n = X_{n-k}`` with ``d -> (-1)^k d`` (the domain-exact sign flip).
        Homology is carried up by ``k``; an odd shift flips every differential's sign,
        preserving ``d.d = 0`` (both differentials flip, their product is unchanged)."""
        dom = self.domain
        k = int(k)
        new_terms = {n + k: M for n, M in self._terms.items()}
        if k % 2 == 0:
            new_dmats = {n + k: mat for n, mat in self._dmats.items()}
        else:
            neg1 = dom.neg(dom.one())
            new_dmats = {n + k: [[dom.mul(neg1, x) for x in row] for row in mat]
                         for n, mat in self._dmats.items()}
        Y = ChainComplex(new_terms, new_dmats, check=False)
        if getattr(self, "_perfect", False):
            Y._perfect = True
        prov = getattr(self, "_proj_vertices", None)
        if prov is not None:                    # X[k]_{n+k} = X_n: carry the vertex lists
            Y._proj_vertices = {n + k: list(vs) for n, vs in prov.items()}
        tagged = getattr(self, "_term_provenance", None)
        if tagged is not None:
            Y._term_provenance = (tagged[0],
                                  {n + k: list(vs) for n, vs in tagged[1].items()})
        return Y

    def truncate(self, lo, hi):
        """Brutal (stupid) truncation to degrees ``[lo, hi]``: keep those terms and
        the differentials strictly inside the window (``lo < n <= hi``); the boundary
        differential ``d_lo`` (out of the window) is dropped."""
        new_terms = {n: M for n, M in self._terms.items() if lo <= n <= hi}
        if not any(M.dim for M in new_terms.values()):
            raise QuiverlabError(
                f"truncate: the window [{lo}, {hi}] misses the support "
                f"{self.degrees()} entirely -- an empty truncation has no home "
                "(devil's-advocate targeted message, 2026-08-05).")
        new_dmats = {n: mat for n, mat in self._dmats.items() if lo < n <= hi}
        out = ChainComplex(new_terms, new_dmats, check=False)
        # brutal truncation of a perfect complex stays perfect
        out._perfect = self._perfect
        # provenance survives, restricted to the window (devil's-advocate LOW)
        prov = getattr(self, "_proj_vertices", None)
        if prov is not None:
            out._proj_vertices = {n: list(vs) for n, vs in prov.items()
                                  if lo <= n <= hi}
        tagged = getattr(self, "_term_provenance", None)
        if tagged is not None:
            out._term_provenance = (tagged[0], {n: list(vs)
                                    for n, vs in tagged[1].items()
                                    if lo <= n <= hi})
        return out

    # -- homology ------------------------------------------------------------ #
    def homology_dims(self, lo=None, hi=None):
        """``{n: dim H_n}`` over ``[lo, hi]`` (default: the support range).
        ``dim H_n = dim ker d_n - rank d_{n+1}`` (absent differentials are zero maps)."""
        dom = self.domain
        degs = self.degrees()
        if not degs:
            return {}
        if lo is None:
            lo = degs[0]
        if hi is None:
            hi = degs[-1]
        out = {}
        for n in range(lo, hi + 1):
            cn = self.term(n).dim
            d_n = self._dmats.get(n)
            d_n1 = self._dmats.get(n + 1)
            r_n = lm.mat_rank(d_n, dom) if (d_n and d_n[0]) else 0
            r_n1 = lm.mat_rank(d_n1, dom) if (d_n1 and d_n1[0]) else 0
            out[n] = (cn - r_n) - r_n1
        return out

    def homology(self, n):
        """``H_n = Z_n / B_n`` as a Module: cycles ``Z_n = ker d_n`` (the whole term
        when ``d_n`` is absent), boundaries ``B_n = im d_{n+1}`` expressed in the cycle
        basis then quotiented out (the ``complex_reps._reps_from_complex`` pattern,
        promoted to genuine sub/quotient modules)."""
        from quiverlab.modules.radtopsoc import submodule, quotient
        dom = self.domain
        n = int(n)
        Cn = self.term(n)
        if Cn.dim == 0:
            return self._zero_module(name=f"H_{n}")
        d_n = self._dmats.get(n)
        if d_n is not None and d_n and d_n[0] and self.term(n - 1).dim:
            cyc = lm.kernel_columns(d_n, dom)
        else:
            ident = lm.identity(Cn.dim, dom)
            cyc = [lm.col(ident, j) for j in range(Cn.dim)]
        if not cyc:
            return self._zero_module(name=f"H_{n}")
        Z = submodule(Cn, cyc, name=f"Z_{n}")
        d_n1 = self._dmats.get(n + 1)
        bdy_in_Z = []
        if d_n1 is not None and d_n1 and d_n1[0]:
            bdyC = [lm.col(d_n1, j) for j in range(len(d_n1[0]))]
            B = lm.cols_to_matrix(cyc)
            coords = lm.solve_columns(B, lm.cols_to_matrix(bdyC), dom)
            if coords is None:
                raise QuiverlabError(
                    f"homology({n}): a boundary column is not a cycle -- "
                    "im d_{n+1} escapes ker d_n (d.d != 0)")
            bdy_in_Z = coords
        return quotient(Z, bdy_in_Z, name=f"H_{n}")

    def is_acyclic(self):
        return all(v == 0 for v in self.homology_dims().values())

    # -- provenance / perfectness ------------------------------------------- #
    def is_perfect(self):
        """Every nonzero term is projective. Fast path: the ``from_projective_resolution``
        / cone / shift provenance flag. Otherwise the EXACT projective-cover-rank
        certificate: ``Q`` is projective iff its projective cover ``P(Q) ->> Q`` (built
        by ``projective_cover``, exact over any Domain) has ``dim P(Q) = dim Q`` (a
        surjection between equal finite dimensions is an iso). This decides in every
        case; it raises only when the algebra has no quiver presentation (so there is
        no projective to compare against)."""
        if getattr(self, "_perfect", False):
            return True
        from quiverlab.modules.resolution import projective_cover
        for n in self.degrees():
            Q = self.term(n)
            if Q.dim == 0:
                continue
            try:
                Q0, _, _ = projective_cover(Q)
            except QuiverlabError as exc:
                raise QuiverlabError(
                    f"is_perfect: cannot certify projectivity of the degree-{n} term: "
                    f"{exc}")
            if Q0.dim != Q.dim:
                return False
        return True

    # -- constructors -------------------------------------------------------- #
    @classmethod
    def stalk(cls, M, degree=0):
        """The one-term complex with ``M`` in ``degree`` (a module viewed as a complex)."""
        return cls({int(degree): M}, {}, check=True)

    @classmethod
    def from_projective_resolution(cls, M, length):
        """The perfect complex ``Q_length -> ... -> Q_0`` (degrees ``length..0``) of the
        minimal projective resolution of ``M``; its homology is ``M`` in degree 0 (by
        exactness). The augmentation ``d_0: Q_0 -> M`` is EXCLUDED (the complex is
        ``Q_*``; ``H_0 = coker d_1 = M``)."""
        from quiverlab.modules.resolution import minimal_resolution
        terms_list, dmats_list = minimal_resolution(M, length)
        terms = {}
        for n, td in enumerate(terms_list):
            if td.module is not None and td.dim > 0:
                terms[n] = td.module
        dmats = {}
        for n in range(1, len(dmats_list)):
            if n in terms and (n - 1) in terms:
                mat = dmats_list[n]
                if mat and mat[0]:
                    dmats[n] = mat
        X = cls(terms, dmats, check=True)
        X._perfect = True
        # Provenance for the derived AR translate (Plan 43 / Task 2): the per-degree
        # projective-summand vertex multiset, in the resolution's summand order (the
        # block order of every dmat). nu is applied termwise off exactly this list.
        X._proj_vertices = {n: list(terms_list[n].vertices)
                            for n in terms if terms_list[n].module is not None}
        return X

    def __repr__(self):
        degs = self.degrees()
        if not degs:
            return "ChainComplex(0)"
        body = " <- ".join(f"C_{n}({self.term(n).dim})" for n in reversed(degs))
        return f"ChainComplex[{self.side} {self.algebra}]: {body}"


# --------------------------------------------------------------------------- #
# Domain-generic helpers.
# --------------------------------------------------------------------------- #
def _get(mat, i, j, dom):
    """Entry ``(i, j)`` of a possibly-shapeless matrix (matmul returns ``[]`` when an
    inner dimension is 0); out-of-range / empty reads as the domain zero."""
    if not mat or i >= len(mat) or not mat[0] or j >= len(mat[0]):
        return dom.zero()
    return mat[i][j]


def _squares_equal(lhs, rhs, rows, cols, dom):
    """Exact equality of two ``rows x cols`` matrices, tolerant of the shapeless ``[]``
    that ``matmul`` returns through a zero-dimensional module."""
    for i in range(rows):
        for j in range(cols):
            if not dom.is_zero(dom.sub(_get(lhs, i, j, dom), _get(rhs, i, j, dom))):
                return False
    return True


def _block_diag_module(Xm, Ym, dom, name="(+)"):
    """``Xm (+) Ym`` (X first) as a block-diagonal Module over the shared algebra.
    Zero-dimensional summands contribute nothing (no action read), so a stalk cone with
    one empty block is exact."""
    dx, dy = Xm.dim, Ym.dim
    n = dx + dy
    A = Xm.algebra if dx else Ym.algebra
    side = Xm.side if dx else Ym.side
    if dx and dy:
        labels = [lab for lab in Xm.action if lab in Ym.action]
    elif dx:
        labels = list(Xm.action)
    else:
        labels = list(Ym.action)
    action = {}
    for lab in labels:
        blk = lm.zeros(n, n, dom)
        if dx:
            Xb = Xm.action[lab]
            for i in range(dx):
                bi, xi = blk[i], Xb[i]
                for j in range(dx):
                    bi[j] = xi[j]
        if dy:
            Yb = Ym.action[lab]
            for i in range(dy):
                bi, yi = blk[dx + i], Yb[i]
                for j in range(dy):
                    bi[dx + j] = yi[j]
        action[lab] = blk
    return Module(A, n, action, name=name, side=side)


class ChainMap:
    """A validated chain map ``f: src -> tgt`` between complexes on one side.

    ``components[n]`` is the module map ``f_n : src_n -> tgt_n`` (rows=target). With
    ``check=True`` every ``f_n`` is validated as an ``A``-module map AND every square
    ``d^tgt_n . f_n == f_{n-1} . d^src_n`` is asserted -- both refuse loudly ("chain
    map ...")."""

    def __init__(self, src, tgt, components, check=True):
        if src.algebra is not tgt.algebra:
            raise QuiverlabError(
                "ChainMap: src and tgt are complexes over different algebras")
        if src.side != tgt.side:
            raise QuiverlabError(
                f"ChainMap: src is {src.side}, tgt is {tgt.side} (a chain map stays "
                "on one side)")
        self.src, self.tgt = src, tgt
        self.domain = src.domain
        self._comps = {int(n): _coerce_matrix(mat, self.domain)
                       for n, mat in components.items()}
        if check:
            self._validate()

    def component(self, n):
        src_n, tgt_n = self.src.term(n), self.tgt.term(n)
        mat = self._comps.get(int(n))
        if mat is None:
            return lm.zeros(tgt_n.dim, src_n.dim, self.domain)
        return mat

    def _validate(self):
        from quiverlab.modules.morphism import ModuleHom
        dom = self.domain
        for n, mat in self._comps.items():
            try:
                ModuleHom(self.src.term(n), self.tgt.term(n), mat, check=True)
            except QuiverlabError as exc:
                raise QuiverlabError(
                    f"chain map: component f_{n} is not a module map ({exc})")
        degs = set(self.src.degrees()) | set(self.tgt.degrees())
        degs |= {n - 1 for n in list(degs)} | {n + 1 for n in list(degs)}
        for n in degs:
            fn = self.component(n)
            fn1 = self.component(n - 1)
            dX = self.src._dmats.get(n)
            dY = self.tgt._dmats.get(n)
            if dX is None:
                dX = lm.zeros(self.src.term(n - 1).dim, self.src.term(n).dim, dom)
            if dY is None:
                dY = lm.zeros(self.tgt.term(n - 1).dim, self.tgt.term(n).dim, dom)
            lhs = lm.matmul(dY, fn, dom)          # X_n -> Y_{n-1}
            rhs = lm.matmul(fn1, dX, dom)         # X_n -> Y_{n-1}
            rows, cols = self.tgt.term(n - 1).dim, self.src.term(n).dim
            if not _squares_equal(lhs, rhs, rows, cols, dom):
                raise QuiverlabError(
                    f"chain map square at degree {n} does not commute: "
                    "d^tgt . f != f . d^src")

    # -- mapping cone / triangle / quasi-iso -------------------------------- #
    def _cone_degrees(self):
        ndegs = set()
        for p in self.src.degrees():
            ndegs.add(p + 1)                      # X_{n-1} nonzero -> n = p+1
        for q in self.tgt.degrees():
            ndegs.add(q)
        return ndegs

    def _cone_diff(self, n):
        """``d^cone_n : X_{n-1} (+) Y_n -> X_{n-2} (+) Y_{n-1}``, the block matrix
        ``[[-d_X, 0], [f, d_Y]]`` (rows=target)."""
        X, Y, dom = self.src, self.tgt, self.domain
        dXm1, dYn = X.term(n - 1).dim, Y.term(n).dim
        dXm2, dYm1 = X.term(n - 2).dim, Y.term(n - 1).dim
        M = lm.zeros(dXm2 + dYm1, dXm1 + dYn, dom)
        dX = X._dmats.get(n - 1)                   # X_{n-1} -> X_{n-2}
        if dX and dX[0]:
            for i in range(dXm2):
                mi, di = M[i], dX[i]
                for j in range(dXm1):
                    mi[j] = dom.neg(di[j])
        f = self.component(n - 1)                  # X_{n-1} -> Y_{n-1}
        for i in range(dYm1):
            mi = M[dXm2 + i]
            for j in range(dXm1):
                mi[j] = _get(f, i, j, dom)
        dY = Y._dmats.get(n)                       # Y_n -> Y_{n-1}
        if dY and dY[0]:
            for i in range(dYm1):
                mi, di = M[dXm2 + i], dY[i]
                for j in range(dYn):
                    mi[dXm1 + j] = di[j]
        return M

    def cone(self):
        """The mapping cone ``cone(f)_n = X_{n-1} (+) Y_n`` with differential
        ``[[-d_X, 0], [f, d_Y]]``. ``d^2 = 0`` is automatic from ``f`` being a chain map
        and is re-asserted here (``check=True``) as the self-certificate."""
        X, Y, dom = self.src, self.tgt, self.domain
        terms = {}
        for n in self._cone_degrees():
            Xm, Ym = X.term(n - 1), Y.term(n)
            if Xm.dim + Ym.dim == 0:
                continue
            terms[n] = _block_diag_module(Xm, Ym, dom, name=f"cone_{n}")
        dmats = {}
        for n in list(terms):
            if (n - 1) in terms:
                dmats[n] = self._cone_diff(n)
        C = ChainComplex(terms, dmats, check=True)
        if getattr(X, "_perfect", False) and getattr(Y, "_perfect", False):
            C._perfect = True
        return C

    def then(self, g):
        """``self`` then ``g`` (left-to-right, mirroring ``ModuleHom.then``):
        ``(self.then(g)).component(n) = g.component(n) @ self.component(n)``. Refuses
        if the middle complexes differ (``g.src is self.tgt``). ``check=False`` -- each
        square is the matmul of two commuting squares; callers re-validate on demand."""
        if g.src is not self.tgt:
            raise QuiverlabError("ChainMap.then: middle complexes differ")
        dom = self.domain
        comps = {}
        degs = set(self.src.degrees()) | set(g.tgt.degrees())
        for n in degs:
            fn, gn = self.component(n), g.component(n)
            if self.tgt.term(n).dim == 0:            # matmul is shapeless through 0
                comps[n] = lm.zeros(g.tgt.term(n).dim, self.src.term(n).dim, dom)
            else:
                comps[n] = lm.matmul(gn, fn, dom)
        return ChainMap(self.src, g.tgt, comps, check=False)

    def is_quasi_iso(self):
        """``f`` is a quasi-isomorphism iff its mapping cone is acyclic."""
        return self.cone().is_acyclic()

    def triangle(self):
        """The distinguished triangle ``X -f-> Y -i-> cone(f) -p-> X[1]``: returns
        ``(X, Y, cone, (i, p))`` with ``i`` the block inclusion ``Y -> cone`` and ``p``
        the block projection ``cone -> X[1]``."""
        X, Y, dom = self.src, self.tgt, self.domain
        C = self.cone()
        # i: Y_n -> cone_n = X_{n-1} (+) Y_n, the second-block inclusion [0; I].
        incl = {}
        for n in Y.degrees():
            dXm1, dYn = X.term(n - 1).dim, Y.term(n).dim
            if dYn == 0:
                continue
            mat = lm.zeros(dXm1 + dYn, dYn, dom)
            for j in range(dYn):
                mat[dXm1 + j][j] = dom.one()
            incl[n] = mat
        i_map = ChainMap(Y, C, incl, check=True)
        # p: cone_n = X_{n-1} (+) Y_n -> X[1]_n = X_{n-1}, the first-block projection.
        Xshift = X.shift(1)
        proj = {}
        for n in C.degrees():
            dXm1, dYn = X.term(n - 1).dim, Y.term(n).dim
            if dXm1 == 0:
                continue
            mat = lm.zeros(dXm1, dXm1 + dYn, dom)
            for i in range(dXm1):
                mat[i][i] = dom.one()
            proj[n] = mat
        p_map = ChainMap(C, Xshift, proj, check=True)
        return (X, Y, C, (i_map, p_map))


def identity_chain_map(X):
    """The identity chain map ``id_X : X -> X``."""
    dom = X.domain
    comps = {n: lm.identity(X.term(n).dim, dom) for n in X.degrees()}
    return ChainMap(X, X, comps, check=True)


# --------------------------------------------------------------------------- #
# Hom total complex (hyper-Hom) for a perfect source.
#
# SIGN / INDEX CONVENTION (Weibel, "An Introduction to Homological Algebra",
# 2.7.4) -- documented verbatim:
#
#   Hom^n(X, Y) = (+)_p Hom_A(X_p, Y_{p-n})
#   (delta f)_p = d^Y_{p-n} . f_p  -  (-1)^n . f_{p-1} . d^X_p              (*)
#
# delta: Hom^n -> Hom^{n+1}. The two blocks of (*) land, for a basis element
# phi in Hom(X_{p0}, Y_{p0-n}), in:
#   * block p0   of Hom^{n+1}: value  d^Y_{p0-n} . phi   (into Y_{p0-n-1});
#   * block p0+1 of Hom^{n+1}: value  -(-1)^n . phi . d^X_{p0+1}.
#
# HONEST NOTE ON THE SIGN (verified, Plan-39 implementation session): with the
# alternating coefficient c_n = eps.(-1)^n on the second block, (delta^2 f)_p
# collects d^Y.d^Y.f_p (=0), f_{p-2}.d^X.d^X (=0), and a single cross term with
# coefficient c_n + c_{n+1} = eps.(-1)^n + eps.(-1)^{n+1} = 0. This vanishes for
# BOTH eps = -1 (Weibel) and eps = +1: the two choices give ISOMORPHIC cochain
# complexes (a degreewise +/-1 rescaling conjugates one to the other), hence the
# SAME homology dimensions. So the arbiter `test_resolution_hyper_hom_computes_ext`
# and the delta.delta = 0 self-check below are BOTH sign-independent -- they pin
# the block indexing and module-map placement, not the sign. The sign is a genuine
# convention; we take Weibel's eps = -1. (Were the sign ever to matter -- e.g. for
# a signed representative, not a dimension -- flip eps here once, never per degree.)
# For a PERFECT X, H^n(Hom^.) computes Hom_{D^b(mod A)}(X, Y[n]) (hyper-Ext); with
# Y = stalk(N) and X = a projective resolution of M this is Ext^n_A(M, N).
# --------------------------------------------------------------------------- #
def _flatten(mat):
    return [x for row in mat for x in row]


def _hom_total_blocks(X, Y, n, dom):
    """Ordered block basis of ``Hom^n(X, Y) = (+)_p Hom_A(X_p, Y_{p-n})``.
    Returns ``(blocks, total_dim)``; each block is
    ``{p, q, homs, offset, count}`` with ``homs`` the ``hom_space`` basis
    (``Y_q.dim x X_p.dim`` matrices) and ``q = p - n``."""
    from quiverlab.modules.hom import hom_space
    blocks, off = [], 0
    for p in X.degrees():
        Xp = X.term(p)
        q = p - n
        Yq = Y.term(q)
        if Xp.dim == 0 or Yq.dim == 0:
            continue
        H = hom_space(Xp, Yq)
        if not H:
            continue
        blocks.append({"p": p, "q": q, "homs": H, "offset": off, "count": len(H)})
        off += len(H)
    return blocks, off


def _place_hom(colvec, tgt_by_p, p, A, dom):
    """Write the coordinates of the module map ``A`` into the ``p``-block of a
    Hom^{n+1} column, solving against that block's ``hom_space`` basis."""
    if not A or all(dom.is_zero(x) for row in A for x in row):
        return
    entry = tgt_by_p.get(p)
    if entry is None:
        raise QuiverlabError(
            "hyper-Hom: a differential image landed outside the Hom total space "
            "(block indexing bug)")
    b, flat = entry
    coords = lm.solve_columns(flat, lm.cols_to_matrix([_flatten(A)]), dom)
    if coords is None:
        raise QuiverlabError(
            "hyper-Hom: a differential image is not expressible in the target Hom "
            "block basis (not a module map)")
    base = b["offset"]
    for i, c in enumerate(coords[0]):
        colvec[base + i] = dom.add(colvec[base + i], c)


def _delta_total(X, Y, n, dom):
    """Matrix of ``delta^n : Hom^n -> Hom^{n+1}`` in the ordered block bases
    (rows = ``dim Hom^{n+1}``, cols = ``dim Hom^n``), by convention (*) above.
    Returns ``(matrix, dim Hom^n, dim Hom^{n+1})``."""
    src_blocks, src_dim = _hom_total_blocks(X, Y, n, dom)
    tgt_blocks, tgt_dim = _hom_total_blocks(X, Y, n + 1, dom)
    tgt_by_p = {}
    for b in tgt_blocks:
        flat = lm.cols_to_matrix([_flatten(h) for h in b["homs"]])
        tgt_by_p[b["p"]] = (b, flat)
    neg_sign = dom.neg(dom.one()) if n % 2 == 0 else dom.one()   # -(-1)^n
    cols = []
    for b in src_blocks:
        p0, q0 = b["p"], b["q"]
        dY = Y._dmats.get(q0)                     # Y_{q0} -> Y_{q0-1}
        dX = X._dmats.get(p0 + 1)                 # X_{p0+1} -> X_{p0}
        for phi in b["homs"]:
            colvec = [dom.zero()] * tgt_dim
            if dY and dY[0]:
                _place_hom(colvec, tgt_by_p, p0, lm.matmul(dY, phi, dom), dom)
            if dX and dX[0]:
                A2 = lm.matmul(phi, dX, dom)
                A2 = [[dom.mul(neg_sign, x) for x in row] for row in A2]
                _place_hom(colvec, tgt_by_p, p0 + 1, A2, dom)
            cols.append(colvec)
    mat = lm.cols_to_matrix(cols) if cols else lm.zeros(tgt_dim, 0, dom)
    return mat, src_dim, tgt_dim


def hyper_hom_dims(X, Y, lo, hi):
    """``{n: dim H^n(Hom^.(X, Y))}`` for ``n`` in ``lo..hi``, using convention (*).

    For ``X`` PERFECT (a bounded complex of projectives) this is
    ``Hom_{D^b(mod A)}(X, Y[n])`` (hyper-Ext); with ``Y = stalk(N)`` and ``X`` a
    projective resolution of ``M`` it is ``Ext^n_A(M, N)`` (the pinned arbiter).
    Raises loudly if ``X`` is not certified perfect (``hyper_ext_dims`` lifts this
    to any bounded source via a projective model)."""
    if not X.is_perfect():
        raise QuiverlabError(
            "hyper_hom_dims: X must be a perfect complex (a bounded complex of "
            "projectives). Use hyper_ext_dims(X, Y, ...) for a general source, "
            "which resolves X to a certified projective model first.")
    dom = X.domain
    deltas = {m: _delta_total(X, Y, m, dom) for m in range(lo - 1, hi + 1)}
    # self-certificate: delta^m . delta^{m-1} == 0. Certifies the block indexing /
    # module-map placement (a mis-indexed block would leave it nonzero); it is
    # sign-independent (see the module header note), so it never flags the sign.
    for m in range(lo, hi + 1):
        dm, dm1 = deltas[m][0], deltas[m - 1][0]
        if dm and dm[0] and dm1 and dm1[0]:
            comp = lm.matmul(dm, dm1, dom)
            if comp and any(not dom.is_zero(x) for row in comp for x in row):
                raise QuiverlabError(
                    "hyper-Hom: delta.delta != 0 -- the Hom total-complex construction "
                    "is inconsistent (block indexing; see the module header)")
    out = {}
    for n in range(lo, hi + 1):
        cn = deltas[n][1]                          # dim Hom^n
        rn = lm.mat_rank(deltas[n][0], dom) if (deltas[n][0] and deltas[n][0][0]) else 0
        rprev = (lm.mat_rank(deltas[n - 1][0], dom)
                 if (deltas[n - 1][0] and deltas[n - 1][0][0]) else 0)
        out[n] = cn - rn - rprev
    return out


# --------------------------------------------------------------------------- #
# Certified projective model + general hyper-Ext.
#
# The model of a bounded complex X is built by iterated MAPPING CONES over the
# stupid filtration (equivalent to the twisted totalization of the termwise
# resolutions, but organized so every step is one Task-2 cone plus a single
# chain-map lift -- no bare higher-correction bookkeeping):
#
#   tau_{<=j} X  =  cone( g_j : stalk(X_j)[j-1] --d^X_j--> tau_{<=j-1} X )
#
# (the connecting map of  0 -> tau_{<=j-1}X -> tau_{<=j}X -> stalk(X_j)[j] -> 0).
# Replace both ends by their perfect models, LIFT g_j through the augmentations
# to a strict-square chain map of models, take the cone: the cone of perfect
# complexes is perfect, and the cone of quasi-isos with a STRICT commuting square
# is a quasi-iso (5-lemma / the LES of the triangle). The augmentation of every
# model built this way is degreewise SURJECTIVE (covers are surjective; a cone of
# surjections is surjective), which is exactly what makes the strict lift exist.
#
# CERTIFICATION: `projective_model` never returns an uncertified model -- it
# asserts P.is_perfect() AND eps.is_quasi_iso() (cone acyclicity) before return.
# A genuine quasi-iso needs each term's resolution to TERMINATE within `length`
# (finite projective dimension); otherwise the truncated model is not acyclic and
# the certificate refuses loudly (honest -- a truncated resolution is not a model).
# The lift-solve is canonicalized via reduce_mod_nullspace (byte-reproducible, the
# CS correction-solve precedent).
# --------------------------------------------------------------------------- #
def _combine_homs(H, coeffs, sdim, tdim, dom):
    out = lm.zeros(tdim, sdim, dom)
    for c, mat in zip(coeffs, H):
        if dom.is_zero(c):
            continue
        for i in range(tdim):
            oi, mi = out[i], mat[i]
            for j in range(sdim):
                oi[j] = dom.add(oi[j], dom.mul(c, mi[j]))
    return out


def _constraints_ok(psi, constraints, dom):
    for (kind, Mmat, rhs) in constraints:
        got = lm.matmul(Mmat, psi, dom) if kind == "left" else lm.matmul(psi, Mmat, dom)
        rows = len(rhs)
        cols = len(rhs[0]) if (rhs and rhs[0]) else 0
        if not _squares_equal(got, rhs, rows, cols, dom):
            return False
    return True


def _solve_module_map(S, T, constraints, dom):
    """The module map ``psi: S -> T`` (rows=target) meeting every linear constraint
    ``("left", L, rhs)`` (``L . psi = rhs``) / ``("right", R, rhs)`` (``psi . R = rhs``),
    solved over a basis of ``Hom_A(S, T)`` so ``psi`` stays a module map. Returns the
    canonical (``reduce_mod_nullspace``) solution or ``None`` if inconsistent."""
    from quiverlab.fields import linalg as flinalg
    from quiverlab.modules.hom import hom_space
    if S.dim == 0:
        return lm.zeros(T.dim, 0, dom)
    if T.dim == 0:
        return lm.zeros(0, S.dim, dom) if _constraints_ok(
            lm.zeros(0, S.dim, dom), constraints, dom) else None
    H = hom_space(S, T)
    if not H:
        Z = lm.zeros(T.dim, S.dim, dom)
        return Z if _constraints_ok(Z, constraints, dom) else None
    contrib = [[] for _ in H]
    tflat = []
    for (kind, Mmat, rhs) in constraints:
        tflat.extend(_flatten(rhs))
        for k, Hk in enumerate(H):
            c = lm.matmul(Mmat, Hk, dom) if kind == "left" else lm.matmul(Hk, Mmat, dom)
            contrib[k].extend(_flatten(c))
    if not tflat:
        return lm.zeros(T.dim, S.dim, dom)
    Amat = [[contrib[k][r] for k in range(len(H))] for r in range(len(tflat))]
    x = flinalg.solve(Amat, tflat, dom)
    if x is None:
        return None
    x = flinalg.reduce_mod_nullspace(x, Amat, dom)
    return _combine_homs(H, x, S.dim, T.dim, dom)


def _stalk_model(M, degree, length, dom):
    """The perfect model of ``stalk(M)[degree]``: the minimal projective resolution
    of ``M`` shifted so its homology sits in ``degree``, plus the augmentation chain
    map ``eps: model -> stalk(M, degree)`` (the cover at ``degree``, zero elsewhere)."""
    from quiverlab.modules.resolution import minimal_resolution
    terms_list, dmats_list = minimal_resolution(M, length)
    terms, dmats = {}, {}
    for n, td in enumerate(terms_list):
        if td.module is not None and td.dim > 0:
            terms[degree + n] = td.module
    for n in range(1, len(dmats_list)):
        if (degree + n) in terms and (degree + n - 1) in terms:
            mat = dmats_list[n]
            if mat and mat[0]:
                dmats[degree + n] = mat
    model = ChainComplex(terms, dmats, check=True)
    model._perfect = True
    stalk = ChainComplex.stalk(M, degree)
    eps_comps = {degree: dmats_list[0]} if degree in terms else {}
    eps = ChainMap(model, stalk, eps_comps, check=True)
    return model, eps, stalk


def _connecting_map(stalk, Xcur, j, X, dom):
    """``g_j: stalk(X_j)[j-1] -> tau_{<=j-1}X`` with the single component ``d^X_j`` at
    degree ``j-1`` (zero -- a disconnected direct sum -- when ``X_{j-1} = 0``)."""
    dXj = X._dmats.get(j)
    comps = {}
    if dXj is not None and dXj and dXj[0] and Xcur.term(j - 1).dim:
        comps[j - 1] = dXj
    return ChainMap(stalk, Xcur, comps, check=True)


def _lift_cone_map(Smodel, epsS, model, eps, g, dom):
    """Lift ``g . epsS: Smodel -> tau_{<=j-1}X`` through the surjective quasi-iso
    ``eps: model -> tau_{<=j-1}X`` to a strict chain map ``gt: Smodel -> model``
    (``eps . gt == g . epsS``), built degree by degree upward (projectivity of
    ``Smodel_m`` + acyclicity of ``cone(eps)`` guarantee each solve)."""
    Sstalk = g.src
    target = {}
    for m in Smodel.degrees():
        xdim, sdim = g.tgt.term(m).dim, Smodel.term(m).dim
        if Sstalk.term(m).dim == 0:
            target[m] = lm.zeros(xdim, sdim, dom)
        else:
            target[m] = lm.matmul(g.component(m), epsS.component(m), dom)
    gt = {}
    for m in sorted(Smodel.degrees()):
        Sm, Tm = Smodel.term(m), model.term(m)
        constraints = [("left", eps.component(m), target[m])]
        dS = Smodel._dmats.get(m)                  # Smodel_m -> Smodel_{m-1}
        if dS is not None and dS and dS[0] and (m - 1) in gt:
            dModel = model._dmats.get(m)           # model_m -> model_{m-1}
            if dModel is None or not (dModel and dModel[0]):
                dModel = lm.zeros(model.term(m - 1).dim, Tm.dim, dom)
            rhs = lm.matmul(gt[m - 1], dS, dom)
            constraints.append(("left", dModel, rhs))
        psi = _solve_module_map(Sm, Tm, constraints, dom)
        if psi is None:
            raise QuiverlabError(
                "projective_model: could not lift the connecting map through the "
                "surjective quasi-iso at degree %d (the term resolutions are "
                "inconsistent -- report this)" % m)
        gt[m] = psi
    return ChainMap(Smodel, model, gt, check=True)


def _cone_of_eps(g_tilde, g, epsS, eps, model_new, Xcur_new, dom):
    """The induced augmentation ``cone(g_tilde) -> cone(g)``, block-diagonal
    ``diag(epsS_{n-1}, eps_n)`` on ``cone_n = (src)_{n-1} (+) (tgt)_n``. It is a chain
    map because the lift square ``eps . g_tilde = g . epsS`` is strict."""
    comps = {}
    for n in model_new.degrees():
        srcS, srcM = g_tilde.src.term(n - 1).dim, g_tilde.tgt.term(n).dim
        tgtS, tgtM = g.src.term(n - 1).dim, g.tgt.term(n).dim
        eS, eM = epsS.component(n - 1), eps.component(n)
        mat = lm.zeros(tgtS + tgtM, srcS + srcM, dom)
        for i in range(tgtS):
            for jj in range(srcS):
                mat[i][jj] = _get(eS, i, jj, dom)
        for i in range(tgtM):
            for jj in range(srcM):
                mat[tgtS + i][srcS + jj] = _get(eM, i, jj, dom)
        comps[n] = mat
    return ChainMap(model_new, Xcur_new, comps, check=True)


def projective_model(X, length):
    """``(P, eps)`` with ``P`` a perfect complex and ``eps: P -> X`` a chain map,
    CERTIFIED: ``P.is_perfect()`` and ``eps.is_quasi_iso()`` are asserted before
    return (never an uncertified model). Built by iterated cones over the stupid
    filtration (see the section header). ``length`` bounds each term's resolution;
    a genuine quasi-iso requires those resolutions to terminate within it (finite
    pd), else the certificate refuses loudly."""
    support = X.degrees()
    if not support:
        raise QuiverlabError(
            "projective_model: the zero complex has no meaningful projective model")
    dom = X.domain
    p0 = support[0]
    model, eps, Xcur = _stalk_model(X.term(p0), p0, length, dom)
    for j in support[1:]:
        Smodel, epsS, _Sstalk = _stalk_model(X.term(j), j - 1, length, dom)
        g = _connecting_map(_Sstalk, Xcur, j, X, dom)
        g_tilde = _lift_cone_map(Smodel, epsS, model, eps, g, dom)
        model_new = g_tilde.cone()
        Xcur_new = g.cone()
        eps = _cone_of_eps(g_tilde, g, epsS, eps, model_new, Xcur_new, dom)
        model, Xcur = model_new, Xcur_new
    # retarget eps onto the ORIGINAL X (Xcur == X structurally: same terms, diffs)
    eps_to_X = ChainMap(model, X,
                        {n: eps.component(n) for n in model.degrees()}, check=True)
    if not model.is_perfect():
        raise QuiverlabError("projective_model failed certification: the model is "
                             "not perfect")
    if not eps_to_X.is_quasi_iso():
        raise QuiverlabError(
            "projective_model failed certification: the augmentation is not a "
            "quasi-isomorphism -- a term's resolution did not terminate within "
            "length=%d (raise length for a finite-pd complex)" % length)
    model._perfect = True
    return model, eps_to_X


def hyper_ext_dims(X, Y, lo, hi, length=8):
    """``{n: dim Hom_{D^b(mod A)}(X, Y[n])}`` (hyper-Ext) for any bounded source
    ``X``: resolve ``X`` to a certified :func:`projective_model` then apply
    :func:`hyper_hom_dims`. Honest window: a model of length ``L`` certifies degrees
    ``hi <= L - 1`` only (the rank formula at degree ``hi`` reads ``d_{hi+1}``, one
    resolution step of slack); a larger ``hi`` refuses loudly rather than report
    silently-truncated dimensions."""
    if hi > length - 1:
        raise QuiverlabError(
            "hyper_ext_dims: requested top degree hi=%d exceeds the model window "
            "(length=%d certifies hyper-Ext degrees <= length-1=%d); raise length"
            % (hi, length, length - 1))
    P, _eps = projective_model(X, length)
    return hyper_hom_dims(P, Y, lo, hi)
