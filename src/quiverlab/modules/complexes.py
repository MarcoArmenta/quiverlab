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
        return Y

    def truncate(self, lo, hi):
        """Brutal (stupid) truncation to degrees ``[lo, hi]``: keep those terms and
        the differentials strictly inside the window (``lo < n <= hi``); the boundary
        differential ``d_lo`` (out of the window) is dropped."""
        new_terms = {n: M for n, M in self._terms.items() if lo <= n <= hi}
        new_dmats = {n: mat for n, mat in self._dmats.items() if lo < n <= hi}
        return ChainComplex(new_terms, new_dmats, check=False)

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
