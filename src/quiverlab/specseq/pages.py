"""SpectralSequence + Page -- the Weibel section-5.4 page engine.

Homological, increasing filtration ``F_p`` (Weibel *An Introduction to
Homological Algebra*, CUP 1994, 5.4.6). Total degree ``n = p + q``; every subspace
is held as a column-span coordinate matrix over ``Tot_n``'s basis, and all algebra
is exact over the Domain (``fields.linalg`` / ``modules.linalg_mod`` / the shared
``_subspace`` helpers).

    Z^r_{p,q} = { x in F_p C_{p+q} : d(x) in F_{p-r} C_{p+q-1} }
    Bdry^r_{p,q} = F_p C_{p+q} ∩ d( F_{p+r-1} C_{p+q+1} )   (= d Z^{r-1}_{p+r-1,q-r+2})
    E^r_{p,q} = Z^r_{p,q} / ( Z^{r-1}_{p-1,q+1} + Bdry^r_{p,q} )
    d^r : E^r_{p,q} -> E^r_{p-r, q+r-1},   induced by d.

ARBITRATION NOTE (Plan-42 implementation, one-time index fix documented here per
the house "flip/correct ONCE" rule): the boundary sub-object of ``E^r`` is
``d Z^{r-1}_{p+r-1,q-r+2} = F_p C_n ∩ d(F_{p+r-1} C_{n+1})`` -- the McCleary/Weibel
5.4.6 form. The plan-brief's ``B^{r-1}_{p,q} = F_p ∩ d(F_{p+r-2} C_{n+1})`` is one
filtration-step short and makes ``E^1_{0,0}`` of the one-step filtration come out
``dim V_0`` instead of ``H_0`` (the boundaries ``im d_1`` vanish); the arbiters
``test_trivial_filtration_collapses_at_E1`` (r=0/r=1 base cases) and the standing
``E_inf totals == total homology`` self-certificate both fail under the short form
and pass under this one. So ``_bdry(p,q,r)`` uses ``F_{p+r-1}``.

Canonical representatives (``Subquotient.reps``): the greedy independent-modulo
filter of ``hochschild/cyclic._capture_generic_reps`` -- columns of ``Z^r`` that
grow the rank over the denominator span, picked in deterministic ``rref`` column
order, so the pick is byte-reproducible (the CS canonicalization mandate applied to
pages). No floats (``src/`` AST gate)."""
from quiverlab.errors import QuiverlabError
from quiverlab.fields import linalg as flinalg
from quiverlab.modules import linalg_mod as lm
from quiverlab.specseq import _subspace as sub


class Subquotient:
    """The module ``E^r_{p,q}`` at one page position: its ``dim`` and its canonical
    coordinate ``reps`` (columns over ``Tot_{p+q}``'s basis)."""

    __slots__ = ("p", "q", "dim", "reps")

    def __init__(self, p, q, reps):
        self.p, self.q = p, q
        self.reps = [list(c) for c in reps]
        self.dim = len(self.reps)

    def __repr__(self):
        return f"E_[{self.p},{self.q}](dim {self.dim})"


class Page:
    """Page ``E^r``: the subquotients at every position plus the ``d^r``
    differentials, and an ``M2 netPage``-style grid."""

    def __init__(self, ss, r, cells):
        self.ss = ss
        self.r = int(r)
        self._cells = cells                      # {(p, q): Subquotient}

    def __getitem__(self, pq):
        p, q = pq
        sq = self._cells.get((int(p), int(q)))
        if sq is None:
            return Subquotient(int(p), int(q), [])
        return sq

    def dim(self, p, q):
        return self[p, q].dim

    @property
    def spots(self):
        return sorted(pq for pq, sq in self._cells.items() if sq.dim > 0)

    def differential(self, p, q):
        """The ``d^r`` matrix ``E^r_{p,q} -> E^r_{p-r, q+r-1}`` (rows = target reps,
        cols = source reps)."""
        return self.ss._dr_matrix(int(p), int(q), self.r)

    def grid(self):
        """``netPage``-style 2-D grid: ``p`` across (columns), ``q`` up (rows), the
        cell ``dim`` in each spot, ``.`` for empty; a fenced monospace block."""
        cells = self._cells
        live = [pq for pq, sq in cells.items() if sq.dim > 0]
        if not live:
            return f"E_{self.r}: (empty)"
        ps = sorted({p for (p, q) in cells})
        qs = sorted({q for (p, q) in cells})
        w = max(2, max(len(str(cells[(p, q)].dim)) for (p, q) in cells
                       if cells[(p, q)].dim > 0) if live else 1)
        lines = [f"E_{self.r} page  (p across, q up):"]
        for q in reversed(qs):
            row = [f"q={q:>2} |"]
            for p in ps:
                sq = cells.get((p, q))
                cell = str(sq.dim) if (sq and sq.dim > 0) else "."
                row.append(cell.rjust(w))
            lines.append(" ".join(row))
        footer = "      " + " ".join(f"p={p}".rjust(w) for p in ps)
        lines.append(footer)
        return "```\n" + "\n".join(lines) + "\n```"

    def __repr__(self):
        return f"Page(E_{self.r}, {len(self.spots)} nonzero spots)"


class SpectralSequence:
    """The spectral sequence of a bounded :class:`FilteredComplex` (Weibel 5.4).

    Pages ``E^r`` (``r >= 0``, memoized) with canonical representatives and induced
    ``d^r`` differentials; converges to the homology of the total complex. The
    standing self-certificate (``E_inf`` totals == total homology) runs at
    construction (attached in Task 4)."""

    def __init__(self, F):
        self.F = F
        self.dom = F.dom
        self._degrees = F.degrees()
        self._levels = F.levels()
        ps = self._levels
        degs = self._degrees
        self._candidates = [(p, n - p) for n in degs for p in ps]
        self.width = len(ps) if ps else 1
        qs = sorted({n - p for n in degs for p in ps}) if (degs and ps) else [0]
        self.height = len(qs) if qs else 1
        self._pages = {}
        from quiverlab.specseq.convergence import certify_convergence
        self.convergence = certify_convergence(self)

    # -- Weibel subspaces ---------------------------------------------------- #
    def _Zr(self, p, q, r):
        """``Z^r_{p,q} = {x in F_p C_n : d(x) in F_{p-r} C_{n-1}}`` as columns over
        ``C_n``."""
        F, dom = self.F, self.dom
        n = p + q
        U = F.piece(n, p)
        if not U:
            return []
        dn = F._dmats.get(n)
        if dn is None or not (dn and dn[0]) or F.dim(n - 1) == 0:
            return U                              # d = 0 => condition vacuous
        Umat = sub.colmat(U)
        A = lm.matmul(dn, Umat, dom)              # d_n |_{F_p C_n}, in C_{n-1} coords
        W = F.piece(n - 1, p - r)                 # F_{p-r} C_{n-1}
        ycols = sub.preimage_selecting(A, W, dom)
        if not ycols:
            return []
        Z = lm.matmul(Umat, sub.colmat(ycols), dom)
        return [lm.col(Z, j) for j in range(len(Z[0]))] if (Z and Z[0]) else []

    def _bdry(self, p, q, r):
        """``Bdry^r_{p,q} = F_p C_n ∩ d(F_{p+r-1} C_{n+1})`` (the McCleary/Weibel
        boundary sub-object; see the module header arbitration note)."""
        F, dom = self.F, self.dom
        n = p + q
        dn1 = F._dmats.get(n + 1)
        if dn1 is None or not (dn1 and dn1[0]):
            return []
        img = sub.image(dn1, F.piece(n + 1, p + r - 1), dom)
        if not img:
            return []
        return sub.intersect(F.piece(n, p), img, dom)

    def _cell(self, p, q, r):
        """``(reps, denom)`` for ``E^r_{p,q}``: canonical representatives (columns of
        ``Z^r`` independent modulo the denominator) and the denominator column set
        ``Z^{r-1}_{p-1,q+1} + Bdry^r_{p,q}`` (both subspaces of ``Z^r``)."""
        Z = self._Zr(p, q, r)
        if not Z:
            return [], []
        denom = self._Zr(p - 1, q + 1, r - 1) + self._bdry(p, q, r)
        idx = lm.independent_modulo(Z, denom, self.dom)
        reps = [list(Z[i]) for i in idx]
        return reps, denom

    # -- pages --------------------------------------------------------------- #
    def page(self, r):
        r = int(r)
        if r < 0:
            raise QuiverlabError(f"spectral sequence: page index r={r} < 0")
        if r in self._pages:
            return self._pages[r]
        cells = {}
        for (p, q) in self._candidates:
            reps, _ = self._cell(p, q, r)
            cells[(p, q)] = Subquotient(p, q, reps)
        page = Page(self, r, cells)
        self._pages[r] = page
        return page

    def _dr_matrix(self, p, q, r):
        """Matrix of ``d^r : E^r_{p,q} -> E^r_{p-r, q+r-1}`` (rows = target reps,
        cols = source reps). For each source rep ``x``, ``d(x)`` lands in
        ``F_{p-r} C_{n-1}`` (guaranteed by ``Z^r``); express it in the target's rep
        basis modulo the target's denominator, canonicalized by
        ``reduce_mod_nullspace``."""
        dom = self.dom
        src_reps, _ = self._cell(p, q, r)
        tp, tq = p - r, q + r - 1
        tgt_reps, tgt_denom = self._cell(tp, tq, r)
        nrows = len(tgt_reps)
        if not src_reps:
            return lm.zeros(nrows, 0, dom)
        n = p + q
        dn = self.F.dmat(n)
        basis = tgt_reps + tgt_denom
        B = sub.colmat(basis)
        cols = []
        for x in src_reps:
            dx = lm.matvec(dn, x, dom) if (dn and dn[0]) else [dom.zero()] * self.F.dim(n - 1)
            if all(dom.is_zero(v) for v in dx):
                cols.append([dom.zero()] * nrows)
                continue
            if not B:                             # target E^r = 0 but d(x) != 0 in C_{n-1}
                cols.append([dom.zero()] * nrows)
                continue
            coeff = flinalg.solve(B, dx, dom)
            if coeff is None:
                raise QuiverlabError(
                    f"spectral sequence: d^{r}(E_[{p},{q}]) image is not a class in "
                    f"E_[{tp},{tq}] -- page bookkeeping bug (report this)")
            coeff = flinalg.reduce_mod_nullspace(coeff, B, dom)
            cols.append(coeff[:nrows])
        return lm.cols_to_matrix(cols) if cols else lm.zeros(nrows, 0, dom)

    # -- convergence surface ------------------------------------------------- #
    @property
    def total_homology_dims(self):
        return self.F.total_homology_dims()

    # set by presets that certify a finite abutment window against an external
    # oracle (devil's-advocate fix, 2026-08-05: E_inf totals of a TRUNCATED
    # double complex are silently wrong past the certified degrees).
    certified_window = None

    def certified_abutment(self, n):
        """The abutment dimension at Ext-degree ``n``, readable ONLY inside the
        preset-certified window -- out-of-window reads refuse loudly instead of
        returning a truncation artifact."""
        from quiverlab.errors import QuiverlabError
        if self.certified_window is None:
            raise QuiverlabError(
                "this spectral sequence carries no certified abutment window; "
                "read total_homology_dims only with an external certificate")
        lo, hi = self.certified_window
        if not (lo <= n <= hi):
            raise QuiverlabError(
                f"abutment degree {n} is outside the certified window [{lo}, {hi}] "
                "-- a truncated double complex is silently wrong out there; "
                "raise p_len/q_len")
        return self.total_homology_dims.get(-n, 0)

    def einf_page(self):
        """The stabilized page ``E_inf`` (``page(e_infinity_page)``)."""
        return self.page(self.convergence.e_infinity_page)

    def __repr__(self):
        return (f"SpectralSequence(width={self.width}, height={self.height}, "
                f"abutment={self.total_homology_dims})")
