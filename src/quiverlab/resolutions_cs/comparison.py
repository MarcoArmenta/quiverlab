"""Comparison morphisms CS <-> bar, and windowed transport of the Tamarkin-Tsygan
cup/bracket classes (spec Plan-04 Task 12).

The Chouhy-Solotar resolution P_* is a retract of the (reduced) normalized bar
resolution B_*: an inclusion Phi: P -> B and a projection Psi: B -> P with
Psi Phi = id.  On generators the inclusion is the block map
    Phi_n(1 (x) sigma (x) 1) = 1 (x) u_0 (x) ... (x) u_{n-1} (x) 1
of the CS left decomposition sigma = u_0 | ... | u_{n-1}, PLUS the lower-order
correction that makes it a chain map: for a relation generator the image is the
whole relation, tip - beta(tip), not just the leading tip.  (For a monomial
presentation beta(tip)=0, so the block term alone is already a chain map; for a
quadratic binomial the correction is the length-2 words of the normal form.)

Applying Hom_{A^e}(-, A) turns the covariant inclusion into the contravariant
comparison COCHAIN map that this module works with,
    Phi# : C^n_bar -> C^n_cs ,   (Phi# g)(sigma) = g(Phi_n sigma),
a genuine cochain map (delta_cs Phi# = Phi# delta_bar) and a quasi-isomorphism.
Its inverse on HH^* transports classes the other way; cup/bracket are computed
on the bar side with engine.tt_calculus (the GF(p) facade) and pulled back through
Phi#.  (Cap products live on the bar engine, engine.tt_calculus, only; there is no
CS transport wrapper for them yet.)  Everything is exact over GF(p); the transport is delivered only inside the
bar-comparison WINDOW (a cup at degree n needs bar cochains up to degree 2n+1).

Bases.  The bar side uses the engine's cochain basis engine.scan3.cochain_basis
(pairs (w, s), interior word w in R^n, output slot s), matching tt_calculus.  The
CS side uses ChouhySolotarResolution._basis(n, "coh") = [(sigma, j)].  The unit
of A is a single basis vector only after unit_adaptation, so A-vectors coming from
the (original-basis) CS arithmetic are pushed through the change of basis P (its
columns are the unit-adapted basis in the original coordinates) before they meet
the bar side.
"""
import numpy as np

from quiverlab.errors import FieldError
from quiverlab.fields.linalg import nullspace, solve
from quiverlab.fields.primefield import PrimeField

_WINDOW_MSG = ("cup/bracket transport is delivered only within the bar-comparison "
               "window; native CS cup is a later phase")


class CSClass:
    """A CS-side Hochschild cohomology class: a degree plus a representative cocycle
    (coordinate vector over ChouhySolotarResolution._basis(degree, "coh"))."""
    __slots__ = ("degree", "vec")

    def __init__(self, degree, vec):
        self.degree = degree
        self.vec = list(vec)


class Comparison:
    def __init__(self, A, window=None, max_cells=4_000_000):
        if not isinstance(A.domain, PrimeField):
            raise FieldError(
                f"the CS<->bar comparison uses engine.tt_calculus (the GF(p) facade); "
                f"this algebra is over {A.domain.name}",
                hint="construct the algebra over GF(p); the domain-generic bar/CS "
                     "(co)homology dims serve any field")
        if A.quiver is None or A.relations is None:
            raise ValueError("Comparison needs an algebra built by Quiver.algebra "
                             "(quiver + relations present)")
        from quiverlab.engine.adapter import to_engine
        from quiverlab.resolutions_cs.build import reduction_system_of
        from quiverlab.resolutions_cs.homology import _require_admissible

        self.A = A
        self.dom = A.domain
        self.p = self.dom.p
        self.max_cells = max_cells
        self.B = A.unit_adapted()
        self.m = self.B.dim
        self.E = to_engine(self.B)
        self._rs = reduction_system_of(A)
        _require_admissible(self._rs)
        self._monomial = all(len(r.tail) == 0 for r in self._rs.rules)
        self._build_change_of_basis()
        self.window = self._default_window() if window is None else window
        # lazily grown CS resolution + degreewise caches
        self._res = None
        self._maxdeg = -1
        self._phi_cache = {}
        self._dcs_cache = {}
        self._dbar_cache = {}
        self._zbar_cache = {}
        self._barbasis_cache = {}

    # -- setup ---------------------------------------------------------------
    def _default_window(self):
        """Largest n with the bar cochains needed for a degree-n cup (up to degree
        2n+1) inside max_cells, under a conservative growth model that floors the
        reduced branching at 2 -- so even algebras with trivial reduced growth
        (k[x]/x^2 has m-1 = 1) get a finite, well-defined comparison window."""
        eff = max(self.m - 1, 2)
        n = 0
        while self.m * eff ** (2 * n + 1) <= self.max_cells:
            n += 1
        return max(n - 1, 0)

    def _build_change_of_basis(self):
        """P: columns = unit-adapted basis in the original coordinates (int mod p);
        Pinv its inverse.  Reproduces Algebra.unit_adapted's construction exactly so
        the CS (original-basis) and bar (unit-adapted-basis) sides agree."""
        dom, m, p = self.dom, self.m, self.p
        if self.A.is_unit_adapted:
            P = [[1 if r == c else 0 for c in range(m)] for r in range(m)]
        else:
            j0 = next(i for i, c in enumerate(self.A.unit) if not dom.is_zero(c))
            P = [[1 if r == c else 0 for c in range(m)] for r in range(m)]
            for r in range(m):
                P[r][j0] = int(self.A.unit[r]) % p
            if j0 != 0:
                for r in range(m):
                    P[r][0], P[r][j0] = P[r][j0], P[r][0]
        Pd = [[dom.coerce(P[i][j]) for j in range(m)] for i in range(m)]
        Pinv = [[0] * m for _ in range(m)]
        for k in range(m):
            ek = [dom.one() if i == k else dom.zero() for i in range(m)]
            x = solve(Pd, ek, dom)
            if x is None:
                raise ValueError("unit-adaptation change of basis is singular")
            for i in range(m):
                Pinv[i][k] = int(x[i]) % p
        self._P = [[int(P[i][j]) % p for j in range(m)] for i in range(m)]
        self._Pinv = Pinv

    def _ensure(self, need):
        """Grow the CS resolution so that S(need) is available; clear caches that
        depend on the resolution instance."""
        if self._res is not None and self._maxdeg >= need:
            return
        from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
        self._res = ChouhySolotarResolution(self.A, self._rs, max_degree=need + 1,
                                            max_cells=self.max_cells)
        self._maxdeg = need + 1
        self._phi_cache.clear()
        self._dcs_cache.clear()

    # -- A-vector helpers ----------------------------------------------------
    def _A_to_Bred(self, vA):
        """Original-basis A-int-vector -> {reduced unit-adapted index: coeff}
        (drops the unit component, index 0 of B)."""
        m, p, Pinv = self.m, self.p, self._Pinv
        out = {}
        for s in range(1, m):
            acc = sum(Pinv[s][k] * (vA[k] % p) for k in range(m)) % p
            if acc:
                out[s] = acc
        return out

    def _block_reduced_index(self, word):
        """A CS block (arrow-name path) as a single reduced unit-adapted basis index.
        Blocks in the CS left decomposition of a quadratic/monomial presentation are
        irreducible single-basis paths, so their normal form is one basis vector."""
        ar = self._res.ar
        vA = [int(x) % self.p for x in ar.path_vec(word)]
        d = self._A_to_Bred(vA)
        if len(d) != 1 or next(iter(d.values())) != 1:
            raise NotImplementedError(
                "comparison block does not reduce to a single basis path; the honest "
                "Phi for this (non-quadratic/non-monomial) presentation needs the "
                "Skoldberg homotopy expansion (a later phase)")
        return next(iter(d))

    # -- the comparison expansion Phi_n(sigma) = sum coeff * (1 (x) word (x) 1) ----
    def _expansion(self, n, sigma):
        """{interior word (reduced unit-adapted indices): coeff} for the honest bar
        image Phi_n(sigma).  Leading block term for every n; for a non-monomial
        presentation at n = 2 also the -beta(tip) relation correction."""
        if n == 0:
            return {(): 1}
        lead = tuple(self._block_reduced_index(b) for b in sigma.blocks)
        out = {lead: 1}
        if not self._monomial and n >= 2:
            nf = self._rs.normal_form(sigma.word)   # beta(tip) = sum lambda_k v_k
            for w, coeff in nf.items():
                if len(w) != len(sigma.word):
                    raise NotImplementedError(
                        "inhomogeneous reduction in the comparison correction; the "
                        "honest Phi past the quadratic window needs the Skoldberg "
                        "homotopy expansion (a later phase)")
                if n != 2:
                    raise NotImplementedError(
                        "the honest comparison correction is delivered through degree "
                        "2 for non-monomial presentations (quadratic window)")
                key = tuple(self._block_reduced_index((letter,)) for letter in w)
                out[key] = (out.get(key, 0) - int(coeff)) % self.p
        return out

    # -- bar-side bases / matrices (engine convention) -----------------------
    def _bar_basis(self, n):
        b = self._barbasis_cache.get(n)
        if b is None:
            from quiverlab.engine.scan3 import cochain_basis
            b = cochain_basis(self.E, n)
            self._barbasis_cache[n] = b
        return b

    def _dbar(self, n):
        """Bar coboundary delta^n : C^n_bar -> C^{n+1}_bar (int mod p)."""
        M = self._dbar_cache.get(n)
        if M is None:
            from quiverlab.engine.scan3 import coboundary_matrix
            bn = self._bar_basis(n)
            idx = {g: i for i, g in enumerate(self._bar_basis(n + 1))}
            D = coboundary_matrix(self.E, n, bn, idx)
            M = [[int(D[i][j]) % self.p for j in range(D.shape[1])]
                 for i in range(D.shape[0])]
            self._dbar_cache[n] = M
        return M

    def _dcs(self, n):
        """CS coboundary delta^n : C^n_cs -> C^{n+1}_cs (int mod p)."""
        M = self._dcs_cache.get(n)
        if M is None:
            self._ensure(n + 1)
            raw = self._res.matrix(n, "coh")
            M = [[int(x) % self.p for x in row] for row in raw]
            self._dcs_cache[n] = M
        return M

    def _zbar(self, n):
        """A basis of bar cocycles Z^n_bar = ker delta^n_bar (int mod p vectors)."""
        Z = self._zbar_cache.get(n)
        if Z is None:
            dom = self.dom
            dbar = self._dbar(n)
            ker = nullspace([[dom.coerce(x) for x in row] for row in dbar], dom)
            Z = [[int(x) % self.p for x in v] for v in ker]
            self._zbar_cache[n] = Z
        return Z

    # -- Phi# : C^n_bar -> C^n_cs (the comparison cochain map) ---------------
    def Phi(self, n):
        """Matrix of the comparison cochain map Phi# : C^n_bar -> C^n_cs
        (rows = CS coh basis, cols = engine bar cochain basis), int mod p."""
        M = self._phi_cache.get(n)
        if M is not None:
            return M
        self._ensure(n)
        p = self.p
        csb = self._res._basis(n, "coh")
        bcb = self._bar_basis(n)
        Pi = self._P
        M = [[0] * len(bcb) for _ in range(len(csb))]
        exps = {}
        for ri, (sigma, j) in enumerate(csb):
            key = sigma.word
            if key not in exps:
                exps[key] = self._expansion(n, sigma)
            ex = exps[key]
            for ci, (w, s) in enumerate(bcb):
                c = ex.get(tuple(w), 0)
                if c:
                    M[ri][ci] = (M[ri][ci] + c * Pi[j][s]) % p
        self._phi_cache[n] = M
        return M

    def Psi(self, n):
        """The inverse comparison on cohomology, HH^n_cs -> HH^n_bar, as a matrix in
        the (CS / bar) cohomology-representative bases.  Psi = (Phi#|_{HH})^{-1}: the
        Hom-dual of the CS projection Psi: B -> P, delivered on HH^* (the resolution-
        level projection is the Skoldberg splitting, a later phase)."""
        Fstar, _, _ = self._induced_on_cohomology(n)
        inv = self._invert_modp(Fstar)
        if inv is None:
            raise AssertionError(
                f"comparison is not invertible on HH^{n}: Phi# fails to be a "
                f"quasi-isomorphism at degree {n} (a bug, never an approximation)")
        return inv

    # -- linear-algebra utilities (int mod p) --------------------------------
    def _matmul(self, Aa, Bb):
        p = self.p
        if not Aa or not Bb or not Bb[0]:
            r = len(Aa)
            c = len(Bb[0]) if Bb and Bb[0] else 0
            return [[0] * c for _ in range(r)]
        r, inn, c = len(Aa), len(Bb), len(Bb[0])
        return [[sum(Aa[i][k] * Bb[k][j] for k in range(inn)) % p for j in range(c)]
                for i in range(r)]

    def _matvec(self, M, v):
        p = self.p
        if not M:
            return []
        return [sum(M[i][k] * v[k] for k in range(len(v))) % p for i in range(len(M))]

    def _cols(self, M):
        return [[row[c] for row in M] for c in range(len(M[0]))] if M and M[0] else []

    def _invert_modp(self, M):
        """Inverse of a square int-mod-p matrix, or None if singular."""
        dom = self.dom
        d = len(M)
        if d == 0:
            return []
        cols = []
        for k in range(d):
            ek = [dom.one() if i == k else dom.zero() for i in range(d)]
            Md = [[dom.coerce(M[i][j]) for j in range(d)] for i in range(d)]
            x = solve(Md, ek, dom)
            if x is None:
                return None
            cols.append([int(x[i]) % self.p for i in range(d)])
        return [[cols[j][i] for j in range(d)] for i in range(d)]

    # -- transport of cocycles ----------------------------------------------
    def transport_class_bar_to_cs(self, bar_cocycle, n):
        """A bar cocycle -> the CS cocycle Phi#(bar_cocycle).  Direct: Phi# is the
        Hom-dual (transpose) of the inclusion chain map."""
        return self._matvec(self.Phi(n), [int(x) % self.p for x in bar_cocycle])

    def transport_cocycle_cs_to_bar(self, cs_cocycle, n):
        """A CS cocycle c -> a bar cocycle g with Phi#(g) cohomologous to c (the
        inverse comparison, realised by solve against the bar cocycles and the CS
        coboundaries)."""
        dom, p = self.dom, self.p
        self._ensure(n)
        c = [int(x) % p for x in cs_cocycle]
        Z = self._zbar(n)
        PhiZ = [self._matvec(self.Phi(n), g) for g in Z]
        Bcs = self._cols(self._dcs(n - 1)) if n >= 1 else []
        columns = PhiZ + Bcs
        nrows = len(c)
        if not columns:
            if any(x % p for x in c):
                raise AssertionError("cs->bar transport: nonzero class in a zero space")
            return [0] * len(self._bar_basis(n))
        Mat = [[dom.coerce(columns[cc][rr]) for cc in range(len(columns))]
               for rr in range(nrows)]
        rhs = [dom.coerce(c[rr]) for rr in range(nrows)]
        x = solve(Mat, rhs, dom)
        if x is None:
            raise AssertionError(
                f"cs->bar transport is inconsistent at degree {n}: the comparison is "
                f"not a quasi-isomorphism here (a bug, never an approximation)")
        a = [int(x[i]) % p for i in range(len(Z))]
        g = [0] * len(self._bar_basis(n))
        for i, ai in enumerate(a):
            if ai:
                Zi = Z[i]
                for k in range(len(g)):
                    g[k] = (g[k] + ai * Zi[k]) % p
        return g

    # -- cohomology plumbing -------------------------------------------------
    def cs_cohomology_basis(self, n):
        """Representative CS cocycles of a basis of HH^n_cs (int mod p vectors)."""
        from quiverlab.resolutions_cs.homology import cs_hh_basis
        reps = cs_hh_basis(self.A, n, "coh", max_cells=self.max_cells)
        return [[int(x) % self.p for x in v] for v in reps]

    def hh_class_cs(self, n, i):
        """The i-th representative class of HH^n on the CS side."""
        reps = self.cs_cohomology_basis(n)
        if i >= len(reps):
            raise IndexError(f"HH^{n}(CS) has {len(reps)} basis classes; no index {i}")
        return CSClass(n, reps[i])

    def _induced_on_cohomology(self, n):
        """Matrix of the induced Phi# : HH^n_bar -> HH^n_cs in chosen representative
        bases, together with those bases.  Returns (Fstar, cs_reps, bar_reps)."""
        cs_reps = self.cs_cohomology_basis(n)
        bar_reps = self._bar_cohomology_basis(n)
        if len(cs_reps) != len(bar_reps):
            raise AssertionError(
                f"HH^{n} dimension mismatch CS={len(cs_reps)} bar={len(bar_reps)}: the "
                f"comparison cannot be an isomorphism (a bug)")
        d = len(cs_reps)
        Fstar = [[0] * d for _ in range(d)]
        for j, g in enumerate(bar_reps):
            img = self._matvec(self.Phi(n), g)         # a CS cocycle
            coords = self._class_coords_cs(n, img, cs_reps)
            for i in range(d):
                Fstar[i][j] = coords[i]
        return Fstar, cs_reps, bar_reps

    def _bar_cohomology_basis(self, n):
        """Representative bar cocycles of a basis of HH^n_bar (int mod p vectors)."""
        from quiverlab.engine import tt_calculus as TT
        H = TT.cohomology_classes(self.E, n, self.p)
        return [[int(H.reps[r, k]) % self.p for r in range(H.reps.shape[0])]
                for k in range(H.dim)]

    def _class_coords_cs(self, n, cocycle, cs_reps):
        """Coordinates of a CS cocycle in the cohomology basis cs_reps (mod the CS
        coboundaries)."""
        dom, p = self.dom, self.p
        Bcs = self._cols(self._dcs(n - 1)) if n >= 1 else []
        columns = [list(r) for r in cs_reps] + Bcs
        nrows = len(cocycle)
        if not columns:
            if any(x % p for x in cocycle):
                raise AssertionError("class coords: nonzero cocycle in a zero space")
            return []
        Mat = [[dom.coerce(columns[cc][rr]) for cc in range(len(columns))]
               for rr in range(nrows)]
        rhs = [dom.coerce(int(cocycle[rr]) % p) for rr in range(nrows)]
        x = solve(Mat, rhs, dom)
        if x is None:
            raise AssertionError(
                "a transported bar cocycle is not a CS cocycle (descent failed -- a bug)")
        return [int(x[i]) % p for i in range(len(cs_reps))]

    def same_cohomology_class(self, x, y, degree):
        """True iff x - y is a CS coboundary at the given degree (fields.linalg.solve
        against im delta_cs^{degree-1})."""
        dom, p = self.dom, self.p
        self._ensure(degree)
        diff = [(int(x[i]) - int(y[i])) % p for i in range(len(x))]
        Bcs = self._cols(self._dcs(degree - 1)) if degree >= 1 else []
        if not Bcs:
            return all(d % p == 0 for d in diff)
        Mat = [[dom.coerce(Bcs[c][r]) for c in range(len(Bcs))] for r in range(len(diff))]
        rhs = [dom.coerce(d) for d in diff]
        return solve(Mat, rhs, dom) is not None

    # -- windowed transport of cup / bracket ---------------------------------
    def _check_window(self, *degrees):
        if any(d > self.window for d in degrees):
            raise NotImplementedError(_WINDOW_MSG)

    def _bar_cochain_of(self, u):
        """CSClass -> its transported bar cocycle (numpy int64) at degree u.degree."""
        g = self.transport_cocycle_cs_to_bar(u.vec, u.degree)
        return np.array(g, dtype=np.int64)

    def cup_of_cs_classes(self, u, v):
        """Psi*( Phi*(u) cup_bar Phi*(v) ): the CS-side cup of two CS classes,
        computed by transporting representatives to the bar, cupping there with the
        Tamarkin-Tsygan calculus (engine.tt_calculus.cup_cochain), and pulling the
        result back through Phi#.  Returns a CS cochain at degree u.degree+v.degree."""
        from quiverlab.engine import tt_calculus as TT
        p, q = u.degree, v.degree
        self._check_window(p, q)
        gu, gv = self._bar_cochain_of(u), self._bar_cochain_of(v)
        cup = TT.cup_cochain(self.E, p, q, gu, gv)          # bar (p+q)-cochain
        self._ensure(p + q)
        return self.transport_class_bar_to_cs([int(x) % self.p for x in cup], p + q)

    def bracket_of_cs_classes(self, u, v):
        """Psi*( [Phi*(u), Phi*(v)]_bar ): the CS-side Gerstenhaber bracket, via
        engine.tt_calculus.gerstenhaber_bracket_cochain.  Returns a CS cochain at
        degree u.degree+v.degree-1."""
        from quiverlab.engine import tt_calculus as TT
        p, q = u.degree, v.degree
        self._check_window(p, q)
        gu, gv = self._bar_cochain_of(u), self._bar_cochain_of(v)
        br = TT.gerstenhaber_bracket_cochain(self.E, p, q, gu, gv)  # bar (p+q-1)-cochain
        self._ensure(p + q - 1)
        return self.transport_class_bar_to_cs([int(x) % self.p for x in br], p + q - 1)

    def transport_then_bar_cup(self, u, v):
        """A second, engine-native route to the transported cup: transport u, v to
        bar CLASSES, cup them with the induced structure-constant tensor
        (tt_calculus.cup_product_matrix), reconstruct a bar cocycle representative
        from the resulting class, and pull it back through Phi#.  Cohomologous to
        cup_of_cs_classes -- a cross-check of the cochain-level and class-level cups."""
        from quiverlab.engine import tt_calculus as TT
        p, q = u.degree, v.degree
        self._check_window(p, q)
        gu, gv = self._bar_cochain_of(u), self._bar_cochain_of(v)
        Hp = TT.cohomology_classes(self.E, p, self.p)
        Hq = TT.cohomology_classes(self.E, q, self.p)
        Hpq = TT.cohomology_classes(self.E, p + q, self.p)
        cu = Hp.coords(gu % self.p)
        cv = Hq.coords(gv % self.p)
        C, dp, dq, dpq = TT.cup_product_matrix(self.E, p, q, self.p)
        cls = np.zeros(dpq, dtype=np.int64)
        for k in range(dpq):
            acc = 0
            for i in range(dp):
                for j in range(dq):
                    acc += int(C[k, i, j]) * int(cu[i]) * int(cv[j])
            cls[k] = acc % self.p
        rep = np.zeros(Hpq.reps.shape[0], dtype=np.int64)
        for k in range(Hpq.dim):
            rep = (rep + int(cls[k]) * Hpq.reps[:, k]) % self.p
        self._ensure(p + q)
        return self.transport_class_bar_to_cs([int(x) % self.p for x in rep], p + q)

    # -- assertions ----------------------------------------------------------
    def assert_chain_map(self, upto=2):
        """The comparison morphism intertwines the two differentials, as matrix
        identities over GF(p):

          (1) delta_cs^n Phi#_n = Phi#_{n+1} delta_bar^n   (Phi# is a cochain map;
              the Hom form of b_bar Phi = Phi d_cs), for n = 0..upto-1; and
          (2) on HH^k the induced Phi* is a two-sided isomorphism -- its inverse Psi*
              (the CS projection on cohomology) satisfies Phi* Psi* = Psi* Phi* = id
              (the Hom form of d_cs Psi = Psi b_bar), for k = 0..upto.

        The identities are asserted only where genuinely non-vacuous: (1) requires a
        nonzero delta on BOTH sides at some tested degree (guarded), and (2) is a
        matrix identity on the actual (nonzero) HH bases.
        """
        saw_nontrivial = False
        for n in range(upto):
            dcs = self._dcs(n)
            dbar = self._dbar(n)
            lhs = self._matmul(dcs, self.Phi(n))            # C^n_bar -> C^{n+1}_cs
            rhs = self._matmul(self.Phi(n + 1), dbar)
            if lhs != rhs:
                raise AssertionError(
                    f"comparison is not a cochain map at degree {n}: "
                    f"delta_cs Phi# != Phi# delta_bar")
            nz_cs = any(v % self.p for row in dcs for v in row)
            nz_bar = any(v % self.p for row in dbar for v in row)
            if nz_cs and nz_bar:
                saw_nontrivial = True
        if upto >= 1 and not saw_nontrivial:
            raise AssertionError(
                "assert_chain_map is vacuous: no tested degree has nonzero "
                "differentials on both sides")
        for k in range(upto + 1):
            Fstar, cs_reps, _ = self._induced_on_cohomology(k)
            inv = self._invert_modp(Fstar)
            if inv is None:
                raise AssertionError(
                    f"comparison is not invertible on HH^{k}: Phi* has no inverse "
                    f"(not a quasi-isomorphism -- a bug)")
            d = len(cs_reps)
            ident = [[1 if i == j else 0 for j in range(d)] for i in range(d)]
            if self._matmul(Fstar, inv) != ident or self._matmul(inv, Fstar) != ident:
                raise AssertionError(
                    f"Phi* Psi* != id on HH^{k}: the CS projection is not a genuine "
                    f"inverse comparison")

    def assert_transport_roundtrip_identity(self, upto=3):
        """Psi*(Phi*(c)) = c on HH^n for every basis class c (n = 0..upto): each CS
        class is transported to the bar and pulled straight back, and must land in
        its own CS cohomology class.  Asserted on genuinely nonzero HH^n."""
        for n in range(upto + 1):
            reps = self.cs_cohomology_basis(n)
            for c in reps:
                g = self.transport_cocycle_cs_to_bar(c, n)         # Phi*(c): to bar
                back = self.transport_class_bar_to_cs(g, n)         # Psi*(...): to CS
                if not self.same_cohomology_class(back, c, degree=n):
                    raise AssertionError(
                        f"transport round-trip Psi*Phi* != id on HH^{n}: a class "
                        f"failed to return to itself (a bug, never an approximation)")
