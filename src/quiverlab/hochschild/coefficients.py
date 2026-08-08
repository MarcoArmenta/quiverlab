"""Exact-Domain A-bimodule coefficients for Hochschild (co)homology (Plan 52).

A :class:`Bimodule` is a light ``(Lact, Ract)``-over-``A`` container (DD1): the
coefficient of ``HH^•(A, M)`` / ``HH_•(A, M)``. It stores, for a
:class:`quiverlab.core.Algebra` ``A`` of dimension ``m`` over an exact Domain:

- ``dim_M`` — ``dim_k M``;
- ``Lact`` — a length-``m`` list of ``dim_M × dim_M`` Domain matrices, ROW-indexed
  by the input: ``Lact[j][s]`` is the M-vector ``b_j · m_s`` (left action by the
  ``j``-th ``A``-basis element). Equivalently ``b_j · v = v @ Lact[j]`` for a row
  vector ``v``;
- ``Ract`` — likewise, ``Ract[j][s]`` is ``m_s · b_j`` (right action);
- ``algebra``, ``domain``, ``name``.

This is the exact-Domain twin of the ported GF(p)/int64
``quiverlab.engine.bimodule`` (a read-only reference oracle). The engine stores its
actions COLUMN-indexed (``Lact[j][:, s] = b_j · m_s``); :meth:`same_actions_mod_p`
compares across the transpose.

Named constructors — :meth:`regular`, :meth:`dual` (``D(A) = Hom_k(A, k)``),
:meth:`twisted` (``{}_φA_ψ``), :meth:`twisted_by_nakayama` (``{}_1A_ν``),
:meth:`mod_socle` (``A/soc_{A^e}A``), :meth:`from_actions` (generic no-code).

Downstream consumers (named for the record):
- **P54 (BV operator)** — :meth:`twisted_by_nakayama` + the chain-level accessor
  :func:`twisted_homology_classes` (the cross-plan contract).
- **P72 (split-extension LES)** — ``B.hochschild_cohomology(top, coefficients=M)``.
- **P74 (skew group algebras)** — :meth:`twisted` with ``phi=g_matrix``.
- generic no-code bimodules — :meth:`from_actions`.
"""

from quiverlab.errors import QuiverlabError
from quiverlab.fields.domain import reject_inexact
from quiverlab.fields.linalg import nullspace, rank, solve

# Allocation guard (mirrors modules' _MAX_MODULE_DIM): a generic from_actions
# coefficient of enormous dim would build dim_M^2 * m Domain entries.
_MAX_COEFF_DIM = 400


def _mm(dom, X, Y):
    """Standard matrix product (X @ Y) over the Domain: (X@Y)[a][c] = sum_b X[a][b] Y[b][c]."""
    n, k = len(X), len(X[0]) if X else 0
    w = len(Y[0]) if Y else 0
    out = [[dom.zero()] * w for _ in range(n)]
    for a in range(n):
        Xa = X[a]
        for b in range(k):
            xab = Xa[b]
            if dom.is_zero(xab):
                continue
            Yb = Y[b]
            outa = out[a]
            for c in range(w):
                yc = Yb[c]
                if not dom.is_zero(yc):
                    outa[c] = dom.add(outa[c], dom.mul(xab, yc))
    return out


def _scaled_sum(dom, coeffs, mats, nrows, ncols):
    """sum_c coeffs[c] * mats[c], a matrix of shape (nrows, ncols)."""
    out = [[dom.zero()] * ncols for _ in range(nrows)]
    for c, cf in enumerate(coeffs):
        if dom.is_zero(cf):
            continue
        Mc = mats[c]
        for a in range(nrows):
            outa, Ma = out[a], Mc[a]
            for b in range(ncols):
                mb = Ma[b]
                if not dom.is_zero(mb):
                    outa[b] = dom.add(outa[b], dom.mul(cf, mb))
    return out


def _eq_mat(dom, X, Y):
    return all(dom.eq(X[a][b], Y[a][b]) for a in range(len(X)) for b in range(len(X[0])))


class Bimodule:
    def __init__(self, algebra, dim_M, Lact, Ract, name="M"):
        self.algebra = algebra
        self.domain = algebra.domain
        self.dim_M = dim_M
        self.Lact = Lact          # list[m] of dim_M x dim_M Domain matrices (row s = b_j . m_s)
        self.Ract = Ract          # list[m] of dim_M x dim_M Domain matrices (row s = m_s . b_j)
        self.name = name

    # -- linear action helpers (row convention: b_j . v = v @ Lact[j]) --------
    def _unit_vec(self, s):
        dom = self.domain
        v = [dom.zero()] * self.dim_M
        v[s] = dom.one()
        return v

    def left_apply(self, b_vec, v):
        """b · v where b = sum_j b_vec[j] b_j (an A-vector) and v an M-vector."""
        dom = self.domain
        out = [dom.zero()] * self.dim_M
        for j, bj in enumerate(b_vec):
            if dom.is_zero(bj):
                continue
            Lj = self.Lact[j]
            for s, vs in enumerate(v):
                if dom.is_zero(vs):
                    continue
                c = dom.mul(bj, vs)
                row = Lj[s]
                for t in range(self.dim_M):
                    w = row[t]
                    if not dom.is_zero(w):
                        out[t] = dom.add(out[t], dom.mul(c, w))
        return out

    def right_apply(self, a_vec, v):
        """v · a where a = sum_j a_vec[j] b_j (an A-vector) and v an M-vector."""
        dom = self.domain
        out = [dom.zero()] * self.dim_M
        for j, aj in enumerate(a_vec):
            if dom.is_zero(aj):
                continue
            Rj = self.Ract[j]
            for s, vs in enumerate(v):
                if dom.is_zero(vs):
                    continue
                c = dom.mul(aj, vs)
                row = Rj[s]
                for t in range(self.dim_M):
                    w = row[t]
                    if not dom.is_zero(w):
                        out[t] = dom.add(out[t], dom.mul(c, w))
        return out

    # -- certification --------------------------------------------------------
    def check(self):
        """Verify the A^e-module axioms (exact, loud). Left action is an algebra
        map, right action an anti-map, the two commute, and the unit acts as the
        identity on both sides. Raises QuiverlabError on any failure."""
        A = self.algebra
        dom = self.domain
        m = A.dim
        dm = self.dim_M
        for j in range(m):
            for mat, side in ((self.Lact, "left"), (self.Ract, "right")):
                if len(mat) != m or len(mat[j]) != dm or any(len(r) != dm for r in mat[j]):
                    raise QuiverlabError(
                        f"bimodule {side} action tables have the wrong shape "
                        f"(need m={m} matrices of size {dm}x{dm})",
                        hint="build actions with the constructors or from_actions")
        # associativity: Lact[j] @ Lact[i] == sum_c (b_i b_j)[c] Lact[c]
        #                Ract[i] @ Ract[j] == sum_c (b_i b_j)[c] Ract[c]
        for i in range(m):
            for j in range(m):
                bij = A.T[i][j]
                lhsL = _mm(dom, self.Lact[j], self.Lact[i])
                rhsL = _scaled_sum(dom, bij, self.Lact, dm, dm)
                if not _eq_mat(dom, lhsL, rhsL):
                    raise QuiverlabError(
                        f"bimodule left action is not associative at (b{i}, b{j})",
                        hint="a·(b·m) must equal (a b)·m — check the action tables")
                lhsR = _mm(dom, self.Ract[i], self.Ract[j])
                rhsR = _scaled_sum(dom, bij, self.Ract, dm, dm)
                if not _eq_mat(dom, lhsR, rhsR):
                    raise QuiverlabError(
                        f"bimodule right action is not associative at (b{i}, b{j})",
                        hint="(m·a)·b must equal m·(a b) — check the action tables")
                # commute: Lact[i] @ Ract[j] == Ract[j] @ Lact[i]
                if not _eq_mat(dom, _mm(dom, self.Lact[i], self.Ract[j]),
                               _mm(dom, self.Ract[j], self.Lact[i])):
                    raise QuiverlabError(
                        f"bimodule left/right actions do not commute at (b{i}, b{j})",
                        hint="(a·m)·b must equal a·(m·b)")
        # unit acts as identity on both sides
        unit = list(A.unit)
        eye = [[dom.one() if a == b else dom.zero() for b in range(dm)] for a in range(dm)]
        Lu = _scaled_sum(dom, unit, self.Lact, dm, dm)
        Ru = _scaled_sum(dom, unit, self.Ract, dm, dm)
        if not _eq_mat(dom, Lu, eye) or not _eq_mat(dom, Ru, eye):
            raise QuiverlabError(
                "the unit of A does not act as the identity on the bimodule",
                hint="1·m = m = m·1 must hold")
        return True

    # -- degree-0 closed forms (Task 6) ---------------------------------------
    def _diff_rows(self):
        """Rows { (Lact[j]-Ract[j])[s] : j, s } — the columns of delta^0 / the
        span [A, M] (both HH^0 and HH_0 read these vectors)."""
        dom = self.domain
        rows = []
        for j in range(self.algebra.dim):
            Lj, Rj = self.Lact[j], self.Ract[j]
            for s in range(self.dim_M):
                rows.append([dom.sub(Lj[s][t], Rj[s][t]) for t in range(self.dim_M)])
        return rows

    def invariants_dim(self):
        """dim M^A = dim ker (m -> (a·m - m·a) over the generators). M^A is the
        left-nullspace of the (Lact-Ract) maps; dim = dim_M - rank of their columns."""
        dom = self.domain
        # constraint (j, t): sum_s m[s] (Lact[j]-Ract[j])[s][t] = 0
        rows = []
        for j in range(self.algebra.dim):
            Lj, Rj = self.Lact[j], self.Ract[j]
            for t in range(self.dim_M):
                rows.append([dom.sub(Lj[s][t], Rj[s][t]) for s in range(self.dim_M)])
        return self.dim_M - (rank(rows, dom) if rows else 0)

    def coinvariants_dim(self):
        """dim M/[A, M] = dim_M - dim [A, M], [A, M] = span{a·m - m·a}."""
        rows = self._diff_rows()
        return self.dim_M - (rank(rows, self.domain) if rows else 0)

    # -- engine bridge --------------------------------------------------------
    def same_actions_mod_p(self, engine_bimodule, p):
        """True iff this Bimodule equals the ported GF(p) ``engine.bimodule.Bimodule``
        mod p. The engine stores actions COLUMN-indexed, so we compare our row-index
        against its column-index (the transpose)."""
        eb = engine_bimodule
        if self.dim_M != int(eb.mu):
            return False
        for j in range(self.algebra.dim):
            Lj, Rj = self.Lact[j], self.Ract[j]
            ebL, ebR = eb.Lact[j], eb.Ract[j]
            for s in range(self.dim_M):
                for t in range(self.dim_M):
                    if int(Lj[s][t]) % p != int(ebL[t][s]) % p:
                        return False
                    if int(Rj[s][t]) % p != int(ebR[t][s]) % p:
                        return False
        return True

    def engine_actions(self, p):
        """Row-convention int64 action tensors reduced mod p: L[j][s][t] = coeff of
        m_t in b_j·m_s (GF(p) only; the minimal A^e engine cross-check)."""
        import numpy as np
        m = self.algebra.dim
        dm = self.dim_M
        L = np.zeros((m, dm, dm), dtype=np.int64)
        R = np.zeros((m, dm, dm), dtype=np.int64)
        for j in range(m):
            Lj, Rj = self.Lact[j], self.Ract[j]
            for s in range(dm):
                for t in range(dm):
                    L[j, s, t] = int(Lj[s][t]) % p
                    R[j, s, t] = int(Rj[s][t]) % p
        return L, R

    def describe(self):
        """Short provenance descriptor for the HHTable / report (e.g. 'D(A)')."""
        return self.name

    # -- change of basis (for the bar unit-adaptation) -----------------------
    def change_of_basis(self, P, newA):
        """Re-express the actions in a new A-basis given ``P`` (columns = new basis
        vectors in old coords): the action of the new basis element ``f_j`` is
        ``sum_c P[c][j] · (old action of b_c)``. M's own space is unchanged; only the
        A-element indexing moves. Used to move a coefficient (built in A's basis)
        into the unit-adapted basis the bar complex needs."""
        dom = self.domain
        m = self.algebra.dim
        dm = self.dim_M
        colP = [[P[c][j] for c in range(m)] for j in range(m)]   # colP[j] = column j of P
        newL = [_scaled_sum(dom, colP[j], self.Lact, dm, dm) for j in range(m)]
        newR = [_scaled_sum(dom, colP[j], self.Ract, dm, dm) for j in range(m)]
        return Bimodule(newA, dm, newL, newR, name=self.name)

    # -- constructors ---------------------------------------------------------
    @classmethod
    def regular(cls, A):
        """The regular bimodule ``M = A``: Lact = Ract = A.T, dim_M = m."""
        dom = A.domain
        m = A.dim
        Lact = [[[A.T[j][s][t] for t in range(m)] for s in range(m)] for j in range(m)]
        Ract = [[[A.T[s][j][t] for t in range(m)] for s in range(m)] for j in range(m)]
        return cls(A, m, Lact, Ract, name="A")

    @classmethod
    def dual(cls, A):
        """The dualizing bimodule ``D(A) = Hom_k(A, k)`` with ``(a·f·b)(x) = f(b x a)``.
        On the dual basis the actions are transposes of A's multiplication."""
        dom = A.domain
        m = A.dim
        # (e_j · phi_s) = sum_t (e_t e_j)[s] phi_t  =>  Lact[j][s][t] = A.T[t][j][s]
        # (phi_s · e_j) = sum_t (e_j e_t)[s] phi_t  =>  Ract[j][s][t] = A.T[j][t][s]
        Lact = [[[A.T[t][j][s] for t in range(m)] for s in range(m)] for j in range(m)]
        Ract = [[[A.T[j][t][s] for t in range(m)] for s in range(m)] for j in range(m)]
        return cls(A, m, Lact, Ract, name="D(A)")

    @classmethod
    def twisted(cls, A, phi=None, psi=None, name=None):
        """The twisted bimodule ``{}_φA_ψ``: underlying space A, ``a·x·b = φ(a) x ψ(b)``.
        ``phi``/``psi`` are algebra automorphisms as ``m×m`` Domain matrices
        (columns = images of the basis) IN A's OWN basis; ``None`` = identity
        (``φ=ψ=id`` is the regular bimodule)."""
        dom = A.domain
        m = A.dim
        phi = _identity_matrix(dom, m) if phi is None else _coerce_matrix(dom, phi, m)
        psi = _identity_matrix(dom, m) if psi is None else _coerce_matrix(dom, psi, m)
        # Lact[j][s] = phi(e_j) · e_s = sum_c phi[c][j] A.T[c][s]
        # Ract[j][s] = e_s · psi(e_j) = sum_c psi[c][j] A.T[s][c]
        Lact = []
        Ract = []
        for j in range(m):
            phicol = [phi[c][j] for c in range(m)]
            psicol = [psi[c][j] for c in range(m)]
            Lj, Rj = [], []
            for s in range(m):
                lvec = [dom.zero()] * m
                rvec = [dom.zero()] * m
                for c in range(m):
                    if not dom.is_zero(phicol[c]):
                        for t in range(m):
                            w = A.T[c][s][t]
                            if not dom.is_zero(w):
                                lvec[t] = dom.add(lvec[t], dom.mul(phicol[c], w))
                    if not dom.is_zero(psicol[c]):
                        for t in range(m):
                            w = A.T[s][c][t]
                            if not dom.is_zero(w):
                                rvec[t] = dom.add(rvec[t], dom.mul(psicol[c], w))
                Lj.append(lvec)
                Rj.append(rvec)
            Lact.append(Lj)
            Ract.append(Rj)
        if name is None:
            name = "twisted {}_φA_ψ"
        return cls(A, m, Lact, Ract, name=name)

    @classmethod
    def twisted_by_nakayama(cls, A):
        """``{}_1A_ν`` with ``ν`` the Nakayama automorphism. For self-injective A,
        ``D(A) ≅ {}_1A_ν``. (P54 entry point.)

        ``ν`` is obtained via the field-general trace-form route
        ``invariants.frobenius.nakayama_automorphism_generic`` in A's OWN basis —
        NOT ``Algebra.nakayama_automorphism()``, which over GF(p) returns the ENGINE
        unit-adapted basis and would be a silent bug against the A.T-built actions
        (DD1b)."""
        from quiverlab.invariants.frobenius import nakayama_automorphism_generic
        nu = nakayama_automorphism_generic(A)
        M = cls.twisted(A, psi=nu, name="twisted _1A_ν")
        return M

    @classmethod
    def mod_socle(cls, A):
        """The quotient ``A / soc_{A^e}A`` with the induced actions
        (``soc_{A^e}A = {x : (rad A)·x = 0 = x·(rad A)}``)."""
        from quiverlab.families.trivial_extension import _bimodule_socle
        from quiverlab.invariants.pathbasis import path_type_basis
        dom = A.domain
        m = A.dim
        try:
            _idem, rad, _src, _tgt = path_type_basis(A, "mod_socle coefficient")
        except QuiverlabError as exc:
            raise QuiverlabError(
                "A/soc coefficient needs a path-type (quiver) presentation to read "
                "off the radical", hint="present A via Quiver(...).algebra(...)") from exc
        socle = _bimodule_socle(A, rad, dom)              # coordinate vectors of soc
        comp, project = _quotient_projection(A, socle)    # complement indices + projector
        dq = len(comp)
        # induced actions on the quotient (row convention):
        #   Lact_M[j][s] = pi(b_j · rep(m_s)) = pi(A.T[j][comp[s]])
        #   Ract_M[j][s] = pi(rep(m_s) · b_j) = pi(A.T[comp[s]][j])
        Lact, Ract = [], []
        for j in range(m):
            Lj, Rj = [], []
            for s in range(dq):
                Lj.append(project([A.T[j][comp[s]][t] for t in range(m)]))
                Rj.append(project([A.T[comp[s]][j][t] for t in range(m)]))
            Lact.append(Lj)
            Ract.append(Rj)
        return cls(A, dq, Lact, Ract, name="A/soc")

    @classmethod
    def from_actions(cls, A, dim_M, left_maps, right_maps, name="M"):
        """Generic no-code bimodule from one exact ``dim_M × dim_M`` matrix per
        GENERATOR (each vertex ``e_v`` and each arrow, keyed by its basis label) per
        side; the actions of every other basis element are FOLDED via A's structure
        constants (``b_i = α · b_k`` ⇒ left is a map, right an anti-map). Floats are
        refused loudly; :meth:`check` certifies the result per instance. (P72
        ``H^*(B,M)``, P74 ``{}_gA``, generic no-code.) The canvas matrix EDITOR for
        this two-sided form is ledger-deferred to P80 (DD5); the library + server
        accept the explicit block here."""
        if dim_M > _MAX_COEFF_DIM:
            raise QuiverlabError(
                f"coefficient dimension {dim_M} exceeds the guard _MAX_COEFF_DIM="
                f"{_MAX_COEFF_DIM}", hint="split the computation or raise the guard")
        if A.quiver is None or not A.basis_labels:
            raise QuiverlabError(
                "from_actions needs a quiver presentation (generator-indexed action maps)",
                hint="present A via Quiver(...).algebra(...)")
        dom = A.domain
        m = A.dim
        labels = A.basis_labels
        # coerce ALL provided matrices first (floats refused regardless of shape)
        Lgiven = {k: _coerce_matrix(dom, v, dim_M) for k, v in left_maps.items()}
        Rgiven = {k: _coerce_matrix(dom, v, dim_M) for k, v in right_maps.items()}
        arrow_names = set(A.quiver.arrows)
        gen_index = {lab: i for i, lab in enumerate(labels)
                     if lab.startswith("e_") or lab in arrow_names}
        Lact = [None] * m
        Ract = [None] * m
        for lab, i in gen_index.items():
            if lab not in Lgiven or lab not in Rgiven:
                raise QuiverlabError(
                    f"from_actions is missing a left/right map for generator {lab!r}",
                    hint="supply one exact matrix per vertex and per arrow, both sides")
            Lact[i], Ract[i] = Lgiven[lab], Rgiven[lab]
        # fold to every remaining basis element via factorization b_i = arrow · b_k
        arrow_idx = [gen_index[a] for a in arrow_names if a in gen_index]
        remaining = [i for i in range(m) if Lact[i] is None]
        progress = True
        while remaining and progress:
            progress = False
            still = []
            for i in remaining:
                ei = A._basis_vec(i)
                done = False
                for ai in arrow_idx:
                    ea = A._basis_vec(ai)
                    for k in range(m):
                        if Lact[k] is None:
                            continue
                        if A.multiply(ea, A._basis_vec(k)) == ei:
                            # Lact_of(α·b_k) = Lact_of(b_k) @ Lact_of(α)  (left is a map)
                            # Ract_of(α·b_k) = Ract_of(α) @ Ract_of(b_k)  (right an anti-map)
                            Lact[i] = _mm(dom, Lact[k], Lact[ai])
                            Ract[i] = _mm(dom, Ract[ai], Ract[k])
                            done = True
                            break
                    if done:
                        break
                if done:
                    progress = True
                else:
                    still.append(i)
            remaining = still
        if remaining:
            raise QuiverlabError(
                "from_actions could not fold the action to every basis element "
                "(the basis is not generated by the vertices and arrows given)",
                hint="supply a matrix for every generator (vertex + arrow)")
        M = cls(A, dim_M, Lact, Ract, name=name)
        M.check()
        return M

    @classmethod
    def inflate(cls, L, B, coeff_over_B, *, vanishing_arrows):
        """The ``L``-bimodule obtained from a ``B``-bimodule ``coeff_over_B`` by
        **inflation along the split projection** ``π: L ↠ B`` (Plan 72, CMRS
        ``math/0102194``). ``L = B ⋉ M`` is a split square-zero extension: its
        generators split into those coming from ``B`` (the shared vertices ``e_v``
        and ``B``'s arrows) and the ``M``-part ``vanishing_arrows`` (the new arrows
        spanning the ideal ``M``, with ``M·M = 0``). For a ``B``-generator the
        ``L``-action copies ``coeff_over_B``'s action at that generator's ``B``-basis
        index; the ``vanishing_arrows`` act as **zero** (``M·M = 0``). The whole
        action is then folded to the full ``L``-basis by :meth:`from_actions` and
        certified by :meth:`check`.

        The ``M``-space (``dim_M = coeff_over_B.dim_M``) is unchanged; only the
        ambient algebra becomes ``L``. For the trivial-extension flagship
        (``coeff_over_B = Bimodule.dual(B)``, ``vanishing_arrows`` = the dual arrows)
        this returns ``D(B)`` viewed as the ideal of ``T(B)``; for
        ``coeff_over_B = Bimodule.regular(B)`` it returns ``B = L/M`` pulled back
        along ``π``. (Reused by P74 skew group algebras — noted for the record.)
        """
        if L.quiver is None or not L.basis_labels:
            raise QuiverlabError(
                "Bimodule.inflate needs a quiver presentation of L to read its "
                "generators", hint="present L via Quiver(...).algebra(...)")
        if B.basis_labels is None:
            raise QuiverlabError(
                "Bimodule.inflate needs B's basis labels to map generators of L "
                "back to B", hint="present B via Quiver(...).algebra(...)")
        dom = L.domain
        dim_M = coeff_over_B.dim_M
        b_index = {lab: i for i, lab in enumerate(B.basis_labels)}
        vanish = set(vanishing_arrows)
        zero = [[dom.zero()] * dim_M for _ in range(dim_M)]
        left_maps, right_maps = {}, {}
        gens = [f"e_{v}" for v in L.quiver.vertices] + list(L.quiver.arrows)
        for lab in gens:
            if lab in vanish:
                left_maps[lab] = [row[:] for row in zero]
                right_maps[lab] = [row[:] for row in zero]
                continue
            if lab not in b_index:
                raise QuiverlabError(
                    f"Bimodule.inflate: generator {lab!r} of L is neither a "
                    f"vanishing (M-part) arrow nor a generator of B",
                    hint="vanishing_arrows must be exactly the M-part arrows of L")
            j = b_index[lab]
            left_maps[lab] = [list(r) for r in coeff_over_B.Lact[j]]
            right_maps[lab] = [list(r) for r in coeff_over_B.Ract[j]]
        return cls.from_actions(L, dim_M, left_maps, right_maps,
                                name=coeff_over_B.name)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _identity_matrix(dom, m):
    return [[dom.one() if a == b else dom.zero() for b in range(m)] for a in range(m)]


def _coerce_matrix(dom, mat, dim_M):
    """Coerce a raw ``dim_M × dim_M`` matrix of exact entries into the Domain,
    refusing floats loudly (floats are rejected BEFORE the shape check)."""
    for row in mat:
        for x in row:
            reject_inexact(x)                 # floats refused regardless of shape
    if len(mat) != dim_M or any(len(row) != dim_M for row in mat):
        raise QuiverlabError(
            f"expected a {dim_M}x{dim_M} matrix", hint="one exact entry per cell")
    return [[dom.coerce(x) for x in row] for row in mat]


def _quotient_projection(A, socle):
    """Given a subspace ``socle`` (a list of coordinate vectors of A, a sub-bimodule),
    return ``(comp, project)`` where ``comp`` is a list of A-basis indices completing
    ``socle`` to a full basis and ``project(v)`` maps ``v`` to its coordinates in the
    quotient basis {images of e_i : i in comp} (dropping the socle part)."""
    dom = A.domain
    m = A.dim
    d = len(socle)
    # greedily pick standard basis vectors e_i extending socle to a basis of A
    base = [list(v) for v in socle]
    comp = []
    for i in range(m):
        ei = A._basis_vec(i)
        if rank(base + [ei], dom) > rank(base, dom):
            base.append(ei)
            comp.append(i)
        if len(comp) == m - d:
            break
    # change-of-basis: columns = socle vectors then complement e_i; Cmat rows.
    cols = [list(v) for v in socle] + [A._basis_vec(i) for i in comp]
    Cmat = [[cols[k][t] for k in range(m)] for t in range(m)]

    def project(v):
        x = solve(Cmat, list(v), dom)
        if x is None:
            raise QuiverlabError(
                "quotient projection failed: image left A (socle not a sub-bimodule?)",
                hint="internal invariant — please report this algebra")
        return x[d:]                          # drop the socle coordinates

    return comp, project


# ---------------------------------------------------------------------------
# Chain-level twisted (co)homology data (P54 cross-plan contract W2)
# ---------------------------------------------------------------------------
def _reps_mod_image(cycles, image, dom):
    """A subset of ``cycles`` linearly independent modulo span(image)."""
    reps, base = [], list(image)
    base_rank = rank(base, dom) if base else 0
    for v in cycles:
        rr = rank(base + [v], dom)
        if rr > base_rank:
            reps.append(v)
            base = base + [v]
            base_rank = rr
    return reps


def _columns(M):
    return [[row[c] for row in M] for c in range(len(M[0]))] if M and M[0] else []


def twisted_homology_classes(A, M, top, max_cells=4_000_000):
    """Chain-level twisted Hochschild HOMOLOGY data (the P54 BV cross-plan
    contract). Returns a list of per-degree dicts ``0..top`` with keys:

    - ``basis`` — the twisted bar-chain basis, the ``(s, J)`` enumeration of
      ``M ⊗ Ā^{⊗n}`` (in the unit-adapted basis the bar complex uses);
    - ``boundary`` — the exact Domain boundary matrix ``b_n : C_n -> C_{n-1}``
      (``None`` at ``n = 0``);
    - ``classes`` — a basis of homology-class representatives, a chosen lift of
      ``ker b_n / im b_{n+1}`` (each a coordinate vector in ``C_n``).

    CONVENTION BLOCK (cross-plan contract with P54, which will mirror it VERBATIM;
    a divergence between this block and P54 §3 is a review defect on whichever plan
    merges second). It states what THIS code does — ``bar.boundary_matrix``:

        *Twist placement:* the coefficient ``M`` occupies the leading tensor factor
        ``M ⊗ Ā^{⊗n}``. The boundary's outer two faces use ``M``'s actions: **term
        0 / the FIRST face** is the RIGHT action ``m·a_1`` (``+``), and **term n /
        the LAST face** is the LEFT action ``a_n·m`` (``(−1)^n``). For the Nakayama
        twist ``{}_1A_ν`` (``Lact`` untwisted, ``Ract`` = ``ψ``-twisted) the twist
        therefore lives in the RIGHT action and enters on the **FIRST face**
        (``m·a_1 = m·ψ(a_1)``); the last face (left action ``a_n·m``) is untwisted.
        *Normalization:* reduced bar (``Ā = A/k·1``), so index 0 is dropped from
        every bar slot. *Basis ordering:* ``M``-basis index ``s`` outermost
        (slowest), then the bar multi-index ``J ∈ {1..m-1}^n`` in lexicographic
        order (fastest last). *Cohomology twin:* ``twisted_cohomology_classes``
        mirrors with ``Hom_{A^e}`` collapse and the transposed convention (the
        twisted right action enters the coboundary's LAST face, ``f(a_1..a_n)·a_{n+1}``).
    """
    from quiverlab.hochschild import bar
    B = A.unit_adapted()
    Mb = bar._coeff_in_unit_basis(A, B, M)
    dom = B.domain
    m = B.dim
    dim_M = m if Mb is None else Mb.dim_M
    bmats = {}
    for n in range(0, top + 2):
        if n == 0:
            bmats[0] = None
            continue
        bmats[n], _, _ = bar.boundary_matrix(B, n, max_cells, coefficients=Mb)
    out = []
    for n in range(top + 1):
        basis = bar._cochain_basis(m, n, dim_M)
        bn = bmats[n]
        bnp1 = bmats[n + 1]
        # cycles = ker b_n (all of C_n when n == 0)
        if n == 0:
            cycles = [[dom.one() if i == j else dom.zero() for j in range(len(basis))]
                      for i in range(len(basis))]
        else:
            cycles = nullspace(bn, dom)
        image = _columns(bnp1) if bnp1 else []
        classes = _reps_mod_image(cycles, image, dom)
        out.append({"basis": basis, "boundary": bn, "classes": classes})
    return out


def twisted_cohomology_classes(A, M, top, max_cells=4_000_000):
    """Chain-level twisted Hochschild COHOMOLOGY data (the cohomology twin of
    :func:`twisted_homology_classes`; see its CONVENTION BLOCK). Returns per-degree
    dicts with ``basis``, ``coboundary`` (``δ^n : C^n -> C^{n+1}``), and
    ``classes`` = a basis of ``ker δ^n / im δ^{n-1}``."""
    from quiverlab.hochschild import bar
    B = A.unit_adapted()
    Mb = bar._coeff_in_unit_basis(A, B, M)
    dom = B.domain
    m = B.dim
    dim_M = m if Mb is None else Mb.dim_M
    dmats = {}
    for n in range(0, top + 1):
        dmats[n], _, _ = bar.coboundary_matrix(B, n, max_cells, coefficients=Mb)
    out = []
    for n in range(top + 1):
        basis = bar._cochain_basis(m, n, dim_M)
        dn = dmats[n]
        # cocycles = ker delta^n
        cycles = nullspace(dn, dom) if dn and dn[0] else \
            [[dom.one() if i == j else dom.zero() for j in range(len(basis))]
             for i in range(len(basis))]
        image = _columns(dmats[n - 1]) if n else []
        classes = _reps_mod_image(cycles, image, dom)
        out.append({"basis": basis, "coboundary": dn, "classes": classes})
    return out
