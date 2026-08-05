"""Recollements from an idempotent (Plan 47). For a vertex subset ``S`` with
``e = sum_{v in S} e_v``:

* the corner algebra ``eAe`` -- NOT the subquiver algebra: a path may leave ``S`` in its
  interior yet satisfy ``e*p*e = p != 0`` (``kA3`` with ``S={1,3}``: the path ``ab: 1->3``
  travels through vertex ``2 not in S`` but ``e(ab)e = ab``, so ``eAe = kA2`` of dim 3, NOT
  the subquiver ``k x k`` of dim 2). Built on the corner path-type basis (the ``A``-basis
  paths with BOTH endpoints in ``S``), presented as ``kQ'/I'`` and certified
  ``dim = sum_{v,w in S} dim e_v A e_w``.
* the quotient ``A/AeA``, presentable on the complement full subquiver (the ``A``-basis
  paths AVOIDING ``S``), certified ``dim = dim A - dim A e_S A`` (interior-aware).
* the six functors ``i_*, i^*, i^!, j^*, j_!, j_*`` of the recollement
  ``(mod A/AeA, mod A, mod eAe)`` as bimodule operations.

Char-clean exact linear algebra (a ``GF(2)`` cell proves it): both algebras are built on
subset path-type bases, never through the char-scoped trace-form radical. The self-certs --
the four adjunction dim identities, the two canonical exact sequences at each joint, and
``j^* j_! ~= id`` / ``j^* j_* ~= id`` -- ARE the test battery (QPA has NO recollement
surface, so the oracle class here is theory pins + internal self-certificates).

Cline-Parshall-Scott (J. Reine Angew. Math. 391 (1988), 85-99); recollement:
Beilinson-Bernstein-Deligne (Faisceaux pervers, Asterisque 100, 1982).
"""
from __future__ import annotations

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.families._present import present_from_pi
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules import radtopsoc
from quiverlab.modules.builders import (_label_vertex_source, _label_vertex_target,
                                        _require_provenance, projective)
from quiverlab.modules.hom import hom_space
from quiverlab.modules.module import Module
from quiverlab.modules.morphism import direct_sum


def _same_span(cols1, cols2, dom):
    """True iff the two column-sets span the same subspace (each rank == the union rank)."""
    r1 = lm.mat_rank(lm.cols_to_matrix(cols1), dom) if cols1 else 0
    r2 = lm.mat_rank(lm.cols_to_matrix(cols2), dom) if cols2 else 0
    both = (cols1 or []) + (cols2 or [])
    ru = lm.mat_rank(lm.cols_to_matrix(both), dom) if both else 0
    return r1 == r2 == ru


def _reduce(cols, dom):
    if not cols:
        return []
    piv = lm.column_space_pivots(lm.cols_to_matrix(cols), dom)
    return [cols[j] for j in piv]


class Recollement:
    """The recollement ``(mod A/AeA, mod A, mod eAe)`` of a presented algebra ``A`` at the
    idempotent ``e = sum_{v in S} e_v`` (``S`` a vertex subset). Builds ``eAe`` (corner,
    certified dim) and ``quotient = A/AeA`` (complement subquiver, certified dim) plus the
    six functors as methods. Char-clean."""

    def __init__(self, A, S):
        _require_provenance(A, "Recollement")
        self.A = A
        self.dom = A.domain
        self.S = [v for v in A.quiver.vertices if v in set(S)]   # keep A's vertex order
        self.Sset = set(S)
        if not self.Sset.issubset(set(A.quiver.vertices)):
            raise QuiverlabError(
                f"Recollement: S={list(S)!r} has vertices not in the quiver "
                f"{list(A.quiver.vertices)!r}")
        # basis-label endpoint vertices + visit classification (path-type, char-clean)
        self._srcv = [_label_vertex_source(A, lab) for lab in A.basis_labels]
        self._tgtv = [_label_vertex_target(A, lab) for lab in A.basis_labels]
        self._corner_idx = [k for k in range(A.dim)
                            if self._srcv[k] in self.Sset and self._tgtv[k] in self.Sset]
        self._visit_idx = [k for k in range(A.dim) if self._visits(A.basis_labels[k])]
        self._avoid_idx = [k for k in range(A.dim) if k not in set(self._visit_idx)]
        self._idem_vec = {}                        # vertex -> A-idempotent vector
        for k, lab in enumerate(A.basis_labels):
            if lab.startswith("e_"):
                self._idem_vec[self._srcv[k]] = A._basis_vec(k)
        self.eAe, self._eae_iota = self._build_corner()
        self.quotient, self._q_iota = self._build_quotient()
        # A.dim x quotient.dim embedding matrix (columns iota_B(label)), for A ->> B projection
        cols = [self._q_iota[lab] for lab in self.quotient.basis_labels]
        self._q_iota_mat = lm.cols_to_matrix(cols) if cols else lm.zeros(A.dim, 0, self.dom)

    # -- classification helpers ------------------------------------------------ #
    def _visits(self, label):
        """True iff the basis path ``label`` visits some vertex in ``S`` (endpoints OR
        interior) -- this is exactly membership in ``A e_S A`` (interior-aware)."""
        A = self.A
        if label.startswith("e_"):
            return _label_vertex_source(A, label) in self.Sset
        word = label.split("*")
        verts = [A.quiver.source(word[0])] + [A.quiver.target(a) for a in word]
        return any(v in self.Sset for v in verts)

    # -- corner algebra eAe ---------------------------------------------------- #
    def _build_corner(self):
        A, dom = self.A, self.dom
        rad_idx = [k for k in self._corner_idx if not A.basis_labels[k].startswith("e_")]
        rad_vecs = [A._basis_vec(k) for k in rad_idx]
        prods = [A.multiply(r1, r2) for r1 in rad_vecs for r2 in rad_vecs]
        rad2 = _reduce(prods, dom)
        keep = lm.independent_modulo(rad_vecs, rad2, dom)         # rad/rad^2 = the arrows
        arrows, img = {}, {}
        for n, ki in enumerate(keep, 1):
            gi = rad_idx[ki]
            name = f"c{n}"
            arrows[name] = (self._srcv[gi], self._tgtv[gi])
            img[name] = rad_vecs[ki]
        Q = Quiver(list(self.S), arrows)
        dim_expected = len(self._corner_idx)
        base_bound = len(rad_idx) + 2
        eAe = present_from_pi(Q, img, A, dom, dim_expected, base_bound, citations=("cps",))
        if eAe.dim != dim_expected:                              # present_from_pi already gates
            raise QuiverlabError(
                f"Recollement: eAe dim {eAe.dim} != corner sum {dim_expected}")
        iota = self._compose_iota(eAe, img)
        return eAe, iota

    # -- quotient A/AeA -------------------------------------------------------- #
    def _build_quotient(self):
        A, dom = self.A, self.dom
        comp_verts = [v for v in A.quiver.vertices if v not in self.Sset]
        comp_arrows, img = {}, {}
        for a, (s, t) in A.quiver.arrows.items():
            if s not in self.Sset and t not in self.Sset:
                comp_arrows[a] = (s, t)
                img[a] = A._basis_vec(A.basis_labels.index(a))
        Q = Quiver(comp_verts, comp_arrows)
        dim_expected = len(self._avoid_idx)
        base_bound = dim_expected + 2
        B = present_from_pi(Q, img, A, dom, dim_expected, base_bound, citations=("cps",))
        iota = self._compose_iota(B, img)
        return B, iota

    def _compose_iota(self, alg, arrow_img):
        """Embedding ``alg -> A``: each ``alg`` basis label -> its ``A``-vector (idempotent
        or product of arrow images along the word)."""
        A = self.A
        iota = {}
        for lab in alg.basis_labels:
            if lab.startswith("e_"):
                v = _label_vertex_source(alg, lab)
                iota[lab] = list(self._idem_vec[v])
            else:
                acc = None
                for a in lab.split("*"):
                    acc = list(arrow_img[a]) if acc is None else A.multiply(acc, arrow_img[a])
                iota[lab] = acc
        return iota

    # -- public dim certificates (asserted in tests) --------------------------- #
    def _corner_dim_sum(self):
        """sum_{v,w in S} dim e_v A e_w (the corner certificate target)."""
        return len(self._corner_idx)

    def _aeA_dim(self):
        """dim A e_S A = number of A-basis paths visiting S (interior-aware)."""
        return len(self._visit_idx)

    # -- action utilities ------------------------------------------------------ #
    def _A_action(self, M, xvec):
        """The right action of the ``A``-element ``xvec`` (A-coords) on the ``A``-module
        ``M``: ``sum_k xvec[k] * M.action[label_k]`` (an ``M.dim x M.dim`` matrix)."""
        dom = self.dom
        n = M.dim
        out = lm.zeros(n, n, dom)
        for k, ck in enumerate(xvec):
            if dom.is_zero(ck):
                continue
            Mk = M.action[self.A.basis_labels[k]]
            for i in range(n):
                oi, mi = out[i], Mk[i]
                for j in range(n):
                    oi[j] = dom.add(oi[j], dom.mul(ck, mi[j]))
        return out

    @staticmethod
    def _alg_action(mod, alg, avec):
        """The action on ``mod`` of the ``alg``-element ``avec`` (alg-coords):
        ``sum_t avec[t] * mod.action[alg.basis_labels[t]]``."""
        dom = mod.domain
        n = mod.dim
        out = lm.zeros(n, n, dom)
        for t, ct in enumerate(avec):
            if dom.is_zero(ct):
                continue
            Mt = mod.action[alg.basis_labels[t]]
            for i in range(n):
                oi, mi = out[i], Mt[i]
                for j in range(n):
                    oi[j] = dom.add(oi[j], dom.mul(ct, mi[j]))
        return out

    def _proj_A_to_B(self, avec):
        """The image of the ``A``-vector ``avec`` in ``B = A/AeA``, as a ``B``-vector: drop
        the visiting-``S`` coordinates (that is ``A e_S A``) and solve against the ``B``
        embedding. ``None``-safe: a fully-visiting element maps to the zero ``B``-vector."""
        dom = self.dom
        masked = [avec[k] if k in set(self._avoid_idx) else dom.zero()
                  for k in range(self.A.dim)]
        B = self.quotient
        if B.dim == 0:
            return []
        if all(dom.is_zero(x) for x in masked):
            return [dom.zero()] * B.dim
        sol = lm.solve_columns(self._q_iota_mat, lm.cols_to_matrix([masked]), dom)
        if sol is None:
            raise QuiverlabError("Recollement: A ->> A/AeA projection failed (bug)")
        return sol[0]

    # -- j-side (corner eAe) --------------------------------------------------- #
    def j_upper_star(self, M):
        """``j^*(M) = M e_S`` with the ``eAe``-action (the corner restriction). Underlying
        space ``sum_{v in S} M e_v``; an ``eAe`` basis element acts by the ``A``-action of
        its embedded ``A``-element restricted to ``M e_S``."""
        A, dom = self.A, self.dom
        eAe = self.eAe
        cols = []
        for v in self.S:
            Ev = M.action[f"e_{v}"]
            cols.extend(lm.col(Ev, j) for j in range(M.dim))
        me = _reduce(cols, dom)
        r = len(me)
        Bm = lm.cols_to_matrix(me) if me else lm.zeros(M.dim, 0, dom)
        action = {}
        for lab in eAe.basis_labels:
            Mx = self._A_action(M, self._eae_iota[lab])
            outcols = []
            for c in me:
                img = lm.matvec(Mx, c, dom)
                sol = lm.solve_columns(Bm, lm.cols_to_matrix([img]), dom)
                if sol is None:
                    raise QuiverlabError("Recollement.j^*: M e_S not eAe-stable (bug)")
                outcols.append(sol[0])
            action[lab] = lm.cols_to_matrix(outcols) if outcols else lm.zeros(r, r, dom)
        return Module(eAe, r, action, name=f"j^*({M.name})", side="right")

    def j_shriek(self, X):
        """``j_!(X) = X (x)_{eAe} eA``, a right ``A``-module (induction, left adjoint of
        ``j^*``). Built as the quotient of the ambient ``X (x)_k eA`` by the balancing
        ``A``-submodule ``<(x.c) (x) p - x (x) (c.p)>``. Self-cert: ``j^* j_! X ~= X``."""
        A, dom = self.A, self.dom
        eAe = self.eAe
        Pmods = [projective(A, v) for v in self.S]
        eA_labels = [lab for P in Pmods for lab in P._pv_basis_labels]
        eA_idx = [A.basis_labels.index(lab) for lab in eA_labels]
        eA_avecs = [A._basis_vec(k) for k in eA_idx]
        eA_dim = len(eA_idx)
        dX = X.dim
        if dX == 0 or eA_dim == 0:
            return self._zero_A_module(f"j_!({X.name})")
        ambient, _inc, _prj = direct_sum(*([eA_module(A, self.S)] * dX))
        Wcols = []
        for t, clab in enumerate(eAe.basis_labels):
            XA = X.action[clab]                     # dX x dX : right eAe action of c
            ce = self._eae_iota[clab]
            for j in range(eA_dim):
                cp = A.multiply(ce, eA_avecs[j])    # left eAe action c.p_j (A-vector)
                cp_ea = [cp[eA_idx[b]] for b in range(eA_dim)]
                for i in range(dX):
                    w = [dom.zero()] * (dX * eA_dim)
                    for a in range(dX):             # (x_i . c) (x) p_j
                        w[a * eA_dim + j] = XA[a][i]
                    for b in range(eA_dim):         # - x_i (x) (c . p_j)
                        w[i * eA_dim + b] = dom.sub(w[i * eA_dim + b], cp_ea[b])
                    Wcols.append(w)
        Wcols = _reduce(Wcols, dom)
        T = radtopsoc.quotient(ambient, Wcols, name=f"j_!({X.name})", side="right")
        return T

    def j_star(self, X):
        """``j_*(X) = Hom_{eAe}(Ae, X)``, a right ``A``-module (coinduction, right adjoint of
        ``j^*``). ``Ae`` is the right ``eAe``-module of ``A``-paths ending in ``S``; the
        residual LEFT ``A``-action on ``Ae`` gives the right ``A``-action on the Hom-space.
        Self-cert: ``j^* j_* X ~= X``."""
        A, dom = self.A, self.dom
        eAe = self.eAe
        ae_idx = [k for k in range(A.dim) if self._tgtv[k] in self.Sset]
        ae_avecs = [A._basis_vec(k) for k in ae_idx]
        d_ae = len(ae_idx)
        if d_ae == 0 or X.dim == 0:
            return self._zero_A_module(f"j_*({X.name})")
        # Ae as a right eAe-module: (m . c) = m * iota(c), expressed in the ae basis.
        ae_action = {}
        for lab in eAe.basis_labels:
            c = self._eae_iota[lab]
            cols = []
            for k in range(d_ae):
                prod = A.multiply(ae_avecs[k], c)
                cols.append([prod[ae_idx[b]] for b in range(d_ae)])
            ae_action[lab] = lm.cols_to_matrix(cols)
        AeMod = Module(eAe, d_ae, ae_action, name="Ae", side="right")
        H = hom_space(AeMod, X)                      # X.dim x d_ae matrices
        m = len(H)
        if m == 0:
            return self._zero_A_module(f"j_*({X.name})")

        def flat(mat):
            return [x for row in mat for x in row]

        basisH = lm.cols_to_matrix([flat(h) for h in H])
        action = {}
        for lab_A in A.basis_labels:
            bvec = A._basis_vec(A.basis_labels.index(lab_A))
            Lb = lm.cols_to_matrix(  # column k = (b . m_k) in ae coords (left A-mult)
                [[A.multiply(bvec, ae_avecs[k])[ae_idx[b]] for b in range(d_ae)]
                 for k in range(d_ae)])
            cols = []
            for h in H:
                comp = lm.matmul(h, Lb, dom)         # f . b = f o Lb : X.dim x d_ae
                sol = lm.solve_columns(basisH, lm.cols_to_matrix([flat(comp)]), dom)
                if sol is None:
                    raise QuiverlabError("Recollement.j_*: f.b left Hom_{eAe}(Ae,X) (bug)")
                cols.append(sol[0])
            action[lab_A] = lm.cols_to_matrix(cols)
        return Module(A, m, action, name=f"j_*({X.name})", side="right")

    # -- i-side (quotient A/AeA) ---------------------------------------------- #
    def i_star(self, N):
        """An ``A/AeA``-module ``N`` as an ``A``-module via ``A ->> A/AeA`` (inflation,
        exact fully faithful): the ``A``-label ``b`` acts as the ``B``-action of its image
        in ``A/AeA`` (any path visiting ``S`` acts as ``0``)."""
        A, dom = self.A, self.dom
        n = N.dim
        action = {}
        for k, lab in enumerate(A.basis_labels):
            bvec = self._proj_A_to_B(A._basis_vec(k))
            action[lab] = (self._alg_action(N, self.quotient, bvec) if bvec
                           else lm.zeros(n, n, dom))
        return Module(A, n, action, name=f"i_*({N.name})", side="right")

    def i_upper_star(self, M):
        """``i^*(M) = M / M(AeA) = M (x)_A A/AeA`` (left adjoint of ``i_*``): quotient ``M``
        by the submodule ``M . (A e_S A)``; viewed as an ``A/AeA``-module."""
        A, dom = self.A, self.dom
        subcols = []
        for k in self._visit_idx:
            Ak = M.action[A.basis_labels[k]]
            subcols.extend(lm.col(Ak, j) for j in range(M.dim))
        subcols = _reduce(subcols, dom)
        Q = radtopsoc.quotient(M, subcols, name=f"i^*({M.name})", side="right")
        return self._restrict_to_B(Q, name=f"i^*({M.name})")

    def i_upper_shriek(self, M):
        """``i^!(M) = Hom_A(A/AeA, M)`` = the largest submodule of ``M`` annihilated by
        ``A e_S A`` (right adjoint of ``i_*``): ``intersection_g ker(M.action[g])``, viewed
        as an ``A/AeA``-module."""
        A, dom = self.A, self.dom
        inter = None
        for k in self._visit_idx:
            ker = lm.kernel_columns(M.action[A.basis_labels[k]], dom)
            inter = ker if inter is None else radtopsoc._intersect(inter, ker, dom)
        if inter is None:                            # no AeA -> whole M annihilated
            inter = [lm.col(lm.identity(M.dim, dom), j) for j in range(M.dim)]
        sub = radtopsoc.submodule(M, inter, name=f"i^!({M.name})", side="right")
        return self._restrict_to_B(sub, name=f"i^!({M.name})")

    def _restrict_to_B(self, Q, name="res_B"):
        """View the ``A``-module ``Q`` (annihilated by ``A e_S A``) as a ``B = A/AeA``-module,
        same underlying space: ``action_B[l] = `` the ``A``-action of ``iota_B(l)`` on Q."""
        B = self.quotient
        action = {lab: self._A_action(Q, self._q_iota[lab]) for lab in B.basis_labels}
        return Module(B, Q.dim, action, name=name, side="right")

    def _zero_A_module(self, name):
        A, dom = self.A, self.dom
        action = {lab: lm.zeros(0, 0, dom) for lab in A.basis_labels}
        return Module(A, 0, action, name=name, side="right")

    # -- certificate helpers (used by the tests) ------------------------------ #
    def _M_AeA_cols(self, M):
        """Columns spanning ``M . (A e_S A) = ker(unit M ->> i_* i^* M)``."""
        dom = self.dom
        cols = []
        for k in self._visit_idx:
            Ak = M.action[self.A.basis_labels[k]]
            cols.extend(lm.col(Ak, j) for j in range(M.dim))
        return _reduce(cols, dom)

    def _MeA_cols(self, M):
        """Columns spanning ``M e A = im(counit j_! j^* M -> M)`` (the A-submodule generated
        by ``M e_S``)."""
        A, dom = self.A, self.dom
        me = []
        for v in self.S:
            Ev = M.action[f"e_{v}"]
            me.extend(lm.col(Ev, j) for j in range(M.dim))
        me = _reduce(me, dom)
        gens = []
        for lab in A.basis_labels:
            Ab = M.action[lab]
            for c in me:
                gens.append(lm.matvec(Ab, c, dom))
        return _reduce(gens, dom)

    def _counit_sequence_exact(self, M):
        """``j_! j^* M --> M --> i_* i^* M --> 0`` exact at each joint: the second map is the
        (surjective) quotient by ``M(AeA)``, and ``im(counit) = M e A = M(AeA) = ker`` -- a
        span identity (both submodules of M coincide)."""
        return _same_span(self._MeA_cols(M), self._M_AeA_cols(M), self.dom)


def eA_module(A, S):
    """``eA = (+)_{v in S} e_v A`` as a right ``A``-module (a projective)."""
    Pmods = [projective(A, v) for v in A.quiver.vertices if v in set(S)]
    if not Pmods:
        from quiverlab.modules.duality import _zero_module
        return _zero_module(A, side="right")
    D, _inc, _prj = direct_sum(*Pmods)
    return D
