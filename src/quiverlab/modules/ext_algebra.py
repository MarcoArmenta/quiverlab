# SPDX-License-Identifier: MIT
"""The Yoneda / Ext-algebra E(A) = Ext^*_A(A/J, A/J) as a graded quiver-with-relations
presentation over the semisimple base R = k^{Q_0} (Plan 27).

Right modules, left-to-right composition (ASS). The locked convention (arbitrated by
explicit lift computations on k[x]/(x^2), kA_2 and kQ/J^2, mirrored by the Plan-05
module engine):

  * Ext^n(S_i, S_j) sits in corner (i, j): e_i E e_j = (+)_n E^n_{ij}.
  * The product is LEFT-TO-RIGHT: g.f for g in E^n_{ij}, f in E^m_{jk} lands in
    E^{n+m}_{ik} (composability target(g) = source(f)). Lift the FIRST factor g
    through the resolution of S_j, then apply f.
  * Hence the Ext-quiver equals Q: dim Ext^1(S_i,S_j) = #{arrows i->j},
    dim Ext^2(S_i,S_k) = #{minimal relations i->k}.

Two facts make the engine a pure exercise in graded Betti numbers + one chain-map
lift, with no coboundary arithmetic anywhere:

  1. Ext dims = graded Betti numbers.  On the Plan-05 MINIMAL resolution P_* -> S_i
     (modules/resolution.py::minimal_resolution), minimality d(P_{n+1}) subset rad P_n
     forces every differential of Hom(P_*, S_j) to vanish, so
     Ext^n(S_i,S_j) = Hom(P_n, S_j) = k^{b_{n,ij}}, b_{n,ij} = # of P_j summands in
     P_n = res.term(n).count(j). A basis cocycle is a summand projection
     P_n ->> P_j ->> top P_j = S_j, addressed as (source i, degree n, summand index).

  2. Yoneda product = chain-map lifting.  For a summand-projection cocycle
     g: P_n^{(i)} -> S_j lift it to a chain map g_t: P_{n+t}^{(i)} -> Q_t^{(j)} over the
     resolution Q_* -> S_j (g_0 covers g; d^Q_{t+1} g_{t+1} = g_t d^P_{n+t+1}, each an
     exact underdetermined solve canonicalised with reduce_mod_nullspace a la Plan 17
     for byte-reproducibility).  The class of g.f is read off f o g_m on the summand
     basis of E^{n+m}_{ik} -- lift-independent ON THE NOSE here (two lifts differ by
     d^Q h + h d^P; f kills the first term as a cocycle, minimality kills the second
     literally on generators, not merely up to coboundary).

Cocycles are represented by their address (source, degree, summand); products return
coordinate vectors over the summand basis of the target corner.  Generators (complement
of the decomposables per corner) and minimal relations (complement of the
lower-relation ideal inside ker of the free-cover multiplication) are exact linear
algebra over R = k^{Q_0} (composable corner pairs only).  A three-valued Koszulity
verdict is delegated to modules/koszul.py.
"""
from dataclasses import dataclass
from typing import Optional

from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.builders import _require_provenance, projective
from quiverlab.modules.ext import GlobalDimension
from quiverlab.modules.resolution import minimal_resolution


# --------------------------------------------------------------------------------------
# per-simple minimal resolution with retained summand structure
# --------------------------------------------------------------------------------------
class _SimpleResolution:
    """The minimal projective resolution of S_i (Plan 05) with the per-degree summand
    layout retained.  Each term P_n = (+)_s P_{v_s}; summand s occupies the basis block
    [off_s : off_s + dim P_{v_s}] and its generator (the trivial path e_{v_s}) is the
    first basis vector of the block, i.e. the column at off_s."""

    def __init__(self, A, i, length):
        self.A = A
        self.i = i
        self.terms, self.dmats = minimal_resolution(A.simple(i), length)
        self._info = [None] * len(self.terms)

    def dim(self, n):
        return self.terms[n].dim

    def module(self, n):
        return self.terms[n].module

    def has_degree(self, n):
        return 0 <= n < len(self.terms)

    def info(self, n):
        """List of (vertex, offset, dim, basis_labels) for the summands of P_n, in the
        direct-sum block order.  Rebuilt from projective(A, v) -- deterministic, so the
        offsets match the block layout _direct_sum produced."""
        if self._info[n] is None:
            out, off = [], 0
            for v in self.terms[n].vertices:
                Pv = projective(self.A, v)
                out.append((v, off, Pv.dim, Pv._pv_basis_labels))
                off += Pv.dim
            self._info[n] = out
        return self._info[n]

    def differential(self, n):
        """d_n: P_n -> P_{n-1} (dim P_{n-1} x dim P_n); column c = d_n(basis vector c)."""
        return self.dmats[n] if n < len(self.dmats) else []


def _build_module_map(src_res, src_deg, tgt_res, tgt_deg, ys, dom):
    """The full matrix (dim Q_{tgt_deg} x dim P_{src_deg}) of the module map
    P_{src_deg}^{(src)} -> Q_{tgt_deg}^{(tgt)} sending each summand generator to the
    prescribed image ys[s] in Q_{tgt_deg}.  A map out of a projective is determined by
    the generator images; the block basis vector (generator_s * path) maps to
    ys[s] * path = action[path] @ ys[s]."""
    info = src_res.info(src_deg)
    Q = tgt_res.module(tgt_deg)
    G = lm.zeros(tgt_res.dim(tgt_deg), src_res.dim(src_deg), dom)
    for (v, off, dv, labels), y in zip(info, ys):
        for k, plabel in enumerate(labels):
            col = lm.matvec(Q.action[plabel], y, dom)      # y * path_k in Q
            for r in range(len(G)):
                G[r][off + k] = col[r]
    return G


def _solve_in_component(dQ, dimQ, Qmod, v, rhs, dom):
    """Solve d^Q y = rhs for a generator image y living in the vertex-v component of Q
    (Q * e_v).  The idempotent action is diagonal on a direct sum of P_w's, so the
    component is the set of basis columns c with action[e_v][c][c] != 0.  The solve is
    canonicalised with reduce_mod_nullspace so the lift is byte-reproducible (Plan 17)."""
    ev = Qmod.action[f"e_{v}"]
    cols_v = [c for c in range(dimQ) if not dom.is_zero(ev[c][c])]
    if not cols_v:
        # empty component: rhs must be 0 (d^Q respects the grading and rhs sits in the
        # v-component of Q_t); the only image is 0.
        return [dom.zero()] * dimQ
    Dred = [[dQ[r][c] for c in cols_v] for r in range(len(dQ))]
    z = linalg.solve(Dred, rhs, dom)
    if z is None:
        raise QuiverlabError(
            "ext_algebra: the chain-map lift solve is inconsistent",
            hint="this should never happen on a genuine minimal resolution (the "
                 "obstruction is a boundary); report the presentation")
    z = linalg.reduce_mod_nullspace(z, Dred, dom)
    y = [dom.zero()] * dimQ
    for c, zc in zip(cols_v, z):
        y[c] = zc
    return y


def _lift(resG, n, s_star, resQ, m, dom):
    """Chain-map lift of the summand-projection cocycle g = (i, n, s_star): P_n^{(i)} ->
    S_j through the resolution Q_* -> S_j.  Returns g_m: P_{n+m}^{(i)} -> Q_m^{(j)}.

    g_0 covers g: it sends the distinguished summand's generator to the generator of
    Q_0 = P_j and every other generator to 0 (a valid lift, ep_j o g_0 = g).  Then
    g_{t+1} is the per-generator solve d^Q_{t+1} g_{t+1} = g_t d^P_{n+t+1}."""
    info0 = resG.info(n)
    dimQ0 = resQ.dim(0)
    ys = []
    for idx, (v, off, dv, labels) in enumerate(info0):
        y = [dom.zero()] * dimQ0
        if idx == s_star:
            y[0] = dom.one()                                # generator of Q_0 = P_j
        ys.append(y)
    g = _build_module_map(resG, n, resQ, 0, ys, dom)
    for t in range(m):
        src_deg = n + t + 1
        dP = resG.differential(src_deg)
        dQ = resQ.differential(t + 1)
        Qnext = resQ.module(t + 1)
        dimQnext = resQ.dim(t + 1)
        ys = []
        for (v, off, dv, labels) in resG.info(src_deg):
            gencol = [dP[r][off] for r in range(len(dP))]   # d^P(generator_s)
            rhs = lm.matvec(g, gencol, dom)                 # g_t(d^P(generator_s))
            ys.append(_solve_in_component(dQ, dimQnext, Qnext, v, rhs, dom))
        g = _build_module_map(resG, src_deg, resQ, t + 1, ys, dom)
    return g


# --------------------------------------------------------------------------------------
# records
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ExtGenerator:
    """A minimal generator of E over R: the cochain (source, degree, summand) is its
    canonical representative, `id` names it in relation words, and (source -> target)
    is its corner."""
    id: int
    degree: int
    source: object
    target: object
    cochain: tuple

    def __repr__(self):
        return f"gen#{self.id}(deg {self.degree}, {self.source}->{self.target})"


@dataclass(frozen=True)
class ExtRelation:
    """A minimal relation of E in a fixed corner/degree: an R-linear combination of
    composable generator-words `terms = ((coeff, (gen_id, ...)), ...)`."""
    degree: int
    source: object
    target: object
    terms: tuple

    def word_length(self):
        """The path length (number of generator letters) of the relation's words. A
        relation is 'quadratic' iff every word has length 2."""
        return sorted({len(word) for _, word in self.terms})

    def __repr__(self):
        return (f"rel(deg {self.degree}, {self.source}->{self.target}, "
                f"{len(self.terms)} terms)")


# --------------------------------------------------------------------------------------
# the engine
# --------------------------------------------------------------------------------------
class _ExtAlgebraEngine:
    """Builds the graded pieces of E: Betti numbers, the Yoneda product on basis
    cochains, and the minimal generators/relations over R = k^{Q_0}."""

    def __init__(self, A, reslen):
        self.A = A
        self.dom = A.domain
        self.verts = list(A.quiver.vertices)
        self.reslen = reslen
        self._res = {}
        self._prod = {}

    # -- resolutions & corner bases ----------------------------------------------------
    def resolution(self, i):
        r = self._res.get(i)
        if r is None:
            r = _SimpleResolution(self.A, i, self.reslen)
            self._res[i] = r
        return r

    def corner_basis(self, i, n, j):
        """Summand indices s of P_n^{(i)} whose vertex is j -- the basis of E^n_{ij}."""
        r = self.resolution(i)
        if not r.has_degree(n):
            return []
        return [s for s, (v, off, dv, lab) in enumerate(r.info(n)) if v == j]

    def ext_dim(self, i, n, j):
        return len(self.corner_basis(i, n, j))

    def target_of(self, cochain):
        i, n, s = cochain
        return self.resolution(i).info(n)[s][0]

    # -- the Yoneda product ------------------------------------------------------------
    def product(self, cg, cf):
        """g.f for basis cochains cg = (i, n, s) and cf = (j, m, s'): returns
        (coords, k) with coords over corner_basis(i, n+m, k), k = target(f).  Requires
        composability target(g) = source(f) = j."""
        key = (cg, cf)
        cached = self._prod.get(key)
        if cached is not None:
            return cached
        i, n, s = cg
        j, m, sp = cf
        resG = self.resolution(i)
        jg = resG.info(n)[s][0]
        if jg != j:
            raise QuiverlabError(
                "ext_algebra: internal product of non-composable cochains "
                f"(target {jg} != source {j})")
        resQ = self.resolution(j)
        kf = resQ.info(m)[sp][0]
        # degree-0 factors are the base idempotents R = k^{Q_0}: two-sided unit.
        if n == 0:                                          # g = e_i, product = f
            out = (self._unit(i, m, sp, kf), kf)
        elif m == 0:                                        # f = e_j, product = g
            out = (self._unit(i, n, s, kf), kf)
        elif not resG.has_degree(n + m) or not self.corner_basis(i, n + m, kf):
            out = ([], kf)                                  # lands in a zero corner
        else:
            gm = _lift(resG, n, s, resQ, m, self.dom)
            r = resQ.info(m)[sp][1]                         # generator row of f in Q_m
            infoNM = resG.info(n + m)
            out = ([gm[r][infoNM[t][1]] for t in self.corner_basis(i, n + m, kf)], kf)
        self._prod[key] = out
        return out

    def _unit(self, i, n, s, k):
        cb = self.corner_basis(i, n, k)
        return [self.dom.one() if t == s else self.dom.zero() for t in cb]

    def _mul_elem_by_gen(self, i, p, l, elem, gen):
        """Multiply an element of E^p_{il} (coords over corner_basis(i,p,l)) on the RIGHT
        by a basis generator gen in E^{q}_{l, l'}; returns coords over
        corner_basis(i, p+q, l')."""
        cb = self.corner_basis(i, p, l)
        lprime = gen.target
        out = [self.dom.zero()] * self.ext_dim(i, p + gen.degree, lprime)
        for pos, coeff in enumerate(elem):
            if self.dom.is_zero(coeff):
                continue
            coords, _ = self.product((i, p, cb[pos]), gen.cochain)
            for a, c in enumerate(coords):
                out[a] = self.dom.add(out[a], self.dom.mul(coeff, c))
        return out, lprime

    def phi(self, word, gen_by_id):
        """Image in E of a generator-word (tuple of generator ids), left-to-right; the
        free-cover multiplication.  Returns (coords, source, target, degree)."""
        gens = [gen_by_id[g] for g in word]
        g0 = gens[0]
        i = g0.source
        elem = self._unit(i, g0.degree, g0.cochain[2], g0.target)
        deg, cur = g0.degree, g0.target
        for gen in gens[1:]:
            elem, cur = self._mul_elem_by_gen(i, deg, cur, elem, gen)
            deg += gen.degree
        return elem, i, cur, deg

    # -- generators --------------------------------------------------------------------
    def decomposable_columns(self, i, n, j):
        """A spanning set (as columns over corner_basis(i,n,j)) of the decomposables
        D^n_{ij} = sum_{a+b=n, a,b>=1} im(E^a_{il} (x) E^b_{lj} -> E^n_{ij})."""
        cols = []
        for a in range(1, n):
            b = n - a
            for l in self.verts:
                Ba = self.corner_basis(i, a, l)
                Bb = self.corner_basis(l, b, j)
                if not Ba or not Bb:
                    continue
                for su in Ba:
                    for sw in Bb:
                        coords, _ = self.product((i, a, su), (l, b, sw))
                        cols.append(coords)
        return cols

    def new_generators(self, i, n, j, next_id):
        """The basis cochains of E^n_{ij} that are NOT decomposable -- one minimal
        generator each.  `next_id` is a mutable [int] counter."""
        cb = self.corner_basis(i, n, j)
        if not cb:
            return []
        dom = self.dom
        ident = lm.identity(len(cb), dom)                   # the standard basis of E^n_{ij}
        cands = [lm.col(ident, c) for c in range(len(cb))]
        keep = lm.independent_modulo(cands, self.decomposable_columns(i, n, j), dom)
        gens = []
        for pos in keep:
            gens.append(ExtGenerator(next_id[0], n, i, j, (i, n, cb[pos])))
            next_id[0] += 1
        return gens


# --------------------------------------------------------------------------------------
# relations over the free R-algebra on the accumulated generators
# --------------------------------------------------------------------------------------
class _RelationBuilder:
    def __init__(self, engine, generators):
        self.eng = engine
        self.dom = engine.dom
        self.gens = generators
        self.by_id = {g.id: g for g in generators}
        self._words = {}
        self.relations_by_degree = {}

    def words(self, n, i, j):
        """All composable generator-words (tuples of ids) of total degree n from i to
        j, length >= 1."""
        key = (n, i, j)
        cached = self._words.get(key)
        if cached is not None:
            return cached
        out = []
        for g in self.gens:
            if g.source != i:
                continue
            if g.degree == n and g.target == j:
                out.append((g.id,))
            elif g.degree < n:
                for rest in self.words(n - g.degree, g.target, j):
                    out.append((g.id,) + rest)
        self._words[key] = out
        return out

    def _words_or_empty(self, deg, i, j):
        if deg == 0:
            return [()] if i == j else []
        return self.words(deg, i, j)

    def compute_degree(self, n):
        rels = []
        for i in self.eng.verts:
            for j in self.eng.verts:
                rels.extend(self._corner(n, i, j))
        self.relations_by_degree[n] = rels
        return rels

    def _corner(self, n, i, j):
        words = self.words(n, i, j)
        if not words:
            return []
        dom = self.dom
        cdim = self.eng.ext_dim(i, n, j)
        # phi matrix: cdim rows, one column per word (its image in E^n_{ij}).
        phi_cols = [self.eng.phi(w, self.by_id)[0] for w in words]
        if cdim == 0:
            # target corner is 0-dimensional: every word maps to 0, so ker = all of F^n.
            ker = [lm.col(lm.identity(len(words), dom), c) for c in range(len(words))]
        else:
            phi_matrix = lm.cols_to_matrix(phi_cols)
            ker = lm.kernel_columns(phi_matrix, dom)
        if not ker:
            return []
        word_index = {w: k for k, w in enumerate(words)}
        lower = self._lower_ideal(n, i, j, words, word_index)
        keep = lm.independent_modulo(ker, lower, dom)
        rels = []
        for idx in keep:
            vec = ker[idx]
            terms = tuple((vec[k], words[k]) for k in range(len(words))
                          if not dom.is_zero(vec[k]))
            rels.append(ExtRelation(n, i, j, terms))
        return rels

    def _lower_ideal(self, n, i, j, words, word_index):
        """Columns (over `words`) spanning the degree-n part of the ideal generated by
        the already-found lower-degree relations: u . r . w for every lower relation r
        of degree p in corner (a,b) and words u: i->a, w: b->j with deg u + p + deg w = n."""
        dom = self.dom
        vecs = []
        for p in range(2, n):
            for rel in self.relations_by_degree.get(p, []):
                a, b = rel.source, rel.target
                for du in range(0, n - p + 1):
                    dw = n - p - du
                    for u in self._words_or_empty(du, i, a):
                        for wtail in self._words_or_empty(dw, b, j):
                            vec = [dom.zero()] * len(words)
                            for coeff, rword in rel.terms:
                                full = u + rword + wtail
                                k = word_index[full]
                                vec[k] = dom.add(vec[k], coeff)
                            vecs.append(vec)
        return vecs


# --------------------------------------------------------------------------------------
# the public result
# --------------------------------------------------------------------------------------
class YonedaPresentation:
    """A graded quiver-with-relations presentation of the Yoneda algebra
    E(A) = Ext^*_A(A/J, A/J) over R = k^{Q_0}, computed through a certified degree
    (Plan 27).  Honest by construction: when gl.dim A is finite (exact) the presentation
    is complete and self-hosting via `as_algebra()`; otherwise it names the truncation
    degree and never claims completeness."""

    def __init__(self, A, engine, degree_gen, degree_rel, gd,
                 generators_by_degree, relations_by_degree, ext_quiver):
        self.algebra = A
        self._eng = engine
        self._gd = gd
        self.certified_through_degree = degree_gen
        self.relation_degree = degree_rel
        self.generators_by_degree = generators_by_degree
        self.relations_by_degree = relations_by_degree
        self.ext_quiver = ext_quiver
        self.is_finite_dimensional = True if gd.exact else None
        # koszul verdict + the derived finiteness flags
        self.koszul, self.koszul_obstruction, self._koszul_reason = self._verdict()
        if gd.exact:
            self.is_finitely_generated_certified = True
        elif self.koszul is True:
            self.is_finitely_generated_certified = True     # Koszul => generated in deg 1
        else:
            self.is_finitely_generated_certified = None

    # -- graded dims -------------------------------------------------------------------
    def hilbert_matrix_through(self, d):
        """The Betti table b_{n,ij} = dim E^n_{ij} for n = 0..d, as
        `M[n][a][b]` indexed by vertex POSITION (a, b) in the fixed vertex order."""
        if d > self.certified_through_degree:
            if not self._gd.exact:
                raise QuiverlabError(
                    f"ext_algebra: dims certified only through degree "
                    f"{self.certified_through_degree}",
                    hint="raise `top`; gl.dim(A) is not finite so higher degrees are "
                         "genuinely open")
            # finite gl.dim: E^n = 0 above gl.dim, so pad honestly with zeros.
        verts = self._eng.verts
        out = []
        for n in range(d + 1):
            out.append([[self._eng.ext_dim(i, n, j) for j in verts] for i in verts])
        return out

    def graded_dims_through(self, d):
        """Total dim E^n for n = 0..d (sum over corners)."""
        verts = self._eng.verts
        return [sum(self._eng.ext_dim(i, n, j) for i in verts for j in verts)
                for n in range(d + 1)]

    # -- Koszulity ---------------------------------------------------------------------
    def _verdict(self):
        """Three-valued Koszul verdict, delegating the quadraticity / G-quadratic /
        Froberg machinery to modules/koszul.py.  If that sibling has not landed the
        verdict is left None with a clear reason (the tests skip accordingly)."""
        try:
            from quiverlab.modules.koszul import (
                froberg_obstruction, g_quadratic_certificate, is_quadratic)
        except ImportError:
            return None, None, "koszul module not available"
        A, d = self.algebra, self.certified_through_degree
        if g_quadratic_certificate(A):
            return True, None, "certified: G-quadratic (Priddy PBW)"
        if not is_quadratic(A):
            return False, (2, "non-quadratic: a defining relation is not length 2"), None
        for deg in range(2, self.certified_through_degree + 1):
            if self.generators_by_degree.get(deg):
                return (False,
                        (deg, f"a new Ext-algebra generator appears in degree {deg}"),
                        None)
        hilbert = self.hilbert_matrix_through(d)
        fob = froberg_obstruction(A, hilbert, d)
        if fob is not None:
            return (False,
                    (fob, f"Froberg identity P(t)C(-t)=I fails in degree {fob}"), None)
        return None, None, f"no obstruction through degree {d}"

    # -- self-hosting ------------------------------------------------------------------
    def as_algebra(self):
        """E(A) as a genuine quiverlab Algebra, via structure constants over its own
        graded basis -- available ONLY when gl.dim A is finite (exact), the case in
        which E is finite-dimensional and every product is known.  Raises with an honest
        truncation message otherwise.

        The table entries are already exact elements of A.domain (products through the
        Yoneda engine), so the algebra is built directly over that domain -- no float
        ever appears, and it works uniformly over GF(p) / CC / any exact Domain."""
        if self.is_finite_dimensional is not True:
            raise QuiverlabError(
                "ext_algebra.as_algebra: E(A) is not certified finite-dimensional -- "
                f"only computed through degree {self.certified_through_degree} "
                "(a truncation, not the whole algebra)",
                hint="as_algebra() needs gl.dim(A) < infinity (exact); the Ext-algebra "
                     "here is infinite-dimensional or its finiteness is unproven. "
                     "Inspect generators_by_degree / relations_by_degree instead")
        eng, dom = self._eng, self._eng.dom
        g = self._gd.value
        # flat graded basis: (i, n, s) for n = 0..g, degree 0 (the R idempotents) first
        basis, index, labels = [], {}, []
        for n in range(g + 1):
            for i in eng.verts:
                for j in eng.verts:
                    for s in eng.corner_basis(i, n, j):
                        index[(i, n, s)] = len(basis)
                        labels.append(f"E{n}_{i}_{j}_{s}")
                        basis.append((i, n, s))
        dimE = len(basis)
        z = [dom.zero()] * dimE
        T = [[list(z) for _ in range(dimE)] for _ in range(dimE)]
        for a, ca in enumerate(basis):
            ia, na, sa = ca
            ja = eng.target_of(ca)
            for b, cb in enumerate(basis):
                ib, nb, sb = cb
                if ja != ib or na + nb > g:                  # non-composable / above top -> 0
                    continue
                coords, kb = eng.product(ca, cb)
                col = list(z)
                for pos, s in enumerate(eng.corner_basis(ia, na + nb, kb)):
                    col[index[(ia, na + nb, s)]] = coords[pos]
                T[a][b] = col
        unit = list(z)
        for i in eng.verts:
            for s in eng.corner_basis(i, 0, i):
                unit[index[(i, 0, s)]] = dom.one()
        from quiverlab.core.algebra import Algebra
        # entries already live in A.domain; build directly and validate (assoc + unit).
        E = Algebra(dom, T, unit, basis_labels=labels)
        E._validate()
        if E.dim != sum(self.graded_dims_through(g)):
            raise QuiverlabError(
                "ext_algebra.as_algebra: built dimension disagrees with the graded "
                "dimensions -- internal inconsistency")
        return E

    def __repr__(self):
        ngen = sum(len(v) for v in self.generators_by_degree.values())
        nrel = sum(len(v) for v in self.relations_by_degree.values())
        tail = "finite-dim" if self.is_finite_dimensional else \
            f"truncated at degree {self.certified_through_degree}"
        return (f"YonedaPresentation(E({self.algebra.__class__.__name__}), "
                f"{ngen} generators, {nrel} relations, {tail}, koszul={self.koszul})")


# --------------------------------------------------------------------------------------
# the factory
# --------------------------------------------------------------------------------------
def _finiteness(A, bound):
    """gl.dim A as a GlobalDimension (value, exact) -- the modules/ext.py::global_dimension
    algorithm (sup over simples of pd S_v), but with a modest resolution depth `bound` and
    growth cap so an algebra with exponentially growing minimal terms (e.g. rad^2=0 on
    several loops -- gl.dim infinite) is reported infinite promptly instead of blowing up.
    A syzygy that overshoots the cap is treated exactly like a simple that fails to resolve
    within `bound`: gl.dim is not finite (exact=False)."""
    best, exact = 0, True
    for v in A.quiver.vertices:
        try:
            res = A.simple(v).projective_resolution(bound, max_term_dim=20000)
        except DepthLimitError:
            exact = False
            best = max(best, bound)
            continue
        pd = res.pd()
        if pd is None:
            exact = False
            best = max(best, bound)
        else:
            best = max(best, pd)
    return GlobalDimension(best, exact)


def ext_algebra(A, top=6):
    """Compute the Yoneda algebra presentation E(A) through degree `top` (or, when
    gl.dim A is finite and exact, complete through gl.dim).  Returns a
    YonedaPresentation (Plan 27)."""
    _require_provenance(A, "ext_algebra")
    if A.quiver is None or A.relations is None:
        raise QuiverlabError(
            "ext_algebra needs the quiver presentation",
            hint="build the algebra via Quiver.algebra(...); structure-constant "
                 "algebras carry no path basis")
    gd = _finiteness(A, bound=max(top, 6))
    if gd.exact:
        degree_gen = gd.value                 # E^{>gl.dim} = 0: generators stop here
        degree_rel = 2 * gd.value             # minimal relations can reach twice that
        reslen = max(gd.value, 1)
    else:
        degree_gen = degree_rel = top
        reslen = top
    eng = _ExtAlgebraEngine(A, reslen)

    # generators degree by degree (deg 0 = the base R, not a generator of the aug ideal)
    generators_by_degree = {}
    all_gens = []
    next_id = [0]
    for n in range(1, degree_gen + 1):
        gens_n = []
        for i in eng.verts:
            for j in eng.verts:
                gens_n.extend(eng.new_generators(i, n, j, next_id))
        if gens_n:
            generators_by_degree[n] = gens_n
            all_gens.extend(gens_n)

    # minimal relations degree by degree over the free R-algebra on the generators
    rb = _RelationBuilder(eng, all_gens)
    relations_by_degree = {}
    for n in range(2, degree_rel + 1):
        rels_n = rb.compute_degree(n)
        if rels_n:
            relations_by_degree[n] = rels_n

    # the Ext-quiver: same vertices as Q, one arrow per E^1 basis element (corner-tagged)
    ext_quiver = _build_ext_quiver(eng, generators_by_degree.get(1, []))

    return YonedaPresentation(A, eng, degree_gen, degree_rel, gd,
                              generators_by_degree, relations_by_degree, ext_quiver)


def _build_ext_quiver(eng, deg1_generators):
    """Quiver with vertices Q_0 and one arrow i->j per basis element of E^1_{ij}
    (= per minimal degree-1 generator).  The Ext-quiver equals Q for admissible A."""
    from quiverlab.combinat.quiver import Quiver
    arrows = {}
    for k, g in enumerate(deg1_generators):
        arrows[f"y{k}"] = (g.source, g.target)
    return Quiver(list(eng.verts), arrows)
