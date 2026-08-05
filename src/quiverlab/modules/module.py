"""The right A-module object and its radical/top/socle (spec §3.6).

RIGHT modules, anti-homomorphism convention (see the package docstring): an element m
is a COLUMN vector in a fixed k-basis of M; the action of an algebra basis element b is
the matrix action[b] with m*b = action[b] @ m, and action[x*y] = action[y] @ action[x].
The vertex subspace M*e_v is the image of action['e_v']; dimension_vector[v] = its rank.
"""
from fractions import Fraction

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm


_SIDES = ("right", "left")


def _other_side(side):
    """The categorical opposite side ("right" <-> "left"), used by the contravariant
    functors D and Tr, which exchange the two sides over the same base algebra."""
    return "left" if side == "right" else "right"


class Module:
    """A finite-dimensional A-module over an exact Domain (Plan 05; Plan 24 sides).

    The action is ALWAYS stored as a RIGHT action of ``self.algebra`` (the
    *representation* algebra). ``self.side`` records how the user reads it:
    ``"right"`` (default) means a right ``self.algebra``-module; ``"left"`` means a
    left ``base_algebra``-module, represented as a right ``self.algebra =
    base_algebra^op``-module (a left A-module IS a right A^op-module). All algorithms
    read only ``(self.algebra, self.action)`` and are blind to ``side`` -- the left
    side reuses every right-module algorithm run over A^op, no duplicated math."""

    def __init__(self, algebra, dim, action, name="M", side="right"):
        if side not in _SIDES:
            raise QuiverlabError(f"Module side must be one of {_SIDES}, got {side!r}")
        self.algebra = algebra
        self.domain = algebra.domain
        self.dim = dim
        self.side = side
        # Action entries may be given as plain int / rational literals for
        # convenience; normalize them into the algebra's exact field (mirroring
        # Algebra.from_structure_constants). Entries already in the field (e.g.
        # sympy MPQ over CC) are used as-is.
        self.action = {label: _coerce_matrix(mat, self.domain)
                       for label, mat in action.items()}
        self.name = name

    @property
    def base_algebra(self):
        """The algebra the user reads this module as a (side)-module over: the
        representation algebra when right, its opposite when left."""
        return self.algebra if self.side == "right" else self.algebra.opposite()

    def with_side(self, side):
        """A re-tagged twin over the SAME representation (same algebra + action),
        read on the given side. This is the exact "side translation" between a right
        A^op-module and a left A-module -- a relabelling, no recomputation."""
        import copy
        twin = copy.copy(self)
        twin.side = side
        twin.action = dict(self.action)      # shallow copy so mutations don't alias
        return twin

    def _idem_label(self, v):
        return f"e_{v}"

    def vertex_projection(self, v):
        return self.action[self._idem_label(v)]

    def dimension_vector(self):
        dom = self.domain
        out = {}
        for v in self.algebra.quiver.vertices:
            out[v] = lm.mat_rank(self.vertex_projection(v), dom)
        return out

    def _arrow_labels(self):
        return list(self.algebra.quiver.arrows)

    def check_module(self, extra_labels=None):
        """Verify the representation is a genuine right A-module and return
        (True, None) or (False, witness). The witness names exactly what failed.

        Checks, in order:
        (i)   the vertex idempotents act as orthogonal projections summing to the
              identity: sum_v P_v == I, P_v @ P_v == P_v, and P_v @ P_w == 0 (v != w);
        (ii)  every stored arrow action a: s(a) -> t(a) is compatible with the vertex
              grading, i.e. action[a] == P_{t(a)} @ action[a] @ P_{s(a)}. This is the
              exact identity forced on a RIGHT (anti-homomorphism) action by
              a = e_{s(a)} * a * e_{t(a)} in A: with m*b = action[b] @ m and
              action[x*y] = action[y] @ action[x] (see _action_of_word), the left
              idempotent e_{s(a)} multiplies on the RIGHT of action[a] and the right
              idempotent e_{t(a)} on the LEFT;
        (iii) every relation sum_w c_w * word is satisfied: sum_w c_w * action[word] == 0;
        (iv)  the action is multiplicative on the composite basis-label products it
              stores: for each stored label word = a1*...*ak, action[word] equals the
              anti-homomorphism composition action[ak] @ ... @ action[a1].

        (Multiplicativity is checked only for the composite path labels actually stored
        in `action`; a module carrying only idempotents and single arrows has no such
        product to check -- this is the honest scope of the stored data, not the full
        multiplication table.)"""
        dom = self.domain
        n = self.dim
        verts = list(self.algebra.quiver.vertices)
        projs = {v: self.vertex_projection(v) for v in verts}
        I = lm.identity(n, dom)
        # (i) idempotents: orthogonal projections summing to the identity
        acc = lm.zeros(n, n, dom)
        for v in verts:
            acc = _add(acc, projs[v], dom)
        if not _mat_eq(acc, I, dom):
            return False, "sum of e_v actions != identity"
        for v in verts:
            Pv = projs[v]
            if not _mat_eq(lm.matmul(Pv, Pv, dom), Pv, dom):
                return False, f"idempotent action e_{v} is not a projection: P_{v} @ P_{v} != P_{v}"
        for a in range(len(verts)):
            for b in range(a + 1, len(verts)):
                v, w = verts[a], verts[b]
                if not _mat_is_zero(lm.matmul(projs[v], projs[w], dom), dom):
                    return False, f"idempotent actions e_{v}, e_{w} are not orthogonal: P_{v} @ P_{w} != 0"
        # (ii) grading: each stored arrow action respects the source/target idempotents
        for a in self._arrow_labels():
            if a not in self.action:
                continue
            Aa = self.action[a]
            s = self.algebra.quiver.source(a)
            t = self.algebra.quiver.target(a)
            graded = lm.matmul(projs[t], lm.matmul(Aa, projs[s], dom), dom)
            if not _mat_eq(Aa, graded, dom):
                return False, (f"arrow action {a!r} ({s} -> {t}) violates the vertex grading: "
                               f"action[{a}] != P_{t} @ action[{a}] @ P_{s}")
        # (iii) relations: for each relation sum c_w * word, sum c_w * action[word] == 0
        for rel in (self.algebra.relations or []):
            M = lm.zeros(n, n, dom)
            for coeff, word in _relation_terms(rel, dom):
                M = _add(M, _scale(self._action_of_word(word), coeff, dom), dom)
            if not _mat_is_zero(M, dom):
                return False, f"relation not satisfied: {rel}"
        # (iv) multiplicativity on the composite basis-label products actually stored
        for label, mat in self.action.items():
            if label.startswith("e_") or "*" not in label:
                continue
            word = tuple(label.split("*"))
            if any(a not in self.action for a in word):
                continue
            if not _mat_eq(mat, self._action_of_word(word), dom):
                return False, (f"action[{label}] is not multiplicative: it differs from the "
                               f"anti-homomorphism composition action[{word[-1]}] @ ... "
                               f"@ action[{word[0]}]")
        return True, None

    def _action_of_word(self, word):
        """action of a path word (tuple of arrow names) by composing arrow actions in
        anti-homomorphism order: action[a1*...*ak] = action[ak] @ ... @ action[a1]."""
        dom = self.domain
        if word == ():
            return lm.identity(self.dim, dom)
        M = None
        for name in word:  # left to right; anti-homo => multiply on the LEFT
            Aa = self.action[name]
            M = Aa if M is None else lm.matmul(Aa, M, dom)
        return M

    @classmethod
    def from_arrow_action(cls, algebra, dimension_vector, arrow_action, name="M"):
        """Build a module from per-arrow action matrices plus the dimension vector.
        The idempotent actions are the block projections implied by dimension_vector
        (in the vertex-ordered basis), and every non-trivial basis-path label's action
        is composed from the arrow actions. Validated before return."""
        dom = algebra.domain
        verts = list(algebra.quiver.vertices)
        dims = [dimension_vector.get(v, 0) for v in verts]
        n = sum(dims)
        # basis ordered by vertex block: build idempotent projections
        action = {}
        offset = 0
        starts = {}
        for v, dv in zip(verts, dims):
            starts[v] = offset
            offset += dv
        for v, dv in zip(verts, dims):
            P = lm.zeros(n, n, dom)
            for i in range(starts[v], starts[v] + dv):
                P[i][i] = dom.one()
            action[f"e_{v}"] = P
        for aname, mat in arrow_action.items():
            action[aname] = mat
        M = cls(algebra, n, action, name=name)
        # fill every algebra basis-label action (paths + idempotents) by composition
        M._extend_to_basis_labels()
        ok, why = M.check_module()
        if not ok:
            raise QuiverlabError(f"from_arrow_action({name}): not a module: {why}",
                                 hint="check that the arrow matrices satisfy the relations")
        return M

    def _extend_to_basis_labels(self):
        """Ensure action[label] exists for every algebra basis label (idempotents and
        path words), computed by composing the stored arrow/idempotent actions."""
        for label in self.algebra.basis_labels:
            if label in self.action:
                continue
            if label.startswith("e_"):
                # already set for genuine vertices; any missing means a semisimple gap
                continue
            word = tuple(label.split("*"))
            self.action[label] = self._action_of_word(word)

    def identity_hom(self):
        """The identity homomorphism id_M as a ModuleHom (Plan 37 / C1)."""
        from quiverlab.modules import linalg_mod as lm
        from quiverlab.modules.morphism import ModuleHom
        return ModuleHom(self, self, lm.identity(self.dim, self.domain), check=False)

    def end_algebra(self):
        """End_A(M) as a structure-constant Algebra (Plan 37 / C1)."""
        from quiverlab.modules.endomorphism import end_algebra
        return end_algebra(self)

    # -- covers / envelopes as maps; radical & socle series (Plan 37 / C1) ----
    def projective_cover_hom(self):
        """The projective cover ``P(M) ->> M`` as a ModuleHom (epi, superfluous
        kernel ker ⊆ rad P(M)) (Plan 37)."""
        from quiverlab.modules.morphism import ModuleHom
        from quiverlab.modules.resolution import projective_cover
        Q0, d0, _ = projective_cover(self)
        return ModuleHom(Q0, self, d0, check=False)

    def injective_envelope_hom(self):
        """The injective envelope ``M >-> E(M)`` as a ModuleHom (mono, essential:
        soc E(M) = soc M) (Plan 37)."""
        from quiverlab.modules.injective import injective_resolution
        from quiverlab.modules.morphism import ModuleHom
        res = injective_resolution(self, 1)
        E0 = res.terms[0]
        iota = res.differential(0)              # iota: M -> E^0 (dual bases)
        return ModuleHom(self, E0, iota, check=False)

    def radical_series(self):
        """The descending radical (lower Loewy) series ``[M, rad M, rad^2 M, ..., 0]``,
        each term a Module, ending at the zero module (Plan 37)."""
        series = [self]
        cur = self
        while cur.dim > 0:
            r = cur.radical()
            if r.dim >= cur.dim:                # f.d.: the radical strictly shrinks
                break
            series.append(r)
            cur = r
        return series

    def socle_series(self):
        """The ascending socle (upper Loewy) series ``[0, soc M, soc^2 M, ..., M]``,
        each term a submodule of M (soc^{k+1}/soc^k = soc(M/soc^k)) (Plan 37)."""
        from quiverlab.modules import linalg_mod as lm
        from quiverlab.modules.radtopsoc import submodule
        from quiverlab.modules.yoneda import _quotient_with_maps
        dom = self.domain
        sub_cols = []
        socs = [submodule(self, [], name=f"soc^0 {self.name}")]
        k = 0
        while len(sub_cols) < self.dim:
            _Q, _proj, lift = _quotient_with_maps(
                self, sub_cols, dom, name=f"{self.name}/soc^{k}")
            soc_cols_Q = _socle_columns(_Q, dom)
            if not soc_cols_Q:                  # a nonzero module has a nonzero socle
                break
            lifted = [lm.matvec(lift, c, dom) for c in soc_cols_Q]
            sub_cols = [list(c) for c in sub_cols] + lifted
            k += 1
            socs.append(submodule(self, sub_cols, name=f"soc^{k} {self.name}"))
        return socs

    def loewy_layers(self):
        """The Loewy (radical) layers, top to bottom: ``[top(M), top(rad M),
        top(rad^2 M), ...]`` as str-keyed, vertex-sorted composition-factor
        multiplicity dicts (Plan 37). This is the public home of the logic
        ``trace.modules._radical_layers`` used to carry; the report renderers
        delegate here byte-for-byte."""
        layers = []
        cur = self
        while cur.dim > 0:
            layers.append(_normalize_dv(cur.top().dimension_vector()))
            r = cur.radical()
            if r.dim >= cur.dim:                # radical strictly shrinks for f.d.
                break
            cur = r
        return layers

    def composition_factors(self):
        """Total composition-factor multiplicities of M as a str-keyed dict
        (summed Loewy layers) (Plan 37)."""
        total = {}
        for layer in self.loewy_layers():
            for k, v in layer.items():
                total[k] = total.get(k, 0) + v
        return total

    def radical(self):
        from quiverlab.modules.radtopsoc import radical as _r
        return _r(self)

    def top(self):
        from quiverlab.modules.radtopsoc import top as _t
        return _t(self)

    def socle(self):
        from quiverlab.modules.radtopsoc import socle as _s
        return _s(self)

    def projective_resolution(self, length, max_term_dim=200000):
        from quiverlab.modules.resolution import minimal_resolution, ProjectiveResolution
        terms, dmats = minimal_resolution(self, length, max_term_dim=max_term_dim)
        return ProjectiveResolution(self, terms, dmats)

    # -- duality, transpose, AR translates (Plan 23; Plan 24 sides) -----------
    def dualize(self):
        """D M = Hom_k(M, k). Contravariant, exchanging the two sides over the SAME
        base algebra (Plan 24): D of a right A-module is a LEFT A-module and vice
        versa. D preserves dimension vectors and D.D = id."""
        from quiverlab.modules.duality import dualize
        return dualize(self)

    def transpose(self):
        """Tr M = coker(Hom(P_0,A) -> Hom(P_1,A)); Hom_A(-,A) lands in the OTHER-side
        modules, so Tr flips the side over the same base algebra (Plan 24)."""
        from quiverlab.modules.duality import transpose_module
        return transpose_module(self)

    def tau(self):
        """Auslander-Reiten translate tau M = D(Tr M). tau(projective) = 0."""
        from quiverlab.modules.duality import tau
        return tau(self)

    def tau_minus(self):
        """inverse AR translate tau^- M = Tr(D M). tau^-(injective) = 0."""
        from quiverlab.modules.duality import tau_minus
        return tau_minus(self)

    def nakayama(self):
        """The Nakayama functor value nu M = D Hom_A(M, A) (Plan 41 / C3). nu(P_v) = I_v;
        ker(nu P_1 -> nu P_0) = tau M ties it to the trusted AR translate."""
        from quiverlab.modules.ar import nakayama_functor
        return nakayama_functor(self)

    def nakayama_minus(self):
        """The inverse Nakayama functor nu^- M = Hom_A(DA, M) (Plan 41 / C3).
        nu^-(I_v) = P_v."""
        from quiverlab.modules.ar import nakayama_functor_minus
        return nakayama_functor_minus(self)

    def is_isomorphic(self, other):
        """True iff self and other are isomorphic right modules (exact certificate)."""
        from quiverlab.modules.hom import is_isomorphic
        return is_isomorphic(self, other)

    # -- Krull-Schmidt decomposition (Plan 30) --------------------------------
    def decompose(self, budget=None):
        """Krull-Schmidt decomposition ``[(M_i, m_i), ...]`` into indecomposable
        summands with multiplicities, ``(+) M_i^{m_i} ~ self`` (each summand certified
        indecomposable; grouped up to iso). Raises loudly when a summand cannot be
        certified within budget -- never a silent wrong answer.

        tau-additivity: ``tau(self) = (+) tau(M_i)^{m_i}`` (the AR translate is additive),
        so a translate of a decomposable module is certified summand-wise."""
        from quiverlab.modules.decompose import decompose
        return decompose(self) if budget is None else decompose(self, budget=budget)

    def is_indecomposable(self, budget=None):
        """True iff self is indecomposable (End local, certified); False iff a Fitting
        split exists. Raises loudly if undecidable within budget (see decompose)."""
        from quiverlab.modules.decompose import is_indecomposable
        return (is_indecomposable(self) if budget is None
                else is_indecomposable(self, budget=budget))

    def injective_resolution(self, length, max_term_dim=200000):
        """Minimal injective coresolution 0 -> M -> E^0 -> E^1 -> ... (Plan 23)."""
        from quiverlab.modules.injective import injective_resolution
        return injective_resolution(self, length, max_term_dim=max_term_dim)

    def injective_dimension(self, bound=32, max_term_dim=200000):
        """inj.dim_A(M) = pd_{A^op}(DM): int, or None if unresolved within bound."""
        from quiverlab.modules.injective import injective_dimension
        return injective_dimension(self, bound=bound, max_term_dim=max_term_dim)

    def __repr__(self):
        dv = self.dimension_vector()
        # side="right" reproduces the pre-Plan-24 string byte-for-byte (base_algebra
        # is self.algebra for right modules); left modules show "left <base algebra>".
        return (f"{self.name}: {self.side} {self.base_algebra} module, "
                f"dim {self.dim}, dimvec {dv}")


def _coerce_matrix(mat, dom):
    """Normalize plain int/rational literals into the domain; entries that are already
    field elements (e.g. sympy MPQ over CC, ints over GF(p)) pass through the coercion
    as a no-op or, for native non-int/Fraction elements, unchanged."""
    return [[dom.coerce(x) if isinstance(x, (int, Fraction)) else x for x in row]
            for row in mat]


def _normalize_dv(dimvec):
    """A dimension vector as a str-keyed, vertex-sorted dict of ints. Byte-identical
    to trace.modules._dv, which now delegates the Loewy layers here (Plan 37)."""
    return {str(v): int(n)
            for v, n in sorted(dimvec.items(), key=lambda kv: str(kv[0]))}


def _socle_columns(Q, dom):
    """The socle basis columns of a module Q in Q's own coordinates: the intersection
    over the arrows of ker(action[arrow]) (mirrors radtopsoc.socle, returning the
    spanning columns rather than the submodule so socle_series can lift them)."""
    from quiverlab.modules.radtopsoc import _intersect
    arrows = list(Q.algebra.quiver.arrows)
    if not arrows:                              # semisimple: soc Q = Q
        ident = lm.identity(Q.dim, dom)
        return [lm.col(ident, j) for j in range(Q.dim)]
    inter = None
    for a in arrows:
        ker = lm.kernel_columns(Q.action[a], dom)
        inter = ker if inter is None else _intersect(inter, ker, dom)
    return inter or []


def _add(A, B, dom):
    return [[dom.add(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def _mat_is_zero(M, dom):
    return all(dom.is_zero(x) for row in M for x in row)


def _mat_eq(A, B, dom):
    """Exact matrix equality through the domain (dom.is_zero(dom.sub(...)) elementwise),
    never a raw `!=` on the underlying field elements."""
    if len(A) != len(B):
        return False
    return all(len(ra) == len(rb) and all(dom.is_zero(dom.sub(x, y)) for x, y in zip(ra, rb))
               for ra, rb in zip(A, B))


def _scale(A, c, dom):
    return [[dom.mul(c, x) for x in row] for row in A]


def _relation_terms(rel, dom):
    """Yield (domain coeff, word) pairs of a Plan-01 Relation, coefficients coerced."""
    for coeff, word in rel.terms:
        yield dom.coerce(coeff), tuple(word)
