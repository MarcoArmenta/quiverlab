"""The right A-module object and its radical/top/socle (spec §3.6).

RIGHT modules, anti-homomorphism convention (see the package docstring): an element m
is a COLUMN vector in a fixed k-basis of M; the action of an algebra basis element b is
the matrix action[b] with m*b = action[b] @ m, and action[x*y] = action[y] @ action[x].
The vertex subspace M*e_v is the image of action['e_v']; dimension_vector[v] = its rank.
"""
from fractions import Fraction

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm


class Module:
    def __init__(self, algebra, dim, action, name="M"):
        self.algebra = algebra
        self.domain = algebra.domain
        self.dim = dim
        # Action entries may be given as plain int / rational literals for
        # convenience; normalize them into the algebra's exact field (mirroring
        # Algebra.from_structure_constants). Entries already in the field (e.g.
        # sympy MPQ over CC) are used as-is.
        self.action = {label: _coerce_matrix(mat, self.domain)
                       for label, mat in action.items()}
        self.name = name

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

    def radical(self):
        from quiverlab.modules.radtopsoc import radical as _r
        return _r(self)

    def top(self):
        from quiverlab.modules.radtopsoc import top as _t
        return _t(self)

    def socle(self):
        from quiverlab.modules.radtopsoc import socle as _s
        return _s(self)

    def projective_resolution(self, length):
        from quiverlab.modules.resolution import minimal_resolution, ProjectiveResolution
        terms, dmats = minimal_resolution(self, length)
        return ProjectiveResolution(self, terms, dmats)

    def __repr__(self):
        dv = self.dimension_vector()
        return f"{self.name}: right {self.algebra} module, dim {self.dim}, dimvec {dv}"


def _coerce_matrix(mat, dom):
    """Normalize plain int/rational literals into the domain; entries that are already
    field elements (e.g. sympy MPQ over CC, ints over GF(p)) pass through the coercion
    as a no-op or, for native non-int/Fraction elements, unchanged."""
    return [[dom.coerce(x) if isinstance(x, (int, Fraction)) else x for x in row]
            for row in mat]


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
