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
        """Verify the representation is a genuine right A-module: (i) the idempotents
        sum to the identity and are orthogonal projections; (ii) every relation is
        satisfied; (iii) action is multiplicative on the basis-label products it is
        given. Returns (True, None) or (False, witness)."""
        dom = self.domain
        n = self.dim
        # (i) sum of idempotent actions == I
        acc = lm.zeros(n, n, dom)
        for v in self.algebra.quiver.vertices:
            P = self.vertex_projection(v)
            acc = _add(acc, P, dom)
        if acc != lm.identity(n, dom):
            return False, "sum of e_v actions != identity"
        # (ii) relations: for each relation sum c_w * word, sum c_w * action[word] == 0
        for rel in (self.algebra.relations or []):
            M = lm.zeros(n, n, dom)
            for coeff, word in _relation_terms(rel, dom):
                M = _add(M, _scale(self._action_of_word(word), coeff, dom), dom)
            if any(not dom.is_zero(x) for row in M for x in row):
                return False, f"relation not satisfied: {rel}"
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


def _scale(A, c, dom):
    return [[dom.mul(c, x) for x in row] for row in A]


def _relation_terms(rel, dom):
    """Yield (domain coeff, word) pairs of a Plan-01 Relation, coefficients coerced."""
    for coeff, word in rel.terms:
        yield dom.coerce(coeff), tuple(word)
