"""Skew group algebras ``A rtimes G`` (the smash product ``A # kG``) as first-class
inputs (Plan 74 / R8).

A finite group ``G`` acting on ``A = kQ/I`` by algebra automorphisms (given
EXPLICITLY, never inferred) yields the smash product

    ``A rtimes G = A (x) kG``,  basis ``{a.g : a in A-basis, g in G}``,
    product ``(a.g)(b.h) = (a . g(b)) . (gh)``,   dim = |G| * dim A.

The action is specified on the QUIVER (primary, no-code): a
:class:`QuiverAutomorphism` triple ``(vertex_perm, arrow_perm, arrow_scalars)``
that induces an algebra automorphism ``phi(e_v)=e_{pi v}``, ``phi(a)=s(a).rho(a)``,
extended multiplicatively; a per-instance :meth:`QuiverAutomorphism.check`
certifies it preserves ``I`` (is an algebra automorphism). A :class:`GroupAction`
closes the generators to a finite group (Cayley table, conjugacy classes,
centralizers).

The constructor is CHARACTERISTIC-AGNOSTIC (the smash structure constants never
divide by ``|G|``); only the Stefan conjugacy-class HH decomposition
(``quiverlab.hochschild.skew_group``) needs ``char k does not divide |G|``.

References: Cibils-Marcos (Proc. AMS 134, 2006, the smash construction + free
action + Galois covering); Shepler-Witherspoon (J. Algebra 351, 2012) and
Stefan (JPAA 103, 1995) for the HH decomposition consumed by the hochschild twin.
"""
from fractions import Fraction

from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields.domain import reject_inexact
from quiverlab.fields.linalg import rank, solve

_CLOSURE_BUDGET = 4096          # max group order the BFS closure will build


# --------------------------------------------------------------------------- #
# small exact matrix helpers (Domain entries)
# --------------------------------------------------------------------------- #
def _eye(dom, m):
    return [[dom.one() if i == j else dom.zero() for j in range(m)] for i in range(m)]


def _mat_mul(dom, X, Y):
    n, k, w = len(X), len(Y), len(Y[0])
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


def _mat_apply(dom, M, v):
    """M @ v where columns of M are images of basis vectors."""
    n = len(M)
    out = [dom.zero()] * n
    for i in range(n):
        Mi = M[i]
        acc = dom.zero()
        for j, vj in enumerate(v):
            if not dom.is_zero(vj) and not dom.is_zero(Mi[j]):
                acc = dom.add(acc, dom.mul(Mi[j], vj))
        out[i] = acc
    return out


def _mat_key(dom, M):
    """A hashable canonical key for a Domain matrix (for group-closure dedup).
    Domain elements are stored reduced, so ``repr`` is a canonical fingerprint."""
    return tuple(tuple(repr(x) for x in row) for row in M)


def _mat_pow(dom, M, n):
    R = _eye(dom, len(M))
    for _ in range(n):
        R = _mat_mul(dom, R, M)
    return R


# --------------------------------------------------------------------------- #
# QuiverAutomorphism
# --------------------------------------------------------------------------- #
class QuiverAutomorphism:
    """A quiver automorphism ``(pi, rho, s)`` inducing an algebra automorphism of
    ``A = kQ/I``: ``phi(e_v) = e_{pi(v)}``, ``phi(alpha) = s(alpha) . rho(alpha)``,
    extended multiplicatively along paths.

    ``vertex_perm`` maps each vertex to a vertex; ``arrow_perm`` maps each arrow
    name to an arrow name; ``arrow_scalars`` (optional) gives a nonzero exact
    scalar per arrow (default ``1``). The scalars are essential: ``x |-> -x`` on
    ``k[x]/(x^2)`` is the identity permutation with ``s(x) = -1``.
    """

    def __init__(self, vertex_perm, arrow_perm, arrow_scalars=None):
        self.vertex_perm = dict(vertex_perm)
        self.arrow_perm = dict(arrow_perm)
        raw = {} if arrow_scalars is None else dict(arrow_scalars)
        for a, sc in raw.items():
            reject_inexact(sc)              # floats refused loudly
        self.arrow_scalars = {a: (raw[a] if a in raw else 1) for a in self.arrow_perm}

    # -- composition (self after other): (self o other)(x) = self(other(x)) ----
    def compose(self, other):
        vp = {v: self.vertex_perm[other.vertex_perm[v]] for v in other.vertex_perm}
        ap = {a: self.arrow_perm[other.arrow_perm[a]] for a in other.arrow_perm}
        sc = {}
        for a in other.arrow_perm:
            s_other = _as_fraction(other.arrow_scalars[a])
            s_self = _as_fraction(self.arrow_scalars[other.arrow_perm[a]])
            sc[a] = s_other * s_self
        return QuiverAutomorphism(vp, ap, sc)

    def is_identity_perm(self):
        return (all(v == w for v, w in self.vertex_perm.items())
                and all(a == b for a, b in self.arrow_perm.items()))

    # -- the induced m x m automorphism matrix in A's own basis -----------------
    def matrix(self, A):
        """The ``m x m`` Domain matrix whose column ``i`` is ``phi(b_i)`` (columns =
        images of ``A``'s basis). The form :meth:`Bimodule.twisted` consumes."""
        self._validate_shape(A)
        dom = A.domain
        m = A.dim
        labels = A.basis_labels
        label_index = {lab: i for i, lab in enumerate(labels)}
        arrow_names = set(A.quiver.arrows)
        images = [None] * m
        # generators: vertices e_v and arrows
        for i, lab in enumerate(labels):
            if isinstance(lab, str) and lab.startswith("e_"):
                v = _vertex_of_label(lab, A)
                tgt = f"e_{self.vertex_perm[v]}"
                if tgt not in label_index:
                    raise QuiverlabError(
                        f"vertex automorphism sends {lab} to missing vertex {tgt}",
                        hint="vertex_perm must be a permutation of the quiver vertices")
                images[i] = A._basis_vec(label_index[tgt])
            elif lab in arrow_names:
                b = self.arrow_perm[lab]
                if b not in label_index:
                    raise QuiverlabError(
                        f"arrow automorphism sends {lab} to missing arrow {b}",
                        hint="arrow_perm must be a permutation of the quiver arrows")
                sc = dom.coerce(self.arrow_scalars[lab])
                if dom.is_zero(sc):
                    raise QuiverlabError(
                        f"arrow scalar for {lab!r} is zero (must be in k^x)",
                        hint="arrow scalars are nonzero units")
                vec = A._basis_vec(label_index[b])
                images[i] = [dom.mul(sc, x) for x in vec]
        # fold to path basis elements via a factorization b_i = arrow . b_k
        arrow_idx = [label_index[a] for a in arrow_names if a in label_index]
        remaining = [i for i in range(m) if images[i] is None]
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
                        if images[k] is None:
                            continue
                        if A.multiply(ea, A._basis_vec(k)) == ei:
                            images[i] = A.multiply(images[ai], images[k])
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
                "quiver automorphism could not be folded to every basis element",
                hint="the basis is not generated by the vertices and arrows given")
        # columns = images
        return [[images[j][t] for j in range(m)] for t in range(m)]

    def check(self, A):
        """Certify (loud) that the induced map is a bijective, unital algebra
        automorphism of ``A`` (equivalently preserves ``I``). Raises on failure."""
        dom = A.domain
        m = A.dim
        M = self.matrix(A)
        if rank(M, dom) != m:
            raise QuiverlabError(
                "the quiver automorphism is not bijective on A",
                hint="an automorphism must be an invertible linear map")
        # unital: phi(1_A) = 1_A
        phi_unit = _mat_apply(dom, M, list(A.unit))
        if phi_unit != list(A.unit):
            raise QuiverlabError(
                "the quiver automorphism does not fix the unit of A",
                hint="phi(1_A) must equal 1_A")
        # multiplicative: phi(b_i b_j) = phi(b_i) phi(b_j)
        for i in range(m):
            phi_i = _mat_apply(dom, M, A._basis_vec(i))
            for j in range(m):
                lhs = _mat_apply(dom, M, A.T[i][j])
                rhs = A.multiply(phi_i, _mat_apply(dom, M, A._basis_vec(j)))
                if lhs != rhs:
                    raise QuiverlabError(
                        f"the quiver automorphism does not preserve the relations of A "
                        f"(fails multiplicativity at basis ({i}, {j}))",
                        hint="phi does not descend to A = kQ/I; it is not an automorphism")
        return True

    def _validate_shape(self, A):
        if A.quiver is None or not A.basis_labels:
            raise QuiverlabError(
                "a quiver automorphism needs a quiver-presented algebra",
                hint="build A via Quiver(...).algebra(...)")
        verts = set(A.quiver.vertices)
        if set(self.vertex_perm) != verts or set(self.vertex_perm.values()) != verts:
            raise QuiverlabError(
                "vertex_perm is not a permutation of the quiver vertices",
                hint=f"vertices are {sorted(verts, key=str)}")
        arrows = set(A.quiver.arrows)
        if set(self.arrow_perm) != arrows or set(self.arrow_perm.values()) != arrows:
            raise QuiverlabError(
                "arrow_perm is not a permutation of the quiver arrows",
                hint=f"arrows are {sorted(arrows)}")
        # compatibility: source(rho a) = pi(source a), target(rho a) = pi(target a)
        for a in arrows:
            s, t = A.quiver.arrows[a]
            b = self.arrow_perm[a]
            bs, bt = A.quiver.arrows[b]
            if bs != self.vertex_perm[s] or bt != self.vertex_perm[t]:
                raise QuiverlabError(
                    f"arrow_perm is incompatible with vertex_perm at arrow {a!r}",
                    hint="need source(rho a)=pi(source a), target(rho a)=pi(target a)")

    def canonical_key(self):
        """A stable hashable key (for canonical-key generator-order normalization)."""
        vp = tuple(sorted((str(v), str(w)) for v, w in self.vertex_perm.items()))
        ap = tuple(sorted((str(a), str(b)) for a, b in self.arrow_perm.items()))
        sc = tuple(sorted((str(a), str(_as_fraction(self.arrow_scalars[a])))
                          for a in self.arrow_scalars))
        return (vp, ap, sc)


def _as_fraction(x):
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, str):
        return Fraction(x)
    # exact Domain element with an int value (GF(p)) -> keep symbolic via repr
    try:
        return Fraction(int(x))
    except (TypeError, ValueError):
        raise QuiverlabError(f"arrow scalar {x!r} is not an exact rational/int",
                             hint="use int, a 'p/q' string, or -1")


def _vertex_of_label(lab, A):
    tail = lab[2:]
    for v in A.quiver.vertices:
        if str(v) == tail:
            return v
    raise QuiverlabError(f"vertex label {lab!r} does not name a quiver vertex",
                         hint="internal invariant")


# --------------------------------------------------------------------------- #
# GroupAction
# --------------------------------------------------------------------------- #
class _MatrixElement:
    """A power-user group element given by an explicit ``m x m`` automorphism
    matrix (no quiver data; not exposed in the no-code GUI)."""

    def __init__(self, matrix):
        self._matrix = matrix
        self.vertex_perm = None

    def matrix(self, A):
        return [list(r) for r in self._matrix]


class GroupAction:
    """A finite group acting on ``A`` by algebra automorphisms, closed from a
    finite list of generators. ``elements`` is the closed set (Cayley-indexed,
    ``elements[0]`` = identity); ``conjugacy_classes`` and :meth:`centralizer`
    are derived combinatorially."""

    def __init__(self, generators, elements, cayley, conjugacy_classes,
                 centralizers, order, *, matrix_mode=False):
        self.generators = list(generators)
        self.elements = elements
        self._cayley = cayley
        self.conjugacy_classes = conjugacy_classes
        self._centralizers = centralizers
        self.order = order
        self.matrix_mode = matrix_mode

    def centralizer(self, g_index):
        return list(self._centralizers[g_index])

    def matrix_of(self, i, A):
        return self.elements[i].matrix(A)

    def product(self, i, j):
        return self._cayley[i][j]

    # -- constructors ---------------------------------------------------------
    @classmethod
    def cyclic(cls, generator, order):
        """The cyclic group ``<generator>`` of the given order (abelian: every
        class is a singleton, every centralizer is the whole group). Built
        abstractly (no field needed); :meth:`check` certifies the order over the
        field when the algebra is supplied."""
        if not (isinstance(order, int) and order >= 1):
            raise QuiverlabError(f"cyclic order {order!r} must be a positive int",
                                 hint="e.g. GroupAction.cyclic(sigma, 2)")
        ident = QuiverAutomorphism(
            {v: v for v in generator.vertex_perm},
            {a: a for a in generator.arrow_perm},
            {a: 1 for a in generator.arrow_perm})
        elements = [ident]
        cur = ident
        for _ in range(1, order):
            cur = cur.compose(generator)
            elements.append(cur)
        cayley = [[(i + j) % order for j in range(order)] for i in range(order)]
        conj = [[i] for i in range(order)]
        cent = [list(range(order)) for _ in range(order)]
        return cls([generator], elements, cayley, conj, cent, order)

    @classmethod
    def from_generators(cls, generators, A, *, budget=_CLOSURE_BUDGET):
        """Close the quiver-automorphism ``generators`` to a finite group by a
        matrix BFS over ``A``'s field (field-correct: the order of an arrow scalar
        depends on the characteristic). Loud on non-closure within ``budget``."""
        dom = A.domain
        ident = QuiverAutomorphism(
            {v: v for v in A.quiver.vertices},
            {a: a for a in A.quiver.arrows},
            {a: 1 for a in A.quiver.arrows})
        return cls._close(generators, ident, A, dom, budget, matrix_mode=False)

    @classmethod
    def from_matrices(cls, matrices, A, *, budget=_CLOSURE_BUDGET):
        """Power-user twin: close explicit ``m x m`` automorphism matrices over
        ``A``'s field. Not exposed in the no-code GUI."""
        dom = A.domain
        gens = [_MatrixElement([[dom.coerce(x) for x in row] for row in M])
                for M in matrices]
        ident = _MatrixElement(_eye(dom, A.dim))
        return cls._close(gens, ident, A, dom, budget, matrix_mode=True)

    @classmethod
    def _close(cls, generators, ident, A, dom, budget, *, matrix_mode):
        gen_mats = [g.matrix(A) for g in generators]
        I = ident.matrix(A)
        elements = [ident]
        mats = [I]
        seen = {_mat_key(dom, I): 0}
        frontier = [(ident, I)]
        while frontier:
            nf = []
            for (elt, M) in frontier:
                for gi, G in enumerate(gen_mats):
                    P = _mat_mul(dom, M, G)
                    key = _mat_key(dom, P)
                    if key not in seen:
                        if len(elements) >= budget:
                            raise QuiverlabError(
                                f"group closure exceeded the budget of {budget} elements "
                                "(is the action infinite?)",
                                hint="check the generator orders / arrow scalars")
                        seen[key] = len(elements)
                        new_elt = (elt.compose(generators[gi]) if not matrix_mode
                                   else _MatrixElement(P))
                        elements.append(new_elt)
                        mats.append(P)
                        nf.append((new_elt, P))
            frontier = nf
        order = len(elements)
        idx = seen
        cayley = [[idx[_mat_key(dom, _mat_mul(dom, mats[i], mats[j]))]
                   for j in range(order)] for i in range(order)]
        inv = [next(j for j in range(order) if cayley[i][j] == 0) for i in range(order)]
        # conjugacy classes + centralizers from the Cayley table
        conj, cent = _classes_and_centralizers(order, cayley, inv)
        return cls(list(generators), elements, cayley, conj, cent, order,
                   matrix_mode=matrix_mode)

    # -- certification --------------------------------------------------------
    def check(self, A):
        """Certify (loud) that every generator is an algebra automorphism of ``A``
        and that the recorded structure is a group of the claimed order over
        ``A``'s field (the Cayley products match the matrices)."""
        dom = A.domain
        for g in self.generators:
            if isinstance(g, QuiverAutomorphism):
                g.check(A)
        mats = [e.matrix(A) for e in self.elements]
        keys = {_mat_key(dom, M): i for i, M in enumerate(mats)}
        if len(keys) != self.order:
            raise QuiverlabError(
                "the recorded group elements are not distinct over this field "
                f"({len(keys)} distinct of {self.order}) -- the claimed order is wrong "
                "for this characteristic",
                hint="over GF(p) an arrow scalar's order can differ; pass the true order")
        for i in range(self.order):
            for j in range(self.order):
                prod = _mat_key(dom, _mat_mul(dom, mats[i], mats[j]))
                if prod not in keys or keys[prod] != self._cayley[i][j]:
                    raise QuiverlabError(
                        "the group is not closed under composition over this field "
                        f"(product of elements {i}, {j} escaped)",
                        hint="the generators do not generate the claimed finite group")
        return True


def _classes_and_centralizers(order, cayley, inv):
    seen = set()
    classes = []
    for g in range(order):
        if g in seen:
            continue
        cls_g = set()
        for k in range(order):
            # k g k^{-1}
            cls_g.add(cayley[cayley[k][g]][inv[k]])
        seen |= cls_g
        classes.append(sorted(cls_g))
    cent = []
    for g in range(order):
        cent.append([k for k in range(order) if cayley[k][g] == cayley[g][k]])
    return classes, cent


# --------------------------------------------------------------------------- #
# the smash-product constructor
# --------------------------------------------------------------------------- #
def skew_group_algebra(A, action, *, field=None):
    """The skew group algebra ``A rtimes G`` in structure constants, dimension
    ``|G| * dim A``, basis labels ``"<a>.<g>"``, product
    ``(a.g)(b.h) = (a . g(b)) . (gh)``. Characteristic-agnostic. Certified per
    instance (associativity via ``_validate`` + the dim-law assert)."""
    if A.quiver is None:
        raise QuiverlabError(
            "skew_group_algebra needs a quiver-presented base A (the action is quiver data)",
            hint="build A via Quiver(...).algebra(...)")
    dom = A.domain
    m = A.dim
    action = _ensure_bound(A, action)
    order = action.order
    mats = [action.matrix_of(g, A) for g in range(order)]
    d = m * order
    zero = dom.zero()

    def n(a, g):
        return a * order + g

    T = [[[zero] * d for _ in range(d)] for _ in range(d)]
    for a in range(m):
        ea = A._basis_vec(a)
        for g in range(order):
            na = n(a, g)
            Trow = T[na]
            for b in range(m):
                gb = _mat_apply(dom, mats[g], A._basis_vec(b))  # g(b)
                prod = A.multiply(ea, gb)                       # a . g(b)
                for h in range(order):
                    gh = action.product(g, h)
                    nb = n(b, h)
                    cell = Trow[nb]
                    for c in range(m):
                        if not dom.is_zero(prod[c]):
                            cell[c * order + gh] = prod[c]
    unit = [zero] * d
    for c in range(m):
        if not dom.is_zero(A.unit[c]):
            unit[n(c, 0)] = A.unit[c]
    base_labels = A.basis_labels or [f"a{i}" for i in range(m)]
    labels = [f"{base_labels[a]}.g{g}" for a in range(m) for g in range(order)]
    S = Algebra(dom, T, unit, basis_labels=labels)
    S._validate()
    if S.dim != order * m:
        raise QuiverlabError(
            f"dim law failed: dim(A rtimes G) = {S.dim} != |G|*dim A = {order * m}",
            hint="internal invariant -- please report this action")
    S._family_citations = ("cibils_marcos_smash",)
    # stamp the base + action so the skew_group_hh decomposition kind can recover
    # them (a smash algebra flows through every other compute kind as a plain
    # structure-constant Algebra; only the Stefan decomposition needs the origin).
    S._skew_base = A
    S._skew_action = action
    return S


def build_skew_group_from_params(params, field):
    """Build ``A rtimes G`` from a flattened GUI/HPC ``SkewGroupAlgebra`` param block:
    ``vertices`` (list of ints), ``arrows`` (name -> [source, target]), ``relations``
    (list of strings), ``generators`` (a list of quiver-automorphism dicts
    ``{"vertex_perm", "arrow_perm", "arrow_scalars"}``). The generator list is
    CANONICALIZED (each generator to a stable form, the list sorted + de-duplicated)
    so two orderings of the same action build the same algebra. Loud
    ``QuiverlabError`` on any ill-typed param; the group is closed field-correctly."""
    from quiverlab.combinat.quiver import Quiver
    verts = params.get("vertices")
    if not (isinstance(verts, (list, tuple)) and verts
            and all(isinstance(v, int) and not isinstance(v, bool) for v in verts)):
        raise QuiverlabError(
            "SkewGroupAlgebra.vertices must be a non-empty list of vertex integers",
            hint="e.g. [1] for the one-loop base of k[x]/(x^2)")
    arrows_p = params.get("arrows") or {}
    if not isinstance(arrows_p, dict):
        raise QuiverlabError("SkewGroupAlgebra.arrows must map names to [source, target]",
                             hint='e.g. {"x": [1, 1]}')
    arrows = {}
    for name, st in arrows_p.items():
        if not (isinstance(st, (list, tuple)) and len(st) == 2):
            raise QuiverlabError(
                f"SkewGroupAlgebra.arrows[{name!r}] must be a [source, target] pair",
                hint='e.g. {"x": [1, 1]}')
        arrows[str(name)] = (st[0], st[1])
    rels = params.get("relations") or []
    if not (isinstance(rels, (list, tuple)) and all(isinstance(r, str) for r in rels)):
        raise QuiverlabError("SkewGroupAlgebra.relations must be a list of strings",
                             hint='e.g. ["x*x"]')
    A = Quiver(vertices=list(verts), arrows=arrows).algebra(
        relations=list(rels), field=field)
    gen_dicts = params.get("generators")
    if not (isinstance(gen_dicts, (list, tuple)) and gen_dicts):
        raise QuiverlabError(
            "SkewGroupAlgebra.generators must be a non-empty list of quiver-automorphism "
            "dicts {vertex_perm, arrow_perm, arrow_scalars}",
            hint='e.g. [{"vertex_perm": {"1": 1}, "arrow_perm": {"x": "x"}, '
                 '"arrow_scalars": {"x": -1}}]')
    gens = [_parse_generator_dict(g, A) for g in gen_dicts]
    # certify each generator is a genuine algebra automorphism BEFORE closing the group
    # (a non-automorphism gives a singular matrix that breaks the closure) -- loud.
    for g in gens:
        g.check(A)
    # canonical generator-order normalization: sort + de-duplicate by canonical key
    seen, canon = set(), []
    for g in sorted(gens, key=lambda q: q.canonical_key()):
        k = g.canonical_key()
        if k not in seen:
            seen.add(k)
            canon.append(g)
    action = GroupAction.from_generators(canon, A)
    return skew_group_algebra(A, action)


def _parse_generator_dict(g, A):
    """Coerce one JSON generator dict into a QuiverAutomorphism, mapping vertex
    labels (possibly strings from JSON) back to the quiver's vertex type."""
    if not isinstance(g, dict):
        raise QuiverlabError("each SkewGroupAlgebra generator must be a dict",
                             hint='{"vertex_perm": ..., "arrow_perm": ..., "arrow_scalars": ...}')
    vlabel = {str(v): v for v in A.quiver.vertices}

    def cv(x):
        if x in vlabel:
            return vlabel[x]
        if str(x) in vlabel:
            return vlabel[str(x)]
        raise QuiverlabError(f"generator names an unknown vertex {x!r}",
                             hint=f"vertices are {sorted(A.quiver.vertices, key=str)}")

    vp_raw = g.get("vertex_perm") or {}
    vertex_perm = {cv(k): cv(v) for k, v in vp_raw.items()}
    # a vertex omitted from vertex_perm is fixed
    for v in A.quiver.vertices:
        vertex_perm.setdefault(v, v)
    ap_raw = g.get("arrow_perm") or {}
    arrow_perm = {str(k): str(v) for k, v in ap_raw.items()}
    for a in A.quiver.arrows:
        arrow_perm.setdefault(a, a)
    scalars = g.get("arrow_scalars") or {}
    arrow_scalars = {str(k): v for k, v in scalars.items()}
    return QuiverAutomorphism(vertex_perm, arrow_perm, arrow_scalars)


def _ensure_bound(A, action):
    """Return an action whose group structure is realized over ``A``'s field. A
    cyclic/abstract action passes through (its matrices realize the order over the
    field, certified by :meth:`GroupAction.check`); anything else is already
    bound."""
    if not isinstance(action, GroupAction):
        raise QuiverlabError(
            "action must be a GroupAction",
            hint="build it via GroupAction.cyclic(...) or GroupAction.from_generators(...)")
    action.check(A)
    return action


def is_free_action(A, action):
    """True iff the ``G``-action is FREE on the vertices ``Q_0``: no non-identity
    element fixes a vertex (a finite check enabling the Galois-covering reading,
    Cibils-Marcos math/0312214). Requires the quiver route (vertex data)."""
    if not isinstance(action, GroupAction) or action.matrix_mode:
        raise QuiverlabError(
            "is_free_action needs the quiver route (vertex permutations)",
            hint="build the action from QuiverAutomorphism generators, not matrices")
    for i, g in enumerate(action.elements):
        if i == 0:
            continue
        vp = g.vertex_perm
        if vp is None:
            raise QuiverlabError(
                "is_free_action needs vertex-permutation data on every element",
                hint="use the quiver route")
        if any(v == w for v, w in vp.items()):
            return False
    return True
