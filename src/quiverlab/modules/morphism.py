"""First-class module homomorphisms (Plan 37 / C1).

A :class:`ModuleHom` is ``(src, tgt, matrix)``: ``matrix`` is dense
``tgt.dim x src.dim`` over the shared :class:`~quiverlab.fields.domain.Domain`,
satisfying ``tgt.action[b] @ matrix == matrix @ src.action[b]`` for every generator
label ``b`` (arrows and vertex idempotents) -- validated at construction
(``check=True``, via the shared predicate :func:`quiverlab.modules.yoneda._is_module_map`).

Composition is written LEFT-TO-RIGHT like paths: ``f.then(g)`` is
``src --f--> tgt --g--> g.tgt`` with matrix ``g.matrix @ f.matrix``. ``*`` is never
overloaded (the house convention). This module also provides kernel / image / cokernel
as Modules-with-maps (Task 2) and ``k``-ary direct sums / ``is_direct_summand``
(Task 4); short exact sequences and pushout/pullback live in
:mod:`quiverlab.modules.ses`.

Float-free: all arithmetic goes through the Domain (the AST gate scans this file).
"""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import _assert_comparable, hom_space
from quiverlab.modules.module import Module, _coerce_matrix
from quiverlab.modules.yoneda import _is_module_map


# --------------------------------------------------------------------------- #
# Single-chokepoint Domain shims (the plan's `add`/`negate`/`zero` adapted to the
# real Domain API: dom.zero()/one() are METHODS, negation is dom.neg). Task 3's
# ses.py imports these so pushout/pullback speak one arithmetic vocabulary.
# --------------------------------------------------------------------------- #
def _dadd(dom, x, y):
    return dom.add(x, y)


def _dneg(dom, x):
    return dom.neg(x)


def _cols_to_mat(cols, nrows, dom):
    """Assemble the ``nrows x len(cols)`` matrix whose columns are ``cols``. Unlike
    :func:`linalg_mod.cols_to_matrix` this keeps the row count when ``cols`` is empty
    (a genuine ``nrows x 0`` matrix, not the shapeless ``[]``) -- the yoneda precedent
    for the empty-column case."""
    return lm.cols_to_matrix(cols) if cols else lm.zeros(nrows, 0, dom)


class ModuleHom:
    """A validated A-module homomorphism ``src -> tgt`` (a dense ``tgt.dim x src.dim``
    matrix over the shared Domain)."""

    def __init__(self, src, tgt, matrix, check=True):
        _assert_comparable(src, tgt, "Hom")
        self.src, self.tgt = src, tgt
        self.domain = src.domain
        self.matrix = _coerce_matrix(matrix, self.domain)
        if len(self.matrix) != tgt.dim or (tgt.dim and any(
                len(r) != src.dim for r in self.matrix)):
            raise QuiverlabError(
                f"not a module map: expected a {tgt.dim}x{src.dim} matrix "
                f"(got {len(self.matrix)}x{len(self.matrix[0]) if self.matrix else 0})")
        if check and not _is_module_map(self.matrix, src, tgt, self.domain):
            raise QuiverlabError(
                "not a module map: the matrix does not intertwine the actions "
                "(tgt.action[b] @ f != f @ src.action[b] for some generator b)")

    def then(self, g: "ModuleHom") -> "ModuleHom":
        """``self`` then ``g`` (left-to-right): ``src --self--> tgt --g--> g.tgt``,
        matrix ``g.matrix @ self.matrix``."""
        if g.src is not self.tgt:
            raise QuiverlabError(
                "cannot compose: middle modules differ "
                f"({self.tgt.name} vs {g.src.name})")
        dom = self.domain
        # matmul returns [] (shapeless) when the inner dimension is 0; composition
        # through the zero module is the zero map of the outer shape -- build it.
        if self.tgt.dim == 0:
            mat = lm.zeros(g.tgt.dim, self.src.dim, dom)
        else:
            mat = lm.matmul(g.matrix, self.matrix, dom)
        return ModuleHom(self.src, g.tgt, mat, check=False)

    def rank(self) -> int:
        m = self.matrix
        return lm.mat_rank(m, self.domain) if (m and m[0]) else 0

    def is_zero(self) -> bool:
        dom = self.domain
        return all(dom.is_zero(x) for row in self.matrix for x in row)

    def is_mono(self) -> bool:
        return self.rank() == self.src.dim

    def is_epi(self) -> bool:
        return self.rank() == self.tgt.dim

    def is_iso(self) -> bool:
        return self.src.dim == self.tgt.dim and self.rank() == self.src.dim

    # -- kernel / image / cokernel as Modules-with-maps (Task 2) ------------- #
    def kernel(self):
        """``(K, iota)`` with ``K = ker(self)`` a submodule of ``src`` and
        ``iota: K -> src`` the (mono) inclusion. ``iota.then(self) == 0``."""
        from quiverlab.modules import radtopsoc
        dom = self.domain
        cols = lm.kernel_columns(self.matrix, dom) if (self.matrix and self.matrix[0]) \
            else ([lm.col(lm.identity(self.src.dim, dom), j)
                   for j in range(self.src.dim)] if self.src.dim else [])
        K = radtopsoc.submodule(self.src, cols,
                                name=f"ker({self.src.name}->{self.tgt.name})")
        iota = ModuleHom(K, self.src, _cols_to_mat(cols, self.src.dim, dom),
                         check=False)
        return K, iota

    def image(self):
        """``(I, epi, mono)``: ``I = im(self)`` as a submodule of ``tgt``, with
        ``epi: src ->> I`` and ``mono: I >-> tgt`` and ``self == epi.then(mono)``."""
        from quiverlab.modules import radtopsoc
        dom = self.domain
        pivots = lm.column_space_pivots(self.matrix, dom)
        cols = [lm.col(self.matrix, j) for j in pivots]
        I = radtopsoc.submodule(self.tgt, cols,
                                name=f"im({self.src.name}->{self.tgt.name})")
        mono = ModuleHom(I, self.tgt, _cols_to_mat(cols, self.tgt.dim, dom),
                         check=False)
        if I.dim == 0:
            epi = ModuleHom(self.src, I, lm.zeros(0, self.src.dim, dom), check=False)
        else:
            # express each column of self (image of a src basis vector) in the chosen
            # image basis (mono's columns) -- solvable because those columns span im.
            coeffs = lm.solve_columns(mono.matrix, self.matrix, dom)
            epi = ModuleHom(self.src, I,
                            [[coeffs[j][i] for j in range(self.src.dim)]
                             for i in range(I.dim)], check=False)
        return I, epi, mono

    def cokernel(self):
        """``(C, proj)`` with ``C = tgt / im(self)`` and ``proj: tgt ->> C`` (epi).
        ``self.then(proj) == 0``."""
        from quiverlab.modules.yoneda import _quotient_with_maps
        dom = self.domain
        img_cols = [lm.col(self.matrix, j) for j in range(self.src.dim)]
        C, proj_mat, _lift = _quotient_with_maps(
            self.tgt, img_cols, dom,
            name=f"coker({self.src.name}->{self.tgt.name})")
        return C, ModuleHom(self.tgt, C, proj_mat, check=False)

    def __repr__(self):
        return f"Hom({self.src.name} -> {self.tgt.name}, rank {self.rank()})"


def hom_basis(M, N):
    """Basis of ``Hom_A(M, N)`` as validated :class:`ModuleHom` objects."""
    _assert_comparable(M, N, "Hom")
    return [ModuleHom(M, N, mat, check=False) for mat in hom_space(M, N)]


def zero_hom(M, N):
    """The zero homomorphism ``M -> N``."""
    _assert_comparable(M, N, "Hom")
    z = N.domain.zero()
    return ModuleHom(M, N, [[z] * M.dim for _ in range(N.dim)], check=False)


# --------------------------------------------------------------------------- #
# k-ary direct sums with inclusion/projection maps; Krull-Schmidt summand test
# (Task 4). Block-diagonal on every shared action label -- the k-ary
# generalization of yoneda._direct_sum2, plus the obvious block inclusion/
# projection matrices (certified by the biproduct identities in the test).
# --------------------------------------------------------------------------- #
def direct_sum(*modules):
    """``D = M_1 (+) ... (+) M_k`` with the biproduct maps: returns
    ``(D, inclusions, projections)`` where ``inclusions[i]: M_i -> D`` and
    ``projections[i]: D -> M_i`` satisfy ``proj_i . incl_i = id_{M_i}`` (as
    ``incl_i.then(proj_i)``) and ``sum_i incl_i . proj_i = id_D``."""
    if not modules:
        raise QuiverlabError("direct_sum needs at least one summand")
    A = modules[0].algebra
    dom = A.domain
    side = modules[0].side
    for M in modules[1:]:
        if M.algebra is not A:
            raise QuiverlabError(
                "direct_sum: the summands are over different algebras")
        if M.side != side:
            raise QuiverlabError(
                f"direct_sum: mixing a {side} and a {M.side} summand")
    dims = [M.dim for M in modules]
    n = sum(dims)
    offs, o = [], 0
    for d in dims:
        offs.append(o)
        o += d
    action = {}
    for label in modules[0].action:
        blk = lm.zeros(n, n, dom)
        for off, M in zip(offs, modules):
            Mb = M.action.get(label)
            if Mb is None:
                raise QuiverlabError(
                    f"direct_sum: summand {M.name!r} is missing the action of "
                    f"{label!r} that the first summand carries")
            for i in range(M.dim):
                bi = blk[off + i]
                Mbi = Mb[i]
                for j in range(M.dim):
                    bi[off + j] = Mbi[j]
        action[label] = blk
    D = Module(A, n, action, name="(+)".join(M.name for M in modules), side=side)
    one = dom.one()
    incls, projs = [], []
    for off, M in zip(offs, modules):
        inc = lm.zeros(n, M.dim, dom)
        prj = lm.zeros(M.dim, n, dom)
        for i in range(M.dim):
            inc[off + i][i] = one
            prj[i][off + i] = one
        incls.append(ModuleHom(M, D, inc, check=False))
        projs.append(ModuleHom(D, M, prj, check=False))
    return D, incls, projs


def is_direct_summand(N, M, budget=512):
    """Krull-Schmidt test: ``N | M`` iff every indecomposable summand of ``N``
    appears in ``M`` with at least its multiplicity. Certified by the Plan-30
    ``decompose`` (each summand provably indecomposable, grouped by the exact
    ``is_isomorphic`` certificate); raises loudly when ``decompose`` cannot certify
    within budget -- never a silent verdict."""
    from quiverlab.modules.decompose import decompose
    from quiverlab.modules.hom import is_isomorphic
    _assert_comparable(N, M, "is_direct_summand")
    if N.dim == 0:
        return True
    dn = decompose(N, budget=budget)
    dm = [[mi, mmult] for (mi, mmult) in decompose(M, budget=budget)]
    for (ni, nmult) in dn:
        hit = None
        for k, (mi, mmult) in enumerate(dm):
            if ni.dim == mi.dim and is_isomorphic(ni, mi):
                hit = k
                break
        if hit is None or dm[hit][1] < nmult:
            return False
        dm[hit][1] -= nmult
    return True
