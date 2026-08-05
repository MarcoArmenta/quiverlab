"""Quasi-hereditary structure: standard/costandard modules, the quasi-heredity test,
good-filtration multiplicities, the characteristic tilting module and its Ringel dual
(Plan 47). Right modules; exact over any Domain. Delta/Nabla/qh/filtration are char-clean
pure linear algebra; only the tilting summand COUNT (P30 decompose) and a PRESENTED Ringel
dual (P44 presented_form) inherit the char 0 / char > dim caveat.

Order convention (pinned once, verbatim below and in every public docstring). ``order`` is
a sequence of the vertices listed **lowest -> highest**; ``rank(v)`` = its index in that
list; ``j > i`` ("j above i") means ``rank(j) > rank(i)``. ``order=None`` is the natural
order ``sorted(A.quiver.vertices)`` (ascending vertex labels). The standard module is

    Delta(i) = P(i) / trace{P(j) : rank(j) > rank(i)}

(Dlab-Ringel, Illinois J. Math. 33 (1989), 280-291). Quasi-heredity is ORDER-DEPENDENT (a
hereditary algebra is quasi-hereditary for every order -- the two-orders oracle).

Hand-derived facts (arbiters; verbatim). ``kA_n = linear_path_algebra(n)`` has forward
arrows ``a_i: i -> i+1``, RIGHT modules, so ``P(v) = e_v A`` = paths STARTING at ``v`` =
the interval module ``[v..n]``, and ``Hom(P(j), P(i)) = e_i A e_j`` = paths ``i -> j``,
nonzero iff ``i <= j``, with image the interval submodule ``[j..n]``.
  * Natural order ``1 < 2 < ... < n`` (rank(v)=v-1): ``trace{P(j):j>i} = [i+1..n] = rad P(i)``
    so ``Delta(i) = S_i``; dually ``Nabla(i) = D(Delta_{A^op}(i)) = I(i) = [1..i]``.
  * Opposite order ``n < ... < 1``: ``trace = 0`` so ``Delta(i) = P(i) = [i..n]``,
    ``Nabla(i) = S_i``.
"""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules import radtopsoc
from quiverlab.modules.builders import _require_provenance
from quiverlab.modules.hom import _assert_comparable, hom_space
from quiverlab.modules.morphism import ModuleHom, direct_sum, hom_basis


def trace_module(sources, M, name="trace"):
    """The trace submodule ``tr_{add S}(M) = sum over N in sources, f in Hom_A(N, M) of
    im(f)``: the smallest ``A``-submodule of ``M`` containing every hom image (``A``-stable
    as a sum of module-map images). Returns ``(T, iota)`` with ``iota: T >-> M`` the (mono)
    inclusion; ``T.dim == 0`` gives the zero submodule + a ``0``-column mono. Every ``N``
    must be comparable to ``M`` (same algebra/side). Char-clean pure linear algebra."""
    dom = M.domain
    cols = []
    for N in sources:
        _assert_comparable(N, M, "trace")
        for f in hom_basis(N, M):
            _I, _epi, mono = f.image()               # mono: im(f) >-> M, cols in M-coords
            for j in range(_I.dim):
                cols.append(lm.col(mono.matrix, j))
    if not cols:
        T = radtopsoc.submodule(M, [], name=name, side=M.side)
        return T, ModuleHom(T, M, lm.zeros(M.dim, 0, dom), check=False)
    piv = lm.column_space_pivots(lm.cols_to_matrix(cols), dom)
    basis = [cols[j] for j in piv]
    T = radtopsoc.submodule(M, basis, name=name, side=M.side)     # A-stable span
    iota = ModuleHom(T, M, lm.cols_to_matrix(basis), check=False)
    return T, iota


# --------------------------------------------------------------------------- #
# Task B: standard & costandard modules Delta(i) / Nabla(i).
# --------------------------------------------------------------------------- #
def _order_ranks(A, order):
    """(rank map vertex->index, order list), validated. ``order=None`` = natural order
    ``sorted(A.quiver.vertices)`` (ints before strings, each ascending)."""
    _require_provenance(A, "standard_modules")
    verts = list(A.quiver.vertices)

    def _key(v):
        return (isinstance(v, str), str(v))

    if order is None:
        order = sorted(verts, key=_key)
    order = list(order)
    if sorted(order, key=_key) != sorted(verts, key=_key):
        raise QuiverlabError(
            f"quasi-hereditary order {order!r} is not a permutation of the vertices "
            f"{verts!r}",
            hint="pass each vertex exactly once, lowest->highest; None = natural order")
    return {v: r for r, v in enumerate(order)}, order


def standard_modules(A, order=None):
    """The standard modules ``Delta(i) = P(i) / trace{P(j) : rank(j) > rank(i)}`` for the
    given ``order`` (dict vertex -> Delta(i), right ``A``-modules). Char-clean. See the
    module docstring for the order convention and the ``kA_n`` hand derivation."""
    ranks, order = _order_ranks(A, order)
    P = {v: A.projective(v) for v in order}
    out = {}
    for i in order:
        higher = [P[j] for j in order if ranks[j] > ranks[i]]
        _T, iota = trace_module(higher, P[i], name=f"tr>{i}")
        subcols = [lm.col(iota.matrix, j) for j in range(_T.dim)]
        Q = radtopsoc.quotient(P[i], subcols, name=f"Delta_{i}", side=P[i].side)
        out[i] = Q
    return out


def costandard_modules(A, order=None):
    """The costandard modules ``Nabla_A(i) = D(Delta_{A^op}(i))`` for the given ``order``
    (dict vertex -> Nabla(i), right ``A``-modules). The D-of-opposite identity IS the
    Dlab-Ringel definition of Nabla; reuses ``standard_modules`` over ``A^op`` -- no dual
    trace is re-implemented. ``dualize`` is side-aware (Plan 24): D of a right ``A^op``-mod
    is a right ``A``-mod, so ``Nabla(v) = I(v)`` / socle ``Nabla(v) = S_v`` arbitrate the
    side bookkeeping."""
    from quiverlab.modules.duality import dualize
    Aop = A.opposite()
    dstd = standard_modules(Aop, order)             # Delta over the opposite algebra
    out = {}
    for i, D in dstd.items():
        # dualize of a right A^op-module has representation algebra A and side tag "left"
        # (a left A^op-module IS a right A-module). Re-tag to "right" so Nabla(i) is a
        # right A-module comparable with A.injective(v) -- the P24 "side is a presentation
        # tag" precedent; the socle Nabla(v) = S_v / Nabla(v) = I(v) oracles arbitrate.
        nab = dualize(D).with_side("right")         # right A^op-mod -> right A-mod
        nab.name = f"Nabla_{i}"
        out[i] = nab
    return out
