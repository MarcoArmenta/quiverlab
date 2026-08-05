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


# --------------------------------------------------------------------------- #
# Task C: is_quasi_hereditary (QHReport) + delta_multiplicities + BGG reciprocity.
# --------------------------------------------------------------------------- #
@dataclass
class QHReport:
    """Three-valued (like ``TiltingReport``): ``bool(report)`` is
    ``is_quasi_hereditary``; ``repr`` names the failing clause. ``per_index`` records, per
    vertex, the ``brick`` (``End Delta(i) = k``) and ``delta_filters_P`` certificates."""
    is_quasi_hereditary: bool
    order: list
    gl_dim: object          # GlobalDimension (qh => finite, Dlab-Ringel)
    per_index: dict         # vertex -> {"brick": bool, "delta_filters_P": bool, "note": str}
    note: str               # names the failing clause, "quasi-hereditary" when ok

    def __bool__(self):
        return self.is_quasi_hereditary

    def __repr__(self):
        if self.is_quasi_hereditary:
            return f"quasi-hereditary (order {self.order})"
        return f"not quasi-hereditary: {self.note}"


def _hom_dim(M, N):
    return len(hom_space(M, N))


def _combine_homs(homs, coeffs, dom):
    """The dom-linear combination ``sum coeffs[k] * homs[k].matrix`` (a raw matrix)."""
    nrow = len(homs[0].matrix)
    ncol = len(homs[0].matrix[0]) if homs[0].matrix else 0
    out = lm.zeros(nrow, ncol, dom)
    for c, f in zip(coeffs, homs):
        if dom.is_zero(c):
            continue
        for i in range(nrow):
            oi, fi = out[i], f.matrix[i]
            for j in range(ncol):
                oi[j] = dom.add(oi[j], dom.mul(c, fi[j]))
    return out


def _find_surjection(M, D):
    """An epi ``M ->> D`` as a raw matrix, or ``None``. A single basis hom already surjects
    for the directed/uniserial oracles; the general case is a bounded search over small
    field-coefficient combinations of ``hom_basis(M, D)`` (never fabricate a surjection --
    ``None`` is the honest 'cannot peel Delta here')."""
    dom = M.domain
    homs = hom_basis(M, D)
    for f in homs:
        if f.is_epi():
            return f.matrix
    if not homs or D.dim == 0:
        return None
    # bounded combination search: units first, then a small integer ladder.
    import itertools
    from quiverlab.fields.primefield import PrimeField
    r = len(homs)
    if isinstance(dom, PrimeField):
        vals = [dom.coerce(i) for i in range(dom.p)]
    else:
        vals = [dom.coerce(i) for i in range(-1, 3)]
    cap = 20000
    if len(vals) ** r > cap:
        return None                                  # too big to search honestly
    for idx in itertools.product(range(len(vals)), repeat=r):
        if all(t == 0 for t in idx):
            continue
        C = _combine_homs(homs, [vals[t] for t in idx], dom)
        if lm.mat_rank(C, dom) == D.dim:
            return C
    return None


def delta_multiplicities(M, deltas, order=None):
    """Greedy top-down Delta-peel: at each step take the HIGHEST-in-``order`` ``Delta(i)``
    whose top ``S_i`` appears in ``top(M)`` and which admits an epi ``M ->> Delta(i)``;
    quotient by its kernel; recurse. Returns ``(mult: vertex -> (M:Delta(i)), certified)``.
    ``certified=False`` (LOUD) when the peel stalls with ``M != 0`` (no Delta-filtration) --
    three-valued honesty like ``TiltingReport``; never a fabricated count."""
    if not deltas:
        return {}, M.dim == 0
    A = deltas[next(iter(deltas))].base_algebra
    ranks, order = _order_ranks(A, order)
    mult = {i: 0 for i in deltas}
    cur = M
    guard = 0
    while cur.dim > 0:
        guard += 1
        if guard > M.dim + 1:                        # cannot exceed dim M genuine layers
            return mult, False
        tops = cur.top().dimension_vector()
        cands = sorted((i for i in deltas if tops.get(i, 0) > 0),
                       key=lambda i: ranks.get(i, -1), reverse=True)
        peeled = False
        for i in cands:
            fmat = _find_surjection(cur, deltas[i])
            if fmat is None:
                continue
            epi = ModuleHom(cur, deltas[i], fmat, check=False)
            K, _iota = epi.kernel()
            mult[i] += 1
            cur = K
            peeled = True
            break
        if not peeled:
            return mult, False                        # uncertified: no Delta peels off the top
    return mult, True


def is_quasi_hereditary(A, order=None):
    """A :class:`QHReport`. ``A`` with order ``>`` is quasi-hereditary iff for every ``i``
    (1) ``End_A(Delta(i)) = k`` (Delta(i) is a brick) and (2) ``P(i)`` has a Delta-filtration
    (top ``Delta(i)``, rest ``Delta(j)``, ``j > i``); a NECESSARY classical consequence is
    ``gl.dim A < infinity`` (Dlab-Ringel). Three-valued: ``True``, or ``False`` with a
    ``note`` naming the first failing clause. Char-clean pure linear algebra."""
    from quiverlab.modules.ext import global_dimension
    ranks, order = _order_ranks(A, order)
    deltas = standard_modules(A, order)
    gld = global_dimension(A)                          # qh => finite (Dlab-Ringel)
    per, ok = {}, True
    for i in order:
        brick = (_hom_dim(deltas[i], deltas[i]) == 1)   # End Delta(i) = k
        _mult, filt = delta_multiplicities(A.projective(i), deltas, order)
        note = ("ok" if (brick and filt) else
                ("End Delta != k" if not brick else "P has no Delta-filtration"))
        per[i] = {"brick": brick, "delta_filters_P": filt, "note": note}
        ok = ok and brick and filt
    ok = ok and gld.exact
    if ok:
        note = "quasi-hereditary"
    elif not gld.exact:
        note = f"gl.dim not finite ({gld!r})"
    else:
        note = "; ".join(f"{i}: {per[i]['note']}" for i in order if per[i]["note"] != "ok")
    return QHReport(ok, order, gld, per, note)
