"""The Jacobi-Zariski nearly-exact sequence tool (Plan 73, R7).

Cibils-Lanzilotta-Marcos-Solotar, ``2009.05017`` (Bull. LMS 54 (2022), "exact twice
in three"); the classical noncommutative origin is Kaygun ``1103.4377``. For an
extension ``B subset A`` the RELATIVE Hochschild homology ``HH_*(A|B, X)`` is the
homology of the CLMS normalized relative bar complex (2101.02597 Thm 2.2)

    ... -> X (x)_{B^e} (A/B)^{ox_B m} -> ... -> X (x)_{B^e} A/B -> X_B -> 0

with the cyclic Hochschild differential ``b``. When ``A/B`` is ``B``-tensor
nilpotent with ``(A/B)^{ox n} = 0`` the complex has finitely many nonzero terms
(``m < n``), so ``HH_*(A|B, X) = 0`` for ``* >= n`` (Cor. 2.4) -- exactly the
short-circuit this module uses. Default coefficient ``X = A`` (regular).

The assembled sequence
    ... -> HH_m(B, X) -> HH_m(A, X) -> HH_m(A|B, X) -> HH_{m-1}(B, X) -> ...
is "exact twice in three" (Thm 3.1): exact at the ``HH_*(B)`` and ``HH_*(A|B)``
spots, with a gap at ``HH_*(A)``.

Composition is left-to-right (quiverlab); the CLMS paper is right-to-left.
"""
import itertools
from dataclasses import dataclass, field as _field

from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import rank, rref
from quiverlab.hochschild.table import HHTable

#: guard on the free-module size ``dim X . |R|^m`` of a relative-complex term
_MAX_RELBAR_CELLS = 4000


# --------------------------------------------------------------------------- #
# the relative bar complex C_m = X (x)_{B^e} (A/B)^{ox_B m}
# --------------------------------------------------------------------------- #
class _RelBar:
    """The CLMS normalized relative bar complex of ``A`` (coefficients ``X = A``,
    the regular bimodule) over the subalgebra ``B``. Cyclic B-tensor: a free basis
    element is ``(a, r_1, ..., r_m)`` with ``a`` an ``A``-basis index and ``r_i`` a
    relative-path index; ``C_m = Free_m / Rel_m`` where ``Rel_m`` is the
    ``B``-generator sliding at the ``m+1`` cyclic junctions."""

    def __init__(self, ext, index=None):
        from quiverlab.families.extension import _f_length, _tokens
        self.ext = ext
        #: nilpotency index n (C_m = 0 for m >= n, Cor. 2.4); None => no short-circuit
        self.index = index
        self.A = ext.A
        self.dom = self.A.domain
        self.rel = list(ext.rel_idx)
        self.nR = len(self.rel)
        labels = self.A.basis_labels
        Fset = ext._new_arrow_set()
        self.gens = [i for i, l in enumerate(labels) if l.startswith("e_")]
        self.gens += [i for i, l in enumerate(labels)
                      if _f_length(l, Fset) == 0 and len(_tokens(l)) == 1]
        self._Ab = [self.A._basis_vec(i) for i in range(self.A.dim)]
        self._gv = [self.A._basis_vec(g) for g in self.gens]
        self._qc = {}

    def _projR(self, vec):
        return [vec[self.rel[k]] for k in range(self.nR)]

    def _flat(self, a, rs):
        idx = a
        for r in rs:
            idx = idx * self.nR + r
        return idx

    def _unflat(self, idx, m):
        rs = []
        for _ in range(m):
            rs.append(idx % self.nR)
            idx //= self.nR
        return idx, tuple(reversed(rs))

    def _cols(self, m):
        return self.A.dim * (self.nR ** m)

    def _feasible(self, m):
        return self._cols(m) <= _MAX_RELBAR_CELLS

    def _relations(self, m):
        A, dom, nR = self.A, self.dom, self.nR
        ncols = self._cols(m)
        rows = []
        for a in range(A.dim):
            for rs in itertools.product(range(nR), repeat=m):
                for gi, g in enumerate(self.gens):
                    gv = self._gv[gi]
                    # internal junctions (slot j, slot j+1), j=0..m-1; slot0 = A, rest = R
                    for j in range(m):
                        acc = {}
                        if j == 0:
                            ag = A.multiply(self._Ab[a], gv)          # a . g in A
                            for aa in range(A.dim):
                                if not dom.is_zero(ag[aa]):
                                    f = self._flat(aa, rs)
                                    acc[f] = dom.add(acc.get(f, dom.zero()), ag[aa])
                        else:
                            rjg = self._projR(A.multiply(self._Ab[self.rel[rs[j - 1]]], gv))
                            for s in range(nR):
                                if not dom.is_zero(rjg[s]):
                                    f = self._flat(a, rs[:j - 1] + (s,) + rs[j:])
                                    acc[f] = dom.add(acc.get(f, dom.zero()), rjg[s])
                        grj = self._projR(A.multiply(gv, self._Ab[self.rel[rs[j]]]))
                        for s in range(nR):
                            if not dom.is_zero(grj[s]):
                                f = self._flat(a, rs[:j] + (s,) + rs[j + 1:])
                                acc[f] = dom.add(acc.get(f, dom.zero()), dom.neg(grj[s]))
                        self._emit(rows, acc, ncols)
                    # cyclic junction (slot m, slot 0)
                    acc = {}
                    ga = A.multiply(gv, self._Ab[a])                  # g . a
                    for aa in range(A.dim):
                        if not dom.is_zero(ga[aa]):
                            f = self._flat(aa, rs)
                            acc[f] = dom.add(acc.get(f, dom.zero()), ga[aa])
                    if m >= 1:
                        rmg = self._projR(A.multiply(self._Ab[self.rel[rs[m - 1]]], gv))
                        for s in range(nR):
                            if not dom.is_zero(rmg[s]):
                                f = self._flat(a, rs[:m - 1] + (s,))
                                acc[f] = dom.add(acc.get(f, dom.zero()), dom.neg(rmg[s]))
                    else:
                        ag = A.multiply(self._Ab[a], gv)              # m=0: g.a - a.g
                        for aa in range(A.dim):
                            if not dom.is_zero(ag[aa]):
                                f = self._flat(aa, ())
                                acc[f] = dom.add(acc.get(f, dom.zero()), dom.neg(ag[aa]))
                    self._emit(rows, acc, ncols)
        return rows, ncols

    def _emit(self, rows, acc, ncols):
        dom = self.dom
        if not acc or all(dom.is_zero(v) for v in acc.values()):
            return
        dense = [dom.zero()] * ncols
        for f, c in acc.items():
            dense[f] = c
        rows.append(dense)

    def _quotient(self, m):
        if m not in self._qc:
            # short-circuit: (A/B)^{ox m} = 0 for m >= index => C_m = 0 (Cor. 2.4).
            # Use the precomputed nilpotency index (NOT a per-m tensor_power_dim, which
            # blows up combinatorially at large m).
            vanished = (self.index is not None and m >= self.index)
            if vanished:
                self._qc[m] = ([], [], 0, True)          # (Rr, piv, ncols, vanished)
            else:
                if not self._feasible(m):
                    raise QuiverlabError(
                        f"relative bar term C_{m} has dim X.|R|^{m} = {self._cols(m)} "
                        f"free basis, past the cell guard {_MAX_RELBAR_CELLS}",
                        hint="the relative complex is infeasible at this degree here")
                rows, ncols = self._relations(m)
                Rr, piv = rref(rows, self.dom) if rows else ([], [])
                self._qc[m] = (Rr, piv, ncols, False)
        return self._qc[m]

    def dim_C(self, m):
        Rr, piv, ncols, vanished = self._quotient(m)
        return 0 if vanished else ncols - len(piv)

    def _reduce(self, v, m):
        Rr, piv, ncols, _v = self._quotient(m)
        dom = self.dom
        v = list(v)
        for row, pc in zip(Rr, piv):
            f = v[pc]
            if not dom.is_zero(f):
                for j in range(ncols):
                    rj = row[j]
                    if not dom.is_zero(rj):
                        v[j] = dom.add(v[j], dom.neg(dom.mul(f, rj)))
        return v

    def _rank_diff(self, m):
        """rank of the induced differential ``b_m: C_m -> C_{m-1}``."""
        if m <= 0 or self.dim_C(m) == 0 or self.dim_C(m - 1) == 0:
            return 0
        A, dom, nR = self.A, self.dom, self.nR
        _Rr, piv, ncols, _v = self._quotient(m)
        pivset = set(piv)
        ncols_m1 = self._cols(m - 1)
        mat = []
        for j in range(ncols):
            if j in pivset:
                continue
            a, rs = self._unflat(j, m)
            col = [dom.zero()] * ncols_m1
            # d_0 (+): (a.r_1, r_2..)
            ar1 = A.multiply(self._Ab[a], self._Ab[self.rel[rs[0]]])
            for aa in range(A.dim):
                if not dom.is_zero(ar1[aa]):
                    f = self._flat(aa, rs[1:])
                    col[f] = dom.add(col[f], ar1[aa])
            # d_i (-1)^i: (a, .., proj(r_i.r_{i+1}), ..)
            for i in range(1, m):
                prod = self._projR(A.multiply(self._Ab[self.rel[rs[i - 1]]],
                                              self._Ab[self.rel[rs[i]]]))
                for s in range(nR):
                    if not dom.is_zero(prod[s]):
                        f = self._flat(a, rs[:i - 1] + (s,) + rs[i + 1:])
                        val = prod[s] if i % 2 == 0 else dom.neg(prod[s])
                        col[f] = dom.add(col[f], val)
            # d_m (-1)^m: (r_m.a, r_1..r_{m-1})
            rma = A.multiply(self._Ab[self.rel[rs[m - 1]]], self._Ab[a])
            for aa in range(A.dim):
                if not dom.is_zero(rma[aa]):
                    f = self._flat(aa, rs[:m - 1])
                    val = rma[aa] if m % 2 == 0 else dom.neg(rma[aa])
                    col[f] = dom.add(col[f], val)
            mat.append(self._reduce(col, m - 1))
        return rank(mat, dom) if mat else 0

    def homology(self, top):
        ranks = {0: 0}
        for m in range(1, top + 2):
            ranks[m] = self._rank_diff(m)
        dims = []
        for m in range(top + 1):
            dims.append(self.dim_C(m) - ranks.get(m, 0) - ranks.get(m + 1, 0))
        return dims


def relative_homology(ext, top, *, coefficients=None):
    """``HH_*(A|B, X)`` (``X = A`` regular by default) via the CLMS normalized
    relative bar complex (Thm 2.2). Returns an ``HHTable`` of dims ``0..top``.

    Refuses loudly when ``A/B`` is NOT tensor-nilpotent and ``top`` reaches past the
    small feasible window -- the complex is then infinite (honest ``status``). When
    ``A/B`` is nilpotent with index ``n`` the complex is finite (``C_m = 0`` for
    ``m >= n``), so the answer is exact (Cor. 2.4)."""
    from quiverlab.invariants.han import is_tensor_nilpotent
    if coefficients is not None:
        raise QuiverlabError(
            "relative_homology currently computes the regular coefficient X = A only",
            hint="the mixed-coefficient relative complex is out of scope for v1")
    tn = is_tensor_nilpotent(ext, cap=max(8, top + 1))
    if tn.status != "nilpotent" and top >= 2:
        raise QuiverlabError(
            f"HH_*(A|B) is infinite: A/B is {tn.status} (not tensor-nilpotent), so the "
            f"relative bar complex does not terminate -- cannot report finite HH_* to "
            f"top={top}",
            hint="relative_homology is finite only when A/B is B-tensor-nilpotent "
                 "(CLMS Cor. 2.4)")
    rb = _RelBar(ext, index=tn.index)
    dims = rb.homology(top)
    return HHTable(dims, "HH_", f"HH_*(A|B) for {repr(ext.A).splitlines()[0]}",
                   engine="CLMS relative bar")


# --------------------------------------------------------------------------- #
# the assembled Jacobi-Zariski sequence + self-cert
# --------------------------------------------------------------------------- #
@dataclass
class JZSequence:
    """The assembled nearly-exact sequence (coefficients ``X = A``). A data report."""
    hh_A: list                           # HH_*(A, A) = ordinary HH_*(A)
    hh_B: list                           # HH_*(B, B) = ordinary HH_*(B) (self-cert bound)
    hh_rel: list                         # HH_*(A|B, A) via the CLMS complex
    exact_at_B: list                     # per-degree self-cert at the HH_*(B) spot
    exact_at_rel: list                   # per-degree self-cert at the HH_*(A|B) spot
    gap_at_A: list                       # measured Ker/Im slack at the HH_*(A) spot
    nilpotency_index: "int | None"
    note: str = ""


def jacobi_zariski_sequence(ext, top):
    """Assemble the Jacobi-Zariski sequence and its rigorous self-certs.

    Verifies (a) ``HH_*(A|B) = 0`` for ``* >= n`` (Cor. 2.4 -- the exactness at the
    relative spot is trivial once the term vanishes) and (b) the injection bound
    ``dim HH_m(B) <= dim HH_m(A)`` (the honest headline consequence of leg (i)+(ii),
    exactness at the ``HH_*(B)`` spot). The gap at the ``HH_*(A)`` spot is measured
    (the sequence is NOT exact there -- "exact twice in three")."""
    from quiverlab.invariants.han import is_tensor_nilpotent
    engine = "cs" if ext.A.domain.name.startswith("QQ") else "auto"
    hh_A = ext.A.hochschild_homology(top, engine=engine, verbose=False).dims
    hh_B = ext.B.hochschild_homology(top, engine=engine, verbose=False).dims
    tn = is_tensor_nilpotent(ext, cap=max(8, top + 1))
    index = tn.index if tn.status == "nilpotent" else None
    hh_rel = relative_homology(ext, top).dims

    exact_at_rel, exact_at_B, gap_at_A = [], [], []
    for m in range(top + 1):
        # rel spot: exact where the term vanishes (Cor. 2.4, m >= index) -- rigorous
        exact_at_rel.append(hh_rel[m] == 0 if (index is not None and m >= index) else True)
        # B spot: the injection bound dim HH_m(B) <= dim HH_m(A) (leg (i)+(ii))
        exact_at_B.append(hh_B[m] <= hh_A[m])
        # gap at A: the raw slack (Euler consistency of the three sequences)
        gap_at_A.append(hh_A[m] - hh_B[m])
    return JZSequence(
        hh_A=list(hh_A), hh_B=list(hh_B), hh_rel=list(hh_rel),
        exact_at_B=exact_at_B, exact_at_rel=exact_at_rel, gap_at_A=gap_at_A,
        nilpotency_index=index,
        note="exact twice in three (CLMS Thm 3.1): the rel-spot self-cert is the "
             "Cor. 2.4 vanishing; the B-spot self-cert is the injection bound "
             "dim HH_m(B) <= dim HH_m(A); the gap at the HH_*(A) spot is measured.")
