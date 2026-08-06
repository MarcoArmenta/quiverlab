"""DoubleComplex -- a bounded homological double complex over any Domain, its
total complex, and the row / column filtrations that feed the page engine.

Homological throughout: ``d_h : D_{p,q} -> D_{p-1,q}`` and
``d_v : D_{p,q} -> D_{p,q-1}`` (rows = target). The sign convention is the one
that makes the total differential square to zero -- the two paths of every square
ANTICOMMUTE:

    d_h[(p,q-1)] . d_v[(p,q)]  +  d_v[(p-1,q)] . d_h[(p,q)]  ==  0

(a strictly commuting bicomplex must be sign-adjusted by the caller -- the standard
``s_{p,q} = (-1)^p`` on one differential -- BEFORE construction; the Hochschild
``(b, B)`` preset's differentials already anticommute by the mixed-complex identity
``bB + Bb = 0``, and the Grothendieck Hom double complex carries the Koszul sign).

``total()`` assembles ``Tot_n = (+)_{p+q=n} D_{p,q}`` with the block-place layout of
``engine/cyclic._total_differential``; both filtrations are genuine subcomplex
filtrations (``FilteredComplex`` re-validates closed + exhaustive). No floats
(``src/`` AST gate)."""
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.specseq.filtered import FilteredComplex


def _nonzero(mat, dom):
    return bool(mat) and any(not dom.is_zero(x) for row in mat for x in row)


class DoubleComplex:
    """A bounded homological double complex over a Domain.

    ``terms`` : ``{(p, q): int}`` -- ``dim D_{p,q}`` (a missing entry is zero).
    ``d_h``   : ``{(p, q): matrix}`` -- ``D_{p,q} -> D_{p-1,q}``, rows=target.
    ``d_v``   : ``{(p, q): matrix}`` -- ``D_{p,q} -> D_{p,q-1}``, rows=target.

    ``check=True`` asserts ``d_h . d_h == 0``, ``d_v . d_v == 0`` and the
    anticommutation at every ``(p, q)`` -- ``QuiverlabError`` naming ``(p, q)``."""

    def __init__(self, terms, d_h, d_v, dom, check=True):
        self.dom = dom
        self._terms = {(int(p), int(q)): int(d)
                       for (p, q), d in terms.items() if int(d) > 0}
        self._d_h = {(int(p), int(q)): [list(r) for r in mat]
                     for (p, q), mat in d_h.items() if mat and mat[0]}
        self._d_v = {(int(p), int(q)): [list(r) for r in mat]
                     for (p, q), mat in d_v.items() if mat and mat[0]}
        if check:
            self._validate()

    def dim(self, p, q):
        return self._terms.get((int(p), int(q)), 0)

    # -- validation ---------------------------------------------------------- #
    def _validate(self):
        dom = self.dom
        for (p, q) in self._terms:
            # d_h . d_h : D_{p,q} -> D_{p-2,q}
            a, b = self._d_h.get((p - 1, q)), self._d_h.get((p, q))
            if a and b and _nonzero(lm.matmul(a, b, dom), dom):
                raise QuiverlabError(
                    f"DoubleComplex: d_h . d_h != 0 at (p={p}, q={q})")
            # d_v . d_v : D_{p,q} -> D_{p,q-2}
            a, b = self._d_v.get((p, q - 1)), self._d_v.get((p, q))
            if a and b and _nonzero(lm.matmul(a, b, dom), dom):
                raise QuiverlabError(
                    f"DoubleComplex: d_v . d_v != 0 at (p={p}, q={q})")
            # anticommutation over the target D_{p-1,q-1}
            rows = self._terms.get((p - 1, q - 1), 0)
            cols = self._terms.get((p, q), 0)
            if not (rows and cols):
                continue
            t1h, t1v = self._d_h.get((p, q - 1)), self._d_v.get((p, q))
            t2v, t2h = self._d_v.get((p - 1, q)), self._d_h.get((p, q))
            S = lm.zeros(rows, cols, dom)
            for term in (lm.matmul(t1h, t1v, dom) if (t1h and t1v) else None,
                         lm.matmul(t2v, t2h, dom) if (t2v and t2h) else None):
                if not term:
                    continue
                for i in range(rows):
                    Si, Ti = S[i], term[i]
                    for j in range(cols):
                        Si[j] = dom.add(Si[j], Ti[j])
            if _nonzero(S, dom):
                raise QuiverlabError(
                    "DoubleComplex: d_h.d_v + d_v.d_h != 0 (anticommutation "
                    f"broken) at (p={p}, q={q}) -- the total differential would "
                    "not square to zero; sign-adjust one differential first")

    # -- total complex ------------------------------------------------------- #
    def _layout(self):
        """``(order, offset, terms_tot)``: per total degree ``n`` the block list
        (descending ``p``), each block's row/col offset within ``Tot_{p+q}``, and
        ``dim Tot_n``."""
        by_n = {}
        for (p, q) in self._terms:
            by_n.setdefault(p + q, []).append((p, q))
        order, offset, terms_tot = {}, {}, {}
        for n, blocks in by_n.items():
            blocks.sort(key=lambda pq: (-pq[0], pq[1]))   # descending p
            order[n] = blocks
            off = 0
            for pq in blocks:
                offset[pq] = off
                off += self._terms[pq]
            terms_tot[n] = off
        return order, offset, terms_tot

    def total(self):
        """``(terms, dmats, dom)``: ``Tot_n = (+)_{p+q=n} D_{p,q}`` with
        ``D = d_h + d_v : Tot_n -> Tot_{n-1}``, each block placed at its target's
        row-offset (the ``engine/cyclic`` layout). ``d.d = 0`` holds because the
        squares anticommute (asserted at construction)."""
        dom = self.dom
        order, offset, terms_tot = self._layout()
        dmats = {}
        for n, blocks in order.items():
            nrows = terms_tot.get(n - 1, 0)
            ncols = terms_tot[n]
            if not (nrows and ncols):
                continue
            D = lm.zeros(nrows, ncols, dom)
            for (p, q) in blocks:
                c0 = offset[(p, q)]
                for (tp, tq), blk in (((p - 1, q), self._d_h.get((p, q))),
                                      ((p, q - 1), self._d_v.get((p, q)))):
                    if blk is None or (tp, tq) not in self._terms:
                        continue
                    r0 = offset[(tp, tq)]
                    for i, rowv in enumerate(blk):
                        Dr = D[r0 + i]
                        for j, val in enumerate(rowv):
                            if not dom.is_zero(val):
                                Dr[c0 + j] = dom.add(Dr[c0 + j], val)
            dmats[n] = D
        return dict(terms_tot), dmats, dom

    # -- filtrations --------------------------------------------------------- #
    def _filtration(self, by_first):
        """Shared builder: ``by_first=True`` filters by ``p`` (column filtration),
        ``False`` by ``q`` (row filtration). ``F_j Tot_n`` = the identity columns of
        the blocks whose filtering index ``<= j``."""
        dom = self.dom
        order, offset, terms_tot = self._layout()
        keys = [(p if by_first else q) for (p, q) in self._terms]
        lo, hi = (min(keys), max(keys)) if keys else (0, 0)
        filt = {}
        for n, blocks in order.items():
            d = terms_tot[n]
            ident = lm.identity(d, dom)
            levels = []
            for j in range(lo, hi + 1):
                cols = []
                for (p, q) in blocks:
                    if (p if by_first else q) <= j:
                        base = offset[(p, q)]
                        cols.extend(lm.col(ident, base + k)
                                    for k in range(self._terms[(p, q)]))
                levels.append(cols)
            filt[n] = levels
        terms_t, dmats_t, _ = self.total()
        return FilteredComplex(terms_t, dmats_t, filt, lo, dom, check=True)

    def column_filtration(self):
        """``F_j Tot_n = (+)_{p<=j} D_{p,n-p}`` -- ``d_h`` lowers ``p``, ``d_v``
        preserves it, so ``d(F_j) <= F_j`` (a subcomplex filtration)."""
        return self._filtration(by_first=True)

    def row_filtration(self):
        """``F_j Tot_n = (+)_{q<=j} D_{n-q,q}`` -- ``d_v`` lowers ``q``, ``d_h``
        preserves it (a subcomplex filtration)."""
        return self._filtration(by_first=False)

    def __repr__(self):
        return (f"DoubleComplex[{len(self._terms)} blocks over {self.dom.name}]")
