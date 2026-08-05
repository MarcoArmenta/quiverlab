"""FilteredComplex -- a bounded homological Domain vector-space complex with an
increasing, exhaustive, Hausdorff filtration (Weibel *An Introduction to
Homological Algebra*, CUP 1994, section 5.4).

Homological convention (pinned, verbatim -- never silently reindexed): every
differential is ``d_n : V_n -> V_{n-1}`` with rows = target ``V_{n-1}``, cols =
source ``V_n`` -- byte-for-byte the ``modules.complexes.ChainComplex._dmats``
layout and ``engine/cyclic``'s total-differential layout. The filtration is given
per degree by an INCREASING chain of column-span bases
``F_lo V_n <= F_{lo+1} V_n <= ... <= V_n`` (filtration degree ``p = lo + j``); it
must be a subcomplex filtration (``d_n(F_p V_n) <= F_p V_{n-1}``) and exhaustive
(top level spans ``V_n``). Hausdorff-ness is the clamp: ``piece(n, p) = []`` below
``lo``. Every gate is loud (:class:`~quiverlab.errors.QuiverlabError`).

All algebra is exact over the Domain (``fields.linalg`` / ``modules.linalg_mod``);
no floats (the ``src/`` AST gate scans this file)."""
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.specseq import _subspace as sub


class FilteredComplex:
    """A bounded homological Domain vector-space complex with an increasing,
    exhaustive, Hausdorff filtration given per degree by column-span bases.

    ``terms`` : ``{n: int}`` -- ``dim V_n`` (a missing degree is the zero space).
    ``dmats`` : ``{n: matrix}`` -- ``d_n : V_n -> V_{n-1}``, rows=target
    (homological).
    ``filt``  : ``{n: list[cols]}`` -- ``filt[n][j]`` is a column-basis (list of
    coordinate columns over ``V_n``'s basis) for ``F_{lo+j} V_n``, an INCREASING
    chain ``F_lo <= F_{lo+1} <= ... <= V_n``; the filtration degree ``p = lo + j``
    (``lo`` may be negative).
    ``lo``    : ``int`` -- the least filtration degree present.

    With ``check=True``: (1) every ``d_n`` is a genuine chain differential
    (``d_{n-1} . d_n == 0``, loud); (2) each ``filt[n]`` is a nested increasing
    chain whose top level spans all of ``V_n`` (exhaustive); (3)
    ``d_n(F_p V_n) <= F_p V_{n-1}`` for every ``p`` (a subcomplex filtration),
    naming ``(n, p)`` otherwise."""

    def __init__(self, terms, dmats, filt, lo, dom, check=True):
        self.dom = dom
        self.lo = int(lo)
        self._terms = {int(n): int(d) for n, d in terms.items() if int(d) > 0}
        self._dmats = {int(n): [list(r) for r in mat] for n, mat in dmats.items()}
        self._filt = {int(n): [[list(c) for c in level] for level in chain]
                      for n, chain in filt.items()}
        if check:
            self._validate()

    # -- constructors -------------------------------------------------------- #
    @classmethod
    def from_chain_complex(cls, X, filt, lo=0, check=True):
        """Forget a P39 :class:`~quiverlab.modules.complexes.ChainComplex` ``X`` to
        its underlying ``(dims, dmats)`` vector-space complex (``term(n).dim`` +
        ``differential(n).matrix``) and attach the caller's per-degree ``filt``
        whose common origin is ``lo`` (``filt[n][j]`` spans ``F_{lo+j} V_n``). The
        module structure is dropped -- homology is unchanged. ``lo`` defaults to
        ``0`` (the one-step / stupid filtrations); the radical-filtration preset
        passes the negative origin of its increasing radical chain."""
        dims = {n: X.term(n).dim for n in X.degrees()}
        dmats = {}
        for n in X.degrees():
            if X.term(n - 1).dim and X.term(n).dim:
                mat = X.differential(n).matrix
                if mat and mat[0]:
                    dmats[n] = mat
        return cls(dims, dmats, filt, lo, X.domain, check=check)

    # -- accessors ----------------------------------------------------------- #
    def degrees(self):
        return sorted(self._terms)

    def dim(self, n):
        return self._terms.get(int(n), 0)

    def dmat(self, n):
        n = int(n)
        mat = self._dmats.get(n)
        if mat is not None and mat and mat[0]:
            return mat
        return lm.zeros(self.dim(n - 1), self.dim(n), self.dom)

    def _top_level(self, n):
        chain = self._filt.get(int(n), [])
        return self.lo + len(chain) - 1

    def levels(self):
        """Filtration degrees ``p`` present, ascending (``lo`` .. global top)."""
        if not self._terms:
            return [self.lo]
        hi = max((self._top_level(n) for n in self._terms), default=self.lo)
        return list(range(self.lo, max(hi, self.lo) + 1))

    def piece(self, n, p):
        """Column-basis of ``F_p V_n``: below ``lo`` -> ``[]`` (Hausdorff clamp);
        at/above the top declared filtration degree -> the raw top level
        ``filt[n][-1]`` (the filtration stabilizes there; it spans ``V_n`` exactly
        when exhaustive -- so page/closedness algebra reads spans honestly, never
        the assumed-full identity); otherwise ``filt[n][p - lo]``. A term with no
        recorded filtration defaults to the whole space (the trivial filtration)."""
        n, p = int(n), int(p)
        if p < self.lo:
            return []
        chain = self._filt.get(n, [])
        if not chain:
            d = self.dim(n)
            ident = lm.identity(d, self.dom)
            return [lm.col(ident, j) for j in range(d)]
        top = self._top_level(n)
        if p >= top:
            return [list(c) for c in chain[-1]]
        return [list(c) for c in chain[p - self.lo]]

    # -- validation ---------------------------------------------------------- #
    def _validate(self):
        dom = self.dom
        # (1) d . d == 0
        for n in sorted(self._dmats):
            dn, dn1 = self._dmats.get(n), self._dmats.get(n + 1)
            if dn is None or dn1 is None or not (dn and dn[0]) or not (dn1 and dn1[0]):
                continue
            comp = lm.matmul(dn, dn1, dom)
            if comp and any(not dom.is_zero(x) for row in comp for x in row):
                raise QuiverlabError(
                    f"FilteredComplex: d_{n} . d_{n + 1} != 0 -- the differentials "
                    "do not form a chain complex",
                    hint="d_n : V_n -> V_{n-1}, rows=target; check the layout")
        # (2) each filt[n] a nested increasing chain
        for n, chain in self._filt.items():
            prev = []
            for j, level in enumerate(chain):
                if not sub.span_contains(prev, level, dom):
                    raise QuiverlabError(
                        f"FilteredComplex: filtration at degree {n} is not "
                        f"increasing -- level {j} does not contain level {j - 1}")
                prev = level
        # (3) closed: d_n(F_p V_n) <= F_p V_{n-1}  (checked BEFORE exhaustiveness so a
        #     filtration that both escapes d AND under-spans reports the subcomplex bug)
        levels = self.levels()
        for n in self.degrees():
            dn = self._dmats.get(n)
            if dn is None or not (dn and dn[0]):
                continue
            for p in levels:
                src = self.piece(n, p)
                if not src:
                    continue
                tgt = self.piece(n - 1, p)
                for c in src:
                    img = lm.matvec(dn, c, dom)
                    if not sub.in_span(img, tgt, dom):
                        raise QuiverlabError(
                            f"FilteredComplex: filtration is NOT closed under d at "
                            f"(n={n}, p={p}): d_{n}(F_{p} V_{n}) escapes "
                            f"F_{p} V_{n - 1} -- not a subcomplex filtration")
        # (4) exhaustive: the top DECLARED level (raw, not the piece() full-clamp)
        #     spans V_n
        for n, d in self._terms.items():
            chain = self._filt.get(n, [])
            top_level = chain[-1] if chain else []
            if sub.span_dim(top_level, dom) != d:
                raise QuiverlabError(
                    f"FilteredComplex: filtration at degree {n} is not exhaustive "
                    f"-- its top level does not span V_{n} (dim {d})")

    # -- homology ------------------------------------------------------------ #
    def total_homology_dims(self):
        """``{n: dim H_n(V_.)}`` over the support range, via the rank formula
        ``dim H_n = dim V_n - rank d_n - rank d_{n+1}``."""
        dom = self.dom
        degs = self.degrees()
        if not degs:
            return {}
        out = {}
        for n in range(degs[0], degs[-1] + 1):
            vn = self.dim(n)
            dn, dn1 = self._dmats.get(n), self._dmats.get(n + 1)
            rn = lm.mat_rank(dn, dom) if (dn and dn[0]) else 0
            rn1 = lm.mat_rank(dn1, dom) if (dn1 and dn1[0]) else 0
            out[n] = vn - rn - rn1
        return out

    def __repr__(self):
        degs = self.degrees()
        body = " <- ".join(f"V_{n}({self.dim(n)})" for n in reversed(degs))
        return f"FilteredComplex[lo={self.lo}]: {body}"
