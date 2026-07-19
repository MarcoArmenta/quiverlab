"""CSResolution: an engine.resolutions.Resolution backend backed by the
Chouhy-Solotar construction over GF(p).

The hanlab engine's Resolution contract wants ordered hashable term generators and
int64 differential matrices expressed in exactly that ordering (homology side:
shape (dim_{n-1}, dim_n); cohomology side: (dim_{n+1}, dim_n)). This facade wraps a
ChouhySolotarResolution (built over the algebra's PrimeField domain) and packs its
domain-generic matrices into numpy.int64. The CS arithmetic is native GF(p), so the
packed entries are the field representatives in [0, p) (there is no wider-ring lift
to keep un-reduced); the matrices are the raw CS differentials, never row-reduced.

The generators are the (Chain, corner-index) pairs of ChouhySolotarResolution._basis
(Chain is a frozen dataclass -> hashable), so building {g: i} indices works directly."""
import numpy as np

from quiverlab.engine.resolutions import Resolution


class CSResolution(Resolution):
    """Test-only surface: gate-free by design. It enforces admissibility and prime-field
    scope but does NOT run the two binding gates (assert_dd_zero / assert_order_condition)
    that the public path in homology.py runs; its certificate is the cross-check batteries
    instead. The public Algebra dispatch always routes to the gated homology.py, never here."""

    def __init__(self, A, max_cells=4_000_000):
        from quiverlab.errors import FieldError
        from quiverlab.fields.primefield import PrimeField
        if not isinstance(A.domain, PrimeField):
            raise FieldError(
                f"the CS int64 engine facade computes over prime fields only; this "
                f"algebra is over {A.domain.name}",
                hint="use the domain-generic cs_(co)homology_dims for any field, or "
                     "construct the algebra over GF(p)")
        self.A = A
        self.max_cells = max_cells
        self._res = None
        self._top = -1

    def _cs(self, need):
        """A ChouhySolotarResolution with max_degree >= need (grown on demand)."""
        from quiverlab.resolutions_cs.build import reduction_system_of
        from quiverlab.resolutions_cs.homology import _require_admissible
        from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
        if self._res is None or self._top < need:
            rs = reduction_system_of(self.A); _require_admissible(rs)
            self._res = ChouhySolotarResolution(self.A, rs, max_degree=need,
                                                max_cells=self.max_cells)
            self._top = need
        return self._res

    @staticmethod
    def _pack(cs, M, rows_basis, index_other, cols_basis):
        col_of = {g: c for c, g in enumerate(cols_basis)}
        out = np.zeros((len(index_other), len(cols_basis)), dtype=np.int64)
        for cj, g in enumerate(cols_basis):
            c = col_of[g]
            for ri, rg in enumerate(rows_basis):
                r = index_other.get(rg)
                if r is not None:
                    out[r, cj] = int(M[ri][c])
        return out

    # -- homology side -----------------------------------------------------------
    def term_basis(self, alg, n):
        return list(self._cs(n)._basis(n, "hom"))

    def differential_matrix(self, alg, n, basis_n, index_nm1):
        cs = self._cs(n)
        return self._pack(cs, cs.matrix(n, "hom"), cs._basis(n - 1, "hom"),
                          index_nm1, basis_n)

    # -- cohomology side ---------------------------------------------------------
    def cochain_basis(self, alg, n):
        return list(self._cs(n + 1)._basis(n, "coh"))

    def coboundary_matrix(self, alg, n, basis_n, index_np1):
        cs = self._cs(n + 1)
        return self._pack(cs, cs.matrix(n, "coh"), cs._basis(n + 1, "coh"),
                          index_np1, basis_n)
