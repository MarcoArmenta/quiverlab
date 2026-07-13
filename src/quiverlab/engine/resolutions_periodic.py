# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
#
# Scope note (quiverlab port): the two lazy engine imports below are package-relative.
# QuantumCIResolution wraps resolutions_cs, which is EXCLUDED from Plan 02 (not ported),
# so that backend is dormant here -- importing it raises ImportError only if the class is
# actually instantiated and used. CyclicNakayamaResolution wraps resolutions_bardzell,
# which lands in Task 10, so it self-heals then. The module itself imports cleanly today
# (only `resolutions.Resolution`, already ported) and is import-closed.
"""Family-specific resolution backends for the eventually-periodic algebras.

Two homology-only `Resolution` backends that let a caller request the small/fast
resolution for a known algebra family by parameter, without hand-building a reduction
system or a monomial presentation:

  * `QuantumCIResolution(c)`      -- quantum complete intersection
                                     A = k<x,y>/(x^2, y^2, yx - c*xy)   (H.quantum_ci(c)),
                                     the NON-monomial headline algebra. Eventually periodic:
                                     P_n has n+1 generators, vs the bar complex's 4*3^n.
  * `CyclicNakayamaResolution(n, ell)` -- self-injective Nakayama A = kZ_n/rad^ell
                                     (H.cyclic_nakayama(n, ell)).

HONEST SCOPE.  These are *thin, verified wrappers*, not independent re-derivations:
`QuantumCIResolution` delegates to the Chouhy-Solotar closed form
(`resolutions_cs.ChouhySolotarResolution`) and `CyclicNakayamaResolution` to the Bardzell
minimal resolution (`resolutions_bardzell.BardzellResolution`). Both underlying engines are
cross-checked exactly against the normalized-bar-complex oracle in their own test files; the
wrappers add a by-family entry point and are themselves cross-checked against the bar oracle
in tests/test_periodic_resolution.py. They provide the SAME numbers as the engines they wrap
(no independent signal) -- use Bardzell-vs-CS-vs-bar for independent confirmation.

CONTRACT (see resolutions.py docstring):
  term_basis(alg, n) -> ordered list of hashable generators (len = dim P_n).
  differential_matrix(alg, n, basis_n, index_nm1) -> int64 matrix of shape
      (len(index_nm1), len(basis_n)) == (dim P_{n-1}, dim P_n), NEVER pre-reduced mod p.
"""

from quiverlab.engine.resolutions import Resolution


class QuantumCIResolution(Resolution):
    """Eventually-periodic backend for the quantum complete intersection
    A = k<x,y>/(x^2, y^2, yx - c*xy). Verified wrapper over the Chouhy-Solotar closed form.

    Usage:
        A = H.quantum_ci(2)
        H.homology_dims(A, N, resolution=QuantumCIResolution(2))
    """

    def __init__(self, c):
        self.c = c
        self._inner = None

    def _ensure(self, alg):
        if self._inner is None:
            from quiverlab.engine.resolutions_cs import (  # excluded from Plan 02 -> dormant
                ChouhySolotarResolution, qci_reduction_system)
            self._inner = ChouhySolotarResolution(qci_reduction_system(self.c), alg)
        return self._inner

    def term_basis(self, alg, n):
        return self._ensure(alg).term_basis(alg, n)

    def differential_matrix(self, alg, n, basis_n, index_nm1):
        return self._ensure(alg).differential_matrix(alg, n, basis_n, index_nm1)


class CyclicNakayamaResolution(Resolution):
    """Backend for the self-injective Nakayama algebra A = kZ_n/rad^ell. Verified wrapper
    over the Bardzell minimal monomial resolution.

    Usage:
        A, _ = H.cyclic_nakayama(3, 2)
        H.homology_dims(A, N, resolution=CyclicNakayamaResolution(3, 2))
    """

    def __init__(self, n, ell):
        self.n = n
        self.ell = ell
        self._inner = None

    def _ensure(self):
        if self._inner is None:
            from quiverlab.engine.resolutions_bardzell import (  # self-heals at Task 10
                BardzellResolution, MonomialPresentation)
            self._inner = BardzellResolution(
                MonomialPresentation.cyclic_nakayama(self.n, self.ell))
        return self._inner

    def term_basis(self, alg, n):
        return self._ensure().term_basis(alg, n)

    def differential_matrix(self, alg, n, basis_n, index_nm1):
        return self._ensure().differential_matrix(alg, n, basis_n, index_nm1)
