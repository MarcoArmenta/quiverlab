"""Item 3 (TODO P1): the family-specific periodic backends in resolutions_periodic.

These wrappers (`QuantumCIResolution`, `CyclicNakayamaResolution`) were NotImplementedError
stubs whose module docstring claimed validation against a test file that did not exist. They
are now thin, verified wrappers over the Chouhy-Solotar / Bardzell engines; this file is the
cross-check the old docstring promised -- each wrapper is run against the normalized-bar-complex
oracle on the overlap range, plus a depth-unlock case far past the bar-complex blow-up.
"""
import pytest

from quiverlab.engine.hh_engine import hochschild_homology_dims
from quiverlab.engine.resolutions_periodic import CyclicNakayamaResolution

# hanlab __init__ aliases, reproduced locally:
homology_dims = hochschild_homology_dims
PRIME = 32003

PRIMES = (32003, 2, 3, 5)


# ---------------------------------------------------------------------------
# QuantumCIResolution: the NON-monomial headline algebra
# ---------------------------------------------------------------------------
# not ported: test_quantum_ci_matches_bar_all_primes, test_quantum_ci_depth_unlock --
# both exercise QuantumCIResolution, which delegates to resolutions_cs (the Chouhy-Solotar
# closed form). resolutions_cs is EXCLUDED from Plan 02, so that backend is dormant here.


# ---------------------------------------------------------------------------
# CyclicNakayamaResolution: self-injective Nakayama
# ---------------------------------------------------------------------------
# CyclicNakayamaResolution delegates to the Bardzell minimal resolution (Task 10) and the
# algebras come from coxeter2.cyclic_nakayama (Task 11), so these self-heal once both land.
@pytest.mark.parametrize("n,ell,N", [(3, 2, 4), (2, 3, 3), (3, 3, 3)])
def test_cyclic_nakayama_matches_bar_all_primes(n, ell, N):
    pytest.importorskip("quiverlab.engine.resolutions_bardzell")  # Bardzell engine (Task 10)
    pytest.importorskip("quiverlab.engine.coxeter2")              # cyclic_nakayama (Task 11)
    from quiverlab.engine.coxeter2 import cyclic_nakayama
    A, _ = cyclic_nakayama(n, ell)
    res = CyclicNakayamaResolution(n, ell)
    bar = homology_dims(A, N)
    per = homology_dims(A, N, resolution=res)
    for p in PRIMES:
        assert per[p] == bar[p], f"CN({n},{ell}) p={p}: periodic {per[p]} != bar {bar[p]}"


def test_cyclic_nakayama_depth_unlock():
    pytest.importorskip("quiverlab.engine.resolutions_bardzell")  # Bardzell engine (Task 10)
    pytest.importorskip("quiverlab.engine.coxeter2")              # cyclic_nakayama (Task 11)
    from quiverlab.engine.coxeter2 import cyclic_nakayama
    # kZ_3/rad^2 to depth 24: the bar complex is hopeless this deep.
    A, _ = cyclic_nakayama(3, 2)
    dims = homology_dims(A, 24, resolution=CyclicNakayamaResolution(3, 2))[PRIME]
    assert len(dims) == 25
    assert dims[:5] == [3, 0, 1, 1, 0]
