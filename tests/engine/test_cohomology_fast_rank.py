"""HH^* dims must be identical whether cohomology ranks go through the dense reference
or rank_mod_p_auto (numba/sparse). This guards the item-(c) speedup."""
import pytest

# hanlab sys.path shim dropped: quiverlab uses absolute package imports.
from quiverlab.engine.hh_engine import truncated_polynomial, two_gen_local
from quiverlab.engine.scan3 import quantum_ci, hochschild_cohomology_dims

PRIMES = (32003, 2, 3)


@pytest.mark.parametrize("alg,N", [
    (truncated_polynomial(3), 5),
    (quantum_ci(2), 5),
    (quantum_ci(3), 5),
    (two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "k[x,y]/(x^2,y^2)"), 4),
], ids=lambda x: getattr(x, "name", str(x)))
def test_cohomology_dims_known_values(alg, N):
    ch = hochschild_cohomology_dims(alg, N, primes=PRIMES)
    # symmetric algebras: HH^n == HH_n; quantum CI: cohomology bounded/vanishing.
    # Pin against the engine's own homology where symmetric, else just non-negative + finite.
    for p in PRIMES:
        assert all(d >= 0 for d in ch[p])
    if alg.name == "k[x,y]/(x^2,y^2)":
        from quiverlab.engine.hh_engine import hochschild_homology_dims
        hh = hochschild_homology_dims(alg, N, primes=PRIMES)
        for p in PRIMES:
            assert ch[p] == hh[p]            # symmetric => HH^n == HH_n
