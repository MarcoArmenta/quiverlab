"""B1 Increment 1: the Resolution protocol.

The bar complex, wrapped as BarResolution and passed via the new `resolution=` kwarg,
must reproduce the default engine EXACTLY (it is the same computation behind an
indirection). RED until hanlab.resolutions.BarResolution and the `resolution=` kwargs
exist; then GREEN with zero numeric change. This is the plug-in seam the Bardzell /
Chouhy-Solotar backends slot into (see notes/B1_chouhy_solotar_plan.md).
"""
import pytest

from quiverlab.engine.hh_engine import truncated_polynomial, hochschild_homology_dims
from quiverlab.engine.scan3 import quantum_ci, hochschild_cohomology_dims
from quiverlab.engine.coxeter import (
    nakayama_quantum_ci, induced_on_HH_homology, is_identity)
from quiverlab.engine.resolutions import BarResolution

# hanlab __init__ aliases, reproduced locally:
homology_dims = hochschild_homology_dims
cohomology_dims = hochschild_cohomology_dims
induced_homology = induced_on_HH_homology
PRIME = 32003

P = PRIME


@pytest.mark.parametrize("alg_fn,N", [
    (lambda: truncated_polynomial(3), 5),
    (lambda: quantum_ci(2), 5),
])
def test_bar_resolution_matches_default_homology(alg_fn, N):
    A = alg_fn()
    assert homology_dims(A, N, resolution=BarResolution())[P] == homology_dims(A, N)[P]


def test_bar_resolution_matches_default_cohomology():
    A = quantum_ci(2)
    assert cohomology_dims(A, 5, resolution=BarResolution())[P] == cohomology_dims(A, 5)[P]


def test_bar_resolution_matches_default_induced_action():
    Q = quantum_ci(2)
    S, Sinv = nakayama_quantum_ci(2, P)
    M_def, d_def = induced_homology(Q, 2, S, P)
    M_res, d_res = induced_homology(Q, 2, S, P, resolution=BarResolution())
    assert d_def == d_res
    assert is_identity(M_def, P) == is_identity(M_res, P)
