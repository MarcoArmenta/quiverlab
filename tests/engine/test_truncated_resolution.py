"""B1 Increment 2: TruncatedPolynomialResolution -- a small backend proving the depth
unlock. P_n = A (rank 1) for all n, vs the bar complex's a*(a-1)^n. Validated against
the bar engine on the overlap range, then run far past the bar-complex depth cap.
"""
import pytest

from quiverlab.engine.hh_engine import truncated_polynomial, hochschild_homology_dims
from quiverlab.engine.resolutions import TruncatedPolynomialResolution

# hanlab __init__ alias, reproduced locally:
homology_dims = hochschild_homology_dims
PRIME = 32003

P = PRIME


@pytest.mark.parametrize("a", [2, 3, 4])
def test_matches_bar_engine_char0(a):
    A = truncated_polynomial(a)
    N = 5
    bar = homology_dims(A, N)[P]
    cs = homology_dims(A, N, resolution=TruncatedPolynomialResolution())[P]
    assert cs == bar == [a] + [a - 1] * N


def test_depth_unlock_far_past_bar_cap():
    # k[x]/(x^4): the bar complex caps near N=6 (dim C_6 = 4*3^6 = 2916); the periodic
    # resolution reaches N=40 instantly with rank-1-per-degree complexes.
    A = truncated_polynomial(4)
    dims = homology_dims(A, 40, resolution=TruncatedPolynomialResolution())[P]
    assert dims == [4] + [3] * 40


def test_char_p_dividing_a_collapses_all_maps():
    # k[x]/(x^3) at p = 3: a = 3 == 0 mod 3, so every induced map vanishes; HH_n = 3.
    A = truncated_polynomial(3)
    dims = homology_dims(A, 20, resolution=TruncatedPolynomialResolution())[3]
    assert dims == [3] * 21
