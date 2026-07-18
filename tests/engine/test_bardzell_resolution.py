"""Bardzell minimal bimodule resolution for monomial algebras: cross-checked exactly
against the normalized-bar-complex oracle (homology_dims) on the overlap range, for
each monomial test algebra and each prime in {32003, 2, 3, 5}, plus a depth-unlock
assertion running the backend to N >= 20 (far past the bar complex's blow-up).

Mirrors tests/engine/test_truncated_resolution.py.
"""
import pytest

from quiverlab.engine.hh_engine import truncated_polynomial, hochschild_homology_dims
from quiverlab.engine.scan2 import local_3gen_radsq
from quiverlab.engine.resolutions_bardzell import BardzellResolution, MonomialPresentation

# hanlab __init__ aliases, reproduced locally:
homology_dims = hochschild_homology_dims
PRIME = 32003

PRIMES = (32003, 2, 3, 5)


# ----------------------------------------------------------------------
# 1. k[x]/(x^a) -- single loop, single relation x^a
# ----------------------------------------------------------------------
@pytest.mark.parametrize("a", [2, 3, 4])
def test_truncated_polynomial_matches_oracle(a):
    A = truncated_polynomial(a)
    pres = MonomialPresentation.truncated_polynomial(a)
    res = BardzellResolution(pres)
    N = 6
    for p in PRIMES:
        bar = homology_dims(A, N)[p]
        bard = homology_dims(A, N, resolution=res)[p]
        assert bard == bar, f"k[x]/(x^{a}) p={p}: {bard} != {bar}"


# ----------------------------------------------------------------------
# 2. cyclic Nakayama kZ_n / rad^ell -- a single cycle; ell > 2 exercises
#    genuinely OVERLAPPING relations (the general overlap-chain differential).
#
# N is capped per dimension because the cross-check routes through the bar
# oracle, whose term blows up as dim C_n = m*(m-1)^n (m = n*ell). For dim-8
# (kZ_4/rad^2) and dim-10 (kZ_5/rad^2) the oracle's b_5 would need a ~289 GiB
# allocation at N=4, so those compare at N=2 -- still a 4-prime cross-check of
# HH_0..HH_2. Deep Bardzell behaviour is covered by test_depth_unlock_* below,
# and the ell>2 OVERLAP differential is cross-checked at full depth N=3 by the
# small-dimension (2,3,3) and (3,3,3) cases.
#
# quiverlab port: cyclic_nakayama (the algebra builder) lives in coxeter2 (Task 11),
# so this cross-check is coxeter2-gated -- it self-heals to PASS once T11 lands. The
# MonomialPresentation.cyclic_nakayama builder is local to resolutions_bardzell.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("n,ell,N", [
    (3, 2, 4),
    (4, 2, 2),   # dim 8: bar oracle infeasible past N~3 (m*(m-1)^n blow-up)
    (5, 2, 2),   # dim 10: bar oracle needs ~289 GiB at N=4; cap the comparison
    (2, 3, 3),   # overlapping relations (ell > 2)
    (3, 3, 3),   # overlapping relations (ell > 2)
])
def test_cyclic_nakayama_matches_oracle(n, ell, N):
    pytest.importorskip("quiverlab.engine.coxeter2")  # cyclic_nakayama (Task 11)
    from quiverlab.engine.coxeter2 import cyclic_nakayama
    A, _ = cyclic_nakayama(n, ell)
    pres = MonomialPresentation.cyclic_nakayama(n, ell)
    res = BardzellResolution(pres)
    for p in PRIMES:
        bar = homology_dims(A, N)[p]
        bard = homology_dims(A, N, resolution=res)[p]
        assert bard == bar, f"CN({n},{ell}) p={p}: {bard} != {bar}"


# ----------------------------------------------------------------------
# 3. local_3gen_radsq -- one vertex, 3 loops, all length-2 relations.
#    The richest monomial test: oracle HH starts [4, 6, 14, 32, ...].
# ----------------------------------------------------------------------
def test_local_3gen_radsq_matches_oracle():
    A = local_3gen_radsq()
    pres = MonomialPresentation.local_radsq(3)
    res = BardzellResolution(pres)
    N = 4
    for p in PRIMES:
        bar = homology_dims(A, N)[p]
        bard = homology_dims(A, N, resolution=res)[p]
        assert bard == bar, f"local_3gen_radsq p={p}: {bard} != {bar}"
    # sanity-check the known oracle values (char 0 proxy and p = 2)
    bard = homology_dims(A, 4, resolution=res)
    assert bard[32003][:4] == [4, 6, 14, 32]
    assert bard[2][:4] == [4, 9, 17, 35]


# ----------------------------------------------------------------------
# 4. DEPTH UNLOCK -- the backend reaches depths the bar complex cannot.
# ----------------------------------------------------------------------
def test_depth_unlock_truncated_polynomial():
    # k[x]/(x^4): bar complex caps near N=6 (dim C_6 = 4*3^6); Bardzell is rank-a
    # per degree and reaches N=40 instantly.
    A = truncated_polynomial(4)
    pres = MonomialPresentation.truncated_polynomial(4)
    dims = homology_dims(A, 40, resolution=BardzellResolution(pres))[PRIME]
    assert dims == [4] + [3] * 40


def test_depth_unlock_cyclic_nakayama():
    pytest.importorskip("quiverlab.engine.coxeter2")  # cyclic_nakayama (Task 11)
    from quiverlab.engine.coxeter2 import cyclic_nakayama
    # kZ_3/rad^2 to depth 24: bar complex is hopeless this deep.
    A, _ = cyclic_nakayama(3, 2)
    pres = MonomialPresentation.cyclic_nakayama(3, 2)
    dims = homology_dims(A, 24, resolution=BardzellResolution(pres))[PRIME]
    assert len(dims) == 25
    assert dims[:5] == [3, 0, 1, 1, 0]


def test_depth_unlock_local_3gen_radsq():
    # local_3gen_radsq to depth 5: the bar complex has dim C_n = 4*9^n (already
    # ~2.3e5 at n=5, infeasible to take ranks of), while Bardzell stays at 4*3^n and
    # reproduces the full HH sequence.
    A = local_3gen_radsq()
    pres = MonomialPresentation.local_radsq(3)
    dims = homology_dims(A, 5, resolution=BardzellResolution(pres))[PRIME]
    assert dims == [4, 6, 14, 32, 72, 170]
    assert all(d > 0 for d in dims)            # complexity >= 1 (Han signature)
