"""Minimal A^e HH(A,M) over GF(p) == bar HH(A,M), every named coefficient.
Plan-16 covariance generalized: hom b·w·a on e_w M e_v, coh a·w·b on the SWAPPED
e_v M e_w. Multi-vertex included (the corner path).

Bar is exponential (cochain basis dim_M*(m-1)^n), so the minimal==bar oracle runs
to the deepest degree bar can CERTIFY (top capped per algebra dim); the minimal
engine itself is deep. (DEVIATION FROM PLAN top=5/6 on multi-vertex: bar's blowup.)
"""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.engine.resolutions_minimal import minimal_homology_dims, minimal_cohomology_dims

pytestmark = pytest.mark.oracle_crossengine


@pytest.mark.parametrize("mk", ["regular", "dual", "twisted_by_nakayama"])
def test_minimal_matches_bar_single_vertex(mk):
    A = ql.truncated_polynomial(3, field=ql.GF(32003))     # dim 3: bar reaches depth 6
    M = getattr(Bimodule, mk)(A)
    N = 6
    assert minimal_homology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_homology(N, engine="bar", coefficients=M).dims
    assert minimal_cohomology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_cohomology(N, engine="bar", coefficients=M).dims


def test_minimal_matches_bar_multivertex():
    # kA3 with a*b=0: exercises the corner path + SWAPPED-tag coh block. dim 5, so
    # bar's cochain basis 5*4^n caps the oracle at top=3 (d^4 is 6.5M > max_cells).
    A = ql.Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=ql.GF(32003))
    M = Bimodule.dual(A)
    N = 3
    assert minimal_cohomology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_cohomology(N, engine="bar", coefficients=M).dims
    assert minimal_homology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_homology(N, engine="bar", coefficients=M).dims


def test_minimal_matches_bar_nontrivial_nu_multivertex():
    # C2: multi-vertex NON-symmetric self-injective over GF(p) -- exercises the
    # twist path (nu != id) AND the corner/SWAPPED-tag block simultaneously, and
    # is the DD1b basis discriminator on the minimal engine. dim 9 => bar caps top=2.
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=ql.GF(32003))
    assert A.is_symmetric() is False
    M = Bimodule.twisted_by_nakayama(A)
    N = 2
    assert minimal_cohomology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_cohomology(N, engine="bar", coefficients=M).dims
    assert minimal_homology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        A.hochschild_homology(N, engine="bar", coefficients=M).dims


def test_minimal_matches_cs_deep_nontrivial_nu():
    # Deeper cross-check where bar cannot reach: minimal == CS to depth 5 on the
    # dim-9 non-symmetric self-injective algebra (both non-bar engines).
    from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=ql.GF(32003))
    M = Bimodule.dual(A)
    N = 5
    assert minimal_cohomology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        cs_cohomology_dims(A, N, coefficients=M).dims
    assert minimal_homology_dims(A, N, primes=(32003,), coefficients=M)[32003] == \
        cs_homology_dims(A, N, coefficients=M).dims


def test_bar_transport_deep_multivertex():
    # Review micro-fix 1: the bar UNIT-ADAPTATION TRANSPORT of a coefficient built
    # in A's own basis, exercised at DEPTH (top=4) on a NON-unit-adapted multi-vertex
    # algebra -- bar(transport) ≡ CS ≡ minimal, all three engines. kA3 (a*b=0, dim 5,
    # unit = e1+e2+e3 so A is NOT unit-adapted); bar's cochain basis stays sparse so
    # the exact rank reaches degree 4 (~5s) with a raised max_cells. This closes the
    # gap where the shipped multi-vertex bar-transport coverage stopped at degree 2/3.
    from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
    A = ql.Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=ql.GF(32003))
    assert A.is_unit_adapted is False           # the transport is genuinely exercised
    M = Bimodule.dual(A)
    N, MC = 4, 20_000_000
    bar_c = A.hochschild_cohomology(N, engine="bar", coefficients=M,
                                    max_cells=MC, verbose=False).dims
    bar_h = A.hochschild_homology(N, engine="bar", coefficients=M,
                                  max_cells=MC, verbose=False).dims
    assert bar_c == cs_cohomology_dims(A, N, coefficients=M).dims == \
        minimal_cohomology_dims(A, N, primes=(32003,), coefficients=M)[32003]
    assert bar_h == cs_homology_dims(A, N, coefficients=M).dims == \
        minimal_homology_dims(A, N, primes=(32003,), coefficients=M)[32003]


def test_minimal_refuses_non_gfp_coefficient():
    from quiverlab.errors import QuiverlabError
    A = ql.truncated_polynomial(3, field=ql.CC)
    M = Bimodule.dual(A)
    with pytest.raises(QuiverlabError):
        minimal_cohomology_dims(A, 3, primes=(0,), coefficients=M)
