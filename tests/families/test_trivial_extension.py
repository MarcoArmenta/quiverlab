"""TrivialExtension(A) = A |x D(A), symmetric. Fixture N10."""
from quiverlab.combinat import Quiver
from quiverlab.families import TrivialExtension
from quiverlab.fields import GF


def _kA2(field):
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=field)


def test_trivial_extension_kA2_dim6():
    A = TrivialExtension(_kA2(GF(32003)))
    assert A.dim == 6                                       # 2 * dim(kA_2) = 2*3 (Fixture N10)
    A._validate()                                          # regression guard: block multiplication stays associative + two-sided unit (raises on failure)


def test_trivial_extension_is_symmetric_hh_duality():
    A = TrivialExtension(_kA2(GF(32003)))                  # symmetric: HH^n dim = HH_n dim
    # T(A) carries no quiver, so HH runs on the exponential normalized-bar complex:
    # dim C^n = m*(m-1)^n with m = 6, so the degree-3 coboundary d^3 is 3750 x 750
    # = 2.8M cells, under the 4M default guard -- no max_cells override needed.
    # Degree 4 was verified once at implementation time (HH^4 == HH_4, values equal) but is
    # excluded here for cost (~61s / ~1.3 GB at max_cells=71M); the degree-3 check pins the
    # same symmetric-duality property HH^n == HH_n.
    co = A.hochschild_cohomology(3).dims
    ho = A.hochschild_homology(3).dims
    assert co == ho
