"""TrivialExtension(A) = A |x D(A), symmetric. Fixture N10."""
from quiverlab.combinat import Quiver
from quiverlab.families import TrivialExtension
from quiverlab.fields import GF


def _kA2(field):
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=field)


def test_trivial_extension_kA2_dim6():
    A = TrivialExtension(_kA2(GF(32003)))
    assert A.dim == 6                                       # 2 * dim(kA_2) = 2*3 (Fixture N10)


def test_trivial_extension_is_symmetric_hh_duality():
    A = TrivialExtension(_kA2(GF(32003)))                  # symmetric: HH^n dim = HH_n dim
    # T(A) carries no quiver, so HH runs on the exponential normalized-bar complex:
    # dim C^n = m*(m-1)^n with m = 6, so the degree-4 coboundary d^4 is 18750 x 3750
    # = 70.3M cells, above the 4M default guard. Raise max_cells to admit degree 4.
    # (FLAG p06-t9: transcribed test used default max_cells and hit DepthLimitError.)
    co = A.hochschild_cohomology(4, max_cells=71_000_000).dims
    ho = A.hochschild_homology(4, max_cells=71_000_000).dims
    assert co == ho
