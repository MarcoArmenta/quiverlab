"""TrivialExtension(A) = A |x D(A): a presented symmetric algebra (Plan 31). Fixture N10."""
import pytest

from quiverlab.combinat import Quiver
from quiverlab.families import TrivialExtension
from quiverlab.fields import GF


def _kA2(field):
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=field)


@pytest.mark.oracle_selfcert
def test_trivial_extension_kA2_dim6():
    A = TrivialExtension(_kA2(GF(32003)))
    assert A.dim == 6                                       # 2 * dim(kA_2) = 2*3 (Fixture N10)
    A._validate()                                          # regression guard: block multiplication stays associative + two-sided unit (raises on failure)


@pytest.mark.oracle_literature
def test_trivial_extension_is_symmetric_hh_duality():
    A = TrivialExtension(_kA2(GF(32003)))                  # symmetric: HH^n dim = HH_n dim
    # Plan 31: T(A) is now a genuine kQ_T/I_T-presented Algebra; over GF(32003) the
    # auto engine is the fast bar-rank engine (int64 mod-p rank on the normalized-bar
    # complex). The values are iso-invariant (unchanged from the retained
    # structure-constant build). Degree 3 pins the symmetric-duality HH^n == HH_n;
    # degree 4 is excluded for cost.
    assert A.is_symmetric() is True
    co = A.hochschild_cohomology(3).dims
    ho = A.hochschild_homology(3).dims
    assert co == ho
