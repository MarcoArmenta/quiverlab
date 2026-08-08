"""The annulus has one hole: pi1^ab = Z. VERIFIED LIVE during authoring:
jacobian_of(annulus_triangulation(2,2)) is 4 vertices, 4 arrows, EMPTY relations
(hereditary) -- so free_rank = 1 is the first Betti number of a 4-cycle (a DERIVED
fact, not a published value): oracle_selfcert."""
import pytest
from quiverlab import annulus_triangulation, jacobian_of
from quiverlab.invariants.coverings import fundamental_group

pytestmark = pytest.mark.oracle_selfcert


def test_annulus_C22_pi1ab_is_Z():
    A = jacobian_of(annulus_triangulation(2, 2))   # 4 vertices, 4 arrows, no relations
    assert A.relations == []                        # hereditary (re-verify the shape)
    g = fundamental_group(A)
    assert g.free_rank == 1 and g.invariant_factors == ()   # Betti(4-cycle) = 1 => Z
