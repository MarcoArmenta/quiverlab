"""QPA cross-check battery over a small zoo (spec §8 ring 3). qpa-marked: skips
locally, mandatory under QUIVERLAB_REQUIRE_QPA=1 in CI."""
import pytest

from quiverlab import Quiver, GF
from quiverlab.fields import QQ          # QQ lives in quiverlab.fields (not re-exported top-level)
from quiverlab.qpa import session

# skip locally (no GAP); under QUIVERLAB_REQUIRE_QPA=1 the predicate is False, so the
# tests RUN and fail naturally if GAP is missing/broken (no silent green skip).
pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


@pytest.mark.parametrize("field", [GF(2), GF(3), QQ])
def test_hochschild_crosscheck_small_zoo(field):
    # commutative square kQ/(a*b - c*d) = kA_2 (x) kA_2 : HH = [1,0,0] (Kunneth)
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=field)
    A.crosscheck("hochschild", 2).assert_agree()


def test_module_self_ext_crosscheck():
    # self-Ext of the simple S1 in kA_2 over GF(2): Ext^0=Hom(S1,S1)=1, Ext^{>=1}=0
    # (no loop at vertex 1). ExtAlgebraGenerators(M, n)[1] gives the dim series.
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(2))
    A.crosscheck("module_ext", A.simple(1), 2).assert_agree()   # -> [1, 0, 0]
