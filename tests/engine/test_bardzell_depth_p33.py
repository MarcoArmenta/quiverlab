"""Plan 33 -- self-injective Nakayama AT SCALE + the Bardzell depth showcase.

The symmetric Brauer stars kZ_n/J^L (Taillefer/Erdmann self-injective cyclic
Nakayama family) and the deg-300 depth regression that the Plan-33 src recursion
guard (D2) unlocks.  Split from the value batteries because the operative oracle
here is the ENGINE reaching depths the bar/CS routes cannot:

* Brauer-star SYMMETRY -- ``oracle_literature``: the Skowronski-Yamagata criterion
  kZ_n/J^L symmetric iff n | (L - 1) (``skowronski_yamagata``).
* Brauer-star Bardzell DEPTH -- ``oracle_crossengine``: the CS engine and the
  Bardzell engine agree degreewise in the shallow window, then Bardzell alone runs
  to degree 120 (self-injective => infinite complexity, the bar oracle is hopeless
  at dim 36 / 55).
* deg-300 regression -- UNMARKED (D2 infrastructure): kZ_20/J^11 (dim 220) HH_* to
  degree 300 completes via the engine-level recursion guard (no local
  ``sys.setrecursionlimit`` bump), pinning both reachability and the leading values.

``tests/engine/`` auto-assigns to the deep CI bucket (``tests/conftest.py``).  Every
value below was recomputed live (Plan 33) before pinning.
"""
import pytest

from quiverlab import GF
from quiverlab.families import NakayamaAlgebra
from quiverlab.resolutions_cs.homology import cs_homology_dims
from quiverlab.engine.hh_engine import hochschild_homology_dims
from quiverlab.engine.resolutions_bardzell import BardzellResolution, MonomialPresentation
from quiverlab.engine.coxeter2 import cyclic_nakayama

pytest.importorskip("quiverlab.groebner")

PRIME = 32003
F = GF(PRIME)


def _bardzell(n, L):
    """(engine algebra, Bardzell resolution) for the cyclic Nakayama kZ_n/J^L."""
    A, _ = cyclic_nakayama(n, L)
    return A, BardzellResolution(MonomialPresentation.cyclic_nakayama(n, L))


# =========================================================================== #
# Brauer stars kZ_4/J^9 (dim 36) and kZ_5/J^11 (dim 55)                          #
# =========================================================================== #
@pytest.mark.oracle_literature
@pytest.mark.parametrize("n,L,dim", [(4, 9, 36), (5, 11, 55)])
def test_brauer_star_symmetry(n, L, dim):
    """The symmetric Brauer stars kZ_4/J^9 (dim 36) and kZ_5/J^11 (dim 55): all four
    self-injectivity booleans are True -- SYMMETRIC (n | (L - 1): 4|8 and 5|10),
    weakly symmetric, Frobenius, self-injective (Skowronski-Yamagata symmetric-
    Nakayama criterion, ``skowronski_yamagata``).  Over GF(32003) (the trace-form
    certifier is loud on small fields)."""
    A = NakayamaAlgebra(n=n, l=L, cyclic=True, field=F)
    assert A.dim == dim
    assert (L - 1) % n == 0                                 # the symmetry criterion
    assert A.is_symmetric() is True
    assert A.is_weakly_symmetric() is True
    assert A.is_frobenius() is True
    assert A.is_selfinjective() is True


@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("n,L,shallow", [(4, 9, [6, 2, 2, 2, 2]), (5, 11, [7, 2, 2, 2, 2])])
def test_brauer_star_bardzell_depth(n, L, shallow):
    """Brauer-star HH_* cross-engine + depth: CS =(deg 4)= Bardzell on kZ_n/J^L
    (HH_* = [6,2,2,2,2] for kZ_4/J^9, [7,2,2,2,2] for kZ_5/J^11), then the Bardzell
    engine alone runs to DEGREE 120 (len 121) -- self-injective => infinite
    complexity, so HH_* never dies and the bar oracle (dim 36 / 55) is hopeless
    (``bardzell``; ``chouhy_solotar``)."""
    A_cs = NakayamaAlgebra(n=n, l=L, cyclic=True, field=F)
    cs4 = cs_homology_dims(A_cs, 4).dims
    Ae, res = _bardzell(n, L)
    bard4 = hochschild_homology_dims(Ae, 4, resolution=res)[PRIME]
    assert cs4 == bard4 == shallow                          # CS vs Bardzell, deg 4
    deep = hochschild_homology_dims(Ae, 120, resolution=res)[PRIME]
    assert len(deep) == 121                                 # Bardzell alone, deg 120
    assert deep[:5] == shallow                              # prefix consistent with CS


# =========================================================================== #
# deg-300 depth regression (D2 recursion guard) -- UNMARKED infrastructure       #
# =========================================================================== #
def test_bardzell_deg300_regression_kZ20_J11():
    """Plan-33 D2 regression: the Bardzell HH_* of kZ_20/J^11 (dim 220) reaches
    DEGREE 300 through the engine-level recursion guard -- WITHOUT any local
    ``sys.setrecursionlimit`` bump (the default Python limit alone stops it near
    deg 160).  This pins both the reachability (len == 301, ~15 s: the cost is the
    dim-220 basis enumeration, not the guard) and the leading values.  Infrastructure
    for the recursion guard, so UNMARKED (no external/second-engine oracle: the bar
    and CS routes are both infeasible at dim 220)."""
    Ae, res = _bardzell(20, 11)
    dims = hochschild_homology_dims(Ae, 300, resolution=res)[PRIME]
    assert len(dims) == 301                                 # reached degree 300 (guard)
    assert dims[:12] == [20, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
