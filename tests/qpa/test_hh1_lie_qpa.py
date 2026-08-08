"""GAP's core Lie library cross-checks HH^1 = Der/Inn (Plan 70 / R11).

QPA itself has NO HH^1-Lie surface, but the [qpa] libgap backend ships GAP's core
Lie machinery (LieAlgebraByStructureConstants, IsLieSolvable, LieDerivedSeries,
SemiSimpleType). Feeding it the computed HH^1 structure constants is a genuine
independent oracle -- an honest skip that FAILS if the Lie functions ever disappear."""
import pytest

pytestmark = pytest.mark.qpa

from quiverlab.qpa.session import gap_available


@pytest.mark.skipif(not gap_available(), reason="[qpa] libgap backend not installed")
def test_hh1_lie_gap_sl2_and_solvable():
    from quiverlab.fields import QQ
    from quiverlab.combinat.quiver import Quiver
    from quiverlab.families import truncated_polynomial
    from quiverlab.qpa.crosscheck import crosscheck_hh1_lie
    # sl2: GAP SemiSimpleType == "A1", IsLieSolvable == false
    r = crosscheck_hh1_lie(Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=QQ))
    r.assert_agree()
    assert r.gap_semisimple_type == "A1" and r.gap_solvable is False
    # k[x]/x^3 over QQ: solvable, derived series dims [2,1,0]
    r2 = crosscheck_hh1_lie(truncated_polynomial(3, field=QQ))
    r2.assert_agree()
    assert r2.gap_solvable is True
    assert r2.gap_derived_dims == [2, 1, 0]


@pytest.mark.skipif(not gap_available(), reason="[qpa] libgap backend not installed")
def test_hh1_lie_gap_trivial_extension_kronecker():
    """T(kK2) = k semidirect sl2: GAP confirms NOT solvable, and the semisimple
    Levi type sl2 -> A1 (SemiSimpleType is read only when our radical is 0, so this
    exercises solvable=False without the type read)."""
    from quiverlab.fields import QQ
    from quiverlab.combinat.quiver import Quiver
    from quiverlab.families import TrivialExtension
    from quiverlab.qpa.crosscheck import crosscheck_hh1_lie
    A = TrivialExtension(Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=QQ))
    r = crosscheck_hh1_lie(A)
    r.assert_agree()
    assert r.gap_solvable is False
