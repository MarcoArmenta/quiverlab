"""engine='ghms': fast HH via the GHMS resolution == the minimal syzygy engine == bar/CS,
degreewise -- the THIRD independent HH oracle class (Plan 75 / R10). Over any Domain.
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.resolutions_minimal import minimal_homology_dims
from quiverlab.errors import QuiverlabError
from quiverlab.families import IncidenceAlgebra
from quiverlab.families.basic import truncated_polynomial
from quiverlab.families.exterior import ExteriorAlgebra
from quiverlab.families.preprojective import PreprojectiveAlgebra
from quiverlab.fields import GF, QQ
from quiverlab.hochschild.koszul_ghms import ghms_cohomology_dims, ghms_homology_dims

xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature

P = 32003


def _kz3(dom):
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}).algebra(
        relations=["a*b", "b*c", "c*a"], field=dom)


def _diamond(dom):
    return IncidenceAlgebra([("b", "x"), ("b", "y"), ("x", "t"), ("y", "t")], field=dom)


# --------------------------------------------- GHMS == the minimal syzygy engine
@xeng
@pytest.mark.parametrize("name,build,top,expect", [
    ("exterior2", lambda d: ExteriorAlgebra(2, field=d), 6, [3, 4, 6, 8, 10, 12, 14]),
    ("exterior3", lambda d: ExteriorAlgebra(3, field=d), 4, [5, 12, 24, 40, 60]),
    ("kz3_radsq", _kz3, 5, [3, 0, 1, 1, 0, 0]),
    ("diamond", _diamond, 4, [4, 0, 0, 0, 0]),
])
def test_ghms_homology_equals_minimal(name, build, top, expect):
    A = build(GF(P))
    got = ghms_homology_dims(A, top)
    assert got == expect                                        # LIVE-VERIFIED
    assert got == list(minimal_homology_dims(to_engine(A), top, primes=(P,))[P])


@xeng
@pytest.mark.parametrize("build,top", [
    (lambda d: ExteriorAlgebra(2, field=d), 3),
    (_kz3, 3),
    (_diamond, 3),
])
def test_ghms_cohomology_equals_cs_over_QQ(build, top):
    """Domain-general: the syzygy engine is GF(p)-only, GHMS runs over QQ directly."""
    A = build(QQ)
    assert ghms_cohomology_dims(A, top) == list(
        A.hochschild_cohomology(top, engine="cs", verbose=False).dims)


@xeng
def test_engine_kwarg_routes_and_agrees():
    A = ExteriorAlgebra(2, field=QQ)
    assert (list(A.hochschild_homology(4, engine="ghms", verbose=False).dims)
            == list(A.hochschild_homology(4, engine="cs", verbose=False).dims))
    assert (list(A.hochschild_cohomology(4, engine="ghms", verbose=False).dims)
            == list(A.hochschild_cohomology(4, engine="cs", verbose=False).dims))


@lit
def test_diamond_is_the_bridge_between_both_halves_of_the_plan():
    """The commutative-square incidence algebra is Koszul AND an incidence algebra, so
    BOTH engines of this plan apply -- and they must not be confused:

      * Hochschild COHOMOLOGY = the order-complex cohomology  = [1, 0, 0]  (R9 theorem)
      * Hochschild HOMOLOGY   = [4, 0, 0, 0, 0]  (HH_0 = dim A/[A,A] = 4, NOT 1)

    Conflating the two is the error the plan's critic caught, so both are pinned here and
    the GHMS route is checked against the order-complex route on the cohomology side.
    """
    A = _diamond(QQ)
    assert ghms_cohomology_dims(A, 2) == A.incidence_cohomology(2).dims == [1, 0, 0]
    assert ghms_homology_dims(A, 4) == [4, 0, 0, 0, 0]


# ------------------------------------------------------------------- refusals
@selfcert
def test_engine_ghms_refuses_non_quadratic():
    with pytest.raises(QuiverlabError, match="[Kk]oszul|quadratic|QUADRATIC"):
        truncated_polynomial(3, field=QQ).hochschild_homology(3, engine="ghms",
                                                              verbose=False)


@lit
def test_engine_ghms_refuses_not_koszul_but_minimal_still_computes():
    """Pi(A_3) is NOT Koszul, so GHMS refuses NAMING the Ext obstruction -- while the
    general syzygy engine still computes its HH. The refusal is a scope boundary, not a
    capability gap."""
    with pytest.raises(QuiverlabError, match="NOT Koszul|obstruction"):
        PreprojectiveAlgebra("A3", field=QQ).hochschild_homology(3, engine="ghms",
                                                                 verbose=False)
    assert list(minimal_homology_dims(
        to_engine(PreprojectiveAlgebra("A3", field=GF(P))), 5, primes=(P,))[P]) \
        == [3, 1, 1, 1, 1, 1]


@selfcert
def test_unknown_engine_still_refuses():
    with pytest.raises(QuiverlabError, match="unknown engine"):
        ExteriorAlgebra(2, field=QQ).hochschild_homology(2, engine="nope", verbose=False)


@selfcert
def test_auto_is_unchanged_by_this_plan():
    """`auto` is deliberately NOT routed to GHMS, so every pre-existing result stays
    byte-identical and no golden can drift."""
    A = ExteriorAlgebra(2, field=QQ)
    auto = A.hochschild_homology(3, engine="auto", verbose=False)
    assert "GHMS" not in auto.engine
