"""PreprojectiveAlgebra Pi(Q). Fixture N8."""
import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.families import PreprojectiveAlgebra
from quiverlab.fields import CC


def test_preprojective_A2_dim4():
    A = PreprojectiveAlgebra("A2", field=CC)
    assert A.dim == 4                                   # monomial mesh relations
    assert A.hochschild_cohomology(0).dims == [1]


def test_preprojective_A3_dim10():
    A = PreprojectiveAlgebra("A3", field=CC)            # tetrahedral number T_3
    assert A.dim == 10                                  # Groebner-certified oracle


def test_preprojective_affine_refused_loudly():
    # Pi(~A2) is infinite-dimensional; refuse up front (do NOT launch a doomed
    # Groebner completion). The message must name the type and say it is infinite.
    with pytest.raises(QuiverlabError) as exc:
        PreprojectiveAlgebra("~A2", field=CC)
    msg = str(exc.value)
    assert "~A2" in msg and "infinite-dimensional" in msg


def test_preprojective_finite_dynkin_An_still_build():
    # The guard must not disturb the finite Dynkin A_n cases.
    assert PreprojectiveAlgebra("A2", field=CC).dim == 4
    assert PreprojectiveAlgebra("A3", field=CC).dim == 10
