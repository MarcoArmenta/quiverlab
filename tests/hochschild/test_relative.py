"""Relative HH over the vertex subalgebra E=kQ0. Separable E => relative ==
absolute degreewise (free strong oracle); byte-identical single-vertex; general B
refused loudly."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.errors import QuiverlabError

pytestmark = pytest.mark.oracle_crossengine


def test_relative_equals_absolute_multivertex():
    A = ql.Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=ql.CC)
    for k in range(6):
        assert A.hochschild_homology(k, relative_to="vertices").dims == \
            A.hochschild_homology(k).dims                         # separability


def test_relative_cohomology_equals_absolute_multivertex():
    # extra (beyond the plan's list): validate the E-relative COBOUNDARY on a
    # multi-vertex algebra, not only single-vertex.
    A = ql.Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=ql.CC)
    for k in range(6):
        assert A.hochschild_cohomology(k, relative_to="vertices").dims == \
            A.hochschild_cohomology(k).dims


@pytest.mark.oracle_selfcert
def test_relative_byte_identical_single_vertex():
    A = ql.truncated_polynomial(3, field=ql.GF(7))                # E = k
    assert A.hochschild_cohomology(5, relative_to="vertices").dims == \
        A.hochschild_cohomology(5).dims
    assert A.hochschild_homology(5, relative_to="vertices").dims == \
        A.hochschild_homology(5).dims


def test_relative_with_coefficients():
    A = ql.Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=ql.CC)
    M = Bimodule.dual(A)
    for k in range(5):
        assert A.hochschild_homology(k, relative_to="vertices", coefficients=M).dims == \
            A.hochschild_homology(k, coefficients=M).dims
        assert A.hochschild_cohomology(k, relative_to="vertices", coefficients=M).dims == \
            A.hochschild_cohomology(k, coefficients=M).dims


def test_general_B_refused():
    A = ql.truncated_polynomial(3, field=ql.CC)
    with pytest.raises(QuiverlabError):
        A.hochschild_homology(3, relative_to="some-subalgebra")
    with pytest.raises(QuiverlabError):
        A.hochschild_cohomology(3, relative_to="some-subalgebra")
