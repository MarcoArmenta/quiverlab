# SPDX-License-Identifier: MIT
"""The public surface of the cluster slice (Plan 79 / R31, Task 6)."""
import pytest

from quiverlab.fields import QQ

selfcert = pytest.mark.oracle_selfcert


@selfcert
def test_cluster_category_is_importable_from_the_top_level():
    from quiverlab import ClusterCategory
    assert ClusterCategory("A2").n == 2


@selfcert
def test_the_algebra_delegate_serves_the_same_object():
    from quiverlab import ClusterCategory, PathAlgebra
    A = PathAlgebra("A3", field=QQ)
    C = A.cluster_category()
    assert isinstance(C, ClusterCategory)
    assert C.num_cluster_tilting()["count"] == 14
    assert C.num_indecomposables() == 9


@selfcert
def test_the_delegate_refuses_a_non_hereditary_algebra():
    from quiverlab import Quiver
    from quiverlab.errors import QuiverlabError
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)
    with pytest.raises(QuiverlabError, match="hereditary"):
        A.cluster_category()
