"""M2 script builders are pure string generation -- exercised without M2."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.m2 import scripts

pytestmark = pytest.mark.fast


def _dual_numbers():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(7))


def test_graded_dims_script_shape():
    s = scripts.graded_dims_script(_dual_numbers(), top=4)
    assert "ZZ/7" in s
    assert "x*x" in s
    assert scripts.SENTINEL in s
    assert "ncBasis" in s


def test_multi_vertex_refused():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))
    with pytest.raises(QuiverlabError, match="single-vertex"):
        scripts.graded_dims_script(A, top=3)


def test_caret_relations_translate():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x^3"], field=GF(5))
    s = scripts.graded_dims_script(A, top=3)
    assert "x^3" in s        # M2 shares the caret power syntax


def test_parse_sentinels_roundtrip():
    out = "junk\n<<QL>> 0 1\n<<QL>> 1 2\n<<QL>> 2 1\nnoise"
    assert scripts.parse_sentinels(out) == [1, 2, 1]


def test_parse_sentinels_rejects_noninteger():
    with pytest.raises(ValueError):
        scripts.parse_sentinels("<<QL>> 0 1.5")


def test_parse_sentinels_rejects_gap_in_degrees():
    with pytest.raises(ValueError):
        scripts.parse_sentinels("<<QL>> 0 1\n<<QL>> 2 1")


def test_commutative_ext_script_shape():
    s = scripts.commutative_ext_script(7, ["x", "y"], ["x^2", "y^2"], top=6)
    assert "ZZ/7[x,y]" in s.replace(" ", "")
    assert "freeResolution" in s and "LengthLimit" in s
    assert scripts.SENTINEL in s
