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


def test_coefficient_and_sign_render_verbatim():
    # Devil's-advocate finding (2026-08-05): the graded-dims oracle is
    # coefficient-blind on quadratic CIs (Hilbert [1,2,1] for ANY q), so the
    # rendered source itself must pin coefficients and signs.
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x-2*x*y"], field=GF(32003))
    s = scripts.graded_dims_script(A, top=4).replace(" ", "")
    assert "y*x-2*x*y" in s
    B = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y+y*x"], field=GF(32003))
    t = scripts.graded_dims_script(B, top=4).replace(" ", "")
    assert "x*y+y*x" in t


def test_non_rational_coefficient_refused_loudly():
    import fractions
    import pytest as _pytest
    from quiverlab.m2.scripts import _relation_to_m2

    class _FakeRel:
        terms = ((object(), ("x", "y")),)     # non-rational coefficient

    with _pytest.raises(QuiverlabError, match="rational"):
        _relation_to_m2(_FakeRel())
