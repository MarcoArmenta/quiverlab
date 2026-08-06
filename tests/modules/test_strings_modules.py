"""String and band modules (Plan 46 / C5).

Self-cert: every materialised module passes check_module (inside from_arrow_action)
and is indecomposable; the simple / projective string modules match the builtin ones.
Loud: a band eigenvalue outside the field. Over QQ / GF(32003) so is_indecomposable
decides (char > dim).

PLAN CORRECTIONS (documented): (1) the plan's peak walk (a, b^-1) over gentle
kA3/(ab) is NOT composable (b^-1 starts at vertex 3, not 2), and gentle kA3/(ab) has
no peaks at all -- a genuine peak string is tested on the quiver 1 -a-> 2 <-b- 3;
(2) the band example uses the Kronecker quiver (the 2-cycle has no band -- see
test_strings_walks); (3) the loud-eigenvalue token is "1/5" (a rational not in
GF(5)) rather than "1/2" (which IS 3 in GF(5))."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.decompose import is_indecomposable
from quiverlab.modules.hom import is_isomorphic
from quiverlab.strings.modules import band_module, string_module
from quiverlab.strings.walks import find_bands

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _gentle_a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)


def _peak():
    # 1 -a-> 2 <-b- 3, no relations: gentle, has the peak string a b^-1 (the "V").
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (3, 2)}).algebra(
        relations=[], field=QQ)


def _kronecker(field=QQ):
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(
        relations=[], field=field)


def _some_band(A):
    bands = find_bands(A, max_length=6)
    assert bands, "expected at least one band"
    return bands[0]


@selfcert
def test_string_modules_are_indecomposable():
    A = _gentle_a3()
    for w in [((None, 1),), (("a", +1),), (("b", +1),)]:
        assert is_indecomposable(string_module(A, w))
    # a genuine peak string a b^-1 (direct then inverse) at vertex 2:
    P = _peak()
    peak = string_module(P, (("a", +1), ("b", -1)))
    assert peak.dimension_vector() == {1: 1, 2: 1, 3: 1}
    assert is_indecomposable(peak)


@xeng
def test_length_one_direct_string_is_the_projective_arm():
    # the walk (a) over gentle kA3/(ab): the module [1;2] = P1 (dim vector {1:1,2:1}).
    A = _gentle_a3()
    M = string_module(A, (("a", +1),))
    assert M.dimension_vector() == {1: 1, 2: 1, 3: 0}
    assert is_isomorphic(M, A.projective(1))


@xeng
def test_trivial_walk_is_the_simple():
    A = _gentle_a3()
    for v in (1, 2, 3):
        assert is_isomorphic(string_module(A, ((None, v),)), A.simple(v))


@selfcert
def test_band_module_dim_and_indecomposable():
    A = _kronecker(field=QQ)
    band = _some_band(A)
    M = band_module(A, band, eigenvalue=1, mult=1)
    assert M.dim == len(band)                          # mult=1
    assert is_indecomposable(M)


@selfcert
def test_band_eigenvalue_must_be_in_the_field():
    A = _kronecker(field=GF(5))
    band = _some_band(A)
    with pytest.raises(QuiverlabError):
        band_module(A, band, eigenvalue="1/5", mult=1)  # 1/5 is not in GF(5)


@selfcert
def test_band_eigenvalue_must_be_nonzero():
    A = _kronecker(field=QQ)
    band = _some_band(A)
    with pytest.raises(QuiverlabError):
        band_module(A, band, eigenvalue=0, mult=1)
