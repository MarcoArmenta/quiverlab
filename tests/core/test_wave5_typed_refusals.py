"""Wave 5 (v1.0.1): raw exceptions on the public surface become typed refusals.

Every finding here was a raw ``KeyError``/``TypeError``/``IndexError``/``AttributeError``
escaping a documented public path -- the bug class Wave 5 fixes. Each test pins the
NEW typed behaviour (a ``QuiverlabError`` subclass) or the newly-accepted input.
"""
import pytest

import quiverlab
from quiverlab import fields
from quiverlab.errors import QuiverlabError

selfcert = pytest.mark.oracle_selfcert


# --------------------------------------------------------------------------- #
# 5a. The five recognizer predicates must not leak a raw KeyError when the
#     ladder honestly refuses (over GF(2), and on local/self-injective CC).
# --------------------------------------------------------------------------- #
_PREDS = ("is_shod", "is_weakly_shod", "is_laura", "is_ada", "is_quasi_tilted")


@selfcert
def test_recognizer_predicates_typed_over_gf2():
    A = quiverlab.linear_path_algebra(3, field=quiverlab.GF(2))
    for m in _PREDS:
        with pytest.raises(QuiverlabError) as ei:
            getattr(A, m)()
        # not a raw KeyError, and the honest ladder note is carried through
        assert not isinstance(ei.value, KeyError)
        assert "char 2" in str(ei.value) or "did not complete" in str(ei.value)


@selfcert
def test_recognizer_predicates_typed_on_selfinjective_cc():
    # dual numbers k[x]/(x^2) over CC: self-injective => status 'unsupported'
    duals = quiverlab.truncated_polynomial(2)
    for m in _PREDS:
        with pytest.raises(QuiverlabError):
            getattr(duals, m)()


@selfcert
def test_recognizer_three_valued_escape_hatch_available():
    # The user can still tell "provably not SHOD" (False) from "cannot decide"
    # via the rich record -- exactly the koszul_profile / YonedaPresentation.koszul
    # house style (three-valued lives on the record, not on the bool predicate).
    A = quiverlab.linear_path_algebra(3, field=quiverlab.GF(2))
    L = A.recognizer_ladder()
    assert L.is_complete is False
    assert L.status in ("error", "unsupported", "budget")
    assert L.note  # a human-readable reason


@selfcert
def test_recognizer_predicates_still_return_bool_when_decidable():
    # DO NOT BREAK: a rep-finite QQ algebra still returns plain booleans.
    A = quiverlab.linear_path_algebra(3, field=fields.QQ)  # hereditary kA_3
    assert A.is_shod() is True
    assert A.is_quasi_tilted() is True
    assert isinstance(A.is_ada(), bool)


# --------------------------------------------------------------------------- #
# 5c(i). fields.QQ must be exported at top level next to CC/GF/QQi/E.
# --------------------------------------------------------------------------- #
@selfcert
def test_QQ_exported_at_top_level():
    assert hasattr(quiverlab, "QQ")
    assert quiverlab.QQ is fields.QQ
    assert "QQ" in quiverlab.__all__


# --------------------------------------------------------------------------- #
# 5d.1. A.crosscheck() with its own documented default what='hochschild' must
#       not leak a raw TypeError about a missing positional argument.
# --------------------------------------------------------------------------- #
@selfcert
def test_crosscheck_missing_top_is_typed():
    A = quiverlab.linear_path_algebra(3)
    with pytest.raises(QuiverlabError) as ei:
        A.crosscheck()  # what defaults to 'hochschild', top omitted
    assert not isinstance(ei.value, TypeError)
    assert "top" in str(ei.value)


# --------------------------------------------------------------------------- #
# 5d.2. Quiver([], {}).algebra() (empty quiver) must be a typed refusal, not
#       a raw IndexError.
# --------------------------------------------------------------------------- #
@selfcert
def test_empty_quiver_algebra_typed():
    with pytest.raises(QuiverlabError) as ei:
        quiverlab.Quiver([], {}).algebra()
    assert not isinstance(ei.value, IndexError)
    assert "vertex" in str(ei.value).lower()


# --------------------------------------------------------------------------- #
# 5d.3. A.module((1,1,1), ...) with a positional/tuple dimension vector must be
#       accepted (vertex-ordered), not a raw AttributeError.
# --------------------------------------------------------------------------- #
@selfcert
def test_module_accepts_tuple_dimvec_single_vertex():
    A = quiverlab.truncated_polynomial(2, field=fields.QQ)  # k[x]/(x^2), 1 vertex
    maps = {"x": [[0, 0], [1, 0]]}
    M_tuple = A.module((2,), maps)
    M_dict = A.module({1: 2}, maps)
    assert M_tuple.dimension_vector() == M_dict.dimension_vector()


@selfcert
def test_module_accepts_tuple_dimvec_multivertex():
    A = quiverlab.linear_path_algebra(2, field=fields.QQ)  # kA_2, verts [1,2], arrow a1:1->2
    maps = {"a1": [[0, 0], [1, 0]]}
    M_tuple = A.module((1, 1), maps)
    M_dict = A.module({1: 1, 2: 1}, maps)
    assert M_tuple.dimension_vector() == M_dict.dimension_vector() == {1: 1, 2: 1}


@selfcert
def test_module_tuple_wrong_length_typed():
    A = quiverlab.linear_path_algebra(2, field=fields.QQ)
    with pytest.raises(QuiverlabError) as ei:
        A.module((1, 1, 1), {"a1": [[0, 0], [1, 0]]})
    assert not isinstance(ei.value, AttributeError)
    assert "vertice" in str(ei.value).lower() or "vertex" in str(ei.value).lower()
