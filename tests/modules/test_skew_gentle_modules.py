"""Skew-gentle module classification via special-string re-gluing (Plan 68; HZZ / GL
Table 1). Self-cert: every materialised module passes check_module (inside
from_arrow_action) and is indecomposable. Cross-engine: on a rep-finite skew-gentle
algebra the number of materialised indecomposables equals the AR-quiver vertex count.
Literature: a special (type-p) string has exactly TWO forms (the +/- split vertices);
an ordinary (u,u) string has one. Over QQ / GF(32003) (decompose char caveat)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.decompose import is_indecomposable
from quiverlab.skewgentle.modules import (classify, skew_gentle_indecomposables,
                                          skew_gentle_module)
from quiverlab.skewgentle.triple import SkewGentleTriple

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


def _t():                              # Q: 1 --a--> 2 , Sp = {2}
    return SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})


@selfcert
def test_materialised_modules_are_indecomposable():
    t = _t()
    for M in skew_gentle_indecomposables(t, field=QQ):
        assert is_indecomposable(M)                  # check_module already passed at build


@lit
def test_special_string_has_two_forms():
    t = _t()
    specials = [s for s in classify(t) if "p" in s.type and not s.is_band]
    assert specials, "expected at least one type-p string over Sp={2}"
    for s in specials:
        assert len(s.forms) == 2                     # the +/- forms = split vertices
    ordinary = [s for s in classify(t) if s.type == ("u", "u")]
    for s in ordinary:
        assert len(s.forms) == 1


@xeng
def test_indecomposable_count_equals_ar_quiver_count():
    # a rep-finite skew-gentle algebra: #materialised indecomposables == #AR vertices.
    from quiverlab.modules.ar import knit_ar_quiver
    from quiverlab.skewgentle.split import SkewGentleAlgebra
    t = _t()
    A = SkewGentleAlgebra(t, field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    mods = skew_gentle_indecomposables(t, field=QQ)
    assert len(mods) == len(ar.vertices)


@selfcert
def test_two_forms_are_non_isomorphic():
    from quiverlab.modules.hom import is_isomorphic
    t = _t()
    s = next(x for x in classify(t) if "p" in x.type and not x.is_band)
    M0 = skew_gentle_module(t, s, form=0, field=QQ)
    M1 = skew_gentle_module(t, s, form=1, field=QQ)
    assert not is_isomorphic(M0, M1)                  # the +/- forms are distinct A-modules
