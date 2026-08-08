"""Skew-gentle module classification via special-string re-gluing (Plan 68; HZZ / GL
Table 1). Self-cert: every materialised module passes check_module (inside
from_arrow_action) and is indecomposable. Cross-engine: the loop-free string census is a
documented STRICT SUBSET of the AR-quiver indecomposables (headline 5 of 6, mesh 8 of 11)
and every census module embeds in the AR list up to iso -- the missing modules are the
loop-traversal / mixed-eigenvalue ones (e.g. the projective P_1). Literature: a special
(type-p) string has exactly TWO forms (the +/- split vertices); an ordinary (u,u) string
has one. The enumeration runs on the char-FREE split model over QQ (a presentation
invariant, M3), so a GF(2) request is answered soundly (the char-2 regression)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import GF, QQ
from quiverlab.modules.decompose import is_indecomposable
from quiverlab.skewgentle.block import skew_gentle_block
from quiverlab.skewgentle.modules import (_string_census_modules, classify,
                                          skew_gentle_indecomposables,
                                          skew_gentle_module)
from quiverlab.skewgentle.triple import SkewGentleTriple

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


def _t():                              # Q: 1 --a--> 2 , Sp = {2}
    return SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})


def _t_mesh():                         # Q: 1 --a--> 2 --b--> 3 , Sp = {2}, a*b in I
    return SkewGentleTriple.make(Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}),
                                 relations=["a*b"], special={2})


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
@pytest.mark.parametrize("triple_fn, n_census, n_ar",
                         [(_t, 5, 6), (_t_mesh, 8, 11)])
def test_string_census_is_a_strict_subset_of_the_ar_indecomposables(triple_fn,
                                                                     n_census, n_ar):
    # De-vacuumed cross-engine oracle (Plan 68 fix round): the OLD test compared the AR
    # route to itself (skew_gentle_indecomposables IS the AR vertices on rep-finite input),
    # so it agreed by construction. The REAL relationship: the loop-free A^g-walk census
    # is a documented STRICT SUBSET of the indecomposables (headline 5 of 6, mesh 8 of 11).
    # The missing modules are the loop-traversal / mixed-eigenvalue ones (e.g. the
    # projective P_1 of the headline) that no loop-free walk produces -- the symmetric-
    # string enumeration that would produce them is not implemented (a backlog item).
    #
    # Two INDEPENDENT engines: the string materialisation (classify + skew_gentle_module)
    # vs the P41 AR knit. The genuine agreement pinned here is the EMBEDDING -- every
    # census module is isomorphic to some AR vertex -- together with the strict-subset
    # counts. Both live over the SAME cached QQ split algebra, so is_isomorphic compares.
    from quiverlab.modules.hom import is_isomorphic
    t = triple_fn()
    ar_mods = skew_gentle_indecomposables(t, field=QQ)   # AR route (complete => authoritative)
    census = _string_census_modules(t)
    assert len(ar_mods) == n_ar                          # AR = authoritative full count
    assert len(census) == n_census                       # loop-free census = strict subset
    assert n_census < n_ar                               # documented STRICT subset
    # every census module embeds in the AR list up to iso (the real cross-engine pin):
    for cm in census:
        assert any(is_isomorphic(cm, am) for am in ar_mods), \
            "a loop-free-census module is missing from the AR enumeration"


def test_char2_block_counts_on_the_char_free_model():
    # M3 REGRESSION (the fix round's MAJOR): over GF(2) the AR knit returns
    # is_complete=False (status='error', the char-2 decompose caveat). The buggy code read
    # that as "rep-infinite" and returned the incomplete loop-free string sample (5), and
    # support_tau_tilting ran the exchange-graph BFS over the UNSOUND GF(2). Fixed: both
    # counts route through the char-FREE split model over QQ (a presentation invariant), so
    # the GF(2) block reports the SAME 6 indecomposables as QQ, WITH an explicit note.
    # (VERIFIED live: the pre-fix GF(2) block reported num_indecomposables == 5.)
    b2 = skew_gentle_block(_t(), field=GF(2))
    bq = skew_gentle_block(_t(), field=QQ)
    cl2 = b2["classification"]
    assert cl2["num_indecomposables"] == 6                      # char-free, NOT the 5-sample
    assert cl2["num_indecomposables"] == bq["classification"]["num_indecomposables"]
    assert b2["tau_tilting"]["num_pairs"] == bq["tau_tilting"]["num_pairs"]   # both = 14
    assert "char" in cl2["note"].lower() and "qq" in cl2["note"].lower()      # honest note
    # the module-level function agrees directly:
    assert len(skew_gentle_indecomposables(_t(), field=GF(2))) == 6


@selfcert
def test_two_forms_are_non_isomorphic():
    from quiverlab.modules.hom import is_isomorphic
    t = _t()
    s = next(x for x in classify(t) if "p" in x.type and not x.is_band)
    M0 = skew_gentle_module(t, s, form=0, field=QQ)
    M1 = skew_gentle_module(t, s, form=1, field=QQ)
    assert not is_isomorphic(M0, M1)                  # the +/- forms are distinct A-modules


@selfcert
def test_named_string_materialises_to_the_expected_module():
    # A direct materialisation self-cert (beyond the pairwise non-iso check): take the
    # named special string M(a^-1) on the headline (type (p, u), two +/- forms) and pin its
    # exact dimension vector on the split quiver {1, 2+, 2-}. Form 0 lands the special end
    # on the 2+ copy, form 1 on 2-; both are 2-dimensional and self-certify (check_module).
    t = _t()
    s = next(x for x in classify(t)
             if len(x.walk) == 1 and "p" in x.type and not x.is_band)   # walk (a, -1)
    M0 = skew_gentle_module(t, s, form=0, field=QQ)
    M1 = skew_gentle_module(t, s, form=1, field=QQ)
    for M in (M0, M1):
        ok, witness = M.check_module()
        assert ok, witness                                              # genuine A-module
        assert M.dim == 2
    assert M0.dimension_vector() == {"1": 1, "2+": 1, "2-": 0}          # + copy
    assert M1.dimension_vector() == {"1": 1, "2+": 0, "2-": 1}          # - copy
