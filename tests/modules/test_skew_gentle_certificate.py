"""brick-finite <=> representation-finite for skew-gentle algebras (Plan 68).
Certificate = DIJ (brick-finite <=> tau-tilting-finite) o Garcia-Lavoue Thm 3.1
(rep-finite <=> brick-finite, char != 2). Cross-engine: the exchange-graph completeness
route (PRIMARY) and the associated-gentle band-census route (a ONE-SIDED no-false-negative
cross-check on rep-infiniteness -- find_bands(A^g) may miss special bands, W4). Literature:
a rep-finite skew-gentle algebra certifies True; a rep-infinite one (a band) certifies
False.

DEVIATION FROM THE PLAN DOC (documented): the plan's `_rep_infinite_triple` -- the
2-cycle 1<=>2 with I = {ab, ba} and Sp = {1} -- is actually REPRESENTATION-FINITE (its
split is AR-complete with 10 indecomposables, verified live) and additionally trips a
P45 `exchange_graph` `status == 'error'`. It is replaced here with the KRONECKER
`1 =>=> 2` (two parallel arrows, I = empty, Sp = empty), the canonical rep-infinite
gentle algebra, whose band `find_bands(A^g)` detects instantly -- a genuinely
rep-infinite skew-gentle triple, exercising the certificate's rep-infinite path."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.skewgentle.certificate import (brick_finite_certificate,
                                              is_representation_finite)
from quiverlab.skewgentle.triple import SkewGentleTriple

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _rep_finite_triple():             # Q: 1 --a--> 2, Sp = {2}: split is hereditary A-type
    return SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})


def _rep_infinite_triple():           # the Kronecker 1 =>=> 2: a band, rep-infinite (tame)
    Q = Quiver([1, 2], {"a": (1, 2), "c": (1, 2)})
    return SkewGentleTriple.make(Q, relations=[], special=set())


@lit
def test_rep_finite_certifies_true():
    assert is_representation_finite(_rep_finite_triple(), field=QQ) is True


@lit
def test_rep_infinite_certifies_false():
    assert is_representation_finite(_rep_infinite_triple(), field=QQ) is False


@xeng
def test_band_route_is_a_sound_one_sided_check():
    # W4: find_bands(A^g) may be INCOMPLETE for SPECIAL bands, so the band-census route is
    # NOT a full cross-engine AGREEMENT oracle. It is a ONE-SIDED, no-false-negative check
    # on rep-INFINITENESS: a band FOUND certifies rep-infinite (sound), but "no-bands"
    # does NOT on its own certify rep-finite. Route 1 (exchange graph) is the PRIMARY,
    # rigorous (char != 2) decider; route 2 may only DOWNGRADE to rep-infinite.
    for t in (_rep_finite_triple(), _rep_infinite_triple()):
        cert = brick_finite_certificate(t, field=QQ)
        if cert["band_route"] == "bands":
            assert cert["rep_finite"] is False        # found band => genuinely rep-infinite
        # per-instance coincidence for these two hand-picked cases (documented, NOT a
        # general guarantee): the census happens to be complete here.
        assert cert["rep_finite"] == (cert["band_route"] == "no-bands")
        assert cert["tau_tilting_finite"] == cert["rep_finite"]   # char != 2 (route 1)


@lit
def test_char2_scope_is_narrowed_not_lying():
    from quiverlab.fields import GF
    cert = brick_finite_certificate(_rep_finite_triple(), field=GF(2))
    assert cert["char"] == 2 and cert["scope"] == "char==2 (narrowed)"
    # M3: the tau-tilting-finiteness verdict is computed on the char-FREE split model over
    # a decompose-rigorous field (QQ) -- NOT over GF(2). The rep-finite UPGRADE (char != 2,
    # Garcia-Lavoue) is WITHHELD -- never a guessed rep-finite over char 2.
    assert cert["tau_tilting_finite"] in (True, False)   # decided on the QQ model
    assert cert["rep_finite"] is None                    # withheld over char 2 (honest)
