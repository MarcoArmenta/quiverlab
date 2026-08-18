"""Public Tate-Hochschild surface (Plan 76 / R3): native GF(p) splice, duality dims
(symmetric, any Domain), positive threshold (any Gorenstein).  Bergh-Jorgensen JNCG 7
(2013).

BUCKET NOTE (plan Global Constraints): `tests/hochschild/` is auto-assigned FAST, so
every test here that performs a real GF(p) complete-resolution build carries an
explicit `@pytest.mark.deep` -- an explicit bucket beats the directory default.
"""
import pytest

from quiverlab import (NakayamaAlgebra, QuantumCI, linear_path_algebra,
                       truncated_polynomial)
from quiverlab.errors import QuiverlabError
from quiverlab.fields import GF, QQ

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


# ---------------------------------------------------------------------------
# native route -- the full Z-graded ring
# ---------------------------------------------------------------------------
@lit
@pytest.mark.deep
def test_kx3_native_full_ring():
    A = truncated_polynomial(3, field=GF(32003))
    r = A.tate_hochschild(5, engine="native")
    assert r.pos_dims == (2, 2, 2, 2, 2, 2)          # HHhat^0..5
    assert r.neg_dims == (2, 2, 2, 2, 2)             # HHhat^-1..-5
    assert r.hat_hh0 == 2 and r.agrees_from == 1
    assert r.scope == "self_injective" and r.engine == "native"
    assert r.ordinary_pos[0] == 3                    # HH^0 = dim Z(A) = 3 != HHhat^0


@selfcert
@pytest.mark.deep
def test_kz3j2_native_nu_twist_witness():
    """The nu-twist discriminating witness (M3): kZ3/J2 is self-injective, NON-symmetric,
    vertex permutation pi = (1 2 3).  Positive degrees hit the threshold; the negatives
    are engine-computed and are NOT the symmetric mirror."""
    from quiverlab.hochschild.tate import nakayama_permutation
    A = NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(32003))
    assert A.is_selfinjective() and not A.is_symmetric()
    assert nakayama_permutation(A) == {1: 2, 2: 3, 3: 1}
    r = A.tate_hochschild(4, engine="native")
    hh = A.hochschild_cohomology(4).dims
    assert tuple(r.pos_dims[1:]) == tuple(hh[1:])                 # threshold, m >= 1
    assert any(r.neg_dims[j - 1] != r.pos_dims[j - 1] for j in range(1, 5))


@xeng
@pytest.mark.deep
def test_native_positive_part_matches_shipped_hh():
    """The NEW splice engine against the SHIPPED bar/fast engine, both directions of
    the zoo (single vertex, symmetric multi-vertex, non-symmetric multi-vertex)."""
    for A in (truncated_polynomial(3, field=GF(32003)),
              NakayamaAlgebra(n=2, l=3, cyclic=True, field=GF(32003)),
              NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(32003))):
        r = A.tate_hochschild(4, engine="native")
        hh = A.hochschild_cohomology(4).dims
        assert tuple(r.pos_dims[1:]) == tuple(hh[1:]), A
        assert r.ordinary_pos == tuple(hh)


@selfcert
@pytest.mark.deep
def test_auto_picks_native_over_gfp_selfinjective():
    A = truncated_polynomial(3, field=GF(32003))
    assert A.tate_hochschild(3).engine == "native"


# ---------------------------------------------------------------------------
# positive route -- universal, and honest about degree 0
# ---------------------------------------------------------------------------
@xeng
def test_positive_threshold_universal_QQ():
    """HHhat^{m>=1} = HH^m for any Gorenstein A over any Domain.  Degree 0 is None on
    the positive route -- HHhat^0 is native-only and is NOT HH^0 (M2)."""
    A = truncated_polynomial(3, field=QQ)
    r = A.tate_hochschild(6, engine="positive")
    hh = A.hochschild_cohomology(6).dims
    assert r.pos_dims[0] is None
    assert r.pos_dims[1:] == tuple(hh[1:])
    assert r.ordinary_pos == tuple(hh)
    assert r.neg_dims is None and r.agrees_from == 1 and r.hat_hh0 is None


@lit
def test_qci_positive_threshold_QQ():
    """QuantumCI(2,2,2) over QQ (q not a root of unity): HHhat^{n>=1} = HH^n =
    [2,1,0,0,0].  Bergh-Jorgensen compute the full ring as 1,2,1 in degrees 0,1,2 and
    0 elsewhere, so HHhat^0 = 1 while HH^0 = 2 -- the surface reports HH^0 only in
    `ordinary_pos` and never claims it as the Tate degree 0."""
    A = QuantumCI(2, 2, 2, field=QQ)
    r = A.tate_hochschild(5, engine="positive")
    assert r.pos_dims[0] is None
    assert r.pos_dims[1:] == (2, 1, 0, 0, 0)
    assert r.ordinary_pos == (2, 2, 1, 0, 0, 0)
    assert r.neg_dims is None


# ---------------------------------------------------------------------------
# duality route
# ---------------------------------------------------------------------------
@xeng
def test_duality_route_degree0_absent_QQ():
    """Symmetric k[x]/(x^3) over QQ: m >= 1 is HH^m and m <= -2 is HH^{|m|-1}; degrees
    0 and -1 are a closed 2-cycle the duality cannot reach -> None, never HH^0."""
    A = truncated_polynomial(3, field=QQ)
    r = A.tate_hochschild(5, engine="duality")
    hh = A.hochschild_cohomology(5).dims                       # [3,2,2,2,2,2]
    assert r.pos_dims[0] is None and r.hat_hh0 is None
    assert r.pos_dims[1:] == tuple(hh[1:])
    assert r.neg_dims[0] is None
    assert r.neg_dims[1:] == tuple(hh[j - 1] for j in range(2, 6))


@xeng
@pytest.mark.deep
def test_duality_agrees_with_native_where_both_speak():
    """Two independent routes on the same symmetric algebra over the same GF(p): the
    duality theorem against the complete-resolution engine, in every degree the
    duality route claims."""
    A = truncated_polynomial(3, field=GF(32003))
    dua = A.tate_hochschild(4, engine="duality")
    nat = A.tate_hochschild(4, engine="native")
    for m in range(1, 5):
        assert dua.pos_dims[m] == nat.pos_dims[m]
    for j in range(2, 5):
        assert dua.neg_dims[j - 1] == nat.neg_dims[j - 1]


@selfcert
def test_qci_duality_refused():
    """QuantumCI is Frobenius but NOT symmetric, so the duality hypothesis fails."""
    A = QuantumCI(2, 2, 2, field=QQ)
    assert not A.is_symmetric()
    with pytest.raises(QuiverlabError):
        A.tate_hochschild(4, engine="duality")


# ---------------------------------------------------------------------------
# refusals
# ---------------------------------------------------------------------------
@selfcert
def test_native_non_selfinjective_deferred():
    """M4: kA2 is Gorenstein (gl.dim 1) but NOT self-injective, so the D(P_n)-dual
    splice cannot build its complete resolution -> loud DEFERRED refusal."""
    A = linear_path_algebra(2, field=GF(32003))
    assert not A.is_selfinjective() and A.is_gorenstein()
    with pytest.raises(QuiverlabError, match="DEFERRED"):
        A.tate_hochschild(4, engine="native")


@selfcert
def test_gorenstein_non_selfinjective_positive_is_honest():
    """The same algebra on `auto` does NOT die -- it falls to the positive route, which
    then declines to name an agreement degree (the Gorenstein dimension of A^e is not
    certified) rather than guessing one."""
    A = linear_path_algebra(2, field=GF(32003))
    r = A.tate_hochschild(3)
    assert r.engine == "positive" and r.scope == "positive_only"
    assert r.agrees_from is None and r.neg_dims is None
    assert all(d is None for d in r.pos_dims)
    assert r.ordinary_pos == tuple(A.hochschild_cohomology(3).dims)
    assert "not certified" in r.note


@selfcert
def test_native_over_QQ_refused():
    A = truncated_polynomial(3, field=QQ)
    with pytest.raises(QuiverlabError, match="GF"):
        A.tate_hochschild(3, engine="native")


@selfcert
def test_unknown_engine_refused():
    A = truncated_polynomial(2, field=GF(32003))
    with pytest.raises(QuiverlabError, match="unknown"):
        A.tate_hochschild(2, engine="tate")


@selfcert
def test_negative_top_refused():
    A = truncated_polynomial(2, field=GF(32003))
    with pytest.raises(QuiverlabError):
        A.tate_hochschild(-1)


# ---------------------------------------------------------------------------
# eventual periodicity
# ---------------------------------------------------------------------------
@selfcert
@pytest.mark.deep
@pytest.mark.parametrize("build,expect", [
    (lambda: truncated_polynomial(3, field=GF(32003)), 2),
    (lambda: truncated_polynomial(2, field=GF(32003)), 2),
    (lambda: NakayamaAlgebra(n=2, l=3, cyclic=True, field=GF(32003)), None),
])
def test_periodicity_certificate(build, expect):
    from quiverlab.hochschild.tate import tate_periodicity
    cert = tate_periodicity(build())
    assert cert.eventually_periodic is True and cert.complexity == 1
    assert cert.period is not None
    assert cert.invertible_element_degree == cert.period
    if expect is not None:
        assert cert.period == expect


@selfcert
@pytest.mark.deep
def test_qci_not_periodic():
    """QuantumCI is self-injective (so it IS served) but NOT eventually periodic --
    complexity 2.  The certificate is an honest None, never a guessed period."""
    from quiverlab.hochschild.tate import tate_periodicity
    A = QuantumCI(2, 2, 2, field=GF(32003))
    cert = tate_periodicity(A)
    assert cert.complexity == 2
    assert cert.period is None and cert.eventually_periodic is False
    assert cert.invertible_element_degree is None


@selfcert
@pytest.mark.deep
def test_native_result_carries_the_period():
    A = truncated_polynomial(3, field=GF(32003))
    r = A.tate_hochschild(3, engine="native")
    assert r.period == 2 and r.periodicity_degree == 2


# ---------------------------------------------------------------------------
# the runner block
# ---------------------------------------------------------------------------
@selfcert
@pytest.mark.deep
def test_block_is_json_shaped():
    from quiverlab.hochschild.tate import tate_hochschild_block
    A = truncated_polynomial(3, field=GF(32003))
    b = tate_hochschild_block(A, 3)
    assert b["kind"] == "tate_hochschild" and b["engine"] == "native"
    assert b["pos_dims"] == [2, 2, 2, 2] and b["neg_dims"] == [2, 2, 2]
    assert b["hat_hh0"] == 2 and b["agrees_from"] == 1
    assert b["ordinary_pos"][0] == 3
    assert "bergh_jorgensen_tate" in b["references"]
    import json
    json.dumps(b)                                    # no exotic types leak out
