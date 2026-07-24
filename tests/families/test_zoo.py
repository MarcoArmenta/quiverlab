"""zoo(dim_max) iterator over the curated exact zoo (open_zoo lift)."""
import importlib

import pytest

from quiverlab import zoo
from quiverlab.families.zoo import load_catalog
from quiverlab.fields import GF


def test_catalog_is_bundled_and_well_formed():
    cat = load_catalog()
    assert len(cat) >= 10
    for rec in cat:
        assert {"name", "ngen", "dim", "rules"} <= set(rec)
        if "arrows" in rec:                             # multi-vertex record (Plan 18)
            assert rec["ngen"] == len(rec["arrows"])
            assert {v for _n, s, t in rec["arrows"] for v in (s, t)} <= set(rec["vertices"])
        else:
            assert rec["ngen"] in (2, 3)
        for lead, tail in rec["rules"]:                 # exact-integer coefficients only
            assert all(isinstance(i, int) for i in lead)
            for c, w in tail:
                assert isinstance(c, int)               # no floats in the bundled catalog
                assert all(isinstance(i, int) for i in w)


def test_zoo_yields_algebras_up_to_dim_max_in_order():
    algs = list(zoo(dim_max=12))
    dims = [A.dim for A in algs]
    assert dims == sorted(dims)
    assert all(d <= 12 for d in dims)
    assert all(A.dim >= 1 for A in algs)


def test_periodic_symmetric_member_present_dim9():
    names = {getattr(A, "zoo_name", None): A for A in zoo(dim_max=12)}
    assert "open2_33_712" in names
    assert names["open2_33_712"].dim == 9


def test_build_from_record_multivertex_schema():
    """Records may carry vertices + arrows (list order = rules index space);
    legacy records (no arrows) keep the one-vertex loop path."""
    from quiverlab.families.zoo import build_from_record
    rec = {"name": "cn_3_2_inline", "ngen": 3, "dim": 6,
           "vertices": [1, 2, 3],
           "arrows": [["a", 1, 2], ["b", 2, 3], ["c", 3, 1]],
           "rules": [[[0, 1], []], [[1, 2], []], [[2, 0], []]]}
    A = build_from_record(rec)
    assert A.dim == 6                                  # kZ_3/rad^2
    assert A.zoo_name == "cn_3_2_inline"
    nonmono = {"name": "sq_inline", "ngen": 4, "dim": 9,
               "vertices": [1, 2, 3, 4],
               "arrows": [["a", 1, 2], ["b", 2, 4], ["c", 1, 3], ["d", 3, 4]],
               "rules": [[[2, 3], [[1, [0, 1]]]]]}     # cd -> ab  (ab - cd = 0)
    B = build_from_record(nonmono)
    assert B.dim == 9                                  # commutative square


def _is_monomial(rec):
    return all(not tail for _lead, tail in rec["rules"])


def _tip_lengths(rec):
    return {len(lead) for lead, _tail in rec["rules"]}


def test_zoo_diversity_gates():
    """The standing zoo must keep the two shapes whose absence hid the 2026-07-22
    bugs: mixed-length-tip MONOMIAL records (Bardzell straddling) and multi-vertex
    records. Curation that drops them fails here."""
    cat = load_catalog()
    straddle_mono = [r for r in cat if _is_monomial(r) and len(_tip_lengths(r)) > 1]
    multivertex = [r for r in cat if "arrows" in r]
    assert len(straddle_mono) >= 2, "zoo lost its straddling monomial records"
    assert len(multivertex) >= 3, "zoo lost its multi-vertex records"


# dim <= 9 records certify vs the bar oracle (bar is exponential in dim, fine
# here); line_abc_cde (dim 16) is past the bar blow-up and certifies by the
# Plan-13 syzygy pin instead: corner Betti == Bardzell chain counts.
NEW_RECORDS_BAR = ("straddle_xx_yy_xyx", "straddle_xx_yy_yxy", "cn_3_2",
                   "comm_square")


@pytest.mark.parametrize("name", NEW_RECORDS_BAR)
def test_new_records_are_live_certified(name):
    """Each Plan-18 record: recorded dim == built dim, and HH_0..3 via the
    minimal/corner engine == the normalized bar complex over two primes.
    Live oracles -- the records are load-bearing, not dead JSON."""
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    from quiverlab.engine.resolutions_minimal import minimal_homology_dims
    from quiverlab.families.zoo import build_from_record
    rec = next(r for r in load_catalog() if r.get("name") == name)
    for p in (32003, 2):
        A = build_from_record(rec, field=GF(p))
        assert A.dim == rec["dim"], f"{name}: catalog dim {rec['dim']} != built {A.dim}"
        # plain to_engine (NO unit_adapted): the path basis is what the minimal
        # engine's radical identification needs; unit-adaptation is a local-algebra
        # convention that breaks the multi-vertex path-type basis (Plan-13 guard)
        E = to_engine(A)
        mh = minimal_homology_dims(E, 3, primes=(p,))[p]
        bh = hochschild_homology_dims(E, 3, primes=(p,))[p]
        assert mh == bh[:len(mh)], f"{name} p={p}: {mh} != {bh}"


def test_line_abc_cde_certified_by_bardzell_chain_counts():
    """line_abc_cde (dim 16, past the bar blow-up): the corner resolution's Betti
    numbers must equal Bardzell's chain counts 6,5,2,1,0 -- the Plan-12 straddle
    chain `abcde` re-derived by pure syzygy linear algebra (the Plan-13 pin),
    promoted to the standing zoo."""
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import minimal_resolution
    from quiverlab.families.zoo import build_from_record
    rec = next(r for r in load_catalog() if r.get("name") == "line_abc_cde")
    A = build_from_record(rec, field=GF(32003))
    assert A.dim == rec["dim"] == 16
    rks, _cols, _e, trunc = minimal_resolution(to_engine(A), 5, 32003)
    assert trunc is None
    assert [rks[n] for n in range(5)] == [6, 5, 2, 1, 0]


@pytest.mark.skipif(
    importlib.util.find_spec("quiverlab.resolutions_cs") is None,
    reason="open-zone HH golden gated on the Plan 04 CS backend, consistent with Tasks 1/13")
def test_zoo_algebra_hh_matches_open_zone_golden():
    # Build needs the Plan-03 Groebner route (a hard Plan-06 prereq); the depth-16
    # HH is produced by the Plan-02 minimal A^e-resolution engine over GF(32003).
    # open_33_0 has cubic tips (x^3, y^3), which lie OUTSIDE the Chouhy-Solotar
    # quadratic-tip certificate (CS raises NotImplementedError here -- spec-6 risk
    # register), so engine="cs" cannot reach this golden; the minimal engine is its
    # true producer, matching the batch open-zone golden (_analyze_open) and the
    # periodic-symmetric sibling. Gated on CS presence as a coarse "full Plan-04
    # stack present" proxy, consistent with the Task-13 ledger tests.
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import minimal_homology_dims
    A = next(A for A in zoo(dim_max=9, field=GF(32003))
             if getattr(A, "zoo_name", "") == "open_33_0")
    hh = minimal_homology_dims(to_engine(A.unit_adapted()), 16, primes=(32003,))
    assert hh[32003] == [6] + [5] * 16               # golden (Fixture Z1)
