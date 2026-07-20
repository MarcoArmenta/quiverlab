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
