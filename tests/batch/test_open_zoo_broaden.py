"""Re-port of the bank tests/test_open_zoo_broaden.py (LEDGER OBLIGATION).
Curated-catalog dimension bands + spec generation. Import-closed subset."""
from quiverlab.batch import depth_for_dim, max_term_dim_for_dim, open_zoo_to_specs
from quiverlab.families.zoo import load_catalog


def test_catalog_round_trips_and_covers_dimension_bands():
    cat = load_catalog()
    for e in cat:
        assert {"name", "ngen", "dim", "rules"} <= set(e)
    dims = {e["dim"] for e in cat}
    assert {9, 12}.issubset(dims) or len(dims) >= 5      # curated band coverage


def test_open_zoo_to_specs_band_and_probe():
    specs = open_zoo_to_specs(load_catalog(), min_dim=9, max_dim=9)
    assert specs
    for s in specs:
        assert s["builder"] == "reduction_system"
        assert s["N"] == depth_for_dim(9)
        assert s["max_term_dim"] == max_term_dim_for_dim(9)
    # limit branch: cap after filtering -> exactly the band's first entry (bank parity, lines 62-64)
    probe = open_zoo_to_specs(load_catalog(), min_dim=9, max_dim=9, limit=1)
    assert len(probe) == 1 and probe[0] == specs[0]


def test_multivertex_record_flows_through_specs_and_analyze():
    """A multi-vertex catalog record survives the spec adapter and analyze():
    the batch scan surface serves the Plan-18 diversity records (specs carry the
    quiver data; _analyze_open must NOT unit-adapt a multi-vertex algebra -- that
    destroys the path-type basis and trips the Plan-13 radical guard)."""
    from quiverlab import GF
    from quiverlab.batch.scan import analyze
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import minimal_homology_dims
    from quiverlab.families.zoo import build_from_record
    rec = next(r for r in load_catalog() if r.get("name") == "cn_3_2")
    specs = open_zoo_to_specs([rec], primes=(32003,))
    assert len(specs) == 1 and specs[0]["builder"] == "reduction_system"
    out = analyze(specs[0])
    assert "error" not in out, out.get("error")
    assert out["dim"] == 6
    E = to_engine(build_from_record(rec, field=GF(32003)))
    ref = minimal_homology_dims(E, out["N"], primes=(32003,))[32003]
    assert out["HH_homology"]["32003"] == ref
