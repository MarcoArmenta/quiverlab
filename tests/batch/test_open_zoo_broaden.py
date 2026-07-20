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
