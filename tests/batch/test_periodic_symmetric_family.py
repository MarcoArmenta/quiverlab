"""Re-port of the bank tests/test_periodic_symmetric_family.py (LEDGER OBLIGATION).
The family k<x,y>/(x^3, y^b - x^2, yx + xy): periodic, symmetric. Fixture Z1."""
import importlib
import os

import numpy as np
import pytest

pytestmark = [
    pytest.mark.skipif(
        importlib.util.find_spec("quiverlab.resolutions_cs") is None,
        reason="periodic-symmetric family needs the Plan 04 CS/reduction-system backend"),
    pytest.mark.skipif(
        os.environ.get("QLAB_RUN_HEAVY_GOLDENS") != "1",
        reason="needs ~2.5GB+ RSS; run with QLAB_RUN_HEAVY_GOLDENS=1 on a machine "
               "with memory headroom"),
]

from quiverlab.engine.coxeter import is_frobenius, nakayama_automorphism
from quiverlab.engine.resolutions_minimal import minimal_homology_dims
from quiverlab.families.zoo import build_from_record, load_catalog


def _algebra(name):
    rec = next(r for r in load_catalog() if r["name"] == name)
    from quiverlab.engine.adapter import to_engine
    from quiverlab.fields import GF
    return to_engine(build_from_record(rec, field=GF(32003)).unit_adapted())


def test_dim9_sibling_periodic_symmetric_and_p2_growth():
    A = _algebra("open2_33_712")                       # k<x,y>/(x^3, y^3 - x^2, yx + xy)
    assert A.m == 9
    assert is_frobenius(A, 32003)
    nu, _ = nakayama_automorphism(A, 32003)
    assert np.array_equal(nu % 32003, np.eye(A.m, dtype=object) % 32003)
    assert minimal_homology_dims(A, 6, primes=(32003,))[32003] == [6, 5, 5, 5, 5, 5, 5]
    hh2 = minimal_homology_dims(A, 6, primes=(2,))[2]
    assert hh2 == [9, 10, 14, 18, 22, 26, 30]          # growing tail at p=2
    assert hh2[-1] > hh2[0]


def test_dim21_headline_is_symmetric():
    A = _algebra("open2_37_19612")                     # k<x,y>/(x^3, y^7 - x^2, yx + xy)
    assert A.m == 21
    assert is_frobenius(A, 32003)
    nu, _ = nakayama_automorphism(A, 32003)
    assert np.array_equal(nu % 32003, np.eye(A.m, dtype=object) % 32003)
