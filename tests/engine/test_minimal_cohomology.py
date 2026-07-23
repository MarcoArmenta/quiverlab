"""Plan 16: HH cohomology from the minimal/corner A^e resolution (Hom-collapse).

The minimal engine was homology-only; Hom_{A^e}(-, A) on the SAME resolution
gives deep HH^. for any f.d. algebra. The coh-side collapse acts a.w.b (the
homology side is b.w.a), and on the corner path the block of a generator tagged
(v, w) is e_v A e_w -- the OPPOSITE corner of the homology target e_w A e_v.

Oracles: the dual normalized bar complex (scan3.hochschild_cohomology_dims,
live), the CS coh side (second deep engine), and the Happel/Kunneth [1,0,0]
pins."""
import pytest

import quiverlab as ql
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.hh_engine import truncated_polynomial
from quiverlab.engine.scan3 import quantum_ci, hochschild_cohomology_dims
from quiverlab.engine.resolutions_minimal import minimal_cohomology_dims

PRIMES = (32003, 2, 3, 5)


def _eng(vertices, arrows, relations, p=32003):
    Q = ql.Quiver(vertices, arrows)
    return to_engine(Q.algebra(relations=relations, field=ql.GF(p)))


def test_local_zoo_matches_bar():
    """k[x]/x^3 and quantum_ci(2) over four primes: minimal HH^. == bar HH^."""
    for make in (lambda: truncated_polynomial(3), lambda: quantum_ci(2)):
        A = make()
        for p in PRIMES:
            mc = minimal_cohomology_dims(A, 4, primes=(p,))[p]
            bc = hochschild_cohomology_dims(A, 4, primes=(p,))[p]
            assert mc == bc[:len(mc)], f"{A.name} p={p}: {mc} != {bc}"


def test_deep_cross_oracle_vs_cs():
    """Past the bar window: minimal coh == CS coh degreewise to depth 8 on the
    quantum CI over GF(3) -- two INDEPENDENT deep engines agreeing."""
    from quiverlab.resolutions_cs.homology import cs_cohomology_dims
    Q = ql.Quiver([1], {"x": (1, 1), "y": (1, 1)})
    A = Q.algebra(relations=["x*x", "y*y", "x*y + y*x"], field=ql.GF(3))
    cs = cs_cohomology_dims(A, 8).dims
    mc = minimal_cohomology_dims(to_engine(A), 8, primes=(3,))[3]
    assert mc == cs


def test_truncation_semantics():
    """A tiny term cap truncates: the dims list is SHORTER and an exact prefix,
    never padded or wrong (delta^t needs the unknown d_{t+1})."""
    A = quantum_ci(2)
    full = minimal_cohomology_dims(A, 6, primes=(32003,))[32003]
    trunc = minimal_cohomology_dims(A, 6, primes=(32003,), max_term_dim=40)[32003]
    assert len(trunc) < len(full)
    assert trunc == full[:len(trunc)]
