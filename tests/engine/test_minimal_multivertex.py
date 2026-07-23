"""Plan 13: the minimal A^e engine on MULTI-VERTEX algebras (corner-typed projective
resolution). Before this plan the engine silently returned a zero resolution
(rks 1,0,0,...) on any multi-vertex input — its radical formula was local-only and the
fake rad(A^e) swallowed the whole kernel.

Oracles: the normalized bar complex (live, never hardcoded — the same engine bar that
tests/engine/test_multivertex_engine.py certifies) plus theory pins: corner Betti
numbers of a monomial algebra are Bardzell's chain counts |AP^n|."""
import numpy as np
import pytest

import quiverlab as ql
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.hh_engine import hochschild_homology_dims
from quiverlab.engine.resolutions_minimal import (
    minimal_homology_dims, minimal_resolution, radical_basis)
from quiverlab.errors import QuiverlabError

PRIMES = (32003, 2, 3, 5)


def _eng(vertices, arrows, relations, p=32003):
    Q = ql.Quiver(vertices, arrows)
    return to_engine(Q.algebra(relations=relations, field=ql.GF(p)))


def test_ka2_matches_bar_and_terminates():
    """kA_2 (1->2, hereditary, dim 3): HH == bar over four primes; the corner
    resolution is FINITE: P_0 = two vertex corners, P_1 = the arrow corner, done."""
    for p in PRIMES:
        E = _eng([1, 2], {"a": (1, 2)}, [], p=p)
        mh = minimal_homology_dims(E, 4, primes=(p,))
        bh = hochschild_homology_dims(E, 4, primes=(p,))
        assert mh[p] == bh[p][:len(mh[p])], f"kA_2 p={p}: {mh[p]} != {bh[p]}"
    E = _eng([1, 2], {"a": (1, 2)}, [])
    rks, _cols, _e, trunc = minimal_resolution(E, 4, 32003)
    assert trunc is None
    assert rks[0] == 2 and rks[1] == 1 and rks[2] == 0     # finite, minimal


def test_commutative_square_matches_bar():
    """kQ/(ab - cd), dim 9 (Plan-04 Fixture B): HH_. = [4, 0, 0]."""
    E = _eng([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)},
             ["a*b - c*d"])
    mh = minimal_homology_dims(E, 2, primes=(32003,))[32003]
    bh = hochschild_homology_dims(E, 2, primes=(32003,))[32003]
    assert mh == bh[:len(mh)]
    assert mh == [4, 0, 0]


def test_line_quiver_betti_equals_bardzell_chain_counts():
    """kQ/(abc, cde) on 1->...->6 (dim 16, monomial): the corner Betti numbers are
    Bardzell's |AP^n| = 6 vertices, 5 arrows, 2 relations, 1 straddling overlap
    (abcde — the Plan-12 chain), 0. Independent syzygy-side re-derivation of the
    straddle chain."""
    E = _eng([1, 2, 3, 4, 5, 6],
             {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 5), "e": (5, 6)},
             ["a*b*c", "c*d*e"])
    rks, _cols, _e, trunc = minimal_resolution(E, 5, 32003)
    assert trunc is None
    assert [rks[n] for n in range(5)] == [6, 5, 2, 1, 0]


def test_cyclic_nakayama_matches_bar():
    """kZ_3/rad^2 (dim 6): multi-vertex WITH nonzero HH in every degree (periodic);
    the strongest multi-vertex dim cross-check."""
    for p in PRIMES:
        E = _eng([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)},
                 ["a*b", "b*c", "c*a"], p=p)
        mh = minimal_homology_dims(E, 3, primes=(p,))
        bh = hochschild_homology_dims(E, 3, primes=(p,))
        assert mh[p] == bh[p][:len(mh[p])], f"CN(3,2) p={p}: {mh[p]} != {bh[p]}"


def test_non_path_basis_raises_loudly():
    """k x k presented with the NON-path basis {1, u=(1,0)}: u is idempotent, so the
    'non-vertex basis vectors span the radical' premise is false. The engine must
    refuse loudly (nilpotent-closure guard), never return a silent wrong resolution."""
    from quiverlab.engine.hh_engine import Algebra as EngineAlgebra
    T = np.zeros((2, 2, 2), dtype=np.int64)
    T[0, 0, 0] = 1                     # 1*1 = 1
    T[0, 1, 1] = T[1, 0, 1] = 1        # 1*u = u*1 = u
    T[1, 1, 1] = 1                     # u*u = u  (idempotent!)
    E = EngineAlgebra(2, T, [1, 0], name="k x k, non-path basis")
    with pytest.raises(QuiverlabError):
        radical_basis(E, 32003)
    with pytest.raises(QuiverlabError):
        minimal_resolution(E, 3, 32003)
