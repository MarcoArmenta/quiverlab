"""Plan 15: corner-mode checkpoints for deepen (multi-vertex resumable driver).

Before this plan `deepen` refused multi-vertex algebras (NotImplementedError, the
Plan-13 boundary): its checkpoint payload and manual resume-state rebuild had no
corner data.  Corner data (_CornerContext, gens0, rad_ab_pairs, the engine) is
deterministic from (A, prime), so the checkpoint persists only what the stepper
accumulates: cur / cur_r / rks / tags (+ the rolling last_gens and the HH /
per_degree records).

Oracle: minimal_homology_dims (the Plan-13-certified corner engine) — live,
never hardcoded."""
import pytest

import quiverlab as ql
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.deepen import deepen, _load_ckpt
from quiverlab.engine.hh_engine import truncated_polynomial
from quiverlab.engine.resolutions_minimal import minimal_homology_dims
from quiverlab.errors import QuiverlabError

PRIMES = (32003, 2, 3, 5)


def _eng(vertices, arrows, relations, p=32003):
    Q = ql.Quiver(vertices, arrows)
    return to_engine(Q.algebra(relations=relations, field=ql.GF(p)))


def _cn32(p=32003):
    """kZ_3/rad^2: periodic, HH nonzero in every degree (the resume fixture)."""
    return _eng([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)},
                ["a*b", "b*c", "c*a"], p=p)


def test_deepen_corner_matches_minimal_homology_dims(tmp_path):
    """CN(3,2) over four primes: deepen's HH_* == the batch corner engine."""
    for p in PRIMES:
        A = _cn32(p=p)
        out = deepen(A, str(tmp_path / ("ck%d" % p)), prime=p, max_degree=4)
        ref = minimal_homology_dims(A, 4, primes=(p,))[p]
        assert out["HH"] == ref, "p=%d: %s != %s" % (p, out["HH"], ref)


def test_deepen_corner_termination_ka2(tmp_path):
    """kA_2 (hereditary): the corner resolution terminates; deepen must report it
    (stop_reason='terminated', hochschild_dim=1) with HH == the batch engine."""
    A = _eng([1, 2], {"a": (1, 2)}, [])
    out = deepen(A, str(tmp_path / "ck"), prime=32003)
    assert out["stop_reason"] == "terminated" and out["terminated"]
    assert out["hochschild_dim"] == 1
    ref = minimal_homology_dims(A, 1, primes=(32003,))[32003]
    assert out["HH"] == ref


def test_deepen_corner_nonmonomial_square(tmp_path):
    """kQ/(ab - cd) (dim 9, non-monomial multi-vertex): HH_0 = 4 (Plan-13 pin)."""
    A = _eng([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)},
             ["a*b - c*d"])
    out = deepen(A, str(tmp_path / "ck"), prime=32003, max_degree=2)
    ref = minimal_homology_dims(A, 2, primes=(32003,))[32003]
    assert out["HH"] == ref[:len(out["HH"])]
    assert out["HH"][0] == 4


def test_deepen_corner_checkpoint_has_tags(tmp_path):
    """The corner checkpoint persists per-degree tags, consistent with rks."""
    A = _cn32()
    ck = str(tmp_path / "ck")
    deepen(A, ck, prime=32003, max_degree=3)
    payload = _load_ckpt(ck)
    assert payload is not None and "tags" in payload
    for n_, tgs in payload["tags"].items():
        assert len(tgs) == payload["rks"][n_]
        assert all(len(tg) == 2 for tg in tgs)


def test_deepen_corner_memory_wall(tmp_path):
    """A 1-byte transient budget stops with stop_reason='memory' at degree 1 and
    the same (empty) exact HH prefix as the batch engine under the same budget."""
    A = _cn32()
    out = deepen(A, str(tmp_path / "ck"), prime=32003, max_transient_bytes=1)
    assert out["stop_reason"] == "memory"
    assert out["wall_radK_bytes"] > 1 and out["wall_degree"] == 1
    ref = minimal_homology_dims(A, 6, primes=(32003,), max_transient_bytes=1)[32003]
    assert out["HH"] == ref


def test_deepen_corner_finalize_only(tmp_path):
    """finalize_only re-emits the summary from the corner checkpoint unchanged."""
    A = _cn32()
    ck = str(tmp_path / "ck")
    full = deepen(A, ck, prime=32003, max_degree=3)
    fin = deepen(A, ck, prime=32003, finalize_only=True)
    assert fin["stop_reason"] == "checkpoint"
    assert fin["HH"] == full["HH"]
    assert fin["max_degree_reached"] == full["max_degree_reached"]


def test_deepen_corner_resume(tmp_path):
    """A stop-early run then a continue run == one full fresh run == the oracle."""
    A = _cn32()
    ck = str(tmp_path / "ck")
    deepen(A, ck, prime=32003, max_degree=2)                # stop early
    out = deepen(A, ck, prime=32003, max_degree=5)          # resume -> continue
    ref = minimal_homology_dims(A, 5, primes=(32003,))[32003]
    assert out["HH"] == ref
    assert out["max_degree_reached"] >= 3                   # actually went past the cap
    fresh = deepen(A, str(tmp_path / "ck2"), prime=32003, max_degree=5)
    assert out["HH"] == fresh["HH"]


def test_deepen_mode_mismatch_refuses(tmp_path):
    """A corner algebra resumed against a local checkpoint dir (and vice versa)
    must refuse loudly -- a silent overlay would corrupt the resolution state."""
    local, corner = truncated_polynomial(3), _cn32()
    ck_local = str(tmp_path / "ck_local")
    deepen(local, ck_local, prime=32003, max_degree=2)
    with pytest.raises(QuiverlabError):
        deepen(corner, ck_local, prime=32003, max_degree=3)
    ck_corner = str(tmp_path / "ck_corner")
    deepen(corner, ck_corner, prime=32003, max_degree=2)
    with pytest.raises(QuiverlabError):
        deepen(local, ck_corner, prime=32003, max_degree=3)
