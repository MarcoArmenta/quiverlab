"""Plan 35 wave 3d: the instant tier skips the explicit-HH-representatives capture.

The plain hh_cohomology / hh_homology reps run the GF(p) route through
``engine.tt_calculus``, whose cold numba JIT is paid fresh in every spawned instant
child and runs tens of seconds -- over the instant wall net -- while the reps only feed
the report / GUI, which the instant tier discards. So ``run_spec(..., capture_reps=
False)`` (what ``webapp.server.instant`` passes) must ship the dims block WITHOUT the
reps, byte-identical to the pre-wave-3d block; the default (queued / saved path) keeps
the reps.
"""
import pathlib

_BODY = {"schema": 2,
         "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                     "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
         "compute": ["hh_cohomology:0..3", "hh_homology:0..3"],
         "artifacts": {"pdf": False, "tikz": False}}
_REPS_KEYS = ("basis_classes", "chain_basis", "differentials", "inner_dims")


def _run(tmp_path, capture_reps):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    return run_spec(ComputeRequest.model_validate(_BODY),
                    pathlib.Path(tmp_path), capture_reps=capture_reps)


def test_instant_omits_hh_reps(tmp_path):
    res = _run(tmp_path, capture_reps=False)
    for kind in ("hh_cohomology", "hh_homology"):
        blk = res["results"][kind]
        assert blk["dims"]                                   # the dims still ship
        for key in _REPS_KEYS:
            assert key not in blk, (kind, key)               # but no reps


def test_full_run_carries_hh_reps(tmp_path):
    res = _run(tmp_path, capture_reps=True)
    for kind in ("hh_cohomology", "hh_homology"):
        blk = res["results"][kind]
        for key in _REPS_KEYS:
            assert key in blk, (kind, key)


def test_dims_are_identical_across_the_flag(tmp_path):
    """The flag changes ONLY whether reps are attached -- never the dims themselves."""
    off = _run(tmp_path / "a", capture_reps=False)
    on = _run(tmp_path / "b", capture_reps=True)
    for kind in ("hh_cohomology", "hh_homology"):
        assert off["results"][kind]["dims"] == on["results"][kind]["dims"]
