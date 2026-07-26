"""Plan 30 -- schema v2 additions: the `tor` compute kind + `tor_target` (a LEFT
A-module, side defaulting to "left"), and the `decompose` kind (no target).

The load-bearing invariant is cache-key discipline (Plan 25/26): an absent
`tor_target` DROPS from `model_dump`, so every pre-Plan-30 request -- family,
quiver, and the Plan-26 module/ext requests -- canonicalizes byte-identically.
A frozen-literal pin below is the regression fence.
"""
import pytest

from webapp.server.cache import canonical_key
from webapp.server.schema import (
    MODULE_KINDS, MODULE_RANGE_KINDS, ComputeRequest,
)

_V = "0.1.0.dev0"

_QUIVER = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
           "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}}
_A2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
       "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}}


def _req(**over):
    body = {"schema": 2, "algebra": _QUIVER, "compute": ["dimension_vector"],
            "artifacts": {"pdf": False, "tikz": False},
            "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}}}
    body.update(over)
    return ComputeRequest.model_validate(body)


# --------------------------------------------------------------------------- #
# The new kinds are registered as module (range) kinds
# --------------------------------------------------------------------------- #

def test_tor_and_decompose_are_module_kinds():
    assert {"tor", "decompose"} <= MODULE_KINDS
    assert "tor" in MODULE_RANGE_KINDS          # tor:0..n
    assert "decompose" not in MODULE_RANGE_KINDS  # scalar


def test_decompose_kind_parses_with_module_only():
    req = _req(compute=["decompose"])
    assert req.tor_target is None and req.ext_target is None


# --------------------------------------------------------------------------- #
# tor_target: a LEFT module (side defaults to "left"; "right" is rejected)
# --------------------------------------------------------------------------- #

def test_tor_target_side_defaults_to_left():
    req = _req(compute=["tor:0..3"], algebra=_A2,
               module={"dims": {"1": 1, "2": 1}, "maps": {"a": [[1]]}},
               tor_target={"dims": {"1": 1}})           # no side given
    assert req.tor_target.side == "left"

    req2 = _req(compute=["tor:0..3"], algebra=_A2,
                module={"dims": {"1": 1, "2": 1}, "maps": {"a": [[1]]}},
                tor_target={"builtin": {"kind": "simple", "vertex": 1}})
    assert req2.tor_target.side == "left"               # builtin, no side -> left


def test_tor_target_explicit_left_is_accepted():
    req = _req(compute=["tor:0..2"], algebra=_A2,
               module={"dims": {"1": 1}}, tor_target={"dims": {"1": 1}, "side": "left"})
    assert req.tor_target.side == "left"


def test_tor_target_right_is_rejected():
    with pytest.raises(ValueError):
        _req(compute=["tor:0..2"], algebra=_A2, module={"dims": {"1": 1}},
             tor_target={"dims": {"1": 1}, "side": "right"})


def test_tor_target_builtin_right_is_rejected():
    with pytest.raises(ValueError):
        _req(compute=["tor:0..2"], algebra=_A2, module={"dims": {"1": 1}},
             tor_target={"builtin": {"kind": "simple", "vertex": 1, "side": "right"}})


def test_tor_requires_a_tor_target():
    with pytest.raises(ValueError):
        _req(compute=["tor:0..2"], algebra=_A2, module={"dims": {"1": 1}})


def test_tor_target_requires_schema_2():
    with pytest.raises(ValueError):
        ComputeRequest.model_validate(
            {"schema": 1, "algebra": _A2, "compute": ["cartan"],
             "tor_target": {"dims": {"1": 1}, "side": "left"}})


# --------------------------------------------------------------------------- #
# Cache-key discipline (Plan 25/26): absent tor_target DROPS; keys are byte-stable
# --------------------------------------------------------------------------- #

# Frozen literals: the canonical_key of a pre-Plan-30 request must be UNCHANGED.
# These were computed under _V with the drop-when-absent model_dump; the delegation
# goldens (frozen BEFORE Plan 30) still passing is the independent cross-check.
_FAMILY = {"schema": 1,
           "algebra": {"kind": "family", "family": "QuantumCI",
                       "params": {"q": 1}, "field": {"kind": "GF", "p": 2, "n": 1}},
           "compute": ["hh_cohomology:0..3"],
           "artifacts": {"pdf": False, "tikz": False}}
_FAMILY_KEY = "a073f92d0f94c7306f827a4552543a4c28c9dbf463c405930c0706d56d36c272"

_MOD_EXT = {"schema": 2,
            "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                        "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
            "compute": ["ext:0..2"], "artifacts": {"pdf": False, "tikz": False},
            "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}},
            "ext_target": {"builtin": {"kind": "simple", "vertex": 1}}}
_MOD_EXT_KEY = "7aeaaaaa5781b3c9a2583e8d5ed98c3a46e1766e5044bad53e2ba40f0c25cc99"


def test_old_style_request_key_is_unchanged_against_frozen_literal():
    d = ComputeRequest.model_validate(_FAMILY).model_dump(by_alias=True)
    assert "tor_target" not in d and "ext_target" not in d and "module" not in d
    assert set(d) == {"schema", "algebra", "compute", "artifacts"}
    assert canonical_key(d, _V) == _FAMILY_KEY


def test_plan26_module_request_key_is_unchanged_against_frozen_literal():
    # A Plan-26 module+ext request must ALSO be byte-stable: tor_target drops.
    d = ComputeRequest.model_validate(_MOD_EXT).model_dump(by_alias=True)
    assert "tor_target" not in d
    assert canonical_key(d, _V) == _MOD_EXT_KEY


def test_present_tor_target_participates_in_the_key():
    # When a tor_target IS present, it enters the key (so distinct N's don't collide).
    base = ComputeRequest.model_validate(
        {"schema": 2, "algebra": _A2, "compute": ["tor:0..2"],
         "module": {"dims": {"1": 1}},
         "tor_target": {"dims": {"1": 1}, "side": "left"}}).model_dump(by_alias=True)
    other = ComputeRequest.model_validate(
        {"schema": 2, "algebra": _A2, "compute": ["tor:0..2"],
         "module": {"dims": {"1": 1}},
         "tor_target": {"dims": {"2": 1}, "side": "left"}}).model_dump(by_alias=True)
    assert "tor_target" in base
    assert canonical_key(base, _V) != canonical_key(other, _V)


def test_tor_target_omitted_side_and_explicit_left_collide():
    # The side-default cache invariant for tor_target (mirrors Plan 26 for `side`).
    omitted = ComputeRequest.model_validate(
        {"schema": 2, "algebra": _A2, "compute": ["tor:0..2"],
         "module": {"dims": {"1": 1}},
         "tor_target": {"dims": {"1": 1}}}).model_dump(by_alias=True)
    explicit = ComputeRequest.model_validate(
        {"schema": 2, "algebra": _A2, "compute": ["tor:0..2"],
         "module": {"dims": {"1": 1}},
         "tor_target": {"dims": {"1": 1}, "side": "left"}}).model_dump(by_alias=True)
    assert canonical_key(omitted, _V) == canonical_key(explicit, _V)
