"""Plan 34 (Marco's release feedback) -- the webapp schema-v2 module surface accepts
EXPLICIT-matrix, non-builtin, non-projective/non-injective/non-simple modules for
BOTH arguments of Ext (Plan 26 ``ext_target``) and Tor (Plan 30 ``tor_target``),
the instant tier returns the SAME dims as the direct library call, and a
relation-violating second argument surfaces as a clean typed 4xx (never a 500).

The standing runner tests (``test_runner_modules.py`` / ``test_runner_plan30.py``)
exercise a builtin or a single explicit module against a builtin/simple target;
here BOTH arguments are explicit interior modules over the hereditary kA_5 (interval
[2,3] against [3,4] -- neither is a projective [i,5], injective [1,j], nor simple).

fast bucket (webapp/ directory); no network (direct ``run_spec``).
"""
import pytest

from webapp.server.runner import RunError, run_spec
from webapp.server.schema import ComputeRequest
from webapp.server.security import is_safe_error_type

from quiverlab import GF, linear_path_algebra
from quiverlab.modules.tor import tor_dims

# kA_5 : 1->2->3->4->5 over GF(5).  Interior interval [2,3] (right M) and [3,4].
_KA5 = {"kind": "quiver", "vertices": [1, 2, 3, 4, 5],
        "arrows": {"a1": [1, 2], "a2": [2, 3], "a3": [3, 4], "a4": [4, 5]},
        "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}}
# kA_4 with the relation a1*a2 = 0, for the relation-violation refusals.
_KA4_REL = {"kind": "quiver", "vertices": [1, 2, 3, 4],
            "arrows": {"a1": [1, 2], "a2": [2, 3], "a3": [3, 4]},
            "relations": ["a1*a2"], "field": {"kind": "GF", "p": 5, "n": 1}}

# Explicit interior modules as schema-v2 block matrices (target x source dims):
_M_23 = {"dims": {"2": 1, "3": 1}, "maps": {"a2": [[1]]}}   # right interval [2,3]
_N_34_right = {"dims": {"3": 1, "4": 1}, "maps": {"a3": [[1]]}}   # right interval [3,4]
_N_34_left = {"dims": {"3": 1, "4": 1}, "maps": {"a3": [[1]]}}    # left target (side->left)
# A module whose a1, a2 both act invertibly, violating a1*a2 = 0:
_BAD = {"dims": {"1": 1, "2": 1, "3": 1}, "maps": {"a1": [[1]], "a2": [[1]]}}


def _run(tmp_path, algebra, compute, module=None, ext_target=None, tor_target=None):
    body = {"schema": 2, "algebra": algebra, "compute": compute,
            "artifacts": {"pdf": False, "tikz": False}}
    if module is not None:
        body["module"] = module
    if ext_target is not None:
        body["ext_target"] = ext_target
    if tor_target is not None:
        body["tor_target"] = tor_target
    return run_spec(ComputeRequest.model_validate(body), tmp_path)


def _interval(A, i, j, n):
    dim = j - i + 1
    dv = {v: (1 if i <= v <= j else 0) for v in range(1, n + 1)}
    maps = {}
    for k in range(1, n):
        blk = [[0] * dim for _ in range(dim)]
        if i <= k <= j - 1:
            blk[(k + 1) - i][k - i] = 1
        maps[f"a{k}"] = blk
    return A.module(dv, maps)


# --------------------------------------------------------------------------- #
# Ext with two explicit non-builtin interior modules == library
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_crossengine
def test_ext_both_explicit_modules_match_library(tmp_path):
    b = _run(tmp_path, _KA5, ["ext:0..3"], _M_23, ext_target=_N_34_right)["results"]["ext"]
    A = linear_path_algebra(5, field=GF(5))
    M, N = _interval(A, 2, 3, 5), _interval(A, 3, 4, 5)
    assert b["dims"] == [A.ext(M, N, n) for n in range(4)] == [0, 1, 0, 0]
    assert b["target"]["dimvec"] == {"1": 0, "2": 0, "3": 1, "4": 1, "5": 0}
    assert b["references"] == ["module_ext"]


# --------------------------------------------------------------------------- #
# Tor with an explicit right module + an explicit LEFT tor_target == library
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_crossengine
def test_tor_explicit_left_target_matches_library(tmp_path):
    b = _run(tmp_path, _KA5, ["tor:0..3"], _M_23, tor_target=_N_34_left)["results"]["tor"]
    assert b["target"]["dimvec"] == {"1": 0, "2": 0, "3": 1, "4": 1, "5": 0}
    A = linear_path_algebra(5, field=GF(5))
    M = _interval(A, 2, 3, 5)
    # The left target: right A^op-module on {v3, v4} with a3 (A^op: 4->3) acting by 1.
    full = {"a1": [[0, 0], [0, 0]], "a2": [[0, 0], [0, 0]],
            "a3": [[0, 1], [0, 0]], "a4": [[0, 0], [0, 0]]}
    N_left = A.module({3: 1, 4: 1}, full, side="left")
    assert N_left.side == "left"
    assert b["dims"] == tor_dims(A, M, N_left, 3) == [0, 1, 0, 0]
    assert b["references"] == ["minimal_resolution", "module_ext"]


# --------------------------------------------------------------------------- #
# A relation-violating SECOND argument is a clean typed refusal (never a 500)
# --------------------------------------------------------------------------- #
def test_ext_target_relation_violation_is_typed_refusal(tmp_path):
    with pytest.raises(RunError) as exc:
        _run(tmp_path, _KA4_REL, ["ext:0..2"], {"dims": {"2": 1}}, ext_target=_BAD)
    assert is_safe_error_type(exc.value.error_type), exc.value.error_type
    assert "relation" in exc.value.message.lower()


def test_tor_target_relation_violation_is_typed_refusal(tmp_path):
    # The bad module is fed as the LEFT tor_target; the violation surfaces over A^op
    # (the relation a1*a2 becomes a2*a1) as the library's loud module error.
    with pytest.raises(RunError) as exc:
        _run(tmp_path, _KA4_REL, ["tor:0..2"], {"dims": {"2": 1}}, tor_target=_BAD)
    assert is_safe_error_type(exc.value.error_type), exc.value.error_type
    assert "relation" in exc.value.message.lower()


# --------------------------------------------------------------------------- #
# Side guards are SCHEMA-level typed refusals (FastAPI maps a request-body
# validation error to a 422 -- a clean 4xx, never a 500).  Ext^n(M, N) pairs a
# right M with a RIGHT N; Tor^A_n(M, N) a right M with a LEFT N.  A wrong-side
# second argument is refused up front at validation, before the runner is reached.
# --------------------------------------------------------------------------- #
def _validate(algebra, compute, module, **extra):
    body = {"schema": 2, "algebra": algebra, "compute": compute, "module": module,
            "artifacts": {"pdf": False, "tikz": False}, **extra}
    return ComputeRequest.model_validate(body)


def test_ext_target_left_side_is_schema_refusal():
    """A LEFT ``ext_target`` is refused UP FRONT at the schema, symmetric to the Tor
    guard -- Ext^n(M, N) pairs a right M with a RIGHT N (previously a left ext_target
    flowed through to the engine's own late refusal; now the schema stops it)."""
    left_ext = {"dims": {"3": 1, "4": 1}, "maps": {"a3": [[1]]}, "side": "left"}
    with pytest.raises(ValueError) as exc:
        _validate(_KA5, ["ext:0..3"], _M_23, ext_target=left_ext)
    assert "ext_target" in str(exc.value) and "RIGHT" in str(exc.value)


def test_tor_target_right_side_is_schema_refusal():
    """A RIGHT ``tor_target`` is refused at the schema (the existing Plan-30 guard) --
    Tor^A_n(M, N) pairs a right M with a LEFT N."""
    right_tor = {"dims": {"3": 1, "4": 1}, "maps": {"a3": [[1]]}, "side": "right"}
    with pytest.raises(ValueError) as exc:
        _validate(_KA5, ["tor:0..3"], _M_23, tor_target=right_tor)
    assert "tor_target" in str(exc.value) and "LEFT" in str(exc.value)
