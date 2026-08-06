"""Wave-2 no-code compute kinds: radical_filtration_ss, ar_quiver, derived_compare,
plus the quasi_hereditary / orbit_geometry block enrichments.

Cross-runner CONTRACT tests (byte-identical server / Pyodide twin), the honest-
refusal paths, estimator sizing, and canonical-key stability. Unmarked: these are
webapp-contract + infrastructure, not mathematical oracles (Plan-32 ruling -- twin
parity is not an oracle class, and the [web]/[hpc] tiers stay unmarked so the
oracle-class count audit is untouched)."""
import importlib.util
import json
import os
import pathlib
import tempfile

import pytest

from quiverlab.hpc import spec
from webapp.server.cache import canonical_key
from webapp.server.config import Config
from webapp.server.estimator import decide_tier, sizing_dim
from webapp.server.schema import ComputeRequest, SchemaError

_ROOT = pathlib.Path(__file__).resolve().parents[2]

_kA3 = {"kind": "quiver", "vertices": [1, 2, 3], "arrows": {"a1": [1, 2], "a2": [2, 3]},
        "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}}
_kA4 = {"kind": "quiver", "vertices": [1, 2, 3, 4],
        "arrows": {"a1": [1, 2], "a2": [2, 3], "a3": [3, 4]},
        "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}}
_kxx = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
        "relations": ["x*x"], "field": {"kind": "GF", "p": 5, "n": 1}}
_kronecker = {"kind": "quiver", "vertices": [1, 2],
              "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
              "field": {"kind": "GF", "p": 32003, "n": 1}}


def _server_block(req, kind):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"][kind]


def _twin_block(req, compute_key):
    path = os.path.join(_ROOT, "docs", "gui", "runner.py")
    s = importlib.util.spec_from_file_location("gui_twin_wave2_test", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    r = json.loads(mod.compute_one(compute_key))
    assert r["ok"], r
    return r["block"]


# --------------------------------------------------------------------------- #
# End-to-end + twin parity
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("req,compute_key,kind", [
    ({"schema": 1, "algebra": _kA3, "compute": ["radical_filtration_ss:0..4"]},
     "radical_filtration_ss:0..4", "radical_filtration_ss"),
    ({"schema": 1, "algebra": _kA3, "compute": ["ar_quiver"]}, "ar_quiver", "ar_quiver"),
    ({"schema": 1, "algebra": _kA3, "compute": ["ar_quiver:16"]}, "ar_quiver:16", "ar_quiver"),
    ({"schema": 1, "algebra": _kA3, "algebra_b": _kA4, "compute": ["derived_compare"]},
     "derived_compare", "derived_compare"),
])
def test_end_to_end_and_twin_parity(req, compute_key, kind):
    sb = _server_block(req, kind)
    tb = _twin_block(req, compute_key)
    assert sb.get("kind") == kind
    assert json.dumps(sb, sort_keys=True, default=str) == \
           json.dumps(tb, sort_keys=True, default=str)


def test_radical_filtration_ss_shape():
    b = _server_block({"schema": 1, "algebra": _kA3,
                       "compute": ["radical_filtration_ss:0..4"]}, "radical_filtration_ss")
    assert b["kind"] == "radical_filtration_ss"
    assert isinstance(b["abutment"], list) and isinstance(b["einf"], list)
    assert isinstance(b["collapse"], bool) and "prose" in b
    assert any(k == "Weibel1994" for k, _ in b["citations"])


def test_ar_quiver_complete_on_rep_finite():
    b = _server_block({"schema": 1, "algebra": _kA3, "compute": ["ar_quiver"]}, "ar_quiver")
    assert b["status"] == "complete" and b["complete"] is True
    assert b["num_vertices"] == 6 and len(b["arrows"]) == b["num_arrows"]
    assert all({"from", "to", "mult"} <= set(a) for a in b["arrows"])
    assert b["tau_orbits"] and any(k == "ARS1995" for k, _ in b["citations"])


def test_derived_compare_distinguishes():
    b = _server_block({"schema": 1, "algebra": _kA3, "algebra_b": _kA4,
                       "compute": ["derived_compare"]}, "derived_compare")
    assert b["kind"] == "derived_compare"
    assert b["distinguished_by"]                       # kA3 vs kA4 differ
    assert b["verdict"] == "distinguished"
    assert b["verdict_text"].startswith("distinguished by ")
    assert "fingerprint_a" in b and "fingerprint_b" in b


def test_derived_compare_not_distinguished_is_honest():
    # A vs itself: never claims equivalence, only "not distinguished by these invariants".
    b = _server_block({"schema": 1, "algebra": _kA3, "algebra_b": _kA3,
                       "compute": ["derived_compare"]}, "derived_compare")
    assert b["distinguished_by"] == []
    assert b["verdict"] == "not distinguished by these invariants"
    assert "equivalent" not in b["verdict_text"].lower()


# --------------------------------------------------------------------------- #
# Honest refusals
# --------------------------------------------------------------------------- #

def test_ar_quiver_self_injective_refusal():
    b = _server_block({"schema": 1, "algebra": _kxx, "compute": ["ar_quiver"]}, "ar_quiver")
    assert b["complete"] is False and b["status"] == "unsupported"
    assert "self-injective" in b["error"]
    assert b["vertices"] == [] and b["arrows"] == []


def test_ar_quiver_budget_cap_is_labelled_partial():
    # Kronecker is representation-infinite: a tiny budget cannot close the knit.
    b = _server_block({"schema": 1, "algebra": _kronecker, "compute": ["ar_quiver:4"]},
                      "ar_quiver")
    assert b["complete"] is False
    assert b["status"] in ("budget", "error", "unsupported")


def test_derived_compare_missing_algebra_b_refuses_spec():
    with pytest.raises(spec.SpecError):
        spec.parse_request({"schema": 1, "algebra": _kA3, "compute": ["derived_compare"]})


def test_derived_compare_missing_algebra_b_refuses_schema():
    # pydantic wraps the model-validator SchemaError (a ValueError) in a ValidationError
    # (itself a ValueError) -- the plan-30 schema-test convention.
    with pytest.raises(ValueError):
        ComputeRequest.model_validate(
            {"schema": 1, "algebra": _kA3, "compute": ["derived_compare"]})


# --------------------------------------------------------------------------- #
# Estimator sizing: oversized requests route off instant
# --------------------------------------------------------------------------- #

def test_ar_quiver_budget_is_not_a_degree(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = ComputeRequest.model_validate({"schema": 1, "algebra": _kA3,
                                         "compute": ["ar_quiver:512"]})
    # the budget 512 must NOT be treated as degree 512 (which would force queue);
    # a small algebra stays instant, a big one routes off (dim^3 sizing).
    assert decide_tier(4, req, cfg) == "instant"
    assert decide_tier(200, req, cfg) != "instant"


def test_derived_compare_sizes_on_the_larger_algebra(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    big_b = {"kind": "quiver", "vertices": list(range(1, 18)),
             "arrows": {f"a{i}": [i, i + 1] for i in range(1, 17)},
             "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}}  # kA17, dim 153
    req = ComputeRequest.model_validate({"schema": 1, "algebra": _kA3,
                                         "algebra_b": big_b, "compute": ["derived_compare"]})
    # A (kA3) has dim 6; B (kA17) has dim 153 -- the job must size on the larger.
    assert sizing_dim(6, req) == 153
    assert decide_tier(sizing_dim(6, req), req, cfg) != "instant"


def test_radical_filtration_ss_degree_sizes_off_instant(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = ComputeRequest.model_validate({"schema": 1, "algebra": _kA3,
                                         "compute": ["radical_filtration_ss:0..20"]})
    assert decide_tier(6, req, cfg) != "instant"       # degree 20 > instant_max_degree


# --------------------------------------------------------------------------- #
# Canonical-key stability: the algebra_b field never changes a request WITHOUT it
# --------------------------------------------------------------------------- #

_V = "0.1.0.dev0"


def test_absent_algebra_b_is_dropped_from_the_canonical_form():
    req = ComputeRequest.model_validate({"schema": 1, "algebra": _kA3,
                                         "compute": ["cartan"]})
    dumped = req.model_dump(by_alias=True)
    assert "algebra_b" not in dumped                   # absent -> serialized away


def test_canonical_key_unchanged_without_algebra_b():
    # The key of a plain request must be the SAME whether or not the schema knows
    # about algebra_b (it is dropped when None, so the pre-wave-2 bytes reproduce).
    body = {"schema": 1, "algebra": _kA3, "compute": ["cartan"]}
    req = ComputeRequest.model_validate(body)
    key_full = canonical_key(req.model_dump(by_alias=True), _V)
    dumped = req.model_dump(by_alias=True)
    dumped.pop("algebra_b", None)                       # what the pre-wave-2 shape had
    key_pre = canonical_key(dumped, _V)
    assert key_full == key_pre


def test_derived_compare_key_depends_on_algebra_b():
    # A genuine derived_compare request DOES carry algebra_b, so its key reflects B.
    r1 = ComputeRequest.model_validate({"schema": 1, "algebra": _kA3,
                                        "algebra_b": _kA4, "compute": ["derived_compare"]})
    r2 = ComputeRequest.model_validate({"schema": 1, "algebra": _kA3,
                                        "algebra_b": _kA3, "compute": ["derived_compare"]})
    assert canonical_key(r1.model_dump(by_alias=True), _V) != \
           canonical_key(r2.model_dump(by_alias=True), _V)


# --------------------------------------------------------------------------- #
# Block enrichments (items 4 + 5)
# --------------------------------------------------------------------------- #

_kA3_big = {**_kA3, "field": {"kind": "GF", "p": 32003, "n": 1}}


def test_quasi_hereditary_enriched_when_qh():
    # over char > dim the Ringel dual presents as kQ/I, so its Cartan matrix is available.
    b = _server_block({"schema": 1, "algebra": _kA3_big, "compute": ["quasi_hereditary"]},
                      "quasi_hereditary")
    assert b["is_quasi_hereditary"] is True
    ct = b["characteristic_tilting"]
    assert ct["dim"] == 6 and "summands" in ct
    rd = b["ringel_dual"]
    assert rd["dim"] == 6 and rd["cartan"]


def test_quasi_hereditary_ringel_dual_degrades_honestly_small_char():
    # over char 5 <= dim 6 the presented form is refused, so the Cartan is captured as
    # an honest per-field error -- the dimension still stands, the block never crashes.
    b = _server_block({"schema": 1, "algebra": _kA3, "compute": ["quasi_hereditary"]},
                      "quasi_hereditary")
    rd = b["ringel_dual"]
    assert rd["dim"] == 6
    assert "cartan" in rd or "cartan_error" in rd     # one or the other, never a crash


def test_quasi_hereditary_not_enriched_when_not_qh():
    b = _server_block({"schema": 1, "algebra": _kxx, "compute": ["quasi_hereditary"]},
                      "quasi_hereditary")
    assert b["is_quasi_hereditary"] is False
    assert "characteristic_tilting" not in b and "ringel_dual" not in b


def test_orbit_geometry_degeneration_absent_for_unique_module():
    # S_2 over kA3 is the UNIQUE module of dim (0,1,0): a trivial (one-class) order,
    # so no degeneration_order field is added (byte-identical to the pre-wave-2 block).
    req = {"schema": 2, "algebra": {**_kA3, "field": {"kind": "GF", "p": 32003, "n": 1}},
           "module": {"builtin": {"kind": "simple", "vertex": 2}},
           "compute": ["orbit_geometry"]}
    b = _server_block(req, "orbit_geometry")
    assert "degeneration_order" not in b


def test_orbit_geometry_degeneration_present_for_degenerable_module():
    # A decomposable S_1 (+) S_2 of dim (1,1,0) sits below the generic I_2 -- a genuine
    # two-class degeneration order, which the block now carries with its Hasse covers.
    req = {"schema": 2, "algebra": {**_kA3, "field": {"kind": "GF", "p": 32003, "n": 1}},
           "module": {"dims": {"1": 1, "2": 1, "3": 0},
                      "maps": {"a1": [[0]], "a2": []}, "side": "right"},
           "compute": ["orbit_geometry"]}
    b = _server_block(req, "orbit_geometry")
    deg = b["degeneration_order"]
    assert deg["complete"] is True and len(deg["vertices"]) == 2
    assert deg["covers"] == [[1, 0]]
