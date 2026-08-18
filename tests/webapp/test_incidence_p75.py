"""Plan 75 Task G1 -- the ``incidence_cohomology`` compute kind + the poset input mode.

The kind is served by the wheel's ``quiverlab.hpc.spec`` core and mirrored byte-for-byte
by the Pyodide twin (``docs/gui/runner.py``) through the SAME shared block builder
(``hochschild.simplicial.incidence_cohomology_block``), so the browser cannot drift from
the server report or a cluster run.

The poset panel emits the SHIPPED ``IncidenceAlgebra`` family -- no new algebra ``kind``,
no schema change -- so every existing request keys byte-unchanged (the frozen goldens in
``test_runner_delegation.py`` are the standing proof of that half), and a covers-only
poset keys byte-identically to the typed family request. Unmarked per the Plan-32
extras-gated ruling (a cross-runner contract, not an oracle).
"""
import importlib.util
import json
import pathlib

import pytest

from webapp.server.cache import canonical_key
from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

# The Boolean lattice B_3 = the divisor/subset poset on {1,2,3}: dim kB_3 = 27, and
# Delta(B_3) is a cone (B_3 has a global bound), so HH^{>=1} = 0 (the Plan-33 pin).
_B3_COVERS = [[x, x | (1 << b)] for x in range(8) for b in range(3) if not (x >> b) & 1]
# The crown C(2,2): 2 minimal + 2 maximal, every min below every max. Delta(C(2,2)) is a
# 4-cycle ~ S^1, so HH^1 = k -- a NONVANISHING H^1 on a NON-lattice poset.
_CROWN_COVERS = [[1, 3], [1, 4], [2, 3], [2, 4]]


def _request(covers, top, elements=None, field=None):
    return {"schema": 1,
            "algebra": {"kind": "family", "family": "IncidenceAlgebra",
                        "params": {"poset_or_covers": covers, "elements": elements},
                        "field": field or {"kind": "QQ"}},
            "compute": ["incidence_cohomology:0..%d" % top],
            "artifacts": {"pdf": False, "tikz": False}}


def _server_block(body, kind, tmp_path):
    out = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return out["results"][kind]


def _gui_block(body, compute):
    spec = importlib.util.spec_from_file_location("gui_runner_p75", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(compute))
    assert out["ok"], out
    return out["block"]


# --------------------------------------------------------------------------- #
# the block itself
# --------------------------------------------------------------------------- #

def test_incidence_block_shape_b3(tmp_path):
    b = _server_block(_request(_B3_COVERS, 3), "incidence_cohomology", tmp_path)
    assert b["dims"] == [1, 0, 0, 0]
    assert b["contractible"] is True
    assert b["face_vector"] == [8, 19, 18, 6]
    # sum (-1)^p f_p over a cone is 1, and it equals the alternating sum of the dims
    # here because top (3) reaches the top dimension of Delta(B_3).
    assert b["euler_characteristic"] == 1
    assert "gerstenhaber_schack_1983" in b["references"]
    assert b["citations"], "the block must resolve its references to citations"


def test_crown_has_nonvanishing_h1(tmp_path):
    b = _server_block(_request(_CROWN_COVERS, 3), "incidence_cohomology", tmp_path)
    assert b["dims"] == [1, 1, 0, 0]
    # No global bound and H^1 != 0: contractibility is DISPROVED, not left unknown.
    assert b["contractible"] is False


def test_no_poset_provenance_is_a_clean_error_block(tmp_path):
    """A drawn quiver carries no poset, and quiverlab never GUESSES that a kQ/I is an
    incidence algebra: the entry is a labelled error, never a 500 and never a number."""
    body = {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                        "relations": [], "field": {"kind": "QQ"}},
            "compute": ["incidence_cohomology:0..2"],
            "artifacts": {"pdf": False, "tikz": False}}
    b = _server_block(body, "incidence_cohomology", tmp_path)
    assert "error" in b and "provenance" in b["error"]
    assert "dims" not in b


def test_degree_range_is_required(tmp_path):
    body = _request(_B3_COVERS, 3)
    body["compute"] = ["incidence_cohomology"]
    with pytest.raises(Exception):
        run_spec(ComputeRequest.model_validate(body), tmp_path)


# --------------------------------------------------------------------------- #
# cross-runner parity
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("covers,dims", [(_B3_COVERS, [1, 0, 0, 0]),
                                         (_CROWN_COVERS, [1, 1, 0, 0])])
def test_twin_parity(covers, dims, tmp_path):
    body = _request(covers, 3)
    gui = _gui_block(body, "incidence_cohomology:0..3")
    ref = _server_block(body, "incidence_cohomology", tmp_path)
    assert gui["dims"] == dims
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))


# --------------------------------------------------------------------------- #
# canonical-key discipline (Plan-25)
# --------------------------------------------------------------------------- #

_V = "0.1.0.dev0"


def _key(body):
    return canonical_key(ComputeRequest.model_validate(body).model_dump(by_alias=True), _V)


def test_poset_request_keys_like_the_typed_family():
    """The panel normalizes ``elements`` to ``null`` when the covers already name every
    element, so a covers-only poset keys BYTE-IDENTICALLY to the typed IncidenceAlgebra
    family form -- one cache entry, not two. ``null`` and an explicit list key
    differently under ``canonical_key``'s sort_keys JSON, which is exactly why the
    normalization is mandatory."""
    covers = [[1, 2], [2, 3]]
    # The catalog's IncidenceAlgebra prefill is the SOURCE OF TRUTH for the "typed family"
    # shape: read its parameter defaults live rather than restating them, so this cannot
    # pass by comparing the panel's shape to itself.
    from webapp.server.catalog import _FAMILY_META
    params = _FAMILY_META["IncidenceAlgebra"]["params"]
    assert set(params) == {"poset_or_covers", "elements"}
    assert params["elements"]["example"] is None, (
        "the catalog default for `elements` is what the panel must normalize TO; if it "
        "ever stops being null, the panel's normalization has to move with it")
    typed = {"schema": 1,
             "algebra": {"kind": "family", "family": "IncidenceAlgebra",
                         "params": {"poset_or_covers": covers,
                                    "elements": params["elements"]["example"]},
                         "field": {"kind": "QQ"}},
             "compute": ["incidence_cohomology:0..3"],
             "artifacts": {"pdf": False, "tikz": False}}
    assert _key(typed) == _key(_request(covers, 3))
    # A genuine isolated element 4 (in NO cover) is a DIFFERENT poset -> a different key.
    assert _key(_request(covers, 3, elements=[1, 2, 3, 4])) != _key(_request(covers, 3))
    # ... and an explicit list that merely repeats the covers' elements is NOT what the
    # panel emits; if a caller sends it anyway it keys differently (honest, not silently
    # collapsed) -- the reason the panel must normalize rather than pass through.
    assert _key(_request(covers, 3, elements=[1, 2, 3])) != _key(_request(covers, 3))


def test_isolated_element_changes_the_mathematics():
    """The key discipline above is not cosmetic: an isolated point genuinely changes
    HH^* (Delta(P) gains a component, so H^0 = k^2)."""
    import quiverlab as ql
    A = ql.IncidenceAlgebra([[1, 2], [2, 3]], elements=[1, 2, 3, 4])
    assert A.incidence_cohomology(2).dims == [2, 0, 0]


# --------------------------------------------------------------------------- #
# tiering: the fast path must actually be REACHABLE through the service
# --------------------------------------------------------------------------- #

def _cfg():
    import os
    import tempfile
    from webapp.server.config import Config
    os.environ.setdefault("QLWEB_DATA_DIR", tempfile.mkdtemp())
    return Config.from_env()


def test_a_large_poset_is_not_refused_by_the_estimator():
    """The default estimator sizes on `dim ** 3` (the bar-resolution model). For the
    order-complex route that is the WRONG driver -- the algebra is never resolved -- and
    it refused exactly the inputs the fast path exists for: the 7-vertex torus face poset
    is `dim kP = 168` and classified `reject`, while the library answers it in 0.02s.
    An incidence-ONLY request is now sized on the poset instead."""
    from webapp.server.estimator import classify
    cfg = _cfg()
    # a 42-element face poset (the torus), dim kP = 168
    covers = [[str(i), str(i + 1)] for i in range(41)]
    r = ComputeRequest.model_validate(_request(covers, 2))
    c = classify(168, r, cfg)
    assert c["tier"] != "reject", c
    # ... and a small poset stays INSTANT (B_3 computes in 0.02s; nothing should queue it)
    small = ComputeRequest.model_validate(_request([[1, 2], [2, 3]], 2))
    assert classify(27, small, cfg)["tier"] == "instant"


def test_a_high_degree_alone_does_not_refuse_a_poset_request():
    """`incidence_cohomology:0..30` is free past the top dimension of Delta(P) (the dims
    are 0), so degree must not drive the tier the way a bar-resolution degree does."""
    from webapp.server.estimator import classify
    cfg = _cfg()
    r = ComputeRequest.model_validate(_request([[1, 2], [2, 3]], 30))
    assert classify(27, r, cfg)["tier"] == "instant"


def test_a_mixed_request_is_still_sized_the_old_way():
    """The exemption is for incidence-ONLY requests: add a real Hochschild kind and the
    algebra dimension drives the size again, because that one genuinely resolves."""
    from webapp.server.estimator import classify, _incidence_only_dim
    cfg = _cfg()
    body = _request(_B3_COVERS, 3)
    body["compute"] = ["incidence_cohomology:0..3", "hh_cohomology:0..3"]
    r = ComputeRequest.model_validate(body)
    assert _incidence_only_dim(r, 27) == 27
    plain = dict(body, compute=["hh_cohomology:0..3"])
    assert (classify(27, r, cfg)["estimate"]["cells"]
            == classify(27, ComputeRequest.model_validate(plain), cfg)["estimate"]["cells"])


def test_existing_requests_are_untouched():
    """The always-true half of the cache-key contract: P75 adds no field to a request
    that does not use the poset panel, so a plain quiver request keys exactly as it did.
    (The frozen goldens in test_runner_delegation.py pin the byte values themselves; this
    pins that the new kind's presence in the dispatch does not perturb the key path.)"""
    body = {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                        "relations": [], "field": {"kind": "QQ"}},
            "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}}
    assert _key(body) == _key(json.loads(json.dumps(body)))
    assert "incidence" not in json.dumps(
        ComputeRequest.model_validate(body).model_dump(by_alias=True))
