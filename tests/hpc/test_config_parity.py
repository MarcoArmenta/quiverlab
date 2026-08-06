"""Config-parity: a GUI-shaped request validates IDENTICALLY (accept/reject) under
the stdlib spec core and under the webapp's pydantic ``ComputeRequest``. Importing
the webapp schema inside the test is fine (the import-boundary contract is about
``import quiverlab.hpc``, not about the tests)."""
import pytest

from quiverlab.hpc.spec import parse_request, SpecError

# The QPA CI job installs without the [web] extra; collection must skip cleanly
# there (same guard pattern as tests/webapp/conftest.py).
pytest.importorskip("pydantic", reason="config-parity test needs the [web] extra")

from webapp.server.schema import ComputeRequest, SchemaError as WSchemaError
from pydantic import ValidationError

_GF2 = {"kind": "GF", "p": 2, "n": 1}

# GUI-shaped valid requests (schema 1 quiver, schema 1 family, schema 2 module).
_VALID = [
    {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                              "relations": ["x*x*x"], "field": _GF2},
     "compute": ["hh_cohomology:0..4", "cartan"], "artifacts": {"pdf": False, "tikz": False}},
    {"schema": 1, "algebra": {"kind": "family", "family": "QuantumCI", "params": {"q": 1},
                              "field": _GF2},
     "compute": ["hh_homology:0..3"], "artifacts": {"pdf": False, "tikz": False}},
    {"schema": 2, "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                              "relations": ["x*x*x"], "field": _GF2},
     "compute": ["dimension_vector"], "artifacts": {"pdf": False, "tikz": False},
     "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}}},
    # QQ (the rational field) is now an accepted exact input field alongside CC/GF.
    {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                              "relations": [], "field": {"kind": "QQ"}},
     "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}},
]

# Requests both validators must REJECT.
_INVALID = [
    # empty compute
    {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1], "arrows": {},
                              "relations": [], "field": _GF2},
     "compute": [], "artifacts": {}},
    # empty vertices
    {"schema": 1, "algebra": {"kind": "quiver", "vertices": [], "arrows": {},
                              "relations": [], "field": _GF2},
     "compute": ["cartan"], "artifacts": {}},
    # unknown field kind (QQ is now valid; RR -- the reals -- is not an exact field)
    {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1], "arrows": {},
                              "relations": [], "field": {"kind": "RR"}},
     "compute": ["cartan"], "artifacts": {}},
    # module block on schema 1
    {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                              "relations": ["x*x*x"], "field": _GF2},
     "compute": ["cartan"], "artifacts": {},
     "module": {"dims": {"1": 1}, "maps": {}}},
    # module compute kind with no module block
    {"schema": 2, "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                              "relations": ["x*x*x"], "field": _GF2},
     "compute": ["dimension_vector"], "artifacts": {}},
    # unparseable compute item
    {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1], "arrows": {},
                              "relations": [], "field": _GF2},
     "compute": ["not a kind!"], "artifacts": {}},
]


@pytest.mark.parametrize("body", _VALID)
def test_valid_requests_accepted_by_both(body):
    req = parse_request(body)
    model = ComputeRequest.model_validate(body)
    assert req.schema_version == model.schema_version
    assert req.algebra.kind == model.algebra.kind
    assert list(req.compute) == list(model.compute)


@pytest.mark.parametrize("body", _INVALID)
def test_invalid_requests_rejected_by_both(body):
    with pytest.raises(SpecError):
        parse_request(body)
    with pytest.raises((ValidationError, WSchemaError)):
        ComputeRequest.model_validate(body)
