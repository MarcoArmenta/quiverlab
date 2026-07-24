import pytest

from webapp.server.schema import ComputeRequest, parse_compute_item, SchemaError


def test_parses_family_request():
    req = ComputeRequest.model_validate({
        "schema": 1,
        "algebra": {"kind": "family", "family": "QuantumCI",
                    "params": {"n": 3}, "field": {"kind": "GF", "p": 5, "n": 1}},
        "compute": ["hh_cohomology:0..6", "coxeter_polynomial"],
        "artifacts": {"pdf": True, "tikz": False},
    })
    assert req.algebra.kind == "family"
    assert req.algebra.family == "QuantumCI"
    assert req.algebra.field.p == 5


def test_parses_quiver_request():
    # The Plan-10 GUI already emits quiver requests; schema v1 must accept them
    # (amendment point 7). Relation strings are opaque data at this layer.
    req = ComputeRequest.model_validate({
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2],
                    "arrows": {"a": [1, 2]}, "relations": [],
                    "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": ["cartan"],
        "artifacts": {"pdf": False, "tikz": False},
    })
    assert req.algebra.kind == "quiver"
    assert req.algebra.vertices == [1, 2]
    assert req.algebra.arrows == {"a": (1, 2)}
    assert req.algebra.relations == []


def test_rejects_unknown_algebra_kind():
    with pytest.raises(ValueError):
        ComputeRequest.model_validate({
            "schema": 1,
            "algebra": {"kind": "matrix", "family": "x", "params": {},
                        "field": {"kind": "CC"}},
            "compute": ["cartan"],
            "artifacts": {"pdf": False, "tikz": False},
        })


def test_compute_item_range():
    item = parse_compute_item("hh_cohomology:0..6")
    assert (item.kind, item.lo, item.hi) == ("hh_cohomology", 0, 6)


def test_compute_item_scalar():
    item = parse_compute_item("coxeter_polynomial")
    assert (item.kind, item.lo, item.hi) == ("coxeter_polynomial", None, None)


def test_compute_item_bad_range():
    with pytest.raises(SchemaError):
        parse_compute_item("hh_cohomology:6..0")


def test_anonymous_request_has_no_email_field():
    # Email belongs ONLY to the big-job endpoint (Task 13); the anonymous
    # compute request must never carry it (pydantic drops the unknown key).
    req = ComputeRequest.model_validate({
        "schema": 1,
        "algebra": {"kind": "family", "family": "QuantumCI",
                    "params": {"n": 3}, "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False},
        "email": "sneaky@example.org"})
    assert not hasattr(req, "email")
