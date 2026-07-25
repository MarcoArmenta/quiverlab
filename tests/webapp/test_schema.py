import pytest

from webapp.server.schema import (
    ComputeRequest, ModuleSpec, parse_compute_item, SchemaError,
)


_QUIVER = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
           "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}}


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


# --------------------------------------------------------------------------- #
# Schema v2: the no-code module block (Plan 26)
# --------------------------------------------------------------------------- #

def _module_req(**over):
    body = {"schema": 2, "algebra": _QUIVER, "compute": ["dimension_vector"],
            "artifacts": {"pdf": False, "tikz": False},
            "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}}}
    body.update(over)
    return ComputeRequest.model_validate(body)


def test_parses_explicit_module_request():
    req = _module_req(compute=["dimension_vector", "rad_top_soc", "tau", "tau_minus"])
    assert req.schema_version == 2
    assert req.module.dims == {"1": 2}
    assert req.module.maps == {"x": [[0, 0], [1, 0]]}
    assert req.module.side == "right"          # default
    assert req.module.builtin is None


def test_parses_builtin_pick_list_with_nested_side():
    # The task's builtin shape carries `side` INSIDE the builtin dict; it is lifted
    # to the top-level canonical `side`.
    req = _module_req(module={"builtin": {"kind": "projective", "vertex": 1,
                                          "side": "left"}})
    assert req.module.builtin.kind == "projective"
    assert req.module.builtin.vertex == 1
    assert req.module.side == "left"           # lifted out of builtin


def test_side_omitted_and_explicit_right_canonicalize_identically():
    # The load-bearing cache invariant (Plan 25): an omitted side and an explicit
    # "right" must dump identically, so the same module typed twice collides.
    omitted = _module_req().model_dump(by_alias=True)
    explicit = _module_req(module={"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]},
                                   "side": "right"}).model_dump(by_alias=True)
    assert omitted == explicit
    assert omitted["module"]["side"] == "right"


def test_non_module_request_dumps_without_module_key():
    # A family/quiver request with no module block must dump byte-identically to the
    # pre-Plan-26 shape (no `module`/`ext_target` keys) so its cache key is unchanged.
    req = ComputeRequest.model_validate(
        {"schema": 1, "algebra": _QUIVER, "compute": ["cartan"]})
    d = req.model_dump(by_alias=True)
    assert "module" not in d and "ext_target" not in d
    assert set(d) == {"schema", "algebra", "compute", "artifacts"}


def test_module_block_requires_schema_2():
    with pytest.raises(ValueError):
        ComputeRequest.model_validate(
            {"schema": 1, "algebra": _QUIVER, "compute": ["tau"],
             "module": {"dims": {"1": 1}, "maps": {"x": [[0]]}}})


def test_module_kind_requires_a_module_block():
    with pytest.raises(ValueError):
        ComputeRequest.model_validate(
            {"schema": 2, "algebra": _QUIVER, "compute": ["tau"]})


def test_ext_requires_ext_target():
    with pytest.raises(ValueError):
        ComputeRequest.model_validate(
            {"schema": 2, "algebra": _QUIVER, "compute": ["ext:0..3"],
             "module": {"dims": {"1": 1}, "maps": {"x": [[0]]}}})


def test_ext_with_two_modules_parses():
    req = ComputeRequest.model_validate(
        {"schema": 2, "algebra": _QUIVER, "compute": ["ext:0..3"],
         "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}},
         "ext_target": {"builtin": {"kind": "simple", "vertex": 1}}})
    assert req.ext_target.builtin.kind == "simple"
    assert "ext_target" in req.model_dump(by_alias=True)


def test_float_matrix_entry_is_rejected():
    with pytest.raises(ValueError):
        ModuleSpec.model_validate({"dims": {"1": 1}, "maps": {"x": [[1.5]]}})


def test_bool_matrix_entry_is_rejected():
    with pytest.raises(ValueError):
        ModuleSpec.model_validate({"dims": {"1": 1}, "maps": {"x": [[True]]}})


def test_exact_string_entries_are_accepted():
    m = ModuleSpec.model_validate({"dims": {"1": 2}, "maps": {"x": [["1/2", "-3"], ["0", "7"]]}})
    assert m.maps["x"] == [["1/2", "-3"], ["0", "7"]]


def test_non_rectangular_matrix_is_rejected():
    with pytest.raises(ValueError):
        ModuleSpec.model_validate({"dims": {"1": 2}, "maps": {"x": [[1, 0], [0]]}})


def test_module_cannot_be_both_explicit_and_builtin():
    with pytest.raises(ValueError):
        ModuleSpec.model_validate(
            {"dims": {"1": 1}, "maps": {}, "builtin": {"kind": "simple", "vertex": 1}})


def test_builtin_conflicting_side_is_rejected():
    with pytest.raises(ValueError):
        ModuleSpec.model_validate(
            {"builtin": {"kind": "simple", "vertex": 1, "side": "left"}, "side": "right"})


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
