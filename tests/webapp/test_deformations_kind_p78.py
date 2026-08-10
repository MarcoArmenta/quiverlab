"""Plan 78 Task 5: the `deformations` compute kind served by the webapp/hpc tier --
grammar (bare + budget suffix + refusal), the block contract, the HH-richness estimator
sizing, and canonical-key stability. Unmarked per the Plan-32 extras-gated ruling (a
cross-runner/server contract, not an oracle)."""
import pytest

from webapp.server.estimator import sizing_dim
from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest, SchemaError, parse_compute_item


def _quiver(vertices, arrows, relations, field, compute):
    return {"schema": 2,
            "algebra": {"kind": "quiver", "vertices": vertices, "arrows": arrows,
                        "relations": relations, "field": field},
            "compute": compute, "artifacts": {"pdf": False, "tikz": False}}


# rad^2 = 0 on a 3-cycle: the CHEAP side of the cost law (HH^2 thin) -- 0.03 s.
_RADSQ = _quiver([1, 2, 3], {"a": [1, 2], "b": [2, 3], "c": [3, 1]},
                 ["a*b", "b*c", "c*a"], {"kind": "QQ"}, ["deformations"])


# --------------------------------------------------------------- grammar
def test_grammar_bare_and_budget_suffix():
    bare = parse_compute_item("deformations")
    assert bare.kind == "deformations" and bare.lo is None and bare.hi is None
    sized = parse_compute_item("deformations:32")
    assert sized.kind == "deformations" and sized.hi == 32
    # the budget is a DIM cap, not a degree range -- the 'name:0..N' grammar must NOT apply
    with pytest.raises(SchemaError, match="positive integer"):
        parse_compute_item("deformations:0..2")


# --------------------------------------------------------------- the block
def test_radsq_block_is_dg_lie(tmp_path):
    block = run_spec(ComputeRequest.model_validate(_RADSQ), tmp_path)["results"]["deformations"]
    assert block["kind"] == "deformations"
    # rad^2 = 0 => the B(A)[1] L-infinity companion is UNCONDITIONALLY dg-Lie: Plan 78's
    # GUARANTEED deliverable (no l_>=3 is ever shipped -- the spike stays frozen).
    assert block["dg_lie"] is True
    assert block["hh2_dim"] is not None and block["hh3_dim"] is not None
    assert block["citations"], "the block must carry resolved citations, not just refs"


# ------------------------------------------------- estimator: HH-richness, not dim
def test_sizing_is_hh_richness_not_algebra_dim():
    """The whole point of `_deformations_dim` (Plan 78 / H1): a SMALL algebra with rich
    HH^2 must size LARGER than a BIGGER algebra with HH^2 = 0, because the obstruction
    bracket -- not the dimension -- is the cost driver. Sizing on A.dim is backwards."""
    rich = ComputeRequest.model_validate(
        _quiver([1], {"x": [1, 1], "y": [1, 1]}, ["x*x", "y*y", "x*y - y*x"],
                {"kind": "QQ"}, ["deformations"]))          # dim 4, HH^2 = 5
    thin = ComputeRequest.model_validate(
        _quiver([1, 2, 3], {"a": [1, 2], "b": [2, 3], "c": [3, 1]},
                ["a*b", "b*c", "c*a"], {"kind": "QQ"}, ["deformations"]))   # dim 6, HH^2 thin
    assert sizing_dim(4, rich) > sizing_dim(6, thin)


def test_sizing_untouched_without_a_deformations_item():
    """Every request that does NOT ask for `deformations` must classify EXACTLY as before
    (the helper returns 0 and never builds/probes anything)."""
    other = ComputeRequest.model_validate(
        _quiver([1], {"x": [1, 1], "y": [1, 1]}, ["x*x", "y*y", "x*y - y*x"],
                {"kind": "QQ"}, ["hh_cohomology:0..2"]))
    assert sizing_dim(4, other) == 4


def test_budget_is_not_a_degree_for_the_tier_estimate():
    """`hi` is a DIM budget; it must not be read as a homological degree by the estimator
    (the hh1_lie / split_extension precedent -- it joins the skip tuple)."""
    from webapp.server.estimator import _max_degree
    req = ComputeRequest.model_validate(
        _quiver([1, 2, 3], {"a": [1, 2], "b": [2, 3], "c": [3, 1]},
                ["a*b", "b*c", "c*a"], {"kind": "QQ"}, ["deformations:32"]))
    assert _max_degree(req) < 32
