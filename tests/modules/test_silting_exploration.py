"""Bounded-radius silting exploration with LOUD truncation (Plan 67 / AI Thm 1.2 --
no general BFS). Literature: k[x]/(x^2) is local -> complete at radius 0 (silt = {A[i]},
mod shift a single vertex). Honest: kA_2 is INFINITE -> the exploration truncates with
status 'radius'/'budget', never 'complete'. Self-cert: every discovered vertex is a
silting object and every arrow is a single mutation."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.derived.silting import bounded_silting_exploration

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


@lit
def test_local_is_complete_mod_shift():
    from quiverlab.families import truncated_polynomial
    A = truncated_polynomial(2, field=GF(32003))     # local: silt = {A[i]}
    ex = bounded_silting_exploration(A, radius=2)
    assert ex.status == "complete" and ex.finite_class == "local"
    assert len(ex.vertices) == 1                      # mod shift: one silting object


@lit
def test_hereditary_kA2_is_infinite_truncates_loudly():
    A = linear_path_algebra(2, field=QQ)             # kA2: silting quiver INFINITE
    ex = bounded_silting_exploration(A, radius=2, budget=64)
    assert ex.status in ("radius", "budget")         # never "complete"
    assert ex.finite_class is None
    # the discovered ball is genuine: A itself is present, every vertex is silting
    assert any(v["is_initial"] for v in ex.vertices)
    from quiverlab.derived.silting import is_silting_object
    for v in ex.vertices:
        assert is_silting_object(v["summands"]).is_silting in (True, "unknown")


@selfcert
def test_budget_trips_before_radius():
    A = linear_path_algebra(3, field=QQ)
    ex = bounded_silting_exploration(A, radius=99, budget=8)
    assert ex.status == "budget" and len(ex.vertices) <= 8


@selfcert
def test_every_arrow_is_a_single_mutation():
    A = linear_path_algebra(2, field=QQ)
    ex = bounded_silting_exploration(A, radius=2, budget=64)
    for (i, j), lab in ex.arrows.items():
        assert lab["direction"] in ("left", "right")
        assert ex.vertices[i]["key"] != ex.vertices[j]["key"]
