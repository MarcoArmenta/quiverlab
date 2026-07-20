"""The QPA GAP-source builders are PURE string generation from an Algebra's
public presentation -- no GAP needed to exercise them. These `fast` tests run
everywhere (they escape the `qpa` bucket via the explicit marker) and pin the
Python-side correctness that the GAP-only path would otherwise only hit in CI:
reading `algebra.quiver`/`algebra.relations` (public attrs, NOT `_quiver`/
`_relations`) and `domain.characteristic` (an int attribute, NOT a method)."""
import pytest

from quiverlab import Quiver, GF
from quiverlab.qpa.scripts import quiver_and_algebra_script

pytestmark = pytest.mark.fast


def test_ka2_gf2_script_no_relations():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=GF(2))
    src = quiver_and_algebra_script(A)
    assert 'Quiver(2, [[1, 2, "a"]])' in src
    assert "PathAlgebra(GF(2), Q)" in src
    assert "A := kQ;;" in src            # no relations -> A is the free path algebra


def test_gfp_monomial_script_includes_relations():
    # a relational presentation must reach the `rels` branch -- the `_relations`
    # bug silently dropped it (getattr(algebra, "_relations", None) is always None).
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x*x"], field=GF(3))
    src = quiver_and_algebra_script(A)
    assert "PathAlgebra(GF(3), Q)" in src
    assert "rels :=" in src and "A := kQ/rels;;" in src


def test_builder_reads_public_attrs_not_private():
    # regression guard: a swap back to algebra._quiver / _relations / .characteristic()
    # would AttributeError/TypeError here, since Algebra has no _quiver/_relations attr
    # and Domain.characteristic is an int, not a callable.
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=GF(5))
    quiver_and_algebra_script(A)  # must not raise
