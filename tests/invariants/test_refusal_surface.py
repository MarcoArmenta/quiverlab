"""Plan 19 refusal contract: off GF(p), only path-basis-needing invariants
on quiver-less algebras refuse — with an honest message (no 'later phase'
promise anywhere in src/). cyclic_homology refuses NOWHERE."""
import pathlib

import pytest

import quiverlab
from quiverlab import CC, truncated_polynomial
from quiverlab.core.algebra import Algebra
from quiverlab.errors import FieldError


def _dual_numbers_sc(field=CC):
    """k[x]/(x^2) via raw structure constants: NO quiver, NO e_ labels."""
    T = [[[1, 0], [0, 1]],
         [[0, 1], [0, 0]]]
    return Algebra.from_structure_constants(T, [1, 0], field=field)


def test_structure_constants_refusals_are_honest():
    A = _dual_numbers_sc(CC)
    for call in (A.is_frobenius, A.is_symmetric, A.nakayama_automorphism,
                 lambda: A.complexity(3)):
        with pytest.raises(FieldError) as ei:
            call()
        msg = str(ei.value)
        assert "later phase" not in msg
        assert "quiver" in msg or "path-type" in msg


def test_cyclic_homology_needs_no_quiver():
    # the bar mixed complex serves ANY unital algebra: same dims from the
    # structure-constants presentation and the quiver presentation
    got = _dual_numbers_sc(CC).cyclic_homology(2)
    want = truncated_polynomial(2, field=CC).cyclic_homology(2)
    assert got.dims == want.dims


def test_no_later_phase_promise_left_in_src():
    root = pathlib.Path(quiverlab.__file__).parent
    hits = [str(p) for p in root.rglob("*.py")
            if "later phase that generalizes" in p.read_text(encoding="utf-8")]
    assert hits == []
