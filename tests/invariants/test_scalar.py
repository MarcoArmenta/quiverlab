"""loewy_length / complexity / center (spec section 3.5). Fixtures A, B; k[x]/(x^n)."""
import pytest
from quiverlab import Quiver, CC, GF, linear_path_algebra, truncated_polynomial
from quiverlab.errors import FieldError


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"])


def test_loewy_length():
    assert linear_path_algebra(2).loewy_length() == 2        # rad^2 = 0
    assert _square().loewy_length() == 3                     # rad^3 = 0, rad^2 = <a*b>
    assert truncated_polynomial(4).loewy_length() == 4       # k[x]/(x^4): rad^4 = 0
    assert truncated_polynomial(2, field=GF(3)).loewy_length() == 2


def test_center_dimension_and_basis():
    d, basis = linear_path_algebra(2).center()
    assert d == 1                                            # connected -> center = k*1
    assert len(basis) == 1
    assert _square().center()[0] == 1
    # k[x]/(x^2) is commutative -> center is the whole algebra
    assert truncated_polynomial(2).center()[0] == 2


def test_center_over_gfp():
    assert linear_path_algebra(2, field=GF(7)).center()[0] == 1


def test_complexity_gfp():
    # finite gl.dim (kA_2) -> minimal A^e resolution terminates -> complexity 0
    assert linear_path_algebra(2, field=GF(32003)).complexity(6) == 0
    # k[x]/(x^2) self-injective -> constant-rank periodic resolution -> complexity 1
    assert truncated_polynomial(2, field=GF(32003)).complexity(6) == 1


def test_complexity_cc_computes():
    # Plan 19: off GF(p) complexity runs on the relative-Tor Betti complex
    assert truncated_polynomial(2, field=CC).complexity(4) == 1


def test_complexity_refuses_a_truncated_resolution(monkeypatch):
    # Honesty guard: a memory/term-dimension-truncated minimal A^e resolution would
    # under-report complexity (an early stop reads as a genuine termination). The build's
    # truncation marker (4th return value) is now consulted -- a non-None marker raises
    # QuiverlabError instead of certifying a wrong number. The default caps never
    # truncate in normal use, so this monkeypatches the engine to simulate a truncated
    # build (cheap + deterministic) and asserts the loud refusal.
    from quiverlab.errors import QuiverlabError
    import quiverlab.engine.resolutions_minimal as rm
    A = truncated_polynomial(2, field=GF(32003))

    real = rm.minimal_resolution

    def _truncated(*args, **kwargs):
        rks, cols, eng, _ = real(*args, **kwargs)
        return rks, cols, eng, 1                 # pretend the build stopped at degree 1
    monkeypatch.setattr(rm, "minimal_resolution", _truncated)
    with pytest.raises(QuiverlabError, match="truncated"):
        A.complexity(6)


def test_complexity_untruncated_value_is_byte_identical(monkeypatch):
    # The guard must NOT touch the normal (untruncated) return: with the marker forced
    # to None (the untruncated case) the value equals the shipped complexity.
    import quiverlab.engine.resolutions_minimal as rm
    A = truncated_polynomial(2, field=GF(32003))
    real = rm.minimal_resolution

    def _untruncated(*args, **kwargs):
        rks, cols, eng, _ = real(*args, **kwargs)
        return rks, cols, eng, None
    monkeypatch.setattr(rm, "minimal_resolution", _untruncated)
    assert A.complexity(6) == 1                   # k[x]/x^2 self-injective -> complexity 1
