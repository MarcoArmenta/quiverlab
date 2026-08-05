"""Derived fingerprint: a necessary-condition invariant tuple + honest comparison.
Literature: D4 vs A4 distinguished by the Coxeter polynomial; the 8-vertex
cospectral trees NOT distinguished (equal fingerprint) -- honest verdict."""
import pytest

from quiverlab import Quiver
from quiverlab.derived.fingerprint import derived_fingerprint, compare_fingerprints

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _d4():
    return Quiver([1, 2, 3, 4], {"a": (2, 1), "b": (3, 1), "c": (4, 1)}).algebra()


def _a4():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra()


def _spider():
    return Quiver(list(range(1, 9)),
                  {"a1": (1, 2), "a2": (2, 3), "a3": (3, 4),
                   "b": (1, 5), "c": (1, 6), "d": (1, 7), "e": (1, 8)}).algebra()


def _double_star():
    return Quiver(list(range(1, 9)),
                  {"m": (1, 2), "a": (1, 3), "b": (1, 4), "c": (1, 5),
                   "d": (2, 6), "e": (2, 7), "f": (2, 8)}).algebra()


@lit
def test_d4_a4_distinguished_by_coxeter():
    fa, fd = derived_fingerprint(_a4()), derived_fingerprint(_d4())
    res = compare_fingerprints(fa, fd)
    assert "coxeter_polynomial" in res["distinguished_by"]
    assert res["verdict"] == "distinguished"


@lit
def test_cospectral_trees_not_distinguished():
    fs = derived_fingerprint(_spider())
    fd = derived_fingerprint(_double_star())
    res = compare_fingerprints(fs, fd)
    assert res["distinguished_by"] == []
    assert res["verdict"] == "not distinguished by these invariants"
    # (these trees are NOT derived equivalent -- the honest-scope demonstration)


@selfcert
def test_fingerprint_never_claims_equivalent():
    res = compare_fingerprints(derived_fingerprint(_a4()), derived_fingerprint(_a4()))
    assert "equivalent" not in res["verdict"]        # necessary-condition language only


@selfcert
def test_singular_or_presentationless_fields_are_honest_errors():
    from quiverlab.families import truncated_polynomial
    from quiverlab.fields import GF
    fp = derived_fingerprint(truncated_polynomial(2, field=GF(7)))   # non-unimodular C
    # coxeter still computes (det C = 2 != 0 -> Phi integral, t+1); center/HH compute;
    # gl_dim is an honest infinite lower bound.
    assert isinstance(fp["cartan_det"], int) and fp["cartan_det"] == 2
    assert "gl.dim" in fp["gl_dim"] or ">=" in fp["gl_dim"]
