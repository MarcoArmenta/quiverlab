"""Weak positivity by the Ovsienko box (cited, branch-and-bound); weak nonnegativity by
the classified hypercritical list (primary) + a sound witness-finder. Oracles: Dynkin =>
weakly positive; Euclidean (extended Dynkin) => weakly nonnegative not weakly positive
(isotropic witness q=0); 3-Kronecker (small witness) and T_{2,3,7} (LIST entry, sincere
10-var witness) => not weakly nonnegative. Refs: Bongartz 1984; von Hohne 1996; BJP 2019."""
import pytest

from quiverlab import CC, Quiver, linear_path_algebra
from quiverlab.invariants.tits import (as_unit_form, is_weakly_nonnegative,
                                       is_weakly_positive)
from quiverlab.invariants.tits import tits_form_combinatorial as titsc

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _kron(m):
    return Quiver([1, 2], {f"a{k}": (1, 2) for k in range(m)}).algebra(field=CC)


def _star(arms):
    verts, ars, nxt = [0], {}, 1
    for ai, length in enumerate(arms):
        prev = 0
        for _ in range(length):
            verts.append(nxt)
            ars[f"e{ai}_{nxt}"] = (prev, nxt)
            prev = nxt
            nxt += 1
    return Quiver(verts, ars).algebra(field=CC)


@lit
def test_dynkin_weakly_positive():
    for A in (linear_path_algebra(5, field=CC), _star([1, 2, 4])):   # A5, E8 (n<=8)
        f = as_unit_form(A)
        assert is_weakly_positive(f).holds is True
        assert is_weakly_nonnegative(f).holds is True                # positive => nonneg


@lit
def test_euclidean_weakly_nonnegative_not_positive():
    # 2-Kronecker = ~A1 (isotropic (1,1)); affine E8 = ~E8 = star arms [1,2,5] (9 verts)
    for A in (_kron(2), _star([1, 2, 5])):
        f = as_unit_form(A)
        wp = is_weakly_positive(f)
        assert wp.holds is False and titsc(A, list(wp.witness)) <= 0   # isotropic root
        assert is_weakly_nonnegative(f).holds is True                  # certified (PSD)
    assert titsc(_kron(2), [1, 1]) == 0                                # the isotropic delta


@lit
def test_wild_not_weakly_nonnegative():
    # 3-Kronecker (n=2, small witness / list) and T_{2,3,7}=star arms [1,2,6]
    # (n=10, decided by the LIST entry -- its witness is SINCERE on all 10 vertices,
    # (12,6,8,4,10,9,7,6,4,2), q=-1: unreachable by any support/box cap).
    for A in (_kron(3), _star([1, 2, 6])):
        f = as_unit_form(A)
        v = is_weakly_nonnegative(f)
        assert v.holds is False
        assert titsc(A, list(v.witness)) < 0                          # strict negative


@selfcert
def test_witness_is_exact_and_nonnegative():
    f = as_unit_form(_kron(3))
    v = is_weakly_nonnegative(f)
    assert all(x >= 0 for x in v.witness) and any(x > 0 for x in v.witness)
    assert v.witness_value == titsc(_kron(3), list(v.witness))


@selfcert
def test_weak_positivity_budget_is_honest_None():
    # a tiny budget must return None (undecided) OR the fast PD certificate's True,
    # never a silent False.
    f = as_unit_form(linear_path_algebra(5, field=CC))
    v = is_weakly_positive(f, budget=1)
    assert v.holds in (True, None)
    if v.holds is None:
        assert v.reason.startswith("budget")


@selfcert
def test_weak_nonneg_outside_coverage_is_None_not_guessed_true():
    # No form gets a GUESSED True. A non-PSD form not caught by the list and starved
    # of the witness-finder budget returns honest None (never True); the same form
    # with a real budget returns a sound False witness. And PSD forms DO get True.
    f5 = as_unit_form(_kron(5))                # 5-Kronecker: indefinite, not listed
    starved = is_weakly_nonnegative(f5, budget=0)
    assert starved.holds is None               # honest -- NEVER a guessed True
    full = is_weakly_nonnegative(f5)           # default budget: finds (1,1), q=-3
    assert full.holds is False and titsc(_kron(5), list(full.witness)) < 0
    # PSD is the sound True certificate (no coverage caveat needed)
    assert is_weakly_nonnegative(as_unit_form(_kron(2))).holds is True
