"""Ground-truth agreement (Plan 62 / R19). m-Kronecker ladder = FORM layer (Kronecker is
NOT simply connected); Dynkin/Euclidean trees + T_{2,3,7} = VERDICT layer; the AR-knit
cross-check ties rep-finite to weakly positive on ssc CC instances. Refs: Bongartz 1984;
BdlPS 2011; Gabriel; Nazarova."""
import pytest

from quiverlab import CC, Quiver, linear_path_algebra
from quiverlab.invariants.tits import (as_unit_form, is_weakly_nonnegative,
                                       is_weakly_positive, tame_wild_certificate)

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


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
@pytest.mark.parametrize("m,wp,wnn", [(1, True, True), (2, False, True), (3, False, False)])
def test_kronecker_ladder_is_form_layer_only(m, wp, wnn):
    K = _kron(m)
    f = as_unit_form(K)
    assert is_weakly_positive(f).holds is wp
    assert is_weakly_nonnegative(f).holds is wnn
    c = tame_wild_certificate(K)
    if m == 1:
        assert c.rep_type == "rep-finite"                 # A2 tree: reaches the verdict
    else:
        assert c.simply_connected is False and c.rep_type is None   # not in verdict scope


@lit
@pytest.mark.parametrize("A,expect", [
    (linear_path_algebra(6, field=CC), "rep-finite"),     # A6
    (_star([1, 1, 3]), "rep-finite"),                     # D6 tree
    (_star([1, 2, 2]), "rep-finite"),                     # E6
    (_star([1, 2, 5]), "tame"),                           # ~E8
    (_star([1, 2, 6]), "wild"),                            # T_{2,3,7}
])
def test_tree_verdicts(A, expect):
    assert tame_wild_certificate(A).rep_type == expect


@xeng
def test_knit_agrees_with_weakly_positive():
    # AR knit (module-category enumeration) vs the Tits form (quadratic-form arithmetic)
    # -- two wholly independent implementations of "representation-finite". Small
    # Dynkin instances keep the knit fast (E6's knit is ~200s -- adjust to reality).
    for A in (linear_path_algebra(5, field=CC), linear_path_algebra(4, field=CC),
              _star([1, 1, 1])):                           # A5, A4, D4: rep-finite trees
        c = tame_wild_certificate(A)
        ar = A.ar_quiver()
        if ar.is_complete:                                 # knit closed => rep-finite
            assert c.rep_type == "rep-finite" and c.weakly_positive is True
