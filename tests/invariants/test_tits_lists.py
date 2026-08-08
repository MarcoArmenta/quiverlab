"""The classified critical (Euclidean) / hypercritical lists. Critical: q(radical)=0,
not weakly positive, weakly nonnegative (cross-oracle). Hypercritical: q(defect)<0, drives
the weak-nonnegativity decision -- including T_{2,3,7} (10 vars, defect (12,6,8,4,10,9,7,6,
4,2), q=-1, verified in exact integer arithmetic). Refs: von Hohne 1996; de la Pena
Banach-26 (1990); BJP 2019."""
import pytest

from quiverlab.invariants._tits_lists import (CRITICAL_FORMS, HYPERCRITICAL_COVERAGE,
                                              HYPERCRITICAL_FORMS,
                                              matches_hypercritical)
from quiverlab.invariants.tits import UnitForm, is_weakly_nonnegative, is_weakly_positive

pytestmark = pytest.mark.oracle_literature


@pytest.mark.parametrize("entry", CRITICAL_FORMS, ids=[c["name"] for c in CRITICAL_FORMS])
def test_critical_radical_isotropic_and_not_weakly_positive(entry):
    f = UnitForm(len(entry["gram"]), entry["gram"], tuple(range(len(entry["gram"]))))
    assert f.evaluate(entry["radical"]) == 0                 # radical: q(z) = 0
    assert is_weakly_positive(f).holds is False              # critical => not weakly pos
    assert is_weakly_nonnegative(f).holds is True            # Euclidean => weakly nonneg


@pytest.mark.parametrize("entry", HYPERCRITICAL_FORMS,
                         ids=[c["name"] for c in HYPERCRITICAL_FORMS] or ["<deferred>"])
def test_hypercritical_defect_negative_and_decided_false(entry):
    f = UnitForm(len(entry["gram"]), entry["gram"], tuple(range(len(entry["gram"]))))
    assert f.evaluate(entry["defect"]) < 0                   # q(defect) < 0
    v = is_weakly_nonnegative(f)
    assert v.holds is False                                  # list-decided False + witness
    assert f.evaluate(v.witness) < 0                         # a real negative witness


def test_T237_is_a_10_variable_hypercritical_entry():
    # THE refutation, encoded so no one reintroduces a <=9 cap: T_{2,3,7} is hypercritical
    # with 10 variables and a SINCERE defect vector.
    names = {c["name"] for c in HYPERCRITICAL_FORMS}
    assert "T237" in names
    e = next(c for c in HYPERCRITICAL_FORMS if c["name"] == "T237")
    assert len(e["gram"]) == 10 and all(x > 0 for x in e["defect"])   # sincere, 10 vars


def test_matcher_is_relabelling_invariant():
    # matches_hypercritical must recognise a permuted copy of a listed form.
    e = next(c for c in HYPERCRITICAL_FORMS if c["name"] == "T334")
    n = len(e["gram"])
    perm = list(range(1, n)) + [0]                       # a cyclic relabelling
    g = e["gram"]
    pg = tuple(tuple(g[perm[i]][perm[j]] for j in range(n)) for i in range(n))
    f = UnitForm(n, pg, tuple(range(n)))
    m = matches_hypercritical(f)
    assert m is not None and m["name"] == "T334"
    assert f.evaluate(m["defect"]) < 0                   # defect remapped to f's indices


def test_coverage_is_documented():
    # honest scope: coverage is an explicit datum the verification page cites.
    assert HYPERCRITICAL_COVERAGE is not None
