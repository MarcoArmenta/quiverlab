"""Result .references plumbing: HHTable carries the engine's citation keys."""
from quiverlab import citations
from quiverlab.combinat import Quiver
from quiverlab.fields import CC, GF
from quiverlab.hochschild.table import HHTable


def test_hhtable_references_default_empty_and_repr():
    t = HHTable([1, 0], "HH^", "A")
    assert t.references == ()
    t2 = HHTable([1, 0], "HH^", "A", references=("bar",))
    assert t2.references == ("bar",) and "bar" in repr(t2)
    assert HHTable([1, 0], "HH^", "A") == HHTable([1, 0], "HH^", "A", references=("bar",))


def test_bar_path_reports_bar_key():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=CC)   # kA_2
    hh = A.hochschild_cohomology(1)
    assert "bar" in hh.references
    for key in hh.references:                                            # all in registry
        citations.reference(key)


def test_result_references_include_family_and_engine():
    """spec §3.9: table.references names engine + family keys (ops added by Plan 04/05)."""
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(5))
    A._family_citations = ("nakayama", "assem_book")                     # simulate a family stamp
    hh = A.hochschild_cohomology(1)
    assert "nakayama" in hh.references and "bar" in hh.references        # family + engine
    for key in A.citations():
        citations.reference(key)                                         # every key valid
