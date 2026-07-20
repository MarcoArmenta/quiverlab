"""Plan-07 freshness gate (spec §3.7-3.9). Pins the exact shape of everything
Plan 07 builds ON TOP OF from merged Plans 03-06. Any drift here means STOP and
reconcile the plan before touching Plan-07 source. This test has NO Plan-07
imports on purpose -- it must be writable and runnable before any Plan-07 code
exists."""
from dataclasses import fields as dc_fields

import pytest

import quiverlab as ql


def _field_names(cls):
    return {f.name for f in dc_fields(cls)}


# --- shipped Groebner events (Plan 03) ---------------------------------------
def test_groebner_events_shape():
    from quiverlab.groebner.events import Dispatch, ReductionStep
    assert _field_names(Dispatch) == {"route", "reason", "n_relations"}
    assert _field_names(ReductionStep) == {"word", "rule_lead", "before", "after"}


# --- Chouhy-Solotar trace events (Plan 04) -----------------------------------
def test_cs_trace_events_shape():
    from quiverlab.resolutions_cs.trace import (
        AmbiguityEvent, ResolutionTerm, DifferentialEvent, LiftStep,
    )
    assert _field_names(AmbiguityEvent) == {"degree", "chain_words"}
    assert _field_names(ResolutionTerm) == {"degree", "n_generators", "collapsed_dim"}
    assert _field_names(DifferentialEvent) == {"degree", "chain", "terms"}
    assert _field_names(LiftStep) == {"degree", "kind", "detail"}


# --- HHTable shape (Plan 01, still frozen) -----------------------------------
def test_hhtable_shape():
    from quiverlab.hochschild.table import HHTable
    t = HHTable([1, 0], "HH^", "an algebra")
    assert t.dims == [1, 0]
    assert hasattr(t, "engine") and isinstance(t.engine, str)
    # Plan 07 will attach `.references`; it must NOT already be a slotted/blocked attr.
    t.references = ("X",)
    assert t.references == ("X",)


# --- Hochschild method contract Plan 07 EXTENDS (must not regress) ------------
def test_hochschild_method_contract_is_frozen():
    """Plan 07 only ADDS `verbose=` to hochschild_cohomology/homology; the merged
    `engine='cs'`/`auto_cs`/`trace` routing, the unknown-engine QuiverlabError, and
    the `.references = self.citations()` (family+engine) union are contracts it must
    PRESERVE (Tasks 8 & 11). Pinned here so a regression (dropping 'cs', raising
    ValueError, clobbering .references engine-only, or un-wrapping the fast list)
    fails at Task 1, not silently at Task 13."""
    import inspect
    from quiverlab.errors import QuiverlabError
    from quiverlab.core.algebra import Algebra
    params = inspect.signature(Algebra.hochschild_cohomology).parameters
    assert "auto_cs" in params and "trace" in params, params
    # unknown engine raises QuiverlabError (NOT a bare ValueError; QuiverlabError is
    # not a ValueError subclass), and its guidance still names 'cs' as a valid engine.
    A = ql.truncated_polynomial(2, field=ql.CC)
    with pytest.raises(QuiverlabError) as ei:
        A.hochschild_cohomology(1, engine="definitely-not-an-engine")
    assert "cs" in str(ei.value)
    # .references is the family+engine union (== citations()), not engine-only: use a
    # stamped catalog family (NakayamaAlgebra stamps ("nakayama","assem_book")).
    N = ql.NakayamaAlgebra(n=3, l=2)
    refs = N.hochschild_cohomology(2).references
    assert refs == N.citations() and "nakayama" in refs and "bar" in refs
    # the fast GF(p) path returns an HHTable whose `.engine` is a string (list wrapped).
    fast = ql.truncated_polynomial(2, field=ql.GF(2)).hochschild_cohomology(2)
    assert isinstance(fast.engine, str) and fast.engine


# --- bar-complex entry points Plan 07 instruments ----------------------------
def test_bar_entry_points_exist():
    from quiverlab.hochschild import bar
    for name in ("coboundary_matrix", "hochschild_cohomology_dims",
                 "boundary_matrix", "hochschild_homology_dims"):
        assert hasattr(bar, name), name


# --- Plan-06 citation registry (bibliography()) ------------------------------
def test_bibliography_registry_has_needed_keys():
    """Plan 06's bibliography() returns a Bibliography dataclass: a `.keys` TUPLE
    (lowercase REGISTRY keys) and iteration (__iter__) yielding entry views with
    .key / .formatted / .bibtex_key (and .doi/.arxiv/.topic/.annotation). There is
    NO .keys() method and NO subscripting. Plan-07's provenance map (Task 11) keys
    off these registry names -- if Plan 06 renames them, update
    quiverlab.trace.provenance.ENGINE_REFERENCES and this set together."""
    assert hasattr(ql, "bibliography"), (
        "Plan 06 must export bibliography() (the citation registry Plan 07 renders)")
    bib = ql.bibliography()
    assert hasattr(bib, "keys") and isinstance(bib.keys, tuple), (
        "bibliography().keys must be a TUPLE attribute (not a method)")
    needed = {"bar", "bardzell", "chouhy_solotar"}
    missing = needed - set(bib.keys)
    assert not missing, f"bibliography() missing registry keys Plan 07 needs: {sorted(missing)}"
    # entry views by registry key (the __iter__ protocol Plan 06 is adding this round)
    by_key = {e.key: e for e in bib}
    assert needed <= set(by_key), "bibliography() iteration must yield the needed keys"
    bar = by_key["bar"]
    for attr in ("key", "formatted", "bibtex_key"):
        assert hasattr(bar, attr), f"entry view missing .{attr}"
    # The `bar` registry key is backed by the genuine Hochschild (1945) reference the
    # golden worked-steps trace displays. STOP if Plan 06 lands without it.
    assert bar.bibtex_key == "Hochschild1945", (
        "Plan 06's `bar` entry must carry bibtex_key == 'Hochschild1945' "
        "(G. Hochschild, Ann. of Math. 46 (1945), 58-67); coordinate its addition")
    assert "Hochschild" in bar.formatted


# --- Plan-06 family catalog --------------------------------------------------
def test_family_catalog_exists():
    assert hasattr(ql, "families"), "Plan 06 must export families() (catalog discovery)"
    cat = ql.families()
    # FamilyListing is always truthy; assert the catalog actually enumerates names.
    assert cat.names(), "families() catalog is empty (FamilyListing.names() returned nothing)"


# --- verbose flag is a plain module attribute Plan 07/09 toggle --------------
def test_verbose_flag_is_settable():
    # Before Plan 07 the attribute may be absent; after Task 2 it exists and defaults True.
    # This assertion documents the contract; it is xfail-until-Task-2.
    if not hasattr(ql, "verbose"):
        pytest.xfail("quiverlab.verbose is introduced in Plan-07 Task 2")
    prev = ql.verbose
    try:
        ql.verbose = False
        assert ql.verbose is False
    finally:
        ql.verbose = prev
