"""Plan-52 interface freshness gate. Pins the HH dispatch, HHTable, bar.py
reduction convention, the ported engine/bimodule reference, and the Nakayama /
bimodule-socle helpers. STOP on any drift."""
import inspect
import quiverlab as ql
from quiverlab.hochschild import bar, table


def test_hh_signatures_have_no_coefficients_yet():
    for name in ("hochschild_cohomology", "hochschild_homology"):
        params = inspect.signature(getattr(ql.Algebra, name)).parameters
        assert {"top", "max_cells", "engine", "auto_cs", "verbose", "trace"} <= set(params)
        assert "coefficients" not in params          # this task adds it (Task 3)
        assert "relative_to" not in params           # Task 7 adds it


def test_hhtable_attribute_set_frozen():
    t = ql.truncated_polynomial(2, field=ql.GF(7)).hochschild_cohomology(2)
    assert isinstance(t, table.HHTable)
    assert set(vars(t)) == {"dims", "kind", "top", "algebra_repr", "engine", "references"}
    #  ^ adding `coefficients` must default None and NOT appear here for the regular case


def test_bar_is_absolute_no_E_composability_filter():
    # SEMANTIC invariant (not a substring match): the reduced bar enumerates ALL
    # (m-1)^n multi-indices over the reduced basis {1..m-1} (unit index 0 dropped),
    # with NO E-composability filter. For a MULTI-VERTEX algebra an E-relative
    # complex (Task 7) would keep only E-composable tuples -- strictly fewer -- so
    # this exact count fails loudly if someone changes the reduction convention or
    # adds an E-filter to bar.py. (If _abar_tuples' signature has drifted, THIS
    # gate STOPs the plan -- which is the point.)
    A = ql.Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=ql.GF(7))
    m = A.dim                                          # kA2: dim 3 (e1, e2, a)
    for n in (1, 2, 3):
        tuples = list(bar._abar_tuples(m, n))
        assert len(tuples) == (m - 1) ** n             # 2, 4, 8 -- absolute, unfiltered
        assert all(len(t) == n and all(1 <= i < m for i in t) for t in tuples)


def test_ported_gfp_bimodule_reference_present():
    from quiverlab.engine import bimodule as eb
    for n in ("Bimodule", "regular_bimodule", "dual_bimodule", "twisted_bimodule",
              "check_bimodule", "hochschild_homology_with_coefficients",
              "hochschild_cohomology_with_coefficients"):
        assert hasattr(eb, n), f"engine.bimodule.{n} missing (reference oracle drift)"


def test_nakayama_generic_route_and_socle_helpers():
    # DD1b: twisted_by_nakayama uses the GENERIC trace-form route
    # (nakayama_automorphism_generic, A's OWN basis) -- NOT
    # Algebra.nakayama_automorphism(), which over GF(p) returns the ENGINE
    # unit-adapted basis (core/algebra.py:896-897). Pin the generic helper's
    # presence + A-basis shape, and the bimodule-socle helper for mod_socle.
    from quiverlab.invariants.frobenius import nakayama_automorphism_generic
    A = ql.truncated_polynomial(3, field=ql.CC)        # Frobenius => nu exists
    nu = nakayama_automorphism_generic(A)              # m x m, columns = images, A-basis
    assert len(nu) == A.dim and len(nu[0]) == A.dim
    from quiverlab.families.trivial_extension import _bimodule_socle
    assert callable(_bimodule_socle)
