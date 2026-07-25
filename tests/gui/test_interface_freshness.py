"""Plan-10 freshness gate: pin the exact quiverlab surface docs/gui/runner.py
consumes (public API + the three sanctioned trace helpers). A failure here means
library drift: STOP and amend the plan; never patch the runner around it."""
import inspect

import quiverlab as ql


def test_quiver_and_algebra_build_surface():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    assert list(inspect.signature(Q.algebra).parameters) == [
        "relations", "field", "degree_bound", "trace"]
    A = Q.algebra(relations=["x*x*x"], field=ql.GF(2))
    assert isinstance(A.dim, int) and A.dim == 3
    assert list(A.quiver.vertices) == [1]
    assert dict(A.quiver.arrows) == {"x": (1, 1)}
    assert [str(r) for r in A.relations] == ["x*x*x"]


def test_invariant_surface():
    A = ql.Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
        relations=["x*x*x"], field=ql.GF(2))
    for name in ("hochschild_cohomology", "hochschild_homology"):
        params = inspect.signature(getattr(A, name)).parameters
        assert "verbose" in params and "trace" in params, name
    ev = []
    t = A.hochschild_cohomology(2, verbose=False, trace=ev)
    assert list(t.dims) == [3, 2, 2]
    assert t.kind == "HH^" and isinstance(t.engine, str)
    assert tuple(t.references) == ("bar",)
    assert len(ev) >= 1  # engines fill the explicit event sink
    m = A.cartan_matrix()
    assert m == [[3]]
    p = A.coxeter_polynomial()
    assert hasattr(p, "as_expr")
    g = A.global_dimension()
    assert hasattr(g, "exact") and hasattr(g, "value") and str(g)
    dim_z, basis = A.center()
    assert isinstance(dim_z, int) and isinstance(basis, list)
    assert isinstance(A.tikz(), str) and A.tikz().startswith(r"\begin{tikzpicture}")
    assert isinstance(A.citations(), tuple)


def test_module_surface():
    """Pin the module surface the no-code panel (Plan 26) consumes: the sided
    builders, A.module, and the Module invariant methods. Drift here means STOP
    and amend the plan -- never patch the runner around it."""
    A = ql.Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
        relations=["x*x*x"], field=ql.GF(2))
    # sided builders + constructor signatures
    for name in ("simple", "projective", "injective"):
        assert "side" in inspect.signature(getattr(A, name)).parameters, name
    assert list(inspect.signature(A.module).parameters) == [
        "dimension_vector", "arrow_action", "side", "name"]
    assert hasattr(A, "opposite") and hasattr(A, "ext")
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]})
    assert M.dimension_vector() == {1: 2} and M.dim == 2 and M.side == "right"
    for meth in ("radical", "top", "socle", "tau", "tau_minus",
                 "dualize", "transpose", "projective_resolution",
                 "injective_resolution", "injective_dimension"):
        assert hasattr(M, meth), meth
    # resolution / dimension shapes the runner reads
    pr = M.projective_resolution(3)
    assert hasattr(pr, "dimension_vectors") and hasattr(pr, "betti") and hasattr(pr, "pd")
    ir = M.injective_resolution(3)
    assert hasattr(ir, "dimension_vectors") and hasattr(ir, "injective_dimension")
    assert M.injective_dimension() is None or isinstance(M.injective_dimension(), int)
    # module Ext + citation step-ids the module blocks attach
    from quiverlab.modules.ext import ext_dims
    assert ext_dims(A, M, A.simple(1), 2) == [1, 1, 1]
    from quiverlab.trace.provenance import resolve_references
    for key in ("assem_book", "minimal_resolution", "module_ext"):
        assert len(resolve_references((key,))[0]) == 2, key


def test_fields_zoo_bibliography():
    assert list(inspect.signature(ql.GF).parameters) == ["q", "modulus"]
    assert str(ql.GF(8)) == "GF(2^3)"  # q = p**n spelling the runner uses
    zoo = list(ql.zoo(dim_max=12))
    assert zoo and all(hasattr(a, "dim") and hasattr(a, "quiver") for a in zoo)
    entries = list(ql.bibliography())
    assert entries and all(
        hasattr(e, "key") and hasattr(e, "bibtex_key") and hasattr(e, "formatted")
        for e in entries)


def test_sanctioned_trace_helpers():
    from quiverlab.trace.provenance import references_for, resolve_references
    from quiverlab.trace.render_html import render_html
    A = ql.Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
        relations=["x*x"], field=ql.CC)
    ev = []
    A.hochschild_cohomology(1, verbose=False, trace=ev)
    keys = references_for(ev)
    pairs = resolve_references(keys)
    assert all(len(p) == 2 for p in pairs)
    html = render_html(list(ev), title="t", references=pairs)
    assert isinstance(html, str) and "<html" in html.lower()
    assert isinstance(render_html([], title="t", references=()), str)


def test_error_types_exist():
    for name in ("QuiverlabError", "FieldError", "RelationError",
                 "NotFiniteDimensionalError", "AdmissibilityError",
                 "DepthLimitError", "ExactnessError"):
        assert isinstance(getattr(ql, name), type), name
