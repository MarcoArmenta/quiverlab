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
