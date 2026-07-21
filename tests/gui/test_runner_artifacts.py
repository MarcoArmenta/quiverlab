import io
import json
import contextlib

KRONECKER_CC = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1, 2],
                "arrows": {"a": [1, 2], "b": [1, 2]},
                "relations": [], "field": {"kind": "CC"}},
    "compute": ["hh_cohomology:0..2", "cartan"],
    "artifacts": {"pdf": True, "tikz": True},
}


def _run_all(runner, req):
    out = json.loads(runner.run_build(json.dumps(req)))
    assert out["ok"], out
    for spec in req["compute"]:
        assert json.loads(runner.compute_one(spec))["ok"]


def test_trace_html_after_hh(runner):
    _run_all(runner, KRONECKER_CC)
    html = runner.trace_html()
    assert "<html" in html.lower() and "References" in html


def test_trace_html_empty_without_hh(runner):
    req = dict(KRONECKER_CC, compute=["cartan"])
    _run_all(runner, req)
    assert runner.trace_html() == ""


def test_tikz(runner):
    _run_all(runner, KRONECKER_CC)
    assert runner.tikz().startswith(r"\begin{tikzpicture}")


def test_python_snippet_reproduces(runner):
    _run_all(runner, KRONECKER_CC)
    snippet = runner.python_snippet()
    assert 'Quiver(vertices=[1, 2], arrows={"a": (1, 2), "b": (1, 2)})' in snippet
    assert "A.hochschild_cohomology(2)" in snippet and "A.cartan_matrix()" in snippet
    # The GUI-to-library bridge must actually RUN.
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        exec(snippet, {"__name__": "__snippet__"})
    assert "4" in buf.getvalue()     # print(A.dim) for the Kronecker algebra


def test_result_bundle(runner):
    _run_all(runner, KRONECKER_CC)
    bundle = json.loads(runner.result_bundle())
    assert bundle["schema"] == 1
    assert bundle["request"]["algebra"]["kind"] == "quiver"
    assert bundle["quiverlab_version"]
    assert [r["invariant"] for r in bundle["results"]] == KRONECKER_CC["compute"]


def test_artifacts_before_build_are_empty(runner):
    assert runner.trace_html() == "" and runner.tikz() == ""
    assert runner.python_snippet() == ""


def test_result_bundle_records_failures(runner):
    req = dict(KRONECKER_CC, compute=["hh_cohomology:0..1", "determinant"])
    out = json.loads(runner.run_build(json.dumps(req)))
    assert out["ok"]
    assert json.loads(runner.compute_one("hh_cohomology:0..1"))["ok"]
    assert not json.loads(runner.compute_one("determinant"))["ok"]
    bundle = json.loads(runner.result_bundle())
    assert [r["invariant"] for r in bundle["results"]] == [
        "hh_cohomology:0..1", "determinant"]
    assert bundle["results"][1]["error"]["type"] == "RequestError"
