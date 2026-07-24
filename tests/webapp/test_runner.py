"""Task 7 — compute runner tests.

Adapted per the 2026-07-24 amendment: the reference family is
``QuantumCI(q=1, field=GF(2))`` (``truncated_polynomial`` is not in
``families()``); ``A.dim`` not ``A.dimension()``; HH tables expose
``kind``/``dims``/``engine``/``references`` (no ``.latex()``); every library
call that accepts ``verbose=`` gets ``verbose=False`` so the server never writes
worked-steps PDFs into ``quiverlab_traces/``. Both ``kind: "family"`` and
``kind: "quiver"`` requests must build and run.
"""
import json

from webapp.server.runner import build_algebra, run_spec, RunError
from webapp.server.schema import ComputeRequest


def _req(field, compute, artifacts=None, params=None):
    return ComputeRequest.model_validate({
        "schema": 1,
        "algebra": {"kind": "family", "family": "QuantumCI",
                    "params": params or {"q": 1}, "field": field},
        "compute": compute,
        "artifacts": artifacts or {"pdf": False, "tikz": False}})


def _quiver_req(field, compute, vertices, arrows, relations, artifacts=None):
    return ComputeRequest.model_validate({
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": vertices, "arrows": arrows,
                    "relations": relations, "field": field},
        "compute": compute,
        "artifacts": artifacts or {"pdf": False, "tikz": False}})


def test_runs_hh_over_gf(tmp_path):
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["hh_cohomology:0..4"])
    result = run_spec(req, tmp_path)
    assert "hh_cohomology" in result["results"]
    block = result["results"]["hh_cohomology"]
    assert isinstance(block["dims"], list)
    # Real HH-table contract (amendment pt 5): kind + engine, no latex.
    assert block["kind"] == "HH^"
    assert "engine" in block
    assert "latex" not in block
    assert (tmp_path / "result.json").exists()
    on_disk = json.loads((tmp_path / "result.json").read_text())
    assert on_disk["quiverlab_version"]


def test_build_algebra_family_uses_dim():
    # Amendment pt 4: the algebra exposes `.dim`, never `.dimension()`.
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["cartan"])
    A = build_algebra(req.algebra)
    assert A.dim == 4
    assert not hasattr(A, "dimension")


def test_runs_quiver_kind(tmp_path):
    # Amendment pt 7: schema v1 accepts `kind: "quiver"`; the runner must build
    # from vertices/arrows/relations exactly as docs/gui/runner.py::run_build.
    # k[x]/(x^3): one loop, cubic relation -> dim 3.
    req = _quiver_req({"kind": "GF", "p": 2, "n": 1},
                      ["hh_cohomology:0..3", "cartan"],
                      vertices=[1], arrows={"x": [1, 1]}, relations=["x*x*x"])
    result = run_spec(req, tmp_path)
    assert result["results"]["hh_cohomology"]["dims"] == [3, 2, 2, 2]
    assert result["results"]["cartan"]["matrix"] == [[3]]
    assert result["results"]["cartan"]["latex"].startswith(r"\begin{pmatrix}")


def test_center_and_global_dimension_shapes(tmp_path):
    # center -> (dim, basis) with str()-shipped basis entries; global_dimension
    # -> object with .exact/.value and a rich str() (amendment pt 5).
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["center", "global_dimension"])
    result = run_spec(req, tmp_path)
    center = result["results"]["center"]
    assert center["dim"] == 4
    assert all(isinstance(x, str) for row in center["basis"] for x in row)
    gdim = result["results"]["global_dimension"]
    assert isinstance(gdim["value"], int)
    assert isinstance(gdim["exact"], bool)
    assert isinstance(gdim["text"], str)


def test_unknown_family_is_run_error(tmp_path):
    req = ComputeRequest.model_validate({
        "schema": 1,
        "algebra": {"kind": "family", "family": "does_not_exist",
                    "params": {}, "field": {"kind": "GF", "p": 2}},
        "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}})
    try:
        run_spec(req, tmp_path)
        assert False, "expected RunError"
    except RunError as exc:
        assert exc.error_type in ("CatalogError", "FieldError")


def test_result_size_cap_refuses_before_writing(tmp_path):
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["hh_cohomology:0..4"])
    try:
        run_spec(req, tmp_path, result_max_bytes=8)   # absurdly small cap
        assert False, "expected RunError"
    except RunError as exc:
        assert exc.error_type == "ResultTooLarge"
    assert not (tmp_path / "result.json").exists()     # nothing written


def test_pdf_unavailable_is_recorded(tmp_path, monkeypatch):
    # Honest worked-steps contract: when pdf is requested with an HH item,
    # write_trace produces EITHER a real trace.pdf (pdflatex/tectonic on PATH) OR
    # a self-contained trace_steps.html fallback, and meta["pdf"] names whichever
    # landed. The old FALSE "trace subsystem absent" label must appear nowhere.
    monkeypatch.chdir(tmp_path)
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["hh_cohomology:0..3"],
               artifacts={"pdf": True, "tikz": False})
    art = tmp_path / "art"
    result = run_spec(req, art)
    meta_pdf = result["meta"]["pdf"]
    pdf_ok = (art / "trace.pdf").exists() and meta_pdf == "trace.pdf"
    html_ok = ("trace_steps.html" in meta_pdf
               and (art / "trace_steps.html").exists())
    assert pdf_ok or html_ok, f"neither trace.pdf nor trace_steps.html: {meta_pdf!r}"
    # The false label is gone -- from the payload and from result.json on disk.
    assert "trace subsystem absent" not in json.dumps(result)
    assert "trace subsystem absent" not in (art / "result.json").read_text()
    # write_trace was given an explicit out_dir, so the cwd is never touched.
    assert not (tmp_path / "quiverlab_traces").exists()


def test_pdf_requested_without_hh_item(tmp_path, monkeypatch):
    # pdf requested but no HH item computed -> honest label, no trace file.
    monkeypatch.chdir(tmp_path)
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["cartan"],
               artifacts={"pdf": True, "tikz": False})
    art = tmp_path / "art"
    result = run_spec(req, art)
    assert result["meta"]["pdf"] == "no traced computation requested (PDF covers HH worked steps)"
    assert not (art / "trace.pdf").exists()
    assert not (art / "trace_steps.html").exists()
    assert not (tmp_path / "quiverlab_traces").exists()


def test_duplicate_compute_kind_refused(tmp_path):
    # Finding 3: results is keyed by kind, so two same-kind items would collide;
    # refuse loudly BEFORE computing anything (no result.json written).
    req = _req({"kind": "GF", "p": 2, "n": 1},
               ["hh_cohomology:0..3", "hh_cohomology:0..5"])
    try:
        run_spec(req, tmp_path)
        assert False, "expected RunError"
    except RunError as exc:
        assert exc.error_type == "DuplicateComputeItem"
        assert "requested twice" in exc.message
    assert not (tmp_path / "result.json").exists()


def test_hh_homology_kind(tmp_path):
    # Finding 4b: hh_homology returns kind "HH_" with a dims list.
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["hh_homology:0..3"])
    result = run_spec(req, tmp_path)
    block = result["results"]["hh_homology"]
    assert block["kind"] == "HH_"
    assert isinstance(block["dims"], list) and len(block["dims"]) == 4


def test_bad_relation_is_relation_error(tmp_path):
    # Finding 4b: a quiver with an unparseable relation surfaces the library's
    # RelationError verbatim -- no traceback text in the payload.
    req = _quiver_req({"kind": "GF", "p": 2, "n": 1}, ["cartan"],
                      vertices=[1], arrows={"x": [1, 1]},
                      relations=["x*x*nonsense"])
    try:
        run_spec(req, tmp_path)
        assert False, "expected RunError"
    except RunError as exc:
        assert exc.error_type == "RelationError"
        assert "nonsense" in exc.message            # verbatim library message
        assert "Traceback" not in exc.message
        assert 'File "' not in exc.message
    assert not (tmp_path / "result.json").exists()


def test_tikz_written_when_requested(tmp_path):
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["cartan"],
               artifacts={"pdf": False, "tikz": True})
    run_spec(req, tmp_path)
    tikz = (tmp_path / "tikz.tex")
    assert tikz.exists()
    assert r"\begin{tikzpicture}" in tikz.read_text()


def test_references_present_in_result(tmp_path):
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["hh_cohomology:0..4"])
    result = run_spec(req, tmp_path)
    # Top-level, library-sourced bibliography entries for what the run used.
    assert isinstance(result["references"], list)
    assert result["references"], "QuantumCI HH cites at least one work"
    # Per-item BibTeX keys are carried through too.
    assert isinstance(result["results"]["hh_cohomology"]["references"], list)
    for entry in result["references"]:
        assert "key" in entry and "formatted" in entry
