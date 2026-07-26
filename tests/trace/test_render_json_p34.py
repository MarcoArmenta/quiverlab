"""Plan 34 -- the JSON worked-steps leg: the COMPLETE machine record.

Marco's contract: every worked-steps computation delivers PDF + HTML + JSON, all
cached/served. This pins the JSON renderer (quiverlab.trace.render_json): a
schema-versioned, deterministic, exact serialization of the recorder's event
stream; the writer's three-artifact output; the hpc ``render --format json`` byte
match; and honest handling of a record-elided matrix.
"""
import json
import pathlib
import shutil

import pytest

import quiverlab
from quiverlab import truncated_polynomial, CC
from quiverlab.modules.qpa_module import _json_entry
from quiverlab.trace.recorder import Trace
from quiverlab.trace import writer as W
from quiverlab.trace.render_json import render_json, TRACE_SCHEMA
from quiverlab.trace.events import RankStep, ModuleDifferential
from quiverlab.trace.provenance import references_for, resolve_references


def _hh_events():
    """HH^*(k[x]/(x^3)) to degree 2 -- Dispatch + ResolutionTerm/RankStep, one of
    which (the degree-2 differential) is a genuinely wide 12-column matrix."""
    A = truncated_polynomial(3, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    refs = resolve_references(references_for(list(tr)))
    return A, list(tr), table, refs


# --------------------------------------------------------------------------- #
# Envelope + version.
# --------------------------------------------------------------------------- #
def test_schema_envelope_and_version():
    A, ev, _table, refs = _hh_events()
    obj = json.loads(render_json(ev, title="HH", references=refs, algebra=A))
    assert obj["quiverlab_trace_schema"] == TRACE_SCHEMA == 1
    assert obj["library_version"] == getattr(quiverlab, "__version__", "unknown")
    assert obj["title"] == "HH"
    assert isinstance(obj["events"], list) and obj["events"]
    # citations carry the (key, formatted) provenance the other renderers print.
    assert obj["citations"] and all(
        set(c) == {"key", "formatted"} for c in obj["citations"])
    # every event is typed by its recorder class name.
    assert {e["type"] for e in obj["events"]} <= {
        "Dispatch", "ResolutionTerm", "RankStep"}


def test_event_order_preserved():
    """The event ORDER is the one semantic ordering -- a list, never reordered by
    sort_keys (which only sorts dict keys)."""
    A, ev, _t, refs = _hh_events()
    obj = json.loads(render_json(ev, references=refs, algebra=A))
    got = [e["type"] for e in obj["events"]]
    want = [type(e).__name__ for e in ev]
    assert got == want


# --------------------------------------------------------------------------- #
# Determinism: two renders of identical input are byte-equal.
# --------------------------------------------------------------------------- #
def test_byte_deterministic():
    A, ev, _t, refs = _hh_events()
    a = render_json(ev, title="HH", references=refs, algebra=A)
    b = render_json(ev, title="HH", references=refs, algebra=A)
    assert a == b


def test_deterministic_independent_of_dict_key_order():
    """sort_keys makes an event's dict field key order irrelevant: the same dimvec
    entered in two orders serializes to identical bytes."""
    from quiverlab.trace.events import ModuleTerm
    a = ModuleTerm(degree=0, kind="projective", sym="P", summands=[1, 2],
                   dim=3, dimvec={"1": 1, "2": 2})
    b = ModuleTerm(degree=0, kind="projective", sym="P", summands=[1, 2],
                   dim=3, dimvec={"2": 2, "1": 1})       # reversed insertion order
    assert render_json([a]) == render_json([b])


# --------------------------------------------------------------------------- #
# Exact entries: round-trip + matrices match the recorded event data.
# --------------------------------------------------------------------------- #
def test_matrices_round_trip_and_match_event_data():
    A, ev, _t, refs = _hh_events()
    obj = json.loads(render_json(ev, references=refs, algebra=A))   # clean json.loads
    ranks = [e for e in ev if isinstance(e, RankStep)]
    js_ranks = [e for e in obj["events"] if e["type"] == "RankStep"]
    assert ranks and len(ranks) == len(js_ranks)
    for src, out in zip(ranks, js_ranks):
        assert src.matrix is not None                     # these are recorded in full
        assert out["nrows"] == src.nrows and out["ncols"] == src.ncols
        assert out["rank"] == src.rank
        # every entry equals the exact _json_entry rendering of the recorded string.
        for i in range(src.nrows):
            for j in range(src.ncols):
                assert out["matrix"][i][j] == _json_entry(src.matrix[i][j])


def test_exact_entry_grammar_int_fraction_and_display_only():
    """Integer strings normalise to JSON ints, ``a/b`` to a fraction STRING, and a
    non-re-enterable entry (a GF(p^n) element rendering) flags the event
    ``display_only`` -- never coerced, never dropped."""
    frac = RankStep(degree=0, side="chain", nrows=1, ncols=3, rank=1, field="CC",
                    matrix=[["1/2", "3", "-2"]])
    exotic = RankStep(degree=1, side="cochain", nrows=1, ncols=2, rank=1,
                      field="GF(4)", matrix=[["1", "x^1"]])
    obj = json.loads(render_json([frac, exotic]))
    e0, e1 = obj["events"]
    assert e0["matrix"] == [["1/2", 3, -2]]              # int / "p/q" grammar
    assert "display_only" not in e0                       # all re-enterable
    assert e1["matrix"] == [[1, "x^1"]]
    assert e1["display_only"] is True                     # x^1 is not re-enterable


def test_module_differential_matrix_serialized():
    """A ModuleDifferential (the module worked-steps carrier) also emits its matrix
    through the exact grammar."""
    md = ModuleDifferential(degree=1, kind="projective", sym="P", symbol="d_1",
                            dom_summands=[1], cod_summands=[1], nrows=2, ncols=2,
                            field="GF(2)", matrix=[["1", "0"], ["0", "1"]])
    obj = json.loads(render_json([md]))
    (e,) = obj["events"]
    assert e["type"] == "ModuleDifferential"
    assert e["matrix"] == [[1, 0], [0, 1]]
    assert e["symbol"] == "d_1" and e["cod_is_module"] is False


# --------------------------------------------------------------------------- #
# The 250k record-elision backstop: shape+rank+note, matrix null, elided true.
# --------------------------------------------------------------------------- #
def test_elided_backstop_serialized_honestly():
    elided = RankStep(degree=5, side="cochain", nrows=600, ncols=600, rank=42,
                      field="GF(2)", matrix=None, elided=True,
                      note="[shape 600x600, rank 42 -- elided]")
    (e,) = json.loads(render_json([elided]))["events"]
    assert e["elided"] is True
    assert e["matrix"] is None                            # body honestly absent
    assert e["nrows"] == 600 and e["ncols"] == 600 and e["rank"] == 42
    assert e["note"] == "[shape 600x600, rank 42 -- elided]"
    assert "display_only" not in e                        # nothing fabricated


def test_recorder_backstop_produces_an_elided_event():
    """The recorder's own helper elides past MATRIX_ELISION_CELLS -- render_json then
    serializes that elided event honestly (matrix null)."""
    from quiverlab.trace.recorder import rankstep, MATRIX_ELISION_CELLS

    class _Dom:
        name = "GF(2)"
    n = MATRIX_ELISION_CELLS + 1
    rs = rankstep(0, "cochain", D=None, nrows=n, ncols=1, rank=0, dom=_Dom())
    assert rs.elided and rs.matrix is None
    (e,) = json.loads(render_json([rs]))["events"]
    assert e["elided"] is True and e["matrix"] is None and e["nrows"] == n


# --------------------------------------------------------------------------- #
# The writer persists all three artifacts (JSON always; tex always; pdf-or-html).
# --------------------------------------------------------------------------- #
def test_writer_writes_json_beside_html_fallback(tmp_path, monkeypatch):
    A, ev, table, _refs = _hh_events()
    monkeypatch.setattr(W, "have_latex", lambda: None)      # no toolchain -> HTML
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    path = W.write_trace(ev, table, algebra=A, kind="HH^", top=2, out_dir=str(tmp_path))
    stem = pathlib.Path(path)
    assert stem.suffix == ".html"
    assert stem.with_suffix(".tex").exists()
    j = stem.with_suffix(".json")
    assert j.exists()
    obj = json.loads(j.read_text())
    assert obj["quiverlab_trace_schema"] == 1 and obj["events"]


def test_writer_writes_json_beside_pdf(tmp_path, monkeypatch):
    A, ev, table, _refs = _hh_events()
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")

    def fake_compile(tex, out_pdf, engine):
        pathlib.Path(out_pdf).write_bytes(b"%PDF-1.5 fake\n%%EOF\n")
        return 1

    monkeypatch.setattr(W, "_compile_pdf", fake_compile)
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    path = W.write_trace(ev, table, algebra=A, kind="HH^", top=2, out_dir=str(tmp_path))
    stem = pathlib.Path(path)
    assert stem.suffix == ".pdf" and stem.exists()
    assert stem.with_suffix(".tex").exists()
    assert stem.with_suffix(".json").exists()
    assert json.loads(stem.with_suffix(".json").read_text())["events"]


def test_writer_json_is_pure_function_of_events(tmp_path, monkeypatch):
    """The JSON sidecar is byte-identical whichever PDF/HTML branch the writer takes
    (a pure function of the events), so the machine record never depends on the
    toolchain."""
    A, ev, table, refs = _hh_events()
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)

    monkeypatch.setattr(W, "have_latex", lambda: None)
    p1 = pathlib.Path(W.write_trace(ev, table, algebra=A, kind="HH^", top=2,
                                    references=refs, out_dir=str(tmp_path / "html")))
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")
    monkeypatch.setattr(W, "_compile_pdf",
                        lambda tex, o, e: (pathlib.Path(o).write_bytes(b"%PDF\n"), 1)[1])
    p2 = pathlib.Path(W.write_trace(ev, table, algebra=A, kind="HH^", top=2,
                                    references=refs, out_dir=str(tmp_path / "pdf")))
    assert (p1.with_suffix(".json").read_text()
            == p2.with_suffix(".json").read_text())


@pytest.mark.skipif(W.have_latex() is None, reason="no LaTeX toolchain on PATH")
def test_writer_real_pdf_plus_json(tmp_path):
    A, ev, table, _refs = _hh_events()
    path = pathlib.Path(W.write_trace(ev, table, algebra=A, kind="HH^", top=2,
                                      out_dir=str(tmp_path)))
    assert path.suffix == ".pdf" and path.read_bytes()[:5] == b"%PDF-"
    assert path.with_suffix(".json").exists() and path.with_suffix(".tex").exists()


# --------------------------------------------------------------------------- #
# The hpc pipeline: spec.run promotes trace.json; `render --format json` emits
# exactly those bytes.
# --------------------------------------------------------------------------- #
_HPC_SPEC = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["hh_cohomology:0..2"],
    "artifacts": {"pdf": True, "tikz": False},
}


def test_spec_run_promotes_trace_json(tmp_path, monkeypatch):
    monkeypatch.setattr(W, "have_latex", lambda: None)
    from quiverlab.hpc.spec import run as spec_run
    art = tmp_path / "art"
    art.mkdir()
    spec_run(_HPC_SPEC, art)
    tj = art / "trace.json"
    assert tj.exists(), "the worked-steps run must promote trace.json"
    obj = json.loads(tj.read_text())
    assert obj["quiverlab_trace_schema"] == 1 and obj["events"]


def test_hpc_render_format_json_byte_matches(tmp_path, monkeypatch):
    monkeypatch.setattr(W, "have_latex", lambda: None)
    from quiverlab.hpc.spec import run as spec_run
    from quiverlab.hpc import report
    art = tmp_path / "art"
    art.mkdir()
    spec_run(_HPC_SPEC, art)
    promoted = (art / "trace.json").read_text()
    out, fmt = report.render(str(art / "result.json"), art / "emitted.json", fmt="json")
    assert fmt == "json"
    assert out.read_text() == promoted            # the SAME bytes, verbatim


def test_hpc_render_json_default_out_name():
    from quiverlab.hpc import report
    assert report.default_out_name("json") == "trace.json"


def test_hpc_render_json_missing_sibling_raises(tmp_path):
    """--format json with no trace.json beside result.json fails LOUDLY (no
    worked-steps event stream was produced)."""
    from quiverlab.hpc import report
    (tmp_path / "result.json").write_text(
        json.dumps({"quiverlab_version": "x", "results": {}}))
    with pytest.raises(report.ReportError):
        report.render(str(tmp_path / "result.json"), tmp_path / "o.json", fmt="json")


def test_hpc_render_json_from_dict_raises(tmp_path):
    """An in-memory result dict has no directory to look in; --format json says so."""
    from quiverlab.hpc import report
    with pytest.raises(report.ReportError):
        report.render({"quiverlab_version": "x", "results": {}},
                      tmp_path / "o.json", fmt="json")


@pytest.mark.skipif(shutil.which("tectonic") is None and shutil.which("pdflatex") is None,
                    reason="no LaTeX toolchain on PATH")
def test_hpc_render_json_alongside_real_pdf(tmp_path):
    from quiverlab.hpc.spec import run as spec_run
    from quiverlab.hpc import report
    art = tmp_path / "art"
    art.mkdir()
    spec_run(_HPC_SPEC, art)
    assert (art / "trace.pdf").exists() and (art / "trace.json").exists()
    out, fmt = report.render(str(art / "result.json"), art / "emitted.json", fmt="json")
    assert out.read_text() == (art / "trace.json").read_text()
