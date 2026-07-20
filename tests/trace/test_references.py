"""Provenance -> References section, resolved through Plan 06's bibliography()."""
import pathlib

import quiverlab
from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace.provenance import (
    references_for, resolve_references, ENGINE_REFERENCES,
)
from quiverlab.trace.events import Dispatch


def test_bar_engine_maps_to_registry_key():
    tr = [Dispatch(route="normalized bar complex", reason="", n_relations=1)]
    assert references_for(tr) == ("bar",)                # REGISTRY key, not a BibTeX id
    assert ENGINE_REFERENCES["bardzell"] == ("bardzell",)


def test_resolve_uses_plan06_registry():
    pairs = resolve_references(("bar",))                 # registry key in
    assert len(pairs) == 1
    bibtex_key, entry = pairs[0]                          # (bibtex_key, formatted) out
    assert bibtex_key == "Hochschild1945"               # backed by Plan 06's `bar` entry
    assert "Hochschild" in entry


def test_resolve_unknown_key_raises_loudly():
    import pytest
    with pytest.raises(KeyError):
        resolve_references(("not_a_registry_key",))


def test_result_object_carries_the_citations_union():
    A = truncated_polynomial(2, field=CC)
    table = A.hochschild_cohomology(2, trace=Trace())
    # .references is the merged family+engine union (== A.citations()); this no-family
    # fixture makes that exactly ("bar",). Task 11 must NOT overwrite it engine-only.
    assert table.references == A.citations() == ("bar",)


def test_cs_engine_yields_chouhy_solotar_citation():
    # A CS (engine="cs") run must emit a Dispatch(route="chouhy-solotar") into the trace
    # so the worked-steps References section cites Chouhy-Solotar (previously silently empty).
    A = truncated_polynomial(2, field=CC)   # admissible algebra the CS engine accepts
    tr = []
    table = A.hochschild_cohomology(2, engine="cs", trace=tr)
    keys = references_for(tr)
    assert "chouhy_solotar" in keys                        # REGISTRY key implied by the CS Dispatch
    pairs = resolve_references(keys)
    assert any(bibtex_key == "ChouhySolotar2015" for bibtex_key, _ in pairs)
    # .references stays the FROZEN family+engine union (== citations()); the CS Dispatch
    # feeds the References SECTION only, it must NOT overwrite the result's .references.
    assert table.references == A.citations()


def test_verbose_run_writes_html_with_references(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Force the HTML path so the citation is greppable (a compiled PDF may compress it);
    # PDF selection is covered separately in test_writer.py.
    from quiverlab.trace import writer as W
    monkeypatch.setattr(W, "have_latex", lambda: None)
    quiverlab.verbose = True
    try:
        A = truncated_polynomial(2, field=CC)
        A.hochschild_cohomology(2)   # no explicit trace -> verbose auto-writes a file
    finally:
        quiverlab.verbose = False
    out = tmp_path / "quiverlab_traces"
    assert out.is_dir()
    files = list(out.glob("HHc_*.html"))
    assert files, "no worked-steps file written"
    # the References section shows the bibtex id resolved from Plan 06's `bar` entry
    assert "Hochschild1945" in files[0].read_text()
