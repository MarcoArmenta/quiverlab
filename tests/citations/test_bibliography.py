"""quiverlab.bibliography(): iterable, grouped-by-topic, annotated (web /literature)."""
from quiverlab import bibliography, citations


def test_full_bibliography_groups_by_topic_and_annotates():
    b = bibliography()
    assert set(b.groups) <= {"Algorithms", "Families", "Finite fields", "Foundations"}
    assert "bardzell" in b.keys and "nakayama" in b.keys
    s = str(b)
    assert "Bardzell" in s                                        # formatted citation
    assert "minimal projective bimodule resolution" in s.lower()  # annotation present


def test_iterates_with_plan09_entry_view_fields():
    """The frozen contract Plan 09's references.entry_view() consumes: each entry
    exposes key/formatted/doi/arxiv/topic/annotation as attributes."""
    entries = list(bibliography())                               # ITERABLE
    assert entries
    by_key = {e.key: e for e in entries}
    for attr in ("key", "formatted", "doi", "arxiv", "topic", "annotation"):
        assert hasattr(by_key["bardzell"], attr)
    cs = by_key["chouhy_solotar"]
    assert cs.arxiv == "1406.2300"                               # parsed from .bib note
    assert cs.doi == "10.1016/j.jalgebra.2015.02.019"
    assert cs.topic == "Algorithms"
    assert by_key["bardzell"].doi == "10.1006/jabr.1996.6813"


def test_subset_bibliography_and_bibtex_dedup():
    b = bibliography(keys=["quantum_ci", "qci_hh_oracle", "quantum_ci"])
    assert b.keys == ("quantum_ci", "qci_hh_oracle")             # order-preserving, deduped
    bt = b.bibtex()
    assert bt.count("@article{BGMS2005") == 1                    # bibtex_key dedup
    assert "@article{BerghErdmann2008" in bt


def test_to_dict_is_json_ready():
    import json
    d = bibliography(keys=["bardzell"]).to_dict()
    json.dumps(d)                                                # no exception
    assert d["groups"]["Algorithms"][0]["bibtex_key"] == "Bardzell1997"
    assert d["groups"]["Algorithms"][0]["annotation"]
