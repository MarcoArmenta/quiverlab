"""Search-first landing (Marco 2026-08-06): a search bar on /draw whose matches
load a small, pre-validated example into the canvas and auto-run it.

Contract under test:
  * the embedded ``SEARCH_INDEX`` in the vendored gui.js is a well-formed set of
    environments, each with a runnable example (quiver + compute [+ module]);
  * the union of example kinds covers the whole v0.2.0 surface, not just HH;
  * every example is a VALID server request (``quiverlab.hpc.spec.parse_request``),
    and the smallest one runs end-to-end with no per-block error -- the same
    dispatch the GUI's worker calls, so a loaded example cannot 4xx/500;
  * ``/draw`` under EVERY configured language prefix carries the four
    ``data-search-*`` i18n attributes gui.js reads;
  * gui.js / gui.css stay byte-identical across their two homes.
"""
import json
import re
import tempfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from quiverlab.hpc import spec
from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.i18n import LANGS

ROOT = Path(__file__).resolve().parents[2]
GUI_JS = ROOT / "webapp" / "static" / "gui" / "gui.js"

# The kinds a newcomer should be able to reach by search -- the whole surface,
# deliberately spanning Hochschild / derived / modules / AR / geometry / gentle.
REQUIRED_KINDS = {
    "hh_cohomology", "ss_hochschild", "tau_tilting", "strings", "quasi_hereditary",
    "derived_fingerprint", "ext", "tau", "projective_resolution", "decompose",
    "orbit_geometry", "recognizers", "ext_algebra",
}


def _load_index():
    src = GUI_JS.read_text(encoding="utf-8")
    m = re.search(
        r"// QLGUI-SEARCH-INDEX-BEGIN\s*var SEARCH_INDEX\s*=\s*(\[.*?\])\s*;\s*"
        r"// QLGUI-SEARCH-INDEX-END",
        src, re.S)
    assert m, "SEARCH_INDEX sentinels/array not found in vendored gui.js"
    return json.loads(m.group(1))


def _to_request(example):
    algebra = {"kind": "quiver", "vertices": example["vertices"],
               "arrows": example["arrows"], "relations": example["relations"],
               "field": example["field"]}
    has_mod = any(k in example for k in ("module", "ext_target", "tor_target"))
    req = {"schema": 2 if has_mod else 1, "algebra": algebra,
           "compute": example["compute"],
           "artifacts": {"pdf": False, "tikz": False}}
    for k in ("module", "ext_target", "tor_target"):
        if k in example:
            req[k] = example[k]
    return req


INDEX = _load_index()


def test_index_shape():
    assert len(INDEX) >= 15, "search index should cover the whole surface"
    ids = set()
    for e in INDEX:
        assert e["id"] and e["id"] not in ids, ("id missing or duplicated", e)
        ids.add(e["id"])
        assert e["title"].strip(), e
        assert e["category"].strip(), e
        assert e["keywords"], e
        ex = e["example"]
        for field in ("vertices", "arrows", "field", "compute"):
            assert ex.get(field) is not None, (e["id"], field)
        assert ex["vertices"], e["id"]
        assert ex["compute"], e["id"]


def test_index_covers_the_surface():
    kinds = set()
    for e in INDEX:
        for item in e["example"]["compute"]:
            kinds.add(item.split(":")[0])
    missing = REQUIRED_KINDS - kinds
    assert not missing, f"search index does not reach: {sorted(missing)}"


def test_keywords_carry_multilingual_synonyms():
    # At least some CJK (Chinese) synonyms must be present so search works in the
    # zh UI even though the labels are English.
    blob = "".join(k for e in INDEX for k in e["keywords"])
    assert re.search(r"[一-鿿]", blob), "no Chinese keywords in the index"


@pytest.mark.parametrize("entry", INDEX, ids=[e["id"] for e in INDEX])
def test_every_example_is_a_valid_request(entry):
    # Validation only (fast): the same parse the webapp runs before dispatch.
    spec.parse_request(_to_request(entry["example"]))


def test_smallest_example_runs_without_error():
    smallest = min(INDEX, key=lambda e: len(e["example"]["compute"]))
    req = _to_request(smallest["example"])
    with tempfile.TemporaryDirectory() as d:
        res = spec.run(req, d, write_result=False)
    errs = {k: v["error"] for k, v in res["results"].items()
            if isinstance(v, dict) and v.get("error")}
    assert not errs, (smallest["id"], errs)


def test_draw_page_carries_search_attributes_in_every_language(tmp_path):
    c = TestClient(create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})))
    for lang in LANGS:
        prefix = "" if lang == "en" else "/" + lang
        html = c.get(f"{prefix}/draw").text
        for attr in ("data-search-placeholder", "data-search-related",
                     "data-search-none", "data-search-load"):
            assert attr in html, (lang, attr)


def test_gui_twins_byte_identical():
    for name in ("gui.js", "gui.css"):
        docs = (ROOT / "docs" / "gui" / name).read_bytes()
        vendored = (ROOT / "webapp" / "static" / "gui" / name).read_bytes()
        assert docs == vendored, f"{name} drifted between docs/gui and webapp/static/gui"
