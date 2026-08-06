"""Plan 50 Task A -- v0.2.0 cross-tier integration reachability audit.

RELEASE / INFRASTRUCTURE gate. UNMARKED (no oracle-class markers) per the Plan-32
ruling for release/infra gates; fast bucket (tests/release/).

The per-plan batteries verify each new compute kind in isolation. This sweep
verifies the RELEASE-level contract across ALL v0.2.0 kinds at once: every kind is
(a) dispatchable in the wheel runner ``quiverlab.hpc.spec`` AND its Pyodide twin
``docs/gui/runner.py``, (b) in the GUI pick-list (both ``gui.js`` copies
byte-identical), (c) rendered by ``gui.js`` and given a heading in
``quiverlab.trace.results_html``, (d) i18n-complete (every referenced ``inv.*``
key exists in both catalogs, en/es key sets equal), (e) ETA-estimated (no kind is
silently "instant"), and (f) runnable end-to-end on ONE small algebra with the
Pyodide twin agreeing on the math.

Data-driven: the ONE canonical kinds list lives at the top; every check
introspects the real dispatch tables / pick-lists / i18n dicts / source, so there
is no hand-maintained duplicate list to drift.

SCOPE NOTE. "Compute kinds" here are the entries of a request's ``compute`` list.
The P44 (constructions) and P48 (marked-surface) presets the metaplan also
mentions are algebra-CONSTRUCTION inputs (they build an ``algebra`` block -- the
webapp algebra ``kind`` is only ``family``/``quiver``, see
``spec._build_algebra``), NOT compute kinds, so they are audited by the family /
surface batteries, not by this compute-kind reachability sweep.
"""
import importlib.util
import inspect
import json
import pathlib
import re
import tempfile

import pytest

from quiverlab.hpc import spec
from quiverlab.trace import results_html
from webapp.server.i18n import catalog
from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest

ROOT = pathlib.Path(__file__).resolve().parents[2]

# --------------------------------------------------------------------------- #
# THE canonical v0.2.0 compute-kind list (verified against the hpc/spec.py
# dispatch by test_canonical_list_is_faithful below).
# --------------------------------------------------------------------------- #
# Algebra-level kinds act on the algebra block (no ``module`` request block).
V020_ALGEBRA_KINDS = (
    "ext_algebra", "recognizers", "homological_profile", "ss_hochschild",
    "derived_fingerprint", "strings", "quasi_hereditary", "tau_tilting",
)
# Module-level kinds require a ``module`` block in the request.
V020_MODULE_KINDS = ("orbit_geometry", "tilting_check", "almost_split")
V020_KINDS = V020_ALGEBRA_KINDS + V020_MODULE_KINDS


# --------------------------------------------------------------------------- #
# Introspection helpers -- read the REAL surfaces, never a second hand list.
# --------------------------------------------------------------------------- #
def _gui_js() -> str:
    return (ROOT / "docs" / "gui" / "gui.js").read_text(encoding="utf-8")


_PYO = None


def _pyodide_runner():
    """Import docs/gui/runner.py (the Pyodide twin) as a module, once."""
    global _PYO
    if _PYO is None:
        path = ROOT / "docs" / "gui" / "runner.py"
        s = importlib.util.spec_from_file_location("_v020_gui_runner_twin", path)
        m = importlib.util.module_from_spec(s)
        s.loader.exec_module(m)
        _PYO = m
    return _PYO


def test_canonical_list_is_wellformed():
    """The single canonical list has no duplicates and algebra/module are disjoint
    -- guards the one hand-maintained artifact this file carries."""
    assert len(set(V020_KINDS)) == len(V020_KINDS), "duplicate in the canonical list"
    assert not (set(V020_ALGEBRA_KINDS) & set(V020_MODULE_KINDS)), \
        "a kind is listed as BOTH algebra- and module-level"


# --------------------------------------------------------------------------- #
# (a) dispatchable in the wheel runner AND its Pyodide twin
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("kind", V020_MODULE_KINDS)
def test_module_kind_in_dispatch_tables(kind):
    # Module kinds are gated by REAL frozenset tables on both tiers -- assert
    # membership directly (no regex).
    assert kind in spec.MODULE_KINDS, \
        f"{kind} missing from quiverlab.hpc.spec.MODULE_KINDS"
    pyo = _pyodide_runner()
    assert kind in pyo._MODULE_KINDS, \
        f"{kind} missing from docs/gui/runner.py _MODULE_KINDS"


@pytest.mark.parametrize("kind", V020_ALGEBRA_KINDS)
def test_algebra_kind_dispatched_in_both_runners(kind):
    # Algebra kinds have no frozenset table -- the dispatch is a ``if kind == ...``
    # branch chain, so scan the ACTUAL dispatch function bodies for the literal.
    server_src = inspect.getsource(spec._dispatch)
    pyo_src = inspect.getsource(_pyodide_runner().compute_one)
    assert ('"%s"' % kind) in server_src, \
        f"{kind} not dispatched in quiverlab.hpc.spec._dispatch"
    assert ('"%s"' % kind) in pyo_src, \
        f"{kind} not dispatched in docs/gui/runner.py compute_one"


def test_canonical_list_is_faithful():
    """Every canonical module kind is a real spec.MODULE_KINDS member, and every
    canonical algebra kind is a real branch -- so the hand-maintained list cannot
    name a kind the dispatch does not actually serve."""
    assert set(V020_MODULE_KINDS) <= set(spec.MODULE_KINDS)
    server_src = inspect.getsource(spec._dispatch)
    for kind in V020_ALGEBRA_KINDS:
        assert ('"%s"' % kind) in server_src


# --------------------------------------------------------------------------- #
# (b) GUI pick-list (both gui.js copies byte-identical)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("kind", V020_KINDS)
def test_kind_in_gui_picklist(kind):
    assert ('id="qlgui-%s"' % kind) in _gui_js(), \
        f"{kind} has no checkbox (id=qlgui-{kind}) in docs/gui/gui.js pick-list"


def test_gui_js_copies_byte_identical():
    a = (ROOT / "docs" / "gui" / "gui.js").read_text(encoding="utf-8")
    b = (ROOT / "webapp" / "static" / "gui" / "gui.js").read_text(encoding="utf-8")
    assert a == b, "docs/gui/gui.js and webapp/static/gui/gui.js must be byte-identical"


# --------------------------------------------------------------------------- #
# (c) rendered by gui.js AND given a heading in trace/results_html.py
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("kind", V020_KINDS)
def test_kind_rendered_in_gui_js(kind):
    assert ('name === "%s"' % kind) in _gui_js(), \
        f"{kind} has no renderBlock branch (name === \"{kind}\") in docs/gui/gui.js"


@pytest.mark.parametrize("kind", V020_KINDS)
def test_kind_has_report_heading(kind):
    assert kind in results_html._HEADINGS, \
        f"{kind} missing from quiverlab.trace.results_html._HEADINGS"


# --------------------------------------------------------------------------- #
# (d) i18n: referenced inv.* keys exist in BOTH catalogs; en/es key sets equal
# --------------------------------------------------------------------------- #
def _referenced_inv_keys():
    """Every ``inv.<name>`` i18n key referenced as a QUOTED string literal in the
    webapp templates or static JS (the ``t("inv.x")`` calls). The quote-anchored
    pattern skips JS member access like ``inv.split(...)``."""
    pat = re.compile(r"""["']inv\.([a-z_]+)["']""")
    keys = set()
    for sub in ("templates", "static"):
        base = ROOT / "webapp" / sub
        if not base.exists():
            continue
        for f in list(base.rglob("*.html")) + list(base.rglob("*.js")):
            for m in pat.finditer(f.read_text(encoding="utf-8")):
                keys.add("inv." + m.group(1))
    return keys


def test_referenced_inv_keys_exist_in_both_catalogs():
    en, es = catalog("en"), catalog("es")
    refs = _referenced_inv_keys()
    assert refs, "no inv.* i18n references found -- the collector regex went stale"
    missing_en = sorted(k for k in refs if k not in en)
    missing_es = sorted(k for k in refs if k not in es)
    assert not missing_en, f"referenced inv.* keys absent from en.json: {missing_en}"
    assert not missing_es, f"referenced inv.* keys absent from es.json: {missing_es}"


def test_en_es_key_sets_equal():
    # Release-level restatement of tests/webapp/test_i18n.test_en_es_key_parity, so
    # this audit is self-contained: a key added to one catalog and not the other
    # would leave a page half-translated.
    en, es = set(catalog("en")), set(catalog("es"))
    assert en == es, "en/es i18n key mismatch: " + str(sorted(en ^ es))


# --------------------------------------------------------------------------- #
# (e) ETA scalars: no v0.2.0 kind is silently "instant"
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("kind", V020_KINDS)
def test_kind_has_eta_estimate(kind):
    # The per-kind ETA scalars table lives ONLY on the Pyodide side
    # (docs/gui/runner.py::ETA_MODEL["scalars"]); the server estimator
    # (webapp/server/estimator.py) sizes ops-based, with no per-kind table (see
    # test_server_estimator_has_no_rival_scalar_table). A v0.2.0 kind missing from
    # the scalars dict falls through to the 0.1 default in ``_units_for`` and is
    # silently under-estimated as "instant" -- exactly the gap this audit caught
    # for tilting_check / orbit_geometry (fixed in the same commit).
    scalars = _pyodide_runner().ETA_MODEL["scalars"]
    assert kind in scalars, (
        f"{kind} missing from docs/gui/runner.py ETA_MODEL['scalars'] -- "
        f"it would be silently estimated at the 0.1 'instant' default")


def test_server_estimator_has_no_rival_scalar_table():
    # Documents WHY (e) checks only the Pyodide table: the server estimator sizes
    # ops-based, with NO per-kind scalar dict. If one is ever added it MUST agree
    # key-for-key with the Pyodide table -- promote (e) to a two-sided parity check
    # then (the P45 merge nearly let the two diverge).
    import webapp.server.estimator as est
    assert not hasattr(est, "ETA_MODEL"), \
        "server estimator grew an ETA_MODEL -- make (e) a two-sided parity check"


# --------------------------------------------------------------------------- #
# (f) end-to-end smoke: ONE small algebra runs EVERY v0.2.0 kind through
#     webapp.server.runner.run_spec, and the Pyodide twin agrees on the math.
# --------------------------------------------------------------------------- #
_KA3_GF7 = {"kind": "quiver", "vertices": [1, 2, 3],
            "arrows": {"a": [1, 2], "b": [2, 3]}, "relations": [],
            "field": {"kind": "GF", "p": 7, "n": 1}}
_SIMPLE_S2 = {"builtin": {"kind": "simple", "vertex": 2}}   # neither proj nor inj in kA3

# Non-math provenance keys: both tiers resolve citations from the SAME references,
# but the server keeps the raw ``references`` list on some blocks and the Pyodide
# twin does not (e.g. tilting_check) -- neither is mathematics, so both are dropped
# before the cross-tier equality check.
_PROVENANCE = ("references", "citations")

# Compute item per algebra kind (range/budget kinds get a small argument).
_ALG_ITEMS = {
    "ext_algebra": "ext_algebra:0..3",
    "recognizers": "recognizers",
    "homological_profile": "homological_profile",
    "ss_hochschild": "ss_hochschild:0..4",
    "derived_fingerprint": "derived_fingerprint",
    "strings": "strings",
    "quasi_hereditary": "quasi_hereditary",
    "tau_tilting": "tau_tilting:256",
}

# ALLOWED honest-refusal on kA3/GF(7). ss_hochschild builds the exponential (b,B)
# bar bicomplex; on kA3 its degree-5 boundary is 70,312,500 cells (> max_cells =
# 4,000,000), so the shared specseq_block returns the library's documented
# honest-refusal block (a *string* ``error`` that keeps kind/top), NEVER a crash.
# This is the deterministic max_cells guard, exercised here on purpose -- the block
# is still byte-identical across both tiers. Every OTHER v0.2.0 kind must compute
# cleanly on this algebra.
_MAY_HONEST_REFUSE = {"ss_hochschild"}


def _math_canon(block):
    """The block's MATH content as a canonical JSON string. Provenance keys are
    dropped; the whole thing is round-tripped through JSON so the server result
    (raw Python dict -- int dict-keys such as ``brick_dimvec`` / dimension vectors)
    and the Pyodide result (already crossed JSON, so those keys are strings)
    compare EQUAL when the mathematics agrees. This is the delegation/twin-test
    idiom (json.dumps(sort_keys=True))."""
    math = {k: v for k, v in block.items() if k not in _PROVENANCE}
    return json.dumps(json.loads(json.dumps(math, default=str)), sort_keys=True)


def _server_result(body):
    with tempfile.TemporaryDirectory() as d:
        return run_spec(ComputeRequest.model_validate(body), d)


def _pyodide_result(body, item):
    pyo = _pyodide_runner()
    assert json.loads(pyo.run_build(json.dumps(body)))["ok"], "Pyodide run_build failed"
    return json.loads(pyo.compute_one(item))


def _assert_block_reachable(kind, block, may_refuse):
    assert isinstance(block, dict) and block, f"{kind} block empty / non-dict"
    # A hard crash surfaces as the runner's {'error': {'type': ...}} wrapper (a
    # DICT error). That is never acceptable.
    assert not isinstance(block.get("error"), dict), \
        f"{kind} hard-crashed: {block.get('error')}"
    # A STRING error is the documented honest-refusal shape -- allowed only for the
    # known-heavy kinds on this algebra.
    if isinstance(block.get("error"), str):
        assert kind in may_refuse, \
            f"{kind} unexpectedly refused on kA3/GF(7): {block['error']}"


@pytest.mark.parametrize("kind", V020_ALGEBRA_KINDS)
def test_smoke_algebra_kind_cross_tier(kind):
    item = _ALG_ITEMS[kind]
    body = {"schema": 2, "algebra": _KA3_GF7, "compute": [item]}
    res = _server_result(body)
    assert "error" not in res, f"top-level failure for {kind}: {res.get('error')}"
    assert kind in res["results"], f"no result block for {kind}"
    block = res["results"][kind]
    _assert_block_reachable(kind, block, _MAY_HONEST_REFUSE)
    # cross-tier: the Pyodide twin produces the SAME mathematics (provenance aside).
    p = _pyodide_result(body, item)
    assert p["ok"], p
    assert _math_canon(block) == _math_canon(p["block"]), \
        f"{kind} server/Pyodide math disagree"


@pytest.mark.parametrize("kind", V020_MODULE_KINDS)
def test_smoke_module_kind_cross_tier(kind):
    body = {"schema": 2, "algebra": _KA3_GF7, "module": _SIMPLE_S2, "compute": [kind]}
    res = _server_result(body)
    assert "error" not in res, f"top-level failure for {kind}: {res.get('error')}"
    assert kind in res["results"], f"no result block for {kind}"
    block = res["results"][kind]
    # S_2 is neither projective nor injective in kA3, so all three module kinds
    # (almost-split, tilting test, orbit geometry) compute cleanly -- no refusal
    # is allowed here.
    _assert_block_reachable(kind, block, may_refuse=set())
    assert "error" not in block, f"{kind} refused on S_2/kA3/GF(7): {block.get('error')}"
    p = _pyodide_result(body, kind)
    assert p["ok"], p
    assert _math_canon(block) == _math_canon(p["block"]), \
        f"{kind} server/Pyodide math disagree"
