"""P80 -- the three-tier compute-grammar parity GATE.

quiverlab parses compute strings in THREE places: ``quiverlab.hpc.spec`` (the CLI and the
worker), ``webapp.server.schema`` (the deployed service's request validation) and
``docs/gui/runner.py`` (the Pyodide twin in the browser). They must agree about what a
valid request IS and what it MEANS -- otherwise the same string is accepted by one tier
and rejected by another, which is what the P80 audit found for ``skew_gentle:<budget>``
(accepted by the CLI and the twin, rejected by the webapp schema).

This turns that one-off audit into a standing gate, so the drift cannot recur silently.
Unmarked (extras-gated dir, the Plan-32 ruling).
"""
import importlib.util
import pathlib
import re

import pytest

from quiverlab.hpc.spec import SpecError
from quiverlab.hpc.spec import parse_compute_item as hpc_parse
from webapp.server.schema import SchemaError
from webapp.server.schema import parse_compute_item as web_parse

ROOT = pathlib.Path(__file__).resolve().parents[2]
_SPEC_SRC = (ROOT / "src" / "quiverlab" / "hpc" / "spec.py").read_text(encoding="utf-8")

# Every kind whose compute string carries a BUDGET (or a radius,budget pair) instead of a
# 'name:0..N' degree range. Derived from the canonical grammar site so a kind added later
# is covered automatically.
_SPECIAL = sorted(set(re.findall(r's == "([a-z_0-9]+)" or s\.startswith\(', _SPEC_SRC)))


def _twin():
    spec = importlib.util.spec_from_file_location(
        "gui_runner_p80", ROOT / "docs" / "gui" / "runner.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_derived_special_form_list_is_not_vacuous():
    # If the regex ever stops matching, every parity assertion below would pass over an
    # EMPTY list -- a green suite proving nothing. Anchor it on known members.
    assert len(_SPECIAL) >= 12, _SPECIAL
    for known in ("tau_tilting", "wall_chamber", "skew_gentle", "cluster_category"):
        assert known in _SPECIAL, (known, _SPECIAL)


@pytest.mark.parametrize("kind", _SPECIAL)
def test_the_bare_form_parses_identically_in_all_three_sites(kind):
    a, b = hpc_parse(kind), web_parse(kind)
    assert a.kind == b.kind == kind
    assert (a.lo, a.hi) == (b.lo, b.hi)
    assert _twin()._parse_compute(kind)[0] == kind


@pytest.mark.parametrize("kind", _SPECIAL)
def test_the_budget_form_parses_identically_in_all_three_sites(kind):
    # `silting` is the one two-argument form (radius,budget); everything else takes a
    # single integer, so it is exercised with its own grammar.
    spec = f"{kind}:3,64" if kind == "silting" else f"{kind}:64"
    a, b = hpc_parse(spec), web_parse(spec)
    assert a.kind == b.kind == kind
    assert (a.lo, a.hi) == (b.lo, b.hi), (kind, (a.lo, a.hi), (b.lo, b.hi))
    twin_kind, twin_top = _twin()._parse_compute(spec)
    assert twin_kind == kind
    expected_top = (a.lo, a.hi) if kind == "silting" else a.hi
    assert twin_top == expected_top, (kind, twin_top, expected_top)


@pytest.mark.parametrize("kind", _SPECIAL)
def test_a_malformed_budget_is_refused_by_both_server_grammars(kind):
    spec = f"{kind}:abc"
    with pytest.raises(SpecError):
        hpc_parse(spec)
    with pytest.raises(SchemaError):
        web_parse(spec)
