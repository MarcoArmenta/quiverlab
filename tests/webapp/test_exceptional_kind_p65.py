"""The ``exceptional_sequences`` algebra-only compute kind (Plan 65 / R27+R28): served by
hpc.spec (the wheel), the shared block dispatches BOTH halves -- classical hereditary
(braid-orbit counts + Dynkin closed form) and tau-exceptional (n!*#sTt) -- each honest
about applicability, with resolved citations. Unmarked (extras-gated dir)."""
import tempfile

from quiverlab.hpc.spec import parse_request, run


def _linear(n, field=None):
    field = field or {"kind": "GF", "p": 7, "n": 1}
    return {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": list(range(1, n + 1)),
                        "arrows": {f"a{i}": [i, i + 1] for i in range(1, n)},
                        "relations": [], "field": field},
            "compute": ["exceptional_sequences:512"],
            "artifacts": {"pdf": False, "tikz": False}}


def _nonhered():
    # k(1<->2)/rad^2 : all length-2 paths zero -- non-hereditary, tau-tilting-finite.
    return {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1, 2],
                        "arrows": {"a": [1, 2], "b": [2, 1]},
                        "relations": ["a*b", "b*a"], "field": {"kind": "GF", "p": 7, "n": 1}},
            "compute": ["exceptional_sequences:512"],
            "artifacts": {"pdf": False, "tikz": False}}


def test_exceptional_block_shape_hereditary():
    out = run(parse_request(_linear(3)), tempfile.mkdtemp())
    b = out["results"]["exceptional_sequences"]
    assert b["kind"] == "exceptional_sequences" and b["n"] == 3
    assert b["hereditary"] is True and b["dynkin_type"] == "A_3"
    assert b["classical"]["count"] == 16
    assert b["classical"]["closed_form_count"] == 16          # n! h^n / |W|
    assert b["classical"]["transitive"] is True
    assert b["tau"]["signed_count"] == 84                     # 3! * 14
    assert b["tau"]["stt_count"] == 14
    assert "buan_marsh_tau_exceptional" in b["references"]
    assert "crawley_boevey_exceptional" in b["references"]
    keys = [k for k, _ in b["citations"]]                     # resolved bibtex keys
    assert "BuanMarsh2021" in keys and "Obaid2013" in keys


def test_exceptional_block_nonhereditary_tau_only():
    out = run(parse_request(_nonhered()), tempfile.mkdtemp())
    b = out["results"]["exceptional_sequences"]
    assert b["hereditary"] is False
    assert b["classical"] is None                             # braid-orbit surface N/A
    assert b["tau"]["signed_count"] == 12                     # 2! * 6
    assert b["tau"]["stt_count"] == 6


def test_exceptional_rep_infinite_classical_refused():
    # 2-Kronecker: hereditary but rep-INFINITE -> classical enumeration refuses (infinite
    # braid orbit), tau-exceptional refuses (tau-tilting-infinite). Both honest, never a count.
    req = {"schema": 1,
           "algebra": {"kind": "quiver", "vertices": [1, 2],
                       "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
                       "field": {"kind": "GF", "p": 7, "n": 1}},
           "compute": ["exceptional_sequences:200"],
           "artifacts": {"pdf": False, "tikz": False}}
    out = run(parse_request(req), tempfile.mkdtemp())
    b = out["results"]["exceptional_sequences"]
    assert b["hereditary"] is True
    assert b["classical"]["complete"] is False                # infinite braid orbit
    assert b["tau"] is None                                   # tau-tilting-infinite
