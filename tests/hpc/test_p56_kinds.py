"""The fundamental_group + simply_connected kinds through the spec runner (Plan 56).

``quiverlab.hpc.spec`` exposes ``parse_request`` + ``run``; this drives the real
``_dispatch`` branches. The webapp runner delegates to the same dispatch (pinned
byte-stable by ``fundamental_group_square_gf7`` / ``simply_connected_zito`` in
``tests/webapp/_runner_goldens.json``), and the Pyodide twin (``docs/gui/runner.py``)
is checked byte-identical in ``tests/webapp/test_coverings_exposure_p56.py``.
"""
from quiverlab.hpc.spec import parse_request
from quiverlab.hpc.spec import run as spec_run

_SQUARE_REL = {"kind": "quiver", "vertices": [1, 2, 3, 4],
               "arrows": {"a": [1, 2], "b": [2, 4], "c": [1, 3], "d": [3, 4]},
               "relations": ["a*b - c*d"], "field": {"kind": "GF", "p": 7, "n": 1}}
_SQUARE_FREE = {"kind": "quiver", "vertices": [1, 2, 3, 4],
                "arrows": {"a": [1, 2], "b": [2, 4], "c": [1, 3], "d": [3, 4]},
                "relations": [], "field": {"kind": "GF", "p": 7, "n": 1}}
_ZITO = {"kind": "quiver", "vertices": [1, 2, 3, 4, 5],
         "arrows": {"al": [1, 2], "be": [2, 3], "ga": [3, 5],
                    "de": [2, 4], "ep": [4, 5]},
         "relations": ["al*be*ga - al*de*ep"], "field": {"kind": "QQ"}}


def _req(algebra, compute):
    return parse_request({"schema": 1, "algebra": algebra, "compute": compute,
                          "artifacts": {"pdf": False, "tikz": False}})


def test_fundamental_group_commutative_square(tmp_path):
    b = spec_run(_req(_SQUARE_REL, ["fundamental_group"]),
                 tmp_path)["results"]["fundamental_group"]
    assert b["abelianization"] == {"free_rank": 0, "invariant_factors": []}
    assert b["components"] == 1
    assert b["hom_to_additive_dim"] == 0
    assert b["presentation_note"].startswith("presentation invariant")
    assert "intrinsic_note" in b and "Z x C_p" in b["intrinsic_note"]
    assert b["references"][0] == "assem_delapena" and b["citations"]


def test_fundamental_group_hereditary_square_is_Z(tmp_path):
    b = spec_run(_req(_SQUARE_FREE, ["fundamental_group"]),
                 tmp_path)["results"]["fundamental_group"]
    assert b["abelianization"] == {"free_rank": 1, "invariant_factors": []}
    assert b["latex"].endswith(r"\mathbb{Z}")            # pi1^ab = Z


def test_simply_connected_zito_true_but_not_strongly(tmp_path):
    b = spec_run(_req(_ZITO, ["simply_connected"]),
                 tmp_path)["results"]["simply_connected"]
    assert b["verdict"] is True                          # Zito: simply connected
    assert b["strongly"]["verdict"] is False             # but NOT strongly (P62 gate)
    assert b["strongly"]["witness"]["vertex"] == 2
    assert b["references"] == ["assem_book", "skowronski_ssc", "le_meur_pi1"]
    assert b["citations"]


def test_simply_connected_hereditary_square_false(tmp_path):
    b = spec_run(_req(_SQUARE_FREE, ["simply_connected"]),
                 tmp_path)["results"]["simply_connected"]
    assert b["verdict"] is False and b["witness"]["kind"] == "nontrivial_pi1ab"
