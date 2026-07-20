"""zoo(dim_max): the curated exact zoo (spec §3.4), lifted from hanlab open_zoo.
Each record is a confluent reduction system over k<x,y[,z]>; rebuilt through the
Plan-03/04 reduction-system -> Algebra path."""
import json
import pathlib

from quiverlab.combinat.quiver import Quiver

_CATALOG = pathlib.Path(__file__).with_name("zoo_catalog.json")
_GENS = ("x", "y", "z")


def load_catalog():
    data = json.loads(_CATALOG.read_text(encoding="utf-8"))
    return [r for r in data if isinstance(r, dict) and "rules" in r]


def _word(idxs):
    return "*".join(_GENS[i] for i in idxs)


def _render_terms(terms):
    """Render a signed sum of (int coeff, generator-index word) as a flat
    relation string with folded signs and no parentheses.

    The reduction-system records encode a rule ``lead = sum c*w``; we hand the
    parser ``lead - sum c*w = 0``. The Plan-03 relation splitter (``_split_terms``)
    does not track parenthesis depth, so a nested/parenthesised coefficient such
    as ``((1)`` or ``(-1)`` breaks parsing -- signs MUST be folded into the joining
    operator here (see the module test-suite and the Task-11 report).
    """
    out = ""
    for c, w in terms:
        if c == 0:
            continue
        mag = abs(c)
        piece = _word(w) if mag == 1 else f"{mag}*{_word(w)}"
        if not out:
            out = f"-{piece}" if c < 0 else piece
        else:
            out += (" - " if c < 0 else " + ") + piece
    return out


def build_from_record(rec, field=None):
    ngen = rec["ngen"]
    gens = _GENS[:ngen]
    Q = Quiver([1], {g: (1, 1) for g in gens})
    rels = []
    for lead, tail in rec["rules"]:
        if not tail:
            rels.append(_word(lead))                            # lead -> 0 (monomial)
        else:
            # lead = sum c*w  <=>  lead - sum c*w = 0
            rels.append(_render_terms([(1, lead)] + [(-c, w) for c, w in tail]))
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("han_conjecture", "chouhy_solotar")
    A.zoo_name = rec["name"]
    return A


def zoo(dim_max=12, field=None):
    recs = [r for r in load_catalog() if r["dim"] <= dim_max]
    recs.sort(key=lambda r: (r["dim"], r["name"]))
    for rec in recs:
        yield build_from_record(rec, field=field)
