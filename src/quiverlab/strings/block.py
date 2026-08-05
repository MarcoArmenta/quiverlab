"""The algebra-only ``strings`` compute block (Plan 46 / C5).

Mirrors ``invariants/recognizers.py::recognizers_block``: returns a plain dict, shared
by both runners (server + Pyodide) so their blocks are byte-identical; each runner
stamps ``block["citations"]`` from ``block["references"]``.

Content: the special-biserial / string / gentle recognizer verdicts, a string census
up to a default length, band presence, an honest rep-type verdict for STRING algebras
(rep-finite iff no bands), and the AG invariant when GENTLE. A non-string algebra
reports null-ish strings/bands/rep_type plus a ``note``; the recognizer verdicts still
populate. Every field is guarded so one failure never crashes the block."""
from quiverlab.invariants.recognizers import (is_gentle, is_special_biserial,
                                              is_string)

_MAX_LENGTH = 8
_SAMPLE = 12


def _walk_repr(walk):
    """JSON-friendly string for a walk: ``e_v`` for a trivial walk, else the letters
    (``a`` direct, ``a^-1`` inverse) separated by spaces."""
    if len(walk) == 1 and walk[0][0] is None:
        return f"e_{walk[0][1]}"
    return " ".join(nm if d > 0 else f"{nm}^-1" for nm, d in walk)


def _guard(fn):
    try:
        return bool(fn())
    except Exception as exc:                       # honest per-flag error, never a 500
        return {"error": str(exc)}


def strings_block(A):
    """The ``strings`` no-code compute block: recognizer verdicts + string census +
    band presence + rep-type + (gentle) AG invariant. SHARED by both runners."""
    from quiverlab.strings.ag import ag_invariant
    from quiverlab.strings.walks import enumerate_strings, find_bands

    recognizers = {
        "is_special_biserial": _guard(lambda: is_special_biserial(A)),
        "is_string": _guard(lambda: is_string(A)),
        "is_gentle": _guard(lambda: is_gentle(A)),
    }

    strings = None
    bands = None
    rep_type = "unknown"
    ag = None
    note = None

    string_ok = recognizers["is_string"] is True
    if string_ok:
        try:
            cen = enumerate_strings(A, max_length=_MAX_LENGTH)
            strings = {
                "count": cen.count,
                "status": cen.status,
                "max_length": cen.max_length,
                "sample": [_walk_repr(w) for w in cen.walks[:_SAMPLE]],
            }
            band_walks = find_bands(A, max_length=_MAX_LENGTH)
            bands = {
                "exist": bool(band_walks),
                "sample": [_walk_repr(b) for b in band_walks[:_SAMPLE]],
            }
            # rep-finite iff no bands AND the census closed (honest contract).
            if band_walks:
                rep_type = "infinite"
            elif cen.status == "complete":
                rep_type = "finite"
            else:
                rep_type = "unknown"
        except Exception as exc:
            note = f"string census unavailable: {exc}"
    else:
        note = ("not a string algebra: strings/bands/rep-type are defined for "
                "special-biserial monomial kQ/I")

    if recognizers["is_gentle"] is True:
        try:
            inv = ag_invariant(A)
            ag = [[n, m] for (n, m) in inv.pairs]
        except Exception as exc:
            ag = {"error": str(exc)}

    block = {
        "recognizers": recognizers,
        "strings": strings,
        "bands": bands,
        "rep_type": rep_type,
        "ag_invariant": ag,
        "references": ["butler_ringel", "avella_geiss", "assem_book"],
    }
    if note is not None:
        block["note"] = note
    return block
