"""JSON worked-steps renderer -- the COMPLETE machine record of a trace (Plan 34).

Marco's contract: every worked-steps computation delivers PDF + HTML + JSON, all
cached/served. This module renders the JSON leg: the recorder's full event stream
serialized EXACTLY -- every step, every matrix entry exact -- deterministically and
schema-versioned, so a downstream tool re-reads precisely what the engine computed.

The calling convention MIRRORS render_text/render_html exactly:
``render_json(events, title="", references=(), algebra=None) -> str`` (whatever
``writer.py`` and the GUI/hpc callers pass). The three human renderers DERIVE a
resulting-dimensions / projectives-injectives view from the events + ``algebra``;
the JSON renderer instead serializes the events THEMSELVES (the machine record), so
``algebra`` is accepted for signature parity but not embedded -- its repr already
rides the ``title``, and every derived view is reproducible from the event stream.
``references`` becomes the ``citations`` list (the same ``(bibtex_key, formatted)``
pairs the other renderers print).

Determinism: a pure function of ``(events, title, references)``. Matrices are
emitted through ``modules.qpa_module._json_entry`` (the exact ``int`` / ``"p/q"`` /
display-string grammar the no-code module input round-trips -- imported, never
duplicated), and a matrix carrying any NON-re-enterable entry (e.g. a GF(p^n)
element) flags its event ``display_only`` (mirroring ``module_blocks``). Dict keys
are sorted (``json.dumps(sort_keys=True)``); the event ORDER -- the only semantic
ordering -- is preserved as a list. No timestamps, no absolute paths, no floats:
two renders of identical input are byte-equal.
"""
import dataclasses
import json

import quiverlab
from quiverlab.modules.qpa_module import _json_entry, _reenterable

# The trace-JSON envelope version. Bump ONLY on an incompatible envelope change; a
# reader keys its parsing off this.
TRACE_SCHEMA = 1


def _matrix_cells(matrix):
    """Serialize a recorded matrix (``list[list[str]]`` of domain-element renderings)
    into exact JSON entries via ``_json_entry`` (int / "p/q" / display string), and
    report whether any entry is NON-re-enterable through the module input grammar (a
    GF(p^n) element, ...). Returns ``(cells, display_only)``. The cells are already
    strings (the recorder called ``str`` on each domain element), so ``dom=None``:
    integer/fraction strings normalise to exact ints / ``"p/q"``, anything else is
    kept verbatim and marks the event ``display_only`` -- never fabricated."""
    display_only = False
    out = []
    for row in matrix:
        out_row = []
        for cell in row:
            val = _json_entry(cell)              # cell is a string; dom=None
            if not _reenterable(val):
                display_only = True
            out_row.append(val)
        out.append(out_row)
    return out, display_only


def _safe(value):
    """Recursively convert an event field value to JSON-native, deterministic data.

    ``bool``/``int``/``str``/``None`` pass through unchanged -- booleans stay
    ``true``/``false`` (NOT the 0/1 that ``_json_entry`` would coerce, so the
    ``elided``/``heading``/``*_is_module`` flags read honestly); tuples and lists
    become arrays (CS chain words, ambiguity chains, resolution term lists); dicts
    become string-keyed objects with their keys stringified (the word-keyed
    linear-combination dicts of a ``ReductionStep`` included); any other leaf -- a
    domain-element coefficient -- goes through the exact ``_json_entry`` grammar."""
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [_safe(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _safe(v) for k, v in value.items()}
    return _json_entry(value)


def _event_dict(event):
    """One event as ``{"type": <ClassName>, <field>: <value>, ...}``.

    The ``matrix`` field (RankStep / ModuleDifferential / ExtDegree) is emitted
    through the exact entry grammar and may add ``display_only``; an
    elided-at-record event keeps ``matrix: null`` (``_safe(None)``) alongside
    ``elided: true`` and its shape/rank ``note`` -- shape+rank survive, the body is
    honestly absent, never fabricated."""
    if not dataclasses.is_dataclass(event):
        # Defensive: an unexpected non-dataclass event is recorded honestly by its
        # repr rather than silently dropped from the machine record.
        return {"type": type(event).__name__, "repr": repr(event)}
    out = {"type": type(event).__name__}
    display_only = False
    for f in dataclasses.fields(event):
        val = getattr(event, f.name)
        if f.name == "matrix" and val is not None:
            out["matrix"], display_only = _matrix_cells(val)
        else:
            out[f.name] = _safe(val)
    if display_only:
        out["display_only"] = True
    return out


def render_json(events, title="", references=(), algebra=None):
    """The trace as a schema-versioned, deterministic JSON string -- the complete
    machine record of the event stream. See the module docstring for the envelope
    and determinism contract."""
    events = list(events)
    from quiverlab.errors import QuiverlabError
    from quiverlab.trace.events import ALL_EVENTS
    for e in events:
        if not isinstance(e, ALL_EVENTS):
            raise QuiverlabError(
                "render_json received a non-event object of type %r -- likely an "
                "unpacked (events, result) tuple from a trace_* helper; the machine "
                "record must serialize only trace events" % type(e).__name__)
    envelope = {
        "quiverlab_trace_schema": TRACE_SCHEMA,
        "library_version": getattr(quiverlab, "__version__", "unknown"),
        "title": title,
        "citations": [{"key": k, "formatted": v} for k, v in references],
        "events": [_event_dict(e) for e in events],
    }
    return json.dumps(envelope, sort_keys=True, ensure_ascii=False, indent=2) + "\n"
