"""Worked-steps output-path contract + the printed one-liner (spec §3.8).

Renders the recorded worked steps to a self-contained, no-JavaScript HTML report
and its JSON machine record. Output:
./quiverlab_traces/HHc_<hash>.<ext> (cohomology) / HHh_<hash>.<ext> (homology)
(Plan 09 collects the newest *.html from this directory -- the glob is
extension-based, so the safe stem does not affect it).

The filename hash is 12 hex chars (48 bits): the plan's original 4-hex (16-bit)
stem would collide by the birthday bound at ~256 distinct traces and silently
overwrite; 12 hex pushes the collision horizon out of practical reach while staying
fully deterministic (no floats)."""
import hashlib
import pathlib

from quiverlab.trace.render_html import render_html
from quiverlab.trace.render_json import render_json

# Filesystem-safe filename stems for the caret-bearing kinds (no "^" in a filename).
_SAFE_STEM = {"HH^": "HHc", "HH_": "HHh"}


def _hash(algebra, kind, top):
    h = hashlib.sha1(("%s|%s|%s" % (repr(algebra), kind, top)).encode("utf-8"))
    return h.hexdigest()[:12]


def _superseding_dispatch_index(events):
    """Index of a SUPERSEDING ``Dispatch`` -- the auto->CS depth fallback
    (``core.algebra._cs_depth_fallback``) records a SECOND ``Dispatch`` (route
    "chouhy-solotar") AFTER the abandoned bar/fast route already filled the recorder with
    its partial per-degree worked steps. ``write_trace`` only ever renders a SINGLE HH
    computation, so a second Dispatch is unambiguously that fallback (the up-front CS route
    and every normal route emit exactly one). Returns the index of the LAST Dispatch when
    more than one is present, else ``None`` (every non-fallback stream)."""
    from quiverlab.trace.events import Dispatch
    idxs = [i for i, e in enumerate(events) if isinstance(e, Dispatch)]
    return idxs[-1] if len(idxs) > 1 else None


def _drop_superseded(events):
    """Correction 2 (auto->CS fallback report coherence): when the bar/fast route hit its
    depth wall and ``engine="auto"`` rerouted to Chouhy-Solotar, the recorder holds BOTH the
    ABANDONED partial bar worked steps AND the full CS worked steps in one stream. Rendering
    both is confusing in a homework document, and it lets the HH drift gate pass only by
    last-wins luck (``derive_dims``'s degree-keyed dicts happen to resolve to all-CS).

    Drop the abandoned worked-step events (everything that is NOT a ``Dispatch``) that
    precede the superseding CS dispatch, keeping the routing STORY -- both ``Dispatch`` lines,
    the second stating WHY it rerouted ("... exceeded max_cells ... reroutes to the
    Chouhy-Solotar resolution ...") -- and the CS worked steps as the single authoritative
    set. A no-op for every non-fallback stream (<=1 Dispatch), so normal reports (and the
    goldens) are byte-unchanged."""
    from quiverlab.trace.events import Dispatch
    cut = _superseding_dispatch_index(events)
    if cut is None:
        return list(events)
    return [e for i, e in enumerate(events) if i >= cut or isinstance(e, Dispatch)]


def _authoritative_result(events, table):
    """Cross-check the event-DERIVED HH dims against the engine's AUTHORITATIVE dims (the
    HH DRIFT GATE -- mirroring the module drift gates in ``trace.modules``) and append a
    ``ResultDims`` event carrying those authoritative dims + variance as the Result source.

    Why: the CS engine used to record only ResolutionTerms, so ``derive_dims`` re-derived
    C_n (the cochain SPACE dims) and mislabelled HH^ as HH_; the report shipped wrong
    numbers with no gate to catch it. Now the renderers show ``ResultDims`` (the engine's
    ``HHTable.dims``/``.kind``) as the Result, and this gate refuses to ship a report whose
    event-derived dims/variance disagree with what the engine returned.

    ``table`` is the HHTable an HH computation returned (``None`` for a module report --
    then nothing is injected and the module Ext/Tor Result stands). When the engine records
    no per-degree worked steps (the fast GF(p) engine emits only a Dispatch) there is
    nothing to cross-check; an honest note rides on the ResultDims so the report still has a
    Result line without fabricating steps.

    Correction 1 (elision vs the gate): when the MAX_EVENTS buffer ELIDED steps the surviving
    event list is TRUNCATED, so ``derive_dims`` reconstructs a short/partial dim list.
    Comparing that to the engine's full ``table.dims`` would HARD-RAISE and defeat the buffer
    cap, whose whole point is to let an over-long report still ship (with an elision note).
    So the derive-vs-table CROSS-CHECK is SKIPPED under elision (detected via the recorder's
    reserved BUFFER_FULL StepNote); the authoritative ResultDims is emitted either way."""
    if table is None:
        return events
    from quiverlab.errors import QuiverlabError
    from quiverlab.trace.events import (
        ResultDims, ResolutionTerm, RankStep, ModuleTerm, ModuleDifferential, ExtDegree,
        StepNote,
    )
    from quiverlab.trace.recorder import BUFFER_FULL_PREFIX
    from quiverlab.trace.render_text import derive_dims, _dims_kind
    elided = any(isinstance(e, StepNote) and BUFFER_FULL_PREFIX in e.text for e in events)
    derived = derive_dims(events)
    if derived and not elided:
        if list(derived) != list(table.dims):
            raise QuiverlabError(
                "HH worked-steps drift: the report's event-derived dimensions %s do not "
                "match the engine's returned dimensions %s -- refusing to ship a report "
                "that misstates the computation" % (list(derived), list(table.dims)))
        if _dims_kind(events) != table.kind:
            raise QuiverlabError(
                "HH worked-steps drift: the report's event-derived variance %r does not "
                "match the engine's %r" % (_dims_kind(events), table.kind))
    has_steps = any(isinstance(e, (ResolutionTerm, RankStep, ModuleTerm,
                                   ModuleDifferential, ExtDegree)) for e in events)
    # An honest "fast engine has no worked steps" note only when there truly are none AND
    # nothing was elided (under elision the reserved BUFFER_FULL StepNote already explains
    # the missing steps, so claiming the fast engine recorded none would be a lie).
    note = ("" if (has_steps or elided) else
            "The fast engine records no per-degree worked steps; the dimensions shown "
            "are the computed result.")
    return events + [ResultDims(kind=table.kind, dims=list(table.dims), note=note)]


def write_trace(events, table, algebra, kind, top, references=(), out_dir=None):
    """Render the worked steps to the self-contained, print-ready HTML report and the
    JSON machine record, print the one-liner, and return the produced HTML path (str).

    Both deliverables are pure functions of the (deduplicated) event stream, so they
    are byte-identical for identical input. The HTML report is print-ready (the GUI /
    a browser can Print -> Save as PDF); the JSON is the complete, schema-versioned
    event stream."""
    events = list(events)
    # Correction 2 (auto->CS depth fallback): if the bar/fast route hit its depth wall and
    # auto rerouted to Chouhy-Solotar, drop the abandoned partial bar worked steps so the
    # report shows a single coherent CS worked-steps set (the two Dispatch lines keep the
    # routing story). Done BEFORE the drift gate so it derives over the authoritative CS
    # steps by construction, not by last-wins luck. A no-op for every non-fallback stream.
    events = _drop_superseded(events)
    # HH drift gate + authoritative Result (Plan 34 fix): cross-check the event-derived
    # dims against the engine's returned dims and record the authoritative dims/variance
    # as a ResultDims event so every renderer shows the engine's numbers with the correct
    # HH^/HH_ label (and the fast GF(p) empty-steps report still gets a Result line).
    events = _authoritative_result(events, table)
    out = pathlib.Path(out_dir) if out_dir is not None else (pathlib.Path.cwd() / "quiverlab_traces")
    out.mkdir(parents=True, exist_ok=True)
    stem = "%s_%s" % (_SAFE_STEM.get(kind, kind), _hash(algebra, kind, top))
    title = "%s of %s" % (kind, repr(algebra).splitlines()[0])
    # The JSON machine record (Plan 34): the complete event stream, deterministic and
    # schema-versioned. A pure function of the events, byte-identical for identical input.
    (out / (stem + ".json")).write_text(
        render_json(events, title=title, references=references, algebra=algebra))
    html = out / (stem + ".html")
    html.write_text(render_html(events, title=title, references=references, algebra=algebra))
    print("Worked steps: %s (HTML, no JavaScript; print to PDF from your browser)"
          % _rel(html))
    return str(html)


def _rel(p):
    try:
        return str(p.relative_to(pathlib.Path.cwd()))
    except ValueError:
        return str(p)
