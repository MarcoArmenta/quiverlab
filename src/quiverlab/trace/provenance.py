"""Trace provenance: which literature a computation used, resolved through
Plan 06's bibliography() registry (spec §3.9; Plan-07 owns the RENDERING, Plan 06
owns the registry data). ENGINE_REFERENCES maps an engine `route` to Plan 06's
lowercase REGISTRY keys (never BibTeX ids); the Task-1 freshness gate asserts every
such key resolves. bibliography() returns a Bibliography whose `.keys` is a tuple
and whose iteration yields entry views (.key / .formatted / .bibtex_key / ...);
there is no subscripting and no .keys() method, so we build a key->entry map by
iterating."""
import quiverlab

ENGINE_REFERENCES = {
    "normalized bar complex": ("bar",),
    "hanlab fast GF(p) rank": ("bar",),
    "bardzell": ("bardzell",),
    "chouhy-solotar": ("chouhy_solotar",),
}

# Cup/bracket (Tamarkin-Tsygan) operation references are wired when the cup/bracket
# trace lands (Plan 04+) against Plan 06's operation REGISTRY keys; empty until then,
# so a verbose run never asks the registry for a key it does not have.
_OPERATION_REFERENCES = ()


def references_for(events):
    """Plan-06 REGISTRY keys implied by a trace: the engine (from every Dispatch
    whose route names a known engine) plus operation keys when LiftSteps are present.
    These registry keys are what a result's `.references` stores (they resolve to
    formatted citations at render time via resolve_references)."""
    from quiverlab.trace.events import Dispatch, LiftStep
    keys = []
    for e in events:
        if isinstance(e, Dispatch):
            for k in ENGINE_REFERENCES.get(e.route, ()):
                if k not in keys:
                    keys.append(k)
    if any(isinstance(e, LiftStep) for e in events):
        for k in _OPERATION_REFERENCES:
            if k not in keys:
                keys.append(k)
    return tuple(keys)


def _entries_by_key():
    """Map registry key -> entry view by iterating bibliography() (no subscripting)."""
    return {e.key: e for e in quiverlab.bibliography()}


def resolve_references(keys):
    """(bibtex_key, formatted) pairs for a tuple of REGISTRY keys, read from Plan 06's
    bibliography(). A registry key absent from the bibliography raises loudly (never
    silently dropped) -- the Task-1 freshness gate keeps this from firing."""
    by_key = _entries_by_key()
    pairs = []
    for k in keys:
        if k not in by_key:
            raise KeyError(
                "citation registry has no key %r (Plan 06 bibliography drift; update "
                "ENGINE_REFERENCES / the freshness gate)" % (k,))
        e = by_key[k]
        pairs.append((getattr(e, "bibtex_key", k), e.formatted))
    return tuple(pairs)
