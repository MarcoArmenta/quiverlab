"""quiverlab GUI runner (Plan 10): executes landing-page GUI requests.

Runs IDENTICALLY under CPython (pytest, tests/gui/) and Pyodide (the docs-site
Web Worker) — this exact file is shipped into the built site and imported in
the browser. Import policy: the public quiverlab surface + the sanctioned trace
helpers (render_html / render_json / references_for / resolve_references);
pinned by tests/gui/test_interface_freshness.py. quiverlab.engine.* is forbidden.

All public functions take and return JSON STRINGS (postMessage-friendly)."""
import json
import os
import random
import traceback
from fractions import Fraction

os.environ.setdefault("MPLBACKEND", "Agg")   # never let matplotlib probe for a display

import quiverlab

SCHEMA_VERSION = 1
# The GUI tags requests carrying module / Ext / Tor blocks as schema 2, and a
# Hochschild coefficient block (Plan 52) as schema 3; this runner's dispatch
# handles all three identically.
ACCEPTED_SCHEMAS = (1, 2, 3)
# The Hochschild compute kinds that may carry a coefficient bimodule (Plan 52).
_HH_COEFFICIENT_KINDS = frozenset({"hh_cohomology", "hh_homology"})
MAX_DEGREE = 10
# Depth to which projective dimension is probed before reporting "infinite"
# (matches the library's injective_dimension(bound=32) default).
_PD_BOUND = 32
# Ceiling on a random module's total dimension -- mirrors spec.py::_MAX_MODULE_DIM
# (the webapp module-block parse cap): each arrow's dense action is n x n cells.
_MAX_RANDOM_MODULE_DIM = 2048

# Module compute kinds (Plan 26; Plan 30 adds `tor`/`decompose`). `ext`/`tor` also
# need a second module (`ext_target`/`tor_target`, the latter a LEFT A-module).
_MODULE_KINDS = frozenset({
    "dimension_vector", "rad_top_soc", "ext", "tor", "tau", "tau_minus",
    "projective_resolution", "injective_resolution",
    "projective_dimension", "injective_dimension", "decompose", "almost_split",
    "tilting_check", "orbit_geometry", "barcode",
})

_state = {"algebra": None, "request": None, "events": None, "results": None,
          "module": None, "ext_target": None, "tor_target": None,
          "algebra_b": None, "coefficients": None, "new_arrows": None}


class RequestError(Exception):
    """Invalid GUI request (schema violation) — reported like a library error."""


def _fail(exc):
    if isinstance(exc, (quiverlab.QuiverlabError, RequestError)):
        return {"ok": False,
                "error": {"type": type(exc).__name__, "message": str(exc)}}
    # Unexpected: generic message to the page, full traceback for the console.
    return {"ok": False,
            "error": {"type": "InternalError",
                      "message": "unexpected engine error — details in the browser console"},
            "detail": traceback.format_exc()}


def _field_from_spec(spec):
    kind = spec.get("kind") if isinstance(spec, dict) else None
    if kind == "CC":
        return quiverlab.CC
    if kind == "QQ":
        # QQ (the rational field) is an exact input field alongside CC/GF; not a
        # top-level export, so import it from quiverlab.fields.
        from quiverlab.fields import QQ
        return QQ
    if kind == "GF":
        p, n = spec.get("p"), spec.get("n", 1)
        if not (isinstance(p, int) and isinstance(n, int) and p >= 2 and n >= 1):
            raise RequestError("field GF needs integers p >= 2 and n >= 1")
        return quiverlab.GF(p ** n)   # FieldError (p not prime, ...) surfaces verbatim
    raise RequestError("unknown field kind %r (expected 'CC', 'GF' or 'QQ')" % (kind,))


def _algebra_from_spec(alg):
    """Build a quiverlab Algebra from a GUI ``algebra`` block (kind 'quiver').
    Shared by :func:`run_build` and :func:`random_module`."""
    kind = alg.get("kind")
    if kind == "family":
        # SkewGentleAlgebra (Plan 68) is the ONE non-scalar constructor the no-code GUI
        # builds client-side: a triple (Q, I, Sp) drawn on the canvas + special-loop
        # picks. Every OTHER family stays the server tier. Byte-identical to the server
        # build (quiverlab.hpc.spec._build_skew_gentle).
        if alg.get("family") == "SkewGentleAlgebra":
            return _skew_gentle_from_spec(alg)
        # SkewGroupAlgebra (Plan 74): a base kQ/I + an explicit finite group action.
        # Byte-identical to the server build (both call the shared library builder
        # quiverlab.families.skew_group.build_skew_group_from_params).
        if alg.get("family") == "SkewGroupAlgebra":
            return _skew_group_from_spec(alg)
        # IncidenceAlgebra (Plan 75): the poset input mode's constructor -- cover pairs
        # (+ optional isolated elements) straight into the shipped library family, which
        # is what stashes the `_poset` provenance the order-complex route requires.
        # Byte-identical to the server build (both call quiverlab.IncidenceAlgebra).
        if alg.get("family") == "IncidenceAlgebra":
            return _incidence_from_spec(alg)
        raise RequestError("algebra kind 'family' is the server tier (Plan 09); "
                           "this GUI submits kind 'quiver' only")
    if kind != "quiver":
        raise RequestError("unknown algebra kind %r (expected 'quiver')" % (kind,))
    vertices = alg.get("vertices")
    if not (isinstance(vertices, list) and vertices
            and all(isinstance(v, int) for v in vertices)):
        raise RequestError("algebra.vertices must be a non-empty list of integers")
    arrows = alg.get("arrows")
    if not (isinstance(arrows, dict) and all(
            isinstance(st, list) and len(st) == 2
            and all(isinstance(x, int) for x in st)
            for st in arrows.values())):
        raise RequestError("algebra.arrows must map names to [source, target] pairs")
    relations = alg.get("relations", [])
    if not (isinstance(relations, list)
            and all(isinstance(r, str) for r in relations)):
        raise RequestError("algebra.relations must be a list of strings")
    potential = alg.get("potential")
    if potential is not None and not isinstance(potential, str):
        raise RequestError("algebra.potential must be a string (a k-linear sum of "
                           "oriented cycles, e.g. 'a*b*c - d*e*f')")
    field = _field_from_spec(alg.get("field"))
    Q = quiverlab.Quiver(vertices=vertices,
                         arrows={k: (s, t) for k, (s, t) in arrows.items()})
    if isinstance(potential, str) and potential.strip():
        # A potential routes through the P44 Jacobian kQ/(d_a W); its relations ARE
        # the cyclic derivatives, so explicit relations are then forbidden. Terms are
        # parsed with the shared relation-term parser (loud RelationError on an
        # unknown arrow / non-composable factor); Potential certifies each is a cycle.
        if relations:
            raise RequestError(
                "a 'potential' and explicit 'relations' cannot both be given: the "
                "Jacobian algebra's relations ARE the cyclic derivatives of the "
                "potential (leave 'relations' empty when a 'potential' is set)")
        from quiverlab.combinat.relations import _parse_term, _split_terms
        terms = [_parse_term(t, Q) for t in _split_terms(potential)]
        W = quiverlab.Potential(Q, terms)
        return quiverlab.JacobianAlgebra(Q, W, field=field)
    return Q.algebra(relations=relations, field=field)


def _skew_gentle_from_spec(alg):
    """Build the split algebra SkewGentleAlgebra(Q, I, Sp) from a GUI
    ``family: SkewGentleAlgebra`` block with flattened triple params (vertices, arrows,
    relations, special). Byte-identical to the server build
    (quiverlab.hpc.spec._build_skew_gentle)."""
    params = alg.get("params") or {}
    verts = params.get("vertices")
    if not (isinstance(verts, list) and verts and all(isinstance(v, int) for v in verts)):
        raise RequestError("SkewGentleAlgebra.vertices must be a non-empty list of ints")
    arrows_p = params.get("arrows") or {}
    if not (isinstance(arrows_p, dict) and all(
            isinstance(st, list) and len(st) == 2 and all(isinstance(x, int) for x in st)
            for st in arrows_p.values())):
        raise RequestError("SkewGentleAlgebra.arrows must map names to [source, target] pairs")
    rels = params.get("relations", [])
    if not (isinstance(rels, list) and all(isinstance(r, str) for r in rels)):
        raise RequestError("SkewGentleAlgebra.relations must be a list of strings")
    special = params.get("special", [])
    if not (isinstance(special, list) and all(isinstance(v, int) for v in special)):
        raise RequestError("SkewGentleAlgebra.special must be a list of vertex integers")
    field = _field_from_spec(alg.get("field"))
    from quiverlab.skewgentle.split import SkewGentleAlgebra
    Q = quiverlab.Quiver(vertices=list(verts),
                         arrows={k: (s, t) for k, (s, t) in arrows_p.items()})
    return SkewGentleAlgebra(quiver=Q, relations=list(rels), special=set(special),
                             field=field)


def _incidence_from_spec(alg):
    """Build the incidence algebra kP from a GUI ``family: IncidenceAlgebra`` block
    (params ``poset_or_covers`` + optional ``elements``) -- the Plan-75 poset input mode.
    Byte-identical to the server build: both call ``quiverlab.IncidenceAlgebra``, which
    validates the cover data (a directed cycle raises, surfaced here as a clean
    RequestError) and stashes the poset provenance."""
    from quiverlab.errors import QuiverlabError
    params = alg.get("params") or {}
    covers = params.get("poset_or_covers")
    if not isinstance(covers, list) or not all(
            isinstance(c, (list, tuple)) and len(c) == 2 for c in covers):
        raise RequestError("IncidenceAlgebra.poset_or_covers must be a list of cover "
                           "pairs [a, b], meaning a < b is a cover")
    elements = params.get("elements")
    if elements is not None and not isinstance(elements, list):
        raise RequestError("IncidenceAlgebra.elements must be a list or null")
    field = _field_from_spec(alg.get("field"))
    try:
        return quiverlab.IncidenceAlgebra([list(c) for c in covers], elements, field=field)
    except QuiverlabError as exc:
        raise RequestError(str(exc)) from exc


def _skew_group_from_spec(alg):
    """Build the skew group algebra A⋊G from a GUI ``family: SkewGroupAlgebra`` block
    (flattened params: vertices, arrows, relations, generators). Byte-identical to the
    server build -- both delegate to the shared library builder
    quiverlab.families.skew_group.build_skew_group_from_params (Plan 74)."""
    from quiverlab.errors import QuiverlabError
    from quiverlab.families.skew_group import build_skew_group_from_params
    params = dict(alg.get("params") or {})
    # canonical generator-order normalization (mirror the server schema): sort the
    # generators list by a stable JSON key so two orderings build the same algebra.
    gens = params.get("generators")
    if isinstance(gens, list):
        import json as _json
        params["generators"] = sorted(
            gens, key=lambda g: _json.dumps(g, sort_keys=True, default=str))
    field = _field_from_spec(alg.get("field"))
    try:
        return build_skew_group_from_params(params, field)
    except QuiverlabError as exc:
        raise RequestError(str(exc)) from exc


def run_build(request_json):
    """Parse + validate a schema-1 request, build the algebra, reset all state."""
    _state.update(algebra=None, request=None, events=[], results=[],
                  module=None, ext_target=None, tor_target=None, algebra_b=None,
                  coefficients=None, new_arrows=None)
    quiverlab.verbose = False   # the GUI renders its own report; never write trace files
    try:
        req = json.loads(request_json)
        if req.get("schema") not in ACCEPTED_SCHEMAS:
            raise RequestError("unsupported schema %r (this GUI speaks schema 1/2/3)"
                               % (req.get("schema"),))
        alg = req.get("algebra") or {}
        A = _algebra_from_spec(alg)
        # For the family: SkewGentleAlgebra path the drawn quiver rides under params;
        # the summary reports the ORIGINAL (Q, arrows), not the split.
        _src = alg.get("params") if alg.get("kind") == "family" else alg
        vertices = _src.get("vertices") or []
        arrows = _src.get("arrows") or {}
        # Plan 75: the IncidenceAlgebra family's params are POSET data (covers), not a
        # quiver, so the summary counts come off the BUILT Hasse quiver instead of the
        # request -- reporting 0 vertices / 0 arrows for a real algebra would be a lie.
        if not vertices and getattr(A, "quiver", None) is not None:
            vertices = list(A.quiver.vertices)
            arrows = dict(A.quiver.arrows)
        # Module blocks (Plan 26) ride alongside the algebra; the module itself is
        # built lazily in compute_one, so a relation-violating matrix surfaces as a
        # per-computation error (rendered on the page), never a build crash.
        # wave 2: the SECOND algebra for derived_compare rides alongside the request
        # (spec dict, built lazily in compute_one -- a bad B surfaces as a per-
        # computation error, like the module blocks, never a build crash).
        _state.update(algebra=A, request=req, module=req.get("module"),
                      ext_target=req.get("ext_target"),
                      tor_target=req.get("tor_target"),
                      algebra_b=req.get("algebra_b"),
                      coefficients=req.get("coefficients"),
                      new_arrows=req.get("new_arrows"))
        out = {"ok": True, "dim": A.dim, "n_vertices": len(vertices),
               "n_arrows": len(arrows), "algebra": repr(A).splitlines()[0]}
    except Exception as exc:
        out = _fail(exc)
    return json.dumps(out)


def _parse_compute(spec):
    name, _, rng = spec.partition(":")
    top = None
    # tau_tilting carries a PAIR BUDGET, not a degree range (Plan 45): 'tau_tilting' or
    # 'tau_tilting:512'. The budget is not a homological degree, so it skips MAX_DEGREE.
    if name == "tau_tilting":
        if rng and not rng.isdigit():
            raise RequestError("tau_tilting budget must be a positive integer (got %r)"
                               % (spec,))
        return "tau_tilting", (int(rng) if rng else None)
    # congruences carries a PAIR BUDGET, not a degree range (Plan 64): 'congruences' or
    # 'congruences:512'. The budget is not a homological degree, so it skips MAX_DEGREE.
    if name == "congruences":
        if rng and not rng.isdigit():
            raise RequestError("congruences budget must be a positive integer (got %r)"
                               % (spec,))
        return "congruences", (int(rng) if rng else None)
    # tau_cluster carries a PAIR BUDGET, not a degree range (Plan 66): 'tau_cluster' or
    # 'tau_cluster:512'. The budget is not a homological degree, so it skips MAX_DEGREE.
    # Server twin: quiverlab.hpc.spec parses the same form.
    if name == "tau_cluster":
        if rng and not rng.isdigit():
            raise RequestError("tau_cluster budget must be a positive integer (got %r)"
                               % (spec,))
        return "tau_cluster", (int(rng) if rng else None)
    # hh1_lie carries a DIM BUDGET, not a degree range (Plan 70): 'hh1_lie' or
    # 'hh1_lie:48'. The budget caps A.dim for the Der solve, not a homological degree,
    # so it skips MAX_DEGREE. Server twin: quiverlab.hpc.spec parses the same form.
    if name == "hh1_lie":
        if rng and not rng.isdigit():
            raise RequestError("hh1_lie budget must be a positive integer (got %r)"
                               % (spec,))
        return "hh1_lie", (int(rng) if rng else None)
    # deformations carries a DIM BUDGET, not a degree range (Plan 78): 'deformations' or
    # 'deformations:32'. The budget caps A.dim for the CS obstruction bracket (a coarse DoS
    # backstop -- the real cost is HH^2/HH^3 richness x resolution size), not a homological
    # degree, so it skips MAX_DEGREE. Server twin: quiverlab.hpc.spec parses the same form.
    if name == "deformations":
        if rng and not rng.isdigit():
            raise RequestError("deformations budget must be a positive integer (got %r)"
                               % (spec,))
        return "deformations", (int(rng) if rng else None)
    # wall_chamber carries a PAIR BUDGET too (Plan 63): 'wall_chamber' or 'wall_chamber:512'
    # -- the exchange-graph pair budget, not a homological degree, so it skips MAX_DEGREE.
    if name == "wall_chamber":
        if rng and not rng.isdigit():
            raise RequestError("wall_chamber budget must be a positive integer (got %r)"
                               % (spec,))
        return "wall_chamber", (int(rng) if rng else None)
    # cluster_category (Plan 79 / R31) carries the EXCHANGE-GRAPH PAIR BUDGET, not a
    # degree: 'cluster_category' or 'cluster_category:512'. Server twin: quiverlab.hpc.spec
    # and webapp.server.schema parse the same special form.
    if name == "cluster_category":
        if rng and not rng.isdigit():
            raise RequestError("cluster_category budget must be a positive integer "
                               "(got %r)" % (spec,))
        return "cluster_category", (int(rng) if rng else None)
    # silting carries a RADIUS,BUDGET pair, not a degree range (Plan 67): 'silting' or
    # 'silting:3,64'. The top is the (radius, budget) tuple; neither is a homological
    # degree, so it skips MAX_DEGREE. Server twin: quiverlab.hpc.spec parses the same form.
    if name == "silting":
        if rng:
            parts = rng.split(",")
            if len(parts) != 2 or not all(p.isdigit() for p in parts):
                raise RequestError("silting suffix must be 'radius,budget' with positive "
                                   "integers (got %r)" % (spec,))
            return "silting", (int(parts[0]), int(parts[1]))
        return "silting", None
    # ar_quiver carries a MODULE BUDGET, not a degree range (wave 2): 'ar_quiver' or
    # 'ar_quiver:512'. The budget is not a homological degree, so it skips MAX_DEGREE.
    if name == "ar_quiver":
        if rng and not rng.isdigit():
            raise RequestError("ar_quiver budget must be a positive integer (got %r)"
                               % (spec,))
        return "ar_quiver", (int(rng) if rng else None)
    # split_extension / arrow_removal (Plan 72) carry a TOP-DEGREE budget, not a lo..hi
    # range: 'split_extension' / 'split_extension:6'. Skips MAX_DEGREE like ar_quiver.
    if name in ("split_extension", "arrow_removal"):
        if rng and not rng.isdigit():
            raise RequestError("%s budget must be a positive integer (got %r)"
                               % (name, spec))
        return name, (int(rng) if rng else None)
    # skew_group_hh carries a TOP-DEGREE budget, not a lo..hi range (Plan 74):
    # 'skew_group_hh' / 'skew_group_hh:3'. Skips MAX_DEGREE like split_extension.
    if name == "skew_group_hh":
        if rng and not rng.isdigit():
            raise RequestError("skew_group_hh budget must be a positive integer (got %r)"
                               % (spec,))
        return "skew_group_hh", (int(rng) if rng else None)
    # exceptional_sequences carries an ENUMERATION BUDGET, not a degree range (Plan 65):
    # 'exceptional_sequences' or 'exceptional_sequences:512'. Skips MAX_DEGREE like tau_tilting.
    if name == "exceptional_sequences":
        if rng and not rng.isdigit():
            raise RequestError(
                "exceptional_sequences budget must be a positive integer (got %r)" % (spec,))
        return "exceptional_sequences", (int(rng) if rng else None)
    # radical_filtration + ar_invariants (Plan 57) carry a MODULE BUDGET, not a degree
    # range (parsed like ar_quiver -- skips MAX_DEGREE). NOTE: 'radical_filtration'
    # (the module-category radical rad^n(X,Y)) is DISTINCT from 'radical_filtration_ss'
    # (the Loewy radical-series spectral sequence, a DIFFERENT object).
    if name in ("radical_filtration", "ar_invariants"):
        if rng and not rng.isdigit():
            raise RequestError("%s budget must be a positive integer (got %r)"
                               % (name, spec))
        return name, (int(rng) if rng else None)
    # left_right_parts carries a MODULE BUDGET, not a degree range (Plan 55): the budget caps
    # the knitted universe, so it skips MAX_DEGREE (like ar_quiver / tau_tilting).
    if name == "left_right_parts":
        if rng and not rng.isdigit():
            raise RequestError("left_right_parts budget must be a positive integer (got %r)"
                               % (spec,))
        return "left_right_parts", (int(rng) if rng else None)
    # tilted_check carries the KNIT budget (budget_modules), not a degree range (Plan 60): it
    # skips MAX_DEGREE like ar_quiver / left_right_parts. budget_sections stays internal.
    if name == "tilted_check":
        if rng and not rng.isdigit():
            raise RequestError("tilted_check budget must be a positive integer (got %r)"
                               % (spec,))
        return "tilted_check", (int(rng) if rng else None)
    # recognizer_ladder carries a MODULE BUDGET, not a degree range (Plan 61): the budget caps
    # the knitted universe, so it skips MAX_DEGREE (like left_right_parts / ar_quiver).
    if name == "recognizer_ladder":
        if rng and not rng.isdigit():
            raise RequestError("recognizer_ladder budget must be a positive integer (got %r)"
                               % (spec,))
        return "recognizer_ladder", (int(rng) if rng else None)
    # skew_gentle carries a tau-tilting PAIR BUDGET, not a degree range (Plan 68):
    # 'skew_gentle' or 'skew_gentle:512'. Skips MAX_DEGREE (like tau_tilting).
    if name == "skew_gentle":
        if rng and not rng.isdigit():
            raise RequestError("skew_gentle budget must be a positive integer (got %r)"
                               % (spec,))
        return "skew_gentle", (int(rng) if rng else None)
    if rng:
        lo, _, hi = rng.partition("..")
        if lo != "0" or not hi.isdigit():
            raise RequestError("bad compute range %r (expected 'name:0..N')" % (spec,))
        top = int(hi)
        if top > MAX_DEGREE:
            raise RequestError("degree cap is %d (got %d)" % (MAX_DEGREE, top))
    return name, top


_TATE_REFERENCES = ("bergh_jorgensen_tate", "wang_singular_hh", "keller_singular_hh",
                    "usui_tate_periodic")


def _citation_pairs(keys):
    from quiverlab.trace.provenance import resolve_references
    return [list(p) for p in resolve_references(tuple(keys))]


def _latex_matrix(rows):
    body = r" \\ ".join(" & ".join(str(x) for x in row) for row in rows)
    return r"\begin{pmatrix} %s \end{pmatrix}" % body


# --- module block (Plan 26): build a module from the request + dispatch ------
# The webapp server tier (webapp/server/runner.py) carries the SAME logic; this
# is the Pyodide/client copy (the two runners already duplicate the algebra
# dispatch -- they cannot import each other).

def _parse_entry(x):
    """Parse an exact matrix entry to int/Fraction. Entries are DATA (never
    evaluated); a float is refused (exactness is the point)."""
    if isinstance(x, bool):
        raise RequestError("module entry %r is not a number" % (x,))
    if isinstance(x, int):
        return x
    if isinstance(x, float):
        raise RequestError("module entry %r is a float; entries must be exact "
                           "(integers or strings like '1/2')" % (x,))
    if isinstance(x, str):
        s = x.strip()
        try:
            return int(s)
        except ValueError:
            pass
        try:
            return Fraction(s)
        except (ValueError, ZeroDivisionError):
            raise RequestError("module entry %r is not an exact integer or "
                               "fraction" % (x,))
    raise RequestError("module entry %r is not a number" % (x,))


def _full_matrices(A, mspec):
    """Expand per-arrow BLOCK matrices (dim[target] x dim[source] in the
    representation quiver) into the full vertex-ordered arrow actions A.module
    consumes. side selects the representation quiver (right=A, left=A^op)."""
    rep = A if mspec.get("side", "right") == "right" else A.opposite()
    verts = list(rep.quiver.vertices)
    by_str = {str(v): v for v in verts}
    dimvec = {v: 0 for v in verts}
    for key, n in (mspec.get("dims") or {}).items():
        if key not in by_str:
            raise RequestError("module: no vertex %r in the algebra" % (key,))
        dimvec[by_str[key]] = int(n)
    start, off = {}, 0
    for v in verts:
        start[v] = off
        off += dimvec[v]
    n = off
    arrow_names = list(rep.quiver.arrows)
    for a in (mspec.get("maps") or {}):
        if a not in arrow_names:
            raise RequestError("module: no arrow %r in the algebra" % (a,))
    action = {}
    for a in arrow_names:
        s, t = rep.quiver.source(a), rep.quiver.target(a)
        rows, cols = dimvec[t], dimvec[s]
        full = [[0] * n for _ in range(n)]
        block = (mspec.get("maps") or {}).get(a)
        if block is not None:
            if len(block) != rows or any(len(r) != cols for r in block):
                raise RequestError("module map %r must be %dx%d (target x source "
                                   "dims)" % (a, rows, cols))
            for i in range(rows):
                for j in range(cols):
                    full[start[t] + i][start[s] + j] = _parse_entry(block[i][j])
        action[a] = full
    return dimvec, action


def _build_module(A, mspec, name):
    if not isinstance(mspec, dict):
        raise RequestError("this computation needs a module block")
    b = mspec.get("builtin")
    if b is not None:
        v, kind, side = b.get("vertex"), b.get("kind"), mspec.get("side", "right")
        builder = {"simple": A.simple, "projective": A.projective,
                   "injective": A.injective}.get(kind)
        if builder is None:
            raise RequestError("unknown builtin module kind %r" % (kind,))
        for vv in A.quiver.vertices:
            if vv == v or str(vv) == str(v):
                return builder(vv, side=side)
        raise RequestError("module: no vertex %r in the algebra" % (v,))
    dimvec, action = _full_matrices(A, mspec)
    return A.module(dimvec, action, side=mspec.get("side", "right"), name=name)


def _build_coefficient(A, cspec):
    """Build a Hochschild coefficient Bimodule from a request coefficient block
    (Plan 52). None => the regular bimodule (no block emitted). Mirrors
    quiverlab.hpc.spec._build_coefficient / _parse_coefficients."""
    if cspec is None:
        return None
    from quiverlab.hochschild.coefficients import Bimodule
    if not isinstance(cspec, dict):
        raise RequestError("coefficients must be an object")
    builtin = cspec.get("builtin")
    if builtin is not None:
        kind = builtin.get("kind") if isinstance(builtin, dict) else None
        builder = {"regular": Bimodule.regular, "dual": Bimodule.dual,
                   "twisted_nakayama": Bimodule.twisted_by_nakayama,
                   "quotient_socle": Bimodule.mod_socle}.get(kind)
        if builder is None:
            raise RequestError("unknown coefficient builtin kind %r" % (kind,))
        return builder(A)
    dim, lm, rm = cspec.get("dim"), cspec.get("left_maps"), cspec.get("right_maps")
    if dim is None or lm is None or rm is None:
        raise RequestError("coefficients: needs a 'builtin' pick-list OR "
                           "'dim'+'left_maps'+'right_maps'")
    return Bimodule.from_actions(A, dim, lm, rm)


def _random_module_core(A, dims, side, seed, tries):
    """Draw a random representation of ``A`` (GitHub #3). Byte-identical twin of
    ``quiverlab.hpc.spec.random_module``: SAME draw order (sorted arrow names,
    row-major) under ``random.Random(seed)``, so both tiers agree for a fixed seed.
    Exact entries only -- GF(p): uniform ``0..p-1`` (the prime subfield, the only
    field elements the entry grammar can name); char 0: integers ``-5..5``. With
    relations the draw is rejection-sampled (the first that builds through the
    library's checked ``A.module`` wins); over char 0 with relations a random point
    never satisfies them exactly, so it refuses with ``{"error": "char0"}``."""
    if side not in ("right", "left"):
        raise RequestError("random module: side must be 'right' or 'left'")
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise RequestError("random module: seed must be an integer")
    if not isinstance(tries, int) or isinstance(tries, bool) or tries < 1:
        raise RequestError("random module: tries must be a positive integer")
    rep = A if side == "right" else A.opposite()
    by_str = {str(v): v for v in rep.quiver.vertices}
    dimvec = {v: 0 for v in rep.quiver.vertices}
    total = 0
    for key, nv in (dims or {}).items():
        if str(key) not in by_str:
            raise RequestError("random module: no vertex %r in the algebra" % (key,))
        if not isinstance(nv, int) or isinstance(nv, bool) or nv < 0:
            raise RequestError("random module: dim[%r] must be a non-negative integer"
                               % (key,))
        dimvec[by_str[str(key)]] = nv
        total += nv
    if total > _MAX_RANDOM_MODULE_DIM:
        raise RequestError("random module: total dimension %d exceeds the %d cap"
                           % (total, _MAX_RANDOM_MODULE_DIM))
    char = A.domain.characteristic
    has_rel = bool(A.relations)
    if has_rel and char == 0:
        return {"error": "char0"}
    dims_str = {str(k): int(v) for k, v in (dims or {}).items()}
    arrow_names = sorted(rep.quiver.arrows)

    def _draw(rng):
        maps = {}
        for a in arrow_names:
            rows = dimvec[rep.quiver.target(a)]
            cols = dimvec[rep.quiver.source(a)]
            if char:
                maps[a] = [[rng.randrange(char) for _ in range(cols)]
                           for _ in range(rows)]
            else:
                maps[a] = [[rng.randint(-5, 5) for _ in range(cols)]
                           for _ in range(rows)]
        return maps

    rng = random.Random(seed)
    for k in range(1, (1 if not has_rel else tries) + 1):
        maps = _draw(rng)
        try:
            _build_module(A, {"dims": dims_str, "maps": maps, "side": side}, "random")
        except quiverlab.QuiverlabError:
            continue                       # relation-violating draw: try again
        return {"maps": maps, "seed": seed, "tries": k}
    return {"error": "budget", "tries": tries}


def random_module(request_json):
    """JSON-string entry for the Pyodide worker's ``random`` verb:
    ``{algebra, dims, side?, seed?, tries?}`` -> ``{maps, seed, tries}`` |
    ``{error}``. The server tier's twin is POST /api/gui/random-module."""
    try:
        req = json.loads(request_json)
        A = _algebra_from_spec(req.get("algebra") or {})
        out = _random_module_core(A, req.get("dims") or {},
                                  req.get("side", "right"),
                                  req.get("seed", 0), req.get("tries", 200))
    except Exception as exc:
        out = _fail(exc)
    return json.dumps(out)


def _tor_target_spec():
    """The Tor second module N (a LEFT A-module). An omitted side defaults to
    ``"left"`` -- mirrors the schema/spec-core default -- so the GUI need not force
    the toggle. A builtin's own nested side is honored as-is."""
    ts = _state.get("tor_target")
    if isinstance(ts, dict):
        b = ts.get("builtin")
        has_side = "side" in ts or (isinstance(b, dict) and "side" in b)
        if not has_side:
            ts = dict(ts, side="left")
    return ts


def _dv(dimvec):
    return {str(v): int(n) for v, n in sorted(dimvec.items(), key=lambda kv: str(kv[0]))}


def _dv_latex(dimvec):
    d = _dv(dimvec)
    return "(" + ",\\, ".join(str(d[k]) for k in d) + ")" if d else "()"


def _homdim_latex(op, value):
    """Display latex for a homological-dimension block (``pd``/``id``), mirroring
    ``quiverlab.hpc.spec._homdim_latex``. An UNRESOLVED probe is not a proof of
    infinity -- the resolution merely did not terminate by ``_PD_BOUND`` -- so it
    states the certified lower bound, exactly as ``global_dimension`` does."""
    if value is None:
        return r"\operatorname{%s} M > %d" % (op, _PD_BOUND)
    return r"\operatorname{%s} M = %d" % (op, value)


_HOMDIM_UNRESOLVED = ("certified lower bound; the resolution did not terminate "
                      "within the probed depth %d" % _PD_BOUND)


def _mod_view(m):
    return {"dimvec": _dv(m.dimension_vector()), "dim": m.dim}


def _mod_repr(m):
    """A module as the no-code INPUT schema ``{"dims": {v: n}, "maps": {arrow: [[..]]}}``
    (Plan 34, Marco): the per-vertex dimension VECTOR + the exact per-arrow action
    matrices, redundant total dim dropped. Routed through the SAME library serializer as
    the server tier (``quiverlab.modules.qpa_module.module_blocks``) so the two runners
    cannot drift."""
    from quiverlab.modules.qpa_module import module_blocks
    return module_blocks(m)


_MOD_REFS = {
    "dimension_vector": ["assem_book"], "rad_top_soc": ["assem_book"],
    "tau": ["assem_book"], "tau_minus": ["assem_book"], "ext": ["module_ext"],
    "tor": ["minimal_resolution", "module_ext"],
    "projective_resolution": ["minimal_resolution"],
    "projective_dimension": ["minimal_resolution"],
    "injective_resolution": ["minimal_resolution", "assem_book"],
    "injective_dimension": ["minimal_resolution", "assem_book"],
    "decompose": ["assem_book"],
    "almost_split": ["assem_book", "ars_book"],
    "tilting_check": ["bongartz_tilting", "assem_book"],
    "orbit_geometry": ["voigt_rigidity", "kac_canonical",
                       "schofield_general_reps", "derksen_weyman_canonical"],
    "barcode": ["escolar_hiraoka", "botnan_crawley_boevey", "gabriel", "assem_book"],
}


# The additive-tau note carried by an AR-translate block whose input decomposes.
# The block ships a stable KEY; each renderer localizes it (app.js via i18n, gui.js
# via its own English map) -- so the block stays language-neutral data.
_TAU_ADDITIVE_KEY = "mod.tau_additive"


def _summands_latex(vertices, letter):
    """A resolution term's summand multiset as LaTeX, e.g. ``P_{1}^{2} \\oplus P_{3}``
    (``letter`` = ``P`` projectives / ``I`` injectives). ``0`` for the zero term."""
    if not vertices:
        return "0"
    counts = {}
    for v in vertices:
        counts[v] = counts.get(v, 0) + 1

    def key(v):
        return (0, v) if isinstance(v, int) and not isinstance(v, bool) else (1, str(v))
    parts = []
    for v in sorted(counts, key=key):
        base = "%s_{%s}" % (letter, v)
        parts.append(base if counts[v] == 1 else "%s^{%d}" % (base, counts[v]))
    return " \\oplus ".join(parts)


# mirrors quiverlab.hpc.spec._differential_blocks (same cap, same shape)
_MAX_DIFF_CELLS = 250_000


def _differential_blocks(res, n_terms):
    """The resolution's maps as exact matrices (rows: target basis, columns:
    source basis). Projective: entry 0 = eps: Q_0 -> M, entry n = d_n: Q_n ->
    Q_{n-1}. Injective: entry 0 = iota: M -> E^0, entry n = d^n: E^{n-1} -> E^n."""
    out = []
    dmats = getattr(res, "dmats", None) or []
    for n in range(min(n_terms, len(dmats))):
        D = res.differential(n)
        nrows = len(D)
        ncols = len(D[0]) if nrows else 0
        if nrows * ncols > _MAX_DIFF_CELLS:
            out.append({"rows": nrows, "cols": ncols, "elided": True})
        else:
            out.append({"rows": nrows, "cols": ncols,
                        "matrix": [[str(x) for x in row] for row in D]})
    return out


def _term_basis_blocks(res, kind, M):
    """The ordered k-basis of each resolution term, as the concatenated path bases of
    its projective / injective summands (Plan 35 UNIT 2). Shape-identical to the server
    tier (``quiverlab.hpc.spec._term_basis_blocks``): ``term_basis[n]`` lists one label
    per basis vector of term n, so ``len(term_basis[n])`` = the term dimension (the
    differential's column count for a projective resolution, row count for an injective
    one). ``None`` (field omitted) when the algebra carries no path basis or the total
    is implausibly large -- renderers tolerate the absence."""
    try:
        from quiverlab.modules.builders import projective
        if kind == "projective_resolution":
            base, sym = M.algebra, "P"
        else:
            from quiverlab.modules.opposite import opposite_algebra
            base, sym = opposite_algebra(M.algebra), "I"
        cache = {}

        def paths(v):
            if v not in cache:
                cache[v] = list(projective(base, v)._pv_basis_labels)
            return cache[v]

        out, total = [], 0
        for n in range(len(res.terms)):
            labels = []
            for v in res.term(n):
                labels.extend("%s_%s: %s" % (sym, v, p) for p in paths(v))
            total += len(labels)
            if total > _MAX_DIFF_CELLS:               # implausible; omit (bloat guard)
                return None
            out.append(labels)
        return out
    except Exception:
        return None


def _summand_view(mod, mult):
    """One indecomposable summand: dimension vector, multiplicity, and either a
    STANDARD name (S_v / P_v / I_v) or its full per-arrow matrices (Marco
    2026-07-29). Same library serializer as the server tier -- no drift."""
    from quiverlab.modules.qpa_module import summand_blocks
    return summand_blocks(mod, mult)


def _decompose_engine():
    """(is_indecomposable, decompose) from the Plan-30 Part-A engine, or ``None``
    when not yet importable (so callers degrade to the pre-Plan-30 block shape)."""
    try:
        from quiverlab.modules.decompose import decompose, is_indecomposable
    except ImportError:
        return None
    return is_indecomposable, decompose


def _input_certificate(M):
    """Certify the INPUT of an AR translate (Marco #1): indecomposable, or its
    Krull-Schmidt decomposition + the additivity note. Best-effort + honest: ``{}``
    (no claim) when the decompose engine is unavailable OR cannot certify within
    budget (it is LOUD, e.g. char <= dim M) -- tau itself stays valid, so the block
    omits the certificate rather than assert what it could not prove."""
    eng = _decompose_engine()
    if eng is None:
        return {}
    is_indec, decompose = eng
    try:
        if is_indec(M):
            return {"indecomposable": True}
        return {"indecomposable": False,
                "decomposition": [_summand_view(s, m) for (s, m) in decompose(M)],
                "note_key": _TAU_ADDITIVE_KEY}
    except quiverlab.QuiverlabError:
        return {}


def _ar_translate(mod, kind, name):
    """One AR translate as a self-contained display payload: the symbol, its
    dimension-vector latex, the FULL representation ({dims, maps}) and the input's
    indecomposability certificate. Mirrors ``quiverlab.hpc.spec._ar_translate`` so
    the two runners ship the same shape."""
    t = mod.tau() if kind == "tau" else mod.tau_minus()
    sym = (r"\tau %s" if kind == "tau" else r"\tau^{-} %s") % name
    latex = ((sym + " = 0") if t.dim == 0
             else (r"\underline{\dim}\, " + sym + " = "
                   + _dv_latex(t.dimension_vector())))
    out = {"name": name, "side": t.side, "is_zero": t.dim == 0,
           "latex": latex, **_mod_view(t)}
    if t.dim > 0:
        out["repr"] = _mod_repr(t)
    out.update(_input_certificate(mod))
    return out


def _target_translates(A, kind):
    """The AR translates of the SECOND module(s) N the request names (the Ext/Tor
    argument), so a tau block covers every module in play, with full matrices
    (Marco, 2026-07-29). A loud refusal on one target becomes an honest error entry
    -- tau M is already computed and stays valid."""
    out = []
    for role, spec in (("ext_target", _state.get("ext_target")),
                       ("tor_target", _tor_target_spec())):
        if not isinstance(spec, dict):
            continue
        try:
            entry = _ar_translate(_build_module(A, spec, "N"), kind, "N")
        except (quiverlab.QuiverlabError, RequestError) as exc:
            entry = {"name": "N", "error": str(exc)}
        entry["role"] = role
        out.append(entry)
    return out


def _module_block(name, top):
    """Dispatch one module compute kind against the built module(s). Blocks carry
    a `latex` display and `citations`, like the algebra invariants."""
    A = _state["algebra"]
    M = _build_module(A, _state.get("module"), "M")
    keys = _MOD_REFS[name]
    cites = _citation_pairs(keys)
    if name == "dimension_vector":
        return {"kind": name, "side": M.side, "references": list(keys),
                "citations": cites, **_mod_view(M),
                "latex": r"\underline{\dim}\, M = " + _dv_latex(M.dimension_vector())}
    if name == "rad_top_soc":
        # "series" = the Loewy (radical) series top-to-bottom (Plan 37), byte-identical
        # to the hpc spec core so the draw page and the report/CLI agree.
        return {"kind": name, "side": M.side, "references": list(keys),
                "citations": cites,
                "radical": _mod_repr(M.radical()), "top": _mod_repr(M.top()),
                "socle": _mod_repr(M.socle()),
                "series": [dict(layer) for layer in M.loewy_layers()]}
    if name in ("tau", "tau_minus"):
        # The translate ships as a full representation ({dims, maps}) -- mirrors the
        # hpc spec core, so both dispatches carry the AR translate's per-arrow
        # matrices -- together with the input certificate (Marco #1).
        entry = _ar_translate(M, name, "M")
        entry.pop("name", None)
        block = {"kind": name, "references": list(keys), "citations": cites, **entry}
        # ... and the same for the SECOND module N when the request names one
        # (Marco, 2026-07-29): a tau block covers every module in the request.
        targets = _target_translates(A, name)
        if targets:
            block["targets"] = targets
        return block
    if name == "decompose":
        eng = _decompose_engine()
        if eng is None:
            raise RequestError("module decomposition needs the Krull-Schmidt engine "
                               "(Plan 30 Part A); not available in this build")
        _, decompose = eng
        summands = [_summand_view(s, m) for (s, m) in decompose(M)]
        return {"kind": name, "side": M.side, "summands": summands,
                "iso_classes": len(summands), "references": list(keys),
                "citations": cites}
    if name == "barcode":
        # The persistence/TDA barcode (Plan 69 / R33). SAME shared core builder as the
        # hpc spec dispatch (quiverlab.modules.barcode.barcode_block) + references ->
        # citations, so the two runners emit byte-identical blocks (a refusal is an
        # {"error": ...} block, never a raise).
        from quiverlab.modules.barcode import barcode_block
        block = barcode_block(A, M)
        block["references"] = list(keys)
        block["citations"] = cites
        return block
    if name == "almost_split":
        # The almost-split sequence 0 -> tau M -> E -> M -> 0 for M indecomposable
        # non-projective (Plan 41). Byte-identical block shape to quiverlab.hpc.spec's
        # almost_split branch (tau = full repr; E's summands via the shared serializer;
        # honest exists:false refusal for projective/decomposable/undecidable input).
        from quiverlab.modules.decompose import decompose
        try:
            ses = M.almost_split_sequence()
        except quiverlab.QuiverlabError as exc:
            return {"kind": name, "exists": False, "reason": str(exc),
                    "references": list(keys), "citations": cites}
        return {"kind": name, "exists": True, "side": M.side,
                "tau": _mod_repr(ses.L),
                "middle": {"summands": [_summand_view(s, m)
                                        for (s, m) in decompose(ses.M)]},
                "M": _mod_view(M), "indecomposable": True,
                "latex": r"0 \to \tau M \to E \to M \to 0",
                "references": list(keys), "citations": cites}
    if name == "ext":
        if top is None:
            raise RequestError("ext needs a range, e.g. 'ext:0..4'")
        from quiverlab.modules.ext import ext_dims
        N = _build_module(A, _state.get("ext_target"), "N")
        # Plan 35 wave 3a: explicit Ext cocycle representatives + self-cert data,
        # additive keys shared byte-for-byte with the hpc spec runner.
        # Plan 35 wave 3c: interpret=True captures the Yoneda exact sequence of each class.
        raw, reps = ext_dims(A, M, N, top, with_reps=True, interpret=True)
        block = {"kind": name, "top": top, "dims": [int(d) for d in raw],
                 "target": _mod_view(N), "references": list(keys), "citations": cites,
                 # Marco 2026-08-03: WHICH module was resolved, by WHICH resolution
                 # (byte-identical to the hpc spec runner's stamp).
                 "resolved": {"module": "M", "side": M.side,
                              "resolution": "minimal projective resolution"}}
        block.update(reps)
        return block
    if name == "tor":
        if top is None:
            raise RequestError("tor needs a range, e.g. 'tor:0..4'")
        try:
            from quiverlab.modules.tor import tor_dims
        except ImportError:
            raise RequestError("Tor requires the Tor engine (Plan 29); not available "
                               "in this build")
        N = _build_module(A, _tor_target_spec(), "N")
        raw, reps = tor_dims(A, M, N, top, with_reps=True)
        block = {"kind": name, "top": top, "dims": [int(d) for d in raw],
                 "target": _mod_view(N), "references": list(keys), "citations": cites,
                 "resolved": {"module": "M", "side": M.side,
                              "resolution": "minimal projective resolution"}}
        block.update(reps)
        return block
    if name in ("projective_resolution", "injective_resolution"):
        if top is None:
            raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
        res = (M.projective_resolution(top) if name == "projective_resolution"
               else M.injective_resolution(top))
        letter = "P" if name == "projective_resolution" else "I"
        terms = [_dv(dv) for dv in res.dimension_vectors()]
        block = {"kind": name, "top": top, "terms": terms,
                 "betti": [res.betti(i) for i in range(len(terms))],
                 "summands": [_summands_latex(res.term(i), letter)
                              for i in range(len(terms))],
                 "differentials": _differential_blocks(res, len(terms)),
                 "references": list(keys), "citations": cites}
        term_basis = _term_basis_blocks(res, name, M)
        if term_basis is not None:
            block["term_basis"] = term_basis
        if name == "projective_resolution":
            block["pd"] = res.pd()
        else:
            block["injective_dimension"] = res.injective_dimension()
        return block
    if name == "projective_dimension":
        pd = M.projective_resolution(_PD_BOUND).pd()
        return {"kind": name, "value": pd, "finite": pd is not None,
                "references": list(keys), "citations": cites,
                "bound": _PD_BOUND, "latex": _homdim_latex("pd", pd),
                **({} if pd is not None else {"note": _HOMDIM_UNRESOLVED})}
    if name == "injective_dimension":
        idim = M.injective_dimension(bound=_PD_BOUND)
        return {"kind": name, "value": idim, "finite": idim is not None,
                "references": list(keys), "citations": cites,
                "bound": _PD_BOUND, "latex": _homdim_latex("id", idim),
                **({} if idim is not None else {"note": _HOMDIM_UNRESOLVED})}
    if name == "tilting_check":
        # The candidate T is the request's module M (schema-2, no second module). top = the
        # optional degree n, default 1. A decompose char-caveat becomes an honest per-block
        # error (mirrors the hpc spec runner byte-for-byte on the math subkeys).
        from quiverlab.errors import QuiverlabError
        from quiverlab.modules.tilting import is_tilting_module
        try:
            rep = is_tilting_module(M, n=(top if top is not None else 1))
        except QuiverlabError as exc:
            return {"kind": name, "error": str(exc),
                    "references": list(keys), "citations": cites}
        return {"kind": name, "is_tilting": rep.is_tilting, "n": rep.n, "pd": rep.pd,
                "self_ext_vanishes": rep.self_ext_vanishes,
                "num_summands": rep.num_summands, "num_vertices": rep.num_vertices,
                "note": rep.note, "references": list(keys), "citations": cites}
    if name == "orbit_geometry":
        # Byte-identical to quiverlab.hpc.spec's orbit_geometry branch: the shared
        # library builder + _with_refs's references(list)+citations(pairs). Orbit
        # dim / rigidity / codim for ANY module; the Kac canonical decomposition is
        # the hereditary-Dynkin extra. A char-scope QuiverlabError becomes an honest
        # per-block error (never fatal).
        from quiverlab.invariants.geometry import orbit_geometry_block
        try:
            block = orbit_geometry_block(M)
        except quiverlab.QuiverlabError as exc:
            block = {"kind": name, "error": str(exc)}
        block["references"] = list(keys)
        block["citations"] = cites
        return block
    raise RequestError("unknown module invariant %r" % (name,))


# Module kinds with a Plan-30 Part-C worked-steps hook (quiverlab.trace.modules):
# when the report is requested (artifacts.pdf), a module compute also emits the
# exhaustive step events, so the HTML/.tex bundle covers modules like HH does.
_MODULE_TRACE_KINDS = frozenset({
    "projective_resolution", "injective_resolution", "ext", "tau", "tau_minus",
})


def _wants_trace():
    req = _state.get("request") or {}
    return bool((req.get("artifacts") or {}).get("pdf"))


def _emit_module_trace(name, top):
    """Append the worked-step events for a traceable module compute into
    ``_state['events']`` (only when the report is requested). Best-effort: the
    compute already succeeded, so a trace re-run failure silently skips the report
    enhancement rather than break the (already-rendered) result."""
    if name not in _MODULE_TRACE_KINDS or not _wants_trace():
        return
    A = _state["algebra"]
    try:
        from quiverlab.trace import modules as tm
        M = _build_module(A, _state.get("module"), "M")
        if name == "projective_resolution":
            events, _ = tm.trace_projective_resolution(M, top)
        elif name == "injective_resolution":
            events, _ = tm.trace_injective_resolution(M, top)
        elif name == "ext":
            N = _build_module(A, _state.get("ext_target"), "N")
            events, _ = tm.trace_ext(A, M, N, top)
        else:                                       # tau / tau_minus
            events, _ = tm.trace_tau(M, kind=name)
        if _state["events"] is not None:
            _state["events"].extend(events)
    except Exception:                               # best-effort report adornment
        pass


def compute_one(spec):
    """Run ONE Plan-09 compute string against the built algebra."""
    A = _state["algebra"]
    try:
        if A is None:
            raise RequestError("no algebra built (run_build first)")
        name, top = _parse_compute(spec)
        if name in _MODULE_KINDS:
            block = _module_block(name, top)
            _emit_module_trace(name, top)
        elif name in ("hh_cohomology", "hh_homology"):
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            method = (A.hochschild_cohomology if name == "hh_cohomology"
                      else A.hochschild_homology)
            coeff = _build_coefficient(A, _state.get("coefficients"))   # Plan 52 (None=regular)
            hh_kwargs = {"coefficients": coeff} if coeff is not None else {}
            table = method(top, verbose=False, trace=_state["events"], **hh_kwargs)
            keys = list(table.references)
            if coeff is not None:              # Plan 52 coefficient provenance
                for ckey in ("chaparro_schroll_solotar", "lindell_rubio_relative"):
                    if ckey not in keys:
                        keys.append(ckey)
            block = {"kind": table.kind, "top": top, "dims": list(table.dims),
                     "engine": table.engine,
                     "citations": _citation_pairs(keys)}
            if coeff is not None:
                block["coefficients"] = coeff.describe()
            # Plan 35 wave 3d: capture the explicit HH^n / HH_n representatives alongside
            # the dims (basis_classes / chain_basis / differentials / inner_dims per
            # degree) from the SAME dims path -- key-for-key identical to the server twin
            # (quiverlab.hpc.spec._dispatch), so the cross-runner contract holds. None
            # (dims-only) when no representative route applies. Skipped with a coefficient
            # (reps are for the regular bimodule).
            if coeff is None:
                from quiverlab.hochschild.hh_reps import hh_reps_blocks
                try:                           # reps are ADDITIVE + best-effort: a
                    reps = hh_reps_blocks(A, name, top, list(table.dims), table.engine)
                except Exception:              # capture must NEVER break the dims block
                    reps = None
                if reps:
                    block.update(reps)
        elif name == "cyclic_homology":
            # Plan-35 follow-up: cyclic homology HC_0..HC_n (Connes (b, B) mixed
            # complex). Range kind; the block is key-for-key identical to the server
            # twin (quiverlab.hpc.spec._dispatch) -- kind/top/dims/engine/references/
            # citations -- so the cross-runner contract holds.
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            # Plan 35 wave 3b: capture the explicit HC representatives alongside the
            # dims (basis_classes / chain_basis / differentials / column_structure) from
            # the SAME (b, B) total complex -- key-for-key identical to the server twin
            # (quiverlab.hpc.spec._dispatch), so the cross-runner contract holds.
            table, reps = A.cyclic_homology(top, with_reps=True)
            keys = ["cyclic"]
            block = {"kind": table.kind, "top": top, "dims": list(table.dims),
                     "engine": table.engine, "references": keys,
                     "citations": _citation_pairs(keys)}
            block.update(reps)
        elif name == "ss_hochschild":
            # Hochschild (b, B) spectral sequence (Plan 42). Range kind on the algebra
            # block; byte-identical to the server twin (quiverlab.hpc.spec._dispatch) --
            # SAME shared builder (specseq.block.specseq_block) + references->citations.
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            from quiverlab.specseq.block import specseq_block
            block = specseq_block(A, top)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "radical_filtration_ss":
            # Radical-filtration (associated-graded) spectral sequence (P42 preset,
            # wave 2). Range kind on the algebra block; byte-identical to the server
            # twin (quiverlab.hpc.spec._dispatch) -- SAME shared builder
            # (specseq.block.radical_filtration_ss_block) + references->citations.
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            from quiverlab.specseq.block import radical_filtration_ss_block
            block = radical_filtration_ss_block(A, top)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "ar_quiver":
            # AR quiver (P41, wave 2): an ALGEBRA-level BUDGET kind (not a degree
            # range). Byte-identical to the server twin (quiverlab.hpc.spec._dispatch):
            # SAME shared builder (modules.ar.ar_quiver_block) + references->citations.
            from quiverlab.modules.ar import ar_quiver_block
            block = ar_quiver_block(A, budget=top if top is not None else 512)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "split_extension":
            # Split-extension LES (Plan 72 / R5): an ALGEBRA-level TOP-DEGREE budget
            # kind. Byte-identical to the server twin (quiverlab.hpc.spec._dispatch):
            # SAME shared builder (split_extension.split_extension_block), which folds
            # its own loud refusals into status='unsupported' + error.
            from quiverlab.hochschild.split_extension import split_extension_block
            block = split_extension_block(A, top=top if top is not None else 6)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "arrow_removal":
            # Certified arrow removal (Plan 72 / R6): an ALGEBRA-level TOP-DEGREE budget
            # kind. Byte-identical to the server twin (quiverlab.hpc.spec._dispatch):
            # SAME shared builder (arrow_removal.arrow_removal_block).
            from quiverlab.hochschild.arrow_removal import arrow_removal_block
            block = arrow_removal_block(A, top=top if top is not None else 6)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "skew_group_hh":
            # Skew-group HH decomposition (Plan 74 / R8): an ALGEBRA-level TOP-DEGREE
            # budget kind on a SkewGroupAlgebra input. Byte-identical to the server twin
            # (quiverlab.hpc.spec._dispatch): SAME shared builder
            # (hochschild.skew_group.skew_group_hh_block); the Stefan conjugacy-class
            # decomposition is cross-checked against the DIRECT engine, the modular case
            # reports direct-only, a bad action -> clean error field (never a crash).
            from quiverlab.hochschild.skew_group import skew_group_hh_block
            block = skew_group_hh_block(A, top if top is not None else 3)
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "han_transport":
            # Han transport (Plan 73 / R7): an ALGEBRA-level CERTIFICATE kind carrying
            # the new-arrow subset F (in _state["new_arrows"], the derived_compare
            # second-input precedent). Byte-identical to the server twin
            # (quiverlab.hpc.spec._dispatch): SAME shared han.han_transport_block.
            na = _state.get("new_arrows")
            if not na:
                raise RequestError("han_transport needs a 'new_arrows' field (the "
                                   "subset F of arrows whose removal from A defines B)")
            from quiverlab.invariants.han import han_transport_block
            block = han_transport_block(A, list(na), hh_top=top)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "radical_filtration":
            # The radical filtration of mod A (Plan 57 / R37): an ALGEBRA-level BUDGET
            # kind. Byte-identical to the server twin (quiverlab.hpc.spec._dispatch):
            # SAME shared builder (modules.radical.radical_filtration_block). A
            # char-scope refusal is caught into an `error` field, never a crash.
            # NOTE: distinct from radical_filtration_ss (the Loewy radical-series
            # spectral sequence, a DIFFERENT object).
            from quiverlab.errors import QuiverlabError
            from quiverlab.modules.radical import radical_filtration_block
            try:
                block = radical_filtration_block(A, budget=top if top is not None else 512)
            except QuiverlabError as exc:
                block = {"kind": "radical_filtration", "error": str(exc)}
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "ar_invariants":
            # The AR-component invariants (Plan 57 / R21): Liu degrees, partition,
            # directing, rep-directed recognizer. Same BUDGET-kind contract + shared
            # builder (modules.ar_invariants.ar_invariants_block); byte-identical twin.
            from quiverlab.errors import QuiverlabError
            from quiverlab.modules.ar_invariants import ar_invariants_block
            try:
                block = ar_invariants_block(A, budget=top if top is not None else 512)
            except QuiverlabError as exc:
                block = {"kind": "ar_invariants", "error": str(exc)}
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "left_right_parts":
            # Left/right parts (P55, wave 2): an ALGEBRA-level BUDGET kind (not a degree
            # range). Byte-identical to the server twin (quiverlab.hpc.spec._dispatch):
            # SAME shared builder (modules.left_right.left_right_parts_block) +
            # references->citations.
            from quiverlab.modules.left_right import left_right_parts_block
            block = left_right_parts_block(A, budget=top if top is not None else 256)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "tilted_check":
            # Tilted-algebra recognizer (Plan 60): an ALGEBRA-level KNIT-BUDGET kind (not a
            # degree range). Byte-identical to the server twin (quiverlab.hpc.spec._dispatch):
            # SAME shared builder (modules.tilted.tilted_check_block) + references->citations.
            # budget_sections keeps its internal default 4096.
            from quiverlab.modules.tilted import tilted_check_block
            block = tilted_check_block(A, budget_modules=top if top is not None else 256)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "recognizer_ladder":
            # The recognizer ladder (P61, wave 2): an ALGEBRA-level BUDGET kind (not a degree
            # range). Byte-identical to the server twin (quiverlab.hpc.spec._dispatch): SAME
            # shared builder (modules.recognizers_ladder.recognizer_ladder_block) +
            # references->citations.
            from quiverlab.modules.recognizers_ladder import recognizer_ladder_block
            block = recognizer_ladder_block(A, budget=top if top is not None else 256)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "cartan":
            # PER-INVARIANT citation keys, matching the server twin
            # (quiverlab.hpc.spec._dispatch) BYTE-FOR-BYTE. NEVER A.citations() here:
            # that set ACCUMULATES the HH bar-resolution key across a run, so the
            # Cartan matrix was citing Hochschild 1945 instead of ASS2006 (the same
            # bug the server fixed for its report; Marco's report-example.pdf).
            m = A.cartan_matrix()
            keys = ["assem_book"]
            block = {"matrix": m, "latex": _latex_matrix(m),
                     "references": keys, "citations": _citation_pairs(keys)}
        elif name == "coxeter_polynomial":
            import sympy
            p = A.coxeter_polynomial()
            keys = ["lenzing_delapena_spectral", "assem_book"]
            block = {"latex": sympy.latex(p.as_expr()), "text": str(p.as_expr()),
                     "references": keys, "citations": _citation_pairs(keys)}
        elif name == "global_dimension":
            g = A.global_dimension()
            keys = ["assem_book"]
            block = {"text": str(g), "exact": bool(g.exact), "value": g.value,
                     "references": keys, "citations": _citation_pairs(keys)}
        elif name == "homological_profile":
            # The C6 homological-dimension family (Plan 40). ONE shared library
            # builder (modules.homdims.homological_profile) drives this Pyodide twin
            # and the server (quiverlab.hpc.spec._dispatch), byte-identical by
            # construction; we only add the resolved citation pairs from `references`
            # (the cross-runner contract asserts key-for-key equality with the server).
            from quiverlab.modules.homdims import homological_profile
            block = homological_profile(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "fractional_cy":
            # Fractional Calabi-Yau dimension of the stable category (Plan 53 / R24): an
            # ALGEBRA-level scalar kind (schema v1). Byte-identical to the server twin
            # (quiverlab.hpc.spec._dispatch): SAME shared builder
            # (modules.fractional_cy.fractional_cy_block) + references->citations. A
            # non-self-injective input is caught INSIDE the builder into {"error": ...}
            # (no `references` on that shape), never a crash.
            from quiverlab.modules.fractional_cy import fractional_cy_block
            block = fractional_cy_block(A)
            if "references" in block:
                block["citations"] = _citation_pairs(block["references"])
        elif name == "center":
            dim_z, basis = A.center()
            # Basis entries are exact ints/rationals (sympy MPQ over CC) — not
            # JSON-serializable; ship them as exact strings. Per-invariant citation
            # keys matching the server twin (Z(A) = HH^0(A), Hochschild's paper);
            # not A.citations() (see the cartan note).
            keys = ["bar"]
            block = {"dim": dim_z,
                     "basis": [[str(x) for x in row] for row in basis],
                     "references": keys, "citations": _citation_pairs(keys)}
        elif name == "dimension":
            # Parity with the server/HPC runner (quiverlab.hpc.spec._dispatch),
            # which serves `dimension` = A.dim; same `value` semantics, GUI block
            # shape (value + references + citations). Per-invariant keys, not
            # A.citations() (which cited Hochschild 1945 here -- see the cartan note).
            keys = ["assem_book"]
            block = {"value": A.dim, "references": keys,
                     "citations": _citation_pairs(keys)}
        elif name == "ext_algebra":
            # Yoneda / Ext-algebra + Koszulity (Plan 38). Byte-identical to the
            # server twin (quiverlab.hpc.spec._dispatch): SAME library block
            # builder, and citations resolved from `references` the same way.
            from quiverlab.modules.ext_algebra import ext_algebra_block
            block = ext_algebra_block(A, top if top is not None else 6)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "koszul":
            # Generalized Koszulity (Plan 77): Berger N-Koszul, Cassidy-Shelton K2
            # with an explicit certified window, Brenner-Butler-King
            # (p,q)-almost-Koszul, Herscovich's scoped multi-Koszul, and the INTERNAL
            # (path-length) generation degrees of Ext(k,k). Byte-identical to the
            # server twin (quiverlab.hpc.spec._dispatch): SAME library block builder
            # (modules.nkoszul.koszul_profile_block), same default top, same
            # `references`->citations, and the SAME QuiverlabError catch so a
            # presentation-less algebra yields a typed error block, never a traceback.
            from quiverlab.errors import QuiverlabError
            from quiverlab.modules.nkoszul import koszul_profile_block
            try:
                block = koszul_profile_block(A, top if top is not None else 8)
            except QuiverlabError as exc:
                block = {"kind": "koszul", "error": str(exc)}
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "recognizers":
            # Recognizer batch + type detection (Plan 38). Byte-identical to the
            # server twin: SAME library block builder + `references`->citations.
            from quiverlab.invariants.recognizers import recognizers_block
            block = recognizers_block(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "coxeter_spectral":
            # Certified Coxeter spectral analysis (Plan 58 / R20). Byte-identical to
            # the server twin (quiverlab.hpc.spec._dispatch): SAME library block
            # builder (invariants.coxeter_spectral.coxeter_spectral_block) +
            # `references`->citations.
            from quiverlab.invariants.coxeter_spectral import coxeter_spectral_block
            block = coxeter_spectral_block(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "derived_fingerprint":
            # Derived fingerprint (Plan 43). Byte-identical to the server twin
            # (quiverlab.hpc.spec._dispatch): SAME library block builder
            # (derived.block.derived_fingerprint_block) + `references`->citations.
            from quiverlab.derived.block import derived_fingerprint_block
            block = derived_fingerprint_block(A, top if top is not None else 4)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "derived_compare":
            # Two-algebra derived-fingerprint comparison (P43, wave 2): needs the SECOND
            # algebra B, carried in _state["algebra_b"] (built lazily here, like the
            # module blocks). Byte-identical to the server twin
            # (quiverlab.hpc.spec._dispatch): SAME shared builder
            # (derived.block.derived_compare_block) + references->citations.
            spec_b = _state.get("algebra_b")
            if spec_b is None:
                raise RequestError("derived_compare needs a second algebra 'algebra_b' "
                                   "(the B in the derived-fingerprint comparison)")
            B = _algebra_from_spec(spec_b)
            from quiverlab.derived.block import derived_compare_block
            block = derived_compare_block(A, B, top if top is not None else 4)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "strings":
            # Gentle / string subsystem (Plan 46): census + bands + rep-type + AG.
            # Byte-identical to the server twin (quiverlab.hpc.spec._dispatch): SAME
            # library block builder + `references`->citations.
            from quiverlab.strings.block import strings_block
            block = strings_block(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "quasi_hereditary":
            # Quasi-hereditary structure (Plan 47): an algebra-scalar kind (schema v1).
            # Byte-identical to the server twin (quiverlab.hpc.spec._dispatch): SAME
            # library block builder (modules.quasihereditary.quasi_hereditary_block,
            # natural order + order-dependence note) + `references`->citations.
            from quiverlab.modules.quasihereditary import quasi_hereditary_block
            block = quasi_hereditary_block(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "fundamental_group":
            # pi1(Q, I) + abelianization (Plan 56): an algebra-scalar kind (schema v1).
            # Byte-identical to the server twin (quiverlab.hpc.spec._dispatch): SAME
            # shared builder (invariants.coverings_block.fundamental_group_block) +
            # `references`->citations.
            from quiverlab.invariants.coverings_block import fundamental_group_block
            block = fundamental_group_block(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "simply_connected":
            # Three-valued simple connectivity + the R16 strongly-simply-connected
            # certificate (Plan 56). Byte-identical to the server twin
            # (quiverlab.hpc.spec._dispatch): SAME shared builder
            # (invariants.coverings_block.simply_connected_block) + references->citations.
            from quiverlab.invariants.coverings_block import simply_connected_block
            block = simply_connected_block(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name in ("cup", "cap", "bracket", "connes_b", "bv_operator"):
            # HH product surface (Plan 35) + the BV operator (Plan 54: bv_operator).
            # Each library method returns a frozen result whose .blocks() IS the block
            # dict (kind/top/engine + tables|matrices + references); we only add the
            # resolved citation pairs, exactly as the server twin does
            # (quiverlab.hpc.spec._dispatch). The block keeps `references` -- the
            # cross-runner contract asserts key-for-key equality with the server.
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            method = {"cup": A.cup_products, "cap": A.cap_products,
                      "bracket": A.gerstenhaber_brackets,
                      "connes_b": A.connes_differentials,
                      "bv_operator": A.bv_operator}[name]
            block = method(top).blocks()
            block["citations"] = _citation_pairs(block["references"])
        elif name == "hh_lie_module":
            # HH^* as a graded Lie module over HH^1 (Plan 71 / R12): a top-carrying HH
            # kind (the gerstenhaber_brackets precedent). Shared block builder ships
            # hh_dims, the module-axiom / inner-zero verdicts, the char-0 weight table
            # and the indecomposable-summand decomposition table (byte-identical twin of
            # quiverlab.hpc.spec._dispatch). Default top = 2 when no range is given.
            from quiverlab.hochschild.lie_module import hh_lie_module_block
            block = hh_lie_module_block(A, top if top is not None else 2)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "tau_tilting":
            # C4 tau-tilting engine (Plan 45): algebra-level, budget (not degree). SAME
            # shared library builder (tautilting.block.tau_tilting_block) + references ->
            # citations as the server twin (quiverlab.hpc.spec._dispatch), so the
            # cross-runner contract holds byte-for-byte. Honest budget cap block.
            from quiverlab.tautilting.block import tau_tilting_block
            block = tau_tilting_block(A, budget=top if top is not None else 512)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "congruences":
            # Torsion-lattice congruences (Plan 64): algebra-level, pair budget (not degree).
            # SAME shared library builder (tautilting.congruence.congruences_block) +
            # references -> citations as the server twin (quiverlab.hpc.spec._dispatch), so the
            # cross-runner contract holds byte-for-byte -- INCLUDING the char-caveat error path:
            # a QuiverlabError refusal (rigorous over char 0 / char > dim) is caught into the
            # SAME {"kind","error","references"} shape spec.py returns (the silting-branch
            # pattern), so the two runners' error paths are byte-identical; a non-QuiverlabError
            # bug surfaces loudly (fail-fast). Honest complete-iff block.
            from quiverlab.tautilting.congruence import _CITATIONS as _CONG_KEYS
            from quiverlab.tautilting.congruence import congruences_block
            try:
                block = congruences_block(A, budget=top if top is not None else 512)
            except quiverlab.QuiverlabError as exc:
                block = {"kind": "congruences", "error": str(exc),
                         "references": list(_CONG_KEYS)}
            block["citations"] = _citation_pairs(block["references"])
        elif name == "tau_cluster":
            # tau-cluster morphism category W(A) + picture group (Plan 66 / R29):
            # algebra-level, pair budget (not degree). SAME shared library builder
            # (tautilting.cluster_morphism.tau_cluster_block) + references -> citations as the
            # server twin (quiverlab.hpc.spec._dispatch), so the cross-runner contract holds
            # byte-for-byte -- INCLUDING the char-caveat / tau-tilting-infinite error path: a
            # QuiverlabError refusal is caught into the SAME {"kind","error","references"} shape
            # spec.py returns (the congruences-branch pattern). Objects = #wide (ties P64), the
            # Hanson-Igusa classifying-space cube complex + K(pi,1) verdict + the picture group,
            # certified complete iff A is tau-tilting-finite. Honest complete-iff block.
            from quiverlab.tautilting.cluster_morphism import _REFERENCES as _TCL_KEYS
            from quiverlab.tautilting.cluster_morphism import tau_cluster_block
            try:
                block = tau_cluster_block(A, budget=top if top is not None else 512)
            except quiverlab.QuiverlabError as exc:
                block = {"kind": "tau_cluster", "error": str(exc),
                         "references": list(_TCL_KEYS)}
            block["citations"] = _citation_pairs(block["references"])
        elif name == "hh1_lie":
            # HH^1 as a Lie algebra (Plan 70 / R11): algebra-level, DIM budget (not
            # degree). SAME shared library builder (invariants.hh1_lie.hh1_lie_block) +
            # references -> citations as the server twin (quiverlab.hpc.spec._dispatch),
            # so the cross-runner contract holds byte-for-byte. Der/Inn + bracket + series
            # + solvable/nilpotent over any exact field; char-0 radical/Levi/sl2-count.
            from quiverlab.invariants.hh1_lie import hh1_lie_block
            block = hh1_lie_block(A, budget=top if top is not None else 48)
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "deformations":
            # Formal deformations / L-infinity / Maurer-Cartan (Plan 78 / R13): algebra-level,
            # DIM budget (not degree), the plan's OWN DEFORM_MAXDIM (32), NOT P70's 48 -- the
            # cost class is different (CS bracket on HH^2/HH^3, not a Der solve). SAME shared
            # library builder (hochschild.deformations.deformations_block) + references ->
            # citations as the server twin (quiverlab.hpc.spec._dispatch), byte-for-byte.
            from quiverlab.hochschild.deformations import (
                DEFORM_MAXDIM, deformations_block)
            block = deformations_block(A, budget=top if top is not None else DEFORM_MAXDIM)
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "wall_chamber":
            # Wall-and-chamber structure via bricks (Plan 63 / R25): algebra-level, budget
            # (not degree). SAME shared library builder
            # (tautilting.wallchamber.wall_chamber_structure) + references -> citations as
            # the server twin (quiverlab.hpc.spec._dispatch), byte-for-byte. Certified
            # complete iff brick-finite <=> tau-tilting-finite, else an honest bounded region
            # (status='budget', no count). The brick/is_isomorphic char caveat -> error block.
            from quiverlab.tautilting.wallchamber import wall_chamber_structure
            try:
                block = wall_chamber_structure(A, budget=top if top is not None else 512)
            except quiverlab.QuiverlabError as exc:
                block = {"kind": "wall_chamber", "error": str(exc)}
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "cluster_category":
            # The Amiot-Keller cluster category, certified acyclic slice (Plan 79 / R31).
            # SAME shared library builder (cluster.category.cluster_category_block) and the
            # same references -> citations as the server twin, byte-for-byte.
            from quiverlab.cluster.category import cluster_category_block
            try:
                block = cluster_category_block(A, budget_pairs=top if top is not None
                                               else 512)
            except quiverlab.QuiverlabError as exc:
                block = {"kind": "cluster_category", "error": str(exc)}
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "silting":
            # Silting theory (Plan 67 / Aihara-Iyama): algebra-level, RADIUS,BUDGET pair
            # (not a degree). SAME shared library builder (derived.block.silting_block) +
            # references -> citations as the server twin (quiverlab.hpc.spec._dispatch), so
            # the cross-runner contract holds byte-for-byte. A QuiverlabError refusal (the
            # char-scope / presentation / verifier-edge path) is caught into an `error`
            # field; a non-QuiverlabError bug is NOT swallowed here -- it surfaces loudly
            # (the fail-fast house rule), so this narrows to "the typed refusals never
            # crash the block", not "never a crash".
            radius, budget = top if top is not None else (3, 64)
            from quiverlab.derived.block import silting_block
            try:
                block = silting_block(A, radius=radius, budget=budget)
            except quiverlab.QuiverlabError as exc:
                block = {"kind": "silting", "error": str(exc)}
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "exceptional_sequences":
            # Exceptional sequences (Plan 65 / R27+R28): algebra-level, enumeration budget
            # (not a degree). SAME shared library builder
            # (tautilting.exceptional.exceptional_sequences_block) + references -> citations
            # as the server twin (quiverlab.hpc.spec._dispatch), so the cross-runner contract
            # holds byte-for-byte. Classical (if hereditary) + tau counts, honest status.
            from quiverlab.tautilting.exceptional import exceptional_sequences_block
            block = exceptional_sequences_block(
                A, budget=top if top is not None else 4096)   # sane DoS cap (Plan 65 H-3)
            block["citations"] = _citation_pairs(block.get("references", []))
        elif name == "string_homological":
            # Homological string-algebra test (Plan 59 / R34, Suarez-Alvarez). Byte-
            # identical to the server twin (quiverlab.hpc.spec._dispatch): SAME library
            # block builder (modules.string_homological.string_homological_block) +
            # `references`->citations. A rep-infinite / self-injective / presentation-
            # less input returns an {"error": ...} block, never a raise.
            from quiverlab.modules.string_homological import string_homological_block
            block = string_homological_block(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "toupie":
            # Toupie structure (Plan 59 / R35). Byte-identical to the server twin
            # (quiverlab.hpc.spec._dispatch): SAME library block builder
            # (families.toupie.toupie_block) -- recognizer + branch/direct-arrow counts +
            # HH + char-0 sl_a lower bound -- + `references`->citations.
            from quiverlab.families.toupie import toupie_block
            block = toupie_block(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "tame_wild":
            # Tits-form tame/wild certificate (Plan 62 / R19). Byte-identical to the
            # server twin (quiverlab.hpc.spec._dispatch): SAME library block builder
            # (invariants.tits_block.tame_wild_block) -- combinatorial Tits form + weak
            # positivity/nonnegativity + the rep-finite/tame/wild verdict gated on the
            # P56 certificate over char 0 -- + `references`->citations. A presentation-
            # less / non-triangular input returns an {"error": ...} block, never a raise.
            from quiverlab.invariants.tits_block import tame_wild_block
            block = tame_wild_block(A)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "skew_gentle":
            # Skew-gentle world (Plan 68 / R32). Byte-identical to the server twin
            # (quiverlab.hpc.spec._dispatch): the SAME library block builder
            # (skewgentle.block.skew_gentle_block) reads the triple off the split
            # algebra's construction marker -- recognizer + split shape + dim law +
            # classification counts + support tau-tilting + rep-type certificate --
            # + `references`->citations.
            from quiverlab.skewgentle.block import skew_gentle_block
            block = skew_gentle_block(A, budget=(top or 512))
            block["citations"] = _citation_pairs(block["references"])
        elif name == "incidence_cohomology":
            # HH^* of an incidence algebra through the order complex (Plan 75 / R9).
            # Byte-identical to the server twin (quiverlab.hpc.spec._dispatch): SAME
            # library block builder (hochschild.simplicial.incidence_cohomology_block)
            # -- the Gerstenhaber-Schack/Cibils theorem HH^n(kP) = H^n(Delta(P); k) via
            # the exact integer SNF simplicial engine, plus the face vector, the
            # integral torsion and the char-dependence verdict -- + `references` ->
            # citations. A DEGREE-RANGE kind. An algebra without poset provenance
            # returns an {"error": ...} block, never a raise.
            from quiverlab.hochschild.simplicial import incidence_cohomology_block
            if top is None:
                raise RequestError(
                    "incidence_cohomology needs a degree range, e.g. "
                    "'incidence_cohomology:0..3'")
            block = incidence_cohomology_block(A, top)
            block["citations"] = _citation_pairs(block["references"])
        elif name == "tate_hochschild":
            # Tate-Hochschild (singular Hochschild) cohomology (Plan 76 / R3).
            # Byte-identical to the server twin (quiverlab.hpc.spec._dispatch): SAME
            # library block builder (hochschild.tate.tate_hochschild_block) -- HHhat^m
            # for every m in [-top, top], NEGATIVE degrees included, off a complete
            # resolution of A over A^e -- + `references` -> citations. A DEGREE-RANGE
            # kind. Degrees a route cannot know are null, never a wrong number
            # (HHhat^0 is the stable centre, NOT HH^0, so it is native-only). A
            # non-Gorenstein algebra, or the DEFERRED Gorenstein-but-not-self-injective
            # native request, returns an {"error": ...} block, never a raise.
            from quiverlab.hochschild.tate import tate_hochschild_block
            if top is None:
                raise RequestError(
                    "tate_hochschild needs a degree range, e.g. "
                    "'tate_hochschild:0..3'")
            from quiverlab.errors import QuiverlabError as _QErr
            try:
                block = tate_hochschild_block(A, top)
            except _QErr as exc:
                block = {"kind": "tate_hochschild", "error": str(exc),
                         "references": list(_TATE_REFERENCES)}
            block["citations"] = _citation_pairs(block["references"])
        else:
            raise RequestError("unknown invariant %r" % (name,))
        _state["results"].append(dict(block, invariant=spec))
        out = {"ok": True, "invariant": spec, "block": block}
    except Exception as exc:
        out = _fail(exc)
        out["invariant"] = spec
        # The downloadable bundle must not silently omit failures.
        # (No bundle exists until run_build resets results to a list.)
        if _state["results"] is not None:
            _state["results"].append({"invariant": spec, "error": out["error"]})
    return json.dumps(out)


def _named_modules():
    """The modules the request named, as ``(label, Module)`` for the report's "The
    modules" section (mirrors ``quiverlab.hpc.spec._named``). Best-effort: a module
    that no longer builds is skipped -- the report describes, it does not verify."""
    A = _state.get("algebra")
    if A is None:
        return []
    specs = [("M", _state.get("module")), ("N", _state.get("ext_target")),
             ("N", _tor_target_spec())]
    named, seen_n = [], 0
    for label, mspec in specs:
        if not isinstance(mspec, dict):
            continue
        if label == "N":
            seen_n += 1
        try:
            named.append((label, _build_module(A, mspec, label)))
        except Exception:
            continue
    if seen_n > 1:                        # both an Ext and a Tor target: tell them apart
        roles = iter(("N (Ext target)", "N (Tor target)"))
        named = [(next(roles) if lbl == "N" else lbl, m) for lbl, m in named]
    return named


def trace_html():
    """The worked-steps report as an HTML string.

    Carries the session's COMPUTED RESULT BLOCKS as well as the worked steps, so the
    saved report is everything the page showed (Marco 2026-07-29) -- and is therefore
    produced even when nothing recorded worked steps (Cartan + centre, say).

    '' when the request did not ask for a report (``artifacts.pdf``), or when there
    is nothing at all to report."""
    events = _state["events"] or []
    results = _state["results"] or []
    if not events and not (results and _wants_trace()):
        return ""
    from quiverlab.trace.provenance import references_for, resolve_references
    from quiverlab.trace.render_html import render_html
    from quiverlab.trace.json_guide import build_json_guide
    title = "Worked steps — %s" % repr(_state["algebra"]).splitlines()[0]
    return render_html(list(events), title=title, algebra=_state["algebra"],
                       references=resolve_references(references_for(events)),
                       results=results, modules=_named_modules(),
                       json_guide=build_json_guide(results))


def trace_json():
    """The worked-steps report as the JSON machine record ('' when nothing was traced).

    Plan 34: the third mandated artifact (PDF/HTML + this) -- the complete event
    stream, deterministic and schema-versioned. A GUI download button hooks this
    accessor beside trace_tex()/trace_html(). Uses the SAME (events, title,
    references) inputs as the other two renderers, so the byte-for-byte machine
    record matches what the writer persists for the same computation."""
    events = _state["events"] or []
    if not events:
        return ""
    from quiverlab.trace.provenance import references_for, resolve_references
    from quiverlab.trace.render_json import render_json
    title = "Worked steps — %s" % repr(_state["algebra"]).splitlines()[0]
    return render_json(list(events), title=title, algebra=_state["algebra"],
                       references=resolve_references(references_for(events)))


def tikz():
    return "" if _state["algebra"] is None else _state["algebra"].tikz()


# The P44/P46/P48 construction families whose real constructors take
# Algebra / Module / BrauerGraph / Triangulation arguments (so the generic
# ql.<Family>(field=..., **params) reproduce form crashes). Kept in sync (by name)
# with quiverlab.hpc.spec._SYNTHETIC_FAMILY_PARAMS and webapp/server/catalog.py.
_SYNTHETIC_FAMILY_NAMES = frozenset({
    "BrauerGraphAlgebra", "OnePointExtension", "CornerAlgebra",
    "OppositeAlgebra", "MarkedSurface"})


def _synthetic_reproduce_lines(family, params, field_ref):
    """Runnable construction lines for one P44/P46/P48 construction family, mirroring
    ``quiverlab.hpc.spec._build_synthetic`` EXACTLY (the emitted script builds a
    byte-identical algebra). ``field_ref`` is a source expression for the field
    ('ql.CC', 'ql.GF(5)' or 'QQ'); the caller has already emitted ``import quiverlab
    as ql`` and, for QQ, its import. BYTE-FOR-BYTE identical to the server twin
    (quiverlab.hpc.spec._synthetic_reproduce_lines)."""
    if family == "OppositeAlgebra":
        return [f"A = ql.PathAlgebra({params.get('base')!r}, field={field_ref})",
                "A = A.opposite()"]
    if family == "CornerAlgebra":
        return [f"A = ql.PathAlgebra({params.get('base')!r}, field={field_ref})",
                f"A = A.corner_algebra({list(params.get('vertices', []))!r})"]
    if family == "OnePointExtension":
        kind = params.get("module_kind", "simple")
        # M0 (not M) so a request module block named M is never clobbered.
        return [f"A = ql.PathAlgebra({params.get('base')!r}, field={field_ref})",
                f"M0 = A.{kind}({params.get('module_vertex')!r})",
                "A = ql.OnePointExtension(A, M0)"]
    if family == "MarkedSurface":
        preset = params.get("preset")
        if preset == "annulus_C22":
            return ["from quiverlab.surfaces import annulus_triangulation, jacobian_of",
                    "T = annulus_triangulation(2, 2)",
                    f"A = jacobian_of(T, field={field_ref})"]
        if preset == "hexagon_internal":
            return ["from quiverlab.surfaces import "
                    "hexagon_with_internal_triangle, jacobian_of",
                    "T = hexagon_with_internal_triangle()",
                    f"A = jacobian_of(T, field={field_ref})"]
        return ["from quiverlab.surfaces import fan_triangulation, jacobian_of",
                "T = fan_triangulation(6)",       # disc fan -> gentle A3
                f"A = jacobian_of(T, field={field_ref})"]
    if family == "BrauerGraphAlgebra":
        # Reconstruct the un-flattened BrauerGraph args (the cosmetic ribbon end-tags
        # are synthesized exactly as _build_brauer does, so the graph is identical).
        edges = tuple((e[0], e[1]) for e in params.get("edges", []))
        cyclic = {vtx: tuple((int(k), f"{vtx}:{pos}") for pos, k in enumerate(idxs))
                  for vtx, idxs in params.get("cyclic_order", [])}
        mult = {vtx: int(m) for vtx, m in params.get("multiplicities", [])}
        return ["from quiverlab.families.brauer import "
                "BrauerGraph, BrauerGraphAlgebra",
                f"G = BrauerGraph(edges={edges!r}, cyclic_order={cyclic!r})",
                f"A = BrauerGraphAlgebra(G, {mult!r}, field={field_ref})"]
    if family == "SkewGroupAlgebra":
        # Rebuild A|xG from the flattened params via the shared library builder (the
        # base quiver + the explicit group action). Byte-identical to the server twin.
        return ["from quiverlab.families.skew_group import build_skew_group_from_params",
                f"A = build_skew_group_from_params({dict(params)!r}, {field_ref})"]
    return [f"# built via the webapp family {family!r}; see docs"]


def python_snippet():
    """Copy-paste reproduction of the GUI computation (the GUI-to-library bridge)."""
    req = _state["request"]
    if req is None:
        return ""
    alg = req["algebra"]
    f = alg["field"]
    if f["kind"] == "CC":
        field_name, field_expr = "CC", "CC"
    elif f["kind"] == "QQ":
        field_name, field_expr = "QQ", "QQ"
    else:
        q = f["p"] ** f.get("n", 1)
        field_name, field_expr = "GF", "GF(%d)" % q
    # QQ is not a top-level export, so it needs its own import line.
    header = ("from quiverlab import Quiver\nfrom quiverlab.fields import QQ"
              if field_name == "QQ" else "from quiverlab import Quiver, %s" % field_name)
    if alg.get("kind") == "family":
        # A `family` block (server tier for most families -- the twin BUILDS only
        # SkewGentle / SkewGroup / Incidence and refuses the rest, but the
        # reproduce string mirrors quiverlab.hpc.spec._snippet's family branch so the
        # cross-runner snippet contract holds). Real constructors, never a crashing
        # ql.<Family>(field=..., **params) form for the five construction families.
        field_ref = ("ql.CC" if f["kind"] == "CC"
                     else "QQ" if f["kind"] == "QQ" else "ql.%s" % field_expr)
        head = ["import quiverlab as ql"]
        if f["kind"] == "QQ":
            head.append("from quiverlab.fields import QQ")
        family, params = alg.get("family"), alg.get("params", {})
        if family in _SYNTHETIC_FAMILY_NAMES:
            lines = head + _synthetic_reproduce_lines(family, params, field_ref)
        else:
            ptxt = "".join(", %s=%r" % (k, v) for k, v in params.items())
            lines = head + ["A = ql.%s(field=%s%s)" % (family, field_ref, ptxt)]
        lines.append("print(A.dim)")
    else:
        arrows = ", ".join('"%s": (%d, %d)' % (k, s, t)
                           for k, (s, t) in alg["arrows"].items())
        q_line = "Q = Quiver(vertices=%r, arrows={%s})" % (alg["vertices"], arrows)
        potential = alg.get("potential")
        if isinstance(potential, str) and potential.strip():
            lines = [header,
                     "from quiverlab import Potential, JacobianAlgebra",
                     "from quiverlab.combinat.relations import _split_terms, _parse_term",
                     "", q_line,
                     "W = Potential(Q, [_parse_term(t, Q) for t in _split_terms(%r)])"
                     % potential,
                     "A = JacobianAlgebra(Q, W, field=%s)" % field_expr,
                     "print(A.dim)"]
        else:
            lines = [header, "", q_line,
                     "A = Q.algebra(relations=%r, field=%s)"
                     % (list(alg.get("relations", [])), field_expr),
                     "print(A.dim)"]
    if req.get("module"):
        lines += _module_snippet_lines(req["module"], "M")
    if req.get("ext_target"):
        lines += _module_snippet_lines(req["ext_target"], "N")
    if req.get("tor_target"):
        lines += _module_snippet_lines(_tor_target_spec(), "N")
    calls = {"hh_cohomology": "A.hochschild_cohomology(%d)",
             "hh_homology": "A.hochschild_homology(%d)",
             "cyclic_homology": "A.cyclic_homology(%d)",
             "ss_hochschild": "A.hochschild_bB_ss(%d)",
             "cartan": "A.cartan_matrix()", "coxeter_polynomial": "A.coxeter_polynomial()",
             "global_dimension": "A.global_dimension()", "center": "A.center()",
             # `dimension` is a scalar invariant compute_one serves (A.dim) -- it MUST
             # have a snippet entry or python_snippet() KeyErrors when it is requested.
             "dimension": "A.dim",
             # Plan 38: Koszulity / Ext-algebra + the recognizer batch. Both are
             # scalar kinds (no %d) so the reproduce snippet never leaves a literal.
             "ext_algebra": "A.ext_algebra()",
             # Plan 77: the generalized-Koszulity profile, a scalar kind (no %d).
             "koszul": "A.koszul_profile()",
             # Plan 79: the cluster-category slice, a scalar kind (no %d).
             "cluster_category": "A.cluster_category()",
             "recognizers": ("[A.is_semisimple(), A.is_hereditary(), A.is_gentle(), "
                             "A.dynkin_type(), A.form_type()]"),
             # Quasi-hereditary structure (Plan 47): a scalar kind, no %d.
             "quasi_hereditary": "A.is_quasi_hereditary()",
             # Plan 59 recognizer batteries: two scalar algebra-only kinds, no %d.
             "string_homological": ("homological_string_test(A)  "
                                    "# from quiverlab.modules.string_homological"),
             "toupie": "(is_toupie(A), toupie_block(A))  # from quiverlab.families.toupie",
             # pi1 + simple connectivity (Plan 56): scalar kinds, no %d.
             "fundamental_group": "A.fundamental_group()",
             "simply_connected": "A.is_simply_connected()",
             # Tits-form tame/wild certificate (Plan 62 / R19): a scalar kind, no %d.
             "tame_wild": "A.tame_wild_certificate()",
             # Plan 75 / R9: the order-complex route to HH^* of an incidence algebra --
             # a degree-range kind (%d = top). Needs the IncidenceAlgebra provenance.
             "incidence_cohomology": "A.incidence_cohomology(%d)",
             # Plan 76 / R3: Tate-Hochschild -- the block reports HHhat^{-top..top}
             # (negative degrees included). A degree-range kind (%d = top).
             "tate_hochschild": "A.tate_hochschild(%d)",
             # Derived fingerprint (Plan 43): a scalar kind, no %d (top defaults to 4).
             "derived_fingerprint": "derived_fingerprint(A)  # from quiverlab.derived",
             # HH product surface (Plan 35): same four calls as the server snippet
             # map (quiverlab.hpc.spec._snippet); each needs a range (%d = top).
             "cup": "A.cup_products(%d)", "cap": "A.cap_products(%d)",
             "bracket": "A.gerstenhaber_brackets(%d)",
             "connes_b": "A.connes_differentials(%d)",
             "bv_operator": "A.bv_operator(%d)",
             # Plan 71: HH^* as a Lie module over HH^1, a top-carrying HH kind (%d = top).
             "hh_lie_module": "A.hh_lie_module(top=%d)",
             # Plan 45: the C4 tau-tilting kind carries a pair budget (%d = budget_pairs).
             "tau_tilting": "A.exchange_graph(budget_pairs=%d)",
             # Plan 64: the congruences kind carries a pair budget (%d = budget).
             "congruences": "A.congruence_lattice(budget=%d)",
             # Plan 66: the tau_cluster kind carries a pair budget (%d = budget); the block
             # also builds A.picture_group(budget=...) with the same budget.
             "tau_cluster": "A.tau_cluster_category(budget=%d)",
             # Plan 70: HH^1 as a Lie algebra, a scalar algebra-only kind, no %d.
             "hh1_lie": "A.hh1_lie_structure()",
             # Plan 78: deformations is a scalar algebra-only kind, no %d.
             "deformations": "A.deformation_structure()",
             # Plan 63: the wall-and-chamber kind carries a pair budget (%d = budget_pairs).
             "wall_chamber": "A.wall_chamber_structure(budget_pairs=%d)",
             # Plan 72: split_extension / arrow_removal carry a top-degree budget (%d = top).
             "split_extension": "A.split_extension_cohomology(%d)",
             "arrow_removal": "A.arrow_removal(top=%d)",
             # Plan 74: skew_group_hh carries a top-degree budget (%d = top). A is the
             # built A|xG; stefan_decomposition takes its base + action.
             "skew_group_hh": ("stefan_decomposition(base, action, %d)  "
                               "# from quiverlab.hochschild.skew_group"),
             # Plan 73: han_transport carries the new-arrow subset F (no %d).
             "han_transport": "A.han_transport(F)  # F = the new-arrow subset defining B",
             # Plan 67: silting carries a RADIUS,BUDGET pair (top = (radius, budget) tuple;
             # tmpl % top fills both %d).
             "silting": "A.silting_exploration(radius=%d, budget=%d)",
             # Plan 65: exceptional_sequences (classical + tau), a scalar kind, no %d.
             "exceptional_sequences": ("(A.exceptional_sequences(), "
                                       "A.tau_exceptional_sequences(want_sequences=False))"),
             # Plan 57: radical_filtration + ar_invariants carry a module budget.
             "radical_filtration": "A.radical_filtration(budget_modules=%d)",
             "ar_invariants": "A.ar_invariants(budget_modules=%d)",
             # Plan 55: the left/right parts kind carries a module budget (%d = budget).
             "left_right_parts": "A.left_right_parts(budget=%d)",
             # Plan 60: the tilted recognizer carries the knit budget (%d = budget_modules).
             "tilted_check": "A.tilted_check(budget_modules=%d)",
             # Plan 61: the recognizer ladder carries a module budget (%d = budget).
             "recognizer_ladder": "A.recognizer_ladder(budget=%d)",
             "dimension_vector": "M.dimension_vector()",
             "rad_top_soc": "(M.radical(), M.top(), M.socle())",
             "tau": "M.tau()", "tau_minus": "M.tau_minus()",
             "ext": "[A.ext(M, N, i) for i in range(%d + 1)]",
             "tor": "tor_dims(A, M, N, %d)  # from quiverlab.modules.tor",
             "decompose": "M.decompose()",
             "barcode": "A.barcode(M)  # from quiverlab.modules.barcode import barcode",
             "almost_split": "M.almost_split_sequence()",
             "projective_resolution": "M.projective_resolution(%d).dimension_vectors()",
             "injective_resolution": "M.injective_resolution(%d).dimension_vectors()",
             "projective_dimension": "M.projective_resolution(%d).pd()" % _PD_BOUND,
             "injective_dimension": "M.injective_dimension()"}
    for spec in req.get("compute", []):
        name, top = _parse_compute(spec)
        tmpl = calls[name]
        call = tmpl % top if ("%d" in tmpl and top is not None) else tmpl
        lines.append("print(%s)" % call)
    return "\n".join(lines) + "\n"


def _fmt_scalar(x):
    if isinstance(x, Fraction):
        return (str(x.numerator) if x.denominator == 1
                else "Fraction(%d, %d)" % (x.numerator, x.denominator))
    return str(x)


def _pymat(mat):
    return ("[" + ", ".join("[" + ", ".join(_fmt_scalar(x) for x in row) + "]"
                            for row in mat) + "]")


def _module_snippet_lines(mspec, varname):
    """Runnable lines rebuilding `varname`: a builtin call, or the FULL exact
    arrow-action matrices A.module consumes."""
    b = mspec.get("builtin")
    if b is not None:
        return ['%s = A.%s(%r, side="%s")'
                % (varname, b.get("kind"), b.get("vertex"), mspec.get("side", "right"))]
    dimvec, action = _full_matrices(_state["algebra"], mspec)
    lines = []
    if any(isinstance(x, Fraction) for m in action.values() for row in m for x in row):
        lines.append("from fractions import Fraction")
    dv_lit = "{" + ", ".join("%r: %d" % (v, n) for v, n in dimvec.items()) + "}"
    maps_lit = "{" + ", ".join('"%s": %s' % (a, _pymat(m)) for a, m in action.items()) + "}"
    lines.append('%s = A.module(%s, %s, side="%s", name="%s")'
                 % (varname, dv_lit, maps_lit, mspec.get("side", "right"), varname))
    return lines


def result_bundle():
    # Marco 2026-07-31 (ADDENDUM 2): mirror the spec envelope's json_guide so the
    # Pyodide result bundle documents how to recover every computed object too.
    from quiverlab.trace.json_guide import build_json_guide
    results = _state["results"] or []
    return json.dumps({"schema": SCHEMA_VERSION, "request": _state["request"],
                       "quiverlab_version": quiverlab.__version__,
                       "results": results,
                       "json_guide": build_json_guide(results)}, indent=1)


# --- Plan 11: wait-time estimation (pure arithmetic; spec 2026-07-22) --------

# Fitted on THIS machine (native, numba BLOCKED to match Pyodide's pure path),
# 2026-07-22, scripts/fit_eta_model.py. Worst-off factor on heavy (>0.3 s)
# grid points: bar 1.46x, fast 3.76x — inside one bucket width; the in-flight
# rescale absorbs the rest. Units are calibrated native-seconds; the browser
# factor comes from calibrate(). Do not hand-tune: rerun the fit script.
ETA_MODEL = {
    "bar":  {"alpha": 1.4622e-07, "p": 1.3},
    "fast": {"alpha": 5.3447e-07, "p": 1.1},
    "scalars": {"cartan": 0.01, "coxeter_polynomial": 0.2,
                # Plan 58: coxeter_spectral is bimodal-but-fast -- real-dominant certifies
                # sub-second (minpoly + Sturm interval), complex-dominant is REFUSED
                # without computing (the deterministic _real_roots_suffice gate, never the
                # measured 121 s minimal_polynomial hang), so 0.5 is honest in both branches.
                "coxeter_spectral": 0.5,
                "center": 0.05, "global_dimension": 0.5,
                # Plan 40: the C6 family aggregates gl.dim + finitistic + dominant +
                # Gorenstein + Igusa-Todorov (several resolutions), so a bit heavier.
                # Plan 53 added phidim/psidim (an AR knit) + LIT to the same block.
                "homological_profile": 2.0,
                # Plan 53: fractional_cy iterates nu/Omega on the simples with a bounded
                # (m,ell) search -- a few small syzygy/Nakayama passes.
                "fractional_cy": 3.0,
                # module kinds (Plan 26): cheap dim-vector reads up to
                # resolution/dimension probes that build syzygies to depth.
                "dimension_vector": 0.02, "rad_top_soc": 0.05,
                "tau": 0.1, "tau_minus": 0.1, "ext": 0.2, "tor": 0.2,
                "decompose": 0.3, "almost_split": 0.3,
                # Plan 69: barcode = decompose (cheap A_n/zigzag) OR a full AR knit
                # (CL). The knit-heavy CL case is pushed off the instant tier by the
                # estimator's classify() upgrade (reason="knit_heavy"), not this weight.
                "barcode": 0.3,
                "projective_resolution": 0.2, "injective_resolution": 0.2,
                "projective_dimension": 0.3, "injective_dimension": 0.3,
                # Plan 44 / 49: single-module homological probes (tilting_check =
                # self-Ext vanishing + summand count; orbit_geometry = Voigt
                # rigidity + Kac canonical decomposition). Same cost class as the
                # other module resolutions/probes -- WITHOUT these keys they fell
                # through to the 0.1 default and were silently under-estimated.
                "tilting_check": 0.3, "orbit_geometry": 0.3,
                # Plan 38: ext_algebra walks a resolution + Yoneda products;
                # recognizers is cheap structural combinatorics + a reduction system.
                "ext_algebra": 2.0, "recognizers": 0.1,
                # Plan 77: koszul runs ext_algebra AND resolves every simple again for
                # the internal degrees, so it sits just above ext_algebra.
                "koszul": 2.2,
                # Plan 79: cluster_category runs the AR knit AND the tau-tilting exchange
                # graph BFS -- the BFS dominates and is NOT algebra-dim-bound (E6 has
                # dim kQ 36 but 833 cluster-tilting objects), so it sits well above the
                # module-level kinds.
                "cluster_category": 6.0,
                # Plan 42: the (b, B) spectral sequence builds the same exponential
                # bar (b, B) bicomplex cyclic homology uses, plus the page algebra.
                "ss_hochschild": 2.0,
                # Plan 43: derived_fingerprint runs Cartan/Coxeter + HH/HC to top=4
                # (the HH pass dominates; cyclic may fall back or error honestly).
                "derived_fingerprint": 1.0,
                # Plan 46: strings = a bounded string/band DFS + AG walk on the
                # reduction system; a touch heavier than recognizers.
                "strings": 0.2,
                # Plan 47: quasi_hereditary builds Delta/Nabla + a gl.dim check +
                # the greedy Delta-peel of each P(v); a few small resolutions.
                "quasi_hereditary": 0.5,
                # Plan 56: fundamental_group = a reduction system + block linear
                # algebra + a small ZZ SNF; simply_connected additionally runs the
                # convex-subset separation sweep (decompose per vertex per subset).
                "fundamental_group": 0.5, "simply_connected": 2.0,
                # Plan 45: the C4 tau-tilting engine BFSes the exchange graph via the
                # 2-term silting mutation (per-pair K^b Hom + minimal approximations);
                # heavier than the string DFS, budget-capped honestly.
                "tau_tilting": 2.0,
                # Plan 64: congruences BFSes the exchange graph (as tau_tilting) and then
                # runs the principal-congruence fixed points + the kappa/CLO build -- a bit
                # heavier than tau_tilting alone.
                "congruences": 3.0,
                # Plan 66: tau_cluster BFSes the exchange graph, builds the wide poset AND a
                # reduction / sub-g-fan per object (the closed-star enumeration + the
                # picture-group Ext scan) -- heavier than congruences.
                "tau_cluster": 4.0,
                # Plan 70: hh1_lie runs the Der/Inn Leibniz null space (d^2 unknowns /
                # d^3 equations, ~ d^5.4 over QQ) + the bracket/series/Killing; budget-
                # capped honestly at dim 48. The same cost class as tau_tilting.
                "hh1_lie": 2.0,
                # Plan 78: deformations runs HH^2 + HH^3 and then the CS Gerstenhaber
                # bracket [alpha, alpha] per basis direction -- the bracket, not the dim,
                # is the driver, and it tracks HH-RICHNESS (a dim-8 HH-rich algebra can
                # cost more than a dim-20 HH-thin one). Weighted well above hh1_lie.
                "deformations": 20.0,
                # Plan 63: wall_chamber runs the tau_tilting exchange-graph BFS PLUS the
                # per-brick submodule enumeration for each D(B) -- just above tau_tilting.
                "wall_chamber": 2.5,
                # Plan 72: split_extension runs the CS Hom-complex of L = T(B)
                # (dim 2*dim B) to degree top+2 and block-partitions it for the snake;
                # arrow_removal runs HH_* + HH^* of A and B. Both are CS-dominated,
                # around the tau_tilting cost class.
                "split_extension": 2.0, "arrow_removal": 1.5,
                # Plan 74: skew_group_hh runs the twisted-bar HH of the base A per
                # conjugacy class + the Z(g)-transport + a DIRECT HH of A|xG (dim
                # |G|*dim A). Bar-dominated on the |G|*dim A algebra, tau_tilting class.
                "skew_group_hh": 2.0,
                # Plan 73: han_transport = the bounded-extension legs (tensor powers +
                # module pd + gl.dim) plus HH_*(A), HH_*(B); resolution-dominated.
                "han_transport": 3.0,
                # Plan 67: silting = a bounded-radius BFS of the silting quiver via K^b
                # Hom + minimal approximations + cone/reduce per step; the hyper-Hom passes
                # dominate. Budget-capped honestly (complete only for local).
                "silting": 1.5,
                # Plan 55: left/right parts = an AR knit + the N^2 Hom predecessor matrix +
                # a pd/id sweep + the two support-algebra End certificates; knit-dominated,
                # the same cost class as tau_tilting.
                "left_right_parts": 2.0,
                # Plan 60: tilted_check = an AR knit + a budget-capped transversal search
                # (faithful + Hom(X,tauY)=0 + tilting/presented-End certificate per candidate);
                # knit- and certificate-dominated, the same cost class as left_right_parts.
                "tilted_check": 2.0,
                # Plan 61: the recognizer ladder reads the P55 atlas + gl.dim + a second AR
                # knit (weakly-shod SCC) + HH^1 (ada/Theorem B); a touch heavier than P55.
                "recognizer_ladder": 2.5,
                # Plan 59: string_homological KNITS the AR quiver + realizes/decomposes
                # extensions (expensive, ar_quiver class); toupie is a small HH + a
                # graph-shape scan (cheap).
                "string_homological": 2.0, "toupie": 0.5,
                # Plan 62: tame_wild = the Tits form (P56 minimal-relation counts) +
                # weak positivity/nonnegativity (box/PSD/list) + the P56 simple/strong-
                # simple-connectivity convex sweep -- the convex sweep dominates
                # (simply_connected class), sized above the cheap scalars.
                "tame_wild": 3.0,
                # Plan 75: incidence_cohomology enumerates the CHAINS of the poset and
                # runs one integer Smith normal form -- combinatorics on a complex far
                # smaller than the enveloping algebra the general HH engines resolve, so
                # it sizes BELOW every hh_* route (this is the whole point of the fast
                # path). Sized with the cheap structural scalars.
                "incidence_cohomology": 0.3,
                # Plan 76: tate_hochschild builds the minimal A^e-resolution to
                # top+2 AND splices its Nakayama-twisted dual, then collapses BOTH
                # halves -- roughly three times the work of the hh_cohomology route
                # it shares the resolution with, and it also runs the periodicity
                # certificate. Sized above hh_cohomology accordingly.
                "tate_hochschild": 3.0,
                # Plan 65: exceptional_sequences KNITS the AR quiver + runs the braid-orbit
                # BFS (classical) and the exchange-graph BFS (tau); knit- and BFS-dominated,
                # the heavier ar_quiver/tame_wild cost class.
                "exceptional_sequences": 3.0},
}
_MAX_CELLS = 4_000_000        # the library's bar guard (frozen contract)
_BUCKETS = (                  # (upper bound in seconds, id, label)
    (15.0, "seconds", "estimated: a few seconds"),
    (75.0, "minute", "estimated: under a minute"),
    (360.0, "minutes", "estimated: a few minutes"),
    (None, "long", "estimated: could be long — Cancel anytime"),
)


def _bar_guard_cells(m, n):
    """rows*cols of the bar coboundary d^n guard: (m(m-1)^{n+1}) * (m(m-1)^n)."""
    return (m * (m - 1) ** (n + 1)) * (m * (m - 1) ** n)


def _cap_degree(m, top):
    """First degree in 0..top whose guard exceeds the library's max_cells."""
    if m <= 2:
        return None               # (m-1) <= 1: sizes stay tiny forever
    for n in range(top + 1):
        if _bar_guard_cells(m, n) > _MAX_CELLS:
            return n
    return None


def _hh_units(m, top, route):
    mdl = ETA_MODEL[route]
    return mdl["alpha"] * sum(_bar_guard_cells(m, n) ** mdl["p"]
                              for n in range(top + 1))


def _units_for(dim, field_spec, compute):
    """(total units, per-invariant breakdown, cap info or None)."""
    route = ("fast" if field_spec.get("kind") == "GF"
             and field_spec.get("n", 1) == 1 else "bar")
    total, breakdown, cap = 0.0, [], None
    for spec in compute:
        name, top = _parse_compute(spec)
        if name in ("hh_cohomology", "hh_homology"):
            k = _cap_degree(dim, top)
            if k is not None:
                # The homology engine's own DepthLimitError names b_{k+1} —
                # boundary indexing is shifted by one vs the coboundary guard —
                # so the SHOWN degree gets +1 for hh_homology (units keep raw k).
                shown = k + 1 if name == "hh_homology" else k
                if cap is None or shown < cap["degree"]:
                    cap = {"degree": shown, "invariant": spec}
            u = _hh_units(dim, min(top, (k - 1) if k is not None else top), route)
        else:
            u = ETA_MODEL["scalars"].get(name, 0.1)
        total += u
        breakdown.append({"invariant": spec, "units": u})
    return total, breakdown, cap


def bucket_for_seconds(seconds):
    for bound, bid, label in _BUCKETS:
        if bound is None or seconds < bound:
            return json.dumps({"bucket": bid, "label": label})
    raise AssertionError("unreachable")


def estimate(factor):
    """Estimate the CURRENT request against the just-built algebra."""
    try:
        A, req = _state["algebra"], _state["request"]
        if A is None or req is None:
            raise RequestError("no algebra built (run_build first)")
        units, breakdown, cap = _units_for(
            A.dim, req["algebra"]["field"], req.get("compute", []))
        seconds = units * float(factor)
        if cap is not None:
            bucket, label = "cap", ("will hit the engine's cell cap near "
                                    "degree %d" % cap["degree"])
        else:
            b = json.loads(bucket_for_seconds(seconds))
            bucket, label = b["bucket"], b["label"]
        out = {"ok": True, "dim": A.dim, "units": units, "seconds": seconds,
               "bucket": bucket, "label": label,
               "cap_degree": cap["degree"] if cap else None,
               "breakdown": breakdown}
    except Exception as exc:
        out = _fail(exc)
    return json.dumps(out)


_CAL_FIELD = {"kind": "GF", "p": 2, "n": 1}
_CAL_COMPUTE = ["hh_cohomology:0..6", "cartan", "center"]


def calibrate():
    """Time a fixed workload; factor = seconds per model unit on THIS machine.
    Builds locally (never via _state) so a visitor's probe state survives."""
    import time
    Q = quiverlab.Quiver(vertices=[1], arrows={"x": (1, 1)})
    t0 = time.monotonic()
    A = Q.algebra(relations=["x*x*x"], field=quiverlab.GF(2))
    A.hochschild_cohomology(6, verbose=False)
    A.cartan_matrix()
    A.center()
    seconds = time.monotonic() - t0
    units, _, _ = _units_for(3, _CAL_FIELD, _CAL_COMPUTE)
    return json.dumps({"seconds": seconds, "units": units,
                       "factor": seconds / units})
