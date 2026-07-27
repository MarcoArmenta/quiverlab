"""Plan 28 -- the ``quiverlab-hpc`` command-line interface (argparse only; no
click). Verbs: ``run``, ``render``, ``sample-config``, ``estimate``, ``version``,
``selftest``, ``gui``.

Progress goes to stderr; stdout stays clean (so ``sample-config`` /
``estimate`` output can be piped). Thread caps (``NUMBA_NUM_THREADS`` /
``OMP_NUM_THREADS``) are set from the detected usable cores unless already set
(``$SLURM_CPUS_PER_TASK`` still wins). Exit codes are BSD sysexits:

  * 0  -- done
  * 75 -- clean checkpoint stop (deepen), resumable (sbatch requeues)
  * 65 -- bad config / relation / data error
  * 64 -- usage error
  * 73 -- cannot write the output
  * 70 -- internal error

The ``gui`` verb imports ``webapp`` LAZILY inside its handler, so importing
``quiverlab.hpc`` never pulls the web stack (pinned by the import-boundary test).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from quiverlab.hpc import resources
from quiverlab.hpc.spec import (
    RESULT_SCHEMA, ComputeError, SpecError, CheckpointStop, load_config,
    parse_request, build_algebra, parse_compute_item,
)

# BSD sysexits (see plan doc "Exit codes").
EX_OK = 0
EX_USAGE = 64
EX_DATAERR = 65
EX_SOFTWARE = 70
EX_CANTCREAT = 73
EX_CHECKPOINT = 75

# Error types that are honest data/config refusals (-> exit 65); anything else is
# an unexpected internal failure (-> exit 70).
_DATA_ERROR_TAGS = frozenset({
    "SchemaError", "CatalogError", "FieldError", "ResultTooLarge",
    "DuplicateComputeItem", "ExactnessError", "RelationError",
    "AdmissibilityError", "NotFiniteDimensionalError", "DepthLimitError",
    "CitationError",
})


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #

class _Parser(argparse.ArgumentParser):
    """ArgumentParser that exits 64 (EX_USAGE) instead of argparse's default 2."""

    def error(self, message):
        self.print_usage(sys.stderr)
        self.exit(EX_USAGE, f"{self.prog}: error: {message}\n")


def _eprint(*args):
    print(*args, file=sys.stderr, flush=True)


def _set_thread_caps(env=None):
    env = os.environ if env is None else env
    cap = str(resources.default_thread_cap(env))
    for var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS"):
        env.setdefault(var, cap)


def _stderr_progress(ev: dict):
    if "deepen" in ev:
        _eprint(ev["deepen"])
    else:
        _eprint(f"[{ev['step'] + 1}/{ev['of']}] {ev['kind']}")


def _compute_error_code(error_type: str) -> int:
    return EX_DATAERR if error_type in _DATA_ERROR_TAGS else EX_SOFTWARE


def _max_degree(items) -> int:
    hi = 0
    for it in items:
        if it.hi is not None:
            hi = max(hi, it.hi)
    return hi


# --------------------------------------------------------------------------- #
# Annotated sample config (parseable + validating YAML)
# --------------------------------------------------------------------------- #

SAMPLE_CONFIG = """\
# quiverlab-hpc run config (Plan 28).
# Usage:  quiverlab-hpc run this-file.yaml -o result.json
#         quiverlab-hpc render result.json -o report.html
schema: 1                 # request schema (1, or 2 to add a `module` block)
algebra:
  kind: quiver            # "quiver" (below) or "family" (a named builder)
  vertices: [1]           # vertex labels (integers)
  arrows:
    x: [1, 1]             # arrow name -> [source, target]
  relations: ["x*x*x"]    # opaque relation strings (k[x]/(x^3) here)
  field:
    kind: GF              # "GF" (finite field) or "CC"
    p: 2                  # the prime (GF only)
    n: 1                  # GF(p^n); 1 for a prime field
compute:                  # each item is "kind" or "kind:lo..hi"
  - hh_cohomology:0..4    # HH^0..HH^4
  - cartan                # the Cartan matrix
artifacts:
  pdf: false              # worked-steps report (trace_steps.html + trace.json)
  tikz: false             # the quiver as tikz.tex
hpc:                      # optional batch/checkpoint knobs (CLI-only)
  checkpoint_dir: null    # a dir on $SCRATCH -> resumable big hh_homology
  time_limit_s: null      # wall budget; a clean checkpoint stop exits 75
  max_mem_bytes: null     # transient-memory guard for the deepen path
  prime: null             # field for checkpointed HH_* (null = the algebra's own)
"""


# --------------------------------------------------------------------------- #
# Verb: run
# --------------------------------------------------------------------------- #

def _cmd_run(args) -> int:
    try:
        cfg = load_config(args.config)
    except SpecError as exc:
        _eprint(f"config error: {exc}")
        return EX_DATAERR
    # CLI overrides fold into the hpc block.
    hpc = dict(cfg.get("hpc") or {})
    if args.checkpoint_dir is not None:
        hpc["checkpoint_dir"] = args.checkpoint_dir
    if args.time_limit is not None:
        hpc["time_limit_s"] = args.time_limit
    if args.max_mem is not None:
        hpc["max_mem_bytes"] = args.max_mem
    if args.prime is not None:
        hpc["prime"] = args.prime
    if args.allow_large:
        hpc["allow_large"] = True
    # Default the deepen transient-memory guard to 4/5 of detected RAM when the
    # user set neither a config value nor --max-mem (integer arithmetic only).
    if hpc.get("checkpoint_dir") and hpc.get("max_mem_bytes") is None:
        mem = resources.detect_mem_bytes()
        if mem is not None:
            hpc["max_mem_bytes"] = mem * 4 // 5
    if hpc:
        cfg["hpc"] = hpc

    out = Path(args.output)
    from quiverlab.hpc.spec import run as _run
    try:
        result = _run(cfg, out.parent, progress_cb=_stderr_progress,
                      result_schema=RESULT_SCHEMA, write_result=False)
    except CheckpointStop as exc:
        _eprint(f"checkpoint stop: {exc}")
        _eprint("resumable: rerun the same command (deepen resumes from the "
                "checkpoint); under SLURM, sbatch requeues.")
        return EX_CHECKPOINT
    except SpecError as exc:
        _eprint(f"config error: {exc}")
        return EX_DATAERR
    except ComputeError as exc:
        _eprint(f"{exc.error_type}: {exc.message}")
        return _compute_error_code(exc.error_type)
    except OSError as exc:                     # e.g. artifact dir not creatable
        _eprint(f"cannot write output: {exc}")
        return EX_CANTCREAT
    try:
        out.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    except OSError as exc:
        _eprint(f"cannot write {out}: {exc}")
        return EX_CANTCREAT
    _eprint(f"wrote {out}")
    return EX_OK


# --------------------------------------------------------------------------- #
# Verb: render
# --------------------------------------------------------------------------- #

def _cmd_render(args) -> int:
    from quiverlab.hpc import report
    try:
        out, actual = report.render(args.result, args.output, fmt=args.format,
                                    on_warn=_eprint)
    except report.ResultSchemaError as exc:
        _eprint(f"cannot render: {exc}")
        return EX_DATAERR
    except report.ReportWriteError as exc:
        _eprint(str(exc))
        return EX_CANTCREAT
    except report.ReportError as exc:
        _eprint(f"cannot render: {exc}")
        return EX_DATAERR
    _eprint(f"wrote {out} ({actual})")
    return EX_OK


# --------------------------------------------------------------------------- #
# Verb: sample-config
# --------------------------------------------------------------------------- #

def _cmd_sample_config(args) -> int:
    sys.stdout.write(SAMPLE_CONFIG)
    return EX_OK


# --------------------------------------------------------------------------- #
# Verb: estimate
# --------------------------------------------------------------------------- #

_FIELD_MULT = {"GF": 1, "CC": 50}
_OPS_PER_MINUTE = 500_000_000


def _estimate_ops(dim, max_degree, field_kind):
    """Try the webapp estimator (if importable), else the identical wheel-side
    formula. Kept lazy so the base wheel works with no webapp present."""
    try:
        from webapp.server import estimator as _wa       # noqa
        return _wa.estimate_ops(dim, max_degree, field_kind)
    except Exception:
        return (dim ** 3) * (max_degree + 1) * _FIELD_MULT.get(field_kind, 50)


def _tier_of(ops, max_deg):
    if ops <= 2_000_000 and max_deg <= 8:
        return "instant"
    if ops <= 500_000_000 and max_deg <= 20:
        return "queued"
    if ops <= 50_000_000_000 and max_deg <= 40:
        return "big"
    return "beyond-big"


def _cmd_estimate(args) -> int:
    try:
        cfg = load_config(args.config)
        req = parse_request(cfg)
        A = build_algebra(req.algebra)
    except SpecError as exc:
        _eprint(f"config error: {exc}")
        return EX_DATAERR
    except ComputeError as exc:
        _eprint(f"{exc.error_type}: {exc.message}")
        return _compute_error_code(exc.error_type)
    items = [parse_compute_item(s) for s in req.compute]
    max_deg = _max_degree(items)
    dim = A.dim
    ops = _estimate_ops(dim, max_deg, req.algebra.field.kind)
    minutes = max(1, -(-ops // _OPS_PER_MINUTE))
    tier = _tier_of(ops, max_deg)
    # Suggested SLURM knobs: 2x the time hint, and an order-of-magnitude memory
    # figure (all integer). A generous floor keeps tiny jobs from asking for 0.
    suggest_time_s = minutes * 60 * 2
    suggest_mem_bytes = max((dim * dim) * (max_deg + 1) * 64, 512 * 1024 ** 2)
    res = resources.detect_resources()

    print("tier:            %s" % tier)
    print("algebra dim:     %d" % dim)
    print("max degree:      %d" % max_deg)
    print("est. cells:      %d" % ops)
    print("est. minutes:    %d" % minutes)
    print("suggest --time:  %d s (%d min)" % (suggest_time_s, -(-suggest_time_s // 60)))
    print("suggest --mem:   %d bytes (%d MiB)" % (suggest_mem_bytes,
                                                  suggest_mem_bytes // (1024 ** 2)))
    print(resources.format_resources(res))
    host_mem = res.get("mem_bytes")
    if host_mem is not None and suggest_mem_bytes > host_mem:
        print("WARNING: the suggested memory exceeds detected host RAM; this host "
              "may be too small -- use a bigger node or narrow the degree range.")
    return EX_OK


# --------------------------------------------------------------------------- #
# Verb: version
# --------------------------------------------------------------------------- #

def _cmd_version(args) -> int:
    import quiverlab
    print("quiverlab %s" % getattr(quiverlab, "__version__", "unknown"))
    print("quiverlab-hpc CLI (result_schema %d)" % RESULT_SCHEMA)
    print(resources.format_resources(resources.detect_resources()))
    return EX_OK


# --------------------------------------------------------------------------- #
# Verb: selftest
# --------------------------------------------------------------------------- #

def _cmd_selftest(args) -> int:
    import tempfile
    from quiverlab.hpc.spec import run as _run
    from quiverlab.hpc import report
    cfg = {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                    "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": ["hh_cohomology:0..3", "cartan", "dimension"],
        "artifacts": {"pdf": False, "tikz": False},
    }
    try:
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            result = _run(cfg, d, result_schema=RESULT_SCHEMA)
            dims = result["results"]["hh_cohomology"]["dims"]
            assert dims == [3, 2, 2, 2], dims
            assert result["results"]["cartan"]["matrix"] == [[3]]
            assert result["results"]["dimension"]["value"] == 3
            assert result["result_schema"] == RESULT_SCHEMA
            out, fmt = report.render(result, d / "report.txt", fmt="txt")
            text = out.read_text(encoding="utf-8")
            assert "quiverlab report" in text and "HH^" in text, text[:200]
    except Exception as exc:  # any failure is an internal error
        _eprint(f"selftest FAILED: {type(exc).__name__}: {exc}")
        return EX_SOFTWARE
    print("selftest OK: k[x]/(x^3) over GF(2) HH^0..3 = [3, 2, 2, 2]; report rendered")
    print(resources.format_resources(resources.detect_resources()))
    return EX_OK


# --------------------------------------------------------------------------- #
# Verb: gui (offline webapp; lazy [web] import)
# --------------------------------------------------------------------------- #

def _cmd_gui(args) -> int:
    try:
        from webapp.server.offline import serve_offline
    except ImportError as exc:
        _eprint("the offline GUI needs the web extras and the webapp package: "
                "install quiverlab[web,hpc] (and run from a source checkout / the "
                "container image). Details: %s" % exc)
        return EX_SOFTWARE
    data_dir = args.data_dir or os.environ.get("QUIVERLAB_DATA")
    try:
        serve_offline(port=args.port, data_dir=data_dir, open_hint=not args.no_open)
    except KeyboardInterrupt:
        return EX_OK
    return EX_OK


# --------------------------------------------------------------------------- #
# Parser
# --------------------------------------------------------------------------- #

def _build_parser() -> _Parser:
    p = _Parser(prog="quiverlab-hpc",
                description="quiverlab HPC / container CLI (Plan 28)")
    sub = p.add_subparsers(dest="verb")

    pr = sub.add_parser("run", help="validate a config, compute, write result.json")
    pr.add_argument("config", help="config file (YAML or JSON); '-' reads stdin")
    pr.add_argument("-o", "--output", default="result.json", help="result JSON path")
    pr.add_argument("--checkpoint-dir", default=None,
                    help="enable resumable big hh_homology via this checkpoint dir")
    pr.add_argument("--time-limit", type=int, default=None,
                    help="wall-time budget in seconds (clean checkpoint stop = exit 75)")
    pr.add_argument("--max-mem", type=int, default=None,
                    help="transient-memory guard in bytes for the deepen path")
    pr.add_argument("--prime", type=int, default=None,
                    help="prime the checkpointed HH_* is computed over (default 32003)")
    pr.add_argument("--allow-large", action="store_true",
                    help="permit large results (relaxes result guards)")
    pr.set_defaults(func=_cmd_run)

    pd = sub.add_parser("render",
                        help="render result.json to HTML/text, or emit trace.json")
    pd.add_argument("result", help="path to result.json")
    pd.add_argument("-o", "--output", default=None, help="report output path")
    pd.add_argument("--format", choices=("auto", "html", "txt", "json"),
                    default="auto",
                    help="output format (default auto = HTML; "
                         "json = worked-steps event stream trace.json)")
    pd.set_defaults(func=_cmd_render)

    ps = sub.add_parser("sample-config", help="print an annotated config YAML")
    ps.set_defaults(func=_cmd_sample_config)

    pe = sub.add_parser("estimate", help="tier + suggested --time/--mem for a config")
    pe.add_argument("config", help="config file (YAML or JSON); '-' reads stdin")
    pe.set_defaults(func=_cmd_estimate)

    pv = sub.add_parser("version", help="print versions + detected host resources")
    pv.set_defaults(func=_cmd_version)

    pt = sub.add_parser("selftest", help="tiny end-to-end compute + render check")
    pt.set_defaults(func=_cmd_selftest)

    pg = sub.add_parser("gui", help="serve the offline local webapp (needs [web])")
    pg.add_argument("--port", type=int, default=8000, help="port (default 8000)")
    pg.add_argument("--data-dir", default=None,
                    help="data directory (default $QUIVERLAB_DATA or ~/.quiverlab)")
    pg.add_argument("--no-open", action="store_true",
                    help="do not print the browser open hint")
    pg.set_defaults(func=_cmd_gui)
    return p


def main(argv=None) -> int:
    _set_thread_caps()
    parser = _build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "func", None) is None:
        parser.print_help(sys.stderr)
        return EX_USAGE
    return args.func(args)


if __name__ == "__main__":     # pragma: no cover
    sys.exit(main())
