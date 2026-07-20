"""Pure per-algebra analysis + scan driver (quiverlab.batch).

Adapted from hanlab's ``labdb.analyze`` / ``_analyze_open`` / ``run_scan`` (Marco
Armenta's HansConjecture, MIT).  The per-algebra computation is a **pure function**
``analyze(spec) -> record``: a ``spec`` is a small JSON-serializable description
``{"builder", "args", "N", "primes"}`` and ``analyze`` builds the algebra from the
builder registry, computes HH_* and HH^* over the requested primes, the apparent
complexities, and the homology/cohomology asymmetry, and returns a JSON record.
``run_scan`` drives a list of specs -- single-worker by default; the same
``analyze`` is the unit of work for a SLURM array on a cluster (one spec per task,
merge the JSON afterwards -- no shared state).

Structural builders route through the exact F_p Hochschild engine
(``Algebra.hochschild_homology`` / ``hochschild_cohomology`` over GF(PRIME)) plus
the apparent-complexity classifier ``scan3.complexity_of``.  The ``reduction_system``
builder routes to the open-zone minimal A^e-resolution path.  The ``scan_open``
policy helpers ``han_verdict`` / ``depth_for_dim`` / ``max_term_dim_for_dim`` are
ported verbatim from hanlab.scan_open below.

All record values are exact ints / strings (dims are ints; complexities are ints
or the sentinel string ``'>=2'``) -- there is no float round-trip anywhere.
"""
from quiverlab.batch.builders import BUILDERS
from quiverlab.fields import GF

PRIME = 32003


# ---------------------------------------------------------------------------
# scan_open policy helpers -- ported verbatim from hanlab.scan_open (MIT)
# ---------------------------------------------------------------------------
def depth_for_dim(d):
    """How deep to push HH_* by dimension (memory- and time-aware on a small host)."""
    if d <= 6:
        return 20
    if d <= 9:
        return 16
    if d <= 12:
        return 12
    return 9


def max_term_dim_for_dim(d):
    # cap on m^2 * r_n (the k-dimension of a single resolution term)
    if d <= 9:
        return 9000
    if d <= 12:
        return 13000
    return 18000


def han_verdict(seq):
    """Classify an HH_* sequence (char-0 proxy) for the Han question.

    Returns one of:
      'nonzero-tail'      : HH_n > 0 for every computed n (consistent with Han);
      'COUNTEREXAMPLE?'   : HH_n = 0 for the last computed n (and some earlier n>0),
                            i.e. homology appears to vanish in high degree -- the
                            shape a Han counterexample would take;
      'all-zero-after'    : trailing run of zeros (subset of the above).
    Plus the growth descriptor: 'bounded' if the tail is eventually constant/period,
    else 'growing'."""
    if not seq:
        return "empty", "?"
    nz = [i for i, v in enumerate(seq) if v != 0]
    last = seq[-1]
    if last == 0:
        verdict = "COUNTEREXAMPLE?"
    else:
        verdict = "nonzero-tail"
    # growth: compare second half max to first-half max
    h = len(seq) // 2
    tail_max = max(seq[h:]) if seq[h:] else 0
    growth = "growing" if (len(seq) >= 6 and seq[-1] > seq[max(0, len(seq) - 4)]) else "bounded"
    return verdict, growth


# ---------------------------------------------------------------------------
# open-zoo -> spec adapter -- ported verbatim from hanlab.cluster.make_shards (MIT)
# ---------------------------------------------------------------------------
def _reduction_system_args(entry):
    """Catalog entry -> the ``reduction_system`` builder's args ``[ngen, rules, name]``.

    Matches ``builders._reduction_system(ngen, rules, name)``; if the builder's arg
    shape ever changes, this is the only line to touch."""
    return [entry["ngen"], entry["rules"], entry["name"]]


def open_zoo_to_specs(catalog, primes=(PRIME,), max_dim=None, min_dim=None, limit=None,
                      max_transient_bytes=None):
    """Convert open-zoo catalog entries to ``reduction_system`` specs, using the same
    depth/budget policy as the scan_open helpers (so a cluster sweep matches the local one).

    ``min_dim``/``max_dim`` keep only algebras with min_dim <= catalog dim <= max_dim (a
    dim band -- the cheap-first curriculum); ``limit`` caps the number of specs (after
    filtering).  ``max_transient_bytes`` (bytes; None = off) is baked into each spec as the
    per-worker peak-memory budget for the minimal-resolution radK transient, so the spec
    list stays a pure function of its inputs (sharded == serial).  When None the key is
    omitted, so existing specs/records are byte-unchanged (back-compat)."""
    specs = []
    for e in catalog:
        d = e["dim"]
        if max_dim is not None and d > max_dim:
            continue
        if min_dim is not None and d < min_dim:
            continue
        spec = {
            "builder": "reduction_system",
            "args": _reduction_system_args(e),
            "N": depth_for_dim(d),
            "primes": list(primes),
            "max_term_dim": max_term_dim_for_dim(d),
        }
        if max_transient_bytes is not None:
            spec["max_transient_bytes"] = int(max_transient_bytes)
        specs.append(spec)
        if limit is not None and len(specs) >= limit:
            break
    return specs


# ---------------------------------------------------------------------------
# The pure unit of work
# ---------------------------------------------------------------------------
def analyze(spec):
    """Pure function spec -> record.  Deterministic; safe to run in any process.

    spec: {"builder": str, "args": list, "N": int, "primes": list[int] (optional)}.
    The "reduction_system" builder routes to the open-zone minimal-resolution path
    (HH_* only, no cohomology); every other builder takes the structural fast-engine
    path (HH_* and HH^* over GF(p))."""
    if spec["builder"] == "reduction_system":
        return _analyze_open(spec)
    from quiverlab.engine.scan3 import complexity_of

    primes = tuple(spec.get("primes", (PRIME, 2, 3)))
    N = spec["N"]
    builder = spec["builder"]
    args = list(spec.get("args", []))
    # rebuild the family over each prime field; the fast engine computes exact F_p HH
    algs = {p: BUILDERS[builder](*args, field=GF(p)) for p in primes}
    alg0 = algs[primes[0]]
    name = spec.get("name") or getattr(alg0, "zoo_name", None) or f"{builder}{tuple(args)}"
    hh = {p: list(map(int, algs[p].hochschild_homology(N, engine="auto").dims))
          for p in primes}
    hc = {p: list(map(int, algs[p].hochschild_cohomology(N, engine="auto").dims))
          for p in primes}
    cx_h = complexity_of(hh[primes[0]])
    cx_c = complexity_of(hc[primes[0]])
    hseq = hh[primes[0]]                          # char-0 proxy drives the flags
    return {
        "name": name,
        "builder": builder,
        "args": args,
        "dim": int(alg0.dim),
        "N": N,
        # every quiverlab family is associative by construction (validated at build)
        "associative": True,
        "HH_homology": {str(p): hh[p] for p in primes},
        "HH_cohomology": {str(p): hc[p] for p in primes},
        "cx_homology": cx_h,
        "cx_cohomology": cx_c,
        "asymmetric": bool(cx_h != cx_c),
        "homology_bounded": bool(cx_h == 0),                # HH_* eventually zero
        "homology_nonzero_tail": bool(hseq[-1] != 0),       # HH_N != 0
    }


def _analyze_open(spec):
    """Open-zone record path: HH_* via the minimal A^e-resolution over F_p.

    The minimal engine has no cohomology, so cohomology-side fields are explicit None
    (`cx_cohomology`, `asymmetric`); the counterexample-query flags are still computed
    from the HH_* tail so `bounded_homology_candidates` / `by_complexity` work on open
    rows.  A non-confluent / infinite / unbuildable reduction system yields a JSON error
    record (never a raise), so a poison spec cannot kill a multiprocessing shard.

    NOTE: this path needs a full open-zone reduction-system -> Algebra build (the Plan-04
    Chouhy-Solotar backend) to reproduce the golden values on non-monomial systems; the
    golden port test is skipped until that backend is on the branch."""
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import (
        minimal_homology_dims, minimal_resolution)
    from quiverlab.engine.scan3 import complexity_of
    from quiverlab.errors import QuiverlabError

    primes = tuple(spec.get("primes", (PRIME,)))
    args = spec["args"]
    name = args[2] if len(args) > 2 else "A"
    rec = {
        "name": name,
        "builder": "reduction_system",
        "args": list(args),
        "kind": "open",                       # marks the open-zone record family
        "N": spec.get("N"),                   # carried even on the error path so the
        # resume-skip / merge-dedup identity (builder, args, N) matches the spec's.
    }
    # build (catch bad-input errors -> error record, never crash a worker)
    try:
        A = BUILDERS["reduction_system"](*args, field=GF(PRIME))
    except (QuiverlabError, ValueError) as e:
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["associative"] = False            # an unbuildable system is not an algebra
        return rec

    dim = int(A.dim)
    N = int(spec.get("N", depth_for_dim(dim)))
    budget = int(spec.get("max_term_dim", max_term_dim_for_dim(dim)))
    # Peak-memory guard (bytes; None = off): bounds the dense radK transient so a high-dim
    # build truncates gracefully instead of OOM-killing the worker.
    mtb = spec.get("max_transient_bytes")
    mtb = int(mtb) if mtb is not None else None
    rec["dim"] = dim
    rec["N"] = N
    rec["max_term_dim"] = budget
    rec["max_transient_bytes"] = mtb
    rec["associative"] = True                 # validated at build

    E = to_engine(A.unit_adapted())           # F_p engine algebra (b_0 = 1_A)
    hh = minimal_homology_dims(E, N, primes=primes, max_term_dim=budget,
                               max_transient_bytes=mtb)
    rks, _cols, _eng, trunc = minimal_resolution(E, N, PRIME, max_term_dim=budget,
                                                 max_transient_bytes=mtb)
    hseq = hh[PRIME]                          # char-0 proxy drives the verdict/flags
    verdict, growth = han_verdict(hseq)

    rec["HH_homology"] = {str(p): list(map(int, hh[p])) for p in primes}
    rec["resolution_ranks"] = [int(rks.get(n, 0)) for n in range(0, N + 2)]
    rec["depth_reached"] = len(hseq) - 1
    rec["truncated_at"] = trunc               # None or int (JSON null / int)
    rec["han_verdict"] = verdict
    rec["growth"] = growth
    # cohomology side does not exist on this engine -> explicit None (not a misleading 0)
    rec["cx_cohomology"] = None
    rec["asymmetric"] = None
    # counterexample-query flags from the HH_* tail (keep the structural queries working)
    rec["cx_homology"] = complexity_of(hseq)
    rec["homology_bounded"] = bool(rec["cx_homology"] == 0)
    rec["homology_nonzero_tail"] = bool(hseq[-1] != 0) if hseq else False
    return rec


def run_scan(specs, n_workers=1):
    """Run `analyze` over a list of specs.  n_workers>1 uses a multiprocessing Pool.

    Returns the list of records in the same order as `specs` (deterministic).
    Default single-worker: the Pool is imported/created only when n_workers > 1
    (guarded -- no process fan-out on a small host)."""
    if n_workers <= 1:
        return [analyze(s) for s in specs]
    from multiprocessing import Pool
    with Pool(processes=n_workers) as pool:
        return pool.map(analyze, specs)
