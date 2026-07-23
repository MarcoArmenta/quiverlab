# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""Checkpointed, resumable, incremental-N driver for the minimal A^e-resolution.

Pushes the resolution one degree at a time (reusing resolutions_minimal's stepper),
finalizing HH_{n-1} after each advance, recording the per-degree memory (predicted radK
transient + actual VmHWM), and writing an ATOMIC checkpoint after every degree.  A re-run
with the same ckpt_dir resumes from the last completed degree -- so a timeout / node
failure costs at most one degree's recompute, and the memory guard turns the OOM wall into
a clean `stop_reason='memory'` record (the last degree + its radK_bytes = the wall).

Both stepper modes are supported (Plan 15): local algebras run the kernel-accelerated
free path; multi-vertex algebras run the corner-typed projective path (Plan 13).  Corner
data (_CornerContext, gens0, rad_ab_pairs, the engine) is deterministic from (A, prime)
and rebuilt on resume -- corner checkpoints persist only the extra per-degree `tags`.
A ckpt_dir belongs to one (algebra, prime) run: resuming across modes refuses loudly."""

import os
import time
import pickle

from quiverlab.engine.resolutions_minimal import (
    _init_resolution, _advance_resolution, _contracted_degree,
    _corner_contracted_degree)
from quiverlab.engine.hh_engine import rank_mod_p
from quiverlab.errors import QuiverlabError

PRIME = 32003

# Predict the next degree's walltime as this multiple of the last completed degree's
# (degrees grow geometrically); stop BEFORE a degree that won't fit the time budget
# rather than starting it and getting SIGKILL'd mid-degree (a hard TIMEOUT = stale JSON).
_TIME_PRED_FACTOR = 2


def _vmhwm_bytes():
    """Peak resident set size in bytes (Linux /proc), or 0 if unavailable."""
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) * 1024     # kB -> bytes
    except OSError:
        pass
    return 0


def _save_ckpt(ckpt_dir, payload):
    """Atomically write payload (temp file + os.replace); keep the last 2 checkpoints.

    Hardened so an accidental second job sharing this ckpt_dir cannot corrupt it: the temp
    file is per-process (two writers never race on a shared `.tmp` in os.replace -- the
    FileNotFoundError that crashed a doubly-launched job), and `latest.txt` only ever
    advances (a slower/stale writer cannot regress the pointer and lose deeper progress)."""
    os.makedirs(ckpt_dir, exist_ok=True)
    n = payload["n"]
    final = os.path.join(ckpt_dir, "ckpt_%d.pkl" % n)
    tmp = final + ".tmp.%d" % os.getpid()                    # per-process: no shared-tmp race
    with open(tmp, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, final)                                   # atomic
    latest_path = os.path.join(ckpt_dir, "latest.txt")
    cur_latest = -1
    if os.path.exists(latest_path):
        try:
            cur_latest = int(open(latest_path).read().strip())
        except (ValueError, OSError):
            cur_latest = -1
    if n >= cur_latest:                                      # monotonic: never regress
        with open(latest_path, "w") as f:
            f.write(str(n))
    for old in sorted(int(fn[5:-4]) for fn in os.listdir(ckpt_dir)
                      if fn.startswith("ckpt_") and fn.endswith(".pkl")):
        if old < n - 1:
            try:
                os.remove(os.path.join(ckpt_dir, "ckpt_%d.pkl" % old))
            except OSError:
                pass


def _load_ckpt(ckpt_dir):
    """Return the latest checkpoint payload, or None if none / unreadable."""
    latest = os.path.join(ckpt_dir, "latest.txt")
    if not os.path.exists(latest):
        return None
    with open(latest) as f:
        n = int(f.read().strip())
    for cand in (n, n - 1):                                  # tolerate a torn newest
        path = os.path.join(ckpt_dir, "ckpt_%d.pkl" % cand)
        if os.path.exists(path):
            try:
                with open(path, "rb") as f:
                    return pickle.load(f)
            except (pickle.UnpicklingError, EOFError):
                continue
    return None


def deepen(A, ckpt_dir, prime=PRIME, max_transient_bytes=None, max_term_dim=10 ** 9,
           max_degree=None, time_limit_s=None, finalize_only=False, log=None):
    """Build/continue the minimal resolution of A over F_prime, checkpointing per degree.

    Stops at the memory wall (guard), resolution termination, the term cap, `max_degree`,
    or `time_limit_s` (clean checkpoint + resumable).  Returns the summary dict described
    in the plan's interface block.  Works for local AND multi-vertex algebras (the
    corner-typed Plan-13 path); a corner checkpoint additionally persists the per-degree
    corner `tags` -- everything else corner-related is rebuilt from (A, prime) on resume.

    `finalize_only=True` does not advance the resolution: it just re-emits the summary from
    the latest checkpoint (`stop_reason='checkpoint'`) -- the recovery path for a job that
    timed out (hard SIGKILL) before it could write a fresh JSON, leaving the checkpoint
    deeper than the stale on-disk summary."""
    t0 = time.time()
    log = log or (lambda *a: None)
    ck = _load_ckpt(ckpt_dir)
    if finalize_only:
        if ck is None:
            log("deepen: finalize_only but no checkpoint -> empty summary")
            return _summary(A, prime, [], [], "checkpoint", False)
        log("deepen: finalize_only -> summary from checkpoint n=%d (HH len %d)"
            % (ck["n"], len(ck["HH"])))
        return _summary(A, prime, list(ck["HH"]), list(ck["per_degree"]), "checkpoint", False)
    st = _init_resolution(A, prime)         # engine, rad pairs (+ corner ctx, gens0):
    corner = st.get("corner") is not None   # deterministic from (A, prime) -- rebuilt, never pickled
    if ck is None:
        last_gens = None                                    # cols[0] is None
        HH = []
        per_degree = []
        log("deepen: fresh start (dim=%d, prime=%d%s)"
            % (A.m, prime, ", corner mode" if corner else ""))
    else:
        if corner != ("tags" in ck):
            raise QuiverlabError(
                "deepen: ckpt_dir holds a %s-mode checkpoint but this algebra needs "
                "%s mode" % ("corner" if "tags" in ck else "local",
                             "corner" if corner else "local"),
                hint="each ckpt_dir belongs to one (algebra, prime) run -- point "
                     "this run at its own directory")
        st.update({"cur": ck["cur"], "cur_r": ck["cur_r"], "rks": ck["rks"],
                   "n": ck["n"], "cols": {}})
        if corner:
            st["tags"] = ck["tags"]
        last_gens = ck["last_gens"]
        HH = list(ck["HH"])
        per_degree = list(ck["per_degree"])
        log("deepen: resumed at degree n=%d (HH so far: %s)" % (ck["n"], HH))

    m = A.m
    # seed the per-degree wall-time predictor from the last checkpointed degree (old
    # checkpoints have no 'secs' -> 0.0 -> the first resumed degree is attempted unguarded)
    last_secs = (per_degree[-1].get("secs", 0) if per_degree else 0)
    while True:
        if max_degree is not None and st["n"] >= max_degree + 1:
            return _summary(A, prime, HH, per_degree, "max_degree", False)
        if time_limit_s is not None and last_secs > 0:
            elapsed = time.time() - t0
            if elapsed + _TIME_PRED_FACTOR * last_secs > time_limit_s:
                log("deepen: PREDICTIVE TIME STOP before degree %d (elapsed %.0fs + est "
                    "%.0fs > limit %.0fs) -> last degree already checkpointed (resumable)"
                    % (st["n"] + 1, elapsed, _TIME_PRED_FACTOR * last_secs, time_limit_s))
                return _summary(A, prime, HH, per_degree, "time", False)
        n_prev = st["n"]
        d_t0 = time.time()
        res = _advance_resolution(st, prime, max_term_dim, max_transient_bytes)
        k = st["n"]                                         # degree just attempted/built
        status = res["status"]
        if status == "memory":
            log("deepen: MEMORY WALL at degree %d (radK ~%.1f GB > budget)"
                % (n_prev + 1, (res["radK_bytes"] or 0) / 10 ** 9))
            return _summary(A, prime, HH, per_degree, "memory", False,
                            wall_degree=n_prev + 1, wall_radK_bytes=res["radK_bytes"])
        # we built (or terminated at) degree k; gens = cols[k]
        gens = st["cols"].get(k, [])
        # finalize HH_{k-1} now that d_k is known (rolling pair last_gens=cols[k-1], gens)
        if k >= 1:
            r_km1 = st["rks"].get(k - 1, 0)
            if corner:
                # contracted blocks are corners e_w A e_v, not full copies of A
                ctx = st["corner"]
                tags = st["tags"]
                dbar_km1 = (_corner_contracted_degree(
                                st["eng"], ctx, last_gens or [], tags.get(k - 1, []),
                                tags.get(k - 2, []), prime)
                            if (k - 1) >= 1 and r_km1 > 0 else None)
                dbar_k = (_corner_contracted_degree(
                              st["eng"], ctx, gens or [], tags.get(k, []),
                              tags.get(k - 1, []), prime)
                          if st["rks"].get(k, 0) > 0 else None)
                dimn = sum(ctx.corner_dim_A(tg) for tg in tags.get(k - 1, []))
            else:
                r_km2 = st["rks"].get(k - 2, 0)
                dbar_km1 = (_contracted_degree(st["eng"], last_gens or [], r_km2, k - 1)
                            if (k - 1) >= 1 and r_km1 > 0 else None)
                dbar_k = (_contracted_degree(st["eng"], gens or [], r_km1, k)
                          if st["rks"].get(k, 0) > 0 else None)
                dimn = m * r_km1
            rn = rank_mod_p(dbar_km1, prime) if dbar_km1 is not None else 0
            rnp1 = rank_mod_p(dbar_k, prime) if dbar_k is not None else 0
            HH.append(int(dimn - rn - rnp1))                # HH_{k-1}
        deg_secs = time.time() - d_t0
        per_degree.append({"n": k, "r_n": res["r_n"], "radK_bytes": res["radK_bytes"],
                           "vmhwm_bytes": _vmhwm_bytes(), "secs": deg_secs})
        last_secs = deg_secs
        log("deepen: degree %d done  r_%d=%s  HH_0..=%s  VmHWM=%.1fGB"
            % (k, k, res["r_n"], HH, _vmhwm_bytes() / 10 ** 9))

        if status == "terminated":
            return _summary(A, prime, HH, per_degree, "terminated", True,
                            hochschild_dim=k - 1)
        if status == "term":
            return _summary(A, prime, HH, per_degree, "term", False)

        last_gens = gens                                    # roll forward
        payload = {"n": st["n"], "cur": st["cur"], "cur_r": st["cur_r"],
                   "rks": st["rks"], "last_gens": last_gens,
                   "HH": HH, "per_degree": per_degree}
        if corner:
            payload["tags"] = st["tags"]        # tiny: one vertex pair per generator
        _save_ckpt(ckpt_dir, payload)
        st["cols"] = {}                                     # free history (rolling)
        if time_limit_s is not None and (time.time() - t0) > time_limit_s:
            log("deepen: TIME LIMIT after degree %d -> checkpoint + exit (resumable)" % k)
            return _summary(A, prime, HH, per_degree, "time", False)


def _summary(A, prime, HH, per_degree, stop_reason, terminated, **extra):
    out = {"dim": int(A.m), "prime": int(prime),
           "max_degree_reached": (len(HH) - 1) if HH else -1,
           "HH": list(HH), "stop_reason": stop_reason, "terminated": bool(terminated),
           "per_degree": per_degree}
    out.update(extra)
    return out
