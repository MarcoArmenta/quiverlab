"""Exact-integer host-resource detection for the whole HPC tier (Plan 28).

``detect_resources()`` reports the cores / memory / GPUs the current process may
actually use -- the MINIMUM over every visible limit (OS count, CPU affinity,
cgroup v1/v2 quotas, SLURM allocation env), so a job never over-subscribes the
slice the scheduler handed it. Stdlib only; all arithmetic is integer (the
``src/`` float gate forbids float literals and ``float()`` calls). No caching, but
at most ONE ``nvidia-smi`` attempt per call (guarded) -- the offline-GUI sibling
imports ``detect_resources`` from webapp code, so the signature stays stable/cheap.

GPUs are REPORT-ONLY: quiverlab's exact engines are CPU-only, so the dict carries
``{"gpus": n, "gpus_used": False}`` and any user-facing print says "detected but
not used (exact CPU engines)".
"""
from __future__ import annotations

import os
import shutil
import subprocess

# cgroup / proc paths (module constants so tests can monkeypatch them to tmp files).
CGROUP_V2_CPU_MAX = "/sys/fs/cgroup/cpu.max"
CGROUP_V1_CFS_QUOTA = "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
CGROUP_V1_CFS_PERIOD = "/sys/fs/cgroup/cpu/cpu.cfs_period_us"
CGROUP_V2_MEM_MAX = "/sys/fs/cgroup/memory.max"
CGROUP_V1_MEM_LIMIT = "/sys/fs/cgroup/memory/memory.limit_in_bytes"
MEMINFO = "/proc/meminfo"

# A cgroup "unlimited" sentinel: memory.limit_in_bytes reports a near-2^63 value
# when uncapped. No real machine has this much RAM, so treat it as no limit.
_MEM_UNLIMITED = 1 << 62
_MiB = 1024 * 1024
_KiB = 1024
# Default cap on the auto-detected thread count when NOT under a SLURM allocation
# (a laptop with 32 logical cores should not spin up 32 numba threads by default).
DEFAULT_CORE_CAP = 8


def _read_text(path: str) -> str | None:
    try:
        with open(path, "r", encoding="utf-8") as f:
            return f.read()
    except OSError:
        return None


def _read_int_file(path: str) -> int | None:
    """A file holding a single integer (cgroup limit files), or None."""
    txt = _read_text(path)
    if txt is None:
        return None
    txt = txt.strip()
    try:
        return int(txt)
    except ValueError:
        return None


def _ceil_div(a: int, b: int) -> int:
    return -(-a // b)


# --------------------------------------------------------------------------- #
# Cores
# --------------------------------------------------------------------------- #

def cgroup_v2_cpu(path: str | None = None) -> int | None:
    """cgroup v2 ``cpu.max`` = "<quota> <period>" (or "max <period>" = uncapped)
    -> ceil(quota/period) cores, or None. ``path`` defaults to the module constant
    at CALL time (so tests can monkeypatch it)."""
    txt = _read_text(CGROUP_V2_CPU_MAX if path is None else path)
    if txt is None:
        return None
    parts = txt.split()
    if len(parts) < 2 or parts[0] == "max":
        return None
    try:
        quota, period = int(parts[0]), int(parts[1])
    except ValueError:
        return None
    if quota <= 0 or period <= 0:
        return None
    return max(1, _ceil_div(quota, period))


def cgroup_v1_cpu(quota_path: str | None = None,
                  period_path: str | None = None) -> int | None:
    """cgroup v1 CFS ``cpu.cfs_quota_us`` / ``cpu.cfs_period_us`` -> cores, or
    None when unset (quota == -1) or unreadable."""
    quota = _read_int_file(CGROUP_V1_CFS_QUOTA if quota_path is None else quota_path)
    period = _read_int_file(CGROUP_V1_CFS_PERIOD if period_path is None else period_path)
    if quota is None or period is None or quota <= 0 or period <= 0:
        return None
    return max(1, _ceil_div(quota, period))


def slurm_cpus(env=None) -> int | None:
    env = os.environ if env is None else env
    raw = env.get("SLURM_CPUS_PER_TASK")
    if raw is None or raw == "":
        return None
    try:
        n = int(raw)
    except ValueError:
        return None
    return n if n > 0 else None


def _affinity_cpus() -> int | None:
    getaff = getattr(os, "sched_getaffinity", None)
    if getaff is None:
        return None
    try:
        return len(getaff(0))
    except OSError:
        return None


def detect_cores(env=None) -> int:
    """The usable core count: the MIN over os.cpu_count, CPU affinity, cgroup v2,
    cgroup v1, and ``$SLURM_CPUS_PER_TASK`` -- always at least 1."""
    env = os.environ if env is None else env
    candidates = []
    cpu_count = os.cpu_count()
    if cpu_count:
        candidates.append(cpu_count)
    for c in (_affinity_cpus(), cgroup_v2_cpu(), cgroup_v1_cpu(), slurm_cpus(env)):
        if c is not None:
            candidates.append(c)
    return max(1, min(candidates)) if candidates else 1


# --------------------------------------------------------------------------- #
# Memory
# --------------------------------------------------------------------------- #

def meminfo_total(path: str | None = None) -> int | None:
    """Linux ``/proc/meminfo`` MemTotal in bytes (the line is in kB), or None."""
    txt = _read_text(MEMINFO if path is None else path)
    if txt is None:
        return None
    for line in txt.splitlines():
        if line.startswith("MemTotal:"):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    return int(parts[1]) * _KiB
                except ValueError:
                    return None
    return None


def sysctl_memsize() -> int | None:
    """macOS physical RAM via ``sysctl -n hw.memsize`` (bytes), or None."""
    if not shutil.which("sysctl"):
        return None
    try:
        out = subprocess.run(["sysctl", "-n", "hw.memsize"],
                             capture_output=True, text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    try:
        return int(out.stdout.strip())
    except ValueError:
        return None


def cgroup_v2_mem(path: str | None = None) -> int | None:
    """cgroup v2 ``memory.max`` in bytes ("max" = unlimited -> None)."""
    txt = _read_text(CGROUP_V2_MEM_MAX if path is None else path)
    if txt is None:
        return None
    txt = txt.strip()
    if txt == "max":
        return None
    try:
        v = int(txt)
    except ValueError:
        return None
    return v if 0 < v < _MEM_UNLIMITED else None


def cgroup_v1_mem(path: str | None = None) -> int | None:
    """cgroup v1 ``memory.limit_in_bytes`` (near-2^63 sentinel = unlimited)."""
    v = _read_int_file(CGROUP_V1_MEM_LIMIT if path is None else path)
    if v is None:
        return None
    return v if 0 < v < _MEM_UNLIMITED else None


def slurm_mem_bytes(env=None, cores: int | None = None) -> int | None:
    """SLURM memory allocation in bytes: ``$SLURM_MEM_PER_NODE`` (MB) or, failing
    that, ``$SLURM_MEM_PER_CPU`` (MB) x cores."""
    env = os.environ if env is None else env
    per_node = env.get("SLURM_MEM_PER_NODE")
    if per_node:
        try:
            return int(per_node) * _MiB
        except ValueError:
            pass
    per_cpu = env.get("SLURM_MEM_PER_CPU")
    if per_cpu:
        try:
            mb = int(per_cpu)
        except ValueError:
            return None
        n = cores if cores is not None else detect_cores(env)
        return mb * _MiB * n
    return None


def detect_mem_bytes(env=None, cores: int | None = None) -> int | None:
    """Usable memory in bytes: the MIN over /proc/meminfo (or macOS sysctl),
    cgroup v2, cgroup v1, and the SLURM allocation env. None if nothing is
    detectable (unknown host)."""
    env = os.environ if env is None else env
    candidates = []
    host = meminfo_total()
    if host is None:
        host = sysctl_memsize()
    for c in (host, cgroup_v2_mem(), cgroup_v1_mem(), slurm_mem_bytes(env, cores)):
        if c is not None:
            candidates.append(c)
    return min(candidates) if candidates else None


# --------------------------------------------------------------------------- #
# GPUs (report-only; the exact engines are CPU-only)
# --------------------------------------------------------------------------- #

def detect_gpus() -> tuple[int, list]:
    """(count, names) from a single guarded ``nvidia-smi -L`` -- (0, []) when the
    binary is absent or fails. NEVER raises: absence is normal."""
    if not shutil.which("nvidia-smi"):
        return 0, []
    try:
        out = subprocess.run(["nvidia-smi", "-L"],
                             capture_output=True, text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return 0, []
    if out.returncode != 0:
        return 0, []
    names = []
    for line in out.stdout.splitlines():
        line = line.strip()
        if not line.startswith("GPU"):
            continue
        # "GPU 0: NVIDIA A100 (UUID: GPU-...)" -> "NVIDIA A100"
        body = line.split(":", 1)[1] if ":" in line else line
        name = body.split("(UUID")[0].strip()
        names.append(name or line)
    return len(names), names


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def detect_resources(env=None) -> dict:
    """Host resources as exact ints:
    ``{"cores", "mem_bytes"(|None), "gpus", "gpu_names", "gpus_used": False}``.
    Cheap and side-effect-free apart from at most one guarded ``nvidia-smi``."""
    env = os.environ if env is None else env
    cores = detect_cores(env)
    mem = detect_mem_bytes(env, cores)
    n_gpu, gpu_names = detect_gpus()
    return {
        "cores": cores,
        "mem_bytes": mem,
        "gpus": n_gpu,
        "gpu_names": gpu_names,
        "gpus_used": False,     # exact CPU engines: GPUs are detected, never used
    }


def default_thread_cap(env=None) -> int:
    """Default thread count for the CLI: ``$SLURM_CPUS_PER_TASK`` verbatim when
    under a SLURM allocation, else detected cores capped at ``DEFAULT_CORE_CAP``
    (a laptop should not default to all logical cores)."""
    env = os.environ if env is None else env
    slurm = slurm_cpus(env)
    if slurm is not None:
        return slurm
    return min(detect_cores(env), DEFAULT_CORE_CAP)


def format_resources(res: dict) -> str:
    """A one-line human summary for ``version`` / ``estimate`` / ``selftest``."""
    mem = res.get("mem_bytes")
    mem_str = ("%d MiB" % (mem // _MiB)) if mem else "unknown"
    line = "host: %d core(s), %s RAM" % (res["cores"], mem_str)
    if res.get("gpus"):
        line += "; %d GPU(s) detected but not used (exact CPU engines)" % res["gpus"]
    return line
