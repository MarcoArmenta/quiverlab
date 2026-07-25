"""Exact-int host detection (Plan 28 addendum): cgroup v1/v2 quotas, SLURM env,
absent nvidia-smi, and the MIN-over-limits precedence. Stdlib only, deterministic
(monkeypatched files/env), no real host probing in the unit assertions."""
from quiverlab.hpc import resources


# --------------------------------------------------------------------------- #
# Cores
# --------------------------------------------------------------------------- #

def test_cgroup_v2_cpu_quota(tmp_path):
    f = tmp_path / "cpu.max"
    f.write_text("200000 100000\n")
    assert resources.cgroup_v2_cpu(str(f)) == 2
    f.write_text("250000 100000\n")           # ceil(2.5) = 3
    assert resources.cgroup_v2_cpu(str(f)) == 3
    f.write_text("max 100000\n")               # uncapped
    assert resources.cgroup_v2_cpu(str(f)) is None


def test_cgroup_v1_cpu_quota(tmp_path):
    q, p = tmp_path / "quota", tmp_path / "period"
    q.write_text("400000")
    p.write_text("100000")
    assert resources.cgroup_v1_cpu(str(q), str(p)) == 4
    q.write_text("-1")                          # unset
    assert resources.cgroup_v1_cpu(str(q), str(p)) is None


def test_slurm_cpus_env():
    assert resources.slurm_cpus({"SLURM_CPUS_PER_TASK": "3"}) == 3
    assert resources.slurm_cpus({}) is None
    assert resources.slurm_cpus({"SLURM_CPUS_PER_TASK": "0"}) is None


def test_detect_cores_takes_the_minimum(monkeypatch, tmp_path):
    # A tiny cgroup v2 quota must win the MIN over the (larger) real cpu_count.
    f = tmp_path / "cpu.max"
    f.write_text("100000 100000\n")            # exactly 1 core
    monkeypatch.setattr(resources, "CGROUP_V2_CPU_MAX", str(f))
    assert resources.detect_cores({}) == 1
    # SLURM even smaller than everything else also wins.
    assert resources.detect_cores({"SLURM_CPUS_PER_TASK": "1"}) == 1


def test_default_thread_cap_precedence(monkeypatch):
    # SLURM_CPUS_PER_TASK wins verbatim, even above the laptop cap.
    assert resources.default_thread_cap({"SLURM_CPUS_PER_TASK": "16"}) == 16
    # No SLURM: detected cores capped at DEFAULT_CORE_CAP.
    cap = resources.default_thread_cap({})
    assert 1 <= cap <= resources.DEFAULT_CORE_CAP


# --------------------------------------------------------------------------- #
# Memory
# --------------------------------------------------------------------------- #

def test_meminfo_total(tmp_path):
    f = tmp_path / "meminfo"
    f.write_text("MemTotal:       16384 kB\nMemFree:  100 kB\n")
    assert resources.meminfo_total(str(f)) == 16384 * 1024
    assert resources.meminfo_total(str(tmp_path / "absent")) is None


def test_cgroup_v2_mem(tmp_path):
    f = tmp_path / "memory.max"
    f.write_text("2147483648\n")
    assert resources.cgroup_v2_mem(str(f)) == 2147483648
    f.write_text("max\n")
    assert resources.cgroup_v2_mem(str(f)) is None


def test_cgroup_v1_mem_sentinel_is_unlimited(tmp_path):
    f = tmp_path / "limit"
    f.write_text("2147483648")
    assert resources.cgroup_v1_mem(str(f)) == 2147483648
    f.write_text(str(1 << 63))                  # near-2^63 sentinel = unlimited
    assert resources.cgroup_v1_mem(str(f)) is None


def test_slurm_mem_bytes():
    assert resources.slurm_mem_bytes({"SLURM_MEM_PER_NODE": "1024"}) == 1024 * 1024 * 1024
    assert resources.slurm_mem_bytes({"SLURM_MEM_PER_CPU": "512"}, cores=2) == \
        512 * 1024 * 1024 * 2
    assert resources.slurm_mem_bytes({}) is None


def test_detect_mem_bytes_takes_the_minimum(monkeypatch, tmp_path):
    # A small SLURM allocation caps the reported memory below the host RAM.
    small = 256 * 1024 * 1024
    got = resources.detect_mem_bytes({"SLURM_MEM_PER_NODE": "256"})
    assert got is not None and got <= small


# --------------------------------------------------------------------------- #
# GPUs (report-only)
# --------------------------------------------------------------------------- #

def test_detect_gpus_absent_binary(monkeypatch):
    monkeypatch.setattr(resources.shutil, "which", lambda name: None)
    assert resources.detect_gpus() == (0, [])


def test_detect_resources_shape_and_gpu_flag():
    r = resources.detect_resources({})
    assert set(r) >= {"cores", "mem_bytes", "gpus", "gpu_names", "gpus_used"}
    assert r["cores"] >= 1
    assert r["gpus_used"] is False               # exact CPU engines never use GPUs
    assert isinstance(r["gpus"], int)


def test_format_resources_mentions_gpu_not_used():
    line = resources.format_resources({"cores": 4, "mem_bytes": 8 * 1024 ** 3,
                                       "gpus": 2, "gpu_names": ["A", "B"],
                                       "gpus_used": False})
    assert "4 core" in line and "not used" in line
    line2 = resources.format_resources({"cores": 1, "mem_bytes": None, "gpus": 0,
                                        "gpu_names": [], "gpus_used": False})
    assert "unknown" in line2
