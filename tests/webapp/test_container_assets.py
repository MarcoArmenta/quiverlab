"""Structural assertions on the Plan-28 container / SLURM / CI assets.

Mirrors tests/webapp/test_deploy_assets.py: text-based checks on committed files so
they run in the fast bucket with no Docker/Apptainer/SLURM. The image build itself
is exercised by .github/workflows/container.yml and container/SMOKE.md.
"""
import pathlib

import yaml

from webapp.server.schema import ComputeRequest

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
CONTAINER = ROOT / "container"
SLURM = ROOT / "slurm"
WF = ROOT / ".github" / "workflows"
SBATCH_FILES = [SLURM / "quiverlab.sbatch", SLURM / "quiverlab-drac.sbatch"]


def _read(p: pathlib.Path) -> str:
    return p.read_text(encoding="utf-8")


# --------------------------------------------------------------------------- #
# Dockerfile
# --------------------------------------------------------------------------- #
def _dockerfile_from() -> str:
    for line in _read(CONTAINER / "Dockerfile").splitlines():
        s = line.strip()
        if s.startswith("FROM "):
            return s[len("FROM "):].strip()
    raise AssertionError("container/Dockerfile has no FROM line")


def _dockerfile_pip_line() -> str:
    for line in _read(CONTAINER / "Dockerfile").splitlines():
        if "pip install" in line:
            return line
    raise AssertionError("container/Dockerfile has no pip install line")


def test_dockerfile_core_directives():
    text = _read(CONTAINER / "Dockerfile")
    # HPC batch is the default surface.
    assert 'ENTRYPOINT ["quiverlab-hpc"]' in text
    # Offline, reproducible, headless runtime env.
    assert "LANG=C.UTF-8" in text
    assert "MPLBACKEND=Agg" in text
    # The exact extras on the install line: numba + HPC CLI + offline webapp, no qpa
    # (check the pip line itself so an explanatory comment cannot spoof it).
    pip = _dockerfile_pip_line()
    assert "[fast,hpc,web]" in pip
    assert "qpa" not in pip
    # Base pinned by digest (linux/amd64).
    assert _dockerfile_from().startswith("python:3.12-slim@sha256:")


def test_dockerfile_has_no_pdf_toolchain():
    # PDF/TeX report output was removed (reports are HTML + JSON; text for console).
    # The image must no longer carry a TeX engine or a PDF text-extraction tool -- the
    # Apptainer def mirrors the Dockerfile, so pin both.
    for name in ("Dockerfile", "quiverlab.def"):
        low = _read(CONTAINER / name).lower()
        assert "tectonic" not in low, name
        assert "pdftotext" not in low, name
        assert "poppler" not in low, name


def test_dockerfile_is_cpu_only_by_design():
    # No CUDA/GPU base (check the FROM ref itself, not comments); explicit CPU-only
    # note (scope addendum).
    text = _read(CONTAINER / "Dockerfile")
    assert "CPU-ONLY" in text
    base = _dockerfile_from().lower()
    assert "cuda" not in base
    assert "nvidia" not in base
    assert base.startswith("python:")


# --------------------------------------------------------------------------- #
# Apptainer definition -- must stay in sync with the Dockerfile
# --------------------------------------------------------------------------- #
def _def_from() -> str:
    for line in _read(CONTAINER / "quiverlab.def").splitlines():
        s = line.strip()
        if s.lower().startswith("from:"):
            return s.split(":", 1)[1].strip()
    raise AssertionError("container/quiverlab.def has no From: line")


def test_def_matches_dockerfile_base_and_runscript():
    text = _read(CONTAINER / "quiverlab.def")
    assert text.strip().splitlines()[0].startswith("#")           # has a header comment
    assert "Bootstrap: docker" in text
    # Same base image (digest-pinned) as the Dockerfile.
    assert _def_from() == _dockerfile_from()
    # Runscript execs the CLI.
    assert "%runscript" in text
    assert "exec quiverlab-hpc" in text
    # Same install extras as the Dockerfile (sync).
    assert "[fast,hpc,web]" in text


# --------------------------------------------------------------------------- #
# sbatch templates
# --------------------------------------------------------------------------- #
def _first_code_line_with(text: str, token: str) -> str | None:
    """First NON-comment, non-blank line mentioning ``token`` (so a `$SCRATCH`
    inside an explanatory comment does not count as executable use)."""
    for line in text.splitlines():
        stripped = line.lstrip()
        if not stripped or stripped.startswith("#"):
            continue
        if token in line:
            return line
    return None


def test_sbatch_signal_and_strict_mode():
    for f in SBATCH_FILES:
        text = _read(f)
        assert text.startswith("#!/bin/bash")
        assert "set -euo pipefail" in text                # strict mode
        assert "#SBATCH --signal" in text                 # checkpoint signal


def test_sbatch_three_way_launcher_fallback():
    for f in SBATCH_FILES:
        text = _read(f)
        assert "QUIVERLAB_LAUNCH" in text                 # env override
        assert "apptainer" in text
        assert "singularity" in text
        assert "python -m quiverlab.hpc" in text          # venv/emulator fallback


def test_sbatch_scratch_defined_before_use():
    # The FIRST executable mention of SCRATCH must be the default-assignment
    # (`${SCRATCH:=...}`), so `set -u` never trips off-cluster.
    for f in SBATCH_FILES:
        first = _first_code_line_with(_read(f), "SCRATCH")
        assert first is not None, f"{f.name}: no SCRATCH usage found"
        assert "${SCRATCH:=" in first, (
            f"{f.name}: first SCRATCH use is not a default-assignment: {first!r}")


def test_sbatch_checkpoint_resume_contract():
    for f in SBATCH_FILES:
        text = _read(f)
        assert "-eq 75" in text                           # exit-75 == clean ckpt stop
        assert "scontrol requeue" in text                 # primary resume
        assert "--checkpoint-dir" in text                 # ckpt dir threaded through


def test_sbatch_is_cpu_only_by_design():
    # Prominent no-GPU guidance (scope addendum).
    for f in SBATCH_FILES:
        text = _read(f)
        assert "CPU-ONLY" in text
        low = text.lower()
        assert "do not request gpus" in low or "do not request a gpu" in low
        assert "estimate" in low                          # points at the estimator


def test_drac_template_header():
    text = _read(SLURM / "quiverlab-drac.sbatch")
    assert "#SBATCH --account=def-CHANGEME" in text
    assert "module load apptainer" in text


# --------------------------------------------------------------------------- #
# CI workflow
# --------------------------------------------------------------------------- #
def test_container_workflow_shape():
    text = _read(WF / "container.yml")
    # Valid YAML (the `on:` -> True 1.1 gotcha is fine; we assert triggers by text).
    assert yaml.safe_load(text) is not None
    # Tag-gated + manual.
    assert 'tags: ["v*"]' in text
    assert "workflow_dispatch" in text
    # Builds + smokes + pushes to GHCR (buildx: amd64 smoked locally, then a
    # multi-arch amd64+arm64 push -- Apple-Silicon laptops run natively).
    assert "docker/build-push-action" in text
    assert "ghcr.io" in text
    assert "push: true" in text
    assert "linux/amd64,linux/arm64" in text
    assert "packages: write" in text
    # Registry paths are lowercase-only; the owner is "MarcoArmenta", so the
    # image name MUST go through the bash lowercase expansion.
    assert "${GITHUB_REPOSITORY_OWNER,,}" in text
    assert "ghcr.io/${{ github.repository_owner }}" not in text
    # The gui verb is smoked THROUGH a published port (the Plan-28 image shipped
    # with gui broken: webapp unimportable + loopback-only bind; CI never ran it).
    assert "gui --no-open" in text
    assert "curl" in text
    # Renders an HTML report inside docker and checks it has real content -- no PDF
    # toolchain (reports are HTML + JSON now).
    assert "render" in text
    assert "report.html" in text
    assert "pdftotext" not in text
    # SIF fallback asset with the < 2 GiB assertion.
    assert "apptainer build" in text
    assert "1024 * 1024 * 1024" in text


# --------------------------------------------------------------------------- #
# CI fixture parity with the request schema
# --------------------------------------------------------------------------- #
def test_ci_tiny_config_validates_as_compute_request():
    data = yaml.safe_load(_read(CONTAINER / "ci-tiny.yaml"))
    # The `hpc:` block is spec-core-only; the HTTP schema ignores it. Drop it and
    # the remainder must validate as a ComputeRequest (the run/render parity pin).
    data.pop("hpc", None)
    req = ComputeRequest.model_validate(data)
    assert req.schema_version == 1
    assert req.compute == ["hh_cohomology:0..2"]
