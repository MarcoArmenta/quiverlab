# quiverlab HPC + offline-GUI image (Plan 28).
#
# ONE image, two user stories:
#   * HPC batch      -- `quiverlab-hpc run config.yaml -o result.json` on a cluster
#                       (Apptainer/Singularity preferred; Docker works too), then
#                       `quiverlab-hpc render result.json -o report.html` locally.
#   * Offline laptop -- `quiverlab-hpc gui` serves the webapp on localhost with an
#                       embedded worker, fully offline (vendored KaTeX,
#                       seeded example cache).
#
# The CLI itself lives in the WHEEL (the `quiverlab-hpc` console script from the
# [hpc] extra); this image is a convenience wrapper. The webapp deploy assets in
# webapp/deploy/ are a SEPARATE artifact (the public server tier) -- do not confuse
# the two.
#
# KEEP THE BUILD STEPS IN SYNC with container/quiverlab.def -- the Apptainer
# definition's %post mirrors every RUN below (same base digest, same
# `pip install -e ".[fast,hpc,web]"`, same guarded seed).
#
# ---------------------------------------------------------------------------
# Base: python:3.12-slim, PINNED BY DIGEST for linux/amd64 (HPC clusters + CI are
# x86-64). Refresh the digest with:
#   docker manifest inspect python:3.12-slim \
#     | jq -r '.manifests[] | select(.platform.architecture=="amd64" and .platform.os=="linux") | .digest'
#
# CPU-ONLY BY DESIGN -- a plain python:slim base, NOT a CUDA/GPU image. quiverlab
# computes with EXACT arithmetic (integers mod p, rationals); there is no floating
# point and no GPU kernel anywhere, so a GPU would sit idle. Do not swap in an
# nvidia/cuda base or request GPUs on the cluster -- ask for cores and RAM instead
# (see slurm/*.sbatch and docs/hpc.md "What resources to request").
# ---------------------------------------------------------------------------
FROM python:3.12-slim@sha256:cab2dbf575e971934a81e4622f5aba17aa7929719bd7e31033a3a83b97fd0464

# Build-time version stamp: CI (container.yml) passes
# --build-arg QUIVERLAB_VERSION=<git tag without the leading v> and asserts the
# running image reports it (see the "Image version" step in the workflow).
ARG QUIVERLAB_VERSION=0.0.0-dev

LABEL org.opencontainers.image.title="quiverlab" \
      org.opencontainers.image.description="Exact Hochschild theory for quivers with relations -- HPC batch + offline GUI container" \
      org.opencontainers.image.source="https://github.com/MarcoArmenta/quiverlab" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="${QUIVERLAB_VERSION}"

# Fully offline, reproducible, headless runtime.
#   QUIVERLAB_SEED_CACHE-> build->runtime handoff path for the seeded Plan-25 cache
#                          DB; the `gui` verb copies it into the data dir on first run.
#   PYTHONPATH=/app     -> `import webapp` for the `gui` verb: webapp/ is not
#                          packaged into the wheel, and a console script does not
#                          put the WORKDIR on sys.path, so without this the gui
#                          verb dies with "No module named 'webapp'".
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    MPLBACKEND=Agg \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    NUMBA_NUM_THREADS=2 \
    OMP_NUM_THREADS=2 \
    PYTHONPATH=/app \
    QUIVERLAB_IMAGE_VERSION=${QUIVERLAB_VERSION} \
    QUIVERLAB_SEED_CACHE=/opt/quiverlab/seed-cache.db

WORKDIR /app

# --- sources + install -----------------------------------------------------
# LICENSE + README back the pyproject metadata (setuptools>=77 aborts without
# LICENSE). webapp/ is the offline GUI surface and is NOT packaged into the wheel,
# so the install is EDITABLE (-e) exactly like webapp/deploy/Dockerfile -- that
# keeps `import webapp` importable for the `gui` verb from the in-image source tree.
# container/ carries the guarded seed hook below. [fast,hpc,web]: numba + the HPC
# CLI + the offline webapp; NO [qpa] (test-only, huge).
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY webapp ./webapp
COPY container ./container
RUN pip install -e ".[fast,hpc,web]"

# --- offline example cache seed (Plan 25 result_cache) ---------------------
# OPTIONAL + GUARDED. The seeding script (container/seed_cache.py) and the manifest
# (webapp/precomputed/manifest.yaml) are delivered by the wheel/webapp side of Plan
# 28. This step is written so the image builds GREEN whether or not they have landed
# yet: the guard skips cleanly when either is absent, and even a present-but-failing
# seed never breaks the build (best-effort `|| true`).
RUN mkdir -p /opt/quiverlab \
    && if [ -f container/seed_cache.py ] && [ -f webapp/precomputed/manifest.yaml ]; then \
         python container/seed_cache.py --manifest webapp/precomputed/manifest.yaml \
             --out "$QUIVERLAB_SEED_CACHE" || true ; \
       else \
         echo "[quiverlab] seed_cache.py/manifest not present yet -- skipping example-cache seed (Plan 28 wheel side)" ; \
       fi

# --- non-root user (defense in depth; mirrors webapp/deploy/Dockerfile) ----
# Apptainer/Singularity run as the calling user regardless; under Docker we drop to
# uid 10001. The user owns /opt/quiverlab so the seed DB is readable.
RUN useradd --create-home --uid 10001 quiverlab \
    && chown -R quiverlab:quiverlab /app /opt/quiverlab
USER quiverlab

# HPC batch is the default verb surface; the offline GUI is `... gui` (port 8000).
ENTRYPOINT ["quiverlab-hpc"]
CMD ["--help"]
