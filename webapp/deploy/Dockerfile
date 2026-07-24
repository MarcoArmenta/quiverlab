# quiverlab-web production image. The FastAPI app and the resource-capped worker
# share ONE image; docker-compose runs them as two services (different CMDs) so a
# rebuild keeps the library, webapp, and vendored assets in lockstep.
#
# Base: python:3.12-slim satisfies pyproject `requires-python = ">=3.10"` and
# matches the repo's supported range (3.10-3.13). The default build fetches
# nothing over the network except PyPI wheels: KaTeX is already vendored under
# webapp/static/katex (self-hosted, no CDN at build OR run time), and a TeX
# engine is intentionally omitted -- see the optional block below.
FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# LICENSE is REQUIRED at build time: pyproject sets `license-files = ["LICENSE"]`
# and setuptools>=77 aborts the build if it is missing. README.md backs the
# `readme = "README.md"` metadata. Copy metadata + sources, then install the
# [web] (FastAPI/uvicorn/jinja/pydantic/ulid) and [fast] (numba) extras.
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY webapp ./webapp
RUN pip install -e ".[web,fast]"

# ---------------------------------------------------------------------------
# OPTIONAL -- worked-steps PDF artifacts (a TeX engine).
#
# Without a TeX engine on PATH the runner falls back HONESTLY to a self-contained
# no-JS HTML transcript (see webapp/server/runner.py::_PDF_HTML_FALLBACK); the
# result is correct, just HTML instead of PDF. Real trace.pdf needs pdflatex or
# tectonic. Tectonic is the lighter choice but still adds ~50 MB and a build-time
# download; the default image ships WITHOUT it (see PROVISIONING.md ->
# "Worked-steps PDF"). To enable PDF, uncomment this block and rebuild:
#
# ARG TECTONIC_VERSION=0.15.0
# ARG TECTONIC_ARCH=x86_64-unknown-linux-gnu
# RUN apt-get update \
#     && apt-get install -y --no-install-recommends curl ca-certificates \
#     && curl -fsSL "https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%40${TECTONIC_VERSION}/tectonic-${TECTONIC_VERSION}-${TECTONIC_ARCH}.tar.gz" \
#        | tar xz -C /usr/local/bin tectonic \
#     && tectonic --version \
#     && apt-get purge -y curl && apt-get autoremove -y \
#     && rm -rf /var/lib/apt/lists/*
# ---------------------------------------------------------------------------

# Thread caps for the app tier (the worker child re-pins these to 1 per job).
# QLWEB_DATA_DIR holds the SQLite DB + artifacts; compose mounts a volume here.
ENV NUMBA_NUM_THREADS=2 \
    OMP_NUM_THREADS=2 \
    QLWEB_DATA_DIR=/data

# Run as a non-root user (defense in depth: the worker executes user-submitted
# compute). Pre-create + chown the /data mountpoint so a FRESH named volume
# inherits writable ownership (Docker seeds an empty named volume from the image
# dir's owner/mode). A PRE-EXISTING root-owned volume needs a one-time host-side
# `chown 10001:10001` -- see PROVISIONING.md "Operational limits".
RUN useradd --create-home --uid 10001 quiverlab \
    && mkdir -p /data \
    && chown -R quiverlab:quiverlab /data /app
USER quiverlab

EXPOSE 8000

# App tier: uvicorn calls the create_app() factory. The worker service overrides
# this via entrypoint-worker.sh (python -m webapp.worker.run_loop).
CMD ["uvicorn", "webapp.server.app:create_app", "--factory", "--host", "0.0.0.0", "--port", "8000"]
