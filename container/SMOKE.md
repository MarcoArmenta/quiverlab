# Container smoke runbooks

Manual, copy-pasteable smoke tests for the Plan-28 image. None of these are a
pytest (they need Docker / Apptainer / a SLURM emulator); the host pytest suite
covers the CLI and the asset shapes. Run these before cutting a release or after
touching `container/**` or `slurm/**`.

The image ships one entrypoint, `quiverlab-hpc`, with verbs
`run / render / sample-config / estimate / gui / version / selftest`.

---

## (a) Docker: build, selftest, tiny run + render

From the repo root:

```bash
# Build (context = repo root; the .dockerignore keeps .venv/.git/docs/tests out).
docker build -f container/Dockerfile --build-arg QUIVERLAB_VERSION=0.1.0.dev0 \
    -t quiverlab:local .

# Sanity: version + selftest.
docker run --rm quiverlab:local version
docker run --rm quiverlab:local selftest

# Tiny run + render. The container writes to /out; make it world-writable so the
# in-image non-root uid (10001) can write to the bind mount.
mkdir -p out && chmod 777 out
docker run --rm -v "$PWD/container:/cfg:ro" -v "$PWD/out:/out" \
    quiverlab:local run /cfg/ci-tiny.yaml -o /out/result.json
docker run --rm -v "$PWD/out:/out" \
    quiverlab:local render /out/result.json -o /out/report.html --format html

# Assert: result envelope + a non-empty rendered HTML report (no PDF toolchain;
# math is KaTeX/MathML, so no network is needed).
python -c "import json; d=json.load(open('out/result.json')); print('keys:', sorted(d)); assert 'result_schema' in d"
test -s out/report.html && echo "report.html OK ($(wc -c < out/report.html) bytes)"
grep -q "quiverlab report" out/report.html && echo "report.html has rendered content"
```

Offline GUI (the laptop app):

```bash
# The image sets QUIVERLAB_GUI_HOST=0.0.0.0 so the published port reaches the
# GUI (loopback inside the container's namespace would not be; only published
# ports are exposed, so this stays local-only). Works on Docker Desktop too.
docker run --rm -p 8000:8000 quiverlab:local gui
# open http://localhost:8000  -- computes locally, shows memory/time estimates,
# ships the seeded example cache. No internet required after the image is pulled.
```

---

## (b) Apptainer: pull from GHCR, run

Rootless, no admin. The PRIMARY path -- no `.def` build needed:

```bash
export APPTAINER_CACHEDIR="${SCRATCH:-$PWD}/.apptainer"   # keep the cache off $HOME quota
apptainer pull quiverlab.sif docker://ghcr.io/marcoarmenta/quiverlab:latest

apptainer run quiverlab.sif version
apptainer run quiverlab.sif selftest

# tiny run + render (bind $PWD so the SIF sees your files)
apptainer run --bind "$PWD" quiverlab.sif run container/ci-tiny.yaml -o result.json
apptainer run --bind "$PWD" quiverlab.sif render result.json -o report.html --format html
test -s report.html && echo "report.html OK"

# named app entrypoints also work:
apptainer run --app gui quiverlab.sif           # offline GUI on :8000
```

Native SIF build (only if you prefer building over pulling):

```bash
apptainer build quiverlab.sif container/quiverlab.def
```

---

## (c) drac-local emulator: the venv-fallback orchestration

The drac-local emulator has **no Apptainer** in its compute images, so this
exercises the `python -m quiverlab.hpc` launcher fallback baked into the sbatch
templates -- i.e. the SLURM orchestration (headers, USR1 handler, checkpoint dir,
resume), NOT the container itself.

```bash
# 1. Make the project installable on the emulator's Linux venv. Add an editable
#    install line to requirements-cluster.txt at the repo root:
printf -- '-e .[fast,hpc]\n' > requirements-cluster.txt

# 2. Build the Linux venv (one-time per project).
slurm-local prepare "$PWD"

# 3. Activate a cluster profile and submit the DRAC template with a tiny config.
#    No apptainer/singularity on PATH -> the launcher auto-falls back to
#    `python -m quiverlab.hpc`.
eval "$(slurm-local activate narval)"
cp container/ci-tiny.yaml tiny.yaml
sbatch slurm/quiverlab-drac.sbatch tiny.yaml result.json

# 4. Watch it, then assert the artifacts appear and render.
slurm-local status
#   ... wait for the job to finish (squeue empty) ...
test -s result.json && echo "result.json OK"
python -m quiverlab.hpc render result.json -o report.html --format html
test -s report.html && echo "report.html OK"
```

Notes:
- `--mem` is **track-only** in drac-local (recorded, not enforced); real OOM
  behaviour is validated by the host `deepen` tests, not here.
- The exit-75 -> requeue resume path uses `scontrol requeue` on a real cluster; in
  the emulator a plain re-`sbatch slurm/quiverlab-drac.sbatch tiny.yaml result.json`
  resumes from the same checkpoint dir under `$SCRATCH`.
- `QUIVERLAB_LAUNCH` overrides the auto-detection entirely, e.g.
  `QUIVERLAB_LAUNCH="python -m quiverlab.hpc" sbatch slurm/quiverlab.sbatch ...`.
