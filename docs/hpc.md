# Run quiverlab on your HPC cluster

quiverlab ships as **one container** you send to a SLURM cluster, run with a tiny
YAML config, and pull a result back from. Very large examples -- the ones a laptop
or the web tier will not finish -- become reachable, with **atomic checkpoints** so
a job that runs out of wall time simply resumes on the next submit.

You need only `ssh`, `scp`, and `sbatch`. No admin rights: Apptainer is rootless.

---

## 5 steps

### 1. Pull the image (once) and copy it up

On a machine with internet, pull the image into a single `.sif` file and `scp` it to
your cluster's `$SCRATCH` (Apptainer is rootless -- no `sudo`, no daemon):

```bash
apptainer pull quiverlab.sif docker://ghcr.io/marcoarmenta/quiverlab:latest
scp quiverlab.sif you@cluster:/scratch/you/
```

On the cluster you can instead pull directly (set `APPTAINER_CACHEDIR` to keep the
layer cache off your small `$HOME` quota):

```bash
export APPTAINER_CACHEDIR=$SCRATCH/.apptainer
apptainer pull $SCRATCH/quiverlab.sif docker://ghcr.io/marcoarmenta/quiverlab:latest
```

### 2. Write a config

A config is the same request shape as the web GUI, as YAML, plus an `hpc:` block.
Get an annotated starting point, or design the example in the GUI and export it:

```bash
apptainer run quiverlab.sif sample-config > my-config.yaml   # annotated template
# ...or, in the browser GUI, click "Config (YAML)" to download exactly your example.
```

Edit `my-config.yaml` to name your algebra, field, and what to compute.
The full config reference -- anatomy, every compute kind, module blocks, and
worked examples from characteristic 0 to explicit quivers with matrices -- is in
the README section [Writing and running config files](https://github.com/MarcoArmenta/quiverlab#writing-and-running-config-files-the-containerized-app)
and in [`container/examples/`](https://github.com/MarcoArmenta/quiverlab/tree/main/container/examples).

### 3. Submit

Copy a template from [`slurm/`](https://github.com/MarcoArmenta/quiverlab/tree/main/slurm),
edit the resource lines (and, on the Alliance, `--account`), then:

```bash
sbatch slurm/quiverlab-drac.sbatch my-config.yaml result.json   # DRAC / Alliance
# or the generic template:
sbatch slurm/quiverlab.sbatch my-config.yaml result.json
```

The template auto-detects the launcher (`apptainer` -> `singularity` ->
`python -m quiverlab.hpc`), binds `$PWD` and `$SCRATCH`, and puts the checkpoint dir
on `$SCRATCH`.

### 4. Pull the result back

```bash
scp you@cluster:/scratch/you/result.json .
```

`result.json` is the exact result envelope: dimensions, references, metadata, and a
`result_schema` stamp. (Heavy payloads like bases are written to a sidecar only when
you pass `--allow-large`.)

### 5. Render a report locally

`result.json` renders to a self-contained, print-ready HTML report (or plain text)
with the **same** container -- pick whichever runtime you have:

```bash
# Apptainer:
apptainer run --bind "$PWD" quiverlab.sif render result.json -o report.html
# Docker:
docker run --rm -v "$PWD:/w" -w /w ghcr.io/marcoarmenta/quiverlab:latest render result.json -o report.html
# Plain pip install (no container):
pip install "quiverlab[hpc]" && quiverlab-hpc render result.json -o report.html
```

The HTML report is print-ready (export to PDF via the browser's Print → Save as
PDF); `--format json` emits the worked-steps event stream (`trace.json`). PDF/TeX
report output has been removed. The output is fully self-contained;
`--format txt` always works.

---

## What resources to request

**quiverlab is CPU-only. Never request a GPU** (`--gres=gpu`, `--gpus`, or
`apptainer run --nv`). Every computation is exact -- integers modulo a prime, or
rationals -- with no floating point and no GPU kernel anywhere, so a GPU allocation
would sit idle and waste your account. Ask for two things instead:

- **cores** -- `--cpus-per-task`. The container caps its thread pools
  (`NUMBA_NUM_THREADS`/`OMP_NUM_THREADS`) to the cores it is given.
- **RAM** -- `--mem`. Deep resolutions are memory-bound; this is usually the
  binding constraint.

Let the estimator size both for a given example:

```bash
apptainer run quiverlab.sif estimate my-config.yaml
# prints the tier and a suggested --time / --cpus-per-task / --mem
```

The container is **cgroup- and SLURM-aware**: it auto-detects the cores and memory
it has actually been granted (from `$SLURM_CPUS_PER_TASK` and the cgroup limits) and
sizes its thread pools and memory budget accordingly. If a GPU happens to be visible
it is reported as *detected but unused* -- never engaged.

---

## Resume: a job that runs out of time

Big jobs checkpoint atomically per degree onto `$SCRATCH`. If a job hits the wall
limit, the template's `--signal=B:USR1@180` handler asks the compute to write a
clean checkpoint and exit **75**; the script then requeues itself. You do not have
to do anything -- but if a job was killed hard, **just submit the same command
again**:

```bash
sbatch slurm/quiverlab-drac.sbatch my-config.yaml result.json   # resumes from $SCRATCH
```

`deepen` reads the checkpoint dir and continues from the last completed degree; the
final result is identical to an uninterrupted run.

---

## Array sweeps

To sweep a family of examples, uncomment the array block in the template and drive
one config per task index:

```bash
#SBATCH --array=0-9
CONFIG="configs/sweep-${SLURM_ARRAY_TASK_ID}.yaml"
RESULT="results/result-${SLURM_ARRAY_TASK_ID}.json"
```

Each task checkpoints independently under `$SCRATCH/quiverlab-ckpt/<job-id>`.

---

## Singularity-only clusters

Older sites ship **Singularity** rather than Apptainer. The commands are identical
with `singularity` in place of `apptainer`, and the sbatch templates auto-detect it.
Pull with `singularity pull quiverlab.sif docker://ghcr.io/marcoarmenta/quiverlab:latest`.

## Air-gapped clusters

If the compute nodes cannot reach GHCR, download the **SIF release asset** attached
to each [GitHub release](https://github.com/MarcoArmenta/quiverlab/releases) on an
internet-connected host and `scp` it up. It is a single self-contained file with the
seeded example cache baked in -- nothing else to fetch.

## `APPTAINER_CACHEDIR`

Always point `APPTAINER_CACHEDIR` at `$SCRATCH` (the templates do this for you):
`$HOME` quotas on shared clusters are small, and Apptainer's layer cache will blow
through them otherwise.

---

## Hosting the web service on DRAC cloud (aside)

This page is about the **batch** container. If instead you want to self-host the
**web tier** (the queued/big-job server) on a DRAC RAS cloud instance, choose a
**CPU/RAM-heavy flavor and never a GPU flavor** -- the same CPU-only reasoning
applies. See [`webapp/deploy/PROVISIONING.md`](https://github.com/MarcoArmenta/quiverlab/blob/main/webapp/deploy/PROVISIONING.md).

## Offline laptop app

The same image is also a **zero-code offline app** -- open the GUI in your browser
with no internet. See [Offline laptop app](offline-app.md).
