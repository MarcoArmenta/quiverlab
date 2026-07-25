#!/bin/bash
# quiverlab SLURM batch template for the Digital Research Alliance of Canada
# (DRAC / Alliance: Narval, Fir, Rorqual, Nibi, ...). Plan 28.
#
# Same behaviour as slurm/quiverlab.sbatch (auto-detected launcher, USR1 ->
# checkpoint -> requeue resume) with the DRAC-specific header: an --account line
# you MUST edit, and `module load apptainer`.
#
#     # one-time: pull the image onto $SCRATCH (rootless, no admin)
#     export APPTAINER_CACHEDIR=$SCRATCH/.apptainer
#     apptainer pull $SCRATCH/quiverlab.sif docker://ghcr.io/MarcoArmenta/quiverlab:latest
#     # then:
#     sbatch slurm/quiverlab-drac.sbatch my-config.yaml result.json
#
# =============================================================================
# CPU-ONLY -- do NOT request GPUs (no --gres=gpu, no --gpus, no --nv). quiverlab
# computes with EXACT arithmetic (integers mod p, rationals) and has no GPU
# kernels, so a GPU allocation would sit idle and burn your RAC/RAS allocation.
# Ask for CORES (--cpus-per-task) and RAM (--mem) instead. Size them with:
#     quiverlab-hpc estimate my-config.yaml
# =============================================================================
#
#SBATCH --account=def-CHANGEME
#SBATCH --job-name=quiverlab
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=quiverlab-%j.out
#SBATCH --error=quiverlab-%j.err
# Deliver USR1 to THIS batch script (B: prefix) 180s before the wall limit so the
# compute can checkpoint and exit 75 (handler below).
#SBATCH --signal=B:USR1@180
#
# ---- Array sweep variant (uncomment; one config per task) -------------------
# #SBATCH --array=0-9
# Then drive one config per index instead of the positional args below:
#     CONFIG="configs/sweep-${SLURM_ARRAY_TASK_ID}.yaml"
#     RESULT="results/result-${SLURM_ARRAY_TASK_ID}.json"
# -----------------------------------------------------------------------------

set -euo pipefail

# DRAC/Alliance clusters expose Apptainer as an environment module.
if command -v module >/dev/null 2>&1; then
    module load apptainer 2>/dev/null || true
fi

# --- scratch (define FIRST) --------------------------------------------------
# $SCRATCH is defined on every Alliance cluster; DEFAULT IT FIRST (`:=`) so
# `set -u` never trips off-cluster (drac-local emulator / laptop) AND so the
# IMAGE default below can dereference it safely.
: "${SCRATCH:=$PWD/.quiverlab-scratch}"
mkdir -p "$SCRATCH"

# --- inputs (positional, with defaults) --------------------------------------
CONFIG="${1:-quiverlab.yaml}"
RESULT="${2:-result.json}"
IMAGE="${QUIVERLAB_IMAGE:-$SCRATCH/quiverlab.sif}"

# --- checkpoint dir ----------------------------------------------------------
CKPT="$SCRATCH/quiverlab-ckpt/${SLURM_JOB_ID:-local}"
mkdir -p "$CKPT"

# Keep Apptainer's layer cache off your (small) $HOME quota.
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$SCRATCH/.apptainer}"

# Thread caps follow the allocation (the CLI does this too; belt and braces).
export NUMBA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

# --- launcher auto-detect ----------------------------------------------------
if [ -n "${QUIVERLAB_LAUNCH:-}" ]; then
    read -r -a LAUNCH <<< "$QUIVERLAB_LAUNCH"     # explicit override
elif command -v apptainer >/dev/null 2>&1; then
    LAUNCH=(apptainer run --bind "$PWD,$SCRATCH" "$IMAGE")
elif command -v singularity >/dev/null 2>&1; then
    LAUNCH=(singularity run --bind "$PWD,$SCRATCH" "$IMAGE")
else
    LAUNCH=(python -m quiverlab.hpc)              # emulator/venv fallback
fi

echo "[quiverlab] launcher: ${LAUNCH[*]}" >&2
echo "[quiverlab] config=$CONFIG result=$RESULT ckpt=$CKPT" >&2

# --- USR1 handler: forward to the compute for a clean checkpoint stop --------
CHILD=""
on_usr1() {
    echo "[quiverlab] USR1 received -- signalling compute to checkpoint" >&2
    [ -n "$CHILD" ] && kill -USR1 "$CHILD" 2>/dev/null || true
}
trap on_usr1 USR1

# --- run in the background so the trap fires promptly, then wait --------------
"${LAUNCH[@]}" run "$CONFIG" -o "$RESULT" --checkpoint-dir "$CKPT" &
CHILD=$!

STATUS=0
while true; do
    if wait "$CHILD"; then STATUS=0; else STATUS=$?; fi
    if kill -0 "$CHILD" 2>/dev/null; then
        continue        # still alive -> interrupted by the trap; keep waiting
    fi
    break               # child gone -> STATUS holds its real exit code
done

# --- exit-75 == clean checkpoint stop -> resume ------------------------------
if [ "$STATUS" -eq 75 ]; then
    echo "[quiverlab] checkpoint written (exit 75); resuming via requeue" >&2
    if [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then
        scontrol requeue "$SLURM_JOB_ID"
    elif command -v sbatch >/dev/null 2>&1; then
        sbatch "$0" "$CONFIG" "$RESULT"     # fallback: resubmit; deepen resumes from $CKPT
    fi
    exit 0
fi

echo "[quiverlab] finished with exit $STATUS" >&2
exit "$STATUS"
