#!/bin/bash
# quiverlab generic SLURM batch template (Plan 28).
#
# Runs one quiverlab HPC job -- `quiverlab-hpc run <config> -o <result>` -- with
# atomic per-degree checkpoints. If the compute checkpoints out (exit 75, e.g. it
# ran out of wall time), this script resumes it: a plain re-`sbatch` picks up from
# the checkpoint dir, and the USR1 handler below requeues automatically.
#
# EDIT the #SBATCH resource lines for your cluster, then:
#     sbatch slurm/quiverlab.sbatch my-config.yaml result.json
#
# Launcher precedence (auto-detected, override with $QUIVERLAB_LAUNCH):
#     $QUIVERLAB_LAUNCH  ->  apptainer  ->  singularity  ->  python -m quiverlab.hpc
# so the SAME script drives a SIF on a cluster and the pip-installed wheel in a
# plain venv (e.g. the drac-local emulator, which ships no apptainer).
#
# =============================================================================
# CPU-ONLY -- do NOT request GPUs (no --gres/--gpus, no --nv). quiverlab computes
# with EXACT arithmetic (integers mod p, rationals) and has no GPU kernels, so a
# GPU allocation would sit idle and waste your account. Ask for CORES
# (--cpus-per-task) and RAM (--mem) instead. Size them with:
#     quiverlab-hpc estimate my-config.yaml
# which prints a suggested --time / --cpus-per-task / --mem for the example.
# =============================================================================
#
#SBATCH --job-name=quiverlab
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=quiverlab-%j.out
#SBATCH --error=quiverlab-%j.err
# Deliver USR1 to THIS batch script (the B: prefix) 180s before the wall-clock
# limit so the compute can write a clean checkpoint and exit 75 (handler below).
#SBATCH --signal=B:USR1@180
#
# ---- Array sweep variant (uncomment; one config per task) -------------------
# #SBATCH --array=0-9
# Then drive one config per index instead of the positional args below:
#     CONFIG="configs/sweep-${SLURM_ARRAY_TASK_ID}.yaml"
#     RESULT="results/result-${SLURM_ARRAY_TASK_ID}.json"
# -----------------------------------------------------------------------------

set -euo pipefail

# --- inputs (positional, with defaults) --------------------------------------
CONFIG="${1:-quiverlab.yaml}"
RESULT="${2:-result.json}"
IMAGE="${QUIVERLAB_IMAGE:-quiverlab.sif}"

# --- scratch + checkpoint dir ------------------------------------------------
# $SCRATCH is defined on real clusters; DEFAULT IT FIRST (`:=`) so `set -u` never
# trips off-cluster (drac-local emulator / laptop). Every later use is now bound.
: "${SCRATCH:=$PWD/.quiverlab-scratch}"
mkdir -p "$SCRATCH"
CKPT="$SCRATCH/quiverlab-ckpt/${SLURM_JOB_ID:-local}"
mkdir -p "$CKPT"

# Keep Apptainer's layer cache off your (small) $HOME quota.
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$SCRATCH/.apptainer}"

# Thread caps follow the allocation (the CLI does this too; belt and braces).
export NUMBA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

# --- launcher auto-detect ----------------------------------------------------
if [ -n "${QUIVERLAB_LAUNCH:-}" ]; then
    # Explicit override, e.g.
    #   QUIVERLAB_LAUNCH="apptainer run /path/quiverlab.sif"
    #   QUIVERLAB_LAUNCH="docker run --rm -v $PWD:$PWD -w $PWD ghcr.io/marcoarmenta/quiverlab:latest"
    read -r -a LAUNCH <<< "$QUIVERLAB_LAUNCH"
elif command -v apptainer >/dev/null 2>&1; then
    LAUNCH=(apptainer run --bind "$PWD,$SCRATCH" "$IMAGE")
elif command -v singularity >/dev/null 2>&1; then
    LAUNCH=(singularity run --bind "$PWD,$SCRATCH" "$IMAGE")
else
    # No container runtime: the pip-installed wheel (the emulator/venv fallback).
    LAUNCH=(python -m quiverlab.hpc)
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

# Re-issue wait until the child truly exits: a trapped USR1 makes `wait` return
# early (128+n) while the child is still alive.
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
