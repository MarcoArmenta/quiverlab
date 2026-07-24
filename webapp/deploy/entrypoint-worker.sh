#!/usr/bin/env bash
# Worker service entrypoint: run cfg.worker_processes concurrent poll loops in
# one container (spawn context). `exec` so SIGTERM (from `compose stop`) reaches
# Python as PID 1 directly: run_loop.main() drains in-flight jobs, then exits.
# The compose `stop_grace_period` bounds that drain; any hard-killed job is
# requeued at the next startup (JobStore.requeue_stale_running).
set -euo pipefail
exec python -m webapp.worker.run_loop
