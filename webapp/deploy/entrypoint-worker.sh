#!/usr/bin/env bash
# Worker service entrypoint: run cfg.worker_processes concurrent poll loops in
# one container (spawn context). exec so signals (SIGTERM from `compose stop`)
# reach Python directly and the daemon loops exit with it.
set -euo pipefail
exec python -m webapp.worker.run_loop
