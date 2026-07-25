"""``python -m quiverlab.hpc`` -- the same entry point as the ``quiverlab-hpc``
console script."""
import sys

from quiverlab.hpc.cli import main

if __name__ == "__main__":
    sys.exit(main())
