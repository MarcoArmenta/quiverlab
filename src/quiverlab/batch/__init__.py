"""quiverlab.batch: SQLite persistence + scan driver for family sweeps.

Expert-facing; deliberately NOT exported at the quiverlab top level -- reach it via
``import quiverlab.batch``.  This package is the quiverlab lift of hanlab's ``labdb``
(Marco Armenta's HansConjecture, MIT): the pure ``analyze(spec) -> record`` unit of
work, the named family builder registry, and the ``ResultsDB`` SQLite store with the
Han counterexample-search queries.

  quiverlab.batch.db.ResultsDB       -- SQLite results store (context manager)
  quiverlab.batch.builders.BUILDERS  -- {name: builder}; build_algebra(spec)
  quiverlab.batch.scan.analyze       -- pure spec -> JSON record; run_scan(specs, n_workers)
  quiverlab.batch.PRIME              -- the char-0 proxy prime, 32003
"""
from quiverlab.batch import builders, db, scan  # noqa: F401
from quiverlab.batch.scan import PRIME  # noqa: F401
