"""Macaulay2 oracle bridge (subprocess). Importing this package does NOT
probe the M2 binary -- everything is lazy, mirroring quiverlab.qpa. The typed
``crosscheck`` dispatcher is added in Task 3 (a lazy wrapper, like qpa)."""
from quiverlab.m2.session import m2_available, require_m2, should_skip_m2

__all__ = ["m2_available", "require_m2", "should_skip_m2"]
