"""Plan 54 Task D -- the bracket recovered from Delta via the BV relation, and the
in-window ARBITER against the independent Gerstenhaber bracket.

(Task B staging: ``attach_bracket_arbiter`` is a no-op until Task D lands the real
arbiter -- the transport calls it so the wiring exists.)
"""


def attach_bracket_arbiter(A, bv, top, max_cells):
    """Attach ``bv.derived_bracket`` (built from Delta via the (BV) relation) and
    ``bv.bracket_check`` (the in-window derived == independent Gerstenhaber
    equality). A mismatch is a loud refusal. (Task D fills this in.)"""
    return bv
