"""Client-IP hashing, HTTP security headers, ULID validation, and error-type
sanitisation.

Raw client addresses are never stored: only a salted SHA-256 goes into the jobs
table, so per-IP rate limiting works without retaining a personal identifier.

Error payloads surface ``QuiverlabError`` type names + messages verbatim (and the
runner's own honest domain tags), but any UNEXPECTED exception type name -- the
runner leaks ``type(exc).__name__`` + ``str(exc)`` into ``RunError`` for
server-log fidelity -- is GENERICISED here so no internal class name or stray
message reaches a client."""
from __future__ import annotations

import hashlib
from functools import lru_cache

# Strict CSP: self-only, no inline script. All JS is vendored/static (Task 11),
# so 'unsafe-inline' is never needed for scripts. KaTeX CSS is self-hosted.
CSP = ("default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; "
       "img-src 'self' data:; font-src 'self'; connect-src 'self'; "
       "object-src 'none'; base-uri 'self'; frame-ancestors 'none'")

SECURITY_HEADERS = {
    "Content-Security-Policy": CSP,
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "DENY",
    "Referrer-Policy": "no-referrer",
}


def hash_ip(ip: str, salt: str) -> str:
    """Salted SHA-256 of a client address, stored instead of the raw address."""
    return hashlib.sha256((salt + ":" + ip).encode("utf-8")).hexdigest()


# --------------------------------------------------------------------------- #
# ULID validation (Crockford base32) -- every job-id route validates the id
# BEFORE touching the store or the filesystem.
# --------------------------------------------------------------------------- #

# Crockford's base32 excludes I, L, O, U to avoid transcription ambiguity.
_CROCKFORD = frozenset("0123456789ABCDEFGHJKMNPQRSTVWXYZ")


def valid_ulid(s: str) -> bool:
    """True iff ``s`` is a 26-char Crockford-base32 string (a well-formed ULID).

    Rejects the wrong length and any out-of-alphabet character -- so a traversal
    payload (``../``, ``%2e``) or an injection attempt never reaches the store or
    an artifact path."""
    return len(s) == 26 and all(c in _CROCKFORD for c in s.upper())


# --------------------------------------------------------------------------- #
# Error-type sanitisation
# --------------------------------------------------------------------------- #

GENERIC_ERROR_TYPE = "InternalError"
GENERIC_ERROR_MESSAGE = "unexpected error; the incident was logged"

# The runner tags a handful of honest, client-safe refusals with its OWN names
# (these are NOT QuiverlabError subclasses); they surface verbatim alongside the
# library error names. Kept in sync with every ``RunError("<tag>", ...)`` the
# runner raises directly.
_RUNNER_TAGS = frozenset({
    "CatalogError", "SchemaError", "FieldError",
    "ResultTooLarge", "DuplicateComputeItem",
})


@lru_cache(maxsize=1)
def _quiverlab_error_names() -> frozenset[str]:
    """Names of every ``QuiverlabError`` subclass (imported lazily so the module
    carries no import-time dependency on quiverlab for the non-error paths)."""
    from quiverlab import errors as qerr

    names: set[str] = set()
    stack = [qerr.QuiverlabError]
    while stack:
        cls = stack.pop()
        names.add(cls.__name__)
        stack.extend(cls.__subclasses__())
    return frozenset(names)


def is_safe_error_type(error_type: str) -> bool:
    """True when ``error_type`` is a ``QuiverlabError`` subclass name or one of
    the runner's own honest tags -- i.e. safe to surface verbatim. Anything else
    is an unexpected internal exception name and must be genericised."""
    return error_type in _RUNNER_TAGS or error_type in _quiverlab_error_names()


def sanitize_error(error_type: str, message: str) -> tuple[str, str]:
    """Return a client-safe ``(error_type, message)``.

    A safe error type passes through verbatim; an unexpected internal error is
    replaced with a generic type + message that leaks neither the class name nor
    the (possibly sensitive) exception string."""
    if is_safe_error_type(error_type):
        return error_type, message
    return GENERIC_ERROR_TYPE, GENERIC_ERROR_MESSAGE


def sanitize_error_string(raw: str | None) -> str | None:
    """Client-safe form of a STORED job error string (``"Type: message"``).

    The worker stores the raw ``f"{type(exc).__name__}: {exc}"`` (and the
    parent's ``"worker error: <Type>: ..."``) for server-side forensics; this is
    the READ boundary that genericises it, so the async/queued path leaks no more
    than the sync path's ``sanitize_error``. A safe type (``QuiverlabError``
    subclass name or a runner tag) passes VERBATIM -- honest refusals like
    ``RelationError``/``ResultTooLarge``/``DepthLimitError`` still show their type
    and message; an unexpected internal type becomes the same generic
    ``InternalError`` type+message the sync path returns.

    The type is the substring before the FIRST ``": "`` (``str.partition`` splits
    once, so a message that itself contains ``": "`` is preserved intact); a
    string with no ``": "`` cannot be vetted and is genericised."""
    if raw is None:
        return None
    error_type, sep, message = raw.partition(": ")
    if not sep:                          # no "Type: message" shape -> cannot vet
        error_type, message = raw, ""
    etype, msg = sanitize_error(error_type, message)
    return f"{etype}: {msg}" if msg else etype
