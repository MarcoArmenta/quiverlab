"""Plan 28 -- the "reproduce on a cluster" YAML.

Serialises a validated compute request (schema 1/2, exactly the shape
``ComputeRequest`` accepts) into a YAML config a user can drop straight into::

    quiverlab-hpc run this-file.yaml -o result.json

Kept in its own module (NOT ``runner.py``) so the render helpers stay import-light
and the runner's dispatch is untouched. The ``docs/gui/gui.js`` ``configYaml()``
emitter MIRRORS this: both must round-trip through a YAML parser back to the same
compute request (the export test pins the round-trip, not byte-identity). Change
one, change the other.
"""
from __future__ import annotations

import json

# The one-line runnable header. app.js/gui.js emit the same instruction verbatim
# (see docs/gui/gui.js el.config handler) so a downloaded file and a server-shown
# block carry identical guidance.
HEADER = "# quiverlab-hpc run this-file.yaml -o result.json\n"


def cluster_config_yaml(spec: dict) -> str:
    """Return the annotated YAML for ``spec`` (a request dict). The spec is first
    normalised through JSON so arrow *tuples* (from ``model_dump``) become lists
    and nothing exotic reaches the YAML dumper; keys are sorted for a stable,
    diff-friendly file. The result parses back (comments ignored) into a
    schema-valid request."""
    import yaml

    normalized = json.loads(json.dumps(spec, default=str))
    body = yaml.safe_dump(normalized, sort_keys=True, default_flow_style=False,
                          allow_unicode=True)
    return HEADER + body
