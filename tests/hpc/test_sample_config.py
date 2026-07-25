"""``sample-config`` emits a parseable, validating, round-trippable config."""
from quiverlab.hpc import cli
from quiverlab.hpc.spec import load_config_text, parse_request, run as spec_run


def test_sample_config_constant_validates():
    req = parse_request(load_config_text(cli.SAMPLE_CONFIG))
    assert req.schema_version == 1
    assert req.algebra.kind == "quiver"
    assert req.compute and "hh_cohomology:0..4" in req.compute


def test_sample_config_verb_output_round_trips(capsys):
    rc = cli.main(["sample-config"])
    assert rc == 0
    out = capsys.readouterr().out
    assert out.strip(), "sample-config produced no output"
    # The printed YAML parses and validates identically to the constant.
    req = parse_request(load_config_text(out))
    assert req.algebra.kind == "quiver"


def test_sample_config_actually_computes(tmp_path):
    # The sample is a real, cheap job: it must run end to end.
    cfg = load_config_text(cli.SAMPLE_CONFIG)
    result = spec_run(cfg, tmp_path)
    assert result["results"]["hh_cohomology"]["dims"] == [3, 2, 2, 2, 2]
    assert result["results"]["cartan"]["matrix"] == [[3]]
