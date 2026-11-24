import json
import subprocess
from pathlib import Path

import pytest

from jury.tracy import basecall, parse_basecalls, quality_window


def test_actual_upstream_json_fields_and_longest_quality_window():
    seq, qualities, bounds = parse_basecalls({
        "primarySeq": "ACGTACGT", "secondarySeq": "ACGTACGT",
        "basecallQual": [10, 30, 30, 30, 5, 30, 30, 10],
        "basecallPos": [12, 24, 36, 48, 60, 72, 84, 96],
    })
    assert seq == "ACGTACGT"
    assert qualities == (10, 30, 30, 30, 5, 30, 30, 10)
    assert bounds == (2, 4)


def test_ambiguous_base_breaks_quality_window():
    assert quality_window("AACNACGT", [30] * 8) == (5, 8)
    assert quality_window("ACGT", [0] * 4) == (1, 0)


@pytest.mark.parametrize("data", [
    {}, {"primarySeq": "ACGT", "basecallQual": [30]},
    {"primarySeq": "ACGT", "basecallQual": [30, 30, 30, "30"]},
    {"primarySeq": "ACGT", "basecallQual": [30, 30, 30, True]},
    {"primarySeq": "AC-T", "basecallQual": [30] * 4},
])
def test_malformed_tracy_json_rejected(data):
    with pytest.raises(ValueError):
        parse_basecalls(data)


def test_subprocess_argument_list_and_temporary_cleanup(tmp_path, monkeypatch):
    trace = tmp_path / "sample with spaces.ab1"
    trace.write_bytes(b"synthetic placeholder")
    outputs = []

    def fake_run(command, **options):
        assert command[:4] == ["tracy", "basecall", "-f", "json"]
        assert command[-1] == str(trace.resolve())
        assert options["check"] is True
        assert options["capture_output"] is True
        assert options["timeout"] == 120
        output = Path(command[command.index("-o") + 1])
        outputs.append(output)
        output.write_text(json.dumps({"primarySeq": "ACGT", "basecallQual": [30] * 4}), encoding="utf-8")

    monkeypatch.setattr("jury.tracy.subprocess.run", fake_run)
    assert basecall(trace, executable="tracy") == ("ACGT", (30, 30, 30, 30), (1, 4))
    assert not outputs[0].exists()
    assert not outputs[0].parent.exists()


def test_subprocess_error_is_clear_and_does_not_echo_private_output(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(2, args[0], stderr="private sequence")

    monkeypatch.setattr("jury.tracy.subprocess.run", fail)
    with pytest.raises(ValueError, match="exit 2") as error:
        basecall(tmp_path / "synthetic.ab1", executable="tracy")
    assert "private sequence" not in str(error.value)
