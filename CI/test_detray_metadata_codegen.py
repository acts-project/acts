"""Check metadata output interfaces and deterministic generation."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parent.parent
ENTRIES = [
    entry
    for entry in json.loads((ROOT / "codegen/manifest.json").read_text())[
        "units"
    ].values()
    if entry["script"].startswith("Detray/detectors/python/")
]


@pytest.mark.parametrize("entry", ENTRIES, ids=lambda entry: entry["arguments"][-1])
def test_metadata_output_interfaces_match(tmp_path, entry):
    detector = Path(entry["output"]).stem.removesuffix("_metadata")
    original = tmp_path / "original"
    original.mkdir()
    subprocess.run(
        [
            sys.executable,
            str(ROOT / entry["script"]),
            *entry["arguments"],
            "--no-format",
            "--output",
            str(original),
        ],
        env={
            **os.environ,
            "PYTHONPATH": str(ROOT / "Detray/python"),
            "PYTHONHASHSEED": "0",
        },
        check=True,
    )
    # The common generator contract accepts any output filename, not just the
    # detector script's implicit <detector>_metadata.hpp name.
    output = tmp_path / "generated" / "renamed.hpp"
    subprocess.run(
        [sys.executable, str(ROOT / entry["script"]), str(output), *entry["arguments"]],
        cwd=tmp_path,
        env={
            **os.environ,
            "PYTHONPATH": str(ROOT / "Detray/python"),
            "PYTHONHASHSEED": "1",
        },
        check=True,
    )
    assert output.read_bytes() == (original / f"{detector}_metadata.hpp").read_bytes()


@pytest.mark.parametrize("function", ["configure", "missing"])
def test_external_metadata(tmp_path, function):
    configuration = tmp_path / "custom.py"
    configuration.write_text(
        "from detray.detectors.configurations.toy import add_toy_types\n"
        "def configure(md):\n"
        "    add_toy_types(md)\n"
    )
    output = tmp_path / "custom_metadata.hpp"
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "Detray/detectors/python/generate_metadata.py"),
            str(output),
            "--detector",
            f"{configuration}:{function}",
        ],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(ROOT / "Detray/python")},
        capture_output=True,
        text=True,
    )
    if function == "missing":
        assert result.returncode != 0
        assert "Invalid detector configuration" in result.stderr
    else:
        assert result.returncode == 0, result.stderr
        assert "custom_metadata" in output.read_text()
