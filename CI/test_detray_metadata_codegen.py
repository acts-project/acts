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
