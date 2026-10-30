"""Check manifest metadata generation preserves the existing detector scripts."""

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


@pytest.mark.parametrize("entry", ENTRIES, ids=lambda entry: Path(entry["script"]).stem)
def test_metadata_matches_existing_script(tmp_path, entry):
    generator = Path(entry["script"]).stem
    original = tmp_path / "original"
    original.mkdir()
    subprocess.run(
        [
            sys.executable,
            str(ROOT / f"Detray/detectors/python/{generator}.py"),
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
        [sys.executable, str(ROOT / entry["script"]), str(output)],
        cwd=tmp_path,
        env={
            **os.environ,
            "PYTHONPATH": str(ROOT / "Detray/python"),
            "PYTHONHASHSEED": "1",
        },
        check=True,
    )
    assert output.read_bytes() == (original / f"{generator}.hpp").read_bytes()
