"""Check metadata output interfaces and deterministic generation."""

import json
import os
import re
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
    if entry["script"] == "Detray/detectors/codegen/generate_metadata.py"
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
            "PYTHONPATH": str(ROOT / "Detray/detectors/codegen"),
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
            "PYTHONPATH": str(ROOT / "Detray/detectors/codegen"),
            "PYTHONHASHSEED": "1",
        },
        check=True,
    )
    assert output.read_bytes() == (original / f"{detector}_metadata.hpp").read_bytes()


@pytest.mark.parametrize("function", ["configure", "missing"])
def test_external_metadata(tmp_path, function):
    configuration = tmp_path / "custom_detector.py"
    configuration.write_text(
        "from detray_metadata.configurations import toy_detector\n"
        "def configure(md):\n"
        "    toy_detector.configure(md)\n"
    )
    output = tmp_path / "custom_metadata.hpp"
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "Detray/detectors/codegen/generate_metadata.py"),
            str(output),
            "--detector",
            f"{configuration}:{function}",
        ],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(ROOT / "Detray/detectors/codegen")},
        capture_output=True,
        text=True,
    )
    if function == "missing":
        assert result.returncode != 0
        assert "Invalid detector configuration" in result.stderr
    else:
        assert result.returncode == 0, result.stderr
        assert "struct custom_metadata" in output.read_text()


def test_shipped_configurations_match_manifest():
    """Every packaged configuration has a manifest unit and is always generated."""
    configurations = ROOT / "Detray/detectors/codegen/detray_metadata/configurations"
    names = sorted(
        path.stem.removesuffix("_detector")
        for path in configurations.glob("*_detector.py")
    )
    assert sorted(entry["arguments"][-1] for entry in ENTRIES) == names
    for entry in ENTRIES:
        assert (
            entry["output"] == f"detray/detectors/{entry['arguments'][-1]}_metadata.hpp"
        )
    cmake = (ROOT / "Detray/detectors/CMakeLists.txt").read_text()
    shipped = re.search(r"set\(_detray_shipped_metadata\s+([^)]*)\)", cmake)
    assert shipped is not None
    assert sorted(shipped.group(1).split()) == names


def test_traccc_itk_metadata_is_external(tmp_path):
    """traccc's ITk configuration is not shipped but generates as external."""
    assert (
        "Detray/detectors/detray/detectors/itk_metadata.hpp"
        not in json.loads((ROOT / "codegen/manifest.json").read_text())["units"]
    )
    output = tmp_path / "itk_metadata.hpp"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "Detray/detectors/codegen/generate_metadata.py"),
            str(output),
            "--detector",
            f"{ROOT / 'Traccc/metadata/itk_detector.py'}:configure",
        ],
        env={**os.environ, "PYTHONPATH": str(ROOT / "Detray/detectors/codegen")},
        check=True,
    )
    assert "struct itk_metadata" in output.read_text()
