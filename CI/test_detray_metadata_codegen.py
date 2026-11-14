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


@pytest.mark.parametrize(
    ("selector", "prebuilt", "error"),
    [
        ("itk_detector.py:configure", True, None),
        ("itk.py:configure", True, None),
        ("other.py:configure", True, "must select an external ITk configuration"),
        ("itk_detector.py:configure", False, "No pre-generated file"),
    ],
)
def test_traccc_metadata_prebuilt_install(tmp_path, selector, prebuilt, error):
    """Exercise the real target, strict lookup, and installation without C++ deps."""
    source = tmp_path / "source"
    source.mkdir()
    prebuilts = tmp_path / "prebuilt"
    prebuilts.mkdir()
    header = prebuilts / "Detray/detectors/detray/detectors/itk_metadata.hpp"
    if prebuilt:
        header.parent.mkdir(parents=True)
        header.write_text("// prebuilt ITk metadata\n")
    (source / "CMakeLists.txt").write_text(f"""cmake_minimum_required(VERSION 3.19)
project(TracccMetadataCheck NONE)
set(CMAKE_INSTALL_INCLUDEDIR include)
set(CMAKE_INSTALL_LIBDIR lib)
set(CMAKE_INSTALL_BINDIR bin)
include("{ROOT}/cmake/ActsCodegen.cmake")
include("{ROOT}/Detray/cmake/detray-functions.cmake")
include("{ROOT}/Traccc/cmake/traccc-functions.cmake")
add_library(detray::detectors INTERFACE IMPORTED)
add_subdirectory("{ROOT}/Traccc/metadata" metadata)
install(EXPORT traccc-exports NAMESPACE traccc:: DESTINATION lib/cmake/traccc)
""")
    build = tmp_path / "build"
    result = subprocess.run(
        [
            "cmake",
            "-G",
            "Ninja",
            "-S",
            str(source),
            "-B",
            str(build),
            f"-DACTS_CODEGEN_PREBUILT_DIR={prebuilts}",
            "-DACTS_CODEGEN_REQUIRE_PREBUILT=ON",
            f"-DTRACCC_ITK_METADATA_CONFIGURATION={tmp_path / selector}",
        ],
        capture_output=True,
        text=True,
    )
    if error:
        assert result.returncode != 0
        assert error in result.stderr
        return
    assert result.returncode == 0, result.stderr
    subprocess.run(["cmake", "--build", str(build)], check=True)
    prefix = tmp_path / "install"
    subprocess.run(
        ["cmake", "--install", str(build), "--prefix", str(prefix)], check=True
    )
    assert (
        prefix / "include/detray/detectors/itk_metadata.hpp"
    ).read_bytes() == header.read_bytes()
    exports = (prefix / "lib/cmake/traccc/traccc-exports.cmake").read_text()
    assert "traccc::detray_metadata" in exports
    assert "detray::detectors" in exports
