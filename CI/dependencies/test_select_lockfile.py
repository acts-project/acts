#!/usr/bin/env python3
"""Self-tests for select_lockfile.py.

Guards asset-name parsing and the selection rules. Picking the wrong lockfile
does not fail here: it fails minutes later in a configure step, on a runner.

Run directly; non-zero exit on failure. Also collectable by pytest.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Optional

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import select_lockfile as sl  # noqa: E402

# Real asset names from v24.0.0, the first release carrying flavors, in the
# shape the GitHub releases API returns them.
V24_ASSETS = [
    "spack_darwin-tahoe-aarch64.lock",
    "spack_darwin-tahoe-aarch64_apple-clang@17.0.0_cxx23.lock",
    "spack_linux-almalinux9-x86_64.lock",
    "spack_linux-almalinux9-x86_64_gcc@15.2.1_cxx23.lock",
    "spack_linux-ubuntu26.04-aarch64_gcc@15.2.0_cxx23.lock",
    "spack_linux-ubuntu26.04-x86_64.lock",
    "spack_linux-ubuntu26.04-x86_64_gcc@15.2.0_cxx20.lock",
    "spack_linux-ubuntu26.04-x86_64_gcc@15.2.0_cxx23.lock",
    "spack_linux-ubuntu26.04-x86_64_gcc@15.2.0_cxx23_cuda13.lock",
    "spack_linux-ubuntu26.04-x86_64_gcc@15.2.0_cxx23_rocm-gfx90a.lock",
    "spack_linux-ubuntu26.04-x86_64_llvm@22.1.8_cxx23.lock",
]

# v23.3.1, pre-flavor, so today's host selections stay covered too.
V23_ASSETS = [
    "spack_linux-ubuntu24.04-x86_64.lock",
    "spack_linux-ubuntu24.04-x86_64_gcc@13.3.0_cxx20.lock",
    "spack_linux-ubuntu24.04-x86_64_llvm@22.1.7_cxx20.lock",
    "spack_linux-ubuntu24.04-x86_64_llvm@22.1.7_cxx23.lock",
]

UBUNTU26 = "linux-ubuntu26.04-x86_64"
UBUNTU24 = "linux-ubuntu24.04-x86_64"


def _release(names: list[str]) -> Dict:
    return {
        "assets": [
            {"name": n, "browser_download_url": f"https://example.invalid/{n}"}
            for n in names
        ]
    }


def _select(
    names: list[str],
    arch: str,
    compiler: Optional[str],
    cxx: str = "cxx20",
    flavor: Optional[str] = None,
) -> str:
    """Select against a synthetic release, returning the chosen asset name."""
    lockfiles = sl.parse_assets(_release(names))
    url = sl.select_lockfile(lockfiles, arch, compiler, cxx, flavor)
    return url.rsplit("/", 1)[-1]


# --- flavor token extraction -----------------------------------------------


def test_extract_flavor_host_specs():
    assert sl.extract_flavor("gcc@15.2.0_cxx23") is None
    assert sl.extract_flavor("gcc@13.3.0") is None, "pre-cxx spec has no flavor"
    assert sl.extract_flavor("default") is None


def test_extract_flavor_accelerator_specs():
    assert sl.extract_flavor("gcc@15.2.0_cxx23_cuda13") == "cuda13"
    assert sl.extract_flavor("gcc@15.2.0_cxx23_rocm-gfx90a") == "rocm-gfx90a"


def test_extract_flavor_keeps_underscores_in_flavor_name():
    # Not shipped today, but the grammar permits it and take-last would truncate.
    assert sl.extract_flavor("gcc@15.2.0_cxx23_cuda13_sm90") == "cuda13_sm90"


def test_normalize_flavor():
    assert sl.normalize_flavor("host") is None, "'host' is the plain CPU stack"
    assert sl.normalize_flavor("") is None
    assert sl.normalize_flavor(None) is None
    assert sl.normalize_flavor("  cuda13 ") == "cuda13"


# --- flavor selection -------------------------------------------------------


def test_flavored_lockfile_is_selected():
    got = _select(V24_ASSETS, UBUNTU26, "gcc@15.2.0", cxx="cxx23", flavor="cuda13")
    assert got == "spack_linux-ubuntu26.04-x86_64_gcc@15.2.0_cxx23_cuda13.lock", got


def test_flavor_name_with_dash_is_selected():
    got = _select(V24_ASSETS, UBUNTU26, "gcc@15.2.0", cxx="cxx23", flavor="rocm-gfx90a")
    assert (
        got == "spack_linux-ubuntu26.04-x86_64_gcc@15.2.0_cxx23_rocm-gfx90a.lock"
    ), got


def test_host_request_never_picks_a_flavored_lockfile():
    # Before flavors were parsed, these were ordinary candidates in the pool.
    for cxx in ("cxx20", "cxx23"):
        got = _select(V24_ASSETS, UBUNTU26, "gcc@15.2.0", cxx=cxx)
        assert sl.extract_flavor(got) is None, f"{cxx} picked flavored {got}"


def test_flavor_falls_back_across_cxx_but_not_across_flavor():
    # Flavors ship at cxx23 only, so cxx20+cuda13 must land on the cxx23 CUDA
    # stack -- never the cxx20 *host* stack, the wrong-but-plausible answer.
    got = _select(V24_ASSETS, UBUNTU26, "gcc@15.2.0", cxx="cxx20", flavor="cuda13")
    assert got == "spack_linux-ubuntu26.04-x86_64_gcc@15.2.0_cxx23_cuda13.lock", got


def test_unknown_flavor_exits_rather_than_falling_back():
    # Returning the host stack would surface as a missing-CUDA configure error
    # much later, pointing nowhere near the cause.
    try:
        _select(V24_ASSETS, UBUNTU26, "gcc@15.2.0", cxx="cxx23", flavor="cuda12")
    except SystemExit as e:
        assert e.code == 1, e.code
    else:
        raise AssertionError("expected a hard failure for an unavailable flavor")


def test_flavor_request_on_arch_without_flavors_exits():
    try:
        _select(V23_ASSETS, UBUNTU24, "gcc@13.3.0", cxx="cxx20", flavor="cuda13")
    except SystemExit as e:
        assert e.code == 1, e.code
    else:
        raise AssertionError("expected a hard failure on a flavorless release")


# --- pre-existing selection rules, unchanged by the flavor axis --------------


def test_exact_compiler_and_cxx_match():
    got = _select(V24_ASSETS, UBUNTU26, "gcc@15.2.0", cxx="cxx20")
    assert got == "spack_linux-ubuntu26.04-x86_64_gcc@15.2.0_cxx20.lock", got


def test_clang_is_an_alias_for_llvm():
    got = _select(V24_ASSETS, UBUNTU26, "clang@22.1.8", cxx="cxx23")
    assert got == "spack_linux-ubuntu26.04-x86_64_llvm@22.1.8_cxx23.lock", got


def test_unknown_compiler_version_uses_highest_of_family():
    got = _select(V23_ASSETS, UBUNTU24, "llvm@22.1.1", cxx="cxx23")
    assert got == "spack_linux-ubuntu24.04-x86_64_llvm@22.1.7_cxx23.lock", got


def test_no_compiler_uses_the_arch_default():
    got = _select(V24_ASSETS, UBUNTU26, None, cxx="cxx20")
    assert got == "spack_linux-ubuntu26.04-x86_64.lock", got


def test_v23_host_selection_is_unchanged():
    # In production today: whatever else moves, this must not.
    got = _select(V23_ASSETS, UBUNTU24, "gcc@13.3.0", cxx="cxx20")
    assert got == "spack_linux-ubuntu24.04-x86_64_gcc@13.3.0_cxx20.lock", got


def _main() -> int:
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok   {name}")
            except AssertionError as e:
                failures += 1
                print(f"FAIL {name}: {e}")
    print(f"\n{failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())
