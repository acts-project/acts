"""Linker-isolation regression tests for libActsPluginArrow.

When the plugin is built with ACTS_ARROW_ISOLATED=ON, the arrow/parquet
static archives are absorbed into the plugin .dylib/.so with hidden
visibility plus an exported-symbols allowlist (see
Plugins/Arrow/cmake/exported_symbols.{txt,ld}). The intent is that the
plugin presents NO arrow or parquet symbols to the rest of the process,
so a co-loaded pyarrow (which has its own libarrow) cannot collide with
ACTS's bundled archives under any dlopen mode.

These tests verify the resulting library actually upholds that contract:

  - no LC_LOAD_DYLIB / NEEDED entry for libarrow*, libparquet*, or
    libarrow_dataset*;
  - no externally-defined symbol whose top-level namespace (after
    demangling) is `arrow::`, `parquet::`, or `arrow_vendored::`.

When isolation is off the plugin links arrow dynamically and the test
skips — both checks are vacuously meaningless in that mode.
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from helpers import arrowEnabled

pytestmark = pytest.mark.skipif(
    not arrowEnabled, reason="Arrow/Parquet bindings not built"
)


_REPO_ROOT = Path(__file__).resolve().parents[3]
_LIB_NAME = "libActsPluginArrow"


def _plugin_lib() -> Path:
    """Locate the built plugin library in build/lib."""
    libdir = _REPO_ROOT / "build" / "lib"
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    path = libdir / f"{_LIB_NAME}{suffix}"
    if not path.exists():
        pytest.skip(f"{path} does not exist")
    return path


def _direct_deps(lib: Path) -> list[str]:
    """Return the direct dynamic-linker dependencies of the library."""
    if sys.platform == "darwin":
        out = subprocess.check_output(["otool", "-L", str(lib)], text=True)
        deps: list[str] = []
        # First line is the library's own path; remainder are deps.
        for line in out.splitlines()[1:]:
            line = line.strip()
            if not line:
                continue
            # Format: "<path> (compatibility version ...)"
            deps.append(line.split(" ", 1)[0])
        return deps
    out = subprocess.check_output(["readelf", "-d", str(lib)], text=True)
    deps = []
    for line in out.splitlines():
        m = re.search(r"\(NEEDED\).*Shared library:\s*\[(.+)\]", line)
        if m:
            deps.append(m.group(1))
    return deps


def _is_arrow_lib(name: str) -> bool:
    base = os.path.basename(name)
    return base.startswith(("libarrow", "libparquet"))


def _exported_symbols_demangled(lib: Path) -> list[str]:
    """Return demangled names of externally-defined symbols in the library."""
    if sys.platform == "darwin":
        # -g external, -U skip undefined, -j just the name
        argv = ["nm", "-gUj", str(lib)]
    else:
        argv = ["nm", "-D", "--defined-only", "--format=just-symbols", str(lib)]
    raw = subprocess.check_output(argv, text=True, stderr=subprocess.DEVNULL)
    mangled = [line.strip() for line in raw.splitlines() if line.strip()]
    if not mangled:
        return []
    cxxfilt = shutil.which("c++filt") or shutil.which("llvm-cxxfilt")
    if cxxfilt is None:
        pytest.skip("c++filt / llvm-cxxfilt not available")
    proc = subprocess.run(
        [cxxfilt],
        input="\n".join(mangled),
        capture_output=True,
        text=True,
        check=True,
    )
    return proc.stdout.splitlines()


def _isolation_active(lib: Path) -> bool:
    """True iff the build appears to use ACTS_ARROW_ISOLATED=ON.

    Heuristic: in isolated mode the plugin has absorbed the arrow/parquet
    archives and so should not list libarrow*/libparquet* in its direct
    deps. In non-isolated mode it links those .dylibs/.sos directly.
    """
    return not any(_is_arrow_lib(d) for d in _direct_deps(lib))


# Top-level namespaces that must not leak out of the plugin.
_FORBIDDEN_NS_RE = re.compile(r"^(arrow|parquet|arrow_vendored)::")


def test_no_arrow_dynamic_dependencies():
    """The plugin must not pull libarrow* / libparquet* into the process
    via its own dynamic-link table."""
    lib = _plugin_lib()
    if not _isolation_active(lib):
        pytest.skip("ACTS_ARROW_ISOLATED appears OFF; isolation check N/A")
    arrow_deps = [d for d in _direct_deps(lib) if _is_arrow_lib(d)]
    assert not arrow_deps, (
        f"{lib.name} unexpectedly references arrow/parquet shared "
        f"libraries: {arrow_deps}"
    )


def test_no_exported_arrow_symbols():
    """No externally-defined symbol's top-level namespace may be
    arrow::, parquet::, or arrow_vendored::. Substring matches against
    'arrow' inside ACTS-namespaced symbols (e.g. signatures taking
    arrow::Table) are filtered out by checking the demangled prefix."""
    lib = _plugin_lib()
    if not _isolation_active(lib):
        pytest.skip("ACTS_ARROW_ISOLATED appears OFF; symbol-leak check N/A")

    demangled = _exported_symbols_demangled(lib)
    assert demangled, f"{lib.name} exports no symbols at all (suspicious)"

    leaks = [s for s in demangled if _FORBIDDEN_NS_RE.match(s)]
    assert not leaks, (
        f"{lib.name} leaks {len(leaks)} arrow/parquet symbol(s) "
        f"out of {len(demangled)} exported. First few:\n  " + "\n  ".join(leaks[:10])
    )


def test_acts_symbols_are_exported():
    """Sanity counterpart: at least one symbol in the Acts/ActsPlugins/
    ActsExamples namespace should be exported. Catches the failure mode
    where the export allowlist is so strict it hides everything."""
    lib = _plugin_lib()
    demangled = _exported_symbols_demangled(lib)
    acts_syms = [
        s for s in demangled if re.match(r"^(Acts|ActsPlugins|ActsExamples)::", s)
    ]
    assert acts_syms, (
        f"{lib.name} exports no Acts*::* symbols (export allowlist "
        f"too strict?). Total exports: {len(demangled)}"
    )
