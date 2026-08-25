"""Runtime paths and environment setup (before ``import acts`` where needed)."""

from __future__ import annotations

import ctypes
import os
import sys
from pathlib import Path
from typing import Optional

_DIGI_CONFIG_REL = Path("Examples/Configs/odd-digi-smearing-config.json")


def preload_dd4hep_for_macos() -> None:
    """Preload DD4hep plugins with RTLD_GLOBAL on macOS before importing acts."""
    if sys.platform != "darwin":
        return
    root = os.environ.get("DD4HEP_ROOT") or os.environ.get("DD4HEP_INSTALL")
    libdir: Optional[Path] = Path(root) / "lib" if root else None
    if libdir is None or not libdir.is_dir():
        try:
            candidate = (
                Path(__file__).resolve().parents[5] / "DD4hep" / "install" / "lib"
            )
        except (IndexError, OSError):
            candidate = None
        if candidate is not None and candidate.is_dir():
            libdir = candidate
    if libdir is None or not libdir.is_dir():
        return
    path_bits = [str(libdir)]
    odd = os.environ.get("ODD_PATH")
    if odd:
        fdir = Path(odd).resolve().parent / "odd-build" / "factory"
        if fdir.is_dir():
            path_bits.append(str(fdir))
    for var in ("DYLD_LIBRARY_PATH", "DD4HEP_LIBRARY_PATH"):
        cur = os.environ.get(var, "")
        merged = []
        for p in path_bits + [x for x in cur.split(":") if x]:
            if p not in merged:
                merged.append(p)
        os.environ[var] = ":".join(merged)
    mode = getattr(ctypes, "RTLD_GLOBAL", 8)
    for names in (
        ("libDD4hepGaudiPluginMgr.1.32.dylib", "libDD4hepGaudiPluginMgr.dylib"),
        ("libDDCorePlugins.1.32.dylib", "libDDCorePlugins.dylib"),
    ):
        for n in names:
            p = libdir / n
            if p.is_file():
                try:
                    ctypes.CDLL(str(p), mode=mode)
                except OSError:
                    break
                break


def setup_runtime_env() -> None:
    """ACTS Examples runtime flags shared by all workflow stages."""
    if os.environ.get("ACTS_SEQUENCER_FAIL_ON_UNMASKED_FPE") is None:
        os.environ["ACTS_SEQUENCER_FAIL_ON_UNMASKED_FPE"] = "0"


def find_acts_source_dir(start: Optional[Path] = None) -> Path:
    """Walk upward from ``start`` (default: this file) to find the ACTS source root."""
    here = (start or Path(__file__)).resolve()
    if here.is_file():
        here = here.parent
    for candidate in (here, *here.parents):
        if (candidate / _DIGI_CONFIG_REL).is_file():
            return candidate
    raise FileNotFoundError(
        f"Could not find {_DIGI_CONFIG_REL} by walking up from {here}"
    )


def default_digi_config(srcdir: Path) -> Path:
    return srcdir / _DIGI_CONFIG_REL
