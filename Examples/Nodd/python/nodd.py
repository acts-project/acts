"""Helpers appended to the generated acts.examples.nodd binding module."""

import os
import sys
from pathlib import Path

import acts
import acts.examples


def getNoddDirectory(nodd_dir=None):
    """Resolve an explicit nODD checkout, or the NODD_PATH environment variable."""
    if nodd_dir is None:
        nodd_dir = os.environ.get("NODD_PATH")
    if not nodd_dir:
        raise RuntimeError("Set NODD_PATH or pass nodd_dir to the nODD checkout")
    directory = Path(nodd_dir).expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"nODD directory not found: {directory}")
    return directory


def getNoddDetector(
    nodd_dir=None,
    build_dir=None,
    compact_file=None,
    logLevel=acts.logging.INFO,
    gen3=True,
):
    """Load the BeamPipe + complete Tracker compact and derive Gen3 geometry.

    build_dir defaults to <nodd_dir>/build/dd4hep; a relative override is relative
    to nodd_dir. An explicit compact_file is relative to the current directory.
    The model must already be exported and built by the nODD repository.
    Sensor transforms and pipe material come from DD4hep. No ODD map is used.
    Pass gen3=False for DD4hep-only Geant4 jobs; explicit pixel controls remain available.
    """
    directory = getNoddDirectory(nodd_dir)
    build = Path(build_dir) if build_dir is not None else Path("build/dd4hep")
    if not build.is_absolute():
        build = directory / build
    plugin_dir = build.expanduser().resolve() / "detector"
    compact = (
        Path(compact_file).expanduser().resolve()
        if compact_file is not None
        else plugin_dir / "tracker/tracker.xml"
    )
    if not compact.is_file():
        raise FileNotFoundError(
            f"nODD tracker compact not found: {compact}. Build the nODD DD4hep model first."
        )

    suffix = {"darwin": ".dylib", "linux": ".so"}.get(sys.platform)
    if suffix is None:
        raise RuntimeError(f"nODD factory loading is unsupported on {sys.platform}")
    libraries = ["nODDPixelBarrel"]
    if compact_file is None or gen3:
        libraries += [
            "nODDBeamPipe",
            "nODDShortStripBarrel",
            "nODDShortStripEndcap",
            "nODDLongStrip",
        ]
    import ROOT

    for name in libraries:
        library = plugin_dir / f"lib{name}{suffix}"
        components = plugin_dir / f"lib{name}.components"
        for path in (library, components):
            if not path.is_file():
                raise FileNotFoundError(f"nODD factory file not found: {path}")
        if ROOT.gSystem.Load(str(library)) < 0:
            raise RuntimeError(f"Failed to load nODD factory library: {library}")
    search_path = os.environ.get("DD4HEP_LIBRARY_PATH", "").split(os.pathsep)
    if str(plugin_dir) not in search_path:
        os.environ["DD4HEP_LIBRARY_PATH"] = os.pathsep.join(
            [str(plugin_dir), *filter(None, search_path)]
        )

    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)
    config = DD4hepNoddDetector.Config(
        xmlFileNames=[str(compact)],
        gen3=gen3,
        logLevel=customLogLevel(),
        dd4hepLogLevel=customLogLevel(minLevel=acts.logging.WARNING),
    )
    return DD4hepNoddDetector(config)
