from pathlib import Path
import sys
import warnings
import pytest_check as check


sys.path += [
    str(Path(__file__).parent / "Examples/Scripts/Python"),
    str(Path(__file__).parent / "Examples/Python/tests" ),
]


import pytest

import acts

try:
    import ROOT

    ROOT.gSystem.ResetSignals()
except ImportError:
    pass

try:
    if acts.logging.getFailureThreshold() != acts.logging.WARNING:
        acts.logging.setFailureThreshold(acts.logging.WARNING)
except RuntimeError:
    # Repackage with different error string
    errtype = (
        "negative"
        if acts.logging.getFailureThreshold() < acts.logging.WARNING
        else "positive"
    )
    warnings.warn(
        "Runtime log failure threshold could not be set. "
        "Compile-time value is probably set via CMake, i.e. "
        f"`ACTS_LOG_FAILURE_THRESHOLD={acts.logging.getFailureThreshold().name}` is set, "
        "or `ACTS_ENABLE_LOG_FAILURE_THRESHOLD=OFF`. "
        f"The pytest test-suite can produce false-{errtype} results in this configuration"
    )


def pytest_addoption(parser):
    parser.addoption("--physmon-output-path", action="store", default=Path.cwd()/"physmon", type=Path)
