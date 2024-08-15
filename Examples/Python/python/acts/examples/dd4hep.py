import subprocess
import sys


# Cannot conveniently catch linker errors, so we launch a suprocess to
# try importing and see if it works in order to provide a useful error message
try:
    subprocess.check_call(
        [sys.executable, "-c", "from acts import ActsPythonBindingsDD4hep"]
    )
except subprocess.CalledProcessError as e:
    print("Error encountered importing DD4hep. Likely you need to set LD_LIBRARY_PATH.")
    sys.exit(1)

from acts._adapter import _patch_config, _detector_create, _patch_detectors
from acts import ActsPythonBindingsDD4hep

_patch_config(ActsPythonBindingsDD4hep)
_patch_detectors(ActsPythonBindingsDD4hep)
ActsPythonBindingsDD4hep.DD4hepDetector.create = _detector_create(
    ActsPythonBindingsDD4hep.DD4hepDetector,
    ActsPythonBindingsDD4hep.DD4hepGeometryService.Config,
)

from acts.ActsPythonBindingsDD4hep import *
