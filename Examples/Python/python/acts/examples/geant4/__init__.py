import subprocess
import sys

# Cannot conveniently catch linker errors, so we launch a subprocess to
# try importing and see if it works in order to provide a useful error message
try:
    subprocess.check_call(
        [sys.executable, "-c", "from acts import ActsExamplesPythonBindingsGeant4"]
    )
except subprocess.CalledProcessError as e:
    print("Error encountered importing DD4hep. Likely you need to set LD_LIBRARY_PATH.")
    sys.exit(1)


from acts._adapter import _patch_config
from acts import ActsExamplesPythonBindingsGeant4

_patch_config(ActsExamplesPythonBindingsGeant4)

from acts.ActsExamplesPythonBindingsGeant4 import *
