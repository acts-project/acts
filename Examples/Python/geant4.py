import subprocess
import sys

# Cannot conveniently catch linker errors, so we launch a subprocess to
# try importing and see if it works in order to provide a useful error message
try:
    subprocess.check_call(
        [
            sys.executable,
            "-c",
            "from acts.plugins.ActsPluginsPythonBindingsGeant4 import *",
        ]
    )
except subprocess.CalledProcessError as e:
    print(
        "Error encountered importing Geant4 as a plugin. Did you build with ACTS_BUILD_PLUGIN_GEANT4=ON ?"
    )
    sys.exit(1)

# from acts._adapter import _patch_config
# from acts import ActsPluginsPythonBindingsGeant4

# _patch_config(ActsPluginsPythonBindingsGeant4)

from .ActsPluginsPythonBindingsGeant4 import *
