from . import ActsPluginsPythonBindings
from acts.plugins.ActsPluginsPythonBindings import *

# Cannot conveniently catch linker errors, so we launch a suprocess to
# try importing and see if it works in order to provide a useful error message
try:
    from . import ActsPluginsPythonBindingsDD4hep as dd4hep
except ImportError as e:
    print("Error encountered importing DD4hep as a plugin.")
    print("Either you did not build with ACTS_BUILD_PLUGIN_DD4HEP=ON,")
    print("or you need to set (DY)LD_LIBRARY_PATH to include the DD4hep librariess.")

# Cannot conveniently catch linker errors, so we launch a suprocess to
# try importing and see if it works in order to provide a useful error message
try:
    from . import ActsPluginsPythonBindingsGeant4 as geant4
except ImportError as e:
    print("Error encountered importing Geant4 as a plugin.")
    print("DID you build with ACTS_BUILD_PLUGIN_GEANT=ON ?")
