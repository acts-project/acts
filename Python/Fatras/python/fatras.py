import warnings

try:
    from .ActsFatrasPythonBindings import *

    from acts._adapter import _patch_config
    from acts import ActsFatrasPythonBindings

    _patch_config(ActsFatrasPythonBindings)

except ImportError:
    warnings.warn("Fatras python bindings not available.")
    raise
