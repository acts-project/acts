from acts._adapter import _patch_config
from acts import ActsPythonBindings

if not hasattr(ActsPythonBindings._examples, "_edm4hep"):
    raise ImportError("ActsPythonBindings._examples._edm4hep not found")

_patch_config(ActsPythonBindings._examples._edm4hep)

from acts.ActsPythonBindings._examples._edm4hep import *
