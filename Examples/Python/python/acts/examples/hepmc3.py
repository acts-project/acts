from acts._adapter import _patch_config
from acts import ActsPythonBindings

if not hasattr(ActsPythonBindings._examples, "_hepmc3"):
    raise ImportError("ActsPythonBindings._examples._hepmc3 not found")

_patch_config(ActsPythonBindings._examples._hepmc3)

from acts.ActsPythonBindings._examples._hepmc3 import *
