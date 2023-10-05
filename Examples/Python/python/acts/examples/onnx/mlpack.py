from acts._adapter import _patch_config
from acts import ActsPythonBindings

if not hasattr(ActsPythonBindings._examples, "_mlpack"):
    raise ImportError("ActsPythonBindings._examples._mlpack not found")

_patch_config(ActsPythonBindings._examples._mlpack)

from acts.ActsPythonBindings._examples._mlpack import *
