from acts._adapter import _patch_config
from acts import ActsPythonBindings

_patch_config(ActsPythonBindings._examples._edm4hep)

from acts.ActsPythonBindings._examples._edm4hep import *
