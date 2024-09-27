from acts._adapter import _patch_config
from acts import ActsPythonBindings

_patch_config(ActsPythonBindings._examples.traccc)

from acts.ActsPythonBindings._examples.traccc import *
