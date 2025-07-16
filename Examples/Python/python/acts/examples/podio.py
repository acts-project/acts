from acts._adapter import _patch_config
from acts import ActsPythonBindingsPodio

_patch_config(ActsPythonBindingsPodio)

from acts.ActsPythonBindingsPodio import *
