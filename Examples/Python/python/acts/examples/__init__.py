import sys, inspect
from pathlib import Path
from typing import Optional, Protocol, Union, List, Dict, Tuple
import os
import re

from acts.ActsPythonBindings._examples import *
from acts import ActsPythonBindings
import acts
from acts._adapter import _patch_config, _patchKwargsConstructor

_propagators = []
_concrete_propagators = []
for stepper in ("Eigen", "Atlas", "StraightLine", "Sympy"):
    _propagators.append(getattr(ActsPythonBindings._propagator, f"{stepper}Propagator"))
    _concrete_propagators.append(
        getattr(
            ActsPythonBindings._propagator,
            f"{stepper}ConcretePropagator",
        )
    )


def ConcretePropagator(propagator):
    for prop, prop_if in zip(_propagators, _concrete_propagators):
        if isinstance(propagator, prop):
            return prop_if(propagator)

    raise TypeError(f"Unknown propagator {type(propagator).__name__}")


_patch_config(ActsPythonBindings._examples)


