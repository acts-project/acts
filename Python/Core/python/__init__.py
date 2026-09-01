from pathlib import Path
from typing import Callable, TypeVar, Union
import os
import warnings
import functools


from .ActsPythonBindings import *
from .ActsPythonBindings import __version__
from .ActsPythonBindings import _demo_histogram1, _demo_profile1, _demo_efficiency1
from .ActsPythonBindings import _testing
from . import ActsPythonBindings
from ._adapter import _patch_config
from .histogram import _patch_histogram_types

# Applied on import, before any logger exists: a logger copies the threshold
# when it is built.
if (_threshold := os.environ.get("ACTS_LOG_FAILURE_THRESHOLD")) is not None:
    _level = getattr(logging, _threshold, None)
    if isinstance(_level, logging.Level):
        logging.setFailureThreshold(_level)
    else:
        _names = ", ".join(lvl.name for lvl in logging.Level.__members__.values())
        _error = (
            f"`ACTS_LOG_FAILURE_THRESHOLD={_threshold}` is not a log level. "
            f"Expected one of: {_names}."
        )
        if "PYTEST_CURRENT_TEST" in os.environ:
            # test environment, fail hard
            raise RuntimeError(_error)
        else:
            warnings.warn(_error + "\nNo failure threshold will be applied.")


def Propagator(stepper, navigator, level=ActsPythonBindings.logging.INFO):
    for prefix in ("Eigen", "Atlas", "StraightLine"):
        _stepper = getattr(ActsPythonBindings, f"{prefix}Stepper")
        if isinstance(stepper, _stepper):
            return getattr(ActsPythonBindings, f"{prefix}Propagator")(
                stepper, navigator, level
            )
    raise TypeError(f"Unknown stepper {type(stepper).__name__}")


_patch_config(ActsPythonBindings)
_patch_histogram_types(ActsPythonBindings)


T = TypeVar("T")


class with_log_threshold:
    def __init__(self, level: logging.Level):
        self.level = level

    def __call__(self, func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args, **kwargs) -> T:
            with logging.ScopedFailureThreshold(self.level):
                return func(*args, **kwargs)

        return wrapper


@staticmethod
def _decoratorFromFile(file: Union[str, Path], **kwargs):
    if isinstance(file, str):
        file = Path(file)

    kwargs.setdefault("level", ActsPythonBindings.logging.INFO)

    if file.suffix in (".json", ".cbor"):
        from .json import MaterialMapJsonConverter, JsonMaterialDecorator

        c = MaterialMapJsonConverter.Config()
        for k in kwargs.keys():
            if hasattr(c, k):
                setattr(c, k, kwargs.pop(k))

        return JsonMaterialDecorator(jFileName=str(file), rConfig=c, **kwargs)
    elif file.suffix == ".root":
        from .root import RootMaterialDecorator

        return RootMaterialDecorator(fileName=str(file), **kwargs)
    else:
        raise ValueError(f"Unknown file type {file.suffix}")


ActsPythonBindings.IMaterialDecorator.fromFile = _decoratorFromFile
