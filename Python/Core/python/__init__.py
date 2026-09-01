from pathlib import Path
from typing import Union


from .ActsPythonBindings import *
from .ActsPythonBindings import __version__
from .ActsPythonBindings import _demo_histogram1, _demo_profile1, _demo_efficiency1
from .ActsPythonBindings import _testing
from . import ActsPythonBindings
from ._adapter import _patch_config
from .histogram import _patch_histogram_types


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
