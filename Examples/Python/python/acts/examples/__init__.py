import sys, inspect
from pathlib import Path
from typing import Optional, Protocol, Union, List, Dict
import os
import re

from acts.ActsPythonBindings._examples import *
from acts import ActsPythonBindings
import acts
from acts._adapter import _patch_config, _patch_detectors, _patchKwargsConstructor

_propagators = []
_concrete_propagators = []
for prefix in ("Eigen", "Atlas", "StraightLine"):
    _propagators.append(getattr(ActsPythonBindings._propagator, f"{prefix}Propagator"))
    _concrete_propagators.append(
        getattr(ActsPythonBindings._propagator, f"{prefix}ConcretePropagator")
    )


def ConcretePropagator(propagator):
    for prop, prop_if in zip(_propagators, _concrete_propagators):
        if isinstance(propagator, prop):
            return prop_if(propagator)

    raise TypeError(f"Unknown propagator {type(propagator).__name__}")


_patch_config(ActsPythonBindings._examples)

_patch_detectors(ActsPythonBindings._examples)

# Manually patch ExaTrkX constructors
# Need to do it this way, since they are not always present
for module in [
    "TorchMetricLearning",
    "OnnxMetricLearning",
    "TorchEdgeClassifier",
    "OnnxEdgeClassifier",
]:
    if hasattr(ActsPythonBindings._examples, module):
        _patchKwargsConstructor(getattr(ActsPythonBindings._examples, module))


def _makeLayerTriplet(*args, **kwargs):

    if len(args) == 1:
        _type = type(args[0])
        negative = central = positive = args[0]
    else:
        negative = kwargs.get("negative")
        central = kwargs.get("central")
        positive = kwargs.get("positive")
        types = []
        if negative is not None:
            types.append(type(negative))
        if central is not None:
            types.append(type(central))
        if positive is not None:
            types.append(type(positive))

        assert all(types[0] == t for t in types), "Inconsistent types"

        _type = types[0]

    def fill(obj):
        if negative is not None:
            obj.negative = negative
        if central is not None:
            obj.central = central
        if positive is not None:
            obj.positive = positive
        return obj

    if _type == bool:
        return fill(TGeoDetector.Config.LayerTripletBool())
    elif _type == list:
        if all(
            all(isinstance(v, str) for v in vv)
            for vv in (negative, central, positive)
            if vv is not None
        ):
            return fill(TGeoDetector.Config.LayerTripletVectorString())
        elif all(
            all(
                (isinstance(v, tuple) or isinstance(v, list))
                and isinstance(v[0], int)
                and isinstance(v[1], inspect.unwrap(TGeoDetector.Config.BinningType))
                for v in vv
            )
            for vv in (negative, central, positive)
            if vv is not None
        ):
            return fill(TGeoDetector.Config.LayerTripletVectorBinning())
        else:
            raise TypeError("Invalid types for list input")
    elif _type == tuple:
        if all(
            all(isinstance(v, float) for v in vv)
            for vv in (negative, central, positive)
            if vv is not None
        ):
            negative = Interval(*negative) if negative is not None else None
            central = Interval(*central) if central is not None else None
            positive = Interval(*positive) if positive is not None else None
            return fill(TGeoDetector.Config.LayerTripletInterval())
        else:
            raise TypeError("Invalid types for tuple input")
    elif _type == Interval:
        return fill(TGeoDetector.Config.LayerTripletInterval())
    elif _type == str:
        return fill(TGeoDetector.Config.LayerTripletString())
    elif _type == float:
        return fill(TGeoDetector.Config.LayerTripletDouble())
    else:
        raise TypeError("Unknown type given")


TGeoDetector.Config.LayerTriplet = _makeLayerTriplet


def _process_volume_intervals(kwargs):
    if len(kwargs) == 0:
        return kwargs  # prevent infinite recursion
    _kwargs = kwargs.copy()

    v = TGeoDetector.Config.Volume()
    for name, value in inspect.getmembers(inspect.unwrap(TGeoDetector.Config.Volume)):
        if not isinstance(getattr(v, name), inspect.unwrap(Interval)):
            continue
        if not name in _kwargs:
            continue
        if not isinstance(_kwargs[name], tuple):
            continue
        _kwargs[name] = Interval(*_kwargs[name])

    return _kwargs


_patchKwargsConstructor(TGeoDetector.Config.Volume, proc=_process_volume_intervals)


def NamedTypeArgs(**namedTypeArgs):
    """Decorator to move args of a named type (eg. `namedtuple` or `Enum`) to kwargs based on type, so user doesn't need to specify the key name.
    Also allows the keyword argument to be converted from a built-in type (eg. `tuple` or `int`)."""

    namedTypeClasses = {c: a for a, c in namedTypeArgs.items()}

    def NamedTypeArgsDecorator(func):
        from functools import wraps

        @wraps(func)
        def NamedTypeArgsWrapper(*args, **kwargs):
            from collections.abc import Iterable

            for k, v in kwargs.items():
                cls = namedTypeArgs.get(k)
                if (
                    cls is not None
                    and v is not None
                    and type(v).__module__ == int.__module__  # is v a 'builtins'?
                ):
                    if issubclass(cls, Iterable):
                        kwargs[k] = cls(*v)
                    else:
                        kwargs[k] = cls(v)

            newargs = []
            for i, a in enumerate(args):
                k = namedTypeClasses.get(type(a))
                if k is None:
                    newargs.append(a)
                    if i > len(newargs):
                        types = [type(a).__name__ for a in args]
                        raise TypeError(
                            f"{func.__name__}() positional argument {i} of type {type(a)} follows named-type arguments, which were converted to keyword arguments. All argument types: {types}"
                        )
                elif k in kwargs:
                    raise TypeError(f"{func.__name__}() keyword argument repeated: {k}")
                else:
                    kwargs[k] = a
            return func(*newargs, **kwargs)

        return NamedTypeArgsWrapper

    return NamedTypeArgsDecorator


def defaultKWArgs(**kwargs) -> dict:
    """Removes keyword arguments that are None or a list of all None (eg. [None,None]).
    This keeps the called function's defaults."""
    from collections.abc import Iterable

    return {
        k: v
        for k, v in kwargs.items()
        if not (
            v is None or (isinstance(v, Iterable) and all([vv is None for vv in v]))
        )
    }


def dump_args(func):
    """
    Decorator to print function call details.
    This includes parameters names and effective values.
    https://stackoverflow.com/questions/6200270/decorator-that-prints-function-call-details-parameters-names-and-effective-valu
    """
    from functools import wraps

    @wraps(func)
    def dump_args_wrapper(*args, **kwargs):
        try:
            func_args = inspect.signature(func).bind(*args, **kwargs).arguments
            func_args_str = ", ".join(
                map("{0[0]} = {0[1]!r}".format, func_args.items())
            )
        except ValueError:
            func_args_str = ", ".join(
                list(map("{0!r}".format, args))
                + list(map("{0[0]} = {0[1]!r}".format, kwargs.items()))
            )
        if not (
            func_args_str == ""
            and any([a == "Config" for a in func.__qualname__.split(".")])
        ):
            print(f"{func.__module__}.{func.__qualname__} ( {func_args_str} )")
        return func(*args, **kwargs)

    # fix up any attributes broken by the wrapping
    for name in dir(func):
        if not name.startswith("__"):
            obj = getattr(func, name)
            wrapped = getattr(dump_args_wrapper, name, None)
            if type(obj) is not type(wrapped):
                setattr(dump_args_wrapper, name, obj)

    return dump_args_wrapper


def dump_args_calls(myLocal=None, mods=None, quiet=False):
    """
    Wrap all Python bindings calls to acts and its submodules in dump_args.
    Specify myLocal=locals() to include imported symbols too.
    """
    import collections

    def _allmods(mod, base, found):
        import types

        mods = [mod]
        found.add(mod)
        for name, obj in sorted(
            vars(mod).items(),
            key=lambda m: (2, m[0])
            if m[0] == "ActsPythonBindings"
            else (1, m[0])
            if m[0].startswith("_")
            else (0, m[0]),
        ):
            if (
                not name.startswith("__")
                and type(obj) is types.ModuleType
                and obj.__name__.startswith(base)
                and f"{mod.__name__}.{name}" in sys.modules
                and obj not in found
            ):
                mods += _allmods(obj, base, found)
        return mods

    if mods is None:
        mods = _allmods(acts, "acts.", set())
    elif not isinstance(mods, list):
        mods = [mods]

    donemods = []
    alldone = 0
    for mod in mods:
        done = 0
        for name in dir(mod):
            # if name in {"Config", "Interval", "IMaterialDecorator"}: continue  # skip here if we don't fix up attributes in dump_args and unwrap classes elsewhere
            obj = getattr(mod, name, None)
            if not (
                not name.startswith("__")
                and isinstance(obj, collections.abc.Callable)
                and hasattr(obj, "__module__")
                and obj.__module__.startswith("acts.ActsPythonBindings")
                and not hasattr(obj, "__wrapped__")
            ):
                continue
            # wrap class's contained methods
            done += dump_args_calls(myLocal, [obj], True)
            wrapped = dump_args(obj)
            setattr(mod, name, wrapped)
            if myLocal and hasattr(myLocal, name):
                setattr(myLocal, name, wrapped)
            done += 1
        if done:
            alldone += done
            donemods.append(f"{mod.__name__}:{done}")
    if not quiet and donemods:
        print("dump_args for module functions:", ", ".join(donemods))
    return alldone


class CustomLogLevel(Protocol):
    def __call__(
        self,
        minLevel: acts.logging.Level = acts.logging.VERBOSE,
        maxLevel: acts.logging.Level = acts.logging.FATAL,
    ) -> acts.logging.Level:
        ...


def defaultLogging(
    s=None,
    logLevel: Optional[acts.logging.Level] = None,
) -> CustomLogLevel:
    """
    Establishes a default logging strategy for the python examples interface.

    Returns a function that determines the log level in the following schema:
    - if `logLevel` is set use it otherwise use the log level of the sequencer `s.config.logLevel`
    - the returned log level is bound between `minLevel` and `maxLevel` provided to `customLogLevel`

    Examples:
    - `customLogLevel(minLevel=acts.logging.INFO)` to get a log level that is INFO or higher
      (depending on the sequencer and `logLevel` param) which is useful to suppress a component which
      produces a bunch of logs below INFO and you are actually more interested in another component
    - `customLogLevel(maxLevel=acts.logging.INFO)` to get a log level that is INFO or lower
      (depending on the sequencer and `logLevel` param) which is useful to get more details from a
      component that will produce logs of interest below the default level
    - in summary `minLevel` defines the maximum amount of logging and `maxLevel` defines the minimum amount of logging
    """

    def customLogLevel(
        minLevel: acts.logging.Level = acts.logging.VERBOSE,
        maxLevel: acts.logging.Level = acts.logging.FATAL,
    ) -> acts.logging.Level:
        l = logLevel if logLevel is not None else s.config.logLevel
        return acts.logging.Level(min(maxLevel.value, max(minLevel.value, l.value)))

    return customLogLevel


class Sequencer(ActsPythonBindings._examples._Sequencer):

    _autoFpeMasks: Optional[List["FpeMask"]] = None

    def __init__(self, *args, **kwargs):

        if "fpeMasks" in kwargs:
            m = kwargs["fpeMasks"]
            if isinstance(m, list) and len(m) > 0 and isinstance(m[0], tuple):
                n = []
                for loc, fpe, count in m:
                    t = _fpe_types_to_enum[fpe] if isinstance(fpe, str) else fpe
                    n.append(self.FpeMask(loc, t, count))
                kwargs["fpeMasks"] = n

            kwargs["fpeMasks"] = kwargs.get("fpeMasks", []) + self._getAutoFpeMasks()

        cfg = self.Config()
        if len(args) == 1 and isinstance(args[0], self.Config):
            cfg = args[0]
            args = args[1:]
        if "config" in kwargs:
            cfg = kwargs.pop("config")

        for k, v in kwargs.items():
            if not hasattr(cfg, k):
                raise ValueError(f"Sequencer.Config does not have field {k}")
            if isinstance(v, Path):
                v = str(v)

            setattr(cfg, k, v)

        super().__init__(cfg)

    class FpeMask(ActsPythonBindings._examples._Sequencer._FpeMask):
        @classmethod
        def fromFile(cls, file: Union[str, Path]) -> List["FpeMask"]:
            if isinstance(file, str):
                file = Path(file)

            if file.suffix in (".yml", ".yaml"):
                try:
                    return cls.fromYaml(file)
                except ImportError:
                    print("FPE mask input file is YAML, but PyYAML is not installed")
                    raise

        @classmethod
        def fromYaml(cls, file: Union[str, Path]) -> List["FpeMask"]:
            import yaml

            with file.open() as fh:
                d = yaml.safe_load(fh)

            return cls.fromDict(d)

        _fpe_types_to_enum = {v.name: v for v in acts.FpeType.values}

        @staticmethod
        def toDict(
            masks: List["FpeMask"],
        ) -> Dict[str, Dict[str, int]]:
            out = {}
            for mask in masks:
                out.setdefault(mask.loc, {})
                out[mask.loc][mask.type.name] = mask.count

                return out

        @classmethod
        def fromDict(cls, d: Dict[str, Dict[str, int]]) -> List["FpeMask"]:
            out = []
            for loc, types in d.items():
                for fpe, count in types.items():
                    out.append(cls(loc, cls._fpe_types_to_enum[fpe], count))
            return out

    @classmethod
    def _getAutoFpeMasks(cls) -> List[FpeMask]:
        if cls._autoFpeMasks is not None:
            return cls._autoFpeMasks

        srcdir = Path(cls._sourceLocation).parent.parent.parent.parent

        cls._autoFpeMasks = []

        for root, _, files in os.walk(srcdir):
            root = Path(root)
            for f in files:
                if (
                    not f.endswith(".hpp")
                    and not f.endswith(".cpp")
                    and not f.endswith(".ipp")
                ):
                    continue
                f = root / f
                #  print(f)
                with f.open("r") as fh:
                    for i, line in enumerate(fh, start=1):
                        if m := re.match(r".*\/\/ ?MARK: ?(fpeMask.*)$", line):
                            exp = m.group(1)
                            for m in re.findall(
                                r"fpeMask\((\w+), ?(\d+) ?, ?issue: ?(\d+)\)", exp
                            ):
                                fpeType, count, _ = m
                                count = int(count)
                                rel = f.relative_to(srcdir)
                                cls._autoFpeMasks.append(
                                    cls.FpeMask(
                                        f"{rel}:{i}",
                                        cls.FpeMask._fpe_types_to_enum[fpeType],
                                        count,
                                    )
                                )

        return cls._autoFpeMasks
