import sys, inspect

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
                and isinstance(v[1], TGeoDetector.Config.BinningType)
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
    for name, value in inspect.getmembers(TGeoDetector.Config.Volume):
        if not isinstance(getattr(v, name), Interval):
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
                        raise TypeError(
                            f"{func.__name__}() positional argument {i} follows named-type arguments, which were converted to keyword arguments"
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
        import inspect

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
        print(f"{func.__module__}.{func.__qualname__} ( {func_args_str} )")
        return func(*args, **kwargs)

    return dump_args_wrapper


def dump_args_calls(
    myLocal=None,
    mod=sys.modules[__name__],
):
    """
    Wrap all calls to acts.examples Python bindings in dump_args.
    Specify myLocal=locals() to include imported symbols too.
    """
    import collections

    for n in dir(mod):
        if n.startswith("_") or n == "Config" or n == "Interval":
            continue
        f = getattr(mod, n)
        if not (
            isinstance(f, collections.abc.Callable)
            and f.__module__.startswith("acts.ActsPythonBindings")
            and not hasattr(f, "__wrapped__")
        ):
            continue
        dump_args_calls(myLocal, f)  # wrap class's contained methods
        w = dump_args(f)
        setattr(mod, n, w)
        if myLocal and hasattr(myLocal, n):
            setattr(myLocal, n, w)
