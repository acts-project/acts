import inspect
import functools
from typing import Optional, Callable, Dict, Any

import acts


def _make_config_adapter(fn):
    @functools.wraps(fn)
    def wrapped(self, *args, **kwargs):
        if len(args) > 0:
            maybe_config = args[0]
            if isinstance(maybe_config, inspect.unwrap(type(self).Config)):
                # is already config, nothing to do here
                fn(self, maybe_config, *args[1:], **kwargs)
                return

        if "config" in kwargs:
            config = kwargs.pop("config")
            fn(self, config, *args, **kwargs)
            return

        cfg = type(self).Config()
        _kwargs = {}
        for k, v in kwargs.items():
            if hasattr(cfg, k):
                setattr(cfg, k, v)
            else:
                _kwargs[k] = v
        try:
            fn(self, cfg, *args, **_kwargs)
        except TypeError as e:
            import textwrap

            print("-" * 80)
            print("Patched config constructor failed for", type(self))
            message = (
                "This is most likely because one of the following kwargs "
                "could not be assigned to the Config object, and the constructor "
                "did not accept it as an additional argument:"
            )
            print("\n".join(textwrap.wrap(message, width=80)))
            print("->", ", ".join(_kwargs.keys()))
            members = inspect.getmembers(type(cfg), lambda a: not inspect.isroutine(a))
            members = [m for m, _ in members if not m.startswith("_")]
            print(type(cfg), "has the following properties:\n->", ", ".join(members))
            print("-" * 80)
            raise e

    return wrapped


def _make_config_constructor(
    cls, proc: Optional[Callable[[Dict[str, Any]], Dict[str, Any]]] = None
):
    fn = cls.__init__

    @functools.wraps(fn)
    def wrapped(self, *args, **kwargs):
        _kwargs = {}
        for k in list(kwargs.keys()):
            if hasattr(cls, k):
                _kwargs[k] = kwargs.pop(k)

        fn(self, *args, **kwargs)

        if proc is not None:
            _kwargs = proc(_kwargs)
        for k, v in _kwargs.items():
            setattr(self, k, v)

    return wrapped


def _patchKwargsConstructor(
    cls, proc: Optional[Callable[[Dict[str, Any]], Dict[str, Any]]] = None
):
    cls.__init__ = _make_config_constructor(cls, proc)


def _patch_config(m):
    for name, cls in inspect.getmembers(m, inspect.isclass):
        if name == "Config":
            _patchKwargsConstructor(cls)

        if name.endswith("Detector"):
            continue

        if hasattr(cls, "Config"):
            cls.__init__ = _make_config_adapter(cls.__init__)
            _patchKwargsConstructor(cls.Config)


def _detector_create(cls, config_class=None):
    def create(*args, mdecorator=None, **kwargs):
        if mdecorator is not None:
            if not isinstance(mdecorator, inspect.unwrap(acts.IMaterialDecorator)):
                raise TypeError("Material decorator is not valid")
        if config_class is None:
            cfg = cls.Config()
        else:
            cfg = config_class()
        _kwargs = {}
        for k, v in kwargs.items():
            try:
                setattr(cfg, k, v)
            except AttributeError:
                _kwargs[k] = v
        det = cls()
        tg, deco = det.finalize(cfg, mdecorator, *args, **_kwargs)
        return det, tg, deco

    return create


def _patch_detectors(m):
    for name, cls in inspect.getmembers(m, inspect.isclass):
        if name.endswith("Detector"):
            cls.create = _detector_create(cls)
