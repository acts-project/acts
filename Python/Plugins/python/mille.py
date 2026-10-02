import warnings

try:
    from .ActsPluginsPythonBindingsMille import *

    from acts._adapter import _patch_config, _patchKwargsConstructor
    from acts import ActsPluginsPythonBindingsMille

    _patch_config(
        ActsPluginsPythonBindingsMille,
        [
            "Config",
            "MillePedeSteeringConfig",
            "MillePedeEqualityConstraint",
            "MillePedeParameterResult",
            "Result",
        ],
    )

    # manually add kwargs constructors for structs not matching


except ModuleNotFoundError:
    warnings.warn(
        "Mille plugin not available. Try building with @PLUGIN_BUILD_FLAG@=ON."
    )
    raise
