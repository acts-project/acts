import warnings
import acts

try:
    if acts.logging.getFailureThreshold() != acts.logging.WARNING:
        acts.logging.setFailureThreshold(acts.logging.WARNING)
except RuntimeError:
    # Repackage with different error string
    errtype = (
        "negative"
        if acts.logging.getFailureThreshold() < acts.logging.WARNING
        else "positive"
    )
    warnings.warn(
        "Runtime log failure threshold could not be set. "
        "Compile-time value is probably set via CMake, i.e. "
        f"`ACTS_LOG_FAILURE_THRESHOLD={acts.logging.getFailureThreshold().name}` is set, "
        "or `ACTS_ENABLE_LOG_FAILURE_THRESHOLD=OFF`. "
        f"The pytest test-suite can produce false-{errtype} results in this configuration"
    )
