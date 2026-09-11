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
        "Log failure threshold could not be set: this build was configured with "
        "`ACTS_ENABLE_LOG_FAILURE_THRESHOLD=OFF`. "
        f"The pytest test-suite can produce false-{errtype} results in this configuration"
    )
