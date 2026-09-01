import warnings
import acts

# setFailureThreshold is a no-op in a build configured with
# ACTS_ENABLE_LOG_FAILURE_THRESHOLD=OFF, so check that it took effect rather
# than relying on it to raise.
acts.logging.setFailureThreshold(acts.logging.WARNING)
if acts.logging.getFailureThreshold() != acts.logging.WARNING:
    errtype = (
        "negative"
        if acts.logging.getFailureThreshold() < acts.logging.WARNING
        else "positive"
    )
    warnings.warn(
        "Log failure threshold could not be set; it is "
        f"`{acts.logging.getFailureThreshold().name}`. This build was probably "
        "configured with `ACTS_ENABLE_LOG_FAILURE_THRESHOLD=OFF`. "
        f"The pytest test-suite can produce false-{errtype} results in this configuration"
    )
