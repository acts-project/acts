import acts
import acts.examples
import pytest


def test_logging_threshold():
    assert acts.examples.getLogFailureThreshold() == acts.logging.WARNING


def test_logging_threshold_context_manager():
    with acts.examples.ScopedFailureThreshold(acts.logging.ERROR):
        assert acts.examples.getLogFailureThreshold() == acts.logging.ERROR
    assert acts.examples.getLogFailureThreshold() == acts.logging.WARNING


def test_logging_threshold_context_manager_exception():
    with pytest.raises(RuntimeError):
        with acts.examples.ScopedFailureThreshold(level=acts.logging.ERROR):
            assert acts.examples.getLogFailureThreshold() == acts.logging.ERROR
            raise RuntimeError("test")
    assert acts.examples.getLogFailureThreshold() == acts.logging.WARNING


def test_default_logger_is_not_armed():
    """Core builds unarmed loggers; only the framework arms what it builds."""
    logger = acts.getDefaultLogger("unarmed", acts.logging.VERBOSE)
    assert logger.failureThreshold == acts.logging.MAX
