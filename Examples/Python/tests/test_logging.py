import acts
import pytest


def test_logging_threshold():
    assert acts.logging.getFailureThreshold() == acts.logging.Level.WARNING


def test_logging_threshold_context_manager():
    with acts.logging.ScopedFailureThreshold(acts.logging.ERROR):
        assert acts.logging.getFailureThreshold() == acts.logging.Level.ERROR
    assert acts.logging.getFailureThreshold() == acts.logging.Level.WARNING


def test_logging_threshold_context_manager_exception():
    with pytest.raises(RuntimeError):
        with acts.logging.ScopedFailureThreshold(level=acts.logging.ERROR):
            assert acts.logging.getFailureThreshold() == acts.logging.Level.ERROR
            raise RuntimeError("test")
    assert acts.logging.getFailureThreshold() == acts.logging.Level.WARNING
