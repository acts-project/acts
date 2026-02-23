import acts
import pytest


def test_get_default_logger_creates_logger():
    logger = acts.getDefaultLogger("test_logger", acts.logging.INFO)
    assert logger is not None
    assert logger.name == "test_logger"
    assert logger.level == acts.logging.INFO


def test_get_default_logger_default_level():
    logger = acts.getDefaultLogger("test_default_level")
    assert logger.level == acts.logging.INFO


@acts.with_log_threshold(acts.logging.MAX)
def test_get_default_logger_different_levels():
    for level in (
        acts.logging.VERBOSE,
        acts.logging.DEBUG,
        acts.logging.WARNING,
        acts.logging.ERROR,
        acts.logging.FATAL,
    ):
        logger = acts.getDefaultLogger(f"logger_{level}", level)
        assert logger.level == level


@acts.with_log_threshold(acts.logging.MAX)
def test_get_default_logger_log_methods(capfd):
    logger = acts.getDefaultLogger("test_log_methods", acts.logging.VERBOSE)
    logger.verbose("verbose message")
    logger.debug("debug message")
    logger.info("info message")
    logger.warning("warning message")
    logger.error("error message")
    logger.fatal("fatal message")

    captured = capfd.readouterr()
    assert "test_log_met" in captured.out
    assert "verbose message" in captured.out
    assert "debug message" in captured.out
    assert "info message" in captured.out
    assert "warning message" in captured.out
    assert "error message" in captured.out
    assert "fatal message" in captured.out


def test_get_default_logger_format_args(capfd):
    logger = acts.getDefaultLogger("test_format", acts.logging.INFO)
    logger.info("value: {}", 42)
    logger.info("two values: {} {}", "foo", "bar")

    captured = capfd.readouterr()
    assert "test_format" in captured.out
    assert "value: 42" in captured.out
    assert "two values: foo bar" in captured.out


def test_get_default_logger_level_filtering(capfd):
    logger = acts.getDefaultLogger("test_filter", acts.logging.INFO)
    logger.debug("filtered out")
    logger.info("should appear")

    captured = capfd.readouterr()
    assert "should appear" in captured.out
    assert "filtered out" not in captured.out


def test_logging_threshold():
    assert acts.logging.getFailureThreshold() == acts.logging.WARNING


def test_logging_threshold_context_manager():
    with acts.logging.ScopedFailureThreshold(acts.logging.ERROR):
        assert acts.logging.getFailureThreshold() == acts.logging.ERROR
    assert acts.logging.getFailureThreshold() == acts.logging.WARNING


def test_logging_threshold_context_manager_exception():
    with pytest.raises(RuntimeError):
        with acts.logging.ScopedFailureThreshold(level=acts.logging.ERROR):
            assert acts.logging.getFailureThreshold() == acts.logging.ERROR
            raise RuntimeError("test")
    assert acts.logging.getFailureThreshold() == acts.logging.WARNING
