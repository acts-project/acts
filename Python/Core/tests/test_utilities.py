import acts
import pytest


def test_logging():
    for l in ("VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL"):
        assert hasattr(acts.logging, l)
        assert hasattr(acts.logging.Level, l)
