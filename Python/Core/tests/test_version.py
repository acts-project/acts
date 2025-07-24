import pytest
import acts


def test_version():
    assert hasattr(acts, "__version__")
    assert hasattr(acts, "version")
    assert hasattr(acts.version, "major")
    assert hasattr(acts.version, "minor")
    assert hasattr(acts.version, "patch")
    assert hasattr(acts.version, "commit_hash")
    assert hasattr(acts.version, "commit_hash_short")
