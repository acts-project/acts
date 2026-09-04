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


def test_logging():
    for l in ("VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL"):
        assert hasattr(acts.logging, l)
        assert hasattr(acts.logging.Level, l)


def test_pgd_particle():
    assert len(acts.PdgParticle.__members__) == 32


def test_algebra():
    v3 = acts.Vector3(1, 2, 3)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2, 3, 4)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2)

    v3 = acts.Vector3([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2])
    with pytest.raises(TypeError):
        acts.Vector3()

    v4 = acts.Vector4(1, 2, 3, 4)
    with pytest.raises(TypeError):
        v4 = acts.Vector4(1, 2, 3)
    v4 = acts.Vector4([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector4([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector4()


def test_geometry_context_factory():
    """Test that GeometryContext factory method works without warnings"""
    import warnings

    # New factory method should not produce warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
        assert gctx is not None


def test_geometry_context_deprecated_constructor():
    """Test that GeometryContext default constructor produces deprecation warning"""
    # Old constructor should produce a deprecation warning
    with pytest.warns(DeprecationWarning, match="GeometryContext.*deprecated"):
        gctx = acts.GeometryContext()
        assert gctx is not None
