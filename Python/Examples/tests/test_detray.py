import acts
import pytest

from helpers import detrayEnabled

u = acts.UnitConstants


@pytest.mark.skipif(not detrayEnabled, reason="Detray plugin not available")
def test_detray_propagator_bindings():
    import acts.examples.detray as exdetray

    assert hasattr(exdetray, "StraightLinePropagatorODD")
    assert hasattr(exdetray, "HelixPropagatorODD")


@pytest.mark.skipif(not detrayEnabled, reason="Detray plugin not available")
def test_detray_surface_grid_policy():
    from acts import detray
    import acts.json

    cfg = detray.DetraySurfaceGridConfig()
    assert cfg.vertexZHalfRange == pytest.approx(150 * u.mm)
    assert cfg.minPt == pytest.approx(0.5 * u.GeV)
    assert cfg.numberOfBinsFactor == pytest.approx(2.0)
    assert cfg.lookupOnLayer
    cfg.minPt = 1 * u.GeV
    cfg.numberOfBinsFactor = None
    assert cfg.numberOfBinsFactor is None

    factory = detray.addDetraySurfaceGridPolicy(
        acts.NavigationPolicyFactory.make().add(acts.TryAllNavigationPolicy), cfg
    )
    assert factory is not None

    payloadCfg = detray.DetrayPayloadConverter.Config()
    assert payloadCfg.surfaceGridFill is not None
    payloadCfg.surfaceGridFill = None
    assert payloadCfg.surfaceGridFill is None
    payloadCfg.surfaceGridFill = cfg
    assert payloadCfg.surfaceGridFill.minPt == pytest.approx(1 * u.GeV)

    jsonCfg = acts.json.TrackingGeometryJsonConverter.Config()
    detray.addDetrayNavigationPolicies(jsonCfg)
    acts.json.TrackingGeometryJsonConverter(jsonCfg)
