import acts
import pytest

from helpers import detrayEnabled

u = acts.UnitConstants


@pytest.mark.skipif(not detrayEnabled, reason="Detray plugin not available")
def test_detray_propagator_bindings():
    import acts.examples.detray as exdetray

    assert hasattr(exdetray, "StraightLinePropagatorODD")
    assert hasattr(exdetray, "HelixPropagatorODD")
