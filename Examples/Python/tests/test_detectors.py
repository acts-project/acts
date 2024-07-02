import pytest

from helpers import dd4hepEnabled

import acts.examples
from acts.examples.odd import getOpenDataDetector


def count_surfaces(geo):
    __tracebackhide__ = True
    nSurfaces = 0

    def visit(srf):
        nonlocal nSurfaces
        nSurfaces += 1

    geo.visitSurfaces(visit)

    return nSurfaces


def check_extra_odd(srf):
    if srf.geometryId().volume() in [28, 30, 23, 25, 16, 18]:
        assert srf.geometryId().extra() != 0
    return


def test_generic_geometry():
    detector, geo, contextDecorators = acts.examples.GenericDetector.create()
    assert detector is not None
    assert geo is not None
    assert contextDecorators is not None

    assert count_surfaces(geo) == 18728


def test_telescope_geometry():
    n_surfaces = 10

    detector, geo, contextDecorators = acts.examples.TelescopeDetector.create(
        bounds=[100, 100],
        positions=[10 * i for i in range(n_surfaces)],
        stereos=[0] * n_surfaces,
        binValue=0,
    )

    assert detector is not None
    assert geo is not None
    assert contextDecorators is not None

    assert count_surfaces(geo) == n_surfaces


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep is not set up")
def test_odd():
    detector, trackingGeometry, decorators = getOpenDataDetector()

    trackingGeometry.visitSurfaces(check_extra_odd)

    assert count_surfaces(trackingGeometry) == 18824


def test_aligned_detector():
    detector, trackingGeometry, decorators = acts.examples.AlignedDetector.create()

    assert detector is not None
    assert trackingGeometry is not None
    assert decorators is not None

    assert count_surfaces(trackingGeometry) == 18728


import itertools


def test_tgeo_config_triplet(monkeypatch):
    from acts.examples import TGeoDetector, Interval

    # monkeypatch the comparison operator
    def eq(self, other):
        return self.lower == other.lower and self.upper == other.upper

    monkeypatch.setattr(Interval, "__eq__", eq)

    LayerTriplet = TGeoDetector.Config.LayerTriplet
    c = TGeoDetector.Config

    def assert_combinations(value, _type):
        t = LayerTriplet(value)
        assert t.negative == value and t.central == value and t.positive == value
        assert isinstance(t, _type)

        keys = ["negative", "central", "positive"]

        combinations = (
            [(k,) for k in keys] + list(itertools.combinations(keys, 2)) + [keys]
        )

        for c in combinations:
            d = {k: value for k in c}

            t = LayerTriplet(**d)
            assert isinstance(t, _type)
            for k in c:
                assert getattr(t, k) == value

    v = ["Some::SensorName"]
    assert_combinations(v, c.LayerTripletVectorString)

    with pytest.raises(TypeError):
        LayerTriplet(["Some::SensorName", 848])

    with pytest.raises(TypeError):
        LayerTriplet(("Some::SensorName", 848))

    for v in (True, False):
        assert_combinations(v, c.LayerTripletBool)

    assert_combinations("hallo", c.LayerTripletString)

    assert_combinations(5.3, c.LayerTripletDouble)

    assert_combinations(Interval(5.0, 9.0), c.LayerTripletInterval)

    with pytest.raises(TypeError):
        LayerTriplet(("a", 9))

    v = (4.4, 2.2)
    t = LayerTriplet(v)
    assert t.negative == Interval(*v)
    assert t.central == Interval(*v)
    assert t.positive == Interval(*v)


def test_tgeo_config_volume(monkeypatch):
    from acts.examples import TGeoDetector, Interval

    # monkeypatch the comparison operator
    def eq(self, other):
        return self.lower == other.lower and self.upper == other.upper

    monkeypatch.setattr(Interval, "__eq__", eq)

    Volume = TGeoDetector.Config.Volume

    v = Volume(name="blubb")
    assert v

    for key in ("binToleranceR", "binToleranceZ", "binTolerancePhi"):
        v = Volume(**{key: Interval(4, 5)})
        assert getattr(v, key) == Interval(4, 5)

        v = Volume(**{key: (4, 5)})
        assert getattr(v, key) == Interval(4, 5)

        v = Volume(**{key: (None, 5)})
        assert getattr(v, key) == Interval(None, 5)

        v = Volume(**{key: (4, None)})
        assert getattr(v, key) == Interval(4, None)
