import pytest
import warnings
from pathlib import Path

from helpers import dd4hepEnabled, geant4Enabled, rootEnabled

import acts.examples


def count_surfaces(geo):
    __tracebackhide__ = True
    nSurfaces = 0

    def visit(srf):
        nonlocal nSurfaces
        nSurfaces += 1

    geo.visitSurfaces(visit)

    return nSurfaces


def check_extra_odd(srf):
    if srf.geometryId.volume in [28, 30, 23, 25, 16, 18]:
        assert srf.geometryId.extra != 0
    return


def test_generic_geometry():
    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()
    contextDecorators = detector.contextDecorators()
    assert detector is not None
    assert trackingGeometry is not None
    assert contextDecorators is not None

    assert count_surfaces(trackingGeometry) == 18728


def test_telescope_geometry():
    n_surfaces = 10

    config = acts.examples.TelescopeDetector.Config(
        bounds=[100, 100],
        positions=[10 * i for i in range(n_surfaces)],
        stereos=[0] * n_surfaces,
        binValue=0,
    )
    detector = acts.examples.TelescopeDetector(config)
    trackingGeometry = detector.trackingGeometry()
    contextDecorators = detector.contextDecorators()

    assert detector is not None
    assert trackingGeometry is not None
    assert contextDecorators is not None

    assert count_surfaces(trackingGeometry) == n_surfaces


@pytest.mark.skipif(not geant4Enabled, reason="Geant4 is not set up")
@pytest.mark.parametrize("binValue", [0, 1, 2])
def test_telescope_geant4_geometry(binValue):
    from acts.examples.geant4 import SensitiveSurfaceMapper

    n_surfaces = 6

    config = acts.examples.TelescopeDetector.Config(
        bounds=[100, 200],
        positions=[30 * (i + 1) for i in range(n_surfaces)],
        stereos=[0] * n_surfaces,
        offsets=[10, -20],
        binValue=binValue,
    )
    detector = acts.examples.TelescopeDetector(config)
    trackingGeometry = detector.trackingGeometry()
    gctx = detector.nominalGeometryContext()

    # every sensitive surface has to be backed by a Geant4 volume in the same place
    smmConfig = SensitiveSurfaceMapper.Config()
    smmConfig.materialMappings = ["Silicon"]
    mapper = SensitiveSurfaceMapper.create(
        smmConfig, acts.logging.INFO, trackingGeometry
    )

    state = SensitiveSurfaceMapper.State()
    mapper.remapSensitiveNames(
        state, gctx, detector, acts.Transform3(acts.Vector3(0, 0, 0))
    )

    assert mapper.checkMapping(state, gctx, False, False)


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep is not set up")
def test_odd(odd_detector):
    detector = odd_detector
    trackingGeometry = detector.trackingGeometry()

    trackingGeometry.visitSurfaces(check_extra_odd)

    assert count_surfaces(trackingGeometry) == 18824


import itertools


def test_tgeo_config_triplet(monkeypatch):

    from acts.examples.tgeo import TGeoDetector, Interval

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
    from acts.examples.tgeo import TGeoDetector, Interval

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


@pytest.mark.skipif(not rootEnabled, reason="ROOT not set up")
def test_tgeo_detector_aligned_element_factory():
    import acts.examples.tgeo as tgeo
    from acts.examples.tgeo import TGeoDetector

    u = acts.UnitConstants
    Volume = TGeoDetector.Config.Volume
    LayerTriplet = TGeoDetector.Config.LayerTriplet

    root_file = (
        Path(__file__).parent.parent.parent.parent / "Tests" / "Data" / "panda.root"
    )
    assert root_file.exists()

    # Select the pixel barrel layer out of the panda.root geometry, same
    # selection as the `TGeoLayerBuilderTests` C++ unit test (`b0Config`).
    volume = Volume(
        name="Pixels",
        layers=LayerTriplet(negative=False, central=True, positive=False),
        subVolumeName=LayerTriplet(central="*"),
        sensitiveNames=LayerTriplet(
            central=[
                "PixelActiveo2",
                "PixelActiveo4",
                "PixelActiveo5",
                "PixelActiveo6",
            ]
        ),
        sensitiveAxes=LayerTriplet(central="XYZ"),
        rRange=LayerTriplet(central=(0.0, 40 * u.mm)),
        zRange=LayerTriplet(central=(-60 * u.mm, 15 * u.mm)),
        splitTolR=LayerTriplet(central=-1.0),
        splitTolZ=LayerTriplet(central=-1.0),
    )

    cfg = TGeoDetector.Config(
        fileName=str(root_file),
        unitScalor=1 * u.cm,  # panda.root stores lengths in ROOT's native cm
        volumes=[volume],
        detectorElementFactory=tgeo.alignedTGeoDetectorElementFactory,
    )

    # The setter dispatches on identity to a native C++ function pointer
    # rather than storing the assigned python object, so the getter always
    # reports None.
    assert cfg.detectorElementFactory is None

    # (Re-)building a TGeoManager from a ROOT file can emit routine,
    # benign ROOT diagnostics (e.g. about registered matrices being
    # replaced) that are not test failures.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        detector = TGeoDetector(cfg)
        trackingGeometry = detector.trackingGeometry()

    assert trackingGeometry is not None
    assert count_surfaces(trackingGeometry) == 14


@pytest.mark.skipif(not rootEnabled, reason="ROOT not set up")
def test_tgeo_detector_element_factory_invalid():
    from acts.examples.tgeo import TGeoDetector

    cfg = TGeoDetector.Config()

    with pytest.raises(ValueError):
        cfg.detectorElementFactory = lambda *args: None

    # None is accepted and leaves the default (non-aligned) factory in place
    cfg.detectorElementFactory = None


def test_coordinate_converter(trk_geo):
    from acts.examples import json

    digiCfg = acts.examples.DigitizationAlgorithm.Config(
        digitizationConfigs=acts.examples.json.readDigiConfigFromJson(
            str(
                Path(__file__).parent.parent.parent.parent
                / "Examples/Configs/generic-digi-smearing-config.json"
            )
        ),
        surfaceByIdentifier=trk_geo.geoIdSurfaceMap(),
    )
    converter = acts.examples.DigitizationCoordinatesConverter(digiCfg)

    def test_surface(surface):
        gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
        geo_id = surface.geometryId.value
        geo_center = surface.center(gctx)
        x, y, z = geo_center[0], geo_center[1], geo_center[2]

        # test if surface center can be reproduced
        assert converter.globalToLocal(geo_id, x, y, z) == (0, 0)
        assert converter.localToGlobal(geo_id, 0, 0) == (x, y, z)

        # test if we can get back to the same local coordinates
        global_shifted = converter.localToGlobal(geo_id, 5, 5)
        local_shifted = converter.globalToLocal(geo_id, *global_shifted)
        assert abs(local_shifted[0] - 5) / 5 < 1e-6
        assert abs(local_shifted[1] - 5) / 5 < 1e-6

    trk_geo.visitSurfaces(test_surface)
