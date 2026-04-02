import pytest
import acts
import functools
from acts.examples import GenericDetector
from acts.examples.odd import getOpenDataDetector
import json

from helpers import dd4hepEnabled


@pytest.mark.parametrize(
    "detectorFactory,aligned,nobj",
    [
        (functools.partial(GenericDetector, gen3=False), True, 2),
        pytest.param(
            functools.partial(GenericDetector, gen3=True),
            True,
            2,  # Gen3 geometry visualiztion produces a single file + materials
        ),
        pytest.param(
            getOpenDataDetector,
            True,
            2,
            marks=[
                pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
                pytest.mark.slow,
                pytest.mark.odd,
            ],
        ),
    ],
    ids=[
        "generic",
        "generic-gen3",
        "odd",
    ],
)
@pytest.mark.slow
def test_geometry_example(detectorFactory, aligned, nobj, tmp_path):
    detector = detectorFactory()
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    from geometry import runGeometry

    json_dir = tmp_path / "json"
    csv_dir = tmp_path / "csv"
    obj_dir = tmp_path / "obj"

    for d in (json_dir, csv_dir, obj_dir):
        d.mkdir()

    events = 5

    kwargs = dict(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        events=events,
        outputDir=tmp_path,
    )

    runGeometry(outputJson=True, **kwargs)
    runGeometry(outputJson=False, **kwargs)

    assert len(list(obj_dir.iterdir())) == nobj

    assert len(list(csv_dir.iterdir())) == 3 * events

    detector_files = [csv_dir / f"event{i:>09}-detectors.csv" for i in range(events)]
    for detector_file in detector_files:
        assert detector_file.exists()
        assert detector_file.stat().st_size > 200

    contents = [f.read_text() for f in detector_files]
    ref = contents[0]
    for c in contents[1:]:
        if aligned:
            assert c == ref, "Detector writeout is expected to be identical"
        else:
            assert c != ref, "Detector writeout is expected to be different"

    if aligned:
        for f in [json_dir / f"event{i:>09}-detector.json" for i in range(events)]:
            assert detector_file.exists()
            with f.open() as fh:
                data = json.load(fh)
                assert data
        material_file = tmp_path / "geometry-map.json"
        assert material_file.exists()
        assert material_file.stat().st_size > 200


class CountingVisitor(acts.TrackingGeometryMutableVisitor):
    def __init__(self):
        super().__init__()
        self.num_surfaces = 0
        self.num_layers = 0
        self.num_volumes = 0
        self.num_portals = 0
        self.num_boundary_surfaces = 0

    def visitSurface(self, surface: acts.Surface):
        self.num_surfaces += 1

    def visitLayer(self, layer: acts.Layer):
        self.num_layers += 1

    def visitVolume(self, volume: acts.Volume):
        self.num_volumes += 1

    def visitPortal(self, portal: acts.Portal):
        self.num_portals += 1

    def visitBoundarySurface(self, boundary: acts.BoundarySurfaceT_TrackingVolume):
        self.num_boundary_surfaces += 1


def test_geometry_visitor(trk_geo):
    visitor = CountingVisitor()
    trk_geo.apply(visitor)

    assert visitor.num_surfaces == 19078
    assert visitor.num_layers == 111
    assert visitor.num_volumes == 18
    assert visitor.num_portals == 0
    assert visitor.num_boundary_surfaces == 67


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.odd
def test_odd_gen1():
    with getOpenDataDetector(gen3=False) as detector:
        trackingGeometry = detector.trackingGeometry()

        visitor = CountingVisitor()
        trackingGeometry.apply(visitor)

        assert visitor.num_surfaces == 19264
        assert visitor.num_layers == 142
        assert visitor.num_volumes == 32
        assert visitor.num_portals == 0  # Gen1: will have no portals
        assert visitor.num_boundary_surfaces == 126  # Gen1: will have boundary surfaces


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.odd
@pytest.mark.parametrize(
    "constructionMethod",
    [
        pytest.param(None, id="default"),
        pytest.param("BarrelEndcap", id="barrel-endcap"),
        pytest.param("DirectLayer", id="direct-layer"),
        pytest.param("DirectLayerGrouped", id="direct-layer-grouped"),
    ],
)
def test_odd_gen3(constructionMethod):
    import acts.examples.dd4hep as dd4hep

    cm = None
    if constructionMethod is not None:
        cm = getattr(
            dd4hep.OpenDataDetector.Config.ConstructionMethod, constructionMethod
        )

    with getOpenDataDetector(gen3=True, constructionMethod=cm) as detector:
        trackingGeometry = detector.trackingGeometry()

        visitor = CountingVisitor()
        trackingGeometry.apply(visitor)

        # Gen3 invariants that hold regardless of construction method
        assert visitor.num_layers == 0  # Gen3: no layers
        assert visitor.num_boundary_surfaces == 0  # Gen3: no boundary surfaces
        assert visitor.num_portals > 0  # Gen3: uses portals instead
        assert visitor.num_surfaces > 0
        assert visitor.num_volumes > 0

        assert visitor.num_surfaces == 19261
        assert visitor.num_volumes == 109
        assert visitor.num_portals == 437
