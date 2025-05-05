import pytest
import acts
import functools
from acts.examples import GenericDetector, AlignedDetector
from acts.examples.odd import getOpenDataDetector
import json

from helpers import dd4hepEnabled


@pytest.mark.parametrize(
    "detectorFactory,aligned,nobj",
    [
        (GenericDetector, True, 450),
        pytest.param(
            getOpenDataDetector,
            True,
            540,
            marks=[
                pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
                pytest.mark.slow,
                pytest.mark.odd,
            ],
        ),
        (functools.partial(AlignedDetector, iovSize=1), False, 450),
    ],
    ids=["generic", "odd", "aligned"],
)
@pytest.mark.slow
def test_geometry_example(detectorFactory, aligned, nobj, tmp_path):
    detector = detectorFactory(decoratorLogLevel=acts.logging.VERBOSE)
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
        outputDir=str(tmp_path),
    )

    runGeometry(outputJson=True, **kwargs)
    runGeometry(outputJson=False, **kwargs)

    assert len(list(obj_dir.iterdir())) == nobj
    assert all(f.stat().st_size > 200 for f in obj_dir.iterdir())

    assert len(list(csv_dir.iterdir())) == 3 * events
    assert all(f.stat().st_size > 200 for f in csv_dir.iterdir())

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


def test_geometry_visitor(trk_geo):
    class Visitor(acts.TrackingGeometryMutableVisitor):
        def __init__(self):
            super().__init__()
            self.num_surfaces = 0
            self.num_layers = 0
            self.num_volumes = 0
            self.num_portals = 0

        def visitSurface(self, surface: acts.Surface):
            self.num_surfaces += 1

        def visitLayer(self, layer: acts.Layer):
            self.num_layers += 1

        def visitVolume(self, volume: acts.Volume):
            self.num_volumes += 1

        def visitPortal(self, portal: acts.Portal):
            self.num_portals += 1

    visitor = Visitor()
    trk_geo.apply(visitor)

    assert visitor.num_surfaces == 19078
    assert visitor.num_layers == 111
    assert visitor.num_volumes == 18
    assert visitor.num_portals == 0


def test_processor_blueprint():
    cfg = acts.Blueprint.Config()
    cfg.envelope[acts.AxisDirection.AxisZ] = (20, 20)
    cfg.envelope[acts.AxisDirection.AxisR] = (1, 2)
    root = acts.Blueprint(cfg)

    # Create a simple cylinder volume
    bounds = acts.CylinderVolumeBounds(10, 20, 30)
    transform = acts.Transform3()

    build_called = False
    connect_called = False
    finalize_called = False

    def on_build(vol):
        print("ON BUILD")
        nonlocal build_called
        build_called = True

    def on_connect(opts, gctx, logger, shell):
        nonlocal connect_called
        connect_called = True
        return shell

    def on_finalize(opts, gctx, parent, logger):
        nonlocal finalize_called
        finalize_called = True

    # Add processor node with callbacks
    processor = root.withProcessor()
    processor.name = "Processor"
    processor.onBuild = on_build
    with pytest.raises(AttributeError):
        print(processor.onBuild)
    processor.onConnect = on_connect
    with pytest.raises(AttributeError):
        print(processor.onConnect)
    processor.onFinalize = on_finalize
    with pytest.raises(AttributeError):
        print(processor.onFinalize)
    processor.addChild(acts.StaticBlueprintNode(transform, bounds, "child"))

    # Construct geometry
    gctx = acts.GeometryContext()
    tracking_geometry = root.construct(
        acts.BlueprintOptions(), gctx, acts.logging.VERBOSE
    )

    assert build_called
    assert connect_called
    assert finalize_called

    # # Test error cases
    # root2 = acts.Blueprint(cfg)
    # processor2 = root2.withProcessor()
    # processor2.setName("Processor")

    # # Should throw when no child is added
    # with pytest.raises(RuntimeError):
    #     root2.construct({}, gctx, logger)

    # # Should throw when multiple children are added
    # vol1 = acts.TrackingVolume(transform, bounds, "child1")
    # vol2 = acts.TrackingVolume(transform, bounds, "child2")
    # processor2.addChild(acts.StaticBlueprintNode(vol1))
    # processor2.addChild(acts.StaticBlueprintNode(vol2))
    # with pytest.raises(RuntimeError):
    #     root2.construct({}, gctx, logger)


#
