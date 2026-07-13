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

    runGeometry(outputSurfacesJson=True, **kwargs)
    runGeometry(outputSurfacesJson=False, **kwargs)

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
        pytest.param("TGeo", id="tgeo"),
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


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.odd
def test_odd_gen3_json_roundtrip(tmp_path):
    from geometry import runGeometry
    from acts.json import TrackingGeometryJsonConverter

    json_dir = tmp_path / "json"
    json_dir.mkdir()

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()

    with getOpenDataDetector(gen3=True) as detector:
        trackingGeometry = detector.trackingGeometry()

        original = CountingVisitor()
        trackingGeometry.apply(original)

        runGeometry(
            trackingGeometry=trackingGeometry,
            decorators=detector.contextDecorators(),
            outputDir=tmp_path,
            events=1,
            outputObj=False,
            outputCsv=False,
            outputSurfacesJson=False,
            serializeGeometryJson=True,
        )

    json_path = json_dir / "tracking-geometry.json"
    assert json_path.exists()
    assert json_path.stat().st_size > 0

    converter = TrackingGeometryJsonConverter()
    rebuilt_geometry = converter.fromJson(gctx, json_path.absolute())

    rebuilt = CountingVisitor()
    rebuilt_geometry.apply(rebuilt)

    assert rebuilt.num_surfaces == original.num_surfaces
    assert rebuilt.num_volumes == original.num_volumes
    assert rebuilt.num_portals == original.num_portals
    assert rebuilt.num_layers == 0
    assert rebuilt.num_boundary_surfaces == 0

    # Propagation comparison: identical seed → identical navigation results
    import numpy as np
    import uproot
    from propagation import runPropagation

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    def run_prop(geo, out_dir):
        out_dir.mkdir(exist_ok=True)
        s = acts.examples.Sequencer(
            events=5, numThreads=1, logLevel=acts.logging.WARNING
        )
        runPropagation(geo, field, str(out_dir), s=s, sterileLogger=False).run()

    orig_dir = tmp_path / "prop_original"
    rebuilt_dir = tmp_path / "prop_rebuilt"
    run_prop(trackingGeometry, orig_dir)
    run_prop(rebuilt_geometry, rebuilt_dir)

    branches = ["nSensitives", "nPortals", "nSuccessfulSteps", "pathLength"]
    with uproot.open(orig_dir / "propagation_summary.root") as f:
        orig = f["propagation_summary"].arrays(branches, library="np")
    with uproot.open(rebuilt_dir / "propagation_summary.root") as f:
        rebuilt_data = f["propagation_summary"].arrays(branches, library="np")

    assert orig["nSensitives"].mean() > 0.0
    assert orig["nPortals"].mean() > 0.0
    assert orig["pathLength"].mean() > 0.0

    np.testing.assert_array_equal(orig["nSensitives"], rebuilt_data["nSensitives"])
    np.testing.assert_array_equal(orig["nPortals"], rebuilt_data["nPortals"])
    np.testing.assert_array_equal(
        orig["nSuccessfulSteps"], rebuilt_data["nSuccessfulSteps"]
    )
    np.testing.assert_allclose(
        orig["pathLength"], rebuilt_data["pathLength"], rtol=1e-5
    )


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.odd
@pytest.mark.json
def test_odd_gen3_json_material(tmp_path):
    """
    Validates that the JSON-only workflow (topology-only geometry JSON +
    JsonMaterialDecorator reading a standalone material-map JSON) reproduces
    the same accumulated material (X0/L0) as the reference DD4hep +
    RootMaterialDecorator workflow, including material assigned to volumes
    (not just surfaces).
    """
    import acts.examples
    from acts.json import (
        TrackingGeometryJsonConverter,
        MaterialMapJsonConverter,
        JsonMaterialDecorator,
    )
    from acts.examples.json import JsonMaterialWriter, JsonFormat
    from material_validation import runMaterialValidation

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()

    # Gen3 geometry, decorated with the default RootMaterialDecorator (the
    # reference workflow) -- this is both the ground truth to compare against
    # and the source geometry for the topology-only JSON below.
    with getOpenDataDetector(gen3=True) as detector:
        trackingGeometry = detector.trackingGeometry()

        # Geometry JSON with no material embedded: geometry description and
        # material are orthogonal, exactly like DD4hep XML + RootMaterialDecorator.
        converter = TrackingGeometryJsonConverter()
        options = TrackingGeometryJsonConverter.Options()
        options.writeMaterial = False
        geometry_json_path = tmp_path / "tracking-geometry.json"
        geometry_json_path.write_text(converter.toJson(gctx, trackingGeometry, options))

        # Material map JSON extracted from the already-decorated reference
        # geometry (includes both surface and volume material).
        material_map_base = tmp_path / "material-map"
        jmw = JsonMaterialWriter(
            level=acts.logging.INFO,
            converterCfg=MaterialMapJsonConverter.Config(
                processSensitives=True,
                processApproaches=True,
                processRepresenting=True,
                processBoundaries=True,
                processVolumes=True,
                processNonMaterial=True,
                context=gctx,
            ),
            fileName=str(material_map_base),
            writeFormat=JsonFormat.Json,
        )
        jmw.write(trackingGeometry)
        material_map_path = tmp_path / "material-map.json"
        assert material_map_path.exists()

        # Rebuild geometry from the topology-only JSON, applying material via
        # JsonMaterialDecorator independently, mirroring RootMaterialDecorator.
        materialDecorator = JsonMaterialDecorator(
            rConfig=MaterialMapJsonConverter.Config(),
            jFileName=str(material_map_path),
            level=acts.logging.INFO,
        )
        rebuilt_geometry = converter.fromJson(
            gctx, geometry_json_path.absolute(), materialDecorator=materialDecorator
        )

    def run_validation(geo, out_dir):
        out_dir.mkdir()
        s = acts.examples.Sequencer(
            events=5, numThreads=1, logLevel=acts.logging.WARNING
        )
        runMaterialValidation(
            surfaces=geo.extractMaterialSurfaces(),
            s=s,
            tracksPerEvent=100,
            outputFileBase=str(out_dir / "material_validation"),
        )
        s.run()

    orig_dir = tmp_path / "val_original"
    rebuilt_dir = tmp_path / "val_rebuilt"
    run_validation(trackingGeometry, orig_dir)
    run_validation(rebuilt_geometry, rebuilt_dir)

    import numpy as np
    import uproot

    with uproot.open(orig_dir / "material_validation.root") as f:
        orig = f["material_tracks"].arrays(["t_X0", "t_L0"], library="np")
    with uproot.open(rebuilt_dir / "material_validation.root") as f:
        rebuilt = f["material_tracks"].arrays(["t_X0", "t_L0"], library="np")

    assert orig["t_X0"].mean() > 0.0
    assert orig["t_L0"].mean() > 0.0

    np.testing.assert_allclose(orig["t_X0"], rebuilt["t_X0"], rtol=1e-9)
    np.testing.assert_allclose(orig["t_L0"], rebuilt["t_L0"], rtol=1e-9)
