import gc
import importlib.util
import subprocess
import sys
import weakref

import pytest
import acts

if importlib.util.find_spec("acts.ActsPluginsPythonBindingsDetray") is None:
    pytest.skip("ACTS detray plugin is not built", allow_module_level=True)


@pytest.mark.parametrize(
    "imports",
    [
        "import detray, detray.io; import acts.detray, acts.vecmem",
        "import acts.vecmem, acts.detray; import detray, detray.io",
    ],
)
def test_canonical_types_in_either_import_order(imports):
    subprocess.run(
        [
            sys.executable,
            "-c",
            imports
            + "; assert acts.detray.Detector is detray.Detector"
            + "; assert acts.vecmem.HostMemoryResource is detray.HostMemoryResource"
            + "; assert acts.vecmem.MemoryResource is detray.MemoryResource",
        ],
        check=True,
    )


@pytest.mark.parametrize("metadata_id", ["default", "odd"])
def test_conversion_io_and_allocator_lifetime(metadata_id, tmp_path):
    import acts.detray
    import detray
    import detray.io

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    root = acts.Blueprint(envelope=acts.ExtentEnvelope(r=[1, 1], z=[1, 1]))
    root.addStaticVolume(
        acts.Transform3.Identity(), acts.CylinderVolumeBounds(10, 20, 30), name="volume"
    )
    geometry = root.construct(acts.BlueprintOptions(), gctx, level=acts.logging.WARNING)
    payload_cfg = acts.detray.DetrayPayloadConverter.Config()
    payload_cfg.beampipeVolume = geometry.findVolumeByName("volume")
    converter_cfg = acts.detray.DetrayGeometryConverter.Config()
    converter_cfg.payloadConverter = acts.detray.DetrayPayloadConverter(
        payload_cfg, acts.logging.WARNING
    )
    converter = acts.detray.DetrayGeometryConverter(converter_cfg, acts.logging.WARNING)
    mr = detray.HostMemoryResource()
    mr_ref = weakref.ref(mr)
    geometry_ref = weakref.ref(geometry)
    result = converter.convert(mr, gctx, geometry, metadataId=metadata_id)
    detector, names = result.detector, result.names
    assert type(detector) is detray.Detector
    assert detector.metadataId == metadata_id
    count = len(detector.surfaces)
    assert count > 0
    del result, mr, geometry
    gc.collect()
    assert mr_ref() is not None
    assert geometry_ref() is not None
    assert len(detector.surfaces) == count

    writer = detray.io.DetectorWriterConfig()
    writer.path = str(tmp_path)
    detray.io.writeDetector(detector, names, writer)
    cfg = detray.io.DetectorReaderConfig().addFile(
        str(tmp_path / f"{detector.name(names)}_geometry.json")
    )
    restored, _ = detray.io.readDetector(detray.HostMemoryResource(), cfg)
    assert restored.metadataId == metadata_id
    assert len(restored.surfaces) == count
    del detector, names
    gc.collect()
    assert mr_ref() is None
    assert geometry_ref() is None
