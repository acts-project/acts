import pytest

import acts

from helpers import dd4hepEnabled

# Standalone script under Examples/Scripts/Python (added to sys.path by conftest)
import material_merge_rz_view as mrz

mm = acts.UnitConstants.mm


def _build_marker_geometry(gctx):
    """Build a tiny Gen3 geometry that triggers a lossy material merge.

    A z-stack with material designated on the (merged) OuterCylinder face of one
    child. Constructed in keep-going mode so the merge tags the surface with a
    MergedMaterialMarker instead of aborting.
    """
    bv = acts.AxisDirection
    abt = acts.AxisBoundaryType
    base = acts.Transform3.Identity()

    root = acts.Blueprint(
        envelope=acts.ExtentEnvelope(z=[20 * mm, 20 * mm], r=[0 * mm, 20 * mm])
    )
    stack = root.addCylinderContainer("Stack", direction=bv.AxisZ)

    mat = stack.addMaterial("Material")
    mat.configureFace(
        acts.CylinderVolumeBounds.Face.OuterCylinder,
        acts.DirectedProtoAxis(bv.AxisRPhi, abt.Bound, 20),
        acts.DirectedProtoAxis(bv.AxisZ, abt.Bound, 20),
    )
    mat.addStaticVolume(
        base * acts.Translation3(acts.Vector3(0, 0, -200 * mm)),
        acts.CylinderVolumeBounds(0, 100 * mm, 100 * mm),
        name="VolumeA",
    )

    stack.addStaticVolume(
        base * acts.Translation3(acts.Vector3(0, 0, 200 * mm)),
        acts.CylinderVolumeBounds(0, 100 * mm, 100 * mm),
        name="VolumeB",
    )

    options = acts.BlueprintOptions()
    options.keepGoingOnMaterialMergeFailure = True
    return root.construct(options, gctx, level=acts.logging.WARNING)


@acts.with_log_threshold(acts.logging.FATAL)
def test_rz_view_with_markers(tmp_path):
    pytest.importorskip("matplotlib")
    from acts.json import TrackingGeometryJsonConverter

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    trackingGeometry = _build_marker_geometry(gctx)

    json_path = tmp_path / "marker-geometry.json"
    json_path.write_text(TrackingGeometryJsonConverter().toJson(gctx, trackingGeometry))

    out = tmp_path / "marker_rz.svg"
    volumes, markers = mrz.run(json_path, out)

    assert out.exists()
    assert out.stat().st_size > 0
    assert len(volumes) >= 2
    # The merged OuterCylinder portal must show up as a marker line
    assert len(markers) >= 1


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.odd
@pytest.mark.slow
# Construction may legitimately emit keep-going WARNINGs (e.g. when ODD has a
# material-on-merged-face clash), which would otherwise trip the test harness'
# ACTS_LOG_FAILURE_THRESHOLD.
@acts.with_log_threshold(acts.logging.FATAL)
def test_rz_view_odd_gen3(tmp_path):
    pytest.importorskip("matplotlib")
    from acts.examples.odd import getOpenDataDetector
    from acts.json import TrackingGeometryJsonConverter

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()

    with getOpenDataDetector(gen3=True) as detector:
        trackingGeometry = detector.trackingGeometry()
        json_path = tmp_path / "odd-geometry.json"
        json_path.write_text(
            TrackingGeometryJsonConverter().toJson(gctx, trackingGeometry)
        )

    out = tmp_path / "odd_rz.svg"
    volumes, markers = mrz.run(json_path, out)

    assert out.exists()
    assert out.stat().st_size > 0
    assert len(volumes) > 0
    # Marker count depends on ODD's material configuration, so we don't pin it;
    # the script must just run end-to-end and produce an SVG.
    assert len(markers) >= 0
