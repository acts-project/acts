import pytest

import acts

mm = acts.UnitConstants.mm
m = acts.UnitConstants.m
degree = acts.UnitConstants.degree

bv = acts.BinningValue

gctx = acts.GeometryContext()
logLevel = acts.logging.VERBOSE


def test_zdirection_container_blueprint(tmp_path):

    def write(root: acts.BlueprintNode, stage: int):
        gz = tmp_path / f"blueprint_{stage}.dot"
        print(gz)
        with gz.open("w") as fh:
            root.graphViz(fh)

    base = acts.Transform3.Identity()

    root = acts.RootBlueprintNode(envelope=acts.ExtentEnvelope(r=[10 * mm, 10 * mm]))

    barrel = root.addCylinderContainer("Barrel", direction=bv.binR)

    r = 25 * mm
    for i in range(1, 3):
        r += 50 * mm
        bounds = acts.CylinderVolumeBounds(r, r + 20 * mm, 200 * mm)
        barrel.addStaticVolume(base, bounds, name=f"Barrel_{i}")

    write(barrel, 1)

    root.clearChildren()

    det = root.addCylinderContainer("Detector", direction=bv.binZ)

    with det.CylinderContainer("nEC", direction=bv.binZ) as ec:
        z = -200
        for i in range(1, 3):
            z -= 200 * mm
            bounds = acts.CylinderVolumeBounds(100 * mm, 150 * mm, 50 * mm)

            trf = base * acts.Translation3(acts.Vector3(0, 0, z))

            ec.addStaticVolume(trf, bounds, name=f"nEC_{i}")

        write(ec, 2)

    det.addChild(barrel)

    write(det, 3)

    with det.CylinderContainer("pEC", direction=bv.binZ) as ec:
        z = 200
        for i in range(1, 3):
            z += 200 * mm
            bounds = acts.CylinderVolumeBounds(100 * mm, 150 * mm, 50 * mm)

            trf = base * acts.Translation3(acts.Vector3(0, 0, z))

            ec.addStaticVolume(trf, bounds, name=f"pEC_{i}")

    write(root, 4)

    # gz = tmp_path / "blueprint.dot"
    # print(gz)
    # with gz.open("w") as fh:
    #     root.graphViz(fh)
    #
    # trackingGeometry = root.construct(
    #     acts.BlueprintNode.Options(), gctx, level=logLevel
    # )
    # assert trackingGeometry is not None


if False:
    with root.CylinderContainer("Detector", direction=bv.binZ) as det:
        det.addStaticVolume(
            base, acts.CylinderVolumeBounds(0, 23 * mm, 3 * m), "BeamPipe"
        )

    trackingGeometry = root.construct(
        acts.BlueprintNode.Options(), gctx, level=logLevel
    )
    assert trackingGeometry is not None
