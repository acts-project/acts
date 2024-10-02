#!/usr/bin/env python3

from pathlib import Path

import acts

mm = acts.UnitConstants.mm
degree = acts.UnitConstants.degree

root = acts.RootBlueprintNode(envelope=acts.ExtentEnvelope(r=[10 * mm, 10 * mm]))


pixel = root.addCylinderContainer(direction=acts.BinningValue.binZ, name="Pixel")
print(repr(pixel))

trf = acts.Transform3.Identity() * acts.Translation3(acts.Vector3(0, 0, 0 * mm))


if True:
    barrel = acts.CylinderContainerBlueprintNode(
        "PixelBarrel",
        acts.BinningValue.binR,
        attachmentStrategy=acts.CylinderVolumeStack.AttachmentStrategy.Gap,
        resizeStrategy=acts.CylinderVolumeStack.ResizeStrategy.Gap,
    )
    pixel.addChild(barrel)

    print("Barrel")
    r = 25 * mm
    for i in range(0, 4):
        r += 50 * mm
        bounds = acts.CylinderVolumeBounds(r, r + 20 * mm, 200 * mm)
        print(bounds)
        brlLayer = barrel.addStaticVolume(trf, bounds, name=f"PixelBarrelLayer{i}")
        assert brlLayer.name == f"PixelBarrelLayer{i}"


if True:

    @pixel.CylinderContainer(acts.BinningValue.binZ)
    def PixelPosEndcap(ec):
        print("Positive Endcap")

        ec.attachmentStrategy = acts.CylinderVolumeStack.AttachmentStrategy.Gap
        ec.resizeStrategy = acts.CylinderVolumeStack.ResizeStrategy.Gap

        z = 200
        for i in range(0, 4):
            z += 200 * mm
            bounds = acts.CylinderVolumeBounds(100 * mm, 150 * mm, 50 * mm)
            print(bounds)

            trf = acts.Transform3.Identity() * acts.Translation3(acts.Vector3(0, 0, z))

            @ec.StaticVolume(trf, bounds, name=f"PixelPosEndcapDisk{i}")
            def PixelEndcapDisk(disk):
                print("Add disk", i)

            assert PixelEndcapDisk.name == f"PixelPosEndcapDisk{i}"


if True:
    with pixel.Material() as mat:
        with mat.CylinderContainer(acts.BinningValue.binZ, name="PixelNegEndcap") as ec:
            ec.attachmentStrategy = acts.CylinderVolumeStack.AttachmentStrategy.Gap

            print("Negative Endcap")

            z = -200
            for i in range(0, 4):
                z -= 200 * mm
                bounds = acts.CylinderVolumeBounds(200 * mm, 300 * mm, 50 * mm)
                print(bounds)

                trf = acts.Transform3.Identity() * acts.Translation3(
                    acts.Vector3(0, 0, z)
                )

                with ec.StaticVolume(
                    trf, bounds, name=f"PixelNegEndcapDisk{i}"
                ) as disk:
                    print("Add disk", i)
                    assert disk.name == f"PixelNegEndcapDisk{i}"


with open("blueprint.dot", "w") as fh:
    root.graphViz(fh)


gctx = acts.GeometryContext()
trackingGeometry = root.construct(
    options=acts.BlueprintNode.Options(), gctx=gctx, level=acts.logging.VERBOSE
)

vis = acts.ObjVisualization3D()
trackingGeometry.visualize(vis, gctx)
with Path("blueprint.obj").open("w") as fh:
    vis.write(fh)
# print("DONE")
