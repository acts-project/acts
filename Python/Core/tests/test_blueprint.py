import pytest

import acts

mm = acts.UnitConstants.mm
m = acts.UnitConstants.m
degree = acts.UnitConstants.degree

bv = acts.AxisDirection

gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
logLevel = acts.logging.VERBOSE


def test_zdirection_container_blueprint(tmp_path):

    def write(root: acts.BlueprintNode, stage: int):
        gz = tmp_path / f"blueprint_{stage}.dot"
        print(gz)
        with gz.open("w") as fh:
            root.graphviz(fh)

    base = acts.Transform3.Identity()

    root = acts.Blueprint(envelope=acts.ExtentEnvelope(r=[10 * mm, 10 * mm]))
    assert root.depth == 0

    barrel = root.addCylinderContainer("Barrel", direction=bv.AxisR)

    assert barrel.depth == 1

    r = 25 * mm
    for i in range(1, 3):
        r += 50 * mm
        bounds = acts.CylinderVolumeBounds(r, r + 20 * mm, 200 * mm)
        vol = barrel.addStaticVolume(base, bounds, name=f"Barrel_{i}")
        assert vol.depth == 2

    write(barrel, 1)

    root.clearChildren()

    assert barrel.depth == 0

    det = root.addCylinderContainer("Detector", direction=bv.AxisZ)

    assert det.depth == 1

    with det.CylinderContainer("nEC", direction=bv.AxisZ) as ec:
        assert ec.depth == 2
        z = -200
        for i in range(1, 3):
            z -= 200 * mm
            bounds = acts.CylinderVolumeBounds(100 * mm, 150 * mm, 50 * mm)

            trf = base * acts.Translation3(acts.Vector3(0, 0, z))

            vol = ec.addStaticVolume(trf, bounds, name=f"nEC_{i}")
            assert vol.depth == 3

        write(ec, 2)

    det.addChild(barrel)
    assert barrel.depth == 2

    write(det, 3)

    with det.CylinderContainer("pEC", direction=bv.AxisZ) as ec:
        assert ec.depth == 2
        z = 200
        for i in range(1, 3):
            z += 200 * mm
            bounds = acts.CylinderVolumeBounds(100 * mm, 150 * mm, 50 * mm)

            trf = base * acts.Translation3(acts.Vector3(0, 0, z))

            vol = ec.addStaticVolume(trf, bounds, name=f"pEC_{i}")
            assert vol.depth == 3

    write(root, 4)


def test_optional_material_keys(tmp_path):
    import json

    acts_json = pytest.importorskip("acts.json")
    root = acts.Blueprint(envelope=acts.ExtentEnvelope(r=[1 * mm, 2 * mm]))
    designator = root.addMaterial("material")
    axis = acts.AxisSpec.DeferredEquidistant(4)
    designator.configureFace(
        acts.CylinderVolumeBounds.Face.OuterCylinder,
        axis,
        axis,
        materialKey="barrel/outer",
    )
    # Existing calls and explicit None remain valid.
    designator.configureFace(acts.CylinderVolumeBounds.Face.NegativeDisc, axis, axis)
    designator.configureFace(
        acts.CylinderVolumeBounds.Face.PositiveDisc, axis, axis, materialKey=None
    )
    with pytest.raises(ValueError, match="already configured"):
        designator.configureFace(
            acts.CylinderVolumeBounds.Face.OuterCylinder, axis, axis, "duplicate"
        )
    designator.addStaticVolume(
        acts.Transform3.Identity(),
        acts.CylinderVolumeBounds(10 * mm, 20 * mm, 30 * mm),
        name="barrel",
    )
    geometry = root.construct(acts.BlueprintOptions(), gctx, level=acts.logging.WARNING)
    keyed = {}

    def collect(surface):
        material = surface.surfaceMaterial
        if (
            isinstance(
                material, (acts.ProtoSurfaceMaterial, acts.ProtoGridSurfaceMaterial)
            )
            and material.materialKey is not None
        ):
            keyed[material.materialKey] = surface

    geometry.visitSurfaces(collect, False)
    assert set(keyed) == {"barrel/outer"}
    assert keyed["barrel/outer"].surfaceMaterial.materialKey == "barrel/outer"

    def legacy_section(name):
        return {
            "acts-geometry-hierarchy-map": {
                "format-version": 0,
                "value-identifier": name,
            },
            "entries": [],
        }

    payload = {
        "type": "homogeneous",
        "mapMaterial": True,
        "data": [[{"material": None, "thickness": 0.0}]],
    }
    document = {
        "Surfaces": legacy_section("Material Surface Map"),
        "Volumes": legacy_section("Material Volume Map"),
        "KeyedSurfaces": [
            {"key": key, "geometry_id": 999, "material": payload}
            for key in ["barrel/outer", "unused/other-detector"]
        ],
    }
    path = tmp_path / "material.json"
    path.write_text(json.dumps(document))
    loader = acts_json.JsonMaterialDecorator(
        acts_json.MaterialMapJsonConverter.Config(), str(path), acts.logging.WARNING
    )
    loader.materialMaps.apply(geometry)
    assert isinstance(
        keyed["barrel/outer"].surfaceMaterial, acts.HomogeneousSurfaceMaterial
    )
