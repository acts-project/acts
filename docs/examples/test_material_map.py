def test_load_material_map():
    #! [Load Material Map]
    import acts
    from acts.examples.odd import getOpenDataDetectorDirectory, getOpenDataDetector

    # Gen1 example: any of .json, .cbor or .root produced by the mapping step
    # works here. Gen3 is JSON/CBOR-only -- see the Gen3 note right after this.
    material_map = getOpenDataDetectorDirectory() / "data/odd-material-maps.root"

    decorator = acts.IMaterialDecorator.fromFile(material_map)
    detector = getOpenDataDetector(materialDecorator=decorator, gen3=False)
    trackingGeometry = detector.trackingGeometry()
    #! [Load Material Map]

    #! [Extract Material Surfaces]
    # The surfaces the mapping writes onto, and reads back when mapping again.
    materialSurfaces = trackingGeometry.extractMaterialSurfaces()
    #! [Extract Material Surfaces]

    assert len(materialSurfaces) > 0
