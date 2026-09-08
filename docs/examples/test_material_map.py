def test_load_material_map():
    #! [Load Material Map]
    import acts
    from acts.examples.odd import getOpenDataDetectorDirectory, getOpenDataDetector

    # Any of .json, .cbor or .root produced by the mapping step works here.
    material_map = getOpenDataDetectorDirectory() / "data/odd-material-maps.root"

    decorator = acts.IMaterialDecorator.fromFile(material_map)
    detector = getOpenDataDetector(materialDecorator=decorator)
    trackingGeometry = detector.trackingGeometry()
    #! [Load Material Map]

    #! [Extract Material Surfaces]
    # The surfaces the mapping writes onto, and reads back when mapping again.
    materialSurfaces = trackingGeometry.extractMaterialSurfaces()
    #! [Extract Material Surfaces]

    assert len(materialSurfaces) > 0
