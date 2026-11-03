def test_load_material_map():
    #! [Load Material Map]
    import acts
    from acts.examples.odd import getOpenDataDetectorDirectory, getOpenDataDetector

    # This bundled map uses ROOT; versioned JSON maps can instead be loaded
    # with TrackingGeometryMaterialJsonConverter and applied to the geometry.
    material_map = getOpenDataDetectorDirectory() / "data/odd-material-maps.root"
    decorator = acts.root.RootMaterialDecorator(
        fileName=str(material_map), level=acts.logging.INFO
    )
    detector = getOpenDataDetector(materialDecorator=decorator)
    trackingGeometry = detector.trackingGeometry()
    #! [Load Material Map]

    #! [Extract Material Surfaces]
    # The surfaces the mapping writes onto, and reads back when mapping again.
    materialSurfaces = trackingGeometry.extractMaterialSurfaces()
    #! [Extract Material Surfaces]

    assert len(materialSurfaces) > 0
