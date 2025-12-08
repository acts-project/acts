def test_basic():
    #! [Basic ODD construction]
    import acts.root
    from acts.examples.odd import getOpenDataDetectorDirectory, getOpenDataDetector

    odd_dir = getOpenDataDetectorDirectory()
    materialDecorator = acts.root.RootMaterialDecorator(
        fileName=str(odd_dir / "data/odd-material-maps.root"),
        level=acts.logging.INFO,
    )
    detector = getOpenDataDetector(materialDecorator=materialDecorator)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    #! [Basic ODD construction]
