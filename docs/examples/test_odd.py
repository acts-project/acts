def test_basic():
    #! [Basic ODD construction]
    from acts.examples.odd import getOpenDataDetector

    # No material decorator: the official ODD repo doesn't yet ship a Gen3
    # material map (data/odd-material-maps.root is Gen1-keyed and
    # incompatible with Gen3 surfaces). See the material mapping how-to for
    # producing your own.
    detector = getOpenDataDetector(gen3=True)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    #! [Basic ODD construction]
