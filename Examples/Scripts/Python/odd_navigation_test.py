#!/usr/bin/env python3
import os, pathlib, acts, acts.examples
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants
geoDir = getOpenDataDetectorDirectory()
outputDir = pathlib.Path.cwd() / "odd_output"

oddMaterialMap = geoDir / "data/odd-material-maps.root"
oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = geoDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector, trackingGeometry, decorators = getOpenDataDetector(
    geoDir, mdecorator=oddMaterialDeco
)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=1,
    numThreads=1,
    outputDir=str(outputDir),
)

s.addAlgorithm(
    acts.examples.NavigationTestAlgorithm(
        trackingGeometry=trackingGeometry,
        randomNumberSvc=rnd,
        magneticField=field,
        ntests=1,
        level=acts.logging.VERBOSE,
    )
)


s.run()
