#!/usr/bin/env python3

import pathlib, acts, acts.examples
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    addCKFTracks,
    addVertexFitting,
    VertexFinder,
)
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector


u = acts.UnitConstants
geoDir = getOpenDataDetectorDirectory()
outputDir = pathlib.Path.cwd() / "tutorial_output"

oddMaterialMap = geoDir / "data/odd-material-maps.root"
oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = geoDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector, trackingGeometry, decorators = getOpenDataDetector(
    geoDir, mdecorator=oddMaterialDeco
)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(events=100, numThreads=1, outputDir=str(outputDir))

addParticleGun(
    s,
    MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
    EtaConfig(-3.0, 3.0, uniform=True),
    ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
    vtxGen=acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
        mean=acts.Vector4(0, 0, 0, 0),
    ),
    multiplicity=5,
    rnd=rnd,
)

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    outputDirRoot=outputDir,
    outputDirCsv=outputDir,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=oddDigiConfig,
    rnd=rnd,
    outputDirRoot=outputDir,
    outputDirCsv=outputDir,
)

addSeeding(
    s,
    trackingGeometry,
    field,
    geoSelectionConfigFile=oddSeedingSel,
    outputDirRoot=outputDir,
)

"""
# tutorial algorithm
s.addAlgorithm(
    acts.examples.TutorialAlgorithm(
        level=acts.logging.VERBOSE,
        message="hello world from python!",
    )
)
"""

addCKFTracks(
    s,
    trackingGeometry,
    field,
    outputDirRoot=outputDir,
    outputDirCsv=outputDir,
)

addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.Iterative,
    outputDirRoot=outputDir,
)

s.run()
