#!/usr/bin/env python3
import argparse
import pathlib, acts, acts.examples
import acts.examples.dd4hep
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector

# acts.examples.dump_args_calls(locals())  # show python binding calls

u = acts.UnitConstants
outputDir = pathlib.Path.cwd() / "odd_output"
outputDir.mkdir(exist_ok=True)

oddDir = getOpenDataDetectorDirectory()

oddMaterialMap = oddDir / "data/odd-material-maps.root"
oddDigiConfig = oddDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = oddDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector, trackingGeometry, decorators = getOpenDataDetector(
    getOpenDataDetectorDirectory(), mdecorator=oddMaterialDeco
)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

from acts.examples.simulation import (
    addParticleGun,
    addPythia8,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    addCKFTracks,
    CKFPerformanceConfig,
    addVertexFitting,
    VertexFinder,
)

"""
s = acts.examples.Sequencer(events=100, skip=1, numThreads=1, logLevel=acts.logging.INFO)
s = addPythia8(
    s,
    npileup=200,
    rnd=rnd,
    outputDirCsv=outputDir,
)
"""
s = acts.examples.Sequencer(events=1, numThreads=1, logLevel=acts.logging.VERBOSE)
s.addReader(
    acts.examples.CsvParticleReader(
        level=acts.logging.INFO,
        inputDir=str(outputDir),
        inputStem="myparticles",
        outputParticles="particles_input",
    )
)
s = addFatras(
    s,
    trackingGeometry,
    field,
    outputDirRoot=outputDir,
    rnd=rnd,
)
"""
s = addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=oddDigiConfig,
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addSeeding(
    s,
    trackingGeometry,
    field,
    geoSelectionConfigFile=oddSeedingSel,
    outputDirRoot=outputDir,
)
s = addCKFTracks(
    s,
    trackingGeometry,
    field,
    CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
    outputDirRoot=outputDir,
)
"""

s.run()
