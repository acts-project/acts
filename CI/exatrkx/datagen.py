#!/usr/bin/env python3

import os
import sys
from pathlib import Path

import acts
import acts.examples
from acts.examples.reconstruction import *
from acts.examples.simulation import *

u = acts.UnitConstants

#######################
# Handle command line #
#######################

def usage():
    print("Usage: {} <target_path> <n_events>".format(sys.argv[0]))

try:
    numberEvents=int(sys.argv[2])
except:
    usage()
    exit(1)

try:
    outputDir = Path(sys.argv[1])
    (outputDir / "train_all").mkdir(parents=True, exist_ok=True)
except:
    usage()
    exit(1)


############################################
# Prepare detector, bfield, random numbers #
############################################

detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()
digiConfig = Path(__file__).parents[1] / "../Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
geoSelection = Path(__file__).parents[1] / "../Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"

assert digiConfig.exists()
assert geoSelection.exists()

field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
rnd = acts.examples.RandomNumbers(seed=1234)


#############################
# Prepare and run sequencer #
#############################

s = acts.examples.Sequencer(
    events=numberEvents,
    numThreads=-1,
    outputDir=str(outputDir),
)

s = addParticleGun(
    s,
    MomentumConfig(1.0 * u.GeV, 5.0 * u.GeV, True),
    EtaConfig(-3.0, 3.0, True),
    ParticleConfig(10, acts.PdgParticle.eMuon, True),
    rnd=rnd,
    multiplicity=50,
    outputDirCsv=str(outputDir/"train_all"),
)

s.addAlgorithm(
    acts.examples.FatrasSimulation(
        level=s.config.logLevel,
        inputParticles="particles_input",
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
    )
)

s = addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=digiConfig,
    outputDirRoot=None,
    outputDirCsv=str(outputDir/"train_all"),
    rnd=rnd,
)

s.addWriter(
    acts.examples.CsvSimHitWriter(
        level=acts.logging.INFO,
        inputSimHits="simhits",
        outputDir=str(outputDir/"train_all"),
        outputStem="truth",
    )
)

s.addWriter(
    acts.examples.CsvMeasurementWriter(
        level=acts.logging.INFO,
        inputMeasurements="measurements",
        inputClusters="clusters",
        inputSimHits="simhits",
        inputMeasurementSimHitsMap="measurement_simhits_map",
        outputDir=str(outputDir/"train_all"),
    )
)

s.addWriter(
    acts.examples.CsvTrackingGeometryWriter(
        level=acts.logging.INFO,
        trackingGeometry=trackingGeometry,
        outputDir=str(outputDir),
        writePerEvent=False,
    )
)
        
s.run()
