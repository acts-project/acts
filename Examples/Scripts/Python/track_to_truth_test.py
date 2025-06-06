#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import argparse

import acts
import acts.examples

actsDir = Path(__file__).parent.parent.parent.parent

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

from acts.examples.simulation import (
    addParticleGun,
    ParticleConfig,
    EtaConfig,
    PhiConfig,
    MomentumConfig,
    TrackToTruthJetConfig,
    addFatras,
    addPythia8,
    addTrackToTruthJetAlg,
    addDigitization,
    ParticleSelectorConfig,
    addDigiParticleSelection,
    addGenParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    addKalmanTracks,
)

parser = argparse.ArgumentParser()
parser.add_argument("--events", "-n", type=int, default=1000)
parser.add_argument("--pileup", "--pu", "-p", type=int, default=0)
parser.add_argument("--hardscatter", "--hs", type=int, default=1)
parser.add_argument("--jobs", "-j", type=int, default=-1)
parser.add_argument("--csv", action="store_true")
args = parser.parse_args()

outputDir = Path.cwd() / "trackToTruth_output"
print(outputDir)
outputDir.mkdir(exist_ok=True)
s = acts.examples.Sequencer(
    events=args.events,
    numThreads=args.jobs,
    logLevel=acts.logging.INFO,
    outputDir=str(outputDir),
)

from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

geoDir = getOpenDataDetectorDirectory()
# acts.examples.dump_args_calls(locals())  # show python binding calls

oddMaterialMap = geoDir / "data/odd-material-maps.root"
assert oddMaterialMap.exists(), f"Material map file {oddMaterialMap} does not exist"

oddDigiConfig = actsDir / "Examples/Configs/odd-digi-smearing-config.json"
assert oddDigiConfig.exists(), f"Digi config file {oddDigiConfig} does not exist"

oddSeedingSel = actsDir / "Examples/Configs/odd-seeding-config.json"
assert oddSeedingSel.exists(), f"Seeding config file {oddSeedingSel} does not exist"

oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector = getOpenDataDetector(
    odd_dir=geoDir,
    materialDecorator=oddMaterialDeco,
)
trackingGeometry = detector.trackingGeometry()

field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


rnd = acts.examples.RandomNumbers(seed=42)

addPythia8(
    s,
    nhard=args.hardscatter,
    npileup=args.pileup,
    hardProcess=["Top:qqbar2ttbar=on"],
    vtxGen=acts.examples.GaussianVertexGenerator(
        mean=acts.Vector4(0, 0, 0, 0),
        stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
    ),
    rnd=rnd,
    outputDirRoot=None,
    outputDirCsv=None,
    writeHepMC3=None,
)

# Effective truth level selection for simulation + track reconstruction
addGenParticleSelection(
    s,
    ParticleSelectorConfig(
        rho=(0.0, 24 * u.mm),
        absZ=(0.0, 1.0 * u.m),
        eta=(-3.0, 3.0),
        pt=(150 * u.MeV, None),
    ),
)

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    enableInteractions=True,
)


addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=oddDigiConfig,
    rnd=rnd,
    logLevel=acts.logging.ERROR,
)

addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        pt=(0.9 * u.GeV, None),
        measurements=(7, None),
        removeNeutral=True,
        removeSecondaries=True,
    ),
)

addSeeding(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    inputParticles="particles_generated",
    seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
    particleHypothesis=acts.ParticleHypothesis.muon,
)

reverseFilteringMomThreshold = 0 * u.GeV

addKalmanTracks(
    s,
    trackingGeometry,
    field,
    reverseFilteringMomThreshold,
    logLevel=acts.logging.FATAL,
)

s.addAlgorithm(
    acts.examples.TrackSelectorAlgorithm(
        level=acts.logging.INFO,
        inputTracks="tracks",
        outputTracks="selected-tracks",
        selectorConfig=acts.TrackSelector.Config(
            minMeasurements=7,
        ),
    )
)
s.addWhiteboardAlias("tracks", "selected-tracks")

# s.addWriter(
#     acts.examples.RootTrackStatesWriter(
#         level=acts.logging.INFO,
#         inputTracks="tracks",
#         inputParticles="particles_selected",
#         inputTrackParticleMatching="track_particle_matching",
#         inputSimHits="simhits",
#         inputMeasurementSimHitsMap="measurement_simhits_map",
#         filePath=str(outputDir / "trackstates_kf.root"),
#     )
# )

s.addWriter(
    acts.examples.RootTrackSummaryWriter(
        level=acts.logging.FATAL,
        inputTracks="tracks",
        inputParticles="particles_selected",
        inputTrackParticleMatching="track_particle_matching",
        filePath=str(outputDir / "tracksummary_kf.root"),
    )
)

s.addWriter(
    acts.examples.TrackFitterPerformanceWriter(
        level=acts.logging.FATAL,
        inputTracks="tracks",
        inputParticles="particles_selected",
        inputTrackParticleMatching="track_particle_matching",
        filePath=str(outputDir / "performance_kf.root"),
    )
)

s.addAlgorithm(
    acts.examples.ParticleSelector(
        level=acts.logging.INFO,
        inputParticles="particles_generated",
        outputParticles="jet_input_particles",
        # eta=(-3.0, 3.0),
        pt=(150 * u.MeV, None),
    )
)

truthJetAlg = acts.examples.TruthJetAlgorithm(
    level=acts.logging.INFO,
    inputTruthParticles="jet_input_particles",
    outputJets="truth_jets",
    jetPtMin=10 * u.GeV,
    inputHepMC3Event="pythia8-event",
    doJetLabeling=True,
    jetLabelingHadronPtMin=5 * u.GeV,
    # if we don't have hard scatter, use all particles, else only use hard scatter particles
    jetLabelingHardScatterHadronsOnly=args.hardscatter != 0,
    clusterHardScatterParticlesOnly=args.hardscatter != 0,
    debugCsvOutput=args.csv,
)

s.addAlgorithm(truthJetAlg)

# addTrackToTruthJetAlg(
#     s,
#     TrackToTruthJetConfig(
#         inputTracks="tracks",
#         inputJets="truth_jets",
#         outputTrackJets="track_jets",
#         maxDeltaR=0.4,
#     ),
#     loglevel=acts.logging.DEBUG,
# )

s.run()
