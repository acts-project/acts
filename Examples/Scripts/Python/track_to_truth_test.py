#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

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
    TruthJetConfig,
    TrackToTruthJetConfig,
    addFatras,
    addPythia8,
    addTruthJetAlg,
    addTrackToTruthJetAlg,
    addDigitization,
    ParticleSelectorConfig,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    addKalmanTracks,
)

s = acts.examples.Sequencer(events=1, numThreads=-1, logLevel=acts.logging.INFO)
outputDir = Path.cwd() / " trackToTruth_output"
outputDir.mkdir(exist_ok=True)

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
    nhard=1,
    npileup=1,
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

s.addWriter(
    acts.examples.RootTrackStatesWriter(
        level=acts.logging.INFO,
        inputTracks="tracks",
        inputParticles="particles_selected",
        inputTrackParticleMatching="track_particle_matching",
        inputSimHits="simhits",
        inputMeasurementSimHitsMap="measurement_simhits_map",
        filePath=str(outputDir / "trackstates_kf.root"),
    )
)

s.addWriter(
    acts.examples.RootTrackSummaryWriter(
        level=acts.logging.INFO,
        inputTracks="tracks",
        inputParticles="particles_selected",
        inputTrackParticleMatching="track_particle_matching",
        filePath=str(outputDir / "tracksummary_kf.root"),
    )
)

s.addWriter(
    acts.examples.TrackFitterPerformanceWriter(
        level=acts.logging.INFO,
        inputTracks="tracks",
        inputParticles="particles_selected",
        inputTrackParticleMatching="track_particle_matching",
        filePath=str(outputDir / "performance_kf.root"),
    )
)

addTruthJetAlg(
    s,
    TruthJetConfig(
        inputTruthParticles="particles_generated_selected",
        outputJets="truth_jets",
        jetPtMin=10 * u.GeV,
        inputHepMC3Event="pythia8-event",
        doJetLabeling=True,
    ),
    loglevel=acts.logging.VERBOSE,
)

addTrackToTruthJetAlg(
    s,
    TrackToTruthJetConfig(
        inputTracks="tracks",
        inputJets="truth_jets",
        outputTrackJets="track_jets",
        maxDeltaR=0.4,
    ),
    loglevel=acts.logging.DEBUG,
)

s.run()
