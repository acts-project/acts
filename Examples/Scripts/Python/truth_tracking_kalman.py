#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

import acts
import acts.examples

u = acts.UnitConstants


def runTruthTrackingKalman(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    inputParticlePath: Optional[Path] = None,
    decorators=[],
    directNavigation=False,
    reverseFilteringMomThreshold=0 * u.GeV,
    s: acts.examples.Sequencer = None,
):
    from acts.examples.simulation import (
        addParticleGun,
        ParticleConfig,
        EtaConfig,
        PhiConfig,
        MomentumConfig,
        addFatras,
        addDigitization,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        TruthSeedRanges,
        addKalmanTracks,
    )

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    if inputParticlePath is None:
        addParticleGun(
            s,
            ParticleConfig(num=1, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
            EtaConfig(-3.0, 3.0, uniform=True),
            MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
            PhiConfig(0.0, 360.0 * u.degree),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=1,
            rnd=rnd,
        )
    else:
        acts.logging.getLogger("Truth tracking example").info(
            "Reading particles from %s", inputParticlePath.resolve()
        )
        assert inputParticlePath.exists()
        s.addReader(
            acts.examples.RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                outputParticles="particles_input",
            )
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
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_input",
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        particleHypothesis=acts.ParticleHypothesis.muon,
        truthSeedRanges=TruthSeedRanges(
            nHits=(7, None),
        ),
    )

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
        directNavigation,
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
            inputParticles="truth_seeds_selected",
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
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf.root"),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_kf.root"),
        )
    )

    return s


if "__main__" == __name__:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # ODD
    from acts.examples.odd import getOpenDataDetector

    detector, trackingGeometry, decorators = getOpenDataDetector()
    digiConfigFile = (
        srcdir / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json"
    )

    ## GenericDetector
    # detector, trackingGeometry, _ = acts.examples.GenericDetector.create()
    # digiConfigFile = (
    #     srcdir
    #     / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    # )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runTruthTrackingKalman(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        outputDir=Path.cwd(),
    ).run()
