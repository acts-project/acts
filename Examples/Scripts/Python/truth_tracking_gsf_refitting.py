#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

import acts
import acts.examples

u = acts.UnitConstants


def runRefittingGsf(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    inputParticlePath: Optional[Path] = None,
    inputSimHitsPath: Optional[Path] = None,
    decorators=[],
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
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        TrackSmearingSigmas,
        addKalmanTracks,
    )

    # Initialize sequencer
    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    # Particle input: read from file or generate electrons
    if inputParticlePath is None:
        addParticleGun(
            s,
            ParticleConfig(num=1, pdg=acts.PdgParticle.eElectron, randomizeCharge=True),
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
        acts.logging.getLogger("GSF Refitting").info(
            "Reading particles from %s", inputParticlePath.resolve()
        )
        assert inputParticlePath.exists()
        s.addReader(
            acts.examples.RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                outputParticles="particles_generated",
            )
        )

    # SimHits input: read from file or run FATRAS
    if inputSimHitsPath is None:
        addFatras(
            s,
            trackingGeometry,
            field,
            rnd=rnd,
            enableInteractions=True,
        )
    else:
        s.addReader(
            acts.examples.RootSimHitReader(
                level=acts.logging.INFO,
                filePath=str(inputSimHitsPath.resolve()),
                outputSimHits="simhits",
            )
        )

    # Digitization
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    # Particle selection (electrons, not muons)
    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(0.9 * u.GeV, None),
            measurements=(7, None),
            removeNeutral=True,
            removeSecondaries=True,
        ),
    )

    # Seeding with electron hypothesis (copied from truth_tracking_gsf.py)
    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_generated",
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        trackSmearingSigmas=TrackSmearingSigmas(
            # zero everything so the GSF has a chance to find the measurements
            loc0=0,
            loc0PtA=0,
            loc0PtB=0,
            loc1=0,
            loc1PtA=0,
            loc1PtB=0,
            time=0,
            phi=0,
            theta=0,
            ptRel=0,
        ),
        particleHypothesis=acts.ParticleHypothesis.electron,
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0 / u.GeV,
            1 * u.ns,
        ],
        initialSigmaQoverPt=0.1 / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1e0, 1e0, 1e0, 1e0, 1e0, 1e0],
    )

    # Initial Kalman fit (produces tracks to refit with GSF)
    addKalmanTracks(
        s,
        trackingGeometry,
        field,
        directNavigation=True,
        reverseFilteringMomThreshold=0.0,
    )

    # Track selection: select good Kalman tracks
    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="kf_tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=7,
            ),
        )
    )

    # NOTE we specify clampToRange as True to silence warnings in the test about
    # queries to the loss distribution outside the specified range, since no dedicated
    # approximation for the ODD is done yet.
    bha = acts.examples.AtlasBetheHeitlerApprox.makeDefault(clampToRange=True)

    gsfOptions = {
        "betheHeitlerApprox": bha,
        "maxComponents": 12,
        "componentMergeMethod": acts.examples.ComponentMergeMethod.maxWeight,
        "mixtureReductionAlgorithm": acts.examples.MixtureReductionAlgorithm.KLDistance,
        "weightCutoff": 1.0e-4,
        "reverseFilteringCovarianceScaling": 100.0,
        "level": acts.logging.INFO,
    }

    s.addAlgorithm(
        acts.examples.RefittingAlgorithm(
            acts.logging.INFO,
            inputTracks="kf_tracks",
            outputTracks="gsf_refit_tracks",
            initialVarInflation=6 * [100.0],
            fit=acts.examples.makeGsfFitterFunction(
                trackingGeometry, field, **gsfOptions
            ),
        )
    )

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=acts.logging.INFO,
            inputTracks="gsf_refit_tracks",
            inputParticles="particles_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="refit_track_particle_matching",
            outputParticleTrackMatching="refit_particle_track_matching",
        )
    )

    s.addWriter(
        acts.examples.RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="gsf_refit_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_gsf_refit.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="gsf_refit_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            filePath=str(outputDir / "tracksummary_gsf_refit.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="gsf_refit_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            filePath=str(outputDir / "performance_gsf_refit.root"),
        )
    )

    return s


if __name__ == "__main__":
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # ODD
    from acts.examples.odd import getOpenDataDetector

    detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()
    digiConfigFile = srcdir / "Examples/Configs/odd-digi-smearing-config.json"

    ## GenericDetector
    # detector = acts.examples.GenericDetector()
    # trackingGeometry = detector.trackingGeometry()
    # digiConfigFile = (
    #     srcdir
    #     / "Examples/Configs/generic-digi-smearing-config.json"
    # )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runRefittingGsf(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        outputDir=Path.cwd(),
    ).run()
