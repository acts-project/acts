#!/usr/bin/env python3

from pathlib import Path
from typing import Optional, Union

from acts.examples import Sequencer, GenericDetector, RootParticleReader

import acts

u = acts.UnitConstants


def runTruthTrackingGsf(
    trackingGeometry,
    digiConfigFile: Path,
    field,
    outputDir: Path,
    outputCsv=True,
    inputParticlePath: Optional[Path] = None,
    decorators=[],
    s=None,
):

    from acts.examples.simulation import (
        addParticleGun,
        EtaConfig,
        PhiConfig,
        ParticleConfig,
        addFatras,
        addDigitization,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        addTruthTrackingGsf,
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
            EtaConfig(-2.0, 2.0),
            ParticleConfig(4, acts.PdgParticle.eElectron, True),
            PhiConfig(0.0, 360.0 * u.degree),
            multiplicity=2,
            rnd=rnd,
        )
    else:
        acts.logging.getLogger("GSF Example").info(
            "Reading particles from %s", inputParticlePath.resolve()
        )
        assert inputParticlePath.exists()
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection="particles_input",
                orderedEvents=False,
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
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
    )

    truthTrkFndAlg = acts.examples.TruthTrackFinder(
        level=acts.logging.INFO,
        inputParticles="truth_seeds_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputProtoTracks="prototracks",
    )

    s.addAlgorithm(truthTrkFndAlg)

    addTruthTrackingGsf(s, trackingGeometry, field)

    # Output
    s.addWriter(
        acts.examples.RootTrajectoryStatesWriter(
            level=acts.logging.INFO,
            inputTrajectories="gsf_trajectories",
            inputParticles="truth_seeds_selected",
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_gsf.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrajectorySummaryWriter(
            level=acts.logging.INFO,
            inputTrajectories="gsf_trajectories",
            inputParticles="truth_seeds_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            filePath=str(outputDir / "tracksummary_gsf.root"),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTrajectories="gsf_trajectories",
            inputParticles="truth_seeds_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            filePath=str(outputDir / "performance_gsf.root"),
        )
    )

    return s


if "__main__" == __name__:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    runTruthTrackingGsf(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        field=field,
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        outputCsv=True,
        inputParticlePath=inputParticlePath,
        outputDir=Path.cwd(),
    ).run()
