#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

import acts
import acts.examples

u = acts.UnitConstants


def runTruthTrackingGx2f(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    outputDir: Path,
    digiConfigFile: Path,
    s: acts.examples.Sequencer = None,
    inputParticlePath: Optional[Path] = None,
):
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
        SeedingAlgorithm,
        TruthSeedRanges,
        addGx2fTracks,
    )

    s = s or acts.examples.Sequencer(
        events=10000, numThreads=-1, logLevel=acts.logging.INFO
    )

    rnd = acts.examples.RandomNumbers()
    outputDir = Path(outputDir)

    if inputParticlePath is None:
        addParticleGun(
            s,
            MomentumConfig(100.0 * u.GeV, 100.0 * u.GeV, transverse=True),
            EtaConfig(-2.0, 2.0),
            ParticleConfig(2, acts.PdgParticle.eMuon, False),
            multiplicity=1,
            rnd=rnd,
            outputDirRoot=outputDir,
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
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        rnd=rnd,
        truthSeedRanges=TruthSeedRanges(
            pt=(1 * u.GeV, None),
            nHits=(9, None),
        ),
    )

    addGx2fTracks(
        s,
        trackingGeometry,
        field,
        nUpdateMax=17,
        relChi2changeCutOff=1e-7,
    )

    # Output
    s.addWriter(
        acts.examples.RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_gx2f.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_gx2f.root"),
            writeGx2fSpecific=True,
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_gx2f.root"),
        )
    )

    return s


if "__main__" == __name__:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # detector, trackingGeometry, _ = getOpenDataDetector()
    detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runTruthTrackingGx2f(
        trackingGeometry=trackingGeometry,
        # decorators=decorators,
        field=field,
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        # "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json",
        # outputCsv=True,
        # inputParticlePath=inputParticlePath,
        outputDir=Path.cwd(),
    ).run()
