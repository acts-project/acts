#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union

import acts
import acts.examples


u = acts.UnitConstants


def addSeeding(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    geoSelectionConfigFile: Union[Path, str],
    outputDir: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> acts.examples.Sequencer:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    geoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection
    outputDir : Path|str, path, None
        the output folder for the Root output, None triggers no output
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    """

    def customLogLevel(custom: acts.logging.Level = acts.logging.INFO):
        """override logging level"""
        if logLevel is None:
            return s.config.logLevel
        return acts.logging.Level(max(custom.value, logLevel.value))

    selAlg = acts.examples.TruthSeedSelector(
        level=customLogLevel(),
        ptMin=1 * u.GeV,
        eta=(-2.5, 2.5),
        nHitsMin=9,
        inputParticles="particles_final",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputParticles="truth_seeds_selected",
    )

    spAlg = acts.examples.SpacePointMaker(
        level=customLogLevel(),
        inputSourceLinks="sourcelinks",
        inputMeasurements="measurements",
        outputSpacePoints="spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
    )

    gridConfig = acts.SpacePointGridConfig(
        bFieldInZ=1.99724 * u.T,
        minPt=500 * u.MeV,
        rMax=200 * u.mm,
        zMax=2000 * u.mm,
        zMin=-2000 * u.mm,
        deltaRMax=60 * u.mm,
        cotThetaMax=7.40627,  # 2.7 eta
    )

    seedFilterConfig = acts.SeedFilterConfig(maxSeedsPerSpM=1, deltaRMin=1 * u.mm)

    seedFinderConfig = acts.SeedfinderConfig(
        rMax=gridConfig.rMax,
        deltaRMin=seedFilterConfig.deltaRMin,
        deltaRMax=gridConfig.deltaRMax,
        deltaRMinTopSP=seedFilterConfig.deltaRMin,
        deltaRMinBottomSP=seedFilterConfig.deltaRMin,
        deltaRMaxTopSP=gridConfig.deltaRMax,
        deltaRMaxBottomSP=gridConfig.deltaRMax,
        collisionRegionMin=-250 * u.mm,
        collisionRegionMax=250 * u.mm,
        zMin=gridConfig.zMin,
        zMax=gridConfig.zMax,
        maxSeedsPerSpM=seedFilterConfig.maxSeedsPerSpM,
        cotThetaMax=gridConfig.cotThetaMax,
        sigmaScattering=50,
        radLengthPerSeed=0.1,
        minPt=gridConfig.minPt,
        bFieldInZ=gridConfig.bFieldInZ,
        beamPos=acts.Vector2(0 * u.mm, 0 * u.mm),
        impactMax=3 * u.mm,
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=customLogLevel(acts.logging.VERBOSE),
        inputSpacePoints=[spAlg.config.outputSpacePoints],
        outputSeeds="seeds",
        outputProtoTracks="prototracks",
        gridConfig=gridConfig,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
    )

    parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
        level=customLogLevel(acts.logging.VERBOSE),
        inputProtoTracks=seedingAlg.config.outputProtoTracks,
        inputSpacePoints=[spAlg.config.outputSpacePoints],
        inputSourceLinks=spAlg.config.inputSourceLinks,
        outputTrackParameters="estimatedparameters",
        outputProtoTracks="prototracks_estimated",
        trackingGeometry=trackingGeometry,
        magneticField=field,
    )

    s.addAlgorithm(selAlg)
    s.addAlgorithm(spAlg)
    s.addAlgorithm(seedingAlg)
    s.addAlgorithm(parEstimateAlg)

    if outputDir is not None:
        s.addWriter(
            acts.examples.TrackFinderPerformanceWriter(
                level=customLogLevel(),
                inputProtoTracks=seedingAlg.config.outputProtoTracks,
                inputParticles=selAlg.config.outputParticles,
                inputMeasurementParticlesMap=selAlg.config.inputMeasurementParticlesMap,
                filePath=str(outputDir / "performance_seeding_trees.root"),
            )
        )

        s.addWriter(
            acts.examples.SeedingPerformanceWriter(
                level=customLogLevel(acts.logging.DEBUG),
                inputProtoTracks=seedingAlg.config.outputProtoTracks,
                inputParticles=selAlg.config.outputParticles,
                inputMeasurementParticlesMap=selAlg.config.inputMeasurementParticlesMap,
                filePath=str(outputDir / "performance_seeding_hists.root"),
            )
        )

        s.addWriter(
            acts.examples.RootTrackParameterWriter(
                level=customLogLevel(acts.logging.VERBOSE),
                inputTrackParameters=parEstimateAlg.config.outputTrackParameters,
                inputProtoTracks=parEstimateAlg.config.outputProtoTracks,
                inputParticles=selAlg.config.inputParticles,
                inputSimHits="simhits",
                inputMeasurementParticlesMap=selAlg.config.inputMeasurementParticlesMap,
                inputMeasurementSimHitsMap="measurement_simhits_map",
                filePath=str(outputDir / "estimatedparams.root"),
                treeName="estimatedparams",
            )
        )

    return s


def runSeeding(trackingGeometry, field, outputDir, s=None):

    from particle_gun import addParticleGun, EtaConfig, PhiConfig, ParticleConfig
    from fatras import addFatras
    from digitization import addDigitization

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )
    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    s = addParticleGun(
        s,
        EtaConfig(-2.0, 2.0),
        ParticleConfig(4, acts.PdgParticle.eMuon, True),
        PhiConfig(0.0, 360.0 * u.degree),
        multiplicity=2,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
    )

    s = addFatras(
        s,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
        preselectParticles=False,
    )

    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    s = addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        rnd=rnd,
    )

    s = addSeeding(
        s,
        trackingGeometry,
        field,
        geoSelectionConfigFile=srcdir
        / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json",
        outputDir=outputDir,
        logLevel=acts.logging.VERBOSE,
    )
    return s


if "__main__" == __name__:
    from common import getOpenDataDetector

    # detector, trackingGeometry, _ = getOpenDataDetector()
    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runSeeding(trackingGeometry, field, outputDir=Path.cwd()).run()
