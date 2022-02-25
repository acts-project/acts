#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union
from enum import Enum

import acts
import acts.examples


u = acts.UnitConstants
SeedingAlgorithm = Enum("SeedingAlgorithm", "Default TruthSmeared TruthEstimated")


def addSeeding(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    geoSelectionConfigFile: Optional[Union[Path, str]] = None,
    seedingAlgorithm: SeedingAlgorithm = SeedingAlgorithm.Default,
    outputDir: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    rnd: Optional[acts.examples.RandomNumbers] = None,
) -> acts.examples.Sequencer:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    geoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection. Not required for SeedingAlgorithm.TruthSmeared.
    outputDir : Path|str, path, None
        the output folder for the Root output, None triggers no output
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    rnd : RandomNumbers, None
        random number generator. Only used by SeedingAlgorithm.TruthSmeared.
    """

    def customLogLevel(custom: acts.logging.Level = acts.logging.INFO):
        """override logging level"""
        if logLevel is None:
            return s.config.logLevel
        return acts.logging.Level(max(custom.value, logLevel.value))

    logger = acts.logging.getLogger("addSeeding")

    selAlg = acts.examples.TruthSeedSelector(
        level=customLogLevel(),
        ptMin=1 * u.GeV,
        eta=(-2.5, 2.5),
        nHitsMin=9,
        inputParticles="particles_final",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputParticles="truth_seeds_selected",
    )
    s.addAlgorithm(selAlg)
    inputParticles = selAlg.config.outputParticles

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if seedingAlgorithm == SeedingAlgorithm.TruthSmeared:
        rnd = rnd or acts.examples.RandomNumbers(seed=42)
        logger.info("Using smeared truth particles for seeding")
        # Run particle smearing
        ptclSmear = acts.examples.ParticleSmearing(
            level=customLogLevel(),
            inputParticles=inputParticles,
            outputTrackParameters="smearedparameters",
            randomNumbers=rnd,
            # gaussian sigmas to smear particle parameters
            sigmaD0=20 * u.um,
            sigmaD0PtA=30 * u.um,
            sigmaD0PtB=0.3 / 1 * u.GeV,
            sigmaZ0=20 * u.um,
            sigmaZ0PtA=30 * u.um,
            sigmaZ0PtB=0.3 / 1 * u.GeV,
            sigmaPhi=1 * u.degree,
            sigmaTheta=1 * u.degree,
            sigmaPRel=0.01,
            sigmaT0=1 * u.ns,
            initialVarInflation=[1, 1, 1, 1, 1, 1],
        )
        s.addAlgorithm(ptclSmear)
    else:

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
        s.addAlgorithm(spAlg)

        # Run either: truth track finding or seeding
        if seedingAlgorithm == SeedingAlgorithm.TruthEstimated:
            logger.info("Using truth track finding from space points for seeding")
            # Use truth tracking
            truthTrackFinder = acts.examples.TruthTrackFinder(
                level=customLogLevel(),
                inputParticles=inputParticles,
                inputMeasurementParticlesMap=selAlg.config.inputMeasurementParticlesMap,
                outputProtoTracks="prototracks",
            )
            s.addAlgorithm(truthTrackFinder)
            inputProtoTracks = truthTrackFinder.config.outputProtoTracks
            inputSeeds = ""
        elif seedingAlgorithm == SeedingAlgorithm.Default:
            logger.info("Using default seeding")
            # Use seeding
            gridConfig = acts.SpacePointGridConfig(
                bFieldInZ=1.99724 * u.T,
                minPt=500 * u.MeV,
                rMax=200 * u.mm,
                zMax=2000 * u.mm,
                zMin=-2000 * u.mm,
                deltaRMax=60 * u.mm,
                cotThetaMax=7.40627,  # 2.7 eta
            )

            seedFilterConfig = acts.SeedFilterConfig(
                maxSeedsPerSpM=1, deltaRMin=1 * u.mm
            )

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
            s.addAlgorithm(seedingAlg)
            inputProtoTracks = seedingAlg.config.outputProtoTracks
            inputSeeds = seedingAlg.config.outputSeeds
        else:
            logger.fatal("unknown seedingAlgorithm", seedingAlgorithm)

        parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
            level=customLogLevel(acts.logging.VERBOSE),
            inputSeeds=inputSeeds,
            inputProtoTracks=inputProtoTracks,
            inputSpacePoints=[spAlg.config.outputSpacePoints],
            inputSourceLinks=spAlg.config.inputSourceLinks,
            outputTrackParameters="estimatedparameters",
            outputProtoTracks="prototracks_estimated",
            trackingGeometry=trackingGeometry,
            magneticField=field,
        )
        s.addAlgorithm(parEstimateAlg)

        if outputDir is not None:
            s.addWriter(
                acts.examples.TrackFinderPerformanceWriter(
                    level=customLogLevel(),
                    inputProtoTracks=inputProtoTracks,
                    inputParticles=inputParticles,  # the original selected particles after digitization
                    inputMeasurementParticlesMap=selAlg.config.inputMeasurementParticlesMap,
                    filePath=str(outputDir / "performance_seeding_trees.root"),
                )
            )

            s.addWriter(
                acts.examples.SeedingPerformanceWriter(
                    level=customLogLevel(acts.logging.DEBUG),
                    inputProtoTracks=inputProtoTracks,
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
        # seedingAlgorithm=SeedingAlgorithm.TruthEstimated,
    )
    return s


if "__main__" == __name__:
    from common import getOpenDataDetector

    # detector, trackingGeometry, _ = getOpenDataDetector()
    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runSeeding(trackingGeometry, field, outputDir=Path.cwd()).run()
