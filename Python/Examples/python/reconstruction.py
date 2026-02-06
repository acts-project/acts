from pathlib import Path
from typing import Optional, Union, List
from enum import Enum
from collections import namedtuple

import acts
import acts.examples

# ROOT might not be available
try:
    from acts.examples.root import (
        RootTrackFinderNTupleWriter,
        RootTrackFinderPerformanceWriter,
        RootTrackFitterPerformanceWriter,
        RootTrackParameterWriter,
        RootTrackStatesWriter,
        RootTrackSummaryWriter,
        RootVertexNTupleWriter,
    )

    ACTS_ROOT_AVAILABLE = True
except ImportError:
    ACTS_ROOT_AVAILABLE = False

u = acts.UnitConstants

SeedingAlgorithm = Enum(
    "SeedingAlgorithm",
    "Default TruthSmeared TruthEstimated Orthogonal HoughTransform AdaptiveHoughTransform Gbts Hashing GridTriplet OrthogonalTriplet",
)

TrackSmearingSigmas = namedtuple(
    "TrackSmearingSigmas",
    [
        "loc0",
        "loc0PtA",
        "loc0PtB",
        "loc1",
        "loc1PtA",
        "loc1PtB",
        "time",
        "phi",
        "theta",
        "ptRel",
    ],
    defaults=[None] * 10,
)

SeedFinderConfigArg = namedtuple(
    "SeedFinderConfig",
    [
        "maxSeedsPerSpM",
        "cotThetaMax",
        "sigmaScattering",
        "radLengthPerSeed",
        "minPt",
        "impactMax",
        "deltaPhiMax",
        "interactionPointCut",
        "deltaZMax",
        "maxPtScattering",
        "zBinEdges",
        "zBinsCustomLooping",
        "rRangeMiddleSP",
        "useVariableMiddleSPRange",
        "binSizeR",
        "seedConfirmation",
        "centralSeedConfirmationRange",
        "forwardSeedConfirmationRange",
        "deltaR",  # (min,max)
        "deltaRBottomSP",  # (min,max)
        "deltaRTopSP",  # (min,max)
        "deltaRMiddleSPRange",  # (min,max)
        "collisionRegion",  # (min,max)
        "r",  # (min,max)
        "z",  # (min,max)
    ],
    defaults=[None] * 18 + [(None, None)] * 7,
)
SeedFinderOptionsArg = namedtuple(
    "SeedFinderOptions", ["beamPos", "bFieldInZ"], defaults=[(None, None), None]
)

SeedFilterConfigArg = namedtuple(
    "SeedFilterConfig",
    [
        "impactWeightFactor",
        "zOriginWeightFactor",
        "compatSeedWeight",
        "compatSeedLimit",
        "numSeedIncrement",
        "seedWeightIncrement",
        "seedConfirmation",
        "maxSeedsPerSpMConf",
        "maxQualitySeedsPerSpMConf",
        "useDeltaRorTopRadius",
        "deltaRMin",
    ],
    defaults=[None] * 11,
)

SpacePointGridConfigArg = namedtuple(
    "SeedGridConfig",
    [
        "rMax",
        "zBinEdges",
        "phiBinDeflectionCoverage",
        "impactMax",
        "deltaRMax",
        "maxPhiBins",
        "phi",  # (min,max)
    ],
    defaults=[None] * 6 + [(None, None)] * 1,
)

SeedingAlgorithmConfigArg = namedtuple(
    "SeedingAlgorithmConfig",
    [
        "allowSeparateRMax",
        "zBinNeighborsTop",
        "zBinNeighborsBottom",
        "numPhiNeighbors",
        "useExtraCuts",
    ],
    defaults=[None] * 5,
)

TruthEstimatedSeedingAlgorithmConfigArg = namedtuple(
    "TruthSeederConfig",
    [
        "deltaR",  # (min,max)
    ],
    defaults=[(None, None)],
)

TrackSelectorConfig = namedtuple(
    "TrackSelectorConfig",
    [
        "loc0",
        "loc1",
        "time",
        "eta",
        "absEta",
        "pt",
        "phi",
        "nMeasurementsMin",
        "maxHoles",
        "maxOutliers",
        "maxHolesAndOutliers",
        "maxSharedHits",
        "maxChi2",
        "nMeasurementsGroupMin",
        "requireReferenceSurface",
    ],
    defaults=[(None, None)] * 7 + [None] * 8,
)


def trackSelectorDefaultKWArgs(c):
    """
    Encapsulate this boilerplate code into a function so different uses do not get out of sync
    """
    return acts.examples.defaultKWArgs(
        loc0Min=c.loc0[0],
        loc0Max=c.loc0[1],
        loc1Min=c.loc1[0],
        loc1Max=c.loc1[1],
        timeMin=c.time[0],
        timeMax=c.time[1],
        phiMin=c.phi[0],
        phiMax=c.phi[1],
        etaMin=c.eta[0],
        etaMax=c.eta[1],
        absEtaMin=c.absEta[0],
        absEtaMax=c.absEta[1],
        ptMin=c.pt[0],
        ptMax=c.pt[1],
        minMeasurements=c.nMeasurementsMin,
        maxHoles=c.maxHoles,
        maxOutliers=c.maxOutliers,
        maxHolesAndOutliers=c.maxHolesAndOutliers,
        maxSharedHits=c.maxSharedHits,
        maxChi2=c.maxChi2,
        measurementCounter=c.nMeasurementsGroupMin,
        requireReferenceSurface=c.requireReferenceSurface,
    )


CkfConfig = namedtuple(
    "CkfConfig",
    [
        "chi2CutOffMeasurement",
        "chi2CutOffOutlier",
        "numMeasurementsCutOff",
        "maxSteps",
        "seedDeduplication",
        "stayOnSeed",
        "pixelVolumes",
        "stripVolumes",
        "maxPixelHoles",
        "maxStripHoles",
        "trimTracks",
        "constrainToVolumes",
        "endOfWorldVolumes",
    ],
    defaults=[
        15.0,
        25.0,
        10,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ],
)

AmbiguityResolutionConfig = namedtuple(
    "AmbiguityResolutionConfig",
    ["maximumSharedHits", "nMeasurementsMin", "maximumIterations"],
    defaults=[None] * 3,
)

ScoreBasedAmbiguityResolutionConfig = namedtuple(
    "ScoreBasedAmbiguityResolutionConfig",
    [
        "minScore",
        "minScoreSharedTracks",
        "maxShared",
        "minUnshared",
        "maxSharedTracksPerMeasurement",
        "useAmbiguityScoring",
    ],
    defaults=[None] * 6,
)

AmbiguityResolutionMLConfig = namedtuple(
    "AmbiguityResolutionMLConfig",
    ["maximumSharedHits", "nMeasurementsMin", "maximumIterations"],
    defaults=[None] * 3,
)

SeedFilterMLDBScanConfig = namedtuple(
    "SeedFilterMLDBScanConfig",
    ["epsilonDBScan", "minPointsDBScan", "minSeedScore"],
    defaults=[None] * 3,
)

HashingTrainingConfigArg = namedtuple(
    "HashingTrainingConfig",
    ["annoySeed", "f"],
    defaults=[None] * 2,
)

HashingAlgorithmConfigArg = namedtuple(
    "HashingAlgorithmConfig",
    ["bucketSize", "zBins", "phiBins"],
    defaults=[None] * 3,
)


class VertexFinder(Enum):
    Truth = (1,)
    AMVF = (2,)
    Iterative = (3,)


@acts.examples.NamedTypeArgs(
    seedingAlgorithm=SeedingAlgorithm,
    trackSmearingSigmas=TrackSmearingSigmas,
    seedFinderConfigArg=SeedFinderConfigArg,
    seedFinderOptionsArg=SeedFinderOptionsArg,
    seedFilterConfigArg=SeedFilterConfigArg,
    spacePointGridConfigArg=SpacePointGridConfigArg,
    seedingAlgorithmConfigArg=SeedingAlgorithmConfigArg,
    hashingTrainingConfigArg=HashingTrainingConfigArg,
    hashingAlgorithmConfigArg=HashingAlgorithmConfigArg,
    truthEstimatedSeedingAlgorithmConfigArg=TruthEstimatedSeedingAlgorithmConfigArg,
    logLevel=acts.logging.Level,
)
def addSeeding(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    geoSelectionConfigFile: Optional[Union[Path, str]] = None,
    stripGeoSelectionConfigFile: Optional[Union[Path, str]] = None,
    layerMappingConfigFile: Optional[Union[Path, str]] = None,
    connectorInputConfigFile: Optional[Union[Path, str]] = None,
    lutInputConfigFile: Optional[Union[Path, str]] = None,
    seedingAlgorithm: SeedingAlgorithm = SeedingAlgorithm.GridTriplet,
    trackSmearingSigmas: TrackSmearingSigmas = TrackSmearingSigmas(),
    initialSigmas: Optional[list] = None,
    initialSigmaQoverPt: Optional[float] = None,
    initialSigmaPtRel: Optional[float] = None,
    initialVarInflation: Optional[list] = None,
    seedFinderConfigArg: SeedFinderConfigArg = SeedFinderConfigArg(),
    seedFinderOptionsArg: SeedFinderOptionsArg = SeedFinderOptionsArg(),
    seedFilterConfigArg: SeedFilterConfigArg = SeedFilterConfigArg(),
    spacePointGridConfigArg: SpacePointGridConfigArg = SpacePointGridConfigArg(),
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(),
    houghTransformConfig: acts.examples.HoughTransformSeeder.Config = acts.examples.HoughTransformSeeder.Config(),
    adaptiveHoughTransformConfig: Optional[
        acts.examples.AdaptiveHoughTransformSeeder.Config
    ] = None,
    hashingTrainingConfigArg: Optional[
        HashingTrainingConfigArg
    ] = HashingTrainingConfigArg(),
    hashingAlgorithmConfigArg: Optional[
        HashingAlgorithmConfigArg
    ] = HashingAlgorithmConfigArg(),
    truthEstimatedSeedingAlgorithmConfigArg: TruthEstimatedSeedingAlgorithmConfigArg = TruthEstimatedSeedingAlgorithmConfigArg(),
    particleHypothesis: Optional[
        acts.ParticleHypothesis
    ] = acts.ParticleHypothesis.pion,
    inputParticles: str = "particles",
    selectedParticles: str = "particles_selected",
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    rnd: Optional[acts.examples.RandomNumbers] = None,
) -> None:
    """This function steers the seeding
    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    geoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection. Not required for SeedingAlgorithm.TruthSmeared.
    stripGeoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection in strips. Needed for SpacePoint making.
    seedingAlgorithm : SeedingAlgorithm, Default
        seeding algorithm to use: one of Default (no truth information used), TruthSmeared, TruthEstimated
    trackSmearingSigmas : TrackSmearingSigmas(loc0, loc0PtA, loc0PtB, loc1, loc1PtA, loc1PtB, time, phi, theta, ptRel)
        TrackSmearing configuration.
        Defaults specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackParameterSmearing.hpp
    initialSigmas : list
        Sets the initial covariance matrix diagonal. This is ignored in case of TruthSmearing.
        Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp
    initialVarInflation : list
        List of 6 scale factors to inflate the initial covariance matrix
        Defaults (all 1) specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackParameterSmearing.hpp
    seedFinderConfigArg : SeedFinderConfigArg(maxSeedsPerSpM, cotThetaMax, sigmaScattering, radLengthPerSeed, minPt, impactMax, deltaPhiMax, interactionPointCut, deltaZMax, maxPtScattering, zBinEdges, zBinsCustomLooping, rRangeMiddleSP, useVariableMiddleSPRange, binSizeR, seedConfirmation, centralSeedConfirmationRange, forwardSeedConfirmationRange, deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z)
        SeedFinderConfig settings. deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z.
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFinderOptionsArg :  SeedFinderOptionsArg(bFieldInZ, beamPos)
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFilterConfigArg : SeedFilterConfigArg(compatSeedWeight, compatSeedLimit, numSeedIncrement, seedWeightIncrement, seedConfirmation, maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf, useDeltaRorTopRadius)
                                Defaults specified in Core/include/Acts/Seeding/SeedFilterConfig.hpp
    spacePointGridConfigArg : SpacePointGridConfigArg(rMax, zBinEdges, phiBinDeflectionCoverage, phi, maxPhiBins, impactMax)
                                SpacePointGridConfigArg settings. phi is specified as a tuple of (min,max).
        Defaults specified in Core/include/Acts/Seeding/SpacePointGrid.hpp
    seedingAlgorithmConfigArg : SeedingAlgorithmConfigArg(allowSeparateRMax, zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors, useExtraCuts)
                                Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/SeedingAlgorithm.hpp
    hashingTrainingConfigArg : HashingTrainingConfigArg(annoySeed, f)
                                Defaults specified in Plugins/Hashing/include/ActsPlugins/Hashing/HashingTrainingConfig.hpp
    hashingAlgorithmConfigArg : HashingAlgorithmConfigArg(bucketSize, zBins, phiBins)
                                Defaults specified in Plugins/Hashing/include/ActsPlugins/Hashing/HashingAlgorithmConfig.hpp
    truthEstimatedSeedingAlgorithmConfigArg : TruthEstimatedSeedingAlgorithmConfigArg(deltaR)
        Currently only deltaR=(min,max) range specified here.
    particleHypothesis : Optional[acts.ParticleHypothesis]
        The hypothesis used for track finding. Defaults to pion.
    inputParticles : str, "particles"
        input particles name in the WhiteBoard
    selectedParticles : str, "particles_selected"
        selected particles name in the WhiteBoard
    outputDirRoot : Path|str, path, None
        the output folder for ROOT output, None triggers no output
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    rnd : RandomNumbers, None
        random number generator. Only used by SeedingAlgorithm.TruthSmeared.
    """

    logLevel = acts.examples.defaultLogging(s, logLevel)()
    logger = acts.logging.getLogger("addSeeding")
    logger.setLevel(logLevel)

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if seedingAlgorithm == SeedingAlgorithm.TruthSmeared:
        logger.info("Using smeared truth particles for seeding")
        addTruthSmearedSeeding(
            s=s,
            rnd=rnd,
            selectedParticles=selectedParticles,
            trackSmearingSigmas=trackSmearingSigmas,
            initialSigmas=initialSigmas,
            initialSigmaQoverPt=initialSigmaQoverPt,
            initialSigmaPtRel=initialSigmaPtRel,
            initialVarInflation=initialVarInflation,
            particleHypothesis=particleHypothesis,
            logLevel=logLevel,
        )
    else:
        spacePoints = addSpacePointsMaking(
            s,
            trackingGeometry,
            geoSelectionConfigFile,
            stripGeoSelectionConfigFile,
            logLevel,
        )
        seeds = None
        perSeedParticleHypothesis = None
        # Run either: truth track finding or seeding
        if seedingAlgorithm == SeedingAlgorithm.TruthEstimated:
            logger.info("Using truth track finding from space points for seeding")
            seeds, perSeedParticleHypothesis = addTruthEstimatedSeeding(
                s,
                spacePoints,
                selectedParticles,
                truthEstimatedSeedingAlgorithmConfigArg,
                particleHypothesis=particleHypothesis,
                logLevel=logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.Default:
            logger.info("Using default seeding")
            seeds = addStandardSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.Orthogonal:
            logger.info("Using orthogonal seeding")
            seeds = addOrthogonalSeeding(
                s,
                spacePoints,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.HoughTransform:
            logger.info("Using Hough Transform seeding")
            houghTransformConfig.inputSpacePoints = [spacePoints]
            houghTransformConfig.inputMeasurements = "measurements"
            houghTransformConfig.outputProtoTracks = "prototracks"
            houghTransformConfig.outputSeeds = "seeds"
            houghTransformConfig.trackingGeometry = trackingGeometry
            seeds = addHoughTransformSeeding(s, houghTransformConfig, logLevel)
        elif seedingAlgorithm == SeedingAlgorithm.AdaptiveHoughTransform:
            logger.info("Using Adaptive Hough Transform seeding")
            adaptiveHoughTransformConfig.inputSpacePoints = [spacePoints]
            adaptiveHoughTransformConfig.outputProtoTracks = "prototracks"
            adaptiveHoughTransformConfig.outputSeeds = "seeds"
            adaptiveHoughTransformConfig.trackingGeometry = trackingGeometry
            adaptiveHoughTransformConfig.threshold = 4
            adaptiveHoughTransformConfig.noiseThreshold = 12
            adaptiveHoughTransformConfig.phiMinBinSize = 3.14 / (2.0 * 257.0)
            adaptiveHoughTransformConfig.qOverPtMinBinSize = 1.1 / (2.0 * 257.0)
            adaptiveHoughTransformConfig.qOverPtMin = 1.1
            adaptiveHoughTransformConfig.doSecondPhase = True
            adaptiveHoughTransformConfig.zMinBinSize = 1 * u.mm
            adaptiveHoughTransformConfig.cotThetaMinBinSize = 0.1
            adaptiveHoughTransformConfig.deduplicate = True
            seeds = addAdaptiveHoughTransformSeeding(
                s, adaptiveHoughTransformConfig, logLevel=logLevel
            )
        elif seedingAlgorithm == SeedingAlgorithm.Gbts:
            logger.info("Using Gbts seeding")
            # output of algs changed, only one output now
            seeds = addGbtsSeeding(
                s,
                spacePoints,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                trackingGeometry,
                logLevel,
                layerMappingConfigFile,
                geoSelectionConfigFile,
                connectorInputConfigFile,
                lutInputConfigFile,
            )
        elif seedingAlgorithm == SeedingAlgorithm.Hashing:
            logger.info("Using Hashing seeding")
            seeds, buckets = addHashingSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                hashingTrainingConfigArg,
                hashingAlgorithmConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.GridTriplet:
            logger.info("Using grid triplet seeding")
            seeds = addGridTripletSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.OrthogonalTriplet:
            logger.info("Using orthogonal triplet seeding")
            seeds = addOrthogonalTripletSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                logLevel,
            )
        else:
            logger.fatal("unknown seedingAlgorithm %s", seedingAlgorithm)

        parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
            level=logLevel,
            inputSeeds=seeds,
            inputParticleHypotheses=perSeedParticleHypothesis,
            outputTrackParameters="estimatedparameters",
            outputSeeds="estimatedseeds",
            trackingGeometry=trackingGeometry,
            magneticField=field,
            **acts.examples.defaultKWArgs(
                initialSigmas=initialSigmas,
                initialSigmaQoverPt=initialSigmaQoverPt,
                initialSigmaPtRel=initialSigmaPtRel,
                initialVarInflation=initialVarInflation,
                particleHypothesis=particleHypothesis,
            ),
        )
        s.addAlgorithm(parEstimateAlg)

        prototracks = "seed-prototracks"
        s.addAlgorithm(
            acts.examples.SeedsToPrototracks(
                level=logLevel,
                inputSeeds="estimatedseeds",
                outputProtoTracks=prototracks,
            )
        )

        tracks = "seed-tracks"
        s.addAlgorithm(
            acts.examples.PrototracksToTracks(
                level=logLevel,
                inputProtoTracks=prototracks,
                inputTrackParameters="estimatedparameters",
                inputMeasurements="measurements",
                outputTracks=tracks,
            )
        )

        s.addAlgorithm(
            acts.examples.TrackTruthMatcher(
                level=logLevel,
                inputTracks=tracks,
                inputParticles=selectedParticles,
                inputMeasurementParticlesMap="measurement_particles_map",
                outputTrackParticleMatching="seed_particle_matching",
                outputParticleTrackMatching="particle_seed_matching",
                matchingRatio=1.0,
                doubleMatching=False,
            )
        )

        if outputDirRoot is not None:
            addSeedPerformanceWriters(
                s,
                outputDirRoot,
                tracks,
                prototracks,
                selectedParticles,
                inputParticles,
                parEstimateAlg.config.outputTrackParameters,
                logLevel,
            )

        if outputDirCsv is not None:
            outputDirCsv = Path(outputDirCsv)

            if not outputDirCsv.exists():
                outputDirCsv.mkdir()

            csvSeedWriter = acts.examples.CsvSeedWriter(
                level=logLevel,
                inputTrackParameters=parEstimateAlg.config.outputTrackParameters,
                inputSimSeeds=seeds,
                inputSimHits="simhits",
                inputMeasurementParticlesMap="measurement_particles_map",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                outputDir=str(outputDirCsv),
                fileName=str(f"seed.csv"),
            )
            s.addWriter(csvSeedWriter)

            if seedingAlgorithm == SeedingAlgorithm.Hashing:
                s.addWriter(
                    acts.examples.CsvSpacePointsBucketWriter(
                        level=logLevel,
                        inputBuckets=buckets,
                        outputDir=str(outputDirCsv),
                    )
                )

    return s


def addTruthSmearedSeeding(
    s: acts.examples.Sequencer,
    rnd: Optional[acts.examples.RandomNumbers],
    selectedParticles: str,
    trackSmearingSigmas: TrackSmearingSigmas,
    initialSigmas: Optional[List[float]],
    initialSigmaQoverPt: Optional[float],
    initialSigmaPtRel: Optional[float],
    initialVarInflation: Optional[List[float]],
    particleHypothesis: Optional[acts.ParticleHypothesis],
    logLevel: acts.logging.Level = None,
):
    """adds algorithm that would mimic detector response uncertainties for truth seeding
    For parameters description see addSeeding
    """

    rnd = rnd or acts.examples.RandomNumbers(seed=42)

    trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
        level=logLevel,
        inputParticles=selectedParticles,
        outputTrackParameters="trueparameters",
    )
    s.addAlgorithm(trkParamExtractor)

    # Smearing track parameters
    trkSmear = acts.examples.TrackParameterSmearing(
        level=logLevel,
        inputTrackParameters=trkParamExtractor.config.outputTrackParameters,
        outputTrackParameters="estimatedparameters",
        randomNumbers=rnd,
        # gaussian sigmas to smear particle parameters
        **acts.examples.defaultKWArgs(
            sigmaLoc0=trackSmearingSigmas.loc0,
            sigmaLoc0PtA=trackSmearingSigmas.loc0PtA,
            sigmaLoc0PtB=trackSmearingSigmas.loc0PtB,
            sigmaLoc1=trackSmearingSigmas.loc1,
            sigmaLoc1PtA=trackSmearingSigmas.loc1PtA,
            sigmaLoc1PtB=trackSmearingSigmas.loc1PtB,
            sigmaTime=trackSmearingSigmas.time,
            sigmaPhi=trackSmearingSigmas.phi,
            sigmaTheta=trackSmearingSigmas.theta,
            sigmaPtRel=trackSmearingSigmas.ptRel,
            initialSigmas=initialSigmas,
            initialSigmaQoverPt=initialSigmaQoverPt,
            initialSigmaPtRel=initialSigmaPtRel,
            initialVarInflation=initialVarInflation,
            particleHypothesis=particleHypothesis,
        ),
    )
    s.addAlgorithm(trkSmear)

    truthTrkFndAlg = acts.examples.TruthTrackFinder(
        level=logLevel,
        inputParticles=selectedParticles,
        inputMeasurements="measurements",
        inputParticleMeasurementsMap="particle_measurements_map",
        inputSimHits="simhits",
        inputMeasurementSimHitsMap="measurement_simhits_map",
        outputProtoTracks="truth_particle_tracks",
    )
    s.addAlgorithm(truthTrkFndAlg)


def addTruthEstimatedSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    inputParticles: str,
    TruthEstimatedSeedingAlgorithmConfigArg: TruthEstimatedSeedingAlgorithmConfigArg,
    particleHypothesis: Optional[acts.ParticleHypothesis] = None,
    logLevel: acts.logging.Level = None,
):
    """adds truth seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()

    truthSeeding = acts.examples.TruthSeedingAlgorithm(
        level=logLevel,
        inputParticles=inputParticles,
        inputParticleMeasurementsMap="particle_measurements_map",
        inputSpacePoints=[spacePoints],
        inputSimHits="simhits",
        inputMeasurementSimHitsMap="measurement_simhits_map",
        outputParticles="truth_seeded_particles",
        outputProtoTracks="truth_particle_tracks",
        outputSeeds="seeds",
        outputParticleHypotheses="seed_particle_hypotheses",
        **acts.examples.defaultKWArgs(
            deltaRMin=TruthEstimatedSeedingAlgorithmConfigArg.deltaR[0],
            deltaRMax=TruthEstimatedSeedingAlgorithmConfigArg.deltaR[1],
            particleHypothesis=particleHypothesis,
        ),
    )
    sequence.addAlgorithm(truthSeeding)

    return truthSeeding.config.outputSeeds, truthSeeding.config.outputParticleHypotheses


def addSpacePointsMaking(
    sequence: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    geoSelectionConfigFile: Union[Path, str],
    stripGeoSelectionConfigFile: Union[Path, str],
    logLevel: acts.logging.Level = None,
):
    """adds space points making
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    spAlg = acts.examples.SpacePointMaker(
        level=logLevel,
        inputMeasurements="measurements",
        outputSpacePoints="spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.json.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
        stripGeometrySelection=(
            acts.examples.json.readJsonGeometryList(str(stripGeoSelectionConfigFile))
            if stripGeoSelectionConfigFile
            else []
        ),
    )
    sequence.addAlgorithm(spAlg)
    return spAlg.config.outputSpacePoints


def addStandardSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    spacePointGridConfigArg: SpacePointGridConfigArg,
    logLevel: acts.logging.Level = None,
    outputSeeds: str = "seeds",
):
    """adds standard seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()

    seedFinderConfig = acts.examples.SeedFinderConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMin=seedFinderConfigArg.deltaR[0],
            deltaRMax=seedFinderConfigArg.deltaR[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            deltaRMiddleMinSPRange=seedFinderConfigArg.deltaRMiddleSPRange[0],
            deltaRMiddleMaxSPRange=seedFinderConfigArg.deltaRMiddleSPRange[1],
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            zBinEdges=seedFinderConfigArg.zBinEdges,
            zBinsCustomLooping=seedFinderConfigArg.zBinsCustomLooping,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            binSizeR=seedFinderConfigArg.binSizeR,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfig.deltaRMin
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfig.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfig.forwardSeedConfirmationRange,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )

    gridConfig = acts.SpacePointGridConfig(
        **acts.examples.defaultKWArgs(
            minPt=seedFinderConfig.minPt,
            rMax=(
                seedFinderConfig.rMax
                if spacePointGridConfigArg.rMax is None
                else spacePointGridConfigArg.rMax
            ),
            zMax=seedFinderConfig.zMax,
            zMin=seedFinderConfig.zMin,
            deltaRMax=(
                seedFinderConfig.deltaRMax
                if spacePointGridConfigArg.deltaRMax is None
                else spacePointGridConfigArg.deltaRMax
            ),
            cotThetaMax=seedFinderConfig.cotThetaMax,
            phiMin=spacePointGridConfigArg.phi[0],
            phiMax=spacePointGridConfigArg.phi[1],
            maxPhiBins=spacePointGridConfigArg.maxPhiBins,
            impactMax=spacePointGridConfigArg.impactMax,
            zBinEdges=spacePointGridConfigArg.zBinEdges,
            phiBinDeflectionCoverage=spacePointGridConfigArg.phiBinDeflectionCoverage,
        )
    )

    gridOptions = acts.SpacePointGridOptions(
        **acts.examples.defaultKWArgs(
            bFieldInZ=seedFinderOptions.bFieldInZ,
        )
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds=outputSeeds,
        **acts.examples.defaultKWArgs(
            allowSeparateRMax=seedingAlgorithmConfigArg.allowSeparateRMax,
            zBinNeighborsTop=seedingAlgorithmConfigArg.zBinNeighborsTop,
            zBinNeighborsBottom=seedingAlgorithmConfigArg.zBinNeighborsBottom,
            numPhiNeighbors=seedingAlgorithmConfigArg.numPhiNeighbors,
            useExtraCuts=seedingAlgorithmConfigArg.useExtraCuts,
        ),
        gridConfig=gridConfig,
        gridOptions=gridOptions,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds


def addGridTripletSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    spacePointGridConfigArg: SpacePointGridConfigArg,
    logLevel: acts.logging.Level = None,
    outputSeeds: str = "seeds",
):
    """adds grid triplet seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()

    seedingAlg = acts.examples.GridTripletSeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=spacePoints,
        outputSeeds=outputSeeds,
        **acts.examples.defaultKWArgs(
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
            minPt=seedFinderConfigArg.minPt,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            impactMax=seedFinderConfigArg.impactMax,
            deltaRMin=seedFinderConfigArg.deltaR[0],
            deltaRMax=seedFinderConfigArg.deltaR[1],
            deltaRMinTop=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTop=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottom=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottom=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            phiMin=spacePointGridConfigArg.phi[0],
            phiMax=spacePointGridConfigArg.phi[1],
            phiBinDeflectionCoverage=spacePointGridConfigArg.phiBinDeflectionCoverage,
            maxPhiBins=spacePointGridConfigArg.maxPhiBins,
            zBinEdges=spacePointGridConfigArg.zBinEdges,
            zBinsCustomLooping=seedFinderConfigArg.zBinsCustomLooping,
            rMinMiddle=None,
            rMaxMiddle=None,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            deltaRMiddleMinSPRange=seedFinderConfigArg.deltaRMiddleSPRange[0],
            deltaRMiddleMaxSPRange=seedFinderConfigArg.deltaRMiddleSPRange[1],
            deltaZMin=None,
            deltaZMax=None,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            helixCutTolerance=None,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            toleranceParam=None,
            deltaInvHelixDiameter=None,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRinsteadOfTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
            useExtraCuts=seedingAlgorithmConfigArg.useExtraCuts,
        ),
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds


def addOrthogonalTripletSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    spacePointGridConfigArg: SpacePointGridConfigArg,
    logLevel: acts.logging.Level = None,
    outputSeeds: str = "seeds",
):
    """adds orthogonal triplet seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()

    seedingAlg = acts.examples.OrthogonalTripletSeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=spacePoints,
        outputSeeds=outputSeeds,
        **acts.examples.defaultKWArgs(
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
            minPt=seedFinderConfigArg.minPt,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            impactMax=seedFinderConfigArg.impactMax,
            deltaRMin=seedFinderConfigArg.deltaR[0],
            deltaRMax=seedFinderConfigArg.deltaR[1],
            deltaRMinTop=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTop=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottom=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottom=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            phiMin=spacePointGridConfigArg.phi[0],
            phiMax=spacePointGridConfigArg.phi[1],
            phiBinDeflectionCoverage=spacePointGridConfigArg.phiBinDeflectionCoverage,
            maxPhiBins=spacePointGridConfigArg.maxPhiBins,
            zBinEdges=spacePointGridConfigArg.zBinEdges,
            zBinsCustomLooping=seedFinderConfigArg.zBinsCustomLooping,
            rMinMiddle=None,
            rMaxMiddle=None,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            deltaRMiddleMinSPRange=seedFinderConfigArg.deltaRMiddleSPRange[0],
            deltaRMiddleMaxSPRange=seedFinderConfigArg.deltaRMiddleSPRange[1],
            deltaZMin=None,
            deltaZMax=None,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            helixCutTolerance=None,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            toleranceParam=None,
            deltaInvHelixDiameter=None,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRinsteadOfTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
            useExtraCuts=seedingAlgorithmConfigArg.useExtraCuts,
        ),
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds


def addOrthogonalSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    logLevel: acts.logging.Level = None,
):
    """adds orthogonal seeding algorithm
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    seedFinderConfig = acts.examples.SeedFinderOrthogonalConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            deltaPhiMax=seedFinderConfigArg.deltaPhiMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfigArg.deltaR[0]
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )
    seedingAlg = acts.examples.SeedingOrthogonalAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds


def addHashingSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    spacePointGridConfigArg: SpacePointGridConfigArg,
    hashingTrainingConfigArg: HashingTrainingConfigArg,
    hashingAlgorithmConfigArg: HashingAlgorithmConfigArg,
    logLevel: acts.logging.Level = None,
):
    """adds Hashing seeding
    For parameters description see addSeeding docstring
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    from acts.examples.hashing import SeedingAlgorithmHashing

    # Same configuration than the standard seeding
    seedFinderConfig = acts.examples.SeedFinderConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMin=seedFinderConfigArg.deltaR[0],
            deltaRMax=seedFinderConfigArg.deltaR[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            deltaRMiddleMinSPRange=seedFinderConfigArg.deltaRMiddleSPRange[0],
            deltaRMiddleMaxSPRange=seedFinderConfigArg.deltaRMiddleSPRange[1],
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            zBinEdges=seedFinderConfigArg.zBinEdges,
            zBinsCustomLooping=seedFinderConfigArg.zBinsCustomLooping,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            binSizeR=seedFinderConfigArg.binSizeR,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfig.deltaRMin
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfig.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfig.forwardSeedConfirmationRange,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )

    gridConfig = acts.SpacePointGridConfig(
        **acts.examples.defaultKWArgs(
            minPt=seedFinderConfig.minPt,
            rMax=(
                seedFinderConfig.rMax
                if spacePointGridConfigArg.rMax is None
                else spacePointGridConfigArg.rMax
            ),
            zMax=seedFinderConfig.zMax,
            zMin=seedFinderConfig.zMin,
            deltaRMax=(
                seedFinderConfig.deltaRMax
                if spacePointGridConfigArg.deltaRMax is None
                else spacePointGridConfigArg.deltaRMax
            ),
            cotThetaMax=seedFinderConfig.cotThetaMax,
            phiMin=spacePointGridConfigArg.phi[0],
            phiMax=spacePointGridConfigArg.phi[1],
            maxPhiBins=spacePointGridConfigArg.maxPhiBins,
            impactMax=spacePointGridConfigArg.impactMax,
            zBinEdges=spacePointGridConfigArg.zBinEdges,
            phiBinDeflectionCoverage=spacePointGridConfigArg.phiBinDeflectionCoverage,
        )
    )

    gridOptions = acts.SpacePointGridOptions(
        **acts.examples.defaultKWArgs(
            bFieldInZ=seedFinderOptions.bFieldInZ,
        )
    )

    # Hashing configuration
    hashingTrainingConfig = acts.hashing.HashingTrainingConfig(
        **acts.examples.defaultKWArgs(
            annoySeed=hashingTrainingConfigArg.annoySeed,
            f=hashingTrainingConfigArg.f,
        ),
    )

    hashingConfig = acts.hashing.HashingAlgorithmConfig(
        **acts.examples.defaultKWArgs(
            bucketSize=hashingAlgorithmConfigArg.bucketSize,
            zBins=hashingAlgorithmConfigArg.zBins,
            phiBins=hashingAlgorithmConfigArg.phiBins,
        ),
    )

    # Seeding algorithm
    seedingAlg = SeedingAlgorithmHashing(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        outputBuckets="buckets",
        **acts.examples.defaultKWArgs(
            allowSeparateRMax=seedingAlgorithmConfigArg.allowSeparateRMax,
            zBinNeighborsTop=seedingAlgorithmConfigArg.zBinNeighborsTop,
            zBinNeighborsBottom=seedingAlgorithmConfigArg.zBinNeighborsBottom,
            numPhiNeighbors=seedingAlgorithmConfigArg.numPhiNeighbors,
        ),
        gridConfig=gridConfig,
        gridOptions=gridOptions,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
        hashingConfig=hashingConfig,
        hashingTrainingConfig=hashingTrainingConfig,
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds, seedingAlg.config.outputBuckets


def addHoughTransformSeeding(
    sequence: acts.examples.Sequencer,
    config: acts.examples.HoughTransformSeeder.Config,
    logLevel: acts.logging.Level = None,
):
    """
    Configures HoughTransform (HT) for seeding, instead of extra proxy config objects it takes
    directly the HT example algorithm config.
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    ht = acts.examples.HoughTransformSeeder(config=config, level=logLevel)
    sequence.addAlgorithm(ht)
    # potentially HT can be extended to also produce seeds, but it is not implemented yet
    # configuration option (outputSeeds) exists
    return ht.config.outputSeeds


def addAdaptiveHoughTransformSeeding(
    sequence: acts.examples.Sequencer,
    config: acts.examples.AdaptiveHoughTransformSeeder.Config,
    logLevel: acts.logging.Level = None,
):
    """
    Configures AdaptiveHoughTransform (AHT) for seeding, instead of extra proxy config objects it takes
    directly the AHT example algorithm config.
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    ht = acts.examples.AdaptiveHoughTransformSeeder(config=config, level=logLevel)
    sequence.addAlgorithm(ht)
    # potentially HT can be extended to also produce proto-tracks, but it is not implemented yet
    # configuration option (outputSeeds) exists and is used
    return ht.config.outputSeeds


def addGbtsSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    trackingGeometry: acts.TrackingGeometry,
    logLevel: acts.logging.Level = None,
    layerMappingConfigFile: Union[Path, str] = None,
    geoSelectionConfigFile: Union[Path, str] = None,
    connectorInputConfigFile: Union[Path, str] = None,
    lutInputConfigFile: Optional[Union[Path, str]] = None,
):
    """Gbts seeding"""

    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    layerMappingFile = str(layerMappingConfigFile)  # turn path into string
    connectorInputFileStr = str(connectorInputConfigFile)
    lutInputConfigFileStr = str(lutInputConfigFile)
    seedFinderConfig = acts.examples.SeedFinderGbtsConfig(
        **acts.examples.defaultKWArgs(
            minPt=seedFinderConfigArg.minPt,
            connectorInputFile=connectorInputFileStr,
            lutInputFile=lutInputConfigFileStr,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )

    seedingAlg = acts.examples.GbtsSeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=spacePoints,
        outputSeeds="seeds",
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
        layerMappingFile=layerMappingFile,
        trackingGeometry=trackingGeometry,
        fill_module_csv=False,
        inputClusters="clusters",
    )

    sequence.addAlgorithm(seedingAlg)
    return seedingAlg.config.outputSeeds


def addSeedPerformanceWriters(
    sequence: acts.examples.Sequencer,
    outputDirRoot: Union[Path, str],
    tracks: str,
    prototracks: str,
    selectedParticles: str,
    inputParticles: str,
    outputTrackParameters: str,
    logLevel: acts.logging.Level = None,
):
    """Writes seeding related performance output"""
    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    assert ACTS_ROOT_AVAILABLE, "ROOT output requested but ROOT is not available"
    outputDirRoot = Path(outputDirRoot)
    if not outputDirRoot.exists():
        outputDirRoot.mkdir()

    sequence.addWriter(
        RootTrackFinderPerformanceWriter(
            level=customLogLevel(),
            inputTracks=tracks,
            inputParticles=selectedParticles,
            inputTrackParticleMatching="seed_particle_matching",
            inputParticleTrackMatching="particle_seed_matching",
            inputParticleMeasurementsMap="particle_measurements_map",
            filePath=str(outputDirRoot / f"performance_seeding.root"),
        )
    )

    sequence.addWriter(
        RootTrackParameterWriter(
            level=customLogLevel(),
            inputTrackParameters=outputTrackParameters,
            inputProtoTracks=prototracks,
            inputParticles=inputParticles,
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDirRoot / "estimatedparams.root"),
            treeName="estimatedparams",
        )
    )


acts.examples.NamedTypeArgs(
    config=SeedFilterMLDBScanConfig,
)


def addSeedFilterML(
    s,
    config: SeedFilterMLDBScanConfig = SeedFilterMLDBScanConfig(),
    onnxModelFile: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)()
    from acts.examples.onnx import SeedFilterMLAlgorithm

    inputParticles = "particles"
    selectedParticles = "particles_selected"
    seeds = "seeds"
    estParams = "estimatedparameters"
    prototracks = "seed-prototracks-ML"
    tracks = "seed-tracks-ML"

    filterML = SeedFilterMLAlgorithm(
        level=customLogLevel,
        inputTrackParameters="estimatedparameters",
        inputSimSeeds="seeds",
        inputSeedFilterNN=onnxModelFile,
        outputTrackParameters="filtered-parameters",
        outputSimSeeds="filtered-seeds",
        **acts.examples.defaultKWArgs(
            epsilonDBScan=config.epsilonDBScan,
            minPointsDBScan=config.minPointsDBScan,
            minSeedScore=config.minSeedScore,
        ),
    )
    s.addAlgorithm(filterML)
    s.addWhiteboardAlias(seeds, "filtered-seeds")
    s.addWhiteboardAlias("estimatedparameters", "filtered-parameters")

    s.addAlgorithm(
        acts.examples.SeedsToPrototracks(
            level=customLogLevel,
            inputSeeds=seeds,
            outputProtoTracks=prototracks,
        )
    )

    s.addAlgorithm(
        acts.examples.PrototracksToTracks(
            level=customLogLevel,
            inputProtoTracks=prototracks,
            inputTrackParameters="estimatedparameters",
            outputTracks=tracks,
        )
    )

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=customLogLevel,
            inputTracks=tracks,
            inputParticles=selectedParticles,
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="seed_particle_matching",
            outputParticleTrackMatching="particle_seed_matching",
            matchingRatio=1.0,
            doubleMatching=False,
        )
    )

    if outputDirRoot is not None:
        addSeedPerformanceWriters(
            s,
            outputDirRoot,
            tracks,
            selectedParticles,
            inputParticles,
            estParams,
            customLogLevel,
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)

        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        csvSeedWriter = acts.examples.CsvSeedWriter(
            level=customLogLevel,
            inputTrackParameters=estParams,
            inputSimSeeds=seeds,
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            outputDir=str(outputDirCsv),
            fileName=str(f"seed.csv"),
        )
        s.addWriter(csvSeedWriter)

    return s


def addKalmanTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    reverseFilteringMomThreshold: float = 0 * u.GeV,
    reverseFilteringCovarianceScaling: float = 100.0,
    inputProtoTracks: str = "truth_particle_tracks",
    multipleScattering: bool = True,
    energyLoss: bool = True,
    clusters: str = None,
    calibrator: acts.examples.MeasurementCalibrator = acts.examples.makePassThroughCalibrator(),
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    kalmanOptions = {
        "multipleScattering": multipleScattering,
        "energyLoss": energyLoss,
        "reverseFilteringMomThreshold": reverseFilteringMomThreshold,
        "reverseFilteringCovarianceScaling": reverseFilteringCovarianceScaling,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "level": customLogLevel(),
        "chi2Cut": float("inf"),
    }

    fitAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        inputClusters=clusters if clusters is not None else "",
        outputTracks="kf_tracks",
        pickTrack=-1,
        fit=acts.examples.makeKalmanFitterFunction(
            trackingGeometry, field, **kalmanOptions
        ),
        calibrator=calibrator,
    )
    s.addAlgorithm(fitAlg)
    s.addWhiteboardAlias("tracks", fitAlg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=fitAlg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="kf_track_particle_matching",
        outputParticleTrackMatching="kf_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    return s


def addTruthTrackingGsf(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    inputProtoTracks: str = "truth_particle_tracks",
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

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
        "level": customLogLevel(),
    }

    gsfAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        outputTracks="gsf_tracks",
        pickTrack=-1,
        fit=acts.examples.makeGsfFitterFunction(trackingGeometry, field, **gsfOptions),
        calibrator=acts.examples.makePassThroughCalibrator(),
    )
    s.addAlgorithm(gsfAlg)
    s.addWhiteboardAlias("tracks", gsfAlg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=gsfAlg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="gsf_track_particle_matching",
        outputParticleTrackMatching="gsf_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    return s


@acts.examples.NamedTypeArgs(
    trackSelectorConfig=TrackSelectorConfig,
    ckfConfig=CkfConfig,
)
def addCKFTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    trackSelectorConfig: Optional[
        Union[TrackSelectorConfig, List[TrackSelectorConfig]]
    ] = None,
    ckfConfig: CkfConfig = CkfConfig(),
    twoWay: bool = True,
    reverseSearch: bool = False,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat=False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for ROOT output, None triggers no output
    trackSelectorConfig : TrackSelectorConfig(loc0, loc1, time, eta, absEta, pt, phi, minMeasurements)
        TrackSelector configuration. Each range is specified as a tuple of (min,max).
        Specify as a list(TrackSelectorConfig) for eta-dependent cuts, with binning specified by absEta[1].
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackSelector.hpp
    writeTrackSummary : bool, True
        write tracksummary_ckf.root ntuple?
    writeTrackStates : bool, False
        write trackstates_ckf.root ntuple? This can be quite large.
    writePerformance : bool, True
        write performance_fitting_ckf.root and performance_finding_ckf.root ntuples?
    writeCovMat : bool, False
        write covaraiance matrices to tracksummary_ckf.root ntuple?
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    tslist = (
        []
        if trackSelectorConfig is None
        else (
            [trackSelectorConfig]
            if type(trackSelectorConfig) is TrackSelectorConfig
            else trackSelectorConfig
        )
    )

    if len(tslist) > 1:
        cutSets = []
        for c in tslist:
            defKW = trackSelectorDefaultKWArgs(c)
            defKW.pop("absEtaMax", None)
            cutSets += [acts.TrackSelector.Config(**(defKW))]
    else:
        cutSets = [
            acts.TrackSelector.Config(**(trackSelectorDefaultKWArgs(c))) for c in tslist
        ]

    if len(tslist) == 0:
        trkSelCfg = None
    elif len(tslist) == 1:
        trkSelCfg = cutSets[0]
    else:
        trkSelCfg = acts.TrackSelector.EtaBinnedConfig(
            cutSets=cutSets,
            absEtaEdges=[cutSets[0].absEtaMin] + [c.absEta[1] for c in tslist],
        )

    # Set up the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=customLogLevel(),
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [
                (
                    acts.GeometryIdentifier(),
                    (
                        [],
                        [ckfConfig.chi2CutOffMeasurement],
                        [ckfConfig.chi2CutOffOutlier],
                        [ckfConfig.numMeasurementsCutOff],
                    ),
                )
            ]
        ),
        inputMeasurements="measurements",
        inputInitialTrackParameters="estimatedparameters",
        inputSeeds=(
            "estimatedseeds"
            if ckfConfig.seedDeduplication or ckfConfig.stayOnSeed
            else ""
        ),
        outputTracks="ckf_tracks",
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field, customLogLevel()
        ),
        **acts.examples.defaultKWArgs(
            trackingGeometry=trackingGeometry,
            magneticField=field,
            trackSelectorCfg=trkSelCfg,
            maxSteps=ckfConfig.maxSteps,
            twoWay=twoWay,
            reverseSearch=reverseSearch,
            seedDeduplication=ckfConfig.seedDeduplication,
            stayOnSeed=ckfConfig.stayOnSeed,
            pixelVolumeIds=ckfConfig.pixelVolumes,
            stripVolumeIds=ckfConfig.stripVolumes,
            maxPixelHoles=ckfConfig.maxPixelHoles,
            maxStripHoles=ckfConfig.maxStripHoles,
            trimTracks=ckfConfig.trimTracks,
            constrainToVolumeIds=ckfConfig.constrainToVolumes,
            endOfWorldVolumeIds=ckfConfig.endOfWorldVolumes,
        ),
    )
    s.addAlgorithm(trackFinder)
    s.addWhiteboardAlias("tracks", trackFinder.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=trackFinder.config.outputTracks,
        inputParticles="particles_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ckf_track_particle_matching",
        outputParticleTrackMatching="ckf_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name="ckf",
        tracks=trackFinder.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=writePerformance,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s


def addGx2fTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    inputProtoTracks: str = "truth_particle_tracks",
    multipleScattering: bool = False,
    energyLoss: bool = False,
    nUpdateMax: int = 5,
    relChi2changeCutOff: float = 1e-7,
    clusters: str = None,
    calibrator: acts.examples.MeasurementCalibrator = acts.examples.makePassThroughCalibrator(),
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    gx2fOptions = {
        "multipleScattering": multipleScattering,
        "energyLoss": energyLoss,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "nUpdateMax": nUpdateMax,
        "relChi2changeCutOff": relChi2changeCutOff,
        "level": customLogLevel(),
    }

    fitAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        inputClusters=clusters if clusters is not None else "",
        outputTracks="gx2f_tracks",
        pickTrack=-1,
        fit=acts.examples.makeGlobalChiSquareFitterFunction(
            trackingGeometry, field, **gx2fOptions
        ),
        calibrator=calibrator,
    )
    s.addAlgorithm(fitAlg)
    s.addWhiteboardAlias("tracks", fitAlg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=fitAlg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="gx2f_track_particle_matching",
        outputParticleTrackMatching="gx2f_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    return s


def addTrackWriters(
    s: acts.examples.Sequencer,
    name: str,
    tracks: str = "tracks",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeSummary: bool = True,
    writeStates: bool = False,
    writeFitterPerformance: bool = False,
    writeFinderPerformance: bool = False,
    logLevel: Optional[acts.logging.Level] = None,
    writeCovMat=False,
):
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if outputDirRoot is not None:
        assert ACTS_ROOT_AVAILABLE, "ROOT output requested but ROOT is not available"
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        if writeSummary:
            trackSummaryWriter = RootTrackSummaryWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                filePath=str(outputDirRoot / f"tracksummary_{name}.root"),
                treeName="tracksummary",
                writeCovMat=writeCovMat,
            )
            s.addWriter(trackSummaryWriter)

        if writeStates:
            trackStatesWriter = RootTrackStatesWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                inputSimHits="simhits",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                filePath=str(outputDirRoot / f"trackstates_{name}.root"),
                treeName="trackstates",
            )
            s.addWriter(trackStatesWriter)

        if writeFitterPerformance:
            trackFitterPerformanceWriter = RootTrackFitterPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                filePath=str(outputDirRoot / f"performance_fitting_{name}.root"),
            )
            s.addWriter(trackFitterPerformanceWriter)

        if writeFinderPerformance:
            trackFinderPerfWriter = RootTrackFinderPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                inputParticleTrackMatching="particle_track_matching",
                inputParticleMeasurementsMap="particle_measurements_map",
                filePath=str(outputDirRoot / f"performance_finding_{name}.root"),
            )
            s.addWriter(trackFinderPerfWriter)

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        if writeSummary:
            csvWriter = acts.examples.CsvTrackWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputMeasurementParticlesMap="measurement_particles_map",
                outputDir=str(outputDirCsv),
                fileName=str(f"tracks_{name}.csv"),
            )
            s.addWriter(csvWriter)

            trackParameterWriter = acts.examples.CsvTrackParameterWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                outputDir=str(outputDirCsv),
                outputStem=str(f"track_parameters_{name}"),
            )
            s.addWriter(trackParameterWriter)


@acts.examples.NamedTypeArgs(
    trackSelectorConfig=TrackSelectorConfig,
)
def addTrackSelection(
    s: acts.examples.Sequencer,
    trackSelectorConfig: TrackSelectorConfig,
    inputTracks: str,
    outputTracks: str,
    logLevel: Optional[acts.logging.Level] = None,
) -> acts.examples.TrackSelectorAlgorithm:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # single cut config for implicit single bin eta configuration
    selectorConfig = acts.TrackSelector.Config(
        **trackSelectorDefaultKWArgs(trackSelectorConfig)
    )

    trackSelector = acts.examples.TrackSelectorAlgorithm(
        level=customLogLevel(),
        inputTracks=inputTracks,
        outputTracks=outputTracks,
        selectorConfig=selectorConfig,
    )

    s.addAlgorithm(trackSelector)

    return trackSelector


GnnBackend = Enum("GnnBackend", "Torch Onnx")


def addGnn(
    s: acts.examples.Sequencer,
    graphConstructor,
    edgeClassifiers: list,
    trackBuilder,
    nodeFeatures: list,
    featureScales: list,
    inputSpacePoints: str = "spacepoints",
    inputClusters: str = "",
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> acts.examples.Sequencer:
    """
    Add GNN track finding with custom stage implementations.

    This is a flexible low-level API that accepts pre-configured GNN stage components.
    For examples of how to configure stages, see gnn_metric_learning.py and gnn_module_map.py.

    Args:
        s: Sequencer to add algorithms to
        graphConstructor: Graph construction stage (TorchMetricLearning, ModuleMapCuda, etc.)
        edgeClassifiers: List of edge classification stages (run sequentially)
        trackBuilder: Track building stage (BoostTrackBuilding, CudaTrackBuilding, etc.)
        nodeFeatures: List of node features to extract from spacepoints/clusters
        featureScales: Scaling factors for each feature
        trackingGeometry: Optional tracking geometry for creating spacepoints
        geometrySelection: Optional geometry selection file for spacepoint creation
        inputSpacePoints: Name of input spacepoint collection (default: "spacepoints")
        inputClusters: Name of input cluster collection (default: "")
        outputDirRoot: Optional output directory for performance ROOT files
        logLevel: Logging level

    Note:
        The trackingGeometry parameter serves two distinct purposes depending on the workflow:
        1. Spacepoint creation: When provided along with geometrySelection, creates spacepoints
           from measurements using SpacePointMaker (typical for simulation workflows)
        2. Module map usage: Some graph constructors (e.g., ModuleMapCuda) require
           trackingGeometry to map module IDs even when using pre-existing spacepoints
    """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Validate that nodeFeatures and featureScales have matching lengths
    if len(nodeFeatures) != len(featureScales):
        raise ValueError(
            f"nodeFeatures and featureScales must have the same length "
            f"(got {len(nodeFeatures)} and {len(featureScales)})"
        )

    # GNN track finding algorithm
    findingAlg = acts.examples.gnn.TrackFindingAlgorithmGnn(
        level=customLogLevel(),
        inputSpacePoints=inputSpacePoints,
        inputClusters=inputClusters,
        outputProtoTracks="gnn_prototracks",
        graphConstructor=graphConstructor,
        edgeClassifiers=edgeClassifiers,
        trackBuilder=trackBuilder,
        nodeFeatures=nodeFeatures,
        featureScales=featureScales,
    )
    s.addAlgorithm(findingAlg)
    s.addWhiteboardAlias("prototracks", findingAlg.config.outputProtoTracks)

    # Convert prototracks to tracks
    s.addAlgorithm(
        acts.examples.PrototracksToTracks(
            level=customLogLevel(),
            inputProtoTracks="prototracks",
            inputMeasurements="measurements",
            outputTracks="tracks",
        )
    )

    # Truth matching
    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks="tracks",
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="gnn_track_particle_matching",
        outputParticleTrackMatching="gnn_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    # Optional performance writer
    if outputDirRoot is not None:
        assert ACTS_ROOT_AVAILABLE, "ROOT output requested but ROOT is not available"
        s.addWriter(
            RootTrackFinderNTupleWriter(
                level=customLogLevel(),
                inputTracks="tracks",
                inputParticles="particles",
                inputParticleMeasurementsMap="particle_measurements_map",
                inputTrackParticleMatching=matchAlg.config.outputTrackParticleMatching,
                filePath=str(Path(outputDirRoot) / "performance_track_finding.root"),
            )
        )

    return s


@acts.examples.NamedTypeArgs(
    config=AmbiguityResolutionConfig,
)
def addAmbiguityResolution(
    s,
    config: AmbiguityResolutionConfig = AmbiguityResolutionConfig(),
    tracks: str = "tracks",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat=False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    from acts.examples import GreedyAmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    alg = GreedyAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=tracks,
        outputTracks="ambi_tracks",
        **acts.examples.defaultKWArgs(
            maximumSharedHits=config.maximumSharedHits,
            nMeasurementsMin=config.nMeasurementsMin,
            maximumIterations=config.maximumIterations,
        ),
    )
    s.addAlgorithm(alg)
    s.addWhiteboardAlias("tracks", alg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=alg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ambi_track_particle_matching",
        outputParticleTrackMatching="ambi_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name="ambi",
        tracks=alg.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=writePerformance,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    config=ScoreBasedAmbiguityResolutionConfig,
)
def addScoreBasedAmbiguityResolution(
    s,
    config: ScoreBasedAmbiguityResolutionConfig = ScoreBasedAmbiguityResolutionConfig(),
    tracks: str = "tracks",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    ambiVolumeFile: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat=False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    from acts.examples import ScoreBasedAmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, acts.logging.INFO)

    algScoreBased = ScoreBasedAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=tracks,
        configFile=ambiVolumeFile,
        outputTracks="ambiTracksScoreBased",
        **acts.examples.defaultKWArgs(
            minScore=config.minScore,
            minScoreSharedTracks=config.minScoreSharedTracks,
            maxShared=config.maxShared,
            minUnshared=config.minUnshared,
            maxSharedTracksPerMeasurement=config.maxSharedTracksPerMeasurement,
            useAmbiguityScoring=config.useAmbiguityScoring,
        ),
    )
    s.addAlgorithm(algScoreBased)
    s.addWhiteboardAlias("tracks", algScoreBased.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=algScoreBased.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ambi_scorebased_track_particle_matching",
        outputParticleTrackMatching="ambi_scorebased_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name="ambi_scorebased",
        tracks=algScoreBased.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=writePerformance,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    config=AmbiguityResolutionMLConfig,
)
def addAmbiguityResolutionML(
    s,
    config: AmbiguityResolutionMLConfig = AmbiguityResolutionMLConfig(),
    tracks: str = "tracks",
    onnxModelFile: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat=False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    from acts.examples.onnx import AmbiguityResolutionMLAlgorithm
    from acts.examples import GreedyAmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, logLevel)
    algML = AmbiguityResolutionMLAlgorithm(
        level=customLogLevel(),
        inputTracks=tracks,
        inputDuplicateNN=onnxModelFile,
        outputTracks="ambiTracksML",
        **acts.examples.defaultKWArgs(
            nMeasurementsMin=config.nMeasurementsMin,
        ),
    )

    algGreedy = GreedyAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=algML.config.outputTracks,
        outputTracks="ambiTracksMLGreedy",
        **acts.examples.defaultKWArgs(
            maximumSharedHits=config.maximumSharedHits,
            nMeasurementsMin=config.nMeasurementsMin,
            maximumIterations=config.maximumIterations,
        ),
    )

    s.addAlgorithm(algML)
    s.addAlgorithm(algGreedy)

    s.addWhiteboardAlias("tracks", algGreedy.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=algGreedy.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ambiML_track_particle_matching",
        outputParticleTrackMatching="ambiML_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name="ambiML",
        tracks=algGreedy.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=writePerformance,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    trackSelectorConfig=TrackSelectorConfig,
)
def addVertexFitting(
    s,
    field,
    tracks: Optional[str] = "tracks",
    trackParameters: Optional[str] = None,
    outputProtoVertices: str = "protovertices",
    outputVertices: str = "fittedVertices",
    vertexFinder: VertexFinder = VertexFinder.Truth,
    maxIterations: Optional[int] = None,
    initialVariances: Optional[List[float]] = None,
    useTime: Optional[bool] = False,
    seeder: Optional[
        acts.examples.VertexSeedFinder
    ] = acts.examples.VertexSeedFinder.GaussianSeeder,
    spatialBinExtent: Optional[float] = None,
    temporalBinExtent: Optional[float] = None,
    simultaneousSeeds: Optional[int] = None,
    trackSelectorConfig: Optional[TrackSelectorConfig] = None,
    writeTrackInfo: bool = False,
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the vertex fitting

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from
        addVertexFitting)
    field : magnetic field
    outputDirRoot : Path|str, path, None
        the output folder for ROOT output, None triggers no output
    outputDirCsv : Path|str, path, None
        the output folder for the CSV output, None triggers no output
    vertexFinder : VertexFinder, Truth
        vertexFinder algorithm: one of Truth, AMVF, Iterative
    seeder : enum member
        determines vertex seeder for AMVF, can be acts.seeder.GaussianSeeder or
        acts.seeder.AdaptiveGridSeeder
    useTime : bool, False
        determines whether time information is used in vertex seeder, finder,
        and fitter
        only implemented for the AMVF and the AdaptiveGridSeeder
    spatialBinExtent : float, None
        spatial bin extent for the AdaptiveGridSeeder
    temporalBinExtent : float, None
        temporal bin extent for the AdaptiveGridSeeder
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    """
    from acts.examples import (
        TruthVertexFinder,
        VertexFitterAlgorithm,
        IterativeVertexFinderAlgorithm,
        AdaptiveMultiVertexFinderAlgorithm,
        CsvVertexWriter,
    )

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if tracks is not None and trackSelectorConfig is not None:
        trackSelector = addTrackSelection(
            s,
            trackSelectorConfig,
            inputTracks=tracks,
            outputTracks="selectedTracksVertexing",
            logLevel=customLogLevel(),
        )
        tracks = trackSelector.config.outputTracks

    if trackParameters is None:
        converter = acts.examples.TracksToParameters(
            level=customLogLevel(),
            inputTracks=tracks,
            outputTrackParameters="selectedTracksParametersVertexing",
        )
        s.addAlgorithm(converter)
        trackParameters = converter.config.outputTrackParameters

    tracks = tracks if tracks is not None else ""
    inputParticles = "particles"
    selectedParticles = "particles_selected"
    inputTruthVertices = "vertices_truth"

    if vertexFinder == VertexFinder.Truth:
        findVertices = TruthVertexFinder(
            level=customLogLevel(),
            inputTracks=tracks,
            inputTrackParticleMatching="track_particle_matching",
            outputProtoVertices=outputProtoVertices,
            excludeSecondaries=True,
        )
        s.addAlgorithm(findVertices)
        fitVertices = VertexFitterAlgorithm(
            level=customLogLevel(),
            inputTrackParameters=trackParameters,
            inputProtoVertices=findVertices.config.outputProtoVertices,
            outputVertices=outputVertices,
            bField=field,
        )
        s.addAlgorithm(fitVertices)
    elif vertexFinder == VertexFinder.Iterative:
        findVertices = IterativeVertexFinderAlgorithm(
            level=customLogLevel(),
            inputTrackParameters=trackParameters,
            outputProtoVertices=outputProtoVertices,
            outputVertices=outputVertices,
            bField=field,
        )
        s.addAlgorithm(findVertices)
    elif vertexFinder == VertexFinder.AMVF:
        findVertices = AdaptiveMultiVertexFinderAlgorithm(
            level=customLogLevel(),
            inputTrackParameters=trackParameters,
            inputTruthParticles=selectedParticles,
            inputTruthVertices=inputTruthVertices,
            outputProtoVertices=outputProtoVertices,
            outputVertices=outputVertices,
            bField=field,
            seedFinder=seeder,
            **acts.examples.defaultKWArgs(
                maxIterations=maxIterations,
                initialVariances=initialVariances,
                useTime=useTime,
                spatialBinExtent=spatialBinExtent,
                temporalBinExtent=temporalBinExtent,
                simultaneousSeeds=simultaneousSeeds,
            ),
        )
        s.addAlgorithm(findVertices)
    else:
        raise RuntimeError("Invalid finder argument")

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()
        s.addWriter(
            CsvVertexWriter(
                level=customLogLevel(),
                inputVertices=outputVertices,
                outputDir=str(outputDirCsv),
                outputStem="vertices",
            )
        )

    if outputDirRoot is not None:
        assert ACTS_ROOT_AVAILABLE, "ROOT output requested but ROOT is not available"
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        s.addWriter(
            RootVertexNTupleWriter(
                level=customLogLevel(),
                inputVertices=outputVertices,
                inputTracks=tracks,
                inputTruthVertices=inputTruthVertices,
                inputParticles=inputParticles,
                inputSelectedParticles=selectedParticles,
                inputTrackParticleMatching="track_particle_matching",
                bField=field,
                writeTrackInfo=writeTrackInfo,
                treeName="vertexing",
                filePath=str(outputDirRoot / "performance_vertexing.root"),
            )
        )

    return s


def addHoughVertexFinding(
    s,
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    inputSpacePoints: Optional[str] = "spacepoints",
    outputVertices: Optional[str] = "fittedHoughVertices",
) -> None:
    from acts.examples import (
        HoughVertexFinderAlgorithm,
        RootVertexNTupleWriter,
    )

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    findHoughVertex = HoughVertexFinderAlgorithm(
        level=customLogLevel(),
        inputSpacepoints=inputSpacePoints,
        outputVertices=outputVertices,
    )
    s.addAlgorithm(findHoughVertex)

    inputParticles = "particles"
    selectedParticles = "particles_selected"
    inputTruthVertices = "vertices_truth"

    if outputDirRoot is not None:
        assert ACTS_ROOT_AVAILABLE, "ROOT output requested but ROOT is not available"
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        s.addWriter(
            RootVertexNTupleWriter(
                level=customLogLevel(),
                inputParticles=inputParticles,
                inputSelectedParticles=selectedParticles,
                inputTracks="",
                inputTrackParticleMatching="",
                writeTrackInfo=False,
                inputTruthVertices=inputTruthVertices,
                inputVertices=outputVertices,
                treeName="houghvertexing",
                filePath=str(outputDirRoot / "performance_houghvertexing.root"),
            )
        )

    return s
