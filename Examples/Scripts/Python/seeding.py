#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union
from enum import Enum
from collections import namedtuple
import argparse

import acts
import acts.examples


# Graciously taken from https://stackoverflow.com/a/60750535/4280680
class EnumAction(argparse.Action):
    """
    Argparse action for handling Enums
    """

    def __init__(self, **kwargs):
        # Pop off the type value
        enum_type = kwargs.pop("enum", None)

        # Ensure an Enum subclass is provided
        if enum_type is None:
            raise ValueError("type must be assigned an Enum when using EnumAction")
        if not issubclass(enum_type, Enum):
            raise TypeError("type must be an Enum when using EnumAction")

        # Generate choices from the Enum
        kwargs.setdefault("choices", tuple(e.name for e in enum_type))

        super(EnumAction, self).__init__(**kwargs)

        self._enum = enum_type

    def __call__(self, parser, namespace, values, option_string=None):
        for e in self._enum:
            if e.name == values:
                setattr(namespace, self.dest, e)
                break
        else:
            raise ValueError("%s is not a validly enumerated algorithm." % values)


u = acts.UnitConstants
SeedingAlgorithm = Enum(
    "SeedingAlgorithm", "Default TruthSmeared TruthEstimated Orthogonal"
)

TruthSeedRanges = namedtuple(
    "TruthSeedRanges",
    ["rho", "z", "phi", "eta", "absEta", "pt", "nHits"],
    defaults=[(None, None)] * 7,
)

ParticleSmearingSigmas = namedtuple(
    "ParticleSmearingSigmas",
    ["d0", "d0PtA", "d0PtB", "z0", "z0PtA", "z0PtB", "t0", "phi", "theta", "pRel"],
    defaults=[None] * 10,
)

SeedfinderConfigArg = namedtuple(
    "SeedfinderConfig",
    [
        "maxSeedsPerSpM",
        "cotThetaMax",
        "sigmaScattering",
        "radLengthPerSeed",
        "minPt",
        "bFieldInZ",
        "impactMax",
        "deltaR",  # (min,max)
        "collisionRegion",  # (min,max)
        "r",  # (min,max)
        "z",  # (min,max)
        "beamPos",  # (x,y)
    ],
    defaults=[None] * 7 + [(None, None)] * 5,
)

TrackParamsEstimationConfig = namedtuple(
    "TrackParamsEstimationConfig",
    [
        "deltaR",  # (min,max)
    ],
    defaults=[(None, None)],
)


@acts.examples.NamedTypeArgs(
    seedingAlgorithm=SeedingAlgorithm,
    truthSeedRanges=TruthSeedRanges,
    particleSmearingSigmas=ParticleSmearingSigmas,
    seedfinderConfigArg=SeedfinderConfigArg,
    trackParamsEstimationConfig=TrackParamsEstimationConfig,
    logLevel=acts.logging.Level,
)
def addSeeding(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    geoSelectionConfigFile: Optional[Union[Path, str]] = None,
    seedingAlgorithm: SeedingAlgorithm = SeedingAlgorithm.Default,
    truthSeedRanges: TruthSeedRanges = TruthSeedRanges(),
    particleSmearingSigmas: ParticleSmearingSigmas = ParticleSmearingSigmas(),
    initialVarInflation: Optional[list] = None,
    seedfinderConfigArg: SeedfinderConfigArg = SeedfinderConfigArg(),
    trackParamsEstimationConfig: TrackParamsEstimationConfig = TrackParamsEstimationConfig(),
    inputParticles: str = "particles_initial",
    outputDirRoot: Optional[Union[Path, str]] = None,
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
    seedingAlgorithm : SeedingAlgorithm, Default
        seeding algorithm to use: one of Default (no truth information used), TruthSmeared, TruthEstimated
    truthSeedRanges : TruthSeedRanges(rho, z, phi, eta, absEta, pt, nHits)
        TruthSeedSelector configuration. Each range is specified as a tuple of (min,max).
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TruthSeedSelector.hpp
    particleSmearingSigmas : ParticleSmearingSigmas(d0, d0PtA, d0PtB, z0, z0PtA, z0PtB, t0, phi, theta, pRel)
        ParticleSmearing configuration.
        Defaults specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/ParticleSmearing.hpp
    initialVarInflation : list
        List of 6 scale factors to inflate the initial covariance matrix
        Defaults (all 1) specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/ParticleSmearing.hpp
    seedfinderConfigArg : SeedfinderConfigArg(maxSeedsPerSpM, cotThetaMax, sigmaScattering, radLengthPerSeed, minPt, bFieldInZ, impactMax, deltaR, collisionRegion, r, z, beamPos)
        SeedfinderConfig settings. deltaR, collisionRegion, r, z are ranges specified as a tuple of (min,max). beamPos is specified as (x,y).
        Defaults specified in Core/include/Acts/Seeding/SeedfinderConfig.hpp
    trackParamsEstimationConfig : TrackParamsEstimationConfig(deltaR)
        TrackParamsEstimationAlgorithm configuration. Currently only deltaR=(min,max) range specified here.
        Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp
    inputParticles : str, "particles_initial"
        input particles name in the WhiteBoard
    outputDirRoot : Path|str, path, None
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

    if int(customLogLevel()) <= int(acts.logging.DEBUG):
        acts.examples.dump_args_calls(locals())
    logger = acts.logging.getLogger("addSeeding")

    if truthSeedRanges is not None:
        selAlg = acts.examples.TruthSeedSelector(
            **acts.examples.defaultKWArgs(
                ptMin=truthSeedRanges.pt[0],
                ptMax=truthSeedRanges.pt[1],
                etaMin=truthSeedRanges.eta[0],
                etaMax=truthSeedRanges.eta[1],
                nHitsMin=truthSeedRanges.nHits[0],
                nHitsMax=truthSeedRanges.nHits[1],
                rhoMin=truthSeedRanges.rho[0],
                rhoMax=truthSeedRanges.rho[1],
                zMin=truthSeedRanges.z[0],
                zMax=truthSeedRanges.z[1],
                phiMin=truthSeedRanges.phi[0],
                phiMax=truthSeedRanges.phi[1],
                absEtaMin=truthSeedRanges.absEta[0],
                absEtaMax=truthSeedRanges.absEta[1],
            ),
            level=customLogLevel(),
            inputParticles=inputParticles,
            inputMeasurementParticlesMap="measurement_particles_map",
            outputParticles="truth_seeds_selected",
        )
        s.addAlgorithm(selAlg)
        selectedParticles = selAlg.config.outputParticles
    else:
        selectedParticles = inputParticles

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if seedingAlgorithm == SeedingAlgorithm.TruthSmeared:
        rnd = rnd or acts.examples.RandomNumbers(seed=42)
        logger.info("Using smeared truth particles for seeding")
        # Run particle smearing
        ptclSmear = acts.examples.ParticleSmearing(
            level=customLogLevel(),
            inputParticles=selectedParticles,
            outputTrackParameters="estimatedparameters",
            randomNumbers=rnd,
            # gaussian sigmas to smear particle parameters
            **acts.examples.defaultKWArgs(
                sigmaD0=particleSmearingSigmas.d0,
                sigmaD0PtA=particleSmearingSigmas.d0PtA,
                sigmaD0PtB=particleSmearingSigmas.d0PtB,
                sigmaZ0=particleSmearingSigmas.z0,
                sigmaZ0PtA=particleSmearingSigmas.z0PtA,
                sigmaZ0PtB=particleSmearingSigmas.z0PtB,
                sigmaT0=particleSmearingSigmas.t0,
                sigmaPhi=particleSmearingSigmas.phi,
                sigmaTheta=particleSmearingSigmas.theta,
                sigmaPRel=particleSmearingSigmas.pRel,
                initialVarInflation=initialVarInflation,
            ),
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
                inputParticles=selectedParticles,
                inputMeasurementParticlesMap=selAlg.config.inputMeasurementParticlesMap,
                outputProtoTracks="prototracks",
            )
            s.addAlgorithm(truthTrackFinder)
            inputProtoTracks = truthTrackFinder.config.outputProtoTracks
            inputSeeds = ""
        elif seedingAlgorithm == SeedingAlgorithm.Default:
            logger.info("Using default seeding")
            # Use seeding
            seedFinderConfig = acts.SeedfinderConfig(
                **acts.examples.defaultKWArgs(
                    rMin=seedfinderConfigArg.r[0],
                    rMax=seedfinderConfigArg.r[1],
                    deltaRMin=seedfinderConfigArg.deltaR[0],
                    deltaRMax=seedfinderConfigArg.deltaR[1],
                    deltaRMinTopSP=seedfinderConfigArg.deltaR[0],
                    deltaRMinBottomSP=seedfinderConfigArg.deltaR[0],
                    deltaRMaxTopSP=seedfinderConfigArg.deltaR[1],
                    deltaRMaxBottomSP=seedfinderConfigArg.deltaR[1],
                    collisionRegionMin=seedfinderConfigArg.collisionRegion[0],
                    collisionRegionMax=seedfinderConfigArg.collisionRegion[1],
                    zMin=seedfinderConfigArg.z[0],
                    zMax=seedfinderConfigArg.z[1],
                    maxSeedsPerSpM=seedfinderConfigArg.maxSeedsPerSpM,
                    cotThetaMax=seedfinderConfigArg.cotThetaMax,
                    sigmaScattering=seedfinderConfigArg.sigmaScattering,
                    radLengthPerSeed=seedfinderConfigArg.radLengthPerSeed,
                    minPt=seedfinderConfigArg.minPt,
                    bFieldInZ=seedfinderConfigArg.bFieldInZ,
                    impactMax=seedfinderConfigArg.impactMax,
                    beamPos=(
                        None
                        if seedfinderConfigArg.beamPos is None
                        or all([x is None for x in seedfinderConfigArg.beamPos])
                        else acts.Vector2(
                            seedfinderConfigArg.beamPos[0] or 0.0,
                            seedfinderConfigArg.beamPos[1] or 0.0,
                        )
                    ),
                ),
            )

            seedFilterConfig = acts.SeedFilterConfig(
                maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
                deltaRMin=seedFinderConfig.deltaRMin,
            )

            gridConfig = acts.SpacePointGridConfig(
                bFieldInZ=seedFinderConfig.bFieldInZ,
                minPt=seedFinderConfig.minPt,
                rMax=seedFinderConfig.rMax,
                zMax=seedFinderConfig.zMax,
                zMin=seedFinderConfig.zMin,
                deltaRMax=seedFinderConfig.deltaRMax,
                cotThetaMax=seedFinderConfig.cotThetaMax,
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
        elif seedingAlgorithm == SeedingAlgorithm.Orthogonal:
            logger.info("Using orthogonal seeding")
            # Use seeding
            seedFinderConfig = acts.SeedFinderOrthogonalConfig(
                **acts.examples.defaultKWArgs(
                    rMin=seedfinderConfigArg.r[0],
                    rMax=seedfinderConfigArg.r[1],
                    deltaRMin=seedfinderConfigArg.deltaR[0],
                    deltaRMax=seedfinderConfigArg.deltaR[1],
                    collisionRegionMin=seedfinderConfigArg.collisionRegion[0],
                    collisionRegionMax=seedfinderConfigArg.collisionRegion[1],
                    zMin=seedfinderConfigArg.z[0],
                    zMax=seedfinderConfigArg.z[1],
                    maxSeedsPerSpM=seedfinderConfigArg.maxSeedsPerSpM,
                    cotThetaMax=seedfinderConfigArg.cotThetaMax,
                    sigmaScattering=seedfinderConfigArg.sigmaScattering,
                    radLengthPerSeed=seedfinderConfigArg.radLengthPerSeed,
                    minPt=seedfinderConfigArg.minPt,
                    bFieldInZ=seedfinderConfigArg.bFieldInZ,
                    impactMax=seedfinderConfigArg.impactMax,
                    beamPos=(
                        None
                        if seedfinderConfigArg.beamPos is None
                        or all([x is None for x in seedfinderConfigArg.beamPos])
                        else acts.Vector2(
                            seedfinderConfigArg.beamPos[0] or 0.0,
                            seedfinderConfigArg.beamPos[1] or 0.0,
                        )
                    ),
                ),
            )

            seedFilterConfig = acts.SeedFilterConfig(
                maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
                deltaRMin=seedFinderConfig.deltaRMin,
            )

            seedingAlg = acts.examples.SeedingOrthogonalAlgorithm(
                level=customLogLevel(acts.logging.VERBOSE),
                inputSpacePoints=[spAlg.config.outputSpacePoints],
                outputSeeds="seeds",
                outputProtoTracks="prototracks",
                seedFilterConfig=seedFilterConfig,
                seedFinderConfig=seedFinderConfig,
            )
            s.addAlgorithm(seedingAlg)
            inputProtoTracks = seedingAlg.config.outputProtoTracks
            inputSeeds = seedingAlg.config.outputSeeds
        else:
            logger.fatal("unknown seedingAlgorithm %s", seedingAlgorithm)

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
            **acts.examples.defaultKWArgs(
                initialVarInflation=initialVarInflation,
                deltaRMin=trackParamsEstimationConfig.deltaR[0],
                deltaRMax=trackParamsEstimationConfig.deltaR[1],
            ),
        )
        s.addAlgorithm(parEstimateAlg)

        if outputDirRoot is not None:
            outputDirRoot = Path(outputDirRoot)
            if not outputDirRoot.exists():
                outputDirRoot.mkdir()
            s.addWriter(
                acts.examples.TrackFinderPerformanceWriter(
                    level=customLogLevel(),
                    inputProtoTracks=inputProtoTracks,
                    inputParticles=selectedParticles,  # the original selected particles after digitization
                    inputMeasurementParticlesMap=selAlg.config.inputMeasurementParticlesMap,
                    filePath=str(outputDirRoot / "performance_seeding_trees.root"),
                )
            )

            s.addWriter(
                acts.examples.SeedingPerformanceWriter(
                    level=customLogLevel(acts.logging.DEBUG),
                    inputProtoTracks=inputProtoTracks,
                    inputParticles=selectedParticles,
                    inputMeasurementParticlesMap=selAlg.config.inputMeasurementParticlesMap,
                    filePath=str(outputDirRoot / "performance_seeding_hists.root"),
                )
            )

            s.addWriter(
                acts.examples.RootTrackParameterWriter(
                    level=customLogLevel(acts.logging.VERBOSE),
                    inputTrackParameters=parEstimateAlg.config.outputTrackParameters,
                    inputProtoTracks=parEstimateAlg.config.outputProtoTracks,
                    inputParticles=inputParticles,
                    inputSimHits="simhits",
                    inputMeasurementParticlesMap=selAlg.config.inputMeasurementParticlesMap,
                    inputMeasurementSimHitsMap="measurement_simhits_map",
                    filePath=str(outputDirRoot / "estimatedparams.root"),
                    treeName="estimatedparams",
                )
            )

    return s


def runSeeding(
    trackingGeometry,
    field,
    outputDir,
    s=None,
    seedingAlgorithm=SeedingAlgorithm.Default,
):

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
        TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-2.5, 2.5), nHits=(9, None)),
        ParticleSmearingSigmas(pRel=0.01),  # only used by SeedingAlgorithm.TruthSmeared
        SeedfinderConfigArg(
            r=(None, 200 * u.mm),  # rMin=default, 33mm
            deltaR=(1 * u.mm, 60 * u.mm),
            collisionRegion=(-250 * u.mm, 250 * u.mm),
            z=(-2000 * u.mm, 2000 * u.mm),
            maxSeedsPerSpM=1,
            sigmaScattering=50,
            radLengthPerSeed=0.1,
            minPt=500 * u.MeV,
            bFieldInZ=1.99724 * u.T,
            impactMax=3 * u.mm,
        ),
        acts.logging.VERBOSE,
        seedingAlgorithm=seedingAlgorithm,
        geoSelectionConfigFile=srcdir
        / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json",
        inputParticles="particles_final",  # use this to reproduce the original root_file_hashes.txt - remove to fix
        outputDirRoot=outputDir,
    )
    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Example script to run seed finding",
    )

    p.add_argument(
        "--algorithm",
        action=EnumAction,
        enum=SeedingAlgorithm,
        default=SeedingAlgorithm.Default,
        help="Select the seeding algorithm to use",
    )

    args = p.parse_args()

    # detector, trackingGeometry, _ = getOpenDataDetector()
    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runSeeding(
        trackingGeometry, field, outputDir=Path.cwd(), seedingAlgorithm=args.algorithm
    ).run()
