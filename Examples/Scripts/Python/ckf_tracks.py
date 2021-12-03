#!/usr/bin/env python3
from pathlib import Path
from typing import Optional
import argparse

from acts.examples import Sequencer, GenericDetector, RootParticleReader, TGeoDetector

import acts

from acts import UnitConstants as u


def runCKFTracks(
    trackingGeometry,
    decorators,
    geometrySelection: Optional[Path],
    digiConfigFile: Path,
    field,
    outputDir: Path,
    truthSmearedSeeded=False,
    truthEstimatedSeeded=False,
    outputCsv=True,
    inputParticlePath: Optional[Path] = None,
    s=None,
    loglevel=acts.logging.INFO,
    events=100,
    numParticles=4,
    numVertices=2,
    atlasDetector=False,
):
    outputDir.mkdir(exist_ok=True)

    s = s or Sequencer(
        events=events, numThreads=-1, outputDir=str(outputDir), logLevel=loglevel
    )

    logger = acts.logging.getLogger("CKFExample")
    logger.setLevel(loglevel)

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    if inputParticlePath is None:
        logger.info("Generating particles using particle gun")
        if atlasDetector:
            particles = acts.examples.ParametricParticleGenerator(
                p=(0.4 * u.GeV, 50 * u.GeV),
                pTransverse=True,
                eta=(-4, 4),
                phi=(0, 360 * u.degree),
                randomizeCharge=True,
                numParticles=numParticles,
            )
        else:
            particles = acts.examples.ParametricParticleGenerator(
                p=(1 * u.GeV, 10 * u.GeV),
                eta=(-2, 2),
                phi=(0, 360 * u.degree),
                randomizeCharge=True,
                numParticles=numParticles,
            )

        evGen = acts.examples.EventGenerator(
            level=loglevel,
            generators=[
                acts.examples.EventGenerator.Generator(
                    multiplicity=acts.examples.FixedMultiplicityGenerator(
                        n=numVertices
                    ),
                    vertex=acts.examples.GaussianVertexGenerator(
                        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                    ),
                    particles=particles,
                )
            ],
            outputParticles="particles_input",
            randomNumbers=rnd,
        )
        s.addReader(evGen)
        inputParticles = evGen.config.outputParticles
    else:
        logger.info("Reading particles from %s", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        inputParticles = "particles_read"
        s.addReader(
            RootParticleReader(
                level=loglevel,
                filePath=str(inputParticlePath.resolve()),
                particleCollection=inputParticles,
                orderedEvents=False,
            )
        )

    # Selector
    selector = acts.examples.ParticleSelector(
        level=loglevel,
        inputParticles=inputParticles,
        outputParticles="particles_selected",
    )
    s.addAlgorithm(selector)

    # Simulation
    simAlg = acts.examples.FatrasSimulation(
        level=loglevel,
        inputParticles=selector.config.outputParticles,
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
    )
    s.addAlgorithm(simAlg)

    # Run the sim hits smearing
    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(str(digiConfigFile)),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits=simAlg.config.outputSimHits,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, loglevel)
    s.addAlgorithm(digiAlg)

    # Run the particle selection
    # The pre-selection will select truth particles satisfying provided criteria
    # from all particles read in by particle reader for further processing. It
    # has no impact on the truth hits themselves
    selAlg = acts.examples.TruthSeedSelector(
        level=loglevel,
        ptMin=500 * u.MeV,
        nHitsMin=9,
        inputParticles=simAlg.config.outputParticlesInitial,
        inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
        outputParticles="particles_seed_selected",
    )
    s.addAlgorithm(selAlg)

    inputParticles = selAlg.config.outputParticles

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if truthSmearedSeeded:
        logger.info("Using smeared truth particles for seeding")
        # Run particle smearing
        ptclSmear = acts.examples.ParticleSmearing(
            level=loglevel,
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
        outputTrackParameters = ptclSmear.config.outputTrackParameters
        s.addAlgorithm(ptclSmear)
    else:
        # Create space points
        spAlg = acts.examples.SpacePointMaker(
            level=loglevel,
            inputSourceLinks=digiAlg.config.outputSourceLinks,
            inputMeasurements=digiAlg.config.outputMeasurements,
            outputSpacePoints="spacepoints",
            trackingGeometry=trackingGeometry,
            geometrySelection=acts.examples.readJsonGeometryList(
                str(geometrySelection)
            ),
        )
        s.addAlgorithm(spAlg)

        # Run either: truth track finding or seeding
        if truthEstimatedSeeded:
            logger.info("Using truth track finding from space points for seeding")
            # Use truth tracking
            truthTrackFinder = acts.examples.TruthTrackFinder(
                level=loglevel,
                inputParticles=inputParticles,
                inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
                outputProtoTracks="prototracks",
            )
            s.addAlgorithm(truthTrackFinder)
            inputProtoTracks = truthTrackFinder.config.outputProtoTracks
            inputSeeds = ""
        else:
            logger.info("Using seeding")
            # Use seeding
            if atlasDetector:
                gridConfig = acts.SpacePointGridConfig(
                    bFieldInZ=1.99724 * u.T,
                    minPt=500 * u.MeV,
                    rMax=300 * u.mm,
                    zMax=3000 * u.mm,
                    zMin=-3000 * u.mm,
                    deltaRMax=100 * u.mm,
                    cotThetaMax=27.3,  # 4.0 eta
                )
            else:
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
                collisionRegionMin=-250 * u.mm,
                collisionRegionMax=250 * u.mm,
                zMin=gridConfig.zMin,
                zMax=gridConfig.zMax,
                maxSeedsPerSpM=seedFilterConfig.maxSeedsPerSpM,
                cotThetaMax=gridConfig.cotThetaMax,
                sigmaScattering=(5 if atlasDetector else 50),
                radLengthPerSeed=(0.5 if atlasDetector else 0.1),
                minPt=gridConfig.minPt,
                bFieldInZ=gridConfig.bFieldInZ,
                beamPos=acts.Vector2(0 * u.mm, 0 * u.mm),
                impactMax=3 * u.mm,
            )
            seeding = acts.examples.SeedingAlgorithm(
                level=loglevel,
                inputSpacePoints=[spAlg.config.outputSpacePoints],
                outputSeeds="seeds",
                outputProtoTracks="prototracks",
                gridConfig=gridConfig,
                seedFilterConfig=seedFilterConfig,
                seedFinderConfig=seedFinderConfig,
            )
            s.addAlgorithm(seeding)
            inputProtoTracks = seeding.config.outputProtoTracks
            inputSeeds = seeding.config.outputSeeds

        # Write truth track finding / seeding performance
        trackFinderPerformanceWriter = acts.examples.TrackFinderPerformanceWriter(
            level=loglevel,
            inputProtoTracks=inputProtoTracks,
            inputParticles=inputParticles,  # the original selected particles after digitization
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            filePath=str(outputDir / "performance_seeding_trees.root"),
        )
        s.addWriter(trackFinderPerformanceWriter)

        # Estimate track parameters from seeds
        paramEstimation = acts.examples.TrackParamsEstimationAlgorithm(
            level=loglevel,
            inputSeeds=inputSeeds,
            inputProtoTracks=inputProtoTracks,
            inputSpacePoints=[spAlg.config.outputSpacePoints],
            inputSourceLinks=digiCfg.outputSourceLinks,
            outputTrackParameters="estimatedparameters",
            outputProtoTracks="prototracks_estimated",
            trackingGeometry=trackingGeometry,
            magneticField=field,
            bFieldMin=0.1 * u.T,
            deltaRMax=100.0 * u.mm,
            deltaRMin=10.0 * u.mm,
            sigmaLoc0=25.0 * u.um,
            sigmaLoc1=100.0 * u.um,
            sigmaPhi=0.02 * u.degree,
            sigmaTheta=0.02 * u.degree,
            sigmaQOverP=0.1 / 1.0 * u.GeV,
            sigmaT0=1400.0 * u.s,
            initialVarInflation=[1, 1, 1, 1, 1, 1],
        )
        s.addAlgorithm(paramEstimation)
        outputTrackParameters = paramEstimation.config.outputTrackParameters

    # Setup the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=loglevel,
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [(acts.GeometryIdentifier(), ([], [15.0], [10]))]
        ),
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputSourceLinks=digiAlg.config.outputSourceLinks,
        inputInitialTrackParameters=outputTrackParameters,
        outputTrajectories="trajectories",
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field
        ),
    )
    s.addAlgorithm(trackFinder)

    # write track states from CKF
    trackStatesWriter = acts.examples.RootTrajectoryStatesWriter(
        level=loglevel,
        inputTrajectories=trackFinder.config.outputTrajectories,
        # @note The full particles collection is used here to avoid lots of warnings
        # since the unselected CKF track might have a majority particle not in the
        # filtered particle collection. This could be avoided when a seperate track
        # selection algorithm is used.
        inputParticles=selector.config.outputParticles,
        inputSimHits=simAlg.config.outputSimHits,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
        filePath=str(outputDir / "trackstates_ckf.root"),
        treeName="trackstates",
    )
    s.addWriter(trackStatesWriter)

    # write track summary from CKF
    trackSummaryWriter = acts.examples.RootTrajectorySummaryWriter(
        level=loglevel,
        inputTrajectories=trackFinder.config.outputTrajectories,
        # @note The full particles collection is used here to avoid lots of warnings
        # since the unselected CKF track might have a majority particle not in the
        # filtered particle collection. This could be avoided when a seperate track
        # selection algorithm is used.
        inputParticles=selector.config.outputParticles,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        filePath=str(outputDir / "tracksummary_ckf.root"),
        treeName="tracksummary",
    )
    s.addWriter(trackSummaryWriter)

    # Write CKF performance data
    ckfPerfWriter = acts.examples.CKFPerformanceWriter(
        level=loglevel,
        inputParticles=inputParticles,
        inputTrajectories=trackFinder.config.outputTrajectories,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        # The bottom seed could be the first, second or third hits on the truth track
        nMeasurementsMin=selAlg.config.nHitsMin - 3,
        ptMin=0.4 * u.GeV,
        filePath=str(outputDir / "performance_ckf.root"),
    )
    s.addWriter(ckfPerfWriter)

    if outputCsv:
        csv_dir = outputDir / "csv"
        csv_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Writing CSV files")
        csvMTJWriter = acts.examples.CsvMultiTrajectoryWriter(
            level=loglevel,
            inputTrajectories=trackFinder.config.outputTrajectories,
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            outputDir=str(csv_dir),
        )
        s.addWriter(csvMTJWriter)

    return s


def getITkGeometry(
    geo_dir: Path,
    mat_input_file: Optional[Path],
    geo_tgeo_jsonconfig: Path,
    geo_tgeo_filename: Path,
    loglevel=acts.logging.INFO,
):
    Volume = TGeoDetector.Config.Volume
    LayerTriplet = TGeoDetector.Config.LayerTriplet

    logger = acts.logging.getLogger("getITkGeometry")
    logger.setLevel(loglevel)

    matDeco = None
    if mat_input_file is not None:
        logger.info("Adding material from %s", mat_input_file.absolute())
        matDeco = acts.IMaterialDecorator.fromFile(
            mat_input_file,
            level=acts.logging.INFO,
        )

    return TGeoDetector.create(
        jsonFile=str(geo_tgeo_jsonconfig),
        fileName=str(geo_tgeo_filename),
        surfaceLogLevel=loglevel,
        layerLogLevel=loglevel,
        volumeLogLevel=loglevel,
        mdecorator=matDeco,
    )


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Example script to run CKF tracking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-l",
        "--loglevel",
        type=int,
        default=2,
        help="The output log level. Please set the wished number (0 = VERBOSE, 1 = DEBUG, 2 = INFO, 3 = WARNING, 4 = ERROR, 5 = FATAL).",
    )
    p.add_argument(
        "-n", "--events", type=int, default=100, help="The number of events to process."
    )
    p.add_argument(
        "--gen-nparticles",
        type=int,
        default=4,
        help="Number of generated particles per event.",
    )
    p.add_argument(
        "--gen-nvertices",
        type=int,
        default=2,
        help="Number of generated vertices.",
    )
    p.add_argument(
        "--generic-detector",
        action="store_true",
        default=True,
        help="Setup Generic Detector.",
    )
    p.add_argument(
        "-a",
        "--atlas",
        action="store_true",
        help="Setup ATLAS Detector using the ITk standalone geometry specified in --geo-dir. Get in touch if you don't have this.",
    )
    p.add_argument(
        "-g",
        "--geo-dir",
        help="Input directory containing geometry.",
    )
    p.add_argument("-i", "--input-particles", type=Path, help="Input particles path.")
    p.add_argument(
        "-o",
        "--output-dir",
        default=Path.cwd(),
        type=Path,
        help="Directory to write outputs to.",
    )
    p.add_argument(
        "--bf-constant-tesla",
        type=float,
        default=2,
        help="Set a constant magnetic field vector in Tesla.",
    )
    p.add_argument(
        "--output-trajectories-csv",
        action="store_true",
        help="output trajectories to CSV files",
    )
    p.add_argument(
        "--no-material",
        action="store_true",
        help="Decorate material to the geometry (--atlas or --build-itk-geometry only)",
    )
    p.add_argument(
        "--mat-input-file",
        type=Path,
        help="Name of the material map input file, supported: '.json', '.cbor' or '.root'.",
    )
    p.add_argument(
        "--ckf-truth-smeared-seeds",
        action="store_true",
        help="Use track parameters smeared from truth particles for steering CKF",
    )
    p.add_argument(
        "--ckf-truth-estimated-seeds",
        action="store_true",
        help="Use track parameters estimated from truth tracks for steering CKF",
    )
    p.add_argument(
        "--geo-selection-config-file",
        type=Path,
        help="Json file for space point geometry selection",
    )
    p.add_argument(
        "--build-itk-geometry",
        action="store_true",
        help="use itk.buildITkGeometry to create TGeo geometry (skips --geo-tgeo-*, defaults --mat-input-file)",
    )
    p.add_argument(
        "--geo-tgeo-jsonconfig",
        type=Path,
        help="TGeo Json config file name (--atlas only).",
    )
    p.add_argument(
        "--geo-tgeo-filename",
        type=Path,
        help="TGeo Root file name (--atlas only).",
    )
    p.add_argument(
        "--digi-config-file",
        type=Path,
        help="Configuration (.json) file for digitization description",
    )

    args = p.parse_args()

    if args.atlas:
        args.generic_detector = False

    if args.generic_detector:
        if args.geo_dir is None:
            args.geo_dir = Path(__file__).resolve().parent.parent.parent.parent
        if not args.ckf_truth_smeared_seeds and args.geo_selection_config_file is None:
            args.geo_selection_config_file = (
                args.geo_dir
                / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
            )
        if args.digi_config_file is None:
            args.digi_config_file = (
                args.geo_dir
                / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
            )
        if not args.build_itk_geometry:
            detector, trackingGeometry, decorators = GenericDetector.create()

    if args.atlas:
        if args.geo_dir is None:
            args.geo_dir = Path("acts-detector-examples")

        if not args.ckf_truth_smeared_seeds and args.geo_selection_config_file is None:
            args.geo_selection_config_file = (
                args.geo_dir / "atlas/itk-hgtd/geoSelection-ITk.json"
            )
        if args.digi_config_file is None:
            args.digi_config_file = (
                args.geo_dir / "atlas/itk-hgtd/itk-smearing-config.json"
            )

        if not args.build_itk_geometry:
            if not args.no_material and args.mat_input_file is None:
                args.mat_input_file = (
                    args.geo_dir / "atlas/itk-hgtd/material-maps-ITk-HGTD.json"
                )
            if args.geo_tgeo_jsonconfig is None:
                args.geo_tgeo_jsonconfig = (
                    args.geo_dir / "atlas/itk-hgtd/tgeo-atlas-itk-hgtd.json"
                )
            if args.geo_tgeo_filename is None:
                args.geo_tgeo_filename = (
                    args.geo_dir / "atlas/itk-hgtd/ATLAS-ITk-HGTD.tgeo.root"
                )
            detector, trackingGeometry, decorators = getITkGeometry(
                args.geo_dir,
                mat_input_file=args.mat_input_file,
                geo_tgeo_jsonconfig=args.geo_tgeo_jsonconfig,
                geo_tgeo_filename=args.geo_tgeo_filename,
                loglevel=acts.logging.Level(args.loglevel),
            )

    if args.build_itk_geometry:
        from itk import buildITkGeometry

        detector, trackingGeometry, decorators = buildITkGeometry(
            args.geo_dir, material=not args.no_material
        )

    field = acts.ConstantBField(acts.Vector3(0, 0, args.bf_constant_tesla * u.T))

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=args.geo_selection_config_file,
        digiConfigFile=args.digi_config_file,
        outputCsv=args.output_trajectories_csv,
        truthSmearedSeeded=args.ckf_truth_smeared_seeds,
        truthEstimatedSeeded=args.ckf_truth_estimated_seeds,
        inputParticlePath=args.input_particles,
        outputDir=args.output_dir,
        loglevel=acts.logging.Level(args.loglevel),
        events=args.events,
        numParticles=args.gen_nparticles,
        numVertices=args.gen_nvertices,
        atlasDetector=args.atlas,
    ).run()
